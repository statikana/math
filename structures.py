from collections.abc import Iterable
from copy import deepcopy
import re
import math
from typing import Any, Union, Generator
from functools import partial


MonoCompat = Union[int, str, "Monomial"]
PolyCompat = Union[int, str, "Monomial", "Polynomial"]

_super_list = ["⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"]
_super = str.maketrans("0123456789", "".join(_super_list))
_super_minus = "⁻"


class Monomial:
    """Represents a monomial expression. A monomial is a product of a number and one or more variables raised to a power.
    Coeffients may be any real number.
    Powers may be any integer."""

    def __init__(self, coeffient: int, **variables: int):
        self.coeffient = coeffient
        self.variables = variables

    def clean(self) -> None:
        delk = []
        for k, v in self.variables.items():
            if v == 0:
                delk.append(k)
        for k in delk:
            del self.variables[k]

    @classmethod
    def read(cls, string: str) -> "Monomial":
        """Reads a monomial from a string. The string should be in the form of a number followed by one or more variables and their respective powers, ie."""
        string = str(string).lower().replace(" ", "")
        regex = re.compile(
            r"(?P<coef>[+-]?[0-9]*)?(?P<variables>([a-zA-Z](\^[0-9]+|((\u00b2|\u00b3|\u00b9|\u2070|\u2074|\u2075|\u2076|\u2077|\u2078|\u2079)+)?))*)?"
        )
        match = regex.fullmatch(string)
        if match is None:
            raise ValueError("Invalid monomial string")

        coef_match = match.group("coef")
        match coef_match:
            case None:
                coef = 1
            case "+":
                coef = 1
            case "-":
                coef = -1
            case _:
                coef = int(coef_match or 1)

        if not match.group("variables"):
            return M(coef)

        variables = {}
        current_var = None
        current_num = 1
        for char in string[len(coef_match or "") :]:
            if re.fullmatch(r"[a-zA-Z]", char):
                if current_var is not None:
                    variables.setdefault(current_var, 0)
                    variables[current_var] += (
                        current_num if current_num is not None else 1
                    )
                current_var = char
                current_num = None
            elif char == "^":
                pass
            elif re.fullmatch(r"[0-9]", char):
                if current_num is None:
                    current_num = 0
                current_num *= 10
                current_num += int(char)
            elif re.fullmatch(
                r"(\u00b2|\u00b3|\u00b9|\u2070|\u2074|\u2075|\u2076|\u2077|\u2078|\u2079)",
                char,
            ):
                if current_num is None:
                    current_num = 0
                current_num *= 10
                current_num += _super_list.index(char)
            else:
                raise
        if current_var:
            variables.setdefault(current_var, 0)
            variables[current_var] += current_num if current_num is not None else 1
        mono = M(coef, **variables)
        mono.clean()
        return mono

    def n_variables(self) -> int:
        """How many unique variables are in this monomial."""
        return len(self.variables)

    @property
    def is_constant(self) -> bool:
        self.clean()
        return not self.variables

    def solve(self, var: str, value: float = 0) -> Generator[complex, None, None]:
        power = self.variables[var]
        omega_n = math.e ** (2j * math.pi / power)

        for k in range(power):
            yield (value / self.coeffient) ** (1 / power) * omega_n**k

    def derivative(self, n: int = 1) -> "Monomial":
        if n == 0:
            return self
        self.clean()

        if self.is_constant:
            return M(0)

        if len(self.variables) > 1:
            raise ValueError()

        var = list(self.variables.keys())[0]
        return Monomial.derivative(
            M(self.coeffient * self.variables[var], **{var: self.variables[var] - 1}),
            n - 1,
        )

    @property
    def variable_list(self) -> list[str]:
        return list(self.variables.keys())

    def __call__(self, **values: float):
        new = deepcopy(self)
        for var in new.variable_list:
            if var in values:
                new.coeffient *= values[var] ** new.variables[var]
                del new.variables[var]

        return new

    def __eq__(self, other: Any) -> bool:
        return (
            isinstance(other, M)
            and self.coeffient == other.coeffient
            and self.variables == other.variables
        )

    def __add__(self, other: MonoCompat) -> "Monomial | Polynomial":
        other = ensure_monomial(other)
        if self.variables == other.variables:
            return M(self.coeffient + other.coeffient, **self.variables)
        return P([self, other])

    def __sub__(self, other: MonoCompat) -> "Monomial | Polynomial":
        other = ensure_monomial(other)
        if self.variables == other.variables:
            return M(self.coeffient - other.coeffient, **self.variables)
        return P([self, other])

    def __mul__(self, other: MonoCompat) -> "Monomial":
        other = ensure_monomial(other)
        var_names = list(self.variables.keys()) + list(other.variables.keys())
        variables = {
            name: self.variables.get(name, 0) + other.variables.get(name, 0)
            for name in var_names
        }
        return M(self.coeffient * other.coeffient, **variables)

    def __div__(self, other: MonoCompat) -> "Monomial":
        other = ensure_monomial(other)
        var_names = list(self.variables.keys()) + list(other.variables.keys())
        variables = {
            name: self.variables.get(name, 0) - other.variables.get(name, 0)
            for name in var_names
        }
        ratio = self.coeffient / other.coeffient
        if not ratio.is_integer():
            raise ValueError("Coefficients not divisible")
        return M(int(ratio), **variables)

    def __truediv__(self, other: "Monomial") -> "Monomial":
        return self.__div__(other)

    def __str__(self) -> str:
        variable_names = sorted(self.variables.keys())
        sorted_vars = {name: self.variables[name] for name in variable_names}
        variables = "".join(
            [
                f"{name}{((_super_minus if power < 0 else "") + str(abs(power)).translate(_super)) if power != 1 else ""}"
                for name, power in sorted_vars.items()
            ]
        )
        coef_str = int(self.coeffient)

        if all(power == 0 for _, power in self.variables.items()):
            variables = ""
        else:
            if self.coeffient == 1:
                coef_str = ""
            elif self.coeffient == -1:
                coef_str = "-"

        return f"{coef_str}{variables}"

    def __repr__(self) -> str:
        return self.__str__()


M = Monomial


class Polynomial:
    def __call__(self, value: int | None = None, **values):
        if self.n_variables > 1 and value is not None:
            raise ValueError(
                "Cannot use position arguments with a polynomial with more than one unique variable"
            )
        elif value is not None:
            var = self.elements[0].variable_list[0]
            values.update({var: value})
        new = deepcopy(self)
        for i, e in enumerate(new.elements):
            new.elements[i] = e(**values)

        new.clean()
        return new
    
    def slope(self, )

    def shell(self):
        if not self.elements:
            raise ValueError("Cannot shell a polynomial with no elements")
        if self.n_variables > 1:
            raise ValueError("Cannot shell a polynomial with more than one variable")

        var: str = list(self.elements[0].variables.keys())[0]

        current = deepcopy(self)
        shelled: list[tuple[M, P, int]] = []

        while len(current.elements) != 1:
            if current.elements[-1].is_constant:
                constant = current.elements[-1].coeffient
                current = P(current.elements[:-1])
            else:
                constant = 0

            factor = current.gc_mono(var)
            if factor.coeffient == 1 and not factor.variables:
                break

            current /= factor

            shelled.append((factor, current, constant))

        return shelled
    
    def solve(self, var: str, value: float = 0, accuracy: int = 10):
        guess = 5



    def gc_mono(self, var: str) -> "Monomial":
        for e in self.elements:
            e.clean()
        coefficients = [e.coeffient for e in self.elements]
        gcd_coef = math.gcd(*coefficients)
        gcv = self.elements[-1].variables[var]
        mono = M(gcd_coef, **{var: gcv})
        mono.clean()
        return mono

    def derivative(self, n: int = 1) -> "Polynomial":
        p = P(map(partial(Monomial.derivative, n=n), self.elements))
        p.clean()
        return p

    def clean(self):
        map(Monomial.clean, self.elements)
        self.elements = Polynomial(self.elements).elements
        self.combine_like_terms()

    def combine_like_terms(self):
        for i, element in enumerate(self.elements):
            element.clean()
            changes = 0
            for j, other in enumerate(self.elements[i + 1 :]):
                if element.variables == other.variables:
                    element.coeffient += other.coeffient
                    del self.elements[j + i + 1 - changes]
                    changes += 1

    @property
    def n_variables(self) -> int:
        """
        How many unique variables are in this polynomial.
        """
        unique = set()
        for element in self.elements:
            unique.update(element.variables.keys())
        return len(unique)

    @classmethod
    def read(cls, string: str) -> "Polynomial":
        # split on + or - to get monomials
        elements = []
        for element in string.split("+"):
            if "-" in element:
                elements.append(M.read(element.split("-")[0]))
                for subelement in element.split("-")[1:]:
                    elements.append(M.read("-" + subelement))
            else:
                elements.append(M.read(element))
        return cls(elements)

    def __init__(self, elements: Iterable[Monomial]):
        self.elements = sorted(
            elements,
            key=lambda e: (
                max(power for _, power in e.variables.items()) if e.variables else 0
            ),
            reverse=True,
        )

    def __str__(self) -> str:
        if not self.elements:
            return "0"
        final = str(self.elements[0])
        for element in self.elements[1:]:
            if element.coeffient > 0:
                final += f" + {element}"
            else:
                final += f" - {Monomial(-element.coeffient, **element.variables)}"

        return final

    def __repr__(self) -> str:
        return self.__str__()

    def __add__(self, other: PolyCompat):
        other = ensure_polynomial(other)
        copy = deepcopy(self)
        for oe in other.elements:
            for se in copy.elements:
                if se.variables == oe.variables:
                    se.coeffient += oe.coeffient
                    break
            else:
                copy.elements.append(oe)

        return copy

    def __sub__(self, other: PolyCompat):
        other = ensure_polynomial(other)
        copy = deepcopy(self)
        for oe in copy.elements:
            for se in other.elements:
                if se.variables == oe.variables:
                    se.coeffient -= oe.coeffient
                    break
            else:
                copy.elements.append(M(-oe.coeffient, **oe.variables))

        return copy

    def __mul__(self, other: PolyCompat):
        other = ensure_polynomial(other)
        copy = deepcopy(self)
        for se in copy.elements:
            for oe in other.elements:
                se *= oe

        return copy

    def __div__(self, other: PolyCompat):
        other = ensure_polynomial(other)
        copy = deepcopy(self)
        for i, se in enumerate(copy.elements):
            for oe in other.elements:
                copy.elements[i] = se / oe
        return copy

    def __truediv__(self, other: PolyCompat):
        return self.__div__(other)


P = Polynomial


def ensure_monomial(item: MonoCompat):
    if isinstance(item, (int, str)):
        return M.read(str(item))
    return item


def ensure_polynomial(item: PolyCompat):
    if isinstance(item, (int, str)):
        return P([M.read(str(item))])
    elif isinstance(item, M):
        return P([item])
    return item


# p = P.read("300x^10 + 3x^3 - 7x^2 + 20x + 2")
# print(p(x=0))
for s in M(10, x=4).solve("x", 10):
    print(s)
