from collections.abc import Iterable
from copy import deepcopy
import re
import math
from typing import Any, Union, Generator
from functools import partial


MonoCompat = Union[int, str, "Monomial"]
PolyCompat = Union[int, str, "Monomial", "Polynomial"]

_super_list = ["⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹", "⁻"]
_super = str.maketrans("0123456789-", "".join(_super_list))


class Monomial:
    """Represents a monomial expression. A monomial is a product of a number and one or more variables raised to a power.
    Coeffients may be any real number.
    Powers may be any integer."""

    def __init__(self, coeffient: int, exponent: int = 1):
        self.coeffient = coeffient
        self.exponent = exponent

    @classmethod
    def read(cls, string: str) -> "Monomial":
        """Reads a monomial from a string. The string should be in the form of a number followed by one or more variables and their respective powers, ie."""
        string = str(string).lower().replace(" ", "")
        regex = re.compile(
            r"(?P<coef>[+-]?[0-9]*)?(?P<variable>([a-zA-Z](\^[0-9]+|((\u00b2|\u00b3|\u00b9|\u2070|\u2074|\u2075|\u2076|\u2077|\u2078|\u2079)+)?))?)?"
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

        var = match.group("variable")
        if not var:
            return M(coef, 0)
        var = var[1:].lstrip("^")
        if not var:
            return M(coef, 1)

        normal = "".join(
            map(lambda c: str(_super_list.index(c)) if c in _super_list else c, var)
        )

        return M(coef, int(normal))

    def solve(self, value: float = 0) -> Generator[complex, None, None]:
        omega_n = math.e ** (2j * math.pi / self.exponent)

        for k in range(self.exponent):
            yield (value / self.coeffient) ** (1 / self.exponent) * omega_n**k

    def derivative(self, n: int = 1) -> "Monomial":
        if n == 0:
            return self

        if self.is_constant:
            return M(0)

        return Monomial.derivative(
            M(self.coeffient * self.exponent, self.exponent - 1),
            n - 1,
        )

    @property
    def is_constant(self) -> bool:
        return self.exponent == 0

    def at(self, value: float):
        return self.coeffient * value**self.exponent

    def __call__(self, value: int):
        return self.at(value)

    def __eq__(self, other: Any) -> bool:
        return (
            isinstance(other, M)
            and self.coeffient == other.coeffient
            and self.exponent == other.exponent
        )

    def __add__(self, other: MonoCompat) -> "Monomial | Polynomial":
        other = ensure_monomial(other)
        if self.exponent == other.exponent:
            return M(self.coeffient + other.coeffient, self.exponent)
        return P([self, other])

    def __sub__(self, other: MonoCompat) -> "Monomial | Polynomial":
        other = ensure_monomial(other)
        if self.exponent == other.exponent:
            return M(self.coeffient - other.coeffient, self.exponent)
        return P([self, M(-other.coeffient, other.exponent)])

    def __mul__(self, other: MonoCompat) -> "Monomial":
        other = ensure_monomial(other)
        return M(self.coeffient * other.coeffient, self.exponent + other.exponent)

    def __div__(self, other: MonoCompat) -> "Monomial":
        other = ensure_monomial(other)
        ratio = self.coeffient / other.coeffient
        if not ratio.is_integer():
            raise ValueError("Coefficients not divisible")
        return M(int(ratio), self.exponent - other.exponent)

    def __truediv__(self, other: "Monomial") -> "Monomial":
        return self.__div__(other)

    def __str__(self) -> str:
        match self.coeffient:
            case 1:
                coef_str = ""
            case -1:
                coef_str = "-"
            case _:
                coef_str = str(self.coeffient)

        match self.exponent:
            case 1:
                exp_str = "x"
            case 0:
                exp_str = ""
            case _:
                exp_str = "x" + str(self.exponent).translate(_super)

        return coef_str + exp_str

    def __repr__(self) -> str:
        return self.__str__()


M = Monomial


class Polynomial:
    def __call__(self, value: float):
        at = partial(Monomial.at, value=value)
        return sum(map(at, self.elements))

    def zero(self, value: float = 0, *, accuracy: int = 100, starting_guess: float = 0):
        # gets a zero of a function using the Newton-Raphson
        guess = starting_guess
        slope = self.derivative()

        for _ in range(accuracy):
            try:
                guess -= self(guess) / slope(guess)
            except ZeroDivisionError:
                break

        return guess

    def solve(self, value: float = 0):
        """Solves a polynomial equation using the Newton-Raphson method and synthetic division."""
        
        zero = self.zero(value=value)
        new = self.synthetic_division(zero)[0]
        print("z", zero, "n", new)
        
        zero2 = new.zero(value=value)
        new2 = new.synthetic_division(zero2)[0]
        print("z", zero2, "n", new2)

        zero3 = new2.zero(value=value)
        new3 = new2.synthetic_division(zero3)[0]
        print("z", zero3, "n", new3)
        

    def synthetic_division(self, value: float) -> tuple["Polynomial", int]:
        """Finds the quotient and remainder of a polynomial division using synthetic division, in the form of P / (x - value)"""
        under = [self.elements[0].coeffient]

        for i, e in enumerate(self.elements[1:], start=1):
            under.append(under[i - 1] * value + e.coeffient)

        return P([M(coef, exp) for coef, exp in zip(under, range(self.degree)[::-1])]), under[-1]  # type: ignore

    @property
    def degree(self) -> int:
        return max(e.exponent for e in self.elements)

    def trim_zero(self):
        for i, e in enumerate(self.elements):
            if e.coeffient == 0:
                del self.elements[i]

    def derivative(self, n: int = 1) -> "Polynomial":
        p = P(map(partial(Monomial.derivative, n=n), self.elements))
        return p

    def slope_line(self, x: float = 0) -> tuple[float, float]:
        y = self(x)
        slope = self.derivative()(x)
        b = y - slope * x
        return slope, b

    @classmethod
    def read(cls, string: str) -> "Polynomial":
        segments = msplit(string, ["+", "-"], retain=["-"])

        elements = list(map(M.read, segments))

        return cls(elements)
    
    def euc_div(self, other: "Polynomial") -> tuple["Polynomial", "Polynomial"]:
        a = deepcopy(self)
        b = deepcopy(other)

        q = P([M(0)])
        r = a
        d = other.degree
        c = other.lc

        while r.degree >= d:
            s =  M(1, r.degree - d) * M(r.lc / c, 1)
            q += s
            print(b, s, b*s)
            r = r - b * s
            print(r.degree, d)
        return q, r



    @property
    def lc(self) -> float:
        # sort by exponent
        self.elements.sort(key=lambda e: e.exponent, reverse=True)
        return self.elements[0].coeffient

    def __init__(self, elements: Iterable[Monomial]):
        self.elements = list(elements)
        self.elements.sort(key=lambda e: e.exponent, reverse=True)

    def __str__(self) -> str:
        if not self.elements:
            return "0"
        final = str(self.elements[0])
        for element in self.elements[1:]:
            if element.coeffient > 0:
                final += f" + {element}"
            else:
                final += f" - {Monomial(-element.coeffient, element.exponent)}"

        return final

    def __repr__(self) -> str:
        return self.__str__()

    def __add__(self, other: PolyCompat):
        other = ensure_polynomial(other)
        copy = deepcopy(self)
        for oe in other.elements:
            for se in copy.elements:
                if se.exponent == oe.exponent:
                    se.coeffient += oe.coeffient
                    break
            else:
                copy.elements.append(oe)

        return copy

    def __sub__(self, other: PolyCompat):
        other = ensure_polynomial(other)
        copy = deepcopy(self)
        for oe in other.elements:
            for se in copy.elements:
                if se.exponent == oe.exponent:
                    se.coeffient -= oe.coeffient
                    break
            else:
                copy.elements.append(M(-oe.coeffient, oe.exponent))

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


def msplit(string: str, delimiters: list[str], retain: list[str] | None = None):
    """Splits a string by multiple delimiters, but retains the delimiters in the right side of the split."""
    if retain is None:
        retain = []

    segments, current = list(), str()
    current = ""
    for char in string:
        if char in delimiters:
            segments.append(current)
            current = char
        elif char in retain:
            current += char
        else:
            current += char

    segments.append(current)
    return segments


# p = M.read("14x^2")
# q = M.read("7x")
# print(p)
# print(q)
# print(p + q)
# print(p - q)
# print(p * q)
# print(p / q)

p = P.read("x^2 - 5x -2")
q = P.read("x - 1")
print(p.euc_div(q))