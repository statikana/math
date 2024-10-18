from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum
import re
from turtle import right
from typing import Generic, TypeVar, Union
from operator import add, sub, mul, truediv

from matplotlib.pyplot import cla
from numpy import var

class Operator(Enum):
    ADD = "+"
    SUB = "-"
    MUL = "*"
    DIV = "/"

_super_list = ["⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"]
_super = str.maketrans("0123456789", "".join(_super_list))

_sub_list = ["₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"]
_sub = str.maketrans("0123456789", "".join(_sub_list))

_super_plus = "⁺"
_super_minus = "⁻"

class Monomial:
    """Represents a monomial expression. A monomial is a product of a number and one or more variables raised to a power.
    Coeffients may be any real number.
    Powers may be any integer."""
    def __init__(self, coeffient: float, **variables: int):
        self.coeffient = coeffient
        self.variables = variables
    
    @classmethod
    def read(cls, string: str) -> "Monomial":
        """Reads a monomial from a string. The string should be in the form of a number followed by one or more variables and their respective powers, ie. """
        string = string.lower().replace(" ", "")
        regex = re.compile(r"(?P<coef>[+-]?[0-9\.]*)?(?P<variables>([a-zA-Z](\^[0-9]+|((\u00b2|\u00b3|\u00b9|\u2070|\u2074|\u2075|\u2076|\u2077|\u2078|\u2079)+)?))*)?")
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
                coef = float(coef_match or 1)

        if not match.group("variables"):
            return M(coef)
        
        variables = {}
        current_var = None
        current_num = 1
        for char in string[len(coef_match or ""):]:
            if re.fullmatch(r"[a-zA-Z]", char):
                if current_var is not None:
                    variables.setdefault(current_var, 0)
                    variables[current_var] += current_num if current_num is not None else 1
                current_var = char
                current_num = None
            elif char == "^":
                pass
            elif re.fullmatch(r"[0-9]", char):
                if current_num is None:
                    current_num = 0
                current_num *= 10
                current_num += int(char)
            elif re.fullmatch(r"(\u00b2|\u00b3|\u00b9|\u2070|\u2074|\u2075|\u2076|\u2077|\u2078|\u2079)", char):
                if current_num is None:
                    current_num = 0
                current_num *= 10
                current_num += _super_list.index(char)
            else:
                raise
        if current_var:
            variables.setdefault(current_var, 0)
            variables[current_var] += current_num if current_num is not None else 1
        return M(coef, **variables)

    def n_variables(self) -> int:
        """How many unique variables are in this monomial."""
        return len(self.variables)
    
    @property
    def is_constant(self) -> bool:
        return not self.variables

    def __add__(self, other: "Monomial") -> "Monomial | Polynomial":
        if self.variables == other.variables:
            return M(self.coeffient + other.coeffient, **self.variables)
        return P([self, other])
    
    def __sub__(self, other: "Monomial") -> "Monomial | Polynomial":
        if self.variables == other.variables:
            return M(self.coeffient - other.coeffient, **self.variables)
        return P([self, other])
    
    def __mul__(self, other: "Monomial") -> "Monomial":
        var_names = list(self.variables.keys()) + list(other.variables.keys())
        variables = {name: self.variables.get(name, 0) + other.variables.get(name, 0) for name in var_names}
        return M(self.coeffient * other.coeffient, **variables)
    
    def __div__(self, other: "Monomial") -> "Monomial":
        var_names = list(self.variables.keys()) + list(other.variables.keys())
        variables = {name: self.variables.get(name, 0) - other.variables.get(name, 0) for name in var_names}
        return M(self.coeffient / other.coeffient, **variables)
    
    def __truediv__(self, other: "Monomial") -> "Monomial":
        return self.__div__(other)
    
    def __str__(self) -> str:
        variable_names = sorted(self.variables.keys())
        sorted_vars = {name: self.variables[name] for name in variable_names}
        variables = "".join([f"{name}{((_super_minus if power < 0 else "") + str(abs(power)).translate(_super)) if power != 1 else ""}" for name, power in sorted_vars.items()])
        coef_str = int(self.coeffient) if self.coeffient.is_integer() else self.coeffient
        if self.coeffient == 1:
            coef_str = ""
        elif self.coeffient == -1:
            coef_str = "-"

        return f"{coef_str}{variables}"

M = Monomial


class Polynomial:
    def solutions(self, value: float = 0) -> list[float]:
        if not self.elements:
            raise ValueError("Cannot solve a polynomial with no elements")
        if self.n_variables > 1:
            raise ValueError("Cannot solve a polynomial with more than one variable")
    
        # at this point, we know there is only one variable
        # separate into factors ie
        # 30x^3 + 10x^2 + 5x
        # 5x(6x^2 + 2x)
        # 5x(2x(3x + 1))

        # 15x^2 + 5x + 2 = 0
        #              ^ ignore
        # 5x(3x + 1) + 2 = 0
        # 5x(3x+1) = -2
        # 5x = -2
        # 3x+1 = -2

        pass




    def first_factor
    
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
        self.elements = sorted(elements, key=lambda e: max(power for _, power in e.variables.items()), reverse=True)

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

P = Polynomial


class Equation:
    def __init__(self, left: Monomial | Polynomial, right: Monomial | Polynomial):
        if isinstance(left, Monomial):
            left = P([left])
        if isinstance(right, Monomial):
            right = P([right])
        self.left = left
        self.right = right
    
    @classmethod
    def read(cls, string: str) -> "Equation":
        left, right = string.split("=")
        return cls(P.read(left), P.read(right))

    def __str__(self) -> str:
        return f"{self.left} = {self.right}"
    

E = Equation

left = M.read("x^2")
right = M.read("x")
print(E(left, right))