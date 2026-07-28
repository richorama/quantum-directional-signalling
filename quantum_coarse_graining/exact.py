"""Small exact matrix algebra over the Gaussian rationals ``Q(i)``."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence, Tuple


@dataclass(frozen=True)
class Gaussian:
    real: Fraction = Fraction(0)
    imag: Fraction = Fraction(0)

    def __post_init__(self) -> None:
        object.__setattr__(self, "real", Fraction(self.real))
        object.__setattr__(self, "imag", Fraction(self.imag))

    def __add__(self, other: "Gaussian") -> "Gaussian":
        return Gaussian(self.real + other.real, self.imag + other.imag)

    def __sub__(self, other: "Gaussian") -> "Gaussian":
        return Gaussian(self.real - other.real, self.imag - other.imag)

    def __mul__(self, other: "Gaussian") -> "Gaussian":
        return Gaussian(
            self.real * other.real - self.imag * other.imag,
            self.real * other.imag + self.imag * other.real,
        )

    def __truediv__(self, other: "Gaussian") -> "Gaussian":
        denominator = other.norm2()
        if denominator == 0:
            raise ZeroDivisionError("division by zero")
        numerator = self * other.conjugate()
        return Gaussian(numerator.real / denominator, numerator.imag / denominator)

    def __neg__(self) -> "Gaussian":
        return Gaussian(-self.real, -self.imag)

    def conjugate(self) -> "Gaussian":
        return Gaussian(self.real, -self.imag)

    def norm2(self) -> Fraction:
        return self.real * self.real + self.imag * self.imag


ZERO = Gaussian()
ONE = Gaussian(1)
I = Gaussian(0, 1)
Matrix = Tuple[Tuple[Gaussian, ...], ...]


def gaussian(value: int | Fraction) -> Gaussian:
    return Gaussian(Fraction(value))


def matrix(rows: Sequence[Sequence[int | Fraction | Gaussian]]) -> Matrix:
    return tuple(
        tuple(value if isinstance(value, Gaussian) else gaussian(value) for value in row)
        for row in rows
    )


def zeros(rows: int, columns: int) -> Matrix:
    return tuple(tuple(ZERO for _ in range(columns)) for _ in range(rows))


def identity(dimension: int) -> Matrix:
    return tuple(
        tuple(ONE if row == column else ZERO for column in range(dimension))
        for row in range(dimension)
    )


def matrix_unit(dimension: int, row: int, column: int) -> Matrix:
    return tuple(
        tuple(
            ONE if (r == row and c == column) else ZERO
            for c in range(dimension)
        )
        for r in range(dimension)
    )


def dagger(value: Matrix) -> Matrix:
    return tuple(
        tuple(value[row][column].conjugate() for row in range(len(value)))
        for column in range(len(value[0]))
    )


def matmul(left: Matrix, right: Matrix) -> Matrix:
    if len(left[0]) != len(right):
        raise ValueError("matrix dimensions do not compose")
    return tuple(
        tuple(
            sum(
                (left[row][k] * right[k][column] for k in range(len(right))),
                ZERO,
            )
            for column in range(len(right[0]))
        )
        for row in range(len(left))
    )


def add(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(left[row][column] + right[row][column] for column in range(len(left[0])))
        for row in range(len(left))
    )


def subtract(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(left[row][column] - right[row][column] for column in range(len(left[0])))
        for row in range(len(left))
    )


def scale(scalar: Gaussian, value: Matrix) -> Matrix:
    return tuple(tuple(scalar * entry for entry in row) for row in value)


def kron(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(
            left[i][j] * right[r][s]
            for j in range(len(left[0]))
            for s in range(len(right[0]))
        )
        for i in range(len(left))
        for r in range(len(right))
    )


def frobenius_norm2(value: Matrix) -> Fraction:
    return sum((entry.norm2() for row in value for entry in row), Fraction(0))


def vectorize(value: Matrix) -> Tuple[Gaussian, ...]:
    return tuple(entry for row in value for entry in row)


def superoperator(images: Sequence[Matrix]) -> Matrix:
    columns = [vectorize(image) for image in images]
    return tuple(
        tuple(columns[column][row] for column in range(len(columns)))
        for row in range(len(columns[0]))
    )


def partial_trace_b(value: Matrix, dim_a: int, dim_b: int) -> Matrix:
    return tuple(
        tuple(
            sum(
                (
                    value[dim_b * a + b][dim_b * c + b]
                    for b in range(dim_b)
                ),
                ZERO,
            )
            for c in range(dim_a)
        )
        for a in range(dim_a)
    )


def conjugate_action(unitary: Matrix, operator: Matrix) -> Matrix:
    return matmul(matmul(unitary, operator), dagger(unitary))


def is_unitary(value: Matrix) -> bool:
    return matmul(dagger(value), value) == identity(len(value))


def permutation_unitary(permutation: Sequence[int]) -> Matrix:
    dimension = len(permutation)
    if sorted(permutation) != list(range(dimension)):
        raise ValueError("not a permutation")
    return tuple(
        tuple(ONE if permutation[column] == row else ZERO for column in range(dimension))
        for row in range(dimension)
    )
