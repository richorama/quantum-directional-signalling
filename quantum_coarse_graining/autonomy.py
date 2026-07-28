"""The Frobenius closure defect is scaled linear operator entanglement."""

from __future__ import annotations

from fractions import Fraction

from .exact import (
    Gaussian,
    Matrix,
    ZERO,
    conjugate_action,
    dagger,
    frobenius_norm2,
    gaussian,
    matmul,
    matrix_unit,
    partial_trace_b,
    scale,
    subtract,
    superoperator,
)


def _trace_and_evolve_superoperators(
    unitary: Matrix, dim_a: int, dim_b: int
) -> tuple[Matrix, Matrix]:
    dimension = dim_a * dim_b
    traced = []
    evolved = []
    for row in range(dimension):
        for column in range(dimension):
            basis = matrix_unit(dimension, row, column)
            traced.append(partial_trace_b(basis, dim_a, dim_b))
            evolved.append(
                partial_trace_b(
                    conjugate_action(unitary, basis), dim_a, dim_b
                )
            )
    return superoperator(traced), superoperator(evolved)


def closure_defect(unitary: Matrix, dim_a: int, dim_b: int) -> Fraction:
    """Squared Frobenius residual of the best unconstrained linear coarse law."""
    traced, evolved = _trace_and_evolve_superoperators(unitary, dim_a, dim_b)
    traced_dagger = dagger(traced)
    best = scale(
        gaussian(Fraction(1, dim_b)),
        matmul(evolved, traced_dagger),
    )
    residual = subtract(evolved, matmul(best, traced))
    return frobenius_norm2(residual)


def operator_reduced_state(unitary: Matrix, dim_a: int, dim_b: int) -> Matrix:
    """Reduced state of the normalized vectorized operator ``|U>`` on ``AA'``."""
    dimension = dim_a * dim_b
    realigned = tuple(
        tuple(
            unitary[dim_b * a + b][dim_b * c + d]
            for b in range(dim_b)
            for d in range(dim_b)
        )
        for a in range(dim_a)
        for c in range(dim_a)
    )
    return scale(
        gaussian(Fraction(1, dimension)),
        matmul(realigned, dagger(realigned)),
    )


def operator_purity(unitary: Matrix, dim_a: int, dim_b: int) -> Fraction:
    reduced = operator_reduced_state(unitary, dim_a, dim_b)
    value = sum(
        (
            reduced[row][column] * reduced[column][row]
            for row in range(len(reduced))
            for column in range(len(reduced))
        ),
        ZERO,
    )
    if value.imag != 0:
        raise ValueError("operator purity was not real")
    return value.real


def linear_operator_entanglement(
    unitary: Matrix, dim_a: int, dim_b: int
) -> Fraction:
    return Fraction(1) - operator_purity(unitary, dim_a, dim_b)


def normalized_closure_defect(
    unitary: Matrix, dim_a: int, dim_b: int
) -> Fraction:
    return closure_defect(unitary, dim_a, dim_b) / (dim_a * dim_a * dim_b)


def operator_entanglement_identity_holds(
    unitary: Matrix, dim_a: int, dim_b: int
) -> bool:
    """Exact identity ``D_A/(d_A^2 d_B) = 1 - Tr(rho_A^2)``."""
    return normalized_closure_defect(
        unitary, dim_a, dim_b
    ) == linear_operator_entanglement(unitary, dim_a, dim_b)
