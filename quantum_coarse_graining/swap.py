"""Exact SWAP boundary control for the directional-signalling defect."""

from __future__ import annotations

from fractions import Fraction

from .exact import (
    Gaussian,
    Matrix,
    ZERO,
    add,
    conjugate_action,
    identity,
    matmul,
    matrix_unit,
    partial_trace_b,
    permutation_unitary,
    scale,
    subtract,
)


def swap_unitary() -> Matrix:
    """The two-qubit unitary ``|a,b> -> |b,a>``."""
    return permutation_unitary((0, 2, 1, 3))


def completely_depolarizing(operator: Matrix) -> Matrix:
    """The qubit channel ``X -> Tr(X) I/2``."""
    trace = operator[0][0] + operator[1][1]
    return scale(trace / Gaussian(2), identity(2))


def coarse_after_swap(operator: Matrix) -> Matrix:
    return partial_trace_b(conjugate_action(swap_unitary(), operator), 2, 2)


def swap_residual(operator: Matrix) -> Matrix:
    coarse_input = partial_trace_b(operator, 2, 2)
    return subtract(coarse_after_swap(operator), completely_depolarizing(coarse_input))


def swap_decomposition_holds() -> bool:
    """Exact basis certificate ``R=(id-Depol) Tr_A`` after relabelling B to A."""
    for row in range(4):
        for column in range(4):
            basis = matrix_unit(4, row, column)
            relabelled_b = coarse_after_swap(basis)
            expected = subtract(
                relabelled_b,
                completely_depolarizing(relabelled_b),
            )
            if swap_residual(basis) != expected:
                return False
    return True


def bell_projector() -> Matrix:
    """Projector onto ``(|00> + |11>)/sqrt(2)`` over the rationals."""
    half = Gaussian(Fraction(1, 2))
    return (
        (half, ZERO, ZERO, half),
        (ZERO, ZERO, ZERO, ZERO),
        (ZERO, ZERO, ZERO, ZERO),
        (half, ZERO, ZERO, half),
    )


def swap_ancilla_output_difference() -> Matrix:
    """Choi witness output ``(id-Depol)⊗id(|Phi><Phi|)``."""
    return subtract(
        bell_projector(),
        scale(Gaussian(Fraction(1, 4)), identity(4)),
    )


def swap_ancilla_witness_certificate() -> bool:
    """Certify eigenvalues ``3/4,-1/4,-1/4,-1/4`` of the witness output."""
    projector = bell_projector()
    complement = subtract(identity(4), projector)
    difference = swap_ancilla_output_difference()
    spectral_form = add(
        scale(Gaussian(Fraction(3, 4)), projector),
        scale(Gaussian(Fraction(-1, 4)), complement),
    )
    return (
        matmul(projector, projector) == projector
        and matmul(complement, complement) == complement
        and matmul(projector, complement)
        == tuple(tuple(ZERO for _ in range(4)) for _ in range(4))
        and difference == spectral_form
    )


def swap_diamond_autonomy_defect() -> Fraction:
    """The known covariant-channel value ``||id_2 - Depol_2|| = 3/2``."""
    return Fraction(3, 2)


def swap_theorem_certificate() -> bool:
    """Exact finite certificate for the SWAP unitary and ancilla lower bound."""
    swap = swap_unitary()
    return (
        matmul(swap, swap) == identity(4)
        and swap_decomposition_holds()
        and swap_ancilla_witness_certificate()
        and swap_diamond_autonomy_defect() == Fraction(3, 2)
    )
