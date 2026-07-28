"""Sharp directional-signalling theorem for a two-qubit Ising interaction."""

from __future__ import annotations

from fractions import Fraction

from .exact import (
    Gaussian,
    Matrix,
    ONE,
    ZERO,
    add,
    conjugate_action,
    identity,
    dagger,
    kron,
    matmul,
    matrix_unit,
    partial_trace_b,
    scale,
    subtract,
)
from .pauli import PAULI_X as _X
from .pauli import PAULI_Y as _Y
from .pauli import PAULI_Z as _Z


def ising_unitary(cosine: Fraction, sine: Fraction) -> Matrix:
    """``exp(-i theta ZxZ) = cos(theta) I - i sin(theta) ZxZ`` over ``Q(i)``."""
    if cosine * cosine + sine * sine != 1:
        raise ValueError("cosine and sine must lie on the unit circle")
    return add(
        scale(Gaussian(cosine), identity(4)),
        scale(Gaussian(0, -sine), kron(_Z, _Z)),
    )


def effective_dephasing(operator: Matrix, cosine: Fraction, sine: Fraction) -> Matrix:
    """The candidate optimal channel ``c^2 id + s^2 Ad_Z``."""
    return add(
        scale(Gaussian(cosine * cosine), operator),
        scale(
            Gaussian(sine * sine),
            matmul(matmul(_Z, operator), _Z),
        ),
    )


def coarse_after_ising(
    operator: Matrix, cosine: Fraction, sine: Fraction
) -> Matrix:
    return partial_trace_b(
        conjugate_action(ising_unitary(cosine, sine), operator), 2, 2
    )


def leakage_map(operator: Matrix) -> Matrix:
    """Unit-diamond-norm map ``-i[Z, Tr_B(Z_B X)]/2``."""
    weighted = matmul(kron(identity(2), _Z), operator)
    moment = partial_trace_b(weighted, 2, 2)
    commutator = subtract(matmul(_Z, moment), matmul(moment, _Z))
    return scale(Gaussian(0, Fraction(-1, 2)), commutator)


def ising_residual(
    operator: Matrix, cosine: Fraction, sine: Fraction
) -> Matrix:
    coarse_input = partial_trace_b(operator, 2, 2)
    return subtract(
        coarse_after_ising(operator, cosine, sine),
        effective_dephasing(coarse_input, cosine, sine),
    )


def ising_decomposition_holds(cosine: Fraction, sine: Fraction) -> bool:
    """Exact basis certificate ``N - E*T = 2cs L``."""
    coefficient = Gaussian(2 * cosine * sine)
    return all(
        ising_residual(matrix_unit(4, row, column), cosine, sine)
        == scale(coefficient, leakage_map(matrix_unit(4, row, column)))
        for row in range(4)
        for column in range(4)
    )


def hidden_correlation_witness() -> Matrix:
    """Trace-norm-one signed operator ``Y_A x Z_B / 4`` with zero A marginal."""
    return scale(Gaussian(Fraction(1, 4)), kron(_Y, _Z))


def product_state_lower_bound_certificate(
    cosine: Fraction, sine: Fraction
) -> bool:
    """Two physical inputs with the same A marginal attain the lower bound."""
    plus_state = scale(
        Gaussian(Fraction(1, 2)),
        add(identity(2), _X),
    )
    b_zero = matrix_unit(2, 0, 0)
    b_one = matrix_unit(2, 1, 1)
    input_zero = kron(plus_state, b_zero)
    input_one = kron(plus_state, b_one)
    output_difference = subtract(
        coarse_after_ising(input_zero, cosine, sine),
        coarse_after_ising(input_one, cosine, sine),
    )
    expected = scale(Gaussian(2 * cosine * sine), _Y)
    return (
        partial_trace_b(input_zero, 2, 2) == plus_state
        and partial_trace_b(input_one, 2, 2) == plus_state
        and output_difference == expected
        and matmul(output_difference, output_difference)
        == scale(
            Gaussian(4 * cosine * cosine * sine * sine),
            identity(2),
        )
    )


def hidden_correlation_witness_holds(
    cosine: Fraction, sine: Fraction
) -> bool:
    """The witness has zero coarse input and output ``-2cs X/2``.

    Since ``X/2`` has trace norm one, this supplies the matching diamond-norm
    induced-norm lower bound ``2|cs|``. A physical-state realization uses a
    classical flag for the positive and negative parts.
    """
    witness = hidden_correlation_witness()
    expected = scale(Gaussian(-cosine * sine), _X)
    return (
        partial_trace_b(witness, 2, 2) == ((ZERO, ZERO), (ZERO, ZERO))
        and coarse_after_ising(witness, cosine, sine) == expected
    )


def witness_trace_norm_certificate(
    cosine: Fraction, sine: Fraction
) -> bool:
    """Certify the trace norms used by the lower-bound witness.

    The Hermitian input obeys ``W^2=I_4/16``, hence has four singular values
    ``1/4`` and trace norm one. Its output is ``-cs X``, whose square is
    ``c^2 s^2 I_2`` and whose trace norm is ``2|cs|``.
    """
    witness = hidden_correlation_witness()
    output = coarse_after_ising(witness, cosine, sine)
    return (
        dagger(witness) == witness
        and matmul(witness, witness)
        == scale(Gaussian(Fraction(1, 16)), identity(4))
        and dagger(output) == output
        and matmul(output, output)
        == scale(Gaussian(cosine * cosine * sine * sine), identity(2))
    )


def flagged_state_witness_certificate(
    cosine: Fraction, sine: Fraction
) -> bool:
    """Certify the positive/negative parts used by the physical flagged state."""
    witness = hidden_correlation_witness()
    plus_projector = add(
        scale(Gaussian(Fraction(1, 2)), identity(4)),
        scale(Gaussian(2), witness),
    )
    minus_projector = subtract(identity(4), plus_projector)
    positive_part = scale(Gaussian(Fraction(1, 4)), plus_projector)
    negative_part = scale(Gaussian(Fraction(1, 4)), minus_projector)
    zero = tuple(tuple(ZERO for _ in range(4)) for _ in range(4))
    return (
        matmul(plus_projector, plus_projector) == plus_projector
        and matmul(minus_projector, minus_projector) == minus_projector
        and matmul(plus_projector, minus_projector) == zero
        and subtract(positive_part, negative_part) == witness
        and sum((positive_part[k][k] for k in range(4)), ZERO).real
        == Fraction(1, 2)
        and sum((negative_part[k][k] for k in range(4)), ZERO).real
        == Fraction(1, 2)
        and subtract(
            ising_residual(positive_part, cosine, sine),
            ising_residual(negative_part, cosine, sine),
        )
        == ising_residual(witness, cosine, sine)
    )


def diamond_autonomy_defect(cosine: Fraction, sine: Fraction) -> Fraction:
    """The sharp autonomy defect, equivalently reverse directional signalling."""
    return abs(2 * cosine * sine)


def ising_theorem_certificate(cosine: Fraction, sine: Fraction) -> bool:
    """Finite exact certificate for the analytic upper and lower bounds."""
    return (
        ising_decomposition_holds(cosine, sine)
        and product_state_lower_bound_certificate(cosine, sine)
        and hidden_correlation_witness_holds(cosine, sine)
        and witness_trace_norm_certificate(cosine, sine)
        and flagged_state_witness_certificate(cosine, sine)
    )
