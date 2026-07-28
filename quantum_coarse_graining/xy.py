"""Exact certificates for the XY/iSWAP directional-signalling family."""

from __future__ import annotations

from fractions import Fraction
from typing import Tuple

from .cartan import CartanFactor, cartan_unitary
from .exact import Gaussian, Matrix, add, identity, matmul, scale
from .pauli import PAULIS

PauliProbabilities = Tuple[Fraction, Fraction, Fraction, Fraction]
QuarticCoefficients = Tuple[Fraction, Fraction, Fraction, Fraction, Fraction]


def _validate_unit_circle(cosine: Fraction, sine: Fraction) -> None:
    if cosine * cosine + sine * sine != 1:
        raise ValueError("cosine and sine must lie on the unit circle")


def xy_unitary(factor: CartanFactor) -> Matrix:
    """Build ``exp[-i theta (XX+YY)]`` exactly."""
    cosine, sine = (Fraction(value) for value in factor)
    _validate_unit_circle(cosine, sine)
    return cartan_unitary(
        (
            (cosine, sine),
            (cosine, sine),
            (Fraction(1), Fraction(0)),
        )
    )


def xy_transfer_parameters(factor: CartanFactor) -> Tuple[Fraction, Fraction]:
    """Return ``(cos(2 theta), sin(2 theta))``."""
    cosine, sine = (Fraction(value) for value in factor)
    _validate_unit_circle(cosine, sine)
    return cosine * cosine - sine * sine, 2 * cosine * sine


def xy_weak_threshold_polynomial(
    transfer_cosine: Fraction,
    transfer_sine: Fraction,
) -> Fraction:
    """Signed weak-threshold equation, nonpositive on the weak branch."""
    cosine, sine = Fraction(transfer_cosine), Fraction(transfer_sine)
    _validate_unit_circle(cosine, sine)
    return (
        2 * sine * cosine
        - (3 * cosine**3 + cosine * cosine - cosine - 1)
    )


def xy_strong_threshold_polynomial(
    transfer_cosine: Fraction,
    transfer_sine: Fraction,
) -> Fraction:
    """The final-channel ``q_Z`` numerator, nonnegative in regime III."""
    cosine, sine = Fraction(transfer_cosine), Fraction(transfer_sine)
    _validate_unit_circle(cosine, sine)
    return 1 + cosine * cosine - 2 * cosine * (1 + sine)


def xy_first_threshold_tangent_polynomial(tangent: Fraction) -> Fraction:
    """Polynomial whose first positive root is the weak/middle threshold."""
    tangent = Fraction(tangent)
    return (
        tangent**6
        - 2 * tangent**5
        - 3 * tangent**4
        + 7 * tangent**2
        + 2 * tangent
        - 1
    )


def xy_second_threshold_tangent_polynomial(tangent: Fraction) -> Fraction:
    """Polynomial whose positive root is the middle/strong threshold."""
    tangent = Fraction(tangent)
    return tangent**3 + tangent**2 - 1


def xy_middle_quartic_coefficients(
    tangent: Fraction,
) -> QuarticCoefficients:
    """Return the conditional middle-branch quartic coefficients.

    If the numerically observed ``q_Z=0`` facet is globally optimal, the
    intermediate defect is the unique physical root of
    ``sum(coefficients[k] * delta**(4-k), k=0..4)``.
    """
    t = Fraction(tangent)
    return (
        t**6 + t**5 + 3 * t**4 + 2 * t**3 + 3 * t**2 + t + 1,
        (
            2 * t**7
            - 46 * t**6
            - 122 * t**5
            - 170 * t**4
            - 138 * t**3
            - 58 * t**2
            - 14 * t
            + 2
        ),
        (
            -84 * t**7
            + 24 * t**6
            + 72 * t**5
            + 96 * t**4
            - 68 * t**3
            - 56 * t**2
            - 32 * t
        ),
        (
            96 * t**7
            + 64 * t**6
            + 160 * t**5
            + 384 * t**4
            + 320 * t**3
            + 128 * t**2
        ),
        64 * t**5 * (2 * t**2 + 2 * t + 1),
    )


def xy_middle_quartic(
    defect: Fraction,
    tangent: Fraction,
) -> Fraction:
    """Evaluate the conditional middle-branch quartic exactly."""
    defect = Fraction(defect)
    coefficients = xy_middle_quartic_coefficients(tangent)
    return sum(
        coefficient * defect ** (4 - index)
        for index, coefficient in enumerate(coefficients)
    )


def xy_weak_probabilities(
    transfer_cosine: Fraction,
    transfer_sine: Fraction,
) -> PauliProbabilities:
    """An optimal Pauli channel on the exact weak branch."""
    cosine, sine = Fraction(transfer_cosine), Fraction(transfer_sine)
    _validate_unit_circle(cosine, sine)
    if (
        cosine < 0
        or sine < 0
        or xy_weak_threshold_polynomial(cosine, sine) > 0
    ):
        raise ValueError("outside the weak XY branch")
    return (
        cosine * cosine,
        sine * sine / 2,
        sine * sine / 2,
        Fraction(0),
    )


def xy_weak_defect(
    transfer_cosine: Fraction,
    transfer_sine: Fraction,
) -> Fraction:
    """Sharp weak-branch value ``sin(4 theta)``."""
    xy_weak_probabilities(transfer_cosine, transfer_sine)
    return 2 * Fraction(transfer_cosine) * Fraction(transfer_sine)


def xy_strong_probabilities(
    transfer_cosine: Fraction,
    transfer_sine: Fraction,
) -> PauliProbabilities:
    """The optimal Pauli channel in the central iSWAP regime."""
    cosine, sine = Fraction(transfer_cosine), Fraction(transfer_sine)
    _validate_unit_circle(cosine, sine)
    if (
        cosine < 0
        or sine < cosine
        or xy_strong_threshold_polynomial(cosine, sine) < 0
    ):
        raise ValueError("outside the strong XY branch")
    common = sine * sine / 4
    coherence = 2 * cosine * (1 + sine)
    return (
        (1 + cosine * cosine + coherence) / 4,
        common,
        common,
        (1 + cosine * cosine - coherence) / 4,
    )


def xy_strong_defect(
    transfer_cosine: Fraction,
    transfer_sine: Fraction,
) -> Fraction:
    """Sharp central-branch value ``sin(2 theta)+sin(2 theta)^2/2``."""
    xy_strong_probabilities(transfer_cosine, transfer_sine)
    sine = Fraction(transfer_sine)
    return sine + sine * sine / 2


def xy_pauli_channel(
    operator: Matrix,
    probabilities: PauliProbabilities,
) -> Matrix:
    """Apply a Pauli channel after validating its probability simplex."""
    probabilities = tuple(Fraction(value) for value in probabilities)
    if len(probabilities) != 4:
        raise ValueError("four Pauli probabilities are required")
    if any(value < 0 for value in probabilities) or sum(probabilities) != 1:
        raise ValueError("Pauli probabilities must be nonnegative and sum to one")
    result = scale(Gaussian(0), identity(2))
    for probability, pauli in zip(probabilities, PAULIS):
        result = add(
            result,
            scale(
                Gaussian(probability),
                matmul(matmul(pauli, operator), pauli),
            ),
        )
    return result


def xy_block_characteristic(
    transfer_cosine: Fraction,
    transfer_sine: Fraction,
    transverse_eigenvalue: Fraction,
    longitudinal_eigenvalue: Fraction,
    even_weight: Fraction,
) -> Tuple[Fraction, Fraction, Fraction, Fraction]:
    """Characteristic coefficients of the invariant ``M_3`` Choi block.

    ``even_weight`` is the parameter ``A`` in ``[0, 1/2]``.  The returned
    tuple is ``(1, p2, p1, p0)`` for
    ``x^3+p2*x^2+p1*x+p0``.
    """
    cosine = Fraction(transfer_cosine)
    sine = Fraction(transfer_sine)
    transverse = Fraction(transverse_eigenvalue)
    longitudinal = Fraction(longitudinal_eigenvalue)
    even = Fraction(even_weight)
    _validate_unit_circle(cosine, sine)
    if not 0 <= even <= Fraction(1, 2):
        raise ValueError("even_weight must lie in [0, 1/2]")
    odd = Fraction(1, 2) - even
    alpha = even * (1 - longitudinal) / 2
    gamma_squared = even * odd / 2
    offset = cosine - transverse
    beta_real = odd * (longitudinal - (cosine * cosine - sine * sine)) / 2
    beta_imag = odd * sine * cosine
    beta_norm = beta_real * beta_real + beta_imag * beta_imag
    p_two = -alpha
    p_one = -(
        2 * gamma_squared * (offset * offset + sine * sine)
        + beta_norm
    )
    p_zero = (
        alpha * beta_norm
        + 2 * gamma_squared * (offset * offset - sine * sine) * beta_real
        - 4 * gamma_squared * offset * sine * beta_imag
    )
    return Fraction(1), p_two, p_one, p_zero


def xy_strong_spectrum_certificate(
    transfer_cosine: Fraction,
    transfer_sine: Fraction,
) -> bool:
    """Certify the rational top root and signs at the maximally mixed witness."""
    cosine, sine = Fraction(transfer_cosine), Fraction(transfer_sine)
    probabilities = xy_strong_probabilities(cosine, sine)
    transverse = probabilities[0] - probabilities[3]
    longitudinal = (
        probabilities[0]
        - probabilities[1]
        - probabilities[2]
        + probabilities[3]
    )
    _, p_two, p_one, p_zero = xy_block_characteristic(
        cosine,
        sine,
        transverse,
        longitudinal,
        Fraction(1, 4),
    )
    top = xy_strong_defect(cosine, sine) / 4
    polynomial = top**3 + p_two * top**2 + p_one * top + p_zero
    other_sum = -p_two - top
    other_product = -p_zero / top
    expected_product = sine * sine * (7 * sine * sine + 2 * sine - 8) / 64
    return (
        polynomial == 0
        and other_sum == -sine / 4
        and other_product == expected_product
        and other_sum < 0
        and other_product > 0
    )


def xy_boundary_certificate() -> bool:
    """Exact weak, strong, and iSWAP controls."""
    weak_cosine, weak_sine = Fraction(12, 13), Fraction(5, 13)
    strong_cosine, strong_sine = Fraction(9, 41), Fraction(40, 41)
    return (
        xy_middle_quartic_coefficients(Fraction(1, 2))
        == (
            Fraction(175, 64),
            Fraction(-3321, 64),
            Fraction(-977, 32),
            Fraction(411, 4),
            Fraction(5),
        )
        and xy_first_threshold_tangent_polynomial(Fraction(1, 2))
        == Fraction(97, 64)
        and xy_second_threshold_tangent_polynomial(Fraction(1, 2))
        == Fraction(-5, 8)
        and xy_weak_defect(weak_cosine, weak_sine) == Fraction(120, 169)
        and xy_weak_probabilities(weak_cosine, weak_sine)
        == (
            Fraction(144, 169),
            Fraction(25, 338),
            Fraction(25, 338),
            Fraction(0),
        )
        and xy_strong_defect(strong_cosine, strong_sine)
        == Fraction(2440, 1681)
        and xy_strong_probabilities(strong_cosine, strong_sine)
        == (
            Fraction(805, 1681),
            Fraction(400, 1681),
            Fraction(400, 1681),
            Fraction(76, 1681),
        )
        and xy_strong_spectrum_certificate(strong_cosine, strong_sine)
        and xy_strong_defect(Fraction(0), Fraction(1)) == Fraction(3, 2)
    )
