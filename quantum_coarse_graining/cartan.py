"""Exact symmetry reduction for two-qubit Cartan interactions."""

from __future__ import annotations

from fractions import Fraction
from typing import Tuple

from .exact import (
    Gaussian,
    Matrix,
    ZERO,
    add,
    conjugate_action,
    identity,
    kron,
    matmul,
    matrix_unit,
    partial_trace_b,
    scale,
    subtract,
)
from .pauli import PAULIS, PAULI_I, PAULI_X, PAULI_Y, PAULI_Z

CartanFactor = Tuple[Fraction, Fraction]
CartanParameters = Tuple[CartanFactor, CartanFactor, CartanFactor]


def _validate(parameters: CartanParameters) -> None:
    for cosine, sine in parameters:
        if cosine * cosine + sine * sine != 1:
            raise ValueError("each Cartan factor must lie on the unit circle")


def cartan_unitary(parameters: CartanParameters) -> Matrix:
    """Build ``exp[-i(a XX+b YY+c ZZ)]`` from exact cosine-sine pairs."""
    _validate(parameters)
    result = identity(4)
    for (cosine, sine), pauli in zip(
        parameters,
        (PAULI_X, PAULI_Y, PAULI_Z),
    ):
        factor = add(
            scale(Gaussian(cosine), identity(4)),
            scale(Gaussian(0, -sine), kron(pauli, pauli)),
        )
        result = matmul(factor, result)
    return result


def cartan_pauli_coefficients(
    parameters: CartanParameters,
) -> Tuple[Gaussian, Gaussian, Gaussian, Gaussian]:
    """Coefficients in ``U=sum_mu a_mu sigma_mu x sigma_mu``."""
    _validate(parameters)
    (cx, sx), (cy, sy), (cz, sz) = parameters
    return (
        Gaussian(cx * cy * cz, -sx * sy * sz),
        Gaussian(cx * sy * sz, -sx * cy * cz),
        Gaussian(sx * cy * sz, -cx * sy * cz),
        Gaussian(sx * sy * cz, -cx * cy * sz),
    )


def cartan_pauli_weights(parameters: CartanParameters) -> Tuple[Fraction, ...]:
    """Probabilities of the maximally-mixed-environment Pauli channel."""
    return tuple(
        coefficient.norm2()
        for coefficient in cartan_pauli_coefficients(parameters)
    )


def cartan_effective_channel(
    operator: Matrix,
    parameters: CartanParameters,
) -> Matrix:
    """The admissible channel ``sum_mu |a_mu|^2 Ad_sigma_mu``."""
    result = tuple(tuple(ZERO for _ in range(2)) for _ in range(2))
    for weight, pauli in zip(cartan_pauli_weights(parameters), PAULIS):
        image = matmul(matmul(pauli, operator), pauli)
        result = add(result, scale(Gaussian(weight), image))
    return result


def cartan_effective_channel_is_cptp(parameters: CartanParameters) -> bool:
    """Random-unitary certificate: all Pauli weights are nonnegative and sum to 1."""
    weights = cartan_pauli_weights(parameters)
    return all(weight >= 0 for weight in weights) and sum(weights) == 1


def coarse_after_cartan(
    operator: Matrix,
    parameters: CartanParameters,
) -> Matrix:
    return partial_trace_b(
        conjugate_action(cartan_unitary(parameters), operator),
        2,
        2,
    )


def cartan_residual(
    operator: Matrix,
    parameters: CartanParameters,
) -> Matrix:
    return subtract(
        coarse_after_cartan(operator, parameters),
        cartan_effective_channel(partial_trace_b(operator, 2, 2), parameters),
    )


def cartan_visible_sector_certificate(parameters: CartanParameters) -> bool:
    """The candidate channel is exact when B is maximally mixed."""
    half_identity = scale(Gaussian(Fraction(1, 2)), PAULI_I)
    for row in range(2):
        for column in range(2):
            operator = matrix_unit(2, row, column)
            extension = kron(operator, half_identity)
            if coarse_after_cartan(
                extension,
                parameters,
            ) != cartan_effective_channel(operator, parameters):
                return False
            if cartan_residual(extension, parameters) != (
                (ZERO, ZERO),
                (ZERO, ZERO),
            ):
                return False
    return True


def cartan_joint_pauli_covariance_certificate(
    parameters: CartanParameters,
) -> bool:
    """Certify ``N Ad_(P x P)=Ad_P N`` and the matching trace covariance."""
    for pauli in PAULIS:
        joint = kron(pauli, pauli)
        for row in range(4):
            for column in range(4):
                operator = matrix_unit(4, row, column)
                transformed = conjugate_action(joint, operator)
                evolved_left = coarse_after_cartan(transformed, parameters)
                evolved_right = conjugate_action(
                    pauli,
                    coarse_after_cartan(operator, parameters),
                )
                traced_left = partial_trace_b(transformed, 2, 2)
                traced_right = conjugate_action(
                    pauli,
                    partial_trace_b(operator, 2, 2),
                )
                if evolved_left != evolved_right or traced_left != traced_right:
                    return False
    return True


def cartan_symmetry_certificate(parameters: CartanParameters) -> bool:
    """Finite exact certificate for the ingredients of the Pauli-twirling lemma."""
    return (
        cartan_effective_channel_is_cptp(parameters)
        and cartan_visible_sector_certificate(parameters)
        and cartan_joint_pauli_covariance_certificate(parameters)
    )
