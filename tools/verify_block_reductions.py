#!/usr/bin/env python3
"""Verify the physical Choi block reductions over exact arithmetic.

This script constructs the unreduced residual Choi operators directly, scales
them by invariant input densities, and compares their characteristic
polynomials with the spectra stated in the manuscript.  It does not import the
certificate package or use the closed-form diamond norms.
"""

from __future__ import annotations

import sys
from typing import Callable, Sequence, Tuple

import sympy as sp


def matrix_unit(dimension: int, row: int, column: int) -> sp.Matrix:
    value = sp.zeros(dimension)
    value[row, column] = 1
    return value


def swap(dimension: int) -> sp.Matrix:
    value = sp.zeros(dimension * dimension)
    for left in range(dimension):
        for right in range(dimension):
            value[right * dimension + left, left * dimension + right] = 1
    return value


def partial_trace_b(
    operator: sp.Matrix, dim_a: int, dim_b: int
) -> sp.Matrix:
    return sp.Matrix(
        dim_a,
        dim_a,
        lambda row, column: sum(
            operator[row * dim_b + hidden, column * dim_b + hidden]
            for hidden in range(dim_b)
        ),
    )


def choi(
    image: Callable[[sp.Matrix], sp.Matrix],
    input_dimension: int,
    output_dimension: int,
) -> sp.Matrix:
    value = sp.zeros(output_dimension * input_dimension)
    for row in range(input_dimension):
        for column in range(input_dimension):
            block = image(matrix_unit(input_dimension, row, column))
            for output_row in range(output_dimension):
                for output_column in range(output_dimension):
                    value[
                        output_row * input_dimension + row,
                        output_column * input_dimension + column,
                    ] = block[output_row, output_column]
    return value


def depolarizing(operator: sp.Matrix, shrinkage: sp.Rational) -> sp.Matrix:
    dimension = operator.rows
    return (
        shrinkage * operator
        + (1 - shrinkage) * sp.trace(operator) * sp.eye(dimension) / dimension
    )


def partial_swap_scaled_choi(
    dimension: int,
    cosine: sp.Rational,
    sine: sp.Rational,
    shrinkage: sp.Rational,
) -> sp.Matrix:
    unitary = cosine * sp.eye(dimension**2) - sp.I * sine * swap(dimension)

    def residual(operator: sp.Matrix) -> sp.Matrix:
        evolved = partial_trace_b(
            unitary * operator * unitary.conjugate().T,
            dimension,
            dimension,
        )
        visible = partial_trace_b(operator, dimension, dimension)
        return sp.simplify(evolved - depolarizing(visible, shrinkage))

    # rho=I/d^2, so sqrt(rho)=I/d and K_rho=J/d^2.
    return choi(residual, dimension**2, dimension) / dimension**2


def partial_swap_predicted_polynomial(
    variable: sp.Symbol,
    dimension: int,
    cosine: sp.Rational,
    sine: sp.Rational,
    shrinkage: sp.Rational,
) -> sp.Expr:
    del cosine
    weight = sp.Rational(1, dimension**2)
    symmetric_weight = sp.Rational(dimension + 1, 2 * dimension)
    antisymmetric_weight = 1 - symmetric_weight
    p_term = sp.Rational((dimension - 1) * (dimension + 2), dimension + 1)
    r_term = sp.Rational((dimension + 1) * (dimension - 2), dimension - 1)
    radical = (1 + shrinkage) ** 2 - 4 * shrinkage * (1 - sine**2)
    root_sum = (
        (1 - shrinkage)
        * (
            p_term * symmetric_weight
            + r_term * antisymmetric_weight
        )
        / dimension**2
    )
    root_product = (
        symmetric_weight
        * antisymmetric_weight
        / dimension**2
        * (
            (1 - shrinkage) ** 2 * p_term * r_term / dimension**2
            - radical
        )
    )
    scalar = -(1 - shrinkage) * weight / dimension
    scalar_multiplicity = dimension * (dimension**2 - 2)
    return sp.expand(
        (variable - scalar) ** scalar_multiplicity
        * (variable**2 - root_sum * variable + root_product) ** dimension
    )


def pauli_channel(
    operator: sp.Matrix,
    transverse: sp.Rational,
    longitudinal: sp.Rational,
) -> sp.Matrix:
    identity = sp.eye(2)
    x = sp.Matrix(((0, 1), (1, 0)))
    y = sp.Matrix(((0, -sp.I), (sp.I, 0)))
    z = sp.diag(1, -1)
    probabilities = (
        (1 + 2 * transverse + longitudinal) / 4,
        (1 - longitudinal) / 4,
        (1 - longitudinal) / 4,
        (1 + longitudinal - 2 * transverse) / 4,
    )
    return sp.simplify(
        sum(
            (
                probability * pauli * operator * pauli
                for probability, pauli in zip(
                    probabilities,
                    (identity, x, y, z),
                )
            ),
            sp.zeros(2),
        )
    )


def xy_scaled_choi(
    cosine: sp.Rational,
    sine: sp.Rational,
    transverse: sp.Rational,
    longitudinal: sp.Rational,
) -> sp.Matrix:
    unitary = sp.Matrix(
        (
            (1, 0, 0, 0),
            (0, cosine, -sp.I * sine, 0),
            (0, -sp.I * sine, cosine, 0),
            (0, 0, 0, 1),
        )
    )

    def residual(operator: sp.Matrix) -> sp.Matrix:
        evolved = partial_trace_b(unitary * operator * unitary.conjugate().T, 2, 2)
        visible = partial_trace_b(operator, 2, 2)
        return sp.simplify(
            evolved - pauli_channel(visible, transverse, longitudinal)
        )

    # A=B=1/4, so sqrt(rho)=I/2 and K_rho=J/4.
    return choi(residual, 4, 2) / 4


def xy_predicted_polynomial(
    variable: sp.Symbol,
    cosine: sp.Rational,
    sine: sp.Rational,
    transverse: sp.Rational,
    longitudinal: sp.Rational,
) -> sp.Expr:
    even = sp.Rational(1, 4)
    odd = sp.Rational(1, 4)
    alpha = even * (1 - longitudinal) / 2
    gamma_squared = even * odd / 2
    offset = cosine - transverse
    beta_real = odd * (
        longitudinal - (cosine**2 - sine**2)
    ) / 2
    beta_imag = odd * sine * cosine
    beta_norm = beta_real**2 + beta_imag**2
    p_two = -alpha
    p_one = -(
        2 * gamma_squared * (offset**2 + sine**2) + beta_norm
    )
    p_zero = (
        alpha * beta_norm
        + 2 * gamma_squared * (offset**2 - sine**2) * beta_real
        - 4 * gamma_squared * offset * sine * beta_imag
    )
    cubic = variable**3 + p_two * variable**2 + p_one * variable + p_zero
    return sp.expand((variable + alpha) ** 2 * cubic**2)


def xy_axial_charge() -> sp.Matrix:
    z = sp.diag(1, -1)
    identity = sp.eye(2)
    return (
        sp.kronecker_product(z, identity, identity)
        - sp.kronecker_product(identity, z, identity)
        - sp.kronecker_product(identity, identity, z)
    ) / 2


def check(name: str, actual: sp.Expr, expected: sp.Expr) -> bool:
    difference = actual - expected
    if isinstance(difference, sp.MatrixBase):
        holds = difference.applyfunc(sp.expand).is_zero_matrix
    else:
        holds = sp.expand(difference) == 0
    print("{0:4s}  {1}".format("PASS" if holds else "FAIL", name))
    return holds


def qubit_stationarity_factorization() -> Tuple[sp.Expr, sp.Expr, sp.Expr]:
    shrinkage, sine = sp.symbols("lambda s")
    radical = (
        (1 + shrinkage) ** 2
        - 4 * shrinkage * (1 - sine**2)
    )
    offset = shrinkage - 1 + 2 * sine**2
    squared_equation = (
        (2 * offset * (1 - shrinkage) + radical) ** 2
        - 9 * offset**2 * radical
    )
    quartic = (
        2 * shrinkage**4
        + (18 * sine**2 - 8) * shrinkage**3
        + (45 * sine**4 - 43 * sine**2 + 12) * shrinkage**2
        + (
            36 * sine**6
            - 54 * sine**4
            + 32 * sine**2
            - 8
        )
        * shrinkage
        + 5 * sine**4
        - 7 * sine**2
        + 2
    )
    endpoint_zero = sp.expand(
        quartic.subs(shrinkage, 0)
        - (sine**2 - 1) * (5 * sine**2 - 2)
    )
    endpoint_one = sp.expand(
        quartic.subs(shrinkage, 1)
        - 4 * sine**4 * (9 * sine**2 - 1)
    )
    endpoint_residual = sp.expand(endpoint_zero**2 + endpoint_one**2)
    return (
        sp.expand(squared_equation),
        sp.expand(-4 * quartic),
        endpoint_residual,
    )


def main(argv: Sequence[str] | None = None) -> int:
    del argv
    variable = sp.Symbol("z")
    failed = 0
    for dimension in (2, 3, 4):
        scaled = partial_swap_scaled_choi(
            dimension,
            sp.Rational(3, 5),
            sp.Rational(4, 5),
            sp.Rational(1, 2),
        )
        actual = scaled.charpoly(variable).as_expr()
        expected = partial_swap_predicted_polynomial(
            variable,
            dimension,
            sp.Rational(3, 5),
            sp.Rational(4, 5),
            sp.Rational(1, 2),
        )
        failed += not check(
            "partial-SWAP physical Choi spectrum, d={0}".format(dimension),
            actual,
            expected,
        )

    scaled_xy = xy_scaled_choi(
        sp.Rational(3, 5),
        sp.Rational(4, 5),
        sp.Rational(1, 3),
        sp.Rational(1, 2),
    )
    failed += not check(
        "XY physical Choi spectrum reduces to two scalar-plus-M3 blocks",
        scaled_xy.charpoly(variable).as_expr(),
        xy_predicted_polynomial(
            variable,
            sp.Rational(3, 5),
            sp.Rational(4, 5),
            sp.Rational(1, 3),
            sp.Rational(1, 2),
        ),
    )
    failed += not check(
        "XY scaled Choi operator commutes with the stated axial charge",
        scaled_xy * xy_axial_charge() - xy_axial_charge() * scaled_xy,
        sp.zeros(8),
    )
    actual, expected, endpoint_residual = qubit_stationarity_factorization()
    failed += not check(
        "qubit strong partial-SWAP stationarity reduces to the stated quartic",
        actual,
        expected,
    )
    failed += not check(
        "qubit stationarity quartic has the stated endpoint factorizations",
        endpoint_residual,
        sp.Integer(0),
    )
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
