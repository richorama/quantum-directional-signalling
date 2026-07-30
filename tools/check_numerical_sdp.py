#!/usr/bin/env python3
"""Independent numerical checks of the unreduced directional-signalling SDP.

This script does not use the covariance reductions or exact certificate
implementation in ``quantum_coarse_graining``.  It optimizes simultaneously
over an arbitrary CPTP effective channel and the dual diamond-norm SDP.
The resulting values are compared with the manuscript formulas.
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

import cvxpy as cp
import numpy as np
from scipy.optimize import brentq, minimize_scalar


@dataclass(frozen=True)
class Check:
    name: str
    unitary: np.ndarray
    dim_a: int
    dim_b: int
    expected: float
    tolerance: float = 8e-4


def matrix_unit(dimension: int, row: int, column: int) -> np.ndarray:
    value = np.zeros((dimension, dimension), dtype=complex)
    value[row, column] = 1
    return value


def partial_trace_b(
    operator: np.ndarray, dim_a: int, dim_b: int
) -> np.ndarray:
    reshaped = operator.reshape(dim_a, dim_b, dim_a, dim_b)
    return np.trace(reshaped, axis1=1, axis2=3)


def reduced_channel_choi(
    unitary: np.ndarray, dim_a: int, dim_b: int
) -> np.ndarray:
    input_dimension = dim_a * dim_b
    choi = np.zeros(
        (dim_a * input_dimension, dim_a * input_dimension),
        dtype=complex,
    )
    for row in range(input_dimension):
        for column in range(input_dimension):
            basis = matrix_unit(input_dimension, row, column)
            image = partial_trace_b(
                unitary @ basis @ unitary.conj().T,
                dim_a,
                dim_b,
            )
            for output_row in range(dim_a):
                for output_column in range(dim_a):
                    choi[
                        output_row * input_dimension + row,
                        output_column * input_dimension + column,
                    ] = image[output_row, output_column]
    return choi


def effective_after_trace_choi(
    channel_choi: cp.Expression, dim_a: int, dim_b: int
) -> cp.Expression:
    input_dimension = dim_a * dim_b
    rows = []
    for output_row in range(dim_a):
        for input_row in range(input_dimension):
            a_row, b_row = divmod(input_row, dim_b)
            entries = []
            for output_column in range(dim_a):
                for input_column in range(input_dimension):
                    a_column, b_column = divmod(input_column, dim_b)
                    if b_row == b_column:
                        entries.append(
                            channel_choi[
                                output_row * dim_a + a_row,
                                output_column * dim_a + a_column,
                            ]
                        )
                    else:
                        entries.append(0)
            rows.append(entries)
    return cp.bmat(rows)


def partial_trace_output(
    operator: cp.Expression, output_dimension: int, input_dimension: int
) -> cp.Expression:
    return cp.bmat(
        [
            [
                sum(
                    operator[
                        output * input_dimension + row,
                        output * input_dimension + column,
                    ]
                    for output in range(output_dimension)
                )
                for column in range(input_dimension)
            ]
            for row in range(input_dimension)
        ]
    )


def signalling_sdp(
    unitary: np.ndarray,
    dim_a: int,
    dim_b: int,
    solver: str = "SCS",
) -> float:
    """Solve ``min_E ||Tr_B Ad_U - E Tr_B||_diamond`` without symmetry."""
    input_dimension = dim_a * dim_b
    reduced_choi = reduced_channel_choi(unitary, dim_a, dim_b)

    channel_choi = cp.Variable(
        (dim_a * dim_a, dim_a * dim_a),
        hermitian=True,
        name="effective_channel_choi",
    )
    dual = cp.Variable(
        (dim_a * input_dimension, dim_a * input_dimension),
        hermitian=True,
        name="diamond_dual",
    )
    bound = cp.Variable(nonneg=True, name="diamond_bound")

    residual = reduced_choi - effective_after_trace_choi(
        channel_choi,
        dim_a,
        dim_b,
    )
    constraints = [
        channel_choi >> 0,
        partial_trace_output(channel_choi, dim_a, dim_a) == np.eye(dim_a),
        dual - residual >> 0,
        dual + residual >> 0,
        bound * np.eye(input_dimension)
        - partial_trace_output(dual, dim_a, input_dimension)
        >> 0,
    ]
    problem = cp.Problem(cp.Minimize(bound), constraints)
    options = {"eps": 2e-6, "max_iters": 100_000} if solver == "SCS" else {}
    value = problem.solve(solver=solver, verbose=False, **options)
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        raise RuntimeError("SDP failed with status {0}".format(problem.status))
    return float(value)


def pauli_x() -> np.ndarray:
    return np.array(((0, 1), (1, 0)), dtype=complex)


def pauli_y() -> np.ndarray:
    return np.array(((0, -1j), (1j, 0)), dtype=complex)


def pauli_z() -> np.ndarray:
    return np.array(((1, 0), (0, -1)), dtype=complex)


def ising_unitary(theta: float) -> np.ndarray:
    zz = np.kron(pauli_z(), pauli_z())
    return math.cos(theta) * np.eye(4) - 1j * math.sin(theta) * zz


def swap_unitary(dimension: int) -> np.ndarray:
    result = np.zeros((dimension * dimension,) * 2, dtype=complex)
    for a in range(dimension):
        for b in range(dimension):
            result[b * dimension + a, a * dimension + b] = 1
    return result


def partial_swap_unitary(dimension: int, phi: float) -> np.ndarray:
    dimension_squared = dimension * dimension
    return (
        math.cos(phi) * np.eye(dimension_squared)
        - 1j * math.sin(phi) * swap_unitary(dimension)
    )


def xy_unitary(theta: float) -> np.ndarray:
    generator = np.kron(pauli_x(), pauli_x()) + np.kron(pauli_y(), pauli_y())
    eigenvalues, eigenvectors = np.linalg.eigh(generator)
    return (
        eigenvectors
        @ np.diag(np.exp(-1j * theta * eigenvalues))
        @ eigenvectors.conj().T
    )


def partial_swap_fixed_defect(
    dimension: int, phi: float, shrinkage: float
) -> float:
    squared_dimension = dimension * dimension
    sine = math.sin(phi)
    a = 1 - shrinkage
    radical = a * a + squared_dimension * shrinkage * sine * sine
    if radical <= 1e-15:
        return 0.0
    h_zero = dimension * (squared_dimension - 3) / (squared_dimension - 1)
    h_one = 2 / (squared_dimension - 1)
    quadratic = a * a * (h_zero * h_zero - h_one * h_one) + 4 * radical
    correction = 1 - a * a * h_one * h_one / (4 * radical)
    interior = (a * h_zero + math.sqrt(quadratic)) / correction
    if a * h_one * interior <= 4 * radical:
        return interior / dimension
    return 2 * a * (h_zero + h_one) / dimension


def partial_swap_expected(dimension: int, phi: float) -> float:
    lower = -1 / (dimension * dimension - 1)
    result = minimize_scalar(
        lambda shrinkage: partial_swap_fixed_defect(
            dimension,
            phi,
            shrinkage,
        ),
        bounds=(lower, 1),
        method="bounded",
        options={"xatol": 1e-12},
    )
    if not result.success:
        raise RuntimeError("scalar partial-SWAP minimization failed")
    return float(result.fun)


def xy_quartic_coefficients(tangent: float) -> Tuple[float, ...]:
    t = tangent
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


def evaluate_polynomial(coefficients: Sequence[float], value: float) -> float:
    result = 0.0
    for coefficient in coefficients:
        result = result * value + coefficient
    return result


THETA_1 = 0.258270520262
THETA_2 = 0.646615513406


def xy_expected(theta: float) -> float:
    sine = math.sin(2 * theta)
    cosine = math.cos(2 * theta)
    if theta <= THETA_1:
        return 2 * sine * cosine
    if theta >= THETA_2:
        return sine + sine * sine / 2
    coefficients = xy_quartic_coefficients(math.tan(theta))
    return brentq(
        lambda defect: evaluate_polynomial(coefficients, defect),
        0.5,
        2,
    )


def checks(full: bool) -> Tuple[Check, ...]:
    cases = [
        Check(
            "Ising interior",
            ising_unitary(0.31),
            2,
            2,
            abs(math.sin(0.62)),
        ),
        Check(
            "partial-SWAP weak",
            partial_swap_unitary(2, math.asin(0.25)),
            2,
            2,
            0.5,
        ),
        Check(
            "partial-SWAP strong",
            partial_swap_unitary(2, math.asin(0.8)),
            2,
            2,
            partial_swap_expected(2, math.asin(0.8)),
        ),
        Check("XY weak", xy_unitary(0.18), 2, 2, xy_expected(0.18)),
        Check("XY middle", xy_unitary(0.45), 2, 2, xy_expected(0.45)),
        Check("XY strong", xy_unitary(0.70), 2, 2, xy_expected(0.70)),
        Check("iSWAP", xy_unitary(math.pi / 4), 2, 2, 1.5),
    ]
    if full:
        cases.extend(
            (
                Check(
                    "qutrit partial-SWAP weak",
                    partial_swap_unitary(3, math.asin(0.6)),
                    3,
                    3,
                    1.2,
                    tolerance=1.5e-3,
                ),
                Check(
                    "qutrit SWAP",
                    swap_unitary(3),
                    3,
                    3,
                    16 / 9,
                    tolerance=1.5e-3,
                ),
            )
        )
    return tuple(cases)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Check manuscript formulas against the unreduced SDP."
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="include the slower qutrit partial-SWAP checks",
    )
    parser.add_argument(
        "--solver",
        default="SCS",
        help="CVXPY solver name (default: SCS)",
    )
    arguments = parser.parse_args(argv)

    failed = 0
    for case in checks(arguments.full):
        observed = signalling_sdp(
            case.unitary,
            case.dim_a,
            case.dim_b,
            solver=arguments.solver,
        )
        error = abs(observed - case.expected)
        status = "PASS" if error <= case.tolerance else "FAIL"
        print(
            "{0:4s}  {1:28s} observed={2:.8f} expected={3:.8f} "
            "error={4:.2e}".format(
                status,
                case.name,
                observed,
                case.expected,
                error,
            )
        )
        failed += status == "FAIL"
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
