#!/usr/bin/env python3
"""Exact SymPy re-derivation of the unconditional intermediate XY theorem.

This is optional tooling.  The core certificate suite in ``tests/`` stays
dependency free; only this script needs SymPy::

    python3 -m pip install "quantum-directional-signalling[proof]"
    python3 tools/derive_xy_middle.py

Every load-bearing exact identity behind the intermediate-branch theorem of
``paper/aqc_signalling.tex`` (and the ``AQC6`` regime II research notes) is
reconstructed from the invariant Choi block ``M_3`` and asserted:

1.  the ``M_3`` characteristic polynomial, its ``q_Z = 0`` facet form
    ``P(lambda;A,u) = lambda^3 - A(1-u)lambda^2 - (ABm+B^2N)lambda - AB^2Theta``
    modulo ``S^2 = 1-C^2``, and the determinant/inertia coefficient identities;
2.  the normal derivative ``P_u = 2 A B G`` taken at fixed ``w``, the tangential
    facet relation ``P_u + 2 P_w``, and the sign of ``G`` at the limiting weak
    saddle;
3.  the zero-multiplier elimination, by a resultant chain *and* by a saturated
    Groebner basis, leaving ``5C^3+C^2+3C-1`` as the only factor that can vanish
    in the open intermediate interval;
4.  the ``t`` parametrization ``F = 64(1+t^2)^2 P(Delta/4; a/2, u)``, its
    ``disc_u F = 4 W_+(a,Delta,t) W_+(a,Delta,-t)`` splitting, the second
    stationary polynomial ``G_8`` built from coefficient derivatives, and
    ``Res_a(W_+,G_8) = -2^36 Delta^12 t^28 (1+t^2)^10 E_1 E_2^2 Q_+``;
5.  the sheet selection ``W_+(a,Delta,-t) < 0``, the quartic discriminant
    factorization, the endpoint substitutions ``Delta = 1/2`` and ``Delta = 2``,
    and the weak/strong threshold joining identities;
6.  an exact Sturm root count for ``Q_+(Delta,1/2)`` on ``(1/2,2)``.

Positivity is certified exactly, never numerically.  Sign-definite quadratics
use their discriminant; box claims use rational Bernstein coefficients with
bisection (the Bernstein basis is nonnegative and sums to one, so a uniform
coefficient sign certifies that sign on the box).

Nothing about the middle branch is transcribed by hand: ``Q_+`` is rebuilt from
``quantum_coarse_graining.xy`` by exact rational interpolation and is then
required to match the eliminated resultant factor, so this script doubles as a
validation of the dependency-free backend.

Each assertion prints one ``PASS`` line.  Any failure prints ``FAIL`` and the
process exits with status 1.
"""

from __future__ import annotations

import argparse
import sys
import time
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Callable, List, Optional, Sequence, Tuple

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from quantum_coarse_graining.xy import (  # noqa: E402
    xy_block_characteristic,
    xy_middle_discriminant,
    xy_middle_quartic_coefficients,
)

# Paper notation.  ``A`` is the even-sector weight of the invariant witness and
# ``B = 1/2 - A`` the odd-sector weight; ``(u, u, w)`` are the Bloch eigenvalues
# of the effective channel; ``C = cos 2theta`` and ``S = sin 2theta``; lambda is
# the top eigenvalue of the Choi block ``M_3``.  The rescaled saddle variables
# are ``t = tan theta``, ``a = 2A`` and ``Delta = 4 lambda``.
A, U, W, C, S, LAM = sp.symbols("A u w C S lambda", real=True)
T, SMALL_A, DELTA = sp.symbols("t a Delta", real=True)
Z = sp.Symbol("z")
B = sp.Rational(1, 2) - A

# Entrywise symbols of the Hermitian block, before any parametrization.
ALPHA, GAMMA, GG, HH, BETA_R, BETA_I = sp.symbols(
    "alpha gamma g h beta_r beta_i", real=True
)
X = sp.Symbol("x")

# Certified windows.  The tangent window brackets the intermediate interval
# (tan theta_1 = 0.2637..., tan theta_2 = 0.7548...) and Delta = delta_A lies in
# the root window (1/2, 2).
T_WINDOW = (sp.Rational(1, 4), sp.Rational(19, 25))
DELTA_WINDOW = (sp.Rational(1, 2), sp.Integer(2))

PSI = U**2 - C * (1 + C) * U + C
NEGATIVE_QUADRATIC = C * U**2 + 2 * C * U - C - 2 * U**2

# Factors removed ("saturated") during the zero-multiplier elimination.  Every
# entry is nonvanishing at an interior facet saddle, where 0 < A < 1/2,
# 0 < u < 1 and 0 < C < 1.  The two sign-definite quadratics are certified by
# their discriminants inside ``check_zero_multiplier_elimination``.
SATURATORS: Tuple[Tuple[str, sp.Expr, str], ...] = (
    ("C", C, "C = cos 2theta > 0 on (0, pi/4)"),
    ("1-C", C - 1, "C < 1 for theta > 0"),
    ("1+C", C + 1, "C > -1"),
    ("A", A, "interior witness: A > 0"),
    ("1-2A", 2 * A - 1, "interior witness: A < 1/2, so 2A-1 = -2B < 0"),
    ("u", U, "interior facet: u > 0"),
    ("1-u", U - 1, "interior facet: u < 1"),
    ("C-u", U - C, "G = 0 forces (C-u)lambda = -B Psi, which is nonzero"),
    ("Psi", PSI, "monic in u with disc C(C-1)(C^2+3C+4) < 0, hence Psi > 0"),
    (
        "Cu^2+2Cu-C-2u^2",
        NEGATIVE_QUADRATIC,
        "leading coefficient C-2 < 0 and disc 8C(C-1) < 0, hence negative",
    ),
)


class Report:
    """Collects PASS/FAIL lines and the final exit status."""

    def __init__(self, stream=sys.stdout) -> None:
        self.stream = stream
        self.passed = 0
        self.failed = 0

    def check(self, claim: str, holds: bool) -> None:
        if holds:
            self.passed += 1
            self.stream.write("PASS  {0}\n".format(claim))
        else:
            self.failed += 1
            self.stream.write("FAIL  {0}\n".format(claim))
        self.stream.flush()

    def note(self, text: str) -> None:
        self.stream.write("      {0}\n".format(text))
        self.stream.flush()

    def section(self, name: str) -> None:
        self.stream.write("--- {0}\n".format(name))
        self.stream.flush()


# ---------------------------------------------------------------------------
# Exact helpers
# ---------------------------------------------------------------------------


def reduce_pythagorean(expression: sp.Expr) -> sp.Expr:
    """Reduce an expression modulo the exact unit-circle relation S^2 = 1-C^2."""
    polynomial = sp.Poly(sp.expand(expression), S)
    reduced = sp.Integer(0)
    for (power,), coefficient in polynomial.terms():
        quotient, remainder = divmod(power, 2)
        reduced += coefficient * (1 - C**2) ** quotient * S**remainder
    return sp.expand(reduced)


def is_zero(expression: sp.Expr) -> bool:
    """Exact zero test for polynomial and rational expressions."""
    return sp.expand(sp.cancel(sp.together(expression))) == 0


def bernstein_coefficients(
    expression: sp.Expr,
    box: Sequence[Tuple[sp.Symbol, sp.Rational, sp.Rational]],
) -> List[sp.Rational]:
    """Bernstein coefficients of a polynomial over a closed rational box.

    The Bernstein basis is nonnegative and forms a partition of unity, so the
    polynomial is bounded by the extreme coefficients.  A uniform strict sign
    therefore certifies that sign on the whole box.
    """
    polynomial = sp.expand(expression)
    for variable, low, high in box:
        polynomial = sp.expand(polynomial.subs(variable, low + (high - low) * variable))
    for variable, _, _ in box:
        degree = sp.Poly(polynomial, variable).degree()
        auxiliary = sp.Dummy("y")
        # (1+y)^n p(y/(1+y)) has the scaled Bernstein coefficients as its
        # ordinary coefficients.
        homogeneous = sp.Poly(
            sp.expand(
                sp.cancel(
                    sp.together(
                        (1 + auxiliary) ** degree
                        * polynomial.subs(variable, auxiliary / (1 + auxiliary))
                    )
                )
            ),
            auxiliary,
        )
        polynomial = sp.expand(
            sum(
                homogeneous.coeff_monomial(auxiliary**power)
                / sp.binomial(degree, power)
                * variable**power
                for power in range(degree + 1)
            )
        )
    return sp.Poly(polynomial, *[variable for variable, _, _ in box]).coeffs()


def certify_sign(
    expression: sp.Expr,
    box: Sequence[Tuple[sp.Symbol, sp.Rational, sp.Rational]],
    positive: bool,
    max_depth: int = 8,
    _depth: int = 0,
) -> bool:
    """Certify a strict sign on a closed rational box by Bernstein bisection."""
    coefficients = bernstein_coefficients(expression, box)
    if all(
        (coefficient > 0) if positive else (coefficient < 0)
        for coefficient in coefficients
    ):
        return True
    if _depth >= max_depth:
        return False
    widths = [high - low for _, low, high in box]
    split = max(range(len(box)), key=lambda index: widths[index])
    variable, low, high = box[split]
    middle = (low + high) / 2
    for bounds in ((low, middle), (middle, high)):
        halved = list(box)
        halved[split] = (variable, bounds[0], bounds[1])
        if not certify_sign(expression, halved, positive, max_depth, _depth + 1):
            return False
    return True


def sturm_root_count(
    expression: sp.Expr,
    variable: sp.Symbol,
    low: sp.Rational,
    high: sp.Rational,
) -> int:
    """Exact number of distinct real roots in ``(low, high]``.

    The Sturm chain is evaluated at the exact rational endpoints and the sign
    variations are counted over the rationals; no floating-point root finder is
    involved.
    """
    chain = [
        sp.Poly(term, variable)
        for term in sp.sturm(sp.Poly(sp.expand(expression), variable))
    ]
    if chain[0].eval(low) == 0 or chain[0].eval(high) == 0:
        raise ValueError("Sturm endpoints must not be roots")

    def variations(point: sp.Rational) -> int:
        signs = [sp.sign(term.eval(point)) for term in chain]
        signs = [sign for sign in signs if sign != 0]
        return sum(
            1 for index in range(len(signs) - 1) if signs[index] != signs[index + 1]
        )

    return int(variations(low) - variations(high))


def cauchy_bound(expression: sp.Expr, variable: sp.Symbol) -> sp.Rational:
    """A rational bound strictly larger than every real root."""
    coefficients = sp.Poly(sp.expand(expression), variable).all_coeffs()
    leading = coefficients[0]
    return 1 + max(abs(coefficient / leading) for coefficient in coefficients[1:])


def saturate(expression: sp.Expr) -> Tuple[sp.Expr, List[str]]:
    """Divide out the documented interior-nonvanishing factors."""
    _, factors = sp.factor_list(sp.expand(expression))
    kept: List[Tuple[sp.Expr, int]] = []
    removed: List[str] = []
    for factor, power in factors:
        name = None
        for candidate, saturator, _ in SATURATORS:
            if is_zero(factor - saturator) or is_zero(factor + saturator):
                name = candidate
                break
        if name is None:
            kept.append((factor, power))
        else:
            removed.append("({0})^{1}".format(name, power) if power > 1 else name)
    core = sp.expand(sp.prod([factor**power for factor, power in kept]))
    return core, removed


# ---------------------------------------------------------------------------
# The invariant Choi block and its characteristic polynomial
# ---------------------------------------------------------------------------


def generic_block() -> sp.Matrix:
    """The Hermitian ``M_3`` of the XY fixed-channel lemma, entrywise symbolic."""
    beta = BETA_R + sp.I * BETA_I
    return sp.Matrix(
        [
            [ALPHA, GAMMA * (GG + sp.I * HH), GAMMA * (-GG + sp.I * HH)],
            [GAMMA * (GG - sp.I * HH), 0, beta],
            [GAMMA * (-GG - sp.I * HH), sp.conjugate(beta), 0],
        ]
    )


def physical_substitution() -> dict:
    """The lemma parametrization of the block entries (gamma enters squared)."""
    return {
        ALPHA: A * (1 - W) / 2,
        GG: C - U,
        HH: S,
        BETA_R: B * (W - (C**2 - S**2)) / 2,
        BETA_I: B * S * C,
    }


@lru_cache(maxsize=None)
def generic_characteristic() -> sp.Expr:
    """``det(x I - M_3)`` with symbolic entries."""
    return sp.expand(sp.det(X * sp.eye(3) - generic_block()))


@lru_cache(maxsize=None)
def general_characteristic() -> sp.Expr:
    """``det(lambda I - M_3)`` in the physical parameters, before the facet."""
    squared = generic_characteristic().subs(GAMMA**2, A * B / 2)
    return sp.expand(squared.subs(physical_substitution()).subs(X, LAM))


@lru_cache(maxsize=None)
def facet_characteristic() -> sp.Expr:
    """The ``q_Z = 0`` characteristic polynomial ``P(lambda;A,u)``."""
    return reduce_pythagorean(general_characteristic().subs(W, 2 * U - 1))


@lru_cache(maxsize=None)
def compact_characteristic() -> sp.Expr:
    """The published compact form of ``P(lambda;A,u)``."""
    m = (U - C) ** 2 + S**2
    n = (U - C**2) ** 2 + C**2 * S**2
    theta = U * (1 - U) * (1 - C) ** 2
    return reduce_pythagorean(
        LAM**3
        - A * (1 - U) * LAM**2
        - (A * B * m + B**2 * n) * LAM
        - A * B**2 * theta
    )


@lru_cache(maxsize=None)
def middle_quartic() -> sp.Expr:
    """``Q_+(Delta,t)`` rebuilt from the dependency-free backend.

    ``xy_middle_quartic_coefficients`` returns rational numbers, so the symbolic
    quartic is recovered by exact Lagrange interpolation through nine rational
    tangents; its coefficients have degree at most seven in ``t``.
    """
    nodes = [sp.Rational(index, 3) for index in range(1, 11)]
    quartic = sp.Integer(0)
    for power in range(5):
        samples = []
        for node in nodes:
            value = xy_middle_quartic_coefficients(Fraction(int(node.p), int(node.q)))[
                4 - power
            ]
            samples.append((node, sp.Rational(value.numerator, value.denominator)))
        quartic += sp.expand(sp.interpolate(samples, T)) * DELTA**power
    return sp.expand(quartic)


def w_plus() -> sp.Expr:
    """The published sheet polynomial ``W_+(a,Delta,t)``."""
    return sp.expand(
        SMALL_A * (DELTA * (1 + T**2) + 4 * (1 - SMALL_A) * T**2) ** 2
        - 2 * DELTA**2 * (1 + T**2) ** 2
        + 8 * DELTA * (1 - SMALL_A) * T * (2 * SMALL_A * T**2 + 1 - T**2)
    )


@lru_cache(maxsize=None)
def scaled_facet_polynomial() -> Tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    """``F = 64(1+t^2)^2 P(Delta/4; a/2, u)`` and its coefficients in ``u``."""
    half_angle = {
        C: (1 - T**2) / (1 + T**2),
        S: 2 * T / (1 + T**2),
        A: SMALL_A / 2,
        LAM: DELTA / 4,
    }
    scaled = sp.expand(
        sp.cancel(
            sp.together(
                64 * (1 + T**2) ** 2 * compact_characteristic().subs(half_angle)
            )
        )
    )
    second, first, zeroth = sp.Poly(scaled, U).all_coeffs()
    return scaled, sp.expand(second), sp.expand(first), sp.expand(zeroth)


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------


def check_m3_characteristic(report: Report) -> None:
    """Item 1: the block characteristic polynomial and its facet form."""
    block = generic_block()
    report.check(
        "m3: M_3 is Hermitian",
        sp.simplify(block - block.conjugate().T) == sp.zeros(3, 3),
    )

    _, second, first, zeroth = sp.Poly(generic_characteristic(), X).all_coeffs()
    modulus = BETA_R**2 + BETA_I**2
    report.check(
        "m3: the spectrum depends on gamma only through gamma^2 = AB/2",
        all(
            set(sp.Poly(coefficient, GAMMA).as_dict().keys()) <= {(0,), (2,)}
            for coefficient in (second, first, zeroth)
        ),
    )
    report.check("m3: p_2 = -alpha", is_zero(second + ALPHA))
    report.check(
        "m3: p_1 = -[2 gamma^2 (g^2+h^2) + |beta|^2]",
        is_zero(first + (2 * GAMMA**2 * (GG**2 + HH**2) + modulus)),
    )
    report.check(
        "m3: p_0 = alpha|beta|^2 + 2 gamma^2 (g^2-h^2) Re beta"
        " - 4 gamma^2 g h Im beta",
        is_zero(
            zeroth
            - (
                ALPHA * modulus
                + 2 * GAMMA**2 * (GG**2 - HH**2) * BETA_R
                - 4 * GAMMA**2 * GG * HH * BETA_I
            )
        ),
    )

    facet = facet_characteristic()
    report.check(
        "facet: P = lambda^3 - A(1-u)lambda^2 - (ABm+B^2N)lambda - AB^2 Theta"
        "  (mod S^2 = 1-C^2)",
        is_zero(facet - compact_characteristic()),
    )

    _, trace, minors, determinant = sp.Poly(facet, LAM).all_coeffs()
    m = reduce_pythagorean((U - C) ** 2 + S**2)
    n = reduce_pythagorean((U - C**2) ** 2 + C**2 * S**2)
    theta = U * (1 - U) * (1 - C) ** 2
    report.check("facet: Tr M_3 = A(1-u) > 0", is_zero(-trace - A * (1 - U)))
    report.check(
        "facet: sum of order-two principal minors = -(ABm + B^2N) < 0",
        is_zero(minors + (A * B * m + B**2 * n)),
    )
    report.check(
        "facet: det M_3 = A B^2 u(1-u)(1-C)^2 > 0",
        is_zero(-determinant - A * B**2 * theta),
    )
    report.check(
        "facet: m and N are strictly positive sums of squares",
        is_zero(m - ((U - C) ** 2 + 1 - C**2))
        and is_zero(n - ((U - C**2) ** 2 + C**2 * (1 - C**2))),
    )

    # Descartes on lambda^3 - e1 lambda^2 + e2 lambda - e3 with e1 > 0, e2 < 0
    # and e3 > 0 gives exactly one positive root; the other two multiply to a
    # positive number, so they are both negative.  An exact Sturm count at an
    # interior rational point of the middle interval confirms the inertia.
    interior = {
        C: sp.Rational(3, 5),
        S: sp.Rational(4, 5),
        A: sp.Rational(1, 4),
        U: sp.Rational(1, 2),
    }
    cubic = sp.expand(facet.subs(interior))
    bound = cauchy_bound(cubic, LAM)
    report.check(
        "facet: inertia (1 positive, 2 negative) by exact Sturm count at"
        " (t,A,u) = (1/2,1/4,1/2)",
        sturm_root_count(cubic, LAM, sp.Integer(0), bound) == 1
        and sturm_root_count(cubic, LAM, -bound, sp.Integer(0)) == 2,
    )

    backend = xy_block_characteristic(
        Fraction(9, 41),
        Fraction(40, 41),
        Fraction(729, 1681),
        Fraction(81, 1681),
        Fraction(1, 4),
    )
    derived = sp.Poly(
        general_characteristic().subs(
            {
                C: sp.Rational(9, 41),
                S: sp.Rational(40, 41),
                U: sp.Rational(729, 1681),
                W: sp.Rational(81, 1681),
                A: sp.Rational(1, 4),
            }
        ),
        LAM,
    ).all_coeffs()
    report.check(
        "backend: xy_block_characteristic reproduces the derived coefficients",
        [sp.Rational(value.numerator, value.denominator) for value in backend]
        == derived,
    )


def check_normal_derivative(report: Report) -> None:
    """Item 2: the normal derivative identity and the tangential relation."""
    general = general_characteristic()
    normal = reduce_pythagorean(sp.diff(general, U).subs(W, 2 * U - 1))
    multiplier = (C - U) * LAM + B * PSI
    report.check(
        "kkt: P_u at fixed w equals 2 A B G, G = (C-u)lambda + B Psi",
        is_zero(normal - 2 * A * B * multiplier),
    )
    report.check(
        "kkt: Psi = u^2 - C(1+C)u + C has disc C(C-1)(C^2+3C+4) < 0, so Psi > 0",
        is_zero(sp.discriminant(PSI, U) - C * (C - 1) * (C**2 + 3 * C + 4)),
    )

    tangential = reduce_pythagorean(
        sp.diff(general, U).subs(W, 2 * U - 1)
        + 2 * sp.diff(general, W).subs(W, 2 * U - 1)
    )
    report.check(
        "kkt: the facet tangent P_u + 2 P_w equals d/du P(lambda;A,u)",
        is_zero(tangential - sp.diff(facet_characteristic(), U)),
    )
    report.check(
        "kkt: G = C(1-C)[(1+C)+SC]/2 > 0 at the limiting weak saddle"
        " (A,u,lambda) = (0,C^2,SC/2)",
        is_zero(
            multiplier.subs({A: 0, U: C**2, LAM: S * C / 2})
            - C * (1 - C) * ((1 + C) + S * C) / 2
        ),
    )


def check_zero_multiplier_elimination(report: Report) -> None:
    """Item 3: eliminate a hypothetical zero KKT multiplier."""
    report.check(
        "saturation: Cu^2+2Cu-C-2u^2 has disc 8C(C-1) < 0 and leading C-2 < 0",
        is_zero(sp.discriminant(NEGATIVE_QUADRATIC, U) - 8 * C * (C - 1)),
    )
    report.note("excluded factors (each nonzero at an interior facet saddle):")
    for name, _, justification in SATURATORS:
        report.note("  {0:16s} {1}".format(name, justification))

    facet = facet_characteristic()
    general = general_characteristic()
    stationary = (
        ("P", facet),
        ("P_A", sp.diff(facet, A)),
        (
            "P_u+2P_w",
            reduce_pythagorean(
                sp.diff(general, U).subs(W, 2 * U - 1)
                + 2 * sp.diff(general, W).subs(W, 2 * U - 1)
            ),
        ),
    )

    def eliminate_lambda(expression: sp.Expr) -> sp.Expr:
        """Insert the unique root of ``G`` and clear its denominator."""
        degree = sp.Poly(expression, LAM).degree()
        return sp.expand(
            sp.cancel(
                sp.together(
                    ((C - U) ** degree * expression).subs(LAM, -B * PSI / (C - U))
                )
            )
        )

    cores = {}
    for name, expression in stationary:
        core, removed = saturate(eliminate_lambda(expression))
        cores[name] = core
        report.check(
            "elim: {0} on G = 0 is linear in A after saturating {1}".format(
                name, "{" + ", ".join(removed) + "}"
            ),
            sp.Poly(core, A).degree() == 1,
        )

    first_core, first_removed = saturate(
        sp.factor(sp.resultant(sp.Poly(cores["P"], A), sp.Poly(cores["P_u+2P_w"], A)))
    )
    second_core, second_removed = saturate(
        sp.factor(sp.resultant(sp.Poly(cores["P"], A), sp.Poly(cores["P_A"], A)))
    )
    report.check(
        "elim: Res_A(P, P_u+2P_w) saturates to u^2 - 2Cu + C^4",
        is_zero(first_core - (U**2 - 2 * C * U + C**4)),
    )
    report.check(
        "elim: Res_A(P, P_A) saturates to a quartic in u",
        sp.Poly(second_core, U).degree() == 4,
    )
    report.note(
        "resultant saturations: {0} and {1}".format(
            "{" + ", ".join(first_removed) + "}", "{" + ", ".join(second_removed) + "}"
        )
    )

    tangent_resultant = sp.factor(
        sp.resultant(sp.Poly(first_core, U), sp.Poly(second_core, U))
    )
    expected = (
        C**4
        * (C - 1) ** 4
        * (C + 1) ** 3
        * (C**2 + 3)
        * (5 * C**3 + C**2 + 3 * C - 1)
    )
    report.check(
        "elim: Res_u = C^4 (C-1)^4 (C+1)^3 (C^2+3) (5C^3+C^2+3C-1)",
        is_zero(tangent_resultant - expected),
    )

    published = (
        C**6
        * (C - 1) ** 14
        * (C + 1) ** 3
        * (C**2 + 3)
        * (5 * C**3 + C**2 + 3 * C - 1)
        * (20 * C**4 + 32 * C**3 + 29 * C**2 + 16 * C + 3)
    )
    quotient, remainder = sp.div(sp.expand(published), sp.expand(tangent_resultant), C)
    report.check(
        "elim: the published product equals this resultant times"
        " C^2 (C-1)^10 (20C^4+32C^3+29C^2+16C+3)",
        sp.expand(remainder) == 0
        and is_zero(
            quotient
            - C**2 * (C - 1) ** 10 * (20 * C**4 + 32 * C**3 + 29 * C**2 + 16 * C + 3)
        ),
    )

    # Sharper cross-check: the saturated elimination ideal itself.
    saturator = sp.prod([expression for _, expression, _ in SATURATORS])
    basis = sp.groebner(
        [cores["P"], cores["P_A"], cores["P_u+2P_w"], Z * saturator - 1],
        Z,
        A,
        U,
        C,
        order="lex",
    )
    eliminated = [
        sp.factor(generator)
        for generator in basis.exprs
        if generator.free_symbols <= {C}
    ]
    report.check(
        "elim: the saturated Groebner elimination ideal in C is <5C^3+C^2+3C-1>",
        len(eliminated) == 1
        and is_zero(eliminated[0] - (5 * C**3 + C**2 + 3 * C - 1)),
    )
    locus = [
        sp.factor(generator)
        for generator in basis.exprs
        if generator.free_symbols <= {U, C} and U in generator.free_symbols
    ]
    report.check(
        "elim: the zero-multiplier locus also forces u = (1+C^2)/2",
        len(locus) == 1
        and (
            is_zero(locus[0] - (2 * U - 1 - C**2))
            or is_zero(locus[0] + (2 * U - 1 - C**2))
        ),
    )
    report.check(
        "elim: u = (1+C^2)/2 is the strong-branch channel u = C(1+S) exactly"
        " when 1+C^2 = 2C(1+S)",
        is_zero(
            (2 * U - 1 - C**2).subs(U, C * (1 + S)) + (1 + C**2 - 2 * C * (1 + S))
        ),
    )

    half_angle = {C: (1 - T**2) / (1 + T**2)}
    threshold = sp.cancel(
        (5 * C**3 + C**2 + 3 * C - 1).subs(half_angle) * (1 + T**2) ** 3
    )
    report.check(
        "elim: 5C^3+C^2+3C-1 = -8(t^3-t^2+1)(t^3+t^2-1)/(1+t^2)^3",
        is_zero(threshold + 8 * (T**3 - T**2 + 1) * (T**3 + T**2 - 1)),
    )
    report.check(
        "elim: t^3-t^2+1 > 0 on [0,1], so only theta_2 (t^3+t^2-1) survives",
        certify_sign(T**3 - T**2 + 1, [(T, sp.Integer(0), sp.Integer(1))], True),
    )


def check_t_parametrization(report: Report) -> None:
    """Item 4: the t parametrization, the W_+ split and Res_a(W_+,G_8)."""
    scaled, second, first, zeroth = scaled_facet_polynomial()
    report.check(
        "t-form: F = 64(1+t^2)^2 P(Delta/4; a/2, u) is polynomial and quadratic in u",
        sp.denom(sp.cancel(scaled)) == 1 and sp.Poly(scaled, U).degree() == 2,
    )

    positive_sheet = w_plus()
    negative_sheet = sp.expand(positive_sheet.subs(T, -T))
    report.check(
        "t-form: disc_u F = 4 W_+(a,Delta,t) W_+(a,Delta,-t)",
        is_zero(first**2 - 4 * second * zeroth - 4 * positive_sheet * negative_sheet),
    )

    eight = sp.expand(
        -(
            sp.diff(second, SMALL_A) * first**2
            - 2 * first * sp.diff(first, SMALL_A) * second
            + 4 * sp.diff(zeroth, SMALL_A) * second**2
        )
        / 16
    )
    report.check(
        "t-form: -16 G_8 = 4 P_2^2 F_a at the F_u = 0 root u = -P_1/(2P_2)",
        is_zero(
            -16 * eight
            - 4 * second**2 * sp.diff(scaled, SMALL_A).subs(U, -first / (2 * second))
        ),
    )
    report.check("t-form: G_8 has degree 8 in a", sp.Poly(eight, SMALL_A).degree() == 8)

    resultant = sp.resultant(sp.Poly(positive_sheet, SMALL_A), sp.Poly(eight, SMALL_A))
    first_cofactor = DELTA**2 * (1 + 3 * T**2) + 4 * (1 - T**2) * (1 + DELTA)
    second_cofactor = DELTA * (1 + T**2) ** 2 + 2 * (1 + T) ** 2 * (1 + 2 * T - T**2)
    report.check(
        "t-form: Res_a(W_+,G_8) = -2^36 Delta^12 t^28 (1+t^2)^10 E_1 E_2^2 Q_+",
        is_zero(
            resultant
            + 2**36
            * DELTA**12
            * T**28
            * (1 + T**2) ** 10
            * first_cofactor
            * second_cofactor**2
            * middle_quartic()
        ),
    )

    for tangent in (sp.Rational(1, 2), sp.Rational(19, 25)):
        backend = xy_middle_quartic_coefficients(
            Fraction(int(tangent.p), int(tangent.q))
        )
        derived = sp.Poly(middle_quartic().subs(T, tangent), DELTA).all_coeffs()
        report.check(
            "backend: xy_middle_quartic_coefficients({0}) matches the"
            " eliminated Q_+".format(tangent),
            [sp.Rational(value.numerator, value.denominator) for value in backend]
            == derived,
        )


def check_sheet_selection(report: Report) -> None:
    """Item 4: the physical saddle lies on W_+ = 0, with positive cofactors."""
    negative_sheet = sp.expand(w_plus().subs(T, -T))
    box = [
        (SMALL_A, sp.Integer(0), sp.Integer(1)),
        (DELTA, DELTA_WINDOW[0], DELTA_WINDOW[1]),
        (T, T_WINDOW[0], T_WINDOW[1]),
    ]
    report.check(
        "sheet: W_+(a,Delta,-t) < 0 on a in [0,1], Delta in [1/2,2], t in [1/4,19/25],"
        " so disc_u F = 0 forces W_+ = 0",
        certify_sign(negative_sheet, box, False),
    )
    report.check(
        "sheet: -W_+(a,Delta,-t)/(1+t^2)^2 >= 171381/972196 on that box",
        certify_sign(
            sp.expand(
                -negative_sheet - sp.Rational(171381, 972196) * (1 + T**2) ** 2
            ),
            box,
            True,
        ),
    )

    window = [(T, T_WINDOW[0], T_WINDOW[1])]
    first_cofactor = sp.Poly(
        sp.expand(DELTA**2 * (1 + 3 * T**2) + 4 * (1 - T**2) * (1 + DELTA)), DELTA
    ).all_coeffs()
    second_cofactor = sp.Poly(
        sp.expand(DELTA * (1 + T**2) ** 2 + 2 * (1 + T) ** 2 * (1 + 2 * T - T**2)),
        DELTA,
    ).all_coeffs()
    report.check(
        "sheet: every Delta-coefficient of E_1 is positive on the window, so E_1 > 0",
        all(certify_sign(coefficient, window, True) for coefficient in first_cofactor),
    )
    report.check(
        "sheet: every Delta-coefficient of E_2 is positive on the window, so E_2 > 0",
        all(certify_sign(coefficient, window, True) for coefficient in second_cofactor),
    )


def check_quartic_discriminant(report: Report) -> None:
    """Item 5: the quartic discriminant factorization and noncollision."""
    quartic = middle_quartic()
    discriminant = sp.factor(sp.discriminant(sp.Poly(quartic, DELTA)))
    ninth = (
        27 * T**9
        + 378 * T**8
        + 1026 * T**7
        + 2106 * T**6
        + 2772 * T**5
        + 2646 * T**4
        + 1658 * T**3
        + 678 * T**2
        + 141 * T
        + 16
    )
    report.check(
        "disc: disc_Delta Q_+ = 2^16 t^7 (1+t^2) (t^3-t^2-t-1)^2 N_9(t)^3",
        is_zero(
            discriminant
            - 2**16 * T**7 * (1 + T**2) * (T**3 - T**2 - T - 1) ** 2 * ninth**3
        ),
    )
    window = [(T, T_WINDOW[0], T_WINDOW[1])]
    report.check(
        "disc: t^3-t^2-t-1 < 0 and N_9 > 0 on the window, so the root count is fixed",
        certify_sign(T**3 - T**2 - T - 1, window, False)
        and certify_sign(ninth, window, True),
    )
    backend = xy_middle_discriminant(Fraction(1, 2))
    report.check(
        "backend: xy_middle_discriminant(1/2) matches the derived discriminant",
        sp.Rational(backend.numerator, backend.denominator)
        == sp.discriminant(sp.Poly(quartic.subs(T, sp.Rational(1, 2)), DELTA)),
    )


def check_root_window(report: Report) -> None:
    """Item 5: the exact endpoint substitutions Delta = 1/2 and Delta = 2."""
    quartic = middle_quartic()
    lower_reference = (
        2484 * T**7
        + 2565 * T**6
        + 2349 * T**5
        + 3119 * T**4
        + 2014 * T**3
        + 687 * T**2
        - 155 * T
        + 5
    )
    upper_reference = 9 * T**5 + 5 * T**4 + 22 * T**3 + 12 * T**2 + 7 * T - 1
    report.check(
        "window: 16 Q_+(1/2,t) = 2484t^7+2565t^6+2349t^5+3119t^4+2014t^3"
        "+687t^2-155t+5",
        is_zero(16 * quartic.subs(DELTA, sp.Rational(1, 2)) - lower_reference),
    )
    report.check(
        "window: Q_+(2,t) = -32(9t^5+5t^4+22t^3+12t^2+7t-1)",
        is_zero(quartic.subs(DELTA, 2) + 32 * upper_reference),
    )
    window = [(T, T_WINDOW[0], T_WINDOW[1])]
    report.check(
        "window: Q_+(1/2,t) > 0 and Q_+(2,t) < 0 on [1/4,19/25]",
        certify_sign(lower_reference, window, True)
        and certify_sign(upper_reference, window, True),
    )


def check_threshold_joining(report: Report) -> None:
    """Item 5: the closed branches join the middle quartic exactly."""
    quartic = middle_quartic()
    weak_value = 4 * T * (1 - T**2) / (1 + T**2) ** 2
    strong_value = 2 * T * (1 + T + T**2) / (1 + T**2) ** 2
    first_threshold = T**6 - 2 * T**5 - 3 * T**4 + 7 * T**2 + 2 * T - 1
    second_threshold = T**3 + T**2 - 1
    weak_positive = T**4 + 2 * T**3 + 4 * T**2 + 2 * T + 1
    strong_positive = (
        8 * T**10
        + 4 * T**9
        + 35 * T**8
        + 43 * T**7
        + 108 * T**6
        + 120 * T**5
        + 143 * T**4
        + 100 * T**3
        + 66 * T**2
        + 23 * T
        + 9
    )

    half_angle = {C: (1 - T**2) / (1 + T**2), S: 2 * T / (1 + T**2)}
    report.check(
        "join: sin 4theta = 4t(1-t^2)/(1+t^2)^2 and S+S^2/2 = 2t(1+t+t^2)/(1+t^2)^2",
        is_zero((2 * S * C).subs(half_angle) - weak_value)
        and is_zero((S + S**2 / 2).subs(half_angle) - strong_value),
    )
    report.check(
        "join: Q_+(sin 4theta, t) = 128 t^3 (t^4+2t^3+4t^2+2t+1)(1+t^2)^-6"
        " (t^6-2t^5-3t^4+7t^2+2t-1)^2",
        is_zero(
            quartic.subs(DELTA, weak_value)
            - 128 * T**3 * weak_positive * first_threshold**2 / (1 + T**2) ** 6
        ),
    )
    report.check(
        "join: Q_+(S+S^2/2, t) = 16 t^3 N_10(t) (1+t^2)^-6 (t^3+t^2-1)^2",
        is_zero(
            quartic.subs(DELTA, strong_value)
            - 16 * T**3 * strong_positive * second_threshold**2 / (1 + T**2) ** 6
        ),
    )
    window = [(T, T_WINDOW[0], T_WINDOW[1])]
    report.check(
        "join: both cofactors are strictly positive on the window",
        certify_sign(weak_positive, window, True)
        and certify_sign(strong_positive, window, True),
    )
    report.check(
        "join: 2SC-(3C^3+C^2-C-1) = 2(t^6-2t^5-3t^4+7t^2+2t-1)/(1+t^2)^3",
        is_zero(
            sp.cancel(
                (2 * S * C - (3 * C**3 + C**2 - C - 1)).subs(half_angle)
                * (1 + T**2) ** 3
            )
            - 2 * first_threshold
        ),
    )
    report.check(
        "join: 1+C^2-2C(1+S) = 4t(t^3+t^2-1)/(1+t^2)^2",
        is_zero(
            sp.cancel((1 + C**2 - 2 * C * (1 + S)).subs(half_angle) * (1 + T**2) ** 2)
            - 4 * T * second_threshold
        ),
    )


def check_sturm_root_count(report: Report) -> None:
    """Item 6: an exact Sturm count, never a floating-point root call."""
    quartic = middle_quartic()
    for tangent in (T_WINDOW[0], sp.Rational(1, 2), T_WINDOW[1]):
        specialized = sp.expand(quartic.subs(T, tangent))
        report.check(
            "sturm: Q_+(Delta,{0}) has exactly one root in (1/2,2)".format(tangent),
            sturm_root_count(specialized, DELTA, DELTA_WINDOW[0], DELTA_WINDOW[1]) == 1,
        )


CHECKS: Tuple[Tuple[str, Callable[[Report], None]], ...] = (
    ("m3-characteristic", check_m3_characteristic),
    ("normal-derivative", check_normal_derivative),
    ("zero-multiplier", check_zero_multiplier_elimination),
    ("t-parametrization", check_t_parametrization),
    ("sheet-selection", check_sheet_selection),
    ("quartic-discriminant", check_quartic_discriminant),
    ("root-window", check_root_window),
    ("threshold-joining", check_threshold_joining),
    ("sturm-count", check_sturm_root_count),
)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Exact SymPy derivation of the intermediate XY theorem."
    )
    parser.add_argument(
        "--list", action="store_true", help="list the available derivations and exit"
    )
    parser.add_argument(
        "--only",
        action="append",
        metavar="NAME",
        help="run only the named derivation (repeatable)",
    )
    arguments = parser.parse_args(argv)

    if arguments.list:
        for name, function in CHECKS:
            print("{0:22s} {1}".format(name, (function.__doc__ or "").splitlines()[0]))
        return 0

    selected = CHECKS
    if arguments.only:
        unknown = set(arguments.only) - {name for name, _ in CHECKS}
        if unknown:
            parser.error("unknown derivation(s): " + ", ".join(sorted(unknown)))
        selected = tuple(item for item in CHECKS if item[0] in set(arguments.only))

    report = Report()
    started = time.time()
    for name, function in selected:
        report.section(name)
        function(report)
    elapsed = time.time() - started
    print(
        "{0} passed, {1} failed in {2:.1f} s".format(
            report.passed, report.failed, elapsed
        )
    )
    return 1 if report.failed else 0


if __name__ == "__main__":
    sys.exit(main())
