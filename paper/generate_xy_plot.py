"""Generate signalling curves and optimizer phase diagrams."""

from __future__ import annotations

import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
THETA_1 = 0.258270520262
THETA_2 = 0.646615513406
GOLDEN_RATIO = (1 + math.sqrt(5)) / 2


def golden_minimum(function, lower: float, upper: float) -> float:
    """Return a deterministic golden-section minimizer."""
    left = upper - (upper - lower) / GOLDEN_RATIO
    right = lower + (upper - lower) / GOLDEN_RATIO
    left_value = function(left)
    right_value = function(right)
    for _ in range(55):
        if left_value <= right_value:
            upper = right
            right = left
            right_value = left_value
            left = upper - (upper - lower) / GOLDEN_RATIO
            left_value = function(left)
        else:
            lower = left
            left = right
            left_value = right_value
            right = lower + (upper - lower) / GOLDEN_RATIO
            right_value = function(right)
    return (lower + upper) / 2


def golden_maximum(function, lower: float, upper: float) -> float:
    """Return a deterministic golden-section maximizer."""
    return golden_minimum(lambda value: -function(value), lower, upper)


def quartic_coefficients(tangent: float) -> tuple[float, ...]:
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


def evaluate_polynomial(coefficients: tuple[float, ...], value: float) -> float:
    result = 0.0
    for coefficient in coefficients:
        result = result * value + coefficient
    return result


def middle_defect(theta: float) -> float:
    coefficients = quartic_coefficients(math.tan(theta))
    lower, upper = 0.5, 2.0
    lower_value = evaluate_polynomial(coefficients, lower)
    for _ in range(80):
        midpoint = (lower + upper) / 2
        midpoint_value = evaluate_polynomial(coefficients, midpoint)
        if midpoint_value == 0:
            return midpoint
        if midpoint_value * lower_value > 0:
            lower = midpoint
            lower_value = midpoint_value
        else:
            upper = midpoint
    return (lower + upper) / 2


def signalling_defect(theta: float) -> tuple[float, str]:
    sine = math.sin(2 * theta)
    cosine = math.cos(2 * theta)
    if theta <= THETA_1:
        return 2 * sine * cosine, "weak"
    if theta < THETA_2:
        return middle_defect(theta), "middle"
    return sine + sine * sine / 2, "strong"


def partial_swap_fixed_defect(sine: float, shrinkage: float) -> float:
    """Qubit fixed-channel norm from the analytic scalar reduction."""
    a = (1 - shrinkage) / 3
    radical = (1 + shrinkage) ** 2 - 4 * shrinkage * (1 - sine * sine)
    if radical <= 4 * a * a:
        return 4 * a
    return radical / (math.sqrt(radical) - a)


def partial_swap_optimizer(phi: float) -> tuple[float, float]:
    """Return optimal shrinkage and symmetric witness weight for qubits."""
    sine = math.sin(phi)
    if sine <= 1 / 3:
        shrinkage = 1.0
    elif phi >= math.pi / 2 - 1e-12:
        shrinkage = 0.0
    else:
        shrinkage = golden_minimum(
            lambda value: partial_swap_fixed_defect(sine, value),
            0.0,
            1.0,
        )

    a = (1 - shrinkage) / 3
    radical = (1 + shrinkage) ** 2 - 4 * shrinkage * math.cos(phi) ** 2
    if radical <= 1e-15:
        witness_coordinate = 0.0
    elif radical <= 4 * a * a:
        witness_coordinate = 1.0
    else:
        root = math.sqrt(radical)
        witness_coordinate = min(1.0, a / (root - a))
    symmetric_weight = (1 + witness_coordinate) / 2
    return shrinkage, symmetric_weight


def xy_characteristic(
    eigenvalue: float,
    theta: float,
    transverse: float,
    even_weight: float,
) -> float:
    """Characteristic polynomial of the middle-facet M3 block."""
    cosine = math.cos(2 * theta)
    sine = math.sin(2 * theta)
    odd_weight = 0.5 - even_weight
    m = (transverse - cosine) ** 2 + sine * sine
    n_term = (
        (transverse - cosine * cosine) ** 2
        + cosine * cosine * sine * sine
    )
    theta_term = transverse * (1 - transverse) * (1 - cosine) ** 2
    return (
        eigenvalue**3
        - even_weight * (1 - transverse) * eigenvalue**2
        - (even_weight * odd_weight * m + odd_weight**2 * n_term)
        * eigenvalue
        - even_weight * odd_weight**2 * theta_term
    )


def xy_largest_eigenvalue(
    theta: float,
    transverse: float,
    even_weight: float,
) -> float:
    """Find the unique positive root of the middle-facet cubic."""
    lower, upper = 0.0, 1.0
    while xy_characteristic(upper, theta, transverse, even_weight) < 0:
        upper *= 2
    for _ in range(64):
        midpoint = (lower + upper) / 2
        if xy_characteristic(midpoint, theta, transverse, even_weight) < 0:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2


def xy_middle_optimizer(theta: float) -> tuple[float, float]:
    """Numerically reconstruct the unique middle-facet saddle."""
    def witness_value(transverse: float) -> tuple[float, float]:
        even_weight = golden_maximum(
            lambda value: xy_largest_eigenvalue(theta, transverse, value),
            0.0,
            0.5,
        )
        return (
            xy_largest_eigenvalue(theta, transverse, even_weight),
            even_weight,
        )

    transverse = golden_minimum(
        lambda value: witness_value(value)[0],
        0.0,
        1.0,
    )
    return transverse, witness_value(transverse)[1]


def xy_optimizer(theta: float) -> tuple[float, float, float, str]:
    """Return transverse/longitudinal channel eigenvalues and witness A."""
    sine = math.sin(2 * theta)
    cosine = math.cos(2 * theta)
    if theta <= THETA_1:
        transverse = cosine * cosine
        return transverse, cosine * cosine - sine * sine, 0.0, "weak"
    if theta < THETA_2:
        transverse, even_weight = xy_middle_optimizer(theta)
        return transverse, 2 * transverse - 1, even_weight, "middle"
    return cosine * (1 + sine), cosine * cosine, 0.25, "strong"


def points(count: int = 241) -> list[tuple[float, float, str]]:
    return [
        (theta, *signalling_defect(theta))
        for theta in (index * math.pi / (4 * (count - 1)) for index in range(count))
    ]


def partial_swap_phase_points(
    count: int = 161,
) -> list[tuple[float, float, float]]:
    return [
        (phi, *partial_swap_optimizer(phi))
        for phi in (
            index * math.pi / (2 * (count - 1))
            for index in range(count)
        )
    ]


def xy_phase_points(
    count: int = 81,
) -> list[tuple[float, float, float, float, str]]:
    return [
        (theta, *xy_optimizer(theta))
        for theta in (
            index * math.pi / (4 * (count - 1))
            for index in range(count)
        )
    ]


def write_data(values: list[tuple[float, float, str]]) -> None:
    lines = ["theta defect regime"]
    lines.extend(f"{theta:.12f} {defect:.12f} {regime}" for theta, defect, regime in values)
    (ROOT / "xy_curve.dat").write_text("\n".join(lines) + "\n", encoding="ascii")


def write_phase_data(
    partial_swap_values: list[tuple[float, float, float]],
    xy_values: list[tuple[float, float, float, float, str]],
) -> None:
    partial_swap_lines = ["phi shrinkage symmetric_weight"]
    partial_swap_lines.extend(
        f"{phi:.12f} {shrinkage:.12f} {weight:.12f}"
        for phi, shrinkage, weight in partial_swap_values
    )
    (ROOT / "partial_swap_phase.dat").write_text(
        "\n".join(partial_swap_lines) + "\n",
        encoding="ascii",
    )

    xy_lines = ["theta transverse longitudinal even_weight regime"]
    xy_lines.extend(
        f"{theta:.12f} {transverse:.12f} {longitudinal:.12f} "
        f"{even_weight:.12f} {regime}"
        for theta, transverse, longitudinal, even_weight, regime in xy_values
    )
    (ROOT / "xy_phase.dat").write_text(
        "\n".join(xy_lines) + "\n",
        encoding="ascii",
    )


def write_svg(values: list[tuple[float, float, str]]) -> None:
    width, height = 760, 430
    left, right, top, bottom = 70, 24, 28, 58
    plot_width = width - left - right
    plot_height = height - top - bottom

    def x(theta: float) -> float:
        return left + theta / (math.pi / 4) * plot_width

    def y(defect: float) -> float:
        return top + (1.6 - defect) / 1.6 * plot_height

    paths = []
    colors = {"weak": "#2878b5", "middle": "#d97b29", "strong": "#3c9d5d"}
    for regime in ("weak", "middle", "strong"):
        regime_points = [(x(theta), y(defect)) for theta, defect, label in values if label == regime]
        coordinates = " ".join(
            f"{'M' if index == 0 else 'L'} {px:.2f} {py:.2f}"
            for index, (px, py) in enumerate(regime_points)
        )
        paths.append(
            f'<path d="{coordinates}" fill="none" stroke="{colors[regime]}" '
            'stroke-width="4" stroke-linecap="round"/>'
        )

    theta_ticks = (
        (0, "0"),
        (THETA_1, "θ₁"),
        (THETA_2, "θ₂"),
        (math.pi / 4, "π/4"),
    )
    x_ticks = "\n".join(
        f'<line x1="{x(value):.2f}" y1="{height-bottom}" x2="{x(value):.2f}" '
        f'y2="{height-bottom+7}" stroke="#333"/>'
        f'<text x="{x(value):.2f}" y="{height-bottom+27}" text-anchor="middle">{label}</text>'
        for value, label in theta_ticks
    )
    y_ticks = "\n".join(
        f'<line x1="{left-7}" y1="{y(value):.2f}" x2="{left}" y2="{y(value):.2f}" stroke="#333"/>'
        f'<text x="{left-12}" y="{y(value)+5:.2f}" text-anchor="end">{value:g}</text>'
        for value in (0, 0.5, 1.0, 1.5)
    )
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}"
viewBox="0 0 {width} {height}" role="img" aria-labelledby="title description">
<title id="title">Exact XY directional signalling curve</title>
<desc id="description">The signalling defect rises through weak, quartic middle, and iSWAP regimes.</desc>
<rect width="100%" height="100%" fill="white"/>
<g font-family="sans-serif" font-size="16" fill="#222">
<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="#333"/>
<line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="#333"/>
{x_ticks}
{y_ticks}
<line x1="{left}" y1="{y(1.5):.2f}" x2="{width-right}" y2="{y(1.5):.2f}"
 stroke="#666" stroke-width="2" stroke-dasharray="3 6"/>
{''.join(paths)}
<text x="{left + plot_width/2:.2f}" y="{height-10}" text-anchor="middle">interaction angle θ</text>
<text x="18" y="{top + plot_height/2:.2f}" text-anchor="middle"
 transform="rotate(-90 18 {top + plot_height/2:.2f})">unavoidable error δₐ</text>
<text x="{x(0.12):.2f}" y="{y(0.25):.2f}" fill="{colors['weak']}">weak</text>
<text x="{x(0.42):.2f}" y="{y(1.22):.2f}" fill="{colors['middle']}">quartic</text>
<text x="{x(0.68):.2f}" y="{y(1.42):.2f}" fill="{colors['strong']}">iSWAP</text>
<text x="{left+8}" y="{y(1.5)-8:.2f}" fill="#555">two-qubit ceiling</text>
</g>
</svg>
"""
    (ROOT / "xy_curve.svg").write_text(svg, encoding="utf-8")


def main() -> None:
    values = points()
    write_data(values)
    write_phase_data(partial_swap_phase_points(), xy_phase_points())
    write_svg(values)


if __name__ == "__main__":
    main()
