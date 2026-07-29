"""Generate the XY signalling curve as data and a lightweight SVG."""

from __future__ import annotations

import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
THETA_1 = 0.258270520262
THETA_2 = 0.646615513406


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


def points(count: int = 241) -> list[tuple[float, float, str]]:
    return [
        (theta, *signalling_defect(theta))
        for theta in (index * math.pi / (4 * (count - 1)) for index in range(count))
    ]


def write_data(values: list[tuple[float, float, str]]) -> None:
    lines = ["theta defect regime"]
    lines.extend(f"{theta:.12f} {defect:.12f} {regime}" for theta, defect, regime in values)
    (ROOT / "xy_curve.dat").write_text("\n".join(lines) + "\n", encoding="ascii")


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
{''.join(paths)}
<text x="{left + plot_width/2:.2f}" y="{height-10}" text-anchor="middle">interaction angle θ</text>
<text x="18" y="{top + plot_height/2:.2f}" text-anchor="middle"
 transform="rotate(-90 18 {top + plot_height/2:.2f})">unavoidable error δₐ</text>
<text x="{x(0.12):.2f}" y="{y(0.25):.2f}" fill="{colors['weak']}">weak</text>
<text x="{x(0.42):.2f}" y="{y(1.22):.2f}" fill="{colors['middle']}">quartic</text>
<text x="{x(0.68):.2f}" y="{y(1.42):.2f}" fill="{colors['strong']}">iSWAP</text>
</g>
</svg>
"""
    (ROOT / "xy_curve.svg").write_text(svg, encoding="utf-8")


def main() -> None:
    values = points()
    write_data(values)
    write_svg(values)


if __name__ == "__main__":
    main()
