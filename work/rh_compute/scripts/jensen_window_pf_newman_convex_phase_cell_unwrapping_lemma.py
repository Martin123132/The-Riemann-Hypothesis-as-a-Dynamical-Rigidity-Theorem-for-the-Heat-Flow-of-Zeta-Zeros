#!/usr/bin/env python3
"""Build the convex-cell phase-unwrapping and exact winding lemma."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Iterable, Sequence

import jensen_window_pf_newman_q208_selected_boundary_pilot as q208


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_convex_phase_cell_unwrapping_lemma"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
Q208_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_q208_selected_boundary_pilot.json"
)

Point = tuple[Fraction, Fraction]
Rectangle = tuple[Fraction, Fraction, Fraction, Fraction]


@dataclass(frozen=True)
class LemmaRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _validate_rectangle(rectangle: Rectangle) -> None:
    x_low, x_high, y_low, y_high = rectangle
    if x_low > x_high or y_low > y_high:
        raise ValueError("invalid rectangle bounds")


def rectangle_contains_origin(rectangle: Rectangle) -> bool:
    _validate_rectangle(rectangle)
    x_low, x_high, y_low, y_high = rectangle
    return x_low <= 0 <= x_high and y_low <= 0 <= y_high


def rectangle_contains_point(rectangle: Rectangle, point: Point) -> bool:
    _validate_rectangle(rectangle)
    x_low, x_high, y_low, y_high = rectangle
    x, y = point
    return x_low <= x <= x_high and y_low <= y <= y_high


def project_zero_to_rectangle(rectangle: Rectangle) -> Point:
    _validate_rectangle(rectangle)

    def project_coordinate(low: Fraction, high: Fraction) -> Fraction:
        if low > 0:
            return low
        if high < 0:
            return high
        return Fraction(0)

    x_low, x_high, y_low, y_high = rectangle
    return (
        project_coordinate(x_low, x_high),
        project_coordinate(y_low, y_high),
    )


def dot(left: Point, right: Point) -> Fraction:
    return left[0] * right[0] + left[1] * right[1]


def cross(left: Point, right: Point) -> Fraction:
    return left[0] * right[1] - left[1] * right[0]


def edge_contains_origin(left: Point, right: Point) -> bool:
    if left == (0, 0) or right == (0, 0):
        return True
    return cross(left, right) == 0 and dot(left, right) <= 0


def polygon_winding(vertices: Sequence[Point]) -> int:
    """Return the exact winding around zero using a half-open positive ray."""
    if len(vertices) < 3:
        raise ValueError("a closed polygon needs at least three vertices")
    winding = 0
    for index, left in enumerate(vertices):
        right = vertices[(index + 1) % len(vertices)]
        if edge_contains_origin(left, right):
            raise ValueError("polygon edge meets the origin")
        determinant = cross(left, right)
        if left[1] <= 0 < right[1] and determinant > 0:
            winding += 1
        elif right[1] <= 0 < left[1] and determinant < 0:
            winding -= 1
    return winding


def certified_polygon_winding(
    cells: Sequence[Rectangle],
    intersection_witnesses: Sequence[Point],
) -> int:
    """Validate a cyclic convex-cell chain and wind its witness polygon."""
    if len(cells) != len(intersection_witnesses):
        raise ValueError("cell and witness counts differ")
    if len(cells) < 3:
        raise ValueError("a cyclic cell chain needs at least three cells")
    for index, cell in enumerate(cells):
        if rectangle_contains_origin(cell):
            raise ValueError(f"cell {index} contains the origin")
        start = intersection_witnesses[index]
        end = intersection_witnesses[(index + 1) % len(cells)]
        previous = cells[(index - 1) % len(cells)]
        if not rectangle_contains_point(previous, start):
            raise ValueError(f"witness {index} misses the previous cell")
        if not rectangle_contains_point(cell, start):
            raise ValueError(f"witness {index} misses the next cell")
        if not rectangle_contains_point(cell, end):
            raise ValueError(f"witness {index + 1} misses cell {index}")
    return polygon_winding(intersection_witnesses)


def edge_rectangles(vertices: Sequence[Point]) -> list[Rectangle]:
    rectangles: list[Rectangle] = []
    for index, left in enumerate(vertices):
        right = vertices[(index + 1) % len(vertices)]
        rectangles.append(
            (
                min(left[0], right[0]),
                max(left[0], right[0]),
                min(left[1], right[1]),
                max(left[1], right[1]),
            )
        )
    return rectangles


def exact_polygon_audits() -> list[dict]:
    polygons: list[tuple[str, list[Point], int]] = [
        (
            "counterclockwise_square",
            [
                (Fraction(1), Fraction(-1)),
                (Fraction(1), Fraction(1)),
                (Fraction(-1), Fraction(1)),
                (Fraction(-1), Fraction(-1)),
            ],
            1,
        ),
        (
            "clockwise_square",
            [
                (Fraction(1), Fraction(-1)),
                (Fraction(-1), Fraction(-1)),
                (Fraction(-1), Fraction(1)),
                (Fraction(1), Fraction(1)),
            ],
            -1,
        ),
        (
            "upper_half_plane_triangle",
            [
                (Fraction(2), Fraction(1)),
                (Fraction(-2), Fraction(1)),
                (Fraction(0), Fraction(3)),
            ],
            0,
        ),
    ]
    audits: list[dict] = []
    for name, vertices, expected in polygons:
        cells = edge_rectangles(vertices)
        actual = certified_polygon_winding(cells, vertices)
        if actual != expected:
            raise RuntimeError(
                f"{name} winding mismatch: expected {expected}, got {actual}"
            )
        audits.append(
            {
                "name": name,
                "vertices": [
                    [str(coordinate) for coordinate in point]
                    for point in vertices
                ],
                "cell_count": len(cells),
                "winding": actual,
            }
        )
    return audits


def rectangle_projection_audits() -> list[dict]:
    fixtures: list[Rectangle] = [
        (Fraction(2), Fraction(5), Fraction(-3), Fraction(4)),
        (Fraction(-5), Fraction(-2), Fraction(-4), Fraction(6)),
        (Fraction(-3), Fraction(7), Fraction(2), Fraction(5)),
        (Fraction(-4), Fraction(8), Fraction(-6), Fraction(-1)),
        (Fraction(2), Fraction(5), Fraction(3), Fraction(7)),
    ]
    audits: list[dict] = []
    for rectangle in fixtures:
        projection = project_zero_to_rectangle(rectangle)
        norm_squared = dot(projection, projection)
        if norm_squared <= 0:
            raise RuntimeError("origin-free fixture has zero projection")
        x_low, x_high, y_low, y_high = rectangle
        vertices = [
            (x_low, y_low),
            (x_low, y_high),
            (x_high, y_low),
            (x_high, y_high),
        ]
        margins = [dot(projection, vertex) - norm_squared for vertex in vertices]
        if min(margins) < 0:
            raise RuntimeError("rectangle projection separation failed")
        audits.append(
            {
                "rectangle": [str(value) for value in rectangle],
                "projection": [str(value) for value in projection],
                "projection_norm_squared": str(norm_squared),
                "minimum_vertex_margin": str(min(margins)),
            }
        )
    return audits


def _fraction(value: str) -> Fraction:
    return Fraction(value)


def _full_derivative_lower(leaf: dict):
    arb = q208.bridge.compact.arb
    return (
        arb(leaf["j_retained_prime_lower"]).lower()
        - arb(leaf["tail_derivative_upper"]).upper()
    )


def _full_interval(leaf: dict, stem: str, tail: str):
    arb = q208.bridge.compact.arb
    retained = arb(leaf[f"{stem}_lower"]).union(
        arb(leaf[f"{stem}_upper"])
    )
    radius = arb(leaf[tail]).upper()
    return retained + arb(0, radius.str(120))


def _contiguous_time_cover(leaves: Iterable[dict]) -> tuple[bool, str, str]:
    intervals = sorted(
        (
            _fraction(leaf["t_low"]),
            _fraction(leaf["t_high"]),
        )
        for leaf in leaves
    )
    if not intervals:
        return False, "", ""
    contiguous = all(
        left[1] == right[0]
        for left, right in zip(intervals, intervals[1:])
    )
    return contiguous, str(intervals[0][0]), str(intervals[-1][1])


def q208_right_edge_audit() -> dict:
    q208.bridge.compact.flint.ctx.prec = 352
    stored = json.loads(Q208_RESULT.read_text(encoding="utf-8"))
    right_records = [
        record
        for record in stored["records"]
        if record["task"]["region"] == "full_time_right_strip"
    ]
    strip_leaves = [
        leaf
        for record in right_records
        for leaf in record["result"]["certified_leaves"]
    ]
    edge_records = [
        record
        for record in right_records
        if record["task"]["x_high"] == "246"
    ]
    edge_leaves = [
        leaf
        for record in edge_records
        for leaf in record["result"]["certified_leaves"]
    ]
    strip_lowers = [_full_derivative_lower(leaf) for leaf in strip_leaves]
    edge_lowers = [_full_derivative_lower(leaf) for leaf in edge_leaves]
    if len(strip_leaves) != 40 or len(edge_leaves) != 20:
        raise RuntimeError("Q208 right-strip or right-edge cell count drifted")
    if not all(lower > 0 for lower in strip_lowers):
        raise RuntimeError("Q208 right strip is not derivative-positive")
    contiguous, time_low, time_high = _contiguous_time_cover(edge_leaves)
    if not contiguous or time_low != "1/1040" or time_high != "1/5":
        raise RuntimeError("Q208 right-edge time cover is not contiguous")
    arb = q208.bridge.compact.arb
    x = arb(246)
    scale = 1 + 1 / x**4
    shear = 4 / x**5
    proxy_derivative_lowers = []
    for leaf in edge_leaves:
        j_box = _full_interval(
            leaf,
            "j_retained",
            "tail_value_upper",
        )
        j_prime_box = _full_interval(
            leaf,
            "j_retained_prime",
            "tail_derivative_upper",
        )
        proxy_derivative_lowers.append(
            (scale * j_prime_box - shear * j_box).lower()
        )
    if not all(lower > 0 for lower in proxy_derivative_lowers):
        raise RuntimeError(
            "axis-safe proxy derivative is not positive on the right edge"
        )
    return {
        "source_kind": stored["kind"],
        "source_sha256": file_hash(Q208_RESULT),
        "right_strip_domain": "[1/1040,1/5]x[245,246]",
        "right_strip_cells": len(strip_leaves),
        "right_edge_domain": "x=246, 1/1040<=t<=1/5",
        "right_edge_cells": len(edge_leaves),
        "right_edge_time_cover_contiguous": contiguous,
        "right_edge_time_low": time_low,
        "right_edge_time_high": time_high,
        "minimum_strip_full_derivative_lower": str(min(strip_lowers)),
        "minimum_edge_full_derivative_lower": str(min(edge_lowers)),
        "proxy": "F_t(x)=16*(1+x^4)*H_t(x)",
        "right_edge_transfer": (
            "F=(1+x^-4)*J and "
            "F'=(1+x^-4)*J'-4*x^-5*J at x=246"
        ),
        "minimum_edge_proxy_derivative_lower": str(
            min(proxy_derivative_lowers)
        ),
        "phase_half_plane": "Im(F_t(246)+i*F_t'(246))>0",
        "phase_conclusion": (
            "The transformed axis-safe right-edge path has a single "
            "principal argument in (0,pi) and cannot contribute a complete "
            "turn by itself."
        ),
    }


def build_exact() -> dict:
    return {
        "projection_separation": (
            "For a compact convex C with 0 not in C, let p be the closest "
            "point of C to 0. Then p dot v >= ||p||^2>0 for every v in C. "
            "Thus C lies in one open half-plane through 0 and admits a "
            "single argument branch of width strictly less than pi."
        ),
        "rectangle_projection": (
            "For C=[a,b]x[c,d], p=(proj_[a,b](0),proj_[c,d](0)); C excludes "
            "0 iff p!=0, and p dot v>=||p||^2 on all of C."
        ),
        "cell_to_chord_homotopy": (
            "If gamma([s_k,s_(k+1)]) is contained in convex C_k and "
            "v_k=gamma(s_k), then (1-u)gamma(s)+u*((1-r)v_k+r*v_(k+1)) "
            "is a homotopy relative endpoints from gamma to its endpoint "
            "chord inside C_k, where r is the affine edge parameter."
        ),
        "cyclic_witness_polygon": (
            "If q_k belongs to C_(k-1) intersect C_k cyclically, the exact "
            "endpoint polygon is homotopic to the rational witness polygon "
            "q_0...q_(n-1) edge by edge inside C_k. Therefore "
            "wind(gamma,0)=wind(q_0...q_(n-1),0)."
        ),
        "exact_ray_crossing": (
            "For each oriented edge p->q, add +1 when p_y<=0<q_y and "
            "det(p,q)>0; add -1 when q_y<=0<p_y and det(p,q)<0. The sum is "
            "the exact winding, provided no edge meets 0."
        ),
        "phase_unwrapping_equivalence": (
            "The closest-point half-plane gives a local argument interval "
            "of width <pi on every C_k. A certified intersection witness "
            "selects the unique 2pi translate shared by adjacent branches. "
            "The resulting lifted phase has the same integer increment as "
            "the exact witness-polygon crossing count."
        ),
        "interval_certificate_contract": (
            "Every cell must enclose the whole continuous path segment, not "
            "sampled points; every cell must exclude 0; each cyclic adjacent "
            "pair must have an exact rational intersection witness; and all "
            "ray-crossing determinants are evaluated exactly. Failure of "
            "any condition is unresolved, never an inferred phase."
        ),
        "axis_safe_proxy": (
            "Set F_t(x)=16*(1+x^4)*H_t(x). Its first jet is "
            "T_x(H_t,H_t') with T_x=16[[1+x^4,0],[4x^3,1+x^4]]. "
            "The homotopy T_(x,s)=[[1-s+16s(1+x^4),0],"
            "[64s*x^3,1-s+16s(1+x^4)]] has positive determinant for "
            "0<=s<=1. Hence F and H have the same contacts and the same "
            "closed-boundary winding, while F is regular at x=0 and its "
            "right-edge jet is rigorously transferred from J."
        ),
        "right_edge_calibration": q208_right_edge_audit(),
        "q208_handoff": (
            "Build a complete fixed-time first-jet cell chain on "
            "t=1/1040, 0<=x<=246, store cyclic endpoint witnesses with the "
            "axis, top, and right edges, and compute one exact closed "
            "witness-polygon winding. Until that chain exists, Q208, "
            "Lambda<=0, RH, and the Clay prize remain open."
        ),
    }


def build_rows(exact: dict) -> list[LemmaRow]:
    return [
        LemmaRow(
            "npcu_01_projection_separation",
            "exact_convex_geometry_lemma",
            "ready_to_apply",
            "Every compact convex origin-free phase cell has a strict separating half-plane.",
            exact["projection_separation"],
            "No smoothness or Xi input is used.",
            rectangle_projection_audits(),
        ),
        LemmaRow(
            "npcu_02_rectangle_projection",
            "exact_interval_geometry",
            "ready_to_apply",
            "Axis-aligned interval cells have an explicit exact separation witness.",
            exact["rectangle_projection"],
            "Endpoint bounds must be outward enclosures.",
            rectangle_projection_audits(),
        ),
        LemmaRow(
            "npcu_03_cell_to_chord",
            "exact_homotopy_lemma",
            "ready_to_apply",
            "A path segment inside one convex origin-free cell may be replaced by its endpoint chord.",
            exact["cell_to_chord_homotopy"],
            "Convexity keeps the full homotopy inside the certified cell.",
        ),
        LemmaRow(
            "npcu_04_witness_polygon",
            "exact_homotopy_composition",
            "ready_to_apply",
            "Cyclic cell intersections reduce continuous-path winding to a rational witness polygon.",
            exact["cyclic_witness_polygon"],
            "Each witness must lie in both adjacent certified cells.",
        ),
        LemmaRow(
            "npcu_05_exact_crossing_count",
            "exact_integer_algorithm",
            "ready_to_apply",
            "The witness-polygon winding is computed without floating-point angles.",
            exact["exact_ray_crossing"],
            "Edges meeting the origin are rejected.",
            exact_polygon_audits(),
        ),
        LemmaRow(
            "npcu_06_phase_unwrapping",
            "exact_phase_lift_equivalence",
            "ready_to_apply",
            "Sequential local argument unwrapping and exact polygon winding are equivalent certificates.",
            exact["phase_unwrapping_equivalence"],
            "The polygon formulation avoids branch-cut and near-pi ambiguity.",
        ),
        LemmaRow(
            "npcu_07_interval_contract",
            "rigorous_numerical_contract",
            "ready_to_apply",
            "The method promotes interval path covers only when all geometric gates are explicit.",
            exact["interval_certificate_contract"],
            "Point sampling alone cannot satisfy this contract.",
        ),
        LemmaRow(
            "npcu_08_axis_safe_proxy",
            "exact_orientation_preserving_proxy",
            "ready_to_apply",
            "One globally regular triangular proxy joins the axis to the high-frequency J cells.",
            exact["axis_safe_proxy"],
            "The positive-determinant homotopy preserves every closed winding integer.",
        ),
        LemmaRow(
            "npcu_09_q208_right_edge",
            "rigorous_source_composition",
            "ready_to_apply",
            "The complete transformed Q208 right edge lies in the open upper proxy half-plane.",
            exact["right_edge_calibration"]["phase_conclusion"],
            "This fixes a phase branch on one open edge, not the closed-boundary winding.",
            exact["right_edge_calibration"],
        ),
        LemmaRow(
            "npcu_10_q208_handoff",
            "open_theorem_target",
            "not_ready_to_apply",
            "The complete Q208 bottom cell chain and closed winding remain to be certified.",
            exact["q208_handoff"],
            "No Q208, cofinal, Lambda<=0, RH, or prize conclusion is promoted.",
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    rows = build_rows(exact)
    return {
        "kind": STEM,
        "date": "2026-07-25",
        "status": (
            "exact convex-cell phase-unwrapping and rational-polygon winding "
            "lemma with rigorous Q208 right-edge calibration; the complete "
            "Q208 bottom chain and closed winding remain open"
        ),
        "proof_boundary": (
            "This artifact proves a reusable topological and exact-arithmetic "
            "certificate and applies its half-plane branch to the stored "
            "Q208 right edge. It does not certify the complete bottom edge, "
            "the closed Q208 boundary winding, Q208, Lambda<=0, RH, or the "
            "Clay prize."
        ),
        "sources": [
            str(Q208_RESULT.relative_to(REPO_ROOT)).replace("\\", "/"),
            (
                "work/rh_compute/scripts/"
                "jensen_window_pf_newman_q208_selected_boundary_pilot.py"
            ),
        ],
        "source_sha256": {
            "q208_result": file_hash(Q208_RESULT),
            "q208_builder": file_hash(Path(q208.__file__).resolve()),
        },
        "exact": exact,
        "projection_audits": rectangle_projection_audits(),
        "polygon_audits": exact_polygon_audits(),
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    right = artifact["exact"]["right_edge_calibration"]
    return "\n".join(
        [
            "# Newman Convex Phase-Cell Unwrapping Lemma",
            "",
            "Date: 2026-07-25",
            "",
            "Status: exact reusable phase/winding certificate with a rigorous",
            "Q208 right-edge calibration. This is not a proof of Q208,",
            "`Lambda<=0`, RH, or the Clay prize.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Convex Cell Lemma",
            "",
            "Let a closed continuous path be split into finitely many pieces,",
            "with the kth piece enclosed by a compact convex set `C_k` that",
            "does not contain the origin. If `p_k` is the closest point of",
            "`C_k` to the origin, then",
            "",
            "```text",
            "p_k dot v >= ||p_k||^2 > 0  for every v in C_k.",
            "```",
            "",
            "Thus every cell lies in a strict open half-plane through the",
            "origin and has one argument branch of width less than `pi`.",
            "",
            "Convexity gives a relative-endpoint homotopy from each path",
            "piece to its endpoint chord. If an exact rational point `q_k`",
            "lies in `C_(k-1) intersect C_k` for every cyclic join, those",
            "chords are in turn homotopic to the rational polygon through",
            "the `q_k`. Therefore the continuous path and witness polygon",
            "have exactly the same winding.",
            "",
            "## Exact Winding",
            "",
            "The polygon integer is evaluated by signed crossings of the",
            "positive real ray:",
            "",
            "```text",
            "+1: p_y <= 0 < q_y and det(p,q) > 0",
            "-1: q_y <= 0 < p_y and det(p,q) < 0.",
            "```",
            "",
            "All coordinates and determinants are rational. No floating-point",
            "`atan2`, branch cut, or tolerance enters the winding integer.",
            "A cell containing the origin, a missing adjacent-cell witness,",
            "or an edge through the origin is rejected as unresolved.",
            "",
            "The checker exercises exact winding `+1`, `-1`, and `0` examples",
            "and independent rejection cases.",
            "",
            "## Q208 Right Edge",
            "",
            "The stored six-term Arb/Taylor theorem has been replayed at the",
            "geometric interface and transferred to the globally regular",
            "proxy `F=16(1+x^4)H`:",
            "",
            "```text",
            f"right strip: {right['right_strip_domain']}",
            f"strip cells: {right['right_strip_cells']}",
            f"right edge: {right['right_edge_domain']}",
            f"edge cells: {right['right_edge_cells']}",
            (
                "minimum full right-edge derivative lower bound: "
                f"{right['minimum_edge_full_derivative_lower']}"
            ),
            (
                "minimum transformed proxy derivative lower bound: "
                f"{right['minimum_edge_proxy_derivative_lower']}"
            ),
            "```",
            "",
            "Every right-edge cell has both `J_t'(246)>0` and",
            "`F_t'(246)>0`, and the 20 time cells form a contiguous cover",
            "from `1/1040` to `1/5`. Hence `F_t(246)+i F_t'(246)` stays in",
            "the open upper half-plane and admits one principal phase branch.",
            "This calibrates the method but does not determine the complete",
            "closed-boundary winding.",
            "",
            "## Live Handoff",
            "",
            "The next certificate must cover the entire fixed-time bottom",
            "edge `t=1/1040`, `0<=x<=246` by origin-free first-jet cells,",
            "store exact rational witnesses at every join, compose them with",
            "the axis, top, and right pieces, and compute one exact cyclic",
            "polygon winding. Until then Q208 remains open.",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman convex phase-cell unwrapping lemma: "
        f"{len(artifact['rows'])} rows, "
        f"{len(artifact['projection_audits'])} projection audits, "
        f"{len(artifact['polygon_audits'])} exact polygon audits, "
        f"{artifact['exact']['right_edge_calibration']['right_edge_cells']} "
        "Q208 right-edge cells, 1 open bottom-edge target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
