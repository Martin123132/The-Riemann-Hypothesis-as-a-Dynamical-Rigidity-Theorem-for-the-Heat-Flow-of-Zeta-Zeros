#!/usr/bin/env python3
"""Compose the complete Q208 boundary and certify its exact winding."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import jensen_window_pf_newman_convex_phase_cell_unwrapping_lemma as phase
import jensen_window_pf_newman_q208_bottom_phase_cell_certificate as bottom
import jensen_window_pf_newman_q208_top_phase_cell_certificate as top


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_q208_closed_boundary_winding_certificate"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
BOTTOM_RESULT = bottom.DEFAULT_OUT
TOP_RESULT = top.DEFAULT_OUT
RIGHT_RESULT = phase.Q208_RESULT
PHASE_RESULT = phase.DEFAULT_OUT
WINDING_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_first_jet_winding_gate.json"
)
BOUNDARY_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_positive_boundary_attainment_lemma.json"
)

compact = bottom.compact
Cell = tuple[object, object]
Point = tuple[Fraction, Fraction]


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_complete(path: Path, flag: str, label: str) -> dict:
    if not path.exists():
        raise RuntimeError(f"{label} result is missing")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("summary", {}).get(flag) is not True:
        raise RuntimeError(f"{label} result is not complete")
    return payload


def horizontal_cells(payload: dict) -> list[Cell]:
    leaves = [
        leaf
        for record in payload["records"]
        for leaf in record["result"]["certified_leaves"]
    ]
    leaves.sort(key=lambda row: Fraction(row["x_low"]))
    if not leaves:
        raise RuntimeError("horizontal phase chain is empty")
    if Fraction(leaves[0]["x_low"]) != bottom.X_LOWER:
        raise RuntimeError("horizontal chain misses x=0")
    if Fraction(leaves[-1]["x_high"]) != bottom.X_UPPER:
        raise RuntimeError("horizontal chain misses x=246")
    for left, right in zip(leaves, leaves[1:]):
        if Fraction(left["x_high"]) != Fraction(right["x_low"]):
            raise RuntimeError("horizontal phase chain has a gap")
    return [
        (
            compact.arb(leaf["full_f"]),
            compact.arb(leaf["full_f_prime"]),
        )
        for leaf in leaves
    ]


def retained_interval(leaf: dict, stem: str, tail: str):
    retained = compact.arb(leaf[f"{stem}_lower"]).union(
        compact.arb(leaf[f"{stem}_upper"])
    )
    radius = compact.arb(leaf[tail]).upper()
    return retained + compact.arb(0, radius.str(120))


def transformed_right_cells() -> tuple[list[Cell], dict]:
    payload = json.loads(RIGHT_RESULT.read_text(encoding="utf-8"))
    records = [
        record
        for record in payload["records"]
        if record["task"]["region"] == "full_time_right_strip"
        and record["task"]["x_high"] == "246"
    ]
    leaves = [
        leaf
        for record in records
        for leaf in record["result"]["certified_leaves"]
    ]
    leaves.sort(key=lambda row: Fraction(row["t_low"]))
    if len(leaves) != 20:
        raise RuntimeError("right edge does not have 20 cells")
    if Fraction(leaves[0]["t_low"]) != Fraction(1, 1040):
        raise RuntimeError("right edge misses t=1/1040")
    if Fraction(leaves[-1]["t_high"]) != Fraction(1, 5):
        raise RuntimeError("right edge misses t=1/5")
    for left, right in zip(leaves, leaves[1:]):
        if Fraction(left["t_high"]) != Fraction(right["t_low"]):
            raise RuntimeError("right-edge time cells have a gap")

    x = compact.arb(246)
    scale = 1 + 1 / x**4
    shear = 4 / x**5
    cells: list[Cell] = []
    derivative_lowers = []
    for leaf in leaves:
        j_box = retained_interval(
            leaf,
            "j_retained",
            "tail_value_upper",
        )
        j_prime_box = retained_interval(
            leaf,
            "j_retained_prime",
            "tail_derivative_upper",
        )
        f_box = scale * j_box
        f_prime_box = scale * j_prime_box - shear * j_box
        if f_prime_box.lower() <= 0:
            raise RuntimeError(
                "transformed right-edge cell leaves the upper half-plane"
            )
        derivative_lowers.append(f_prime_box.lower())
        cells.append((f_box, f_prime_box))
    return cells, {
        "cell_count": len(cells),
        "time_domain": ["1/1040", "1/5"],
        "minimum_proxy_derivative_lower": str(min(derivative_lowers)),
    }


def cell_contains_origin(cell: Cell) -> bool:
    return cell[0].contains(compact.arb(0)) and cell[1].contains(
        compact.arb(0)
    )


def intersection_witness(previous: Cell, following: Cell) -> Point:
    try:
        real = previous[0].intersection(following[0])
        imaginary = previous[1].intersection(following[1])
    except ValueError as exc:
        raise RuntimeError("adjacent boundary cells do not intersect") from exc
    witness = (
        bottom.interval_midpoint_fraction(real),
        bottom.interval_midpoint_fraction(imaginary),
    )
    if not bottom.point_in_cell(witness, previous):
        raise RuntimeError("join witness misses its previous cell")
    if not bottom.point_in_cell(witness, following):
        raise RuntimeError("join witness misses its following cell")
    return witness


def compose_cells(
    bottom_cells: list[Cell],
    right_cells: list[Cell],
    top_cells_forward: list[Cell],
) -> tuple[list[Cell], list[str]]:
    top_cells_reverse = list(reversed(top_cells_forward))
    bottom_first = bottom_cells[0]
    top_first = top_cells_forward[0]
    if bottom_first[0].lower() <= 0 or top_first[0].lower() <= 0:
        raise RuntimeError("axis endpoint cell is not value-positive")
    if not bottom_first[1].contains(compact.arb(0)):
        raise RuntimeError("bottom axis cell misses F'=0")
    if not top_first[1].contains(compact.arb(0)):
        raise RuntimeError("top axis cell misses F'=0")
    axis_value = bottom_first[0].union(top_first[0])
    if axis_value.lower() <= 0:
        raise RuntimeError("axis hull reaches the origin")
    axis_cell: Cell = (axis_value, compact.arb(0))
    cells = [
        *bottom_cells,
        *right_cells,
        *top_cells_reverse,
        axis_cell,
    ]
    labels = [
        *[f"bottom:{index}" for index in range(len(bottom_cells))],
        *[f"right:{index}" for index in range(len(right_cells))],
        *[
            f"top_reverse:{index}"
            for index in range(len(top_cells_reverse))
        ],
        "axis",
    ]
    if any(cell_contains_origin(cell) for cell in cells):
        raise RuntimeError("composed boundary has an origin-containing cell")
    return cells, labels


def build_certificate() -> dict:
    compact.flint.ctx.prec = bottom.bridge.PRECISION_BITS
    bottom_payload = load_complete(
        BOTTOM_RESULT,
        "complete_bottom_edge_certificate",
        "bottom",
    )
    top_payload = load_complete(
        TOP_RESULT,
        "complete_top_edge_certificate",
        "top",
    )
    bottom_cells = horizontal_cells(bottom_payload)
    top_cells = horizontal_cells(top_payload)
    right_cells, right_audit = transformed_right_cells()
    cells, labels = compose_cells(
        bottom_cells,
        right_cells,
        top_cells,
    )
    witnesses = [
        intersection_witness(cells[index - 1], cell)
        for index, cell in enumerate(cells)
    ]
    winding = phase.polygon_winding(witnesses)
    if winding != 0:
        raise RuntimeError(
            f"complete Q208 witness polygon has winding {winding}, not zero"
        )
    edge_types = {
        "bottom": len(bottom_cells),
        "right": len(right_cells),
        "top_reverse": len(top_cells),
        "axis": 1,
    }
    return {
        "kind": STEM,
        "date": "2026-07-26",
        "status": (
            "rigorous finite Q208 closed-boundary zero-winding and "
            "no-contact certificate"
        ),
        "proof_boundary": (
            "The exact convex-cell composition certifies the finite rectangle "
            "[1/1040,1/5]x[0,246]. The published Lambda<=1/5 bound excludes "
            "contacts for 1/5<t<=1/4, so Q208 is contact-free and has zero "
            "first-jet winding. This proves no stage j>=209, no cofinal "
            "termination theorem, no Lambda<=0, no RH, and no Clay prize."
        ),
        "sources": {
            "bottom": str(BOTTOM_RESULT.relative_to(REPO_ROOT)).replace(
                "\\", "/"
            ),
            "top": str(TOP_RESULT.relative_to(REPO_ROOT)).replace("\\", "/"),
            "right": str(RIGHT_RESULT.relative_to(REPO_ROOT)).replace(
                "\\", "/"
            ),
            "phase": str(PHASE_RESULT.relative_to(REPO_ROOT)).replace(
                "\\", "/"
            ),
            "winding": str(WINDING_RESULT.relative_to(REPO_ROOT)).replace(
                "\\", "/"
            ),
            "positive_boundary": str(
                BOUNDARY_RESULT.relative_to(REPO_ROOT)
            ).replace("\\", "/"),
        },
        "source_sha256": {
            "bottom": file_hash(BOTTOM_RESULT),
            "top": file_hash(TOP_RESULT),
            "right": file_hash(RIGHT_RESULT),
            "phase": file_hash(PHASE_RESULT),
            "winding": file_hash(WINDING_RESULT),
            "positive_boundary": file_hash(BOUNDARY_RESULT),
        },
        "proxy": {
            "definition": "F_t(x)=16*(1+x^4)*H_t(x)",
            "first_jet_homotopy": (
                "T_(x,s)=[[1-s+16s(1+x^4),0],"
                "[64s*x^3,1-s+16s(1+x^4)]] in GL+(2,R)"
            ),
            "contact_and_winding_transfer": True,
        },
        "domain": {
            "certified_low_rectangle": "[1/1040,1/5]x[0,246]",
            "q208": "[1/1040,1/4]x[0,246]",
            "high_time_composition": (
                "Lambda<=1/5 gives no contact for 1/5<t<=1/4"
            ),
        },
        "cell_counts": edge_types,
        "total_cells": len(cells),
        "total_witnesses": len(witnesses),
        "right_edge_audit": right_audit,
        "witness_cell_labels": labels,
        "witnesses": [
            [str(coordinate) for coordinate in point]
            for point in witnesses
        ],
        "exact_winding": winding,
        "conclusion": (
            "Q_208 is contact-free and has zero first-jet winding; "
            "therefore Q_1 through Q_208 are certified by containment."
        ),
    }


def render_note(artifact: dict) -> str:
    counts = artifact["cell_counts"]
    return "\n".join(
        [
            "# Newman Q208 Closed-Boundary Winding Certificate",
            "",
            "Date: 2026-07-26",
            "",
            "Status: rigorous finite Q208 zero-winding and no-contact",
            "certificate. This is not a proof of `Lambda<=0`, RH, or the",
            "Clay prize.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Exact Boundary",
            "",
            "The positively oriented low rectangle is composed as",
            "",
            "```text",
            "bottom: t=1/1040, x=0 -> 246",
            "right:  x=246, t=1/1040 -> 1/5",
            "top:    t=1/5, x=246 -> 0",
            "axis:   x=0, t=1/5 -> 1/1040.",
            "```",
            "",
            "Every numerical cell encloses a whole path arc and excludes the",
            "origin. The axis image is positive real because",
            "`F_t(0)=16H_t(0)>0`, `F_t'(0)=0`, and `H_t(0)` increases with",
            "`t`. Exact dyadic witnesses lie in every cyclic adjacent-cell",
            "intersection.",
            "",
            "```text",
            f"bottom cells={counts['bottom']}",
            f"right cells={counts['right']}",
            f"reversed top cells={counts['top_reverse']}",
            f"axis cells={counts['axis']}",
            f"total cells={artifact['total_cells']}",
            f"exact witness-polygon winding={artifact['exact_winding']}",
            "```",
            "",
            "The proxy `F=16(1+x^4)H` is joined to the ordinary first jet by",
            "a positive-determinant triangular homotopy. Zero proxy winding",
            "therefore equals zero `H+iH_x` winding. Since every heat contact",
            "has strictly positive local index in the standard `(x,t)`",
            "orientation, the low rectangle contains",
            "no contact.",
            "",
            "The independent `Lambda<=1/5` theorem excludes contacts for",
            "`1/5<t<=1/4`; the explicitly certified top handles equality.",
            "Thus",
            "",
            "```text",
            artifact["conclusion"],
            "```",
            "",
            "## Scope",
            "",
            "This is one finite successor to Q207. It supplies no stage",
            "`j>=209`, no parameter-uniform termination theorem, no",
            "`Lambda<=0`, and no RH conclusion.",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_certificate()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Q208 closed-boundary winding certificate: "
        f"{artifact['total_cells']} cells, "
        f"{artifact['total_witnesses']} exact witnesses, "
        f"winding={artifact['exact_winding']}, "
        "Q1..Q208 certified, 0 cofinal promotions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
