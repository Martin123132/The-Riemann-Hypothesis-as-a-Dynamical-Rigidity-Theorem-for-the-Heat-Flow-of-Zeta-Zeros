#!/usr/bin/env python3
"""Compose the Q207-to-Q208 successor from Q207 old-edge cells."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import jensen_window_pf_newman_theta_forward_six_term_finite_bridge_interval_certificate as bridge


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_q207_q208_forward_adiabatic_successor_certificate"
RESULT_DIR = REPO_ROOT / "work/rh_compute/results"
DEFAULT_OUT = RESULT_DIR / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"

Q207_STEM = (
    "jensen_window_pf_newman_theta_forward_six_term_"
    "finite_bridge_interval_certificate"
)
COARSE_STEM = (
    "jensen_window_pf_newman_q207_q208_adiabatic_"
    "bottom_collar_certificate"
)
REFINED_STEM = (
    "jensen_window_pf_newman_q207_q208_adiabatic_"
    "bottom_collar_refined_tail_certificate"
)
CORE_STEM = (
    "jensen_window_pf_newman_theta_compact_"
    "transversality_interval_certificate"
)
RIGHT_STEM = "jensen_window_pf_newman_q208_selected_boundary_pilot"
SUCCESSOR_STEM = "jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma"

Q207_CACHE = RESULT_DIR / f"{Q207_STEM}.jsonl"
Q207_RESULT = RESULT_DIR / f"{Q207_STEM}.json"
COARSE_CACHE = RESULT_DIR / f"{COARSE_STEM}.jsonl"
COARSE_RESULT = RESULT_DIR / f"{COARSE_STEM}.json"
REFINED_CACHE = RESULT_DIR / f"{REFINED_STEM}.jsonl"
REFINED_RESULT = RESULT_DIR / f"{REFINED_STEM}.json"
CORE_RESULT = RESULT_DIR / f"{CORE_STEM}.json"
RIGHT_RESULT = RESULT_DIR / f"{RIGHT_STEM}.json"
SUCCESSOR_RESULT = RESULT_DIR / f"{SUCCESSOR_STEM}.json"

DATE = "2026-07-25"
PRECISION_BITS = 352
ENDPOINT_DIGITS = 105
OLD_TIME = Fraction(1, 1035)
NEW_TIME = Fraction(1, 1040)
TIME_STEP = OLD_TIME - NEW_TIME
CORE_END = Fraction(38)
OLD_RADIUS = Fraction(245)
NEW_RADIUS = Fraction(246)
COARSE_END = Fraction(379, 2)

compact = bridge.compact


def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_jsonl(path: Path) -> list[dict]:
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def fraction_to_arb(value: Fraction):
    return compact.arb(value.numerator) / value.denominator


def serialize_interval(value) -> str:
    if not value.is_finite():
        raise RuntimeError(f"nonfinite interval: {value}")
    return value.str(ENDPOINT_DIGITS, more=True)


def source_hashes() -> dict[str, str]:
    return {
        "q207_cache": file_hash(Q207_CACHE),
        "q207_result": file_hash(Q207_RESULT),
        "coarse_cache": file_hash(COARSE_CACHE),
        "coarse_result": file_hash(COARSE_RESULT),
        "refined_cache": file_hash(REFINED_CACHE),
        "refined_result": file_hash(REFINED_RESULT),
        "compact_core_result": file_hash(CORE_RESULT),
        "right_strip_result": file_hash(RIGHT_RESULT),
        "successor_lemma": file_hash(SUCCESSOR_RESULT),
    }


def q207_old_bottom_cells() -> dict[tuple[Fraction, Fraction], dict]:
    cells: dict[tuple[Fraction, Fraction], dict] = {}
    for record in load_jsonl(Q207_CACHE):
        leaves = [
            leaf
            for leaf in record["result"]["certified_leaves"]
            if Fraction(leaf["t_low"]) == OLD_TIME
        ]
        if len(leaves) != 1:
            raise RuntimeError(
                "each Q207 panel must have exactly one old-bottom cell"
            )
        leaf = leaves[0]
        key = (Fraction(leaf["x_low"]), Fraction(leaf["x_high"]))
        if key in cells:
            raise RuntimeError(f"duplicate Q207 old-bottom cell {key}")
        cells[key] = leaf
    return cells


def full_j_boxes(leaf: dict):
    value_tail = compact.arb(leaf["tail_value_upper"]).upper()
    derivative_tail = compact.arb(
        leaf["tail_derivative_upper"]
    ).upper()
    j_box = compact.arb(leaf["j_retained_lower"]).union(
        compact.arb(leaf["j_retained_upper"])
    ) + compact.arb(0, value_tail.str(120))
    j_prime_box = compact.arb(
        leaf["j_retained_prime_lower"]
    ).union(
        compact.arb(leaf["j_retained_prime_upper"])
    ) + compact.arb(0, derivative_tail.str(120))
    return j_box, j_prime_box


def old_f_cell_distance(leaf: dict):
    """Transform J=16*x^4*H to F=16*(1+x^4)*H."""
    x_low = Fraction(leaf["x_low"])
    x_high = Fraction(leaf["x_high"])
    x_box = fraction_to_arb(x_low).union(fraction_to_arb(x_high))
    j_box, j_prime_box = full_j_boxes(leaf)
    scale = 1 + 1 / x_box**4
    shear = 4 / x_box**5
    f_box = scale * j_box
    f_prime_box = scale * j_prime_box - shear * j_box
    distance = (
        f_box.abs_lower() ** 2 + f_prime_box.abs_lower() ** 2
    ).sqrt().lower()
    if distance <= 0:
        raise RuntimeError(
            f"old Q207 transformed cell meets the origin on "
            f"[{x_low},{x_high}]"
        )
    return f_box, f_prime_box, distance


def enclosing_q207_panel(
    x_low: Fraction,
    x_high: Fraction,
    source: str,
) -> tuple[Fraction, Fraction]:
    if source == "coarse":
        return x_low, x_high
    doubled = 2 * x_low
    parent_index = doubled.numerator // doubled.denominator
    parent_low = Fraction(parent_index, 2)
    parent_high = parent_low + Fraction(1, 2)
    if not (parent_low <= x_low < x_high <= parent_high):
        raise RuntimeError(
            f"refined panel [{x_low},{x_high}] has invalid parent "
            f"[{parent_low},{parent_high}]"
        )
    return parent_low, parent_high


def transport_rows() -> list[tuple[str, dict]]:
    rows: list[tuple[str, dict]] = []
    for record in load_jsonl(COARSE_CACHE):
        x_low = Fraction(record["task"]["x_low"])
        x_high = Fraction(record["task"]["x_high"])
        if CORE_END <= x_low and x_high <= COARSE_END:
            rows.append(("coarse", record))
    for record in load_jsonl(REFINED_CACHE):
        x_low = Fraction(record["task"]["x_low"])
        x_high = Fraction(record["task"]["x_high"])
        if COARSE_END <= x_low and x_high <= OLD_RADIUS:
            rows.append(("refined", record))
    rows.sort(key=lambda item: Fraction(item[1]["task"]["x_low"]))
    return rows


def build_records() -> tuple[list[dict], dict]:
    cells = q207_old_bottom_cells()
    rows = transport_rows()
    records: list[dict] = []
    numeric: list[tuple[object, object]] = []
    previous = CORE_END
    used_parents: set[tuple[Fraction, Fraction]] = set()

    for sequence, (source, record) in enumerate(rows, start=1):
        result = record["result"]
        if result.get("status") != "certified":
            raise RuntimeError(f"uncertified transport row {sequence}")
        x_low = Fraction(record["task"]["x_low"])
        x_high = Fraction(record["task"]["x_high"])
        if x_low != previous:
            raise RuntimeError(
                f"transport cover gap: expected {previous}, got {x_low}"
            )
        previous = x_high

        parent = enclosing_q207_panel(x_low, x_high, source)
        if parent not in cells:
            raise RuntimeError(f"missing Q207 parent cell {parent}")
        used_parents.add(parent)
        f_box, f_prime_box, distance = old_f_cell_distance(cells[parent])
        displacement = compact.arb(
            result["transport"]["displacement_upper"]
        ).upper()
        ratio = (displacement / distance).upper()
        if ratio >= 1:
            raise RuntimeError(
                f"forward transport gate failed on [{x_low},{x_high}]: "
                f"{ratio}"
            )
        numeric.append((ratio, distance))
        records.append(
            {
                "sequence": sequence,
                "source": source,
                "x_low": str(x_low),
                "x_high": str(x_high),
                "q207_parent_x_low": str(parent[0]),
                "q207_parent_x_high": str(parent[1]),
                "q207_old_cell_f": serialize_interval(f_box),
                "q207_old_cell_f_prime": serialize_interval(
                    f_prime_box
                ),
                "q207_old_cell_distance_lower": serialize_interval(
                    distance
                ),
                "transport_displacement_upper": serialize_interval(
                    displacement
                ),
                "transport_ratio_upper": serialize_interval(ratio),
                "strict_gate": True,
            }
        )

    if previous != OLD_RADIUS:
        raise RuntimeError(
            f"transport cover ends at {previous}, not {OLD_RADIUS}"
        )
    if used_parents != set(cells):
        missing = sorted(set(cells) - used_parents)
        extra = sorted(used_parents - set(cells))
        raise RuntimeError(
            f"Q207 parent coverage mismatch: missing={missing}, extra={extra}"
        )

    max_index = max(
        range(len(records)),
        key=lambda index: numeric[index][0],
    )
    min_index = min(
        range(len(records)),
        key=lambda index: numeric[index][1],
    )
    summary = {
        "q207_old_parent_cells": len(cells),
        "transport_panels": len(records),
        "certified_transport_panels": sum(
            row["strict_gate"] for row in records
        ),
        "coarse_transport_panels": sum(
            row["source"] == "coarse" for row in records
        ),
        "refined_transport_panels": sum(
            row["source"] == "refined" for row in records
        ),
        "transport_cover": [str(CORE_END), str(OLD_RADIUS)],
        "maximum_transport_ratio_upper": records[max_index][
            "transport_ratio_upper"
        ],
        "maximum_transport_ratio_x": [
            records[max_index]["x_low"],
            records[max_index]["x_high"],
        ],
        "minimum_old_cell_distance_lower": records[min_index][
            "q207_old_cell_distance_lower"
        ],
        "minimum_old_cell_distance_parent_x": [
            records[min_index]["q207_parent_x_low"],
            records[min_index]["q207_parent_x_high"],
        ],
        "complete_forward_collar": True,
        "finite_successor_instances": 1,
        "all_j_promotions": 0,
    }
    return records, summary


def validate_source_conclusions() -> None:
    q207 = load_json(Q207_RESULT)
    core = load_json(CORE_RESULT)
    right = load_json(RIGHT_RESULT)
    successor = load_json(SUCCESSOR_RESULT)
    coarse = load_json(COARSE_RESULT)
    refined = load_json(REFINED_RESULT)
    if q207.get("status") != "complete rigorous finite Q207 six-term interval theorem":
        raise RuntimeError("Q207 source theorem status drifted")
    if not q207.get("summary", {}).get("theorem_ready"):
        raise RuntimeError("Q207 source is not theorem-ready")
    if "|x|<=38" not in core.get("status", ""):
        raise RuntimeError("compact-core source status drifted")
    right_summary = right.get("summary", {}).get("right_strip", {})
    if not right_summary.get("full_derivative_positive"):
        raise RuntimeError("Q208 right-strip source is not derivative-positive")
    if (
        successor.get("status")
        != "exact conditional adiabatic phase-cell successor lemma with finite Q208 calibration; Xi all-j antecedent open"
    ):
        raise RuntimeError("successor lemma status drifted")
    if coarse.get("summary", {}).get("certified_panels") != 379:
        raise RuntimeError("coarse transport prefix source drifted")
    if not refined.get("summary", {}).get(
        "complete_collar_certificate"
    ):
        raise RuntimeError("refined transport source is incomplete")


def render_note(artifact: dict) -> str:
    summary = artifact["summary"]
    return f"""# Q207-to-Q208 Forward Adiabatic Successor Certificate

Date: {DATE}

Status: rigorous finite forward-successor certificate from the Q207 old bottom edge; not a proof of `Lambda<=0` or RH, and not an all-`j` theorem.

## Direction Audit

The earlier hybrid collar transported from cells on the new Q208
bottom edge. That proves the complete collar, but it is a
reverse-endpoint calibration. The conditional induction lemma instead
asks whether cells already known on the old Q207 edge absorb the
downward heat displacement.

For `x>=38`, the Q207 source stores origin-free cells for

```text
J=16*x^4*H, J_x.
```

On each old cell this certificate applies the exact
positive-determinant transform

```text
F=(1+x^(-4))*J,
F_x=(1+x^(-4))*J_x-4*x^(-5)*J,                  (1)
```

so that `F=16*(1+x^4)*H`, the regular proxy used by the collar
derivative bounds.

## Forward Transport

The old and new times are

```text
t_207=1/1035, t_208=1/1040,
delta_207=1/215280.
```

The {summary["q207_old_parent_cells"]} old Q207 parent cells cover
`38<=x<=245`. The stored rigorous heat-jet displacement bounds are
partitioned into {summary["coarse_transport_panels"]} half-unit panels
and {summary["refined_transport_panels"]} quarter-unit panels. Every
one of the {summary["transport_panels"]} transport panels satisfies

```text
delta_207 sup||(F_t,partial_t F_x)||_2
    < dist(0,C_207,k).                            (2)
```

The largest rigorous ratio is
`{summary["maximum_transport_ratio_upper"]}` on the panel
`[{summary["maximum_transport_ratio_x"][0]},{summary["maximum_transport_ratio_x"][1]}]`.

## Finite Composition

The independent compact transversality certificate already covers
`0<=x<=38` for the complete time interval. Equation (2) transports
the Q207 old edge across the new collar on `38<=x<=245`, and the
independent derivative-positive theorem covers the new right strip

```text
[1/1040,1/5]x[245,246].
```

Together with the contact-free Q207 base, these are exactly the old-edge
collar and right-strip hypotheses of the successor lemma. They give a
genuine forward finite Q207-to-Q208 successor and preserve the zero
first-jet degree. This independently recovers the finite Q208
contact-free theorem.

## Proof Boundary

This certificate proves one finite successor using the Q207 old-edge
cells. It does not provide a Q208-to-Q209 estimate, a uniform all-`j`
collar bound, a uniform right-strip cone, a cofinal boundary theorem,
`Lambda<=0`, RH, PF-infinity, or a Clay-prize proof.

Machine-audited files:

```text
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```
"""


def main() -> int:
    compact.flint.ctx.prec = PRECISION_BITS
    validate_source_conclusions()
    records, summary = build_records()
    artifact = {
        "kind": STEM,
        "date": DATE,
        "status": (
            "rigorous Q207-old-edge forward adiabatic successor "
            "certificate for Q208"
        ),
        "proof_boundary": (
            "This artifact proves one finite forward Q207-to-Q208 "
            "successor from Q207 old-edge cells. It proves no Q209 "
            "stage, no uniform all-j collar or right-strip theorem, "
            "no cofinal boundary theorem, no Lambda<=0, and no RH "
            "result."
        ),
        "contract": {
            "schema": "newman_q207_q208_forward_adiabatic_successor_v1",
            "old_time": str(OLD_TIME),
            "new_time": str(NEW_TIME),
            "time_step": str(TIME_STEP),
            "old_radius": str(OLD_RADIUS),
            "new_radius": str(NEW_RADIUS),
            "compact_core": [str(Fraction(0)), str(CORE_END)],
            "transported_collar": [str(CORE_END), str(OLD_RADIUS)],
            "new_right_strip": [str(OLD_RADIUS), str(NEW_RADIUS)],
            "old_proxy": "J=16*x^4*H",
            "transport_proxy": "F=16*(1+x^4)*H",
            "proxy_transform": (
                "F=(1+x^-4)J; "
                "F_x=(1+x^-4)J_x-4*x^-5*J"
            ),
            "strict_gate": (
                "delta_207*sup||(F_t,partial_t F_x)||"
                "<dist(0,C_207,k)"
            ),
            "source_sha256": source_hashes(),
        },
        "builder_sha256": file_hash(Path(__file__)),
        "summary": summary,
        "composition": {
            "base": "Q207 is contact-free with zero first-jet degree.",
            "compact_core": (
                "[0,1/5]x[0,38] is independently contact-free."
            ),
            "old_edge_transport": (
                "All 525 panels on [38,245] satisfy the old-edge "
                "adiabatic thickening gate."
            ),
            "right_strip": (
                "J_x>0 on [1/1040,1/5]x[245,246]."
            ),
            "conclusion": (
                "The exact successor lemma gives Q208 contact-free "
                "with the same zero first-jet degree."
            ),
        },
        "records": records,
    }
    DEFAULT_OUT.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    DEFAULT_NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built forward Q207-Q208 adiabatic successor: "
        f"{summary['q207_old_parent_cells']} old cells, "
        f"{summary['certified_transport_panels']}/"
        f"{summary['transport_panels']} transport panels, "
        "1 finite successor, 0 all-j promotions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
