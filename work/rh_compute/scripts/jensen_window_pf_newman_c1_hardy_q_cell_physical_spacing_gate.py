#!/usr/bin/env python3
"""Test whether adjacent physical pivots are joined by the saved q cells."""

from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
Q_CELLS = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_q_selector_cell_gate.json"
RAE_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_rae_coefficient_roster_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_q_cell_physical_spacing_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_q_cell_physical_spacing_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_q_cell_physical_spacing_gate.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def binary128_fraction(payload: str) -> Fraction:
    require(isinstance(payload, str) and len(payload) == 32, "invalid binary128 payload")
    bits = int(payload, 16)
    sign = -1 if bits >> 127 else 1
    exponent = (bits >> 112) & 0x7FFF
    fraction = bits & ((1 << 112) - 1)
    require(exponent != 0x7FFF, "non-finite binary128 payload")
    if exponent == 0:
        if fraction == 0:
            return Fraction(0)
        mantissa = fraction
        binary_exponent = 1 - 16383 - 112
    else:
        mantissa = (1 << 112) + fraction
        binary_exponent = exponent - 16383 - 112
    if binary_exponent >= 0:
        return Fraction(sign * mantissa * (2**binary_exponent), 1)
    return Fraction(sign * mantissa, 2 ** (-binary_exponent))


def exact_fraction(record: dict[str, int]) -> Fraction:
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def fraction_record(value: Fraction) -> dict[str, int]:
    return {"numerator": value.numerator, "denominator": value.denominator}


def decimal_text(value: Fraction, digits: int = 70) -> str:
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(value.numerator) / Decimal(value.denominator), ".55E")


def metric_record(value: Fraction) -> dict[str, Any]:
    return {"exact": fraction_record(value), "decimal": decimal_text(value)}


def circular_distance(left: Fraction, right: Fraction) -> tuple[Fraction, int]:
    delta = right - left
    candidates = [(abs(delta + shift), shift) for shift in (-1, 0, 1)]
    distance, shift = min(candidates, key=lambda item: (item[0], abs(item[1])))
    require(distance <= Fraction(1, 2), "circle-distance reduction failed")
    return distance, shift


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy q-cell physical-spacing obstruction gate

Date: 2026-08-06

Status: rigorous adjacent endpoint-cell obstruction; not a proof and no continuous coverage claim

## Test

The physical adapter fixes the first 32 block-20 pivots at

```text
rae_j = 657064 + 420*(j-1).
```

For each branch and each adjacent pair, this gate compares the exact binary128
`a1` centers with the exact rational `a1` radii of the two saved q-selector
cells.  Because `a1` is reduced modulo one, the comparison uses the shorter
distance on the circle.  This is generous: each ordinary real interval is
embedded into a circular arc, so disjoint circular arcs imply disjoint original
cell projections.

For every one of the `{aggregate['adjacent_pair_count']}` branch-adjacent pairs,

```text
circular_distance(a1_j,a1_(j+1)) > radius_j + radius_(j+1).
```

The minimum certified projected gap is
`{aggregate['minimum_phase_projected_gap']['decimal']}`.  The smallest ratio of
center distance to the sum of endpoint radii is
`{aggregate['minimum_phase_distance_to_radius_sum_ratio']['decimal']}`.

The raw `a1` formula is smooth for `rae/a>1`; modulo one it is a continuous map
to the circle.  A connected path joining two centers cannot stay in the union
of two disjoint endpoint arcs.  Therefore the two saved endpoint q cells cannot
by themselves certify the intervening continuous physical segment.  No
assumption about `fracL(rae)` is needed for this falsification because failure
in one coordinate already defeats four-dimensional box containment.

## What this changes

The local q cells are valid, but they are islands around discrete source calls,
not an adjacent interval cover.  Trying to stretch them across a physical step
of 420 is the wrong scaling strategy.

The source evaluator itself uses a discrete pivot roster, so continuous coverage
between consecutive pivots is not required merely to reproduce that finite
computation.  The economical next route is to prove the selector and analytic
margins directly on the discrete index `j`, first for all 212 pivots of block 20,
using modular/floor segmentation.  Adaptive intermediate q cells remain a valid
fallback only if a later theorem genuinely needs a continuous `rae` interval.

## Boundary

This gate does not invalidate the q-cell certificates.  It invalidates only the
proposal that adjacent endpoint cells already form a continuous physical cover.
It does not prove the remaining source pivots, special-function values, the
Riemann Hypothesis, or any prize-level conclusion.
"""


def main() -> int:
    for path in (Q_CELLS, RAE_GATE, CHECKER):
        require(path.is_file(), f"missing physical-spacing dependency: {path}")
    q_artifact = json.loads(Q_CELLS.read_text(encoding="utf-8"))
    rae_artifact = json.loads(RAE_GATE.read_text(encoding="utf-8"))
    require(q_artifact["aggregate"]["cell_count"] == 64, "q-cell count drift")
    require(rae_artifact["aggregate"]["pair_count"] == 32, "RAE roster count drift")
    require(rae_artifact["fixture"]["spacing"] == 420, "physical spacing drift")

    q_rows = {int(row["chain"]): row for row in q_artifact["rows"]}
    rae_rows = {int(row["sum_index"]): row for row in rae_artifact["rows"]}
    require(sorted(q_rows) == list(range(1, 65)), "q-cell chain sequence drift")

    centers: dict[tuple[int, int], dict[str, Any]] = {}
    for sum_index in range(1, 33):
        for branch in (1, 2):
            chain = 2 * (sum_index - 1) + branch
            q_row = q_rows[chain]
            rae_branch = rae_rows[sum_index]["branches"][branch - 1]
            require(int(q_row["sum_index"]) == sum_index, "q-cell sum-index drift")
            require(int(q_row["branch"]) == branch, "q-cell branch drift")
            saved = rae_branch["initial_coefficients_hex"]
            a1 = binary128_fraction(q_row["center"]["a1_hex"])
            y = binary128_fraction(q_row["center"]["xr_hex"])
            a3 = binary128_fraction(q_row["center"]["a3_hex"])
            require(a1 == binary128_fraction(saved[0]), "q-cell/physical a1 center mismatch")
            require(y == 2 * binary128_fraction(saved[1]), "q-cell/physical y center mismatch")
            require(a3 == binary128_fraction(saved[2]), "q-cell/physical a3 center mismatch")
            centers[(sum_index, branch)] = {
                "chain": chain,
                "rae": int(rae_rows[sum_index]["target_rae"]),
                "values": {
                    "a1_mod_one": a1,
                    "y_equals_2a2": y,
                    "a3": a3,
                },
                "radii": {
                    "a1_mod_one": exact_fraction(q_row["radii"]["a1"]["exact"]),
                    "y_equals_2a2": exact_fraction(q_row["radii"]["y"]["exact"]),
                    "a3": exact_fraction(q_row["radii"]["a3"]["exact"]),
                },
            }

    rows: list[dict[str, Any]] = []
    phase_gaps: list[tuple[Fraction, int, int]] = []
    phase_ratios: list[tuple[Fraction, int, int]] = []
    coordinate_disjoint_counts = {"a1_mod_one": 0, "y_equals_2a2": 0, "a3": 0}
    for branch in (1, 2):
        for sum_index in range(1, 32):
            left = centers[(sum_index, branch)]
            right = centers[(sum_index + 1, branch)]
            phase_distance, circle_shift = circular_distance(
                left["values"]["a1_mod_one"], right["values"]["a1_mod_one"]
            )
            phase_radius_sum = left["radii"]["a1_mod_one"] + right["radii"]["a1_mod_one"]
            phase_gap = phase_distance - phase_radius_sum
            require(phase_gap > 0, f"adjacent phase arcs overlap at branch {branch}, index {sum_index}")
            phase_ratio = phase_distance / phase_radius_sum
            phase_gaps.append((phase_gap, sum_index, branch))
            phase_ratios.append((phase_ratio, sum_index, branch))

            coordinate_rows: dict[str, Any] = {}
            for coordinate in ("a1_mod_one", "y_equals_2a2", "a3"):
                if coordinate == "a1_mod_one":
                    distance = phase_distance
                else:
                    distance = abs(left["values"][coordinate] - right["values"][coordinate])
                radius_sum = left["radii"][coordinate] + right["radii"][coordinate]
                signed_gap = distance - radius_sum
                disjoint = signed_gap > 0
                coordinate_disjoint_counts[coordinate] += int(disjoint)
                coordinate_rows[coordinate] = {
                    "center_distance": metric_record(distance),
                    "endpoint_radius_sum": metric_record(radius_sum),
                    "signed_projected_gap": metric_record(signed_gap),
                    "endpoint_projections_disjoint": disjoint,
                    "distance_to_radius_sum_ratio": metric_record(distance / radius_sum),
                }
            rows.append(
                {
                    "branch": branch,
                    "left_sum_index": sum_index,
                    "right_sum_index": sum_index + 1,
                    "left_chain": left["chain"],
                    "right_chain": right["chain"],
                    "left_rae": left["rae"],
                    "right_rae": right["rae"],
                    "rae_step": right["rae"] - left["rae"],
                    "short_circle_integer_shift": circle_shift,
                    "coordinates": coordinate_rows,
                    "endpoint_q_cells_form_continuous_cover": False,
                }
            )

    minimum_gap = min(phase_gaps, key=lambda item: item[0])
    minimum_ratio = min(phase_ratios, key=lambda item: item[0])
    require(coordinate_disjoint_counts["a1_mod_one"] == 62, "phase obstruction is not universal")
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_q_cell_physical_spacing_gate",
        "status": "rigorous_adjacent_endpoint_q_cell_obstruction_with_discrete_route_open",
        "test_contract": {
            "physical_roster": "rae_j=657064+420*(j-1), j=1,...,32",
            "branches": [1, 2],
            "tested_adjacent_pairs_per_branch": 31,
            "phase_metric": "min_{n in Z}|a1_right-a1_left+n|",
            "topological_reason": "Disjoint circular a1 projections make the union of two endpoint boxes disconnected between their centers.",
            "generous_enlargement": "Real a1 intervals are embedded as circle arcs before testing disjointness.",
        },
        "aggregate": {
            "adjacent_pair_count": 62,
            "phase_projected_disjoint_count": coordinate_disjoint_counts["a1_mod_one"],
            "y_projected_disjoint_count": coordinate_disjoint_counts["y_equals_2a2"],
            "a3_projected_disjoint_count": coordinate_disjoint_counts["a3"],
            "endpoint_continuous_cover_count": 0,
            "minimum_phase_projected_gap": {
                "left_sum_index": minimum_gap[1],
                "branch": minimum_gap[2],
                **metric_record(minimum_gap[0]),
            },
            "minimum_phase_distance_to_radius_sum_ratio": {
                "left_sum_index": minimum_ratio[1],
                "branch": minimum_ratio[2],
                **metric_record(minimum_ratio[0]),
            },
        },
        "rows": rows,
        "route_decision": {
            "rejected": "Treating the two saved endpoint q cells as a continuous cover of an adjacent physical step.",
            "recommended": "Exploit the evaluator's discrete pivot index: derive selector/margin segments over j and cover all 212 block-20 calls without filling irrelevant continuous gaps.",
            "fallback": "If a later theorem needs continuous rae coverage, generate adaptive intermediate q cells and certify the full four-coordinate map including fracL.",
            "falsification_rule": "A future claimed adjacent cover must close the exact positive phase gaps recorded here; endpoint validity alone is insufficient.",
        },
        "sources": {
            "q_cells": {"path": relative(Q_CELLS), "sha256": file_hash(Q_CELLS)},
            "rae_gate": {"path": relative(RAE_GATE), "sha256": file_hash(RAE_GATE)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "proof_boundary": (
            "This rigorously rejects only adjacent continuous coverage by the saved endpoint q cells. "
            "It does not reject the local cells, does not yet cover all discrete source calls, and does "
            "not prove special-function accuracy, RH, or a prize-level theorem."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print(
        "built Hardy q-cell physical-spacing gate: "
        f"{artifact['aggregate']['phase_projected_disjoint_count']}/62 phase projections disjoint, "
        "discrete-index route selected"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
