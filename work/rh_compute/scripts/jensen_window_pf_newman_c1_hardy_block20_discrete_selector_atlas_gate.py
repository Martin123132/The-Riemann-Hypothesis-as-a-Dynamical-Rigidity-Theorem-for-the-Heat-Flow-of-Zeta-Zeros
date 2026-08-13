#!/usr/bin/env python3
"""Build the complete discrete selector atlas for Hardy block 20."""

from __future__ import annotations

from collections import Counter
from decimal import Decimal
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as local_cells


FIXTURE_RESULT = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10_block20_full/fixture_result.json"
)
TELEMETRY = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10_block20_full/enabled/chain_telemetry.jsonl"
)
PREFIX_Q_CELLS = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_q_selector_cell_gate.json"
RAE_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_rae_coefficient_roster_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_discrete_selector_atlas_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_discrete_selector_atlas_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_discrete_selector_atlas_gate.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def load_chains() -> dict[int, dict[str, Any]]:
    chains: dict[int, dict[str, Any]] = {}
    for line in TELEMETRY.read_text(encoding="utf-8").splitlines():
        record = json.loads(line)
        chain_id = int(record.get("chain", 0))
        if record["type"] == "chain":
            chains[chain_id] = {"header": record, "levels": {}, "q_terms": {}, "recurrences": {}}
        elif record["type"] == "level":
            chains[chain_id]["levels"][int(record["level"])] = record
        elif record["type"] == "q_terms":
            chains[chain_id]["q_terms"][int(record["nit"])] = record
        elif record["type"] == "recurrence":
            chains[chain_id]["recurrences"][int(record["nit"])] = record
    require(sorted(chains) == list(range(1, 425)), "full selector chain sequence drift")
    return chains


def selector_word(row: dict[str, Any]) -> tuple[Any, ...]:
    if row["route"] == "direct_kernel_mit1":
        return ("MIT1", 104)
    cell = row["q_cell"]
    child = cell["child_selectors"]
    q_branch = cell["q_branches"]
    return (
        "MIT2",
        row["child_length"],
        child["n1"],
        child["ix"],
        child["linear_residual_sign"],
        child["child_xr_sign"],
        child["source_conjugates"],
        q_branch["a1_sign"],
        q_branch["jbot"],
        q_branch["ip1"],
        q_branch["newton_used"],
    )


def maximal_segments(rows: list[dict[str, Any]], branch: int) -> list[dict[str, Any]]:
    selected = [row for row in rows if row["branch"] == branch]
    require([row["sum_index"] for row in selected] == list(range(1, 213)), "branch index gap")
    segments: list[dict[str, Any]] = []
    start = 0
    current = selector_word(selected[0])
    for index in range(1, len(selected)):
        candidate = selector_word(selected[index])
        if candidate != current:
            segments.append(
                {
                    "start_sum_index": selected[start]["sum_index"],
                    "end_sum_index": selected[index - 1]["sum_index"],
                    "start_rae": selected[start]["rae"],
                    "end_rae": selected[index - 1]["rae"],
                    "selector_word": list(current),
                }
            )
            start = index
            current = candidate
    segments.append(
        {
            "start_sum_index": selected[start]["sum_index"],
            "end_sum_index": selected[-1]["sum_index"],
            "start_rae": selected[start]["rae"],
            "end_rae": selected[-1]["rae"],
            "selector_word": list(current),
        }
    )
    return segments


def minimum_decimal(rows: list[dict[str, Any]], path: tuple[str, ...]) -> str:
    values: list[Decimal] = []
    for row in rows:
        value: Any = row
        for key in path:
            value = value[key]
        values.append(Decimal(value))
    return str(min(values))


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 discrete selector atlas

Date: 2026-08-06

Status: complete finite selector atlas for all 212 block-20 pivots; not a proof and analytic q remains open

## Source roster

The separate full-block telemetry fixture records both cubic branches at every
source pivot

```text
rae_j = 657064 + 420*(j-1),  1 <= j <= 212.
```

All 424 chains end in the same exact evaluator checkpoint state and the same
fifteen displayed Hardy values as the accepted 64-chain fixture.  The observer
therefore exposes source decisions without changing the finite calculation.

## The transition the prefix missed

The first 64 chains all had `MIT=2`.  Across the full block, however,

```text
MIT=1 direct-kernel chains : {aggregate['mit1_chain_count']}
MIT=2 one-step chains      : {aggregate['mit2_chain_count']}.
```

The 50 direct cases are scattered through the later roster rather than forming
one terminal tail.  They have one saved level, kernel length 104, and no q
recurrence.  They are recorded as their own exact selector word and are not
forced into the earlier `MIT=2` model.

For each of the 374 recursive chains, the Section 11.243 Arb constructor was
replayed at the exact saved center.  Every chain admits a nonzero-radius box
preserving the child integer/parity/orientation decisions, q sign and loop
branches, all six `PSI` paths, both `ERF` paths, frac endpoints, denominator,
and stationary discriminant.  The maximum shrink was
`{aggregate['maximum_q_cell_halvings']}` halvings.

The complete discrete route therefore contains

```text
424 / 424 classified branch calls,
374 rigorous local q-selector cells,
 50 exact direct-kernel calls,
  0 unclassified calls.
```

Branch 1 has `{aggregate['branch_1_selector_segment_count']}` maximal constant
selector-word segments and branch 2 has
`{aggregate['branch_2_selector_segment_count']}`.  This high fragmentation is
why extrapolating the first 32 pivots was unsafe.

## What is and is not covered

This is the all-212 **discrete** selector atlas required by the actual evaluator.
It does not claim continuous coverage between pivots; Section 11.245 proves that
the saved endpoint cells cannot provide such a cover.  It also does not yet
replace the source `PSI`/`ERF` approximations with rigorous values or enclose the
analytic `W1`--`W5` combination and `q`.  Direct `MIT=1` classification is not an
external-error theorem.

The next target is route-specific: rigorously evaluate special functions and the
correlated q correction on the 374 occupied cells, while separately auditing the
50 direct kernels.  Recurrence accumulation and the outer Hardy representation
remain independent obligations.  No physical carrier, determinant sign,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion follows here.
"""


def main() -> int:
    for path in (FIXTURE_RESULT, TELEMETRY, PREFIX_Q_CELLS, RAE_GATE, CHECKER):
        require(path.is_file(), f"missing selector-atlas dependency: {path}")
    fixture = json.loads(FIXTURE_RESULT.read_text(encoding="utf-8"))
    require(fixture["branch_count"] == 424, "full fixture branch count drift")
    require(fixture["equivalence"]["checkpoint_state_equal_ignoring_run_id"] is True, "fixture state drift")
    require(fixture["equivalence"]["displayed_hardy_values_exact"] is True, "fixture value drift")
    prefix = json.loads(PREFIX_Q_CELLS.read_text(encoding="utf-8"))
    prefix_rows = {int(row["chain"]): row for row in prefix["rows"]}
    chains = load_chains()

    rows: list[dict[str, Any]] = []
    q_cells: list[dict[str, Any]] = []
    route_counter: Counter[str] = Counter()
    child_lengths: Counter[int] = Counter()
    mit1_indices: dict[str, list[int]] = {"1": [], "2": []}
    for chain_id in range(1, 425):
        data = chains[chain_id]
        header = data["header"]
        sum_index = int(header["sum_index"])
        branch = int(header["branch"])
        require(chain_id == 2 * (sum_index - 1) + branch, "physical chain ordering drift")
        mit = int(header["mit"])
        require(mit in (1, 2), f"chain {chain_id} unexpected MIT")
        require(sorted(data["levels"]) == list(range(1, mit + 1)), f"chain {chain_id} level gap")
        require(len(data["recurrences"]) == mit - 1, f"chain {chain_id} recurrence count drift")
        require(len(data["q_terms"]) == mit - 1, f"chain {chain_id} q count drift")
        rae = 657064 + 420 * (sum_index - 1)
        if mit == 1:
            route = "direct_kernel_mit1"
            require(int(header["kernel_length"]) == 104, f"chain {chain_id} direct length drift")
            mit1_indices[str(branch)].append(sum_index)
            row = {
                "chain": chain_id,
                "sum_index": sum_index,
                "rae": rae,
                "branch": branch,
                "route": route,
                "mit": 1,
                "initial_length": 104,
                "child_length": None,
                "initial_coefficients_hex": header["initial_coefficients_hex"],
                "level_one_fracL_hex": data["levels"][1]["frac_length_hex"],
                "q_cell": None,
            }
        else:
            route = "recursive_mit2_q_cell"
            cell = local_cells.selector_cell(
                chain_id, header, data["levels"][1], data["levels"][2]
            )
            child_length = int(header["kernel_length"])
            require(child_length == int(data["levels"][2]["length"]), "child length mismatch")
            child_lengths[child_length] += 1
            q_cells.append(cell)
            if chain_id <= 64:
                require(cell == prefix_rows[chain_id], f"chain {chain_id} prefix q-cell replay drift")
            row = {
                "chain": chain_id,
                "sum_index": sum_index,
                "rae": rae,
                "branch": branch,
                "route": route,
                "mit": 2,
                "initial_length": 104,
                "child_length": child_length,
                "initial_coefficients_hex": header["initial_coefficients_hex"],
                "level_one_fracL_hex": data["levels"][1]["frac_length_hex"],
                "q_cell": cell,
            }
        route_counter[route] += 1
        rows.append(row)

    require(route_counter == {"recursive_mit2_q_cell": 374, "direct_kernel_mit1": 50}, "route count drift")
    halving_histogram = Counter(int(cell["halvings"]) for cell in q_cells)
    selector_word_histogram = Counter(json.dumps(selector_word(row), separators=(",", ":")) for row in rows)
    branch_segments = {str(branch): maximal_segments(rows, branch) for branch in (1, 2)}
    recursive_rows = [row["q_cell"] for row in rows if row["q_cell"] is not None]
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_discrete_selector_atlas_gate",
        "status": "complete_424_call_discrete_selector_atlas_with_analytic_q_open",
        "roster": {
            "block": 20,
            "formula": "rae_j=657064+420*(j-1)",
            "sum_count": 212,
            "branch_count": 424,
            "first_rae": 657064,
            "last_rae": 745684,
        },
        "selector_word_schema": [
            "route",
            "child_length",
            "child_n1",
            "child_ix",
            "child_linear_sign",
            "child_xr_sign",
            "source_conjugates",
            "q_a1_sign",
            "q_jbot",
            "q_ip1",
            "q_newton_used",
        ],
        "aggregate": {
            "classified_call_count": 424,
            "unclassified_call_count": 0,
            "mit1_chain_count": route_counter["direct_kernel_mit1"],
            "mit2_chain_count": route_counter["recursive_mit2_q_cell"],
            "child_length_histogram": {str(key): value for key, value in sorted(child_lengths.items())},
            "q_cell_count": len(q_cells),
            "maximum_q_cell_halvings": max(halving_histogram),
            "q_cell_halving_histogram": {str(key): value for key, value in sorted(halving_histogram.items())},
            "distinct_selector_word_count": len(selector_word_histogram),
            "branch_1_selector_segment_count": len(branch_segments["1"]),
            "branch_2_selector_segment_count": len(branch_segments["2"]),
            "minimum_q_cell_radii": {
                name: min((cell["radii"][name]["decimal"] for cell in q_cells), key=Decimal)
                for name in ("a1", "y", "a3", "fracL")
            },
            "minimum_psi_boundary_gap": minimum_decimal(recursive_rows, ("q_branches", "minimum_psi_boundary_gap")),
            "minimum_erf_boundary_gap": minimum_decimal(recursive_rows, ("q_branches", "minimum_erf_boundary_gap")),
            "minimum_denominator_abs_lower": minimum_decimal(recursive_rows, ("analytic_margins", "minimum_denominator_abs_lower")),
            "minimum_stationary_discriminant_lower": minimum_decimal(recursive_rows, ("analytic_margins", "minimum_stationary_discriminant_lower")),
        },
        "mit1_sum_indices_by_branch": mit1_indices,
        "selector_word_histogram": dict(sorted(selector_word_histogram.items())),
        "selector_segments_by_branch": branch_segments,
        "rows": rows,
        "next_handoff": {
            "recursive_route": "Replace source PSI/ERF approximations by rigorous values on the 374 occupied q cells and enclose correlated t1 through t5 and q.",
            "direct_route": "Audit the 50 MIT=1 kernel evaluations separately; no recurrence q correction is present on those calls.",
            "falsification_rule": "Any analytic enclosure must preserve the exact route and selector word recorded for its discrete call; do not interpolate across MIT transitions.",
        },
        "sources": {
            "fixture_result": {"path": relative(FIXTURE_RESULT), "sha256": file_hash(FIXTURE_RESULT)},
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "prefix_q_cells": {"path": relative(PREFIX_Q_CELLS), "sha256": file_hash(PREFIX_Q_CELLS)},
            "rae_gate": {"path": relative(RAE_GATE), "sha256": file_hash(RAE_GATE)},
            "local_cell_builder": {
                "path": relative(Path(local_cells.__file__).resolve()),
                "sha256": file_hash(Path(local_cells.__file__).resolve()),
            },
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "proof_boundary": (
            "Complete finite discrete selector atlas for one low-height 212-pivot block only. "
            "It proves neither continuous coverage nor source special-function accuracy, analytic q, "
            "recurrence/outer representation bounds, physical-height scaling, RH, or a prize-level theorem."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print(
        "built Hardy block-20 discrete selector atlas: "
        "424 classified calls, 374 q cells, 50 direct kernels, 0 unclassified"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
