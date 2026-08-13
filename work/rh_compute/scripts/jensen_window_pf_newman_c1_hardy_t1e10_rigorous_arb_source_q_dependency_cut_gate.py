#!/usr/bin/env python3
"""Cut emitted source-q payloads out of the finite corrected-model path."""

from __future__ import annotations

import ast
from fractions import Fraction
import hashlib
import inspect
import json
import os
from pathlib import Path
import sys
import textwrap
import time
from typing import Any, Callable


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import jensen_window_pf_newman_c1_hardy_block20_coefficient_neighborhood_transport_gate as block20_cells
import jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_full_gate as log_full_builder
import jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_pilot_gate as log_pilot
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
import jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_full_gate as full_cells
import jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_pilot_gate as evaluator
import jensen_window_pf_newman_c1_hardy_t1e10_crossblock_direct_kernel_output_gate as direct_builder
import jensen_window_pf_newman_c1_hardy_t1e10_crossblock_recursive_output_budget_gate as transport
import jensen_window_pf_newman_c1_hardy_t1e10_later_complete_recurrence_cell_full_gate as tail_builder
from flint import acb


FULL_CELLS = full_cells.RESULT
FULL_TAILS = tail_builder.RESULT
LOG_FULL = log_full_builder.RESULT
DIRECT_OUTPUT = direct_builder.RESULT
TELEMETRY = full_cells.TELEMETRY
WEIGHTS = transport.WEIGHTS
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_rigorous_arb_source_q_dependency_cut_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_rigorous_arb_source_q_dependency_cut_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_t1e10_rigorous_arb_source_q_dependency_cut_gate.py"

EXPECTED_CALLS = 1040
STRESS_CHAINS = (626, 653, 3406, 3423)
REPLAY_PRECISION = 110
REQUESTED_SCALE = Fraction(5, 1000)
SERIALIZATION_PAD = Fraction(1, 10**50)
FORBIDDEN_PAYLOAD_KEYS = {
    "q_terms",
    "steps",
    "qq_hex",
    "t1_hex",
    "t2_hex",
    "t3_hex",
    "t4_hex",
    "t5_hex",
    "state_before_hex",
    "multiplier_hex",
    "psi_values",
    "erf_values",
    "erfc_values",
}
FORBIDDEN_SOURCE_TOKENS = tuple(sorted(FORBIDDEN_PAYLOAD_KEYS))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def source_audit(name: str, function: Callable[..., Any]) -> dict[str, Any]:
    source = textwrap.dedent(inspect.getsource(function))
    tree = ast.parse(source)
    data_keys: set[str] = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Subscript):
            continue
        if not isinstance(node.value, ast.Name) or node.value.id != "data":
            continue
        if isinstance(node.slice, ast.Constant) and isinstance(node.slice.value, str):
            data_keys.add(node.slice.value)
    forbidden = [token for token in FORBIDDEN_SOURCE_TOKENS if token in source]
    require(not forbidden, f"{name} imports forbidden payload tokens: {forbidden}")
    require(data_keys <= {"levels"}, f"{name} reads non-level data keys: {sorted(data_keys)}")
    return {
        "function": name,
        "sha256": hashlib.sha256(source.encode("utf-8")).hexdigest(),
        "data_keys": sorted(data_keys),
        "forbidden_payload_token_hits": forbidden,
    }


def evaluator_manifest() -> list[dict[str, Any]]:
    functions = (
        ("block20_cells.coefficient_box", block20_cells.coefficient_box),
        ("evaluator.coefficient_box", evaluator.coefficient_box),
        ("evaluator.coefficient_balls", evaluator.coefficient_balls),
        ("evaluator.cell_sector_contract", evaluator.cell_sector_contract),
        ("evaluator.ray_integral_cell", evaluator.ray_integral_cell),
        ("evaluator.finite_zero_mode_cell", evaluator.finite_zero_mode_cell),
        ("evaluator.exact_nonsaddle_cell", evaluator.exact_nonsaddle_cell),
        ("evaluator.paper_components_cell", evaluator.paper_components_cell),
        ("evaluator.evaluate_cell", evaluator.evaluate_cell),
        ("log_pilot.direction", log_pilot.direction),
        ("log_pilot.kernel", log_pilot.kernel),
        ("log_pilot.error_box", log_pilot.error_box),
        ("log_pilot.harmonic", log_pilot.harmonic),
        ("log_full_builder.breaks_for_cutoff", log_full_builder.breaks_for_cutoff),
    )
    return [source_audit(name, function) for name, function in functions]


def load_stripped_levels(selected: set[int]) -> dict[int, dict[str, Any]]:
    levels = {chain: {} for chain in selected}
    with TELEMETRY.open(encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            chain = int(record.get("chain", 0))
            if chain in selected and record.get("type") == "level":
                levels[chain][int(record["level"])] = record
    for chain, rows in levels.items():
        require(set(rows) >= {1, 2}, f"stripped replay chain {chain} level roster missing")
    return {chain: {"levels": rows} for chain, rows in levels.items()}


def stored_correction(row: dict[str, Any]) -> acb:
    record = row["high_precision_correction"]
    return acb(
        evaluator.interval_from_record(record["real"]),
        evaluator.interval_from_record(record["imag"]),
    )


def replay_stress_cells(cell_rows: dict[int, dict[str, Any]]) -> list[dict[str, Any]]:
    stripped = load_stripped_levels(set(STRESS_CHAINS))
    rows: list[dict[str, Any]] = []
    for chain in STRESS_CHAINS:
        source = cell_rows[chain]
        data = stripped[chain]
        require(set(data) == {"levels"}, f"chain {chain} stripped input gained a top-level key")
        box = evaluator.coefficient_box(data, source["refined_cell"])
        replay = evaluator.evaluate_cell(box, data, REPLAY_PRECISION)
        saved = stored_correction(source)
        overlap = replay["correction"].overlaps(saved)
        require(overlap, f"chain {chain} stripped evaluator replay missed saved certificate")
        rows.append(
            {
                "chain": chain,
                "block": int(source["block"]),
                "supplied_top_level_keys": sorted(data),
                "precision_decimal_digits": REPLAY_PRECISION,
                "selector": replay["selector"],
                "saved_correction_overlap": overlap,
                "replayed_correction_magnitude_upper": decimal(evaluator.upper_abs(replay["correction"])),
            }
        )
    return rows


def forbidden_projection_keys(value: Any) -> set[str]:
    hits: set[str] = set()
    if isinstance(value, dict):
        for key, child in value.items():
            if str(key) in FORBIDDEN_PAYLOAD_KEYS:
                hits.add(str(key))
            hits.update(forbidden_projection_keys(child))
    elif isinstance(value, list):
        for child in value:
            hits.update(forbidden_projection_keys(child))
    return hits


def projection_digest(
    cell_rows: list[dict[str, Any]],
    tail_rows: dict[int, dict[str, Any]],
) -> tuple[str, int]:
    projection: list[dict[str, Any]] = []
    for row in cell_rows:
        chain = int(row["chain"])
        tail = tail_rows[chain]
        projected = {
            "chain": chain,
            "block": int(row["block"]),
            "sum_index": int(row["sum_index"]),
            "branch": int(row["branch"]),
            "refined_cell": row["refined_cell"],
            "high_precision_correction": row["high_precision_correction"],
            "cell_correction_magnitude_upper": row["cell_correction_magnitude_upper"],
            "cell_complete_tail_majorant_upper": tail["cell_complete_tail_majorant_upper"],
            "minimum_stationary_discriminant_lower": tail["minimum_stationary_discriminant_lower"],
        }
        hits = forbidden_projection_keys(projected)
        require(not hits, f"chain {chain} analytic projection contains forbidden keys: {sorted(hits)}")
        projection.append(projected)
    material = json.dumps(projection, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(material).hexdigest(), len(material)


def output_certificate(
    cell_rows: list[dict[str, Any]],
    tail_rows: dict[int, dict[str, Any]],
    full_cell_outputs: dict[int, dict[str, Any]],
    full_tail_outputs: dict[int, dict[str, Any]],
) -> list[dict[str, Any]]:
    factors = transport.load_weight_factors()
    log_full = json.loads(LOG_FULL.read_text(encoding="utf-8"))
    direct = json.loads(DIRECT_OUTPUT.read_text(encoding="utf-8"))
    require(
        log_full["status"] == "rigorous_all_374_logarithmic_endpoint_ray_nonsaddle_enclosures_validated",
        "block-20 analytic dependency failed",
    )
    require(
        direct["status"] == "rigorous_t1e10_all_direct_kernels_and_complete_exact_point_output_budget",
        "direct-kernel dependency failed",
    )
    require(len(log_full["rows"]) == 374, "block-20 analytic roster drift")
    require(len(direct["output_rows"]) == 15, "direct-kernel output roster drift")
    direct_outputs = {int(row["output_index"]): row for row in direct["output_rows"]}

    result: list[dict[str, Any]] = []
    for output in range(1, 16):
        block20 = Fraction(0)
        for row in log_full["rows"]:
            factor = factors[(20, int(row["sum_index"]), output, int(row["branch"]))]
            block20 += factor * (Fraction(row["exact_minus_paper_magnitude_upper"]) + SERIALIZATION_PAD)

        cells_column = Fraction(0)
        tails_column = Fraction(0)
        for row in cell_rows:
            chain = int(row["chain"])
            factor = factors[(int(row["block"]), int(row["sum_index"]), output, int(row["branch"]))]
            require(factor >= 0, f"output {output} has a negative outer factor")
            cells_column += factor * (Fraction(row["cell_correction_magnitude_upper"]) + SERIALIZATION_PAD)
            tails_column += factor * (Fraction(tail_rows[chain]["cell_complete_tail_majorant_upper"]) + SERIALIZATION_PAD)

        direct_column = Fraction(direct_outputs[output]["direct_kernel_majorant_upper"]) + SERIALIZATION_PAD
        source_q_free_cell = block20 + direct_column + cells_column
        source_q_free_complete = source_q_free_cell + tails_column
        stored_cell = Fraction(full_cell_outputs[output]["all_recursive_cell_complete_majorant_upper"])
        stored_partial = Fraction(full_tail_outputs[output]["later_complete_recurrence_partial_majorant_upper"])
        require(source_q_free_cell >= stored_cell, f"output {output} source-q-free cell recombination lost the accepted bound")
        require(source_q_free_complete >= stored_partial, f"output {output} source-q-free complete recombination lost the accepted bound")
        require(source_q_free_complete < REQUESTED_SCALE, f"output {output} source-q-free recombination exceeds requested scale")
        result.append(
            {
                "output_index": output,
                "output_label": full_cell_outputs[output]["output_label"],
                "block20_analytic_correction_upper": decimal(block20),
                "direct_kernel_non_q_upper": decimal(direct_column),
                "later_analytic_cell_correction_upper": decimal(cells_column),
                "later_analytic_complete_tail_upper": decimal(tails_column),
                "source_q_free_cell_majorant_upper": decimal(source_q_free_cell),
                "source_q_free_complete_majorant_upper": decimal(source_q_free_complete),
                "dominance_over_stored_cell_upper": decimal(source_q_free_cell - stored_cell),
                "dominance_over_stored_partial_upper": decimal(source_q_free_complete - stored_partial),
                "margin_below_requested_scale": decimal(REQUESTED_SCALE - source_q_free_complete),
                "within_requested_scale": True,
            }
        )
    return result


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Rigorous Arb source-q dependency cut

Date: 2026-08-09

Status: finite saved-height source-q dependency cut closed; not a proof of the outer Hardy theorem, height uniformity, or RH

The corrected-model cell evaluator was audited function by function.  Its
only access to the chain record is through the exact coefficient `levels`.
It does not read emitted `q_terms`, recurrence `steps`, `qq_hex`, `t1`--`t5`,
saved PSI/ERF payloads, source recurrence states, or source multipliers.
Mathematical pi enters through `arb.pi()`, while `erf` and `digamma` are
evaluated independently by Arb/FLINT rather than imported from the source.

As an executable falsification check, chains {', '.join(map(str, STRESS_CHAINS))}
were replayed at {REPLAY_PRECISION} decimal digits with a dictionary whose only
top-level key was `levels`.  All {len(artifact['stripped_replay_rows'])} replayed
correction balls overlap the saved independent cell certificates.

For all {EXPECTED_CALLS} later recursive cells, a source-q-free projection of
the accepted cell, Arb correction, and analytic complete-tail fields has SHA-256
`{artifact['projection']['sha256']}`.  Exact nonnegative outer transport was
then rebuilt directly as

```text
block-20 analytic correction
  + direct nonrecursive kernel column
  + later Arb cell correction
  + later analytic complete tail.
```

No source-q point column or signed cross-call cancellation occurs in this
recombination.  The largest complete output is
`{aggregate['maximum_source_q_free_complete_majorant_upper']}` and the minimum
margin below `0.005` is `{aggregate['minimum_margin_below_requested_scale']}`.
All 15 outputs close.

This proves a dependency classification, not RH.  Emitted source `q`, PSI/ERF
values, and recurrence states are not mathematical antecedents of this finite
corrected-model Arb bound.  They remain relevant only to the separate task of
reproducing or proving the legacy executable itself.  The direct-kernel column
still audits nonrecursive source arithmetic, and the outer Hardy representation
and remainder, coefficient coverage beyond the saved height, a height-uniform
route theorem, `Lambda<=0`, PF-infinity, RH, and a prize-level proof remain open.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    manifest = evaluator_manifest()

    full = json.loads(FULL_CELLS.read_text(encoding="utf-8"))
    tails = json.loads(FULL_TAILS.read_text(encoding="utf-8"))
    require(bool(full["passed"]) and bool(tails["passed"]), "accepted cell/tail dependency failed")
    require(len(full["rows"]) == EXPECTED_CALLS and len(tails["rows"]) == EXPECTED_CALLS, "later roster drift")
    cell_rows = sorted(full["rows"], key=lambda row: int(row["chain"]))
    cell_by_chain = {int(row["chain"]): row for row in cell_rows}
    tail_by_chain = {int(row["chain"]): row for row in tails["rows"]}
    require(set(cell_by_chain) == set(tail_by_chain), "cell/tail chain roster mismatch")

    replay = replay_stress_cells(cell_by_chain)
    projection_sha, projection_bytes = projection_digest(cell_rows, tail_by_chain)
    full_cell_outputs = {int(row["output_index"]): row for row in full["output_rows"]}
    full_tail_outputs = {int(row["output_index"]): row for row in tails["output_rows"]}
    outputs = output_certificate(cell_rows, tail_by_chain, full_cell_outputs, full_tail_outputs)

    maximum = max((Fraction(row["source_q_free_complete_majorant_upper"]), int(row["output_index"])) for row in outputs)
    minimum_margin = min((Fraction(row["margin_below_requested_scale"]), int(row["output_index"])) for row in outputs)
    aggregate = {
        "evaluator_function_count": len(manifest),
        "forbidden_payload_token_hit_count": sum(len(row["forbidden_payload_token_hits"]) for row in manifest),
        "stripped_replay_count": len(replay),
        "stripped_replay_overlap_count": sum(bool(row["saved_correction_overlap"]) for row in replay),
        "projected_later_cell_count": len(cell_rows),
        "output_count": len(outputs),
        "outputs_within_requested_scale": sum(bool(row["within_requested_scale"]) for row in outputs),
        "uses_source_q_payload": False,
        "uses_signed_cross_call_cancellation": False,
        "maximum_source_q_free_complete_majorant_upper": decimal(maximum[0]),
        "maximum_source_q_free_complete_witness_output": maximum[1],
        "minimum_margin_below_requested_scale": decimal(minimum_margin[0]),
        "minimum_margin_witness_output": minimum_margin[1],
    }
    passed = (
        aggregate["forbidden_payload_token_hit_count"] == 0
        and aggregate["stripped_replay_overlap_count"] == len(STRESS_CHAINS)
        and aggregate["projected_later_cell_count"] == EXPECTED_CALLS
        and aggregate["outputs_within_requested_scale"] == 15
    )
    require(passed, "source-q dependency cut did not close")

    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_rigorous_arb_source_q_dependency_cut_gate",
        "status": "finite_saved_height_corrected_model_arb_path_has_no_emitted_source_q_payload_dependency",
        "passed": True,
        "scope": {
            "height": "t=10^10 saved atlas",
            "later_recursive_cells": EXPECTED_CALLS,
            "outputs": 15,
            "stripped_replay_chains": list(STRESS_CHAINS),
        },
        "decision": {
            "source_q_payload_is_corrected_model_antecedent": False,
            "source_q_payload_remains_legacy_reproduction_obligation": True,
            "independent_special_functions": ["Arb erf", "Arb digamma", "Arb pi"],
        },
        "evaluator_manifest": manifest,
        "stripped_replay_rows": replay,
        "projection": {
            "row_count": len(cell_rows),
            "canonical_json_bytes": projection_bytes,
            "sha256": projection_sha,
            "forbidden_payload_keys": sorted(FORBIDDEN_PAYLOAD_KEYS),
            "forbidden_key_hit_count": 0,
        },
        "output_rows": outputs,
        "aggregate": aggregate,
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "process_priority": priority,
            "workers": 1,
            "flint_threads": 1,
        },
        "dependencies": {
            "accepted_later_cells": {"path": relative(FULL_CELLS), "sha256": file_hash(FULL_CELLS)},
            "accepted_later_tails": {"path": relative(FULL_TAILS), "sha256": file_hash(FULL_TAILS)},
            "block20_analytic_corrections": {"path": relative(LOG_FULL), "sha256": file_hash(LOG_FULL)},
            "direct_nonrecursive_column": {"path": relative(DIRECT_OUTPUT), "sha256": file_hash(DIRECT_OUTPUT)},
            "outer_weights": {"path": relative(WEIGHTS), "sha256": file_hash(WEIGHTS)},
            "coefficient_telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "cell_evaluator": {"path": relative(Path(evaluator.__file__).resolve()), "sha256": file_hash(Path(evaluator.__file__).resolve())},
            "coefficient_box": {"path": relative(Path(block20_cells.__file__).resolve()), "sha256": file_hash(Path(block20_cells.__file__).resolve())},
            "contour_pilot": {"path": relative(Path(log_pilot.__file__).resolve()), "sha256": file_hash(Path(log_pilot.__file__).resolve())},
            "contour_full": {"path": relative(Path(log_full_builder.__file__).resolve()), "sha256": file_hash(Path(log_full_builder.__file__).resolve())},
        },
        "next_obligation": "Certify the outer Hardy representation and remainder, then replace the finite saved-height coefficient atlas by a height-uniform selector, transition, and accumulation theorem. Keep source-q fidelity on a separate legacy-reproduction branch.",
        "proof_boundary": "Finite saved-height dependency cut and corrected-model majorant only. It does not prove the legacy source q, the outer Hardy representation/remainder, height uniformity, Lambda<=0, PF-infinity, RH, or a prize-level conclusion.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "closed rigorous Arb source-q dependency cut: "
        f"{EXPECTED_CALLS} cells, 4 stripped replays, 15 outputs, max={aggregate['maximum_source_q_free_complete_majorant_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
