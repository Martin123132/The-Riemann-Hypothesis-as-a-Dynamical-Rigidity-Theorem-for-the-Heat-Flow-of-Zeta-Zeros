#!/usr/bin/env python3
"""Validate the finite Hardy chain telemetry route-comparison scout."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_chain_telemetry_scout.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_chain_telemetry_scout.md"
FIXTURE = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10/fixture_result.json"
)
TELEMETRY = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10/enabled/chain_telemetry.jsonl"
)
ANALYZER = REPO_ROOT / "work/rh_compute/scripts/analyze_hardy_chain_telemetry_routes.py"
SOURCE_CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_hardy_chain_telemetry_source_contract.py"
ACCEPTED = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
)
ACCEPTED_SHA256 = "0fb64f090c21b185f27edc1c9194254cf471e74bfd6af3b88cb7f0cb40ec0c3d"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def number(value: object) -> Decimal:
    result = Decimal(str(value))
    require(result.is_finite(), f"non-finite diagnostic value: {value}")
    return result


def require_summary(summary: dict, label: str) -> None:
    require(set(summary) == {"min", "median", "p90", "max"}, f"{label} shape drift")
    values = [number(summary[key]) for key in ("min", "median", "p90", "max")]
    require(values == sorted(values), f"{label} order drift")
    require(values[0] >= 0, f"{label} became negative")


def main() -> int:
    for path in (RESULT, NOTE, FIXTURE, TELEMETRY, ANALYZER, SOURCE_CHECKER, ACCEPTED):
        require(path.is_file(), f"missing telemetry scout artifact: {path}")
    require(file_hash(ACCEPTED) == ACCEPTED_SHA256, "accepted evaluator hash drift")

    result = json.loads(RESULT.read_text(encoding="utf-8"))
    fixture = json.loads(FIXTURE.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")

    require(
        result["kind"] == "jensen_window_pf_newman_c1_hardy_chain_telemetry_scout",
        "result kind drift",
    )
    require(result["status"] == "finite_low_height_route_comparison", "result status drift")
    require(result["precision_ladder_decimal_digits"] == [50, 80, 120], "precision ladder drift")
    require(
        result["phase_constant_provenance"]["accepted_source_definition"]
        == "p=4*ATAN(1); tpp=2*p; tpm=-tpp",
        "phase-constant provenance drift",
    )

    for field in (
        "enabled_vs_disabled_journal_exact",
        "accepted_reference_journal_exact_except_run_id",
        "enabled_vs_disabled_output_exact",
        "accepted_reference_output_exact",
    ):
        require(fixture[field] is True, f"fixture equivalence failure: {field}")
    require(fixture["accepted_source_sha256"] == ACCEPTED_SHA256, "fixture accepted hash drift")
    require(fixture["telemetry"]["chain_count"] == 64, "fixture chain count drift")
    require(fixture["telemetry"]["recurrence_count"] == 64, "fixture recurrence count drift")
    require(
        number(fixture["telemetry"]["max_q_identity_scaled_error"]) < Decimal("1e-32"),
        "fixture q identity rounding guard failed",
    )

    aggregate = result["aggregate"]
    require(aggregate["chain_count"] == 64, "analysis chain count drift")
    require(aggregate["recurrence_count"] == 64, "analysis recurrence count drift")
    require(aggregate["branches"] == [1, 2], "branch inventory drift")
    require(aggregate["initial_length_range"] == [104, 104], "initial length drift")
    require(aggregate["kernel_length_range"] == [1, 2], "kernel length drift")
    require(aggregate["first_parent_shell_term_range"] == [105, 105], "shell cost drift")
    require(aggregate["first_parent_exact_elimination_count"] == 64, "shell identity drift")

    for field in (
        "local_defect_norms",
        "source_root_error_norms",
        "kernel_adapter_error_norms",
        "w1_t5_max_per_chain_norms",
        "w1_to_local_defect_ratios",
        "local_defect_to_root_error_ratios",
        "kernel_to_root_error_ratios",
        "q_component_cancellation_ratios",
        "first_parent_root_residual_norms",
        "first_parent_error_reduction_factors",
    ):
        require_summary(aggregate[field], field)

    require(
        number(aggregate["max_precision_ladder_drift"]) < Decimal("1e-70"),
        "precision ladder did not converge",
    )
    require(
        number(aggregate["max_accumulation_identity_error"]) < Decimal("1e-30"),
        "conjugation-aware accumulation identity drift",
    )
    require(
        number(aggregate["kernel_adapter_error_norms"]["max"])
        < number(aggregate["local_defect_norms"]["min"]),
        "kernel/local-defect dominance observation drift",
    )
    require(
        number(aggregate["w1_to_local_defect_ratios"]["median"]) > Decimal(5),
        "observed W1/direct-defect cancellation signal drift",
    )
    require(
        number(aggregate["local_defect_to_root_error_ratios"]["median"])
        > Decimal("0.95"),
        "first local-defect dominance signal drift",
    )
    require(
        number(aggregate["kernel_to_root_error_ratios"]["max"]) < Decimal("0.01"),
        "kernel error unexpectedly dominates a saved root error",
    )

    chains = result["chains"]
    require([chain["chain"] for chain in chains] == list(range(1, 65)), "chain sequence gap")
    require(sum(chain["branch"] == 1 for chain in chains) == 32, "branch-one count drift")
    require(sum(chain["branch"] == 2 for chain in chains) == 32, "branch-two count drift")
    for chain in chains:
        require(chain["mit"] == 2, f"chain {chain['chain']} MIT drift")
        require(chain["initial_length"] == 104, f"chain {chain['chain']} length drift")
        require(chain["first_parent_shell"] is not None, f"chain {chain['chain']} shell missing")
        require(
            chain["first_parent_reduction_factor"] == "exact_elimination",
            f"chain {chain['chain']} shell elimination drift",
        )
        require(
            number(chain["accumulation_identity_error"]) < Decimal("1e-30"),
            f"chain {chain['chain']} accumulation drift",
        )

    routes = result["route_comparison"]
    require(
        routes["provisional_order"]
        == ["direct_local_defect", "exact_shell", "termwise_w1"],
        "route order drift",
    )
    require(routes["direct_local_defect"]["status"] == "best_next_certificate_candidate", "direct route status drift")
    require(routes["exact_shell"]["exact_elimination_count"] == 64, "shell route count drift")
    require(routes["termwise_w1"]["status"] == "diagnostic_only", "W1 status drift")

    require(file_hash(FIXTURE) == result["evaluator_equivalence"]["fixture_sha256"], "fixture hash drift")
    require(file_hash(TELEMETRY) == result["sources"]["telemetry"]["sha256"], "telemetry hash drift")
    require(file_hash(ANALYZER) == result["sources"]["builder"]["sha256"], "analyzer hash drift")
    require(
        file_hash(Path(__file__).resolve()) == result["sources"]["checker"]["sha256"],
        "checker hash drift",
    )

    for token in (
        "Date: 2026-08-06",
        "Status: finite low-height diagnostic; not a proof and no external error theorem",
        "p = 4*atan(1)",
        "`MIT=2`",
        "does not prove a transformed-level interval adapter",
        "prize-level conclusion",
    ):
        require(token in note, f"note missing proof-boundary token: {token}")

    print(
        "validated Hardy chain telemetry scout: 64 chains, 64 direct defects, "
        "2 branches, exact evaluator state, 3 diagnostic routes and 0 uniform bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
