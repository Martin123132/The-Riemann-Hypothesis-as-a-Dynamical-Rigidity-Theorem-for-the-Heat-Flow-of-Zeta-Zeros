#!/usr/bin/env python3
"""Build explicit asymptotic-remainder targets for the degree-40 collar."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESIDUAL_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.json"
)
DEFAULT_FORMAL_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.md"

SAFE_COMMON_MULTIPLIER = flint.arb(1000)


@dataclass(frozen=True)
class FirstOmittedMultiplierRow:
    index: int
    value_first_omitted_scaled: str
    derivative_first_omitted_scaled: str
    value_budget_multiplier_limit: str
    derivative_budget_multiplier_limit: str
    common_budget_multiplier_limit: str
    safe_multiplier_value_budget_fraction: str
    safe_multiplier_derivative_budget_fraction: str
    limiting_channel: str
    proof_boundary: str


@dataclass(frozen=True)
class OptimizedWindowMultiplierRow:
    index: int
    value_window_sum_j21_to_j120: str
    derivative_window_sum_j21_to_j120: str
    value_window_budget_fraction: str
    derivative_window_budget_fraction: str
    least_value_term_j: int
    least_value_term: str
    least_derivative_term_j: int
    least_derivative_term: str
    remaining_value_budget_after_window: str
    remaining_derivative_budget_after_window: str
    least_value_multiplier_limit_after_window: str
    least_derivative_multiplier_limit_after_window: str
    common_least_multiplier_limit_after_window: str
    proof_boundary: str


def arb_text(value: flint.arb, digits: int = 18) -> str:
    return value.str(digits)


def sci_from_arb(value: flint.arb) -> str:
    return f"{float(value):.18E}"


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_paths(residual_path: Path, formal_path: Path) -> dict[str, str]:
    return {
        "residual_budget_json": residual_path.relative_to(REPO_ROOT).as_posix(),
        "formal_tail_obstruction_json": formal_path.relative_to(REPO_ROOT).as_posix(),
        "residual_budget_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md",
        "formal_tail_obstruction_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md",
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


def build_diagnostics(residual: dict, formal: dict) -> dict:
    residual_summary = residual["summary"]
    formal_diagnostics = formal["matrix_rows"][1]["diagnostics"]
    value_budget = flint.arb(residual_summary["value_residual_half_safety_budget_A"])
    derivative_budget = flint.arb(residual_summary["derivative_residual_half_safety_budget_B"])

    first_rows: list[FirstOmittedMultiplierRow] = []
    optimized_rows: list[OptimizedWindowMultiplierRow] = []
    for row in formal_diagnostics["profile_rows"]:
        value_first = flint.arb(row["value_term_j21"])
        derivative_first = flint.arb(row["derivative_term_j21"])
        value_limit = value_budget / value_first
        derivative_limit = derivative_budget / derivative_first
        common_limit = min(value_limit, derivative_limit)
        value_fraction = SAFE_COMMON_MULTIPLIER * value_first / value_budget
        derivative_fraction = SAFE_COMMON_MULTIPLIER * derivative_first / derivative_budget
        first_rows.append(
            FirstOmittedMultiplierRow(
                index=int(row["index"]),
                value_first_omitted_scaled=arb_text(value_first),
                derivative_first_omitted_scaled=arb_text(derivative_first),
                value_budget_multiplier_limit=arb_text(value_limit),
                derivative_budget_multiplier_limit=arb_text(derivative_limit),
                common_budget_multiplier_limit=arb_text(common_limit),
                safe_multiplier_value_budget_fraction=arb_text(value_fraction),
                safe_multiplier_derivative_budget_fraction=arb_text(derivative_fraction),
                limiting_channel="value" if value_limit < derivative_limit else "derivative",
                proof_boundary=(
                    "Sufficient first-omitted-term multiplier only; the actual remainder bound is not proved."
                ),
            )
        )

        value_window = (
            flint.arb(row["value_sum_j21_to_j40"])
            + flint.arb(row["value_sum_j41_to_j80"])
            + flint.arb(row["value_sum_j81_to_j120"])
        )
        derivative_window = (
            flint.arb(row["derivative_sum_j21_to_j40"])
            + flint.arb(row["derivative_sum_j41_to_j80"])
            + flint.arb(row["derivative_sum_j81_to_j120"])
        )
        value_remaining = value_budget - value_window
        derivative_remaining = derivative_budget - derivative_window
        least_value = flint.arb(row["least_value_term"])
        least_derivative = flint.arb(row["least_derivative_term"])
        least_value_limit = value_remaining / least_value
        least_derivative_limit = derivative_remaining / least_derivative
        optimized_rows.append(
            OptimizedWindowMultiplierRow(
                index=int(row["index"]),
                value_window_sum_j21_to_j120=arb_text(value_window),
                derivative_window_sum_j21_to_j120=arb_text(derivative_window),
                value_window_budget_fraction=arb_text(value_window / value_budget),
                derivative_window_budget_fraction=arb_text(derivative_window / derivative_budget),
                least_value_term_j=int(row["least_value_term_j"]),
                least_value_term=arb_text(least_value),
                least_derivative_term_j=int(row["least_derivative_term_j"]),
                least_derivative_term=arb_text(least_derivative),
                remaining_value_budget_after_window=arb_text(value_remaining),
                remaining_derivative_budget_after_window=arb_text(derivative_remaining),
                least_value_multiplier_limit_after_window=arb_text(least_value_limit),
                least_derivative_multiplier_limit_after_window=arb_text(least_derivative_limit),
                common_least_multiplier_limit_after_window=arb_text(min(least_value_limit, least_derivative_limit)),
                proof_boundary=(
                    "Conditional optimized-window budget only; it assumes a future theorem controls the actual remainder "
                    "after a finite formal window by a least-term multiplier."
                ),
            )
        )

    common_first_limit = min(flint.arb(row.common_budget_multiplier_limit) for row in first_rows)
    max_safe_value_fraction = max(flint.arb(row.safe_multiplier_value_budget_fraction) for row in first_rows)
    max_safe_derivative_fraction = max(flint.arb(row.safe_multiplier_derivative_budget_fraction) for row in first_rows)
    common_least_limit = min(flint.arb(row.common_least_multiplier_limit_after_window) for row in optimized_rows)
    max_window_value_fraction = max(flint.arb(row.value_window_budget_fraction) for row in optimized_rows)
    max_window_derivative_fraction = max(flint.arb(row.derivative_window_budget_fraction) for row in optimized_rows)

    return {
        "parameters": {
            "tail_start_k": 22,
            "indices": [21, 22, 23, 24],
            "collar_interval_u": "[0, 1/1156]",
            "degree40_residual_first_j": 21,
            "formal_window_j_range": "21..120",
            "least_term_j": 103,
            "value_budget_A": residual_summary["value_residual_half_safety_budget_A"],
            "derivative_budget_B": residual_summary["derivative_residual_half_safety_budget_B"],
            "safe_common_first_omitted_multiplier": "1000",
        },
        "first_omitted_multiplier_rows": [asdict(row) for row in first_rows],
        "optimized_window_multiplier_rows": [asdict(row) for row in optimized_rows],
        "first_omitted_row_count": len(first_rows),
        "optimized_window_row_count": len(optimized_rows),
        "common_first_omitted_multiplier_limit": arb_text(common_first_limit),
        "safe_common_multiplier_value_fraction": arb_text(max_safe_value_fraction),
        "safe_common_multiplier_derivative_fraction": arb_text(max_safe_derivative_fraction),
        "safe_common_multiplier_closes_first_omitted_target": bool(
            max_safe_value_fraction < 1 and max_safe_derivative_fraction < 1
        ),
        "common_least_multiplier_limit_after_window": arb_text(common_least_limit),
        "max_window_value_budget_fraction": arb_text(max_window_value_fraction),
        "max_window_derivative_budget_fraction": arb_text(max_window_derivative_fraction),
        "proof_boundary_note": (
            "The constants are sufficient targets only. They do not prove that the actual zeta heat-flow "
            "remainder is bounded by a first omitted term or by a least-term optimized truncation."
        ),
    }


def build_artifact(residual_path: Path, formal_path: Path) -> dict:
    residual = load_json(residual_path)
    formal = load_json(formal_path)
    diagnostics = build_diagnostics(residual, formal)
    paths = source_paths(residual_path, formal_path)
    rows = [
        {
            "id": "nlrgart_01_budget_import",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "Import the degree-40 residual budgets A and B and the formal-tail obstruction profile on 0<=u<=1/1156.",
            "diagnostics": diagnostics,
            "source_artifacts": [
                paths["residual_budget_note"],
                paths["formal_tail_obstruction_note"],
            ],
            "proof_boundary": "Exact budget bookkeeping only; no analytic remainder theorem is proved.",
        },
        {
            "id": "nlrgart_02_first_omitted_multiplier_target",
            "role": "exact_sufficient_condition",
            "readiness": "not_ready_to_apply",
            "claim": "A common first-omitted-term remainder theorem with multiplier at most the recorded common limit would close the fixed-k degree-40 residual target.",
            "formula": "|R_i(u)| <= C*|a_(i,21)|*u^21 and |R_i'(u)| <= C*21*|a_(i,21)|*u^20",
            "proof_boundary": "Sufficient condition only; the first-omitted-term theorem remains open.",
        },
        {
            "id": "nlrgart_03_safe_1000x_first_omitted_budget",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "A 1000x first-omitted-term theorem would use less than the available half-safety value and derivative budgets in all four index rows.",
            "proof_boundary": "Finite target calibration only; it does not prove the 1000x theorem.",
        },
        {
            "id": "nlrgart_04_optimized_window_target",
            "role": "conditional_handoff",
            "readiness": "not_ready_to_apply",
            "claim": "If a future asymptotic theorem justifies summing the formal window j=21..120 and bounding the remaining actual error by a least-term multiplier, the recorded post-window multiplier limits are sufficient.",
            "proof_boundary": "Conditional handoff only; it does not justify the formal window or prove the least-term remainder bound.",
        },
        {
            "id": "nlrgart_05_bad_infinite_sum_rejected",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "claim": "The constants in this target allow promotion of the infinite formal Taylor tail as a convergent absolute series.",
            "gap": "The formal-tail obstruction shows the formal terms turn around after the least-term region; this artifact only calibrates finite/conditional asymptotic remainder targets.",
            "proof_boundary": "Rejected shortcut only; not a proof of any residual estimate.",
        },
        {
            "id": "nlrgart_06_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "A promoted proof must identify the actual analytic remainder, the truncation rule, uniformity on the collar, derivative control, and interval-safe comparison with the recorded constants.",
            "source_artifacts": [
                paths["uniform_remainder_target"],
                paths["dependency_graph"],
            ],
            "proof_boundary": "Proof-hygiene gate only; not scaled-curvature monotonicity, cone entry, RH, or Lambda <= 0.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "first_omitted_rows": diagnostics["first_omitted_row_count"],
        "optimized_window_rows": diagnostics["optimized_window_row_count"],
        "common_first_omitted_multiplier_limit": diagnostics["common_first_omitted_multiplier_limit"],
        "safe_common_multiplier": "1000",
        "safe_common_multiplier_value_fraction": diagnostics["safe_common_multiplier_value_fraction"],
        "safe_common_multiplier_derivative_fraction": diagnostics["safe_common_multiplier_derivative_fraction"],
        "safe_common_multiplier_closes_first_omitted_target": diagnostics[
            "safe_common_multiplier_closes_first_omitted_target"
        ],
        "common_least_multiplier_limit_after_window": diagnostics["common_least_multiplier_limit_after_window"],
        "max_window_value_budget_fraction": diagnostics["max_window_value_budget_fraction"],
        "max_window_derivative_budget_fraction": diagnostics["max_window_derivative_budget_fraction"],
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The degree-40 residual budget can be rephrased as a concrete asymptotic-remainder target: "
            "a common 1000x first-omitted-term theorem would still fit inside the half-safety value and "
            "derivative budgets for F_21..F_24 on 0<=u<=1/1156. A finite optimized-window route through "
            "j=120 leaves still larger least-term multiplier slack, but both routes require a real analytic "
            "remainder theorem rather than an infinite formal-tail sum."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target",
        "date": "2026-07-07",
        "status": "exact theorem-search diagnostic",
        "source_degree40_residual_tail_budget": paths["residual_budget_note"],
        "source_formal_tail_obstruction": paths["formal_tail_obstruction_note"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.py",
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.py"
        ),
        "proof_boundary": (
            "Exact theorem-search diagnostic only. It calibrates sufficient first-omitted-term and least-term "
            "asymptotic remainder targets, but it does not prove those remainder estimates, does not prove "
            "scaled-curvature monotonicity, does not prove cone entry, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "The 1000x first-omitted-term condition is a sufficient target, not a proved theorem.",
            "The optimized-window route is conditional on an actual analytic remainder theorem.",
            "The infinite formal Taylor tail is not promoted to a convergent proof object.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][0]["diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian asymptotic remainder target: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['first_omitted_rows']} first-omitted rows, "
        f"{summary['optimized_window_rows']} optimized-window rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Asymptotic Remainder Target",
        "",
        "Date: 2026-07-07",
        "",
        "Status: exact theorem-search diagnostic. This is not a proof",
        "of an analytic residual estimate, scaled-curvature monotonicity,",
        "cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target`.",
        "",
        "Proof boundary: this artifact converts the degree-40 residual budget",
        "and the formal-tail obstruction into explicit sufficient asymptotic",
        "remainder targets. It does not prove those targets.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## First-Omitted-Term Target",
        "",
        "```text",
        f"common multiplier limit: {summary['common_first_omitted_multiplier_limit']}",
        f"safe common multiplier tested: {summary['safe_common_multiplier']}",
        f"safe value budget fraction: {summary['safe_common_multiplier_value_fraction']}",
        f"safe derivative budget fraction: {summary['safe_common_multiplier_derivative_fraction']}",
        f"safe multiplier closes target: {summary['safe_common_multiplier_closes_first_omitted_target']}",
        "```",
        "",
        "Per-index first-omitted rows:",
        "",
        "```text",
    ]
    for row in diagnostics["first_omitted_multiplier_rows"]:
        lines.append(
            f"F_{row['index']}: value term={row['value_first_omitted_scaled']}, "
            f"derivative term={row['derivative_first_omitted_scaled']}, "
            f"common multiplier limit={row['common_budget_multiplier_limit']}, "
            f"1000x fractions=({row['safe_multiplier_value_budget_fraction']}, "
            f"{row['safe_multiplier_derivative_budget_fraction']}), "
            f"limiting={row['limiting_channel']}"
        )
    lines.extend(
        [
            "```",
            "",
            "## Optimized-Window Target",
            "",
            "This conditional route is different from summing the infinite formal",
            "tail: it would need a theorem justifying a finite formal window and",
            "an actual remainder bound after truncation near the least term.",
            "",
            "```text",
            f"formal window: {diagnostics['parameters']['formal_window_j_range']}",
            f"least term j: {diagnostics['parameters']['least_term_j']}",
            f"max value window budget fraction: {summary['max_window_value_budget_fraction']}",
            f"max derivative window budget fraction: {summary['max_window_derivative_budget_fraction']}",
            f"common least-term multiplier limit after window: {summary['common_least_multiplier_limit_after_window']}",
            "```",
            "",
            "Per-index optimized-window rows:",
            "",
            "```text",
        ]
    )
    for row in diagnostics["optimized_window_multiplier_rows"]:
        lines.append(
            f"F_{row['index']}: window fractions=({row['value_window_budget_fraction']}, "
            f"{row['derivative_window_budget_fraction']}), least terms j=({row['least_value_term_j']}, "
            f"{row['least_derivative_term_j']}), common least-term multiplier="
            f"{row['common_least_multiplier_limit_after_window']}"
        )
    lines.extend(
        [
            "```",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_degree40_residual_tail_budget"],
            artifact["source_formal_tail_obstruction"],
            artifact["source_uniform_remainder_target"],
            artifact["source_dependency_graph"],
            "```",
            "",
            "Summary:",
            "",
            summary["main_finding"],
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--residual-json", type=Path, default=DEFAULT_RESIDUAL_JSON)
    parser.add_argument("--formal-json", type=Path, default=DEFAULT_FORMAL_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    residual_path = args.residual_json if args.residual_json.is_absolute() else REPO_ROOT / args.residual_json
    formal_path = args.formal_json if args.formal_json.is_absolute() else REPO_ROOT / args.formal_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(residual_path, formal_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian asymptotic remainder target: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
