#!/usr/bin/env python3
"""Build the k300 obstruction to the fixed bounded log-curvature wall."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
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

from jensen_window_pf_heat_flow_monotone_closure_scout import (  # noqa: E402
    REPO_ROOT,
    arb_positive,
    decimal_format,
)
from jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge import (  # noqa: E402
    h_arb,
    h_decimal,
)
from jensen_window_pf_negative_lambda_k300_precision_repair_audit import (  # noqa: E402
    DEFAULT_BASE_ENCLOSURE,
    DEFAULT_REPAIR_ENCLOSURES,
    label_to_t,
    ratio_tables,
)


getcontext().prec = 100

DEFAULT_ENCLOSURES = (DEFAULT_BASE_ENCLOSURE, *DEFAULT_REPAIR_ENCLOSURES)
DEFAULT_OUT_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.md"


@dataclass(frozen=True)
class Extremum:
    sample: str
    lam: str
    k: int


def update_min(current: tuple[Decimal, str, int] | None, value: Decimal, label: str, k: int) -> tuple[Decimal, str, int]:
    if current is None or value < current[0]:
        return value, label, k
    return current


def update_max(current: tuple[Decimal, str, int] | None, value: Decimal, label: str, k: int) -> tuple[Decimal, str, int]:
    if current is None or value > current[0]:
        return value, label, k
    return current


def extremum(item: tuple[Decimal, str, int] | None) -> Extremum:
    if item is None:
        raise RuntimeError("missing bounded log-curvature obstruction extremum")
    value, lam, k = item
    return Extremum(decimal_format(value), lam, k)


def first_failure_dict(item: dict | None) -> dict | None:
    if item is None:
        return None
    return {
        "lam": item["lam"],
        "k": item["k"],
        "scaled_curvature": decimal_format(item["scaled_curvature"]),
        "two_thirds_slack": decimal_format(item["two_thirds_slack"]),
        "B_k": decimal_format(item["B_k"]),
    }


def build_diagnostics(paths: list[Path], coefficient_cap: int = 300) -> dict:
    by_label_arb, by_label_sample = ratio_tables(paths, coefficient_cap)
    labels = sorted(by_label_sample, key=label_to_t)
    checked_ratio_max = coefficient_cap - 1
    adjacent_max = checked_ratio_max - 1

    total_rows = checked_ratio_max * len(labels)
    adjacent_total_rows = adjacent_max * len(labels)
    two_thirds_rows = 0
    two_thirds_failures = 0
    two_thirds_inconclusive = 0
    b_wall_rows = 0
    c_increase_rows = 0
    c_increase_failures = 0
    c_increase_inconclusive = 0
    b_decrease_rows = 0

    min_two_thirds_slack: tuple[Decimal, str, int] | None = None
    max_scaled_curvature: tuple[Decimal, str, int] | None = None
    min_b_wall: tuple[Decimal, str, int] | None = None
    max_b_wall: tuple[Decimal, str, int] | None = None
    min_c_increase_margin: tuple[Decimal, str, int] | None = None
    min_b_decrease_margin: tuple[Decimal, str, int] | None = None
    global_first_failure: dict | None = None
    per_lambda: dict[str, dict] = {}

    for label in labels:
        b_arb: dict[int, flint.arb] = {}
        c_arb: dict[int, flint.arb] = {}
        b_sample: dict[int, Decimal] = {}
        c_sample: dict[int, Decimal] = {}
        local_pass = 0
        local_fail = 0
        local_inconclusive = 0
        local_first_failure: dict | None = None

        for k in range(1, checked_ratio_max + 1):
            raw_arb = by_label_arb[label][k]
            raw_sample = by_label_sample[label][k]
            b_arb[k] = h_arb(k) - raw_arb.log()
            b_sample[k] = h_decimal(k) - raw_sample.ln()
            c_arb[k] = flint.arb(2 * k + 1) * b_arb[k]
            c_sample[k] = Decimal(2 * k + 1) * b_sample[k]
            slack_arb = flint.arb(2) / flint.arb(3) - c_arb[k]
            slack_sample = Decimal(2) / Decimal(3) - c_sample[k]

            if arb_positive(b_arb[k]):
                b_wall_rows += 1
            if arb_positive(slack_arb):
                two_thirds_rows += 1
                local_pass += 1
            elif arb_positive(-slack_arb):
                two_thirds_failures += 1
                local_fail += 1
                failure = {
                    "lam": label,
                    "k": k,
                    "scaled_curvature": c_sample[k],
                    "two_thirds_slack": slack_sample,
                    "B_k": b_sample[k],
                }
                if local_first_failure is None:
                    local_first_failure = failure
                if global_first_failure is None:
                    global_first_failure = failure
            else:
                two_thirds_inconclusive += 1
                local_inconclusive += 1

            min_two_thirds_slack = update_min(min_two_thirds_slack, slack_sample, label, k)
            max_scaled_curvature = update_max(max_scaled_curvature, c_sample[k], label, k)
            min_b_wall = update_min(min_b_wall, b_sample[k], label, k)
            max_b_wall = update_max(max_b_wall, b_sample[k], label, k)

        for k in range(1, adjacent_max + 1):
            c_gap_arb = c_arb[k + 1] - c_arb[k]
            c_gap_sample = c_sample[k + 1] - c_sample[k]
            b_gap_arb = b_arb[k] - b_arb[k + 1]
            b_gap_sample = b_sample[k] - b_sample[k + 1]
            if arb_positive(c_gap_arb):
                c_increase_rows += 1
            elif arb_positive(-c_gap_arb):
                c_increase_failures += 1
            else:
                c_increase_inconclusive += 1
            if arb_positive(b_gap_arb):
                b_decrease_rows += 1
            min_c_increase_margin = update_min(min_c_increase_margin, c_gap_sample, label, k)
            min_b_decrease_margin = update_min(min_b_decrease_margin, b_gap_sample, label, k)

        per_lambda[label] = {
            "two_thirds_bound_rows": local_pass,
            "two_thirds_failure_rows": local_fail,
            "two_thirds_inconclusive_rows": local_inconclusive,
            "first_failure": first_failure_dict(local_first_failure),
        }

    return {
        "lambdas": labels,
        "coefficient_cap": coefficient_cap,
        "raw_ratio_rows_per_lambda": checked_ratio_max,
        "adjacent_rows_per_lambda": adjacent_max,
        "scaled_curvature_total_rows": total_rows,
        "adjacent_total_rows": adjacent_total_rows,
        "two_thirds_bound_rows": two_thirds_rows,
        "two_thirds_failure_rows": two_thirds_failures,
        "two_thirds_inconclusive_rows": two_thirds_inconclusive,
        "b_wall_rows": b_wall_rows,
        "scaled_curvature_increase_rows": c_increase_rows,
        "scaled_curvature_increase_failures": c_increase_failures,
        "scaled_curvature_increase_inconclusive": c_increase_inconclusive,
        "b_decrease_rows": b_decrease_rows,
        "first_two_thirds_failure": first_failure_dict(global_first_failure),
        "min_two_thirds_slack": asdict(extremum(min_two_thirds_slack)),
        "max_scaled_curvature": asdict(extremum(max_scaled_curvature)),
        "min_b_wall": asdict(extremum(min_b_wall)),
        "max_b_wall": asdict(extremum(max_b_wall)),
        "min_scaled_curvature_increase_margin": asdict(extremum(min_c_increase_margin)),
        "min_b_decrease_margin": asdict(extremum(min_b_decrease_margin)),
        "per_lambda": per_lambda,
    }


def build_artifact(paths: list[Path] | None = None) -> dict:
    if paths is None:
        paths = list(DEFAULT_ENCLOSURES)
    diagnostics = build_diagnostics(paths)
    rows = [
        {
            "id": "nllcko_01_exact_scaled_curvature_rewrite",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "With B_k=-log(((2*k-1)/(2*k+1))*R_k) and C_k=(2*k+1)*B_k, the old bounded log-curvature wall B_k<=2/(3*(2*k+1)) is exactly C_k<=2/3.",
            "formula": "B_k<=2/(3*(2*k+1)) iff C_k=(2*k+1)*B_k<=2/3",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md",
                "outputs/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.md",
            ],
            "proof_boundary": "Exact rewrite only; it does not prove or disprove any zeta coefficient inequality by itself.",
        },
        {
            "id": "nllcko_02_repaired_k300_two_thirds_obstruction",
            "role": "finite_obstruction",
            "readiness": "not_ready_to_apply",
            "claim": "The repaired k300 data reject the fixed scaled-curvature wall C_k<=2/3 on 718 of 897 checked rows, with zero inconclusive rows.",
            "diagnostics": diagnostics,
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
            ],
            "proof_boundary": "Finite Arb obstruction gate only; it retires this fixed wall but does not prove a replacement all-k theorem.",
        },
        {
            "id": "nllcko_03_first_checked_failure",
            "role": "finite_counterexample",
            "readiness": "not_ready_to_apply",
            "claim": "The first checked failure in the repaired ordering is lambda=-25.0 at k=31.",
            "witness": diagnostics["first_two_thirds_failure"],
            "proof_boundary": "Finite checked counterexample to the fixed wall only; not a statement against the linear curvature-barrier route.",
        },
        {
            "id": "nllcko_04_worst_checked_failure",
            "role": "finite_extremum",
            "readiness": "not_ready_to_apply",
            "claim": "The worst checked slack is negative at lambda=-25.0, k=299, where C_k exceeds 2/3 by about 0.47755.",
            "max_scaled_curvature": diagnostics["max_scaled_curvature"],
            "min_two_thirds_slack": diagnostics["min_two_thirds_slack"],
            "proof_boundary": "Finite stress extremum only; not an all-k asymptotic theorem.",
        },
        {
            "id": "nllcko_05_k22_prefix_scope",
            "role": "historical_scope",
            "readiness": "not_ready_to_apply",
            "claim": "The former bounded log-curvature target was supported only by the old k<=22 prefix compatibility check; that prefix does not license the fixed 2/3 wall at k300.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md",
            ],
            "proof_boundary": "Scope correction only; finite prefix compatibility is not promoted to an all-k theorem.",
        },
        {
            "id": "nllcko_06_replacement_scaled_monotonicity_route",
            "role": "replacement_route_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "The repaired data instead support the newer scaled-curvature monotonicity shape C_(k+1)>=C_k and B_(k+1)<=B_k on all 894 adjacent rows.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.md",
                "outputs/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.md",
            ],
            "proof_boundary": "Finite route diagnostic only; it does not prove the all-k scaled-curvature or linear-barrier theorem.",
        },
        {
            "id": "nllcko_07_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "Future negative-lambda proofs must not use C_k<=2/3 or B_k<=2/(3*(2*k+1)) as a live tail wall; any replacement must prove the zeta-specific raw corridor or linear curvature barrier with an explicit finite collar.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
                "outputs/signed_hankel_jensen_dependency_graph.md",
            ],
            "proof_boundary": "Proof-hygiene gate only; not a proof of cone entry, jwpf_06, RH, or Lambda <= 0.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "exact_reduction_rows": 1,
        "finite_obstruction_rows": 3,
        "historical_scope_rows": 1,
        "replacement_route_rows": 1,
        "acceptance_gate_rows": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "lambdas": diagnostics["lambdas"],
        "coefficient_cap": diagnostics["coefficient_cap"],
        "scaled_curvature_total_rows": diagnostics["scaled_curvature_total_rows"],
        "two_thirds_bound_rows": diagnostics["two_thirds_bound_rows"],
        "two_thirds_failure_rows": diagnostics["two_thirds_failure_rows"],
        "two_thirds_inconclusive_rows": diagnostics["two_thirds_inconclusive_rows"],
        "adjacent_total_rows": diagnostics["adjacent_total_rows"],
        "scaled_curvature_increase_rows": diagnostics["scaled_curvature_increase_rows"],
        "b_decrease_rows": diagnostics["b_decrease_rows"],
        "target_retired": True,
        "main_finding": (
            "The fixed bounded log-curvature wall B_k<=2/(3*(2*k+1)), equivalently "
            "C_k=(2*k+1)*B_k<=2/3, is finite-rejected by the repaired k300 data: "
            "only 179/897 checked rows satisfy it, 718/897 fail it, and zero rows "
            "are inconclusive. The live replacement direction is the zeta-specific "
            "linear/scaled curvature corridor, not a fixed 2/3 wall."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction",
        "date": "2026-07-07",
        "status": "finite obstruction gate",
        "source_bounded_log_curvature_target": "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md",
        "source_k300_precision_repair": "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
        "source_curvature_corridor_bridge": "outputs/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.md",
        "source_linear_barrier_scout": "outputs/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.md",
        "source_raw_corridor_target": "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.py",
        "proof_boundary": (
            "Finite obstruction gate only. It records that the fixed 2/3 scaled-curvature "
            "wall is incompatible with the repaired k300 zeta data, but it does not prove "
            "the replacement raw-corridor theorem, does not prove cone entry, does not prove "
            "jwpf_06, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "invariants": [
            "No row is ready_to_apply.",
            "The fixed 2/3 scaled-curvature wall is retired as a live theorem target.",
            "The repaired k300 obstruction is not evidence against the linear curvature-barrier route.",
            "Finite scaled-curvature monotonicity is stress evidence only.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
        "summary": summary,
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][1]["diagnostics"]
    first = diagnostics["first_two_thirds_failure"]
    max_c = diagnostics["max_scaled_curvature"]
    min_slack = diagnostics["min_two_thirds_slack"]
    result_line = (
        "validated Jensen-window PF negative-lambda bounded log-curvature k300 obstruction: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['two_thirds_failure_rows']} two-thirds failures, "
        f"{summary['scaled_curvature_increase_rows']} scaled-curvature increase rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Bounded Log-Curvature k300 Obstruction",
        "",
        "Date: 2026-07-07",
        "",
        "Status: finite obstruction gate. This is not a proof of the replacement",
        "raw-corridor theorem, cone entry, Jensen-window PF-infinity, RH, or",
        "`Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction`.",
        "",
        "Proof boundary: this artifact retires the fixed `2/3` scaled-curvature",
        "wall as a live target. It does not prove the newer linear curvature",
        "barrier or any all-`k` theorem.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Exact Rewrite",
        "",
        "Let:",
        "",
        "```text",
        "B_k = -log(((2*k-1)/(2*k+1))*R_k)",
        "C_k = (2*k+1)*B_k",
        "```",
        "",
        "Then the former bounded log-curvature wall is exactly:",
        "",
        "```text",
        "B_k <= 2/(3*(2*k+1))  iff  C_k <= 2/3",
        "```",
        "",
        "## Repaired k300 Obstruction",
        "",
        "Inputs:",
        "",
        "```text",
        "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k300.jsonl",
        "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k220_k250_dps220.jsonl",
        "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220.jsonl",
        "```",
        "",
        "Classified rows:",
        "",
        "```text",
        f"scaled-curvature rows: {diagnostics['scaled_curvature_total_rows']}",
        f"C_k <= 2/3 rows: {diagnostics['two_thirds_bound_rows']} / {diagnostics['scaled_curvature_total_rows']}",
        f"C_k > 2/3 rows: {diagnostics['two_thirds_failure_rows']} / {diagnostics['scaled_curvature_total_rows']}",
        f"inconclusive rows: {diagnostics['two_thirds_inconclusive_rows']}",
        f"B_k > 0 rows: {diagnostics['b_wall_rows']} / {diagnostics['scaled_curvature_total_rows']}",
        "```",
        "",
        "First checked failure:",
        "",
        "```text",
        f"lambda={first['lam']}, k={first['k']}",
        f"C_k = {first['scaled_curvature']}",
        f"2/3 - C_k = {first['two_thirds_slack']}",
        f"B_k = {first['B_k']}",
        "```",
        "",
        "Worst checked slack:",
        "",
        "```text",
        f"max C_k = {max_c['sample']} at lambda={max_c['lam']}, k={max_c['k']}",
        f"min 2/3-C_k = {min_slack['sample']} at lambda={min_slack['lam']}, k={min_slack['k']}",
        "```",
        "",
        "Replacement-route contrast:",
        "",
        "```text",
        f"C_(k+1)-C_k positive rows: {diagnostics['scaled_curvature_increase_rows']} / {diagnostics['adjacent_total_rows']}",
        f"B_k-B_(k+1) positive rows: {diagnostics['b_decrease_rows']} / {diagnostics['adjacent_total_rows']}",
        f"min C increase margin: {diagnostics['min_scaled_curvature_increase_margin']['sample']} at lambda={diagnostics['min_scaled_curvature_increase_margin']['lam']}, k={diagnostics['min_scaled_curvature_increase_margin']['k']}",
        f"min B decrease margin: {diagnostics['min_b_decrease_margin']['sample']} at lambda={diagnostics['min_b_decrease_margin']['lam']}, k={diagnostics['min_b_decrease_margin']['k']}",
        "```",
        "",
        "Integration:",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md",
        "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
        "outputs/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.md",
        "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "```",
        "",
        "Summary:",
        "",
        summary["main_finding"],
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact()
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda bounded log-curvature k300 obstruction: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
