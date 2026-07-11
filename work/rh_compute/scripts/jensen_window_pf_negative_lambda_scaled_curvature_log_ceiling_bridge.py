#!/usr/bin/env python3
"""Build a log-ceiling bridge for scaled-curvature monotonicity."""

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
from jensen_window_pf_negative_lambda_raw_log_decrement_bridge import (  # noqa: E402
    log_lower_bound_arb,
    log_lower_bound_decimal,
    log_upper_bound_arb,
    log_upper_bound_decimal,
)


getcontext().prec = 100

DEFAULT_ENCLOSURES = (DEFAULT_BASE_ENCLOSURE, *DEFAULT_REPAIR_ENCLOSURES)
DEFAULT_OUT_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.md"


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
        raise RuntimeError("missing scaled-curvature log-ceiling extremum")
    value, lam, k = item
    return Extremum(decimal_format(value), lam, k)


def alpha_arb(k: int) -> flint.arb:
    return flint.arb(2 * k + 1) / flint.arb(2 * k + 3)


def alpha_decimal(k: int) -> Decimal:
    return Decimal(2 * k + 1) / Decimal(2 * k + 3)


def scaled_log_ceiling_arb(k: int, p_value: flint.arb) -> flint.arb:
    return h_arb(k + 1) - alpha_arb(k) * h_arb(k) - flint.arb(2) * p_value / flint.arb(2 * k + 3)


def scaled_log_ceiling_decimal(k: int, p_value: Decimal) -> Decimal:
    return h_decimal(k + 1) - alpha_decimal(k) * h_decimal(k) - Decimal(2) * p_value / Decimal(2 * k + 3)


def build_diagnostics(paths: list[Path], coefficient_cap: int = 300) -> dict:
    by_label_arb, by_label_sample = ratio_tables(paths, coefficient_cap)
    labels = sorted(by_label_sample, key=label_to_t)
    checked_ratio_max = coefficient_cap - 1
    adjacent_max = checked_ratio_max - 1

    p_lower_rows = 0
    p_upper_rows = 0
    scaled_ceiling_rows = 0
    raw_upper_rows = 0
    ceiling_dominance_rows = 0
    log_lower_rows = 0
    scaled_log_corridor_rows = 0
    scaled_width_rows = 0

    min_p_lower_margin: tuple[Decimal, str, int] | None = None
    min_p_upper_margin: tuple[Decimal, str, int] | None = None
    min_scaled_ceiling_margin: tuple[Decimal, str, int] | None = None
    min_raw_upper_margin: tuple[Decimal, str, int] | None = None
    min_ceiling_to_raw_upper_slack: tuple[Decimal, str, int] | None = None
    min_log_lower_margin: tuple[Decimal, str, int] | None = None
    min_scaled_log_width: tuple[Decimal, str, int] | None = None
    min_delta: tuple[Decimal, str, int] | None = None
    max_delta: tuple[Decimal, str, int] | None = None
    min_ceiling: tuple[Decimal, str, int] | None = None
    max_ceiling: tuple[Decimal, str, int] | None = None

    for label in labels:
        p_arb: dict[int, flint.arb] = {}
        p_sample: dict[int, Decimal] = {}
        for k in range(1, checked_ratio_max + 1):
            p_arb[k] = by_label_arb[label][k].log()
            p_sample[k] = by_label_sample[label][k].ln()
            p_upper_margin_arb = h_arb(k) - p_arb[k]
            p_upper_margin_sample = h_decimal(k) - p_sample[k]
            if arb_positive(p_arb[k]):
                p_lower_rows += 1
            if arb_positive(p_upper_margin_arb):
                p_upper_rows += 1
            min_p_lower_margin = update_min(min_p_lower_margin, p_sample[k], label, k)
            min_p_upper_margin = update_min(min_p_upper_margin, p_upper_margin_sample, label, k)

        for k in range(1, adjacent_max + 1):
            delta_arb = (by_label_arb[label][k + 1] / by_label_arb[label][k]).log()
            delta_sample = (by_label_sample[label][k + 1] / by_label_sample[label][k]).ln()
            ceiling_arb = scaled_log_ceiling_arb(k, p_arb[k])
            ceiling_sample = scaled_log_ceiling_decimal(k, p_sample[k])
            raw_upper_arb = log_upper_bound_arb(k, by_label_arb[label][k])
            raw_upper_sample = log_upper_bound_decimal(k, by_label_sample[label][k])
            raw_lower_arb = log_lower_bound_arb(k)
            raw_lower_sample = log_lower_bound_decimal(k)

            scaled_margin_arb = ceiling_arb - delta_arb
            raw_upper_margin_arb = raw_upper_arb - delta_arb
            dominance_arb = raw_upper_arb - ceiling_arb
            lower_margin_arb = delta_arb - raw_lower_arb
            width_arb = ceiling_arb - raw_lower_arb

            scaled_margin_sample = ceiling_sample - delta_sample
            raw_upper_margin_sample = raw_upper_sample - delta_sample
            dominance_sample = raw_upper_sample - ceiling_sample
            lower_margin_sample = delta_sample - raw_lower_sample
            width_sample = ceiling_sample - raw_lower_sample

            if arb_positive(scaled_margin_arb):
                scaled_ceiling_rows += 1
            if arb_positive(raw_upper_margin_arb):
                raw_upper_rows += 1
            if arb_positive(dominance_arb):
                ceiling_dominance_rows += 1
            if arb_positive(lower_margin_arb):
                log_lower_rows += 1
            if arb_positive(scaled_margin_arb) and arb_positive(lower_margin_arb):
                scaled_log_corridor_rows += 1
            if arb_positive(width_arb):
                scaled_width_rows += 1

            min_scaled_ceiling_margin = update_min(min_scaled_ceiling_margin, scaled_margin_sample, label, k)
            min_raw_upper_margin = update_min(min_raw_upper_margin, raw_upper_margin_sample, label, k)
            min_ceiling_to_raw_upper_slack = update_min(min_ceiling_to_raw_upper_slack, dominance_sample, label, k)
            min_log_lower_margin = update_min(min_log_lower_margin, lower_margin_sample, label, k)
            min_scaled_log_width = update_min(min_scaled_log_width, width_sample, label, k)
            min_delta = update_min(min_delta, delta_sample, label, k)
            max_delta = update_max(max_delta, delta_sample, label, k)
            min_ceiling = update_min(min_ceiling, ceiling_sample, label, k)
            max_ceiling = update_max(max_ceiling, ceiling_sample, label, k)

    return {
        "lambdas": labels,
        "coefficient_cap": coefficient_cap,
        "raw_ratio_rows_per_lambda": checked_ratio_max,
        "adjacent_rows_per_lambda": adjacent_max,
        "raw_log_total_rows": checked_ratio_max * len(labels),
        "adjacent_total_rows": adjacent_max * len(labels),
        "p_lower_rows": p_lower_rows,
        "p_upper_rows": p_upper_rows,
        "scaled_ceiling_rows": scaled_ceiling_rows,
        "raw_upper_rows": raw_upper_rows,
        "ceiling_dominance_rows": ceiling_dominance_rows,
        "log_lower_rows": log_lower_rows,
        "scaled_log_corridor_rows": scaled_log_corridor_rows,
        "scaled_width_rows": scaled_width_rows,
        "min_p_lower_margin": asdict(extremum(min_p_lower_margin)),
        "min_p_upper_margin": asdict(extremum(min_p_upper_margin)),
        "min_scaled_ceiling_margin": asdict(extremum(min_scaled_ceiling_margin)),
        "min_raw_upper_margin": asdict(extremum(min_raw_upper_margin)),
        "min_ceiling_to_raw_upper_slack": asdict(extremum(min_ceiling_to_raw_upper_slack)),
        "min_log_lower_margin": asdict(extremum(min_log_lower_margin)),
        "min_scaled_log_width": asdict(extremum(min_scaled_log_width)),
        "min_delta": asdict(extremum(min_delta)),
        "max_delta": asdict(extremum(max_delta)),
        "min_ceiling": asdict(extremum(min_ceiling)),
        "max_ceiling": asdict(extremum(max_ceiling)),
    }


def build_artifact(paths: list[Path] | None = None) -> dict:
    if paths is None:
        paths = list(DEFAULT_ENCLOSURES)
    diagnostics = build_diagnostics(paths)
    rows = [
        {
            "id": "nlsclcb_01_log_curvature_coordinate",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "Set p_k=log(R_k), delta_k=p_(k+1)-p_k, h_k=log((2*k+1)/(2*k-1)), B_k=h_k-p_k, and C_k=(2*k+1)*B_k.",
            "formula": "C_k=(2*k+1)*(h_k-p_k)",
            "proof_boundary": "Exact change of variables only; no zeta inequality is proved.",
        },
        {
            "id": "nlsclcb_02_affine_ceiling_equivalence",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The scaled-curvature monotonicity target C_(k+1)>=C_k is exactly an affine upper ceiling for delta_k.",
            "formula": "delta_k <= h_(k+1)-((2*k+1)/(2*k+3))*h_k-2*p_k/(2*k+3)",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.md",
            ],
            "proof_boundary": "Exact algebraic equivalence only; the affine ceiling remains open for all k.",
        },
        {
            "id": "nlsclcb_03_ceiling_dominates_raw_log_upper",
            "role": "exact_sufficient_condition",
            "readiness": "available_exact",
            "claim": "On the raw upper wall B_k>=0, the affine scaled-curvature ceiling is stronger than the nonlinear raw-log upper corridor wall.",
            "formula": "raw_log_upper(p_k)-scaled_ceiling(p_k)=alpha_k*B_k-L_k(B_k)>=0",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.md",
                "outputs/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.md",
            ],
            "proof_boundary": "Exact sufficient comparison only; it does not prove the affine ceiling for zeta coefficients.",
        },
        {
            "id": "nlsclcb_04_repaired_k300_affine_ceiling",
            "role": "finite_stress",
            "readiness": "not_ready_to_apply",
            "claim": "The repaired k300 data validate the p-wall, affine scaled-curvature ceiling, lower log wall, and scaled log corridor on all checked rows.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
                "outputs/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.md",
            ],
            "diagnostics": diagnostics,
            "proof_boundary": "Finite repaired stress evidence only; not an all-k affine log-ceiling theorem.",
        },
        {
            "id": "nlsclcb_05_sharpness_warning",
            "role": "finite_sharpness",
            "readiness": "not_ready_to_apply",
            "claim": "The affine ceiling is very close to the nonlinear raw-log upper wall at the high-k edge, so a promoted proof needs sharp constants.",
            "min_ceiling_to_raw_upper_slack": diagnostics["min_ceiling_to_raw_upper_slack"],
            "proof_boundary": "Finite sharpness diagnostic only; not an asymptotic theorem.",
        },
        {
            "id": "nlsclcb_06_scaled_log_corridor_target",
            "role": "open_dependency",
            "readiness": "not_ready_to_apply",
            "claim": "A sufficient log-recurrence target is lower_log_bound(k)<=delta_k<=scaled_ceiling(k,p_k), plus the p-wall and the companion B-monotone upper side.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.md",
                "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
            ],
            "proof_boundary": "Open dependency only; it is not a substitute for the companion monotone-contraction theorem.",
        },
        {
            "id": "nlsclcb_07_live_affine_recurrence_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim": "A zeta-specific proof may target the affine delta ceiling directly as a log-ratio recurrence after a finite collar.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
                "outputs/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.md",
            ],
            "proof_boundary": "Live theorem-search route only; no all-k affine recurrence is proved here.",
        },
        {
            "id": "nlsclcb_08_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "A promoted proof must state the tail start, lambda range, finite collar, p-wall, lower log wall, affine ceiling, companion B-monotone side, and forbidden assumptions.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
                "outputs/signed_hankel_jensen_dependency_graph.md",
            ],
            "proof_boundary": "Proof-hygiene gate only; not a proof of scaled-curvature monotonicity, raw-corridor occupancy, cone entry, jwpf_06, RH, or Lambda <= 0.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "exact_available_rows": 3,
        "finite_stress_rows": 1,
        "finite_sharpness_rows": 1,
        "open_dependency_rows": 1,
        "live_routes": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "lambdas": diagnostics["lambdas"],
        "coefficient_cap": diagnostics["coefficient_cap"],
        "raw_log_total_rows": diagnostics["raw_log_total_rows"],
        "adjacent_total_rows": diagnostics["adjacent_total_rows"],
        "p_lower_rows": diagnostics["p_lower_rows"],
        "p_upper_rows": diagnostics["p_upper_rows"],
        "scaled_ceiling_rows": diagnostics["scaled_ceiling_rows"],
        "raw_upper_rows": diagnostics["raw_upper_rows"],
        "ceiling_dominance_rows": diagnostics["ceiling_dominance_rows"],
        "log_lower_rows": diagnostics["log_lower_rows"],
        "scaled_log_corridor_rows": diagnostics["scaled_log_corridor_rows"],
        "scaled_width_rows": diagnostics["scaled_width_rows"],
        "main_finding": (
            "Scaled-curvature monotonicity is exactly equivalent to the affine log-ratio ceiling "
            "delta_k<=h_(k+1)-((2*k+1)/(2*k+3))*h_k-2*p_k/(2*k+3). On repaired k300 data "
            "the p-wall holds on 897/897 rows, the affine ceiling holds on 894/894 adjacent rows, "
            "and the ceiling is stronger than the nonlinear raw-log upper wall on 894/894 rows."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge",
        "date": "2026-07-07",
        "status": "exact finite theorem-search diagnostic",
        "source_scaled_curvature_target": "outputs/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.md",
        "source_raw_log_bridge": "outputs/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.md",
        "source_linear_barrier_scout": "outputs/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.md",
        "source_raw_corridor_target": "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
        "source_monotone_contraction_target": "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.py",
        "proof_boundary": (
            "Exact finite theorem-search diagnostic only. It rewrites scaled-curvature "
            "monotonicity as an affine log-ratio ceiling and checks repaired k300 evidence, "
            "but it does not prove the all-k affine recurrence, does not prove the companion "
            "monotone-contraction theorem, does not prove cone entry, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "The affine log-ratio ceiling remains an open zeta-specific all-k theorem.",
            "Repaired k300 evidence is finite stress evidence only.",
            "The nonlinear raw-log upper wall is weaker than the affine ceiling.",
            "The companion monotone-contraction upper side remains a separate open theorem target.",
            "Generic moment positivity, endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][3]["diagnostics"]
    sharpness = artifact["matrix_rows"][4]["min_ceiling_to_raw_upper_slack"]
    result_line = (
        "validated Jensen-window PF negative-lambda scaled-curvature log-ceiling bridge: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['scaled_ceiling_rows']} scaled-ceiling rows, "
        f"{summary['scaled_log_corridor_rows']} scaled-log-corridor rows, "
        f"{summary['ceiling_dominance_rows']} ceiling-dominance rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Scaled-Curvature Log-Ceiling Bridge",
        "",
        "Date: 2026-07-07",
        "",
        "Status: exact finite theorem-search diagnostic. This is not a proof",
        "of scaled-curvature monotonicity, the raw-corridor theorem, cone entry,",
        "RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge`.",
        "",
        "Proof boundary: this artifact rewrites scaled-curvature monotonicity",
        "as an affine log-ratio ceiling and checks repaired k300 evidence. It",
        "does not prove any all-`k` theorem.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Exact Affine Ceiling",
        "",
        "Let:",
        "",
        "```text",
        "p_k = log(R_k)",
        "delta_k = p_(k+1)-p_k",
        "h_k = log((2*k+1)/(2*k-1))",
        "B_k = h_k-p_k",
        "C_k = (2*k+1)*B_k",
        "```",
        "",
        "Then:",
        "",
        "```text",
        "C_(k+1)>=C_k",
        "<=> delta_k <= h_(k+1)-((2*k+1)/(2*k+3))*h_k-2*p_k/(2*k+3)",
        "```",
        "",
        "On the raw upper wall `B_k>=0`, this affine ceiling is stronger than",
        "the nonlinear raw-log upper wall:",
        "",
        "```text",
        "raw_log_upper(p_k)-scaled_ceiling(p_k)=alpha_k*B_k-L_k(B_k)>=0",
        "```",
        "",
        "## Repaired k300 Stress",
        "",
        "```text",
        f"p-wall rows: {diagnostics['p_lower_rows']} lower and {diagnostics['p_upper_rows']} upper / {diagnostics['raw_log_total_rows']}",
        f"scaled-ceiling rows: {diagnostics['scaled_ceiling_rows']} / {diagnostics['adjacent_total_rows']}",
        f"log-lower rows: {diagnostics['log_lower_rows']} / {diagnostics['adjacent_total_rows']}",
        f"scaled-log-corridor rows: {diagnostics['scaled_log_corridor_rows']} / {diagnostics['adjacent_total_rows']}",
        f"ceiling-dominance rows: {diagnostics['ceiling_dominance_rows']} / {diagnostics['adjacent_total_rows']}",
        "```",
        "",
        "Sharpness:",
        "",
        "```text",
        f"min scaled ceiling margin: {diagnostics['min_scaled_ceiling_margin']['sample']} at lambda={diagnostics['min_scaled_ceiling_margin']['lam']}, k={diagnostics['min_scaled_ceiling_margin']['k']}",
        f"min raw upper margin: {diagnostics['min_raw_upper_margin']['sample']} at lambda={diagnostics['min_raw_upper_margin']['lam']}, k={diagnostics['min_raw_upper_margin']['k']}",
        f"min raw-upper minus scaled-ceiling slack: {sharpness['sample']} at lambda={sharpness['lam']}, k={sharpness['k']}",
        f"min scaled log width: {diagnostics['min_scaled_log_width']['sample']} at lambda={diagnostics['min_scaled_log_width']['lam']}, k={diagnostics['min_scaled_log_width']['k']}",
        "```",
        "",
        "Live recurrence target:",
        "",
        "```text",
        "log(1-4/(2*k+1)^2) <= delta_k",
        "delta_k <= h_(k+1)-((2*k+1)/(2*k+3))*h_k-2*p_k/(2*k+3)",
        "```",
        "",
        "Integration:",
        "",
        "```text",
        artifact["source_scaled_curvature_target"],
        artifact["source_raw_log_bridge"],
        artifact["source_linear_barrier_scout"],
        artifact["source_raw_corridor_target"],
        artifact["source_monotone_contraction_target"],
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
        "wrote Jensen-window PF negative-lambda scaled-curvature log-ceiling bridge: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
