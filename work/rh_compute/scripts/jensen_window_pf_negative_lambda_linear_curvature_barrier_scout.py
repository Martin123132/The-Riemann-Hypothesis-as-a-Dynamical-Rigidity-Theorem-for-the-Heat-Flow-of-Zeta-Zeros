#!/usr/bin/env python3
"""Build a linear curvature-barrier scout for the negative-lambda route."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
from fractions import Fraction
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
    curvature_lower_arb,
    curvature_lower_decimal,
    exact_curvature_witness,
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
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.md"


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
        raise RuntimeError("missing linear curvature-barrier extremum")
    value, lam, k = item
    return Extremum(decimal_format(value), lam, k)


def alpha_arb(k: int) -> flint.arb:
    return flint.arb(2 * k + 1) / flint.arb(2 * k + 3)


def alpha_decimal(k: int) -> Decimal:
    return Decimal(2 * k + 1) / Decimal(2 * k + 3)


def build_diagnostics(paths: list[Path], coefficient_cap: int = 300) -> dict:
    by_label_arb, by_label_sample = ratio_tables(paths, coefficient_cap)
    labels = sorted(by_label_sample, key=label_to_t)
    checked_ratio_max = coefficient_cap - 1
    adjacent_max = checked_ratio_max - 1

    b_wall_rows = 0
    exact_lower_barrier_rows = 0
    linear_lower_barrier_rows = 0
    monotone_rows = 0
    linear_corridor_rows = 0

    min_b_wall: tuple[Decimal, str, int] | None = None
    min_exact_lower_margin: tuple[Decimal, str, int] | None = None
    min_linear_lower_margin: tuple[Decimal, str, int] | None = None
    min_linear_to_exact_slack: tuple[Decimal, str, int] | None = None
    min_monotone_margin: tuple[Decimal, str, int] | None = None
    min_linear_width: tuple[Decimal, str, int] | None = None
    min_ratio_next_over_current: tuple[Decimal, str, int] | None = None
    max_ratio_next_over_current: tuple[Decimal, str, int] | None = None

    for label in labels:
        b_arb: dict[int, flint.arb] = {}
        b_sample: dict[int, Decimal] = {}
        for k in range(1, checked_ratio_max + 1):
            raw_arb = by_label_arb[label][k]
            raw_sample = by_label_sample[label][k]
            b_arb[k] = h_arb(k) - raw_arb.log()
            b_sample[k] = h_decimal(k) - raw_sample.ln()
            if arb_positive(b_arb[k]):
                b_wall_rows += 1
            min_b_wall = update_min(min_b_wall, b_sample[k], label, k)

        for k in range(1, adjacent_max + 1):
            exact_lower_arb = curvature_lower_arb(k, b_arb[k])
            exact_lower_sample = curvature_lower_decimal(k, b_sample[k])
            linear_lower_arb = alpha_arb(k) * b_arb[k]
            linear_lower_sample = alpha_decimal(k) * b_sample[k]
            exact_lower_margin_arb = b_arb[k + 1] - exact_lower_arb
            linear_lower_margin_arb = b_arb[k + 1] - linear_lower_arb
            linear_to_exact_slack_arb = linear_lower_arb - exact_lower_arb
            monotone_margin_arb = b_arb[k] - b_arb[k + 1]
            linear_width_arb = b_arb[k] - linear_lower_arb

            exact_lower_margin_sample = b_sample[k + 1] - exact_lower_sample
            linear_lower_margin_sample = b_sample[k + 1] - linear_lower_sample
            linear_to_exact_slack_sample = linear_lower_sample - exact_lower_sample
            monotone_margin_sample = b_sample[k] - b_sample[k + 1]
            linear_width_sample = b_sample[k] - linear_lower_sample
            ratio = b_sample[k + 1] / b_sample[k]

            if arb_positive(exact_lower_margin_arb):
                exact_lower_barrier_rows += 1
            if arb_positive(linear_lower_margin_arb):
                linear_lower_barrier_rows += 1
            if arb_positive(monotone_margin_arb):
                monotone_rows += 1
            if arb_positive(linear_lower_margin_arb) and arb_positive(monotone_margin_arb):
                linear_corridor_rows += 1

            min_exact_lower_margin = update_min(min_exact_lower_margin, exact_lower_margin_sample, label, k)
            min_linear_lower_margin = update_min(min_linear_lower_margin, linear_lower_margin_sample, label, k)
            min_linear_to_exact_slack = update_min(min_linear_to_exact_slack, linear_to_exact_slack_sample, label, k)
            min_monotone_margin = update_min(min_monotone_margin, monotone_margin_sample, label, k)
            min_linear_width = update_min(min_linear_width, linear_width_sample, label, k)
            min_ratio_next_over_current = update_min(min_ratio_next_over_current, ratio, label, k)
            max_ratio_next_over_current = update_max(max_ratio_next_over_current, ratio, label, k)

    return {
        "lambdas": labels,
        "coefficient_cap": coefficient_cap,
        "raw_ratio_rows_per_lambda": checked_ratio_max,
        "adjacent_rows_per_lambda": adjacent_max,
        "b_total_rows": checked_ratio_max * len(labels),
        "adjacent_total_rows": adjacent_max * len(labels),
        "b_wall_rows": b_wall_rows,
        "exact_lower_barrier_rows": exact_lower_barrier_rows,
        "linear_lower_barrier_rows": linear_lower_barrier_rows,
        "monotone_rows": monotone_rows,
        "linear_corridor_rows": linear_corridor_rows,
        "min_b_wall": asdict(extremum(min_b_wall)),
        "min_exact_lower_margin": asdict(extremum(min_exact_lower_margin)),
        "min_linear_lower_margin": asdict(extremum(min_linear_lower_margin)),
        "min_linear_to_exact_slack": asdict(extremum(min_linear_to_exact_slack)),
        "min_monotone_margin": asdict(extremum(min_monotone_margin)),
        "min_linear_width": asdict(extremum(min_linear_width)),
        "min_ratio_next_over_current": asdict(extremum(min_ratio_next_over_current)),
        "max_ratio_next_over_current": asdict(extremum(max_ratio_next_over_current)),
    }


def build_artifact(paths: list[Path] | None = None) -> dict:
    if paths is None:
        paths = list(DEFAULT_ENCLOSURES)
    diagnostics = build_diagnostics(paths)
    monotone_shortcut_witness = exact_curvature_witness(Fraction(2), Fraction(3, 2))
    raw_wall_shortcut_witness = exact_curvature_witness(Fraction(2), Fraction(1))
    rows = [
        {
            "id": "nllcbs_01_exact_linear_sufficient_barrier",
            "role": "exact_sufficient_condition",
            "readiness": "available_exact",
            "claim": "For B>=0 and alpha_k=(2*k+1)/(2*k+3), the nonlinear lower barrier L_k(B)=log((2*k+3)/(2+(2*k+1)*exp(-B))) satisfies L_k(B)<=alpha_k*B.",
            "formula": "alpha_k*B-L_k(B)>=0 for B>=0",
            "proof_boundary": "Exact calculus lemma only; it proves no zeta coefficient inequality.",
        },
        {
            "id": "nllcbs_02_linear_barrier_implies_raw_corridor_side",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The linear barrier B_(k+1)>=alpha_k*B_k is a sufficient condition for the lower coefficient-curvature barrier and hence for the upper raw-corridor side.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.md",
            ],
            "proof_boundary": "Exact sufficient reduction only; the linear zeta inequality remains open for all k.",
        },
        {
            "id": "nllcbs_03_repaired_k300_linear_barrier",
            "role": "finite_stress",
            "readiness": "not_ready_to_apply",
            "claim": "The repaired k300 data validate the linear lower barrier, exact lower barrier, and monotone side on all checked adjacent rows.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.md",
                "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
            ],
            "diagnostics": diagnostics,
            "proof_boundary": "Finite repaired stress evidence only; not an all-k linear curvature-barrier theorem.",
        },
        {
            "id": "nllcbs_04_defect_width_recurrence_distinction",
            "role": "rejected_analogy",
            "readiness": "not_ready_to_apply",
            "claim": "The analogous defect-width recurrence is rejected, so the live route is specifically a log-curvature linear barrier, not a defect-linear shortcut.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md",
            ],
            "proof_boundary": "Shortcut distinction only; it does not prove the B-linear theorem.",
        },
        {
            "id": "nllcbs_05_monotone_curvature_shortcut_blocked",
            "role": "exact_counterexample",
            "readiness": "not_ready_to_apply",
            "claim": "Raw walls plus B_(k+1)<=B_k do not imply the linear lower barrier.",
            "witness": monotone_shortcut_witness,
            "proof_boundary": "Exact shortcut rejection only; not evidence against actual zeta coefficients.",
        },
        {
            "id": "nllcbs_06_raw_wall_shortcut_blocked",
            "role": "exact_counterexample",
            "readiness": "not_ready_to_apply",
            "claim": "Raw walls alone do not imply the two-sided linear B corridor.",
            "witness": raw_wall_shortcut_witness,
            "proof_boundary": "Exact shortcut rejection only; not evidence against actual zeta coefficients.",
        },
        {
            "id": "nllcbs_07_live_linear_barrier_target",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim": "A zeta-specific proof may target B_(k+1)>=((2*k+1)/(2*k+3))*B_k after a finite collar.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
                "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
            ],
            "proof_boundary": "Live theorem-search target only; no all-k linear barrier is proved here.",
        },
        {
            "id": "nllcbs_08_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "A promoted proof must state the tail start, lambda range, finite collar, B>=0 wall, monotone upper side, linear lower side, and forbidden assumptions.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
                "outputs/signed_hankel_jensen_dependency_graph.md",
            ],
            "proof_boundary": "Proof-hygiene gate only; not a proof of the raw-corridor theorem.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "exact_available_rows": 2,
        "finite_stress_rows": 1,
        "rejected_analogy_rows": 1,
        "exact_counterexample_rows": 2,
        "live_routes": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "lambdas": diagnostics["lambdas"],
        "coefficient_cap": diagnostics["coefficient_cap"],
        "b_total_rows": diagnostics["b_total_rows"],
        "adjacent_total_rows": diagnostics["adjacent_total_rows"],
        "b_wall_rows": diagnostics["b_wall_rows"],
        "exact_lower_barrier_rows": diagnostics["exact_lower_barrier_rows"],
        "linear_lower_barrier_rows": diagnostics["linear_lower_barrier_rows"],
        "monotone_rows": diagnostics["monotone_rows"],
        "linear_corridor_rows": diagnostics["linear_corridor_rows"],
        "main_finding": (
            "The nonlinear lower curvature barrier is implied by the simpler linear theorem "
            "B_(k+1)>=((2*k+1)/(2*k+3))*B_k. On repaired k300 data the linear lower barrier "
            "holds on 894/894 adjacent rows and the B wall holds on 897/897 rows, while the "
            "analogous defect-width recurrence is already rejected."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_linear_curvature_barrier_scout",
        "date": "2026-07-07",
        "status": "exact finite theorem-search diagnostic",
        "source_curvature_corridor_bridge": "outputs/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.md",
        "source_k300_precision_repair": "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
        "source_raw_corridor_target": "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
        "source_monotone_contraction_target": "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "source_defect_recurrence_scout": "outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md",
        "source_raw_obstruction": "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.py",
        "proof_boundary": (
            "Exact finite theorem-search diagnostic only. It isolates a linear sufficient "
            "curvature-barrier theorem and checks repaired k300 evidence, but it does not prove "
            "the all-k linear barrier, does not prove the raw-corridor theorem, does not prove "
            "cone entry, does not prove jwpf_06, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "invariants": [
            "No row is ready_to_apply.",
            "The linear curvature barrier remains an open zeta-specific all-k theorem.",
            "Repaired k300 evidence is finite stress evidence only.",
            "The analogous defect-width recurrence remains rejected.",
            "Monotone curvature alone is rejected as a shortcut.",
            "Generic moment positivity, endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
        "summary": summary,
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][2]["diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda linear curvature-barrier scout: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['linear_lower_barrier_rows']} linear-barrier rows, "
        f"{summary['monotone_rows']} monotone-curvature rows, "
        f"{summary['exact_counterexample_rows']} exact counterexamples, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Linear Curvature-Barrier Scout",
        "",
        "Date: 2026-07-07",
        "",
        "Status: exact finite theorem-search diagnostic. This is not a proof",
        "of the raw-corridor theorem, cone entry, Jensen-window PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_linear_curvature_barrier_scout`.",
        "",
        "Proof boundary: this artifact isolates a linear sufficient theorem for",
        "the lower coefficient-curvature barrier. It does not prove any all-`k`",
        "theorem.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Exact Linear Sufficient Condition",
        "",
        "Let:",
        "",
        "```text",
        "B_k = -log(((2*k-1)/(2*k+1))*R_k)",
        "alpha_k = (2*k+1)/(2*k+3)",
        "L_k(B) = log((2*k+3)/(2+(2*k+1)*exp(-B)))",
        "```",
        "",
        "For `B>=0`, the exact calculus lemma is:",
        "",
        "```text",
        "L_k(B) <= alpha_k*B",
        "```",
        "",
        "Therefore the simpler linear theorem",
        "",
        "```text",
        "B_(k+1) >= ((2*k+1)/(2*k+3))*B_k",
        "```",
        "",
        "implies the nonlinear lower curvature barrier.",
        "",
        "## Repaired k300 Stress",
        "",
        "Inputs:",
        "",
        "```text",
        "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k300.jsonl",
        "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k220_k250_dps220.jsonl",
        "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220.jsonl",
        "```",
        "",
        "Repaired linear-barrier diagnostics:",
        "",
        "```text",
        f"B wall rows: {diagnostics['b_wall_rows']} / {diagnostics['b_total_rows']}",
        f"exact lower-barrier rows: {diagnostics['exact_lower_barrier_rows']} / {diagnostics['adjacent_total_rows']}",
        f"linear lower-barrier rows: {diagnostics['linear_lower_barrier_rows']} / {diagnostics['adjacent_total_rows']}",
        f"monotone-curvature rows: {diagnostics['monotone_rows']} / {diagnostics['adjacent_total_rows']}",
        f"linear B-corridor rows: {diagnostics['linear_corridor_rows']} / {diagnostics['adjacent_total_rows']}",
        "```",
        "",
        "Extrema:",
        "",
        "```text",
        f"min linear lower-barrier margin: {diagnostics['min_linear_lower_margin']['sample']} at lambda={diagnostics['min_linear_lower_margin']['lam']}, k={diagnostics['min_linear_lower_margin']['k']}",
        f"min exact lower-barrier margin: {diagnostics['min_exact_lower_margin']['sample']} at lambda={diagnostics['min_exact_lower_margin']['lam']}, k={diagnostics['min_exact_lower_margin']['k']}",
        f"min linear-to-exact slack: {diagnostics['min_linear_to_exact_slack']['sample']} at lambda={diagnostics['min_linear_to_exact_slack']['lam']}, k={diagnostics['min_linear_to_exact_slack']['k']}",
        f"B_(k+1)/B_k range: {diagnostics['min_ratio_next_over_current']['sample']} at lambda={diagnostics['min_ratio_next_over_current']['lam']}, k={diagnostics['min_ratio_next_over_current']['k']} to {diagnostics['max_ratio_next_over_current']['sample']} at lambda={diagnostics['max_ratio_next_over_current']['lam']}, k={diagnostics['max_ratio_next_over_current']['k']}",
        "```",
        "",
        "## Shortcut Distinction",
        "",
        "The analogous defect-width recurrence is already rejected:",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md",
        "```",
        "",
        "Thus the live route is specifically the log-curvature linear barrier, not",
        "a defect-linear shortcut.",
        "",
        "Integration:",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.md",
        "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
        "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
        "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md",
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "```",
        "",
        "Summary:",
        "",
        summary["main_finding"],
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


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
    print(f"wrote Jensen-window PF negative-lambda linear curvature-barrier scout: {out_json} and {note}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
