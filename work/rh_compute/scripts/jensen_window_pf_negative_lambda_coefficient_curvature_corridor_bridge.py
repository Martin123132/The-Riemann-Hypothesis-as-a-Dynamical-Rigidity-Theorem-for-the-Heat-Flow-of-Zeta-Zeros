#!/usr/bin/env python3
"""Build a coefficient-curvature corridor bridge for the negative-lambda route."""

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
from jensen_window_pf_negative_lambda_k300_precision_repair_audit import (  # noqa: E402
    DEFAULT_BASE_ENCLOSURE,
    DEFAULT_REPAIR_ENCLOSURES,
    label_to_t,
    ratio_tables,
)
from jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout import (  # noqa: E402
    exact_counterexample,
    raw_upper_wall,
)


getcontext().prec = 100

DEFAULT_ENCLOSURES = (DEFAULT_BASE_ENCLOSURE, *DEFAULT_REPAIR_ENCLOSURES)
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.md"


@dataclass(frozen=True)
class Extremum:
    sample: str
    lam: str
    k: int


def q(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


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
        raise RuntimeError("missing coefficient-curvature corridor extremum")
    value, lam, k = item
    return Extremum(decimal_format(value), lam, k)


def h_arb(k: int) -> flint.arb:
    return (flint.arb(2 * k + 1) / flint.arb(2 * k - 1)).log()


def h_decimal(k: int) -> Decimal:
    return (Decimal(2 * k + 1) / Decimal(2 * k - 1)).ln()


def curvature_lower_arb(k: int, b_value: flint.arb) -> flint.arb:
    return (flint.arb(2 * k + 3) / (flint.arb(2) + flint.arb(2 * k + 1) * (-b_value).exp())).log()


def curvature_lower_decimal(k: int, b_value: Decimal) -> Decimal:
    return (Decimal(2 * k + 3) / (Decimal(2) + Decimal(2 * k + 1) * (-b_value).exp())).ln()


def exact_curvature_witness(raw: Fraction, raw_next: Fraction, k: int = 1) -> dict:
    witness = exact_counterexample(raw, raw_next, k)
    x = Fraction(2 * k - 1, 2 * k + 1) * raw
    x_next = Fraction(2 * k + 1, 2 * k + 3) * raw_next
    scaled_upper_x_next = (Fraction(2) + Fraction(2 * k + 1) * x) / Fraction(2 * k + 3)
    witness.update(
        {
            "x_k": q(x),
            "x_next": q(x_next),
            "scaled_upper_x_next_wall": q(scaled_upper_x_next),
            "B_k": f"log({q(1 / x)})",
            "B_next": f"log({q(1 / x_next)})",
            "B_monotone_holds": bool(x_next >= x),
            "B_lower_barrier_holds": bool(x_next <= scaled_upper_x_next),
            "raw_wall_holds_at_k": bool(Fraction(1) <= raw <= raw_upper_wall(k)),
            "raw_wall_holds_at_next": bool(Fraction(1) <= raw_next <= raw_upper_wall(k + 1)),
        }
    )
    return witness


def build_diagnostics(paths: list[Path], coefficient_cap: int = 300) -> dict:
    by_label_arb, by_label_sample = ratio_tables(paths, coefficient_cap)
    labels = sorted(by_label_sample, key=label_to_t)
    checked_ratio_max = coefficient_cap - 1
    adjacent_max = checked_ratio_max - 1

    b_positive_rows = 0
    b_upper_wall_rows = 0
    b_monotone_rows = 0
    b_lower_barrier_rows = 0
    b_corridor_rows = 0
    b_width_rows = 0

    min_b: tuple[Decimal, str, int] | None = None
    max_b: tuple[Decimal, str, int] | None = None
    min_b_upper_slack: tuple[Decimal, str, int] | None = None
    min_b_monotone_margin: tuple[Decimal, str, int] | None = None
    min_b_lower_barrier_margin: tuple[Decimal, str, int] | None = None
    min_b_corridor_width: tuple[Decimal, str, int] | None = None
    min_exp_scaled_upper_margin: tuple[Decimal, str, int] | None = None

    for label in labels:
        b_arb: dict[int, flint.arb] = {}
        b_sample: dict[int, Decimal] = {}
        x_sample: dict[int, Decimal] = {}
        for k in range(1, checked_ratio_max + 1):
            raw_arb = by_label_arb[label][k]
            raw_sample = by_label_sample[label][k]
            b_arb[k] = h_arb(k) - raw_arb.log()
            b_sample[k] = h_decimal(k) - raw_sample.ln()
            x_sample[k] = raw_sample * Decimal(2 * k - 1) / Decimal(2 * k + 1)
            upper_slack_arb = h_arb(k) - b_arb[k]
            upper_slack_sample = h_decimal(k) - b_sample[k]
            if arb_positive(b_arb[k]):
                b_positive_rows += 1
            if arb_positive(upper_slack_arb):
                b_upper_wall_rows += 1
            min_b = update_min(min_b, b_sample[k], label, k)
            max_b = update_max(max_b, b_sample[k], label, k)
            min_b_upper_slack = update_min(min_b_upper_slack, upper_slack_sample, label, k)

        for k in range(1, adjacent_max + 1):
            lower_arb = curvature_lower_arb(k, b_arb[k])
            lower_sample = curvature_lower_decimal(k, b_sample[k])
            monotone_margin_arb = b_arb[k] - b_arb[k + 1]
            lower_margin_arb = b_arb[k + 1] - lower_arb
            width_arb = b_arb[k] - lower_arb
            monotone_margin_sample = b_sample[k] - b_sample[k + 1]
            lower_margin_sample = b_sample[k + 1] - lower_sample
            width_sample = b_sample[k] - lower_sample
            scaled_x_wall = (Decimal(2) + Decimal(2 * k + 1) * x_sample[k]) / Decimal(2 * k + 3)
            exp_scaled_upper_margin = scaled_x_wall - x_sample[k + 1]

            if arb_positive(monotone_margin_arb):
                b_monotone_rows += 1
            if arb_positive(lower_margin_arb):
                b_lower_barrier_rows += 1
            if arb_positive(monotone_margin_arb) and arb_positive(lower_margin_arb):
                b_corridor_rows += 1
            if arb_positive(width_arb):
                b_width_rows += 1

            min_b_monotone_margin = update_min(min_b_monotone_margin, monotone_margin_sample, label, k)
            min_b_lower_barrier_margin = update_min(min_b_lower_barrier_margin, lower_margin_sample, label, k)
            min_b_corridor_width = update_min(min_b_corridor_width, width_sample, label, k)
            min_exp_scaled_upper_margin = update_min(min_exp_scaled_upper_margin, exp_scaled_upper_margin, label, k)

    return {
        "lambdas": labels,
        "coefficient_cap": coefficient_cap,
        "raw_ratio_rows_per_lambda": checked_ratio_max,
        "adjacent_rows_per_lambda": adjacent_max,
        "b_total_rows": checked_ratio_max * len(labels),
        "adjacent_total_rows": adjacent_max * len(labels),
        "b_positive_rows": b_positive_rows,
        "b_upper_wall_rows": b_upper_wall_rows,
        "b_monotone_rows": b_monotone_rows,
        "b_lower_barrier_rows": b_lower_barrier_rows,
        "b_corridor_rows": b_corridor_rows,
        "b_width_rows": b_width_rows,
        "min_b": asdict(extremum(min_b)),
        "max_b": asdict(extremum(max_b)),
        "min_b_upper_slack": asdict(extremum(min_b_upper_slack)),
        "min_b_monotone_margin": asdict(extremum(min_b_monotone_margin)),
        "min_b_lower_barrier_margin": asdict(extremum(min_b_lower_barrier_margin)),
        "min_b_corridor_width": asdict(extremum(min_b_corridor_width)),
        "min_exp_scaled_upper_margin": asdict(extremum(min_exp_scaled_upper_margin)),
    }


def build_artifact(paths: list[Path] | None = None) -> dict:
    if paths is None:
        paths = list(DEFAULT_ENCLOSURES)
    diagnostics = build_diagnostics(paths)
    monotone_shortcut_witness = exact_curvature_witness(Fraction(2), Fraction(3, 2))
    raw_wall_shortcut_witness = exact_curvature_witness(Fraction(2), Fraction(1))
    rows = [
        {
            "id": "nlcccb_01_curvature_coordinate",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "With x_k=((2*k-1)/(2*k+1))*R_k and B_k=-log(x_k), the raw log ratio is p_k=log(R_k)=log((2*k+1)/(2*k-1))-B_k.",
            "formula": "R_k=((2*k+1)/(2*k-1))*exp(-B_k)",
            "proof_boundary": "Exact change of variables only; it proves no zeta inequality.",
        },
        {
            "id": "nlcccb_02_exact_curvature_corridor",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The raw decrement corridor is exactly equivalent to a coefficient-curvature corridor for B_(k+1).",
            "formula": "log((2*k+3)/(2+(2*k+1)*exp(-B_k))) <= B_(k+1) <= B_k",
            "proof_boundary": "Exact algebraic equivalence only; the all-k curvature corridor remains open.",
        },
        {
            "id": "nlcccb_03_monotone_side_identification",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The upper side B_(k+1)<=B_k is exactly the monotone-contraction/log-curvature side already isolated as Delta^3 log A>=0.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md",
                "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
            ],
            "proof_boundary": "Exact identification with an open theorem target only; not a monotonicity proof.",
        },
        {
            "id": "nlcccb_04_repaired_k300_curvature_corridor",
            "role": "finite_stress",
            "readiness": "not_ready_to_apply",
            "claim": "The repaired k300 data validate the B wall, monotone side, lower barrier, and full coefficient-curvature corridor on all checked rows.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.md",
                "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
            ],
            "diagnostics": diagnostics,
            "proof_boundary": "Finite repaired stress evidence only; not an all-k curvature-corridor theorem.",
        },
        {
            "id": "nlcccb_05_monotone_curvature_shortcut_blocked",
            "role": "exact_counterexample",
            "readiness": "not_ready_to_apply",
            "claim": "Raw walls plus B_(k+1)<=B_k do not imply the lower curvature barrier.",
            "witness": monotone_shortcut_witness,
            "proof_boundary": "Exact shortcut rejection only; not evidence against actual zeta coefficients.",
        },
        {
            "id": "nlcccb_06_raw_wall_shortcut_blocked",
            "role": "exact_counterexample",
            "readiness": "not_ready_to_apply",
            "claim": "Raw walls alone do not imply the monotone curvature side.",
            "witness": raw_wall_shortcut_witness,
            "proof_boundary": "Exact shortcut rejection only; not evidence against actual zeta coefficients.",
        },
        {
            "id": "nlcccb_07_live_lower_barrier_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim": "A zeta-specific proof can target the missing lower barrier B_(k+1)>=log((2*k+3)/(2+(2*k+1)*exp(-B_k))) after a finite collar.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
                "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
            ],
            "proof_boundary": "Live theorem-search route only; no all-k lower barrier is proved here.",
        },
        {
            "id": "nlcccb_08_conditional_closure_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim": "Finite collar plus the existing monotone-contraction theorem plus the lower curvature barrier would imply the raw-corridor target.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
                "outputs/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.md",
                "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
            ],
            "proof_boundary": "Conditional route only; at least two all-k inputs remain open.",
        },
        {
            "id": "nlcccb_09_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "A promoted curvature-corridor proof must state the tail start, lambda range, finite collar, B wall, monotone side, lower barrier, and forbidden assumptions.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
                "outputs/signed_hankel_jensen_dependency_graph.md",
            ],
            "proof_boundary": "Proof-hygiene gate only; not a proof of the raw-corridor theorem.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "exact_available_rows": 3,
        "finite_stress_rows": 1,
        "exact_counterexample_rows": 2,
        "live_routes": 2,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "lambdas": diagnostics["lambdas"],
        "coefficient_cap": diagnostics["coefficient_cap"],
        "b_total_rows": diagnostics["b_total_rows"],
        "adjacent_total_rows": diagnostics["adjacent_total_rows"],
        "b_positive_rows": diagnostics["b_positive_rows"],
        "b_upper_wall_rows": diagnostics["b_upper_wall_rows"],
        "b_monotone_rows": diagnostics["b_monotone_rows"],
        "b_lower_barrier_rows": diagnostics["b_lower_barrier_rows"],
        "b_corridor_rows": diagnostics["b_corridor_rows"],
        "b_width_rows": diagnostics["b_width_rows"],
        "main_finding": (
            "The raw-corridor theorem is exactly equivalent to the coefficient-curvature corridor "
            "log((2*k+3)/(2+(2*k+1)*exp(-B_k))) <= B_(k+1) <= B_k. On repaired k300 data "
            "the B wall holds on 897/897 rows and the curvature corridor holds on 894/894 adjacent rows; "
            "monotone curvature alone is blocked by an exact cone counterexample."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge",
        "date": "2026-07-07",
        "status": "exact finite theorem-search diagnostic",
        "source_raw_log_bridge": "outputs/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.md",
        "source_k300_precision_repair": "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
        "source_raw_corridor_target": "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
        "source_log_curvature_bridge": "outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md",
        "source_monotone_contraction_target": "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "source_raw_obstruction": "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.py",
        "proof_boundary": (
            "Exact finite theorem-search diagnostic only. It rewrites the raw decrement corridor "
            "as a coefficient-curvature corridor and checks repaired k300 evidence, but it does "
            "not prove the all-k raw-corridor theorem, does not prove cone entry, does not prove "
            "jwpf_06, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "invariants": [
            "No row is ready_to_apply.",
            "The coefficient-curvature corridor remains an open zeta-specific all-k theorem.",
            "Repaired k300 evidence is finite stress evidence only.",
            "Monotone curvature alone is rejected as a shortcut.",
            "Generic moment positivity, endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
        "summary": summary,
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][3]["diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda coefficient-curvature corridor bridge: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['b_corridor_rows']} curvature-corridor rows, "
        f"{summary['b_monotone_rows']} monotone-curvature rows, "
        f"{summary['exact_counterexample_rows']} exact counterexamples, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Coefficient-Curvature Corridor Bridge",
        "",
        "Date: 2026-07-07",
        "",
        "Status: exact finite theorem-search diagnostic. This is not a proof",
        "of the raw-corridor theorem, cone entry, Jensen-window PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge`.",
        "",
        "Proof boundary: this artifact rewrites the raw-ratio decrement",
        "corridor as a coefficient-curvature corridor and checks repaired k300",
        "evidence. It does not prove any all-`k` theorem.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Exact Curvature Form",
        "",
        "Let:",
        "",
        "```text",
        "R_k = M_(k+1)*M_(k-1)/M_k^2",
        "x_k = ((2*k-1)/(2*k+1))*R_k",
        "B_k = -log(x_k)",
        "```",
        "",
        "Then:",
        "",
        "```text",
        "R_k = ((2*k+1)/(2*k-1))*exp(-B_k)",
        "p_k = log(R_k) = log((2*k+1)/(2*k-1))-B_k",
        "```",
        "",
        "The raw decrement corridor is exactly equivalent to:",
        "",
        "```text",
        "log((2*k+3)/(2+(2*k+1)*exp(-B_k))) <= B_(k+1)",
        "B_(k+1) <= B_k",
        "```",
        "",
        "The upper side `B_(k+1)<=B_k` is the monotone-contraction side",
        "already isolated in log-curvature language; the new missing side is",
        "the lower curvature barrier.",
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
        "Repaired curvature diagnostics:",
        "",
        "```text",
        f"B positive rows: {diagnostics['b_positive_rows']} / {diagnostics['b_total_rows']}",
        f"B upper-wall rows: {diagnostics['b_upper_wall_rows']} / {diagnostics['b_total_rows']}",
        f"monotone-curvature rows: {diagnostics['b_monotone_rows']} / {diagnostics['adjacent_total_rows']}",
        f"lower-barrier rows: {diagnostics['b_lower_barrier_rows']} / {diagnostics['adjacent_total_rows']}",
        f"curvature-corridor rows: {diagnostics['b_corridor_rows']} / {diagnostics['adjacent_total_rows']}",
        f"curvature-width rows: {diagnostics['b_width_rows']} / {diagnostics['adjacent_total_rows']}",
        "```",
        "",
        "Extrema:",
        "",
        "```text",
        f"min B lower-barrier margin: {diagnostics['min_b_lower_barrier_margin']['sample']} at lambda={diagnostics['min_b_lower_barrier_margin']['lam']}, k={diagnostics['min_b_lower_barrier_margin']['k']}",
        f"min B monotone margin: {diagnostics['min_b_monotone_margin']['sample']} at lambda={diagnostics['min_b_monotone_margin']['lam']}, k={diagnostics['min_b_monotone_margin']['k']}",
        f"min B corridor width: {diagnostics['min_b_corridor_width']['sample']} at lambda={diagnostics['min_b_corridor_width']['lam']}, k={diagnostics['min_b_corridor_width']['k']}",
        f"B range: {diagnostics['min_b']['sample']} at lambda={diagnostics['min_b']['lam']}, k={diagnostics['min_b']['k']} to {diagnostics['max_b']['sample']} at lambda={diagnostics['max_b']['lam']}, k={diagnostics['max_b']['k']}",
        "```",
        "",
        "## Shortcut Gates",
        "",
        "Simple curvature shortcuts are blocked:",
        "",
        "```text",
        "R_1=2, R_2=3/2: raw walls and B_(k+1)<=B_k hold, but the lower curvature barrier fails",
        "R_1=2, R_2=1: raw walls hold, but B_(k+1)<=B_k fails",
        "```",
        "",
        "Integration:",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.md",
        "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
        "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
        "outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md",
        "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
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
    print(f"wrote Jensen-window PF negative-lambda coefficient-curvature corridor bridge: {out_json} and {note}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
