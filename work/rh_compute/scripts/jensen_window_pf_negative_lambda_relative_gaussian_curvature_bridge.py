#!/usr/bin/env python3
"""Build the relative-Gaussian curvature bridge for the negative-lambda tail."""

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
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.md"


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
        raise RuntimeError("missing relative-Gaussian curvature extremum")
    value, lam, k = item
    return Extremum(decimal_format(value), lam, k)


def build_diagnostics(paths: list[Path], coefficient_cap: int = 300) -> dict:
    by_label_arb, by_label_sample = ratio_tables(paths, coefficient_cap)
    labels = sorted(by_label_sample, key=label_to_t)
    checked_ratio_max = coefficient_cap - 1
    adjacent_max = checked_ratio_max - 1

    b_positive_rows = 0
    b_decrease_rows = 0
    c_increase_rows = 0
    c_lambda_order_rows = 0

    min_b: tuple[Decimal, str, int] | None = None
    min_b_decrease: tuple[Decimal, str, int] | None = None
    min_c_increase: tuple[Decimal, str, int] | None = None
    min_c_lambda_gap: tuple[Decimal, str, int] | None = None
    min_c: tuple[Decimal, str, int] | None = None
    max_c: tuple[Decimal, str, int] | None = None

    b_arb_by_label: dict[str, dict[int, flint.arb]] = {}
    b_sample_by_label: dict[str, dict[int, Decimal]] = {}
    c_arb_by_label: dict[str, dict[int, flint.arb]] = {}
    c_sample_by_label: dict[str, dict[int, Decimal]] = {}

    for label in labels:
        b_arb_by_label[label] = {}
        b_sample_by_label[label] = {}
        c_arb_by_label[label] = {}
        c_sample_by_label[label] = {}
        for k in range(1, checked_ratio_max + 1):
            b_arb = h_arb(k) - by_label_arb[label][k].log()
            b_sample = h_decimal(k) - by_label_sample[label][k].ln()
            c_arb = flint.arb(2 * k + 1) * b_arb
            c_sample = Decimal(2 * k + 1) * b_sample
            b_arb_by_label[label][k] = b_arb
            b_sample_by_label[label][k] = b_sample
            c_arb_by_label[label][k] = c_arb
            c_sample_by_label[label][k] = c_sample
            if arb_positive(b_arb):
                b_positive_rows += 1
            min_b = update_min(min_b, b_sample, label, k)
            min_c = update_min(min_c, c_sample, label, k)
            max_c = update_max(max_c, c_sample, label, k)

        for k in range(1, adjacent_max + 1):
            b_drop_arb = b_arb_by_label[label][k] - b_arb_by_label[label][k + 1]
            b_drop_sample = b_sample_by_label[label][k] - b_sample_by_label[label][k + 1]
            c_gap_arb = c_arb_by_label[label][k + 1] - c_arb_by_label[label][k]
            c_gap_sample = c_sample_by_label[label][k + 1] - c_sample_by_label[label][k]
            if arb_positive(b_drop_arb):
                b_decrease_rows += 1
            if arb_positive(c_gap_arb):
                c_increase_rows += 1
            min_b_decrease = update_min(min_b_decrease, b_drop_sample, label, k)
            min_c_increase = update_min(min_c_increase, c_gap_sample, label, k)

    for left, right in zip(labels, labels[1:]):
        for k in range(1, checked_ratio_max + 1):
            gap_arb = c_arb_by_label[left][k] - c_arb_by_label[right][k]
            gap_sample = c_sample_by_label[left][k] - c_sample_by_label[right][k]
            if arb_positive(gap_arb):
                c_lambda_order_rows += 1
            min_c_lambda_gap = update_min(min_c_lambda_gap, gap_sample, f"{left} to {right}", k)

    return {
        "lambdas": labels,
        "coefficient_cap": coefficient_cap,
        "raw_log_total_rows": checked_ratio_max * len(labels),
        "adjacent_total_rows": adjacent_max * len(labels),
        "lambda_order_total_rows": max(0, len(labels) - 1) * checked_ratio_max,
        "b_positive_rows": b_positive_rows,
        "b_decrease_rows": b_decrease_rows,
        "c_increase_rows": c_increase_rows,
        "c_lambda_order_rows": c_lambda_order_rows,
        "min_b": asdict(extremum(min_b)),
        "min_b_decrease": asdict(extremum(min_b_decrease)),
        "min_c_increase": asdict(extremum(min_c_increase)),
        "min_c_lambda_gap": asdict(extremum(min_c_lambda_gap)),
        "min_c": asdict(extremum(min_c)),
        "max_c": asdict(extremum(max_c)),
    }


def build_artifact(paths: list[Path] | None = None) -> dict:
    if paths is None:
        paths = list(DEFAULT_ENCLOSURES)
    diagnostics = build_diagnostics(paths)
    rows = [
        {
            "id": "nlrgcb_01_relative_gaussian_log_coordinate",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "Let f_k be the log of the zeta moment sequence after subtracting any Gaussian raw-moment baseline; affine changes in f_k do not affect the second differences below.",
            "formula": "f_k=log(M_k/G_k), with G_(k+1)*G_(k-1)/G_k^2=(2*k+1)/(2*k-1)",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.md",
            ],
            "proof_boundary": "Exact coordinate choice only; no zeta inequality is proved.",
        },
        {
            "id": "nlrgcb_02_deficit_as_negative_second_difference",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The curvature deficit B_k is the negative second difference of the relative Gaussian log moment f_k.",
            "formula": "B_k=h_k-log(R_k)=2*f_k-f_(k-1)-f_(k+1)=-Delta^2 f_(k-1)",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.md",
            ],
            "proof_boundary": "Exact identity only; positivity of B_k remains a zeta-specific input.",
        },
        {
            "id": "nlrgcb_03_scaled_curvature_weighted_third_gap",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The scaled-curvature target C_(k+1)>=C_k is exactly a weighted four-point inequality for f_k.",
            "formula": "C_(k+1)-C_k=(2*k+1)*f_(k-1)-(6*k+5)*f_k+(6*k+7)*f_(k+1)-(2*k+3)*f_(k+2)",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.md",
            ],
            "proof_boundary": "Exact bridge only; the weighted four-point inequality is not proved.",
        },
        {
            "id": "nlrgcb_04_companion_upper_side_third_difference",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The companion monotone upper side B_(k+1)<=B_k is an ordinary third-difference sign condition in the same f_k coordinate.",
            "formula": "B_k-B_(k+1)=-f_(k-1)+3*f_k-3*f_(k+1)+f_(k+2)",
            "source_artifacts": [
                "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
                "outputs/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.md",
            ],
            "proof_boundary": "Exact translation only; it does not prove the companion monotone-contraction theorem.",
        },
        {
            "id": "nlrgcb_05_repaired_k300_curvature_ladder",
            "role": "finite_stress",
            "readiness": "not_ready_to_apply",
            "claim": "The repaired k300 data validate the finite curvature ladder B_k>0, B_(k+1)<=B_k, and C_(k+1)>=C_k on every checked row.",
            "diagnostics": diagnostics,
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
                "outputs/jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.md",
            ],
            "proof_boundary": "Finite repaired stress evidence only; not an all-k weighted curvature theorem.",
        },
        {
            "id": "nlrgcb_06_lambda_order_finite_pattern",
            "role": "finite_order_pattern",
            "readiness": "not_ready_to_apply",
            "claim": "On the checked grid, C_k is ordered downward as |lambda| grows across -25, -50, and -100.",
            "diagnostics": {
                "c_lambda_order_rows": diagnostics["c_lambda_order_rows"],
                "lambda_order_total_rows": diagnostics["lambda_order_total_rows"],
                "min_c_lambda_gap": diagnostics["min_c_lambda_gap"],
            },
            "proof_boundary": "Finite lambda-order diagnostic only; not a continuous-lambda theorem or tail theorem.",
        },
        {
            "id": "nlrgcb_07_live_tail_profile_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim": "A promoted tail proof may now target the relative-Gaussian weighted four-point inequality for f_k, with uniform remainders and a finite collar.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md",
                "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
                "outputs/jensen_window_pf_negative_lambda_high_order_taylor_scout.md",
            ],
            "gap": "Current signed/Taylor diagnostics do not yet provide a uniform all-k remainder theorem for the weighted four-point inequality.",
            "proof_boundary": "Live theorem-search route only.",
        },
        {
            "id": "nlrgcb_08_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "A promoted proof must state the relative Gaussian baseline, finite collar, lambda range, B-wall input, weighted C-monotonicity inequality, companion B-monotone side, and forbidden endpoint assumptions.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
                "outputs/signed_hankel_jensen_dependency_graph.md",
            ],
            "proof_boundary": "Proof-hygiene gate only; not a proof of scaled-curvature monotonicity, raw-corridor occupancy, cone entry, jwpf_06, RH, or Lambda <= 0.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "exact_available_rows": 4,
        "finite_stress_rows": 1,
        "finite_order_rows": 1,
        "live_routes": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "lambdas": diagnostics["lambdas"],
        "coefficient_cap": diagnostics["coefficient_cap"],
        "raw_log_total_rows": diagnostics["raw_log_total_rows"],
        "adjacent_total_rows": diagnostics["adjacent_total_rows"],
        "lambda_order_total_rows": diagnostics["lambda_order_total_rows"],
        "b_positive_rows": diagnostics["b_positive_rows"],
        "b_decrease_rows": diagnostics["b_decrease_rows"],
        "c_increase_rows": diagnostics["c_increase_rows"],
        "c_lambda_order_rows": diagnostics["c_lambda_order_rows"],
        "main_finding": (
            "The scaled-curvature target can be rewritten as a weighted four-point "
            "discrete-curvature inequality for the log moment sequence relative to the "
            "Gaussian baseline. Repaired k300 data validate B_k>0 on 897/897 rows, "
            "B_(k+1)<=B_k on 894/894 rows, C_(k+1)>=C_k on 894/894 rows, and "
            "C_k lambda-ordering on 598/598 rows."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge",
        "date": "2026-07-07",
        "status": "exact finite theorem-search diagnostic",
        "source_gaussian_curvature_matrix": "outputs/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.md",
        "source_scaled_curvature_target": "outputs/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.md",
        "source_scaled_log_ceiling_bridge": "outputs/jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.md",
        "source_monotone_contraction_target": "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "source_raw_corridor_target": "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.py",
        "proof_boundary": (
            "Exact finite theorem-search diagnostic only. It rewrites the scaled-curvature "
            "target in relative-Gaussian log-moment coordinates and checks repaired k300 "
            "finite stress, but it does not prove the all-k weighted four-point inequality, "
            "does not prove the companion monotone-contraction theorem, does not prove cone "
            "entry, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "The weighted four-point relative-Gaussian inequality remains an open zeta-specific all-k theorem.",
            "Repaired k300 evidence is finite stress evidence only.",
            "The companion monotone-contraction upper side remains a separate open theorem target.",
            "Generic moment positivity, endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][4]["diagnostics"]
    order = artifact["matrix_rows"][5]["diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian curvature bridge: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['b_positive_rows']} B-positive rows, "
        f"{summary['b_decrease_rows']} B-decrease rows, "
        f"{summary['c_increase_rows']} C-increase rows, "
        f"{summary['c_lambda_order_rows']} C-lambda-order rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Curvature Bridge",
        "",
        "Date: 2026-07-07",
        "",
        "Status: exact finite theorem-search diagnostic. This is not a proof",
        "of scaled-curvature monotonicity, the raw-corridor theorem, cone entry,",
        "RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge`.",
        "",
        "Proof boundary: this artifact rewrites the scaled-curvature target in",
        "relative-Gaussian log-moment coordinates and checks repaired k300",
        "evidence. It does not prove any all-`k` theorem.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Exact Relative-Gaussian Form",
        "",
        "Let `f_k` be the log moment sequence after subtracting any Gaussian",
        "raw-moment baseline with:",
        "",
        "```text",
        "G_(k+1)*G_(k-1)/G_k^2 = (2*k+1)/(2*k-1)",
        "```",
        "",
        "Then:",
        "",
        "```text",
        "B_k = h_k-log(R_k) = 2*f_k-f_(k-1)-f_(k+1)",
        "C_k = (2*k+1)*B_k",
        "```",
        "",
        "The scaled-curvature target is exactly:",
        "",
        "```text",
        "C_(k+1)-C_k",
        " = (2*k+1)*f_(k-1)-(6*k+5)*f_k+(6*k+7)*f_(k+1)-(2*k+3)*f_(k+2)",
        " >= 0",
        "```",
        "",
        "The companion monotone upper side is:",
        "",
        "```text",
        "B_k-B_(k+1) = -f_(k-1)+3*f_k-3*f_(k+1)+f_(k+2) >= 0",
        "```",
        "",
        "## Repaired k300 Stress",
        "",
        "```text",
        f"B-positive rows: {diagnostics['b_positive_rows']} / {diagnostics['raw_log_total_rows']}",
        f"B-decrease rows: {diagnostics['b_decrease_rows']} / {diagnostics['adjacent_total_rows']}",
        f"C-increase rows: {diagnostics['c_increase_rows']} / {diagnostics['adjacent_total_rows']}",
        f"C lambda-order rows: {order['c_lambda_order_rows']} / {order['lambda_order_total_rows']}",
        "```",
        "",
        "Sharp finite margins:",
        "",
        "```text",
        f"min B: {diagnostics['min_b']['sample']} at lambda={diagnostics['min_b']['lam']}, k={diagnostics['min_b']['k']}",
        f"min B decrease: {diagnostics['min_b_decrease']['sample']} at lambda={diagnostics['min_b_decrease']['lam']}, k={diagnostics['min_b_decrease']['k']}",
        f"min C increase: {diagnostics['min_c_increase']['sample']} at lambda={diagnostics['min_c_increase']['lam']}, k={diagnostics['min_c_increase']['k']}",
        f"min C lambda gap: {order['min_c_lambda_gap']['sample']} at lambda-pair={order['min_c_lambda_gap']['lam']}, k={order['min_c_lambda_gap']['k']}",
        f"C range: {diagnostics['min_c']['sample']} to {diagnostics['max_c']['sample']}",
        "```",
        "",
        "Integration:",
        "",
        "```text",
        artifact["source_gaussian_curvature_matrix"],
        artifact["source_scaled_curvature_target"],
        artifact["source_scaled_log_ceiling_bridge"],
        artifact["source_monotone_contraction_target"],
        artifact["source_raw_corridor_target"],
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
        "wrote Jensen-window PF negative-lambda relative-Gaussian curvature bridge: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
