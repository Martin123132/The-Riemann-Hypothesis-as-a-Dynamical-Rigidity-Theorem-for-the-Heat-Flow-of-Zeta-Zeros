#!/usr/bin/env python3
"""Build the negative-lambda raw-moment bridge matrix for the adaptive route."""

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
    contraction,
    decimal_format,
    load_enclosures,
)


getcontext().prec = 100

DEFAULT_ENCLOSURE_JSONL = (
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k200.jsonl",
)
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.md"


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
        raise RuntimeError("missing extremum")
    value, lam, k = item
    return Extremum(decimal_format(value), lam, k)


def label_to_t(label: str) -> Decimal:
    return abs(Decimal(label))


def raw_ratio_factor_arb(k: int) -> flint.arb:
    return flint.arb(2 * k + 1) / flint.arb(2 * k - 1)


def raw_ratio_factor_decimal(k: int) -> Decimal:
    return Decimal(2 * k + 1) / Decimal(2 * k - 1)


def raw_half_width_threshold_arb(k: int) -> flint.arb:
    return flint.arb(2 * k) / flint.arb(2 * k - 1)


def raw_half_width_threshold_decimal(k: int) -> Decimal:
    return Decimal(2 * k) / Decimal(2 * k - 1)


def raw_one_third_threshold_arb(k: int) -> flint.arb:
    return flint.arb(6 * k + 1) / flint.arb(3 * (2 * k - 1))


def raw_one_third_threshold_decimal(k: int) -> Decimal:
    return Decimal(6 * k + 1) / Decimal(3 * (2 * k - 1))


def bridge_lower_factor_arb(k: int) -> flint.arb:
    return flint.arb((2 * k - 1) * (2 * k + 3)) / flint.arb((2 * k + 1) * (2 * k + 1))


def bridge_lower_factor_decimal(k: int) -> Decimal:
    return Decimal((2 * k - 1) * (2 * k + 3)) / Decimal((2 * k + 1) * (2 * k + 1))


def build_raw_ratios(paths: list[Path]) -> tuple[dict[str, dict[int, flint.arb]], dict[str, dict[int, Decimal]], int]:
    balls, samples, labels = load_enclosures(paths)
    if not labels:
        raise RuntimeError("no coefficient enclosures found")
    coefficient_k_max = max(index for (_lam, index) in balls)
    checked_ratio_max = coefficient_k_max - 1
    by_label_arb: dict[str, dict[int, flint.arb]] = {}
    by_label_sample: dict[str, dict[int, Decimal]] = {}

    for lam in sorted(labels, key=lambda value: label_to_t(labels[value])):
        label = labels[lam]
        by_label_arb[label] = {}
        by_label_sample[label] = {}
        for k in range(1, checked_ratio_max + 1):
            x_arb = contraction({idx: balls[(lam, idx)] for idx in (k - 1, k, k + 1)}, k)
            x_sample = contraction({idx: samples[(lam, idx)] for idx in (k - 1, k, k + 1)}, k)
            by_label_arb[label][k] = x_arb * raw_ratio_factor_arb(k)
            by_label_sample[label][k] = x_sample * raw_ratio_factor_decimal(k)
    return by_label_arb, by_label_sample, checked_ratio_max


def build_matrix(paths: list[Path]) -> dict:
    by_label_arb, by_label_sample, checked_ratio_max = build_raw_ratios(paths)
    labels = sorted(by_label_sample, key=label_to_t)
    coefficient_k_max = checked_ratio_max + 1
    raw_rows = checked_ratio_max * len(labels)

    raw_lower_positive_rows = 0
    raw_upper_positive_rows = 0
    half_width_positive_rows = 0
    one_third_positive_rows = 0
    min_raw_lower_margin: tuple[Decimal, str, int] | None = None
    min_raw_upper_slack: tuple[Decimal, str, int] | None = None
    min_raw_cone_margin: tuple[Decimal, str, int] | None = None
    max_raw_ratio: tuple[Decimal, str, int] | None = None

    for label in labels:
        for k in range(1, checked_ratio_max + 1):
            raw_arb = by_label_arb[label][k]
            raw_sample = by_label_sample[label][k]
            upper_wall_arb = raw_ratio_factor_arb(k)
            upper_wall_sample = raw_ratio_factor_decimal(k)
            lower_margin_arb = raw_arb - flint.arb(1)
            upper_slack_arb = upper_wall_arb - raw_arb
            lower_margin_sample = raw_sample - Decimal(1)
            upper_slack_sample = upper_wall_sample - raw_sample
            if arb_positive(lower_margin_arb):
                raw_lower_positive_rows += 1
            if arb_positive(upper_slack_arb):
                raw_upper_positive_rows += 1
            if arb_positive(raw_arb - raw_half_width_threshold_arb(k)):
                half_width_positive_rows += 1
            if arb_positive(raw_arb - raw_one_third_threshold_arb(k)):
                one_third_positive_rows += 1
            min_raw_lower_margin = update_min(min_raw_lower_margin, lower_margin_sample, label, k)
            min_raw_upper_slack = update_min(min_raw_upper_slack, upper_slack_sample, label, k)
            min_raw_cone_margin = update_min(
                min_raw_cone_margin,
                min(lower_margin_sample, upper_slack_sample),
                label,
                k,
            )
            max_raw_ratio = update_max(max_raw_ratio, raw_sample, label, k)

    corridor_rows = 0
    bridge_lower_positive_rows = 0
    scaled_upper_positive_rows = 0
    corridor_occupancy_positive_rows = 0
    corridor_width_positive_rows = 0
    min_bridge_lower_margin: tuple[Decimal, str, int] | None = None
    min_scaled_upper_margin: tuple[Decimal, str, int] | None = None
    min_corridor_occupancy_margin: tuple[Decimal, str, int] | None = None
    min_corridor_width: tuple[Decimal, str, int] | None = None

    for label in labels:
        for k in range(1, checked_ratio_max):
            corridor_rows += 1
            raw_k_arb = by_label_arb[label][k]
            raw_next_arb = by_label_arb[label][k + 1]
            raw_k_sample = by_label_sample[label][k]
            raw_next_sample = by_label_sample[label][k + 1]

            lower_factor_arb = bridge_lower_factor_arb(k)
            lower_factor_sample = bridge_lower_factor_decimal(k)
            lower_bound_arb = lower_factor_arb * raw_k_arb
            lower_bound_sample = lower_factor_sample * raw_k_sample
            upper_bound_arb = (flint.arb(2) + flint.arb(2 * k - 1) * raw_k_arb) / flint.arb(2 * k + 1)
            upper_bound_sample = (Decimal(2) + Decimal(2 * k - 1) * raw_k_sample) / Decimal(2 * k + 1)

            bridge_lower_margin_arb = raw_next_arb - lower_bound_arb
            bridge_lower_margin_sample = raw_next_sample - lower_bound_sample
            scaled_upper_margin_arb = upper_bound_arb - raw_next_arb
            scaled_upper_margin_sample = upper_bound_sample - raw_next_sample
            corridor_width_arb = upper_bound_arb - lower_bound_arb
            corridor_width_sample = upper_bound_sample - lower_bound_sample

            if arb_positive(bridge_lower_margin_arb):
                bridge_lower_positive_rows += 1
            if arb_positive(scaled_upper_margin_arb):
                scaled_upper_positive_rows += 1
            if arb_positive(bridge_lower_margin_arb) and arb_positive(scaled_upper_margin_arb):
                corridor_occupancy_positive_rows += 1
            if arb_positive(corridor_width_arb):
                corridor_width_positive_rows += 1
            min_bridge_lower_margin = update_min(min_bridge_lower_margin, bridge_lower_margin_sample, label, k)
            min_scaled_upper_margin = update_min(min_scaled_upper_margin, scaled_upper_margin_sample, label, k)
            min_corridor_occupancy_margin = update_min(
                min_corridor_occupancy_margin,
                min(bridge_lower_margin_sample, scaled_upper_margin_sample),
                label,
                k,
            )
            min_corridor_width = update_min(min_corridor_width, corridor_width_sample, label, k)

    half_width_failure_rows = raw_rows - half_width_positive_rows
    one_third_failure_rows = raw_rows - one_third_positive_rows
    summary = {
        "matrix_rows": 8,
        "ready_to_apply_rows": 0,
        "coefficient_k_max": coefficient_k_max,
        "checked_ratio_max": checked_ratio_max,
        "lambdas": labels,
        "raw_ratio_rows": raw_rows,
        "raw_lower_positive_rows": raw_lower_positive_rows,
        "raw_upper_positive_rows": raw_upper_positive_rows,
        "raw_cone_positive_rows": min(raw_lower_positive_rows, raw_upper_positive_rows),
        "half_width_positive_rows": half_width_positive_rows,
        "half_width_failure_rows": half_width_failure_rows,
        "one_third_positive_rows": one_third_positive_rows,
        "one_third_failure_rows": one_third_failure_rows,
        "corridor_rows": corridor_rows,
        "bridge_lower_positive_rows": bridge_lower_positive_rows,
        "scaled_upper_positive_rows": scaled_upper_positive_rows,
        "corridor_occupancy_positive_rows": corridor_occupancy_positive_rows,
        "corridor_width_positive_rows": corridor_width_positive_rows,
        "max_raw_ratio": asdict(extremum(max_raw_ratio)),
        "min_raw_lower_margin": asdict(extremum(min_raw_lower_margin)),
        "min_raw_upper_slack": asdict(extremum(min_raw_upper_slack)),
        "min_raw_cone_margin": asdict(extremum(min_raw_cone_margin)),
        "min_bridge_lower_margin": asdict(extremum(min_bridge_lower_margin)),
        "min_scaled_upper_margin": asdict(extremum(min_scaled_upper_margin)),
        "min_corridor_occupancy_margin": asdict(extremum(min_corridor_occupancy_margin)),
        "min_corridor_width": asdict(extremum(min_corridor_width)),
        "target_closing": False,
        "main_finding": (
            "The adaptive exact-cone route becomes a raw-moment ratio corridor. On the k200 "
            f"negative-lambda prefix, the raw cone holds on {raw_lower_positive_rows}/{raw_rows} lower "
            f"rows and {raw_upper_positive_rows}/{raw_rows} upper rows, while the adaptive corridor "
            f"contains R_(k+1) on {corridor_occupancy_positive_rows}/{corridor_rows} adjacent rows. "
            "The exact corridor is nonempty precisely under the upper raw wall, so the remaining all-k "
            "burden is a zeta-specific upper raw moment-growth theorem plus a corridor occupancy theorem."
        ),
    }
    rows = [
        {
            "id": "nlrmb_01_raw_ratio_reindexing",
            "role": "exact_reindexing",
            "readiness": "available_exact",
            "claim": "With M_k=mu_{2k}, A_k=M_k*k!/(2*k)!, and R_k=M_(k+1)*M_(k-1)/M_k^2, one has x_k=((2*k-1)/(2*k+1))*R_k.",
            "formula": "x_k=((2*k-1)/(2*k+1))*R_k",
            "proof_boundary": "Exact reindexing only; it does not prove any raw moment bound.",
        },
        {
            "id": "nlrmb_02_exact_cone_translation",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The adaptive exact cone 0<=s_k<=1 is equivalent to 1<=R_k<=(2*k+1)/(2*k-1).",
            "formula": "s_k=((2*k+1)-(2*k-1)*R_k)/2",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_adaptive_envelope_obligations.md",
                "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
            ],
            "proof_boundary": "Exact translation only; the upper raw wall remains an open all-k estimate.",
        },
        {
            "id": "nlrmb_03_adaptive_corridor_identity",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The monotone bridge and scaled k-increase are exactly a lower/upper corridor for R_(k+1), and the corridor width is positive exactly under the upper raw wall.",
            "formula": "((2*k-1)*(2*k+3)/(2*k+1)^2)*R_k <= R_(k+1) <= (2+(2*k-1)*R_k)/(2*k+1)",
            "proof_boundary": "Exact algebraic corridor only; no all-k corridor occupancy is proved.",
        },
        {
            "id": "nlrmb_04_raw_cone_prefix",
            "role": "finite_certificate",
            "readiness": "not_ready_to_apply",
            "claim": f"The checked raw ratios satisfy 1<=R_k<=(2*k+1)/(2*k-1) on {min(raw_lower_positive_rows, raw_upper_positive_rows)}/{raw_rows} rows.",
            "proof_boundary": "Finite prefix certificate only; not an all-k raw moment-growth theorem.",
        },
        {
            "id": "nlrmb_05_raw_corridor_prefix",
            "role": "finite_pattern",
            "readiness": "not_ready_to_apply",
            "claim": f"The checked adjacent raw ratios occupy the adaptive corridor on {corridor_occupancy_positive_rows}/{corridor_rows} rows.",
            "proof_boundary": "Finite corridor pattern only; it does not prove the adaptive envelope.",
        },
        {
            "id": "nlrmb_06_fixed_buffer_thresholds",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "claim": (
                f"The half-width threshold R_k>=2*k/(2*k-1) fails on {half_width_failure_rows} rows, "
                f"and the one-third threshold R_k>=(6*k+1)/(3*(2*k-1)) fails on {one_third_failure_rows} rows."
            ),
            "proof_boundary": "Finite rejection of fixed-buffer shortcuts only; not a proof of the adaptive replacement.",
        },
        {
            "id": "nlrmb_07_live_theorem_burden",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim": "A promoted route must prove the zeta-specific upper raw wall and all-k corridor occupancy without using endpoint PF, RH, or Lambda <= 0.",
            "proof_boundary": "Theorem-search burden only; no all-k theorem is supplied here.",
        },
        {
            "id": "nlrmb_08_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "Any proof upgrade must preserve the fixed-buffer rejection gates and keep exact, finite, conditional, and conjectural claims separated.",
            "proof_boundary": "Proof-hygiene gate only; not a mathematical proof.",
        },
    ]
    return {
        "kind": "jensen_window_pf_negative_lambda_raw_moment_bridge_matrix",
        "date": "2026-07-06",
        "status": "exact finite theorem-search diagnostic",
        "source_adaptive_obligations": "outputs/jensen_window_pf_negative_lambda_adaptive_envelope_obligations.md",
        "source_adaptive_matrix": "outputs/jensen_window_pf_negative_lambda_adaptive_envelope_matrix.md",
        "source_cone_prefix": "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_k200_scout.md",
        "source_boundary_threshold": "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
        "source_monotone_contraction_target": "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "source_bounded_log_curvature_target": "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.py",
        "proof_boundary": (
            "Exact algebra plus finite theorem-search diagnostics only. The artifact translates the "
            "adaptive envelope route into raw moment-ratio inequalities and checks the finite "
            "negative-lambda prefix, but it does not prove all-k raw moment-growth bounds, corridor "
            "occupancy, cone entry, jwpf_06, RH, or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "Only exact algebraic rows are available_exact.",
            "No finite row is ready_to_apply.",
            "The lower raw wall R_k>=1 is not the bottleneck; the upper raw wall remains open all-k.",
            "The adaptive corridor occupancy remains open all-k.",
            "The half-width and one-third fixed buffers remain finite-rejected.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    result_line = (
        "validated Jensen-window PF negative-lambda raw-moment bridge matrix: "
        f"{summary['matrix_rows']} matrix rows, 0 issues, "
        f"{summary['raw_cone_positive_rows']} raw-cone rows, "
        f"{summary['corridor_occupancy_positive_rows']} corridor rows, "
        f"{summary['half_width_failure_rows']} half-width failures"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Raw-Moment Bridge Matrix",
        "",
        "Date: 2026-07-06",
        "",
        "Status: exact finite theorem-search diagnostic. This is not a proof of",
        "the adaptive scaled-defect target, cone entry, Jensen-window PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_raw_moment_bridge_matrix`.",
        "",
        "Proof boundary: this artifact translates the adaptive route into raw",
        "moment-ratio inequalities and checks the finite negative-lambda prefix.",
        "It does not prove the all-`k` raw moment-growth or corridor theorem.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.json",
        "```",
        "",
        "Input enclosures:",
        "",
        "```text",
        *artifact["enclosure_jsonl"],
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Raw-Moment Translation",
        "",
        "Let:",
        "",
        "```text",
        "M_k = mu_{2k}",
        "A_k = M_k*k!/(2*k)!",
        "R_k = M_(k+1)*M_(k-1)/M_k^2",
        "x_k = ((2*k-1)/(2*k+1))*R_k",
        "s_k = ((2*k+1)-(2*k-1)*R_k)/2",
        "```",
        "",
        "Then the exact cone is:",
        "",
        "```text",
        "0 <= s_k <= 1 iff 1 <= R_k <= (2*k+1)/(2*k-1)",
        "```",
        "",
        "The adaptive corridor is:",
        "",
        "```text",
        "((2*k-1)*(2*k+3)/(2*k+1)^2)*R_k <= R_(k+1) <= (2+(2*k-1)*R_k)/(2*k+1)",
        "corridor width = 2*(2*k+1-(2*k-1)*R_k)/(2*k+1)^2",
        "```",
        "",
        "Finite diagnostics:",
        "",
        "```text",
        f"lambdas: {', '.join(summary['lambdas'])}",
        f"coefficient range: A_0..A_{summary['coefficient_k_max']}",
        f"checked raw ratios: R_1..R_{summary['checked_ratio_max']}",
        f"raw lower wall rows: {summary['raw_lower_positive_rows']} / {summary['raw_ratio_rows']}",
        f"raw upper wall rows: {summary['raw_upper_positive_rows']} / {summary['raw_ratio_rows']}",
        f"corridor occupancy rows: {summary['corridor_occupancy_positive_rows']} / {summary['corridor_rows']}",
        f"corridor width rows: {summary['corridor_width_positive_rows']} / {summary['corridor_rows']}",
        f"half-width failure rows: {summary['half_width_failure_rows']}",
        f"one-third failure rows: {summary['one_third_failure_rows']}",
        "```",
        "",
        "Extrema:",
        "",
        "```text",
        (
            "max raw ratio: "
            f"{summary['max_raw_ratio']['sample']} at lambda={summary['max_raw_ratio']['lam']}, "
            f"k={summary['max_raw_ratio']['k']}"
        ),
        (
            "min raw lower margin: "
            f"{summary['min_raw_lower_margin']['sample']} at lambda={summary['min_raw_lower_margin']['lam']}, "
            f"k={summary['min_raw_lower_margin']['k']}"
        ),
        (
            "min raw upper slack: "
            f"{summary['min_raw_upper_slack']['sample']} at lambda={summary['min_raw_upper_slack']['lam']}, "
            f"k={summary['min_raw_upper_slack']['k']}"
        ),
        (
            "min bridge lower margin: "
            f"{summary['min_bridge_lower_margin']['sample']} at lambda={summary['min_bridge_lower_margin']['lam']}, "
            f"k={summary['min_bridge_lower_margin']['k']}"
        ),
        (
            "min scaled upper margin: "
            f"{summary['min_scaled_upper_margin']['sample']} at lambda={summary['min_scaled_upper_margin']['lam']}, "
            f"k={summary['min_scaled_upper_margin']['k']}"
        ),
        (
            "min corridor width: "
            f"{summary['min_corridor_width']['sample']} at lambda={summary['min_corridor_width']['lam']}, "
            f"k={summary['min_corridor_width']['k']}"
        ),
        "```",
        "",
        "Fixed-buffer translations:",
        "",
        "```text",
        "s_k <= 1/2 iff R_k >= 2*k/(2*k-1)",
        "s_k <= 1/3 iff R_k >= (6*k+1)/(3*(2*k-1))",
        "```",
        "",
        "Integration:",
        "",
        "```text",
        artifact["source_adaptive_obligations"],
        artifact["source_adaptive_matrix"],
        artifact["source_cone_prefix"],
        artifact["source_boundary_threshold"],
        artifact["source_monotone_contraction_target"],
        artifact["source_bounded_log_curvature_target"],
        "outputs/jensen_window_pf_negative_lambda_half_width_tail_target.md",
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
    parser.add_argument("--enclosures", type=Path, nargs="+", default=list(DEFAULT_ENCLOSURE_JSONL))
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    paths = [path if path.is_absolute() else REPO_ROOT / path for path in args.enclosures]
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_matrix(paths)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda raw-moment bridge matrix: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
