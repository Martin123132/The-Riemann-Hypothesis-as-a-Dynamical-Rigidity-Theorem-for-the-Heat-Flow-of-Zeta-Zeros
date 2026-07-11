#!/usr/bin/env python3
"""Build a raw-ratio decrement-corridor scout for the negative-lambda route."""

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
from jensen_window_pf_negative_lambda_raw_moment_bridge_matrix import (  # noqa: E402
    DEFAULT_ENCLOSURE_JSONL,
    build_raw_ratios,
)


getcontext().prec = 100

DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.md"


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
        raise RuntimeError("missing decrement-corridor extremum")
    value, lam, k = item
    return Extremum(decimal_format(value), lam, k)


def label_to_t(label: str) -> Decimal:
    return abs(Decimal(label))


def decrement_lower_arb(k: int, raw: flint.arb) -> flint.arb:
    return flint.arb(2) * (raw - flint.arb(1)) / flint.arb(2 * k + 1)


def decrement_upper_arb(k: int, raw: flint.arb) -> flint.arb:
    return flint.arb(4) * raw / flint.arb((2 * k + 1) ** 2)


def decrement_lower_decimal(k: int, raw: Decimal) -> Decimal:
    return Decimal(2) * (raw - Decimal(1)) / Decimal(2 * k + 1)


def decrement_upper_decimal(k: int, raw: Decimal) -> Decimal:
    return Decimal(4) * raw / Decimal((2 * k + 1) ** 2)


def theta_arb(k: int, raw: flint.arb, raw_next: flint.arb) -> flint.arb:
    decrement = raw - raw_next
    lower = decrement_lower_arb(k, raw)
    upper = decrement_upper_arb(k, raw)
    return (decrement - lower) / (upper - lower)


def theta_decimal(k: int, raw: Decimal, raw_next: Decimal) -> Decimal:
    decrement = raw - raw_next
    lower = decrement_lower_decimal(k, raw)
    upper = decrement_upper_decimal(k, raw)
    return (decrement - lower) / (upper - lower)


def q(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def raw_upper_wall(k: int) -> Fraction:
    return Fraction(2 * k + 1, 2 * k - 1)


def exact_counterexample(raw: Fraction, raw_next: Fraction, k: int = 1) -> dict:
    decrement = raw - raw_next
    lower = Fraction(2) * (raw - 1) / Fraction(2 * k + 1)
    upper = Fraction(4) * raw / Fraction((2 * k + 1) ** 2)
    corridor_lower = raw - upper
    corridor_upper = raw - lower
    return {
        "k": k,
        "R_k": q(raw),
        "R_next": q(raw_next),
        "R_k_upper_wall": q(raw_upper_wall(k)),
        "R_next_upper_wall": q(raw_upper_wall(k + 1)),
        "D_k": q(decrement),
        "decrement_lower": q(lower),
        "decrement_upper": q(upper),
        "lower_margin_D_minus_lower": q(decrement - lower),
        "upper_margin_upper_minus_D": q(upper - decrement),
        "corridor_lower": q(corridor_lower),
        "corridor_upper": q(corridor_upper),
        "raw_cone_holds_at_k_and_next": bool(
            Fraction(1) <= raw <= raw_upper_wall(k)
            and Fraction(1) <= raw_next <= raw_upper_wall(k + 1)
        ),
        "raw_decrease_holds": bool(raw_next < raw),
    }


def build_diagnostics(paths: list[Path]) -> dict:
    by_label_arb, by_label_sample, checked_ratio_max = build_raw_ratios(paths)
    labels = sorted(by_label_sample, key=label_to_t)
    adjacent_rows = (checked_ratio_max - 1) * len(labels)

    decrease_rows = 0
    lower_decrement_rows = 0
    upper_decrement_rows = 0
    corridor_rows = 0
    theta_rows = 0

    min_decrement: tuple[Decimal, str, int] | None = None
    max_decrement: tuple[Decimal, str, int] | None = None
    min_lower_margin: tuple[Decimal, str, int] | None = None
    min_upper_margin: tuple[Decimal, str, int] | None = None
    min_width: tuple[Decimal, str, int] | None = None
    min_theta: tuple[Decimal, str, int] | None = None
    max_theta: tuple[Decimal, str, int] | None = None
    min_theta_k_drop: tuple[Decimal, str, int] | None = None
    min_theta_lambda_gap: tuple[Decimal, str, int] | None = None

    theta_by_label_arb: dict[str, dict[int, flint.arb]] = {}
    theta_by_label_sample: dict[str, dict[int, Decimal]] = {}

    for label in labels:
        theta_by_label_arb[label] = {}
        theta_by_label_sample[label] = {}
        for k in range(1, checked_ratio_max):
            raw_arb = by_label_arb[label][k]
            next_arb = by_label_arb[label][k + 1]
            raw_sample = by_label_sample[label][k]
            next_sample = by_label_sample[label][k + 1]

            decrement_arb = raw_arb - next_arb
            lower_arb = decrement_lower_arb(k, raw_arb)
            upper_arb = decrement_upper_arb(k, raw_arb)
            lower_margin_arb = decrement_arb - lower_arb
            upper_margin_arb = upper_arb - decrement_arb
            width_arb = upper_arb - lower_arb
            theta_value_arb = (decrement_arb - lower_arb) / width_arb

            decrement_sample = raw_sample - next_sample
            lower_sample = decrement_lower_decimal(k, raw_sample)
            upper_sample = decrement_upper_decimal(k, raw_sample)
            lower_margin_sample = decrement_sample - lower_sample
            upper_margin_sample = upper_sample - decrement_sample
            width_sample = upper_sample - lower_sample
            theta_value_sample = (decrement_sample - lower_sample) / width_sample

            theta_by_label_arb[label][k] = theta_value_arb
            theta_by_label_sample[label][k] = theta_value_sample

            if arb_positive(decrement_arb):
                decrease_rows += 1
            if arb_positive(lower_margin_arb):
                lower_decrement_rows += 1
            if arb_positive(upper_margin_arb):
                upper_decrement_rows += 1
            if arb_positive(lower_margin_arb) and arb_positive(upper_margin_arb):
                corridor_rows += 1
            if arb_positive(theta_value_arb) and arb_positive(flint.arb(1) - theta_value_arb):
                theta_rows += 1

            min_decrement = update_min(min_decrement, decrement_sample, label, k)
            max_decrement = update_max(max_decrement, decrement_sample, label, k)
            min_lower_margin = update_min(min_lower_margin, lower_margin_sample, label, k)
            min_upper_margin = update_min(min_upper_margin, upper_margin_sample, label, k)
            min_width = update_min(min_width, width_sample, label, k)
            min_theta = update_min(min_theta, theta_value_sample, label, k)
            max_theta = update_max(max_theta, theta_value_sample, label, k)

    theta_k_monotone_rows = 0
    for label in labels:
        for k in range(1, checked_ratio_max - 1):
            drop_arb = theta_by_label_arb[label][k] - theta_by_label_arb[label][k + 1]
            drop_sample = theta_by_label_sample[label][k] - theta_by_label_sample[label][k + 1]
            if arb_positive(drop_arb):
                theta_k_monotone_rows += 1
            min_theta_k_drop = update_min(min_theta_k_drop, drop_sample, label, k)

    theta_lambda_order_rows = 0
    for left, right in zip(labels, labels[1:]):
        pair_label = f"{left}->{right}"
        for k in range(1, checked_ratio_max):
            gap_arb = theta_by_label_arb[right][k] - theta_by_label_arb[left][k]
            gap_sample = theta_by_label_sample[right][k] - theta_by_label_sample[left][k]
            if arb_positive(gap_arb):
                theta_lambda_order_rows += 1
            min_theta_lambda_gap = update_min(min_theta_lambda_gap, gap_sample, pair_label, k)

    return {
        "lambdas": labels,
        "coefficient_k_max": checked_ratio_max + 1,
        "checked_ratio_max": checked_ratio_max,
        "adjacent_rows": adjacent_rows,
        "raw_decrease_rows": decrease_rows,
        "lower_decrement_rows": lower_decrement_rows,
        "upper_decrement_rows": upper_decrement_rows,
        "decrement_corridor_rows": corridor_rows,
        "theta_unit_rows": theta_rows,
        "theta_k_monotone_rows": theta_k_monotone_rows,
        "theta_k_monotone_total_rows": (checked_ratio_max - 2) * len(labels),
        "theta_lambda_order_rows": theta_lambda_order_rows,
        "theta_lambda_order_total_rows": (len(labels) - 1) * (checked_ratio_max - 1),
        "min_decrement": asdict(extremum(min_decrement)),
        "max_decrement": asdict(extremum(max_decrement)),
        "min_lower_decrement_margin": asdict(extremum(min_lower_margin)),
        "min_upper_decrement_margin": asdict(extremum(min_upper_margin)),
        "min_decrement_corridor_width": asdict(extremum(min_width)),
        "min_theta": asdict(extremum(min_theta)),
        "max_theta": asdict(extremum(max_theta)),
        "min_theta_k_drop": asdict(extremum(min_theta_k_drop)),
        "min_theta_lambda_gap": asdict(extremum(min_theta_lambda_gap)),
    }


def build_artifact(paths: list[Path] | None = None) -> dict:
    if paths is None:
        paths = list(DEFAULT_ENCLOSURE_JSONL)
    diagnostics = build_diagnostics(paths)
    upper_side_failure = exact_counterexample(Fraction(2), Fraction(3, 2))
    lower_side_failure = exact_counterexample(Fraction(2), Fraction(1))
    rows = [
        {
            "id": "nlrdc_01_exact_decrement_corridor",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The raw adaptive corridor is exactly a two-sided decrement bound for D_k=R_k-R_(k+1).",
            "formula": "2*(R_k-1)/(2*k+1) <= D_k <= 4*R_k/(2*k+1)^2",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.md",
                "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
            ],
            "proof_boundary": "Exact algebra only; no zeta-specific inequality is proved.",
        },
        {
            "id": "nlrdc_02_exact_width_equivalence",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The decrement corridor width is positive exactly under the upper raw wall.",
            "formula": "4*R_k/(2*k+1)^2 - 2*(R_k-1)/(2*k+1) = 2*(2*k+1-(2*k-1)*R_k)/(2*k+1)^2",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.md",
                "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
            ],
            "proof_boundary": "Exact algebra only; proving the upper raw wall for zeta remains open.",
        },
        {
            "id": "nlrdc_03_k200_decrement_corridor",
            "role": "finite_pattern",
            "readiness": "not_ready_to_apply",
            "claim": (
                "The k200 negative-lambda prefix satisfies raw decrease, both decrement walls, "
                f"and decrement-corridor occupancy on {diagnostics['decrement_corridor_rows']}/{diagnostics['adjacent_rows']} adjacent rows."
            ),
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.md",
                "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_k200_scout.md",
            ],
            "proof_boundary": "Finite prefix evidence only; not an all-k recurrence theorem.",
        },
        {
            "id": "nlrdc_04_theta_coordinate_shape",
            "role": "finite_pattern",
            "readiness": "not_ready_to_apply",
            "claim": (
                "The normalized corridor coordinate theta_k=(D_k-L_k)/(U_k-L_k) stays in (0,1), "
                "decreases with k on each checked lambda, and increases as lambda becomes more negative."
            ),
            "formula": "theta_k=(R_k-R_(k+1)-2*(R_k-1)/(2*k+1))/(4*R_k/(2*k+1)^2-2*(R_k-1)/(2*k+1))",
            "proof_boundary": "Finite shape diagnostic only; not a monotone theta theorem.",
        },
        {
            "id": "nlrdc_05_monotone_decrease_shortcut_blocked",
            "role": "exact_counterexample",
            "readiness": "not_ready_to_apply",
            "claim": "Raw cone occupancy plus R_(k+1)<R_k does not imply the scaled-upper side of the decrement corridor.",
            "witness": upper_side_failure,
            "proof_boundary": "Exact counterexample gate only; not evidence against the zeta prefix.",
        },
        {
            "id": "nlrdc_06_overfast_decrease_shortcut_blocked",
            "role": "exact_counterexample",
            "readiness": "not_ready_to_apply",
            "claim": "Raw cone occupancy plus R_(k+1)<R_k also does not imply the monotone-bridge lower side of the decrement corridor.",
            "witness": lower_side_failure,
            "proof_boundary": "Exact counterexample gate only; not evidence against the zeta prefix.",
        },
        {
            "id": "nlrdc_07_live_theta_recurrence_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim": "A zeta-specific proof could promote the finite theta-shape into a recurrence or comparison principle that forces 0<theta_k<1 on the tail.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md",
                "outputs/jensen_window_pf_negative_lambda_adaptive_envelope_obligations.md",
            ],
            "gap": "No all-k formula or interval-safe comparison for theta_k is proved.",
            "proof_boundary": "Live theorem-search route only.",
        },
        {
            "id": "nlrdc_08_open_requirements",
            "role": "open_requirement",
            "readiness": "not_ready_to_apply",
            "claim": "A promoted route must prove the upper raw wall and the two decrement inequalities for actual zeta moments, with a finite/collar handoff and no endpoint PF or RH input.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
                "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
            ],
            "proof_boundary": "Open analytic requirement only.",
        },
        {
            "id": "nlrdc_09_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "Any recurrence upgrade must preserve the exact counterexamples: raw cone, Stieltjes positivity, and monotone decrease are insufficient by themselves.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
                "outputs/signed_hankel_jensen_dependency_graph.md",
            ],
            "proof_boundary": "Proof-hygiene gate only; not a proof.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "exact_available_rows": 2,
        "finite_pattern_rows": 2,
        "exact_counterexample_rows": 2,
        "live_routes": 1,
        "open_requirement_rows": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        **diagnostics,
        "main_finding": (
            "The zeta-specific raw-corridor target can be attacked as a decrement recurrence: "
            "prove 2*(R_k-1)/(2*k+1) <= R_k-R_(k+1) <= 4*R_k/(2*k+1)^2. "
            f"The k200 prefix satisfies this on {diagnostics['decrement_corridor_rows']} adjacent rows, "
            f"and the normalized theta coordinate is k-monotone on {diagnostics['theta_k_monotone_rows']} rows "
            f"and lambda-ordered on {diagnostics['theta_lambda_order_rows']} rows. Exact raw-cone monotone "
            "counterexamples show that the recurrence must use zeta-specific structure."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout",
        "date": "2026-07-07",
        "status": "exact finite theorem-search diagnostic",
        "source_raw_bridge": "outputs/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.md",
        "source_raw_corridor_target": "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
        "source_raw_obstruction": "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
        "source_adaptive_obligations": "outputs/jensen_window_pf_negative_lambda_adaptive_envelope_obligations.md",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.py",
        "proof_boundary": (
            "Exact algebra plus finite theorem-search diagnostics only. The artifact rewrites the "
            "zeta-specific raw corridor as a decrement recurrence and records finite theta-shape "
            "evidence plus exact counterexamples, but it does not prove the raw-corridor theorem, "
            "does not prove cone entry, does not prove jwpf_06, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "Only exact algebra rows are available_exact.",
            "No finite row is ready_to_apply.",
            "The decrement-corridor recurrence remains zeta-specific and open all-k.",
            "Raw cone occupancy plus raw monotone decrease is explicitly blocked as a shortcut.",
            "Theta monotonicity is finite shape evidence only.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    result_line = (
        "validated Jensen-window PF negative-lambda raw-ratio decrement-corridor scout: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['decrement_corridor_rows']} decrement-corridor rows, "
        f"{summary['theta_k_monotone_rows']} theta-k-monotone rows, "
        f"{summary['exact_counterexample_rows']} exact counterexamples, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Raw-Ratio Decrement-Corridor Scout",
        "",
        "Date: 2026-07-07",
        "",
        "Status: exact finite theorem-search diagnostic. This is not a proof",
        "of the raw-corridor theorem, cone entry, Jensen-window PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout`.",
        "",
        "Proof boundary: this artifact rewrites the zeta-specific raw corridor",
        "as a decrement recurrence and validates finite shape evidence. It does",
        "not prove an all-`k` recurrence.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.json",
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
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Exact Decrement Form",
        "",
        "Let:",
        "",
        "```text",
        "M_k = mu_{2k}",
        "R_k = M_(k+1)*M_(k-1)/M_k^2",
        "D_k = R_k - R_(k+1)",
        "```",
        "",
        "The adaptive raw corridor is exactly:",
        "",
        "```text",
        "2*(R_k-1)/(2*k+1) <= D_k <= 4*R_k/(2*k+1)^2",
        "```",
        "",
        "Its width is:",
        "",
        "```text",
        "4*R_k/(2*k+1)^2 - 2*(R_k-1)/(2*k+1)",
        "= 2*(2*k+1-(2*k-1)*R_k)/(2*k+1)^2",
        "```",
        "",
        "Finite diagnostics:",
        "",
        "```text",
        f"lambdas: {', '.join(summary['lambdas'])}",
        f"coefficient range: A_0..A_{summary['coefficient_k_max']}",
        f"checked raw ratios: R_1..R_{summary['checked_ratio_max']}",
        f"raw decrease rows: {summary['raw_decrease_rows']} / {summary['adjacent_rows']}",
        f"lower decrement rows: {summary['lower_decrement_rows']} / {summary['adjacent_rows']}",
        f"upper decrement rows: {summary['upper_decrement_rows']} / {summary['adjacent_rows']}",
        f"decrement-corridor rows: {summary['decrement_corridor_rows']} / {summary['adjacent_rows']}",
        f"theta unit rows: {summary['theta_unit_rows']} / {summary['adjacent_rows']}",
        f"theta-k monotone rows: {summary['theta_k_monotone_rows']} / {summary['theta_k_monotone_total_rows']}",
        f"theta lambda-order rows: {summary['theta_lambda_order_rows']} / {summary['theta_lambda_order_total_rows']}",
        "```",
        "",
        "Extrema:",
        "",
        "```text",
        (
            "min lower decrement margin: "
            f"{summary['min_lower_decrement_margin']['sample']} at "
            f"lambda={summary['min_lower_decrement_margin']['lam']}, "
            f"k={summary['min_lower_decrement_margin']['k']}"
        ),
        (
            "min upper decrement margin: "
            f"{summary['min_upper_decrement_margin']['sample']} at "
            f"lambda={summary['min_upper_decrement_margin']['lam']}, "
            f"k={summary['min_upper_decrement_margin']['k']}"
        ),
        (
            "min decrement corridor width: "
            f"{summary['min_decrement_corridor_width']['sample']} at "
            f"lambda={summary['min_decrement_corridor_width']['lam']}, "
            f"k={summary['min_decrement_corridor_width']['k']}"
        ),
        (
            "theta range: "
            f"{summary['min_theta']['sample']} at lambda={summary['min_theta']['lam']}, k={summary['min_theta']['k']} "
            f"to {summary['max_theta']['sample']} at lambda={summary['max_theta']['lam']}, k={summary['max_theta']['k']}"
        ),
        (
            "min theta k-drop: "
            f"{summary['min_theta_k_drop']['sample']} at "
            f"lambda={summary['min_theta_k_drop']['lam']}, k={summary['min_theta_k_drop']['k']}"
        ),
        (
            "min theta lambda gap: "
            f"{summary['min_theta_lambda_gap']['sample']} at "
            f"pair={summary['min_theta_lambda_gap']['lam']}, k={summary['min_theta_lambda_gap']['k']}"
        ),
        "```",
        "",
        "Exact shortcut blockers:",
        "",
        "```text",
        "R_1=2, R_2=3/2: raw cone and decrease hold, but the lower decrement wall fails",
        "R_1=2, R_2=1: raw cone and decrease hold, but the upper decrement wall fails",
        "```",
        "",
        "Integration:",
        "",
        "```text",
        artifact["source_raw_bridge"],
        artifact["source_raw_corridor_target"],
        artifact["source_raw_obstruction"],
        artifact["source_adaptive_obligations"],
        "outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md",
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
    parser.add_argument("--enclosures", type=Path, nargs="+", default=list(DEFAULT_ENCLOSURE_JSONL))
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    paths = [path if path.is_absolute() else REPO_ROOT / path for path in args.enclosures]
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(paths)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda raw-ratio decrement-corridor scout: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
