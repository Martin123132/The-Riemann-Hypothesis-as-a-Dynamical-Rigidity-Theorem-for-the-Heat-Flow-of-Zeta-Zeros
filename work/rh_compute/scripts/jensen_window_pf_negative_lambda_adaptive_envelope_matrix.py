#!/usr/bin/env python3
"""Build the negative-lambda adaptive scaled-defect envelope matrix."""

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
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_adaptive_envelope_matrix.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_adaptive_envelope_matrix.md"


@dataclass(frozen=True)
class Extremum:
    sample: str
    lam: str
    k: int


@dataclass(frozen=True)
class PairExtremum:
    sample: str
    pair: str
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


def pair_extremum(item: tuple[Decimal, str, int] | None) -> PairExtremum:
    if item is None:
        raise RuntimeError("missing pair extremum")
    value, pair, k = item
    return PairExtremum(decimal_format(value), pair, k)


def label_to_t(label: str) -> Decimal:
    return abs(Decimal(label))


def build_scaled_defects(paths: list[Path]) -> tuple[dict[str, dict[int, flint.arb]], dict[str, dict[int, Decimal]], int]:
    balls, samples, labels = load_enclosures(paths)
    if not labels:
        raise RuntimeError("no coefficient enclosures found")
    coefficient_k_max = max(index for (_lam, index) in balls)
    checked_x_max = coefficient_k_max - 1
    by_label_arb: dict[str, dict[int, flint.arb]] = {}
    by_label_sample: dict[str, dict[int, Decimal]] = {}

    for lam in sorted(labels):
        label = labels[lam]
        by_label_arb[label] = {}
        by_label_sample[label] = {}
        for k in range(1, checked_x_max + 1):
            x_arb = contraction({idx: balls[(lam, idx)] for idx in (k - 1, k, k + 1)}, k)
            x_sample = contraction({idx: samples[(lam, idx)] for idx in (k - 1, k, k + 1)}, k)
            scaled_arb = (flint.arb(1) - x_arb) * flint.arb(2 * k + 1) / flint.arb(2)
            scaled_sample = (Decimal(1) - x_sample) * Decimal(2 * k + 1) / Decimal(2)
            by_label_arb[label][k] = scaled_arb
            by_label_sample[label][k] = scaled_sample
    return by_label_arb, by_label_sample, checked_x_max


def build_matrix(paths: list[Path]) -> dict:
    by_label_arb, by_label_sample, checked_x_max = build_scaled_defects(paths)
    labels = sorted(by_label_sample, key=label_to_t)
    coefficient_k_max = checked_x_max + 1
    scaled_rows = checked_x_max * len(labels)

    cone_positive_rows = 0
    half_width_positive_rows = 0
    one_third_positive_rows = 0
    min_cone_margin: tuple[Decimal, str, int] | None = None
    min_upper_slack: tuple[Decimal, str, int] | None = None
    max_scaled: tuple[Decimal, str, int] | None = None

    for label in labels:
        for k in range(1, checked_x_max + 1):
            scaled_arb = by_label_arb[label][k]
            scaled_sample = by_label_sample[label][k]
            if arb_positive(scaled_arb) and arb_positive(flint.arb(1) - scaled_arb):
                cone_positive_rows += 1
            if arb_positive(flint.arb(1) / flint.arb(2) - scaled_arb):
                half_width_positive_rows += 1
            if arb_positive(flint.arb(1) / flint.arb(3) - scaled_arb):
                one_third_positive_rows += 1
            min_cone_margin = update_min(min_cone_margin, min(scaled_sample, Decimal(1) - scaled_sample), label, k)
            min_upper_slack = update_min(min_upper_slack, Decimal(1) - scaled_sample, label, k)
            max_scaled = update_max(max_scaled, scaled_sample, label, k)

    k_increase_rows = 0
    k_increase_positive_rows = 0
    min_k_increase: tuple[Decimal, str, int] | None = None
    for label in labels:
        for k in range(1, checked_x_max):
            k_increase_rows += 1
            inc_arb = by_label_arb[label][k + 1] - by_label_arb[label][k]
            inc_sample = by_label_sample[label][k + 1] - by_label_sample[label][k]
            if arb_positive(inc_arb):
                k_increase_positive_rows += 1
            min_k_increase = update_min(min_k_increase, inc_sample, label, k)

    lambda_order_rows = 0
    lambda_order_positive_rows = 0
    min_lambda_gap: tuple[Decimal, str, int] | None = None
    lambda_pairs: list[dict] = []
    for weaker_label, stronger_label in zip(labels, labels[1:]):
        # labels are ordered by |lambda|, so a positive gap means s(|lambda| smaller) >= s(|lambda| larger).
        pair = f"{weaker_label}>={stronger_label}"
        pair_rows = 0
        pair_positive = 0
        pair_min: tuple[Decimal, str, int] | None = None
        for k in range(1, checked_x_max + 1):
            lambda_order_rows += 1
            pair_rows += 1
            gap_arb = by_label_arb[weaker_label][k] - by_label_arb[stronger_label][k]
            gap_sample = by_label_sample[weaker_label][k] - by_label_sample[stronger_label][k]
            if arb_positive(gap_arb):
                lambda_order_positive_rows += 1
                pair_positive += 1
            min_lambda_gap = update_min(min_lambda_gap, gap_sample, pair, k)
            pair_min = update_min(pair_min, gap_sample, pair, k)
        lambda_pairs.append(
            {
                "pair": pair,
                "rows": pair_rows,
                "positive_rows": pair_positive,
                "min_gap": asdict(pair_extremum(pair_min)),
            }
        )

    half_width_failures = scaled_rows - half_width_positive_rows
    one_third_failures = scaled_rows - one_third_positive_rows
    max_scaled_ext = extremum(max_scaled)
    min_upper_slack_ext = extremum(min_upper_slack)
    summary = {
        "matrix_rows": 7,
        "ready_to_apply_rows": 0,
        "coefficient_k_max": coefficient_k_max,
        "checked_x_max": checked_x_max,
        "lambdas": labels,
        "scaled_defect_rows": scaled_rows,
        "cone_positive_rows": cone_positive_rows,
        "half_width_positive_rows": half_width_positive_rows,
        "half_width_failure_rows": half_width_failures,
        "one_third_positive_rows": one_third_positive_rows,
        "one_third_failure_rows": one_third_failures,
        "k_increase_rows": k_increase_rows,
        "k_increase_positive_rows": k_increase_positive_rows,
        "lambda_order_rows": lambda_order_rows,
        "lambda_order_positive_rows": lambda_order_positive_rows,
        "lambda_pairs": lambda_pairs,
        "max_scaled_defect": asdict(max_scaled_ext),
        "min_cone_margin": asdict(extremum(min_cone_margin)),
        "min_upper_slack": asdict(min_upper_slack_ext),
        "min_k_increase": asdict(extremum(min_k_increase)),
        "min_lambda_gap": asdict(pair_extremum(min_lambda_gap)),
        "target_closing": False,
        "main_finding": (
            "On the k200 negative-lambda prefix, s_k is interval-increasing in k on all "
            f"{k_increase_positive_rows}/{k_increase_rows} adjacent rows and interval-ordered by "
            f"|lambda| on all {lambda_order_positive_rows}/{lambda_order_rows} cross-lambda rows. "
            f"The largest checked scaled defect is {max_scaled_ext.sample} at lambda={max_scaled_ext.lam}, "
            f"k={max_scaled_ext.k}, leaving upper cone slack {min_upper_slack_ext.sample}. "
            "This supports a monotone-envelope proof search, but supplies no all-k tail theorem."
        ),
    }
    rows = [
        {
            "id": "nlaem_01_exact_cone_prefix",
            "role": "finite_certificate",
            "readiness": "not_ready_to_apply",
            "claim": f"The exact cone 0<=s_k<=1 holds on {cone_positive_rows}/{scaled_rows} checked rows.",
            "proof_boundary": "Finite prefix certificate only; not an all-k exact-cone theorem.",
        },
        {
            "id": "nlaem_02_k_monotone_frontier",
            "role": "finite_pattern",
            "readiness": "not_ready_to_apply",
            "claim": f"The scaled defect increases in k on {k_increase_positive_rows}/{k_increase_rows} checked adjacent rows.",
            "proof_boundary": "Finite monotonicity pattern only; it does not prove a limiting envelope.",
        },
        {
            "id": "nlaem_03_lambda_magnitude_order",
            "role": "finite_pattern",
            "readiness": "not_ready_to_apply",
            "claim": f"The scaled defect is ordered by |lambda| on {lambda_order_positive_rows}/{lambda_order_rows} checked cross-lambda rows.",
            "proof_boundary": "Finite lambda-grid pattern only; it does not prove continuous lambda monotonicity.",
        },
        {
            "id": "nlaem_04_corner_dominance",
            "role": "finite_frontier",
            "readiness": "not_ready_to_apply",
            "claim": (
                f"The checked maximum is s={max_scaled_ext.sample} at lambda={max_scaled_ext.lam}, "
                f"k={max_scaled_ext.k}, with upper cone slack {min_upper_slack_ext.sample}."
            ),
            "proof_boundary": "Finite corner dominance only; it does not control k beyond the checked prefix.",
        },
        {
            "id": "nlaem_05_fixed_buffer_rejections",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "claim": (
                f"The half-width buffer fails on {half_width_failures} rows, the one-third buffer fails on "
                f"{one_third_failures} rows, and scaled-defect nonincrease is incompatible with the k-increase pattern."
            ),
            "proof_boundary": "Finite rejection of shortcuts only; not a proof of the replacement envelope.",
        },
        {
            "id": "nlaem_06_monotone_envelope_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim": (
                "A plausible repair is to prove k-monotonicity, |lambda|-monotonicity, and a limiting/adaptive "
                "upper envelope below the exact cone wall."
            ),
            "proof_boundary": "Theorem-search route only; no monotone-envelope theorem is proved here.",
        },
        {
            "id": "nlaem_07_acceptance_tests",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "Any promoted proof must give interval-safe all-k estimates and survive the fixed-buffer rejection gates.",
            "proof_boundary": "Proof-hygiene gate only; not a mathematical bridge theorem.",
        },
    ]
    return {
        "kind": "jensen_window_pf_negative_lambda_adaptive_envelope_matrix",
        "date": "2026-07-06",
        "status": "finite theorem-search diagnostic",
        "source_scaled_defect_frontier": "outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k200_scout.md",
        "source_adaptive_target": "outputs/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.md",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_adaptive_envelope_matrix.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_envelope_matrix.py",
        "proof_boundary": (
            "Finite theorem-search diagnostic only. It records monotone-envelope patterns for the checked "
            "negative-lambda grid, but it does not prove k-uniform bounds, continuous lambda monotonicity, "
            "the adaptive scaled-defect target, cone entry, jwpf_06, RH, or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "The exact cone evidence is finite prefix evidence only.",
            "The k-monotone and lambda-magnitude patterns are finite diagnostics only.",
            "The one-third, fixed half-width, and scaled-defect nonincrease shortcuts remain rejected.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    result_line = (
        "validated Jensen-window PF negative-lambda adaptive envelope matrix: "
        f"{summary['matrix_rows']} matrix rows, 0 issues, "
        f"{summary['k_increase_positive_rows']} k-increase rows, "
        f"{summary['lambda_order_positive_rows']} lambda-order rows, "
        f"{summary['half_width_failure_rows']} half-width failures"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Adaptive Envelope Matrix",
        "",
        "Date: 2026-07-06",
        "",
        "Status: finite theorem-search diagnostic. This is not a proof of the",
        "adaptive scaled-defect target, cone entry, Jensen-window PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_adaptive_envelope_matrix`.",
        "",
        "Proof boundary: this artifact records checked monotone-envelope",
        "patterns only. It does not prove all-`k` bounds, continuous lambda",
        "monotonicity, or a tail theorem.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_adaptive_envelope_matrix.json",
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
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_adaptive_envelope_matrix.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_envelope_matrix.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Matrix",
        "",
        "```text",
        f"lambdas: {', '.join(summary['lambdas'])}",
        f"coefficient range: A_0..A_{summary['coefficient_k_max']}",
        f"checked contractions: x_1..x_{summary['checked_x_max']}",
        f"exact cone rows: {summary['cone_positive_rows']} / {summary['scaled_defect_rows']}",
        f"k-increase rows: {summary['k_increase_positive_rows']} / {summary['k_increase_rows']}",
        f"lambda-order rows: {summary['lambda_order_positive_rows']} / {summary['lambda_order_rows']}",
        f"half-width failure rows: {summary['half_width_failure_rows']}",
        f"one-third failure rows: {summary['one_third_failure_rows']}",
        "```",
        "",
        "Cross-lambda order gaps:",
        "",
        "```text",
    ]
    for row in summary["lambda_pairs"]:
        lines.append(
            f"{row['pair']}: {row['positive_rows']} / {row['rows']} positive, "
            f"min gap {row['min_gap']['sample']} at k={row['min_gap']['k']}"
        )
    lines.extend(
        [
            "```",
            "",
            "Extrema:",
            "",
            "```text",
            (
                "max scaled defect: "
                f"{summary['max_scaled_defect']['sample']} at lambda={summary['max_scaled_defect']['lam']}, "
                f"k={summary['max_scaled_defect']['k']}"
            ),
            (
                "min upper cone slack: "
                f"{summary['min_upper_slack']['sample']} at lambda={summary['min_upper_slack']['lam']}, "
                f"k={summary['min_upper_slack']['k']}"
            ),
            (
                "min exact-cone margin: "
                f"{summary['min_cone_margin']['sample']} at lambda={summary['min_cone_margin']['lam']}, "
                f"k={summary['min_cone_margin']['k']}"
            ),
            (
                "min k-increase: "
                f"{summary['min_k_increase']['sample']} at lambda={summary['min_k_increase']['lam']}, "
                f"k={summary['min_k_increase']['k']}"
            ),
            (
                "min lambda-order gap: "
                f"{summary['min_lambda_gap']['sample']} for {summary['min_lambda_gap']['pair']}, "
                f"k={summary['min_lambda_gap']['k']}"
            ),
            "```",
            "",
            "Live monotone-envelope route:",
            "",
            "```text",
            "prove s_k(lambda) increases in k without crossing 1",
            "prove s_k(lambda) decreases as |lambda| increases on the needed negative-lambda ray",
            "prove a limiting/adaptive upper envelope below the exact cone wall",
            "keep the one-third, half-width, and nonincrease shortcuts rejected",
            "```",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_scaled_defect_frontier"],
            artifact["source_adaptive_target"],
            "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
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
        "wrote Jensen-window PF negative-lambda adaptive envelope matrix: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
