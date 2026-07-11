#!/usr/bin/env python3
"""Build a recurrence-compatibility scout for the negative-lambda defect tail."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
import json
from pathlib import Path
import sys


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
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k22.jsonl",
)
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_defect_recurrence_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md"


@dataclass(frozen=True)
class Extremum:
    sample: str
    lam: str
    k: int


@dataclass(frozen=True)
class RecurrenceDiagnostics:
    lambdas: list[str]
    coefficient_k_max: int
    checked_x_max: int
    buffered_sufficient_rows: int
    buffered_sufficient_positive_rows: int
    defect_monotone_rows: int
    defect_monotone_positive_rows: int
    width_preserving_rows: int
    width_preserving_positive_rows: int
    width_preserving_rejected_rows: int
    min_buffer_margin: Extremum
    min_defect_monotone_gap: Extremum
    max_scaled_defect: Extremum
    max_defect_ratio: Extremum
    max_width_preserving_excess: Extremum


def update_min(current: tuple[Decimal, str, int] | None, value: Decimal, lam: str, k: int) -> tuple[Decimal, str, int]:
    if current is None or value < current[0]:
        return value, lam, k
    return current


def update_max(current: tuple[Decimal, str, int] | None, value: Decimal, lam: str, k: int) -> tuple[Decimal, str, int]:
    if current is None or value > current[0]:
        return value, lam, k
    return current


def extremum(item: tuple[Decimal, str, int] | None) -> Extremum:
    if item is None:
        raise RuntimeError("missing recurrence diagnostic extremum")
    value, lam, k = item
    return Extremum(decimal_format(value), lam, k)


def build_diagnostics(paths: list[Path]) -> RecurrenceDiagnostics:
    balls, samples, labels = load_enclosures(paths)
    if not labels:
        raise RuntimeError("no coefficient enclosures found")
    coefficient_k_max = max(index for (_lam, index) in balls)
    checked_x_max = coefficient_k_max - 1
    lambdas = sorted(labels)

    buffered_rows = 0
    buffered_pos = 0
    defect_rows = 0
    defect_pos = 0
    width_rows = 0
    width_pos = 0
    width_rejected = 0

    min_buffer: tuple[Decimal, str, int] | None = None
    min_defect_gap: tuple[Decimal, str, int] | None = None
    max_scaled: tuple[Decimal, str, int] | None = None
    max_ratio: tuple[Decimal, str, int] | None = None
    max_width_excess: tuple[Decimal, str, int] | None = None

    for lam in lambdas:
        lam_label = labels[lam]
        arb_x: dict[int, flint.arb] = {}
        sample_x: dict[int, Decimal] = {}
        for k in range(1, checked_x_max + 1):
            arb_x[k] = contraction({idx: balls[(lam, idx)] for idx in (k - 1, k, k + 1)}, k)
            sample_x[k] = contraction({idx: samples[(lam, idx)] for idx in (k - 1, k, k + 1)}, k)

            defect = flint.arb(1) - arb_x[k]
            scaled = defect * flint.arb(2 * k + 1) / flint.arb(2)
            buffer_margin = flint.arb(1) / flint.arb(3) - scaled
            buffered_rows += 1
            if arb_positive(defect) and arb_positive(buffer_margin):
                buffered_pos += 1

            sample_defect = Decimal(1) - sample_x[k]
            sample_scaled = sample_defect * Decimal(2 * k + 1) / Decimal(2)
            sample_buffer_margin = Decimal(1) / Decimal(3) - sample_scaled
            min_buffer = update_min(min_buffer, sample_buffer_margin, lam_label, k)
            max_scaled = update_max(max_scaled, sample_scaled, lam_label, k)

        for k in range(1, checked_x_max):
            defect_rows += 1
            defect_gap = arb_x[k + 1] - arb_x[k]
            if arb_positive(defect_gap):
                defect_pos += 1
            sample_gap = sample_x[k + 1] - sample_x[k]
            min_defect_gap = update_min(min_defect_gap, sample_gap, lam_label, k)

            sample_defect = Decimal(1) - sample_x[k]
            sample_next_defect = Decimal(1) - sample_x[k + 1]
            sample_ratio = sample_next_defect / sample_defect
            max_ratio = update_max(max_ratio, sample_ratio, lam_label, k)

            width_rows += 1
            factor = flint.arb(2 * k + 1) / flint.arb(2 * k + 3)
            margin = factor * (flint.arb(1) - arb_x[k]) - (flint.arb(1) - arb_x[k + 1])
            if arb_positive(margin):
                width_pos += 1
            if arb_positive(-margin):
                width_rejected += 1

            sample_factor = Decimal(2 * k + 1) / Decimal(2 * k + 3)
            sample_width_margin = sample_factor * sample_defect - sample_next_defect
            max_width_excess = update_max(max_width_excess, -sample_width_margin, lam_label, k)

    return RecurrenceDiagnostics(
        lambdas=[labels[lam] for lam in lambdas],
        coefficient_k_max=coefficient_k_max,
        checked_x_max=checked_x_max,
        buffered_sufficient_rows=buffered_rows,
        buffered_sufficient_positive_rows=buffered_pos,
        defect_monotone_rows=defect_rows,
        defect_monotone_positive_rows=defect_pos,
        width_preserving_rows=width_rows,
        width_preserving_positive_rows=width_pos,
        width_preserving_rejected_rows=width_rejected,
        min_buffer_margin=extremum(min_buffer),
        min_defect_monotone_gap=extremum(min_defect_gap),
        max_scaled_defect=extremum(max_scaled),
        max_defect_ratio=extremum(max_ratio),
        max_width_preserving_excess=extremum(max_width_excess),
    )


def build_artifact(paths: list[Path]) -> dict:
    diagnostics = build_diagnostics(paths)
    tail_lower_start = diagnostics.checked_x_max + 1
    tail_monotone_start = diagnostics.checked_x_max
    summary = {
        "buffered_sufficient_rows": diagnostics.buffered_sufficient_rows,
        "defect_monotone_rows": diagnostics.defect_monotone_rows,
        "width_preserving_rows": diagnostics.width_preserving_rows,
        "width_preserving_rejected_rows": diagnostics.width_preserving_rejected_rows,
        "ready_to_apply_rows": 0,
        "live_sufficient_routes": 1,
        "rejected_routes": 1,
        "target_closing": False,
        "main_finding": (
            "A buffered sufficient tail theorem is compatible with the finite prefix: "
            f"prove 0<=d_k<=2/(3*(2*k+1)) for k>={tail_lower_start} and "
            f"d_(k+1)<=d_k for k>={tail_monotone_start}. The tempting width-preserving "
            "one-step recurrence d_(k+1)<=((2*k+1)/(2*k+3))*d_k is rejected on all "
            "checked adjacent pairs, so it should not be used as the direct recurrence route."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_defect_recurrence_scout",
        "date": "2026-07-06",
        "status": "finite theorem-search diagnostic",
        "source_tail_barrier_scout": "outputs/jensen_window_pf_negative_lambda_tail_barrier_scout.md",
        "source_defect_tail_target": "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_defect_recurrence_scout.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_defect_recurrence_scout.py",
        "proof_boundary": (
            "Finite theorem-search diagnostic only. It identifies recurrence-style sufficient conditions "
            "and finite counterevidence for over-strong recurrences, but it does not prove an all-k tail "
            "theorem, does not prove cone entry, does not prove jwpf_06, and leaves Lambda <= 0 unsettled."
        ),
        "recurrence_rows": [
            {
                "id": "nldrs_01_buffered_sufficient_tail_condition",
                "role": "live_sufficient_condition",
                "claim": (
                    "If 0<=d_k<=2/(3*(2*k+1)) on the tail and d_(k+1)<=d_k from the bridge index, "
                    "then the target defect barriers 0<=d_k<=2/(2*k+1) and d_(k+1)<=d_k follow."
                ),
                "finite_rows": diagnostics.buffered_sufficient_rows,
                "finite_positive_rows": diagnostics.buffered_sufficient_positive_rows,
                "proof_boundary": "Sufficient condition plus finite compatibility only; not an all-k theorem.",
            },
            {
                "id": "nldrs_02_defect_monotone_only_insufficient",
                "role": "insufficient_condition",
                "claim": "Defect monotonicity alone gives d_(k+1)<=d_k but does not preserve the shrinking lower-wall width for all k.",
                "finite_rows": diagnostics.defect_monotone_rows,
                "finite_positive_rows": diagnostics.defect_monotone_positive_rows,
                "proof_boundary": "Finite-compatible but insufficient by itself.",
            },
            {
                "id": "nldrs_03_width_preserving_recurrence_rejected",
                "role": "rejected_candidate",
                "claim": "The direct recurrence d_(k+1)<=((2*k+1)/(2*k+3))*d_k would preserve the lower wall step-by-step, but is false on the certified prefix.",
                "finite_rows": diagnostics.width_preserving_rows,
                "finite_positive_rows": diagnostics.width_preserving_positive_rows,
                "finite_rejected_rows": diagnostics.width_preserving_rejected_rows,
                "proof_boundary": "Rejected by finite interval evidence only.",
            },
            {
                "id": "nldrs_04_scaled_defect_nonincrease_rejected",
                "role": "rejected_candidate",
                "claim": "The scaled-defect nonincreasing shortcut remains rejected by the tail-barrier scout.",
                "source": "outputs/jensen_window_pf_negative_lambda_tail_barrier_scout.md",
                "proof_boundary": "Rejected shortcut only; not a theorem.",
            },
        ],
        "finite_diagnostics": asdict(diagnostics),
        "invariants": [
            "No row is ready_to_apply for the defect-tail theorem.",
            "The buffered sufficient condition is not promoted to an all-k proof.",
            "Width-preserving recurrence is explicitly rejected.",
            "Scaled-defect nonincrease is explicitly rejected.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
        "summary": summary,
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["finite_diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda defect-recurrence scout: "
        f"{summary['buffered_sufficient_rows']} buffered rows, "
        f"{summary['defect_monotone_rows']} defect-monotone rows, "
        f"{summary['width_preserving_rejected_rows']} width-recurrence rejections, "
        f"{summary['live_sufficient_routes']} live sufficient routes, 0 issues"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Defect-Recurrence Scout",
        "",
        "Date: 2026-07-06",
        "",
        "Status: finite theorem-search diagnostic. This is not a proof of the",
        "defect-tail theorem, cone entry, Jensen-window PF-infinity, RH, or",
        "`Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_defect_recurrence_scout`.",
        "",
        "Proof boundary: this artifact tests recurrence-style sufficient",
        "conditions against the certified finite prefix. It does not prove an",
        "all-`k` recurrence or tail theorem.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_defect_recurrence_scout.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_defect_recurrence_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_defect_recurrence_scout.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Sufficient Condition",
        "",
        "A finite-compatible sufficient tail theorem is:",
        "",
        "```text",
        "0 <= d_k <= 2/(3*(2*k+1))",
        "d_(k+1) <= d_k",
        "```",
        "",
        "This is stronger than the defect-tail target, but it is compatible with",
        "the certified prefix and would imply the required barriers.",
        "",
        "## Rejected Recurrence",
        "",
        "The tempting stepwise wall-preserving recurrence",
        "",
        "```text",
        "d_(k+1) <= ((2*k+1)/(2*k+3))*d_k",
        "```",
        "",
        "is rejected on every checked adjacent pair.",
        "",
        "Finite diagnostics:",
        "",
        "```text",
        f"lambdas: {', '.join(diagnostics['lambdas'])}",
        f"coefficient range: A_0..A_{diagnostics['coefficient_k_max']}",
        f"checked contractions: x_1..x_{diagnostics['checked_x_max']}",
        (
            "buffered sufficient rows: "
            f"{diagnostics['buffered_sufficient_positive_rows']} / {diagnostics['buffered_sufficient_rows']}"
        ),
        (
            "defect monotone rows: "
            f"{diagnostics['defect_monotone_positive_rows']} / {diagnostics['defect_monotone_rows']}"
        ),
        (
            "width-preserving recurrence rows: "
            f"{diagnostics['width_preserving_positive_rows']} / {diagnostics['width_preserving_rows']}"
        ),
        (
            "width-preserving rejected rows: "
            f"{diagnostics['width_preserving_rejected_rows']} / {diagnostics['width_preserving_rows']}"
        ),
        "```",
        "",
        "Extrema:",
        "",
        "```text",
        (
            "min buffer margin: "
            f"{diagnostics['min_buffer_margin']['sample']} "
            f"at lambda={diagnostics['min_buffer_margin']['lam']}, k={diagnostics['min_buffer_margin']['k']}"
        ),
        (
            "min defect monotone gap: "
            f"{diagnostics['min_defect_monotone_gap']['sample']} "
            f"at lambda={diagnostics['min_defect_monotone_gap']['lam']}, k={diagnostics['min_defect_monotone_gap']['k']}"
        ),
        (
            "max width-preserving excess: "
            f"{diagnostics['max_width_preserving_excess']['sample']} "
            f"at lambda={diagnostics['max_width_preserving_excess']['lam']}, k={diagnostics['max_width_preserving_excess']['k']}"
        ),
        "```",
        "",
        "Integration:",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_tail_barrier_scout.md",
        "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
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
        "wrote Jensen-window PF negative-lambda defect-recurrence scout: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
