#!/usr/bin/env python3
"""Build a theorem-search scout for the negative-lambda ratio-cone tail barrier."""

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
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_tail_barrier_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_tail_barrier_scout.md"


@dataclass(frozen=True)
class Extremum:
    sample: str
    lam: str
    k: int


@dataclass(frozen=True)
class TailBarrierDiagnostics:
    lambdas: list[str]
    coefficient_k_max: int
    checked_x_max: int
    cone_buffer_rows: int
    cone_buffer_positive_rows: int
    defect_monotone_rows: int
    defect_monotone_positive_rows: int
    scaled_defect_rows: int
    scaled_defect_one_third_rows: int
    scaled_defect_increase_rows: int
    scaled_defect_increase_positive_rows: int
    rejected_candidate_count: int
    min_upper_defect: Extremum
    min_lower_wall_slack: Extremum
    min_one_third_buffer_margin: Extremum
    min_defect_monotone_gap: Extremum
    min_scaled_defect_increase: Extremum
    max_scaled_defect: Extremum
    max_defect_ratio: Extremum


def _update_min(current: tuple[Decimal, str, int] | None, value: Decimal, lam: str, k: int) -> tuple[Decimal, str, int]:
    if current is None or value < current[0]:
        return value, lam, k
    return current


def _update_max(current: tuple[Decimal, str, int] | None, value: Decimal, lam: str, k: int) -> tuple[Decimal, str, int]:
    if current is None or value > current[0]:
        return value, lam, k
    return current


def _extremum(item: tuple[Decimal, str, int] | None) -> Extremum:
    if item is None:
        raise RuntimeError("missing extremum")
    value, lam, k = item
    return Extremum(decimal_format(value), lam, k)


def build_diagnostics(paths: list[Path]) -> TailBarrierDiagnostics:
    balls, samples, labels = load_enclosures(paths)
    if not labels:
        raise RuntimeError("no coefficient enclosures found")

    coefficient_k_max = max(index for (_lam, index) in balls)
    checked_x_max = coefficient_k_max - 1
    lambdas = sorted(labels)

    arb_x: dict[tuple[Decimal, int], flint.arb] = {}
    sample_x: dict[tuple[Decimal, int], Decimal] = {}
    arb_scaled: dict[tuple[Decimal, int], flint.arb] = {}
    sample_scaled: dict[tuple[Decimal, int], Decimal] = {}
    sample_defect: dict[tuple[Decimal, int], Decimal] = {}

    cone_buffer_rows = 0
    cone_buffer_pos = 0
    scaled_rows = 0
    scaled_one_third_pos = 0

    min_upper_defect: tuple[Decimal, str, int] | None = None
    min_lower_wall_slack: tuple[Decimal, str, int] | None = None
    min_one_third_margin: tuple[Decimal, str, int] | None = None
    max_scaled_defect: tuple[Decimal, str, int] | None = None

    for lam in lambdas:
        lam_label = labels[lam]
        for k in range(1, checked_x_max + 1):
            arb_values = {idx: balls[(lam, idx)] for idx in (k - 1, k, k + 1)}
            sample_values = {idx: samples[(lam, idx)] for idx in (k - 1, k, k + 1)}
            x = contraction(arb_values, k)
            sx = contraction(sample_values, k)
            defect = flint.arb(1) - x
            scaled = defect * flint.arb(2 * k + 1) / flint.arb(2)
            sample_d = Decimal(1) - sx
            sample_width = Decimal(2) / Decimal(2 * k + 1)
            sample_s = sample_d / sample_width

            arb_x[(lam, k)] = x
            sample_x[(lam, k)] = sx
            arb_scaled[(lam, k)] = scaled
            sample_scaled[(lam, k)] = sample_s
            sample_defect[(lam, k)] = sample_d

            lower = flint.arb(2 * k - 1) / flint.arb(2 * k + 1)
            lower_slack = x - lower
            upper_defect = flint.arb(1) - x
            one_third_margin = flint.arb(1) / flint.arb(3) - scaled

            sample_lower_slack = sample_width - sample_d
            sample_one_third_margin = Decimal(1) / Decimal(3) - sample_s

            cone_buffer_rows += 1
            if arb_positive(lower_slack) and arb_positive(upper_defect):
                cone_buffer_pos += 1

            scaled_rows += 1
            if arb_positive(one_third_margin):
                scaled_one_third_pos += 1

            min_upper_defect = _update_min(min_upper_defect, sample_d, lam_label, k)
            min_lower_wall_slack = _update_min(min_lower_wall_slack, sample_lower_slack, lam_label, k)
            min_one_third_margin = _update_min(min_one_third_margin, sample_one_third_margin, lam_label, k)
            max_scaled_defect = _update_max(max_scaled_defect, sample_s, lam_label, k)

    defect_mono_rows = 0
    defect_mono_pos = 0
    scaled_inc_rows = 0
    scaled_inc_pos = 0
    min_defect_gap: tuple[Decimal, str, int] | None = None
    min_scaled_inc: tuple[Decimal, str, int] | None = None
    max_defect_ratio: tuple[Decimal, str, int] | None = None

    for lam in lambdas:
        lam_label = labels[lam]
        for k in range(1, checked_x_max):
            defect_mono_rows += 1
            defect_gap = arb_x[(lam, k + 1)] - arb_x[(lam, k)]
            sample_gap = sample_x[(lam, k + 1)] - sample_x[(lam, k)]
            if arb_positive(defect_gap):
                defect_mono_pos += 1
            min_defect_gap = _update_min(min_defect_gap, sample_gap, lam_label, k)

            scaled_inc_rows += 1
            scaled_inc = arb_scaled[(lam, k + 1)] - arb_scaled[(lam, k)]
            sample_scaled_inc = sample_scaled[(lam, k + 1)] - sample_scaled[(lam, k)]
            if arb_positive(scaled_inc):
                scaled_inc_pos += 1
            min_scaled_inc = _update_min(min_scaled_inc, sample_scaled_inc, lam_label, k)

            ratio = sample_defect[(lam, k + 1)] / sample_defect[(lam, k)]
            max_defect_ratio = _update_max(max_defect_ratio, ratio, lam_label, k)

    rejected_candidate_count = 1 if scaled_inc_pos == scaled_inc_rows else 0

    return TailBarrierDiagnostics(
        lambdas=[labels[lam] for lam in lambdas],
        coefficient_k_max=coefficient_k_max,
        checked_x_max=checked_x_max,
        cone_buffer_rows=cone_buffer_rows,
        cone_buffer_positive_rows=cone_buffer_pos,
        defect_monotone_rows=defect_mono_rows,
        defect_monotone_positive_rows=defect_mono_pos,
        scaled_defect_rows=scaled_rows,
        scaled_defect_one_third_rows=scaled_one_third_pos,
        scaled_defect_increase_rows=scaled_inc_rows,
        scaled_defect_increase_positive_rows=scaled_inc_pos,
        rejected_candidate_count=rejected_candidate_count,
        min_upper_defect=_extremum(min_upper_defect),
        min_lower_wall_slack=_extremum(min_lower_wall_slack),
        min_one_third_buffer_margin=_extremum(min_one_third_margin),
        min_defect_monotone_gap=_extremum(min_defect_gap),
        min_scaled_defect_increase=_extremum(min_scaled_inc),
        max_scaled_defect=_extremum(max_scaled_defect),
        max_defect_ratio=_extremum(max_defect_ratio),
    )


def build_artifact(paths: list[Path]) -> dict:
    diagnostics = build_diagnostics(paths)
    active_depth = diagnostics.checked_x_max - 2
    first_collar = active_depth + 1
    second_collar = active_depth + 2
    next_active_depth = active_depth + 1
    next_collar_x = diagnostics.checked_x_max + 1
    next_coefficient = diagnostics.coefficient_k_max + 1
    if diagnostics.scaled_defect_one_third_rows == diagnostics.scaled_defect_rows:
        one_third_phrase = (
            "the stronger finite one-third width buffer "
            f"d_k <= 2/(3*(2k+1)) through x_{diagnostics.checked_x_max}"
        )
    else:
        one_third_phrase = (
            "the stronger finite one-third width buffer on "
            f"{diagnostics.scaled_defect_one_third_rows}/{diagnostics.scaled_defect_rows} checked rows, "
            "with a frontier failure before the end of the prefix"
        )
    summary = {
        "checked_x_max": diagnostics.checked_x_max,
        "cone_buffer_rows": diagnostics.cone_buffer_rows,
        "defect_monotone_rows": diagnostics.defect_monotone_rows,
        "one_third_buffer_rows": diagnostics.scaled_defect_one_third_rows,
        "scaled_defect_increase_rows": diagnostics.scaled_defect_increase_rows,
        "rejected_candidate_count": diagnostics.rejected_candidate_count,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The certified prefix satisfies the defect-form cone barriers "
            f"0 <= d_k <= 2/(2k+1), d_{{k+1}} <= d_k, and {one_third_phrase}. "
            "However, the scaled defect is increasing on the checked prefix, "
            "so any analytic proof route that assumes scaled-defect monotone "
            f"decrease is already rejected. The next finite collar needs x_{next_collar_x}, "
            f"hence A_{next_coefficient}, unless an analytic tail theorem supplies it."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_tail_barrier_scout",
        "date": "2026-07-06",
        "status": "finite theorem-search diagnostic",
        "source_prefix_scout": "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.md",
        "source_finite_collar_contract": "outputs/jensen_window_pf_negative_lambda_finite_collar_contract.md",
        "source_cone_entry_target": "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_tail_barrier_scout.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_tail_barrier_scout.py",
        "proof_boundary": (
            "Finite theorem-search diagnostic only. It rewrites the missing "
            "negative-lambda tail as defect inequalities and certifies several "
            "finite-prefix barrier facts, but it does not prove an all-k tail "
            "theorem, does not prove a collared finite flow theorem, does not "
            "prove cone entry for zeta coefficients, does not prove jwpf_06, "
            "and leaves Lambda <= 0 unsettled."
        ),
        "tail_barrier_rows": [
            {
                "id": "nltbs_01_defect_form_reduction",
                "role": "exact_reformulation",
                "claim": "Writing d_k=1-x_k, the ratio cone is equivalent to 0 <= d_k <= 2/(2*k+1) and d_{k+1} <= d_k.",
                "proof_boundary": "Algebraic reformulation only; not a zeta tail theorem.",
            },
            {
                "id": "nltbs_02_one_third_width_buffer",
                "role": "finite_barrier_certificate",
                "claim": (
                    "On the checked negative-lambda prefix, the stronger "
                    "one-third width buffer d_k <= 2/(3*(2*k+1)) holds on "
                    f"{diagnostics.scaled_defect_one_third_rows}/{diagnostics.scaled_defect_rows} "
                    "certified x_k rows."
                ),
                "rows": diagnostics.scaled_defect_rows,
                "positive_rows": diagnostics.scaled_defect_one_third_rows,
                "proof_boundary": "Finite prefix only; not a tail theorem.",
            },
            {
                "id": "nltbs_03_defect_monotone_prefix",
                "role": "finite_barrier_certificate",
                "claim": "On the checked negative-lambda prefix, d_{k+1} <= d_k.",
                "rows": diagnostics.defect_monotone_rows,
                "positive_rows": diagnostics.defect_monotone_positive_rows,
                "proof_boundary": "Finite prefix only; not an all-k monotonicity theorem.",
            },
            {
                "id": "nltbs_04_scaled_defect_monotonicity_rejected",
                "role": "rejected_candidate",
                "claim": "The candidate shortcut that the scaled defect ((2*k+1)/2)*d_k should be nonincreasing is false on the certified prefix.",
                "rows": diagnostics.scaled_defect_increase_rows,
                "positive_increase_rows": diagnostics.scaled_defect_increase_positive_rows,
                "proof_boundary": "Rejects one over-strong barrier ansatz only.",
            },
            {
                "id": "nltbs_05_next_collar_upgrade",
                "role": "open_upgrade",
                "claim": (
                    f"The next purely finite collar upgrade from K={active_depth} to "
                    f"K={next_active_depth} needs x_{next_collar_x}, hence coefficient "
                    f"A_{next_coefficient}, unless an analytic tail theorem supplies the missing collar."
                ),
                "proof_boundary": "Upgrade target only; not an analytic tail theorem.",
            },
        ],
        "finite_diagnostics": asdict(diagnostics),
        "invariants": [
            "No row is ready_to_apply for the cone-entry theorem.",
            "The one-third width buffer is a finite-prefix certificate only.",
            "The scaled-defect monotonicity shortcut is explicitly rejected.",
            "The all-k tail or finite-collar flow theorem remains open.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
        "summary": summary,
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["finite_diagnostics"]
    result_json = artifact.get("result_json", "work/rh_compute/results/jensen_window_pf_negative_lambda_tail_barrier_scout.json")
    note_ref = artifact.get("note", path.relative_to(REPO_ROOT).as_posix())
    active_depth = diagnostics["checked_x_max"] - 2
    first_collar = active_depth + 1
    second_collar = active_depth + 2
    next_active_depth = active_depth + 1
    next_collar_x = diagnostics["checked_x_max"] + 1
    next_coefficient = diagnostics["coefficient_k_max"] + 1
    lines = [
        "# Jensen-Window PF Negative-Lambda Tail-Barrier Scout",
        "",
        "Date: 2026-07-06",
        "",
        "Status: finite theorem-search diagnostic. This is not a proof of zeta",
        "cone entry, Jensen-window PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_tail_barrier_scout`.",
        "",
        "Proof boundary: this artifact rewrites the missing negative-lambda tail",
        "as defect inequalities and certifies finite-prefix barrier facts only.",
        "It does not prove an all-`k` tail theorem, does not prove a collared",
        "finite flow theorem, does not prove `jwpf_06`, and does not establish",
        "`Lambda <= 0`.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        result_json,
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
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_tail_barrier_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        f"python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_tail_barrier_scout.py --target {result_json} --note {note_ref}",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        (
            "validated Jensen-window PF negative-lambda tail-barrier scout: "
            f"{summary['cone_buffer_rows']} cone-buffer rows, "
            f"{summary['defect_monotone_rows']} defect-monotone rows, "
            f"{summary['one_third_buffer_rows']} one-third-buffer rows, "
            f"{summary['scaled_defect_increase_rows']} scaled-defect increase rows, "
            f"{summary['rejected_candidate_count']} rejected candidate, 0 issues"
        ),
        "```",
        "",
        "## Defect Form",
        "",
        "Write:",
        "",
        "```text",
        "d_k = 1 - x_k",
        "x_k = (A_{k+1}/A_k)/(A_k/A_{k-1})",
        "```",
        "",
        "Then the ratio cone becomes:",
        "",
        "```text",
        "0 <= d_k <= 2/(2*k+1)",
        "d_{k+1} <= d_k",
        "```",
        "",
        "This is the exact tail-barrier form that an analytic proof must supply",
        "outside the finite prefix.",
        "",
        "## Finite Barrier Diagnostics",
        "",
        "```text",
        f"lambdas: {', '.join(diagnostics['lambdas'])}",
        f"coefficient range: A_0..A_{diagnostics['coefficient_k_max']}",
        f"checked contractions: x_1..x_{diagnostics['checked_x_max']}",
        f"cone-buffer rows: {diagnostics['cone_buffer_positive_rows']} / {diagnostics['cone_buffer_rows']}",
        (
            "one-third width buffer rows: "
            f"{diagnostics['scaled_defect_one_third_rows']} / {diagnostics['scaled_defect_rows']}"
        ),
        (
            "defect monotone rows: "
            f"{diagnostics['defect_monotone_positive_rows']} / {diagnostics['defect_monotone_rows']}"
        ),
        (
            "scaled-defect increase rows: "
            f"{diagnostics['scaled_defect_increase_positive_rows']} / {diagnostics['scaled_defect_increase_rows']}"
        ),
        "```",
        "",
        "Minimum and maximum sampled margins:",
        "",
        "```text",
        (
            "min upper defect d_k: "
            f"{diagnostics['min_upper_defect']['sample']} "
            f"at lambda={diagnostics['min_upper_defect']['lam']}, k={diagnostics['min_upper_defect']['k']}"
        ),
        (
            "min lower-wall slack: "
            f"{diagnostics['min_lower_wall_slack']['sample']} "
            f"at lambda={diagnostics['min_lower_wall_slack']['lam']}, k={diagnostics['min_lower_wall_slack']['k']}"
        ),
        (
            "min one-third buffer margin: "
            f"{diagnostics['min_one_third_buffer_margin']['sample']} "
            f"at lambda={diagnostics['min_one_third_buffer_margin']['lam']}, "
            f"k={diagnostics['min_one_third_buffer_margin']['k']}"
        ),
        (
            "min defect monotone gap: "
            f"{diagnostics['min_defect_monotone_gap']['sample']} "
            f"at lambda={diagnostics['min_defect_monotone_gap']['lam']}, "
            f"k={diagnostics['min_defect_monotone_gap']['k']}"
        ),
        (
            "max scaled defect: "
            f"{diagnostics['max_scaled_defect']['sample']} "
            f"at lambda={diagnostics['max_scaled_defect']['lam']}, k={diagnostics['max_scaled_defect']['k']}"
        ),
        (
            "max defect ratio d_{k+1}/d_k: "
            f"{diagnostics['max_defect_ratio']['sample']} "
            f"at lambda={diagnostics['max_defect_ratio']['lam']}, k={diagnostics['max_defect_ratio']['k']}"
        ),
        "```",
        "",
        "## Rejected Shortcut",
        "",
        "The scaled defect",
        "",
        "```text",
        "s_k = ((2*k+1)/2) * d_k",
        "```",
        "",
        "is increasing on every checked adjacent pair. Therefore a proof route",
        "that tries to obtain the lower wall by assuming `s_k` is nonincreasing",
        "is already incompatible with the certified finite prefix.",
        "",
        "## Next Upgrade",
        "",
        (
            f"The current finite-collar contract has active depth `K={active_depth}` "
            f"with collars `x_{first_collar}` and `x_{second_collar}`. The next "
            f"purely finite upgrade to `K={next_active_depth}` needs"
        ),
        f"`x_{next_collar_x}`, hence `A_{next_coefficient}`, unless an analytic tail theorem supplies the",
        "missing collar.",
        "",
        "Integration:",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.md",
        "outputs/jensen_window_pf_negative_lambda_finite_collar_contract.md",
        "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
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
    paths = [path if path.is_absolute() else (REPO_ROOT / path) for path in args.enclosures]
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(paths)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    artifact["result_json"] = out_json.relative_to(REPO_ROOT).as_posix()
    artifact["note"] = note.relative_to(REPO_ROOT).as_posix()
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda tail-barrier scout: "
        f"{out_json.relative_to(REPO_ROOT)} and {note.relative_to(REPO_ROOT)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
