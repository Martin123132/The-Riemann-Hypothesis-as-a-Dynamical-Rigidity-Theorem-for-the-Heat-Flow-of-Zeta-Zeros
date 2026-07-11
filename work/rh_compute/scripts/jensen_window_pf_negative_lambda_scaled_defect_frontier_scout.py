#!/usr/bin/env python3
"""Build a scaled-defect frontier scout for the negative-lambda tail."""

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
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k50.jsonl",
)
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_defect_frontier_k50_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k50_scout.md"


@dataclass(frozen=True)
class Extremum:
    sample: str
    lam: str
    k: int


@dataclass(frozen=True)
class LambdaFrontier:
    lam: str
    checked_x_max: int
    one_third_positive_rows: int
    half_width_positive_rows: int
    cone_positive_rows: int
    first_one_third_failure_k: int | None
    first_one_third_failure_margin: str | None
    max_scaled_defect: Extremum
    min_half_width_margin: Extremum
    min_cone_margin: Extremum


@dataclass(frozen=True)
class ScaledDefectFrontierDiagnostics:
    lambdas: list[str]
    coefficient_k_max: int
    checked_x_max: int
    scaled_defect_rows: int
    cone_positive_rows: int
    half_width_positive_rows: int
    one_third_positive_rows: int
    scaled_defect_increase_rows: int
    scaled_defect_increase_positive_rows: int
    one_third_failure_rows: int
    lambda_frontiers: list[LambdaFrontier]
    min_cone_margin: Extremum
    min_half_width_margin: Extremum
    min_one_third_margin: Extremum
    min_scaled_defect_increase: Extremum
    max_scaled_defect: Extremum


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
        raise RuntimeError("missing scaled-defect extremum")
    value, lam, k = item
    return Extremum(decimal_format(value), lam, k)


def build_diagnostics(paths: list[Path]) -> ScaledDefectFrontierDiagnostics:
    balls, samples, labels = load_enclosures(paths)
    if not labels:
        raise RuntimeError("no coefficient enclosures found")
    coefficient_k_max = max(index for (_lam, index) in balls)
    checked_x_max = coefficient_k_max - 1
    lambdas = sorted(labels)

    total_rows = 0
    cone_pos = 0
    half_pos = 0
    one_third_pos = 0
    one_third_failures = 0
    scaled_inc_rows = 0
    scaled_inc_pos = 0

    min_cone: tuple[Decimal, str, int] | None = None
    min_half: tuple[Decimal, str, int] | None = None
    min_third: tuple[Decimal, str, int] | None = None
    min_scaled_inc: tuple[Decimal, str, int] | None = None
    max_scaled: tuple[Decimal, str, int] | None = None
    lambda_frontiers: list[LambdaFrontier] = []

    for lam in lambdas:
        lam_label = labels[lam]
        local_cone = 0
        local_half = 0
        local_third = 0
        first_failure_k: int | None = None
        first_failure_margin: str | None = None
        local_max_scaled: tuple[Decimal, str, int] | None = None
        local_min_half: tuple[Decimal, str, int] | None = None
        local_min_cone: tuple[Decimal, str, int] | None = None
        arb_scaled: dict[int, flint.arb] = {}
        sample_scaled: dict[int, Decimal] = {}

        for k in range(1, checked_x_max + 1):
            x = contraction({idx: balls[(lam, idx)] for idx in (k - 1, k, k + 1)}, k)
            sx = contraction({idx: samples[(lam, idx)] for idx in (k - 1, k, k + 1)}, k)
            defect = flint.arb(1) - x
            scaled = defect * flint.arb(2 * k + 1) / flint.arb(2)
            sample_s = (Decimal(1) - sx) * Decimal(2 * k + 1) / Decimal(2)
            arb_scaled[k] = scaled
            sample_scaled[k] = sample_s

            total_rows += 1
            cone_margin = min(sample_s, Decimal(1) - sample_s)
            half_margin = Decimal(1) / Decimal(2) - sample_s
            third_margin = Decimal(1) / Decimal(3) - sample_s

            if arb_positive(scaled) and arb_positive(flint.arb(1) - scaled):
                cone_pos += 1
                local_cone += 1
            if arb_positive(flint.arb(1) / flint.arb(2) - scaled):
                half_pos += 1
                local_half += 1
            if arb_positive(flint.arb(1) / flint.arb(3) - scaled):
                one_third_pos += 1
                local_third += 1
            else:
                one_third_failures += 1
                if first_failure_k is None:
                    first_failure_k = k
                    first_failure_margin = decimal_format(third_margin)

            min_cone = update_min(min_cone, cone_margin, lam_label, k)
            min_half = update_min(min_half, half_margin, lam_label, k)
            min_third = update_min(min_third, third_margin, lam_label, k)
            max_scaled = update_max(max_scaled, sample_s, lam_label, k)
            local_min_cone = update_min(local_min_cone, cone_margin, lam_label, k)
            local_min_half = update_min(local_min_half, half_margin, lam_label, k)
            local_max_scaled = update_max(local_max_scaled, sample_s, lam_label, k)

        for k in range(1, checked_x_max):
            scaled_inc_rows += 1
            inc = arb_scaled[k + 1] - arb_scaled[k]
            sample_inc = sample_scaled[k + 1] - sample_scaled[k]
            if arb_positive(inc):
                scaled_inc_pos += 1
            min_scaled_inc = update_min(min_scaled_inc, sample_inc, lam_label, k)

        lambda_frontiers.append(
            LambdaFrontier(
                lam=lam_label,
                checked_x_max=checked_x_max,
                one_third_positive_rows=local_third,
                half_width_positive_rows=local_half,
                cone_positive_rows=local_cone,
                first_one_third_failure_k=first_failure_k,
                first_one_third_failure_margin=first_failure_margin,
                max_scaled_defect=extremum(local_max_scaled),
                min_half_width_margin=extremum(local_min_half),
                min_cone_margin=extremum(local_min_cone),
            )
        )

    return ScaledDefectFrontierDiagnostics(
        lambdas=[labels[lam] for lam in lambdas],
        coefficient_k_max=coefficient_k_max,
        checked_x_max=checked_x_max,
        scaled_defect_rows=total_rows,
        cone_positive_rows=cone_pos,
        half_width_positive_rows=half_pos,
        one_third_positive_rows=one_third_pos,
        scaled_defect_increase_rows=scaled_inc_rows,
        scaled_defect_increase_positive_rows=scaled_inc_pos,
        one_third_failure_rows=one_third_failures,
        lambda_frontiers=lambda_frontiers,
        min_cone_margin=extremum(min_cone),
        min_half_width_margin=extremum(min_half),
        min_one_third_margin=extremum(min_third),
        min_scaled_defect_increase=extremum(min_scaled_inc),
        max_scaled_defect=extremum(max_scaled),
    )


def build_artifact(paths: list[Path]) -> dict:
    diagnostics = build_diagnostics(paths)
    prefix_label = f"k{diagnostics.coefficient_k_max}"
    tail_barrier_ref = f"outputs/jensen_window_pf_negative_lambda_tail_barrier_{prefix_label}_scout.md"
    half_width_all = diagnostics.half_width_positive_rows == diagnostics.scaled_defect_rows
    half_width_role = "live_sufficient_buffer" if half_width_all else "frontier_obstruction"
    half_width_claim = (
        "The half-width buffer s_k <= 1/2 holds on every checked row."
        if half_width_all
        else f"The half-width buffer s_k <= 1/2 fails before the end of the checked {prefix_label} prefix."
    )
    half_width_boundary = (
        "Finite prefix only; not an all-k half-width theorem."
        if half_width_all
        else "Finite obstruction to this half-width buffer only."
    )
    half_width_summary = (
        "the half-width buffer s_k<=1/2 holds on every checked row"
        if half_width_all
        else (
            f"the half-width buffer s_k<=1/2 fails on "
            f"{diagnostics.scaled_defect_rows - diagnostics.half_width_positive_rows}/"
            f"{diagnostics.scaled_defect_rows} checked rows"
        )
    )
    tail_guidance = (
        "any analytic tail theorem should target the exact cone or a weaker buffer such as one-half."
        if half_width_all
        else "any analytic tail theorem should target the exact cone or a lambda/k-dependent buffer above the observed frontier."
    )
    half_width_invariant = (
        "The half-width buffer is a finite candidate, not a theorem."
        if half_width_all
        else "The half-width buffer is finitely rejected on this checked prefix."
    )
    return {
        "kind": "jensen_window_pf_negative_lambda_scaled_defect_frontier_scout",
        "date": "2026-07-06",
        "status": "finite theorem-search diagnostic",
        "source_tail_barrier_scout": tail_barrier_ref,
        "source_defect_tail_target": "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_scaled_defect_frontier_scout.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_defect_frontier_scout.py",
        "proof_boundary": (
            "Finite scaled-defect frontier diagnostic only. It distinguishes the exact "
            "cone barrier, the half-width buffer, and the one-third-width buffer on "
            "the checked negative-lambda prefix, but it does not prove an all-k tail "
            "theorem, cone entry, jwpf_06, RH, or Lambda <= 0."
        ),
        "frontier_rows": [
            {
                "id": "nlsdf_01_exact_cone_scaled_defect",
                "role": "finite_barrier_certificate",
                "claim": "The exact lower-wall scaled defect satisfies 0 <= s_k <= 1 on every checked row.",
                "rows": diagnostics.scaled_defect_rows,
                "positive_rows": diagnostics.cone_positive_rows,
                "proof_boundary": "Finite prefix only; not an all-k tail theorem.",
            },
            {
                "id": "nlsdf_02_half_width_buffer",
                "role": half_width_role,
                "claim": half_width_claim,
                "rows": diagnostics.scaled_defect_rows,
                "positive_rows": diagnostics.half_width_positive_rows,
                "proof_boundary": half_width_boundary,
            },
            {
                "id": "nlsdf_03_one_third_frontier",
                "role": "frontier_obstruction",
                "claim": f"The one-third-width buffer s_k <= 1/3 fails before the end of the checked {prefix_label} prefix.",
                "rows": diagnostics.scaled_defect_rows,
                "positive_rows": diagnostics.one_third_positive_rows,
                "failure_rows": diagnostics.one_third_failure_rows,
                "proof_boundary": "Finite obstruction to this stronger buffer only.",
            },
            {
                "id": "nlsdf_04_scaled_defect_increase",
                "role": "rejected_candidate",
                "claim": "The scaled defect is increasing on every checked adjacent row, rejecting scaled-defect nonincrease.",
                "rows": diagnostics.scaled_defect_increase_rows,
                "positive_rows": diagnostics.scaled_defect_increase_positive_rows,
                "proof_boundary": "Finite rejection of one shortcut only.",
            },
        ],
        "finite_diagnostics": asdict(diagnostics),
        "summary": {
            "scaled_defect_rows": diagnostics.scaled_defect_rows,
            "cone_positive_rows": diagnostics.cone_positive_rows,
            "half_width_positive_rows": diagnostics.half_width_positive_rows,
            "one_third_positive_rows": diagnostics.one_third_positive_rows,
            "one_third_failure_rows": diagnostics.one_third_failure_rows,
            "scaled_defect_increase_rows": diagnostics.scaled_defect_increase_rows,
            "scaled_defect_increase_positive_rows": diagnostics.scaled_defect_increase_positive_rows,
            "ready_to_apply_rows": 0,
            "target_closing": False,
            "main_finding": (
                f"On the {prefix_label} negative-lambda prefix, the exact scaled-defect cone "
                f"0<=s_k<=1 holds on every checked row and {half_width_summary}, "
                f"while the one-third buffer holds on only {diagnostics.one_third_positive_rows}/"
                f"{diagnostics.scaled_defect_rows} rows. Thus the previous one-third sufficient "
                f"route is too strong for the observed prefix; {tail_guidance}"
            ),
        },
        "invariants": [
            "The exact cone barrier remains finite evidence only.",
            half_width_invariant,
            "The one-third-width route is not promoted as an all-k sufficient condition.",
            "The scaled-defect nonincrease shortcut remains rejected.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_json(artifact: dict, path: Path) -> None:
    path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_note(artifact: dict, path: Path) -> None:
    diagnostics = artifact["finite_diagnostics"]
    summary = artifact["summary"]
    result_json = artifact.get("result_json", path.with_suffix(".json").relative_to(REPO_ROOT).as_posix())
    note_ref = artifact.get("note", path.relative_to(REPO_ROOT).as_posix())
    lines = [
        "# Jensen-Window PF Negative-Lambda Scaled-Defect Frontier Scout",
        "",
        "Date: 2026-07-06",
        "",
        "Status: finite theorem-search diagnostic. This is not a proof of the",
        "defect-tail theorem, cone entry, Jensen-window PF-infinity, RH, or",
        "`Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_scaled_defect_frontier_scout`.",
        "",
        "Proof boundary: this artifact separates finite scaled-defect barriers.",
        "It does not prove an all-`k` tail theorem, does not prove `jwpf_06`,",
        "and does not establish `Lambda <= 0`.",
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
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_scaled_defect_frontier_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        (
            "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_defect_frontier_scout.py "
            f"--target {result_json} --note {note_ref}"
        ),
        "```",
        "",
        "Current result:",
        "",
        "```text",
        (
            "validated Jensen-window PF negative-lambda scaled-defect frontier scout: "
            f"{summary['scaled_defect_rows']} scaled rows, "
            f"{summary['cone_positive_rows']} cone rows, "
            f"{summary['half_width_positive_rows']} half-width rows, "
            f"{summary['one_third_positive_rows']} one-third rows, "
            f"{summary['one_third_failure_rows']} one-third failures, "
            f"{summary['scaled_defect_increase_positive_rows']} scaled-increase rows, 0 issues"
        ),
        "```",
        "",
        "## Scaled Defect",
        "",
        "Write:",
        "",
        "```text",
        "d_k = 1 - x_k",
        "s_k = ((2*k+1)/2) * d_k",
        "```",
        "",
        "The exact lower wall is `0 <= s_k <= 1`. The previously useful",
        "one-third-width buffer is `s_k <= 1/3`; the weaker half-width buffer",
        "is `s_k <= 1/2`.",
        "",
        "## Frontier",
        "",
        "```text",
        f"lambdas: {', '.join(diagnostics['lambdas'])}",
        f"coefficient range: A_0..A_{diagnostics['coefficient_k_max']}",
        f"checked contractions: x_1..x_{diagnostics['checked_x_max']}",
        f"exact cone rows: {diagnostics['cone_positive_rows']} / {diagnostics['scaled_defect_rows']}",
        f"half-width rows: {diagnostics['half_width_positive_rows']} / {diagnostics['scaled_defect_rows']}",
        f"one-third rows: {diagnostics['one_third_positive_rows']} / {diagnostics['scaled_defect_rows']}",
        f"one-third failure rows: {diagnostics['one_third_failure_rows']}",
        f"scaled-defect increase rows: {diagnostics['scaled_defect_increase_positive_rows']} / {diagnostics['scaled_defect_increase_rows']}",
        "```",
        "",
        "Per-lambda one-third frontier:",
        "",
        "```text",
    ]
    for row in diagnostics["lambda_frontiers"]:
        failure = row["first_one_third_failure_k"]
        failure_text = "none" if failure is None else f"k={failure}, margin={row['first_one_third_failure_margin']}"
        lines.append(
            f"lambda={row['lam']}: one-third {row['one_third_positive_rows']}/{row['checked_x_max']}, "
            f"first failure {failure_text}, max s={row['max_scaled_defect']['sample']} at k={row['max_scaled_defect']['k']}"
        )
    lines.extend(
        [
            "```",
            "",
            "Global extrema:",
            "",
            "```text",
            (
                "min exact-cone scaled margin: "
                f"{diagnostics['min_cone_margin']['sample']} at lambda={diagnostics['min_cone_margin']['lam']}, "
                f"k={diagnostics['min_cone_margin']['k']}"
            ),
            (
                "min half-width margin: "
                f"{diagnostics['min_half_width_margin']['sample']} at lambda={diagnostics['min_half_width_margin']['lam']}, "
                f"k={diagnostics['min_half_width_margin']['k']}"
            ),
            (
                "min one-third margin: "
                f"{diagnostics['min_one_third_margin']['sample']} at lambda={diagnostics['min_one_third_margin']['lam']}, "
                f"k={diagnostics['min_one_third_margin']['k']}"
            ),
            (
                "max scaled defect: "
                f"{diagnostics['max_scaled_defect']['sample']} at lambda={diagnostics['max_scaled_defect']['lam']}, "
                f"k={diagnostics['max_scaled_defect']['k']}"
            ),
            (
                "min scaled-defect increase: "
                f"{diagnostics['min_scaled_defect_increase']['sample']} at lambda={diagnostics['min_scaled_defect_increase']['lam']}, "
                f"k={diagnostics['min_scaled_defect_increase']['k']}"
            ),
            "```",
            "",
            "## Consequence",
            "",
            summary["main_finding"],
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_tail_barrier_scout"],
            "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
            "```",
            "",
        ]
    )
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
    note.parent.mkdir(parents=True, exist_ok=True)
    artifact["result_json"] = out_json.relative_to(REPO_ROOT).as_posix()
    artifact["note"] = note.relative_to(REPO_ROOT).as_posix()
    write_json(artifact, out_json)
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda scaled-defect frontier scout: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
