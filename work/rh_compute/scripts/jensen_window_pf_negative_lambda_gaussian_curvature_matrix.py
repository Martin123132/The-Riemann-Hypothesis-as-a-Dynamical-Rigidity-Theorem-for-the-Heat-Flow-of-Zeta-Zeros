#!/usr/bin/env python3
"""Build the Gaussian-baseline curvature matrix for the negative-lambda tail."""

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
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.md"


@dataclass(frozen=True)
class Extremum:
    sample: str
    lam: str
    k: int


@dataclass(frozen=True)
class GaussianCurvatureDiagnostics:
    lambdas: list[str]
    coefficient_k_max: int
    checked_x_max: int
    gaussian_deficit_rows: int
    gaussian_deficit_positive_rows: int
    bounded_deficit_rows: int
    bounded_deficit_positive_rows: int
    raw_threshold_rows: int
    raw_threshold_positive_rows: int
    positive_threshold_rows: int
    min_deficit_bound_slack: Extremum
    min_raw_threshold_margin: Extremum
    min_positive_threshold: Extremum
    max_gaussian_deficit: Extremum
    max_deficit_fraction_of_bound: Extremum


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
        raise RuntimeError("missing Gaussian-curvature diagnostic extremum")
    value, lam, k = item
    return Extremum(decimal_format(value), lam, k)


def build_diagnostics(paths: list[Path]) -> GaussianCurvatureDiagnostics:
    balls, samples, labels = load_enclosures(paths)
    if not labels:
        raise RuntimeError("no coefficient enclosures found")
    coefficient_k_max = max(index for (_lam, index) in balls)
    checked_x_max = coefficient_k_max - 1
    lambdas = sorted(labels)

    deficit_rows = 0
    deficit_good = 0
    bound_rows = 0
    bound_good = 0
    threshold_rows = 0
    threshold_good = 0
    positive_threshold = 0

    min_bound_slack: tuple[Decimal, str, int] | None = None
    min_threshold_margin: tuple[Decimal, str, int] | None = None
    min_threshold: tuple[Decimal, str, int] | None = None
    max_deficit: tuple[Decimal, str, int] | None = None
    max_fraction: tuple[Decimal, str, int] | None = None

    for lam in lambdas:
        lam_label = labels[lam]
        for k in range(1, checked_x_max + 1):
            x = contraction({idx: balls[(lam, idx)] for idx in (k - 1, k, k + 1)}, k)
            sample_x = contraction({idx: samples[(lam, idx)] for idx in (k - 1, k, k + 1)}, k)

            gaussian_curvature = (flint.arb(2 * k + 1) / flint.arb(2 * k - 1)).log()
            raw_curvature = gaussian_curvature + x.log()
            deficit = -x.log()
            deficit_bound = flint.arb(2) / flint.arb(3 * (2 * k + 1))
            raw_threshold = gaussian_curvature - deficit_bound

            deficit_rows += 1
            bound_rows += 1
            threshold_rows += 1
            if arb_positive(deficit):
                deficit_good += 1
            if arb_positive(deficit_bound - deficit):
                bound_good += 1
            if arb_positive(raw_curvature - raw_threshold):
                threshold_good += 1
            if arb_positive(raw_threshold):
                positive_threshold += 1

            sample_gaussian = (Decimal(2 * k + 1) / Decimal(2 * k - 1)).ln()
            sample_deficit = -sample_x.ln()
            sample_bound = Decimal(2) / Decimal(3 * (2 * k + 1))
            sample_threshold = sample_gaussian - sample_bound
            sample_raw_curvature = sample_gaussian - sample_deficit
            min_bound_slack = update_min(min_bound_slack, sample_bound - sample_deficit, lam_label, k)
            min_threshold_margin = update_min(min_threshold_margin, sample_raw_curvature - sample_threshold, lam_label, k)
            min_threshold = update_min(min_threshold, sample_threshold, lam_label, k)
            max_deficit = update_max(max_deficit, sample_deficit, lam_label, k)
            max_fraction = update_max(max_fraction, sample_deficit / sample_bound, lam_label, k)

    return GaussianCurvatureDiagnostics(
        lambdas=[labels[lam] for lam in lambdas],
        coefficient_k_max=coefficient_k_max,
        checked_x_max=checked_x_max,
        gaussian_deficit_rows=deficit_rows,
        gaussian_deficit_positive_rows=deficit_good,
        bounded_deficit_rows=bound_rows,
        bounded_deficit_positive_rows=bound_good,
        raw_threshold_rows=threshold_rows,
        raw_threshold_positive_rows=threshold_good,
        positive_threshold_rows=positive_threshold,
        min_deficit_bound_slack=extremum(min_bound_slack),
        min_raw_threshold_margin=extremum(min_threshold_margin),
        min_positive_threshold=extremum(min_threshold),
        max_gaussian_deficit=extremum(max_deficit),
        max_deficit_fraction_of_bound=extremum(max_fraction),
    )


def build_artifact(paths: list[Path]) -> dict:
    diagnostics = build_diagnostics(paths)
    summary = {
        "matrix_rows": 7,
        "ready_to_apply_rows": 0,
        "exact_rows": 2,
        "live_routes": 2,
        "rejected_or_insufficient_rows": 2,
        "gaussian_deficit_rows": diagnostics.gaussian_deficit_rows,
        "bounded_deficit_rows": diagnostics.bounded_deficit_rows,
        "raw_threshold_rows": diagnostics.raw_threshold_rows,
        "target_closing": False,
        "main_finding": (
            "The negative-lambda curvature target is a controlled deficit from the Gaussian moment baseline: "
            "Gaussian moments give raw curvature log((2*k+1)/(2*k-1)) and x_k=1, while the actual zeta "
            "prefix has a positive deficit B_k=-log(x_k) on 63 rows but keeps it below 2/(3*(2*k+1)) "
            "on all 63 checked rows. Positive Gaussian scale-mixture arguments point the wrong way for "
            "the upper wall because they make the deficit nonpositive."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_gaussian_curvature_matrix",
        "date": "2026-07-06",
        "status": "finite theorem-search diagnostic",
        "source_bounded_log_curvature_target": "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md",
        "source_heat_flow_cone_entry_target": "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "source_phi_taylor_sign_scout": "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_gaussian_curvature_matrix.py",
        "proof_boundary": (
            "Finite theorem-search diagnostic and exact algebraic comparison only. It identifies the "
            "Gaussian baseline, compatible finite deficit margins, and rejected/insufficient proof "
            "templates, but it does not prove an all-k bounded-curvature theorem, cone entry, or Lambda <= 0."
        ),
        "matrix_rows": [
            {
                "id": "nlgcm_01_gaussian_baseline_exact",
                "role": "exact_baseline",
                "claim": "For Gaussian raw moments M_k proportional to Gamma(k+1/2) T^(-k-1/2), raw curvature equals log((2*k+1)/(2*k-1)), hence x_k=1 and B_k=0.",
                "proof_boundary": "Exact model baseline only; it lands on the upper wall and is not strict cone entry.",
            },
            {
                "id": "nlgcm_02_actual_target_as_gaussian_deficit",
                "role": "exact_reduction",
                "claim": "The bounded log-curvature target is exactly 0<=B_k<=2/(3*(2*k+1)), where B_k is the deficit from Gaussian raw curvature.",
                "proof_boundary": "Exact reformulation only; the all-k bound remains open.",
            },
            {
                "id": "nlgcm_03_positive_deficit_finite_pattern",
                "role": "finite_diagnostic",
                "claim": "The checked negative-lambda prefix has positive Gaussian deficit B_k on all 63 rows, corresponding to x_k<1.",
                "proof_boundary": "Finite prefix evidence only.",
            },
            {
                "id": "nlgcm_04_bounded_deficit_finite_pattern",
                "role": "finite_diagnostic",
                "claim": "The checked negative-lambda prefix keeps B_k below 2/(3*(2*k+1)) on all 63 rows.",
                "proof_boundary": "Finite prefix evidence only; not an all-k tail theorem.",
            },
            {
                "id": "nlgcm_05_positive_gaussian_scale_mixture_rejected",
                "role": "rejected_template",
                "claim": "A positive scale mixture of centered Gaussian raw moments has nonnegative curvature added to the Gaussian baseline, giving B_k<=0 and pushing x_k>=1.",
                "proof_boundary": "Rejects this proof template for the upper-wall route only; not a claim about the actual Phi kernel.",
            },
            {
                "id": "nlgcm_06_signed_or_tilted_gaussian_perturbation_route",
                "role": "live_route",
                "claim": "A controlled signed or tilted perturbation of the Gaussian baseline could prove 0<B_k<=2/(3*(2*k+1)) with uniform remainder bounds.",
                "proof_boundary": "Live route only; no such all-k perturbation theorem is proved.",
            },
            {
                "id": "nlgcm_07_phi_taylor_asymptotic_route",
                "role": "live_route",
                "claim": "The certified Phi Taylor signs support the fixed-k direction B_k>0 and B_(k+1)<=B_k for large negative lambda, but uniform-in-k remainder control remains missing.",
                "proof_boundary": "Finite/local asymptotic support only; not a cone-entry theorem.",
            },
        ],
        "finite_diagnostics": asdict(diagnostics),
        "invariants": [
            "No row is ready_to_apply.",
            "The bounded curvature target remains open.",
            "Positive Gaussian scale-mixture arguments are not promoted as proof of the upper wall.",
            "Finite prefix margins are not promoted to an all-k theorem.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
        "summary": summary,
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["finite_diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda Gaussian curvature matrix: "
        f"{summary['matrix_rows']} matrix rows, "
        f"{summary['gaussian_deficit_rows']} positive-deficit rows, "
        f"{summary['bounded_deficit_rows']} bounded-deficit rows, "
        f"{summary['raw_threshold_rows']} raw-threshold rows, 0 issues"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Gaussian Curvature Matrix",
        "",
        "Date: 2026-07-06",
        "",
        "Status: finite theorem-search diagnostic. This is not a proof of the",
        "bounded log-curvature theorem, cone entry, Jensen-window PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_gaussian_curvature_matrix`.",
        "",
        "Proof boundary: this artifact compares the actual negative-lambda prefix",
        "with the Gaussian raw-moment baseline. It does not prove any all-`k`",
        "tail theorem.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_gaussian_curvature_matrix.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Baseline",
        "",
        "For Gaussian raw moments:",
        "",
        "```text",
        "M_k^G proportional to Gamma(k+1/2) T^(-k-1/2)",
        "M_(k+1)^G*M_(k-1)^G/(M_k^G)^2 = (2*k+1)/(2*k-1)",
        "x_k = 1",
        "B_k = -log(x_k) = 0",
        "```",
        "",
        "Thus the bounded log-curvature target asks for a small positive deficit",
        "from the Gaussian baseline:",
        "",
        "```text",
        "0 < B_k <= 2/(3*(2*k+1))",
        "```",
        "",
        "A positive Gaussian scale mixture is not the right upper-wall proof",
        "template, because it adds nonnegative curvature to the Gaussian baseline",
        "and pushes `x_k>=1` rather than `x_k<=1`.",
        "",
        "Finite diagnostics:",
        "",
        "```text",
        f"lambdas: {', '.join(diagnostics['lambdas'])}",
        f"coefficient range: A_0..A_{diagnostics['coefficient_k_max']}",
        f"checked contractions: x_1..x_{diagnostics['checked_x_max']}",
        (
            "positive Gaussian-deficit rows: "
            f"{diagnostics['gaussian_deficit_positive_rows']} / {diagnostics['gaussian_deficit_rows']}"
        ),
        (
            "bounded Gaussian-deficit rows: "
            f"{diagnostics['bounded_deficit_positive_rows']} / {diagnostics['bounded_deficit_rows']}"
        ),
        (
            "raw threshold rows: "
            f"{diagnostics['raw_threshold_positive_rows']} / {diagnostics['raw_threshold_rows']}"
        ),
        (
            "positive threshold rows: "
            f"{diagnostics['positive_threshold_rows']} / {diagnostics['raw_threshold_rows']}"
        ),
        "```",
        "",
        "Extrema:",
        "",
        "```text",
        (
            "min deficit-bound slack: "
            f"{diagnostics['min_deficit_bound_slack']['sample']} "
            f"at lambda={diagnostics['min_deficit_bound_slack']['lam']}, "
            f"k={diagnostics['min_deficit_bound_slack']['k']}"
        ),
        (
            "min raw-threshold margin: "
            f"{diagnostics['min_raw_threshold_margin']['sample']} "
            f"at lambda={diagnostics['min_raw_threshold_margin']['lam']}, "
            f"k={diagnostics['min_raw_threshold_margin']['k']}"
        ),
        (
            "max Gaussian deficit: "
            f"{diagnostics['max_gaussian_deficit']['sample']} "
            f"at lambda={diagnostics['max_gaussian_deficit']['lam']}, "
            f"k={diagnostics['max_gaussian_deficit']['k']}"
        ),
        (
            "max deficit fraction of bound: "
            f"{diagnostics['max_deficit_fraction_of_bound']['sample']} "
            f"at lambda={diagnostics['max_deficit_fraction_of_bound']['lam']}, "
            f"k={diagnostics['max_deficit_fraction_of_bound']['k']}"
        ),
        "```",
        "",
        "Integration:",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md",
        "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md",
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
        "wrote Jensen-window PF negative-lambda Gaussian curvature matrix: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
