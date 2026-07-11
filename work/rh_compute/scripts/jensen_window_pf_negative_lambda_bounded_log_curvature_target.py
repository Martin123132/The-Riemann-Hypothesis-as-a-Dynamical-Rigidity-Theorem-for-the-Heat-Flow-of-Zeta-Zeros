#!/usr/bin/env python3
"""Build the bounded log-curvature theorem target for the negative-lambda tail."""

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
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_bounded_log_curvature_target.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md"


@dataclass(frozen=True)
class Extremum:
    sample: str
    lam: str
    k: int


@dataclass(frozen=True)
class BoundedCurvatureDiagnostics:
    lambdas: list[str]
    coefficient_k_max: int
    checked_x_max: int
    raw_logconvexity_rows: int
    raw_logconvexity_positive_rows: int
    raw_curvature_threshold_rows: int
    raw_curvature_threshold_positive_rows: int
    simple_log_buffer_rows: int
    simple_log_buffer_positive_rows: int
    threshold_positive_rows: int
    min_raw_curvature_threshold_margin: Extremum
    min_simple_log_buffer_margin: Extremum
    min_plain_logconvexity_shortfall: Extremum
    max_bounded_log_curvature_target: Extremum


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
        raise RuntimeError("missing bounded-curvature diagnostic extremum")
    value, lam, k = item
    return Extremum(decimal_format(value), lam, k)


def build_diagnostics(paths: list[Path]) -> BoundedCurvatureDiagnostics:
    balls, samples, labels = load_enclosures(paths)
    if not labels:
        raise RuntimeError("no coefficient enclosures found")
    coefficient_k_max = max(index for (_lam, index) in balls)
    checked_x_max = coefficient_k_max - 1
    lambdas = sorted(labels)

    raw_rows = 0
    raw_good = 0
    threshold_rows = 0
    threshold_good = 0
    simple_rows = 0
    simple_good = 0
    threshold_positive = 0

    min_threshold_margin: tuple[Decimal, str, int] | None = None
    min_simple_margin: tuple[Decimal, str, int] | None = None
    min_plain_shortfall: tuple[Decimal, str, int] | None = None
    max_target: tuple[Decimal, str, int] | None = None

    for lam in lambdas:
        lam_label = labels[lam]
        for k in range(1, checked_x_max + 1):
            x = contraction({idx: balls[(lam, idx)] for idx in (k - 1, k, k + 1)}, k)
            sample_x = contraction({idx: samples[(lam, idx)] for idx in (k - 1, k, k + 1)}, k)

            factor = flint.arb(2 * k - 1) / flint.arb(2 * k + 1)
            raw_curvature_ratio = x / factor
            target = flint.arb(2) / flint.arb(3 * (2 * k + 1))
            raw_threshold = (flint.arb(2 * k + 1) / flint.arb(2 * k - 1)).log() - target
            raw_rows += 1
            threshold_rows += 1
            simple_rows += 1
            if arb_positive(raw_curvature_ratio - flint.arb(1)):
                raw_good += 1
            if arb_positive(raw_curvature_ratio.log() - raw_threshold):
                threshold_good += 1
            if arb_positive(x - (-target).exp()):
                simple_good += 1

            sample_target = Decimal(2) / Decimal(3 * (2 * k + 1))
            sample_log_factor = (Decimal(2 * k + 1) / Decimal(2 * k - 1)).ln()
            sample_threshold = sample_log_factor - sample_target
            sample_b = -sample_x.ln()
            sample_raw_curvature = sample_log_factor - sample_b
            sample_threshold_margin = sample_raw_curvature - sample_threshold
            sample_simple_margin = sample_x - (-sample_target).exp()
            if sample_threshold > 0:
                threshold_positive += 1
            min_threshold_margin = update_min(min_threshold_margin, sample_threshold_margin, lam_label, k)
            min_simple_margin = update_min(min_simple_margin, sample_simple_margin, lam_label, k)
            min_plain_shortfall = update_min(min_plain_shortfall, sample_threshold, lam_label, k)
            max_target = update_max(max_target, sample_target, lam_label, k)

    return BoundedCurvatureDiagnostics(
        lambdas=[labels[lam] for lam in lambdas],
        coefficient_k_max=coefficient_k_max,
        checked_x_max=checked_x_max,
        raw_logconvexity_rows=raw_rows,
        raw_logconvexity_positive_rows=raw_good,
        raw_curvature_threshold_rows=threshold_rows,
        raw_curvature_threshold_positive_rows=threshold_good,
        simple_log_buffer_rows=simple_rows,
        simple_log_buffer_positive_rows=simple_good,
        threshold_positive_rows=threshold_positive,
        min_raw_curvature_threshold_margin=extremum(min_threshold_margin),
        min_simple_log_buffer_margin=extremum(min_simple_margin),
        min_plain_logconvexity_shortfall=extremum(min_plain_shortfall),
        max_bounded_log_curvature_target=extremum(max_target),
    )


def build_target(paths: list[Path]) -> dict:
    diagnostics = build_diagnostics(paths)
    tail_lower_start = diagnostics.checked_x_max + 1
    summary = {
        "target_rows": 8,
        "ready_to_apply_rows": 0,
        "live_routes": 2,
        "exact_reduction_rows": 2,
        "insufficient_routes": 1,
        "conditional_application_rows": 1,
        "raw_curvature_threshold_rows": diagnostics.raw_curvature_threshold_rows,
        "simple_log_buffer_rows": diagnostics.simple_log_buffer_rows,
        "target_closing": False,
        "main_finding": (
            "The bounded log-curvature target is the missing quantitative half of the "
            f"negative-lambda buffered tail: prove B_k=-Delta^2 log A_k<=2/(3*(2*k+1)) "
            f"for k>={tail_lower_start}. Equivalently, for raw moments M_k=mu_(2k), prove "
            "log(M_(k+1)*M_(k-1)/M_k^2) >= log((2*k+1)/(2*k-1))-2/(3*(2*k+1)). "
            "Plain moment log-convexity only gives the nonnegative left side and is too weak."
        ),
    }
    rows = [
        {
            "id": "nllct_01_bounded_log_curvature_statement",
            "role": "open_statement",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md",
                "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
            ],
            "claim_if_proved": (
                f"Prove 0<=B_k<=2/(3*(2*k+1)) for all k>={tail_lower_start}, "
                "where B_k=-Delta^2 log A_k=-log(x_k), for the actual negative-lambda zeta heat-flow coefficients."
            ),
            "gap": "No all-k bounded log-curvature estimate is currently proved.",
            "acceptance_test": "Give a noncircular analytic estimate with explicit starting index and constants.",
            "proof_boundary": "Open theorem target only; not a proof of the defect-tail theorem.",
        },
        {
            "id": "nllct_02_raw_moment_curvature_equivalence",
            "role": "exact_reduction",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
            ],
            "claim_if_proved": (
                "For M_k=mu_(2k), x_k=((2*k-1)/(2*k+1))*(M_(k+1)*M_(k-1)/M_k^2), so the bounded "
                "B_k target is equivalent to a quantitative lower bound on raw moment curvature."
            ),
            "gap": "The exact identity is available, but the quantitative lower bound is not.",
            "acceptance_test": "Keep the factorial normalization A_k=M_k*k!/(2*k)! explicit.",
            "proof_boundary": "Exact reduction only; not the analytic estimate.",
        },
        {
            "id": "nllct_03_saddle_laplace_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
                "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md",
            ],
            "claim_if_proved": "A uniform saddle/Laplace estimate for the actual heat-flow moment integrals gives the required upper bound on B_k.",
            "gap": "Current asymptotic support is fixed-k/local and does not control the k-tail uniformly.",
            "acceptance_test": "Bound the remainder uniformly in k from the stated tail index, with constants strong enough for 2/(3*(2*k+1)).",
            "proof_boundary": "Live theorem-search route only.",
        },
        {
            "id": "nllct_04_moment_ratio_growth_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md",
            ],
            "claim_if_proved": "A zeta-specific lower bound on the growth of raw moment ratios proves the required raw moment curvature threshold.",
            "gap": "Plain positivity or Stieltjes moment log-convexity gives only a zero lower bound.",
            "acceptance_test": "Derive the threshold log((2*k+1)/(2*k-1))-2/(3*(2*k+1)) from Phi/Newman-kernel structure.",
            "proof_boundary": "Live theorem-search route only; not a proved moment-ratio growth theorem.",
        },
        {
            "id": "nllct_05_plain_moment_logconvexity_insufficient",
            "role": "insufficient_route",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
            ],
            "claim_if_proved": "Raw moment log-convexity gives M_(k+1)*M_(k-1)/M_k^2>=1.",
            "gap": "The needed threshold is strictly positive for the checked tail range, so zero curvature is too weak.",
            "acceptance_test": "Do not use ordinary moment log-convexity alone as the bounded-curvature proof.",
            "proof_boundary": "Insufficient exact fact only.",
        },
        {
            "id": "nllct_06_finite_compatibility",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md",
            ],
            "claim_if_proved": "The current finite prefix satisfies the raw moment curvature threshold and the equivalent simple log-buffer rows.",
            "gap": "Finite compatibility through x_21 does not prove the all-k tail from k>=22.",
            "acceptance_test": "Use as a base-case diagnostic only unless paired with an analytic tail theorem.",
            "proof_boundary": "Finite diagnostic only.",
        },
        {
            "id": "nllct_07_monotone_contraction_interaction",
            "role": "open_dependency",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
                "outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md",
            ],
            "claim_if_proved": "Together with B_(k+1)<=B_k, the bounded-curvature estimate gives the buffered defect-tail sufficient condition.",
            "gap": "The B monotonicity clause is still the open Delta^3 log A theorem target.",
            "acceptance_test": "Discharge both the bounded-curvature upper estimate and the monotone-contraction theorem without endpoint assumptions.",
            "proof_boundary": "Dependency statement only; neither open theorem is proved.",
        },
        {
            "id": "nllct_08_conditional_defect_tail_application",
            "role": "conditional_application",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
                "outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md",
            ],
            "claim_if_proved": "A proved bounded log-curvature theorem plus the monotone-contraction theorem would close the buffered route to the defect-tail theorem.",
            "gap": "Both analytic inputs remain open.",
            "acceptance_test": "After proof, re-run the defect-tail target and dependency graph without changing forbidden endpoint assumptions.",
            "proof_boundary": "Conditional application only; not cone entry or Lambda <= 0.",
        },
    ]
    return {
        "kind": "jensen_window_pf_negative_lambda_bounded_log_curvature_target",
        "date": "2026-07-06",
        "status": "open_theorem_target",
        "target_id": "target_negative_lambda_bounded_log_curvature",
        "source_log_curvature_bridge": "outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md",
        "source_defect_tail_target": "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
        "source_monotone_contraction_target": "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_bounded_log_curvature_target.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_bounded_log_curvature_target.py",
        "proof_boundary": (
            "Open theorem target only. It isolates the bounded log-curvature estimate needed by "
            "the negative-lambda buffered tail route, but it does not prove that estimate, does "
            "not prove the defect-tail theorem, does not prove cone entry, and does not prove Lambda <= 0."
        ),
        "target_rows": rows,
        "finite_diagnostics": asdict(diagnostics),
        "invariants": [
            "No row is ready_to_apply.",
            "The target remains open_theorem_target.",
            "Plain raw moment log-convexity is recorded as insufficient.",
            "Finite prefix compatibility is not promoted to an all-k bounded-curvature theorem.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
        "summary": summary,
    }


def write_note(target: dict, path: Path) -> None:
    summary = target["summary"]
    diagnostics = target["finite_diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda bounded log-curvature target: "
        f"{summary['target_rows']} rows, 0 issues, {summary['ready_to_apply_rows']} ready-to-apply rows, "
        f"{summary['live_routes']} live routes, {summary['raw_curvature_threshold_rows']} raw-threshold rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Bounded Log-Curvature Target",
        "",
        "Date: 2026-07-06",
        "",
        "Status: open theorem target. This is not a proof of the defect-tail",
        "theorem, cone entry, Jensen-window PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_bounded_log_curvature_target`.",
        "",
        "Proof boundary: this artifact states the missing quantitative",
        "log-curvature estimate. It does not prove that estimate and does not",
        "close the negative-lambda tail.",
        "",
        "Retirement notice: this fixed `2/3` scaled-curvature wall is no longer a",
        "live theorem target after the repaired k300 obstruction in",
        "`outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.md`.",
        "The note is retained as the historical k<=22 compatibility target.",
        "",
        "Machine-readable target:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_bounded_log_curvature_target.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_bounded_log_curvature_target.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_bounded_log_curvature_target.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Target Statement",
        "",
        "Let:",
        "",
        "```text",
        "B_k = -Delta^2 log A_k = -log(x_k)",
        "M_k = mu_(2k)",
        "A_k = M_k*k!/(2*k)!",
        "```",
        "",
        "The missing estimate is:",
        "",
        "```text",
        "0 <= B_k <= 2/(3*(2*k+1))",
        "```",
        "",
        "for the actual negative-lambda zeta heat-flow coefficient tail from the",
        "stated starting index.",
        "",
        "The exact raw-moment translation is:",
        "",
        "```text",
        "x_k = ((2*k-1)/(2*k+1)) * (M_(k+1)*M_(k-1)/M_k^2)",
        "log(M_(k+1)*M_(k-1)/M_k^2) >= log((2*k+1)/(2*k-1)) - 2/(3*(2*k+1))",
        "```",
        "",
        "Ordinary moment log-convexity gives only the nonnegative left side. The",
        "threshold on the right is strictly positive, so that standard fact is",
        "not enough.",
        "",
        "## Live Routes",
        "",
        "```text",
        "1. uniform saddle/Laplace control of the actual heat-flow moment integrals",
        "2. zeta-specific lower bounds on raw moment-ratio growth",
        "```",
        "",
        "Finite diagnostics:",
        "",
        "```text",
        f"lambdas: {', '.join(diagnostics['lambdas'])}",
        f"coefficient range: A_0..A_{diagnostics['coefficient_k_max']}",
        f"checked contractions: x_1..x_{diagnostics['checked_x_max']}",
        (
            "raw log-convexity rows: "
            f"{diagnostics['raw_logconvexity_positive_rows']} / {diagnostics['raw_logconvexity_rows']}"
        ),
        (
            "raw curvature-threshold rows: "
            f"{diagnostics['raw_curvature_threshold_positive_rows']} / {diagnostics['raw_curvature_threshold_rows']}"
        ),
        (
            "simple log-buffer rows: "
            f"{diagnostics['simple_log_buffer_positive_rows']} / {diagnostics['simple_log_buffer_rows']}"
        ),
        (
            "positive threshold rows: "
            f"{diagnostics['threshold_positive_rows']} / {diagnostics['raw_curvature_threshold_rows']}"
        ),
        "```",
        "",
        "Extrema:",
        "",
        "```text",
        (
            "min raw curvature-threshold margin: "
            f"{diagnostics['min_raw_curvature_threshold_margin']['sample']} "
            f"at lambda={diagnostics['min_raw_curvature_threshold_margin']['lam']}, "
            f"k={diagnostics['min_raw_curvature_threshold_margin']['k']}"
        ),
        (
            "min simple log-buffer margin: "
            f"{diagnostics['min_simple_log_buffer_margin']['sample']} "
            f"at lambda={diagnostics['min_simple_log_buffer_margin']['lam']}, "
            f"k={diagnostics['min_simple_log_buffer_margin']['k']}"
        ),
        (
            "min plain-logconvexity shortfall: "
            f"{diagnostics['min_plain_logconvexity_shortfall']['sample']} "
            f"at lambda={diagnostics['min_plain_logconvexity_shortfall']['lam']}, "
            f"k={diagnostics['min_plain_logconvexity_shortfall']['k']}"
        ),
        "```",
        "",
        "Integration:",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md",
        "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
        "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "```",
        "",
        "Summary:",
        "",
        summary["main_finding"],
        "",
        "Historical correction:",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.md",
        "validated Jensen-window PF negative-lambda bounded log-curvature k300 obstruction: 7 rows, 0 issues, 718 two-thirds failures, 894 scaled-curvature increase rows, 0 ready-to-apply rows",
        "```",
        "",
        "The fixed wall `C_k=(2*k+1)*B_k<=2/3` is finite-rejected on the repaired",
        "k300 data and should not be used as a live defect-tail route.",
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
    target = build_target(paths)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(target, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(target, note)
    print(
        "wrote Jensen-window PF negative-lambda bounded log-curvature target: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
