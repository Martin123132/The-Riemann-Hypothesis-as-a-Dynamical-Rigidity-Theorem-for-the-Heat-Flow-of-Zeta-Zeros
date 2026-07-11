#!/usr/bin/env python3
"""Build a log-curvature bridge for the negative-lambda defect tail."""

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
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_log_curvature_bridge.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md"


@dataclass(frozen=True)
class Extremum:
    sample: str
    lam: str
    k: int


@dataclass(frozen=True)
class LogCurvatureDiagnostics:
    lambdas: list[str]
    coefficient_k_max: int
    checked_x_max: int
    positive_curvature_rows: int
    positive_curvature_positive_rows: int
    exact_defect_buffer_rows: int
    exact_defect_buffer_positive_rows: int
    simple_log_buffer_rows: int
    simple_log_buffer_positive_rows: int
    curvature_monotone_rows: int
    curvature_monotone_positive_rows: int
    min_exact_defect_buffer_margin: Extremum
    min_simple_log_buffer_margin: Extremum
    min_simple_log_slack: Extremum
    min_exact_log_slack: Extremum
    min_curvature_monotone_gap: Extremum
    max_log_curvature: Extremum


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
        raise RuntimeError("missing log-curvature diagnostic extremum")
    value, lam, k = item
    return Extremum(decimal_format(value), lam, k)


def build_diagnostics(paths: list[Path]) -> LogCurvatureDiagnostics:
    balls, samples, labels = load_enclosures(paths)
    if not labels:
        raise RuntimeError("no coefficient enclosures found")
    coefficient_k_max = max(index for (_lam, index) in balls)
    checked_x_max = coefficient_k_max - 1
    lambdas = sorted(labels)

    positive_rows = 0
    positive_good = 0
    exact_rows = 0
    exact_good = 0
    simple_rows = 0
    simple_good = 0
    monotone_rows = 0
    monotone_good = 0

    min_exact_margin: tuple[Decimal, str, int] | None = None
    min_simple_margin: tuple[Decimal, str, int] | None = None
    min_simple_slack: tuple[Decimal, str, int] | None = None
    min_exact_slack: tuple[Decimal, str, int] | None = None
    min_monotone_gap: tuple[Decimal, str, int] | None = None
    max_log_curvature: tuple[Decimal, str, int] | None = None

    for lam in lambdas:
        lam_label = labels[lam]
        arb_x: dict[int, flint.arb] = {}
        sample_x: dict[int, Decimal] = {}
        for k in range(1, checked_x_max + 1):
            arb_x[k] = contraction({idx: balls[(lam, idx)] for idx in (k - 1, k, k + 1)}, k)
            sample_x[k] = contraction({idx: samples[(lam, idx)] for idx in (k - 1, k, k + 1)}, k)

            c = flint.arb(2) / flint.arb(3 * (2 * k + 1))
            exact_margin = arb_x[k] - (flint.arb(1) - c)
            simple_margin = arb_x[k] - (-c).exp()
            positive_rows += 1
            exact_rows += 1
            simple_rows += 1
            if arb_positive(arb_x[k]) and arb_positive(flint.arb(1) - arb_x[k]):
                positive_good += 1
            if arb_positive(exact_margin):
                exact_good += 1
            if arb_positive(simple_margin):
                simple_good += 1

            sample_c = Decimal(2) / Decimal(3 * (2 * k + 1))
            sample_log_curvature = -sample_x[k].ln()
            sample_exact_margin = sample_x[k] - (Decimal(1) - sample_c)
            sample_simple_margin = sample_x[k] - (-sample_c).exp()
            sample_simple_slack = sample_c - sample_log_curvature
            sample_exact_slack = -(Decimal(1) - sample_c).ln() - sample_log_curvature
            min_exact_margin = update_min(min_exact_margin, sample_exact_margin, lam_label, k)
            min_simple_margin = update_min(min_simple_margin, sample_simple_margin, lam_label, k)
            min_simple_slack = update_min(min_simple_slack, sample_simple_slack, lam_label, k)
            min_exact_slack = update_min(min_exact_slack, sample_exact_slack, lam_label, k)
            max_log_curvature = update_max(max_log_curvature, sample_log_curvature, lam_label, k)

        for k in range(1, checked_x_max):
            monotone_rows += 1
            gap = arb_x[k + 1] - arb_x[k]
            if arb_positive(gap):
                monotone_good += 1
            min_monotone_gap = update_min(min_monotone_gap, sample_x[k + 1] - sample_x[k], lam_label, k)

    return LogCurvatureDiagnostics(
        lambdas=[labels[lam] for lam in lambdas],
        coefficient_k_max=coefficient_k_max,
        checked_x_max=checked_x_max,
        positive_curvature_rows=positive_rows,
        positive_curvature_positive_rows=positive_good,
        exact_defect_buffer_rows=exact_rows,
        exact_defect_buffer_positive_rows=exact_good,
        simple_log_buffer_rows=simple_rows,
        simple_log_buffer_positive_rows=simple_good,
        curvature_monotone_rows=monotone_rows,
        curvature_monotone_positive_rows=monotone_good,
        min_exact_defect_buffer_margin=extremum(min_exact_margin),
        min_simple_log_buffer_margin=extremum(min_simple_margin),
        min_simple_log_slack=extremum(min_simple_slack),
        min_exact_log_slack=extremum(min_exact_slack),
        min_curvature_monotone_gap=extremum(min_monotone_gap),
        max_log_curvature=extremum(max_log_curvature),
    )


def build_artifact(paths: list[Path]) -> dict:
    diagnostics = build_diagnostics(paths)
    tail_lower_start = diagnostics.checked_x_max + 1
    monotone_bridge_k = diagnostics.checked_x_max
    summary = {
        "bridge_rows": 5,
        "positive_curvature_rows": diagnostics.positive_curvature_rows,
        "exact_defect_buffer_rows": diagnostics.exact_defect_buffer_rows,
        "simple_log_buffer_rows": diagnostics.simple_log_buffer_rows,
        "curvature_monotone_rows": diagnostics.curvature_monotone_rows,
        "ready_to_apply_rows": 0,
        "live_reduction_rows": 2,
        "target_closing": False,
        "main_finding": (
            "The buffered defect-tail condition reduces to a cleaner log-curvature target: "
            f"prove 0<=B_k<=2/(3*(2*k+1)) for k>={tail_lower_start}, where "
            "B_k=-Delta^2 log A_k=-log(x_k), and prove B_(k+1)<=B_k from "
            f"k>={monotone_bridge_k}. The second clause is exactly the Delta^3 log A "
            "monotone-contraction target, while the bounded-curvature upper estimate remains open."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_log_curvature_bridge",
        "date": "2026-07-06",
        "status": "finite theorem-search diagnostic",
        "source_defect_recurrence_scout": "outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md",
        "source_defect_tail_target": "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
        "source_monotone_contraction_target": "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_log_curvature_bridge.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_log_curvature_bridge.py",
        "proof_boundary": (
            "Finite theorem-search diagnostic and exact algebraic reduction only. It shows how a "
            "bounded log-curvature theorem would imply the buffered defect-tail condition, but it "
            "does not prove the all-k bound, does not prove the monotone-contraction theorem, does "
            "not prove cone entry, and does not prove Lambda <= 0."
        ),
        "bridge_rows": [
            {
                "id": "nllcb_01_defect_log_equivalence",
                "role": "exact_bridge",
                "claim": "For 0<c<1 and 0<x_k<=1, d_k=1-x_k<=c is equivalent to B_k=-log(x_k)<=-log(1-c).",
                "proof_boundary": "Exact algebra only; it does not prove any zeta tail bound.",
            },
            {
                "id": "nllcb_02_simple_log_curvature_sufficient",
                "role": "live_sufficient_condition",
                "claim": "The stronger bound 0<=B_k<=2/(3*(2*k+1)) implies the buffered defect condition d_k<=2/(3*(2*k+1)).",
                "proof_boundary": "Sufficient condition only; the all-k bound remains unproved.",
            },
            {
                "id": "nllcb_03_monotone_defect_log_curvature",
                "role": "exact_bridge",
                "claim": "For 0<x_k<=1, d_(k+1)<=d_k is equivalent to B_(k+1)<=B_k.",
                "proof_boundary": "Exact translation only; not a proof of monotonicity for zeta coefficients.",
            },
            {
                "id": "nllcb_04_delta3_reduction",
                "role": "live_reduction",
                "claim": "B_(k+1)<=B_k is equivalent to Delta^3 log A_{k-1}>=0, matching the existing monotone-contraction theorem target.",
                "source": "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
                "proof_boundary": "Reduction to an open theorem target only.",
            },
            {
                "id": "nllcb_05_finite_compatibility",
                "role": "finite_diagnostic",
                "claim": "The current certified prefix validates the exact defect buffer, the stronger simple log-buffer, and curvature monotonicity on all checked rows.",
                "proof_boundary": "Finite compatibility only; not an all-k theorem.",
            },
        ],
        "finite_diagnostics": asdict(diagnostics),
        "invariants": [
            "No row is ready_to_apply for the defect-tail theorem.",
            "The bounded log-curvature estimate is not proved for all k.",
            "The Delta^3 log A monotone-contraction theorem target remains open.",
            "Finite prefix compatibility is not promoted to cone entry.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
        "summary": summary,
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["finite_diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda log-curvature bridge: "
        f"{summary['simple_log_buffer_rows']} simple log-buffer rows, "
        f"{summary['exact_defect_buffer_rows']} exact defect-buffer rows, "
        f"{summary['curvature_monotone_rows']} curvature-monotone rows, "
        f"{summary['bridge_rows']} bridge rows, 0 issues"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Log-Curvature Bridge",
        "",
        "Date: 2026-07-06",
        "",
        "Status: finite theorem-search diagnostic. This is not a proof of the",
        "defect-tail theorem, cone entry, Jensen-window PF-infinity, RH, or",
        "`Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_log_curvature_bridge`.",
        "",
        "Proof boundary: this artifact translates the buffered defect tail into",
        "log-curvature inequalities and checks finite compatibility. It does",
        "not prove the required all-`k` log-curvature theorem.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_log_curvature_bridge.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_log_curvature_bridge.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_log_curvature_bridge.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Log-Curvature Translation",
        "",
        "Let:",
        "",
        "```text",
        "B_k = -Delta^2 log A_k = -log(x_k)",
        "x_k = (A_{k+1}/A_k)/(A_k/A_{k-1})",
        "d_k = 1 - x_k",
        "```",
        "",
        "A sufficient all-tail theorem for the buffered route is:",
        "",
        "```text",
        "0 <= B_k <= 2/(3*(2*k+1))",
        "B_(k+1) <= B_k",
        "```",
        "",
        "The second line is the existing monotone-contraction sign target in",
        "log-curvature language:",
        "",
        "```text",
        "Delta^3 log A_{k-1} >= 0",
        "```",
        "",
        "Thus the remaining new analytic burden is the bounded-curvature upper",
        "estimate `B_k <= 2/(3*(2*k+1))` on the actual zeta heat-flow tail.",
        "",
        "Finite diagnostics:",
        "",
        "```text",
        f"lambdas: {', '.join(diagnostics['lambdas'])}",
        f"coefficient range: A_0..A_{diagnostics['coefficient_k_max']}",
        f"checked contractions: x_1..x_{diagnostics['checked_x_max']}",
        (
            "positive curvature rows: "
            f"{diagnostics['positive_curvature_positive_rows']} / {diagnostics['positive_curvature_rows']}"
        ),
        (
            "exact defect-buffer rows: "
            f"{diagnostics['exact_defect_buffer_positive_rows']} / {diagnostics['exact_defect_buffer_rows']}"
        ),
        (
            "simple log-buffer rows: "
            f"{diagnostics['simple_log_buffer_positive_rows']} / {diagnostics['simple_log_buffer_rows']}"
        ),
        (
            "curvature-monotone rows: "
            f"{diagnostics['curvature_monotone_positive_rows']} / {diagnostics['curvature_monotone_rows']}"
        ),
        "```",
        "",
        "Extrema:",
        "",
        "```text",
        (
            "min simple log slack: "
            f"{diagnostics['min_simple_log_slack']['sample']} "
            f"at lambda={diagnostics['min_simple_log_slack']['lam']}, k={diagnostics['min_simple_log_slack']['k']}"
        ),
        (
            "min exact log slack: "
            f"{diagnostics['min_exact_log_slack']['sample']} "
            f"at lambda={diagnostics['min_exact_log_slack']['lam']}, k={diagnostics['min_exact_log_slack']['k']}"
        ),
        (
            "max log curvature: "
            f"{diagnostics['max_log_curvature']['sample']} "
            f"at lambda={diagnostics['max_log_curvature']['lam']}, k={diagnostics['max_log_curvature']['k']}"
        ),
        (
            "min curvature-monotone gap: "
            f"{diagnostics['min_curvature_monotone_gap']['sample']} "
            f"at lambda={diagnostics['min_curvature_monotone_gap']['lam']}, k={diagnostics['min_curvature_monotone_gap']['k']}"
        ),
        "```",
        "",
        "Integration:",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md",
        "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
        "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
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
        "wrote Jensen-window PF negative-lambda log-curvature bridge: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
