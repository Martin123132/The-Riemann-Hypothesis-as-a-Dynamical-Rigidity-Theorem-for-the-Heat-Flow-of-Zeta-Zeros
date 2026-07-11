#!/usr/bin/env python3
"""Validate the worst-row full relative-Gaussian expectation certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate import (  # noqa: E402
    DEFAULT_COMPACT_JSON,
    DEFAULT_FAR_TAIL_JSON,
    DEFAULT_FIRST_OMITTED_JSON,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
    result_line,
)


REQUIRED_ROW_IDS = {
    "nlrgwrfec_01_compact_source_import",
    "nlrgwrfec_02_far_tail_source_import",
    "nlrgwrfec_03_global_n_tail_certificate",
    "nlrgwrfec_04_normalization_correction",
    "nlrgwrfec_05_complete_expectation_composition",
    "nlrgwrfec_06_below_first_omitted_certificate",
    "nlrgwrfec_07_all_row_handoff",
    "nlrgwrfec_08_acceptance_gate",
}

ALLOWED_ROLES = {
    "interval_certificate_import",
    "tail_certificate_import",
    "analytic_tail_certificate",
    "normalization_certificate",
    "full_expectation_certificate",
    "ratio_certificate",
    "intervalization_handoff",
    "acceptance_gate",
}

ALLOWED_READINESS = {
    "available_interval_certificate",
    "available_tail_certificate",
    "available_normalization_certificate",
    "available_worst_row_expectation_certificate",
    "available_worst_row_ratio_certificate",
    "not_ready_to_apply",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Full-Expectation Certificate",
    "Status: worst-row full-expectation certificate",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate.py",
    (
        "validated Jensen-window PF negative-lambda relative-Gaussian worst-row full-expectation certificate: "
        "8 rows, 0 issues, 3 composed sources, 2 below-one full ratios, "
        "0 open worst-row integral sources, 0 ready-to-apply rows"
    ),
    "tail start n: 31",
    "global extension certified: True",
    "both c0 lower bounds certified: True",
    "full value certified negative: True",
    "full derivative certified negative: True",
    "value ratio / first omitted upper: 0.9707100590203297651200076500707905236731092924001689841866982789573976",
    "derivative ratio / first omitted upper: 0.9693567774758049120828364647202533809659101015524784238953305230805692",
    "both full ratios below one: True",
    "complete worst-row expectation certified: True",
    "quadrature needed for worst-row expectation: False",
    "Apply the same direct compact/global-tail certificate to the other recorded T/index rows.",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class FullExpectationIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> FullExpectationIssue:
    return FullExpectationIssue(section=section, issue=name, detail=detail)


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("e", "E"))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[FullExpectationIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[FullExpectationIssue]:
    issues: list[FullExpectationIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate"
    ):
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "worst-row full-expectation certificate":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_compact_x_moment_taylor_certificate",
        "source_compact_x_moment_taylor_json",
        "source_endpoint_x_moment_taylor_certificate",
        "source_far_tail_split_certificate",
        "source_far_tail_split_json",
        "source_first_omitted_denominator_certificate",
        "source_first_omitted_denominator_json",
        "source_coefficient_core_certificate",
        "source_phi_tail_grid_certificate",
        "source_node_c0_range_certificate",
        "source_finite_plus_tail_budget_certificate",
        "source_quadrature_remainder_route_matrix",
        "source_intervalization_target",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "source_formal_tail_obstruction",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "complete true relative-gaussian expectation certificate",
        "single t=10000, f_21 worst row",
        "finite n<=30 compact and far-tail integrals",
        "global n>=31",
        "both first-omitted ratios below one",
        "does not certify the other finite-grid rows",
        "finite-grid-to-collar",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[FullExpectationIssue]:
    try:
        recomputed = build_artifact(
            DEFAULT_COMPACT_JSON,
            DEFAULT_FAR_TAIL_JSON,
            DEFAULT_FIRST_OMITTED_JSON,
        )
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[FullExpectationIssue] = []
    for key in ("diagnostics", "matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputation"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[FullExpectationIssue], int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[FullExpectationIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "not-list", repr(type(rows)))], 0, 0, 0
    seen = {row.get("id") for row in rows if isinstance(row, dict)}
    missing = REQUIRED_ROW_IDS - seen
    extra = seen - REQUIRED_ROW_IDS
    if missing:
        issues.append(issue("matrix_rows", "missing-row-ids", repr(sorted(missing))))
    if extra:
        issues.append(issue("matrix_rows", "extra-row-ids", repr(sorted(extra))))
    expectation_rows = 0
    ratio_rows = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("matrix_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "claim", "proof_boundary", "source_artifacts"):
            if key not in row or row[key] in (None, ""):
                issues.append(issue(row_id, "missing-field", key))
        if row.get("role") not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(row.get("role"))))
        if row.get("readiness") not in ALLOWED_READINESS:
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        if row.get("readiness") == "available_worst_row_expectation_certificate":
            expectation_rows += 1
        if row.get("readiness") == "available_worst_row_ratio_certificate":
            ratio_rows += 1
        if not any(
            marker in str(row.get("proof_boundary", "")).lower()
            for marker in ("only", "not ", "does not")
        ):
            issues.append(issue(row_id, "weak-proof-boundary", repr(row.get("proof_boundary"))))
        refs = row.get("source_artifacts", [])
        if not isinstance(refs, list) or not refs:
            issues.append(issue(row_id, "missing-source-artifacts", repr(refs)))
        else:
            for ref in refs:
                issues.extend(validate_ref(row_id, ref))
    if expectation_rows != 1:
        issues.append(issue("matrix_rows", "bad-expectation-row-count", str(expectation_rows)))
    if ratio_rows != 1:
        issues.append(issue("matrix_rows", "bad-ratio-row-count", str(ratio_rows)))
    return issues, len(rows), expectation_rows, ratio_rows


def validate_summary(
    artifact: dict,
    row_count: int,
    expectation_rows: int,
    ratio_rows: int,
) -> list[FullExpectationIssue]:
    summary = artifact.get("summary", {})
    diagnostics = artifact.get("diagnostics", {})
    issues: list[FullExpectationIssue] = []
    expected = {
        "matrix_rows": 8,
        "T": 10000,
        "index": 21,
        "precision_bits": 8192,
        "tail_start_n": 31,
        "global_n_tail_extension_certified": True,
        "finite_and_full_c0_lowers_certified": True,
        "full_value_certified_negative": True,
        "full_derivative_certified_negative": True,
        "both_full_expectation_ratios_below_one": True,
        "composed_source_count": 3,
        "below_one_ratio_count": 2,
        "complete_worst_row_expectation_certified": True,
        "quadrature_needed_for_worst_row_expectation": False,
        "remaining_worst_row_integral_source_count": 0,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, expected_value in expected.items():
        if summary.get(key) != expected_value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {expected_value!r}"))
    if row_count != summary.get("matrix_rows"):
        issues.append(issue("summary", "row-count-mismatch", f"{row_count} != {summary.get('matrix_rows')!r}"))
    if expectation_rows != 1 or ratio_rows != 1:
        issues.append(issue("summary", "certificate-row-mismatch", f"{expectation_rows}, {ratio_rows}"))
    thresholds = (
        ("value_n_tail_global_bound_upper", None, Decimal("1.1e-1300")),
        ("derivative_n_tail_global_bound_upper", None, Decimal("7e-1295")),
        ("c0_n_tail_bound_upper", None, Decimal("1.3e-1304")),
        ("value_normalized_n_tail_correction_upper", None, Decimal("1.2e-1297")),
        ("derivative_normalized_n_tail_correction_upper", None, Decimal("7.5e-1295")),
        ("value_full_ratio_to_first_omitted_upper", Decimal("0.970"), Decimal("0.971")),
        ("derivative_full_ratio_to_first_omitted_upper", Decimal("0.969"), Decimal("0.970")),
        ("value_remaining_margin_below_one_lower", Decimal("0.029"), Decimal("0.030")),
        ("derivative_remaining_margin_below_one_lower", Decimal("0.030"), Decimal("0.031")),
    )
    for key, lower, upper in thresholds:
        try:
            value = dec(summary[key])
        except Exception as exc:
            issues.append(issue("summary", f"bad-numeric-{key}", f"{type(exc).__name__}: {exc}"))
            continue
        if lower is not None and value <= lower:
            issues.append(issue("summary", f"too-small-{key}", str(value)))
        if upper is not None and value >= upper:
            issues.append(issue("summary", f"too-large-{key}", str(value)))
    if dec(diagnostics.get("value_x_ge_1_monotonicity_margin_lower", "0")) <= Decimal("600000"):
        issues.append(issue("diagnostics", "weak-value-monotonicity-margin", repr(diagnostics.get("value_x_ge_1_monotonicity_margin_lower"))))
    if dec(diagnostics.get("derivative_x_ge_1_monotonicity_margin_lower", "0")) <= Decimal("600000"):
        issues.append(issue("diagnostics", "weak-derivative-monotonicity-margin", repr(diagnostics.get("derivative_x_ge_1_monotonicity_margin_lower"))))
    if dec(diagnostics.get("finite_c0_lower", "0")) <= Decimal("0.44"):
        issues.append(issue("diagnostics", "finite-c0-lower-too-small", repr(diagnostics.get("finite_c0_lower"))))
    if dec(diagnostics.get("full_c0_lower", "0")) <= Decimal("0.44"):
        issues.append(issue("diagnostics", "full-c0-lower-too-small", repr(diagnostics.get("full_c0_lower"))))
    if diagnostics.get("remaining_worst_row_integral_sources") != []:
        issues.append(issue("diagnostics", "remaining-integral-sources", repr(diagnostics.get("remaining_worst_row_integral_sources"))))
    tail_argument = str(diagnostics.get("global_tail_argument", "")).lower()
    for required in ("0<=x<=1", "x>=1", "decreasing", "n>=31"):
        if required not in tail_argument:
            issues.append(issue("diagnostics", "weak-global-tail-argument", required))
    normalization = str(diagnostics.get("normalization_correction_formula", ""))
    for required in ("Phi/Phi0", "Phi30/Phi30_0", "Phi0_lower", "Phi30_0_lower"):
        if required not in normalization:
            issues.append(issue("diagnostics", "weak-normalization-formula", required))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "three rigorous sources",
        "global n>=31 normalization correction",
        "complete true value and derivative expectations",
        "first-omitted ratio uppers below 0.971",
        "margins above 0.029",
        "removes generalized gauss-laguerre quadrature",
        "one-row certificate",
        "finite-grid-to-collar theorem remain open",
    ):
        if required not in finding:
            issues.append(issue("summary", "weak-main-finding", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in (
        "only to t=10000, f_21",
        "global n>=31 normalization correction",
        "x>=1 monotonicity margins",
        "strictly below one first omitted term",
        "quadrature is not used",
        "no claim about the other finite-grid rows",
        "no row is ready_to_apply",
        "lambda <= 0",
    ):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[FullExpectationIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[FullExpectationIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-required-string", required))
    lowered = text.lower()
    for forbidden in (
        "all-row finite-grid certificate is complete",
        "uniform collar theorem is proved",
        "scaled-curvature monotonicity is proved",
        "cone entry is proved",
        "rh is proved",
        "lambda <= 0 is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-promotion-language", forbidden))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact_path = args.artifact if args.artifact.is_absolute() else REPO_ROOT / args.artifact
    note_path = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = load_json(artifact_path)
    issues: list[FullExpectationIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, expectation_rows, ratio_rows = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, expectation_rows, ratio_rows))
    issues.extend(validate_note(note_path))
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"FULL-EXPECTATION {item.section} [{item.issue}] {item.detail}")
        print(result_line(artifact).replace("0 issues", f"{len(issues)} issues"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
