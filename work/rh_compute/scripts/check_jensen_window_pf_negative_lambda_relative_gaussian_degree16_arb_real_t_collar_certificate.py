#!/usr/bin/env python3
"""Validate the Arb finite-degree real-T collar certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlrgd16arb_01_coefficient_ball_setup",
    "nlrgd16arb_02_normalizer_bernstein_certificate",
    "nlrgd16arb_03_B_product_bernstein_certificate",
    "nlrgd16arb_04_companion_product_bernstein_certificate",
    "nlrgd16arb_05_weighted_gap_derivative_bernstein_certificate",
    "nlrgd16arb_06_real_t_interval_surrogate_upgrade",
    "nlrgd16arb_07_live_infinite_tail_upgrade",
    "nlrgd16arb_08_interval_surrogate_promotion_rejected",
}

ALLOWED_ROLES = {
    "finite_interval_surrogate",
    "finite_interval_certificate",
    "finite_diagnostic",
    "live_route",
    "rejected_route",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Degree-16 Arb Real-T Collar Certificate",
    "Status: finite theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate`",
    "degree-16 finite surrogate at fixed `k=22`, using Arb coefficient-ratio",
    "It does not bound the infinite residual Taylor tail.",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian degree-16 Arb real-T collar certificate: 8 rows, 0 issues, 4 positive normalizer rows, 3 certified stencil rows, 0 ready-to-apply rows",
    "real T interval: [1156, infinity)",
    "u interval: [0, 1/1156]",
    "ratio ball rows: 9",
    "positive normalizer rows: 4",
    "certified stencil rows: 3",
    "failed Bernstein rows: 0",
    "Bernstein subdivisions: 1",
    "finite-degree Arb collar certified: True",
    "c16/c0: [-237753170.949415650656506482065 +/- 4.12e-22]",
    "F_21: degree=8, positive Bernstein coefficients=9/9, certified=True",
    "F_24: degree=8, positive Bernstein coefficients=9/9, certified=True",
    "B_product: degree=16, zero order=2, Bernstein degree=14, positive coefficients=15/15, certified=True",
    "companion_product: degree=32, zero order=3, Bernstein degree=29, positive coefficients=30/30, certified=True",
    "weighted_gap_derivative: degree=31, zero order=1, Bernstein degree=30, positive coefficients=31/31, certified=True",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class ArbRealTCollarIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> ArbRealTCollarIssue:
    return ArbRealTCollarIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[ArbRealTCollarIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[ArbRealTCollarIssue]:
    issues: list[ArbRealTCollarIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "finite theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_real_t_collar_scout",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "finite",
        "arb coefficient-ratio balls",
        "fixed k=22",
        "does not bound",
        "infinite residual taylor tail",
        "scaled-curvature",
        "cone entry",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[ArbRealTCollarIssue]:
    params = artifact.get("matrix_rows", [{}])[0].get("diagnostics", {}).get("parameters", {})
    try:
        recomputed = build_artifact(
            int(params.get("coefficient_max_taylor_degree", 16)),
            int(params.get("tail_cutoff_n", 80)),
            int(params.get("precision_bits", 256)),
            int(params.get("tail_start_k", 22)),
            int(params.get("continuation_M", 8)),
            int(params.get("collar_start_T", 1156)),
        )
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[ArbRealTCollarIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_dict(section: str, row: dict, expected: dict) -> list[ArbRealTCollarIssue]:
    issues: list[ArbRealTCollarIssue] = []
    for key, value in expected.items():
        if row.get(key) != value:
            issues.append(issue(section, f"bad-{key}", f"{row.get(key)!r} != {value!r}"))
    return issues


def validate_diagnostics(artifact: dict) -> list[ArbRealTCollarIssue]:
    diagnostics = artifact.get("matrix_rows", [{}])[0].get("diagnostics", {})
    issues: list[ArbRealTCollarIssue] = []
    expected_params = {
        "coefficient_max_taylor_degree": 16,
        "tail_cutoff_n": 80,
        "precision_bits": 256,
        "tail_start_k": 22,
        "continuation_M": 8,
        "collar_start_T": 1156,
        "real_interval_u": "[0, 1/1156]",
        "real_interval_T": "[1156, infinity)",
    }
    if diagnostics.get("parameters") != expected_params:
        issues.append(issue("parameters", "bad-parameters", repr(diagnostics.get("parameters"))))
    expected_scalars = {
        "ratio_ball_rows_count": 9,
        "positive_normalizer_rows": 4,
        "certified_stencil_rows": 3,
        "failed_bernstein_rows": 0,
        "bernstein_subdivisions": 1,
        "finite_degree_arb_collar_certified": True,
    }
    for key, value in expected_scalars.items():
        if diagnostics.get(key) != value:
            issues.append(issue("diagnostics", f"bad-{key}", f"{diagnostics.get(key)!r} != {value!r}"))

    ratio_rows = diagnostics.get("ratio_ball_rows", [])
    if [row.get("degree") for row in ratio_rows] != [0, 2, 4, 6, 8, 10, 12, 14, 16]:
        issues.append(issue("ratio_rows", "bad-degree-order", repr([row.get("degree") for row in ratio_rows])))
    expected_last = {
        12: {
            "ratio_to_c0": "[-3304359.71281471153631522353332 +/- 4.76e-24]",
            "radius": "[1.674061151e-65 +/- 1.96e-75]",
            "sign": "negative",
        },
        14: {
            "ratio_to_c0": "[30122498.1078104844171418342933 +/- 4.05e-23]",
            "radius": "[9.034890952e-64 +/- 4.91e-74]",
            "sign": "positive",
        },
        16: {
            "ratio_to_c0": "[-237753170.949415650656506482065 +/- 4.12e-22]",
            "radius": "[4.114298378e-62 +/- 4.67e-74]",
            "sign": "negative",
        },
    }
    for row in ratio_rows[-3:]:
        issues.extend(validate_dict(f"c{row.get('degree')}", row, expected_last.get(row.get("degree"), {})))

    normalizers = diagnostics.get("normalizer_rows", [])
    if [row.get("name") for row in normalizers] != ["F_21", "F_22", "F_23", "F_24"]:
        issues.append(issue("normalizers", "bad-order", repr([row.get("name") for row in normalizers])))
    for row in normalizers:
        expected = {
            "polynomial_degree": 8,
            "stripped_zero_order": 0,
            "bernstein_degree": 8,
            "bernstein_subdivisions": 1,
            "positive_bernstein_coefficients": 9,
            "total_bernstein_coefficients": 9,
            "certified_positive_on_interval": True,
        }
        issues.extend(validate_dict(str(row.get("name")), row, expected))

    stencils = diagnostics.get("stencil_rows", [])
    if [row.get("name") for row in stencils] != ["B_product", "companion_product", "weighted_gap_derivative"]:
        issues.append(issue("stencil_rows", "bad-order", repr([row.get("name") for row in stencils])))
    expected_stencils = {
        "B_product": {
            "polynomial_degree": 16,
            "stripped_zero_order": 2,
            "bernstein_degree": 14,
            "positive_bernstein_coefficients": 15,
            "total_bernstein_coefficients": 15,
            "endpoint_value": "[37.4095664624940562618495050575 +/- 1.08e-29]",
            "certified_positive_on_interval": True,
        },
        "companion_product": {
            "polynomial_degree": 32,
            "stripped_zero_order": 3,
            "bernstein_degree": 29,
            "positive_bernstein_coefficients": 30,
            "total_bernstein_coefficients": 30,
            "endpoint_value": "[0.0630023568050961538569390725332 +/- 2.64e-33]",
            "certified_positive_on_interval": True,
        },
        "weighted_gap_derivative": {
            "polynomial_degree": 31,
            "stripped_zero_order": 1,
            "bernstein_degree": 30,
            "positive_bernstein_coefficients": 31,
            "total_bernstein_coefficients": 31,
            "endpoint_value": "[35.2027878277016019205889550443 +/- 2.30e-29]",
            "certified_positive_on_interval": True,
        },
    }
    for row in stencils:
        issues.extend(validate_dict(str(row.get("name")), row, expected_stencils.get(row.get("name"), {})))
    note = str(diagnostics.get("proof_boundary_note", "")).lower()
    for required in ("degree-16 finite surrogate", "arb coefficient balls", "does not bound", "infinite residual"):
        if required not in note:
            issues.append(issue("diagnostics", "weak-proof-boundary-note", required))
    return issues


def validate_rows(artifact: dict) -> tuple[list[ArbRealTCollarIssue], int, int, int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[ArbRealTCollarIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0, 0, 0, 0
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
    surrogate_rows = 0
    certificate_rows = 0
    diagnostic_rows = 0
    live_routes = 0
    rejected_routes = 0
    ready_to_apply = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("matrix_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "claim", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        role = row.get("role")
        if role not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(role)))
        if role == "finite_interval_surrogate":
            surrogate_rows += 1
        if role == "finite_interval_certificate":
            certificate_rows += 1
        if role == "finite_diagnostic":
            diagnostic_rows += 1
        if role == "live_route":
            live_routes += 1
        if role == "rejected_route":
            rejected_routes += 1
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
        elif row.get("readiness") != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("finite", "not", "only", "live", "rejected", "surrogate")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    if ready_to_apply:
        issues.append(issue("matrix_rows", "unexpected-ready-to-apply", str(ready_to_apply)))
    return issues, len(rows), surrogate_rows, certificate_rows, diagnostic_rows, live_routes + rejected_routes


def validate_summary(
    artifact: dict,
    row_count: int,
    surrogate_rows: int,
    certificate_rows: int,
    diagnostic_rows: int,
    route_rows: int,
) -> list[ArbRealTCollarIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 8,
        "ratio_ball_rows": 9,
        "positive_normalizer_rows": 4,
        "certified_stencil_rows": 3,
        "failed_bernstein_rows": 0,
        "bernstein_subdivisions": 1,
        "real_T_arb_finite_degree_collar_start_T": "1156",
        "finite_degree_arb_collar_certified": True,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    issues: list[ArbRealTCollarIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 8:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if surrogate_rows != 1:
        issues.append(issue("summary", "bad-surrogate-row-count", str(surrogate_rows)))
    if certificate_rows != 4:
        issues.append(issue("summary", "bad-certificate-row-count", str(certificate_rows)))
    if diagnostic_rows != 1:
        issues.append(issue("summary", "bad-diagnostic-row-count", str(diagnostic_rows)))
    if route_rows != 2:
        issues.append(issue("summary", "bad-route-row-count", str(route_rows)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "arb coefficient-ratio balls",
        "finite surrogate",
        "bernstein coefficients certify",
        "f_21..f_24",
        "stripped companion product",
        "weighted-gap derivative numerator",
        "0<=u<=1/1156",
        "infinite residual taylor-tail theorem open",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "degree 16", "infinite residual", "k=22", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[ArbRealTCollarIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ArbRealTCollarIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "infinite taylor-tail theorem is proved",
        "uniform taylor-tail theorem is proved",
        "scaled-curvature monotonicity is proved",
        "cone entry is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[ArbRealTCollarIssue], dict]:
    artifact = load_json(target_path)
    issues: list[ArbRealTCollarIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    issues.extend(validate_diagnostics(artifact))
    row_issues, row_count, surrogate_rows, certificate_rows, diagnostic_rows, route_rows = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, surrogate_rows, certificate_rows, diagnostic_rows, route_rows))
    issues.extend(validate_note(note_path))
    return issues, artifact.get("summary", {})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    target = args.target if args.target.is_absolute() else REPO_ROOT / args.target
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    issues, summary = validate(target, note)
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "summary": summary, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-REL-GAUSS-DEG16-ARB-REAL-T {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda relative-Gaussian degree-16 Arb real-T collar certificate: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('positive_normalizer_rows')} positive normalizer rows, "
            f"{summary.get('certified_stencil_rows')} certified stencil rows, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
