#!/usr/bin/env python3
"""Validate the degree-16 real-T collar scout for the relative-Gaussian surrogate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlrgd16rt_01_surrogate_setup",
    "nlrgd16rt_02_normalizer_positivity",
    "nlrgd16rt_03_B_product_collar",
    "nlrgd16rt_04_companion_product_collar",
    "nlrgd16rt_05_weighted_gap_derivative_collar",
    "nlrgd16rt_06_integer_scan_to_real_surrogate",
    "nlrgd16rt_07_live_tail_upgrade",
    "nlrgd16rt_08_surrogate_promotion_rejected",
}

ALLOWED_ROLES = {
    "finite_surrogate",
    "finite_surrogate_certificate",
    "finite_diagnostic",
    "live_route",
    "rejected_route",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Degree-16 Real-T Collar Scout",
    "Status: finite theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout`",
    "certifies a real-`T` collar only for the",
    "rationalized degree-16 finite surrogate at fixed `k=22`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian degree-16 real-T collar scout: 8 rows, 0 issues, 4 positive normalizer rows, 3 certified surrogate stencil rows, 0 ready-to-apply rows",
    "real T interval: [1156, infinity)",
    "u interval: (0, 1/1156]",
    "positive normalizer rows: 4",
    "certified product stencil rows: 2",
    "certified derivative stencil rows: 1",
    "certified surrogate stencil rows: 3",
    "root-count failures: 0",
    "real-T surrogate collar certified: True",
    "B: 1.259416498980593755E-4",
    "companion: 8.851471227548859301E-10",
    "weighted gap: 2.518416978813492713E-4",
    "F_21: degree=8, roots=0, endpoint sign=positive, certified=True",
    "F_24: degree=8, roots=0, endpoint sign=positive, certified=True",
    "B: degree=16, zero order=2, stripped degree=14, roots=0, certified=True",
    "companion: degree=32, zero order=3, stripped degree=29, roots=0, certified=True",
    "weighted_gap derivative: degree=30, zero order=1, stripped degree=29, roots=0, certified=True",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class RealTCollarIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> RealTCollarIssue:
    return RealTCollarIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[RealTCollarIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[RealTCollarIssue]:
    issues: list[RealTCollarIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "finite theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_degree16_collar_scan",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "finite",
        "surrogate",
        "fixed k=22",
        "does not prove",
        "arb coefficient-ball",
        "infinite taylor-tail",
        "scaled-curvature",
        "cone entry",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[RealTCollarIssue]:
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
    issues: list[RealTCollarIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_row_dict(section: str, row: dict, expected: dict) -> list[RealTCollarIssue]:
    issues: list[RealTCollarIssue] = []
    for key, value in expected.items():
        if row.get(key) != value:
            issues.append(issue(section, f"bad-{key}", f"{row.get(key)!r} != {value!r}"))
    return issues


def validate_diagnostics(artifact: dict) -> list[RealTCollarIssue]:
    diagnostics = artifact.get("matrix_rows", [{}])[0].get("diagnostics", {})
    issues: list[RealTCollarIssue] = []
    expected_scalars = {
        "certified_positive_normalizer_rows": 4,
        "certified_product_stencil_rows": 2,
        "certified_derivative_stencil_rows": 1,
        "certified_surrogate_stencil_rows": 3,
        "root_count_failures": 0,
        "real_T_surrogate_collar_certified": True,
    }
    for key, value in expected_scalars.items():
        if diagnostics.get(key) != value:
            issues.append(issue("diagnostics", f"bad-{key}", f"{diagnostics.get(key)!r} != {value!r}"))
    expected_params = {
        "coefficient_max_taylor_degree": 16,
        "tail_cutoff_n": 80,
        "precision_bits": 256,
        "tail_start_k": 22,
        "continuation_M": 8,
        "collar_start_T": 1156,
        "real_interval_u": "(0, 1/1156]",
        "real_interval_T": "[1156, infinity)",
    }
    if diagnostics.get("parameters") != expected_params:
        issues.append(issue("parameters", "bad-parameters", repr(diagnostics.get("parameters"))))

    normalizers = diagnostics.get("normalizer_rows", [])
    if [row.get("index") for row in normalizers] != [21, 22, 23, 24]:
        issues.append(issue("normalizers", "bad-index-order", repr([row.get("index") for row in normalizers])))
    normalizer_values = {
        21: "4.881704699409963404E-1",
        22: "4.714791458554841687E-1",
        23: "4.553011796707405967E-1",
        24: "4.396229616097434235E-1",
    }
    for row in normalizers:
        index = row.get("index")
        expected = {
            "polynomial_degree": 8,
            "zero_order_at_u0": 0,
            "stripped_degree": 8,
            "roots_in_open_interval": 0,
            "sign_at_u0_after_stripping": "positive",
            "sign_at_endpoint": "positive",
            "value_at_endpoint": normalizer_values.get(index),
            "certified_positive_on_interval": True,
        }
        issues.extend(validate_row_dict(f"F_{index}", row, expected))

    product_rows = diagnostics.get("product_stencil_rows", [])
    if [row.get("name") for row in product_rows] != ["B", "companion"]:
        issues.append(issue("product_rows", "bad-product-order", repr([row.get("name") for row in product_rows])))
    product_expected = {
        "B": {
            "polynomial_degree": 16,
            "zero_order_at_u0": 2,
            "stripped_degree": 14,
            "roots_in_open_interval": 0,
            "sign_at_u0_after_stripping": "positive",
            "sign_at_endpoint": "positive",
            "value_at_endpoint": "2.799413206146811600E-5",
            "certified_positive_on_interval": True,
        },
        "companion": {
            "polynomial_degree": 32,
            "zero_order_at_u0": 3,
            "stripped_degree": 29,
            "roots_in_open_interval": 0,
            "sign_at_u0_after_stripping": "positive",
            "sign_at_endpoint": "positive",
            "value_at_endpoint": "4.078338730298926972E-11",
            "certified_positive_on_interval": True,
        },
    }
    for row in product_rows:
        issues.extend(validate_row_dict(str(row.get("name")), row, product_expected.get(row.get("name"), {})))

    derivative_rows = diagnostics.get("derivative_stencil_rows", [])
    if len(derivative_rows) != 1:
        issues.append(issue("derivative_rows", "bad-count", str(len(derivative_rows))))
    else:
        expected = {
            "name": "weighted_gap",
            "numerator_degree": 30,
            "zero_order_at_u0": 1,
            "stripped_degree": 29,
            "roots_in_open_interval": 0,
            "sign_at_u0_after_stripping": "positive",
            "sign_at_endpoint": "positive",
            "numerator_value_at_endpoint": "3.045223860527820235E-2",
            "certified_positive_derivative_on_interval": True,
            "certified_positive_stencil_on_interval": True,
        }
        issues.extend(validate_row_dict("weighted_gap_derivative", derivative_rows[0], expected))

    endpoint = diagnostics.get("endpoint_log_stencils", {})
    expected_endpoint = {
        "T": "1156",
        "B": "1.259416498980593755E-4",
        "companion": "8.851471227548859301E-10",
        "weighted_gap": "2.518416978813492713E-4",
    }
    if endpoint != expected_endpoint:
        issues.append(issue("endpoint", "bad-endpoint-log-stencils", repr(endpoint)))
    surrogate_note = str(diagnostics.get("surrogate_note", "")).lower()
    for required in ("midpoint", "finite surrogate", "not an arb coefficient-ball theorem"):
        if required not in surrogate_note:
            issues.append(issue("diagnostics", "weak-surrogate-note", required))
    return issues


def validate_rows(artifact: dict) -> tuple[list[RealTCollarIssue], int, int, int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[RealTCollarIssue] = []
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
        if role == "finite_surrogate":
            surrogate_rows += 1
        if role == "finite_surrogate_certificate":
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
        if not any(marker in boundary for marker in ("finite", "surrogate", "not", "only", "live", "rejected")):
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
) -> list[RealTCollarIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 8,
        "positive_normalizer_rows": 4,
        "certified_product_stencil_rows": 2,
        "certified_derivative_stencil_rows": 1,
        "certified_surrogate_stencil_rows": 3,
        "root_count_failures": 0,
        "real_T_surrogate_collar_start_T": "1156",
        "real_T_surrogate_collar_certified": True,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    issues: list[RealTCollarIssue] = []
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
        "rationalized degree-16 finite surrogate",
        "t=1156",
        "t>=1156",
        "four normalizers are positive",
        "no stripped roots",
        "weighted-gap derivative numerator is positive",
        "not an infinite taylor-tail theorem",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "rationalized", "infinite residual", "k=22", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[RealTCollarIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[RealTCollarIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "analytic zeta-tail theorem is proved",
        "infinite taylor-tail theorem is proved",
        "uniform taylor-tail theorem is proved",
        "scaled-curvature monotonicity is proved",
        "cone entry is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[RealTCollarIssue], dict]:
    artifact = load_json(target_path)
    issues: list[RealTCollarIssue] = []
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
            print(f"JWPF-NEG-LAMBDA-REL-GAUSS-DEG16-REAL-T {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda relative-Gaussian degree-16 real-T collar scout: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('positive_normalizer_rows')} positive normalizer rows, "
            f"{summary.get('certified_surrogate_stencil_rows')} certified surrogate stencil rows, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
