#!/usr/bin/env python3
"""Validate the Arb degree-40 collar ladder stress artifact."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlrgd40acl_01_degree_ladder_setup",
    "nlrgd40acl_02_all_ladder_levels_bernstein_positive",
    "nlrgd40acl_03_first_omitted_terms_survive",
    "nlrgd40acl_04_degree40_collar_endpoint_stability",
    "nlrgd40acl_05_companion_stencil_bottleneck",
    "nlrgd40acl_06_live_residual_tail_upgrade",
    "nlrgd40acl_07_finite_ladder_promotion_rejected",
    "nlrgd40acl_08_acceptance_gate",
}

ALLOWED_ROLES = {
    "finite_interval_ladder",
    "finite_interval_certificate",
    "finite_diagnostic",
    "live_route",
    "rejected_route",
    "acceptance_gate",
}

REQUIRED_DEGREES = [16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40]

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Degree-40 Arb Collar Ladder Stress",
    "Status: finite theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress`",
    "bound the infinite residual Taylor tail.",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian degree-40 Arb collar ladder stress: 8 rows, 0 issues, 13 degree levels, max degree 40, 0 failed Bernstein rows, 0 ready-to-apply rows",
    "real T interval: [1156, infinity)",
    "u interval: [0, 1/1156]",
    "degree values: [16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40]",
    "all degree levels certified: True",
    "total failed Bernstein rows: 0",
    "weakest stencil lower: degree 16 companion_product = [0.0630023568050961538569390725332 +/- 2.64e-33]",
    "degree=40, M=20: normalizers=4/4, stencils=3/3, failed Bernstein=0",
    "highest ratio c40/c0: [4944469183781380656.77812493750 +/- 3.86e-13] (positive)",
    "B_product endpoint: [37.2033572903651204156102209610 +/- 2.67e-29]",
    "companion_product endpoint: [17.2250223756550133023672457518 +/- 2.10e-29]",
    "weighted_gap_derivative endpoint: [27.4244024269680055723947152492 +/- 3.51e-29]",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class LadderStressIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> LadderStressIssue:
    return LadderStressIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[LadderStressIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[LadderStressIssue]:
    issues: list[LadderStressIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "finite theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_degree16_arb_collar_certificate",
        "source_uniform_remainder_target",
        "source_stencil_remainder_obligations",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "finite",
        "degree 40",
        "fixed k=22",
        "t>=1156",
        "does not bound",
        "infinite residual taylor tail",
        "scaled-curvature",
        "cone entry",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[LadderStressIssue]:
    params = artifact.get("matrix_rows", [{}])[0].get("diagnostics", {}).get("parameters", {})
    try:
        recomputed = build_artifact(
            list(params.get("degree_values", REQUIRED_DEGREES)),
            int(params.get("tail_cutoff_n", 80)),
            int(params.get("precision_bits", 384)),
            int(params.get("tail_start_k", 22)),
            int(params.get("collar_start_T", 1156)),
        )
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[LadderStressIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_diagnostics(artifact: dict) -> list[LadderStressIssue]:
    diagnostics = artifact.get("matrix_rows", [{}])[0].get("diagnostics", {})
    params = diagnostics.get("parameters", {})
    issues: list[LadderStressIssue] = []
    expected_params = {
        "min_taylor_degree": 16,
        "max_taylor_degree": 40,
        "degree_step": 2,
        "degree_values": REQUIRED_DEGREES,
        "tail_cutoff_n": 80,
        "precision_bits": 384,
        "tail_start_k": 22,
        "collar_start_T": 1156,
        "real_interval_u": "[0, 1/1156]",
        "real_interval_T": "[1156, infinity)",
    }
    if params != expected_params:
        issues.append(issue("parameters", "bad-parameters", repr(params)))
    if diagnostics.get("degree_ladder_row_count") != 13:
        issues.append(issue("diagnostics", "bad-degree-ladder-row-count", repr(diagnostics.get("degree_ladder_row_count"))))
    if diagnostics.get("all_degree_ladder_rows_certified") is not True:
        issues.append(issue("diagnostics", "ladder-not-certified", repr(diagnostics.get("all_degree_ladder_rows_certified"))))
    if diagnostics.get("total_failed_bernstein_rows") != 0:
        issues.append(issue("diagnostics", "unexpected-bernstein-failures", repr(diagnostics.get("total_failed_bernstein_rows"))))

    rows = diagnostics.get("degree_ladder_rows", [])
    if [row.get("max_taylor_degree") for row in rows] != REQUIRED_DEGREES:
        issues.append(issue("degree_ladder_rows", "bad-degree-order", repr([row.get("max_taylor_degree") for row in rows])))
    for row in rows:
        degree = row.get("max_taylor_degree")
        expected_ratio_rows = (int(degree) // 2) + 1
        expected_scalars = {
            "continuation_M": int(degree) // 2,
            "ratio_ball_rows": expected_ratio_rows,
            "highest_ratio_degree": degree,
            "positive_normalizer_rows": 4,
            "certified_stencil_rows": 3,
            "failed_bernstein_rows": 0,
            "finite_degree_arb_collar_certified": True,
            "weakest_stencil_name": "companion_product",
        }
        for key, value in expected_scalars.items():
            if row.get(key) != value:
                issues.append(issue(f"degree_{degree}", f"bad-{key}", f"{row.get(key)!r} != {value!r}"))
        if set(row.get("normalizer_min_lowers", {})) != {"F_21", "F_22", "F_23", "F_24"}:
            issues.append(issue(f"degree_{degree}", "bad-normalizer-keys", repr(row.get("normalizer_min_lowers"))))
        if set(row.get("stencil_min_lowers", {})) != {"B_product", "companion_product", "weighted_gap_derivative"}:
            issues.append(issue(f"degree_{degree}", "bad-stencil-keys", repr(row.get("stencil_min_lowers"))))
    degree40 = diagnostics.get("degree40_row", {})
    expected_degree40 = {
        "max_taylor_degree": 40,
        "continuation_M": 20,
        "highest_ratio_degree": 40,
        "highest_ratio_sign": "positive",
        "highest_ratio_to_c0": "[4944469183781380656.77812493750 +/- 3.86e-13]",
        "positive_normalizer_rows": 4,
        "certified_stencil_rows": 3,
        "failed_bernstein_rows": 0,
        "finite_degree_arb_collar_certified": True,
    }
    for key, value in expected_degree40.items():
        if degree40.get(key) != value:
            issues.append(issue("degree40_row", f"bad-{key}", f"{degree40.get(key)!r} != {value!r}"))
    expected_endpoints = {
        "B_product": "[37.2033572903651204156102209610 +/- 2.67e-29]",
        "companion_product": "[17.2250223756550133023672457518 +/- 2.10e-29]",
        "weighted_gap_derivative": "[27.4244024269680055723947152492 +/- 3.51e-29]",
    }
    if degree40.get("stencil_endpoints") != expected_endpoints:
        issues.append(issue("degree40_row", "bad-stencil-endpoints", repr(degree40.get("stencil_endpoints"))))
    weakest = diagnostics.get("weakest_stencil_across_ladder", {})
    if weakest != {
        "degree": 16,
        "name": "companion_product",
        "sample": "[0.0630023568050961538569390725332 +/- 2.64e-33]",
    }:
        issues.append(issue("diagnostics", "bad-weakest-stencil", repr(weakest)))
    note = str(diagnostics.get("proof_boundary_note", "")).lower()
    for required in ("finite-surrogate", "degree 40", "does not estimate", "infinite residual"):
        if required not in note:
            issues.append(issue("diagnostics", "weak-proof-boundary-note", required))
    return issues


def validate_rows(artifact: dict) -> tuple[list[LadderStressIssue], int, int, int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[LadderStressIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0, 0, 0, 0
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
    ladder_rows = 0
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
        if role == "finite_interval_ladder":
            ladder_rows += 1
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
        if not any(marker in boundary for marker in ("finite", "not", "only", "live", "rejected", "proof-hygiene")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    if ready_to_apply:
        issues.append(issue("matrix_rows", "unexpected-ready-to-apply", str(ready_to_apply)))
    return issues, len(rows), ladder_rows, certificate_rows, diagnostic_rows, live_routes + rejected_routes


def validate_summary(
    artifact: dict,
    row_count: int,
    ladder_rows: int,
    certificate_rows: int,
    diagnostic_rows: int,
    route_rows: int,
) -> list[LadderStressIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 8,
        "degree_ladder_rows": 13,
        "min_taylor_degree": 16,
        "max_taylor_degree": 40,
        "collar_start_T": "1156",
        "degree40_positive_normalizer_rows": 4,
        "degree40_certified_stencil_rows": 3,
        "total_failed_bernstein_rows": 0,
        "all_degree_ladder_rows_certified": True,
        "weakest_stencil_across_ladder": {
            "degree": 16,
            "name": "companion_product",
            "sample": "[0.0630023568050961538569390725332 +/- 2.64e-33]",
        },
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    issues: list[LadderStressIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 8:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if ladder_rows != 1:
        issues.append(issue("summary", "bad-ladder-row-count", str(ladder_rows)))
    if certificate_rows != 1:
        issues.append(issue("summary", "bad-certificate-row-count", str(certificate_rows)))
    if diagnostic_rows != 3:
        issues.append(issue("summary", "bad-diagnostic-row-count", str(diagnostic_rows)))
    if route_rows != 2:
        issues.append(issue("summary", "bad-route-row-count", str(route_rows)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "t>=1156",
        "degree 16 through 40",
        "fixed k=22",
        "all four normalizers",
        "weighted-gap stencils",
        "zero bernstein failures",
        "residual-tail bound beyond the last tested degree",
        "fixed-k limitation",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "degree 40", "infinite residual", "k=22", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[LadderStressIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[LadderStressIssue] = []
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


def validate(target_path: Path, note_path: Path) -> tuple[list[LadderStressIssue], dict]:
    artifact = load_json(target_path)
    issues: list[LadderStressIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    issues.extend(validate_diagnostics(artifact))
    row_issues, row_count, ladder_rows, certificate_rows, diagnostic_rows, route_rows = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, ladder_rows, certificate_rows, diagnostic_rows, route_rows))
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
            print(f"JWPF-NEG-LAMBDA-REL-GAUSS-DEG40-ARB-COLLAR-LADDER {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda relative-Gaussian degree-40 Arb collar ladder stress: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('degree_ladder_rows')} degree levels, "
            f"max degree {summary.get('max_taylor_degree')}, "
            f"{summary.get('total_failed_bernstein_rows')} failed Bernstein rows, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
