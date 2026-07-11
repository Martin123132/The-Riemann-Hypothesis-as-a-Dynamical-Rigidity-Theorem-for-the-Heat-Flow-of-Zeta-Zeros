#!/usr/bin/env python3
"""Validate the endpoint x-moment Taylor certificate."""

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

from jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate import (  # noqa: E402
    DEFAULT_ENDPOINT_ROUTE_JSON,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_PHI_TAIL_JSON,
    REPO_ROOT,
    build_artifact,
    result_line,
)


REQUIRED_ROW_IDS = {
    "nlrgexmtc_01_route_scope_import",
    "nlrgexmtc_02_exact_transformed_moments",
    "nlrgexmtc_03_complex_disk_majorant",
    "nlrgexmtc_04_taylor_moment_interval",
    "nlrgexmtc_05_first_panel_remainder_certificate",
    "nlrgexmtc_06_tail_and_splice_handoff",
    "nlrgexmtc_07_acceptance_gate",
}

ALLOWED_ROLES = {
    "scope_reduction",
    "exact_moment_reduction",
    "analytic_remainder_certificate",
    "arb_interval_certificate",
    "first_panel_certificate",
    "compact_interval_handoff",
    "acceptance_gate",
}

ALLOWED_READINESS = {
    "available_scope_reduction",
    "available_exact_reduction",
    "available_analytic_certificate",
    "available_interval_certificate",
    "available_first_panel_certificate",
    "not_ready_to_apply",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Endpoint X-Moment Taylor Certificate",
    "Status: endpoint x-moment Taylor certificate",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.py",
    (
        "validated Jensen-window PF negative-lambda relative-Gaussian endpoint x-moment Taylor certificate: "
        "7 rows, 0 issues, 65 exact-moment rows, 1 certified first panel, 5 open compact panels, "
        "0 ready-to-apply rows"
    ),
    "mu_k=T^(-k/2)*lower_gamma(alpha+1+k/2,1)/Gamma(alpha+1)",
    "disk radius R: 0.2",
    "4R<pi/2 certified: True",
    "finite Phi_30 disk majorant upper: 95.0575638764905004952154227153911180127486223883434865867288",
    "value certified negative: True",
    "derivative certified negative: True",
    "both channels below caps: True",
    "Compose the already separate n>30 Phi-tail source without changing normalization.",
    "Certify the five compact panels 1<=y<=200",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class EndpointXMomentIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> EndpointXMomentIssue:
    return EndpointXMomentIssue(section=section, issue=name, detail=detail)


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("e", "E"))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[EndpointXMomentIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[EndpointXMomentIssue]:
    issues: list[EndpointXMomentIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate"
    ):
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "endpoint x-moment Taylor certificate":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_endpoint_x_panel_route_matrix",
        "source_endpoint_x_panel_route_json",
        "source_phi_tail_grid_certificate",
        "source_phi_tail_grid_json",
        "source_endpoint_parity_repair_matrix",
        "source_interpolation_remainder_route_matrix",
        "source_finite_part_weighted_sum_interval_certificate",
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
        "finite n<=30 first-panel core only",
        "t=10000",
        "f_21",
        "exact transformed gamma moments",
        "complex-disk cauchy remainder",
        "does not compose the n>30 tail",
        "does not certify the five panels",
        "does not prove the full compact-interval integral",
        "uniform collar",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[EndpointXMomentIssue]:
    try:
        recomputed = build_artifact(DEFAULT_ENDPOINT_ROUTE_JSON, DEFAULT_PHI_TAIL_JSON)
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[EndpointXMomentIssue] = []
    for key in ("diagnostics", "matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputation"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[EndpointXMomentIssue], int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[EndpointXMomentIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "not-list", repr(type(rows)))], 0, 0
    seen = {row.get("id") for row in rows if isinstance(row, dict)}
    missing = REQUIRED_ROW_IDS - seen
    extra = seen - REQUIRED_ROW_IDS
    if missing:
        issues.append(issue("matrix_rows", "missing-row-ids", repr(sorted(missing))))
    if extra:
        issues.append(issue("matrix_rows", "extra-row-ids", repr(sorted(extra))))
    ready_to_apply = 0
    first_panel_rows = 0
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
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
        if row.get("readiness") == "available_first_panel_certificate":
            first_panel_rows += 1
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("only", "not ", "does not")):
            issues.append(issue(row_id, "weak-proof-boundary", repr(row.get("proof_boundary"))))
        refs = row.get("source_artifacts", [])
        if not isinstance(refs, list) or not refs:
            issues.append(issue(row_id, "missing-source-artifacts", repr(refs)))
        else:
            for ref in refs:
                issues.extend(validate_ref(row_id, ref))
    if ready_to_apply:
        issues.append(issue("matrix_rows", "unexpected-ready-to-apply", str(ready_to_apply)))
    if first_panel_rows != 1:
        issues.append(issue("matrix_rows", "bad-first-panel-certificate-count", str(first_panel_rows)))
    return issues, len(rows), first_panel_rows


def validate_moment_rows(artifact: dict) -> list[EndpointXMomentIssue]:
    diagnostics = artifact.get("diagnostics", {})
    rows = diagnostics.get("taylor_moment_rows", [])
    issues: list[EndpointXMomentIssue] = []
    if not isinstance(rows, list) or len(rows) != 65:
        return [issue("moment_rows", "bad-row-count", repr(len(rows) if isinstance(rows, list) else rows))]
    degrees = [row.get("degree") for row in rows if isinstance(row, dict)]
    if degrees != list(range(65)):
        issues.append(issue("moment_rows", "bad-degree-ladder", repr(degrees)))
    required_fields = {
        "degree",
        "phi_taylor_coefficient_ball",
        "subtracted_even_ratio_ball",
        "value_core_coefficient_ball",
        "derivative_core_coefficient_ball",
        "transformed_gamma_moment_ball",
        "value_moment_contribution_ball",
        "derivative_moment_contribution_ball",
        "proof_boundary",
    }
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("moment_rows", "bad-row", repr(row)))
            continue
        degree = row.get("degree")
        missing = required_fields - set(row)
        if missing:
            issues.append(issue(f"moment_{degree}", "missing-fields", repr(sorted(missing))))
        should_subtract = isinstance(degree, int) and degree % 2 == 0 and degree <= 40
        has_subtraction = row.get("subtracted_even_ratio_ball") is not None
        if should_subtract != has_subtraction:
            issues.append(
                issue(
                    f"moment_{degree}",
                    "bad-even-ratio-subtraction",
                    f"expected {should_subtract}, got {has_subtraction}",
                )
            )
        if "finite n<=30 core only" not in str(row.get("proof_boundary", "")):
            issues.append(issue(f"moment_{degree}", "weak-proof-boundary", repr(row.get("proof_boundary"))))
    selected = diagnostics.get("selected_taylor_moment_rows", [])
    selected_degrees = {row.get("degree") for row in selected if isinstance(row, dict)}
    if selected_degrees != {0, 1, 2, 20, 40, 41, 42, 48, 56, 64}:
        issues.append(issue("moment_rows", "bad-selected-degrees", repr(sorted(selected_degrees))))
    return issues


def validate_summary(artifact: dict, row_count: int, first_panel_rows: int) -> list[EndpointXMomentIssue]:
    summary = artifact.get("summary", {})
    issues: list[EndpointXMomentIssue] = []
    expected = {
        "matrix_rows": 7,
        "T": 10000,
        "index": 21,
        "alpha": "20.5",
        "y_interval": "0<=y<=1",
        "x_interval": "0<=x<=0.01",
        "precision_bits": 8192,
        "phi_term_count": 30,
        "polynomial_M": 20,
        "taylor_degree": 64,
        "first_tail_degree": 65,
        "taylor_moment_row_count": 65,
        "disk_radius": "0.2",
        "disk_radius_admissible": True,
        "x_to_disk_radius_ratio": "0.05000000000000000000000000000000000000000",
        "value_first_panel_certified_negative": True,
        "derivative_first_panel_certified_negative": True,
        "both_channels_below_caps": True,
        "finite_core_first_panel_certified": True,
        "first_panel_certificate_rows": 1,
        "full_n_tail_composed_here": False,
        "remaining_compact_panel_count": 5,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, expected_value in expected.items():
        if summary.get(key) != expected_value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {expected_value!r}"))
    if row_count != summary.get("matrix_rows"):
        issues.append(issue("summary", "row-count-mismatch", f"{row_count} != {summary.get('matrix_rows')!r}"))
    if first_panel_rows != summary.get("first_panel_certificate_rows"):
        issues.append(
            issue(
                "summary",
                "first-panel-row-count-mismatch",
                f"{first_panel_rows} != {summary.get('first_panel_certificate_rows')!r}",
            )
        )

    numeric_thresholds = (
        ("real_exp4_lower_on_disk", Decimal("0.3"), None),
        ("finite_phi_complex_disk_majorant_upper", None, Decimal("100")),
        ("first_panel_mass_upper", None, Decimal("2e-21")),
        ("value_cauchy_tail_integral_radius_upper", None, Decimal("1e-102")),
        ("derivative_cauchy_tail_integral_radius_upper", None, Decimal("1e-100")),
        ("value_first_panel_abs_upper", None, Decimal("3e-86")),
        ("derivative_first_panel_abs_upper", None, Decimal("6e-85")),
        ("value_abs_to_cap_ratio_upper", None, Decimal("5e-47")),
        ("derivative_abs_to_cap_ratio_upper", None, Decimal("5e-47")),
    )
    for key, lower, upper in numeric_thresholds:
        try:
            value = dec(summary[key])
        except Exception as exc:
            issues.append(issue("summary", f"bad-numeric-{key}", f"{type(exc).__name__}: {exc}"))
            continue
        if lower is not None and value <= lower:
            issues.append(issue("summary", f"too-small-{key}", str(value)))
        if upper is not None and value >= upper:
            issues.append(issue("summary", f"too-large-{key}", str(value)))

    try:
        value_tail = dec(summary["value_cauchy_tail_integral_radius_upper"])
        derivative_tail = dec(summary["derivative_cauchy_tail_integral_radius_upper"])
        value_total = dec(summary["value_first_panel_abs_upper"])
        derivative_total = dec(summary["derivative_first_panel_abs_upper"])
        if value_tail >= value_total:
            issues.append(issue("summary", "value-tail-not-resolved", f"{value_tail} >= {value_total}"))
        if derivative_tail >= derivative_total:
            issues.append(
                issue("summary", "derivative-tail-not-resolved", f"{derivative_tail} >= {derivative_total}")
            )
    except Exception as exc:
        issues.append(issue("summary", "tail-total-comparison-failed", f"{type(exc).__name__}: {exc}"))

    diagnostics = artifact.get("diagnostics", {})
    if diagnostics.get("source_first_panel_mass_upper") != (
        "1.6155560239464103501984383230311864446832165667162E-21"
    ):
        issues.append(
            issue(
                "diagnostics",
                "bad-source-first-panel-mass",
                repr(diagnostics.get("source_first_panel_mass_upper")),
            )
        )
    if diagnostics.get("source_tail_certificate_rows") != 3:
        issues.append(issue("diagnostics", "bad-tail-row-count", repr(diagnostics.get("source_tail_certificate_rows"))))
    if diagnostics.get("source_tail_certified") is not True:
        issues.append(issue("diagnostics", "tail-source-not-certified", repr(diagnostics.get("source_tail_certified"))))
    if diagnostics.get("full_n_tail_composed_here") is not False:
        issues.append(issue("diagnostics", "tail-silently-composed", repr(diagnostics.get("full_n_tail_composed_here"))))
    if len(diagnostics.get("remaining_compact_panels", [])) != 5:
        issues.append(issue("diagnostics", "bad-remaining-panel-list", repr(diagnostics.get("remaining_compact_panels"))))

    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "closes the finite n<=30 first-panel integration problem",
        "degree-64 arb taylor model",
        "exact transformed gamma moments",
        "|z|<=0.2 cauchy majorant",
        "certified negative",
        "endpoint branch obstruction",
        "finite core only",
        "separate n>30 tail",
        "five y>=1 panels",
        "uniform collar theorem remain open",
    ):
        if required not in finding:
            issues.append(issue("summary", "weak-main-finding", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in (
        "finite n<=30",
        "degree 0 through 64",
        "4r<pi/2",
        "full post-degree-64 cauchy tail",
        "not silently absorbed",
        "five panels",
        "no row is ready_to_apply",
        "lambda <= 0",
    ):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[EndpointXMomentIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[EndpointXMomentIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-required-string", required))
    lowered = text.lower()
    for forbidden in (
        "full compact interval is certified",
        "all six compact panels are certified",
        "n>30 tail is composed here",
        "finite-grid interval certificate is complete",
        "uniform collar theorem is proved",
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
    issues: list[EndpointXMomentIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, first_panel_rows = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_moment_rows(artifact))
    issues.extend(validate_summary(artifact, row_count, first_panel_rows))
    issues.extend(validate_note(note_path))

    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"ENDPOINT-X-MOMENT {item.section} [{item.issue}] {item.detail}")
        print(result_line(artifact).replace("0 issues", f"{len(issues)} issues"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
