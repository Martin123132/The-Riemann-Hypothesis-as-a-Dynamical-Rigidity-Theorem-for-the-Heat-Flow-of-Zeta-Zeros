#!/usr/bin/env python3
"""Validate the relative-Gaussian first-omitted denominator certificate."""

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

VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402

from jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate import (  # noqa: E402
    DEFAULT_GRID_JSON,
    DEFAULT_INTERVAL_JSON,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_QUADRATURE_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlrgfodc_01_first_omitted_formula",
    "nlrgfodc_02_r21_arb_sign_magnitude",
    "nlrgfodc_03_grid_denominator_lowers",
    "nlrgfodc_04_ratio_cap_absolute_radius_handoff",
    "nlrgfodc_05_denominator_promotion_rejected",
    "nlrgfodc_06_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_reduction",
    "arb_coefficient_certificate",
    "arb_denominator_certificate",
    "exact_budget_handoff",
    "rejected_route",
    "acceptance_gate",
}

ALLOWED_READINESS = {
    "available_exact",
    "available_arb",
    "available_for_intervalization",
    "not_ready_to_apply",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian First-Omitted Denominator Certificate",
    "Status: first-omitted denominator lower certificate",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian first-omitted denominator certificate: 6 rows, 0 issues, 20 denominator rows, 2 ratio-cap rows, 0 ready-to-apply rows",
    "D_value(i,T)=|r_21|*(i+1/2)_21*T^(-18)",
    "D_derivative(i,T)=21*|r_21|*(i+1/2)_21*T^(-19)",
    "minimum value denominator lower: 6.782032247872604818E-22 at T=10000, F_21",
    "minimum derivative denominator lower: 1.424226772053247012E-24 at T=10000, F_21",
    "ratio cap 1.000000000000000000E-6",
    "ratio cap 2.000000000000000000E-3",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.json",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class DenominatorIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> DenominatorIssue:
    return DenominatorIssue(section=section, issue=name, detail=detail)


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("e", "E"))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[DenominatorIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def arb_positive_from_text(text: object) -> bool:
    try:
        value = flint.arb(str(text))
    except Exception:
        return False
    return bool(value > 0 and not value.contains(0))


def validate_top_level(artifact: dict) -> list[DenominatorIssue]:
    issues: list[DenominatorIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "first-omitted denominator lower certificate":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_cancellation_reduced_grid_scout",
        "source_cancellation_reduced_grid_json",
        "source_intervalization_target",
        "source_intervalization_target_json",
        "source_quadrature_ladder_scout",
        "source_quadrature_ladder_json",
        "source_degree16_certificate",
        "source_degree40_residual_budget",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "denominator lower certificate only",
        "divisor side",
        "does not enclose",
        "quadrature-remainder",
        "uniform collar",
        "scaled-curvature",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(
    artifact: dict,
    grid_path: Path,
    interval_path: Path,
    quadrature_path: Path,
) -> list[DenominatorIssue]:
    try:
        recomputed = build_artifact(grid_path, interval_path, quadrature_path)
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[DenominatorIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_diagnostics(artifact: dict) -> list[DenominatorIssue]:
    rows_by_id = {row.get("id"): row for row in artifact.get("matrix_rows", []) if isinstance(row, dict)}
    issues: list[DenominatorIssue] = []
    coeff = rows_by_id.get("nlrgfodc_02_r21_arb_sign_magnitude", {}).get("diagnostics", {})
    expected_coeff = {
        "coefficient_index": 21,
        "coefficient_degree": 42,
        "ratio_cutoff_n": 80,
        "precision_bits": 384,
        "ratio_to_c0_sign": "negative",
    }
    for key, value in expected_coeff.items():
        if coeff.get(key) != value:
            issues.append(issue("coefficient", f"bad-{key}", f"{coeff.get(key)!r} != {value!r}"))
    if dec(coeff.get("absolute_ratio_to_c0_lower", "0")) <= Decimal("3.47E19"):
        issues.append(issue("coefficient", "abs-r21-lower-too-small", repr(coeff.get("absolute_ratio_to_c0_lower"))))
    if not str(coeff.get("proof_boundary", "")).lower().startswith("arb sign"):
        issues.append(issue("coefficient", "weak-proof-boundary", repr(coeff.get("proof_boundary"))))

    diag = rows_by_id.get("nlrgfodc_03_grid_denominator_lowers", {}).get("diagnostics", {})
    denom_summary = diag.get("denominator_summary", {})
    expected_summary = {
        "denominator_rows": 20,
        "t_grid": [1156, 1500, 2000, 5000, 10000],
        "indices": [21, 22, 23, 24],
        "first_j": 21,
        "value_denominator_formula": "|r_21|*(i+1/2)_21*T^(-18)",
        "derivative_denominator_formula": "21*|r_21|*(i+1/2)_21*T^(-19)",
        "minimum_value_denominator_lower": "6.782032247872604818E-22",
        "minimum_value_denominator_location": "T=10000, F_21",
        "minimum_derivative_denominator_lower": "1.424226772053247012E-24",
        "minimum_derivative_denominator_location": "T=10000, F_21",
        "all_denominators_positive": True,
    }
    for key, value in expected_summary.items():
        if denom_summary.get(key) != value:
            issues.append(issue("denominator_summary", f"bad-{key}", f"{denom_summary.get(key)!r} != {value!r}"))
    denom_rows = diag.get("denominator_rows", [])
    if len(denom_rows) != 20:
        issues.append(issue("denominator_rows", "bad-row-count", str(len(denom_rows))))
    seen = {(row.get("T"), row.get("index")) for row in denom_rows if isinstance(row, dict)}
    expected_seen = {(T, index) for T in [1156, 1500, 2000, 5000, 10000] for index in [21, 22, 23, 24]}
    if seen != expected_seen:
        issues.append(issue("denominator_rows", "bad-grid", repr(sorted(seen))))
    for row in denom_rows:
        if not isinstance(row, dict):
            issues.append(issue("denominator_rows", "bad-row", repr(row)))
            continue
        row_id = f"T={row.get('T')},F={row.get('index')}"
        for key in (
            "first_j",
            "rising_factor_exact",
            "value_denominator_ball",
            "value_denominator_lower",
            "derivative_denominator_ball",
            "derivative_denominator_lower",
            "derivative_to_value_denominator_ratio",
            "proof_boundary",
        ):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if row.get("first_j") != 21:
            issues.append(issue(row_id, "bad-first-j", repr(row.get("first_j"))))
        if dec(row.get("value_denominator_lower", "0")) <= 0:
            issues.append(issue(row_id, "value-lower-not-positive", repr(row.get("value_denominator_lower"))))
        if dec(row.get("derivative_denominator_lower", "0")) <= 0:
            issues.append(issue(row_id, "derivative-lower-not-positive", repr(row.get("derivative_denominator_lower"))))
        expected_ratio = f"21/{row.get('T')}"
        if row.get("derivative_to_value_denominator_ratio") != expected_ratio:
            issues.append(issue(row_id, "bad-derivative-value-ratio", repr(row.get("derivative_to_value_denominator_ratio"))))
        boundary = str(row.get("proof_boundary", "")).lower()
        for required in ("denominator only", "numerator", "quadrature"):
            if required not in boundary:
                issues.append(issue(row_id, "weak-row-boundary", required))

    cap_rows = rows_by_id.get("nlrgfodc_04_ratio_cap_absolute_radius_handoff", {}).get("diagnostics", {}).get(
        "ratio_cap_rows", []
    )
    expected_caps = {
        "1.000000000000000000E-6": (
            "6.782032247872604818E-28",
            "1.424226772053247012E-30",
        ),
        "2.000000000000000000E-3": (
            "1.356406449574520964E-24",
            "2.848453544106494024E-27",
        ),
    }
    if len(cap_rows) != 2:
        issues.append(issue("ratio_cap_rows", "bad-row-count", str(len(cap_rows))))
    for row in cap_rows:
        cap = row.get("ratio_error_cap")
        if cap not in expected_caps:
            issues.append(issue("ratio_cap_rows", "bad-cap", repr(cap)))
            continue
        value_cap, derivative_cap = expected_caps[cap]
        if row.get("value_scaled_absolute_radius_cap_from_min_denominator") != value_cap:
            issues.append(issue("ratio_cap_rows", "bad-value-cap", repr(row)))
        if row.get("derivative_scaled_absolute_radius_cap_from_min_denominator") != derivative_cap:
            issues.append(issue("ratio_cap_rows", "bad-derivative-cap", repr(row)))
        if row.get("limiting_value_location") != "T=10000, F_21":
            issues.append(issue("ratio_cap_rows", "bad-value-location", repr(row)))
        if row.get("limiting_derivative_location") != "T=10000, F_21":
            issues.append(issue("ratio_cap_rows", "bad-derivative-location", repr(row)))
    return issues


def validate_rows(artifact: dict) -> tuple[list[DenominatorIssue], int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[DenominatorIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
    ready_to_apply = 0
    available_for_intervalization = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("matrix_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "claim", "proof_boundary", "source_artifacts"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if row.get("role") not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(row.get("role"))))
        if row.get("readiness") not in ALLOWED_READINESS:
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
        if row.get("readiness") == "available_for_intervalization":
            available_for_intervalization += 1
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("only", "not", "rejected")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    if ready_to_apply:
        issues.append(issue("matrix_rows", "unexpected-ready-to-apply", str(ready_to_apply)))
    return issues, len(rows), available_for_intervalization


def validate_summary(artifact: dict, row_count: int, available_for_intervalization: int) -> list[DenominatorIssue]:
    summary = artifact.get("summary", {})
    issues: list[DenominatorIssue] = []
    expected = {
        "matrix_rows": 6,
        "denominator_rows": 20,
        "ratio_cap_rows": 2,
        "certified_denominator_conditions": 2,
        "source_grid_rows": 20,
        "source_interval_obligations": 8,
        "source_quadrature_ladder_rows": 7,
        "first_j": 21,
        "r21_sign_certified": True,
        "minimum_value_denominator_lower": "6.782032247872604818E-22",
        "minimum_value_denominator_location": "T=10000, F_21",
        "minimum_derivative_denominator_lower": "1.424226772053247012E-24",
        "minimum_derivative_denominator_location": "T=10000, F_21",
        "finest_ratio_cap": "1.000000000000000000E-6",
        "finest_value_scaled_absolute_radius_cap": "6.782032247872604818E-28",
        "finest_derivative_scaled_absolute_radius_cap": "1.424226772053247012E-30",
        "all_denominators_positive": True,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 6:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if available_for_intervalization != 2:
        issues.append(issue("summary", "bad-available-for-intervalization-count", str(available_for_intervalization)))
    if dec(summary.get("absolute_r21_lower", "0")) <= Decimal("3.47E19"):
        issues.append(issue("summary", "abs-r21-lower-too-small", repr(summary.get("absolute_r21_lower"))))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "r_21",
        "denominator",
        "t=10000, f_21",
        "ratio-radius",
        "1e-6",
        "numerator",
        "quadrature",
        "grid-to-collar",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "denominator", "residual numerator", "uniform collar", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[DenominatorIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[DenominatorIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "finite-grid interval certificate is complete",
        "uniform residual estimate is proved",
        "residual-tail theorem is proved",
        "quadrature-remainder theorem is proved",
        "first-omitted residual comparison is proved",
        "scaled-curvature monotonicity is proved",
        "cone entry is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(
    target_path: Path,
    note_path: Path,
    grid_path: Path,
    interval_path: Path,
    quadrature_path: Path,
) -> tuple[list[DenominatorIssue], dict]:
    artifact = load_json(target_path)
    issues: list[DenominatorIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact, grid_path, interval_path, quadrature_path))
    issues.extend(validate_diagnostics(artifact))
    row_issues, row_count, available_for_intervalization = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, available_for_intervalization))
    issues.extend(validate_note(note_path))
    return issues, artifact.get("summary", {})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--grid-json", type=Path, default=DEFAULT_GRID_JSON)
    parser.add_argument("--interval-json", type=Path, default=DEFAULT_INTERVAL_JSON)
    parser.add_argument("--quadrature-json", type=Path, default=DEFAULT_QUADRATURE_JSON)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    target = args.target if args.target.is_absolute() else REPO_ROOT / args.target
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    grid_path = args.grid_json if args.grid_json.is_absolute() else REPO_ROOT / args.grid_json
    interval_path = args.interval_json if args.interval_json.is_absolute() else REPO_ROOT / args.interval_json
    quadrature_path = args.quadrature_json if args.quadrature_json.is_absolute() else REPO_ROOT / args.quadrature_json
    issues, summary = validate(target, note, grid_path, interval_path, quadrature_path)
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "summary": summary, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-REL-GAUSS-FIRST-OMITTED-DENOM {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda relative-Gaussian first-omitted denominator certificate: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('denominator_rows')} denominator rows, "
            f"{summary.get('ratio_cap_rows')} ratio-cap rows, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
