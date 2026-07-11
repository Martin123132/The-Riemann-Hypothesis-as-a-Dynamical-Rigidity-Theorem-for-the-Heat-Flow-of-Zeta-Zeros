#!/usr/bin/env python3
"""Validate the worst-row Arb Chebyshev interpolant-moment scout."""

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

from jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout import (  # noqa: E402
    DEFAULT_FLOATING_SCOUT_JSON,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_QUADRATURE_ROUTE_JSON,
    REPO_ROOT,
    build_artifact,
    result_line,
)


REQUIRED_ROW_IDS = {
    "nlrgwracims_01_floating_scout_import",
    "nlrgwracims_02_arb_interpolant_moment_arithmetic",
    "nlrgwracims_03_arb_degree_ladder",
    "nlrgwracims_04_cauchy_delta_vs_caps",
    "nlrgwracims_05_reference_panel_contributions",
    "nlrgwracims_06_remainder_gap",
    "nlrgwracims_07_acceptance_gate",
}

ALLOWED_ROLES = {
    "scope_reduction",
    "arb_interpolant_diagnostic",
    "rejected_route",
    "acceptance_gate",
}

ALLOWED_READINESS = {
    "diagnostic_only",
    "not_ready_to_apply",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Arb Chebyshev Interpolant-Moment Scout",
    "Status: worst-row Arb Chebyshev interpolant-moment scout",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian worst-row Arb Chebyshev interpolant-moment scout: 7 rows, 0 issues, 4 degrees, 3 Cauchy pairs, 3 cap-safe pairs, 0 ready-to-apply rows",
    "compact interval: 0<=y<=200",
    "degrees: [16, 20, 24, 32]",
    "first cap-safe pair: 16->20",
    "cap-safe pair count: 3",
    "target closing: False",
    "does not bound the",
    "interpolation remainder",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class ArbInterpolantIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> ArbInterpolantIssue:
    return ArbInterpolantIssue(section=section, issue=name, detail=detail)


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("e", "E"))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[ArbInterpolantIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[ArbInterpolantIssue]:
    issues: list[ArbInterpolantIssue] = []
    if (
        artifact.get("kind")
        != "jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout"
    ):
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "worst-row Arb Chebyshev interpolant-moment scout":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_floating_chebyshev_panel_moment_scout",
        "source_floating_chebyshev_panel_moment_scout_json",
        "source_quadrature_remainder_route_matrix",
        "source_quadrature_remainder_route_json",
        "source_compact_interval_integration_scout",
        "source_far_tail_split_certificate",
        "source_finite_part_weighted_sum_interval_certificate",
        "source_intervalization_target",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "arb chebyshev interpolant-moment scout only",
        "interpolant coefficients",
        "does not bound the interpolation remainder",
        "does not prove the true compact integral",
        "quadrature-remainder theorem",
        "all rows",
        "finite-grid interval certificate",
        "uniform collar",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[ArbInterpolantIssue]:
    try:
        recomputed = build_artifact(DEFAULT_FLOATING_SCOUT_JSON, DEFAULT_QUADRATURE_ROUTE_JSON)
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[ArbInterpolantIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[ArbInterpolantIssue], int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[ArbInterpolantIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "not-list", repr(type(rows)))], 0, 0
    seen = {row.get("id") for row in rows if isinstance(row, dict)}
    missing = REQUIRED_ROW_IDS - seen
    extra = seen - REQUIRED_ROW_IDS
    if missing:
        issues.append(issue("matrix_rows", "missing-row-ids", repr(sorted(missing))))
    if extra:
        issues.append(issue("matrix_rows", "extra-row-ids", repr(sorted(extra))))
    diagnostic_rows = 0
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
            issues.append(issue(row_id, "forbidden-ready-to-apply", row_id))
        if row.get("readiness") == "diagnostic_only":
            diagnostic_rows += 1
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("only", "not ", "does not", "rejected")):
            issues.append(issue(row_id, "weak-proof-boundary", repr(row.get("proof_boundary"))))
        refs = row.get("source_artifacts", [])
        if not isinstance(refs, list) or not refs:
            issues.append(issue(row_id, "missing-source-artifacts", repr(refs)))
        else:
            for ref in refs:
                issues.extend(validate_ref(row_id, ref))
    by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    rejected = by_id.get("nlrgwracims_06_remainder_gap", {})
    if "true cancellation-reduced core" not in str(rejected.get("gap", "")).lower():
        issues.append(issue("nlrgwracims_06_remainder_gap", "weak-gap", repr(rejected.get("gap"))))
    gate = by_id.get("nlrgwracims_07_acceptance_gate", {}).get("diagnostics", {})
    if "true cancellation-reduced core" not in str(gate.get("recommended_upgrade", "")).lower():
        issues.append(issue("nlrgwracims_07_acceptance_gate", "weak-upgrade", repr(gate)))
    return issues, len(rows), diagnostic_rows


def validate_summary(artifact: dict, row_count: int, diagnostic_rows: int) -> list[ArbInterpolantIssue]:
    summary = artifact.get("summary", {})
    issues: list[ArbInterpolantIssue] = []
    expected = {
        "matrix_rows": 7,
        "T": 10000,
        "index": 21,
        "u": "1/10000",
        "alpha": "20.5",
        "precision_bits": 1024,
        "polynomial_M": 20,
        "ratio_cutoff_n": 80,
        "phi_term_count": 30,
        "panel_count": 6,
        "degree_count": 4,
        "degrees": [16, 20, 24, 32],
        "cauchy_pair_count": 3,
        "cap_safe_pair_count": 3,
        "first_cap_safe_pair": "16->20",
        "reference_degree": 32,
        "compact_interval": "0<=y<=200",
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, expected_value in expected.items():
        if summary.get(key) != expected_value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {expected_value!r}"))
    if row_count != summary.get("matrix_rows"):
        issues.append(issue("summary", "row-count-mismatch", f"{row_count} != {summary.get('matrix_rows')!r}"))
    if diagnostic_rows != 5:
        issues.append(issue("summary", "bad-diagnostic-row-count", str(diagnostic_rows)))
    cauchy_rows = None
    for row in artifact.get("matrix_rows", []):
        if isinstance(row, dict) and row.get("id") == "nlrgwracims_04_cauchy_delta_vs_caps":
            cauchy_rows = row.get("diagnostics", {}).get("cauchy_rows")
    if not isinstance(cauchy_rows, list) or len(cauchy_rows) != 3:
        issues.append(issue("summary", "bad-cauchy-rows", repr(cauchy_rows)))
        cauchy_rows = []
    by_pair = {row.get("pair"): row for row in cauchy_rows if isinstance(row, dict)}
    for pair in ("16->20", "20->24", "24->32"):
        if by_pair.get(pair, {}).get("both_deltas_below_caps") is not True:
            issues.append(issue("summary", "bad-cap-safe-pair", repr(by_pair.get(pair))))
    if dec(by_pair.get("16->20", {}).get("value_delta_to_cap_ratio_upper", "1")) >= Decimal("0.02"):
        issues.append(issue("summary", "value-16-20-ratio-too-large", repr(by_pair.get("16->20"))))
    if dec(by_pair.get("16->20", {}).get("derivative_delta_to_cap_ratio_upper", "1")) >= Decimal("0.02"):
        issues.append(issue("summary", "derivative-16-20-ratio-too-large", repr(by_pair.get("16->20"))))
    if dec(by_pair.get("20->24", {}).get("value_delta_to_cap_ratio_upper", "1")) >= Decimal("1e-8"):
        issues.append(issue("summary", "value-20-24-ratio-too-large", repr(by_pair.get("20->24"))))
    if dec(by_pair.get("20->24", {}).get("derivative_delta_to_cap_ratio_upper", "1")) >= Decimal("1e-8"):
        issues.append(issue("summary", "derivative-20-24-ratio-too-large", repr(by_pair.get("20->24"))))
    if abs(dec(summary.get("reference_value_interpolant_mid", "0"))) <= Decimal("1e-35"):
        issues.append(issue("summary", "reference-value-too-small", repr(summary.get("reference_value_interpolant_mid"))))
    if abs(dec(summary.get("reference_derivative_interpolant_mid", "0"))) <= Decimal("1e-33"):
        issues.append(
            issue("summary", "reference-derivative-too-small", repr(summary.get("reference_derivative_interpolant_mid")))
        )
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "arb chebyshev interpolant-moment route",
        "upgrades the floating scout",
        "all consecutive cauchy deltas",
        "below the unscaled quadrature caps",
        "does not close",
        "interpolation remainders remain unbounded",
    ):
        if required not in finding:
            issues.append(issue("summary", "weak-main-finding", required))
    upgrade = str(summary.get("recommended_upgrade", "")).lower()
    for required in ("difference", "arb chebyshev interpolant", "true cancellation-reduced core", "panel"):
        if required not in upgrade:
            issues.append(issue("summary", "weak-recommended-upgrade", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "interpolant arithmetic", "panel interpolation remainders", "compact interval", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[ArbInterpolantIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ArbInterpolantIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-required-string", required))
    lowered = text.lower()
    for forbidden in (
        "compact interval-integration certificate is complete",
        "quadrature-remainder theorem is proved",
        "finite-grid interval certificate is complete",
        "uniform collar theorem is proved",
        "rh is proved",
        "lambda <= 0 is proved",
        "arb chebyshev interpolant-moment ladder proves",
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
    issues: list[ArbInterpolantIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, diagnostic_rows = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, diagnostic_rows))
    issues.extend(validate_note(note_path))

    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"ARB-INTERPOLANT-SCOUT {item.section} [{item.issue}] {item.detail}")
        print(result_line(artifact).replace("0 issues", f"{len(issues)} issues"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
