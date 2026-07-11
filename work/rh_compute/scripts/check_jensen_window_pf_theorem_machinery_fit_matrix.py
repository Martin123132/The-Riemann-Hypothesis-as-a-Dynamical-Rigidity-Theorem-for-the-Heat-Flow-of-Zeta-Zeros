#!/usr/bin/env python3
"""Validate the Jensen-window PF theorem-machinery fit matrix."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MATRIX = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_theorem_machinery_fit_matrix.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_theorem_machinery_fit_matrix.md"


ALLOWED_VERDICTS = {
    "endpoint_equivalence_only",
    "possible_if_new_kernel_representation",
    "possible_if_preserver_hypotheses_are_proved",
    "route_mismatch",
    "conditional_downstream_only",
    "rejected_circular",
}

REQUIRED_IDS = {
    "tm_01_asw_edrei_pf_sequence_characterization",
    "tm_02_schoenberg_pf_functions_variation_diminishing",
    "tm_03_karlin_basic_composition_cauchy_binet",
    "tm_04_polya_schur_multiplier_preservers",
    "tm_05_gantmacher_krein_sign_regular_matrices",
    "tm_06_laguerre_polya_jensen_limit",
    "tm_07_finite_grid_or_rh_assuming_shortcuts",
}

NON_CLOSING_VERDICTS = {
    "endpoint_equivalence_only",
    "route_mismatch",
    "conditional_downstream_only",
    "rejected_circular",
}

REQUIRED_FEATURES = {
    "outputs every binomially weighted Jensen-window Toeplitz matrix for all d,n",
    "handles binomial weights binom(d,j)",
    "handles all shifts n uniformly",
    "starts from proved noncircular hypotheses about the actual A_k(0)",
    "does not assume Jensen hyperbolicity, Laguerre-Polya membership, RH, or Lambda <= 0",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Theorem Machinery Fit Matrix",
    "Status: theorem-search fit matrix",
    "This is not a proof of Jensen-window",
    "jwpf_06_sign_regular_to_jensen_pf_conversion",
    "work/rh_compute/results/jensen_window_pf_theorem_machinery_fit_matrix.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_theorem_machinery_fit_matrix.py",
    "validated Jensen-window PF theorem machinery fit matrix: 7 rows, 0 issues, 0 ready-to-apply rows",
    "tm_01_asw_edrei_pf_sequence_characterization",
    "tm_02_schoenberg_pf_functions_variation_diminishing",
    "tm_03_karlin_basic_composition_cauchy_binet",
    "tm_04_polya_schur_multiplier_preservers",
    "tm_05_gantmacher_krein_sign_regular_matrices",
    "tm_06_laguerre_polya_jensen_limit",
    "tm_07_finite_grid_or_rh_assuming_shortcuts",
    "no `ready_to_apply` row",
    "degree-3 countermodel",
    "Structural Ansatz Workbench",
    "outputs/jensen_window_pf_structural_ansatz_matrix.md",
    "work/rh_compute/results/jensen_window_pf_structural_ansatz_matrix.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py",
    "validated Jensen-window PF structural ansatz matrix: 6 ansatz rows, 0 issues, 0 ready-to-apply rows",
    "outputs/jensen_window_pf_schur_shape_contract.md",
    "work/rh_compute/results/jensen_window_pf_schur_shape_contract.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_schur_shape_contract.py",
    "finite-band shape obligations",
    "outputs/jensen_window_pf_cauchy_binet_low_degree_scout.md",
    "work/rh_compute/results/jensen_window_pf_cauchy_binet_low_degree_scout.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_cauchy_binet_low_degree_scout.py",
    "`15` formula rows with nonnegative Bernstein coefficients",
    "`0` kernel identities found",
)


@dataclass(frozen=True)
class MatrixIssue:
    row_id: str
    issue: str
    detail: str


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_source(row_id: str, source: str) -> list[MatrixIssue]:
    if source.startswith(("http://", "https://")):
        return []
    if (REPO_ROOT / source).exists():
        return []
    return [MatrixIssue(row_id, "missing-source", source)]


def validate_row(row: dict) -> list[MatrixIssue]:
    issues: list[MatrixIssue] = []
    row_id = str(row.get("id", "<missing-id>"))
    required_fields = (
        "id",
        "theorem_family",
        "primary_sources",
        "fit_to_jwpf06",
        "required_hypotheses",
        "current_evidence",
        "fatal_gap",
        "next_action",
        "verdict",
        "hypotheses_verified",
        "would_close_jwpf06_if_verified",
        "proof_boundary",
    )
    for key in required_fields:
        if key not in row:
            issues.append(MatrixIssue(row_id, "missing-field", key))

    verdict = row.get("verdict")
    if verdict not in ALLOWED_VERDICTS:
        issues.append(MatrixIssue(row_id, "bad-verdict", repr(verdict)))
    if row.get("hypotheses_verified") is not False:
        issues.append(MatrixIssue(row_id, "hypotheses-should-not-be-verified", repr(row.get("hypotheses_verified"))))
    if verdict in NON_CLOSING_VERDICTS and row.get("would_close_jwpf06_if_verified") is True:
        issues.append(MatrixIssue(row_id, "nonclosing-verdict-closes", str(verdict)))
    if row.get("would_close_jwpf06_if_verified") is True and verdict not in {
        "possible_if_new_kernel_representation",
        "possible_if_preserver_hypotheses_are_proved",
    }:
        issues.append(MatrixIssue(row_id, "bad-closing-verdict", str(verdict)))

    for key in ("primary_sources", "required_hypotheses", "current_evidence"):
        value = row.get(key)
        if not isinstance(value, list) or not value:
            issues.append(MatrixIssue(row_id, f"bad-{key}", repr(value)))
    for source in row.get("primary_sources", []):
        if isinstance(source, str):
            issues.extend(validate_source(row_id, source))
        else:
            issues.append(MatrixIssue(row_id, "bad-source", repr(source)))

    if not str(row.get("fatal_gap", "")).strip():
        issues.append(MatrixIssue(row_id, "missing-fatal-gap", row_id))
    if not str(row.get("next_action", "")).strip():
        issues.append(MatrixIssue(row_id, "missing-next-action", row_id))
    boundary = str(row.get("proof_boundary", "")).lower()
    if not any(marker in boundary for marker in ("not ", "only", "missing", "rejected", "conditional", "candidate")):
        issues.append(MatrixIssue(row_id, "weak-proof-boundary", row.get("proof_boundary", "")))

    combined = " ".join(str(row.get(key, "")) for key in ("fit_to_jwpf06", "fatal_gap", "proof_boundary")).lower()
    for forbidden in ("therefore rh", "proves lambda <= 0", "bridge is proved", "ready_to_apply"):
        if forbidden in combined:
            issues.append(MatrixIssue(row_id, "forbidden-overclaim", forbidden))
    return issues


def validate_note(path: Path) -> list[MatrixIssue]:
    if not path.exists():
        return [MatrixIssue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[MatrixIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(MatrixIssue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "the bridge is proved"):
        if forbidden in lowered:
            issues.append(MatrixIssue("note", "forbidden-text", forbidden))
    return issues


def validate(matrix_path: Path, note_path: Path) -> tuple[list[MatrixIssue], int, int]:
    matrix = load_json(matrix_path)
    issues: list[MatrixIssue] = []
    if matrix.get("kind") != "jensen_window_pf_theorem_machinery_fit_matrix":
        issues.append(MatrixIssue("<matrix>", "bad-kind", repr(matrix.get("kind"))))
    if matrix.get("target_obligation") != "jwpf_06_sign_regular_to_jensen_pf_conversion":
        issues.append(MatrixIssue("<matrix>", "bad-target-obligation", repr(matrix.get("target_obligation"))))

    boundary = str(matrix.get("proof_boundary", "")).lower()
    if "not a proof" not in boundary or "lambda <= 0" not in boundary:
        issues.append(MatrixIssue("<matrix>", "weak-proof-boundary", matrix.get("proof_boundary", "")))

    features = set(matrix.get("required_bridge_features", []))
    missing_features = REQUIRED_FEATURES - features
    if missing_features:
        issues.append(MatrixIssue("<matrix>", "missing-required-features", repr(sorted(missing_features))))

    allowed = set(matrix.get("allowed_verdicts", []))
    if not ALLOWED_VERDICTS.issubset(allowed):
        issues.append(MatrixIssue("<matrix>", "missing-allowed-verdicts", repr(sorted(ALLOWED_VERDICTS - allowed))))

    rows = matrix.get("rows", [])
    if not isinstance(rows, list) or not rows:
        issues.append(MatrixIssue("<matrix>", "missing-rows", "rows must be a nonempty list"))
        rows = []

    seen: set[str] = set()
    ready_count = 0
    closing_candidate_count = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(MatrixIssue("<matrix>", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        if row_id in seen:
            issues.append(MatrixIssue(row_id, "duplicate-id", row_id))
        seen.add(row_id)
        if row.get("verdict") == "ready_to_apply":
            ready_count += 1
        if row.get("would_close_jwpf06_if_verified") is True:
            closing_candidate_count += 1
        issues.extend(validate_row(row))

    for missing in sorted(REQUIRED_IDS - seen):
        issues.append(MatrixIssue(missing, "missing-required-row", missing))
    if ready_count != 0:
        issues.append(MatrixIssue("<matrix>", "ready-to-apply-row-present", str(ready_count)))
    if closing_candidate_count < 2:
        issues.append(MatrixIssue("<matrix>", "too-few-structural-candidates", str(closing_candidate_count)))
    if "tm_07_finite_grid_or_rh_assuming_shortcuts" in seen:
        shortcut = next(row for row in rows if isinstance(row, dict) and row.get("id") == "tm_07_finite_grid_or_rh_assuming_shortcuts")
        if shortcut.get("verdict") != "rejected_circular":
            issues.append(MatrixIssue("tm_07_finite_grid_or_rh_assuming_shortcuts", "shortcut-not-rejected", str(shortcut.get("verdict"))))

    issues.extend(validate_note(note_path))
    return issues, len(rows), ready_count


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, row_count, ready_count = validate(args.matrix, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "rows": row_count,
                    "ready_to_apply_rows": ready_count,
                    "issues": [asdict(issue) for issue in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for issue in issues:
            print(f"JWPF-THEOREM-MACHINERY {issue.row_id} [{issue.issue}] {issue.detail}")
        print(
            "validated Jensen-window PF theorem machinery fit matrix: "
            f"{row_count} rows, {len(issues)} issues, {ready_count} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
