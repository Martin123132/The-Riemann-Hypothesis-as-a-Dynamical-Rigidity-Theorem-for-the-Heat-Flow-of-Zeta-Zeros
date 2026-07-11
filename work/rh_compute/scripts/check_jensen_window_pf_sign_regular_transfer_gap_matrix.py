#!/usr/bin/env python3
"""Validate the Jensen-window PF sign-regular transfer gap matrix."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MATRIX = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_sign_regular_transfer_gap_matrix.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_sign_regular_transfer_gap_matrix.md"


EXPECTED_ROW_IDS = (
    "srgt_01_degree2_exact_contact",
    "srgt_02_degree3_new_obligation",
    "srgt_03_finite_reshaped_hankel_countermodel",
    "srgt_04_selected_toeplitz_countermodel",
    "srgt_05_all_order_antecedent_requirement",
    "srgt_06_xi_phi_specific_transfer_requirement",
    "srgt_07_binomial_shift_uniformity_requirement",
    "srgt_08_acceptable_theorem_shape",
    "srgt_09_forbidden_transfer_shortcuts",
)

ALLOWED_ROLES = {
    "exact_contact",
    "exact_obstruction",
    "rejected_shortcut",
    "open_requirement",
    "live_contract",
}

REQUIRED_SOURCE_PATHS = (
    "outputs/jensen_hankel_bridge_algebra.md",
    "outputs/jensen_window_pf_obligation_algebra.md",
    "outputs/signed_hankel_jensen_bridge_target.md",
    "outputs/jensen_window_pf_bridge_obligations.md",
    "work/rh_compute/results/jensen_hankel_bridge_algebra.json",
    "work/rh_compute/results/jensen_window_pf_obligation_algebra.json",
)

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Sign-Regular Transfer Gap Matrix",
    "Status: finite theorem-search diagnostic",
    "Artifact kind: `jensen_window_pf_sign_regular_transfer_gap_matrix`.",
    "work/rh_compute/results/jensen_window_pf_sign_regular_transfer_gap_matrix.json",
    "python work/rh_compute/scripts/jensen_window_pf_sign_regular_transfer_gap_matrix.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_sign_regular_transfer_gap_matrix.py",
    "validated Jensen-window PF sign-regular transfer gap matrix: 9 transfer rows, 2 countermodel gates, 3 open requirements, 3 rejected shortcuts, 0 ready-to-apply rows, 0 issues",
    "A_0..A_4 = 1, 33/40, 429/640, 4719/12800, 4719/35840",
    "degree-3 Jensen discriminant: -2476526481/3276800000",
    "d=3 first negative contiguous size: 8",
    "d=4 quartic discriminant: -668519580649275927/359661568000000000",
    "d=4 first negative contiguous size: 6",
    "1. all-order, all-shift reshaped-Hankel sign consistency for the actual A_k(0)",
    "2. noncircular Xi/Phi-specific structure absent from arbitrary positive sequences",
    "3. binomial-weight and shift uniformity for every Jensen-window Toeplitz minor",
    "4. no endpoint PF, Jensen hyperbolicity, Laguerre-Polya, RH, or Lambda <= 0 assumptions",
    "srgt_01_degree2_exact_contact",
    "srgt_09_forbidden_transfer_shortcuts",
)

FORBIDDEN_TEXT = (
    "therefore rh",
    "we have proved lambda <= 0",
    "proves lambda <= 0",
    "bridge is proved",
    "sign-regular-to-jensen conversion is proved",
    "jensen-window pf-infinity is proved",
)


@dataclass(frozen=True)
class TransferGapIssue:
    row_id: str
    issue: str
    detail: str


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def require(
    condition: bool,
    issues: list[TransferGapIssue],
    row_id: str,
    issue: str,
    detail: str,
) -> None:
    if not condition:
        issues.append(TransferGapIssue(row_id, issue, detail))


def validate_paths(matrix: dict) -> list[TransferGapIssue]:
    issues: list[TransferGapIssue] = []
    for path in REQUIRED_SOURCE_PATHS:
        require((REPO_ROOT / path).exists(), issues, "<paths>", "missing-required-source", path)
    for key in ("bridge_json", "obligation_json", "generator", "checker"):
        value = matrix.get(key)
        require(isinstance(value, str) and (REPO_ROOT / value).exists(), issues, "<paths>", f"missing-{key}", repr(value))
    return issues


def validate_countermodel(matrix: dict) -> list[TransferGapIssue]:
    issues: list[TransferGapIssue] = []
    data = matrix.get("countermodel_data", {})
    bridge = load_json(REPO_ROOT / str(matrix.get("bridge_json", "")))
    obligation = load_json(REPO_ROOT / str(matrix.get("obligation_json", "")))
    bridge_counter = bridge.get("finite_countermodel", {})
    obligation_counter = obligation.get("finite_countermodel", {})

    expected = {
        "sequence_A0_to_A4": ["1", "33/40", "429/640", "4719/12800", "4719/35840"],
        "finite_reshaped_signs_pass": True,
        "degree3_jensen_discriminant": "-2476526481/3276800000",
        "degree3_jensen_hyperbolicity_breaks": True,
        "d3_selected_toeplitz_minors_positive": True,
        "d3_first_negative_contiguous_size": 8,
        "d3_first_negative_contiguous_determinant": "-435846079534239/104857600000000",
        "d4_selected_toeplitz_minors_positive": True,
        "d4_quartic_discriminant": "-668519580649275927/359661568000000000",
        "d4_first_negative_contiguous_size": 6,
        "d4_first_negative_contiguous_determinant": "-229760849637/28672000000",
    }
    for key, value in expected.items():
        require(data.get(key) == value, issues, "countermodel", f"bad-{key}", repr(data.get(key)))

    require(
        data.get("sequence_A0_to_A4") == bridge_counter.get("sequence_A0_to_A4"),
        issues,
        "countermodel",
        "bridge-sequence-mismatch",
        repr((data.get("sequence_A0_to_A4"), bridge_counter.get("sequence_A0_to_A4"))),
    )
    require(
        data.get("degree3_jensen_discriminant") == bridge_counter.get("degree3_jensen_discriminant"),
        issues,
        "countermodel",
        "bridge-discriminant-mismatch",
        repr((data.get("degree3_jensen_discriminant"), bridge_counter.get("degree3_jensen_discriminant"))),
    )
    require(
        data.get("d3_first_negative_contiguous_size")
        == obligation_counter.get("d3_first_negative_contiguous_toeplitz_minor", {}).get("size"),
        issues,
        "countermodel",
        "obligation-d3-size-mismatch",
        repr(obligation_counter.get("d3_first_negative_contiguous_toeplitz_minor")),
    )
    require(
        data.get("d4_first_negative_contiguous_determinant")
        == obligation_counter.get("d4_first_negative_contiguous_toeplitz_minor", {}).get("determinant"),
        issues,
        "countermodel",
        "obligation-d4-det-mismatch",
        repr(obligation_counter.get("d4_first_negative_contiguous_toeplitz_minor")),
    )
    return issues


def validate_rows(matrix: dict) -> tuple[list[TransferGapIssue], int, int, int, int]:
    issues: list[TransferGapIssue] = []
    rows = matrix.get("transfer_rows", [])
    require(isinstance(rows, list), issues, "<rows>", "bad-rows", repr(type(rows)))
    if not isinstance(rows, list):
        return issues, 0, 0, 0, 0

    ids = [row.get("id") for row in rows if isinstance(row, dict)]
    require(tuple(ids) == EXPECTED_ROW_IDS, issues, "<rows>", "bad-row-ids", repr(ids))

    ready_count = 0
    open_requirement_count = 0
    rejected_shortcut_count = 0
    countermodel_gate_count = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(TransferGapIssue("<rows>", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        role = row.get("role")
        require(role in ALLOWED_ROLES, issues, row_id, "bad-role", repr(role))
        if row.get("readiness") == "ready_to_apply":
            ready_count += 1
        require(row.get("readiness") == "not_ready_to_apply", issues, row_id, "bad-readiness", repr(row.get("readiness")))
        for key in ("claim", "acceptance_test", "proof_boundary"):
            require(bool(str(row.get(key, "")).strip()), issues, row_id, f"missing-{key}", repr(row.get(key)))
        if role == "open_requirement":
            open_requirement_count += 1
        if role == "rejected_shortcut":
            rejected_shortcut_count += 1
        if "countermodel" in row_id:
            countermodel_gate_count += 1
        for source in row.get("source_artifacts", []):
            require((REPO_ROOT / source).exists(), issues, row_id, "missing-source-artifact", str(source))

        combined = " ".join(str(row.get(key, "")) for key in ("claim", "acceptance_test", "proof_boundary")).lower()
        for forbidden in FORBIDDEN_TEXT:
            if forbidden in combined:
                issues.append(TransferGapIssue(row_id, "forbidden-overclaim", forbidden))
    return issues, ready_count, open_requirement_count, rejected_shortcut_count, countermodel_gate_count


def validate_note(path: Path) -> list[TransferGapIssue]:
    if not path.exists():
        return [TransferGapIssue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[TransferGapIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        require(required in text, issues, "note", "missing-text", required)
    lowered = text.lower()
    for forbidden in FORBIDDEN_TEXT:
        if forbidden in lowered:
            issues.append(TransferGapIssue("note", "forbidden-text", forbidden))
    return issues


def validate(matrix_path: Path, note_path: Path) -> tuple[list[TransferGapIssue], dict]:
    matrix = load_json(matrix_path)
    issues: list[TransferGapIssue] = []
    require(
        matrix.get("kind") == "jensen_window_pf_sign_regular_transfer_gap_matrix",
        issues,
        "<matrix>",
        "bad-kind",
        repr(matrix.get("kind")),
    )
    require(
        matrix.get("status") == "finite_theorem_search_diagnostic",
        issues,
        "<matrix>",
        "bad-status",
        repr(matrix.get("status")),
    )
    boundary = str(matrix.get("proof_boundary", "")).lower()
    for required in ("finite theorem-search diagnostic", "does not prove", "lambda <= 0"):
        require(required in boundary, issues, "<matrix>", "weak-proof-boundary", matrix.get("proof_boundary", ""))

    issues.extend(validate_paths(matrix))
    issues.extend(validate_countermodel(matrix))
    row_issues, ready_count, open_count, rejected_count, countermodel_count = validate_rows(matrix)
    issues.extend(row_issues)

    summary = matrix.get("summary", {})
    expected_summary = {
        "transfer_rows": 9,
        "countermodel_gates": 2,
        "open_requirement_rows": 3,
        "rejected_shortcut_rows": 3,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, value in expected_summary.items():
        require(summary.get(key) == value, issues, "summary", f"bad-{key}", repr(summary.get(key)))
    require(summary.get("transfer_rows") == len(matrix.get("transfer_rows", [])), issues, "summary", "row-count-mismatch", repr(summary))
    require(summary.get("ready_to_apply_rows") == ready_count, issues, "summary", "ready-count-mismatch", repr(ready_count))
    require(summary.get("open_requirement_rows") == open_count, issues, "summary", "open-count-mismatch", repr(open_count))
    require(summary.get("rejected_shortcut_rows") == rejected_count, issues, "summary", "rejected-count-mismatch", repr(rejected_count))
    require(summary.get("countermodel_gates") == countermodel_count, issues, "summary", "countermodel-count-mismatch", repr(countermodel_count))

    invariants = "\n".join(str(item) for item in matrix.get("invariants", []))
    for required in (
        "No row is ready_to_apply.",
        "Degree-2 contact is not promoted to all-degree transfer.",
        "Finite reshaped-Hankel signs are not promoted to all-order sign consistency.",
        "Selected low-order Toeplitz minors are not promoted to Jensen-window PF-infinity.",
        "Endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, and Lambda <= 0 are forbidden as inputs.",
    ):
        require(required in invariants, issues, "invariants", "missing-invariant", required)

    issues.extend(validate_note(note_path))
    return issues, matrix


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, matrix = validate(args.matrix, args.note)
    summary = matrix.get("summary", {})
    if args.json:
        print(json.dumps({"ok": not issues, "issues": [asdict(issue) for issue in issues]}, indent=2, sort_keys=True))
    else:
        for issue in issues:
            print(f"JENSEN-WINDOW-PF-SIGN-REGULAR-TRANSFER-GAP {issue.row_id} [{issue.issue}] {issue.detail}")
        print(
            "validated Jensen-window PF sign-regular transfer gap matrix: "
            f"{summary.get('transfer_rows', 0)} transfer rows, "
            f"{summary.get('countermodel_gates', 0)} countermodel gates, "
            f"{summary.get('open_requirement_rows', 0)} open requirements, "
            f"{summary.get('rejected_shortcut_rows', 0)} rejected shortcuts, "
            f"{summary.get('ready_to_apply_rows', 0)} ready-to-apply rows, "
            f"{len(issues)} issues"
        )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
