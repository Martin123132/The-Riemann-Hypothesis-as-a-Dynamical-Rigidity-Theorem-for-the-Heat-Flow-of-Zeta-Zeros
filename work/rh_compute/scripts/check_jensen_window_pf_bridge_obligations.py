#!/usr/bin/env python3
"""Validate the Jensen-window PF bridge obligation ledger."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_LEDGER = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_bridge_obligations.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_bridge_obligations.md"


ALLOWED_STATUSES = {
    "available_exact",
    "finite_evidence",
    "guard_validated",
    "open_obligation",
    "conditional_consequence",
    "rejected_by_countermodel",
    "route_mismatch",
}

ALLOWED_ROLES = {
    "exact_equivalence",
    "exact_contact",
    "algebraic_obstruction",
    "finite_evidence",
    "open_antecedent",
    "open_bridge_theorem",
    "open_uniformity_obligation",
    "conditional_consequence",
    "rejection_test",
    "route_separation",
}

ROLE_STATUSES = {
    "exact_equivalence": {"available_exact"},
    "exact_contact": {"guard_validated"},
    "algebraic_obstruction": {"guard_validated"},
    "finite_evidence": {"finite_evidence"},
    "open_antecedent": {"open_obligation"},
    "open_bridge_theorem": {"open_obligation"},
    "open_uniformity_obligation": {"open_obligation"},
    "conditional_consequence": {"conditional_consequence"},
    "rejection_test": {"rejected_by_countermodel"},
    "route_separation": {"route_mismatch"},
}

NON_CLOSING_ROLES = {
    "exact_equivalence",
    "exact_contact",
    "algebraic_obstruction",
    "finite_evidence",
    "conditional_consequence",
    "rejection_test",
    "route_separation",
}

REQUIRED_IDS = {
    "jwpf_01_window_pf_jensen_equivalence",
    "jwpf_02_degree2_signed_hankel_contact",
    "jwpf_03_low_degree_extra_toeplitz_obligations",
    "jwpf_04_current_finite_pf_sturm_evidence",
    "jwpf_05_all_order_shifted_sign_consistency",
    "jwpf_05b_weaker_xi_specific_antecedent",
    "jwpf_06_sign_regular_to_jensen_pf_conversion",
    "jwpf_07_binomial_weight_and_shift_uniformity",
    "jwpf_08_jensen_to_laguerre_polya_limit",
    "jwpf_09_finite_rectangle_promotion_rejected",
    "jwpf_10_ordinary_coefficient_pf_route_separated",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Bridge Obligations",
    "Status: theorem-obligation ledger",
    "This is not a proof of PF-infinity",
    "work/rh_compute/results/jensen_window_pf_bridge_obligations.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_bridge_obligations.py",
    "validated Jensen-window PF bridge obligations: 11 obligations, 0 issues, 3 open obligations",
    "jwpf_06_sign_regular_to_jensen_pf_conversion",
    "jwpf_05b_weaker_xi_specific_antecedent",
    "Q_(10,n)(-100)<0",
    "outputs/jensen_window_pf_endpoint_order10_counterexample.md",
    "would_close_target=true",
    "outputs/jensen_window_pf_bridge_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
    "outputs/jensen_window_pf_theorem_machinery_fit_matrix.md",
    "work/rh_compute/results/proof_claim_ledger.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_theorem_machinery_fit_matrix.py",
    "validated Jensen-window PF theorem machinery fit matrix: 7 rows, 0 issues, 0 ready-to-apply rows",
    "outputs/jensen_window_pf_structural_ansatz_matrix.md",
    "work/rh_compute/results/jensen_window_pf_structural_ansatz_matrix.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py",
    "validated Jensen-window PF structural ansatz matrix: 6 ansatz rows, 0 issues, 0 ready-to-apply rows",
)


@dataclass(frozen=True)
class ObligationIssue:
    obligation_id: str
    issue: str
    detail: str


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(obligation_id: str, ref: str) -> list[ObligationIssue]:
    if ref.startswith(("http://", "https://")):
        return []
    if (REPO_ROOT / ref).exists():
        return []
    return [ObligationIssue(obligation_id, "missing-ref", ref)]


def has_boundary(text: str) -> bool:
    lowered = text.lower()
    return any(
        marker in lowered
        for marker in (
            "not ",
            "not a ",
            "only",
            "open",
            "finite",
            "conditional",
            "rejected",
            "separate",
        )
    )


def validate_obligation(row: dict) -> list[ObligationIssue]:
    issues: list[ObligationIssue] = []
    obligation_id = str(row.get("id", "<missing-id>"))
    for key in ("id", "role", "status", "refs", "claim", "needed_upgrade", "would_close_target", "proof_boundary"):
        if key not in row:
            issues.append(ObligationIssue(obligation_id, "missing-field", key))

    role = row.get("role")
    status = row.get("status")
    if role not in ALLOWED_ROLES:
        issues.append(ObligationIssue(obligation_id, "bad-role", repr(role)))
    if status not in ALLOWED_STATUSES:
        issues.append(ObligationIssue(obligation_id, "bad-status", repr(status)))
    if role in ROLE_STATUSES and status not in ROLE_STATUSES[role]:
        issues.append(ObligationIssue(obligation_id, "role-status-mismatch", f"{role!r} -> {status!r}"))

    would_close = row.get("would_close_target")
    if not isinstance(would_close, bool):
        issues.append(ObligationIssue(obligation_id, "bad-would-close-target", repr(would_close)))
    if role in NON_CLOSING_ROLES and would_close is True:
        issues.append(ObligationIssue(obligation_id, "nonclosing-row-closes-target", role))
    if would_close is True and status != "open_obligation":
        issues.append(ObligationIssue(obligation_id, "closing-row-not-open", repr(status)))

    refs = row.get("refs")
    if not isinstance(refs, list) or not refs:
        issues.append(ObligationIssue(obligation_id, "missing-refs", "refs must be a nonempty list"))
    else:
        for ref in refs:
            if not isinstance(ref, str):
                issues.append(ObligationIssue(obligation_id, "bad-ref", repr(ref)))
            else:
                issues.extend(validate_ref(obligation_id, ref))

    boundary = str(row.get("proof_boundary", ""))
    if not has_boundary(boundary):
        issues.append(ObligationIssue(obligation_id, "weak-proof-boundary", boundary))

    text = f"{row.get('claim', '')} {row.get('needed_upgrade', '')} {boundary}".lower()
    if row.get("status") in {"finite_evidence", "guard_validated", "rejected_by_countermodel", "route_mismatch"}:
        for forbidden in ("therefore rh", "therefore lambda <= 0", "proves lambda <= 0"):
            if forbidden in text:
                issues.append(ObligationIssue(obligation_id, "forbidden-overclaim", forbidden))
    return issues


def validate_note(path: Path) -> list[ObligationIssue]:
    if not path.exists():
        return [ObligationIssue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ObligationIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(ObligationIssue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "the bridge is proved"):
        if forbidden in lowered:
            issues.append(ObligationIssue("note", "forbidden-text", forbidden))
    return issues


def validate(ledger_path: Path, note_path: Path) -> tuple[list[ObligationIssue], int, int]:
    ledger = load_json(ledger_path)
    issues: list[ObligationIssue] = []
    if ledger.get("kind") != "jensen_window_pf_bridge_obligation_ledger":
        issues.append(ObligationIssue("<ledger>", "bad-kind", repr(ledger.get("kind"))))
    if ledger.get("target") != "target_jensen_window_pf_bridge":
        issues.append(ObligationIssue("<ledger>", "bad-target", repr(ledger.get("target"))))
    boundary = str(ledger.get("proof_boundary", "")).lower()
    if "not a proof" not in boundary or "lambda <= 0" not in boundary:
        issues.append(ObligationIssue("<ledger>", "weak-proof-boundary", ledger.get("proof_boundary", "")))

    obligations = ledger.get("obligations", [])
    if not isinstance(obligations, list) or not obligations:
        issues.append(ObligationIssue("<ledger>", "missing-obligations", "obligations must be a nonempty list"))
        obligations = []

    seen: set[str] = set()
    open_count = 0
    closing_count = 0
    for row in obligations:
        if not isinstance(row, dict):
            issues.append(ObligationIssue("<ledger>", "bad-obligation", repr(row)))
            continue
        obligation_id = str(row.get("id", "<missing-id>"))
        if obligation_id in seen:
            issues.append(ObligationIssue(obligation_id, "duplicate-id", obligation_id))
        seen.add(obligation_id)
        if row.get("status") == "open_obligation":
            open_count += 1
        if row.get("would_close_target") is True:
            closing_count += 1
        issues.extend(validate_obligation(row))

    for missing in sorted(REQUIRED_IDS - seen):
        issues.append(ObligationIssue(missing, "missing-required-obligation", missing))
    if open_count < 3:
        issues.append(ObligationIssue("<ledger>", "too-few-open-obligations", str(open_count)))
    if closing_count != 1:
        issues.append(ObligationIssue("<ledger>", "bad-closing-obligation-count", str(closing_count)))
    if "jwpf_04_current_finite_pf_sturm_evidence" in seen:
        finite_row = next(row for row in obligations if isinstance(row, dict) and row.get("id") == "jwpf_04_current_finite_pf_sturm_evidence")
        boundary_text = f"{finite_row.get('claim', '')} {finite_row.get('proof_boundary', '')}".lower()
        if "not all-minor jensen-window pf-infinity" not in boundary_text:
            issues.append(
                ObligationIssue(
                    "jwpf_04_current_finite_pf_sturm_evidence",
                    "missing-finite-boundary",
                    "must say not all-minor Jensen-window PF-infinity",
                )
            )

    issues.extend(validate_note(note_path))
    return issues, len(obligations), open_count


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, obligation_count, open_count = validate(args.ledger, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "obligations": obligation_count,
                    "open_obligations": open_count,
                    "issues": [asdict(issue) for issue in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JENSEN-WINDOW-PF-OBLIGATION {item.obligation_id} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF bridge obligations: "
            f"{obligation_count} obligations, {len(issues)} issues, {open_count} open obligations"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
