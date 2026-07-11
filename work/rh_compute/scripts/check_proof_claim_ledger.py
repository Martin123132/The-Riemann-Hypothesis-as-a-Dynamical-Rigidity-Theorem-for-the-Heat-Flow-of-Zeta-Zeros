#!/usr/bin/env python3
"""Validate the proof-programme claim ledger.

The ledger is a machine-readable separation of exact lemmas, finite
certificates, diagnostics, theorem targets, countermodel gates, forbidden
promotions, and hygiene gates.  This checker enforces that open bridge targets
are not accidentally recorded as proved and that referenced artifacts exist.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_LEDGER = REPO_ROOT / "work/rh_compute/results/proof_claim_ledger.json"


ALLOWED_CATEGORIES = {
    "exact_lemma",
    "finite_certificate",
    "interval_certificate",
    "diagnostic",
    "algebraic_reindexing",
    "countermodel_gate",
    "theorem_target",
    "forbidden_promotion",
    "hygiene_gate",
}


ALLOWED_STATUSES = {
    "available_exact",
    "finite_validated",
    "interval_validated",
    "diagnostic_validated",
    "guard_validated",
    "open_target",
    "rejected_by_countermodel",
}


CATEGORY_STATUSES = {
    "exact_lemma": {"available_exact"},
    "finite_certificate": {"finite_validated"},
    "interval_certificate": {"interval_validated"},
    "diagnostic": {"diagnostic_validated"},
    "algebraic_reindexing": {"guard_validated"},
    "countermodel_gate": {"guard_validated"},
    "theorem_target": {"open_target"},
    "forbidden_promotion": {"rejected_by_countermodel"},
    "hygiene_gate": {"guard_validated"},
}


@dataclass(frozen=True)
class LedgerIssue:
    claim_id: str
    issue: str
    detail: str


def load_ledger(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def has_boundary(text: str) -> bool:
    lowered = text.lower()
    return any(
        marker in lowered
        for marker in (
            "not ",
            "not a ",
            "not currently proved",
            "only",
            "open theorem target",
            "rejected",
            "finite",
            "proof-safety",
        )
    )


def validate_refs(claim: dict) -> list[LedgerIssue]:
    issues: list[LedgerIssue] = []
    claim_id = str(claim.get("id", "<missing-id>"))
    refs = claim.get("refs")
    if not isinstance(refs, list) or not refs:
        return [LedgerIssue(claim_id, "missing-refs", "claim must have at least one reference")]
    for ref in refs:
        if not isinstance(ref, str):
            issues.append(LedgerIssue(claim_id, "bad-ref", f"non-string ref {ref!r}"))
            continue
        if ref.startswith("http://") or ref.startswith("https://"):
            continue
        if not (REPO_ROOT / ref).exists():
            issues.append(LedgerIssue(claim_id, "missing-ref", ref))
    return issues


def validate_claim(claim: dict) -> list[LedgerIssue]:
    issues: list[LedgerIssue] = []
    claim_id = str(claim.get("id", "<missing-id>"))

    for key in ("id", "title", "category", "status", "claim", "proof_boundary"):
        if key not in claim or not claim[key]:
            issues.append(LedgerIssue(claim_id, "missing-field", key))

    category = claim.get("category")
    status = claim.get("status")
    if category not in ALLOWED_CATEGORIES:
        issues.append(LedgerIssue(claim_id, "bad-category", repr(category)))
    if status not in ALLOWED_STATUSES:
        issues.append(LedgerIssue(claim_id, "bad-status", repr(status)))
    if category in CATEGORY_STATUSES and status not in CATEGORY_STATUSES[category]:
        issues.append(
            LedgerIssue(
                claim_id,
                "category-status-mismatch",
                f"{category!r} cannot have status {status!r}",
            )
        )

    boundary = str(claim.get("proof_boundary", ""))
    if not has_boundary(boundary):
        issues.append(LedgerIssue(claim_id, "weak-proof-boundary", boundary))

    if category in {"finite_certificate", "interval_certificate", "diagnostic", "algebraic_reindexing", "countermodel_gate", "hygiene_gate"}:
        if "validation_command" not in claim:
            issues.append(LedgerIssue(claim_id, "missing-validation-command", "validated claim needs executable gate"))

    if category == "theorem_target":
        if status != "open_target":
            issues.append(LedgerIssue(claim_id, "target-not-open", f"theorem target has status {status!r}"))
        for key in ("current_blocker", "required_upgrade"):
            if key not in claim or not claim[key]:
                issues.append(LedgerIssue(claim_id, "missing-target-field", key))
        forbidden_words = ("proved", "validated", "complete")
        title_status = f"{claim.get('title', '')} {claim.get('status', '')}".lower()
        if any(word in title_status for word in forbidden_words if word != "open"):
            issues.append(LedgerIssue(claim_id, "target-overpromoted", title_status))

    if category == "forbidden_promotion" and status != "rejected_by_countermodel":
        issues.append(LedgerIssue(claim_id, "forbidden-promotion-not-rejected", str(status)))

    issues.extend(validate_refs(claim))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    ledger = load_ledger(args.ledger)
    claims = ledger.get("claims", [])
    issues: list[LedgerIssue] = []

    if ledger.get("kind") != "rh_proof_programme_claim_ledger":
        issues.append(LedgerIssue("<ledger>", "bad-kind", repr(ledger.get("kind"))))
    if not isinstance(claims, list) or not claims:
        issues.append(LedgerIssue("<ledger>", "missing-claims", "claims must be a nonempty list"))
        claims = []

    seen_ids: set[str] = set()
    category_counts: dict[str, int] = {}
    status_counts: dict[str, int] = {}
    for claim in claims:
        claim_id = str(claim.get("id", "<missing-id>"))
        if claim_id in seen_ids:
            issues.append(LedgerIssue(claim_id, "duplicate-id", claim_id))
        seen_ids.add(claim_id)
        category_counts[str(claim.get("category"))] = category_counts.get(str(claim.get("category")), 0) + 1
        status_counts[str(claim.get("status"))] = status_counts.get(str(claim.get("status")), 0) + 1
        issues.extend(validate_claim(claim))

    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "claims": len(claims),
                    "category_counts": category_counts,
                    "status_counts": status_counts,
                    "issues": [asdict(issue) for issue in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for issue in issues:
            print(f"CLAIM {issue.claim_id} [{issue.issue}] {issue.detail}")
        print(
            "validated proof-claim ledger: "
            f"{len(claims)} claims, {len(issues)} issues, "
            f"{category_counts.get('theorem_target', 0)} open theorem targets"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
