#!/usr/bin/env python3
"""Validate the Jensen-window PF column recurrence contract."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CONTRACT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_column_recurrence_contract.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_column_recurrence_contract.md"

EXPECTED_DEGREE_TERMS = {
    2: ["+ h1 * C[m-1]", "- h0*h2 * C[m-2]"],
    3: ["+ h1 * C[m-1]", "- h0*h2 * C[m-2]", "+ h0**2*h3 * C[m-3]"],
    4: ["+ h1 * C[m-1]", "- h0*h2 * C[m-2]", "+ h0**2*h3 * C[m-3]", "- h0**3*h4 * C[m-4]"],
    5: [
        "+ h1 * C[m-1]",
        "- h0*h2 * C[m-2]",
        "+ h0**2*h3 * C[m-3]",
        "- h0**3*h4 * C[m-4]",
        "+ h0**4*h5 * C[m-5]",
    ],
}

EXPECTED_HARD_ROWS = {
    "d3_column_recurrence_m8": {
        "source_shape_row": "d3_column_shape_m8",
        "countermodel_value": "-435846079534239/104857600000000",
        "recurrence_polynomial_h": "-3*h0**5*h2*h3**2 + 6*h0**4*h1**2*h3**2 + 12*h0**4*h1*h2**2*h3 + h0**4*h2**4 - 20*h0**3*h1**3*h2*h3 - 10*h0**3*h1**2*h2**3 + 6*h0**2*h1**5*h3 + 15*h0**2*h1**4*h2**2 - 7*h0*h1**6*h2 + h1**8",
    },
    "d4_column_recurrence_m6": {
        "source_shape_row": "d4_column_shape_m6",
        "countermodel_value": "-229760849637/28672000000",
        "recurrence_polynomial_h": "2*h0**4*h2*h4 + h0**4*h3**2 - 3*h0**3*h1**2*h4 - 6*h0**3*h1*h2*h3 - h0**3*h2**3 + 4*h0**2*h1**3*h3 + 6*h0**2*h1**2*h2**2 - 5*h0*h1**4*h2 + h1**6",
    },
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Column Recurrence Contract",
    "Status: column-shape recurrence diagnostic",
    "This is not a proof",
    "work/rh_compute/results/jensen_window_pf_column_recurrence_contract.json",
    "python work/rh_compute/scripts/jensen_window_pf_column_recurrence_contract.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_column_recurrence_contract.py",
    "validated Jensen-window PF column recurrence contract: 4 degree rows, 0 issues, 2 hard frontier rows",
    "C_m = h_0^m * e_m",
    "E(t) = 1 / H(-t)",
    "C[m] = sum_{j=1..min(d,m)}",
    "d3_column_recurrence_m8",
    "d4_column_recurrence_m6",
    "-435846079534239/104857600000000",
    "-229760849637/28672000000",
    "necessary Schur-positivity condition",
    "outputs/jensen_window_pf_column_recurrence_finite_coverage.md",
    "work/rh_compute/results/jensen_window_pf_column_recurrence_finite_coverage.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_column_recurrence_finite_coverage.py",
    "validated Jensen-window PF column recurrence finite coverage: 1470 direct positive rows, 210 hard recurrence rows, 315 Sturm/PF windows, 0 issues",
    "finite evidence only",
)


@dataclass(frozen=True)
class ColumnIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> ColumnIssue:
    return ColumnIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_refs(contract: dict) -> list[ColumnIssue]:
    issues: list[ColumnIssue] = []
    for key in ("source_algebra", "source_schur_shape_contract"):
        ref = contract.get(key)
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            issues.append(issue("<contract>", f"missing-{key}", repr(ref)))
    return issues


def validate_top_level(contract: dict) -> list[ColumnIssue]:
    issues: list[ColumnIssue] = []
    if contract.get("kind") != "jensen_window_pf_column_recurrence_contract":
        issues.append(issue("<contract>", "bad-kind", repr(contract.get("kind"))))
    if contract.get("target_obligation") != "jwpf_06_sign_regular_to_jensen_pf_conversion":
        issues.append(issue("<contract>", "bad-target-obligation", repr(contract.get("target_obligation"))))
    boundary = str(contract.get("proof_boundary", "")).lower()
    if "not a proof" not in boundary or "lambda <= 0" not in boundary:
        issues.append(issue("<contract>", "weak-proof-boundary", str(contract.get("proof_boundary", ""))))

    summary = contract.get("summary", {})
    expected = {
        "degree_rows": 4,
        "hard_frontier_rows": 2,
        "negative_countermodel_rows": 2,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("<summary>", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("elementary-symmetric", "generic log-concavity", "negative"):
        if required not in finding:
            issues.append(issue("<summary>", "missing-main-finding-text", required))
    return issues


def validate_normalization(contract: dict) -> list[ColumnIssue]:
    normalization = contract.get("normalization", {})
    expected = {
        "h_j": "binom(d,j) * A_{n+j}",
        "g_j": "h_j / h_0",
        "C_m": "det(h_{1+col-row}) for rows=0..m-1 and cols=1..m",
        "C_m_relation_to_elementary": "C_m = h_0^m * e_m(g_1,...,g_d)",
        "generating_identity": "sum_m e_m t^m = 1 / (1 - g_1*t + g_2*t^2 - ... + (-1)^d*g_d*t^d)",
    }
    issues: list[ColumnIssue] = []
    for key, value in expected.items():
        if normalization.get(key) != value:
            issues.append(issue("normalization", f"bad-{key}", repr(normalization.get(key))))
    return issues


def validate_degree_rows(contract: dict) -> list[ColumnIssue]:
    issues: list[ColumnIssue] = []
    rows = contract.get("degree_contract_rows", [])
    by_degree = {row.get("degree"): row for row in rows if isinstance(row, dict)}
    for degree, expected_terms in EXPECTED_DEGREE_TERMS.items():
        row = by_degree.get(degree)
        if row is None:
            issues.append(issue("degree_contract_rows", "missing-degree", str(degree)))
            continue
        terms = [term.get("term") for term in row.get("recurrence_terms", []) if isinstance(term, dict)]
        if terms != expected_terms:
            issues.append(issue(f"degree-{degree}", "bad-recurrence-terms", repr(terms)))
        if "C[0]=1" not in str(row.get("unnormalized_column_recurrence", "")):
            issues.append(issue(f"degree-{degree}", "missing-initial-condition", str(row.get("unnormalized_column_recurrence", ""))))
        if row.get("necessary_condition") != "C[m] >= 0 for every m >= 0 and every shift n":
            issues.append(issue(f"degree-{degree}", "bad-necessary-condition", repr(row.get("necessary_condition"))))
    return issues


def validate_hard_rows(contract: dict) -> list[ColumnIssue]:
    issues: list[ColumnIssue] = []
    rows = contract.get("hard_frontier_recurrence_rows", [])
    by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for row_id, expected in EXPECTED_HARD_ROWS.items():
        row = by_id.get(row_id)
        if row is None:
            issues.append(issue("hard_frontier_recurrence_rows", "missing-row", row_id))
            continue
        for key, value in expected.items():
            if row.get(key) != value:
                issues.append(issue(row_id, f"bad-{key}", f"{row.get(key)!r} != {value!r}"))
        if row.get("direct_determinant_polynomial_h") != row.get("recurrence_polynomial_h"):
            issues.append(issue(row_id, "recurrence-determinant-mismatch", repr(row)))
        if row.get("recurrence_matches_determinant") is not True:
            issues.append(issue(row_id, "missing-recurrence-match-flag", repr(row.get("recurrence_matches_determinant"))))
        if row.get("countermodel_negative") is not True:
            issues.append(issue(row_id, "countermodel-not-negative", repr(row.get("countermodel_negative"))))
        if row.get("matches_obligation_algebra_countermodel") is not True:
            issues.append(issue(row_id, "countermodel-not-matched", repr(row.get("matches_obligation_algebra_countermodel"))))
        if "not from generic log-concavity" not in str(row.get("structural_obligation", "")).lower():
            issues.append(issue(row_id, "missing-structural-obligation-warning", str(row.get("structural_obligation", ""))))
    return issues


def validate_note(path: Path) -> list[ColumnIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ColumnIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "schur positivity is proved", "pf-infinity is proved"):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(contract_path: Path, note_path: Path) -> tuple[list[ColumnIssue], int, int]:
    contract = load_json(contract_path)
    issues: list[ColumnIssue] = []
    issues.extend(validate_top_level(contract))
    issues.extend(validate_refs(contract))
    issues.extend(validate_normalization(contract))
    issues.extend(validate_degree_rows(contract))
    issues.extend(validate_hard_rows(contract))
    issues.extend(validate_note(note_path))
    summary = contract.get("summary", {})
    return issues, int(summary.get("degree_rows", 0)), int(summary.get("hard_frontier_rows", 0))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, degree_rows, hard_rows = validate(args.contract, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "degree_rows": degree_rows,
                    "hard_frontier_rows": hard_rows,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-COLUMN-RECURRENCE {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF column recurrence contract: "
            f"{degree_rows} degree rows, {len(issues)} issues, {hard_rows} hard frontier rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
