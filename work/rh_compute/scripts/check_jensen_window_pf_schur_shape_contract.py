#!/usr/bin/env python3
"""Validate the Jensen-window PF Schur shape contract."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CONTRACT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_schur_shape_contract.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_schur_shape_contract.md"

EXPECTED_GRID = {
    2: {
        "tested_minors": 12020,
        "upper_toeplitz_structural_nonzero": 4488,
        "upper_toeplitz_structural_zero": 7532,
        "finite_band_structural_nonzero": 2907,
        "finite_band_zero_after_upper_support": 1581,
        "unique_finite_band_shapes": 1848,
    },
    3: {
        "tested_minors": 12020,
        "upper_toeplitz_structural_nonzero": 4488,
        "upper_toeplitz_structural_zero": 7532,
        "finite_band_structural_nonzero": 3960,
        "finite_band_zero_after_upper_support": 528,
        "unique_finite_band_shapes": 2654,
    },
    4: {
        "tested_minors": 12020,
        "upper_toeplitz_structural_nonzero": 4488,
        "upper_toeplitz_structural_zero": 7532,
        "finite_band_structural_nonzero": 4370,
        "finite_band_zero_after_upper_support": 118,
        "unique_finite_band_shapes": 2989,
    },
    5: {
        "tested_minors": 12020,
        "upper_toeplitz_structural_nonzero": 4488,
        "upper_toeplitz_structural_zero": 7532,
        "finite_band_structural_nonzero": 4472,
        "finite_band_zero_after_upper_support": 16,
        "unique_finite_band_shapes": 3078,
    },
}

EXPECTED_FRONTIERS = {
    "d3_column_shape_m8": {
        "term_count": 10,
        "lambda": [1, 1, 1, 1, 1, 1, 1, 1],
        "mu": [0, 0, 0, 0, 0, 0, 0, 0],
        "determinant_polynomial_h": "-3*h0**5*h2*h3**2 + 6*h0**4*h1**2*h3**2 + 12*h0**4*h1*h2**2*h3 + h0**4*h2**4 - 20*h0**3*h1**3*h2*h3 - 10*h0**3*h1**2*h2**3 + 6*h0**2*h1**5*h3 + 15*h0**2*h1**4*h2**2 - 7*h0*h1**6*h2 + h1**8",
    },
    "d4_column_shape_m6": {
        "term_count": 9,
        "lambda": [1, 1, 1, 1, 1, 1],
        "mu": [0, 0, 0, 0, 0, 0],
        "determinant_polynomial_h": "2*h0**4*h2*h4 + h0**4*h3**2 - 3*h0**3*h1**2*h4 - 6*h0**3*h1*h2*h3 - h0**3*h2**3 + 4*h0**2*h1**3*h3 + 6*h0**2*h1**2*h2**2 - 5*h0*h1**4*h2 + h1**6",
    },
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Schur Shape Contract",
    "Status: exact shape-contract diagnostic",
    "This is not a proof",
    "work/rh_compute/results/jensen_window_pf_schur_shape_contract.json",
    "python work/rh_compute/scripts/jensen_window_pf_schur_shape_contract.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_schur_shape_contract.py",
    "validated Jensen-window PF Schur shape contract: 4 grid rows, 0 issues, 2 frontier rows",
    "h_j = binom(d,j) * A_{n+j}",
    "finite-band structural nonzero",
    "d3_column_shape_m8",
    "d4_column_shape_m6",
    "mixed-sign h-monomial expansions",
    "positive Schur, network, or determinant-integral route",
    "outputs/jensen_window_pf_column_recurrence_contract.md",
    "work/rh_compute/results/jensen_window_pf_column_recurrence_contract.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_column_recurrence_contract.py",
    "validated Jensen-window PF column recurrence contract: 4 degree rows, 0 issues, 2 hard frontier rows",
    "C_m = h_0^m * e_m",
)


@dataclass(frozen=True)
class ShapeIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> ShapeIssue:
    return ShapeIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_top_level(contract: dict) -> list[ShapeIssue]:
    issues: list[ShapeIssue] = []
    if contract.get("kind") != "jensen_window_pf_schur_shape_contract":
        issues.append(issue("<contract>", "bad-kind", repr(contract.get("kind"))))
    if contract.get("target_obligation") != "jwpf_06_sign_regular_to_jensen_pf_conversion":
        issues.append(issue("<contract>", "bad-target-obligation", repr(contract.get("target_obligation"))))
    boundary = str(contract.get("proof_boundary", "")).lower()
    if "not a proof" not in boundary or "lambda <= 0" not in boundary:
        issues.append(issue("<contract>", "weak-proof-boundary", str(contract.get("proof_boundary", ""))))

    summary = contract.get("summary", {})
    expected = {
        "grid_rows": 4,
        "frontier_rows": 2,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("<summary>", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("schur positivity", "finite-support", "mixed-sign"):
        if required not in finding:
            issues.append(issue("<summary>", "missing-main-finding-text", required))
    return issues


def validate_specialization(contract: dict) -> list[ShapeIssue]:
    spec = contract.get("finite_window_specialization", {})
    expected = {
        "h_j": "binom(d,j) * A_{n+j}",
        "support": "h_j = 0 for j < 0 or j > d",
        "toeplitz_minor": "det(h_{q_j-r_i})",
        "jacobi_trudi_form": "s_{lambda/mu}(h_0,...,h_d)",
    }
    issues: list[ShapeIssue] = []
    for key, value in expected.items():
        if spec.get(key) != value:
            issues.append(issue("finite_window_specialization", f"bad-{key}", repr(spec.get(key))))
    return issues


def validate_grid(contract: dict) -> list[ShapeIssue]:
    issues: list[ShapeIssue] = []
    grid = contract.get("bounded_grid", {})
    if grid.get("matrix_size") != 8:
        issues.append(issue("bounded_grid", "bad-matrix-size", repr(grid.get("matrix_size"))))
    if grid.get("max_order") != 5:
        issues.append(issue("bounded_grid", "bad-max-order", repr(grid.get("max_order"))))
    if grid.get("total_tested_minors") != 48080:
        issues.append(issue("bounded_grid", "bad-total-tested-minors", repr(grid.get("total_tested_minors"))))
    if grid.get("total_finite_band_nonzero") != 15709:
        issues.append(issue("bounded_grid", "bad-total-finite-band-nonzero", repr(grid.get("total_finite_band_nonzero"))))

    rows = grid.get("rows", [])
    by_degree = {row.get("degree"): row for row in rows if isinstance(row, dict)}
    for degree, expected in EXPECTED_GRID.items():
        row = by_degree.get(degree)
        if row is None:
            issues.append(issue("bounded_grid", "missing-degree-row", str(degree)))
            continue
        for key, value in expected.items():
            if row.get(key) != value:
                issues.append(issue(f"degree-{degree}", f"bad-{key}", f"{row.get(key)!r} != {value!r}"))
        example = row.get("first_noncontiguous_shape_example", {})
        if example.get("rows") != [0, 1] or example.get("cols") != [0, 2]:
            issues.append(issue(f"degree-{degree}", "bad-first-noncontiguous-example", repr(example)))
    return issues


def validate_frontiers(contract: dict) -> list[ShapeIssue]:
    issues: list[ShapeIssue] = []
    rows = contract.get("hard_frontier_column_shapes", [])
    by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for row_id, expected in EXPECTED_FRONTIERS.items():
        row = by_id.get(row_id)
        if row is None:
            issues.append(issue("hard_frontier_column_shapes", "missing-row", row_id))
            continue
        for key, value in expected.items():
            if row.get(key) != value:
                issues.append(issue(row_id, f"bad-{key}", f"{row.get(key)!r} != {value!r}"))
        if row.get("has_mixed_monomial_signs") is not True:
            issues.append(issue(row_id, "missing-mixed-sign-flag", repr(row.get("has_mixed_monomial_signs"))))
        if row.get("substitution") != "h_j = binom(d,j) * A_{n+j}":
            issues.append(issue(row_id, "bad-substitution", repr(row.get("substitution"))))
    return issues


def validate_note(path: Path) -> list[ShapeIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ShapeIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "schur positivity is proved", "pf-infinity is proved"):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(contract_path: Path, note_path: Path) -> tuple[list[ShapeIssue], int, int]:
    contract = load_json(contract_path)
    issues: list[ShapeIssue] = []
    issues.extend(validate_top_level(contract))
    issues.extend(validate_specialization(contract))
    issues.extend(validate_grid(contract))
    issues.extend(validate_frontiers(contract))
    issues.extend(validate_note(note_path))
    summary = contract.get("summary", {})
    return issues, int(summary.get("grid_rows", 0)), int(summary.get("frontier_rows", 0))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, grid_rows, frontier_rows = validate(args.contract, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "grid_rows": grid_rows,
                    "frontier_rows": frontier_rows,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-SCHUR-SHAPE {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF Schur shape contract: "
            f"{grid_rows} grid rows, {len(issues)} issues, {frontier_rows} frontier rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
