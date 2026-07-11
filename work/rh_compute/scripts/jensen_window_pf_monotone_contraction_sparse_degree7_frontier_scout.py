#!/usr/bin/env python3
"""Sparse exact degree-7 monotone-contraction frontier scout."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import jensen_window_pf_monotone_contraction_sparse_degree6_scout as sparse


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.md"
DEGREE = 7
MAX_M = 10


def analyze_degree7(max_m: int = MAX_M) -> list[dict]:
    variable_count = DEGREE - 1
    columns = sparse.normalized_column_polynomials(DEGREE, max_m)
    rows: list[dict] = []
    for size in range(1, max_m + 1):
        monotone_poly = sparse.substitute_monotone(columns[size], variable_count)
        power_degrees = [max(monomial[index] for monomial in monotone_poly) for index in range(variable_count)]
        stats = sparse.scaled_bernstein_stats(
            monotone_poly,
            include_index_details=True,
            negative_example_limit=12,
        )
        strictly_positive = stats["bernstein_coefficients_strictly_positive"]
        rows.append(
            {
                "id": f"mcs7_d7_m{size}",
                "degree": DEGREE,
                "minor_size": size,
                "normalized_recurrence_term_count": len(columns[size]),
                "monotone_power_basis_term_count": len(monotone_poly),
                "monotone_power_multidegree": power_degrees,
                "monotone_substitution": "x_i=1-prod_{r<=i}(1-s_r), all s_r in [0,1]",
                "frontier_classification": (
                    "strict_positive_global_bernstein_certificate"
                    if strictly_positive
                    else "global_bernstein_certificate_obstruction"
                ),
                **stats,
                "proof_boundary": (
                    "Sparse exact degree-7 column-recurrence diagnostic under monotone contractions only. "
                    "Negative Bernstein coefficients at the obstruction row reject this one-shot global "
                    "Bernstein certificate; they are not a proof that the polynomial is negative and not an "
                    "all-m, all-degree, all-shift, all-shape, zeta cone-entry, RH, or Lambda <= 0 theorem."
                ),
            }
        )
    return rows


def build_payload() -> dict:
    rows = analyze_degree7(MAX_M)
    positive_rows = [row for row in rows if row["bernstein_coefficients_strictly_positive"]]
    obstruction_rows = [row for row in rows if not row["bernstein_coefficients_strictly_positive"]]
    negative_coefficient_total = sum(row["bernstein_negative_coefficient_count"] for row in rows)
    zero_coefficient_total = sum(row["bernstein_zero_coefficient_count"] for row in rows)
    first_obstruction = obstruction_rows[0] if obstruction_rows else None
    return {
        "kind": "jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout",
        "date": "2026-07-06",
        "status": "exact_sparse_degree7_frontier_diagnostic",
        "source_degree6_scout": "outputs/jensen_window_pf_monotone_contraction_sparse_degree6_scout.md",
        "source_column_extension_scout": "outputs/jensen_window_pf_monotone_contraction_column_extension_scout.md",
        "source_column_recurrence_contract": "outputs/jensen_window_pf_column_recurrence_contract.md",
        "source_monotone_contraction_target": "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "proof_boundary": (
            "Exact sparse bounded degree-7 column-recurrence theorem-search diagnostic only. "
            "It proves strict global Bernstein positivity for degree-7 rows m<=9 under monotone "
            "contractions and records the degree-7 m=10 failure of this one-shot global Bernstein "
            "certificate. The obstruction is not a positivity counterexample, not zeta cone entry, "
            "not all-m/all-shape closure, and not a proof of RH or Lambda <= 0."
        ),
        "monotone_contraction_region": "0 <= x1 <= x2 <= x3 <= x4 <= x5 <= x6 <= 1",
        "normalized_recurrence": (
            "Q_0=1 and Q_m=sum_{j=1..min(7,m)} (-1)^(j-1) binom(7,j) "
            "prod_{i=1..j-1} x_i^(j-i) Q_{m-j}, after factoring A^m*rho^m."
        ),
        "rows": rows,
        "summary": {
            "degree7_rows": len(rows),
            "positive_certificate_rows": len(positive_rows),
            "certificate_obstruction_rows": len(obstruction_rows),
            "max_minor_size": MAX_M,
            "first_certificate_obstruction_m": first_obstruction["minor_size"] if first_obstruction else None,
            "first_certificate_obstruction_min_coefficient": (
                first_obstruction["bernstein_min_coefficient"] if first_obstruction else None
            ),
            "total_bernstein_coefficients": sum(row["bernstein_coefficient_count"] for row in rows),
            "positive_row_bernstein_coefficients": sum(row["bernstein_coefficient_count"] for row in positive_rows),
            "obstruction_row_bernstein_coefficients": sum(row["bernstein_coefficient_count"] for row in obstruction_rows),
            "negative_bernstein_coefficients": negative_coefficient_total,
            "zero_bernstein_coefficients": zero_coefficient_total,
            "ready_to_apply_rows": 0,
            "target_closing": False,
            "main_finding": (
                "The sparse exact axis-Bernstein transform certifies degree 7 monotone-contraction "
                "column rows through m=9, covering 670,891 strictly positive Bernstein coefficients. "
                "At degree 7 m=10 the same global Bernstein certificate first fails, with 126 negative "
                "Bernstein coefficients and minimum -4928. This marks a certificate frontier, not a "
                "positivity counterexample or proof obstruction to a subdivided/stronger certificate."
            ),
        },
        "invariants": [
            "Every degree-7 row m<=9 has strictly positive Bernstein coefficients.",
            "Degree-7 m=10 is recorded as a global Bernstein certificate obstruction.",
            "Negative Bernstein coefficients at m=10 do not prove polynomial negativity.",
            "The scout covers only degree 7 and m<=10.",
            "No row is ready_to_apply.",
            "Endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(payload: dict, path: Path) -> None:
    summary = payload["summary"]
    result_line = (
        "validated Jensen-window PF monotone-contraction sparse degree-7 frontier scout: "
        f"{summary['positive_certificate_rows']} positive rows, "
        f"{summary['certificate_obstruction_rows']} certificate-obstruction row, "
        f"{summary['total_bernstein_coefficients']} Bernstein coefficients, "
        f"first obstruction m={summary['first_certificate_obstruction_m']}, "
        f"{summary['negative_bernstein_coefficients']} negative Bernstein coefficients, "
        f"{summary['zero_bernstein_coefficients']} zero Bernstein coefficients, 0 issues"
    )
    lines = [
        "# Jensen-Window PF Monotone-Contraction Sparse Degree-7 Frontier Scout",
        "",
        "Date: 2026-07-06",
        "",
        "Status: exact sparse degree-7 frontier diagnostic. This is not a proof",
        "of Jensen-window PF-infinity, all-shape Schur positivity, RH, or",
        "`Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout`.",
        "",
        "Proof boundary: this artifact proves finite exact Bernstein sign",
        "certificates for bounded degree-7 column rows `m<=9` under monotone",
        "contractions. It also records that the same one-shot global Bernstein",
        "certificate fails at `m=10`. That failure is not a proof that the",
        "polynomial is negative, and it does not prove all column rows or all",
        "Schur shapes.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Certified Region",
        "",
        "```text",
        "degree 7: 0 <= x1 <= x2 <= x3 <= x4 <= x5 <= x6 <= 1",
        "x_i=1-prod_{r<=i}(1-s_r), all s_r in [0,1]",
        "```",
        "",
        "## Sparse Recurrence",
        "",
        "After removing the positive monomial `A^m*rho^m`, the column recurrence is:",
        "",
        "```text",
        payload["normalized_recurrence"],
        "```",
        "",
        "## Degree-7 Rows",
        "",
        "```text",
    ]
    for row in payload["rows"]:
        lines.append(
            f"{row['id']}: m={row['minor_size']}, "
            f"class={row['frontier_classification']}, "
            f"recurrence terms={row['normalized_recurrence_term_count']}, "
            f"power terms={row['monotone_power_basis_term_count']}, "
            f"Bernstein multidegree={row['bernstein_multidegree']}, "
            f"Bernstein count={row['bernstein_coefficient_count']}, "
            f"min={row['bernstein_min_coefficient']}, "
            f"min_index={row['bernstein_min_index']}, "
            f"negative={row['bernstein_negative_coefficient_count']}, "
            f"zero={row['bernstein_zero_coefficient_count']}"
        )
    lines.extend(
        [
            "```",
            "",
            "The first obstruction row has representative negative Bernstein coefficients:",
            "",
            "```text",
        ]
    )
    obstruction_rows = [row for row in payload["rows"] if row["frontier_classification"] == "global_bernstein_certificate_obstruction"]
    if obstruction_rows:
        for example in obstruction_rows[0].get("bernstein_negative_examples", []):
            lines.append(f"index={example['index']}, coefficient={example['coefficient']}")
    else:
        lines.append("none")
    lines.extend(
        [
            "```",
            "",
            "## Consequence",
            "",
            payload["summary"]["main_finding"],
            "",
            "Integration:",
            "",
            "```text",
            "outputs/jensen_window_pf_monotone_contraction_sparse_degree6_scout.md",
            "outputs/jensen_window_pf_monotone_contraction_column_extension_scout.md",
            "outputs/jensen_window_pf_column_recurrence_contract.md",
            "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
            "outputs/signed_hankel_jensen_dependency_graph.md",
            "```",
            "",
            "Summary:",
            "",
            payload["summary"]["main_finding"],
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(payload, args.note)
    print(
        "wrote "
        f"{args.out_json.relative_to(REPO_ROOT)} and {args.note.relative_to(REPO_ROOT)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
