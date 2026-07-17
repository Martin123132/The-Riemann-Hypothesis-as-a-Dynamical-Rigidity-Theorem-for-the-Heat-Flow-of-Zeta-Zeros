#!/usr/bin/env python3
"""Prove coarse normalized eleventh- and twelfth-cumulant corridors."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_compound_order6_high_cumulant_coarse_corridor as order6


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_high_cumulant_coarse_corridor.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order7_high_cumulant_coarse_corridor.md"
)
TARGET_ORDERS = (11, 12)
FINITE_FORMAL_CAP = 14000
RAY_FORMAL_CAP = 700000
FINITE_EXACT_CAP = 14001
RAY_EXACT_CAP = 700001
EXACT_CORRIDOR_CAP = 1000000


@dataclass(frozen=True)
class CorridorRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def exact_audit() -> dict:
    expressions = order6.formal_scaled_expressions(TARGET_ORDERS)
    finite_caps = {
        3: Fraction(6, 5),
        4: Fraction(3, 2),
        5: Fraction(2),
        6: Fraction(3),
        7: Fraction(9, 2),
        8: Fraction(7),
        9: Fraction(12),
        10: Fraction(21),
        11: Fraction(38),
        12: Fraction(71),
    }
    ray_caps = {order: Fraction(2) for order in range(3, 13)}
    symbols = [sp.symbols(f"L_{order}") for order in range(3, 13)]
    epsilon = sp.symbols("epsilon")
    rows = {}
    for order, expression in expressions.items():
        finite_bound = order6.termwise_bound(
            expression,
            finite_caps,
            order6.FINITE_Q_INVERSE_CAP,
        )
        ray_bound = order6.termwise_bound(
            expression,
            ray_caps,
            order6.RAY_Q_INVERSE_CAP,
        )
        if finite_bound >= FINITE_FORMAL_CAP:
            raise RuntimeError(f"finite formal kappa_{order} cap failed")
        if ray_bound >= RAY_FORMAL_CAP:
            raise RuntimeError(f"ray formal kappa_{order} cap failed")
        rows[str(order)] = {
            "scaled_expression": sp.sstr(expression),
            "term_count": len(sp.Poly(expression, *symbols, epsilon).terms()),
            "finite_termwise_bound": str(finite_bound),
            "finite_cap": FINITE_FORMAL_CAP,
            "ray_termwise_bound": str(ray_bound),
            "ray_cap": RAY_FORMAL_CAP,
        }

    cauchy_factor = 12 * 11
    finite_residual = cauchy_factor * Fraction(19, 123750)
    ray_residual = cauchy_factor * Fraction(1099, 9900000)
    if finite_residual >= 1 or ray_residual >= 1:
        raise RuntimeError("order-seven Cauchy residual cap failed")
    if FINITE_FORMAL_CAP + 1 > FINITE_EXACT_CAP:
        raise RuntimeError("finite exact cumulant reserve failed")
    if RAY_FORMAL_CAP + 1 > RAY_EXACT_CAP:
        raise RuntimeError("ray exact cumulant reserve failed")
    if max(FINITE_FORMAL_CAP, RAY_FORMAL_CAP) + 1 >= EXACT_CORRIDOR_CAP:
        raise RuntimeError("order-seven exact cumulant reserve failed")
    return {
        "formal_rows": rows,
        "formal_orders": list(TARGET_ORDERS),
        "formal_term_counts": sum(row["term_count"] for row in rows.values()),
        "cauchy_factor": cauchy_factor,
        "cauchy_factor_identity": "r!/(r-2)!=r*(r-1)<=132 for r=11,12",
        "finite_scaled_residual_bound": str(finite_residual),
        "ray_scaled_residual_bound": str(ray_residual) + "/u",
        "finite_formal_cap": FINITE_FORMAL_CAP,
        "ray_formal_cap": RAY_FORMAL_CAP,
        "finite_exact_corridor_cap": FINITE_EXACT_CAP,
        "finite_exact_corridor": (
            "|kappa_r|*q^(r/2-1)/(r-2)!<14001, r=11,12, 2<=u<=20"
        ),
        "ray_exact_corridor_cap": RAY_EXACT_CAP,
        "ray_exact_corridor": (
            "|kappa_r|*q^(r/2-1)/(r-2)!<700001, r=11,12, u>=20"
        ),
        "exact_corridor_cap": EXACT_CORRIDOR_CAP,
        "exact_corridor": (
            "|kappa_r|*q^(r/2-1)/(r-2)!<1000000, r=11,12, u>=2"
        ),
    }


def build_artifact() -> dict:
    sources = order6.source_contract()
    exact = exact_audit()
    rows = [
        CorridorRow(
            "co7hccc_01_formal_extension",
            "exact_formal_algebra",
            "ready_to_apply",
            "The epsilon-ten partition recurrence gives exact scaled formal cumulants eleven and twelve.",
            "scaled kappa_r^[10], r=11,12, through epsilon^10",
            "Exact finite symbolic algebra only.",
            {
                "term_counts": {
                    key: value["term_count"]
                    for key, value in exact["formal_rows"].items()
                }
            },
        ),
        CorridorRow(
            "co7hccc_02_finite_formal_bound",
            "exact_inequality",
            "ready_to_apply",
            "Termwise substitution of the proved finite normalized-potential boxes gives a coarse strict formal bound.",
            "|scaled kappa_r^[10]|<14000 on 2<=u<=20, r=11,12",
            "Rational termwise majorization; no sampling.",
        ),
        CorridorRow(
            "co7hccc_03_ray_formal_bound",
            "exact_inequality",
            "ready_to_apply",
            "The asymptotic normalized-potential box and q floor give a coarse strict ray bound.",
            "|scaled kappa_r^[10]|<700000 on u>=20, r=11,12",
            "Rational termwise majorization on the complete ray.",
        ),
        CorridorRow(
            "co7hccc_04_cauchy_extension",
            "exact_analytic_consequence",
            "ready_to_apply",
            "The unit-disk partition residual controls derivatives eleven and twelve with the exact factor 132.",
            exact["cauchy_factor_identity"],
            "Same proved complex-disk residual; exact Cauchy arithmetic only.",
        ),
        CorridorRow(
            "co7hccc_05_exact_corridor",
            "exact_analytic_theorem",
            "ready_to_apply",
            "Formal and exact-minus-formal bounds compose to finite, asymptotic, and global coarse corridors for both new cumulants.",
            exact["exact_corridor"],
            "Exact first-summand tilted cumulants on u>=2.",
            {
                "finite": exact["finite_exact_corridor"],
                "ray": exact["ray_exact_corridor"],
            },
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order7_high_cumulant_coarse_corridor",
        "date": "2026-07-14",
        "status": "global coarse exact eleventh- and twelfth-cumulant corridor theorem",
        "proof_boundary": (
            "This artifact proves only a deliberately coarse absolute corridor "
            "for the normalized eleventh and twelfth tilted cumulants. It does "
            "not prove order-seven curvature, entry, PF-infinity, RH, or Lambda<=0."
        ),
        "source_contract": sources,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": len(rows),
            "formal_orders": 2,
            "formal_terms": exact["formal_term_counts"],
            "cauchy_extensions": 1,
            "global_coarse_corridors": 2,
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order7_high_cumulant_coarse_corridor.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order7_high_cumulant_coarse_corridor.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Jensen-Window PF Compound Order-Seven High-Cumulant Coarse Corridor",
        "",
        "Date: 2026-07-14",
        "",
        "Status: global coarse exact eleventh- and twelfth-cumulant corridor",
        "theorem. This is not a proof of order-seven curvature, PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "## Formal Bounds",
        "",
        "The exact epsilon-ten tilted-Gaussian recurrence needs no potential",
        "orders beyond the already-certified `L_3,...,L_12` box. Rational",
        "termwise substitution gives",
        "",
        "```text",
        "|scaled kappa_r^[10]|<14000 for 2<=u<=20,",
        "|scaled kappa_r^[10]|<700000 for u>=20,",
        f"formal terms={exact['formal_term_counts']}.",
        "```",
        "",
        "The exact Cauchy factor is `12*11=132`; both scaled residual budgets",
        "are below one. Hence the finite and asymptotic source boxes give",
        "",
        "```text",
        exact["finite_exact_corridor"],
        exact["ray_exact_corridor"],
        exact["exact_corridor"],
        "```",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order7_high_cumulant_coarse_corridor.json",
        "outputs/jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    print(
        "wrote order-seven high-cumulant corridor: "
        f"{artifact['summary']['formal_terms']} formal terms, "
        f"{artifact['summary']['global_coarse_corridors']} exact corridors"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
