#!/usr/bin/env python3
"""Prove coarse normalized ninth- and tenth-cumulant corridors for order six."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
import math
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.md"
)
SOURCE_FINITE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.json"
)
SOURCE_RAY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate.json"
)
SOURCE_DISK = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.json"
)
MAX_EPSILON_ORDER = 10
MAX_POTENTIAL_ORDER = 12
TARGET_ORDERS = (9, 10)
FINITE_Q_INVERSE_CAP = Fraction(1, 9000)
RAY_Q_INVERSE_CAP = Fraction(1, 10**35)
FINITE_FORMAL_CAP = 1600
RAY_FORMAL_CAP = 36000
EXACT_CORRIDOR_CAP = 50000


@dataclass(frozen=True)
class CorridorRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_contract() -> dict:
    finite = load_json(SOURCE_FINITE)
    ray = load_json(SOURCE_RAY)
    disk = load_json(SOURCE_DISK)
    caps = finite.get("normalized_jet_caps", {})
    expected_caps = {
        "3": "6/5",
        "4": "3/2",
        "5": "2",
        "6": "3",
        "7": "9/2",
        "8": "7",
        "9": "12",
        "10": "21",
        "11": "38",
        "12": "71",
    }
    if caps != expected_caps:
        raise RuntimeError("finite normalized-potential cap source changed")
    if finite.get("parameters", {}).get("mode_interval") != ["2", "20"]:
        raise RuntimeError("finite normalized-potential range changed")
    geometry = ray.get("scalar_geometry", {})
    if geometry.get("full_normalized_domain") != "|L_r|<2, r=3,...,12":
        raise RuntimeError("asymptotic normalized-potential source changed")
    if geometry.get("q_floor") != "100000000000000000000000000000000000":
        raise RuntimeError("asymptotic q floor changed")
    cauchy = disk.get("cauchy_budgets", {})
    if cauchy.get("finite", {}).get("total_log_constant") != "19/123750":
        raise RuntimeError("finite complex-disk constant changed")
    if cauchy.get("ray", {}).get("total_log_constant") != "1099/9900000":
        raise RuntimeError("ray complex-disk constant changed")
    return {
        "finite_caps": expected_caps,
        "finite_mode": ["2", "20"],
        "ray_caps": "|L_r|<2 for r=3,...,12",
        "finite_q_inverse_cap": str(FINITE_Q_INVERSE_CAP),
        "ray_q_inverse_cap": str(RAY_Q_INVERSE_CAP),
        "cauchy_formula": cauchy["cauchy_formula"],
        "finite_log_constant": cauchy["finite"]["total_log_constant"],
        "ray_log_constant": cauchy["ray"]["total_log_constant"],
    }


def formal_scaled_expressions(
    target_orders: tuple[int, ...] = TARGET_ORDERS,
) -> dict[int, sp.Expr]:
    z, y, epsilon = sp.symbols("z y epsilon")
    symbols = {
        order: sp.symbols(f"L_{order}")
        for order in range(3, MAX_POTENTIAL_ORDER + 1)
    }
    perturbation = [sp.Integer(0) for _ in range(MAX_EPSILON_ORDER + 1)]
    for order in range(3, MAX_POTENTIAL_ORDER + 1):
        perturbation[order - 2] = (
            symbols[order] * y**order / sp.factorial(order)
        )

    exponential = [sp.Integer(0) for _ in range(MAX_EPSILON_ORDER + 1)]
    exponential[0] = sp.Integer(1)
    for degree in range(1, MAX_EPSILON_ORDER + 1):
        exponential[degree] = sp.expand(
            -sum(
                index * perturbation[index] * exponential[degree - index]
                for index in range(1, degree + 1)
            )
            / degree
        )
    maximum_y_degree = max(sp.Poly(value, y).degree() for value in exponential)
    gaussian_moments = [sp.Integer(1)]
    for _ in range(maximum_y_degree):
        gaussian_moments.append(
            sp.expand(z * gaussian_moments[-1] + sp.diff(gaussian_moments[-1], z))
        )

    def tilted_expectation(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sum(
                coefficient * gaussian_moments[monomial[0]]
                for monomial, coefficient in sp.Poly(polynomial, y).terms()
            )
        )

    partition = [tilted_expectation(value) for value in exponential]
    logarithm = [sp.Integer(0) for _ in range(MAX_EPSILON_ORDER + 1)]
    for degree in range(1, MAX_EPSILON_ORDER + 1):
        logarithm[degree] = sp.expand(
            partition[degree]
            - sum(
                index * logarithm[index] * partition[degree - index]
                for index in range(1, degree)
            )
            / degree
        )

    expressions = {}
    for order in target_orders:
        scaled = sum(
            sp.diff(logarithm[degree], z, order).subs(z, 0)
            * epsilon ** (degree - (order - 2))
            for degree in range(order - 2, MAX_EPSILON_ORDER + 1, 2)
        )
        expressions[order] = sp.factor(
            sp.Rational((-1) ** order, math.factorial(order - 2)) * scaled
        )
    return expressions


def termwise_bound(
    expression: sp.Expr,
    potential_caps: dict[int, Fraction],
    q_inverse_cap: Fraction,
) -> Fraction:
    symbols = [sp.symbols(f"L_{order}") for order in range(3, 13)]
    epsilon = sp.symbols("epsilon")
    total = Fraction(0)
    for powers, coefficient in sp.Poly(expression, *symbols, epsilon).terms():
        rational = sp.Rational(coefficient)
        term = Fraction(abs(int(rational.p)), int(rational.q))
        for order, power in zip(range(3, 13), powers[:-1]):
            term *= potential_caps[order] ** power
        epsilon_power = powers[-1]
        if epsilon_power % 2:
            raise RuntimeError("scaled cumulant expression has odd epsilon power")
        term *= q_inverse_cap ** (epsilon_power // 2)
        total += term
    return total


def exact_audit() -> dict:
    expressions = formal_scaled_expressions()
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
    rows = {}
    for order, expression in expressions.items():
        finite_bound = termwise_bound(
            expression, finite_caps, FINITE_Q_INVERSE_CAP
        )
        ray_bound = termwise_bound(expression, ray_caps, RAY_Q_INVERSE_CAP)
        if finite_bound >= FINITE_FORMAL_CAP:
            raise RuntimeError(f"finite formal kappa_{order} cap failed")
        if ray_bound >= RAY_FORMAL_CAP:
            raise RuntimeError(f"ray formal kappa_{order} cap failed")
        rows[str(order)] = {
            "scaled_expression": sp.sstr(expression),
            "term_count": len(
                sp.Poly(
                    expression,
                    *[sp.symbols(f"L_{r}") for r in range(3, 13)],
                    sp.symbols("epsilon"),
                ).terms()
            ),
            "finite_termwise_bound": str(finite_bound),
            "finite_cap": FINITE_FORMAL_CAP,
            "ray_termwise_bound": str(ray_bound),
            "ray_cap": RAY_FORMAL_CAP,
        }

    cauchy_factor = 10 * 9
    finite_residual = cauchy_factor * Fraction(19, 123750)
    ray_residual = cauchy_factor * Fraction(1099, 9900000)
    if finite_residual >= 1 or ray_residual >= 1:
        raise RuntimeError("extended Cauchy residual cap failed")
    if max(FINITE_FORMAL_CAP, RAY_FORMAL_CAP) + 1 >= EXACT_CORRIDOR_CAP:
        raise RuntimeError("coarse exact cumulant reserve failed")
    return {
        "formal_rows": rows,
        "formal_orders": list(TARGET_ORDERS),
        "formal_term_counts": sum(row["term_count"] for row in rows.values()),
        "cauchy_factor": cauchy_factor,
        "cauchy_factor_identity": "r!/(r-2)!=r*(r-1)<=90 for r=9,10",
        "finite_scaled_residual_bound": str(finite_residual),
        "ray_scaled_residual_bound": str(ray_residual) + "/u",
        "finite_formal_cap": FINITE_FORMAL_CAP,
        "ray_formal_cap": RAY_FORMAL_CAP,
        "exact_corridor_cap": EXACT_CORRIDOR_CAP,
        "exact_corridor": (
            "|kappa_r|*q^(r/2-1)/(r-2)!<50000, r=9,10, u>=2"
        ),
    }


def build_artifact() -> dict:
    sources = source_contract()
    exact = exact_audit()
    rows = [
        CorridorRow(
            "co6hccc_01_formal_extension",
            "exact_formal_algebra",
            "ready_to_apply",
            "The epsilon-ten partition recurrence gives exact scaled formal cumulants nine and ten.",
            "scaled kappa_r^[10], r=9,10, through epsilon^10",
            "Exact finite symbolic algebra only.",
            {"term_counts": {k: v["term_count"] for k, v in exact["formal_rows"].items()}},
        ),
        CorridorRow(
            "co6hccc_02_finite_formal_bound",
            "exact_inequality",
            "ready_to_apply",
            "Termwise substitution of the proved finite normalized-potential boxes gives a coarse strict formal bound.",
            "|scaled kappa_r^[10]|<1600 on 2<=u<=20, r=9,10",
            "Rational termwise majorization; no sampling.",
        ),
        CorridorRow(
            "co6hccc_03_ray_formal_bound",
            "exact_inequality",
            "ready_to_apply",
            "The asymptotic normalized-potential box and q floor give a coarse strict ray bound.",
            "|scaled kappa_r^[10]|<36000 on u>=20, r=9,10",
            "Rational termwise majorization on the complete ray.",
        ),
        CorridorRow(
            "co6hccc_04_cauchy_extension",
            "exact_analytic_consequence",
            "ready_to_apply",
            "The existing unit-disk partition residual controls derivatives nine and ten after replacing 56 by the exact factor 90.",
            exact["cauchy_factor_identity"],
            "Same proved complex-disk residual; exact Cauchy arithmetic only.",
        ),
        CorridorRow(
            "co6hccc_05_exact_corridor",
            "exact_analytic_theorem",
            "ready_to_apply",
            "Formal and exact-minus-formal bounds compose to a global coarse corridor for both new cumulants.",
            exact["exact_corridor"],
            "Exact first-summand tilted cumulants on u>=2.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order6_high_cumulant_coarse_corridor",
        "date": "2026-07-13",
        "status": "global coarse exact ninth- and tenth-cumulant corridor theorem",
        "proof_boundary": (
            "This artifact proves only a deliberately coarse absolute corridor for "
            "the normalized ninth and tenth tilted cumulants. It does not prove the "
            "order-six curvature theorem by itself, order-six entry, PF-infinity, "
            "RH, or Lambda<=0."
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
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Jensen-Window PF Compound Order-Six High-Cumulant Coarse Corridor",
        "",
        "Date: 2026-07-13",
        "",
        "Status: global coarse exact ninth- and tenth-cumulant corridor theorem.",
        "This is not a proof of order-six curvature, PF-infinity, RH, or",
        "`Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.py",
        "```",
        "",
        "## Formal Bounds",
        "",
        "The exact tilted-Gaussian partition recurrence is extended only far",
        "enough to extract the scaled epsilon-ten cumulants `r=9,10`. Rational",
        "termwise substitution gives",
        "",
        "```text",
        "|scaled kappa_r^[10]|<1600 for 2<=u<=20,",
        "|scaled kappa_r^[10]|<36000 for u>=20.",
        "```",
        "",
        "These constants are intentionally loose; the order-six stable hierarchy",
        "suppresses the new derivative uncertainty by several powers of `q`.",
        "",
        "## Exact-Density Transfer",
        "",
        "The existing complex-disk theorem already controls one analytic",
        "partition residual on `|z|<=1`. Cauchy's estimate is not limited to",
        "order eight. For the two new derivatives,",
        "",
        "```text",
        exact["cauchy_factor_identity"],
        "finite scaled residual < " + exact["finite_scaled_residual_bound"],
        "ray scaled residual < " + exact["ray_scaled_residual_bound"],
        "```",
        "",
        "so the complete exact corridor is",
        "",
        "```text",
        exact["exact_corridor"] + ".",
        "```",
        "",
        "No sign is asserted or needed for these two coarse corridors.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.md",
        "outputs/jensen_window_pf_compound_order6_first_summand_curvature_bridge.md",
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
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-six high-cumulant corridor: "
        f"{summary['formal_orders']} formal orders, "
        f"{summary['formal_terms']} terms, "
        f"{summary['global_coarse_corridors']} exact corridors"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
