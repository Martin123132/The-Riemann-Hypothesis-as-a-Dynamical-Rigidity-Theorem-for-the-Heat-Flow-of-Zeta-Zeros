#!/usr/bin/env python3
"""Prove coarse normalized twenty-first/twenty-second-cumulant corridors."""

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
    "jensen_window_pf_compound_order12_high_cumulant_coarse_corridor.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order12_high_cumulant_coarse_corridor.md"
)
SOURCE_EXACT_PARTITION = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_central_residual_certificate.json"
)
TARGET_ORDERS = (21, 22)
EXACT_CORRIDOR_CAP = 1


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
    if any(sp.expand(expression) != 0 for expression in expressions.values()):
        raise RuntimeError("epsilon-ten order-twelve high-cumulant vanishing failed")
    cauchy_factor = max(order * (order - 1) for order in TARGET_ORDERS)
    finite_residual = cauchy_factor * Fraction(19, 123750)
    ray_residual = cauchy_factor * Fraction(1099, 9900000)
    if finite_residual >= EXACT_CORRIDOR_CAP or ray_residual >= EXACT_CORRIDOR_CAP:
        raise RuntimeError("twenty-first/twenty-second Cauchy budget failed")
    return {
        "formal_orders": list(TARGET_ORDERS),
        "formal_rows": {
            str(order): {
                "scaled_expression": sp.sstr(expressions[order]),
                "term_count": 0,
            }
            for order in TARGET_ORDERS
        },
        "formal_term_counts": 0,
        "cauchy_factor": cauchy_factor,
        "cauchy_factor_identity": (
            "r!/(r-2)!=r*(r-1)<=22*21=462 for r=21,22"
        ),
        "finite_scaled_residual_bound": str(finite_residual),
        "ray_scaled_residual_bound": f"{ray_residual}/u",
        "finite_exact_corridor": (
            "|kappa_r|*q^(r/2-1)/(r-2)!<1, r=21,22, 2<=u<=20"
        ),
        "ray_exact_corridor": (
            "|kappa_r|*q^(r/2-1)/(r-2)!<1, r=21,22, u>=20"
        ),
        "exact_corridor": (
            "|kappa_r|*q^(r/2-1)/(r-2)!<1, r=21,22, u>=2"
        ),
    }


def build_artifact() -> dict:
    exact_partition = json.loads(SOURCE_EXACT_PARTITION.read_text(encoding="utf-8"))
    if (
        exact_partition.get("status")
        != "global exact central residual and unit-disk partition theorem"
        or exact_partition.get("summary", {}).get("global_partition_residual_closed")
        is not True
    ):
        raise RuntimeError("exact unit-disk partition source changed")
    sources = {
        **order6.source_contract(),
        "exact_partition_residual": (
            "global unit-disk exact-minus-epsilon-ten partition theorem"
        ),
        "exact_partition_source": SOURCE_EXACT_PARTITION.relative_to(
            REPO_ROOT
        ).as_posix(),
    }
    exact = exact_audit()
    rows = [
        CorridorRow(
            "co12hccc_01_formal_vanishing",
            "exact_formal_algebra",
            "ready_to_apply",
            "The epsilon-ten formal log partition has zero normalized cumulants in orders twenty-one and twenty-two.",
            "scaled kappa_21^[10]=scaled kappa_22^[10]=0",
            "Exact finite symbolic algebra only.",
        ),
        CorridorRow(
            "co12hccc_02_cauchy_extension",
            "exact_analytic_consequence",
            "ready_to_apply",
            "The proved unit-disk log-partition residual controls both new derivatives with Cauchy factor 462.",
            exact["cauchy_factor_identity"],
            "Same exact-minus-epsilon-ten complex-disk theorem; no extrapolation.",
            {
                "finite_scaled_residual_bound": exact[
                    "finite_scaled_residual_bound"
                ],
                "ray_scaled_residual_bound": exact["ray_scaled_residual_bound"],
            },
        ),
        CorridorRow(
            "co12hccc_03_finite_corridor",
            "exact_analytic_theorem",
            "ready_to_apply",
            "Formal vanishing and the finite Cauchy residual give a unit corridor for both new cumulants.",
            exact["finite_exact_corridor"],
            "Exact first-summand tilted cumulants on 2<=u<=20.",
        ),
        CorridorRow(
            "co12hccc_04_global_corridor",
            "exact_analytic_theorem",
            "ready_to_apply",
            "The asymptotic residual is smaller still, so the same unit corridor holds on the complete saddle ray.",
            exact["exact_corridor"],
            "Exact first-summand tilted cumulants on u>=2.",
            {"ray": exact["ray_exact_corridor"]},
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order12_high_cumulant_coarse_corridor",
        "date": "2026-07-22",
        "status": (
            "global coarse exact twenty-first- and twenty-second-cumulant corridor theorem"
        ),
        "proof_boundary": (
            "This artifact proves only an absolute normalized corridor for the "
            "twenty-first and twenty-second tilted cumulants. It does not prove "
            "order-twelve curvature, entry, PF-infinity, RH, or Lambda<=0."
        ),
        "source_contract": sources,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": len(rows),
            "formal_orders": 2,
            "formal_terms": 0,
            "cauchy_extensions": 1,
            "global_coarse_corridors": 2,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order12_high_cumulant_coarse_corridor.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order12_high_cumulant_coarse_corridor.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Jensen-Window PF Compound Order-Twelve High-Cumulant Coarse Corridor",
        "",
        "Date: 2026-07-22",
        "",
        "Status: global coarse exact twenty-first- and twenty-second-cumulant",
        "corridor theorem. This is not a proof of order-twelve curvature,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "## Vanishing And Cauchy Budget",
        "",
        "```text",
        "scaled kappa_21^[10]=scaled kappa_22^[10]=0",
        exact["cauchy_factor_identity"],
        f"finite residual < {exact['finite_scaled_residual_bound']}",
        f"ray residual < {exact['ray_scaled_residual_bound']}",
        exact["finite_exact_corridor"],
        exact["ray_exact_corridor"],
        exact["exact_corridor"],
        "```",
        "",
        "The bound uses the same proved unit-disk exact-minus-epsilon-ten",
        "partition residual as the order-eleven extension. No numerical",
        "extrapolation or RH assumption enters this corridor.",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print(
        "wrote order-twelve high-cumulant corridor: "
        f"{artifact['summary']['formal_terms']} formal terms, "
        f"{artifact['summary']['global_coarse_corridors']} exact corridors"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
