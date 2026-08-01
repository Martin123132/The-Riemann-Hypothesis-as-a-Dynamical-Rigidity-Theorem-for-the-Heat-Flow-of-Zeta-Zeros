#!/usr/bin/env python3
"""Independently validate the arbitrary-N stable remainder theorem."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_theta_arbitrary_n_stable_remainder_gate"
DEFAULT_ARTIFACT = (
    REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
)
MODULAR_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_modular_blend_gate.json"
)
FINITE_PREDECESSOR = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_forward_remainder_tail_gate.json"
)
EXPECTED_IDS = [
    "ntansrg_01_all_n_modular_identity",
    "ntansrg_02_arbitrary_cap_split",
    "ntansrg_03_cap_free_split",
    "ntansrg_04_forward_derivative_recurrence",
    "ntansrg_05_uniform_endpoint_domination",
    "ntansrg_06_uniform_geometric_ratio",
    "ntansrg_07_explicit_arbitrary_tail",
    "ntansrg_08_adaptive_cofinal_handoff",
]
MAX_ORDER = 9
WITNESS_K = [2, 3, 4, 5, 8, 11, 13, 21, 33, 65]
X_SYMBOL = sp.symbols("X", nonnegative=True)


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def parse_exact(text: str) -> sp.Expr:
    return sp.sympify(
        text,
        locals={"X": X_SYMBOL, "pi": sp.pi, "exp": sp.exp},
    )


def validate(path: Path) -> list[str]:
    issues: list[str] = []
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("kind") != STEM:
        issues.append("artifact kind mismatch")

    rows = artifact.get("rows", [])
    ids = [row.get("id") for row in rows]
    if ids != EXPECTED_IDS:
        issues.append(f"row id/order mismatch: {ids}")
    expected_readiness = ["proved"] * 7 + ["not_ready_to_apply"]
    if [row.get("readiness") for row in rows] != expected_readiness:
        issues.append("row readiness mismatch")

    parameters = artifact.get("parameters", {})
    expected_parameters = {
        "max_derivative_order": MAX_ORDER,
        "minimum_retained_count": 1,
        "minimum_first_omitted": 2,
        "retained_count_upper_bound": None,
        "cap_upper_bound": None,
        "witness_K": WITNESS_K,
    }
    if parameters != expected_parameters:
        issues.append(f"parameter mismatch: {parameters}")

    modular = json.loads(MODULAR_SOURCE.read_text(encoding="utf-8"))
    exact = modular.get("exact", {})
    identity = exact.get("decaying_tail_enclosure", {}).get("tail_kernel")
    series = exact.get("theta_summand", {}).get("kernel_series", "")
    if identity != (
        "r_N(u)=sum_(n>N)b_n(u)=Phi(u)-sum_(n<=N)b_n(u)"
    ):
        issues.append("upstream modular remainder identity mismatch")
    if "Phi(u)=sum_(n>=1)phi_n(u)" not in series:
        issues.append("upstream forward theta series mismatch")

    audit = artifact.get("source_audit", {})
    modular_audit = audit.get("modular_source", {})
    predecessor_audit = audit.get("finite_predecessor", {})
    if modular_audit.get("sha256") != file_hash(MODULAR_SOURCE):
        issues.append("modular source hash mismatch")
    if modular_audit.get("tail_identity") != identity:
        issues.append("stored modular identity mismatch")
    if modular_audit.get("forward_series") != series:
        issues.append("stored forward series mismatch")
    if predecessor_audit.get("sha256") != file_hash(FINITE_PREDECESSOR):
        issues.append("finite predecessor hash mismatch")
    predecessor = json.loads(
        FINITE_PREDECESSOR.read_text(encoding="utf-8")
    )
    if predecessor_audit.get("kind") != predecessor.get("kind"):
        issues.append("finite predecessor kind mismatch")
    if predecessor_audit.get("first_omitted") != 13:
        issues.append("finite predecessor first-omitted mismatch")

    polynomials = artifact.get("polynomials", [])
    witnesses = artifact.get("witness_bounds", [])
    if len(polynomials) != MAX_ORDER + 1:
        issues.append("polynomial row count mismatch")
    if len(witnesses) != len(WITNESS_K) * (MAX_ORDER + 1):
        issues.append("witness row count mismatch")

    x = X_SYMBOL
    current = 2 * x - 3
    coefficient_norms: list[int] = []
    for q in range(MAX_ORDER + 1):
        if q >= len(polynomials):
            continue
        polynomial = sp.Poly(sp.expand(current), x)
        stored = polynomials[q]
        coefficient_norm = int(
            sum(abs(value) for value in polynomial.all_coeffs())
        )
        coefficient_norms.append(coefficient_norm)
        if stored.get("q") != q:
            issues.append(f"polynomial order mismatch at q={q}")
        if stored.get("degree") != q + 1:
            issues.append(f"polynomial degree mismatch at q={q}")
        if stored.get("A_q") != coefficient_norm:
            issues.append(f"coefficient norm mismatch at q={q}")
        try:
            if sp.expand(parse_exact(stored["P_q"]) - polynomial.as_expr()) != 0:
                issues.append(f"polynomial recurrence drift at q={q}")
        except Exception as exc:
            issues.append(f"polynomial parse failed at q={q}: {exc}")
        current = sp.expand(
            (5 - 4 * x) * polynomial.as_expr()
            + 4 * x * sp.diff(polynomial.as_expr(), x)
        )

    # Exact endpoint reductions valid for all K>=2 and q<=9.
    if not bool(16 * sp.pi - 45 > 0):
        issues.append("uniform u-monotonicity endpoint failed")
    if not bool(11 - 5 * sp.pi < -4):
        issues.append("uniform ratio endpoint failed")
    if artifact.get("endpoint_extrema") != {
        "u_monotonicity_worst_case": "16*pi-45>48-45>0",
        "ratio_exponent_worst_case": "11-5*pi<11-15=-4",
    }:
        issues.append("stored endpoint extrema mismatch")

    for index, witness in enumerate(witnesses):
        expected_k = WITNESS_K[index // (MAX_ORDER + 1)]
        expected_q = index % (MAX_ORDER + 1)
        if witness.get("K") != expected_k or witness.get("q") != expected_q:
            issues.append(f"witness ordering mismatch at row {index}")
            continue
        q = expected_q
        k = expected_k
        d = 2 * q + 4
        rho = sp.exp(sp.Rational(d, k) - sp.pi * (2 * k + 1))
        bound = (
            coefficient_norms[q]
            * sp.pi ** (q + 2)
            * k**d
            * sp.exp(-sp.pi * k**2)
            / (1 - rho)
        )
        try:
            if sp.simplify(parse_exact(witness["rho_exact"]) - rho) != 0:
                issues.append(f"rho drift at K={k}, q={q}")
            if sp.simplify(parse_exact(witness["B_qK_exact"]) - bound) != 0:
                issues.append(f"tail bound drift at K={k}, q={q}")
            decimal = sp.Float(witness["B_qK_decimal"], 80)
            relative_error = abs(sp.N(decimal / bound - 1, 70))
            if not bool(relative_error < sp.Float("1e-28")):
                issues.append(f"tail decimal drift at K={k}, q={q}")
            if not bool(sp.N(rho, 80) < 1):
                issues.append(f"rho is not below one at K={k}, q={q}")
        except Exception as exc:
            issues.append(f"witness parse failed at K={k}, q={q}: {exc}")

    # Recheck the chain-rule recurrence independently of stored polynomials.
    u = sp.symbols("u", real=True)
    p = sp.Function("P")
    x_of_u = sp.exp(4 * u)
    lhs = sp.diff(sp.exp(5 * u - x_of_u) * p(x_of_u), u)
    rhs = sp.exp(5 * u - x_of_u) * (
        (5 - 4 * x_of_u) * p(x_of_u)
        + 4 * x_of_u * sp.Subs(sp.diff(p(x), x), x, x_of_u)
    )
    if sp.simplify(lhs.doit() - rhs.doit()) != 0:
        issues.append("generic derivative chain rule failed")

    stable = artifact.get("stable_remainder", "")
    for marker in (
        "sum_(n=1)^N D^q delta_n(u)",
        "sum_(n=N+1)^M D^q phi_n(u)",
        "T_(M+1,q)(u)",
        "N>=1, M>=N, K=M+1",
    ):
        if marker not in stable:
            issues.append(f"arbitrary-cap marker missing: {marker}")
    cap_free = artifact.get("cap_free_remainder", "")
    for marker in ("sum_(n=1)^N D^q delta_n(u)", "T_(N+1,q)(u)"):
        if marker not in cap_free:
            issues.append(f"cap-free marker missing: {marker}")

    expected_tail_formula = (
        "B_(q,K)=A_q*pi^(q+2)*K^(2q+4)*exp(-pi*K^2)"
        "/(1-exp((2q+4)/K-pi*(2K+1)))"
    )
    if artifact.get("tail_bound_formula") != expected_tail_formula:
        issues.append("tail formula mismatch")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "switch-defect",
        "full d0/d1 budgets",
        "first-jet separation",
        "x>38",
        "Lambda<=0",
        "RH",
        "Clay-prize",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary guard missing: {marker}")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    args = parser.parse_args()
    issues = validate(args.artifact)
    if issues:
        for issue in issues:
            print(f"ISSUE: {issue}")
        raise SystemExit(1)
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    print(
        "validated Newman theta arbitrary-N stable remainder gate: "
        f"{len(artifact['rows'])} rows, 0 issues, "
        f"{len(artifact['polynomials'])} derivative polynomials, "
        f"{len(artifact['witness_bounds'])} witness bounds, "
        "1 exact all-N split, 1 cofinal handoff"
    )


if __name__ == "__main__":
    main()
