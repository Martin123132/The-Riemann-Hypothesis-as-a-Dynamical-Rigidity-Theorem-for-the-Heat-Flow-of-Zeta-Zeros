#!/usr/bin/env python3
"""Independently validate the stable forward-remainder tail gate."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_theta_forward_remainder_tail_gate"
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
EXPECTED_IDS = [
    "ntfrtg_01_stable_modular_remainder",
    "ntfrtg_02_forward_derivative_recurrence",
    "ntfrtg_03_uniform_forward_tail",
    "ntfrtg_04_geometric_arithmetic_sum",
    "ntfrtg_05_compact_quadratic_correction",
    "ntfrtg_06_stable_matrix_handoff",
]
MAX_ORDER = 9
FIRST_OMITTED = 13
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
    expected_readiness = ["proved"] * 5 + ["not_ready_to_apply"]
    readiness = [row.get("readiness") for row in rows]
    if readiness != expected_readiness:
        issues.append(f"readiness mismatch: {readiness}")
    if any(row.get("readiness") == "promoted" for row in rows):
        issues.append("unexpected promoted row")

    parameters = artifact.get("parameters", {})
    if parameters.get("max_derivative_order") != MAX_ORDER:
        issues.append("maximum derivative order mismatch")
    if parameters.get("forward_cap") != FIRST_OMITTED - 1:
        issues.append("forward cap mismatch")
    if parameters.get("first_omitted") != FIRST_OMITTED:
        issues.append("first omitted index mismatch")
    if parameters.get("retained_range") != ["1", "12"]:
        issues.append("retained range mismatch")
    if parameters.get("compact_interval") != ["0", "11/5"]:
        issues.append("compact interval mismatch")

    modular = json.loads(MODULAR_SOURCE.read_text(encoding="utf-8"))
    exact = modular.get("exact", {})
    identity = exact.get("decaying_tail_enclosure", {}).get("tail_kernel")
    series = exact.get("theta_summand", {}).get("kernel_series", "")
    expected_identity = (
        "r_N(u)=sum_(n>N)b_n(u)=Phi(u)-sum_(n<=N)b_n(u)"
    )
    if identity != expected_identity:
        issues.append("upstream modular remainder identity mismatch")
    if "Phi(u)=sum_(n>=1)phi_n(u)" not in series:
        issues.append("upstream forward theta series mismatch")
    audit = artifact.get("source_audit", {})
    if audit.get("sha256") != file_hash(MODULAR_SOURCE):
        issues.append("upstream source hash mismatch")
    if audit.get("tail_identity") != identity:
        issues.append("stored source identity mismatch")
    if audit.get("forward_series") != series:
        issues.append("stored source series mismatch")

    x = X_SYMBOL
    current = 2 * x - 3
    stored_polynomials = artifact.get("polynomials", [])
    stored_tails = artifact.get("tail_bounds", [])
    if len(stored_polynomials) != MAX_ORDER + 1:
        issues.append("polynomial row count mismatch")
    if len(stored_tails) != MAX_ORDER + 1:
        issues.append("tail row count mismatch")

    for q in range(MAX_ORDER + 1):
        if q >= len(stored_polynomials) or q >= len(stored_tails):
            continue
        polynomial = sp.Poly(sp.expand(current), x)
        stored_polynomial = stored_polynomials[q]
        if stored_polynomial.get("q") != q:
            issues.append(f"polynomial order mismatch at q={q}")
        try:
            stored_expression = parse_exact(stored_polynomial["P_q"])
        except Exception as exc:
            issues.append(f"polynomial parse failed at q={q}: {exc}")
            stored_expression = sp.nan
        if sp.expand(stored_expression - polynomial.as_expr()) != 0:
            issues.append(f"polynomial recurrence drift at q={q}")
        coefficient_norm = int(
            sum(abs(value) for value in polynomial.all_coeffs())
        )
        if stored_polynomial.get("degree") != q + 1:
            issues.append(f"polynomial degree mismatch at q={q}")
        if stored_polynomial.get("A_q") != coefficient_norm:
            issues.append(f"coefficient norm mismatch at q={q}")

        tail = stored_tails[q]
        d = 2 * q + 4
        a = 4 * q + 9
        if tail.get("q") != q or tail.get("d_q") != d:
            issues.append(f"tail exponent mismatch at q={q}")
        if tail.get("a_q") != a:
            issues.append(f"u exponent mismatch at q={q}")
        if not bool(4 * sp.pi * FIRST_OMITTED**2 - a > 0):
            issues.append(f"u monotonicity failed at q={q}")

        rho = sp.exp(
            sp.Rational(d, FIRST_OMITTED)
            - sp.pi * (2 * FIRST_OMITTED + 1)
        )
        try:
            stored_rho = parse_exact(tail["rho_exact"])
            if sp.simplify(stored_rho - rho) != 0:
                issues.append(f"geometric ratio drift at q={q}")
        except Exception as exc:
            issues.append(f"rho parse failed at q={q}: {exc}")
        if not bool(sp.N(rho, 80) < 1):
            issues.append(f"geometric ratio is not below one at q={q}")

        expected_bound = (
            coefficient_norm
            * sp.pi ** (q + 2)
            * FIRST_OMITTED ** (2 * q + 4)
            * sp.exp(-sp.pi * FIRST_OMITTED**2)
            / (1 - rho)
        )
        try:
            stored_bound = parse_exact(tail["B_q_exact"])
            if sp.simplify(stored_bound - expected_bound) != 0:
                issues.append(f"tail bound drift at q={q}")
            stored_decimal = sp.Float(tail["B_q_decimal"], 80)
            relative_error = abs(
                sp.N(stored_decimal / expected_bound - 1, 70)
            )
            if not bool(relative_error < sp.Float("1e-28")):
                issues.append(f"tail decimal drift at q={q}")
        except Exception as exc:
            issues.append(f"tail bound parse failed at q={q}: {exc}")

        current = sp.expand(
            (5 - 4 * x) * polynomial.as_expr()
            + 4 * x * sp.diff(polynomial.as_expr(), x)
        )

    # Verify the generic chain-rule step independently of the stored rows.
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
        "erfc(3*sinh(4u))/2",
        "sum_(n=N+1)^12 phi_n(u)",
        "sum_(n=1)^N delta_n(u)",
        "sum_(n>=13)D^q phi_n(u)",
        "1<=N<=12",
    ):
        if marker not in stable:
            issues.append(f"stable remainder marker missing: {marker}")
    if artifact.get("compact_correction") != (
        "sqrt(M_(p,b)Q_f)+B_q M_(p,b)"
    ):
        issues.append("compact correction mismatch")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "outer interval u>11/5",
        "first-jet",
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
        "validated Newman theta forward-remainder tail gate: "
        f"{len(artifact['rows'])} rows, 0 issues, "
        f"{len(artifact['polynomials'])} exact derivative polynomials, "
        f"{len(artifact['tail_bounds'])} explicit arithmetic-tail bounds, "
        "1 open stable-matrix handoff"
    )


if __name__ == "__main__":
    main()
