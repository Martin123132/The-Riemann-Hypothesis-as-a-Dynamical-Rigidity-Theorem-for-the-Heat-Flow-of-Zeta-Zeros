#!/usr/bin/env python3
"""Independently validate the Newman Gasper fake-Xi remainder gate."""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import sympy as sp

import jensen_window_pf_newman_gasper_fake_xi_remainder_gate as target


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_gasper_fake_xi_remainder_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_gasper_fake_xi_remainder_gate.md"

EXPECTED_IDS = [
    "gfrr_01_fake_xi_normalization",
    "gfrr_02_gasper_real_zero_benchmark",
    "gfrr_03_laguerre_remainder_identity",
    "gfrr_04_scalar_budget_minimum",
    "gfrr_05_scalar_budget_witnesses",
    "gfrr_06_kernel_ratio_endpoints",
    "gfrr_07_positive_cosh_obstruction",
    "gfrr_08_sign_aware_handoff",
]


def relative_error(left: mp.mpf, right: mp.mpf) -> mp.mpf:
    return abs(left - right) / max(abs(right), mp.mpf("1e-300"))


def validate() -> list[str]:
    issues: list[str] = []
    if not RESULT.exists():
        return ["missing stored result"]
    if not NOTE.exists():
        return ["missing rendered note"]

    stored = json.loads(RESULT.read_text(encoding="utf-8"))
    rebuilt = target.build_payload()
    if stored != rebuilt:
        issues.append("stored Gasper fake-Xi payload differs from reconstruction")
    if stored.get("status") != (
        "exact fake-Xi comparison with scalar and positive-convolution obstructions"
    ):
        issues.append("status drifted")

    rows = stored.get("rows", [])
    if len(rows) != 8:
        issues.append("expected 8 rows")
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    expected_roles = {
        "exact_identity": 2,
        "established_theorem": 1,
        "exact_theorem": 2,
        "interval_counter_certificate": 1,
        "exact_route_obstruction": 1,
        "open_handoff": 1,
    }
    for role, expected in expected_roles.items():
        if sum(row.get("role") == role for row in rows) != expected:
            issues.append(f"expected {expected} row(s) with role {role}")

    exact = stored.get("exact", {})

    # Re-derive the fake-Xi kernel and Bessel scaling under v=4u.
    u = sp.symbols("u", real=True)
    psi = 4 * sp.pi**2 * sp.exp(-2 * sp.pi * sp.cosh(4 * u)) * sp.cosh(9 * u)
    expanded = (
        2
        * sp.pi**2
        * sp.exp(9 * u - sp.pi * sp.exp(4 * u) - sp.pi * sp.exp(-4 * u))
        * (1 + sp.exp(-18 * u))
    )
    if sp.simplify(psi.rewrite(sp.exp) - expanded) != 0:
        issues.append("fake-Xi kernel expansion failed")
    benchmark = exact.get("fake_xi_benchmark", {})
    if benchmark.get("normalization") != "P_0(z)=Xi_star(z/2)/8":
        issues.append("fake-Xi normalization drifted")
    if "pi^2/2" not in benchmark.get("bessel_form", ""):
        issues.append("Bessel normalization missing")
    if "only real zeros" not in benchmark.get("established_result", ""):
        issues.append("Gasper benchmark theorem missing")

    # Re-derive the quadratic self/cross/residual identities.
    h, hp, hpp, p, pp, ppp, alpha = sp.symbols(
        "h hp hpp p pp ppp alpha", real=True
    )
    l_h = hp**2 - h * hpp
    l_p = pp**2 - p * ppp
    b_ph = 2 * pp * hp - p * hpp - h * ppp
    e, ep, epp = h - alpha * p, hp - alpha * pp, hpp - alpha * ppp
    l_e = ep**2 - e * epp
    b_ape = 2 * alpha * pp * ep - alpha * p * epp - e * alpha * ppp
    if sp.simplify(l_e - (l_h - alpha * b_ph + alpha**2 * l_p)) != 0:
        issues.append("residual Laguerre identity failed")
    if sp.simplify(b_ape - (alpha * b_ph - 2 * alpha**2 * l_p)) != 0:
        issues.append("mixed Laguerre identity failed")
    if sp.simplify(l_h - (alpha**2 * l_p + b_ape + l_e)) != 0:
        issues.append("full Laguerre decomposition failed")

    # Re-derive the exact scalar minimization in its rho<2 sign chamber.
    a, b, lp, lh = sp.symbols("a b lp lh", positive=True)
    rho = b**2 / (lp * lh)
    branch = 3 - 2 * b / (a * lp) + lh / (a**2 * lp)
    square = (a * b - lh) ** 2 / (a**2 * lp * lh)
    if sp.simplify(branch - (3 - rho) - square) != 0:
        issues.append("scalar minimum square completion failed")
    left_boundary = 4 * lh * lp / b**2 - 1
    left_gap = (b**2 - 2 * lh * lp) ** 2 / (b**2 * lp * lh)
    if sp.simplify(left_boundary - (3 - rho) - left_gap) != 0:
        issues.append("scalar alternate-branch comparison failed")
    scalar = exact.get("scalar_absolute_budget", {})
    if "inf_(alpha>0) R_alpha=3-rho>1" not in scalar.get("theorem", ""):
        issues.append("stored scalar theorem drifted")
    if scalar.get("optimizer") != "alpha_star=L[H]/B[P,H]":
        issues.append("stored scalar optimizer drifted")
    if "rho<2 makes L[E_alpha]>0" not in scalar.get("branch_check", ""):
        issues.append("stored scalar branch proof missing")

    # Validate the rigorous endpoint-time interval certificates.
    interval = stored.get("numerics", {}).get("interval_certificate", {})
    if interval.get("rigorous") is not True:
        issues.append("scalar interval certificate is not marked rigorous")
    if interval.get("cell_count") != 2048 or interval.get("theta_summands") != 8:
        issues.append("scalar interval certificate geometry drifted")
    if interval.get("omitted_tail_radius_per_jet") != "1e-45":
        issues.append("scalar interval tail cap drifted")
    interval_rows = interval.get("rows", [])
    if [row.get("t") for row in interval_rows] != ["0", "1/2"]:
        issues.append("scalar interval endpoint times drifted")
    target.ctx.dps = 70
    for row in interval_rows:
        try:
            l_h_ball = target.arb(row.get("L_H"))
            l_p_ball = target.arb(row.get("L_P"))
            b_ball = target.arb(row.get("B_PH"))
            rho_ball = target.arb(row.get("rho"))
            minimum_ball = target.arb(row.get("minimum_absolute_budget"))
        except (TypeError, ValueError):
            issues.append(f"unparseable scalar interval row at t={row.get('t')}")
            continue
        if not (l_h_ball > 0 and l_p_ball > 0 and b_ball > 0):
            issues.append(f"scalar interval sign chamber failed at t={row.get('t')}")
        if not (rho_ball > 0 and rho_ball < 2):
            issues.append(f"scalar interval rho enclosure failed at t={row.get('t')}")
        if not minimum_ball > 1:
            issues.append(f"scalar interval budget failed at t={row.get('t')}")
        if not all(row.get("certified_signs", {}).values()):
            issues.append(f"stored scalar interval flags failed at t={row.get('t')}")
    tail_proof = interval.get("tail_proof", "")
    for phrase in ("n>=9", "u=1+v>=1", "exp(-587v)", "1e-60", "1e-45"):
        if phrase not in tail_proof:
            issues.append(f"scalar interval tail proof missing: {phrase}")

    # Independently recompute the selected jets with a third quadrature order.
    independent = target.numerical_rows(520)
    stored_rows = stored.get("numerics", {}).get("rows", [])
    if len(stored_rows) != 2 or len(independent) != 2:
        issues.append("expected two endpoint-time numerical rows")
    fields = ("L_H", "L_P", "B_PH", "rho", "alpha_opt", "minimum_absolute_budget")
    for observed, check in zip(stored_rows, independent, strict=False):
        if observed.get("t") != check.get("t") or observed.get("x") != "25":
            issues.append("numerical witness coordinate drifted")
            continue
        for field in fields:
            left = mp.mpf(observed.get(field, "nan"))
            right = mp.mpf(check.get(field, "nan"))
            if relative_error(left, right) > mp.mpf("2e-10"):
                issues.append(f"independent quadrature mismatch: t={observed.get('t')} {field}")
        rho_value = mp.mpf(observed.get("rho", "nan"))
        minimum = mp.mpf(observed.get("minimum_absolute_budget", "nan"))
        if not (0 < rho_value < 2):
            issues.append(f"rho sign chamber failed at t={observed.get('t')}")
        if minimum <= mp.mpf("2.5"):
            issues.append(f"scalar obstruction margin failed at t={observed.get('t')}")

    # Compare the generated fast quadrature against separately computed 80-digit values.
    references = stored.get("numerics", {}).get("independent_mpmath_80_digit_references", {})
    by_time = {"0": "t=0", "0.5": "t=1/2"}
    for observed in stored_rows:
        reference = references.get(by_time.get(observed.get("t", ""), ""), {})
        if not reference:
            issues.append(f"missing mpmath reference for t={observed.get('t')}")
            continue
        for field in fields:
            expected = mp.mpf(reference.get(field, "nan"))
            actual = mp.mpf(observed.get(field, "nan"))
            if relative_error(actual, expected) > mp.mpf("5e-11"):
                issues.append(f"mpmath reference mismatch: t={observed.get('t')} {field}")

    # Check the exact endpoint comparison behind the cosh-mixture obstruction.
    mp.mp.dps = 60
    phi_one_zero = mp.pi * (2 * mp.pi - 3) * mp.e ** (-mp.pi)
    psi_zero = 4 * mp.pi**2 * mp.e ** (-2 * mp.pi)
    if not phi_one_zero > psi_zero:
        issues.append("origin kernel inequality failed")
    ratio = exact.get("kernel_ratio", {})
    if "M(0)>1" not in ratio.get("origin", ""):
        issues.append("origin ratio theorem missing")
    if ratio.get("limit") != "lim_(u->infinity) M(u)=1":
        issues.append("kernel-ratio tail limit drifted")
    obstruction = exact.get("positive_cosh_obstruction", {})
    for phrase in ("finite limit", "mu(R\\{0})=0", "constant"):
        if phrase not in obstruction.get("boundedness_lemma", ""):
            issues.append(f"positive-cosh boundedness lemma missing: {phrase}")
    if "direct Cardon/Polya" not in obstruction.get("conclusion", ""):
        issues.append("positive-convolution rejection missing")
    if "signed measures" not in obstruction.get("scope", ""):
        issues.append("positive-convolution proof boundary missing")

    note = NOTE.read_text(encoding="utf-8")
    required_note = (
        "inf R_alpha=2.5852035",
        "Arb midpoint",
        "M(0)>1",
        "lim_(u->infinity) M(u)=1",
        "sign-aware",
        "This is not a proof of RH",
    )
    for phrase in required_note:
        if phrase not in note:
            issues.append(f"rendered note missing: {phrase}")
    for forbidden in ("therefore RH", "we have proved RH", "Lambda <= 0 is proved"):
        if forbidden.lower() in note.lower():
            issues.append(f"rendered note contains forbidden promotion: {forbidden}")
    return issues


def main() -> int:
    issues = validate()
    for issue in issues:
        print(f"GASPER-FAKE-XI [{issue}]")
    print(
        "validated Jensen-window PF Newman Gasper fake-Xi remainder gate: "
        f"8 rows, {len(issues)} issues, 2 exact transform identities, "
        "1 established real-zero benchmark, 1 scalar algebra theorem, "
        "2 interval scalar witnesses, 2 high-precision cross-checks, "
        "1 exact positive-convolution obstruction, "
        "1 sign-aware handoff"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
