#!/usr/bin/env python3
"""Independently validate the positive-time strong-log-concavity gate."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_positive_time_strong_logconcavity_gate as target
import flint


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_positive_time_strong_logconcavity_gate.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_positive_time_strong_logconcavity_gate.md"
)

EXPECTED_IDS = [
    "npslc_01_published_root_concavity",
    "npslc_02_uniform_phi_curvature",
    "npslc_03_origin_curvature_certificate",
    "npslc_04_positive_time_admissibility",
    "npslc_05_correlation_marginal",
    "npslc_06_zeroth_square_transform",
    "npslc_07_xi_shape_nonpromotion",
    "npslc_08_target_window_margin",
    "npslc_09_weighted_hierarchy_handoff",
]


def parse_ball(text: str) -> flint.arb:
    return flint.arb(text.replace("E", "e"))


def validate() -> list[str]:
    issues: list[str] = []
    if not RESULT.exists():
        return ["missing stored result"]
    if not NOTE.exists():
        return ["missing rendered note"]

    stored = json.loads(RESULT.read_text(encoding="utf-8"))
    rebuilt = target.build_payload()
    if stored != rebuilt:
        issues.append("stored payload differs from reconstruction")
    if stored.get("kind") != (
        "jensen_window_pf_newman_positive_time_strong_logconcavity_gate"
    ):
        issues.append("kind drifted")
    if stored.get("status") != (
        "exact positive-time strong-log-concavity theorem with Xi-specific nonpromotion gate"
    ):
        issues.append("status drifted")

    rows = stored.get("rows", [])
    if len(rows) != 9:
        issues.append("expected nine rows")
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    expected_roles = {
        "published_theorem": 1,
        "exact_lemma": 1,
        "interval_certificate": 1,
        "exact_theorem": 2,
        "published_theorem_composition": 1,
        "exact_identity": 1,
        "nonpromotion_gate": 1,
        "open_handoff": 1,
    }
    for role, count in expected_roles.items():
        if sum(row.get("role") == role for row in rows) != count:
            issues.append(f"expected {count} row(s) with role {role}")

    # Independently derive the two theta-summand origin formulas.
    u, a = sp.symbols("u a", positive=True)
    summand = (
        2 * a**2 * sp.exp(9 * u) - 3 * a * sp.exp(5 * u)
    ) * sp.exp(-a * sp.exp(4 * u))
    expected0 = a * (2 * a - 3) * sp.exp(-a)
    expected2 = a * (32 * a**3 - 224 * a**2 + 330 * a - 75) * sp.exp(-a)
    if sp.simplify(summand.subs(u, 0) - expected0) != 0:
        issues.append("Phi summand value formula failed")
    if sp.simplify(sp.diff(summand, u, 2).subs(u, 0) - expected2) != 0:
        issues.append("Phi summand second-derivative formula failed")

    # Recheck the calculus transfer from q(r)=log Phi(sqrt(r)).
    ell1, ell2 = sp.symbols("ell1 ell2", real=True)
    q_second = (u * ell2 - ell1) / (4 * u**3)
    if sp.factor(4 * u**3 * q_second) != u * ell2 - ell1:
        issues.append("root-variable curvature identity failed")

    certificate = stored.get("certificate", {})
    fresh_certificate = target.curvature_certificate()
    if certificate != fresh_certificate:
        issues.append("stored Arb certificate differs from fresh computation")
    for field in (
        "phi0_ball",
        "phi2_ball",
        "kappa_ball",
        "threshold_kappa_over_2_ball",
        "threshold_margin_ball",
        "target_strong_curvature_lower_ball",
        "correlation_strong_curvature_lower_ball",
    ):
        if field not in certificate:
            issues.append(f"certificate missing {field}")
    if certificate.get("passed") is not True:
        issues.append("Arb certificate not marked passed")
    if certificate.get("tail_start") != 9:
        issues.append("theta-tail start drifted")
    if parse_ball(certificate.get("phi0_ball", "0")) <= 0:
        issues.append("Phi(0) did not certify positive")
    if parse_ball(certificate.get("phi2_ball", "0")) >= 0:
        issues.append("Phi''(0) did not certify negative")
    if parse_ball(certificate.get("threshold_margin_ball", "0")) <= 0:
        issues.append("kappa/2 did not certify beyond 1/5")
    if parse_ball(certificate.get("target_strong_curvature_lower_ball", "0")) <= 0:
        issues.append("target kernel curvature margin did not certify")
    if parse_ball(certificate.get("correlation_strong_curvature_lower_ball", "0")) <= 0:
        issues.append("target correlation curvature margin did not certify")

    # Check the joint Hessian decomposition used before Prekopa marginalization.
    p, q, alpha, beta, n, s = sp.symbols(
        "p q alpha beta n s", real=True, positive=True
    )
    matrix_form = (
        (alpha + beta - 2 * n / s**2) * p**2
        + 2 * (alpha - beta) * p * q
        + (alpha + beta) * q**2
    )
    pair_form = alpha * (p + q) ** 2 + beta * (p - q) ** 2 - 2 * n * p**2 / s**2
    if sp.expand(matrix_form - pair_form) != 0:
        issues.append("correlation Hessian pair decomposition failed")

    exact = stored.get("exact", {})
    published = exact.get("published_input", {})
    if published.get("primary_source") != "https://doi.org/10.1007/BF02075457":
        issues.append("published curvature source drifted")
    correlation = exact.get("correlation", {})
    if correlation.get("marginal_source") != (
        "https://doi.org/10.1016/0022-1236(76)90004-5"
    ):
        issues.append("Prekopa source drifted")
    for phrase in ("-2*(kappa-2*t)", "every n>=0"):
        if phrase not in correlation.get("marginal_theorem", ""):
            issues.append(f"correlation theorem missing {phrase}")

    guard = exact.get("zeroth_correlation_guard", {})
    if guard.get("identity") != "Fourier[K_(0,t)](xi)=2*H_t(xi/2)^2":
        issues.append("zeroth-correlation normalization drifted")
    if "Hardy's theorem" not in guard.get("zero_source", ""):
        issues.append("unconditional real-zero source missing")
    decision = exact.get("nonpromotion_decision", "")
    for phrase in ("K_(0,0)", "s^2 weighting", "do not force Fourier zero-freeness"):
        if phrase not in decision:
            issues.append(f"nonpromotion decision missing {phrase}")
    handoff = exact.get("open_handoff", "")
    for phrase in ("K_(1,t)", "K_(0,t)", "hierarchy coercivity"):
        if phrase not in handoff:
            issues.append(f"weighted hierarchy handoff missing {phrase}")

    note = NOTE.read_text(encoding="utf-8")
    required_note = [
        "not a proof of RH",
        "Uniform Curvature",
        "Arb Threshold",
        "Correlation Propagation",
        "Nonpromotion Gate",
        "K_(0,0)",
        "s^2 weighting",
        "https://doi.org/10.1007/BF02075457",
    ]
    for phrase in required_note:
        if phrase not in note:
            issues.append(f"note missing required phrase: {phrase}")
    return issues


def main() -> int:
    issues = validate()
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated Jensen-window PF Newman positive-time strong-log-concavity gate: "
        "9 rows, 0 issues, 1 published input, 2 exact curvature/admissibility "
        "theorems, 1 Arb threshold certificate, 1 Prekopa correlation theorem, "
        "1 square-transform identity, 1 Xi-specific nonpromotion gate, "
        "1 target-window margin, 1 weighted-hierarchy handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
