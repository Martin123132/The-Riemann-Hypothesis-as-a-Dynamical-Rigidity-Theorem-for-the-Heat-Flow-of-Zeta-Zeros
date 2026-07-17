#!/usr/bin/env python3
"""Independently validate the Newman strict-Laguerre correlation target."""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import sympy as sp

import jensen_window_pf_newman_strict_laguerre_correlation_target as target


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_strict_laguerre_correlation_target.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_strict_laguerre_correlation_target.md"


EXPECTED_IDS = [
    "nslc_01_boundary_laguerre_zero",
    "nslc_02_strict_laguerre_equivalence",
    "nslc_03_correlation_identity",
    "nslc_04_wiener_density",
    "nslc_05_density_rh_equivalence",
    "nslc_06_positive_definite_scope",
    "nslc_07_strict_logconcavity_guard",
    "nslc_08_generic_kernel_nonpromotion",
    "nslc_09_first_correlation_target",
    "nslc_10_xi_density_handoff",
]


def validate() -> list[str]:
    issues: list[str] = []
    if not RESULT.exists():
        return ["missing stored result"]
    if not NOTE.exists():
        return ["missing rendered note"]

    stored = json.loads(RESULT.read_text(encoding="utf-8"))
    rebuilt = target.build_payload()
    if stored != rebuilt:
        issues.append("stored strict-Laguerre payload differs from reconstruction")
    if stored.get("status") != (
        "exact strict-Laguerre/Wiener equivalence with generic guards and a "
        "retired strict-monotonicity subroute"
    ):
        issues.append("status drifted")

    rows = stored.get("rows", [])
    if len(rows) != 10:
        issues.append("expected 10 target rows")
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    for role, expected in [
        ("exact_equivalence", 2),
        ("exact_identity", 1),
        ("published_theorem", 1),
        ("exact_countermodel", 1),
        ("nonpromotion_gate", 2),
        ("open_theorem_target", 1),
        ("open_handoff", 1),
    ]:
        if sum(row.get("role") == role for row in rows) != expected:
            issues.append(f"expected {expected} row(s) with role {role}")

    exact = stored.get("exact", {})
    definitions = exact.get("definitions", {})
    if "H_t'(x)^2-H_t(x)*H_t''(x)" not in definitions.get("first_laguerre", ""):
        issues.append("first Laguerre definition drifted")
    if "0<t<=1/5" not in exact.get("strict_laguerre_equivalence", ""):
        issues.append("positive-time interval missing")
    if "Platt-Trudgian Corollary 2" not in exact.get(
        "published_time_window_reduction", ""
    ):
        issues.append("published time-window reduction missing")
    if "finite multiple zero" not in exact.get("strict_laguerre_equivalence", ""):
        issues.append("boundary-attainment dependency missing")
    if "L_t(0)>0" not in exact.get("origin_sign", ""):
        issues.append("origin sign anchor missing")

    # Re-derive the midpoint/difference normalization of the correlation identity.
    a, b, u, v = sp.symbols("a b u v", real=True)
    laguerre_pair_weight = (a - b) ** 2 / 2
    transformed_weight = sp.expand(
        laguerre_pair_weight.subs({a: u + v, b: u - v}) * 2
    )
    if sp.simplify(transformed_weight - 4 * v**2) != 0:
        issues.append("midpoint/difference Jacobian normalization failed")
    if exact.get("correlation_normalization", {}).get("full_laguerre") != "L_1[F_t](x)=4*L_t(x)":
        issues.append("full-to-half transform Laguerre factor failed")
    if "Fourier[K_(1,t)](xi)=L_t(xi/2)" not in exact.get("correlation_identity", ""):
        issues.append("correlation Fourier scaling missing")

    wiener = exact.get("wiener_density", {})
    if wiener.get("source") != "https://arxiv.org/abs/1606.05011":
        issues.append("Wiener/correlation source drifted")
    if "if and only if" not in wiener.get("statement", ""):
        issues.append("Wiener equivalence missing")
    if "L_t(0)>0" not in wiener.get("application", ""):
        issues.append("zero-free to strict-positive step missing")
    if "if and only if" not in exact.get("density_rh_equivalence", ""):
        issues.append("density RH equivalence missing")
    published = exact.get("published_correlation_source", {})
    if published.get("source") != "https://arxiv.org/abs/1309.0055":
        issues.append("Csordas correlation source drifted")
    if "Theorem 3.7" not in published.get("theorem", ""):
        issues.append("Csordas theorem number missing")

    # Independently verify the exact triangular-Gaussian Fourier guard.
    a_s, sigma_s, xi = sp.symbols("a sigma xi", positive=True)
    u_transform = 2 * sp.sin(a_s * xi / 2) / xi
    triangle_transform = sp.expand(u_transform**2)
    gaussian_transform = sp.sqrt(2 * sp.pi) * sigma_s * sp.exp(
        -sigma_s**2 * xi**2 / 2
    )
    guard_transform = sp.factor(triangle_transform * gaussian_transform)
    stored_guard = exact.get("strict_logconcavity_guard", {})
    parsed_guard_transform = sp.sympify(
        stored_guard.get("fourier_transform", "nan"),
        locals={"a": a_s, "sigma": sigma_s, "xi": xi},
    )
    if sp.simplify(
        parsed_guard_transform - guard_transform
    ) != 0:
        issues.append("stored guard Fourier transform failed")

    concrete = sp.simplify(guard_transform.subs({a_s: 1, sigma_s: 2}))
    for n in (1, 2, 3):
        zero = 2 * sp.pi * n
        if sp.simplify(concrete.subs(xi, zero)) != 0:
            issues.append(f"guard value failed at n={n}")
        if sp.simplify(sp.diff(concrete, xi).subs(xi, zero)) != 0:
            issues.append(f"guard slope failed at n={n}")
        second = sp.simplify(sp.limit(sp.diff(concrete, xi, 2), xi, zero))
        if second == 0:
            issues.append(f"guard zero is not exactly double at n={n}")
    if sp.limit(concrete, xi, 0, dir="+") != 2 * sp.sqrt(2 * sp.pi):
        issues.append("guard transform origin limit failed")

    # Check the conditional-variance log-curvature identity numerically.
    mp.mp.dps = 60
    sigma_value = mp.mpf(2)
    curvature_upper = -mp.mpf(3) / 16
    curvature_rows: list[dict] = []
    for x_value in [mp.mpf(k) / 2 for k in range(-12, 13)]:
        def unnormalized(y: mp.mpf) -> mp.mpf:
            return (1 - abs(y)) * mp.exp(-(x_value - y) ** 2 / 8)

        mass = mp.quad(unnormalized, [-1, 0, 1])
        mean = mp.quad(lambda y: y * unnormalized(y), [-1, 0, 1]) / mass
        second_moment = (
            mp.quad(lambda y: y**2 * unnormalized(y), [-1, 0, 1]) / mass
        )
        variance = second_moment - mean**2
        curvature = variance / sigma_value**4 - 1 / sigma_value**2
        if curvature > curvature_upper + mp.mpf("1e-45"):
            issues.append(f"guard curvature bound failed at x={x_value}")
        curvature_rows.append(
            {
                "x": str(x_value),
                "variance": mp.nstr(variance, 20),
                "curvature": mp.nstr(curvature, 20),
            }
        )
    if stored_guard.get("curvature_bound_a1_sigma2") != "-3/16":
        issues.append("stored concrete curvature bound failed")
    if "positive definite" not in stored_guard.get("positive_definite", ""):
        issues.append("guard positive-definiteness conclusion missing")
    if "non-dense" not in stored_guard.get("translation_failure", ""):
        issues.append("guard Wiener failure missing")
    if "does not reproduce the Xi tail" not in stored_guard.get("scope_warning", ""):
        issues.append("countermodel tail boundary missing")

    decision = exact.get("nonpromotion_decision", "")
    for phrase in (
        "strict log-concavity",
        "positive definiteness",
        "theta-type",
        "arithmetic or modular",
        "total-positivity",
    ):
        if phrase not in decision:
            issues.append(f"nonpromotion decision missing: {phrase}")
    if "every 0<t<=1/5" not in exact.get("open_handoff", ""):
        issues.append("uniform time handoff missing")
    for phrase in (
        "M_t(x)=-L_t'(x)>0",
        "M_0(1401016343/100000)<0",
        "corrected C1 double-zero transversality",
        "do not impose global monotonicity",
    ):
        if phrase not in exact.get("open_handoff", ""):
            issues.append(f"monotonicity handoff missing: {phrase}")

    note = NOTE.read_text(encoding="utf-8")
    required = [
        "not a proof of RH or `Lambda <= 0`",
        "There is no missing normalization factor",
        "Wiener's theorem",
        "positive definite",
        "strictly log-concave",
        "zeros are double",
        "not dense",
        "Xi-specific tail",
        "theta-tail weighted countermodel",
        "strict_laguerre_monotonicity_scout.md",
        "rigorous Xi counterexample",
        "https://arxiv.org/abs/1606.05011",
    ]
    for phrase in required:
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
        "validated Jensen-window PF Newman strict-Laguerre correlation target: "
        "10 rows, 0 issues, 1 strict-Laguerre equivalence, "
        "1 exact correlation identity, 1 Wiener-density equivalence, "
        "1 RH-equivalent density target, "
        "1 exact strict-log-concavity/positive-definiteness countermodel, "
        "2 non-promotion gates, 1 open Xi handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
