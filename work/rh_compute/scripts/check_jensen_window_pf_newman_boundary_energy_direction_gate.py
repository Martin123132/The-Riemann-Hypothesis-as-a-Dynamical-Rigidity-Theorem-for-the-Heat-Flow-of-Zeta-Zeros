#!/usr/bin/env python3
"""Independently validate the Newman boundary-energy direction gate."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_boundary_energy_direction_gate as gate


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_boundary_energy_direction_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_boundary_energy_direction_gate.md"


EXPECTED_IDS = [
    "nbeng_01_pair_gap_ode",
    "nbeng_02_stiffness_boundary_jet",
    "nbeng_03_exact_higher_jet_birth",
    "nbeng_04_renormalized_potential",
    "nbeng_05_collision_energy_asymptotic",
    "nbeng_06_forward_relaxation_guard",
    "nbeng_07_boundary_integrability_criterion",
    "nbeng_08_published_energy_scope",
    "nbeng_09_directionality_gate",
    "nbeng_10_xi_boundary_energy_handoff",
]


def validate() -> list[str]:
    issues: list[str] = []
    if not RESULT.exists():
        return ["missing stored result"]
    if not NOTE.exists():
        return ["missing rendered note"]

    stored = json.loads(RESULT.read_text(encoding="utf-8"))
    rebuilt = gate.build_payload()
    if stored != rebuilt:
        issues.append("stored boundary-energy payload differs from reconstruction")
    if stored.get("status") != "exact boundary-energy singularity and published-directionality gate":
        issues.append("status drifted")

    rows = stored.get("rows", [])
    if len(rows) != 10:
        issues.append("expected 10 theorem rows")
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    for role, expected in [
        ("exact_countermodel", 1),
        ("conditional_criterion", 1),
        ("circularity_gate", 1),
        ("open_handoff", 1),
    ]:
        if sum(row.get("role") == role for row in rows) != expected:
            issues.append(f"expected {expected} row with role {role}")

    # Re-derive the squared-gap equation directly from the two root ODEs.
    x_plus, x_minus, outsider = sp.symbols("x_plus x_minus outsider", real=True)
    gap = x_plus - x_minus
    outsider_difference = sp.factor(
        1 / (x_plus - outsider) - 1 / (x_minus - outsider)
    )
    expected_difference = -gap / (
        (x_plus - outsider) * (x_minus - outsider)
    )
    if sp.simplify(outsider_difference - expected_difference) != 0:
        issues.append("independent outsider gap identity failed")
    singular_gap_velocity = 4 / gap
    q_velocity = sp.factor(2 * gap * singular_gap_velocity)
    if q_velocity != 8:
        issues.append("universal squared-gap velocity failed")

    # Reconstruct the exact even birth model and all displayed branch coefficients.
    z, tau = sp.symbols("z tau", real=True)
    a_squared = 1 + 16 / (8 + sp.pi)
    initial = sp.expand((z**2 - 1) ** 2 * (z**2 - a_squared))
    heat = sp.expand(
        sum(
            (-tau) ** order
            / sp.factorial(order)
            * sp.diff(initial, z, 2 * order)
            for order in range(4)
        )
    )
    if sp.simplify(sp.diff(heat, tau) + sp.diff(heat, z, 2)) != 0:
        issues.append("independent polynomial heat equation failed")
    cofactor = sp.cancel(initial / (z - 1) ** 2)
    field = sp.simplify(sp.diff(cofactor, z).subs(z, 1) / cofactor.subs(z, 1))
    if field != -sp.pi / 8:
        issues.append("independent classical field failed")

    stiffness = sp.factor(
        sp.Rational(1, 2)
        + 2 * (1 + a_squared) / (1 - a_squared) ** 2
    )
    expected_stiffness = (sp.pi**2 + 24 * sp.pi + 160) / 64
    if sp.simplify(stiffness - expected_stiffness) != 0:
        issues.append("independent stiffness failed")
    if not bool(sp.N(stiffness, 30) > 0):
        issues.append("model stiffness is not positive")

    eps = sp.symbols("eps", positive=True)
    solved_branches: dict[int, sp.Expr] = {}
    expected_p4 = (sp.pi + 12) * (sp.pi**2 + 24 * sp.pi + 96) / 128
    expected_l5 = (
        15 * sp.pi**4
        + 720 * sp.pi**3
        + 11904 * sp.pi**2
        + 82432 * sp.pi
        + 211968
    ) / 8192
    for sign in (1, -1):
        m, n, p4, r5 = sp.symbols(f"m_{sign} n_{sign} p4_{sign} r5_{sign}")
        trial = (
            1
            + sign * sp.sqrt(2) * eps
            + m * eps**2
            + n * eps**3
            + p4 * eps**4
            + r5 * eps**5
        )
        residual = sp.series(
            heat.subs({tau: eps**2, z: trial}), eps, 0, 7
        ).removeO().expand()
        equations = [
            sp.Eq(sp.simplify(residual.coeff(eps, power)), 0)
            for power in range(3, 7)
        ]
        solutions = sp.solve(equations, (m, n, p4, r5), dict=True, simplify=True)
        if len(solutions) != 1:
            issues.append(f"branch coefficient solve failed for sign={sign}")
            continue
        solution = solutions[0]
        expected = {
            m: -sp.pi / 4,
            n: -sign * sp.sqrt(2) * expected_stiffness,
            p4: expected_p4,
            r5: sign * sp.sqrt(2) * expected_l5,
        }
        for coefficient, value in expected.items():
            if sp.simplify(solution[coefficient] - value) != 0:
                issues.append(
                    f"branch coefficient {coefficient} failed for sign={sign}"
                )
        solved_branches[sign] = sp.expand(trial.subs(solution))

    expected_c3 = (
        sp.pi**4
        + 48 * sp.pi**3
        + 800 * sp.pi**2
        + 5632 * sp.pi
        + 14848
    ) / 32
    if len(solved_branches) == 2:
        q_series = sp.series(
            (solved_branches[1] - solved_branches[-1]) ** 2,
            eps,
            0,
            8,
        ).removeO().expand()
        expected_q = (
            8 * eps**2
            - 16 * expected_stiffness * eps**4
            + expected_c3 * eps**6
        )
        if sp.simplify(q_series - expected_q) != 0:
            issues.append("independent cubic squared-gap jet failed")
    if not bool(sp.N(expected_c3, 30) > 0):
        issues.append("cubic gap coefficient is not positive")

    # Check the renormalized potential and its collision asymptotics independently.
    ratio, reference_gap = sp.symbols("ratio Delta", positive=True)
    potential = ratio ** -2 + 2 * ratio - 3
    if potential.subs(ratio, 1) != 0:
        issues.append("renormalized potential value at equilibrium failed")
    if sp.diff(potential, ratio).subs(ratio, 1) != 0:
        issues.append("renormalized potential slope at equilibrium failed")
    if sp.simplify(sp.diff(potential, ratio, 2) - 6 / ratio**4) != 0:
        issues.append("renormalized potential convexity failed")

    q_model = (
        8 * tau
        - 16 * expected_stiffness * tau**2
        + expected_c3 * tau**3
    )
    gap_model = sp.sqrt(q_model)
    pair_energy = 1 / q_model + 2 * gap_model / reference_gap**3 - 3 / reference_gap**2
    if sp.limit(tau * pair_energy, tau, 0, dir="+") != sp.Rational(1, 8):
        issues.append("collision energy leading coefficient failed")
    derivative = sp.diff(pair_energy, tau)
    if sp.limit(tau**2 * derivative, tau, 0, dir="+") != -sp.Rational(1, 8):
        issues.append("collision energy derivative leading coefficient failed")
    cutoff = sp.symbols("cutoff", positive=True)
    divergent_integral = sp.integrate(1 / (8 * tau), (tau, cutoff, 1))
    if sp.limit(divergent_integral, cutoff, 0, dir="+") != sp.oo:
        issues.append("logarithmic boundary divergence failed")

    exact = stored.get("exact", {})
    scope = exact.get("rodgers_tao_scope", {})
    if scope.get("source") != "https://arxiv.org/abs/1801.05914":
        issues.append("Rodgers-Tao source drifted")
    if "Lambda<0" not in scope.get("contradiction_regime", ""):
        issues.append("negative-Lambda assumption warning missing")
    if "Lambda/4" not in scope.get("integrated_energy_theorem", ""):
        issues.append("Theorem 17 time interval missing")
    if "3*abs(Lambda)/4" not in scope.get("boundary_separation", ""):
        issues.append("boundary separation calculation missing")
    if "unavailable" not in scope.get("directionality", ""):
        issues.append("opposite-direction warning missing")
    if "finite local integrated" not in exact.get("open_handoff", ""):
        issues.append("boundary-energy handoff missing")

    note = NOTE.read_text(encoding="utf-8")
    required = [
        "not a proof of RH or `Lambda <= 0`",
        "first two nontrivial gap-jet signs",
        "nonintegrable trace",
        "Nonnegative forward relaxation",
        "`(Lambda,Lambda+epsilon)`",
        "hypothetical regime `Lambda>0`",
        "boundary-trace theorem",
        "arXiv:1801.05914v5",
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
        "validated Jensen-window PF Newman boundary-energy direction gate: "
        "10 rows, 0 issues, 1 universal gap law, 1 exact higher-jet birth model, "
        "1 nonintegrable collision-energy asymptotic, 1 conditional exclusion criterion, "
        "1 published-scope audit, 1 open Xi boundary-energy handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
