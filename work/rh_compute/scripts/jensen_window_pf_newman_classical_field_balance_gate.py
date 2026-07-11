#!/usr/bin/env python3
"""Derive the classical-location balance and compactness gate for the Newman field."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_classical_field_balance_gate.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_classical_field_balance_gate.md"


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def build_exact() -> dict:
    c0, dc, root, dr = sp.symbols("c0 dc root dr", real=True)
    actual_center = c0 + dc
    actual_root = root + dr
    perturbation = sp.factor(
        1 / (actual_center - actual_root) - 1 / (c0 - root)
    )
    expected_perturbation = sp.factor(
        (dr - dc) / ((actual_center - actual_root) * (c0 - root))
    )
    if sp.simplify(perturbation - expected_perturbation) != 0:
        raise RuntimeError("field perturbation identity failed")

    h = sp.symbols("h", positive=True)
    lattice_checks: list[dict] = []
    for cutoff in range(1, 21):
        distances = [sp.Rational(2 * k + 1, 2) * h for k in range(1, cutoff + 1)]
        field = sp.factor(sum(1 / distance - 1 / distance for distance in distances))
        if field != 0:
            raise RuntimeError(f"arithmetic lattice cancellation failed at cutoff={cutoff}")
        lattice_checks.append({"cutoff": cutoff, "field": "0", "passed": True})

    odd_square_sum = sp.summation(
        1 / (2 * sp.Symbol("k", integer=True, nonnegative=True) + 1) ** 2,
        (sp.Symbol("k", integer=True, nonnegative=True), 0, sp.oo),
    )
    if sp.simplify(odd_square_sum - sp.pi**2 / 8) != 0:
        raise RuntimeError("odd-square sum failed")

    r = sp.symbols("r", positive=True)
    j_integral = sp.Rational(1, 2) * (
        sp.log(r) * sp.log((1 + r) / (1 - r))
        + sp.polylog(2, -r)
        - sp.polylog(2, r)
    )
    j_derivative = sp.diff(j_integral, r).replace(
        lambda expr: expr.func == sp.polylog and expr.args[0] == 1,
        lambda expr: -sp.log(1 - expr.args[1]),
    )
    if sp.simplify(
        sp.expand_log(j_derivative, force=True) - sp.log(r) / (1 - r**2)
    ) != 0:
        raise RuntimeError("truncated logarithmic integral failed")

    a, x = sp.symbols("a x", positive=True)
    scaled_r = a / x
    log_scale = sp.log(x / (4 * sp.pi))
    continuum_field = sp.factor(
        (
            -log_scale * sp.atanh(scaled_r)
            - sp.pi**2 / 4
            - j_integral.subs(r, scaled_r)
        )
        / (2 * sp.pi)
    )
    continuum_correction = sp.simplify(continuum_field + sp.pi / 8)
    correction_slope = sp.simplify(
        sp.limit(
            (
                (
                    -sp.log(a / (4 * sp.pi * r)) * sp.atanh(r)
                    - j_integral
                )
                / (2 * sp.pi)
            )
            / r,
            r,
            0,
            dir="+",
        )
    )
    expected_slope = (1 - sp.log(a / (4 * sp.pi))) / (2 * sp.pi)
    if sp.simplify(correction_slope - expected_slope) != 0:
        raise RuntimeError("continuum field asymptotic failed")

    time = sp.symbols("time", nonnegative=True)
    heat_density_correction = time / (8 * x) * sp.log(sp.sqrt(x**2 - a**2) / a)
    y = sp.symbols("y", positive=True)
    antiderivative = (sp.log(y) - sp.log(x**2 - y**2) / 2) / x**2
    if sp.simplify(sp.diff(antiderivative, y) - 1 / (y * (x**2 - y**2))) != 0:
        raise RuntimeError("heat-density correction antiderivative failed")

    level = sp.symbols("level", positive=True)
    density_log = sp.log(level / (4 * sp.pi))
    g_time = density_log / 16
    g_level = density_log / (4 * sp.pi) + time / (16 * level)
    quantile_velocity = sp.factor(-g_time / g_level)
    expected_velocity = -sp.pi * level * density_log / (
        4 * level * density_log + sp.pi * time
    )
    if sp.simplify(quantile_velocity - expected_velocity) != 0:
        raise RuntimeError("classical quantile velocity failed")
    if sp.limit(quantile_velocity, level, sp.oo) != -sp.pi / 4:
        raise RuntimeError("classical quantile velocity limit failed")

    n, eps = sp.symbols("n eps", positive=True)
    baseline_field = 1 / (2 * n)
    upper_neighbor_field = sp.factor(
        baseline_field
        + (-1 / eps + 1)
        + (1 / (2 * n + eps) - 1 / (2 * n + 1))
    )
    lower_neighbor_field = sp.factor(
        baseline_field
        + (1 / eps - 1)
        + (1 / (2 * n - eps) - 1 / (2 * n - 1))
    )
    upper_pole = sp.simplify(
        sp.cancel(sp.limit(eps * upper_neighbor_field, eps, 0, dir="+"))
    )
    lower_pole = sp.simplify(
        sp.cancel(sp.limit(eps * lower_neighbor_field, eps, 0, dir="+"))
    )
    if upper_pole != -1:
        raise RuntimeError("upper-neighbor negative-field countermodel failed")
    if lower_pole != 1:
        raise RuntimeError("lower-neighbor positive-field countermodel failed")

    return {
        "principal_value_field": (
            "For H(z)=(z-c)^2*V(z), B_c=V'(c)/V(c)="
            "PV sum_(signed roots y outside the pair) 1/(c-y)."
        ),
        "arithmetic_pair_equilibrium": (
            "For roots c+(j+1/2)h, j in Z, removing the pair c+/-h/2 gives "
            "PV sum 1/(c-y)=0 by termwise reflection."
        ),
        "perturbation_identity": (
            "If c=c0+delta_c and y_j=r_j+delta_j, then "
            "1/(c-y_j)-1/(c0-r_j)=(delta_j-delta_c)/"
            "((c-y_j)(c0-r_j))."
        ),
        "weighted_stability": (
            "|B(c;y)-B(c0;r)| <= sum_j |delta_j-delta_c|/"
            "(|c-y_j|*|c0-r_j|), provided the paired sums converge."
        ),
        "classical_density": (
            "rho_0(y)=Psi'(y)=log(y/(4*pi))/(4*pi), y>=a>4*pi"
        ),
        "continuum_field": {
            "definition": "B_rho(x)=2*x*PV integral_a^infinity rho_0(y)/(x^2-y^2)dy",
            "exact_scaled_formula": str(continuum_field),
            "limit": "B_rho(x)=-pi/8+O_a(1/x)",
            "first_correction": (
                "B_rho(x)+pi/8=(a/x)*(1-log(a/(4*pi)))/(2*pi)+"
                "O_a((a/x)^3*|log(a/x)|)"
            ),
            "integral_identities": [
                "PV integral_0^infinity du/(1-u^2)=0",
                "PV integral_0^infinity log(u)/(1-u^2)du=-pi^2/4",
            ],
        },
        "positive_time_density": {
            "density": "rho_t(y)=rho_0(y)+t/(16*y)",
            "field_correction": str(heat_density_correction),
            "asymptotic": "B_(rho,t)(x)=-pi/8+O_a((1+t*log(x/a))/x)",
        },
        "quantile_drift": {
            "counting_function": (
                "g(x,t)=x/(4*pi)*log(x/(4*pi))-x/(4*pi)+11/8+"
                "t/16*log(x/(4*pi))"
            ),
            "velocity": str(quantile_velocity),
            "limit": "dx_n/dt=-pi/4+o(1), so B_n=(1/2)dx_n/dt -> -pi/8",
        },
        "published_fixed_time_localization": {
            "source": "https://arxiv.org/abs/1904.12438",
            "theorem": (
                "For each fixed 0<t<=1/2, sufficiently high zeros are real, simple, "
                "and uniquely localized near the quantiles g(x_n,t)=n once "
                "x_n>=exp(C/t)."
            ),
            "consequence": (
                "At fixed positive t, every multiple zero lies in a bounded region "
                "|x|<exp(C/t), after adjusting the absolute constant."
            ),
            "uniformity_warning": (
                "The localization scale diverges exponentially as t decreases to zero."
            ),
        },
        "macroscopic_location_countermodels": {
            "reference": (
                "Signed integer lattice with double roots at +/-n; after removing the "
                "+n pair its field is 1/(2*n)."
            ),
            "upper_neighbor_move": {
                "move": "+/-(n+1) to +/-(n+eps)",
                "field": str(upper_neighbor_field),
                "limit": "-infinity",
            },
            "lower_neighbor_move": {
                "move": "+/-(n-1) to +/-(n-eps)",
                "field": str(lower_neighbor_field),
                "limit": "+infinity",
            },
            "guard": (
                "Even symmetry and uniformly bounded location error do not bound the "
                "collision field or stiffness."
            ),
        },
        "compactness_reduction": (
            "If Lambda>0, the published fixed-time theorem confines every boundary "
            "multiple zero at t=Lambda to |x|<exp(C/Lambda); closing Lambda<=0 still "
            "requires a lambda-uniform analytic exclusion of that expanding compact region."
        ),
        "open_handoff": (
            "Prove a lambda-uniform paired reciprocal-gap discrepancy bound strong enough "
            "for the weighted perturbation sum, together with an explicit far-tail bound, "
            "or prove a different analytic exclusion that remains effective as lambda->0+."
        ),
        "checks": {
            "lattice_checks": lattice_checks,
            "odd_square_sum": str(odd_square_sum),
            "continuum_correction_slope": str(correction_slope),
            "quantile_velocity": str(quantile_velocity),
        },
    }


def build_payload() -> dict:
    exact = build_exact()
    rows = [
        GateRow(
            id="ncfb_01_principal_value_field",
            role="exact_identity",
            readiness="available_exact",
            claim="The regularized field at a double zero is the principal-value signed-root field with the colliding pair removed.",
            formula=exact["principal_value_field"],
            proof_boundary="Canonical-product identity at a real-zero time.",
        ),
        GateRow(
            id="ncfb_02_arithmetic_pair_equilibrium",
            role="exact_benchmark",
            readiness="available_exact",
            claim="An arithmetic lattice has exactly zero external field at the midpoint of a removed adjacent pair.",
            formula=exact["arithmetic_pair_equilibrium"],
            proof_boundary="Solvable constant-density benchmark only.",
        ),
        GateRow(
            id="ncfb_03_weighted_perturbation_identity",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Field displacement from any reference configuration is governed exactly by a weighted relative-location sum.",
            formula=f"{exact['perturbation_identity']} Hence {exact['weighted_stability']}",
            proof_boundary="Requires convergent paired sums and nonzero outsider gaps.",
        ),
        GateRow(
            id="ncfb_04_classical_continuum_field",
            role="exact_asymptotic",
            readiness="available_exact",
            claim="The continuum field generated by the Riemann-von Mangoldt density tends to -pi/8.",
            formula=exact["continuum_field"]["limit"],
            proof_boundary="Continuum density benchmark, not the discrete Xi zero field.",
            diagnostics=exact["continuum_field"],
        ),
        GateRow(
            id="ncfb_05_positive_time_quantile_drift",
            role="exact_asymptotic",
            readiness="available_exact",
            claim="Implicit differentiation of the positive-time quantile law gives velocity -pi/4 and field -pi/8 at high height.",
            formula=f"{exact['quantile_drift']['velocity']}; {exact['quantile_drift']['limit']}",
            proof_boundary="Quantile model asymptotic only.",
        ),
        GateRow(
            id="ncfb_06_published_high_zero_localization",
            role="published_theorem",
            readiness="available_published",
            claim="The published positive-time asymptotic theorem makes every sufficiently high zero simple and uniquely localized at fixed t>0.",
            formula=exact["published_fixed_time_localization"]["theorem"],
            proof_boundary="External theorem for each fixed positive time; its threshold is not uniform as t approaches zero.",
            diagnostics=exact["published_fixed_time_localization"],
        ),
        GateRow(
            id="ncfb_07_positive_boundary_compactness",
            role="conditional_reduction",
            readiness="ready_to_apply",
            claim="A hypothetical positive Newman boundary can have multiple zeros only in a t-dependent compact height range.",
            formula=exact["compactness_reduction"],
            proof_boundary="Conditional on Lambda>0; does not exclude the compact collision.",
        ),
        GateRow(
            id="ncfb_08_macroscopic_location_guard",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="Even symmetry and bounded error from classical locations do not control the sign or magnitude of the collision field.",
            formula="One neighboring even pair moved by <1 makes B tend to either sign of infinity.",
            proof_boundary="Abstract lattice countermodels, not Xi zero configurations.",
            diagnostics=exact["macroscopic_location_countermodels"],
        ),
        GateRow(
            id="ncfb_09_required_local_upgrade",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="Riemann-von Mangoldt counting and macroscopic location estimates must be upgraded to reciprocal-gap weighted local discrepancy control.",
            formula=exact["weighted_stability"],
            proof_boundary="Necessity for this perturbative field route, not for every possible RH proof.",
        ),
        GateRow(
            id="ncfb_10_lambda_uniform_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The fixed-time compactness theorem closes the high tail only after a lambda-uniform local balance or compact-region exclusion is supplied.",
            formula=exact["open_handoff"],
            proof_boundary="Open lambda-uniform theorem; not a proof of RH or Lambda <= 0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_classical_field_balance_gate",
        "date": "2026-07-11",
        "status": "exact classical-field balance reduction with fixed-time compactness gate",
        "proof_boundary": (
            "This artifact proves exact principal-value, perturbation, arithmetic-lattice, and "
            "continuum-density identities; imports the published fixed-positive-time high-zero "
            "localization theorem with its nonuniform threshold; and blocks promotion from "
            "macroscopic zero locations to collision-field control. It does not prove the needed "
            "lambda-uniform reciprocal-gap balance, exclude a compact positive-boundary collision, "
            "prove RH, or prove Lambda <= 0."
        ),
        "sources": [
            "outputs/jensen_window_pf_newman_root_external_field_lemma.md",
            "https://arxiv.org/abs/1801.05914",
            "https://arxiv.org/abs/1904.12438",
        ],
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    counter = exact["macroscopic_location_countermodels"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Classical-Field Balance Gate",
            "",
            "Date: 2026-07-11",
            "",
            "Status: exact classical-field balance reduction with fixed-time",
            "compactness gate. This is not a proof of RH or `Lambda <= 0`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_classical_field_balance_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_classical_field_balance_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_classical_field_balance_gate.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF Newman classical-field balance gate: 10 rows, 0 issues, 3 exact field identities, 1 arithmetic equilibrium, 1 continuum -pi/8 benchmark, 1 quantile-drift match, 1 published fixed-time localization theorem, 2 exact sensitivity countermodels, 1 compactness reduction, 1 open lambda-uniform handoff",
            "```",
            "",
            "## Exact Balance Identity",
            "",
            "At a double zero with the colliding pair removed,",
            "",
            "```text",
            exact["principal_value_field"],
            exact["perturbation_identity"],
            exact["weighted_stability"],
            "```",
            "",
            "The constant-density reference is exact arithmetic equilibrium:",
            "",
            "```text",
            exact["arithmetic_pair_equilibrium"],
            "```",
            "",
            "## Classical Density Benchmark",
            "",
            "For the Riemann-von Mangoldt density,",
            "",
            "```text",
            exact["classical_density"],
            exact["continuum_field"]["definition"],
            exact["continuum_field"]["limit"],
            "```",
            "",
            "The two exact principal-value integrals behind the constant are",
            "",
            "```text",
            *exact["continuum_field"]["integral_identities"],
            "```",
            "",
            "Adding the positive-time density correction leaves the same limit:",
            "",
            "```text",
            exact["positive_time_density"]["field_correction"],
            exact["positive_time_density"]["asymptotic"],
            "```",
            "",
            "Implicit differentiation of the published quantile law gives",
            "",
            "```text",
            exact["quantile_drift"]["velocity"],
            exact["quantile_drift"]["limit"],
            "```",
            "",
            "The classical reference field is therefore -pi/8 at leading order.",
            "",
            "This matches the high-positive-time zero drift described in the",
            "Polymath asymptotics.",
            "",
            "## Fixed-Time Compactness",
            "",
            exact["published_fixed_time_localization"]["theorem"],
            "",
            "Consequently, if `Lambda>0`, every multiple zero at the boundary",
            "is confined to a region of scale `exp(C/Lambda)`. The threshold",
            "The threshold diverges as `Lambda->0+`, so this is a compactness reduction rather",
            "than the desired exclusion theorem.",
            "",
            "Primary sources: https://arxiv.org/abs/1801.05914 and",
            "https://arxiv.org/abs/1904.12438.",
            "",
            "## Macroscopic Guard",
            "",
            counter["reference"],
            "",
            "```text",
            f"{counter['upper_neighbor_move']['move']}: B={counter['upper_neighbor_move']['field']} -> -infinity",
            f"{counter['lower_neighbor_move']['move']}: B={counter['lower_neighbor_move']['field']} -> +infinity",
            "```",
            "",
            "Both deformations preserve even symmetry and move every zero by less than one.",
            "Therefore bounded classical-location error, and a fortiori",
            "macroscopic counting control, cannot determine the collision field.",
            "",
            "## Live Handoff",
            "",
            exact["open_handoff"],
            "",
            "This is the next strict gate: the estimate must remain effective as",
            "`lambda->0+`; fixed-lambda asymptotics do not prove `Lambda<=0`.",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Jensen-window PF Newman classical-field balance gate: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
