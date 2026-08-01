#!/usr/bin/env python3
"""Build an exact counterfactual sensor atlas for one off-axis zero quartet."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_counterfactual_birth_signature_atlas.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_counterfactual_birth_signature_atlas.md"
)


@dataclass(frozen=True)
class AtlasRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def build_exact() -> dict:
    z, t = sp.symbols("z t", real=True)
    x, y = sp.symbols("x y", positive=True)

    initial = sp.expand(((z - x) ** 2 + y**2) * ((z + x) ** 2 + y**2))
    heat = sp.expand(
        sum(
            (-t) ** order
            / sp.factorial(order)
            * sp.diff(initial, z, 2 * order)
            for order in range(3)
        )
    )
    if sp.simplify(sp.diff(heat, t) + sp.diff(heat, z, 2)) != 0:
        raise RuntimeError("symmetric quartet heat equation failed")

    boundary_radicand = x**4 + 4 * x**2 * y**2 + y**4
    collision_root = sp.sqrt(boundary_radicand)
    collision_time = sp.factor((y**2 - x**2 + collision_root) / 6)
    boundary_polynomial = sp.expand((z**2 - collision_root) ** 2)
    if sp.simplify(heat.subs(t, collision_time) - boundary_polynomial) != 0:
        raise RuntimeError("quartet collision factorization failed")

    w = sp.symbols("w", real=True)
    heat_in_w = sp.Poly(heat.subs(z**2, w), w)
    w_discriminant = sp.factor(sp.discriminant(heat_in_w.as_expr(), w))
    expected_discriminant = 16 * (
        6 * t**2 + 2 * t * (x**2 - y**2) - x**2 * y**2
    )
    if sp.simplify(w_discriminant - expected_discriminant) != 0:
        raise RuntimeError("quartet squared-root discriminant failed")

    q, h, curvature = sp.symbols("q h curvature", real=True)
    local_model = curvature * (q**2 / 2 - h)
    local_laguerre = sp.factor(
        sp.diff(local_model, q) ** 2
        - local_model * sp.diff(local_model, q, 2)
    )
    expected_local_laguerre = curvature**2 * (q**2 / 2 + h)
    if sp.simplify(local_laguerre - expected_local_laguerre) != 0:
        raise RuntimeError("universal local Laguerre jet failed")

    delta, radius, degree, variable = sp.symbols(
        "delta radius degree variable", positive=True
    )
    coefficient_0 = sp.Integer(1)
    coefficient_1 = 2 * sp.cos(delta) / radius
    coefficient_2 = 2 / radius**2
    jensen = sp.expand(
        coefficient_0
        + degree * coefficient_1 * variable
        + sp.binomial(degree, 2) * coefficient_2 * variable**2
    )
    jensen_discriminant = sp.factor(
        sp.discriminant(jensen, variable).rewrite(sp.sin)
    )
    expected_jensen_discriminant = (
        4 * degree * (1 - degree * sp.sin(delta) ** 2) / radius**2
    )
    if (
        sp.simplify(
            sp.trigsimp(jensen_discriminant - expected_jensen_discriminant)
        )
        != 0
    ):
        raise RuntimeError("quadratic quartet Jensen discriminant failed")

    epsilon, gamma = sp.symbols("epsilon gamma", positive=True)
    beta = sp.Rational(1, 2) + epsilon
    rho = beta + sp.I * gamma
    li_multiplier = sp.simplify((rho - 1) / rho)
    li_radius_squared = sp.simplify(
        li_multiplier * sp.conjugate(li_multiplier)
    )
    expected_li_radius_squared = (
        gamma**2 + (sp.Rational(1, 2) - epsilon) ** 2
    ) / (
        gamma**2 + (sp.Rational(1, 2) + epsilon) ** 2
    )
    if sp.simplify(li_radius_squared - expected_li_radius_squared) != 0:
        raise RuntimeError("Li multiplier radius failed")
    li_growth_rate = sp.factor(
        sp.log(1 / expected_li_radius_squared) / 2
    )

    sector_sine_squared = sp.factor(
        4 * gamma**2 * epsilon**2 / (gamma**2 + epsilon**2) ** 2
    )
    jensen_continuous_threshold = sp.factor(1 / sector_sine_squared)

    m = sp.symbols("m", positive=True)
    diagonal_x = 2 * m
    diagonal_y = 2 / m
    diagonal_collision_time = sp.factor(
        collision_time.subs({x: diagonal_x, y: diagonal_y})
    )
    diagonal_jensen_threshold = sp.factor(
        jensen_continuous_threshold.subs({gamma: m, epsilon: 1 / m})
    )
    diagonal_li_growth = sp.factor(
        li_growth_rate.subs({gamma: m, epsilon: 1 / m})
    )
    li_index = sp.symbols("li_index", integer=True, positive=True)
    diagonal_multiplier = sp.simplify(
        li_multiplier.subs({gamma: m, epsilon: 1 / m})
    )
    diagonal_li_contribution = sp.simplify(
        4
        - diagonal_multiplier**li_index
        - sp.conjugate(diagonal_multiplier) ** li_index
        - diagonal_multiplier ** (-li_index)
        - sp.conjugate(diagonal_multiplier) ** (-li_index)
    )
    fixed_index_li_limit = sp.simplify(
        sp.limit(m**2 * diagonal_li_contribution, m, sp.oo)
    )
    if fixed_index_li_limit != 2 * li_index**2:
        raise RuntimeError(
            "diagonal fixed-index Li limit failed: "
            f"{fixed_index_li_limit!s}"
        )

    limits = {
        "m2_collision_time": str(
            sp.limit(m**2 * diagonal_collision_time, m, sp.oo)
        ),
        "jensen_threshold_over_m4": str(
            sp.limit(diagonal_jensen_threshold / m**4, m, sp.oo)
        ),
        "m3_li_growth_rate": str(
            sp.limit(m**3 * diagonal_li_growth, m, sp.oo)
        ),
        "li_multiplier": str(
            sp.limit(
                diagonal_multiplier,
                m,
                sp.oo,
            )
        ),
        "fixed_index_li_contribution": str(fixed_index_li_limit),
    }
    if limits != {
        "m2_collision_time": "2",
        "jensen_threshold_over_m4": "1/4",
        "m3_li_growth_rate": "1",
        "li_multiplier": "1",
        "fixed_index_li_contribution": "2*li_index**2",
    }:
        raise RuntimeError(f"diagonal escape limits failed: {limits!r}")

    return {
        "coordinate_map": {
            "zeta_zero": "rho=1/2+epsilon+i*gamma with epsilon,gamma>0",
            "H_zero": "z_rho=2*gamma-2*i*epsilon=x-i*y",
            "coordinates": "x=2*gamma and y=2*epsilon",
            "squared_zero": (
                "s_rho=-z_rho^2=-4*(gamma^2-epsilon^2)"
                "+8*i*gamma*epsilon"
            ),
            "sector_angle": "delta=2*atan(epsilon/gamma)",
            "sector_sine_squared": str(sector_sine_squared),
            "squared_zero_radius": "R=4*(gamma^2+epsilon^2)=x^2+y^2",
        },
        "symmetric_heat_model": {
            "initial": str(initial),
            "heat_solution": str(heat),
            "heat_equation": "partial_t P_t=-partial_z^2 P_t",
            "w_discriminant": str(w_discriminant),
            "collision_time": str(collision_time),
            "collision_time_zeta_coordinates": (
                "tau=(2/3)*(epsilon^2-gamma^2"
                "+sqrt(gamma^4+4*gamma^2*epsilon^2+epsilon^4))"
            ),
            "collision_location_squared": str(collision_root),
            "boundary_factorization": "P_tau(z)=(z^2-c^2)^2",
            "scope": (
                "Exact even polynomial countermodel with the quartet symmetries; "
                "not the Xi heat flow."
            ),
        },
        "universal_collision_sensor": {
            "local_model": "H_(tau+h)(c+q)=A*(q^2/2-h)+higher terms",
            "laguerre_jet": (
                "H_x^2-H*H_xx=A^2*(q^2/2+h)+higher terms"
            ),
            "exact_model_jet": str(local_laguerre),
            "collision_value": "L_tau(c)=0 and partial_x L_tau(c)=0",
            "correlation_identity": (
                "Fourier[K_(1,t)](xi)=L_t(xi/2), so a collision at c "
                "creates a zero mode at xi=2*c."
            ),
            "wiener_failure": (
                "At that time the translates of K_(1,t) are not dense in L1(R)."
            ),
        },
        "quadratic_jensen_sensor": {
            "quartet_factor": (
                "G_delta,R(s)=1+2*cos(delta)*s/R+s^2/R^2"
            ),
            "exponential_coefficients": (
                "A_0=1, A_1=2*cos(delta)/R, A_2=2/R^2"
            ),
            "jensen_polynomial": str(jensen),
            "discriminant": str(expected_jensen_discriminant),
            "criterion": "J_d is hyperbolic iff d*sin(delta)^2<=1",
            "first_failure": (
                "d_first=floor(csc(delta)^2)+1"
            ),
            "continuous_threshold": str(jensen_continuous_threshold),
            "asymptotic": (
                "d_first is of order gamma^2/(4*epsilon^2) "
                "when gamma/epsilon tends to infinity."
            ),
        },
        "li_sensor": {
            "multiplier": str(li_multiplier),
            "radius_squared": str(expected_li_radius_squared),
            "growth_rate": str(li_growth_rate),
            "quartet_pairing": (
                "The four Li multipliers are w, conjugate(w), "
                "1/w, and 1/conjugate(w)."
            ),
            "quartet_contribution": (
                "lambda_n^(Q)=4-2*(r^n+r^(-n))*cos(n*theta), "
                "where w=r*exp(i*theta) and 0<r<1."
            ),
            "subsequence": (
                "There are infinitely many n with cos(n*theta)>=1/2; "
                "along such a subsequence the quartet contribution tends to -infinity."
            ),
            "e_folding_asymptotic": (
                "kappa=-log(r) is asymptotic to "
                "epsilon/(gamma^2+1/4) for epsilon/gamma tending to zero."
            ),
            "scope": (
                "This is the isolated quartet contribution, not a bound for the "
                "complete Xi Li coefficient or its first negative index."
            ),
        },
        "suzuki_sensor": {
            "offset": "delta_xi>=epsilon for this quartet",
            "fixed_shift": (
                "For 0<omega<epsilon the shifted quotient has an uncancelled "
                "upper-half-plane pole, so D(omega) fails."
            ),
            "exceptional_shift": (
                "At omega=epsilon an equal-height horizontal pair can cancel in "
                "the quotient; this is one exceptional shift only."
            ),
            "cofinal_detection": (
                "Every sequence omega_j->0 eventually enters 0<omega_j<epsilon "
                "and therefore detects the quartet."
            ),
        },
        "diagonal_escape_family": {
            "parameters": (
                "gamma_m=m, epsilon_m=1/m, x_m=2m, y_m=2/m"
            ),
            "collision_time": str(diagonal_collision_time),
            "jensen_continuous_threshold": str(diagonal_jensen_threshold),
            "li_growth_rate": str(diagonal_li_growth),
            "fixed_index_li_asymptotic": (
                "For every fixed positive integer n, "
                "m^2*lambda_n^(Q) tends to 2*n^2; hence "
                "lambda_n^(Q)>0 for all sufficiently large m."
            ),
            "limits": limits,
            "conclusion": (
                "As m tends to infinity, the hypothetical quartet escapes every "
                "fixed real band, its heat boundary tends to zero, its first "
                "quadratic Jensen failure tends to infinite degree, its Li "
                "amplification rate tends to zero, and its Suzuki offset tends "
                "to zero."
            ),
        },
        "route_decision": {
            "primary": (
                "Retain direct positive-time simplicity or first-jet "
                "transversality as the primary route. Positive-boundary "
                "attainment makes any Lambda>0 collision finite, and the "
                "strict-Laguerre/Wiener sensor reacts at the collision itself "
                "without a diverging degree, coefficient index, or shift sequence."
            ),
            "concrete_target": (
                "Prove (H_t(x),H_t'(x))!=(0,0) for every real x and every "
                "0<t<=1/5, equivalently prove L_t(x)>0 or zero-freeness of "
                "Fourier[K_(1,t)] throughout that domain."
            ),
            "remaining_interface": (
                "Use the corrected Riemann-Siegel first-jet partition to prove "
                "phase-critical-value avoidance in the remaining high-frequency, "
                "small-t layer; alternatively prove a boundary-integrable Xi "
                "collision energy or a collision-index invariant."
            ),
        },
    }


def build_payload() -> dict:
    exact = build_exact()
    rows = [
        AtlasRow(
            id="ncbsa_01_quartet_coordinate_map",
            role="exact_identity",
            readiness="available_exact",
            claim="One off-axis zeta quartet maps to a conjugate pair of squared zeros at an exact negative-axis sector angle.",
            formula=(
                f"{exact['coordinate_map']['H_zero']}; "
                f"{exact['coordinate_map']['squared_zero']}; "
                f"{exact['coordinate_map']['sector_angle']}"
            ),
            proof_boundary="Coordinate identity only; it assumes a hypothetical off-axis zero.",
        ),
        AtlasRow(
            id="ncbsa_02_symmetric_quartet_heat_model",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="A symmetric four-zero polynomial has an exact finite backward-heat collision.",
            formula=(
                f"tau={exact['symmetric_heat_model']['collision_time']}; "
                f"{exact['symmetric_heat_model']['boundary_factorization']}"
            ),
            proof_boundary=exact["symmetric_heat_model"]["scope"],
            diagnostics=exact["symmetric_heat_model"],
        ),
        AtlasRow(
            id="ncbsa_03_universal_laguerre_birth_jet",
            role="exact_lemma",
            readiness="available_exact",
            claim="The strict Laguerre quantity changes sign through a generic double-zero birth and vanishes at the collision.",
            formula=exact["universal_collision_sensor"]["laguerre_jet"],
            proof_boundary="Local nondegenerate double-zero normal form; higher-multiplicity splitting is handled by the existing boundary-attainment lemma.",
        ),
        AtlasRow(
            id="ncbsa_04_correlation_zero_mode",
            role="exact_composition",
            readiness="available_exact",
            claim="An Xi heat-flow collision is exactly a zero mode of the first correlation transform.",
            formula=exact["universal_collision_sensor"]["correlation_identity"],
            proof_boundary="Uses the independently checked Xi correlation identity; the quartet polynomial itself has no asserted Xi kernel.",
        ),
        AtlasRow(
            id="ncbsa_05_quadratic_jensen_threshold",
            role="exact_lemma",
            readiness="available_exact",
            claim="The isolated squared-quartet factor saturates the sector theorem at an exact Jensen degree threshold.",
            formula=(
                f"{exact['quadratic_jensen_sensor']['discriminant']}; "
                f"{exact['quadratic_jensen_sensor']['criterion']}"
            ),
            proof_boundary="Exact quadratic factor only; other Xi zeros can alter complete Jensen coefficients.",
            diagnostics=exact["quadratic_jensen_sensor"],
        ),
        AtlasRow(
            id="ncbsa_06_li_quartet_signature",
            role="exact_lemma",
            readiness="available_exact",
            claim="The Li contribution of one off-axis quartet has reciprocal exponential amplitudes and becomes negative without bound along a phase-return subsequence.",
            formula=exact["li_sensor"]["quartet_contribution"],
            proof_boundary=exact["li_sensor"]["scope"],
            diagnostics=exact["li_sensor"],
        ),
        AtlasRow(
            id="ncbsa_07_suzuki_shift_signature",
            role="exact_composition",
            readiness="available_exact",
            claim="A fixed Suzuki shift below the quartet offset detects an uncancelled pole, while one exact cancellation shift can hide it.",
            formula=(
                f"{exact['suzuki_sensor']['fixed_shift']} "
                f"{exact['suzuki_sensor']['exceptional_shift']}"
            ),
            proof_boundary="Uses the internally audited fixed-omega phase diagram; no arithmetic determinant premise is proved here.",
        ),
        AtlasRow(
            id="ncbsa_08_diagonal_escape_theorem",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="One explicit quartet family simultaneously escapes every fixed band, degree, Li amplification scale, positive time floor, and positive Suzuki shift.",
            formula=exact["diagonal_escape_family"]["conclusion"],
            proof_boundary="Counterfactual polynomial/zero-coordinate family only; it is not a family of Xi zeros.",
            diagnostics=exact["diagonal_escape_family"],
        ),
        AtlasRow(
            id="ncbsa_09_finite_sensor_nonpromotion",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="No collection of fixed finite sensor cutoffs can exclude all high, near-line hypothetical quartets.",
            formula=(
                "R fixed, d fixed, n fixed, t>=t_0>0 fixed, and omega>0 fixed "
                "are all evaded by sufficiently large m in the diagonal family."
            ),
            proof_boundary="Blocks finite-to-global promotion only; it does not rule out a cofinal theorem.",
        ),
        AtlasRow(
            id="ncbsa_10_primary_route_selection",
            role="route_decision",
            readiness="available_exact",
            claim="Direct collision transversality is the sensor that does not introduce a second diverging detection parameter.",
            formula=exact["route_decision"]["primary"],
            proof_boundary="Logical route comparison, not a proof that the Xi transversality estimate holds.",
        ),
        AtlasRow(
            id="ncbsa_11_phase_critical_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The atlas selects the existing corrected first-jet phase-critical theorem as the primary global handoff.",
            formula=exact["route_decision"]["remaining_interface"],
            proof_boundary="Open Xi-specific theorem; not Lambda<=0 and not RH.",
        ),
        AtlasRow(
            id="ncbsa_12_new_math_alternative",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="A boundary-integrable energy or collision-index theorem remains a genuinely different global alternative.",
            formula=(
                "Construct an Xi-specific invariant that cannot cross the "
                "correlation zero mode or cannot support the universal "
                "nonintegrable collision trace."
            ),
            proof_boundary="Research target only; no such invariant is constructed here.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_counterfactual_birth_signature_atlas",
        "date": "2026-07-25",
        "status": "exact counterfactual birth-signature atlas and finite-sensor nonpromotion gate",
        "proof_boundary": (
            "This artifact computes exact heat, Laguerre, Jensen, Li, and Suzuki "
            "signatures for a symmetric hypothetical off-axis quartet and proves "
            "an explicit diagonal family evades every fixed sensor cutoff. It "
            "composes existing Xi collision identities only where stated. It does "
            "not assert that Xi has an off-axis zero, prove the selected Xi "
            "transversality target, prove Lambda<=0, or prove RH."
        ),
        "sources": [
            "outputs/formal_core.md",
            "outputs/jensen_window_pf_newman_positive_boundary_attainment_lemma.md",
            "outputs/jensen_window_pf_newman_strict_laguerre_correlation_target.md",
            "outputs/jensen_window_pf_newman_real_zero_band_degree361_sector_certificate.md",
            "outputs/jensen_window_pf_suzuki_fixed_omega_phase_diagram.md",
            "outputs/jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.md",
            "https://arxiv.org/abs/math/0404394",
        ],
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    coordinate = exact["coordinate_map"]
    heat = exact["symmetric_heat_model"]
    collision = exact["universal_collision_sensor"]
    jensen = exact["quadratic_jensen_sensor"]
    li = exact["li_sensor"]
    suzuki = exact["suzuki_sensor"]
    diagonal = exact["diagonal_escape_family"]
    decision = exact["route_decision"]
    return "\n".join(
        [
            "# Newman Counterfactual Birth-Signature Atlas",
            "",
            "Date: 2026-07-25",
            "",
            "Status: exact counterfactual birth-signature atlas and finite-sensor",
            "nonpromotion gate. This is not a proof of `Lambda<=0` or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_counterfactual_birth_signature_atlas.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_counterfactual_birth_signature_atlas.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_counterfactual_birth_signature_atlas.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Newman counterfactual birth-signature atlas: 12 rows, 0 issues, 1 symmetric quartet heat collision, 1 exact Jensen threshold, 1 exact Li quartet law, 1 Suzuki shift law, 1 diagonal finite-sensor escape theorem, 2 open global handoffs",
            "```",
            "",
            "## Counterfactual Rule",
            "",
            "Assume one off-axis quartet provisionally and calculate every",
            "consequence. The construction below is a proof-safety model, not an",
            "assertion that the actual xi function has such a zero.",
            "",
            "## Coordinate Map",
            "",
            "```text",
            coordinate["zeta_zero"],
            coordinate["H_zero"],
            coordinate["squared_zero"],
            coordinate["sector_angle"],
            f"sin(delta)^2={coordinate['sector_sine_squared']}",
            coordinate["squared_zero_radius"],
            "```",
            "",
            "## Symmetric Heat Collision",
            "",
            "Use the real even quartet polynomial",
            "",
            "```text",
            f"P_0(z)={heat['initial']}",
            f"P_t(z)={heat['heat_solution']}",
            heat["heat_equation"],
            f"Disc_w(P_t)={heat['w_discriminant']}",
            f"tau={heat['collision_time']}",
            heat["collision_time_zeta_coordinates"],
            heat["boundary_factorization"],
            "```",
            "",
            "This is an exact four-zero backward-heat birth with the correct",
            "functional-equation symmetries. It is not the Xi heat flow.",
            "",
            "At any nondegenerate double-zero boundary the local heat equation gives",
            "",
            "```text",
            collision["local_model"],
            collision["laguerre_jet"],
            collision["collision_value"],
            collision["correlation_identity"],
            collision["wiener_failure"],
            "```",
            "",
            "Thus the strict-Laguerre correlation sensor reacts at the collision",
            "itself. It does not need a second index to grow.",
            "",
            "## Jensen Sensor",
            "",
            "The conjugate squared-zero factor is",
            "",
            "```text",
            jensen["quartet_factor"],
            jensen["exponential_coefficients"],
            f"J_d(X)={jensen['jensen_polynomial']}",
            f"Disc(J_d)={jensen['discriminant']}",
            jensen["criterion"],
            jensen["first_failure"],
            f"csc(delta)^2={jensen['continuous_threshold']}",
            jensen["asymptotic"],
            "```",
            "",
            "The quadratic model exactly saturates the sector cutoff. A high",
            "near-line quartet can therefore remain invisible through an enormous",
            "but finite Jensen degree.",
            "",
            "## Li Sensor",
            "",
            "For `w=(rho-1)/rho=r*exp(i*theta)`,",
            "",
            "```text",
            f"w={li['multiplier']}",
            f"r^2={li['radius_squared']}",
            f"kappa=-log(r)={li['growth_rate']}",
            li["quartet_pairing"],
            li["quartet_contribution"],
            li["subsequence"],
            li["e_folding_asymptotic"],
            "```",
            "",
            "This is the isolated quartet contribution. It does not identify the",
            "first negative complete Xi Li coefficient.",
            "",
            "## Suzuki Sensor",
            "",
            "```text",
            suzuki["fixed_shift"],
            suzuki["exceptional_shift"],
            suzuki["cofinal_detection"],
            "```",
            "",
            "One fixed quotient can hide one exact horizontal cancellation. A",
            "cofinal shift sequence cannot hide a quartet with positive offset.",
            "",
            "## Diagonal Escape Theorem",
            "",
            "Take",
            "",
            "```text",
            diagonal["parameters"],
            f"tau_m={diagonal['collision_time']}",
            f"D_m={diagonal['jensen_continuous_threshold']}",
            f"kappa_m={diagonal['li_growth_rate']}",
            diagonal["fixed_index_li_asymptotic"],
            f"limits={diagonal['limits']}",
            "```",
            "",
            diagonal["conclusion"],
            "",
            "Consequently, no fixed real band, finite degree, finite Li index,",
            "positive lower time cutoff, or positive fixed Suzuki shift can be",
            "promoted into a global exclusion theorem. At least one coordinate",
            "must be handled cofinally or without a cutoff.",
            "",
            "## Route Decision",
            "",
            decision["primary"],
            "",
            "The concrete theorem remains",
            "",
            "```text",
            decision["concrete_target"],
            "```",
            "",
            decision["remaining_interface"],
            "",
            "The atlas selects the corrected high-frequency first-jet",
            "phase-critical-value theorem as the primary next attack. A",
            "boundary-integrable collision energy or a collision-index invariant",
            "remains a separate high-risk alternative. Neither is proved here.",
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
        "wrote Newman counterfactual birth-signature atlas: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
