#!/usr/bin/env python3
"""Reduce corrected critical transversality to a thin-collar Wronskian target."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import mpmath as mp
import sympy as sp

import jensen_window_pf_newman_polymath15_critical_scaled_coercivity_scout as base


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction.md"
)
COARSE_DPS = 55
FINE_DPS = 75


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
    x, y, u, v, ell = sp.symbols("x y u v ell", real=True, nonzero=True)
    amplitude, radial, phase_velocity, theta = sp.symbols(
        "amplitude radial phase_velocity theta", real=True, nonzero=True
    )
    wronskian = sp.expand(sp.im((u + sp.I * v) * (x - sp.I * y)))
    if sp.simplify(wronskian - (v * x - u * y)) != 0:
        raise RuntimeError("complex Wronskian identity failed")

    polar_value = amplitude * sp.cos(theta)
    polar_first = amplitude * (
        radial * sp.cos(theta) - phase_velocity * sp.sin(theta)
    )
    polar_transversality = sp.expand(
        4 * (polar_value**2 + (polar_first / ell) ** 2)
    )
    expected = 4 * amplitude**2 * (
        sp.cos(theta) ** 2
        + (
            radial * sp.cos(theta) - phase_velocity * sp.sin(theta)
        )
        ** 2
        / ell**2
    )
    if sp.simplify(polar_transversality - expected) != 0:
        raise RuntimeError("polar first-jet identity failed")

    at_crossing = sp.simplify(
        expected.subs(sp.cos(theta), 0).subs(sp.sin(theta) ** 2, 1)
    )
    expected_crossing = 4 * amplitude**2 * phase_velocity**2 / ell**2
    if sp.simplify(at_crossing - expected_crossing) != 0:
        raise RuntimeError("crossing Wronskian identity failed")

    return {
        "complex_split": (
            "Choose the sharp endpoint half q_(N,t) so Q_(N,t)=2*Re(q_(N,t)) "
            "on the real axis, and set E_(N,t)=exp(i*beta_t)D_(N,t)-q_(N,t); "
            "then J_(N,t)=2*Re(E_(N,t)) and J_(N,t)'=2*Re(E_(N,t)')"
        ),
        "cartesian": (
            "For E=X+iY and E'=U+iV, T_L[J]=4*(X^2+(U/L)^2)"
        ),
        "wronskian": (
            "W_E=Im(E'*conj(E))=V*X-U*Y"
        ),
        "polar": (
            "If E=a*exp(i*theta), u=(log(a))', and v=theta', then "
            "T_L[J]=4*a^2*(cos(theta)^2+"
            "(u*cos(theta)-v*sin(theta))^2/L^2) and W_E=a^2*v"
        ),
        "crossing": (
            "At Re(E)=0 with E!=0, T_L[J]=4*W_E^2/(|E|^2*L^2); "
            "J=J'=0 if and only if theta'=0"
        ),
        "collision_wronskian_cap": (
            "If Z=J+r=Z'=0, |r|<=epsilon_0, and |r'|<=L*epsilon_1, "
            "then 2*|W_E|<=epsilon_0*|E'|+L*epsilon_1*|E|"
        ),
        "thin_collar_target": (
            "For every critical point, either |Re(E)|>epsilon_0/2 or "
            "2*|W_E|>epsilon_0*|E'|+L*epsilon_1*|E|"
        ),
        "target_domains": (
            "The explicit formulas remain valid on the original domain L>=50, "
            "0<tL<=25. After the oscillatory-zeta theorem, the live asymptotic "
            "high-frequency layer is 0<tL<=c_*+o(1), with "
            "c_*=4911678521/1933561194; uncovered bounded L remains a separate "
            "compact-certificate obligation"
        ),
        "explicit_target": (
            "With epsilon_0=2500*exp(-3L/4) and epsilon_1=5000*exp(-3L/4), "
            "either |Re(E)|>1250*exp(-3L/4) or "
            "2*|W_E|>exp(-3L/4)*(2500*|E'|+5000*L*|E|)"
        ),
        "critical_value_form": (
            "Where E!=0, collision means theta is pi/2 modulo pi and theta'=0; "
            "the open arithmetic theorem is quantitative avoidance of those phase critical levels"
        ),
        "generator_bridge": (
            "At J=0 with E=iY and Y!=0, W_E=-Y*J'/2; along a moving simple "
            "zero x_*(t), x_*'=-J_t/J'=J_t*Y/(2*W_E)"
        ),
        "pair_wronskian": (
            "If J(x)=(x-alpha)*(x-beta)*G(x), then "
            "|W_E(alpha)|=|E(alpha)|*(beta-alpha)*|G(alpha)|/2; "
            "the Wronskian margin is linearly close-pair sensitive"
        ),
    }


def corrected_complex_main(x: mp.mpf, time: mp.mpf, cutoff: int) -> mp.mpc:
    s = (1 - mp.j * x) / 2
    alpha = base.alpha(s)
    shifted_s = s + time * alpha / 2
    finite_sum = mp.fsum(
        mp.exp(
            time * mp.log(n) ** 2 / 4 - shifted_s * mp.log(n)
        )
        for n in range(1, cutoff + 1)
    )
    normalizer = base.m_time(s, time)
    saddle_time = x / 2 + mp.pi * time / 8
    saddle = mp.sqrt(saddle_time / (2 * mp.pi))
    p = 1 - 2 * (saddle - cutoff)
    u_rs = mp.exp(
        -mp.j
        * (
            saddle_time / 2 * mp.log(saddle_time / (2 * mp.pi))
            - saddle_time / 2
            - mp.pi / 8
        )
    )
    endpoint_half = (
        (-1) ** cutoff
        * mp.exp(time * mp.pi**2 / 64)
        * base.m_zero(mp.j * saddle_time)
        * base.c_zero(p)
        * u_rs
        * mp.exp(mp.j * mp.pi / 8)
        / abs(normalizer)
    )
    return normalizer / abs(normalizer) * finite_sum - endpoint_half


def phase_velocity(x: mp.mpf, time: mp.mpf, cutoff: int) -> tuple[mp.mpc, mp.mpc, mp.mpf]:
    function = lambda value: corrected_complex_main(value, time, cutoff)
    value = function(x)
    first = mp.diff(function, x)
    velocity = mp.im(first * mp.conj(value)) / abs(value) ** 2
    return value, first, velocity


def diagnostics(dps: int) -> dict:
    mp.mp.dps = dps
    time = mp.mpf(0)
    sign_rows = []
    for x_text, cutoff in (("200", 3), ("520", 6)):
        x = mp.mpf(x_text)
        value, _, velocity = phase_velocity(x, time, cutoff)
        sign_rows.append(
            {
                "x": x_text,
                "cutoff": cutoff,
                "corrected_value": mp.nstr(2 * mp.re(value), 35, strip_zeros=False),
                "phase_velocity": mp.nstr(velocity, 35, strip_zeros=False),
            }
        )

    root_rows = []
    for bracket in (("517.1", "517.3"), ("519.6", "519.9")):
        function = lambda value: 2 * mp.re(corrected_complex_main(value, time, 6))
        root = mp.findroot(function, tuple(mp.mpf(value) for value in bracket))
        value, _, velocity = phase_velocity(root, time, 6)
        root_rows.append(
            {
                "bracket": list(bracket),
                "root": mp.nstr(root, 40, strip_zeros=False),
                "abs_corrected_value_residual": mp.nstr(
                    abs(2 * mp.re(value)), 25, strip_zeros=False
                ),
                "abs_E": mp.nstr(abs(value), 35, strip_zeros=False),
                "phase_velocity": mp.nstr(velocity, 35, strip_zeros=False),
            }
        )
    if not (
        mp.mpf(sign_rows[0]["phase_velocity"]) < 0
        and mp.mpf(sign_rows[1]["phase_velocity"]) > 0
    ):
        raise RuntimeError("corrected phase sign-change diagnostic failed")
    if not (
        mp.mpf(root_rows[0]["phase_velocity"]) < 0
        and mp.mpf(root_rows[1]["phase_velocity"]) > 0
    ):
        raise RuntimeError("corrected crossing phase signs failed")
    if any(
        mp.mpf(row["abs_corrected_value_residual"]) >= mp.mpf("1e-45")
        for row in root_rows
    ):
        raise RuntimeError("corrected crossing root residual failed")

    lehmer_rows = []
    for bracket in (("14010.12", "14010.13"), ("14010.19", "14010.21")):
        function = lambda value: 2 * mp.re(corrected_complex_main(value, time, 33))
        root = mp.findroot(function, tuple(mp.mpf(value) for value in bracket))
        value, first, velocity = phase_velocity(root, time, 33)
        lehmer_rows.append(
            {
                "bracket": list(bracket),
                "root": mp.nstr(root, 40, strip_zeros=False),
                "abs_corrected_value_residual": mp.nstr(
                    abs(2 * mp.re(value)), 25, strip_zeros=False
                ),
                "abs_E": mp.nstr(abs(value), 35, strip_zeros=False),
                "abs_W": mp.nstr(
                    abs(mp.im(first * mp.conj(value))), 35, strip_zeros=False
                ),
                "phase_velocity": mp.nstr(velocity, 35, strip_zeros=False),
            }
        )
    lehmer_gap = mp.mpf(lehmer_rows[1]["root"]) - mp.mpf(
        lehmer_rows[0]["root"]
    )
    if not lehmer_gap > 0:
        raise RuntimeError("corrected Lehmer pair ordering failed")
    if any(
        mp.mpf(row["abs_corrected_value_residual"]) >= mp.mpf("1e-45")
        for row in lehmer_rows
    ):
        raise RuntimeError("corrected Lehmer crossing residual failed")
    return {
        "role": "finite_route_shaping_diagnostics_only",
        "proof_boundary": (
            "Moderate-frequency t=0 points show that neither a global nor a "
            "crossing-only one-sided phase sign follows formally. They do not "
            "establish the L>=50 thin-collar target."
        ),
        "dps": dps,
        "sign_rows": sign_rows,
        "root_rows": root_rows,
        "lehmer_pair": {
            "roots": lehmer_rows,
            "gap": mp.nstr(lehmer_gap, 35, strip_zeros=False),
            "min_abs_phase_velocity": mp.nstr(
                min(abs(mp.mpf(row["phase_velocity"])) for row in lehmer_rows),
                30,
                strip_zeros=False,
            ),
        },
    }


def compare_diagnostics(coarse: dict, fine: dict) -> dict:
    root_delta = mp.mpf(0)
    velocity_delta = mp.mpf(0)
    for left, right in zip(coarse["root_rows"], fine["root_rows"], strict=True):
        root_delta = max(
            root_delta,
            abs(mp.mpf(left["root"]) - mp.mpf(right["root"])),
        )
        velocity_delta = max(
            velocity_delta,
            abs(
                mp.mpf(left["phase_velocity"])
                - mp.mpf(right["phase_velocity"])
            ),
        )
    for left, right in zip(
        coarse["lehmer_pair"]["roots"],
        fine["lehmer_pair"]["roots"],
        strict=True,
    ):
        root_delta = max(
            root_delta,
            abs(mp.mpf(left["root"]) - mp.mpf(right["root"])),
        )
        velocity_delta = max(
            velocity_delta,
            abs(
                mp.mpf(left["phase_velocity"])
                - mp.mpf(right["phase_velocity"])
            ),
        )
    if root_delta >= mp.mpf("1e-35") or velocity_delta >= mp.mpf("1e-30"):
        raise RuntimeError("Wronskian diagnostics are precision-unstable")
    return {
        "max_abs_root_delta": mp.nstr(root_delta, 25, strip_zeros=False),
        "max_abs_phase_velocity_delta": mp.nstr(
            velocity_delta, 25, strip_zeros=False
        ),
    }


def build_artifact() -> dict:
    exact = build_exact()
    coarse = diagnostics(COARSE_DPS)
    fine = diagnostics(FINE_DPS)
    convergence = compare_diagnostics(coarse, fine)
    rows = [
        GateRow(
            id="np15cwpr_01_complex_split",
            role="exact_definition",
            readiness="ready_to_apply",
            claim="The corrected real main is the real part of one explicit complex main.",
            formula=exact["complex_split"],
            proof_boundary="Real axis with one fixed analytic cutoff lift.",
        ),
        GateRow(
            id="np15cwpr_02_cartesian_jet",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The corrected C1 target is the Euclidean first jet of the real part of E.",
            formula=exact["cartesian"],
            proof_boundary="No division by E.",
        ),
        GateRow(
            id="np15cwpr_03_polar_jet",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Away from zeros of E, the first jet separates radial and phase motion exactly.",
            formula=exact["polar"],
            proof_boundary="Only the polar rewrite assumes E is nonzero.",
        ),
        GateRow(
            id="np15cwpr_04_crossing_wronskian",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="At a corrected crossing, simplicity is exactly nonvanishing phase velocity.",
            formula=f"{exact['wronskian']}; {exact['crossing']}",
            proof_boundary="E nonzero at the crossing.",
        ),
        GateRow(
            id="np15cwpr_05_collision_cap",
            role="exact_exclusion_lemma",
            readiness="ready_to_apply",
            claim="A true double zero forces an explicit upper bound on the complex Wronskian.",
            formula=exact["collision_wronskian_cap"],
            proof_boundary="Uses only the certified normalized C1 remainder caps.",
        ),
        GateRow(
            id="np15cwpr_06_thin_collar",
            role="strategy_reduction",
            readiness="conditional_ready",
            claim="It is sufficient to prove Wronskian separation only in the exponentially thin corrected-value collar.",
            formula=exact["thin_collar_target"],
            proof_boundary="Sufficient, not asserted as established.",
        ),
        GateRow(
            id="np15cwpr_07_explicit_target",
            role="open_theorem_target",
            readiness="not_ready_to_apply",
            claim="The global remainder constants give a concrete crossing-localized arithmetic target.",
            formula=exact["explicit_target"],
            proof_boundary=(
                "Open on the original L>=50, 0<tL<=25 domain. Asymptotically "
                "the live high-frequency layer is 0<tL<=c_*+o(1); bounded-L "
                "residuals remain separate compact obligations."
            ),
        ),
        GateRow(
            id="np15cwpr_08_phase_critical_values",
            role="dynamical_reformulation",
            readiness="conditional_ready",
            claim="The surviving phase-critical obstruction is exactly the collapsing denominator of the zero generator and is linearly close-pair sensitive.",
            formula=(
                f"{exact['critical_value_form']}; {exact['generator_bridge']}; "
                f"{exact['pair_wronskian']}"
            ),
            proof_boundary="Quantitative control is still missing.",
        ),
        GateRow(
            id="np15cwpr_09_sign_diagnostic",
            role="route_rejection_diagnostic",
            readiness="diagnostic_only",
            claim="Finite corrected examples have both phase signs, including at corrected crossings.",
            formula="one-sided phase sign is not an algebraic consequence of the corrected representation",
            proof_boundary="Moderate-frequency finite evidence only.",
            diagnostics={
                "convergence": convergence,
                "sign_rows": fine["sign_rows"],
                "root_rows": fine["root_rows"],
                "lehmer_pair": fine["lehmer_pair"],
            },
        ),
        GateRow(
            id="np15cwpr_10_nonpromotion",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="The Wronskian reduction and finite crossing diagnostics do not prove phase-level avoidance.",
            formula="exact reduction plus finite diagnostics != corrected transversality",
            proof_boundary="Not Lambda<=0 or RH.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction",
        "date": "2026-07-17",
        "status": (
            "exact thin-collar Wronskian reduction of corrected critical C1 "
            "transversality; open arithmetic separation, not Lambda<=0 or RH"
        ),
        "proof_boundary": (
            "This artifact proves exact complex, polar, Wronskian, and robust "
            "collision identities. The numerical rows only reject an unjustified "
            "one-sided phase-sign shortcut. It does not prove the L>=50 thin-collar "
            "Wronskian target on either the live asymptotic layer or the bounded-L "
            "residual domains, corrected transversality, Lambda<=0, or RH."
        ),
        "exact": exact,
        "diagnostics": fine,
        "convergence": convergence,
        "rows": [asdict(row) for row in rows],
        "sources": [
            "outputs/jensen_window_pf_newman_polymath15_critical_transversality_target.md",
            "outputs/jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate.md",
            "outputs/jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem.md",
            "outputs/jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate.md",
            "https://arxiv.org/abs/1904.12438",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    diagnostics_payload = artifact["diagnostics"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Polymath-15 Critical Wronskian Phase Reduction",
            "",
            "Date: 2026-07-17",
            "",
            "Status: exact crossing-localized reduction of the corrected `C1`",
            "target. The arithmetic separation remains open; this is not a proof",
            "of `Lambda <= 0` or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction.py",
            "```",
            "",
            "## Complex Main",
            "",
            "Absorb one half of the explicit endpoint lift into the complex finite",
            "sum:",
            "",
            "```text",
            exact["complex_split"],
            exact["cartesian"],
            "```",
            "",
            "For `E=a*exp(i*theta)`, the exact polar form is",
            "",
            "```text",
            exact["polar"],
            exact["crossing"],
            "```",
            "",
            "Thus a corrected real-part crossing is multiple precisely when its",
            "complex phase is stationary, apart from the separately visible case",
            "`E=0`.",
            "",
            "## Robust Collision Gate",
            "",
            "The certified remainder caps imply",
            "",
            "```text",
            exact["collision_wronskian_cap"],
            "```",
            "",
            "Consequently it is enough to establish the disjunction",
            "",
            "```text",
            exact["thin_collar_target"],
            exact["explicit_target"],
            "```",
            "",
            "The target formulas remain valid on their original broader domain,",
            "while the current proof partition is",
            "",
            "```text",
            exact["target_domains"],
            "```",
            "",
            "This localizes the arithmetic theorem to an exponentially thin collar",
            "around corrected real-part crossings. It does not require a global sign",
            "for the phase derivative.",
            "",
            "## Sign Shortcut Rejected",
            "",
            "At moderate frequency and `t=0`, direct corrected evaluations give",
            "both signs of `theta'`. The same occurs at two certified-to-high-precision",
            "corrected crossings:",
            "",
            *[
                "- `x={root}`: `theta'={velocity}`, residual `{residual}`".format(
                    root=row["root"],
                    velocity=row["phase_velocity"],
                    residual=row["abs_corrected_value_residual"],
                )
                for row in diagnostics_payload["root_rows"]
            ],
            "",
            "These finite diagnostics only shape the theorem search. They show that",
            "one-sided phase monotonicity is not supplied by the algebra alone; they",
            "do not prove or disprove the `L>=50` thin-collar separation.",
            "",
            "## Dynamical Bridge",
            "",
            "At a corrected crossing the phase Wronskian is also the exact",
            "denominator of the zero-flow generator:",
            "",
            "```text",
            exact["generator_bridge"],
            exact["pair_wronskian"],
            "```",
            "",
            "The corrected Lehmer-pair diagnostic has gap",
            f"`{diagnostics_payload['lehmer_pair']['gap']}` and minimum absolute",
            "phase velocity",
            f"`{diagnostics_payload['lehmer_pair']['min_abs_phase_velocity']}`.",
            "This is the finite-main signature of the same square-root velocity",
            "blow-up identified by the original dynamical-rigidity programme.",
            "",
            "## Live Target",
            "",
            "```text",
            exact["critical_value_form"],
            "```",
            "",
            "The problem is now a quantitative phase-critical-value avoidance theorem",
            "for the corrected Riemann-Siegel complex main. That theorem remains",
            "RH-level and is not established here.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman Polymath-15 critical Wronskian phase reduction: "
        f"{len(artifact['rows'])} rows, 6 exact identities, 1 robust collision gate, "
        "2 precision-stable crossing diagnostics, 1 open thin-collar target"
    )


if __name__ == "__main__":
    main()
