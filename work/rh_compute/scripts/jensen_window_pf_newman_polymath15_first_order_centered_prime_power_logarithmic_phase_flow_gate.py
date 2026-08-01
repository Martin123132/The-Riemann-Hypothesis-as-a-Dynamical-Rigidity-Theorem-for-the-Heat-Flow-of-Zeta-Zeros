#!/usr/bin/env python3
"""Build the prime-power logarithmic phase-flow gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_centered_prime_power_"
    "logarithmic_phase_flow_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "complete_chain": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_complete_prime_power_"
        "chain_occupation_gate.json"
    ),
    "ray_flow": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_ray_bottom_logarithmic_flow_reduction.json"
    ),
    "geometric_phase": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_geometric_prime_power_"
        "phase_monotonicity_gate.json"
    ),
    "prime_power_heat": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_prime_power_heat_block_"
        "composition_guard.json"
    ),
    "critical_frame": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_real_residual_reduction.json"
    ),
    "normalized_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_normalized_prefix_phase_flux_reduction.json"
    ),
}


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCES.items()
    }


def audit_sources(payloads: dict[str, dict]) -> dict:
    chain = payloads["complete_chain"].get("chain", {})
    if "Q_(m,p)(w)=sum_(k=0)^R_m" not in chain.get(
        "chain_polynomial", ""
    ):
        raise RuntimeError("complete-chain polynomial drifted")

    ray = payloads["ray_flow"].get("exact", {})
    if "D_j=d/dL" not in ray.get("logarithmic_derivative", {}).get(
        "operator", ""
    ):
        raise RuntimeError("ray logarithmic derivative drifted")

    geometric = payloads["geometric_phase"].get("exact", {})
    if "22/15-sqrt(2)" not in geometric.get("thresholds", {}).get(
        "dyadic", ""
    ):
        raise RuntimeError("geometric dyadic margin drifted")

    heat = payloads["prime_power_heat"].get("exact", {})
    if "p^(-12/25)" not in heat.get("prime_power_contraction", ""):
        raise RuntimeError("prime-power contraction drifted")

    frame = payloads["critical_frame"].get("exact", {}).get(
        "critical_frame", {}
    )
    if "s_*'=t*D_x/4" not in frame.get("saddle_derivative", ""):
        raise RuntimeError("critical-frame saddle derivative drifted")

    prefix = payloads["normalized_prefix"].get("exact", {})
    if "|epsilon_n|<8446/x^2" not in prefix.get(
        "q_ge_1_inward", ""
    ):
        raise RuntimeError("normalized correction-current bound drifted")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCES.items()
        },
        "source_kinds": {
            key: payload.get("kind") for key, payload in payloads.items()
        },
    }


def numeric_audit() -> dict:
    mp.mp.dps = 100
    p = mp.mpf(2)
    h = mp.log(p)
    length = 8
    t = mp.mpf(1) / 4
    base = 562_538_277
    cutoff = 128 * base
    left_scale = mp.mpf(cutoff) ** 2 - mp.mpf(1) / 64
    left_L = mp.log(left_scale)
    left_x = 4 * mp.pi * left_scale

    geometric = [p ** (-mp.mpf(k) / 2) for k in range(length)]
    heated = [
        mp.exp(
            -mp.mpf(k) * h / 2
            + (mp.mpf(k) ** 2 - 14 * mp.mpf(k)) * h**2 / 16
        )
        for k in range(length)
    ]

    def aligned_speed(coefficients: list[mp.mpf]) -> mp.mpf:
        mass = mp.fsum(coefficients)
        moment = mp.fsum(
            mp.mpf(k) * coefficient
            for k, coefficient in enumerate(coefficients)
        )
        return 1 + moment / mass

    geometric_speed = aligned_speed(geometric)
    heated_speed = aligned_speed(heated)
    margin = mp.mpf(22) / 15 - mp.sqrt(2)
    return {
        "arithmetic": (
            "100-digit evaluation; the checker independently proves "
            "the displayed coarse inequalities with rational enclosures"
        ),
        "prime": 2,
        "block_length": length,
        "t": "1/4",
        "prime_free_base": base,
        "cutoff": cutoff,
        "L_left": mp.nstr(left_L, 60),
        "tL_left": mp.nstr(t * left_L, 60),
        "x_left": mp.nstr(left_x, 60),
        "u_cell_upper": mp.nstr(mp.log(1 + mp.mpf(1) / cutoff), 60),
        "geometric_margin_mu_8": mp.nstr(margin, 60),
        "geometric_aligned_speed": mp.nstr(geometric_speed, 60),
        "heated_aligned_speed": mp.nstr(heated_speed, 60),
        "aligned_speed_difference": mp.nstr(
            geometric_speed - heated_speed, 60
        ),
        "difference_over_mu_8": mp.nstr(
            (geometric_speed - heated_speed) / margin, 60
        ),
        "heated_coefficients": [
            mp.nstr(value, 50) for value in heated
        ],
    }


def build_rows() -> list[GateRow]:
    return [
        GateRow(
            "pplpf_01_chart",
            "fixed-chart domain",
            "proved",
            "Fix one complete prime-power chain inside one cutoff chart.",
            "D_j=x*partial_x at fixed t,N; p does not divide m; "
            "n_k=m*p^k, 0<=k<=M-1.",
            "Cutoff equality is owned by the canonical new-N chart.",
        ),
        GateRow(
            "pplpf_02_factorization",
            "exact chain algebra",
            "proved",
            "The actual heat-corrected chain has one polynomial factorization.",
            "Z_(m,p)=B_m*Q, Q=sum_k A_k, "
            "A_k=a_(m,k)*w_(m,p)^k.",
            "The d_n factors remain inside A_k.",
        ),
        GateRow(
            "pplpf_03_currents",
            "exact logarithmic derivative",
            "proved",
            "Each internal coefficient has an exact fixed-ray current.",
            "epsilon_n=d_(n,x)/(1+d_n); "
            "D_j A_k=x*(epsilon_(m*p^k)-k*h*s_*')*A_k.",
            "Here h=log(p) and t is fixed.",
        ),
        GateRow(
            "pplpf_04_moments",
            "exact definitions",
            "proved",
            "Two chain moments contain the full internal derivative.",
            "H=sum_k k*A_k; E=sum_k epsilon_(m*p^k)*A_k; "
            "D_j Q=x*E-x*h*s_*'*H.",
            "No argument branch or division is used.",
        ),
        GateRow(
            "pplpf_05_clock",
            "exact phase clock",
            "proved",
            "The moving base phase supplies a positive angular clock.",
            "s_*'=c+i*b; tau=-x*h*b>0; rho=c/(-b)>=0.",
            "The certified critical-frame formula has b<-1/2.",
        ),
        GateRow(
            "pplpf_06_reduced_derivative",
            "exact logarithmic derivative",
            "proved",
            "The chain derivative separates angular, radial, and correction currents.",
            "D_j Q=tau*[(i-rho)*H+E/(h*(-b))].",
            "The identity is exact on every fixed chart.",
        ),
        GateRow(
            "pplpf_07_shifted_speed",
            "exact phase derivative",
            "proved",
            "The shifted internal block P=wQ has an exact normalized ray speed.",
            "D_j arg(P)/tau=J_ray/|Q|^2, where "
            "J_ray=|Q|^2+Re(H*conj(Q))-rho*Im(H*conj(Q))"
            "+Im(E*conj(Q))/(h*(-b)).",
            "The quotient statement is used only when Q is nonzero; "
            "J_ray itself is division-free.",
        ),
        GateRow(
            "pplpf_08_frozen_current",
            "exact phase derivative",
            "proved",
            "Freezing the coefficient shape removes precisely two ray terms.",
            "J_ang=|Q|^2+Re(H*conj(Q)); "
            "partial_theta arg[z*Q_theta]=J_ang/|Q_theta|^2.",
            "A frozen angular theorem and an L-flow theorem are not identical.",
        ),
        GateRow(
            "pplpf_09_external_phase",
            "exact phase bookkeeping",
            "proved",
            "The p-free external factor has its own moving phase.",
            "D_j arg(B_m)=x*(beta_t'-b*log(m)); "
            "D_j arg(Z_(m,p))/tau=(beta_t'-b*log(m))/(-h*b)"
            "+J_ray/|Q|^2-1.",
            "An internal one-turn theorem is not a joined-chain theorem.",
        ),
        GateRow(
            "pplpf_10_heat_profile",
            "exact dimensionless reduction",
            "proved",
            "The correction-free heat block is a compact q-deformation of the geometric block.",
            "q=exp(-t*h^2/4), v=2*(u-delta_a)/h; "
            "C_k^0/C_0^0=p^(-k/2)*q^[k*(2M-2-k+v)].",
            "Here u=log(a/(m*p^(M-1))) and "
            "delta_a=log(a)-Re(alpha).",
        ),
        GateRow(
            "pplpf_11_geometric_boundary",
            "boundary identification",
            "proved",
            "The previous geometric model is exactly the q=1 boundary.",
            "t=0 => q=1 and C_k^0/C_0^0=p^(-k/2).",
            "For q<1 the multiplier need not be uniformly close to one.",
        ),
        GateRow(
            "pplpf_12_witness_chart",
            "complete chart witness",
            "proved",
            "A complete dyadic eight-level chain exists just above L=50.",
            "t=1/4, p=2, m=562538277, N=128*m=72004899456; "
            "a=N at L=log(N^2-1/64)>50 and tL<25.",
            "The base m is odd and the chain has exactly M=8 levels.",
        ),
        GateRow(
            "pplpf_13_phase_sweep",
            "exact phase-attainment theorem",
            "proved",
            "The base phase attains theta=0 modulo 2*pi inside the witness cell.",
            "Delta theta>2*pi*log(2)*(2N+1)>2*pi on "
            "[lambda_N,lambda_(N+1)).",
            "At the attained point 0<=u<log(1+1/N)<1/N.",
        ),
        GateRow(
            "pplpf_14_heat_separation",
            "rigorous coarse enclosure",
            "proved",
            "At aligned phase the heat-only block is far from the geometric speed.",
            "Theta_geo(0)>14/5; Theta_heat(0)<11/5; "
            "Theta_geo(0)-Theta_heat(0)>3/5.",
            "The checker proves these inequalities with exact rational "
            "enclosures for log(2), sqrt(2), and exp.",
        ),
        GateRow(
            "pplpf_15_actual_stability",
            "rigorous perturbation enclosure",
            "proved",
            "The actual d_n and radial-current terms cannot erase the witness separation.",
            "|Theta_actual_flow(0)-Theta_heat(0)|<10^(-8).",
            "This uses u<1/N, delta_a<1/(4x), |d_n|<2189/x, "
            "|epsilon_n|<8446/x^2, and x>12*(N^2-1).",
        ),
        GateRow(
            "pplpf_16_C1_rejection",
            "counterexample gate",
            "proved",
            "The old absolute geometric C1 perturbation target is false.",
            "|Theta_actual_flow'-Theta_geometric'|>1/2 "
            "while mu_8=22/15-sqrt(2)<1/15.",
            "The failure is caused by the retained heat profile, not by a numerical accident.",
        ),
        GateRow(
            "pplpf_17_replacement",
            "new theorem target",
            "open",
            "Replace absolute closeness by direct heat-block starlikeness.",
            "J_heat=|Q_heat|^2+Re(H_heat*conj(Q_heat))>0 "
            "on |z|=1 for the admissible long families.",
            "This is a one-sided inequality and may hold while the absolute "
            "difference from the geometric speed is order one.",
        ),
        GateRow(
            "pplpf_18_modulus_payoff",
            "conditional theorem",
            "open",
            "A quantitative starlike margin supplies its own modulus budget.",
            "Re[z*P_heat'/P_heat]>=kappa>0, P_heat'(0)=1 "
            "=> a normalized starlike block; combine its distortion "
            "margin with the d_n C1 bounds.",
            "A uniform kappa for all admissible p,M,q,v is not yet proved.",
        ),
        GateRow(
            "pplpf_19_family_split",
            "route decision",
            "open",
            "Test the compact heat current at the three sharp geometric thresholds.",
            "Start with (p,M)=(2,8),(3,4),(5,2), then prove length "
            "propagation or an analytic tail theorem.",
            "Short p=2,3 chains and singleton chains remain explicit.",
        ),
        GateRow(
            "pplpf_20_rejoin",
            "nonpromotion guard",
            "proved",
            "Even a proved internal heat-block theorem must be rejoined before counting.",
            "Retain B_m, all p-free bases, short and singleton chains, "
            "the recurrent endpoint, q=1, cutoff homotopies, and "
            "mathsf_X*(L*D_j mathcal_C_N-mathcal_C_N).",
            "No joined Abel gap or successor winding bound is promoted.",
        ),
    ]


def build_payload() -> dict:
    payloads = load_sources()
    source_audit = audit_sources(payloads)
    rows = build_rows()
    audit = numeric_audit()
    proof_boundary = (
        "This artifact proves the exact fixed-ray logarithmic derivative "
        "of an actual heat-corrected prime-power chain, separates frozen "
        "angular, radial, correction, and external-phase currents, gives "
        "the exact dimensionless heat profile, and rigorously rejects the "
        "previous absolute geometric C1 perturbation target with a complete "
        "dyadic cutoff-cell witness. It replaces that failed route by a "
        "direct heat-block starlikeness target. It does not prove that "
        "target uniformly, a corrected all-family phase margin, a joined "
        "p-free or endpoint theorem, the Xi Abel gap, H_j<3*pi/2, contact "
        "exclusion, Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion."
    )
    return {
        "kind": STEM,
        "date": "2026-07-30",
        "status": (
            "exact_prime_power_logarithmic_phase_flow_and_"
            "absolute_C1_perturbation_counterexample"
        ),
        "proof_boundary": proof_boundary,
        "source_audit": source_audit,
        "exact": {
            "domain": (
                "L>=50, 0<tL<=25, fixed t,N,p,m; p does not divide m; "
                "n_k=m*p^k, 0<=k<=M-1"
            ),
            "chain": {
                "factorization": (
                    "Z_(m,p)=B_m*Q, Q=sum_k A_k, "
                    "A_k=a_(m,k)*w_(m,p)^k"
                ),
                "moments": (
                    "H=sum_k k*A_k; "
                    "E=sum_k epsilon_(m*p^k)*A_k"
                ),
                "coefficient_current": (
                    "epsilon_n=d_(n,x)/(1+d_n); "
                    "D_j A_k=x*(epsilon_(m*p^k)-k*h*s_*')*A_k"
                ),
                "sum_current": "D_j Q=x*E-x*h*s_*'*H",
            },
            "phase_flow": {
                "clock": (
                    "s_*'=c+i*b; tau=-x*h*b>0; rho=c/(-b)>=0"
                ),
                "reduced_derivative": (
                    "D_j Q=tau*[(i-rho)*H+E/(h*(-b))]"
                ),
                "ray_current": (
                    "J_ray=|Q|^2+Re(H*conj(Q))"
                    "-rho*Im(H*conj(Q))"
                    "+Im(E*conj(Q))/(h*(-b))"
                ),
                "shifted_speed": "D_j arg(wQ)/tau=J_ray/|Q|^2",
                "frozen_current": (
                    "J_ang=|Q|^2+Re(H*conj(Q)); "
                    "partial_theta arg[z*Q_theta]=J_ang/|Q_theta|^2"
                ),
                "external_phase": (
                    "D_j arg(B_m)=x*(beta_t'-b*log(m)); "
                    "D_j arg(Z_(m,p))/tau="
                    "(beta_t'-b*log(m))/(-h*b)+J_ray/|Q|^2-1"
                ),
            },
            "heat_profile": {
                "parameters": (
                    "q=exp(-t*h^2/4), "
                    "v=2*(u-delta_a)/h, "
                    "u=log(a/(m*p^(M-1))), "
                    "delta_a=log(a)-Re(alpha)"
                ),
                "coefficient": (
                    "C_k^0/C_0^0="
                    "p^(-k/2)*q^[k*(2M-2-k+v)]"
                ),
                "geometric_boundary": (
                    "t=0 => q=1 and C_k^0/C_0^0=p^(-k/2)"
                ),
                "actual_correction": (
                    "C_k/C_0=(C_k^0/C_0^0)"
                    "*(1+d_(m*p^k))/(1+d_m)"
                ),
            },
            "dyadic_counterexample": {
                "chart": (
                    "t=1/4, p=2, m=562538277, "
                    "N=128*m=72004899456, M=8"
                ),
                "left_boundary": (
                    "a=N, L=lambda_N=log(N^2-1/64)>50, tL<25"
                ),
                "phase_sweep": (
                    "Delta theta>2*pi*log(2)*(2N+1)>2*pi; "
                    "there is theta=0 mod 2*pi in the fixed-N cell"
                ),
                "cell_shape": (
                    "At the aligned point 0<=u<log(1+1/N)<1/N"
                ),
                "ideal_coefficients": (
                    "c_k=exp[-k*log(2)/2"
                    "+(k^2-14k)*log(2)^2/16], 0<=k<=7"
                ),
                "coarse_speeds": (
                    "Theta_geo(0)>14/5; Theta_heat(0)<11/5; "
                    "Theta_geo(0)-Theta_heat(0)>3/5"
                ),
                "actual_perturbation": (
                    "|Theta_actual_flow(0)-Theta_heat(0)|<10^(-8)"
                ),
                "rejection": (
                    "|Theta_actual_flow'-Theta_geometric'|>1/2 "
                    "but mu_8=22/15-sqrt(2)<1/15"
                ),
            },
            "replacement_target": {
                "current": (
                    "J_heat=|Q_heat|^2"
                    "+Re(H_heat*conj(Q_heat))>0 on |z|=1"
                ),
                "families": (
                    "Begin with (p,M)=(2,8),(3,4),(5,2), "
                    "then prove length propagation or an analytic tail theorem"
                ),
                "payoff": (
                    "A quantitative normalized starlike-block theorem "
                    "supplies both one-sided phase speed and a modulus "
                    "budget for the d_n C1 perturbation"
                ),
                "rejoin": (
                    "Insert any proved long-block estimate into the full "
                    "p-free, short, singleton, recurrent-endpoint, cutoff, "
                    "and ray-defect composition"
                ),
            },
            "numeric_audit": audit,
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_logarithmic_phase_identities": 7,
            "dimensionless_heat_profiles": 1,
            "complete_dyadic_counterexamples": 1,
            "rejected_absolute_C1_targets": 1,
            "replacement_starlikeness_targets": 1,
            "proved_joined_abel_gaps": 0,
            "proved_successor_winding_bounds": 0,
        },
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    audit = exact["numeric_audit"]
    lines = [
        "# Prime-Power Logarithmic Phase-Flow Gate",
        "",
        "Date: 2026-07-30",
        "",
        "Status: exact fixed-ray chain current and a rigorous rejection of",
        "the absolute geometric C1 perturbation route. The replacement",
        "heat-block starlikeness theorem remains open. This gate does not",
        "promote RH or Lambda<=0.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"python work/rh_compute/scripts/{STEM}.py",
        f"python work/rh_compute/scripts/check_{STEM}.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        (
            "validated prime-power logarithmic phase-flow gate: "
            "20 rows, 7 exact logarithmic-phase identities, "
            "1 dimensionless heat profile, 1 complete dyadic "
            "counterexample, 1 rejected absolute C1 target, "
            "1 replacement starlikeness target, 0 joined Abel gaps, "
            "0 successor winding bounds"
        ),
        "```",
        "",
        "## Exact Chain Flow",
        "",
        "On one fixed cutoff chart put `h=log(p)` and",
        "",
        "```text",
        exact["chain"]["factorization"],
        exact["chain"]["moments"],
        exact["chain"]["coefficient_current"],
        exact["chain"]["sum_current"],
        "```",
        "",
        "With the positive moving angular clock,",
        "",
        "```text",
        exact["phase_flow"]["clock"],
        exact["phase_flow"]["reduced_derivative"],
        exact["phase_flow"]["ray_current"],
        exact["phase_flow"]["shifted_speed"],
        "```",
        "",
        "The frozen-circle current is instead",
        "",
        "```text",
        exact["phase_flow"]["frozen_current"],
        "```",
        "",
        "so a frozen angular theorem cannot be promoted silently to an",
        "`L`-flow theorem. The p-free factor is also moving:",
        "",
        "```text",
        exact["phase_flow"]["external_phase"],
        "```",
        "",
        "## Heat Coordinates",
        "",
        "After normalizing the correction-free coefficient at `k=0`,",
        "",
        "```text",
        exact["heat_profile"]["parameters"],
        exact["heat_profile"]["coefficient"],
        exact["heat_profile"]["geometric_boundary"],
        exact["heat_profile"]["actual_correction"],
        "```",
        "",
        "The exponent contains `k(2M-2-k+v)`. It can be order one or",
        "larger across a long block even when the radial drift per phase",
        "turn is tiny. Heat therefore need not be an absolute C1-small",
        "perturbation of the `t=0` geometric polynomial.",
        "",
        "## Complete Counterexample",
        "",
        "The obstruction occurs in an actual complete fixed-cutoff chain:",
        "",
        "```text",
        exact["dyadic_counterexample"]["chart"],
        exact["dyadic_counterexample"]["left_boundary"],
        exact["dyadic_counterexample"]["phase_sweep"],
        exact["dyadic_counterexample"]["cell_shape"],
        exact["dyadic_counterexample"]["ideal_coefficients"],
        "```",
        "",
        "At the aligned phase, exact rational enclosures prove",
        "",
        "```text",
        exact["dyadic_counterexample"]["coarse_speeds"],
        exact["dyadic_counterexample"]["actual_perturbation"],
        exact["dyadic_counterexample"]["rejection"],
        "```",
        "",
        "High-precision values, used only as a readable audit, are",
        "",
        "```text",
        f"L_left={audit['L_left']}",
        f"mu_8={audit['geometric_margin_mu_8']}",
        f"Theta_geo(0)={audit['geometric_aligned_speed']}",
        f"Theta_heat(0)={audit['heated_aligned_speed']}",
        f"difference={audit['aligned_speed_difference']}",
        f"difference/mu_8={audit['difference_over_mu_8']}",
        "```",
        "",
        "Thus the old absolute target is false by a wide margin. The",
        "heat deformation is not harmful in this witness: its aligned",
        "speed remains positive. What fails is the choice of norm.",
        "",
        "## Replacement Target",
        "",
        "The surviving one-sided theorem is",
        "",
        "```text",
        exact["replacement_target"]["current"],
        "```",
        "",
        exact["replacement_target"]["families"] + ".",
        "",
        exact["replacement_target"]["payoff"] + ".",
        "",
        "Even after that theorem is proved,",
        "",
        "```text",
        exact["replacement_target"]["rejoin"],
        "```",
        "",
        "The `pi` in the cutoff cell comes from the completed-zeta saddle",
        "normalization. The `2*pi` in the phase sweep is the period of",
        "`exp(i*theta)`. No prime polygon or block defines a new value of",
        "`pi`.",
        "",
        "## Boundary",
        "",
        payload["proof_boundary"],
        "",
    ]
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "built prime-power logarithmic phase-flow gate: "
        "20 rows, 7 exact logarithmic-phase identities, "
        "1 dimensionless heat profile, 1 complete dyadic counterexample, "
        "1 rejected absolute C1 target, 1 replacement starlikeness target, "
        "0 joined Abel gaps, 0 successor winding bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
