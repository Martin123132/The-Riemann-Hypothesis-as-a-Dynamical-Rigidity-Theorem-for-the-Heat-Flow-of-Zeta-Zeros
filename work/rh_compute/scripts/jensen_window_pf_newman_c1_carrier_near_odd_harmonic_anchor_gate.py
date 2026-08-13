#!/usr/bin/env python3
"""Extract the odd-harmonic endpoint anchor in the live C1 ideal package."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_carrier_near_odd_harmonic_anchor_gate"
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
DATE = "2026-08-04"

SOURCES = {
    "finite_band": (
        "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_finite_band_"
        "dirichlet_cell_reduction_gate.json"
    ),
    "symmetric_outer": (
        "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_ideal_cubic_"
        "symmetric_outer_pairing_gate.json"
    ),
    "terminal_lift": (
        "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_terminal_"
        "nonterminal_relative_lift_gate.json"
    ),
    "physical_chart": (
        "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_five_moment_q1_saddle_phase_variance_"
        "reduction.json"
    ),
    "frontier": (
        "work/rh_compute/results/"
        "jensen_window_pf_newman_c1_remainder_supersession_frontier_gate.json"
    ),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_path(key: str) -> Path:
    return REPO_ROOT / SOURCES[key]


def load_source(key: str) -> dict:
    return json.loads(source_path(key).read_text(encoding="utf-8"))


def contains(payload: dict, needle: str) -> bool:
    return needle in json.dumps(payload, sort_keys=True)


def odd_reciprocal_sum(m_value: int, n_value: int) -> sp.Rational:
    return sum(
        (sp.Rational(1, value) for value in range(m_value, n_value + 1) if value % 2),
        sp.Rational(0),
    )


def lower_half_mode(mode: int) -> sp.Expr:
    if mode == 0:
        return sp.Rational(1, 2)
    if mode % 2 == 0:
        return sp.Integer(0)
    return sp.I / (sp.pi * mode)


def upper_half_mode(mode: int) -> sp.Expr:
    if mode == 0:
        return sp.Rational(1, 2)
    if mode % 2 == 0:
        return sp.Integer(0)
    return -sp.I / (sp.pi * mode)


def half_cell_checks() -> int:
    cases = 0
    for m_value in (1, 2, 3, 4, 7, 10):
        for width in (0, 1, 2, 3, 5, 8, 13, 21):
            n_value = m_value + width
            modes = range(-(m_value - 1), n_value + 1)
            d_value = odd_reciprocal_sum(m_value, n_value)
            lower = sp.simplify(sum(lower_half_mode(mode) for mode in modes))
            upper = sp.simplify(sum(upper_half_mode(mode) for mode in modes))
            require(
                sp.simplify(lower - (sp.Rational(1, 2) + sp.I * d_value / sp.pi)) == 0,
                f"lower half-cell identity failed at m={m_value}, n={n_value}",
            )
            require(
                sp.simplify(upper - (sp.Rational(1, 2) - sp.I * d_value / sp.pi)) == 0,
                f"upper half-cell identity failed at m={m_value}, n={n_value}",
            )
            cases += 1
    return cases


def physical_roster_checks() -> int:
    cases = 0
    for n_cutoff in (4, 5, 8, 16, 50):
        for theta in (Fraction(0), Fraction(1, 3), Fraction(1)):
            a_value = Fraction(n_cutoff) + theta
            for epsilon in (Fraction(1, 4), Fraction(3, 4)):
                alpha = a_value * a_value - epsilon
                quotient = alpha / Fraction(2 * n_cutoff + 1, 2)
                m_value = quotient.numerator // quotient.denominator + 1
                doubled = 2 * alpha
                n_value = doubled.numerator // doubled.denominator
                d_value = sum(
                    (Fraction(1, value) for value in range(m_value, n_value + 1) if value % 2),
                    Fraction(0),
                )
                require(Fraction(n_value + 1, m_value + 1) > a_value, "roster ratio stress failed")
                require(float(d_value) > math.log(float(a_value)) / 2, "odd-harmonic stress failed")
                cases += 1
    return cases


def symbolic_anchor_checks() -> dict[str, str]:
    u, u_n, u_x, ell, d_odd, e_n, s_1 = sp.symbols(
        "u u_N u_x ell D_N E_N S_1", real=True
    )
    b = ell - u_n
    d = ell + u_n
    p_h = sp.I * (u - u_n) * (u + u_n) ** 2 / 4
    p_t = sp.I * (u + u_n) * (u - u_n) ** 2 / 4 - u_x * (u + u_n)

    p_h_1 = sp.expand(p_h.subs(u, ell))
    p_h_n = sp.expand(p_h.subs(u, u_n))
    p_t_1 = sp.expand(p_t.subs(u, ell))
    p_t_n = sp.expand(p_t.subs(u, u_n))
    require(sp.simplify(p_h_1 - sp.I * b * d**2 / 4) == 0, "Hermitian lower endpoint drifted")
    require(p_h_n == 0, "Hermitian terminal endpoint drifted")
    require(sp.simplify(p_t_1 - (sp.I * d * b**2 / 4 - u_x * d)) == 0, "transpose lower endpoint drifted")
    require(sp.simplify(p_t_n + 2 * u_x * u_n) == 0, "transpose terminal endpoint drifted")

    h_endpoint = sp.expand(sp.I * d_odd * (p_h_1 - e_n * p_h_n) / sp.pi)
    h_expected = -d_odd * b * d**2 / (4 * sp.pi)
    require(sp.simplify(h_endpoint - h_expected) == 0, "Hermitian anchor drifted")

    kappa = 1 / (2 * sp.pi * sp.I)
    b_t_n = -sp.I * u_x
    terminal = sp.simplify(2 * sp.I * u_n * e_n * kappa * b_t_n * s_1)
    t_endpoint = sp.expand(
        sp.I * d_odd * (p_t_1 - e_n * p_t_n) / sp.pi + terminal
    )
    t_expected = (
        -d_odd * d * b**2 / (4 * sp.pi)
        - sp.I * d_odd * u_x * d / sp.pi
        + sp.I * e_n * u_x * u_n * (2 * d_odd - s_1) / sp.pi
    )
    require(sp.simplify(t_endpoint - t_expected) == 0, "transpose anchor drifted")

    return {
        "P_H_lower": "i*log(N)*log(a^2/N)^2/4",
        "P_H_terminal": "0",
        "P_T_lower": "i*log(a^2/N)*log(N)^2/4-u_(N,x)*log(a^2/N)",
        "P_T_terminal": "-2*u_N*u_(N,x)",
        "hermitian_endpoint": "-D_N*log(N)*log(a^2/N)^2/(4*pi)",
        "transpose_endpoint": (
            "-D_N*log(a^2/N)*log(N)^2/(4*pi)"
            "-i*D_N*u_(N,x)*log(a^2/N)/pi"
            "+i*E_N*u_N*u_(N,x)*(2D_N-S_(1,N))/pi"
        ),
    }


def build_payload() -> dict:
    finite = load_source("finite_band")
    outer = load_source("symmetric_outer")
    terminal = load_source("terminal_lift")
    physical = load_source("physical_chart")
    frontier = load_source("frontier")

    require(
        contains(finite, "H_N:=integral_0^(1/2)K_N(v)dv=2kappa"),
        "finite-band half kernel drifted",
    )
    require(
        contains(finite, "mathscr M_N[P]-mathscr D_N[P]=2mathscr M_N[P]"),
        "finite-band join drifted",
    )
    require(
        outer["counts"]["signed_complete_ideal_bounds"] == 0
        and contains(outer, "K_near,N(u)=sum_(j=-(m_N-1))^n_N e(ju)"),
        "symmetric outer frontier drifted",
    )
    require(
        contains(terminal, "T_mu[P]=e(alpha_Plog mu)kappa"),
        "terminal conditional operator drifted",
    )
    require(
        contains(physical, "a^2=exp(L)+1/(32L^2)"),
        "physical q=1 scale drifted",
    )
    require(
        frontier["summary"]["required_c2_boundary_bounds"] == 0
        and frontier["summary"]["open_rows"] == 3,
        "C1 frontier drifted",
    )

    half_cases = half_cell_checks()
    roster_cases = physical_roster_checks()
    anchors = symbolic_anchor_checks()

    payloads = {
        "finite_band": finite,
        "symmetric_outer": outer,
        "terminal_lift": terminal,
        "physical_chart": physical,
        "frontier": frontier,
    }
    source_audit = {}
    for key, source in payloads.items():
        path = source_path(key)
        source_audit[key] = {
            "path": SOURCES[key],
            "sha256": sha256_path(path),
            "kind": source["kind"],
            "status": source["status"],
        }

    rows = [
        {
            "id": "cnha_01_sources",
            "readiness": "proved",
            "claim": "Five parent certificates are hash-pinned and replay-compatible.",
        },
        {
            "id": "cnha_02_odd_band",
            "readiness": "proved",
            "claim": "D_N is exactly the odd reciprocal sum over m_N<=r<=n_N.",
        },
        {
            "id": "cnha_03_lower_half",
            "readiness": "proved",
            "claim": "The lower near half-cell has integral 1/2+iD_N/pi.",
        },
        {
            "id": "cnha_04_upper_half",
            "readiness": "proved",
            "claim": "The upper near half-cell has integral 1/2-iD_N/pi.",
        },
        {
            "id": "cnha_05_near_decomposition",
            "readiness": "proved",
            "claim": "The near sum equals the starred carrier plus iD_N(f(1)-f(N))/pi plus one zero-at-cell-centres variation.",
        },
        {
            "id": "cnha_06_complete_recomposition",
            "readiness": "proved",
            "claim": "After the paired remote correction, the complete join is 2M_N+iD_N(f(1)-f(N))/pi-V_N.",
        },
        {
            "id": "cnha_07_roster_ratio",
            "readiness": "proved",
            "claim": "On the physical chart, (n_N+1)/(m_N+1)>a.",
        },
        {
            "id": "cnha_08_odd_lower_bound",
            "readiness": "proved",
            "claim": "On the fixed physical q=2tL^2=1 chart, D_N>(1/2)log(a)>L/4.",
        },
        {
            "id": "cnha_09_log_geometry",
            "readiness": "proved",
            "claim": "log N>L/2-1 and log(a^2/N)>L/2 on L>=50.",
        },
        {
            "id": "cnha_10_hermitian_endpoints",
            "readiness": "proved",
            "claim": "The ideal Hermitian lower endpoint is positive imaginary and its terminal endpoint vanishes.",
        },
        {
            "id": "cnha_11_hermitian_anchor",
            "readiness": "proved",
            "claim": "The raw Hermitian join is -A_H+R_H with A_H=D_N log(N)log(a^2/N)^2/(4pi)>0 and R_H=2M_H-V_H.",
        },
        {
            "id": "cnha_12_transpose_endpoints",
            "readiness": "proved",
            "claim": "Both ideal transpose endpoint values and the conditional terminal survivor remain explicit.",
        },
        {
            "id": "cnha_13_transpose_anchor",
            "readiness": "proved",
            "claim": "The raw transpose join is -A_T+R_T+Q_T with an explicit h-suppressed complex Q_T.",
        },
        {
            "id": "cnha_14_anchor_scale",
            "readiness": "proved",
            "claim": "On q=1, A_H>L^3(L-2)/(128pi)>14900 and A_T>L^2(L-2)^2/(128pi)>14300 at L>=50.",
        },
        {
            "id": "cnha_15_magnitude_handoff",
            "readiness": "open",
            "claim": "A stationary cancellation estimate for R_H and R_T+Q_T would turn the raw anchors into phase-independent modulus gaps.",
        },
        {
            "id": "cnha_16_phase_guard",
            "readiness": "guard_validated",
            "claim": "The terminal multiplier rotates either negative raw anchor by the terminal phase, so no physical-current sign follows without a phase sector or the complete Hermitian/transpose projection.",
        },
        {
            "id": "cnha_17_cutoff_guard",
            "readiness": "guard_validated",
            "claim": "D_N jumps at odd roster transfers; only the complete mode-transfer invariant package may be used across adjacent cutoffs.",
        },
        {
            "id": "cnha_18_boundary",
            "readiness": "guard_validated",
            "claim": "The signed stationary residual and every RH-level conclusion remain open.",
        },
    ]

    return {
        "kind": STEM,
        "schema_version": 1,
        "date": DATE,
        "status": (
            "exact odd-harmonic endpoint anchors and fixed-q1 quantitative floors "
            "proved with physical phase/cutoff guards; the stationary "
            "carrier-variation estimate remains open"
        ),
        "source_sha256": sha256_path(Path(__file__).resolve()),
        "source_audit": source_audit,
        "exact": {
            "roster": {
                "m_N": "floor(alpha_P/(N+1/2))+1",
                "n_N": "floor(2alpha_P)",
                "D_N": "sum_(m_N<=r<=n_N, r odd)1/r",
            },
            "near_half_cells": {
                "lower": "integral_0^(1/2)K_(near,N)(v)dv=1/2+iD_N/pi",
                "upper": "integral_(-1/2)^0K_(near,N)(v)dv=1/2-iD_N/pi",
            },
            "near_decomposition": (
                "mathscr N_N[P]=mathscr M_N[P]+iD_N{f_P(1)-f_P(N)}/pi+V_(near,N)[P]"
            ),
            "near_variation": (
                "V_(near,N)[P]=integral_0^(1/2){f_P(1+v)-f_P(1)}K_(near,N)(v)dv"
                "+sum_(q=2)^(N-1)integral_(-1/2)^(1/2){f_P(q+v)-f_P(q)}K_(near,N)(v)dv"
                "+integral_(-1/2)^0{f_P(N+v)-f_P(N)}K_(near,N)(v)dv"
            ),
            "paired_variation": "V_(near,N)[P]+mathscr C_N[P]=-V_N[P]",
            "complete_join": (
                "mathcal J_P=2mathscr M_N[P]+iD_N{f_P(1)-f_P(N)}/pi-V_N[P]"
            ),
            "terminal_factor": "E_N=e(alpha_P*log(N))*exp(S(log(N)))",
            "physical_bound": {
                "scope": "fixed q=2tL^2=1 physical chart, L>=50",
                "roster_ratio": "(n_N+1)/(m_N+1)>a",
                "integral_bound": "D_N>=(1/2)log((n_N+1)/(m_N+1))",
                "conclusion": "D_N>(1/2)log(a)>L/4",
            },
            "anchors": anchors,
            "hermitian_split": {
                "A_H": "D_N*log(N)*log(a^2/N)^2/(4*pi)",
                "R_H": "2mathscr M_N[P_H^0]-V_N[P_H^0]",
                "identity": "mathcal J_H^0=-A_H+R_H",
                "lower_bound": "A_H>L^3*(L-2)/(128*pi)>14900 on q=1, L>=50",
            },
            "transpose_split": {
                "A_T": "D_N*log(a^2/N)*log(N)^2/(4*pi)",
                "R_T": "2mathscr M_N[P_T^0]-V_N[P_T^0]",
                "Q_T": (
                    "-iD_N*u_(N,x)*log(a^2/N)/pi"
                    "+iE_N*u_N*u_(N,x)*(2D_N-S_(1,N))/pi"
                ),
                "identity": "mathcal J_T^0=-A_T+R_T+Q_T",
                "lower_bound": "A_T>L^2*(L-2)^2/(128*pi)>14300 on q=1, L>=50",
            },
            "magnitude_handoff": {
                "hermitian": "|R_H|<=A_H-delta implies |mathcal J_H^0|>=delta",
                "transpose": "|R_T|+|Q_T|<=A_T-delta implies |mathcal J_T^0|>=delta",
            },
        },
        "analytic_proof": {
            "roster": (
                "From a^2-1<alpha_P<=a^2, a=N+theta, and 0<=theta<=1: "
                "n_N+1>2(a^2-1) and m_N+1<=a^2/(a-1/2)+2. For a>=4, "
                "a^3-3a^2-a+1>0, so their ratio exceeds a."
            ),
            "odd_integral": (
                "Writing the odd integers as 2k+1 and comparing the decreasing "
                "summand with its unit-interval integral gives "
                "D_N>=(1/2)log((n_N+1)/(m_N+1))."
            ),
            "logs": (
                "Since a^2>exp(L), log(a)>L/2. Also N>=a-1>a/e and N<=a, "
                "so log N>L/2-1 and log(a^2/N)>=log(a)>L/2."
            ),
            "rational_floor": (
                "Using pi<22/7 at L=50 gives A_H>328125/22>14900 and "
                "A_T>157500/11>14300; both polynomial lower bounds increase for L>=50."
            ),
        },
        "diagnostics": {
            "half_cell_cases": half_cases,
            "half_cell_identities": 2 * half_cases,
            "physical_roster_stress_cases": roster_cases,
        },
        "rows": rows,
        "summary": {
            "rows": len(rows),
            "half_cell_cases": half_cases,
            "half_cell_identities": 2 * half_cases,
            "physical_roster_stress_cases": roster_cases,
            "explicit_raw_anchors": 2,
            "phase_rotation_guards": 1,
            "cutoff_transfer_guards": 1,
            "open_stationary_residual_targets": 1,
        },
        "next_action": (
            "Use the exact negative anchors only as centers of phase-independent "
            "modulus disks. Derive a reciprocal-stationary decomposition of "
            "R_H=2M_H-V_H and R_T=2M_T-V_T that preserves diagonal/adjacent "
            "cancellation and the contact-band condition. First test on the fixed "
            "q=1 chart whether "
            "|R_H|<A_H or |R_T|+|Q_T|<A_T is even scale-compatible. Do not infer "
            "a physical-current sign from the unrotated anchor and do not freeze "
            "D_N across an odd roster transfer."
        ),
        "proof_boundary": (
            "This gate proves exact endpoint-anchor decompositions and quantitative "
            "raw anchor lower bounds on the fixed q=1 chart. It proves no extension "
            "of those floors to other parabolic frequencies, stationary carrier-variation "
            "bound, complete physical-current sign, contact-box exclusion, ordinary "
            "retained-main separation, Q209 shell, cofinal successor, Lambda<=0, "
            "RH, or prize-level conclusion."
        ),
    }


def success_line(payload: dict) -> str:
    summary = payload["summary"]
    return (
        "validated Newman C1 carrier-near odd-harmonic anchor gate: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['half_cell_identities']} half-cell identities, "
        "q=1 D_N>L/4, 2 explicit raw anchors, 1 phase-rotation guard, "
        "1 live stationary residual"
    )


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    lines = [
        "# Newman C1 Carrier-Near Odd-Harmonic Anchor",
        "",
        f"Date: {DATE}",
        "",
        "Status: exact-lemma gate extracting two raw endpoint anchors and their phase/cutoff guards; this is not a proof of RH.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"python work/rh_compute/scripts/{STEM}.py",
        f"python work/rh_compute/scripts/check_{STEM}.py",
        "```",
        "",
        "## Odd-Harmonic Coordinate",
        "",
        "Put",
        "",
        "```text",
        "D_N=sum_(m_N<=r<=n_N, r odd)1/r.",
        "```",
        "",
        "For the contiguous near kernel, direct integration of each Fourier",
        "character gives",
        "",
        "```text",
        exact["near_half_cells"]["lower"],
        exact["near_half_cells"]["upper"],
        "```",
        "",
        "Thus the starred endpoint halves and the near kernel compose exactly as",
        "",
        "```text",
        exact["near_variation"],
        exact["near_decomposition"],
        exact["paired_variation"],
        exact["complete_join"],
        "```",
        "",
        "No remote one-sided series has been separated.",
        "",
        "## Physical Size",
        "",
        "On the fixed physical q=2tL^2=1 chart with L>=50, the floor",
        "inequalities and an integral comparison over odd integers prove",
        "",
        "```text",
        exact["physical_bound"]["roster_ratio"],
        exact["physical_bound"]["integral_bound"],
        exact["physical_bound"]["conclusion"],
        "```",
        "",
        "## Hermitian Anchor",
        "",
        "The lower ideal endpoint is purely imaginary and the terminal endpoint",
        "vanishes. Therefore",
        "",
        "```text",
        exact["hermitian_split"]["identity"],
        f"A_H={exact['hermitian_split']['A_H']}",
        f"R_H={exact['hermitian_split']['R_H']}",
        exact["hermitian_split"]["lower_bound"],
        "```",
        "",
        "## Transpose Anchor",
        "",
        "Keeping the terminal conditional survivor gives",
        "",
        "```text",
        exact["terminal_factor"],
        exact["transpose_split"]["identity"],
        f"A_T={exact['transpose_split']['A_T']}",
        f"R_T={exact['transpose_split']['R_T']}",
        f"Q_T={exact['transpose_split']['Q_T']}",
        exact["transpose_split"]["lower_bound"],
        "```",
        "",
        "## What This Buys",
        "",
        "The anchors are useful without choosing a phase if the stationary",
        "remainders can be placed in smaller disks:",
        "",
        "```text",
        exact["magnitude_handoff"]["hermitian"],
        exact["magnitude_handoff"]["transpose"],
        "```",
        "",
        "They are not standalone sign terms. The terminal multiplier rotates the",
        "negative real centers, and D_N changes when an odd reciprocal mode crosses",
        "a roster boundary. The complete Hermitian/transpose projection and exact",
        "mode transfer must remain attached.",
        "",
        "## Next Stage",
        "",
        payload["next_action"],
        "",
        "## Pi Provenance",
        "",
        "The factor pi comes only from e(y)=exp(2pi i y): integrating an odd",
        "Fourier mode across a half-cell gives i/(pi r). The rational estimate",
        "pi<22/7 is used only to display conservative decimal-free anchor floors.",
        "No circle, polygon, or fitted plotting constant is introduced.",
        "",
        "## Proof Boundary",
        "",
        payload["proof_boundary"],
        "",
        success_line(payload),
        "",
    ]
    return "\n".join(lines)


def write_text_atomic(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(content, encoding="utf-8")
    temporary.replace(path)


def main() -> int:
    payload = build_payload()
    write_text_atomic(RESULT, json.dumps(payload, indent=2) + "\n")
    write_text_atomic(NOTE, render_note(payload))
    print(success_line(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
