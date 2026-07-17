#!/usr/bin/env python3
"""Derive and scout the adjacent-cutoff endpoint-correction recurrence."""

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
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_endpoint_C0_recurrence_transition_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_endpoint_C0_recurrence_transition_gate.md"
)
N_VALUES = (10, 100, 1000, 1_000_000, 100_000_000_000)
C_VALUES = ("0", "0.1", "1", "4", "25")


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
    p = sp.symbols("p", real=True)

    def c0(q):
        return (
            sp.exp(sp.pi * sp.I * (q**2 / 2 + sp.Rational(3, 8)))
            - sp.I * sp.sqrt(2) * sp.cos(sp.pi * q / 2)
        ) / (2 * sp.cos(sp.pi * q))

    recurrence = sp.trigsimp(
        sp.expand_trig(
            c0(p + 2)
            + c0(p)
            - sp.exp(
                sp.pi
                * sp.I
                * (p**2 / 2 + p + sp.Rational(3, 8))
            )
        ).rewrite(sp.exp)
    )
    if sp.simplify(recurrence) != 0:
        raise RuntimeError("C0 adjacent recurrence failed")

    return {
        "recurrence": (
            "C0(p+2)+C0(p)=exp(pi*i*(p^2/2+p+3/8))"
        ),
        "cutoff_parameters": (
            "p_N=1-2*(a-N); replacing N by N+1 sends p_N to p_N+2"
        ),
        "signed_jump": (
            "(-1)^(N+1)C0(p+2)-(-1)^N C0(p)="
            "(-1)^(N+1)*exp(pi*i*(p^2/2+p+3/8))"
        ),
        "corrected_jump": (
            "J_(N+1)-J_N equals one added two-saddle Dirichlet block minus "
            "the explicit signed endpoint-recurrence block"
        ),
        "boundary": (
            "At a=N+1, p_N=-1 and p_(N+1)=1"
        ),
        "t0_ratio": (
            "At t=0, T=2*pi*n^2: main/endpoint="
            "n^(-1/2)*M0(1/2+i*T)/(exp(pi*i/8)*M0(i*T))"
        ),
        "observed_scale": (
            "Delta J=O(exp(-5L/4-c*L/16)) and "
            "Delta J'/L=O(exp(-5L/4-c*L/16))"
        ),
        "live_bound": (
            "Prove |Delta J|+|Delta J'|/L<=exp(-5L/4) on every "
            "L>=50 radius-1/L cutoff collar with 0<tL<=25"
        ),
        "absorption": (
            "The live bound is absorbed by the existing 2500*exp(-3L/4) "
            "fixed-cell C1 collar budget"
        ),
    }


def fixed_correction(x: mp.mpf, t: mp.mpf, cutoff: int) -> mp.mpf:
    s = (1 - mp.j * x) / 2
    amplitude = abs(base.m_time(s, t))
    t_prime = x / 2 + mp.pi * t / 8
    saddle = mp.sqrt(t_prime / (2 * mp.pi))
    p = 1 - 2 * (saddle - cutoff)
    u = mp.exp(
        -mp.j
        * (
            t_prime / 2 * mp.log(t_prime / (2 * mp.pi))
            - t_prime / 2
            - mp.pi / 8
        )
    )
    raw = (
        2
        * ((-1) ** cutoff)
        * mp.exp(t * mp.pi**2 / 64)
        * mp.re(
            base.m_zero(mp.j * t_prime)
            * base.c_zero(p)
            * u
            * mp.exp(mp.j * mp.pi / 8)
        )
    )
    return raw / amplitude


def boundary_time(n: int, c: mp.mpf) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    if c == 0:
        t = mp.mpf(0)
    else:
        initial = c / (2 * mp.log(n))
        t = mp.findroot(
            lambda value: value * mp.log(n * n - value / 16) - c,
            initial,
        )
    x = 4 * mp.pi * n * n - mp.pi * t / 4
    ell = mp.log(x / (4 * mp.pi))
    return x, ell, t


def jump_point(n: int, c_text: str, dps: int) -> dict | None:
    mp.mp.dps = dps
    c = mp.mpf(c_text)
    x, ell, t = boundary_time(n, c)
    if t > mp.mpf("0.5"):
        return None

    s = (1 - mp.j * x) / 2
    alpha = base.alpha(s)
    alpha_first = base.alpha_prime(s)
    alpha_second = base.alpha_second(s)
    s_star = s + t * alpha / 2
    s_star_first = -mp.j / 2 - mp.j * t * alpha_first / 4
    s_star_second = -t * alpha_second / 8
    log_m_first = alpha + t * alpha * alpha_first / 2
    log_m_second = alpha_first + t * (
        alpha_first**2 + alpha * alpha_second
    ) / 2
    beta_first = -mp.re(log_m_first) / 2
    beta_second = -mp.im(log_m_second) / 4

    log_n = mp.log(n)
    term = mp.exp(t * log_n**2 / 4 - s_star * log_n)
    term_first = -s_star_first * log_n * term
    term_second = (
        s_star_first**2 * log_n**2 - s_star_second * log_n
    ) * term
    normalizer = base.m_time(s, t)
    phase = normalizer / abs(normalizer)
    main0 = 2 * mp.re(phase * term)
    main1 = 2 * mp.re(
        phase * (term_first + mp.j * beta_first * term)
    )
    main2 = 2 * mp.re(
        phase
        * (
            term_second
            + 2 * mp.j * beta_first * term_first
            + (mp.j * beta_second - beta_first**2) * term
        )
    )

    new_q = lambda value: fixed_correction(value, t, n)
    old_q = lambda value: fixed_correction(value, t, n - 1)
    correction0 = new_q(x) - old_q(x)
    correction1 = mp.diff(new_q, x) - mp.diff(old_q, x)
    correction2 = mp.diff(new_q, x, 2) - mp.diff(old_q, x, 2)
    jump0 = main0 - correction0
    jump1 = main1 - correction1
    jump2 = main2 - correction2
    scale = mp.exp(5 * ell / 4 + c * ell / 16)
    return {
        "n": n,
        "c": c_text,
        "L": mp.nstr(ell, 40),
        "t": mp.nstr(t, 40),
        "main_jump": mp.nstr(main0, 45),
        "correction_jump": mp.nstr(correction0, 45),
        "corrected_jump": mp.nstr(jump0, 45),
        "corrected_first_over_L": mp.nstr(jump1 / ell, 45),
        "corrected_second_over_L2": mp.nstr(jump2 / ell**2, 45),
        "scaled_value": mp.nstr(abs(jump0) * scale, 35),
        "scaled_first": mp.nstr(abs(jump1 / ell) * scale, 35),
    }


def diagnostics(coarse_dps: int, fine_dps: int) -> dict:
    coarse = [
        row
        for n in N_VALUES
        for c in C_VALUES
        if (row := jump_point(n, c, coarse_dps)) is not None
    ]
    fine = [
        row
        for n in N_VALUES
        for c in C_VALUES
        if (row := jump_point(n, c, fine_dps)) is not None
    ]
    if [(row["n"], row["c"]) for row in coarse] != [
        (row["n"], row["c"]) for row in fine
    ]:
        raise RuntimeError("transition scout key mismatch")
    max_relative = mp.mpf(0)
    for left, right in zip(coarse, fine, strict=True):
        for key in ("corrected_jump", "corrected_first_over_L"):
            lv = mp.mpf(left[key])
            rv = mp.mpf(right[key])
            if rv:
                max_relative = max(max_relative, abs(lv - rv) / abs(rv))
    if max_relative >= mp.mpf("1e-25"):
        raise RuntimeError("transition scout precision failed")
    max_scaled_value = max(mp.mpf(row["scaled_value"]) for row in fine)
    max_scaled_first = max(mp.mpf(row["scaled_first"]) for row in fine)
    if max_scaled_value >= 1 or max_scaled_first >= 1:
        raise RuntimeError("transition scout scale guard failed")
    return {
        "coarse_dps": coarse_dps,
        "fine_dps": fine_dps,
        "row_count": len(fine),
        "max_relative_delta": mp.nstr(max_relative, 25),
        "max_scaled_value": mp.nstr(max_scaled_value, 25),
        "max_scaled_first": mp.nstr(max_scaled_first, 25),
        "rows": fine,
    }


def build_artifact(coarse_dps: int, fine_dps: int) -> dict:
    exact = build_exact()
    diag = diagnostics(coarse_dps, fine_dps)
    rows = [
        GateRow(
            id="np15c0rtg_01_recurrence",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="C0 obeys an exact two-step recurrence across adjacent Riemann-Siegel cutoffs.",
            formula=exact["recurrence"],
            proof_boundary="Exact away from apparent poles and extended through them by removability.",
        ),
        GateRow(
            id="np15c0rtg_02_signed_jump",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The alternating endpoint correction converts the recurrence into one explicit signed jump block.",
            formula=f"{exact['cutoff_parameters']}; {exact['signed_jump']}",
            proof_boundary="Fixed adjacent integers N and N+1.",
        ),
        GateRow(
            id="np15c0rtg_03_corrected_jump",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim="The adjacent corrected-main mismatch is one Dirichlet block minus one explicit endpoint block.",
            formula=exact["corrected_jump"],
            proof_boundary="No long finite-sum subtraction remains.",
        ),
        GateRow(
            id="np15c0rtg_04_boundary_ratio",
            role="exact_specialization",
            readiness="ready_to_apply",
            claim="At t=0 and an exact cutoff, the complex block ratio is an elementary M0 quotient tending to one.",
            formula=f"{exact['boundary']}; {exact['t0_ratio']}",
            proof_boundary="t=0 boundary specialization only.",
        ),
        GateRow(
            id="np15c0rtg_05_transition_scout",
            role="high_precision_diagnostic",
            readiness="diagnostic_only",
            claim="Independent precision runs show cancellation through the first normalized jet over eleven decades of cutoff index.",
            formula=exact["observed_scale"],
            proof_boundary="Point diagnostics at exact cutoff boundaries.",
            diagnostics=diag,
        ),
        GateRow(
            id="np15c0rtg_06_scale_guard",
            role="diagnostic_gate",
            readiness="guard_validated",
            claim="Every sampled adjacent mismatch is below the proposed next-saddle scale with constant one.",
            formula="scaled value and first-jet maxima are both below 1",
            proof_boundary="Finite sample only.",
            diagnostics={
                "value": diag["max_scaled_value"],
                "first": diag["max_scaled_first"],
            },
        ),
        GateRow(
            id="np15c0rtg_07_live_bound",
            role="open_theorem_target",
            readiness="not_ready_to_apply",
            claim="The cutoff gap is reduced to one explicit complex adjacent-block estimate.",
            formula=exact["live_bound"],
            proof_boundary="Must hold on full complex collars, not only at real boundaries.",
        ),
        GateRow(
            id="np15c0rtg_08_absorption",
            role="conditional_closure",
            readiness="conditional_ready",
            claim="The target estimate is more than sufficient for the existing fixed-cell C1 budget.",
            formula=exact["absorption"],
            proof_boundary="Conditional on the open complex-collar estimate.",
        ),
        GateRow(
            id="np15c0rtg_09_no_large_subtraction",
            role="strategy_improvement",
            readiness="ready_to_apply",
            claim="Future interval work can evaluate the small adjacent ratio directly instead of subtracting two full corrected sums.",
            formula="one block minus one recurrence block",
            proof_boundary="Computational conditioning improvement, not a sign theorem.",
        ),
        GateRow(
            id="np15c0rtg_10_nonpromotion",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="The recurrence and point scout do not close cutoff collars or corrected transversality.",
            formula="finite boundary scout != complex-collar theorem",
            proof_boundary="Not Lambda<=0 or RH.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_polymath15_endpoint_C0_recurrence_transition_gate",
        "date": "2026-07-17",
        "status": (
            "exact adjacent C0 recurrence and high-precision corrected-jump scout; "
            "not a complex-collar transition theorem, Lambda<=0, or RH"
        ),
        "proof_boundary": (
            "This artifact proves the exact C0 recurrence and adjacent corrected-jump "
            "reduction. Its numerical rows are point diagnostics. It does not prove "
            "the uniform complex-collar jump estimate, close all cutoff transitions, "
            "prove corrected first-jet transversality, exclude a positive Newman "
            "boundary, prove Lambda<=0, or prove RH."
        ),
        "exact": exact,
        "diagnostics": diag,
        "rows": [asdict(row) for row in rows],
        "sources": [
            "outputs/jensen_window_pf_newman_polymath15_critical_C1_cell_remainder_certificate.md",
            "https://arxiv.org/abs/1904.12438",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    diag = artifact["diagnostics"]
    table = [
        "| n | c=tL | L | Delta J | Delta J'/L | scaled value | scaled first |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in diag["rows"]:
        table.append(
            "| {n} | {c} | {L} | {j} | {j1} | {s0} | {s1} |".format(
                n=row["n"],
                c=row["c"],
                L=mp.nstr(mp.mpf(row["L"]), 7),
                j=mp.nstr(mp.mpf(row["corrected_jump"]), 8),
                j1=mp.nstr(mp.mpf(row["corrected_first_over_L"]), 8),
                s0=mp.nstr(mp.mpf(row["scaled_value"]), 7),
                s1=mp.nstr(mp.mpf(row["scaled_first"]), 7),
            )
        )
    return "\n".join(
        [
            "# Jensen-Window PF Newman Polymath-15 Endpoint C0 Recurrence Transition Gate",
            "",
            "Date: 2026-07-17",
            "",
            "Status: exact adjacent-cutoff recurrence and corrected-jump diagnostics.",
            "This is not a proof of `Lambda <= 0` or RH; the complex-collar",
            "transition theorem remains open.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_polymath15_endpoint_C0_recurrence_transition_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_endpoint_C0_recurrence_transition_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_endpoint_C0_recurrence_transition_gate.py",
            "```",
            "",
            "## Exact Recurrence",
            "",
            "Direct simplification of the explicit C0 formula gives",
            "",
            "```text",
            exact["recurrence"],
            exact["cutoff_parameters"],
            exact["signed_jump"],
            "```",
            "",
            "Therefore",
            "",
            "```text",
            exact["corrected_jump"],
            "```",
            "",
            "The large adjacent corrected sums never need to be subtracted.",
            "",
            "## Boundary Scout",
            "",
            "At exact cutoff boundaries, independent precision runs give",
            "",
            *table,
            "",
            "Here the final columns multiply by `exp(5L/4+cL/16)`. Their maxima",
            f"are `{diag['max_scaled_value']}` for the value and",
            f"`{diag['max_scaled_first']}` for the first jet. This strongly supports",
            "the next-saddle scale but remains finite point evidence.",
            "",
            "## Live Estimate",
            "",
            "The cutoff problem is now the explicit statement",
            "",
            "```text",
            exact["live_bound"],
            "```",
            "",
            "If established,",
            "",
            "```text",
            exact["absorption"],
            "```",
            "",
            "The remaining difficulty is a uniform complex-ratio estimate for one",
            "block. Neither the recurrence nor the point scout proves that estimate,",
            "corrected first-jet transversality, `Lambda <= 0`, or RH.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--coarse-dps", type=int, default=80)
    parser.add_argument("--fine-dps", type=int, default=110)
    args = parser.parse_args()
    artifact = build_artifact(args.coarse_dps, args.fine_dps)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman Polymath-15 endpoint C0 recurrence transition gate: "
        f"{len(artifact['rows'])} rows, 1 exact recurrence, "
        f"{artifact['diagnostics']['row_count']} transition diagnostics, "
        "1 open complex-collar bound"
    )


if __name__ == "__main__":
    main()
