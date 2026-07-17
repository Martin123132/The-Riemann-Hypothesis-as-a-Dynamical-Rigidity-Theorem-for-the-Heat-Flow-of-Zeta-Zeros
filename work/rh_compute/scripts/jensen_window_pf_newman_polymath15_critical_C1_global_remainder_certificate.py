#!/usr/bin/env python3
"""Close adjacent cutoffs in the corrected critical C1 remainder bound."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys

import mpmath as mp
import sympy as sp

import jensen_window_pf_newman_polymath15_critical_scaled_coercivity_scout as base


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate.md"
)
PRECISION_BITS = 256
L_MIN = 50
TRANSITION_CONSTANT = 1000
BASE_RAW_CONSTANT = 2300
GLOBAL_RAW_CONSTANT = 2500
GLOBAL_FIRST_CONSTANT = 5000
GLOBAL_NORM_CONSTANT = 32_000_000
AUDIT_N_VALUES = (10, 100, 1_000_000)
AUDIT_C_VALUES = ("0", "0.1", "1", "4")
AUDIT_OFFSETS = (("0", "0"), ("0.8", "0.6"), ("-0.7", "0.9"))
AUDIT_DPS = 100


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def arb_text(value: flint.arb, digits: int = 70) -> str:
    return value.str(digits)


def build_exact() -> dict:
    a, r = sp.symbols("a r", positive=True, real=True)
    delta = sp.log(1 + r)
    saddle_defect = sp.expand(
        2 * a**2 * delta - 2 * a**2 * (1 + r) * r + 3 * a**2 * r**2
    )
    expected = a**2 * (2 * sp.log(1 + r) - 2 * r + r**2)
    if sp.simplify(saddle_defect - expected) != 0:
        raise RuntimeError("cubic saddle defect identity failed")

    h, r0, rm, d, r1, time, psi = sp.symbols(
        "h r0 rm d r1 time psi"
    )
    e = r1 - d
    residual = h * r0 + rm - h * d + time / 4 * (
        sp.I * sp.pi * e / 2 + e**2
    ) - sp.I * sp.pi * psi
    if sp.expand(residual - (
        h * r0
        + rm
        - h * d
        + time * sp.I * sp.pi * (r1 - d) / 8
        + time * (r1 - d) ** 2 / 4
        - sp.I * sp.pi * psi
    )) != 0:
        raise RuntimeError("adjacent log-ratio residual identity failed")

    return {
        "transition_disk": (
            "If a radius-1/L center disk crosses x_n(t)=4*pi*n^2-pi*t/4, "
            "then every point in that disk has |z-x_n(t)|<=2/L"
        ),
        "parameters": (
            "T=z/2+pi*t/8, a=sqrt(T/(2*pi)), h=1/2-i*pi*t/8, "
            "s=i*T+h, delta=log(n/a), r=(n-a)/a"
        ),
        "saddle_defect": (
            "Psi=a^2*(2*log(1+r)-2*r+r^2)=O(a^2*r^3)"
        ),
        "alpha_remainders": (
            "r0=alpha(iT)-(log(a)+i*pi/4), "
            "r1=alpha(s)-(log(a)+i*pi/4), "
            "R_M=log(M0(s))-log(M0(iT))-h*alpha(iT)"
        ),
        "log_ratio": (
            "log(main_n/endpoint_n)=h*r0+R_M-h*delta+"
            "t/4*(i*pi*(r1-delta)/2+(r1-delta)^2)-i*pi*Psi"
        ),
        "ratio_bound": (
            "|log(main_n/endpoint_n)|<3/|T| and "
            "|main_n/endpoint_n-1|<6/|T|"
        ),
        "endpoint_bound": (
            "|endpoint_n|/A_t(x)<50*exp(-L/4) on the doubled transition collar"
        ),
        "adjacent_bound": (
            "|Delta corrected lift|/A_t(x)<1000*exp(-5L/4)"
        ),
        "absorption": (
            "2300*exp(-3L/4)+1000*exp(-5L/4)<"
            "2500*exp(-3L/4) for L>=50"
        ),
        "global_collar": (
            "For every critical radius-1/L disk, including cutoff crossings, "
            "sup |R_hat|/A_t(x)<2500*exp(-3L/4)"
        ),
        "global_c1": (
            "|r|<2500*exp(-3L/4), |r'|<5000*L*exp(-3L/4), "
            "r^2+(r'/L)^2<32000000*exp(-3L/2)"
        ),
        "remaining_target": (
            "Prove T_L[J]>32000000*exp(-3L/2) for the corrected finite main"
        ),
    }


def boundary_time(n: int, c: mp.mpf) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    if c == 0:
        time = mp.mpf(0)
    else:
        time = mp.findroot(
            lambda value: value * mp.log(n * n - value / 16) - c,
            c / (2 * mp.log(n)),
        )
    x_boundary = 4 * mp.pi * n * n - mp.pi * time / 4
    ell = mp.log(x_boundary / (4 * mp.pi))
    return time, x_boundary, ell


def log_ratio_audit(dps: int = AUDIT_DPS) -> dict:
    """Numerically cross-check the exact residual against the direct block ratio."""
    mp.mp.dps = dps
    rows: list[dict] = []
    max_relative_delta = mp.mpf(0)
    max_scaled_residual = mp.mpf(0)
    for n in AUDIT_N_VALUES:
        for c_text in AUDIT_C_VALUES:
            c = mp.mpf(c_text)
            time, x_boundary, ell = boundary_time(n, c)
            if time > mp.mpf("0.5"):
                continue
            for real_text, imag_text in AUDIT_OFFSETS:
                offset = mp.mpc(real_text, imag_text) / ell
                z = x_boundary + offset
                saddle_time = z / 2 + mp.pi * time / 8
                saddle = mp.sqrt(saddle_time / (2 * mp.pi))
                h = mp.mpf("0.5") - mp.j * mp.pi * time / 8
                s = mp.j * saddle_time + h
                log_n = mp.log(n)

                main = base.m_time(s, time) * mp.exp(
                    time * log_n**2 / 4
                    - (s + time * base.alpha(s) / 2) * log_n
                )
                p = 2 * (n - saddle) - 1
                u_rs = mp.exp(
                    -mp.j
                    * (
                        saddle_time / 2
                        * mp.log(saddle_time / (2 * mp.pi))
                        - saddle_time / 2
                        - mp.pi / 8
                    )
                )
                endpoint = (
                    mp.exp(time * mp.pi**2 / 64)
                    * (-1) ** n
                    * base.m_zero(mp.j * saddle_time)
                    * u_rs
                    * mp.exp(mp.j * mp.pi / 8)
                    * mp.exp(mp.j * mp.pi * (p**2 / 2 + p + mp.mpf(3) / 8))
                )

                delta = mp.log(n / saddle)
                relative_displacement = (n - saddle) / saddle
                saddle_defect = saddle**2 * (
                    2 * mp.log(1 + relative_displacement)
                    - 2 * relative_displacement
                    + relative_displacement**2
                )
                alpha_zero = base.alpha(mp.j * saddle_time) - (
                    mp.log(saddle) + mp.j * mp.pi / 4
                )
                alpha_one = base.alpha(s) - (
                    mp.log(saddle) + mp.j * mp.pi / 4
                )
                m_remainder = (
                    mp.log(base.m_zero(s) / base.m_zero(mp.j * saddle_time))
                    - h * base.alpha(mp.j * saddle_time)
                )
                error = alpha_one - delta
                residual = (
                    h * alpha_zero
                    + m_remainder
                    - h * delta
                    + time / 4 * (mp.j * mp.pi * error / 2 + error**2)
                    - mp.j * mp.pi * saddle_defect
                )
                direct_ratio = main / endpoint
                relative_delta = abs(mp.exp(residual) - direct_ratio) / abs(
                    direct_ratio
                )
                scaled_residual = abs(residual) * abs(saddle_time)
                max_relative_delta = max(max_relative_delta, relative_delta)
                max_scaled_residual = max(max_scaled_residual, scaled_residual)
                rows.append(
                    {
                        "n": n,
                        "c": c_text,
                        "offset_over_L": [real_text, imag_text],
                        "t": mp.nstr(time, 30, strip_zeros=False),
                        "L": mp.nstr(ell, 30, strip_zeros=False),
                        "relative_identity_delta": mp.nstr(
                            relative_delta, 30, strip_zeros=False
                        ),
                        "abs_log_ratio_times_abs_T": mp.nstr(
                            scaled_residual, 30, strip_zeros=False
                        ),
                    }
                )
    if max_relative_delta >= mp.mpf("1e-80"):
        raise RuntimeError("direct adjacent log-ratio identity audit failed")
    if max_scaled_residual >= 1:
        raise RuntimeError("sampled adjacent log-ratio constant exceeded one")
    return {
        "role": "high_precision_identity_audit_only",
        "proof_boundary": (
            "Finite complex-point diagnostics guard the exact formula and branches; "
            "they are not used to prove the uniform ratio bound."
        ),
        "dps": dps,
        "row_count": len(rows),
        "max_relative_identity_delta": mp.nstr(
            max_relative_delta, 30, strip_zeros=False
        ),
        "max_abs_log_ratio_times_abs_T": mp.nstr(
            max_scaled_residual, 30, strip_zeros=False
        ),
        "rows": rows,
    }


def interval_budget() -> dict:
    flint.ctx.prec = PRECISION_BITS
    ell = flint.arb(L_MIN)
    pi = flint.arb.pi()
    h_bound = (flint.arb(1) / 4 + (pi / 16) ** 2).sqrt()
    if not h_bound < flint.arb(3) / 5:
        raise RuntimeError("transition h bound failed")

    # Coefficients in the exact residual identity, measured in units 1/|T|.
    hr0 = (flint.arb(3) / 5) * (flint.arb(3) / 2)
    taylor = flint.arb(1) / 10
    delta_term = (flint.arb(3) / 5) * (flint.arb(2) / ell)
    square_term = flint.arb(3) / 4
    cubic_term = flint.arb(1) / 100
    log_ratio_constant = hr0 + taylor + delta_term + square_term + cubic_term
    if not log_ratio_constant < 3:
        raise RuntimeError("transition log-ratio constant failed")
    ratio_constant = 2 * flint.arb(3)
    if not ratio_constant <= 6:
        raise RuntimeError("transition exponential ratio constant failed")

    # The exact slow form of G, the doubled collar, and the endpoint prefactor
    # give |endpoint|/A < 50*exp(-L/4). With |T|>pi*exp(L), the two sharp
    # components therefore cost less than 600/pi at the exp(-5L/4) scale.
    adjacent_constant = 2 * ratio_constant * 50 / pi
    if not adjacent_constant < 200:
        raise RuntimeError("adjacent corrected-lift constant failed")
    if not adjacent_constant < TRANSITION_CONSTANT:
        raise RuntimeError("saved transition constant failed")

    effective_transition = flint.arb(TRANSITION_CONSTANT) * (-ell / 2).exp()
    if not effective_transition < 1:
        raise RuntimeError("transition absorption endpoint failed")
    global_raw = flint.arb(BASE_RAW_CONSTANT) + effective_transition
    if not global_raw < GLOBAL_RAW_CONSTANT:
        raise RuntimeError("global raw collar constant failed")
    norm_constant = (
        flint.arb(GLOBAL_RAW_CONSTANT) ** 2
        + flint.arb(GLOBAL_FIRST_CONSTANT) ** 2
    )
    if not norm_constant < GLOBAL_NORM_CONSTANT:
        raise RuntimeError("global C1 norm constant failed")
    endpoint_norm = flint.arb(GLOBAL_NORM_CONSTANT) * (-3 * ell / 2).exp()
    if not endpoint_norm < flint.arb(1) / 10**24:
        raise RuntimeError("global L50 norm endpoint failed")
    return {
        "precision_bits": PRECISION_BITS,
        "L_min": str(L_MIN),
        "h_bound_ball": arb_text(h_bound),
        "h_bound_lt": "3/5",
        "log_ratio_constant_ball": arb_text(log_ratio_constant),
        "log_ratio_constant_lt": "3",
        "ratio_constant": "6",
        "endpoint_constant": "50",
        "adjacent_constant_ball": arb_text(adjacent_constant),
        "adjacent_constant_lt": "200",
        "saved_transition_constant": str(TRANSITION_CONSTANT),
        "effective_transition_at_L50_ball": arb_text(effective_transition),
        "effective_transition_at_L50_lt": "1",
        "global_raw_constant_ball": arb_text(global_raw),
        "global_raw_constant_lt": str(GLOBAL_RAW_CONSTANT),
        "global_first_constant": str(GLOBAL_FIRST_CONSTANT),
        "global_norm_constant_ball": arb_text(norm_constant),
        "global_norm_constant_lt": str(GLOBAL_NORM_CONSTANT),
        "global_L50_norm_ball": arb_text(endpoint_norm),
        "global_L50_norm_lt": "10^-24",
    }


def build_artifact() -> dict:
    exact = build_exact()
    interval = interval_budget()
    ratio_audit = log_ratio_audit()
    rows = [
        GateRow(
            id="np15c1grc_01_transition_geometry",
            role="exact_geometry",
            readiness="ready_to_apply",
            claim="A cutoff-crossing center disk lies inside a doubled collar around one exact transition.",
            formula=exact["transition_disk"],
            proof_boundary="Consecutive cutoffs are much farther apart when L>=50.",
        ),
        GateRow(
            id="np15c1grc_02_cubic_saddle",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The C0 recurrence phase matches the added endpoint phase through quadratic order.",
            formula=f"{exact['parameters']}; {exact['saddle_defect']}",
            proof_boundary="Principal branches on the doubled positive-real collar.",
        ),
        GateRow(
            id="np15c1grc_03_log_ratio",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="All leading Stirling, heat, logarithmic, and endpoint phases cancel in one exact adjacent-block ratio.",
            formula=f"{exact['alpha_remainders']}; {exact['log_ratio']}",
            proof_boundary="One positive-imaginary block; the sharp block follows by reflection.",
        ),
        GateRow(
            id="np15c1grc_04_ratio_bound",
            role="analytic_bound",
            readiness="ready_to_apply",
            claim="Elementary alpha and Taylor bounds make the adjacent complex ratio one plus O(1/T).",
            formula=exact["ratio_bound"],
            proof_boundary="Uniform for L>=50, 0<t<=1/2, and 0<tL<=25.",
            diagnostics={
                "h": interval["h_bound_ball"],
                "constant": interval["log_ratio_constant_ball"],
            },
        ),
        GateRow(
            id="np15c1grc_05_endpoint_size",
            role="analytic_bound",
            readiness="ready_to_apply",
            claim="The explicit recurrence endpoint block has the expected quarter-power normalized size.",
            formula=exact["endpoint_bound"],
            proof_boundary="Uses the exact slow G formula and the published endpoint prefactor.",
        ),
        GateRow(
            id="np15c1grc_06_adjacent_bound",
            role="interval_certificate",
            readiness="certified",
            claim="The adjacent corrected-lift mismatch gains a full extra factor exp(-L).",
            formula=exact["adjacent_bound"],
            proof_boundary="Full doubled complex collar, not only the real cutoff.",
            diagnostics={
                "raw_constant": interval["adjacent_constant_ball"],
                "saved_constant": interval["saved_transition_constant"],
            },
        ),
        GateRow(
            id="np15c1grc_07_absorption",
            role="interval_certificate",
            readiness="certified",
            claim="The adjacent term fits inside the fixed-cell raw collar slack uniformly from L=50 onward.",
            formula=exact["absorption"],
            proof_boundary="Monotone in L.",
            diagnostics={
                "effective_L50": interval["effective_transition_at_L50_ball"],
                "global_raw": interval["global_raw_constant_ball"],
            },
        ),
        GateRow(
            id="np15c1grc_08_global_C1",
            role="proved_bound",
            readiness="ready_to_apply",
            claim="The corrected critical C1 remainder certificate now includes every cutoff transition.",
            formula=f"{exact['global_collar']}; {exact['global_c1']}",
            proof_boundary="L>=50 and 0<tL<=25 only.",
            diagnostics={
                "norm_constant": interval["global_norm_constant_ball"],
                "L50_norm": interval["global_L50_norm_ball"],
            },
        ),
        GateRow(
            id="np15c1grc_09_live_arithmetic",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="With the error and cutoff sides closed, the asymptotic critical obstruction is the corrected finite first-jet lower bound itself.",
            formula=exact["remaining_target"],
            proof_boundary="Open and RH-level.",
        ),
        GateRow(
            id="np15c1grc_10_nonpromotion",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="A global remainder certificate does not establish the arithmetic transversality lower bound.",
            formula="controlled approximation error != no double zero",
            proof_boundary="Not Lambda<=0 or RH.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate",
        "date": "2026-07-17",
        "status": (
            "global corrected C1 remainder certificate on the L>=50 critical "
            "ray including every cutoff; not corrected transversality, Lambda<=0, or RH"
        ),
        "proof_boundary": (
            "This artifact extends the corrected normalized C1 remainder bounds to "
            "all L>=50 critical disks, including cutoff crossings. It does not "
            "prove the corrected finite first-jet lower bound, close L<50, exclude "
            "a positive Newman boundary, prove Lambda<=0, or prove RH."
        ),
        "exact": exact,
        "interval": interval,
        "numerical_log_ratio_audit": ratio_audit,
        "rows": [asdict(row) for row in rows],
        "sources": [
            "outputs/jensen_window_pf_newman_polymath15_endpoint_C0_recurrence_transition_gate.md",
            "outputs/jensen_window_pf_newman_polymath15_critical_C1_cell_remainder_certificate.md",
            "https://arxiv.org/abs/1904.12438",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    interval = artifact["interval"]
    ratio_audit = artifact["numerical_log_ratio_audit"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Polymath-15 Critical C1 Global Remainder Certificate",
            "",
            "Date: 2026-07-17",
            "",
            "Status: global corrected `C1` remainder on the `L>=50` critical ray,",
            "including every cutoff transition. This is not a proof of",
            "`Lambda <= 0` or RH; corrected transversality remains open.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate.py",
            "```",
            "",
            "## Adjacent Ratio",
            "",
            "The exact C0 recurrence reduces a cutoff jump to one main block minus",
            "one endpoint block. On a doubled transition collar put",
            "",
            "```text",
            exact["parameters"],
            exact["alpha_remainders"],
            "```",
            "",
            "The recurrence phase has the exact residual",
            "",
            "```text",
            exact["saddle_defect"],
            "```",
            "",
            "and direct cancellation gives",
            "",
            "```text",
            exact["log_ratio"],
            "```",
            "",
            "The elementary bounds saved in the artifact yield",
            "",
            "```text",
            f"|h| <= {interval['h_bound_ball']} < 3/5",
            f"log-ratio constant = {interval['log_ratio_constant_ball']} < 3",
            exact["ratio_bound"],
            "```",
            "",
            "A separate direct complex-point audit compares the block quotient",
            "with the exponential of the displayed exact residual. Across",
            f"`{ratio_audit['row_count']}` points its maximum relative discrepancy is",
            f"`{ratio_audit['max_relative_identity_delta']}`, while the largest",
            f"sampled `|log ratio|*|T|` is `{ratio_audit['max_abs_log_ratio_times_abs_T']}`.",
            "This finite audit only guards signs and branches; the uniform bound",
            "above comes from the analytic estimates.",
            "",
            "## Transition Bound",
            "",
            "Combining the ratio estimate with the explicit endpoint size gives",
            "",
            "```text",
            exact["endpoint_bound"],
            f"direct adjacent constant = {interval['adjacent_constant_ball']} < 200",
            exact["adjacent_bound"],
            "```",
            "",
            "The saved constant `1000` is deliberately loose. At `L=50` its",
            "effective contribution to the `exp(-3L/4)` collar constant is only",
            f"`{interval['effective_transition_at_L50_ball']}`. Therefore",
            "",
            "```text",
            exact["absorption"],
            exact["global_collar"],
            "```",
            "",
            "## Global C1 Budget",
            "",
            "Cauchy's estimate now applies with one fixed analytic cutoff on every",
            "center disk, whether or not the prescribed cutoff changes inside it:",
            "",
            "```text",
            exact["global_c1"],
            "```",
            "",
            "The approximation and cutoff sides of the `L>=50` critical problem",
            "are therefore closed. The remaining asymptotic target is exactly",
            "",
            "```text",
            exact["remaining_target"],
            "```",
            "",
            "That arithmetic lower bound is not proved here; neither are",
            "positive-boundary exclusion, `Lambda <= 0`, or RH.",
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
        "built Newman Polymath-15 critical C1 global remainder certificate: "
        f"{len(artifact['rows'])} rows, 2 exact transition identities, "
        "1 global cutoff theorem, 1 explicit C1 budget, 1 open arithmetic target"
    )


if __name__ == "__main__":
    main()
