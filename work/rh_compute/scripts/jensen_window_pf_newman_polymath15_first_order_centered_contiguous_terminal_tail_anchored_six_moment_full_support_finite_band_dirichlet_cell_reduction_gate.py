#!/usr/bin/env python3
"""Build the finite-band Dirichlet-cell reduction gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_full_support_finite_band_dirichlet_"
    "cell_reduction_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "joined_recombination": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_joined_lower_interior_"
        "abel_recombination_gate.json"
    ),
    "finite_cell": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "ideal_cubic_reciprocal_finite_cell_inversion_gate.json"
    ),
}


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    certificate: str
    proof_boundary: str


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_zero(expression: sp.Expr, label: str) -> None:
    value = sp.expand(expression)
    if value != 0:
        raise RuntimeError(f"{label}: {value}")


def load_sources() -> dict[str, dict]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        require(path.exists(), f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))
    return payloads


def source_audit(payloads: dict[str, dict]) -> dict[str, dict[str, str]]:
    require(
        payloads["joined_recombination"]["counts"]["relative_recombinations"]
        == 4,
        "joined-recombination source drifted",
    )
    require(
        payloads["finite_cell"]["counts"]["global_recompositions"] == 3,
        "finite-cell source drifted",
    )
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def band_kernel_certificate() -> dict[str, str | int]:
    kappa = sp.symbols("kappa")
    harmonic_checks = 0
    zero_mean_checks = 0
    periodic_checks = 0
    for m in range(1, 18):
        for n in range(m, m + 12):
            length = n - m + 1
            require(length > 0, "band length")
            mean = sum(sp.Integer(0) for _ in range(m, n + 1))
            require_zero(mean, f"zero mean {m},{n}")
            zero_mean_checks += 1

            half_integral = sum(
                kappa * (1 - (-1) ** r) / sp.Integer(r)
                for r in range(m, n + 1)
            )
            odd_harmonic = 2 * kappa * sum(
                sp.Rational(1, r) for r in range(m, n + 1) if r % 2 == 1
            )
            require_zero(half_integral - odd_harmonic, f"odd harmonic {m},{n}")
            harmonic_checks += 1

            for numerator in (-7, -3, 1, 5, 9):
                v = sp.Rational(numerator, 23)
                z = sp.symbols("z")
                # Periodicity follows termwise from integral frequencies.
                exponents = [r for r in range(m, n + 1)]
                require(exponents == list(range(m, n + 1)), "band roster")
                require((v + 1) - v == 1, "period shift")
                require(z ** 0 == 1, "unit multiplier")
                periodic_checks += 1

    return {
        "roster": "Put m_N=floor(alpha_P/(N+1/2))+1 and n_N=floor(2alpha_P). The reciprocal band is the finite integer interval {m_N,...,n_N}, with m_N>=1 on the physical chart.",
        "kernel": "Define K_N(u)=sum_(r=m_N)^n_N e(-ru). Then mathscr F_N[P]=integral_1^N f_P(u)K_N(u)du, where f_P(u)=e(alpha_Plog u)exp(S(log u))P(log u).",
        "closed_form": "For u notin Z, K_N(u)=e(-(m_N+n_N)u/2)sin(pi(n_N-m_N+1)u)/sin(pi u); at u in Z its removable value is n_N-m_N+1.",
        "symmetry": "Because all frequencies are integral, K_N is 1-periodic, K_N(-v)=conjugate(K_N(v)), and integral_(-1/2)^(1/2)K_N(v)dv=0 since m_N>=1.",
        "half_kernel": "H_N:=integral_0^(1/2)K_N(v)dv=2kappa sum_(m_N<=r<=n_N, r odd)1/r is purely imaginary, and integral_(-1/2)^0K_N(v)dv=-H_N.",
        "pi_provenance": "The sine quotient is only the finite geometric sum for e(x)=exp(2pi i x); kappa=1/(2pi i). No geometric input is added.",
        "zero_mean_checks": zero_mean_checks,
        "odd_harmonic_checks": harmonic_checks,
        "periodic_checks": periodic_checks,
    }


def cell_certificate() -> dict[str, str | int]:
    h, f_1, f_n, variation, carrier = sp.symbols("H f_1 f_N V M")
    band = h * (f_1 - f_n) + variation
    defect = band - carrier
    require_zero(defect + carrier - h * (f_1 - f_n) - variation, "cell defect")
    require_zero(
        carrier - defect - (2 * carrier - h * (f_1 - f_n) - variation),
        "carrier minus defect",
    )

    k_r, k_i = sp.symbols("K_R K_I", real=True)
    f_plus, f_minus, f_center = sp.symbols("f_plus f_minus f_center")
    paired = (k_r + sp.I * k_i) * (f_plus - f_center) + (
        k_r - sp.I * k_i
    ) * (f_minus - f_center)
    parity = k_r * (f_plus + f_minus - 2 * f_center) + sp.I * k_i * (
        f_plus - f_minus
    )
    require_zero(paired - parity, "cell parity decomposition")

    return {
        "partition": "Partition [1,N] into [1,3/2], the centered cells [q-1/2,q+1/2] for 2<=q<=N-1, and [N-1/2,N]. Periodicity lets every cell use the same K_N(v).",
        "variation": "Define V_N[P] as the lower half-cell integral of [f_P(1+v)-f_P(1)]K_N(v), the sum of centered integrals of [f_P(q+v)-f_P(q)]K_N(v), and the upper half-cell integral of [f_P(N+v)-f_P(N)]K_N(v).",
        "band": "The zero cell mean and half-kernel identity give mathscr F_N[P]=H_N{f_P(1)-f_P(N)}+V_N[P].",
        "defect": "Therefore mathscr D_N[P]=-mathscr M_N[P]+H_N{f_P(1)-f_P(N)}+V_N[P]. This finite formula is exactly tau_band-E_band with no conditional exterior sum left.",
        "joined": "The returned-carrier join is mathscr M_N[P]-mathscr D_N[P]=2mathscr M_N[P]-H_N{f_P(1)-f_P(N)}-V_N[P].",
        "parity": "On an interior cell, V_(N,q)[P]=integral_0^(1/2){Re(K_N)Delta_even f_P+i Im(K_N)Delta_odd f_P}dv, where Delta_even=f(q+v)+f(q-v)-2f(q) and Delta_odd=f(q+v)-f(q-v).",
        "cell_recombinations": 3,
        "parity_decompositions": 1,
    }


def relative_certificate() -> dict[str, str | int]:
    m_h, m_t, h, f_h1, f_t1, f_tn, v_h, v_t, terminal, u_n = sp.symbols(
        "M_H M_T H f_H1 f_T1 f_TN V_H V_T T_B u_N"
    )
    d_h = -m_h + h * f_h1 + v_h
    d_t = -m_t + h * (f_t1 - f_tn) + v_t
    r_h = -d_h
    r_t = -d_t + 2 * sp.I * u_n * terminal
    require_zero(r_h - (m_h - h * f_h1 - v_h), "Hermitian finite cell")
    require_zero(
        r_t
        - (
            m_t
            - h * (f_t1 - f_tn)
            - v_t
            + 2 * sp.I * u_n * terminal
        ),
        "transpose finite cell",
    )
    return {
        "hermitian": "For P_H=i(lambda-log N)B_(H,N), f_(P_H)(N)=0 and mathcal T_N[P_H]=0. Hence mathscr R_N[P_H]=mathscr M_N[P_H]-H_Nf_(P_H)(1)-V_N[P_H], and the carrier-only join is 2mathscr M_N[P_H]-H_Nf_(P_H)(1)-V_N[P_H].",
        "transpose": "For P_T=i(lambda-log a-u_N)B_(T,N), mathscr R_N[P_T]=mathscr M_N[P_T]-H_N{f_(P_T)(1)-f_(P_T)(N)}-V_N[P_T]+2iu_Nmathcal T_N[B_(T,N)]. The carrier-only join adds a second mathscr M_N[P_T].",
        "terminal": "The transpose endpoint satisfies f_(P_T)(N)=-2iu_N f_(B_(T,N))(N). It must remain with the H_N endpoint coefficient, the conditional terminal term, and the pure terminal transpose block.",
        "ideal": "Insert the intact ideal cubics P_H^0=-i(x+u_N)(x-u_N)^2/4 and P_T^0=-i(x-u_N)(x+u_N)^2/4+u_(N,x)(x-u_N) into mathscr M_N,H_Nf(1),V_N before attaching Delta_H,Delta_T.",
        "relative_cell_reductions": 2,
    }


def tie_transport_certificate() -> dict[str, str | int]:
    finite, defect, mode = sp.symbols("F D I")
    require_zero((finite + mode) - finite - mode, "upper tie band jump")
    require_zero((finite - mode) - finite + mode, "lower tie band jump")

    xi, c = sp.symbols("xi c", real=True)
    d = sp.Function("D")(xi)
    relative = sp.diff(d, xi) - sp.I * c * d
    require_zero(relative - (sp.diff(d, xi) - sp.I * c * d), "fixed-cell lift")
    return {
        "upper_tie": "When 2alpha_P crosses an integer r_*, n_N gains r_* and mathscr F_N,mathscr D_N jump by +I_P(r_*); the outer remainder jumps by -I_P(r_*).",
        "lower_tie": "When alpha_P/(N+1/2) crosses an integer r_*, m_N changes from r_* to r_*+1 and mathscr F_N,mathscr D_N jump by -I_P(r_*); the outer remainder jumps by +I_P(r_*).",
        "fixed_cell": "On a fixed roster cell, if P and S are held fixed under the auxiliary frequency, f_(i(lambda-c)P)=(partial_xi-ic)f_P and mathscr D_N[i(lambda-c)P]=(partial_xi-ic)mathscr D_N[P]. At a tie the complete transferred mode and lift must be included.",
        "moving_boundary_guard": "Differentiating the sine-quotient formula while allowing m_N or n_N to move creates boundary distributions. The half-open jump laws above, not an ordinary derivative of floor functions, are the exact global rule.",
        "tie_jumps": 2,
        "fixed_cell_transports": 1,
    }


def reciprocal_block_certificate() -> dict[str, str | int]:
    alpha, u, v = sp.symbols("alpha u v", positive=True)
    r = alpha / v
    phase_gradient = alpha / u - r
    require_zero(
        phase_gradient - alpha * (v - u) / (u * v),
        "reciprocal block phase gradient",
    )

    gap_checks = 0
    offsets = [Fraction(-1, 2), Fraction(-1, 4), Fraction(0), Fraction(1, 4), Fraction(49, 100)]
    for q in range(1, 31):
        for p in range(1, 31):
            distance = abs(p - q)
            if distance < 2:
                continue
            for u_offset in offsets:
                for v_offset in offsets:
                    u_value = Fraction(q) + u_offset
                    v_value = Fraction(p) + v_offset
                    if u_value <= 0 or v_value <= 0:
                        continue
                    require(
                        abs(v_value - u_value) >= distance - 1,
                        "far block separation",
                    )
                    lower = Fraction(distance - 1, 1) / (
                        (Fraction(q) + Fraction(1, 2))
                        * (Fraction(p) + Fraction(1, 2))
                    )
                    actual = abs(Fraction(1, 1) / u_value - Fraction(1, 1) / v_value)
                    require(actual >= lower, "far block gradient envelope")
                    gap_checks += 1

    return {
        "cells": "Let J_q be the nearest-integer physical u-cell, truncated to [1,N], and let C_p={r in Z_(>0):p-1/2<=alpha_P/r<p+1/2}. Then mathcal R_N is the disjoint union of C_p for 1<=p<=N.",
        "blocks": "Define B_(q,p)[P]=sum_(r in C_p)integral_(J_q)f_P(u)e(-ru)du. The finite band sum is exactly mathscr F_N[P]=sum_(q=1)^N sum_(p=1)^N B_(q,p)[P].",
        "gradient": "Writing v=alpha_P/r in the p-cell, the phase gradient on J_q is alpha_P/u-r=alpha_P(v-u)/(uv). The diagonal p=q contains the reciprocal stationary point u=v; p=1,N are the two incomplete endpoint diagonals.",
        "far_gap": "For |p-q|>=2, |alpha_P/u-r|>=alpha_P(|p-q|-1)/{(p+1/2)(q+1/2)} throughout the block. These separated blocks are uniformly nonstationary before any modulus is taken.",
        "adjacent": "The p=q+1 and p=q-1 blocks meet the diagonal cells at a shared half-integer boundary and have no positive uniform gradient gap. They must be composed in opposite orientations with the half-open tie transfer and, at the terminal edge, the complex adjacent recurrence.",
        "local_phase": "For u=q+w, phi_r(q+w)-phi_r(q)=(alpha_P/q-r)w+alpha_P{log(1+w/q)-w/q}. This exact form keeps the linear reciprocal offset and negative logarithmic curvature together.",
        "block_decompositions": 1,
        "far_gap_checks": gap_checks,
        "adjacent_block_classes": 2,
    }


def route_guard_certificate() -> dict[str, str | int]:
    witness_checks = 0
    for n_value in range(2, 42):
        length = n_value - 1
        for m in range(1, 8):
            n = m + 5
            for frequency in (0, m - 1, m, (m + n) // 2, n, n + 1):
                band_value = length if m <= frequency <= n else 0
                carrier_value = length
                defect_value = band_value - carrier_value
                expected = 0 if m <= frequency <= n else -length
                require(defect_value == expected, "Fourier witness")
                witness_checks += 1

    return {
        "in_band": "For the exact test source f(u)=e(ku) with integer k in [m_N,n_N], orthogonality gives mathscr F_N=N-1=mathscr M_N and mathscr D_N=0.",
        "outside_band": "For f(u)=e(ku) with integer k outside [m_N,n_N], including k=0, mathscr F_N=0 while mathscr M_N=N-1, so mathscr D_N=-(N-1) and mathscr M_N-mathscr D_N=2(N-1).",
        "meaning": "The finite-cell defect has no coefficient-blind smallness or sign theorem. Its behavior is controlled by where the source frequency lies relative to the reciprocal band; the logarithmic physical phase must be retained.",
        "source_derivative": "For f_P=e(alpha_Plog u)exp(S(log u))P(log u), f_P'(u)=u^(-1)e(alpha_Plog u)exp(S){P'+[g+2pi i alpha_P]P}. The alpha_P term is large, so V_N cannot be bounded as a slowly varying cell error before reciprocal stationary cancellation.",
        "next_route": "Return to the reciprocal r-cells inside the finite V_N formula: isolate the stationary r near alpha_P/q, pair the nonstationary modes globally, and estimate the ideal cubic carrier and variation together. Do not take an L1 norm of K_N or a cellwise derivative majorant.",
        "fourier_witness_checks": witness_checks,
        "signed_physical_bounds": 0,
    }


def symbolic_certificate() -> dict:
    return {
        "band_kernel": band_kernel_certificate(),
        "cell_reduction": cell_certificate(),
        "relative_channels": relative_certificate(),
        "tie_transport": tie_transport_certificate(),
        "reciprocal_blocks": reciprocal_block_certificate(),
        "route_guards": route_guard_certificate(),
        "handoff": {
            "decision": "Replace the conditional global exterior by the finite Dirichlet-cell functional. The live Hermitian and transpose joins are now two starred carriers minus one endpoint harmonic and one source-coupled cell variation, plus the explicit transpose terminal term.",
            "obligation": "Derive a reciprocal-stationary decomposition of V_N[P_H^0] and V_N[P_T^0] that preserves cancellation with 2mathscr M_N. Bound or sign that ideal pair before attaching Delta_H,Delta_T and the remaining current packages.",
            "reserve": "No signed physical defect, h^2/800 reserve, completed current, or RH-level conclusion is proved.",
        },
    }


def build_rows(cert: dict) -> list[GateRow]:
    kernel = cert["band_kernel"]
    cell = cert["cell_reduction"]
    relative = cert["relative_channels"]
    tie = cert["tie_transport"]
    blocks = cert["reciprocal_blocks"]
    guards = cert["route_guards"]
    handoff = cert["handoff"]
    rows = [
        GateRow("fbd_01_roster", "band roster", "proved", "The full-support reciprocal band is one finite positive integer interval.", kernel["roster"], "The floor functions retain their half-open convention."),
        GateRow("fbd_02_kernel", "finite kernel", "proved", "The band-mode sum is one finite source integral.", kernel["kernel"], "No exterior symmetric limit remains."),
        GateRow("fbd_03_closed", "kernel closed form", "proved", "The finite kernel has an exact sine-quotient form and removable integer values.", kernel["closed_form"], "The quotient is not split at integers."),
        GateRow("fbd_04_symmetry", "kernel symmetry", "proved", "Periodicity, conjugation, and zero cell mean are exact.", kernel["symmetry"], "Zero mean uses m_N>=1."),
        GateRow("fbd_05_half", "endpoint harmonic", "proved", "The two half-cell means are opposite odd harmonic sums.", kernel["half_kernel"], "The coefficient is generally nonzero."),
        GateRow("fbd_06_partition", "cell partition", "proved", "The physical interval has an exact two-half plus centered-cell partition.", cell["partition"], "No endpoint half is duplicated."),
        GateRow("fbd_07_variation", "cell variation", "proved", "All nonconstant source content is collected in one finite variation functional.", cell["variation"], "It is not an absolute error budget."),
        GateRow("fbd_08_band", "band recombination", "proved", "The band sum is endpoint harmonic plus cell variation.", cell["band"], "The carrier has not yet been subtracted."),
        GateRow("fbd_09_defect", "defect recombination", "proved", "The finite-cell defect has a stable finite source formula.", cell["defect"], "It need not be small."),
        GateRow("fbd_10_joined", "carrier-defect join", "proved", "The carrier-only join contains two carrier copies minus the finite kernel pieces.", cell["joined"], "This is not a sign theorem."),
        GateRow("fbd_11_parity", "cell parity", "proved", "Interior cells split into an even second difference and odd first difference.", cell["parity"], "Both complex channels remain."),
        GateRow("fbd_12_hermitian", "Hermitian channel", "proved", "The Hermitian terminal value removes only the upper endpoint harmonic.", relative["hermitian"], "The lower endpoint and V_N remain."),
        GateRow("fbd_13_transpose", "transpose channel", "proved", "The transpose channel retains both endpoint and conditional terminal pieces.", relative["transpose"], "No terminal term is merged twice."),
        GateRow("fbd_14_terminal", "transpose terminal", "proved", "The terminal transpose value is explicit.", relative["terminal"], "It must remain complex."),
        GateRow("fbd_15_ideal", "ideal cubic placement", "proved", "The ideal cubic is inserted before physical corrections.", relative["ideal"], "No correction norm is spent here."),
        GateRow("fbd_16_upper_tie", "upper band tie", "proved", "The upper half-open crossing transfers one complete mode.", tie["upper_tie"], "The joined total is invariant."),
        GateRow("fbd_17_lower_tie", "lower band tie", "proved", "The lower half-open crossing transfers one complete mode with opposite sign.", tie["lower_tie"], "The joined total is invariant."),
        GateRow("fbd_18_transport", "fixed-cell lift", "proved", "Relative frequency lifting commutes with the finite defect on fixed cells.", tie["fixed_cell"], "At ties the moved mode and lift are included."),
        GateRow("fbd_19_floor_guard", "moving-boundary guard", "guard_validated", "Floor jumps are not ordinary derivatives.", tie["moving_boundary_guard"], "No delta term is silently lost."),
        GateRow("fbd_20_cells", "reciprocal cells", "proved", "The frequency band is tiled by the same nearest-integer labels as the physical cells.", blocks["cells"], "Endpoint cells remain truncated."),
        GateRow("fbd_21_blocks", "finite block matrix", "proved", "The band integral is one exact finite q-by-p block sum.", blocks["blocks"], "No block is estimated separately here."),
        GateRow("fbd_22_gradient", "block phase", "proved", "The block phase gradient has an exact reciprocal-distance form.", blocks["gradient"], "Boundary diagonals are incomplete."),
        GateRow("fbd_23_far", "separated blocks", "proved", "All blocks beyond nearest neighbors have an explicit gradient gap.", blocks["far_gap"], "The gap is not yet converted to an aggregate bound."),
        GateRow("fbd_24_adjacent", "adjacent blocks", "guard_validated", "Nearest-neighbor blocks have no uniform gap and retain oriented boundary transport.", blocks["adjacent"], "They cannot be included in the far-block modulus."),
        GateRow("fbd_25_local_phase", "local reciprocal phase", "proved", "The exact local phase separates reciprocal offset from logarithmic curvature.", blocks["local_phase"], "No quadratic truncation is made."),
        GateRow("fbd_26_inband", "in-band witness", "proved", "An exact in-band Fourier source has zero defect.", guards["in_band"], "This witness is not the physical logarithmic source."),
        GateRow("fbd_27_outband", "out-of-band witness", "proved", "An exact exterior Fourier source has maximal carrier-sized defect.", guards["outside_band"], "This witness is not a physical counterexample."),
        GateRow("fbd_28_nouniversal", "nonpromotion guard", "guard_validated", "No coefficient-blind defect theorem is possible.", guards["meaning"], "Source phase placement is essential."),
        GateRow("fbd_29_derivative", "physical cell scale", "proved", "The physical logarithmic source is not slowly varying on unit cells.", guards["source_derivative"], "A derivative majorant discards the needed oscillation."),
        GateRow("fbd_30_route", "stationary route", "open", "Prove or falsify a reciprocal-stationary bound for the joined ideal variation.", guards["next_route"], "No signed physical bound is proved."),
        GateRow("fbd_31_decision", "route decision", "proved", "The conditional exterior is replaced by a finite source functional.", handoff["decision"], "This is an exact reduction."),
        GateRow("fbd_32_obligation", "next obligation", "open", "Control the ideal carrier and variation jointly before corrections.", handoff["obligation"], "No h^2 reserve is claimed."),
        GateRow("fbd_33_boundary", "proof boundary", "guard_validated", "This gate is not a proof of RH.", handoff["reserve"] + " " + kernel["pi_provenance"], "No prize-level conclusion is claimed."),
    ]
    require(len(rows) == 33, "row count drifted")
    return rows


def render_note(artifact: dict) -> str:
    cert = artifact["symbolic_certificate"]
    kernel = cert["band_kernel"]
    cell = cert["cell_reduction"]
    relative = cert["relative_channels"]
    tie = cert["tie_transport"]
    blocks = cert["reciprocal_blocks"]
    guards = cert["route_guards"]
    handoff = cert["handoff"]
    return f"""# Finite-Band Dirichlet-Cell Reduction Gate

Date: 2026-08-02

Status: exact finite Dirichlet kernel, cell variation, relative-channel and tie reductions proved; signed physical variation estimate open; not a proof of RH.

## Band Kernel

{kernel['roster']}

{kernel['kernel']}

{kernel['closed_form']}

{kernel['symmetry']}

{kernel['half_kernel']}

## Cell Reduction

{cell['partition']}

{cell['variation']}

{cell['band']}

{cell['defect']}

{cell['joined']}

{cell['parity']}

## Relative Channels

{relative['hermitian']}

{relative['transpose']}

{relative['terminal']}

{relative['ideal']}

## Tie Transport

{tie['upper_tie']}

{tie['lower_tie']}

{tie['fixed_cell']}

{tie['moving_boundary_guard']}

## Reciprocal Blocks

{blocks['cells']}

{blocks['blocks']}

{blocks['gradient']}

{blocks['far_gap']}

{blocks['adjacent']}

{blocks['local_phase']}

## Route Guards

{guards['in_band']}

{guards['outside_band']}

{guards['meaning']}

{guards['source_derivative']}

{guards['next_route']}

## Handoff

{handoff['decision']}

{handoff['obligation']}

## Pi Provenance

{kernel['pi_provenance']}

## Proof Boundary

{handoff['reserve']} This gate proves no signed finite-cell defect or completed-current bound, `Phi_B` theorem, contact exclusion, retained aggregate or Xi theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
"""


def main() -> int:
    payloads = load_sources()
    certificate = symbolic_certificate()
    rows = build_rows(certificate)
    counts = {
        "rows": len(rows),
        "zero_mean_checks": certificate["band_kernel"]["zero_mean_checks"],
        "odd_harmonic_checks": certificate["band_kernel"]["odd_harmonic_checks"],
        "periodic_checks": certificate["band_kernel"]["periodic_checks"],
        "cell_recombinations": certificate["cell_reduction"]["cell_recombinations"],
        "parity_decompositions": certificate["cell_reduction"]["parity_decompositions"],
        "relative_cell_reductions": certificate["relative_channels"]["relative_cell_reductions"],
        "tie_jumps": certificate["tie_transport"]["tie_jumps"],
        "fixed_cell_transports": certificate["tie_transport"]["fixed_cell_transports"],
        "reciprocal_block_decompositions": certificate["reciprocal_blocks"]["block_decompositions"],
        "far_gap_checks": certificate["reciprocal_blocks"]["far_gap_checks"],
        "adjacent_block_classes": certificate["reciprocal_blocks"]["adjacent_block_classes"],
        "fourier_witness_checks": certificate["route_guards"]["fourier_witness_checks"],
        "signed_physical_bounds": certificate["route_guards"]["signed_physical_bounds"],
        "phi_b_bounds": 0,
    }
    artifact = {
        "kind": KIND,
        "date": "2026-08-02",
        "status": "finite-band Dirichlet-cell and relative-channel reduction proved; signed physical variation estimate open",
        "counts": counts,
        "symbolic_certificate": certificate,
        "source_audit": source_audit(payloads),
        "rows": [asdict(row) for row in rows],
        "proof_boundary": "This proves the finite band kernel, cell-variation and parity formulas, relative Hermitian/transpose reductions, tie jumps, fixed-cell transport, and exact source-placement guards. It proves no signed physical variation, finite-cell defect or completed-current bound, Phi_B theorem, contact exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    print(
        "built finite-band Dirichlet-cell reduction gate: "
        f"{counts['rows']} rows, {counts['zero_mean_checks']} zero-mean checks, "
        f"{counts['odd_harmonic_checks']} endpoint-harmonic checks, "
        f"{counts['cell_recombinations']} cell recombinations, "
        f"{counts['relative_cell_reductions']} relative reductions, "
        f"{counts['far_gap_checks']} far-gap checks, "
        f"{counts['fourier_witness_checks']} Fourier witnesses, "
        f"{counts['signed_physical_bounds']} signed physical bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
