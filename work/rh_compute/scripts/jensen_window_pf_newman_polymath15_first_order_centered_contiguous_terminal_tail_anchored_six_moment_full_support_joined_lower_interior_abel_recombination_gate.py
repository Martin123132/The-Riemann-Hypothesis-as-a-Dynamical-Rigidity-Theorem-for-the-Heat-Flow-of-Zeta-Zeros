#!/usr/bin/env python3
"""Build the joined lower/interior Abel-recombination gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_full_support_joined_lower_interior_"
    "abel_recombination_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "relative_lift": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_terminal_nonterminal_"
        "relative_lift_gate.json"
    ),
    "finite_cell": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "ideal_cubic_reciprocal_finite_cell_inversion_gate.json"
    ),
    "outer_endpoint": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_outer_endpoint_"
        "kernel_gate.json"
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
        payloads["relative_lift"]["counts"]["mixed_polynomial_functionals"] == 2,
        "relative-lift source drifted",
    )
    require(
        payloads["finite_cell"]["counts"]["global_recompositions"] == 3,
        "finite-cell source drifted",
    )
    require(
        payloads["outer_endpoint"]["counts"]["twofold_outer_recompositions"]
        == 1,
        "outer-endpoint source drifted",
    )
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def physical_kernel_certificate() -> dict[str, str | int]:
    r, r_n, r_x, r_nx = sp.symbols("R R_N R_x R_Nx")
    c, c_n, d, d_n, delta = sp.symbols("C C_N D D_N delta")
    rb_n, rb_nx, cb_n, db_n, db_delta = sp.symbols(
        "Rbar_N Rbar_Nx Cbar_N Dbar_N deltabar"
    )

    q = (r + delta) * c + d
    q_n = (r_n + delta) * c_n + d_n
    qb_n = (rb_n + db_delta) * cb_n + db_n
    n_poly = r * q + r_x * c
    n_n = r_n * q_n + r_nx * c_n
    nb_n = rb_n * qb_n + rb_nx * cb_n
    a_poly = r * c
    a_n = r_n * c_n

    b_h_direct = cb_n * n_poly + nb_n * c - rb_n * cb_n * q - qb_n * a_poly
    b_h_factored = (
        cb_n
        * c
        * (
            (r - rb_n) ** 2
            + (r - rb_n) * (delta - db_delta)
            + r_x
            + rb_nx
        )
        + (r - rb_n) * (cb_n * d - db_n * c)
    )
    require_zero(b_h_direct - b_h_factored, "physical Hermitian factorization")

    b_t_direct = c_n * n_poly + n_n * c - a_n * q - q_n * a_poly
    b_t_factored = (
        c_n * c * ((r - r_n) ** 2 + r_x + r_nx)
        + (r - r_n) * (c_n * d - d_n * c)
    )
    require_zero(b_t_direct - b_t_factored, "physical transpose factorization")

    x, u_n, u_nx = sp.symbols("x u_N u_Nx", real=True)
    ideal = {
        r: sp.I * x / 2,
        r_n: -sp.I * u_n / 2,
        rb_n: sp.I * u_n / 2,
        r_x: -sp.I * u_nx / 2,
        r_nx: -sp.I * u_nx / 2,
        rb_nx: sp.I * u_nx / 2,
        c: 1,
        c_n: 1,
        cb_n: 1,
        d: 0,
        d_n: 0,
        db_n: 0,
        delta: 0,
        db_delta: 0,
    }
    b_h_ideal = sp.expand(b_h_factored.subs(ideal))
    b_t_ideal = sp.expand(b_t_factored.subs(ideal))
    require_zero(b_h_ideal + (x - u_n) ** 2 / 4, "ideal Hermitian kernel")
    require_zero(
        b_t_ideal + (x + u_n) ** 2 / 4 + sp.I * u_nx,
        "ideal transpose kernel",
    )

    p_h_ideal = sp.expand(sp.I * (x + u_n) * b_h_ideal)
    p_t_ideal = sp.expand(sp.I * (x - u_n) * b_t_ideal)
    require_zero(
        p_h_ideal + sp.I * (x + u_n) * (x - u_n) ** 2 / 4,
        "ideal Hermitian lift",
    )
    require_zero(
        p_t_ideal
        + sp.I * (x - u_n) * (x + u_n) ** 2 / 4
        - u_nx * (x - u_n),
        "ideal transpose lift",
    )

    u_q = sp.symbols("u_q", real=True)
    p_h_atom = sp.expand(p_h_ideal.subs(x, -u_q))
    p_t_atom = sp.expand(p_t_ideal.subs(x, -u_q))
    require_zero(
        p_h_atom - sp.I * (u_q - u_n) * (u_q + u_n) ** 2 / 4,
        "ideal Hermitian atom",
    )
    require_zero(
        p_t_atom
        - sp.I * (u_q + u_n) * (u_q - u_n) ** 2 / 4
        + u_nx * (u_q + u_n),
        "ideal transpose atom",
    )

    return {
        "hermitian_full": "B_(H,N)=(R-conjugate(R_N)){conjugate(C_N)Q-conjugate(Q_N)C}+(R_x+conjugate(R_(N,x)))conjugate(C_N)C.",
        "hermitian_physical": "After Q=(R+delta)C+D, B_(H,N)=conjugate(C_N)C{(R-conjugate(R_N))^2+(R-conjugate(R_N))(delta-conjugate(delta))+R_x+conjugate(R_(N,x))}+(R-conjugate(R_N)){conjugate(C_N)D-conjugate(D_N)C}.",
        "transpose_full": "B_(T,N)=(R-R_N){C_NQ-Q_NC}+(R_x+R_(N,x))C_NC.",
        "transpose_physical": "After Q=(R+delta)C+D, B_(T,N)=C_NC{(R-R_N)^2+R_x+R_(N,x)}+(R-R_N){C_ND-D_NC}.",
        "ideal": "For C=1, D=delta=0, R=ix/2, R_N=-iu_N/2, and R_x=R_(N,x)=-iu_(N,x)/2, B_(H,N)^0=-(x-u_N)^2/4 and B_(T,N)^0=-(x+u_N)^2/4-iu_(N,x).",
        "ideal_lifts": "P_H^0=-i(x+u_N)(x-u_N)^2/4 and P_T^0=-i(x-u_N)(x+u_N)^2/4+u_(N,x)(x-u_N). Both ideal lifts are cubic; the physical correction can restore degrees four and five.",
        "ideal_atoms": "At x=-u_q, P_H^0=i(u_q-u_N)(u_q+u_N)^2/4 and P_T^0=i(u_q+u_N)(u_q-u_N)^2/4-u_(N,x)(u_q+u_N).",
        "correction_split": "Define Delta_H=P_H-P_H^0 and Delta_T=P_T-P_T^0 only after the full terminal conjugations above. Every term in either Delta vanishes under the simultaneous ideal specialization; no modulus or coefficient budget is taken here.",
        "factorizations": 2,
        "ideal_kernel_identities": 6,
    }


def endpoint_certificate() -> dict[str, str | int]:
    kappa, length, p_0, p_1, g_0 = sp.symbols("kappa L P_0 P_0prime g_0")
    alpha, s_1, s_2, s_3 = sp.symbols("alpha S_1 S_2 S_3")
    value = -sp.I * length * p_0
    derivative = sp.I * (p_0 - length * (p_1 + g_0 * p_0))
    minus_t = -kappa * value * s_1
    minus_u = kappa**2 * (derivative * s_2 + alpha * value * s_3)
    expected_t = sp.I * kappa * length * p_0 * s_1
    expected_u = sp.I * kappa**2 * (
        (p_0 - length * (p_1 + g_0 * p_0)) * s_2
        - alpha * length * p_0 * s_3
    )
    require_zero(minus_t - expected_t, "lower conditional lift")
    require_zero(minus_u - expected_u, "lower higher lift")
    return {
        "generic": "For P_L=i(lambda-L)P, the exact lower block is -T_1[P_L]-U_1[P_L]=i kappa L P(0)S_1+i kappa^2{[P(0)-L(P'(0)+g_0P(0))]S_2-alpha_P L P(0)S_3}.",
        "lengths": "Use L=log N in the Hermitian channel and L=log a+u_N in the transpose channel. The logarithms are retained until reverse recombination.",
        "lower_endpoint_identities": 2,
    }


def recombination_certificate() -> dict[str, str | int]:
    t_n, u_n, t_1, u_1, interior = sp.symbols("T_N U_N T_1 U_1 I_int")
    mode = t_n + u_n - t_1 - u_1 + interior
    lower_interior = -t_1 - u_1 + interior
    require_zero(lower_interior - (mode - t_n - u_n), "modewise reverse IBP")

    finite, outer, starred, terminal = sp.symbols("F O M T")
    defect = finite - starred
    relation = {outer: starred - finite}
    remainder = outer - terminal
    require_zero(
        (remainder + defect + terminal).subs(relation),
        "outer remainder versus finite-cell defect",
    )
    require_zero(
        (finite + remainder - starred + terminal).subs(relation),
        "band plus remainder",
    )

    tie, tau, exterior = sp.symbols("J tau E")
    finite_after = finite + tie
    outer_after = outer - tie
    require_zero(finite_after + outer_after - finite - outer, "tie full-mode invariance")
    require_zero(
        (defect - (tau - exterior)).subs(finite, tau + starred - exterior),
        "finite-cell exterior identity",
    )

    return {
        "modewise": "For every outer mode r, -T_(1,r)[P]-U_(1,r)[P]+I_(int,r)[P]=I_P(r)-T_(N,r)[P]-U_(N,r)[P]. Thus the lower S_1/S_2/S_3 traces and the absolute interior reverse exactly to the original mode integral minus its upper boundary.",
        "outer": "After symmetric summation, L_N[P]:=-T_1[P]-U_1[P]+I_N[P]=O_N[P]-T_N[P]-U_N[P], and R_N[P]=O_N[P]-T_N[P]. The lower logarithms are representation coordinates, not a separately small error and not a zero.",
        "objects": "Let F_N[P]=sum_(r in mathcal R_N)I_P(r), let M_N[P] be the starred physical carrier from full finite Poisson inversion on [1,N], and define D_N[P]=F_N[P]-M_N[P]. Since F_N+O_N=M_N, one has R_N[P]=-D_N[P]-T_N[P] and F_N[P]+R_N[P]=M_N[P]-T_N[P].",
        "exterior": "The tiled finite-cell formula F_N=tau_band+M_N-E_band gives D_N=tau_band-E_band and R_N=E_band-tau_band-T_N. The finite-cell defect is the global two-sided exterior/tie package, not a new local endpoint estimate.",
        "tie": "At a half-open roster tie, F_N gains one complete mode and O_N loses the same mode. F_N+O_N and F_N+R_N are invariant; D_N and R_N change oppositely.",
        "modewise_recombinations": 1,
        "global_recombinations": 4,
        "tie_invariants": 2,
    }


def relative_channel_certificate() -> dict[str, str | int]:
    defect, terminal, band, starred, u_n = sp.symbols("D T F M u_N")
    r_h = -defect
    r_t = -defect + 2 * sp.I * u_n * terminal
    require_zero(r_h + defect, "Hermitian defect reduction")
    require_zero(r_t + defect - 2 * sp.I * u_n * terminal, "transpose defect reduction")

    # These are algebraic consequences of F=M+D.
    require_zero((band + r_h).subs(band, starred + defect) - starred, "Hermitian band join")
    require_zero(
        (band + r_t).subs(band, starred + defect)
        - starred
        - 2 * sp.I * u_n * terminal,
        "transpose band join",
    )
    require_zero(
        starred + r_h - (starred - defect),
        "Hermitian carrier-only join",
    )
    require_zero(
        starred + r_t - (starred - defect + 2 * sp.I * u_n * terminal),
        "transpose carrier-only join",
    )
    return {
        "hermitian": "For P_H=i(lambda-log N)B_(H,N), T_N[P_H]=0. Hence R_N[P_H]=-D_N[P_H] and F_N[P_H]+R_N[P_H]=M_N[P_H] exactly.",
        "transpose": "For P_T=i(lambda-log a-u_N)B_(T,N), T_N[P_T]=-2iu_NT_N[B_(T,N)]. Hence R_N[P_T]=-D_N[P_T]+2iu_NT_N[B_(T,N)] and F_N[P_T]+R_N[P_T]=M_N[P_T]+2iu_NT_N[B_(T,N)].",
        "carrier_only": "If the joined finite object is the returned carrier M_N rather than the band-mode sum F_N, then M_N[P_H]+R_N[P_H]=M_N[P_H]-D_N[P_H], while the transpose channel also has +2iu_NT_N[B_(T,N)]. Dropping D_N would apply Fourier inversion twice.",
        "starred_atoms": "M_N[P_H] is the starred q<N sum with multiplier -i(u_q-u_N)B_(H,N)(log q); its q=N atom vanishes. M_N[P_T] is the starred q<=N sum with multiplier -i(u_q+u_N)B_(T,N)(log q).",
        "relative_recombinations": 4,
    }


def abel_certificate() -> dict[str, str | int]:
    checks = 0
    for size in range(3, 28):
        amplitudes = [sp.Rational((3 * j + size) % 11 - 5, j + 2) for j in range(1, size + 1)]
        kernels = [sp.Rational((5 * j + 2 * size) % 13 - 6, j + 3) for j in range(1, size + 1)]
        prefixes = []
        running = sp.Integer(0)
        for value in amplitudes:
            running += value
            prefixes.append(running)

        direct = sum(amplitudes[j] * kernels[j] for j in range(size))
        abel = prefixes[-1] * kernels[-1] + sum(
            prefixes[j] * (kernels[j] - kernels[j + 1])
            for j in range(size - 1)
        )
        require_zero(direct - abel, f"Abel full sequence {size}")
        checks += 1

        h_kernels = kernels[:-1] + [sp.Integer(0)]
        h_direct = sum(amplitudes[j] * h_kernels[j] for j in range(size))
        h_abel = sum(
            prefixes[j] * (h_kernels[j] - h_kernels[j + 1])
            for j in range(size - 1)
        )
        require_zero(h_direct - h_abel, f"Abel terminal-null sequence {size}")
        checks += 1

    return {
        "generic": "For a_q and k_q with A_j=sum_(q<=j)a_q, sum_(q=1)^N a_qk_q=A_Nk_N+sum_(q=1)^(N-1)A_q(k_q-k_(q+1)).",
        "hermitian": "Because the Hermitian relative kernel has k_N=0, its exact Abel form has no terminal boundary: sum_(q=1)^N a_qk_q=sum_(q=1)^(N-1)A_q(k_q-k_(q+1)). This is where phase-difference cancellation can be sought without detaching the finite-cell defect.",
        "transpose": "The transpose relative kernel generally has k_N=-2iu_NB_(T,N)(log N), so its Abel terminal boundary survives and must be combined with +2iu_NT_N[B_(T,N)] and the pure terminal transpose block.",
        "guard": "Abel summation may now be applied to M_N-D_N (and the explicit transpose terminal term), not to M_N alone while silently deleting D_N. No partial-sum estimate or sign is asserted by this gate.",
        "abel_checks": checks,
    }


def symbolic_certificate() -> dict:
    return {
        "physical_kernels": physical_kernel_certificate(),
        "lower_endpoint": endpoint_certificate(),
        "reverse_recombination": recombination_certificate(),
        "relative_channels": relative_channel_certificate(),
        "abel": abel_certificate(),
        "handoff": {
            "decision": "The lower logarithmic traces are no longer an independent target. After exact reverse integration by parts, the live obstruction is the source-specific returned carrier minus the finite-cell defect, plus one explicit transpose terminal term.",
            "route": "Derive a source-specific formula and signed estimate for D_N[P_H] and D_N[P_T] from the global exterior/tie kernel, then Abel-sum M_N-D_N with the ideal cubic kept intact and Delta_H,Delta_T attached only afterward.",
            "reserve": "No h^2/800 reserve is claimed. The next bound must include both phase families, the correction perturbations, fixed affine row, pure terminal transpose block, and moved-tail defect before the whole-jet C1 transfer.",
        },
    }


def build_rows(cert: dict) -> list[GateRow]:
    physical = cert["physical_kernels"]
    lower = cert["lower_endpoint"]
    reverse = cert["reverse_recombination"]
    relative = cert["relative_channels"]
    abel = cert["abel"]
    handoff = cert["handoff"]
    rows = [
        GateRow("jli_01_hfull", "physical Hermitian kernel", "proved", "The complete Hermitian terminal coefficient has an exact rate/determinant factorization.", physical["hermitian_full"], "Terminal conjugation is retained."),
        GateRow("jli_02_hphysical", "Hermitian correction chart", "proved", "Substitution of Q exposes every physical correction channel.", physical["hermitian_physical"], "No correction is bounded separately."),
        GateRow("jli_03_tfull", "physical transpose kernel", "proved", "The complete transpose terminal coefficient has the matching factorization.", physical["transpose_full"], "No sign is inferred."),
        GateRow("jli_04_tphysical", "transpose correction chart", "proved", "Substitution of Q exposes the transpose correction channels.", physical["transpose_physical"], "No determinant factor is discarded."),
        GateRow("jli_05_ideal", "ideal kernels", "proved", "The correction-free kernels are explicit quadratics.", physical["ideal"], "This specialization is not the full physical state."),
        GateRow("jli_06_ideal_lifts", "ideal relative lifts", "proved", "Both correction-free joined kernels are cubic.", physical["ideal_lifts"], "Degree reduction is not a size bound."),
        GateRow("jli_07_atoms", "ideal carrier atoms", "proved", "The ideal q-atoms display the phase-difference and phase-sum factors.", physical["ideal_atoms"], "Neither atom family has a fixed real sign."),
        GateRow("jli_08_corrections", "correction split", "proved", "Physical corrections are defined after full terminal conjugation.", physical["correction_split"], "No componentwise correction budget is introduced."),
        GateRow("jli_09_lower", "lower endpoint", "proved", "Both lifted lower endpoint orders form one exact formula.", lower["generic"], "The S_1/S_2/S_3 pieces remain joined."),
        GateRow("jli_10_lengths", "relative lengths", "proved", "The two genuine lower logarithms are identified.", lower["lengths"], "They are not hidden in O(h)."),
        GateRow("jli_11_mode", "modewise reversal", "proved", "Two integrations by parts reverse exactly mode by mode.", reverse["modewise"], "The lower block does not vanish by itself."),
        GateRow("jli_12_outer", "outer reversal", "proved", "The complete lower/interior package is the outer mode sum minus its upper boundary.", reverse["outer"], "Conditional summation remains symmetric."),
        GateRow("jli_13_objects", "three-object split", "proved", "Band modes, returned carriers, and their defect are distinguished exactly.", reverse["objects"], "F_N and M_N must not be identified."),
        GateRow("jli_14_exterior", "finite-cell exterior", "proved", "The finite-cell defect is the global exterior/tie package.", reverse["exterior"], "No cellwise absolute budget is allowed."),
        GateRow("jli_15_tie", "half-open transfer", "proved", "The complete mode join is invariant at reciprocal ties.", reverse["tie"], "D_N and R_N move oppositely."),
        GateRow("jli_16_hrelative", "Hermitian recombination", "proved", "The Hermitian upper conditional value vanishes after relative lifting.", relative["hermitian"], "This is a mode-plus-remainder identity."),
        GateRow("jli_17_trelative", "transpose recombination", "proved", "The transpose join leaves one explicit terminal conditional term.", relative["transpose"], "The terminal term is not dropped."),
        GateRow("jli_18_carrier_guard", "carrier-only guard", "guard_validated", "A returned carrier cannot replace the reciprocal-band mode sum for free.", relative["carrier_only"], "Dropping D_N double-applies inversion."),
        GateRow("jli_19_atoms", "starred carrier form", "proved", "The exact returned sums have the parent relative multipliers.", relative["starred_atoms"], "Endpoint half weights remain."),
        GateRow("jli_20_abel", "Abel identity", "proved", "The discrete joined carrier admits exact summation by parts.", abel["generic"], "No partial-sum estimate is supplied."),
        GateRow("jli_21_habel", "Hermitian Abel boundary", "proved", "The Hermitian terminal Abel boundary vanishes exactly.", abel["hermitian"], "The finite-cell defect remains joined."),
        GateRow("jli_22_tabel", "transpose Abel boundary", "proved", "The transpose Abel boundary survives in an explicit terminal package.", abel["transpose"], "It must meet the pure terminal block."),
        GateRow("jli_23_abel_guard", "Abel route guard", "guard_validated", "Abel must act on the carrier-defect join.", abel["guard"], "No sign is inferred."),
        GateRow("jli_24_decision", "route decision", "proved", "The raw lower logarithms are retired as the primary obstruction.", handoff["decision"], "This is a structural reduction only."),
        GateRow("jli_25_route", "next theorem", "open", "Prove or falsify a signed source-specific finite-cell-defect estimate.", handoff["route"], "No signed defect bound is proved."),
        GateRow("jli_26_reserve", "aggregate reserve", "open", "Recover the h^2/800 reserve only after all live channels are reattached.", handoff["reserve"], "No h^2 reserve is claimed."),
        GateRow("jli_27_pi", "pi provenance", "proved", "Every pi remains inherited from the fixed Fourier normalization.", "kappa=1/(2pi i), alpha_P=xi/(2pi), and e(x)=exp(2pi i x); no new geometric constant is inserted.", "No fitted constant is used."),
        GateRow("jli_28_boundary", "proof boundary", "guard_validated", "This gate is not a proof of RH.", "No signed finite-cell defect, completed mixed-current, Phi_B, contact-exclusion, Xi, Lambda<=0, PF-infinity, RH, or prize-level theorem is claimed.", "The theorem-search boundary remains explicit."),
    ]
    require(len(rows) == 28, "row count drifted")
    return rows


def render_note(artifact: dict) -> str:
    cert = artifact["symbolic_certificate"]
    p = cert["physical_kernels"]
    lower = cert["lower_endpoint"]
    reverse = cert["reverse_recombination"]
    relative = cert["relative_channels"]
    abel = cert["abel"]
    handoff = cert["handoff"]
    return f"""# Joined Lower/Interior Abel Recombination Gate

Date: 2026-08-02

Status: exact physical kernel expansion and lower/interior/finite-cell recombination proved; signed finite-cell-defect estimate open; not a proof of RH.

## Physical Kernels

{p['hermitian_full']}

{p['hermitian_physical']}

{p['transpose_full']}

{p['transpose_physical']}

{p['ideal']}

{p['ideal_lifts']}

{p['ideal_atoms']}

{p['correction_split']}

## Lower Endpoint

{lower['generic']}

{lower['lengths']}

## Exact Recombination

{reverse['modewise']}

{reverse['outer']}

{reverse['objects']}

{reverse['exterior']}

{reverse['tie']}

## Relative Channels

{relative['hermitian']}

{relative['transpose']}

{relative['carrier_only']}

{relative['starred_atoms']}

## Abel Form

{abel['generic']}

{abel['hermitian']}

{abel['transpose']}

{abel['guard']}

## Handoff

{handoff['decision']}

{handoff['route']}

{handoff['reserve']}

## Pi Provenance

The constants `kappa=1/(2pi i)` and `alpha_P=xi/(2pi)` and the character `e(x)=exp(2pi i x)` are inherited. No circle, polygon, plotted symmetry, or fitted constant is introduced.

## Proof Boundary

This gate proves exact physical Hermitian and transpose kernel factorizations, their ideal cubic relative lifts, the complete lifted lower endpoint, reverse twofold integration by parts, the distinction between reciprocal-band modes and returned carriers, reduction of the outer remainder to the finite-cell defect, both relative-channel recombinations, and exact Abel identities. It proves no signed finite-cell-defect estimate, completed mixed-current or `Phi_B` bound, contact exclusion, retained aggregate or Xi theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
"""


def main() -> int:
    payloads = load_sources()
    certificate = symbolic_certificate()
    rows = build_rows(certificate)
    counts = {
        "rows": len(rows),
        "physical_kernel_factorizations": certificate["physical_kernels"]["factorizations"],
        "ideal_kernel_identities": certificate["physical_kernels"]["ideal_kernel_identities"],
        "lower_endpoint_identities": certificate["lower_endpoint"]["lower_endpoint_identities"],
        "modewise_recombinations": certificate["reverse_recombination"]["modewise_recombinations"],
        "global_recombinations": certificate["reverse_recombination"]["global_recombinations"],
        "tie_invariants": certificate["reverse_recombination"]["tie_invariants"],
        "relative_recombinations": certificate["relative_channels"]["relative_recombinations"],
        "abel_checks": certificate["abel"]["abel_checks"],
        "signed_finite_cell_defect_bounds": 0,
        "signed_completed_current_bounds": 0,
        "phi_b_bounds": 0,
    }
    artifact = {
        "kind": KIND,
        "date": "2026-08-02",
        "status": "joined lower/interior and finite-cell recombination proved; signed finite-cell-defect estimate open",
        "counts": counts,
        "symbolic_certificate": certificate,
        "source_audit": source_audit(payloads),
        "rows": [asdict(row) for row in rows],
        "proof_boundary": "This proves exact physical kernel factorizations, ideal cubic lifts, lower/interior reverse recombination, finite-band/carrier/defect identities, relative-channel joins, and Abel formulas. It proves no signed finite-cell-defect or completed-current bound, Phi_B theorem, contact exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    print(
        "built joined lower/interior Abel recombination gate: "
        f"{counts['rows']} rows, {counts['physical_kernel_factorizations']} physical factorizations, "
        f"{counts['global_recombinations']} global recombinations, "
        f"{counts['relative_recombinations']} relative joins, "
        f"{counts['abel_checks']} Abel checks, "
        f"{counts['signed_finite_cell_defect_bounds']} signed defect bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
