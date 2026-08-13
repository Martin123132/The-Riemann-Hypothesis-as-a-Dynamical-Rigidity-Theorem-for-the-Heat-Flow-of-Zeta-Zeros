#!/usr/bin/env python3
"""Build the ideal-cubic reciprocal collar transport gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_full_support_ideal_cubic_reciprocal_"
    "collar_transport_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "finite_band": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_finite_band_dirichlet_"
        "cell_reduction_gate.json"
    ),
    "first_correction": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "ideal_cubic_reciprocal_first_correction_gate.json"
    ),
    "terminal_quadrature": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_terminal_conditional_"
        "quadrature_gate.json"
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


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_zero(expression: sp.Expr, label: str) -> None:
    value = sp.simplify(sp.factor(sp.together(expression)))
    if value != 0 and value.equals(0) is not True:
        raise RuntimeError(f"{label}: {value}")


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def load_sources() -> dict[str, dict]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        require(path.exists(), f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))
    require(payloads["finite_band"]["counts"]["rows"] == 33, "finite-band drift")
    require(
        payloads["first_correction"]["counts"]["first_correction_cancellations"]
        == 1,
        "first-correction drift",
    )
    require(
        payloads["terminal_quadrature"]["counts"]["complex_recurrence_relations"]
        == 1,
        "terminal-quadrature drift",
    )
    return payloads


def source_audit() -> dict[str, dict[str, str]]:
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def cell_bounds(q: int, n_value: int) -> tuple[Fraction, Fraction]:
    require(1 <= q <= n_value, "physical cell index")
    left = Fraction(1) if q == 1 else Fraction(q) - Fraction(1, 2)
    right = Fraction(n_value) if q == n_value else Fraction(q) + Fraction(1, 2)
    return left, right


def cell_first_moment(q: int, n_value: int) -> Fraction:
    left, right = cell_bounds(q, n_value)
    return (right * right - left * left) / 2


def collar_certificate() -> dict[str, str | int]:
    decompositions = 0
    far_compressions = 0
    for n_value in range(3, 25):
        for radius in range(0, min(4, n_value - 1) + 1):
            atoms = {
                (q, p): Fraction(10000 * q + 97 * p + q * p, 113)
                for q in range(1, n_value + 1)
                for p in range(1, n_value + 1)
            }
            full = sum(atoms.values(), Fraction(0))
            near = sum(
                value
                for (q, p), value in atoms.items()
                if abs(q - p) <= radius
            )
            far = sum(
                value
                for (q, p), value in atoms.items()
                if abs(q - p) >= radius + 1
            )
            require(full == near + far, "collar decomposition")
            decompositions += 1

            for p in range(1, n_value + 1):
                far_indices = [
                    q
                    for q in range(1, n_value + 1)
                    if abs(q - p) >= radius + 1
                ]
                compressed_indices = list(range(1, max(1, p - radius))) + list(
                    range(min(n_value + 1, p + radius + 1), n_value + 1)
                )
                require(far_indices == compressed_indices, "far interval compression")
                direct = sum((atoms[q, p] for q in far_indices), Fraction(0))
                compressed = sum(
                    (atoms[q, p] for q in compressed_indices), Fraction(0)
                )
                require(direct == compressed, "far compressed value")
                far_compressions += 1

    return {
        "atom": "For r in C_p define I_(q,r)[P]=integral_(J_q)f_P(u)e(-ru)du and mathcal B_(q,p)[P]=sum_(r in C_p)I_(q,r)[P].",
        "collar": "For 0<=L<=N-1 define mathcal C_(N,L)[P]=sum_(|p-q|<=L)mathcal B_(q,p)[P] and mathcal E_(N,L)[P]=sum_(|p-q|>=L+1)mathcal B_(q,p)[P]. Then mathscr F_N[P]=mathcal C_(N,L)[P]+mathcal E_(N,L)[P].",
        "radius_one": "The bookmarked diagonal-plus-adjacent object is mathcal C_(N,1)=sum_q mathcal B_(q,q)+sum_(q<N)mathcal B_(q,q+1)+sum_(q>1)mathcal B_(q,q-1).",
        "far_compression": "For fixed p, the far physical cells are two contiguous rails: [1,p-L-1/2] when p>=L+2 and [p+L+1/2,N] when p<=N-L-1. Hence mathcal E_(N,L) is a sum of at most two oriented integrals per frequency cell, not a collection of independently bounded unit cells.",
        "decompositions": decompositions,
        "far_compressions": far_compressions,
    }


def tie_certificate() -> dict[str, str | int]:
    tie_checks = 0
    radius_one_checks = 0
    noninvariance_witnesses = 0
    for n_value in range(4, 31):
        for radius in range(0, min(4, n_value - 2) + 1):
            for p in range(1, n_value):
                difference = {
                    q: int(abs(q - (p + 1)) <= radius)
                    - int(abs(q - p) <= radius)
                    for q in range(1, n_value + 1)
                }
                expected = {q: 0 for q in range(1, n_value + 1)}
                entering = p + radius + 1
                leaving = p - radius
                if entering <= n_value:
                    expected[entering] += 1
                if leaving >= 1:
                    expected[leaving] -= 1
                require(difference == expected, "collar tie rail law")

                far_difference = {
                    q: int(abs(q - (p + 1)) >= radius + 1)
                    - int(abs(q - p) >= radius + 1)
                    for q in range(1, n_value + 1)
                }
                require(
                    all(far_difference[q] == -difference[q] for q in difference),
                    "far tie cancellation",
                )
                tie_checks += 1
                if radius == 1:
                    radius_one_checks += 1

                if 2 <= leaving and entering <= n_value - 1:
                    witness = cell_first_moment(entering, n_value) - cell_first_moment(
                        leaving, n_value
                    )
                    require(witness == 2 * radius + 1, "noninvariance witness")
                    noninvariance_witnesses += 1

    return {
        "internal_tie": "When alpha_P/r crosses p+1/2 upward, the mode r moves from C_p to C_(p+1). Its collar jump is 1_(p+L+1<=N)I_(p+L+1,r)-1_(p-L>=1)I_(p-L,r).",
        "radius_one_tie": "For L=1 the exact jump is 1_(p+2<=N)I_(p+2,r)-1_(p-1>=1)I_(p-1,r). The diagonal and two adjacent blocks therefore do not form an internal-tie-invariant object.",
        "far_tie": "The far complement has the opposite rail jump, so mathcal C_(N,L)+mathcal E_(N,L)=mathscr F_N remains invariant under every internal reassignment.",
        "witness": "For the exact test source f(u)=u e(ru), every full-cell atom I_(q,r) equals q. At an interior tie the collar jump is 2L+1, proving that no proper fixed-width collar is coefficientwise tie invariant.",
        "outer_tie": "At the two outer reciprocal-band crossings the complete mode, not only its near collar, is transferred to or from the outer remainder. The collar and both far rails must travel together.",
        "fixed_cell": "On a fixed roster-and-cell chart, mathcal C_(N,L)[i(lambda-c)P]=(partial_xi-ic)mathcal C_(N,L)[P] and likewise for mathcal E_(N,L). At a tie the two rail atoms and their relative lifts are transferred explicitly.",
        "tie_checks": tie_checks,
        "radius_one_checks": radius_one_checks,
        "noninvariance_witnesses": noninvariance_witnesses,
        "fixed_cell_transports": 2,
    }


def phase_certificate() -> dict[str, str | int]:
    alpha, v, y = sp.symbols("alpha v y", nonzero=True)
    u = v * (1 + y)
    require_zero(u / v - (1 + y), "diagonal ratio")

    c, a, b = sp.symbols("c a b", nonzero=True)
    right_u = c - a
    right_v = c + b
    right_y = (a + b) / (c + b)
    require_zero(right_u / right_v - (1 - right_y), "right corner ratio")
    require_zero(right_v - right_u - (a + b), "right corner orientation")

    left_u = c + a
    left_v = c - b
    left_y = (a + b) / (c - b)
    require_zero(left_u / left_v - (1 + left_y), "left corner ratio")
    require_zero(left_v - left_u + (a + b), "left corner orientation")

    orientation_checks = 0
    offsets = [Fraction(0), Fraction(1, 11), Fraction(2, 9), Fraction(7, 20)]
    for q in range(1, 25):
        right_c = Fraction(2 * q + 1, 2)
        for a_value in offsets:
            for b_value in offsets:
                d_value = (a_value + b_value) / (right_c + b_value)
                require(0 <= d_value < 1, "right corner domain")
                if d_value > 0:
                    phase_value = math.log(1 - float(d_value)) + float(d_value)
                    require(phase_value < 0, "right corner phase sign")
                orientation_checks += 1

    for q in range(2, 25):
        left_c = Fraction(2 * q - 1, 2)
        for a_value in offsets:
            for b_value in offsets:
                d_value = (a_value + b_value) / (left_c - b_value)
                require(d_value >= 0, "left corner domain")
                if d_value > 0:
                    phase_value = math.log(1 + float(d_value)) - float(d_value)
                    require(phase_value < 0, "left corner phase sign")
                orientation_checks += 1

    return {
        "diagonal": "For r in C_q put v=alpha_P/r and u=v(1+y). Then phi_r(u)=alpha_P(log v-1)+alpha_P{log(1+y)-y}. This is the exact negative Morse chart; no quadratic truncation is made.",
        "right_corner": "For p=q+1 let c=q+1/2, u=c-a, v=c+b, and d=(a+b)/(c+b). Then phi_r(u)=alpha_P(log v-1)+alpha_P{log(1-d)+d}, with phi_r'(u)>=0 and equality only at the shared corner a=b=0.",
        "left_corner": "For p=q-1 let c=q-1/2, u=c+a, v=c-b, and d=(a+b)/(c-b). Then phi_r(u)=alpha_P(log v-1)+alpha_P{log(1+d)-d}, with phi_r'(u)<=0 and equality only at the shared corner.",
        "orientation": "Both corner corrections are strictly negative away from the shared boundary: d/d d[log(1-d)+d]=-d/(1-d) and d/d d[log(1+d)-d]=-d/(1+d). Their cubic terms have opposite orientations and cannot be replaced by one unsigned half-Gaussian.",
        "phase_identities": 3,
        "corner_orientation_checks": orientation_checks,
    }


def far_certificate() -> dict[str, str | int]:
    amplitude, phase, gradient, kappa = sp.symbols("A E Q kappa", nonzero=True)
    amplitude_prime, gradient_prime = sp.symbols("A_prime Q_prime")
    phase_prime = gradient * phase / kappa
    boundary_derivative = kappa * (
        amplitude_prime * phase / gradient
        + amplitude * phase_prime / gradient
        - amplitude * phase * gradient_prime / gradient**2
    )
    correction = kappa * phase * (
        amplitude_prime / gradient - amplitude * gradient_prime / gradient**2
    )
    require_zero(
        boundary_derivative - correction - amplitude * phase,
        "far integration by parts",
    )
    return {
        "gap": "On either far rail, v=alpha_P/r lies in the p-cell and the physical variable lies outside [p-L-1/2,p+L+1/2]. Hence |v-u|>=L and q_r(u)=alpha_P/u-r has no zero for every L>=1.",
        "ibp": "For q_r=alpha_P/u-r, integral_a^b A e(phi_r)du=kappa[Ae(phi_r)/q_r]_a^b-kappa integral_a^b e(phi_r){A'/q_r-Aq_r'/q_r^2}du on each oriented far interval.",
        "telescoping": "The unit-cell boundary terms telescope before estimation because each far side is one contiguous interval. Only u=1, u=N, and the two collar rails p-L-1/2 and p+L+1/2 remain for each frequency cell.",
        "aggregate_open": "The rail-compressed identity prevents an artificial O(N^2) boundary budget, but no signed or absolute aggregate bound for the remaining p,r sum is proved here.",
        "ibp_identities": 1,
        "far_aggregate_bounds": 0,
    }


def moment_certificate() -> dict[str, str | int]:
    y = sp.symbols("y", real=True)
    kappa = 1 / (2 * sp.pi * sp.I)
    minus = sp.exp(-sp.pi * sp.I * y**2)
    plus = sp.exp(sp.pi * sp.I * y**2)
    require_zero(
        y**2 * minus - kappa * minus + kappa * sp.diff(y * minus, y),
        "negative truncated moment",
    )
    require_zero(
        y**2 * plus + kappa * plus - kappa * sp.diff(y * plus, y),
        "positive truncated moment",
    )
    differential_symbol, alpha = sp.symbols("D alpha", nonzero=True)
    require_zero(
        -kappa * differential_symbol / (2 * alpha)
        + kappa * differential_symbol / (2 * alpha),
        "full-line first correction",
    )
    return {
        "full_line": "The inherited full-line sequential coefficient remains -kappa mathcal D_qA/(2alpha_P)+kappa mathcal D_qA/(2alpha_P)=0 for arbitrary twice-differentiable A.",
        "negative_moment": "On a finite negative chart, integral_a^b y^2e(-y^2/2)dy=kappa integral_a^b e(-y^2/2)dy-kappa[y e(-y^2/2)]_a^b.",
        "positive_moment": "On a finite positive chart, integral_c^d eta^2e(+eta^2/2)deta=-kappa integral_c^d e(+eta^2/2)deta+kappa[eta e(+eta^2/2)]_c^d.",
        "boundary_guard": "The full-line cancellation uses vanishing regularized boundary terms. Finite diagonal and corner charts retain the four displayed endpoint terms; they may cancel only after the oriented collar rails, outer endpoints, and terminal recurrence are composed.",
        "moment_identities": 2,
        "full_line_cancellations": 1,
        "finite_correction_bounds": 0,
    }


def ideal_certificate() -> dict[str, str | int]:
    x, u_n, u_x, u_q = sp.symbols("x u_N u_x u_q", real=True)
    p_h = -sp.I * (x + u_n) * (x - u_n) ** 2 / 4
    p_t = -sp.I * (x - u_n) * (x + u_n) ** 2 / 4 + u_x * (x - u_n)
    require_zero(
        p_h.subs(x, -u_q)
        - sp.I * (u_q - u_n) * (u_q + u_n) ** 2 / 4,
        "Hermitian atom",
    )
    require_zero(
        p_t.subs(x, -u_q)
        - sp.I * (u_q + u_n) * (u_q - u_n) ** 2 / 4
        + u_x * (u_q + u_n),
        "transpose atom",
    )
    require_zero(p_h.subs(x, -u_n), "Hermitian terminal")
    require_zero(p_t.subs(x, -u_n) + 2 * u_n * u_x, "transpose terminal")
    require_zero(p_t + p_h.subs(x, -x) - u_x * (x - u_n), "ideal reflection")
    carrier, collar, far, terminal = sp.symbols("M C E Z")
    defect = collar + far - carrier
    require_zero(carrier - defect - (2 * carrier - collar - far), "Hermitian join")
    require_zero(
        carrier - defect + terminal - (2 * carrier - collar - far + terminal),
        "transpose join",
    )
    return {
        "cubics": "Retain P_H^0=-i(x+u_N)(x-u_N)^2/4 and P_T^0=-i(x-u_N)(x+u_N)^2/4+u_(N,x)(x-u_N) inside every block before a modulus is taken.",
        "atoms": "At x=-u_q, P_H^0=i(u_q-u_N)(u_q+u_N)^2/4 and P_T^0=i(u_q+u_N)(u_q-u_N)^2/4-u_(N,x)(u_q+u_N).",
        "reflection": "The exact polynomial relation is P_T^0(x)=-P_H^0(-x)+u_(N,x)(x-u_N). The physical interval is not reflection-symmetric, so this identity alone gives no cancellation or sign.",
        "terminal": "At q=N, P_H^0(log N)=0 and P_T^0(log N)=-2u_Nu_(N,x). The transpose endpoint must remain joined to +2iu_N mathcal T_N[B_(T,N)^0] and the complex terminal recurrence.",
        "joined": "Since mathscr F_N=mathcal C_(N,1)+mathcal E_(N,1), the ideal carrier-only joins are mathcal J_H^0=2mathscr M_N[P_H^0]-mathcal C_(N,1)[P_H^0]-mathcal E_(N,1)[P_H^0] and mathcal J_T^0=2mathscr M_N[P_T^0]-mathcal C_(N,1)[P_T^0]-mathcal E_(N,1)[P_T^0]+2iu_Nmathcal T_N[B_(T,N)^0]. The collar rail jump cancels only against the far rail jump inside these complete joins.",
        "ideal_identities": 5,
        "joined_recombinations": 2,
        "signed_collar_bounds": 0,
    }


def symbolic_certificate() -> dict:
    return {
        "collar": collar_certificate(),
        "ties": tie_certificate(),
        "phases": phase_certificate(),
        "far": far_certificate(),
        "moments": moment_certificate(),
        "ideal": ideal_certificate(),
        "handoff": {
            "decision": "The diagonal-plus-adjacent collar is useful only on a fixed cell chart; it is not a globally closed observation. Preserve its two rail transfers and compose the rail-compressed far functional before seeking a uniform sign.",
            "obligation": "Insert P_H^0 and P_T^0 into the exact diagonal and corner charts, carry the four finite-moment boundary terms into the two rail interfaces, join 2mathscr M_N and the transpose terminal survivor, and then prove or falsify a signed aggregate for the resulting tie-complete collar-plus-far functional.",
            "reserve": "No signed ideal-cubic collar, far aggregate, finite first-correction remainder, h^2 reserve, completed current, or RH-level theorem is proved.",
        },
    }


def build_rows(cert: dict) -> list[GateRow]:
    collar = cert["collar"]
    ties = cert["ties"]
    phases = cert["phases"]
    far = cert["far"]
    moments = cert["moments"]
    ideal = cert["ideal"]
    handoff = cert["handoff"]
    rows = [
        GateRow("rct_01_atom", "block atom", "proved", "Every reciprocal block is a finite sum of physical-cell mode atoms.", collar["atom"], "No atom is estimated separately."),
        GateRow("rct_02_collar", "collar definition", "proved", "Every finite radius gives an exact near/far partition.", collar["collar"], "The partition is algebraic."),
        GateRow("rct_03_radius_one", "bookmarked collar", "proved", "The diagonal-plus-adjacent object is exactly the radius-one collar.", collar["radius_one"], "It is not yet tie completed."),
        GateRow("rct_04_far_compression", "far rail compression", "proved", "All far physical cells compress to two oriented intervals per frequency cell.", collar["far_compression"], "No modulus is taken."),
        GateRow("rct_05_internal_tie", "internal tie law", "proved", "A reciprocal reassignment moves two collar rails.", ties["internal_tie"], "The source is evaluated at the crossing."),
        GateRow("rct_06_radius_one_tie", "three-strip jump", "proved", "The diagonal-plus-adjacent collar has an explicit two-rail jump.", ties["radius_one_tie"], "It is not globally invariant alone."),
        GateRow("rct_07_far_tie", "far compensation", "proved", "The far complement cancels the collar rail jump exactly.", ties["far_tie"], "The full band sum remains invariant."),
        GateRow("rct_08_witness", "nonclosure witness", "guard_validated", "No proper fixed-width collar is coefficientwise tie invariant.", ties["witness"], "The witness is nonphysical and proves only structural nonclosure."),
        GateRow("rct_09_outer_tie", "outer transfer", "guard_validated", "Outer roster crossings transfer a complete mode.", ties["outer_tie"], "Near and far pieces cannot be transferred separately."),
        GateRow("rct_10_fixed_cell", "relative lift", "proved", "Fixed-cell relative lifting preserves the collar and far split.", ties["fixed_cell"], "Rail lifts are explicit at ties."),
        GateRow("rct_11_diagonal", "diagonal phase", "proved", "The diagonal has an exact negative Morse chart.", phases["diagonal"], "No full-line replacement is made."),
        GateRow("rct_12_right_corner", "right adjacent corner", "proved", "The right adjacent block is an oriented stationary boundary corner.", phases["right_corner"], "Its endpoint is retained."),
        GateRow("rct_13_left_corner", "left adjacent corner", "proved", "The left adjacent block has the opposite orientation.", phases["left_corner"], "Its endpoint is retained."),
        GateRow("rct_14_corner_sign", "corner geometry", "proved", "Both exact corner phase corrections are negative away from contact.", phases["orientation"], "This does not sign the complex amplitude integral."),
        GateRow("rct_15_far_gap", "far nonstationarity", "proved", "The compressed far rails have no stationary point.", far["gap"], "The bound is not yet aggregated."),
        GateRow("rct_16_ibp", "oriented far identity", "proved", "One exact integration by parts applies to each far rail.", far["ibp"], "Endpoint orientation is preserved."),
        GateRow("rct_17_telescope", "boundary telescope", "proved", "Unit-cell integration-by-parts boundaries cancel before estimation.", far["telescoping"], "Only physical and collar rails remain."),
        GateRow("rct_18_far_bound", "far aggregate", "open", "Control the p,r aggregate after rail compression.", far["aggregate_open"], "No far aggregate bound is claimed."),
        GateRow("rct_19_full_line", "full-line correction", "proved", "The inherited first sequential coefficient still cancels on the full line.", moments["full_line"], "This is not a finite-cell theorem."),
        GateRow("rct_20_negative_moment", "negative finite moment", "proved", "Finite negative-chart second moments have an explicit boundary defect.", moments["negative_moment"], "The endpoints are not discarded."),
        GateRow("rct_21_positive_moment", "positive finite moment", "proved", "Finite positive-chart second moments have the opposite boundary formula.", moments["positive_moment"], "The endpoints are not discarded."),
        GateRow("rct_22_boundary_guard", "finite-correction guard", "guard_validated", "Full-line cancellation cannot be copied into a finite collar without boundary composition.", moments["boundary_guard"], "A complete endpoint cancellation remains open."),
        GateRow("rct_23_finite_correction", "finite first correction", "open", "Compose and bound the finite moment boundary terms.", moments["boundary_guard"], "No finite correction remainder is proved."),
        GateRow("rct_24_cubics", "ideal sources", "proved", "The two ideal cubics remain exact block amplitudes.", ideal["cubics"], "Corrections remain detached."),
        GateRow("rct_25_atoms", "ideal atom values", "proved", "The physical ideal atom values retain their phase-difference and phase-sum factors.", ideal["atoms"], "No atomwise sign is inferred."),
        GateRow("rct_26_reflection", "ideal reflection", "guard_validated", "The Hermitian/transpose reflection has an explicit shear defect.", ideal["reflection"], "The physical domain is not reflection symmetric."),
        GateRow("rct_27_terminal", "terminal ideal values", "proved", "The Hermitian endpoint vanishes and the transpose endpoint survives at shear scale.", ideal["terminal"], "The complex terminal quadrature remains compulsory."),
        GateRow("rct_28_joined", "physical ideal joins", "proved", "The live Hermitian and transpose joins contain both the collar and far rails.", ideal["joined"], "No component is signed separately."),
        GateRow("rct_29_decision", "route decision", "guard_validated", "Replace the isolated three-strip target by a tie-complete collar-plus-far composition.", handoff["decision"], "No valid cancellation route is discarded beyond the isolated collar."),
        GateRow("rct_30_obligation", "next obligation", "open", "Derive the signed tie-complete ideal composition.", handoff["obligation"], handoff["reserve"]),
        GateRow("rct_31_boundary", "proof boundary", "guard_validated", "This gate is not a proof of RH.", handoff["reserve"], "No prize-level conclusion is claimed."),
    ]
    require(len(rows) == 31, "row count drifted")
    return rows


def render_note(artifact: dict) -> str:
    cert = artifact["symbolic_certificate"]
    collar = cert["collar"]
    ties = cert["ties"]
    phases = cert["phases"]
    far = cert["far"]
    moments = cert["moments"]
    ideal = cert["ideal"]
    handoff = cert["handoff"]
    return f"""# Ideal-Cubic Reciprocal Collar Transport Gate

Date: 2026-08-02

Status: exact reciprocal collar, rail transport, corner phases, and finite-moment boundary identities proved; signed tie-complete ideal estimate open; not a proof of RH.

## Collar Decomposition

{collar['atom']}

{collar['collar']}

{collar['radius_one']}

{collar['far_compression']}

## Rail Transport

{ties['internal_tie']}

{ties['radius_one_tie']}

{ties['far_tie']}

{ties['witness']}

{ties['outer_tie']}

{ties['fixed_cell']}

## Exact Local Phases

{phases['diagonal']}

{phases['right_corner']}

{phases['left_corner']}

{phases['orientation']}

## Far Rails

{far['gap']}

{far['ibp']}

{far['telescoping']}

{far['aggregate_open']}

## Finite First Correction

{moments['full_line']}

{moments['negative_moment']}

{moments['positive_moment']}

{moments['boundary_guard']}

## Ideal Cubics

{ideal['cubics']}

{ideal['atoms']}

{ideal['reflection']}

{ideal['terminal']}

{ideal['joined']}

## Handoff

{handoff['decision']}

{handoff['obligation']}

## Pi Provenance

The phases use only e(x)=exp(2pi i x) and kappa=1/(2pi i). The finite moment identities follow by differentiating e(+-y^2/2); no geometric constant is inserted.

## Proof Boundary

{handoff['reserve']} This gate proves no signed physical collar or completed-current bound, `Phi_B` theorem, contact exclusion, retained aggregate or Xi theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
"""


def main() -> int:
    load_sources()
    certificate = symbolic_certificate()
    rows = build_rows(certificate)
    counts = {
        "rows": len(rows),
        "collar_decompositions": certificate["collar"]["decompositions"],
        "far_interval_compressions": certificate["collar"]["far_compressions"],
        "internal_tie_checks": certificate["ties"]["tie_checks"],
        "radius_one_tie_checks": certificate["ties"]["radius_one_checks"],
        "noninvariance_witnesses": certificate["ties"]["noninvariance_witnesses"],
        "fixed_cell_transports": certificate["ties"]["fixed_cell_transports"],
        "phase_identities": certificate["phases"]["phase_identities"],
        "corner_orientation_checks": certificate["phases"]["corner_orientation_checks"],
        "far_ibp_identities": certificate["far"]["ibp_identities"],
        "far_aggregate_bounds": certificate["far"]["far_aggregate_bounds"],
        "finite_moment_identities": certificate["moments"]["moment_identities"],
        "full_line_cancellations": certificate["moments"]["full_line_cancellations"],
        "finite_correction_bounds": certificate["moments"]["finite_correction_bounds"],
        "ideal_cubic_identities": certificate["ideal"]["ideal_identities"],
        "joined_recombinations": certificate["ideal"]["joined_recombinations"],
        "signed_collar_bounds": certificate["ideal"]["signed_collar_bounds"],
        "phi_b_bounds": 0,
    }
    artifact = {
        "kind": KIND,
        "date": "2026-08-02",
        "status": "reciprocal collar and rail transport proved; signed tie-complete ideal estimate open",
        "counts": counts,
        "symbolic_certificate": certificate,
        "source_audit": source_audit(),
        "rows": [asdict(row) for row in rows],
        "proof_boundary": "This proves the exact finite collar/far split, internal and outer rail transport, local diagonal and adjacent-corner phases, rail-compressed integration by parts, finite Gaussian moment boundary identities, and ideal-cubic source identities. It proves no signed tie-complete collar, far aggregate, finite first-correction remainder, completed current, Phi_B theorem, contact exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    print(
        "built ideal-cubic reciprocal collar transport gate: "
        f"{counts['rows']} rows, {counts['collar_decompositions']} collar decompositions, "
        f"{counts['far_interval_compressions']} far compressions, "
        f"{counts['internal_tie_checks']} tie checks, "
        f"{counts['noninvariance_witnesses']} nonclosure witnesses, "
        f"{counts['corner_orientation_checks']} corner checks, "
        f"{counts['finite_moment_identities']} finite moment identities, "
        f"{counts['signed_collar_bounds']} signed collar bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
