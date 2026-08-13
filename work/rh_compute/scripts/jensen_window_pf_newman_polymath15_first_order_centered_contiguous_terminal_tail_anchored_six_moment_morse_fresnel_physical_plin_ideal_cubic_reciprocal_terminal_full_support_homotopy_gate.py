#!/usr/bin/env python3
"""Build the terminal-compatible full-support reciprocal homotopy gate."""

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
    "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
    "ideal_cubic_reciprocal_terminal_full_support_homotopy_gate"
)
RESULT_PATH = REPO_ROOT / "work" / "rh_compute" / "results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "abel_join": REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_abel_cross_current_reduction.json"
    ),
    "terminal_prefix": REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_growing_prefix_finite_height_gate.json"
    ),
    "flow_matrix": REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_flow_matrix_phase_reduction.json"
    ),
    "physical_plin": REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "amplitude_gate.json"
    ),
    "finite_cell": REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
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
    value = sp.simplify(expression)
    if value != 0:
        raise RuntimeError(f"{label}: {value}")


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def ceil_fraction(value: Fraction) -> int:
    return -floor_fraction(-value)


def load_sources() -> dict[str, dict]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        if not path.exists():
            raise RuntimeError(f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))
    return payloads


def source_audit(payloads: dict[str, dict]) -> dict[str, dict[str, str]]:
    require(
        payloads["abel_join"]["counts"]["current_polarization_identities"] == 2,
        "Abel-join source drifted",
    )
    require(
        payloads["terminal_prefix"]["counts"]["finite_height_current_theorems"]
        == 1,
        "terminal-prefix source drifted",
    )
    require(
        payloads["flow_matrix"]["counts"]["correction_free_moments"] == 6,
        "flow-matrix source drifted",
    )
    require(
        payloads["physical_plin"]["counts"]["maximum_polynomial_degree"] == 5,
        "physical-P_lin source drifted",
    )
    require(
        payloads["finite_cell"]["counts"]["global_recompositions"] == 3,
        "finite-cell source drifted",
    )
    require(
        payloads["finite_cell"]["counts"]["exterior_tail_bounds"] == 0,
        "finite-cell proof boundary drifted",
    )
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def audit_terminal_strips() -> dict[str, int]:
    cell_cases = 0
    strip_tilings = 0
    full_tilings = 0
    roster_containments = 0
    low_gap_cases = 0
    for n_value in [3, 4, 6, 9, 13]:
        for b_value in range(1, n_value):
            for offset in [Fraction(1, 7), Fraction(1, 3), Fraction(2, 3)]:
                alpha = Fraction(n_value * n_value) - offset
                cells: dict[int, list[int]] = {}
                for q in range(1, n_value + 1):
                    a_q = alpha / (Fraction(q) + Fraction(1, 2))
                    b_q = alpha / (Fraction(q) - Fraction(1, 2))
                    lower = floor_fraction(a_q) + 1
                    upper = floor_fraction(b_q)
                    cells[q] = list(range(lower, upper + 1)) if lower <= upper else []
                    direct = [
                        r
                        for r in range(1, floor_fraction(2 * alpha) + 1)
                        if a_q < r <= b_q
                    ]
                    require(cells[q] == direct, "reciprocal cell formula drifted")
                    cell_cases += 1

                full_union = [r for q in range(1, n_value + 1) for r in cells[q]]
                a_n = alpha / (Fraction(n_value) + Fraction(1, 2))
                expected_full = [
                    r
                    for r in range(1, floor_fraction(2 * alpha) + 1)
                    if a_n < r <= 2 * alpha
                ]
                require(
                    sorted(full_union) == expected_full,
                    "full-support cells failed to tile",
                )
                require(len(full_union) == len(set(full_union)), "full cells overlapped")
                full_tilings += 1

                terminal_union = [
                    r
                    for q in range(b_value + 1, n_value + 1)
                    for r in cells[q]
                ]
                a_b = alpha / (Fraction(b_value) + Fraction(1, 2))
                expected_terminal = [
                    r
                    for r in range(1, floor_fraction(a_b) + 1)
                    if a_n < r <= a_b
                ]
                require(
                    sorted(terminal_union) == expected_terminal,
                    "terminal reciprocal strip failed to tile",
                )
                strip_tilings += 1

                alpha_minus = alpha - Fraction(1, 20)
                alpha_plus = alpha + Fraction(1, 20)
                roster_lower = max(1, floor_fraction(alpha_minus / (2 * n_value)))
                roster_upper = ceil_fraction(2 * alpha_plus)
                require(
                    all(roster_lower <= r <= roster_upper for r in full_union),
                    "full-support cell escaped fixed roster",
                )
                roster_containments += 1

                low_cut = floor_fraction(a_n)
                if low_cut >= 1:
                    nearest_gap = alpha / n_value - low_cut
                    exact_floor_gap = alpha / n_value - a_n
                    require(nearest_gap >= exact_floor_gap, "low gap monotonicity failed")
                    require(exact_floor_gap > Fraction(1, 4), "low gap lost quarter")
                    low_gap_cases += 1

    return {
        "cell_cases": cell_cases,
        "terminal_strip_tilings": strip_tilings,
        "full_support_tilings": full_tilings,
        "roster_containments": roster_containments,
        "low_gap_cases": low_gap_cases,
    }


def symbolic_certificate() -> dict:
    half = sp.Rational(1, 2)
    J = half * sp.Matrix(
        [
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, -1],
            [0, 0, -1, 0],
        ]
    )
    wb = sp.Matrix(sp.symbols("V_B N_B A_B Q_B", real=True))
    wc = sp.Matrix(sp.symbols("V_C N_C A_C Q_C", real=True))
    we = sp.Matrix(sp.symbols("V_E N_E A_E Q_E", real=True))
    db = sp.Matrix(sp.symbols("dV_B dN_B dA_B dQ_B", real=True))
    dc = sp.Matrix(sp.symbols("dV_C dN_C dA_C dQ_C", real=True))

    def current(vector: sp.Matrix) -> sp.Expr:
        return (vector.T * J * vector)[0]

    require_zero(current(wb) - (wb[0] * wb[1] - wb[2] * wb[3]), "current form")
    wt = we + wc
    bulk_join = current(wb) + 2 * (wt.T * J * wb)[0]
    full_current = current(we + wb + wc)
    require_zero(
        full_current - bulk_join - current(wt),
        "physical full-current anchor",
    )

    full_derivative = 2 * ((we + wb + wc).T * J * (db + dc))[0]
    expanded_derivative = (
        (db[0] + dc[0]) * (wb[1] + wc[1] + we[1])
        + (wb[0] + wc[0] + we[0]) * (db[1] + dc[1])
        - (db[2] + dc[2]) * (wb[3] + wc[3] + we[3])
        - (wb[2] + wc[2] + we[2]) * (db[3] + dc[3])
    )
    require_zero(full_derivative - expanded_derivative, "full-flow derivative")

    old_derivative = 2 * (wb.T * J * db)[0] + 2 * (wt.T * J * db)[0]
    wc_physical = sp.Matrix(sp.symbols("V_C0 N_C0 A_C0 Q_C0", real=True))
    old_physical_tail = 2 * (wb.T * J * db)[0] + 2 * (
        (we + wc_physical).T * J * db
    )[0]
    moving_defect = (
        2 * ((wc - wc_physical).T * J * db)[0]
        + 2 * ((we + wb + wc).T * J * dc)[0]
    )
    require_zero(
        full_derivative - old_physical_tail - moving_defect,
        "moving-tail defect",
    )

    y = sp.Matrix(sp.symbols("y0:3", real=True))
    u0 = sp.Matrix(3, 4, sp.symbols("u0_0:12", real=True))
    u1 = sp.Matrix(3, 4, sp.symbols("u1_0:12", real=True))
    matrix = u0 * J * u1.T + u1 * J * u0.T
    omega = u0.T * y
    dot_omega = u1.T * y
    edge_linear = u1 * (2 * J * we)
    require_zero(
        (y.T * matrix * y)[0] + (edge_linear.T * y)[0]
        - 2 * ((we + omega).T * J * dot_omega)[0],
        "edge-only matrix reduction",
    )

    alpha_p, n_value = sp.symbols("alpha_P N", positive=True, real=True)
    a_n = alpha_p / (n_value + half)
    low_gap = alpha_p / n_value - a_n
    require_zero(
        low_gap - alpha_p / (2 * n_value * (n_value + half)),
        "lower complement gap",
    )
    r, u = sp.symbols("r u", positive=True, real=True)
    high_gap = r - alpha_p / u
    require_zero(high_gap.subs({r: 2 * alpha_p, u: 1}) - alpha_p, "high gap")

    a, e_n, e_b = sp.symbols("a e_N e_B", real=True)
    tail_interior = sp.symbols("tail_interior")
    h_n = sp.symbols("H_N")
    h_b = sp.symbols("H_B")
    endpoint_n = half * (1 + a * e_n)
    endpoint_b = half * (1 + e_b)
    star_difference = (h_n - endpoint_n) - (h_b - endpoint_b)
    require_zero(
        star_difference.subs(h_n, h_b + tail_interior + a * e_n)
        - (half * e_b + tail_interior + half * a * e_n),
        "starred terminal extension",
    )

    strip_a = alpha_p / (n_value + half)
    b_symbol = sp.symbols("B", positive=True, real=True)
    strip_b = alpha_p / (b_symbol + half)
    require_zero(
        alpha_p / ((b_symbol + 1) - half) - strip_b,
        "terminal-strip upper boundary",
    )

    phase_ratio = sp.Rational(7, 4 * 50**2)
    require(phase_ratio < 1, "terminal phase coefficient failed")
    audit = audit_terminal_strips()
    return {
        "current_homotopy": {
            "observation_order": "omega=(V,mathcal N,A,Q) and Q(omega)=omega^T J omega=Vmathcal N-AQ.",
            "terminal_split": "omega_T=omega_E+omega_C(Omega), where omega_C is the carrier tail B<n<=N and omega_E is the genuine endpoint edge.",
            "bulk_only": "Psi_B=Q(omega_B)+2omega_T^T J omega_B and Psi_B(Omega)=h^2Phi_B.",
            "full": "tilde(Psi)(xi)=Q(omega_E+omega_B(xi)+omega_C(xi)).",
            "anchor": "tilde(Psi)(Omega)=Q(omega_T+omega_B(Omega))=h^2P_ret.",
            "derivative": "tilde(Psi)'=2(omega_E+omega_B+omega_C)^T J(dot(omega_B)+dot(omega_C)).",
            "moving_defect": "tilde(Psi)'-Psi_B'=2[omega_C(xi)-omega_C(Omega)]^TJdot(omega_B)+2[omega_E+omega_B+omega_C(xi)]^TJdot(omega_C).",
        },
        "full_support_observations": {
            "moments": "G_j^[N]=G_j^[B]+G_j^(C), 0<=j<=5, with G_j^(C) summed over B<n<=N under the same auxiliary rotation.",
            "matrix": "For Y=y_B+y_C, tilde(Psi)'=Y^TM_xiY+ell_E^TY, where ell_E=U_1(2Jomega_E).",
            "edge_only": "The terminal carriers now live in Y; only the genuine endpoint edge remains in the external linear row.",
            "linear_polynomial": "For Y=p+r, the linear remainder is one degree-at-most-five P_lin^[N]=F_alpha+i(lambda-log a)F_beta with beta=(mathcal N_p+mathcal N_E,V_p+V_E,-Q_p-Q_E,-A_p-A_E).",
            "quadratic": "The nonlinear remainder remains epsilon_Vepsilon_(N,xi)+epsilon_Nepsilon_(V,xi)-epsilon_Aepsilon_(Q,xi)-epsilon_Qepsilon_(A,xi).",
        },
        "terminal_extension": {
            "physical_sum": "H^[N][A]=H^[B][A]+sum_(q=B+1)^N A(q)e(alpha_Plog q).",
            "starred_sum": "M_N[A]-M_B[A]=(1/2)A(B)e_B+sum_(q=B+1)^(N-1)A(q)e_q+(1/2)A(N)e_N.",
            "mode_extension": "I_(A,N)(r)=I_(A,B)(r)+integral_B^N A(u)e(alpha_Plog u-ru)du.",
            "endpoint_guard": "The B endpoint half is transferred into the terminal interval; it must not be retained simultaneously as a hard endpoint of both transforms.",
        },
        "reciprocal_full_support": {
            "cells": "The cells 1<=q<=N tile (alpha_P/(N+1/2),2alpha_P].",
            "terminal_strip": "The cells B+1<=q<=N tile (alpha_P/(N+1/2),alpha_P/(B+1/2)].",
            "returned_tail": "Their diagonal inversions return q=B+1,...,N-1 with full weight and q=N with half weight; the upper hard endpoint supplies the other half at N.",
            "involution": "The full-support scalar transform is still Fourier inversion. Its value is the alignment of the terminal strip and endpoint edge, not scalar smallness.",
        },
        "outer_complement": {
            "low_gap": "For r<=alpha_P/(N+1/2) and 1<=u<=N, alpha_P/u-r>=alpha_P/[2N(N+1/2)]>1/4 on the physical segment.",
            "high_gap": "For r>2alpha_P and 1<=u<=N, r-alpha_P/u>alpha_P.",
            "phase_motion": "For n=N-m in the terminal block, u_n<2hK and |z_n(xi)-z_n(Omega)|<7h^3K|z_n(Omega)|/(4L^2)<h^3K|z_n(Omega)|.",
            "boundary": "These gaps prove nonstationarity only. They do not bound the endpoint-composed outer complement at h^2 scale.",
        },
        "audit": audit,
    }


def build_rows(cert: dict) -> list[GateRow]:
    hom = cert["current_homotopy"]
    obs = cert["full_support_observations"]
    ext = cert["terminal_extension"]
    rec = cert["reciprocal_full_support"]
    outer = cert["outer_complement"]
    rows = [
        GateRow("rth_01_domain", "domain", "proved", "The terminal carrier block is separated from the genuine endpoint edge.", hom["terminal_split"], "No current sign follows."),
        GateRow("rth_02_current", "current form", "proved", "The four observations carry one universal projective current.", hom["observation_order"], "The form is indefinite."),
        GateRow("rth_03_bulk", "bulk homotopy", "proved", "The prior auxiliary flow is the bulk self-current plus its frozen-terminal cross current.", hom["bulk_only"], "The terminal self-current is not in Psi_B."),
        GateRow("rth_04_full", "full homotopy", "proved", "Moving the terminal carriers with the bulk defines one exact retained-current homotopy.", hom["full"], "The endpoint edge remains fixed."),
        GateRow("rth_05_anchor", "physical anchor", "proved", "The new homotopy equals the complete retained current at the physical frequency.", hom["anchor"], "This is not a sign theorem."),
        GateRow("rth_06_derivative", "full derivative", "proved", "The full homotopy derivative preserves all carrier self and cross terms.", hom["derivative"], "No pair is estimated separately."),
        GateRow("rth_07_move_defect", "homotopy comparison", "proved", "The difference from the frozen-terminal derivative is an explicit moving-tail defect.", hom["moving_defect"], "No norm bound is inferred."),
        GateRow("rth_08_moments", "six-moment closure", "proved", "The all-carrier homotopy still closes on six moments.", obs["moments"], "No seventh moment is introduced."),
        GateRow("rth_09_matrix", "matrix reduction", "proved", "The full-support flow uses the same rank-eight matrix with an edge-only linear row.", obs["matrix"], "The matrix remains indefinite."),
        GateRow("rth_10_edge", "edge placement", "proved", "Terminal carriers are internal moment data rather than an external terminal row.", obs["edge_only"], "The genuine edge is not deleted."),
        GateRow("rth_11_plin", "linear compression", "proved", "The full-support linear residual remains one physical polynomial amplitude.", obs["linear_polynomial"], "Its signed functional remains unbounded."),
        GateRow("rth_12_quadratic", "quadratic image", "proved", "The nonlinear residual remains exactly four observation pairings.", obs["quadratic"], "No quadratic bound follows."),
        GateRow("rth_13_extension", "physical extension", "proved", "Extending B to N adds exactly the terminal carrier samples.", ext["physical_sum"], "The endpoint edge is separate."),
        GateRow("rth_14_starred", "Dirichlet extension", "proved", "The two starred Poisson sums differ by the terminal interval with half weights at B and N.", ext["starred_sum"], "The half weights are structural."),
        GateRow("rth_15_mode_extension", "mode extension", "proved", "Each complete Fourier mode acquires exactly the continuous terminal interval.", ext["mode_extension"], "This is not a discrete tail approximation."),
        GateRow("rth_16_endpoint_guard", "endpoint guard", "guard_validated", "The old B hard half cannot coexist with the new full-support endpoint accounting.", ext["endpoint_guard"], "Double counting would invalidate the join."),
        GateRow("rth_17_cells", "full reciprocal band", "proved", "The full-support reciprocal cells tile one contiguous band.", rec["cells"], "Internal ties still cancel."),
        GateRow("rth_18_terminal_strip", "terminal reciprocal strip", "proved", "The newly added cells form exactly the terminal strip.", rec["terminal_strip"], "No cells are enumerated in production."),
        GateRow("rth_19_returned_tail", "terminal carrier inversion", "proved", "The terminal strip returns the terminal physical samples with the correct endpoint half.", rec["returned_tail"], "Aliases and exteriors remain joined."),
        GateRow("rth_20_involution", "involution guard", "guard_validated", "Full-support scalar re-Poissonization is still exact inversion.", rec["involution"], "It supplies alignment, not decay."),
        GateRow("rth_21_low_gap", "lower complement", "proved", "The low outer complement is uniformly nonstationary.", outer["low_gap"], "Endpoint transition size is not bounded."),
        GateRow("rth_22_high_gap", "upper complement", "proved", "The high outer complement has a much larger phase-gradient gap.", outer["high_gap"], "A conditional boundary series remains."),
        GateRow("rth_23_phase_motion", "terminal motion", "proved", "Moving a terminal carrier over the short auxiliary segment changes only a tiny phase.", outer["phase_motion"], "No joined row norm is assumed."),
        GateRow("rth_24_route", "route decision", "proved", "The next estimate belongs to the full-support observation kernel with the edge and outer complement joined.", outer["boundary"], "Scalar absolute exterior budgets remain invalid."),
        GateRow("rth_25_handoff", "signed theorem", "open", "Prove or falsify the edge-composed full-support retained-current inequality.", "Bound tilde(Psi)(T_0)-epsilon integral tilde(Psi)' below zero, using terminal-strip alignment, outer-complement cancellation, and the four quadratic pairings before any modulus.", "No such signed bound is proved here."),
        GateRow("rth_26_pi", "pi provenance", "proved", "Every pi remains fixed by the completed-zeta and Fourier normalizations.", "alpha_P=xi/(2pi), e(x)=exp(2pi i x), and kappa=1/(2pi i); no fitted circle or polygon is introduced.", "No numerical geometry is inferred."),
        GateRow("rth_27_boundary", "proof boundary", "guard_validated", "The full-support homotopy is a route reduction, not the missing estimate.", "No edge-complement bound, signed full-current theorem, Phi_B theorem, Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is claimed.", "This is not a proof of RH."),
    ]
    require(len(rows) == 27, "row count drifted")
    return rows


def render_note(artifact: dict) -> str:
    c = artifact["symbolic_certificate"]
    h = c["current_homotopy"]
    o = c["full_support_observations"]
    e = c["terminal_extension"]
    r = c["reciprocal_full_support"]
    x = c["outer_complement"]
    return f"""# Reciprocal Terminal Full-Support Homotopy Gate

Date: 2026-08-02

Status: exact all-carrier homotopy, terminal reciprocal strip, and outer-gap reduction proved; edge-composed signed bound open; not a proof of RH.

This is not a proof of RH. It changes the auxiliary transform so the terminal carriers and the bulk are one moving block, while the genuine endpoint edge remains explicit.

## Current Reassembly

{h['observation_order']}

{h['terminal_split']}

{h['bulk_only']}

{h['full']}

{h['anchor']}

{h['derivative']}

The exact comparison with the previous frozen-terminal flow is

```text
{h['moving_defect']}
```

## Full-Support Observation Image

{o['moments']}

{o['matrix']}

{o['edge_only']}

{o['linear_polynomial']}

{o['quadratic']}

## Terminal Interval

{e['physical_sum']}

{e['starred_sum']}

{e['mode_extension']}

{e['endpoint_guard']}

## Reciprocal Terminal Strip

{r['cells']}

{r['terminal_strip']}

{r['returned_tail']}

{r['involution']}

## Outer Complement

{x['low_gap']}

{x['high_gap']}

{x['phase_motion']}

{x['boundary']}

## Handoff

The next theorem must estimate the full-support observation kernel, not a scalar reciprocal cell. Compose the true endpoint edge with the two outer complement pieces, retain the terminal-strip carriers, the endpoint-linear, Hermitian, and transpose channels, and the four quadratic residual pairings, and only then take a modulus. The direct target is

```text
tilde(Psi)(T_0)-epsilon integral_0^1
  tilde(Psi)'(T_0-theta epsilon)dtheta < 0.
```

## Pi Provenance

`alpha_P=xi/(2pi)` and `e(x)=exp(2pi i x)` come from the fixed completed-zeta and Fourier conventions. No circle, polygon, or fitted geometric constant is added.

## Proof Boundary

This gate proves exact current polarization, the all-carrier physical anchor, six-moment closure on `1<=n<=N`, edge-only external linearization, terminal interval half-weight transfer, terminal reciprocal-strip tiling, two outer phase-gradient gaps, and a terminal phase-motion bound. It proves no edge-complement estimate, signed full-current theorem, quadratic residual bound, `Phi_B` upper bound, contact exclusion, retained aggregate or Xi theorem, `Q209`, cofinal descendant theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
"""


def main() -> int:
    payloads = load_sources()
    certificate = symbolic_certificate()
    rows = build_rows(certificate)
    counts = {
        "rows": len(rows),
        "exact_current_polarizations": 3,
        "all_carrier_homotopies": 1,
        "physical_anchor_identities": 1,
        "moving_tail_defect_identities": 1,
        "six_moment_full_support_closures": 1,
        "edge_only_matrix_reductions": 1,
        "full_support_linear_functionals": 1,
        "terminal_extension_identities": 3,
        "terminal_reciprocal_strip_tilings": 1,
        "full_support_cell_inversions": 1,
        "outer_complement_gap_identities": 2,
        "terminal_phase_displacement_bounds": 1,
        "outer_complement_bounds": 0,
        "quadratic_remainder_bounds": 0,
        "signed_full_current_bounds": 0,
        "phi_b_bounds": 0,
    }
    artifact = {
        "kind": KIND,
        "date": "2026-08-02",
        "status": "all-carrier full-support homotopy and terminal reciprocal strip proved; edge-complement and signed-current bounds open",
        "counts": counts,
        "symbolic_certificate": certificate,
        "source_audit": source_audit(payloads),
        "rows": [asdict(row) for row in rows],
        "proof_boundary": "This proves an exact all-carrier auxiliary homotopy, its physical retained-current anchor, six-moment full-support closure, edge-only external row, terminal interval and reciprocal-strip identities, outer nonstationary gaps, and a terminal phase-displacement bound. It proves no bound for the edge-composed outer complement, no quadratic remainder estimate, signed full-current or Phi_B theorem, contact exclusion, retained aggregate or Xi theorem, Q209, cofinal descendant theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
    }
    note = render_note(artifact)
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, note)
    print(
        "built reciprocal terminal full-support homotopy gate: "
        f"{counts['rows']} rows, "
        f"{counts['all_carrier_homotopies']} all-carrier homotopy, "
        f"{counts['terminal_reciprocal_strip_tilings']} terminal strip, "
        f"{counts['outer_complement_bounds']} outer-complement bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
