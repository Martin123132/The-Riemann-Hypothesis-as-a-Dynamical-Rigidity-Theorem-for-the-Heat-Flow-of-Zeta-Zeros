#!/usr/bin/env python3
"""Build the C1 terminal-phase winding and covariant two-channel gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_terminal_phase_winding_covariant_two_channel_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "q1_phase": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_five_moment_q1_saddle_phase_variance_reduction.json"
    ),
    "relative_lift": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_terminal_nonterminal_"
        "relative_lift_gate.json"
    ),
    "anchor": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_carrier_near_odd_harmonic_anchor_gate.json",
    "disk": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_reciprocal_stationary_disk_geometry_gate.json",
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
    expanded = sp.expand(expression)
    if expanded != 0:
        raise RuntimeError(f"{label}: {expanded}")


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
        require(path.is_file(), f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))
    return payloads


def audit_sources(payloads: dict[str, dict]) -> dict[str, dict[str, str]]:
    require(
        payloads["q1_phase"]["counts"]["physical_q1_parameter_laws"] == 4,
        "q1 phase source drifted",
    )
    require(
        payloads["q1_phase"]["counts"]["common_phase_laws"] == 1,
        "q1 common-phase source drifted",
    )
    require(
        payloads["relative_lift"]["counts"]["mixed_current_compressions"] == 2,
        "relative-lift source drifted",
    )
    require(
        payloads["anchor"]["summary"]["explicit_raw_anchors"] == 2,
        "anchor source drifted",
    )
    require(
        payloads["disk"]["summary"]["live_signed_stationary_targets"] == 1,
        "disk frontier drifted",
    )
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def phase_cell_certificate() -> dict[str, str | int]:
    ell = sp.symbols("L", positive=True, real=True)
    a_squared = sp.exp(ell) + sp.Rational(1, 32) / ell**2
    derivative = sp.exp(ell) - sp.Rational(1, 16) / ell**3
    require_zero(sp.diff(a_squared, ell) - derivative, "cutoff derivative")

    n, eps_left, eps_right, log_n, pi = sp.symbols(
        "N eps_left eps_right log_N pi", positive=True, real=True
    )
    omega_left = 2 * pi * n**2 - eps_left
    omega_right = 2 * pi * (n + 1) ** 2 - eps_right
    delta_omega = 2 * pi * (2 * n + 1) + eps_left - eps_right
    require_zero(omega_right - omega_left - delta_omega, "endpoint frequency")
    delta_phi = delta_omega * log_n
    require_zero(
        delta_phi
        - (2 * pi * (2 * n + 1) + eps_left - eps_right) * log_n,
        "endpoint phase",
    )

    # The physical epsilon bound is far below one already with x>1.
    require(Fraction(7, 8 * 50**2) < 1, "epsilon below one")
    # Using pi>3, log(2)>1/2, and pi<22/7 gives a rational full-turn check.
    full_turn_lower = Fraction(4 * 3 * 2, 2)
    full_turn_upper = Fraction(44, 7)
    require(full_turn_lower > full_turn_upper, "full-turn rational check")

    return {
        "cutoff_law": "On q=1, a(L)^2=exp(L)+1/(32L^2), T_0=2pi a(L)^2, Omega=T_0-epsilon, and 0<epsilon<7/(8L^2x), with x=4pi exp(L).",
        "cutoff_monotonicity": "The derivative of a(L)^2 is exp(L)-1/(16L^3)>0 for L>=50. Thus a is continuous, strictly increasing, and unbounded.",
        "complete_cells": "For every integer N>=ceil(a(50)), there are unique L_N and L_(N+1) with a(L_N)=N and a(L_(N+1))=N+1. This is a complete fixed-N cutoff cell.",
        "terminal_phase": "At the physical frequency xi=Omega, E_N=|E_N|exp(i phi_N) with phi_N(L)=Omega(L)log N.",
        "endpoint_identity": "Across a complete fixed-N cell, Delta phi_N=[2pi(2N+1)+epsilon(L_N)-epsilon(L_(N+1))]log N.",
        "epsilon_guard": "For L>=50, x>1 makes 0<epsilon<7/(8L^2x)<1. Hence Delta phi_N>[2pi(2N+1)-1]log N>4pi N log N.",
        "full_turn": "Since N>=2 and log N>=log 2>1/2, Delta phi_N>4pi N log N>2pi. The last comparison can be checked using pi>3 and pi<22/7; the actual physical N is vastly larger.",
        "circle_coverage": "A continuous real phase whose endpoint increment exceeds 2pi assumes a representative of every residue class modulo 2pi in the cell interior: choose k=floor((phi_left-theta)/(2pi))+1, so phi_left<theta+2pi k<=phi_left+2pi<phi_right. Therefore exp(i phi_N) covers the whole unit circle independently of the half-open endpoint convention.",
        "sector_guard": "No fixed proper phase sector contains the terminal phase on a complete fixed-N cell. In particular sin(phi_N) takes both signs and vanishes, so the raw negative anchors in Section 11.213 cannot have a cell-uniform physical projection by a fixed terminal sector.",
        "complete_cell_winding_theorems": 1,
        "circle_surjectivity_theorems": 1,
        "fixed_sector_obstructions": 1,
    }


def covariance_certificate() -> dict[str, str | int]:
    h, u, v, x, y, lam = sp.symbols("h u v x y lambda", real=True)
    matrix = sp.Matrix([[h + u, -v], [-v, h - u]]) / 2
    vector = sp.Matrix([x, y])
    quadratic = (vector.T * matrix * vector)[0]
    expected = (
        h * (x**2 + y**2) + u * (x**2 - y**2) - 2 * v * x * y
    ) / 2
    require_zero(quadratic - expected, "two-channel quadratic form")
    require_zero(sp.trace(matrix) - h, "matrix trace")
    require_zero(
        matrix.det() - (h**2 - u**2 - v**2) / 4,
        "matrix determinant",
    )
    characteristic = sp.expand((lam * sp.eye(2) - matrix).det())
    expected_characteristic = sp.expand(
        (lam - h / 2) ** 2 - (u**2 + v**2) / 4
    )
    require_zero(characteristic - expected_characteristic, "matrix spectrum")

    a_h, a_t, rho, s = sp.symbols("A_H A_T rho s", positive=True, real=True)
    h_anchor = -a_h * rho * s
    t_anchor_modulus = a_t * rho
    require_zero(
        h_anchor.subs(s, -1) + t_anchor_modulus - rho * (a_h + a_t),
        "anchor center obstruction",
    )

    return {
        "mixed_current": "Let h=Re(conjugate(tau_0)mathcal J_H^0), t=tau_0 mathcal J_T^0, and z=nu c_xi. The ideal joined mixed current is R_id(z)=(1/2){|z|^2 h+Re(z^2 t)}.",
        "matrix_form": "Writing z=x+iy and t=u+iv gives R_id=[x y]Q(h,t)[x y]^T with Q(h,t)=(1/2)[[h+u,-v],[-v,h-u]].",
        "spectral_form": "The trace and determinant are h and (h^2-|t|^2)/4, so the eigenvalues are (h-|t|)/2 and (h+|t|)/2.",
        "actual_phase_criterion": "For z=r exp(i gamma), R_id=(r^2/2){h+Re(exp(2i gamma)t)}. Negativity for the actual source phase is exactly h+Re(exp(2i gamma)t)<0.",
        "uniform_phase_criterion": "Negativity for every nonzero source phase is equivalent to negative definiteness of Q, hence to the single phase-invariant inequality h+|t|<0. Nonpositivity is equivalent to h+|t|<=0.",
        "anchor_pair": "For tau_0=rho exp(i(phi_N-pi/2)) and the raw centres mathcal J_H^0=-A_H, mathcal J_T^0=-A_T, one has h_0=-rho A_H sin(phi_N) and |t_0|=rho A_T.",
        "center_obstruction": "Because every complete cell contains phi_N=3pi/2 modulo 2pi, the centre-only phase-uniform margin h_0+|t_0| equals rho(A_H+A_T)>0 there. Static negative centres cannot prove a source-phase-uniform current theorem on a complete cell.",
        "surviving_target": "Estimate the complete tie-invariant Hermitian/transpose pair before projection. Either prove the actual-phase inequality h+Re(exp(2i gamma_phys)t) plus all correction rows is below the required negative reserve, or prove the stronger invariant inequality h+|t| plus those rows is below it. Residual phase correlation, not two independent anchor disks, is now the live arithmetic content.",
        "proof_guard": "The phase-uniform criterion is stronger than the theorem at one known physical source phase. Its failure for the bare centres does not disprove the complete current: the signed reciprocal residuals, physical corrections, fixed affine row, and moved-tail defect have not been estimated.",
        "covariant_matrix_identities": 1,
        "source_phase_uniform_criteria": 1,
        "center_only_obstructions": 1,
    }


def build_rows(phase: dict, covariance: dict) -> list[GateRow]:
    return [
        GateRow("tpw_01_q1_law", "physical phase law", "proved", "The q=1 cutoff and physical frequency have an exact one-sided epsilon law.", phase["cutoff_law"], "Only the fixed physical q=1 chart is used."),
        GateRow("tpw_02_monotone", "cutoff monotonicity", "proved", "The saddle cutoff is strictly increasing for L>=50.", phase["cutoff_monotonicity"], "No integer roster is frozen across a crossing."),
        GateRow("tpw_03_cells", "complete cells", "proved", "Every sufficiently large integer cutoff has a unique complete L-cell.", phase["complete_cells"], "The first partial cell at L=50 is excluded."),
        GateRow("tpw_04_phase", "terminal phase", "proved", "The physical terminal phase is Omega log N on a fixed-N cell.", phase["terminal_phase"], "This is the phase of E_N before source projection."),
        GateRow("tpw_05_endpoint", "endpoint identity", "proved", "The exact cell endpoint phase increment is explicit.", phase["endpoint_identity"], "No phase samples or fitted unwrapping are used."),
        GateRow("tpw_06_epsilon", "epsilon guard", "proved", "The one-sided physical correction cannot remove the macroscopic phase increment.", phase["epsilon_guard"], "The bound is deliberately coarse but uniform."),
        GateRow("tpw_07_turn", "full-turn theorem", "proved", "Every complete fixed-N cell gains more than one full terminal turn.", phase["full_turn"], "This is a terminal-phase theorem, not a current sign."),
        GateRow("tpw_08_circle", "circle coverage", "proved", "The terminal phase covers every direction on every complete cell.", phase["circle_coverage"], "Coverage follows from continuity and endpoint variation."),
        GateRow("tpw_09_sector", "sector obstruction", "proved", "A fixed proper terminal phase sector cannot be uniform on a complete cell.", phase["sector_guard"], "Adaptive phase localization is not ruled out."),
        GateRow("tpw_10_current", "two-channel current", "proved", "The ideal joined current has one Hermitian scalar and one transpose spin-two scalar.", covariance["mixed_current"], "Correction and affine rows remain outside this ideal block."),
        GateRow("tpw_11_matrix", "covariant matrix", "proved", "The source phase is represented by an exact real symmetric 2x2 quadratic form.", covariance["matrix_form"], "No positivity of that matrix is asserted."),
        GateRow("tpw_12_spectrum", "matrix spectrum", "proved", "The two invariant eigenvalues are determined by h and |t|.", covariance["spectral_form"], "This is algebraic compression only."),
        GateRow("tpw_13_actual", "actual-phase criterion", "proved", "The physical source phase has an exact scalar sign criterion.", covariance["actual_phase_criterion"], "The inequality itself remains unproved."),
        GateRow("tpw_14_uniform", "phase-uniform criterion", "proved", "A source-phase-uniform theorem is equivalent to h+|t|<0.", covariance["uniform_phase_criterion"], "This criterion is stronger than one fixed physical phase."),
        GateRow("tpw_15_centres", "centre-only obstruction", "proved", "The raw anchor centres cannot make the covariant matrix negative on a complete cell.", covariance["center_obstruction"], "This does not obstruct signed residual rotation."),
        GateRow("tpw_16_target", "covariant current target", "open", "The live target is a joined phase-correlated Hermitian/transpose current estimate.", covariance["surviving_target"], "No signed current, contact exclusion, Lambda<=0, or RH conclusion is claimed."),
    ]


def render_note(payload: dict) -> str:
    phase = payload["exact"]["phase_cell"]
    covariance = payload["exact"]["covariant_current"]
    lines = [
        "# Newman C1 Terminal-Phase Winding and Covariant Two-Channel Gate",
        "",
        "Date: 2026-08-04",
        "",
        "Status: complete fixed-N terminal-phase winding and the exact two-channel source-phase matrix proved; joined physical-current inequality open; not a proof of RH.",
        "",
        "## Complete-cell phase winding",
        "",
        phase["cutoff_law"],
        "",
        phase["cutoff_monotonicity"],
        "",
        "```text",
        "Delta phi_N",
        " =[2pi(2N+1)+epsilon(L_N)-epsilon(L_(N+1))]log N",
        " >4pi N log N",
        " >2pi.",
        "```",
        "",
        phase["circle_coverage"],
        "",
        phase["sector_guard"],
        "",
        "## Covariant two-channel form",
        "",
        covariance["mixed_current"],
        "",
        "```text",
        "Q(h,t)=(1/2)[[h+Re(t), -Im(t)],",
        "             [-Im(t), h-Re(t)]].",
        "",
        "eigenvalues(Q)=(h-|t|)/2, (h+|t|)/2.",
        "```",
        "",
        covariance["actual_phase_criterion"],
        "",
        covariance["uniform_phase_criterion"],
        "",
        covariance["center_obstruction"],
        "",
        covariance["proof_guard"],
        "",
        "## Gate rows",
        "",
        "| id | role | state | claim |",
        "|---|---|---|---|",
    ]
    for row in payload["rows"]:
        lines.append(
            f"| {row['id']} | {row['role']} | {row['readiness']} | {row['claim']} |"
        )
    lines.extend(
        [
            "",
            "## Next action",
            "",
            payload["next_action"],
            "",
            "## Pi provenance",
            "",
            payload["pi_provenance"],
            "",
            "## Proof boundary",
            "",
            payload["proof_boundary"],
            "",
            payload["success"],
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    sources = load_sources()
    source_audit = audit_sources(sources)
    phase = phase_cell_certificate()
    covariance = covariance_certificate()
    rows = build_rows(phase, covariance)
    require(len(rows) == 16, "row count")
    require(sum(row.readiness == "open" for row in rows) == 1, "open row count")

    success = (
        "built Newman C1 terminal-phase winding/covariant two-channel gate: "
        "16 rows, 0 issues, 1 full-cell winding theorem, 1 circle-surjectivity "
        "theorem, 1 phase-uniform criterion, 1 centre-only obstruction, "
        "1 live covariant current target"
    )
    payload = {
        "kind": KIND,
        "schema_version": 1,
        "date": "2026-08-04",
        "status": "complete fixed-N terminal-phase winding and exact covariant two-channel matrix proved; joined physical-current inequality open",
        "source_sha256": {key: item["sha256"] for key, item in source_audit.items()},
        "source_audit": source_audit,
        "exact": {
            "phase_cell": phase,
            "covariant_current": covariance,
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "complete_cell_winding_theorems": phase["complete_cell_winding_theorems"],
            "circle_surjectivity_theorems": phase["circle_surjectivity_theorems"],
            "fixed_sector_obstructions": phase["fixed_sector_obstructions"],
            "covariant_matrix_identities": covariance["covariant_matrix_identities"],
            "source_phase_uniform_criteria": covariance["source_phase_uniform_criteria"],
            "center_only_obstructions": covariance["center_only_obstructions"],
            "live_covariant_current_targets": 1,
        },
        "next_action": "Work with the complete physical pair h=Re(conjugate(tau_0)mathcal J_H) and t=tau_0 mathcal J_T, including the correction polynomials, fixed affine row, and moved-tail defect before projection. First derive the exact actual-source-phase scalar h+Re(exp(2i gamma_phys)t) in the contact-box normalization. Then test whether reciprocal stationary recomposition can prove its signed reserve. Use h+|t| only if a genuinely source-phase-uniform theorem is sought; do not estimate the two channels or their anchor disks independently.",
        "pi_provenance": "The factor 2pi in T_0=2pi a^2 and alpha_P=xi/(2pi) comes from the completed-zeta saddle and the Fourier character e(y)=exp(2pi i y). The threshold 2pi is one period of exp(i phi_N). The rational bounds pi>3 and pi<22/7 only certify a deliberately coarse full-turn inequality. No circle, polygon, or fitted plotting constant is inserted into the model.",
        "proof_boundary": "This gate proves complete-cell winding of the physical q=1 terminal phase, failure of any fixed proper terminal phase sector on such a cell, the exact Hermitian/transpose source-phase quadratic matrix, its invariant eigenvalues, and a centre-only obstruction to source-phase-uniform negativity. It proves no signed reciprocal residual estimate, actual-source-phase current inequality, correction-row bound, all-q transport, contact-box exclusion, Q209 shell, cofinal successor, Lambda<=0, RH, or prize-level conclusion.",
        "success": success,
    }

    atomic_write(RESULT_PATH, json.dumps(payload, indent=2) + "\n")
    atomic_write(NOTE_PATH, render_note(payload))
    print(success)


if __name__ == "__main__":
    main()
