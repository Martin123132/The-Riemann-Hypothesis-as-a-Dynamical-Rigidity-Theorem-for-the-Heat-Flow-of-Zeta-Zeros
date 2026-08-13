#!/usr/bin/env python3
"""Build the actual-source saddle-phase lock and centre-symmetry gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_saddle_phase_lock_center_symmetry_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "normalizer_lock": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_reciprocal_normalizer_phase_reinforcement_guard.json",
    "absolute_anchor": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_absolute_phase_anchor_reduction.json",
    "q1_phase": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_five_moment_q1_saddle_phase_variance_reduction.json"
    ),
    "raw_anchors": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_carrier_near_odd_harmonic_anchor_gate.json",
    "relative_lift": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_terminal_nonterminal_"
        "relative_lift_gate.json"
    ),
    "terminal_covariance": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_terminal_phase_winding_covariant_two_channel_gate.json",
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
    expanded = sp.expand_trig(sp.expand(expression))
    if sp.simplify(expanded) != 0:
        raise RuntimeError(f"{label}: {expanded}")


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(content, encoding="utf-8")
    os.replace(temporary, path)


def load_and_audit_sources() -> tuple[dict[str, dict], dict[str, dict[str, str]]]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        require(path.is_file(), f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))

    normalizer_text = json.dumps(payloads["normalizer_lock"], sort_keys=True)
    require("|psi_t|<15/(16x)<1/x" in normalizer_text, "normalizer phase bound drifted")
    require("Q(v)=exp(i*[v*log(v/(2*pi))-v-pi/4])" in normalizer_text, "Q phase drifted")
    require(payloads["absolute_anchor"]["constants"]["d_absolute"] == 2189, "d_1 bound drifted")
    require("eta=f_1/|f_1|" in payloads["absolute_anchor"]["exact"]["unit_anchor"], "unit anchor drifted")
    require(payloads["q1_phase"]["counts"]["common_phase_laws"] == 1, "q1 phase drifted")
    require(payloads["raw_anchors"]["summary"]["explicit_raw_anchors"] == 2, "raw anchors drifted")
    require(payloads["relative_lift"]["counts"]["mixed_current_compressions"] == 2, "mixed current drifted")
    require(payloads["terminal_covariance"]["summary"]["source_phase_uniform_criteria"] == 1, "covariance source drifted")

    audit = {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }
    return payloads, audit


def source_phase_certificate() -> dict[str, str | int]:
    return {
        "normalizer": "The common source normalizer is nu=eta/S_a=omega_c/(B_0T_0), where omega_c=(-1)^N eta and eta=phi(1+d_1)/|1+d_1|. At xi=Omega one has c_xi=1, |nu|=(B_0T_0)^(-1), and exp(2i gamma_phys)=eta^2; cutoff parity cancels from the transpose phase.",
        "continuous_phase": "Because M_t(s) is nonzero and |d_1|<1/2, choose continuous unwrapped beta_t=arg M_t(s), delta_1=arg(1+d_1), theta_eta=beta_t+delta_1, and chi_N=phi_N+theta_eta on each complete fixed-N q=1 cell.",
        "stationary_lock": "For Q(v)=exp(i[v log(v/(2pi))-v-pi/4]) and gamma_t=conjugate(M_t)/M_t=exp(-2i beta_t), Section 11.119 supplies the continuous psi=arg(Q(Omega)/gamma_t), with |psi|<1/x. Hence arg Q(Omega)+2beta_t=psi exactly on the chosen lift.",
        "source_normalizer_identifications": 1,
        "continuous_phase_lifts": 1,
    }


def reverse_turn_certificate() -> dict[str, str | int]:
    u = sp.symbols("u", positive=True, real=True)
    lower_log = u - u**2 / 2 + u**3 / 3 - u**4 / 4
    polynomial_upper = sp.cancel(
        ((1 + u) ** 2 * (1 - 2 * lower_log) - 1) / u**2
    )
    expected_upper = -2 - 2 * u / 3 + u**2 / 6 + u**3 / 3 + u**4 / 2
    require_zero(polynomial_upper - expected_upper, "alternating-log polynomial")

    # For 0<u<=1/2, the three positive terms are bounded by
    # u/12+u/12+u/16=11u/48, leaving a 7u/16 strict margin.
    require(
        Fraction(1, 12) + Fraction(1, 12) + Fraction(1, 16)
        == Fraction(11, 48),
        "positive tail budget",
    )
    require(Fraction(2, 3) - Fraction(11, 48) == Fraction(7, 16), "turn margin")

    # The physical endpoint perturbation is below 1947/N^2. The actual
    # cells have N>2^25, while N>742 already preserves the ideal margin.
    require(Fraction(1947, 2**25) < Fraction(21, 8), "physical margin audit")

    return {
        "lock_identity": "Let K_N=2pi N^2 and F_N(v)=v[1+log(K_N/v)]. The stationary lock gives the exact unwrapped identity 2chi_N=F_N(Omega)+pi/4+psi+2delta_1.",
        "endpoint_identity": "At the ideal cutoff endpoints Omega=2pi N^2 and 2pi(N+1)^2, the increment of 2chi_N is 2pi G_N, where G_N=(N+1)^2[1-2log(1+1/N)]-N^2.",
        "ideal_decrement": "Put u=1/N. The alternating lower bound log(1+u)>=u-u^2/2+u^3/3-u^4/4 yields G_N<=-2-(7/16)u. Thus the ideal increment obeys Delta(2chi_N)<=-4pi-7pi/(8N).",
        "endpoint_scale": "At an endpoint with a=m>=2, x=4pi[m^2-1/(32L^2)]>9m^2, using pi>3. Also N>=ceil(a(50))>2^25 because e>2.",
        "physical_error": "Replacing the ideal endpoint frequencies by Omega=2pi a^2-epsilon contributes less than 1/N^2 to Delta F_N. The lock phases contribute |Delta psi|<2/(9N^2). Since |arg(1+d_1)|<2|d_1|<4378/x, the term 2Delta delta_1 contributes less than 17512/(9N^2). The total perturbation of Delta(2chi_N) is therefore less than 1947/N^2.",
        "reverse_turn": "Since 1947/N^2<7pi/(8N) for N>2^25, every complete physical q=1 cell satisfies Delta chi_N<-2pi. Moreover Delta chi_N tends to -2pi as N tends to infinity.",
        "absolute_circle": "The continuous absolute terminal carrier eta exp(i phi_N)=exp(i chi_N) covers the full unit circle on every complete cell, now in the reverse direction. This is source-inclusive and is distinct from the much faster positive winding of exp(i phi_N).",
        "phase_inventory": "The previous Delta phi_N>2pi and the new Delta chi_N<-2pi imply Delta theta_eta=Delta chi_N-Delta phi_N<-4pi. The transpose phase phi_N+2theta_eta=2chi_N-phi_N has increment below -6pi. Thus the terminal, source, absolute-terminal, and transpose phases all cover the circle; their correlations, not a fixed sector, are the remaining content.",
        "reverse_full_turn_theorems": 1,
        "source_inclusive_circle_theorems": 1,
        "phase_inventory_corollaries": 2,
        "physical_error_coefficient": 1947,
    }


def centre_symmetry_certificate() -> dict[str, str | int]:
    a_t, b, u_n, phi, theta = sp.symbols(
        "A_T b u_N phi theta", positive=True, real=True
    )
    a_h = a_t * (1 + 2 * u_n / b)
    original = a_h * sp.sin(phi) + a_t * sp.sin(phi + 2 * theta)
    factored = 2 * a_t * (
        sp.cos(theta) * sp.sin(phi + theta) + (u_n / b) * sp.sin(phi)
    )
    require_zero(original - factored, "centre symmetry")

    return {
        "anchor_ratio": "For b=log N, d=log(a^2/N)=b+2u_N, the raw anchors satisfy A_H/A_T=d/b=1+2u_N/b and A_H-A_T=D_N b d u_N/(2pi)>=0.",
        "actual_current": "At xi=Omega, write tau_0=rho exp(i[phi_N-pi/2]), rho>0. The actual-source raw-centre current is R_center=-(rho|nu|^2/2)[A_H sin(phi_N)+A_T sin(phi_N+2theta_eta)].",
        "symmetry_factorization": "The two-channel centre has the exact saddle symmetry R_center=-rho|nu|^2 A_T[cos(theta_eta)sin(phi_N+theta_eta)+(u_N/b)sin(phi_N)]. The first term couples the source phase to the source-inclusive absolute terminal phase chi_N; the second is the exact Hermitian/transpose anchor mismatch.",
        "branch_free": "Equivalently, R_center=-rho|nu|^2 A_T[Re(eta)Im(eta exp(i phi_N))+(u_N/b)Im(exp(i phi_N))]. Choosing omega_c=(-1)^Neta instead gives the same product, so this form is parity invariant and branch free.",
        "mismatch_bound": "Across a complete cell, 0<=u_N<=log(1+1/N)<1/N and b=log N>1/2, hence 0<=u_N/b<2/N. The raw centres are therefore an extremely close equal-anchor pair, but this small coefficient does not control the current near zeros of Re(eta) or Im(eta exp(i phi_N)).",
        "tie_guard": "D_N changes at internal odd-roster ties. The ratio and factorization hold on both adjacent raw-centre charts, but the centre is not itself a tie-invariant physical observable. Any sign theorem must restore the reciprocal residual, terminal survivor, correction rows, affine row, and moved-tail defect before crossing a tie.",
        "open_correlation": "The next falsifiable centre-level question is whether the Xi-specific joint path (Re eta, Im[eta exp(i phi_N)]) has a signed product budget strong enough to dominate the O(u_N/log N) mismatch on roster interiors. Circle coverage of either factor separately cannot answer it.",
        "exact_center_factorizations": 2,
        "anchor_mismatch_bounds": 1,
        "live_joint_phase_targets": 1,
    }


def build_rows(source: dict, turn: dict, centre: dict) -> list[GateRow]:
    return [
        GateRow("asl_01_normalizer", "source normalizer", "proved", "The physical source phase and modulus are fixed exactly.", source["normalizer"], "This is the q=1 source normalization."),
        GateRow("asl_02_lift", "phase lift", "proved", "The source and terminal phases admit continuous unwrapped lifts.", source["continuous_phase"], "Only complete fixed-N cells are used."),
        GateRow("asl_03_lock", "stationary lock", "proved", "The reciprocal stationary phase is conjugate-locked to the normalizer.", source["stationary_lock"], "The lock is imported with its O(1/x) error."),
        GateRow("asl_04_identity", "absolute phase identity", "proved", "The source-inclusive terminal phase has an exact scalar normal form.", turn["lock_identity"], "No phase sampling is used."),
        GateRow("asl_05_endpoint", "ideal endpoint", "proved", "The ideal complete-cell increment is one elementary logarithmic expression.", turn["endpoint_identity"], "Physical epsilon corrections are added later."),
        GateRow("asl_06_decrement", "ideal decrement", "proved", "The ideal absolute phase decreases by more than one full turn.", turn["ideal_decrement"], "The alternating logarithm bound is uniform for N>=2."),
        GateRow("asl_07_scale", "physical scale", "proved", "The actual q=1 cells make every endpoint correction tiny relative to the turn margin.", turn["endpoint_scale"], "The coarse N>2^25 floor is sufficient."),
        GateRow("asl_08_error", "physical perturbation", "proved", "All physical phase-lock and first-anchor corrections fit one O(N^-2) budget.", turn["physical_error"], "The budget is intentionally coarse."),
        GateRow("asl_09_reverse", "reverse full turn", "proved", "The actual source-inclusive terminal phase loses more than 2pi per complete cell.", turn["reverse_turn"], "This is a phase theorem, not a current sign."),
        GateRow("asl_10_circle", "absolute circle", "proved", "The actual absolute terminal carrier covers every unit direction.", turn["absolute_circle"], "Coverage follows from continuity and endpoint decrement."),
        GateRow("asl_11_inventory", "phase inventory", "proved", "The source and transpose phases also make complete turns.", turn["phase_inventory"], "Separate winding does not determine their joint product."),
        GateRow("asl_12_ratio", "anchor ratio", "proved", "The Hermitian and transpose raw anchors differ by exactly 2u_N/log N.", centre["anchor_ratio"], "D_N remains roster dependent."),
        GateRow("asl_13_current", "actual centre", "proved", "The physical source can be inserted into the bare two-channel centre exactly.", centre["actual_current"], "Residual and correction rows are excluded here."),
        GateRow("asl_14_symmetry", "centre symmetry", "proved", "The bare centre factors into the source and absolute-terminal phases plus one mismatch.", centre["symmetry_factorization"], "No sign of the phase product is asserted."),
        GateRow("asl_15_branch", "branch-free form", "proved", "The symmetry has a parity-invariant branch-free expression.", centre["branch_free"], "A continuous lift is used only for winding counts."),
        GateRow("asl_16_mismatch", "mismatch bound", "proved", "The unequal-anchor term is below 2/N relative to A_T.", centre["mismatch_bound"], "Smallness alone fails at phase-product zeros."),
        GateRow("asl_17_ties", "roster guard", "proved", "The raw-centre factorization cannot be promoted through ties by itself.", centre["tie_guard"], "The complete joined current is the invariant object."),
        GateRow("asl_18_joint", "joint phase target", "open", "A signed Xi-specific joint-phase product estimate remains open.", centre["open_correlation"], "Marginal circle coverage is insufficient."),
        GateRow("asl_19_complete", "complete current", "open", "The joint phase theorem must be recomposed with every retained residual and correction row.", "Prove the actual-phase inequality for the complete tie-invariant pair (h,t), using the new locked phase coordinate before any absolute values.", "No contact-box, all-q, Lambda<=0, or RH theorem is claimed."),
    ]


def render_note(payload: dict) -> str:
    source = payload["exact"]["source_phase"]
    turn = payload["exact"]["reverse_turn"]
    centre = payload["exact"]["centre_symmetry"]
    lines = [
        "# Newman C1 Actual-Source Saddle-Phase Lock and Centre Symmetry Gate",
        "",
        "Date: 2026-08-04",
        "",
        "Status: exact actual-source normalization, reverse full-turn theorem, and bare-centre symmetry proved; joined current sign open; not a proof of RH.",
        "",
        "## Source phase",
        "",
        source["normalizer"],
        "",
        source["stationary_lock"],
        "",
        "## Reverse turn",
        "",
        "```text",
        "2chi_N=F_N(Omega)+pi/4+psi+2delta_1,",
        "F_N(v)=v[1+log(2pi N^2/v)],",
        "Delta chi_N<-2pi.",
        "```",
        "",
        turn["ideal_decrement"],
        "",
        turn["physical_error"],
        "",
        turn["absolute_circle"],
        "",
        turn["phase_inventory"],
        "",
        "## Centre symmetry",
        "",
        centre["actual_current"],
        "",
        "```text",
        "R_center=-rho|nu|^2 A_T[",
        "  cos(theta_eta)sin(phi_N+theta_eta)",
        "  +(u_N/log N)sin(phi_N)].",
        "```",
        "",
        centre["branch_free"],
        "",
        centre["mismatch_bound"],
        "",
        centre["tie_guard"],
        "",
        "## Gate rows",
        "",
        "| id | role | state | claim |",
        "|---|---|---|---|",
    ]
    for row in payload["rows"]:
        lines.append(f"| {row['id']} | {row['role']} | {row['readiness']} | {row['claim']} |")
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
    _, source_audit = load_and_audit_sources()
    source = source_phase_certificate()
    turn = reverse_turn_certificate()
    centre = centre_symmetry_certificate()
    rows = build_rows(source, turn, centre)
    require(len(rows) == 19, "row count")
    require(sum(row.readiness == "open" for row in rows) == 2, "open row count")

    success = (
        "built Newman C1 actual-source saddle-phase lock/centre-symmetry gate: "
        "19 rows, 0 issues, 1 reverse full-turn theorem, 1 source-inclusive "
        "circle theorem, 2 exact centre factorizations, 1 live joint-phase target"
    )
    payload = {
        "kind": KIND,
        "schema_version": 1,
        "date": "2026-08-04",
        "status": "actual-source normalization, reverse complete-cell turn, and exact bare-centre symmetry proved; joined current sign open",
        "source_sha256": {key: item["sha256"] for key, item in source_audit.items()},
        "source_audit": source_audit,
        "exact": {
            "source_phase": source,
            "reverse_turn": turn,
            "centre_symmetry": centre,
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "source_normalizer_identifications": source["source_normalizer_identifications"],
            "continuous_phase_lifts": source["continuous_phase_lifts"],
            "reverse_full_turn_theorems": turn["reverse_full_turn_theorems"],
            "source_inclusive_circle_theorems": turn["source_inclusive_circle_theorems"],
            "phase_inventory_corollaries": turn["phase_inventory_corollaries"],
            "physical_error_coefficient": turn["physical_error_coefficient"],
            "exact_center_factorizations": centre["exact_center_factorizations"],
            "anchor_mismatch_bounds": centre["anchor_mismatch_bounds"],
            "live_joint_phase_targets": centre["live_joint_phase_targets"],
            "open_rows": 2,
        },
        "next_action": "Derive a local scale-separation theorem for the joint path Re(eta)*Im(eta exp(i phi_N)) on roster interiors. Use the exact reverse-turn coordinate chi_N and the much faster terminal phase phi_N to determine whether the dominant product must take both signs, while retaining the O(u_N/log N) mismatch. If that centre-level sign oscillates, promote it as a route obstruction and move directly to the complete reciprocal residual; if a one-sided Xi correlation survives, insert it into the tie-invariant complete (h,t) current before estimating corrections.",
        "pi_provenance": "The pi in Q(v)=exp(i[v log(v/(2pi))-v-pi/4]) is the completed-zeta reciprocal stationary phase, with pi/4 its negative-curvature signature. The factor 2pi in K_N=2pi N^2 and T_0=2pi a^2 is the same completed-zeta/Fourier normalization. The threshold 2pi is one period of exp(i chi_N). Only elementary pi>3 is used in the coarse endpoint-error audit.",
        "proof_boundary": "This gate proves the exact physical source normalizer, a source-inclusive reverse full turn on every complete q=1 cutoff cell, the induced source and transpose winding inventory, and the exact parity-invariant bare-centre symmetry with its O(u_N/log N) anchor mismatch. It proves no sign for the joint phase product, raw-centre sign theorem, signed reciprocal residual estimate, complete actual-source current inequality, correction-row bound, tie-spliced current theorem, all-q transport, contact-box exclusion, Q209 shell, cofinal successor, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
        "success": success,
    }

    atomic_write(RESULT_PATH, json.dumps(payload, indent=2) + "\n")
    atomic_write(NOTE_PATH, render_note(payload))
    print(success)


if __name__ == "__main__":
    main()
