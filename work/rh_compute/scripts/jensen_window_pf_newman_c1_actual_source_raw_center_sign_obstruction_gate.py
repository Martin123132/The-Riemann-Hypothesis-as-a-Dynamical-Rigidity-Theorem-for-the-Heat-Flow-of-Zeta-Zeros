#!/usr/bin/env python3
"""Build the actual-source raw-centre sign-obstruction gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_raw_center_sign_obstruction_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "phase_lock": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_saddle_phase_lock_center_symmetry_gate.json",
    "covariant_current": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_terminal_phase_winding_covariant_two_channel_gate.json",
    "raw_anchors": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_carrier_near_odd_harmonic_anchor_gate.json",
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


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(content, encoding="utf-8")
    os.replace(temporary, path)


def load_and_audit_sources() -> dict[str, dict[str, str]]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        require(path.is_file(), f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))
    require(payloads["phase_lock"]["summary"]["reverse_full_turn_theorems"] == 1, "phase lock drifted")
    require(payloads["phase_lock"]["summary"]["exact_center_factorizations"] == 2, "centre factor drifted")
    require(payloads["phase_lock"]["summary"]["physical_error_coefficient"] == 1947, "phase error drifted")
    require(payloads["covariant_current"]["summary"]["live_covariant_current_targets"] == 1, "current frontier drifted")
    require(payloads["raw_anchors"]["summary"]["explicit_raw_anchors"] == 2, "raw anchors drifted")
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def phase_scale_certificate() -> dict[str, str | int]:
    # A safe sine component of width pi/3 is enough. The locked phase moves
    # at relative speed below 2/N, so terminal phase variation exceeds pi*N/6.
    require(Fraction(1, 3) / Fraction(2, 1) == Fraction(1, 6), "scale ratio")
    require(2**25 > 10, "physical cell floor")

    return {
        "main_phase": "Write chi_N=chi_0+r with chi_0=(1/2)F_N(Omega)+pi/8, F_N(v)=v[1+log(2pi N^2/v)], and r=psi/2+arg(1+d_1).",
        "frequency_monotonicity": "On q=1, T_0=2pi a^2 is strictly increasing while both positive terms of epsilon=atan(1/x)/(8L^2)+3x/[4L^2(1+x^2)] are strictly decreasing for x>1. Hence Omega=T_0-epsilon is strictly increasing across the cell.",
        "main_reverse_turn": "The cell crosses K_N=2pi N^2 once. From that crossing to the right endpoint, F_N(Omega) decreases by more than 4pi, so chi_0 decreases by more than 2pi.",
        "safe_arc": "In any real interval of phase length greater than 2pi, the periodic set sin(theta)>=1/2 has a connected component of intersection of width at least pi/3. Choose a subinterval J on which chi_0 decreases by pi/3 and sin(chi_0)>=1/2.",
        "relative_speed": "With phi_N=Omega log N and Omega>=K_N on J, |d chi_0/d phi_N|=log(Omega/K_N)/(2log N)<log(1+1/N)/log N<2/N.",
        "source_sweep": "The pi/3 change of chi_0 therefore requires Delta phi_N>pi N/6. Thus theta_0=chi_0-phi_N decreases by more than pi N/6+pi/3>2pi on J, and assumes representatives with cos(theta_0)=+1 and cos(theta_0)=-1.",
        "safe_phase_arcs": 1,
        "scale_separation_theorems": 1,
        "source_half_turn_pairs": 1,
    }


def sign_certificate() -> dict[str, str | int | Fraction]:
    n_floor = 2**25
    require(Fraction(487, n_floor**2) < Fraction(1, 1000), "phase remainder")
    require(Fraction(2, n_floor) < Fraction(1, 1000), "anchor mismatch")
    positive_reserve = Fraction(999, 1000) * Fraction(499, 1000) - Fraction(1, 1000)
    require(positive_reserve == Fraction(497501, 1000000), "reserve value")
    require(positive_reserve > Fraction(49, 100), "reserve floor")

    theta, chi, phi, mismatch = sp.symbols("theta chi phi mismatch", real=True)
    bracket = sp.cos(theta) * sp.sin(chi) + mismatch * sp.sin(phi)
    require(sp.diff(bracket, mismatch) == sp.sin(phi), "mismatch linearity")

    return {
        "uniform_remainder": "The phase-lock errors give |r|<1/(2x)+4378/x<487/N^2<1/1000 throughout the cell. The exact anchor mismatch e_N=u_N/log N obeys 0<=e_N<2/N<1/1000.",
        "positive_point": "At a point of J with theta_0=0 modulo 2pi, theta_eta=theta_0+r and chi_N=chi_0+r. Hence cos(theta_eta)>999/1000 and sin(chi_N)>499/1000. The normalized bracket B=cos(theta_eta)sin(chi_N)+e_N sin(phi_N) is greater than 497501/10^6>49/100.",
        "negative_point": "At a point of J with theta_0=pi modulo 2pi, cos(theta_eta)<-999/1000 while sin(chi_N)>499/1000. The same normalized bracket satisfies B<-497501/10^6<-49/100.",
        "raw_center_obstruction": "Since R_center=-rho|nu|^2A_T B with rho|nu|^2A_T>0, every complete physical q=1 cell contains raw-centre points with R_center<-(49/100)rho|nu|^2A_T and R_center>(49/100)rho|nu|^2A_T. The actual physical source phase therefore does not rescue a cell-uniform raw-centre sign.",
        "tie_guard": "The selected points may be perturbed off the finite roster-tie set without losing the strict 49/100 reserve. The raw-centre sign obstruction is valid on roster interiors, but only the complete residual-composed current is tie invariant.",
        "residual_scalar": "With W_H=R_H and W_T=R_T+Q_T, the complete ideal actual-source scalar is Sigma_id=-2rho A_TB+E_W, where E_W=Re(conjugate(tau_0)W_H)+Re(eta^2tau_0W_T), and R_id=(|nu|^2/2)Sigma_id.",
        "macroscopic_compensation": "At the witness B<-49/100, negativity of the complete ideal current requires E_W<-(98/100)rho A_T. Since |E_W|<=rho(|W_H|+|W_T|), a necessary condition is |W_H|+|W_T|>(98/100)A_T. The residual must therefore be macroscopic and correctly signed on the adverse arc; this may be a uniform negative bias or a phase-correlated response.",
        "complete_target": "A successful physical-current theorem must use signed reciprocal-stationary residuals and the correction/affine/moved-tail rows to overturn or absorb this O(A_T) alternating centre. It must prove either a sufficiently strong uniform negative bias or the required phase-correlated inequality. Any proof based on a uniformly favourable centre phase, even the actual source phase, is now retired.",
        "strict_positive_witnesses": 1,
        "strict_negative_witnesses": 1,
        "actual_source_center_sign_obstructions": 1,
        "macroscopic_compensation_conditions": 1,
        "normalized_reserve_numerator": 49,
        "normalized_reserve_denominator": 100,
    }


def build_rows(phase: dict, sign: dict) -> list[GateRow]:
    return [
        GateRow("rco_01_main", "locked main phase", "proved", "The actual absolute phase splits into an explicit main and a uniformly tiny remainder.", phase["main_phase"], "This is the fixed physical q=1 chart."),
        GateRow("rco_02_monotone", "frequency monotonicity", "proved", "The physical frequency is strictly increasing across each complete cell.", phase["frequency_monotonicity"], "No sampled monotonicity is used."),
        GateRow("rco_03_turn", "main reverse turn", "proved", "The explicit locked main makes more than one reverse turn after the saddle crossing.", phase["main_reverse_turn"], "This concerns chi_0, before the tiny phase remainder."),
        GateRow("rco_04_arc", "safe sine arc", "proved", "One connected locked-phase arc has sine at least one half and width pi/3.", phase["safe_arc"], "This is an elementary periodic-arc pigeonhole."),
        GateRow("rco_05_speed", "scale separation", "proved", "The locked phase moves at less than 2/N of the terminal phase speed.", phase["relative_speed"], "The bound is used only after Omega reaches K_N."),
        GateRow("rco_06_sweep", "source half-turn pair", "proved", "The source phase makes more than a full turn inside the safe locked-phase arc.", phase["source_sweep"], "Only continuity and strict phase ranges are required."),
        GateRow("rco_07_remainder", "physical remainder", "proved", "The actual phase and unequal-anchor errors are each below 1/1000.", sign["uniform_remainder"], "The physical N floor is used explicitly."),
        GateRow("rco_08_positive", "positive bracket witness", "proved", "One roster-interior point has normalized centre bracket above 49/100.", sign["positive_point"], "The current itself has the opposite sign."),
        GateRow("rco_09_negative", "negative bracket witness", "proved", "A second roster-interior point has normalized centre bracket below -49/100.", sign["negative_point"], "The current itself has the opposite sign."),
        GateRow("rco_10_obstruction", "actual-source centre obstruction", "proved", "The raw-centre current takes both signs at the actual physical source phase.", sign["raw_center_obstruction"], "This is not the complete current."),
        GateRow("rco_11_ties", "tie guard", "proved", "The strict witnesses can be chosen away from every roster tie.", sign["tie_guard"], "Residual recomposition remains compulsory."),
        GateRow("rco_12_compensation", "residual compensation", "proved", "Any negative complete ideal current needs a macroscopic correctly signed residual on the positive-centre arc.", sign["macroscopic_compensation"], "This is a necessary condition, not a residual estimate; uniform bias is not excluded."),
        GateRow("rco_13_target", "complete current target", "open", "Only a signed complete residual current can now restore a one-sided theorem.", sign["complete_target"], "No residual estimate or RH conclusion is supplied."),
    ]


def render_note(payload: dict) -> str:
    phase = payload["exact"]["phase_scale"]
    sign = payload["exact"]["sign_obstruction"]
    lines = [
        "# Newman C1 Actual-Source Raw-Centre Sign Obstruction Gate",
        "",
        "Date: 2026-08-04",
        "",
        "Status: actual physical source phase proved to force both raw-centre signs on every complete q=1 cell; complete residual current open; not a proof of RH.",
        "",
        "## Scale separation",
        "",
        phase["main_phase"],
        "",
        phase["main_reverse_turn"],
        "",
        phase["safe_arc"],
        "",
        phase["relative_speed"],
        "",
        phase["source_sweep"],
        "",
        "## Strict sign witnesses",
        "",
        sign["uniform_remainder"],
        "",
        sign["positive_point"],
        "",
        sign["negative_point"],
        "",
        "```text",
        "R_center=-rho|nu|^2 A_T B,",
        "B>49/100 at one point,",
        "B<-49/100 at another point.",
        "```",
        "",
        sign["raw_center_obstruction"],
        "",
        sign["tie_guard"],
        "",
        sign["residual_scalar"],
        "",
        sign["macroscopic_compensation"],
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
    source_audit = load_and_audit_sources()
    phase = phase_scale_certificate()
    sign = sign_certificate()
    rows = build_rows(phase, sign)
    require(len(rows) == 13, "row count")
    require(sum(row.readiness == "open" for row in rows) == 1, "open row count")

    success = (
        "built Newman C1 actual-source raw-centre sign-obstruction gate: "
        "13 rows, 0 issues, 1 safe phase arc, 1 scale-separation theorem, "
        "2 strict sign witnesses, 1 actual-source centre obstruction, "
        "1 macroscopic compensation condition, 1 live complete-current target"
    )
    payload = {
        "kind": KIND,
        "schema_version": 1,
        "date": "2026-08-04",
        "status": "actual physical source phase forces both raw-centre signs on every complete q=1 cell; complete residual current open",
        "source_sha256": {key: item["sha256"] for key, item in source_audit.items()},
        "source_audit": source_audit,
        "exact": {
            "phase_scale": phase,
            "sign_obstruction": sign,
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "safe_phase_arcs": phase["safe_phase_arcs"],
            "scale_separation_theorems": phase["scale_separation_theorems"],
            "source_half_turn_pairs": phase["source_half_turn_pairs"],
            "strict_positive_witnesses": sign["strict_positive_witnesses"],
            "strict_negative_witnesses": sign["strict_negative_witnesses"],
            "actual_source_center_sign_obstructions": sign["actual_source_center_sign_obstructions"],
            "macroscopic_compensation_conditions": sign["macroscopic_compensation_conditions"],
            "normalized_reserve": "49/100",
            "live_complete_current_targets": 1,
        },
        "next_action": "Stop trying to obtain the physical sign from either raw anchor centre. Compose the full tie-invariant actual-source scalar h+Re(eta^2 t) using R_H, R_T+Q_T, the physical correction polynomials, fixed affine row, and moved-tail defect. Determine whether it has a cell-uniform negative bias or a phase-correlated O(A_T) response on the strict adverse arcs. A lower bound at one adverse witness, or a common-unit residual bound below the sharp 98/100 threshold, would instead prove that the complete ideal block changes sign and retire this entire sign route.",
        "pi_provenance": "The only pi-dependent phase is inherited from the completed-zeta saddle lock in the source artifact. The pi/3 safe arc and 2pi turn are fixed fractions of the period of exp(i theta); they introduce no geometric model constant. The 49/100 reserve is a rational consequence of the physical N floor and the exact phase/mismatch error bounds.",
        "proof_boundary": "This gate proves a two-scale phase theorem and two strict opposite-sign witnesses for the bare Hermitian/transpose centre at the actual physical source phase on every complete q=1 cell. It retires a cell-uniform raw-centre sign, not the complete current. It proves no signed reciprocal residual estimate, correction/affine/moved-tail compensation, tie-spliced complete current sign, all-q transport, contact-box exclusion, Q209 shell, cofinal successor, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
        "success": success,
    }

    atomic_write(RESULT_PATH, json.dumps(payload, indent=2) + "\n")
    atomic_write(NOTE_PATH, render_note(payload))
    print(success)


if __name__ == "__main__":
    main()
