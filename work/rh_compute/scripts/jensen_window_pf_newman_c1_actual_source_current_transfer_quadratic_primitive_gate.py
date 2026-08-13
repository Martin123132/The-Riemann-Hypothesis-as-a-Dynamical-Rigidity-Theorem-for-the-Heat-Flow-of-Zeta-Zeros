#!/usr/bin/env python3
"""Build the current/transfer ownership and quadratic-primitive gate."""

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
KIND = "jensen_window_pf_newman_c1_actual_source_current_transfer_quadratic_primitive_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "homotopy": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_ideal_cubic_reciprocal_terminal_full_support_homotopy_gate.json",
    "observation_image": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_observation_image_compression_gate.json",
    "terminal_split": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_terminal_bulk_phase_completion_gate.json",
    "ledger": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_complete_compensation_ledger_gate.json",
    "edge_remote": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_edge_affine_remote_budget_gate.json",
    "moving_transfer": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_moving_tail_budget_gate.json",
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


def rows_by_id(payload: dict) -> dict[str, dict]:
    return {row["id"]: row for row in payload["rows"]}


def load_and_audit_sources() -> dict[str, dict[str, str]]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        require(path.is_file(), f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))

    homotopy_rows = rows_by_id(payloads["homotopy"])
    require(
        "tilde(Psi)'" in homotopy_rows["rth_07_move_defect"]["certificate"],
        "moving-current identity drifted",
    )
    require(
        "omega_C(xi)-omega_C(Omega)" in homotopy_rows["rth_07_move_defect"]["certificate"],
        "moving displacement drifted",
    )

    observation_rows = rows_by_id(payloads["observation_image"])
    require(
        "x_Vz_N+x_Nz_V-x_Az_Q-x_Qz_A" in observation_rows["mfic_14_quadratic"]["certificate"],
        "quadratic observation identity drifted",
    )

    split = payloads["terminal_split"]
    split_text = json.dumps(split, sort_keys=True)
    require("epsilon=a+rho" in split_text, "terminal/nonterminal split drifted")
    require("Re(H_aa)=0" in split_text, "pure-terminal ownership drifted")

    ledger = payloads["ledger"]
    require(
        "E_X=2R_X/|nu|^2" in ledger["exact"]["physical_partition"]["remaining"],
        "scalar normalization drifted",
    )
    require(
        "Sigma_full" in ledger["exact"]["normalization"]["scalar"],
        "current ledger drifted",
    )

    edge = payloads["edge_remote"]
    require(edge["summary"]["closed_secondary_fraction"] == "280759/858000", "five-package budget drifted")
    require(edge["exact"]["source_ratio"] == "B_0T_0/rho<198/329<61/100.", "source ratio drifted")

    moving = payloads["moving_transfer"]
    require(
        "-varepsilon integral_0^1D_move" in moving["exact"]["definition"],
        "integrated moving package drifted",
    )
    require(moving["summary"]["moving_tail_budget_denominator"] == 1_000_000, "transfer budget drifted")

    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def symbolic_certificate() -> dict[str, str | int]:
    x = sp.symbols("x0:4", real=True)
    y = sp.symbols("y0:4", real=True)
    xv = sp.Matrix(x)
    yv = sp.Matrix(y)
    j = sp.Rational(1, 2) * sp.Matrix(
        [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]]
    )
    determinant = x[0] * x[1] - x[2] * x[3]
    require(sp.expand((xv.T * j * xv)[0]) == determinant, "J determinant")
    derivative = sum(sp.diff(determinant, x[index]) * y[index] for index in range(4))
    expected = x[0] * y[1] + x[1] * y[0] - x[2] * y[3] - x[3] * y[2]
    require(sp.expand(derivative) == expected, "determinant derivative")
    require(sp.expand(2 * (xv.T * j * yv)[0]) == expected, "J derivative")

    h, r, k = sp.symbols("h r k", positive=True, real=True)
    p_star = r**2 / 2
    s0 = sp.Rational(16, 7) * h ** sp.Rational(-1, 2)
    t1 = 4 * h ** sp.Rational(3, 2) * k**2
    carrier = sp.simplify(8 * p_star**2 * s0 * t1)
    edge = sp.simplify(sp.Rational(24, 5) * p_star * t1)
    require(carrier == sp.Rational(128, 7) * h * r**4 * k**2, "pointwise carrier reduction")
    require(edge == sp.Rational(48, 5) * h ** sp.Rational(3, 2) * r**2 * k**2, "pointwise edge reduction")

    vr, vi, nr, ni, ar, ai, qr, qi = sp.symbols("vr vi nr ni ar ai qr qi", real=True)
    zv = vr + sp.I * vi
    zn = nr + sp.I * ni
    za = ar + sp.I * ai
    zq = qr + sp.I * qi
    polarized = sp.re(zv * sp.conjugate(zn) - za * sp.conjugate(zq) + zv * zn - za * zq).expand(complex=True)
    require(sp.simplify(polarized - 2 * (vr * nr - ar * qr)) == 0, "complex determinant polarization")

    return {
        "observation_order": "e=(e_V,e_mathcalN,e_A,e_Q)^T and J=(1/2)[[0,1,0,0],[1,0,0,0],[0,0,0,-1],[0,0,-1,0]].",
        "determinant": "D_e=e^TJe=e_Ve_mathcalN-e_Ae_Q.",
        "primitive": "D_e'=e_Vdot(e_mathcalN)+e_mathcalNdot(e_V)-e_Adot(e_Q)-e_Qdot(e_A)=R_quad^cur.",
        "pointwise_move": "At xi=Omega, Delta omega_C=0 and D_move^cur=2(omega_E+omega_B+omega_C)^TJdot(omega_C).",
        "pointwise_reductions": "8P_*^2S_0T_1<(128/7)hR_*^4K^2 and (24/5)(B_0T_0)P_*T_1<(48/5)(B_0T_0)h^(3/2)R_*^2K^2.",
        "complex_primitive": "For U_X=(nu c_xi/|nu|)Z_X, Q_quad=Re{U_V conjugate(U_mathcalN)-U_A conjugate(U_Q)+U_VU_mathcalN-U_AU_Q}.",
        "symbolic_audits": 6,
    }


def exact_certificate() -> dict[str, str | int]:
    h0 = Fraction(1, 72_000_000_000)
    l0 = Fraction(50)
    r0 = Fraction(53, 2)

    exp25_lower = sum(Fraction(25) ** index / math.factorial(index) for index in range(61))
    require(exp25_lower > 72_000_000_000, "exp(25) floor")
    require(72_000_000_000**7 > 13_000**18, "fractional-power floor")

    current_main = Fraction(422_400, 7 * 329 * 14_300) * r0**4 / (l0 * 13_000)
    current_edge = Fraction(2_928, 5 * 100 * 14_300) * h0 * r0**2
    current_total = current_main + current_edge
    require(current_total < Fraction(1, 100), "pointwise moving budget")

    closed = Fraction(280_759, 858_000) + Fraction(1, 100)
    reserve = Fraction(98, 100) - closed
    sufficient = Fraction(98, 100) + closed
    require(closed == Fraction(289_339, 858_000), "corrected closed fraction")
    require(reserve == Fraction(551_501, 858_000), "corrected necessary reserve")
    require(sufficient == Fraction(1_130_179, 858_000), "corrected sufficient threshold")

    return {
        "ownership_correction": "D_move(xi)=tilde(Psi)'(xi)-Psi_B'(xi) is a pointwise current defect. R_move^tr=-varepsilon integral_0^1D_move(T_0-theta varepsilon)dtheta is its integrated transfer. The two quantities must have distinct names and budgets.",
        "prior_gate_boundary": "The 1/1000000 theorem remains valid for E_move^tr=2R_move^tr/|nu|^2. It does not bound the pointwise E_move^cur used in Sigma_full.",
        "physical_current": "At the physical point xi=Omega, Delta omega_C(Omega)=0, so R_move^cur=D_move(Omega)=2[omega_E+omega_B(Omega)+omega_C(Omega)]^TJdot(omega_C)(Omega).",
        "pointwise_bound": "|E_move^cur|<8P_*^2S_0T_1+(24/5)(B_0T_0)P_*T_1<(128/7)hR_*^4K^2+(48/5)(B_0T_0)h^(3/2)R_*^2K^2.",
        "normalized_bound": "After K^2<=h^(-1/9), rho>(329/3300)Lexp(-L/4), exp(L/4)<h^(-1/2), B_0T_0/rho<61/100, and A_T>14300, the two ratios are below [422400/(7*329*14300)]h^(7/18)R_*^4/L and [2928/(5*100*14300)]h^(25/18)R_*^2.",
        "monotonicity": "After h<exp(-L/2), the logarithmic derivatives are -7/36+4/(L+3)-1/L and -25/36+2/(L+3), both negative for L>=50.",
        "endpoint_check": "At L=50, h<1/72000000000, (72000000000)^7>13000^18 gives h^(7/18)<1/13000, and h^(25/18)<h. The resulting exact rational sum is below 1/100.",
        "current_budget": "|E_move^cur|<rho A_T/100.",
        "corrected_closed_budget": "|E_R+E_term+E_aa+E_Delta+E_aff^C+E_move^cur|<(289339/858000)rho A_T.",
        "quadratic_constituents": "E_quad^cur=(2/|nu|^2){e_Ve_(mathcalN,xi)+e_mathcalNe_(V,xi)-e_Ae_(Q,xi)-e_Qe_(A,xi)}, with e containing only the nonterminal remainder after epsilon=a+rho_obs.",
        "quadratic_primitive": "Q_quad=(2/|nu|^2)(e_Ve_mathcalN-e_Ae_Q), and E_quad^cur=dQ_quad/dxi exactly on every fixed chart.",
        "quadratic_polarization": "For normalized complex nonterminal amplitudes U_X, Q_quad=Re{U_V conjugate(U_mathcalN)-U_A conjugate(U_Q)+U_VU_mathcalN-U_AU_Q}; differentiation gives the Hermitian phase-difference and transpose phase-sum currents without discarding either phase.",
        "quadratic_transfer": "For xi_theta=T_0-theta varepsilon and Omega=T_0-varepsilon, -varepsilon integral_0^1E_quad^cur(xi_theta)dtheta=Q_quad(Omega)-Q_quad(T_0).",
        "ownership": "E_quad contains no terminal a factor. E_aa is the a-a block; E_term and the finite/remote/correction channels contain terminal-nonterminal terms; E_aff is linear; E_move^cur is the terminal-motion current defect.",
        "signed_reclassification": "Define E_main^cur=E_FAN+E_quad^cur. With the corrected absolute budget, adverse-witness negativity necessarily requires E_main^cur<-(551501/858000)rho A_T. The stronger E_main^cur<-(1130179/858000)rho A_T is sufficient independently of the absolute package signs.",
        "live_obligation": "Choose one coherent route: prove the signed pointwise inequality for E_FAN+E_quad^cur, or use the transfer primitive and bound the endpoint determinant difference jointly with tilde(Psi)(T_0). No separate |E_quad| budget is assumed.",
        "current_main_numerator": current_main.numerator,
        "current_main_denominator": current_main.denominator,
        "current_edge_numerator": current_edge.numerator,
        "current_edge_denominator": current_edge.denominator,
        "closed_numerator": closed.numerator,
        "closed_denominator": closed.denominator,
        "reserve_numerator": reserve.numerator,
        "reserve_denominator": reserve.denominator,
        "sufficient_numerator": sufficient.numerator,
        "sufficient_denominator": sufficient.denominator,
    }


def build_rows(exact: dict[str, str | int]) -> list[GateRow]:
    return [
        GateRow("ctq_01_domain", "domain", "proved", "Pointwise current and integrated transfer are distinct objects.", str(exact["ownership_correction"]), "No numerical estimate is used in this distinction."),
        GateRow("ctq_02_prior", "corrigendum", "guard_validated", "The prior millionth-scale theorem is retained only in its valid transfer role.", str(exact["prior_gate_boundary"]), "Its insertion into the pointwise current ledger is superseded."),
        GateRow("ctq_03_current", "current ownership", "proved", "The physical moving current simplifies at xi=Omega.", str(exact["physical_current"]), "The vanished displacement term is not discarded away from Omega."),
        GateRow("ctq_04_bound", "pointwise envelope", "proved", "The physical moving current has two source-normalized terms.", str(exact["pointwise_bound"]), "No transfer factor varepsilon is inserted."),
        GateRow("ctq_05_ratio", "anchor ratio", "proved", "Both pointwise terms admit explicit anchor ratios.", str(exact["normalized_bound"]), "The growing-prefix and source floors are retained."),
        GateRow("ctq_06_monotone", "monotonicity", "proved", "The two normalized majorants decrease from L=50.", str(exact["monotonicity"]), "This is a uniform fixed-q1 statement."),
        GateRow("ctq_07_endpoint", "rational endpoint", "proved", "One exact endpoint comparison closes the pointwise budget.", str(exact["endpoint_check"]), "No floating-point inequality is used."),
        GateRow("ctq_08_move", "current budget", "proved", "The pointwise moving current costs less than one percent of the anchor.", str(exact["current_budget"]), "The stronger transfer budget remains separate."),
        GateRow("ctq_09_closed", "corrected ledger", "proved", "Six pointwise absolute packages have a corrected common budget.", str(exact["corrected_closed_budget"]), "The old mixed current/transfer fraction is superseded."),
        GateRow("ctq_10_four", "quadratic ownership", "proved", "The nonterminal quadratic current has exactly four normalized constituents.", str(exact["quadratic_constituents"]), "No terminal, affine, or moving-tail term is repeated."),
        GateRow("ctq_11_primitive", "determinant primitive", "proved", "The four constituents are one exact frequency derivative.", str(exact["quadratic_primitive"]), "This identity alone supplies no pointwise sign."),
        GateRow("ctq_12_phase", "phase polarization", "proved", "The primitive preserves both Hermitian and transpose phases.", str(exact["quadratic_polarization"]), "No real or imaginary component is dropped."),
        GateRow("ctq_13_transfer", "transfer telescope", "proved", "The integrated quadratic current is an endpoint determinant difference.", str(exact["quadratic_transfer"]), "The endpoint difference is not yet bounded."),
        GateRow("ctq_14_exclusion", "no-double-counting", "guard_validated", "Every neighbouring package remains excluded from E_quad.", str(exact["ownership"]), "The split epsilon=a+rho_obs is mandatory."),
        GateRow("ctq_15_reclass", "signed reclassification", "proved", "The unsuppressed quadratic belongs with the signed near-affine main.", str(exact["signed_reclassification"]), "The necessary and sufficient thresholds are not the same."),
        GateRow("ctq_16_handoff", "live theorem", "open", "The next theorem is a signed pointwise join or an endpoint-primitive transfer bound.", str(exact["live_obligation"]), "Neither alternative is proved here."),
        GateRow("ctq_17_boundary", "proof boundary", "guard_validated", "This correction and primitive do not prove the current sign or RH.", "No signed main inequality, endpoint determinant bound, all-q transport, contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion is claimed.", "The result is exact bookkeeping and route reduction."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    lines = [
        "# Newman C1 Current/Transfer Ownership and Quadratic Primitive Gate",
        "",
        "Date: 2026-08-05",
        "",
        "Status: exact ownership correction and determinant primitive; not a proof of the signed current or RH.",
        "",
        "## Ownership corrigendum",
        "",
        str(exact["ownership_correction"]),
        "",
        str(exact["prior_gate_boundary"]),
        "",
        "At the physical point,",
        "",
        "```text",
        str(exact["physical_current"]),
        str(exact["pointwise_bound"]),
        str(exact["current_budget"]),
        "```",
        "",
        "The corrected six-package current budget is",
        "",
        "```text",
        str(exact["corrected_closed_budget"]),
        "```",
        "",
        "## Quadratic primitive",
        "",
        "```text",
        str(exact["quadratic_constituents"]),
        str(exact["quadratic_primitive"]),
        str(exact["quadratic_transfer"]),
        "```",
        "",
        str(exact["quadratic_polarization"]),
        "",
        str(exact["ownership"]),
        "",
        "## Revised frontier",
        "",
        str(exact["signed_reclassification"]),
        "",
        str(exact["live_obligation"]),
        "",
        "This is not a proof of a complete-current sign, Lambda<=0, PF-infinity, RH, or a prize-level conclusion.",
        "",
        payload["success"],
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    source_audit = load_and_audit_sources()
    symbolic = symbolic_certificate()
    exact = exact_certificate()
    rows = build_rows(exact)
    payload = {
        "kind": KIND,
        "schema_version": 1,
        "date": "2026-08-05",
        "status": "current/transfer moving-tail ownership corrected; pointwise moving budget closed; quadratic determinant primitive proved; signed main open",
        "source_sha256": {key: value["sha256"] for key, value in source_audit.items()},
        "source_audit": source_audit,
        "symbolic": symbolic,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ownership_corrections": 1,
            "valid_transfer_budgets_retained": 1,
            "pointwise_moving_budget_denominator": 100,
            "quadratic_constituents": 4,
            "quadratic_primitives": 1,
            "corrected_closed_fraction": "289339/858000",
            "necessary_signed_fraction": "551501/858000",
            "sufficient_signed_fraction": "1130179/858000",
            "open_rows": sum(row.readiness == "open" for row in rows),
        },
        "next_action": str(exact["live_obligation"]),
        "pi_provenance": "No new pi is introduced. All pi factors remain those of the completed-zeta saddle, the Fourier character, kappa=1/(2pi i), and the existing anchor normalization.",
        "proof_boundary": "This gate corrects current/transfer ownership, proves a fixed-q1 pointwise moving-current budget, and identifies the exact nonterminal quadratic primitive. It proves no signed main inequality, endpoint determinant bound, complete-current sign, all-q transport, contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
        "success": "built Newman C1 current/transfer quadratic-primitive gate: 17 rows, 0 issues, pointwise |E_move^cur|<rho*A_T/100, 1 determinant primitive, 1 signed-main reclassification, 1 open theorem row",
    }
    atomic_write(RESULT_PATH, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(payload))
    print(payload["success"])


if __name__ == "__main__":
    main()
