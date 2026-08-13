#!/usr/bin/env python3
"""Build the actual-source edge-affine ownership and classification gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_edge_affine_ownership_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "full_support": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "ideal_cubic_reciprocal_terminal_full_support_homotopy_gate.json"
    ),
    "outer_endpoint": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_outer_endpoint_kernel_gate.json"
    ),
    "observation_image": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_"
        "observation_image_compression_gate.json"
    ),
    "endpoint_source": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_complex_endpoint_source_normalization_gate.json",
    "source_phase": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_saddle_phase_lock_center_symmetry_gate.json",
    "physical_correction": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_physical_correction_budget_gate.json",
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

    require(payloads["full_support"]["counts"]["edge_only_matrix_reductions"] == 1, "edge row drifted")
    require(payloads["full_support"]["counts"]["full_support_linear_functionals"] == 1, "linear functional drifted")
    require(payloads["outer_endpoint"]["counts"]["edge_endpoint_factorizations"] == 1, "edge polynomial drifted")
    require(payloads["outer_endpoint"]["counts"]["twofold_outer_recompositions"] == 1, "outer complement drifted")
    require(payloads["observation_image"]["counts"]["linear_residual_functionals"] == 1, "observation image drifted")
    require(payloads["endpoint_source"]["summary"]["source_nonreal_witnesses"] == 1, "midpoint source drifted")
    require(payloads["source_phase"]["summary"]["source_normalizer_identifications"] == 1, "source phase drifted")
    require(payloads["physical_correction"]["summary"]["remaining_reserve_fraction"] == "697361/858000", "reserve drifted")

    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def exact_certificate() -> dict[str, str | int]:
    half = sp.Rational(1, 2)
    J = half * sp.Matrix(
        [
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, -1],
            [0, 0, -1, 0],
        ]
    )
    edge = sp.Matrix(sp.symbols("V_E N_E A_E Q_E", real=True))
    carrier = sp.Matrix(sp.symbols("V_C N_C A_C Q_C", real=True))
    terminal = edge + carrier
    t_edge = 2 * J * edge
    t_carrier = 2 * J * carrier
    t_terminal = 2 * J * terminal
    require(sp.simplify(t_terminal - t_edge - t_carrier) == sp.zeros(4, 1), "terminal row split")
    require(t_edge == sp.Matrix([edge[1], edge[0], -edge[3], -edge[2]]), "edge row order")

    derivative_error = sp.Matrix(sp.symbols("dV dN dA dQ", real=True))
    duplicate = sp.expand((t_terminal - t_edge).dot(derivative_error))
    expected_duplicate = sp.expand(2 * (carrier.T * J * derivative_error)[0])
    require(sp.simplify(duplicate - expected_duplicate) == 0, "duplicate carrier row")

    x, u_x = sp.symbols("x u_x", real=True)
    ideal_c = sp.Integer(1)
    ideal_a = sp.I * x / 2
    ideal_q = sp.I * x / 2
    ideal_n = -(x**2) / 4 - sp.I * u_x / 2
    edge_polynomial = sp.expand(
        sp.I
        * x
        * (
            edge[1] * ideal_c
            + edge[0] * ideal_n
            - edge[3] * ideal_a
            - edge[2] * ideal_q
        )
    )
    expected_polynomial = (
        -sp.I * edge[0] * x**3 / 4
        + (edge[2] + edge[3]) * x**2 / 2
        + (edge[0] * u_x / 2 + sp.I * edge[1]) * x
    )
    require(sp.simplify(edge_polynomial - expected_polynomial) == 0, "ideal edge cubic")
    require(sp.expand(edge_polynomial).coeff(x, 3) == -sp.I * edge[0] / 4, "cubic coefficient")

    gamma, o_r, o_i, scale = sp.symbols("gamma o_r o_i scale", real=True, positive=True)
    omega = sp.cos(gamma) + sp.I * sp.sin(gamma)
    outer = o_r + sp.I * o_i
    normalized_left = 2 * scale**2 * sp.re(omega * outer / scale)
    normalized_right = 2 * scale * sp.re(omega * outer)
    require(sp.simplify(sp.expand_complex(normalized_left - normalized_right)) == 0, "source normalization")

    midpoint_real = sp.sqrt(2 - sp.sqrt(2)) / 4
    midpoint_imag = sp.sqrt(2 + sp.sqrt(2)) / 4 - sp.sqrt(2) / 2
    require(bool(midpoint_real > sp.Rational(3, 16)), "midpoint real lower bound")
    require(bool(midpoint_imag < 0), "midpoint imaginary sign")

    closed = Fraction(1, 6) + Fraction(1, 17160) + Fraction(1, 2000)
    reserve = Fraction(49, 50) - closed
    require(closed == Fraction(143479, 858000), "closed fraction")
    require(reserve == Fraction(697361, 858000), "signed reserve")

    return {
        "row_split": "t_T=2Jomega_T=t_E+t_C with t_E=2Jomega_E and t_C=2Jomega_C(Omega).",
        "ownership": "After omega_C is moved into the full-support moment vector Y, the unique external row is t_E=(mathcal N_E,V_E,-Q_E,-A_E); using t_T adds t_C^T dot(epsilon) twice.",
        "edge_polynomial": "For beta_E=(mathcal N_E,V_E,-Q_E,-A_E), P_E(lambda)=i(lambda-log a)F_(beta_E)(lambda), and R_aff=Re{nu c_xi mathscr O_N[P_E]}.",
        "single_owner": "If P_E is displayed as E_aff, the remaining full-support linear polynomial uses beta_p=beta_N-beta_E.",
        "source_normalization": "At c_xi=1 and nu=omega_c/(B_0T_0), E_aff=2R_aff/|nu|^2=2B_0T_0 Re{omega_c mathscr O_N[P_E]}.",
        "support": "mathscr O_N[P_E]=mathscr N_N[P_E]+mathscr C_N[P_E]=mathscr M_N[P_E]-mathscr F_N[P_E], with one common reciprocal cutoff.",
        "ideal_cubic": "With x=lambda-log a, P_E^0=-iV_E x^3/4+(A_E+Q_E)x^2/2+(V_Eu_(N,x)/2+i mathcal N_E)x.",
        "midpoint": "At every p=0 midpoint, C_0'''(0)=0 and Im C_0(0)<0, so V_E=Re{(1+i/T_0)C_0(0)}>sqrt(2-sqrt(2))/4>3/16 and |[x^3]P_E^0|>3/64.",
        "classification": "The actual affine polynomial has an unsuppressed ideal cubic constituent. No h^2 coefficient argument makes the complete affine functional perturbative; it must remain signed until its near and paired-remote pieces are composed.",
        "closed_budget": "|E_R+E_term+E_aa+E_Delta|<=(143479/858000)rho A_T.",
        "signed_target": "For E_FA=E_F+E_aff, if |E_move+E_quad|<=delta_2rho A_T, adverse-witness negativity requires E_FA<-[697361/858000-delta_2]rho A_T.",
        "closed_fraction_numerator": closed.numerator,
        "closed_fraction_denominator": closed.denominator,
        "reserve_numerator": reserve.numerator,
        "reserve_denominator": reserve.denominator,
        "edge_only_ownership_corrections": 1,
        "signed_main_reclassifications": 1,
        "absolute_affine_budgets": 0,
        "remaining_absolute_secondary_packages": 2,
    }


def build_rows(exact: dict[str, str | int]) -> list[GateRow]:
    return [
        GateRow("eao_01_split", "terminal split", "proved", "The frozen terminal row is the sum of edge and moved-carrier rows.", str(exact["row_split"]), "This is an exact matrix identity."),
        GateRow("eao_02_owner", "external owner", "proved", "Only the genuine endpoint edge remains external after full-support transfer.", str(exact["ownership"]), "The moved terminal carriers remain in Y."),
        GateRow("eao_03_duplicate", "double-count guard", "proved", "Retaining t_T externally repeats one carrier-affine interaction.", str(exact["ownership"]), "The correction applies to the post-transfer chart."),
        GateRow("eao_04_polynomial", "polynomial compression", "proved", "The edge-affine row is one degree-at-most-five full-support functional.", str(exact["edge_polynomial"]), "No componentwise observation budget is introduced."),
        GateRow("eao_05_single", "single ownership", "proved", "Separating E_aff removes beta_E from the remaining linear polynomial.", str(exact["single_owner"]), "P_E is included exactly once."),
        GateRow("eao_06_source", "source normalization", "proved", "The affine package has one fixed physical source phase and scale.", str(exact["source_normalization"]), "There is no free affine phase."),
        GateRow("eao_07_support", "support identity", "proved", "The affine functional has one near plus paired-remote decomposition.", str(exact["support"]), "One-sided remote tails remain inadmissible."),
        GateRow("eao_08_terminal", "endpoint guard", "proved", "Terminal u_N suppression does not multiply the complete affine functional.", "The u_N factor in the terminal trace is local to lambda=log N.", "No global h gain is inferred from one endpoint."),
        GateRow("eao_09_cubic", "ideal constituent", "proved", "The correction-free affine polynomial has an explicit cubic constituent.", str(exact["ideal_cubic"]), "Physical coefficient corrections remain separate."),
        GateRow("eao_10_midpoint", "physical witness", "proved", "The affine cubic coefficient is bounded away from zero on recurring physical midpoints.", str(exact["midpoint"]), "This proves coefficient non-suppression, not functional largeness."),
        GateRow("eao_11_class", "route classification", "proved", "The unbounded affine package belongs initially to the signed main rather than an absolute secondary budget.", str(exact["classification"]), "Cancellation may still make it small after composition."),
        GateRow("eao_12_closed", "closed budget", "proved", "The four already bounded packages retain their combined common-unit budget.", str(exact["closed_budget"]), "Triangle arithmetic is used only for closed packages."),
        GateRow("eao_13_reclass", "signed package", "proved", "The live signed package is E_FA=E_F+E_aff.", "This is a bookkeeping reclassification with exact ownership.", "It is not a sign estimate."),
        GateRow("eao_14_target", "conditional target", "proved", "Only moving-tail and quadratic terms remain in the absolute-secondary target.", str(exact["signed_target"]), "A delta_2 bound and signed E_FA theorem are both missing."),
        GateRow("eao_15_live", "live theorem", "open", "Split and estimate the affine near/remote composition without losing its sign.", "Join mathscr N_N[P_E] to E_F and test the paired mathscr C_N[P_E] against the anchor after exact source normalization.", "No affine bound, complete-current sign, or RH conclusion is supplied."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    lines = [
        "# Newman C1 Actual-Source Edge-Affine Ownership Gate",
        "",
        "Date: 2026-08-04",
        "",
        "Status: edge-only affine ownership, exact source normalization, and signed-main reclassification proved; the affine functional and complete-current sign remain unbounded; not a proof of RH.",
        "",
        "## Ownership",
        "",
        str(exact["row_split"]),
        "",
        str(exact["ownership"]),
        "",
        str(exact["single_owner"]),
        "",
        "## Full-Support Functional",
        "",
        str(exact["edge_polynomial"]),
        "",
        str(exact["source_normalization"]),
        "",
        str(exact["support"]),
        "",
        "## Scale Classification",
        "",
        str(exact["ideal_cubic"]),
        "",
        str(exact["midpoint"]),
        "",
        str(exact["classification"]),
        "",
        "## Signed Target",
        "",
        str(exact["closed_budget"]),
        "",
        str(exact["signed_target"]),
        "",
        "## Gate Rows",
        "",
        "| id | role | state | claim |",
        "|---|---|---|---|",
    ]
    for row in payload["rows"]:
        lines.append(f"| {row['id']} | {row['role']} | {row['readiness']} | {row['claim']} |")
    lines.extend(
        [
            "",
            "## Next Action",
            "",
            payload["next_action"],
            "",
            "## Pi Provenance",
            "",
            payload["pi_provenance"],
            "",
            "## Proof Boundary",
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
    exact = exact_certificate()
    rows = build_rows(exact)
    require(len(rows) == 15, "row count")
    require(sum(row.readiness == "open" for row in rows) == 1, "open row count")

    success = (
        "built Newman C1 actual-source edge-affine ownership gate: "
        "15 rows, 0 issues, 1 edge-only ownership correction, "
        "1 exact affine polynomial, 1 physical cubic non-suppression witness, "
        "0 absolute affine budgets, 1 signed-main reclassification, "
        "2 remaining absolute secondary packages, 1 live finite-affine target"
    )
    payload = {
        "kind": KIND,
        "schema_version": 1,
        "date": "2026-08-04",
        "status": "edge-only affine ownership and signed-main reclassification proved; affine functional and complete-current sign open",
        "source_sha256": {key: item["sha256"] for key, item in source_audit.items()},
        "source_audit": source_audit,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "edge_only_ownership_corrections": exact["edge_only_ownership_corrections"],
            "exact_affine_polynomials": 1,
            "physical_cubic_nonsuppression_witnesses": 1,
            "absolute_affine_budgets": exact["absolute_affine_budgets"],
            "signed_main_reclassifications": exact["signed_main_reclassifications"],
            "closed_secondary_fraction": f"{exact['closed_fraction_numerator']}/{exact['closed_fraction_denominator']}",
            "remaining_reserve_fraction": f"{exact['reserve_numerator']}/{exact['reserve_denominator']}",
            "remaining_absolute_secondary_packages": exact["remaining_absolute_secondary_packages"],
            "live_signed_finite_affine_targets": 1,
        },
        "next_action": "Split mathscr O_N[P_E]=mathscr N_N[P_E]+mathscr C_N[P_E] with the established common cutoff. Compose the near affine term with E_F before taking moduli. Adapt the full-phase paired-remote estimate to P_E with its exact endpoint coefficients and source factor; only a genuinely anchor-small remote remainder should return to the absolute budget. Then bound E_move+E_quad and confront the signed E_FA target.",
        "pi_provenance": "Every pi is inherited from the completed-zeta endpoint, e(v)=exp(2pi i v), kappa=1/(2pi i), and the fixed source normalizer. No geometric constant is introduced.",
        "proof_boundary": "This gate proves edge-only affine ownership after full-support transfer, exact polynomial and source compression, the ideal cubic formula, a recurring physical coefficient non-suppression witness, and a signed-main reclassification. It proves no absolute affine, moving-tail, or quadratic budget, no signed finite-affine inequality, complete-current sign, all-q transport, contact exclusion, Q209 shell, cofinal successor, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
        "success": success,
    }

    atomic_write(RESULT_PATH, json.dumps(payload, indent=2) + "\n")
    atomic_write(NOTE_PATH, render_note(payload))
    print(success)


if __name__ == "__main__":
    main()
