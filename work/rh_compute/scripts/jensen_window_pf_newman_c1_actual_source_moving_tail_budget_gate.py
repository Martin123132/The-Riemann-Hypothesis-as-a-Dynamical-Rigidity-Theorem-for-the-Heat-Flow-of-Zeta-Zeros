#!/usr/bin/env python3
"""Build the actual-source moving-tail common-unit budget gate."""

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
KIND = "jensen_window_pf_newman_c1_actual_source_moving_tail_budget_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "homotopy": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_ideal_cubic_reciprocal_terminal_full_support_homotopy_gate.json",
    "growing_prefix": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_growing_prefix_finite_height_gate.json",
    "terminal_budget": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_terminal_secondary_budget_gate.json",
    "physical_correction": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_physical_correction_budget_gate.json",
    "joined_coefficients": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_joined_coefficient_envelope_gate.json",
    "edge_remote": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_edge_affine_remote_budget_gate.json",
    "anchor": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_carrier_near_odd_harmonic_anchor_gate.json",
    "ledger": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_complete_compensation_ledger_gate.json",
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

    homotopy = payloads["homotopy"]
    require(homotopy["counts"]["all_carrier_homotopies"] == 1, "homotopy drifted")
    require(homotopy["counts"]["moving_tail_defect_identities"] == 1, "moving defect drifted")
    homotopy_rows = rows_by_id(homotopy)
    require(
        "omega_C(xi)-omega_C(Omega)" in homotopy_rows["rth_07_move_defect"]["certificate"],
        "moving defect formula drifted",
    )

    prefix = payloads["growing_prefix"]
    require("h*K^18<=1" in prefix["exact"]["prefix"], "prefix exponent drifted")
    require(prefix["exact_integer_budget"]["h_max"] == "1/72000000000", "height floor drifted")

    terminal = payloads["terminal_budget"]
    require("sigma>1/2" in terminal["exact"]["terminal_carrier"], "carrier exponent drifted")
    require("log a<(L+1)/2" in terminal["exact"]["rho"], "logarithmic range drifted")

    correction = payloads["physical_correction"]
    require("epsilon=h^2(L+70)^2" in correction["exact"]["field_box"], "field box drifted")
    require(correction["summary"]["physical_correction_budget_denominator"] == 2000, "correction gate drifted")

    joined_rows = rows_by_id(payloads["joined_coefficients"])
    require("|delta|<4h^2" in joined_rows["pjce_03_rate_box"]["certificate"], "delta box drifted")

    edge = payloads["edge_remote"]
    require(edge["exact"]["source_ratio"] == "B_0T_0/rho<198/329<61/100.", "source ratio drifted")
    require("|V_E|<3/5" in edge["exact"]["edge_rows"], "edge row drifted")
    require(edge["summary"]["closed_secondary_fraction"] == "280759/858000", "closed budget drifted")

    anchor = payloads["anchor"]
    require("A_T>" in anchor["exact"]["transpose_split"]["lower_bound"], "anchor floor drifted")
    require("14300" in anchor["exact"]["transpose_split"]["lower_bound"], "anchor constant drifted")

    ledger = payloads["ledger"]
    require(
        "moving-tail defect is attached once" in ledger["exact"]["physical_partition"]["ownership"],
        "moving ownership drifted",
    )
    require("E_X=2R_X/|nu|^2" in ledger["exact"]["physical_partition"]["remaining"], "normalization drifted")

    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def symbolic_certificate() -> dict[str, str | int]:
    h, L, R, K = sp.symbols("h L R K", positive=True, real=True)
    eps = sp.Rational(7, 8) * h**2 / L**2
    p_star = R**2 / 2
    s0 = sp.Rational(16, 7) * h ** sp.Rational(-1, 2)
    s1 = sp.Rational(16, 7) * R * h ** sp.Rational(-1, 2)
    t1 = 4 * h ** sp.Rational(3, 2) * K**2

    carrier = sp.simplify(8 * eps * p_star**2 * s0 * t1)
    displacement = sp.simplify(8 * eps**2 * p_star**2 * t1 * s1)
    edge = sp.simplify(sp.Rational(24, 5) * eps * p_star * t1)
    require(carrier == 16 * h**3 * R**4 * K**2 / L**2, "carrier reduction")
    require(displacement == 14 * h**5 * R**5 * K**2 / L**4, "displacement reduction")
    require(edge == sp.Rational(42, 5) * h ** sp.Rational(7, 2) * R**2 * K**2 / L**2, "edge reduction")

    x = sp.symbols("x0:4", real=True)
    y = sp.symbols("y0:4", real=True)
    xv = sp.Matrix(x)
    yv = sp.Matrix(y)
    j = sp.Rational(1, 2) * sp.Matrix(
        [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]]
    )
    require(sp.expand((xv.T * j * xv)[0]) == x[0] * x[1] - x[2] * x[3], "J quadratic")
    bilinear = sp.expand(2 * (xv.T * j * yv)[0])
    expected = x[0] * y[1] + x[1] * y[0] - x[2] * y[3] - x[3] * y[2]
    require(bilinear == expected, "J bilinear")

    return {
        "observation_order": "omega=(V,mathcal N,A,Q)^T and Q(omega)=omega^T J omega=Vmathcal N-AQ.",
        "moving_defect": "D_move(xi)=2Delta omega_C(xi)^T J dot(omega_B)(xi)+2[omega_E+omega_B(xi)+omega_C(xi)]^T J dot(omega_C)(xi).",
        "moving_package": "For varepsilon=T_0-Omega and xi_theta=T_0-theta varepsilon, R_move=-varepsilon integral_0^1 D_move(xi_theta)dtheta and E_move=2R_move/|nu|^2.",
        "bilinear_bound": "For real four-vectors, |2x^TJy|<=4||x||_infinity||y||_infinity.",
        "carrier_reduction": "8varepsilon P_*^2S_0T_1<16h^3R_*^4K^2/L^2.",
        "displacement_reduction": "8varepsilon^2P_*^2T_1S_1<14h^5R_*^5K^2/L^4.",
        "edge_reduction": "(24/5)varepsilon P_*(B_0T_0)T_1<(42/5)(B_0T_0)h^(7/2)R_*^2K^2/L^2.",
        "symbolic_audits": 5,
    }


def exact_certificate() -> dict[str, str | int]:
    h0 = Fraction(1, 72_000_000_000)
    l0 = Fraction(50)
    r0 = Fraction(53, 2)

    require(Fraction(4) < Fraction(120**2, 1000), "delta absorption")
    require(Fraction(1001, 2000) < Fraction(3, 5), "A envelope")
    q_upper = (r0 / 2 + Fraction(1, 1_000_000)) * Fraction(1001, 1000) + Fraction(1, 1_000_000)
    require(q_upper < Fraction(3, 5) * r0, "Q envelope")
    n_upper = r0 * Fraction(3, 5) * r0 / 2 + Fraction(1, 500_000) * Fraction(1001, 1000)
    require(n_upper < r0**2 / 2, "mathcal N envelope")
    require(Fraction(1001, 1000) < r0**2 / 2, "V envelope")

    exp25_lower = sum(Fraction(25) ** k / math.factorial(k) for k in range(61))
    require(exp25_lower > 72_000_000_000, "exp(25) rational floor")

    main_endpoint = Fraction(52800, 329 * 14300) * h0**2 * r0**4 / l0**3
    displacement_endpoint = Fraction(46200, 329 * 14300) * h0**4 * r0**5 / l0**5
    edge_endpoint = Fraction(2562, 500 * 14300) * h0**3 * r0**2 / l0**2
    endpoint_sum = main_endpoint + displacement_endpoint + edge_endpoint
    require(endpoint_sum < Fraction(1, 1_000_000), "moving budget endpoint")

    closed = Fraction(280759, 858000) + Fraction(1, 1_000_000)
    reserve = Fraction(98, 100) - closed
    require(closed == Fraction(140379929, 429000000), "closed budget")
    require(reserve == Fraction(280040071, 429000000), "remaining reserve")

    return {
        "ownership": "The all-carrier transfer differs from the frozen-terminal transfer only through the integral of D_move; this observation-level term is outside the later relative-lift and quadratic packages and is inserted exactly once.",
        "definition": "R_move=-varepsilon integral_0^1D_move(T_0-theta varepsilon)dtheta and E_move=2R_move/|nu|^2, where varepsilon=T_0-Omega.",
        "weight": "For q_n=exp[t(log n)^2/4-sigma log n], q=1 gives sigma>1/2 and t(log n)^2/4<1/8, hence q_n<(8/7)n^(-1/2).",
        "observation_box": "For every 1<=n<=N, |P_V(log n)|,|P_mathcalN(log n)|,|P_A(log n)|,|P_Q(log n)|<P_*:=R_*^2/2, R_*=(L+3)/2.",
        "all_carrier_mass": "S_0=sum_(n<=N)q_n<(16/7)h^(-1/2).",
        "all_carrier_moment": "S_1=sum_(n<=N)u_nq_n<(16/7)R_*h^(-1/2).",
        "terminal_moment": "For T_1=sum_(B<n<=N)u_nq_n, K^18<=a and K>=2 give n>=a/2, q_n<2h^(1/2), and T_1<4h^(3/2)K^2.",
        "observation_envelopes": "Uniformly on the auxiliary segment, ||omega_B+omega_C||_infinity<=|nu|P_*S_0, ||dot(omega_B)||_infinity<=|nu|P_*S_1, ||dot(omega_C)||_infinity<=|nu|P_*T_1, and ||Delta omega_C||_infinity<=varepsilon|nu|P_*T_1.",
        "pointwise_defect": "|D_move|<4varepsilon|nu|^2P_*^2T_1S_1+4|nu|^2P_*^2S_0T_1+(12/5)|nu|P_*T_1.",
        "integrated_bound": "|E_move|<8varepsilon P_*^2(varepsilon T_1S_1+S_0T_1)+(24/5)varepsilon P_*(B_0T_0)T_1.",
        "raw_terms": "Using varepsilon<7h^2/(8L^2), |E_move|<16h^3R_*^4K^2/L^2+14h^5R_*^5K^2/L^4+(42/5)(B_0T_0)h^(7/2)R_*^2K^2/L^2.",
        "normalized_terms": "Using K^2<=h^(-1/9), B_0T_0/rho<61/100, rho>(329/3300)Lexp(-L/4), A_T>14300, and exp(L/4)<h^(-1/2), the three ratios are below [52800/(329*14300)]h^(43/18)R_*^4/L^3, [46200/(329*14300)]h^(79/18)R_*^5/L^5, and [2562/(500*14300)]h^(61/18)R_*^2/L^2.",
        "monotonicity": "After h<exp(-L/2), the logarithmic derivatives are -43/36+4/(L+3)-3/L, -79/36+5/(L+3)-5/L, and -61/36+2/(L+3)-2/L; all are negative for L>=50.",
        "endpoint_check": "At L=50, exp(-25)<1/72000000000 follows from the 61-term positive Taylor floor for exp(25). Replacing powers 43/18,79/18,61/18 by 2,4,3 gives an exact rational sum below 1/1000000.",
        "moving_budget": "|E_move|<rho A_T/1000000.",
        "closed_budget": "|E_R+E_term+E_aa+E_Delta+E_aff^C+E_move|<(140379929/429000000)rho A_T.",
        "signed_target": "If |E_quad|<=delta_1rho A_T, adverse-witness negativity requires E_FAN<-[280040071/429000000-delta_1]rho A_T.",
        "endpoint_main_numerator": main_endpoint.numerator,
        "endpoint_main_denominator": main_endpoint.denominator,
        "endpoint_displacement_numerator": displacement_endpoint.numerator,
        "endpoint_displacement_denominator": displacement_endpoint.denominator,
        "endpoint_edge_numerator": edge_endpoint.numerator,
        "endpoint_edge_denominator": edge_endpoint.denominator,
        "closed_numerator": closed.numerator,
        "closed_denominator": closed.denominator,
        "reserve_numerator": reserve.numerator,
        "reserve_denominator": reserve.denominator,
        "remaining_absolute_secondary_packages": 1,
    }


def build_rows(exact: dict[str, str | int]) -> list[GateRow]:
    return [
        GateRow("mtb_01_ownership", "package ownership", "proved", "The moving-tail defect is a single observation-level transfer package.", str(exact["ownership"]), "No terminal or nonterminal quadratic row is repeated."),
        GateRow("mtb_02_definition", "exact definition", "proved", "The package has a fixed sign and scalar normalization.", str(exact["definition"]), "The absolute estimate does not depend on the displayed sign."),
        GateRow("mtb_03_weight", "source-free weight", "proved", "Every correction-free carrier weight has a square-root majorant.", str(exact["weight"]), "This is restricted to the physical q=1 chart."),
        GateRow("mtb_04_observation", "observation box", "proved", "All four carrier observation polynomials share one source-free envelope.", str(exact["observation_box"]), "The polynomial coefficients are fixed during the auxiliary rotation."),
        GateRow("mtb_05_mass", "all-carrier mass", "proved", "The complete carrier mass is O(h^(-1/2)).", str(exact["all_carrier_mass"]), "No oscillatory cancellation is used."),
        GateRow("mtb_06_moment", "all-carrier moment", "proved", "The complete first distance moment is O(R_*h^(-1/2)).", str(exact["all_carrier_moment"]), "The crude logarithmic radius remains explicit."),
        GateRow("mtb_07_terminal", "terminal moment", "proved", "The moving terminal derivative has an h^(3/2)K^2 envelope.", str(exact["terminal_moment"]), "The growing-prefix condition is used before summing."),
        GateRow("mtb_08_displacement", "terminal displacement", "proved", "The terminal observation displacement gains one more varepsilon.", str(exact["observation_envelopes"]), "This follows from the exact auxiliary phase derivative."),
        GateRow("mtb_09_bilinear", "J bilinear", "proved", "The four-coordinate current has a dimension-fixed sup-norm bound.", "|2x^TJy|<=4||x||_infinity||y||_infinity.", "No dimension-dependent hidden constant is used."),
        GateRow("mtb_10_pointwise", "pointwise defect", "proved", "Both moving-tail terms have complete pointwise source envelopes.", str(exact["pointwise_defect"]), "The genuine edge is kept separate from the carrier mass."),
        GateRow("mtb_11_integrated", "integrated defect", "proved", "Physical integration and source normalization give three explicit positive terms.", str(exact["integrated_bound"]), "The factor 1/|nu|=B_0T_0 appears only in the genuine-edge term."),
        GateRow("mtb_12_main", "carrier ratio", "proved", "The leading all-carrier/terminal interaction is anchor-negligible.", str(exact["normalized_terms"]), "The rho lower bound is used only after the source factors cancel."),
        GateRow("mtb_13_delta", "displacement ratio", "proved", "The terminal-displacement interaction is anchor-negligible.", str(exact["monotonicity"]), "The extra varepsilon is retained."),
        GateRow("mtb_14_edge", "edge ratio", "proved", "The genuine-edge/moving-terminal interaction is source-normalized and anchor-negligible.", str(exact["endpoint_check"]), "The edge coefficient is O(1), so B_0T_0/rho is essential."),
        GateRow("mtb_15_close", "moving budget", "proved", "The complete moving-tail package costs less than one millionth of the anchor.", str(exact["moving_budget"]), "No claim about E_quad is included."),
        GateRow("mtb_16_target", "live theorem", "open", "Only the nonterminal quadratic package remains in the absolute-secondary frontier.", str(exact["signed_target"]), "A complete-current sign or RH conclusion still requires the quadratic and signed-main theorems."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    lines = [
        "# Newman C1 Actual-Source Moving-Tail Budget Gate",
        "",
        "Date: 2026-08-04",
        "",
        "Status: the complete observation-level moving-tail package is below rho A_T/1000000; the nonterminal quadratic package and signed main remain open; not a proof of RH.",
        "",
        "## Ownership and definition",
        "",
        str(exact["ownership"]),
        "",
        str(exact["definition"]),
        "",
        "## Source-free envelopes",
        "",
        str(exact["weight"]),
        "",
        str(exact["observation_box"]),
        "",
        str(exact["all_carrier_mass"]),
        "",
        str(exact["all_carrier_moment"]),
        "",
        str(exact["terminal_moment"]),
        "",
        str(exact["observation_envelopes"]),
        "",
        "## Common-unit budget",
        "",
        str(exact["pointwise_defect"]),
        "",
        str(exact["integrated_bound"]),
        "",
        str(exact["raw_terms"]),
        "",
        str(exact["normalized_terms"]),
        "",
        str(exact["monotonicity"]),
        "",
        str(exact["endpoint_check"]),
        "",
        str(exact["moving_budget"]),
        "",
        str(exact["closed_budget"]),
        "",
        str(exact["signed_target"]),
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
    symbolic = symbolic_certificate()
    exact = exact_certificate()
    rows = build_rows(exact)
    require(len(rows) == 16, "row count")
    require(sum(row.readiness == "open" for row in rows) == 1, "open row count")

    success = (
        "built Newman C1 actual-source moving-tail budget gate: "
        "16 rows, 0 issues, |E_move|<rho*A_T/1000000, "
        "1 remaining absolute secondary package, 1 live signed near-affine target"
    )
    payload = {
        "kind": KIND,
        "schema_version": 1,
        "date": "2026-08-04",
        "status": "moving-tail package bounded in anchor units; nonterminal quadratic package and signed near-affine theorem open",
        "source_sha256": {key: item["sha256"] for key, item in source_audit.items()},
        "source_audit": source_audit,
        "symbolic": symbolic,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_moving_tail_definitions": 1,
            "source_free_observation_envelopes": 4,
            "source_free_moment_envelopes": 3,
            "monotone_normalized_terms": 3,
            "moving_tail_budget_denominator": 1_000_000,
            "closed_secondary_fraction": "140379929/429000000",
            "remaining_reserve_fraction": "280040071/429000000",
            "remaining_absolute_secondary_packages": 1,
            "live_signed_near_affine_targets": 1,
        },
        "next_action": "Retain the exact four nonterminal quadratic pairings and their Hermitian/transpose phase structure. Derive a common-unit E_quad bound without componentwise row-norm promotion; then prove or falsify E_FAN<-[280040071/429000000-delta_1]rho A_T on every adverse actual-source witness arc.",
        "pi_provenance": "No new geometric pi is introduced. The only pi-dependent rate used in the observation box is u_(N,x)=h^2/(8pi), inherited from the completed-zeta/Riemann-Siegel saddle coordinate; pi>3 supplies a rational upper bound.",
        "proof_boundary": "This gate proves a fixed-q1 source-free observation envelope and the common-unit budget |E_move|<rho A_T/1000000 for the complete observation-level moving-tail defect. It proves no nonterminal quadratic budget, signed finite-near-affine inequality, complete-current sign, all-q transport, contact exclusion, Q209 shell, cofinal successor, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
        "success": success,
    }

    atomic_write(RESULT_PATH, json.dumps(payload, indent=2) + "\n")
    atomic_write(NOTE_PATH, render_note(payload))
    print(success)


if __name__ == "__main__":
    main()
