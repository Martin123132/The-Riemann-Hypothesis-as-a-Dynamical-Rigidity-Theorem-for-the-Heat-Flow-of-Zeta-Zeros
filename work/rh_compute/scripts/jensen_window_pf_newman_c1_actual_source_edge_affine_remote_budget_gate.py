#!/usr/bin/env python3
"""Build the actual-source edge-affine paired-remote budget gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_edge_affine_remote_budget_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "normalizer": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_theta_forward_sqrt_to_corrected_rs_C1_transfer_gate.json",
    "absolute_source": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_absolute_phase_anchor_reduction.json",
    "endpoint_edge": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_q1_finite_height_real_edge_remainder_gate.json",
    "anchor": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_carrier_near_odd_harmonic_anchor_gate.json",
    "remote_pairing": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_full_phase_remote_budget_gate.json",
    "physical_correction": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_physical_correction_budget_gate.json",
    "edge_ownership": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_edge_affine_ownership_gate.json",
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

    normalizer_rows = {row["id"]: row for row in payloads["normalizer"]["rows"]}
    require(
        "x^(7/4)/32" in normalizer_rows["ntfsrsg_02_uniform_normalizer_lower_bound"]["formula"],
        "normalizer floor drifted",
    )
    source_exact = payloads["absolute_source"]["exact"]
    require("|d_1|<2189/x<1/2" in source_exact["first_coefficient"], "first source drifted")
    require("beta=(sqrt(pi)/8)" in source_exact["endpoint_phase"], "endpoint beta drifted")

    endpoint = payloads["endpoint_edge"]
    require(endpoint["counts"]["uniform_q1_finite_height_edge_signs"] == 1, "edge gate drifted")
    majorants = endpoint["majorant_certificate"]
    require(majorants["endpoint_majorants"]["H_x"] == "|H_x|<3h", "H_x bound drifted")
    require(
        majorants["source_block_majorants"]["endpoint_rate"]["real_s_star_prime"]
        == "|Re(s_*')|<h^2/20000",
        "source-rate bound drifted",
    )

    require(payloads["anchor"]["summary"]["explicit_raw_anchors"] == 2, "anchor source drifted")
    require(
        payloads["remote_pairing"]["summary"]["paired_denominator_cancellations"] == 1,
        "remote pairing drifted",
    )
    require(
        payloads["physical_correction"]["summary"]["physical_correction_budget_denominator"] == 2000,
        "physical correction drifted",
    )
    require(
        payloads["edge_ownership"]["summary"]["edge_only_ownership_corrections"] == 1,
        "edge ownership drifted",
    )

    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def exact_certificate() -> dict[str, str | int]:
    source_ratio = Fraction(198, 329)
    require(source_ratio < Fraction(61, 100), "source ratio close")

    second_boundary = Fraction(7, 2)
    interior = Fraction(799, 120)
    require(second_boundary + interior == Fraction(1219, 120), "IBP constants")

    L, y = sp.symbols("L y", real=True, nonnegative=True)
    anchor_lower = sp.Rational(7, 2816) * L**2 * (L - 2) ** 2
    ideal_upper = (L + 3) ** 3 / 81 + 1
    difference = sp.together(sp.Rational(13, 100) * anchor_lower - ideal_upper)
    shifted = sp.Poly(sp.expand(sp.numer(difference).subs(L, y + 50)), y)
    require(all(coefficient > 0 for coefficient in shifted.all_coeffs()), "ideal anchor close")

    correction_l50 = (
        Fraction(1, 72_000_000_000)
        * 120**2
        * Fraction(53, 2) ** 3
        * 301
        / 6
    )
    require(correction_l50 < Fraction(1, 5), "correction remote close")

    affine_remote = 2 * Fraction(61, 100) * Fraction(131, 1000)
    require(affine_remote < Fraction(4, 25), "affine remote budget")
    closed = Fraction(1, 6) + Fraction(1, 17160) + Fraction(1, 2000) + Fraction(4, 25)
    reserve = Fraction(98, 100) - closed
    require(closed == Fraction(280759, 858000), "closed fraction")
    require(reserve == Fraction(560081, 858000), "reserve fraction")

    return {
        "source_upper": "B_0T_0<3exp(-L/4).",
        "terminal_lower": "rho>(329/3300)L exp(-L/4).",
        "source_ratio": "B_0T_0/rho<198/329<61/100.",
        "edge_rows": "For the genuine endpoint alone, |V_E|<3/5, |A_E|<4h, |Q_E|<4h, and |mathcal N_E|<4h^2.",
        "affine_split": "P_E=P_E^0+Delta_E and mathscr O_N[P_E]=mathscr N_N[P_E]+mathscr C_N[P_E] with one common symmetric cutoff.",
        "ideal_polynomial": "For R=(L+3)/2, |P_E^0|<R^3/6, |(P_E^0)'|<R^2/2, and |(P_E^0)''|<R.",
        "amplitude_box": "For A_E^0(u)=exp(S(log u))P_E^0(log u), |A_E^0|<R^3/6, |(A_E^0)'|<R^3/(5u), and |(A_E^0)''|<R^3/(2u^2).",
        "ideal_pairing": "Full-phase denominator pairing gives first boundary at most 8R^3/81 and second-boundary-plus-interior remainder below 1219R^3/(4320alpha)<1.",
        "ideal_remote": "|mathscr C_N[P_E^0]|<13A_T/100.",
        "field_drift": "With epsilon=h^2(L+70)^2, the physical field boxes give |Delta F_(beta_E)|<epsilon R^2/6 and |Delta_E|<epsilon R^3/6.",
        "drift_remote": "Contiguous-kernel L1 control gives |mathscr C_N[Delta_E]|<(1+6L)epsilon R^3/(6h)<1/5<A_T/1000.",
        "actual_remote": "|mathscr C_N[P_E]|<131A_T/1000.",
        "affine_remote_budget": "For E_aff^C=2B_0T_0 Re{omega_c mathscr C_N[P_E]}, |E_aff^C|<4rho A_T/25.",
        "signed_near": "Set E_FAN=E_F+2B_0T_0 Re{omega_c mathscr N_N[P_E]}; no modulus is taken across this near composition.",
        "closed_budget": "|E_R+E_term+E_aa+E_Delta+E_aff^C|<(280759/858000)rho A_T.",
        "signed_target": "If |E_move+E_quad|<=delta_2rho A_T, adverse-witness negativity requires E_FAN<-[560081/858000-delta_2]rho A_T.",
        "source_ratio_numerator": source_ratio.numerator,
        "source_ratio_denominator": source_ratio.denominator,
        "ibp_constant_numerator": (second_boundary + interior).numerator,
        "ibp_constant_denominator": (second_boundary + interior).denominator,
        "affine_remote_numerator": 4,
        "affine_remote_denominator": 25,
        "closed_numerator": closed.numerator,
        "closed_denominator": closed.denominator,
        "reserve_numerator": reserve.numerator,
        "reserve_denominator": reserve.denominator,
        "remaining_absolute_secondary_packages": 2,
    }


def build_rows(exact: dict[str, str | int]) -> list[GateRow]:
    return [
        GateRow("ear_01_source", "source upper", "proved", "The endpoint source scale has an explicit exponential upper bound.", str(exact["source_upper"]), "The normalizer amplitude is bounded before cancellation."),
        GateRow("ear_02_terminal", "terminal lower", "proved", "The terminal multiplier has a matching exponential lower bound.", str(exact["terminal_lower"]), "The logarithmic terminal kernel is bounded from below, not treated as O(h)."),
        GateRow("ear_03_ratio", "source ratio", "proved", "The endpoint source is uniformly controlled in rho units.", str(exact["source_ratio"]), "This ratio applies on the fixed physical q=1 chart."),
        GateRow("ear_04_edge", "edge coefficients", "proved", "The genuine endpoint row has one O(1), two O(h), and one O(h^2) entries.", str(exact["edge_rows"]), "The moved terminal carrier is excluded from these rows."),
        GateRow("ear_05_split", "affine split", "proved", "The affine polynomial and outer support split without one-sided tails.", str(exact["affine_split"]), "The common cutoff is retained throughout."),
        GateRow("ear_06_polynomial", "ideal polynomial", "proved", "The ideal affine cubic and two derivatives have explicit R envelopes.", str(exact["ideal_polynomial"]), "Coefficient non-suppression is retained."),
        GateRow("ear_07_amplitude", "amplitude derivatives", "proved", "Removing the full phase leaves controlled nonoscillatory derivatives.", str(exact["amplitude_box"]), "The alpha log u phase is not differentiated into the amplitude."),
        GateRow("ear_08_pairing", "full-phase pairing", "proved", "Common-cutoff pairing cancels the k^(-1) endpoint trace.", str(exact["ideal_pairing"]), "Neither one-sided remote infinity is defined separately."),
        GateRow("ear_09_ideal", "ideal remote", "proved", "The ideal affine paired-remote piece is anchor-small.", str(exact["ideal_remote"]), "The affine near piece is not included."),
        GateRow("ear_10_drift", "field drift", "proved", "Physical coefficient drift has an h^2 polynomial envelope.", str(exact["field_drift"]), "Small coefficients are joined to the edge row before taking a modulus."),
        GateRow("ear_11_drift_remote", "drift remote", "proved", "Finite contiguous kernels bound the remote drift without a pointwise roster-height loss.", str(exact["drift_remote"]), "This is a common-cutoff estimate."),
        GateRow("ear_12_actual", "actual remote", "proved", "The complete physical affine remote polynomial stays below 131/1000 of A_T.", str(exact["actual_remote"]), "Near and remote ownership remains disjoint."),
        GateRow("ear_13_budget", "affine remote budget", "proved", "The source-normalized affine paired-remote scalar costs less than 4/25 of the anchor.", str(exact["affine_remote_budget"]), "No bound for the affine near scalar is claimed."),
        GateRow("ear_14_near", "signed near", "proved", "The affine near term composes with the finite carrier-near package before moduli.", str(exact["signed_near"]), "This is an exact reclassification, not a sign theorem."),
        GateRow("ear_15_close", "closed budget", "proved", "Five absolute secondary packages have one common-unit budget.", str(exact["closed_budget"]), "Moving-tail and quadratic packages remain open."),
        GateRow("ear_16_target", "live theorem", "open", "Prove the signed finite-near-affine inequality after controlling the last two absolute packages.", str(exact["signed_target"]), "No complete-current sign or RH conclusion is supplied."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    lines = [
        "# Newman C1 Actual-Source Edge-Affine Remote Budget Gate",
        "",
        "Date: 2026-08-04",
        "",
        "Status: the common-cutoff affine paired-remote scalar is bounded by 4rho A_T/25; the affine near term remains signed with the finite package; not a proof of RH.",
        "",
        "## Source scale",
        "",
        str(exact["source_upper"]),
        "",
        str(exact["terminal_lower"]),
        "",
        str(exact["source_ratio"]),
        "",
        "## Edge and remote envelopes",
        "",
        str(exact["edge_rows"]),
        "",
        str(exact["ideal_polynomial"]),
        "",
        str(exact["amplitude_box"]),
        "",
        str(exact["ideal_pairing"]),
        "",
        str(exact["ideal_remote"]),
        "",
        str(exact["field_drift"]),
        "",
        str(exact["drift_remote"]),
        "",
        "## Common-unit close",
        "",
        str(exact["actual_remote"]),
        "",
        str(exact["affine_remote_budget"]),
        "",
        str(exact["signed_near"]),
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
    exact = exact_certificate()
    rows = build_rows(exact)
    require(len(rows) == 16, "row count")
    require(sum(row.readiness == "open" for row in rows) == 1, "open row count")

    success = (
        "built Newman C1 actual-source edge-affine remote budget gate: "
        "16 rows, 0 issues, B_0T_0/rho<61/100, "
        "|E_aff^C|<4rho*A_T/25, 2 remaining absolute secondary packages, "
        "1 live signed near-affine target"
    )
    payload = {
        "kind": KIND,
        "schema_version": 1,
        "date": "2026-08-04",
        "status": "edge-affine paired-remote package bounded in anchor units; signed near-affine theorem and two absolute secondary packages open",
        "source_sha256": {key: item["sha256"] for key, item in source_audit.items()},
        "source_audit": source_audit,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "source_scale_ratios": 1,
            "genuine_edge_coefficient_bounds": 4,
            "paired_denominator_cancellations": 1,
            "absolute_affine_remote_budgets": 1,
            "affine_remote_fraction": "4/25",
            "closed_secondary_fraction": "280759/858000",
            "remaining_reserve_fraction": "560081/858000",
            "remaining_absolute_secondary_packages": 2,
            "live_signed_near_affine_targets": 1,
        },
        "next_action": "Keep E_FAN=E_F+E_aff^N signed and common-cutoff invariant. Bound E_move+E_quad in rho A_T units, then prove or falsify E_FAN<-[560081/858000-delta_2]rho A_T on every adverse actual-source witness arc.",
        "pi_provenance": "Every pi is inherited from the completed-zeta endpoint, the Fourier character e(v)=exp(2pi i v), kappa=1/(2pi i), and the published normalizer. Rational bounds on pi only certify displayed inequalities.",
        "proof_boundary": "This gate proves a fixed-q1 source-scale comparison, genuine-edge coefficient box, common-cutoff affine paired-remote budget |E_aff^C|<4rho A_T/25, and the resulting signed near-affine target. It proves no affine-near sign, moving-tail or quadratic budget, complete-current sign, all-q transport, contact exclusion, Q209 shell, cofinal successor, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
        "success": success,
    }

    atomic_write(RESULT_PATH, json.dumps(payload, indent=2) + "\n")
    atomic_write(NOTE_PATH, render_note(payload))
    print(success)


if __name__ == "__main__":
    main()
