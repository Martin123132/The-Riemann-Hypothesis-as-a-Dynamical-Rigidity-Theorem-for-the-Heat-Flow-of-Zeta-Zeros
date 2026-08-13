#!/usr/bin/env python3
"""Build the actual-source full-phase paired-remote budget gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_full_phase_remote_budget_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "compensation_ledger": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_complete_compensation_ledger_gate.json",
    "terminal_budget": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_terminal_secondary_budget_gate.json",
    "outer_pairing": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_ideal_cubic_symmetric_outer_pairing_gate.json",
    "relative_lift": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_terminal_nonterminal_relative_lift_gate.json",
    "raw_anchors": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_carrier_near_odd_harmonic_anchor_gate.json",
    "aggregate_scaling": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_aggregate_scaling_gate.json",
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

    require(payloads["compensation_ledger"]["summary"]["disjoint_scalar_packages"] == 8, "ledger drifted")
    require(payloads["terminal_budget"]["summary"]["combined_terminal_denominator"] == 17160, "terminal budget drifted")
    require(payloads["outer_pairing"]["counts"]["remote_pair_formula_checks"] == 5, "remote pairing drifted")
    require(payloads["outer_pairing"]["counts"]["one_sided_divergence_witnesses"] == 1, "cutoff guard drifted")
    require(payloads["relative_lift"]["counts"]["outer_remainder_components"] == 4, "IBP source drifted")
    require("14300" in payloads["raw_anchors"]["exact"]["transpose_split"]["lower_bound"], "anchor floor drifted")
    require(payloads["aggregate_scaling"]["counts"]["explicit_outside_roster_tail_bounds"] == 6, "physical scale source drifted")

    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def exact_certificate() -> dict[str, str | int]:
    x, k = sp.symbols("x k", positive=True, real=True)
    paired = 1 / (x - k) + 1 / (x + k)
    require(sp.simplify(paired - 2 * x / (x**2 - k**2)) == 0, "denominator pairing")

    endpoint_factor = Fraction(8, 27)
    endpoint_amplitude = Fraction(6503, 3000)
    anchor_conversion = Fraction(44, 175)
    endpoint_ratio = endpoint_factor * endpoint_amplitude * anchor_conversion
    remainder_ratio = Fraction(1, 14300)
    require(endpoint_ratio + remainder_ratio < Fraction(1, 6), "remote budget")

    known_secondary = Fraction(1, 6) + Fraction(1, 17160)
    require(known_secondary == Fraction(2861, 17160), "known secondary sum")

    return {
        "remote_phase": "For k>n_N=floor(2alpha), retain q_-=alpha/u-k and q_+=alpha/u+k; then k>2alpha, |q_-|>k/2, and q_+>k.",
        "denominator_pair": "1/q_-+1/q_+=2(alpha/u)/{(alpha/u)^2-k^2}, whose modulus is at most 8alpha/(3u k^2).",
        "tail_sum": "Since n_N>3alpha/2 and sum_(k>n_N)k^(-2)<1/n_N, one has alpha sum_(k>n_N)k^(-2)<2/3.",
        "first_boundary": "With |kappa|<1/6, both first endpoints contribute at most (8/27)(|A(1)|+|A(N)|).",
        "amplitude_box": "For R=(L+3)/2 and P=P_H^0 or P_T^0, |P|<R^3/2, |P'|<R^2, |P''|<2R, hence |A|<R^3/2, |A'|<R^3/u, and |A''|<2R^3/u^2.",
        "second_boundary": "The paired second-boundary coefficient is bounded by k^(-2){5|A'|+(9/2)|A|/u^2}, giving less than (29/2)R^3 after both endpoints.",
        "interior": "The paired C_r interior is bounded by k^(-2){5|A''|+(27/2)|A'|/u+(9/u^3+51/(4u^4))|A|}, whose integral is below (223/8)R^3.",
        "ibp_remainder": "After |kappa|^2<1/36 and the k^(-2) sum, one polynomial contributes below 339R^3/(288alpha); both together contribute below 339R^3/(144alpha)<1.",
        "endpoint_anchor": "The endpoint amplitudes obey c_H/c_T<7/6 and u_(N,x)d+2u_Nu_(N,x)<c_T/1000, while c_T/A_T=pi/D_N<44/175.",
        "rational_close": "(8/27)(6503/3000)(44/175)+1/14300<1/6.",
        "remote_budget": "|mathscr C_N[P_H^0]|+|mathscr C_N[P_T^0]|<A_T/6 and therefore |E_R|<rho A_T/6.",
        "reduced_secondary": "E_four=E_Delta+E_aff+E_move+E_quad contains the four secondary packages still lacking common-unit bounds.",
        "finite_target": "If |E_four|<=delta_4rho A_T, adverse-witness negativity forces E_F<-[98/100-delta_4-1/6-1/17160]rho A_T.",
        "endpoint_ratio_numerator": endpoint_ratio.numerator,
        "endpoint_ratio_denominator": endpoint_ratio.denominator,
        "known_secondary_numerator": known_secondary.numerator,
        "known_secondary_denominator": known_secondary.denominator,
        "closed_remote_packages": 1,
        "remaining_unbounded_secondary_packages": 4,
    }


def build_rows(exact: dict[str, str | int]) -> list[GateRow]:
    return [
        GateRow("frb_01_phase", "remote phase", "proved", "Every remote paired mode is uniformly nonstationary in both full phases.", str(exact["remote_phase"]), "This is the high symmetric remote tail only."),
        GateRow("frb_02_denominator", "endpoint pairing", "proved", "The first endpoint denominators cancel from order k^(-1) to k^(-2).", str(exact["denominator_pair"]), "The common cutoff is essential."),
        GateRow("frb_03_tail", "tail mass", "proved", "The paired denominator tail has alpha-weight below two thirds.", str(exact["tail_sum"]), "No one-sided infinite sum is formed."),
        GateRow("frb_04_first", "first boundary", "proved", "The first boundaries cost at most 8/27 of their endpoint amplitudes.", str(exact["first_boundary"]), "Both integer endpoints are retained."),
        GateRow("frb_05_amplitude", "amplitude box", "proved", "Both ideal cubics and two amplitude derivatives have one uniform R^3 box.", str(exact["amplitude_box"]), "The oscillatory alpha log u phase remains in the phase, not the amplitude."),
        GateRow("frb_06_second", "second boundary", "proved", "The paired second boundary is absolutely summable with an explicit constant.", str(exact["second_boundary"]), "This follows after the exact first pairing."),
        GateRow("frb_07_interior", "IBP interior", "proved", "The twofold full-phase interior is absolutely summable with an explicit constant.", str(exact["interior"]), "No differentiated oscillatory amplitude is used."),
        GateRow("frb_08_remainder", "IBP remainder", "proved", "Both second-boundary and interior remainders together are below one.", str(exact["ibp_remainder"]), "The physical exponential relation a^2>exp L supplies the scale."),
        GateRow("frb_09_anchor", "anchor conversion", "proved", "All ideal endpoint amplitudes convert to the transpose anchor with rational loss.", str(exact["endpoint_anchor"]), "The A_T>14300 floor absorbs the sub-unit IBP remainder."),
        GateRow("frb_10_close", "rational close", "proved", "The endpoint and remainder constants close below one sixth of A_T.", str(exact["rational_close"]), "This is an exact rational comparison."),
        GateRow("frb_11_remote", "remote budget", "proved", "The complete ideal paired-remote scalar is below rho A_T/6.", str(exact["remote_budget"]), "The result does not apply to physical correction polynomials."),
        GateRow("frb_12_cutoff", "cutoff guard", "proved", "The proof preserves the common symmetric cutoff throughout.", "The k^(-1) endpoint cancellation occurs before every modulus.", "Separate one-sided remote tails remain divergent proof objects."),
        GateRow("frb_13_reduced", "reduced frontier", "proved", "Only four secondary packages remain without common-unit bounds.", str(exact["reduced_secondary"]), "The finite package is not counted as secondary."),
        GateRow("frb_14_transfer", "finite target", "proved", "The adverse finite-package threshold now includes the remote and terminal budgets.", str(exact["finite_target"]), "A delta_4 estimate remains missing."),
        GateRow("frb_15_target", "live theorem", "open", "Bound the correction, affine, moving-tail, and quadratic sum in common units.", "Join every h-small coefficient to its multiplying observation row before taking a modulus.", "No complete-current sign or RH conclusion is supplied."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    lines = [
        "# Newman C1 Actual-Source Full-Phase Remote Budget Gate",
        "",
        "Date: 2026-08-04",
        "",
        "Status: the common-cutoff ideal paired-remote scalar is bounded by rho A_T/6; four secondary packages and the complete-current sign remain open; not a proof of RH.",
        "",
        "## Full-phase pairing",
        "",
        str(exact["remote_phase"]),
        "",
        str(exact["denominator_pair"]),
        "",
        str(exact["tail_sum"]),
        "",
        str(exact["first_boundary"]),
        "",
        "## Remainder control",
        "",
        str(exact["amplitude_box"]),
        "",
        str(exact["second_boundary"]),
        "",
        str(exact["interior"]),
        "",
        str(exact["ibp_remainder"]),
        "",
        "## Anchor comparison",
        "",
        str(exact["endpoint_anchor"]),
        "",
        str(exact["rational_close"]),
        "",
        str(exact["remote_budget"]),
        "",
        str(exact["finite_target"]),
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
    require(len(rows) == 15, "row count")
    require(sum(row.readiness == "open" for row in rows) == 1, "open row count")

    success = (
        "built Newman C1 actual-source full-phase remote budget gate: "
        "15 rows, 0 issues, 1 paired denominator cancellation, "
        "2 ideal remote channels closed, |E_R|<rho*A_T/6, "
        "4 unbounded secondary packages, 1 live common-unit four-package target"
    )
    payload = {
        "kind": KIND,
        "schema_version": 1,
        "date": "2026-08-04",
        "status": "ideal paired-remote package bounded in anchor units; four secondary packages and complete-current sign open",
        "source_sha256": {key: item["sha256"] for key, item in source_audit.items()},
        "source_audit": source_audit,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "paired_denominator_cancellations": 1,
            "closed_ideal_remote_channels": 2,
            "remote_budget_denominator": 6,
            "known_secondary_fraction": f"{exact['known_secondary_numerator']}/{exact['known_secondary_denominator']}",
            "remaining_unbounded_secondary_packages": exact["remaining_unbounded_secondary_packages"],
            "live_common_unit_four_package_targets": 1,
        },
        "next_action": "Bound E_four=E_Delta+E_aff+E_move+E_quad in common rho A_T units. The physical correction and moving-tail terms must be joined to their observation rows before exploiting h powers; the four quadratic pairings should remain signed. If |E_four|<=delta_4rho A_T with delta_4<98/100-1/6-1/17160, apply the sharpened finite-package target.",
        "pi_provenance": "All pi factors are inherited from the Fourier character, kappa=1/(2pi i), and the established anchor. The new constants come from exact paired denominators and rational inequalities.",
        "proof_boundary": "This gate proves a common-unit rho A_T/6 bound for the two ideal symmetrically paired remote channels and reduces the unbounded secondary frontier to four packages. It proves no physical-correction, fixed-affine, moving-tail, quadratic, signed finite carrier-near, complete-current, all-q, contact-exclusion, Q209 shell, cofinal-successor, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
        "success": success,
    }

    atomic_write(RESULT_PATH, json.dumps(payload, indent=2) + "\n")
    atomic_write(NOTE_PATH, render_note(payload))
    print(success)


if __name__ == "__main__":
    main()
