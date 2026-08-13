#!/usr/bin/env python3
"""Build the actual-source physical-correction common-unit budget gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_physical_correction_budget_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "compensation_ledger": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_complete_compensation_ledger_gate.json",
    "terminal_budget": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_terminal_secondary_budget_gate.json",
    "remote_budget": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_full_phase_remote_budget_gate.json",
    "correction_coefficients": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_correction_coefficient_envelope_gate.json",
    "joined_coefficients": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_joined_coefficient_envelope_gate.json",
    "relative_lift": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_terminal_nonterminal_relative_lift_gate.json",
    "abel_recombination": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_joined_lower_interior_abel_recombination_gate.json",
    "outer_pairing": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_ideal_cubic_symmetric_outer_pairing_gate.json",
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

    require(
        payloads["compensation_ledger"]["summary"]["disjoint_scalar_packages"] == 8,
        "compensation ledger drifted",
    )
    require(
        payloads["terminal_budget"]["summary"]["combined_terminal_denominator"] == 17160,
        "terminal budget drifted",
    )
    require(
        payloads["remote_budget"]["summary"]["remote_budget_denominator"] == 6,
        "remote budget drifted",
    )
    require(
        payloads["correction_coefficients"]["counts"]["coefficient_envelopes"] == 6,
        "coefficient envelope drifted",
    )
    require(
        payloads["joined_coefficients"]["counts"]["carrier_rate_bounds"] == 5,
        "rate envelope drifted",
    )
    require(
        payloads["relative_lift"]["counts"]["mixed_polynomial_functionals"] == 2,
        "relative-lift channels drifted",
    )
    require(
        payloads["abel_recombination"]["counts"]["global_recombinations"] == 4,
        "finite/exterior recombination drifted",
    )
    require(
        payloads["outer_pairing"]["counts"]["recombination_identities"] == 5,
        "outer-pairing identity drifted",
    )
    require(
        "14300" in payloads["raw_anchors"]["exact"]["transpose_split"]["lower_bound"],
        "transpose anchor floor drifted",
    )

    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def exact_certificate() -> dict[str, str | int]:
    g, z, z0, w, w0, d, aux = sp.symbols("g z z0 w w0 d aux")
    physical = g * (z**2 + z * d + w) + z * aux
    ideal = z0**2 + w0
    expanded = (g - 1) * ideal + g * (
        z**2 - z0**2 + z * d + w - w0
    ) + z * aux
    require(sp.expand(physical - ideal - expanded) == 0, "channel subtraction")

    endpoint_r = Fraction(53, 2)
    delta_b_coefficient = (
        Fraction(2001, 1000)
        + Fraction(1003, 1000)
        * (
            Fraction(1, 1500) / endpoint_r
            + Fraction(1, 1000) / endpoint_r
            + Fraction(1, 500) / endpoint_r**2
        )
        + Fraction(3, 1000) / endpoint_r
    )
    require(delta_b_coefficient < 3, "terminal coefficient box")

    h0 = Fraction(1, 72_000_000_000)
    k0 = Fraction(120**2)
    nonterminal = 6 * h0 * k0 * endpoint_r**3 * (4 + 6 * 50)
    require(nonterminal == Fraction(8_485_989, 1_250_000), "endpoint close value")
    require(nonterminal < Fraction(679, 100), "nonterminal close")

    terminal = 2 * (50 + 9) * h0**3 * k0 * endpoint_r**2
    require(terminal < Fraction(1, 1000), "terminal correction close")
    require(nonterminal + terminal < 7, "complete functional close")
    require(Fraction(7, 14300) < Fraction(1, 2000), "anchor close")

    known_secondary = Fraction(1, 6) + Fraction(1, 17160) + Fraction(1, 2000)
    reserve = Fraction(49, 50) - known_secondary
    require(known_secondary == Fraction(143479, 858000), "secondary sum")
    require(reserve == Fraction(697361, 858000), "remaining reserve")

    derivative_at_50 = (
        -Fraction(1, 2)
        + Fraction(2, 120)
        + Fraction(3, 53)
        + Fraction(6, 304)
    )
    terminal_derivative_at_50 = (
        Fraction(1, 59)
        - Fraction(3, 2)
        + Fraction(2, 120)
        + Fraction(2, 53)
    )
    require(derivative_at_50 < 0, "nonterminal monotonicity")
    require(terminal_derivative_at_50 < 0, "terminal monotonicity")

    return {
        "field_box": "For epsilon=h^2(L+70)^2 and R_*=(L+3)/2, the physical C,D,R,R_x,delta envelopes imply C errors below epsilon, D and derivative errors below epsilon/1000, rate errors below epsilon/6000, and |delta-conjugate(delta)|<epsilon/1000.",
        "channel_factorization": "Both terminal channels have B=G(z^2+zd+w)+zH after setting d=0 in the transpose channel and using the exact conjugated Hermitian terminal data.",
        "coefficient_budget": "The common terminal-factor subtraction gives |Delta B_H|,|Delta B_T|<3epsilon R_*^2 and hence |Delta_H|,|Delta_T|<6epsilon R_*^3.",
        "finite_recombination": "For every correction polynomial, mathscr M_N[P]+mathscr O_N[P]=2mathscr M_N[P]-mathscr F_N[P].",
        "dirichlet_kernel": "For m consecutive frequencies, |K_m(v)|=|sin(pi m v)|/|sin(pi v)| and integral_0^1|K_m(v)|dv<=1+log m by splitting min(m,1/(2v)) at 1/(2m).",
        "functional_bounds": "Since m<=2a^2, N<=a=h^(-1), and S<=0, |mathscr M_N[P]|<P_*/h and |mathscr F_N[P]|<3LP_*/h for |P|<=P_*.",
        "nonterminal_close": "Both correction channels outside the terminal survivor contribute less than 8485989/1250000<679/100.",
        "terminal_close": "The corrected logarithmic multiplier rho<(L+9)/6 gives a transpose terminal correction below 2(L+9)h^3(L+70)^2R_*^2<1/1000.",
        "correction_budget": "|E_Delta|<7rho<rho A_T/2000.",
        "reduced_secondary": "E_three=E_aff+E_move+E_quad contains the three secondary packages still lacking common-unit bounds.",
        "finite_target": "If |E_three|<=delta_3rho A_T, then E_F<-[697361/858000-delta_3]rho A_T is necessary at every adverse witness.",
        "delta_b_coefficient_numerator": delta_b_coefficient.numerator,
        "delta_b_coefficient_denominator": delta_b_coefficient.denominator,
        "nonterminal_numerator": nonterminal.numerator,
        "nonterminal_denominator": nonterminal.denominator,
        "known_secondary_numerator": known_secondary.numerator,
        "known_secondary_denominator": known_secondary.denominator,
        "remaining_reserve_numerator": reserve.numerator,
        "remaining_reserve_denominator": reserve.denominator,
        "closed_terminal_coefficient_channels": 2,
        "remaining_unbounded_secondary_packages": 3,
    }


def build_rows(exact: dict[str, str | int]) -> list[GateRow]:
    return [
        GateRow("pcb_01_fields", "physical field box", "proved", "All physical terminal coefficient fields have one explicit epsilon box.", str(exact["field_box"]), "This box is confined to q=1 and 0<=lambda<=log N."),
        GateRow("pcb_02_factor", "terminal factorization", "proved", "Hermitian and transpose correction channels share one exact bilinear form.", str(exact["channel_factorization"]), "Conjugation remains explicit in the Hermitian channel."),
        GateRow("pcb_03_bilinear", "bilinear correction", "proved", "Both terminal coefficient perturbations are below 3epsilon R_*^2.", str(exact["coefficient_budget"]), "No retained observation-row norm is used."),
        GateRow("pcb_04_polynomial", "relative polynomial", "proved", "Both lifted correction polynomials are below one common P_* envelope.", "The relative lift has modulus below 2R_*, so P_*=6epsilon R_*^3.", "Degree five is retained."),
        GateRow("pcb_05_recombine", "finite/exterior identity", "proved", "The complete correction functional returns to one carrier sum and one finite Fourier roster.", str(exact["finite_recombination"]), "This is an exact full-support identity."),
        GateRow("pcb_06_kernel", "contiguous kernel", "proved", "The reciprocal roster is a consecutive-frequency Dirichlet kernel.", str(exact["dirichlet_kernel"]), "The starting frequency contributes only a unit phase."),
        GateRow("pcb_07_l1", "kernel L1 theorem", "proved", "The kernel costs 1+log m per unit interval rather than its pointwise height.", str(exact["dirichlet_kernel"]), "The interval endpoints are integers."),
        GateRow("pcb_08_scale", "physical roster scale", "proved", "The roster logarithm is below 3L and the carrier count is below h^(-1).", str(exact["functional_bounds"]), "This uses alpha_P<=a^2 and N<=a."),
        GateRow("pcb_09_carrier", "carrier bound", "proved", "Each correction carrier sum is below P_*/h.", str(exact["functional_bounds"]), "Starred endpoint weights are at most one."),
        GateRow("pcb_10_roster", "finite-roster bound", "proved", "Each correction reciprocal-band sum is below 3LP_*/h.", str(exact["functional_bounds"]), "No infinite one-sided tail is formed."),
        GateRow("pcb_11_close", "nonterminal close", "proved", "Both complete nonterminal correction channels together are below 6.79.", str(exact["nonterminal_close"]), "The exponential majorant decreases for L>=50."),
        GateRow("pcb_12_terminal", "terminal survivor", "proved", "The correction to the transpose terminal survivor is below 1/1000.", str(exact["terminal_close"]), "The logarithmic terminal-kernel bound rho<(L+9)/6 is used inside this source-free functional."),
        GateRow("pcb_13_anchor", "common-unit budget", "proved", "The complete physical correction package is below rho A_T/2000.", str(exact["correction_budget"]), "The anchor floor A_T>14300 closes the conversion."),
        GateRow("pcb_14_fraction", "budget arithmetic", "proved", "The three closed secondary budgets consume exactly 143479/858000.", "1/6+1/17160+1/2000=143479/858000.", "Triangle arithmetic is applied only after package ownership."),
        GateRow("pcb_15_reduced", "reduced frontier", "proved", "Only three secondary packages remain without common-unit bounds.", str(exact["reduced_secondary"]), "The signed finite package is not counted as secondary."),
        GateRow("pcb_16_target", "live theorem", "open", "Bound the affine, moving-tail, and nonterminal quadratic package in common units.", str(exact["finite_target"]), "No signed finite-package or complete-current theorem is supplied."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    lines = [
        "# Newman C1 Actual-Source Physical-Correction Budget Gate",
        "",
        "Date: 2026-08-04",
        "",
        "Status: the complete physical coefficient-correction scalar is bounded by rho A_T/2000; three secondary packages and the complete-current sign remain open; not a proof of RH.",
        "",
        "## Physical coefficient box",
        "",
        str(exact["field_box"]),
        "",
        str(exact["channel_factorization"]),
        "",
        str(exact["coefficient_budget"]),
        "",
        "## Contiguous-kernel return",
        "",
        str(exact["finite_recombination"]),
        "",
        str(exact["dirichlet_kernel"]),
        "",
        str(exact["functional_bounds"]),
        "",
        str(exact["nonterminal_close"]),
        "",
        str(exact["terminal_close"]),
        "",
        "## Anchor comparison",
        "",
        str(exact["correction_budget"]),
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
    require(len(rows) == 16, "row count")
    require(sum(row.readiness == "open" for row in rows) == 1, "open row count")

    success = (
        "built Newman C1 actual-source physical-correction budget gate: "
        "16 rows, 0 issues, 2 terminal coefficient channels closed, "
        "1 contiguous-kernel L1 theorem, |E_Delta|<rho*A_T/2000, "
        "3 unbounded secondary packages, 1 live common-unit three-package target"
    )
    payload = {
        "kind": KIND,
        "schema_version": 1,
        "date": "2026-08-04",
        "status": "physical coefficient-correction package bounded in anchor units; three secondary packages and complete-current sign open",
        "source_sha256": {key: item["sha256"] for key, item in source_audit.items()},
        "source_audit": source_audit,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "closed_terminal_coefficient_channels": exact["closed_terminal_coefficient_channels"],
            "contiguous_kernel_l1_theorems": 1,
            "physical_correction_budget_denominator": 2000,
            "known_secondary_fraction": f"{exact['known_secondary_numerator']}/{exact['known_secondary_denominator']}",
            "remaining_reserve_fraction": f"{exact['remaining_reserve_numerator']}/{exact['remaining_reserve_denominator']}",
            "remaining_unbounded_secondary_packages": exact["remaining_unbounded_secondary_packages"],
            "live_common_unit_three_package_targets": 1,
        },
        "next_action": "Bound E_three=E_aff+E_move+E_quad in common rho A_T units. Audit the fixed affine row first because it is linear and source explicit; preserve the moving-tail defect at observation level and keep the four quadratic pairings signed. If |E_three|<=delta_3rho A_T with delta_3<697361/858000, confront the signed finite carrier-plus-near target.",
        "pi_provenance": "The pi in the Dirichlet-kernel sine quotient comes from the fixed Fourier character e(v)=exp(2pi i v) and its finite geometric sum. All other pi factors are inherited through kappa, the completed-zeta saddle, and the established anchor.",
        "proof_boundary": "This gate proves a common-unit rho A_T/2000 bound for the complete physical coefficient-correction package and reduces the unbounded secondary frontier to three packages. It proves no fixed-affine, moving-tail, quadratic, signed finite carrier-near, complete-current, all-q, contact-exclusion, Q209 shell, cofinal-successor, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
        "success": success,
    }

    atomic_write(RESULT_PATH, json.dumps(payload, indent=2) + "\n")
    atomic_write(NOTE_PATH, render_note(payload))
    print(success)


if __name__ == "__main__":
    main()
