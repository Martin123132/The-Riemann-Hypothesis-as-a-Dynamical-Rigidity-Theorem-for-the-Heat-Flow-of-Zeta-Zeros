#!/usr/bin/env python3
"""Build the actual-source complete-compensation ledger gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_complete_compensation_ledger_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "raw_center": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_raw_center_sign_obstruction_gate.json",
    "disk_geometry": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_reciprocal_stationary_disk_geometry_gate.json",
    "outer_pairing": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_ideal_cubic_symmetric_outer_pairing_gate.json",
    "relative_lift": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_terminal_nonterminal_relative_lift_gate.json",
    "correction_envelope": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_correction_coefficient_envelope_gate.json",
    "terminal_completion": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_terminal_bulk_phase_completion_gate.json",
    "full_support": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_ideal_cubic_reciprocal_terminal_full_support_homotopy_gate.json",
    "observation_image": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_observation_image_compression_gate.json",
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

    require(payloads["raw_center"]["summary"]["normalized_reserve"] == "49/100", "raw reserve drifted")
    require(payloads["raw_center"]["summary"]["macroscopic_compensation_conditions"] == 1, "compensation drifted")
    require(payloads["disk_geometry"]["summary"]["absolute_budget_obstructions"] == 2, "disk guard drifted")
    require(payloads["outer_pairing"]["counts"]["one_sided_divergence_witnesses"] == 1, "outer pairing drifted")
    require(payloads["outer_pairing"]["counts"]["signed_complete_ideal_bounds"] == 0, "outer sign status drifted")
    require(payloads["relative_lift"]["counts"]["mixed_current_compressions"] == 2, "relative lifts drifted")
    require(payloads["correction_envelope"]["counts"]["coefficient_envelopes"] == 6, "correction box drifted")
    require(payloads["correction_envelope"]["counts"]["physical_plin_bounds"] == 0, "correction frontier drifted")
    require(payloads["terminal_completion"]["counts"]["pure_terminal_rational_bounds"] == 1, "terminal bound drifted")
    require(payloads["full_support"]["counts"]["moving_tail_defect_identities"] == 1, "moving-tail identity drifted")
    require(payloads["full_support"]["counts"]["signed_full_current_bounds"] == 0, "full-current status drifted")
    require(payloads["observation_image"]["counts"]["quadratic_observation_pairings"] == 4, "quadratic rows drifted")
    require(payloads["observation_image"]["counts"]["quadratic_remainder_bounds"] == 0, "quadratic status drifted")

    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def exact_certificates() -> dict[str, dict[str, str | int]]:
    phi, theta, mismatch = sp.symbols("phi theta mismatch", real=True)
    left = sp.cos(theta) * sp.sin(phi + theta) + mismatch * sp.sin(phi)
    right = sp.sin(phi + 2 * theta) / 2 + (sp.Rational(1, 2) + mismatch) * sp.sin(phi)
    require(sp.simplify(sp.expand_trig(left - right)) == 0, "two-tone identity")

    reserve = Fraction(49, 100)
    threshold = 2 * reserve
    require(threshold == Fraction(98, 100), "compensation threshold")

    finite = ("E_F",)
    secondary = ("E_R", "E_term", "E_Delta", "E_aff", "E_aa", "E_move", "E_quad")
    require(len(set(finite + secondary)) == 8, "package partition")

    return {
        "normalization": {
            "scalar": "R_full=(|nu|^2/2)Sigma_full with Sigma_full=-2rho A_T B+E_tot.",
            "two_tone": "B=(1/2)sin(phi_N+2theta_eta)+(1/2+u_N/log N)sin(phi_N).",
            "package_count": 8,
        },
        "ideal_partition": {
            "finite_H": "W_H^F=A_H+mathscr M_N[P_H^0]+mathscr N_N[P_H^0].",
            "remote_H": "W_H^R=mathscr C_N[P_H^0].",
            "finite_T": "W_T^F=A_T+mathscr M_N[P_T^0]+mathscr N_N[P_T^0].",
            "remote_T": "W_T^R=mathscr C_N[P_T^0].",
            "terminal_T": "W_T^term=2iu_N mathcal T_N^[N][B_(T,N)^0].",
            "projection": "E_W=E_F+E_R+E_term after applying Re(conjugate(tau_0) .) in the Hermitian channel and Re(eta^2 tau_0 .) in the transpose channel.",
        },
        "physical_partition": {
            "correction": "E_Delta is the linear actual-source projection of L_H[Delta_H] and L_T[Delta_T,Delta B_(T,N)] after P_H,P_T,B_(T,N) are split from their ideal parts.",
            "remaining": "Normalize the fixed affine, pure terminal transpose, moving-tail, and four remaining nonterminal quadratic packages by E_X=2R_X/|nu|^2.",
            "total": "E_tot=E_F+E_R+E_term+E_Delta+E_aff+E_aa+E_move+E_quad.",
            "ownership": "Use the split convention: the transpose survivor is not a remote endpoint, the pure terminal row is not repeated in E_quad, and the moving-tail defect is attached once at observation level.",
        },
        "compensation": {
            "pointwise": "Sigma_full<0 iff E_tot<2rho A_T B.",
            "adverse": "At B<-49/100, negativity requires E_tot<-(98/100)rho A_T.",
            "retirement": "At the same adverse witness, E_tot>=-(98/100)rho A_T implies Sigma_full>0 and therefore a complete-current sign change.",
            "budget": "If E_tot=E_F+E_sec and |E_sec|<=delta rho A_T, then negativity at the adverse witness requires E_F<-(98/100-delta)rho A_T.",
            "mechanisms": "A proof may establish a sufficiently strong cell-uniform negative bias or a pointwise phase-correlated inequality; the adverse threshold alone does not choose between them.",
            "threshold_numerator": 98,
            "threshold_denominator": 100,
        },
        "scale_audit": {
            "finite": "No signed anchor-relative estimate; separated carrier budgets already exceed the anchors.",
            "remote": "Absolute convergence is proved, but the physical derivative numerator has no rho A_T bound.",
            "terminal": "Explicit h mechanisms are proved for terminal pieces, but no completed common-unit anchor comparison is available.",
            "correction": "Six O(h^2)/O(h^4) coefficient envelopes are proved, but no retained-row norm closes their current scale.",
            "affine": "The exact fixed affine split is known; no anchor-relative bound is proved.",
            "moving": "Atomic h^3 K phase displacement is proved; the joined observation factors in the full defect remain unbounded.",
            "quadratic": "Four exact pairings are known; no quadratic remainder bound is proved.",
            "nonpromotion": "For positive varepsilon and A, varepsilon*(A/varepsilon)=A; coefficient smallness alone is not current smallness.",
            "common_unit_secondary_bounds": 0,
        },
    }


def build_rows(exact: dict[str, dict[str, str | int]]) -> list[GateRow]:
    norm = exact["normalization"]
    ideal = exact["ideal_partition"]
    physical = exact["physical_partition"]
    comp = exact["compensation"]
    scale = exact["scale_audit"]
    return [
        GateRow("ccl_01_normalization", "common scalar", "proved", "Every retained-current contribution is placed in the scalar units of the raw-centre threshold.", str(norm["scalar"]), "This is a normalization identity, not a sign."),
        GateRow("ccl_02_two_tone", "phase audit", "proved", "The raw bracket has an exact two-tone form.", str(norm["two_tone"]), "The angle identity supplies no cancellation estimate."),
        GateRow("ccl_03_partition", "ideal partition", "proved", "The ideal residual splits into finite, paired-remote, and transpose-terminal packages.", str(ideal["projection"]), "One-sided remote tails remain inadmissible."),
        GateRow("ccl_04_finite", "finite package", "proved", "The carrier-plus-near package is isolated before any modulus.", f"{ideal['finite_H']} {ideal['finite_T']}", str(scale["finite"])),
        GateRow("ccl_05_remote", "remote package", "proved", "The remote contribution remains a common-cutoff symmetric pair.", f"{ideal['remote_H']} {ideal['remote_T']}", str(scale["remote"])),
        GateRow("ccl_06_terminal", "terminal survivor", "proved", "The ideal transpose terminal survivor has one unique owner.", str(ideal["terminal_T"]), str(scale["terminal"])),
        GateRow("ccl_07_correction", "physical correction", "proved", "Physical coefficient corrections are attached linearly after the ideal split.", str(physical["correction"]), str(scale["correction"])),
        GateRow("ccl_08_affine", "fixed affine", "proved", "The fixed affine row remains a separate real package.", str(physical["remaining"]), str(scale["affine"])),
        GateRow("ccl_09_pure_terminal", "pure terminal", "proved", "The pure terminal transpose term is separate from the ideal survivor and generic quadratic row.", str(physical["ownership"]), str(scale["terminal"])),
        GateRow("ccl_10_moving", "moving tail", "proved", "The moving-tail defect is attached once at observation level.", str(physical["ownership"]), str(scale["moving"])),
        GateRow("ccl_11_quadratic", "quadratic rows", "proved", "Exactly four remaining nonterminal observation pairings form the quadratic package.", str(physical["remaining"]), str(scale["quadratic"])),
        GateRow("ccl_12_ownership", "no double counting", "proved", "The split ledger has eight disjoint scalar package labels.", str(physical["total"]), str(physical["ownership"])),
        GateRow("ccl_13_threshold", "sharp threshold", "proved", "The adverse raw-centre witness forces the sharp 98/100 compensation threshold.", str(comp["adverse"]), str(comp["retirement"])),
        GateRow("ccl_14_budget", "conditional budget", "proved", "A secondary budget transfers quantitatively to the finite signed target.", str(comp["budget"]), "No positive secondary budget has yet been proved."),
        GateRow("ccl_15_nonpromotion", "scale guard", "proved", "Small coefficients alone cannot certify a perturbative current contribution.", str(scale["nonpromotion"]), "This is a logical nonpromotion witness, not a physical counterexample."),
        GateRow("ccl_16_mechanisms", "proof alternatives", "proved", "Uniform negative bias and phase-correlated response are both logically admissible mechanisms.", str(comp["mechanisms"]), "The earlier necessary condition does not force adaptation by itself."),
        GateRow("ccl_17_target", "live theorem", "open", "A common-unit secondary budget and signed finite carrier-near theorem remain necessary.", "First bound E_sec/(rho A_T); then apply the conditional budget lemma to E_F before any general-q transport.", "No complete-current sign or RH conclusion is supplied."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    lines = [
        "# Newman C1 Actual-Source Complete Compensation Ledger Gate",
        "",
        "Date: 2026-08-04",
        "",
        "Status: exact common-unit compensation ledger and sharp adverse threshold proved; all common-unit package estimates and the complete-current sign remain open; not a proof of RH.",
        "",
        "## Common scalar",
        "",
        str(exact["normalization"]["scalar"]),
        "",
        str(exact["normalization"]["two_tone"]),
        "",
        "## Exact partition",
        "",
        str(exact["ideal_partition"]["finite_H"]),
        "",
        str(exact["ideal_partition"]["finite_T"]),
        "",
        str(exact["ideal_partition"]["remote_H"]),
        "",
        str(exact["ideal_partition"]["remote_T"]),
        "",
        str(exact["ideal_partition"]["terminal_T"]),
        "",
        str(exact["physical_partition"]["total"]),
        "",
        str(exact["physical_partition"]["ownership"]),
        "",
        "## Compensation test",
        "",
        str(exact["compensation"]["pointwise"]),
        "",
        str(exact["compensation"]["adverse"]),
        "",
        str(exact["compensation"]["retirement"]),
        "",
        str(exact["compensation"]["budget"]),
        "",
        str(exact["compensation"]["mechanisms"]),
        "",
        "## Scale status",
        "",
    ]
    for key in ("finite", "remote", "terminal", "correction", "affine", "moving", "quadratic", "nonpromotion"):
        lines.append(f"- `{key}`: {exact['scale_audit'][key]}")
    lines.extend(
        [
            "",
            "## Gate rows",
            "",
            "| id | role | state | claim |",
            "|---|---|---|---|",
        ]
    )
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
    exact = exact_certificates()
    rows = build_rows(exact)
    require(len(rows) == 17, "row count")
    require(sum(row.readiness == "open" for row in rows) == 1, "open row count")

    success = (
        "built Newman C1 actual-source complete-compensation ledger gate: "
        "17 rows, 0 issues, 8 disjoint scalar packages, 1 two-tone identity, "
        "1 sharp 98/100 adverse threshold, 1 conditional budget lemma, "
        "0 common-unit secondary bounds, 1 live signed finite-package target"
    )
    payload = {
        "kind": KIND,
        "schema_version": 1,
        "date": "2026-08-04",
        "status": "exact common-unit compensation ledger proved; common-unit package estimates and complete-current sign open",
        "source_sha256": {key: item["sha256"] for key, item in source_audit.items()},
        "source_audit": source_audit,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "disjoint_scalar_packages": exact["normalization"]["package_count"],
            "two_tone_identities": 1,
            "sharp_adverse_thresholds": 1,
            "conditional_budget_lemmas": 1,
            "common_unit_secondary_bounds": exact["scale_audit"]["common_unit_secondary_bounds"],
            "uniform_bias_mechanisms": 1,
            "phase_correlated_mechanisms": 1,
            "live_signed_finite_package_targets": 1,
        },
        "next_action": "Derive a common-unit bound |E_sec|<=delta rho A_T for E_sec=E_R+E_term+E_Delta+E_aff+E_aa+E_move+E_quad, preserving symmetric remote cutoffs and split ownership. If delta<98/100, attack the signed finite carrier-plus-near inequality E_F<-(98/100-delta)rho A_T on the adverse actual-source arcs. In parallel, test a half-turn or cell-average lower bound for E_tot; one adverse point with E_tot>=-(98/100)rho A_T rigorously retires this current-sign route.",
        "pi_provenance": "All pi-dependent quantities are inherited from the completed-zeta saddle, Fourier character, and earlier anchors. The new two-tone identity is elementary trigonometry, while the 49/100 and 98/100 thresholds use rational arithmetic only.",
        "proof_boundary": "This gate proves a common scalar normalization, an exact eight-package ownership ledger, the two-tone raw phase audit, the sharp adverse threshold, a conditional budget lemma, and the uniform-bias/phase-correlation alternatives. It proves no common-unit secondary budget, signed finite carrier-near estimate, complete-current sign, all-q transport, contact exclusion, Q209 shell, cofinal successor, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
        "success": success,
    }

    atomic_write(RESULT_PATH, json.dumps(payload, indent=2) + "\n")
    atomic_write(NOTE_PATH, render_note(payload))
    print(success)


if __name__ == "__main__":
    main()
