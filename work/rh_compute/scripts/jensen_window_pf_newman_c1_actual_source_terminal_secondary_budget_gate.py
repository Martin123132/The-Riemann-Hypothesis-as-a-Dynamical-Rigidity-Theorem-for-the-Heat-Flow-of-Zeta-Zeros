#!/usr/bin/env python3
"""Build the actual-source terminal-secondary common-unit budget gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_terminal_secondary_budget_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "compensation_ledger": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_complete_compensation_ledger_gate.json",
    "raw_anchors": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_carrier_near_odd_harmonic_anchor_gate.json",
    "terminal_completion": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_terminal_bulk_phase_completion_gate.json",
    "relative_lift": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_terminal_nonterminal_relative_lift_gate.json",
    "source_phase": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_saddle_phase_lock_center_symmetry_gate.json",
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
    require(payloads["compensation_ledger"]["summary"]["sharp_adverse_thresholds"] == 1, "threshold drifted")
    require(payloads["raw_anchors"]["summary"]["explicit_raw_anchors"] == 2, "anchor drifted")
    require("14300" in payloads["raw_anchors"]["exact"]["transpose_split"]["lower_bound"], "transpose floor drifted")
    require(payloads["terminal_completion"]["counts"]["pure_terminal_rational_bounds"] == 1, "pure terminal drifted")
    require(payloads["terminal_completion"]["counts"]["pure_terminal_transpose_collapses"] == 1, "terminal collapse drifted")
    require(payloads["relative_lift"]["counts"]["lower_conditional_relative_traces"] == 2, "S1 source drifted")
    require(payloads["source_phase"]["summary"]["source_normalizer_identifications"] == 1, "source normalization drifted")

    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def exact_certificate() -> dict[str, str | int]:
    rho_log_coefficient = Fraction(1, 6)
    term_prefactor = (
        2
        * rho_log_coefficient
        * 2
        * Fraction(1, 24)
        * Fraction(1, 14300)
    )
    require(term_prefactor == Fraction(1, 514800), "terminal survivor prefactor")
    term_ratio = Fraction(1, 68640)
    require(term_prefactor < term_ratio, "terminal survivor target")

    aa_prefactor = (
        Fraction(10, 20)
        * rho_log_coefficient
        * Fraction(1, 14300)
    )
    require(aa_prefactor == Fraction(1, 171600), "pure terminal prefactor")
    aa_ratio = Fraction(1, 22880)
    require(aa_prefactor < aa_ratio, "pure terminal target")

    combined = term_ratio + aa_ratio
    require(combined == Fraction(1, 17160), "combined terminal ratio")
    require(59 * Fraction(1, 72_000_000_000) ** 3 < 1, "endpoint decay factor")

    return {
        "terminal_carrier": "On q=1, t=1/(2L^2), sigma>1/2, log N<(L+1)/2, and log N>L/2-1 imply S(log N)<0 and |E_N|<1 for L>=50.",
        "rho": "For the terminal kernel S_(1,N)=S_1^[N](alpha_P/N), the logarithmic bound (11.184.11), N+1<a+1<2a, log a<(L+1)/2, and log 36<4 give 0<S_(1,N)<L+9. Thus rho=|kappa|S_(1,N)|E_N|<(L+9)/6.",
        "decay_factor": "Since h<exp(-L/2), (L+9)h^3 is bounded by the decreasing function (L+9)exp(-3L/2); the L=50 rational check 59/(72000000000)^3<1 proves (L+9)h^3<1.",
        "survivor_identity": "W_T^term=-iE_N u_N u_(N,x)S_(1,N)/pi and S_(1,N)|E_N|/pi=2rho, so |E_term|<=2rho^2u_Nu_(N,x).",
        "survivor_budget": "Using u_N<2h, u_(N,x)<h^2/24, rho<(L+9)/6, A_T>14300, and (L+9)h^3<1 gives |E_term|<rho A_T/68640.",
        "source_relation": "At xi=Omega, w_N=nu E_N and S_(1,N)|E_N|=2pi rho.",
        "pure_identity": "The proved |R_aa|<h^3S_(1,N)^2|w_N|^2/160 gives |E_aa|<h^3pi^2rho^2/20 after E_aa=2R_aa/|nu|^2.",
        "pure_budget": "Using pi^2<10, rho<(L+9)/6, A_T>14300, and (L+9)h^3<1 gives |E_aa|<rho A_T/22880.",
        "combined_budget": "|E_term+E_aa|<rho A_T/17160.",
        "remaining_secondary": "E_sec^*=E_R+E_Delta+E_aff+E_move+E_quad contains exactly five still-unbounded secondary packages.",
        "finite_target": "If |E_sec^*|<=delta_*rho A_T, adverse-witness negativity forces E_F<-[98/100-delta_*-1/17160]rho A_T.",
        "survivor_prefactor_denominator": term_prefactor.denominator,
        "pure_prefactor_denominator": aa_prefactor.denominator,
        "term_denominator": term_ratio.denominator,
        "aa_denominator": aa_ratio.denominator,
        "combined_denominator": combined.denominator,
        "closed_terminal_packages": 2,
        "remaining_unbounded_secondary_packages": 5,
    }


def build_rows(exact: dict[str, str | int]) -> list[GateRow]:
    return [
        GateRow("tsb_01_carrier", "terminal carrier", "proved", "The source-free terminal carrier has modulus below one.", str(exact["terminal_carrier"]), "This uses the fixed physical q=1 chart."),
        GateRow("tsb_02_rho", "terminal multiplier", "proved", "The common terminal multiplier satisfies rho<(L+9)/6.", str(exact["rho"]), "This uses the logarithmic terminal kernel at alpha_P/N, not the lower kernel at alpha_P."),
        GateRow("tsb_03_survivor_identity", "ideal survivor", "proved", "The ideal transpose survivor reduces to 2rho^2u_Nu_(N,x) in modulus.", str(exact["survivor_identity"]), "This is the ideal terminal/nonterminal survivor, not the pure terminal row."),
        GateRow("tsb_04_survivor_budget", "survivor budget", "proved", "The ideal transpose survivor is bounded in rho A_T units.", str(exact["survivor_budget"]), "No other secondary package is included."),
        GateRow("tsb_05_source_relation", "source normalization", "proved", "The source-inclusive terminal carrier factors as nu times E_N.", str(exact["source_relation"]), "At the physical frequency c_xi=1."),
        GateRow("tsb_06_pure_identity", "pure terminal", "proved", "The pure terminal transpose envelope converts exactly to the common scalar units.", str(exact["pure_identity"]), "The Hermitian pure terminal real current is already zero."),
        GateRow("tsb_07_pure_budget", "pure budget", "proved", "The pure terminal transpose row is bounded in rho A_T units.", str(exact["pure_budget"]), "This is added once under the split convention."),
        GateRow("tsb_08_combined", "combined budget", "proved", "The two terminal secondary packages have a joint 1/17160 budget.", str(exact["combined_budget"]), "Triangle inequality is applied only after exact package ownership."),
        GateRow("tsb_09_ownership", "ownership guard", "proved", "The two bounded terminal packages are distinct and not double counted.", "E_term is the ideal transpose survivor; E_aa is the pure terminal quadratic row.", "The unsplit quadratic convention must not add E_aa again."),
        GateRow("tsb_10_remaining", "remaining secondary", "proved", "Exactly five secondary packages remain without common-unit bounds.", str(exact["remaining_secondary"]), "Remote, correction, affine, moving-tail, and quadratic packages remain open."),
        GateRow("tsb_11_transfer", "budget transfer", "proved", "The finite-package threshold improves by the certified terminal budget.", str(exact["finite_target"]), "A delta_* estimate is still missing."),
        GateRow("tsb_12_target", "live theorem", "open", "Bound the five-package secondary sum before attacking the signed finite package.", "Derive delta_* in common rho A_T units without one-sided tails or coefficient-only promotion.", "No complete-current sign or RH conclusion is supplied."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    lines = [
        "# Newman C1 Actual-Source Terminal Secondary Budget Gate",
        "",
        "Date: 2026-08-04",
        "",
        "Status: the ideal transpose survivor and pure terminal transpose row are jointly bounded by rho A_T/17160; five secondary packages and the complete-current sign remain open; not a proof of RH.",
        "",
        "## Terminal normalization",
        "",
        str(exact["terminal_carrier"]),
        "",
        str(exact["rho"]),
        "",
        str(exact["decay_factor"]),
        "",
        "## Common-unit bounds",
        "",
        str(exact["survivor_identity"]),
        "",
        str(exact["survivor_budget"]),
        "",
        str(exact["pure_identity"]),
        "",
        str(exact["pure_budget"]),
        "",
        str(exact["combined_budget"]),
        "",
        "## Reduced frontier",
        "",
        str(exact["remaining_secondary"]),
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
    require(len(rows) == 12, "row count")
    require(sum(row.readiness == "open" for row in rows) == 1, "open row count")

    success = (
        "built Newman C1 actual-source terminal-secondary budget gate: "
        "12 rows, 0 issues, rho<(L+9)/6, 2 terminal packages closed, "
        "joint ratio 1/17160, 5 unbounded secondary packages, "
        "1 live common-unit secondary target"
    )
    payload = {
        "kind": KIND,
        "schema_version": 1,
        "date": "2026-08-04",
        "status": "two terminal secondary packages bounded in anchor units; five secondary packages and complete-current sign open",
        "source_sha256": {key: item["sha256"] for key, item in source_audit.items()},
        "source_audit": source_audit,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "rho_upper_bound": "(L+9)/6",
            "terminal_survivor_prefactor_denominator": exact["survivor_prefactor_denominator"],
            "pure_terminal_prefactor_denominator": exact["pure_prefactor_denominator"],
            "terminal_survivor_denominator": exact["term_denominator"],
            "pure_terminal_denominator": exact["aa_denominator"],
            "combined_terminal_denominator": exact["combined_denominator"],
            "closed_terminal_packages": exact["closed_terminal_packages"],
            "remaining_unbounded_secondary_packages": exact["remaining_unbounded_secondary_packages"],
            "live_common_unit_secondary_targets": 1,
        },
        "next_action": "Bound E_sec^*=E_R+E_Delta+E_aff+E_move+E_quad in common rho A_T units. Start with the symmetrically paired remote numerator and retain its common cutoff; then join correction and moving-tail factors to their multiplying observation rows before using the h powers. If the resulting delta_* is below 98/100-1/17160, apply the exact finite-package threshold from this gate.",
        "pi_provenance": "The only pi factors are inherited through kappa=1/(2pi i), the terminal carrier, and the established anchor. The proof uses pi>3 and pi^2<10 only for rational upper bounds.",
        "proof_boundary": "This gate proves common-unit bounds for the ideal transpose terminal survivor and pure terminal transpose row, jointly below rho A_T/17160, using the logarithmic terminal-kernel bound, and sharpens the conditional finite-package target. It proves no paired-remote, correction, affine, moving-tail, or quadratic budget, signed finite carrier-near estimate, complete-current sign, all-q transport, contact exclusion, Q209 shell, cofinal successor, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
        "success": success,
    }

    atomic_write(RESULT_PATH, json.dumps(payload, indent=2) + "\n")
    atomic_write(NOTE_PATH, render_note(payload))
    print(success)


if __name__ == "__main__":
    main()
