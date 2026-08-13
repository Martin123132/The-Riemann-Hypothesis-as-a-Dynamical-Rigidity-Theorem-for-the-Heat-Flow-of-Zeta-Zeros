#!/usr/bin/env python3
"""Build the C1 reciprocal-stationary anchor-disk geometry gate."""

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
KIND = "jensen_window_pf_newman_c1_reciprocal_stationary_disk_geometry_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "anchor": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_carrier_near_odd_harmonic_anchor_gate.json",
    "finite_band": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_finite_band_"
        "dirichlet_cell_reduction_gate.json"
    ),
    "ridge": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_ideal_cubic_double_"
        "morse_ridge_carrier_inversion_gate.json"
    ),
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
    expanded = sp.expand(expression)
    if expanded != 0:
        raise RuntimeError(f"{label}: {expanded}")


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def load_sources() -> dict[str, dict]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        require(path.is_file(), f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))
    return payloads


def audit_sources(payloads: dict[str, dict]) -> dict[str, dict[str, str]]:
    require(
        payloads["anchor"]["summary"]["explicit_raw_anchors"] == 2,
        "anchor source drifted",
    )
    require(
        payloads["anchor"]["summary"]["open_stationary_residual_targets"] == 1,
        "anchor frontier drifted",
    )
    require(
        payloads["finite_band"]["counts"]["reciprocal_block_decompositions"]
        == 1,
        "finite-band source drifted",
    )
    require(
        payloads["ridge"]["counts"]["inversion_identities"] == 5,
        "ridge source drifted",
    )
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def algebra_certificate() -> dict[str, str | int]:
    a_h, a_t = sp.symbols("A_H A_T", positive=True, real=True)
    r_h, r_t, q_t = sp.symbols("R_H R_T Q_T")
    j_h = -a_h + r_h
    j_t = -a_t + r_t + q_t
    w_h = r_h
    w_t = r_t + q_t
    require_zero(w_h - (a_h + j_h), "Hermitian combined residual")
    require_zero(w_t - (a_t + j_t), "transpose combined residual")

    m, tau, exterior, alias = sp.symbols("M tau E Lambda")
    band_cell = tau + m - exterior + alias
    cell_join = 2 * m - band_cell
    require_zero(
        cell_join - (m - tau + exterior - alias),
        "per-cell reciprocal decomposition",
    )

    x, y, delta = sp.symbols("x y delta", real=True)
    disk_square = (a_h + x) ** 2 + y**2 - a_h**2
    require_zero(disk_square - (x**2 + y**2 + 2 * a_h * x), "disk square")
    contracted_square = (a_h - delta) ** 2 - a_h**2
    require_zero(contracted_square - (-2 * a_h * delta + delta**2), "delta disk")

    return {
        "hermitian_join": "mathcal J_H^0=-A_H+R_H, hence W_H:=R_H=A_H+mathcal J_H^0.",
        "transpose_join": "mathcal J_T^0=-A_T+R_T+Q_T, hence W_T:=R_T+Q_T=A_T+mathcal J_T^0.",
        "transpose_improvement": "The exact sufficient condition |R_T+Q_T|<=A_T-delta is weaker than the previous separated condition |R_T|+|Q_T|<=A_T-delta and still implies |mathcal J_T^0|>=delta.",
        "cell_join": "For S_q=tau_q+M_q-E_q^(ext)+Lambda_q^(alias), one has 2M_q-S_q=M_q-tau_q+E_q^(ext)-Lambda_q^(alias).",
        "global_join": "Summing cells gives mathcal J_H^0=mathscr M_N[P_H^0]-tau_band[P_H^0]+E_band[P_H^0], while mathcal J_T^0 has the same expression for P_T^0 plus the exact transpose terminal survivor.",
        "anchored_stationary": "Therefore W_H and W_T are A_H and A_T plus the complete tie-invariant reciprocal joins. Internal aliases, exterior pieces, and tie halves must recompose before a modulus is taken.",
        "disk_identity": "For A>0 and J=x+iy, |A+J|^2-A^2=|J|^2+2A Re(J).",
        "contracted_disk": "For 0<delta<=A, |A+J|<=A-delta iff |J|^2+2A Re(J)<=-2A delta+delta^2; the reverse triangle inequality then gives |J|>=delta.",
        "signed_necessity": "Any nonzero J in the closed tangent disk |A+J|<=A satisfies Re(J)<=-|J|^2/(2A)<0. A magnitude-only estimate for J cannot establish this directional condition.",
        "combined_residuals": 2,
        "optimal_disk_identities": 2,
    }


def budget_certificate() -> dict[str, str | int | float]:
    # The exact lower comparison at L=50 uses a positive Taylor truncation.
    exp_lower = sum(
        Fraction(25, 2) ** k / math.factorial(k)
        for k in range(41)
    )
    anchor_upper_50 = Fraction(57 * 51 * 53 * 53, 192)
    budget_lower_50 = Fraction(243, 520) * exp_lower
    require(budget_lower_50 > anchor_upper_50, "L=50 budget comparison")
    require(Fraction(4, 51) < Fraction(1, 4), "monotone comparison")

    # All constants in the source-specific interval argument are rational.
    count_factor = Fraction(1, 65)
    exponent_factor = Fraction(9, 10)
    square_root_factor = Fraction(5, 1)
    polynomial_factor = Fraction(27, 4)
    two_carriers = Fraction(2, 1)
    carrier_coefficient = (
        count_factor
        * exponent_factor
        * square_root_factor
        * polynomial_factor
        * two_carriers
    )
    require(carrier_coefficient == Fraction(243, 260), "carrier coefficient")

    return {
        "test_block": "Let I_N={q in Z:ceil(N/64)<=q<=floor(N/32)}. On L>=50, N>4160 and #I_N>=N/65; every such q is an interior carrier with full starred weight.",
        "amplitude_floor": "Since sigma<1/2+1/(8L), log(a)<(L+1)/2, and (L+1)/(16L)<=51/800<log(10/9), exp(S(log q))>9/(10 sqrt(q))>9/(2 sqrt(N)) on I_N.",
        "polynomial_floor": "For q in I_N, u_q-u_N=log(N/q)>=log 32>3 and u_q+u_N>=u_q-u_N. Thus |P_H^0(log q)|>27/4. The imaginary cubic part of P_T^0 has the same floor and is orthogonal to its real affine part, so |P_T^0(log q)|>27/4.",
        "carrier_budget": "For C=H,T define B_C=2 sum_(q=1)^N omega_N(q)|f_(P_C^0)(q)|. The preceding block gives B_C>(243/260)sqrt(N)>(243/520)exp(L/4).",
        "anchor_upper": "The floor definitions give D_N<1+log(2N+1)<(L+7)/2, while b<(L+1)/2 and d<(L+3)/2. Since pi>3, both A_H and A_T are below U(L)=(L+7)(L+1)(L+3)^2/192.",
        "uniform_comparison": "At L=50 the 41-term positive Taylor lower sum for exp(25/2) gives (243/520)exp(25/2)>125396 while U(50)<42531. The logarithmic derivative of U is below 4/(L+1)<=4/51<1/4, so the strict budget gap increases for all L>=50.",
        "route_guard": "Therefore B_H>A_H and B_T>A_T throughout the q=1 chart. A proof that first applies the triangle inequality separately to carrier atoms, and then adds any nonnegative variation budget, cannot close either anchor disk. This does not lower-bound |mathscr M_N| or disprove signed cancellation.",
        "carrier_coefficient": str(carrier_coefficient),
        "exp_taylor_terms": 41,
        "budget_lower_50": float(budget_lower_50),
        "anchor_upper_50": float(anchor_upper_50),
        "budget_to_anchor_ratio_50": float(budget_lower_50 / anchor_upper_50),
        "absolute_budget_obstructions": 2,
    }


def build_rows(algebra: dict, budget: dict) -> list[GateRow]:
    return [
        GateRow("rsd_01_hermitian", "combined residual", "proved", "The Hermitian disk residual is A_H plus the complete ideal join.", algebra["hermitian_join"], "This is an identity, not a disk estimate."),
        GateRow("rsd_02_transpose", "combined residual", "proved", "The transpose correction belongs inside one combined residual.", algebra["transpose_join"], "R_T and Q_T are not estimated separately."),
        GateRow("rsd_03_improvement", "optimal handoff", "proved", "The combined transpose disk strictly improves the separated triangle handoff.", algebra["transpose_improvement"], "No combined disk bound is claimed."),
        GateRow("rsd_04_cell", "reciprocal cell", "proved", "Each finite reciprocal cell has an exact one-carrier join.", algebra["cell_join"], "Aliases retain symmetric ordering."),
        GateRow("rsd_05_global", "global recomposition", "proved", "Internal ties and aliases recompose into one global exterior.", algebra["global_join"], "The global exterior is not bounded here."),
        GateRow("rsd_06_anchor_stationary", "anchored decomposition", "proved", "Both disk residuals have a complete tie-invariant reciprocal decomposition.", algebra["anchored_stationary"], "No proper fixed-width collar is frozen."),
        GateRow("rsd_07_disk", "disk geometry", "proved", "The tangent anchor disk has an exact signed quadratic equation.", algebra["disk_identity"], "The disk is tangent to zero."),
        GateRow("rsd_08_delta", "quantitative disk", "proved", "The contracted disk gives the exact delta margin.", algebra["contracted_disk"], "This is sufficient, not necessary, for nonvanishing."),
        GateRow("rsd_09_sign", "sign necessity", "proved", "Every nonzero point in the anchor disk has negative raw real part.", algebra["signed_necessity"], "Raw sign is not yet the physical-current sign."),
        GateRow("rsd_10_block", "carrier block", "proved", "A macroscopic interior carrier block has a rational counting floor.", budget["test_block"], "The block is used only for an absolute-budget guard."),
        GateRow("rsd_11_amplitude", "amplitude floor", "proved", "The physical q=1 weight has an explicit lower floor on that block.", budget["amplitude_floor"], "No phase is removed from the true sum."),
        GateRow("rsd_12_hermitian_poly", "Hermitian floor", "proved", "The Hermitian ideal cubic is uniformly nonzero on the test block.", budget["polynomial_floor"], "This is a modulus statement for a route guard."),
        GateRow("rsd_13_transpose_poly", "transpose floor", "proved", "The transpose ideal cubic has the same modulus floor.", budget["polynomial_floor"], "Its real affine part cannot cancel its imaginary cubic part in modulus."),
        GateRow("rsd_14_budget", "absolute carrier budget", "proved", "Both fully termwise carrier budgets grow at least like exp(L/4).", budget["carrier_budget"], "This does not lower-bound the signed carrier sum."),
        GateRow("rsd_15_anchor_upper", "anchor upper bound", "proved", "Both anchors have a polynomial q=1 upper envelope.", budget["anchor_upper"], "The sharper lower anchors remain valid."),
        GateRow("rsd_16_comparison", "uniform mismatch", "proved", "The termwise budgets already exceed the anchors for every L>=50.", budget["uniform_comparison"], "The comparison is a route obstruction only."),
        GateRow("rsd_17_guard", "nonpromotion guard", "proved", "Separated termwise absolute estimates cannot prove either disk inclusion.", budget["route_guard"], "Signed reciprocal cancellation may still close the theorem."),
        GateRow("rsd_18_target", "signed stationary target", "open", "The live target is a signed complete-join estimate in the tangent half-plane, preferably after physical phase composition.", "Bound the tie-invariant M-tau+E package before moduli, or prove a phase-weighted Hermitian/transpose alternative that bypasses the raw anchor disk.", "No stationary sign, physical-current sign, contact exclusion, Lambda<=0, or RH conclusion is claimed."),
    ]


def render_note(payload: dict) -> str:
    algebra = payload["exact"]["algebra"]
    budget = payload["analytic_proof"]["absolute_budget_guard"]
    rows = payload["rows"]
    lines = [
        "# Newman C1 Reciprocal-Stationary Disk Geometry Gate",
        "",
        "Date: 2026-08-04",
        "",
        "Status: exact anchored reciprocal decomposition and disk geometry proved; separated termwise absolute budgets obstructed; signed complete-join estimate open; not a proof of RH.",
        "",
        "## Exact combined residuals",
        "",
        "```text",
        "W_H:=R_H=A_H+mathcal J_H^0,",
        "W_T:=R_T+Q_T=A_T+mathcal J_T^0.",
        "```",
        "",
        algebra["transpose_improvement"],
        "",
        "For each reciprocal cell,",
        "",
        "```text",
        "2M_q-S_q=M_q-tau_q+E_q^(ext)-Lambda_q^(alias).",
        "```",
        "",
        algebra["anchored_stationary"],
        "",
        "## Tangent-disk geometry",
        "",
        "```text",
        "|A+J|^2-A^2=|J|^2+2A Re(J).",
        "",
        "|A+J|<=A-delta",
        "iff |J|^2+2A Re(J)<=-2A delta+delta^2.",
        "```",
        "",
        algebra["signed_necessity"],
        "",
        "## Absolute-budget obstruction",
        "",
        budget["test_block"],
        "",
        budget["amplitude_floor"],
        "",
        budget["polynomial_floor"],
        "",
        "```text",
        "B_H>A_H,   B_T>A_T,                  L>=50, q=1.",
        "```",
        "",
        budget["uniform_comparison"],
        "",
        budget["route_guard"],
        "",
        "## Gate rows",
        "",
        "| id | role | state | claim |",
        "|---|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['id']} | {row['role']} | {row['readiness']} | {row['claim']} |"
        )
    lines.extend(
        [
            "",
            "## Next action",
            "",
            payload["next_action"],
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
    sources = load_sources()
    source_audit = audit_sources(sources)
    algebra = algebra_certificate()
    budget = budget_certificate()
    rows = build_rows(algebra, budget)
    require(len(rows) == 18, "row count")
    require(sum(row.readiness == "open" for row in rows) == 1, "open row count")

    success = (
        "built Newman C1 reciprocal-stationary disk geometry gate: "
        "18 rows, 0 issues, 2 combined residuals, 2 optimal disk identities, "
        "2 absolute-budget obstructions, 1 live signed stationary target"
    )
    payload = {
        "kind": KIND,
        "schema_version": 1,
        "date": "2026-08-04",
        "status": "exact anchored reciprocal decomposition and tangent-disk geometry proved; separated termwise absolute budgets fail; signed complete-join estimate open",
        "source_sha256": {key: item["sha256"] for key, item in source_audit.items()},
        "source_audit": source_audit,
        "exact": {"algebra": algebra},
        "analytic_proof": {"absolute_budget_guard": budget},
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "combined_residuals": algebra["combined_residuals"],
            "optimal_disk_identities": algebra["optimal_disk_identities"],
            "absolute_budget_obstructions": budget["absolute_budget_obstructions"],
            "live_signed_stationary_targets": 1,
        },
        "next_action": "Use the exact global M-tau+E decomposition to seek a signed estimate for the complete Hermitian or transpose join, with aliases, exterior, tie halves, and terminal phase composed before moduli. The raw disk route requires a negative-real tangent-half-plane estimate; alternatively derive a phase-weighted two-channel inequality that proves the physical current directly. Do not return to separated carrier/variation triangle budgets.",
        "proof_boundary": "This gate proves exact combined anchor residuals, an improved transpose handoff, tangent-disk geometry, and a q=1 source-specific obstruction to separated termwise absolute budgets. It proves no signed reciprocal join, disk inclusion, physical-current sign, all-q transport, contact-box exclusion, Q209 shell, cofinal successor, Lambda<=0, RH, or prize-level conclusion.",
        "success": success,
    }

    atomic_write(RESULT_PATH, json.dumps(payload, indent=2) + "\n")
    atomic_write(NOTE_PATH, render_note(payload))
    print(success)


if __name__ == "__main__":
    main()
