#!/usr/bin/env python3
"""Build a high-precision calibrated roster-interior phase-root scout."""

from __future__ import annotations

import ctypes
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_calibrated_roster_interior_phase_root_scout"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "raw_center": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_raw_center_sign_obstruction_gate.json",
    "phase_lock": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_saddle_phase_lock_center_symmetry_gate.json",
    "normalizer": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_reciprocal_normalizer_phase_reinforcement_guard.json",
    "absolute_anchor": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_first_order_centered_absolute_phase_anchor_reduction.json",
    "q1_phase": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_five_moment_q1_saddle_phase_variance_reduction.json"
    ),
    "determinant_chart": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_finite_determinant_reciprocal_block_symmetry_gate.json",
}

L_MIN = mp.mpf(50)
LEVEL = Fraction(-49, 100)
PRECISION_LADDER = (90, 130, 170)


@dataclass(frozen=True)
class ScoutRow:
    id: str
    role: str
    readiness: str
    claim: str
    certificate: str
    proof_boundary: str


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def set_below_normal_priority() -> None:
    if os.name != "nt":
        return
    try:
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.GetCurrentProcess.restype = ctypes.c_void_p
        kernel32.SetPriorityClass.argtypes = (ctypes.c_void_p, ctypes.c_uint32)
        kernel32.SetPriorityClass.restype = ctypes.c_int
        kernel32.SetPriorityClass(kernel32.GetCurrentProcess(), 0x00004000)
    except (AttributeError, OSError):
        pass


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
        payloads["raw_center"]["summary"]["strict_negative_witnesses"] == 1,
        "raw-center witness drifted",
    )
    require(
        payloads["phase_lock"]["summary"]["exact_center_factorizations"] == 2,
        "phase-lock factorization drifted",
    )
    require(
        "M_0(s)=sqrt(2*pi)"
        in payloads["normalizer"]["exact"]["source_coordinate"],
        "normalizer formula drifted",
    )
    require(
        "eta=f_1/|f_1|"
        in payloads["absolute_anchor"]["exact"]["unit_anchor"],
        "absolute anchor drifted",
    )
    require(
        "Omega=-Im(s_*)"
        in payloads["q1_phase"]["physical_q1_certificate"]["frequency"],
        "q1 frequency drifted",
    )
    require(
        payloads["determinant_chart"]["summary"]["reciprocal_block_classes"] == 3,
        "determinant block chart drifted",
    )
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def alpha(s: mp.mpc) -> mp.mpc:
    return 1 / (2 * s) + 1 / (s - 1) + mp.log(s / (2 * mp.pi)) / 2


def alpha_prime(s: mp.mpc) -> mp.mpc:
    return -1 / (2 * s**2) - 1 / (s - 1) ** 2 + 1 / (2 * s)


def m_zero(s: mp.mpc) -> mp.mpc:
    return (
        s
        * (s - 1)
        / 16
        * mp.pi ** (-s / 2)
        * mp.sqrt(2 * mp.pi)
        * mp.exp((s / 2 - mp.mpf("0.5")) * mp.log(s / 2) - s / 2)
    )


def a_squared(ell: mp.mpf) -> mp.mpf:
    return mp.exp(ell) + 1 / (32 * ell**2)


def omega_explicit(ell: mp.mpf) -> mp.mpf:
    x = 4 * mp.pi * mp.exp(ell)
    return (
        x / 2
        + mp.atan(x) / (8 * ell**2)
        - 3 * x / (4 * ell**2 * (1 + x**2))
    )


def bisect_monotone(
    function,
    target: mp.mpf,
    lower: mp.mpf,
    upper: mp.mpf,
    *,
    increasing: bool,
    tolerance: mp.mpf,
) -> tuple[mp.mpf, mp.mpf]:
    lower_value = function(lower) - target
    upper_value = function(upper) - target
    if increasing:
        require(lower_value <= 0 <= upper_value, "increasing bisection bracket")
    else:
        require(lower_value >= 0 >= upper_value, "decreasing bisection bracket")
    iterations = 0
    while upper - lower > tolerance:
        middle = (lower + upper) / 2
        middle_value = function(middle) - target
        if (increasing and middle_value < 0) or (
            not increasing and middle_value > 0
        ):
            lower = middle
        else:
            upper = middle
        iterations += 1
        require(iterations < 2_000, "bisection iteration cap")
    return lower, upper


def cell_endpoint(integer_a: int, tolerance: mp.mpf) -> mp.mpf:
    guess = 2 * mp.log(integer_a)
    lower, upper = bisect_monotone(
        a_squared,
        mp.mpf(integer_a) ** 2,
        guess - mp.mpf("1e-8"),
        guess + mp.mpf("1e-8"),
        increasing=True,
        tolerance=tolerance,
    )
    return (lower + upper) / 2


def phase_geometry(ell: mp.mpf, n_value: int) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    omega = omega_explicit(ell)
    k_n = 2 * mp.pi * n_value**2
    chi_0 = omega * (1 + mp.log(k_n / omega)) / 2 + mp.pi / 8
    theta_0 = chi_0 - omega * mp.log(n_value)
    return omega, chi_0, theta_0


def physical_state(ell: mp.mpf, expected_n: int) -> dict[str, mp.mpf | mp.mpc | int]:
    time = 1 / (2 * ell**2)
    x = 4 * mp.pi * mp.exp(ell)
    s = (1 - mp.j * x) / 2
    alpha_value = alpha(s)
    alpha_one = alpha_prime(s)
    s_star = s + time * alpha_value / 2
    omega_alpha = -mp.im(s_star)
    omega = omega_explicit(ell)
    a_value = mp.sqrt(a_squared(ell))
    n_value = int(mp.floor(a_value))
    require(n_value == expected_n, "root left fixed-N cell")

    d_1 = 1 / (6 * s) + alpha_one * (
        time / 4 + time**2 * alpha_value**2 / 8
    )
    normalizer = mp.exp(time * alpha_value**2 / 4) * m_zero(s)
    phi = normalizer / abs(normalizer)
    eta = phi * (1 + d_1) / abs(1 + d_1)
    phi_n = omega * mp.log(n_value)
    terminal_phase = mp.exp(mp.j * phi_n)
    u_n = mp.log(a_value / n_value)
    mismatch = u_n / mp.log(n_value)
    bracket = (
        mp.re(eta) * mp.im(eta * terminal_phase)
        + mismatch * mp.im(terminal_phase)
    )
    theta_eta = mp.arg(eta)
    bracket_trig = (
        mp.cos(theta_eta) * mp.sin(phi_n + theta_eta)
        + mismatch * mp.sin(phi_n)
    )

    k_n = 2 * mp.pi * n_value**2
    chi_0 = omega * (1 + mp.log(k_n / omega)) / 2 + mp.pi / 8
    theta_0 = chi_0 - phi_n
    q_phase = mp.exp(
        mp.j * (omega * mp.log(omega / (2 * mp.pi)) - omega - mp.pi / 4)
    )
    gamma = mp.conj(normalizer) / normalizer
    psi = mp.arg(q_phase / gamma)
    delta_1 = mp.arg(1 + d_1)
    remainder = psi / 2 + delta_1
    squared_lock_delta = abs(
        eta**2 - mp.exp(2 * mp.j * (theta_0 + remainder))
    )
    return {
        "L": ell,
        "t": time,
        "x": x,
        "s": s,
        "s_star": s_star,
        "a": a_value,
        "N": n_value,
        "Omega": omega,
        "omega_formula_delta": abs(omega - omega_alpha),
        "d_1": d_1,
        "eta": eta,
        "phi_N": phi_n,
        "u_N": u_n,
        "mismatch": mismatch,
        "B": bracket,
        "B_trig": bracket_trig,
        "B_chart_delta": abs(bracket - bracket_trig),
        "chi_0": chi_0,
        "theta_0": theta_0,
        "psi": psi,
        "delta_1": delta_1,
        "phase_remainder": remainder,
        "squared_lock_delta": squared_lock_delta,
    }


def text(value: mp.mpf | mp.mpc, digits: int) -> str:
    return mp.nstr(value, max(30, digits), strip_zeros=False)


def state_payload(state: dict[str, mp.mpf | mp.mpc | int], digits: int) -> dict:
    result: dict[str, str | int] = {}
    for key, value in state.items():
        if isinstance(value, int):
            result[key] = value
        else:
            result[key] = text(value, digits)
    return result


def compute_witness(dps: int, cell_offset: int = 0) -> dict:
    mp.mp.dps = dps
    tolerance = mp.power(10, -(dps - 32))
    level = mp.mpf(LEVEL.numerator) / LEVEL.denominator

    a_at_floor = mp.sqrt(a_squared(mp.mpf(50)))
    n_value = int(mp.floor(a_at_floor)) + 1 + cell_offset
    left = cell_endpoint(n_value, tolerance)
    right = cell_endpoint(n_value + 1, tolerance)
    k_n = 2 * mp.pi * n_value**2
    cross_lo, cross_hi = bisect_monotone(
        omega_explicit,
        k_n,
        left,
        right,
        increasing=True,
        tolerance=tolerance,
    )
    crossing = (cross_lo + cross_hi) / 2

    chi_cross = mp.pi * n_value**2 + mp.pi / 8
    chi_turn = int(mp.floor((chi_cross - mp.pi / 2) / (2 * mp.pi)))
    chi_target = mp.pi / 2 + 2 * mp.pi * chi_turn
    chi_function = lambda ell: phase_geometry(ell, n_value)[1]
    chi_lo, chi_hi = bisect_monotone(
        chi_function,
        chi_target,
        crossing,
        right,
        increasing=False,
        tolerance=tolerance,
    )
    chi_point = (chi_lo + chi_hi) / 2

    theta_function = lambda ell: phase_geometry(ell, n_value)[2]
    theta_at_chi = theta_function(chi_point)
    theta_turn = int(mp.nint((theta_at_chi - mp.pi) / (2 * mp.pi)))
    negative_target = mp.pi + 2 * mp.pi * theta_turn
    neg_lo, neg_hi = bisect_monotone(
        theta_function,
        negative_target,
        crossing,
        right,
        increasing=False,
        tolerance=tolerance,
    )
    negative_point = (neg_lo + neg_hi) / 2

    early_target = 3 * mp.pi / 2 + 2 * mp.pi * theta_turn
    early_lo, early_hi = bisect_monotone(
        theta_function,
        early_target,
        crossing,
        negative_point,
        increasing=False,
        tolerance=tolerance,
    )
    early_point = (early_lo + early_hi) / 2
    early_state = physical_state(early_point, n_value)
    negative_state = physical_state(negative_point, n_value)
    require(early_state["B"] > level, "early endpoint did not exceed level")
    require(negative_state["B"] < level, "negative endpoint did not lie below level")

    root_lower, root_upper = early_point, negative_point
    root_lower_value = early_state["B"] - level
    root_upper_value = negative_state["B"] - level
    root_iterations = 0
    while root_upper - root_lower > tolerance:
        middle = (root_lower + root_upper) / 2
        middle_value = physical_state(middle, n_value)["B"] - level
        if middle_value > 0:
            root_lower = middle
            root_lower_value = middle_value
        else:
            root_upper = middle
            root_upper_value = middle_value
        root_iterations += 1
        require(root_iterations < 2_000, "level-root bisection iteration cap")
    root = (root_lower + root_upper) / 2
    root_state = physical_state(root, n_value)
    root_value = root_state["B"] - level

    require(root_lower_value >= 0 >= root_upper_value, "final level bracket")
    require(abs(root_state["omega_formula_delta"]) < mp.power(10, -(dps - 35)), "Omega formula")
    require(abs(root_state["B_chart_delta"]) < mp.power(10, -(dps - 35)), "B chart")
    require(abs(root_state["squared_lock_delta"]) < mp.power(10, -(dps - 35)), "phase lock")
    require(abs(root_state["psi"]) < 1 / root_state["x"], "stationary phase bound")
    require(root_state["a"] - n_value > mp.mpf("0.1"), "lower roster margin")
    require(n_value + 1 - root_state["a"] > mp.mpf("0.01"), "upper roster margin")

    digits = dps - 12
    return {
        "dps": dps,
        "N": n_value,
        "chi_turn": chi_turn,
        "theta_turn": theta_turn,
        "cell": {
            "L_left": text(left, digits),
            "L_right": text(right, digits),
            "L_width": text(right - left, digits),
            "Omega_crossing_L": text(crossing, digits),
            "safe_chi_L": text(chi_point, digits),
            "safe_chi_target": text(chi_target, digits),
        },
        "bracket": {
            "L_early": text(early_point, digits),
            "B_early": text(early_state["B"], digits),
            "L_negative": text(negative_point, digits),
            "B_negative": text(negative_state["B"], digits),
            "theta_early_target": text(early_target, digits),
            "theta_negative_target": text(negative_target, digits),
        },
        "root_bracket": {
            "L_lower": text(root_lower, digits),
            "B_minus_level_lower": text(root_lower_value, digits),
            "L_upper": text(root_upper, digits),
            "B_minus_level_upper": text(root_upper_value, digits),
            "width": text(root_upper - root_lower, digits),
            "iterations": root_iterations,
        },
        "root": state_payload(root_state, digits),
        "root_level_residual": text(root_value, digits),
        "roster_margins": {
            "a_minus_N": text(root_state["a"] - n_value, digits),
            "N_plus_1_minus_a": text(n_value + 1 - root_state["a"], digits),
            "L_minus_left": text(root - left, digits),
            "right_minus_L": text(right - root, digits),
        },
    }


def build_convergence(ladder: list[dict]) -> dict:
    mp.mp.dps = max(PRECISION_LADDER) + 20
    comparisons = []
    for lower, upper in zip(ladder, ladder[1:]):
        root_delta = abs(mp.mpf(lower["root"]["L"]) - mp.mpf(upper["root"]["L"]))
        a_delta = abs(mp.mpf(lower["root"]["a"]) - mp.mpf(upper["root"]["a"]))
        require(root_delta < mp.power(10, -(lower["dps"] - 38)), "root convergence")
        # da/dL is approximately a/2 here, so absolute a errors lose eleven
        # decimal places relative to the L root without losing roster digits.
        require(a_delta < mp.power(10, -(lower["dps"] - 50)), "a convergence")
        comparisons.append(
            {
                "dps_pair": [lower["dps"], upper["dps"]],
                "L_delta": text(root_delta, 50),
                "a_delta": text(a_delta, 50),
            }
        )
    return {
        "comparisons": comparisons,
        "stable_N": len({row["N"] for row in ladder}) == 1,
        "precision_ladders": len(ladder),
    }


def exact_certificate(
    high: dict,
    additional_roots: list[dict],
) -> dict[str, str | int | list[str]]:
    additional_summary = "; ".join(
        f"N={row['N']}, L={row['root']['L']}, B={row['root']['B']}"
        for row in additional_roots
    )
    return {
        "physical_chart": "On q=1 use t=1/(2L^2), x=4pi exp(L), s=(1-ix)/2, s_*=s+t alpha(s)/2, Omega=-Im(s_*), and a^2=exp(L)+1/(32L^2). The first complete cell above L=50 has N=floor(a)=72004899338.",
        "source_phase": "Use the exact Polymath-15 M_t(s)=exp[t alpha(s)^2/4]M_0(s), d_1=1/(6s)+alpha'(s){t/4+t^2alpha(s)^2/8}, and eta=(M_t/|M_t|)(1+d_1)/|1+d_1|.",
        "branch_free_bracket": "With phi_N=Omega log N and u_N=log(a/N), evaluate B=Re(eta)Im(eta exp(i phi_N))+(u_N/log N)Im(exp(i phi_N)). This equals the trigonometric centre factor without choosing an argument branch.",
        "root_construction": "After Omega crosses K_N=2pi N^2, solve chi_0=F_N(Omega)/2+pi/8 at the congruence chi_0=pi/2 mod 2pi. In that safe arc solve theta_0=chi_0-phi_N first at 3pi/2 and then at pi modulo 2pi. The resulting endpoint values straddle B=-49/100, so bisection constructs a roster-interior numerical level point.",
        "root_summary": f"At {high['dps']} decimal digits the root lies at L={high['root']['L']}, with N={high['N']}, B={high['root']['B']}, a-N={high['roster_margins']['a_minus_N']}, and N+1-a={high['roster_margins']['N_plus_1_minus_a']}.",
        "additional_root_summary": "The next two complete cells give " + additional_summary + ".",
        "available_fields": [
            "L,t,x,s,s_*,Omega,a,N",
            "M_t(s),d_1,eta,phi_N,u_N",
            "B and its safe-arc phase coordinates",
            "fixed-cell margins and precision-ladder root brackets",
        ],
        "missing_determinant_fields": [
            "omega_E and dot(omega_E)",
            "retained observation p and dot(p)",
            "finite block observations f^(0),f^(1),f^(2) and their tangents",
            "independent near observation n and common-cutoff paired-remote observation c",
            "det(B_mat)=h^2P_ret and every determinant current increment",
        ],
        "scale_guard": "The calibrated cell has N=72004899338. Direct enumeration of N carriers or an N by N reciprocal block table is not an admissible next algorithm. The missing evaluator must use the certified stationary diagonal, adjacent recurrence, and far paired/compressed representations.",
        "next_action": "Use this saved root as the immutable input row for a determinant-evaluator contract. First implement omega_E, p, and the diagonal-plus-adjacent observation through the exact physical coefficient and Morse/adjacent formulas; cross-check it against the near-first chart. Add the far sector only through a resumable paired or Type-I/II compression. Do not infer a sign theorem from this phase scout.",
        "pi_provenance": "Every pi here is inherited from the completed-zeta normalizer, T_0=2pi a^2, Q(v)=exp(i[v log(v/(2pi))-v-pi/4]), and the period of the complex exponential. No fitted circle or polygon supplies pi.",
        "proof_boundary": "This is reproducible high-precision diagnostic evidence for calibrated physical phase levels in three consecutive fixed q=1 roster interiors. It is not interval arithmetic and proves no uniqueness of any B root, no all-cell statement, determinant-row value, P_ret nonsingularity, geometric-current sign or bound, complete-current sign, all-q transport, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
        "calibrated_roots": 1 + len(additional_roots),
        "physical_phase_rows": 1 + len(additional_roots),
        "evaluated_determinant_rows": 0,
    }


def build_rows(exact: dict, convergence: dict) -> list[ScoutRow]:
    return [
        ScoutRow("cpr_01_chart", "physical chart", "diagnostic_validated", "The q=1 chart is evaluated from its exact L-dependent formulas.", str(exact["physical_chart"]), "No asymptotic replacement of Omega or a is used."),
        ScoutRow("cpr_02_source", "source phase", "diagnostic_validated", "The actual normalizer and first correction determine eta.", str(exact["source_phase"]), "The source phase is not fitted or freely rotated."),
        ScoutRow("cpr_03_bracket", "branch-free bracket", "diagnostic_validated", "B is evaluated without an argument branch.", str(exact["branch_free_bracket"]), "The equivalent trigonometric chart is checked numerically."),
        ScoutRow("cpr_04_safe", "safe phase arc", "diagnostic_validated", "An explicit safe chi_0 point selects the local source half-turn.", str(exact["root_construction"]), "This selects one of many possible roots."),
        ScoutRow("cpr_05_straddle", "level straddle", "diagnostic_validated", "Explicit early and negative phase points straddle -49/100.", str(exact["root_construction"]), "Floating-point signs are diagnostic, not interval certificates."),
        ScoutRow("cpr_06_root", "calibrated roots", "diagnostic_validated", "B=-49/100 is resolved in three consecutive physical roster interiors.", str(exact["root_summary"]) + " " + str(exact["additional_root_summary"]), "No uniqueness in any cell is asserted."),
        ScoutRow("cpr_07_ladder", "precision ladder", "diagnostic_validated", "Three independently recomputed decimal precisions converge to the same row.", json.dumps(convergence, sort_keys=True), "Cross-precision agreement does not turn mpmath into interval arithmetic."),
        ScoutRow("cpr_08_roster", "roster margin", "diagnostic_validated", "The selected point is separated from both fixed-N cell boundaries.", str(exact["root_summary"]), "This is a numerical margin at the selected row."),
        ScoutRow("cpr_09_available", "available inputs", "diagnostic_validated", "The phase-level input contract is now executable.", "; ".join(exact["available_fields"]), "Observation matrices are not included."),
        ScoutRow("cpr_10_missing", "determinant gap", "open", "The physical determinant-row evaluator remains absent.", "; ".join(exact["missing_determinant_fields"]), str(exact["scale_guard"])),
        ScoutRow("cpr_11_scale", "algorithmic scale", "guard_validated", "Direct carrier or block-table enumeration is rejected at the physical N.", str(exact["scale_guard"]), "A compressed arithmetic representation is required."),
        ScoutRow("cpr_12_boundary", "proof boundary", "guard_validated", "A calibrated numerical row is not a current theorem.", str(exact["proof_boundary"]), "All prize-level conclusions remain open."),
    ]


def render_note(payload: dict) -> str:
    high = payload["precision_ladder"][-1]
    exact = payload["exact"]
    lines = [
        "# Newman C1 Calibrated Roster-Interior Phase-Root Scout",
        "",
        "Date: 2026-08-05",
        "",
        "Status: high-precision physical phase diagnostic validated; determinant observations remain open; not a proof of RH.",
        "",
        "## Physical Row",
        "",
        exact["physical_chart"],
        "",
        exact["source_phase"],
        "",
        exact["branch_free_bracket"],
        "",
        "```text",
        f"N = {high['N']}",
        f"L = {high['root']['L']}",
        f"a-N = {high['roster_margins']['a_minus_N']}",
        f"N+1-a = {high['roster_margins']['N_plus_1_minus_a']}",
        f"B = {high['root']['B']}",
        f"B+49/100 = {high['root_level_residual']}",
        "```",
        "",
        "## Construction",
        "",
        exact["root_construction"],
        "",
        f"The early bracket value is `{high['bracket']['B_early']}` and the negative bracket value is `{high['bracket']['B_negative']}`. The final L-bracket width is `{high['root_bracket']['width']}`.",
        "",
        "## Precision Ladder",
        "",
        "The row was rebuilt from scratch at 90, 130, and 170 decimal digits. Consecutive L differences are:",
        "",
    ]
    for comparison in payload["convergence"]["comparisons"]:
        lines.append(
            f"- `{comparison['dps_pair'][0]} -> {comparison['dps_pair'][1]}`: `{comparison['L_delta']}`"
        )
    lines.extend(
        [
            "",
            "## Adjacent Cells",
            "",
            "The same branch-free construction was independently run at 170 digits in the next two complete cells:",
            "",
        ]
    )
    for row in payload["additional_roots"]:
        lines.append(
            f"- `N={row['N']}`: `L={row['root']['L']}`, `B={row['root']['B']}`, `a-N={row['roster_margins']['a_minus_N']}`"
        )
    lines.extend(
        [
            "",
            "This is a convergence diagnostic, not an interval enclosure.",
            "",
            "## Determinant Input Gap",
            "",
            "The following fields are still missing:",
            "",
        ]
    )
    lines.extend(f"- `{item}`" for item in exact["missing_determinant_fields"])
    lines.extend(
        [
            "",
            exact["scale_guard"],
            "",
            "## Next Action",
            "",
            exact["next_action"],
            "",
            "## Pi Provenance",
            "",
            exact["pi_provenance"],
            "",
            "## Proof Boundary",
            "",
            exact["proof_boundary"],
            "",
            payload["success"],
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    set_below_normal_priority()
    source_audit = load_and_audit_sources()
    ladder = [compute_witness(dps) for dps in PRECISION_LADDER]
    additional_roots = [
        compute_witness(max(PRECISION_LADDER), cell_offset)
        for cell_offset in (1, 2)
    ]
    convergence = build_convergence(ladder)
    require(convergence["stable_N"], "N changed across precision ladder")
    require(
        [row["N"] for row in additional_roots]
        == [ladder[-1]["N"] + 1, ladder[-1]["N"] + 2],
        "adjacent-cell sequence",
    )
    exact = exact_certificate(ladder[-1], additional_roots)
    rows = build_rows(exact, convergence)
    require(len(rows) == 12, "row count")
    require(sum(row.readiness == "open" for row in rows) == 1, "open row count")

    success = (
        "built Newman C1 calibrated roster-interior phase-root scout: "
        "12 rows, 0 issues, 6 source audits, 3 precision ladders, "
        "3 physical q=1 cells, 3 calibrated B=-49/100 roots, "
        "3 fixed-roster margins, 0 determinant rows, 1 open evaluator row"
    )
    payload = {
        "kind": KIND,
        "schema_version": 1,
        "date": "2026-08-05",
        "status": "calibrated physical phase root located; determinant-row evaluator open",
        "source_sha256": {key: item["sha256"] for key, item in source_audit.items()},
        "source_audit": source_audit,
        "level": {"numerator": LEVEL.numerator, "denominator": LEVEL.denominator},
        "precision_ladder": ladder,
        "additional_roots": additional_roots,
        "convergence": convergence,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "issues": 0,
            "source_audits": len(source_audit),
            "precision_ladders": len(ladder),
            "physical_q1_cells": 1 + len(additional_roots),
            "calibrated_roots": exact["calibrated_roots"],
            "fixed_roster_margins": 1 + len(additional_roots),
            "evaluated_determinant_rows": exact["evaluated_determinant_rows"],
            "open_rows": sum(row.readiness == "open" for row in rows),
        },
        "next_action": exact["next_action"],
        "pi_provenance": exact["pi_provenance"],
        "proof_boundary": exact["proof_boundary"],
        "success": success,
    }
    atomic_write(RESULT_PATH, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(payload))
    print(success)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
