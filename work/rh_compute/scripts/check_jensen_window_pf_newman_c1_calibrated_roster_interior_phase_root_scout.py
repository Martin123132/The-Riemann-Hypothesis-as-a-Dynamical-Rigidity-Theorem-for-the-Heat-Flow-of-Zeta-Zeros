#!/usr/bin/env python3
"""Independently check the calibrated roster-interior phase-root scout."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_calibrated_roster_interior_phase_root_scout"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"
CHECK_DPS = 145


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def alpha(s: mp.mpc) -> mp.mpc:
    return 1 / (2 * s) + 1 / (s - 1) + mp.log(s / (2 * mp.pi)) / 2


def alpha_prime(s: mp.mpc) -> mp.mpc:
    return -1 / (2 * s**2) - 1 / (s - 1) ** 2 + 1 / (2 * s)


def independent_state(ell: mp.mpf, expected_n: int) -> dict[str, mp.mpf | mp.mpc | int]:
    time = 1 / (2 * ell**2)
    x = 4 * mp.pi * mp.exp(ell)
    s = (1 - mp.j * x) / 2
    alpha_value = alpha(s)
    alpha_one = alpha_prime(s)
    s_star = s + time * alpha_value / 2
    omega_from_alpha = -mp.im(s_star)
    omega = (
        x / 2
        + mp.atan(x) / (8 * ell**2)
        - 3 * x / (4 * ell**2 * (1 + x**2))
    )
    a_value = mp.sqrt(mp.exp(ell) + 1 / (32 * ell**2))
    n_value = int(mp.floor(a_value))
    require(n_value == expected_n, "independent fixed-N check")

    d_1 = 1 / (6 * s) + alpha_one * (
        time / 4 + time**2 * alpha_value**2 / 8
    )
    log_m_zero = (
        mp.log(s * (s - 1) / 16)
        + mp.log(2 * mp.pi) / 2
        - s * mp.log(mp.pi) / 2
        + (s / 2 - mp.mpf("0.5")) * mp.log(s / 2)
        - s / 2
    )
    log_m_time = time * alpha_value**2 / 4 + log_m_zero
    phi = mp.exp(mp.j * mp.im(log_m_time))
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
    gamma = mp.exp(-2 * mp.j * mp.im(log_m_time))
    psi = mp.arg(q_phase / gamma)
    remainder = psi / 2 + mp.arg(1 + d_1)
    squared_lock_delta = abs(
        eta**2 - mp.exp(2 * mp.j * (theta_0 + remainder))
    )
    return {
        "N": n_value,
        "a": a_value,
        "Omega": omega,
        "omega_delta": abs(omega - omega_from_alpha),
        "B": bracket,
        "B_chart_delta": abs(bracket - bracket_trig),
        "d_1": d_1,
        "psi": psi,
        "x": x,
        "chi_0": chi_0,
        "theta_0": theta_0,
        "squared_lock_delta": squared_lock_delta,
    }


def check_sources(payload: dict) -> int:
    require(set(payload["source_sha256"]) == set(payload["source_audit"]), "source keys")
    require(len(payload["source_audit"]) == 6, "source count")
    for key, item in payload["source_audit"].items():
        path = REPO_ROOT / item["path"]
        require(path.is_file(), f"missing source {key}")
        digest = file_hash(path)
        require(digest == item["sha256"], f"stored source hash {key}")
        require(digest == payload["source_sha256"][key], f"source hash map {key}")
    return 6


def check_ladder(payload: dict) -> int:
    ladder = payload["precision_ladder"]
    require([row["dps"] for row in ladder] == [90, 130, 170], "precision ladder")
    require(len({row["N"] for row in ladder}) == 1, "ladder N")
    require(ladder[-1]["N"] == 72_004_899_338, "physical N")
    for row in ladder:
        level = mp.mpf(-49) / 100
        require(mp.mpf(row["bracket"]["B_early"]) > level, "early sign")
        require(mp.mpf(row["bracket"]["B_negative"]) < level, "negative sign")
        require(mp.mpf(row["root_bracket"]["B_minus_level_lower"]) >= 0, "root lower sign")
        require(mp.mpf(row["root_bracket"]["B_minus_level_upper"]) <= 0, "root upper sign")
        require(mp.mpf(row["roster_margins"]["a_minus_N"]) > mp.mpf("0.8"), "lower roster margin")
        require(mp.mpf(row["roster_margins"]["N_plus_1_minus_a"]) > mp.mpf("0.05"), "upper roster margin")

    comparisons = payload["convergence"]["comparisons"]
    require(len(comparisons) == 2, "convergence comparisons")
    require(mp.mpf(comparisons[0]["L_delta"]) < mp.mpf("1e-50"), "90/130 convergence")
    require(mp.mpf(comparisons[1]["L_delta"]) < mp.mpf("1e-90"), "130/170 convergence")
    require(payload["convergence"]["stable_N"], "stored N stability")
    return len(ladder)


def check_high_row(payload: dict) -> int:
    mp.mp.dps = CHECK_DPS
    high = payload["precision_ladder"][-1]
    n_value = high["N"]
    ell = mp.mpf(high["root"]["L"])
    state = independent_state(ell, n_value)
    level = mp.mpf(-49) / 100
    require(abs(state["B"] - level) < mp.mpf("1e-90"), "independent root residual")
    require(state["omega_delta"] < mp.mpf("1e-115"), "independent Omega identity")
    require(state["B_chart_delta"] < mp.mpf("1e-115"), "independent B chart")
    require(state["squared_lock_delta"] < mp.mpf("1e-110"), "independent phase lock")
    require(abs(state["psi"]) < 1 / state["x"], "independent psi bound")
    require(abs(state["d_1"]) < mp.mpf("1e-22"), "independent d1 scale")
    require(state["a"] - n_value > mp.mpf("0.8"), "independent lower margin")
    require(n_value + 1 - state["a"] > mp.mpf("0.05"), "independent upper margin")

    early = independent_state(mp.mpf(high["bracket"]["L_early"]), n_value)
    negative = independent_state(mp.mpf(high["bracket"]["L_negative"]), n_value)
    require(early["B"] > mp.mpf("-0.01"), "independent early bracket")
    require(negative["B"] < mp.mpf("-0.99"), "independent negative bracket")
    require(
        abs(early["theta_0"] - mp.mpf(high["bracket"]["theta_early_target"]))
        < mp.mpf("1e-90"),
        "early phase target",
    )
    require(
        abs(negative["theta_0"] - mp.mpf(high["bracket"]["theta_negative_target"]))
        < mp.mpf("1e-90"),
        "negative phase target",
    )
    additional = payload["additional_roots"]
    require(
        [row["N"] for row in additional] == [n_value + 1, n_value + 2],
        "adjacent-cell sequence",
    )
    require({row["N"] % 2 for row in [high, *additional]} == {0, 1}, "parity coverage")
    for row in additional:
        row_state = independent_state(mp.mpf(row["root"]["L"]), row["N"])
        require(abs(row_state["B"] - level) < mp.mpf("1e-90"), "adjacent root residual")
        row_early = independent_state(mp.mpf(row["bracket"]["L_early"]), row["N"])
        row_negative = independent_state(mp.mpf(row["bracket"]["L_negative"]), row["N"])
        require(row_early["B"] > mp.mpf("-0.01"), "adjacent early bracket")
        require(row_negative["B"] < mp.mpf("-0.99"), "adjacent negative bracket")
        require(row_state["a"] - row["N"] > mp.mpf("0.05"), "adjacent lower margin")
        require(row["N"] + 1 - row_state["a"] > mp.mpf("0.05"), "adjacent upper margin")
        require(row_state["B_chart_delta"] < mp.mpf("1e-110"), "adjacent B chart")
        require(row_state["squared_lock_delta"] < mp.mpf("1e-105"), "adjacent lock")
        require(abs(row_state["psi"]) < 1 / row_state["x"], "adjacent psi bound")
    return 16


def check_rows(payload: dict) -> None:
    rows = payload["rows"]
    require(len(rows) == 12, "row count")
    require(len({row["id"] for row in rows}) == 12, "row ids")
    require(sum(row["readiness"] == "open" for row in rows) == 1, "open rows")
    require(sum(row["readiness"] == "diagnostic_validated" for row in rows) == 9, "diagnostic rows")
    require(rows[-1]["readiness"] == "guard_validated", "boundary row")
    require(all(row["claim"] and row["certificate"] and row["proof_boundary"] for row in rows), "row fields")


def main() -> int:
    require(RESULT_PATH.is_file(), "missing result")
    require(NOTE_PATH.is_file(), "missing note")
    payload = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    require(payload["kind"] == KIND, "kind")
    require(payload["schema_version"] == 1, "schema")
    require(payload["level"] == {"numerator": -49, "denominator": 100}, "level")
    source_audits = check_sources(payload)
    precision_ladders = check_ladder(payload)
    numerical_audits = check_high_row(payload)
    check_rows(payload)

    summary = payload["summary"]
    require(summary["rows"] == 12 and summary["issues"] == 0, "summary rows")
    require(summary["source_audits"] == 6, "summary sources")
    require(summary["precision_ladders"] == 3, "summary ladder")
    require(summary["physical_q1_cells"] == 3, "summary cell")
    require(summary["calibrated_roots"] == 3, "summary root")
    require(summary["fixed_roster_margins"] == 3, "summary margin")
    require(summary["evaluated_determinant_rows"] == 0, "summary determinant boundary")
    require(summary["open_rows"] == 1, "summary open")

    note = NOTE_PATH.read_text(encoding="utf-8")
    require(payload["success"] in note, "note success")
    require("72004899338" in note and "B = -0.49" in note, "root note")
    require("not an interval enclosure" in note, "interval boundary")
    require("0 determinant rows" in note, "determinant boundary")
    require("not a proof of RH" in note, "RH boundary")
    require("prize-level conclusion" in note, "prize boundary")

    success = (
        "validated Newman C1 calibrated roster-interior phase-root scout: "
        f"12 rows, 0 issues, {source_audits} source audits, "
        f"{precision_ladders} precision ladders, {numerical_audits} independent numerical audits, "
        "3 physical q=1 cells, 3 calibrated B=-49/100 roots, "
        "3 fixed-roster margins, 0 determinant rows, 1 open evaluator row"
    )
    print(success)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
