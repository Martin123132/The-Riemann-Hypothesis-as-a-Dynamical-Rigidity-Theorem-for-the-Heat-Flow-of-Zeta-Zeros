#!/usr/bin/env python3
"""Independently check the physical scalar and genuine-edge adapter."""

from __future__ import annotations

import ctypes
import hashlib
import json
import os
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_physical_scalar_edge_adapter_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"
CHECK_DPS = 155
ROW_ORDER = ("P_V", "P_N", "P_A", "P_Q")


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


def decode(value: dict[str, str]) -> mp.mpc:
    return mp.mpc(mp.mpf(value["re"]), mp.mpf(value["im"]))


def alpha(s: mp.mpc) -> mp.mpc:
    return 1 / (2 * s) + 1 / (s - 1) + mp.log(s / (2 * mp.pi)) / 2


def alpha_prime(s: mp.mpc) -> mp.mpc:
    return -1 / (2 * s**2) - 1 / (s - 1) ** 2 + 1 / (2 * s)


def chart(x: mp.mpf | mp.mpc, time: mp.mpf, n_value: int) -> dict[str, mp.mpf | mp.mpc]:
    s = (1 - mp.j * x) / 2
    a = mp.sqrt(x / (4 * mp.pi) + time / 16)
    h = 1 / a
    theta = a - n_value
    return {
        "s": s,
        "a": a,
        "h": h,
        "theta": theta,
        "p": 1 - 2 * theta,
        "T_0": 2 * mp.pi / h**2,
        "u_N": mp.log(a / n_value),
    }


def scaled_diff(function, x: mp.mpf, order: int = 1) -> mp.mpc:
    return mp.diff(lambda z: function(x * (1 + z)), 0, order, addprec=55) / x**order


def correction_value(x: mp.mpf | mp.mpc, time: mp.mpf, log_n: mp.mpf) -> mp.mpc:
    s = (1 - mp.j * x) / 2
    a = alpha(s)
    a_one = alpha_prime(s)
    return 1 / (6 * s) + a_one * (
        time / 4 + time**2 * (a - log_n) ** 2 / 8
    )


def normalized_correction(
    x: mp.mpf | mp.mpc, time: mp.mpf, log_value: mp.mpf
) -> mp.mpc:
    return (1 + correction_value(x, time, log_value)) / (
        1 + correction_value(x, time, mp.mpf(0))
    )


def c0(p: mp.mpf | mp.mpc) -> mp.mpc:
    return (
        mp.exp(mp.j * mp.pi * (p**2 / 2 + mp.mpf(3) / 8))
        - mp.j * mp.sqrt(2) * mp.cos(mp.pi * p / 2)
    ) / (2 * mp.cos(mp.pi * p))


def endpoint_h(x: mp.mpf | mp.mpc, time: mp.mpf, n_value: int) -> mp.mpc:
    state = chart(x, time, n_value)
    p = state["p"]
    c_1 = mp.diff(c0, p, 3, addprec=55) / (12 * mp.pi**2)
    return c0(p) + state["h"] * c_1


def s_star(x: mp.mpf | mp.mpc, time: mp.mpf) -> mp.mpc:
    s = (1 - mp.j * x) / 2
    return s + time * alpha(s) / 2


def endpoint_mu(x: mp.mpf | mp.mpc, time: mp.mpf, n_value: int) -> mp.mpc:
    state = chart(x, time, n_value)
    s = state["s"]
    a_value = alpha(s)
    a_one = alpha_prime(s)
    chi = a_value - mp.log(state["a"])
    t_0 = state["T_0"]
    k_rate = -mp.pi / 8 + 1 / (4 * t_0) + 1 / (2 * (t_0 + mp.j))
    return k_rate + mp.j * chi / 2 + mp.j * (time / 4) * a_one * chi


def endpoint_j(x: mp.mpf, time: mp.mpf, n_value: int) -> mp.mpc:
    h_value = endpoint_h(x, time, n_value)
    h_x = scaled_diff(lambda xx: endpoint_h(xx, time, n_value), x)
    return h_x + endpoint_mu(x, time, n_value) * h_value


def edge_value(x: mp.mpf, time: mp.mpf, n_value: int) -> mp.mpf:
    state = chart(x, time, n_value)
    tau = 1 + mp.j / state["T_0"]
    return mp.re(tau * endpoint_h(x, time, n_value))


def edge_scalar(x: mp.mpf, time: mp.mpf, n_value: int) -> mp.mpf:
    state = chart(x, time, n_value)
    tau = 1 + mp.j / state["T_0"]
    s_prime = scaled_diff(lambda xx: s_star(xx, time), x)
    alpha_c = mp.re(s_prime) * state["u_N"]
    return mp.re(tau * endpoint_j(x, time, n_value)) - alpha_c * edge_value(
        x, time, n_value
    )


def terminal_chi(x: mp.mpf, time: mp.mpf, n_value: int) -> mp.mpc:
    state = chart(x, time, n_value)
    s = state["s"]
    a_value = alpha(s)
    a_one = alpha_prime(s)
    h = state["h"]
    theta = state["theta"]
    t_0 = state["T_0"]
    log_n = mp.log(n_value)
    d_1 = correction_value(x, time, mp.mpf(0))
    d_n = correction_value(x, time, log_n)
    d_1_x = scaled_diff(lambda xx: correction_value(xx, time, mp.mpf(0)), x)
    d_n_x = scaled_diff(lambda xx: correction_value(xx, time, log_n), x)
    stable_ell = mp.log(n_value / state["a"])
    saddle_rho = time * h**2 / 16
    saddle_r = 1 - saddle_rho
    y = h**2 / (4 * mp.pi * saddle_r)
    chi_r = mp.log(saddle_r) / 2 + mp.log(1 + y**2) / 4 - y**2 / (1 + y**2)
    eta_small = mp.atan(y) / 2 + 3 * y / (1 + y**2)
    w = a_value - log_n
    w_x = (
        eta_small / 2
        - 3 / (4 * t_0)
        + mp.j * (stable_ell + h * theta - chi_r) / 2
        - mp.j * (time / 4) * a_one * w
        + d_n_x / (1 + d_n)
    )
    delta_n = d_n_x / (1 + d_n) - d_1_x / (1 + d_1)
    return -mp.j * h * theta / 2 + w_x - delta_n


def polynomial_add(left: list[mp.mpc], right: list[mp.mpc]) -> list[mp.mpc]:
    size = max(len(left), len(right))
    return [
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(size)
    ]


def polynomial_scale(scale: mp.mpc, values: list[mp.mpc]) -> list[mp.mpc]:
    return [scale * value for value in values]


def polynomial_multiply(left: list[mp.mpc], right: list[mp.mpc]) -> list[mp.mpc]:
    result = [mp.mpc(0) for _ in range(len(left) + len(right) - 1)]
    for i, a_value in enumerate(left):
        for j, b_value in enumerate(right):
            result[i + j] += a_value * b_value
    return result


def pad(values: list[mp.mpc], size: int) -> list[mp.mpc]:
    require(len(values) <= size, "polynomial overflow")
    return [*values, *([mp.mpc(0)] * (size - len(values)))]


def lift(values: list[mp.mpc], log_a: mp.mpf) -> list[mp.mpc]:
    result = [mp.mpc(0) for _ in range(len(values) + 1)]
    for index, value in enumerate(values):
        result[index] -= mp.j * log_a * value
        result[index + 1] += mp.j * value
    return result


def independent_rows(table: dict[str, mp.mpf | mp.mpc]) -> dict[str, list[mp.mpc]]:
    c_poly = [1, table["rho_1"], table["rho_2"]]
    d_poly = [0, table["rho_1_x"], table["rho_2_x"]]
    r_poly = [
        table["s_prime"] * table["log_N"] + mp.j * table["b"] * table["u_N"],
        -table["s_prime"],
    ]
    r_x_poly = [
        table["s_second"] * table["log_N"]
        + mp.j * (table["b_x"] * table["u_N"] + table["b"] * table["u_N_x"]),
        -table["s_second"],
    ]
    q_rate = [table["chi_rate"] + table["s_prime"] * table["log_N"], -table["s_prime"]]
    q_poly = polynomial_add(polynomial_multiply(q_rate, c_poly), d_poly)
    a_poly = polynomial_multiply(r_poly, c_poly)
    n_poly = polynomial_add(
        polynomial_multiply(r_poly, q_poly), polynomial_multiply(r_x_poly, c_poly)
    )
    return {
        "P_V": pad(c_poly, 5),
        "P_N": pad(n_poly, 5),
        "P_A": pad(a_poly, 5),
        "P_Q": pad(q_poly, 5),
    }


def maximum_delta(left: list[mp.mpc], right: list[mp.mpc]) -> mp.mpf:
    require(len(left) == len(right), "vector length")
    return max(abs(a_value - b_value) for a_value, b_value in zip(left, right))


def check_sources(payload: dict) -> int:
    require(set(payload["source_sha256"]) == set(payload["source_audit"]), "source keys")
    for key, item in payload["source_audit"].items():
        path = REPO_ROOT / item["path"]
        require(path.is_file(), f"missing source {key}")
        digest = file_hash(path)
        require(digest == item["sha256"], f"source hash {key}")
        require(digest == payload["source_sha256"][key], f"source map {key}")
    return len(payload["source_audit"])


def check_row(row: dict) -> dict[str, mp.mpf]:
    n_value = int(row["N"])
    ell = mp.mpf(row["chart"]["L"])
    time = 1 / (2 * ell**2)
    x = 4 * mp.pi * mp.exp(ell)
    state = chart(x, time, n_value)
    require(int(mp.floor(mp.re(state["a"]))) == n_value, "fixed-N chart")

    c_at_1 = normalized_correction(x, time, mp.mpf(1))
    c_at_2 = normalized_correction(x, time, mp.mpf(2))
    rho_2 = (c_at_2 - 2 * c_at_1 + 1) / 2
    rho_1 = c_at_1 - 1 - rho_2
    d_at_1 = scaled_diff(lambda xx: normalized_correction(xx, time, mp.mpf(1)), x)
    d_at_2 = scaled_diff(lambda xx: normalized_correction(xx, time, mp.mpf(2)), x)
    rho_2_x = (d_at_2 - 2 * d_at_1) / 2
    rho_1_x = d_at_1 - rho_2_x
    s_prime = scaled_diff(lambda xx: s_star(xx, time), x)
    s_second = scaled_diff(lambda xx: s_star(xx, time), x, 2)
    u_n_x = scaled_diff(lambda xx: mp.log(chart(xx, time, n_value)["a"] / n_value), x)
    chi_rate = terminal_chi(x, time, n_value)

    inputs = row["compiler_inputs"]
    tolerance = mp.mpf("1e-100")
    scalar_deltas = [
        abs(rho_1 - decode(inputs["rho_1"])),
        abs(rho_2 - decode(inputs["rho_2"])),
        abs(rho_1_x - decode(inputs["rho_1_x"])),
        abs(rho_2_x - decode(inputs["rho_2_x"])),
        abs(mp.re(s_prime) - mp.mpf(inputs["c"])),
        abs(mp.im(s_prime) - mp.mpf(inputs["b"])),
        abs(mp.re(s_second) - mp.mpf(inputs["c_x"])),
        abs(mp.im(s_second) - mp.mpf(inputs["b_x"])),
        abs(u_n_x - mp.mpf(inputs["u_N_x"])),
        abs(chi_rate - decode(inputs["chi_rate"])),
    ]
    require(max(scalar_deltas) < tolerance, f"scalar adapter row N={n_value}")

    v_e = edge_value(x, time, n_value)
    q_e = scaled_diff(lambda xx: edge_value(xx, time, n_value), x)
    a_e = edge_scalar(x, time, n_value)
    n_e = scaled_diff(lambda xx: edge_scalar(xx, time, n_value), x)
    edge = row["omega_E"]
    edge_deltas = [
        abs(v_e - mp.mpf(edge["V_E"])),
        abs(n_e - mp.mpf(edge["N_E"])),
        abs(a_e - mp.mpf(edge["A_E"])),
        abs(q_e - mp.mpf(edge["Q_E"])),
    ]
    require(max(edge_deltas) < mp.mpf("1e-95"), f"direct edge derivatives N={n_value}")

    table: dict[str, mp.mpf | mp.mpc] = {
        "log_N": mp.log(n_value),
        "log_a": mp.log(state["a"]),
        "u_N": mp.log(state["a"] / n_value),
        "u_N_x": u_n_x,
        "rho_1": rho_1,
        "rho_2": rho_2,
        "rho_1_x": rho_1_x,
        "rho_2_x": rho_2_x,
        "s_prime": s_prime,
        "s_second": s_second,
        "b": mp.im(s_prime),
        "b_x": mp.im(s_second),
        "chi_rate": chi_rate,
    }
    rows = independent_rows(table)
    maximum_polynomial = mp.mpf(0)
    maximum_lift = mp.mpf(0)
    stored_base = row["polynomials"]["base_rows"]
    stored_tangent = row["polynomials"]["tangent_rows"]
    tangents: dict[str, list[mp.mpc]] = {}
    for name in ROW_ORDER:
        actual = rows[name]
        expected = [decode(value) for value in stored_base[name]]
        maximum_polynomial = max(maximum_polynomial, maximum_delta(actual, expected))
        tangents[name] = lift(actual, table["log_a"])
        expected_lift = [decode(value) for value in stored_tangent[name]]
        maximum_lift = max(maximum_lift, maximum_delta(tangents[name], expected_lift))
    require(maximum_polynomial < mp.mpf("1e-90"), f"polynomial rows N={n_value}")
    require(maximum_lift < mp.mpf("1e-89"), f"tangent rows N={n_value}")

    beta_e = (n_e, v_e, -q_e, -a_e)
    p_e = [mp.mpc(0) for _ in range(6)]
    for weight, name in zip(beta_e, ROW_ORDER):
        p_e = polynomial_add(p_e, polynomial_scale(weight, tangents[name]))
    expected_p_e = [decode(value) for value in row["polynomials"]["P_E"]]
    p_e_delta = maximum_delta(p_e, expected_p_e)
    require(p_e_delta < mp.mpf("1e-88"), f"edge-affine polynomial N={n_value}")
    require(row["polynomials"]["P_E_degree"] == 5, "edge-affine degree")

    for value in row["edge_bound_ratios"].values():
        require(mp.mpf(value) < 1, f"edge bound ratio N={n_value}")
    require(abs(chi_rate + mp.j * state["h"] * state["theta"] / 2) < 3 * state["h"] ** 2, "physical chi bound")
    return {
        "scalar": max(scalar_deltas),
        "edge": max(edge_deltas),
        "polynomial": maximum_polynomial,
        "lift": maximum_lift,
        "P_E": p_e_delta,
    }


def check_note(payload: dict) -> None:
    note = NOTE_PATH.read_text(encoding="utf-8")
    markers = (
        "fixed heat time t and fixed cutoff N",
        "genuine endpoint alone",
        "degree-five edge-affine polynomial",
        "Every pi is inherited",
        "not interval enclosures",
    )
    for marker in markers:
        require(marker in note, f"missing note marker: {marker}")
    require(payload["proof_boundary"] in note, "proof boundary note")


def main() -> int:
    set_below_normal_priority()
    mp.mp.dps = CHECK_DPS
    require(RESULT_PATH.is_file(), "missing result")
    require(NOTE_PATH.is_file(), "missing note")
    payload = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    require(payload["kind"] == KIND, "kind")
    require(payload["schema_version"] == 1, "schema version")
    source_audits = check_sources(payload)
    ladder = payload["precision_ladder"]
    require([item["dps"] for item in ladder] == [90, 130, 170], "precision ladder")
    require(all(item["physical_rows"] == 3 for item in ladder), "ladder row count")
    require(len(payload["convergence"]["comparisons"]) == 2, "convergence count")
    for comparison in payload["convergence"]["comparisons"]:
        require(
            mp.mpf(comparison["maximum_absolute_delta"])
            < mp.mpf(comparison["required_lt"]),
            "stored convergence",
        )
    require(payload["convergence"]["stable_roster"], "stable roster")

    physical_rows = payload["physical_rows"]
    require([row["N"] for row in physical_rows] == [72_004_899_338, 72_004_899_339, 72_004_899_340], "physical roster")
    diagnostics = [check_row(row) for row in physical_rows]
    check_note(payload)

    rows = payload["rows"]
    require(len(rows) == 15, "gate row count")
    require(len({row["id"] for row in rows}) == 15, "gate row ids")
    require(sum(row["readiness"] == "open" for row in rows) == 1, "open gate rows")
    require(all(row["claim"] and row["certificate"] and row["proof_boundary"] for row in rows), "gate row fields")

    summary = payload["summary"]
    require(summary["rows"] == 15 and summary["issues"] == 0, "summary rows")
    require(summary["source_audits"] == source_audits == 6, "summary sources")
    require(summary["physical_roots"] == 3, "summary roots")
    require(summary["physical_scalar_rows"] == 3, "summary scalar rows")
    require(summary["genuine_edge_rows"] == 3, "summary edge rows")
    require(summary["compiled_base_polynomials"] == 12, "summary base rows")
    require(summary["compiled_tangent_polynomials"] == 12, "summary tangent rows")
    require(summary["edge_affine_polynomials"] == 3, "summary P_E")
    require(summary["physical_alpha_rows"] == 0, "summary alpha boundary")
    require(summary["physical_beta_p_rows"] == 0, "summary beta boundary")
    require(summary["signed_determinant_rows"] == 0, "summary sign boundary")
    require(payload["status"].startswith("physical scalar and genuine-edge adapter validated"), "status")
    require("not interval enclosures" in payload["proof_boundary"], "proof boundary")

    maxima = {
        key: max(row[key] for row in diagnostics)
        for key in ("scalar", "edge", "polynomial", "lift", "P_E")
    }
    print(
        "validated Newman C1 physical scalar and genuine-edge adapter gate: "
        f"{len(rows)} rows, {source_audits} source audits, 3 physical roots, "
        "10 scalar inputs per root, 2 independent chi routes, "
        "4 direct edge-derivative checks per root, 12 base and 12 tangent rows, "
        "3 degree-five edge-affine polynomials, "
        f"max scalar delta {mp.nstr(maxima['scalar'], 8)}, "
        f"max edge delta {mp.nstr(maxima['edge'], 8)}, "
        "0 retained-observation rows and 0 signed determinant rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
