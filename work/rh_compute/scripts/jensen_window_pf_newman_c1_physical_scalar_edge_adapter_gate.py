#!/usr/bin/env python3
"""Build physical scalar and genuine-edge rows at the calibrated C1 roots."""

from __future__ import annotations

import ctypes
from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_physical_scalar_edge_adapter_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"
PRECISION_LADDER = (90, 130, 170)
ROW_ORDER = ("P_V", "P_N", "P_A", "P_Q")

SOURCE_PATHS = {
    "calibrated_roots": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_calibrated_roster_interior_phase_root_scout.json",
    "finite_height_edge": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "q1_finite_height_real_edge_remainder_gate.json"
    ),
    "correction_envelope": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "correction_coefficient_envelope_gate.json"
    ),
    "polynomial_compiler": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_physical_observation_polynomial_compiler_gate.json",
    "edge_ownership": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_edge_affine_ownership_gate.json",
    "complex_endpoint": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "complex_endpoint_source_normalization_gate.json"
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


def load_and_audit_sources() -> tuple[dict[str, dict], dict[str, dict[str, str]]]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        require(path.is_file(), f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))

    require(payloads["calibrated_roots"]["summary"]["calibrated_roots"] == 3, "root count")
    require(
        payloads["finite_height_edge"]["counts"]["uniform_q1_finite_height_edge_signs"] == 1,
        "finite-height edge source drifted",
    )
    require(
        payloads["correction_envelope"]["counts"]["centered_coefficients"] == 6,
        "correction envelope drifted",
    )
    require(
        payloads["polynomial_compiler"]["summary"]["base_polynomials"] == 4,
        "compiler row count drifted",
    )
    require(
        payloads["edge_ownership"]["summary"]["edge_only_ownership_corrections"] == 1,
        "edge ownership drifted",
    )
    require(
        payloads["complex_endpoint"]["summary"]["corrected_cartesian_projections"] == 4,
        "complex endpoint source drifted",
    )
    audit = {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }
    return payloads, audit


def alpha_jet(s: mp.mpc) -> tuple[mp.mpc, mp.mpc, mp.mpc]:
    alpha = 1 / (2 * s) + 1 / (s - 1) + mp.log(s / (2 * mp.pi)) / 2
    alpha_one = -1 / (2 * s**2) - 1 / (s - 1) ** 2 + 1 / (2 * s)
    alpha_two = 1 / s**3 + 2 / (s - 1) ** 3 - 1 / (2 * s**2)
    return alpha, alpha_one, alpha_two


def correction(
    s: mp.mpc,
    time: mp.mpf,
    alpha: mp.mpc,
    alpha_one: mp.mpc,
    alpha_two: mp.mpc,
    log_n: mp.mpf,
) -> tuple[mp.mpc, mp.mpc]:
    w = alpha - log_n
    block = time / 4 + time**2 * w**2 / 8
    value = 1 / (6 * s) + alpha_one * block
    derivative = (-mp.j / 2) * (
        -1 / (6 * s**2)
        + alpha_two * block
        + time**2 * w * alpha_one**2 / 4
    )
    return value, derivative


def c0(p: mp.mpf | mp.mpc) -> mp.mpc:
    return (
        mp.exp(mp.j * mp.pi * (p**2 / 2 + mp.mpf(3) / 8))
        - mp.j * mp.sqrt(2) * mp.cos(mp.pi * p / 2)
    ) / (2 * mp.cos(mp.pi * p))


def c0_jet(p: mp.mpf) -> list[mp.mpc]:
    return [mp.diff(c0, p, order, addprec=50) for order in range(6)]


def vector_add(*vectors: list[mp.mpc]) -> list[mp.mpc]:
    return [sum(entries) for entries in zip(*vectors)]


def vector_scale(scale: mp.mpc, vector: list[mp.mpc]) -> list[mp.mpc]:
    return [scale * entry for entry in vector]


def lift(coefficients: list[mp.mpc], log_a: mp.mpf) -> list[mp.mpc]:
    result = [mp.mpc(0) for _ in range(len(coefficients) + 1)]
    for index, coefficient in enumerate(coefficients):
        result[index] -= mp.j * log_a * coefficient
        result[index + 1] += mp.j * coefficient
    return result


def trim_degree(coefficients: list[mp.mpc], tolerance: mp.mpf) -> int:
    for index in range(len(coefficients) - 1, -1, -1):
        if abs(coefficients[index]) > tolerance:
            return index
    return -1


def compile_rows(state: dict[str, mp.mpf | mp.mpc]) -> dict[str, list[mp.mpc]]:
    log_n = state["log_N"]
    rho_1 = state["rho_1"]
    rho_2 = state["rho_2"]
    rho_1_x = state["rho_1_x"]
    rho_2_x = state["rho_2_x"]
    s_prime = state["s_prime"]
    s_second = state["s_second"]
    chi_rate = state["chi_rate"]
    b = state["b"]
    b_x = state["b_x"]
    u_n = state["u_N"]
    u_n_x = state["u_N_x"]

    v_0 = [1, rho_1, rho_2, 0, 0]
    v_1 = [
        log_n,
        log_n * rho_1 - 1,
        log_n * rho_2 - rho_1,
        -rho_2,
        0,
    ]
    v_2 = [
        log_n**2,
        log_n**2 * rho_1 - 2 * log_n,
        log_n**2 * rho_2 - 2 * log_n * rho_1 + 1,
        rho_1 - 2 * log_n * rho_2,
        rho_2,
    ]
    e_0 = [0, rho_1_x, rho_2_x, 0, 0]
    e_1 = [0, log_n * rho_1_x, log_n * rho_2_x - rho_1_x, -rho_2_x, 0]
    q_0 = vector_add(
        vector_scale(chi_rate, v_0), vector_scale(s_prime, v_1), e_0
    )
    q_1 = vector_add(
        vector_scale(chi_rate, v_1), vector_scale(s_prime, v_2), e_1
    )
    a_0 = vector_add(
        vector_scale(s_prime, v_1), vector_scale(mp.j * b * u_n, v_0)
    )
    n_0 = vector_add(
        vector_scale(s_second, v_1),
        vector_scale(s_prime, q_1),
        vector_scale(mp.j * (b_x * u_n + b * u_n_x), v_0),
        vector_scale(mp.j * b * u_n, q_0),
    )
    return {"P_V": v_0, "P_N": n_0, "P_A": a_0, "P_Q": q_0}


def physical_state(ell: mp.mpf, expected_n: int) -> dict[str, object]:
    time = 1 / (2 * ell**2)
    x = 4 * mp.pi * mp.exp(ell)
    s = (1 - mp.j * x) / 2
    a = mp.sqrt(mp.exp(ell) + time / 16)
    n_value = int(mp.floor(a))
    require(n_value == expected_n, "adapter root left fixed-N cell")
    h = 1 / a
    theta = a - n_value
    p = 1 - 2 * theta
    log_n = mp.log(n_value)
    log_a = mp.log(a)
    u_n = log_a - log_n
    t_0 = 2 * mp.pi / h**2

    alpha, alpha_one, alpha_two = alpha_jet(s)
    d_1, d_1_x = correction(s, time, alpha, alpha_one, alpha_two, mp.mpf(0))
    d_n, d_n_x = correction(s, time, alpha, alpha_one, alpha_two, log_n)
    a_0 = 1 + d_1
    correction_b = alpha_one * time**2 / 8
    correction_b_x = -mp.j * alpha_two * time**2 / 16
    rho_2 = correction_b / a_0
    rho_1 = -2 * alpha * rho_2
    rho_2_x = correction_b_x / a_0 - correction_b * d_1_x / a_0**2
    alpha_x = -mp.j * alpha_one / 2
    rho_1_x = -2 * alpha_x * rho_2 - 2 * alpha * rho_2_x

    s_prime = -mp.j / 2 - mp.j * time * alpha_one / 4
    s_second = -time * alpha_two / 8
    c = mp.re(s_prime)
    b = mp.im(s_prime)
    c_x = mp.re(s_second)
    b_x = mp.im(s_second)
    u_n_x = h**2 / (8 * mp.pi)

    m_rate = (-mp.j / 2) * (alpha + time * alpha * alpha_one / 2)
    beta_rate = -mp.pi / 8 + 1 / (4 * t_0)
    t_0_rate = 1 / (2 * t_0)
    chi_rate = (
        m_rate
        + d_1_x / (1 + d_1)
        - beta_rate
        - t_0_rate
        - s_prime * log_n
    )

    stable_ell = mp.log(1 - h * theta)
    saddle_rho = time * h**2 / 16
    saddle_r = 1 - saddle_rho
    y = h**2 / (4 * mp.pi * saddle_r)
    chi_r = mp.log(saddle_r) / 2 + mp.log(1 + y**2) / 4 - y**2 / (1 + y**2)
    eta_small = mp.atan(y) / 2 + 3 * y / (1 + y**2)
    w = alpha - log_a - stable_ell
    w_x = (
        eta_small / 2
        - 3 / (4 * t_0)
        + mp.j * (stable_ell + h * theta - chi_r) / 2
        - mp.j * (time / 4) * alpha_one * w
        + d_n_x / (1 + d_n)
    )
    delta_n = d_n_x / (1 + d_n) - d_1_x / (1 + d_1)
    chi_terminal = -mp.j * h * theta / 2 + w_x - delta_n

    f = c0_jet(p)
    c_1 = f[3] / (12 * mp.pi**2)
    c_1_p = f[4] / (12 * mp.pi**2)
    c_1_pp = f[5] / (12 * mp.pi**2)
    p_x = -h / (4 * mp.pi)
    h_x = -h**3 / (8 * mp.pi)
    p_xx = h**3 / (32 * mp.pi**2)
    h_xx = 3 * h**5 / (64 * mp.pi**2)
    endpoint_h = f[0] + h * c_1
    endpoint_h_x = p_x * (f[1] + h * c_1_p) + h_x * c_1
    endpoint_h_xx = (
        p_xx * (f[1] + h * c_1_p)
        + p_x**2 * (f[2] + h * c_1_pp)
        + 2 * p_x * h_x * c_1_p
        + h_xx * c_1
    )

    saddle_chi = alpha - log_a
    k_rate = -mp.pi / 8 + 1 / (4 * t_0) + 1 / (2 * (t_0 + mp.j))
    mu = k_rate + mp.j * saddle_chi / 2 + mp.j * (time / 4) * alpha_one * saddle_chi
    k_rate_x = -1 / (8 * t_0**2) - 1 / (4 * (t_0 + mp.j) ** 2)
    saddle_chi_x = -mp.j * alpha_one / 2 - h**2 / (8 * mp.pi)
    alpha_one_x = -mp.j * alpha_two / 2
    mu_x = (
        k_rate_x
        + mp.j * saddle_chi_x / 2
        + mp.j
        * (time / 4)
        * (alpha_one_x * saddle_chi + alpha_one * saddle_chi_x)
    )
    endpoint_j = endpoint_h_x + mu * endpoint_h
    endpoint_j_x = endpoint_h_xx + mu_x * endpoint_h + mu * endpoint_h_x

    tau = 1 + mp.j / t_0
    tau_x = -mp.j / (2 * t_0**2)
    alpha_c = c * u_n
    alpha_c_x = c_x * u_n + c * u_n_x
    v_e = mp.re(tau * endpoint_h)
    q_e = mp.re(tau_x * endpoint_h + tau * endpoint_h_x)
    a_e = mp.re(tau * endpoint_j) - alpha_c * v_e
    n_e = (
        mp.re(tau_x * endpoint_j + tau * endpoint_j_x)
        - alpha_c_x * v_e
        - alpha_c * q_e
    )

    compiler_state: dict[str, mp.mpf | mp.mpc] = {
        "log_N": log_n,
        "log_a": log_a,
        "u_N": u_n,
        "u_N_x": u_n_x,
        "rho_1": rho_1,
        "rho_2": rho_2,
        "rho_1_x": rho_1_x,
        "rho_2_x": rho_2_x,
        "s_prime": s_prime,
        "s_second": s_second,
        "c": c,
        "b": b,
        "c_x": c_x,
        "b_x": b_x,
        "chi_rate": chi_rate,
    }
    rows = compile_rows(compiler_state)
    tangent_rows = {name: lift(row, log_a) for name, row in rows.items()}
    beta_e = (n_e, v_e, -q_e, -a_e)
    p_e = [mp.mpc(0) for _ in range(6)]
    for weight, name in zip(beta_e, ROW_ORDER):
        p_e = vector_add(p_e, vector_scale(weight, tangent_rows[name]))

    return {
        "N": n_value,
        "L": ell,
        "t": time,
        "x": x,
        "s": s,
        "a": a,
        "h": h,
        "theta": theta,
        "p": p,
        "T_0": t_0,
        "log_N": log_n,
        "log_a": log_a,
        "u_N": u_n,
        "u_N_x": u_n_x,
        "alpha": alpha,
        "alpha_prime": alpha_one,
        "alpha_second": alpha_two,
        "d_1": d_1,
        "d_1_x": d_1_x,
        "d_N": d_n,
        "d_N_x": d_n_x,
        "rho_1": rho_1,
        "rho_2": rho_2,
        "rho_1_x": rho_1_x,
        "rho_2_x": rho_2_x,
        "s_prime": s_prime,
        "s_second": s_second,
        "c": c,
        "b": b,
        "c_x": c_x,
        "b_x": b_x,
        "chi_rate": chi_rate,
        "chi_terminal": chi_terminal,
        "chi_delta": abs(chi_rate - chi_terminal),
        "endpoint_H": endpoint_h,
        "endpoint_H_x": endpoint_h_x,
        "endpoint_H_xx": endpoint_h_xx,
        "endpoint_J": endpoint_j,
        "endpoint_J_x": endpoint_j_x,
        "V_E": v_e,
        "N_E": n_e,
        "A_E": a_e,
        "Q_E": q_e,
        "rows": rows,
        "tangent_rows": tangent_rows,
        "P_E": p_e,
    }


def number_text(value: mp.mpf, digits: int) -> str:
    return mp.nstr(value, max(35, digits), strip_zeros=False)


def complex_payload(value: mp.mpc, digits: int) -> dict[str, str]:
    return {
        "re": number_text(mp.re(value), digits),
        "im": number_text(mp.im(value), digits),
    }


def encode_vector(values: list[mp.mpc], digits: int) -> list[dict[str, str]]:
    return [complex_payload(value, digits) for value in values]


def encode_state(state: dict[str, object], digits: int) -> dict:
    complex_names = (
        "s",
        "alpha",
        "alpha_prime",
        "alpha_second",
        "d_1",
        "d_1_x",
        "d_N",
        "d_N_x",
        "rho_1",
        "rho_2",
        "rho_1_x",
        "rho_2_x",
        "s_prime",
        "s_second",
        "chi_rate",
        "chi_terminal",
        "endpoint_H",
        "endpoint_H_x",
        "endpoint_H_xx",
        "endpoint_J",
        "endpoint_J_x",
    )
    real_names = (
        "L",
        "t",
        "x",
        "a",
        "h",
        "theta",
        "p",
        "T_0",
        "log_N",
        "log_a",
        "u_N",
        "u_N_x",
        "c",
        "b",
        "c_x",
        "b_x",
        "chi_delta",
        "V_E",
        "N_E",
        "A_E",
        "Q_E",
    )
    result = {"N": state["N"]}
    result["chart"] = {
        name: number_text(state[name], digits)
        for name in (
            "L",
            "t",
            "x",
            "a",
            "h",
            "theta",
            "p",
            "T_0",
            "log_N",
            "log_a",
            "u_N",
        )
    }
    result["compiler_inputs"] = {
        name: complex_payload(state[name], digits)
        for name in ("rho_1", "rho_2", "rho_1_x", "rho_2_x", "chi_rate")
    }
    result["compiler_inputs"].update(
        {
            name: number_text(state[name], digits)
            for name in ("c", "b", "c_x", "b_x", "u_N_x")
        }
    )
    result["rates"] = {
        "s_prime": complex_payload(state["s_prime"], digits),
        "s_second": complex_payload(state["s_second"], digits),
        "chi_carrier_route": complex_payload(state["chi_rate"], digits),
        "chi_terminal_route": complex_payload(state["chi_terminal"], digits),
        "chi_route_delta": number_text(state["chi_delta"], digits),
        "chi_route_delta_over_h2": number_text(
            state["chi_delta"] / state["h"] ** 2, digits
        ),
    }
    result["correction"] = {
        name: complex_payload(state[name], digits)
        for name in ("d_1", "d_1_x", "d_N", "d_N_x")
    }
    result["endpoint_jets"] = {
        name: complex_payload(state[name], digits)
        for name in (
            "endpoint_H",
            "endpoint_H_x",
            "endpoint_H_xx",
            "endpoint_J",
            "endpoint_J_x",
        )
    }
    result["omega_E"] = {
        "order": ["V_E", "N_E", "A_E", "Q_E"],
        **{
            name: number_text(state[name], digits)
            for name in ("V_E", "N_E", "A_E", "Q_E")
        },
    }
    result["edge_bound_ratios"] = {
        "abs_V_E_over_3_5": number_text(abs(state["V_E"]) / mp.mpf("0.6"), digits),
        "abs_A_E_over_4h": number_text(abs(state["A_E"]) / (4 * state["h"]), digits),
        "abs_Q_E_over_4h": number_text(abs(state["Q_E"]) / (4 * state["h"]), digits),
        "abs_N_E_over_4h2": number_text(abs(state["N_E"]) / (4 * state["h"] ** 2), digits),
    }
    result["polynomials"] = {
        "coefficient_order": "ascending powers of ell=log n",
        "observation_order": ["V", "mathcal N", "A", "Q"],
        "base_rows": {
            name: encode_vector(state["rows"][name], digits) for name in ROW_ORDER
        },
        "tangent_rows": {
            name: encode_vector(state["tangent_rows"][name], digits)
            for name in ROW_ORDER
        },
        "P_E": encode_vector(state["P_E"], digits),
        "P_E_degree": trim_degree(state["P_E"], mp.power(10, -(digits - 25))),
    }
    require(set(complex_names).issubset(state), "complex serialization fields")
    require(set(real_names).issubset(state), "real serialization fields")
    return result


def root_inputs(root_payload: dict) -> list[tuple[int, str]]:
    high = root_payload["precision_ladder"][-1]
    rows = [high, *root_payload["additional_roots"]]
    result = sorted((int(row["N"]), str(row["root"]["L"])) for row in rows)
    require([n for n, _ in result] == [72_004_899_338, 72_004_899_339, 72_004_899_340], "root roster")
    return result


def flatten_diagnostic(state: dict[str, object]) -> list[mp.mpc]:
    values: list[mp.mpc] = []
    for name in (
        "rho_1",
        "rho_2",
        "rho_1_x",
        "rho_2_x",
        "s_prime",
        "s_second",
        "chi_rate",
        "chi_terminal",
        "V_E",
        "N_E",
        "A_E",
        "Q_E",
    ):
        values.append(mp.mpc(state[name]))
    for name in ROW_ORDER:
        values.extend(state["rows"][name])
        values.extend(state["tangent_rows"][name])
    values.extend(state["P_E"])
    return values


def compute_ladder(
    roots: list[tuple[int, str]],
) -> tuple[list[dict], list[dict], dict]:
    numerical: dict[int, list[dict[str, object]]] = {}
    ladder_payload: list[dict] = []
    for dps in PRECISION_LADDER:
        mp.mp.dps = dps
        states = [physical_state(mp.mpf(ell), n_value) for n_value, ell in roots]
        numerical[dps] = states
        max_chi = max(state["chi_delta"] / state["h"] ** 2 for state in states)
        max_edge_ratio = max(
            max(
                abs(state["V_E"]) / mp.mpf("0.6"),
                abs(state["A_E"]) / (4 * state["h"]),
                abs(state["Q_E"]) / (4 * state["h"]),
                abs(state["N_E"]) / (4 * state["h"] ** 2),
            )
            for state in states
        )
        ladder_payload.append(
            {
                "dps": dps,
                "physical_rows": len(states),
                "N_values": [state["N"] for state in states],
                "maximum_chi_route_delta_over_h2": number_text(max_chi, dps - 10),
                "maximum_edge_bound_ratio": number_text(max_edge_ratio, dps - 10),
            }
        )

    convergence: list[dict] = []
    for lower, upper in zip(PRECISION_LADDER, PRECISION_LADDER[1:]):
        mp.mp.dps = upper + 20
        maximum = mp.mpf(0)
        for lower_state, upper_state in zip(numerical[lower], numerical[upper]):
            lower_values = flatten_diagnostic(lower_state)
            upper_values = flatten_diagnostic(upper_state)
            require(len(lower_values) == len(upper_values), "diagnostic vector length")
            maximum = max(
                maximum,
                max(abs(left - right) for left, right in zip(lower_values, upper_values)),
            )
        threshold = mp.power(10, -(lower - 35))
        require(maximum < threshold, f"{lower}/{upper} precision convergence")
        convergence.append(
            {
                "dps_pair": [lower, upper],
                "maximum_absolute_delta": number_text(maximum, 55),
                "required_lt": number_text(threshold, 55),
            }
        )

    high_states = numerical[PRECISION_LADDER[-1]]
    mp.mp.dps = PRECISION_LADDER[-1]
    for state in high_states:
        require(state["chi_delta"] < mp.power(10, -120), "chi route agreement")
        require(abs(state["V_E"]) < mp.mpf("0.6"), "V_E bound")
        require(abs(state["A_E"]) < 4 * state["h"], "A_E bound")
        require(abs(state["Q_E"]) < 4 * state["h"], "Q_E bound")
        require(abs(state["N_E"]) < 4 * state["h"] ** 2, "N_E bound")
        require(trim_degree(state["P_E"], mp.mpf("1e-130")) == 5, "edge polynomial degree")
    encoded = [encode_state(state, 145) for state in high_states]
    return encoded, ladder_payload, {"comparisons": convergence, "stable_roster": True}


def exact_contract() -> dict[str, str]:
    return {
        "fixed_derivative": "Every x derivative is taken at fixed heat time t and fixed cutoff N. The immutable L root selects the state but is not differentiated along the q=1 root curve.",
        "correction": "With B=alpha'(s)t^2/8 and a_0=1+d_1, rho_2=B/a_0, rho_1=-2alpha rho_2, rho_(2,x)=B_x/a_0-B d_(1,x)/a_0^2, and rho_(1,x)=-2alpha_xrho_2-2alpha rho_(2,x).",
        "rates": "s_*'=-i/2-it alpha'/4, s_*''=-t alpha''/8, and u_(N,x)=h^2/(8pi).",
        "chi": "The carrier-normalizer route chi_rate=(M_(t,x)/M_t)+d_(1,x)/(1+d_1)-beta_x/beta-T_(0,x)/T_0-s_*'log N agrees with the terminal route -ih theta/2+W_x-delta_N.",
        "edge": "After the terminal carrier is moved into the full-support vector, omega_E is the genuine endpoint alone: V_E=Re(tau H), Q_E=Re(tau_xH+tau H_x), A_E=Re(tau J)-alpha_cV_E, and N_E=Re(tau_xJ+tau J_x)-alpha_(c,x)V_E-alpha_cQ_E.",
        "compiler": "The twelve physical scalar values compile P_V,P_N,P_A,P_Q in ascending ell powers and their exact tangent lifts. beta_E=(N_E,V_E,-Q_E,-A_E) then compiles the degree-five edge-affine polynomial P_E.",
        "pi_provenance": "Every pi is inherited from the completed-zeta/Riemann-Siegel source: x=4pi exp(L), T_0=2pi a^2, the C_0 phase, p_x=-h/(4pi), and u_(N,x)=h^2/(8pi). The adapter introduces no polygonal or fitted pi.",
        "missing": "The retained carrier observations alpha and beta_p remain unevaluated. They require the stationary diagonal, oriented adjacent recurrence, and common-cutoff near/paired-remote compression; direct enumeration to N is still rejected.",
        "proof_boundary": "These are reproducible high-precision point values with cross-precision and independent-formula checks, not interval enclosures. They prove no determinant sign, all-cell statement, complete-current inequality, Lambda<=0, RH, or prize-level conclusion.",
    }


def build_rows(contract: dict[str, str], physical_rows: list[dict]) -> list[GateRow]:
    n_values = [row["N"] for row in physical_rows]
    max_chi = max(mp.mpf(row["rates"]["chi_route_delta_over_h2"]) for row in physical_rows)
    max_edge = max(
        max(mp.mpf(value) for value in row["edge_bound_ratios"].values())
        for row in physical_rows
    )
    return [
        GateRow("pse_01_roots", "immutable roots", "diagnostic_validated", "All three calibrated physical roots are consumed without re-solving them.", f"N={n_values}.", "The root scout remains numerical, not interval-certified."),
        GateRow("pse_02_derivative", "derivative convention", "guard_validated", "The adapter uses the source fixed-t, fixed-N x derivative.", contract["fixed_derivative"], "No derivative along the q=1 L curve is substituted."),
        GateRow("pse_03_correction", "correction quotient", "diagnostic_validated", "The four normalized correction coefficients are evaluated from the exact quotient.", contract["correction"], "No asymptotic replacement of a_0 is made."),
        GateRow("pse_04_rates", "carrier rates", "diagnostic_validated", "The physical first and second carrier rates and u_(N,x) are explicit.", contract["rates"], "The printed phase angle chi_N is not used."),
        GateRow("pse_05_chi_carrier", "carrier chi route", "diagnostic_validated", "chi_rate is evaluated from the absolute source phase and real normalizer.", contract["chi"], "Large terms cancel at high precision before serialization."),
        GateRow("pse_06_chi_terminal", "terminal chi route", "diagnostic_validated", "The stable terminal W_x identity independently reproduces chi_rate.", f"Maximum |Delta chi|/h^2={number_text(max_chi, 35)}.", "Agreement of floating values is not an interval identity."),
        GateRow("pse_07_endpoint", "complex endpoint jets", "diagnostic_validated", "C_0 through its fifth derivative supplies H,H_x,H_xx,J,J_x.", "The complete first endpoint correction and every source rate are retained.", "No real-H specialization is made."),
        GateRow("pse_08_owner", "genuine edge owner", "guard_validated", "Only the genuine endpoint remains in omega_E.", contract["edge"], "The moved terminal carrier is not counted twice."),
        GateRow("pse_09_edge", "physical edge rows", "diagnostic_validated", "All four genuine-edge coordinates are evaluated at all three roots.", f"Maximum certified-box ratio observed={number_text(max_edge, 35)}<1.", "The source bounds are inherited; point values are not interval enclosures."),
        GateRow("pse_10_compile", "physical base rows", "diagnostic_validated", "Four physical observation polynomials are compiled at every root.", contract["compiler"], "Retained carrier moments have not yet been projected."),
        GateRow("pse_11_lift", "physical tangent rows", "diagnostic_validated", "The four physical rows pass the exact tangent-lift interface.", "Each root stores four length-five base rows and four length-six lifts.", "Auxiliary-frequency differentiation holds source coefficients fixed."),
        GateRow("pse_12_affine", "edge-affine polynomial", "diagnostic_validated", "The single-owner beta_E row produces a degree-five physical P_E at every root.", "P_E=i(ell-log a)F_(beta_E) is stored in ascending powers.", "No signed arithmetic functional is inferred from its coefficients."),
        GateRow("pse_13_pi", "pi provenance", "guard_validated", "Every occurrence of pi has an explicit source origin.", contract["pi_provenance"], "No geometric analogy supplies a constant."),
        GateRow("pse_14_remaining", "retained observations", "open", "Compute alpha and beta_p through compressed actual-source observations.", contract["missing"], "Direct enumeration of roughly 7.2e10 carriers is inadmissible."),
        GateRow("pse_15_boundary", "proof boundary", "guard_validated", "The adapter is a physical numerical interface, not a determinant theorem.", contract["proof_boundary"], "All prize-level implications remain open."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    summary = payload["summary"]
    return "\n".join(
        [
            "# Newman C1 Physical Scalar and Genuine-Edge Adapter Gate",
            "",
            "Date: 2026-08-05",
            "",
            "Status: three physical scalar rows, genuine-edge rows, and edge-affine polynomials evaluated with independent route checks; retained observations and determinant sign remain open; not a proof of RH.",
            "",
            "## Derivative Convention",
            "",
            exact["fixed_derivative"],
            "",
            "## Scalar Adapter",
            "",
            exact["correction"],
            "",
            exact["rates"],
            "",
            exact["chi"],
            "",
            "## Genuine Edge",
            "",
            exact["edge"],
            "",
            "## Compiler Handoff",
            "",
            exact["compiler"],
            "",
            "## Validation",
            "",
            f"- Physical roots: `{summary['physical_roots']}`.",
            f"- Precision levels: `{summary['precision_levels']}`.",
            f"- Physical scalar rows: `{summary['physical_scalar_rows']}`.",
            f"- Genuine-edge rows: `{summary['genuine_edge_rows']}`.",
            f"- Compiled base polynomials: `{summary['compiled_base_polynomials']}`.",
            f"- Compiled tangent polynomials: `{summary['compiled_tangent_polynomials']}`.",
            f"- Degree-five edge-affine polynomials: `{summary['edge_affine_polynomials']}`.",
            "",
            "## Next Layer",
            "",
            exact["missing"],
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


def main() -> int:
    set_below_normal_priority()
    payloads, source_audit = load_and_audit_sources()
    roots = root_inputs(payloads["calibrated_roots"])
    physical_rows, measured_ladder, convergence = compute_ladder(roots)
    exact = exact_contract()
    rows = build_rows(exact, physical_rows)
    summary = {
        "rows": len(rows),
        "issues": 0,
        "source_audits": len(source_audit),
        "physical_roots": len(physical_rows),
        "precision_levels": len(PRECISION_LADDER),
        "physical_scalar_rows": len(physical_rows),
        "genuine_edge_rows": len(physical_rows),
        "compiled_base_polynomials": len(physical_rows) * len(ROW_ORDER),
        "compiled_tangent_polynomials": len(physical_rows) * len(ROW_ORDER),
        "edge_affine_polynomials": len(physical_rows),
        "independent_chi_routes": 2,
        "open_rows": sum(row.readiness == "open" for row in rows),
        "physical_alpha_rows": 0,
        "physical_beta_p_rows": 0,
        "signed_determinant_rows": 0,
    }
    success = (
        "built Newman C1 physical scalar and genuine-edge adapter gate: "
        f"{len(rows)} rows, 0 issues, {len(physical_rows)} physical roots, "
        f"{len(PRECISION_LADDER)} precision levels, {len(physical_rows)} scalar rows, "
        f"{len(physical_rows)} genuine-edge rows, "
        f"{len(physical_rows) * len(ROW_ORDER)} base and "
        f"{len(physical_rows) * len(ROW_ORDER)} tangent polynomials, "
        f"{len(physical_rows)} degree-five edge-affine polynomials, "
        "0 retained-observation rows and 0 signed determinant rows"
    )
    payload = {
        "kind": KIND,
        "schema_version": 1,
        "date": "2026-08-05",
        "status": "physical scalar and genuine-edge adapter validated; retained observations and determinant sign open",
        "source_sha256": {key: item["sha256"] for key, item in source_audit.items()},
        "source_audit": source_audit,
        "schema": {
            "x_derivative": "fixed t and fixed N",
            "polynomial_variable": "ell=log n",
            "coefficient_order": "ascending powers",
            "observation_order": ["V", "mathcal N", "A", "Q"],
            "edge_order": ["V_E", "N_E", "A_E", "Q_E"],
            "compiler_inputs": [
                "L_N",
                "log_a",
                "rho_1",
                "rho_2",
                "rho_1_x",
                "rho_2_x",
                "c",
                "b",
                "c_x",
                "b_x",
                "chi_rate",
                "u_N_x",
            ],
        },
        "exact": exact,
        "precision_ladder": [
            {
                **item,
                "root_source": source_audit["calibrated_roots"]["path"],
            }
            for item in measured_ladder
        ],
        "convergence": convergence,
        "physical_rows": physical_rows,
        "rows": [asdict(row) for row in rows],
        "summary": summary,
        "next_action": exact["missing"],
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
