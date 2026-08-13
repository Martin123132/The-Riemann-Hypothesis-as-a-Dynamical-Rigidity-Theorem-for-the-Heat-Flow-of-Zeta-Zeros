#!/usr/bin/env python3
"""Independently check the C1 physical-observation polynomial compiler gate."""

from __future__ import annotations

import ctypes
import hashlib
import json
import math
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_physical_observation_polynomial_compiler_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"
LOW_N_PATH = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_fixed_cell_low_n_determinant_evaluator_gate.json"
ROW_ORDER = ("P_V", "P_N", "P_A", "P_Q")
FUNCTIONALS = ("carrier", "band_blocks", "near", "remote")
TOLERANCE = 2.0e-12


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


def symbol_table() -> dict[str, sp.Symbol | sp.Expr]:
    ell, x = sp.symbols("ell x", real=True)
    log_n, log_a, u_n_x = sp.symbols("L_N log_a u_N_x", real=True)
    c, b, c_x, b_x = sp.symbols("c b c_x b_x", real=True)
    rho_1, rho_2, rho_1_x, rho_2_x, chi_rate = sp.symbols(
        "rho_1 rho_2 rho_1_x rho_2_x chi_rate"
    )
    names = (
        "N_p_xi",
        "V_p_xi",
        "Q_p_xi",
        "A_p_xi",
        "N_p",
        "V_p",
        "Q_p",
        "A_p",
        "N_E",
        "V_E",
        "Q_E",
        "A_E",
    )
    table: dict[str, sp.Symbol | sp.Expr] = {
        "ell": ell,
        "x": x,
        "L_N": log_n,
        "log_a": log_a,
        "u_N_x": u_n_x,
        "c": c,
        "b": b,
        "c_x": c_x,
        "b_x": b_x,
        "rho_1": rho_1,
        "rho_2": rho_2,
        "rho_1_x": rho_1_x,
        "rho_2_x": rho_2_x,
        "chi_rate": chi_rate,
        "I": sp.I,
    }
    table.update({name: sp.symbols(name, real=True) for name in names})
    return table


def pad(values: list[sp.Expr], size: int) -> list[sp.Expr]:
    require(len(values) <= size, "polynomial overflow")
    return [*values, *([sp.S.Zero] * (size - len(values)))]


def add(*rows: list[sp.Expr], size: int | None = None) -> list[sp.Expr]:
    length = size if size is not None else max(len(row) for row in rows)
    result = [sp.S.Zero] * length
    for row in rows:
        require(len(row) <= length, "addition overflow")
        for index, value in enumerate(row):
            result[index] = sp.expand(result[index] + value)
    return result


def scale(row: list[sp.Expr], scalar: sp.Expr, size: int | None = None) -> list[sp.Expr]:
    values = [sp.expand(scalar * value) for value in row]
    return pad(values, size) if size is not None else values


def multiply(
    left: list[sp.Expr], right: list[sp.Expr], size: int | None = None
) -> list[sp.Expr]:
    natural = len(left) + len(right) - 1
    length = natural if size is None else size
    require(natural <= length, "multiplication overflow")
    result = [sp.S.Zero] * length
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            result[left_index + right_index] = sp.expand(
                result[left_index + right_index] + left_value * right_value
            )
    return result


def rows_equal(left: list[sp.Expr], right: list[sp.Expr]) -> bool:
    return len(left) == len(right) and all(
        sp.expand(a_value - b_value) == 0 for a_value, b_value in zip(left, right)
    )


def direct_rows(table: dict[str, sp.Symbol | sp.Expr]) -> dict[str, list[sp.Expr]]:
    log_n = table["L_N"]
    log_a = table["log_a"]
    u_n = log_a - log_n
    c = table["c"]
    b = table["b"]
    c_x = table["c_x"]
    b_x = table["b_x"]
    s_prime = c + sp.I * b
    s_second = c_x + sp.I * b_x
    correction = [1, table["rho_1"], table["rho_2"]]
    correction_x = [0, table["rho_1_x"], table["rho_2_x"]]
    log_difference = [log_n, -1]
    rate = add(
        scale(log_difference, s_prime),
        [sp.I * b * u_n],
    )
    rate_x = add(
        scale(log_difference, s_second),
        [sp.I * (b_x * u_n + b * table["u_N_x"])],
    )
    q_factor = add([table["chi_rate"]], scale(log_difference, s_prime))
    q_row = add(multiply(q_factor, correction), correction_x, size=4)
    a_row = multiply(rate, correction, size=4)
    n_row = add(
        multiply(rate, q_row),
        multiply(rate_x, correction),
        size=5,
    )
    return {
        "P_V": pad(correction, 5),
        "P_N": n_row,
        "P_A": pad(a_row, 5),
        "P_Q": pad(q_row, 5),
    }


def moment_rows(table: dict[str, sp.Symbol | sp.Expr]) -> dict[str, list[sp.Expr]]:
    log_n = table["L_N"]
    log_a = table["log_a"]
    u_n = log_a - log_n
    c = table["c"]
    b = table["b"]
    c_x = table["c_x"]
    b_x = table["b_x"]
    s_prime = c + sp.I * b
    s_second = c_x + sp.I * b_x
    v_0 = [1, table["rho_1"], table["rho_2"]]
    e_0 = [0, table["rho_1_x"], table["rho_2_x"]]
    multiplier = [log_n, -1]
    v_1 = multiply(multiplier, v_0)
    v_2 = multiply(multiplier, v_1)
    e_1 = multiply(multiplier, e_0)
    q_0 = add(
        scale(v_0, table["chi_rate"]),
        scale(v_1, s_prime),
        e_0,
        size=5,
    )
    q_1 = add(
        scale(v_1, table["chi_rate"]),
        scale(v_2, s_prime),
        e_1,
        size=5,
    )
    a_0 = add(
        scale(v_1, s_prime),
        scale(v_0, sp.I * b * u_n),
        size=5,
    )
    n_0 = add(
        scale(v_1, s_second),
        scale(q_1, s_prime),
        scale(v_0, sp.I * (b_x * u_n + b * table["u_N_x"])),
        scale(q_0, sp.I * b * u_n),
        size=5,
    )
    return {
        "P_V": pad(v_0, 5),
        "P_N": n_0,
        "P_A": a_0,
        "P_Q": q_0,
    }


def lift(row: list[sp.Expr], log_a: sp.Expr) -> list[sp.Expr]:
    result = [sp.S.Zero] * (len(row) + 1)
    for index, value in enumerate(row):
        result[index] = sp.expand(result[index] - sp.I * log_a * value)
        result[index + 1] = sp.expand(result[index + 1] + sp.I * value)
    return result


def translate(row: list[sp.Expr], offset: sp.Expr) -> list[sp.Expr]:
    result = [sp.S.Zero] * len(row)
    for power, coefficient in enumerate(row):
        for centered_power in range(power + 1):
            result[centered_power] = sp.expand(
                result[centered_power]
                + coefficient
                * sp.binomial(power, centered_power)
                * offset ** (power - centered_power)
            )
    return result


def centered_rows(table: dict[str, sp.Symbol | sp.Expr]) -> dict[str, list[sp.Expr]]:
    log_n = table["L_N"]
    log_a = table["log_a"]
    u_n = log_a - log_n
    c = table["c"]
    b = table["b"]
    c_x = table["c_x"]
    b_x = table["b_x"]
    s_prime = c + sp.I * b
    s_second = c_x + sp.I * b_x
    correction = [
        1 + table["rho_1"] * log_a + table["rho_2"] * log_a**2,
        table["rho_1"] + 2 * table["rho_2"] * log_a,
        table["rho_2"],
    ]
    correction_x = [
        table["rho_1_x"] * log_a + table["rho_2_x"] * log_a**2,
        table["rho_1_x"] + 2 * table["rho_2_x"] * log_a,
        table["rho_2_x"],
    ]
    rate = [-c * u_n, -s_prime]
    rate_x = [-c_x * u_n + sp.I * b * table["u_N_x"], -s_second]
    delta = table["chi_rate"] - sp.I * b * u_n
    q_row = add(multiply(add(rate, [delta]), correction), correction_x, size=4)
    a_row = multiply(rate, correction, size=4)
    n_row = add(multiply(rate, q_row), multiply(rate_x, correction), size=5)
    return {
        "P_V": pad(correction, 5),
        "P_N": n_row,
        "P_A": pad(a_row, 5),
        "P_Q": pad(q_row, 5),
    }


def parse_row(values: list[str], table: dict[str, sp.Symbol | sp.Expr]) -> list[sp.Expr]:
    return [sp.sympify(value, locals=table) for value in values]


def parse_substitution(
    values: dict[str, str], table: dict[str, sp.Symbol | sp.Expr]
) -> dict[sp.Expr, sp.Expr]:
    return {
        sp.sympify(key, locals=table): sp.sympify(value, locals=table)
        for key, value in values.items()
    }


def weight_vectors(table: dict[str, sp.Symbol | sp.Expr]) -> dict[str, list[sp.Expr]]:
    return {
        "alpha": [table["N_p_xi"], table["V_p_xi"], -table["Q_p_xi"], -table["A_p_xi"]],
        "beta_p": [table["N_p"], table["V_p"], -table["Q_p"], -table["A_p"]],
        "beta_E": [table["N_E"], table["V_E"], -table["Q_E"], -table["A_E"]],
    }


def joined_coefficients(
    rows: dict[str, list[sp.Expr]],
    table: dict[str, sp.Symbol | sp.Expr],
) -> list[sp.Expr]:
    weights = weight_vectors(table)
    beta = [
        left + right for left, right in zip(weights["beta_p"], weights["beta_E"])
    ]
    result = [sp.S.Zero] * 6
    for row_index, name in enumerate(ROW_ORDER):
        result = add(
            result,
            scale(rows[name], weights["alpha"][row_index], size=6),
            scale(lift(rows[name], table["log_a"]), beta[row_index]),
            size=6,
        )
    return result


def decode_complex(value: dict[str, str]) -> complex:
    return complex(float(value["re"]), float(value["im"]))


def dot(coefficients: list[complex], basis: list[complex]) -> complex:
    require(len(coefficients) <= len(basis), "basis too short")
    return sum(coefficient * basis[index] for index, coefficient in enumerate(coefficients))


def compare_complex(actual: complex, encoded: dict[str, str], message: str) -> float:
    error = abs(actual - decode_complex(encoded))
    require(error < TOLERANCE, f"{message}: {error}")
    return error


def check_symbolic(payload: dict, table: dict[str, sp.Symbol | sp.Expr]) -> dict:
    direct = direct_rows(table)
    moments = moment_rows(table)
    centered = centered_rows(table)
    maximum_nonzero = 0
    for name in ROW_ORDER:
        require(rows_equal(direct[name], moments[name]), f"independent dual derivation failed: {name}")
        stored_direct = parse_row(payload["symbolic"]["polynomial_law_rows"][name], table)
        stored_moment = parse_row(payload["symbolic"]["moment_vector_rows"][name], table)
        stored_centered = parse_row(payload["symbolic"]["centered_rows"][name], table)
        stored_tangent = parse_row(payload["symbolic"]["tangent_rows"][name], table)
        require(rows_equal(direct[name], stored_direct), f"stored direct row drifted: {name}")
        require(rows_equal(moments[name], stored_moment), f"stored moment row drifted: {name}")
        require(rows_equal(translate(direct[name], table["log_a"]), centered[name]), f"binomial centered translation failed: {name}")
        require(rows_equal(centered[name], stored_centered), f"stored centered row drifted: {name}")
        require(rows_equal(lift(direct[name], table["log_a"]), stored_tangent), f"stored tangent row drifted: {name}")

    identities = payload["symbolic"]["identities"]
    for group in ("moment_vector_vs_polynomial_law", "centered_translation", "tangent_lift"):
        require(set(identities[group]) == set(ROW_ORDER), f"identity group drifted: {group}")
        require(all(value == "0" for value in identities[group].values()), f"nonzero identity: {group}")
    for name in ("rate_cancellation", "edge_additive_ownership", "lift_zero_endpoint"):
        require(sp.sympify(identities[name], locals=table) == 0, f"nonzero identity: {name}")
    return {"exact_symbolic_rows": 16, "maximum_nonzero": maximum_nonzero}


def check_exact_fixture(payload: dict, table: dict[str, sp.Symbol | sp.Expr]) -> dict:
    substitution = parse_substitution(payload["exact_fixture"]["substitution"], table)
    direct = direct_rows(table)
    centered = centered_rows(table)
    for name in ROW_ORDER:
        expected = [
            sp.simplify(sp.sympify(value).subs(substitution)) for value in direct[name]
        ]
        expected_centered = [
            sp.simplify(sp.sympify(value).subs(substitution))
            for value in centered[name]
        ]
        expected_tangent = [
            sp.simplify(value.subs(substitution))
            for value in lift(direct[name], table["log_a"])
        ]
        require(rows_equal(expected, parse_row(payload["exact_fixture"]["rows"][name], table)), f"fixture row drifted: {name}")
        require(rows_equal(expected_centered, parse_row(payload["exact_fixture"]["centered_rows"][name], table)), f"fixture centered row drifted: {name}")
        require(rows_equal(expected_tangent, parse_row(payload["exact_fixture"]["tangent_rows"][name], table)), f"fixture tangent row drifted: {name}")

    joined = [sp.simplify(value.subs(substitution)) for value in joined_coefficients(direct, table)]
    stored_joined = parse_row(payload["exact_fixture"]["P_lin_coefficients"], table)
    require(rows_equal(joined, stored_joined), "fixture joined polynomial drifted")
    degree = max(index for index, value in enumerate(joined) if value != 0)
    require(degree == payload["exact_fixture"]["P_lin_degree"] == 5, "fixture degree drifted")
    return {"degree": degree, "nonzero_coefficients": sum(value != 0 for value in joined)}


def check_low_n(
    payload: dict,
    low_n_payload: dict,
    table: dict[str, sp.Symbol | sp.Expr],
) -> dict:
    substitution = parse_substitution(payload["exact_fixture"]["substitution"], table)
    generic_rows = direct_rows(table)
    generic_joined = joined_coefficients(generic_rows, table)
    source_by_id = {fixture["id"]: fixture for fixture in low_n_payload["fixtures"]}
    maximum_projection_error = 0.0
    maximum_lift_identity_error = 0.0

    for recorded in payload["low_n_interface"]["fixtures"]:
        source = source_by_id[recorded["id"]]
        local = dict(substitution)
        local[table["L_N"]] = sp.Float(math.log(float(source["parameters"]["N"])), 18)
        local[table["log_a"]] = sp.Float(math.log(float(source["parameters"]["a"])), 18)
        rows = {
            name: [
                complex(sp.N(sp.sympify(value).subs(local), 18))
                for value in generic_rows[name]
            ]
            for name in ROW_ORDER
        }
        tangents = {
            name: [
                complex(sp.N(value.subs(local), 18))
                for value in lift(generic_rows[name], table["log_a"])
            ]
            for name in ROW_ORDER
        }
        joined = [complex(sp.N(value.subs(local), 18)) for value in generic_joined]

        for functional in FUNCTIONALS:
            basis = [decode_complex(value) for value in source["basis"][functional]]
            stored = recorded["projections"][functional]
            for name in ROW_ORDER:
                base_value = dot(rows[name], basis)
                tangent_value = dot(tangents[name], basis)
                maximum_projection_error = max(
                    maximum_projection_error,
                    compare_complex(base_value, stored["base"][name], f"{recorded['id']} {functional} {name} base"),
                    compare_complex(tangent_value, stored["lifted"][name], f"{recorded['id']} {functional} {name} tangent"),
                )
                multiplication = dot([0j, *rows[name]], basis)
                lift_expected = 1j * (
                    multiplication - float(local[table["log_a"]]) * base_value
                )
                maximum_lift_identity_error = max(
                    maximum_lift_identity_error, abs(tangent_value - lift_expected)
                )
            maximum_projection_error = max(
                maximum_projection_error,
                compare_complex(dot(joined, basis), stored["P_lin"], f"{recorded['id']} {functional} P_lin"),
            )

    require(maximum_lift_identity_error < TOLERANCE, "low-N lift identity drifted")
    require(
        abs(maximum_lift_identity_error - payload["low_n_interface"]["maximum_lift_projection_error"])
        < TOLERANCE,
        "stored low-N lift maximum drifted",
    )
    return {
        "fixtures": len(payload["low_n_interface"]["fixtures"]),
        "maximum_projection_replay_error": maximum_projection_error,
        "maximum_lift_identity_error": maximum_lift_identity_error,
    }


def main() -> int:
    set_below_normal_priority()
    require(RESULT_PATH.is_file(), "missing compiler result")
    require(NOTE_PATH.is_file(), "missing compiler note")
    require(LOW_N_PATH.is_file(), "missing low-N basis")
    payload = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    low_n_payload = json.loads(LOW_N_PATH.read_text(encoding="utf-8"))
    require(payload["kind"] == KIND, "kind drifted")
    require(
        payload["status"] == "c1_physical_observation_polynomial_compiler_validated",
        "status drifted",
    )
    require(payload["schema"]["coefficient_order"] == "ascending powers", "coefficient order drifted")
    require(payload["schema"]["row_order"] == list(ROW_ORDER), "row order drifted")
    require("not the later phase angle" in payload["schema"]["symbol_guard"], "chi symbol guard missing")

    for source in payload["source_audit"].values():
        path = REPO_ROOT / source["path"]
        require(path.is_file(), f"missing audited source: {path}")
        require(file_hash(path) == source["sha256"], f"source hash drifted: {path}")

    table = symbol_table()
    symbolic = check_symbolic(payload, table)
    fixture = check_exact_fixture(payload, table)
    low_n = check_low_n(payload, low_n_payload, table)
    rows = payload["rows"]
    require(len(rows) == payload["summary"]["rows"] == 20, "gate row count drifted")
    require(len({row["id"] for row in rows}) == len(rows), "duplicate gate row id")
    require(sum(row["readiness"] == "open" for row in rows) == 1, "open-row count drifted")
    require(payload["summary"]["physical_numeric_rows"] == 0, "physical row boundary drifted")
    require(payload["summary"]["physical_signed_bounds"] == 0, "signed-bound boundary drifted")
    note = NOTE_PATH.read_text(encoding="utf-8")
    require(payload["success"] in note, "note lost success line")
    require("not a proof of RH" in note, "note lost proof boundary")
    print(
        "independent physical-observation compiler check passed: "
        f"{symbolic['exact_symbolic_rows']} exact symbolic rows, "
        f"degree-{fixture['degree']} joined fixture, "
        f"{low_n['fixtures']} low-N fixtures, "
        f"max projection replay error {low_n['maximum_projection_replay_error']:.3e}, "
        f"max lift identity error {low_n['maximum_lift_identity_error']:.3e}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
