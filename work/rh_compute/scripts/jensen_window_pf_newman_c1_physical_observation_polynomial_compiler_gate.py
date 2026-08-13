#!/usr/bin/env python3
"""Build the exact C1 physical-observation polynomial compiler gate."""

from __future__ import annotations

import ctypes
from dataclasses import asdict, dataclass
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

SOURCE_PATHS = {
    "moment_rows": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_five_moment_explicit_phi_quadratic_reduction.json"
    ),
    "observation_image": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_observation_image_"
        "compression_gate.json"
    ),
    "centered_amplitude": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "amplitude_gate.json"
    ),
    "edge_ownership": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_edge_affine_ownership_gate.json",
    "low_n_basis": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_fixed_cell_low_n_determinant_evaluator_gate.json",
}

ROW_ORDER = ("P_V", "P_N", "P_A", "P_Q")
BASE_DEGREES = {"P_V": 2, "P_N": 4, "P_A": 3, "P_Q": 3}
LIFTED_DEGREES = {name: degree + 1 for name, degree in BASE_DEGREES.items()}
LOW_N_FUNCTIONALS = ("carrier", "band_blocks", "near", "remote")


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

    require(
        payloads["moment_rows"]["counts"]["real_bulk_observations"] == 4,
        "moment-row source lost four observations",
    )
    require(
        payloads["moment_rows"]["counts"]["complex_base_coefficient_vectors"] == 5,
        "moment-row source lost five coefficient vectors",
    )
    require(
        payloads["observation_image"]["counts"]["maximum_polynomial_degree"] == 5,
        "observation-image degree ceiling drifted",
    )
    require(
        payloads["centered_amplitude"]["counts"]["physical_base_polynomials"] == 4,
        "centered-amplitude source lost four base polynomials",
    )
    require(
        payloads["edge_ownership"]["summary"]["exact_affine_polynomials"] == 1,
        "edge affine ownership drifted",
    )
    require(
        payloads["low_n_basis"]["summary"]["basis_degree"] == 6,
        "low-N basis degree drifted",
    )
    require(
        payloads["low_n_basis"]["summary"]["fixtures"] == 3,
        "low-N fixture count drifted",
    )
    audit = {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }
    return payloads, audit


def symbols() -> dict[str, sp.Symbol | sp.Expr]:
    ell, x = sp.symbols("ell x", real=True)
    log_n, log_a, u_n_x = sp.symbols("L_N log_a u_N_x", real=True)
    c, b, c_x, b_x = sp.symbols("c b c_x b_x", real=True)
    rho_1, rho_2, rho_1_x, rho_2_x, chi_rate = sp.symbols(
        "rho_1 rho_2 rho_1_x rho_2_x chi_rate"
    )
    return {
        "ell": ell,
        "x": x,
        "L_N": log_n,
        "log_a": log_a,
        "u_N": log_a - log_n,
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
        "s_prime": c + sp.I * b,
        "s_second": c_x + sp.I * b_x,
    }


def coefficient_vector(expression: sp.Expr, variable: sp.Symbol, size: int) -> list[sp.Expr]:
    polynomial = sp.Poly(sp.expand(expression), variable)
    degree = polynomial.degree()
    require(degree is sp.S.NegativeInfinity or degree < size, "coefficient vector overflow")
    return [sp.expand(polynomial.nth(index)) for index in range(size)]


def vector_rows(values: dict[str, sp.Symbol | sp.Expr]) -> dict[str, list[sp.Expr]]:
    log_n = values["L_N"]
    rho_1 = values["rho_1"]
    rho_2 = values["rho_2"]
    rho_1_x = values["rho_1_x"]
    rho_2_x = values["rho_2_x"]
    s_prime = values["s_prime"]
    s_second = values["s_second"]
    chi_rate = values["chi_rate"]
    b = values["b"]
    b_x = values["b_x"]
    u_n = values["u_N"]
    u_n_x = values["u_N_x"]

    v_0 = sp.Matrix([1, rho_1, rho_2, 0, 0])
    v_1 = sp.Matrix(
        [
            log_n,
            log_n * rho_1 - 1,
            log_n * rho_2 - rho_1,
            -rho_2,
            0,
        ]
    )
    v_2 = sp.Matrix(
        [
            log_n**2,
            log_n**2 * rho_1 - 2 * log_n,
            log_n**2 * rho_2 - 2 * log_n * rho_1 + 1,
            rho_1 - 2 * log_n * rho_2,
            rho_2,
        ]
    )
    e_0 = sp.Matrix([0, rho_1_x, rho_2_x, 0, 0])
    e_1 = sp.Matrix(
        [
            0,
            log_n * rho_1_x,
            log_n * rho_2_x - rho_1_x,
            -rho_2_x,
            0,
        ]
    )
    q_0 = chi_rate * v_0 + s_prime * v_1 + e_0
    q_1 = chi_rate * v_1 + s_prime * v_2 + e_1
    a_0 = s_prime * v_1 + sp.I * b * u_n * v_0
    n_0 = (
        s_second * v_1
        + s_prime * q_1
        + sp.I * (b_x * u_n + b * u_n_x) * v_0
        + sp.I * b * u_n * q_0
    )
    return {
        "P_V": [sp.expand(entry) for entry in v_0],
        "P_N": [sp.expand(entry) for entry in n_0],
        "P_A": [sp.expand(entry) for entry in a_0],
        "P_Q": [sp.expand(entry) for entry in q_0],
    }


def polynomial_laws(
    values: dict[str, sp.Symbol | sp.Expr], variable: sp.Symbol
) -> tuple[dict[str, sp.Expr], dict[str, sp.Expr]]:
    rho_1 = values["rho_1"]
    rho_2 = values["rho_2"]
    rho_1_x = values["rho_1_x"]
    rho_2_x = values["rho_2_x"]
    log_n = values["L_N"]
    s_prime = values["s_prime"]
    s_second = values["s_second"]
    chi_rate = values["chi_rate"]
    b = values["b"]
    b_x = values["b_x"]
    u_n = values["u_N"]
    u_n_x = values["u_N_x"]

    correction = 1 + rho_1 * variable + rho_2 * variable**2
    correction_x = rho_1_x * variable + rho_2_x * variable**2
    rate = s_prime * (log_n - variable) + sp.I * b * u_n
    rate_x = s_second * (log_n - variable) + sp.I * (b_x * u_n + b * u_n_x)
    q_polynomial = (chi_rate + s_prime * (log_n - variable)) * correction + correction_x
    a_polynomial = rate * correction
    n_polynomial = rate * q_polynomial + rate_x * correction
    laws = {
        "C": sp.expand(correction),
        "D": sp.expand(correction_x),
        "R": sp.expand(rate),
        "R_x": sp.expand(rate_x),
        "delta": sp.expand(chi_rate - sp.I * b * u_n),
    }
    rows = {
        "P_V": sp.expand(correction),
        "P_N": sp.expand(n_polynomial),
        "P_A": sp.expand(a_polynomial),
        "P_Q": sp.expand(q_polynomial),
    }
    return laws, rows


def centered_laws(
    values: dict[str, sp.Symbol | sp.Expr], variable: sp.Symbol
) -> tuple[dict[str, sp.Expr], dict[str, sp.Expr]]:
    log_a = values["log_a"]
    u_n = values["u_N"]
    rho_1 = values["rho_1"]
    rho_2 = values["rho_2"]
    rho_1_x = values["rho_1_x"]
    rho_2_x = values["rho_2_x"]
    c = values["c"]
    b = values["b"]
    c_x = values["c_x"]
    s_prime = values["s_prime"]
    s_second = values["s_second"]
    chi_rate = values["chi_rate"]
    u_n_x = values["u_N_x"]

    correction = 1 + rho_1 * (log_a + variable) + rho_2 * (log_a + variable) ** 2
    correction_x = rho_1_x * (log_a + variable) + rho_2_x * (log_a + variable) ** 2
    rate = -s_prime * variable - c * u_n
    rate_x = -s_second * variable - c_x * u_n + sp.I * b * u_n_x
    delta = chi_rate - sp.I * b * u_n
    q_polynomial = (rate + delta) * correction + correction_x
    rows = {
        "P_V": sp.expand(correction),
        "P_N": sp.expand(rate * q_polynomial + rate_x * correction),
        "P_A": sp.expand(rate * correction),
        "P_Q": sp.expand(q_polynomial),
    }
    laws = {
        "C": sp.expand(correction),
        "D": sp.expand(correction_x),
        "R": sp.expand(rate),
        "R_x": sp.expand(rate_x),
        "delta": sp.expand(delta),
    }
    return laws, rows


def lift_coefficients(coefficients: list[sp.Expr], log_a: sp.Expr) -> list[sp.Expr]:
    result = [sp.S.Zero] * (len(coefficients) + 1)
    for index, coefficient in enumerate(coefficients):
        result[index] = sp.expand(result[index] - sp.I * log_a * coefficient)
        result[index + 1] = sp.expand(result[index + 1] + sp.I * coefficient)
    return result


def expressions_equal(left: list[sp.Expr], right: list[sp.Expr]) -> bool:
    return len(left) == len(right) and all(
        sp.expand(a_value - b_value) == 0 for a_value, b_value in zip(left, right)
    )


def expression_list(values: list[sp.Expr]) -> list[str]:
    return [str(sp.expand(value)) for value in values]


def fixture_substitution(values: dict[str, sp.Symbol | sp.Expr]) -> dict[sp.Expr, sp.Expr]:
    return {
        values["L_N"]: sp.Rational(7, 3),
        values["log_a"]: sp.Rational(19, 7),
        values["u_N_x"]: sp.Rational(2, 47),
        values["c"]: -sp.Rational(1, 2),
        values["b"]: sp.Rational(3, 5),
        values["c_x"]: sp.Rational(1, 31),
        values["b_x"]: -sp.Rational(2, 37),
        values["rho_1"]: sp.Rational(2, 5) + sp.I * sp.Rational(1, 7),
        values["rho_2"]: -sp.Rational(1, 11) + sp.I * sp.Rational(2, 13),
        values["rho_1_x"]: sp.Rational(3, 17) - sp.I * sp.Rational(1, 19),
        values["rho_2_x"]: -sp.Rational(2, 23) + sp.I * sp.Rational(1, 29),
        values["chi_rate"]: sp.Rational(1, 41) - sp.I * sp.Rational(3, 43),
    }


def exact_substituted(values: list[sp.Expr], substitution: dict[sp.Expr, sp.Expr]) -> list[sp.Expr]:
    return [sp.simplify(value.subs(substitution)) for value in values]


def complex_payload(value: complex) -> dict[str, str]:
    return {"re": f"{value.real:.17e}", "im": f"{value.imag:.17e}"}


def decode_complex(value: dict[str, str]) -> complex:
    return complex(float(value["re"]), float(value["im"]))


def dot(coefficients: list[complex], basis: list[complex]) -> complex:
    require(len(coefficients) <= len(basis), "basis too short for coefficient vector")
    return sum(coefficient * basis[index] for index, coefficient in enumerate(coefficients))


def shifted_coefficients(coefficients: list[complex]) -> list[complex]:
    return [0j, *coefficients]


def numeric_row_coefficients(
    rows: dict[str, list[sp.Expr]], substitution: dict[sp.Expr, sp.Expr]
) -> dict[str, list[complex]]:
    return {
        name: [complex(sp.N(value.subs(substitution), 18)) for value in coefficients]
        for name, coefficients in rows.items()
    }


def physical_weight_symbols() -> dict[str, sp.Symbol]:
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
    return {name: sp.symbols(name, real=True) for name in names}


def physical_weight_vectors(weights: dict[str, sp.Symbol]) -> dict[str, list[sp.Expr]]:
    return {
        "alpha": [
            weights["N_p_xi"],
            weights["V_p_xi"],
            -weights["Q_p_xi"],
            -weights["A_p_xi"],
        ],
        "beta_p": [
            weights["N_p"],
            weights["V_p"],
            -weights["Q_p"],
            -weights["A_p"],
        ],
        "beta_E": [
            weights["N_E"],
            weights["V_E"],
            -weights["Q_E"],
            -weights["A_E"],
        ],
    }


def fixture_weight_substitution(weights: dict[str, sp.Symbol]) -> dict[sp.Expr, sp.Expr]:
    numerators = (2, -3, 5, -7, 11, -13, 17, -19, 23, -29, 31, -37)
    denominators = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)
    return {
        weights[name]: sp.Rational(numerator, denominator)
        for name, numerator, denominator in zip(weights, numerators, denominators)
    }


def linear_combination(expressions: list[sp.Expr], weights: list[sp.Expr]) -> sp.Expr:
    return sp.expand(sum(weight * expression for weight, expression in zip(weights, expressions)))


def low_n_interface(
    low_n_payload: dict,
    values: dict[str, sp.Symbol | sp.Expr],
    polynomial_rows: dict[str, list[sp.Expr]],
    weights: dict[str, sp.Symbol],
    weight_vectors: dict[str, list[sp.Expr]],
) -> dict:
    base_substitution = fixture_substitution(values)
    weight_substitution = fixture_weight_substitution(weights)
    alpha = [complex(sp.N(value.subs(weight_substitution), 18)) for value in weight_vectors["alpha"]]
    beta = [
        complex(sp.N((left + right).subs(weight_substitution), 18))
        for left, right in zip(weight_vectors["beta_p"], weight_vectors["beta_E"])
    ]
    fixture_rows = []
    maximum_lift_error = 0.0
    maximum_combination_error = 0.0

    for fixture in low_n_payload["fixtures"]:
        parameters = fixture["parameters"]
        substitution = dict(base_substitution)
        substitution[values["L_N"]] = sp.Float(math.log(float(parameters["N"])), 18)
        substitution[values["log_a"]] = sp.Float(math.log(float(parameters["a"])), 18)
        base_rows = numeric_row_coefficients(polynomial_rows, substitution)
        lifted_rows = {
            name: [
                complex(sp.N(value.subs(substitution), 18))
                for value in lift_coefficients(coefficients, values["log_a"])
            ]
            for name, coefficients in polynomial_rows.items()
        }

        projections: dict[str, dict] = {}
        lift_errors = []
        combination_errors = []
        for functional in LOW_N_FUNCTIONALS:
            basis = [decode_complex(value) for value in fixture["basis"][functional]]
            base_projection = {
                name: dot(base_rows[name], basis) for name in ROW_ORDER
            }
            lifted_projection = {
                name: dot(lifted_rows[name], basis) for name in ROW_ORDER
            }
            for name in ROW_ORDER:
                multiplication_projection = dot(shifted_coefficients(base_rows[name]), basis)
                expected = 1j * (
                    multiplication_projection
                    - float(substitution[values["log_a"]]) * base_projection[name]
                )
                lift_errors.append(abs(lifted_projection[name] - expected))

            combined_coefficients = [0j] * 6
            for row_index, name in enumerate(ROW_ORDER):
                for degree, coefficient in enumerate(base_rows[name]):
                    combined_coefficients[degree] += alpha[row_index] * coefficient
                for degree, coefficient in enumerate(lifted_rows[name]):
                    combined_coefficients[degree] += beta[row_index] * coefficient
            combined_direct = dot(combined_coefficients, basis)
            combined_rows = sum(
                alpha[index] * base_projection[name]
                + beta[index] * lifted_projection[name]
                for index, name in enumerate(ROW_ORDER)
            )
            combination_errors.append(abs(combined_direct - combined_rows))
            projections[functional] = {
                "base": {name: complex_payload(base_projection[name]) for name in ROW_ORDER},
                "lifted": {
                    name: complex_payload(lifted_projection[name]) for name in ROW_ORDER
                },
                "P_lin": complex_payload(combined_direct),
            }

        fixture_lift_error = max(lift_errors)
        fixture_combination_error = max(combination_errors)
        maximum_lift_error = max(maximum_lift_error, fixture_lift_error)
        maximum_combination_error = max(maximum_combination_error, fixture_combination_error)
        fixture_rows.append(
            {
                "id": fixture["id"],
                "parameters": {
                    "N": int(parameters["N"]),
                    "a": float(parameters["a"]),
                    "L_N": f"{float(substitution[values['L_N']]):.17e}",
                    "log_a": f"{float(substitution[values['log_a']]):.17e}",
                    "u_N": f"{float(substitution[values['log_a']] - substitution[values['L_N']]):.17e}",
                },
                "projections": projections,
                "maximum_lift_projection_error": fixture_lift_error,
                "maximum_P_lin_combination_error": fixture_combination_error,
            }
        )

    require(maximum_lift_error < 2.0e-12, "low-N lift projection mismatch")
    require(maximum_combination_error < 2.0e-12, "low-N P_lin combination mismatch")
    return {
        "coefficient_fixture": "synthetic exact rates with each low-N fixture's own L_N=log N and log_a=log a",
        "functional_order": list(LOW_N_FUNCTIONALS),
        "fixtures": fixture_rows,
        "maximum_lift_projection_error": maximum_lift_error,
        "maximum_P_lin_combination_error": maximum_combination_error,
    }


def build_certificate(low_n_payload: dict) -> dict:
    values = symbols()
    ell = values["ell"]
    x = values["x"]
    laws, polynomial_expressions = polynomial_laws(values, ell)
    polynomial_rows = {
        name: coefficient_vector(polynomial_expressions[name], ell, 5)
        for name in ROW_ORDER
    }
    moment_rows = vector_rows(values)
    for name in ROW_ORDER:
        require(
            expressions_equal(polynomial_rows[name], moment_rows[name]),
            f"moment-vector/polynomial-law mismatch: {name}",
        )

    lifted_rows = {
        name: lift_coefficients(polynomial_rows[name], values["log_a"])
        for name in ROW_ORDER
    }
    for name in ROW_ORDER:
        direct = coefficient_vector(
            sp.I * (ell - values["log_a"]) * polynomial_expressions[name], ell, 6
        )
        require(expressions_equal(lifted_rows[name], direct), f"lift mismatch: {name}")

    centered_rate_laws, centered_expressions = centered_laws(values, x)
    centered_rows = {
        name: coefficient_vector(centered_expressions[name], x, 5)
        for name in ROW_ORDER
    }
    for name in ROW_ORDER:
        shifted = coefficient_vector(
            polynomial_expressions[name].subs(ell, values["log_a"] + x), x, 5
        )
        require(expressions_equal(centered_rows[name], shifted), f"centered mismatch: {name}")

    rate_cancellation = sp.expand(
        polynomial_expressions["P_A"]
        - polynomial_expressions["P_Q"]
        - ((sp.I * values["b"] * values["u_N"] - values["chi_rate"]) * laws["C"] - laws["D"])
    )
    require(rate_cancellation == 0, "rate cancellation failed")

    weights = physical_weight_symbols()
    weight_vectors = physical_weight_vectors(weights)
    row_expression_list = [polynomial_expressions[name] for name in ROW_ORDER]
    f_alpha = linear_combination(row_expression_list, weight_vectors["alpha"])
    f_beta_p = linear_combination(row_expression_list, weight_vectors["beta_p"])
    f_beta_e = linear_combination(row_expression_list, weight_vectors["beta_E"])
    p_bulk = sp.expand(f_alpha + sp.I * (ell - values["log_a"]) * f_beta_p)
    p_edge = sp.expand(sp.I * (ell - values["log_a"]) * f_beta_e)
    p_total = sp.expand(p_bulk + p_edge)
    p_joined = sp.expand(
        f_alpha
        + sp.I
        * (ell - values["log_a"])
        * linear_combination(
            row_expression_list,
            [left + right for left, right in zip(weight_vectors["beta_p"], weight_vectors["beta_E"])],
        )
    )
    require(sp.expand(p_total - p_joined) == 0, "edge ownership failed")
    require(sp.degree(p_total, ell) <= 5, "joined polynomial exceeded degree five")
    require(sp.expand(p_total.subs(ell, values["log_a"]) - f_alpha.subs(ell, values["log_a"])) == 0, "lift-zero endpoint failed")

    substitution = fixture_substitution(values)
    weight_substitution = fixture_weight_substitution(weights)
    complete_substitution = {**substitution, **weight_substitution}
    fixture_rows = {
        name: expression_list(exact_substituted(polynomial_rows[name], substitution))
        for name in ROW_ORDER
    }
    fixture_lifted = {
        name: expression_list(exact_substituted(lifted_rows[name], substitution))
        for name in ROW_ORDER
    }
    fixture_centered = {
        name: expression_list(exact_substituted(centered_rows[name], substitution))
        for name in ROW_ORDER
    }
    fixture_p_lin = coefficient_vector(p_total.subs(complete_substitution), ell, 6)
    require(sp.degree(sp.Poly(sum(value * ell**index for index, value in enumerate(fixture_p_lin)), ell)) == 5, "fixture did not exercise degree five")

    low_n = low_n_interface(
        low_n_payload, values, polynomial_rows, weights, weight_vectors
    )

    identity_zeros = {
        "moment_vector_vs_polynomial_law": {name: "0" for name in ROW_ORDER},
        "centered_translation": {name: "0" for name in ROW_ORDER},
        "tangent_lift": {name: "0" for name in ROW_ORDER},
        "rate_cancellation": str(rate_cancellation),
        "edge_additive_ownership": str(sp.expand(p_total - p_joined)),
        "lift_zero_endpoint": str(
            sp.expand(
                p_total.subs(ell, values["log_a"])
                - f_alpha.subs(ell, values["log_a"])
            )
        ),
    }
    return {
        "schema": {
            "polynomial_variable": "ell=log n",
            "centered_variable": "x=ell-log a",
            "coefficient_order": "ascending powers",
            "observation_order": ["V", "mathcal N", "A", "Q"],
            "row_order": list(ROW_ORDER),
            "inputs": [
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
            "derived": [
                "u_N=log_a-L_N",
                "s_prime=c+i*b",
                "s_second=c_x+i*b_x",
                "delta=chi_rate-i*b*u_N",
            ],
            "symbol_guard": "chi_rate is the complex rate chi_N of Sections 11.157-11.173, not the later phase angle also denoted chi_N",
            "base_degree_ceilings": BASE_DEGREES,
            "lifted_degree_ceilings": LIFTED_DEGREES,
        },
        "laws": {
            "ell_chart": {name: str(expression) for name, expression in laws.items()},
            "centered_chart": {
                name: str(expression) for name, expression in centered_rate_laws.items()
            },
        },
        "symbolic": {
            "moment_vector_rows": {
                name: expression_list(moment_rows[name]) for name in ROW_ORDER
            },
            "polynomial_law_rows": {
                name: expression_list(polynomial_rows[name]) for name in ROW_ORDER
            },
            "centered_rows": {
                name: expression_list(centered_rows[name]) for name in ROW_ORDER
            },
            "tangent_rows": {
                name: expression_list(lifted_rows[name]) for name in ROW_ORDER
            },
            "identities": identity_zeros,
        },
        "exact_fixture": {
            "substitution": {
                str(key): str(value) for key, value in complete_substitution.items()
            },
            "rows": fixture_rows,
            "centered_rows": fixture_centered,
            "tangent_rows": fixture_lifted,
            "P_lin_coefficients": expression_list(fixture_p_lin),
            "P_lin_degree": int(sp.degree(p_total.subs(complete_substitution), ell)),
            "alpha": expression_list(
                [value.subs(weight_substitution) for value in weight_vectors["alpha"]]
            ),
            "beta_p": expression_list(
                [value.subs(weight_substitution) for value in weight_vectors["beta_p"]]
            ),
            "beta_E": expression_list(
                [value.subs(weight_substitution) for value in weight_vectors["beta_E"]]
            ),
        },
        "low_n_interface": low_n,
    }


def exact_contract() -> dict[str, str]:
    return {
        "compiler": "The compiler maps the twelve explicit scalar inputs into ascending coefficient arrays for P_V=C, P_N=mathcal N, P_A=RC, and P_Q=Q. It derives u_N=log_a-L_N internally, preventing an inconsistent duplicate input.",
        "dual_derivation": "Each base row is constructed both from the five moment vectors v_0,v_1,v_2,e_0,e_1 of Section 11.159 and from the carrier laws C,D,R,R_x,Q,mathcal N of Section 11.160. All twenty symbolic coefficient differences vanish.",
        "centered": "Translation ell=log_a+x exactly reproduces R=-s_prime*x-c*u_N, R_x=-s_second*x-c_x*u_N+i*b*u_N_x, Q=(R+delta)C+D, and mathcal N=RQ+R_xC with delta=chi_rate-i*b*u_N.",
        "lift": "The tangent compiler applies dot(P)(ell)=i(ell-log_a)P(ell), giving the exact ascending recurrence dot(p)_j=i[p_(j-1)-log_a*p_j].",
        "ownership": "For alpha=(N_(p,xi),V_(p,xi),-Q_(p,xi),-A_(p,xi)), beta_p=(N_p,V_p,-Q_p,-A_p), and beta_E=(N_E,V_E,-Q_E,-A_E), P_lin=F_alpha+i(ell-log_a)F_(beta_p+beta_E) splits exactly into one bulk term plus the single edge lift. The edge is not inserted into alpha or counted twice.",
        "low_n": "The exact fixture is recompiled at each validated low-N chart using that chart's L_N=log N and log_a=log a. Base, tangent, and joined P_lin arrays project through the stored degree-six carrier, finite-band, near, and paired-remote basis with the recorded machine-roundoff identities.",
        "missing_inputs": "A physical numerical row now requires explicit certified values of chi_rate, rho_1,rho_2 and their x-derivatives, c,b,c_x,b_x,u_(N,x), plus the genuine endpoint edge observations and retained carrier observations used in alpha and beta. The compiler deliberately does not invent them.",
        "pi_provenance": "No new pi is introduced by the coefficient algebra. Any pi inside u_(N,x), calibrated rates, or the low-N Fourier basis is inherited from the completed-zeta/Fourier normalizations already audited by the source gates.",
        "proof_boundary": "This gate proves exact symbolic compatibility, coefficient ownership, degree ceilings, and a low-N software interface. Its low-N rates are synthetic regression data. It does not evaluate the physical chi_rate or endpoint jets, enumerate the physical carrier, prove interval enclosures, a determinant sign, a flow inequality, Lambda<=0, RH, or a prize-level conclusion.",
    }


def build_rows(contract: dict[str, str], certificate: dict) -> list[GateRow]:
    low_n = certificate["low_n_interface"]
    return [
        GateRow("poc_01_schema", "compiler schema", "validated", "The physical observation compiler has an explicit variable, order, and input contract.", contract["compiler"], "Complex rates remain explicit inputs."),
        GateRow("poc_02_symbol", "symbol collision guard", "guard_validated", "The complex chi rate is not confused with the later phase angle.", certificate["schema"]["symbol_guard"], "Notation is disambiguated in code and output."),
        GateRow("poc_03_moment", "moment-vector derivation", "validated", "The Section 11.159 vectors compile all four observation rows.", "Rows (v_0,n_0,a_0,q_0) are serialized in ascending ell degree.", "No real-part projection is applied at coefficient level."),
        GateRow("poc_04_carrier", "carrier-law derivation", "validated", "The Section 11.160 carrier laws independently compile the same rows.", contract["dual_derivation"], "The equality is coefficient-by-coefficient."),
        GateRow("poc_05_value", "value row", "validated", "P_V=C has degree at most two.", "P_V coefficients agree in both derivations.", "No numerical correction bound is inferred."),
        GateRow("poc_06_scalar", "centered-scalar row", "validated", "P_A=RC has degree at most three.", "P_A coefficients agree in both derivations.", "This is the bulk scalar, not the total endpoint scalar."),
        GateRow("poc_07_value_rate", "value-rate row", "validated", "P_Q=Q has degree at most three.", "P_Q coefficients agree in both derivations.", "chi_rate remains unspecialized."),
        GateRow("poc_08_scalar_rate", "scalar-rate row", "validated", "P_N=RQ+R_xC has degree at most four.", "P_N coefficients agree in both derivations.", "No sign follows from degree four."),
        GateRow("poc_09_centered", "centered chart", "validated", "Large log N and log a terms cancel before norms are taken.", contract["centered"], "The identity uses u_N=log_a-L_N exactly."),
        GateRow("poc_10_defect", "rate cancellation", "validated", "P_A-P_Q=(i*b*u_N-chi_rate)C-D.", "The expanded symbolic difference is zero.", "This cancellation is algebraic, not a bound on delta."),
        GateRow("poc_11_lift", "tangent lift", "validated", "All four auxiliary frequency derivatives use one exact coefficient lift.", contract["lift"], "The source polynomial is held fixed on the auxiliary chart."),
        GateRow("poc_12_degree", "degree ceiling", "validated", "Base degrees (2,4,3,3) lift to (3,5,4,4).", "No moment above ell^5 enters P_lin; the evaluator retains degree six as margin.", "Degree control alone gives no amplitude bound."),
        GateRow("poc_13_alpha", "alpha channel", "validated", "The non-lifted physical channel compiles as F_alpha.", "alpha=(N_(p,xi),V_(p,xi),-Q_(p,xi),-A_(p,xi)).", "The observations still require physical evaluation."),
        GateRow("poc_14_beta", "beta channel", "validated", "The lifted physical channel compiles as i(ell-log_a)F_beta.", "beta=beta_p+beta_E is assembled before projection.", "Componentwise absolute values are not introduced."),
        GateRow("poc_15_edge", "edge ownership", "validated", "The genuine endpoint edge appears exactly once in beta_E.", contract["ownership"], "No edge value is inserted into alpha."),
        GateRow("poc_16_endpoint", "lift-zero endpoint", "validated", "At ell=log a the entire beta lift vanishes.", "P_lin(log a)=F_alpha(log a) exactly.", "The full polynomial need not contain a universal ell-log a factor."),
        GateRow("poc_17_fixture", "exact fixture", "validated", "A rational-complex fixture exercises every input and reaches degree five.", "All base, centered, tangent, and joined coefficients are serialized exactly.", "The fixture is not asserted to be Xi data."),
        GateRow("poc_18_low_n", "low-N interface", "validated", "The compiler output passes through the existing degree-six functional API.", f"3 fixtures; max lift error={low_n['maximum_lift_projection_error']:.3e}; max joined error={low_n['maximum_P_lin_combination_error']:.3e}.", "The basis values are floating-point regression data."),
        GateRow("poc_19_inputs", "physical scalar adapter", "open", "Supply certified calibrated-root values for every compiler input.", contract["missing_inputs"], "No missing rate is fitted or guessed."),
        GateRow("poc_20_boundary", "proof boundary", "guard_validated", "The compiler is infrastructure, not a signed theorem.", contract["proof_boundary"], "All prize-level implications remain open."),
    ]


def render_note(payload: dict) -> str:
    contract = payload["exact"]
    summary = payload["summary"]
    return "\n".join(
        [
            "# Newman C1 Physical Observation Polynomial Compiler Gate",
            "",
            "Date: 2026-08-05",
            "",
            "Status: exact coefficient compiler and low-N interface validated; physical scalar calibration and signed determinant remain open; not a proof of RH.",
            "",
            "## Compiler Contract",
            "",
            contract["compiler"],
            "",
            contract["dual_derivation"],
            "",
            contract["centered"],
            "",
            contract["lift"],
            "",
            "## Physical Ownership",
            "",
            contract["ownership"],
            "",
            "## Validation",
            "",
            f"- Symbolic row identities: `{summary['symbolic_row_identities']}`.",
            f"- Centered translation identities: `{summary['centered_identities']}`.",
            f"- Tangent identities: `{summary['tangent_identities']}`.",
            f"- Exact fixture degree: `{summary['exact_fixture_P_lin_degree']}`.",
            f"- Low-N fixtures: `{summary['low_n_fixtures']}`.",
            f"- Largest low-N lift projection discrepancy: `{summary['maximum_low_n_lift_error']:.3e}`.",
            f"- Largest low-N joined projection discrepancy: `{summary['maximum_low_n_joined_error']:.3e}`.",
            "",
            "## Next Input Layer",
            "",
            contract["missing_inputs"],
            "",
            "The compiler calls the complex coefficient rate `chi_rate`. This is intentionally distinct from the later phase angle that reuses the printed symbol `chi_N` elsewhere in the corpus.",
            "",
            "## Pi Provenance",
            "",
            contract["pi_provenance"],
            "",
            "## Proof Boundary",
            "",
            contract["proof_boundary"],
            "",
            payload["success"],
            "",
        ]
    )


def main() -> int:
    set_below_normal_priority()
    payloads, source_audit = load_and_audit_sources()
    certificate = build_certificate(payloads["low_n_basis"])
    contract = exact_contract()
    rows = build_rows(contract, certificate)
    low_n = certificate["low_n_interface"]
    summary = {
        "rows": len(rows),
        "base_polynomials": 4,
        "symbolic_row_identities": 4,
        "centered_identities": 4,
        "tangent_identities": 4,
        "edge_ownership_identities": 1,
        "exact_fixture_P_lin_degree": certificate["exact_fixture"]["P_lin_degree"],
        "low_n_fixtures": len(low_n["fixtures"]),
        "low_n_functionals_per_fixture": len(LOW_N_FUNCTIONALS),
        "maximum_low_n_lift_error": low_n["maximum_lift_projection_error"],
        "maximum_low_n_joined_error": low_n["maximum_P_lin_combination_error"],
        "physical_numeric_rows": 0,
        "physical_signed_bounds": 0,
    }
    success = (
        "built C1 physical-observation polynomial compiler gate: "
        "4 base rows from 2 exact derivations, 4 centered translations, "
        "4 tangent lifts, single edge ownership, 3 low-N interface fixtures, "
        "0 physical numerical rows and 0 signed bounds"
    )
    payload = {
        "kind": KIND,
        "date": "2026-08-05",
        "status": "c1_physical_observation_polynomial_compiler_validated",
        "source_audit": source_audit,
        "exact": contract,
        **certificate,
        "rows": [asdict(row) for row in rows],
        "summary": summary,
        "proof_boundary": contract["proof_boundary"],
        "success": success,
    }
    atomic_write(RESULT_PATH, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(payload))
    print(success)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
