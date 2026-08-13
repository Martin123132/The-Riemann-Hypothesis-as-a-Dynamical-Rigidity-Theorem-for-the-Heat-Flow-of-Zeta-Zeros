#!/usr/bin/env python3
"""Independently check the quadratic-residual two-carrier kernel gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_quadratic_residual_two_carrier_kernel_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def complex_symbol(prefix: str) -> sp.Expr:
    real, imag = sp.symbols(f"{prefix}_R {prefix}_I", real=True)
    return real + sp.I * imag


def is_zero(expression: sp.Expr) -> bool:
    return sp.simplify(sp.expand_complex(expression)) == 0


def check_sources(payload: dict) -> int:
    require(set(payload["source_sha256"]) == set(payload["source_audit"]), "source keys")
    require(len(payload["source_audit"]) == 4, "source count")
    for key, item in payload["source_audit"].items():
        path = REPO_ROOT / item["path"]
        require(path.is_file(), f"missing source {key}")
        digest = file_hash(path)
        require(digest == item["sha256"], f"source hash drift {key}")
        require(digest == payload["source_sha256"][key], f"source map drift {key}")
    return len(payload["source_audit"])


def check_symbolics() -> int:
    audits = 0
    c_l = complex_symbol("check_c_l")
    c_m = complex_symbol("check_c_m")
    q_l = complex_symbol("check_q_l")
    q_m = complex_symbol("check_q_m")
    r_l = complex_symbol("check_r_l")
    r_m = complex_symbol("check_r_m")
    rx_l = complex_symbol("check_rx_l")
    rx_m = complex_symbol("check_rx_m")
    a_l, a_m = r_l * c_l, r_m * c_m
    n_l, n_m = r_l * q_l + rx_l * c_l, r_m * q_m + rx_m * c_m

    direct_t = c_l * n_m + n_l * c_m - a_l * q_m - q_l * a_m
    factored_t = (r_m - r_l) * (c_l * q_m - c_m * q_l) + (rx_l + rx_m) * c_l * c_m
    require(is_zero(direct_t - factored_t), "transpose factorization")
    audits += 1

    delta_l = complex_symbol("check_delta_l")
    delta_m = complex_symbol("check_delta_m")
    d_l = complex_symbol("check_d_l")
    d_m = complex_symbol("check_d_m")
    physical_q_l = (r_l + delta_l) * c_l + d_l
    physical_q_m = (r_m + delta_m) * c_m + d_m
    direct_wedge = c_l * physical_q_m - c_m * physical_q_l
    factored_wedge = (
        c_l * c_m * ((r_m - r_l) + (delta_m - delta_l))
        + c_l * d_m
        - c_m * d_l
    )
    require(is_zero(direct_wedge - factored_wedge), "physical wedge")
    audits += 1

    direct_h = (
        c_l * sp.conjugate(n_m)
        + n_l * sp.conjugate(c_m)
        - a_l * sp.conjugate(q_m)
        - q_l * sp.conjugate(a_m)
    )
    factored_h = (
        (sp.conjugate(r_m) - r_l)
        * (c_l * sp.conjugate(q_m) - q_l * sp.conjugate(c_m))
        + (rx_l + sp.conjugate(rx_m)) * c_l * sp.conjugate(c_m)
    )
    require(is_zero(direct_h - factored_h), "Hermitian factorization")
    audits += 1

    x_l, x_m, u_x = sp.symbols("check_x_l check_x_m check_u_x", real=True)
    ideal_r_l, ideal_r_m = sp.I * x_l / 2, sp.I * x_m / 2
    ideal_rx = -sp.I * u_x / 2
    ideal_h = (sp.conjugate(ideal_r_m) - ideal_r_l) ** 2 + ideal_rx + sp.conjugate(ideal_rx)
    ideal_t = (ideal_r_m - ideal_r_l) ** 2 + 2 * ideal_rx
    require(is_zero(ideal_h + (x_l + x_m) ** 2 / 4), "ideal H")
    audits += 1
    require(is_zero(ideal_t + (x_m - x_l) ** 2 / 4 + sp.I * u_x), "ideal T")
    audits += 1

    u_v = complex_symbol("check_u_v")
    u_n = complex_symbol("check_u_n")
    u_a = complex_symbol("check_u_a")
    u_q = complex_symbol("check_u_q")
    direct_primitive = 2 * (sp.re(u_v) * sp.re(u_n) - sp.re(u_a) * sp.re(u_q))
    h_projection = (
        u_v * sp.conjugate(u_n)
        + u_n * sp.conjugate(u_v)
        - u_a * sp.conjugate(u_q)
        - u_q * sp.conjugate(u_a)
    )
    t_projection = u_v * u_n + u_n * u_v - u_a * u_q - u_q * u_a
    require(is_zero(direct_primitive - sp.re(h_projection + t_projection) / 2), "polarization")
    audits += 1

    lambda_l, lambda_m, log_a = sp.symbols("check_lambda_l check_lambda_m check_log_a", real=True)
    lift_l, lift_m = lambda_l - log_a, lambda_m - log_a
    require(
        is_zero(sp.I * lift_l + sp.conjugate(sp.I * lift_m) - sp.I * (lambda_l - lambda_m)),
        "Hermitian multiplier",
    )
    audits += 1
    require(
        is_zero(sp.I * lift_l + sp.I * lift_m - sp.I * (lambda_l + lambda_m - 2 * log_a)),
        "transpose multiplier",
    )
    audits += 1
    return audits


def check_finite_functional(payload: dict) -> int:
    i = sp.I
    xs = (sp.Rational(-2), sp.Rational(1, 3), sp.Rational(5, 2))
    weights = (1 + i / 2, -sp.Rational(2, 3) + i / 5, sp.Rational(3, 7) - 2 * i / 9)
    c = (1 + i / 7, sp.Rational(2, 3) - i / 4, -sp.Rational(1, 5) + 2 * i / 3)
    q = (sp.Rational(4, 5) - i / 6, -sp.Rational(3, 8) + i / 2, sp.Rational(7, 9) + i / 10)
    r = (-sp.Rational(1, 3) + i / 2, sp.Rational(5, 7) - i / 9, -sp.Rational(2, 5) - 3 * i / 8)
    rx = (i / 11, -sp.Rational(1, 13) - i / 17, sp.Rational(2, 19) - i / 23)
    a = tuple(r_j * c_j for r_j, c_j in zip(r, c))
    n = tuple(r_j * q_j + rx_j * c_j for r_j, q_j, rx_j, c_j in zip(r, q, rx, c))

    def apply(values: tuple[sp.Expr, ...], lifted: bool = False) -> sp.Expr:
        factors = (i * x for x in xs) if lifted else (sp.Integer(1) for _ in xs)
        return sum(factor * weight * value for factor, weight, value in zip(factors, weights, values))

    u_c, u_n, u_a, u_q = (apply(values) for values in (c, n, a, q))
    du_c, du_n, du_a, du_q = (apply(values, lifted=True) for values in (c, n, a, q))
    direct_primitive = 2 * (sp.re(u_c) * sp.re(u_n) - sp.re(u_a) * sp.re(u_q))
    direct_current = 2 * (
        sp.re(du_c) * sp.re(u_n)
        + sp.re(u_c) * sp.re(du_n)
        - sp.re(du_a) * sp.re(u_q)
        - sp.re(u_a) * sp.re(du_q)
    )

    kernel_primitive = 0
    kernel_current = 0
    diagonal_h = 0
    diagonal_t = 0
    for left in range(3):
        for right in range(3):
            k_h = (
                c[left] * sp.conjugate(n[right])
                + n[left] * sp.conjugate(c[right])
                - a[left] * sp.conjugate(q[right])
                - q[left] * sp.conjugate(a[right])
            )
            k_t = c[left] * n[right] + n[left] * c[right] - a[left] * q[right] - q[left] * a[right]
            h_weight = weights[left] * sp.conjugate(weights[right])
            t_weight = weights[left] * weights[right]
            kernel_primitive += (k_h * h_weight + k_t * t_weight) / 2
            kernel_current += (
                i * (xs[left] - xs[right]) * k_h * h_weight
                + i * (xs[left] + xs[right]) * k_t * t_weight
            ) / 2
            if left == right:
                diagonal_h += i * (xs[left] - xs[right]) * k_h * h_weight / 2
                diagonal_t += i * (xs[left] + xs[right]) * k_t * t_weight / 2

    primitive_difference = sp.simplify(direct_primitive - sp.re(kernel_primitive))
    current_difference = sp.simplify(direct_current - sp.re(kernel_current))
    diagonal_h = sp.simplify(sp.re(diagonal_h))
    diagonal_t = sp.simplify(sp.re(diagonal_t))
    require(primitive_difference == 0, "finite primitive")
    require(current_difference == 0, "finite current")
    require(diagonal_h == 0, "finite H diagonal")
    require(diagonal_t != 0, "finite T diagonal guard")

    stored = payload["finite_functional"]
    require(stored["carrier_count"] == 3, "stored carrier count")
    require(stored["primitive_difference"] == str(primitive_difference), "stored primitive")
    require(stored["current_difference"] == str(current_difference), "stored current")
    require(stored["hermitian_diagonal_current"] == str(diagonal_h), "stored H diagonal")
    require(stored["transpose_diagonal_current"] == str(diagonal_t), "stored T diagonal")
    return 4


def check_rows(payload: dict) -> None:
    suffixes = (
        "domain",
        "functional",
        "rows",
        "transpose",
        "wedge",
        "hermitian",
        "primitive",
        "current",
        "diagonal",
        "ideal",
        "match",
        "phases",
        "grouping",
        "join",
        "target",
        "boundary",
    )
    expected = [f"qrtk_{index:02d}_{suffix}" for index, suffix in enumerate(suffixes, start=1)]
    require([row["id"] for row in payload["rows"]] == expected, "row ids")
    require(sum(row["readiness"] == "open" for row in payload["rows"]) == 1, "open rows")
    require(payload["rows"][-2]["readiness"] == "open", "frontier row")


def main() -> None:
    require(RESULT_PATH.is_file(), "missing result")
    require(NOTE_PATH.is_file(), "missing note")
    payload = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    require(payload["kind"] == KIND, "kind")
    require(payload["schema_version"] == 1, "schema")
    source_audits = check_sources(payload)
    symbolic_audits = check_symbolics()
    finite_audits = check_finite_functional(payload)
    check_rows(payload)

    summary = payload["summary"]
    require(summary["rows"] == 16, "summary rows")
    require(summary["source_audits"] == 4, "source audits")
    require(summary["symbolic_audits"] == 8, "symbolic audits")
    require(summary["finite_functional_audits"] == 4, "finite audits")
    require(summary["kernel_families"] == 2, "kernel families")
    require(summary["new_phase_families"] == 0, "new phase families")
    require(summary["hermitian_diagonal_current"] == 0, "H diagonal")
    require(summary["transpose_diagonal_guards"] == 1, "T diagonal")
    require(summary["open_rows"] == 1, "open summary")

    note = NOTE_PATH.read_text(encoding="utf-8")
    require(payload["success"] in note, "note success")
    require("K_H" in note and "K_T" in note, "kernel notation")
    require("phase differences" in note and "phase sums" in note, "phase closure")
    require("not a proof" in note.lower(), "boundary note")
    require("prize-level conclusion" in note, "prize boundary")

    print(
        "validated Newman C1 quadratic residual two-carrier kernel gate: "
        f"16 rows, 0 issues, {symbolic_audits} symbolic audits, "
        f"{finite_audits} finite-functional audits, {source_audits} source audits, "
        "2 exact kernel families, 1 open signed theorem row"
    )


if __name__ == "__main__":
    main()
