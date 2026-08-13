#!/usr/bin/env python3
"""Independently check the modular-primitive/contact nonclosure gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_modular_primitive_contact_nonclosure_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def expression_hash(expression: sp.Expr) -> str:
    return hashlib.sha256(sp.srepr(sp.expand(expression)).encode("utf-8")).hexdigest()


def derivative_of(function: sp.Expr, x: sp.Symbol, order: int) -> sp.Expr:
    return function if order == 0 else sp.diff(function, x, order)


def heat_monomial(m: int, z: sp.Symbol, direction: int) -> sp.Expr:
    return sp.expand(
        sum(
            direction**order * sp.diff(z**m, z, 2 * order) / sp.factorial(order)
            for order in range(m // 2 + 1)
        )
    )


def main() -> int:
    require(RESULT_PATH.is_file(), "missing result artifact")
    require(NOTE_PATH.is_file(), "missing theorem note")
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    require(artifact.get("kind") == KIND, "kind drift")

    sources = artifact.get("source_audit", {})
    require(len(sources) == 4, "source count drift")
    for source in sources.values():
        path = REPO_ROOT / source["path"]
        require(path.is_file(), f"missing source {path}")
        require(file_hash(path) == source["sha256"], f"source hash drift {path}")

    expected_ids = [f"mpc_{index:02d}_{suffix}" for index, suffix in enumerate(
        (
            "sources", "theta_coordinate", "field", "operator_jets", "summands",
            "bessel_density", "bessel_transform", "nonvanishing", "model",
            "field_model", "phase", "shape", "heat", "embedding", "primitive",
            "reflection", "probability", "nonclosure", "arithmetic", "uniform",
            "boundary",
        ),
        start=1,
    )]
    require([row["id"] for row in artifact.get("rows", [])] == expected_ids, "row roster drift")

    x, t = sp.symbols("x t", real=True)
    cfun = sp.Function("C")(x)
    operator = -4 * t**2 * sp.diff(cfun, x, 2) + 4 * t * x * sp.diff(cfun, x) + (2 * t - 1 - x**2) * cfun
    stored_operator = artifact.get("operator_jet_audit", [])
    require(len(stored_operator) == 17, "operator audit count drift")
    for j, stored in enumerate(stored_operator):
        direct = sp.expand(sp.diff(operator, x, j))
        expected = (
            -4 * t**2 * derivative_of(cfun, x, j + 2)
            + 4 * t * x * derivative_of(cfun, x, j + 1)
            + (4 * t * j + 2 * t - 1 - x**2) * derivative_of(cfun, x, j)
        )
        if j >= 1:
            expected -= 2 * j * x * derivative_of(cfun, x, j - 1)
        if j >= 2:
            expected -= j * (j - 1) * derivative_of(cfun, x, j - 2)
        require(sp.simplify(direct - expected) == 0, f"independent operator drift j={j}")
        require(expression_hash(direct) == stored["sha256"], f"stored operator hash drift j={j}")

    h, alpha, beta = sp.symbols("h alpha beta", nonzero=True)
    stored_contacts = artifact.get("contact_audit", [])
    require(len(stored_contacts) == 15, "contact audit count drift")
    for m, stored in zip(range(2, 17), stored_contacts, strict=True):
        theta_jet = -sp.Rational(1, 2) + h**m * (1 + alpha * h + beta * h**2)
        heat_jet = (theta_jet + sp.Rational(1, 2)) / 8
        theta_ratio = sp.factor(sp.diff(theta_jet, h, m + 1).subs(h, 0) / ((m + 1) * sp.diff(theta_jet, h, m).subs(h, 0)))
        heat_ratio = sp.factor(sp.diff(heat_jet, h, m + 1).subs(h, 0) / ((m + 1) * sp.diff(heat_jet, h, m).subs(h, 0)))
        require(theta_ratio == alpha and heat_ratio == alpha, f"independent contact drift m={m}")
        require(stored["multiplicity"] == m, f"stored contact multiplicity drift m={m}")

    stored_models = artifact.get("model_audit", [])
    require(len(stored_models) == 15, "model audit count drift")
    phases = (sp.Rational(5, 4) * sp.pi, sp.Rational(3, 2) * sp.pi, sp.Rational(7, 4) * sp.pi)
    cursor = 0
    r = sp.symbols("r", real=True)
    for m in range(2, 7):
        for phase in phases:
            stored = stored_models[cursor]
            cursor += 1
            field = sp.simplify(phase * sp.cot(phase) - m - 1 + r)
            ell = sp.simplify((2 * field - m) / 4)
            require(stored["multiplicity"] == m and sp.sympify(stored["phase"]) == phase, "stored model index drift")
            stored_field = sp.sympify(stored["field"], locals={"r": r})
            stored_ell = sp.sympify(stored["ell"], locals={"r": r})
            require(sp.simplify(stored_field - field) == 0, f"stored field drift m={m}, phase={phase}")
            require(sp.simplify(stored_ell - ell) == 0, f"stored ell drift m={m}, phase={phase}")

    mp.mp.dps = 80
    stored_bessel = artifact.get("bessel_audit", [])
    require(len(stored_bessel) == 3, "Bessel audit count drift")
    for stored in stored_bessel:
        c = mp.mpf(stored["c"])
        a = mp.mpf(stored["a"])

        def transform(value: mp.mpf) -> mp.mpf:
            return mp.re(mp.besselk(1j * value / 4, a)) / mp.besselk(0, a)

        value = transform(c)
        log_derivative = mp.diff(transform, c) / value
        require(value > mp.mpf(31) / 32, f"independent Bessel lower bound failed c={c}")
        require(abs(value - mp.mpf(stored["transform_at_c"])) < mp.mpf("1e-50"), f"Bessel value drift c={c}")
        require(abs(log_derivative - mp.mpf(stored["log_derivative_at_c"])) < mp.mpf("1e-50"), f"Bessel derivative drift c={c}")

    u, v = sp.symbols("u v", real=True)
    kernel = 8 * sp.sinh(v - u)
    stored_kernel = artifact.get("primitive_kernel_audit", {})
    require(sp.simplify(kernel.subs(v, u)) == 0, "kernel diagonal drift")
    require(sp.simplify(sp.diff(kernel, u).subs(v, u)) == -8, "kernel derivative drift")
    require(sp.simplify(sp.diff(kernel, u, 2) - kernel) == 0, "kernel homogeneous equation drift")
    require(expression_hash(kernel) == stored_kernel["sha256"], "kernel hash drift")

    phi = sp.symbols("q0:16")
    jets: list[sp.Expr] = [sp.Symbol("s0"), -sp.Rational(1, 2)]
    for order in range(14):
        jets.append(sp.expand(jets[order] + (0 if order % 2 else 8 * phi[order])))
    stored_modular = artifact.get("modular_jet_audit", [])
    require(len(stored_modular) == 8, "modular jet count drift")
    for stored in stored_modular:
        order = stored["order"]
        require(order % 2 == 1 and jets[order] == -sp.Rational(1, 2), f"modular odd jet drift order={order}")
        require(stored["jet"] == "-1/2", f"stored modular odd jet drift order={order}")

    z = sp.symbols("z")
    stored_hermite = artifact.get("hermite_audit", [])
    require(len(stored_hermite) == 15, "Hermite audit count drift")
    for m, stored in zip(range(2, 17), stored_hermite, strict=True):
        forward = heat_monomial(m, z, -1)
        backward = heat_monomial(m, z, 1)
        require(int(sp.Poly(forward, z).count_roots(-sp.oo, sp.oo)) == m, f"forward root drift m={m}")
        require(int(sp.Poly(backward, z).count_roots(-sp.oo, sp.oo)) == m % 2, f"backward root drift m={m}")
        require(expression_hash(forward) == stored["forward_sha256"], f"forward hash drift m={m}")
        require(expression_hash(backward) == stored["backward_sha256"], f"backward hash drift m={m}")

    counts = artifact.get("counts", {})
    expected_counts = {
        "rows": 21,
        "sources": 4,
        "operator_jets": 17,
        "contact_multiplicities": 15,
        "sinc_bessel_models": 15,
        "bessel_diagnostics": 3,
        "modular_odd_jets": 8,
        "hermite_checks": 15,
        "open_handoffs": 2,
        "actual_xi_field_bounds": 0,
        "degree_uniform_bounds": 0,
        "rh_conclusions": 0,
    }
    require(counts == expected_counts, "count contract drift")

    note = NOTE_PATH.read_text(encoding="utf-8")
    for phrase in (
        "R_*(-u)-R_*(u)=sinh(u)",
        "n^(-1/2)",
        "does not reproduce the discrete Jacobi-theta roster",
        "not a proof of RH",
    ):
        require(phrase in note, f"note boundary phrase missing: {phrase}")
    require("does not reproduce the discrete Jacobi-theta summand roster" in artifact.get("proof_boundary", ""), "artifact proof boundary drift")

    print(
        "validated modular-primitive contact nonclosure gate: "
        "21 rows, 4 sources, 17 operator jets, 15 contact multiplicities, "
        "15 sinc-Bessel models, 3 Bessel diagnostics, 8 modular odd jets, "
        "15 Hermite checks, 2 open handoffs, 0 Xi field bounds, "
        "0 degree-uniform bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
