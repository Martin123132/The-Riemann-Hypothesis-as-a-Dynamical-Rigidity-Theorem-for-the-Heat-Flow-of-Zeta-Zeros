#!/usr/bin/env python3
"""Independently check the arbitrary-multiplicity Jensen boundary layer."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_arbitrary_multiplicity_jensen_boundary_layer_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"
MIN_MULTIPLICITY = 2
MAX_AUDIT_MULTIPLICITY = 16

EXPECTED_COUNTS = {
    "gate_rows": 15,
    "source_artifacts": 7,
    "multiplicity_rows": 15,
    "universal_coefficient_checks": 79,
    "hermite_good_side_checks": 15,
    "imaginary_bad_side_checks": 15,
    "orientation_checks": 15,
    "exact_polynomial_layer_checks": 14,
    "fixed_jensen_shifts": 1,
    "uniform_in_degree_bounds": 0,
    "rh_conclusions": 0,
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rebuild_polynomial(m: int, eta: sp.Symbol, a: sp.Symbol) -> sp.Expr:
    terms = []
    for q in range(m // 2 + 1):
        derivative = sp.diff(eta**m, eta, 2 * q)
        terms.append(a**q * derivative / sp.factorial(q))
    return sp.expand(sum(terms))


def check_source_audit(artifact: dict) -> None:
    audit = artifact.get("source_audit", {})
    require(len(audit) == EXPECTED_COUNTS["source_artifacts"], "source count drifted")
    for key, row in audit.items():
        path = REPO_ROOT / row.get("path", "")
        require(path.is_file(), f"missing source path: {key}")
        require(file_hash(path) == row.get("sha256"), f"source hash drift: {key}")


def check_rows(artifact: dict) -> None:
    rows = artifact.get("rows", [])
    expected_ids = [f"amb_{index:02d}_{suffix}" for index, suffix in enumerate(
        [
            "sources", "fixed_shift", "jensen_operator", "jensen_limit", "heat_limit",
            "combined", "coefficients", "good_side", "bad_side", "collision",
            "boundary", "degree361", "hankel", "route", "proof_boundary",
        ],
        1,
    )]
    require([row.get("id") for row in rows] == expected_ids, "gate row ids drifted")
    require(sum(row.get("readiness") == "open" for row in rows) == 1, "open row drifted")
    require(rows[13].get("id") == "amb_14_route", "route row drifted")
    for row in rows:
        for field in ("role", "readiness", "claim", "certificate", "proof_boundary"):
            require(bool(row.get(field)), f"empty {field}: {row.get('id')}")


def check_multiplicities(artifact: dict) -> None:
    eta, a = sp.symbols("eta a")
    rows = artifact.get("multiplicity_audit", [])
    require(
        [row.get("multiplicity") for row in rows]
        == list(range(MIN_MULTIPLICITY, MAX_AUDIT_MULTIPLICITY + 1)),
        "multiplicity range drifted",
    )
    coefficient_checks = 0
    for row in rows:
        m = row["multiplicity"]
        polynomial = rebuild_polynomial(m, eta, a)
        require(str(polynomial) == row.get("universal_polynomial"), f"polynomial drift at m={m}")
        coefficients = row.get("coefficients", [])
        require(len(coefficients) == m // 2 + 1, f"term count drift at m={m}")
        require(row.get("term_count") == len(coefficients), f"stored term count drift at m={m}")
        for q, coefficient_row in enumerate(coefficients):
            expected = sp.factorial(m) // (
                sp.factorial(q) * sp.factorial(m - 2 * q)
            )
            require(coefficient_row.get("q") == q, f"q drift at m={m}")
            require(coefficient_row.get("eta_power") == m - 2 * q, f"power drift at m={m}")
            require(coefficient_row.get("integer_coefficient") == int(expected), f"coefficient drift at m={m}")
            coefficient_checks += 1
        require(
            sp.expand(polynomial.subs(a, -sp.Rational(1, 2)) - sp.hermite_prob(m, eta)) == 0,
            f"good-side Hermite drift at m={m}",
        )
        require(
            sp.expand(
                polynomial.subs(a, sp.Rational(1, 2))
                - sp.I**m * sp.hermite_prob(m, eta / sp.I)
            )
            == 0,
            f"bad-side Hermite drift at m={m}",
        )
        require(row.get("good_side_real_simple_roots") == m, f"good roots drift at m={m}")
        require(row.get("bad_side_nonreal_roots") == m - (m % 2), f"bad roots drift at m={m}")
        require(row.get("bad_side_real_roots") == m % 2, f"bad real roots drift at m={m}")
        require(row.get("collision_parameter_scale") == "tau=rho/8", f"scale drift at m={m}")
    require(coefficient_checks == EXPECTED_COUNTS["universal_coefficient_checks"], "coefficient count drifted")


def check_exact_polynomial_layers(artifact: dict) -> None:
    h, x, z, eta, rho, tau = sp.symbols("h x z eta rho tau")
    stored = artifact.get("exact_polynomial_layer_audit", [])
    require(len(stored) == EXPECTED_COUNTS["exact_polynomial_layer_checks"], "exact layer count drifted")
    position = 0
    for n in (0, 3):
        for m in range(2, 9):
            row = stored[position]
            position += 1
            require(row.get("shift") == n and row.get("multiplicity") == m, "exact layer index drifted")
            iterate = (x - rho) ** m
            heat_polynomial = sp.Integer(0)
            for order in range(m + 1):
                heat_polynomial += (tau * h**2) ** order * iterate / sp.factorial(order)
                iterate = sp.expand(
                    4 * x * sp.diff(iterate, x, 2)
                    + (4 * n + 2) * sp.diff(iterate, x)
                )
            scaled_jensen = sp.Integer(0)
            for (degree,), coefficient in sp.Poly(sp.expand(heat_polynomial), x).terms():
                ratio = sp.prod(1 - index * h**2 for index in range(degree))
                scaled_jensen += coefficient * ratio * z**degree
            exact_limit = sp.expand(
                sp.limit(
                    sp.expand(scaled_jensen.subs(z, rho + eta * h) / h**m),
                    h,
                    0,
                )
            )
            target = rebuild_polynomial(m, eta, 4 * rho * tau - rho**2 / 2)
            require(sp.expand(exact_limit - target) == 0, f"exact layer mismatch at n={n}, m={m}")
            require(row.get("matched") is True, f"stored exact layer status drifted at n={n}, m={m}")
            require(row.get("limit") == str(exact_limit), f"stored exact limit drifted at n={n}, m={m}")
            require(row.get("target") == str(target), f"stored exact target drifted at n={n}, m={m}")


def check_certificate(artifact: dict) -> None:
    certificate = artifact.get("symbolic_certificate", {})
    expected_keys = {
        "newman_map", "jensen_operator", "jensen_local_limit", "heat_local_limit",
        "combined_layer", "coefficient_formula", "good_side", "bad_side",
        "collision_scale", "positive_boundary", "fixed_order_guard", "open_target",
    }
    require(set(certificate) == expected_keys, "certificate key drifted")
    require("rho=-c^2" in certificate["newman_map"], "Newman map drifted")
    require("n=0" in certificate["positive_boundary"], "fixed shift missing")
    require("rho/(8D_j)" in certificate["collision_scale"], "collision scale missing")
    require("K_(m,a)" in certificate["combined_layer"], "universal layer missing")
    require("k=D" in certificate["fixed_order_guard"], "growing Hankel order missing")
    require("uniform in degree" in certificate["open_target"], "uniform target missing")


def check_note(note: str) -> None:
    markers = [
        "# Newman Arbitrary-Multiplicity Jensen Boundary-Layer Gate",
        "## Positive-Boundary Coordinate",
        "single shift `n=0`",
        "## Exact Jensen Operator",
        "## Universal Multiplicity Layer",
        "K_(m,a)(eta)",
        "For `m=2`",
        "t_j=Lambda-c^2/(8D_j)+o(D_j^-1)",
        "## Quantifier Guard",
        "## Pi Provenance",
        "## Proof Boundary",
    ]
    for marker in markers:
        require(marker in note, f"note marker missing: {marker}")


def main() -> int:
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    note = NOTE_PATH.read_text(encoding="utf-8")
    require(artifact.get("kind") == KIND, "kind drifted")
    require("one degree-uniform Xi handoff open" in artifact.get("status", ""), "status drifted")
    require(artifact.get("counts") == EXPECTED_COUNTS, "counts drifted")
    require("no degree-uniform exclusion" in artifact.get("proof_boundary", ""), "proof boundary drifted")
    check_source_audit(artifact)
    check_rows(artifact)
    check_multiplicities(artifact)
    check_exact_polynomial_layers(artifact)
    check_certificate(artifact)
    check_note(note)
    print(
        "validated arbitrary-multiplicity Jensen boundary-layer gate: "
        f"{EXPECTED_COUNTS['gate_rows']} rows, {EXPECTED_COUNTS['source_artifacts']} sources, "
        f"{EXPECTED_COUNTS['multiplicity_rows']} multiplicities, "
        f"{EXPECTED_COUNTS['universal_coefficient_checks']} universal coefficients, "
        f"{EXPECTED_COUNTS['hermite_good_side_checks']} Hermite checks, "
        f"{EXPECTED_COUNTS['imaginary_bad_side_checks']} imaginary-side checks, "
        f"{EXPECTED_COUNTS['orientation_checks']} orientation checks, "
        f"{EXPECTED_COUNTS['exact_polynomial_layer_checks']} exact polynomial layers, "
        f"{EXPECTED_COUNTS['fixed_jensen_shifts']} fixed shift, "
        f"{EXPECTED_COUNTS['uniform_in_degree_bounds']} uniform-degree bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
