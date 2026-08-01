#!/usr/bin/env python3
"""Independently validate the stable-remainder outer-tail theorem."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import comb, factorial
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR_DIR = REPO_ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR_DIR) not in sys.path:
    sys.path.insert(0, str(VENDOR_DIR))

import flint
from flint import arb
import sympy as sp


STEM = "jensen_window_pf_newman_theta_stable_remainder_outer_tail_gate"
DEFAULT_ARTIFACT = (
    REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
)
TAIL_GATE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_forward_remainder_tail_gate.json"
)
MODULAR_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_modular_blend_gate.json"
)
PRECISION_BITS = 256
M_ORDER = 9
N_VALUES = tuple(range(4, 11))
EXPECTED_IDS = [
    "ntsrotg_01_stable_outer_split",
    "ntsrotg_02_forward_outer_envelope",
    "ntsrotg_03_switch_derivative_envelope",
    "ntsrotg_04_reflected_kernel_envelope",
    "ntsrotg_05_heat_absorption",
    "ntsrotg_06_outer_matrix",
    "ntsrotg_07_compact_handoff",
]


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def derivative_data() -> list[dict]:
    x = sp.symbols("X", nonnegative=True)
    current = 2 * x - 3
    rows: list[dict] = []
    for q in range(M_ORDER + 1):
        polynomial = sp.Poly(sp.expand(current), x)
        rows.append(
            {
                "q": q,
                "P_q": str(polynomial.as_expr()),
                "A_q": int(
                    sum(
                        abs(value)
                        for value in polynomial.all_coeffs()
                    )
                ),
            }
        )
        current = sp.expand(
            (5 - 4 * x) * polynomial.as_expr()
            + 4 * x * sp.diff(polynomial.as_expr(), x)
        )
    return rows


def laurent_norm(expression: sp.Expr, y: sp.Symbol) -> tuple[sp.Expr, int]:
    norm = sp.Rational(0)
    maximum: int | None = None
    for term in sp.Add.make_args(sp.expand(expression)):
        coefficient, exponent = term.as_coeff_exponent(y)
        if coefficient.has(y) or not exponent.is_integer:
            raise RuntimeError(f"non-Laurent switch term: {term}")
        norm += abs(coefficient)
        maximum = (
            int(exponent)
            if maximum is None
            else max(maximum, int(exponent))
        )
    if maximum is None:
        raise RuntimeError("empty switch Laurent polynomial")
    return norm, maximum


def switch_data() -> list[dict]:
    y = sp.symbols("y", positive=True)
    z = sp.Rational(3, 2) * (y - y**-1)
    z_prime = 4 * y * sp.diff(z, y)
    current = sp.expand(-z_prime)
    rows = [
        {
            "a": 0,
            "laurent_polynomial": None,
            "coefficient_norm": "1/2",
            "maximum_degree": 0,
            "normalization": "1",
            "bound": (
                "0<s(u)<=0.5*exp(-(9/16)*exp(8u))"
            ),
        }
    ]
    for order in range(1, M_ORDER + 1):
        norm, maximum = laurent_norm(current, y)
        rows.append(
            {
                "a": order,
                "laurent_polynomial": str(current),
                "coefficient_norm": str(norm),
                "maximum_degree": maximum,
                "normalization": "1/sqrt(pi)",
                "bound": (
                    "|s^(a)(u)|<=C_a*exp(4*h_a*u)"
                    "*exp(-(9/16)*exp(8u))/sqrt(pi)"
                ),
            }
        )
        current = sp.expand(
            4 * y * sp.diff(current, y)
            - 2 * z * z_prime * current
        )
    return rows


def specs() -> list[dict]:
    rows: list[dict] = []
    for b in range(M_ORDER + 1):
        rows.extend(
            [
                {
                    "component": "d0",
                    "p": 0,
                    "b": b,
                    "q": M_ORDER - b,
                    "coefficient": comb(M_ORDER, b),
                },
                {
                    "component": "d1_u",
                    "p": 1,
                    "b": b,
                    "q": M_ORDER - b,
                    "coefficient": comb(M_ORDER, b),
                },
            ]
        )
    for b in range(M_ORDER):
        rows.append(
            {
                "component": "d1_lower",
                "p": 0,
                "b": b,
                "q": M_ORDER - 1 - b,
                "coefficient": M_ORDER * comb(M_ORDER - 1, b),
            }
        )
    return rows


def rational_ball(text: str) -> arb:
    value = sp.Rational(text)
    return arb(int(value.p)) / int(value.q)


def heat_polynomial(order: int, t: arb, u: arb) -> arb:
    return sum(
        (
            arb(factorial(order))
            * (2 * u) ** (order - 2 * ell)
            * t ** (order - ell)
            / (factorial(ell) * factorial(order - 2 * ell))
            for ell in range(order // 2 + 1)
        ),
        arb(0),
    )


def serialize(value: arb) -> dict:
    log10_value = value.log() / arb(10).log()
    return {
        "enclosure": value.str(65, more=True),
        "upper": value.upper().str(65),
        "log10_enclosure": log10_value.str(50, more=True),
        "relative_accuracy_bits": value.rel_accuracy_bits(),
    }


def rebuild_entry(
    retained: int,
    spec: dict,
    polynomials: list[dict],
    switches: list[dict],
) -> tuple[dict, arb]:
    u0 = arb(11) / 5
    t = arb(1) / 5
    y4 = (4 * u0).exp()
    y8 = (8 * u0).exp()
    pi = arb.pi()
    p, b, q = spec["p"], spec["b"], spec["q"]
    first = retained + 1
    degree = 2 * q + 4
    linear = 4 * q + 9
    rho = (
        arb(degree) / first - pi * y4 * (2 * first + 1)
    ).exp()
    first_value = (
        polynomials[q]["A_q"]
        * pi ** (q + 2)
        * first**degree
        * (linear * u0 - pi * first**2 * y4).exp()
    )
    weight = (
        u0**p
        * heat_polynomial(b, t, u0)
        * (t * u0**2).exp()
    )
    forward_decay = (
        4 * pi * first**2 * y4
        - 2 * t * u0
        - linear
        - arb(p + b) / u0
    )
    if rho.upper() >= 1 or forward_decay.lower() <= 0:
        raise RuntimeError("independent forward decay check failed")
    forward = weight * first_value / (1 - rho) / forward_decay

    c = arb(9) / 16
    defect = arb(0)
    term_count = 0
    for index in range(1, retained + 1):
        if (pi * index**2 / y4).upper() >= 1:
            raise RuntimeError("independent reflected-X check failed")
        for switch_order in range(q + 1):
            kernel_order = q - switch_order
            a_k = polynomials[kernel_order]["A_q"]
            kernel_linear = 4 * kernel_order + 9
            forward_at_split = (
                a_k
                * pi ** (kernel_order + 2)
                * index ** (2 * kernel_order + 4)
                * (
                    kernel_linear * u0
                    - pi * index**2 * y4
                ).exp()
            )
            raw_decay = 4 * pi * index**2 * y4 - kernel_linear
            if raw_decay.lower() <= 5:
                raise RuntimeError(
                    "independent kernel decay check failed"
                )
            difference = (
                a_k * pi * index**2
                + forward_at_split * (5 * u0).exp()
            )
            switch = switches[switch_order]
            h = int(switch["maximum_degree"])
            switch_constant = rational_ball(
                switch["coefficient_norm"]
            )
            if switch_order:
                switch_constant /= pi.sqrt()
            switch_linear = 4 * h - 5
            initial = (
                comb(q, switch_order)
                * switch_constant
                * difference
                * u0**p
                * heat_polynomial(b, t, u0)
                * (
                    t * u0**2
                    + switch_linear * u0
                    - c * y8
                ).exp()
            )
            decay = (
                8 * c * y8
                - 2 * t * u0
                - switch_linear
                - arb(p + b) / u0
            )
            if decay.lower() <= 0:
                raise RuntimeError(
                    "independent switch decay check failed"
                )
            defect += initial / decay
            term_count += 1
    total = forward + defect
    entry = {
        **spec,
        "N": retained,
        "forward_outer": serialize(forward),
        "retained_defect_outer": serialize(defect),
        "total_outer": serialize(total),
        "forward_diagnostic": {
            "first_forward_index": first,
            "degree": degree,
            "linear_exponent": linear,
            "rho_enclosure": rho.str(45, more=True),
            "decay_rate_enclosure": forward_decay.str(45, more=True),
        },
        "defect_term_count": term_count,
    }
    return entry, total


def rebuild_numerics() -> tuple[list[dict], list[dict]]:
    polynomials = derivative_data()
    switches = switch_data()
    entries: list[dict] = []
    aggregates: list[dict] = []
    for retained in N_VALUES:
        components = {
            "d0": arb(0),
            "d1_u": arb(0),
            "d1_lower": arb(0),
        }
        for spec in specs():
            entry, total = rebuild_entry(
                retained, spec, polynomials, switches
            )
            entries.append(entry)
            components[spec["component"]] += (
                spec["coefficient"] * total
            )
        d1 = components["d1_u"] + components["d1_lower"]
        aggregates.append(
            {
                "N": retained,
                "outer_d0": serialize(components["d0"]),
                "outer_d1_u": serialize(components["d1_u"]),
                "outer_d1_lower": serialize(
                    components["d1_lower"]
                ),
                "outer_d1": serialize(d1),
            }
        )
    return entries, aggregates


def validate(path: Path) -> list[str]:
    flint.ctx.prec = PRECISION_BITS
    issues: list[str] = []
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("kind") != STEM:
        issues.append("artifact kind mismatch")
    parameters = artifact.get("parameters", {})
    expected_parameters = {
        "precision_bits": PRECISION_BITS,
        "m_order": M_ORDER,
        "n_values": list(N_VALUES),
        "u_split_exact": "11/5",
        "t_cap_exact": "1/5",
        "forward_phase": "pi*n^2*exp(4u)",
        "switch_phase_lower": "(9/16)*exp(8u)",
    }
    if parameters != expected_parameters:
        issues.append("parameter mismatch")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    if [row.get("readiness") for row in rows] != (
        ["proved"] * 6 + ["not_ready_to_apply"]
    ):
        issues.append("row readiness mismatch")

    audit = artifact.get("source_audit", {})
    if audit.get("tail_gate_sha256") != file_hash(TAIL_GATE):
        issues.append("tail-gate source hash mismatch")
    if audit.get("modular_source_sha256") != file_hash(
        MODULAR_SOURCE
    ):
        issues.append("modular source hash mismatch")

    rebuilt_polynomials = derivative_data()
    if artifact.get("derivative_polynomials") != rebuilt_polynomials:
        issues.append("derivative-polynomial recurrence mismatch")
    rebuilt_switches = switch_data()
    if artifact.get("switch_derivatives") != rebuilt_switches:
        issues.append("switch Laurent recurrence mismatch")

    y = sp.symbols("y", positive=True)
    phase_slack = sp.factor(
        9 * ((y - y**-1) / 2) ** 2
        - sp.Rational(9, 16) * y**2
    )
    expected_slack = (
        9 * (3 * y**2 - 2) * (y**2 - 2) / (16 * y**2)
    )
    if sp.simplify(phase_slack - expected_slack) != 0:
        issues.append("switch phase factorization failed")

    rebuilt_entries, rebuilt_aggregates = rebuild_numerics()
    if artifact.get("outer_entries") != rebuilt_entries:
        issues.append("independent outer-entry reconstruction mismatch")
    if artifact.get("outer_aggregates") != rebuilt_aggregates:
        issues.append("independent outer aggregation mismatch")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "complete arithmetic remainder",
        "u>11/5",
        "compact matrix",
        "first-jet",
        "transition-cell",
        "Lambda<=0",
        "RH",
        "Clay-prize",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary guard missing: {marker}")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    args = parser.parse_args()
    issues = validate(args.artifact)
    if issues:
        for issue in issues:
            print(f"ISSUE: {issue}")
        raise SystemExit(1)
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    print(
        "validated Newman theta stable-remainder outer-tail gate: "
        f"{len(artifact['rows'])} rows, "
        f"{len(artifact['outer_entries'])} directed entries, "
        f"{len(artifact['outer_aggregates'])} retained counts, 0 issues"
    )


if __name__ == "__main__":
    main()
