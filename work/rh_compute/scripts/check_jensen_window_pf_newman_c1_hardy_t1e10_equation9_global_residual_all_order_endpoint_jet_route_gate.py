#!/usr/bin/env python3
"""Independently check the all-order endpoint-jet route gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_all_order_endpoint_jet_route_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")

A = 159_577
B = 5_122_421
T_LO = 622
T_HI = 39_894


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def independent_pair_check() -> None:
    d = sp.symbols("d", nonzero=True)
    jumps = sp.symbols("j0:10")
    r_plus, r_minus = sp.symbols("r_plus r_minus")
    for order in (2, 4, 6, 8, 10):
        plus = -sum(jumps[j] / d ** (j + 1) for j in range(order)) + r_plus / d**order
        minus = -sum(jumps[j] / (-d) ** (j + 1) for j in range(order)) + r_minus / (-d) ** order
        target = -2 * sum(jumps[j] / d ** (j + 1) for j in range(1, order, 2))
        target += (r_plus + r_minus) / d**order
        require(sp.cancel(plus + minus - target) == 0, f"independent pair failure at {order}")

    y, q = sp.symbols("y q")
    p = y
    for n in range(11):
        poly = sp.Poly(p, y)
        require(poly.degree() == n + 1, f"independent degree failure at {n}")
        require(sp.expand(poly.LC() - q**n) == 0, f"independent leading coefficient failure at {n}")
        p = sp.expand(2 * sp.diff(p, y) + q * y * p)


def independent_scout_row(mode: int, x: mp.mpf, endpoint: int) -> list[mp.mpf]:
    mp.mp.dps = 90
    q = mp.j * mp.pi * x
    coefficients = [mp.mpc(0), mp.mpc(1)]
    polynomial_rows = [coefficients]
    for _ in range(7):
        derivative = [(j + 1) * coefficients[j + 1] for j in range(len(coefficients) - 1)]
        following = [mp.mpc(0)] * max(len(derivative), len(coefficients) + 1)
        for j, value in enumerate(derivative):
            following[j] += 2 * value
        for j, value in enumerate(coefficients):
            following[j + 1] += q * value
        coefficients = following
        polynomial_rows.append(coefficients)
    fourier_d = 2 * mp.pi * mp.j * mode
    values: list[mp.mpf] = []
    for n in (1, 3, 5, 7):
        value = mp.mpc(0)
        for coefficient in reversed(polynomial_rows[n]):
            value = value * endpoint + coefficient
        value *= mp.exp(q * endpoint**2 / 4)
        values.append(abs(2 * value / fourier_d ** (n + 1)))
    return values


def main() -> int:
    mp.mp.dps = 90
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True and artifact.get("kind") == STEM, "artifact identity drift")
    decisions = artifact["decision"]
    for key in (
        "all_order_signed_pair_identity_proved",
        "only_odd_endpoint_jumps_survive_pairing",
        "endpoint_derivative_recurrence_proved",
        "highest_jet_ratio_proved",
        "B_target_half_cell_ratio_above_32",
        "A_half_boundary_bracketed_by_39894_39895",
        "phase_adapted_transition_charts_still_required",
    ):
        require(decisions.get(key) is True, f"missing decision: {key}")
    for key in (
        "raw_higher_jet_absolute_hierarchy_uniformly_descending",
        "cancellation_aware_resummation_excluded",
        "quantitative_R_after_A_bound_proved",
        "rh_implication",
    ):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    independent_pair_check()
    require(Fraction(B, 4 * T_HI) > 32, "independent B ratio failure")
    require(Fraction(A, 4 * T_HI) > 1, "independent A target ratio failure")
    require(Fraction(A, 4 * (T_HI + 1)) < 1, "independent A outer ratio failure")
    require(Fraction(2 * T_LO, B) * B == 2 * T_LO, "independent B crossing failure")

    rows = artifact["diagnostic_order_scout"]["rows"]
    cases = (
        (T_LO, mp.mpf(2 * T_LO) / B, B),
        (T_HI, mp.mpf("0.5"), A),
        (T_HI + 1, mp.mpf("0.5"), A),
    )
    require(len(rows) == len(cases), "diagnostic row count drift")
    for row, (mode, x, endpoint) in zip(rows, cases, strict=True):
        expected = independent_scout_row(mode, x, endpoint)
        saved = [mp.mpf(term["endpoint_term_modulus"]) for term in row["terms"]]
        require(len(saved) == 4, "diagnostic term count drift")
        for observed, reference in zip(saved, expected, strict=True):
            require(abs(observed / reference - 1) < mp.mpf("1e-20"), "diagnostic replay drift")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(BUILDER.is_file() and NOTE.is_file(), "builder or note missing")
    note = NOTE.read_text(encoding="utf-8")
    require("every even endpoint derivative cancels" in note, "pair cancellation statement missing")
    require("endpoint jumps remain" in note, "odd-jump statement missing")
    require("cannot give a uniform descending absolute hierarchy" in note, "route decision missing")
    require("No lower bound on the true signed" in note, "proof boundary missing")
    print("independently checked all-order endpoint-jet route guard", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
