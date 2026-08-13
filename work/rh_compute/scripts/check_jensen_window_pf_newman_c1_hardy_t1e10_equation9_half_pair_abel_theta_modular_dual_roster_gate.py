#!/usr/bin/env python3
"""Independent replay of the Abel-theta modular dual-roster gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_modular_dual_roster_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()

A = 159_577
B = 5_122_421
L = (B - A) // 2


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def theta(tau: mp.mpc, cutoff: int) -> mp.mpc:
    terms = [mp.mpf(1)]
    terms.extend(2 * mp.exp(mp.j * mp.pi * tau * n * n) for n in range(1, cutoff + 1))
    return mp.fsum(terms)


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")

    # Replay the saddle calculation in face-distance y=x-z coordinates,
    # independently of the production z-coordinate differentiation.
    x, y, endpoint, n, height = sp.symbols("x y D n t", positive=True, real=True)
    phase = (
        sp.pi * endpoint**2 * (x - y) / 4
        + height * sp.log((1 - x) / x) / 2
        - sp.pi * n**2 * x * (x - y) / y
    )
    y_star = 2 * n * x / endpoint
    require(sp.simplify(sp.diff(phase, y).subs(y, y_star)) == 0, "face-distance saddle replay failed")
    alpha = endpoint - 2 * n
    tangential = sp.factor(sp.simplify(sp.diff(phase, x).subs(y, y_star)))
    expected = sp.pi * alpha**2 / 4 - height / (2 * x * (1 - x))
    require(sp.simplify(tangential - expected) == 0, "dual tangential saddle replay failed")

    require(B - 2 * L == A, "B roster does not terminate at A")
    roster = range(B, A - 1, -2)
    require(len(roster) == L + 1, "B roster cardinality drift")
    require(roster[0] == B and roster[-1] == A, "B roster endpoint drift")
    for k in (0, 1, 97, 10_000):
        require(B - 2 * (L + k) == A - 2 * k, "continuation-label identity failed")

    # Use a different point and truncation from the production witness.
    mp.mp.dps = 110
    tau = mp.mpc("-0.219", "0.311")
    direct = theta(tau, 36)
    modular = (-mp.j * tau) ** (-mp.mpf("0.5")) * theta(-1 / tau, 36)
    require(abs(direct - modular) < mp.mpf("1e-95"), "independent modular witness failed")

    decision = artifact["decision"]
    require(decision["nonzero_pairs_reassembled_as_one_Abel_theta_triangle"] is True, "theta reduction drift")
    require(decision["endpoint_difference_retained_before_majorization"] is True, "endpoint grouping drift")
    require(decision["B_dual_indices_recover_complete_odd_source_roster"] is True, "dual roster drift")
    require(decision["continuation_endpoints_cancel_termwise"] is False, "unproved continuation cancellation")
    require(decision["modular_transform_alone_proves_quantitative_gain"] is False, "modular overclaim")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    print("checked Abel-theta modular roster in independent face-distance coordinates", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
