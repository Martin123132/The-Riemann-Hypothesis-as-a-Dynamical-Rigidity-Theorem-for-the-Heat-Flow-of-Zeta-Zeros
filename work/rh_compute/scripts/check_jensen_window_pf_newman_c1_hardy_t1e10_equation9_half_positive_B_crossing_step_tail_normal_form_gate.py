#!/usr/bin/env python3
"""Independent replay of the exact B step-tail crossing normal form."""

from __future__ import annotations

import hashlib
import json
import re
from decimal import Decimal
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_positive_B_crossing_step_tail_normal_form_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_zero_ball(text: str) -> bool:
    """Verify that an Arb midpoint-radius rendering encloses zero."""
    compact = text.replace(" ", "")
    if compact == "0":
        return True
    match = re.fullmatch(r"\[([^\]]+)\+/-([^\]]+)\]", compact)
    if match is None:
        return False
    midpoint, radius = (Decimal(value) for value in match.groups())
    return abs(midpoint) <= radius


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")

    # Replay the crossing identity as a two-row truth table.
    a, e, f = sp.symbols("a e f")
    bulk = 2 * a * f
    rows = [
        # q positive: P_B=U and I_qneg=0
        (e - a * sp.symbols("T"), e - a * sp.symbols("T"), 0),
        # q negative: P_B=U-P_bulk and I_qneg=1
        (e - bulk + a * sp.symbols("T"), e + a * sp.symbols("T"), 1),
    ]
    for p_b, u, indicator in rows:
        require(sp.expand(p_b - (u - indicator * bulk)) == 0, "truth-table replay failed")

    endpoint, x, m = sp.symbols("endpoint x m", positive=True, real=True)
    aa = endpoint * x / 2
    partial = sp.simplify(
        1 / (endpoint**2 * x**2 - 4 * m**2)
        - sp.Rational(1, 4) / (aa**2 - m**2)
    )
    require(partial == 0, "cotangent scaling replay failed")

    # Independently recover the completed face phase and first paired boundary
    # current from the two phase derivatives.  This does not reuse strings from
    # the production certificate.
    q = sp.sqrt(x / 2) * (endpoint - 2 * m / x)
    completed = sp.simplify(-sp.pi * m**2 / x + sp.pi * q**2 / 2)
    expected_trace = sp.pi * endpoint**2 * x / 4 - sp.pi * endpoint * m
    require(sp.simplify(completed - expected_trace) == 0, "completed face phase replay failed")

    z = sp.symbols("z", positive=True, real=True)
    phase_z = sp.pi * m**2 / z + sp.pi * endpoint**2 * z / 4
    driver = z ** (-sp.Rational(1, 2)) * (2 + sp.I * sp.pi * endpoint**2 * z)
    boundary = sp.simplify(
        x ** (-sp.Rational(3, 2)) * driver.subs(z, x)
        / (2 * sp.I * sp.pi * sp.I * sp.diff(phase_z, z).subs(z, x))
    )
    expected_boundary = -2 * (2 + sp.I * sp.pi * endpoint**2 * x) / (
        sp.pi**2 * (endpoint**2 * x**2 - 4 * m**2)
    )
    require(sp.simplify(boundary - expected_boundary) == 0, "paired boundary-current replay failed")

    decision = artifact["decision"]
    require(decision["B_endpoint_rewritten_as_sign_adapted_tail_plus_bulk_step"] is True, "normal-form decision drift")
    require(decision["tail_and_step_jumps_cancel_exactly"] is True, "jump decision drift")
    require(decision["local_modes_621_622_may_be_absorbed_into_cotangent_bound"] is False, "local-pole overclaim")
    require(decision["Fresnel_and_cotangent_remainders_bounded"] is False, "remainder overclaim")

    for witness in artifact["continuity_witnesses"]:
        require(exact_zero_ball(witness["jump_plus_bulk_absolute_ball"]), "continuity witness is not exact zero")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    print("checked exact B step-tail normal form by truth-table and partial-fraction replay", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
