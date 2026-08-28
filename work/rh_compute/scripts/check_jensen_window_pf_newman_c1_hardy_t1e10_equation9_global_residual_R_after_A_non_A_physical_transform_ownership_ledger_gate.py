#!/usr/bin/env python3
"""Independently check the non-A physical-transform ownership ledger."""

from __future__ import annotations

from decimal import Decimal, getcontext
import hashlib
import json
import os
from pathlib import Path
import re
import sys
from typing import Any


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
from mpmath import mp


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_physical_transform_ownership_ledger_gate"
)
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(Path(__file__).name.removeprefix("check_"))
BALL_PATTERN = re.compile(
    r"^\[?\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
    r"(?:\s*\+/-\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?))?\s*\]?$"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def bounds(text: str) -> tuple[Decimal, Decimal]:
    match = BALL_PATTERN.match(text.strip())
    require(match is not None, f"cannot parse ball: {text}")
    midpoint = Decimal(match.group(1))
    radius = Decimal(match.group(2) or "0")
    require(radius >= 0, "negative radius")
    return midpoint - radius, midpoint + radius


def absolute_upper(text: str) -> Decimal:
    lower, upper = bounds(text)
    return max(abs(lower), abs(upper))


def independent_b_trace(dependencies: dict[str, dict[str, Any]]) -> tuple[Decimal, Decimal]:
    local = dependencies["B_local"]["interval_certificate"]
    window = dependencies["B_window"]["interval_certificate"]
    principal_lower, principal_upper = bounds(
        local["combined_principal_physical_signed_projection_ball"]
    )
    signed_lower, signed_upper = bounds(window["signed_nonlocal_first_current_ball"])
    radius = sum(
        (
            absolute_upper(local["positive_normal_remainder_absolute_allowance_ball"]),
            absolute_upper(local["negative_partner_remainder_absolute_allowance_ball"]),
            absolute_upper(window["nonlocal_higher_current_absolute_allowance_ball"]),
            absolute_upper(window["nonlocal_post_third_current_absolute_allowance_ball"]),
            absolute_upper(window["Fresnel_boundary_dictionary_absolute_allowance_ball"]),
        ),
        Decimal(0),
    )
    return principal_lower + signed_lower - radius, principal_upper + signed_upper + radius


def check_kummer_euler_identity() -> str:
    mp.dps = 80
    t = mp.mpf("7.25")
    alpha = mp.mpf(5)
    a = mp.mpf(3) / 4 - mp.j * t / 2
    b = mp.mpf(3) / 4 + mp.j * t / 2
    z = mp.j * mp.pi * alpha**2 / 4
    integrand = lambda x: x ** (a - 1) * (1 - x) ** (b - 1) * mp.exp(z * x)
    direct = mp.quad(integrand, [0, mp.mpf("0.25"), mp.mpf("0.5"), mp.mpf("0.75"), 1])
    closed = mp.beta(a, b) * mp.hyp1f1(a, mp.mpf("1.5"), z)
    relative = abs(direct - closed) / abs(closed)
    require(relative < mp.mpf("1e-55"), "independent Euler-Kummer spot check failed")
    return mp.nstr(relative, 35)


def main() -> int:
    getcontext().prec = 130
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True, "builder artifact failed")
    require(artifact["kind"] == STEM, "artifact kind drift")
    decision = artifact["decision"]
    require(decision["complete_source_transform_identified_exactly_as_Q_K"] is True, "Q_K drift")
    require(decision["single_joined_J_Z_transform_isolated"] is True, "J_Z drift")
    require(decision["source_hybrid_telemetry_transfer_permitted"] is False, "telemetry guard lost")
    require(decision["J_Z_interval_enclosed"] is False, "J_Z overclaim")
    require(decision["non_A_bound_proved"] is False, "non-A overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")

    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(Path(__file__)) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    dependencies: dict[str, dict[str, Any]] = {}
    for name, record in artifact["dependencies"].items():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {name}")
        dependencies[name] = load_json(path)

    b_lower, b_upper = independent_b_trace(dependencies)
    stored_lower, stored_upper = bounds(
        artifact["B_trace_two_sided_certificate"]["derived_two_sided_ball"]
    )
    require(stored_lower <= b_lower <= b_upper <= stored_upper, "stored B interval misses Decimal replay")
    require(b_lower > Decimal("-0.000134948"), "B lower scale guard failed")
    require(b_upper < Decimal("-0.00013198"), "B upper scale guard failed")

    a_text = dependencies["A_transition"]["certificate"]["components"][
        "projector_completed_exact_A_transition_ball"
    ]
    a_lower, a_upper = bounds(a_text)
    outer = Decimal("8e-10")
    face = Decimal("0.00364")
    target = Decimal("0.0331747039947")
    known_lower = b_lower - outer - face
    known_upper = b_upper + outer + face
    exact_j_lower = known_upper - target
    exact_j_upper = known_lower + target
    require(Decimal("-0.02966668") > exact_j_lower, "published J_Z lower endpoint failed")
    require(Decimal("0.02939975") < exact_j_upper, "published J_Z upper endpoint failed")
    require(-Decimal("0.02939975") > exact_j_lower, "symmetric J_Z lower endpoint failed")

    rkg_lower = a_upper + known_upper - target
    rkg_upper = a_lower + known_lower + target
    require(Decimal("-0.0663408") > rkg_lower, "published R_KGamma lower endpoint failed")
    require(Decimal("-0.00727437") < rkg_upper, "published R_KGamma upper endpoint failed")
    gamma = Decimal("5e-20")
    require(Decimal("-0.0663408") > rkg_lower + gamma, "Q_K-T lower correction failed")
    require(Decimal("-0.00727437") < rkg_upper - gamma, "Q_K-T upper correction failed")

    identity = artifact["identity_certificate"]
    require(identity["joined_source_owned_transform"] == "J_Z=P_t[Z]=Q_K-G-A_transition=R_KGamma-A_transition", "J_Z identity text drift")
    require(identity["non_A_identity"] == "R_nonA=J_Z-E_Btr,win-E_outer-I_(A,42)", "non-A sign drift")
    relative = check_kummer_euler_identity()
    print(
        "independently checked non-A physical-transform ledger; "
        f"Euler-Kummer surrogate relative discrepancy {relative}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
