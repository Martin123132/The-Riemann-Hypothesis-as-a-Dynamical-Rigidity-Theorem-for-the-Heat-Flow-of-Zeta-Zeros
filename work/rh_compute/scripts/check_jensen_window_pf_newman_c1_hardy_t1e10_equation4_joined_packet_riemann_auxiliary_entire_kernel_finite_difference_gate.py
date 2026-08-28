#!/usr/bin/env python3
"""Independently check the entire auxiliary-kernel finite difference."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_riemann_auxiliary_entire_kernel_finite_difference_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fresh_kappa(z: mp.mpc) -> mp.mpc:
    numerator = mp.exp(-1j * mp.pi * z) - mp.exp(-1j * mp.pi * z * z)
    return numerator / (2j * mp.sin(mp.pi * z))


def fresh_prefix(z: mp.mpc, terms: int) -> mp.mpc:
    source = mp.fsum(mp.exp(1j * mp.pi * (2 * r + 1) * z) for r in range(terms))
    return mp.exp(-1j * mp.pi * z * z) * source


def main() -> int:
    priority = gate.set_low_priority()
    require(priority == "below_normal_one_cpu", f"checker resource cap unavailable: {priority}")
    require(gate.RESULT.is_file(), "missing entire-kernel result")
    require(gate.NOTE.is_file(), "missing entire-kernel note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("kind") == gate.STEM, "kind drift")
    require(artifact.get("status") == "riemann_auxiliary_entire_kernel_finite_difference_certified", "status drift")
    require(artifact.get("passed") is True, "gate did not pass")
    for name, path in gate.DEPENDENCIES.items():
        require(digest(path) == artifact["dependencies"][name]["sha256"], f"dependency hash drift: {name}")
    require(digest(gate.BUILDER) == artifact["source_hashes"]["builder"], "builder hash drift")
    require(digest(Path(__file__)) == artifact["source_hashes"]["checker"], "checker hash drift")

    mp.mp.dps = 100
    for z in (
        mp.mpc("0.17", "0.29"),
        mp.mpc("0.91", "-0.23"),
        mp.mpc("3.2", "0.14"),
        mp.mpc("-1.13", "0.37"),
    ):
        alternate = (
            mp.exp(-0.5j * mp.pi * (z * z + z))
            * mp.sin(0.5 * mp.pi * (z * z - z))
            / mp.sin(mp.pi * z)
        )
        require(abs(fresh_kappa(z) - alternate) < mp.mpf("1e-90"), "fresh auxiliary form failed")
        for terms in (1, 3, 6, 11):
            require(
                abs(fresh_prefix(z, terms) - (fresh_kappa(z) - fresh_kappa(z - terms)))
                < mp.mpf("1e-88"),
                "fresh prefix finite difference failed",
            )

        direct_window = fresh_prefix(z, 11) - fresh_prefix(z, 3)
        entire_window = fresh_kappa(z - 3) - fresh_kappa(z - 11)
        require(abs(direct_window - entire_window) < mp.mpf("1e-87"), "fresh window identity failed")

    # Changed integers and one-sided Richardson limits independently test removability.
    for integer in (-9, -1, 2, 6, 14):
        expected = mp.mpf(integer) - mp.mpf("0.5")
        eps = mp.mpf("1e-20")
        one = fresh_kappa(mp.mpf(integer) + eps)
        half = fresh_kappa(mp.mpf(integer) + eps / 2)
        richardson = 2 * half - one
        require(abs(richardson - expected) < mp.mpf("1e-35"), "fresh integer removable value failed")

    scope = artifact["scope"]
    n_minus = scope["n_minus"]
    n_plus = scope["n_plus"]
    require(n_plus - n_minus == scope["label_count"], "actual count drift")
    removable_zero = (-mp.mpf(n_minus) - mp.mpf("0.5")) - (-mp.mpf(n_plus) - mp.mpf("0.5"))
    require(removable_zero == scope["label_count"], "actual zero extension failed")

    identities = artifact["finite_difference_identities"]
    require(identities["all_csc_integer_poles_removed_before_integration"] is True, "pole-removal drift")
    require(artifact["common_functional"]["common_base_kappa_cancels_before_norms"] is True, "common-base drift")
    decision = artifact["decision"]
    require(decision["joined_packet_enclosed"] is False, "joined-packet overclaim")
    require(decision["J_Z_enclosed"] is False, "J_Z overclaim")
    require(decision["D_K_enclosed"] is False, "D_K overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print("independently checked entire auxiliary-kernel finite difference for the joined packet")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
