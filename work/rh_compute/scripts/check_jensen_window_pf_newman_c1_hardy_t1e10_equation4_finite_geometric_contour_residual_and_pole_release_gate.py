#!/usr/bin/env python3
"""Independently check the finite RSI residual and pole-release identity."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_geometric_contour_residual_and_pole_release_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def direct_kernel_check(z: mp.mpc, terms: int) -> mp.mpc:
    # Use the equivalent finite geometric polynomial in reverse order.
    prefix = mp.fsum(
        mp.e ** (1j * mp.pi * (2 * k + 1) * z)
        for k in range(terms - 1, -1, -1)
    )
    return mp.csc(mp.pi * z) + 2j * prefix - mp.e ** (2j * mp.pi * terms * z) * mp.csc(mp.pi * z)


def independent_contour(center: mp.mpf, terms: int, s: mp.mpc) -> mp.mpc:
    rotation = (1 - 1j) / mp.sqrt(2)

    def integrand(q: mp.mpf) -> mp.mpc:
        u = center + q * rotation
        return (
            mp.e ** (-1j * mp.pi * u**2)
            * mp.power(u + terms, -s)
            / mp.sin(mp.pi * u)
            * rotation
            / (2j)
        )

    return mp.quad(integrand, [-mp.inf, -5, -2, 0, 2, 5, mp.inf])


def main() -> int:
    require(gate.RESULT.is_file(), "missing finite residual result")
    require(gate.NOTE.is_file(), "missing finite residual note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == gate.STEM, "kind drift")
    require(
        artifact["status"]
        == "finite_RSI_geometric_prefix_residual_translation_and_pole_release_identity_certified",
        "status drift",
    )
    require(artifact.get("passed") is True, "finite residual gate failed")

    mp.mp.dps = 75
    for z in (mp.mpc("0.19", "0.41"), mp.mpc("0.73", "-0.16")):
        for terms in (1, 4, 9):
            require(abs(direct_kernel_check(z, terms)) < mp.mpf("1e-68"), "fresh kernel identity failed")

    s = mp.mpf("0.5") + mp.mpf("1.75") * 1j
    for terms in (2, 4):
        left = independent_contour(mp.mpf("0.5") - terms, terms, s)
        right = independent_contour(mp.mpf("0.5"), terms, s)
        dirichlet = mp.fsum(mp.power(n, -s) for n in range(1, terms + 1))
        require(abs(left - right - dirichlet) < mp.mpf("1e-60"), "fresh pole-release orientation failed")

    translated = artifact["translated_residual_and_poles"]
    require(translated["crossed_poles"] == "m=1-N,...,0", "pole roster drift")
    require(translated["residue_at_m"] == "(m+N)^(-s)/pi", "residue formula drift")
    require(translated["branch_point_minus_N_outside_deformation_strip"] is True, "branch guard lost")

    window = artifact["actual_window_identity"]
    require(window["prefix_before_window"] == (window["alpha_min"] - 1) // 2, "lower prefix arithmetic failed")
    require(window["prefix_through_window"] == (window["alpha_max"] + 1) // 2, "upper prefix arithmetic failed")
    require(window["prefix_through_window"] - window["prefix_before_window"] == window["window_label_count"], "window count failed")
    require(window["window_label_count"] == 2_481_423, "actual label count drift")
    require(window["released_Dirichlet_block"] == [79_789, 2_561_211], "released block drift")

    decision = artifact["decision"]
    require(decision["A21_infinite_interchange_required_for_finite_prefix"] is False, "A21 was reintroduced")
    require(decision["finite_prefix_plus_residual_exact"] is True, "finite residual identity lost")
    require(decision["integer_translation_releases_Dirichlet_prefix_exactly"] is True, "pole release lost")
    require(decision["physical_QK_prefactor_and_branch_alignment_complete"] is False, "Q_K alignment overclaim")
    require(decision["released_pole_block_joined_to_G_and_A_transition"] is False, "pole block overclaim")
    require(decision["D_K_enclosed"] is False, "D_K overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")

    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file {path}")
            require(file_hash(path) == record["sha256"], f"hash drift {path}")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in (
        "not a proof of the continuation defect",
        "infinite interchange is unnecessary",
        "Exactly the poles",
        "infinite remainder vanishes is needed",
        "No geometric circle construction",
        "matched branch-for-branch",
    ):
        require(token in note, f"note boundary token missing: {token}")

    print("independently checked finite RSI residual, pole release, and actual-window arithmetic", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
