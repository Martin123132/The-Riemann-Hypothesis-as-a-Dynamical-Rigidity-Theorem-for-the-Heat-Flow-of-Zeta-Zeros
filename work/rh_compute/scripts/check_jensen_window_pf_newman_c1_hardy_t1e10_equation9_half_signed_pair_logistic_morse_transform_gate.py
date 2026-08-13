#!/usr/bin/env python3
"""Independent replay of the signed-pair half-domain Morse transform."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_signed_pair_logistic_morse_transform_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")

    # Replay square completion with a signed integer variable rather than the
    # builder's two explicit branches.
    alpha, x, n = sp.symbols("alpha x n", real=True, nonzero=True)
    q = sp.sqrt(x / 2) * (alpha - 2 * n / x)
    direct = sp.pi * x * alpha**2 / 4 - sp.pi * n * alpha
    completed = -sp.pi * n**2 / x + sp.pi * q**2 / 2
    require(sp.simplify(direct - completed) == 0, "signed square completion failed")
    require(sp.simplify(completed.subs(n, -n) + sp.pi * n**2 / x - sp.pi * (sp.sqrt(x / 2) * (alpha + 2 * n / x)) ** 2 / 2) == 0, "negative branch replay failed")

    mp.mp.dps = 90
    t = mp.mpf(10) ** 10
    nt = mp.sqrt(t / (2 * mp.pi))
    require(mp.mpf(39_894) < nt < mp.mpf(39_895), "threshold interval failed")
    for row in artifact["finite_height_ledger"]["rows"]:
        mode = mp.mpf(row["mode"])
        v_half = 2 * mp.pi * mode**2 / t
        s_half = mp.sign(v_half - 1) * mp.sqrt(2 * (v_half - 1 - mp.log(v_half)))
        require(abs(v_half - mp.mpf(row["v_H"])) < mp.mpf("1e-68"), "v_H replay drift")
        require(abs(s_half - mp.mpf(row["s_H"])) < mp.mpf("1e-68"), "s_H replay drift")

    decision = artifact["decision"]
    require(decision["positive_and_negative_modes_share_outer_Morse_phase"] is True, "common phase drift")
    require(decision["signed_Fresnel_coordinates_are_identical"] is False, "signed-coordinate overclaim")
    require(decision["signed_pair_formed_inside_one_half_domain_Gaussian_integral"] is True, "pair transform drift")
    require(decision["negative_partner_is_second_stationary_branch"] is False, "stationary-branch overclaim")
    require(decision["pair_amplitude_quantitatively_bounded"] is False, "quantitative overclaim")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    print("checked signed-pair half-domain Morse transform by signed-variable replay", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
