#!/usr/bin/env python3
"""Independent replay of the half-Kummer reflection and branch reduction."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_kummer_reflection_branch_reduction_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def integrand(alpha: int, t: mp.mpf, x: mp.mpf) -> mp.mpc:
    return (
        alpha
        * mp.exp(mp.j * mp.pi * alpha * alpha * x / 4)
        * mp.exp(mp.j * t * mp.log((1 - x) / x) / 2)
        / (x * (1 - x)) ** 0.25
    )


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "half-reflection artifact is not passed")
    require(artifact["sources"]["builder"]["sha256"] == file_hash(BUILDER), "builder hash drift")
    require(artifact["sources"]["checker"]["sha256"] == file_hash(CHECKER), "checker hash drift")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(dependency["sha256"] == file_hash(path), f"dependency hash drift: {path}")

    mp.mp.dps = 80
    for alpha in (1, 3, 159_577, 5_122_421):
        for t in (mp.mpf(17), mp.mpf("1234.5")):
            for x in (mp.mpf("0.071"), mp.mpf("0.223"), mp.mpf("0.499")):
                left = integrand(alpha, t, 1 - x)
                right = mp.exp(mp.j * mp.pi / 4) * mp.conj(integrand(alpha, t, x))
                scale = max(mp.mpf(1), abs(left), abs(right))
                require(abs(left - right) < mp.mpf("1e-60") * scale, f"numeric reflection failed at alpha={alpha}, x={x}")

    t = 10_000_000_000.0
    nt = math.sqrt(t / (2 * math.pi))
    x_mode = lambda mode: 2 * math.pi * mode * mode / (t + 2 * math.pi * mode * mode)
    alpha_mode = lambda mode: 2 * mode + t / (math.pi * mode)
    require(39_894 < nt < 39_895, "N_t replay failed")
    require(x_mode(39_894) < 0.5 < x_mode(39_895), "half-domain split replay failed")
    require(alpha_mode(621) > 5_122_421 > alpha_mode(622), "B endpoint crossing replay failed")
    require(alpha_mode(39_852) > 159_577 > alpha_mode(39_853), "A endpoint crossing replay failed")

    for mode in (622, 1000, 10_000, 39_894):
        partner = t / (2 * math.pi * mode)
        require(abs((2 * partner + t / (math.pi * partner)) - alpha_mode(mode)) < 2e-8, "alpha involution replay failed")
        require(abs(x_mode(partner) - (1 - x_mode(mode))) < 2e-15, "x involution replay failed")

    decision = artifact["decision"]
    require(decision["odd_alpha_half_kummer_reflection_proved"] is True, "reflection decision drift")
    require(decision["half_domain_has_one_positive_stationary_branch"] is True, "branch decision drift")
    require(decision["reflected_upper_branch_requires_separate_stationary_estimate"] is False, "reflected-branch overclaim")
    require(decision["half_domain_nonstationary_completion_bounded"] is False, "quantitative completion overclaim")
    require(decision["complete_T_upper_proved"] is False, "T_upper overclaim")
    print("checked odd-alpha reflection and one-branch half-domain saddle ledger", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
