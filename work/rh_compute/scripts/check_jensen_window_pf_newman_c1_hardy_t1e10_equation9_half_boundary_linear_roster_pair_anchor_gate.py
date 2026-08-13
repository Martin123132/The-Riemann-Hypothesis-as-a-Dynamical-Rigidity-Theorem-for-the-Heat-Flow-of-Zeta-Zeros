#!/usr/bin/env python3
"""Independent replay of the midpoint interpolation guard."""

from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_boundary_linear_roster_pair_anchor_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"
CHECKER = Path(__file__).resolve()
PRECISION = 130


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "midpoint guard artifact is not passed")
    require(artifact["sources"]["builder"]["sha256"] == file_hash(BUILDER), "builder hash drift")
    require(artifact["sources"]["checker"]["sha256"] == file_hash(CHECKER), "checker hash drift")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(dependency["sha256"] == file_hash(path), f"dependency hash drift: {path}")

    a, b = artifact["scope"]["source_endpoints"]
    count = artifact["scope"]["alpha_count"]
    require(all(alpha * alpha % 8 == 1 for alpha in range(a, b + 1, 2)), "discrete parity replay failed")
    require(sum(range(a, b + 1, 2)) == count * (a + b) // 2, "direct sum replay failed")

    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    imaginary_unit = acb(0, 1)
    full = (1 + imaginary_unit) / 2

    def direct_fresnel(q: arb) -> acb:
        argument = (-imaginary_unit * pi / 4).exp() * (pi / 2).sqrt() * q
        return full - full * argument.erfc()

    def mode_value(mode: int) -> acb:
        m = arb(mode)
        sign = -1 if abs(mode) % 2 else 1
        return 4 * m * sign * (direct_fresnel((arb(b) - 4 * m) / 2) - direct_fresnel((arb(a) - 4 * m) / 2))

    for row in reversed(artifact["witnesses"]):
        mode = row["mode"]
        positive = mode_value(mode)
        negative = mode_value(-mode)
        pair = positive + negative
        require(positive.overlaps(parse_complex(row["I_plus_ball"])), f"positive witness mismatch at {mode}")
        require(negative.overlaps(parse_complex(row["I_minus_ball"])), f"negative witness mismatch at {mode}")
        require(pair.overlaps(parse_complex(row["pair_ball"])), f"pair witness mismatch at {mode}")
        require(not pair.real.contains(0) or not pair.imag.contains(0), f"counterexample lost at {mode}")

    decision = artifact["decision"]
    require(decision["direct_discrete_midpoint_roster_simplifies"] is True, "direct parity drift")
    require(decision["canonical_continuous_interpolation_is_linear"] is False, "linear interpolation overclaim")
    require(decision["all_nonzero_symmetric_Poisson_pairs_cancel_at_midpoint"] is False, "pair cancellation overclaim")
    require(decision["prior_linear_interpolation_anchor_rejected"] is True, "rejection decision drift")
    require(decision["complete_T_upper_proved"] is False, "T_upper overclaim")
    print("checked midpoint interpolation guard at 130 digits using direct erfc Fresnel values", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
