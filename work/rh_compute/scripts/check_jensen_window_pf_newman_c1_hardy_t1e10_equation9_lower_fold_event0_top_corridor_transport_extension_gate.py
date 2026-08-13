#!/usr/bin/env python3
"""Independently check event-0 transport through the selector top corridor."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, acb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event0_top_corridor_transport_extension_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
PRECISION = 110


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def main() -> None:
    require(RESULT.is_file(), "missing event-0 extension artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "event-0 extension artifact did not pass")
    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file() and file_hash(path) == record["sha256"], f"hash drift: {path}")
    ctx.dps = PRECISION
    atlas_path = REPO_ROOT / artifact["dependencies"]["continuous_event_atlas"]["path"]
    atlas = json.loads(atlas_path.read_text(encoding="utf-8"))
    row = atlas["certificate"]["all_rows"][0]
    coefficients = [parse_complex(record["transport_coefficient_ball"]) for record in row["moments"]]
    mass = arb(row["finite_contour_absolute_mass_bound"])
    x = arb(row["maximum_transport_argument"])
    far0 = arb(row["event_far_tail_absolute_bound"])
    zmax = (arb(21) * 21 + 1).sqrt()
    polynomial = sum((arb(2) ** degree * abs(value).upper() for degree, value in enumerate(coefficients)), arb(0))
    remainder = mass * (2 * x).exp() * (2 * x) ** 7 / arb(math.factorial(7))
    far = (2 * x / zmax).exp() * far0
    total = polynomial + remainder + far
    physical = arb(2).sqrt() / arb.pi() * total
    require(total < arb("2.15e-6"), "independent theta<=2 bound failed")
    require(total <= arb(artifact["certificate"]["uniform_normalized_exact_minus_beta4_bound_ball"]).upper(), "independent bound exceeds saved enclosure")
    require(physical <= arb(artifact["certificate"]["uniform_physical_exact_minus_beta4_bound_ball"]).upper(), "independent physical bound exceeds saved enclosure")
    pi = arb.pi()
    c = arb(159_577)
    m = arb(39_894)
    tstar = pi * c * c / 8
    tau = pi * m * (c - 2 * m)
    require((tau + 2 * pi / 16 - tstar).contains(0), "independent theta endpoint identity failed")
    require(artifact["decision"]["remaining_398_mode_grouped_estimate_proved"] is False, "remaining modes overpromoted")
    print("independently validated event-0 theta<=2 transport extension", flush=True)


if __name__ == "__main__":
    main()
