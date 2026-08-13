#!/usr/bin/env python3
"""Independently check the exact Fresnel endpoint normal form."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fresnel_endpoint_normal_form_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fresnel_endpoint_normal_form_gate.md"
T = 10_000_000_000
A = 159577
B = 5122421


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def z_star(mode: int, endpoint: int) -> arb:
    m = arb(mode)
    t = arb(T)
    pi = arb.pi()
    return m * (pi / (t + 2 * pi * m**2)).sqrt() * (arb(endpoint) - 2 * m - t / (pi * m))


def main() -> None:
    artifact = load(RESULT)
    decision = artifact["decision"]
    require(artifact["passed"] is True, "artifact is not passed")
    require(decision["exact_alpha_fresnel_normal_form_proved"] is True, "normal-form decision drift")
    require(decision["all_84_A_endpoint_modes_share_one_compact_chart"] is True, "compact-chart decision drift")
    require(decision["alpha_integral_Airy_normal_form_required"] is False, "alpha Airy guard drift")
    require(decision["composite_x_endpoint_fold_uniformization_required"] == "open_near_characteristic", "composite fold guard drift")
    require(decision["reduced_x_saddle_nondegenerate"] is True, "x-saddle decision drift")
    require(decision["x_stationary_remainder_proved"] is False, "x-remainder overpromotion")

    alpha, x, m, t, C, z = sp.symbols("alpha x m t C z", positive=True)
    z_alpha = sp.sqrt(x / 2) * (alpha - 2 * m / x)
    require(sp.simplify(x * (alpha - 2 * m / x) ** 2 / 4 - z_alpha**2 / 2) == 0, "independent Fresnel scaling failed")
    k = sp.exp(-sp.I * sp.pi / 4) * sp.sqrt(sp.pi / 2)
    primitive = sp.exp(sp.I * sp.pi / 4) * sp.erf(k * z) / sp.sqrt(2)
    require(sp.simplify(sp.diff(primitive, z) - sp.exp(sp.I * sp.pi * z**2 / 2)) == 0, "independent Fresnel primitive failed")
    x_star = 2 * sp.pi * m**2 / (t + 2 * sp.pi * m**2)
    alpha_star = 2 * m + t / (sp.pi * m)
    zc = sp.sqrt(x / 2) * (C - 2 * m / x)
    expected_z = m * sp.sqrt(sp.pi / (t + 2 * sp.pi * m**2)) * (C - alpha_star)
    require(sp.simplify(zc.subs(x, x_star) - expected_z) == 0, "independent endpoint coordinate failed")
    psi = -sp.pi * m**2 / x + t * sp.log((1 - x) / x) / 2
    expected_curvature = -(t + 2 * sp.pi * m**2) ** 4 / (8 * sp.pi**2 * m**4 * t)
    require(sp.simplify(sp.diff(psi, x).subs(x, x_star)) == 0, "independent reduced saddle failed")
    require(sp.simplify(sp.diff(psi, x, 2).subs(x, x_star) - expected_curvature) == 0, "independent reduced curvature failed")

    ctx.dps = 130
    rows = [(mode, z_star(mode, A), z_star(mode, B)) for mode in range(39853, 39937)]
    min_a = min(rows, key=lambda row: row[1])
    max_a = max(rows, key=lambda row: row[1])
    min_b = min(rows, key=lambda row: row[2])
    max_b = max(rows, key=lambda row: row[2])
    error = max_a[1] + 2 / (arb.pi() * min_b[2])
    saved = artifact["certified_atlas"]
    require(saved["A_transition_mode_count"] == len(rows) == 84, "transition count drift")
    require(saved["zA_transition_min_mode"] == min_a[0] == 39936, "zA minimum mode drift")
    require(saved["zA_transition_max_mode"] == max_a[0] == 39894, "zA maximum mode drift")
    require(arb(saved["zA_transition_min_ball"]).overlaps(min_a[1]), "zA minimum ball drift")
    require(arb(saved["zA_transition_max_ball"]).overlaps(max_a[1]), "zA maximum ball drift")
    require(arb(saved["zB_transition_min_ball"]).overlaps(min_b[2]), "zB minimum ball drift")
    require(arb(saved["zB_transition_max_ball"]).overlaps(max_b[2]), "zB maximum ball drift")
    require(arb(saved["half_saddle_factor_error_bound_ball"]).overlaps(error), "half-saddle error ball drift")
    require(arb(0) < min_a[1] <= max_a[1] < arb("0.044"), "compact A chart failed")
    require(min_b[2] > arb("2480000") and error < arb("0.044"), "far-endpoint chart failed")

    neighbors = {
        "zA_mode_39852": z_star(39852, A),
        "zA_mode_39937": z_star(39937, A),
        "zB_mode_621": z_star(621, B),
        "zB_mode_622": z_star(622, B),
        "zB_mode_2560588": z_star(2560588, B),
        "zB_mode_2560589": z_star(2560589, B),
    }
    for name, value in neighbors.items():
        require(arb(saved["key_neighbor_balls"][name]).overlaps(value), f"neighbor ball drift: {name}")
    require(neighbors["zA_mode_39852"] < 0 and neighbors["zA_mode_39937"] < 0, "A-neighbor signs failed")
    require(neighbors["zB_mode_621"] < 0 < neighbors["zB_mode_622"], "lower B-neighbor signs failed")
    require(neighbors["zB_mode_2560589"] < 0 < neighbors["zB_mode_2560588"], "upper B-neighbor signs failed")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {dependency['path']}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {dependency['path']}")
    for source in artifact["sources"].values():
        path = REPO_ROOT / source["path"]
        require(file_hash(path) == source["sha256"], f"source hash drift: {source['path']}")

    note = NOTE.read_text(encoding="utf-8")
    for token in (
        "No asymptotic expansion has been used",
        "alpha integral is a smooth Fresnel endpoint",
        "characteristic boundary tangent",
        "every one of these",
        "half-saddle Fresnel contributions",
        "six positive ranges",
    ):
        require(token in note, f"note token missing: {token}")

    print("validated Fresnel endpoint normal form independently: modes=84, chart<0.044, ranges=6")


if __name__ == "__main__":
    main()
