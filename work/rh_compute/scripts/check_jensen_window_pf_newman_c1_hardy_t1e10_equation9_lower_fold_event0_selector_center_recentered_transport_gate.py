#!/usr/bin/env python3
"""Independently check selector-centered event-zero transport."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))

import flint
from flint import acb, arb


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event0_selector_center_recentered_transport_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
HEIGHT_BUILDER = REPO_ROOT / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_second_order_all_event_continuous_height_gate.py"
MODE = 39_894


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def main() -> None:
    require(RESULT.is_file(), "missing recentered event-zero artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "recentered artifact did not pass")
    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file() and file_hash(path) == record["sha256"], f"hash drift: {path}")

    height = load_module(HEIGHT_BUILDER, "event0_recentered_check_height")
    atlas, second, first, module, face_class, _ = height.load_sources("event0_recentered_check")
    module.SERIES_ORDER = 82
    module.INNER_ORDER = 112
    module.QUADRATIC_ORDER = 34
    flint.ctx.cap = module.SERIES_ORDER + 8
    flint.ctx.dps = 100
    face = first.event_face(face_class, module, 0, MODE)
    face.t = face.tstar
    face.eta = arb(1)
    face.lam = arb(0)
    order = 9
    moments, mass, _ = height.certified_moments(atlas, module, face, order)
    kappa = face.pi / (16 * face.beta)
    coefficients = [
        (face.i * kappa) ** degree / arb(math.factorial(degree)) * moment
        for degree, moment in enumerate(moments)
    ]
    zmax = (module.HORIZONTAL_END**2 + module.CONTOUR_HEIGHT**2).sqrt()
    size = kappa * zmax
    remainder = mass * size.exp() * size ** (order + 1) / arb(math.factorial(order + 1))
    far = kappa.exp() * second.second_order_far_bound(first, module, face)
    total = sum((abs(value).upper() for value in coefficients), arb(0)) + remainder + far
    require(total < arb("1e-8"), f"independent recentered bound exceeds 1e-8: {total}")
    saved = arb(artifact["certificate"]["uniform_normalized_exact_minus_beta4_bound_ball"])
    require(saved < arb("1e-8"), "saved recentered bound exceeds 1e-8")
    direct_center = height.integrate_height_event(atlas, second, first, module, face_class, 0, 2)
    direct_lower = height.integrate_height_event(atlas, second, first, module, face_class, 0, 1)
    saved_center = parse_complex(artifact["certificate"]["selector_center_ball"])
    saved_lower = parse_complex(artifact["certificate"]["top_corridor_lower_face_ball"])
    require(direct_center.real.overlaps(saved_center.real) and direct_center.imag.overlaps(saved_center.imag), "independent center overlap failed")
    require(direct_lower.real.overlaps(saved_lower.real) and direct_lower.imag.overlaps(saved_lower.imag), "independent lower overlap failed")
    require(artifact["decision"]["remaining_398_mode_finite_integral_splice_proved"] is False, "remaining modes overpromoted")
    print("independently validated selector-centered event-zero transport", flush=True)


if __name__ == "__main__":
    main()
