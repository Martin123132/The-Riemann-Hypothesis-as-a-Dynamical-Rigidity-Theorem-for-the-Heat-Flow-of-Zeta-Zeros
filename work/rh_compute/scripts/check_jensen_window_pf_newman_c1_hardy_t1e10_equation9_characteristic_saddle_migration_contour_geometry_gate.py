#!/usr/bin/env python3
"""Independently check the characteristic contour-geometry gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))

from flint import arb, ctx


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_saddle_migration_contour_geometry_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_saddle_migration_contour_geometry_gate.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    require(RESULT.is_file() and NOTE.is_file(), "missing result or note")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "gate is not passed")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file() and file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file() and file_hash(path) == record["sha256"], f"source hash drift: {path}")

    ctx.dps = 130
    t = arb(10_000_000_000); c = arb(159577); pi = arb.pi()
    eta = pi*c**2/(8*t); beta = (t*eta)**(arb(1)/3); lam = (eta-1)*t/beta
    old = arb(16)-lam; new = arb(81)-lam; maximum = (lam+64).sqrt()
    margin = arb(81)-(lam+64); decay = margin/2
    phase_error = arb(64)/15*arb(9)**5/beta**2 + arb(64)*arb(9)**3/(3*beta**2) + arb(9)*arb(64)**2/(4*beta**2)
    amp = max(4*beta/(pi*c)*arb(64)/c, 1-(arb(9)/beta).cosh()**(-arb(3)/2))
    values = artifact["certified_geometry"]
    checks = {
        "old_radius_saddle_crossing_y_ball": old,
        "new_radius_saddle_crossing_y_ball": new,
        "maximum_saddle_radius_ball": maximum,
        "new_radius_saddle_clearance_ball": arb(9)-maximum,
        "outgoing_ray_quadratic_margin_ball": margin,
        "outgoing_ray_linear_decay_coefficient_ball": decay,
        "expanded_core_total_phase_error_ball": phase_error,
        "expanded_core_relative_amplitude_error_ball": amp,
    }
    for key, value in checks.items():
        require(arb(values[key]).contains(value), f"independent value escaped recorded ball: {key}")
    require(old < 64 < new, "saddle migration classification failed")
    require(maximum < 9 and decay > arb("5.9"), "R9 contour margin failed")
    require(phase_error < arb("0.06") and amp < arb("1.4e-5"), "expanded core admissibility failed")
    require(artifact["decision"]["old_abs_z_gt_4_term_is_pure_nonstationary_tail"] is False, "migration guard drift")
    require(artifact["decision"]["contour_tail_integral_bound_proved"] is False, "proof boundary drift")
    note = NOTE.read_text(encoding="utf-8")
    for token in ("16-lambda", "exp(i*pi/6)", "<0.06", "not a proof"):
        require(token in note, f"note token missing: {token}")
    print("validated contour geometry independently: crossing=10.8901, max saddle<9, ray decay>5.9")


if __name__ == "__main__":
    main()
