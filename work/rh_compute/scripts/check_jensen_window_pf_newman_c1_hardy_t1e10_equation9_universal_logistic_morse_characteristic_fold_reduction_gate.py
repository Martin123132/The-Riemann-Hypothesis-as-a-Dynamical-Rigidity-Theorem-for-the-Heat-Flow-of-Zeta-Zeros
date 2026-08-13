#!/usr/bin/env python3
"""Independently check the universal Morse characteristic-fold reduction."""

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


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_universal_logistic_morse_characteristic_fold_reduction_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_universal_logistic_morse_characteristic_fold_reduction_gate.md"
T = 10_000_000_000
A = 159577


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def kappa_and_defect(mode: int) -> tuple[arb, arb]:
    m = arb(mode)
    t = arb(T)
    pi = arb.pi()
    kappa = -t.sqrt() * (pi * arb(A) * m + 2 * pi * m**2 + t) / (pi.sqrt() * (2 * pi * m**2 + t) ** (arb(3) / 2))
    return kappa, 1 - pi * kappa**2 / 2


def main() -> None:
    artifact = load(RESULT)
    decision = artifact["decision"]
    require(artifact["passed"] is True, "artifact is not passed")
    require(decision["universal_global_morse_phase_proved"] is True, "Morse decision drift")
    require(decision["grouped_boundary_fresnel_current_preserved"] is True, "grouped-current decision drift")
    require(decision["classical_inverse_sqrt_mode_carrier_recovered"] is True, "classical carrier drift")
    require(decision["full_two_variable_hessian_degenerate"] is False, "Hessian guard drift")
    require(decision["composite_boundary_tangent_characteristic"] is True, "characteristic decision drift")
    require(decision["plain_frozen_fresnel_stationary_remainder_uniform"] is False, "uniformity guard drift")
    require(decision["fold_uniform_remainder_proved"] is False, "fold theorem overpromotion")

    m, t, v, C, q = sp.symbols("m t v C q", positive=True)
    pi = sp.pi
    r = t / (2 * pi * m**2)
    x = 1 / (1 + r * v)
    psi = -pi * m**2 / x + t * sp.log(r * v) / 2
    psi_star = sp.simplify(psi.subs(v, 1))
    require(sp.simplify(psi - psi_star + t * (v - 1 - sp.log(v)) / 2) == 0, "independent universal phase failed")
    measure = sp.simplify((-sp.diff(x, v)) * (x * (1 - x)) ** (-sp.Rational(1, 4)))
    expected_measure = r ** sp.Rational(3, 4) * v ** (-sp.Rational(1, 4)) * (1 + r * v) ** (-sp.Rational(3, 2))
    require(sp.simplify(measure - expected_measure) == 0, "independent measure failed")

    alpha_q = 2 * m / x + sp.sqrt(2 / x) * q
    density = sp.simplify(alpha_q * sp.diff(alpha_q, q) / 2)
    require(sp.simplify(density - m * sp.sqrt(2) * x ** (-sp.Rational(3, 2)) - q / x) == 0, "independent grouped current failed")

    q_c = sp.sqrt(x / 2) * (C - 2 * m / x)
    kappa = sp.factor(2 * sp.diff(q_c, v).subs(v, 1) / sp.sqrt(t))
    defect = sp.factor(1 - pi * kappa**2 / 2)
    alpha_star = 2 * m + t / (pi * m)
    expected_boundary = (2 * pi * m**2 - t) / (2 * pi * m**2 + t)
    require(sp.simplify(defect.subs(C, alpha_star) - expected_boundary) == 0, "independent boundary defect failed")
    nt = sp.sqrt(t / (2 * pi))
    a = sp.sqrt(8 * t / pi)
    require(sp.simplify(kappa.subs({m: nt, C: a}) + sp.sqrt(2 / pi)) == 0, "independent coalescent slope failed")
    require(sp.simplify(defect.subs({m: nt, C: a})) == 0, "independent coalescent defect failed")

    x_star = sp.simplify(x.subs(v, 1))
    grouped = 2 * m * x_star ** (-sp.Rational(3, 2)) * sp.exp(sp.I * pi / 4)
    measure_star = r ** sp.Rational(3, 4) * (1 + r) ** (-sp.Rational(3, 2))
    gaussian = 2 * sp.sqrt(pi / t) * sp.exp(-sp.I * pi / 4)
    raw = sp.powsimp(grouped * measure_star * gaussian, force=True)
    expected_raw = 2 ** sp.Rational(5, 4) * (t / pi) ** sp.Rational(1, 4) / sp.sqrt(m)
    require(sp.simplify(sp.powsimp(raw / expected_raw, force=True) - 1) == 0, "independent raw main failed")
    normalization = (pi / (32 * t)) ** sp.Rational(1, 4)
    require(sp.simplify(sp.powsimp(raw * normalization * sp.sqrt(m), force=True) - 1) == 0, "independent paper normalization failed")

    ctx.dps = 130
    rows = [(mode, *kappa_and_defect(mode)) for mode in range(39853, 39937)]
    require(all(row[2] < 0 for row in rows[:42]), "independent lower defect signs failed")
    require(all(row[2] > 0 for row in rows[42:]), "independent upper defect signs failed")
    closest = min(rows, key=lambda row: abs(float(row[2])))
    reciprocal = 1 / abs(closest[2])
    saved = artifact["certified_fold_atlas"]
    require(saved["transition_count"] == len(rows) == 84, "transition count drift")
    require(saved["lower_defect_negative_count"] == saved["upper_defect_positive_count"] == 42, "42+42 split drift")
    require(saved["closest_integer_mode"] == closest[0] == 39894, "closest mode drift")
    require(arb(saved["closest_defect_ball"]).overlaps(closest[2]), "closest defect ball drift")
    require(arb(saved["closest_reciprocal_ball"]).overlaps(reciprocal), "reciprocal ball drift")
    require(abs(closest[2]) < arb("6.27e-6") and reciprocal > arb("159000"), "near-characteristic scale drift")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {dependency['path']}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {dependency['path']}")
    for source in artifact["sources"].values():
        path = REPO_ROOT / source["path"]
        require(file_hash(path) == source["sha256"], f"source hash drift: {source['path']}")

    note = NOTE.read_text(encoding="utf-8")
    for token in (
        "independent of `m`",
        "kept grouped without separating",
        "leaving `1/sqrt(m)`",
        "boundary tangent is characteristic",
        "cannot be the desired uniform theorem",
    ):
        require(token in note, f"note token missing: {token}")

    print("validated universal logistic Morse fold reduction independently: modes=84, split=42+42, closest<6.27e-6")


if __name__ == "__main__":
    main()
