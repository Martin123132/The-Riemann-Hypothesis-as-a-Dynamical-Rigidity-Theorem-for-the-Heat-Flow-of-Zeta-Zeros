#!/usr/bin/env python3
"""Independently check the finite-Poisson portcullis saddle reduction."""

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


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_finite_poisson_portcullis_saddle_reduction_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_finite_poisson_portcullis_saddle_reduction_gate.md"
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


def inverse_roots(alpha_value: int) -> tuple[arb, arb]:
    alpha = arb(alpha_value)
    disc = (1 - 8 * arb(T) / (arb.pi() * alpha**2)).sqrt()
    return alpha * (1 - disc) / 4, alpha * (1 + disc) / 4


def main() -> None:
    artifact = load(RESULT)
    require(artifact["passed"] is True, "artifact is not passed")
    require(artifact["decision"]["equation124_is_joint_poisson_kummer_saddle_map"] is True, "equation-(124) decision drift")
    require(artifact["decision"]["projected_alpha_turning_point_is_full_saddle_degeneracy"] is False, "turning guard drift")
    require(artifact["decision"]["reflected_involution_preserves_integer_lattice"] is False, "lattice guard drift")
    require(artifact["decision"]["height_uniform_remainder_proved"] is False, "uniform theorem overpromotion")

    alpha, x, m, t = sp.symbols("alpha x m t", positive=True)
    phase = sp.pi * alpha**2 * x / 4 - sp.pi * m * alpha + t * sp.log((1 - x) / x) / 2
    x_star = 2 * sp.pi * m**2 / (t + 2 * sp.pi * m**2)
    alpha_star = 2 * m + t / (sp.pi * m)
    require(sp.simplify(sp.diff(phase, alpha).subs({alpha: alpha_star, x: x_star})) == 0, "independent alpha saddle failed")
    require(sp.simplify(sp.diff(phase, x).subs({alpha: alpha_star, x: x_star})) == 0, "independent x saddle failed")
    det = sp.factor(sp.simplify(sp.hessian(phase, (alpha, x)).det().subs({alpha: alpha_star, x: x_star})))
    expected = -(2 * sp.pi * m**2 + t) ** 3 / (8 * m**2 * t)
    require(sp.simplify(det - expected) == 0, "independent Hessian failed")
    require(sp.simplify(alpha_star - 2 * m - t / (sp.pi * m)) == 0, "independent equation-(124) failed")
    partner = t / (2 * sp.pi * m)
    require(sp.simplify(x_star.subs(m, partner) - (1 - x_star)) == 0, "independent reflection failed")

    phi = x * alpha**2 / 4 - m * alpha
    exponential = sp.exp(sp.I * sp.pi * phi)
    decomposition = sp.simplify(
        alpha * exponential
        - (2 / x) * (sp.diff(exponential, alpha) / (sp.I * sp.pi) + m * exponential)
    )
    require(decomposition == 0, "independent mode decomposition failed")

    ctx.dps = 130
    lower_a, upper_a = inverse_roots(A)
    lower_b, upper_b = inverse_roots(B)
    nt = (arb(T) / (2 * arb.pi())).sqrt()
    require(arb(39852) < lower_a < arb(39853), "independent lower-A classification failed")
    require(arb(39936) < upper_a < arb(39937), "independent upper-A classification failed")
    require(arb(621) < lower_b < arb(622), "independent lower-B classification failed")
    require(arb(2560588) < upper_b < arb(2560589), "independent upper-B classification failed")
    require(arb(39894) < nt < arb(39895), "independent Nt classification failed")
    saved = artifact["stationary_classification"]
    for field, value in (
        ("lower_root_at_first_alpha_ball", lower_a),
        ("upper_root_at_first_alpha_ball", upper_a),
        ("lower_root_at_last_alpha_ball", lower_b),
        ("upper_root_at_last_alpha_ball", upper_b),
        ("riemann_siegel_nt_ball", nt),
    ):
        require(arb(saved[field]).overlaps(value), f"saved interval drift: {field}")
    require(saved["lower_branch_interior_count"] == 39231, "lower interior count drift")
    require(saved["lower_turning_endpoint_count"] == 42, "lower turning count drift")
    require(saved["upper_turning_endpoint_count"] == 42, "upper turning count drift")
    require(saved["turning_gap_count"] == 84, "turning gap drift")
    require(saved["classical_interior_term_count"] + saved["classical_lower_endpoint_fresnel_term_count"] == saved["classical_source_upper_term_count"] == 39273, "classical roster closure failed")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {dependency['path']}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {dependency['path']}")
    for source in artifact["sources"].values():
        path = REPO_ROOT / source["path"]
        require(file_hash(path) == source["sha256"], f"source hash drift: {source['path']}")

    note = NOTE.read_text(encoding="utf-8")
    for token in (
        "Equation (P1) is exactly paper equation (124)",
        "global two-variable saddle itself does not degenerate",
        "exactly 42 on each side",
        "the integer mode lattice in general",
        "interchange theorem for the Kummer endpoint singularities",
    ):
        require(token in note, f"note token missing: {token}")

    print(
        "validated finite-Poisson portcullis saddle reduction independently: "
        "interior=39231, turning=84, endpoint-classical=42"
    )


if __name__ == "__main__":
    main()
