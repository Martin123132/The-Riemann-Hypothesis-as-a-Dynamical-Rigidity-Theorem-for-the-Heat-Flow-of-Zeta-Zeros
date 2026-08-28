#!/usr/bin/env python3
"""Certify the normal/tangential geometry created by the post-A edge shift."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_normal_tangential_split_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "extended_notch": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_normal_form_gate.json",
    "extended_pole_free": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_pole_free_kernel_gate.json",
    "stationary_geometry": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_separable_triangle_geometry_gate.json",
}

A = 159_577
T = 10_000_000_000
OLD_END = 39_894
NEW_END = 39_936
OLD_EDGE = Fraction(2 * OLD_END + 1, 2)
NEW_EDGE = Fraction(2 * NEW_END + 1, 2)
SHIFT_START = OLD_END + 1
SHIFT_END = NEW_END
PRECISION = 100


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def translation_certificate() -> dict[str, Any]:
    require(NEW_EDGE - OLD_EDGE == 42, "edge shift drift")
    require(SHIFT_END - SHIFT_START + 1 == 42, "shifted finite block drift")
    old = set(range(622, OLD_END + 1))
    new = set(range(622, NEW_END + 1))
    shifted = set(range(SHIFT_START, SHIFT_END + 1))
    require(new == old | shifted and old.isdisjoint(shifted), "edge translation set identity failed")

    return {
        "old_target": [622, OLD_END],
        "new_target": [622, NEW_END],
        "shifted_block": [SHIFT_START, SHIFT_END],
        "shifted_mode_count": len(shifted),
        "edge_shift": str(NEW_EDGE - OLD_EDGE),
        "fixed_regulator_identity": "C_U,epsilon(s)-C_T,epsilon(s)=-sum_(m=39895)^39936 exp(-pi*epsilon*m^2)exp(-2*pi*i*m*s)",
        "interpolant_identity": "h_U,epsilon(y)-h_T,epsilon(y)=-q_epsilon(y)1_(39894.5<y<39936.5), up to endpoint conventions of measure zero",
        "dual_identity": "sum_(k in Z)[hat h_U,epsilon(k+s)-hat h_T,epsilon(k+s)]=-sum_(m=39895)^39936 w_m exp(-2*pi*i*m*s)",
        "zero_regulator_polynomial": "C_U,0(s)-C_T,0(s)=-exp(-pi*i*79831*s)sin(42*pi*s)/sin(pi*s), with continuous value -42 at s in Z",
        "ownership_guard": "The improved upper-edge geometry is exactly balanced by the explicit 42-mode block 39895..39936; it is not an independent small remainder.",
    }


def phase_geometry_certificate() -> dict[str, Any]:
    alpha, x, m = sp.symbols("alpha x m", positive=True)
    phase = sp.pi * alpha**2 * x / 4 - sp.pi * m * alpha + sp.Rational(T, 2) * sp.log((1 - x) / x)
    normal = sp.simplify(sp.diff(phase, alpha))
    tangential = sp.simplify(sp.diff(phase.subs(alpha, A), x))
    require(sp.simplify(normal - sp.pi * (alpha * x / 2 - m)) == 0, "normal derivative drift")
    require(sp.simplify(tangential - (sp.pi * A**2 / 4 - sp.Rational(T, 2) / (x * (1 - x)))) == 0, "tangential derivative drift")

    old_corner_gap = OLD_EDGE - Fraction(A, 4)
    new_corner_gap = NEW_EDGE - Fraction(A, 4)
    require(old_corner_gap == Fraction(1, 4), "old corner gap drift")
    require(new_corner_gap == Fraction(169, 4), "new corner gap drift")
    require(new_corner_gap / old_corner_gap == 169, "normal-gap gain drift")

    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    discriminant = 1 - arb(8 * T) / (pi * arb(A) ** 2)
    x_star = (1 - discriminant.sqrt()) / 2
    old_face_gap = arb(OLD_EDGE.numerator) / OLD_EDGE.denominator - arb(A) * x_star / 2
    new_face_gap = arb(NEW_EDGE.numerator) / NEW_EDGE.denominator - arb(A) * x_star / 2
    require(old_face_gap.lower() > arb(42), "old face-saddle normal gap drift")
    require(new_face_gap.lower() > arb(84), "new face-saddle normal gap drift")
    require((new_face_gap - old_face_gap).contains(42), "face-saddle gap shift drift")

    curvature = arb(T) * (1 - 2 * x_star) / (2 * (x_star * (1 - x_star)) ** 2)
    require(curvature.lower() > arb(0), "tangential Morse curvature drift")

    return {
        "joint_phase": "Phi(alpha,x;m)=pi*alpha^2*x/4-pi*m*alpha+(t/2)log((1-x)/x)",
        "A_face_normal_derivative": "partial_alpha Phi(A,x;m)=pi*(A*x/2-m)",
        "A_face_tangential_derivative": "partial_x Phi(A,x;m)=pi*A^2/4-t/[2*x*(1-x)]",
        "uniform_half_domain_normal_gaps": {
            "old_edge": "|partial_alpha Phi(A,x;39894.5)|>=pi/4 for 0<=x<=1/2",
            "new_edge": "|partial_alpha Phi(A,x;39936.5)|>=169*pi/4 for 0<=x<=1/2",
            "exact_gain_factor": 169,
        },
        "corner_normal_coordinates": {
            "old_q_A": "-1/2",
            "new_q_A": "-169/2",
        },
        "tangential_stationary_point": "x_*=[1-sqrt(1-8*t/(pi*A^2))]/2",
        "x_star_ball": x_star.str(80, more=True),
        "old_edge_gap_at_x_star_ball": old_face_gap.str(80, more=True),
        "new_edge_gap_at_x_star_ball": new_face_gap.str(80, more=True),
        "tangential_curvature_ball": curvature.str(80, more=True),
        "classification": "The new upper edge is uniformly nonstationary in the A-face normal direction but retains a nondegenerate tangential Morse stationary point.",
        "licensed_order": "Perform the signed normal endpoint reduction first, then a tangential Morse/Fresnel estimate; do not integrate by parts through x_*.",
    }


def render_note(artifact: dict[str, Any]) -> str:
    geometry = artifact["phase_geometry_certificate"]
    return f"""# Post-A extended-notch normal/tangential split

Date: 2026-08-23

Status: exact edge-translation identity and mixed endpoint/Morse geometry
certified; quantitative signed assembly remains open

Let `T={{622,...,39894}}`, `U={{622,...,39936}}`, and retain the common
Gaussian regulator.  Since `U minus T={{39895,...,39936}}`, exactly

```text
C_U,epsilon(s)-C_T,epsilon(s)
 =-sum_(m=39895)^39936
    exp(-pi*epsilon*m^2)exp(-2*pi*i*m*s).             (NT1)
```

On the continuous interpolants this is the 42-mode edge translation

```text
h_U,epsilon(y)-h_T,epsilon(y)
 =-q_epsilon(y)1_(39894.5<y<39936.5),                (NT2)
```

up to measure-zero endpoint conventions.  Poisson summation therefore gives

```text
sum_(k in Z)[hat h_U,epsilon(k+s)-hat h_T,epsilon(k+s)]
 =-sum_(m=39895)^39936 w_m exp(-2*pi*i*m*s).          (NT3)
```

At zero regulator the right side is the finite polynomial

```text
-exp(-pi*i*79831*s)sin(42*pi*s)/sin(pi*s),            (NT4)
```

with continuous value `-42` at integer `s`.  Thus the apparent analytical
gain from moving the upper edge is exactly balanced by the explicit 42-mode
outside-target half of the A window.  The edge and that block must stay in
one signed assembly.

For the exact joint phase

```text
Phi(alpha,x;m)
 =pi*alpha^2*x/4-pi*m*alpha+(t/2)log((1-x)/x),

partial_alpha Phi(A,x;m)=pi(A*x/2-m),
partial_x Phi(A,x;m)=pi*A^2/4-t/[2*x(1-x)].           (NT5)
```

On `0<=x<=1/2`, the old and new upper edges obey

```text
|partial_alpha Phi(A,x;39894.5)| >= pi/4,
|partial_alpha Phi(A,x;39936.5)| >= 169*pi/4.         (NT6)
```

The edge shift therefore improves the uniform A-face normal denominator by
the exact factor `169`.  At the tangential stationary point

```text
x_*=[1-sqrt(1-8t/(pi*A^2))]/2
   ={geometry['x_star_ball']},                        (NT7)
```

the mode-unit normal gaps are

```text
old: {geometry['old_edge_gap_at_x_star_ball']}
new: {geometry['new_edge_gap_at_x_star_ball']}.       (NT8)
```

But `partial_x Phi(A,x_*;m)=0` for every `m`, and the positive tangential
curvature is `{geometry['tangential_curvature_ball']}`.  The new edge is
therefore uniformly nonstationary in the normal direction, not in the
tangential direction.  The licensed order is a signed normal endpoint
reduction followed by a tangential Morse/Fresnel estimate.  Integrating by
parts through `x_*` is invalid.

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Pi provenance: every `pi` in (NT1)--(NT8) is inherited from the Gaussian
Abel/Fourier kernel and the exact Kummer phase.  No fitted or geometric
occurrence is introduced.

Proof boundary: exact fixed-regulator edge translation, finite-block
ownership, uniform normal gap, and tangential Morse classification only.
No endpoint-expansion remainder, signed 42-mode cancellation bound,
quantitative `R_after_A`, `R_Dir`, or `Q_K-T` estimate, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["extended_notch"]["decision"]["six_row_mask_compressed_coefficientwise"] is True, "extended-notch compression drift")
    require(dependencies["extended_pole_free"]["decision"]["pole_free_closed_centred_cell_formula_certified"] is True, "extended pole-free kernel drift")
    require(dependencies["stationary_geometry"]["decision"]["endpoint_characteristic_is_triangle_face_crossing"] is True, "stationary geometry drift")

    artifact = {
        "kind": STEM,
        "status": "exact_post_A_edge_translation_and_A_face_normal_tangential_geometry_certified_signed_estimate_open",
        "passed": True,
        "scope": {
            "height": T,
            "A": A,
            "old_upper_edge": str(OLD_EDGE),
            "new_upper_edge": str(NEW_EDGE),
            "shifted_modes": [SHIFT_START, SHIFT_END],
        },
        "translation_certificate": translation_certificate(),
        "phase_geometry_certificate": phase_geometry_certificate(),
        "decision": {
            "edge_translation_equals_explicit_42_mode_block": True,
            "zero_regulator_translation_is_finite_trigonometric_polynomial": True,
            "new_upper_edge_uniformly_nonstationary_in_A_face_normal_direction": True,
            "normal_denominator_gain_factor_169_certified": True,
            "tangential_A_Morse_stationary_point_persists": True,
            "normal_then_tangential_analysis_licensed": True,
            "blanket_x_integration_by_parts_licensed": False,
            "signed_42_mode_cancellation_bound_proved": False,
            "R_after_A_bound_proved": False,
            "R_Dir_bound_proved": False,
            "QK_minus_T_bound_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "arb_threads": 1,
            "process_priority": priority,
            "precision_decimal_digits": PRECISION,
        },
        "next_obligation": "Derive the first two signed A-face normal endpoint currents for the translated upper edge and combine them algebraically with the explicit 42-mode block before taking norms. Put the resulting one-dimensional amplitude into the existing exact A-face Morse coordinate and run a bounded floating scout to choose the certified core/tail partition.",
        "proof_boundary": "Exact fixed-regulator edge translation, finite-block ownership, uniform normal gap, and tangential Morse classification only. No endpoint-expansion remainder, signed 42-mode cancellation bound, quantitative R_after_A, R_Dir, or Q_K-T estimate, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified post-A extended-notch normal/tangential split", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
