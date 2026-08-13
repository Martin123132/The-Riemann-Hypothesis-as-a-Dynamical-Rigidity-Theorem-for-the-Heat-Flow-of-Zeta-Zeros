#!/usr/bin/env python3
"""Certify a uniform discrete tangential/normal gradient gap on the B face."""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_discrete_joint_gradient_gap_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "face_trace": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_face_trace_lattice_deghosting_gate.json",
    "B_window": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_gaussian_window_fresnel_remainder_gate.json",
    "step_tail": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_positive_B_crossing_step_tail_normal_form_gate.json",
}

PRECISION = 100
T = 10_000_000_000
B = 5_122_421
FLOOR_MODE = 621
CEIL_MODE = 622


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


def symbolic_certificate() -> dict[str, str]:
    x, endpoint, mode, height, hessian = sp.symbols("x B m t H", positive=True, real=True)
    trace_derivative = sp.pi * endpoint**2 / 4 - height / (2 * x * (1 - x))
    normal_coordinate = (endpoint * x - 2 * mode) / sp.sqrt(2 * x)
    g = trace_derivative / sp.sqrt(hessian)
    require(
        sp.simplify(sp.diff(trace_derivative, x) - height * (1 - 2 * x) / (2 * x**2 * (1 - x) ** 2)) == 0,
        "trace monotonicity derivative failed",
    )
    require(
        sp.simplify(sp.diff(normal_coordinate, x) - (endpoint * x + 2 * mode) / (2 * sp.sqrt(2) * x ** sp.Rational(3, 2))) == 0,
        "normal-coordinate monotonicity derivative failed",
    )

    return {
        "tangential_coordinate": "G(x)=Phi_B'(x)/sqrt(H_B)",
        "normal_coordinate": "q_m(x)=(B*x-2m)/sqrt(2x)",
        "monotonicity": "G'(x)>0 and q_m'(x)>0 on 0<x<1/2",
        "trace_root": "G(x_B^-)=0",
        "normal_root": "q_m(2m/B)=0",
        "left_balance": "For m<=floor(m_B^-), min_x max(|G|,|q_m|) is the unique solution -G=q_m between 2m/B and x_B^-.",
        "right_balance": "For m>=ceil(m_B^-), min_x max(|G|,|q_m|) is the unique solution G=-q_m between x_B^- and 2m/B, or before the half boundary if the crossing lies outside.",
        "lattice_monotonicity": "The left minimax margin increases as m decreases; the right margin increases as m increases. Hence only floor/ceil modes need comparison.",
        "smooth_partition": "chi_T=G^2/(G^2+q_m^2), chi_N=q_m^2/(G^2+q_m^2); the certified gap bounds every denominator away from zero.",
        "scope_guard": "This is a face-residual gradient gap after the sharp interior saddle allocation, not a claim that the full triangle has no interior saddle.",
        "unused_symbolic_G": sp.sstr(g),
    }


def interval_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    t = arb(T)
    endpoint = arb(B)
    pi = arb.pi()
    x0 = (1 - (1 - 8 * t / (pi * endpoint**2)).sqrt()) / 2
    hessian = t * (1 - 2 * x0) / (2 * x0**2 * (1 - x0) ** 2)
    root_h = hessian.sqrt()

    def g(x: arb) -> arb:
        return (pi * endpoint**2 / 4 - t / (2 * x * (1 - x))) / root_h

    def q(mode: int, x: arb) -> arb:
        return (endpoint * x - 2 * arb(mode)) / (2 * x).sqrt()

    left_low = arb("0.00024258422528363")
    left_high = arb("0.00024258422528364")
    left_h_low = -g(left_low) - q(FLOOR_MODE, left_low)
    left_h_high = -g(left_high) - q(FLOOR_MODE, left_high)
    require(left_h_low > 0 and left_h_high < 0, "mode-621 balance bracket failed")
    left_margin_lower = -g(left_high)
    require(left_margin_lower > arb("28.08"), "mode-621 minimax margin below 28.08")

    right_low = arb("0.00024275750696095")
    right_high = arb("0.00024275750696097")
    right_h_low = g(right_low) + q(CEIL_MODE, right_low)
    right_h_high = g(right_high) + q(CEIL_MODE, right_high)
    require(right_h_low < 0 and right_h_high > 0, "mode-622 balance bracket failed")
    right_margin_lower = g(right_low)
    require(right_margin_lower > arb("22.40"), "mode-622 minimax margin below 22.40")
    require(left_margin_lower > right_margin_lower, "nearest-side margin ordering failed")

    left_xi = ((left_low + left_high) / 2 - x0) * root_h
    right_xi = ((right_low + right_high) / 2 - x0) * root_h
    require(-arb("28.08") < left_xi < -arb("28.06"), "mode-621 balance xi drift")
    require(arb("22.41") < right_xi < arb("22.43"), "mode-622 balance xi drift")

    gap_square_lower = arb("22.40") ** 2
    require(gap_square_lower > arb(501), "partition denominator square below 501")

    return {
        "height": T,
        "endpoint": B,
        "trace_root_ball": x0.str(PRECISION, more=True),
        "trace_hessian_ball": hessian.str(PRECISION, more=True),
        "floor_mode": FLOOR_MODE,
        "ceil_mode": CEIL_MODE,
        "mode_621_balance_x_bracket": [str(left_low), str(left_high)],
        "mode_621_balance_sign_balls": [left_h_low.str(PRECISION, more=True), left_h_high.str(PRECISION, more=True)],
        "mode_621_balance_xi_ball": left_xi.str(PRECISION, more=True),
        "mode_621_minimax_margin_lower_ball": left_margin_lower.str(PRECISION, more=True),
        "mode_622_balance_x_bracket": [str(right_low), str(right_high)],
        "mode_622_balance_sign_balls": [right_h_low.str(PRECISION, more=True), right_h_high.str(PRECISION, more=True)],
        "mode_622_balance_xi_ball": right_xi.str(PRECISION, more=True),
        "mode_622_minimax_margin_lower_ball": right_margin_lower.str(PRECISION, more=True),
        "all_positive_integer_face_minimax_margin_lower_bound": "22.40",
        "partition_denominator_G2_plus_q2_lower_bound": "501",
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    return f"""# Discrete joint-gradient gap on the B face

Date: 2026-08-13

Status: exact monotonic reduction plus saved-height interval certificate; not
a proof of the complete B-face estimate

After the sharp interior-saddle allocation, measure the two face directions
by

```text
G(x)=Phi_B'(x)/sqrt(H_B),
q_m(x)=(B*x-2m)/sqrt(2x).                             (DG1)
```

On `0<x<1/2`, both are strictly increasing:

```text
G'(x)=t(1-2x)/[2sqrt(H_B)x^2(1-x)^2]>0,
q_m'(x)=(Bx+2m)/[2sqrt(2)x^(3/2)]>0.                 (DG2)
```

The tangential root is `x_B^-`; the normal root is `2m/B`.  If a crossing
lies to the left of the trace root, the minimum of
`max(|G|,|q_m|)` is the unique balance `-G=q_m` between them.  If it lies to
the right, the balance is `G=-q_m`.  Monotonicity in the integer mode proves
that only the adjacent lattice modes 621 and 622 can minimize the gap.

For mode 621, Arb brackets the balance at

```text
x in [{c['mode_621_balance_x_bracket'][0]},
      {c['mode_621_balance_x_bracket'][1]}],
xi={c['mode_621_balance_xi_ball']},
min max(|G|,|q_621|)>{c['mode_621_minimax_margin_lower_ball']}>28.08. (DG3)
```

For mode 622,

```text
x in [{c['mode_622_balance_x_bracket'][0]},
      {c['mode_622_balance_x_bracket'][1]}],
xi={c['mode_622_balance_xi_ball']},
min max(|G|,|q_622|)>{c['mode_622_minimax_margin_lower_ball']}>22.40. (DG4)
```

Consequently, for every positive integer mode and every B-face point,

```text
max(|G(x)|,|q_m(x)|)>22.40,
G(x)^2+q_m(x)^2>501.                                 (DG5)
```

This is the discrete gain missed by the continuum tangent model.  The B
trace saddle lies between two integer crossings, but in the natural
Gaussian/Fresnel units it is more than 22 units from simultaneous
criticality.  One may therefore use the smooth exact partition

```text
chi_T=G^2/(G^2+q_m^2),
chi_N=q_m^2/(G^2+q_m^2)                              (DG6)
```

and integrate tangentially where `G` dominates and normally where `q_m`
dominates.  No hard crossing collar or singular cotangent split is required
at the theorem level.

The scope is important: (DG5) is a face-residual statement after allocation
of the interior saddle.  For mode 622 the complete triangle still has its
ordinary interior saddle; this certificate does not erase or double-count
that target contribution.

Pi provenance: `pi` comes from the exact equation-(9) face phase and its
Fresnel normal coordinate.  No fitted constant is used.

Proof boundary: exact discrete face noncriticality and a denominator margin
for a two-direction partition only.  No derivative bounds for the partition
or amplitude, completed integration-by-parts estimate, A-fold splice,
complete paired residual, `T_upper`, `Lambda<=0`, PF-infinity, RH, or a
prize-level conclusion is established.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(
        dependencies["face_trace"]["decision"]["B_face_local_integer_roster_is_621_622"] is True,
        "face-roster dependency drift",
    )
    require(
        dependencies["B_window"]["decision"]["xi_70_window_contains_exactly_B_crossings_621_622"] is True,
        "B-window dependency drift",
    )
    require(
        dependencies["step_tail"]["decision"]["tail_and_step_jumps_cancel_exactly"] is True,
        "step-tail continuity dependency drift",
    )

    artifact = {
        "kind": STEM,
        "status": "B_face_discrete_tangential_normal_minimax_gap_above_22p40_certified",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(),
        "decision": {
            "B_face_tangential_and_normal_coordinates_strictly_increasing": True,
            "adjacent_integer_modes_control_global_discrete_face_gap": True,
            "all_positive_integer_face_minimax_margin_above_22p40": True,
            "smooth_two_direction_partition_denominator_above_501": True,
            "continuum_tangent_degeneracy_is_hit_by_an_integer_mode": False,
            "full_triangle_has_no_interior_saddle": False,
            "two_direction_IBP_constant_completed": False,
            "complete_B_face_estimate_proved": False,
            "complete_T_upper_proved": False,
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
            "flint_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Apply the smooth chi_T/chi_N partition to the exact grouped B face residual. Bound its first two derivatives using G^2+q^2>501, integrate the chi_T channel in the trace direction and the chi_N channel in the Fresnel-normal direction, and retain the interior mode-622 saddle as the already allocated target term.",
        "proof_boundary": "Exact discrete B-face gradient gap only. No partition derivative/amplitude constants, complete B estimate, A-fold splice, complete paired residual, complete T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified discrete B-face joint-gradient gap above 22.40", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
