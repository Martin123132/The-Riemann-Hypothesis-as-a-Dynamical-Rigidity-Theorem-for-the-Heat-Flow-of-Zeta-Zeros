#!/usr/bin/env python3
"""Certify the exact face trace and isolate its genuine lattice-local modes."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_face_trace_lattice_deghosting_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "triangle": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_separable_triangle_geometry_gate.json",
    "bi_Morse": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_exact_bimorse_face_fold_gate.json",
    "B_tangent": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_kummer_B_crossing_tangent_profile_barrier_gate.json",
    "half_reflection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_kummer_reflection_branch_reduction_gate.json",
}

PRECISION = 100
T = 10_000_000_000
A = 159_577
B = 5_122_421


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
    z, x, m, endpoint, t = sp.symbols("z x m D t", positive=True, real=True)
    phase = (
        sp.pi * m**2 / z
        + sp.pi * endpoint**2 * z / 4
        - sp.pi * m**2 / x
        + t * sp.log((1 - x) / x) / 2
    )
    trace = sp.simplify(phase.subs(z, x))
    expected_trace = sp.pi * endpoint**2 * x / 4 + t * sp.log((1 - x) / x) / 2
    require(sp.simplify(trace - expected_trace) == 0, "mode cancellation on face failed")
    require(sp.simplify(sp.diff(trace, m)) == 0, "face trace retained mode dependence")

    discriminant = sp.sqrt(1 - 8 * t / (sp.pi * endpoint**2))
    x_minus = (1 - discriminant) / 2
    x_plus = (1 + discriminant) / 2
    trace_derivative = sp.diff(trace, x)
    require(sp.simplify(trace_derivative.subs(x, x_minus)) == 0, "lower trace root failed")
    require(sp.simplify(trace_derivative.subs(x, x_plus)) == 0, "upper trace root failed")

    m_minus = sp.simplify(endpoint * x_minus / 2)
    m_plus = sp.simplify(endpoint * x_plus / 2)
    alpha_m = 2 * m + t / (sp.pi * m)
    require(sp.simplify(alpha_m.subs(m, m_minus) - endpoint) == 0, "lower inverse root failed")
    require(sp.simplify(alpha_m.subs(m, m_plus) - endpoint) == 0, "upper inverse root failed")

    normal_derivative = sp.simplify(sp.diff(phase, z).subs({z: x_minus, x: x_minus}))
    expected_normal = sp.pi * endpoint**2 / 4 * (1 - (m / m_minus) ** 2)
    require(sp.simplify(normal_derivative - expected_normal) == 0, "normal detuning failed")

    tangential_hessian = sp.simplify(sp.diff(trace, x, 2).subs(x, x_minus))
    expected_hessian = t * (1 - 2 * x_minus) / (2 * x_minus**2 * (1 - x_minus) ** 2)
    require(sp.simplify(tangential_hessian - expected_hessian) == 0, "trace Hessian failed")

    return {
        "face_trace": "Phi_D(x)=pi*D^2*x/4+(t/2)log((1-x)/x), independent of m",
        "trace_roots": "x_D^pm=[1 pm sqrt(1-8t/(pi D^2))]/2",
        "lattice_centers": "m_D^pm=D*x_D^pm/2=[D pm sqrt(D^2-8t/pi)]/4",
        "inverse_root_identity": "alpha_(m_D^pm)=D for alpha_m=2m+t/(pi*m)",
        "normal_detuning": "partial_z Phi_(m,D)|(z=x=x_D^-)=pi*D^2[1-(m/m_D^-)^2]/4",
        "trace_hessian": "Phi_D''(x_D^-)=t(1-2x_D^-)/[2x_D^-^2(1-x_D^-)^2]>0",
        "genuine_face_criticality": "Both trace and normal derivatives vanish only at the noninteger lattice center m=m_D^-.",
        "tangent_deghosting": "Delta_B near zero at mode 257 is a null tangent at the mode-centered outer saddle, not a joint or face critical point; the exact normal derivative at the B trace saddle is nonzero.",
    }


def endpoint_record(endpoint_value: int) -> dict[str, Any]:
    t = arb(T)
    pi = arb.pi()
    endpoint = arb(endpoint_value)
    discriminant = (1 - 8 * t / (pi * endpoint**2)).sqrt()
    x_minus = (1 - discriminant) / 2
    x_plus = (1 + discriminant) / 2
    m_minus = endpoint * x_minus / 2
    m_plus = endpoint * x_plus / 2
    hessian = t * (1 - 2 * x_minus) / (2 * x_minus**2 * (1 - x_minus) ** 2)

    floor_mode = int(float(m_minus.mid()))
    require(arb(floor_mode) < m_minus < arb(floor_mode + 1), "lattice-center floor is not certified")
    ceil_mode = floor_mode + 1

    def normal(mode: int) -> arb:
        return pi * endpoint**2 / 4 * (1 - (arb(mode) / m_minus) ** 2)

    return {
        "endpoint": endpoint_value,
        "discriminant_ball": discriminant.str(PRECISION, more=True),
        "x_minus_ball": x_minus.str(PRECISION, more=True),
        "x_plus_ball": x_plus.str(PRECISION, more=True),
        "m_minus_ball": m_minus.str(PRECISION, more=True),
        "m_plus_ball": m_plus.str(PRECISION, more=True),
        "adjacent_integer_modes": [floor_mode, ceil_mode],
        "distance_to_floor_ball": (m_minus - floor_mode).str(PRECISION, more=True),
        "distance_to_ceil_ball": (ceil_mode - m_minus).str(PRECISION, more=True),
        "lower_trace_hessian_ball": hessian.str(PRECISION, more=True),
        "normal_derivative_floor_ball": normal(floor_mode).str(PRECISION, more=True),
        "normal_derivative_ceil_ball": normal(ceil_mode).str(PRECISION, more=True),
    }


def interval_ledger() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    b = endpoint_record(B)
    a = endpoint_record(A)
    pi = arb.pi()
    endpoint = arb(B)
    m_b = arb(b["m_minus_ball"])

    def normal_b(mode: int) -> arb:
        return pi * endpoint**2 / 4 * (1 - (arb(mode) / m_b) ** 2)

    x_257 = 2 * pi * arb(257) ** 2 / (arb(T) + 2 * pi * arb(257) ** 2)
    x_b = arb(b["x_minus_ball"])
    require(arb("621.55") < m_b < arb("621.56"), "B lattice center drift")
    require(arb(b["distance_to_floor_ball"]) > arb("0.55"), "B floor margin below 0.55")
    require(arb(b["distance_to_ceil_ball"]) > arb("0.44"), "B ceil margin below 0.44")
    require(abs(normal_b(621)) > arb("3.68e10"), "B mode-621 normal margin too small")
    require(abs(normal_b(622)) > arb("2.94e10"), "B mode-622 normal margin too small")
    require(abs(normal_b(257)) > arb("1.70e13"), "mode-257 normal margin too small")
    require(x_b - x_257 > arb("0.0002"), "mode-257 outer saddle too close to B trace saddle")
    require(arb("39852.39") < arb(a["m_minus_ball"]) < arb("39852.40"), "A lower center drift")
    require(arb("39936.10") < arb(a["m_plus_ball"]) < arb("39936.11"), "A upper center drift")
    return {
        "height": T,
        "B": b,
        "A": a,
        "B_nonadjacent_lattice_distance_lower_bound": "1.44",
        "B_mode_257_normal_derivative_ball": normal_b(257).str(PRECISION, more=True),
        "B_mode_257_outer_saddle_ball": x_257.str(PRECISION, more=True),
        "B_trace_saddle_minus_mode257_outer_saddle_ball": (x_b - x_257).str(PRECISION, more=True),
        "local_chart_ledger": {
            "B_lower_face": [621, 622],
            "A_lower_face": [39_852, 39_853],
            "A_reflected_upper_face": [39_936, 39_937],
            "half_boundary_outer_saddle": [39_894, 39_895],
        },
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_ledger"]
    b, a = c["B"], c["A"]
    return f"""# Face-trace lattice geometry and tangent deghosting

Date: 2026-08-13

Status: interval certificate; not a proof of the B face remainder bound

On the triangular face `z=x`, the two mode phases cancel exactly:

```text
Phi_(m,D)(x,x)
 =pi*D^2*x/4+(t/2)log((1-x)/x)=Phi_D(x).              (FT1)
```

The trace is independent of `m`.  Its two stationary points are

```text
x_D^pm=[1 pm sqrt(1-8t/(pi D^2))]/2,                  (FT2)
```

and their lattice centers are

```text
m_D^pm=D*x_D^pm/2
      =[D pm sqrt(D^2-8t/pi)]/4.                      (FT3)
```

These are exactly the two inverse roots of `alpha_m=D`.  At the lower trace
saddle, the normal derivative is

```text
partial_z Phi_(m,D)
 =pi*D^2[1-(m/m_D^-)^2]/4.                            (FT4)
```

Hence a genuine face critical point requires `m=m_D^-`, not merely a null
tangent in a mode-centered approximation.

For `B=5122421`, rigorous intervals give

```text
x_B^-={b['x_minus_ball']},
m_B^-={b['m_minus_ball']}.                             (FT5)
```

Thus only modes `621` and `622` are adjacent to the B face event.  Their
lattice distances exceed `0.55` and `0.44`, and even their exact normal
derivatives have magnitudes above `3.68e10` and `2.94e10`.  Every other
integer mode is at least `1.44` lattice units from the event.

This resolves the apparent tangent singularity near mode `257`.  Although
the old mode-centered tangent defect is close to zero there, the exact B-face
normal derivative is

```text
{c['B_mode_257_normal_derivative_ball']},              (FT6)
```

with magnitude above `1.70e13`; its outer saddle is also more than `0.0002`
away from `x_B^-`.  Mode `257` is not a joint or face critical point.  Its
large tangent correction came from extending a local line to a remote
stationary region.

For the lower A face,

```text
m_A^-={a['m_minus_ball']},                             (FT7)
```

so modes `39852,39853` are the adjacent lower-face pair.  The reflected root
is `{a['m_plus_ball']}` and is already represented by half-Kummer
conjugation.  The separate `39894,39895` pair belongs to the `x=1/2` outer
saddle corner.

The quantitative theorem can therefore use a two-mode B face chart plus
normal integration by parts for all other B modes.  The old 51-mode tangent
collar remains a diagnostic partition, not the natural exact local roster.

Pi provenance: all `pi` factors come from the exact equation-(9) pair phase
and its boundary trace.  No fitted constant is used.

Proof boundary: exact trace cancellation, stationary roots, and saved-height
lattice/normal margins only.  No normal integration-by-parts constants,
two-mode B face estimate, A-fold splice, complete paired residual,
`T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion follows.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(
        dependencies["triangle"]["decision"]["endpoint_characteristic_is_triangle_face_crossing"] is True,
        "triangle dependency drift",
    )
    require(
        dependencies["bi_Morse"]["decision"]["pair_triangle_phase_is_exactly_hyperbolic_quadratic"] is True,
        "bi-Morse dependency drift",
    )
    require(
        dependencies["B_tangent"]["certificate"]["minimum_integer_abs_Delta_B_mode"] == 257,
        "tangent witness dependency drift",
    )
    require(
        dependencies["half_reflection"]["decision"]["reflected_upper_branch_requires_separate_stationary_estimate"] is False,
        "reflection dependency drift",
    )

    artifact = {
        "kind": STEM,
        "status": "exact_mode_independent_face_trace_and_saved_height_lattice_margins_proved",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_ledger": interval_ledger(),
        "decision": {
            "triangle_face_trace_is_mode_independent": True,
            "face_lattice_centers_equal_portcullis_inverse_roots": True,
            "B_face_local_integer_roster_is_621_622": True,
            "mode_257_is_true_face_critical_point": False,
            "old_tangent_Delta_near_zero_is_exact_residual_singularity": False,
            "B_nonlocal_normal_IBP_bound_proved": False,
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
        "next_obligation": "Derive the normal integration-by-parts operator around the exact B trace saddle for all m outside {621,622}, using the explicit normal detuning in (FT4). Evaluate modes 621 and 622 with the exact curved bi-Morse face chart and complete transformed amplitude, then aggregate before absolute values.",
        "proof_boundary": "Exact face trace and saved-height lattice margins only. No normal IBP constants, two-mode B face estimate, A-fold splice, complete paired residual, complete T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified face trace: B local roster={621,622}, mode257 deghosted", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
