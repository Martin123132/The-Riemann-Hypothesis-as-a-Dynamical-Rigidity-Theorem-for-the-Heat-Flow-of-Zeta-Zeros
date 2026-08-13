#!/usr/bin/env python3
"""Recenter event-zero beta^-4 transport at the selector center."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import flint
from flint import acb, arb


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event0_selector_center_recentered_transport_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
HEIGHT_BUILDER = REPO_ROOT / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_second_order_all_event_continuous_height_gate.py"
DEPENDENCIES = {
    "old_event0_extension": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event0_top_corridor_transport_extension_gate.json",
    "completed_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.json",
}

MODE = 39_894
MOMENT_ORDER = 8
PRECISION = 90


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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


def complex_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": value.real.str(PRECISION, more=True),
        "imag_ball": value.imag.str(PRECISION, more=True),
    }


def configure(module) -> None:
    module.SERIES_ORDER = 70
    module.INNER_ORDER = 96
    module.QUADRATIC_ORDER = 30
    flint.ctx.cap = module.SERIES_ORDER + 8


def calculate() -> dict[str, Any]:
    height = load_module(HEIGHT_BUILDER, "event0_recentered_height")
    atlas, second, first, module, face_class, _ = height.load_sources("event0_recentered")
    configure(module)
    flint.ctx.dps = PRECISION
    face = first.event_face(face_class, module, 0, MODE)
    face.t = face.tstar
    face.eta = arb(1)
    face.lam = arb(0)

    moments, mass, mass_by_piece = height.certified_moments(atlas, module, face, MOMENT_ORDER)
    kappa = face.pi / (16 * face.beta)
    coefficients = [
        (face.i * kappa) ** degree / arb(math.factorial(degree)) * moment
        for degree, moment in enumerate(moments)
    ]
    z_max = (module.HORIZONTAL_END**2 + module.CONTOUR_HEIGHT**2).sqrt()
    size = kappa * z_max
    remainder = mass * size.exp() * size ** (MOMENT_ORDER + 1) / arb(math.factorial(MOMENT_ORDER + 1))
    center_far = second.second_order_far_bound(first, module, face)
    transported_far = kappa.exp() * center_far
    uniform = sum((abs(value).upper() for value in coefficients), arb(0)) + remainder + transported_far
    physical = arb(2).sqrt() / face.pi * uniform
    center = module.add_error(moments[0], center_far)
    lower = module.add_error(height.polynomial_value(coefficients, -1), remainder + transported_far)
    direct_center = height.integrate_height_event(atlas, second, first, module, face_class, 0, 2)
    direct_lower = height.integrate_height_event(atlas, second, first, module, face_class, 0, 1)
    require(center.real.overlaps(direct_center.real) and center.imag.overlaps(direct_center.imag), "center direct overlap failed")
    require(lower.real.overlaps(direct_lower.real) and lower.imag.overlaps(direct_lower.imag), "lower-face direct overlap failed")
    require(uniform < arb("1e-8"), f"recentered top-corridor bound exceeds 1e-8: {uniform}")
    return {
        "transport_center": "t*",
        "transport_interval": "t*-pi/16<=t<=t*",
        "transport_parameter": "h=(pi/16)theta, -1<=theta<=0",
        "moment_order": MOMENT_ORDER,
        "series_order": module.SERIES_ORDER,
        "inner_order": module.INNER_ORDER,
        "quadratic_order": module.QUADRATIC_ORDER,
        "finite_contour_absolute_mass_bound_ball": mass.str(PRECISION, more=True),
        "absolute_mass_by_piece": {name: value.str(PRECISION, more=True) for name, value in sorted(mass_by_piece.items())},
        "maximum_transport_argument_ball": size.str(PRECISION, more=True),
        "finite_exponential_remainder_bound_ball": remainder.str(PRECISION, more=True),
        "center_far_tail_bound_ball": center_far.str(PRECISION, more=True),
        "transported_far_tail_bound_ball": transported_far.str(PRECISION, more=True),
        "uniform_normalized_exact_minus_beta4_bound_ball": uniform.str(PRECISION, more=True),
        "uniform_physical_exact_minus_beta4_bound_ball": physical.str(PRECISION, more=True),
        "selector_center_ball": complex_record(center),
        "top_corridor_lower_face_ball": complex_record(lower),
        "direct_selector_center_ball": complex_record(direct_center),
        "direct_top_corridor_lower_face_ball": complex_record(direct_lower),
        "moments": [
            {"degree": degree, "moment_ball": complex_record(moment), "transport_coefficient_ball": complex_record(coefficients[degree])}
            for degree, moment in enumerate(moments)
        ],
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Event-zero transport recentered at the selector

Date: 2026-08-13

Status: rigorous event-zero exact-minus-beta-four top-corridor enclosure;
not a complete selected-branch or `T_upper` splice

The earlier extension transported seven wide event-centered moment balls over
`0<=theta<=2`.  Its `2.141e-6` majorant was dominated by lost interval
correlation even though the selector-center midpoint was near `10^-13`.

Recenter the exact transport identity at `t*` instead:

```text
D_(0,t*+h)(z)=exp(ihz/beta)D_(0,t*)(z),
h=(pi/16)theta,       -1<=theta<=0.                  (RC1)
```

The `pi/16` width is inherited from the exact selector/event geometry.  A
fresh contour calculation through moment degree `{c['moment_order']}`, rather
than reparsing the old moment balls, gives

```text
sup_(t*-pi/16<=t<=t*) |Delta I_2(0,t)|
 < {c['uniform_normalized_exact_minus_beta4_bound_ball']} <1e-8,
physical < {c['uniform_physical_exact_minus_beta4_bound_ball']}. (RC2)
```

Direct contour integrations at both faces overlap the recentered transport
balls.  This replaces the old top-corridor majorant for event zero; it does
not alter or invalidate the original all-event cells.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No remaining-398-mode finite-integral theorem, completed-remainder bound,
all-corridor continuation, complete `Q_K-T` or `T_upper`, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), HEIGHT_BUILDER, CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    require(dependencies["old_event0_extension"]["decision"]["event0_fold_chart_covers_complete_top_corridor"] is True, "old event-zero dependency drift")
    require(dependencies["completed_projection"].get("passed") is True, "completed projection dependency failed")
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "event_zero_top_corridor_recentered_at_selector_center_with_fresh_moments",
        "passed": True,
        "certificate": calculate(),
        "decision": {
            "old_event_centered_2_141e_minus_6_majorant_superseded_on_top_corridor": True,
            "fresh_selector_centered_moments_used": True,
            "direct_face_integrations_overlap_transport_balls": True,
            "remaining_398_mode_finite_integral_splice_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "height_builder": {"path": relative(HEIGHT_BUILDER), "sha256": file_hash(HEIGHT_BUILDER)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
        "next_action": (
            "Combine the sharpened event-zero collar with the exact event-ordered 398-term projection. "
            "Derive the finite-integral coherent/leakage kernel before spending the remaining top-corridor budget."
        ),
        "proof_boundary": (
            "One exact-minus-beta4 event mode on t*-pi/16<=t<=t* only. No remaining-mode finite-integral theorem, "
            "completed-remainder estimate, all-corridor continuation, complete Q_K-T or T_upper theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("recentered event zero at selector center with fresh top-corridor moments", flush=True)


if __name__ == "__main__":
    main()
