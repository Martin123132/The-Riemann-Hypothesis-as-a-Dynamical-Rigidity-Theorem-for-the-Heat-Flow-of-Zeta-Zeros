#!/usr/bin/env python3
"""Certify the prototype full-line join throughout one complete height cell."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
from typing import Any, Iterator


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import flint
from flint import acb, acb_series, arb


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_continuous_height_cell_gate"
FULL_STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_full_line_join_gate"
FULL_BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{FULL_STEM}.py"
FULL_RESULT = REPO_ROOT / f"work/rh_compute/results/{FULL_STEM}.json"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 85
MOMENT_ORDER = 12
TARGET_NORMALIZED = arb("0.000019")
TARGET_PHYSICAL = arb("0.0000086")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def load_full_module():
    spec = importlib.util.spec_from_file_location("full_line_transport_source", FULL_BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load full-line builder")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def complex_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": value.real.str(PRECISION, more=True),
        "imag_ball": value.imag.str(PRECISION, more=True),
    }


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def contour_panels(module, face) -> Iterator[tuple[str, acb, arb, acb, int]]:
    for index in range(module.REAL_PANELS):
        center = -module.CONTOUR_RADIUS + module.REAL_HALF_WIDTH + 2 * module.REAL_HALF_WIDTH * index
        yield "real_core", acb(center), module.REAL_HALF_WIDTH, acb(1), 1
    for side, orientation in ((1, 1), (-1, -1)):
        for index in range(module.VERTICAL_PANELS):
            center_y = module.VERTICAL_HALF_WIDTH + 2 * module.VERTICAL_HALF_WIDTH * index
            yield "right_connector" if side > 0 else "left_connector", acb(side * module.CONTOUR_RADIUS, center_y), module.VERTICAL_HALF_WIDTH, face.i, orientation
    for side in (-1, 1):
        start = -module.HORIZONTAL_END if side < 0 else module.CONTOUR_RADIUS
        for index in range(module.HORIZONTAL_PANELS_PER_SIDE):
            center_x = start + module.REAL_HALF_WIDTH + 2 * module.REAL_HALF_WIDTH * index
            yield "horizontal_near", acb(center_x, module.CONTOUR_HEIGHT), module.REAL_HALF_WIDTH, acb(1), 1


def segment_ball(center: acb, half_width: arb, direction: acb) -> acb:
    if direction.imag.contains(0):
        return acb(arb(center.real, half_width), center.imag)
    return acb(center.real, arb(center.imag, half_width))


def moment_panel_error(module, face, center: acb, half_width: arb, direction: acb, degree: int) -> arb:
    is_lifted_horizontal = direction.imag.contains(0) and center.imag.lower() > arb("0.9")
    disk_radius = arb("1.2" if is_lifted_horizontal else "1.8") * half_width
    z_disk = acb(arb(center.real, disk_radius), arb(center.imag, disk_radius))
    q0 = center / face.beta
    tanh0 = -face.i * (face.i * q0).tan()
    exact_b0 = -(face.detuning + face.beta * tanh0)
    canonical_b0 = -(center + face.detuning)
    use_removable = min(abs(exact_b0).upper(), abs(canonical_b0).upper()) < arb("0.6")
    inner_error = face.inner_truncation_bound(z_disk, use_removable)
    disk_power = abs(z_disk).upper() ** degree
    maximum = (face.true_integrand_bound(z_disk) + inner_error) * disk_power
    ratio = half_width / disk_radius
    outer_tail = (
        2
        * maximum
        * half_width
        * ratio**module.SERIES_ORDER
        / ((module.SERIES_ORDER + 1) * (1 - ratio))
    )
    path_power = abs(segment_ball(center, half_width, direction)).upper() ** degree
    return outer_tail + 2 * half_width * path_power * inner_error


def certified_moments(module, face, order: int) -> tuple[list[acb], arb, dict[str, arb]]:
    moments = [acb(0) for _ in range(order + 1)]
    errors = [arb(0) for _ in range(order + 1)]
    absolute_mass = arb(0)
    mass_by_piece: dict[str, arb] = {}
    for piece, center, half_width, direction, orientation in contour_panels(module, face):
        base = face.difference_series(center, direction)
        variable = acb_series([center, direction], prec=module.SERIES_ORDER)
        power = acb_series([1], prec=module.SERIES_ORDER)
        for degree in range(order + 1):
            moments[degree] += orientation * module.integrate_symmetric_series(base * power, half_width)
            errors[degree] += moment_panel_error(module, face, center, half_width, direction, degree)
            power *= variable
        segment = segment_ball(center, half_width, direction)
        panel_mass = 2 * half_width * face.true_integrand_bound(segment)
        absolute_mass += panel_mass
        mass_by_piece[piece] = mass_by_piece.get(piece, arb(0)) + panel_mass
    return [module.add_error(value, error) for value, error in zip(moments, errors)], absolute_mass, mass_by_piece


def polynomial_value(coefficients: list[acb], theta: int) -> acb:
    total = acb(0)
    power = arb(1)
    for coefficient in coefficients:
        total += power * coefficient
        power *= theta
    return total


def render_note(artifact: dict[str, Any]) -> str:
    return f"""# Continuous-height prototype full-line join

Date: 2026-08-12

Status: one prototype crossing mode is certified throughout its full
`pi/8` height cell; this is not a proof of 399-event propagation or RH

Let `tau=pi*m*(C-2m)`, `t*=pi*C^2/8`, `beta=(t*)^(1/3)`, and
`h=t-tau`.  Since `eta=t*/t`, the two phases in the common full-line
chart simplify exactly to

```text
Theta_exact(t)=t*z/beta-t* tanh(z/beta)+normal terms,
Theta_Airy(t)=z^3/3+(t-t*)z/beta+normal terms.       (CH1)
```

Every remaining amplitude and normal Fresnel factor is independent of `t`.
Consequently their signed difference satisfies the exact transport identity

```text
D_t(z,y)=exp(i*h*z/beta) D_tau(z,y).                 (CH2)
```

Write `h=(pi/16) theta`, `|theta|<=1`, and
`kappa=pi/(16 beta)`.  On the finite deformed contour we certified the
moments `M_n=int z^n D_tau(z) dz` through degree
{artifact['method']['moment_order']} and used

```text
Delta I(theta)=sum_(n=0)^N (i*kappa*theta)^n M_n/n!+R_N(theta).
```

The exponential Taylor remainder is bounded against a direct Arb enclosure
of the finite-contour absolute mass.  The already certified tails on
`Im z=1` acquire at most the exact factor `exp(kappa)`.  This gives uniformly
for every real `t` in `[tau-pi/16,tau+pi/16]`

```text
|Delta I(t)| <= {artifact['certification']['uniform_normalized_absolute_bound']},
(sqrt(2)/pi)|Delta I(t)| <= {artifact['certification']['uniform_physical_absolute_bound']}.
```

Both are below `0.000019` and `0.0000086`, respectively.  The interval
reconstruction also overlaps the separately certified lower face, event,
and upper face balls.

## Boundary

This closes continuous height only for the present one-mode prototype and
its present event parameters.  It does not yet make the constants uniform
in the event index, treat all 399 event cells, control the remaining modes,
establish complete `T_upper`, prove `Lambda<=0`, prove RH, or establish a
prize-level conclusion.
"""


def main() -> None:
    flint.ctx.dps = PRECISION
    priority = set_low_priority()
    require(FULL_BUILDER.is_file(), "missing full-line builder")
    require(FULL_RESULT.is_file(), "missing full-line result")
    require(CHECKER.is_file(), "missing independent checker")
    source = json.loads(FULL_RESULT.read_text(encoding="utf-8"))
    require(source.get("status") == "full_line_one_mode_hyperbolic_Morse_Fresnel_Airy_join_certified", "bad full-line dependency")

    module = load_full_module()
    flint.ctx.cap = module.SERIES_ORDER + 8
    face = module.FullLineFace("event", 0)
    moments, absolute_mass, mass_by_piece = certified_moments(module, face, MOMENT_ORDER)

    kappa = face.pi / (16 * face.beta)
    coefficients = [(face.i * kappa) ** degree / arb(math.factorial(degree)) * moment for degree, moment in enumerate(moments)]
    z_max = (module.HORIZONTAL_END**2 + module.CONTOUR_HEIGHT**2).sqrt()
    expansion_size = kappa * z_max
    finite_remainder = absolute_mass * expansion_size.exp() * expansion_size ** (MOMENT_ORDER + 1) / arb(math.factorial(MOMENT_ORDER + 1))
    far_at_event = module.horizontal_far_bound(face)
    uniform_far = kappa.exp() * far_at_event
    uniform_bound = sum((abs(value).upper() for value in coefficients), arb(0)) + finite_remainder + uniform_far
    physical_bound = arb(2).sqrt() / face.pi * uniform_bound
    require(uniform_bound < TARGET_NORMALIZED, "continuous-height normalized target failed")
    require(physical_bound < TARGET_PHYSICAL, "continuous-height physical target failed")

    stored = {row["label"]: row for row in source["certified_faces"]}
    endpoint_checks: list[dict[str, Any]] = []
    for label, theta in (("lower_face", -1), ("event", 0), ("upper_face", 1)):
        reconstructed = module.add_error(polynomial_value(coefficients, theta), finite_remainder + uniform_far)
        recorded = parse_complex(stored[label]["normalized_full_line_difference_ball"])
        real_overlap = reconstructed.real.overlaps(recorded.real)
        imag_overlap = reconstructed.imag.overlaps(recorded.imag)
        require(real_overlap and imag_overlap, f"{label} transport reconstruction does not overlap")
        endpoint_checks.append({
            "label": label,
            "theta": theta,
            "transport_reconstruction_ball": complex_record(reconstructed),
            "stored_full_line_ball": complex_record(recorded),
            "real_overlap": real_overlap,
            "imag_overlap": imag_overlap,
        })

    artifact: dict[str, Any] = {
        "kind": STEM,
        "date": "2026-08-12",
        "status": "one_mode_continuous_height_cell_transport_certified",
        "precision_decimal_digits": PRECISION,
        "resource_policy": {"worker_cap": 1, "active_workers": 1, "process_priority": priority},
        "height_cell": {
            "center_tau_ball": face.tau.str(PRECISION, more=True),
            "half_width_pi_over_16_ball": (face.pi / 16).str(PRECISION, more=True),
            "theta_interval": "[-1,1]",
            "beta_ball": face.beta.str(PRECISION, more=True),
            "kappa_pi_over_16beta_ball": kappa.str(PRECISION, more=True),
        },
        "method": {
            "exact_transport_identity": "D_t(z,y)=exp(i*(t-tau)*z/beta)*D_tau(z,y)",
            "moment_order": MOMENT_ORDER,
            "moment_definition": "M_n=int_Gamma z^n D_tau(z) dz on the finite deformed contour",
            "finite_contour_z_absolute_bound": z_max.str(PRECISION, more=True),
            "maximum_transport_argument": expansion_size.str(PRECISION, more=True),
            "far_tail_transport_factor": "exp(kappa) on Im(z)=1",
        },
        "certification": {
            "moments": [
                {"degree": degree, "moment_ball": complex_record(moment), "transport_coefficient_ball": complex_record(coefficients[degree])}
                for degree, moment in enumerate(moments)
            ],
            "finite_contour_absolute_mass_bound": absolute_mass.str(PRECISION, more=True),
            "absolute_mass_by_piece": {key: value.str(PRECISION, more=True) for key, value in sorted(mass_by_piece.items())},
            "uniform_finite_exponential_remainder_bound": finite_remainder.str(PRECISION, more=True),
            "event_far_tail_absolute_bound": far_at_event.str(PRECISION, more=True),
            "uniform_transported_far_tail_bound": uniform_far.str(PRECISION, more=True),
            "uniform_normalized_absolute_bound": uniform_bound.str(PRECISION, more=True),
            "uniform_physical_absolute_bound": physical_bound.str(PRECISION, more=True),
            "endpoint_overlap_checks": endpoint_checks,
        },
        "claims": {
            "exact_common_height_transport_identity_proved": True,
            "one_mode_full_height_cell_uniformly_certified": True,
            "normalized_difference_below_0_000019_throughout_cell": True,
            "physical_difference_below_0_0000086_throughout_cell": True,
            "all_399_event_cells_proved": False,
            "complete_T_upper_proved": False,
            "RH_proved": False,
        },
        "proof_boundary": (
            "One prototype mode throughout one fixed-parameter height cell only. Event-index uniformity, all 399 cells, "
            "remaining modes, complete T_upper, Lambda<=0, RH, and a prize-level conclusion remain open."
        ),
        "next_action": (
            "Expose every contour and normal-integral constant as a function of the event detuning, then prove monotonic "
            "or interval-atlas control across the finite 399-event parameter range."
        ),
        "dependencies": {
            "full_line_gate": {"path": relative(FULL_RESULT), "sha256": file_hash(FULL_RESULT)},
            "full_line_builder": {"path": relative(FULL_BUILDER), "sha256": file_hash(FULL_BUILDER)},
        },
        "artifacts": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    artifact["artifacts"]["note"] = {"path": relative(NOTE), "sha256": file_hash(NOTE)}
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"certified continuous height cell: |Delta I|<={uniform_bound}; "
        f"physical<={physical_bound}; priority={priority}",
        flush=True,
    )


if __name__ == "__main__":
    main()
