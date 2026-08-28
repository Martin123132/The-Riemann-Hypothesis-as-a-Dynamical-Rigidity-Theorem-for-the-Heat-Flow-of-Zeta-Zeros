#!/usr/bin/env python3
"""Certify the rational four-term A-exterior common-phase integral."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time
from typing import Any, Sequence


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, acb_series, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_exact_exterior_rational_common_phase_contour_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "reduction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_exact_exterior_common_phase_reduction_gate.json",
    "endpoint_remainder": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_exact_exterior_endpoint_tail_remainder_gate.json",
}

HEIGHT = 10_000_000_000
A = 159_577
LOWER_START = 39_853
UPPER_END = 39_936
MODES = tuple(range(LOWER_START, UPPER_END + 1))
Y0_TEXT = "0.0037"
FAR_SPLICE_Q_TEXT = "0.0083"
CONTOUR_HEIGHT_TEXT = "0.001"
PRODUCTION_PRECISION = 60
PRODUCTION_SERIES_ORDER = 36
PRODUCTION_COMPACT_MAX_WIDTH_TEXT = "0.00002"
PRODUCTION_VERTICAL_MAX_WIDTH_TEXT = "0.00001"


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


def complex_record(value: acb, digits: int = 55) -> dict[str, str]:
    return {
        "real_ball": value.real.str(digits, more=True),
        "imag_ball": value.imag.str(digits, more=True),
        "absolute_ball": abs(value).str(digits, more=True),
    }


def add_complex_error(value: acb, radius: arb) -> acb:
    return acb(arb(value.real, radius), arb(value.imag, radius))


def endpoint_terms(x: Any, mode: arb, endpoint: arb) -> list[Any]:
    pi = arb.pi()
    imaginary = acb(0, 1)
    denominator = (endpoint * x - 2 * mode) * (endpoint * x + 2 * mode)
    b0 = 4 * x ** arb("1.5") * (pi * endpoint**2 * x - 2 * imaginary) / (pi * denominator)
    b1 = (
        8 * imaginary * x ** arb("2.5")
        * (
            pi * endpoint**4 * x**3
            - 20 * pi * endpoint**2 * mode**2 * x
            + 2 * imaginary * endpoint**2 * x**2
            + 24 * imaginary * mode**2
        )
        / (pi**2 * denominator**3)
    )
    b2 = (
        16 * x ** arb("3.5")
        * (
            pi * endpoint**6 * x**5
            - 56 * pi * endpoint**4 * mode**2 * x**3
            + 6 * imaginary * endpoint**4 * x**4
            - 560 * pi * endpoint**2 * mode**4 * x
            + 240 * imaginary * endpoint**2 * mode**2 * x**2
            + 480 * imaginary * mode**4
        )
        / (pi**3 * denominator**5)
    )
    b3 = (
        -96 * imaginary * x ** arb("4.5")
        * (
            pi * endpoint**8 * x**7
            - 108 * pi * endpoint**6 * mode**2 * x**5
            + 10 * imaginary * endpoint**6 * x**6
            - 3024 * pi * endpoint**4 * mode**4 * x**3
            + 840 * imaginary * endpoint**4 * mode**2 * x**4
            - 6720 * pi * endpoint**2 * mode**6 * x
            + 5600 * imaginary * endpoint**2 * mode**4 * x**2
            + 4480 * imaginary * mode**6
        )
        / (pi**4 * denominator**7)
    )
    return [b0, b1, b2, b3]


def parameters() -> dict[str, Any]:
    pi, t, endpoint = arb.pi(), arb(HEIGHT), arb(A)
    c = pi * endpoint**2 / 4
    kappa = (1 - 8 * t / (pi * endpoint**2)).sqrt()
    x_star = (1 - kappa) / 2
    q_star = ((1 - x_star) / x_star).log()
    phase_star = c * x_star + t * q_star / 2
    return {
        "pi": pi,
        "t": t,
        "endpoint": endpoint,
        "c": c,
        "x_star": x_star,
        "q_star": q_star,
        "phase_star": phase_star,
        "q_split": arb(FAR_SPLICE_Q_TEXT),
        "eta": arb(CONTOUR_HEIGHT_TEXT),
        "y0": arb(Y0_TEXT),
        "i": acb(0, 1),
    }


def cutoff_q(mode_int: int, p: dict[str, Any]) -> arb:
    mode = arb(mode_int)
    return (p["t"] * (1 + p["y0"]) / (2 * p["pi"] * mode**2)).log()


def pole_q(mode_int: int, p: dict[str, Any]) -> arb:
    mode = arb(mode_int)
    return ((p["endpoint"] - 2 * mode) / (2 * mode)).log()


def endpoint_sum(q: Any, modes: Sequence[int], p: dict[str, Any]) -> Any:
    x = 1 / (1 + q.exp())
    total = q * 0
    for mode_int in modes:
        total += sum(endpoint_terms(x, arb(mode_int), p["endpoint"]), q * 0)
    return total


def transformed_amplitude(q: Any, modes: Sequence[int], p: dict[str, Any]) -> Any:
    return -((arb("0.75") * q).exp()) * endpoint_sum(q, modes, p) / (2 * p["pi"] * p["i"])


def relative_phase(q: Any, p: dict[str, Any]) -> Any:
    x = 1 / (1 + q.exp())
    return p["c"] * x + p["t"] * q / 2 - p["phase_star"]


def phase_prime(q: Any, p: dict[str, Any]) -> Any:
    return p["t"] / 2 - p["c"] / (4 * (q / 2).cosh() ** 2)


def phase_second(q: Any, p: dict[str, Any]) -> Any:
    return p["c"] / 4 / (q / 2).cosh() ** 2 * (q / 2).tanh()


def integrand_series(
    center: acb,
    direction: acb,
    modes: Sequence[int],
    p: dict[str, Any],
    order: int,
) -> acb_series:
    q = acb_series([center, direction], prec=order)
    series = direction * transformed_amplitude(q, modes, p) * (p["i"] * relative_phase(q, p)).exp()
    require(series.prec == order, f"Taylor series cap drift: expected {order}, got {series.prec}")
    return series


def integrate_symmetric_series(series: acb_series, half_width: arb) -> acb:
    total = acb(0)
    for degree, coefficient in enumerate(series.coeffs()):
        if degree % 2 == 0:
            total += coefficient * 2 * half_width ** (degree + 1) / (degree + 1)
    return total


def panel_certificate(
    center: acb,
    direction: acb,
    half_width: arb,
    modes: Sequence[int],
    p: dict[str, Any],
    order: int,
) -> tuple[acb, arb, arb, arb]:
    disk_radius = 2 * half_width
    series = integrand_series(center, direction, modes, p, order)
    polynomial = integrate_symmetric_series(series, half_width)

    q_disk = acb(
        arb(center.real, disk_radius),
        arb(center.imag, disk_radius),
    )
    amplitude_bound = abs(transformed_amplitude(q_disk, modes, p)).upper()
    derivative_bound = (
        abs(phase_prime(center, p)).upper()
        + disk_radius * abs(phase_second(q_disk, p)).upper()
    )
    center_phase_modulus = abs((p["i"] * relative_phase(center, p)).exp()).upper()
    maximum = amplitude_bound * center_phase_modulus * (disk_radius * derivative_bound).exp()
    ratio = half_width / disk_radius
    remainder = (
        2 * maximum * half_width * ratio**order
        / ((order + 1) * (1 - ratio))
    )
    require(remainder.is_finite(), "nonfinite Cauchy panel remainder")
    return add_complex_error(polynomial, remainder), remainder, maximum, disk_radius * derivative_bound


def panel_count(width: arb, max_width: arb) -> int:
    require(width.lower() > 0, "nonpositive integration width")
    return max(1, math.ceil(float(width.upper() / max_width.lower())))


def integrate_real_band(
    left: arb,
    right: arb,
    modes: Sequence[int],
    p: dict[str, Any],
    order: int,
    max_width: arb,
) -> tuple[acb, dict[str, Any]]:
    width = right - left
    panels = panel_count(width, max_width)
    total = acb(0)
    cauchy_error = arb(0)
    max_disk_phase = arb(0)
    max_integrand = arb(0)
    for index in range(panels):
        panel_left = left + width * index / panels
        panel_right = left + width * (index + 1) / panels
        center = acb((panel_left + panel_right) / 2)
        half_width = (panel_right - panel_left) / 2
        value, error, maximum, disk_phase = panel_certificate(
            center, acb(1), half_width, modes, p, order
        )
        total += value
        cauchy_error += error
        max_integrand = max(max_integrand, maximum)
        max_disk_phase = max(max_disk_phase, disk_phase)
    return total, {
        "panel_count": panels,
        "cauchy_error_ball": cauchy_error.str(45, more=True),
        "max_disk_integrand_ball": max_integrand.str(45, more=True),
        "max_disk_phase_variation_ball": max_disk_phase.str(45, more=True),
    }


def compact_cutoff_sweep(
    p: dict[str, Any], order: int, max_width: arb
) -> tuple[acb, list[dict[str, Any]], int, arb]:
    cutoffs = sorted(
        ((cutoff_q(mode, p), mode) for mode in MODES),
        key=lambda row: float(row[0].mid()),
    )
    require([mode for _, mode in cutoffs] == list(range(UPPER_END, LOWER_START - 1, -1)), "cutoff order drift")
    require(cutoffs[-1][0].upper() < p["q_split"].lower(), "far splice does not clear all cutoffs")

    active: list[int] = []
    total = acb(0)
    rows: list[dict[str, Any]] = []
    panel_total = 0
    cauchy_total = arb(0)
    for index, (left, mode) in enumerate(cutoffs):
        active.append(mode)
        right = cutoffs[index + 1][0] if index + 1 < len(cutoffs) else p["q_split"]
        value, diagnostics = integrate_real_band(left, right, tuple(active), p, order, max_width)
        total += value
        panel_total += diagnostics["panel_count"]
        cauchy_total += arb(diagnostics["cauchy_error_ball"])
        rows.append({
            "band_index": index,
            "new_mode": mode,
            "active_mode_count": len(active),
            "q_left_ball": left.str(45, more=True),
            "q_right_ball": right.str(45, more=True),
            "relative_contribution": complex_record(value, 45),
            **diagnostics,
        })
    return total, rows, panel_total, cauchy_total


def vertical_leg(
    p: dict[str, Any], order: int, max_width: arb
) -> tuple[acb, dict[str, Any]]:
    panels = panel_count(p["eta"], max_width)
    total = acb(0)
    cauchy_error = arb(0)
    max_disk_phase = arb(0)
    max_integrand = arb(0)
    for index in range(panels):
        left = p["eta"] * index / panels
        right = p["eta"] * (index + 1) / panels
        center = acb(p["q_split"], (left + right) / 2)
        half_width = (right - left) / 2
        value, error, maximum, disk_phase = panel_certificate(
            center, p["i"], half_width, MODES, p, order
        )
        total += value
        cauchy_error += error
        max_integrand = max(max_integrand, maximum)
        max_disk_phase = max(max_disk_phase, disk_phase)
    return total, {
        "panel_count": panels,
        "cauchy_error_ball": cauchy_error.str(45, more=True),
        "max_disk_integrand_ball": max_integrand.str(45, more=True),
        "max_disk_phase_variation_ball": max_disk_phase.str(45, more=True),
        "relative_contribution": complex_record(total, 50),
    }


def horizontal_amplitude_mode_majorant(mode_int: int, x_bound: arb, decay_rate: arb, p: dict[str, Any]) -> arb:
    pi, endpoint, mode = p["pi"], p["endpoint"], arb(mode_int)
    denominator = 4 * mode**2 - endpoint**2 * x_bound**2
    require(denominator.lower() > 0, f"horizontal denominator lost at mode {mode_int}")
    groups: list[tuple[arb, arb]] = []
    groups.extend((
        (4 * endpoint**2 / denominator, arb("2.5")),
        (8 / (pi * denominator), arb("1.5")),
    ))
    prefactor = 8 / (pi**2 * denominator**3)
    groups.extend((
        (prefactor * pi * endpoint**4, arb("5.5")),
        (prefactor * 20 * pi * endpoint**2 * mode**2, arb("3.5")),
        (prefactor * 2 * endpoint**2, arb("4.5")),
        (prefactor * 24 * mode**2, arb("2.5")),
    ))
    prefactor = 16 / (pi**3 * denominator**5)
    groups.extend((
        (prefactor * pi * endpoint**6, arb("8.5")),
        (prefactor * 56 * pi * endpoint**4 * mode**2, arb("6.5")),
        (prefactor * 6 * endpoint**4, arb("7.5")),
        (prefactor * 560 * pi * endpoint**2 * mode**4, arb("4.5")),
        (prefactor * 240 * endpoint**2 * mode**2, arb("5.5")),
        (prefactor * 480 * mode**4, arb("3.5")),
    ))
    prefactor = 96 / (pi**4 * denominator**7)
    groups.extend((
        (prefactor * pi * endpoint**8, arb("11.5")),
        (prefactor * 108 * pi * endpoint**6 * mode**2, arb("9.5")),
        (prefactor * 10 * endpoint**6, arb("10.5")),
        (prefactor * 3024 * pi * endpoint**4 * mode**4, arb("7.5")),
        (prefactor * 840 * endpoint**4 * mode**2, arb("8.5")),
        (prefactor * 6720 * pi * endpoint**2 * mode**6, arb("5.5")),
        (prefactor * 5600 * endpoint**2 * mode**4, arb("6.5")),
        (prefactor * 4480 * mode**6, arb("4.5")),
    ))
    total = arb(0)
    for coefficient, power in groups:
        exponent = power * decay_rate - arb("0.75")
        require(exponent.lower() > 0, "horizontal amplitude decay exponent lost")
        total += coefficient * x_bound**power / exponent
    return (arb("0.75") * p["q_split"]).exp() * total / (2 * pi)


def horizontal_ray_bound(p: dict[str, Any]) -> dict[str, str]:
    q, eta, c, t = p["q_split"], p["eta"], p["c"], p["t"]
    exponential = q.exp()
    cosine, sine = eta.cos(), eta.sin()
    x_bound = 1 / (1 + exponential * cosine)
    decay_rate = exponential * cosine / (1 + exponential * cosine)
    phase_damping = (
        t * eta / 2
        - c * exponential * sine
        / (1 + 2 * exponential * cosine + exponential**2)
    )
    require(phase_damping.lower() > 80, "horizontal phase damping lost")
    amplitude_majorant = sum(
        (horizontal_amplitude_mode_majorant(mode, x_bound, decay_rate, p) for mode in MODES),
        arb(0),
    )
    integral_bound = (-phase_damping).exp() * amplitude_majorant
    require(integral_bound.upper() < arb("2e-31"), "horizontal ray bound misses cap")
    return {
        "x_modulus_envelope_at_splice_ball": x_bound.str(50, more=True),
        "exponential_decay_rate_ball": decay_rate.str(50, more=True),
        "phase_damping_exponent_ball": phase_damping.str(50, more=True),
        "undamped_amplitude_integral_majorant_ball": amplitude_majorant.str(50, more=True),
        "relative_horizontal_integral_bound_ball": integral_bound.str(50, more=True),
    }


def certificate(
    precision: int,
    series_order: int,
    compact_max_width_text: str,
    vertical_max_width_text: str,
) -> dict[str, Any]:
    ctx.dps = precision
    ctx.threads = 1
    ctx.cap = series_order
    p = parameters()
    compact_width = arb(compact_max_width_text)
    vertical_width = arb(vertical_max_width_text)

    cutoff_rows = [(mode, cutoff_q(mode, p), pole_q(mode, p)) for mode in MODES]
    pole_margin = min((q0 - pole).lower() for _, q0, pole in cutoff_rows)
    require(pole_margin > arb("0.00369"), "cutoff/pole analytic margin lost")
    require(p["eta"].upper() < arb.pi().lower(), "contour reaches logistic singularity")

    compact, bands, compact_panels, compact_error = compact_cutoff_sweep(
        p, series_order, compact_width
    )
    vertical, vertical_diagnostics = vertical_leg(p, series_order, vertical_width)
    horizontal = horizontal_ray_bound(p)
    horizontal_error = arb(horizontal["relative_horizontal_integral_bound_ball"])
    relative_rational = add_complex_error(compact + vertical, horizontal_error)

    carrier = (p["i"] * p["phase_star"]).exp()
    canonical_rational = carrier * relative_rational
    normalizer = (p["pi"] / (32 * p["t"])) ** arb("0.25")
    paper_rotation = (-p["i"] * p["pi"] / 8).exp()
    physical_rational = 2 * normalizer * (paper_rotation * canonical_rational).real

    replacement = load_json(DEPENDENCIES["endpoint_remainder"])
    replacement_error = arb(replacement["certificate"]["physical_complete_exterior_replacement_error_ball"]).upper()
    physical_exact = arb(physical_rational, replacement_error)
    require(physical_exact.is_finite(), "full exact exterior enclosure is nonfinite")
    require(physical_exact.upper() < arb("-0.0023"), "exact exterior sign/scale guard lost")
    require(physical_exact.lower() > arb("-0.0025"), "exact exterior lower scale guard lost")

    phase_at_split = relative_phase(p["q_split"], p)
    return {
        "height": HEIGHT,
        "endpoint": A,
        "mode_range": [LOWER_START, UPPER_END],
        "mode_count": len(MODES),
        "y0": Y0_TEXT,
        "q_star_ball": p["q_star"].str(55, more=True),
        "q_split": FAR_SPLICE_Q_TEXT,
        "phase_at_split_ball": phase_at_split.str(55, more=True),
        "contour_height": CONTOUR_HEIGHT_TEXT,
        "minimum_cutoff_to_endpoint_pole_q_margin_ball": pole_margin.str(50, more=True),
        "series_order": series_order,
        "compact_max_panel_width": compact_max_width_text,
        "vertical_max_panel_width": vertical_max_width_text,
        "compact_panel_count": compact_panels,
        "compact_cauchy_error_sum_ball": compact_error.str(50, more=True),
        "compact_relative_integral_ball": complex_record(compact),
        "vertical_leg": vertical_diagnostics,
        "horizontal_ray": horizontal,
        "relative_rational_integral_ball": complex_record(relative_rational),
        "canonical_rational_integral_ball": complex_record(canonical_rational),
        "physical_rational_integral_ball": physical_rational.str(55, more=True),
        "endpoint_replacement_physical_error_ball": replacement_error.str(55, more=True),
        "physical_full_exact_exterior_ball": physical_exact.str(55, more=True),
        "cutoff_bands": bands,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    h = c["horizontal_ray"]
    v = c["vertical_leg"]
    return f"""# Rational common-phase contour certificate for the exact A exterior

Date: 2026-08-14

Status: rigorous four-term rational common-phase value and full exact-exterior
enclosure after the certified endpoint replacement; not a complete A theorem

Put `q=log((1-x)/x)`.  The four-term replacement from Section 11.438 gives

```text
I_4=sum_m integral_(q0_m)^infinity G_m(q)exp(i(Phi(q)-Phi_*))dq,
G_m(q)=-exp(3q/4)(b_0+b_1+b_2+b_3)/(2pi i).          (RC1)
```

The compact interval ends at the exact rational splice `q={FAR_SPLICE_Q_TEXT}`,
where `Phi-Phi_*={c['phase_at_split_ball']}`.  Splitting at all 84 exact
cutoffs and summing the active modes before integration, {c['compact_panel_count']}
Taylor panels of order {c['series_order']} give

```text
I_compact={c['compact_relative_integral_ball']['real_ball']}
          +i {c['compact_relative_integral_ball']['imag_ball']}.             (RC2)
```

Every panel integrates its Arb power series exactly.  Its omitted analytic
tail is bounded by Cauchy's estimate on a disk of twice the panel radius;
the total compact Cauchy error is
`{c['compact_cauchy_error_sum_ball']}`.

For the far interval, shift to `q=s+i*{CONTOUR_HEIGHT_TEXT}`.  No endpoint pole
or logistic singularity lies in the contour strip: the minimum real
cutoff-to-pole margin is `{c['minimum_cutoff_to_endpoint_pole_q_margin_ball']}`.
The finite vertical leg gives

```text
I_vertical={v['relative_contribution']['real_ball']}
           +i {v['relative_contribution']['imag_ball']}.                     (RC3)
```

On the horizontal ray, monotonicity of
`R/(1+2R cos(eta)+R^2)` and the exact rational amplitude imply

```text
Im(Phi(s+i eta)-Phi_*) >= {h['phase_damping_exponent_ball']},
integral_(q_split)^infinity |G(s+i eta)exp(i(Phi-Phi_*))|ds
 <= {h['relative_horizontal_integral_bound_ball']}.                          (RC4)
```

Restoring the common phase, paper rotation, and equation-(9) normalization
gives

```text
four-term rational exterior = {c['physical_rational_integral_ball']},
endpoint replacement error <= {c['endpoint_replacement_physical_error_ball']},
full exact A exterior        = {c['physical_full_exact_exterior_ball']}.      (RC5)
```

An independent checker increases precision and Taylor order and decreases
both panel widths.  Its physical and canonical enclosures overlap the saved
balls and have smaller radii.

Pi provenance: every `pi` in (RC1)--(RC5) comes from the exact endpoint
hierarchy, common outer phase, paper rotation, and equation-(9) normalization.
No fitted constant is introduced.

Proof boundary: the full exact exterior for the finite A-transition roster
39853..39936 at t=10^10 only.  No compact exact-minus-affine transformed-
amplitude value, complete A endpoint block, `R_Dir` estimate, complete
`Q_K-T` or `T_upper`, all-height theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["reduction"]["decision"]["full_exact_exterior_common_phase_reduction_proved"] is True, "common-phase dependency drift")
    require(dependencies["endpoint_remainder"]["decision"]["four_term_endpoint_tail_used"] is True, "endpoint-remainder dependency drift")

    certified = certificate(
        PRODUCTION_PRECISION,
        PRODUCTION_SERIES_ORDER,
        PRODUCTION_COMPACT_MAX_WIDTH_TEXT,
        PRODUCTION_VERTICAL_MAX_WIDTH_TEXT,
    )
    artifact = {
        "kind": STEM,
        "status": "full_exact_A_exterior_rigorously_enclosed_by_rational_common_phase_contour",
        "passed": True,
        "certificate": certified,
        "decision": {
            "logistic_contour_strip_analytic_for_all_84_modes": True,
            "grouped_compact_common_phase_integral_certified": True,
            "finite_vertical_contour_leg_certified": True,
            "horizontal_ray_below_2e_minus_31_canonical": True,
            "four_term_rational_exterior_value_certified": True,
            "full_exact_exterior_value_after_endpoint_replacement_certified": True,
            "compact_exact_minus_affine_transformed_amplitude_certified": False,
            "complete_A_endpoint_block_proved": False,
            "R_Dir_bound_proved": False,
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
            "precision_decimal_digits": PRODUCTION_PRECISION,
        },
        "next_obligation": "Certify the compact exact-minus-affine transformed-amplitude integral on the localized exact A domain and combine it with the finite-box affine carrier and this exact exterior value.",
        "proof_boundary": "Full exact exterior value for the 84-mode A-transition roster at t=10^10 only. No compact exact-minus-affine transformed-amplitude value, complete A endpoint theorem, R_Dir estimate, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified rational common-phase contour and full exact A exterior", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
