#!/usr/bin/env python3
"""Transport the block-20 endpoint model over every occupied coefficient cell."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any, Callable


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_block20_retained_decay_endpoint_majorant_gate as retained
import jensen_window_pf_newman_c1_hardy_block20_transported_endpoint_budget_gate as transport
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import acb, arb, ctx


ATLAS = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_discrete_selector_atlas_gate.json"
CENTER_ENDPOINT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_retained_decay_endpoint_majorant_gate.json"
WEIGHT_FIXTURE = transport.WEIGHT_FIXTURE
WEIGHTS = transport.WEIGHTS
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_coefficient_neighborhood_transport_gate.json"
CACHE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_coefficient_neighborhood_transport_gate_cache.jsonl"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_coefficient_neighborhood_transport_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_block20_coefficient_neighborhood_transport_gate.py"
PRECISIONS = (60, 90)
REQUESTED_ERROR_SCALE = Fraction(5, 1000)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def exact_fraction(record: dict[str, Any]) -> Fraction:
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def upper(value: arb) -> Fraction:
    return cells.bound_fraction(value.upper())


def lower(value: arb) -> Fraction:
    return cells.bound_fraction(value.lower())


def upper_abs(value: arb | acb) -> Fraction:
    return upper(abs(value))


def bounds(value: arb) -> tuple[Fraction, Fraction]:
    return cells.bounds(value)


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover - platform fallback
        return f"unavailable:{type(exc).__name__}"


def interval(center: Fraction, radius: Fraction) -> arb:
    return cells.interval_ball(center, radius)


def acb_bounds_record(value: acb) -> dict[str, Any]:
    return point_balls.acb_record(value, 45)


def arb_bounds_record(value: arb) -> dict[str, Any]:
    return point_balls.arb_record(value, 45)


def absolute_interval_extrema(center: Fraction, radius: Fraction) -> tuple[Fraction, Fraction]:
    low = center - radius
    high = center + radius
    require(not (low <= 0 <= high), "coefficient sign crosses zero inside selector cell")
    return min(abs(low), abs(high)), max(abs(low), abs(high))


def integrate_envelope(function: Callable[[acb, bool], acb], dps: int) -> arb:
    return retained.integrate_real(function, retained.ORIGIN, retained.CUTOFF, dps)


def generic_envelope(
    quadratic_low: Fraction,
    quadratic_high: Fraction,
    cubic_low: Fraction,
    cubic_high: Fraction,
    gap_low: Fraction,
    family: str,
    cubic_sign: int,
    dps: int,
) -> arb:
    require(
        0 < quadratic_low <= quadratic_high and 0 <= cubic_low <= cubic_high and gap_low > 0,
        "generic coefficient envelope domain failure",
    )
    ctx.dps = dps
    ctx.threads = 1
    sine_one, sine_two, sine_three, _ = retained.ray_sines(family, cubic_sign)
    pi = arb.pi()
    alpha = 2 * pi * sine_one
    beta_low = 2 * pi * cells.fraction_ball(quadratic_low) * sine_two
    gamma_low = 2 * pi * cells.fraction_ball(cubic_low) * sine_three
    quadratic_high_ball = cells.fraction_ball(quadratic_high)
    cubic_high_ball = cells.fraction_ball(cubic_high)
    gap_low_ball = cells.fraction_ball(gap_low)

    def integrand(radius: acb, _analytic: bool) -> acb:
        minimum_decay = beta_low * radius**2 + gamma_low * radius**3
        homotopy_factor = (1 - (-minimum_decay).exp()) / minimum_decay
        geometric_denominator = 1 - (-alpha * radius).exp()
        return (
            2
            * pi
            * (quadratic_high_ball * radius**2 + cubic_high_ball * radius**3)
            * (-alpha * gap_low_ball * radius).exp()
            * homotopy_factor
            / geometric_denominator
        )

    central = integrate_envelope(integrand, dps)
    epsilon = cells.fraction_ball(retained.ORIGIN)
    cutoff = cells.fraction_ball(retained.CUTOFF)
    origin_radius = upper(
        2
        * pi
        * (alpha * epsilon).exp()
        / alpha
        * (quadratic_high_ball * epsilon**2 / 2 + cubic_high_ball * epsilon**3 / 3)
    )
    rate = alpha * gap_low_ball
    denominator_floor = 1 - (-alpha * cutoff).exp()
    tail_radius = upper(
        2
        * pi
        / denominator_floor
        * (
            quadratic_high_ball * retained.exponential_moment_two(rate, cutoff)
            + cubic_high_ball * retained.exponential_moment_three(rate, cutoff)
        )
    )
    return central + retained.error_ball(origin_radius + tail_radius)


def exceptional_b_envelope(
    cubic_high: Fraction,
    length: int,
    gap_low: Fraction,
    quadratic_floor_low: Fraction,
    cubic_sign: int,
    dps: int,
) -> arb:
    require(cubic_high > 0 and gap_low >= 0 and quadratic_floor_low > 0, "W2 envelope domain failure")
    ctx.dps = dps
    ctx.threads = 1
    sine_one, _, _, _ = retained.ray_sines("b", cubic_sign)
    sine_two = cells.fraction_ball(Fraction(1, 2) if cubic_sign >= 0 else Fraction(1))
    pi = arb.pi()
    cubic = cells.fraction_ball(cubic_high)
    linear_rate = 2 * pi * sine_one * cells.fraction_ball(gap_low)
    quadratic_rate = 2 * pi * cells.fraction_ball(quadratic_floor_low) * sine_two

    def integrand(radius: acb, _analytic: bool) -> acb:
        return (
            2
            * pi
            * cubic
            * (3 * length * radius**2 + radius**3)
            * (-linear_rate * radius - quadratic_rate * radius**2).exp()
        )

    central = retained.integrate_real(integrand, Fraction(0), retained.CUTOFF, dps)
    cutoff = cells.fraction_ball(retained.CUTOFF)
    tail = upper(
        2
        * pi
        * cubic
        * (
            3 * length * retained.gaussian_moment_two(quadratic_rate, cutoff)
            + retained.gaussian_moment_three(quadratic_rate, cutoff)
        )
    )
    return central + retained.error_ball(tail)


def exceptional_c_envelope(
    cubic_high: Fraction,
    gap_low: Fraction,
    phi2_low: Fraction,
    cubic_sign: int,
    dps: int,
) -> arb:
    require(cubic_high > 0 and gap_low >= 0 and phi2_low > 0, "lower exceptional envelope domain failure")
    ctx.dps = dps
    ctx.threads = 1
    sine_one, _, _, _ = retained.ray_sines("c", cubic_sign)
    pi = arb.pi()
    cubic = cells.fraction_ball(cubic_high)
    linear_rate = 2 * pi * sine_one * cells.fraction_ball(gap_low)
    quadratic_rate = pi * cells.fraction_ball(phi2_low)

    def integrand(radius: acb, _analytic: bool) -> acb:
        return 2 * pi * cubic * radius**3 * (-linear_rate * radius - quadratic_rate * radius**2).exp()

    central = retained.integrate_real(integrand, Fraction(0), retained.CUTOFF, dps)
    cutoff = cells.fraction_ball(retained.CUTOFF)
    tail = upper(2 * pi * cubic * retained.gaussian_moment_three(quadratic_rate, cutoff))
    return central + retained.error_ball(tail)


def coefficient_box(data: dict[str, Any], atlas_row: dict[str, Any]) -> dict[str, Any]:
    parent = data["levels"][1]
    cell = atlas_row["q_cell"]
    a1, phi2, a3 = (cells.binary128_fraction(value) for value in parent["coefficients_hex"])
    y = 2 * phi2
    frac = cells.binary128_fraction(parent["frac_length_hex"])
    require(parent["coefficients_hex"][0] == cell["center"]["a1_hex"], "a1 center drift")
    require(parent["xr_hex"] == cell["center"]["xr_hex"], "xr center drift")
    require(parent["coefficients_hex"][2] == cell["center"]["a3_hex"], "a3 center drift")
    require(parent["frac_length_hex"] == cell["center"]["fracL_hex"], "fracL center drift")
    radii = {name: exact_fraction(cell["radii"][name]["exact"]) for name in ("a1", "y", "a3", "fracL")}
    a1_bounds = (a1 - radii["a1"], a1 + radii["a1"])
    y_bounds = (y - radii["y"], y + radii["y"])
    a3_bounds = (a3 - radii["a3"], a3 + radii["a3"])
    frac_bounds = (frac - radii["fracL"], frac + radii["fracL"])
    require(y_bounds[0] > 0 and 0 < frac_bounds[0] < frac_bounds[1] < 1, "basic cell domain drift")
    require(not (a1_bounds[0] <= 0 <= a1_bounds[1]), "a1 selector changes sign")
    require(not (a3_bounds[0] <= 0 <= a3_bounds[1]), "a3 ray selector changes sign")
    return {
        "center": {"a1": a1, "y": y, "a3": a3, "fracL": frac},
        "radii": radii,
        "bounds": {"a1": a1_bounds, "y": y_bounds, "a3": a3_bounds, "fracL": frac_bounds},
    }


def structural_transport(data: dict[str, Any], atlas_row: dict[str, Any], box: dict[str, Any], dps: int) -> dict[str, Any]:
    ctx.dps = dps
    ctx.threads = 1
    center = box["center"]
    radii = box["radii"]
    a1 = interval(center["a1"], radii["a1"])
    y = interval(center["y"], radii["y"])
    a3 = interval(center["a3"], radii["a3"])
    frac = interval(center["fracL"], radii["fracL"])
    selectors = atlas_row["q_cell"]["child_selectors"]
    n1 = int(selectors["n1"])
    ix = int(selectors["ix"])
    sigma = int(selectors["linear_residual_sign"])
    orientation = 1 if bool(selectors["source_conjugates"]) else -1

    x = -a1
    s1 = 1 / (2 * y)
    s2 = -a3 / y**3
    raw0 = s1 * x**2 + s2 * x**3
    raw1 = 2 * s1 * x + 3 * s2 * x**2
    raw2 = s1 + 3 * s2 * x
    raw3 = s2
    linear = raw1 - n1
    if ix % 2:
        linear -= cells.fraction_ball(Fraction(sigma, 2))
    child_xr = 2 * raw2 - ix
    stored = [orientation * linear, orientation * child_xr / 2, orientation * raw3]

    child = data["levels"][2]
    logged_coefficients = [point_balls.binary128_ball(value) for value in child["coefficients_hex"]]
    for index, logged in enumerate(logged_coefficients):
        require(stored[index].contains(logged), f"child coefficient {index + 1} escaped reconstructed cell")

    h1 = 3 * a3 / y**2
    child_kernel = acb(0)
    for index in range(int(child["length"]) + 1):
        phase = arb(0)
        for coefficient in reversed(stored):
            phase = (phase + coefficient) * index
        u = index - a1
        denominator = 1 + h1 * u - (arb(3) / 2) * h1**2 * u**2
        require(not denominator.contains(0), "child denominator interval contains zero")
        tpm = point_balls.binary128_ball(data["header"]["tpm_hex"])
        child_kernel += acb(arb(0), tpm * phase).exp() / denominator
    if child["subtract_one"]:
        child_kernel -= 1
    if child["conjugate"]:
        child_kernel = child_kernel.conjugate()
    source_child = point_balls.binary128_complex(data["steps"][1]["state_before_hex"])
    require(child_kernel.contains(source_child), "source child kernel escaped interval transport")

    y_center = cells.fraction_ball(center["y"])
    a1_center = cells.fraction_ball(center["a1"])
    a3_center = cells.fraction_ball(center["a3"])
    x_center = -a1_center
    raw0_center = x_center**2 / (2 * y_center) + (-a3_center / y_center**3) * x_center**3
    source_multiplier = point_balls.binary128_complex(data["steps"][1]["multiplier_hex"])
    tpp = -point_balls.binary128_ball(data["header"]["tpm_hex"])
    anchored_ratio = (y_center / y).sqrt() * acb(arb(0), tpp * (raw0 - raw0_center)).exp()
    multiplier = source_multiplier * anchored_ratio
    require(multiplier.contains(source_multiplier), "source multiplier escaped anchored interval transport")

    parent_length = int(data["levels"][1]["length"])
    child_length = int(child["length"])
    jbot = 1 if center["a1"] > 0 else 0
    integral_count = child_length - jbot + 1
    phase_lipschitz = -point_balls.binary128_ball(data["header"]["tpm_hex"]) * (
        cells.fraction_ball(radii["a1"]) * parent_length
        + cells.fraction_ball(radii["y"] / 2) * parent_length**2
        + cells.fraction_ball(radii["a3"]) * parent_length**3
    )
    saddle_lipschitz = -point_balls.binary128_ball(data["header"]["tpm_hex"]) * integral_count * (
        cells.fraction_ball(radii["a1"]) * parent_length**2 / 2
        + cells.fraction_ball(radii["y"] / 2) * parent_length**3 / 3
        + cells.fraction_ball(radii["a3"]) * parent_length**4 / 4
    )
    return {
        "child_coefficients": stored,
        "child_kernel": child_kernel,
        "source_child": source_child,
        "multiplier": multiplier,
        "source_multiplier": source_multiplier,
        "phase_lipschitz": phase_lipschitz,
        "saddle_lipschitz": saddle_lipschitz,
        "frac": frac,
    }


def endpoint_envelope(data: dict[str, Any], box: dict[str, Any], dps: int) -> dict[str, Any]:
    bounds_map = box["bounds"]
    a1_low, a1_high = bounds_map["a1"]
    y_low, y_high = bounds_map["y"]
    a3_low, a3_high = bounds_map["a3"]
    phi2_low, phi2_high = y_low / 2, y_high / 2
    cubic_low, cubic_high = absolute_interval_extrema(box["center"]["a3"], box["radii"]["a3"])
    cubic_sign = 1 if box["center"]["a3"] > 0 else -1
    length = int(data["levels"][1]["length"])
    child_length = int(data["levels"][2]["length"])
    mode_b = child_length + 1
    xi_low = a1_low + y_low * length + 3 * a3_low * length**2
    xi_high = a1_high + y_high * length + 3 * a3_high * length**2
    endpoint_quadratic_low = phi2_low + 3 * a3_low * length
    endpoint_quadratic_high = phi2_high + 3 * a3_high * length
    require(endpoint_quadratic_low > 0, "endpoint quadratic floor changed sign")

    b_gap_low = Fraction(mode_b) - xi_high
    b_zero_gap_low = Fraction(mode_b) - a1_high
    b_length_gap_low = Fraction(mode_b + 1) - xi_high
    b_floor_low = phi2_low if cubic_sign >= 0 else endpoint_quadratic_low
    b_exceptional = exceptional_b_envelope(
        cubic_high, length, b_gap_low, b_floor_low, cubic_sign, dps
    )
    b_zero = generic_envelope(
        phi2_low, phi2_high, cubic_low, cubic_high, b_zero_gap_low, "b", cubic_sign, dps
    )
    b_length = generic_envelope(
        endpoint_quadratic_low,
        endpoint_quadratic_high,
        cubic_low,
        cubic_high,
        b_length_gap_low,
        "b",
        cubic_sign,
        dps,
    )

    if a1_high < 0:
        selector = "W3"
        c_exceptional_gap_low = 1 + a1_low
        c_zero_gap_low = 2 + a1_low
        c_length_gap_low = 1 + xi_low
    elif a1_low > 0:
        selector = "W4"
        c_exceptional_gap_low = a1_low
        c_zero_gap_low = 1 + a1_low
        c_length_gap_low = xi_low
    else:  # pragma: no cover - guarded by coefficient_box
        raise RuntimeError("lower selector changes inside cell")
    c_exceptional = exceptional_c_envelope(
        cubic_high, c_exceptional_gap_low, phi2_low, cubic_sign, dps
    )
    c_zero = generic_envelope(
        phi2_low, phi2_high, cubic_low, cubic_high, c_zero_gap_low, "c", cubic_sign, dps
    )
    c_length = generic_envelope(
        endpoint_quadratic_low,
        endpoint_quadratic_high,
        cubic_low,
        cubic_high,
        c_length_gap_low,
        "c",
        cubic_sign,
        dps,
    )
    upper_family = b_exceptional + b_zero + b_length
    lower_family = c_exceptional + c_zero + c_length
    return {
        "selector": selector,
        "phi2_low": phi2_low,
        "endpoint_quadratic_low": endpoint_quadratic_low,
        "b_gap_low": b_gap_low,
        "c_exceptional_gap_low": c_exceptional_gap_low,
        "xi_low": xi_low,
        "xi_high": xi_high,
        "upper_family": upper_family,
        "lower_family": lower_family,
        "complete": upper_family + lower_family,
    }


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 coefficient-neighborhood transport gate

Date: 2026-08-09

Status: {artifact['status']}; finite block-20 theorem only, not a proof of the complete evaluator or RH

## Cellwise contract

Each of the {aggregate['cell_count']} occupied recursive selector cells is used
as an exact box in `(Phi1, xr=2*Phi2, Phi3, fracL)`.  The fixed selector word is
replayed through the cubic Legendre child map.  The transformed child kernel,
the source-anchored recurrence multiplier, the parent phase, the exact
equation-(69) saddle family, and the retained-decay endpoint functional are
recorded in separate columns.

For the saddle family, `|exp(iu)-exp(iv)| <= |u-v|` gives the source-normalized
explicit bound

```text
|Delta sum I_n| <= tpp_src*J*(r1*N^2/2+r2*N^3/3+r3*N^4/4),
```

where `J` is the fixed number of admitted saddle indices and `tpp_src` is the
accepted binary128 value produced by the evaluator's `8*atan(1)`.  Mathematical
`2*pi` and its independently bounded source-constant transport remain distinct.

For the endpoint functional, every positive numerator coefficient is replaced
by its cell maximum, while every gap and every retained quadratic/cubic decay
coefficient is replaced by its cell minimum.  This gives a rigorous upper
envelope without assuming monotonicity of the full composite expression.

## Finite result

```text
maximum center endpoint budget          {aggregate['maximum_center_endpoint_majorant_upper']}
maximum cell endpoint budget            {aggregate['maximum_cell_endpoint_majorant_upper']}
maximum cell/center inflation            {aggregate['maximum_cell_endpoint_inflation_upper']}
outputs below the 0.005 endpoint scale   {aggregate['outputs_within_requested_scale']} / 15
maximum parent phase Lipschitz budget    {aggregate['maximum_parent_phase_lipschitz_upper']}
maximum saddle-family Lipschitz budget   {aggregate['maximum_saddle_lipschitz_upper']}
maximum child-kernel variation           {aggregate['maximum_child_kernel_variation_upper']}
maximum multiplier relative variation    {aggregate['maximum_multiplier_relative_variation_upper']}
```

The endpoint figure uses the already certified exact output weights.  The
child, multiplier, saddle, and `fracL` arithmetic columns are not silently
added to that endpoint-only budget; they identify the remaining complete
recurrence obligations.  Each is a separate open budget.

## Boundary

{artifact['proof_boundary']}
"""


def load_cache(cache: Path, fingerprint: dict[str, Any]) -> dict[int, dict[str, Any]]:
    if not cache.is_file():
        return {}
    cached: dict[int, dict[str, Any]] = {}
    for line_number, line in enumerate(cache.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        record = json.loads(line)
        require(record.get("kind") == "hardy_block20_coefficient_neighborhood_cache_row", f"cache kind drift at line {line_number}")
        if record.get("fingerprint") != fingerprint:
            continue
        chain = int(record["chain"])
        require(chain == int(record["row"]["chain"]), f"cache chain drift at line {line_number}")
        cached[chain] = record["row"]
    return cached


def append_cache(cache: Path, fingerprint: dict[str, Any], row: dict[str, Any]) -> None:
    cache.parent.mkdir(parents=True, exist_ok=True)
    record = {
        "kind": "hardy_block20_coefficient_neighborhood_cache_row",
        "fingerprint": fingerprint,
        "chain": int(row["chain"]),
        "row": row,
    }
    with cache.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(record, separators=(",", ":")) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def build(limit: int | None, precisions: tuple[int, ...], cache: Path) -> dict[str, Any]:
    for path in (ATLAS, CENTER_ENDPOINT, WEIGHT_FIXTURE, WEIGHTS, CHECKER):
        require(path.is_file(), f"missing coefficient-neighborhood dependency: {path}")
    priority = set_low_priority()
    atlas = json.loads(ATLAS.read_text(encoding="utf-8"))
    center_endpoint = json.loads(CENTER_ENDPOINT.read_text(encoding="utf-8"))
    atlas_rows = {int(row["chain"]): row for row in atlas["rows"] if int(row["mit"]) == 2}
    center_rows = {int(row["chain"]): row for row in center_endpoint["rows"]}
    recursive = native_q.load_recursive_chains()
    require(set(atlas_rows) == set(center_rows) == set(recursive), "coefficient-neighborhood roster mismatch")
    selected_chains = sorted(recursive)
    if limit is not None:
        require(1 <= limit <= len(selected_chains), "limit outside recursive roster")
        selected_chains = selected_chains[:limit]

    fingerprint = {
        "schema": 1,
        "precisions": list(precisions),
        "atlas_sha256": file_hash(ATLAS),
        "center_endpoint_sha256": file_hash(CENTER_ENDPOINT),
        "telemetry_sha256": file_hash(native_q.TELEMETRY),
        "builder_sha256": file_hash(BUILDER),
    }
    cached = load_cache(cache, fingerprint)
    rows: list[dict[str, Any]] = []
    evaluated: dict[int, dict[str, Any]] = {}
    for chain in selected_chains:
        if chain in cached:
            row = cached[chain]
            rows.append(row)
            evaluated[chain] = {"cell_upper": Fraction(row["endpoint_column"]["cell_majorant_upper"])}
            continue
        box = coefficient_box(recursive[chain], atlas_rows[chain])
        structural = structural_transport(recursive[chain], atlas_rows[chain], box, max(precisions))
        endpoint_values = [endpoint_envelope(recursive[chain], box, dps) for dps in precisions]
        high = endpoint_values[-1]
        for low_value in endpoint_values[:-1]:
            for name in ("upper_family", "lower_family", "complete"):
                require(low_value[name].overlaps(high[name]), f"chain {chain} {name} precision nonoverlap")
        center = center_rows[chain]
        center_upper = Fraction(center["complete_majorant_upper"])
        cell_upper = upper(high["complete"])
        require(center_upper <= cell_upper, f"chain {chain} center endpoint escaped cell envelope")
        require(high["selector"] == center["selector"], f"chain {chain} lower selector drift")
        evaluated[chain] = {"cell_upper": cell_upper}

        child_delta = upper_abs(structural["child_kernel"] - structural["source_child"])
        multiplier_delta = upper_abs(structural["multiplier"] - structural["source_multiplier"])
        multiplier_center_lower = lower(abs(structural["source_multiplier"]))
        require(multiplier_center_lower > 0, "zero recurrence multiplier center")
        multiplier_relative = multiplier_delta / multiplier_center_lower
        child_coefficient_widths = [upper_abs(value - logged) for value, logged in zip(
            structural["child_coefficients"],
            [point_balls.binary128_ball(value) for value in recursive[chain]["levels"][2]["coefficients_hex"]],
        )]
        row = {
                "chain": chain,
                "sum_index": int(recursive[chain]["header"]["sum_index"]),
                "branch": int(recursive[chain]["header"]["branch"]),
                "child_length": int(recursive[chain]["levels"][2]["length"]),
                "selector": high["selector"],
                "selector_halvings": int(atlas_rows[chain]["q_cell"]["halvings"]),
                "parent_radii": {name: decimal(value) for name, value in box["radii"].items()},
                "minimum_positive_parameters": {
                    "phi2": decimal(high["phi2_low"]),
                    "endpoint_quadratic": decimal(high["endpoint_quadratic_low"]),
                    "upper_exceptional_gap": decimal(high["b_gap_low"]),
                    "lower_exceptional_gap": decimal(high["c_exceptional_gap_low"]),
                    "xi_lower": decimal(high["xi_low"]),
                    "xi_upper": decimal(high["xi_high"]),
                },
                "phase_column": {
                    "parent_phase_lipschitz_upper": decimal(upper(structural["phase_lipschitz"])),
                },
                "child_column": {
                    "coefficient_variation_uppers": [decimal(value) for value in child_coefficient_widths],
                    "kernel_variation_upper": decimal(child_delta),
                    "source_kernel_contained": True,
                },
                "multiplier_column": {
                    "variation_upper": decimal(multiplier_delta),
                    "relative_variation_upper": decimal(multiplier_relative),
                    "source_multiplier_contained": True,
                },
                "saddle_column": {
                    "equation_69_lipschitz_upper": decimal(upper(structural["saddle_lipschitz"])),
                },
                "endpoint_column": {
                    "center_majorant_upper": decimal(center_upper),
                    "cell_majorant_lower": decimal(lower(high["complete"])),
                    "cell_majorant_upper": decimal(cell_upper),
                    "cell_to_center_inflation_upper": decimal(cell_upper / center_upper),
                    "precision_overlap": True,
                },
                "arithmetic_column": {
                    "fracL_radius": decimal(box["radii"]["fracL"]),
                    "status": "open_not_allocated_to_endpoint_budget",
                },
            }
        rows.append(row)
        append_cache(cache, fingerprint, row)

    weights = transport.load_weights()
    output_rows: list[dict[str, Any]] = []
    if len(selected_chains) == 374:
        for output_index in range(1, 16):
            center_total = Fraction(0)
            cell_total = Fraction(0)
            for row in rows:
                weight = weights[(int(row["sum_index"]), output_index, int(row["branch"]))]
                factor = weight["amplitude_fraction"] * weight["phase_abs_upper"]
                center_total += factor * Fraction(row["endpoint_column"]["center_majorant_upper"])
                cell_total += factor * Fraction(row["endpoint_column"]["cell_majorant_upper"])
            output_rows.append(
                {
                    "output_index": output_index,
                    "output_label": transport.output_label(output_index),
                    "center_endpoint_majorant_upper": decimal(center_total),
                    "cell_endpoint_majorant_upper": decimal(cell_total),
                    "cell_to_center_inflation_upper": decimal(cell_total / center_total),
                    "cell_endpoint_to_requested_scale_ratio": decimal(cell_total / REQUESTED_ERROR_SCALE),
                    "within_requested_scale": cell_total < REQUESTED_ERROR_SCALE,
                }
            )

    def maximum(row_path: tuple[str, ...]) -> tuple[Fraction, int]:
        values: list[tuple[Fraction, int]] = []
        for row in rows:
            value: Any = row
            for key in row_path:
                value = value[key]
            values.append((Fraction(value), int(row["chain"])))
        return max(values)

    max_center = maximum(("endpoint_column", "center_majorant_upper"))
    max_cell = maximum(("endpoint_column", "cell_majorant_upper"))
    max_inflation = maximum(("endpoint_column", "cell_to_center_inflation_upper"))
    max_phase = maximum(("phase_column", "parent_phase_lipschitz_upper"))
    max_child = maximum(("child_column", "kernel_variation_upper"))
    max_multiplier = maximum(("multiplier_column", "relative_variation_upper"))
    max_saddle = maximum(("saddle_column", "equation_69_lipschitz_upper"))
    full = len(selected_chains) == 374
    outputs_within = sum(bool(row["within_requested_scale"]) for row in output_rows)
    status = (
        "rigorous_full_cell_endpoint_transport_closes_finite_block20_scale_other_columns_open"
        if full and outputs_within == 15
        else "rigorous_full_cell_endpoint_transport_exceeds_finite_block20_scale_other_columns_open"
        if full
        else "rigorous_coefficient_neighborhood_transport_pilot"
    )
    aggregate = {
        "cell_count": len(rows),
        "full_roster": full,
        "precision_ladder_decimal_digits": list(precisions),
        "precision_overlap_count": len(rows),
        "w3_count": sum(row["selector"] == "W3" for row in rows),
        "w4_count": sum(row["selector"] == "W4" for row in rows),
        "maximum_center_endpoint_majorant_upper": decimal(max_center[0]),
        "maximum_center_endpoint_majorant_witness": max_center[1],
        "maximum_cell_endpoint_majorant_upper": decimal(max_cell[0]),
        "maximum_cell_endpoint_majorant_witness": max_cell[1],
        "maximum_cell_endpoint_inflation_upper": decimal(max_inflation[0]),
        "maximum_cell_endpoint_inflation_witness": max_inflation[1],
        "maximum_parent_phase_lipschitz_upper": decimal(max_phase[0]),
        "maximum_parent_phase_lipschitz_witness": max_phase[1],
        "maximum_child_kernel_variation_upper": decimal(max_child[0]),
        "maximum_child_kernel_variation_witness": max_child[1],
        "maximum_multiplier_relative_variation_upper": decimal(max_multiplier[0]),
        "maximum_multiplier_relative_variation_witness": max_multiplier[1],
        "maximum_saddle_lipschitz_upper": decimal(max_saddle[0]),
        "maximum_saddle_lipschitz_witness": max_saddle[1],
        "requested_endpoint_error_scale": decimal(REQUESTED_ERROR_SCALE),
        "outputs_within_requested_scale": outputs_within,
        "maximum_transported_cell_endpoint_majorant_upper": (
            decimal(max(Fraction(row["cell_endpoint_majorant_upper"]) for row in output_rows))
            if output_rows
            else None
        ),
    }
    return {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_coefficient_neighborhood_transport_gate",
        "status": status,
        "scope": {
            "block": 20,
            "selected_recursive_calls": len(rows),
            "full_recursive_calls": 374,
            "cell_coordinates": ["Phi1", "xr=2*Phi2", "Phi3", "fracL"],
            "process_priority": priority,
            "resumable_cache_rows_reused": sum(chain in cached for chain in selected_chains),
        },
        "transport_theorem": {
            "phase": "|Delta phase(n)| <= tpp_src*(r1*n+r2*n^2+r3*n^3)",
            "saddle": "|Delta sum I_n| <= tpp_src*J*(r1*N^2/2+r2*N^3/3+r3*N^4/4)",
            "child": "exact interval cubic Legendre map with fixed n1, ix, parity sign, orientation, and nonzero centered denominator",
            "multiplier": "source-anchored M=M0*sqrt(y0/y)*exp(i*tpp*(r0-r0_center))",
            "endpoint": "positive numerator maxima, gap minima, and retained-decay coefficient minima give a pointwise integral upper envelope",
            "arithmetic": "fracL radius recorded but PSI/ERF and recurrence arithmetic remain a separate open budget",
        },
        "rows": rows,
        "transported_outputs": output_rows,
        "aggregate": aggregate,
        "dependencies": {
            "selector_atlas": {"path": relative(ATLAS), "sha256": file_hash(ATLAS)},
            "center_endpoint": {"path": relative(CENTER_ENDPOINT), "sha256": file_hash(CENTER_ENDPOINT)},
            "weight_fixture": {"path": relative(WEIGHT_FIXTURE), "sha256": file_hash(WEIGHT_FIXTURE)},
            "weights": {"path": relative(WEIGHTS), "sha256": file_hash(WEIGHTS)},
            "telemetry": {"path": relative(native_q.TELEMETRY), "sha256": file_hash(native_q.TELEMETRY)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "cache": {"path": relative(cache), "sha256": file_hash(cache)},
        },
        "next_obligation": (
            "Use the separate child, multiplier, saddle, and fracL arithmetic columns to choose accuracy-refined "
            "subcells and certify the complete recurrence and Legendre-tail budget without borrowing the endpoint allocation."
        ),
        "proof_boundary": (
            "Rigorous coefficient-cell transport for the finite block-20 recursive roster. The fixed-weight endpoint "
            "column is separate from child-kernel, multiplier, saddle-family, PSI/ERF arithmetic, Legendre-tail, "
            "other-block, outer-Hardy, and height-uniform obligations. No complete evaluator theorem, Lambda <= 0, "
            "PF-infinity, RH, or prize-level conclusion follows."
        ),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--precisions", type=int, nargs="+", default=list(PRECISIONS))
    parser.add_argument("--output", type=Path, default=RESULT)
    parser.add_argument("--note", type=Path, default=NOTE)
    parser.add_argument("--cache", type=Path, default=CACHE)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    precisions = tuple(int(value) for value in args.precisions)
    require(precisions and list(precisions) == sorted(set(precisions)), "precision ladder must be strictly increasing")
    cache = args.cache if args.cache.is_absolute() else REPO_ROOT / args.cache
    artifact = build(args.limit, precisions, cache)
    output = args.output if args.output.is_absolute() else REPO_ROOT / args.output
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    output.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    output_tmp = output.with_suffix(output.suffix + ".tmp")
    note_tmp = note.with_suffix(note.suffix + ".tmp")
    output_tmp.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    note_tmp.write_text(build_note(artifact), encoding="utf-8")
    output_tmp.replace(output)
    note_tmp.replace(note)
    print(
        "built Hardy block-20 coefficient-neighborhood transport: "
        f"{artifact['aggregate']['cell_count']} cells, "
        f"{artifact['aggregate']['outputs_within_requested_scale']}/15 endpoint outputs below 0.005"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
