#!/usr/bin/env python3
"""Enclose the signed corrected-model bridge on two later-block coefficient cells."""

from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_full_gate as log_full
import jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_pilot_gate as log_pilot
import jensen_window_pf_newman_c1_hardy_block20_coefficient_neighborhood_transport_gate as block20_cells
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
import jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_pilot_gate as point_bridge
import jensen_window_pf_newman_c1_hardy_t1e10_crossblock_structural_atlas_gate as crossblock
from flint import acb, arb, ctx


TELEMETRY = crossblock.TELEMETRY
BLOCK21_POINTS = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_full_gate.json"
LATER_POINTS = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_blocks22_28_source_to_corrected_model_bridge_full_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_pilot_gate.json"
CACHE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_pilot_gate_cache.jsonl"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_pilot_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_pilot_gate.py"
WITNESS_CHAINS = (780, 1456)
PRECISIONS = (70, 110)
TAIL_TARGET = Fraction(1, 10**30)
MAXIMUM_CELL_TO_POINT_INFLATION = Fraction(5, 4)
FINITE_EXACT_POINT_MARGIN = Fraction(
    "4.55601112596361510673046545569964374799791938123474E-3"
)
CPU_PARK_THRESHOLD = 75.0
I = acb(0, 1)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def fraction_record(value: Fraction) -> dict[str, Any]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": decimal(value),
    }


def exact_fraction(record: dict[str, Any]) -> Fraction:
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def upper(value: arb) -> Fraction:
    return cells.bound_fraction(value.upper())


def lower(value: arb) -> Fraction:
    return cells.bound_fraction(value.lower())


def upper_abs(value: acb | arb) -> Fraction:
    return upper(abs(value))


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


def cpu_percent() -> float | None:
    try:
        import psutil

        return float(psutil.cpu_percent(interval=1.0))
    except Exception:  # pragma: no cover - optional monitor
        return None


def load_witnesses() -> dict[int, dict[str, Any]]:
    selected = set(WITNESS_CHAINS)
    chains: dict[int, dict[str, Any]] = {}
    with TELEMETRY.open(encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            chain = int(record.get("chain", 0))
            if chain not in selected:
                continue
            if record["type"] == "chain":
                chains[chain] = {"header": record, "levels": {}, "q_terms": {}, "steps": {}}
            elif record["type"] == "level":
                chains[chain]["levels"][int(record["level"])] = record
            elif record["type"] == "q_terms":
                chains[chain]["q_terms"][int(record["nit"])] = record
            elif record["type"] == "recurrence":
                chains[chain]["steps"][int(record["nit"])] = record
    require(set(chains) == selected, "coefficient-cell pilot witness roster missing")
    for chain, data in chains.items():
        require(int(data["header"]["mit"]) == 2, f"chain {chain} is not recursive")
        require(set(data["levels"]) >= {1, 2}, f"chain {chain} level roster missing")
        require(int(data["levels"][2]["length"]) == 1, f"chain {chain} child length drift")
    return chains


def load_point_rows() -> dict[int, dict[str, Any]]:
    rows: dict[int, dict[str, Any]] = {}
    for path in (BLOCK21_POINTS, LATER_POINTS):
        artifact = json.loads(path.read_text(encoding="utf-8"))
        require(bool(artifact["passed"]), f"failed point dependency: {path.name}")
        for row in artifact["rows"]:
            chain = int(row["chain"])
            if chain in WITNESS_CHAINS:
                require(chain not in rows, f"duplicate point witness {chain}")
                rows[chain] = row
    require(set(rows) == set(WITNESS_CHAINS), "point witness rows missing")
    return rows


def algebra_audit() -> int:
    checks = 0
    fixtures = (
        (1, 2, 3, 5, 7, 11, 13, 17, 19),
        (-2, 3, -5, 7, -11, 13, -17, 19, -23),
        (Fraction(1, 3), Fraction(-2, 5), Fraction(3, 7), Fraction(5, 11),
         Fraction(-7, 13), Fraction(11, 17), Fraction(13, 19), Fraction(-17, 23), Fraction(19, 29)),
    )
    for ps, isrc, mc, qs, t5, qextra, pp, ip, qp in fixtures:
        ps, isrc, mc, qs, t5, qextra, pp, ip, qp = map(
            Fraction, (ps, isrc, mc, qs, t5, qextra, pp, ip, qp)
        )
        qdrop = qs - t5 - qextra
        dsrc = ps - mc - qs
        d1 = dsrc - (isrc - mc - t5)
        d2 = d1 + qextra
        dcorr = d2 + (pp - ps) - (ip - isrc) - (qp - qdrop)
        require(d1 == ps - isrc - (qs - t5), "signed W1 bridge algebra failure")
        checks += 1
        require(d2 == ps - isrc - qdrop, "signed t2 bridge algebra failure")
        checks += 1
        require(dcorr == pp - ip - qp, "signed corrected bridge algebra failure")
        checks += 1
    return checks


def base_selector_cell(chain: int, data: dict[str, Any]) -> dict[str, Any]:
    return cells.selector_cell(
        chain,
        data["header"],
        data["levels"][1],
        data["levels"][2],
    )


def refine_selector_cell(base: dict[str, Any], extra_halvings: int) -> dict[str, Any]:
    require(extra_halvings >= 0, "negative extra cell halvings")
    refined = deepcopy(base)
    factor = Fraction(1, 2**extra_halvings)
    for name, record in refined["radii"].items():
        radius = exact_fraction(record["exact"]) * factor
        record["exact"] = fraction_record(radius)
        record["decimal"] = decimal(radius)
    total_halvings = int(base["halvings"]) + extra_halvings
    refined["halvings"] = total_halvings
    refined["scale"] = fraction_record(Fraction(1, 2**total_halvings))
    refined["inherited_selector_cell_halvings"] = int(base["halvings"])
    refined["transport_extra_halvings"] = extra_halvings
    return refined


def coefficient_box(data: dict[str, Any], cell: dict[str, Any]) -> dict[str, Any]:
    return block20_cells.coefficient_box(data, {"q_cell": cell})


def initial_extra_halvings(box: dict[str, Any], data: dict[str, Any], point_upper: Fraction) -> int:
    radii = box["radii"]
    length = int(data["levels"][1]["length"])
    # The integrated phase variation is a rigorous first-order envelope for
    # |exp(i*f)-exp(i*g)| after integrating over y in [0,L].
    integrated_phase = 2 * Fraction(355, 113) * (
        radii["a1"] * length**2 / 2
        + (radii["y"] / 2) * length**3 / 3
        + radii["a3"] * length**4 / 4
    )
    target = point_upper / 32
    if integrated_phase <= target:
        return 0
    return max(0, math.ceil(math.log2(float(integrated_phase / target))))


def coefficient_balls(box: dict[str, Any]) -> tuple[arb, arb, arb]:
    center = box["center"]
    radii = box["radii"]
    return (
        cells.interval_ball(center["a1"], radii["a1"]),
        cells.interval_ball(center["y"], radii["y"]),
        cells.interval_ball(center["a3"], radii["a3"]),
    )


def cell_sector_contract(box: dict[str, Any], data: dict[str, Any]) -> dict[str, Any]:
    bounds = box["bounds"]
    a1_low, a1_high = bounds["a1"]
    y_low, y_high = bounds["y"]
    a3_low, a3_high = bounds["a3"]
    length = int(data["levels"][1]["length"])
    child_length = int(data["levels"][2]["length"])
    xi_low = a1_low + y_low * length + 3 * a3_low * length**2
    xi_high = a1_high + y_high * length + 3 * a3_high * length**2
    require(Fraction(child_length) < xi_low < xi_high < Fraction(child_length + 1), "dual selector escaped refined cell")
    curvature_low = min(y_low, y_low + 6 * a3_low * length)
    require(curvature_low > 0, "cell curvature floor is nonpositive")
    cubic_sign = 1 if box["center"]["a3"] > 0 else -1
    cubic_low = min(abs(a3_low), abs(a3_high))
    require(cubic_low > 0, "cell cubic coefficient crosses zero")
    if cubic_sign < 0:
        b = {
            "angle": "5*pi/6",
            "linear": (Fraction(child_length + 1) - xi_high) / 2,
            "quadratic": None,
            "cubic": cubic_low,
        }
        c = {"angle": "pi/2", "linear": 1 + a1_low, "quadratic": Fraction(0), "cubic": cubic_low}
    else:
        b = {"angle": "pi/2", "linear": Fraction(child_length + 1) - xi_high, "quadratic": Fraction(0), "cubic": cubic_low}
        c = {"angle": "pi/6", "linear": (1 + a1_low) / 2, "quadratic": None, "cubic": cubic_low}
    require(b["linear"] > 0 and c["linear"] > 0, "cell ray linear decay is nonpositive")
    return {
        "cubic_sign": cubic_sign,
        "curvature_low": curvature_low,
        "b": b,
        "c": c,
    }


def ray_integral_cell(
    box: dict[str, Any],
    data: dict[str, Any],
    sector: dict[str, Any],
    endpoint: int,
    family: str,
    dps: int,
    cutoff: Fraction,
) -> dict[str, Any]:
    ctx.dps = dps
    ctx.threads = 1
    a1, y, a3 = coefficient_balls(box)
    a2 = y / 2
    contract = sector[family]
    direction, sine = log_pilot.direction(contract["angle"])
    phase_sign = -1 if family == "b" else 1
    start = int(data["levels"][2]["length"]) + 1 if family == "b" else 1
    endpoint_ball = arb(endpoint)
    pi = arb.pi()

    def phase(z: acb) -> acb:
        return a1 * z + a2 * z**2 + a3 * z**3

    def derivative(z: acb) -> acb:
        return a1 + 2 * a2 * z + 3 * a3 * z**2

    def make_integrand(hypergeometric: bool):
        def integrand(radius: acb, _analytic: bool) -> acb:
            u = direction * radius
            z = endpoint_ball + u
            exponential = (2 * pi * I * u).exp()
            kernel = log_pilot.kernel(exponential, start, hypergeometric)
            return derivative(z) * (phase_sign * 2 * pi * I * phase(z)).exp() * kernel * direction

        return integrand

    breaks = log_full.breaks_for_cutoff(cutoff)
    tolerance = arb(f"1e-{dps - 25}")
    central = acb(0)
    for lo, hi in zip(breaks[:-1], breaks[1:]):
        central += acb.integral(
            make_integrand(lo >= Fraction(1, 2)),
            acb(cells.fraction_ball(lo)),
            acb(cells.fraction_ball(hi)),
            abs_tol=tolerance,
            rel_tol=tolerance,
            deg_limit=64,
            eval_limit=200000,
            depth_limit=64,
        )

    epsilon = cells.fraction_ball(log_pilot.ORIGIN_CUTOFF)
    far = cells.fraction_ball(cutoff)
    angular_decay = 2 * pi * cells.fraction_ball(sine)
    endpoint_derivative = derivative(endpoint_ball)
    endpoint_second = 2 * a2 + 6 * a3 * endpoint_ball
    amplitude_origin = abs(endpoint_derivative) + abs(endpoint_second) * epsilon + 3 * abs(a3) * epsilon**2
    kernel_origin_integral = (angular_decay * epsilon).exp() * (
        epsilon * (1 - (angular_decay * epsilon).log()) + angular_decay * epsilon**2 / 2
    )
    origin_radius = abs(direction) * amplitude_origin * kernel_origin_integral

    linear_decay = cells.fraction_ball(contract["linear"])
    if contract["quadratic"] is None:
        quadratic_decay = arb(3).sqrt() * cells.fraction_ball(sector["curvature_low"]) / 4
    else:
        quadratic_decay = cells.fraction_ball(contract["quadratic"])
    cubic_decay = cells.fraction_ball(contract["cubic"])
    decay_rate = 2 * pi * (linear_decay + quadratic_decay * far + cubic_decay * far**2)
    require(lower(decay_rate) > 0, "cell ray tail decay is nonpositive")
    inverse = 1 / decay_rate
    p0, p1, p2 = abs(endpoint_derivative), abs(endpoint_second), 3 * abs(a3)
    polynomial_tail = (-decay_rate * far).exp() * (
        p0 * inverse
        + p1 * (far * inverse + inverse**2)
        + p2 * (far**2 * inverse + 2 * far * inverse**2 + 2 * inverse**3)
    )
    denominator = start * (1 - (-angular_decay * far).exp())
    tail_radius = abs(direction) * polynomial_tail / denominator
    origin_upper, tail_upper = upper(origin_radius), upper(tail_radius)
    return {
        "total": central + log_pilot.error_box(origin_upper + tail_upper),
        "origin_upper": origin_upper,
        "tail_upper": tail_upper,
    }


def finite_zero_mode_cell(box: dict[str, Any], data: dict[str, Any], dps: int) -> acb:
    ctx.dps = dps
    ctx.threads = 1
    a1, y, a3 = coefficient_balls(box)
    a2 = y / 2
    length = int(data["levels"][1]["length"])
    pi = arb.pi()
    tolerance = arb(f"1e-{dps - 25}")

    def integrand(value: acb, _analytic: bool) -> acb:
        phase = a1 * value + a2 * value**2 + a3 * value**3
        return (-2 * pi * I * phase).exp()

    total = acb(0)
    lo = 0
    while lo < length:
        hi = min(length, lo + 13)
        total += acb.integral(
            integrand,
            acb(lo),
            acb(hi),
            abs_tol=tolerance,
            rel_tol=tolerance,
            deg_limit=64,
            eval_limit=200000,
            depth_limit=64,
        )
        lo = hi
    return total


def exact_nonsaddle_cell(box: dict[str, Any], data: dict[str, Any], dps: int) -> dict[str, Any]:
    ctx.dps = dps
    a1, y, a3 = coefficient_balls(box)
    a2 = y / 2
    length = int(data["levels"][1]["length"])
    child_length = int(data["levels"][2]["length"])
    sector = cell_sector_contract(box, data)
    cutoff = log_full.INITIAL_CUTOFF
    while True:
        rays = [
            ray_integral_cell(box, data, sector, endpoint, family, dps, cutoff)
            for family in ("b", "c")
            for endpoint in (0, length)
        ]
        if max(ray["tail_upper"] for ray in rays) < TAIL_TARGET:
            break
        cutoff *= 2
        require(cutoff <= log_full.MAXIMUM_CUTOFF, "cell ray cutoff exceeded safety ceiling")
    b_source = rays[0]["total"] - rays[1]["total"]
    c_source = (rays[3]["total"] - rays[2]["total"]).conjugate()
    pi = arb.pi()
    endpoint_phase = a1 * length + a2 * length**2 + a3 * length**3
    endpoint = (-2 * pi * I * endpoint_phase).exp()
    endpoint_half = (1 + endpoint) / 2
    endpoint_harmonic = I * (endpoint - 1) * cells.fraction_ball(log_pilot.harmonic(child_length)) / (2 * pi)
    zero_mode = finite_zero_mode_cell(box, data, dps) if box["bounds"]["a1"][0] > 0 else acb(0)
    return {
        "value": endpoint_half + endpoint_harmonic + b_source + c_source + zero_mode,
        "maximum_ray_origin_upper": max(ray["origin_upper"] for ray in rays),
        "maximum_ray_tail_upper": max(ray["tail_upper"] for ray in rays),
        "selected_ray_cutoff": cutoff,
    }


def paper_components_cell(box: dict[str, Any], data: dict[str, Any], dps: int) -> dict[str, Any]:
    ctx.dps = dps
    a1, y, a3 = coefficient_balls(box)
    a2 = y / 2
    length = int(data["levels"][1]["length"])
    child_length = int(data["levels"][2]["length"])
    xi = a1 + y * length + 3 * a3 * length**2
    require(lower(xi) > child_length and upper(xi) < child_length + 1, "paper dual selector escaped cell")
    delta = child_length + 1 - xi
    sqrt_y = y.sqrt()
    pi = arb.pi()
    root_two = arb(2).sqrt()
    epi_minus = acb(1 / root_two, -1 / root_two)
    endpoint_phase = a1 * length + a2 * length**2 + a3 * length**3
    endpoint = acb(0, -2 * pi * endpoint_phase).exp()
    internal_endpoint = endpoint.conjugate()

    def fresnel_boundary(argument: arb) -> acb:
        z = epi_minus * (pi / y).sqrt() * argument
        return I * epi_minus * acb(0, -pi * argument**2 / y).exp() * (1 - z.erf()) / (2 * sqrt_y)

    w2 = (
        I
        * internal_endpoint
        * (
            epi_minus
            * acb(0, -pi * delta**2 / y).exp()
            * (1 - (epi_minus * (pi / y).sqrt() * delta).erf())
            / (2 * sqrt_y)
            - 1 / (2 * pi * delta)
        )
    ).conjugate()
    if box["bounds"]["a1"][1] < 0:
        selector = "W3"
        argument = a1 + 1
        w34 = (fresnel_boundary(argument) - I / (2 * pi * argument)).conjugate()
    else:
        require(box["bounds"]["a1"][0] > 0, "paper W3/W4 selector changes inside cell")
        selector = "W4"
        w34 = (fresnel_boundary(a1) - I * internal_endpoint / (2 * pi * xi)).conjugate()
    w5 = (
        I
        / (2 * pi)
        * (
            internal_endpoint * ((xi + 1).digamma() - delta.digamma())
            - ((a1 + 1).digamma() - (child_length + 1 - a1).digamma())
        )
    ).conjugate()
    half_sum = ((1 + internal_endpoint) / 2).conjugate()
    return {"selector": selector, "value": w2 + w34 + w5 + half_sum}


def dyadic(pair: list[int]) -> Fraction:
    mantissa, exponent = int(pair[0]), int(pair[1])
    return Fraction(mantissa * 2**exponent) if exponent >= 0 else Fraction(mantissa, 2 ** (-exponent))


def interval_from_record(record: dict[str, Any]) -> arb:
    lo, hi = dyadic(record["lower_dyadic"]), dyadic(record["upper_dyadic"])
    require(lo < hi, "degenerate stored point interval")
    return cells.interval_ball((lo + hi) / 2, (hi - lo) / 2)


def point_correction_box(row: dict[str, Any]) -> acb:
    record = row["exact_minus_paper"]
    return acb(interval_from_record(record["real"]), interval_from_record(record["imag"]))


def evaluate_cell(box: dict[str, Any], data: dict[str, Any], dps: int) -> dict[str, Any]:
    exact = exact_nonsaddle_cell(box, data, dps)
    paper = paper_components_cell(box, data, dps)
    require(paper["selector"] in ("W3", "W4"), "cell selector drift")
    correction = exact["value"] - paper["value"]
    return {
        "selector": paper["selector"],
        "correction": correction,
        "exact_nonsaddle": exact["value"],
        "paper_nonsaddle": paper["value"],
        "maximum_ray_origin_upper": exact["maximum_ray_origin_upper"],
        "maximum_ray_tail_upper": exact["maximum_ray_tail_upper"],
        "selected_ray_cutoff": exact["selected_ray_cutoff"],
    }


def cache_fingerprint() -> str:
    payload = {
        "builder": file_hash(BUILDER),
        "checker": file_hash(CHECKER),
        "telemetry": file_hash(TELEMETRY),
        "block21_points": file_hash(BLOCK21_POINTS),
        "later_points": file_hash(LATER_POINTS),
        "selector_builder": file_hash(Path(cells.__file__).resolve()),
        "witnesses": list(WITNESS_CHAINS),
        "precisions": list(PRECISIONS),
        "maximum_inflation": [MAXIMUM_CELL_TO_POINT_INFLATION.numerator, MAXIMUM_CELL_TO_POINT_INFLATION.denominator],
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def load_cache(fingerprint: str) -> dict[int, dict[str, Any]]:
    if not CACHE.exists():
        CACHE.parent.mkdir(parents=True, exist_ok=True)
        CACHE.write_text(json.dumps({"kind": "later_cell_bridge_pilot_cache", "fingerprint": fingerprint}, sort_keys=True) + "\n", encoding="utf-8")
        return {}
    lines = [line for line in CACHE.read_text(encoding="utf-8").splitlines() if line.strip()]
    require(lines, "empty coefficient-cell pilot cache")
    header = json.loads(lines[0])
    if header.get("fingerprint") != fingerprint:
        stale = CACHE.with_name(f"{CACHE.stem}.stale.{header.get('fingerprint', 'unknown')[:12]}{CACHE.suffix}")
        require(not stale.exists(), f"stale cache destination exists: {stale}")
        CACHE.replace(stale)
        CACHE.write_text(json.dumps({"kind": "later_cell_bridge_pilot_cache", "fingerprint": fingerprint}, sort_keys=True) + "\n", encoding="utf-8")
        return {}
    rows: dict[int, dict[str, Any]] = {}
    for line in lines[1:]:
        row = json.loads(line)
        chain = int(row["chain"])
        require(chain not in rows, f"duplicate cached witness {chain}")
        rows[chain] = row
    return rows


def append_cache(row: dict[str, Any]) -> None:
    with CACHE.open("a", encoding="utf-8", buffering=1) as handle:
        handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def build_row(chain: int, data: dict[str, Any], point_row: dict[str, Any]) -> dict[str, Any]:
    base = base_selector_cell(chain, data)
    base_box = coefficient_box(data, base)
    point_upper = Fraction(point_row["exact_correction_magnitude_upper"])
    extra = initial_extra_halvings(base_box, data, point_upper)
    point_box = point_correction_box(point_row)
    accepted: tuple[dict[str, Any], dict[str, Any], dict[str, Any], int] | None = None
    for extra_halvings in range(extra, extra + 13):
        cell = refine_selector_cell(base, extra_halvings)
        box = coefficient_box(data, cell)
        low = evaluate_cell(box, data, PRECISIONS[0])
        cell_upper = upper_abs(low["correction"])
        if cell_upper > MAXIMUM_CELL_TO_POINT_INFLATION * point_upper:
            continue
        high = evaluate_cell(box, data, PRECISIONS[1])
        if not low["correction"].overlaps(high["correction"]):
            continue
        if not high["correction"].overlaps(point_box):
            continue
        high_upper = upper_abs(high["correction"])
        if high_upper > MAXIMUM_CELL_TO_POINT_INFLATION * point_upper:
            continue
        accepted = (cell, low, high, extra_halvings)
        break
    require(accepted is not None, f"chain {chain} failed to close a refined signed cell")
    cell, low, high, extra_halvings = accepted
    box = coefficient_box(data, cell)
    high_upper = upper_abs(high["correction"])
    retained = Fraction(point_row["retained_local_majorant_upper"])
    weight = Fraction(point_row["outer_weight_factor"])
    inflation = high_upper / point_upper
    transported = weight * high_upper
    point_transported = weight * point_upper
    return {
        "chain": chain,
        "block": int(data["header"]["block"]),
        "sum_index": int(data["header"]["sum_index"]),
        "branch": int(data["header"]["branch"]),
        "parent_length": int(data["levels"][1]["length"]),
        "child_length": int(data["levels"][2]["length"]),
        "selector": high["selector"],
        "phi3_sign": "positive" if box["center"]["a3"] > 0 else "negative",
        "base_selector_cell": base,
        "refined_cell": cell,
        "transport_extra_halvings": extra_halvings,
        "radii_nonzero": all(exact_fraction(record["exact"]) > 0 for record in cell["radii"].values()),
        "low_precision_correction": point_balls.acb_record(low["correction"], 45),
        "high_precision_correction": point_balls.acb_record(high["correction"], 45),
        "precision_overlap": low["correction"].overlaps(high["correction"]),
        "point_correction_overlap": high["correction"].overlaps(point_box),
        "point_correction_magnitude_upper": decimal(point_upper),
        "cell_correction_magnitude_upper": decimal(high_upper),
        "cell_to_point_inflation_upper": decimal(inflation),
        "maximum_allowed_inflation": decimal(MAXIMUM_CELL_TO_POINT_INFLATION),
        "retained_local_majorant_upper": decimal(retained),
        "cell_correction_within_retained_majorant": high_upper < retained,
        "outer_weight_factor": decimal(weight),
        "point_transported_correction_upper": decimal(point_transported),
        "cell_transported_correction_upper": decimal(transported),
        "transported_inflation_upper": decimal(transported - point_transported),
        "maximum_ray_origin_upper": decimal(high["maximum_ray_origin_upper"]),
        "maximum_ray_tail_upper": decimal(high["maximum_ray_tail_upper"]),
        "selected_ray_cutoff": decimal(high["selected_ray_cutoff"]),
        "signed_bridge_closed_on_analytic_cell": True,
    }


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy t=1e10 later-block signed coefficient-cell bridge pilot

Date: 2026-08-09

Status: rigorous two-witness corrected-model coefficient-cell pilot; not a full recurrence, height-uniform theorem, or RH proof

The worst exact/retained-ratio witnesses from block 21 and the later recursive
blocks are promoted from exact coefficient points to nonzero-radius,
selector-stable boxes.  The source bridge is kept signed algebraically:

```text
Dsrc = Psrc - MC - qsrc
D1   = Dsrc - (I69src - MC - t5src)
D2   = D1 + qextra = Psrc - I69src - Qsrc,drop
Dcorr = D2 + (Ppi-Psrc) - (I69pi-I69src) - (Qpaper-Qsrc,drop)
      = Ppi - I69pi - Qpaper = Qexact - Qpaper.
```

The last expression is evaluated as one complex Arb interval on each cell;
the exact contour and paper terms are not separately converted to absolute
majorants before subtraction.  Both 70- and 110-digit enclosures overlap the
saved exact-point correction.

Witness cells closed: `{aggregate['closed_cell_count']}/{aggregate['cell_count']}`

Maximum cell/point correction inflation: `{aggregate['maximum_cell_to_point_inflation_upper']}`

Maximum transported two-cell inflation: `{aggregate['total_transported_inflation_upper']}`

Remaining finite exact-point margin: `{aggregate['finite_exact_point_margin']}`

Boundary: these are analytic corrected-model cells.  The theorem does not yet
transport source binary128 rounding, the complete cubic recurrence tail, all
1,040 later recursive calls, the outer Hardy representation/remainder, or any
height-uniform conclusion.  Those remain explicit obligations.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime-limit-seconds", type=float, default=10800.0)
    args = parser.parse_args()
    for path in (TELEMETRY, BLOCK21_POINTS, LATER_POINTS, CHECKER):
        require(path.is_file(), f"missing coefficient-cell pilot dependency: {path}")
    started = time.monotonic()
    priority = set_low_priority()
    ctx.threads = 1
    witnesses = load_witnesses()
    points = load_point_rows()
    algebra_checks = algebra_audit()
    fingerprint = cache_fingerprint()
    cached = load_cache(fingerprint)
    require(set(cached).issubset(WITNESS_CHAINS), "unexpected coefficient-cell cache witness")
    high_cpu_streak = 0
    for chain in WITNESS_CHAINS:
        if chain in cached:
            continue
        if time.monotonic() - started >= args.runtime_limit_seconds:
            print(f"parked later coefficient-cell pilot after {len(cached)}/{len(WITNESS_CHAINS)} witnesses")
            return 75
        row = build_row(chain, witnesses[chain], points[chain])
        append_cache(row)
        cached[chain] = row
        print(f"closed later coefficient-cell witness {chain} ({len(cached)}/{len(WITNESS_CHAINS)})", flush=True)
        sampled = cpu_percent()
        if sampled is not None:
            high_cpu_streak = high_cpu_streak + 1 if sampled > CPU_PARK_THRESHOLD else 0
            if high_cpu_streak >= 2 and len(cached) < len(WITNESS_CHAINS):
                print(f"parked later coefficient-cell pilot after {len(cached)}/{len(WITNESS_CHAINS)} witnesses (sustained CPU)")
                return 75
    rows = [cached[chain] for chain in WITNESS_CHAINS]
    maximum_inflation = max(Fraction(row["cell_to_point_inflation_upper"]) for row in rows)
    total_transported_inflation = sum(
        (Fraction(row["transported_inflation_upper"]) for row in rows), Fraction(0)
    )
    aggregate = {
        "cell_count": len(rows),
        "closed_cell_count": sum(bool(row["signed_bridge_closed_on_analytic_cell"]) for row in rows),
        "precision_overlap_count": sum(bool(row["precision_overlap"]) for row in rows),
        "point_overlap_count": sum(bool(row["point_correction_overlap"]) for row in rows),
        "nonzero_radius_cell_count": sum(bool(row["radii_nonzero"]) for row in rows),
        "within_retained_majorant_count": sum(bool(row["cell_correction_within_retained_majorant"]) for row in rows),
        "exact_algebra_checks": algebra_checks,
        "maximum_cell_to_point_inflation_upper": decimal(maximum_inflation),
        "maximum_allowed_inflation": decimal(MAXIMUM_CELL_TO_POINT_INFLATION),
        "total_transported_inflation_upper": decimal(total_transported_inflation),
        "finite_exact_point_margin": decimal(FINITE_EXACT_POINT_MARGIN),
        "transported_inflation_within_finite_margin": total_transported_inflation < FINITE_EXACT_POINT_MARGIN,
        "remaining_later_recursive_calls": 1038,
    }
    passed = (
        aggregate["closed_cell_count"] == len(WITNESS_CHAINS)
        and aggregate["precision_overlap_count"] == len(WITNESS_CHAINS)
        and aggregate["point_overlap_count"] == len(WITNESS_CHAINS)
        and aggregate["nonzero_radius_cell_count"] == len(WITNESS_CHAINS)
        and aggregate["within_retained_majorant_count"] == len(WITNESS_CHAINS)
        and maximum_inflation <= MAXIMUM_CELL_TO_POINT_INFLATION
        and aggregate["transported_inflation_within_finite_margin"]
    )
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_pilot_gate",
        "status": "two_worst_ratio_later_recursive_signed_coefficient_cells_closed" if passed else "later_recursive_signed_coefficient_cell_pilot_falsified",
        "scope": "Block-21 chain 780 and block-23 chain 1456 at the saved t=1e10 fixture.",
        "passed": passed,
        "precisions_decimal_digits": list(PRECISIONS),
        "signed_bridge_identity": "Dcorr=Ppi-I69pi-Qpaper=Qexact-Qpaper",
        "rows": rows,
        "aggregate": aggregate,
        "runtime": {
            "priority": priority,
            "flint_threads": 1,
            "worker_count": 1,
            "elapsed_seconds": time.monotonic() - started,
            "resumable_cache": relative(CACHE),
            "cache_fingerprint": fingerprint,
            "cpu_park_threshold_percent": CPU_PARK_THRESHOLD,
        },
        "dependencies": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "block21_points": {"path": relative(BLOCK21_POINTS), "sha256": file_hash(BLOCK21_POINTS)},
            "later_points": {"path": relative(LATER_POINTS), "sha256": file_hash(LATER_POINTS)},
            "selector_builder": {"path": relative(Path(cells.__file__).resolve()), "sha256": file_hash(Path(cells.__file__).resolve())},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": "Scale the correlated cell correction to all 1,040 later recursive calls, then add complete recurrence-tail and source-arithmetic transport.",
        "proof_boundary": "Two analytic corrected-model coefficient cells at one saved height; no complete later-block recurrence, outer Hardy remainder, height uniformity, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
    }
    require(passed, "later coefficient-cell signed bridge pilot did not close")
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    tmp = RESULT.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp.replace(RESULT)
    note_tmp = NOTE.with_suffix(".md.tmp")
    note_tmp.write_text(render_note(artifact), encoding="utf-8")
    note_tmp.replace(NOTE)
    print(
        "built later coefficient-cell signed bridge pilot: "
        f"{aggregate['closed_cell_count']}/{aggregate['cell_count']} cells, "
        f"max inflation {aggregate['maximum_cell_to_point_inflation_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
