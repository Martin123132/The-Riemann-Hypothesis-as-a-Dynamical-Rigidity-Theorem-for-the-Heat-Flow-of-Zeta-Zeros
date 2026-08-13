#!/usr/bin/env python3
"""Test the signed source-to-corrected-model bridge on block-21 witnesses."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
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
import jensen_window_pf_newman_c1_hardy_block20_exact_eq69_saddle_family_gate as eq69
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
import jensen_window_pf_newman_c1_hardy_t1e10_crossblock_recursive_output_budget_gate as transport
import jensen_window_pf_newman_c1_hardy_t1e10_crossblock_structural_atlas_gate as atlas
from flint import acb, arb, ctx


TELEMETRY = atlas.TELEMETRY
ENDPOINTS = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_retained_endpoint_gate.json"
WEIGHTS = transport.WEIGHTS
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_pilot_gate.json"
CACHE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_pilot_gate_cache.jsonl"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_pilot_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_pilot_gate.py"
WITNESS_CHAINS = (432, 443, 448, 489)
PRECISIONS = (70, 110)
OUTPUT_INDEX = 15
TAIL_TARGET = Fraction(1, 10**30)
FORMULA_GAP_REPORTING_CEILING = Fraction(5, 10**14)
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


def exact_arb(payload: str) -> arb:
    return cells.fraction_ball(cells.binary128_fraction(payload))


def exact_acb(payload: list[str]) -> acb:
    return point_balls.binary128_complex(payload)


def upper(value: arb) -> Fraction:
    return cells.bound_fraction(value.upper())


def upper_abs(value: acb) -> Fraction:
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
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def load_witnesses() -> dict[int, dict[str, Any]]:
    chains: dict[int, dict[str, Any]] = {}
    selected = set(WITNESS_CHAINS)
    for line in TELEMETRY.read_text(encoding="utf-8").splitlines():
        record = json.loads(line)
        chain = int(record.get("chain", 0))
        if chain not in selected:
            continue
        kind = record["type"]
        if kind == "chain":
            chains[chain] = {"header": record, "levels": {}, "q_terms": {}, "steps": {}}
        elif kind == "level":
            chains[chain]["levels"][int(record["level"])] = record
        elif kind == "q_terms":
            chains[chain]["q_terms"][int(record["nit"])] = record
        elif kind == "recurrence":
            chains[chain]["steps"][int(record["nit"])] = record
    require(set(chains) == selected, "block-21 witness roster missing")
    for chain, data in chains.items():
        require(int(data["header"]["block"]) == 21, f"chain {chain} block drift")
        require(int(data["header"]["mit"]) == 2, f"chain {chain} route drift")
        require(set(data["levels"]) >= {1, 2}, f"chain {chain} level roster drift")
        require(1 in data["q_terms"] and 1 in data["steps"], f"chain {chain} recurrence payload missing")
    return chains


def fraction_record(value: Fraction) -> dict[str, Any]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": decimal(value),
    }


def sector_contract(data: dict[str, Any], dps: int) -> dict[str, Any]:
    ctx.dps = dps
    parent = data["levels"][1]
    child = data["levels"][2]
    require(int(parent["degree"]) == 3, "block-21 sector parent degree drift")
    phi1, phi2, phi3 = (cells.binary128_fraction(value) for value in parent["coefficients_hex"])
    require(phi3 != 0, "block-21 zero cubic coefficient")
    length = int(parent["length"])
    xi = phi1 + 2 * phi2 * length + 3 * phi3 * length**2
    child_length = int(child["length"])
    require(xi.numerator // xi.denominator == child_length, "block-21 xi selector drift")
    delta = child_length + 1 - xi
    second_minimum = min(2 * phi2, 2 * phi2 + 6 * phi3 * length)
    c_gap = 1 + phi1
    require(0 < delta < 1 and second_minimum > 0 and c_gap > Fraction(1, 2), "block-21 sector decay drift")
    slanted = arb(3).sqrt() * cells.fraction_ball(second_minimum) / 4
    if phi3 < 0:
        values = {
            "b": ("5*pi/6", delta / 2, slanted, -phi3),
            "c": ("pi/2", c_gap, arb(0), -phi3),
        }
    else:
        values = {
            "b": ("pi/2", delta, arb(0), phi3),
            "c": ("pi/6", c_gap / 2, slanted, phi3),
        }
    result: dict[str, Any] = {
        "phi3_sign": "negative" if phi3 < 0 else "positive",
        "length": length,
        "child_length": child_length,
    }
    for family, (angle, linear, quadratic, cubic) in values.items():
        require(linear > 0 and cubic > 0, f"block-21 {family} decay drift")
        result[f"{family}_contract"] = {
            "angle": angle,
            "linear_decay": fraction_record(linear),
            "quadratic_decay": point_balls.arb_record(quadratic, 45),
            "cubic_decay": fraction_record(cubic),
        }
    return result


def saddle_integral_sum(data: dict[str, Any], dps: int, mathematical_pi: bool) -> dict[str, Any]:
    ctx.dps = dps
    ctx.threads = 1
    parent = data["levels"][1]
    child = data["levels"][2]
    require(int(parent["degree"]) == 3, "block-21 equation-(69) parent degree drift")
    a1, a2, a3 = (exact_arb(value) for value in parent["coefficients_hex"])
    a1_fraction, a2_fraction, a3_fraction = (
        cells.binary128_fraction(value) for value in parent["coefficients_hex"]
    )
    length = int(parent["length"])
    require(a2_fraction > 0, "block-21 equation-(69) quadratic drift")
    require(min(2 * a2_fraction, 2 * a2_fraction + 6 * a3_fraction * length) > 0, "block-21 curvature drift")
    child_length = int(child["length"])
    jbot = eq69.ceil_fraction(a1_fraction)
    require(jbot in (0, 1), "block-21 equation-(69) starting index drift")
    require(bool(child["subtract_one"]) == (jbot == 1), "block-21 child-start adapter drift")
    source_two_pi = -exact_arb(data["header"]["tpm_hex"])
    phase = 2 * arb.pi() if mathematical_pi else source_two_pi
    tolerance = arb(f"1e-{dps - 20}")
    total = acb(0)
    minimum_accuracy = 10**9
    for index in range(jbot, child_length + 1):
        n = arb(index)

        def integrand(y: acb, _analytic: bool) -> acb:
            g = n * y - a1 * y - a2 * y**2 - a3 * y**3
            return (I * phase * g).exp()

        value = acb.integral(
            integrand,
            arb(0),
            arb(length),
            rel_tol=tolerance,
            abs_tol=tolerance,
            eval_limit=200000,
            depth_limit=30,
        )
        require(value.is_finite(), f"block-21 equation-(69) nonfinite integral n={index}")
        minimum_accuracy = min(minimum_accuracy, int(value.rel_accuracy_bits()))
        total += value
    return {
        "indices": list(range(jbot, child_length + 1)),
        "sum": total,
        "minimum_relative_accuracy_bits": minimum_accuracy,
        "source_two_pi": source_two_pi,
    }


def paper_components(data: dict[str, Any], dps: int) -> dict[str, Any]:
    ctx.dps = dps
    parent = data["levels"][1]
    child = data["levels"][2]
    coefficients = [cells.binary128_fraction(value) for value in parent["coefficients_hex"]]
    a1_fraction, a2_fraction, a3_fraction = coefficients
    require(a2_fraction > 0, "block-21 W2--W5 quadratic drift")
    a1, a2, a3 = (cells.fraction_ball(value) for value in coefficients)
    length = int(parent["length"])
    xi_fraction = a1_fraction + 2 * a2_fraction * length + 3 * a3_fraction * length**2
    child_length = int(child["length"])
    require(xi_fraction.numerator // xi_fraction.denominator == child_length, "block-21 W2--W5 xi drift")
    frac_fraction = xi_fraction - child_length
    require(0 < frac_fraction < 1, "block-21 W2--W5 boundary drift")
    xi = cells.fraction_ball(xi_fraction)
    delta = cells.fraction_ball(1 - frac_fraction)
    xr = 2 * a2
    sqrt_xr = xr.sqrt()
    pi = arb.pi()
    sqrt_two = arb(2).sqrt()
    epi_minus = acb(1 / sqrt_two, -1 / sqrt_two)
    phase_at_endpoint = a1 * length + a2 * length**2 + a3 * length**3
    endpoint = acb(0, -2 * pi * phase_at_endpoint).exp()
    internal_endpoint = endpoint.conjugate()

    def fresnel_boundary(argument: arb) -> acb:
        z = epi_minus * (pi / xr).sqrt() * argument
        return I * epi_minus * acb(0, -pi * argument**2 / xr).exp() * (1 - z.erf()) / (2 * sqrt_xr)

    w2 = (
        I
        * internal_endpoint
        * (
            epi_minus
            * acb(0, -pi * delta**2 / xr).exp()
            * (1 - (epi_minus * (pi / xr).sqrt() * delta).erf())
            / (2 * sqrt_xr)
            - 1 / (2 * pi * delta)
        )
    ).conjugate()
    if a1_fraction < 0:
        selector = "W3"
        argument = a1 + 1
        w34 = (fresnel_boundary(argument) - I / (2 * pi * argument)).conjugate()
    else:
        require(a1_fraction > 0, "block-21 W3/W4 selector on zero")
        selector = "W4"
        argument = a1
        w34 = (fresnel_boundary(argument) - I * internal_endpoint / (2 * pi * xi)).conjugate()
    w5 = (
        I
        / (2 * pi)
        * (
            internal_endpoint * ((xi + 1).digamma() - delta.digamma())
            - ((a1 + 1).digamma() - (child_length + 1 - a1).digamma())
        )
    ).conjugate()
    half_sum = ((1 + internal_endpoint) / 2).conjugate()
    return {
        "selector": selector,
        "endpoint": endpoint,
        "w2": w2,
        "w34": w34,
        "w5": w5,
        "half_sum": half_sum,
        "aggregate": w2 + w34 + w5 + half_sum,
    }


def select_cutoff(data: dict[str, Any], sector: dict[str, Any]) -> tuple[Fraction, Fraction]:
    cutoff = log_full.INITIAL_CUTOFF
    length = int(data["levels"][1]["length"])
    while True:
        worst = max(
            log_full.tail_bound(data, sector, endpoint, family, cutoff)
            for family in ("b", "c")
            for endpoint in (0, length)
        )
        if worst < TAIL_TARGET:
            return cutoff, worst
        cutoff *= 2
        require(cutoff <= log_full.MAXIMUM_CUTOFF, "block-21 ray cutoff exceeded safety ceiling")


def finite_zero_mode(data: dict[str, Any], dps: int) -> acb:
    ctx.dps = dps
    parent = data["levels"][1]
    phi1, phi2, phi3 = (exact_arb(value) for value in parent["coefficients_hex"])
    length = int(parent["length"])
    pi = arb.pi()

    def integrand(y: acb, _analytic: bool) -> acb:
        phase = phi1 * y + phi2 * y**2 + phi3 * y**3
        return (-2 * pi * I * phase).exp()

    tolerance = arb(f"1e-{dps - 25}")
    total = acb(0)
    lower_break = 0
    while lower_break < length:
        upper_break = min(length, lower_break + 13)
        total += acb.integral(
            integrand,
            acb(lower_break),
            acb(upper_break),
            abs_tol=tolerance,
            rel_tol=tolerance,
            deg_limit=64,
            eval_limit=200000,
            depth_limit=64,
        )
        lower_break = upper_break
    return total


def exact_contour_nonsaddle(
    data: dict[str, Any], sector: dict[str, Any], dps: int, cutoff: Fraction
) -> dict[str, Any]:
    ctx.dps = dps
    parent = data["levels"][1]
    length = int(parent["length"])
    child_length = int(data["levels"][2]["length"])
    breaks = log_full.breaks_for_cutoff(cutoff)
    b_zero = log_pilot.ray_integral(data, sector, 0, "b", dps, cutoff, breaks)
    b_length = log_pilot.ray_integral(data, sector, length, "b", dps, cutoff, breaks)
    c_zero = log_pilot.ray_integral(data, sector, 0, "c", dps, cutoff, breaks)
    c_length = log_pilot.ray_integral(data, sector, length, "c", dps, cutoff, breaks)
    b_source = b_zero["total"] - b_length["total"]
    c_source = (c_length["total"] - c_zero["total"]).conjugate()
    phi1, phi2, phi3 = (exact_arb(value) for value in parent["coefficients_hex"])
    pi = arb.pi()
    endpoint_phase = phi1 * length + phi2 * length**2 + phi3 * length**3
    endpoint = (-2 * pi * I * endpoint_phase).exp()
    endpoint_half = (1 + endpoint) / 2
    a_endpoint = I * (endpoint - 1) * cells.fraction_ball(log_pilot.harmonic(child_length)) / (2 * pi)
    zero_mode = finite_zero_mode(data, dps) if cells.binary128_fraction(parent["coefficients_hex"][0]) > 0 else acb(0)
    return {
        "exact_nonsaddle": endpoint_half + a_endpoint + b_source + c_source + zero_mode,
        "maximum_ray_origin_upper": max(ray["origin_upper"] for ray in (b_zero, b_length, c_zero, c_length)),
        "maximum_ray_tail_upper": max(ray["tail_upper"] for ray in (b_zero, b_length, c_zero, c_length)),
    }


def evaluate(data: dict[str, Any], dps: int, cutoff: Fraction) -> dict[str, Any]:
    ctx.dps = dps
    ctx.threads = 1
    sector = sector_contract(data, dps)
    source_integral = saddle_integral_sum(data, dps, mathematical_pi=False)
    paper_integral = saddle_integral_sum(data, dps, mathematical_pi=True)
    require(source_integral["indices"] == paper_integral["indices"], "block-21 integral roster drift")
    paper = paper_components(data, dps)
    require(paper["selector"] in ("W3", "W4"), "block-21 selector drift")
    contour = exact_contour_nonsaddle(data, sector, dps, cutoff)

    parent = data["levels"][1]
    child = data["levels"][2]
    q_terms = data["q_terms"][1]
    step = data["steps"][1]
    tpm = exact_arb(data["header"]["tpm_hex"])
    source_two_pi = -tpm
    parent_source = point_balls.direct_sum(parent, tpm)
    parent_paper = point_balls.direct_sum(parent, -2 * arb.pi())
    child_raw = point_balls.direct_sum(child, tpm)
    child_adapted = point_balls.adapt(child, child_raw)
    saddle_model = exact_acb(step["multiplier_hex"]) * child_adapted
    source_q = exact_acb(q_terms["qq_hex"])
    source_t5 = exact_acb(q_terms["t5_hex"])
    endpoint = exact_acb(q_terms["endpoint_hex"])
    q_extra = I * endpoint / (source_two_pi * (int(child["length"]) + 1))
    source_nonsaddle = source_q - source_t5 - q_extra

    d_source = parent_source - saddle_model - source_q
    exact_w1_replacement = source_integral["sum"] - saddle_model - source_t5
    d_after_w1 = d_source - exact_w1_replacement
    d_after_t2 = d_after_w1 + q_extra
    expected_after_w1 = parent_source - source_integral["sum"] - (source_q - source_t5)
    expected_after_t2 = parent_source - source_integral["sum"] - source_nonsaddle
    parent_transport = parent_paper - parent_source
    integral_transport = paper_integral["sum"] - source_integral["sum"]
    nonsaddle_transport = paper["aggregate"] - source_nonsaddle
    d_corrected_transport = d_after_t2 + parent_transport - integral_transport - nonsaddle_transport
    d_corrected_direct = parent_paper - paper_integral["sum"] - paper["aggregate"]
    independent_target = parent_paper - paper_integral["sum"]
    formula_target_gap = contour["exact_nonsaddle"] - independent_target
    exact_minus_paper = contour["exact_nonsaddle"] - paper["aggregate"]

    return {
        "selector": paper["selector"],
        "phi3_sign": sector["phi3_sign"],
        "source_integral": source_integral["sum"],
        "paper_integral": paper_integral["sum"],
        "paper_nonsaddle": paper["aggregate"],
        "exact_nonsaddle": contour["exact_nonsaddle"],
        "independent_target": independent_target,
        "formula_target_gap": formula_target_gap,
        "d_source": d_source,
        "d_after_w1": d_after_w1,
        "d_after_t2": d_after_t2,
        "d_corrected_transport": d_corrected_transport,
        "d_corrected_direct": d_corrected_direct,
        "exact_minus_paper": exact_minus_paper,
        "w1_identity_gap": d_after_w1 - expected_after_w1,
        "t2_identity_gap": d_after_t2 - expected_after_t2,
        "corrected_identity_gap": d_corrected_transport - d_corrected_direct,
        "endpoint_closure_gap": d_corrected_direct - exact_minus_paper,
        "q_extra": q_extra,
        "parent_transport": parent_transport,
        "integral_transport": integral_transport,
        "nonsaddle_transport": nonsaddle_transport,
        "maximum_ray_origin_upper": contour["maximum_ray_origin_upper"],
        "maximum_ray_tail_upper": contour["maximum_ray_tail_upper"],
        "minimum_integral_accuracy_bits": min(
            source_integral["minimum_relative_accuracy_bits"],
            paper_integral["minimum_relative_accuracy_bits"],
        ),
    }


def cache_fingerprint() -> str:
    payload = {
        "builder": file_hash(BUILDER),
        "checker": file_hash(CHECKER),
        "telemetry": file_hash(TELEMETRY),
        "endpoints": file_hash(ENDPOINTS),
        "weights": file_hash(WEIGHTS),
        "dependencies": {
            "log_pilot": file_hash(Path(log_pilot.__file__).resolve()),
            "log_full": file_hash(Path(log_full.__file__).resolve()),
            "eq69": file_hash(Path(eq69.__file__).resolve()),
            "point_balls": file_hash(Path(point_balls.__file__).resolve()),
        },
        "witnesses": list(WITNESS_CHAINS),
        "precisions": list(PRECISIONS),
        "tail_target": [TAIL_TARGET.numerator, TAIL_TARGET.denominator],
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def load_cache(fingerprint: str) -> dict[int, dict[str, Any]]:
    if not CACHE.exists():
        CACHE.write_text(json.dumps({"kind": "block21_bridge_pilot_cache", "fingerprint": fingerprint}, sort_keys=True) + "\n", encoding="utf-8")
        return {}
    lines = [line for line in CACHE.read_text(encoding="utf-8").splitlines() if line.strip()]
    require(lines, "block-21 bridge cache empty")
    header = json.loads(lines[0])
    if header.get("fingerprint") != fingerprint:
        stale = CACHE.with_name(f"{CACHE.stem}.stale.{header.get('fingerprint', 'unknown')[:12]}{CACHE.suffix}")
        require(not stale.exists(), f"stale cache destination already exists: {stale}")
        CACHE.replace(stale)
        CACHE.write_text(json.dumps({"kind": "block21_bridge_pilot_cache", "fingerprint": fingerprint}, sort_keys=True) + "\n", encoding="utf-8")
        return {}
    rows: dict[int, dict[str, Any]] = {}
    for line in lines[1:]:
        row = json.loads(line)
        chain = int(row["chain"])
        require(chain not in rows, f"duplicate cached bridge chain {chain}")
        rows[chain] = row
    return rows


def append_cache(row: dict[str, Any]) -> None:
    with CACHE.open("a", encoding="utf-8", buffering=1) as handle:
        handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def endpoint_rows() -> dict[int, dict[str, Any]]:
    artifact = json.loads(ENDPOINTS.read_text(encoding="utf-8"))
    rows = {int(row["chain"]): row for row in artifact["rows"] if int(row["chain"]) in WITNESS_CHAINS}
    require(set(rows) == set(WITNESS_CHAINS), "block-21 endpoint witness rows missing")
    return rows


def build_row(
    chain: int,
    data: dict[str, Any],
    endpoint_row: dict[str, Any],
    weight_factor: Fraction,
) -> dict[str, Any]:
    sector = sector_contract(data, PRECISIONS[1])
    cutoff, selected_tail = select_cutoff(data, sector)
    low = evaluate(data, PRECISIONS[0], cutoff)
    high = evaluate(data, PRECISIONS[1], cutoff)
    complex_fields = (
        "source_integral",
        "paper_integral",
        "paper_nonsaddle",
        "exact_nonsaddle",
        "independent_target",
        "formula_target_gap",
        "d_source",
        "d_after_w1",
        "d_after_t2",
        "d_corrected_transport",
        "d_corrected_direct",
        "exact_minus_paper",
        "w1_identity_gap",
        "t2_identity_gap",
        "corrected_identity_gap",
        "endpoint_closure_gap",
        "q_extra",
        "parent_transport",
        "integral_transport",
        "nonsaddle_transport",
    )
    for name in complex_fields:
        require(low[name].overlaps(high[name]), f"chain {chain} precision nonoverlap: {name}")
    require(low["selector"] == high["selector"] == endpoint_row["selector"], f"chain {chain} selector mismatch")
    require(low["phi3_sign"] == high["phi3_sign"], f"chain {chain} Phi3 sign mismatch")
    zero = acb(0)
    for name in ("formula_target_gap", "w1_identity_gap", "t2_identity_gap", "corrected_identity_gap", "endpoint_closure_gap"):
        require(high[name].overlaps(zero), f"chain {chain} bridge closure failed: {name}")
    retained = Fraction(endpoint_row["complete_majorant_upper"])
    correction_upper = upper_abs(high["exact_minus_paper"])
    return {
        "chain": chain,
        "block": int(endpoint_row["block"]),
        "sum_index": int(endpoint_row["sum_index"]),
        "branch": int(endpoint_row["branch"]),
        "parent_length": int(endpoint_row["parent_length"]),
        "child_length": int(endpoint_row["child_length"]),
        "selector": high["selector"],
        "phi3_sign": high["phi3_sign"],
        "output_index": OUTPUT_INDEX,
        "outer_weight_factor": decimal(weight_factor),
        "retained_local_majorant_upper": decimal(retained),
        "retained_transported_contribution_upper": decimal(weight_factor * retained),
        "exact_correction_magnitude_upper": decimal(correction_upper),
        "exact_to_retained_ratio_upper": decimal(correction_upper / retained),
        "exact_correction_within_retained_majorant": correction_upper <= retained,
        "transported_exact_correction_upper": decimal(weight_factor * correction_upper),
        "selected_ray_cutoff": decimal(cutoff),
        "selected_tail_bound_upper": decimal(selected_tail),
        "maximum_ray_origin_upper": decimal(high["maximum_ray_origin_upper"]),
        "maximum_ray_tail_upper": decimal(high["maximum_ray_tail_upper"]),
        "minimum_integral_accuracy_bits": int(high["minimum_integral_accuracy_bits"]),
        "formula_target_gap_upper": decimal(upper_abs(high["formula_target_gap"])),
        "w1_identity_gap_upper": decimal(upper_abs(high["w1_identity_gap"])),
        "t2_identity_gap_upper": decimal(upper_abs(high["t2_identity_gap"])),
        "corrected_identity_gap_upper": decimal(upper_abs(high["corrected_identity_gap"])),
        "endpoint_closure_gap_upper": decimal(upper_abs(high["endpoint_closure_gap"])),
        "d_source": point_balls.acb_record(high["d_source"], 45),
        "d_after_w1": point_balls.acb_record(high["d_after_w1"], 45),
        "d_after_t2": point_balls.acb_record(high["d_after_t2"], 45),
        "d_corrected": point_balls.acb_record(high["d_corrected_direct"], 45),
        "exact_minus_paper": point_balls.acb_record(high["exact_minus_paper"], 45),
        "precision_overlap": True,
        "signed_bridge_closed": True,
    }


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    roster = ", ".join(str(chain) for chain in WITNESS_CHAINS)
    return f"""# Hardy t=1e10 block-21 signed bridge pilot

Date: 2026-08-09

Status: representative finite signed-bridge gate closed; not a proof of a coefficient neighborhood, arbitrary height, or RH

This pilot tests chains `{roster}`.  They are selected before exact contour
evaluation to cover the largest output-15 block-21 contribution and every
occupied `(W3/W4, sign(Phi3))` branch.  All have parent length 119, so none is
a relabelled block-20 length-104 call.

For each witness, exact source algebra is transported in four signed stages:

```text
Dsrc = Psrc - MC - qsrc
D1   = Dsrc - (I69src - MC - t5src)
D2   = D1 + qextra = Psrc - I69src - Qsrc,drop
Dcorr = D2 + (Ppi-Psrc) - (I69pi-I69src) - (Qpaper-Qsrc,drop)
      = Ppi - I69pi - Qpaper.
```

An independent sign-aware endpoint-contour calculation supplies
`Qexact = Ppi-I69pi`.  Thus `Dcorr=Qexact-Qpaper` closes without a fitted
correction or cancellation between calls.

Witnesses closed: `{aggregate['closed_witness_count']}/{aggregate['witness_count']}`

Maximum formula-to-independent-target gap: `{aggregate['maximum_formula_target_gap_upper']}`

Maximum signed bridge identity gap: `{aggregate['maximum_signed_bridge_gap_upper']}`

Exact corrections inside the earlier retained majorant:
`{aggregate['exact_corrections_within_retained_count']}/{aggregate['witness_count']}`

The result is a representative finite bridge, not a proof for all 212
block-21 pivots, blocks 22--28, a coefficient neighborhood, arbitrary height,
or RH.  Its next valid use is to promote the length-generic machinery to the
complete block-21 roster with an append-only cache.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime-limit-seconds", type=float, default=10800.0)
    args = parser.parse_args()
    for path in (TELEMETRY, ENDPOINTS, WEIGHTS, CHECKER):
        require(path.is_file(), f"missing block-21 bridge dependency: {path}")
    priority = set_low_priority()
    ctx.threads = 1
    started = time.monotonic()
    witnesses = load_witnesses()
    endpoints = endpoint_rows()
    factors = transport.load_weight_factors()
    fingerprint = cache_fingerprint()
    cached = load_cache(fingerprint)
    require(set(cached).issubset(WITNESS_CHAINS), "unexpected block-21 bridge cache row")
    for chain in WITNESS_CHAINS:
        if chain in cached:
            continue
        elapsed = time.monotonic() - started
        if elapsed >= args.runtime_limit_seconds:
            print(f"parked block-21 bridge pilot after {len(cached)}/{len(WITNESS_CHAINS)} witnesses")
            return 75
        row = endpoints[chain]
        weight = factors[(21, int(row["sum_index"]), OUTPUT_INDEX, int(row["branch"]))]
        built = build_row(chain, witnesses[chain], row, weight)
        append_cache(built)
        cached[chain] = built
        print(f"closed block-21 bridge witness {chain} ({len(cached)}/{len(WITNESS_CHAINS)})", flush=True)

    rows = [cached[chain] for chain in WITNESS_CHAINS]
    require({(row["selector"], row["phi3_sign"]) for row in rows} == {
        ("W3", "negative"), ("W3", "positive"), ("W4", "negative"), ("W4", "positive")
    }, "block-21 witness branch coverage drift")
    maximum_identity = max(
        Fraction(row[name])
        for row in rows
        for name in ("w1_identity_gap_upper", "t2_identity_gap_upper", "corrected_identity_gap_upper", "endpoint_closure_gap_upper")
    )
    aggregate = {
        "witness_count": len(rows),
        "closed_witness_count": sum(bool(row["signed_bridge_closed"]) for row in rows),
        "branch_coverage": sorted({f"{row['selector']}:{row['phi3_sign']}" for row in rows}),
        "maximum_formula_target_gap_upper": decimal(max(Fraction(row["formula_target_gap_upper"]) for row in rows)),
        "maximum_signed_bridge_gap_upper": decimal(maximum_identity),
        "exact_corrections_within_retained_count": sum(bool(row["exact_correction_within_retained_majorant"]) for row in rows),
        "maximum_exact_to_retained_ratio_upper": decimal(max(Fraction(row["exact_to_retained_ratio_upper"]) for row in rows)),
        "maximum_transported_exact_correction_upper": decimal(max(Fraction(row["transported_exact_correction_upper"]) for row in rows)),
        "minimum_integral_accuracy_bits": min(int(row["minimum_integral_accuracy_bits"]) for row in rows),
    }
    passed = (
        aggregate["closed_witness_count"] == len(rows)
        and Fraction(aggregate["maximum_formula_target_gap_upper"]) < FORMULA_GAP_REPORTING_CEILING
        and aggregate["exact_corrections_within_retained_count"] == len(rows)
    )
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_pilot_gate",
        "status": "representative_block21_signed_source_bridge_closed" if passed else "representative_block21_signed_source_bridge_falsified",
        "scope": "Four exact block-21 witnesses covering the maximum output-15 contribution and every W3/W4 by Phi3-sign branch.",
        "method": {
            "selection": "Maximum retained transported output-15 contribution in each occupied (W3/W4, sign(Phi3)) class.",
            "bridge": "Exact source W1 replacement, source-only t2 deletion, pi transport, and total source-to-paper W2--W5 transport.",
            "independent_endpoint": "Sign-aware logarithmic endpoint rays with rigorous origin and far-tail boxes plus an adaptive finite zero mode over parent length 119.",
            "closure_criterion": "Every identity-gap interval overlaps zero; 5e-14 is only a regression ceiling for the deliberate origin enclosure, matching the admitted block-20 contour scale.",
            "no_special_function_split_needed": True,
            "reason": "The old special-function probe only partitions the total nonsaddle transport; the signed total is reconstructed directly from emitted q terms and published W2--W5.",
        },
        "precisions_decimal_digits": list(PRECISIONS),
        "runtime": {
            "priority": priority,
            "flint_threads": 1,
            "worker_count": 1,
            "elapsed_seconds": time.monotonic() - started,
            "resumable_cache": relative(CACHE),
            "cache_fingerprint": fingerprint,
        },
        "provenance": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "endpoints": {"path": relative(ENDPOINTS), "sha256": file_hash(ENDPOINTS)},
            "weights": {"path": relative(WEIGHTS), "sha256": file_hash(WEIGHTS)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "rows": rows,
        "aggregate": aggregate,
        "passed": passed,
        "boundary": {
            "proved": "The full signed exact-point bridge and independent corrected-model endpoint identity close on the four preselected block-21 witnesses.",
            "not_proved": "No claim is made for the remaining 208 block-21 calls, blocks 22--28, coefficient cells, arbitrary height, or RH.",
            "next_action": "Promote the length-generic bridge to all 212 block-21 calls using the append-only witness cache pattern, then transport exact corrections through all 15 outer outputs.",
        },
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print(
        "built block-21 signed bridge pilot: "
        f"passed={passed}, max_formula_gap={aggregate['maximum_formula_target_gap_upper']}, "
        f"retained={aggregate['exact_corrections_within_retained_count']}/{len(rows)}"
    )
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
