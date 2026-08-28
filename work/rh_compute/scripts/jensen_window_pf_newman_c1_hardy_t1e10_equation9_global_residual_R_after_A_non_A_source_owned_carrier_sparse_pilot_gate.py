#!/usr/bin/env python3
"""Run a sparse altered-representation pilot for the joined non-A source block."""

from __future__ import annotations

from bisect import bisect_left
from fractions import Fraction
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time
from typing import Any


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

REPO_ROOT = Path(__file__).resolve().parents[3]
REPO_VENDOR = REPO_ROOT / "work/rh_compute/vendor"
import sys

sys.path.insert(0, str(REPO_VENDOR))
from mpmath import mp


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_source_owned_carrier_sparse_pilot_gate"
)
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
ROWS = REPO_ROOT / f"work/rh_compute/results/{STEM}_rows.jsonl"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

REASSEMBLY_STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_folded_abel_source_roster_reassembly_gate"
)
REFERENCE_STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_reference_evaluator_gate"
)
DEPENDENCIES = {
    "source_roster_reassembly": REPO_ROOT / f"work/rh_compute/results/{REASSEMBLY_STEM}.json",
    "source_reference_evaluator": REPO_ROOT / f"work/rh_compute/results/{REFERENCE_STEM}.json",
}
REFERENCE_MODULE = REPO_ROOT / f"work/rh_compute/scripts/{REFERENCE_STEM}.py"

A = 159_577
B = 5_122_421
L = 2_481_422
HEIGHT = 10_000_000_000
FIRST = 622
ORDINARY_LAST = 39_852
A_WINDOW_FIRST = 39_853
A_WINDOW_LAST = 39_936
WORKERS = 1
MP_DPS = 80
PHASE_BLOCK = 127
TAIL_SWITCH = 16
TAIL_MAX_TERMS = 72
TAIL_RELATIVE_GOAL = "1e-62"
ALTERED_RELATIVE_ALLOWANCE = mp.mpf("1e-42")
SOURCE_MASS_ALLOWANCE = mp.mpf("2e-10")


def x_b(mode: int) -> Fraction:
    return Fraction(2 * mode, B)


def x_a(mode: int) -> Fraction:
    return Fraction(2 * mode, A)


PILOT_SPECS = (
    ("first_B_entry", x_b(622), True, False),
    ("dyadic_1_over_4096", Fraction(1, 4096), False, True),
    ("decimal_1_over_1000", Fraction(1, 1000), False, False),
    ("first_A_exit", x_a(622), True, False),
    ("closest_pair_A_630", x_a(630), True, False),
    ("closest_pair_B_20223", x_b(20_223), True, False),
    ("last_B_entry", x_b(39_936), True, False),
    ("dyadic_1_over_64", Fraction(1, 64), False, True),
    ("dyadic_1_over_32", Fraction(1, 32), False, True),
    ("dyadic_1_over_8", Fraction(1, 8), False, True),
    ("quarter", Fraction(1, 4), True, True),
    ("third", Fraction(1, 3), False, True),
    ("two_fifths", Fraction(2, 5), False, True),
    ("seven_sixteenths", Fraction(7, 16), False, True),
    ("last_interior_A_exit", x_a(39_852), False, False),
    ("first_A_window_exit", x_a(39_853), True, False),
    ("last_lower_A_exit", x_a(39_894), False, False),
    ("midpoint", Fraction(1, 2), True, True),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


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


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency {relative(path)}")
    return json.loads(path.read_text(encoding="utf-8"))


def load_reference_module():
    spec = importlib.util.spec_from_file_location("non_a_reference_evaluator", REFERENCE_MODULE)
    require(spec is not None and spec.loader is not None, "cannot load source reference evaluator")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def mp_fraction(value: Fraction) -> Any:
    return mp.mpf(value.numerator) / value.denominator


def mp_cis_pi_fraction(value: Fraction) -> Any:
    reduced = value % 2
    return mp.exp(mp.j * mp.pi * reduced.numerator / reduced.denominator)


def complex_record(value: Any, digits: int = 55) -> dict[str, str]:
    return {
        "real": mp.nstr(mp.re(value), digits),
        "imag": mp.nstr(mp.im(value), digits),
    }


def real_text(value: Any, digits: int = 45) -> str:
    return mp.nstr(value, digits)


def fresnel_tail_positive(q: Any) -> tuple[Any, Any, int]:
    """Return the Abel tail integral from positive q to infinity."""

    require(q > 0, "Fresnel tail requires q>0")
    term = mp.j / (mp.pi * q)
    series = term
    relative_goal = mp.mpf(TAIL_RELATIVE_GOAL)
    used = 1
    for index in range(TAIL_MAX_TERMS - 1):
        next_term = term * (2 * index + 1) / (mp.j * mp.pi * q * q)
        if abs(next_term) > abs(term):
            break
        series += next_term
        term = next_term
        used += 1
        if abs(term) <= relative_goal * max(mp.mpf(1), abs(series)):
            break
    tail = mp.exp(mp.j * mp.pi * q * q / 2) * series
    return tail, abs(term), used


def central_fresnel(q: Any) -> Any:
    return mp.fresnelc(q) + mp.j * mp.fresnels(q)


def lower_fresnel(q: Any, diagnostics: dict[str, Any] | None = None) -> Any:
    """Return the Abel integral from minus infinity to real q."""

    half_total = (1 + mp.j) / 2
    if abs(q) <= TAIL_SWITCH:
        return central_fresnel(q) + half_total
    tail, last_term, used = fresnel_tail_positive(abs(q))
    if diagnostics is not None:
        diagnostics["tail_calls"] += 1
        diagnostics["max_tail_terms"] = max(diagnostics["max_tail_terms"], used)
        diagnostics["max_last_tail_term"] = max(diagnostics["max_last_tail_term"], last_term)
    return (1 + mp.j) - tail if q > 0 else tail


def upper_fresnel(q: Any, diagnostics: dict[str, Any] | None = None) -> Any:
    """Return the Abel integral from real q to plus infinity."""

    if abs(q) <= TAIL_SWITCH:
        return (1 + mp.j) / 2 - central_fresnel(q)
    tail, last_term, used = fresnel_tail_positive(abs(q))
    if diagnostics is not None:
        diagnostics["tail_calls"] += 1
        diagnostics["max_tail_terms"] = max(diagnostics["max_tail_terms"], used)
        diagnostics["max_last_tail_term"] = max(diagnostics["max_last_tail_term"], last_term)
    return tail if q > 0 else (1 + mp.j) - tail


def phase_value(mode: int, x: Fraction) -> Any:
    return mp_cis_pi_fraction(Fraction(mode, 1) - Fraction(mode * mode, 1) / x)


def weighted_phase_sum(x: Fraction, first: int, last: int, block_size: int = PHASE_BLOCK) -> Any:
    """Return sum m*(-1)^m*exp(-i*pi*m^2/x) with rational block reseeding."""

    pieces: list[Any] = []
    second_ratio = mp_cis_pi_fraction(Fraction(-2, 1) / x)
    for start in range(first, last + 1, block_size):
        stop = min(last + 1, start + block_size)
        phase = phase_value(start, x)
        ratio = mp_cis_pi_fraction(Fraction(1, 1) - Fraction(2 * start + 1, 1) / x)
        local: list[Any] = []
        for mode in range(start, stop):
            local.append(mode * phase)
            phase *= ratio
            ratio *= second_ratio
        pieces.append(mp.fsum(local))
    return mp.fsum(pieces)


def source_value(reference: Any, x: Fraction) -> tuple[Any, dict[str, Any]]:
    observed = reference.recurrence_theta_current(x, Fraction(0), L + 1, block_size=256)[2]
    value = mp.mpc(repr(observed.real), repr(observed.imag))
    record: dict[str, Any] = {
        "method": "O(L)_rational_phase_block_recurrence",
        "block_size": 256,
        "value": complex_record(value),
    }
    return value, record


def periodic_source_check(reference: Any, x: Fraction, observed: Any) -> dict[str, Any]:
    period = reference.minimal_period(x, Fraction(0))
    expected = reference.periodic_high_precision(x, Fraction(0), L + 1, period)[2]
    mass = mp.mpf((L + 1) * (A + B) // 2)
    error = abs(observed - expected)
    return {
        "minimal_period": period,
        "periodic_high_precision": complex_record(expected),
        "absolute_error": real_text(error),
        "source_mass_normalized_error": real_text(error / mass),
        "passed": error / mass <= SOURCE_MASS_ALLOWANCE,
    }


def selected_carriers(x: Fraction, diagnostics: dict[str, Any]) -> dict[str, Any]:
    x_mp = mp_fraction(x)
    root = mp.sqrt(2 / x_mp)
    total = 1 + mp.j
    scale = root / x_mp
    endpoint_a = mp_cis_pi_fraction(x * A * A / 4)

    ordinary_weighted = weighted_phase_sum(x, FIRST, ORDINARY_LAST)
    ordinary_p = scale * total * ordinary_weighted

    high_p_terms: list[Any] = []
    high_a_pair_terms: list[Any] = []
    for mode in range(A_WINDOW_FIRST, A_WINDOW_LAST + 1):
        phase = phase_value(mode, x)
        coefficient = phase * mode * scale
        q_minus = mp.sqrt(x_mp / 2) * (A - 2 * mode / x_mp)
        q_plus = mp.sqrt(x_mp / 2) * (A + 2 * mode / x_mp)
        lower_minus = lower_fresnel(q_minus, diagnostics)
        lower_plus = lower_fresnel(q_plus, diagnostics)
        high_p_terms.append(coefficient * total)
        high_a_pair_terms.append(
            -2 * endpoint_a / (mp.j * mp.pi * x_mp)
            + coefficient * (lower_plus - lower_minus)
        )

    high_p = mp.fsum(high_p_terms)
    high_a_pair = mp.fsum(high_a_pair_terms)
    full_p = ordinary_p + high_p
    carrier = -full_p - high_a_pair
    return {
        "ordinary_P_sum": ordinary_p,
        "A_window_P_sum": high_p,
        "A_window_A_pair_sum": high_a_pair,
        "full_P_sum": full_p,
        "joined_carrier": carrier,
    }


def raw_carriers(x: Fraction, diagnostics: dict[str, Any]) -> Any:
    """Evaluate -I plus the two original endpoint rows directly."""

    x_mp = mp_fraction(x)
    root = mp.sqrt(2 / x_mp)
    scale = root / x_mp
    endpoint_scale = 1 / (mp.j * mp.pi * x_mp)
    endpoint_a = mp_cis_pi_fraction(x * A * A / 4)
    endpoint_b = mp_cis_pi_fraction(x * B * B / 4)
    endpoint_difference = (endpoint_b - endpoint_a) * endpoint_scale
    terms: list[Any] = []

    for mode in range(FIRST, A_WINDOW_LAST + 1):
        phase = phase_value(mode, x)
        coefficient = phase * mode * scale
        q_a = mp.sqrt(x_mp / 2) * (A - 2 * mode / x_mp)
        q_b = mp.sqrt(x_mp / 2) * (B - 2 * mode / x_mp)
        lower_a = lower_fresnel(q_a, diagnostics)
        lower_b = lower_fresnel(q_b, diagnostics)
        upper_b = upper_fresnel(q_b, diagnostics)
        finite_mode = endpoint_difference + coefficient * (lower_b - lower_a)

        if mode <= ORDINARY_LAST:
            endpoint_row = endpoint_difference - coefficient * (lower_a + upper_b)
        else:
            q_a_negative_mode = mp.sqrt(x_mp / 2) * (A + 2 * mode / x_mp)
            a_negative = (
                -endpoint_a * endpoint_scale
                + coefficient * lower_fresnel(q_a_negative_mode, diagnostics)
            )
            b_positive = endpoint_b * endpoint_scale - coefficient * upper_b
            endpoint_row = b_positive - a_negative
        terms.append(-finite_mode + endpoint_row)

    return mp.fsum(terms)


def physical_density(x: Fraction, value: Any) -> Any:
    x_mp = mp_fraction(x)
    normalization = 2 * (mp.pi / (32 * HEIGHT)) ** mp.mpf("0.25")
    weight = mp.exp(mp.j * HEIGHT / 2 * mp.log((1 - x_mp) / x_mp))
    weight /= (x_mp * (1 - x_mp)) ** mp.mpf("0.25")
    return normalization * mp.re(mp.exp(-mp.j * mp.pi / 8) * weight * value)


def transition_roster() -> list[tuple[Fraction, str, int]]:
    events = [(x_b(mode), "B_entry", mode) for mode in range(FIRST, A_WINDOW_LAST + 1)]
    events.extend((x_a(mode), "A_exit", mode) for mode in range(FIRST, 39_894 + 1))
    events.sort(key=lambda item: item[0])
    require(len(events) == 78_588, "transition roster count drift")
    require(all(events[index][0] < events[index + 1][0] for index in range(len(events) - 1)), "event collision")
    return events


def transition_record(x: Fraction, events: list[tuple[Fraction, str, int]]) -> dict[str, Any]:
    coordinates = [item[0] for item in events]
    index = bisect_left(coordinates, x)
    candidates = []
    if index < len(events):
        candidates.append(events[index])
    if index > 0:
        candidates.append(events[index - 1])
    nearest = min(candidates, key=lambda item: abs(item[0] - x))
    distance = abs(nearest[0] - x)
    b_last = min(A_WINDOW_LAST, math.floor(B * x / 2))
    a_last = min(39_894, math.floor(A * x / 2))
    return {
        "events_strictly_below": index,
        "B_entries_reached": max(0, b_last - FIRST + 1),
        "A_exits_reached": max(0, a_last - FIRST + 1),
        "nearest_event": {
            "kind": nearest[1],
            "mode": nearest[2],
            "x": str(nearest[0]),
            "distance": str(distance),
            "distance_decimal": format(float(distance), ".17g"),
        },
    }


def config_signature() -> str:
    payload = {
        "points": [(label, str(x), raw, periodic) for label, x, raw, periodic in PILOT_SPECS],
        "mp_dps": MP_DPS,
        "phase_block": PHASE_BLOCK,
        "tail_switch": TAIL_SWITCH,
        "tail_max_terms": TAIL_MAX_TERMS,
        "dependencies": {name: file_hash(path) for name, path in DEPENDENCIES.items()},
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()


def load_cache(signature: str) -> dict[str, dict[str, Any]]:
    cached: dict[str, dict[str, Any]] = {}
    if not ROWS.is_file():
        return cached
    for line in ROWS.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        if row.get("config_signature") == signature and row.get("passed") is True:
            cached[row["label"]] = row
    return cached


def append_cache(row: dict[str, Any]) -> None:
    ROWS.parent.mkdir(parents=True, exist_ok=True)
    with ROWS.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(row, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def evaluate_row(
    reference: Any,
    events: list[tuple[Fraction, str, int]],
    label: str,
    x: Fraction,
    raw_check: bool,
    periodic_check: bool,
    signature: str,
) -> dict[str, Any]:
    started = time.time()
    diagnostics = {"tail_calls": 0, "max_tail_terms": 0, "max_last_tail_term": mp.mpf(0)}
    source, source_record = source_value(reference, x)
    selected = selected_carriers(x, diagnostics)
    assembly = source + selected["joined_carrier"]
    component_mass = abs(source) + abs(selected["joined_carrier"])
    cancellation_ratio = abs(assembly) / max(mp.mpf(1), component_mass)

    raw_record = None
    if raw_check:
        raw = raw_carriers(x, diagnostics)
        discrepancy = abs(raw - selected["joined_carrier"])
        denominator = max(mp.mpf(1), abs(raw), abs(selected["joined_carrier"]))
        relative = discrepancy / denominator
        raw_record = {
            "joined_carrier": complex_record(raw),
            "absolute_discrepancy": real_text(discrepancy),
            "relative_discrepancy": real_text(relative),
            "passed": relative <= ALTERED_RELATIVE_ALLOWANCE,
        }

    periodic_record = periodic_source_check(reference, x, source) if periodic_check else None
    passed = (raw_record is None or raw_record["passed"]) and (
        periodic_record is None or periodic_record["passed"]
    )
    row = {
        "config_signature": signature,
        "label": label,
        "x": str(x),
        "x_decimal": format(float(x), ".17g"),
        "source": source_record,
        "selected_representation": {
            "ordinary_P_sum": complex_record(selected["ordinary_P_sum"]),
            "A_window_P_sum": complex_record(selected["A_window_P_sum"]),
            "A_window_A_pair_sum": complex_record(selected["A_window_A_pair_sum"]),
            "full_P_sum": complex_record(selected["full_P_sum"]),
            "joined_carrier": complex_record(selected["joined_carrier"]),
            "source_plus_joined_carrier": complex_record(assembly),
            "component_mass": real_text(component_mass),
            "cancellation_ratio": real_text(cancellation_ratio),
            "normalized_physical_density": real_text(physical_density(x, assembly)),
        },
        "altered_endpoint_Fresnel_cross_check": raw_record,
        "periodic_source_cross_check": periodic_record,
        "transition_position": transition_record(x, events),
        "tail_diagnostics": {
            "tail_calls": diagnostics["tail_calls"],
            "maximum_terms_used": diagnostics["max_tail_terms"],
            "maximum_last_retained_term": real_text(diagnostics["max_last_tail_term"]),
            "rigorous_remainder_bound": False,
        },
        "runtime_seconds": round(time.time() - started, 3),
        "passed": passed,
    }
    require(passed, f"pilot row failed at {label}")
    return row


def render_note(artifact: dict[str, Any]) -> str:
    summary = artifact["pilot_summary"]
    lines = [
        "# Sparse source-minus-owned-carriers altered-representation pilot",
        "",
        "Date: 2026-08-26",
        "",
        "Status: **diagnostic sparse pilot passed; interval physical quadrature remains open**.",
        "",
        "## Common mathematical level",
        "",
        "This gate evaluates only the pointwise source/ownership block from Core 11.476:",
        "",
        "```text",
        "S(x)=sum_(n=0)^L f_x(n)",
        "     -sum_(m=622)^39936 P_m(x)",
        "     -sum_(m=39853)^39936[A_m(x)+A_-m(x)].",
        "```",
        "",
        "The saved B-trace, B-outer, and translated A-face artifacts are physical integrals or bounds.",
        "They are deliberately not inserted as pointwise constants.  Their ownership is unchanged.",
        "",
        "The primary evaluator joins `P_m+A_m+A_-m` on the 84 A-window modes before",
        "summation.  Selected rows are independently replayed as `-I_m` plus the original",
        "ordinary and A-window endpoint rows.",
        "",
        "## Results",
        "",
        f"- Rational sample rows: `{summary['row_count']}`.",
        f"- Full 39,315-mode altered-form rows: `{summary['altered_cross_check_count']}`.",
        f"- Maximum altered-form relative discrepancy: `{summary['maximum_altered_relative_discrepancy']}`.",
        f"- Periodic full-source oracle rows: `{summary['periodic_source_cross_check_count']}`.",
        f"- Maximum source-mass-normalized recurrence error: `{summary['maximum_source_mass_normalized_error']}`.",
        f"- Cancellation-ratio range: `{summary['minimum_cancellation_ratio']}` to `{summary['maximum_cancellation_ratio']}`.",
        f"- Maximum sampled normalized physical-density magnitude: `{summary['maximum_abs_normalized_physical_density']}`.",
        "",
        "| label | x | nearest transition | cancellation ratio | physical density |",
        "|---|---:|---|---:|---:|",
    ]
    for row in artifact["rows"]:
        nearest = row["transition_position"]["nearest_event"]
        lines.append(
            f"| {row['label']} | {row['x']} | {nearest['kind']} {nearest['mode']} "
            f"(distance {nearest['distance']}) | "
            f"{row['selected_representation']['cancellation_ratio']} | "
            f"{row['selected_representation']['normalized_physical_density']} |"
        )
    lines.extend(
        [
            "",
            "## Decision",
            "",
            "The joined pointwise representation is numerically coherent on the sparse roster and is the",
            "selected evaluator for the next bounded integration experiment.  The raw endpoint-plus-Fresnel",
            "form remains a checker only.  This gate does not turn the sparse samples or the asymptotic",
            "Fresnel-tail diagnostics into an interval error theorem.",
            "",
            "The next obligation is a physical-transform ownership ledger for the complete source, finite",
            "`P/A` carriers, and unchanged physical channels.  Reuse certified transformed integrals where",
            "they exist; launch event-aware panels only for a term that has no such transform.  Only a joined",
            "interval result below the inherited non-A allowance may be promoted.",
            "",
            "## Pi provenance",
            "",
            "Every `pi` comes from the original quadratic Kummer character, its integer Fourier transform,",
            "the canonical Fresnel primitive, or the inherited equation-(9) physical normalization.",
            "No fitted or geometric occurrence is introduced.",
            "",
            "## Boundary",
            "",
            artifact["proof_boundary"],
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    mp.dps = MP_DPS
    require(CHECKER.is_file(), "missing independent checker")
    loaded = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in loaded.values()), "a dependency gate failed")
    require(
        loaded["source_roster_reassembly"]["decision"][
            "source_roster_minus_owned_bulk_and_A_face_route_selected"
        ] is True,
        "source-minus-carriers route dependency drift",
    )
    require(
        loaded["source_reference_evaluator"]["decision"][
            "reference_evaluator_ready_as_fast_evaluator_oracle"
        ] is True,
        "source reference evaluator dependency drift",
    )
    require(
        file_hash(REFERENCE_MODULE)
        == loaded["source_reference_evaluator"]["sources"]["builder"]["sha256"],
        "source reference implementation hash drift",
    )

    specs = sorted(PILOT_SPECS, key=lambda item: item[1])
    require(len({item[0] for item in specs}) == len(specs), "duplicate pilot label")
    require(len({item[1] for item in specs}) == len(specs), "duplicate pilot coordinate")
    require(all(Fraction(0) < item[1] <= Fraction(1, 2) for item in specs), "pilot point outside half-domain")

    reference = load_reference_module()
    events = transition_roster()
    signature = config_signature()
    cached = load_cache(signature)
    rows: list[dict[str, Any]] = []
    for label, x, raw_check, periodic_check in specs:
        if label in cached:
            row = cached[label]
        else:
            row = evaluate_row(reference, events, label, x, raw_check, periodic_check, signature)
            append_cache(row)
        rows.append(row)
        print(f"completed {label} x={x}", flush=True)

    require(all(row["passed"] is True for row in rows), "a sparse pilot row failed")
    altered_errors = [
        mp.mpf(row["altered_endpoint_Fresnel_cross_check"]["relative_discrepancy"])
        for row in rows
        if row["altered_endpoint_Fresnel_cross_check"] is not None
    ]
    source_errors = [
        mp.mpf(row["periodic_source_cross_check"]["source_mass_normalized_error"])
        for row in rows
        if row["periodic_source_cross_check"] is not None
    ]
    cancellation = [mp.mpf(row["selected_representation"]["cancellation_ratio"]) for row in rows]
    densities = [abs(mp.mpf(row["selected_representation"]["normalized_physical_density"])) for row in rows]
    require(max(altered_errors) <= ALTERED_RELATIVE_ALLOWANCE, "altered representation mismatch")
    require(max(source_errors) <= SOURCE_MASS_ALLOWANCE, "periodic source oracle mismatch")

    proof_boundary = (
        "A deterministic floating sparse-point pilot at t=10^10 for the complete source minus its owned "
        "positive full-line P block and 84 paired A atoms, including rational-phase block reseeding, exact "
        "transition labels, selected full-roster periodic source checks, and altered endpoint-plus-Fresnel "
        "replays. The Fresnel asymptotic stopping diagnostic is not a rigorous remainder enclosure, the O(L) "
        "source recurrence has no uniform floating error theorem, and no interval x quadrature is performed. "
        "The B-trace, B-outer, and translated A-face physical channels are not reevaluated or mixed into this "
        "pointwise gate. No non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, "
        "PF-infinity, RH, or prize-level conclusion is proved."
    )
    artifact = {
        "kind": STEM,
        "status": "sparse_joined_source_owned_carrier_pilot_passed_physical_transform_ledger_selected",
        "passed": True,
        "scope": {
            "height": HEIGHT,
            "x_interval": ["0", "1/2"],
            "source_odd_roster": [A, B],
            "source_count": L + 1,
            "positive_P_mode_range": [FIRST, A_WINDOW_LAST],
            "paired_A_mode_range": [A_WINDOW_FIRST, A_WINDOW_LAST],
        },
        "mathematical_level_audit": {
            "evaluated_here": "pointwise_x_integrand_source_and_owned_P/A_carriers",
            "left_unchanged": [
                "physical_B_trace_window",
                "physical_completed_B_outer_channel",
                "physical_translated_A_face_channel",
            ],
            "forbidden_mixing_avoided": True,
            "reason": "the saved unchanged channels are physical integrals or bounds, not pointwise x values",
        },
        "evaluator": {
            "source": "validated O(L) exact-rational phase-reseeded recurrence on the L+1 endpoint-complete roster",
            "ordinary_modes": "joined contribution -P_m",
            "A_window_modes": "joined contribution -P_m-A_m-A_-m before summation",
            "altered_cross_check": "-I_m plus A_m+B_m or B_m-A_-m from canonical Fresnel tails",
            "phase_block": PHASE_BLOCK,
            "mp_decimal_digits": MP_DPS,
            "workers": WORKERS,
        },
        "transition_schedule": {
            "event_count": len(events),
            "first": {"x": str(events[0][0]), "kind": events[0][1], "mode": events[0][2]},
            "last": {"x": str(events[-1][0]), "kind": events[-1][1], "mode": events[-1][2]},
            "minimum_gap": "882/817420575917",
            "all_events_distinct": True,
        },
        "pilot_summary": {
            "row_count": len(rows),
            "altered_cross_check_count": len(altered_errors),
            "periodic_source_cross_check_count": len(source_errors),
            "maximum_altered_relative_discrepancy": real_text(max(altered_errors)),
            "maximum_source_mass_normalized_error": real_text(max(source_errors)),
            "minimum_cancellation_ratio": real_text(min(cancellation)),
            "maximum_cancellation_ratio": real_text(max(cancellation)),
            "maximum_abs_normalized_physical_density": real_text(max(densities)),
            "all_rows_passed": True,
        },
        "rows": rows,
        "decision": {
            "common_pointwise_mathematical_level_respected": True,
            "joined_source_P_A_evaluator_sparse_pilot_passed": True,
            "endpoint_plus_39315_Fresnel_form_sparse_cross_check_passed": True,
            "raw_endpoint_Fresnel_form_selected_for_production": False,
            "joined_source_owned_carrier_form_selected_for_next_pilot": True,
            "naive_value_space_panel_route_selected": False,
            "physical_transform_ownership_ledger_selected_next": True,
            "event_aware_panel_quadrature_retained_as_fallback": True,
            "sparse_values_promoted_to_interval_bound": False,
            "uniform_source_floating_error_theorem_proved": False,
            "Fresnel_tail_interval_remainder_proved": False,
            "physical_x_quadrature_completed": False,
            "non_A_bound_proved": False,
            "rh_implication": False,
        },
        "next_obligation": (
            "Build an exact physical-transform ownership ledger for the complete source, finite P block, "
            "paired A block, B trace, completed B outer channel, and translated A-face channel. Identify which "
            "terms already have certified physical values or bounds and isolate only genuinely unevaluated "
            "transforms. Use event-aware x panels only as a fallback for those terms, with a cumulative interval "
            "error budget below the inherited non-A allowance."
        ),
        "proof_boundary": proof_boundary,
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "rows": {"path": relative(ROWS), "sha256": file_hash(ROWS)},
        },
        "runtime": {
            "workers": WORKERS,
            "threads": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.time() - started, 3),
        },
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("completed sparse joined source-owned-carrier altered-representation pilot", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
