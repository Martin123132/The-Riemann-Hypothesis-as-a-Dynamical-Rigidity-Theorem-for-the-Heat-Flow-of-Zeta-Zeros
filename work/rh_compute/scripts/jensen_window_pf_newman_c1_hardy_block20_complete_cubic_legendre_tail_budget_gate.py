#!/usr/bin/env python3
"""Bound the complete saved-point cubic Legendre truncation at block 20."""

from __future__ import annotations

from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import jensen_window_pf_newman_c1_hardy_block20_coefficient_neighborhood_transport_gate as coefficient_transport
import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_block20_transported_endpoint_budget_gate as transport
import jensen_window_pf_newman_c1_hardy_cubic_child_legendre_adapter_gate as child_adapter
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import acb, arb, ctx


COEFFICIENT_TRANSPORT = coefficient_transport.RESULT
WEIGHT_FIXTURE = transport.WEIGHT_FIXTURE
WEIGHTS = transport.WEIGHTS
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_complete_cubic_legendre_tail_budget_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_complete_cubic_legendre_tail_budget_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_block20_complete_cubic_legendre_tail_budget_gate.py"
PRECISIONS = (70, 110)
WEIGHT_PRECISION = 140
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


def upper_decimal(value: Fraction) -> str:
    with localcontext() as context:
        context.prec = 100
        context.rounding = ROUND_CEILING
        return format(Decimal(value.numerator) / Decimal(value.denominator), ".80E")


def lower_decimal(value: Fraction) -> str:
    with localcontext() as context:
        context.prec = 100
        context.rounding = ROUND_FLOOR
        return format(Decimal(value.numerator) / Decimal(value.denominator), ".80E")


def upper(value: arb) -> Fraction:
    return cells.bound_fraction(value.upper())


def lower(value: arb) -> Fraction:
    return cells.bound_fraction(value.lower())


def upper_abs(value: arb | acb) -> Fraction:
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


def source_phase(coefficients: list[arb], index: int) -> arb:
    value = arb(0)
    for coefficient in reversed(coefficients):
        value = (value + coefficient) * index
    return value


def evaluate_chain(data: dict[str, Any], dps: int) -> dict[str, Any]:
    """Enclose all omitted Legendre orders while preserving source normalization."""
    ctx.dps = dps
    ctx.threads = 1
    header = data["header"]
    parent_level = data["levels"][1]
    child_level = data["levels"][2]
    step = data["steps"][1]

    parent = [cells.binary128_fraction(value) for value in parent_level["coefficients_hex"]]
    require(len(parent) == 3, "cubic parent required")
    a1_f, a2_f, a3_f = parent
    y_f = 2 * a2_f
    require(y_f > 0, "nonpositive parent quadratic coefficient")
    a1, y, a3 = (cells.fraction_ball(value) for value in (a1_f, y_f, a3_f))
    tpm = point_balls.binary128_ball(header["tpm_hex"])
    stored = [point_balls.binary128_ball(value) for value in child_level["coefficients_hex"]]
    multiplier = point_balls.binary128_complex(step["multiplier_hex"])
    orientation = 1 if bool(child_level["conjugate"]) else -1

    delta_sum = acb(0)
    direct_sum_upper = Fraction(0)
    elementary_sum_upper = Fraction(0)
    maximum_tail = Fraction(0)
    maximum_z = Fraction(0)
    minimum_discriminant: Fraction | None = None
    minimum_exact_denominator: Fraction | None = None
    minimum_source_denominator: Fraction | None = None
    term_rows: list[dict[str, Any]] = []

    for index in range(int(child_level["length"]) + 1):
        u = cells.fraction_ball(Fraction(index) - a1_f)
        z = 12 * a3 * u / y**2
        discriminant = 1 + z
        require(lower(discriminant) > 0, "stationary discriminant is not positive")
        root = discriminant.sqrt()
        stationary_x = 2 * u / (y * (1 + root))
        exact_h = u * stationary_x - (y / 2) * stationary_x**2 - a3 * stationary_x**3
        truncated_h = u**2 / (2 * y) - a3 * u**3 / y**3
        tail = exact_h - truncated_h

        exact_denominator = root.sqrt()
        source_denominator_f, centered_denominator_f = child_adapter.denominator_forms(parent, index)
        require(source_denominator_f == centered_denominator_f, "source denominator identity drift")
        source_denominator = cells.fraction_ball(source_denominator_f)
        taylor_denominator = 1 + z / 4 - 3 * z**2 / 32
        require(source_denominator.overlaps(taylor_denominator), "fourth-root Taylor denominator drift")
        require(lower(source_denominator) > 0, "source denominator is not positive")

        angle = source_phase(stored, index)
        exact_angle = angle + orientation * tail
        exact_term = acb(arb(0), tpm * exact_angle).exp() / exact_denominator
        source_term = acb(arb(0), tpm * angle).exp() / source_denominator
        term_delta = exact_term - source_term
        direct_upper = upper_abs(term_delta)
        amplitude_upper = upper_abs(1 / exact_denominator - 1 / source_denominator)
        phase_upper = upper_abs(tpm) * upper_abs(tail) / lower(source_denominator)
        elementary_upper = amplitude_upper + phase_upper
        require(direct_upper <= elementary_upper, "direct term enclosure escaped elementary majorant")

        delta_sum += term_delta
        direct_sum_upper += direct_upper
        elementary_sum_upper += elementary_upper
        maximum_tail = max(maximum_tail, upper_abs(tail))
        maximum_z = max(maximum_z, upper_abs(z))
        discriminant_low = lower(discriminant)
        exact_denominator_low = lower(exact_denominator)
        source_denominator_low = lower(source_denominator)
        minimum_discriminant = (
            discriminant_low if minimum_discriminant is None else min(minimum_discriminant, discriminant_low)
        )
        minimum_exact_denominator = (
            exact_denominator_low
            if minimum_exact_denominator is None
            else min(minimum_exact_denominator, exact_denominator_low)
        )
        minimum_source_denominator = (
            source_denominator_low
            if minimum_source_denominator is None
            else min(minimum_source_denominator, source_denominator_low)
        )
        term_rows.append(
            {
                "index": index,
                "z_abs_upper": upper_decimal(upper_abs(z)),
                "discriminant_lower": lower_decimal(discriminant_low),
                "exact_denominator_lower": lower_decimal(exact_denominator_low),
                "source_denominator_lower": lower_decimal(source_denominator_low),
                "legendre_phase_tail_abs_upper": upper_decimal(upper_abs(tail)),
                "direct_term_gap_abs_upper": upper_decimal(direct_upper),
                "elementary_term_majorant_upper": upper_decimal(elementary_upper),
            }
        )

    if bool(child_level["conjugate"]):
        delta_sum = delta_sum.conjugate()
    local_delta = multiplier * delta_sum
    multiplier_upper = upper_abs(multiplier)
    local_majorant = multiplier_upper * direct_sum_upper
    local_elementary_majorant = multiplier_upper * elementary_sum_upper
    require(upper_abs(local_delta) <= local_majorant, "local tail escaped termwise majorant")
    return {
        "actual_delta": local_delta,
        "multiplier_abs_upper": multiplier_upper,
        "direct_sum_upper": direct_sum_upper,
        "elementary_sum_upper": elementary_sum_upper,
        "local_majorant": local_majorant,
        "local_elementary_majorant": local_elementary_majorant,
        "maximum_tail": maximum_tail,
        "maximum_z": maximum_z,
        "minimum_discriminant": minimum_discriminant,
        "minimum_exact_denominator": minimum_exact_denominator,
        "minimum_source_denominator": minimum_source_denominator,
        "term_rows": term_rows,
    }


def build() -> dict[str, Any]:
    for path in (COEFFICIENT_TRANSPORT, WEIGHT_FIXTURE, WEIGHTS, CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    priority = set_low_priority()
    coefficient_result = json.loads(COEFFICIENT_TRANSPORT.read_text(encoding="utf-8"))
    endpoint_outputs = {
        int(row["output_index"]): row for row in coefficient_result["transported_outputs"]
    }
    require(set(endpoint_outputs) == set(range(1, 16)), "endpoint output roster drift")

    recursive = native_q.load_recursive_chains()
    require(len(recursive) == 374, "recursive roster drift")
    ctx.dps = WEIGHT_PRECISION
    ctx.threads = 1
    weights = transport.load_weights()
    rows: list[dict[str, Any]] = []
    local_deltas: dict[int, acb] = {}
    for chain in sorted(recursive):
        low = evaluate_chain(recursive[chain], PRECISIONS[0])
        high = evaluate_chain(recursive[chain], PRECISIONS[1])
        require(low["actual_delta"].overlaps(high["actual_delta"]), f"chain {chain} precision drift")
        require(low["local_majorant"] >= upper_abs(high["actual_delta"]), f"chain {chain} low-precision escape")
        header = recursive[chain]["header"]
        child = recursive[chain]["levels"][2]
        step = recursive[chain]["steps"][1]
        local_deltas[chain] = high["actual_delta"].conjugate() if bool(step["conjugate"]) else high["actual_delta"]
        rows.append(
            {
                "chain": chain,
                "sum_index": int(header["sum_index"]),
                "branch": int(header["branch"]),
                "child_length": int(child["length"]),
                "source_child_conjugated": bool(child["conjugate"]),
                "outer_recurrence_conjugated": bool(step["conjugate"]),
                "precision_overlap": True,
                "maximum_dimensionless_z_abs_upper": upper_decimal(high["maximum_z"]),
                "minimum_stationary_discriminant_lower": lower_decimal(high["minimum_discriminant"]),
                "minimum_exact_fourth_root_denominator_lower": lower_decimal(high["minimum_exact_denominator"]),
                "minimum_source_taylor_denominator_lower": lower_decimal(high["minimum_source_denominator"]),
                "maximum_exact_legendre_phase_tail_abs_upper": upper_decimal(high["maximum_tail"]),
                "multiplier_abs_upper": upper_decimal(high["multiplier_abs_upper"]),
                "child_termwise_direct_gap_upper": upper_decimal(high["direct_sum_upper"]),
                "child_termwise_elementary_majorant_upper": upper_decimal(high["elementary_sum_upper"]),
                "local_recurrence_tail_majorant_upper": upper_decimal(high["local_majorant"]),
                "local_recurrence_tail_elementary_majorant_upper": upper_decimal(high["local_elementary_majorant"]),
                "local_recurrence_tail_actual": point_balls.acb_record(high["actual_delta"], 48),
                "local_recurrence_tail_actual_abs_upper": upper_decimal(upper_abs(high["actual_delta"])),
                "terms": high["term_rows"],
            }
        )

    row_by_chain = {int(row["chain"]): row for row in rows}
    transported_outputs: list[dict[str, Any]] = []
    for output_index in range(1, 16):
        tail_majorant = Fraction(0)
        signed_tail = arb(0)
        for chain in sorted(recursive):
            row = row_by_chain[chain]
            weight = weights[(int(row["sum_index"]), output_index, int(row["branch"]))]
            tail_majorant += (
                weight["amplitude_fraction"]
                * weight["phase_abs_upper"]
                * Fraction(row["local_recurrence_tail_majorant_upper"])
            )
            signed_tail += transport.projected(weight, local_deltas[chain])
        endpoint = Fraction(endpoint_outputs[output_index]["cell_endpoint_majorant_upper"])
        combined = endpoint + tail_majorant
        transported_outputs.append(
            {
                "output_index": output_index,
                "output_label": transport.output_label(output_index),
                "coefficient_cell_endpoint_majorant_upper": upper_decimal(endpoint),
                "complete_cubic_legendre_tail_majorant_upper": upper_decimal(tail_majorant),
                "endpoint_plus_legendre_tail_majorant_upper": upper_decimal(combined),
                "combined_to_requested_scale_ratio": upper_decimal(combined / REQUESTED_ERROR_SCALE),
                "combined_within_requested_scale": combined < REQUESTED_ERROR_SCALE,
                "signed_legendre_tail_diagnostic": point_balls.arb_record(signed_tail, 48),
                "signed_legendre_tail_diagnostic_abs_upper": upper_decimal(upper_abs(signed_tail)),
            }
        )

    def maximum(field: str) -> tuple[Fraction, int]:
        return max((Fraction(row[field]), int(row["chain"])) for row in rows)

    def minimum(field: str) -> tuple[Fraction, int]:
        return min((Fraction(row[field]), int(row["chain"])) for row in rows)

    def output_maximum(field: str) -> tuple[Fraction, int]:
        return max((Fraction(row[field]), int(row["output_index"])) for row in transported_outputs)

    max_z = maximum("maximum_dimensionless_z_abs_upper")
    min_disc = minimum("minimum_stationary_discriminant_lower")
    max_phase_tail = maximum("maximum_exact_legendre_phase_tail_abs_upper")
    max_local = maximum("local_recurrence_tail_majorant_upper")
    max_output_tail = output_maximum("complete_cubic_legendre_tail_majorant_upper")
    max_combined = output_maximum("endpoint_plus_legendre_tail_majorant_upper")
    closes = all(bool(row["combined_within_requested_scale"]) for row in transported_outputs)
    aggregate = {
        "recursive_call_count": len(rows),
        "child_index_count": sum(len(row["terms"]) for row in rows),
        "child_length_histogram": {
            str(length): sum(int(row["child_length"]) == length for row in rows) for length in (1, 2)
        },
        "precision_ladder_decimal_digits": list(PRECISIONS),
        "precision_overlap_count": sum(bool(row["precision_overlap"]) for row in rows),
        "maximum_dimensionless_z_abs_upper": upper_decimal(max_z[0]),
        "maximum_dimensionless_z_abs_witness": max_z[1],
        "minimum_stationary_discriminant_lower": lower_decimal(min_disc[0]),
        "minimum_stationary_discriminant_witness": min_disc[1],
        "maximum_exact_legendre_phase_tail_abs_upper": upper_decimal(max_phase_tail[0]),
        "maximum_exact_legendre_phase_tail_witness": max_phase_tail[1],
        "maximum_local_recurrence_tail_majorant_upper": upper_decimal(max_local[0]),
        "maximum_local_recurrence_tail_majorant_witness": max_local[1],
        "maximum_transported_legendre_tail_majorant_upper": upper_decimal(max_output_tail[0]),
        "maximum_transported_legendre_tail_majorant_witness": max_output_tail[1],
        "maximum_endpoint_plus_legendre_tail_majorant_upper": upper_decimal(max_combined[0]),
        "maximum_endpoint_plus_legendre_tail_majorant_witness": max_combined[1],
        "requested_error_scale": decimal(REQUESTED_ERROR_SCALE),
        "outputs_with_combined_majorant_within_requested_scale": sum(
            bool(row["combined_within_requested_scale"]) for row in transported_outputs
        ),
    }
    status = (
        "complete_saved_point_cubic_legendre_tail_and_endpoint_budget_close"
        if closes
        else "complete_saved_point_cubic_legendre_tail_breaks_endpoint_budget"
    )
    return {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_complete_cubic_legendre_tail_budget_gate",
        "status": status,
        "scope": "The exact 374 saved recursive calls and their 753 child indices at block 20.",
        "theorem": {
            "stationary_coordinate": "z=12*a3*(k-a1)/y^2 and F''(X)=y*sqrt(1+z)",
            "exact_amplitude_denominator": "A(z)=(1+z)^(1/4)",
            "source_amplitude_denominator": "D2(z)=1+z/4-3*z^2/32",
            "exact_phase": "H(u)=u*X-(y/2)*X^2-a3*X^3, X=2*u/[y*(1+sqrt(1+z))]",
            "source_phase": "H3(u)=u^2/(2*y)-a3*u^3/y^3",
            "term_bound": "|exp(i*tpm*H)/A-exp(i*tpm*H3)/D2| <= |1/A-1/D2|+|tpm|*|H-H3|/D2",
            "normalization": "The exact tail is added after the fixed integer/parity/orientation normalization, so no exact-2*pi shortcut is used.",
            "local_transport": "Multiply the termwise child majorant by the exact saved recurrence-multiplier magnitude.",
            "output_transport": "Multiply each local majorant by the exact saved amplitude and the outward bound for the saved complex branch phase.",
        },
        "rows": rows,
        "transported_outputs": transported_outputs,
        "aggregate": aggregate,
        "dependencies": {
            "coefficient_cell_endpoint": {"path": relative(COEFFICIENT_TRANSPORT), "sha256": file_hash(COEFFICIENT_TRANSPORT)},
            "weight_fixture": {"path": relative(WEIGHT_FIXTURE), "sha256": file_hash(WEIGHT_FIXTURE)},
            "weight_rows": {"path": relative(WEIGHTS), "sha256": file_hash(WEIGHTS)},
        },
        "runtime": {"priority": priority, "flint_threads": 1},
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": (
            "Extend the joint phase/amplitude tail bound over accuracy-refined coefficient subcells, then certify "
            "the remaining source arithmetic, other blocks, cross-block accumulation, and outer Hardy remainder."
        ),
        "proof_boundary": (
            "This is a complete all-orders cubic Legendre truncation bound only at the finite saved child indices. "
            "It does not prove a coefficient-neighborhood tail theorem, source-roundoff closure, a height-uniform "
            "recurrence, other-block or outer-Hardy control, Lambda<=0, PF-infinity, RH, or a prize-level conclusion."
        ),
    }


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 complete cubic Legendre-tail budget

Date: 2026-08-09
Status: {artifact['status']}; finite saved-point theorem only, not a proof of the complete evaluator or RH

## Exact Joint Remainder

For `u=k-a1` and `z=12*a3*u/y^2`, the exact cubic stationary map gives

```text
F''(X)=y*sqrt(1+z),
A(z)=(1+z)^(1/4),
D2(z)=1+z/4-3*z^2/32.
```

Thus the omitted transformed-child term has both a phase part `H-H3` and an
amplitude part `1/A-1/D2`.  The gate evaluates the exact radicals rather than
stopping at the first omitted Taylor coefficient.  It preserves the source's
binary `tpm` and fixed integer/parity/orientation normalization, so no exact
`2*pi` periodicity is assumed.

## Certified Roster

```text
recursive calls                         {aggregate['recursive_call_count']}
saved child indices                     {aggregate['child_index_count']}
precision overlaps                      {aggregate['precision_overlap_count']} / {aggregate['recursive_call_count']}
maximum |z|                             {aggregate['maximum_dimensionless_z_abs_upper']}
minimum 1+z                             {aggregate['minimum_stationary_discriminant_lower']}
maximum |H-H3|                          {aggregate['maximum_exact_legendre_phase_tail_abs_upper']}
maximum local recurrence tail           {aggregate['maximum_local_recurrence_tail_majorant_upper']}
maximum transported tail                {aggregate['maximum_transported_legendre_tail_majorant_upper']}
maximum endpoint plus tail              {aggregate['maximum_endpoint_plus_legendre_tail_majorant_upper']}
outputs below 0.005                     {aggregate['outputs_with_combined_majorant_within_requested_scale']} / 15
```

The endpoint and Legendre-tail columns remain separately recorded.  The
combined number uses absolute call-by-call transport and does not use the
signed output diagnostic.

## Boundary

This closes the all-orders cubic Legendre phase-and-amplitude truncation only
at the exact saved block-20 child indices.  Tail transport over coefficient
subcells, source arithmetic, other blocks, cross-block accumulation, the outer
Hardy representation, `Lambda<=0`, PF-infinity, RH, and a prize-level theorem
remain open.
"""


def main() -> int:
    artifact = build()
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built complete cubic Legendre-tail budget: "
        f"{artifact['aggregate']['recursive_call_count']} calls, "
        f"{artifact['aggregate']['outputs_with_combined_majorant_within_requested_scale']}/15 outputs below 0.005"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
