#!/usr/bin/env python3
"""Stress-test later coefficient cells against the complete cubic recurrence tail."""

from __future__ import annotations

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

import jensen_window_pf_newman_c1_hardy_block20_coefficient_cell_complete_cubic_legendre_tail_budget_gate as cell_tail
import jensen_window_pf_newman_c1_hardy_block20_coefficient_neighborhood_transport_gate as coefficient_transport
import jensen_window_pf_newman_c1_hardy_block20_complete_cubic_legendre_tail_budget_gate as point_tail
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
import jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_full_gate as later_cells
import jensen_window_pf_newman_c1_hardy_t1e10_crossblock_recursive_output_budget_gate as transport
from flint import acb, arb, ctx


FULL_CELLS = later_cells.RESULT
TELEMETRY = later_cells.TELEMETRY
WEIGHTS = transport.WEIGHTS
SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/telemetry/zeta14cubicmult_telemetry.f90"
ROUNDING_MODE_PROBE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/binary128_rounding_mode_probe/probe_result.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_later_complete_recurrence_stress_cell_pilot_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_later_complete_recurrence_stress_cell_pilot_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_t1e10_later_complete_recurrence_stress_cell_pilot_gate.py"
WITNESS_CHAINS = (626, 653, 3406, 3423)
EXPECTED_BLOCKS = {626: 21, 653: 21, 3406: 28, 3423: 28}
PRECISIONS = (70, 110)
REQUESTED_SCALE = Fraction(5, 1000)
BINARY128_UNIT_ROUNDOFF = Fraction(1, 2**113)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


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
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def exact_fraction(record: dict[str, Any]) -> Fraction:
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def binary128_rounding_preimage(payload: str) -> dict[str, Any]:
    """Return the exact round-to-nearest Voronoi cell of a normal payload."""
    bits = int(payload, 16)
    exponent = (bits >> 112) & 0x7FFF
    require(0 < exponent < 0x7FFF, f"non-normal binary128 payload {payload}")
    center = cells.binary128_fraction(payload)
    if bits >> 127:
        predecessor = cells.binary128_fraction(f"{bits + 1:032X}")
        successor = cells.binary128_fraction(f"{bits - 1:032X}")
    else:
        predecessor = cells.binary128_fraction(f"{bits - 1:032X}")
        successor = cells.binary128_fraction(f"{bits + 1:032X}")
    require(predecessor < center < successor, "binary128 neighbor order drift")
    lower_bound = (predecessor + center) / 2
    upper_bound = (center + successor) / 2
    return {
        "hex": payload,
        "center": center,
        "lower": lower_bound,
        "upper": upper_bound,
        "lower_radius": center - lower_bound,
        "upper_radius": upper_bound - center,
    }


def rounding_record(payload: str) -> dict[str, Any]:
    record = binary128_rounding_preimage(payload)
    return {
        "hex": payload,
        "center": decimal(record["center"]),
        "lower": decimal(record["lower"]),
        "upper": decimal(record["upper"]),
        "lower_radius": decimal(record["lower_radius"]),
        "upper_radius": decimal(record["upper_radius"]),
    }


def fraction_interval_inside(lo: Fraction, hi: Fraction, outer_lo: Fraction, outer_hi: Fraction) -> bool:
    return outer_lo <= lo <= hi <= outer_hi


def rounding_preimage_inside_arb(payload: str, enclosure: arb) -> bool:
    record = binary128_rounding_preimage(payload)
    return lower(enclosure) <= record["lower"] and record["upper"] <= upper(enclosure)


def maximum_abs_interval(value: arb) -> Fraction:
    return max(abs(lower(value)), abs(upper(value)))


def complex_multiply_add_roundoff_upper(multiplier: acb, child: acb, source_q: acb) -> Fraction:
    """Uniform RN error for two products and two additions per component."""

    def component(a: Fraction, b: Fraction, q: Fraction) -> Fraction:
        products = a + b
        product_and_first_add = BINARY128_UNIT_ROUNDOFF * (2 + BINARY128_UNIT_ROUNDOFF) * products
        return product_and_first_add + BINARY128_UNIT_ROUNDOFF * (
            products + q + product_and_first_add
        )

    mr = maximum_abs_interval(multiplier.real)
    mi = maximum_abs_interval(multiplier.imag)
    cr = maximum_abs_interval(child.real)
    ci = maximum_abs_interval(child.imag)
    qr = maximum_abs_interval(source_q.real)
    qi = maximum_abs_interval(source_q.imag)
    real_error = component(mr * cr, mi * ci, qr)
    imag_error = component(mr * ci, mi * cr, qi)
    return real_error + imag_error


def exact_complex_multiply_add(data: dict[str, Any]) -> tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]:
    step = data["steps"][1]
    m_re, m_im = (cells.binary128_fraction(value) for value in step["multiplier_hex"])
    c_re, c_im = (cells.binary128_fraction(value) for value in step["state_before_hex"])
    q_re, q_im = (cells.binary128_fraction(value) for value in step["qq_hex"])
    exact = (m_re * c_re - m_im * c_im + q_re, m_re * c_im + m_im * c_re + q_im)
    emitted = tuple(cells.binary128_fraction(value) for value in step["state_pre_transform_hex"])
    return exact, emitted


def parent_rounding_guards(data: dict[str, Any], box: dict[str, Any]) -> dict[str, Any]:
    parent = data["levels"][1]
    a1_hex, phi2_hex, a3_hex = parent["coefficients_hex"]
    frac_hex = parent["frac_length_hex"]
    bounds = box["bounds"]
    radii = box["radii"]
    a1 = binary128_rounding_preimage(a1_hex)
    phi2 = binary128_rounding_preimage(phi2_hex)
    a3 = binary128_rounding_preimage(a3_hex)
    frac = binary128_rounding_preimage(frac_hex)
    guards = {
        "a1": fraction_interval_inside(a1["lower"], a1["upper"], *bounds["a1"]),
        "y_from_twice_phi2": fraction_interval_inside(2 * phi2["lower"], 2 * phi2["upper"], *bounds["y"]),
        "a3": fraction_interval_inside(a3["lower"], a3["upper"], *bounds["a3"]),
        "fracL": fraction_interval_inside(frac["lower"], frac["upper"], *bounds["fracL"]),
    }
    require(all(guards.values()), "parent binary128 rounding preimage escaped accepted cell")
    maximum_ratio = max(
        max(a1["lower_radius"], a1["upper_radius"]) / radii["a1"],
        2 * max(phi2["lower_radius"], phi2["upper_radius"]) / radii["y"],
        max(a3["lower_radius"], a3["upper_radius"]) / radii["a3"],
        max(frac["lower_radius"], frac["upper_radius"]) / radii["fracL"],
    )
    return {
        "guards": guards,
        "all_preimages_inside_cell": all(guards.values()),
        "maximum_rounding_radius_to_cell_radius_ratio": decimal(maximum_ratio),
        "a1": rounding_record(a1_hex),
        "phi2": rounding_record(phi2_hex),
        "a3": rounding_record(a3_hex),
        "fracL": rounding_record(frac_hex),
    }


def build() -> dict[str, Any]:
    for path in (FULL_CELLS, TELEMETRY, WEIGHTS, SOURCE, ROUNDING_MODE_PROBE, CHECKER):
        require(path.is_file(), f"missing stress-pilot dependency: {path}")
    priority = set_low_priority()
    ctx.threads = 1
    full = json.loads(FULL_CELLS.read_text(encoding="utf-8"))
    require(bool(full["passed"]), "later coefficient-cell dependency failed")
    rounding_mode = json.loads(ROUNDING_MODE_PROBE.read_text(encoding="utf-8"))
    require(bool(rounding_mode["passed"]), "binary128 rounding-mode dependency failed")
    require(rounding_mode["probe_values"]["rounding_mode"] == "nearest", "native rounding mode is not nearest")
    require(rounding_mode["probe_values"]["digits"] == "113", "native binary128 precision drift")
    full_rows = {int(row["chain"]): row for row in full["rows"]}
    require(set(WITNESS_CHAINS).issubset(full_rows), "stress-cell roster missing")
    chains = later_cells.load_chains(set(WITNESS_CHAINS))

    rows: list[dict[str, Any]] = []
    for chain in WITNESS_CHAINS:
        data = chains[chain]
        source_row = full_rows[chain]
        require(int(source_row["block"]) == EXPECTED_BLOCKS[chain], f"chain {chain} block drift")
        require(int(source_row["child_length"]) == 1, f"chain {chain} child-length drift")
        cell = source_row["refined_cell"]
        atlas_row = {"q_cell": cell}
        box = coefficient_transport.coefficient_box(data, atlas_row)

        low = cell_tail.evaluate_cell(data, atlas_row, PRECISIONS[0])
        high = cell_tail.evaluate_cell(data, atlas_row, PRECISIONS[1])
        require(len(low["term_gaps"]) == len(high["term_gaps"]) == 2, "stress tail term roster drift")
        require(all(a.overlaps(b) for a, b in zip(low["term_gaps"], high["term_gaps"])), f"chain {chain} tail precision nonoverlap")
        point = point_tail.evaluate_chain(data, PRECISIONS[1])
        require(point["local_majorant"] <= high["local_majorant"], f"chain {chain} point tail escaped cell tail")

        structural = coefficient_transport.structural_transport(data, atlas_row, box, PRECISIONS[1])
        parent_rounding = parent_rounding_guards(data, box)
        child_hex = data["levels"][2]["coefficients_hex"]
        child_rounding_guards = [
            rounding_preimage_inside_arb(payload, enclosure)
            for payload, enclosure in zip(child_hex, structural["child_coefficients"])
        ]
        require(all(child_rounding_guards), f"chain {chain} child rounding preimage escaped reconstructed cell")
        multiplier_hex = data["steps"][1]["multiplier_hex"]
        multiplier_rounding_guards = [
            rounding_preimage_inside_arb(multiplier_hex[0], structural["multiplier"].real),
            rounding_preimage_inside_arb(multiplier_hex[1], structural["multiplier"].imag),
        ]
        require(all(multiplier_rounding_guards), f"chain {chain} multiplier rounding preimage escaped anchored cell")
        state_hex = data["steps"][1]["state_before_hex"]
        state_rounding_guards = [
            rounding_preimage_inside_arb(state_hex[0], structural["child_kernel"].real),
            rounding_preimage_inside_arb(state_hex[1], structural["child_kernel"].imag),
        ]
        require(all(state_rounding_guards), f"chain {chain} child-state rounding preimage escaped cell kernel")

        source_q = point_balls.binary128_complex(data["steps"][1]["qq_hex"])
        operation_roundoff = complex_multiply_add_roundoff_upper(
            structural["multiplier"], structural["child_kernel"], source_q
        )
        exact_operation, emitted_operation = exact_complex_multiply_add(data)
        observed_operation_gap = abs(emitted_operation[0] - exact_operation[0]) + abs(
            emitted_operation[1] - exact_operation[1]
        )
        point_multiplier = point_balls.binary128_complex(multiplier_hex)
        point_child = point_balls.binary128_complex(state_hex)
        point_operation_roundoff = complex_multiply_add_roundoff_upper(point_multiplier, point_child, source_q)
        require(observed_operation_gap <= point_operation_roundoff, f"chain {chain} source multiply-add escaped RN envelope")

        tpm_rounding = binary128_rounding_preimage(data["header"]["tpm_hex"])
        ctx.dps = PRECISIONS[1]
        minus_two_pi = -2 * arb.pi()
        tpm_contains_mathematical_value = (
            tpm_rounding["lower"] <= lower(minus_two_pi)
            and upper(minus_two_pi) <= tpm_rounding["upper"]
        )
        require(tpm_contains_mathematical_value, f"chain {chain} tpm rounding cell missed -2*pi")

        rows.append(
            {
                "chain": chain,
                "block": int(source_row["block"]),
                "sum_index": int(source_row["sum_index"]),
                "branch": int(source_row["branch"]),
                "parent_length": int(source_row["parent_length"]),
                "child_length": int(source_row["child_length"]),
                "selector": source_row["selector"],
                "phi3_sign": source_row["phi3_sign"],
                "precision_overlap": True,
                "point_tail_contained_in_cell_tail": True,
                "maximum_dimensionless_z_abs_upper": point_tail.upper_decimal(high["maximum_z"]),
                "minimum_stationary_discriminant_lower": point_tail.lower_decimal(high["minimum_discriminant"]),
                "maximum_factored_legendre_phase_tail_upper": point_tail.upper_decimal(high["maximum_tail"]),
                "point_complete_tail_majorant_upper": point_tail.upper_decimal(point["local_majorant"]),
                "cell_complete_tail_majorant_upper": point_tail.upper_decimal(high["local_majorant"]),
                "cell_to_point_tail_inflation_upper": point_tail.upper_decimal(high["local_majorant"] / point["local_majorant"]),
                "term_gap_abs_uppers": [point_tail.upper_decimal(upper(value)) for value in high["term_gaps"]],
                "parent_binary128_rounding": parent_rounding,
                "child_binary128_rounding_preimages_inside_reconstructed_cell": child_rounding_guards,
                "multiplier_binary128_rounding_preimages_inside_anchored_cell": multiplier_rounding_guards,
                "child_state_binary128_rounding_preimages_inside_cell_kernel": state_rounding_guards,
                "tpm_rounding_preimage_contains_minus_two_pi": tpm_contains_mathematical_value,
                "source_complex_multiply_add_observed_gap": decimal(observed_operation_gap),
                "source_complex_multiply_add_point_roundoff_upper": decimal(point_operation_roundoff),
                "source_complex_multiply_add_cell_roundoff_upper": decimal(operation_roundoff),
            }
        )

    factors = transport.load_weight_factors()
    prior_outputs = {int(row["output_index"]): row for row in full["output_rows"]}
    output_rows: list[dict[str, Any]] = []
    for output_index in range(1, 16):
        tail = Fraction(0)
        roundoff = Fraction(0)
        for row in rows:
            key = (int(row["block"]), int(row["sum_index"]), output_index, int(row["branch"]))
            factor = factors[key]
            tail += factor * Fraction(row["cell_complete_tail_majorant_upper"])
            roundoff += factor * Fraction(row["source_complex_multiply_add_cell_roundoff_upper"])
        prior = Fraction(prior_outputs[output_index]["all_recursive_cell_complete_majorant_upper"])
        augmented = prior + tail + roundoff
        output_rows.append(
            {
                "output_index": output_index,
                "output_label": prior_outputs[output_index]["output_label"],
                "prior_all_recursive_cell_complete_majorant_upper": decimal(prior),
                "four_stress_cell_complete_tail_upper": decimal(tail),
                "four_stress_cell_source_roundoff_upper": decimal(roundoff),
                "stress_augmented_complete_majorant_upper": decimal(augmented),
                "margin_below_requested_scale": decimal(REQUESTED_SCALE - augmented),
                "within_requested_scale": augmented < REQUESTED_SCALE,
            }
        )

    maximum_tail = max((Fraction(row["cell_complete_tail_majorant_upper"]), int(row["chain"])) for row in rows)
    maximum_inflation = max((Fraction(row["cell_to_point_tail_inflation_upper"]), int(row["chain"])) for row in rows)
    maximum_rounding_ratio = max(
        (Fraction(row["parent_binary128_rounding"]["maximum_rounding_radius_to_cell_radius_ratio"]), int(row["chain"]))
        for row in rows
    )
    maximum_output = max((Fraction(row["stress_augmented_complete_majorant_upper"]), int(row["output_index"])) for row in output_rows)
    minimum_margin = min((Fraction(row["margin_below_requested_scale"]), int(row["output_index"])) for row in output_rows)
    aggregate = {
        "stress_cell_count": len(rows),
        "precision_ladder_decimal_digits": list(PRECISIONS),
        "precision_overlap_count": sum(bool(row["precision_overlap"]) for row in rows),
        "point_tail_containment_count": sum(bool(row["point_tail_contained_in_cell_tail"]) for row in rows),
        "parent_rounding_preimage_containment_count": sum(bool(row["parent_binary128_rounding"]["all_preimages_inside_cell"]) for row in rows),
        "child_rounding_preimage_containment_count": sum(all(row["child_binary128_rounding_preimages_inside_reconstructed_cell"]) for row in rows),
        "multiplier_rounding_preimage_containment_count": sum(all(row["multiplier_binary128_rounding_preimages_inside_anchored_cell"]) for row in rows),
        "child_state_rounding_preimage_containment_count": sum(all(row["child_state_binary128_rounding_preimages_inside_cell_kernel"]) for row in rows),
        "tpm_rounding_containment_count": sum(bool(row["tpm_rounding_preimage_contains_minus_two_pi"]) for row in rows),
        "source_multiply_add_roundoff_containment_count": sum(
            Fraction(row["source_complex_multiply_add_observed_gap"])
            <= Fraction(row["source_complex_multiply_add_point_roundoff_upper"])
            for row in rows
        ),
        "maximum_local_cell_complete_tail_upper": decimal(maximum_tail[0]),
        "maximum_local_cell_complete_tail_witness": maximum_tail[1],
        "maximum_cell_to_point_tail_inflation_upper": decimal(maximum_inflation[0]),
        "maximum_cell_to_point_tail_inflation_witness": maximum_inflation[1],
        "maximum_parent_rounding_radius_to_cell_radius_ratio": decimal(maximum_rounding_ratio[0]),
        "maximum_parent_rounding_radius_to_cell_radius_witness": maximum_rounding_ratio[1],
        "outputs_with_stress_augmented_majorant_within_requested_scale": sum(bool(row["within_requested_scale"]) for row in output_rows),
        "maximum_stress_augmented_complete_majorant_upper": decimal(maximum_output[0]),
        "maximum_stress_augmented_complete_witness_output": maximum_output[1],
        "minimum_margin_below_requested_scale": decimal(minimum_margin[0]),
        "minimum_margin_witness_output": minimum_margin[1],
    }
    passed = all(
        aggregate[key] == len(rows)
        for key in (
            "precision_overlap_count",
            "point_tail_containment_count",
            "parent_rounding_preimage_containment_count",
            "child_rounding_preimage_containment_count",
            "multiplier_rounding_preimage_containment_count",
            "child_state_rounding_preimage_containment_count",
            "tpm_rounding_containment_count",
            "source_multiply_add_roundoff_containment_count",
        )
    ) and aggregate["outputs_with_stress_augmented_majorant_within_requested_scale"] == 15
    return {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_later_complete_recurrence_stress_cell_pilot_gate",
        "status": "four_later_stress_cells_close_complete_tail_and_pinned_binary128_rounding_pilot" if passed else "later_complete_recurrence_stress_cell_pilot_falsified",
        "scope": "Chains 626, 653, 3406, and 3423 at the saved t=1e10 fixture.",
        "passed": passed,
        "theorem": {
            "factored_phase_remainder": "H-H3=(u^2/y)*z^2*(s^2+4*s+1)/(12*(1+s)^4), s=sqrt(1+z)",
            "common_phase": "The normalized transformed-child phase is factored out at unit modulus before termwise absolute values.",
            "joint_covariance": "With exact W69(M,C)=I69-M*C, M*C+W69 is I69; child and multiplier cell variation is not an independent corrected-model error column.",
            "binary128_preimage": "For round-to-nearest ties-to-even, each finite normal payload represents the exact Voronoi interval between adjacent binary128 midpoints.",
            "source_operation": "The emitted complex recurrence multiply-add is bounded as two rounded products and two rounded additions per real component with unit roundoff 2^-113.",
        },
        "rows": rows,
        "output_rows": output_rows,
        "aggregate": aggregate,
        "runtime": {"priority": priority, "flint_threads": 1, "worker_count": 1},
        "dependencies": {
            "full_later_cells": {"path": relative(FULL_CELLS), "sha256": file_hash(FULL_CELLS)},
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "weights": {"path": relative(WEIGHTS), "sha256": file_hash(WEIGHTS)},
            "block20_cell_tail_builder": {"path": relative(Path(cell_tail.__file__).resolve()), "sha256": file_hash(Path(cell_tail.__file__).resolve())},
            "block20_point_tail_builder": {"path": relative(Path(point_tail.__file__).resolve()), "sha256": file_hash(Path(point_tail.__file__).resolve())},
            "coefficient_transport_builder": {"path": relative(Path(coefficient_transport.__file__).resolve()), "sha256": file_hash(Path(coefficient_transport.__file__).resolve())},
            "instrumented_source": {"path": relative(SOURCE), "sha256": file_hash(SOURCE)},
            "binary128_rounding_mode_probe": {"path": relative(ROUNDING_MODE_PROBE), "sha256": file_hash(ROUNDING_MODE_PROBE)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": "Enclose the source q/special-function cell column, then scale the complete tail and source-operation roundoff to all 1,040 later cells.",
        "proof_boundary": "A four-cell finite saved-height pilot in the pinned source/compiler environment. The source q and special-function variation over cells, the other 1,036 tails, outer Hardy representation, height uniformity, Lambda<=0, PF-infinity, RH, and a prize-level conclusion remain open.",
    }


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy t=1e10 later complete-recurrence stress-cell pilot

Date: 2026-08-09

Status: {artifact['status']}; four-cell finite pilot only, not a proof of the complete recurrence or RH

The maximum-inflation cell (chain 626), narrowest cell (chain 653), and both
longest block-28 cells (chains 3406 and 3423) preserve the complete factored
cubic Legendre remainder at 70 and 110 decimal digits.  The common transformed
child phase is removed at unit modulus before absolute values.

Every parent, reconstructed child, anchored multiplier, child state, and
`tpm=-2*pi` round-to-nearest binary128 preimage lies inside its corresponding
analytic interval.  The emitted source complex multiply-add also lies inside
the explicit `u=2^-113` operation envelope on all four saved calls.

```text
stress cells                                      {aggregate['stress_cell_count']}
precision overlaps                                {aggregate['precision_overlap_count']} / {aggregate['stress_cell_count']}
maximum local complete tail                       {aggregate['maximum_local_cell_complete_tail_upper']}
maximum cell/point tail inflation                 {aggregate['maximum_cell_to_point_tail_inflation_upper']}
maximum parent rounding-radius/cell-radius ratio  {aggregate['maximum_parent_rounding_radius_to_cell_radius_ratio']}
maximum stress-augmented output                   {aggregate['maximum_stress_augmented_complete_majorant_upper']}
minimum remaining 0.005 margin                    {aggregate['minimum_margin_below_requested_scale']}
outputs below 0.005                               {aggregate['outputs_with_stress_augmented_majorant_within_requested_scale']} / 15
```

The exact corrected recurrence still uses the covariance
`W69(M,C)=I69-M*C`; the child and multiplier are not charged twice.

## Boundary

The source-derived probe pins round-to-nearest ties-to-even binary128 in the
accepted Docker/compiler environment.  Source `q` and special-function variation over
the cells, the other 1036 later tails, the outer Hardy representation, height uniformity,
`Lambda<=0`, PF-infinity, RH, and a prize-level conclusion remain open.
"""


def main() -> int:
    artifact = build()
    require(bool(artifact["passed"]), "later complete-recurrence stress-cell pilot did not close")
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    tmp = RESULT.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp.replace(RESULT)
    note_tmp = NOTE.with_suffix(".md.tmp")
    note_tmp.write_text(render_note(artifact), encoding="utf-8")
    note_tmp.replace(NOTE)
    print(
        "built later complete-recurrence stress-cell pilot: "
        f"{artifact['aggregate']['stress_cell_count']} cells, "
        f"{artifact['aggregate']['outputs_with_stress_augmented_majorant_within_requested_scale']}/15 outputs"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
