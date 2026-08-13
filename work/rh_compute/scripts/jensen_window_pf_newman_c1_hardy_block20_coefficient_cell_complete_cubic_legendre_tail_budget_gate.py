#!/usr/bin/env python3
"""Transport the complete cubic Legendre tail over all block-20 selector cells."""

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

import jensen_window_pf_newman_c1_hardy_block20_coefficient_neighborhood_transport_gate as coefficient_transport
import jensen_window_pf_newman_c1_hardy_block20_complete_cubic_legendre_tail_budget_gate as point_tail
import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_block20_transported_endpoint_budget_gate as transport
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import acb, arb, ctx


ATLAS = coefficient_transport.ATLAS
ENDPOINT_CELLS = coefficient_transport.RESULT
POINT_TAIL = point_tail.RESULT
WEIGHT_FIXTURE = transport.WEIGHT_FIXTURE
WEIGHTS = transport.WEIGHTS
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_coefficient_cell_complete_cubic_legendre_tail_budget_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_coefficient_cell_complete_cubic_legendre_tail_budget_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_block20_coefficient_cell_complete_cubic_legendre_tail_budget_gate.py"
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


def algebra_audit() -> int:
    checks = 0
    for s in (Fraction(1, 2), Fraction(2, 3), Fraction(1), Fraction(4, 3), Fraction(2)):
        exact_factor = 2 * (1 + 2 * s) / (3 * (1 + s) ** 2)
        cubic_factor = Fraction(7, 12) - s**2 / 12
        factored = (s - 1) ** 2 * (s**2 + 4 * s + 1) / (12 * (1 + s) ** 2)
        require(exact_factor - cubic_factor == factored, "phase-remainder factorization drift")
        z = s**2 - 1
        z_factored = z**2 * (s**2 + 4 * s + 1) / (12 * (1 + s) ** 4)
        require(factored == z_factored, "z-factorization drift")
        checks += 2
    return checks


def evaluate_cell(data: dict[str, Any], atlas_row: dict[str, Any], dps: int) -> dict[str, Any]:
    ctx.dps = dps
    ctx.threads = 1
    box = coefficient_transport.coefficient_box(data, atlas_row)
    structural = coefficient_transport.structural_transport(data, atlas_row, box, dps)
    center = box["center"]
    radii = box["radii"]
    a1 = coefficient_transport.interval(center["a1"], radii["a1"])
    y = coefficient_transport.interval(center["y"], radii["y"])
    a3 = coefficient_transport.interval(center["a3"], radii["a3"])
    require(lower(y) > 0, "cell quadratic coefficient is not positive")
    tpm = point_balls.binary128_ball(data["header"]["tpm_hex"])

    term_gaps: list[arb] = []
    elementary_gaps: list[arb] = []
    tail_values: list[arb] = []
    z_values: list[arb] = []
    discriminants: list[arb] = []
    exact_denominators: list[arb] = []
    source_denominators: list[arb] = []
    for index in range(int(data["levels"][2]["length"]) + 1):
        u = index - a1
        z = 12 * a3 * u / y**2
        discriminant = 1 + z
        require(lower(discriminant) > 0, "cell stationary discriminant is not positive")
        s = discriminant.sqrt()
        tail = u**2 / y * z**2 * (s**2 + 4 * s + 1) / (12 * (1 + s) ** 4)
        require(lower(tail) >= 0, "factored phase remainder lost positivity")
        exact_denominator = s.sqrt()
        source_denominator = 1 + z / 4 - 3 * z**2 / 32
        require(lower(exact_denominator) > 0, "exact amplitude denominator is not positive")
        require(lower(source_denominator) > 0, "source amplitude denominator is not positive")
        reduced_gap = acb(arb(0), tpm * tail).exp() / exact_denominator - 1 / source_denominator
        direct_gap = abs(reduced_gap)
        amplitude_gap = abs(1 / exact_denominator - 1 / source_denominator)
        elementary_gap = amplitude_gap + abs(tpm) * tail / source_denominator
        require(upper(direct_gap) <= upper(elementary_gap), "cell direct term escaped elementary majorant")
        term_gaps.append(direct_gap)
        elementary_gaps.append(elementary_gap)
        tail_values.append(tail)
        z_values.append(abs(z))
        discriminants.append(discriminant)
        exact_denominators.append(exact_denominator)
        source_denominators.append(source_denominator)

    child_majorant = sum((upper(value) for value in term_gaps), Fraction(0))
    elementary_majorant = sum((upper(value) for value in elementary_gaps), Fraction(0))
    multiplier_upper = upper_abs(structural["multiplier"])
    return {
        "term_gaps": term_gaps,
        "child_majorant": child_majorant,
        "elementary_majorant": elementary_majorant,
        "multiplier_upper": multiplier_upper,
        "local_majorant": multiplier_upper * child_majorant,
        "local_elementary_majorant": multiplier_upper * elementary_majorant,
        "maximum_tail": max(upper(value) for value in tail_values),
        "maximum_z": max(upper(value) for value in z_values),
        "minimum_discriminant": min(lower(value) for value in discriminants),
        "minimum_exact_denominator": min(lower(value) for value in exact_denominators),
        "minimum_source_denominator": min(lower(value) for value in source_denominators),
    }


def build() -> dict[str, Any]:
    for path in (ATLAS, ENDPOINT_CELLS, POINT_TAIL, WEIGHT_FIXTURE, WEIGHTS, CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    priority = set_low_priority()
    exact_algebra_checks = algebra_audit()
    atlas = json.loads(ATLAS.read_text(encoding="utf-8"))
    atlas_rows = {int(row["chain"]): row for row in atlas["rows"] if int(row["mit"]) == 2}
    endpoint_result = json.loads(ENDPOINT_CELLS.read_text(encoding="utf-8"))
    endpoint_outputs = {int(row["output_index"]): row for row in endpoint_result["transported_outputs"]}
    point_result = json.loads(POINT_TAIL.read_text(encoding="utf-8"))
    point_rows = {int(row["chain"]): row for row in point_result["rows"]}
    recursive = native_q.load_recursive_chains()
    require(set(atlas_rows) == set(point_rows) == set(recursive), "cell-tail roster mismatch")

    rows: list[dict[str, Any]] = []
    for chain in sorted(recursive):
        low = evaluate_cell(recursive[chain], atlas_rows[chain], PRECISIONS[0])
        high = evaluate_cell(recursive[chain], atlas_rows[chain], PRECISIONS[1])
        require(len(low["term_gaps"]) == len(high["term_gaps"]), "term roster precision drift")
        require(all(a.overlaps(b) for a, b in zip(low["term_gaps"], high["term_gaps"])), f"chain {chain} precision nonoverlap")
        point_upper = Fraction(point_rows[chain]["local_recurrence_tail_majorant_upper"])
        require(point_upper <= high["local_majorant"], f"chain {chain} saved point escaped cell tail")
        header = recursive[chain]["header"]
        rows.append(
            {
                "chain": chain,
                "sum_index": int(header["sum_index"]),
                "branch": int(header["branch"]),
                "child_length": int(recursive[chain]["levels"][2]["length"]),
                "precision_overlap": True,
                "maximum_dimensionless_z_abs_upper": point_tail.upper_decimal(high["maximum_z"]),
                "minimum_stationary_discriminant_lower": point_tail.lower_decimal(high["minimum_discriminant"]),
                "minimum_exact_fourth_root_denominator_lower": point_tail.lower_decimal(high["minimum_exact_denominator"]),
                "minimum_source_taylor_denominator_lower": point_tail.lower_decimal(high["minimum_source_denominator"]),
                "maximum_factored_legendre_phase_tail_upper": point_tail.upper_decimal(high["maximum_tail"]),
                "multiplier_cell_abs_upper": point_tail.upper_decimal(high["multiplier_upper"]),
                "child_cell_termwise_direct_gap_upper": point_tail.upper_decimal(high["child_majorant"]),
                "child_cell_termwise_elementary_majorant_upper": point_tail.upper_decimal(high["elementary_majorant"]),
                "local_point_tail_majorant_upper": point_tail.upper_decimal(point_upper),
                "local_cell_tail_majorant_upper": point_tail.upper_decimal(high["local_majorant"]),
                "local_cell_tail_elementary_majorant_upper": point_tail.upper_decimal(high["local_elementary_majorant"]),
                "cell_to_point_tail_inflation_upper": point_tail.upper_decimal(high["local_majorant"] / point_upper),
                "term_gap_abs_uppers": [point_tail.upper_decimal(upper(value)) for value in high["term_gaps"]],
            }
        )

    ctx.dps = WEIGHT_PRECISION
    ctx.threads = 1
    weights = transport.load_weights()
    transported_outputs: list[dict[str, Any]] = []
    for output_index in range(1, 16):
        tail_majorant = Fraction(0)
        for row in rows:
            weight = weights[(int(row["sum_index"]), output_index, int(row["branch"]))]
            tail_majorant += (
                weight["amplitude_fraction"]
                * weight["phase_abs_upper"]
                * Fraction(row["local_cell_tail_majorant_upper"])
            )
        endpoint = Fraction(endpoint_outputs[output_index]["cell_endpoint_majorant_upper"])
        combined = endpoint + tail_majorant
        transported_outputs.append(
            {
                "output_index": output_index,
                "output_label": transport.output_label(output_index),
                "coefficient_cell_endpoint_majorant_upper": point_tail.upper_decimal(endpoint),
                "coefficient_cell_legendre_tail_majorant_upper": point_tail.upper_decimal(tail_majorant),
                "combined_cell_majorant_upper": point_tail.upper_decimal(combined),
                "combined_to_requested_scale_ratio": point_tail.upper_decimal(combined / REQUESTED_ERROR_SCALE),
                "combined_within_requested_scale": combined < REQUESTED_ERROR_SCALE,
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
    max_tail = maximum("maximum_factored_legendre_phase_tail_upper")
    max_inflation = maximum("cell_to_point_tail_inflation_upper")
    max_local = maximum("local_cell_tail_majorant_upper")
    max_output_tail = output_maximum("coefficient_cell_legendre_tail_majorant_upper")
    max_combined = output_maximum("combined_cell_majorant_upper")
    closures = sum(bool(row["combined_within_requested_scale"]) for row in transported_outputs)
    aggregate = {
        "exact_algebra_checks": exact_algebra_checks,
        "cell_count": len(rows),
        "child_index_count": sum(len(row["term_gap_abs_uppers"]) for row in rows),
        "precision_ladder_decimal_digits": list(PRECISIONS),
        "precision_overlap_count": sum(bool(row["precision_overlap"]) for row in rows),
        "maximum_dimensionless_z_abs_upper": point_tail.upper_decimal(max_z[0]),
        "maximum_dimensionless_z_abs_witness": max_z[1],
        "minimum_stationary_discriminant_lower": point_tail.lower_decimal(min_disc[0]),
        "minimum_stationary_discriminant_witness": min_disc[1],
        "maximum_factored_legendre_phase_tail_upper": point_tail.upper_decimal(max_tail[0]),
        "maximum_factored_legendre_phase_tail_witness": max_tail[1],
        "maximum_cell_to_point_tail_inflation_upper": point_tail.upper_decimal(max_inflation[0]),
        "maximum_cell_to_point_tail_inflation_witness": max_inflation[1],
        "maximum_local_cell_tail_majorant_upper": point_tail.upper_decimal(max_local[0]),
        "maximum_local_cell_tail_majorant_witness": max_local[1],
        "maximum_transported_cell_tail_majorant_upper": point_tail.upper_decimal(max_output_tail[0]),
        "maximum_transported_cell_tail_majorant_witness": max_output_tail[1],
        "maximum_combined_cell_majorant_upper": point_tail.upper_decimal(max_combined[0]),
        "maximum_combined_cell_majorant_witness": max_combined[1],
        "requested_error_scale": cells.decimal_text(REQUESTED_ERROR_SCALE),
        "outputs_with_combined_cell_majorant_within_requested_scale": closures,
    }
    status = (
        "complete_cubic_legendre_tail_transported_over_all_selector_cells_and_combined_budget_closes"
        if closures == 15
        else "complete_cubic_legendre_tail_cell_transport_breaks_combined_budget"
    )
    return {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_coefficient_cell_complete_cubic_legendre_tail_budget_gate",
        "status": status,
        "scope": "All 374 occupied block-20 selector cells and their 753 transformed-child indices.",
        "theorem": {
            "factored_phase_remainder": "H-H3=(u^2/y)*z^2*(s^2+4*s+1)/(12*(1+s)^4), s=sqrt(1+z)",
            "positivity": "The factored phase remainder is nonnegative whenever y>0 and 1+z>0.",
            "phase_cancellation": "The common normalized child phase has unit modulus and factors out before bounding each term.",
            "reduced_term": "|exp(i*tpm*(H-H3))/(1+z)^(1/4)-1/(1+z/4-3*z^2/32)|",
            "multiplier": "The full-cell source-anchored recurrence multiplier enclosure is applied after the termwise child sum.",
        },
        "rows": rows,
        "transported_outputs": transported_outputs,
        "aggregate": aggregate,
        "dependencies": {
            "selector_atlas": {"path": relative(ATLAS), "sha256": file_hash(ATLAS)},
            "endpoint_cells": {"path": relative(ENDPOINT_CELLS), "sha256": file_hash(ENDPOINT_CELLS)},
            "point_tail": {"path": relative(POINT_TAIL), "sha256": file_hash(POINT_TAIL)},
            "weight_fixture": {"path": relative(WEIGHT_FIXTURE), "sha256": file_hash(WEIGHT_FIXTURE)},
            "weight_rows": {"path": relative(WEIGHTS), "sha256": file_hash(WEIGHTS)},
        },
        "runtime": {"priority": priority, "flint_threads": 1},
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": (
            "Certify the remaining parent/saddle and source-arithmetic recurrence columns on accuracy-refined "
            "subcells, then extend the endpoint-plus-tail theorem to other blocks and the outer Hardy remainder."
        ),
        "proof_boundary": (
            "Rigorous finite block-20 coefficient-cell theorem for the endpoint and complete cubic Legendre-tail "
            "columns only. It does not close parent/saddle recurrence arithmetic, source rounding, other blocks, "
            "cross-block accumulation, the outer Hardy representation, Lambda<=0, PF-infinity, RH, or a prize-level conclusion."
        ),
    }


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 coefficient-cell complete cubic Legendre-tail budget

Date: 2026-08-09
Status: {artifact['status']}; finite block-20 theorem only, not a proof of the complete evaluator or RH

## Stable Exact Remainder

For `s=sqrt(1+z)`, exact algebra factors the complete phase remainder as

```text
H-H3=(u^2/y) z^2 (s^2+4s+1)/[12(1+s)^4].
```

This form is nonnegative on every admitted cell and avoids subtracting two
nearly equal interval phases.  The common normalized child phase has unit
modulus and is removed before estimation.  Each term is therefore enclosed as

```text
|exp(i*tpm*(H-H3))/(1+z)^(1/4)
 -1/(1+z/4-3z^2/32)|.
```

The full-cell recurrence multiplier is applied only after the termwise child
sum, and the exact output weights are applied only after the local majorant.

## Result

```text
selector cells                           {aggregate['cell_count']}
child indices                            {aggregate['child_index_count']}
precision overlaps                       {aggregate['precision_overlap_count']} / {aggregate['cell_count']}
maximum cell |z|                         {aggregate['maximum_dimensionless_z_abs_upper']}
minimum cell 1+z                         {aggregate['minimum_stationary_discriminant_lower']}
maximum cell phase tail                  {aggregate['maximum_factored_legendre_phase_tail_upper']}
maximum cell/point tail inflation        {aggregate['maximum_cell_to_point_tail_inflation_upper']}
maximum transported cell tail            {aggregate['maximum_transported_cell_tail_majorant_upper']}
maximum endpoint plus cell tail          {aggregate['maximum_combined_cell_majorant_upper']}
outputs below 0.005                      {aggregate['outputs_with_combined_cell_majorant_within_requested_scale']} / 15
```

No signed output cancellation is used.

## Boundary

This closes only the coefficient-cell endpoint and complete cubic Legendre-tail
columns at finite block 20.  Parent/saddle recurrence arithmetic, source
rounding, other blocks, cross-block accumulation, the outer Hardy
representation, `Lambda<=0`, PF-infinity, RH, and a prize-level theorem remain
open.
"""


def main() -> int:
    artifact = build()
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built coefficient-cell complete cubic Legendre-tail budget: "
        f"{artifact['aggregate']['cell_count']} cells, "
        f"{artifact['aggregate']['outputs_with_combined_cell_majorant_within_requested_scale']}/15 outputs below 0.005"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
