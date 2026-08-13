#!/usr/bin/env python3
"""Prove the corrected joint recurrence reduces exactly to the endpoint error."""

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
import jensen_window_pf_newman_c1_hardy_block20_exact_eq69_saddle_family_gate as eq69
import jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_full_gate as log_rays
import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_block20_transported_endpoint_budget_gate as transport
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells


ATLAS = coefficient_transport.ATLAS
ENDPOINT_CELLS = coefficient_transport.RESULT
EQ69 = eq69.RESULT
LOG_RAYS = log_rays.RESULT
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_coefficient_cell_joint_recurrence_identity_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_coefficient_cell_joint_recurrence_identity_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_block20_coefficient_cell_joint_recurrence_identity_gate.py"
REQUESTED_ERROR_SCALE = Fraction(5, 1000)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def exact_fraction(record: dict[str, Any]) -> Fraction:
    return Fraction(int(record["numerator"]), int(record["denominator"]))


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
    fixtures = (
        (1, 2, 3, 5, 7, 11),
        (-2, 3, -5, 7, -11, 13),
        (3, -4, 5, -6, 7, -8),
        (Fraction(1, 3), Fraction(2, 5), Fraction(-3, 7), Fraction(5, 11), Fraction(7, 13), Fraction(-11, 17)),
    )
    for parent, integral, multiplier, child, q_exact, q_paper in fixtures:
        parent = Fraction(parent)
        integral = Fraction(integral)
        multiplier = Fraction(multiplier)
        child = Fraction(child)
        q_exact = Fraction(q_exact)
        q_paper = Fraction(q_paper)
        parent = integral + q_exact
        w69 = integral - multiplier * child
        corrected_model = multiplier * child + w69 + q_paper
        require(parent - corrected_model == q_exact - q_paper, "joint recurrence identity drift")
        delta_child = Fraction(7, 19)
        shifted_w69 = integral - multiplier * (child + delta_child)
        shifted_model = multiplier * (child + delta_child) + shifted_w69 + q_paper
        require(shifted_model == corrected_model, "child/W69 covariance drift")
        delta_multiplier = Fraction(-5, 23)
        shifted_multiplier = multiplier + delta_multiplier
        shifted_w69 = integral - shifted_multiplier * child
        shifted_model = shifted_multiplier * child + shifted_w69 + q_paper
        require(shifted_model == corrected_model, "multiplier/W69 covariance drift")
        checks += 3
    return checks


def cell_margins(data: dict[str, Any], atlas_row: dict[str, Any]) -> dict[str, Any]:
    box = coefficient_transport.coefficient_box(data, atlas_row)
    a1_low, a1_high = box["bounds"]["a1"]
    y_low, y_high = box["bounds"]["y"]
    a3_low, a3_high = box["bounds"]["a3"]
    frac_low, frac_high = box["bounds"]["fracL"]
    length = int(data["levels"][1]["length"])
    child_length = int(data["levels"][2]["length"])
    jbot = int(atlas_row["q_cell"]["q_branches"]["jbot"])
    require(length == 104, "parent length drift")
    require(child_length in (1, 2), "child length drift")

    # The atlas is a product box in derived inputs.  The recurrence theorem is
    # used on its physical slice F'(N)=L+fracL; independent coefficient extrema
    # need not themselves satisfy that defining relation.
    xi_high = a1_high + y_high * length + 3 * a3_high * length**2
    require(Fraction(0) < frac_low < frac_high < Fraction(1), "physical-slice dual selector drift")
    if jbot == 1:
        require(a1_low > 0, "positive jbot cell crossed zero")
    else:
        require(a1_high < 0, "zero jbot cell crossed zero")
    require(a1_low > Fraction(-1, 2) and a1_high < Fraction(1, 2), "a1 fundamental cell drift")

    if a3_low > 0:
        curvature_floor = y_low
        cubic_sign = 1
    elif a3_high < 0:
        curvature_floor = y_low + 6 * a3_low * length
        cubic_sign = -1
    else:
        raise RuntimeError("cubic-sign cell drift")
    require(curvature_floor > 0, "curvature floor is not positive")
    upper_gap = Fraction(child_length + 1) - xi_high
    lower_gap = 1 + a1_low
    require(upper_gap > 0 and lower_gap > 0, "nonsaddle mode gap is not positive")
    return {
        "jbot": jbot,
        "child_length": child_length,
        "cubic_sign": cubic_sign,
        "a1_lower_margin": a1_low + Fraction(1, 2),
        "a1_upper_margin": Fraction(1, 2) - a1_high,
        "dual_lower_margin": frac_low,
        "dual_upper_margin": 1 - frac_high,
        "upper_nonsaddle_gap": upper_gap,
        "lower_nonsaddle_gap": lower_gap,
        "curvature_floor": curvature_floor,
    }


def build() -> dict[str, Any]:
    for path in (ATLAS, ENDPOINT_CELLS, EQ69, LOG_RAYS, CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    priority = set_low_priority()
    exact_algebra_checks = algebra_audit()
    atlas = json.loads(ATLAS.read_text(encoding="utf-8"))
    endpoint = json.loads(ENDPOINT_CELLS.read_text(encoding="utf-8"))
    eq69_result = json.loads(EQ69.read_text(encoding="utf-8"))
    ray_result = json.loads(LOG_RAYS.read_text(encoding="utf-8"))
    recursive = native_q.load_recursive_chains()
    atlas_rows = {int(row["chain"]): row for row in atlas["rows"] if int(row["mit"]) == 2}
    eq69_rows = {int(row["chain"]): row for row in eq69_result["rows"]}
    ray_rows = {int(row["chain"]): row for row in ray_result["rows"]}
    require(set(recursive) == set(atlas_rows) == set(eq69_rows) == set(ray_rows), "joint recurrence roster drift")

    rows: list[dict[str, Any]] = []
    for chain in sorted(recursive):
        margins = cell_margins(recursive[chain], atlas_rows[chain])
        require(bool(eq69_rows[chain]["precision_overlap"]), f"chain {chain} equation-(69) precision drift")
        require(bool(ray_rows[chain]["precision_overlap"]), f"chain {chain} ray precision drift")
        require(bool(ray_rows[chain]["formula_target_identity_contains_zero"]), f"chain {chain} point Poisson identity drift")
        rows.append(
            {
                "chain": chain,
                "sum_index": int(recursive[chain]["header"]["sum_index"]),
                "branch": int(recursive[chain]["header"]["branch"]),
                "jbot": margins["jbot"],
                "child_length": margins["child_length"],
                "cubic_sign": margins["cubic_sign"],
                "a1_lower_fundamental_margin": cells.decimal_text(margins["a1_lower_margin"]),
                "a1_upper_fundamental_margin": cells.decimal_text(margins["a1_upper_margin"]),
                "dual_lower_selector_margin": cells.decimal_text(margins["dual_lower_margin"]),
                "dual_upper_selector_margin": cells.decimal_text(margins["dual_upper_margin"]),
                "upper_nonsaddle_gap_lower": cells.decimal_text(margins["upper_nonsaddle_gap"]),
                "lower_nonsaddle_gap_lower": cells.decimal_text(margins["lower_nonsaddle_gap"]),
                "curvature_floor_lower": cells.decimal_text(margins["curvature_floor"]),
                "point_equation_69_precision_overlap": True,
                "point_poisson_identity_contains_zero": True,
            }
        )

    outputs: list[dict[str, Any]] = []
    for source in endpoint["transported_outputs"]:
        bound = Fraction(source["cell_endpoint_majorant_upper"])
        outputs.append(
            {
                "output_index": int(source["output_index"]),
                "output_label": source["output_label"],
                "joint_corrected_recurrence_majorant_upper": source["cell_endpoint_majorant_upper"],
                "majorant_to_requested_scale_ratio": cells.decimal_text(bound / REQUESTED_ERROR_SCALE),
                "within_requested_scale": bound < REQUESTED_ERROR_SCALE,
            }
        )

    def minimum(field: str) -> tuple[Fraction, int]:
        return min((Fraction(row[field]), int(row["chain"])) for row in rows)

    def maximum_output(field: str) -> tuple[Fraction, int]:
        return max((Fraction(row[field]), int(row["output_index"])) for row in outputs)

    min_a1_low = minimum("a1_lower_fundamental_margin")
    min_a1_high = minimum("a1_upper_fundamental_margin")
    min_dual_low = minimum("dual_lower_selector_margin")
    min_dual_high = minimum("dual_upper_selector_margin")
    min_upper_gap = minimum("upper_nonsaddle_gap_lower")
    min_lower_gap = minimum("lower_nonsaddle_gap_lower")
    min_curvature = minimum("curvature_floor_lower")
    max_output = maximum_output("joint_corrected_recurrence_majorant_upper")
    closures = sum(bool(row["within_requested_scale"]) for row in outputs)
    aggregate = {
        "exact_algebra_checks": exact_algebra_checks,
        "cell_count": len(rows),
        "point_equation_69_overlap_count": sum(bool(row["point_equation_69_precision_overlap"]) for row in rows),
        "point_poisson_identity_count": sum(bool(row["point_poisson_identity_contains_zero"]) for row in rows),
        "jbot_histogram": {str(value): sum(int(row["jbot"]) == value for row in rows) for value in (0, 1)},
        "child_length_histogram": {str(value): sum(int(row["child_length"]) == value for row in rows) for value in (1, 2)},
        "cubic_sign_histogram": {str(value): sum(int(row["cubic_sign"]) == value for row in rows) for value in (-1, 1)},
        "minimum_a1_lower_fundamental_margin": cells.decimal_text(min_a1_low[0]),
        "minimum_a1_lower_fundamental_margin_witness": min_a1_low[1],
        "minimum_a1_upper_fundamental_margin": cells.decimal_text(min_a1_high[0]),
        "minimum_a1_upper_fundamental_margin_witness": min_a1_high[1],
        "minimum_dual_lower_selector_margin": cells.decimal_text(min_dual_low[0]),
        "minimum_dual_lower_selector_margin_witness": min_dual_low[1],
        "minimum_dual_upper_selector_margin": cells.decimal_text(min_dual_high[0]),
        "minimum_dual_upper_selector_margin_witness": min_dual_high[1],
        "minimum_upper_nonsaddle_gap": cells.decimal_text(min_upper_gap[0]),
        "minimum_upper_nonsaddle_gap_witness": min_upper_gap[1],
        "minimum_lower_nonsaddle_gap": cells.decimal_text(min_lower_gap[0]),
        "minimum_lower_nonsaddle_gap_witness": min_lower_gap[1],
        "minimum_curvature_floor": cells.decimal_text(min_curvature[0]),
        "minimum_curvature_floor_witness": min_curvature[1],
        "maximum_joint_corrected_recurrence_majorant_upper": cells.decimal_text(max_output[0]),
        "maximum_joint_corrected_recurrence_majorant_witness": max_output[1],
        "remaining_requested_scale_slack_lower": cells.decimal_text(REQUESTED_ERROR_SCALE - max_output[0]),
        "outputs_within_requested_scale": closures,
    }
    require(closures == 15, "joint corrected recurrence budget did not close")
    return {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_coefficient_cell_joint_recurrence_identity_gate",
        "status": "joint_corrected_recurrence_reduces_exactly_to_cell_endpoint_error_on_all_374_physical_cell_slices",
        "scope": (
            "The mathematically corrected exact-W1 recurrence on the physical slices "
            "F'(N)=L+fracL of all 374 occupied block-20 coefficient cells."
        ),
        "exact_identity": {
            "poisson_split": "P=I69+Q_exact",
            "exact_w1": "W69(M,C)=I69-M*C",
            "corrected_model": "R_corr=M*C+W69(M,C)+Q_paper=I69+Q_paper",
            "error": "P-R_corr=Q_exact-Q_paper",
            "covariance": "For arbitrary delta C or delta M, the induced change in M*C is cancelled exactly by W69=I69-M*C.",
            "consequence": "Parent, transformed-child, multiplier, and equation-(69) saddle variations are not independent error columns in the corrected joint model.",
        },
        "analytic_contract": {
            "finite_poisson": "Endpoint-complete finite Poisson summation for the entire cubic exponential on integer [0,N].",
            "mode_split": "jbot=ceil(Phi1), L=floor(F'(N)); saddle modes jbot..L and nonsaddle complements are fixed on each cell.",
            "nonsaddle_reduction": "The complements are the four logarithmic endpoint rays plus the explicit half-endpoint and optional zero mode of Formal Core 11.264.",
            "uniformity": (
                "On each physical slice F'(N)=L+fracL, positive fracL margins fix the dual selector; "
                "the independent coefficient bounds certify the endpoint majorant, nonsaddle gaps, curvature floor, and cubic sign."
            ),
            "normalization": "This identity uses mathematical 2*pi. Source binary normalization and implementation arithmetic remain separate."
        },
        "rows": rows,
        "transported_outputs": outputs,
        "aggregate": aggregate,
        "runtime": {"priority": priority, "active_workers": 1},
        "dependencies": {
            "selector_atlas": {"path": relative(ATLAS), "sha256": file_hash(ATLAS)},
            "endpoint_cells": {"path": relative(ENDPOINT_CELLS), "sha256": file_hash(ENDPOINT_CELLS)},
            "equation_69": {"path": relative(EQ69), "sha256": file_hash(EQ69)},
            "logarithmic_rays": {"path": relative(LOG_RAYS), "sha256": file_hash(LOG_RAYS)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": (
            "Transport source-to-corrected-model arithmetic separately: exact W1 replacement, the source-only "
            "t2 L+1 deletion, special-function normalization, binary tpm, and floating recurrence arithmetic."
        ),
        "proof_boundary": (
            "Exact corrected-model identity and finite block-20 coefficient-cell endpoint budget only. It does not "
            "prove that the unmodified source computes exact W1 or Q_paper, does not close source arithmetic, other "
            "blocks, height uniformity, cross-block accumulation, the outer Hardy representation, Lambda<=0, "
            "PF-infinity, RH, or a prize-level conclusion."
        ),
    }


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 coefficient-cell joint recurrence identity

Date: 2026-08-09
Status: {artifact['status']}; corrected finite model only, not a proof of the source evaluator or RH

## Correlated Identity

Let `P` be the exact parent sum, `I69` the exact equation-(69) saddle family,
`M*C` any transformed-child model, and

```text
W69(M,C)=I69-M*C.
```

Endpoint-complete finite Poisson summation gives `P=I69+Q_exact`.  Therefore

```text
M*C+W69(M,C)+Q_paper = I69+Q_paper,
P-[M*C+W69(M,C)+Q_paper] = Q_exact-Q_paper.
```

This is invariant under arbitrary changes of `M` or `C` when `W69` is changed
with them.  The broad parent, child, multiplier, and saddle Lipschitz columns
must therefore not be added independently in the corrected joint model.

## Cell Guards

On the physical slice `F'(N)=L+fracL`, all {aggregate['cell_count']} occupied
cells preserve `jbot`, `L`, cubic sign, positive curvature, and positive
upper/lower nonsaddle gaps.  The enclosing product boxes are deliberately
wider in their four derived coordinates; their independent coefficient bounds
certify the endpoint majorant.  The 374 saved equation-(69) integrations
overlap across precision, and all 374 independent point Poisson identities
contain zero.

```text
minimum dual lower selector margin      {aggregate['minimum_dual_lower_selector_margin']}
minimum dual upper selector margin      {aggregate['minimum_dual_upper_selector_margin']}
minimum upper nonsaddle gap             {aggregate['minimum_upper_nonsaddle_gap']}
minimum lower nonsaddle gap             {aggregate['minimum_lower_nonsaddle_gap']}
minimum curvature floor                 {aggregate['minimum_curvature_floor']}
maximum corrected recurrence budget     {aggregate['maximum_joint_corrected_recurrence_majorant_upper']}
remaining 0.005 slack                   {aggregate['remaining_requested_scale_slack_lower']}
outputs below 0.005                     {aggregate['outputs_within_requested_scale']} / 15
```

The corrected joint recurrence budget is exactly the already certified
coefficient-cell endpoint budget.  No Legendre-tail column is added because a
change in `M*C` is cancelled by the defining exact W1 correction.

## Boundary

The unmodified source does not compute exact `W69` or exact `Q_paper` by
definition.  Its W1 replacement, source-only `t2` term, special functions,
binary normalization, and floating arithmetic must now be transported as a
separate source-to-corrected-model column.  Other blocks, uniform height,
cross-block accumulation, the outer Hardy representation, `Lambda<=0`,
PF-infinity, RH, and a prize-level theorem remain open.
"""


def main() -> int:
    artifact = build()
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built coefficient-cell joint recurrence identity: "
        f"{artifact['aggregate']['cell_count']} cells, "
        f"{artifact['aggregate']['outputs_within_requested_scale']}/15 outputs below 0.005"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
