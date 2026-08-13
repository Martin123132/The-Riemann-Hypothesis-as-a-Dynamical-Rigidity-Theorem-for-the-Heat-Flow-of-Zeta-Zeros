#!/usr/bin/env python3
"""Assemble the exact finite source-aligned Hardy hybrid-error ledger."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx


RESULT_ROOT = REPO_ROOT / "work/rh_compute/results"
CALIBRATION = RESULT_ROOT / "jensen_window_pf_newman_c1_hardy_t1e10_exact_hardy_arb_calibration_gate.json"
COMPONENT_SPLIT = RESULT_ROOT / "jensen_window_pf_newman_c1_hardy_t1e10_lower_rs_hybrid_upper_component_split_gate.json"
CLASSICAL_SPLIT = RESULT_ROOT / "jensen_window_pf_newman_c1_hardy_t1e10_classical_upper_main_remainder_split_gate.json"
BLOCK_PARTITION = RESULT_ROOT / "jensen_window_pf_newman_c1_hardy_t1e10_equation124_classical_block_partition_gate.json"
BLOCK_ZERO_IDENTITY = RESULT_ROOT / "jensen_window_pf_newman_c1_hardy_t1e10_block_zero_source_formula_identity_gate.json"
OUTER_LEDGER = RESULT_ROOT / "jensen_window_pf_newman_c1_hardy_t1e10_outer_hardy_representation_dependency_ledger_gate.json"
CORRECTED_INTERNAL = RESULT_ROOT / "jensen_window_pf_newman_c1_hardy_t1e10_crossblock_hybrid_exact_output_gate.json"
SOURCE_Q_CUT = RESULT_ROOT / "jensen_window_pf_newman_c1_hardy_t1e10_rigorous_arb_source_q_dependency_cut_gate.json"
PAPER = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_Brereton_2607.15310.pdf"
SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/telemetry/zeta14cubicmult_telemetry.f90"
RESULT = RESULT_ROOT / "jensen_window_pf_newman_c1_hardy_t1e10_source_aligned_hybrid_error_ledger_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_source_aligned_hybrid_error_ledger_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 110
EXPECTED_OUTPUTS = 15
EXPECTED_BLOCKS = 36
EXPECTED_BLOCK_ROWS = EXPECTED_OUTPUTS * EXPECTED_BLOCKS
TOLERANCE = arb("0.005")
FIRST_OPEN_ID = "source_aligned_upper_block_sum"

EXPECTED_KINDS = {
    CALIBRATION: "jensen_window_pf_newman_c1_hardy_t1e10_exact_hardy_arb_calibration_gate",
    COMPONENT_SPLIT: "jensen_window_pf_newman_c1_hardy_t1e10_lower_rs_hybrid_upper_component_split_gate",
    CLASSICAL_SPLIT: "jensen_window_pf_newman_c1_hardy_t1e10_classical_upper_main_remainder_split_gate",
    BLOCK_PARTITION: "jensen_window_pf_newman_c1_hardy_t1e10_equation124_classical_block_partition_gate",
    BLOCK_ZERO_IDENTITY: "jensen_window_pf_newman_c1_hardy_t1e10_block_zero_source_formula_identity_gate",
    OUTER_LEDGER: "jensen_window_pf_newman_c1_hardy_t1e10_outer_hardy_representation_dependency_ledger_gate",
    CORRECTED_INTERNAL: "jensen_window_pf_newman_c1_hardy_t1e10_crossblock_hybrid_exact_output_gate",
    SOURCE_Q_CUT: "jensen_window_pf_newman_c1_hardy_t1e10_rigorous_arb_source_q_dependency_cut_gate",
}

SOURCE_ASSERTIONS = (
    {
        "id": "source_hybrid_split",
        "lines": (656, 669),
        "tokens": ("ZP(t)+RSremainder(t).", "variable zsum", "variable rszsum"),
    },
    {
        "id": "block_zero_transition_and_direct_terms",
        "lines": (675, 702),
        "tokens": ("call start(a,t6,yphase,transit,M2)", "call ter(ial,1/aar(i),yphasear(i),tar(i),term)", "zsum(i)=zsum(i)*sqrt(8/aar(i))"),
    },
    {
        "id": "cubic_blocks_one_through_totblock",
        "lines": (778, 850),
        "tokens": ("do iblock=rh_block_start,totblock", "call alphasum(MTM,rae,pcsum)", "zsum=zsum+zsum2"),
    },
    {
        "id": "source_endpoint_cutoff",
        "lines": (892, 924),
        "tokens": ("raecutoff=RN1", "CE=raecutoff*(1-sqrt(1-(a/(raecutoff*C))**2))/CE", "NC=aint(CE,dp1)"),
    },
    {
        "id": "lower_rs_sum_and_leading_endpoint_correction",
        "lines": (921, 1035),
        "tokens": ("CALL pari_calc", "RSN=aint(CE,dp1)", "rszsum(i)=2*rszsum(i)+((-1)**(RSN-1))*CE*C"),
    },
    {
        "id": "final_source_addition",
        "lines": (1062, 1068),
        "tokens": ("zsum(i)=zsum(i)+rszsum(i)", "rh_stage_complete"),
    },
    {
        "id": "direct_or_recursive_cubic_sums",
        "lines": (1521, 1565),
        "tokens": ("if (MTM.lt.K) then", "call computarrs", "call mgausssum", "pcsum(i)=amp(i)*(real(c1*pcgsum(1)+c2*pcgsum(2)))"),
    },
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        result = json.load(handle)
    require(result.get("kind") == EXPECTED_KINDS[path], f"kind drift: {relative(path)}")
    if path == CORRECTED_INTERNAL:
        require(
            result.get("status") == "rigorous_t1e10_complete_6784_call_endpoint_direct_exact_point_column_below_0_005",
            "legacy corrected-internal status drift",
        )
        require(result["aggregate"].get("outputs_within_requested_scale") == EXPECTED_OUTPUTS, "legacy corrected-internal aggregate drift")
    else:
        require(result.get("passed") is True, f"upstream gate is not passed: {relative(path)}")
    return result


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


def sign_name(value: arb) -> str:
    if value > 0:
        return "positive"
    if value < 0:
        return "negative"
    return "contains_zero"


def minimum(values: list[arb]) -> arb:
    return min(values, key=lambda value: value.lower())


def maximum(values: list[arb]) -> arb:
    return max(values, key=lambda value: value.upper())


def audit_source() -> list[dict[str, Any]]:
    source_lines = SOURCE.read_text(encoding="utf-8", errors="replace").splitlines()
    rows: list[dict[str, Any]] = []
    for assertion in SOURCE_ASSERTIONS:
        audit_id = assertion["id"]
        start, end = assertion["lines"]
        excerpt = "\n".join(source_lines[start - 1 : end])
        missing = [token for token in assertion["tokens"] if token not in excerpt]
        require(not missing, f"source audit {audit_id} missing tokens: {missing}")
        digest = hashlib.sha256(excerpt.encode("utf-8")).hexdigest()
        rows.append(
            {
                "id": audit_id,
                "line_start": start,
                "line_end": end,
                "required_tokens": list(assertion["tokens"]),
                "missing_token_count": 0,
                "excerpt_sha256": digest,
            }
        )
    return rows


def channel_ledger(corrected: dict[str, Any], q_cut: dict[str, Any]) -> list[dict[str, Any]]:
    corrected_bound = corrected["aggregate"]["maximum_hybrid_complete_majorant_upper"]
    q_cut_bound = q_cut["aggregate"]["maximum_source_q_free_complete_majorant_upper"]
    return [
        {
            "id": "published_hybrid_identity_and_source_control_flow",
            "state": "source_audited",
            "scope": "Equations (124)--(127), the endpoint cutoff, lower sum, hybrid upper blocks, endpoint correction, and final addition.",
            "gate": relative(OUTER_LEDGER),
            "limitation": "The paper explicitly says the hybrid representation is not exact and does not supply usable constants for every error term here.",
        },
        {
            "id": "source_endpoint_cutoff_partition",
            "state": "interval_certified_saved_heights",
            "scope": "All 36 equation-(124) endpoint cutoffs and all 39,273 classical upper-main terms at the fifteen saved heights.",
            "gate": relative(BLOCK_PARTITION),
            "limitation": "Finite saved-height partition, not a height-uniform error theorem.",
        },
        {
            "id": "exact_hardy_target",
            "state": "interval_certified_saved_heights",
            "scope": "Independent Arb Hardy Z values at all fifteen saved heights.",
            "gate": relative(CALIBRATION),
            "limitation": "Finite calibration only.",
        },
        {
            "id": "lower_rs_source_error",
            "state": "interval_certified_saved_heights",
            "scope": "E_lower=L_source-L_exact, including the source leading endpoint correction.",
            "gate": relative(COMPONENT_SPLIT),
            "limitation": "No uniform source-arithmetic and Riemann--Siegel remainder theorem has been promoted.",
        },
        {
            "id": "exact_upper_remaining_correction",
            "state": "interval_certified_saved_heights",
            "scope": "R_upper=U_exact-T_upper.",
            "gate": relative(CLASSICAL_SPLIT),
            "limitation": "Its small saved-height size is not a uniform bound.",
        },
        {
            "id": "final_source_addition_rounding",
            "state": "interval_certified_saved_heights",
            "scope": "rho_add=Z_source-(L_source+Q_source).",
            "gate": relative(COMPONENT_SPLIT),
            "limitation": "Finite saved source values only.",
        },
        {
            "id": "source_aligned_block_increments",
            "state": "interval_certified_saved_heights",
            "scope": "Every e_k=Delta Q_k-Delta T_k for k=0,...,35, with exact telescoping to Q_source-T_upper.",
            "gate": relative(BLOCK_PARTITION),
            "limitation": "The 540 finite intervals do not prove a signed estimate for arbitrary height.",
        },
        {
            "id": "corrected_internal_model_column",
            "state": "interval_certified_but_not_source_aligned_outer_error",
            "scope": f"A corrected internal endpoint/direct column is below {corrected_bound}; its source-q-free successor is below {q_cut_bound}.",
            "gate": relative(SOURCE_Q_CUT),
            "limitation": "Neither result is the legacy source-aligned equations-(126)--(127) error and neither may close this ledger.",
        },
        {
            "id": FIRST_OPEN_ID,
            "state": "open_height_uniform_signed_theorem",
            "scope": "E_hyb_main(t)=sum_k e_k(t)=Q_source(t)-T_upper(t), with the source endpoint convention retained.",
            "required_upgrade": "Derive a constant-bearing height-uniform signed estimate before taking blockwise absolute values; preserve transition, endpoint, direct/recursive, and cross-block cancellation.",
        },
        {
            "id": "source_aligned_subordinate_channels",
            "state": "open_height_uniform_bounds",
            "scope": "E_lower(t), R_upper(t), and rho_add(t).",
            "required_upgrade": "Promote explicit height-uniform Riemann--Siegel, source arithmetic, and final-rounding bounds in the same parameter domain.",
        },
        {
            "id": "assembled_source_aligned_hardy_error",
            "state": "open_height_uniform_theorem",
            "scope": "E_total(t)=E_lower(t)+E_hyb_main(t)-R_upper(t)+rho_add(t).",
            "required_upgrade": "Combine the signed hybrid theorem and subordinate uniform bounds without importing the corrected-model column as a source-error estimate.",
        },
    ]


def render_note(artifact: dict[str, Any]) -> str:
    agg = artifact["aggregate"]
    center = artifact["output_rows"][7]
    return f"""# Source-aligned hybrid-error ledger at t=10^10

Date: 2026-08-10

Status: exact finite interval identity validated; not a proof of the height-uniform signed theorem

The primary-source route uses the endpoint cutoff from equations (124)--(127)
and the Fortran control flow.  The diagnostic nearest-root midpoint is not used
in this ledger.

Write `Z_src=L_src+Q_src+rho_add`, `Z_exact=L_exact+U_exact`, and let
`T_upper` be the exact classical upper-main sum induced by the 36 source
cutoffs.  For block increments

```text
e_k = (Q_k-Q_(k-1)) - (T_k-T_(k-1)),  k=0,...,35,
```

the independently checked finite identity is

```text
Z_src-Z_exact
  = (L_src-L_exact) + sum_(k=0)^35 e_k
    - (U_exact-T_upper) + rho_add.
```

All {agg['output_count']} output identities and all {agg['block_row_count']} block
rows are interval certified.  At the central output,

```text
E_total    = {center['total_source_minus_exact_ball']}
E_lower    = {center['lower_source_minus_exact_ball']}
sum e_k    = {center['source_aligned_upper_block_sum_ball']}
R_upper    = {center['exact_upper_remaining_correction_ball']}
rho_add    = {center['final_addition_rounding_ball']}
sum |e_k|  = {center['independent_block_triangle_ball']}
```

The triangle sum rejects `0.005` on all {agg['block_triangle_rejects_0p005_count']}
outputs.  Its range is

```text
{agg['minimum_independent_block_triangle_ball']}
to
{agg['maximum_independent_block_triangle_ball']}.
```

Only {agg['minimum_retained_fraction_ball']} to
{agg['maximum_retained_fraction_ball']} of that absolute mass remains in the
signed block sum.  Equivalently, the certified cancellation fraction is
{agg['minimum_cancellation_fraction_ball']} to
{agg['maximum_cancellation_fraction_ball']}.  A proof that takes absolute
values block by block therefore destroys the mechanism the data says matters.

The first genuinely unbounded obligation is

```text
E_hyb_main(t)=sum_k e_k(t)=Q_src(t)-T_upper(t).
```

It needs a constant-bearing, height-uniform signed estimate using the published
endpoint convention.  The existing corrected internal model is useful
diagnostic machinery, but it is not this source-aligned outer error and cannot
be substituted for it.

Proof boundary: rigorous finite decomposition for one fifteen-point window
near `t=10^10`.  It proves no height-uniform hybrid estimate, no cofinal
Jensen/PF theorem, no `Lambda<=0`, no RH, and no prize-level conclusion.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    ctx.dps = PRECISION
    ctx.threads = 1

    calibration = load_json(CALIBRATION)
    component = load_json(COMPONENT_SPLIT)
    classical = load_json(CLASSICAL_SPLIT)
    partition = load_json(BLOCK_PARTITION)
    block_zero = load_json(BLOCK_ZERO_IDENTITY)
    outer = load_json(OUTER_LEDGER)
    corrected = load_json(CORRECTED_INTERNAL)
    q_cut = load_json(SOURCE_Q_CUT)

    require(outer["paper_audit"]["paper_declares_hybrid_not_exact"] is True, "paper non-exact audit drift")
    require(outer["paper_audit"]["missing_token_count"] == 0, "paper token audit drift")
    require(outer["paper_audit"]["equations"] == [126, 127], "paper equation audit drift")
    require(partition["aggregate"]["cutoff_roster"][0] == 37946, "source endpoint cutoff drift")
    require(partition["aggregate"]["partition_total_terms"] == 39273, "partition size drift")
    require(block_zero["aggregate"]["transition_zero_count"] == EXPECTED_OUTPUTS, "block-zero transition audit drift")
    require(block_zero["decision"]["diagnostic_midpoint_is_a_published_cutoff_rule"] is False, "midpoint attribution drift")
    source_audit = audit_source()

    artifacts = (calibration, component, classical, partition)
    require(all(len(artifact["output_rows"]) == EXPECTED_OUTPUTS for artifact in artifacts), "output count drift")

    output_rows: list[dict[str, Any]] = []
    totals: list[arb] = []
    main_sums: list[arb] = []
    triangles: list[arb] = []
    retained_fractions: list[arb] = []
    cancellation_fractions: list[arb] = []
    subordinate_triangles: list[arb] = []

    for index, rows in enumerate(zip(*(artifact["output_rows"] for artifact in artifacts)), start=1):
        cal_row, component_row, classical_row, partition_row = rows
        require(all(int(row["output_index"]) == index for row in rows), f"output order drift at {index}")
        targets = [arb(row["target_t"]) for row in rows]
        require(all(targets[0].overlaps(target) for target in targets[1:]), f"target height drift at {index}")

        total = arb(cal_row["signed_source_minus_hardy_ball"])
        direct_total = arb(component_row["direct_source_minus_hardy_ball"])
        lower = arb(component_row["lower_source_minus_exact_ball"])
        upper = arb(component_row["hybrid_upper_source_minus_exact_ball"])
        addition = arb(component_row["final_addition_rounding_ball"])
        main_gap = arb(classical_row["hybrid_minus_classical_main_ball"])
        correction = arb(classical_row["exact_remaining_correction_ball"])
        classical_upper = arb(classical_row["hybrid_minus_exact_upper_ball"])
        block_rows = partition_row["block_rows"]
        require(len(block_rows) == EXPECTED_BLOCKS, f"block count drift at {index}")
        increments = [arb(row["increment_residual_ball"]) for row in block_rows]
        block_sum = sum(increments, arb(0))
        triangle = sum((abs(value) for value in increments), arb(0))
        block_zero_value = increments[0]
        later_net = sum(increments[1:], arb(0))
        retained = abs(block_sum) / triangle
        cancellation = 1 - retained
        subordinate = abs(lower) + abs(correction) + abs(addition)

        component_identity = lower + upper + addition - total
        upper_identity = main_gap - correction - upper
        block_identity = block_sum - main_gap
        total_identity = lower + block_sum - correction + addition - total
        require(total.overlaps(direct_total), f"direct total mismatch at {index}")
        require(upper.overlaps(classical_upper), f"upper split mismatch at {index}")
        require(component_identity.contains(0), f"component identity failed at {index}")
        require(upper_identity.contains(0), f"upper identity failed at {index}")
        require(block_identity.contains(0), f"block identity failed at {index}")
        require(total_identity.contains(0), f"total source-aligned identity failed at {index}")
        require(arb(block_rows[-1]["cumulative_residual_ball"]).overlaps(block_sum), f"telescoping failed at {index}")
        require(total < 0 and main_gap < 0 and block_zero_value < 0 and later_net > 0, f"signed pattern drift at {index}")
        require(triangle > TOLERANCE, f"block triangle no longer rejects tolerance at {index}")
        require(retained > 0 and retained < 1 and cancellation > 0 and cancellation < 1, f"cancellation ratio drift at {index}")

        totals.append(abs(total))
        main_sums.append(abs(block_sum))
        triangles.append(triangle)
        retained_fractions.append(retained)
        cancellation_fractions.append(cancellation)
        subordinate_triangles.append(subordinate)
        output_rows.append(
            {
                "output_index": index,
                "output_label": cal_row["output_label"],
                "target_t": cal_row["target_t"],
                "total_source_minus_exact_ball": total.str(PRECISION, more=True),
                "lower_source_minus_exact_ball": lower.str(PRECISION, more=True),
                "hybrid_upper_source_minus_exact_ball": upper.str(PRECISION, more=True),
                "source_aligned_upper_block_sum_ball": block_sum.str(PRECISION, more=True),
                "exact_upper_remaining_correction_ball": correction.str(PRECISION, more=True),
                "final_addition_rounding_ball": addition.str(PRECISION, more=True),
                "block_zero_increment_residual_ball": block_zero_value.str(PRECISION, more=True),
                "later_blocks_net_increment_residual_ball": later_net.str(PRECISION, more=True),
                "independent_block_triangle_ball": triangle.str(PRECISION, more=True),
                "retained_fraction_ball": retained.str(PRECISION, more=True),
                "cancellation_fraction_ball": cancellation.str(PRECISION, more=True),
                "subordinate_channel_triangle_ball": subordinate.str(PRECISION, more=True),
                "component_identity_residual_ball": component_identity.str(PRECISION, more=True),
                "upper_identity_residual_ball": upper_identity.str(PRECISION, more=True),
                "block_telescoping_identity_residual_ball": block_identity.str(PRECISION, more=True),
                "source_aligned_total_identity_residual_ball": total_identity.str(PRECISION, more=True),
                "sign_pattern": {
                    "total": sign_name(total),
                    "source_aligned_upper_block_sum": sign_name(block_sum),
                    "block_zero": sign_name(block_zero_value),
                    "later_blocks_net": sign_name(later_net),
                },
                "block_count": len(block_rows),
            }
        )

    aggregate = {
        "output_count": EXPECTED_OUTPUTS,
        "block_count": EXPECTED_BLOCKS,
        "block_row_count": EXPECTED_BLOCK_ROWS,
        "source_audit_count": len(source_audit),
        "source_audit_missing_token_count": sum(row["missing_token_count"] for row in source_audit),
        "source_endpoint_first_cutoff": partition["aggregate"]["cutoff_roster"][0],
        "source_endpoint_last_cutoff": partition["aggregate"]["cutoff_roster"][-1],
        "classical_upper_partition_term_count": partition["aggregate"]["partition_total_terms"],
        "total_identity_contains_zero_count": sum(arb(row["source_aligned_total_identity_residual_ball"]).contains(0) for row in output_rows),
        "block_telescoping_identity_contains_zero_count": sum(arb(row["block_telescoping_identity_residual_ball"]).contains(0) for row in output_rows),
        "strict_negative_total_count": sum(row["sign_pattern"]["total"] == "negative" for row in output_rows),
        "strict_negative_source_aligned_upper_block_sum_count": sum(row["sign_pattern"]["source_aligned_upper_block_sum"] == "negative" for row in output_rows),
        "strict_negative_block_zero_count": sum(row["sign_pattern"]["block_zero"] == "negative" for row in output_rows),
        "strict_positive_later_net_count": sum(row["sign_pattern"]["later_blocks_net"] == "positive" for row in output_rows),
        "block_triangle_rejects_0p005_count": sum(arb(row["independent_block_triangle_ball"]) > TOLERANCE for row in output_rows),
        "minimum_total_absolute_error_ball": minimum(totals).str(PRECISION, more=True),
        "maximum_total_absolute_error_ball": maximum(totals).str(PRECISION, more=True),
        "minimum_source_aligned_upper_block_sum_absolute_ball": minimum(main_sums).str(PRECISION, more=True),
        "maximum_source_aligned_upper_block_sum_absolute_ball": maximum(main_sums).str(PRECISION, more=True),
        "minimum_independent_block_triangle_ball": minimum(triangles).str(PRECISION, more=True),
        "maximum_independent_block_triangle_ball": maximum(triangles).str(PRECISION, more=True),
        "minimum_retained_fraction_ball": minimum(retained_fractions).str(PRECISION, more=True),
        "maximum_retained_fraction_ball": maximum(retained_fractions).str(PRECISION, more=True),
        "minimum_cancellation_fraction_ball": minimum(cancellation_fractions).str(PRECISION, more=True),
        "maximum_cancellation_fraction_ball": maximum(cancellation_fractions).str(PRECISION, more=True),
        "maximum_subordinate_channel_triangle_ball": maximum(subordinate_triangles).str(PRECISION, more=True),
    }
    require(aggregate["total_identity_contains_zero_count"] == EXPECTED_OUTPUTS, "total identity aggregate drift")
    require(aggregate["block_telescoping_identity_contains_zero_count"] == EXPECTED_OUTPUTS, "block identity aggregate drift")
    require(aggregate["block_triangle_rejects_0p005_count"] == EXPECTED_OUTPUTS, "triangle aggregate drift")

    dependencies = {
        path.stem: {"path": relative(path), "sha256": file_hash(path)}
        for path in EXPECTED_KINDS
    }
    dependencies.update(
        {
            "paper": {"path": relative(PAPER), "sha256": file_hash(PAPER)},
            "source_file": {"path": relative(SOURCE), "sha256": file_hash(SOURCE)},
        }
    )

    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_source_aligned_hybrid_error_ledger_gate",
        "status": "finite_source_aligned_equations_124_127_error_identity_certified_height_uniform_signed_block_sum_open",
        "passed": True,
        "scope": {
            "height_center": "1e10",
            "outputs": EXPECTED_OUTPUTS,
            "blocks": "0-35",
            "precision_decimal_digits": PRECISION,
            "cutoff_convention": "published_and_source_endpoint",
            "diagnostic_midpoint_used": False,
        },
        "source_aligned_identity": {
            "definitions": {
                "E_total": "Z_source-Z_exact",
                "E_lower": "L_source-L_exact",
                "E_hyb_main": "Q_source-T_upper=sum_(k=0)^35 e_k",
                "R_upper": "U_exact-T_upper",
                "rho_add": "Z_source-(L_source+Q_source)",
                "e_k": "(Q_k-Q_(k-1))-(T_k-T_(k-1))",
            },
            "identity": "E_total=E_lower+sum_(k=0)^35 e_k-R_upper+rho_add",
            "cutoff_coordinate": "alpha(N;t)=2N+t/(pi*N)",
            "source_tail_start_rule": "N_k=AINT(lower_root(alpha_k,t))+1",
            "first_tail_start_n_at_center": 37946,
        },
        "paper_audit": outer["paper_audit"],
        "source_audit": source_audit,
        "channel_ledger": channel_ledger(corrected, q_cut),
        "output_rows": output_rows,
        "aggregate": aggregate,
        "decision": {
            "finite_source_aligned_total_identity_certified": True,
            "source_endpoint_is_primary": True,
            "diagnostic_midpoint_subtraction_used": False,
            "blockwise_absolute_values_destroy_observed_cancellation": True,
            "corrected_internal_model_closes_source_aligned_outer_error": False,
            "first_unbounded_signed_combination_id": FIRST_OPEN_ID,
            "height_uniform_source_aligned_error_bound_proved": False,
            "rh_implication": False,
        },
        "next_obligation": "Derive a constant-bearing height-uniform signed estimate for E_hyb_main(t)=sum_k e_k(t) with the published/source endpoint cutoff, preserving cross-block cancellation; then promote uniform bounds for E_lower, R_upper, and rho_add and assemble the full source error.",
        "proof_boundary": "Rigorous finite source-aligned decomposition for one fifteen-point window near t=10^10. It proves no height-uniform equations-(126)--(127) remainder, no cofinal Jensen/PF theorem, no Lambda<=0, no RH, and no prize-level conclusion.",
        "dependencies": dependencies,
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "flint_threads": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
    }
    note = render_note(artifact)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(note, encoding="utf-8")
    print(
        "built source-aligned hybrid error ledger: "
        f"outputs={EXPECTED_OUTPUTS}, blocks={EXPECTED_BLOCK_ROWS}, open-first={FIRST_OPEN_ID}"
    )


if __name__ == "__main__":
    main()
