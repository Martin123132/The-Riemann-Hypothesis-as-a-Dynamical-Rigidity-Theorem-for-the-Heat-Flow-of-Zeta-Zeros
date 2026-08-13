#!/usr/bin/env python3
"""Compose the finite source-to-corrected-model bridge for Hardy block 20."""

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

import jensen_window_pf_newman_c1_hardy_block20_coefficient_cell_joint_recurrence_identity_gate as joint
import jensen_window_pf_newman_c1_hardy_block20_corrected_recurrence_defect_gate as corrected
import jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_full_gate as endpoint
import jensen_window_pf_newman_c1_hardy_block20_exact_eq69_saddle_family_gate as eq69
import jensen_window_pf_newman_c1_hardy_block20_special_function_residual_gate as special
import jensen_window_pf_newman_c1_hardy_block20_t2_lplus1_mismatch_gate as t2
import jensen_window_pf_newman_c1_hardy_block20_w2_w5_correspondence_gate as w2w5
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import acb, arb, ctx


EQ69 = eq69.RESULT
T2 = t2.RESULT
SPECIAL = special.RESULT
W2W5 = w2w5.RESULT
ENDPOINT = endpoint.RESULT
JOINT = joint.RESULT
CORRECTED_DEFECT = corrected.RESULT
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_source_to_corrected_model_bridge_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_source_to_corrected_model_bridge_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_block20_source_to_corrected_model_bridge_gate.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def upper_abs(value: arb | acb) -> Fraction:
    return cells.bound_fraction(abs(value).upper())


def lower_abs(value: acb) -> Fraction:
    return cells.bound_fraction(abs(value).lower())


def complex_record(value: acb) -> dict[str, Any]:
    return point_balls.acb_record(value, 55)


def algebra_audit() -> int:
    checks = 0
    fixtures = (
        (1, 2, 3, 5, 7, 11, 13, 17),
        (-2, 3, -5, 7, -11, 13, -17, 19),
        (Fraction(1, 3), Fraction(-2, 5), Fraction(3, 7), Fraction(5, 11), Fraction(-7, 13), Fraction(11, 17), Fraction(13, 19), Fraction(-17, 23)),
        (5, -4, 3, -2, 1, -6, 7, -8),
    )
    for ps, integral, model, t5_value, q_drop, q_extra, delta_parent, delta_integral in fixtures:
        ps, integral, model, t5_value = map(Fraction, (ps, integral, model, t5_value))
        q_drop, q_extra, delta_parent, delta_integral = map(
            Fraction, (q_drop, q_extra, delta_parent, delta_integral)
        )
        q_paper = q_drop + Fraction(5, 29)
        delta_q = q_paper - q_drop
        source_q = t5_value + q_drop + q_extra
        source_defect = ps - model - source_q
        delta_w1 = integral - model - t5_value
        after_w1 = source_defect - delta_w1
        after_deletion = after_w1 + q_extra
        final_transport = after_deletion + delta_parent - delta_integral - delta_q
        direct = ps + delta_parent - (integral + delta_integral) - q_paper
        require(after_w1 == ps - integral - q_drop - q_extra, "exact-W1 bridge algebra drift")
        require(after_deletion == ps - integral - q_drop, "t2 deletion bridge algebra drift")
        require(final_transport == direct, "source-to-corrected bridge algebra drift")
        checks += 3
    return checks


def load_rows(path: Path) -> tuple[dict[str, Any], dict[int, dict[str, Any]]]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    return artifact, {int(row["chain"]): row for row in artifact["rows"]}


def build() -> dict[str, Any]:
    for path in (EQ69, T2, SPECIAL, W2W5, ENDPOINT, JOINT, CORRECTED_DEFECT, CHECKER):
        require(path.is_file(), f"missing bridge dependency: {path}")
    ctx.dps = 180
    ctx.threads = 1
    exact_checks = algebra_audit()
    eq_artifact, eq_rows = load_rows(EQ69)
    t2_artifact, t2_rows = load_rows(T2)
    special_artifact, special_rows = load_rows(SPECIAL)
    w2_artifact, w2_rows = load_rows(W2W5)
    endpoint_artifact, endpoint_rows = load_rows(ENDPOINT)
    joint_artifact = json.loads(JOINT.read_text(encoding="utf-8"))
    corrected_artifact = json.loads(CORRECTED_DEFECT.read_text(encoding="utf-8"))
    roster = set(eq_rows)
    require(
        roster == set(t2_rows) == set(special_rows) == set(w2_rows) == set(endpoint_rows)
        and len(roster) == 374,
        "source-to-corrected bridge roster drift",
    )

    rows: list[dict[str, Any]] = []
    for chain in sorted(roster):
        eq_row = eq_rows[chain]
        t2_row = t2_rows[chain]
        special_row = special_rows[chain]
        w2_row = w2_rows[chain]
        endpoint_row = endpoint_rows[chain]
        exact_w1 = corrected.complex_from_record(eq_row["exact_w1_replacement"])
        post_w1 = corrected.complex_from_record(eq_row["residual_after_exact_w1"])
        t2_prior = corrected.complex_from_record(t2_row["prior_exact_w1_residual"])
        deletion = corrected.complex_from_record(t2_row["deletion_correction"])
        post_deletion = corrected.complex_from_record(t2_row["residual_after_deletion"])
        w2_prior = corrected.complex_from_record(w2_row["prior_post_lplus1_residual"])
        special_delta = corrected.complex_from_record(special_row["correlated_replacement"]["total_q_delta"])
        nonsaddle_transport = corrected.complex_from_record(w2_row["complete_nonsaddle_transport"])
        normalization_remainder = nonsaddle_transport - special_delta
        paper_residual = corrected.complex_from_record(w2_row["paper_residual"])
        endpoint_correction = corrected.complex_from_record(endpoint_row["exact_correction"])

        require(post_w1.overlaps(t2_prior), f"chain {chain} exact-W1/t2 handoff drift")
        require(post_deletion.overlaps(w2_prior), f"chain {chain} t2/W2-W5 handoff drift")
        require(paper_residual.overlaps(endpoint_correction), f"chain {chain} endpoint correction bridge drift")
        require(bool(eq_row["precision_overlap"]), f"chain {chain} equation-(69) precision drift")
        require(bool(t2_row["precision_overlap"]), f"chain {chain} t2 precision drift")
        require(bool(w2_row["precision_overlap"]), f"chain {chain} W2-W5 precision drift")
        require(
            Fraction(w2_row["residual_two_construction_gap_upper"]) < Fraction(1, 10**70),
            f"chain {chain} W2-W5 bridge gap drift",
        )

        rows.append(
            {
                "chain": chain,
                "sum_index": int(eq_row["sum_index"]),
                "branch": int(eq_row["branch"]),
                "child_length": int(eq_row["child_length"]),
                "exact_w1_replacement": complex_record(exact_w1),
                "exact_w1_replacement_magnitude_upper": decimal(upper_abs(exact_w1)),
                "post_exact_w1_residual": complex_record(post_w1),
                "t2_deletion_correction": complex_record(deletion),
                "t2_deletion_correction_magnitude_upper": decimal(upper_abs(deletion)),
                "post_t2_deletion_residual": complex_record(post_deletion),
                "special_function_q_replacement": complex_record(special_delta),
                "special_function_q_replacement_magnitude_upper": decimal(upper_abs(special_delta)),
                "complete_nonsaddle_transport": complex_record(nonsaddle_transport),
                "complete_nonsaddle_transport_magnitude_upper": decimal(upper_abs(nonsaddle_transport)),
                "post_special_normalization_remainder": complex_record(normalization_remainder),
                "post_special_normalization_remainder_magnitude_upper": decimal(upper_abs(normalization_remainder)),
                "parent_source_to_pi_transport_upper": w2_row["parent_source_to_pi_transport_upper"],
                "eq69_source_to_pi_transport_upper": w2_row["eq69_source_to_pi_transport_upper"],
                "final_paper_residual": complex_record(paper_residual),
                "final_paper_residual_magnitude_lower": decimal(lower_abs(paper_residual)),
                "final_paper_residual_magnitude_upper": decimal(upper_abs(paper_residual)),
                "endpoint_decomposition_overlap": True,
                "precision_and_identity_handoffs": True,
            }
        )

    def maximum(field: str) -> tuple[Fraction, int]:
        return max((Fraction(row[field]), int(row["chain"])) for row in rows)

    def minimum(field: str) -> tuple[Fraction, int]:
        return min((Fraction(row[field]), int(row["chain"])) for row in rows)

    max_w1 = maximum("exact_w1_replacement_magnitude_upper")
    max_t2 = maximum("t2_deletion_correction_magnitude_upper")
    max_special = maximum("special_function_q_replacement_magnitude_upper")
    max_nonsaddle = maximum("complete_nonsaddle_transport_magnitude_upper")
    max_normalization = maximum("post_special_normalization_remainder_magnitude_upper")
    max_parent_pi = maximum("parent_source_to_pi_transport_upper")
    max_eq69_pi = maximum("eq69_source_to_pi_transport_upper")
    min_final = minimum("final_paper_residual_magnitude_lower")
    max_final = maximum("final_paper_residual_magnitude_upper")

    outputs = [
        {
            "output_index": int(row["output_index"]),
            "output_label": row["output_label"],
            "corrected_model_majorant_upper": row["joint_corrected_recurrence_majorant_upper"],
            "within_requested_scale": bool(row["within_requested_scale"]),
        }
        for row in joint_artifact["transported_outputs"]
    ]
    require(len(outputs) == 15 and all(row["within_requested_scale"] for row in outputs), "bridge output budget drift")
    max_output = max((Fraction(row["corrected_model_majorant_upper"]), int(row["output_index"])) for row in outputs)

    aggregate = {
        "exact_algebra_checks": exact_checks,
        "call_count": len(rows),
        "identity_handoff_count": sum(bool(row["precision_and_identity_handoffs"]) for row in rows),
        "endpoint_decomposition_overlap_count": sum(bool(row["endpoint_decomposition_overlap"]) for row in rows),
        "maximum_exact_w1_replacement_magnitude_upper": decimal(max_w1[0]),
        "maximum_exact_w1_replacement_witness": max_w1[1],
        "maximum_t2_deletion_correction_magnitude_upper": decimal(max_t2[0]),
        "maximum_t2_deletion_correction_witness": max_t2[1],
        "maximum_special_function_q_replacement_magnitude_upper": decimal(max_special[0]),
        "maximum_special_function_q_replacement_witness": max_special[1],
        "maximum_complete_nonsaddle_transport_magnitude_upper": decimal(max_nonsaddle[0]),
        "maximum_complete_nonsaddle_transport_witness": max_nonsaddle[1],
        "maximum_post_special_normalization_remainder_magnitude_upper": decimal(max_normalization[0]),
        "maximum_post_special_normalization_remainder_witness": max_normalization[1],
        "maximum_parent_source_to_pi_transport_upper": decimal(max_parent_pi[0]),
        "maximum_parent_source_to_pi_transport_witness": max_parent_pi[1],
        "maximum_eq69_source_to_pi_transport_upper": decimal(max_eq69_pi[0]),
        "maximum_eq69_source_to_pi_transport_witness": max_eq69_pi[1],
        "minimum_final_paper_residual_magnitude_lower": decimal(min_final[0]),
        "minimum_final_paper_residual_witness": min_final[1],
        "maximum_final_paper_residual_magnitude_upper": decimal(max_final[0]),
        "maximum_final_paper_residual_witness": max_final[1],
        "maximum_source_recurrence_reconstruction_roundoff_upper": corrected_artifact["aggregate"]["maximum_source_roundoff_gap_abs_upper"],
        "maximum_corrected_model_output_majorant_upper": decimal(max_output[0]),
        "maximum_corrected_model_output_majorant_witness": max_output[1],
        "outputs_within_requested_scale": sum(bool(row["within_requested_scale"]) for row in outputs),
    }
    return {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_source_to_corrected_model_bridge_gate",
        "status": "finite_source_to_corrected_model_bridge_closes_at_374_saved_points_and_15_outputs",
        "scope": "The signed correction map from the accepted source recurrence to the mathematical corrected model at the 374 saved recursive block-20 calls.",
        "bridge_identity": {
            "source_defect": "D_src=P_src-M*C-q_src",
            "exact_w1_step": "D_1=D_src-(I69_src-M*C-t5_src)",
            "t2_deletion_step": "D_2=D_1+q_extra=P_src-I69_src-Q_src_drop",
            "normalization_step": "D_corr=D_2+(P_pi-P_src)-(I69_pi-I69_src)-(Q_paper-Q_src_drop)",
            "direct_target": "D_corr=P_pi-I69_pi-Q_paper=Q_exact-Q_paper",
            "special_split": "Q_paper-Q_src_drop=delta_q_special+delta_q_post_special_normalization",
        },
        "rows": rows,
        "transported_outputs": outputs,
        "aggregate": aggregate,
        "dependencies": {
            "exact_equation_69": {"path": relative(EQ69), "sha256": file_hash(EQ69)},
            "t2_lplus1": {"path": relative(T2), "sha256": file_hash(T2)},
            "special_functions": {"path": relative(SPECIAL), "sha256": file_hash(SPECIAL)},
            "w2_w5_correspondence": {"path": relative(W2W5), "sha256": file_hash(W2W5)},
            "endpoint_decomposition": {"path": relative(ENDPOINT), "sha256": file_hash(ENDPOINT)},
            "joint_recurrence": {"path": relative(JOINT), "sha256": file_hash(JOINT)},
            "source_roundoff": {"path": relative(CORRECTED_DEFECT), "sha256": file_hash(CORRECTED_DEFECT)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": (
            "Implement and source-equivalence-test the correction map, then transport its W1, t2, special-function, "
            "normalization, and floating-arithmetic pieces over physical coefficient-cell slices before scaling to other blocks."
        ),
        "proof_boundary": (
            "Rigorous composition of previously certified finite block-20 exact-point corrections and output endpoint "
            "majorants only. The accepted source is not modified and the correction map is not yet enclosed over whole "
            "physical cells, other blocks, or height. This does not close cross-block accumulation, the outer Hardy "
            "representation, Lambda<=0, PF-infinity, RH, or a prize-level conclusion."
        ),
    }


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 source-to-corrected-model bridge

Date: 2026-08-09
Status: {artifact['status']}; finite correction theorem only, not a proof of RH

## Signed bridge

The admitted source recurrence and the corrected mathematical model are joined
without adding unrelated absolute maxima.  With `Q_src_drop` denoting the
source nonsaddle correction after deleting the source-only `t2` term,

```text
D_src = P_src-M*C-q_src,
D_1   = D_src-(I69_src-M*C-t5_src),
D_2   = D_1+q_extra = P_src-I69_src-Q_src_drop,
D_corr= D_2+(P_pi-P_src)-(I69_pi-I69_src)
              -(Q_paper-Q_src_drop)
      = P_pi-I69_pi-Q_paper
      = Q_exact-Q_paper.
```

All 374 exact-W1 to `t2`, `t2` to W2--W5, and final endpoint-decomposition
handoffs overlap.  The complete nonsaddle transport also splits exactly into
the independently measured PSI/ERF replacement and a post-special-function
normalization remainder.

```text
maximum exact-W1 replacement             {aggregate['maximum_exact_w1_replacement_magnitude_upper']}
maximum t2 deletion correction           {aggregate['maximum_t2_deletion_correction_magnitude_upper']}
maximum PSI/ERF q replacement             {aggregate['maximum_special_function_q_replacement_magnitude_upper']}
maximum complete nonsaddle transport      {aggregate['maximum_complete_nonsaddle_transport_magnitude_upper']}
maximum post-special normalization        {aggregate['maximum_post_special_normalization_remainder_magnitude_upper']}
maximum parent source-to-2*pi transport   {aggregate['maximum_parent_source_to_pi_transport_upper']}
maximum I69 source-to-2*pi transport      {aggregate['maximum_eq69_source_to_pi_transport_upper']}
maximum source reconstruction roundoff    {aggregate['maximum_source_recurrence_reconstruction_roundoff_upper']}
maximum corrected output majorant         {aggregate['maximum_corrected_model_output_majorant_upper']}
outputs below 0.005                       {aggregate['outputs_within_requested_scale']} / 15
```

## Boundary

This proves the signed correction map at the saved finite inputs and connects
its final residual to the independently certified endpoint correction.  It
does not alter the accepted source, enclose the correction map over entire
physical cells, treat other blocks or height, control the outer Hardy
representation, or prove `Lambda<=0`, PF-infinity, RH, or a prize theorem.
"""


def main() -> int:
    artifact = build()
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built source-to-corrected-model bridge: "
        f"{artifact['aggregate']['call_count']} calls, "
        f"{artifact['aggregate']['outputs_within_requested_scale']}/15 outputs below 0.005"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
