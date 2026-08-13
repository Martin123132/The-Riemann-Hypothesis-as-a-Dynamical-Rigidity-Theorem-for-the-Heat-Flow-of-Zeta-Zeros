#!/usr/bin/env python3
"""Independently check the saved block-zero source-formula identity gate."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

from flint import arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_block_zero_source_formula_identity_gate as gate
import jensen_window_pf_newman_c1_hardy_t1e10_equation124_classical_block_partition_gate as partition


CHECK_PRECISION = 130


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(gate.RESULT.is_file(), "missing block-zero identity result")
    require(gate.NOTE.is_file(), "missing block-zero identity note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_block_zero_source_formula_identity_gate", "kind drift")
    require(artifact["status"] == "finite_block_zero_equation62_approximation_error_localized", "status drift")
    require(bool(artifact["passed"]), "block-zero identity gate failed")
    require(artifact["scope"]["precisions_decimal_digits"] == list(gate.PRECISIONS), "precision roster drift")
    require(artifact["scope"]["alpha_terms_per_output"] == gate.ALPHA_COUNT, "alpha count drift")

    gate.audit_source()
    algebra = gate.exact_algebra_audit()
    require(bool(algebra["all_polynomial_remainders_zero"]), "fresh exact algebra audit failed")
    require(algebra == artifact["equation62_source_identity"], "saved exact algebra audit drift")

    transition = gate.central_parameters(CHECK_PRECISION)
    saved_transition = artifact["transition_certificate"]["high_precision"]
    for key, saved_key in (("a", "a_ball"), ("t6", "t6_ball"), ("g", "g_ball"), ("gm", "gm_ball"), ("yphase", "yphase_ball")):
        require(transition[key].overlaps(arb(saved_transition[saved_key])), f"transition ball drift: {key}")
    require(transition["g"] > arb("3.2") and transition["gm"] > arb("3.2"), "transition exclusion lost")

    parent = json.loads(partition.RESULT.read_text(encoding="utf-8"))
    source_block_zero = partition.load_stage_zero()[0]
    rows = artifact["output_rows"]
    require(len(rows) == gate.EXPECTED_OUTPUTS, "output roster drift")
    source_model_abs: list[arb] = []
    shift_effect_abs: list[arb] = []
    model_target_abs: list[arb] = []
    source_target_abs: list[arb] = []
    for output_index, row in enumerate(rows, start=1):
        offset = Decimal(output_index - 8) / Decimal(100)
        target_t = format(Decimal("1e10") + offset, "f")
        require(row["output_index"] == output_index, "output index drift")
        require(row["output_label"] == f"t{offset:+.2f}", "output label drift")
        require(row["target_t"] == target_t, "target height drift")
        require(row["alpha_count"] == gate.ALPHA_COUNT, "row alpha count drift")

        model = gate.ideal_block_zero(offset, CHECK_PRECISION, source_literal=True)
        intended = gate.ideal_block_zero(offset, CHECK_PRECISION, source_literal=False)
        saved_high = row["high_precision"]
        for key, saved_key in (
            ("a", "shifted_a_ball"),
            ("yphase", "shifted_yphase_ball"),
            ("raw_sum", "raw_ter_sum_ball"),
            ("normalization", "normalization_ball"),
            ("normalized_sum", "normalized_equation62_sum_ball"),
        ):
            require(model[key].overlaps(arb(saved_high[saved_key])), f"fresh model misses saved ball at output {output_index}: {key}")
        require(intended["normalized_sum"].overlaps(arb(saved_high["intended_exact_shift_normalized_sum_ball"])), f"fresh intended-shift model misses saved ball at output {output_index}")

        cutoff = int(parent["output_rows"][output_index - 1]["block_rows"][0]["cutoff_n"])
        require(row["classical_cutoff_n"] == cutoff, "classical cutoff drift")
        target = partition.cumulative_targets(target_t, [cutoff], CHECK_PRECISION)[cutoff]
        require(target.overlaps(arb(saved_high["exact_classical_target_ball"])), f"fresh target misses saved ball at output {output_index}")

        ctx.dps = CHECK_PRECISION
        source = arb(source_block_zero["zsum"][output_index - 1])
        source_minus_model = source - model["normalized_sum"]
        intended_minus_source_model = intended["normalized_sum"] - model["normalized_sum"]
        model_minus_target = model["normalized_sum"] - target
        source_minus_target = source - target
        decomposition = source_minus_model + model_minus_target - source_minus_target
        require(decomposition.contains(0), f"fresh decomposition excludes zero at output {output_index}")
        require(abs(source_minus_model) < gate.SOURCE_MODEL_GUARD, f"fresh source/model guard failed at output {output_index}")
        require(abs(intended_minus_source_model) < gate.SHIFT_LITERAL_GUARD, f"fresh default-real shift guard failed at output {output_index}")
        require(model_minus_target < -gate.REQUESTED_TOLERANCE, f"fresh model/target guard failed at output {output_index}")
        require(source_minus_model.overlaps(arb(row["source_minus_model_ball"])), "saved source/model residual drift")
        require(intended_minus_source_model.overlaps(arb(row["intended_minus_source_default_real_model_ball"])), "saved default-real shift residual drift")
        require(model_minus_target.overlaps(arb(row["model_minus_target_ball"])), "saved model/target residual drift")
        require(source_minus_target.overlaps(arb(row["source_minus_target_ball"])), "saved source/target residual drift")
        source_model_abs.append(abs(source_minus_model))
        shift_effect_abs.append(abs(intended_minus_source_model))
        model_target_abs.append(abs(model_minus_target))
        source_target_abs.append(abs(source_minus_target))

    aggregate = artifact["aggregate"]
    require(aggregate["output_count"] == gate.EXPECTED_OUTPUTS, "output aggregate drift")
    require(aggregate["alpha_terms_per_output"] == gate.ALPHA_COUNT, "alpha aggregate drift")
    require(aggregate["total_reconstructed_alpha_terms"] == gate.EXPECTED_OUTPUTS * gate.ALPHA_COUNT, "term aggregate drift")
    require(aggregate["transition_zero_count"] == gate.EXPECTED_OUTPUTS, "transition aggregate drift")
    require(aggregate["source_model_below_1e_minus_20_count"] == gate.EXPECTED_OUTPUTS, "source/model aggregate drift")
    require(aggregate["default_real_shift_effect_below_2e_minus_12_count"] == gate.EXPECTED_OUTPUTS, "default-real shift aggregate drift")
    require(aggregate["model_target_rejects_0p005_count"] == gate.EXPECTED_OUTPUTS, "model/target aggregate drift")
    require(max(source_model_abs, key=lambda x: x.upper()).overlaps(arb(aggregate["maximum_source_model_residual_absolute_ball"])), "maximum source/model aggregate drift")
    require(max(shift_effect_abs, key=lambda x: x.upper()).overlaps(arb(aggregate["maximum_default_real_shift_effect_absolute_ball"])), "maximum default-real shift aggregate drift")
    require(min(model_target_abs, key=lambda x: x.lower()).overlaps(arb(aggregate["minimum_model_target_residual_absolute_ball"])), "minimum model/target aggregate drift")
    require(max(model_target_abs, key=lambda x: x.upper()).overlaps(arb(aggregate["maximum_model_target_residual_absolute_ball"])), "maximum model/target aggregate drift")
    require(max(source_target_abs, key=lambda x: x.upper()).overlaps(arb(aggregate["maximum_source_target_residual_absolute_ball"])), "maximum source/target aggregate drift")

    context = artifact["published_remainder_context"]
    ctx.dps = CHECK_PRECISION
    equation57_bound = arb("6.15") * arb("1e10") ** (-arb(1) / 12)
    require(context["paper_equation"] == "57", "published remainder equation drift")
    require(equation57_bound.overlaps(arb(context["equation57_bound_at_1e10_ball"])), "published equation-57 bound drift")
    require((equation57_bound / gate.REQUESTED_TOLERANCE).overlaps(arb(context["bound_to_0p005_ratio_ball"])), "published bound ratio drift")
    require(context["strict_threshold_integer"] == str(1230**12), "published bound threshold drift")
    require(not bool(context["is_signed_block_zero_remainder"]), "global bound promoted to signed block remainder")
    require(not bool(context["closes_0p005_at_1e10"]), "global bound falsely closes saved tolerance")

    decision = artifact["decision"]
    for key in (
        "transition_creates_block_zero_defect",
        "amplitude_normalization_mismatch_creates_block_zero_defect",
        "phase_translation_mismatch_creates_block_zero_defect",
        "observable_source_ideal_departure_creates_block_zero_defect",
        "default_real_shift_literal_mismatch_creates_block_zero_defect",
    ):
        require(not bool(decision[key]), f"false cause promoted: {key}")
    require(bool(decision["source_aligned_equation62_hybrid_residual_rejects_0p005_at_saved_outputs"]), "source-aligned residual localization lost")
    require(not bool(decision["diagnostic_midpoint_is_a_published_cutoff_rule"]), "diagnostic midpoint promoted to published cutoff")
    require(not bool(decision["height_uniform_equation62_remainder_bound_proved"]), "height-uniform overpromotion")
    require(not bool(decision["finite_gate_has_rh_implication"]), "RH overpromotion")

    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file {path}")
            require(file_hash(path) == record["sha256"], f"hash drift {path}")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("not a proof of RH", "Source-attribution notice", "equation (62)", "below `1e-20`", "equation (57)", "not a signed", "no RH implication"):
        require(token in note, f"note boundary token missing: {token}")
    print(
        "validated block-zero source identity: "
        f"{aggregate['total_reconstructed_alpha_terms']} terms, transition-zero={aggregate['transition_zero_count']}, "
        f"source-model<1e-20={aggregate['source_model_below_1e_minus_20_count']}, "
        f"eq62-rejects-0.005={aggregate['model_target_rejects_0p005_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
