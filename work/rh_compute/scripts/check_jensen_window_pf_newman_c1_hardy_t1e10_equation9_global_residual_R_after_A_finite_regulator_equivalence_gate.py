#!/usr/bin/env python3
"""Independently check the finite-regulator R_after_A equivalence gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_finite_regulator_equivalence_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
B = 5_122_421
TARGET_START = 622
TARGET_END = 39_894
A_START = 39_853
A_END = 39_936
ATOM_ORDER = ("P_plus", "P_minus", "A_plus", "A_minus", "B_plus", "B_minus")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def vector(*values: int) -> dict[str, int]:
    return dict(zip(ATOM_ORDER, values, strict=True))


def independent_vectors(mode: int) -> tuple[dict[str, int], dict[str, int], dict[str, int]]:
    chi = int(TARGET_START <= mode <= TARGET_END)
    pre = vector(1 - chi, 1, 1, 1, 1, 1)
    notch = vector(0, 0, 0, 0, 0, 0)
    if A_START <= mode <= A_END:
        notch = vector(1 - chi, 0, 1, 1, 0, 0)
    post = {atom: pre[atom] - notch[atom] for atom in ATOM_ORDER}
    return pre, notch, post


def main() -> int:
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True and artifact.get("kind") == STEM, "artifact identity drift")
    decisions = artifact["decision"]
    for key in (
        "finite_R_after_A_defined_on_one_common_regulator",
        "common_kernel_equals_six_class_post_A_mask",
        "analytic_B_outer_cutoff_condition_M_ge_B_enforced",
        "A_endpoint_and_outside_target_positive_bulk_deleted_algebraically_before_norms",
        "B_window_and_outer_subtractions_identical_in_both_forms",
        "common_physical_limit_order_preserved",
        "zero_half_negative_endpoint_and_remote_positive_sectors_remain_joined",
    ):
        require(decisions.get(key) is True, f"missing decision: {key}")
    for key in (
        "quantitative_R_after_A_upper_bound_proved",
        "complete_R_Dir_proved",
        "complete_Q_K_minus_T_proved",
        "rh_implication",
    ):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    require(artifact["scope"]["finite_cutoff"] == f"integer M>=B={B}", "finite-cutoff condition drift")
    require(artifact["scope"]["limit_order"] == "M to infinity first at fixed epsilon; epsilon down to zero second", "limit-order drift")

    Pp, Pm, Ap, Am, Bp, Bm, chi = sp.symbols("Pp Pm Ap Am Bp Bm chi")
    pair = (Pp + Ap + Bp) + (Pm + Am + Bm) - chi * Pp
    mask = (1 - chi) * Pp + Pm + Ap + Am + Bp + Bm
    require(sp.expand(pair - mask) == 0, "independent pair identity failed")
    H, D, IT, PT, QT = sp.symbols("H D IT PT QT")
    require(sp.expand((H + D - IT + QT).subs(IT, PT + QT) - (H + D - PT)) == 0, "independent common-kernel identity failed")

    certificate = artifact["finite_equivalence_certificate"]
    rows = certificate["positive_class_rows"]
    expected_ids = [
        "low_positive_pairs",
        "ordinary_target_pairs_before_A_window",
        "A_window_target_pairs",
        "A_window_outer_pairs",
        "remote_positive_pairs",
    ]
    require([row["sector_id"] for row in rows] == expected_ids, "class roster drift")
    expected_representatives = (1, 622, 39_853, 39_895, 39_937)
    for row, mode in zip(rows, expected_representatives, strict=True):
        _, _, post = independent_vectors(mode)
        require(row["post_A_coefficients"] == post, f"class coefficient drift: {row['sector_id']}")

    stored_samples = {sample["mode"]: sample for sample in certificate["boundary_samples"]}
    expected_samples = {1, 621, 622, 39_852, 39_853, 39_894, 39_895, 39_936, 39_937, B - 1, B, B + 17}
    require(set(stored_samples) == expected_samples, "boundary-sample roster drift")
    for mode, sample in stored_samples.items():
        pre, notch, post = independent_vectors(mode)
        require(sample["pre_A_coefficients"] == pre, f"pre-A coefficient drift at {mode}")
        require(sample["A_notch_coefficients"] == notch, f"A-notch coefficient drift at {mode}")
        require(sample["post_A_coefficients"] == post, f"post-A coefficient drift at {mode}")

    for cutoff_record in certificate["count_checks"]:
        cutoff = cutoff_record["cutoff"]
        require(cutoff >= B, "count check below B cutoff")
        counts = [621, 39_231, 42, 42, cutoff - 39_936]
        require(cutoff_record["class_counts"] == counts, f"class counts drift at M={cutoff}")
        require(sum(counts) == cutoff_record["positive_mode_total"] == cutoff, f"positive total drift at M={cutoff}")
        require(cutoff_record["symmetric_mode_total_with_zero"] == 2 * cutoff + 1, f"symmetric total drift at M={cutoff}")

    require(certificate["zero_sector"]["summand"] == "H_x+I_0", "zero/half-current join drift")
    require(certificate["zero_sector"]["separate_norm_admissible"] is False, "zero-sector norm guard lost")
    require(certificate["B_subtractions_common_to_both_forms"] is True, "B common-subtraction guard lost")
    require(set(certificate["independent_norm_forbidden_sectors"]) == {
        "endpoint_half_current",
        "zero_mode",
        "negative_bulk",
        "remaining_endpoint_atoms",
        "remote_positive_pairs",
    }, "joined-sector roster drift")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(BUILDER.is_file() and NOTE.is_file(), "builder or note missing")
    note = NOTE.read_text(encoding="utf-8")
    require("M>=B=5122421" in note, "strong cutoff missing from note")
    require("M to infinity" in note and "epsilon down 0" in note, "limit order missing from note")
    require("No numerical or analytic upper bound" in note, "proof boundary missing from note")
    print("independently checked finite-regulator R_after_A equivalence", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
