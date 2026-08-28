#!/usr/bin/env python3
"""Independently check the global R_Dir ownership ledger."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_Dir_projector_ownership_ledger_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
ATOM_ORDER = ("P_plus", "P_minus", "A_plus", "A_minus", "B_plus", "B_minus")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def add(left: dict[str, int], right: dict[str, int]) -> dict[str, int]:
    return {atom: left[atom] + right[atom] for atom in ATOM_ORDER}


def main() -> int:
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True and artifact.get("kind") == STEM, "artifact identity drift")
    decisions = artifact["decision"]
    for key in (
        "all_finite_mode_classes_partitioned_without_gap_or_overlap",
        "A_endpoint_and_positive_bulk_coefficients_transferred_exactly",
        "B_window_A_transition_and_analytic_B_outer_no_overlap_guards_certified",
        "zero_half_negative_and_remote_positive_sectors_remain_joined",
        "post_A_one_sided_targets_derived",
    ):
        require(decisions.get(key) is True, f"missing decision: {key}")
    for key in ("joined_R_after_A_bound_proved", "complete_R_Dir_proved", "complete_Q_K_minus_T_proved", "rh_implication"):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    Pp, Pm, Ap, Am, Bp, Bm, chi = sp.symbols("Pp Pm Ap Am Bp Bm chi")
    pair = (Pp + Ap + Bp) + (Pm + Am + Bm) - chi * Pp
    projected = (1 - chi) * Pp + Pm + Ap + Am + Bp + Bm
    require(sp.expand(pair - projected) == 0, "independent paired-projector check failed")
    M = sp.symbols("M", integer=True, positive=True)
    require(sp.expand(621 + 39_231 + 42 + 42 + (M - 39_936) - M) == 0, "independent partition count failed")

    rows = artifact["mode_ownership_rows"]
    require([row["sector_id"] for row in rows] == [
        "zero_mode",
        "low_positive_pairs",
        "ordinary_target_pairs_before_A_window",
        "A_window_target_pairs",
        "A_window_outer_pairs",
        "remote_positive_pairs",
    ], "mode-row order or identity drift")
    for row in rows[1:]:
        require(add(row["A_transition_coefficients"], row["post_A_coefficients"]) == row["pre_A_coefficients"], f"row split failed: {row['sector_id']}")
        require(row["separate_norm_admissible"] is False, f"separate-norm guard lost: {row['sector_id']}")
    require(rows[3]["A_transition_coefficients"] == dict(zip(ATOM_ORDER, (0, 0, 1, 1, 0, 0), strict=True)), "A target transfer drift")
    require(rows[4]["A_transition_coefficients"] == dict(zip(ATOM_ORDER, (1, 0, 1, 1, 0, 0), strict=True)), "A outer transfer drift")

    ctx.dps = 120
    ctx.threads = 1
    A_dep_record = artifact["dependencies"]["A_transition"]
    A_dep = load_json(REPO_ROOT / A_dep_record["path"])
    A_ball = arb(A_dep["certificate"]["components"]["projector_completed_exact_A_transition_ball"])
    cert = artifact["extraction_certificate"]
    require(A_ball.overlaps(arb(cert["A_transition_ball"])), "A transition interval drift")
    work_room = arb("0.00014057919999999995") - A_ball
    negative_room = arb("0.00013197919999999995") - A_ball
    require(arb(cert["working_R_after_A_threshold_ball"]).overlaps(work_room), "working threshold misses")
    require(arb(cert["negative_R_after_A_threshold_ball"]).overlaps(negative_room), "negative threshold misses")
    require(arb(cert["working_R_after_A_safe_upper"]) <= arb(work_room.lower()), "unsafe working threshold")
    require(arb(cert["negative_R_after_A_safe_upper"]) <= arb(negative_room.lower()), "unsafe negative threshold")
    require(arb(cert["working_R_after_A_safe_upper"]) > arb("0.0368"), "working room scale drift")

    extractions = {row["sector_id"]: row for row in cert["extractions"]}
    require(set(extractions) == {"B_trace_window", "analytic_B_outer", "projector_completed_A_transition", "joined_R_after_A"}, "extraction roster drift")
    require(extractions["B_trace_window"]["owner"] == "outside_R_Dir", "B window ownership drift")
    require(extractions["analytic_B_outer"]["owner"] == "outside_R_Dir", "B outer ownership drift")
    require(extractions["projector_completed_A_transition"]["owner"] == "inside_R_Dir_extracted_signed", "A ownership drift")
    require(extractions["joined_R_after_A"]["bound_type"] == "open_one_sided_join", "joined-bound status drift")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(NOTE.is_file(), "missing note")
    note = NOTE.read_text(encoding="utf-8")
    require("zero mode, negative full-line bulk" in note, "joined-sector guard missing")
    require("No bound for" in note and "R_after_A" in note, "proof boundary missing")
    print("independently checked global R_Dir ownership ledger", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
