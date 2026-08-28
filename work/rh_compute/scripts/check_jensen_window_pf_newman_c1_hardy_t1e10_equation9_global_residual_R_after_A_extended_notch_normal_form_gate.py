#!/usr/bin/env python3
"""Independently check the post-A extended-notch normal form gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_normal_form_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")

A = 159_577
T_LO = 622
T_HI = 39_894
W_LO = 39_853
W_HI = 39_936


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True and artifact.get("kind") == STEM, "artifact identity drift")
    decisions = artifact["decision"]
    for key in (
        "post_A_positive_bulk_mask_is_one_contiguous_extended_notch",
        "six_row_mask_compressed_coefficientwise",
        "post_A_endpoint_roster_reduced_to_two_rows",
        "extended_notch_one_cell_fold_and_modular_transform_inherited",
        "A_upper_projector_edge_moved_out_of_half_boundary_collar",
        "B_subtractions_and_negative_bulk_ownership_preserved",
    ):
        require(decisions.get(key) is True, f"missing decision: {key}")
    for key in (
        "extended_notch_alone_proves_joined_bound",
        "quantitative_R_after_A_bound_proved",
        "rh_implication",
    ):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    target = set(range(T_LO, T_HI + 1))
    window = set(range(W_LO, W_HI + 1))
    extended = set(range(T_LO, W_HI + 1))
    require(target | window == extended, "independent union failure")
    require((len(target), len(window), len(target & window), len(window - target), len(extended)) == (39_273, 84, 42, 42, 39_315), "independent count failure")
    for mode in range(1, W_HI + 2):
        chi_t = int(mode in target)
        chi_w = int(mode in window)
        chi_u = int(mode in extended)
        require(chi_t + chi_w * (1 - chi_t) == chi_u, f"independent mask failure at {mode}")

    source_all, source_u, endpoint_u, bulk_u, a_pair = sp.symbols("S S_U Q_U P_U A_pair")
    expression = source_all - source_u + endpoint_u - a_pair
    expression = expression.subs(endpoint_u, source_u - bulk_u)
    require(sp.expand(expression - (source_all - bulk_u - a_pair)) == 0, "independent folded identity failure")

    old_edge = Fraction(2 * T_HI + 1, 2)
    new_edge = Fraction(2 * W_HI + 1, 2)
    half_boundary = Fraction(A, 4)
    require(new_edge - old_edge == 42, "independent edge shift failure")
    require(old_edge - half_boundary == Fraction(1, 4), "independent old gap failure")
    require(new_edge - half_boundary == Fraction(169, 4), "independent new gap failure")
    q_half = lambda mode: (A - 4 * mode) / 2
    require(q_half(Fraction(W_LO)) == Fraction(165, 2), "independent first deleted q failure")
    require(q_half(Fraction(W_HI)) == Fraction(-167, 2), "independent last deleted q failure")
    require(q_half(Fraction(W_LO - 1)) == Fraction(169, 2), "independent lower survivor q failure")
    require(q_half(Fraction(W_HI + 1)) == Fraction(-171, 2), "independent upper survivor q failure")
    require(q_half(new_edge) == Fraction(-169, 2), "independent new-edge q failure")

    modular = artifact["modular_certificate"]
    require("d_A=39936.5" in modular["edges"]["upper"], "extended upper edge missing")
    require("(-1)^k" in modular["edge_phase_guard"], "half-integer phase guard missing")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")

    require(BUILDER.is_file() and NOTE.is_file(), "builder or note missing")
    note = NOTE.read_text(encoding="utf-8")
    require("contiguous extended notch" in note, "contiguous-notch statement missing")
    require("half-integer edge has `q_A=-84.5`" in note, "normal-gap statement missing")
    require("No B atom or negative Gamma bulk" in note, "ownership guard missing")
    require("No quantitative bound for the remaining joined kernel" in note, "proof boundary missing")
    print("independently checked post-A extended-notch normal form", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
