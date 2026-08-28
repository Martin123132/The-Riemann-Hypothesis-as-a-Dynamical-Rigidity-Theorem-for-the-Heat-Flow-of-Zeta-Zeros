#!/usr/bin/env python3
"""Independently check the A-corner carrier and projector orientation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_corner_projector_orientation_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
HEIGHT = 10_000_000_000
A = 159_577


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def main() -> int:
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True and artifact.get("kind") == STEM, "artifact identity drift")
    decisions = artifact["decision"]
    for key in (
        "odd_endpoint_phase_reduced_to_common_classical_carrier",
        "A_endpoint_sign_identified_as_epsilon_A_minus_one",
        "paired_endpoint_coefficient_continuous_across_target_edge",
        "positive_full_line_bulk_projector_jump_identified",
        "equation9_normalization_and_half_projection_restored",
        "two_mode_affine_A_tangent_wedge_projection_certified",
        "two_mode_affine_correction_is_target_scale_significant",
    ):
        require(decisions.get(key) is True, f"missing decision: {key}")
    for key in ("exact_A_corner_remainder_bound_proved", "complete_A_roster_sum_proved", "R_Dir_bound_proved", "rh_implication"):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    corner_path = REPO_ROOT / artifact["dependencies"]["corner"]["path"]
    corner = load_json(corner_path)
    corner_rows = {row["mode"]: row for row in corner["rows"]}
    saved = artifact["interval_certificate"]
    saved_rows = {row["mode"]: row for row in saved["rows"]}
    require(set(saved_rows) == {39_894, 39_895}, "mode roster drift")

    ctx.dps = 125
    ctx.threads = 1
    pi, t, endpoint = arb.pi(), arb(HEIGHT), arb(A)
    imaginary = acb(0, 1)
    normalizer = (pi / (32 * t)) ** (arb(1) / 4)
    rotation = (-imaginary * pi / 8).exp()
    fresh_sum = arb(0)
    fresh_correction_sum = arb(0)
    for mode_int in (39_894, 39_895):
        mode = arb(mode_int)
        source = corner_rows[mode_int]
        affine = parse_complex(source["affine_wedge_C_aff"])
        scalar = parse_complex(source["canonical_wedge_C"])
        K = pi * endpoint * mode
        r = t / (2 * pi * mode**2)
        W0 = acb(K, -1) * (4 * r ** arb("0.75") / (endpoint * (pi * t).sqrt()))

        x_m = 2 * pi * mode**2 / (t + 2 * pi * mode**2)
        psi_star = -pi * mode**2 / x_m + t * ((1 - x_m) / x_m).log() / 2
        unsimplified_phase = pi * mode * endpoint + psi_star
        fresh_carrier = (imaginary * unsimplified_phase).exp()
        raw_affine = -fresh_carrier * W0 * affine
        raw_correction = -fresh_carrier * W0 * (affine - scalar)
        fresh = 2 * normalizer * (rotation * raw_affine).real
        fresh_correction = 2 * normalizer * (rotation * raw_correction).real
        require(fresh.overlaps(arb(saved_rows[mode_int]["physical_affine_A_tangent_wedge_ball"])), f"physical affine value misses at {mode_int}")
        require(fresh_correction.overlaps(arb(saved_rows[mode_int]["physical_affine_correction_ball"])), f"physical affine correction misses at {mode_int}")
        require(saved_rows[mode_int]["paired_endpoint_projector_coefficient"] == 1, f"endpoint coefficient drift at {mode_int}")
        require(saved_rows[mode_int]["target_indicator"] == (1 if mode_int == 39_894 else 0), f"target indicator drift at {mode_int}")
        require(saved_rows[mode_int]["positive_full_line_bulk_projector_coefficient"] == (0 if mode_int == 39_894 else 1), f"bulk projector coefficient drift at {mode_int}")
        fresh_sum += fresh
        fresh_correction_sum += fresh_correction

    require(fresh_sum.overlaps(arb(saved["two_mode_physical_affine_A_tangent_wedge_ball"])), "two-mode affine sum misses")
    require(fresh_correction_sum.overlaps(arb(saved["two_mode_physical_affine_correction_ball"])), "two-mode correction sum misses")
    require(fresh_sum.upper() < arb("-0.0077"), "two-mode sign drift")
    require(fresh_correction_sum.upper() < arb("-6.3e-6"), "affine correction scale drift")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(NOTE.is_file(), "missing note")
    note = NOTE.read_text(encoding="utf-8")
    require("coefficient `+1` on both" in note, "projector orientation missing from note")
    require("No fitted constant" in note, "pi provenance missing from note")
    require("No\nexact A-corner remainder" in note, "proof boundary missing from note")
    print("independently checked affine A-corner carrier and projector orientation", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
