#!/usr/bin/env python3
"""Independently check corrected-selected source projection orientation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
from flint import arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_corrected_selected_source_projection_orientation_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
C = 159_577


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    require(RESULT.is_file(), "missing result artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "saved gate did not pass")
    decision = artifact["decision"]
    require(decision["corrected_selected_augmented_strip_projection_has_negative_real_sign"] is True, "signed strip decision drift")
    require(decision["corrected_selected_component_embedded_in_fold_owned_paired_residual"] is False, "unproved paired-residual embedding was promoted")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file() and file_hash(path) == record["sha256"], f"dependency hash drift: {path}")

    require(C % 8 == 1 and ((C * C - 1) // 8) % 2 == 0, "independent parity check failed")
    # Reconstruct more sharply than the builder's saved 80-digit balls so
    # overlap tests compare the mathematics, not two unequal outward roundings.
    ctx.dps = 100
    splice_path = REPO_ROOT / artifact["dependencies"]["corrected_selected_splice"]["path"]
    splice = json.loads(splice_path.read_text(encoding="utf-8"))["certificate"]
    factor = 4 * arb(2).sqrt()
    real_upper = (factor * arb(splice["corrected_selected_uniform_real_part_upper"])).upper()
    modulus = (factor * arb(splice["corrected_selected_uniform_modulus_bound"])).upper()
    saved = artifact["certificate"]
    require(arb(saved["source_projected_corrected_real_part_upper"]).overlaps(real_upper), "source real projection arithmetic drift")
    require(arb(saved["source_projected_corrected_modulus_bound"]).overlaps(modulus), "source modulus arithmetic drift")
    require(real_upper < -arb("1.263e-4"), "independent source-projected sign failed")
    print("PASS: independent parity, normalization, and source-projected corrected sign check")


if __name__ == "__main__":
    main()
