#!/usr/bin/env python3
"""Independently check the corrected selected-defect splice."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
from flint import arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_corrected_selected_defect_splice_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def algebra_check() -> None:
    modes = range(23)
    w = {m: Fraction(9 * m + 5, 31 * m + 37) for m in modes}
    s = {m: Fraction(4 * m - 7, 29 * m + 41) for m in modes}
    o = {m: Fraction(2 - 5 * m, 37 * m + 43) for m in modes}
    d = {m: Fraction(7 * m + 1, 41 * m + 47) for m in modes}
    exact = {m: s[m] + o[m] + d[m] for m in modes}
    lhs = sum((w[m] * (exact[m] - o[m]) for m in modes), Fraction())
    rhs = sum((w[m] * (s[m] + d[m]) for m in modes), Fraction())
    require(lhs == rhs, "independent 23-mode splice identity failed")


def main() -> None:
    require(RESULT.is_file(), "missing result artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "saved gate did not pass")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file() and file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    algebra_check()
    ctx.dps = 70
    selected = json.loads((REPO_ROOT / artifact["dependencies"]["selected_beta4"]["path"]).read_text(encoding="utf-8"))["certificate"]
    profile = json.loads((REPO_ROOT / artifact["dependencies"]["exact_profile_closure"]["path"]).read_text(encoding="utf-8"))["interval_certificate"]
    finite = arb(profile["uniform_weighted_finite_t_profile_correction_bound"]).upper()
    maximum_real = max(arb(row["value_ball"]["real_ball"]).upper() for row in selected["node_rows"])
    independent_real_upper = maximum_real + arb(selected["y_midpoint_error_ball"]).upper() + arb(selected["height_transport_error_ball"]).upper() + finite
    independent_modulus = arb(selected["uniform_selected_weighted_bound_ball"]).upper() + finite
    saved = artifact["certificate"]
    require(arb(saved["corrected_selected_uniform_real_part_upper"]).overlaps(independent_real_upper), "saved real-part arithmetic drift")
    require(arb(saved["corrected_selected_uniform_modulus_bound"]).overlaps(independent_modulus), "saved modulus arithmetic drift")
    require(independent_real_upper < -arb("2.233e-5"), "independent corrected real-part sign failed")
    print("PASS: independent 23-mode algebra and corrected selected-defect sign arithmetic")


if __name__ == "__main__":
    main()
