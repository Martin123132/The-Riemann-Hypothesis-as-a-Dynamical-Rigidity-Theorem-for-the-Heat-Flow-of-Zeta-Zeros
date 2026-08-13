#!/usr/bin/env python3
"""Validate the t=1e10 cross-block hybrid exact output gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_hybrid_exact_output_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_hybrid_exact_output_gate.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_crossblock_hybrid_exact_output_gate",
        "kind drift",
    )
    require(len(artifact["output_rows"]) == 15, "output count drift")
    aggregate = artifact["aggregate"]
    require(aggregate["total_call_count"] == 6784, "total call count drift")
    require(aggregate["exact_block20_call_count"] == 374, "exact block-20 count drift")
    require(aggregate["retained_later_call_count"] == 1040, "retained later count drift")
    require(aggregate["uses_signed_cross_call_cancellation"] is False, "signed cancellation was used")
    require(aggregate["block20_complex_intervals_decoded_from_dyadic_endpoints"], "dyadic decode flag drift")
    require(aggregate["outputs_within_requested_scale"] == 15, "not all outputs close")
    for row in artifact["output_rows"]:
        block20 = Fraction(row["block20_exact_correction_majorant_upper"])
        later = Fraction(row["blocks21_28_retained_majorant_upper"])
        direct = Fraction(row["direct_kernel_majorant_upper"])
        combined = Fraction(row["hybrid_complete_majorant_upper"])
        require(
            abs(block20 + later + direct - combined) <= Fraction(1, 10**49),
            f"hybrid decomposition drift at output {row['output_index']}",
        )
        require(combined < Fraction(5, 1000), f"hybrid scale failure at output {row['output_index']}")
        require(Fraction(row["margin_below_requested_scale"]) > 0, "nonpositive margin")
    for record in (*artifact["dependencies"].values(), *artifact["sources"].values()):
        path = REPO_ROOT / record["path"]
        require(path.is_file(), f"missing artifact: {path}")
        require(file_hash(path) == record["sha256"], f"artifact hash drift: {path}")
    require(NOTE.is_file(), "missing markdown note")
    print(
        "validated t=1e10 cross-block hybrid exact output gate: "
        f"6784 calls, max {aggregate['maximum_hybrid_complete_majorant_upper']}, 15/15 below 0.005"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
