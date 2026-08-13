"""Validate the block-20 t2 L+1 source/paper mismatch certificate."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_t2_lplus1_mismatch_gate.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(RESULT.is_file(), "missing t2 L+1 mismatch artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_block20_t2_lplus1_mismatch_gate", "t2 L+1 artifact kind drift")
    require(artifact["status"] == "exact_source_paper_t2_lplus1_mismatch_isolated_and_finitely_corrected", "t2 L+1 artifact status drift")
    scope = artifact["scope"]
    aggregate = artifact["aggregate"]
    rows = artifact["rows"]
    require(scope["recursive_call_count"] == len(rows) == 374, "t2 L+1 row count drift")
    require(scope["precision_overlap_count"] == 374 and scope["precision_ladder_decimal_digits"] == [90, 150], "t2 L+1 precision drift")
    require(aggregate["improved_call_count"] == 374, "t2 L+1 improvement count drift")
    require(aggregate["worsened_call_count"] == aggregate["indeterminate_call_count"] == 0, "t2 L+1 classification drift")
    require(aggregate["residual_after_excluding_zero_count"] == 374, "t2 L+1 residual exclusion drift")
    require(Decimal(aggregate["minimum_residual_after_magnitude_lower"]) > Decimal("1e-7"), "t2 L+1 minimum residual drift")
    require(Decimal(aggregate["maximum_residual_after_magnitude_upper"]) < Decimal("0.0016"), "t2 L+1 maximum residual drift")
    require(Decimal(aggregate["maximum_residual_ratio_to_prior_upper"]) < Decimal("0.02"), "t2 L+1 reduction ratio drift")
    require(Decimal(aggregate["maximum_pi_transport_bound_upper"]) < Decimal("1e-30"), "t2 L+1 pi transport drift")
    require(Decimal(aggregate["maximum_derivation_identity_gap_upper"]) < Decimal("1e-120"), "t2 L+1 algebra identity drift")
    require(all(row["precision_overlap"] and row["magnitude_classification"] == "improved" for row in rows), "t2 L+1 row classification drift")
    require(all(row["residual_after_excludes_zero"] for row in rows), "t2 L+1 row residual drift")
    require(all(row["child_length"] in (1, 2) for row in rows), "t2 L+1 child roster drift")
    for group in (artifact["source"], artifact["paper"]["pdf"], *artifact["sources"].values()):
        if "path" not in group or "sha256" not in group:
            continue
        path = REPO_ROOT / group["path"]
        require(path.is_file(), f"t2 L+1 dependency missing: {group['path']}")
        require(file_hash(path) == group["sha256"], f"t2 L+1 dependency hash drift: {group['path']}")
    require("not a proof of the remaining nonsaddle estimates" in artifact["proof_boundary"], "t2 L+1 proof boundary drift")
    print(
        "validated Hardy block-20 t2 L+1 mismatch gate: "
        "374 improved, max residual ratio below 0.02, 374 residuals exclude zero"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
