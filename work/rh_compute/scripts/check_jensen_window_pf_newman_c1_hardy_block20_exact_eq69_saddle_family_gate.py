"""Validate the block-20 exact equation-(69) saddle-family certificate."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_exact_eq69_saddle_family_gate.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(RESULT.is_file(), "missing exact equation-(69) saddle-family artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_block20_exact_eq69_saddle_family_gate", "equation-(69) artifact kind drift")
    require(artifact["status"] == "rigorous_finite_exact_equation_69_saddle_family_enclosed", "equation-(69) artifact status drift")
    scope = artifact["scope"]
    aggregate = artifact["aggregate"]
    rows = artifact["rows"]
    require(scope["recursive_call_count"] == 374 and len(rows) == 374, "equation-(69) row count drift")
    require(scope["exact_integral_count"] == aggregate["exact_integral_count"] == 544, "equation-(69) integral count drift")
    require(scope["precision_overlap_count"] == aggregate["precision_overlap_count"] == 374, "equation-(69) precision overlap drift")
    require(scope["precision_ladder_decimal_digits"] == [70, 110], "equation-(69) precision ladder drift")
    require(aggregate["residual_after_excluding_zero_count"] == 374, "equation-(69) residual exclusion drift")
    require(
        aggregate["improved_call_count"] + aggregate["worsened_call_count"] + aggregate["indeterminate_call_count"] == 374,
        "equation-(69) classification count drift",
    )
    require(Decimal(aggregate["minimum_residual_after_magnitude_lower"]) > 0, "equation-(69) zero residual admitted")
    require(Decimal(aggregate["maximum_phase_transport_bound_upper"]) < Decimal("1e-25"), "equation-(69) pi transport scale drift")
    require(Decimal(aggregate["maximum_cancellation_identity_gap_upper"]) < Decimal("1e-80"), "equation-(69) cancellation identity drift")
    require(all(row["precision_overlap"] for row in rows), "equation-(69) row precision drift")
    require(all(row["residual_after_excludes_zero"] for row in rows), "equation-(69) row residual drift")
    require(sum(len(row["indices"]) for row in rows) == 544, "equation-(69) saved index count drift")
    require(all(row["indices"] in ([0, 1], [1], [1, 2]) for row in rows), "equation-(69) index roster drift")
    require(all(bool(row["child_subtract_one"]) == (row["indices"][0] == 1) for row in rows), "equation-(69) start adapter drift")
    for group in (artifact["source"], artifact["paper"]["pdf"], *artifact["sources"].values()):
        if "path" not in group or "sha256" not in group:
            continue
        path = REPO_ROOT / group["path"]
        require(path.is_file(), f"equation-(69) dependency missing: {group['path']}")
        require(file_hash(path) == group["sha256"], f"equation-(69) dependency hash drift: {group['path']}")
    require("not a height-uniform recurrence theorem" in artifact["proof_boundary"], "equation-(69) proof boundary drift")
    print(
        "validated Hardy block-20 exact equation-(69) saddle-family gate: "
        f"544 integrals, {aggregate['improved_call_count']} improved, 374 residuals exclude zero"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
