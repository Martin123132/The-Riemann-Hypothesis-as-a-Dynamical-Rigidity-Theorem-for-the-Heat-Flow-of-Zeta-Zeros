#!/usr/bin/env python3
"""Validate the displayed-W1 finite-ip completion artifact."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_w1_finite_ip_completion_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_w1_finite_ip_completion_gate.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(RESULT.is_file() and NOTE.is_file(), "missing finite-ip gate output")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact["status"] == "rigorous_displayed_w1_finite_index_completion_enclosed", "finite-ip status drift")
    require(artifact["scope"]["recursive_call_count"] == 374, "finite-ip call count drift")
    require(artifact["scope"]["precision_overlap_count"] == 374, "finite-ip precision count drift")
    require(len(artifact["rows"]) == 374, "finite-ip row count drift")
    aggregate = artifact["aggregate"]
    require(
        aggregate["fully_covered_call_count"] + aggregate["calls_with_missing_indices"] == 374,
        "finite-ip coverage partition drift",
    )
    require(
        aggregate["missing_lower_term_count"] + aggregate["missing_upper_term_count"]
        == aggregate["missing_term_count"],
        "finite-ip term partition drift",
    )
    require(
        aggregate["improved_call_count"] + aggregate["worsened_call_count"]
        + aggregate["indeterminate_call_count"] == 374,
        "finite-ip classification partition drift",
    )
    require(all(row["precision_overlap"] for row in artifact["rows"]), "finite-ip precision flag drift")
    require(
        sum(bool(row["missing_lower_indices"] or row["missing_upper_indices"]) for row in artifact["rows"])
        == aggregate["calls_with_missing_indices"],
        "finite-ip missing-row recount drift",
    )
    require(
        sum(row["residual_after_excludes_zero"] for row in artifact["rows"])
        == aggregate["residual_after_excluding_zero_count"],
        "finite-ip nonzero recount drift",
    )

    source = artifact["source"]
    source_path = REPO_ROOT / source["path"]
    require(source_path.is_file() and file_hash(source_path) == source["sha256"], "finite-ip source hash drift")
    for name, record in artifact["sources"].items():
        if not isinstance(record, dict) or "path" not in record:
            continue
        path = REPO_ROOT / record["path"]
        require(path.is_file(), f"missing finite-ip dependency: {name}")
        require(file_hash(path) == record["sha256"], f"finite-ip dependency hash drift: {name}")

    note = " ".join(NOTE.read_text(encoding="utf-8").split())
    for token in (
        "equation (69)",
        "equation (81)",
        "ip1=1",
        "R_after = R_q - T_ip",
        "not used as an analytic contour-error bound",
        "does not turn equations (72), (75), or (77) into equalities",
    ):
        require(token in note, f"finite-ip note token missing: {token}")
    print(
        "validated Hardy block-20 displayed-W1 finite-ip completion gate: "
        f"{aggregate['fully_covered_call_count']} fully covered, "
        f"{aggregate['missing_term_count']} omitted terms, "
        f"{aggregate['residual_after_excluding_zero_count']} residuals exclude zero"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
