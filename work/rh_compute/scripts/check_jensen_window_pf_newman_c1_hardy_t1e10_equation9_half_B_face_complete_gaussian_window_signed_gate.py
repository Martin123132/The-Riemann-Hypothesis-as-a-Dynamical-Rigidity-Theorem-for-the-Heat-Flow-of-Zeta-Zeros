#!/usr/bin/env python3
"""Independent accounting replay of the complete B trace package."""

from __future__ import annotations

import hashlib
import json
import re
from decimal import Decimal, getcontext
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_complete_gaussian_window_signed_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bounds(text: str) -> tuple[Decimal, Decimal]:
    match = re.fullmatch(r"\[([^\]]+) \+/- ([^\]]+)\]", text)
    require(match is not None, f"unparsed Arb ball: {text}")
    midpoint, radius = (Decimal(value) for value in match.groups())
    return midpoint - radius, midpoint + radius


def upper(text: str) -> Decimal:
    return bounds(text)[1]


def main() -> int:
    getcontext().prec = 90
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    deps = {
        name: json.loads((REPO_ROOT / record["path"]).read_text(encoding="utf-8"))
        for name, record in artifact["dependencies"].items()
    }
    first = upper(deps["signed_nonlocal_first"]["numerical_certificate"]["physically_normalized_signed_projection_ball"])
    higher = upper(deps["nonlocal_higher"]["interval_certificate"]["combined_higher_current_physical_absolute_integral_ball"])
    partial = deps["partial_budget"]["interval_certificate"]
    remainder = upper(partial["nonlocal_post_third_current_physical_ball"])
    dictionary = upper(partial["truncation_dictionary_physical_absolute_ball"])
    local = upper(deps["local_exact"]["interval_certificate"]["exact_local_B_window_physical_upper_ball"])
    independent_upper = first + higher + remainder + dictionary + local
    require(independent_upper < Decimal("-1.3198e-4"), "independent signed join failed")

    certificate = artifact["interval_certificate"]
    require(certificate["local_normal_remainder_count"] == 1, "local remainder multiplicity drift")
    require(upper(certificate["complete_B_trace_package_window_physical_upper_ball"]) < Decimal("-1.3198e-4"), "stored signed upper bound failed")
    decision = artifact["decision"]
    require(decision["local_normal_remainder_counted_exactly_once"] is True, "double-count guard drift")
    require(decision["nonlocal_constant_bulk_step_completion_included"] is False, "nonlocal step completion silently included")
    require(decision["complete_saved_height_endpoint_B_sector_upper_bound_proved"] is False, "endpoint B-sector overclaim")
    require(decision["outside_window_grouped_trace_tails_bounded"] is False, "outside-window overclaim")
    require(decision["complete_B_face_estimate_proved"] is False, "complete-B overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print(f"independently checked complete B trace-package window upper={independent_upper}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
