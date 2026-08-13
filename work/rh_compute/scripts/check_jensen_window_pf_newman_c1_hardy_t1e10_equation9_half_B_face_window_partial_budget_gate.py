#!/usr/bin/env python3
"""Independent replay of the saved-height partial B-window budget."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_window_partial_budget_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()

T = 10_000_000_000
B = 5_122_421


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ball_bounds(text: str) -> tuple[mp.mpf, mp.mpf]:
    match = re.fullmatch(r"\[([^\]]+) \+/- ([^\]]+)\]", text)
    require(match is not None, f"unparsed Arb ball: {text}")
    midpoint, radius = (mp.mpf(value) for value in match.groups())
    return midpoint - radius, midpoint + radius


def ball_upper(text: str) -> mp.mpf:
    return ball_bounds(text)[1]


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    dependencies = {
        name: json.loads((REPO_ROOT / record["path"]).read_text(encoding="utf-8"))
        for name, record in artifact["dependencies"].items()
    }
    mp.mp.dps = 100
    t, endpoint = mp.mpf(T), mp.mpf(B)
    x0 = (1 - mp.sqrt(1 - 8 * t / (mp.pi * endpoint**2))) / 2
    hessian = t * (1 - 2 * x0) / (2 * x0**2 * (1 - x0) ** 2)
    root_hessian = mp.sqrt(hessian)
    x_low = x0 - 70 / root_hessian
    x_high = x0 + 70 / root_hessian
    require(0 < x_low < x_high < mp.mpf("0.5"), "independent window reconstruction failed")

    remainder_upper = ball_upper(
        dependencies["nonlocal_remainder"]["interval_certificate"]
        ["infinite_nonlocal_completed_remainder_ball"]
    )
    physical_factor = 2 * (mp.pi / (32 * t)) ** mp.mpf("0.25")
    weight_max = (x_low * (1 - x_low)) ** mp.mpf("-0.25")
    nonlocal_physical = physical_factor * (x_high - x_low) * weight_max * remainder_upper
    require(nonlocal_physical < mp.mpf("1.03e-13"), "independent physical remainder bound failed")

    first_upper = ball_upper(
        dependencies["signed_first_current"]["numerical_certificate"]
        ["physically_normalized_signed_projection_ball"]
    )
    higher_upper = ball_upper(
        dependencies["higher_currents"]["interval_certificate"]
        ["combined_higher_current_physical_absolute_integral_ball"]
    )
    local_upper = ball_upper(
        dependencies["local_normal_remainder"]["interval_certificate"]
        ["combined_physically_normalized_remainder_ball"]
    )
    dictionary_upper = ball_upper(
        dependencies["truncation_dictionary"]["interval_certificate"]
        ["physical_dictionary_correction_ball"]
    )
    independent_partial_upper = first_upper + higher_upper + local_upper + nonlocal_physical + dictionary_upper
    independent_remaining_lower = mp.mpf("8.6e-6") - independent_partial_upper
    require(independent_partial_upper < mp.mpf("6.276e-6"), "independent partial budget failed")
    require(independent_remaining_lower > mp.mpf("2.324e-6"), "independent remaining budget failed")

    certificate = artifact["interval_certificate"]
    require(certificate["kummer_weight_max_location"] == "x_window_low", "stored weight endpoint drift")
    require(ball_upper(certificate["nonlocal_post_third_current_physical_ball"]) < mp.mpf("1.03e-13"), "stored physical remainder too large")
    require(ball_upper(certificate["certified_partial_B_window_upper_ball"]) < mp.mpf("6.276e-6"), "stored partial budget too large")
    require(ball_bounds(certificate["remaining_triangle_budget_ball"])[0] > mp.mpf("2.324e-6"), "stored remaining budget too small")

    decision = artifact["decision"]
    require(decision["signed_first_current_preserved_before_absolute_value"] is True, "signed-channel guard drift")
    require(decision["complementary_outer_safe_local_currents_bounded"] is False, "local-current overclaim")
    require(decision["outside_window_grouped_trace_tails_bounded"] is False, "outside-window overclaim")
    require(decision["complete_B_face_estimate_proved"] is False, "complete-B overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print(
        "independently checked partial B-window budget: "
        f"remainder={mp.nstr(nonlocal_physical, 16)}, "
        f"partial={mp.nstr(independent_partial_upper, 16)}, "
        f"remaining={mp.nstr(independent_remaining_lower, 16)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
