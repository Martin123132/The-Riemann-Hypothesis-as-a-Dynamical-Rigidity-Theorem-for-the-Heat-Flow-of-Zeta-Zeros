#!/usr/bin/env python3
"""Check the event-parameter common-contour geometry independently."""

from __future__ import annotations

import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))

from flint import arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_contour_geometry_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
C = 159_577
Q = (C - 1) // 4


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def mode(index: int) -> int:
    return Q - index // 2 if index % 2 == 0 else Q + (index + 1) // 2


def main() -> None:
    ctx.dps = 110
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("status") == "event_parameter_range_and_common_contour_geometry_certified", "bad status")
    c = artifact["certificate"]
    pi = arb.pi()
    beta = (pi * arb(C) ** 2 / 8) ** (arb(1) / 3)
    kappa = pi / (16 * beta)
    failures = []
    detunings = []
    for index, row in enumerate(c["roster"]):
        m = mode(index)
        signed = 4 * m - C
        require(row["mode"] == m and row["signed_odd"] == signed, f"roster mismatch at {index}")
        require(signed == (-1 if index % 2 == 0 else 1) * (2 * index + 1), f"odd law mismatch at {index}")
        d = beta * arb(signed) / C
        recorded = arb(row["detuning_ball"])
        require(d.overlaps(recorded), f"detuning mismatch at {index}")
        saddle = (d**2 + kappa + 64).sqrt()
        require(saddle.overlaps(arb(row["canonical_outer_saddle_ball"])), f"saddle mismatch at {index}")
        if saddle.lower() > 9:
            failures.append(index)
        detunings.append(abs(d).upper())

    require(len(failures) == c["old_R9_failure_count"], "R9 failure count mismatch")
    dmax = max(detunings, key=float)
    lammax = dmax**2 + kappa
    canmax = (lammax + 64).sqrt()
    ratio = (beta**2 + 64 + arb(64) ** 2 / (4 * beta**2)) / (beta**2 - dmax**2 - kappa)
    exactmax = beta * ratio.sqrt().acosh()
    require(canmax.overlaps(arb(c["maximum_canonical_saddle_abs_ball"])), "canonical envelope mismatch")
    require(exactmax.overlaps(arb(c["maximum_exact_saddle_abs_ball"])), "exact envelope mismatch")
    require(canmax < 14 and exactmax < 14, "R14 saddle containment failed")

    contour = c["common_contour"]
    for key in (
        "canonical_saddle_clearance_ball",
        "exact_saddle_clearance_ball",
        "canonical_lifted_path_bracket_ball",
        "exact_lifted_path_bracket_ball",
        "minimum_lifted_Cauchy_bracket_ball",
        "exact_far_exponent_margin_ball",
        "canonical_far_exponent_margin_ball",
    ):
        require(arb(contour[key]) > 0, f"nonpositive contour margin: {key}")
    require(arb(contour["maximum_local_tanh_q_radius_ball"]) < arb("0.02"), "tanh radius failed")
    require(artifact["claims"]["uniform_full_line_remainder_proved"] is False, "remainder boundary promoted")
    require(artifact["claims"]["RH_proved"] is False, "RH boundary promoted")
    note = (REPO_ROOT / artifact["artifacts"]["note"]["path"]).read_text(encoding="utf-8")
    for token in ("prototype `R=9` contour is not atlas-uniform", "fixed contour", "does not show"):
        require(token in note, f"missing note token: {token}")
    print(f"validated event contour geometry independently: events=399, R9 failures={len(failures)}, R14 passes")


if __name__ == "__main__":
    main()
