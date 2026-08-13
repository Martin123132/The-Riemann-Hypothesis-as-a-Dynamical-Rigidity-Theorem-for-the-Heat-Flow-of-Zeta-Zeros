#!/usr/bin/env python3
"""Join the certified saved-height B-face window channels into one budget."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_window_partial_budget_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "signed_first_current": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_pole_subtracted_cotangent_signed_window_gate.json",
    "higher_currents": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_higher_boundary_current_absolute_window_gate.json",
    "nonlocal_remainder": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_gaussian_window_fresnel_remainder_gate.json",
    "truncation_dictionary": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_fresnel_boundary_truncation_dictionary_gate.json",
    "local_normal_remainder": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_local_normal_safe_fresnel_remainder_gate.json",
}

PRECISION = 100
T = 10_000_000_000
B = 5_122_421
WINDOW_XI = 70
REFERENCE_TARGET = "8.6e-6"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def interval_certificate(dependencies: dict[str, dict[str, Any]]) -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1

    remainder_certificate = dependencies["nonlocal_remainder"]["interval_certificate"]
    x_low = arb(remainder_certificate["x_window_low_ball"])
    x_high = arb(remainder_certificate["x_window_high_ball"])
    remainder_density = arb(remainder_certificate["infinite_nonlocal_completed_remainder_ball"])
    require(x_low > 0 and x_high < arb(1) / 2 and x_low < x_high, "B window geometry drift")

    # The Kummer weight decreases on 0<x<1/2, so its window maximum is at x_low.
    width = x_high - x_low
    weight_max = (x_low * (1 - x_low)) ** (-arb(1) / 4)
    physical_factor = 2 * (arb.pi() / (32 * arb(T))) ** (arb(1) / 4)
    nonlocal_remainder_physical = physical_factor * width * weight_max * remainder_density

    first = arb(
        dependencies["signed_first_current"]["numerical_certificate"]
        ["physically_normalized_signed_projection_ball"]
    )
    higher = arb(
        dependencies["higher_currents"]["interval_certificate"]
        ["combined_higher_current_physical_absolute_integral_ball"]
    )
    local_normal = arb(
        dependencies["local_normal_remainder"]["interval_certificate"]
        ["combined_physically_normalized_remainder_ball"]
    )
    dictionary_correction = arb(
        dependencies["truncation_dictionary"]["interval_certificate"]
        ["physical_dictionary_correction_ball"]
    )
    partial = first + higher + local_normal + nonlocal_remainder_physical + dictionary_correction
    target = arb(REFERENCE_TARGET)
    remaining = target - partial

    require(nonlocal_remainder_physical < arb("1.03e-13"), "physical nonlocal remainder exceeds 1.03e-13")
    require(partial < arb("6.276e-6"), "certified partial B-window budget exceeds 6.276e-6")
    require(remaining > arb("2.324e-6"), "remaining reference allowance fell below 2.324e-6")

    return {
        "height": T,
        "endpoint": B,
        "window_xi": [-WINDOW_XI, WINDOW_XI],
        "x_window": [x_low.str(PRECISION, more=True), x_high.str(PRECISION, more=True)],
        "window_width_ball": width.str(PRECISION, more=True),
        "kummer_weight_max_ball": weight_max.str(PRECISION, more=True),
        "kummer_weight_max_location": "x_window_low",
        "physical_projection_factor_ball": physical_factor.str(PRECISION, more=True),
        "nonlocal_post_third_current_density_ball": remainder_density.str(PRECISION, more=True),
        "nonlocal_post_third_current_physical_ball": nonlocal_remainder_physical.str(PRECISION, more=True),
        "signed_first_current_physical_ball": first.str(PRECISION, more=True),
        "higher_currents_physical_absolute_ball": higher.str(PRECISION, more=True),
        "local_normal_remainder_physical_absolute_ball": local_normal.str(PRECISION, more=True),
        "truncation_dictionary_physical_absolute_ball": dictionary_correction.str(PRECISION, more=True),
        "certified_partial_B_window_upper_ball": partial.str(PRECISION, more=True),
        "certified_partial_B_window_upper_bound": "6.276e-6",
        "reference_physical_target": REFERENCE_TARGET,
        "remaining_triangle_budget_ball": remaining.str(PRECISION, more=True),
        "remaining_triangle_budget_lower_bound": "2.324e-6",
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    return f"""# Saved-height B-face window partial budget

Date: 2026-08-13

Status: rigorous join of four already certified window channels; not a
complete B-face estimate

Let `R_3(x)` denote the absolutely summed nonlocal remainder left after the
third B-face integration-by-parts current.  The earlier Fresnel gate proves
uniformly on `|xi|<=70` that

```text
|R_3(x)| <= {c['nonlocal_post_third_current_density_ball']}.       (BP1)
```

The physical equation-(9) projection contributes the factor

```text
N_t=2(pi/(32t))^(1/4),
W_t(x)=[x(1-x)]^(-1/4).                              (BP2)
```

This window lies below `x=1/2`.  Therefore `x(1-x)` is increasing and
`W_t` is decreasing, so its maximum is attained at `x_low`, not `x_high`.
Using the certified width and `W_t(x_low)` gives

```text
N_t <= {c['physical_projection_factor_ball']},
sup W_t <= {c['kummer_weight_max_ball']},
E_rem <= {c['nonlocal_post_third_current_physical_ball']}
      < 1.03e-13.                                    (BP3)
```

This closes the physical normalization of the nonlocal post-third-current
remainder on the Gaussian window.  It is then joined to the signed first
current, the absolute second/third currents, and the normal-safe local
Fresnel remainders without replacing the signed first current by its much
larger complex modulus:

```text
signed first current                 {c['signed_first_current_physical_ball']}
higher-current absolute allowance    {c['higher_currents_physical_absolute_ball']}
local normal-safe remainder          {c['local_normal_remainder_physical_absolute_ball']}
nonlocal post-third remainder        {c['nonlocal_post_third_current_physical_ball']}
Fresnel/boundary dictionary          {c['truncation_dictionary_physical_absolute_ball']}
--------------------------------------------------------------------------
certified partial upper allowance    {c['certified_partial_B_window_upper_ball']}
                                   < 6.276e-6.        (BP4)
```

Against the working physical target `8.6e-6`, an absolute-value treatment
of every still-open channel therefore has at least

```text
{c['remaining_triangle_budget_ball']} > 2.324e-6      (BP5)
```

available.  This is only accounting headroom.  The exact complementary
outer-safe currents for modes 621 and 622 and both `|xi|>70` grouped trace
tails remain unbounded, so BP4 is not an upper bound for the complete B face.

Pi provenance: `pi` in BP2 is inherited directly from the equation-(9)
normalization.  The Fresnel remainder in BP1 inherits its powers of `pi`
from repeated integration by parts of `exp(i*pi*u^2/2)`.  No fitted or
geometric surrogate for `pi` is introduced.

Proof boundary: four specified saved-height `|xi|<=70` B-face channels
only.  No complementary local current, outside-window tail, complete B
estimate, A-fold splice, complete paired residual, complete `T_upper`,
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["signed_first_current"]["decision"]["first_pole_subtracted_current_integrated_before_absolute_value"] is True, "signed-first-current dependency drift")
    require(dependencies["higher_currents"]["decision"]["combined_higher_window_current_below_1p491e_minus_6"] is True, "higher-current dependency drift")
    require(dependencies["higher_currents"]["interval_certificate"].get("kummer_weight_max_location") == "x_window_low", "higher-current weight endpoint drift")
    require(dependencies["nonlocal_remainder"]["decision"]["nonlocal_completed_remainder_absolutely_summable"] is True, "nonlocal-remainder dependency drift")
    require(dependencies["truncation_dictionary"]["decision"]["exact_rational_dictionary_correction_derived"] is True, "truncation-dictionary dependency drift")
    require(dependencies["local_normal_remainder"]["decision"]["combined_physical_remainder_below_1p97e_minus_10"] is True, "local-remainder dependency drift")

    artifact = {
        "kind": STEM,
        "status": "four_certified_B_window_channels_below_6p276e_minus_6_with_more_than_2p324e_minus_6_reference_allowance_remaining",
        "passed": True,
        "symbolic_certificate": {
            "physical_remainder_bound": "E_rem<=2(pi/(32t))^(1/4)(x_high-x_low)[x_low(1-x_low)]^(-1/4)sup|R_3|",
            "weight_monotonicity": "d[x(1-x)]/dx=1-2x>0 on the saved window, hence [x(1-x)]^(-1/4) is maximized at x_low",
            "partial_budget": "E_partial=E_C+|E_1|+|E_2|+|E_local_normal|+|E_rem|",
            "signed_channel_guard": "E_C is inserted as its certified real-projection interval, not as the complex modulus",
        },
        "interval_certificate": interval_certificate(dependencies),
        "decision": {
            "post_third_current_nonlocal_window_remainder_physically_bounded": True,
            "weight_maximum_assigned_to_x_window_low": True,
            "signed_first_current_preserved_before_absolute_value": True,
            "four_certified_window_channels_joined": True,
            "Fresnel_boundary_truncation_dictionary_correction_included": True,
            "certified_partial_B_window_budget_below_6p276e_minus_6": True,
            "remaining_reference_triangle_budget_above_2p324e_minus_6": True,
            "complementary_outer_safe_local_currents_bounded": False,
            "outside_window_grouped_trace_tails_bounded": False,
            "complete_B_face_estimate_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "flint_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Bound the exact complementary outer-safe 621/622 local currents on |xi|<=70, preserving their step-tail cancellation, then bound both |xi|>70 grouped trace tails. The sum of absolute allowances for those open channels must fit below the certified remaining 2.324e-6 reference budget.",
        "proof_boundary": "Only the signed first current, absolute higher currents, normal-safe local Fresnel remainders, and physically normalized nonlocal post-third-current remainder on the saved-height |xi|<=70 B window. No complementary local current, outside-window tail, complete B estimate, A-fold splice, complete paired residual, complete T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level theorem is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified partial B-window allowance below 6.276e-6", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
