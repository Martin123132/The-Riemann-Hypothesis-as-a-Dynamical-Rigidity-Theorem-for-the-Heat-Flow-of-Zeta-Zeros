#!/usr/bin/env python3
"""Join the certified B trace-current package on the saved |xi|<=70 window."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_complete_gaussian_window_signed_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "signed_nonlocal_first": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_pole_subtracted_cotangent_signed_window_gate.json",
    "nonlocal_higher": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_higher_boundary_current_absolute_window_gate.json",
    "partial_budget": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_window_partial_budget_gate.json",
    "local_exact": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_local_outer_safe_exact_current_gate.json",
}

PRECISION = 90
T = 10_000_000_000
B = 5_122_421


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
    first = arb(
        dependencies["signed_nonlocal_first"]["numerical_certificate"]
        ["physically_normalized_signed_projection_ball"]
    )
    higher = arb(
        dependencies["nonlocal_higher"]["interval_certificate"]
        ["combined_higher_current_physical_absolute_integral_ball"]
    )
    partial = dependencies["partial_budget"]["interval_certificate"]
    nonlocal_remainder = arb(partial["nonlocal_post_third_current_physical_ball"])
    dictionary = arb(partial["truncation_dictionary_physical_absolute_ball"])
    local = arb(
        dependencies["local_exact"]["interval_certificate"]
        ["exact_local_B_window_physical_upper_ball"]
    )

    # The local theorem already includes the complete positive normal-safe
    # remainder allowance, so that term from the earlier partial ledger is
    # deliberately not charged a second time here.
    upper = first + higher + nonlocal_remainder + dictionary + local
    positive_nonlocal_allowance = first + higher + nonlocal_remainder + dictionary
    require(positive_nonlocal_allowance < arb("6.276e-6"), "nonlocal window allowance exceeds 6.276e-6")
    require(upper < arb("-1.3198e-4"), "complete B trace-package window bound did not close")
    return {
        "height": T,
        "endpoint": B,
        "window_xi": [-70, 70],
        "signed_nonlocal_first_current_ball": first.str(PRECISION, more=True),
        "nonlocal_higher_current_absolute_allowance_ball": higher.str(PRECISION, more=True),
        "nonlocal_post_third_current_absolute_allowance_ball": nonlocal_remainder.str(PRECISION, more=True),
        "Fresnel_boundary_dictionary_absolute_allowance_ball": dictionary.str(PRECISION, more=True),
        "exact_local_621_622_upper_ball": local.str(PRECISION, more=True),
        "positive_nonlocal_window_allowance_ball": positive_nonlocal_allowance.str(PRECISION, more=True),
        "complete_B_trace_package_window_physical_upper_ball": upper.str(PRECISION, more=True),
        "complete_B_trace_package_window_physical_upper_bound": "-1.3198e-4",
        "local_normal_remainder_count": 1,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    return f"""# Complete signed B trace-current package on the Gaussian window

Date: 2026-08-13

Status: rigorous saved-height `|xi|<=70` B trace-current package; not a proof
of the complete endpoint B sector or RH; the nonlocal constant bulk-step
completion and all other global sectors remain open

The nonlocal B modes are partitioned into the signed first pole-subtracted
current, absolute second/third-current allowance, post-third-current
Fresnel remainder, and the exact Fresnel-to-boundary truncation correction.
The local modes 621 and 622 are supplied by the cancellation-safe exact
step-tail theorem.  Thus the complete trace-package upper estimate is

```text
E_Btr,win^+
 = E_nonlocal,0
   + |E_nonlocal,1:2| + |E_nonlocal,rem| + |E_dictionary|
   + E_local,621:622^+ .                              (BW1)
```

The certified terms are

```text
signed nonlocal first       {c['signed_nonlocal_first_current_ball']}
nonlocal higher allowance   {c['nonlocal_higher_current_absolute_allowance_ball']}
nonlocal Fresnel remainder  {c['nonlocal_post_third_current_absolute_allowance_ball']}
truncation dictionary       {c['Fresnel_boundary_dictionary_absolute_allowance_ball']}
exact local upper estimate  {c['exact_local_621_622_upper_ball']}
--------------------------------------------------------------------------
complete trace-package upper {c['complete_B_trace_package_window_physical_upper_ball']}
                          < -1.3198e-4.               (BW2)
```

There is no double counting in BW2.  The exact local theorem already charges
the full positive normal-safe local Fresnel remainder adversely, so the
same remainder from the earlier partial ledger is not added again.  The
nonlocal Fresnel remainder excludes modes 621 and 622, and the cotangent and
higher-current sums have those modes removed algebraically.

For `tau_m=1_(622<=m<=39894)`, the complete positive B endpoint sector also
contains the step coefficient

```text
sigma_m(x)=(1-tau_m)-1_(q_B,m(x)<0).                  (BW3)
```

BW1 includes `sigma_m P_bulk,m` for the two local crossing modes 621 and
622, because their jumps must remain attached to the exact local tails.  For
every nonlocal mode it contains the sign-adapted trace tail `U_B,m` but not
the constant-on-window term `sigma_m P_bulk,m`.  Those nonlocal steps remain
in the global completion.  Therefore the rigorous negative conclusion is
for the complete trace-current package on the window, not for the complete
endpoint B sector.  The exterior trace package and the omitted nonlocal
steps must be handled in the global grouped remainder.

Pi provenance: all terms retain the equation-(9) normalization and phases;
the `pi/8` projection is the exact odd-square half-Kummer reflection.  No
fitted value of `pi` or empirical phase is introduced.

Proof boundary: the complete saved-height B trace-current package restricted
to `|xi|<=70`, including the exact local crossing steps only.  No bound for
the nonlocal constant bulk-step completion, exterior trace package, complete
endpoint B sector, A-fold splice, complete `Q_K-T` or `T_upper`,
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["signed_nonlocal_first"]["decision"]["first_pole_subtracted_current_integrated_before_absolute_value"] is True, "signed-current dependency drift")
    require(dependencies["nonlocal_higher"]["decision"]["local_621_622_poles_excluded_before_summation"] is True, "higher-current local exclusion drift")
    require(dependencies["partial_budget"]["decision"]["Fresnel_boundary_truncation_dictionary_correction_included"] is True, "partial-budget dictionary drift")
    require(dependencies["local_exact"]["decision"]["exact_local_621_622_window_projection_below_minus_1p369e_minus_4"] is True, "local exact dependency drift")
    artifact = {
        "kind": STEM,
        "status": "complete_saved_height_B_trace_package_on_abs_xi_at_most_70_below_minus_1p3198e_minus_4",
        "passed": True,
        "symbolic_certificate": {
            "step_coefficient": "sigma_m=(1-tau_m)-1_(q_B,m<0), tau_m=1_(622<=m<=39894)",
            "nonlocal_partition": "exact nonlocal B trace package = C_0_hat+C_1_hat+C_2_hat+R_F-D_dictionary; its constant sigma_m P_bulk,m completion remains global",
            "local_partition": "modes 621,622 use the exact step-tail/direct-Fresnel splice, including their local sigma_m P_bulk,m step, with the remainder charged once",
            "disjointness": "all nonlocal sums exclude 621,622; the local theorem contains exactly those two modes",
            "upper_join": "E_Btr,win<=E_nonlocal,0+|E_nonlocal,1:2|+|R_nonlocal|+|D_dictionary|+E_local^+",
        },
        "interval_certificate": interval_certificate(dependencies),
        "decision": {
            "all_nonlocal_B_trace_current_window_channels_covered": True,
            "exact_local_621_622_channels_covered": True,
            "local_normal_remainder_counted_exactly_once": True,
            "complete_saved_height_B_trace_package_window_upper_bound_proved": True,
            "complete_B_trace_package_window_projection_below_minus_1p3198e_minus_4": True,
            "nonlocal_constant_bulk_step_completion_included": False,
            "complete_saved_height_endpoint_B_sector_upper_bound_proved": False,
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
        "next_obligation": "Embed this negative trace-current package in the exact global residual while retaining every nonlocal sigma_m P_bulk,m step in the complement. Bound the exterior trace package and the resulting global joined remainder without assigning the omitted steps an independent norm.",
        "proof_boundary": "Only the complete saved-height B trace-current package on |xi|<=70, including exact local crossing steps. No nonlocal constant bulk-step completion, exterior trace package, complete endpoint B sector, A-fold splice, complete Q_K-T or T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level theorem is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified complete B trace-package window below -1.3198e-4", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
