#!/usr/bin/env python3
"""Guard the B trace package from promotion to the complete endpoint sector."""

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
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_trace_package_bulk_step_scope_guard_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "step_tail": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_positive_B_crossing_step_tail_normal_form_gate.json",
    "window_geometry": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_gaussian_window_fresnel_remainder_gate.json",
    "trace_window": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_complete_gaussian_window_signed_gate.json",
    "global_extraction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_B_window_extraction_target_gate.json",
}


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


def symbolic_certificate() -> dict[str, Any]:
    bulk, A_plus, U_plus, I_minus, B_minus = sp.symbols(
        "P_bulk P_A_plus U_B_plus I_minus P_B_minus"
    )
    tau, local, q_negative = sp.symbols("tau ell qneg")
    sigma = 1 - tau - q_negative
    pair = A_plus + I_minus + U_plus + sigma * bulk
    trace = U_plus + B_minus + local * sigma * bulk
    grouped = A_plus + I_minus - B_minus + (1 - local) * sigma * bulk
    require(sp.expand(pair - trace - grouped) == 0, "trace/grouped identity failed")

    rows = [
        {"mode_range": "1..620", "tau": 0, "ell": 0, "q_negative": 0, "sigma": 1, "owner": "grouped_remainder"},
        {"mode_range": "621 at left/right window sides", "tau": 0, "ell": 1, "q_negative": "1/0", "sigma": "0/1", "owner": "local_trace_package"},
        {"mode_range": "622 at left/right window sides", "tau": 1, "ell": 1, "q_negative": "1/0", "sigma": "-1/0", "owner": "local_trace_package"},
        {"mode_range": "623..39894", "tau": 1, "ell": 0, "q_negative": 1, "sigma": -1, "owner": "grouped_remainder"},
        {"mode_range": "39895 and above", "tau": 0, "ell": 0, "q_negative": 1, "sigma": 0, "owner": "zero"},
    ]
    return {
        "step_coefficient": "sigma_m=(1-tau_m)-1_(q_B,m<0), tau_m=1_(622<=m<=39894), ell_m=1_(m in {621,622})",
        "pair_identity": "I_m+I_-m-tau_m P_bulk=[U_B,m+P_B,-m+ell_m sigma_m P_bulk]+[P_A,m+I_-m-P_B,-m+(1-ell_m)sigma_m P_bulk]",
        "trace_package": "Btr_M=sum_[U_B,m+P_B,-m+ell_m sigma_m P_bulk,m]",
        "nonlocal_window_step_block": "sum_(m=1)^620 P_bulk,m-sum_(m=623)^39894 P_bulk,m",
        "window_case_table": rows,
    }


def interval_certificate(dependencies: dict[str, dict[str, Any]]) -> dict[str, str]:
    ctx.dps = 80
    ctx.threads = 1
    geometry = dependencies["window_geometry"]["interval_certificate"]
    c_low = arb(geometry["c_window_low_ball"])
    c_high = arb(geometry["c_window_high_ball"])
    require(c_low > arb(620) and c_low < arb(621), "left window lattice position drift")
    require(c_high > arb(622) and c_high < arb(623), "right window lattice position drift")
    return {
        "c_window_low_ball": c_low.str(80, more=True),
        "c_window_high_ball": c_high.str(80, more=True),
        "certified_nonlocal_positive_q_range": "m<=620",
        "certified_local_crossing_modes": "621,622",
        "certified_nonlocal_negative_q_range": "m>=623",
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    return f"""# B trace-package bulk-step scope guard

Date: 2026-08-13

Status: exact ownership audit; not a proof of the complete endpoint-B window
or RH; the negative trace-package bound survives while that promotion is
rejected

With

```text
tau_m=1_(622<=m<=39894),
ell_m=1_(m in {{621,622}}),
sigma_m=(1-tau_m)-1_(q_B,m<0),
```

the exact finite-mode identity is

```text
I_m+I_-m-tau_m P_bulk,m
 =[U_B,m+P_B,-m+ell_m sigma_m P_bulk,m]
  +[P_A,m+I_-m-P_B,-m
    +(1-ell_m)sigma_m P_bulk,m].                     (SG1)
```

The first bracket is the certified trace package.  The second bracket keeps
every nonlocal bulk step in the grouped remainder.  The window geometry

```text
c_low={c['c_window_low_ball']},
c_high={c['c_window_high_ball']}
```

places exactly modes 621 and 622 across a crossing.  Hence on the complete
window the omitted nonlocal step block is exactly

```text
sum_(m=1)^620 P_bulk,m - sum_(m=623)^39894 P_bulk,m. (SG2)
```

Modes 621 and 622 retain their variable steps in the local trace package;
modes at least 39895 have `sigma_m=0`.  Therefore the already-certified
bound

```text
E_Btr,win<-1.3198e-4
```

is preserved, but it is not a bound for the complete endpoint-B sector.
No estimate of (SG2) is asserted; it stays in `R_group` under the common
symmetric Abel prescription.

Proof boundary: exact finite-mode ownership, window sign roster, and a
fail-closed non-promotion guard only.  No independent bound for the omitted
bulk-step block, grouped remainder, complete endpoint B sector, complete
`Q_K-T` or `T_upper`, height-uniform theorem, `Lambda<=0`, PF-infinity, RH,
or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["step_tail"]["decision"]["B_endpoint_rewritten_as_sign_adapted_tail_plus_bulk_step"] is True, "step-tail identity drift")
    require(dependencies["window_geometry"]["decision"]["xi_70_window_contains_exactly_B_crossings_621_622"] is True, "window crossing roster drift")
    require(dependencies["trace_window"]["decision"]["complete_saved_height_B_trace_package_window_upper_bound_proved"] is True, "trace-package bound drift")
    require(dependencies["trace_window"]["decision"]["nonlocal_constant_bulk_step_completion_included"] is False, "trace-package scope drift")
    require(dependencies["global_extraction"]["decision"]["nonlocal_constant_bulk_steps_retained_in_grouped_remainder"] is True, "global ownership drift")
    require(dependencies["global_extraction"]["decision"]["complete_endpoint_B_window_extracted"] is False, "endpoint-sector overclaim drift")

    artifact = {
        "kind": STEM,
        "status": "B_trace_package_scope_guard_with_nonlocal_bulk_steps_retained_globally",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(dependencies),
        "decision": {
            "finite_mode_trace_grouped_identity_proved": True,
            "window_nonlocal_sigma_ranges_proved": True,
            "local_621_622_crossing_steps_owned_by_trace_package": True,
            "nonlocal_bulk_steps_owned_by_grouped_remainder": True,
            "negative_trace_package_bound_preserved": True,
            "independent_nonlocal_bulk_step_bound_proved": False,
            "complete_endpoint_B_window_bound_proved": False,
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
        "next_obligation": "Bound the common-regulator joined remainder while retaining sum_(m=1)^620 P_bulk,m-sum_(m=623)^39894 P_bulk,m inside it. Do not promote the negative B trace package to the complete endpoint B sector without an exact estimate of that completion.",
        "proof_boundary": "Exact B trace-package ownership and non-promotion guard only. No independent omitted-step or grouped-remainder bound, complete endpoint B sector, complete Q_K-T or T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified B trace-package scope; nonlocal bulk steps remain grouped", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
