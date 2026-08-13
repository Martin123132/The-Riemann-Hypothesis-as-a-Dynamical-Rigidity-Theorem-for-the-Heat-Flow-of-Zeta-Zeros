#!/usr/bin/env python3
"""Certify source-roster recovery and same-index theta-tail reassembly."""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_weber_source_roster_tail_reassembly_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "Weber_completion": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_weber_completion_cancellation_gate.json",
    "shift_defect": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_cutoff_shift_defect_vanishing_gate.json",
    "Abel_theta_roster": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_modular_dual_roster_gate.json",
    "symmetric_Poisson": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate.json",
}

T = 10_000_000_000
A = 159_577
B = 5_122_421
L = (B - A) // 2


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


def symbolic_certificate() -> dict[str, str]:
    x, alpha = sp.symbols("x alpha", positive=True, real=True)
    outer_amplitude = x ** (-sp.Rational(7, 4)) * (1 - x) ** (-sp.Rational(1, 4))
    completed_kernel = 2 * sp.I * sp.pi * alpha * x ** sp.Rational(3, 2)
    normalized = sp.simplify((2 / (4 * sp.I * sp.pi)) * outer_amplitude * completed_kernel)
    expected = alpha * x ** (-sp.Rational(1, 4)) * (1 - x) ** (-sp.Rational(1, 4))
    require(sp.simplify(normalized - expected) == 0, "source normalization failed")

    endpoint, mode, s = sp.symbols("D n s", positive=True, real=True)
    p = sp.sqrt(x * s)
    tail_phase = (
        endpoint**2 * x / 4
        + mode**2 * x
        - endpoint**2 * p**2 / 4
        - mode**2 * x**2 / p**2
    )
    expected_tail_phase = x * (endpoint**2 * (1 - s) / 4 + mode**2 * (1 - 1 / s))
    require(sp.simplify(tail_phase - expected_tail_phase) == 0, "tail phase failed")

    driver = sp.exp(sp.I * sp.pi * endpoint**2 * x * (1 - s) / 4) * (
        2 + sp.I * sp.pi * endpoint**2 * x * (1 - s)
    )
    require(sp.simplify(driver.subs(s, 1) - 2) == 0, "tail driver face value failed")
    d_driver = sp.simplify(sp.diff(driver, s).subs(s, 1))
    require(d_driver == -3 * sp.I * sp.pi * endpoint**2 * x / 2, "tail driver face derivative failed")

    require(B - 2 == A + 2 * (L - 1), "completed source roster endpoints failed")
    return {
        "positive_dual_identity": "sum_(n>=1)(B_n-A_n)=sum_(n=1)^L B_n+sum_(k>=1)(B_(L+k)-A_k)",
        "completed_finite_B_block": "2/(4i*pi) times K_B^[0,infinity](alpha) times a_t(x)=alpha[x(1-x)]^(-1/4)e^(i*pi*alpha^2*x/4)",
        "recovered_source_labels": "alpha=B-2,B-4,...,A (L labels); the upper label B remains in the zero/direct exceptional sector",
        "same_index_tail_identity": "sum_(n>=1)(B_n-A_n)=sum_(n=1)^L Full_B,n+sum_(n>=1)(Tail_A,n-Tail_B,n)",
        "tail_coordinate": "s=p^2/x>=1, tau_x(s)=x(1-1/s)",
        "tail_driver": "g_D(x,s)=exp[i*pi*D^2*x(1-s)/4]{2+i*pi*D^2*x(1-s)}",
        "summed_tail": "2sum_(n>=1)(Tail_A,n-Tail_B,n)=e^(i*pi/4)x int_1^infinity s^(-1/2)[g_A-g_B][Theta(tau_x(s))-1]ds",
        "tail_face_cancellation": "g_A(x,1)-g_B(x,1)=0 and partial_s(g_A-g_B)(x,1)=(3i*pi/2)(B^2-A^2)x",
        "exceptional_sector": "dual n=0 endpoint difference and the direct -1 in Theta(tau)-1 stay grouped with the missing upper source label B",
    }


def render_note() -> str:
    return f"""# Weber source-roster and same-index theta-tail reassembly

Date: 2026-08-13

Status: exact fixed-height reassembly and source normalization; not a proof
of a bound for the surviving theta tail or zero-mode sector

Write `B_n` and `A_n` for the positive modular-dual endpoint kernels and
`Full_(D,n)`, `Tail_(D,n)` for their Weber completion and `p>=sqrt(x)`
tail.  The licensed cutoff shift gives

```text
sum_(n>=1)(B_n-A_n)
 =sum_(n=1)^L B_n+sum_(k>=1)(B_(L+k)-A_k).           (WR1)
```

Every fixed continuation completion cancels, while
`B_n=Full_(B,n)-Tail_(B,n)`.  Therefore (WR1) becomes exactly

```text
sum_(n>=1)(B_n-A_n)
 =sum_(n=1)^L Full_(B,n)
  +sum_(n>=1)[Tail_(A,n)-Tail_(B,n)].                (WR2)
```

For `alpha=B-2n`, the completed kernel is

```text
Full_(B,n)=2*i*pi*alpha*x^(3/2)
              exp(i*pi*alpha^2*x/4).                 (WR3)
```

The Abel-triangle prefactor is `1/(4i*pi)`, positive dual modes have
multiplicity two, and `a_t(x)=x^(-7/4)(1-x)^(-1/4)` before the common
logistic phase.  Hence

```text
[2/(4i*pi)] a_t(x) Full_(B,n)
 =alpha[x(1-x)]^(-1/4)exp(i*pi*alpha^2*x/4),         (WR4)
```

which is exactly the equation-(9) source integrand.  Thus the finite block
in (WR2) recovers the `L={L}` labels

```text
B-2,B-4,...,A.                                       (WR5)
```

There is no normalization or phase discrepancy.  The missing upper label
`B`, the dual `n=0` endpoint difference, and the direct `-1` in
`Theta(tau)-1` remain in one exceptional sector; this gate does not split
them.

The remainder in (WR2) can be summed at equal dual index.  Put

```text
s=p^2/x>=1,
tau_x(s)=x(1-1/s),
g_D(x,s)=exp[i*pi*D^2*x(1-s)/4]
           {{2+i*pi*D^2*x(1-s)}}.                   (WR6)
```

Then, in the inherited cutoff/Abel sense,

```text
2sum_(n>=1)[Tail_(A,n)-Tail_(B,n)]
 =e^(i*pi/4)x int_1^infinity s^(-1/2)
   [g_A(x,s)-g_B(x,s)]
   [Theta(tau_x(s))-1] ds.                           (WR7)
```

This restores the endpoint difference before any absolute value.  At the
new face,

```text
g_A(x,1)-g_B(x,1)=0,
partial_s(g_A-g_B)(x,1)
  =(3i*pi/2)(B^2-A^2)x.                              (WR8)
```

So the apparent termwise tail boundary is absent from the summed object.
The remaining challenge is quantitative: control (WR7) jointly in `s` and
`x` and close the exceptional zero/direct sector against the existing
face/fold atlas.

Pi provenance: the `pi` factors are inherited from equation (9), the Jacobi
modular factor, and the Weber boundary value.  Equation (WR4) explicitly
shows their cancellation into the original source normalization.

Proof boundary: exact source recovery for labels (WR5), same-index tail
reassembly, and face cancellation (WR8) only.  No estimate for (WR7), no
identification of the exceptional sector with the missing `B` source term,
no complete source-minus-target estimate, `T_upper`, height-uniform theorem,
`Lambda<=0`, PF-infinity, RH, and any prize-level conclusion remain open.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    weber = dependencies["Weber_completion"]["decision"]
    require(weber["completed_B_minus_A_kernel_cancels_for_each_fixed_k_at_least_1"] is True, "Weber cancellation drift")
    shift = dependencies["shift_defect"]["decision"]
    require(shift["fixed_height_symmetric_dual_cutoff_reindexing_justified"] is True, "cutoff reindexing drift")
    roster = dependencies["Abel_theta_roster"]["decision"]
    require(roster["B_dual_indices_recover_complete_odd_source_roster"] is True, "dual roster drift")
    poisson = dependencies["symmetric_Poisson"]["symbolic_reduction"]
    require(
        poisson["alpha_roster_integrand"] == "f_x(u)=(A+2u) exp(i*pi*x*(A+2u)^2/4)",
        "source normalization dependency drift",
    )

    artifact = {
        "kind": STEM,
        "status": "completed_B_block_recovers_source_roster_except_B_remaining_positive_dual_sector_is_one_grouped_theta_tail",
        "passed": True,
        "scope": {
            "height": T,
            "lower_endpoint": A,
            "upper_endpoint": B,
            "roster_length": L,
            "recovered_completed_labels": L,
        },
        "symbolic_certificate": symbolic_certificate(),
        "decision": {
            "positive_dual_sector_reindexed_with_cutoff_guard": True,
            "completed_finite_B_block_has_exact_equation9_source_normalization": True,
            "completed_B_block_recovers_labels_B_minus_2_through_A": True,
            "upper_source_label_B_recovered_by_this_gate": False,
            "all_remaining_positive_dual_tails_reassembled_at_same_index": True,
            "summed_tail_retains_endpoint_difference_before_absolute_values": True,
            "summed_tail_driver_vanishes_at_s_equal_1": True,
            "summed_theta_tail_quantitatively_bounded": False,
            "zero_mode_direct_exceptional_sector_closed": False,
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
            "process_priority": priority,
        },
        "next_obligation": "Analyze WR7 as one endpoint-grouped theta tail. Use the linear s=1 zero before majorization, seek a second modular or contour representation for s>1, and isolate the x=0/s=1 corner. In parallel, combine the dual n=0 and direct -1 terms without separating endpoints and test whether they supply the missing upper source label B plus the certified face/fold currents.",
        "proof_boundary": "Exact fixed-height source-roster recovery except B and same-index theta-tail reassembly only. No tail bound, exceptional-sector closure, complete source-minus-target estimate, T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(), encoding="utf-8")
    print("certified Weber source-roster recovery and same-index theta-tail reassembly", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
