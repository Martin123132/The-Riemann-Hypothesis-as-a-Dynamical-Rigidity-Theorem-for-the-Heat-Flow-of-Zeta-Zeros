#!/usr/bin/env python3
"""Certify exact Weber completion cancellation for matched dual continuations."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_weber_completion_cancellation_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "continuation_partition": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_continuation_stationary_cancellation_gate.json",
}

T = 10_000_000_000
A = 159_577
B = 5_122_421
L = (B - A) // 2
K = (A - 1) // 2


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
    endpoint, mode, x, p = sp.symbols("D n x p", positive=True, real=True)
    alpha = endpoint - 2 * mode
    u = endpoint * p / 2 - mode * x / p
    phase = endpoint**2 * x / 4 + mode**2 * x - endpoint**2 * p**2 / 4 - mode**2 * x**2 / p**2
    require(sp.factor(phase - (alpha**2 * x / 4 - u**2)) == 0, "exact Gaussian phase failed")

    du_dp = sp.diff(u, p)
    expected_du_dp = endpoint / 2 + mode * x / p**2
    require(sp.simplify(du_dp - expected_du_dp) == 0, "Gaussian Jacobian failed")
    require(sp.simplify(u.subs(p, sp.sqrt(x)) - alpha * sp.sqrt(x) / 2) == 0, "tail endpoint failed")

    transformed_amplitude = sp.factor(
        2 * sp.sqrt(x) * (2 + sp.I * sp.pi * endpoint**2 * (x - p**2)) / du_dp
    )
    expected_amplitude = (
        4 * sp.sqrt(x) * p**2
        * (2 + sp.I * sp.pi * endpoint**2 * (x - p**2))
        / (endpoint * p**2 + 2 * mode * x)
    )
    require(sp.simplify(transformed_amplitude - expected_amplitude) == 0, "transformed amplitude failed")

    # Differentiate the convergent Weber formula before taking its Abel
    # boundary value a->i*pi*D^2/4, b->i*pi*n^2*x^2.
    a, b = sp.symbols("a b", positive=True)
    i0 = sp.sqrt(sp.pi) * sp.exp(-2 * sp.sqrt(a * b)) / (2 * sp.sqrt(a))
    i2 = sp.factor(-sp.diff(i0, a))
    expected_i2 = sp.factor(i0 * (1 / (2 * a) + sp.sqrt(b / a)))
    require(sp.simplify(i2 - expected_i2) == 0, "Weber moment derivative failed")

    i0_boundary = sp.exp(-sp.I * sp.pi / 4) * sp.exp(-sp.I * sp.pi * endpoint * mode * x) / endpoint
    i2_boundary = i0_boundary * (2 / (sp.I * sp.pi * endpoint**2) + 2 * mode * x / endpoint)
    bracket_integral = sp.factor(
        (2 + sp.I * sp.pi * endpoint**2 * x) * i0_boundary
        - sp.I * sp.pi * endpoint**2 * i2_boundary
    )
    expected_bracket = sp.I * sp.pi * alpha * endpoint * x * i0_boundary
    require(sp.simplify(bracket_integral - expected_bracket) == 0, "completed bracket failed")

    exponent_identity = sp.expand(endpoint**2 / 4 + mode**2 - endpoint * mode - alpha**2 / 4)
    require(exponent_identity == 0, "completed phase identity failed")

    return {
        "matched_labels": "alpha=A-2k with n_B=L+k and n_A=k for every integer k>=1",
        "face_square_root": "p=sqrt(x-z), so z=x-p^2 and 0<p<sqrt(x)",
        "exact_Gaussian_coordinate": "u_D=Dp/2-nx/p",
        "exact_phase": "D^2*x/4+n^2*x-D^2*p^2/4-n^2*x^2/p^2=alpha^2*x/4-u_D^2",
        "common_tail_endpoint": "u_D(sqrt(x))=alpha*sqrt(x)/2",
        "Jacobian": "du_D/dp=D/2+nx/p^2>0",
        "Gaussian_amplitude": "J_D=4sqrt(x)p_D(u)^2[2+i*pi*D^2(x-p_D(u)^2)]/[D*p_D(u)^2+2nx]",
        "Weber_formula": "I0(a,b)=int_0^infinity exp(-a p^2-b/p^2)dp=sqrt(pi)exp(-2sqrt(ab))/(2sqrt(a))",
        "Weber_second_moment": "I2=-partial_a I0=I0[1/(2a)+sqrt(b/a)]",
        "completed_kernel": "K_D^[0,infinity]=2*i*pi*alpha*x^(3/2)*exp(i*pi*alpha^2*x/4)",
        "completed_pair_cancellation": "K_B^[0,infinity]-K_A^[0,infinity]=0 exactly at fixed alpha",
        "truncated_pair_tail": "K_B^[0,sqrt(x)]-K_A^[0,sqrt(x)]=-e^(i*pi/4+i*pi*alpha^2*x/4)int_(alpha*sqrt(x)/2)^infinity (J_B-J_A)e^(-i*pi*u^2)du",
        "finite_cutoff_identity": "sum_(n=0)^N(B_n-A_n)=sum_(n=0)^L B_n-A_0+sum_(k=1)^(N-L)(B_(L+k)-A_k)-sum_(k=N-L+1)^N A_k for N>L",
        "shift_defect": "the final A block has exactly L terms and must be shown to vanish in the original cutoff/Abel prescription before N->infinity",
    }


def render_note() -> str:
    return f"""# Exact Weber completion of the modular continuation

Date: 2026-08-13

Status: exact Abel-regularized completion and cancellation; not a proof of a
bound for the surviving common Gaussian tails

Fix any matched continuation label

```text
alpha=A-2k,
n_B=L+k, n_A=k, k>=1.                                (WC1)
```

For endpoint `D` and `n=(D-alpha)/2`, put `p=sqrt(x-z)`.  Including the
principal modular factor `exp(i*pi/4)`, its dual endpoint kernel becomes

```text
K_D=exp(i*pi/4) 2sqrt(x) exp(i*pi*x(D^2/4+n^2))
 int_0^sqrt(x) [2+i*pi*D^2(x-p^2)]
 exp(-i*pi[D^2*p^2/4+n^2*x^2/p^2]) dp.              (WC2)
```

The exact monotone coordinate

```text
u_D=Dp/2-nx/p                                         (WC3)
```

satisfies

```text
x(D^2/4+n^2)-D^2*p^2/4-n^2*x^2/p^2
   =alpha^2*x/4-u_D^2,
u_D(0+)=-infinity,
u_D(sqrt(x))=alpha*sqrt(x)/2.                         (WC4)
```

Thus both endpoints have exactly the same Gaussian phase and interval.  No
stationary-phase approximation has been used.

Now complete the `p` integral in (WC2) to infinity in the Abel sense.  For
`Re(a),Re(b)>0`,

```text
I0(a,b)=int_0^infinity exp(-a*p^2-b/p^2)dp
       =sqrt(pi)/(2sqrt(a)) exp(-2sqrt(a*b)),
I2=-partial_a I0=I0[1/(2a)+sqrt(b/a)].                (WC5)
```

Taking the principal boundary values
`a -> i*pi*D^2/4`, `b -> i*pi*n^2*x^2` gives the exact completed kernel

```text
K_D^[0,infinity]
 =2*i*pi*alpha*x^(3/2)*exp(i*pi*alpha^2*x/4).         (WC6)
```

The right side is independent of `D`.  Hence the B and A completed kernels
cancel exactly for every fixed label in (WC1), including negative `alpha`.
Returning to the original finite interval leaves only

```text
K_B-K_A=-exp(i*pi/4+i*pi*alpha^2*x/4)
 int_(alpha*sqrt(x)/2)^infinity
       [J_B(u)-J_A(u)] exp(-i*pi*u^2) du,             (WC7)

J_D=4sqrt(x)p_D(u)^2[2+i*pi*D^2(x-p_D(u)^2)]
       /[D*p_D(u)^2+2nx].                             (WC8)
```

Equation (WC7) is cancellation preserving: one phase, one interval, and one
amplitude difference.  It also explains why the scalar term in the leading
stationary coefficient was not the final continuation residue; the exact
completion includes every stationary correction and cancels the completed
bulk for each fixed pair.

There is a necessary cutoff guard.  With the same dual cutoff `N>L` on both
endpoint sums, exact finite reindexing gives

```text
sum_(n=0)^N (B_n-A_n)
 =sum_(n=0)^L B_n-A_0
  +sum_(k=1)^(N-L)(B_(L+k)-A_k)
  -sum_(k=N-L+1)^N A_k.                              (WC9)
```

The last shift-defect block has exactly `L={L}` terms.  It cannot be dropped
merely because every fixed high-index term tends to zero.  Vanishing of this
block in the original cutoff or Abel prescription is an explicit remaining
obligation before (WC7) may be summed through `k=infinity`.

The exceptional `k=0` term is outside this certificate because the A term
has `n_A=0` and is a face contribution, not an interior Weber saddle.  For
`k>={K + 1}` the lower tail endpoint in (WC7) is negative; those terms are
already included in the common tail series and are not a separate bulk.

Pi provenance: every `pi` in (WC2)-(WC8) comes from the equation-(9)
Fourier/Kummer phase and the principal Jacobi modular factor.  Formula (WC5)
is first used in its absolutely convergent half-plane and only then taken to
the stated Abel boundary.

Proof boundary: the exact completion and completed-pair cancellation for
every fixed `k>=1`, together with the finite-cutoff identity (WC9), are
proved.  Infinite-cutoff reindexing, vanishing of the shift defect, a bound
for the resulting common-tail series, the exceptional `k=0` face/fold
splice, the complete source-minus-target estimate, `T_upper`, `Lambda<=0`,
PF-infinity, RH, and any prize-level conclusion remain open.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    decision = dependencies["continuation_partition"]["decision"]
    require(decision["matched_positive_continuation_pairs"] == K, "continuation count drift")
    require(decision["all_matched_positive_continuation_reduced_phases_are_x_nonstationary"] is True, "x-phase dependency drift")
    require(decision["exceptional_k_0_face_fold_term_cancelled"] is False, "k=0 dependency overclaim")

    artifact = {
        "kind": STEM,
        "status": "each_fixed_matched_dual_continuation_Weber_bulk_cancels_exactly_finite_cutoff_shift_defect_open",
        "passed": True,
        "scope": {
            "height": T,
            "lower_endpoint": A,
            "upper_endpoint": B,
            "B_roster_offset": L,
            "matched_positive_continuation_pairs": K,
            "fixed_matched_continuation_labels": "every integer k>=1",
        },
        "symbolic_certificate": symbolic_certificate(),
        "decision": {
            "matched_endpoint_dual_phases_share_exact_Gaussian_coordinate": True,
            "matched_endpoint_dual_domains_share_exact_Gaussian_interval": True,
            "Abel_completed_Weber_kernel_evaluated_exactly": True,
            "completed_B_minus_A_kernel_cancels_for_each_fixed_k_at_least_1": True,
            "every_finite_matched_prefix_reduced_to_common_Gaussian_tails": True,
            "same_cutoff_reindexing_has_L_term_A_shift_defect": True,
            "infinite_continuation_Abel_reindexing_justified": False,
            "shift_defect_vanishing_proved": False,
            "stationary_expansion_remainder_needed_for_completed_bulk": False,
            "exceptional_k_0_face_fold_term_cancelled": False,
            "common_Gaussian_tail_sum_bounded": False,
            "complete_modular_continuation_bounded": False,
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
        "next_obligation": "First prove that the L-term high-index A shift defect in WC9 vanishes in the original dual cutoff or Abel prescription, with a bound uniform enough to pass through the x,z integral. Only then form the all-k common Gaussian tail series, retain the x=0 endpoint grouping, and splice k=0 to the face/fold atlas.",
        "proof_boundary": "Exact Abel-Weber completion for each fixed matched k and exact finite-cutoff bookkeeping only. Infinite reindexing, shift-defect vanishing, the summed common-tail bound, k=0 fold splice, complete source-minus-target estimate, T_upper, Lambda<=0, PF-infinity, RH, and any prize-level conclusion remain open.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(), encoding="utf-8")
    print("certified exact Weber completion cancellation for matched dual continuations", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
