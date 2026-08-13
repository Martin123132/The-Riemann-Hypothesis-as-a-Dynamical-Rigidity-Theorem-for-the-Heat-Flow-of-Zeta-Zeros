#!/usr/bin/env python3
"""Certify the leading dual-continuation cancellation and saddle partition."""

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

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_continuation_stationary_cancellation_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "Abel_theta_roster": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_modular_dual_roster_gate.json",
}

PRECISION = 100
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
    x, alpha, gap = sp.symbols("x alpha gap", positive=True, real=True)
    endpoint = alpha + gap
    mode = gap / 2
    z = sp.symbols("z", positive=True, real=True)

    r = 1 / z - 1 / x
    phase = (
        sp.pi * endpoint**2 * z / 4
        - sp.pi * mode**2 / r
        + sp.symbols("t", positive=True, real=True) * sp.log((1 - x) / x) / 2
    )
    z_star = x * alpha / endpoint
    require(sp.simplify(sp.diff(phase, z).subs(z, z_star)) == 0, "dual z saddle failed")

    reduced_phase = sp.factor(sp.simplify(phase.subs(z, z_star)))
    expected_reduced = sp.pi * alpha**2 * x / 4 + sp.symbols("t", positive=True, real=True) * sp.log((1 - x) / x) / 2
    require(sp.simplify(reduced_phase - expected_reduced) == 0, "reduced phase failed")

    hessian = sp.factor(sp.simplify(sp.diff(phase, z, 2).subs(z, z_star)))
    expected_hessian = -sp.pi * endpoint**3 / (2 * gap * x)
    require(sp.simplify(hessian - expected_hessian) == 0, "dual z Hessian failed")

    r_star = gap / (x * alpha)
    hessian_abs = -expected_hessian
    driver = z_star ** (-sp.Rational(1, 2)) * (2 + sp.I * sp.pi * endpoint**2 * z_star)
    stationary_multiplier = (
        r_star ** (-sp.Rational(1, 2))
        * driver
        * sp.sqrt(2 * sp.pi / hessian_abs)
    )
    expected_multiplier = 4 * sp.sqrt(x) / endpoint + 2 * sp.I * sp.pi * alpha * x ** sp.Rational(3, 2)
    replay = sp.powsimp(sp.powdenest(stationary_multiplier, force=True), force=True)
    require(sp.simplify(replay - expected_multiplier) == 0, "stationary multiplier failed")

    d1, d2 = sp.symbols("D_1 D_2", positive=True, real=True)
    multiplier_1 = 4 * sp.sqrt(x) / d1 + 2 * sp.I * sp.pi * alpha * x ** sp.Rational(3, 2)
    multiplier_2 = 4 * sp.sqrt(x) / d2 + 2 * sp.I * sp.pi * alpha * x ** sp.Rational(3, 2)
    require(
        sp.simplify(multiplier_1 - multiplier_2 - 4 * sp.sqrt(x) * (1 / d1 - 1 / d2)) == 0,
        "paired leading-channel cancellation failed",
    )

    return {
        "continuation_label": "alpha=A-2k=alpha_(B,L+k)=alpha_(A,k)",
        "dual_z_saddle": "z_*(D,alpha)=x*alpha/D for 0<alpha<D",
        "reduced_phase": "psi_alpha(x)=pi*alpha^2*x/4+(t/2)log((1-x)/x)",
        "z_Hessian": "Psi_zz(z_*)=-pi*D^3/[2(D-alpha)x]",
        "modular_coordinate": "r_*=(D-alpha)/(x*alpha)",
        "endpoint_driver": "h_D(z)=z^(-1/2)(2+i*pi*D^2*z)",
        "leading_stationary_multiplier": "M_D(alpha,x)=4*sqrt(x)/D+2*i*pi*alpha*x^(3/2)",
        "paired_difference": "M_B-M_A=4*sqrt(x)(1/B-1/A)",
        "cancelled_channel": "the universal 2*i*pi*alpha*x^(3/2) channel cancels exactly",
        "outer_amplitude_residual": "a_t(x)(M_B-M_A)=4(1/B-1/A)x^(-5/4)(1-x)^(-1/4)",
        "x_phase_derivative": "psi_alpha'(x)=pi*alpha^2/4-t/[2x(1-x)]",
    }


def interval_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    t = arb(T)
    pi = arb.pi()
    threshold = (8 * t / pi).sqrt()
    lower_label = arb(A - 2)
    upper_label = arb(A)
    require(lower_label < threshold < upper_label, "continuation threshold separation failed")

    tangential_margin = 2 * t - pi * lower_label**2 / 4
    require(tangential_margin > arb("479304.70"), "tangential derivative margin drift")

    scalar_coefficient = 4 * (arb(1) / B - arb(1) / A)
    require(scalar_coefficient < 0, "scalar residual sign drift")

    require(B - 2 * L == A, "B roster endpoint drift")
    require(A - 2 * K == 1, "last positive continuation label drift")
    require(A - 2 * (K + 1) == -1, "first negative continuation label drift")

    return {
        "height": T,
        "lower_endpoint": A,
        "upper_endpoint": B,
        "B_roster_offset": L,
        "positive_continuation_pairs": K,
        "threshold_sqrt_8t_over_pi_ball": threshold.str(PRECISION, more=True),
        "A_minus_threshold_ball": (upper_label - threshold).str(PRECISION, more=True),
        "threshold_minus_A_minus_2_ball": (threshold - lower_label).str(PRECISION, more=True),
        "uniform_x_derivative_magnitude_lower_ball": tangential_margin.str(PRECISION, more=True),
        "paired_scalar_coefficient_ball": scalar_coefficient.str(PRECISION, more=True),
        "paired_scalar_coefficient_exact": f"4*({A}-{B})/({A}*{B})",
        "continuation_partition": {
            "k_0": "B has the alpha=A interior z saddle while A has its n=0 face term; no cancellation is asserted",
            "k_1_through_K": f"{K} matched positive-alpha interior-z saddle pairs; all are x-nonstationary",
            "k_at_least_K_plus_1": f"k>={K + 1} gives alpha<=-1, so the formal z saddle lies outside 0<z<x",
        },
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    return f"""# Modular-dual continuation stationary cancellation

Date: 2026-08-13

Status: exact leading-stationary cancellation and complete saddle-roster
partition; not a proof of a full modular-continuation bound

For a positive dual index, write

```text
alpha=D-2n,  z_*=x alpha/D,  r_*=(D-alpha)/(x alpha). (MC1)
```

At the `z` saddle the dual phase and Hessian are exactly

```text
psi_alpha(x)=pi alpha^2 x/4+(t/2)log((1-x)/x),
Psi_zz(z_*)=-pi D^3/[2(D-alpha)x].                    (MC2)
```

Include the modular factor `r^(-1/2)`, the endpoint driver
`z^(-1/2)(2+i*pi*D^2*z)`, and the leading stationary scale
`sqrt(2*pi/|Psi_zz|)`.  Their product simplifies without an asymptotic
approximation to

```text
M_D(alpha,x)=4 sqrt(x)/D+2 i*pi*alpha*x^(3/2).        (MC3)
```

For each matched continuation pair

```text
alpha_(B,L+k)=alpha_(A,k)=A-2k,                      (MC4)
```

the common stationary phase, Hessian signature, dual `+/-n` multiplicity,
and universal second term in (MC3) agree.  With the compulsory endpoint
signs, the leading difference is therefore exactly

```text
M_B-M_A=4 sqrt(x)(1/B-1/A).                          (MC5)
```

This cancels the larger `2 i*pi*alpha*x^(3/2)` channel.  It does not cancel
the scalar endpoint channel, the stationary-expansion remainder, or any
face and corner boundary term.

The arithmetic threshold is rigorously separated by Arb:

```text
A-2 < sqrt(8t/pi) < A,
sqrt(8t/pi)={c['threshold_sqrt_8t_over_pi_ball']},
A-sqrt(8t/pi)={c['A_minus_threshold_ball']},
sqrt(8t/pi)-(A-2)={c['threshold_minus_A_minus_2_ball']}. (MC6)
```

Consequently the dual continuation has the complete partition

```text
k=0: exceptional A-face/B-interior fold chart;
1<=k<={K}: matched interior-z saddles, but no real x saddle;
k>={K + 1}: alpha<=-1, hence no interior z saddle.             (MC7)
```

For the middle block,

```text
|psi_alpha'(x)| >= 2t-pi(A-2)^2/4
                      > {c['uniform_x_derivative_magnitude_lower_ball']}
```

throughout `0<x<=1/2`.  Thus its reduced leading phase is uniformly
nonstationary.  The surviving scalar coefficient in (MC5) is
`{c['paired_scalar_coefficient_exact']}` with Arb ball
`{c['paired_scalar_coefficient_ball']}`.

Pi provenance: `pi` is inherited from the equation-(9) Fourier/Kummer phase
and the canonical Jacobi transformation.  The checker independently bounds
it from Machin's identity using exact rational alternating sums.

Proof boundary: the identity (MC5), the threshold separation, and the
continuation saddle partition are proved only for the leading interior
`z`-stationary contribution.  No uniform stationary-expansion remainder,
summed scalar-channel estimate, exceptional `k=0` fold splice, complete
modular continuation bound, source-minus-target estimate, `T_upper`,
`Lambda<=0`, PF-infinity, RH, and any prize-level conclusion remain open.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    roster = dependencies["Abel_theta_roster"]["decision"]
    require(roster["Jacobi_modular_dual_phase_identified"] is True, "modular-phase dependency drift")
    require(roster["B_dual_indices_recover_complete_odd_source_roster"] is True, "dual-roster dependency drift")
    require(roster["continuation_endpoints_cancel_termwise"] is False, "dependency overclaims cancellation")

    artifact = {
        "kind": STEM,
        "status": "leading_dual_continuation_stationary_channel_cancels_and_remaining_positive_labels_are_x_nonstationary",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(),
        "decision": {
            "dual_continuation_saddle_roster_partition_complete": True,
            "matched_positive_continuation_pairs": K,
            "leading_universal_z_stationary_channel_cancels_for_k_1_through_K": True,
            "scalar_endpoint_channel_survives": True,
            "all_matched_positive_continuation_reduced_phases_are_x_nonstationary": True,
            "all_k_at_least_K_plus_1_are_z_nonstationary": True,
            "exceptional_k_0_face_fold_term_cancelled": False,
            "uniform_z_stationary_expansion_remainder_bounded": False,
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
            "flint_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Derive a uniform paired z-stationary expansion for k=1..K with the B-minus-A driver grouped, bound its remainder before summing, integrate the surviving scalar channel in x using the certified derivative gap, treat k=0 in the existing face/fold chart, and bound k>=K+1 by z-nonstationary integration.",
        "proof_boundary": "Leading z-stationary cancellation and exact continuation saddle classification only. No uniform stationary remainder, summed continuation bound, k=0 fold splice, complete source-minus-target estimate, T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified modular-dual continuation leading cancellation and saddle partition", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
