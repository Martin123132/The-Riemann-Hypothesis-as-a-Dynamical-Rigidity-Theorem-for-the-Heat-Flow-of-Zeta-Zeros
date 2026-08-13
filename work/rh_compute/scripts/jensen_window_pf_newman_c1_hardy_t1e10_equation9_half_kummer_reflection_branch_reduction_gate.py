#!/usr/bin/env python3
"""Certify exact half-Kummer reflection and the one-branch Poisson saddle ledger."""

from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path
from typing import Any

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_kummer_reflection_branch_reduction_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "portcullis": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_finite_poisson_portcullis_saddle_reduction_gate.json",
    "symmetric_poisson": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate.json",
    "symmetric_endpoint_residual": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_symmetric_endpoint_tail_reassembly_gate.json",
    "coverage_ledger": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_ordinary_mode_coverage_ledger_gate.json",
}

T = 10_000_000_000
A = 159_577
B = 5_122_421
TARGET_START = 622
TARGET_END = 39_894


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
    alpha, x, t, m = sp.symbols("alpha x t m", positive=True, real=True)
    phase = sp.pi * alpha**2 * x / 4 + t * sp.log((1 - x) / x) / 2
    # On 0<x<1, log(x/(1-x))=-log((1-x)/x) on the real branch.
    reflected_phase = sp.pi * alpha**2 * (1 - x) / 4 - t * sp.log((1 - x) / x) / 2
    conjugate_phase = -phase
    require(
        sp.simplify(reflected_phase - conjugate_phase - sp.pi * alpha**2 / 4) == 0,
        "Kummer reflection phase failed",
    )

    k = sp.symbols("k", integer=True)
    odd_alpha = 2 * k + 1
    odd_square_character = sp.expand(odd_alpha**2 - 1) / 8
    require(sp.simplify(odd_square_character - k * (k + 1) / 2) == 0, "odd-square residue failed")
    require(all((2 * value + 1) ** 2 % 8 == 1 for value in range(-20, 21)), "odd-square residue witnesses failed")

    x_m = 2 * sp.pi * m**2 / (t + 2 * sp.pi * m**2)
    nt = sp.sqrt(t / (2 * sp.pi))
    require(sp.simplify(x_m.subs(m, nt) - sp.Rational(1, 2)) == 0, "half-line saddle threshold failed")
    x_derivative = sp.factor(sp.diff(x_m, m))
    expected_x_derivative = 4 * sp.pi * m * t / (2 * sp.pi * m**2 + t) ** 2
    require(sp.simplify(x_derivative - expected_x_derivative) == 0, "x saddle derivative failed")

    alpha_m = 2 * m + t / (sp.pi * m)
    require(sp.simplify(alpha_m.subs(m, t / (2 * sp.pi * m)) - alpha_m) == 0, "continuous reflection failed")
    require(sp.simplify(x_m.subs(m, t / (2 * sp.pi * m)) - (1 - x_m)) == 0, "x reflection failed")

    alpha_derivative = sp.diff(
        sp.pi * alpha**2 * x / 4 - sp.pi * m * alpha + t * sp.log((1 - x) / x) / 2,
        alpha,
    )
    require(sp.simplify(alpha_derivative - sp.pi * (alpha * x / 2 - m)) == 0, "alpha derivative failed")

    return {
        "source_integrand": "K_(alpha,t;x)=alpha exp(i*pi*alpha^2*x/4) exp(i*t*log((1-x)/x)/2)/[x(1-x)]^(1/4)",
        "odd_reflection": "K_(alpha,t;1-x)=exp(i*pi/4) conjugate(K_(alpha,t;x)) for every odd integer alpha",
        "half_integral_identity": "exp(-i*pi/8) integral_0^1 K dx=2 Re[exp(-i*pi/8) integral_0^(1/2) K dx]",
        "finite_roster_identity": "The same identity holds after the finite odd-alpha sum, including its two endpoint half-weights.",
        "half_poisson_formula": "Q_half=integral_0^(1/2) W_t(x){H_x+lim_(M->infinity)sum_(m=-M)^M I_m(x)}dx",
        "joint_saddle": "x_m=2*pi*m^2/(t+2*pi*m^2), alpha_m=2m+t/(pi*m)",
        "half_domain_threshold": "x_m<1/2 iff 0<m<sqrt(t/(2*pi)); x_m=1/2 at m=sqrt(t/(2*pi))",
        "continuous_reflection": "m*=t/(2*pi*m) gives alpha_(m*)=alpha_m and x_(m*)=1-x_m",
        "branch_consequence": "The m>sqrt(t/(2*pi)) stationary branch lies in the reflected x>1/2 half and is already supplied by exact conjugation; it is not a second stationary family to estimate in Q_half.",
        "nonpositive_guard": "For m<=0 and alpha>0, x>0, partial_alpha Phi=pi(alpha*x/2-m)>0, so no joint saddle exists.",
        "pi_provenance": "pi comes from the equation-(9) Kummer phase, odd-square character, and Fourier-Poisson mode phase.",
    }


def finite_height_ledger() -> dict[str, Any]:
    import mpmath as mp

    mp.mp.dps = 80
    pi = mp.pi
    nt = mp.sqrt(T / (2 * pi))
    x = lambda mode: 2 * pi * mode**2 / (T + 2 * pi * mode**2)
    alpha = lambda mode: 2 * mode + T / (pi * mode)
    require(mp.mpf(TARGET_END) < nt < mp.mpf(TARGET_END + 1), "N_t interval drift")
    require(x(TARGET_END) < mp.mpf("0.5") < x(TARGET_END + 1), "half-domain mode split drift")
    require(alpha(621) > B and alpha(622) < B, "B crossing drift")
    require(alpha(39852) > A and alpha(39853) < A, "A crossing drift")
    return {
        "height": T,
        "source_endpoints": [A, B],
        "N_t": mp.nstr(nt, 80),
        "N_t_interval": [TARGET_END, TARGET_END + 1],
        "x_39894": mp.nstr(x(TARGET_END), 80),
        "x_39895": mp.nstr(x(TARGET_END + 1), 80),
        "half_domain_stationary_partition": {
            "nonpositive_modes": "no alpha stationary point",
            "low_B_endpoint_side": [1, 621],
            "lower_interior": [622, 39852],
            "lower_A_endpoint_side": [39853, 39894],
            "half_boundary_and_outer_positive_completion": "m>=39895 has no interior x saddle on 0<x<=1/2",
        },
        "source_target": [TARGET_START, TARGET_END],
        "source_target_count": TARGET_END - TARGET_START + 1,
        "removed_full_domain_stationary_ranges": {
            "reflected_A_transition": [39895, 39936],
            "upper_interior": [39937, 2560588],
            "interpretation": "These are stationary only in x>1/2 and are represented exactly by conjugation of the half integral.",
        },
        "remaining_half_domain_transitions": [
            "B endpoint crossing 621|622",
            "A endpoint/half-boundary crossing 39894|39895",
        ],
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["finite_height_ledger"]
    return f"""# Exact half-Kummer reflection and one-branch reduction

Date: 2026-08-13

Status: exact-lemma certificate; not a proof of the quantitative half-domain remainder

For an odd source index `alpha`, write the equation-(9) integrand as

```text
K_alpha(x)=alpha exp(i*pi*alpha^2*x/4)
 exp[i(t/2)log((1-x)/x)]/[x(1-x)]^(1/4).
```

Odd squares satisfy `alpha^2=1 (mod 8)`.  Therefore

```text
K_alpha(1-x)=exp(i*pi/4) conjugate(K_alpha(x)),         (HR1)

exp(-i*pi/8) integral_0^1 K_alpha(x)dx
 =2 Re[exp(-i*pi/8) integral_0^(1/2)K_alpha(x)dx].     (HR2)
```

Equation (HR2) is termwise exact, so it remains exact after the finite odd
source-roster sum and its endpoint half-weights.  It also survives symmetric
finite Poisson summation and the already certified interchange:

```text
Q_main=2 Re[exp(-i*pi/8) Q_half],

Q_half=integral_0^(1/2) W_t(x)
 {{H_x+lim_(M->infinity)sum_(m=-M)^M I_m(x)}}dx.       (HR3)
```

This changes the stationary ledger decisively.  A positive Poisson mode has

```text
x_m=2*pi*m^2/(t+2*pi*m^2),
alpha_m=2m+t/(pi*m).                                   (HR4)
```

Since `x_m` is strictly increasing,

```text
x_m<1/2 iff m<sqrt(t/(2*pi)).                          (HR5)
```

At `t=10^10`,

```text
sqrt(t/(2*pi))={c['N_t']},
39894<sqrt(t/(2*pi))<39895.                            (HR6)
```

Consequently the half-integral has only one positive stationary branch:

```text
1..621:          B-endpoint side (alpha_m>B),
622..39852:      lower interior,
39853..39894:    A-endpoint side (alpha_m<A),
m>=39895:        no interior x saddle on 0<x<=1/2.     (HR7)
```

For `m<=0`, `partial_alpha Phi=pi(alpha*x/2-m)>0`, so no joint
saddle exists either.  The formerly listed reflected transition
`39895..39936` and upper interior `39937..2560588` live in `x>1/2`.
They are already represented exactly by the conjugate in (HR2), and must not
be counted as a second stationary source family in (HR3).

The continuous involution explains the geometry:

```text
m*=t/(2*pi*m),  alpha_(m*)=alpha_m,  x_(m*)=1-x_m.     (HR8)
```

It need not preserve integers because (HR1)--(HR3), not rounded reciprocal
pairing, performs the exact discrete reassembly.

The quantitative problem is now smaller.  On the half-domain there are only
two endpoint transitions: the `B` crossing `621|622` and the combined `A`
endpoint/`x=1/2` crossing `39894|39895`.  The latter is the characteristic
fold already represented by the certified atlas.  The former is
noncharacteristic and should be treated as one bulk-plus-Fresnel crossing,
together with the nonpositive and outer-positive nonstationary completion.

Pi provenance: every `pi` in (HR1)--(HR8) is inherited from the equation-(9)
Kummer phase, the odd-square character, or the Fourier-Poisson phase.  No
geometric fit is introduced.

Proof boundary: exact source reflection, half-domain Poisson representation,
and stationary classification only.  No half-domain nonstationary bound,
`B` crossing theorem, quantitative `A` splice, complete `T_upper`,
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["symmetric_poisson"]["decision"]["symmetric_poisson_interchange_proved"] is True, "interchange dependency drift")
    require(dependencies["symmetric_endpoint_residual"]["decision"]["symmetric_residual_limit_proved"] is True, "residual dependency drift")
    coverage = dependencies["coverage_ledger"]["mode_coverage"]["source_classical_target"]
    require(coverage == {"range": [TARGET_START, TARGET_END], "count": 39_273}, "coverage drift")

    artifact = {
        "kind": STEM,
        "status": "exact_half_kummer_reflection_and_single_stationary_branch_reduction_proved",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "finite_height_ledger": finite_height_ledger(),
        "decision": {
            "odd_alpha_half_kummer_reflection_proved": True,
            "full_source_main_equals_twice_real_half_integral": True,
            "symmetric_poisson_may_be_restricted_to_half_domain": True,
            "half_domain_has_one_positive_stationary_branch": True,
            "reflected_upper_branch_requires_separate_stationary_estimate": False,
            "rounded_reciprocal_integer_pairing_used": False,
            "half_domain_nonstationary_completion_bounded": False,
            "B_endpoint_crossing_bounded": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)} for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "sympy_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Derive the noncharacteristic B bulk-plus-Fresnel crossing at 621|622 on the half-domain, then bound the nonpositive and m>=39895 completion jointly while retaining the A/half-boundary fold atlas.",
        "proof_boundary": "Exact odd-alpha half-Kummer reflection, half-domain Poisson representation, and one-branch stationary classification only. No quantitative half-domain complement, B crossing theorem, A splice, complete T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified exact half-Kummer reflection: one positive stationary branch on x<=1/2", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
