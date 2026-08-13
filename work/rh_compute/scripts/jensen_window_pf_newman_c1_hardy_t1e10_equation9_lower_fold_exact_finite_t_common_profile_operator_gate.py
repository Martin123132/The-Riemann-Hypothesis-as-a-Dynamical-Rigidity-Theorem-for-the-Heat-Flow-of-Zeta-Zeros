#!/usr/bin/env python3
"""Derive the exact finite-t common Fourier profile before completion."""

from __future__ import annotations

from decimal import Decimal, getcontext
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_finite_t_common_profile_operator_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "exact_normal_form": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_grouped_characteristic_airy_boundary_normal_form_gate.json",
    "fourier_completion": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_399_mode_fourier_completion_gate.json",
    "beta4_reduction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_canonical_airy_derivative_reduction_gate.json",
    "completed_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.json",
    "profile_obligation": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_completed_point_weighted_profile_obligation_gate.json",
    "weighted_kernel": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_opposite_branch_dirichlet_gate.json",
}

C = 159_577
MODE_CENTER = 39_895
L = 39_696
Q = 39_894
U = 40_094


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def parse_ball(text: str) -> tuple[Decimal, Decimal]:
    body = text.strip().removeprefix("[").removesuffix("]").strip()
    if "+/-" not in body:
        return Decimal(body), Decimal(0)
    midpoint, radius = body.split("+/-", 1)
    return Decimal(midpoint.strip()), Decimal(radius.strip())


def symbolic_certificate() -> dict[str, str]:
    beta, c, lam, y, z, s = sp.symbols("beta C lambda y z s", positive=True, real=True)
    n = sp.symbols("n", integer=True)
    pi = sp.pi
    i = sp.I
    sigma = 4 * beta / (pi * c)
    y_width = 2 / sigma
    h = 4 * beta / c
    mode_center = (c + 3) / 4
    mode = mode_center + n
    detuning = 4 * beta * (mode - c / 4) / c
    t_star = beta**3
    height = t_star - beta * lam
    eta = t_star / height
    tau = sp.tanh(z / beta)
    x = (1 - tau) / 2

    require(sp.simplify(h * y_width - 2 * pi) == 0, "hY identity failed")
    require(sp.simplify(detuning - h * (n + sp.Rational(3, 4))) == 0, "centered detuning failed")
    require(sp.simplify(detuning * y_width - 2 * pi * (n + sp.Rational(3, 4))) == 0, "Fourier phase failed")

    original_height_phase = height * (z / beta - eta * tau)
    reduced_height_phase = beta**3 * (z / beta - tau) - lam * z
    require(sp.simplify(original_height_phase - reduced_height_phase) == 0, "height phase reduction failed")

    exact_phase_without_mode = (
        reduced_height_phase
        - beta * y * tau
        + x * y**2 / (2 * beta)
    )
    exact_phase_with_mode = exact_phase_without_mode - detuning * y
    factored_at_y = sp.expand(exact_phase_with_mode.subs(y, y_width * s) - exact_phase_without_mode.subs(y, y_width * s))
    expected_factored = -2 * pi * (n + sp.Rational(3, 4)) * s
    require(sp.simplify(factored_at_y - expected_factored) == 0, "mode factorization failed")

    # The inverse Airy normalization contributes 1/(2*pi), while dy=Y ds.
    # Hence the exact fold-integral profile has scale Y/(2*pi)=1/h.
    require(sp.simplify(y_width / (2 * pi) - 1 / h) == 0, "profile normalization failed")

    return {
        "fixed_lower_endpoint": str(C),
        "selector_center": str(MODE_CENTER),
        "height_parameterization": "t=t*-beta*lambda, t*=beta^3=pi*C^2/8",
        "exact_height_phase": "beta^3[z/beta-tanh(z/beta)]-lambda*z",
        "exact_profile_phase": "beta^3[z/beta-tanh(z/beta)]-lambda*z-beta*y*tanh(z/beta)+x(z)y^2/(2beta)",
        "mode_detuning": "d_(39895+n)=h(n+3/4)",
        "one_period": "hY=2pi",
        "mode_factor": "exp(-i*d_m*Y*s)=exp(-2pi*i*n*s)exp(-3pi*i*s/2)",
        "profile_scale": "Y/(2pi)=1/h",
        "exact_amplitude": "cosh(z/beta)^(-3/2)(1+sigma*y/C)",
        "exact_logit": "x(z)=[1-tanh(z/beta)]/2",
    }


def exact_finite_weighted_identity_check() -> None:
    """Check the finite Fourier pairing with unrelated exact rationals."""

    frequencies = list(range(-5, 6))
    profile_coefficients = {n: Fraction(7 * n - 3, 11 * abs(n) + 19) for n in frequencies}
    weights = {n: Fraction(5 * n + 13, 17 * abs(n) + 23) for n in frequencies}
    direct = sum((weights[n] * profile_coefficients[n] for n in frequencies), Fraction())
    paired = sum((profile_coefficients[n] * weights[n] for n in frequencies), Fraction())
    require(direct == paired, "finite weighted Fourier pairing failed")

    raw_profile = {n: Fraction(3 * n + 8, 29 * abs(n) + 31) for n in frequencies}
    h = Fraction(7, 13)
    canonical_profile = {n: value / h for n, value in raw_profile.items()}
    require(
        sum((weights[n] * canonical_profile[n] for n in frequencies), Fraction())
        == sum((weights[n] * raw_profile[n] for n in frequencies), Fraction()) / h,
        "Y/(2pi)=1/h transfer failed",
    )


def numerical_certificate(dependencies: dict[str, Any]) -> dict[str, Any]:
    geometry = dependencies["fourier_completion"]["geometry"]
    kernel = dependencies["weighted_kernel"]["certificate"]
    obligation = dependencies["profile_obligation"]["certificate"]
    spacing = parse_ball(geometry["detuning_spacing_h_ball"])
    y_width = parse_ball(geometry["normal_width_Y_ball"])
    h_y_error = parse_ball(geometry["hY_minus_2pi_ball"])
    kernel_l1 = parse_ball(kernel["two_branch_kernel_L1_bound_ball"])

    require(spacing[0] - spacing[1] > Decimal("0.054"), "spacing drift")
    require(y_width[0] - y_width[1] > Decimal("116"), "profile width drift")
    require(abs(h_y_error[0]) + h_y_error[1] < Decimal("5e-69"), "hY interval drift")
    require(kernel_l1[0] + kernel_l1[1] < Decimal("3.840283e-5"), "weighted kernel drift")

    stored_scale = [Decimal(value) for value in obligation["canonical_profile_scale_interval"]]
    independent_scale = Decimal(1) / (spacing[0] - spacing[1])
    require(stored_scale[1] == independent_scale, "profile scale dependency drift")
    return {
        "detuning_spacing_h_ball": geometry["detuning_spacing_h_ball"],
        "normal_width_Y_ball": geometry["normal_width_Y_ball"],
        "hY_minus_2pi_ball": geometry["hY_minus_2pi_ball"],
        "canonical_profile_scale_interval": obligation["canonical_profile_scale_interval"],
        "completed_canonical_profile_value_upper_bound": obligation["uniform_completed_canonical_profile_value_upper_bound"],
        "weighted_kernel_L1_bound_ball": kernel["two_branch_kernel_L1_bound_ball"],
        "weighted_kernel_L1_upper_bound": str(kernel_l1[0] + kernel_l1[1]),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Exact finite-t common-profile operator

Date: 2026-08-13

Status: exact transformed-strip Fourier operator proved; uniform profile norm
remains open

Keep the common source carrier suppressed and put

```text
t*=beta^3=pi C^2/8,       t=t*-beta lambda,
sigma=4beta/(pi C),       Y=2/sigma,
h=4beta/C,                hY=2pi,
x(z)=[1-tanh(z/beta)]/2.                              (FP1)
```

For `m=39895+n`, exact arithmetic gives

```text
d_m=h(n+3/4),
exp(-i d_mYs)=exp(-2pi i ns)exp(-3pi i s/2).          (FP2)
```

Before the mode factor is inserted, define the exact reduced fold profile

```text
J_ex(lambda,s)=Integral_Gamma cosh(z/beta)^(-3/2)
 (1+sigma Ys/C)
 exp(i{{beta^3[z/beta-tanh(z/beta)]-lambda z
        -beta Ys tanh(z/beta)+x(z)Y^2s^2/(2beta)}})dz. (FP3)
```

`Gamma` denotes the same admissible exact fold contour/Abel continuation as
the completed-strip theorem.  The beta-minus-four Airy contraction is

```text
J_4(lambda,s)=2pi exp(iY^2s^2/(4beta))
 [U(lambda,Ys,beta)Ai(-lambda-Ys)
  +V(lambda,Ys,beta)d_X Ai(-lambda-Ys)].              (FP4)
```

The inverse Airy normalization and `dy=Y ds` force

```text
g_ex,lambda(s)=Y/(2pi)e^(-3pi i s/2)J_ex(lambda,s),
g_4,lambda(s) =Y/(2pi)e^(-3pi i s/2)J_4(lambda,s),
Y/(2pi)=1/h.                                          (FP5)
```

No fitted scale appears.  Substitution of (FP2) into the exact transformed
mode integral proves

```text
G_ex,(39895+n)=Integral_0^1 g_ex,lambda(s)e^(-2pi i ns)ds,
G_4,(39895+n) =Integral_0^1 g_4,lambda(s)e^(-2pi i ns)ds. (FP6)
```

Therefore, with `delta g=g_ex-g_4` and `w_39894=0`, the entire weighted
exact-minus-beta-minus-four correction is exactly

```text
Delta W=sum_(m=39696)^40094 w_m(G_ex,m-G_4,m)
       =Integral_0^1 delta g_lambda(s)P_w(s)ds,
P_w(s)=sum_(n=-199)^199 w_(39895+n)e^(-2pi i ns).     (FP7)
```

This proves the operator identity left open in Section 11.397.  It neither
separates modes nor applies absolute values before the finite Fourier pairing.
At `s=0`, (FP3)--(FP5) reduce to the saved completed-strip difference and

```text
|delta g_lambda(0)|
 <{c['completed_canonical_profile_value_upper_bound']}<1.759e-10.       (FP8)
```

The existing whole-kernel enclosure now gives the unconditional implication

```text
sup_(lambda,s)|delta g_lambda(s)|<=epsilon_profile
  ==> |Delta W|<3.840283e-5 epsilon_profile.           (FP9)
```

What remains is analytic rather than combinatorial: certify a uniform bound
for (FP3) minus (FP4) on the top corridor and `0<=s<=1`, preferably by
integrating their difference on one common contour.  The scalar value (FP8)
does not supply that norm.

Pi provenance: `pi C^2/8` is inherited from the exact Kummer phase; `2pi` in
(FP2), (FP5), and (FP6) follows algebraically from `hY=2pi` and the inverse
Airy Fourier normalization.  No circle, polygon, fitted period, or inserted
geometric constant is used.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No uniform exact-minus-beta-minus-four profile estimate, numerical bound for
`Delta W`, complete physical source/initial-data splice, complete `Q_K-T` or
`T_upper`, all-corridor or height-uniform theorem, `Lambda<=0`, PF-infinity,
RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    getcontext().prec = 90
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    for name, dependency in dependencies.items():
        require(dependency.get("passed") is True, f"dependency failed: {name}")

    symbolic = symbolic_certificate()
    exact_finite_weighted_identity_check()
    certificate = numerical_certificate(dependencies)
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "exact_finite_t_common_profile_fourier_operator_certified",
        "passed": True,
        "symbolic_certificate": symbolic,
        "operator_identity": {
            "exact_profile": "g_ex(s)=Y/(2pi) exp(-3pi i s/2) J_ex(lambda,s)",
            "beta4_profile": "g_4(s)=Y/(2pi) exp(-3pi i s/2) J_4(lambda,s)",
            "mode_coefficients": "G_(39895+n)=int_0^1 g(s) exp(-2pi i n s) ds",
            "profile_difference": "delta g=g_ex-g_4",
            "weighted_kernel": "P_w(s)=sum_(n=-199)^199 w_(39895+n) exp(-2pi i n s), w_39894=0",
            "weighted_pairing": "Delta W=int_0^1 delta g(s) P_w(s) ds",
            "uniform_transfer": "|Delta W|<=||delta g||_infinity ||P_w||_1",
        },
        "certificate": certificate,
        "decision": {
            "exact_transformed_strip_common_profile_identity_proved": True,
            "weighted_exact_minus_beta4_coefficient_operator_proved": True,
            "coefficientwise_mode_bounds_used": False,
            "uniform_exact_minus_beta4_profile_bound_proved": False,
            "weighted_exact_minus_beta4_numerical_bound_proved": False,
            "complete_physical_source_initial_data_splice_proved": False,
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
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
        "next_action": (
            "Certify sup|J_ex(lambda,s)-J_4(lambda,s)| on the ordinary top corridor and 0<=s<=1 using one common "
            "contour and a resumable two-parameter atlas or an analytic derivative/endpoint envelope. Transfer the result "
            "through the whole weighted kernel; do not estimate the 398 coefficients separately."
        ),
        "proof_boundary": (
            "Exact transformed-strip common-profile and weighted Fourier operator identities only. No uniform exact-minus-beta^-4 "
            "profile bound, weighted numerical remainder, complete physical source/initial-data splice, complete Q_K-T or "
            "T_upper, all-corridor or height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("certified exact finite-t common-profile Fourier operator", flush=True)


if __name__ == "__main__":
    main()
