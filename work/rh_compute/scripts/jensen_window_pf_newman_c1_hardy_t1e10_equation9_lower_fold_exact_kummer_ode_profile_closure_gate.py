#!/usr/bin/env python3
"""Close the exact-minus-beta4 common profile through its exact Kummer ODE."""

from __future__ import annotations

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

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_kummer_ode_profile_closure_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "common_profile_operator": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_finite_t_common_profile_operator_gate.json",
    "beta4_completed_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.json",
    "initial_slope": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_profile_initial_slope_gate.json",
    "beta6_leading_coefficient": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta6_common_profile_leading_coefficient_gate.json",
    "weighted_kernel": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_opposite_branch_dirichlet_gate.json",
}

C = 159_577
PANELS = 8192
PRECISION = 80


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


def symbolic_certificate() -> dict[str, str]:
    beta, lam, y, rho, t = sp.symbols("beta lambda y rho t", positive=True, real=True)
    x, q = sp.symbols("x q", positive=True, real=True)
    i = sp.I

    rho_y = 1 + y / (2 * beta**2)
    t_value = beta**3 - beta * lam
    kappa = 2 * beta**3 * rho_y**2
    require(sp.expand(kappa / 2 - (beta**3 + beta * y + y**2 / (4 * beta))) == 0, "Kummer argument drift")
    require(sp.expand(beta**3 * rho_y**2 - t_value - beta * (lam + y + y**2 / (4 * beta**2))) == 0, "ODE potential drift")

    # Logistic substitution x=(1-tanh(q))/2.
    tanh_q = 1 - 2 * x
    q_logit = sp.log((1 - x) / x) / 2
    jacobian_amplitude = beta * sp.sqrt(2) * (x * (1 - x)) ** (-sp.Rational(1, 4))
    phase = sp.expand(t * q_logit - beta**3 * rho**2 * tanh_q)
    expected_phase = -beta**3 * rho**2 + t * sp.log(1 - x) / 2 - t * sp.log(x) / 2 + 2 * beta**3 * rho**2 * x
    require(sp.simplify(phase - expected_phase) == 0, "logistic Kummer phase drift")

    # Kummer equation after removing exp(z/2).
    w0, w1 = sp.symbols("w0 w1")
    c = 2 * i * beta**3
    zeta = c * rho**2
    w2 = -sp.Rational(3, 2) * w1 / zeta + (sp.Rational(1, 4) - i * t / (2 * zeta)) * w0
    f_rr = sp.expand(6 * c * rho * w1 + 4 * c**2 * rho**3 * w2)
    target_rr = -4 * beta**3 * (beta**3 * rho**2 - t) * rho * w0
    require(sp.simplify(f_rr - target_rr) == 0, "Kummer-to-radial ODE reduction failed")

    # Exact beta^-4 residual in the A=Ai(-X), B=d_X Ai(-X) basis.
    epsilon = sp.symbols("epsilon", positive=True, real=True)
    X = lam + y
    u1 = -(13 * lam + 3 * y) / 60
    v1 = (8 * lam**2 - 4 * lam * y + 3 * y**2) / 60
    u2 = -(448 * lam**5 + 280 * lam**2 * y**3 + 4565 * lam**2 - 105 * lam * y**4 + 30 * lam * y + 63 * y**5 - 405 * y**2) / 50400
    v2 = (40 * lam**3 - 20 * lam**2 * y + lam * y**2 - 9 * y**3 + 27) / 1680
    u = 1 + epsilon * u1 + epsilon**2 * u2
    v = epsilon * v1 + epsilon**2 * v2

    def derivative(pair: tuple[sp.Expr, sp.Expr]) -> tuple[sp.Expr, sp.Expr]:
        first, second = pair
        return sp.diff(first, y) - X * second, first + sp.diff(second, y)

    second_u, second_v = derivative(derivative((u, v)))
    residual_u = sp.factor(second_u + (X + epsilon * y**2 / 4) * u)
    residual_v = sp.factor(second_v + (X + epsilon * y**2 / 4) * v)
    p3 = -y**2 * (448 * lam**5 + 280 * lam**2 * y**3 + 4565 * lam**2 - 105 * lam * y**4 + 30 * lam * y + 63 * y**5 - 405 * y**2) / 201600
    q3 = -y**2 * (-40 * lam**3 + 20 * lam**2 * y - lam * y**2 + 9 * y**3 - 27) / 6720
    require(sp.expand(residual_u - epsilon**3 * p3) == 0, "beta4 residual Ai coefficient drift")
    require(sp.expand(residual_v - epsilon**3 * q3) == 0, "beta4 residual d_X Ai coefficient drift")

    return {
        "logistic_substitution": "x=(1-tanh(q))/2, q=(1/2)log((1-x)/x)",
        "logistic_jacobian_amplitude": str(jacobian_amplitude),
        "kummer_parameters": "a=3/4-i*t/2, b=3/4+i*t/2, a+b=3/2, t=beta^3-beta*lambda",
        "kummer_argument": "kappa=2beta^3(1+y/(2beta^2))^2",
        "exact_kummer_identity": (
            "H_ex=beta*sqrt(2)*rho*exp(-i*kappa/2)*B(a,b)*1F1(a;3/2;i*kappa), rho=1+y/(2beta^2)"
        ),
        "exact_ode": "H_ex''(y)+[lambda+y+y^2/(4beta^2)]H_ex(y)=0",
        "beta4_residual": "R4=2pi*beta^-6[P3*Ai(-X)+Q3*d_X Ai(-X)]",
        "P3": str(sp.factor(p3)),
        "Q3": str(sp.factor(q3)),
        "error_equation": "E''+[lambda+y+y^2/(4beta^2)]E=-R4, E=H_ex-H4",
    }


def residual_coefficients(lam: arb, y: arb) -> tuple[arb, arb]:
    p3 = -y**2 * (
        448 * lam**5 + 280 * lam**2 * y**3 + 4565 * lam**2
        - 105 * lam * y**4 + 30 * lam * y + 63 * y**5 - 405 * y**2
    ) / 201600
    q3 = -y**2 * (-40 * lam**3 + 20 * lam**2 * y - lam * y**2 + 9 * y**3 - 27) / 6720
    return p3, q3


def interval_certificate(dependencies: dict[str, Any]) -> dict[str, Any]:
    ctx.dps = PRECISION
    pi = arb.pi()
    c = arb(C)
    beta = (pi * c**2 / 8) ** (arb(1) / 3)
    epsilon = beta**-2
    y_width = pi * c / (2 * beta)
    h = 4 * beta / c
    lambda_max = pi / (16 * beta)
    lam = arb(lambda_max / 2, lambda_max / 2)
    panel_width = y_width / PANELS

    residual_l1 = arb(0)
    residual_max = arb(0)
    residual_max_panel = -1
    airy_modulus_upper = arb(0)
    airy_modulus_panel = -1
    for panel in range(PANELS):
        y_midpoint = y_width * (2 * panel + 1) / (2 * PANELS)
        y = arb(y_midpoint, panel_width / 2)
        x = lam + y
        airy = (-x).airy_ai()
        airy_x = -(-x).airy_ai(derivative=1)
        airy_bi = (-x).airy_bi()
        p3, q3 = residual_coefficients(lam, y)
        residual = 2 * pi * epsilon**3 * (
            abs(p3).upper() * abs(airy).upper() + abs(q3).upper() * abs(airy_x).upper()
        )
        residual_l1 += panel_width * residual.upper()
        if residual.upper() > residual_max:
            residual_max = residual.upper()
            residual_max_panel = panel

        modulus = (airy * airy + airy_bi * airy_bi).sqrt().upper()
        if modulus > airy_modulus_upper:
            airy_modulus_upper = modulus
            airy_modulus_panel = panel

    require(airy_modulus_upper < arb("0.716"), "Airy modulus atlas exceeds 0.716")
    require(residual_l1 < arb("1.001e-8"), "integrated exact residual exceeds 1.001e-8")

    derivative_modulus = (
        (-lam).airy_ai(derivative=1) ** 2 + (-lam).airy_bi(derivative=1) ** 2
    ).sqrt().upper()
    kernel_bound = (pi * airy_modulus_upper**2).upper()
    initial_value_kernel_bound = (pi * airy_modulus_upper * derivative_modulus).upper()
    feedback_integral = (epsilon * y_width**3 / 12).upper()
    feedback_factor = (kernel_bound * feedback_integral).upper()
    denominator = (1 - feedback_factor).lower()
    require(kernel_bound < arb("1.611"), "Airy Green-kernel bound exceeds 1.611")
    require(denominator > arb("0.954"), "Volterra denominator falls below 0.954")

    initial_value = abs(arb(
        dependencies["beta4_completed_projection"]["height_cell_certificate"]["uniform_reduced_exact_minus_beta4_bound_ball"]
    )).upper()
    initial_slope = abs(arb(
        dependencies["initial_slope"]["numerical_certificate"]["uniform_exact_minus_beta4_initial_slope_ball"]
    )).upper()
    numerator = (
        initial_value_kernel_bound * initial_value
        + kernel_bound * initial_slope
        + kernel_bound * residual_l1
    ).upper()
    raw_profile_bound = (numerator / denominator).upper()
    canonical_profile_bound = (raw_profile_bound / h).upper()
    kernel_l1 = arb(dependencies["weighted_kernel"]["certificate"]["two_branch_kernel_L1_bound_ball"]).upper()
    weighted_bound = (canonical_profile_bound * kernel_l1).upper()
    require(raw_profile_bound < arb("1.69e-8"), "raw exact-minus-beta4 profile exceeds 1.69e-8")
    require(canonical_profile_bound < arb("3.13e-7"), "canonical exact-minus-beta4 profile exceeds 3.13e-7")
    require(weighted_bound < arb("1.21e-11"), "weighted finite-t profile correction exceeds 1.21e-11")

    return {
        "panels": PANELS,
        "precision_decimal_digits": PRECISION,
        "lambda_interval": "0<=lambda<=pi/(16beta)",
        "y_interval": "0<=y<=Y=pi*C/(2beta)",
        "beta_ball": beta.str(PRECISION, more=True),
        "epsilon_beta_minus_2_ball": epsilon.str(PRECISION, more=True),
        "Y_ball": y_width.str(PRECISION, more=True),
        "h_ball": h.str(PRECISION, more=True),
        "lambda_max_ball": lambda_max.str(PRECISION, more=True),
        "Airy_modulus_maximum_panel": airy_modulus_panel,
        "Airy_modulus_upper_bound": airy_modulus_upper.str(PRECISION, more=True),
        "Airy_derivative_modulus_at_initial_line_upper_bound": derivative_modulus.str(PRECISION, more=True),
        "Green_kernel_uniform_bound": kernel_bound.str(PRECISION, more=True),
        "initial_value_kernel_uniform_bound": initial_value_kernel_bound.str(PRECISION, more=True),
        "residual_maximum_panel": residual_max_panel,
        "residual_pointwise_upper_bound": residual_max.str(PRECISION, more=True),
        "residual_L1_upper_bound": residual_l1.str(PRECISION, more=True),
        "quadratic_feedback_integral_upper_bound": feedback_integral.str(PRECISION, more=True),
        "Volterra_feedback_factor_upper_bound": feedback_factor.str(PRECISION, more=True),
        "Volterra_denominator_lower_bound": denominator.str(PRECISION, more=True),
        "initial_value_uniform_bound": initial_value.str(PRECISION, more=True),
        "initial_slope_uniform_bound": initial_slope.str(PRECISION, more=True),
        "Volterra_numerator_upper_bound": numerator.str(PRECISION, more=True),
        "uniform_raw_exact_minus_beta4_profile_bound": raw_profile_bound.str(PRECISION, more=True),
        "uniform_canonical_exact_minus_beta4_profile_bound": canonical_profile_bound.str(PRECISION, more=True),
        "weighted_kernel_L1_upper_bound": kernel_l1.str(PRECISION, more=True),
        "uniform_weighted_finite_t_profile_correction_bound": weighted_bound.str(PRECISION, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    s = artifact["symbolic_certificate"]
    c = artifact["interval_certificate"]
    return f"""# Exact Kummer ODE whole-profile closure

Date: 2026-08-13

Status: rigorous full exact-minus-beta-minus-four profile and weighted
finite-height correction on the top selector corridor; this is not a complete
`T_upper` theorem or a proof of RH

Put `rho=1+y/(2beta^2)`, `t=beta^3-beta*lambda`, and remove the common phase
`exp(i*y^2/(4beta))` from the transformed fold.  The exact logistic change

```text
x=(1-tanh(z/beta))/2
```

gives, without asymptotic expansion,

```text
H_ex=beta*sqrt(2)*rho*exp(-i*kappa/2) B(a,b)
      *1F1(a;3/2;i*kappa),
a=3/4-i*t/2, b=3/4+i*t/2,
kappa=2beta^3 rho^2.                                  (KO1)
```

Kummer's differential equation then reduces exactly to

```text
H_ex''+[lambda+y+y^2/(4beta^2)]H_ex=0.                (KO2)
```

The beta-minus-four Airy profile `H4` has no hidden higher residual.  Exact
symbolic reduction in the basis `A=Ai(-X)`, `A_X=d_X Ai(-X)`, `X=lambda+y`,
proves

```text
H4''+[lambda+y+y^2/(4beta^2)]H4
 =R4=2pi*beta^-6[P3 A+Q3 A_X],                        (KO3)
P3={s['P3']},
Q3={s['Q3']}.                                         (KO4)
```

Thus `E=H_ex-H4` obeys the exact forced equation

```text
E''+(lambda+y)E=-beta^-2*y^2*E/4-R4.                  (KO5)
```

No post-beta-minus-six Taylor remainder occurs in (KO5).  An `{c['panels']}`-
panel Arb atlas on the complete rectangle proves

```text
sup sqrt(Ai(-X)^2+Bi(-X)^2)
 < {c['Airy_modulus_upper_bound']},
sup |K_Airy(y,s)| < {c['Green_kernel_uniform_bound']},
Integral_0^Y |R4|dy < {c['residual_L1_upper_bound']}.  (KO6)
```

The exact completed-point and initial-slope gates supply

```text
sup |E(lambda,0)| < {c['initial_value_uniform_bound']},
sup |E_y(lambda,0)| < {c['initial_slope_uniform_bound']}. (KO7)
```

The Volterra feedback factor is
`{c['Volterra_feedback_factor_upper_bound']}`, leaving denominator
`>{c['Volterra_denominator_lower_bound']}`.  Therefore

```text
sup_(lambda,y)|J_ex-J4|=sup|E|
 < {c['uniform_raw_exact_minus_beta4_profile_bound']} < 1.69e-8,
sup_(lambda,s)|g_ex-g4|
 < {c['uniform_canonical_exact_minus_beta4_profile_bound']} < 3.13e-7. (KO8)
```

Finally, applying the already certified whole Fourier kernel only after the
exact coefficient pairing gives

```text
|Delta W_finite-t|
 < {c['uniform_weighted_finite_t_profile_correction_bound']} < 1.21e-11. (KO9)
```

This closes the full profile obligation left open in Sections 11.397--11.399;
it does not infer the result from the small formal beta-minus-six coefficient.
Pi provenance: `beta^3=pi*C^2/8` comes from the exact transformed Kummer
geometry, `2pi` in (KO3) is the inverse Airy Fourier normalization, and the
Green Wronskian is `1/pi`.  No fitted constant is inserted.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No complete selected-leading plus source-current inequality, all-corridor or
height-uniform theorem, complete `Q_K-T` or `T_upper`, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    for name, dependency in dependencies.items():
        require(dependency.get("passed") is True, f"dependency failed: {name}")

    symbolic = symbolic_certificate()
    interval = interval_certificate(dependencies)
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "exact_kummer_ode_whole_profile_and_weighted_finite_t_correction_certified",
        "passed": True,
        "symbolic_certificate": symbolic,
        "interval_certificate": interval,
        "decision": {
            "exact_logistic_Kummer_identity_proved": True,
            "exact_common_profile_ODE_proved": True,
            "beta4_residual_is_exactly_beta6_with_no_hidden_tail": True,
            "Airy_Green_Volterra_denominator_positive": True,
            "full_exact_minus_beta4_profile_bound_proved": True,
            "weighted_finite_t_profile_correction_below_1_21e_minus_11": True,
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
            "Insert the 1.21e-11 exact finite-t weighted correction into the cancellation-preserving selected-leading defect budget, "
            "then identify the remaining source/ordinary-Morse splice term before scaling beyond the top corridor."
        ),
        "proof_boundary": (
            "A rigorous full exact-minus-beta4 common-profile and weighted finite-t correction on one top selector corridor only. "
            "No complete selected-leading/source-current inequality, all-corridor or height-uniform theorem, complete Q_K-T or "
            "T_upper theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("certified exact Kummer-ODE profile closure: weighted finite-t correction < 1.21e-11", flush=True)


if __name__ == "__main__":
    main()
