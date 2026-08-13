#!/usr/bin/env python3
"""Complete the selector strip in the beta^-4 Airy basis and split its branch block."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_finite_t_completed_strip_error_gate as point_gate


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "leading_completed_strip": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_nonzero_height_completed_strip_cell_gate.json",
    "point_completed_strip": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_finite_t_completed_strip_error_gate.json",
    "boundary_fourier_completion": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_399_mode_fourier_completion_gate.json",
    "beta4_airy_reduction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_canonical_airy_derivative_reduction_gate.json",
    "branch_reflection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_opposite_branch_roster_reflection_pairing_gate.json",
    "selector_reassembly": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_strip_poisson_reassembly_gate.json",
}

LOWER_ALPHA = 159_577
FIXED_B = 5_122_423
BLOCK_LO = 39_696
BLOCK_HI = 40_094
PRECISION = 100


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
    lam, y, epsilon = sp.symbols("lambda y epsilon", real=True)
    u1 = -(13 * lam + 3 * y) / 60
    v1 = (8 * lam**2 - 4 * lam * y + 3 * y**2) / 60
    u2 = -(448 * lam**5 + 280 * lam**2 * y**3 + 4565 * lam**2 - 105 * lam * y**4 + 30 * lam * y + 63 * y**5 - 405 * y**2) / 50400
    v2 = (40 * lam**3 - 20 * lam**2 * y + lam * y**2 - 9 * y**3 + 27) / 1680
    u = 1 + epsilon * u1 + epsilon**2 * u2
    v = epsilon * v1 + epsilon**2 * v2
    u0 = sp.expand(u.subs(y, 0))
    v0 = sp.expand(v.subs(y, 0))
    require(sp.expand(u0.subs(lam, 0) - 1) == 0, "beta4 U endpoint value drift")
    require(sp.expand(v0.subs(lam, 0) - sp.Rational(9, 560) * epsilon**2) == 0, "beta4 V endpoint value drift")
    require(sp.expand(sp.diff(u0, lam).subs(lam, 0) + sp.Rational(13, 60) * epsilon) == 0, "beta4 U endpoint derivative drift")
    require(sp.expand(sp.diff(v0, lam).subs(lam, 0)) == 0, "beta4 V endpoint derivative drift")

    q = sp.symbols("q", nonzero=True)
    c_plus, c_minus = sp.symbols("C_plus C_minus")
    negative_phases = sum(q ** (4 * j + 1) for j in range(199))
    positive_phases = sum(q ** (-(4 * j + 3)) for j in range(200))
    selected = c_minus * negative_phases + c_plus * positive_phases
    opposite = c_plus * negative_phases + c_minus * positive_phases
    full = (c_plus + c_minus) * (negative_phases + positive_phases)
    require(sp.expand(selected + opposite - full) == 0, "selected/opposite block split failed")

    selected_reflection = c_plus * q**-799
    selected_reflection += sum(
        c_minus * q ** (4 * j + 1) + c_plus * q ** (-(4 * j + 3))
        for j in range(199)
    )
    require(sp.expand(selected - selected_reflection) == 0, "selected reflection pairing failed")

    return {
        "beta4_endpoint_U": str(u0),
        "beta4_endpoint_V": str(v0),
        "completed_beta4_profile": "g4_lambda(s)=Y*exp(i*Y^2*s^2/(4*beta)-3*pi*i*s/2)*(U*A+V*A_X)(lambda,Y*s)",
        "boundary_beta4_coefficient": "G4_(39895+k)=Integral_0^1 g4_lambda(s)*exp(-2*pi*i*k*s) ds, -199<=k<=199",
        "symmetric_midpoint": "M4=(g4_lambda(0)+g4_lambda(1))/2",
        "endpoint_half_current": "H4=(g4_lambda(0)-g4_lambda(1))/2",
        "outer_complement": "C4=M4-sum_(k=-199)^199 G4_(39895+k)",
        "completed_identity": "H4+sum_block G4+C4=g4_lambda(0)",
        "boundary_labels": "4m-A=-793,-789,...,799 for m=39696,...,40094",
        "selected_branch_sum": "B_sel=C_-*sum_(j=0)^198 exp(i*(4j+1)*x)+C_+*sum_(j=0)^199 exp(-i*(4j+3)*x)",
        "selected_reflection_form": "B_sel=exp(-i*799*x)C_+ + exp(-i*x)*sin(398*x)/sin(2*x)*(exp(i*398*x)C_-+exp(-i*398*x)C_+)",
        "opposite_branch_sum": "B_opp=C_+*sum_(j=0)^198 exp(i*(4j+1)*x)+C_-*sum_(j=0)^199 exp(-i*(4j+3)*x)",
        "block_split": "sum_block G4=Integral_0^Y exp(i*y^2/(4*beta))*(B_sel+B_opp)dy",
        "completed_remainder": "R_comp=H4+C4+Integral_0^Y exp(i*y^2/(4*beta))*B_opp dy=g4_lambda(0)-Integral_0^Y exp(i*y^2/(4*beta))*B_sel dy",
        "removable_rule": "The sine quotient is its original finite exponential sum when sin(2*x)=0.",
    }


def beta4_correction_second_derivative(beta: arb, lambda_radius: arb, pi: arb) -> arb:
    epsilon = 1 / (beta * beta)
    lam = arb(0, lambda_radius)
    lam2 = lam * lam
    lam3 = lam2 * lam
    lam4 = lam3 * lam
    lam5 = lam4 * lam

    delta_u = epsilon * (-arb(13) * lam / 60)
    delta_u += epsilon * epsilon * (-(arb(448) * lam5 + arb(4565) * lam2) / 50400)
    v = epsilon * (arb(2) * lam2 / 15)
    v += epsilon * epsilon * ((arb(40) * lam3 + 27) / 1680)
    delta_u_prime = epsilon * (-arb(13) / 60)
    delta_u_prime += epsilon * epsilon * (-(arb(2240) * lam4 + arb(9130) * lam) / 50400)
    delta_u_second = epsilon * epsilon * (-(arb(8960) * lam3 + arb(9130)) / 50400)
    v_prime = epsilon * (arb(4) * lam / 15) + epsilon * epsilon * (arb(120) * lam2 / 1680)
    v_second = epsilon * (arb(4) / 15) + epsilon * epsilon * (arb(240) * lam / 1680)

    p = delta_u_prime - lam * v
    q = delta_u + v_prime
    p_prime = delta_u_second - v - lam * v_prime
    q_prime = delta_u_prime + v_second
    a_coefficient = p_prime - lam * q
    ax_coefficient = p + q_prime
    airy = (-lam).airy_ai()
    airy_x = -(-lam).airy_ai(derivative=1)
    return abs(2 * pi * (a_coefficient * airy + ax_coefficient * airy_x))


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["center_certificate"]
    h = artifact["height_cell_certificate"]
    return f"""# Beta^-4 completed selector strip and branch projection

Date: 2026-08-13

Status: rigorous beta-minus-four all-mode completion and selector-height error
cell; not an ordinary-Morse splice or a proof of `T_upper` or RH

For `A={LOWER_ALPHA}`, `B={FIXED_B}`, `beta^3=pi*A^2/8`, and
`lambda=(t*-t)/beta`, the beta-minus-four Airy reduction has endpoint weights

```text
U(lambda,0)=1-13*lambda/(60*beta^2)
 -[448*lambda^5+4565*lambda^2]/(50400*beta^4),
V(lambda,0)=2*lambda^2/(15*beta^2)
 +[40*lambda^3+27]/(1680*beta^4).                  (BP1)
```

The shared-boundary coefficients are the Fourier coefficients of the smooth
profile `g4_lambda`.  Dirichlet-Jordan and the matching endpoint half-current
therefore give the exact completed identity

```text
H4 + sum_(m=39696)^40094 G4_m + C4 = g4_lambda(0),   (BP2)
```

where `C4` is the zero, negative, and outer-positive complement retained as
one object.  No outer coefficient is bounded or discarded separately.

At `lambda=0`, (BP1) gives

```text
J4-J0=-(9/560)*2*pi*Ai'(0)/beta^4
 ={c['beta4_minus_leading_reduced_ball']},
J_exact-J4={c['exact_minus_beta4_reduced_ball']}.       (BP3)
```

Thus the beta-minus-four term cancels the previously certified exact-minus-
leading discrepancy by a factor greater than
`{c['center_reduced_improvement_factor_lower_ball']}`.  The residual is not
set to zero; it remains an explicitly enclosed finite-height error.

The same subtraction is rigorous on the complete selector cell.  The exact
and beta-minus-four center derivatives differ by
`{h['center_first_derivative_residual_ball']}`.  Using the old exact-minus-
leading second-derivative majorant together with an interval bound for the
beta-minus-four correction proves

```text
t*-pi/16 <= t <= t*+pi/16,
|J_exact(lambda)-J4(lambda)|
 < {h['uniform_reduced_exact_minus_beta4_bound_ball']},
equation-(9) normalized error
 < {h['uniform_equation9_normalized_exact_minus_beta4_bound_ball']} < 4.9e-16. (BP4)
```

On the ordinary-side half `t*-pi/16<=t<=t*`, the positive-`X` branch basis is
valid away from its regular endpoint limit.  In the old-strip coordinate
`x=beta*y/A`, the 399 block splits exactly into selected and opposite pieces,

```text
B_sel=e^(-i*799*x)C_+
 +e^(-i*x) sin(398*x)/sin(2*x)
   [e^(i*398*x)C_-+e^(-i*398*x)C_+],                 (BP5)

sum_block G4 = Integral exp(i*y^2/(4*beta))
                    (B_sel+B_opp)dy.                 (BP6)
```

The admissible completed remainder is consequently

```text
R_comp=H4+C4+Integral exp(i*y^2/(4*beta))B_opp dy
      =g4_lambda(0)-Integral exp(i*y^2/(4*beta))B_sel dy. (BP7)
```

Equation (BP7), rather than `B_opp` alone, is the cancellation-preserving
object.  It inserts the reflection pairing into the completed Poisson theorem
without reallocating the complement or endpoint half-current.  The adjacent
`A+2` event chart still has different detunings and is not identified here.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No selected-branch/logistic amplitude match, all-corridor estimate, complete
`Q_K-T` or `T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion
is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    require(dependencies["leading_completed_strip"].get("passed") is True, "nonzero completed-strip dependency failed")
    require(dependencies["point_completed_strip"].get("passed") is True, "point completed-strip dependency failed")
    require(dependencies["beta4_airy_reduction"]["decision"]["all_derivatives_reduced_to_Ai_and_first_derivative"] is True, "beta4 Airy reduction drift")
    require(dependencies["branch_reflection"]["decision"]["shared_boundary_399_mode_opposite_branch_is_exact_mirror_pairing"] is True, "branch reflection drift")

    ctx.dps = PRECISION
    algebra = symbolic_certificate()
    p = point_gate.parameters()
    beta = p["beta"]
    pi = p["pi"]
    epsilon = 1 / (beta * beta)
    lambda_radius = pi / (16 * beta)
    airy_zero = arb(0).airy_ai()
    airy_prime_zero = arb(0).airy_ai(derivative=1)
    airy_x_zero = -airy_prime_zero

    point = dependencies["point_completed_strip"]
    cell = dependencies["leading_completed_strip"]
    exact_minus_leading = arb(point["numerical_certificate"]["reduced_exact_minus_canonical_fold_ball"])
    beta4_minus_leading = 2 * pi * arb(9) / 560 * epsilon * epsilon * airy_x_zero
    center_residual = exact_minus_leading - beta4_minus_leading
    require(abs(center_residual) < arb("6.6e-22"), "center exact-minus-beta4 residual exceeds 6.6e-22")

    exact_minus_leading_derivative = arb(cell["center_derivative"]["rigorous_derivative_ball"])
    beta4_minus_leading_derivative = -arb(13) / 60 * 2 * pi * airy_zero * epsilon
    derivative_residual = exact_minus_leading_derivative - beta4_minus_leading_derivative
    require(abs(derivative_residual) < arb("2.6e-22"), "center derivative residual exceeds 2.6e-22")

    exact_minus_leading_second = arb(cell["second_derivative_majorant"]["two_ray_second_derivative_bound"])
    beta4_second = beta4_correction_second_derivative(beta, lambda_radius, pi)
    require(beta4_second < arb("5.85e-8"), "beta4 correction second derivative exceeds 5.85e-8")
    residual_second = exact_minus_leading_second + beta4_second
    uniform_reduced = abs(center_residual) + lambda_radius * abs(derivative_residual) + lambda_radius * lambda_radius * residual_second / 2
    require(uniform_reduced < arb("5.3e-15"), "uniform exact-minus-beta4 reduced error exceeds 5.3e-15")

    strip_scale = arb(2).sqrt() / pi * p["Y"]
    height_lower = p["height"] - pi / 16
    maximum_normalizer = (pi / (32 * height_lower)) ** (arb(1) / 4)
    center_physical = strip_scale * center_residual
    center_normalized = (pi / (32 * p["height"])) ** (arb(1) / 4) * center_physical
    uniform_physical = strip_scale * uniform_reduced
    uniform_normalized = maximum_normalizer * uniform_physical
    require(uniform_normalized < arb("4.9e-16"), "uniform normalized exact-minus-beta4 error exceeds 4.9e-16")

    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "beta4_all_mode_selector_completion_and_branch_projection_certified",
        "passed": True,
        "parameters": {
            "lower_endpoint_A": LOWER_ALPHA,
            "fixed_upper_endpoint_B": FIXED_B,
            "boundary_mode_roster": [BLOCK_LO, BLOCK_HI],
            "boundary_detuning_labels": [-793, 799, 4],
            "beta_ball": beta.str(PRECISION, more=True),
            "lambda_radius_ball": lambda_radius.str(PRECISION, more=True),
            "height_interval": "t*-pi/16 <= t <= t*+pi/16",
            "ordinary_side_half_cell": "t*-pi/16 <= t <= t* (equivalently 0<=lambda<=pi/(16beta))",
        },
        "exact_algebra": algebra,
        "center_certificate": {
            "exact_minus_leading_reduced_ball": exact_minus_leading.str(PRECISION, more=True),
            "beta4_minus_leading_reduced_ball": beta4_minus_leading.str(PRECISION, more=True),
            "exact_minus_beta4_reduced_ball": center_residual.str(PRECISION, more=True),
            "exact_minus_beta4_physical_ball": center_physical.str(PRECISION, more=True),
            "exact_minus_beta4_equation9_normalized_ball": center_normalized.str(PRECISION, more=True),
            "center_reduced_improvement_factor_lower_ball": (abs(exact_minus_leading) / abs(center_residual)).lower().str(PRECISION, more=True),
            "beta4_endpoint_identity": "J4(0)-J0(0)=-(9/560)*2*pi*Ai'(0)/beta^4",
        },
        "height_cell_certificate": {
            "center_exact_minus_leading_derivative_ball": exact_minus_leading_derivative.str(PRECISION, more=True),
            "center_beta4_minus_leading_derivative_ball": beta4_minus_leading_derivative.str(PRECISION, more=True),
            "center_first_derivative_residual_ball": derivative_residual.str(PRECISION, more=True),
            "exact_minus_leading_second_derivative_bound_ball": exact_minus_leading_second.str(PRECISION, more=True),
            "beta4_correction_second_derivative_bound_ball": beta4_second.str(PRECISION, more=True),
            "exact_minus_beta4_second_derivative_bound_ball": residual_second.str(PRECISION, more=True),
            "Taylor_bound": "|r(lambda)|<=|r(0)|+L|r'(0)|+L^2*(sup|delta''|+sup|c4''|)/2",
            "uniform_reduced_exact_minus_beta4_bound_ball": uniform_reduced.str(PRECISION, more=True),
            "uniform_physical_exact_minus_beta4_bound_ball": uniform_physical.str(PRECISION, more=True),
            "maximum_equation9_normalizer_ball": maximum_normalizer.str(PRECISION, more=True),
            "uniform_equation9_normalized_exact_minus_beta4_bound_ball": uniform_normalized.str(PRECISION, more=True),
            "uniform_improvement_over_leading_cell_lower_ball": (arb(cell["uniform_error"]["uniform_reduced_completed_strip_error_bound_ball"]) / uniform_reduced).lower().str(PRECISION, more=True),
        },
        "decision": {
            "beta4_all_mode_Fourier_completion_proved": True,
            "beta4_endpoint_half_current_retained": True,
            "zero_negative_outer_positive_complement_retained_as_one_object": True,
            "selected_plus_opposite_branch_equals_boundary_block": True,
            "completed_opposite_remainder_defined_without_termwise_outer_bounds": True,
            "uniform_exact_minus_beta4_selector_cell_below_4_9e_minus_16_normalized": True,
            "positive_X_branch_projection_scoped_only_to_ordinary_side_half_cell": True,
            "adjacent_selector_amplitude_identification_proved": False,
            "selected_branch_logistic_Gamma_match_proved": False,
            "ordinary_Morse_corridor_splice_proved": False,
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
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
        "next_action": (
            "Use B_sel as the sole branch-continuation candidate and R_comp as the inseparable completed remainder. "
            "Derive the exact finite-height selected-branch/logistic-Morse amplitude map on the top ordinary corridor, "
            "then prove or reject its endpoint-uniform remainder before attempting all 400 corridors."
        ),
        "proof_boundary": (
            "Exact beta^-4 all-mode selector completion, cancellation-preserving branch projection, and a rigorous "
            "exact-minus-beta4 error on one selector cell only. No adjacent-chart amplitude identification, selected-branch/logistic "
            "match, ordinary-corridor splice, complete Q_K-T or T_upper theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("built beta^-4 completed selector branch projection: normalized cell residual < 4.9e-16", flush=True)


if __name__ == "__main__":
    main()
