#!/usr/bin/env python3
"""Certify the finite-Poisson joint saddle that produces equation (124)."""

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


CURRENT_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation10_global_theta_current_reduction_gate.json"
SOURCE_LEDGER = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_source_aligned_hybrid_error_ledger_gate.json"
PAPER = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_Brereton_2607.15310.pdf"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_finite_poisson_portcullis_saddle_reduction_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_finite_poisson_portcullis_saddle_reduction_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 110
T = 10_000_000_000
A = 159577
B = 5122421


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


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


def symbolic_reduction() -> dict[str, Any]:
    alpha, x, m, t = sp.symbols("alpha x m t", positive=True)
    pi = sp.pi
    phase = pi * alpha**2 * x / 4 - pi * m * alpha + t * sp.log((1 - x) / x) / 2
    phase_alpha = sp.simplify(sp.diff(phase, alpha))
    phase_x = sp.simplify(sp.diff(phase, x))
    x_star = 2 * pi * m**2 / (t + 2 * pi * m**2)
    alpha_star = 2 * m + t / (pi * m)
    stationary_alpha_residual = sp.simplify(phase_alpha.subs({alpha: alpha_star, x: x_star}))
    stationary_x_residual = sp.simplify(phase_x.subs({alpha: alpha_star, x: x_star}))
    require(stationary_alpha_residual == 0 and stationary_x_residual == 0, "joint saddle equations failed")

    hessian = sp.hessian(phase, (alpha, x))
    hessian_det = sp.factor(sp.simplify(hessian.det().subs({alpha: alpha_star, x: x_star})))
    expected_det = -(2 * pi * m**2 + t) ** 3 / (8 * m**2 * t)
    require(sp.simplify(hessian_det - expected_det) == 0, "joint saddle Hessian failed")

    equation124_residual = sp.simplify(alpha_star - 2 * m * (t / (2 * pi * m**2) + 1))
    require(equation124_residual == 0, "equation-(124) saddle identity failed")
    discriminant = sp.sqrt(1 - 8 * t / (pi * alpha**2))
    n_minus = alpha * (1 - discriminant) / 4
    n_plus = alpha * (1 + discriminant) / 4
    quadratic = 2 * m**2 - alpha * m + t / pi
    require(sp.simplify(quadratic.subs(m, n_minus)) == 0, "lower inverse root failed")
    require(sp.simplify(quadratic.subs(m, n_plus)) == 0, "upper inverse root failed")
    root_sum = sp.simplify(n_minus + n_plus)
    root_product = sp.simplify(n_minus * n_plus)
    require(root_sum == alpha / 2 and root_product == t / (2 * pi), "root involution failed")

    partner = t / (2 * pi * m)
    x_partner = sp.simplify(x_star.subs(m, partner))
    require(sp.simplify(x_partner - (1 - x_star)) == 0, "reflected saddle identity failed")

    phi_one = x * alpha**2 / 4 - m * alpha
    completed_square = x * (alpha - 2 * m / x) ** 2 / 4 - m**2 / x
    square_residual = sp.simplify(phi_one - completed_square)
    require(square_residual == 0, "mode completed square failed")
    exponential = sp.exp(sp.I * pi * phi_one)
    integrand_decomposition = sp.simplify(
        alpha * exponential
        - (2 / x) * (sp.diff(exponential, alpha) / (sp.I * pi) + m * exponential)
    )
    require(integrand_decomposition == 0, "mode integral decomposition failed")

    return {
        "finite_poisson_formula": "sum_(n=0)^L f(n)=(f(0)+f(L))/2+lim_(M->infinity) sum_(q=-M)^M integral_0^L f(u)e^(-2*pi*i*q*u)du",
        "finite_poisson_assumptions": "L is a positive integer and f is C^1 with piecewise C^2 derivative; the Fourier sum is interpreted symmetrically via the Dirichlet kernel",
        "odd_roster_mode_integral": "J_m=(-1)^m/2 integral_A^B alpha exp(i*pi*(x*alpha^2/4-m*alpha))dalpha",
        "mode_integral_identity": "J_m=(-1)^m/x*{[E_m(B)-E_m(A)]/(i*pi)+m*F_m(A,B;x)}",
        "mode_exponential": "E_m(alpha)=exp(i*pi*(x*alpha^2/4-m*alpha))",
        "mode_fresnel_integral": "F_m=exp(-i*pi*m^2/x) integral_A^B exp(i*pi*x*(alpha-2m/x)^2/4)dalpha",
        "completed_square_residual": str(square_residual),
        "mode_integrand_decomposition_residual": str(integrand_decomposition),
        "joint_phase": "Phi=pi*alpha^2*x/4-pi*m*alpha+(t/2)log((1-x)/x)",
        "stationary_x": "x_m=2*pi*m^2/(t+2*pi*m^2)",
        "stationary_alpha": "alpha_m=2*m+t/(pi*m)",
        "stationary_alpha_residual": str(stationary_alpha_residual),
        "stationary_x_residual": str(stationary_x_residual),
        "equation124_residual": str(equation124_residual),
        "hessian_determinant": str(hessian_det),
        "hessian_sign": "strictly_negative_for_m_t_positive",
        "lower_inverse_root": "N_-(alpha;t)=alpha/4*(1-sqrt(1-8*t/(pi*alpha^2)))",
        "upper_inverse_root": "N_+(alpha;t)=alpha/4*(1+sqrt(1-8*t/(pi*alpha^2)))",
        "root_sum": str(root_sum),
        "root_product": str(root_product),
        "reflection_involution": "N*=t/(2*pi*N) and x_(N*)=1-x_N",
        "integer_lattice_guard": "The continuous reflection N->t/(2*pi*N) does not generally preserve integer dual modes; reflected branches must be recombined by the exact theta/Poisson formula, not paired by rounding.",
        "pi_provenance": "pi is inherited from the Kummer quadratic phase exp(i*pi*alpha^2*x/4) and the Fourier factor exp(-2*pi*i*m*u); their joint saddle produces t/(pi*m).",
    }


def inverse_roots(alpha_value: int, precision: int) -> tuple[arb, arb]:
    ctx.dps = precision
    alpha = arb(alpha_value)
    t = arb(T)
    disc = (1 - 8 * t / (arb.pi() * alpha**2)).sqrt()
    return alpha * (1 - disc) / 4, alpha * (1 + disc) / 4


def stationary_classification(precision: int) -> dict[str, Any]:
    ctx.dps = precision
    lower_a, upper_a = inverse_roots(A, precision)
    lower_b, upper_b = inverse_roots(B, precision)
    nt = (arb(T) / (2 * arb.pi())).sqrt()
    require(arb(39852) < lower_a < arb(39853), "lower A root classification failed")
    require(arb(39936) < upper_a < arb(39937), "upper A root classification failed")
    require(arb(621) < lower_b < arb(622), "lower B root classification failed")
    require(arb(2560588) < upper_b < arb(2560589), "upper B root classification failed")
    require(arb(39894) < nt < arb(39895), "N_t classification failed")
    lower_interior_count = 39852 - 622 + 1
    upper_interior_count = 2560588 - 39937 + 1
    lower_turning_count = 39894 - 39853 + 1
    upper_turning_count = 39936 - 39895 + 1
    require(lower_turning_count == upper_turning_count == 42, "turning split drift")
    return {
        "precision_decimal_digits": precision,
        "lower_root_at_first_alpha_ball": lower_a.str(precision, more=True),
        "upper_root_at_first_alpha_ball": upper_a.str(precision, more=True),
        "lower_root_at_last_alpha_ball": lower_b.str(precision, more=True),
        "upper_root_at_last_alpha_ball": upper_b.str(precision, more=True),
        "riemann_siegel_nt_ball": nt.str(precision, more=True),
        "lower_branch_interior_modes": "622..39852",
        "lower_branch_interior_count": lower_interior_count,
        "lower_turning_endpoint_modes": "39853..39894",
        "lower_turning_endpoint_count": lower_turning_count,
        "upper_turning_endpoint_modes": "39895..39936",
        "upper_turning_endpoint_count": upper_turning_count,
        "upper_branch_interior_modes": "39937..2560588",
        "upper_branch_interior_count": upper_interior_count,
        "turning_gap_modes": "39853..39936",
        "turning_gap_count": lower_turning_count + upper_turning_count,
        "classical_source_upper_range": "622..39894",
        "classical_source_upper_term_count": 39894 - 622 + 1,
        "classical_interior_term_count": lower_interior_count,
        "classical_lower_endpoint_fresnel_term_count": lower_turning_count,
        "full_two_variable_saddle_nondegenerate": True,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["stationary_classification"]
    s = artifact["symbolic_reduction"]
    return f"""# Finite-Poisson portcullis saddle reduction

Date: 2026-08-10

Status: exact finite-Poisson and joint-saddle reduction validated; not a proof of the uniform remainder

Insert the exact incomplete-theta derivative from the equation-(10) current
gate into the equation-(9) Kummer integral.  Finite Poisson summation on the
contiguous odd-alpha roster gives dual mode `m` with joint phase

```text
Phi(alpha,x;m,t)
 =pi*alpha^2*x/4-pi*m*alpha+(t/2)log((1-x)/x).
```

The two stationary equations have the unique positive solution

```text
x_m     =2*pi*m^2/(t+2*pi*m^2),
alpha_m =2*m+t/(pi*m).                                  (P1)
```

Equation (P1) is exactly paper equation (124) with `m=N`.  This derives the
appearance of `pi`: it comes from the quadratic Kummer phase and the Fourier
dual frequency, not from a fitted circle or polygon.

The full Hessian is nondegenerate:

```text
det Hess(Phi)= {s['hessian_determinant']} < 0.           (P2)
```

Thus the turning point of the one-dimensional map `alpha(N)` is a projection
effect.  The global two-variable saddle itself does not degenerate.

At `t=10^10`, interval arithmetic certifies

```text
N_-(159577) in (39852,39853),
N_+(159577) in (39936,39937),
N_t            in (39894,39895),
N_-(5122421) in (621,622),
N_+(5122421) in (2560588,2560589).
```

The interior lower branch therefore contains modes
`{c['lower_branch_interior_modes']}`.  The lower-alpha boundary gap has 84
modes: `{c['lower_turning_endpoint_modes']}` below `N_t` and
`{c['upper_turning_endpoint_modes']}` above it, exactly 42 on each side.  The
classical source upper range `622..39894` consists of {c['classical_interior_term_count']}
interior modes plus those 42 lower endpoint-Fresnel modes.

For each Poisson mode, the alpha integral has the exact grouped form

```text
J_m=(-1)^m/x*{{[E_m(B)-E_m(A)]/(i*pi)+m*F_m}}.
```

The boundary current and Fresnel integral must remain paired for each `m`.
Their two series need not converge separately.  Likewise, the continuous
reflection `N*=t/(2*pi*N)` satisfies `x_(N*)=1-x_N`, but it does not preserve
the integer mode lattice in general.  The reflected `x` branches must be
recombined by the exact finite-Poisson/theta formula, not by rounding a
reciprocal partner.

The next analytic target is now precise: prove an endpoint-uniform finite
Poisson/stationary-phase theorem for the weighted Kummer integral, retain the
two reflected branches and all 84 turning modes, and bound symmetric
nonstationary mode tails.  Then compare the resulting exact dual expression
with `T_upper` and only afterward insert the source cubic-saddle truncation and
Gaussian evaluator defects.

Proof boundary: exact algebra, symmetric finite-Poisson identity under its
stated smoothness assumptions, and finite endpoint classification only.  No
interchange theorem for the Kummer endpoint singularities, explicit stationary
remainder, source-error bound, `Lambda<=0`, RH, or prize-level conclusion is
proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    current = load_json(CURRENT_GATE)
    ledger = load_json(SOURCE_LEDGER)
    require(current.get("passed") is True, "current gate is not passed")
    require(ledger.get("passed") is True, "source ledger is not passed")
    require(current["aggregate"]["first_alpha"] == A and current["aggregate"]["last_alpha"] == B, "alpha endpoints drift")
    require(ledger["aggregate"]["classical_upper_partition_term_count"] == 39273, "classical term count drift")

    symbolic = symbolic_reduction()
    classification = stationary_classification(PRECISION)
    require(classification["classical_source_upper_term_count"] == 39273, "stationary/classical roster mismatch")
    require(
        classification["classical_interior_term_count"] + classification["classical_lower_endpoint_fresnel_term_count"] == 39273,
        "classical interior plus endpoint roster does not close",
    )
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_finite_poisson_portcullis_saddle_reduction_gate",
        "status": "exact_finite_poisson_joint_saddle_derives_equation124_and_resolves_projected_turning_geometry",
        "passed": True,
        "scope": {
            "height_center": "1e10",
            "alpha_roster": f"{A}..{B} odd",
            "poisson_mode": "m integer",
            "diagnostic_midpoint_used": False,
        },
        "symbolic_reduction": symbolic,
        "stationary_classification": classification,
        "decision": {
            "equation124_is_joint_poisson_kummer_saddle_map": True,
            "pi_provenance_is_explicit": True,
            "projected_alpha_turning_point_is_full_saddle_degeneracy": False,
            "lower_endpoint_requires_42_classical_fresnel_modes": True,
            "reflected_continuous_saddle_involution_is_exact": True,
            "reflected_involution_preserves_integer_lattice": False,
            "mode_boundary_and_fresnel_series_may_be_split_before_symmetric_summation": False,
            "endpoint_complete_two_branch_route_is_viable": True,
            "height_uniform_remainder_proved": False,
            "rh_implication": False,
        },
        "next_obligation": "Prove an endpoint-uniform finite-Poisson/stationary-phase interchange for the equation-(9) Kummer weight, keep each mode's boundary current paired with its Fresnel integral, recombine x and 1-x branches on the integer lattice, explicitly cover the 84 turning modes, and bound symmetric nonstationary tails before comparing with T_upper.",
        "proof_boundary": "Exact finite-Poisson algebra and finite endpoint classification only. No endpoint-singular interchange theorem, explicit stationary remainder, source-aligned height-uniform Hardy-error bound, Lambda<=0, RH, or prize-level conclusion is proved.",
        "dependencies": {
            "equation10_current_gate": {"path": relative(CURRENT_GATE), "sha256": file_hash(CURRENT_GATE)},
            "source_aligned_ledger": {"path": relative(SOURCE_LEDGER), "sha256": file_hash(SOURCE_LEDGER)},
            "paper": {"path": relative(PAPER), "sha256": file_hash(PAPER)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "flint_threads": 1,
            "sympy_threads": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built finite-Poisson portcullis saddle reduction: "
        f"interior={classification['classical_interior_term_count']}, turning={classification['turning_gap_count']}, endpoint-classical={classification['classical_lower_endpoint_fresnel_term_count']}"
    )


if __name__ == "__main__":
    main()
