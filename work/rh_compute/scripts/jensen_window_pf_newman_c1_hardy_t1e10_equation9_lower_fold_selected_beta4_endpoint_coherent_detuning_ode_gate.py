#!/usr/bin/env python3
"""Derive the endpoint-coherent detuning ODE for the selected beta^-4 branch."""

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

from flint import acb, arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_beta4_endpoint_coherent_detuning_ode_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "beta4_reduction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_canonical_airy_derivative_reduction_gate.json",
    "completed_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.json",
    "event_pairing": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_event_ordered_398_pairing_multiplier_gate.json",
}

C = 159_577
Q = (C - 1) // 4
PAIR_COUNT = 198
PRECISION = 90


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


def complex_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": value.real.str(PRECISION, more=True),
        "imag_ball": value.imag.str(PRECISION, more=True),
    }


def symbolic_certificate() -> dict[str, str]:
    lam, y, epsilon = sp.symbols("lambda y epsilon", real=True)
    x = lam + y
    u = 1 + epsilon * (-(13 * lam + 3 * y) / 60)
    u += epsilon**2 * -(
        448 * lam**5 + 280 * lam**2 * y**3 + 4565 * lam**2
        - 105 * lam * y**4 + 30 * lam * y + 63 * y**5 - 405 * y**2
    ) / 50400
    v = epsilon * (8 * lam**2 - 4 * lam * y + 3 * y**2) / 60
    v += epsilon**2 * (40 * lam**3 - 20 * lam**2 * y + lam * y**2 - 9 * y**3 + 27) / 1680
    forcing_w = sp.factor(sp.diff(u, y, 2) - v - 2 * x * sp.diff(v, y))
    forcing_wx = sp.factor(2 * sp.diff(u, y) + sp.diff(v, y, 2))
    expected_w = -epsilon * y**2 / 4 + epsilon**2 * (13 * lam * y**2 / 240 + y**3 / 80)
    expected_wx = -epsilon**2 * y**2 * (8 * lam**2 - 4 * lam * y + 3 * y**2) / 240
    require(sp.expand(forcing_w - expected_w) == 0, "W forcing identity failed")
    require(sp.expand(forcing_wx - expected_wx) == 0, "W_X forcing identity failed")

    h, a, b, ap, bp = sp.symbols("h a b ap bp", real=True)
    i = sp.I
    f_minus = a + i * b
    f_plus = a - i * b
    fp_minus = ap + i * bp
    fp_plus = ap - i * bp
    lower = 0
    upper_inner = 0
    for j in range(1, PAIR_COUNT + 1):
        d_minus = -h * (sp.Rational(1, 4) + j)
        lower += fp_minus + i * d_minus * f_minus
        upper_inner += i * (h / 2 - d_minus) * f_minus - fp_minus
    for j in range(1, PAIR_COUNT + 3):
        d_plus = h * (j - sp.Rational(1, 4))
        lower += fp_plus + i * d_plus * f_plus
        upper_inner += i * (h / 2 - d_plus) * f_plus - fp_plus
    lower_closed = 398 * ap + sp.Rational(79601, 2) * h * b
    lower_closed += i * (-2 * bp + sp.Rational(599, 2) * h * a)
    upper_closed = -398 * ap - sp.Rational(79599, 2) * h * b
    upper_closed += i * (2 * bp - sp.Rational(201, 2) * h * a)
    require(sp.expand(lower - lower_closed) == 0, "lower endpoint source sum failed")
    require(sp.expand(upper_inner - upper_closed) == 0, "upper endpoint source sum failed")

    return {
        "branch": "C_sigma=(U W_sigma+V W_sigma,X)/2; W_sigma,XX=-X W_sigma",
        "finite_transform": "G_sigma(d)=Integral_0^Y exp(i[y^2/(4beta)-dy]) C_sigma(lambda,y)dy",
        "detuning_ode": (
            "G_sigma''/(4beta^2)+i(1+d/beta)G_sigma'"
            "+(lambda-d^2+i/(2beta))G_sigma=R_0+R_Y+Q_sigma"
        ),
        "lower_source": "R_0=C_sigma,y(lambda,0)+i d C_sigma(lambda,0)",
        "upper_source": (
            "R_Y=exp(i[Y^2/(4beta)-dY])"
            "{i[Y/(2beta)-d]C_sigma(lambda,Y)-C_sigma,y(lambda,Y)}"
        ),
        "forcing_W": str(sp.factor(expected_w)),
        "forcing_WX": str(sp.factor(expected_wx)),
        "forcing_integral": "Q_sigma=Integral exp(i[y^2/(4beta)-dy])*(E_0 W_sigma+E_1 W_sigma,X)/2 dy",
        "remaining_lower_source": "398a'+(79601/2)h b+i[-2b'+(599/2)h a]",
        "remaining_upper_source_before_common_phase": "-398a'-(79599/2)h b+i[2b'-(201/2)h a]",
        "remaining_forcing_projection": (
            "2 exp(ix) sin(396x)/sin(2x) Re(Q_- exp(i398x))"
            "+2 Q_+ exp(-i797x)cos(2x)"
        ),
    }


def weights_and_derivatives(lam: arb, y: arb, beta: arb) -> tuple[arb, arb, arb, arb]:
    epsilon = 1 / beta**2
    u1 = -(13 * lam + 3 * y) / 60
    v1 = (8 * lam**2 - 4 * lam * y + 3 * y**2) / 60
    n_u2 = (
        448 * lam**5 + 280 * lam**2 * y**3 + 4565 * lam**2
        - 105 * lam * y**4 + 30 * lam * y + 63 * y**5 - 405 * y**2
    )
    n_v2 = 40 * lam**3 - 20 * lam**2 * y + lam * y**2 - 9 * y**3 + 27
    u = 1 + epsilon * u1 - epsilon**2 * n_u2 / 50400
    v = epsilon * v1 + epsilon**2 * n_v2 / 1680
    n_u2_y = 840 * lam**2 * y**2 - 420 * lam * y**3 + 30 * lam + 315 * y**4 - 810 * y
    n_v2_y = -20 * lam**2 + 2 * lam * y - 27 * y**2
    u_y = -epsilon / 20 - epsilon**2 * n_u2_y / 50400
    v_y = epsilon * (-4 * lam + 6 * y) / 60 + epsilon**2 * n_v2_y / 1680
    return u, v, u_y, v_y


def branch_and_derivative(lam: arb, y: arb, beta: arb, sigma: int) -> tuple[acb, acb]:
    x = lam + y
    ai, aip, bi, bip = (-x).airy()
    w = acb(ai, -sigma * bi)
    w_x = acb(-aip, sigma * bip)
    u, v, u_y, v_y = weights_and_derivatives(lam, y, beta)
    value = (u * w + v * w_x) / 2
    derivative = ((u_y - x * v) * w + (u + v_y) * w_x) / 2
    return value, derivative


def interval_certificate() -> dict[str, Any]:
    pi = arb.pi()
    c = arb(C)
    beta = (pi * c**2 / 8) ** (arb(1) / 3)
    h = 4 * beta / c
    y_max = pi * c / (2 * beta)
    lambda_max = pi / (16 * beta)
    lam = arb(lambda_max / 2, lambda_max / 2)
    y = arb(y_max / 2, y_max / 2)
    epsilon = 1 / beta**2
    forcing_w = -epsilon * y**2 / 4 + epsilon**2 * (13 * lam * y**2 / 240 + y**3 / 80)
    forcing_wx = -epsilon**2 * y**2 * (8 * lam**2 - 4 * lam * y + 3 * y**2) / 240
    require(abs(forcing_w) < arb("7.31e-4"), "W forcing coefficient exceeds 7.31e-4")
    require(abs(forcing_wx) < arb("1.08e-7"), "W_X forcing coefficient exceeds 1.08e-7")

    lower_value, lower_derivative = branch_and_derivative(lam, arb(0), beta, -1)
    upper_value, upper_derivative = branch_and_derivative(lam, y_max, beta, -1)

    def source(value: acb, derivative: acb, upper: bool) -> acb:
        a, b = value.real, value.imag
        ap, bp = derivative.real, derivative.imag
        if not upper:
            return acb(
                398 * ap + arb(79601) * h * b / 2,
                -2 * bp + arb(599) * h * a / 2,
            )
        inner = acb(
            -398 * ap - arb(79599) * h * b / 2,
            2 * bp - arb(201) * h * a / 2,
        )
        common_phase = acb(0, y_max**2 / (4 * beta) + pi / 2).exp()
        return common_phase * inner

    lower_source = source(lower_value, lower_derivative, False)
    upper_source = source(upper_value, upper_derivative, True)
    require((h * y_max - 2 * pi).contains(0), "hY=2pi identity lost")
    require((y_max / (2 * beta) - h / 2).contains(0), "upper derivative scale identity lost")
    return {
        "height_interval": "t*-pi/16<=t<=t*",
        "lambda_interval": "0<=lambda<=pi/(16beta)",
        "beta_ball": beta.str(PRECISION, more=True),
        "h_ball": h.str(PRECISION, more=True),
        "Y_ball": y_max.str(PRECISION, more=True),
        "lambda_max_ball": lambda_max.str(PRECISION, more=True),
        "hY_identity": "hY=2pi",
        "upper_endpoint_phase": "exp(-i d_m Y)=i for every selected mode because C=1 mod 4",
        "maximum_forcing_W_coefficient_ball": abs(forcing_w).str(PRECISION, more=True),
        "maximum_forcing_WX_coefficient_ball": abs(forcing_wx).str(PRECISION, more=True),
        "remaining_398_lower_endpoint_source_ball": complex_record(lower_source),
        "remaining_398_lower_endpoint_source_absolute_ball": abs(lower_source).str(PRECISION, more=True),
        "remaining_398_upper_endpoint_source_ball": complex_record(upper_source),
        "remaining_398_upper_endpoint_source_absolute_ball": abs(upper_source).str(PRECISION, more=True),
        "lower_branch_value_ball": complex_record(lower_value),
        "lower_branch_derivative_ball": complex_record(lower_derivative),
        "upper_branch_value_ball": complex_record(upper_value),
        "upper_branch_derivative_ball": complex_record(upper_derivative),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    return f"""# Endpoint-coherent detuning ODE for the selected beta-minus-four branch

Date: 2026-08-13

Status: exact finite-transform ODE, endpoint-source compression, and grouped
forcing structure proved; no quantitative finite-integral splice theorem

Let `X=lambda+y`, `epsilon=beta^-2`, and

```text
C_sigma=(U W_sigma+V W_sigma,X)/2,
W_sigma,XX=-X W_sigma.
```

For the finite selected-branch transform

```text
G_sigma(d)=Integral_0^Y exp(i[y^2/(4beta)-d y])
                         C_sigma(lambda,y)dy,          (ED1)
```

two integrations by parts give the exact detuning equation

```text
G_sigma''/(4beta^2)+i(1+d/beta)G_sigma'
 +(lambda-d^2+i/(2beta))G_sigma
 =R_0,sigma+R_Y,sigma+Q_sigma,                        (ED2)

R_0,sigma=C_sigma,y(lambda,0)+i d C_sigma(lambda,0), (ED3)

R_Y,sigma=e^(i[Y^2/(4beta)-dY])
 [i(Y/(2beta)-d)C_sigma(lambda,Y)-C_sigma,y(lambda,Y)]. (ED4)
```

The beta-minus-four correction does not leave an opaque bulk error.  Exact
Airy reduction gives

```text
Q_sigma=Integral_0^Y e^(i[y^2/(4beta)-dy])
                    [E_0 W_sigma+E_1 W_sigma,X]/2 dy,

E_0=-epsilon*y^2/4
    +epsilon^2(13lambda*y^2/240+y^3/80),
E_1=-epsilon^2*y^2(8lambda^2-4lambda*y+3y^2)/240.    (ED5)
```

Both forcing channels vanish quadratically at the lower endpoint.  On the
whole selector top corridor and `0<=y<=Y`, interval arithmetic proves

```text
|E_0| <= {c['maximum_forcing_W_coefficient_ball']} <7.31e-4,
|E_1| <= {c['maximum_forcing_WX_coefficient_ball']} <1.08e-7. (ED6)
```

Now remove event zero and put `h=4beta/C`, `d_0=-h/4`.  The remaining roster
is `d_0-hj`, `1<=j<=198`, and `d_0+hj`, `1<=j<=200`.  If
`C_-=a+ib` and `C_-,y=a'+ib'`, direct finite summation gives

```text
sum R_0=398a'+(79601/2)h b
        +i[-2b'+(599/2)h a],                          (ED7)

sum e^(-i common phase)R_Y
 =-398a'-(79599/2)h b
   +i[2b'-(201/2)h a].                               (ED8)
```

Here (ED7) uses the values at `y=0` and (ED8) the values at `y=Y`.  The
upper phase is genuinely common: `hY=2pi`, `C=1 mod 4`, and therefore
`exp(-i d_m Y)=i` for every selected mode.  The interval source aggregates
are

```text
lower: {c['remaining_398_lower_endpoint_source_absolute_ball']},
upper: {c['remaining_398_upper_endpoint_source_absolute_ball']}. (ED9)
```

They are coherent boundary channels, not 398 unrelated errors.  The interior
forcing also retains the exact event-ordered projection (11.385.4), with
`C_sigma` replaced by `(E_0 W_sigma+E_1 W_sigma,X)/2`.  Thus no termwise
absolute value is required anywhere in (ED2)--(ED8).

This explains structurally why a bare full-saddle comparison can fail even
when individual interior modes look accurate: it omits coherent finite
endpoint currents.  Equations (ED2)--(ED8) do not yet prove that those
currents cancel the completed Poisson complement, nor do they bound the
variation-of-constants kernel.

Pi provenance: `Y=pi C/(2beta)` and `beta^3=pi C^2/8` come from the exact
Kummer/Fourier selector geometry.  They force `hY=2pi` and
`Y/(2beta)=h/2`; no independent circle or fitted period is inserted.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No completed endpoint cancellation, uniform finite-integral kernel bound,
grouped 398-mode splice, all-corridor continuation, complete `Q_K-T` or
`T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    require(dependencies["beta4_reduction"].get("passed") is True, "beta4 reduction dependency failed")
    require(dependencies["completed_projection"].get("passed") is True, "completed projection dependency failed")
    require(dependencies["event_pairing"].get("passed") is True, "event pairing dependency failed")

    ctx.dps = PRECISION
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "selected_beta4_finite_transform_endpoint_sources_and_forcing_grouped_exactly",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(),
        "decision": {
            "corrected_selected_branch_detuning_ODE_proved": True,
            "remaining_398_lower_endpoint_sources_compressed_exactly": True,
            "remaining_398_upper_endpoint_sources_compressed_exactly": True,
            "upper_endpoint_lattice_phase_is_common": True,
            "interior_forcing_retains_event_ordered_Dirichlet_projection": True,
            "forcing_coefficients_uniformly_bounded": True,
            "completed_endpoint_cancellation_proved": False,
            "finite_integral_kernel_bound_proved": False,
            "grouped_398_mode_splice_proved": False,
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
            "Construct the variation-of-constants Green kernel for (ED2) in the centered detuning d=d0+hj. "
            "Compose the two coherent endpoint-source aggregates with g4_lambda(0), the endpoint half-current, and the "
            "zero/negative/outer-positive completion before bounding the grouped forcing projection."
        ),
        "proof_boundary": (
            "Exact corrected-branch finite-transform ODE, endpoint-source compression, common upper lattice phase, and "
            "uniform forcing-coefficient bounds only. No completed endpoint cancellation, finite-integral kernel bound, "
            "grouped 398-mode splice, complete Q_K-T or T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("derived endpoint-coherent beta^-4 detuning ODE for the remaining 398 selected terms", flush=True)


if __name__ == "__main__":
    main()
