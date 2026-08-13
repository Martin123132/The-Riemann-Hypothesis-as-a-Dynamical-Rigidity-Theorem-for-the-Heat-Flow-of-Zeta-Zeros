#!/usr/bin/env python3
"""Derive the exact selected-fold/logistic stationary-action bridge."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_logistic_stationary_action_bridge_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "beta4_completed_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.json",
    "event_selection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_hankel_carrier_event_selection_gate.json",
    "exact_logistic_chart": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_full_line_join_gate.json",
    "gamma_bulk": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_global_morse_gamma_bulk_gate.json",
    "ownership_ledger": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_ordinary_mode_coverage_ledger_gate.json",
}

C = 159_577
MODE_LO = 39_696
MODE_HI = 40_094
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


def symbolic_certificate() -> dict[str, Any]:
    s, mu = sp.symbols("s mu", positive=True)
    r = (1 - mu) / s**2
    q = sp.log(r) / 2
    tanh_q = (r - 1) / (r + 1)
    x = (1 - tanh_q) / 2
    y_exact = sp.simplify(((s - 1) + tanh_q) / x)
    theta_exact = sp.simplify((1 - mu) * q - tanh_q - x * y_exact**2 / 2)

    z = 1 - sp.sqrt(2 * s + mu - 1)
    y_canonical = sp.simplify(z**2 - mu)
    theta_canonical = sp.simplify(z**3 / 3 - mu * z - y_canonical**2 / 4)
    defect = sp.factor(sp.simplify(theta_exact - theta_canonical))
    expected = (
        3 * s**2
        - 8 * s * sp.sqrt(2 * s + mu - 1)
        + 12 * s
        - 4 * mu * sp.sqrt(2 * s + mu - 1)
        + 6 * mu * sp.log(s)
        - 3 * mu * sp.log(1 - mu)
        + 9 * mu
        + 4 * sp.sqrt(2 * s + mu - 1)
        - 6 * sp.log(s)
        + 3 * sp.log(1 - mu)
        - 11
    ) / 6
    require(sp.simplify(defect - expected) == 0, "stationary-action defect simplification failed")
    require(sp.factor(y_exact - ((s - 1) ** 2 - mu) / s) == 0, "exact stationary y simplification failed")
    require(sp.factor(y_canonical - 2 * (s - sp.sqrt(2 * s + mu - 1))) == 0, "canonical stationary y simplification failed")

    center = sp.simplify(defect.subs(mu, 0))
    center_derivatives = [sp.simplify(sp.diff(center, s, order).subs(s, 1)) for order in range(7)]
    require(center_derivatives[:5] == [0] * 5, "center action lost fifth-order tangency")
    require(center_derivatives[5] == 6 and center_derivatives[6] == -90, "center tangency coefficients drift")
    height_derivative = sp.simplify(sp.diff(defect, mu).subs(mu, 0))
    height_derivatives = [sp.simplify(sp.diff(height_derivative, s, order).subs(s, 1)) for order in range(6)]
    require(height_derivatives[:3] == [0] * 3, "height derivative lost cubic tangency")
    require(height_derivatives[3:] == [-1, 9, -81], "height tangency coefficients drift")

    q_symbol, n_symbol = sp.symbols("q n", real=True)
    tanh_symbol = sp.tanh(q_symbol)
    x_symbol = (1 - tanh_symbol) / 2
    exact_phase = (
        (1 - mu) * q_symbol
        - tanh_symbol
        - ((s - 1) + tanh_symbol) * n_symbol
        + x_symbol * n_symbol**2 / 2
    )
    exact_hessian = sp.hessian(exact_phase, (q_symbol, n_symbol))
    exact_det = sp.factor(
        sp.simplify(exact_hessian.det().subs({q_symbol: q, n_symbol: y_exact}).rewrite(sp.exp))
    )
    canonical_hessian = sp.Matrix([[2 * z, -1], [-1, sp.Rational(1, 2)]])
    canonical_det = sp.factor(canonical_hessian.det())
    expected_exact_det = -2 * s**2 * (1 - mu) / (s**2 + 1 - mu)
    expected_canonical_det = -sp.sqrt(2 * s + mu - 1)
    require(sp.simplify(exact_det - expected_exact_det) == 0, "exact stationary Hessian determinant drift")
    require(sp.simplify(canonical_det - expected_canonical_det) == 0, "canonical stationary Hessian determinant drift")

    exact_scalar_amplitude = sp.cosh(q_symbol) ** (-sp.Rational(3, 2)) * (1 + n_symbol / 2)
    exact_scalar_at_saddle = sp.factor(
        sp.simplify(exact_scalar_amplitude.subs({q_symbol: q, n_symbol: y_exact}).rewrite(sp.exp))
    )
    stationary_amplitude_ratio = sp.factor(
        sp.simplify(exact_scalar_at_saddle * sp.sqrt(expected_canonical_det / expected_exact_det))
    )
    expected_amplitude_ratio = ((1 - mu) * (2 * s + mu - 1)) ** sp.Rational(1, 4) / sp.sqrt(s)
    require(
        sp.simplify(stationary_amplitude_ratio**4 - expected_amplitude_ratio**4) == 0,
        "stationary amplitude ratio fourth-power drift",
    )
    center_amplitude = sp.simplify(expected_amplitude_ratio.subs(mu, 0))
    amplitude_derivatives = [sp.simplify(sp.diff(center_amplitude, s, order).subs(s, 1)) for order in range(5)]
    require(
        amplitude_derivatives == [1, 0, -sp.Rational(1, 2), 3, -sp.Rational(81, 4)],
        "stationary amplitude tangency drift",
    )

    return {
        "scaled_variables": "s=4m/C, mu=(t*-t)/t*=lambda/beta^2, beta^3=t*",
        "exact_logistic_stationary_coordinate": "q_e=(1/2)log((1-mu)/s^2)",
        "exact_normal_stationary_coordinate": "y_e/beta^2=((s-1)^2-mu)/s",
        "canonical_selected_stationary_coordinate": "z_c/beta=1-sqrt(2s+mu-1)",
        "canonical_normal_stationary_coordinate": "y_c/beta^2=2(s-sqrt(2s+mu-1))",
        "exact_stationary_action": "Theta_e/beta^3=(1-mu)q_e-tanh(q_e)-x_e*(y_e/beta^2)^2/2",
        "canonical_stationary_action": "Theta_c/beta^3=w^3/3-mu*w-(y_c/beta^2)^2/4, w=z_c/beta",
        "action_defect": str(expected),
        "exact_bridge": "Theta_e-Theta_c=beta^3*F(s,mu); exact logistic carrier=exp(i*beta^3*F)*canonical selected carrier",
        "center_tangency": {
            "derivatives_order_0_through_6": [str(value) for value in center_derivatives],
            "series_lead": "F(s,0)=(s-1)^5/20-(s-1)^6/8+13(s-1)^7/56+...",
        },
        "height_tangency": {
            "d_dmu_derivatives_order_0_through_5": [str(value) for value in height_derivatives],
            "series_lead": "F_mu(s,0)=-(s-1)^3/6+3(s-1)^4/8-27(s-1)^5/40+...",
        },
        "exact_stationary_Hessian_determinant": str(expected_exact_det),
        "canonical_stationary_Hessian_determinant": str(expected_canonical_det),
        "exact_scalar_amplitude_at_saddle": str(exact_scalar_at_saddle),
        "stationary_amplitude_ratio": str(expected_amplitude_ratio),
        "stationary_amplitude_center_tangency": {
            "derivatives_order_0_through_4": [str(value) for value in amplitude_derivatives],
            "series_lead": "R_A(s,0)=1-(s-1)^2/4+(s-1)^3/2-27(s-1)^4/32+...",
        },
        "combined_leading_carrier_bridge": "L_exact/L_canonical=R_A(s,mu)*exp(i*beta^3*F(s,mu))",
    }


def f_center(s: arb) -> arb:
    root = (2 * s - 1).sqrt()
    return (3 * s * s - 8 * s * root + 12 * s + 4 * root - 6 * s.log() - 11) / 6


def f_mu(s: arb, mu: arb) -> arb:
    root = (2 * s + mu - 1).sqrt()
    numerator = 4 * s + 2 * mu - 2 * root * s.log() + root * (1 - mu).log() - 2 * root - 2
    return -numerator / (2 * root)


def row_certificate(mode: int, c: arb, beta: arb, mu_radius: arb) -> dict[str, Any]:
    s = 4 * arb(mode) / c
    mu = arb(mu_radius / 2, mu_radius / 2)
    center_action = beta**3 * f_center(s)
    derivative_enclosure = f_mu(s, mu)
    variation = beta**3 * mu_radius * abs(derivative_enclosure)
    action_bound = abs(center_action) + variation
    amplitude_ratio = ((1 - mu) * (2 * s + mu - 1)) ** (arb(1) / 4) / s.sqrt()
    amplitude_ratio_defect = abs(amplitude_ratio - 1)
    combined_leading_defect = amplitude_ratio_defect + amplitude_ratio.upper() * action_bound
    carrier_scaled_defect = combined_leading_defect / arb(mode).sqrt()
    physical_carrier_scaled_defect = arb(2).sqrt() / arb.pi() * carrier_scaled_defect
    delta = s - 1
    exact_y = beta * beta * (delta * delta - mu) / s
    canonical_y = 2 * beta * beta * (s - (2 * s + mu - 1).sqrt())
    require(exact_y.lower() > 0 and canonical_y.lower() > 0, f"mode {mode} leaves the positive normal saddle")
    return {
        "mode": mode,
        "detuning_label": 4 * mode - C,
        "s_ball": s.str(PRECISION, more=True),
        "center_action_defect_ball": center_action.str(PRECISION, more=True),
        "height_derivative_F_mu_ball": derivative_enclosure.str(PRECISION, more=True),
        "height_variation_bound_ball": variation.str(PRECISION, more=True),
        "uniform_action_defect_absolute_bound_ball": action_bound.str(PRECISION, more=True),
        "uniform_carrier_multiplier_defect_bound_ball": action_bound.str(PRECISION, more=True),
        "stationary_amplitude_ratio_ball": amplitude_ratio.str(PRECISION, more=True),
        "stationary_amplitude_ratio_defect_bound_ball": amplitude_ratio_defect.str(PRECISION, more=True),
        "combined_leading_carrier_defect_bound_ball": combined_leading_defect.str(PRECISION, more=True),
        "carrier_scaled_combined_leading_defect_bound_ball": carrier_scaled_defect.str(PRECISION, more=True),
        "physical_carrier_scaled_combined_leading_defect_bound_ball": physical_carrier_scaled_defect.str(PRECISION, more=True),
        "exact_normal_saddle_y_ball": exact_y.str(PRECISION, more=True),
        "canonical_normal_saddle_y_ball": canonical_y.str(PRECISION, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Selected fold/logistic stationary-action bridge

Date: 2026-08-13

Status: exact all-mode top-corridor stationary carrier map certified; the
uniform finite-integral amplitude match, complete ordinary splice, `T_upper`,
and RH remain open

Put

```text
s=4m/C,  mu=(t*-t)/t*=lambda/beta^2,  beta^3=t*.       (SA1)
```

Solving the exact logistic and selected canonical joint stationary equations
gives

```text
q_e=(1/2)log((1-mu)/s^2),
y_e/beta^2=((s-1)^2-mu)/s,
z_c/beta=1-sqrt(2s+mu-1),
y_c/beta^2=2[s-sqrt(2s+mu-1)].                         (SA2)
```

Substitution into the two stationary actions simplifies exactly to

```text
Theta_e-Theta_c=beta^3 F(s,mu),                        (SA3)

6F=3s^2-8sR+12s-4mu R+6mu log s-3mu log(1-mu)
   +9mu+4R-6log s+3log(1-mu)-11,
R=sqrt(2s+mu-1).                                       (SA4)
```

This is the finite-height carrier bridge: the exact logistic carrier is the
selected canonical carrier multiplied by `exp(i*beta^3*F)`.  It is not an
asymptotic phase declaration.

At the selector center, symbolic differentiation proves

```text
F(s,0)=(s-1)^5/20-(s-1)^6/8+13(s-1)^7/56+...;          (SA5)
```

the derivatives through order four vanish and the fifth derivative is `6`.
The first height derivative has the lower-order but still rigid tangency

```text
F_mu(s,0)=-(s-1)^3/6+3(s-1)^4/8-27(s-1)^5/40+... .    (SA6)
```

All `{c['mode_count']}` modes `m={MODE_LO},...,{MODE_HI}` were enclosed on
the full top corridor `t*-pi/16<=t<=t*`.  Both stationary normal coordinates
remain inside the one-cell strip, with

```text
0 < y_e <= {c['maximum_exact_normal_saddle_y_upper_ball']},
0 < y_c <= {c['maximum_canonical_normal_saddle_y_upper_ball']},
Y={c['strip_width_Y_ball']}.                            (SA7)
```

The uniform action and multiplier bounds are

```text
max_m sup_t |Theta_e-Theta_c|
 < {c['maximum_uniform_action_defect_bound_ball']} < 0.001555,
|exp(i[Theta_e-Theta_c])-1| <= |Theta_e-Theta_c|.       (SA8)
```

The maximum occurs at mode `{c['maximum_mode']}`.  Its center defect is
`{c['maximum_mode_center_action_defect_ball']}` and is rigorously nonzero, so
the multiplier cannot be replaced by an identity.  It must remain inside the
199-pair projection before summation or absolute values.

The exact scalar amplitude and the two Hessian determinants simplify at the
same saddles.  Their leading stationary-coefficient ratio is

```text
R_A(s,mu)=[(1-mu)(2s+mu-1)]^(1/4)/sqrt(s),             (SA9)
L_exact/L_canonical=R_A exp(i*beta^3 F).               (SA10)
```

At `mu=0`, `R_A=1-(s-1)^2/4+(s-1)^3/2+...`.  Across the
full top corridor, the exact inequality

```text
|R_A exp(i*Delta)-1| <= |R_A-1|+R_A|Delta|
```

gives the carrier-scaled bounds

```text
max_m sup_t |R_A exp(i*Delta)-1|/sqrt(m)
 < {c['maximum_carrier_scaled_combined_leading_defect_bound_ball']} < 7.8e-6,
physical < {c['maximum_physical_carrier_scaled_combined_leading_defect_bound_ball']}.
                                                               (SA11)
```

This consumes less than the existing local normalized headroom
`{c['local_normalized_headroom_ball']}` and leaves at least
`{c['remaining_local_normalized_headroom_lower_ball']}` for the uniform
endpoint/Fresnel and amplitude-variation remainder.  These are local
one-carrier diagnostics, not a 399-term triangle estimate and not an aggregate
`T_upper` budget.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No uniform finite-integral selected-branch amplitude theorem, grouped 199-pair error bound, all-400-
corridor continuation, complete `Q_K-T` or `T_upper`, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    require(dependencies["beta4_completed_projection"].get("passed") is True, "beta4 completed projection failed")
    require(dependencies["event_selection"]["decision"]["selected_extracted_carrier_crosses_at_event_center"] is True, "selected carrier dependency drift")
    require(dependencies["gamma_bulk"]["decision"]["ordinary_full_line_bulk_closed_at_saved_height"] is True, "Gamma bulk dependency drift")

    ctx.dps = PRECISION
    algebra = symbolic_certificate()
    pi = arb.pi()
    c = arb(C)
    tstar = pi * c * c / 8
    beta = tstar ** (arb(1) / 3)
    lambda_radius = pi / (16 * beta)
    mu_radius = lambda_radius / (beta * beta)
    require((mu_radius - pi / (16 * tstar)).contains(0), "height scaling identity drift")
    y_width = pi * c / (2 * beta)

    rows = [row_certificate(mode, c, beta, mu_radius) for mode in range(MODE_LO, MODE_HI + 1)]
    require(len(rows) == 399, "shared-boundary mode count drift")
    maximum_row = max(rows, key=lambda row: arb(row["uniform_action_defect_absolute_bound_ball"]).upper())
    maximum_bound = arb(maximum_row["uniform_action_defect_absolute_bound_ball"])
    require(maximum_row["mode"] == MODE_HI, "maximum action defect mode drift")
    require(maximum_bound < arb("0.001555"), "top-corridor action defect exceeds 0.001555")
    require(abs(arb(maximum_row["center_action_defect_ball"])) > arb("0.001553"), "nontrivial phase defect guard failed")
    maximum_carrier_row = max(
        rows,
        key=lambda row: arb(row["carrier_scaled_combined_leading_defect_bound_ball"]).upper(),
    )
    maximum_carrier_scaled = arb(maximum_carrier_row["carrier_scaled_combined_leading_defect_bound_ball"])
    maximum_physical_carrier_scaled = arb(
        maximum_carrier_row["physical_carrier_scaled_combined_leading_defect_bound_ball"]
    )
    require(maximum_carrier_row["mode"] == MODE_HI, "maximum carrier-scaled mode drift")
    require(maximum_carrier_scaled < arb("7.8e-6"), "combined leading carrier defect exceeds 7.8e-6")
    ownership = dependencies["ownership_ledger"]
    local_headroom = arb(ownership["quantitative_frontier"]["local_normalized_headroom_ball"])
    remaining_headroom = local_headroom - maximum_carrier_scaled
    require(remaining_headroom.lower() > arb("9e-6"), "remaining local headroom fell below 9e-6")
    maximum_exact_y = max(arb(row["exact_normal_saddle_y_ball"]).upper() for row in rows)
    maximum_canonical_y = max(arb(row["canonical_normal_saddle_y_ball"]).upper() for row in rows)
    minimum_exact_y = min(arb(row["exact_normal_saddle_y_ball"]).lower() for row in rows)
    minimum_canonical_y = min(arb(row["canonical_normal_saddle_y_ball"]).lower() for row in rows)
    require(minimum_exact_y > arb("9e-5") and minimum_canonical_y > arb("9e-5"), "minimum top-corridor saddle clearance failed")
    require(maximum_exact_y < y_width and maximum_canonical_y < y_width, "stationary saddle leaves selector strip")

    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "exact_selected_fold_logistic_stationary_action_and_leading_amplitude_bridge_certified_on_top_corridor",
        "passed": True,
        "exact_algebra": algebra,
        "certificate": {
            "C": C,
            "mode_roster": [MODE_LO, MODE_HI],
            "mode_count": len(rows),
            "height_interval": "t*-pi/16 <= t <= t*",
            "lambda_interval": "0 <= lambda <= pi/(16beta)",
            "mu_interval": "0 <= mu <= pi/(16t*)",
            "mu_radius_ball": mu_radius.str(PRECISION, more=True),
            "strip_width_Y_ball": y_width.str(PRECISION, more=True),
            "minimum_exact_normal_saddle_y_lower_ball": minimum_exact_y.str(PRECISION, more=True),
            "minimum_canonical_normal_saddle_y_lower_ball": minimum_canonical_y.str(PRECISION, more=True),
            "maximum_exact_normal_saddle_y_upper_ball": maximum_exact_y.str(PRECISION, more=True),
            "maximum_canonical_normal_saddle_y_upper_ball": maximum_canonical_y.str(PRECISION, more=True),
            "maximum_mode": maximum_row["mode"],
            "maximum_mode_center_action_defect_ball": maximum_row["center_action_defect_ball"],
            "maximum_uniform_action_defect_bound_ball": maximum_bound.str(PRECISION, more=True),
            "carrier_multiplier_bound": "|exp(i*Delta)-1|<=|Delta|<0.001555",
            "maximum_carrier_scaled_mode": maximum_carrier_row["mode"],
            "maximum_carrier_scaled_combined_leading_defect_bound_ball": maximum_carrier_scaled.str(PRECISION, more=True),
            "maximum_physical_carrier_scaled_combined_leading_defect_bound_ball": maximum_physical_carrier_scaled.str(PRECISION, more=True),
            "local_normalized_headroom_ball": local_headroom.str(PRECISION, more=True),
            "remaining_local_normalized_headroom_lower_ball": remaining_headroom.lower().str(PRECISION, more=True),
            "rows": rows,
        },
        "decision": {
            "exact_logistic_and_canonical_joint_stationary_coordinates_solved": True,
            "exact_stationary_action_defect_formula_proved": True,
            "selector_center_action_has_fifth_order_tangency": True,
            "first_height_derivative_has_third_order_tangency": True,
            "all_399_top_corridor_stationary_saddles_remain_inside_strip": True,
            "uniform_action_defect_below_0_001555": True,
            "unit_phase_identification_is_not_exact": True,
            "phase_multiplier_must_remain_inside_grouped_projection": True,
            "exact_leading_stationary_amplitude_ratio_proved": True,
            "carrier_scaled_combined_leading_defect_below_7_8e_minus_6": True,
            "combined_leading_carrier_defect_fits_local_headroom": True,
            "selected_branch_amplitude_match_proved": False,
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
            "Retain R_A(s,mu)*exp(i*beta^3*F(s,mu)) inside each reflected selected pair. Derive a uniform endpoint/Fresnel remainder "
            "below the remaining 9e-6 local headroom, then bound finite differences across the 199-pair roster before any Dirichlet/Abel estimate."
        ),
        "proof_boundary": (
            "Exact joint-stationary coordinate, action, and leading Hessian-amplitude bridge on the 399-mode top corridor only. "
            "No uniform finite-integral endpoint/Fresnel remainder, grouped-pair error estimate, all-corridor continuation, complete Q_K-T or T_upper theorem, "
            "Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("built selected/logistic stationary-action bridge: 399 top-corridor modes, |Delta| < 0.001555", flush=True)


if __name__ == "__main__":
    main()
