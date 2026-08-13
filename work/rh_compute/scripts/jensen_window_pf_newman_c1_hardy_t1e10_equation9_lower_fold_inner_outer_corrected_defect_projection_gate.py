#!/usr/bin/env python3
"""Certify the inner/outer split of the corrected weighted fold defect."""

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

from flint import arb, acb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_inner_outer_corrected_defect_projection_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "selected_total": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_selected_branch_grouped_integral_gate.json",
    "event_pairing": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_event_ordered_398_pairing_multiplier_gate.json",
    "branch_envelope": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_scaled_airy_branch_amplitude_envelope_gate.json",
    "weighted_kernel": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_opposite_branch_dirichlet_gate.json",
    "exact_profile": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_kummer_ode_profile_closure_gate.json",
    "source_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_corrected_selected_source_projection_orientation_gate.json",
    "companion_decomposition": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_defect_raw_strip_companion_decomposition_gate.json",
}

C = 159_577
L = 39_696
Q = 39_894
U = 40_094
PANELS = 65_536
HEIGHT_NODES = 5
PRECISION = 70


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


def action_defect(s: arb, mu: arb) -> arb:
    root = (2 * s + mu - 1).sqrt()
    return (
        3 * s * s
        - 8 * s * root
        + 12 * s
        - 4 * mu * root
        + 6 * mu * s.log()
        - 3 * mu * (1 - mu).log()
        + 9 * mu
        + 4 * root
        - 6 * s.log()
        + 3 * (1 - mu).log()
        - 11
    ) / 6


def height_derivative_defect(s: arb, mu: arb) -> arb:
    root = (2 * s + mu - 1).sqrt()
    numerator = 4 * s + 2 * mu - 2 * root * s.log() + root * (1 - mu).log() - 2 * root - 2
    return -numerator / (2 * root)


def leading_weight(mode: int, c: arb, beta: arb, mu: arb) -> acb:
    s = 4 * arb(mode) / c
    ratio = ((1 - mu) * (2 * s + mu - 1)) ** (arb(1) / 4) / s.sqrt()
    phase = beta**3 * action_defect(s, mu)
    return (ratio * acb(0, phase).exp() - 1) / arb(mode).sqrt()


def weights_and_derivatives(lam: arb, y: arb, beta: arb) -> tuple[arb, arb, arb, arb, arb, arb]:
    epsilon = beta**-2
    u1 = -(13 * lam + 3 * y) / 60
    v1 = (8 * lam**2 - 4 * lam * y + 3 * y**2) / 60
    n_u2 = 448 * lam**5 + 280 * lam**2 * y**3 + 4565 * lam**2 - 105 * lam * y**4 + 30 * lam * y + 63 * y**5 - 405 * y**2
    n_v2 = 40 * lam**3 - 20 * lam**2 * y + lam * y**2 - 9 * y**3 + 27
    u = 1 + epsilon * u1 - epsilon**2 * n_u2 / 50400
    v = epsilon * v1 + epsilon**2 * n_v2 / 1680
    n_u2_y = 840 * lam**2 * y**2 - 420 * lam * y**3 + 30 * lam + 315 * y**4 - 810 * y
    n_v2_y = -20 * lam**2 + 2 * lam * y - 27 * y**2
    u_y = -epsilon / 20 - epsilon**2 * n_u2_y / 50400
    v_y = epsilon * (-4 * lam + 6 * y) / 60 + epsilon**2 * n_v2_y / 1680
    n_u2_yy = 1680 * lam**2 * y - 1260 * lam * y**2 + 1260 * y**3 - 810
    n_v2_yy = 2 * lam - 54 * y
    u_yy = -epsilon**2 * n_u2_yy / 50400
    v_yy = epsilon / 10 + epsilon**2 * n_v2_yy / 1680
    return u, v, u_y, v_y, u_yy, v_yy


def branch_value(lam: arb, y: arb, beta: arb, sigma: int) -> acb:
    x = lam + y
    ai, ai_prime, bi, bi_prime = acb(-x).airy()
    w = ai - acb(0, sigma) * bi
    w_x = -ai_prime + acb(0, sigma) * bi_prime
    u, v, _, _, _, _ = weights_and_derivatives(lam, y, beta)
    return (u * w + v * w_x) / 2


def grouped_components(mu: arb, c: arb, beta: arb, y_max: arb, panels: int = PANELS) -> tuple[acb, acb]:
    lam = beta**2 * mu
    minus_input = [acb(0) for _ in range(panels)]
    plus_input = [acb(0) for _ in range(panels)]
    for mode in range(L, U + 1):
        if mode == Q:
            continue
        coefficient = leading_weight(mode, c, beta, mu)
        offset = mode - L
        shifted = coefficient * acb(0, -arb.pi() * offset / panels).exp()
        if mode < Q:
            minus_input[offset] = shifted
        else:
            plus_input[offset] = shifted
    minus_dft = acb.dft(minus_input)
    plus_dft = acb.dft(plus_input)
    d_l = beta * (4 * arb(L) / c - 1)
    inner = acb(0)
    outer = acb(0)
    for index in range(panels):
        y = y_max * arb(2 * index + 1) / (2 * panels)
        common = acb(0, y**2 / (4 * beta) - d_l * y).exp()
        minus_branch = branch_value(lam, y, beta, -1)
        inner += common * minus_branch * minus_dft[index]
        outer += common * minus_branch.conjugate() * plus_dft[index]
    scale = y_max / panels
    return scale * inner, scale * outer


def component_error_bound(
    mode_start: int,
    mode_end: int,
    event_pairing: dict[str, Any],
    branch_envelope: dict[str, Any],
    c: arb,
    beta: arb,
    y_max: arb,
    lambda_max: arb,
    panels: int = PANELS,
) -> dict[str, arb]:
    modes = [mode for mode in range(mode_start, mode_end + 1) if mode != Q]
    pairing = event_pairing["leading_multiplier_certificate"]
    maximum_weight = arb(pairing["maximum_normalized_leading_defect_ball"]).upper()
    weight_l1 = maximum_weight * len(modes)
    detuning_max = max(abs(beta * (4 * arb(mode) / c - 1)).upper() for mode in modes)
    envelope = branch_envelope["certificate"]
    w_bound = arb(envelope["W_modulus_bound_ball"]).upper()
    wx_bound = arb(envelope["W_X_modulus_bound_ball"]).upper()
    atlas_x_min = arb(envelope["X_min_ball"]).upper()
    near_zero_x = arb(atlas_x_min / 2, atlas_x_min / 2)
    ai0, ai_prime0, bi0, bi_prime0 = acb(-near_zero_x).airy()
    w_bound = max(w_bound, abs(ai0 - acb(0, 1) * bi0).upper())
    wx_bound = max(wx_bound, abs(-ai_prime0 + acb(0, 1) * bi_prime0).upper())
    x_max = y_max + lambda_max
    lam = arb(lambda_max / 2, lambda_max / 2)
    y = arb(y_max / 2, y_max / 2)
    u, v, u_y, v_y, u_yy, v_yy = weights_and_derivatives(lam, y, beta)
    c0 = (abs(u).upper() * w_bound + abs(v).upper() * wx_bound) / 2
    c1 = ((abs(u_y).upper() + x_max.upper() * abs(v).upper()) * w_bound + (abs(u).upper() + abs(v_y).upper()) * wx_bound) / 2
    c2 = (
        (abs(u_yy).upper() + abs(v).upper() + x_max.upper() * abs(u).upper() + 2 * x_max.upper() * abs(v_y).upper()) * w_bound
        + (2 * abs(u_y).upper() + x_max.upper() * abs(v).upper() + abs(v_yy).upper()) * wx_bound
    ) / 2
    phase_one = y_max / (2 * beta) + detuning_max
    phase_two = 1 / (2 * beta)
    integrand_second = weight_l1 * (c2 + 2 * phase_one * c1 + (phase_two + phase_one**2) * c0)
    midpoint_error = y_max**3 * integrand_second / (24 * panels**2)

    mu_ball = arb(arb.pi() / (32 * beta**3), arb.pi() / (32 * beta**3))
    weight_lambda_derivative_l1 = arb(0)
    for mode in modes:
        s = 4 * arb(mode) / c
        ratio = ((1 - mu_ball) * (2 * s + mu_ball - 1)) ** (arb(1) / 4) / s.sqrt()
        log_ratio_mu = (-1 / (1 - mu_ball) + 1 / (2 * s + mu_ball - 1)) / 4
        weight_lambda_derivative_l1 += ratio.upper() * (
            abs(log_ratio_mu).upper() / beta**2 + beta * abs(height_derivative_defect(s, mu_ball)).upper()
        ) / arb(mode).sqrt()

    epsilon = beta**-2
    u_lambda = -13 * epsilon / 60 - epsilon**2 * (2240 * lam**4 + 560 * lam * y**3 + 9130 * lam - 105 * y**4 + 30 * y) / 50400
    v_lambda = epsilon * (16 * lam - 4 * y) / 60 + epsilon**2 * (120 * lam**2 - 40 * lam * y + y**2) / 1680
    c_lambda = ((abs(u_lambda).upper() + x_max.upper() * abs(v).upper()) * w_bound + (abs(u).upper() + abs(v_lambda).upper()) * wx_bound) / 2
    height_lipschitz = y_max * (weight_lambda_derivative_l1 * c0 + weight_l1 * c_lambda)
    return {
        "mode_count": len(modes),
        "weight_l1_bound": weight_l1,
        "midpoint_error": midpoint_error,
        "height_lipschitz": height_lipschitz,
        "height_transport_error": height_lipschitz * lambda_max / 8,
    }


def component_certificate(
    name: str,
    values: list[acb],
    errors: dict[str, arb],
    exact_profile_correction: arb,
) -> dict[str, Any]:
    adverse = errors["midpoint_error"].upper() + errors["height_transport_error"].upper() + exact_profile_correction.upper()
    real_upper = max(value.real.upper() for value in values) + adverse
    real_lower = min(value.real.lower() for value in values) - adverse
    modulus_upper = max(abs(value).upper() for value in values) + adverse
    scale = 4 * arb(2).sqrt()
    physical_real_upper = (scale * real_upper).upper()
    physical_real_lower = (scale * real_lower).lower()
    return {
        "name": name,
        "mode_count": errors["mode_count"],
        "midpoint_error_ball": errors["midpoint_error"].str(PRECISION, more=True),
        "height_transport_error_ball": errors["height_transport_error"].str(PRECISION, more=True),
        "exact_finite_t_profile_correction_bound_ball": exact_profile_correction.str(PRECISION, more=True),
        "corrected_canonical_real_lower_ball": real_lower.str(PRECISION, more=True),
        "corrected_canonical_real_upper_ball": real_upper.str(PRECISION, more=True),
        "corrected_canonical_modulus_upper_ball": modulus_upper.str(PRECISION, more=True),
        "physical_projected_real_lower_ball": physical_real_lower.str(PRECISION, more=True),
        "physical_projected_real_upper_ball": physical_real_upper.str(PRECISION, more=True),
        "uniformly_negative": physical_real_upper < 0,
        "uniformly_positive": physical_real_lower > 0,
    }


def render_note(artifact: dict[str, Any]) -> str:
    inner = artifact["certificate"]["inner"]
    outer = artifact["certificate"]["outer"]
    target = artifact["certificate"]["sufficient_paired_companion_real_upper_target"]
    return f"""# Inner/outer corrected-defect projection

Date: 2026-08-13

Status: rigorous branch-separated corrected-defect enclosure on the ordinary
top corridor; this is not yet a bound for the paired companion

Split the carrier-suppressed corrected statistic at event zero:

```text
W_corr=W_inner+W_outer,
W_inner=sum_(m=39696)^39894 a_m(G_ex,m-G_opp,m),
W_outer=sum_(m=39895)^40094 a_m(G_ex,m-G_opp,m),       (IO1)
a_39894=0.
```

The beta-minus-four parts were evaluated as separate 198- and 200-mode
Fourier polynomials on `{PANELS}` midpoint panels at five height nodes.  No
modewise integral or post-summation triangle was used.  Separate derivative
majorants give

```text
inner y error  <{inner['midpoint_error_ball']},
inner t error  <{inner['height_transport_error_ball']},
outer y error  <{outer['midpoint_error_ball']},
outer t error  <{outer['height_transport_error_ball']}.                 (IO2)
```

The exact Kummer-ODE profile correction is split with the two certified
whole-kernel norms, giving bounds

```text
inner exact finite-t correction <{inner['exact_finite_t_profile_correction_bound_ball']},
outer exact finite-t correction <{outer['exact_finite_t_profile_correction_bound_ball']}. (IO3)
```

After restoring the common carrier and paired projection, the rigorous
component intervals are

```text
{inner['physical_projected_real_lower_ball']} < R_inner
 < {inner['physical_projected_real_upper_ball']},
{outer['physical_projected_real_lower_ball']} < R_outer
 < {outer['physical_projected_real_upper_ball']}.      (IO4)
```

Their sum is the previously certified `R_corr`.  The outer component cancels
exactly in Section 11.406, so the operative fold increment is

```text
Delta Q_fold=R_inner+R_pair,
R_pair=2Re[e^(-i*pi/8)J_pair].                        (IO5)
```

In fact (IO4) certifies `R_inner<0` uniformly.  Consequently the explicit
sufficient target for a negative fold increment is

```text
R_pair<{target}.                                      (IO6)
```

No sign for the fold increment is claimed until (IO6) is proved on the same
corridor.

Pi provenance: the projection factor and corridor width come from
`beta^3=pi*C^2/8`, `hY=2pi`, inverse Airy Fourier normalization, and exact
odd-square half-domain reflection.  No fitted constant is introduced.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No paired-companion bound, signed fold increment, fixed-state residual,
all-corridor telescope, complete `Q_K-T` or `T_upper`, height-uniform theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    ctx.dps = PRECISION
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {
        name: json.loads(path.read_text(encoding="utf-8"))
        for name, path in DEPENDENCIES.items()
    }
    require(all(item.get("passed") is True for item in dependencies.values()), "dependency failure")
    require(
        dependencies["companion_decomposition"]["decision"]["outer_corrected_block_cancels_from_collapsed_fold_increment"] is True,
        "companion collapse drift",
    )

    c = arb(C)
    tstar = arb.pi() * c * c / 8
    beta = tstar ** (arb(1) / 3)
    mu_radius = arb.pi() / (16 * tstar)
    lambda_max = beta**2 * mu_radius
    y_max = arb.pi() * c / (2 * beta)
    rows = []
    inner_values: list[acb] = []
    outer_values: list[acb] = []
    old_rows = dependencies["selected_total"]["certificate"]["node_rows"]
    for index in range(HEIGHT_NODES):
        mu = mu_radius * index / (HEIGHT_NODES - 1)
        inner, outer = grouped_components(mu, c, beta, y_max)
        total = inner + outer
        stored = acb(
            arb(old_rows[index]["value_ball"]["real_ball"]),
            arb(old_rows[index]["value_ball"]["imag_ball"]),
        )
        require(abs(total - stored).upper() < arb("1e-55"), f"total component reconstruction drift at node {index}")
        inner_values.append(inner)
        outer_values.append(outer)
        rows.append({
            "height_fraction_index_over_4": index,
            "mu_ball": mu.str(PRECISION, more=True),
            "inner_beta4_value_ball": complex_record(inner),
            "outer_beta4_value_ball": complex_record(outer),
            "reconstructed_total_beta4_value_ball": complex_record(total),
        })

    inner_errors = component_error_bound(L, Q - 1, dependencies["event_pairing"], dependencies["branch_envelope"], c, beta, y_max, lambda_max)
    outer_errors = component_error_bound(Q + 1, U, dependencies["event_pairing"], dependencies["branch_envelope"], c, beta, y_max, lambda_max)
    profile_sup = arb(dependencies["exact_profile"]["interval_certificate"]["uniform_canonical_exact_minus_beta4_profile_bound"]).upper()
    kernel = dependencies["weighted_kernel"]["certificate"]
    inner_kernel_l1 = arb(kernel["minus_branch_kernel"]["L1_bound_ball"]).upper()
    outer_kernel_l1 = arb(kernel["plus_branch_kernel"]["L1_bound_ball"]).upper()
    inner_correction = (profile_sup * inner_kernel_l1).upper()
    outer_correction = (profile_sup * outer_kernel_l1).upper()
    inner_certificate = component_certificate("inner", inner_values, inner_errors, inner_correction)
    outer_certificate = component_certificate("outer", outer_values, outer_errors, outer_correction)
    require(inner_certificate["uniformly_negative"] is True, "inner corrected projection sign did not close")
    paired_target = (-arb(inner_certificate["physical_projected_real_upper_ball"]).upper()).lower()
    require(paired_target > arb("6.84e-5"), "paired companion target lost certified margin")

    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "rigorous_inner_outer_corrected_defect_projection_enclosure",
        "passed": True,
        "scope": {
            "height_interval": "t*-pi/16<=t<=t*",
            "inner_modes": [L, Q],
            "inner_nonzero_weight_count": Q - L,
            "event_zero_mode": Q,
            "outer_modes": [Q + 1, U],
            "outer_nonzero_weight_count": U - Q,
            "midpoint_panels": PANELS,
            "height_nodes": HEIGHT_NODES,
        },
        "node_rows": rows,
        "certificate": {
            "inner": inner_certificate,
            "outer": outer_certificate,
            "sum_identity": "W_corr=W_inner+W_outer",
            "operative_increment": "D_fold=C_inner+J_pair",
            "sufficient_paired_companion_real_upper_target": paired_target.str(PRECISION, more=True),
        },
        "decision": {
            "inner_outer_beta4_grouped_integrals_certified": True,
            "exact_finite_t_profile_correction_split_by_whole_kernel_norm": True,
            "inner_plus_outer_reconstructs_prior_corrected_statistic": True,
            "inner_corrected_projection_uniformly_negative": inner_certificate["uniformly_negative"],
            "inner_corrected_projection_uniformly_positive": inner_certificate["uniformly_positive"],
            "outer_corrected_projection_uniformly_negative": outer_certificate["uniformly_negative"],
            "outer_corrected_projection_uniformly_positive": outer_certificate["uniformly_positive"],
            "total_corrected_sign_transferred_without_split": False,
            "paired_companion_projection_bounded": False,
            "signed_fold_increment_proved": False,
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
            "Prove R_pair=2Re[e^(-i*pi/8)J_pair] is below the saved sufficient target "
            f"{paired_target.str(PRECISION, more=True)} on the ordinary top corridor, using a grouped paired-current "
            "representation without separating the endpoint/Fresnel, negative-mode, or event-zero channels."
        ),
        "proof_boundary": (
            "Rigorous inner/outer corrected weighted-defect projection enclosures on one ordinary top corridor only. "
            "No paired-companion bound, signed fold increment, fixed-state residual, all-corridor telescope, complete "
            "Q_K-T or T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print(
        "certified inner/outer corrected-defect projection: "
        f"inner_negative={inner_certificate['uniformly_negative']}, outer_negative={outer_certificate['uniformly_negative']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
