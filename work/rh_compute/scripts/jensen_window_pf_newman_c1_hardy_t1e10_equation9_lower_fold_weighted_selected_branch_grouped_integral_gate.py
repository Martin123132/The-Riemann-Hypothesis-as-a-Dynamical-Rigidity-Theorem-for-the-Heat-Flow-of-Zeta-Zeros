#!/usr/bin/env python3
"""Certify the leading-defect selected beta^-4 grouped integral directly."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_selected_branch_grouped_integral_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "event_pairing": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_event_ordered_398_pairing_multiplier_gate.json",
    "detuning_ode": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_beta4_endpoint_coherent_detuning_ode_gate.json",
    "branch_envelope": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_scaled_airy_branch_amplitude_envelope_gate.json",
    "weighted_opposite": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_opposite_branch_dirichlet_gate.json",
}

C = 159_577
Q = 39_894
L = 39_696
U = 40_094
Y_PANELS = 65_536
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
    return {"real_ball": value.real.str(PRECISION, more=True), "imag_ball": value.imag.str(PRECISION, more=True)}


def action_defect(s: arb, mu: arb) -> arb:
    root = (2 * s + mu - 1).sqrt()
    return (
        3 * s * s - 8 * s * root + 12 * s - 4 * mu * root
        + 6 * mu * s.log() - 3 * mu * (1 - mu).log() + 9 * mu
        + 4 * root - 6 * s.log() + 3 * (1 - mu).log() - 11
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


def grouped_point(mu: arb, c: arb, beta: arb, tstar: arb, y_max: arb) -> acb:
    lam = beta**2 * mu
    coefficients = [acb(0) for _ in range(U - L + 1)]
    for mode in range(L, U + 1):
        if mode != Q:
            coefficients[mode - L] = leading_weight(mode, c, beta, mu)
    minus_input = [acb(0) for _ in range(Y_PANELS)]
    plus_input = [acb(0) for _ in range(Y_PANELS)]
    for offset, coefficient in enumerate(coefficients):
        shifted = coefficient * acb(0, -arb.pi() * offset / Y_PANELS).exp()
        if L + offset < Q:
            minus_input[offset] = shifted
        elif L + offset > Q:
            plus_input[offset] = shifted
    minus_dft = acb.dft(minus_input)
    plus_dft = acb.dft(plus_input)
    h = 4 * beta / c
    d_l = beta * (4 * arb(L) / c - 1)
    total = acb(0)
    for index in range(Y_PANELS):
        y = y_max * arb(2 * index + 1) / (2 * Y_PANELS)
        minus_kernel = minus_dft[index]
        plus_kernel = plus_dft[index]
        common = acb(0, y**2 / (4 * beta) - d_l * y).exp()
        minus_branch = branch_value(lam, y, beta, -1)
        total += common * (
            minus_branch * minus_kernel
            + minus_branch.conjugate() * plus_kernel
        )
    return y_max * total / Y_PANELS


def derivative_error_bound(pairing: dict[str, Any], branch_envelope: dict[str, Any], c: arb, beta: arb, y_max: arb, lambda_max: arb) -> dict[str, arb]:
    certificate = pairing["leading_multiplier_certificate"]
    maximum_weight = arb(certificate["maximum_normalized_leading_defect_ball"]).upper()
    weight_l1 = maximum_weight * 398
    detuning_max = max(abs(beta * (4 * arb(L) / c - 1)).upper(), abs(beta * (4 * arb(U) / c - 1)).upper())
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
    midpoint_error = y_max**3 * integrand_second / (24 * Y_PANELS**2)

    mu = arb(arb.pi() / (32 * beta**3), arb.pi() / (32 * beta**3))
    weight_lambda_derivative_l1 = arb(0)
    for mode in range(L, U + 1):
        if mode == Q:
            continue
        s = 4 * arb(mode) / c
        ratio = ((1 - mu) * (2 * s + mu - 1)) ** (arb(1) / 4) / s.sqrt()
        log_ratio_mu = (-1 / (1 - mu) + 1 / (2 * s + mu - 1)) / 4
        weight_lambda_derivative_l1 += ratio.upper() * (
            abs(log_ratio_mu).upper() / beta**2 + beta * abs(height_derivative_defect(s, mu)).upper()
        ) / arb(mode).sqrt()

    epsilon = beta**-2
    u_lambda = -13 * epsilon / 60 - epsilon**2 * (
        2240 * lam**4 + 560 * lam * y**3 + 9130 * lam - 105 * y**4 + 30 * y
    ) / 50400
    v_lambda = epsilon * (16 * lam - 4 * y) / 60 + epsilon**2 * (
        120 * lam**2 - 40 * lam * y + y**2
    ) / 1680
    c_lambda = (
        (abs(u_lambda).upper() + x_max.upper() * abs(v).upper()) * w_bound
        + (abs(u).upper() + abs(v_lambda).upper()) * wx_bound
    ) / 2
    height_lipschitz = y_max * (weight_lambda_derivative_l1 * c0 + weight_l1 * c_lambda)
    return {
        "weight_l1_bound": weight_l1,
        "integrand_second_derivative_bound": integrand_second,
        "composite_midpoint_error_bound": midpoint_error,
        "lambda_derivative_bound": height_lipschitz,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Grouped selected-branch leading-defect integral

Date: 2026-08-13

Status: rigorous direct grouped midpoint/height enclosure; valid but too
coarse for the local normalized headroom

The selected leading-defect correction is kept as one finite integral,

```text
W_sel(mu)=Integral_0^Y e^(iy^2/(4beta))
 sum_(m!=q) a_m(mu)e^(-id_my)C_sel,m(lambda,y)dy.    (SG1)
```

At each height the two mode halves are evaluated as whole 198- and 200-term
Fourier polynomials.  No modewise integral or absolute-value sum is taken.
Composite midpoint quadrature on `{Y_PANELS}` panels is enclosed using a
uniform second-derivative bound; five height nodes are transported over their
adjacent subintervals by an analytic `lambda`-derivative majorant.

The resulting ordinary-top-corridor certificate is

```text
max node modulus <{c['maximum_node_absolute_ball']},
y-quadrature error <{c['y_midpoint_error_ball']},
height transport error <{c['height_transport_error_ball']},
sup_mu |W_sel(mu)|<{c['uniform_selected_weighted_bound_ball']}.       (SG2)
```

The direct grouped values remain near `3.2e-5`, but the present analytic
derivative envelope is deliberately conservative.  The bound is therefore
a certified interface, not yet a closure of the previous local headroom.
Its next sharpening should exploit the exact Airy equation inside the
grouped kernel rather than a global `sum|a_m|` derivative majorant.

Pi provenance: the height radius is `pi/(16t*)`, the fold scale obeys
`beta^3=pi C^2/8`, and the mode phases use `hY=2pi`.  No unrelated circle
normalization enters.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No exact finite-integral amplitude remainder, local-headroom closure,
completed source/initial-data splice, complete `Q_K-T` or `T_upper`,
all-corridor theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion
is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    for name, dependency in dependencies.items():
        require(dependency.get("passed") is True, f"dependency failed: {name}")
    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    c = arb(C)
    tstar = pi * c**2 / 8
    beta = tstar ** (arb(1) / 3)
    y_max = pi * c / (2 * beta)
    mu_max = pi / (16 * tstar)
    lambda_max = beta**2 * mu_max
    node_rows = []
    maximum_node = arb(0)
    for index in range(5):
        mu = mu_max * index / 4
        value = grouped_point(mu, c, beta, tstar, y_max)
        maximum_node = max(maximum_node, abs(value).upper())
        node_rows.append({"height_fraction_index_over_4": index, "mu_ball": mu.str(PRECISION, more=True), "value_ball": complex_record(value), "absolute_ball": abs(value).str(PRECISION, more=True)})
    errors = derivative_error_bound(dependencies["event_pairing"], dependencies["branch_envelope"], c, beta, y_max, lambda_max)
    height_transport = errors["lambda_derivative_bound"] * lambda_max / 8
    uniform = maximum_node + errors["composite_midpoint_error_bound"] + height_transport
    require(maximum_node < arb("3.3e-5"), "grouped node value exceeded 3.3e-5")
    require(errors["composite_midpoint_error_bound"] < arb("5e-6"), "midpoint error exceeded 5e-6")
    require(uniform < arb("0.001"), "uniform selected weighted bound exceeded fail-safe cap")
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "grouped_selected_leading_defect_integral_certified_with_conservative_height_transport",
        "passed": True,
        "certificate": {
            "height_interval": "t*-pi/16<=t<=t*",
            "y_panel_count": Y_PANELS,
            "height_node_count": 5,
            "node_rows": node_rows,
            "maximum_node_absolute_ball": maximum_node.str(PRECISION, more=True),
            "y_midpoint_error_ball": errors["composite_midpoint_error_bound"].str(PRECISION, more=True),
            "lambda_derivative_bound_ball": errors["lambda_derivative_bound"].str(PRECISION, more=True),
            "height_transport_error_ball": height_transport.str(PRECISION, more=True),
            "uniform_selected_weighted_bound_ball": uniform.str(PRECISION, more=True),
            "integrand_second_derivative_bound_ball": errors["integrand_second_derivative_bound"].str(PRECISION, more=True),
        },
        "decision": {
            "grouped_selected_leading_defect_integral_bound_proved": True,
            "local_normalized_headroom_closed": False,
            "exact_finite_integral_amplitude_remainder_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {name: {"path": relative(path), "sha256": file_hash(path)} for name, path in DEPENDENCIES.items()},
        "sources": {"builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)}, "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)}},
        "runtime": {"workers": 1, "process_priority": priority, "elapsed_seconds": round(time.perf_counter() - started, 3)},
        "next_action": "Replace the global height-derivative triangle majorant by a grouped Airy-ODE/weighted-Fourier derivative packet, or prove a direct interval enclosure on many height panels. Then address the exact finite-integral-minus-beta^-4 remainder.",
        "proof_boundary": "Direct grouped leading-defect selected-branch integral bound with conservative height transport only. No local-headroom closure, exact finite-integral amplitude remainder, completed source/initial-data splice, complete Q_K-T or T_upper, all-corridor theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("certified grouped selected-branch leading-defect integral", flush=True)


if __name__ == "__main__":
    main()
