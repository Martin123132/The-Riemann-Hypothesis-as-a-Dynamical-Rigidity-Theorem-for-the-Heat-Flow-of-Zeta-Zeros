#!/usr/bin/env python3
"""Certify exact y primitives and explicit envelopes on the three y=64 charts."""

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


PARTITION_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_exact_normal_saddle_partition_gate.json"
FRESNEL_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fresnel_endpoint_normal_form_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_three_chart_primitive_envelope_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_three_chart_primitive_envelope_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 100
T = 10_000_000_000
C = 159577
B = 5122421
MODE_LO = 39853
MODE_HI = 39936
Y0 = 64
RADIUS = 9
FRESNEL_WIDTHS = 4


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


def symbolic_primitive() -> dict[str, str]:
    y, ystar, k, c_amp, h0, e_u, e_l = sp.symbols(
        "y ystar k c H0 E_U E_L", real=True, positive=True
    )
    h1 = ystar * h0 - sp.I * (e_u - e_l) / k
    affine = sp.expand(h0 + c_amp * h1)
    expected = (1 + c_amp * ystar) * h0 - sp.I * c_amp * (e_u - e_l) / k
    require(sp.simplify(affine - expected) == 0, "affine primitive identity failed")

    beta, mode, endpoint_c, x = sp.symbols("beta m C x", positive=True, real=True)
    sigma = endpoint_c / (2 * beta**2)
    c_scaled = sigma / endpoint_c
    detuning = 4 * beta * mode / endpoint_c - beta
    tau = 1 - 2 * x
    p = sp.simplify(detuning + beta * tau)
    exact_ystar = sp.simplify(beta * p / x)
    source_ystar = sp.simplify((2 * mode / x - endpoint_c) / sigma)
    require(sp.simplify(exact_ystar - source_ystar) == 0, "source saddle normalization failed")
    require(
        sp.simplify(endpoint_c * sigma * (1 + c_scaled * exact_ystar) - 2 * mode * sigma / x) == 0,
        "Fresnel-current coefficient failed",
    )

    pi = sp.pi
    endpoint_coefficient = sp.simplify(endpoint_c * sigma * c_scaled / (x / beta))
    relation_residual = sp.together(endpoint_coefficient - 2 / (pi * x))
    relation_residual = sp.simplify(relation_residual.subs(endpoint_c**2, 8 * beta**3 / pi))
    require(relation_residual == 0, "endpoint-current coefficient failed")

    y, k, p0 = sp.symbols("y k p0", positive=True, real=True)
    completed_phase = k * (y - p0 / k) ** 2 / 2 - p0**2 / (2 * k)
    require(
        sp.simplify(k * y**2 / 2 - p0 * y - completed_phase) == 0,
        "quadratic completion carrier failed",
    )

    return {
        "normal_phase": "phi(y)=k*y^2/2-p*y with k=x(z)/beta>0",
        "normal_saddle": "y_*=p/k=beta[d_m+beta*tanh(z/beta)]/x(z)",
        "zeroth_primitive": "H0(L,U)=int_L^U exp(i phi(y))dy",
        "first_moment": "H1(L,U)=y_* H0-i[E(U)-E(L)]/k",
        "affine_primitive": "K(L,U)=(1+c*y_*)H0-i(c/k)[E(U)-E(L)]",
        "source_saddle": "y_*=(2m/x-C)/sigma",
        "source_Fresnel_coefficient": "C*sigma*(1+c*y_*)=2m*sigma/x",
        "source_endpoint_coefficient": "C*sigma*c/k=2/(pi*x)",
        "source_identity": "C*sigma*K=2m*sigma*H0/x+2[E(U)-E(L)]/(i*pi*x)",
        "split_telescope": "K(0,y_B)=K(0,64)+K(64,y_B); the E(64) endpoint cancels exactly",
        "scale_relations": "pi*C^2=8beta^3, sigma=C/(2beta^2), c=sigma/C",
        "full_gaussian": "H_full=exp(-i*p^2/(2*k))*sqrt(2*pi/k)*exp(i*pi/4)",
    }


def certified_envelopes() -> dict[str, Any]:
    ctx.dps = PRECISION
    t = arb(T)
    endpoint_c = arb(C)
    pi = arb.pi()
    eta = pi * endpoint_c**2 / (8 * t)
    beta = (t * eta) ** (arb(1) / 3)
    sigma = 4 * beta / (pi * endpoint_c)
    c_amp = sigma / endpoint_c
    h = 4 * beta / endpoint_c
    detunings = [h * (arb(mode) - endpoint_c / 4) for mode in range(MODE_LO, MODE_HI + 1)]

    def tau(z: arb) -> arb:
        return (z / beta).tanh()

    def x(z: arb) -> arb:
        return (1 - tau(z)) / 2

    x_min = x(arb(RADIUS))
    x_max = x(arb(-RADIUS))
    k_min = x_min / beta
    k_max = x_max / beta
    y_b = arb(B - C) / sigma
    span = y_b - arb(Y0)
    delta = arb(FRESNEL_WIDTHS) * (x_max / beta).sqrt()

    g_global_min = detunings[0] + beta * tau(arb(-RADIUS)) - arb(Y0) * x_max / beta
    d0_max = -g_global_min
    d_b_min = delta + k_min * span
    d_b_max = d0_max + k_max * span
    amplitude_y0 = 1 + c_amp * arb(Y0)
    amplitude_yb = 1 + c_amp * y_b
    coefficient_margin = amplitude_y0 - c_amp * d0_max / k_min

    nonstationary_bound = (
        2 * amplitude_y0 / delta
        + amplitude_yb / d_b_min
        + 2 * c_amp / k_min * (d_b_max / delta).log()
    )

    partition = json.loads(PARTITION_GATE.read_text(encoding="utf-8"))["certified_partition"]
    max_saddle = arb(partition["maximum_normal_saddle_on_R9_ball"])
    upper_q_min = arb(partition["minimum_upper_endpoint_normalized_clearance_ball"])
    inverse_root_k = 1 / k_min.sqrt()
    max_morse_amplitude = 1 + c_amp * max_saddle
    morse_fresnel_remainder = max_morse_amplitude * inverse_root_k * (
        arb(2) / FRESNEL_WIDTHS + arb(2) / upper_q_min
    )
    morse_grouped_main = (
        max_morse_amplitude * (2 * pi / k_min).sqrt()
        + 2 * c_amp / k_min
    )

    q_transition_max = delta * (beta / x_min).sqrt()
    max_transition_saddle = arb(Y0) + beta * delta / x_min
    max_transition_amplitude = 1 + c_amp * max_transition_saddle
    transition_h0_bound = inverse_root_k * (
        2 * q_transition_max
        + 2 / q_transition_max
        + 2 / upper_q_min
    )
    transition_affine_bound = (
        max_transition_amplitude * transition_h0_bound
        + 2 * c_amp / k_min
    )

    require(coefficient_margin > arb("0.99"), "nonstationary coefficient sign margin failed")
    require(nonstationary_bound < arb(33), "nonstationary envelope exceeds 33")
    require(q_transition_max < arb("4.02"), "transition q range exceeds 4.02")
    require(transition_affine_bound < arb(600), "transition envelope exceeds 600")
    require(morse_fresnel_remainder < arb(34), "Morse Fresnel remainder exceeds 34")
    require(morse_grouped_main < arb(200), "grouped Morse main exceeds 200")
    require(abs(amplitude_yb - arb(B) / endpoint_c) < arb("1e-90"), "upper amplitude identity failed")

    return {
        "height": T,
        "mode_range": [MODE_LO, MODE_HI],
        "z_range": [-RADIUS, RADIUS],
        "outer_y_range": [Y0, "y_B"],
        "beta_ball": beta.str(PRECISION, more=True),
        "sigma_ball": sigma.str(PRECISION, more=True),
        "affine_coefficient_c_ball": c_amp.str(PRECISION, more=True),
        "x_min_ball": x_min.str(PRECISION, more=True),
        "x_max_ball": x_max.str(PRECISION, more=True),
        "k_min_ball": k_min.str(PRECISION, more=True),
        "k_max_ball": k_max.str(PRECISION, more=True),
        "y_B_ball": y_b.str(PRECISION, more=True),
        "uniform_g_buffer_delta_ball": delta.str(PRECISION, more=True),
        "nonstationary_boundary_derivative_max_ball": d0_max.str(PRECISION, more=True),
        "nonstationary_upper_derivative_min_ball": d_b_min.str(PRECISION, more=True),
        "nonstationary_upper_derivative_max_ball": d_b_max.str(PRECISION, more=True),
        "nonstationary_coefficient_margin_ball": coefficient_margin.str(PRECISION, more=True),
        "per_mode_nonstationary_outer_bound_ball": nonstationary_bound.str(PRECISION, more=True),
        "transition_normalized_q_max_ball": q_transition_max.str(PRECISION, more=True),
        "transition_H0_bound_ball": transition_h0_bound.str(PRECISION, more=True),
        "per_mode_transition_affine_bound_ball": transition_affine_bound.str(PRECISION, more=True),
        "maximum_transition_saddle_ball": max_transition_saddle.str(PRECISION, more=True),
        "per_mode_grouped_Morse_main_bound_ball": morse_grouped_main.str(PRECISION, more=True),
        "per_mode_Morse_Fresnel_remainder_bound_ball": morse_fresnel_remainder.str(PRECISION, more=True),
        "maximum_Morse_amplitude_ball": max_morse_amplitude.str(PRECISION, more=True),
        "upper_normalized_q_min_ball": upper_q_min.str(PRECISION, more=True),
        "upper_affine_amplitude_ball": amplitude_yb.str(PRECISION, more=True),
        "bound_derivation": {
            "nonstationary": "One integration by parts with phi'(64)>=delta and phi''=k>0; the affine derivative and curvature integrals are evaluated in logarithmic closed form before uniformization.",
            "transition": "Complete the square exactly; split the normalized Fresnel interval at +/-q_max and use |int_Q^infinity exp(iq^2/2)dq|<=2/Q.",
            "Morse": "Keep the exact endpoint current in the grouped main and replace only H0 by its full Gaussian; the two omitted normalized Fresnel tails are bounded by 2/|q_0|+2/q_B.",
        },
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certified_envelopes"]
    return f"""# Exact y=64 three-chart primitives and envelopes

Date: 2026-08-11

Status: exact affine primitive and explicit saved-height chart envelopes
validated; not a proof of the remaining z integral or grouped mode sum

For fixed `(z,m)`, put

```text
k=x(z)/beta,        p=d_m+beta*tanh(z/beta),
y_*=p/k,            c=sigma/C,
E(y)=exp(i[k y^2/2-py]),
H0(L,U)=integral_L^U E(y)dy.                            (Y3E1)
```

Direct differentiation gives the exact affine primitive

```text
K(L,U)=integral_L^U(1+cy)E(y)dy
      =(1+c y_*)H0(L,U)-i(c/k)[E(U)-E(L)].              (Y3E2)
```

Using `pi*C^2=8beta^3` and `sigma=C/(2beta^2)`, (Y3E2) becomes

```text
C sigma K
 =2m sigma H0/x+2[E(U)-E(L)]/(i pi x).                 (Y3E3)
```

Thus (Y3E3) is exactly the source's paired Fresnel and endpoint current.  In
particular

```text
K(0,y_B)=K(0,64)+K(64,y_B),                            (Y3E4)
```

and the artificial `E(64)` endpoint cancels algebraically.  The split creates
no new source term and no double counting.

On the nonstationary chart `g_m(z)<=-delta`, one integration by parts uses
`phi'(64)>=delta` and `phi''=k>0`.  Evaluating the affine derivative and
curvature integrals before taking uniform bounds gives

```text
|K(64,y_B)|
 <= {c['per_mode_nonstationary_outer_bound_ball']} <33. (Y3E5)
```

On the endpoint/Fresnel chart, the normalized lower endpoint remains inside
`|q_0|<={c['transition_normalized_q_max_ball']}`.  Completing the square and
using `|int_Q^infinity exp(iq^2/2)dq|<=2/Q` outside that compact interval
gives

```text
|K(64,y_B)|
 <= {c['per_mode_transition_affine_bound_ball']} <600. (Y3E6)
```

This is an envelope only; the incomplete Fresnel function remains exact in
the working representation.

On the ordinary Morse chart, define the grouped main by replacing only `H0`
in (Y3E2) with the full Gaussian

```text
H_full=exp(-i p^2/(2k))sqrt(2pi/k)exp(i*pi/4),          (Y3E7)
```

while retaining the exact endpoint current in (Y3E2).  Then

```text
|K-K_Morse,grouped|
 <= {c['per_mode_Morse_Fresnel_remainder_bound_ball']} <34. (Y3E8)
```

and the grouped main itself is below

```text
{c['per_mode_grouped_Morse_main_bound_ball']} <200.    (Y3E9)
```

The constants in (Y3E5)--(Y3E9) are per-mode absolute envelopes.  They are
not summed over 84 modes and are not claimed to be small enough for the final
source error.  Their purpose is to close the local chart analysis without
losing the exact endpoint/Fresnel pairing; the next estimate must recover
cancellation in the remaining `z` integral and mode sum.

Pi provenance: every pi in (Y3E3) and (Y3E7) comes from the Kummer quadratic
phase, Fourier--Poisson character, and standard Fresnel Gaussian.  No fitted
geometric normalization is introduced.

Proof boundary: exact source-normalized primitives and explicit per-mode
saved-height envelopes only.  No completed z integral, grouped 84-mode sum,
join to all lower-interior modes, `T_upper` theorem, `Lambda<=0`, RH, or
prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    require(PARTITION_GATE.is_file() and FRESNEL_GATE.is_file() and CHECKER.is_file(), "missing dependency or checker")
    symbolic = symbolic_primitive()
    certified = certified_envelopes()
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_three_chart_primitive_envelope_gate",
        "status": "exact_y64_affine_primitive_source_pairing_and_three_chart_envelopes_complete",
        "passed": True,
        "symbolic_primitive": symbolic,
        "certified_envelopes": certified,
        "decision": {
            "exact_affine_normal_primitive_proved": True,
            "source_endpoint_Fresnel_pairing_recovered": True,
            "artificial_y64_endpoint_cancels_exactly": True,
            "nonstationary_per_mode_envelope_proved": True,
            "transition_per_mode_envelope_proved": True,
            "grouped_Morse_remainder_per_mode_envelope_proved": True,
            "z_integral_completed": False,
            "grouped_84_mode_bound_proved": False,
            "lower_interior_join_proved": False,
            "rh_implication": False,
        },
        "next_obligation": "Insert the exact three-chart primitives into the z integral, preserve the 84-mode Dirichlet cancellation through the mode-dependent cutoffs, and match the grouped Morse carrier to modes 622..39852 before taking absolute values.",
        "proof_boundary": "Exact source-normalized primitives and explicit per-mode saved-height envelopes only. No completed z integral, grouped 84-mode sum, join to all lower-interior modes, T_upper theorem, Lambda<=0, RH, or prize-level conclusion is proved.",
        "dependencies": {
            "partition_gate": {"path": relative(PARTITION_GATE), "sha256": file_hash(PARTITION_GATE)},
            "fresnel_gate": {"path": relative(FRESNEL_GATE), "sha256": file_hash(FRESNEL_GATE)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {"workers": 1, "process_priority": priority, "elapsed_seconds": round(time.perf_counter() - started, 3)},
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("built y=64 three-chart primitive envelopes: nonstat<33, transition<600, Morse remainder<34")


if __name__ == "__main__":
    main()
