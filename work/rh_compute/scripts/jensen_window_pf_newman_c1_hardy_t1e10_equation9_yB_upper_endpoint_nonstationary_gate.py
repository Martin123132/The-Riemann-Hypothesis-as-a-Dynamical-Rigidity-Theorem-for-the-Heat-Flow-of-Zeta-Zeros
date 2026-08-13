#!/usr/bin/env python3
"""Certify the grouped B-endpoint current by nonstationary z integration."""

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


REASSEMBLY_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_endpoint_roster_reassembly_gate.json"
NORMAL_FORM_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_grouped_characteristic_airy_boundary_normal_form_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_yB_upper_endpoint_nonstationary_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_yB_upper_endpoint_nonstationary_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 100
T = 10_000_000_000
C = 159577
B = 5122421
MODE_COUNT = 84
RADIUS = 9


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


def symbolic_reduction() -> dict[str, str]:
    beta, endpoint_c, endpoint_b, t, x = sp.symbols("beta C B t x", positive=True, real=True)
    pi = sp.pi
    tau = 1 - 2 * x
    eta = pi * endpoint_c**2 / (8 * t)
    sigma = 4 * beta / (pi * endpoint_c)
    y_b = (endpoint_b - endpoint_c) / sigma
    transformed = -t * eta * tau - beta * y_b * tau + x * y_b**2 / (2 * beta)
    difference = sp.factor(transformed - pi * endpoint_b**2 * x / 4)
    reduced = sp.factor(difference.subs(beta**3, pi * endpoint_c**2 / 8))
    require(reduced == -pi * endpoint_c * (2 * endpoint_b - endpoint_c) / 8, "B-phase reduction failed")
    return {
        "transformed_common_phase": "t*z/beta-t*eta*tanh(z/beta)-beta*y_B*tanh(z/beta)+x(z)*y_B^2/(2beta)",
        "source_phase_reduction": "P_B(z)=t*z/beta+pi*B^2*x(z)/4 plus a z-independent constant",
        "phase_derivative": "P_B'(z)=[t-(pi*B^2/8)sech(z/beta)^2]/beta",
        "phase_second_derivative": "P_B''(z)=(pi*B^2/(4beta^2))sech(z/beta)^2*tanh(z/beta)",
        "grouped_amplitude": "A_B(z)=(sqrt(2)/pi)*84*(sigma/C)*(beta/x(z))*cosh(z/beta)^(-3/2)",
        "integration_by_parts": "|int A_B exp(iP_B)|<=2Amax/mu+2R[A1max/mu+Amax*P2max/mu^2]",
    }


def certified_bound() -> dict[str, Any]:
    ctx.dps = PRECISION
    pi = arb.pi()
    t = arb(T)
    endpoint_c = arb(C)
    endpoint_b = arb(B)
    radius = arb(RADIUS)
    beta = (pi * endpoint_c**2 / 8) ** (arb(1) / 3)
    sigma = 4 * beta / (pi * endpoint_c)
    u = radius / beta
    tau_r = u.tanh()
    sech2_r = u.cosh() ** -2
    x_min = (1 - tau_r) / 2

    prefactor = arb(2).sqrt() / pi
    amplitude_max = (
        prefactor
        * MODE_COUNT
        * (sigma / endpoint_c)
        * (beta / x_min)
        * u.cosh() ** (-arb(3) / 2)
    )
    log_derivative_max = (1 + tau_r / 2) / beta
    amplitude_derivative_max = amplitude_max * log_derivative_max
    derivative_min = (pi * endpoint_b**2 * sech2_r / 8 - t) / beta
    second_derivative_max = pi * endpoint_b**2 * sech2_r * tau_r / (4 * beta**2)
    ibp_bound = (
        2 * amplitude_max / derivative_min
        + 2 * radius * (
            amplitude_derivative_max / derivative_min
            + amplitude_max * second_derivative_max / derivative_min**2
        )
    )

    require(derivative_min > arb("4.77e9"), "B-endpoint phase derivative is too small")
    require(amplitude_max < arb("0.0177"), "B-endpoint amplitude exceeds 0.0177")
    require(ibp_bound < arb("7.5e-12"), "B-endpoint IBP bound exceeds 7.5e-12")

    return {
        "height": T,
        "source_endpoints": [C, B],
        "mode_count": MODE_COUNT,
        "z_range": [-RADIUS, RADIUS],
        "beta_ball": beta.str(PRECISION, more=True),
        "sigma_ball": sigma.str(PRECISION, more=True),
        "tanh_R_over_beta_ball": tau_r.str(PRECISION, more=True),
        "sech2_R_over_beta_ball": sech2_r.str(PRECISION, more=True),
        "minimum_x_ball": x_min.str(PRECISION, more=True),
        "grouped_amplitude_max_ball": amplitude_max.str(PRECISION, more=True),
        "grouped_amplitude_derivative_max_ball": amplitude_derivative_max.str(PRECISION, more=True),
        "phase_derivative_absolute_min_ball": derivative_min.str(PRECISION, more=True),
        "phase_second_derivative_max_ball": second_derivative_max.str(PRECISION, more=True),
        "grouped_B_endpoint_integral_bound_ball": ibp_bound.str(PRECISION, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certified_bound"]
    return f"""# Grouped B-endpoint nonstationary enclosure

Date: 2026-08-11

Status: exact upper-endpoint phase reduction and explicit R=9
nonstationary bound validated; not a proof of the lower characteristic fold

After the 84 endpoint characters are reassembled as in the y=64 endpoint
gate, the transformed common phase at `y=y_B` reduces, up to a
`z`-independent constant, to the original source endpoint phase

```text
P_B(z)=t*z/beta+pi*B^2*x(z)/4.                         (BE1)
```

The reduction uses only

```text
x=(1-tanh(z/beta))/2,
sigma=4beta/(pi C),
beta^3=pi C^2/8.                                       (BE2)
```

Thus

```text
P_B'(z)=[t-(pi B^2/8)sech(z/beta)^2]/beta,             (BE3)
P_B''(z)=(pi B^2/(4beta^2))sech(z/beta)^2
          tanh(z/beta).                                (BE4)
```

There is no `B`-endpoint stationary point on `|z|<=9`.  Interval arithmetic
gives

```text
|P_B'| >= {c['phase_derivative_absolute_min_ball']} >4.77e9. (BE5)
```

The complete grouped endpoint amplitude, including `D_84(y_B)=-84`, obeys

```text
|A_B| <= {c['grouped_amplitude_max_ball']} <0.0177,    (BE6)
|A_B'| <= {c['grouped_amplitude_derivative_max_ball']}. (BE7)
```

One integration by parts on `[-9,9]` therefore gives

```text
|integral_-9^9 A_B(z)exp(iP_B(z))dz|
 <=2Amax/mu+18[A1max/mu+Amax*P2max/mu^2]
 <= {c['grouped_B_endpoint_integral_bound_ball']}
 <7.5e-12.                                             (BE8)
```

This closes the physical `B` endpoint current for the 84 transition modes
on the saved R=9 chart.  It does not estimate the lower `C` endpoint/Fresnel
fold, whose phase is characteristic and must remain grouped.

Pi provenance: the `pi` in (BE1)--(BE4) is inherited from the source Kummer
quadratic and Fourier-Poisson character through the exact scale identities;
no fitted geometric normalization is used.

Proof boundary: a saved-height grouped `B`-endpoint enclosure on `|z|<=9`
only.  No lower characteristic-fold bound, lower-interior join, complete
`T_upper`, height-uniform source theorem, `Lambda<=0`, RH, or prize-level
conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for dependency in (REASSEMBLY_GATE, NORMAL_FORM_GATE):
        require(dependency.is_file(), f"missing dependency: {dependency}")
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_yB_upper_endpoint_nonstationary_gate",
        "status": "grouped_yB_upper_endpoint_R9_nonstationary_bound_complete",
        "passed": True,
        "symbolic": symbolic_reduction(),
        "certified_bound": certified_bound(),
        "decision": {
            "exact_transformed_to_source_B_phase_reduction_proved": True,
            "B_endpoint_nonstationary_on_R9_proved": True,
            "grouped_B_endpoint_below_7_5e_minus_12": True,
            "lower_characteristic_fold_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            "reassembly_gate": {"path": relative(REASSEMBLY_GATE), "sha256": file_hash(REASSEMBLY_GATE)},
            "normal_form_gate": {"path": relative(NORMAL_FORM_GATE), "sha256": file_hash(NORMAL_FORM_GATE)},
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
        "next_obligation": "Apply a characteristic-fold-uniform estimate to the lower C endpoint/Fresnel current, then join its ordinary Morse side to modes 622..39852 and compare the combined carrier with T_upper.",
        "proof_boundary": "Saved-height grouped B-endpoint nonstationary enclosure only. No lower fold bound, lower-interior join, complete T_upper theorem, Lambda<=0, RH, or prize-level conclusion is proved.",
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("built grouped B-endpoint enclosure: |P'|>4.77e9, integral<7.5e-12")


if __name__ == "__main__":
    main()
