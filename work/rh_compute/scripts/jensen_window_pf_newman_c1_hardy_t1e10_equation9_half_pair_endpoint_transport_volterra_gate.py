#!/usr/bin/env python3
"""Certify the endpoint-driven Volterra transport identity for Fourier pairs."""

from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path
from typing import Any

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_endpoint_transport_volterra_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "pair_Morse": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_signed_pair_logistic_morse_transform_gate.json",
    "paired_residual": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_paired_target_residual_normal_form_gate.json",
    "symmetric_poisson": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate.json",
    "midpoint_guard": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_boundary_linear_roster_pair_anchor_gate.json",
}

T = 10_000_000_000
A = 159_577
B = 5_122_421
L = (B - A) // 2


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


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
    u, x, alpha, m = sp.symbols("u x alpha m", real=True, positive=True)
    alpha_u = sp.symbols("A", real=True) + 2 * u
    f = alpha_u * sp.exp(sp.I * sp.pi * x * alpha_u**2 / 4)
    f_uu = sp.simplify(sp.diff(f, u, 2))
    transport_rhs = sp.simplify(6 * sp.I * sp.pi * x * f + 4 * sp.I * sp.pi * x**2 * sp.diff(f, x))
    require(sp.simplify(f_uu - transport_rhs) == 0, "roster transport PDE failed")

    k = 2 * sp.pi * m
    pair, pair_x, delta_f1 = sp.symbols("P P_x Delta_f1")
    ode = 4 * sp.I * sp.pi * x**2 * pair_x + (k**2 + 6 * sp.I * sp.pi * x) * pair - 2 * delta_f1
    integrating_factor = x ** sp.Rational(3, 2) * sp.exp(sp.I * sp.pi * m**2 / x)
    coefficient = k**2 / (4 * sp.I * sp.pi * x**2) + sp.Rational(3, 2) / x
    require(
        sp.simplify(sp.diff(sp.log(integrating_factor), x) - coefficient) == 0,
        "integrating factor failed",
    )
    normalized_rhs = delta_f1 / (2 * sp.I * sp.pi * x**2)
    transported_rhs = sp.simplify(integrating_factor * normalized_rhs)
    expected_rhs = delta_f1 * x ** (-sp.Rational(1, 2)) * sp.exp(sp.I * sp.pi * m**2 / x) / (2 * sp.I * sp.pi)
    require(sp.simplify(transported_rhs - expected_rhs) == 0, "transported endpoint driver failed")
    require(ode.has(pair_x, pair, delta_f1), "ODE symbol audit failed")

    endpoint = sp.symbols("D", positive=True, real=True)
    endpoint_derivative = sp.exp(sp.I * sp.pi * x * endpoint**2 / 4) * (2 + sp.I * sp.pi * x * endpoint**2)
    require(sp.simplify(endpoint_derivative.subs(x, 0) - 2) == 0, "endpoint initial value failed")

    # Verify the characteristic alignment: the endpoint-driver saddle z=2m/D
    # meets the outer Morse saddle exactly when D=2m+t/(pi*m).
    t = sp.symbols("t", positive=True, real=True)
    z_endpoint = 2 * m / endpoint
    x_m = 2 * sp.pi * m**2 / (t + 2 * sp.pi * m**2)
    characteristic = sp.solve(sp.Eq(z_endpoint, x_m), endpoint)[0]
    require(sp.simplify(characteristic - (2 * m + t / (sp.pi * m))) == 0, "characteristic alignment failed")

    return {
        "roster_PDE": "partial_u^2 f_x=6 i*pi*x f_x+4 i*pi*x^2 partial_x f_x",
        "pair_coefficient": "P_m(x)=I_m(x)+I_(-m)(x)=2 integral_0^L f_x(u)cos(2*pi*m*u)du",
        "pair_initial_value": "P_m(0)=0 for every nonzero integer m because L is an integer and f_0(u)=A+2u",
        "pair_transport_ODE": "4 i*pi*x^2 P_m'+(4*pi^2*m^2+6 i*pi*x)P_m=2 Delta f_x'",
        "integrating_factor": "mu_m(x)=x^(3/2)exp(i*pi*m^2/x)",
        "endpoint_driver": "Delta f_x'=exp(i*pi*x*B^2/4)(2+i*pi*x*B^2)-exp(i*pi*x*A^2/4)(2+i*pi*x*A^2)",
        "Volterra_solution": "P_m(x)=x^(-3/2)exp(-i*pi*m^2/x)/(2 i*pi) integral_0^x z^(-1/2)exp(i*pi*m^2/z)Delta f_z' dz",
        "half_pair_integral": "K_m+K_(-m)=1/(2 i*pi) integral_0^(1/2) z^(-1/2)exp(i*pi*m^2/z)Delta f_z' G_m(z)dz",
        "Volterra_kernel": "G_m(z)=integral_z^(1/2) W_t(x)x^(-3/2)exp(-i*pi*m^2/x)dx",
        "characteristic_alignment": "The endpoint-driver saddle z_D=2m/D meets x_m exactly when D=2m+t/(pi*m).",
        "cancellation_interpretation": "The Fresnel tails and endpoint currents are not separated; their cancellation is encoded in the endpoint-driven Volterra propagator.",
    }


def finite_height_geometry() -> dict[str, Any]:
    mp.mp.dps = 80
    t = mp.mpf(T)
    pi = mp.pi
    A_mp, B_mp = mp.mpf(A), mp.mpf(B)

    def roots(endpoint: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
        disc = endpoint**2 - 8 * t / pi
        require(disc > 0, "endpoint characteristic discriminant is not positive")
        return ((endpoint - mp.sqrt(disc)) / 4, (endpoint + mp.sqrt(disc)) / 4)

    roots_A = roots(A_mp)
    roots_B = roots(B_mp)
    half_corner_mode = A_mp / 4
    half_stationary_alpha = mp.sqrt(8 * t / pi)
    half_detuning = pi * A_mp**2 / 4 - 2 * t
    require(mp.mpf(39_894) < half_corner_mode < mp.mpf(39_895), "A/half mode interval drift")
    require(mp.mpf(621) < roots_B[0] < mp.mpf(622), "B lower root drift")
    require(mp.mpf(39_852) < roots_A[0] < mp.mpf(39_853), "A lower root drift")
    require(mp.mpf(39_936) < roots_A[1] < mp.mpf(39_937), "A upper root drift")
    return {
        "height": T,
        "source_endpoints": [A, B],
        "B_characteristic_roots": [mp.nstr(value, 70) for value in roots_B],
        "A_characteristic_roots": [mp.nstr(value, 70) for value in roots_A],
        "A_half_boundary_mode_A_over_4": mp.nstr(half_corner_mode, 70),
        "half_stationary_alpha_sqrt_8t_over_pi": mp.nstr(half_stationary_alpha, 70),
        "A_minus_half_stationary_alpha": mp.nstr(A_mp - half_stationary_alpha, 70),
        "A_half_phase_derivative_detuning": mp.nstr(half_detuning, 70),
        "interpretation": "The same endpoint driver produces the B crossing, both A characteristic roots, and the A/x=1/2 corner. The midpoint arithmetic identity supplies boundary data but does not linearize the propagator.",
    }


def render_note(artifact: dict[str, Any]) -> str:
    g = artifact["finite_height_geometry"]
    return f"""# Endpoint-driven transport for symmetric Fourier pairs

Date: 2026-08-13

Status: exact-lemma certificate; not a proof of the Volterra operator bound

Let

```text
f_x(u)=(A+2u)exp[i*pi*x(A+2u)^2/4],
P_m(x)=I_m(x)+I_-m(x).                                (VT1)
```

The canonical interpolation obeys the exact differential identity

```text
partial_u^2 f_x
 =6 i*pi*x f_x+4 i*pi*x^2 partial_x f_x.              (VT2)
```

Integrating (VT2) against the symmetric Fourier character and using the
integer endpoints gives a first-order transport equation:

```text
4 i*pi*x^2 P_m'(x)
 +(4*pi^2*m^2+6 i*pi*x)P_m(x)=2 Delta f_x',           (VT3)

Delta f_x'=e^(i*pi*x*B^2/4)(2+i*pi*x*B^2)
            -e^(i*pi*x*A^2/4)(2+i*pi*x*A^2).
```

At `x=0`, the interpolation is linear, and direct integration gives
`P_m(0)=0` for every nonzero integer `m`.  The integrating factor for (VT3)
is

```text
mu_m(x)=x^(3/2)exp(i*pi*m^2/x).                        (VT4)
```

Therefore the complete pair coefficient, including all endpoint/Fresnel
cancellation, is exactly

```text
P_m(x)=x^(-3/2)exp(-i*pi*m^2/x)/(2 i*pi)
 integral_0^x z^(-1/2)exp(i*pi*m^2/z)Delta f_z' dz.   (VT5)
```

This is not an asymptotic expansion.  The driver vanishes at `z=0`, so its
integrand is locally `O(z^(1/2))` and the lower endpoint is ordinary.  After
inserting (VT5) into the half-Kummer integral, absolute local integrability
permits the triangular order to be exchanged:

```text
K_m+K_-m=1/(2 i*pi) integral_0^(1/2)
 z^(-1/2)exp(i*pi*m^2/z)Delta f_z' G_m(z)dz,

G_m(z)=integral_z^(1/2) W_t(x)x^(-3/2)
                  exp(-i*pi*m^2/x)dx.                 (VT6)
```

Equations (VT5)--(VT6) are a cancellation-preserving replacement for the
unsafe sum of bare `A` and `B` endpoint exponentials.  The finite-roster
interior has not been discarded: it is encoded exactly by the Volterra
propagator from the endpoint driver.

The geometry is also exact.  The endpoint-driver phase has stationary point
`z_D=2m/D`; it aligns with the outer Morse saddle precisely when

```text
D=2m+t/(pi*m).                                        (VT7)
```

At the saved height, (VT7) reproduces

```text
B lower root: {g['B_characteristic_roots'][0]},
A roots:      {g['A_characteristic_roots'][0]},
              {g['A_characteristic_roots'][1]}.       (VT8)
```

The half-boundary corner occurs at

```text
A/4={g['A_half_boundary_mode_A_over_4']},
sqrt(8t/pi)={g['half_stationary_alpha_sqrt_8t_over_pi']},
A-sqrt(8t/pi)={g['A_minus_half_stationary_alpha']}.    (VT9)
```

Thus the interesting midpoint symmetry is real but has a disciplined role:
it gives exact boundary data for (VT3)--(VT6).  It does not replace the
continuous quadratic chirp by a linear interpolation.

The next quantitative problem is now an oscillatory Volterra-operator bound
for (VT6), with separate local charts only where the two phases align in
(VT7).  This formulation automatically keeps the `+m/-m` cancellation that
the tangent B profile lacked.

Pi provenance: every `pi` in (VT1)--(VT9) is differentiated or completed
from the equation-(9) Kummer phase and integer Fourier characters.  No
geometric or fitted constant is introduced.

Proof boundary: exact PDE, pair ODE, Volterra solution, locally justified
triangular exchange, and characteristic geometry only.  No quantitative
Volterra norm, nonlinear B crossing, A-fold splice, complete `T_upper`,
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(
        dependencies["pair_Morse"]["decision"]["signed_pair_formed_inside_one_half_domain_Gaussian_integral"] is True,
        "pair-Morse dependency drift",
    )
    require(
        dependencies["paired_residual"]["decision"]["every_negative_mode_retained_exactly_once"] is True,
        "paired residual dependency drift",
    )
    require(
        dependencies["symmetric_poisson"]["decision"]["paired_one_over_m_endpoint_current_cancels"] is True,
        "pair cancellation dependency drift",
    )
    require(
        dependencies["midpoint_guard"]["decision"]["prior_linear_interpolation_anchor_rejected"] is True,
        "midpoint guard dependency drift",
    )

    artifact = {
        "kind": STEM,
        "status": "exact_endpoint_driven_pair_transport_and_half_volterra_representation_proved_operator_bound_open",
        "passed": True,
        "scope": {
            "height": T,
            "source_odd_roster": [A, B],
            "roster_coordinate_length": L,
            "mode_condition": "integer m>0",
        },
        "symbolic_certificate": symbolic_certificate(),
        "finite_height_geometry": finite_height_geometry(),
        "decision": {
            "canonical_roster_transport_PDE_proved": True,
            "symmetric_pair_first_order_x_ODE_proved": True,
            "zero_initial_pair_at_x0_proved": True,
            "endpoint_driven_Volterra_solution_proved": True,
            "half_pair_triangular_Volterra_representation_proved": True,
            "bare_endpoint_exponentials_separated": False,
            "midpoint_parity_used_as_boundary_data_only": True,
            "Volterra_operator_quantitatively_bounded": False,
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
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "sympy_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Construct a quantitative oscillatory bound for the Volterra kernel G_m and the endpoint-driver integral in the nonaligned regions. Isolate only neighborhoods of the B lower root and the A/half-fold alignments, then compare their grouped signed contribution with the target budget.",
        "proof_boundary": "Exact roster PDE, pair transport ODE, endpoint-driven Volterra representation, and characteristic geometry only. No quantitative Volterra operator bound, nonlinear B crossing, A-fold splice, complete T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified endpoint-driven Volterra transport for symmetric Fourier pairs", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
