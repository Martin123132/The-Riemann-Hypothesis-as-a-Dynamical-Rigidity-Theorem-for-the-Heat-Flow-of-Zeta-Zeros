#!/usr/bin/env python3
"""Reduce the full exact A exterior to one common-phase integral."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_exact_exterior_common_phase_reduction_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "triangle": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_separable_triangle_geometry_gate.json",
    "orientation": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_corner_projector_orientation_gate.json",
    "regrouping": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_A_localized_exact_domain_regrouping_gate.json",
    "cancellation_guard": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_exterior_asymptotic_cancellation_guard_gate.json",
}

HEIGHT = 10_000_000_000
A = 159_577
LOWER_START = 39_853
UPPER_END = 39_936
MODES = tuple(range(LOWER_START, UPPER_END + 1))
Y0_TEXT = "0.0037"


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


def complex_record(value: acb, digits: int = 70) -> dict[str, str]:
    return {
        "real_ball": value.real.str(digits, more=True),
        "imag_ball": value.imag.str(digits, more=True),
        "absolute_ball": abs(value).str(digits, more=True),
    }


def symbolic_certificate() -> dict[str, str]:
    alpha, beta, X = sp.symbols("alpha beta X", positive=True)
    w_minus = beta / X - alpha * X
    w_plus = beta / X + alpha * X
    gaussian_primitive = sp.sqrt(sp.pi) / (4 * alpha) * (
        sp.exp(-2 * alpha * beta) * sp.erfc(w_minus)
        - sp.exp(2 * alpha * beta) * sp.erfc(w_plus)
    )
    require(
        sp.simplify(
            sp.diff(gaussian_primitive, X)
            - sp.exp(-alpha**2 * X**2 - beta**2 / X**2)
        ) == 0,
        "inverse-Gaussian primitive failed",
    )
    F0 = 2 * gaussian_primitive
    driver = sp.simplify(2 * sp.diff(alpha * F0, alpha))
    expected = (
        -2 * sp.sqrt(sp.pi) * beta * (
            sp.exp(-2 * alpha * beta) * sp.erfc(w_minus)
            + sp.exp(2 * alpha * beta) * sp.erfc(w_plus)
        )
        + 2 * X * (
            sp.exp(-2 * alpha * beta - w_minus**2)
            + sp.exp(2 * alpha * beta - w_plus**2)
        )
    )
    require(sp.simplify(driver - expected) == 0, "endpoint driver primitive failed")

    x, m, endpoint = sp.symbols("x m A", positive=True, real=True)
    phase_endpoint = sp.pi * m**2 / x + sp.pi * endpoint**2 * x / 4
    phase_outer = -sp.pi * m**2 / x + sp.symbols("t", positive=True) / 2 * sp.log((1 - x) / x)
    common = sp.simplify(phase_endpoint + phase_outer)
    require(not common.has(m), "mode did not cancel from exact face phase")

    return {
        "base_primitive": "int_0^X exp(-alpha^2 q^2-beta^2/q^2)dq=sqrt(pi)/(4alpha){exp(-2alpha beta)erfc(beta/X-alpha X)-exp(2alpha beta)erfc(beta/X+alpha X)}",
        "physical_parameters": "alpha=exp(-i*pi/4)sqrt(pi)A/2, beta=exp(-i*pi/4)sqrt(pi)m, X=sqrt(x)",
        "endpoint_driver_primitive": "H_m(x)=4sqrt(x)exp(iF_m(x))-2sqrt(pi)beta(-1)^(mA){erfc(w_-)+erfc(w_+)}",
        "endpoint_phase": "F_m(x)=pi*m^2/x+pi*A^2*x/4",
        "phase_stripped_ODE": "B_m'+iF_m'B_m=h_A, B_m=exp(-iF_m)H_m, h_A=x^(-1/2)(2+i*pi*A^2*x)",
        "endpoint_IBP_recurrence": "b_0=h_A/(iF_m'), b_(j+1)=-b_j'/(iF_m'), B_m=sum_(j<N)b_j-exp(-iF_m)int_0^x b_(N-1)'exp(iF_m)",
        "common_face_phase": "Phi_A(x)=pi*A^2*x/4+(t/2)log((1-x)/x)",
        "exact_exterior": "C_A,ext=-1/(2pi i)sum_m int_0^(x0_m) x^(-7/4)(1-x)^(-1/4)B_m(x)exp(iPhi_A(x))dx",
        "cutoff": "x0_m={1+[t/(2pi*m^2)](1+y0)}^(-1)",
    }


def endpoint_terms(x: arb, mode: arb, endpoint: arb) -> list[acb]:
    pi = arb.pi()
    imaginary = acb(0, 1)
    denominator = (endpoint * x - 2 * mode) * (endpoint * x + 2 * mode)
    b0 = 4 * x ** arb("1.5") * (pi * endpoint**2 * x - 2 * imaginary) / (pi * denominator)
    b1 = (
        8 * imaginary * x ** arb("2.5")
        * (
            pi * endpoint**4 * x**3
            - 20 * pi * endpoint**2 * mode**2 * x
            + 2 * imaginary * endpoint**2 * x**2
            + 24 * imaginary * mode**2
        )
        / (pi**2 * denominator**3)
    )
    b2 = (
        16 * x ** arb("3.5")
        * (
            pi * endpoint**6 * x**5
            - 56 * pi * endpoint**4 * mode**2 * x**3
            + 6 * imaginary * endpoint**4 * x**4
            - 560 * pi * endpoint**2 * mode**4 * x
            + 240 * imaginary * endpoint**2 * mode**2 * x**2
            + 480 * imaginary * mode**4
        )
        / (pi**3 * denominator**5)
    )
    b3 = (
        -96 * imaginary * x ** arb("4.5")
        * (
            pi * endpoint**8 * x**7
            - 108 * pi * endpoint**6 * mode**2 * x**5
            + 10 * imaginary * endpoint**6 * x**6
            - 3024 * pi * endpoint**4 * mode**4 * x**3
            + 840 * imaginary * endpoint**4 * mode**2 * x**4
            - 6720 * pi * endpoint**2 * mode**6 * x
            + 5600 * imaginary * endpoint**2 * mode**4 * x**2
            + 4480 * imaginary * mode**6
        )
        / (pi**4 * denominator**7)
    )
    return [b0, b1, b2, b3]


def roster_certificate(precision: int = 100) -> dict[str, Any]:
    ctx.dps = precision
    ctx.threads = 1
    pi, t, endpoint = arb.pi(), arb(HEIGHT), arb(A)
    y0 = arb(Y0_TEXT)
    kappa = (1 - 8 * t / (pi * endpoint**2)).sqrt()
    x_star = (1 - kappa) / 2
    phase = lambda x: pi * endpoint**2 * x / 4 + t / 2 * ((1 - x) / x).log()
    rows: list[dict[str, Any]] = []
    maxima = [arb(0) for _ in range(4)]
    positive_modes: list[int] = []

    for mode_int in MODES:
        mode = arb(mode_int)
        r = t / (2 * pi * mode**2)
        x0 = 1 / (1 + r * (1 + y0))
        delta = 2 * (phase(x0) - phase(x_star))
        require(delta.lower() > 0, f"Morse phase gap not positive at mode {mode_int}")
        w_abs = delta.sqrt()
        if x0.lower() > x_star.upper():
            w0 = w_abs
            positive_modes.append(mode_int)
        else:
            require(x0.upper() < x_star.lower(), f"cutoff straddles saddle at mode {mode_int}")
            w0 = -w_abs
        terms = endpoint_terms(x0, mode, endpoint)
        for index, term in enumerate(terms):
            maxima[index] = max(maxima[index], abs(term).upper())
        rows.append(
            {
                "mode": mode_int,
                "x0_ball": x0.str(75, more=True),
                "stationary_Morse_endpoint_w0_ball": w0.str(70, more=True),
                "endpoint_IBP_terms_at_x0": [complex_record(term) for term in terms],
            }
        )

    require(positive_modes == list(range(39_927, 39_937)), "positive stationary roster drift")
    require(maxima[0] < arb("765"), "b0 cutoff cap lost")
    require(maxima[1] < arb("0.012"), "b1 cutoff cap lost")
    require(maxima[2] < arb("5e-7"), "b2 cutoff cap lost")
    require(maxima[3] < arb("3.7e-11"), "b3 cutoff cap lost")
    return {
        "height": HEIGHT,
        "endpoint": A,
        "mode_range": [LOWER_START, UPPER_END],
        "mode_count": len(rows),
        "y0": Y0_TEXT,
        "x_star_ball": x_star.str(75, more=True),
        "positive_stationary_endpoint_modes": positive_modes,
        "w0_min_ball": arb(rows[0]["stationary_Morse_endpoint_w0_ball"]).str(70, more=True),
        "w0_max_ball": arb(rows[-1]["stationary_Morse_endpoint_w0_ball"]).str(70, more=True),
        "max_abs_endpoint_IBP_terms_at_cutoff": [value.str(70, more=True) for value in maxima],
        "rows": rows,
    }


def diagnostic_pilot(dps: int, gauss_nodes: int, outer_terms: int) -> dict[str, Any]:
    import mpmath as mp
    import numpy as np

    mp.mp.dps = dps
    t, endpoint, pi = mp.mpf(HEIGHT), mp.mpf(A), mp.pi
    y0 = mp.mpf(Y0_TEXT)
    c = pi * endpoint**2 / 4
    kappa = mp.sqrt(1 - 8 * t / (pi * endpoint**2))
    x_star = (1 - kappa) / 2
    phase_star = c * x_star + t / 2 * mp.log((1 - x_star) / x_star)

    def phase_delta(x):
        return c * (x - x_star) + t / 2 * (
            mp.log((1 - x) / (1 - x_star)) - mp.log(x / x_star)
        )

    def phase_prime(x):
        return c - t / (2 * x * (1 - x))

    x_left = mp.findroot(
        lambda x: phase_delta(x) - 50,
        (x_star - mp.mpf("0.0015"), x_star - mp.mpf("0.0007")),
    )

    def endpoint_tail(x, mode):
        denominator = (endpoint * x - 2 * mode) * (endpoint * x + 2 * mode)
        b0 = 4 * x**mp.mpf("1.5") * (pi * endpoint**2 * x - 2j) / (pi * denominator)
        b1 = (
            8j * x**mp.mpf("2.5")
            * (
                pi * endpoint**4 * x**3
                - 20 * pi * endpoint**2 * mode**2 * x
                + 2j * endpoint**2 * x**2
                + 24j * mode**2
            )
            / (pi**2 * denominator**3)
        )
        b2 = (
            16 * x**mp.mpf("3.5")
            * (
                pi * endpoint**6 * x**5
                - 56 * pi * endpoint**4 * mode**2 * x**3
                + 6j * endpoint**4 * x**4
                - 560 * pi * endpoint**2 * mode**4 * x
                + 240j * endpoint**2 * mode**2 * x**2
                + 480j * mode**4
            )
            / (pi**3 * denominator**5)
        )
        return b0 + b1 + b2

    def amplitude(x, mode):
        return (
            -1 / (2j * pi)
            * x ** (-mp.mpf(7) / 4)
            * (1 - x) ** (-mp.mpf(1) / 4)
            * endpoint_tail(x, mode)
        )

    nodes, weights = np.polynomial.legendre.leggauss(gauss_nodes)
    near_relative = 0j
    for mode_int in MODES:
        mode = mp.mpf(mode_int)
        r = t / (2 * pi * mode**2)
        x0 = 1 / (1 + r * (1 + y0))
        midpoint = (x_left + x0) / 2
        half_width = (x0 - x_left) / 2
        subtotal = 0j
        for node, weight in zip(nodes, weights):
            x = midpoint + half_width * mp.mpf(str(node))
            subtotal += (
                mp.mpf(str(weight))
                * amplitude(x, mode)
                * mp.exp(1j * phase_delta(x))
            )
        near_relative += half_width * subtotal

    def grouped_amplitude(x):
        return sum(amplitude(x, mp.mpf(mode)) for mode in MODES)

    hierarchy = [lambda x: grouped_amplitude(x) / (1j * phase_prime(x))]
    for _ in range(1, outer_terms):
        previous = hierarchy[-1]
        hierarchy.append(
            lambda x, prior=previous: -mp.diff(prior, x) / (1j * phase_prime(x))
        )
    boundary_terms = [function(x_left) for function in hierarchy]
    far_partial = []
    running = 0j
    for value in boundary_terms:
        running += value
        far_partial.append(mp.exp(50j) * running)

    normalizer = (pi / (32 * t)) ** mp.mpf("0.25")
    paper_rotation = mp.exp(-1j * pi / 8)
    carrier = mp.exp(1j * phase_star)
    physical = [
        2 * normalizer * mp.re(paper_rotation * carrier * (near_relative + far))
        for far in far_partial
    ]
    return {
        "rigorous_certificate": False,
        "purpose": "route-selection telemetry only; Gauss quadrature, endpoint-tail truncation, and outer-IBP remainder are not interval bounded",
        "decimal_digits": dps,
        "gauss_legendre_nodes_per_mode": gauss_nodes,
        "outer_boundary_terms": outer_terms,
        "x_left_w_minus_10": mp.nstr(x_left, 45),
        "near_relative_complex": [mp.nstr(mp.re(near_relative), 45), mp.nstr(mp.im(near_relative), 45)],
        "outer_boundary_term_absolute_values": [mp.nstr(abs(value), 35) for value in boundary_terms],
        "physical_partial_estimates": [mp.nstr(value, 35) for value in physical],
        "final_physical_scale_pilot": mp.nstr(physical[-1], 35),
    }


def render_note(artifact: dict[str, Any]) -> str:
    r = artifact["roster_certificate"]
    p = artifact["diagnostic_pilot"]
    return f"""# Full exact A-exterior common-phase reduction

Date: 2026-08-14

Status: exact one-dimensional reduction and certified cutoff geometry, with
nonrigorous scale telemetry; not an exterior value or complete A bound

For

```text
F_m(z)=pi m^2/z+pi A^2 z/4,
h_A(z)=z^(-1/2)(2+i pi A^2 z),
H_m(x)=integral_0^x h_A(z)exp(iF_m(z))dz,             (XR1)
```

put `alpha=exp(-i*pi/4)sqrt(pi)A/2`,
`beta=exp(-i*pi/4)sqrt(pi)m`, `X=sqrt(x)`, and
`w_+/-=beta/X +/- alpha X`.  Differentiating the exact inverse-Gaussian
erfc primitive gives

```text
H_m(x)=4sqrt(x)exp(iF_m(x))
 -2sqrt(pi)beta(-1)^(mA)[erfc(w_-)+erfc(w_+)].        (XR2)
```

The phase-stripped endpoint tail `B_m=exp(-iF_m)H_m` obeys

```text
B_m'+iF_m'B_m=h_A.                                   (XR3)
```

Since every exterior cutoff lies before the endpoint saddle, define

```text
b_0=h_A/(iF_m'),
b_(j+1)=-b_j'/(iF_m').                               (XR4)
```

This gives an exact finite hierarchy plus one explicit oscillatory
remainder.  At the 84 cutoff points, the maxima of `|b_0|,...,|b_3|` are

```text
{r['max_abs_endpoint_IBP_terms_at_cutoff'][0]},
{r['max_abs_endpoint_IBP_terms_at_cutoff'][1]},
{r['max_abs_endpoint_IBP_terms_at_cutoff'][2]},
{r['max_abs_endpoint_IBP_terms_at_cutoff'][3]}.       (XR5)
```

The endpoint phase cancels the mode-dependent outer term exactly.  Hence the
full convergent exterior is the single common-phase family

```text
C_A,ext^exact=-1/(2pi i) sum_(m=39853)^39936
 integral_0^(x0_m) x^(-7/4)(1-x)^(-1/4)B_m(x)
 exp(i Phi_A(x))dx,                                  (XR6)

Phi_A(x)=pi A^2 x/4+(t/2)log((1-x)/x),
x0_m={{1+[t/(2pi m^2)](1+0.0037)}}^(-1).             (XR7)
```

Its unique stationary point is
`x_*={r['x_star_ball']}`.  The signed Morse endpoint
`w0=sgn(x0-x_*)sqrt(2[Phi_A(x0)-Phi_A(x_*)])` runs from
`{r['w0_min_ball']}` to `{r['w0_max_ball']}` and is positive exactly for
modes `39927..39936`.

Diagnostic only: a 180-node-per-mode Gauss pilot on `w>=-10`, a three-term
endpoint hierarchy, and {p['outer_boundary_terms']} outer boundary terms
give physical partial estimates

```text
{chr(10).join(p['physical_partial_estimates'])}       (XR8)
```

The final displayed scale, `{p['final_physical_scale_pilot']}`, is not an
interval enclosure.  It only shows that the exact exterior is likely a
signed carrier component at a much larger scale than the final `R_Dir`
budget.  Promotion requires interval bounds for the endpoint-hierarchy
remainder, stationary-region quadrature, and the final outer-IBP remainder.

Pi provenance: every `pi` in (XR1)--(XR8) comes from the original
equation-(9) triangle, exact endpoint/outer phases, and paper normalization.
No fitted constant is introduced.

Proof boundary: exact erfc primitive, phase-stripped ODE, one-dimensional
common-phase reduction, cutoff/Morse roster, and pointwise cutoff hierarchy
only.  The diagnostic pilot is not a certificate.  No exact exterior value,
uniform hierarchy remainder, outer-tail bound, compact nonlinear-amplitude
bound, complete A theorem, `R_Dir` estimate, RH, or prize-level conclusion is
proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["triangle"]["decision"]["pair_Volterra_integral_reduced_to_separable_triangle"] is True, "triangle dependency drift")
    require(dependencies["orientation"]["decision"]["A_endpoint_sign_identified_as_epsilon_A_minus_one"] is True, "orientation dependency drift")
    require(dependencies["regrouping"]["decision"]["ten_true_exact_face_stationary_transitions_remain_in_exterior"] is True, "stationary roster dependency drift")
    require(dependencies["cancellation_guard"]["decision"]["common_cutoff_recombination_required_before_exterior_limit"] is True, "cancellation dependency drift")

    artifact = {
        "kind": STEM,
        "status": "full_exact_A_exterior_reduced_to_common_phase_endpoint_tail_family_rigorous_remainder_and_value_open",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "roster_certificate": roster_certificate(),
        "diagnostic_pilot": diagnostic_pilot(dps=60, gauss_nodes=180, outer_terms=5),
        "decision": {
            "exact_endpoint_driver_erfc_primitive_proved": True,
            "phase_stripped_endpoint_tail_ODE_proved": True,
            "full_exact_exterior_common_phase_reduction_proved": True,
            "mode_dependent_face_phase_cancelled_exactly": True,
            "ten_mode_incomplete_stationary_roster_certified": True,
            "endpoint_IBP_hierarchy_evaluated_at_all_cutoffs": True,
            "diagnostic_exterior_scale_is_rigorous": False,
            "uniform_endpoint_hierarchy_remainder_bound_proved": False,
            "outer_IBP_remainder_bound_proved": False,
            "exact_exterior_value_proved": False,
            "complete_A_endpoint_block_proved": False,
            "R_Dir_bound_proved": False,
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
            "arb_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Certify the endpoint-tail hierarchy remainder uniformly on the stationary window, perform interval quadrature in the exact Morse coordinate, and bound the grouped far current after enough outer integrations by parts to make its final remainder smaller than the A-assembly budget.",
        "proof_boundary": "Exact erfc primitive, phase-stripped ODE, common-phase reduction, cutoff/Morse roster, and pointwise cutoff hierarchy only. Diagnostic quadrature is not a certificate. No exact exterior value, uniform remainder, compact nonlinear-amplitude bound, complete A endpoint theorem, R_Dir estimate, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified full exact A-exterior common-phase reduction", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
