#!/usr/bin/env python3
"""Bound two derivatives of the exact completed-B outer Fresnel remainder."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_outer_fresnel_remainder_derivative_tail_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "rational_tail": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_outer_three_current_derivative_tail_gate.json",
    "derivative_transport": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_derivative_transport_gate.json",
    "fresnel_remainder": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_gaussian_window_fresnel_remainder_gate.json",
}

PRECISION = 100
T = 10_000_000_000
B = 5_122_421
M = B
DELTA_TEXT = "1e-4"


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
    q = sp.symbols("q", nonzero=True, real=True)
    pi = sp.pi
    truncated = sp.I / (pi * q) + 1 / (pi**2 * q**3) - 3 * sp.I / (pi**3 * q**5)
    forcing = sp.simplify(-1 - sp.I * pi * q * truncated - sp.diff(truncated, q))
    require(forcing == -15 * sp.I / (pi**3 * q**6), "phase-stripped remainder forcing drift")

    remainder = sp.Function("r")(q)
    first = -sp.I * pi * q * remainder + forcing
    second = sp.simplify(
        -sp.I * pi * remainder
        - sp.I * pi * q * first
        + sp.diff(forcing, q)
    )
    expected_second = sp.simplify(
        -(sp.I * pi + pi**2 * q**2) * remainder
        - 15 / (pi**2 * q**5)
        + 90 * sp.I / (pi**3 * q**7)
    )
    require(sp.simplify(second - expected_second) == 0, "second remainder transport drift")

    x = sp.symbols("x", positive=True, real=True)
    f = sp.Function("f")(x)
    qx = sp.Function("q")(x)
    rx = sp.Function("r")
    product = f * rx(qx)
    first_chain = sp.diff(product, x)
    second_chain = sp.diff(product, x, 2)
    return {
        "phase_stripped_tail": "z(q)=exp(-i*pi*q^2/2)T(q), z'=-1-i*pi*q*z",
        "truncation": "P_3=i/(pi q)+1/(pi^2 q^3)-3i/(pi^3 q^5), r=z-P_3",
        "remainder_transport": "r'=-i*pi*q*r-15i/(pi^3 q^6)",
        "second_remainder_transport": "r''=-(i*pi+pi^2 q^2)r-15/(pi^2 q^5)+90i/(pi^3 q^7)",
        "bounds": "|r|<=30/(pi^4|q|^7), |r'|<=45/(pi^3|q|^6), |r''|<=45/(pi^2|q|^5)+120/(pi^3|q|^7)",
        "completed_remainder": "R_s(x)=f(x)r(q_s(x)), f=m sqrt(2)x^(-3/2), s in {+,-}",
        "first_chain_rule": str(first_chain),
        "second_chain_rule": str(second_chain),
        "outer_geometry": "|q_s|>=alpha sqrt(2/x)m, |q_s'|<=beta sqrt(2/x)m/(2x), |q_s''|<=gamma sqrt(2/x)m/(4x^2), with alpha=3/4, beta=5/4, gamma=13/4",
    }


def series_tail(exponent: int) -> arb:
    start = arb(M)
    return start ** (-exponent) + start ** (1 - exponent) / (exponent - 1)


def interval_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    alpha, beta, gamma = arb(3) / 4, arb(5) / 4, arb(13) / 4
    delta, x_max = arb(DELTA_TEXT), arb(1) / 2

    # One-sign bounds after inserting the uniform outer geometry.  The two
    # completed endpoint signs are finally combined by a factor of two.
    a0 = arb(15) / (4 * pi**4 * alpha**7)
    r0 = a0 * x_max**2 * series_tail(6)

    a11 = arb(45) * beta / (8 * pi**3 * alpha**6)
    r1 = (arb(3) * a0 / 2) * x_max * series_tail(6) + a11 * series_tail(4)

    r2 = (
        (arb(15) * a0 / 4) * series_tail(6)
        + (3 * a11) * delta**-1 * series_tail(4)
        + (arb(45) * beta**2 / (8 * pi**2 * alpha**5)) * delta**-2 * series_tail(2)
        + (arb(15) * beta**2 / (2 * pi**3 * alpha**7)) * delta**-1 * series_tail(4)
        + (arb(45) * gamma / (16 * pi**3 * alpha**6)) * delta**-1 * series_tail(4)
    )
    raw = {0: 2 * r0, 1: 2 * r1, 2: 2 * r2}

    endpoint, t = arb(B), arb(T)
    x0 = (1 - (1 - 8 * t / (pi * endpoint**2)).sqrt()) / 2
    hessian = t * (1 - 2 * x0) / (2 * x0**2 * (1 - x0) ** 2)
    root_hessian = hessian.sqrt()
    weight = (delta * (1 - delta)) ** (-arb(1) / 4)
    log_first = (1 - 2 * delta) / (4 * delta * (1 - delta))
    second_ratio = (12 * delta**2 - 12 * delta + 5) / (16 * delta**2 * (1 - delta) ** 2)
    weighted = {
        0: weight * raw[0],
        1: weight * (raw[1] + log_first * raw[0]),
        2: weight * (raw[2] + 2 * log_first * raw[1] + second_ratio * raw[0]),
    }
    normalized = {
        0: weighted[0],
        1: weighted[1] / root_hessian,
        2: weighted[2] / hessian,
    }
    require(all(value.is_finite() and value > 0 for value in normalized.values()), "nonfinite remainder derivative bound")

    encode = lambda values: {str(order): value.str(PRECISION, more=True) for order, value in values.items()}
    return {
        "height": T,
        "endpoint": B,
        "outer_mode_start": M,
        "x_domain": [DELTA_TEXT, "1/2"],
        "outer_geometry_constants": {
            "alpha": alpha.str(PRECISION, more=True),
            "beta": beta.str(PRECISION, more=True),
            "gamma": gamma.str(PRECISION, more=True),
        },
        "two_sign_remainder_x_derivative_tail_balls": encode(raw),
        "weighted_x_derivative_tail_balls": encode(weighted),
        "weighted_normalized_derivative_tail_balls": encode(normalized),
        "sqrt_trace_hessian_ball": root_hessian.str(PRECISION, more=True),
        "derivative_order_legend": {"0": "remainder amplitude", "1": "first derivative", "2": "second derivative"},
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    raw = c["two_sign_remainder_x_derivative_tail_balls"]
    norm = c["weighted_normalized_derivative_tail_balls"]
    return f"""# Outer differentiated Fresnel remainder for the completed B trace

Date: 2026-08-13

Status: exact remainder-transport and uniform outer derivative-tail
certificate; not a proof of the complete completed-B estimate

Write the sign-adapted Fresnel tail as

```text
z(q)=exp(-i*pi*q^2/2)T(q),
P_3(q)=i/(pi q)+1/(pi^2 q^3)-3i/(pi^3 q^5),
r(q)=z(q)-P_3(q).                                      (OR1)
```

The already certified contour remainder gives
`|r|<=30/(pi^4|q|^7)`.  Differentiating the exact tail equation before
taking norms gives

```text
r' =-i*pi*q*r-15i/(pi^3 q^6),
r''=-(i*pi+pi^2 q^2)r-15/(pi^2 q^5)
     +90i/(pi^3 q^7),                                  (OR2)

|r'| <=45/(pi^3|q|^6),
|r''|<=45/(pi^2|q|^5)+120/(pi^3|q|^7).                (OR3)
```

For either completed outer sign, `R_s=f r(q_s)` with
`f=m sqrt(2)x^(-3/2)`.  On `m>=B={B}` and
`{DELTA_TEXT}<=x<=1/2`,

```text
|q_s|  >=(3/4)sqrt(2/x)m,
|q_s'| <=(5/4)sqrt(2/x)m/(2x),
|q_s''|<=(13/4)sqrt(2/x)m/(4x^2).                     (OR4)
```

Applying the product/chain rules and summing both signs by the decreasing
series integral test gives

```text
sum_s sum_(m>=B)|R_s|    <= {raw['0']}
sum_s sum_(m>=B)|R_s'|   <= {raw['1']}
sum_s sum_(m>=B)|R_s''|  <= {raw['2']}.                (OR5)
```

After the Kummer weight and normalized coordinate are included,

```text
|w R|          <= {norm['0']}
|d_xi(w R)|    <= {norm['1']}
|d_xixi(w R)|  <= {norm['2']}.                         (OR6)
```

Thus the asymptotic remainder is not being dropped: its complete outer tail
is uniformly controlled through the two derivatives needed for tangential
integration by parts.  The companion rational-current gate and this gate do
not yet close the analytic block: the exact Fresnel-to-boundary dictionary
term `D_m` from Section 11.411 must also be differentiated and summed.

The dictionary term, finite block `m<B`, local phase-coupled tangential
integration, window and endpoint boundary terms, and the common-regulator
`R_join` combination remain open.  This gate establishes none of a complete
B estimate, `T_upper`, a height-uniform theorem, `Lambda<=0`, PF-infinity,
RH, or a prize-level conclusion.

Pi provenance: `pi` is the same equation-(9) Fresnel phase constant used in
the exact tail differential equation; no fitted or decorative value enters.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(
        dependencies["rational_tail"]["decision"]["outer_rational_amplitude_first_second_derivative_bounds_uniform_in_cutoff_and_Abel_weight"] is True,
        "rational-tail dependency drift",
    )
    require(
        dependencies["fresnel_remainder"]["decision"]["three_term_Fresnel_tail_with_explicit_remainder_proved"] is True,
        "Fresnel-remainder dependency drift",
    )

    artifact = {
        "kind": STEM,
        "status": "uniform_completed_B_outer_exact_Fresnel_remainder_tail_bounded_through_two_normalized_derivatives",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(),
        "decision": {
            "exact_phase_stripped_remainder_first_second_transport_proved": True,
            "two_sign_outer_Fresnel_remainder_tail_analytically_summed": True,
            "outer_remainder_amplitude_first_second_derivative_bounds_uniform_in_cutoff_and_Abel_weight": True,
            "Fresnel_to_boundary_dictionary_derivative_tail_bounded": False,
            "complete_outer_mode_derivative_block_bounded": False,
            "finite_mode_block_below_outer_start_bounded": False,
            "completed_exterior_current_bound_proved": False,
            "joined_remainder_below_1p4058e_minus_4_proved": False,
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
            "flint_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Differentiate and sum the exact Fresnel-to-boundary dictionary D_m on m>=B. Only then combine the complete outer block with the local tangential phase rather than a global amplitude supremum. In parallel replace the finite m<B derivative roster by the original finite-Poisson or Dirichlet representation before taking norms, retaining every window and endpoint boundary term.",
        "proof_boundary": "Uniform analytic two-sign outer Fresnel-remainder tail bounds through two derivatives on m>=B and delta<=x<=1/2, with the rational block available separately. No differentiated Fresnel-to-boundary dictionary tail, complete outer derivative block, finite m<B block, local phase-coupled exterior integral, complete exterior-B or joined-remainder bound, complete T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified completed-B outer exact Fresnel remainder derivative tail", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
