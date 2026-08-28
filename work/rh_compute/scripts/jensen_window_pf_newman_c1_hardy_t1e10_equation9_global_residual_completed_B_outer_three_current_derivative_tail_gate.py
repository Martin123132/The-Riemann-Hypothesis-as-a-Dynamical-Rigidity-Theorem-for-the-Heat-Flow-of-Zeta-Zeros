#!/usr/bin/env python3
"""Bound two derivatives of the completed-B outer three-current rational tail."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_outer_three_current_derivative_tail_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "derivative_transport": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_derivative_transport_gate.json",
    "higher_currents": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_higher_boundary_current_absolute_window_gate.json",
    "fresnel_remainder": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_gaussian_window_fresnel_remainder_gate.json",
    "tangential_cutoff": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_tangential_exterior_cutoff_gate.json",
}

PRECISION = 100
T = 10_000_000_000
B = 5_122_421
DELTA_TEXT = "1e-4"
OUTER_START = B

# Each row is (coefficient numerator, pi denominator power, B power,
# c power, denominator power) for |alpha B^q c^p(c^2-m^2)^(-k)|.
CURRENT_TERMS: dict[str, tuple[tuple[int, int, int, int, int], ...]] = {
    "C0": (
        (1, 2, 0, 0, 1),
        (1, 1, 1, 1, 1),
    ),
    "C1": (
        (4, 2, 0, 4, 3),
        (4, 3, -1, 3, 3),
        (5, 2, 0, 2, 2),
        (3, 3, -1, 1, 2),
    ),
    "C2": (
        (48, 3, -1, 7, 5),
        (48, 4, -2, 6, 5),
        (84, 3, -1, 5, 4),
        (60, 4, -2, 4, 4),
        (35, 3, -1, 3, 3),
        (15, 4, -2, 2, 3),
    ),
}


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


def symbolic_certificate() -> dict[str, Any]:
    c, endpoint, mode = sp.symbols("c B m", positive=True, real=True)
    denominator = c**2 - mode**2
    pi = sp.pi

    c0 = -(1 + sp.I * pi * endpoint * c) / (pi**2 * denominator)
    c1 = c / (pi**3 * endpoint) * (
        (-4 * pi * endpoint * c**3 + 4 * sp.I * c**2) * denominator**-3
        + (5 * pi * endpoint * c - 3 * sp.I) * denominator**-2
    )
    c2 = -sp.I * c**2 / (pi**4 * endpoint**2) * (
        (-48 * pi * endpoint * c**5 + 48 * sp.I * c**4) * denominator**-5
        + (84 * pi * endpoint * c**3 - 60 * sp.I * c**2) * denominator**-4
        + (-35 * pi * endpoint * c + 15 * sp.I) * denominator**-3
    )

    signed_terms = {
        "C0": -denominator**-1 / pi**2 - sp.I * endpoint * c * denominator**-1 / pi,
        "C1": (
            -4 * c**4 * denominator**-3 / pi**2
            + 4 * sp.I * c**3 * denominator**-3 / (pi**3 * endpoint)
            + 5 * c**2 * denominator**-2 / pi**2
            - 3 * sp.I * c * denominator**-2 / (pi**3 * endpoint)
        ),
        "C2": (
            48 * sp.I * c**7 * denominator**-5 / (pi**3 * endpoint)
            + 48 * c**6 * denominator**-5 / (pi**4 * endpoint**2)
            - 84 * sp.I * c**5 * denominator**-4 / (pi**3 * endpoint)
            - 60 * c**4 * denominator**-4 / (pi**4 * endpoint**2)
            + 35 * sp.I * c**3 * denominator**-3 / (pi**3 * endpoint)
            + 15 * c**2 * denominator**-3 / (pi**4 * endpoint**2)
        ),
    }
    for name, exact in (("C0", c0), ("C1", c1), ("C2", c2)):
        require(sp.simplify(exact - signed_terms[name]) == 0, f"{name} monomial decomposition drift")

    p, k = sp.symbols("p k", integer=True, nonnegative=True)
    f = c**p * denominator**-k
    first = sp.diff(f, c)
    second = sp.diff(f, c, 2)
    expected_first = p * c ** (p - 1) * denominator**-k - 2 * k * c ** (p + 1) * denominator ** (-k - 1)
    expected_second = (
        p * (p - 1) * c ** (p - 2) * denominator**-k
        - 2 * k * (2 * p + 1) * c**p * denominator ** (-k - 1)
        + 4 * k * (k + 1) * c ** (p + 2) * denominator ** (-k - 2)
    )
    require(sp.simplify(first - expected_first) == 0, "first monomial derivative drift")
    require(sp.simplify(second - expected_second) == 0, "second monomial derivative drift")

    x = sp.symbols("x", positive=True, real=True)
    weight = (x * (1 - x)) ** (-sp.Rational(1, 4))
    log_first = sp.factor(sp.diff(weight, x) / weight)
    second_ratio = sp.factor(sp.diff(weight, x, 2) / weight)
    return {
        "current_formulas": {"C0": str(c0), "C1": str(c1), "C2": str(c2)},
        "monomial_derivatives": {
            "first": "d[c^p D^-k]/dc=p c^(p-1)D^-k-2k c^(p+1)D^(-k-1)",
            "second": "d2[c^p D^-k]/dc2=p(p-1)c^(p-2)D^-k-2k(2p+1)c^pD^(-k-1)+4k(k+1)c^(p+2)D^(-k-2)",
        },
        "x_conversion": "c=Bx/2, so d/dx=(B/2)d/dc and d2/dx2=(B/2)^2 d2/dc2",
        "outer_margin": "m>=B and 0<=c<=B/4 imply m^2-c^2>=(15/16)m^2",
        "tail_sum": "sum_(m=M)^infinity m^(-2k)<=M^(-2k)+M^(1-2k)/(2k-1)",
        "weight_log_first": str(log_first),
        "weight_second_ratio": str(second_ratio),
    }


def derivative_components(order: int, p: int, k: int) -> tuple[tuple[int, int, int], ...]:
    if order == 0:
        return ((1, p, k),)
    if order == 1:
        rows: list[tuple[int, int, int]] = []
        if p >= 1:
            rows.append((p, p - 1, k))
        rows.append((2 * k, p + 1, k + 1))
        return tuple(rows)
    require(order == 2, "only derivative orders zero through two are supported")
    rows = []
    if p >= 2:
        rows.append((p * (p - 1), p - 2, k))
    rows.append((2 * k * (2 * p + 1), p, k + 1))
    rows.append((4 * k * (k + 1), p + 2, k + 2))
    return tuple(rows)


def power(base: arb, exponent: int) -> arb:
    return base**exponent if exponent >= 0 else arb(1) / base ** (-exponent)


def zeta_tail_majorant(start: arb, denominator_power: int) -> arb:
    exponent = 2 * denominator_power
    return start ** (-exponent) + start ** (1 - exponent) / (exponent - 1)


def monomial_bound(
    coefficient_numerator: int,
    pi_power: int,
    endpoint_power: int,
    c_power: int,
    denominator_power: int,
    derivative_order: int,
    endpoint: arb,
    c_max: arb,
    rho: arb,
    start: arb,
) -> arb:
    coefficient = arb(coefficient_numerator) / arb.pi() ** pi_power
    total = arb(0)
    for factor, new_c_power, new_denominator_power in derivative_components(
        derivative_order, c_power, denominator_power
    ):
        total += (
            coefficient
            * power(endpoint, endpoint_power)
            * factor
            * c_max**new_c_power
            * rho ** (-new_denominator_power)
            * zeta_tail_majorant(start, new_denominator_power)
        )
    return total * (endpoint / 2) ** derivative_order


def interval_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    endpoint = arb(B)
    start = arb(OUTER_START)
    delta = arb(DELTA_TEXT)
    c_max = endpoint / 4
    rho = 1 - (c_max / start) ** 2
    require(rho == arb(15) / 16, "outer denominator margin drift")

    grouped: dict[str, dict[str, arb]] = {}
    totals = {order: arb(0) for order in range(3)}
    for current, terms in CURRENT_TERMS.items():
        current_rows: dict[str, arb] = {}
        for order in range(3):
            value = arb(0)
            for term in terms:
                value += monomial_bound(*term, order, endpoint, c_max, rho, start)
            current_rows[str(order)] = value
            totals[order] += value
        grouped[current] = current_rows

    t, pi = arb(T), arb.pi()
    x0 = (1 - (1 - 8 * t / (pi * endpoint**2)).sqrt()) / 2
    hessian = t * (1 - 2 * x0) / (2 * x0**2 * (1 - x0) ** 2)
    root_hessian = hessian.sqrt()

    weight_max = (delta * (1 - delta)) ** (-arb(1) / 4)
    weight_log_first_max = (1 - 2 * delta) / (4 * delta * (1 - delta))
    weight_second_ratio_max = (
        12 * delta**2 - 12 * delta + 5
    ) / (16 * delta**2 * (1 - delta) ** 2)
    weighted = {
        0: weight_max * totals[0],
        1: weight_max * (totals[1] + weight_log_first_max * totals[0]),
        2: weight_max * (
            totals[2]
            + 2 * weight_log_first_max * totals[1]
            + weight_second_ratio_max * totals[0]
        ),
    }
    normalized = {
        0: weighted[0],
        1: weighted[1] / root_hessian,
        2: weighted[2] / hessian,
    }
    require(all(value.is_finite() and value > 0 for value in normalized.values()), "nonfinite derivative tail bound")

    def encode(values: dict[Any, arb]) -> dict[str, str]:
        return {str(key): value.str(PRECISION, more=True) for key, value in values.items()}

    return {
        "height": T,
        "endpoint": B,
        "outer_mode_start": OUTER_START,
        "x_domain": [DELTA_TEXT, "1/2"],
        "c_max_ball": c_max.str(PRECISION, more=True),
        "c_over_m_max_ball": (c_max / start).str(PRECISION, more=True),
        "denominator_rho_ball": rho.str(PRECISION, more=True),
        "trace_root_ball": x0.str(PRECISION, more=True),
        "sqrt_trace_hessian_ball": root_hessian.str(PRECISION, more=True),
        "current_x_derivative_tail_balls": {
            name: encode(values) for name, values in grouped.items()
        },
        "summed_three_current_x_derivative_tail_balls": encode(totals),
        "Kummer_weight_max_ball": weight_max.str(PRECISION, more=True),
        "Kummer_log_first_abs_max_ball": weight_log_first_max.str(PRECISION, more=True),
        "Kummer_second_ratio_max_ball": weight_second_ratio_max.str(PRECISION, more=True),
        "weighted_x_derivative_tail_balls": encode(weighted),
        "weighted_normalized_derivative_tail_balls": encode(normalized),
        "derivative_order_legend": {"0": "amplitude", "1": "first derivative", "2": "second derivative"},
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    raw = c["summed_three_current_x_derivative_tail_balls"]
    normalized = c["weighted_normalized_derivative_tail_balls"]
    return f"""# Outer three-current derivative tail for the completed B trace

Date: 2026-08-13

Status: uniform analytic rational-tail derivative certificate; not a proof
of the complete completed-B estimate

Let `c=Bx/2`, `D=c^2-m^2`, and retain the first three paired rational
currents `C_0,C_1,C_2` from the sign-adapted B-endpoint Fresnel expansion.
The exact formulas are

```text
C_0=-(1+i*pi*B*c)/(pi^2 D),

C_1=c/(pi^3 B){{(-4*pi*B*c^3+4*i*c^2)D^-3
                  +(5*pi*B*c-3*i)D^-2}},

C_2=-i*c^2/(pi^4 B^2){{(-48*pi*B*c^5+48*i*c^4)D^-5
                         +(84*pi*B*c^3-60*i*c^2)D^-4
                         +(-35*pi*B*c+15*i)D^-3}}.       (OT1)
```

For the analytic outer tail choose `m>=B={B}`.  On the entire extracted
domain `delta={DELTA_TEXT}<=x<=1/2`,

```text
0<=c<=B/4,       c/m<=1/4,
|D|=m^2-c^2 >= (15/16)m^2.                              (OT2)
```

Each term in (OT1) is a monomial `alpha B^q c^p D^-k`.  It is differentiated
symbolically through order two before absolute values.  Every resulting
lattice tail is then enclosed by

```text
sum_(m=M)^infinity m^(-2k)
 <= M^(-2k)+M^(1-2k)/(2k-1).                            (OT3)
```

Consequently every finite cutoff and every positive Abel weight bounded by
one obeys the following uniform bounds for
`S=sum_(m>=B)(C_0+C_1+C_2)`:

```text
|S|    <= {raw['0']}
|S_x|  <= {raw['1']}
|S_xx| <= {raw['2']}                                    (OT4)
```

After multiplication by `w=[x(1-x)]^(-1/4)` and conversion to
`xi=sqrt(H_B)(x-x_B)`, the corresponding bounds are

```text
|w S|          <= {normalized['0']}
|d_xi(w S)|    <= {normalized['1']}
|d_xixi(w S)|  <= {normalized['2']}.                    (OT5)
```

The large unnormalized derivatives in (OT4) are expected near the artificial
corner `x=delta`; (OT5) records the actual normalized scale needed by the
tangential argument.  No mode enumeration is used.

What this closes: the first three *rational* boundary currents have a
cutoff-uniform, Abel-uniform analytic outer tail through two derivatives.

What remains open: differentiated bounds for the exact Fresnel expansion
remainder, the finite block `m<B`, both cutoff/window edge terms, and the
common-regulator combination with `R_join`.  Thus this gate proves no complete
B estimate, `T_upper`, height-uniform theorem, `Lambda<=0`, PF-infinity, RH,
or prize-level conclusion.

Pi provenance: every `pi` in (OT1) is inherited from the equation-(9) Fresnel
phase and its endpoint integrations by parts.  No geometric or fitted value
of `pi` is introduced in this certificate.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(
        dependencies["derivative_transport"]["decision"]["exact_first_second_pair_derivatives_proved"] is True,
        "derivative transport drift",
    )
    require(
        dependencies["higher_currents"]["decision"]["second_and_third_face_current_formulas_derived_exactly"] is True,
        "higher-current formula drift",
    )
    require(
        dependencies["tangential_cutoff"]["decision"]["exterior_B_trace_uniformly_tangentially_nonstationary"] is True,
        "tangential exterior drift",
    )

    artifact = {
        "kind": STEM,
        "status": "uniform_completed_B_outer_three_rational_current_tail_bounded_through_two_normalized_derivatives",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(),
        "decision": {
            "outer_mode_denominator_margin_proved_on_delta_to_half": True,
            "first_three_rational_current_outer_tail_analytically_summed": True,
            "outer_rational_amplitude_first_second_derivative_bounds_uniform_in_cutoff_and_Abel_weight": True,
            "differentiated_Fresnel_remainder_tail_bounded": False,
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
        "next_obligation": "Differentiate the exact sign-adapted Fresnel remainder with enough retained currents to obtain a summable outer tail through two x derivatives. Then combine it with this rational certificate and replace the remaining m<B roster by a finite-Poisson or Dirichlet representation before taking norms.",
        "proof_boundary": "Uniform analytic tail bounds through two derivatives for C0+C1+C2 on m>=B and delta<=x<=1/2 only. No differentiated exact expansion remainder, finite m<B block, complete exterior-B or joined-remainder bound, complete T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified completed-B outer three-current derivative tail", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
