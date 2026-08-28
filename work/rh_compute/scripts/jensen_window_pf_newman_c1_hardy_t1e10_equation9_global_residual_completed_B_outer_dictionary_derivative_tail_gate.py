#!/usr/bin/env python3
"""Bound two derivatives of the outer Fresnel-to-boundary dictionary tail."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_outer_dictionary_derivative_tail_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "outer_remainder": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_outer_fresnel_remainder_derivative_tail_gate.json",
    "dictionary": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_fresnel_boundary_truncation_dictionary_gate.json",
}

PRECISION = 100
T = 10_000_000_000
B = 5_122_421
M = B
DELTA_TEXT = "1e-4"

# (integer coefficient, c power, m power, denominator power) after the
# common factor 3/(pi^4 B^2).
DICTIONARY_TERMS = ((1, 6, 0, 5), (10, 4, 2, 5), (5, 2, 4, 5))


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


def derivative_components(order: int, p: int, k: int) -> tuple[tuple[int, int, int], ...]:
    if order == 0:
        return ((1, p, k),)
    if order == 1:
        return tuple(row for row in ((p, p - 1, k), (2 * k, p + 1, k + 1)) if row[0])
    require(order == 2, "derivative order outside zero through two")
    return tuple(
        row
        for row in (
            (p * (p - 1), p - 2, k),
            (2 * k * (2 * p + 1), p, k + 1),
            (4 * k * (k + 1), p + 2, k + 2),
        )
        if row[0]
    )


def series_tail(exponent: int) -> arb:
    require(exponent > 1, "dictionary lattice tail is not summable")
    start = arb(M)
    return start ** (-exponent) + start ** (1 - exponent) / (exponent - 1)


def symbolic_certificate() -> dict[str, str]:
    x, endpoint, mode, c = sp.symbols("x B m c", positive=True, real=True)
    raw = 48 * x**2 * (
        endpoint**4 * x**4 + 40 * endpoint**2 * mode**2 * x**2 + 80 * mode**4
    ) / (sp.pi**4 * (endpoint * x - 2 * mode) ** 5 * (endpoint * x + 2 * mode) ** 5)
    reduced = sp.factor(raw.subs(x, 2 * c / endpoint))
    expected = 3 * c**2 * (c**4 + 10 * mode**2 * c**2 + 5 * mode**4) / (
        sp.pi**4 * endpoint**2 * (c**2 - mode**2) ** 5
    )
    require(sp.simplify(reduced - expected) == 0, "dictionary c-reduction drift")
    for order in range(3):
        require(sp.diff(expected, c, order).has(mode), "dictionary derivative lost mode dependence")
    return {
        "original": "D_m=48x^2(B^4x^4+40B^2m^2x^2+80m^4)/[pi^4(Bx-2m)^5(Bx+2m)^5]",
        "reduced": "D_m=3c^2(c^4+10m^2c^2+5m^4)/[pi^4 B^2(c^2-m^2)^5], c=Bx/2",
        "derivative_rule": "Differentiate each c^p(c^2-m^2)^(-k) monomial before absolute values; d/dx=(B/2)d/dc",
        "outer_margin": "m>=B and c<=B/4 imply |c^2-m^2|>=(15/16)m^2",
    }


def interval_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    endpoint, pi = arb(B), arb.pi()
    cmax, rho = endpoint / 4, arb(15) / 16
    common = arb(3) / (pi**4 * endpoint**2)
    raw = {order: arb(0) for order in range(3)}
    for coefficient, p, mode_power, k in DICTIONARY_TERMS:
        for order in range(3):
            for factor, cp, kp in derivative_components(order, p, k):
                exponent = 2 * kp - mode_power
                raw[order] += (
                    common
                    * coefficient
                    * factor
                    * cmax**cp
                    * rho**(-kp)
                    * series_tail(exponent)
                    * (endpoint / 2) ** order
                )

    delta, t = arb(DELTA_TEXT), arb(T)
    x0 = (1 - (1 - 8 * t / (pi * endpoint**2)).sqrt()) / 2
    hessian = t * (1 - 2 * x0) / (2 * x0**2 * (1 - x0) ** 2)
    weight = (delta * (1 - delta)) ** (-arb(1) / 4)
    l1 = (1 - 2 * delta) / (4 * delta * (1 - delta))
    w2 = (12 * delta**2 - 12 * delta + 5) / (16 * delta**2 * (1 - delta) ** 2)
    weighted = {
        0: weight * raw[0],
        1: weight * (raw[1] + l1 * raw[0]),
        2: weight * (raw[2] + 2 * l1 * raw[1] + w2 * raw[0]),
    }
    normalized = {0: weighted[0], 1: weighted[1] / hessian.sqrt(), 2: weighted[2] / hessian}
    require(all(value.is_finite() and value > 0 for value in normalized.values()), "nonfinite dictionary bound")
    encode = lambda values: {str(order): value.str(PRECISION, more=True) for order, value in values.items()}
    return {
        "height": T,
        "endpoint": B,
        "outer_mode_start": M,
        "x_domain": [DELTA_TEXT, "1/2"],
        "dictionary_x_derivative_tail_balls": encode(raw),
        "weighted_x_derivative_tail_balls": encode(weighted),
        "weighted_normalized_derivative_tail_balls": encode(normalized),
        "derivative_order_legend": {"0": "dictionary amplitude", "1": "first derivative", "2": "second derivative"},
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    raw = c["dictionary_x_derivative_tail_balls"]
    norm = c["weighted_normalized_derivative_tail_balls"]
    return f"""# Outer Fresnel-to-boundary dictionary derivative tail

Date: 2026-08-13

Status: uniform analytic dictionary-tail derivative certificate; not a proof
of the complete completed-B estimate

The direct three-term Fresnel expansion and the first three endpoint
boundary currents differ exactly by

```text
D_m=48x^2(B^4x^4+40B^2m^2x^2+80m^4)
    /[pi^4(Bx-2m)^5(Bx+2m)^5].                       (OD1)
```

With `c=Bx/2`, this is

```text
D_m=3c^2(c^4+10m^2c^2+5m^4)
    /[pi^4 B^2(c^2-m^2)^5].                          (OD2)
```

Every monomial in (OD2) is differentiated through order two before norms.
For `m>=B={B}` and `{DELTA_TEXT}<=x<=1/2`, the same
`|c^2-m^2|>=(15/16)m^2` margin and decreasing-series integral test give

```text
sum_(m>=B)|D_m|    <= {raw['0']}
sum_(m>=B)|D_m'|   <= {raw['1']}
sum_(m>=B)|D_m''|  <= {raw['2']}.                     (OD3)
```

After the Kummer weight and normalized coordinate,

```text
|w D|          <= {norm['0']}
|d_xi(w D)|    <= {norm['1']}
|d_xixi(w D)|  <= {norm['2']}.                        (OD4)
```

Together with the rational-current and exact Fresnel-remainder gates, this
closes the complete analytic outer block `m>=B` through two normalized
derivatives.  The finite block `m<B`, phase-coupled exterior integration,
all indicator/cutoff boundary terms, and the common-regulator joined
remainder remain open.  No complete B estimate, `T_upper`, height-uniform
theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.

Pi provenance: `pi` is inherited from the exact equation-(9) Fresnel phase
and the endpoint integration-by-parts dictionary; it is not inserted as a
geometric or fitted parameter.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(
        dependencies["outer_remainder"]["decision"]["complete_outer_mode_derivative_block_bounded"] is False,
        "outer-remainder proof boundary drift",
    )
    require(
        dependencies["dictionary"]["decision"]["exact_rational_dictionary_correction_derived"] is True,
        "dictionary dependency drift",
    )

    artifact = {
        "kind": STEM,
        "status": "uniform_completed_B_outer_dictionary_tail_bounded_through_two_normalized_derivatives",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(),
        "decision": {
            "exact_dictionary_reduced_to_outer_rational_monomials": True,
            "dictionary_outer_tail_analytically_summed_through_two_derivatives": True,
            "complete_outer_mode_derivative_block_bounded": True,
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
        "next_obligation": "Use x-dependent outer majorants together with the exact tangential phase G_B, rather than a global amplitude supremum, to bound the complete analytic m>=B exterior contribution including x_low, x_high, and x=1/2 boundary terms. Separately compress the finite m<B block before norms.",
        "proof_boundary": "Uniform analytic outer dictionary-tail bounds through two derivatives, completing the rational plus remainder plus dictionary outer derivative block on m>=B and delta<=x<=1/2. No finite m<B block, phase-coupled exterior integral, complete exterior-B or joined-remainder bound, complete T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified completed-B outer dictionary derivative tail", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
