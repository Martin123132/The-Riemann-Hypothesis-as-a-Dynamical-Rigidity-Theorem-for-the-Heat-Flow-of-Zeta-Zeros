#!/usr/bin/env python3
"""Audit the unequal-truncation parameter map for the joined equation-(4) packet."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import time
from typing import Any


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_unequal_truncation_saddle_map_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCY = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_PW_unowned_cell_pole_safe_two_boundary_contour_join_gate.json"
)

HEIGHT = 10_000_000_000
ALPHA_MIN = 159_577
ALPHA_MAX = 5_122_421
LOWER_BOUNDARY = mp.mpf("621.5")
UPPER_BOUNDARY = mp.mpf("39936.5")


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
            process.cpu_affinity([process.cpu_affinity()[0]])
            return "below_normal_one_cpu"
        process.nice(10)
        process.cpu_affinity([process.cpu_affinity()[0]])
        return "nice_10_one_cpu"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def text(value: mp.mpf | mp.mpc, digits: int = 45) -> str:
    return mp.nstr(value, digits)


def phase_saddle_roots(t: mp.mpf, label: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    discriminant = label * label - 8 * t / mp.pi
    require(discriminant >= 0, "label is below the real-saddle transition")
    root = mp.sqrt(discriminant)
    return (label - root) / 4, (label + root) / 4


def unequal_map(t: mp.mpf, y: mp.mpf) -> dict[str, mp.mpf]:
    p = t / (2 * mp.pi)
    partner = p / y
    larger = max(y, partner)
    smaller = min(y, partner)
    lam = mp.sqrt(larger / smaller)
    return {
        "alpha": larger,
        "beta": smaller,
        "lambda": lam,
        "label": 2 * (y + partner),
        "product": y * partner,
    }


def completed_ratio(t: mp.mpf, lam: mp.mpf) -> mp.mpf:
    """Ratio of successive displayed Theorem 1.5 remainder scales at sigma=1/2."""
    return mp.power(lam + 1 / lam, 3) / mp.sqrt(t)


def completed_scale(t: mp.mpf, lam: mp.mpf, order: int) -> mp.mpf:
    return mp.power(lam + 1 / lam, 3 * order + mp.mpf("0.5")) / mp.power(
        t, mp.mpf(order) / 2 + mp.mpf("0.25")
    )


def intermediate_scale(t: mp.mpf, lam: mp.mpf, order: int) -> mp.mpf:
    """Displayed Theorem 2.1 error scale at sigma=1/2, without its constant."""
    return mp.sqrt(lam) * mp.power(t, -mp.mpf("0.25") - mp.mpf(order) / 6)


def first_half_integer_not_below(value: mp.mpf) -> mp.mpf:
    return mp.mpf(math.ceil(float(value - mp.mpf("0.5")))) + mp.mpf("0.5")


def rational_root_fixture(beta: mp.mpf, alpha: mp.mpf) -> dict[str, str]:
    t = 2 * mp.pi * alpha * beta
    label = 2 * (alpha + beta)
    lower, upper = phase_saddle_roots(t, label)
    mapped = unequal_map(t, beta)
    return {
        "t_over_pi": text(t / mp.pi),
        "label": text(label),
        "beta": text(beta),
        "alpha": text(alpha),
        "root_lower": text(lower),
        "root_upper": text(upper),
        "product_discrepancy": text(abs(mapped["product"] - t / (2 * mp.pi))),
        "label_discrepancy": text(abs(mapped["label"] - label)),
        "root_discrepancy": text(max(abs(lower - beta), abs(upper - alpha))),
    }


def main() -> int:
    started = time.perf_counter()
    priority = set_low_priority()
    require(priority == "below_normal_one_cpu", f"resource cap unavailable: {priority}")
    require(DEPENDENCY.is_file(), "missing pole-safe join dependency")
    require(CHECKER.is_file(), "missing independent checker")

    mp.mp.dps = 90
    t = mp.mpf(HEIGHT)
    p = t / (2 * mp.pi)
    sqrt_p = mp.sqrt(p)

    # The exact displayed-scale crossover solves
    # (lambda+lambda^-1)^3=sqrt(t), not merely lambda=t^(1/6).
    c = mp.power(t, mp.mpf(1) / 6)
    lambda_star = (c + mp.sqrt(c * c - 4)) / 2
    y_star = sqrt_p / lambda_star
    split = first_half_integer_not_below(y_star)
    prior_half_integer = split - 1

    require(split == mp.mpf("860.5"), "actual-height half-cell split drift")
    require(completed_ratio(t, unequal_map(t, split)["lambda"]) < 1, "central side does not contract")
    require(
        completed_ratio(t, unequal_map(t, prior_half_integer)["lambda"]) > 1,
        "preceding half-cell unexpectedly contracts",
    )

    endpoint_first = int(LOWER_BOUNDARY + mp.mpf("0.5"))
    endpoint_last = int(split - mp.mpf("0.5"))
    central_first = endpoint_last + 1
    central_last = int(UPPER_BOUNDARY - mp.mpf("0.5"))
    endpoint_count = endpoint_last - endpoint_first + 1
    central_count = central_last - central_first + 1
    require((endpoint_first, endpoint_last, endpoint_count) == (622, 860, 239), "endpoint roster drift")
    require((central_first, central_last, central_count) == (861, 39936, 39076), "central roster drift")
    require(endpoint_count + central_count == 39315, "owned-cell count drift")

    label_rows: list[dict[str, str | int]] = []
    for label in (ALPHA_MIN, ALPHA_MAX):
        lower, upper = phase_saddle_roots(t, mp.mpf(label))
        mapped = unequal_map(t, lower)
        label_rows.append(
            {
                "odd_label": label,
                "lower_root": text(lower),
                "upper_root": text(upper),
                "root_sum_discrepancy": text(abs(lower + upper - mp.mpf(label) / 2)),
                "root_product_discrepancy": text(abs(lower * upper - p)),
                "mapped_label_discrepancy": text(abs(mapped["label"] - label)),
                "lambda": text(mapped["lambda"]),
                "completed_successive_scale_ratio": text(completed_ratio(t, mapped["lambda"])),
            }
        )

    sample_points: list[dict[str, Any]] = []
    for name, y in (
        ("lower_owned_boundary", LOWER_BOUNDARY),
        ("last_noncontracting_half_boundary", prior_half_integer),
        ("first_contracting_half_boundary", split),
        ("balanced_saddle", sqrt_p),
        ("upper_owned_boundary", UPPER_BOUNDARY),
    ):
        mapped = unequal_map(t, y)
        sample_points.append(
            {
                "name": name,
                "y": text(y),
                "partner": text(p / y),
                "stationary_label": text(mapped["label"]),
                "lambda": text(mapped["lambda"]),
                "product_discrepancy": text(abs(mapped["product"] - p)),
                "completed_successive_scale_ratio": text(completed_ratio(t, mapped["lambda"])),
                "theorem_1_5_displayed_scales_without_constant": [
                    text(completed_scale(t, mapped["lambda"], order)) for order in range(4)
                ],
                "theorem_2_1_displayed_scales_without_constant": [
                    text(intermediate_scale(t, mapped["lambda"], order)) for order in range(4)
                ],
            }
        )

    fixtures = [
        rational_root_fixture(mp.mpf("1.5"), mp.mpf("3")),
        rational_root_fixture(mp.mpf("2.5"), mp.mpf("4")),
    ]
    for fixture in fixtures:
        require(mp.mpf(fixture["root_discrepancy"]) < mp.mpf("1e-80"), "altered root fixture failed")

    artifact: dict[str, Any] = {
        "kind": STEM,
        "status": "unequal_truncation_saddle_map_and_route_split_certified",
        "passed": True,
        "scope": {
            "height": HEIGHT,
            "s": "1/2+i*t",
            "joined_target": "S_W-C_G*sum_(m=622)^39936 m^(-s)-mathcal_A_A^nat",
        },
        "dependency": {"path": relative(DEPENDENCY), "sha256": file_hash(DEPENDENCY)},
        "exact_parameter_map": {
            "p": "t/(2*pi)",
            "phase_stationarity": "ell=2*y+t/(pi*y)=2*(y+p/y)",
            "unequal_truncations": "alpha=max(y,p/y), beta=min(y,p/y)",
            "product": "t=2*pi*alpha*beta",
            "lambda": "sqrt(alpha/beta)",
            "stationary_roots": "y_+/-=(ell+/-sqrt(ell^2-8*t/pi))/4",
            "label_rows": label_rows,
            "altered_rational_root_fixtures": fixtures,
        },
        "actual_height_route_split": {
            "sqrt_t_over_2pi": text(sqrt_p),
            "theorem_1_5_successive_displayed_scale_ratio": "(lambda+lambda^(-1))^3/sqrt(t)",
            "lambda_crossover": text(lambda_star),
            "y_crossover": text(y_star),
            "first_half_integer_not_below_crossover": text(split),
            "preceding_half_integer": text(prior_half_integer),
            "endpoint_cells": {"first": endpoint_first, "last": endpoint_last, "count": endpoint_count},
            "central_cells": {"first": central_first, "last": central_last, "count": central_count},
            "owned_cell_count": endpoint_count + central_count,
            "sample_points": sample_points,
        },
        "literature_audit": {
            "primary_source": "https://arxiv.org/abs/1811.01130",
            "theorem_1_5_object": "classical symmetric zeta remainder R(s;alpha,beta)",
            "theorem_2_1_object": "classical zeta remainder before completed symmetric re-expansion",
            "direct_identity_with_joined_heat_flow_packet_proved": False,
            "explicit_numerical_error_constants_supplied_by_paper": False,
            "theorem_1_5_scales_are_certified_bounds_here": False,
            "theorem_2_1_scales_are_certified_bounds_here": False,
            "usable_as_exact_saddle_and_mordell_blueprint": True,
        },
        "decision": {
            "exact_saddle_map_certified": True,
            "completed_series_core_starts_at_cell": central_first,
            "completed_series_endpoint_quarantine_cells": endpoint_count,
            "preferred_donor_formula": "O'Sullivan Theorem 2.1 / equation (2.42) before completed re-expansion",
            "reason": "its sigma=1/2 displayed remainder scale contracts by t^(-1/6) and its stated lambda>=1 constant is lambda-independent, but constants must still be made explicit",
            "next_action": "derive an exact common-contour identity from the fully reassembled joined packet to a finite Mordell main plus an explicit remainder; retain Gamma and A subtraction before any norm",
            "J_Z_enclosed": False,
            "D_K_enclosed": False,
            "rh_implication": False,
        },
        "runtime": {
            "seconds": time.perf_counter() - started,
            "active_compute_workers": 1,
            "priority": priority,
            "thread_caps": 1,
        },
        "source_hashes": {"builder": file_hash(BUILDER), "checker": file_hash(CHECKER)},
        "proof_boundary": (
            "Exact stationary-phase/unequal-truncation parameter map, actual-height displayed-scale crossover, "
            "and a cancellation-preserving route split only. No identity between O'Sullivan's classical zeta "
            "remainder and J_Z, no explicit remainder constant, no numerical enclosure of the joined packet, "
            "J_Z, or D_K, and no non-A, all-height, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }

    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    note = f"""# Unequal-truncation saddle map for the joined equation-(4) packet

Date: 2026-08-27

Status: exact parameter-map and route-split certificate; quantitative joined-packet enclosure remains open.

## Exact map

Put `p=t/(2*pi)`.  The real phase of one Fresnel label is stationary when

```text
ell=2*y+t/(pi*y)=2*(y+p/y).
```

Therefore `alpha=max(y,p/y)` and `beta=min(y,p/y)` satisfy

```text
t=2*pi*alpha*beta,
lambda=sqrt(alpha/beta),
ell=2*(alpha+beta).
```

This is exactly the unequal-truncation hyperbola used by O'Sullivan, not a fitted analogy.

## Actual-height split

For Theorem 1.5, the ratio of successive displayed remainder scales at `sigma=1/2` is

```text
rho_5=(lambda+lambda^(-1))^3/sqrt(t).
```

At `t=10^10`, `rho_5=1` occurs at

```text
lambda={text(lambda_star)},
y={text(y_star)}.
```

The first half-integer boundary not below this crossover is `860.5`; the previous one, `859.5`, still has ratio above one.  Hence a completed-series implementation has the exact cell split

```text
endpoint quarantine: 622..860   (239 cells),
contracting core:     861..39936 (39076 cells).
```

This split is a route audit, not an error bound.  The paper's constants are implicit.

## Route decision

The better donor is Theorem 2.1 / equation (2.42), before the completed symmetric re-expansion.  Its displayed `sigma=1/2` error scale is `lambda^(1/2)t^(-1/4-N/6)` and the stated implied constant is independent of `lambda` for `lambda>=1`.  We must still derive an identity for the fully reassembled heat-flow packet and make every constant explicit.  O'Sullivan's `R(s;alpha,beta)` is not identified here with `J_Z`.

Primary source: https://arxiv.org/abs/1811.01130

Pi provenance: every `pi` above comes from the fixed quadratic Riemann-Siegel/Fresnel phase and the standard relation `t=2*pi*alpha*beta`; none is fitted or inserted geometrically.

## Proof boundary

{artifact['proof_boundary']}
"""
    NOTE.write_text(note, encoding="utf-8")
    print("certified unequal-truncation saddle map and actual-height route split")
    print(f"result: {relative(RESULT)}")
    print(f"note: {relative(NOTE)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
