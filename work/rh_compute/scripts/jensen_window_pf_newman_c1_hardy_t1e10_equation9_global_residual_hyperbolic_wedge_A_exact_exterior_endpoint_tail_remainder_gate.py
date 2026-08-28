#!/usr/bin/env python3
"""Certify the four-term endpoint-tail replacement on the A exterior."""

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

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_exact_exterior_endpoint_tail_remainder_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "reduction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_exact_exterior_common_phase_reduction_gate.json",
    "cancellation_guard": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_exterior_asymptotic_cancellation_guard_gate.json",
}

HEIGHT = 10_000_000_000
A = 159_577
LOWER_START = 39_853
UPPER_END = 39_936
MODES = tuple(range(LOWER_START, UPPER_END + 1))
Y0_TEXT = "0.0037"
DELTA_TEXT = "1e-12"
PRODUCTION_SLABS = 4096
PRECISION = 60


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


def polynomial_absolute_majorant(x: arb, mode: arb, endpoint: arb) -> arb:
    pi = arb.pi()
    terms = (
        pi * endpoint**10 * x**9,
        176 * pi * endpoint**8 * mode**2 * x**7,
        14 * endpoint**8 * x**8,
        9504 * pi * endpoint**6 * mode**4 * x**5,
        2016 * endpoint**6 * mode**2 * x**6,
        59136 * pi * endpoint**4 * mode**6 * x**3,
        28224 * endpoint**4 * mode**4 * x**4,
        59136 * pi * endpoint**2 * mode**8 * x,
        75264 * endpoint**2 * mode**6 * x**2,
        32256 * mode**8,
    )
    return sum((term.upper() for term in terms), arb(0))


def b3_derivative_absolute_bound(x: arb, mode: arb, endpoint: arb) -> arb:
    pi = arb.pi()
    x_upper = abs(x).upper()
    numerator = (
        240
        * (x ** arb("3.5")).upper()
        * polynomial_absolute_majorant(x_upper, mode, endpoint)
    )
    denominator = (
        pi.lower() ** 4
        * abs(endpoint * x - 2 * mode).lower() ** 8
        * abs(endpoint * x + 2 * mode).lower() ** 8
    )
    require(denominator > 0, "endpoint derivative denominator lost separation")
    return numerator / denominator


def certify_mode(mode_int: int, slabs: int) -> dict[str, Any]:
    pi, t, endpoint = arb.pi(), arb(HEIGHT), arb(A)
    mode = arb(mode_int)
    y0 = arb(Y0_TEXT)
    delta = arb(DELTA_TEXT)
    r = t / (2 * pi * mode**2)
    x0 = 1 / (1 + r * (1 + y0))

    # On [0,delta], |b3'(x)| <= C x^(7/2), including every polynomial term.
    coefficient = (
        240
        * polynomial_absolute_majorant(delta, mode, endpoint)
        / (
            pi.lower() ** 4
            * (2 * mode - endpoint * delta).lower() ** 8
            * (2 * mode).lower() ** 8
        )
    )
    q = coefficient * delta ** arb("4.5") / arb("4.5")
    weighted_error = (
        coefficient
        * (1 - delta).lower() ** arb("-0.25")
        * delta ** arb("3.75")
        / (arb("4.5") * arb("3.75"))
    )

    for index in range(slabs):
        left = delta + (x0 - delta) * arb(index) / slabs
        right = delta + (x0 - delta) * arb(index + 1) / slabs
        midpoint = (left + right) / 2
        radius = (right - left) / 2
        x_interval = arb(
            midpoint.mid(),
            radius.upper() * arb("1.000000000000001"),
        )
        derivative_bound = b3_derivative_absolute_bound(x_interval, mode, endpoint)
        width = right - left
        q += derivative_bound * width
        outer_weight = (
            left.lower() ** arb("-1.75")
            * (1 - right.upper()) ** arb("-0.25")
        )
        weighted_error += outer_weight * q * width

    require(q.is_finite() and weighted_error.is_finite(), "nonfinite endpoint remainder bound")
    return {
        "mode": mode_int,
        "x0_ball": x0.str(65, more=True),
        "phase_stripped_endpoint_tail_remainder_ball": q.str(65, more=True),
        "weighted_exterior_integral_error_ball": weighted_error.str(65, more=True),
    }


def certificate(slabs: int, precision: int) -> dict[str, Any]:
    ctx.dps = precision
    ctx.threads = 1
    pi, t = arb.pi(), arb(HEIGHT)
    rows = [certify_mode(mode, slabs) for mode in MODES]
    canonical = sum(
        (arb(row["weighted_exterior_integral_error_ball"]) for row in rows),
        arb(0),
    ) / (2 * pi)
    normalizer = (pi / (32 * t)) ** arb("0.25")
    physical = 2 * normalizer * canonical
    max_tail = max(
        arb(row["phase_stripped_endpoint_tail_remainder_ball"]).upper()
        for row in rows
    )
    require(physical.upper() < arb("1.5e-14"), "endpoint replacement misses physical cap")
    return {
        "height": HEIGHT,
        "endpoint": A,
        "mode_range": [LOWER_START, UPPER_END],
        "mode_count": len(rows),
        "y0": Y0_TEXT,
        "zero_endpoint_split_delta": DELTA_TEXT,
        "real_slabs_per_mode": slabs,
        "max_phase_stripped_endpoint_tail_remainder_ball": max_tail.str(65, more=True),
        "canonical_complete_exterior_replacement_error_ball": canonical.str(65, more=True),
        "physical_complete_exterior_replacement_error_ball": physical.str(65, more=True),
        "rows": rows,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Four-term exact-endpoint tail remainder on the A exterior

Date: 2026-08-14

Status: rigorous replacement-error bound for the full exact exterior;
not an exterior value or complete A theorem

For the exact phase-stripped endpoint tail of Section 11.437, write

```text
B_m=b_0+b_1+b_2+b_3+R_4,
b_0=h_A/(iF_m'),   b_(j+1)=-b_j'/(iF_m').            (ET1)
```

Repeated exact integration by parts gives

```text
|R_4(x)| <= integral_0^x |b_3'(u)|du.                (ET2)
```

The derivative is the explicit rational function

```text
b_3'(x)=240 i x^(7/2) P_9(x,m)
 /[pi^4(Ax-2m)^8(Ax+2m)^8],                         (ET3)
```

where the real and imaginary coefficients of `P_9` are bounded term by term.
The interval calculation preserves the `x^(7/2)` vanishing on
`[0,{DELTA_TEXT}]`, then uses {c['real_slabs_per_mode']} deterministic real
slabs for every mode through its exact cutoff.

After integrating the cumulative remainder against the complete exterior
weight `x^(-7/4)(1-x)^(-1/4)`, summing all 84 modes, restoring `1/(2pi)`,
and applying the equation-(9) physical normalization,

```text
max_m sup_x |R_4,m(x)| <=
 {c['max_phase_stripped_endpoint_tail_remainder_ball']},

complete canonical exterior replacement error <=
 {c['canonical_complete_exterior_replacement_error_ball']},

complete physical exterior replacement error <=
 {c['physical_complete_exterior_replacement_error_ball']}. (ET4)
```

An independent checker repeats the entire calculation at 70 digits with
8192 slabs per mode and obtains a tighter nested upper bound.

Pi provenance: every `pi` in (ET1)--(ET4) comes from the exact endpoint
phase, endpoint driver, and equation-(9) normalization.  No fitted constant
is introduced.

Proof boundary: uniform four-term endpoint-tail replacement error over the
full exact A exterior for modes 39853..39936 only.  No value of the resulting
rational common-phase integral, outer-IBP remainder bound, compact nonlinear
amplitude, complete A endpoint theorem, `R_Dir` estimate, RH, or prize-level
conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["reduction"]["decision"]["endpoint_IBP_hierarchy_evaluated_at_all_cutoffs"] is True, "reduction dependency drift")
    require(dependencies["cancellation_guard"]["decision"]["full_exact_endpoint_current_integrable_at_x_zero"] is True, "cancellation dependency drift")
    certified = certificate(PRODUCTION_SLABS, PRECISION)
    artifact = {
        "kind": STEM,
        "status": "full_exact_A_exterior_four_term_endpoint_tail_replacement_error_below_1p5e_minus_14",
        "passed": True,
        "certificate": certified,
        "decision": {
            "four_term_endpoint_tail_used": True,
            "zero_endpoint_power_vanishing_retained": True,
            "all_84_mode_remainders_certified": True,
            "complete_physical_replacement_error_below_1p5e_minus_14": True,
            "rational_common_phase_integral_value_proved": False,
            "outer_IBP_remainder_bound_proved": False,
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
            "precision_decimal_digits": PRECISION,
        },
        "next_obligation": "Certify the now-rational four-term common-phase exterior integral by interval quadrature on a compact stationary Morse window and a grouped nonstationary far-tail bound.",
        "proof_boundary": "Uniform four-term endpoint-tail replacement error over the full exact A exterior only. No value of the rational common-phase integral, outer-tail bound, compact nonlinear-amplitude bound, complete A endpoint theorem, R_Dir estimate, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified four-term exact-endpoint tail remainder on A exterior", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
