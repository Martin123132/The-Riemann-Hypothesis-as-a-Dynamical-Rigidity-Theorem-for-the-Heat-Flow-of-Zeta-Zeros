#!/usr/bin/env python3
"""Reconcile direct Fresnel and z-boundary three-term B expansions."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_fresnel_boundary_truncation_dictionary_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "B_window": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_gaussian_window_fresnel_remainder_gate.json",
    "higher_currents": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_higher_boundary_current_absolute_window_gate.json",
}

PRECISION = 100
T = 10_000_000_000
B = 5_122_421
CUTOFF = 100_000


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
    x, endpoint, mode = sp.symbols("x B m", positive=True, real=True)
    pi, imaginary_unit = sp.pi, sp.I
    c = endpoint * x / 2
    denominator = c**2 - mode**2

    positive = (
        endpoint / (imaginary_unit * pi * (endpoint * x - 2 * mode))
        - mode / (2 * pi**2 * (c - mode) ** 3)
        + 3 * imaginary_unit * mode * x / (4 * pi**3 * (c - mode) ** 5)
    )
    negative = (
        endpoint / (imaginary_unit * pi * (endpoint * x + 2 * mode))
        + mode / (2 * pi**2 * (c + mode) ** 3)
        - 3 * imaginary_unit * mode * x / (4 * pi**3 * (c + mode) ** 5)
    )
    current0 = -2 * (2 + imaginary_unit * pi * endpoint**2 * x) / (
        pi**2 * (endpoint**2 * x**2 - 4 * mode**2)
    )
    current1 = c / (pi**3 * endpoint) * (
        (-4 * pi * endpoint * c**3 + 4 * imaginary_unit * c**2) / denominator**3
        + (5 * pi * endpoint * c - 3 * imaginary_unit) / denominator**2
    )
    current2 = -imaginary_unit * c**2 / (pi**4 * endpoint**2) * (
        (-48 * pi * endpoint * c**5 + 48 * imaginary_unit * c**4) / denominator**5
        + (84 * pi * endpoint * c**3 - 60 * imaginary_unit * c**2) / denominator**4
        + (-35 * pi * endpoint * c + 15 * imaginary_unit) / denominator**3
    )
    correction = 48 * x**2 * (
        endpoint**4 * x**4 + 40 * endpoint**2 * mode**2 * x**2 + 80 * mode**4
    ) / (pi**4 * (endpoint * x - 2 * mode) ** 5 * (endpoint * x + 2 * mode) ** 5)
    require(
        sp.factor(sp.together(current0 + current1 + current2 - positive - negative - correction)) == 0,
        "Fresnel/boundary truncation dictionary failed",
    )
    return {
        "direct_Fresnel_terms": "F_3,m=P_0,+ + P_1,+ + P_2,+ + P_0,- + P_1,- + P_2,-",
        "boundary_currents": "C_3,m=C_0,m+C_1,m+C_2,m",
        "dictionary": "C_3,m=F_3,m+D_m",
        "correction": "D_m=48*x^2*(B^4*x^4+40*B^2*m^2*x^2+80*m^4)/(pi^4*(B*x-2m)^5*(B*x+2m)^5)",
        "remainder_transfer": "exact-C_3=(exact-F_3)-D_m, hence |exact-C_3|<=|R_F,m|+|D_m|",
    }


def interval_certificate(window: dict[str, Any]) -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    pi, endpoint = arb.pi(), arb(B)
    x_low = arb(window["x_window_low_ball"])
    x_high = arb(window["x_window_high_ball"])
    require(0 < x_low < x_high < arb(1) / 2, "window geometry drift")

    numerator_common = endpoint**4 * x_high**4
    correction_sum = arb(0)
    for mode in range(1, 621):
        m = arb(mode)
        numerator = numerator_common + 40 * endpoint**2 * m**2 * x_high**2 + 80 * m**4
        left = endpoint * x_low - 2 * m
        right = endpoint * x_low + 2 * m
        require(left > 0, "left-side denominator lost positivity")
        correction_sum += 48 * x_high**2 * numerator / (pi**4 * left**5 * right**5)
    for mode in range(623, CUTOFF + 1):
        m = arb(mode)
        numerator = numerator_common + 40 * endpoint**2 * m**2 * x_high**2 + 80 * m**4
        left = 2 * m - endpoint * x_high
        right = 2 * m + endpoint * x_low
        require(left > 0, "right-side denominator lost positivity")
        correction_sum += 48 * x_high**2 * numerator / (pi**4 * left**5 * right**5)

    first_tail_mode = arb(CUTOFF + 1)
    rho = 1 - endpoint * x_high / (2 * first_tail_mode)
    ratio = endpoint * x_high / first_tail_mode
    polynomial = ratio**4 + 40 * ratio**2 + 80
    tail = (
        48
        * x_high**2
        * polynomial
        / (pi**4 * arb(2) ** 10 * rho**5)
        * arb(CUTOFF) ** (-5)
        / 5
    )
    correction_sum += tail

    width = x_high - x_low
    weight_max = (x_low * (1 - x_low)) ** (-arb(1) / 4)
    physical_factor = 2 * (pi / (32 * arb(T))) ** (arb(1) / 4)
    physical = physical_factor * width * weight_max * correction_sum
    require(physical < arb("2.1e-20"), "dictionary correction exceeds 2.1e-20")
    return {
        "height": T,
        "endpoint": B,
        "window_xi": [-70, 70],
        "sum_cutoff": CUTOFF,
        "uniform_nonlocal_dictionary_correction_ball": correction_sum.str(PRECISION, more=True),
        "analytic_tail_ball": tail.str(PRECISION, more=True),
        "physical_dictionary_correction_ball": physical.str(PRECISION, more=True),
        "physical_dictionary_correction_upper_bound": "2.1e-20",
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    return f"""# Fresnel-to-boundary B-face truncation dictionary

Date: 2026-08-13

Status: exact algebra plus saved-height interval correction; not a complete
B-face estimate

The three direct sign-adapted Fresnel terms and the first three repeated
`z`-boundary currents are two nearby, but not identical, truncations.  For
one paired positive index `m`, exact rational algebra gives

```text
C_0,m+C_1,m+C_2,m = F_0:2,m + D_m,                  (TD1)

D_m=48*x^2*(B^4*x^4+40*B^2*m^2*x^2+80*m^4)
    /[pi^4*(B*x-2m)^5*(B*x+2m)^5].                  (TD2)
```

Therefore a direct Fresnel remainder cannot be attached to the boundary
truncation as though `D_m=0`; instead

```text
|exact-(C_0+C_1+C_2)| <= |R_F,m|+|D_m|.             (TD3)
```

Modes 621 and 622 are removed before summation.  Endpoint-monotone
denominators, a finite Arb sum through `{c['sum_cutoff']}`, and an analytic
`m^-6` tail prove uniformly on `|xi|<=70`

```text
sum_(m notin {{621,622}})|D_m|
 <= {c['uniform_nonlocal_dictionary_correction_ball']}.             (TD4)
```

After the exact equation-(9) projection factor, window width, and maximum
Kummer weight are restored, this costs only

```text
{c['physical_dictionary_correction_ball']} < 2.1e-20.               (TD5)
```

TD5 is numerically negligible but logically required.  It repairs the
dictionary between the direct Fresnel remainder and the boundary-current
budget without changing the existing `1.491e-6` higher-current ceiling.

Pi provenance: every power of `pi` in TD1--TD5 comes from the exact Fresnel
phase, its integration-by-parts recurrence, or the equation-(9) physical
normalization.  No fitted constant is introduced.

Proof boundary: the nonlocal truncation-dictionary correction on the one
saved-height B window only.  No local 621/622 replacement, outside-window
tail, complete B estimate, A-fold splice, complete `T_upper`, height-uniform
theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["B_window"]["decision"]["nonlocal_completed_remainder_absolutely_summable"] is True, "Fresnel dependency drift")
    require(dependencies["higher_currents"]["decision"]["second_and_third_face_current_formulas_derived_exactly"] is True, "boundary-current dependency drift")
    artifact = {
        "kind": STEM,
        "status": "exact_Fresnel_boundary_truncation_dictionary_with_physical_correction_below_2p1e_minus_20",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(dependencies["B_window"]["interval_certificate"]),
        "decision": {
            "direct_Fresnel_and_boundary_current_truncations_identical": False,
            "exact_rational_dictionary_correction_derived": True,
            "nonlocal_dictionary_correction_absolutely_summed": True,
            "physical_dictionary_correction_below_2p1e_minus_20": True,
            "complete_B_face_estimate_proved": False,
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
        "next_obligation": "Insert this correction into the B-window partial budget, then use the direct Fresnel truncation on normal-safe local arcs and the exact continuous step-tail replacement on complementary outer-safe arcs.",
        "proof_boundary": "Only the exact three-term truncation dictionary and its nonlocal saved-height |xi|<=70 correction. No local replacement, outside-window tail, complete B estimate, A-fold splice, complete T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level theorem is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified Fresnel/boundary truncation dictionary correction below 2.1e-20", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
