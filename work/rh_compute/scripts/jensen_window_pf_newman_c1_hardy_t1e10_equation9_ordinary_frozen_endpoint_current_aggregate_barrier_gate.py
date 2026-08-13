#!/usr/bin/env python3
"""Certify the source-wide obstruction to a frozen endpoint-current Morse main."""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_frozen_endpoint_current_aggregate_barrier_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "coverage_ledger": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_ordinary_mode_coverage_ledger_gate.json",
    "interface_barrier": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_full_saddle_interface_barrier_gate.json",
    "universal_morse": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_universal_logistic_morse_characteristic_fold_reduction_gate.json",
}

PRECISION = 90
T = 10_000_000_000
A = 159_577
B = 5_122_421
LOWER_START = 622
UPPER_END = 39_894
SEGMENTS = (
    ("B_collar", 622, 999),
    ("low_bulk", 1_000, 9_999),
    ("middle_bulk", 10_000, 29_999),
    ("upper_bulk", 30_000, 38_999),
    ("A_collar", 39_000, 39_694),
    ("fold_overlap", 39_695, 39_852),
    ("lower_transition", 39_853, 39_894),
)


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


def complex_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": value.real.str(PRECISION, more=True),
        "imag_ball": value.imag.str(PRECISION, more=True),
    }


def fresnel(z: arb, pi: arb, imaginary_unit: acb) -> acb:
    scale = (-imaginary_unit * pi / 4).exp() * (pi / 2).sqrt()
    return (imaginary_unit * pi / 4).exp() / arb(2).sqrt() * (scale * z).erf()


def frozen_ratio(mode: int, pi: arb, imaginary_unit: acb) -> acb:
    m = arb(mode)
    t = arb(T)
    x = 2 * pi * m**2 / (t + 2 * pi * m**2)
    q_a = (x / 2).sqrt() * (arb(A) - 2 * m / x)
    q_b = (x / 2).sqrt() * (arb(B) - 2 * m / x)
    i0 = fresnel(q_b, pi, imaginary_unit) - fresnel(q_a, pi, imaginary_unit)
    i1 = (
        (imaginary_unit * pi * q_b**2 / 2).exp()
        - (imaginary_unit * pi * q_a**2 / 2).exp()
    ) / (imaginary_unit * pi)
    return i0 / (1 + imaginary_unit) + x.sqrt() * i1 / (
        m * arb(2).sqrt() * (1 + imaginary_unit)
    )


def summarize_segment(
    name: str,
    start: int,
    end: int,
    theta_zero: arb,
    theta_exact: arb,
    pi: arb,
    imaginary_unit: acb,
) -> dict[str, Any]:
    frozen_minus_exact = acb(0)
    endpoint_only = acb(0)
    phase_only = acb(0)
    termwise_triangle = arb(0)
    for mode in range(start, end + 1):
        m = arb(mode)
        ratio = frozen_ratio(mode, pi, imaginary_unit)
        carrier_zero = (imaginary_unit * (theta_zero - arb(T) * m.log())).exp() / m.sqrt()
        carrier_exact = (imaginary_unit * (theta_exact - arb(T) * m.log())).exp() / m.sqrt()
        defect = carrier_zero * ratio - carrier_exact
        frozen_minus_exact += defect
        endpoint_only += carrier_zero * (ratio - 1)
        phase_only += carrier_zero - carrier_exact
        termwise_triangle += abs(defect)

    return {
        "name": name,
        "range": [start, end],
        "count": end - start + 1,
        "complex_frozen_minus_exact_ball": complex_record(frozen_minus_exact),
        "two_real_frozen_minus_exact_ball": (2 * frozen_minus_exact.real).str(PRECISION, more=True),
        "complex_endpoint_current_only_ball": complex_record(endpoint_only),
        "two_real_endpoint_current_only_ball": (2 * endpoint_only.real).str(PRECISION, more=True),
        "complex_phase_only_ball": complex_record(phase_only),
        "two_real_phase_only_ball": (2 * phase_only.real).str(PRECISION, more=True),
        "termwise_complex_triangle_ball": termwise_triangle.str(PRECISION, more=True),
    }


def combine(rows: list[dict[str, Any]]) -> dict[str, Any]:
    total = acb(0)
    endpoint = acb(0)
    phase = acb(0)
    triangle = arb(0)
    for row in rows:
        total += acb(
            arb(row["complex_frozen_minus_exact_ball"]["real_ball"]),
            arb(row["complex_frozen_minus_exact_ball"]["imag_ball"]),
        )
        endpoint += acb(
            arb(row["complex_endpoint_current_only_ball"]["real_ball"]),
            arb(row["complex_endpoint_current_only_ball"]["imag_ball"]),
        )
        phase += acb(
            arb(row["complex_phase_only_ball"]["real_ball"]),
            arb(row["complex_phase_only_ball"]["imag_ball"]),
        )
        triangle += arb(row["termwise_complex_triangle_ball"])
    physical_scale = arb(2).sqrt() / arb.pi()
    return {
        "range": [LOWER_START, UPPER_END],
        "count": UPPER_END - LOWER_START + 1,
        "complex_frozen_minus_exact_ball": complex_record(total),
        "two_real_frozen_minus_exact_ball": (2 * total.real).str(PRECISION, more=True),
        "complex_endpoint_current_only_ball": complex_record(endpoint),
        "two_real_endpoint_current_only_ball": (2 * endpoint.real).str(PRECISION, more=True),
        "complex_phase_only_ball": complex_record(phase),
        "two_real_phase_only_ball": (2 * phase.real).str(PRECISION, more=True),
        "complex_discrepancy_absolute_ball": abs(total).str(PRECISION, more=True),
        "equation9_physical_scale_absolute_ball": (physical_scale * abs(total)).str(PRECISION, more=True),
        "termwise_complex_triangle_ball": triangle.str(PRECISION, more=True),
        "signed_to_triangle_ratio_ball": (abs(total) / triangle).str(PRECISION, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    total = artifact["certificate"]["aggregate"]
    phase = artifact["certificate"]["phase"]
    return f"""# Frozen endpoint-current aggregate barrier

Date: 2026-08-13

Status: source-wide frozen-current obstruction certified; not a proof of the non-frozen Morse remainder

For every source target mode `622<=m<=39894`, let `R_m` be the exact
incomplete Fresnel plus endpoint-current ratio evaluated at the outer Morse
saddle.  The frozen endpoint-retained model is

```text
M_frozen=2 Re sum_m exp(i(theta_0-t log m)) R_m/sqrt(m),
T_upper =2 Re sum_m exp(i(theta-t log m))/sqrt(m).
```

Here `theta` is the exact Riemann--Siegel theta from `log Gamma`, while
`theta_0=t/2[log(t/(2*pi))-1]-pi/8`.  Rigorous Arb summation over all
{total['count']} modes gives

```text
M_frozen-T_upper={total['two_real_frozen_minus_exact_ball']},
sum_m (frozen-exact)={total['complex_frozen_minus_exact_ball']['real_ball']}
                     + i {total['complex_frozen_minus_exact_ball']['imag_ball']},
equation-(9) physical scale={total['equation9_physical_scale_absolute_ball']}.
```

The carrier-phase difference is only

```text
theta-theta_0={phase['theta_exact_minus_theta_zero_ball']},
two-real phase-only aggregate={total['two_real_phase_only_ball']}.
```

It is therefore far too small to explain the discrepancy.  Signed
aggregation does reduce the termwise triangle
`{total['termwise_complex_triangle_ball']}`, but does not cancel the frozen
defect.

This rejects one specific approximation: replacing the exact global Morse
amplitude by its saddle value, even after retaining the incomplete endpoint
current there.  It does **not** reject the endpoint-retained Morse method and
does not lower-bound the exact equation-(9) error.  The omitted integral

```text
integral exp(-i*t*s^2/4) [A_m(s)-A_m(0)] ds
```

can be the same size near the fold and may supply the missing correction.  It
must be derived and bounded before any comparison with `T_upper` is valid.

Pi provenance: every `pi` here comes from the equation-(9) Kummer/Fourier
phase, the Gaussian Morse normalization, or the Riemann--Siegel theta
normalization.  No geometric fit or free polygonal constant is inserted.

Proof boundary: a rigorous finite saved-height sum for the frozen leading
model only.  No bound on the non-frozen Morse-amplitude integral, complete
`T_upper` theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is
proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    ctx.dps = PRECISION
    ctx.threads = 1

    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    coverage = dependencies["coverage_ledger"]["mode_coverage"]
    require(coverage["source_classical_target"] == {"range": [622, 39894], "count": 39273}, "source target drift")

    pi = arb.pi()
    imaginary_unit = acb(0, 1)
    t = arb(T)
    theta_exact = acb(arb("0.25"), t / 2).lgamma().imag - t * pi.log() / 2
    theta_zero = t / 2 * ((t / (2 * pi)).log() - 1) - pi / 8
    theta_correction = theta_exact - theta_zero

    rows = [
        summarize_segment(name, start, end, theta_zero, theta_exact, pi, imaginary_unit)
        for name, start, end in SEGMENTS
    ]
    aggregate = combine(rows)
    total_two_real = arb(aggregate["two_real_frozen_minus_exact_ball"])
    phase_two_real = arb(aggregate["two_real_phase_only_ball"])
    physical = arb(aggregate["equation9_physical_scale_absolute_ball"])
    triangle = arb(aggregate["termwise_complex_triangle_ball"])

    require(sum(row["count"] for row in rows) == 39_273, "segment count drift")
    require(total_two_real < arb("-0.096"), "frozen-current aggregate barrier failed")
    require(abs(phase_two_real) < arb("0.00000001"), "phase-only correction unexpectedly large")
    require(theta_correction > arb("2.0e-12") and theta_correction < arb("2.2e-12"), "theta correction drift")
    require(physical > arb("0.028"), "physical-scale obstruction failed")
    require(triangle > arb("0.61"), "termwise triangle drift")

    artifact = {
        "kind": STEM,
        "status": "frozen_endpoint_retained_morse_main_rigorously_rejected_as_complete_source_comparison",
        "passed": True,
        "model": {
            "frozen": "M_frozen=2 Re sum_(m=622)^39894 exp(i(theta_0-t log m)) R_m/sqrt(m)",
            "target": "T_upper=2 Re sum_(m=622)^39894 exp(i(theta-t log m))/sqrt(m)",
            "ratio": "R_m=I0/(1+i)+sqrt(x_m) I1/[m*sqrt(2)(1+i)]",
            "omitted_term": "sum_m carrier_m integral_R exp(-i*t*s^2/4)[A_m(s)-A_m(0)]ds",
        },
        "certificate": {
            "height": T,
            "alpha_endpoints": [A, B],
            "precision_decimal_digits": PRECISION,
            "phase": {
                "theta_exact_ball": theta_exact.str(PRECISION, more=True),
                "theta_zero_ball": theta_zero.str(PRECISION, more=True),
                "theta_exact_minus_theta_zero_ball": theta_correction.str(PRECISION, more=True),
            },
            "segments": rows,
            "aggregate": aggregate,
        },
        "decision": {
            "frozen_endpoint_current_is_complete_ordinary_join": False,
            "endpoint_retained_nonfrozen_morse_method_rejected": False,
            "nonfrozen_amplitude_variation_must_be_retained": True,
            "complete_x_integrated_error_lower_bound_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)} for name, path in DEPENDENCIES.items()
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
        "next_obligation": "Write the exact global Morse amplitude A_m(s), retain its endpoint dependence q_A(s), q_B(s), exploit odd cancellation on the full real s-line, and certify the signed aggregate of its first nonzero even correction before attempting a remainder bound.",
        "proof_boundary": "Rigorous source-wide finite sum of the frozen endpoint-current leading model at t=10^10 only. No bound on the non-frozen Morse-amplitude integral, complete T_upper theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }

    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "certified frozen endpoint-current aggregate barrier: "
        f"two-real={total_two_real}, physical={physical}, modes={aggregate['count']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
