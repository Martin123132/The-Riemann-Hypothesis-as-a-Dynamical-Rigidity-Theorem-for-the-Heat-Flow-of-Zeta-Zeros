#!/usr/bin/env python3
"""Validate a stable O(L) reference evaluator for the non-A theta current."""

from __future__ import annotations

import cmath
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time
from typing import Any


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
from mpmath import mp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_reference_evaluator_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCY = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_common_kernel_theta_current_gate.json"

A = 159_577
B = 5_122_421
L = 2_481_422
DEFAULT_BLOCK = 256
SHORT_DPS = 90
FULL_DPS = 100
NORMALIZED_ALLOWANCE = 2.0e-10

SHORT_CASES = (
    {"x": "0.000000000001", "s": "0", "length": 17},
    {"x": "0.00000001", "s": "1/3", "length": 257},
    {"x": "0.1732050807568877", "s": "49/100", "length": 2_048},
    {"x": "0.99999999", "s": "7/8", "length": 4_096},
)

FULL_PERIODIC_CASES = (
    {"x": "0", "s": "1/3"},
    {"x": "1", "s": "0"},
    {"x": "1/2", "s": "0"},
    {"x": "1/3", "s": "1/2"},
    {"x": "2/5", "s": "1/3"},
    {"x": "7/16", "s": "3/8"},
    {"x": "1/4096", "s": "1/2"},
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_fraction(text: str) -> Fraction:
    return Fraction(text)


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


def reduce_mod_two(value: Fraction) -> Fraction:
    return value % 2


def cis_pi_fraction(value: Fraction) -> complex:
    reduced = reduce_mod_two(value)
    angle = math.pi * (reduced.numerator / reduced.denominator)
    return complex(math.cos(angle), math.sin(angle))


def mp_cis_pi_fraction(value: Fraction) -> Any:
    reduced = reduce_mod_two(value)
    angle = mp.pi * mp.mpf(reduced.numerator) / reduced.denominator
    return mp.exp(mp.j * angle)


def complex_record(value: complex) -> dict[str, float]:
    return {"real": float(value.real), "imag": float(value.imag)}


def mp_complex_record(value: Any, digits: int = 75) -> dict[str, str]:
    return {"real": mp.nstr(mp.re(value), digits), "imag": mp.nstr(mp.im(value), digits)}


def absolute_mass(length: int, s: Fraction) -> dict[str, float]:
    s0 = float(length)
    s1 = float(length * (length - 1) // 2)
    current = float(Fraction(length) * (A + 2 * s + length - 1))
    return {"S0": max(1.0, s0), "S1": max(1.0, s1), "F": max(1.0, current)}


def normalized_errors(
    observed: tuple[complex, complex, complex],
    expected: tuple[Any, Any, Any],
    length: int,
    s: Fraction,
) -> dict[str, float]:
    masses = absolute_mass(length, s)
    labels = ("S0", "S1", "F")
    return {
        label: float(abs(mp.mpc(observed[index]) - expected[index])) / masses[label]
        for index, label in enumerate(labels)
    }


def recurrence_theta_current(
    x: Fraction,
    s: Fraction,
    length: int,
    block_size: int = DEFAULT_BLOCK,
) -> tuple[complex, complex, complex]:
    """Evaluate S0, S1 and F with exact-rational block phase reseeding."""

    require(length >= 1, "theta roster must be nonempty")
    require(block_size >= 1, "block size must be positive")
    if x == 0:
        s0 = complex(length, 0.0)
        s1 = complex(length * (length - 1) // 2, 0.0)
        current = complex(float(Fraction(length) * (A + 2 * s + length - 1)), 0.0)
        return s0, s1, current

    q = cis_pi_fraction(2 * x)
    block_s0_real: list[float] = []
    block_s0_imag: list[float] = []
    block_s1_real: list[float] = []
    block_s1_imag: list[float] = []

    for start in range(0, length, block_size):
        stop = min(length, start + block_size)
        z_phase = x * (start * start + start * (A + 2 * s))
        ratio_phase = x * (A + 2 * s + 2 * start + 1)
        z = cis_pi_fraction(z_phase)
        ratio = cis_pi_fraction(ratio_phase)
        local_s0_real: list[float] = []
        local_s0_imag: list[float] = []
        local_s1_real: list[float] = []
        local_s1_imag: list[float] = []

        for n in range(start, stop):
            local_s0_real.append(z.real)
            local_s0_imag.append(z.imag)
            local_s1_real.append(n * z.real)
            local_s1_imag.append(n * z.imag)
            z *= ratio
            ratio *= q

        block_s0_real.append(math.fsum(local_s0_real))
        block_s0_imag.append(math.fsum(local_s0_imag))
        block_s1_real.append(math.fsum(local_s1_real))
        block_s1_imag.append(math.fsum(local_s1_imag))

    s0 = complex(math.fsum(block_s0_real), math.fsum(block_s0_imag))
    s1 = complex(math.fsum(block_s1_real), math.fsum(block_s1_imag))
    base = cis_pi_fraction(x * (A + 2 * s) ** 2 / 4)
    current = base * (float(A + 2 * s) * s0 + 2.0 * s1)
    return s0, s1, current


def direct_high_precision(x: Fraction, s: Fraction, length: int) -> tuple[Any, Any, Any]:
    mp.dps = SHORT_DPS
    s0 = mp.mpc(0)
    s1 = mp.mpc(0)
    current = mp.mpc(0)
    for n in range(length):
        reduced_phase = x * (n * n + n * (A + 2 * s))
        z = mp_cis_pi_fraction(reduced_phase)
        s0 += z
        s1 += n * z
        alpha = A + 2 * n + 2 * s
        current += mp.mpf(alpha.numerator) / alpha.denominator * mp_cis_pi_fraction(x * alpha**2 / 4)
    return s0, s1, current


def is_even_integer(value: Fraction) -> bool:
    return value.denominator == 1 and value.numerator % 2 == 0


def minimal_period(x: Fraction, s: Fraction, search_limit: int = 131_072) -> int:
    """Find P with z_(n+P)=z_n from two exact congruences."""

    if x == 0:
        return 1
    for period in range(1, search_limit + 1):
        xp = x * period
        if xp.denominator != 1:
            continue
        if is_even_integer(xp * (period + A + 2 * s)):
            return period
    raise RuntimeError(f"no quadratic period found below {search_limit} for x={x}, s={s}")


def periodic_high_precision(
    x: Fraction,
    s: Fraction,
    length: int,
    period: int,
) -> tuple[Any, Any, Any]:
    """Collapse a full roster to one exact phase cycle and arithmetic moments."""

    mp.dps = FULL_DPS
    cycles, remainder = divmod(length, period)
    z_cycle: list[Any] = []
    w_cycle: list[Any] = []
    alpha_cycle: list[Any] = []
    for j in range(period):
        reduced_phase = x * (j * j + j * (A + 2 * s))
        z_cycle.append(mp_cis_pi_fraction(reduced_phase))
        alpha = A + 2 * j + 2 * s
        alpha_mp = mp.mpf(alpha.numerator) / alpha.denominator
        alpha_cycle.append(alpha_mp)
        w_cycle.append(mp_cis_pi_fraction(x * alpha**2 / 4))

    z0 = mp.fsum(z_cycle)
    z1 = mp.fsum(j * z_cycle[j] for j in range(period))
    zr0 = mp.fsum(z_cycle[:remainder])
    zr1 = mp.fsum(j * z_cycle[j] for j in range(remainder))
    s0 = cycles * z0 + zr0
    s1 = period * cycles * (cycles - 1) * z0 / 2 + cycles * z1
    s1 += cycles * period * zr0 + zr1

    w0 = mp.fsum(w_cycle)
    wa = mp.fsum(alpha_cycle[j] * w_cycle[j] for j in range(period))
    wr0 = mp.fsum(w_cycle[:remainder])
    wra = mp.fsum(alpha_cycle[j] * w_cycle[j] for j in range(remainder))
    current = cycles * wa + period * cycles * (cycles - 1) * w0
    current += 2 * cycles * period * wr0 + wra
    return s0, s1, current


def run_short_rows() -> list[dict[str, Any]]:
    rows = []
    for case in SHORT_CASES:
        x = parse_fraction(case["x"])
        s = parse_fraction(case["s"])
        length = int(case["length"])
        observed = recurrence_theta_current(x, s, length)
        expected = direct_high_precision(x, s, length)
        errors = normalized_errors(observed, expected, length, s)
        require(max(errors.values()) <= NORMALIZED_ALLOWANCE, f"short-roster mismatch: {case}")
        rows.append(
            {
                **case,
                "block_size": DEFAULT_BLOCK,
                "recurrence": {
                    "S0": complex_record(observed[0]),
                    "S1": complex_record(observed[1]),
                    "F": complex_record(observed[2]),
                },
                "direct_high_precision": {
                    "S0": mp_complex_record(expected[0]),
                    "S1": mp_complex_record(expected[1]),
                    "F": mp_complex_record(expected[2]),
                },
                "mass_normalized_errors": errors,
                "passed": True,
            }
        )
    return rows


def run_periodic_rows() -> list[dict[str, Any]]:
    rows = []
    for case in FULL_PERIODIC_CASES:
        x = parse_fraction(case["x"])
        s = parse_fraction(case["s"])
        period = minimal_period(x, s)
        observed = recurrence_theta_current(x, s, L)
        expected = periodic_high_precision(x, s, L, period)
        errors = normalized_errors(observed, expected, L, s)
        require(max(errors.values()) <= NORMALIZED_ALLOWANCE, f"full periodic mismatch: {case}")
        xp = x * period
        endpoint_phase = xp * (period + A + 2 * s)
        rows.append(
            {
                **case,
                "length": L,
                "minimal_period": period,
                "period_congruences": {
                    "x_times_period": str(xp),
                    "xP_times_P_plus_A_plus_2s": str(endpoint_phase),
                    "xP_is_integer": xp.denominator == 1,
                    "endpoint_phase_is_even_integer": is_even_integer(endpoint_phase),
                },
                "block_size": DEFAULT_BLOCK,
                "recurrence": {
                    "S0": complex_record(observed[0]),
                    "S1": complex_record(observed[1]),
                    "F": complex_record(observed[2]),
                },
                "periodic_high_precision": {
                    "S0": mp_complex_record(expected[0]),
                    "S1": mp_complex_record(expected[1]),
                    "F": mp_complex_record(expected[2]),
                },
                "mass_normalized_errors": errors,
                "passed": True,
            }
        )
    return rows


def render_note(artifact: dict[str, Any]) -> str:
    summary = artifact["summary"]
    return f"""# Non-A quadratic-theta reference evaluator

Date: 2026-08-24

Status: one-worker floating reference evaluator validated against independent
high-precision and exact-periodic oracles; interval error theorem and physical
quadrature open

The exact non-A current from Section 11.460 is

```text
S_0(x,s)=sum_(n=0)^(L-1) z_n,
S_1(x,s)=sum_(n=0)^(L-1) n z_n,
z_n=exp(i*pi*x[n^2+n(A+2s)]),

F_x(s)=exp(i*pi*x(A+2s)^2/4)[(A+2s)S_0+2S_1],
L={L}.                                               (RE1)
```

The reference evaluator uses

```text
z_(n+1)=z_n r_n,
r_(n+1)=r_n exp(2*pi*i*x).                          (RE2)
```

It never sends a large unreduced phase to the platform trigonometric
functions.  At the start of every {DEFAULT_BLOCK}-term block, `z_n` and `r_n`
are reseeded from the exact rational input phases reduced modulo two.  Real
and imaginary parts of `S_0` and `S_1` are accumulated by blockwise
`math.fsum`.  The implementation is O(L), uses one below-normal worker, and
is the audit oracle for a later fast incomplete-Gauss evaluator.

The short-roster gate compares four nonperiodic decimal/rational phase rows,
through length 4096, against direct {SHORT_DPS}-digit term evaluation.  The
full-roster gate uses an independent exact compression.  For rational `x,s`,

```text
z_(n+P)/z_n
 =exp(i*pi[2xPn+xP(P+A+2s)]).                       (RE3)
```

Consequently `P` is a period whenever

```text
xP is an integer,
xP(P+A+2s) is an even integer.                      (RE4)
```

Writing `L=kP+R`, one cycle gives exactly

```text
S_0=k Z_0+Z_(0,R),
S_1=P k(k-1)Z_0/2+k Z_1+kP Z_(0,R)+Z_(1,R),        (RE5)
```

with the analogous affine-cycle formula applied directly to the original
weighted current `F_x(s)`.  Thus the seven full-roster rows are checked
without replaying 2,481,422 terms in the independent oracle.  Their periods
range from {summary['minimum_period']} to {summary['maximum_period']}.

All {summary['row_count']} rows pass the declared mass-normalized floating
allowance `{NORMALIZED_ALLOWANCE:.1e}`.  The largest observed normalized
errors are

```text
S_0: {summary['maximum_mass_normalized_errors']['S0']:.6e}
S_1: {summary['maximum_mass_normalized_errors']['S1']:.6e}
F:   {summary['maximum_mass_normalized_errors']['F']:.6e}.       (RE6)
```

These denominators are the corresponding absolute term masses, not a claim
that cancellation is absent.  Agreement at rational test phases validates
the implementation and catches long-run phase drift; it does not bound the
roundoff error uniformly for every physical `x,s`.

The next stage is to build a genuinely sublinear incomplete quadratic-Gauss
evaluator, cross-check it against this O(L) oracle on deterministic points,
and only then compose it with the pole-free `C_U,epsilon` cell kernel.  A
rigorous physical computation additionally needs interval phase reduction,
summation-error enclosures, `x,s` quadrature remainders, and the ordered
`M`-then-`epsilon` limit.

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Pi provenance: `pi` in (RE1)--(RE3) is inherited from the original Kummer
quadratic phase.  Rational periodicity only reduces that existing phase
modulo its ordinary `2*pi` period; it introduces no fitted or geometric
constant.

Proof boundary: floating implementation validation and exact rational-phase
periodicity/cycle identities only.  No uniform floating-error theorem,
interval-certified theta evaluator, fast incomplete-Gauss algorithm,
physical `x,s` quadrature, non-A bound, joined `R_after_A`, `R_Dir`, `Q_K-T`,
all-height theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is
proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependency = json.loads(DEPENDENCY.read_text(encoding="utf-8"))
    require(dependency.get("passed") is True, "non-A theta-current dependency is not passed")
    require(dependency["theta_current_certificate"]["roster"]["L"] == L, "source roster length drift")
    require(dependency["decision"]["reference_term_recurrence_proved"] is True, "term recurrence dependency drift")
    require(CHECKER.is_file(), "independent checker missing")

    short_rows = run_short_rows()
    periodic_rows = run_periodic_rows()
    all_rows = short_rows + periodic_rows
    maxima = {
        label: max(row["mass_normalized_errors"][label] for row in all_rows)
        for label in ("S0", "S1", "F")
    }
    require(max(maxima.values()) <= NORMALIZED_ALLOWANCE, "global reference-evaluator allowance failed")

    artifact = {
        "kind": STEM,
        "status": "floating_O_L_quadratic_theta_reference_evaluator_validated_against_high_precision_and_exact_periodic_oracles_interval_theorem_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "source_odd_roster": [A, B],
            "source_cell_count": L,
            "theta_cell": "0<=s<=1",
            "workers": 1,
        },
        "algorithm": {
            "recurrences": [
                "z_(n+1)=z_n*r_n",
                "r_(n+1)=r_n*exp(2*pi*i*x)",
            ],
            "phase_input": "exact Fraction parsed from deterministic decimal/rational text",
            "phase_reduction": "reduce every reseed exponent modulo 2 before one cos/sin evaluation",
            "block_size": DEFAULT_BLOCK,
            "block_reseeding": "z_n and r_n recomputed from exact rational phases at every block start",
            "summation": "separate real/imaginary blockwise math.fsum for S0 and S1",
            "complexity": "O(L) time and O(block_size) memory",
        },
        "periodic_oracle": {
            "period_identity": "z_(n+P)/z_n=exp(i*pi*(2*x*P*n+x*P*(P+A+2*s)))",
            "period_conditions": [
                "x*P is an integer",
                "x*P*(P+A+2*s) is an even integer",
            ],
            "S0_cycle_formula": "L=kP+R => S0=k*Z0+Z0_R",
            "S1_cycle_formula": "S1=P*k*(k-1)*Z0/2+k*Z1+k*P*Z0_R+Z1_R",
            "current_oracle": "apply the affine cycle formula directly to alpha_n*exp(i*pi*x*alpha_n^2/4)",
        },
        "short_high_precision_rows": short_rows,
        "full_roster_periodic_rows": periodic_rows,
        "summary": {
            "row_count": len(all_rows),
            "short_row_count": len(short_rows),
            "full_roster_row_count": len(periodic_rows),
            "minimum_period": min(row["minimal_period"] for row in periodic_rows),
            "maximum_period": max(row["minimal_period"] for row in periodic_rows),
            "mass_normalized_allowance": NORMALIZED_ALLOWANCE,
            "maximum_mass_normalized_errors": maxima,
        },
        "decision": {
            "short_high_precision_rows_passed": True,
            "full_roster_exact_periodic_rows_passed": True,
            "large_phase_trigonometric_calls_eliminated": True,
            "reference_evaluator_ready_as_fast_evaluator_oracle": True,
            "uniform_floating_error_theorem_proved": False,
            "interval_certified": False,
            "fast_incomplete_quadratic_Gauss_evaluator_built": False,
            "physical_x_s_quadrature_completed": False,
            "non_A_bound_proved": False,
            "R_after_A_bound_proved": False,
            "rh_implication": False,
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "blas_threads": 1,
            "process_priority": priority,
        },
        "dependencies": {
            "non_A_common_kernel_theta_current": {"path": relative(DEPENDENCY), "sha256": file_hash(DEPENDENCY)}
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "proof_boundary": "Floating O(L) implementation validation and exact rational-phase period/cycle identities only. No uniform floating-error theorem, interval-certified theta evaluator, fast incomplete-Gauss algorithm, physical x,s quadrature, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"validated {len(all_rows)} quadratic-theta reference rows; "
        f"max normalized error {max(maxima.values()):.3e}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
