#!/usr/bin/env python3
"""Specialize the truncated-theta Poisson recursion to the non-A current."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import importlib.util
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
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_poisson_current_self_duality_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
REFERENCE_BUILDER = BUILDER.with_name(
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_reference_evaluator_gate.py"
)
DEPENDENCIES = {
    "common_kernel_theta_current": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_common_kernel_theta_current_gate.json",
    "reference_evaluator": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_reference_evaluator_gate.json",
}

A = 159_577
B = 5_122_421
L = 2_481_422
K = L - 1
DIRECT_STOP = 64
DUAL_BLOCK = 256

TRACE_CASES = (
    {"x": "1/2", "s": "0"},
    {"x": "1/2", "s": "1/2"},
    {"x": "1/3", "s": "1/3"},
    {"x": "2/5", "s": "3/8"},
    {"x": "1/8", "s": "7/10"},
    {"x": "1/4096", "s": "1/2"},
    {"x": "123457/1000000", "s": "49/100"},
    {"x": "1/10000000", "s": "1/3"},
)

PILOT_CASES = (
    {"x": "1/2", "s": "0"},
    {"x": "2/5", "s": "1/3"},
    {"x": "1/3", "s": "1/2"},
    {"x": "1/8", "s": "3/8"},
    {"x": "1/4096", "s": "1/2"},
    {"x": "1/1000000", "s": "49/100"},
)


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
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def load_reference_module():
    spec = importlib.util.spec_from_file_location("non_a_theta_reference", REFERENCE_BUILDER)
    require(spec is not None and spec.loader is not None, "reference module loader failed")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def normalize_pair(a: Fraction, b: Fraction) -> tuple[Fraction, Fraction, bool, str]:
    """Use the exact theta periodicities to put b in [0,1/4]."""

    beta = b % 1
    conjugated = False
    if beta <= Fraction(1, 4):
        a_out, b_out, branch = a, beta, "unit_period"
    elif beta < Fraction(1, 2):
        a_out, b_out = Fraction(1, 2) - a, Fraction(1, 2) - beta
        conjugated, branch = True, "half_shift_then_conjugate"
    elif beta <= Fraction(3, 4):
        a_out, b_out, branch = a - Fraction(1, 2), beta - Fraction(1, 2), "half_shift"
    else:
        a_out, b_out = -a, 1 - beta
        conjugated, branch = True, "unit_shift_then_conjugate"
    a_out %= 1
    require(0 <= a_out < 1 and 0 <= b_out <= Fraction(1, 4), "normalization range failed")
    return a_out, b_out, conjugated, branch


def symbolic_certificate() -> dict[str, str]:
    n, v, x, c = sp.symbols("n v x c", real=True, positive=True)
    phase = x * (c + 2 * n) ** 2 / 4 - 2 * v * n
    stationary = sp.simplify(v / x - c / 2)
    require(sp.simplify(sp.diff(phase, n).subs(n, stationary)) == 0, "stationary point drift")
    require(sp.simplify((c + 2 * n).subs(n, stationary) - 2 * v / x) == 0, "current self-duality drift")
    require(sp.simplify(phase.subs(n, stationary) - (c * v - v**2 / x)) == 0, "dual phase drift")
    require(sp.simplify(sp.diff(phase, n, 2) - 2 * x) == 0, "quadratic curvature drift")

    a, b, q, kk = sp.symbols("a b q K", real=True, positive=True)
    c0 = sp.symbols("C_0")
    f0, f1, r0, r1 = sp.symbols("F_0 F_1 R_0 R_1")
    transformed_zero = c0 * f0 + r0
    transformed_first = c0 * (q * f1 - a * f0) / (2 * b * kk) + r1
    source = sp.expand(c * transformed_zero + 2 * kk * transformed_first)
    expected = sp.expand(c0 * ((c - a / b) * f0 + q * f1 / b) + c * r0 + 2 * kk * r1)
    require(sp.simplify(source - expected) == 0, "normalized current coefficient transport drift")
    return {
        "physical_phase": "phi(n;v)=x(c+2n)^2/4-2vn, c=A+2s",
        "stationary_point": "n_v=v/x-c/2",
        "stationary_current": "c+2n_v=2v/x",
        "stationary_phase": "phi(n_v;v)=cv-v^2/x",
        "gaussian_factor": "exp(i*pi/4)/sqrt(x)",
        "pure_dual_current": "F_saddle(x,s)=2*exp(i*pi/4)*x^(-3/2)*sum_(ceil(xc/2)<=v<=floor(xc/2+xK)) v*exp(i*pi*(cv-v^2/x))",
        "normalized_F0": "F(K,0)=C_0 F(q,0;a*,b*)+R_0",
        "normalized_F1": "F(K,1)=C_0[q F(q,1;a*,b*)-a F(q,0;a*,b*)]/(2bK)+R_1",
        "normalized_source": "cS_0+2S_1=C_0[(c-a/b)F(q,0)+(q/b)F(q,1)]+cR_0+2KR_1",
    }


def normalization_identity_certificate() -> list[dict[str, Any]]:
    cases = (
        (Fraction(7, 19), Fraction(3, 20)),
        (Fraction(7, 19), Fraction(7, 20)),
        (Fraction(7, 19), Fraction(13, 20)),
        (Fraction(7, 19), Fraction(17, 20)),
    )
    rows = []
    for a, b in cases:
        a0, b0, conjugated, branch = normalize_pair(a, b)
        for n in range(24):
            original = (a * n + b * n * n) % 1
            normalized = (a0 * n + b0 * n * n) % 1
            if conjugated:
                require((original + normalized).denominator == 1, f"conjugate phase normalization failed: {branch}")
            else:
                require((original - normalized).denominator == 1, f"phase normalization failed: {branch}")
        rows.append(
            {
                "a": str(a),
                "b": str(b),
                "a_normalized": str(a0),
                "b_normalized": str(b0),
                "conjugated": conjugated,
                "branch": branch,
                "integer_phase_checks": 24,
            }
        )
    return rows


def trace_case(x: Fraction, s: Fraction) -> dict[str, Any]:
    c = A + 2 * s
    a, b, conjugated, branch = normalize_pair(x * c / 2, x / 2)
    current_k = K
    rows = []
    for iteration in range(32):
        p = ceil_fraction(a)
        q = floor_fraction(a + 2 * b * current_k)
        row: dict[str, Any] = {
            "iteration": iteration,
            "K": current_k,
            "a": str(a),
            "b": str(b),
            "p": p,
            "q": q,
            "normalization_branch": branch,
            "conjugated": conjugated,
        }
        if current_k <= DIRECT_STOP:
            row["stop"] = "direct_short_sum"
            rows.append(row)
            break
        if b == 0 or q <= p:
            row["stop"] = "Euler_Maclaurin_or_linear_boundary"
            rows.append(row)
            break
        require(q <= (current_k + 1) // 2, "normalized recursion failed to contract by one half")
        row["contracted_K"] = q
        row["contraction_ratio"] = q / current_k
        rows.append(row)
        a_star = a / (2 * b)
        b_star = -Fraction(1, 1) / (4 * b)
        a, b, conjugated, branch = normalize_pair(a_star, b_star)
        current_k = q
    else:
        raise RuntimeError("theta recursion trace exceeded 32 steps")
    return {
        "x": str(x),
        "s": str(s),
        "initial_a": str(x * c / 2),
        "initial_b": str(x / 2),
        "iterations": len(rows),
        "main_steps": sum("contracted_K" in row for row in rows),
        "terminal_K": rows[-1]["K"],
        "terminal_reason": rows[-1]["stop"],
        "maximum_contraction_ratio": max((row.get("contraction_ratio", 0.0) for row in rows), default=0.0),
        "rows": rows,
    }


def dual_saddle_current(reference: Any, x: Fraction, s: Fraction) -> tuple[complex, int, int]:
    c = A + 2 * s
    a = x * c / 2
    p = ceil_fraction(a)
    q = floor_fraction(a + x * K)
    if q < p:
        return 0j, p, q
    q_ratio = reference.cis_pi_fraction(-Fraction(2, 1) / x)
    block_real: list[float] = []
    block_imag: list[float] = []
    for start in range(p, q + 1, DUAL_BLOCK):
        stop = min(q + 1, start + DUAL_BLOCK)
        z = reference.cis_pi_fraction(c * start - Fraction(start * start, 1) / x)
        ratio = reference.cis_pi_fraction(c - Fraction(2 * start + 1, 1) / x)
        local_real: list[float] = []
        local_imag: list[float] = []
        for v in range(start, stop):
            local_real.append(v * z.real)
            local_imag.append(v * z.imag)
            z *= ratio
            ratio *= q_ratio
        block_real.append(math.fsum(local_real))
        block_imag.append(math.fsum(local_imag))
    weighted = complex(math.fsum(block_real), math.fsum(block_imag))
    prefactor = 2.0 * complex(math.cos(math.pi / 4), math.sin(math.pi / 4)) / float(x) ** 1.5
    return prefactor * weighted, p, q


def run_pilot(reference: Any) -> list[dict[str, Any]]:
    rows = []
    for case in PILOT_CASES:
        x = Fraction(case["x"])
        s = Fraction(case["s"])
        full = reference.recurrence_theta_current(x, s, L)[2]
        saddle, p, q = dual_saddle_current(reference, x, s)
        remainder = full - saddle
        mass = float(Fraction(L) * (A + 2 * s + L - 1))
        rows.append(
            {
                **case,
                "stationary_integer_range": [p, q],
                "stationary_term_count": max(0, q - p + 1),
                "compression_ratio": max(0, q - p + 1) / L,
                "reference_current": {"real": full.real, "imag": full.imag},
                "completed_saddle_current": {"real": saddle.real, "imag": saddle.imag},
                "endpoint_and_nonstationary_remainder": {"real": remainder.real, "imag": remainder.imag},
                "remainder_absolute": abs(remainder),
                "remainder_over_absolute_term_mass": abs(remainder) / mass,
                "remainder_over_reference_magnitude": abs(remainder) / max(1.0, abs(full)),
                "diagnostic_only": True,
            }
        )
    return rows


def render_note(artifact: dict[str, Any]) -> str:
    summary = artifact["pilot_summary"]
    return f"""# Quadratic-theta Poisson current self-duality

Date: 2026-08-24

Status: exact truncated-theta specialization, normalization, contraction, and
pure first-moment saddle transform certified; endpoint-complete fast evaluator
open

Put `K=L-1={K}`, `c=A+2s`,

```text
a=xc/2,  b=x/2,

F(K,j;a,b)=K^(-j) sum_(n=0)^K n^j exp(2pi i[an+bn^2]). (PS1)
```

Then exactly `S_0=F(K,0;a,b)` and `S_1=K F(K,1;a,b)`.  This puts the
non-A current in the `j=0,1` case of G. A. Hiary's
[truncated-theta algorithm](https://arxiv.org/abs/0711.5002v4).  The source is
used for the normalization identities and endpoint-complete recursive
architecture; none of its unnamed constants is imported into a corpus bound.

On the physical interval `0<x<=1/2`, `0<b<=1/4` already.  Reducing only `a`
modulo one gives

```text
q=floor(a+2bK)=floor(a+xK)<=(K+1)/2.               (PS2)
```

Every main recursion step therefore contracts by at least one half after the
exact unit/half-period and conjugation normalization.  The eight exact
rational traces terminate in at most {artifact['trace_summary']['maximum_iterations']}
recorded rows, with at most {artifact['trace_summary']['maximum_main_steps']}
main Poisson steps.  This is a parameter-map result; the required remainders
have not been implemented.

The special current has a stronger symmetry than a generic `j=0,1` pair.
Apply Poisson frequency `v` directly to its full quadratic phase:

```text
phi(n;v)=x(c+2n)^2/4-2vn,
n_v=v/x-c/2.                                        (PS3)
```

At the stationary point,

```text
c+2n_v=2v/x,
phi(n_v;v)=cv-v^2/x.                                (PS4)
```

Thus the huge constant phase cancels and the completed stationary channel is
the pure absolute-dual first moment

```text
F_saddle(x,s)
 =2 exp(i*pi/4)x^(-3/2)
   sum_(ceil(xc/2)<=v<=floor(xc/2+xK))
      v exp(i*pi[cv-v^2/x]).                        (PS5)
```

There is no independent zeroth-order dual current in (PS5).  If `v=p+r`, it
reappears only in the tied combination `p F_0+N F_1`, where `N=q-p`; splitting
those two terms before estimating would destroy the derivative-current
cancellation.

In the normalized Hiary coordinate the same fact is

```text
F(K,0)=C_0 F(q,0;a*,b*)+R_0,

F(K,1)=C_0[qF(q,1;a*,b*)-aF(q,0;a*,b*)]/(2bK)+R_1,

cS_0+2S_1
 =C_0[(c-a/b)F(q,0)+(q/b)F(q,1)]
  +cR_0+2KR_1.                                     (PS6)
```

The six floating rows compare (PS5) alone with the validated full-roster
reference current.  The stationary roster uses between
{summary['minimum_stationary_terms']} and {summary['maximum_stationary_terms']}
terms.  The observed omitted endpoint-plus-nonstationary remainder, divided
by absolute source term mass, ranges from
`{summary['minimum_mass_normalized_remainder']:.6e}` to
`{summary['maximum_mass_normalized_remainder']:.6e}`.  These values are route
diagnostics only.  In particular, (PS5) is not promoted as an approximation
until `R_0,R_1` are evaluated with explicit error.

The next implementation stage is narrow: specialize the endpoint-complete
`R_0,R_1` formulas to `j<=1`, preserve the tied current `cR_0+2KR_1`, and
cross-check one complete recursion step against the O(`L`) reference oracle.
The small-`b` branch must retain its Euler--Maclaurin integral and correction
terms rather than silently returning the saddle sum.

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Pi provenance: all `pi` factors in (PS1)--(PS6) come from the inherited
Kummer quadratic phase and the ordinary integer Fourier character in Poisson
summation.  No fitted or geometric occurrence is introduced.

Proof boundary: exact parameter mapping, theta normalization identities,
one-half recursion contraction, and completed stationary-current self-duality
only.  The six remainder values are floating diagnostics.  No
endpoint/nonstationary remainder evaluator, complete fast theta algorithm,
uniform or interval error theorem, physical quadrature, non-A bound, joined
`R_after_A`, `R_Dir`, `Q_K-T`, all-height theorem, `Lambda<=0`, PF-infinity,
RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a self-duality dependency is not passed")
    require(dependencies["common_kernel_theta_current"]["scope"]["source_cell_count"] == L, "source length drift")
    require(dependencies["reference_evaluator"]["decision"]["reference_evaluator_ready_as_fast_evaluator_oracle"] is True, "reference oracle drift")
    require(CHECKER.is_file() and REFERENCE_BUILDER.is_file(), "checker or reference builder missing")

    reference = load_reference_module()
    symbolic = symbolic_certificate()
    normalizations = normalization_identity_certificate()
    traces = [trace_case(Fraction(case["x"]), Fraction(case["s"])) for case in TRACE_CASES]
    pilot = run_pilot(reference)
    trace_summary = {
        "row_count": len(traces),
        "maximum_iterations": max(row["iterations"] for row in traces),
        "maximum_main_steps": max(row["main_steps"] for row in traces),
        "maximum_observed_contraction_ratio": max(row["maximum_contraction_ratio"] for row in traces),
        "all_terminal_reasons": sorted({row["terminal_reason"] for row in traces}),
    }
    require(trace_summary["maximum_observed_contraction_ratio"] <= (K + 1) / (2 * K), "trace contraction summary drift")
    pilot_summary = {
        "row_count": len(pilot),
        "minimum_stationary_terms": min(row["stationary_term_count"] for row in pilot),
        "maximum_stationary_terms": max(row["stationary_term_count"] for row in pilot),
        "minimum_compression_ratio": min(row["compression_ratio"] for row in pilot),
        "maximum_compression_ratio": max(row["compression_ratio"] for row in pilot),
        "minimum_mass_normalized_remainder": min(row["remainder_over_absolute_term_mass"] for row in pilot),
        "maximum_mass_normalized_remainder": max(row["remainder_over_absolute_term_mass"] for row in pilot),
    }

    artifact = {
        "kind": STEM,
        "status": "exact_quadratic_theta_poisson_current_self_duality_and_recursion_contraction_certified_endpoint_complete_remainder_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "A": A,
            "B": B,
            "L": L,
            "K": K,
            "physical_x": "0<=x<=1/2",
            "physical_s": "0<=s<=1",
            "Hiary_moments": [0, 1],
        },
        "primary_source": {
            "author": "Ghaith Ayesh Hiary",
            "title": "A nearly-optimal method to compute the truncated theta function, its derivatives, and integrals",
            "arxiv": "0711.5002v4",
            "url": "https://arxiv.org/abs/0711.5002v4",
            "used_for": "Lemma 4.1 normalization identities and the endpoint-complete recursive architecture for F(K,j;a,b)",
            "imported_numerical_constants": False,
        },
        "theta_mapping": {
            "definition": "F(K,j;a,b)=K^(-j) sum_(n=0)^K n^j exp(2*pi*i*(a*n+b*n^2))",
            "parameters": "K=L-1, c=A+2s, a=x*c/2, b=x/2",
            "S0": "S_0=F(K,0;a,b)",
            "S1": "S_1=K*F(K,1;a,b)",
            "physical_current": "F_x(s)=exp(i*pi*x*c^2/4)*(c*S_0+2*S_1)",
        },
        "symbolic_certificate": symbolic,
        "normalization_identity_rows": normalizations,
        "recursion_traces": traces,
        "trace_summary": trace_summary,
        "floating_saddle_only_rows": pilot,
        "pilot_summary": pilot_summary,
        "decision": {
            "non_A_theta_current_is_Hiary_j0_j1_instance": True,
            "physical_quadratic_parameter_already_in_contracting_range": True,
            "every_main_step_contracts_length_by_at_least_one_half": True,
            "completed_saddle_current_is_pure_absolute_dual_first_moment": True,
            "large_constant_phase_cancels_in_completed_saddle_current": True,
            "zeroth_and_first_reindexed_dual_moments_must_remain_tied": True,
            "endpoint_complete_R0_R1_evaluator_built": False,
            "small_b_Euler_Maclaurin_branch_built": False,
            "complete_fast_theta_evaluator_built": False,
            "uniform_error_theorem_proved": False,
            "physical_quadrature_completed": False,
            "non_A_bound_proved": False,
            "R_after_A_bound_proved": False,
            "rh_implication": False,
        },
        "next_obligation": "Specialize the endpoint-complete R_0 and R_1 terms to j<=1, preserve cR_0+2KR_1 as one current, implement the small-b Euler--Maclaurin boundary branch, and cross-check one complete modular step against the validated O(L) oracle before recursion.",
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "blas_threads": 1,
            "process_priority": priority,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)} for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "reference_builder": {"path": relative(REFERENCE_BUILDER), "sha256": file_hash(REFERENCE_BUILDER)},
        },
        "proof_boundary": "Exact parameter mapping, theta normalization identities, one-half recursion contraction, and completed stationary-current self-duality only. The saddle-only remainder rows are floating diagnostics. No endpoint/nonstationary remainder evaluator, small-b Euler--Maclaurin branch, complete fast theta algorithm, uniform or interval error theorem, physical quadrature, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"certified quadratic-theta Poisson current self-duality; "
        f"{len(traces)} exact traces, {len(pilot)} diagnostic saddle rows",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
