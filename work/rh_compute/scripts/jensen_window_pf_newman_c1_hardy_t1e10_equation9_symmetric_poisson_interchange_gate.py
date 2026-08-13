#!/usr/bin/env python3
"""Certify symmetric finite-Poisson interchange for the Kummer alpha roster."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


SADDLE_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_finite_poisson_portcullis_saddle_reduction_gate.json"
PAPER = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_Brereton_2607.15310.pdf"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 110
A = 159577
B = 5122421
L = (B - A) // 2


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


def symbolic_reduction() -> dict[str, Any]:
    u, x, alpha, A_sym, m = sp.symbols("u x alpha A m", real=True)
    alpha_u = A_sym + 2 * u
    f = alpha_u * sp.exp(sp.I * sp.pi * x * alpha_u**2 / 4)
    f1 = sp.factor(sp.diff(f, u) / sp.exp(sp.I * sp.pi * x * alpha_u**2 / 4))
    f2 = sp.factor(sp.diff(f, u, 2) / sp.exp(sp.I * sp.pi * x * alpha_u**2 / 4))
    expected_f1 = 2 + sp.I * sp.pi * x * alpha_u**2
    expected_f2 = 6 * sp.I * sp.pi * x * alpha_u - sp.pi**2 * x**2 * alpha_u**3
    require(sp.simplify(f1 - expected_f1) == 0, "first alpha-roster derivative failed")
    require(sp.simplify(f2 - expected_f2) == 0, "second alpha-roster derivative failed")

    D = 2 * sp.pi * sp.I * m
    delta_f, delta_f1, r_plus, r_minus = sp.symbols("Delta_f Delta_f1 R_plus R_minus")
    i_plus = -delta_f / D - delta_f1 / D**2 + r_plus / D**2
    i_minus = i_plus.subs({m: -m, r_plus: r_minus})
    paired = sp.simplify(i_plus + i_minus)
    expected_pair = sp.simplify((-2 * delta_f1 + r_plus + r_minus) / D**2)
    require(sp.simplify(paired - expected_pair) == 0, "paired integration-by-parts identity failed")
    require(sp.simplify(sp.diff(paired, delta_f)) == 0, "1/m endpoint current did not cancel")

    beta0 = sp.beta(sp.Rational(3, 4), sp.Rational(3, 4))
    beta1 = sp.beta(sp.Rational(7, 4), sp.Rational(3, 4))
    beta2 = sp.beta(sp.Rational(11, 4), sp.Rational(3, 4))
    require(sp.simplify(sp.expand_func(beta1 - beta0 / 2)) == 0, "first beta recurrence failed")
    require(sp.simplify(sp.expand_func(beta2 - 7 * beta0 / 20)) == 0, "second beta recurrence failed")

    return {
        "odd_roster_coordinate": "alpha=A+2u, 0<=u<=L=(B-A)/2",
        "alpha_roster_integrand": "f_x(u)=(A+2u) exp(i*pi*x*(A+2u)^2/4)",
        "first_derivative": "f_x'(u)=exp(i*pi*x*alpha^2/4)[2+i*pi*x*alpha^2]",
        "second_derivative": "f_x''(u)=exp(i*pi*x*alpha^2/4)[6*i*pi*x*alpha-pi^2*x^2*alpha^3]",
        "first_derivative_residual": str(sp.simplify(f1 - expected_f1)),
        "second_derivative_residual": str(sp.simplify(f2 - expected_f2)),
        "symmetric_partial_sum": "K_(t,M)=int_0^1 W_t(x){[f_x(0)+f_x(L)]/2+sum_(m=-M)^M I_m(x)}dx",
        "mode_integral": "I_m(x)=int_0^L f_x(u) exp(-2*pi*i*m*u)du",
        "paired_ibp_identity": "I_m+I_(-m)=[-2*Delta f_x'+R_m+R_(-m)]/(2*pi*i*m)^2",
        "paired_endpoint_cancellation": "The Delta f_x/(2*pi*i*m) terms cancel exactly before absolute values.",
        "pair_bound": "|I_m+I_(-m)|<=(|Delta f_x'|+int_0^L|f_x''(u)|du)/(2*pi^2*m^2)",
        "pointwise_derivative_majorant": "C_AB(x)=4+pi*x*(5*B^2-A^2)/2+pi^2*x^2*(B^4-A^4)/8",
        "beta_recurrences": "B(7/4,3/4)=B(3/4,3/4)/2 and B(11/4,3/4)=7*B(3/4,3/4)/20",
        "integrated_constant": "C_AB=B(3/4,3/4)[4+pi*(5*B^2-A^2)/4+7*pi^2*(B^4-A^4)/160]",
        "uniform_tail_bound": "|K_t-K_(t,M)|<=C_AB/(2*pi^2*M), uniformly for every real t and M>=1",
        "kummer_weight_absolute_value": "|W_t(x)|=[x(1-x)]^(-1/4)",
        "interchange_reason": "The Kummer endpoint weight is absolutely integrable and the paired Fourier tail has an integrable O(M^-1) majorant independent of t.",
    }


def certified_constants() -> dict[str, Any]:
    ctx.dps = PRECISION
    q = arb(3) / 4
    beta0 = q.gamma() ** 2 / (arb(3) / 2).gamma()
    a = arb(A)
    b = arb(B)
    pi = arb.pi()
    bracket = 4 + pi * (5 * b**2 - a**2) / 4 + 7 * pi**2 * (b**4 - a**4) / 160
    c_ab = beta0 * bracket
    tail_coefficient = c_ab / (2 * pi**2)
    cutoff_for_005 = tail_coefficient / arb("0.005")
    require(beta0 > 0 and c_ab > 0 and tail_coefficient > 0, "certified constants are not positive")
    require(tail_coefficient > arb("1e25"), "raw tail coefficient unexpectedly small")
    return {
        "precision_decimal_digits": PRECISION,
        "A": A,
        "B": B,
        "L": L,
        "alpha_count": L + 1,
        "beta_3_4_3_4_ball": beta0.str(PRECISION, more=True),
        "integrated_C_AB_ball": c_ab.str(PRECISION, more=True),
        "tail_coefficient_C_AB_over_2pi2_ball": tail_coefficient.str(PRECISION, more=True),
        "raw_cutoff_required_for_bound_below_0_005_ball": cutoff_for_005.str(PRECISION, more=True),
        "scale_decision": "The C2 interchange proof is rigorous and height-uniform, but its absolute tail coefficient is too large for a quantitative source-error bound. The stationary/nonstationary phase structure must be used before absolute values.",
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certified_constants"]
    return f"""# Symmetric finite-Poisson interchange for the Kummer roster

Date: 2026-08-10

Status: exact height-uniform interchange theorem validated; not a proof of the stationary-phase remainder

Write the contiguous odd roster as `alpha=A+2u`, `u=0..L`, and set

```text
f_x(u)=(A+2u) exp(i*pi*x*(A+2u)^2/4),
W_t(x)=exp[i(t/2)log((1-x)/x)]/[x(1-x)]^(1/4).
```

For a symmetric Fourier cutoff `M`, finite Poisson gives the partial
reconstruction

```text
K_(t,M)=int_0^1 W_t(x)
  {{[f_x(0)+f_x(L)]/2+sum_(m=-M)^M I_m(x)}} dx,
I_m(x)=int_0^L f_x(u)exp(-2*pi*i*m*u)du.
```

Two integrations by parts, followed by pairing `m` and `-m`, give

```text
I_m+I_(-m)=[-2 Delta f_x'+R_m+R_(-m)]/(2*pi*i*m)^2.  (SP1)
```

The apparent `Delta f_x/(2*pi*i*m)` endpoint current cancels exactly in
(SP1).  Consequently

```text
|I_m+I_(-m)|
 <= (|Delta f_x'|+int_0^L |f_x''(u)|du)/(2*pi^2*m^2).
```

The exact derivatives are

```text
f_x'=exp(i*pi*x*alpha^2/4)[2+i*pi*x*alpha^2],
f_x''=exp(i*pi*x*alpha^2/4)
       [6*i*pi*x*alpha-pi^2*x^2*alpha^3].
```

Since `|W_t(x)|=[x(1-x)]^(-1/4)`, beta integration yields, uniformly for
every real `t` and every integer `M>=1`,

```text
|K_t-K_(t,M)| <= C_AB/(2*pi^2*M),                       (SP2)

C_AB=B(3/4,3/4)
 [4+pi(5B^2-A^2)/4+7pi^2(B^4-A^4)/160].                (SP3)
```

Thus no auxiliary endpoint regulator is required: the Kummer endpoint
singularity is absolutely integrable, and the paired Fourier tail supplies
an integrable majorant independent of height.  This closes the interchange
part of the previous obligation.

For the source roster `A={A}`, `B={B}`, `L={L}`, interval arithmetic gives

```text
C_AB/(2*pi^2) = {c['tail_coefficient_C_AB_over_2pi2_ball']}.
```

This last number is a structural warning.  The raw `C^2` triangle bound would
need a cutoff larger than

```text
{c['raw_cutoff_required_for_bound_below_0_005_ball']}
```

merely to fall below `0.005`.  It proves convergence and interchange, but is
not a viable quantitative Hardy-error estimate.  The next step must exploit
the joint Kummer/Poisson phase, retain the modewise boundary/Fresnel pairing,
and obtain phase-adapted bounds for the interior, turning, and nonstationary
mode ranges before taking absolute values.

Proof boundary: exact finite-roster Poisson reconstruction and a uniform
interchange/tail theorem only.  No useful stationary-phase constant,
source-aligned hybrid remainder, `Lambda<=0`, RH, or prize-level conclusion is
proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    require((B - A) % 2 == 0 and A % 2 == B % 2 == 1, "invalid odd roster")
    require(SADDLE_GATE.is_file() and PAPER.is_file() and CHECKER.is_file(), "missing dependency or checker")
    symbolic = symbolic_reduction()
    constants = certified_constants()
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate",
        "status": "exact_height_uniform_symmetric_poisson_interchange_complete_phase_adapted_remainder_open",
        "passed": True,
        "scope": {
            "height": "every real t for a fixed finite odd-alpha roster",
            "source_roster": f"{A}..{B} odd",
            "diagnostic_midpoint_used": False,
        },
        "symbolic_reduction": symbolic,
        "certified_constants": constants,
        "decision": {
            "symmetric_poisson_interchange_proved": True,
            "auxiliary_endpoint_regulator_required": False,
            "paired_one_over_m_endpoint_current_cancels": True,
            "tail_bound_uniform_in_height": True,
            "raw_c2_tail_bound_quantitatively_viable": False,
            "stationary_phase_remainder_proved": False,
            "rh_implication": False,
        },
        "next_obligation": "Derive phase-adapted bounds for the joint Kummer/Poisson integral: interior saddles, all 84 lower-alpha turning modes, and symmetric nonstationary tails, while retaining each mode's grouped boundary/Fresnel current.",
        "proof_boundary": "Exact symmetric finite-Poisson reconstruction and height-uniform interchange for a fixed finite odd-alpha roster. No quantitatively useful stationary-phase remainder, source-aligned hybrid bound, Lambda<=0, RH, or prize-level conclusion is proved.",
        "dependencies": {
            "finite_poisson_saddle_gate": {"path": relative(SADDLE_GATE), "sha256": file_hash(SADDLE_GATE)},
            "paper": {"path": relative(PAPER), "sha256": file_hash(PAPER)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "sympy_threads": 1,
            "flint_threads": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        "built symmetric finite-Poisson interchange: "
        f"alpha={constants['alpha_count']}, uniform_tail=True, quantitative=False"
    )


if __name__ == "__main__":
    main()
