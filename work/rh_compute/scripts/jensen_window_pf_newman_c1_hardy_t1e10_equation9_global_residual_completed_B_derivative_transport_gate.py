#!/usr/bin/env python3
"""Derive exact first and second derivative transport for the completed B current."""

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

import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_derivative_transport_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "endpoint_decomposition": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_finite_endpoint_tail_decomposition_gate.json",
    "crossing_completion": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_endpoint_safe_completed_B_extraction_gate.json",
    "tangential_cutoff": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_tangential_exterior_cutoff_gate.json",
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
    x, mode, endpoint = sp.symbols("x m B", positive=True)
    pi = sp.pi
    imaginary_unit = sp.I
    rho = 1 / (imaginary_unit * pi)
    a = mode * sp.sqrt(2 / x)
    b = endpoint * sp.sqrt(x / 2)
    r_value = sp.symbols("R")

    sign_rows: list[dict[str, str]] = []
    for sign in (1, -1):
        a_sign = sign * a
        q = b - sign * a
        q_first = sp.simplify(sp.diff(q, x))
        q_second = sp.simplify(sp.diff(q_first, x))
        q_first_expected = (b + sign * a) / (2 * x)
        q_second_expected = (-b - 3 * sign * a) / (4 * x**2)
        require(sp.simplify(q_first - q_first_expected) == 0, "q first derivative failed")
        require(sp.simplify(q_second - q_second_expected) == 0, "q second derivative failed")

        kappa = sp.simplify(1 / (2 * x) + imaginary_unit * pi * q * q_first)
        kappa_first = sp.simplify(sp.diff(kappa, x))
        kappa_first_expected = sp.simplify(
            -1 / (2 * x**2) + imaginary_unit * pi * (q_first**2 + q * q_second)
        )
        require(sp.simplify(kappa_first - kappa_first_expected) == 0, "K first derivative failed")

        y = rho - a_sign * r_value
        r_q = -imaginary_unit * pi * q * r_value - 1
        y_first_chain = sp.simplify(-sp.diff(a_sign, x) * r_value - a_sign * r_q * q_first)
        y_first_transport = sp.simplify(a_sign * q_first + (rho - y) * kappa)
        require(sp.simplify(y_first_chain - y_first_transport) == 0, "Y transport failed")

        y_symbol = sp.symbols(f"Y_{'plus' if sign == 1 else 'minus'}")
        y_first_symbol = a_sign * q_first + (rho - y_symbol) * kappa
        y_second_transport = sp.simplify(
            a_sign * (q_second - q_first / (2 * x))
            - y_first_symbol * kappa
            + (rho - y_symbol) * kappa_first
        )
        sign_rows.append(
            {
                "sign": "+" if sign == 1 else "-",
                "a_sign": str(a_sign),
                "q": str(q),
                "q_first": str(q_first_expected),
                "q_second": str(q_second_expected),
                "K": str(kappa),
                "K_first": str(kappa_first_expected),
                "Y_first": str(y_first_symbol),
                "Y_second": str(y_second_transport),
            }
        )

    # The non-target bulk completion is a homogeneous solution V'=-K_+ V.
    q_plus = b - a
    q_plus_first = sp.diff(q_plus, x)
    completion = a * (1 + imaginary_unit) * sp.exp(-imaginary_unit * pi * q_plus**2 / 2)
    kappa_plus = 1 / (2 * x) + imaginary_unit * pi * q_plus * q_plus_first
    require(sp.simplify(sp.diff(completion, x) + kappa_plus * completion) == 0, "bulk homogeneous transport failed")

    y_plus = sp.Function("Y_plus")(x)
    y_minus = sp.Function("Y_minus")(x)
    pair = (y_plus + y_minus) / x
    pair_first = sp.diff(pair, x)
    pair_second = sp.diff(pair, x, 2)
    pair_first_expected = (sp.diff(y_plus, x) + sp.diff(y_minus, x)) / x - (y_plus + y_minus) / x**2
    pair_second_expected = (
        (sp.diff(y_plus, x, 2) + sp.diff(y_minus, x, 2)) / x
        - 2 * (sp.diff(y_plus, x) + sp.diff(y_minus, x)) / x**2
        + 2 * (y_plus + y_minus) / x**3
    )
    require(sp.simplify(pair_first - pair_first_expected) == 0, "pair first derivative failed")
    require(sp.simplify(pair_second - pair_second_expected) == 0, "pair second derivative failed")

    weight = (x * (1 - x)) ** sp.Rational(-1, 4)
    log_first = sp.simplify(sp.diff(weight, x) / weight)
    weight_second_ratio = sp.simplify(sp.diff(weight, x, 2) / weight)
    require(log_first == (2 * x - 1) / (4 * x * (1 - x)), "weight first derivative failed")
    require(
        sp.simplify(weight_second_ratio - (12 * x**2 - 12 * x + 5) / (16 * x**2 * (1 - x) ** 2)) == 0,
        "weight second derivative failed",
    )

    return {
        "phase_stripped_modes": {
            "Y_plus_tau": "exp(-i*pi*q_+^2/2)[P_B,+ +(1-tau_m)P_bulk,+]",
            "Y_minus": "exp(-i*pi*q_-^2/2)P_B,-",
            "q_sign": "q_s=B sqrt(x/2)-s m sqrt(2/x), s in {+1,-1}",
            "a_sign": "a_s=s m sqrt(2/x)",
            "rho": "1/(i*pi)",
        },
        "transport": "Y_s'=a_s q_s'+(rho-Y_s)K_s, K_s=1/(2x)+i*pi*q_s*q_s'",
        "second_transport": "Y_s''=a_s(q_s''-q_s'/(2x))-Y_s'K_s+(rho-Y_s)K_s'",
        "geometry_derivatives": "q_s'=(b+s a)/(2x), q_s''=(-b-3s a)/(4x^2), K_s'=-1/(2x^2)+i*pi[(q_s')^2+q_s q_s'']",
        "bulk_completion_homogeneous_equation": "V_+'=-K_+ V_+, V_+=(1-tau_m)a(1+i)exp(-i*pi*q_+^2/2)",
        "pair_coefficient": "C_m=(Y_plus_tau+Y_minus)/x",
        "pair_derivatives": {
            "first": "C_m'=(Y_+'+Y_-')/x-(Y_++Y_-)/x^2",
            "second": "C_m''=(Y_+''+Y_-'')/x-2(Y_+'+Y_-')/x^2+2(Y_++Y_-)/x^3",
        },
        "Kummer_weight": "w=[x(1-x)]^(-1/4), w'/w=(2x-1)/[4x(1-x)], w''/w=(12x^2-12x+5)/[16x^2(1-x)^2]",
        "finite_sum_amplitude": "A_(M,epsilon)=w sum_(m=1)^M w_(m,epsilon) C_m; differentiate the finite sum using the displayed transports before taking norms",
        "normalized_coordinate": "partial_xi=H_B^(-1/2)partial_x and partial_xixi=H_B^(-1)partial_xx",
        "sign_rows": sign_rows,
    }


def render_note() -> str:
    return """# Exact derivative transport for the crossing-completed B current

Date: 2026-08-13

Status: exact derivative reduction; not a proof artifact

No quantitative derivative norm or grouped-remainder bound is asserted.

After odd-endpoint parity, factor the common B phase.  For `s=+1,-1` put

```text
a=m sqrt(2/x),                 b=B sqrt(x/2),
a_s=s a,                       q_s=b-s a,
rho=1/(i*pi),

Y_+^tau=e^(-i*pi*q_+^2/2)[P_B,++(1-tau_m)P_bulk,+],
Y_-    =e^(-i*pi*q_-^2/2)P_B,-.                       (DT1)
```

The added non-target bulk term in `Y_+^tau` is precisely the crossing
completion from Section 11.417.  Let

```text
K_s=1/(2x)+i*pi*q_s q_s'.                             (DT2)
```

Using only `F'(q)=exp(i*pi*q^2/2)`, the phase-stripped endpoint tail obeys
the closed transport equation

```text
Y_s'=a_s q_s'+(rho-Y_s)K_s,                           (DT3)

q_s' =(b+s a)/(2x),
q_s''=(-b-3s a)/(4x^2),
K_s' =-1/(2x^2)+i*pi[(q_s')^2+q_s q_s''].             (DT4)
```

The bulk completion is not an extra forcing term.  Its phase-stripped value

```text
V_+=(1-tau_m)a(1+i)e^(-i*pi*q_+^2/2)
```

satisfies exactly `V_+'=-K_+V_+`; therefore the completed `Y_+^tau` obeys
the same inhomogeneous equation (DT3).  Differentiating once more gives

```text
Y_s''=a_s(q_s''-q_s'/(2x))-Y_s'K_s+(rho-Y_s)K_s'.     (DT5)
```

The complete phase-stripped mode pair is

```text
C_m=(Y_+^tau+Y_-)/x,
C_m'=(Y_+'+Y_-')/x-(Y_++Y_-)/x^2,
C_m''=(Y_+''+Y_-'')/x-2(Y_+'+Y_-')/x^2
       +2(Y_++Y_-)/x^3.                               (DT6)
```

Finally, for `w(x)=[x(1-x)]^(-1/4)`,

```text
w'/w =(2x-1)/[4x(1-x)],
w''/w=(12x^2-12x+5)/[16x^2(1-x)^2].                  (DT7)
```

Thus, at every common finite cutoff and regulator,

```text
A_(M,epsilon)=w sum_(m=1)^M w_(m,epsilon)C_m
```

has explicit first and second derivatives obtained from (DT3)--(DT7).
The normalized derivatives are `A_xi=A_x/sqrt(H_B)` and
`A_xixi=A_xx/H_B`.  No numerical differencing of Fresnel tails and no
distributional crossing term is required.

The remaining theorem is quantitative: obtain sum-first bounds for these
finite derivative expressions that are uniform in `M` and `epsilon`, include
the two B-window edge terms, and combine the tangential estimate with the
common-regulator joined remainder below `1.4058e-4`.

Pi provenance: `pi` in (DT1)--(DT5) is inherited from the exact equation-(9)
Fresnel phase.  Equations (DT6)--(DT7) are calculus identities and introduce
no fitted geometric constant.

Proof boundary: exact phase-stripped first/second derivative transport and
finite-sum differentiation only.  This gate does not prove uniform derivative
bounds, a completed exterior-current estimate, a joined-remainder estimate,
complete `Q_K-T` or `T_upper`, or a height-uniform theorem.  It makes no claim
of `Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["endpoint_decomposition"]["decision"]["finite_current_decomposed_exactly"] is True, "endpoint decomposition drift")
    require(dependencies["crossing_completion"]["decision"]["crossing_completed_exterior_B_current_smooth_at_finite_cutoff"] is True, "crossing completion drift")
    require(dependencies["crossing_completion"]["decision"]["completed_exterior_first_second_derivative_bounds_proved"] is False, "derivative-bound scope drift")
    require(dependencies["tangential_cutoff"]["decision"]["exterior_B_trace_uniformly_tangentially_nonstationary"] is True, "tangential phase drift")

    artifact = {
        "kind": STEM,
        "status": "exact_crossing_completed_B_first_second_derivative_transport_reduced",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "decision": {
            "phase_stripped_completed_mode_transport_proved": True,
            "bulk_completion_is_homogeneous_transport_solution": True,
            "exact_first_second_pair_derivatives_proved": True,
            "exact_Kummer_weight_derivatives_proved": True,
            "finite_cutoff_sum_may_be_differentiated_without_crossing_distributions": True,
            "uniform_M_epsilon_derivative_bounds_proved": False,
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
            "process_priority": priority,
        },
        "next_obligation": "Use the exact Y_s transport to derive sum-first enclosures for A, A_xi, and A_xixi on the two completed B exterior intervals, uniform in the finite cutoff and Abel regulator. Split the proof into the finite crossing roster and an analytic outer-mode tail, retain both B-window edge terms, and combine with R_join,comp under one common prescription. The joined target remains below 1.4058e-4.",
        "proof_boundary": "Exact completed-mode derivative transport and finite-sum calculus only. No uniform derivative norm, completed exterior-current or joined-remainder bound, complete Q_K-T or T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(), encoding="utf-8")
    print("certified exact first/second derivative transport for the crossing-completed B current", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
