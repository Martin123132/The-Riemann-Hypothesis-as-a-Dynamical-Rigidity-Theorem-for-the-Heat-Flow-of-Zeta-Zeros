#!/usr/bin/env python3
"""Factor the corrected Airy carrier into exact outgoing Hankel branches."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import mpmath as mp
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_canonical_exact_hankel_branch_factorization_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "beta4_Airy_reduction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_canonical_airy_derivative_reduction_gate.json",
    "ordinary_fold_ownership": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_ordinary_mode_coverage_ledger_gate.json",
}


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


def exact_connection_coefficients() -> dict[str, str]:
    i = sp.I
    sqrt3 = sp.sqrt(3)

    # J_-nu=cos(pi*nu)J_nu-sin(pi*nu)Y_nu and
    # J=(H1+H2)/2, Y=(H1-H2)/(2i).
    def coefficients(nu: sp.Rational, sign: int) -> tuple[sp.Expr, sp.Expr]:
        j_coefficient = sp.cos(sp.pi * nu) + sign
        y_coefficient = -sp.sin(sp.pi * nu)
        h1 = sp.simplify(j_coefficient / 2 + y_coefficient / (2 * i))
        h2 = sp.simplify(j_coefficient / 2 - y_coefficient / (2 * i))
        return h1, h2

    a_h1, a_h2 = coefficients(sp.Rational(1, 3), 1)
    d_h1, d_h2 = coefficients(sp.Rational(2, 3), -1)
    require(sp.simplify(a_h1 - sp.expand_complex(sqrt3 * sp.exp(i * sp.pi / 6) / 2)) == 0, "Ai H1 coefficient drift")
    require(sp.simplify(a_h2 - sp.expand_complex(sqrt3 * sp.exp(-i * sp.pi / 6) / 2)) == 0, "Ai H2 coefficient drift")
    require(sp.simplify(d_h1 + sp.expand_complex(sqrt3 * sp.exp(-i * sp.pi / 6) / 2)) == 0, "A_X H1 coefficient drift")
    require(sp.simplify(d_h2 + sp.expand_complex(sqrt3 * sp.exp(i * sp.pi / 6) / 2)) == 0, "A_X H2 coefficient drift")
    return {
        "J_-1/3_plus_J_1/3_H1": str(a_h1),
        "J_-1/3_plus_J_1/3_H2": str(a_h2),
        "J_-2/3_minus_J_2/3_H1": str(d_h1),
        "J_-2/3_minus_J_2/3_H2": str(d_h2),
    }


def numerical_replay() -> dict[str, Any]:
    mp.mp.dps = 90
    rows = []
    maximum = mp.mpf(0)
    for x in (mp.mpf("0.0003"), mp.mpf("0.17"), mp.mpf("2.75"), mp.mpf("37")):
        xi = mp.mpf(2) * x**mp.mpf("1.5") / 3
        airy = mp.airyai(-x)
        airy_x = -mp.airyai(-x, 1)
        airy_hankel = mp.sqrt(x) / (2 * mp.sqrt(3)) * (
            mp.exp(1j * mp.pi / 6) * mp.hankel1(mp.mpf(1) / 3, xi)
            + mp.exp(-1j * mp.pi / 6) * mp.hankel2(mp.mpf(1) / 3, xi)
        )
        derivative_hankel = -x / (2 * mp.sqrt(3)) * (
            mp.exp(-1j * mp.pi / 6) * mp.hankel1(mp.mpf(2) / 3, xi)
            + mp.exp(1j * mp.pi / 6) * mp.hankel2(mp.mpf(2) / 3, xi)
        )
        error = max(abs(airy - airy_hankel), abs(airy_x - derivative_hankel))
        maximum = max(maximum, error)
        rows.append({"X": mp.nstr(x, 20), "max_absolute_error": mp.nstr(error, 25)})
    require(maximum < mp.mpf("1e-80"), "high-precision Hankel replay failed")
    return {"digits": mp.mp.dps, "rows": rows, "maximum_absolute_error": mp.nstr(maximum, 25)}


def render_note(_: dict[str, Any]) -> str:
    return """# Beta^-4 canonical fold: exact Hankel branch factorization

Date: 2026-08-13
Status: exact branch identity; not a proof of the ordinary-carrier splice

For `X>0`, put

```text
xi=2X^(3/2)/3,       A(X)=Ai(-X),       A_X=dA/dX.
```

The exact Bessel connection formulae give

```text
A(X)=sqrt(X)/(2sqrt(3))
 {e^(i*pi/6) H^(1)_(1/3)(xi)
  +e^(-i*pi/6)H^(2)_(1/3)(xi)},                       (HB1)

A_X(X)=-X/(2sqrt(3))
 {e^(-i*pi/6)H^(1)_(2/3)(xi)
  +e^(i*pi/6) H^(2)_(2/3)(xi)}.                       (HB2)
```

The `pi/6` phases are forced by the Bessel reflection identity at orders
`1/3` and `2/3`; they are not fitted phases.  If `U,V` are the real rational
weights from the beta-minus-four Airy reduction, define

```text
C_1(X)=1/(2sqrt(3)){
 U sqrt(X)e^(i*pi/6)H^(1)_(1/3)(xi)
 -V X e^(-i*pi/6)H^(1)_(2/3)(xi)},

C_2(X)=1/(2sqrt(3)){
 U sqrt(X)e^(-i*pi/6)H^(2)_(1/3)(xi)
 -V X e^(i*pi/6)H^(2)_(2/3)(xi)}.                    (HB3)
```

Then, exactly,

```text
U A+V A_X=C_1+C_2,             C_2=conj(C_1)          (HB4)
```

for real `X,U,V`.  The two Hankel functions are independent because

```text
W(H^(1)_nu,H^(2)_nu)=-4i/(pi*xi).                     (HB5)
```

Their standard leading forms identify the phase contract, without yet
asserting a quantitative remainder:

```text
C_1 ~[U X^(-1/4)+i V X^(1/4)]/(2sqrt(pi))
      *exp(i[xi-pi/4]),
C_2 ~[U X^(-1/4)-i V X^(1/4)]/(2sqrt(pi))
      *exp(-i[xi-pi/4]).                              (HB6)
```

After restoring the normal factor `exp(i[y^2/(4beta)-d*y])`, the two phase
derivatives are therefore

```text
y/(2beta)-d+sqrt(lambda+y),
y/(2beta)-d-sqrt(lambda+y).                           (HB7)
```

These are exactly the two ordinary-saddle equations already visible in the
cubic fold geometry.  Unlike a large-`X` WKB replacement, (HB1)--(HB4) remain
exact for every `X>0`.  The raw Hankel functions diverge at `xi=0`, but the
`sqrt(X)` and `X` prefactors in (HB3) give finite one-sided branch limits.  At
`X=0` the written products need either the regular Airy basis or separately
certified limiting values.  This distinction matters for a stable fold collar.

The remaining obligation is to choose and certify that collar, match the
appropriate exact Hankel branch to the endpoint-retaining logistic/Gamma
carrier on each ordinary corridor, and bound the branch remainder and
overlaps.  This gate proves no such quantitative splice, no `Q_K-T` bound,
no complete `T_upper`, no height-uniform theorem, no `Lambda<=0`, no
PF-infinity, no RH, and no prize-level conclusion.
"""


def main() -> int:
    started = time.perf_counter()
    priority = set_low_priority()
    require(CHECKER.is_file(), "missing independent checker")
    dependencies = {}
    for name, path in DEPENDENCIES.items():
        require(path.is_file(), f"missing dependency: {name}")
        dependencies[name] = json.loads(path.read_text(encoding="utf-8"))
    require(dependencies["beta4_Airy_reduction"]["decision"]["all_derivatives_reduced_to_Ai_and_first_derivative"] is True, "Airy reduction drift")
    require(dependencies["ordinary_fold_ownership"]["decision"]["quantitative_splice_proved"] is False, "ownership boundary drift")

    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "exact_beta_minus_4_Airy_carrier_Hankel_branch_factorization_proved",
        "passed": True,
        "scope": {"X": "positive real", "U_V": "real corrected Airy weights"},
        "connection_coefficients": exact_connection_coefficients(),
        "exact_identities": {
            "A": "sqrt(X)/(2*sqrt(3))*(exp(i*pi/6)*H1_(1/3)(xi)+exp(-i*pi/6)*H2_(1/3)(xi))",
            "A_X": "-X/(2*sqrt(3))*(exp(-i*pi/6)*H1_(2/3)(xi)+exp(i*pi/6)*H2_(2/3)(xi))",
            "xi": "2*X^(3/2)/3",
            "branch_sum": "U*A+V*A_X=C1+C2",
            "real_conjugation": "C2=conjugate(C1) for real X,U,V",
            "Hankel_Wronskian": "W(H1_nu,H2_nu)=-4*i/(pi*xi)",
        },
        "leading_phase_contract": {
            "C1": "(U*X^(-1/4)+i*V*X^(1/4))/(2*sqrt(pi))*exp(i*(xi-pi/4))",
            "C2": "(U*X^(-1/4)-i*V*X^(1/4))/(2*sqrt(pi))*exp(-i*(xi-pi/4))",
            "normal_phase_derivatives": ["y/(2*beta)-d+sqrt(lambda+y)", "y/(2*beta)-d-sqrt(lambda+y)"],
            "quantitative_remainder_proved": False,
        },
        "numerical_replay": numerical_replay(),
        "decision": {
            "exact_Hankel_branch_factorization_proved": True,
            "branch_basis_nonsingular_for_positive_X": True,
            "raw_Hankel_functions_singular_at_X_zero": True,
            "prefactored_branches_have_finite_positive_limits": True,
            "Airy_or_explicit_limit_required_at_X_zero": True,
            "ordinary_logistic_branch_identification_proved": False,
            "corridor_remainder_bound_proved": False,
            "complete_Q_K_minus_T_bound_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "next_action": "Choose a rigorous Airy collar and compare the exact outgoing Hankel branch with the endpoint-retaining logistic/Gamma carrier on each owned ordinary corridor, including overlap remainders.",
        "proof_boundary": "Exact positive-X branch factorization and leading phase contract only. No quantitative Hankel asymptotic remainder, ordinary-carrier identification, corridor splice, Q_K-T bound, complete T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
        "dependencies": {name: {"path": relative(path), "sha256": file_hash(path)} for name, path in DEPENDENCIES.items()},
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "resource_policy": {"workers": 1, "process_priority": priority},
        "runtime": {"elapsed_seconds": round(time.perf_counter() - started, 3)},
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    artifact["sources"]["note"] = {"path": relative(NOTE), "sha256": file_hash(NOTE)}
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print("proved exact positive-X Hankel branch factorization: two conjugate branches, 0 corridor-splice claims", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
