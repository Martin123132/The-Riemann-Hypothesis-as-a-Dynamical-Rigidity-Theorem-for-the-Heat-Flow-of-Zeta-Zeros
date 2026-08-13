#!/usr/bin/env python3
"""Reduce the beta^-4 corrected fold canonical model to Airy derivatives."""

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

import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_canonical_airy_derivative_reduction_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "leading_Airy_reduction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_canonical_core_quadrature_gate.json",
    "beta4_normal_form": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_extreme_second_order_normal_form_gate.json",
    "all_event_cells": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_second_order_all_event_continuous_height_gate.json",
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


def polynomial_data() -> dict[str, Any]:
    z, y, x, lam, eps = sp.symbols("z y X lambda epsilon", real=True)
    i = sp.I
    a1 = y / 2 - 3 * z**2 / 4
    r1 = -2 * z**5 / 15 + z**3 * y / 3 - z * y**2 / 4
    a2 = 13 * z**4 / 32 - 3 * z**2 * y / 8
    r2 = 17 * z**7 / 315 - 2 * z**5 * y / 15 + z**3 * y**2 / 12
    c1 = sp.expand(a1 + i * r1)
    c2 = sp.expand(a2 + i * r2 + i * a1 * r1 - r1**2 / 2)
    multiplier = sp.expand(1 + eps * c1 + eps**2 * c2)
    require(sp.degree(c1, z) == 5, "beta^-2 z degree drift")
    require(sp.degree(c2, z) == 10, "beta^-4 z degree drift")

    p = [sp.Integer(1)]
    q = [sp.Integer(0)]
    for _ in range(10):
        old_p, old_q = p[-1], q[-1]
        p.append(sp.expand(sp.diff(old_p, x) - x * old_q))
        q.append(sp.expand(old_p + sp.diff(old_q, x)))

    u = sp.Integer(0)
    v = sp.Integer(0)
    for k in range(11):
        coefficient = sp.expand(multiplier).coeff(z, k)
        u += coefficient * i**k * p[k]
        v += coefficient * i**k * q[k]
    u = sp.collect(sp.expand(u), eps)
    v = sp.collect(sp.expand(v), eps)

    u1 = -(13 * x - 10 * y) / 60
    v1 = (8 * x**2 - 20 * x * y + 15 * y**2) / 60
    u2 = -(
        448 * x**5 - 2240 * x**4 * y + 4480 * x**3 * y**2
        - 4200 * x**2 * y**3 + 4565 * x**2 + 1575 * x * y**4
        - 9100 * x * y + 4130 * y**2
    ) / 50400
    v2 = (40 * x**3 - 140 * x**2 * y + 161 * x * y**2 - 70 * y**3 + 27) / 1680
    require(sp.expand(u - (1 + eps * u1 + eps**2 * u2)) == 0, "U reduction failed")
    require(sp.expand(v - (eps * v1 + eps**2 * v2)) == 0, "V reduction failed")

    substitution = {x: lam + y}
    u1_lam = sp.factor(u1.subs(substitution))
    v1_lam = sp.factor(v1.subs(substitution))
    u2_lam = sp.factor(u2.subs(substitution))
    v2_lam = sp.factor(v2.subs(substitution))
    require(sp.expand(u1_lam + (13 * lam + 3 * y) / 60) == 0, "lambda U1 simplification failed")
    require(sp.expand(v1_lam - (8 * lam**2 - 4 * lam * y + 3 * y**2) / 60) == 0, "lambda V1 simplification failed")

    return {
        "symbols": {"z": z, "y": y, "X": x, "lambda": lam, "epsilon": eps},
        "normal_form": {"a1": a1, "r1": r1, "a2": a2, "r2": r2, "C1": c1, "C2": c2},
        "multiplier": multiplier,
        "recurrence": list(zip(p, q)),
        "weights_X": {"U1": u1, "V1": v1, "U2": u2, "V2": v2},
        "weights_lambda": {"U1": u1_lam, "V1": v1_lam, "U2": u2_lam, "V2": v2_lam},
    }


def render_note(artifact: dict[str, Any]) -> str:
    w = artifact["reduced_weights_lambda_y"]
    return f"""# Beta^-4 canonical fold: exact Airy-derivative reduction

Date: 2026-08-13
Status: exact Abel/Fourier reduction; not a proof of the ordinary-carrier bound

Put `X=lambda+y`, `epsilon=beta^-2`, and

```text
A(X)=Ai(-X),                    A_X=dA/dX,
D_M(y)=sum_(m in M) exp(-i*d_m*y).
```

For the corrected event normal form use

```text
M_4(z,y)=1+epsilon(a_1+i*r_1)
 +epsilon^2(a_2+i*r_2+i*a_1*r_1-r_1^2/2),            (AR1)

a_1=y/2-3z^2/4,
r_1=-2z^5/15+z^3y/3-zy^2/4,
a_2=13z^4/32-3z^2y/8,
r_2=17z^7/315-2z^5y/15+z^3y^2/12.
```

The real-line `z` integrals are symmetric Abel limits.  The Fourier
normalization of the Airy function gives, for `0<=k<=10`,

```text
integral_R^Abel z^k exp(i[z^3/3-Xz])dz
   =2*pi*i^k d^k A(X)/dX^k.                            (AR2)
```

This is the provenance of `pi`: it is exactly the `2*pi` in the Airy Fourier
transform, not a fitted constant or an arbitrary circle.  Since `A_XX=-X*A`,
write `A^(k)=P_k(X)A+Q_k(X)A_X`, where

```text
P_0=1, Q_0=0,
P_(k+1)=P_k'-X*Q_k,       Q_(k+1)=P_k+Q_k'.            (AR3)
```

Substitution through degree ten collapses the complete corrected two-variable
model to

```text
I_M^[4](beta,lambda,Y)
 =2*pi integral_0^Y exp(i*y^2/(4beta))*D_M(y)
    *[U(lambda,y,beta)A(lambda+y)
      +V(lambda,y,beta)A_X(lambda+y)]dy,               (AR4)

U=1+epsilon*U_1+epsilon^2*U_2,
V=  epsilon*V_1+epsilon^2*V_2,

U_1={w['U1']},
V_1={w['V1']},
U_2={w['U2']},
V_2={w['V2']}.                                        (AR5)
```

Thus the beta^-4 fold model has only two scalar carrier channels, not an
unresolved two-dimensional oscillatory integral.  The leading formula in
Section 11.315 is recovered exactly by setting `epsilon=0`.

The remaining wall is now explicit.  One must prove an Airy branch expansion
with remainder in the ordinary corridors, match both phase branches and their
normalizations to the classical/Gamma carrier, retain the exact Airy form near
`X=0`, and prove overlap ownership without omission or duplication.  This gate
does not perform that asymptotic match and proves no bound for `Q_K-T`, no
complete `T_upper`, no height-uniform theorem, no `Lambda<=0`, no PF-infinity,
no RH, and no prize-level conclusion.
"""


def main() -> int:
    started = time.perf_counter()
    priority = set_low_priority()
    require(CHECKER.is_file(), "missing independent checker")
    dependencies = {}
    for name, path in DEPENDENCIES.items():
        require(path.is_file(), f"missing dependency: {name}")
        dependencies[name] = json.loads(path.read_text(encoding="utf-8"))
    require(dependencies["leading_Airy_reduction"]["decision"]["two_one_dimensional_representations_proved"] is True, "leading Airy reduction drift")
    require(dependencies["beta4_normal_form"]["claims"]["complete_beta_minus_4_multiplier_derived"] is True, "beta^-4 multiplier drift")
    require(dependencies["all_event_cells"]["claims"]["all_399_exact_event_cells_uniformly_certified"] is True, "event atlas drift")

    data = polynomial_data()
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "exact_beta_minus_4_fold_canonical_Airy_derivative_reduction_proved",
        "passed": True,
        "scope": {
            "beta": "positive real",
            "lambda": "real",
            "Y": "finite nonnegative",
            "mode_roster": "arbitrary finite real detuning roster",
            "z_integral": "symmetric Abel oscillatory limit on the real line",
        },
        "normal_form": {name: str(value) for name, value in data["normal_form"].items()},
        "moment_identity": "Integral_Abel z^k exp(i*(z^3/3-X*z)) dz = 2*pi*i^k*d_X^k Ai(-X), 0<=k<=10",
        "Airy_ODE": "A(X)=Ai(-X); A''(X)=-X*A(X)",
        "derivative_recurrence": [
            {"k": k, "P_k": str(p), "Q_k": str(q)}
            for k, (p, q) in enumerate(data["recurrence"])
        ],
        "reduced_weights_X_y": {name: str(value) for name, value in data["weights_X"].items()},
        "reduced_weights_lambda_y": {name: str(value) for name, value in data["weights_lambda"].items()},
        "exact_reduction": "2*pi*Integral_0^Y exp(i*y^2/(4*beta))*D_M(y)*(U*Ai(-lambda-y)+V*d_X Ai(-X)|_(X=lambda+y))dy",
        "pi_provenance": "Ai(-X)=(2*pi)^(-1) Integral_Abel exp(i*(z^3/3-X*z))dz",
        "decision": {
            "beta_minus_4_two_dimensional_model_reduced_exactly": True,
            "highest_required_Airy_derivative": 10,
            "all_derivatives_reduced_to_Ai_and_first_derivative": True,
            "leading_Section_11_315_formula_recovered_at_epsilon_zero": True,
            "ordinary_carrier_branch_match_proved": False,
            "complete_Q_K_minus_T_bound_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "next_action": "Prove a remainder-controlled Ai(-X), d_X Ai(-X) branch expansion on each ordinary corridor and identify it with the classical/Gamma carrier, while retaining the exact Airy form at fold overlap.",
        "proof_boundary": "Exact corrected canonical-model reduction only. No Airy-to-ordinary carrier remainder theorem, corridor splice, Q_K-T bound, complete T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
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
    print(
        "proved exact beta^-4 Airy reduction: z degree 10 -> Ai(-X) and d_X Ai(-X), 0 ordinary-carrier claims",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
