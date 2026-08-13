#!/usr/bin/env python3
"""Derive and bound the beta^-6 common-profile coefficient."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta6_common_profile_leading_coefficient_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "common_profile_operator": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_finite_t_common_profile_operator_gate.json",
    "beta4_reduction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_canonical_airy_derivative_reduction_gate.json",
    "weighted_kernel": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_opposite_branch_dirichlet_gate.json",
}

C = 159_577
PANELS = 4096
PRECISION = 80


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
    q, z, y, x, lam = sp.symbols("q z y X lambda", real=True)
    i = sp.I
    amplitude = sp.cosh(q * z) ** (-sp.Rational(3, 2)) * (1 + q**2 * y / 2)
    phase_remainder = (
        (q * z - sp.tanh(q * z)) / q**3
        - z**3 / 3
        + y * (z - sp.tanh(q * z) / q)
        - q * y**2 * sp.tanh(q * z) / 4
    )
    ratio = sp.series(amplitude * sp.exp(i * phase_remainder), q, 0, 8).removeO().expand()
    c3 = sp.expand(ratio.coeff(q, 6))

    a3 = z**4 * (390 * y - 379 * z**2) / 1920
    r3 = -z**5 * (189 * y**2 - 306 * y * z**2 + 124 * z**4) / 5670
    a1 = y / 2 - 3 * z**2 / 4
    r1 = -2 * z**5 / 15 + z**3 * y / 3 - z * y**2 / 4
    a2 = 13 * z**4 / 32 - 3 * z**2 * y / 8
    r2 = 17 * z**7 / 315 - 2 * z**5 * y / 15 + z**3 * y**2 / 12
    expected_c3 = sp.expand(
        a3
        + i * r3
        - r1 * r2
        - i * r1**3 / 6
        + a1 * (i * r2 - r1**2 / 2)
        + i * a2 * r1
    )
    require(sp.expand(c3 - expected_c3) == 0, "beta^-6 multiplier coefficient drift")
    require(sp.degree(c3, z) == 15, "beta^-6 multiplier degree drift")

    p, r = sp.Integer(1), sp.Integer(0)
    u3, v3 = sp.Integer(0), sp.Integer(0)
    for degree in range(16):
        coefficient = c3.coeff(z, degree)
        u3 += coefficient * i**degree * p
        v3 += coefficient * i**degree * r
        p, r = sp.expand(sp.diff(p, x) - x * r), sp.expand(p + sp.diff(r, x))
    u3 = sp.factor(sp.simplify(u3.subs(x, lam + y)))
    v3 = sp.factor(sp.simplify(v3.subs(x, lam + y)))

    expected_u3 = (
        -49856 * lam**6 + 1344 * lam**5 * y + 3360 * lam**4 * y**2
        - 12920 * lam**3 * y**3 - 516535 * lam**3 + 7785 * lam**2 * y**4
        + 5445 * lam**2 * y - 1836 * lam * y**5 + 21195 * lam * y**2
        + 1674 * y**6 - 21285 * y**3 - 265860
    ) / 9072000
    expected_v3 = -(
        3584 * lam**7 - 1792 * lam**6 * y + 1344 * lam**5 * y**2
        + 2240 * lam**4 * y**3 - 43880 * lam**4 - 1960 * lam**3 * y**4
        + 18580 * lam**3 * y + 1764 * lam**2 * y**5 + 5445 * lam**2 * y**2
        - 567 * lam * y**6 + 14130 * lam * y**3 - 127530 * lam
        + 189 * y**7 - 7605 * y**4 + 42570 * y
    ) / 9072000
    require(sp.expand(u3 - expected_u3) == 0, "U3 contraction drift")
    require(sp.expand(v3 - expected_v3) == 0, "V3 contraction drift")
    return {
        "a3": a3,
        "r3": r3,
        "C3": c3,
        "U3": u3,
        "V3": v3,
    }


def interval_uv3(lam: arb, y: arb) -> tuple[arb, arb]:
    u3 = (
        -49856 * lam**6 + 1344 * lam**5 * y + 3360 * lam**4 * y**2
        - 12920 * lam**3 * y**3 - 516535 * lam**3 + 7785 * lam**2 * y**4
        + 5445 * lam**2 * y - 1836 * lam * y**5 + 21195 * lam * y**2
        + 1674 * y**6 - 21285 * y**3 - 265860
    ) / 9072000
    v3 = -(
        3584 * lam**7 - 1792 * lam**6 * y + 1344 * lam**5 * y**2
        + 2240 * lam**4 * y**3 - 43880 * lam**4 - 1960 * lam**3 * y**4
        + 18580 * lam**3 * y + 1764 * lam**2 * y**5 + 5445 * lam**2 * y**2
        - 567 * lam * y**6 + 14130 * lam * y**3 - 127530 * lam
        + 189 * y**7 - 7605 * y**4 + 42570 * y
    ) / 9072000
    return u3, v3


def interval_certificate(kernel: dict[str, Any]) -> dict[str, Any]:
    ctx.dps = PRECISION
    pi = arb.pi()
    c = arb(C)
    beta = (pi * c**2 / 8) ** (arb(1) / 3)
    epsilon = beta**-2
    y_width = pi * c / (2 * beta)
    lambda_max = pi / (16 * beta)
    h = 4 * beta / c
    lam = arb(lambda_max / 2, lambda_max / 2)
    maximum = arb(0)
    maximum_panel = -1

    for panel in range(PANELS):
        midpoint = y_width * (2 * panel + 1) / (2 * PANELS)
        radius = y_width / (2 * PANELS)
        y = arb(midpoint, radius)
        u3, v3 = interval_uv3(lam, y)
        x = lam + y
        airy = (-x).airy_ai()
        airy_x = (-x).airy_ai(derivative=1)
        bound = 2 * pi * epsilon**3 * (
            abs(u3).upper() * abs(airy).upper()
            + abs(v3).upper() * abs(airy_x).upper()
        )
        if bound.upper() > maximum.upper():
            maximum = bound
            maximum_panel = panel

    kernel_l1 = arb(kernel["certificate"]["two_branch_kernel_L1_bound_ball"])
    profile_bound = maximum / h
    weighted_bound = profile_bound * kernel_l1
    require(maximum < arb("6.951e-10"), "beta^-6 reduced coefficient exceeds 6.951e-10")
    require(profile_bound < arb("1.288e-8"), "beta^-6 canonical profile exceeds 1.288e-8")
    require(weighted_bound < arb("4.95e-13"), "beta^-6 weighted transfer exceeds 4.95e-13")
    return {
        "panels": PANELS,
        "precision_decimal_digits": PRECISION,
        "height_interval": "t*-pi/16<=t<=t*",
        "lambda_interval": "0<=lambda<=pi/(16beta)",
        "y_interval": "0<=y<=Y=pi*C/(2beta)",
        "maximum_panel": maximum_panel,
        "maximum_panel_y_midpoint_ball": (y_width * (2 * maximum_panel + 1) / (2 * PANELS)).str(PRECISION, more=True),
        "beta_ball": beta.str(PRECISION, more=True),
        "epsilon_beta_minus_2_ball": epsilon.str(PRECISION, more=True),
        "Y_ball": y_width.str(PRECISION, more=True),
        "lambda_max_ball": lambda_max.str(PRECISION, more=True),
        "h_ball": h.str(PRECISION, more=True),
        "uniform_reduced_beta6_coefficient_bound_ball": maximum.str(PRECISION, more=True),
        "uniform_canonical_beta6_profile_bound_ball": profile_bound.str(PRECISION, more=True),
        "weighted_kernel_L1_bound_ball": kernel["certificate"]["two_branch_kernel_L1_bound_ball"],
        "uniform_weighted_beta6_coefficient_bound_ball": weighted_bound.str(PRECISION, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    p = artifact["polynomial_certificate"]
    c = artifact["interval_certificate"]
    return f"""# Beta^-6 common-profile leading coefficient

Date: 2026-08-13

Status: exact formal coefficient and rigorous whole-strip coefficient bound;
the exact post-beta-minus-six Taylor remainder remains open

Put `q=beta^-1`, `X=lambda+y`, and factor the exact transformed integrand by
the cubic carrier `exp(i[z^3/3-Xz])`.  Its exact ratio has the expansion

```text
A(q,z,y)exp(iR(q,z,y))
 =1+q^2 C1(z,y)+q^4 C2(z,y)+q^6 C3(z,y)+O(q^8).       (B6.1)
```

The third amplitude and phase coefficients are

```text
a3={p['a3']},
r3={p['r3']}.                                         (B6.2)
```

Exact multiplication gives

```text
C3=a3+i r3-r1 r2-i r1^3/6+a1(i r2-r1^2/2)+i a2 r1.  (B6.3)
```

`C3` has degree fifteen in `z`.  Contracting every monomial with

```text
integral_R^Abel z^k e^(i[z^3/3-Xz])dz
 =2pi i^k d_X^k Ai(-X)                               (B6.4)
```

and reducing by `A_XX=-X A` yields exactly

```text
J6(lambda,y)=2pi beta^-6
 [U3(lambda,y)Ai(-lambda-y)
  +V3(lambda,y)d_X Ai(-lambda-y)],                    (B6.5)

U3={p['U3']},
V3={p['V3']}.                                         (B6.6)
```

A `4096`-panel Arb interval atlas encloses the complete top corridor
`0<=lambda<=pi/(16beta)` and `0<=y<=Y`.  It proves

```text
sup |J6(lambda,y)|
 <{c['uniform_reduced_beta6_coefficient_bound_ball']}<6.951e-10,

sup |g6(lambda,s)|=sup |J6|/h
 <{c['uniform_canonical_beta6_profile_bound_ball']}<1.288e-8.           (B6.7)
```

Passing this coefficient through the already certified whole weighted
kernel gives

```text
|Delta W6|
 <{c['uniform_weighted_beta6_coefficient_bound_ball']}<4.95e-13.       (B6.8)
```

Equations (B6.7)--(B6.8) rigorously bound the coefficient of the first
omitted asymptotic order; they are not a bound on the entire exact-minus-
beta-minus-four profile until the `O(q^8)` remainder in (B6.1) is enclosed.
No extrapolation from coefficient size is promoted.

Pi provenance: `2pi` in (B6.4) is the inverse Airy Fourier normalization;
`beta^3=pi C^2/8`, `Y=pi C/(2beta)`, and `hY=2pi` are inherited from the
exact Kummer/Fourier selector geometry.  No fitted or geometric constant is
introduced.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No exact Taylor-tail estimate, full exact-minus-beta-minus-four profile
bound, weighted finite-`t` remainder, complete source/initial-data splice,
complete `Q_K-T` or `T_upper`, all-corridor or height-uniform theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    for name, dependency in dependencies.items():
        require(dependency.get("passed") is True, f"dependency failed: {name}")

    polynomials = polynomial_data()
    interval = interval_certificate(dependencies["weighted_kernel"])
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "beta6_common_profile_leading_coefficient_uniformly_bounded",
        "passed": True,
        "polynomial_certificate": {name: str(value) for name, value in polynomials.items()},
        "interval_certificate": interval,
        "decision": {
            "beta6_multiplier_coefficient_exactly_derived": True,
            "beta6_Airy_two_channel_contraction_exact": True,
            "beta6_coefficient_uniformly_bounded_on_top_corridor": True,
            "weighted_beta6_coefficient_bound_proved": True,
            "exact_post_beta6_Taylor_remainder_bound_proved": False,
            "full_exact_minus_beta4_profile_bound_proved": False,
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
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
        "next_action": (
            "Enclose the exact post-beta^-6 ratio remainder uniformly on a common contour, or derive beta^-8 and a rigorous "
            "analytic tail majorant. Do not infer the full remainder from the small beta^-6 coefficient alone."
        ),
        "proof_boundary": (
            "Exact beta^-6 coefficient algebra and a rigorous uniform coefficient bound only. No exact Taylor-tail estimate, "
            "full exact-minus-beta^-4 profile bound, weighted finite-t remainder, complete source/initial-data splice, complete "
            "Q_K-T or T_upper, all-corridor or height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("certified beta^-6 common-profile leading coefficient", flush=True)


if __name__ == "__main__":
    main()
