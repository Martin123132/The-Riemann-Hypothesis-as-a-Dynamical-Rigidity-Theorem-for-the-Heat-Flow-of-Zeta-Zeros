#!/usr/bin/env python3
"""Compress the finite physical Q_K roster to one geometric half-line integral."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import mpmath as mp

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_QK_geometric_half_line_and_saddle_route_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "branch_alignment": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_label_branch_correction_QK_alignment_gate.json",
}

HEIGHT = 10_000_000_000
SOURCE_ALPHA_MIN = 159_577
SOURCE_ALPHA_MAX = 5_122_421
SOURCE_COUNT = 2_481_423
PRECISION = 100
SURROGATE_T = mp.mpf("5")
SURROGATE_ALPHA_MIN = 1
SURROGATE_ALPHA_MAX = 7


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {relative(path)}")
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


def source_label(t: mp.mpf, alpha: int) -> mp.mpc:
    pi = mp.pi
    a = mp.mpf("0.25") - 0.5j * t
    b = mp.mpf("0.75") - 0.5j * t
    xi = 1j * pi * alpha * alpha / 4
    c = pi * alpha * (1 + 1j) / mp.sqrt(2)
    A = mp.gamma(a) * pi ** (-a) * mp.hyp1f1(a, mp.mpf("0.5"), xi) / 2
    B = c * mp.gamma(b) * pi ** (-b) * mp.hyp1f1(b, mp.mpf("1.5"), xi) / 2
    return mp.e ** (3 * pi * t / 4 + 3j * pi / 8) * (A - B)


def geometric_factor(x: mp.mpf, alpha_min: int, count: int) -> mp.mpc:
    if x == 0:
        return mp.mpf(count)
    lam = mp.pi * (1 + 1j) / mp.sqrt(2)
    return (
        mp.e ** (-lam * alpha_min * x)
        * mp.expm1(-2 * lam * count * x)
        / mp.expm1(-2 * lam * x)
    )


def compressed_integral(t: mp.mpf, alpha_min: int, alpha_max: int) -> mp.mpc:
    count = (alpha_max - alpha_min) // 2 + 1
    s = mp.mpf("0.5") + 1j * t

    def integrand(x: mp.mpf) -> mp.mpc:
        return x ** (-s) * mp.e ** (-mp.pi * x * x) * geometric_factor(x, alpha_min, count)

    integral = mp.quad(integrand, [0, mp.mpf("0.01"), mp.mpf("0.1"), 1, 4, mp.inf])
    return mp.e ** (3 * mp.pi * t / 4 + 3j * mp.pi / 8) * integral


def saddle_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    t = arb(HEIGHT)
    threshold = (8 * t / pi).sqrt()
    require(arb(SOURCE_ALPHA_MIN) > threshold, "lower source label is not above the saddle threshold")
    require(threshold > arb(SOURCE_ALPHA_MIN - 2), "lower source label is not the first odd above threshold")

    def row(alpha_value: int) -> dict[str, Any]:
        alpha = arb(alpha_value)
        discriminant = acb(-4 * pi, pi * (pi * alpha * alpha - 8 * t))
        root = discriminant.sqrt()
        lam_alpha = acb(pi * alpha / arb(2).sqrt(), pi * alpha / arb(2).sqrt())
        saddle_plus = (-lam_alpha + root) / (4 * pi)
        saddle_minus = (-lam_alpha - root) / (4 * pi)
        require(saddle_plus.real.upper() < 0 and saddle_plus.imag.upper() < 0, "plus saddle left quadrant III")
        require(saddle_minus.real.upper() < 0 and saddle_minus.imag.upper() < 0, "minus saddle left quadrant III")
        return {
            "alpha": alpha_value,
            "discriminant_real_ball": discriminant.real.str(45, more=True),
            "discriminant_imaginary_ball": discriminant.imag.str(55, more=True),
            "saddle_plus_ball": saddle_plus.str(55, more=True),
            "saddle_minus_ball": saddle_minus.str(55, more=True),
            "both_saddles_in_open_quadrant_III": True,
        }

    return {
        "transition_alpha_ball": threshold.str(70, more=True),
        "previous_odd_label": SOURCE_ALPHA_MIN - 2,
        "first_odd_label_above_transition": SOURCE_ALPHA_MIN,
        "first_label_minus_transition_ball": (arb(SOURCE_ALPHA_MIN) - threshold).str(55, more=True),
        "rows": [row(SOURCE_ALPHA_MIN - 2), row(SOURCE_ALPHA_MIN), row(SOURCE_ALPHA_MIN + 2), row(SOURCE_ALPHA_MAX)],
    }


def render_note(artifact: dict[str, Any]) -> str:
    saddle = artifact["saddle_geometry"]
    surrogate = artifact["surrogate_compression_check"]
    return f"""# Finite Q_K geometric half-line compression and saddle route

Date: 2026-08-27

Status: exact finite compression certified; saddle-ray deformation open

Section 11.484 writes every explicit A33 label as one convergent half-line
integral.  Because the physical roster is finite, those integrals can be
summed without any interchange theorem.  Put

```text
s=1/2+it,
lambda=pi(1+i)/sqrt(2),
K_0=exp(3pi*t/4+i3pi/8),
M=(B-A)/2+1=2481423.                                (GH1)
```

For one odd label,

```text
S_alpha=K_0 integral_0^infinity
        x^(-s)exp(-pi*x^2-lambda*alpha*x)dx.         (GH2)
```

Therefore, exactly,

```text
S_W=K_0 integral_0^infinity x^(-s)exp(-pi*x^2)G_W(x)dx,

G_W(x)=sum_(j=0)^(M-1)exp[-lambda(A+2j)x]
      =exp(-lambda*A*x)
       (1-exp(-2lambda*M*x))/(1-exp(-2lambda*x)),    (GH3)

G_W(0)=M.                                           (GH4)
```

For real `x>0`, `|exp(-2lambda*x)|<1`, so the denominator in (GH3) is
nonzero.  The apparent endpoint singularity is removable by (GH4).  With
`Hardy_t[X]=2Re[exp(i theta(t))X]`, the corrected physical roster is

```text
H(t)Q_K=Hardy_t[S_W].                               (GH5)
```

Thus 2,481,423 ill-conditioned direct Kummer calls have been replaced by one
exact scalar integral and a stable finite geometric quotient.  At the
surrogate `t={surrogate['t']}`, roster `{surrogate['alpha_min']}..{surrogate['alpha_max']}`, direct summation of the
individual Kummer labels and independent quadrature of (GH3) differ by only
`{surrogate['absolute_discrepancy']}`.

The one-label exponent in (GH2) is

```text
F_alpha(x)=-pi*x^2-lambda*alpha*x-s Log(x).          (GH6)
```

Its saddles solve

```text
2pi*x^2+lambda*alpha*x+s=0,
D_alpha=lambda^2 alpha^2-8pi*s
       =-4pi+i*pi*(pi*alpha^2-8t).                  (GH7)
```

The transition is therefore exactly

```text
alpha_0=sqrt(8t/pi).
```

At `t=10^10`, Arb gives

```text
alpha_0={saddle['transition_alpha_ball']},
159577-alpha_0={saddle['first_label_minus_transition_ball']}. (GH8)
```

Hence `159577` is the first odd label above the saddle transition; `159575`
is below it.  Arb also places both roots of (GH7) in the open third quadrant
at the previous odd label, the first two roster labels, and the upper roster
endpoint.

There is a revealing algebraic candidate contour.  On the ray
`x=exp(-3i*pi/4)y`, using `Log(x)=log(y)-3i*pi/4`, the whole large prefactor
cancels pointwise:

```text
K_0 x^(-s)exp(-pi*x^2-lambda*alpha*x)dx
 =y^(-s)exp[i*pi*(alpha*y-y^2)]dy.                  (GH9)
```

The real stationary equation of the right-hand side is

```text
2pi*y^2-pi*alpha*y+t=0,                             (GH10)
```

with the same threshold `alpha_0`.  Equation (GH9) explains the source's
transition geometry and removes the exponentially large normalization at
the integrand level.  It is not yet permission to rotate the integral:
`exp(-pi*x^2)` grows in intervening Stokes sectors, and the branch at zero
must be connected with the correct multiplier.  A valid next proof must
derive that connection by a finite steepest-descent contour or an exact
parabolic-cylinder connection formula, including every large-arc and branch
contribution, before using (GH9) as an integral identity.

Pi provenance: all occurrences come from the RSI Gaussian/sine kernel,
quarter-turn contour, and gamma/Kummer normalization.  The transition value
`sqrt(8t/pi)` is derived from the discriminant (GH7), not inserted or fitted.

Proof boundary: exact finite geometric compression, removable endpoint,
surrogate quadrature check, exact saddle polynomial, and actual threshold
arithmetic only.  The saddle-ray integrand identity is algebraic, but its
contour deformation and Stokes multiplier remain open.  No actual-height
quadrature, `B_W`, `Q_K`, `D_K`, `Delta_KU`, `J_Z`, non-A, all-height,
`Lambda<=0`, PF-infinity, RH, or prize-level enclosure is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "dependency failure")
    branch = dependencies["branch_alignment"]
    require(
        branch["decision"]["physical_QK_prefactor_and_branch_alignment_complete_with_correction"] is True,
        "branch-corrected Q_K alignment drift",
    )
    window = branch["actual_window_augmented_identity"]
    require(window["alpha_min"] == SOURCE_ALPHA_MIN, "source lower endpoint drift")
    require(window["alpha_max"] == SOURCE_ALPHA_MAX, "source upper endpoint drift")
    require(window["window_label_count"] == SOURCE_COUNT, "source count drift")

    mp.mp.dps = 65
    direct = mp.fsum(
        source_label(SURROGATE_T, alpha)
        for alpha in range(SURROGATE_ALPHA_MIN, SURROGATE_ALPHA_MAX + 1, 2)
    )
    compressed = compressed_integral(SURROGATE_T, SURROGATE_ALPHA_MIN, SURROGATE_ALPHA_MAX)
    discrepancy = abs(direct - compressed)
    require(discrepancy < mp.mpf("1e-27"), "surrogate geometric compression failed")

    saddle = saddle_certificate()
    artifact = {
        "kind": STEM,
        "status": "exact_finite_QK_geometric_half_line_compression_and_guarded_saddle_route_certified",
        "passed": True,
        "scope": {
            "height": HEIGHT,
            "alpha_min": SOURCE_ALPHA_MIN,
            "alpha_max": SOURCE_ALPHA_MAX,
            "window_label_count": SOURCE_COUNT,
        },
        "exact_geometric_half_line": {
            "one_label": "S_alpha=K_0*integral_0^infinity x^(-s)exp(-pi*x^2-lambda*alpha*x)dx",
            "parameters": "lambda=pi*(1+i)/sqrt(2); K_0=exp(3*pi*t/4+i*3*pi/8)",
            "finite_source": "G_W(x)=sum_(j=0)^(M-1)exp(-lambda*(A+2*j)*x)=exp(-lambda*A*x)*(1-exp(-2*lambda*M*x))/(1-exp(-2*lambda*x))",
            "removable_endpoint": "G_W(0)=M=2481423",
            "denominator_nonzero_for_real_x_positive": True,
            "compressed_source": "S_W=K_0*integral_0^infinity x^(-s)exp(-pi*x^2)G_W(x)dx",
            "physical_projection": "H(t)*Q_K=Hardy_t[S_W]",
            "finite_sum_integral_interchange_requires_limit_theorem": False,
        },
        "surrogate_compression_check": {
            "decimal_digits": 65,
            "t": str(SURROGATE_T),
            "alpha_min": SURROGATE_ALPHA_MIN,
            "alpha_max": SURROGATE_ALPHA_MAX,
            "direct_label_sum": mp.nstr(direct, 52),
            "compressed_integral": mp.nstr(compressed, 52),
            "absolute_discrepancy": mp.nstr(discrepancy, 18),
        },
        "saddle_geometry": saddle,
        "candidate_saddle_ray": {
            "substitution": "x=exp(-3*i*pi/4)*y with Log(x)=log(y)-3*i*pi/4",
            "integrand_identity": "K_0*x^(-s)*exp(-pi*x^2-lambda*alpha*x)dx=y^(-s)*exp(i*pi*(alpha*y-y^2))dy",
            "stationary_equation": "2*pi*y^2-pi*alpha*y+t=0",
            "same_transition_alpha": "sqrt(8*t/pi)",
            "integrand_normalization_cancellation_exact": True,
            "positive_half_line_to_saddle_ray_contour_deformation_certified": False,
            "Stokes_multiplier_and_large_arc_contributions_owned": False,
        },
        "decision": {
            "finite_QK_roster_compressed_to_one_exact_integral": True,
            "direct_Kummer_roster_evaluation_required": False,
            "actual_lower_label_is_first_odd_above_saddle_transition": True,
            "candidate_saddle_ray_selected_for_connection_audit": True,
            "candidate_saddle_ray_permitted_as_integral_identity": False,
            "actual_height_QK_enclosed": False,
            "D_K_enclosed": False,
            "non_A_bound_proved": False,
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
            "numerical_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Prove an exact connection formula from the positive half-line in (GH3) to finite steepest-descent contours through the two roots of 2*pi*x^2+lambda*alpha*x+s=0, uniformly over the actual finite roster. The proof must retain the zero-branch monodromy, Stokes multipliers, and every large-arc segment. Only then use the common geometric factor to seek a signed actual-height enclosure.",
        "proof_boundary": "Exact finite geometric half-line compression, removable endpoint, surrogate quadrature, exact saddle polynomial, actual threshold arithmetic, and an algebraic candidate-ray integrand identity only. No saddle-ray contour deformation, Stokes multiplier, actual-height Q_K, D_K, Delta_KU, J_Z, non-A, all-height, Lambda<=0, PF-infinity, RH, or prize-level enclosure is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified exact finite Q_K geometric half-line compression and guarded saddle route", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
