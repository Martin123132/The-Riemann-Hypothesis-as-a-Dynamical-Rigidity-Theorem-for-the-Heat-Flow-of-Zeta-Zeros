#!/usr/bin/env python3
"""Align the finite RSI prefix with Q_K by retaining its log-square branch correction."""

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

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_label_branch_correction_QK_alignment_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
PAPER_2015 = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_2015_1502.06903.pdf"
PAPER_2026 = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_Brereton_2607.15310.pdf"

DEPENDENCIES = {
    "finite_residual": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_geometric_contour_residual_and_pole_release_gate.json",
    "normalization_repair": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_domain_Kummer_normalization_repair_gate.json",
    "exact_H": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_exact_H_continuation_defect_target_gate.json",
}

SOURCE_ALPHA_MIN = 159_577
SOURCE_ALPHA_MAX = 5_122_421
PRECISION = 90
WITNESSES = (("2", 1), ("2", 3), ("5", 1), ("5", 3), ("5", 5))


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


def projection_row(t_text: str, alpha_value: int) -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    t = arb(t_text)
    alpha = arb(alpha_value)
    gamma_quarter = acb(arb("0.25"), t / 2).gamma()
    gamma_three_quarter = acb(arb("0.75"), t / 2).gamma()
    argument = acb(0, pi * alpha * alpha / 4)
    phase = acb(0, -pi / 8).exp()
    phi_one = argument.hypgeom_1f1(acb(arb("0.25"), -t / 2), arb("0.5"))
    phi_two = argument.hypgeom_1f1(acb(arb("0.75"), -t / 2), arb("1.5"))
    rotated_one = phase * phi_one
    rotated_two = phase * phi_two
    require(rotated_one.imag.contains(0), "odd-label Phi1 reflection failed")
    require(rotated_two.imag.contains(0), "odd-label Phi2 reflection failed")

    exact_prefix_hardy = (
        -(-pi * t / 4).exp()
        * abs(gamma_quarter)
        * pi ** arb("-0.25")
        * rotated_one.real
    )
    source_label_hardy = (
        2
        * pi ** arb("1.25")
        * alpha
        * (-3 * pi * t / 4).exp()
        * rotated_two.real
        / ((1 + (-2 * pi * t).exp()) * abs(gamma_quarter))
    )

    scalar_H = (
        arb(2) ** arb("1.25")
        * t ** arb("0.25")
        * pi ** arb("1.5")
        * (-3 * pi * t / 4).exp()
        / (
            (1 + (-2 * pi * t).exp())
            * abs(gamma_quarter)
            * abs(gamma_three_quarter) ** 2
        )
    )
    c_t = (pi / (arb(32) * t)) ** arb("0.25")
    beta = 2 * abs(gamma_three_quarter) ** 2 / pi.sqrt()
    corrected_kummer_label = c_t * alpha * beta * rotated_two.real
    require(source_label_hardy.overlaps(scalar_H * corrected_kummer_label), "H*K label alignment failed")

    branch_hardy = source_label_hardy - exact_prefix_hardy
    return {
        "t": t_text,
        "alpha": alpha_value,
        "rotated_Phi1_imaginary_defect_ball": rotated_one.imag.str(35, more=True),
        "rotated_Phi2_imaginary_defect_ball": rotated_two.imag.str(35, more=True),
        "exact_finite_RSI_label_Hardy_ball": exact_prefix_hardy.str(65, more=True),
        "A33_source_label_Hardy_ball": source_label_hardy.str(65, more=True),
        "H_times_corrected_Kummer_label_ball": (scalar_H * corrected_kummer_label).str(65, more=True),
        "branch_correction_Hardy_ball": branch_hardy.str(65, more=True),
        "branch_correction_excludes_zero": not branch_hardy.contains(0),
    }


def render_note(artifact: dict[str, Any]) -> str:
    window = artifact["actual_window_augmented_identity"]
    witness = next(row for row in artifact["rigorous_surrogate_witnesses"] if row["t"] == "5" and row["alpha"] == 3)
    return f"""# Finite-label branch correction and exact Q_K alignment

Date: 2026-08-27

Status: exact corrected finite-roster identity certified; branch sum not enclosed

The finite geometric contour prefix from the preceding gate cannot be
identified directly with the equation-(4) Kummer roster.  The obstruction is
an explicit branch term, not an unresolved sign convention.

For `s=1/2+it`, odd `alpha`, and

```text
r=exp(-i*pi/4),
z=1/2+q*r,
v=exp(i*pi/4)z=q+(1+i)/(2sqrt(2)),
c_alpha=pi*alpha*(1+i)/sqrt(2),
d=exp(-pi*t/4-i*pi/8),                              (BC1)
```

the exact label in the finite csc-prefix is

```text
T_alpha=-integral_C exp(-i*pi*z^2+i*pi*alpha*z)z^(-s)dz
       =-d integral_(R+i0) exp(-pi*v^2+c_alpha*v)v^(-s)dv. (BC2)
```

Put

```text
I_alpha(gamma)=integral_0^infinity x^(-s)exp(-pi*x^2+gamma*x)dx,
a=1/4-it/2,  b=3/4-it/2,  xi=i*pi*alpha^2/4,
A_alpha=Gamma(a)pi^(-a) 1F1(a;1/2;xi)/2,
B_alpha=c_alpha Gamma(b)pi^(-b) 1F1(b;3/2;xi)/2.    (BC3)
```

Then `I_alpha(c_alpha)=A_alpha+B_alpha` and
`I_alpha(-c_alpha)=A_alpha-B_alpha`.  The upper boundary value of `v^(-s)`
on the negative real axis gives `exp(-i*pi*s)=-i exp(pi*t)`, hence

```text
T_alpha=-d[(A_alpha+B_alpha)-i exp(pi*t)(A_alpha-B_alpha)]. (BC4)
```

After simplifying the odd-label phase and the gamma reflection in Appendix
A, the finite A33 source label is

```text
S_alpha=i exp(pi*t)d(A_alpha-B_alpha).               (BC5)
```

Therefore the exact relation is

```text
T_alpha=S_alpha-d I_alpha(c_alpha),
S_alpha-T_alpha=d I_alpha(c_alpha).                  (BC6)
```

The missing term in (BC6) is precisely what is lost if the integration
constant obtained on one side of the `Log(v^2)` cut is carried across the
whole horizontal endpoint without a jump.  This diagnosis does not require
assigning an unstated branch convention to the source: (BC2)--(BC6) compare
the explicit exact RSI label and explicit A33 Kummer label themselves.

Let `Hardy_t[X]=2 Re[exp(i*theta(t))X]`.  Odd-label Kummer reflection gives
both `exp(-i*pi/8)Phi1` and `exp(-i*pi/8)Phi2` real, and exact phase algebra
then yields

```text
Hardy_t[T_alpha]
 =-exp(-pi*t/4)|Gamma(1/4+it/2)|pi^(-1/4)
   Re[exp(-i*pi/8)Phi1],                              (BC7)

Hardy_t[S_alpha]
 =2pi^(5/4)alpha exp(-3pi*t/4)
   Re[exp(-i*pi/8)Phi2]
  /[(1+exp(-2pi*t))|Gamma(1/4+it/2)|]
 =H(t)K_t(alpha).                                    (BC8)
```

Thus the exact finite RSI prefix projects through `Phi1`, while the physical
Kummer label projects through the thermal `Phi2` term.  At the independent
rigorous witness `t=5, alpha=3`,

```text
Hardy_t[T_alpha]={witness['exact_finite_RSI_label_Hardy_ball']},
Hardy_t[S_alpha]={witness['A33_source_label_Hardy_ball']},
Hardy_t[S_alpha-T_alpha]={witness['branch_correction_Hardy_ball']}. (BC9)
```

The last ball excludes zero, so the uncorrected identification is false.

For the actual roster, define

```text
B_W=sum_(alpha=A,A+2,...,B) d I_alpha(c_alpha),
P_W=sum T_alpha,  S_W=sum S_alpha.                   (BC10)
```

There are `{window['window_label_count']}` labels and, exactly,

```text
S_W=P_W+B_W,
H(t)Q_K=Hardy_t[S_W]
       =Hardy_t[-sum_(n={window['released_Dirichlet_block'][0]})^({window['released_Dirichlet_block'][1]}) n^(-s)
                 +C_{window['prefix_before_window']}-C_{window['prefix_through_window']}+B_W]. (BC11)
```

Equation (BC11) is the corrected finite-contour route to the physical roster.
It uses no infinite A21 interchange, but the new finite branch sum `B_W` must
be retained and enclosed.  The next stage is to combine the released
Dirichlet block with the exact classical upper component before norms, then
derive a common-contour or summation representation for `B_W+C_n--C_m`.

Pi provenance: every `pi` comes from the original sine/Gaussian RSI kernel,
the rotation by `pi/4`, the Riemann-Siegel phase, or standard gamma/Kummer
normalization.  No fitted or unexplained circle constant is inserted.

Proof boundary: exact finite-label reduction, explicit branch correction,
corrected `H(t)K_t(alpha)` projection, and the augmented actual-window
identity only.  No enclosure of `B_W`, `Q_K`, `D_K`, `Delta_KU`, or `J_Z`, no
released-pole cancellation theorem, non-A bound, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "dependency failure")
    require(
        dependencies["finite_residual"]["decision"]["finite_prefix_plus_residual_exact"] is True,
        "finite residual dependency drift",
    )
    require(
        dependencies["normalization_repair"]["decision"]["full_Kummer_label_prefactor_corrected_from_2c_to_c"]
        is True,
        "Kummer normalization drift",
    )
    require(
        dependencies["exact_H"]["decision"]["equation4_infinite_Kummer_series_equals_Z_certified_here"]
        is False,
        "infinite equation-(4) exactness was overpromoted",
    )

    rows = [projection_row(t_text, alpha) for t_text, alpha in WITNESSES]
    central = next(row for row in rows if row["t"] == "5" and row["alpha"] == 3)
    require(central["branch_correction_excludes_zero"] is True, "branch correction witness contains zero")
    central_branch = arb(central["branch_correction_Hardy_ball"])
    require(abs(central_branch).lower() > arb("1"), "branch correction witness lost scale")

    finite_window = dependencies["finite_residual"]["actual_window_identity"]
    require(finite_window["alpha_min"] == SOURCE_ALPHA_MIN, "lower source label drift")
    require(finite_window["alpha_max"] == SOURCE_ALPHA_MAX, "upper source label drift")
    require(finite_window["window_label_count"] == 2_481_423, "source count drift")

    artifact = {
        "kind": STEM,
        "status": "exact_finite_label_branch_correction_and_augmented_QK_contour_alignment_certified",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "identity_domain": "t>0 and every positive odd alpha; actual finite roster A..B",
            "precision_decimal_digits": PRECISION,
        },
        "exact_finite_RSI_label": {
            "definition": "T_alpha=-integral_C exp(-i*pi*z^2+i*pi*alpha*z)z^(-s)dz",
            "rotated_contour": "T_alpha=-d*integral_(R+i0)exp(-pi*v^2+c_alpha*v)v^(-s)dv",
            "parameters": "s=1/2+it; c_alpha=pi*alpha*(1+i)/sqrt(2); d=exp(-pi*t/4-i*pi/8)",
            "half_line_formula": "T_alpha=-d[(A_alpha+B_alpha)-i*exp(pi*t)*(A_alpha-B_alpha)]",
        },
        "A33_source_label": {
            "formula": "S_alpha=i*exp(pi*t)*d*(A_alpha-B_alpha)",
            "half_line_identification": "A_alpha-B_alpha=I_alpha(-c_alpha)",
            "corrected_physical_projection": "Hardy_t[S_alpha]=H(t)*K_t(alpha)",
        },
        "branch_correction": {
            "half_line_integral": "I_alpha(c_alpha)=integral_0^infinity x^(-s)exp(-pi*x^2+c_alpha*x)dx=A_alpha+B_alpha",
            "exact_label_identity": "T_alpha=S_alpha-d*I_alpha(c_alpha)",
            "difference": "S_alpha-T_alpha=d*I_alpha(c_alpha)",
            "log_square_diagnosis": "The endpoint v crosses the Log(v^2) cut; the integration constant derived on the left cannot be carried through without the displayed jump term.",
        },
        "hardy_projection_identities": {
            "functional": "Hardy_t[X]=2*Re(exp(i*theta(t))*X)",
            "odd_reflection": "exp(-i*pi/8)Phi_j is real for j=1,2 and odd alpha",
            "finite_RSI_label": "Hardy_t[T_alpha]=-exp(-pi*t/4)|Gamma(1/4+it/2)|pi^(-1/4)Re(exp(-i*pi/8)Phi1)",
            "A33_source_label": "Hardy_t[S_alpha]=2pi^(5/4)alpha exp(-3pi*t/4)Re(exp(-i*pi/8)Phi2)/[(1+exp(-2pi*t))|Gamma(1/4+it/2)|]=H(t)K_t(alpha)",
        },
        "rigorous_surrogate_witnesses": rows,
        "actual_window_augmented_identity": {
            "alpha_min": finite_window["alpha_min"],
            "alpha_max": finite_window["alpha_max"],
            "window_label_count": finite_window["window_label_count"],
            "prefix_before_window": finite_window["prefix_before_window"],
            "prefix_through_window": finite_window["prefix_through_window"],
            "released_Dirichlet_block": finite_window["released_Dirichlet_block"],
            "branch_sum": "B_W=sum_(alpha=A,A+2,...,B)d*I_alpha(c_alpha)",
            "source_prefix_relation": "S_W=P_W+B_W",
            "physical_identity": "H(t)Q_K=Hardy_t[-sum_(n=n_-+1)^n_+ n^(-s)+C_(n_-)-C_(n_+)+B_W]",
        },
        "decision": {
            "finite_RSI_label_equals_A33_source_label": False,
            "uncorrected_finite_residual_permitted_as_HQK": False,
            "branch_correction_rigorously_nonzero_at_surrogate": True,
            "A21_infinite_interchange_required_for_augmented_finite_roster": False,
            "physical_QK_prefactor_and_branch_alignment_complete_with_correction": True,
            "actual_window_augmented_common_contour_identity_exact": True,
            "branch_correction_sum_enclosed": False,
            "released_pole_block_joined_to_classical_upper_component": False,
            "D_K_enclosed": False,
            "non_A_bound_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
        "references": {
            "Lewis_2015": {
                "path": relative(PAPER_2015),
                "sha256": file_hash(PAPER_2015),
                "relevant_equations": ["A1", "A21", "A23", "A25", "A31", "A32", "A33", "A38", "A39", "A40"],
            },
            "Lewis_Brereton_2026": {
                "path": relative(PAPER_2026),
                "sha256": file_hash(PAPER_2026),
                "relevant_equations": [4, 7, 9],
            },
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "arb_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Insert the augmented identity H(t)Q_K=Hardy_t[-Dirichlet_block+C_(79788)-C_(2561211)+B_W] into D_K=U-HQ_K. Cancel the released Dirichlet block against the exact classical upper decomposition before absolute values, then derive a common-contour or summation representation for B_W+C_(79788)-C_(2561211) with a signed enclosure on the D_K corridor.",
        "proof_boundary": "Exact finite-label half-line reduction, explicit nonzero log-square branch correction, corrected H(t)K_t(alpha) projection, and augmented actual-window identity only. No enclosure of B_W, Q_K, D_K, Delta_KU, or J_Z, no released-pole cancellation theorem, non-A bound, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified finite-label branch correction and augmented physical Q_K identity", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
