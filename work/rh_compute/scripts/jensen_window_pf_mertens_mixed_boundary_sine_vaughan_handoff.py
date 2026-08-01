#!/usr/bin/env python3
"""Build the mixed-boundary sine and modewise Vaughan handoff."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.md"
)


def row(
    row_id: str,
    role: str,
    status: str,
    statement: str,
    proof_boundary: str,
) -> dict[str, str]:
    return {
        "id": row_id,
        "role": role,
        "status": status,
        "statement": statement,
        "proof_boundary": proof_boundary,
    }


def build_payload() -> dict:
    rows = [
        row(
            "mbsv_01_vector_input",
            "exact_equivalence",
            "available_exact",
            "Lemma 11.22Z gives R_alpha<infinity iff "
            "E_alpha:=sum_(K dyadic)K^(-2-alpha)||B_K||_2^2"
            "<infinity.",
            "The ordinary feature-vector handoff is independently "
            "validated and is the sole input here.",
        ),
        row(
            "mbsv_02_collapsed_coefficient",
            "exact_definition",
            "available_exact",
            "Set x_(K,i)=mu(K+i) for 1<=i<K and x_(K,K)=beta_K.",
            "The future interval remains inside the exact boundary "
            "coefficient beta_K.",
        ),
        row(
            "mbsv_03_suffix_map",
            "exact_identity",
            "available_exact",
            "B_(K,t)=sum_(i=t)^K x_(K,i), 1<=t<=K.",
            "The terminal row is B_(K,K)=beta_K.",
        ),
        row(
            "mbsv_04_brownian_gram",
            "exact_identity",
            "available_exact",
            "||B_K||_2^2=x_K^T C_K x_K with "
            "C_K(i,j)=min(i,j).",
            "The suffix incidence matrix has Gram C_K.",
        ),
        row(
            "mbsv_05_inverse_laplacian",
            "exact_identity",
            "available_exact",
            "C_K^(-1)=Q_K, where Q_K has off-diagonal -1, diagonal "
            "2 for i<K, and final diagonal 1.",
            "Direct multiplication gives the identity, including K=1.",
        ),
        row(
            "mbsv_06_mixed_boundary",
            "exact_identity",
            "available_exact",
            "Q_K is the discrete Laplacian with v_0=0 and "
            "v_(K+1)=v_K.",
            "These are Dirichlet and half-step Neumann boundary "
            "conditions.",
        ),
        row(
            "mbsv_07_half_odd_angles",
            "exact_definition",
            "available_exact",
            "For 1<=r<=K put theta_(K,r)=(2r-1)pi/(2K+1).",
            "The terminal boundary equation is "
            "cos((2K+1)theta/2)=0.",
        ),
        row(
            "mbsv_08_sine_eigenvectors",
            "exact_definition",
            "available_exact",
            "phi_(K,r)(i)=2(2K+1)^(-1/2)sin(i*theta_(K,r)).",
            "The normalization is chosen for an orthonormal basis of "
            "R^K.",
        ),
        row(
            "mbsv_09_orthonormality",
            "exact_identity",
            "available_exact",
            "sum_(i=1)^K sin(i*theta_(K,r))"
            "sin(i*theta_(K,q))=(2K+1)delta_(r,q)/4.",
            "This is the finite mixed-boundary sine orthogonality "
            "identity.",
        ),
        row(
            "mbsv_10_eigenvalues",
            "exact_identity",
            "available_exact",
            "C_K phi_(K,r)=lambda_(K,r)phi_(K,r), where "
            "lambda_(K,r)=[4sin^2(theta_(K,r)/2)]^(-1).",
            "Equivalently Q_K has eigenvalue "
            "4sin^2(theta_(K,r)/2).",
        ),
        row(
            "mbsv_11_spectral_coefficient",
            "exact_definition",
            "available_exact",
            "X_(K,r):=<x_K,phi_(K,r)>.",
            "These are normalized ordinary-Mobius half-odd sine modes.",
        ),
        row(
            "mbsv_12_anchor_substitution",
            "exact_identity",
            "available_exact",
            "X_(K,r)=2(2K+1)^(-1/2)[sum_(i=1)^(K-1)"
            "mu(K+i)sin(i theta_(K,r))+beta_K sin(K theta_(K,r))].",
            "The future anchor appears only in the mixed-Neumann "
            "boundary amplitude.",
        ),
        row(
            "mbsv_13_spectral_feature",
            "exact_definition",
            "available_exact",
            "Define g_(K,r)(K+i)=phi_(K,r)(i) and "
            "g_(K,r)(n)=b_K(n)phi_(K,r)(K) for 2K<=n<4K.",
            "This unfolds beta_K back into ordinary Mobius values.",
        ),
        row(
            "mbsv_14_spectral_mobius_sum",
            "exact_identity",
            "available_exact",
            "X_(K,r)=sum_(K<n<4K)mu(n)g_(K,r)(n).",
            "Every spectral mode is one compact ordinary-Mobius test.",
        ),
        row(
            "mbsv_15_feature_envelope",
            "exact_bound",
            "available_exact",
            "|g_(K,r)(n)|<=2/sqrt(2K+1) for all n,r.",
            "The test coefficients are uniformly bounded and supported "
            "on K<n<4K.",
        ),
        row(
            "mbsv_16_spectral_energy",
            "exact_identity",
            "available_exact",
            "||B_K||_2^2=sum_(r=1)^K"
            "lambda_(K,r)|X_(K,r)|^2.",
            "This is the exact diagonalization of the min(i,j) Gram.",
        ),
        row(
            "mbsv_17_eigenvalue_lower",
            "exact_bound",
            "available_exact",
            "K^(-2)lambda_(K,r)>=4/[pi^2(2r-1)^2].",
            "Use sin(y)<=y and (2K+1)/K>=2.",
        ),
        row(
            "mbsv_18_eigenvalue_upper",
            "exact_bound",
            "available_exact",
            "K^(-2)lambda_(K,r)<=9/[4(2r-1)^2].",
            "Use sin(y)>=2y/pi on [0,pi/2] and "
            "(2K+1)/K<=3.",
        ),
        row(
            "mbsv_19_odd_frequency_square",
            "exact_definition",
            "available_exact",
            "S_alpha:=sum_(K dyadic)K^(-alpha)"
            "sum_(r=1)^K|X_(K,r)|^2/(2r-1)^2.",
            "This is the mixed-boundary additive-twist square function.",
        ),
        row(
            "mbsv_20_spectral_equivalence",
            "exact_equivalence",
            "available_exact",
            "(4/pi^2)S_alpha<=E_alpha<=(9/4)S_alpha; hence "
            "R_alpha<infinity iff S_alpha<infinity.",
            "The constants are uniform in K and alpha.",
        ),
        row(
            "mbsv_21_absolute_barrier",
            "proof_guard",
            "guard_validated",
            "The feature envelope gives |X_(K,r)|=O(sqrt(K)) and "
            "sum_r|X_(K,r)|^2/(2r-1)^2=O(K).",
            "The normalized scale is O(K^(1-alpha)), one full power "
            "short.",
        ),
        row(
            "mbsv_22_davenport_guard",
            "literature_guard",
            "source_backed",
            "Davenport plus partial summation gives uniformly in r "
            "|X_(K,r)|<<_A sqrt(K)log^(-A)K.",
            "The odd-frequency weight then gives only "
            "K^(1-alpha)log^(-2A)K, never a power saving.",
        ),
        row(
            "mbsv_23_square_root_calibration",
            "conditional_calibration",
            "conditional",
            "If the unnormalized bracket defining X_(K,r) is "
            "O_epsilon(K^(1/2+epsilon)) uniformly in r, then "
            "S_alpha<infinity whenever 2epsilon<alpha.",
            "This sharp square-root-plus-epsilon estimate is not proved "
            "unconditionally.",
        ),
        row(
            "mbsv_24_parseval_guard",
            "proof_guard",
            "guard_validated",
            "Orthogonality alone controls the unweighted mode square but "
            "does not provide Mobius cancellation in the odd-weighted "
            "square.",
            "Generic Parseval or large-sieve bounds reproduce the "
            "one-power loss rather than close it.",
        ),
        row(
            "mbsv_25_vaughan_split",
            "exact_definition",
            "available_exact",
            "Split K<n<4K into X<n<=2X for X=K,2K and set "
            "U_X=V_X=floor(sqrt(X)).",
            "The zero test value at n=4K makes the cover exact.",
        ),
        row(
            "mbsv_26_type_i_coefficients",
            "exact_definition",
            "available_exact",
            "A_X(d)=sum_(bc=d,b<=U_X,c<=V_X)mu(b)mu(c).",
            "These are Green-Tao's finite Type I coefficients.",
        ),
        row(
            "mbsv_27_type_ii_coefficients",
            "exact_definition",
            "available_exact",
            "D_X(d)=sum_(c|d,c>V_X)mu(c).",
            "These are the finite Type II divisor coefficients.",
        ),
        row(
            "mbsv_28_type_i_mode",
            "exact_definition",
            "available_exact",
            "T_(I,K,r,X)=sum_(d<=U_XV_X)A_X(d)"
            "sum_(X/d<w<=2X/d)g_(K,r)(dw).",
            "The inner test is deterministic and bounded.",
        ),
        row(
            "mbsv_29_type_ii_mode",
            "exact_definition",
            "available_exact",
            "T_(II,K,r,X)=sum_(V_X<d<=2X/U_X)D_X(d)"
            "sum_(max(U_X,X/d)<w<=2X/d)mu(w)g_(K,r)(dw).",
            "The bilinear Mobius factor remains inside the signed mode.",
        ),
        row(
            "mbsv_30_modewise_vaughan",
            "source_backed_identity",
            "source_backed",
            "With T_(q,K,r)=T_(q,K,r,K)+T_(q,K,r,2K), "
            "X_(K,r)=-T_(I,K,r)+T_(II,K,r).",
            "Apply the finite Vaughan identity directly to each scalar "
            "spectral feature before squaring.",
        ),
        row(
            "mbsv_31_modewise_criterion",
            "exact_equivalence",
            "available_exact",
            "R_alpha<infinity iff sum_K K^(-alpha)sum_r"
            "|-T_(I,K,r)+T_(II,K,r)|^2/(2r-1)^2<infinity.",
            "This is the exact additive-twist Type I/II square-function "
            "target.",
        ),
        row(
            "mbsv_32_separate_mode_guard",
            "proof_guard",
            "guard_validated",
            "Taking |T_I|+|T_II| mode by mode is not equivalent to the "
            "signed criterion.",
            "A viable bilinear estimate must retain cancellation inside "
            "each half-odd mode and across the weighted mode family.",
        ),
        row(
            "mbsv_33_all_interval_guard",
            "literature_guard",
            "source_backed",
            "All-interval H*log^(-A)X cancellation remains scalar and "
            "gives the same K^(1-alpha)log^(-2A)K mode scale.",
            "It supplies no square-root power or correlated spectral "
            "square estimate.",
        ),
        row(
            "mbsv_34_finite_validation",
            "finite_reproduction",
            "finite_check_passed",
            "Independent finite checks verify C_K^(-1), orthonormality, "
            "eigenpairs, spectral reconstruction, feature unfolding, "
            "and both interval Vaughan modes.",
            "Finite checks validate algebra only and are not evidence "
            "for the open asymptotic estimate.",
        ),
        row(
            "mbsv_35_open_spectral_gate",
            "open_theorem_target",
            "open",
            "Prove the modewise criterion in row 31 for every member of "
            "one fixed cofinal sequence alpha_j->0 without assuming RH.",
            "No square-root additive-twist bound, signed Type I/II "
            "power gain, full Burnol bound, RH, PF-infinity, or "
            "Lambda<=0 conclusion is supplied.",
        ),
    ]
    return {
        "kind": (
            "jensen_window_pf_mertens_"
            "mixed_boundary_sine_vaughan_handoff"
        ),
        "date": "2026-07-23",
        "status": (
            "exact mixed-boundary sine and modewise Vaughan reduction "
            "with one open additive-twist gate"
        ),
        "parameters": {
            "alpha_range": "0<alpha<1",
            "dyadic_blocks": "K=2^j",
            "angles": "theta_(K,r)=(2r-1)pi/(2K+1)",
            "boundary_conditions": "Dirichlet at 0, Neumann at K+1/2",
            "mode_weight": "(2r-1)^(-2)",
        },
        "source_anchors": [
            "outputs/jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.md",
            "https://doi.org/10.1093/qmath/os-8.1.313",
            "https://doi.org/10.5802/aif.2401",
            "https://doi.org/10.1017/fmp.2023.28",
        ],
        "rows": rows,
        "audit": {
            "row_count": 35,
            "exact_reduction_count": 27,
            "literature_guard_count": 2,
            "proof_guard_count": 3,
            "conditional_calibration_count": 1,
            "finite_validation_count": 1,
            "open_signed_gate_count": 1,
            "mixed_boundary_spectrum_proved": True,
            "spectral_mobius_features_proved": True,
            "odd_frequency_equivalence_proved": True,
            "modewise_vaughan_proved": True,
            "square_root_additive_twist_proved": False,
            "signed_type_i_ii_gain_proved": False,
            "full_burnol_bound_proved": False,
            "rh_proved": False,
            "pf_infinity_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note() -> str:
    return r"""# Jensen-Window PF Mertens Mixed-Boundary Sine Vaughan Handoff

Date: 2026-07-23

Status: exact mixed-boundary sine and modewise Vaughan reduction with
one open additive-twist gate. This is not a proof of RH, PF-infinity,
or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.json
python work/rh_compute/scripts/jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.py
```

## Mixed-Boundary Sine Diagonalization

Fix `0<alpha<1`. Lemma 11.22Z gives

```text
R_alpha<infinity
 iff
E_alpha
 :=sum_(K dyadic)K^(-2-alpha)||B_K||_2^2
 <infinity.                                           (MBSV.1)
```

Collapse the future anchor into the terminal coefficient:

```text
x_(K,i):=mu(K+i),       1<=i<K,
x_(K,K):=beta_K,
B_(K,t)=sum_(i=t)^K x_(K,i).                         (MBSV.2)
```

If `C_K(i,j)=min(i,j)`, then

```text
||B_K||_2^2=x_K^T C_K x_K.                           (MBSV.3)
```

The inverse `Q_K=C_K^(-1)` is tridiagonal, with off-diagonal `-1`,
diagonal `2` for `i<K`, and final diagonal `1`. It is the discrete
Laplacian with mixed boundary conditions

```text
v_0=0,                  v_(K+1)=v_K.                 (MBSV.4)
```

For `1<=r<=K`, put

```text
theta_(K,r):=(2r-1)pi/(2K+1),

phi_(K,r)(i)
 :=2/sqrt(2K+1)*sin(i*theta_(K,r)).                  (MBSV.5)
```

The finite sine identity

```text
sum_(i=1)^K
 sin(i*theta_(K,r))*sin(i*theta_(K,q))
 =(2K+1)delta_(r,q)/4
```

shows that the `phi_(K,r)` form an orthonormal basis. Direct substitution
in the interior recurrence and the terminal row gives

```text
C_K phi_(K,r)=lambda_(K,r)phi_(K,r),

lambda_(K,r)
 :=1/[4sin^2(theta_(K,r)/2)].                        (MBSV.6)
```

Define

```text
X_(K,r):=<x_K,phi_(K,r)>.
```

Then

```text
||B_K||_2^2
 =sum_(r=1)^K lambda_(K,r)|X_(K,r)|^2.               (MBSV.7)
```

## Ordinary-Mobius Spectral Features

Write

```text
b_K(n):=((2K)/n)^(1+alpha)*(4K-n)/(2K),
                                      2K<=n<4K.
```

Substituting the exact formula for `beta_K` gives

```text
X_(K,r)
 =2/sqrt(2K+1)*
  [sum_(i=1)^(K-1)mu(K+i)sin(i*theta_(K,r))
   +beta_K*sin(K*theta_(K,r))]

 =sum_(K<n<4K)mu(n)g_(K,r)(n),                       (MBSV.8)
```

where

```text
g_(K,r)(K+i):=phi_(K,r)(i),              1<=i<K,
g_(K,r)(n):=b_K(n)phi_(K,r)(K),          2K<=n<4K.
                                                               (MBSV.9)
```

Every test is compact and satisfies

```text
|g_(K,r)(n)|<=2/sqrt(2K+1).                          (MBSV.10)
```

Thus the spectral modes are normalized ordinary-Mobius additive sine
tests, with the future interval entering only through the terminal
mixed-Neumann amplitude.

## Odd-Frequency Square Function

For

```text
y_(K,r):=(2r-1)pi/(4K+2),
```

the elementary inequalities `2y/pi<=sin(y)<=y` give

```text
4/[pi^2(2r-1)^2]
 <=K^(-2)lambda_(K,r)
 <=9/[4(2r-1)^2].                                   (MBSV.11)
```

Consequently, if

```text
S_alpha
 :=sum_(K dyadic)K^(-alpha)
   sum_(r=1)^K |X_(K,r)|^2/(2r-1)^2,                (MBSV.12)
```

then

```text
(4/pi^2)S_alpha<=E_alpha<=(9/4)S_alpha.              (MBSV.13)
```

In particular,

```text
R_alpha<infinity iff S_alpha<infinity.               (MBSV.14)
```

This is the exact half-odd additive-twist square function.

## Sharp Calibration And Guards

The feature envelope gives `|X_(K,r)|=O(sqrt(K))`; since
`sum_r(2r-1)^(-2)<infinity`, it gives only

```text
K^(-alpha)sum_r |X_(K,r)|^2/(2r-1)^2
 =O(K^(1-alpha)).                                    (MBSV.15)
```

This is one full power short. Davenport's uniform estimate for
Mobius exponential sums, with partial summation for the smooth future
weight, improves this only to

```text
O_A(K^(1-alpha)log^(-2A)K).                          (MBSV.16)
```

No fixed logarithmic saving is dyadically summable here. The
all-interval `H*log^(-A)X` theorem has the same scalar limitation.
Generic Parseval or large-sieve bounds likewise do not create
Mobius-specific cancellation.

The sharp simple sufficient condition is transparent. If, uniformly in
`1<=r<=K`, the bracket in (MBSV.8) satisfies

```text
|sum_(i=1)^(K-1)mu(K+i)sin(i*theta_(K,r))
 +beta_K*sin(K*theta_(K,r))|
 =O_epsilon(K^(1/2+epsilon)),                        (MBSV.17)
```

then `|X_(K,r)|=O_epsilon(K^epsilon)`, and (MBSV.12)
converges for every `2epsilon<alpha`. This square-root-plus-epsilon
bound is a conditional calibration, not an unconditional theorem.

## Vaughan Mode By Mode

Split the support into `X<n<=2X` for `X=K,2K`, and put

```text
U_X=V_X=floor(sqrt(X)),

A_X(d)
 :=sum_(bc=d,b<=U_X,c<=V_X)mu(b)mu(c),

D_X(d)
 :=sum_(c|d,c>V_X)mu(c).                             (MBSV.18)
```

Define

```text
T_(I,K,r,X)
 :=sum_(d<=U_X*V_X)A_X(d)
   sum_(X/d<w<=2X/d)g_(K,r)(dw),

T_(II,K,r,X)
 :=sum_(V_X<d<=2X/U_X)D_X(d)
   sum_(max(U_X,X/d)<w<=2X/d)mu(w)g_(K,r)(dw).
                                                               (MBSV.19)
```

Green and Tao's finite Vaughan identity gives, after summing the two
intervals,

```text
X_(K,r)=-T_(I,K,r)+T_(II,K,r),                       (MBSV.20)
```

where `T_(q,K,r)` is the sum of the `X=K` and `X=2K` pieces. Therefore
the live criterion is exactly

```text
sum_(K dyadic)K^(-alpha)
 sum_(r=1)^K
 |-T_(I,K,r)+T_(II,K,r)|^2/(2r-1)^2
 <infinity.                                          (MBSV.21)
```

Vaughan is applied before the mode square. Replacing the signed
difference by separate absolute Type I and Type II bounds is not an
equivalent step.

## Open Gate

Prove (MBSV.21) for every member of one fixed cofinal sequence
`alpha_j->0` without assuming RH. A viable estimate must deliver
Mobius-specific square-root cancellation on average in the half-odd
mode measure, or an equally strong joint Type I/II bilinear bound.
The square-root estimate, signed power gain, full Burnol bound, RH,
PF-infinity, and `Lambda<=0` all remain open.

Primary sources:

- Davenport, *On Some Infinite Series Involving Arithmetical
  Functions (II)*:
  https://doi.org/10.1093/qmath/os-8.1.313
- Green and Tao, *Quadratic Uniformity of the Mobius Function*,
  Lemma 4.1:
  https://doi.org/10.5802/aif.2401
- Matomaki, Shao, Tao, and Teravainen, *Higher Uniformity of
  Arithmetic Functions in Short Intervals I. All Intervals*:
  https://doi.org/10.1017/fmp.2023.28
"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.result.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.result.write_text(
        json.dumps(payload, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(), encoding="utf-8")
    print(
        "built Mertens mixed-boundary sine/Vaughan handoff: "
        f"{len(payload['rows'])} rows, "
        f"{payload['audit']['exact_reduction_count']} exact reductions, "
        f"{payload['audit']['proof_guard_count']} proof guards, "
        f"{payload['audit']['open_signed_gate_count']} open gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
