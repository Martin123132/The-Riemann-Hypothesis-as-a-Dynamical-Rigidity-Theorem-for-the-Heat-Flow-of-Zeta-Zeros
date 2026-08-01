#!/usr/bin/env python3
"""Build the spectral anchor/high-mode reduction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.md"
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
            "sahr_01_spectral_input",
            "exact_equivalence",
            "available_exact",
            "Corollary 11.22Z.1 gives R_alpha<infinity iff "
            "S_alpha:=sum_K K^(-alpha)sum_r"
            "|X_(K,r)|^2/(2r-1)^2<infinity.",
            "The mixed-boundary spectral handoff is independently "
            "validated and is the sole input.",
        ),
        row(
            "sahr_02_cutoff",
            "exact_definition",
            "available_exact",
            "For 1<=R<=K split the mode sum into r<=R and r>R.",
            "The cutoff may depend on K and the fixed positive alpha.",
        ),
        row(
            "sahr_03_current_anchor_split",
            "exact_identity",
            "available_exact",
            "X_(K,r)=Y_(K,r)+beta_K*h_(K,r), where "
            "Y_(K,r)=sum_(i=1)^(K-1)mu(K+i)phi_(K,r)(i) and "
            "h_(K,r)=phi_(K,r)(K).",
            "The current block and future affine anchor are separated "
            "without estimating either one.",
        ),
        row(
            "sahr_04_endpoint_phase",
            "exact_identity",
            "available_exact",
            "h_(K,r)=2(-1)^(r-1)(2K+1)^(-1/2)"
            "cos(theta_(K,r)/2).",
            "Use K*theta=(2r-1)pi/2-theta/2.",
        ),
        row(
            "sahr_05_endpoint_envelope",
            "exact_bound",
            "available_exact",
            "|h_(K,r)|^2<=4/(2K+1)<=2/K.",
            "The future anchor is spread over all modes at boundary "
            "amplitude O(K^(-1/2)).",
        ),
        row(
            "sahr_06_current_parseval",
            "exact_bound",
            "available_exact",
            "sum_(r=1)^K|Y_(K,r)|^2="
            "sum_(i=1)^(K-1)mu(K+i)^2<=K.",
            "This uses the exact sine orthonormality and no Mobius "
            "cancellation.",
        ),
        row(
            "sahr_07_odd_tail",
            "exact_bound",
            "available_exact",
            "For R>=1, sum_(r>R)(2r-1)^(-2)<=1/(2R+1).",
            "Compare the decreasing tail with its elementary integral.",
        ),
        row(
            "sahr_08_current_high_modes",
            "exact_bound",
            "available_exact",
            "sum_(r>R)|Y_(K,r)|^2/(2r-1)^2"
            "<=K/(2R+1)^2.",
            "Use Parseval and the smallest high-mode denominator.",
        ),
        row(
            "sahr_09_anchor_high_modes",
            "exact_bound",
            "available_exact",
            "sum_(r>R)|beta_K*h_(K,r)|^2/(2r-1)^2"
            "<=2|beta_K|^2/[K(2R+1)].",
            "The anchor loses one power of R rather than two because it "
            "occupies every boundary mode.",
        ),
        row(
            "sahr_10_total_high_modes",
            "exact_bound",
            "available_exact",
            "sum_(r>R)|X_(K,r)|^2/(2r-1)^2"
            "<=2K/(2R+1)^2+4|beta_K|^2/[K(2R+1)].",
            "Apply |a+b|^2<=2|a|^2+2|b|^2 after the exact split.",
        ),
        row(
            "sahr_11_weighted_high_scale",
            "exact_bound",
            "available_exact",
            "The normalized high scale is at most "
            "2K^(1-alpha)/(2R+1)^2+"
            "4K^(-1-alpha)|beta_K|^2/(2R+1).",
            "This displays the distinct current and anchor thresholds.",
        ),
        row(
            "sahr_12_trivial_anchor",
            "exact_bound",
            "available_exact",
            "|beta_K|<=2K.",
            "There are 2K future coefficients and every weight lies in "
            "[0,1].",
        ),
        row(
            "sahr_13_unconditional_cutoff",
            "exact_definition",
            "available_exact",
            "Set R_(alpha,K)=ceil(K^(1-alpha/2)).",
            "For fixed alpha>0 this is sublinear and at most K.",
        ),
        row(
            "sahr_14_unconditional_high_bound",
            "exact_bound",
            "available_exact",
            "At R=R_(alpha,K), the normalized high scale is at most "
            "2K^(-1)+8K^(-alpha/2).",
            "Insert the trivial anchor bound and R>=K^(1-alpha/2).",
        ),
        row(
            "sahr_15_high_summability",
            "exact_bound",
            "available_exact",
            "sum_(K dyadic)K^(-alpha)"
            "sum_(r>R_(alpha,K))|X_(K,r)|^2/(2r-1)^2<infinity.",
            "Both K^(-1) and K^(-alpha/2) are dyadically summable.",
        ),
        row(
            "sahr_16_low_mode_equivalence",
            "exact_equivalence",
            "available_exact",
            "R_alpha<infinity iff sum_(K dyadic)K^(-alpha)"
            "sum_(r<=R_(alpha,K))|X_(K,r)|^2/(2r-1)^2<infinity.",
            "The removed high-mode series is unconditionally finite.",
        ),
        row(
            "sahr_17_positive_alpha_guard",
            "proof_guard",
            "guard_validated",
            "The cutoff and high-tail summability are pointwise in each "
            "fixed alpha>0.",
            "As alpha->0 the exponent approaches one; no alpha=0 or "
            "uniform-in-alpha reduction follows.",
        ),
        row(
            "sahr_18_square_root_anchor",
            "conditional_calibration",
            "conditional",
            "Assume |beta_K|=O_epsilon(K^(1/2+epsilon)) with "
            "2epsilon<alpha.",
            "This power-scale anchor cancellation is not known "
            "unconditionally.",
        ),
        row(
            "sahr_19_conditional_cutoff",
            "conditional_calibration",
            "conditional",
            "For any eta>0 with (1-alpha)/2+eta<1, set "
            "R=ceil(K^((1-alpha)/2+eta)).",
            "The current high term is O(K^(-2eta)); the conditional "
            "anchor term is also dyadically summable.",
        ),
        row(
            "sahr_20_conditional_low_modes",
            "conditional_calibration",
            "conditional",
            "Under row 18, only O(K^((1-alpha)/2+eta)) low modes remain.",
            "This stronger reduction depends on the unproved "
            "square-root anchor estimate.",
        ),
        row(
            "sahr_21_davenport_anchor_guard",
            "literature_guard",
            "source_backed",
            "Davenport and partial summation give beta_K<<_A "
            "K*log^(-A)K.",
            "Arbitrary logarithmic savings do not change the power "
            "threshold R>K^(1-alpha).",
        ),
        row(
            "sahr_22_generic_anchor_guard",
            "proof_guard",
            "guard_validated",
            "For a bounded synthetic sequence equal to one on "
            "[2K,3K], beta_K>=2K/9.",
            "No coefficient-envelope or generic bounded-sequence "
            "argument can assume square-root anchor cancellation.",
        ),
        row(
            "sahr_23_finite_validation",
            "finite_reproduction",
            "finite_check_passed",
            "Arbitrary-vector and Mobius checks verify the endpoint "
            "phase, Parseval split, odd tail, current/anchor inequalities, "
            "and both cutoff calibrations.",
            "Finite validation checks the inequalities only and supplies "
            "no asymptotic cancellation.",
        ),
        row(
            "sahr_24_open_low_mode_gate",
            "open_theorem_target",
            "open",
            "Prove the reduced low-mode criterion in row 16 for every "
            "member of one fixed cofinal sequence alpha_j->0 without "
            "assuming RH.",
            "No low-mode square-root gain, anchor power cancellation, "
            "full Burnol bound, RH, PF-infinity, or Lambda<=0 conclusion "
            "is supplied.",
        ),
    ]
    return {
        "kind": (
            "jensen_window_pf_mertens_"
            "spectral_anchor_high_mode_reduction"
        ),
        "date": "2026-07-23",
        "status": (
            "exact spectral anchor separation and summable high-mode "
            "reduction with one open low-mode gate"
        ),
        "parameters": {
            "alpha_range": "0<alpha<1",
            "dyadic_blocks": "K=2^j",
            "unconditional_cutoff": "ceil(K^(1-alpha/2))",
            "current_tail_decay": "R^(-2)",
            "anchor_tail_decay": "R^(-1)",
        },
        "source_anchors": [
            "outputs/jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.md",
            "https://doi.org/10.1093/qmath/os-8.1.313",
        ],
        "rows": rows,
        "audit": {
            "row_count": 24,
            "exact_reduction_count": 16,
            "literature_guard_count": 1,
            "proof_guard_count": 2,
            "conditional_calibration_count": 3,
            "finite_validation_count": 1,
            "open_signed_gate_count": 1,
            "current_anchor_split_proved": True,
            "unconditional_high_modes_summable": True,
            "unconditional_low_mode_equivalence_proved": True,
            "square_root_anchor_proved": False,
            "low_mode_gain_proved": False,
            "full_burnol_bound_proved": False,
            "rh_proved": False,
            "pf_infinity_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note() -> str:
    return r"""# Jensen-Window PF Mertens Spectral Anchor/High-Mode Reduction

Date: 2026-07-23

Status: exact spectral anchor separation and summable high-mode
reduction with one open low-mode gate. This is not a proof of RH,
PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.json
python work/rh_compute/scripts/jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.py
```

## Exact Current/Anchor Split

Fix `0<alpha<1`. Corollary 11.22Z.1 proves

```text
R_alpha<infinity
 iff
S_alpha
 :=sum_(K dyadic)K^(-alpha)
   sum_(r=1)^K |X_(K,r)|^2/(2r-1)^2
 <infinity.                                           (SAHR.1)
```

Write

```text
X_(K,r)=Y_(K,r)+beta_K*h_(K,r),

Y_(K,r)
 :=sum_(i=1)^(K-1)mu(K+i)phi_(K,r)(i),

h_(K,r):=phi_(K,r)(K).                               (SAHR.2)
```

Since

```text
K*theta_(K,r)
 =(2r-1)pi/2-theta_(K,r)/2,
```

the boundary amplitude is exactly

```text
h_(K,r)
 =2(-1)^(r-1)/sqrt(2K+1)
  *cos(theta_(K,r)/2),

|h_(K,r)|^2<=4/(2K+1)<=2/K.                         (SAHR.3)
```

The current block has full Parseval control:

```text
sum_(r=1)^K|Y_(K,r)|^2
 =sum_(i=1)^(K-1)mu(K+i)^2
 <=K.                                                (SAHR.4)
```

## High-Mode Inequality

For an integer `1<=R<=K`, the odd tail satisfies

```text
sum_(r>R)1/(2r-1)^2<=1/(2R+1).                      (SAHR.5)
```

Therefore

```text
sum_(r>R)|Y_(K,r)|^2/(2r-1)^2
 <=K/(2R+1)^2,

sum_(r>R)|beta_K*h_(K,r)|^2/(2r-1)^2
 <=2|beta_K|^2/[K(2R+1)].                           (SAHR.6)
```

Using `|a+b|^2<=2|a|^2+2|b|^2` gives the exact working bound

```text
sum_(r>R)|X_(K,r)|^2/(2r-1)^2
 <=2K/(2R+1)^2
   +4|beta_K|^2/[K(2R+1)].                          (SAHR.7)
```

The two terms have different origins. The current oscillatory block
decays like `R^(-2)`, while the affine future anchor decays only like
`R^(-1)`.

## Unconditional Summable Tail

Every coefficient in `beta_K` lies in `[0,1]`, so

```text
|beta_K|<=2K.                                        (SAHR.8)
```

Choose

```text
R_(alpha,K):=ceil(K^(1-alpha/2)).                    (SAHR.9)
```

Multiplying (SAHR.7) by `K^(-alpha)` gives

```text
K^(-alpha)
sum_(r>R_(alpha,K))|X_(K,r)|^2/(2r-1)^2
 <=2K^(-1)+8K^(-alpha/2).                           (SAHR.10)
```

Both terms are dyadically summable. Hence the high modes are closed
unconditionally, and

```text
R_alpha<infinity
 iff
sum_(K dyadic)K^(-alpha)
 sum_(1<=r<=R_(alpha,K))
 |X_(K,r)|^2/(2r-1)^2
 <infinity.                                          (SAHR.11)
```

For every fixed positive `alpha`, only
`O(K^(1-alpha/2))` modes remain. This statement is pointwise in
`alpha`; it gives neither an `alpha=0` reduction nor uniformity as
`alpha->0`.

## Conditional Stronger Cutoff

Suppose, conditionally, that

```text
|beta_K|=O_epsilon(K^(1/2+epsilon)),
                  2epsilon<alpha.                   (SAHR.12)
```

For any `eta>0` satisfying `(1-alpha)/2+eta<1`, choose

```text
R=ceil(K^((1-alpha)/2+eta)).
```

Then the current contribution in (SAHR.7), after multiplication by
`K^(-alpha)`, is `O(K^(-2eta))`; the boundary contribution is also
dyadically summable because `2epsilon<alpha`. Under (SAHR.12), only

```text
O(K^((1-alpha)/2+eta))                               (SAHR.13)
```

low modes remain. This stronger cutoff is conditional.

## Guards

Davenport and partial summation give
`beta_K<<_A K*log^(-A)K`, but logarithmic savings do not change the
power threshold in (SAHR.7). A generic bounded-sequence argument cannot
do better: if the synthetic future coefficients equal one on
`[2K,3K]`, then every corresponding weight is at least `2/9`, so
`beta_K>=2K/9`.

Thus the unconditional low-mode count cannot be improved by coefficient
envelopes or generic orthogonality. A better cutoff requires genuine
Mobius power cancellation in the anchor, while closing the remaining
low modes requires the signed additive-twist estimate from Corollary
11.22Z.1.

## Open Gate

For every member of one fixed cofinal sequence `alpha_j->0`, prove the
low-mode series in (SAHR.11) without assuming RH. The high modes are now
an exact closed component. Low-mode square-root cancellation, anchor
power cancellation, the full Burnol bound, RH, PF-infinity, and
`Lambda<=0` all remain open.

Primary source:

- Davenport, *On Some Infinite Series Involving Arithmetical
  Functions (II)*:
  https://doi.org/10.1093/qmath/os-8.1.313
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
        "built Mertens spectral anchor/high-mode reduction: "
        f"{len(payload['rows'])} rows, "
        f"{payload['audit']['exact_reduction_count']} exact reductions, "
        f"{payload['audit']['proof_guard_count']} proof guards, "
        f"{payload['audit']['open_signed_gate_count']} open gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
