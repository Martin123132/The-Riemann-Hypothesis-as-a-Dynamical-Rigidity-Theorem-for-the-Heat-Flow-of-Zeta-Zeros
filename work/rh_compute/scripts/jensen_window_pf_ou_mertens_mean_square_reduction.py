#!/usr/bin/env python3
"""Build the OU-tail/Mertens mean-square and anchored-correlation reduction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_ou_mertens_mean_square_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_ou_mertens_mean_square_reduction.md"
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
            "ommsr_01_definitions",
            "exact_definition",
            "available_exact",
            "For 0<alpha<1, set r_k=sum_(n>=k)mu(n)/n, "
            "M(k)=sum_(n<=k)mu(n), B_alpha=sum k^(-alpha)|r_k|^2, "
            "and M_alpha=sum M(k)^2/k^(2+alpha).",
            "The prime number theorem supplies convergence and r_1=0.",
        ),
        row(
            "ommsr_02_pnt_tail",
            "classical_theorem_step",
            "source_backed",
            "The prime number theorem gives sum_(n>=1)mu(n)/n=0 and "
            "M(N)/N->0.",
            "Classical input only; it supplies no RH-strength power saving.",
        ),
        row(
            "ommsr_03_tail_difference",
            "exact_identity",
            "available_exact",
            "r_k-r_(k+1)=mu(k)/k.",
            "Immediate from the reciprocal-Mobius tail definition.",
        ),
        row(
            "ommsr_04_mertens_from_tail",
            "exact_identity",
            "available_exact",
            "M(k)=sum_(n=1)^k r_n-k*r_(k+1), hence "
            "M(k)/k=H(r)(k)-r_(k+1).",
            "Telescope sum n(r_n-r_(n+1)); H is the prefix Hardy operator.",
        ),
        row(
            "ommsr_05_tail_from_mertens",
            "exact_identity",
            "available_exact",
            "r_(k+1)=-M(k)/(k+1)+sum_(n>=k+1)M(n)/(n(n+1)).",
            "Abel summation and M(N)/N->0.",
        ),
        row(
            "ommsr_06_weighted_hardy",
            "classical_theorem_step",
            "source_backed",
            "On l2(k^(-alpha)), H has norm at most "
            "h_alpha=sqrt(2(3+alpha))/(1+alpha).",
            "A nonsharp explicit Schur-test form of the weighted Hardy "
            "inequality; no Mobius cancellation is used.",
        ),
        row(
            "ommsr_07_weighted_copson",
            "classical_theorem_step",
            "source_backed",
            "The tail operator Cb(k)=sum_(n>=k)b_n/n has norm at most "
            "c_alpha=sqrt(2(3-alpha))/(1-alpha) on l2(k^(-alpha)).",
            "A nonsharp explicit Schur-test form of the weighted Copson "
            "inequality; the constant blows up as alpha->1.",
        ),
        row(
            "ommsr_08_tail_mertens_norm_equivalence",
            "exact_inequality",
            "available_exact",
            "M_alpha^(1/2)<=(h_alpha+2^(alpha/2))B_alpha^(1/2) "
            "and B_alpha^(1/2)<=(1+c_alpha)M_alpha^(1/2).",
            "Combine rows 4-7 with the one-step weighted shift bounds.",
        ),
        row(
            "ommsr_09_ou_weight_comparison",
            "exact_inequality",
            "available_exact",
            "For w_(alpha,k)=k^(1-alpha)-(k-1)^(1-alpha), "
            "(1-alpha)k^(-alpha)<=w_(alpha,k)<=k^(-alpha).",
            "Concavity and 1-alpha<=2^(-alpha); valid for every k>=1.",
        ),
        row(
            "ommsr_10_infinite_ou_mertens_equivalence",
            "exact_equivalence",
            "available_exact",
            "The limiting OU tail energy sum w_(alpha,k)|r_k|^2 is "
            "finite iff M_alpha is finite.",
            "Rows 8 and 9 give explicit two-sided norm constants.",
        ),
        row(
            "ommsr_11_continuous_mertens_cell_identity",
            "exact_identity",
            "available_exact",
            "I_alpha=integral_1^infinity M(x)^2*x^(-2-alpha)dx "
            "equals (1/(1+alpha))*sum M(k)^2*"
            "[k^(-1-alpha)-(k+1)^(-1-alpha)].",
            "M(x)=M(k) on each half-open cell [k,k+1).",
        ),
        row(
            "ommsr_12_discrete_continuous_comparison",
            "exact_inequality",
            "available_exact",
            "2^(-2-alpha)M_alpha<=I_alpha<=M_alpha.",
            "Bound x^(-2-alpha) above and below on each unit cell.",
        ),
        row(
            "ommsr_13_mertens_mellin_identity",
            "exact_identity",
            "available_exact",
            "1/zeta(s)=s*integral_1^infinity M(x)x^(-s-1)dx "
            "initially for Re(s)>1.",
            "Classical Abel/Mellin identity for the Mertens function.",
        ),
        row(
            "ommsr_14_mellin_plancherel",
            "classical_theorem_step",
            "source_backed",
            "If I_alpha<infinity and beta=(1+alpha)/2, then "
            "I_alpha=(1/(2pi))*integral_R "
            "dt/[(beta^2+t^2)|zeta(beta+it)|^2].",
            "Mellin-Plancherel plus analytic continuation from Re(s)>1; "
            "finiteness rules out a boundary pole.",
        ),
        row(
            "ommsr_15_dyadic_block_comparison",
            "exact_inequality",
            "available_exact",
            "For D_j=2^(-2j)sum_(2^j<=k<2^(j+1))M(k)^2, "
            "2^(-2-alpha)sum 2^(-alpha*j)D_j<=M_alpha<="
            "sum 2^(-alpha*j)D_j.",
            "Compare k with 2^j on each dyadic block.",
        ),
        row(
            "ommsr_16_root_growth_equivalence",
            "exact_equivalence",
            "available_exact",
            "M_alpha<infinity for every alpha>0 iff "
            "limsup_(j->infinity)D_j^(1/j)<=1.",
            "Root test for the nonnegative generating series "
            "sum 2^(-alpha*j)D_j.",
        ),
        row(
            "ommsr_17_dyadic_mean_square_rh_criterion",
            "exact_equivalence",
            "available_exact",
            "RH iff for every epsilon>0, "
            "sum_(K<=k<2K)M(k)^2=O_epsilon(K^(2+epsilon)) "
            "on dyadic K.",
            "Rows 10, 15, and 16 composed with the cofinal OU criterion; "
            "this records an equivalent antecedent, not its proof.",
        ),
        row(
            "ommsr_18_prefix_pair_correlation",
            "exact_identity",
            "available_exact",
            "sum_(k<=X)M(k)^2=sum_(n<=X)mu(n)^2(X-n+1)+"
            "2sum_(h<X)sum_(m<=X-h)mu(m)mu(m+h)(X-m-h+1).",
            "Expand M(k)^2 and count the k containing each ordered pair.",
        ),
        row(
            "ommsr_19_cumulative_correlation",
            "exact_identity",
            "available_exact",
            "With C_h(Y)=sum_(m<=Y)mu(m)mu(m+h), the off-diagonal "
            "in row 18 equals 2sum_(h<X)sum_(Y<=X-h)C_h(Y).",
            "Reverse the m,Y summations exactly.",
        ),
        row(
            "ommsr_20_anchored_correlation_criterion",
            "exact_equivalence",
            "available_exact",
            "The dyadic RH criterion is equivalent to "
            "sum_(h<X)sum_(Y<=X-h)C_h(Y)=O_epsilon(X^(2+epsilon)); "
            "the squarefree diagonal is O(X^2).",
            "Rows 17-19; the correlation sum is signed and anchored at m=1.",
        ),
        row(
            "ommsr_21_generic_hardy_barrier",
            "proof_guard",
            "guard_validated",
            "Generic weighted Hardy applied to M(k)=sum_(n<=k)mu(n) "
            "costs sum_(n<=N)mu(n)^2*n^(-alpha), of order N^(1-alpha).",
            "The same divergent scale as the absolute OU and remainder "
            "bounds; a generic sequence inequality cannot close RH.",
        ),
        row(
            "ommsr_22_averaged_chowla_source",
            "literature_guard",
            "source_backed",
            "Matomaki-Radziwill-Tao prove an o(H^k X) average over all "
            "shift k-tuples, and state that the result also holds for mu.",
            "Published averaged-Chowla input only; no anchored cumulative "
            "Mertens estimate is attributed to that theorem.",
        ),
        row(
            "ommsr_23_averaged_chowla_mismatch",
            "proof_guard",
            "guard_validated",
            "The published theorem averages every base shift, whereas "
            "C_h(Y) fixes the origin; even an anchored terminal bound "
            "sum_h|C_h(Y)|=o(Y^2), integrated naively in Y, has cubic "
            "rather than X^(2+epsilon) scale.",
            "Blocks promotion of averaged or logarithmically averaged "
            "Chowla into the required origin-anchored cumulative theorem.",
        ),
        row(
            "ommsr_24_endpoint_slice_guard",
            "proof_guard",
            "guard_validated",
            "Almost-all short-interval cancellation leaves the single "
            "prefix M(X) and the h_1=0 correlation slice uncontrolled.",
            "An exceptional anchored slice can carry the RH-strength "
            "obligation even when its density among shifts tends to zero.",
        ),
        row(
            "ommsr_25_open_anchored_mean_square_gate",
            "open_target",
            "open_target",
            "Prove D_j=2^(o(j)), equivalently the dyadic Mertens "
            "mean-square or signed anchored integrated-correlation bound, "
            "without assuming a zero-free half-plane or RH.",
            "This is RH-equivalent. No OU bound, full Burnol bound, RH, "
            "PF-infinity, or Lambda<=0 conclusion is recorded.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_ou_mertens_mean_square_reduction",
        "date": "2026-07-23",
        "status": (
            "exact OU-tail/Mertens/dyadic/anchored-correlation reduction "
            "with one open RH-equivalent mean-square gate"
        ),
        "parameters": {
            "alpha_range": "0<alpha<1",
            "beta": "(1+alpha)/2",
            "dyadic_K": "2^j",
        },
        "source_anchors": [
            "https://arxiv.org/abs/2208.06141",
            "https://arxiv.org/abs/2508.00388",
            "https://arxiv.org/abs/1503.05121",
            "https://annals.math.princeton.edu/2016/183-3/p06",
        ],
        "rows": rows,
        "audit": {
            "row_count": 25,
            "exact_reduction_count": 16,
            "classical_theorem_step_count": 4,
            "proof_guard_count": 4,
            "literature_guard_count": 1,
            "open_mean_square_gate_count": 1,
            "dyadic_mean_square_rh_equivalence_proved": True,
            "averaged_chowla_closes_gate": False,
            "ou_uniform_bound_proved": False,
            "full_burnol_bound_proved": False,
            "rh_proved": False,
            "pf_infinity_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note() -> str:
    return r"""# Jensen-Window PF OU/Mertens Mean-Square Reduction

Date: 2026-07-23

Status: exact OU-tail/Mertens/dyadic/anchored-correlation reduction with one
open RH-equivalent mean-square gate. This is not a proof of RH, PF-infinity,
or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_ou_mertens_mean_square_reduction.json
python work/rh_compute/scripts/jensen_window_pf_ou_mertens_mean_square_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_ou_mertens_mean_square_reduction.py
```

## Tail And Mertens Coordinates

Fix `0<alpha<1`. The prime number theorem gives

```text
sum_(n>=1)mu(n)/n=0,
M(N)/N->0,
```

so define

```text
r_k:=sum_(n>=k)mu(n)/n,
M(k):=sum_(n<=k)mu(n),

B_alpha:=sum_(k>=1)k^(-alpha)|r_k|^2,
M_alpha:=sum_(k>=1)M(k)^2/k^(2+alpha).                       (OMMSR.1)
```

The difference and two Abel identities are

```text
r_k-r_(k+1)=mu(k)/k,                                         (OMMSR.2)

M(k)=sum_(n=1)^k r_n-k*r_(k+1),                              (OMMSR.3)

r_(k+1)
 =-M(k)/(k+1)+sum_(n>=k+1)M(n)/(n(n+1)).                    (OMMSR.4)
```

Let

```text
Hf(k)=(1/k)sum_(n<=k)f(n),
Cf(k)=sum_(n>=k)f(n)/n.
```

A direct Schur test on `l2(k^(-alpha))` gives the nonsharp explicit
constants

```text
||H||<=h_alpha:=sqrt(2(3+alpha))/(1+alpha),
||C||<=c_alpha:=sqrt(2(3-alpha))/(1-alpha).                  (OMMSR.5)
```

Equations (OMMSR.3)-(OMMSR.5), including the one-step shift, imply

```text
M_alpha^(1/2)
 <=(h_alpha+2^(alpha/2))*B_alpha^(1/2),

B_alpha^(1/2)
 <=(1+c_alpha)*M_alpha^(1/2).                                (OMMSR.6)
```

Thus the reciprocal-Mobius tail norm is finite exactly when the weighted
Mertens mean square is finite.

## OU Energy Is The Same Barrier

For

```text
w_(alpha,k):=k^(1-alpha)-(k-1)^(1-alpha),
```

concavity gives

```text
(1-alpha)k^(-alpha)
 <=w_(alpha,k)
 <=k^(-alpha).                                               (OMMSR.7)
```

Therefore the limiting OU energy from Lemma 11.22Q satisfies

```text
(1-alpha)B_alpha
 <=sum_(k>=1)w_(alpha,k)|r_k|^2
 <=B_alpha,                                                  (OMMSR.8)
```

and is norm-equivalent to `M_alpha`. The OU extraction supplies useful
positive structure, but it does not bypass the classical Mertens
mean-square barrier.

The continuous version is exact:

```text
I_alpha
 :=integral_1^infinity M(x)^2*x^(-2-alpha)dx

 =(1/(1+alpha))*sum_(k>=1)M(k)^2
   [k^(-1-alpha)-(k+1)^(-1-alpha)],                          (OMMSR.9)

2^(-2-alpha)M_alpha<=I_alpha<=M_alpha.                       (OMMSR.10)
```

The classical Mellin identity

```text
1/zeta(s)
 =s*integral_1^infinity M(x)x^(-s-1)dx                      (OMMSR.11)
```

starts in `Re(s)>1`. If `I_alpha<infinity`, Cauchy-Schwarz continues the
Mellin transform to `Re(s)>(1+alpha)/2`, and Mellin-Plancherel gives, with
`beta=(1+alpha)/2`,

```text
I_alpha
 =(1/(2*pi))*integral_R
   dt/[(beta^2+t^2)|zeta(beta+it)|^2].                       (OMMSR.12)
```

Finiteness excludes a reciprocal-zeta pole on that boundary. This is a
reciprocal-zeta Hardy energy, not an unconditional estimate for it.

## Exact Dyadic Criterion

Define

```text
D_j
 :=2^(-2j)sum_(2^j<=k<2^(j+1))M(k)^2.                       (OMMSR.13)
```

Block comparison gives

```text
2^(-2-alpha)sum_(j>=0)2^(-alpha*j)D_j
 <=M_alpha
 <=sum_(j>=0)2^(-alpha*j)D_j.                               (OMMSR.14)
```

The root test therefore gives

```text
M_alpha<infinity for every alpha>0
 iff
limsup_(j->infinity)D_j^(1/j)<=1.                            (OMMSR.15)
```

Composed with the cofinal OU criterion:

```text
RH
 iff
for every epsilon>0,
sum_(K<=k<2K)M(k)^2
 =O_epsilon(K^(2+epsilon))
on dyadic K.                                                 (OMMSR.16)
```

This is an exact RH-equivalent average-cancellation target. Equation
(OMMSR.16) has not been proved here.

## Anchored Correlation Form

Put

```text
S(X):=sum_(k<=X)M(k)^2,
C_h(Y):=sum_(m<=Y)mu(m)mu(m+h).
```

Expanding every prefix square and counting the containing prefixes gives

```text
S(X)
 =sum_(n<=X)mu(n)^2(X-n+1)

 +2sum_(h=1)^(X-1)sum_(m=1)^(X-h)
   mu(m)mu(m+h)(X-m-h+1)                                    (OMMSR.17)

 =sum_(n<=X)mu(n)^2(X-n+1)
  +2sum_(h=1)^(X-1)sum_(Y=1)^(X-h)C_h(Y).                   (OMMSR.18)
```

The squarefree diagonal is `O(X^2)`. Hence (OMMSR.16) is equivalently the
signed, origin-anchored integrated-correlation estimate

```text
sum_(h=1)^(X-1)sum_(Y=1)^(X-h)C_h(Y)
 =O_epsilon(X^(2+epsilon)).                                  (OMMSR.19)
```

The signs and the cumulative `Y` summation are essential. Taking absolute
values creates a much stronger target.

## Two Nonpromotion Guards

First, a generic weighted Hardy inequality applied to
`M(k)=sum_(n<=k)mu(n)` gives only

```text
sum_(k<=N)M(k)^2/k^(2+alpha)
 <=C_alpha*sum_(n<=N)mu(n)^2*n^(-alpha)
 =O_alpha(N^(1-alpha)).                                      (OMMSR.20)
```

This reproduces the divergent absolute Gram scale. Hardy/Copson machinery
transfers coordinates but supplies no Mobius cancellation.

Second, Matomaki, Radziwill, and Tao prove an averaged Chowla theorem of
the form

```text
sum_(h_1,...,h_k<=H)
 |sum_(n<=X)mu(n+h_1)...mu(n+h_k)|
 =o(H^k*X),                                                  (OMMSR.21)
```

with the Mobius extension stated in the paper. This averages all base
shifts. By contrast, `C_h(Y)` fixes the exceptional origin slice
`h_1=0` and then integrates every terminal length `Y`. Even a hypothetical
anchored terminal estimate

```text
sum_(h<=Y)|C_h(Y)|=o(Y^2)
```

used naively for every `Y<=X` gives cubic scale after summing over `Y`,
not the `X^(2+epsilon)` target. Almost-all short-interval cancellation
likewise leaves the single prefix `M(X)` uncontrolled. Therefore averaged
Chowla, logarithmically averaged Chowla, and almost-all short-interval
cancellation cannot be promoted into (OMMSR.19) without a new anchored
transfer theorem and an additional full power of cancellation.

## Open Handoff

The exact live target is any one of the equivalent statements

```text
D_j=2^(o(j)),

sum_(K<=k<2K)M(k)^2=O_epsilon(K^(2+epsilon)),

sum_(h<X)sum_(Y<=X-h)C_h(Y)=O_epsilon(X^(2+epsilon)),
```

proved without assuming a zero-free half-plane, RH, or the desired OU
bound. The reduction does not prove this target, the full Burnol norm,
RH, PF-infinity, or `Lambda <= 0`.

## Source Boundary

- [Lee and Leong](https://arxiv.org/abs/2208.06141) record modern explicit
  Mertens and reciprocal-zeta bounds and the standard Mellin connection.
  The exact tail/Mertens and dyadic reductions above are derived here.
- [Das and Manna](https://arxiv.org/abs/2508.00388) review and sharpen
  discrete Hardy/Copson inequalities. Only generic operator bounds are
  used here, and they do not supply arithmetic cancellation.
- [Matomaki, Radziwill, and Tao](https://arxiv.org/abs/1503.05121) prove
  averaged Chowla and explicitly note the Mobius extension. Their theorem
  does not state the anchored cumulative estimate (OMMSR.19).
- [Matomaki and Radziwill](https://annals.math.princeton.edu/2016/183-3/p06)
  prove cancellation for Mobius in almost all short intervals. The
  exceptional origin-anchored prefix remains outside that conclusion.
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
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(), encoding="utf-8")
    print(
        "wrote OU/Mertens mean-square reduction: "
        f"{len(payload['rows'])} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
