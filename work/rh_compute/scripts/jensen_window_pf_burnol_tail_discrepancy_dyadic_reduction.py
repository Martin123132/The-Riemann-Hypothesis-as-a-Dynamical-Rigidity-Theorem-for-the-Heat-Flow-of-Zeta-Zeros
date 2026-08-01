#!/usr/bin/env python3
"""Build the Burnol tail-discrepancy and dyadic square-function reduction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.md"
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
            "btdr_01_weighted_mobius_data",
            "exact_definition",
            "available_exact",
            "alpha=2omega, a_omega(d)=mu(d)d^(-alpha), "
            "A_(omega,N)=sum_(d<=N)a_omega(d)/d.",
            "Definitions for one fixed 0<omega<1/2.",
        ),
        row(
            "btdr_02_limiting_residue",
            "exact_identity",
            "available_exact",
            "A_(omega,infinity)=1/zeta(1+alpha) and "
            "r_(omega,N)=A_(omega,infinity)-A_(omega,N).",
            "Absolute convergence at 1+alpha>1.",
        ),
        row(
            "btdr_03_finite_and_limiting_discrepancies",
            "exact_definition",
            "available_exact",
            "S_(omega,N)(k)=A_(omega,N)k-B_(omega,N)(k), "
            "with the analogous limiting S_(omega,infinity).",
            "Finite divisor sums; B_infinity(k) is finite for each k.",
        ),
        row(
            "btdr_04_tail_divisor_sum",
            "exact_definition",
            "available_exact",
            "H_(omega,N)(k)=B_(omega,infinity)(k)-B_(omega,N)(k)="
            "sum_(N<d<=k)a_omega(d)floor(k/d).",
            "Exact finite omitted-divisor sum.",
        ),
        row(
            "btdr_05_discrepancy_difference",
            "exact_identity",
            "available_exact",
            "D_(omega,N)(k)=S_(omega,N)(k)-"
            "S_(omega,infinity)(k)=H_(omega,N)(k)-k*r_(omega,N).",
            "Exact subtraction of finite and limiting discrepancies.",
        ),
        row(
            "btdr_06_omitted_fractional_tail",
            "exact_identity",
            "available_exact",
            "T_(omega,N)(k)=sum_(d>N)a_omega(d){k/d}="
            "k*r_(omega,N)-H_(omega,N)(k)=-D_(omega,N)(k).",
            "The infinite tail is absolute after d>k and finite before it.",
        ),
        row(
            "btdr_07_abel_tail_operator",
            "exact_identity",
            "available_exact",
            "With g_k(d)=d{k/d}, T_N(k)=r_N*g_k(N+1)+"
            "sum_(d=N+1)^k r_d*(g_k(d+1)-g_k(d)).",
            "Finite summation by parts; g_k is constant for d>k.",
        ),
        row(
            "btdr_08_energy_definitions",
            "exact_definition",
            "available_exact",
            "Q_N=sum_(k>=1)k^(alpha-2)S_N(k)^2, "
            "Q_infinity analogously, and "
            "R_N=sum_(k>=1)k^(alpha-2)D_N(k)^2.",
            "Nonnegative extended-real weighted energies.",
        ),
        row(
            "btdr_09_hilbert_tail_split",
            "exact_equivalence",
            "available_exact",
            "sup_N Q_N<infinity iff Q_infinity<infinity and "
            "sup_N R_N<infinity.",
            "Pointwise convergence plus Fatou in one direction and the "
            "Hilbert triangle inequality in both directions.",
        ),
        row(
            "btdr_10_stable_prefix_energy",
            "exact_identity",
            "available_exact",
            "For k<=N, H_N(k)=0 and D_N(k)=-k*r_N, so "
            "P_N=sum_(k<=N)k^(alpha-2)D_N(k)^2="
            "r_N^2 sum_(k<=N)k^alpha.",
            "Exact stable-divisor prefix.",
        ),
        row(
            "btdr_11_scalar_tail_equivalence",
            "exact_equivalence",
            "available_exact",
            "sup_N P_N<infinity iff "
            "r_N=O(N^(-(1+alpha)/2))="
            "O(N^(-1/2-omega)).",
            "Uses elementary two-sided power-sum bounds.",
        ),
        row(
            "btdr_12_post_prefix_energy",
            "exact_identity",
            "available_exact",
            "U_N=sum_(k>N)k^(alpha-2)|H_N(k)-k*r_N|^2 and "
            "R_N=P_N+U_N.",
            "Exact disjoint split at k=N.",
        ),
        row(
            "btdr_13_three_gate_criterion",
            "exact_equivalence",
            "available_exact",
            "sup_N Q_N<infinity iff Q_infinity<infinity, "
            "r_N=O(N^(-1/2-omega)), and sup_N U_N<infinity.",
            "Composition of the Hilbert split, stable-prefix identity, "
            "and nonnegativity.",
        ),
        row(
            "btdr_14_quotient_reindexing",
            "exact_identity",
            "available_exact",
            "H_N(k)=sum_(m<=floor(k/(N+1)))"
            "[M_alpha(floor(k/m))-M_alpha(N)].",
            "Swap floor(k/d) with its count of multiples.",
        ),
        row(
            "btdr_15_short_interval_reindexing",
            "exact_identity",
            "available_exact",
            "H_N(k)=sum_(q<=floor(k/(N+1)))q*"
            "[M_alpha(floor(k/q))-"
            "M_alpha(max(N,floor(k/(q+1))))].",
            "Partition d by q=floor(k/d); the intervals have ratio "
            "(q+1)/q.",
        ),
        row(
            "btdr_16_dyadic_block_energy",
            "exact_definition",
            "available_exact",
            "V_N(K)=sum_(K<k<=2K)|H_N(k)-k*r_N|^2, "
            "K_j=2^j N.",
            "Dyadic post-prefix block definition.",
        ),
        row(
            "btdr_17_dyadic_equivalence",
            "exact_equivalence",
            "available_exact",
            "2^(alpha-2)sum_j K_j^(alpha-2)V_N(K_j)"
            "<=U_N<=sum_j K_j^(alpha-2)V_N(K_j).",
            "The power weight is monotone on each dyadic block.",
        ),
        row(
            "btdr_18_summable_envelope",
            "exact_implication",
            "available_exact",
            "If V_N(2^jN)<=C*(2^jN)^(2-alpha)*2^(-eta*j) "
            "uniformly for some eta>0, then sup_N U_N<infinity.",
            "A sufficient geometric envelope, not a proved estimate.",
        ),
        row(
            "btdr_19_literature_nonpromotion_guards",
            "literature_guard",
            "source_backed_guard",
            "Baez-Duarte's unweighted natural-approximation lower bounds "
            "and the fractional-part autocorrelation formula motivate but "
            "do not prove the weighted uniform dyadic estimate.",
            "Source comparison only; the weights and finite-tail operator "
            "are different.",
        ),
        row(
            "btdr_20_open_dyadic_square_function",
            "open_target",
            "open_target",
            "Prove Q_infinity<infinity, the reciprocal-zeta tail rate, "
            "and a uniform summable dyadic bound for V_N along one "
            "explicit cofinal omega sequence.",
            "This is RH-strength and remains open; no RH, PF-infinity, "
            "or Lambda<=0 conclusion is recorded.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction",
        "date": "2026-07-23",
        "status": (
            "exact finite-tail and dyadic square-function reduction with "
            "one open arithmetic gate"
        ),
        "parameters": {
            "omega_range": "0<omega<1/2",
            "alpha": "2omega",
            "weight": "k^(alpha-2)",
        },
        "source_anchors": [
            "https://arxiv.org/abs/math/0202166",
            "https://arxiv.org/abs/math/0011254",
            "https://arxiv.org/abs/math/0306251",
        ],
        "rows": rows,
        "audit": {
            "row_count": 20,
            "exact_identity_count": 12,
            "equivalence_step_count": 4,
            "literature_guard_count": 2,
            "open_dyadic_square_function_gate_count": 1,
            "uniform_dyadic_bound_proved": False,
            "rh_proved": False,
            "pf_infinity_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note() -> str:
    return r"""# Jensen-Window PF Burnol Tail-Discrepancy/Dyadic Reduction

Date: 2026-07-23

Status: exact finite-tail and dyadic square-function reduction with one open
arithmetic gate. This is not a proof of RH, PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.json
python work/rh_compute/scripts/jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.py
```

## Finite Versus Limiting Discrepancy

Fix `0<omega<1/2`, set `alpha=2omega`, and write

```text
a_omega(d)=mu(d)d^(-alpha),

A_(omega,N)=sum_(d<=N)a_omega(d)/d,

A_(omega,infinity)=1/zeta(1+alpha),

r_(omega,N)=A_(omega,infinity)-A_(omega,N).         (BTDR.1)
```

For integer `k>=1`, define

```text
B_(omega,N)(k)
 =sum_(d<=N)a_omega(d)floor(k/d),

S_(omega,N)(k)=A_(omega,N)k-B_(omega,N)(k).
```

The limiting divisor sum is finite at every `k`:

```text
B_(omega,infinity)(k)
 =sum_(d<=k)a_omega(d)floor(k/d),

S_(omega,infinity)(k)
 =A_(omega,infinity)k-B_(omega,infinity)(k).
```

Put

```text
H_(omega,N)(k)
 =B_(omega,infinity)(k)-B_(omega,N)(k)
 =sum_(N<d<=k)a_omega(d)floor(k/d).                 (BTDR.2)
```

Direct subtraction gives the exact finite-tail discrepancy

```text
D_(omega,N)(k)
 =S_(omega,N)(k)-S_(omega,infinity)(k)
 =H_(omega,N)(k)-k*r_(omega,N).                    (BTDR.3)
```

Equivalently, the omitted fractional-part tail is

```text
T_(omega,N)(k)
 =sum_(d>N)a_omega(d){k/d}
 =k*r_(omega,N)-H_(omega,N)(k)
 =-D_(omega,N)(k).                                 (BTDR.4)
```

The sum in (BTDR.4) is finite up to `d=k`; beyond `k`,
`{k/d}=k/d`, so the remainder is absolutely convergent.

There is also an exact summation-by-parts form. Let

```text
g_k(d)=d*{k/d},
r_(omega,d)=sum_(m>d)a_omega(m)/m.
```

Since `g_k(d)=k` for `d>k`,

```text
T_(omega,N)(k)
 =r_(omega,N)g_k(N+1)
  +sum_(d=N+1)^k r_(omega,d)
   [g_k(d+1)-g_k(d)].                              (BTDR.5)
```

Thus the remaining operator is driven by jumps of the divisor sawtooth,
not by a generic bounded Hardy operator.

## Exact Three-Gate Split

Define the discrete energies

```text
Q_(omega,N)
 =sum_(k>=1)k^(alpha-2)S_(omega,N)(k)^2,

Q_(omega,infinity)
 =sum_(k>=1)k^(alpha-2)S_(omega,infinity)(k)^2,

R_(omega,N)
 =sum_(k>=1)k^(alpha-2)D_(omega,N)(k)^2.            (BTDR.6)
```

For each fixed `k`, `S_(omega,N)(k)` tends to
`S_(omega,infinity)(k)`. Fatou and the Hilbert triangle inequality give

```text
sup_N Q_(omega,N)<infinity

iff

Q_(omega,infinity)<infinity
and
sup_N R_(omega,N)<infinity.                         (BTDR.7)
```

On the stable prefix `k<=N`, `H_(omega,N)(k)=0`. Hence

```text
P_(omega,N)
 :=sum_(k<=N)k^(alpha-2)D_(omega,N)(k)^2
 =r_(omega,N)^2 sum_(k<=N)k^alpha.                 (BTDR.8)
```

The elementary bounds

```text
N^(1+alpha)/(1+alpha)
 <=sum_(k<=N)k^alpha
 <=N^(1+alpha)
```

show the exact scalar equivalence

```text
sup_N P_(omega,N)<infinity

iff

r_(omega,N)=O(N^(-(1+alpha)/2))
            =O(N^(-1/2-omega)).                    (BTDR.9)
```

The post-prefix energy is

```text
U_(omega,N)
 =sum_(k>N)k^(alpha-2)
  |H_(omega,N)(k)-k*r_(omega,N)|^2,

R_(omega,N)=P_(omega,N)+U_(omega,N).                (BTDR.10)
```

Combining (BTDR.7)-(BTDR.10) yields the exact three-gate criterion

```text
sup_N Q_(omega,N)<infinity

iff

  Q_(omega,infinity)<infinity,
  r_(omega,N)=O(N^(-1/2-omega)),
  sup_N U_(omega,N)<infinity.                       (BTDR.11)
```

The first line is the limiting Jordan-error energy gate. The second is the
forced reciprocal-zeta tail rate. The third is the genuinely finite-tail
square-function gate. None is discarded or hidden inside a generic norm
statement.

## Multiplicative-Interval Form

Let

```text
M_alpha(x)=sum_(d<=floor(x))a_omega(d).
```

Counting `floor(k/d)` by multiples gives

```text
H_(omega,N)(k)
 =sum_(m<=floor(k/(N+1)))
  [M_alpha(floor(k/m))-M_alpha(N)].                 (BTDR.12)
```

Alternatively, partitioning by `q=floor(k/d)` gives

```text
H_(omega,N)(k)
 =sum_(q<=floor(k/(N+1))) q*
  [M_alpha(floor(k/q))
   -M_alpha(max(N,floor(k/(q+1))))].                (BTDR.13)
```

Each bracket in (BTDR.13) is a weighted Mobius increment on the
multiplicative interval

```text
(k/(q+1), k/q],
```

intersected with `(N,infinity)`.

## Dyadic Square Function

For `K>=N`, define

```text
V_(omega,N)(K)
 =sum_(K<k<=2K)
  |H_(omega,N)(k)-k*r_(omega,N)|^2,

K_j=2^j*N.
```

Since `alpha-2<0`, dyadic weight comparison gives

```text
2^(alpha-2) sum_(j>=0) K_j^(alpha-2)V_(omega,N)(K_j)
 <=U_(omega,N)
 <=sum_(j>=0) K_j^(alpha-2)V_(omega,N)(K_j).        (BTDR.14)
```

Therefore the remaining post-prefix condition is exactly a uniform
summable dyadic square function. For example, the sufficient envelope

```text
V_(omega,N)(2^j*N)
 <=C*(2^j*N)^(2-alpha)*2^(-eta*j)                  (BTDR.15)
```

for some `eta>0`, uniformly in `N` and `j`, would close this component.
The critical estimate without the summable factor `2^(-eta*j)` is not
enough by itself.

## Literature Boundary

Burnol records the weighted infinite function and RH-conditional convergence
of its finite natural approximants:

```text
https://arxiv.org/abs/math/0202166
```

Báez-Duarte proves analogous scalar lower-bound obstructions and divergence
results for unweighted natural approximants:

```text
https://arxiv.org/abs/math/0011254
```

The multiplicative fractional-part autocorrelation gives an exact Gram-kernel
framework:

```text
https://arxiv.org/abs/math/0306251
```

Those results motivate (BTDR.5) and the dyadic route, but they do not state or
prove the weighted uniform estimate (BTDR.14). In particular, an unweighted
natural-approximation theorem cannot be transferred by deleting the
`d^(-2omega)` and `k^(2omega-2)` weights.

## Open Gate

For every shift in one explicit cofinal sequence, prove all three lines of
(BTDR.11). The new finite-tail obligation is

```text
sup_N sum_(j>=0)(2^j*N)^(2omega-2)
 V_(omega,N)(2^j*N)<infinity.                      (BTDR.16)
```

Equation (BTDR.16) is open. The exact decomposition does not establish its
hypotheses, RH, PF-infinity, or `Lambda <= 0`.
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
        "wrote Burnol tail-discrepancy/dyadic reduction: "
        f"{len(payload['rows'])} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
