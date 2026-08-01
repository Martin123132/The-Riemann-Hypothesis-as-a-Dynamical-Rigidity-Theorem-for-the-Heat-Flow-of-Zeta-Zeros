#!/usr/bin/env python3
"""Build the Burnol cell-energy and reciprocal-zeta tail obstruction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_burnol_cell_energy_tail_obstruction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_burnol_cell_energy_tail_obstruction.md"
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
            "bcet_01_partial_residue",
            "exact_definition",
            "available_exact",
            "A_(omega,N)=sum_(d<=N)mu(d)d^(-(1+2omega)).",
            "Finite scalar definition.",
        ),
        row(
            "bcet_02_truncated_convolution",
            "exact_definition",
            "available_exact",
            "b_(omega,N)(n)=sum_(d|n,d<=N)mu(d)d^(-2omega).",
            "Finite divisor-sum definition.",
        ),
        row(
            "bcet_03_summatory_floor_identity",
            "exact_identity",
            "available_exact",
            "B_(omega,N)(k)=sum_(n<=k)b_(omega,N)(n)="
            "sum_(d<=N)mu(d)d^(-2omega)floor(k/d).",
            "Exact finite divisor reindexing.",
        ),
        row(
            "bcet_04_reciprocal_fractional_part",
            "exact_identity",
            "available_exact",
            "f_(2omega,N)(1/x)=A_(omega,N)x-"
            "B_(omega,N)(floor(x)).",
            "Exact for noninteger x, with endpoint values irrelevant to L2.",
        ),
        row(
            "bcet_05_weighted_burnol_norm",
            "exact_identity",
            "available_exact",
            "||t^(-omega)f_(2omega,N)||_2^2="
            "integral_0^infinity x^(2omega-2)"
            "|A_(omega,N)x-B_(omega,N)(floor(x))|^2dx.",
            "Exact reciprocal substitution.",
        ),
        row(
            "bcet_06_cell_decomposition",
            "exact_identity",
            "available_exact",
            "The norm is the sum over k>=0 of integral_k^(k+1) "
            "x^(2omega-2)|A_N x-B_N(k)|^2dx.",
            "Nonnegative cell decomposition; Tonelli is immediate.",
        ),
        row(
            "bcet_07_zero_cell",
            "exact_identity",
            "available_exact",
            "The k=0 cell equals A_N^2/(1+2omega).",
            "Exact elementary integral.",
        ),
        row(
            "bcet_08_closed_cell_formula",
            "exact_identity",
            "available_exact",
            "For k>=1 the cell is A_N^2 Delta_(1+2omega)/(1+2omega)"
            "-2A_NB_N(k)Delta_(2omega)/(2omega)"
            "+B_N(k)^2Delta_(2omega-1)/(2omega-1).",
            "Here Delta_a=(k+1)^a-k^a; exact elementary integration.",
        ),
        row(
            "bcet_09_discrete_discrepancy",
            "exact_identity",
            "available_exact",
            "S_(omega,N)(k)=A_(omega,N)k-B_(omega,N)(k)="
            "sum_(d<=N)mu(d)d^(-2omega){k/d}.",
            "Exact finite fractional-part identity.",
        ),
        row(
            "bcet_10_uniform_residue_bound",
            "exact_bound",
            "available_exact",
            "sup_N|A_(omega,N)|<=zeta(1+2omega).",
            "Absolute convergence only; this is not the needed tail rate.",
        ),
        row(
            "bcet_11_norm_implies_discrete_energy",
            "exact_implication",
            "available_exact",
            "Uniform Burnol L2 norm implies uniform "
            "Q_(omega,N)=sum_(k>=1)k^(2omega-2)S_(omega,N)(k)^2.",
            "Uses cell-weight comparability and the uniform A_N bound.",
        ),
        row(
            "bcet_12_discrete_energy_implies_norm",
            "exact_implication",
            "available_exact",
            "Uniform Q_(omega,N) implies the uniform Burnol L2 norm.",
            "Uses |S+A*y|^2<=2S^2+2A^2 and "
            "sum k^(2omega-2)<infinity.",
        ),
        row(
            "bcet_13_exact_discrete_criterion",
            "exact_corollary",
            "available_exact",
            "sup_N||t^(-omega)f_(2omega,N)||_2<infinity iff "
            "sup_N Q_(omega,N)<infinity.",
            "Exact fixed-omega equivalence with explicit elementary constants.",
        ),
        row(
            "bcet_14_stable_prefix",
            "exact_identity",
            "available_exact",
            "For k<=N, B_(omega,N)(k)=B_(omega,infinity)(k), while "
            "S_(omega,N)(k)-S_(omega,infinity)(k)="
            "k(A_(omega,N)-A_(omega,infinity)).",
            "All divisors of n<=k are already present.",
        ),
        row(
            "bcet_15_limit_discrete_energy",
            "exact_implication",
            "available_exact",
            "Uniform Q_(omega,N) implies Q_(omega,infinity)<infinity by "
            "pointwise convergence and Fatou.",
            "One-way compactness statement only.",
        ),
        row(
            "bcet_16_reciprocal_zeta_tail_rate",
            "exact_implication",
            "available_exact",
            "Uniform Q_(omega,N)<=C forces "
            "|A_(omega,N)-1/zeta(1+2omega)|<="
            "2*sqrt(C*(1+2omega))*N^(-1/2-omega).",
            "Weighted l2 triangle inequality on the stable prefix; constants "
            "can be enlarged when starting from the continuous norm.",
        ),
        row(
            "bcet_17_zero_free_consequence",
            "exact_implication",
            "available_exact",
            "The reciprocal-zeta tail rate analytically continues "
            "sum mu(n)n^(-s)=1/zeta(s) to Re(s)>1/2+omega, so that "
            "half-plane is zero-free.",
            "Abel summation gives the continuation in the open half-plane; "
            "no boundary zero claim is made.",
        ),
        row(
            "bcet_18_open_short_interval_gate",
            "open_target",
            "open_target",
            "Prove the uniform weighted discrepancy energy Q_(omega,N) on "
            "one cofinal omega sequence; this requires both the scalar "
            "reciprocal-zeta tail rate and short-multiplicative-interval "
            "cancellation.",
            "The reduction exposes an RH-strength necessary condition but "
            "does not prove it, RH, or Lambda <= 0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_burnol_cell_energy_tail_obstruction",
        "date": "2026-07-23",
        "status": (
            "exact cell-energy/discrete-discrepancy reduction with one "
            "reciprocal-zeta tail obstruction and one open short-interval gate"
        ),
        "rows": rows,
        "refs": [
            "outputs/jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.md",
            "outputs/jensen_window_pf_jordan_muntz_causal_energy_bridge.md",
            "https://arxiv.org/abs/math/0202166",
            "https://arxiv.org/abs/math/0505453",
        ],
        "audit": {
            "row_count": 18,
            "exact_identity_or_bound_count": 10,
            "norm_reduction_step_count": 4,
            "reciprocal_zeta_tail_obstruction_count": 1,
            "zero_free_consequence_count": 1,
            "open_short_interval_gate_count": 1,
            "uniform_discrepancy_proved": False,
            "rh_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note() -> str:
    return r"""# Jensen-Window PF Burnol Cell-Energy/Tail Obstruction

Date: 2026-07-23

Status: exact cell-energy and discrete-discrepancy reduction with one
reciprocal-zeta tail obstruction and one open short-interval gate. This is
not a proof of RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_burnol_cell_energy_tail_obstruction.json
python work/rh_compute/scripts/jensen_window_pf_burnol_cell_energy_tail_obstruction.py
python work/rh_compute/scripts/check_jensen_window_pf_burnol_cell_energy_tail_obstruction.py
```

## Arithmetic Cells

Fix `0<omega<1/2` and define

```text
A_(omega,N)
 =sum_(d<=N)mu(d)d^(-(1+2omega)),

b_(omega,N)(n)
 =sum_(d|n,d<=N)mu(d)d^(-2omega),

B_(omega,N)(k)
 =sum_(n<=k)b_(omega,N)(n)
 =sum_(d<=N)mu(d)d^(-2omega)floor(k/d).             (BCET.1)
```

For Burnol's finite fractional-part function,

```text
f_(2omega,N)(t)
 =sum_(d<=N)mu(d)d^(-2omega){1/(d*t)},
```

the reciprocal coordinate `x=1/t` gives

```text
f_(2omega,N)(1/x)
 =A_(omega,N)x-B_(omega,N)(floor(x)).               (BCET.2)
```

Consequently,

```text
||t^(-omega)f_(2omega,N)||_2^2
 =sum_(k>=0) integral_k^(k+1)
  x^(2omega-2)|A_(omega,N)x-B_(omega,N)(k)|^2 dx.   (BCET.3)
```

The zero cell is exactly

```text
A_(omega,N)^2/(1+2omega).
```

For `k>=1`, write `Delta_a(k)=(k+1)^a-k^a`. Direct integration gives

```text
J_(omega,N)(k)
 =A_N^2*Delta_(1+2omega)(k)/(1+2omega)
  -2*A_N*B_N(k)*Delta_(2omega)(k)/(2omega)
  +B_N(k)^2*Delta_(2omega-1)(k)/(2omega-1).         (BCET.4)
```

Every `J_(omega,N)(k)` is a nonnegative cell integral even though the three
expanded terms need not be nonnegative separately.

## Discrete Energy

Define the summatory discrepancy

```text
S_(omega,N)(k)
 =A_(omega,N)k-B_(omega,N)(k)
 =sum_(d<=N)mu(d)d^(-2omega){k/d}.                  (BCET.5)
```

On `x=k+y`, the cell error is simply

```text
A_N*x-B_N(k)=S_N(k)+A_N*y,  0<=y<1.
```

Because `2omega-2` lies in `(-2,-1)`, the cell weight is comparable to
`k^(2omega-2)`, uniformly in `k>=1`. Also

```text
sup_N |A_(omega,N)|<=zeta(1+2omega).
```

The elementary inequalities

```text
|S+A*y|^2<=2S^2+2A^2,

S^2<=2|S+A*y|^2+2A^2*y^2
```

therefore prove the exact uniform-boundedness equivalence

```text
sup_N ||t^(-omega)f_(2omega,N)||_2<infinity

iff

sup_N Q_(omega,N)<infinity,                         (BCET.6)

Q_(omega,N)
 =sum_(k>=1)k^(2omega-2)S_(omega,N)(k)^2.
```

This is the discrete short-multiplicative-interval form of the natural
mollifier gate.

## Forced Reciprocal-Zeta Rate

Absolute convergence gives

```text
A_(omega,infinity)=1/zeta(1+2omega).
```

For every `k<=N`, all divisors of every `n<=k` already occur, so

```text
B_(omega,N)(k)=B_(omega,infinity)(k)
```

and hence

```text
S_(omega,N)(k)-S_(omega,infinity)(k)
 =k*(A_(omega,N)-A_(omega,infinity)).               (BCET.7)
```

If `Q_(omega,N)<=C` uniformly, Fatou first gives
`Q_(omega,infinity)<=C`. Taking the weighted `l2` norm of (BCET.7) on
`1<=k<=N` and using the triangle inequality gives

```text
|A_(omega,N)-1/zeta(1+2omega)|
 *sqrt(sum_(k<=N)k^(2omega))
 <=2*sqrt(C).
```

Since

```text
sum_(k<=N)k^(2omega)
 >=N^(1+2omega)/(1+2omega),
```

the uniform norm forces

```text
|A_(omega,N)-1/zeta(1+2omega)|
 <=2*sqrt(C*(1+2omega))*N^(-1/2-omega).             (BCET.8)
```

Thus the norm target already contains a square-root reciprocal-zeta tail
estimate before the remaining short-interval cancellation is considered.

## Zero-Free Consequence

Put `s_0=1+2omega` and

```text
r_N=A_(omega,infinity)-A_(omega,N).
```

The bound `r_N=O(N^(-1/2-omega))`, followed by Abel summation, makes

```text
sum_(n>=1)mu(n)n^(-s)
```

converge analytically for

```text
Re(s)>s_0-(1/2+omega)=1/2+omega.                    (BCET.9)
```

It agrees with `1/zeta(s)` where `Re(s)>1`, so analytic continuation makes
the open half-plane in (BCET.9) zero-free. No boundary-line assertion is
made. On a cofinal `omega_j->0` sequence, these open half-planes already
force RH.

This gives a direct arithmetic proof of the hard direction of the uniform
partial-norm criterion. It also shows why a generic Hilbert-space or
total-positivity estimate cannot close the problem: it would have to imply
the explicit reciprocal-zeta rate (BCET.8).

## Open Gate

The remaining target can be stated entirely discretely:

```text
sup_N sum_(k>=1) k^(2omega-2)
 |sum_(d<=N)mu(d)d^(-2omega){k/d}|^2
 <infinity                                           (BCET.10)
```

for every shift in one explicit cofinal sequence. A proof must supply both
the scalar tail rate (BCET.8) and cancellation in the weighted
short-multiplicative-interval discrepancies. Neither component is proved
here.

## Proof Boundary

Equations (BCET.1)-(BCET.9) are finite reindexings, elementary integrals,
weighted `l2` estimates, or Abel-summation consequences.
Equation (BCET.10) is open. The artifact exposes an RH-strength necessary condition;
it does not establish that condition, RH, PF-infinity, or `Lambda <= 0`.
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
        "wrote Burnol cell-energy/tail obstruction: "
        f"{len(payload['rows'])} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
