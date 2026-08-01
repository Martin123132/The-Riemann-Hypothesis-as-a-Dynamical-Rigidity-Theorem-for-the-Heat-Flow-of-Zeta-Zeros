#!/usr/bin/env python3
"""Build the Mertens anchor/logarithmic tail-energy reduction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.md"
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
            "malter_01_definitions",
            "exact_definition",
            "available_exact",
            "For 0<alpha<1 set a_n=mu(n)n^(-alpha), "
            "A_alpha(N)=sum_(n<=N)a_n, "
            "r_N=sum_(n>N)a_n/n, "
            "E_alpha=sum N^(alpha-2)|A_alpha(N)|^2, and "
            "R_alpha=sum N^alpha|r_N|^2.",
            "Here r_0=1/zeta(1+alpha); A_alpha is the weighted Mobius "
            "prefix, not the reciprocal partial sum denoted A_(omega,N) "
            "in the Burnol tail-discrepancy note.",
        ),
        row(
            "malter_02_adjacent_tail_difference",
            "exact_identity",
            "available_exact",
            "r_(N-1)-r_N=a_N/N for every N>=1.",
            "Subtract two adjacent absolutely convergent reciprocal tails.",
        ),
        row(
            "malter_03_prefix_from_all_tails",
            "exact_identity",
            "available_exact",
            "A_alpha(N)=sum_(j=0)^(N-1)r_j-N*r_N.",
            "Multiply row 2 by N and telescope. This couples all earlier "
            "cutoffs and restores the constant prefix mode.",
        ),
        row(
            "malter_04_normalized_hardy_identity",
            "exact_identity",
            "available_exact",
            "With b_N=A_alpha(N)/N and s_N=r_(N-1), "
            "b_N=H(s)(N)-r_N, where H(s)(N)=N^(-1)sum_(j=1)^N s_j.",
            "This is row 3 divided by N in weighted l2 coordinates.",
        ),
        row(
            "malter_05_tail_from_prefix",
            "exact_identity",
            "available_exact",
            "r_N=-A_alpha(N)/(N+1)+"
            "sum_(n>N)A_alpha(n)/(n(n+1)).",
            "Discrete Abel summation; A_alpha(n)/n tends to zero "
            "unconditionally for alpha>0.",
        ),
        row(
            "malter_06_weighted_hardy",
            "classical_theorem_step",
            "source_backed",
            "On l2(N^alpha), H has norm at most "
            "h_alpha=sqrt(2(3-alpha))/(1-alpha).",
            "A nonsharp explicit weighted Hardy bound obtained by the "
            "same Schur test used in the checked OU/Mertens reduction.",
        ),
        row(
            "malter_07_weighted_copson",
            "classical_theorem_step",
            "source_backed",
            "The tail operator Cb(N)=sum_(n>=N)b_n/n has norm at most "
            "c_alpha=sqrt(2(3+alpha))/(1+alpha) on l2(N^alpha).",
            "A nonsharp weighted Copson bound; the shifted kernel "
            "1/(n+1) in row 5 is pointwise dominated by this positive "
            "kernel after shifting the output index.",
        ),
        row(
            "malter_08_prefix_tail_norm_equivalence",
            "exact_equivalence",
            "available_exact",
            "E_alpha<infinity iff R_alpha<infinity. Quantitatively, "
            "R_alpha^(1/2)<=(1+c_alpha)E_alpha^(1/2), while "
            "E_alpha^(1/2)<=(1+2^(alpha/2)h_alpha)"
            "(|r_0|^2+R_alpha)^(1/2).",
            "Rows 4-7 and the one-step weight shift; r_0 is finite for "
            "every alpha>0.",
        ),
        row(
            "malter_09_stable_prefix_energy",
            "exact_identity",
            "available_exact",
            "For alpha=2omega, the stable-prefix Burnol energy is "
            "P_(omega,N)=r_N^2*sum_(k<=N)k^alpha.",
            "This is BTDR.8 with precisely the same reciprocal tail r_N.",
        ),
        row(
            "malter_10_power_sum_comparison",
            "exact_inequality",
            "available_exact",
            "(1/(1+alpha))*N^alpha*r_N^2"
            "<=P_(omega,N)/N<=N^alpha*r_N^2.",
            "Use N^(1+alpha)/(1+alpha)<=sum_(k<=N)k^alpha"
            "<=N^(1+alpha).",
        ),
        row(
            "malter_11_logarithmic_stable_prefix_equivalence",
            "exact_equivalence",
            "available_exact",
            "R_alpha<infinity iff "
            "sum_(N>=1)P_(alpha/2,N)/N<infinity.",
            "Termwise comparison in row 10; no limiting interchange or "
            "arithmetic estimate is used.",
        ),
        row(
            "malter_12_three_way_energy_equivalence",
            "exact_equivalence",
            "available_exact",
            "E_alpha<infinity iff R_alpha<infinity iff "
            "sum_(N>=1)P_(alpha/2,N)/N<infinity.",
            "Compose rows 8 and 11. The global stable-prefix family "
            "contains exactly the missing weighted-prefix anchor energy.",
        ),
        row(
            "malter_13_cofinal_rh_criterion",
            "exact_equivalence",
            "available_exact",
            "For any fixed cofinal sequence alpha_j->0 in (0,1), RH iff "
            "sum_(N>=1)P_(alpha_j/2,N)/N<infinity for every j.",
            "Compose row 12 with the already checked weighted-prefix/"
            "Mertens cofinal RH criterion. This is an equivalent target, "
            "not a proof that any series converges.",
        ),
        row(
            "malter_14_all_cutoff_anchor_resolution",
            "exact_consequence",
            "available_exact",
            "The fixed-cutoff constant anchor is not an independent mode "
            "once the complete adjacent-cutoff tail sequence is retained: "
            "row 3 reconstructs it exactly.",
            "The recovery is global in N and must not be confused with a "
            "bound from one tail value or one q-block.",
        ),
        row(
            "malter_15_fixed_cutoff_gauge_guard",
            "proof_guard",
            "guard_validated",
            "At one fixed cutoff, adding a constant to the observed prefix "
            "block still leaves its post-prefix discrepancy unchanged.",
            "The earlier constant-tail countermodel remains valid locally; "
            "the new result adds cross-cutoff compatibility rather than "
            "invalidating that guard.",
        ),
        row(
            "malter_16_critical_pointwise_consequence",
            "exact_consequence",
            "available_exact",
            "The critical pointwise rate "
            "r_N=O(N^(-(1+alpha)/2)) yields at best "
            "A_alpha(N)=O(N^((1-alpha)/2)) from row 3.",
            "Summing the critical tail envelope gives the boundary power. "
            "Its contribution to E_alpha is of harmonic size.",
        ),
        row(
            "malter_17_uniform_prefix_one_log_short",
            "proof_guard",
            "guard_validated",
            "sup_N P_(omega,N)<infinity does not imply "
            "sum_N P_(omega,N)/N<infinity.",
            "Uniform stable-prefix control supplies the critical pointwise "
            "tail rate, whereas the anchor energy requires a logarithmic "
            "mean-square gain across cutoffs.",
        ),
        row(
            "malter_18_critical_scalar_model",
            "exact_definition",
            "available_exact",
            "Let gamma=(1+alpha)/2, r_N=(N+1)^(-gamma), and "
            "a_N=N[r_(N-1)-r_N]. Then r_N=sum_(n>N)a_n/n.",
            "A positive scalar/operator model with the same adjacent-tail "
            "algebra; it is not a model of the Mobius signs.",
        ),
        row(
            "malter_19_scalar_coefficient_envelope",
            "exact_inequality",
            "available_exact",
            "The scalar model satisfies "
            "0<a_N<=gamma*N^(-gamma)<=gamma*N^(-alpha).",
            "The mean-value theorem and gamma>alpha for 0<alpha<1.",
        ),
        row(
            "malter_20_scalar_prefix_energy_bound",
            "exact_inequality",
            "available_exact",
            "For the scalar model, "
            "2^(-(1+alpha))/(1+alpha)<=P_N<=1 for every N>=1.",
            "Insert r_N=(N+1)^(-gamma) into row 9 and use row 10.",
        ),
        row(
            "malter_21_scalar_logarithmic_divergence",
            "exact_consequence",
            "available_exact",
            "For the scalar model, sum_N P_N/N diverges, and "
            "A(N)~[gamma/(1-gamma)]N^(1-gamma), so "
            "N^(alpha-2)A(N)^2 is asymptotic to a positive constant/N.",
            "This realizes the logarithmic boundary exactly.",
        ),
        row(
            "malter_22_countermodel_scope_guard",
            "proof_guard",
            "guard_validated",
            "The scalar model proves that no generic operator argument "
            "using only the critical pointwise tail rate and coefficient "
            "size can close the logarithmic series.",
            "It does not disprove a Mobius-specific cancellation theorem "
            "and makes no claim about the actual value of P_(omega,N).",
        ),
        row(
            "malter_23_full_q_distinction_guard",
            "proof_guard",
            "guard_validated",
            "The one-log obstruction concerns the stable-prefix component "
            "P alone. A cofinal uniform bound for the full Burnol energy "
            "Q would already imply RH through the existing three-gate "
            "criterion.",
            "Do not advertise the logarithmic P criterion as an extra "
            "condition that must be added after the full Q theorem.",
        ),
        row(
            "malter_24_literature_context",
            "literature_guard",
            "source_backed",
            "Burnol and Baez-Duarte supply the natural Mobius "
            "Nyman-Beurling framework; weighted Dirichlet and Hardy-space "
            "work supplies nearby zero-free and approximation criteria.",
            "Those sources motivate the coordinate but are not cited as "
            "proving rows 2-13.",
        ),
        row(
            "malter_25_focused_source_search_guard",
            "literature_guard",
            "source_backed",
            "A focused search of the cited primary sources found no exact "
            "displayed match for sum_N P_(omega,N)/N and the three-way "
            "equivalence in row 12.",
            "This is a search report, not a novelty or priority claim; a "
            "referee-level literature review remains appropriate.",
        ),
        row(
            "malter_26_open_logarithmic_gain_gate",
            "open_target",
            "open_target",
            "Prove for the actual Mobius tail, on a fixed cofinal sequence "
            "alpha_j->0, the logarithmic estimate "
            "sum_N P_(alpha_j/2,N)/N<infinity without assuming RH.",
            "Equivalently prove the weighted-prefix/Mertens energy. The "
            "reduction supplies no arithmetic cancellation and proves "
            "neither RH, PF-infinity, nor Lambda<=0.",
        ),
    ]
    return {
        "kind": (
            "jensen_window_pf_mertens_anchor_"
            "logarithmic_tail_energy_reduction"
        ),
        "date": "2026-07-23",
        "status": (
            "exact all-cutoff anchor/logarithmic stable-prefix reduction "
            "with one open RH-equivalent arithmetic gate"
        ),
        "parameters": {
            "alpha_range": "0<alpha<1",
            "omega_match": "alpha=2omega",
            "gamma": "(1+alpha)/2",
        },
        "source_anchors": [
            "https://arxiv.org/abs/math/0202166",
            "https://arxiv.org/abs/math/0011254",
            "https://arxiv.org/abs/2508.00388",
            "https://doi.org/10.1007/s11785-025-01661-2",
            "https://arxiv.org/abs/2606.16097",
        ],
        "rows": rows,
        "audit": {
            "row_count": 26,
            "exact_reduction_count": 17,
            "classical_theorem_step_count": 2,
            "literature_guard_count": 2,
            "proof_guard_count": 4,
            "open_logarithmic_gain_gate_count": 1,
            "prefix_tail_norm_equivalence_proved": True,
            "logarithmic_stable_prefix_equivalence_proved": True,
            "cofinal_rh_equivalence_proved": True,
            "uniform_stable_prefix_closes_log_gate": False,
            "logarithmic_gain_proved": False,
            "full_burnol_bound_proved": False,
            "rh_proved": False,
            "pf_infinity_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note() -> str:
    return r"""# Jensen-Window PF Mertens Anchor/Logarithmic Tail-Energy Reduction

Date: 2026-07-23

Status: exact all-cutoff anchor/logarithmic stable-prefix reduction with
one open RH-equivalent arithmetic gate. This is not a proof of RH,
PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.json
python work/rh_compute/scripts/jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.py
```

## Adjacent Tails Recover The Anchor

Fix `0<alpha<1` and write

```text
a_n:=mu(n)n^(-alpha),
A_alpha(N):=sum_(n<=N)a_n,

r_N:=sum_(n>N)a_n/n,
r_0=1/zeta(1+alpha).                                  (MALTER.1)
```

Here `A_alpha` is the weighted Mobius prefix. It is not the reciprocal
partial sum called `A_(omega,N)` in the Burnol tail-discrepancy note.

Adjacent tails satisfy

```text
r_(N-1)-r_N=a_N/N.                                    (MALTER.2)
```

Multiplying by `N` and telescoping gives the exact reconstruction

```text
A_alpha(N)
 =sum_(j=0)^(N-1)r_j-N*r_N.                           (MALTER.3)
```

Thus the constant anchor that is invisible in one post-prefix q-block is
not independent once the compatible tail sequence at every cutoff is
retained. The recovery is global in `N`; one tail value still does not
determine the anchor.

The inverse identity is the previously checked Abel transform

```text
r_N
 =-A_alpha(N)/(N+1)
  +sum_(n>N)A_alpha(n)/(n(n+1)).                       (MALTER.4)
```

No zero-free region is needed for either identity.

## Weighted Hardy/Copson Equivalence

Set

```text
E_alpha:=sum_(N>=1)N^(alpha-2)|A_alpha(N)|^2,
R_alpha:=sum_(N>=1)N^alpha|r_N|^2,
b_N:=A_alpha(N)/N,
s_N:=r_(N-1).
```

Then (MALTER.3) says

```text
b_N=H(s)(N)-r_N,
H(s)(N):=N^(-1)sum_(j=1)^N s_j.                       (MALTER.5)
```

On `l2(N^alpha)`, the weighted Hardy and Copson operators admit the
nonsharp Schur bounds

```text
||H||<=h_alpha:=sqrt(2(3-alpha))/(1-alpha),
||C||<=c_alpha:=sqrt(2(3+alpha))/(1+alpha),            (MALTER.6)

(Cb)(N):=sum_(n>=N)b_n/n.
```

The one-step weight shifts and (MALTER.4)-(MALTER.6) give

```text
R_alpha^(1/2)
 <=(1+c_alpha)E_alpha^(1/2),

E_alpha^(1/2)
 <=(1+2^(alpha/2)h_alpha)
   (|r_0|^2+R_alpha)^(1/2).                            (MALTER.7)
```

Since `r_0` is finite for `alpha>0`,

```text
E_alpha<infinity
 iff
R_alpha<infinity.                                     (MALTER.8)
```

This is an exact finiteness equivalence with explicit nonsharp operator
constants. It does not estimate either arithmetic sequence.

## Stable-Prefix Energy Across All Cutoffs

Put `alpha=2omega`. The stable-prefix term in the Burnol
tail-discrepancy split is

```text
P_(omega,N)
 =r_N^2*sum_(k<=N)k^alpha.                             (MALTER.9)
```

The elementary power-sum comparison yields

```text
(1/(1+alpha))*N^alpha*r_N^2
 <=P_(omega,N)/N
 <=N^alpha*r_N^2.                                     (MALTER.10)
```

Therefore

```text
E_alpha<infinity
 iff
R_alpha<infinity
 iff
sum_(N>=1)P_(alpha/2,N)/N<infinity.                    (MALTER.11)
```

This is the precise all-cutoff form of the weighted-prefix anchor. A
uniform bound on each `P_(omega,N)` is one logarithm weaker than the
required series.

Composing (MALTER.11) with the checked weighted-prefix/Mertens criterion
gives, for every fixed cofinal sequence `alpha_j->0`,

```text
RH
 iff
sum_(N>=1)P_(alpha_j/2,N)/N<infinity
for every j.                                           (MALTER.12)
```

This is an RH-equivalent target, not a proof of convergence.

## The Exact Logarithmic Barrier

The pointwise stable-prefix bound gives only

```text
r_N=O(N^(-(1+alpha)/2)).
```

Inserting this envelope into (MALTER.3) gives at best

```text
A_alpha(N)=O(N^((1-alpha)/2)),
```

whose contribution to `E_alpha` is of order `1/N`. The loss is genuinely
logarithmic, as the following scalar model shows.

Let

```text
gamma:=(1+alpha)/2,
r_N:=(N+1)^(-gamma),
a_N:=N[r_(N-1)-r_N].                                  (MALTER.13)
```

Then

```text
r_N=sum_(n>N)a_n/n,
0<a_N<=gamma*N^(-gamma)<=gamma*N^(-alpha).             (MALTER.14)
```

The associated stable-prefix energy obeys

```text
2^(-(1+alpha))/(1+alpha)<=P_N<=1,                      (MALTER.15)
```

so `sum_N P_N/N` diverges. Moreover,

```text
A(N)
 ~[gamma/(1-gamma)]N^(1-gamma),

N^(alpha-2)A(N)^2
 ~[gamma/(1-gamma)]^2/N.                               (MALTER.16)
```

This is a scalar/operator countermodel, not a surrogate for `mu`. It
proves that coefficient size plus the critical pointwise tail rate cannot
close the logarithmic series by generic analysis. It leaves open a
Mobius-specific mean-square cancellation theorem.

## Relation To The Earlier Local Guard

The earlier fixed-cutoff constant-tail countermodel remains correct. At
one cutoff the q-block observation is invariant under a constant prefix
shift. Equations (MALTER.2)-(MALTER.3) add the missing compatibility
between adjacent cutoffs. In short:

```text
one cutoff:       anchor is invisible;
all cutoff tails: anchor is reconstructed;
critical sup P:   reconstruction is one logarithm nonsummable;
sum P_N/N:        exact weighted-prefix/Mertens energy.
```

The obstruction concerns the stable-prefix component `P` alone. A
cofinal uniform theorem for the full Burnol energy `Q` would already
imply RH through the existing three-gate criterion; (MALTER.12) is a
cross-coordinate sharpening, not an extra condition after full `Q`.

## Literature Guard

The primary context checked was:

- Burnol's analytic estimate and strengthened Nyman-Beurling argument:
  https://arxiv.org/abs/math/0202166
- Baez-Duarte's arithmetical Nyman-Beurling reformulation:
  https://arxiv.org/abs/math/0011254
- weighted Hardy/Copson context:
  https://arxiv.org/abs/2508.00388
- weighted Dirichlet-space zero-free criteria:
  https://doi.org/10.1007/s11785-025-01661-2
- the current Hardy-space approximation formulation:
  https://arxiv.org/abs/2606.16097

A focused search found nearby Mobius approximation and zero-free
criteria, but no exact displayed match for (MALTER.11). That is a search
report, not a novelty or priority claim.

## Open Gate

The sharpened arithmetic target is:

```text
for every alpha in one fixed cofinal sequence alpha_j->0,

sum_(N>=1) P_(alpha/2,N)/N < infinity                  (MALTER.17)
```

without RH or a zero-free hypothesis. Equivalently, prove the
weighted-prefix/Mertens energy. The present reduction supplies no Mobius
cancellation and proves neither RH, PF-infinity, nor `Lambda <= 0`.
"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    args.result.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    payload = build_payload()
    args.result.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(), encoding="utf-8")
    print(
        "wrote Mertens anchor/logarithmic tail-energy reduction: "
        f"{len(payload['rows'])} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
