#!/usr/bin/env python3
"""Build the planar boundary suffix-energy handoff."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.md"
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


def build_rows() -> list[dict[str, str]]:
    return [
        row(
            "pbseh_01_inherited_data",
            "exact_definition",
            "available_exact",
            (
                "Retain N=2K+1, q_r=2r-1, the truncated kernel "
                "P_(K,R), current coefficients u_i=mu(K+i), future "
                "scalars beta and Q, endpoint column p_i=P(i,K), and "
                "c=P(K,K) from Corollary 11.22Z.8."
            ),
            "The range is 1<=R<=K; no Mobius estimate is assumed.",
        ),
        row(
            "pbseh_02_suffixes",
            "exact_definition",
            "available_exact",
            (
                "Set A_i=sum_(j=i)^(K-1)u_j and "
                "H_R^odd=sum_(r=1)^R 1/q_r."
            ),
            "A_i is an anchored current-block Mobius suffix.",
        ),
        row(
            "pbseh_03_boundary_weights",
            "exact_definition",
            "available_exact",
            (
                "Put w_0=0, w_i=P(i,K-1) for 1<=i<=K-1, "
                "w_K=c, and delta_i=w_i-p_i for 1<=i<K."
            ),
            "The final value w_K is defined as c, not P(K,K-1).",
        ),
        row(
            "pbseh_04_shifted_prefix_relation",
            "exact_identity",
            "available_exact",
            (
                "For 2<=i<=K-2, S_i^bd=S_(i-1)^sh+"
                "u_i A_(i+1), with the same formula at i=1 after "
                "setting S_0^sh=0."
            ),
            "This is a finite reindexing of current-current edges.",
        ),
        row(
            "pbseh_05_boundary_split",
            "exact_identity",
            "available_exact",
            (
                "B_CC-A_CF=sum_(i=1)^(K-2)w_i u_i A_(i+1)+"
                "sum_(i=2)^(K-2)delta_i S_(i-1)^sh."
            ),
            "The boundary and shoulder are joined before the split.",
        ),
        row(
            "pbseh_06_edge_kernel",
            "exact_identity",
            "available_exact",
            (
                "An edge (a,a+g) has joined boundary kernel "
                "P(a,K-1)+sum_(i=a+1)^(K-1-g)"
                "[P(i,K-1)-P(i,K)]."
            ),
            "The sum is empty when a+g=K-1.",
        ),
        row(
            "pbseh_07_mode_edge_kernel",
            "exact_identity",
            "available_exact",
            (
                "Modewise the edge kernel numerator is "
                "phi_r(K-1)phi_r(a)+[phi_r(K-1)-e_r]"
                "sum_(i=a+1)^(K-1-g)phi_r(i)."
            ),
            "This is the cancellation-preserving form of B_CC-A_CF.",
        ),
        row(
            "pbseh_08_endpoint_delta",
            "exact_identity",
            "available_exact",
            (
                "delta_i=-(pi^2/N^2)sum_(r<=R)"
                "phi_r(i)e_r sinc^2(theta_r/2)."
            ),
            "The endpoint jump cancels q_r^2 exactly.",
        ),
        row(
            "pbseh_09_delta_interval_bound",
            "exact_bound",
            "available_exact",
            (
                "Every consecutive interval I satisfies "
                "|sum_(i in I)delta_i|<="
                "8pi H_R^odd/N^2."
            ),
            "This uses the finite sine-sum formula and no Mobius input.",
        ),
        row(
            "pbseh_10_shoulder_correction_bound",
            "exact_bound",
            "available_exact",
            (
                "|sum_(i=2)^(K-2)delta_i S_(i-1)^sh|"
                "<pi H_R^odd for |u_i|<=1."
            ),
            "The proof swaps to current edges before taking absolute values.",
        ),
        row(
            "pbseh_11_future_scalar_bound",
            "exact_bound",
            "available_exact",
            "|beta|<=K+1/2=N/2.",
            "Only |mu|<=1 and 0<=b_m<=(4K-m)/(2K) are used.",
        ),
        row(
            "pbseh_12_anchor_delta_bound",
            "exact_bound",
            "available_exact",
            (
                "|beta sum_(i=1)^(K-1)u_i delta_i|"
                "<pi^2/2."
            ),
            "The pointwise delta bound is 4pi^2 R/N^3.",
        ),
        row(
            "pbseh_13_positive_increments",
            "exact_identity",
            "available_exact",
            (
                "w_i-w_(i-1)>0 for 1<=i<=K, where w_K=c."
            ),
            "Strict positivity follows from S_R(u)>0 on 0<u<pi.",
        ),
        row(
            "pbseh_14_increment_bounds",
            "exact_bound",
            "available_exact",
            (
                "With C_sin=1/2+2pi, the uniform odd sine-sum bound "
                "0<S_R(u)<=C_sin gives "
                "0<w_i-w_(i-1)<=4pi C_sin/N^2 for every 1<=i<=K."
            ),
            "The constant is independent of K and R.",
        ),
        row(
            "pbseh_15_extended_energy",
            "exact_definition",
            "available_exact",
            (
                "E_B=sum_(i=1)^(K-1)(w_i-w_(i-1))"
                "(A_i+beta)^2+(c-w_(K-1))beta^2."
            ),
            "Every coefficient is positive, so E_B is nonnegative.",
        ),
        row(
            "pbseh_16_extended_diagonal",
            "exact_definition",
            "available_exact",
            (
                "D_B=sum_(i=1)^(K-1)w_i u_i^2+cQ and "
                "C_delta=sum_(i=2)^(K-2)delta_i S_(i-1)^sh-"
                "beta sum_(i=1)^(K-1)u_i delta_i."
            ),
            "C_delta retains the endpoint-smoothed shoulder and anchor terms.",
        ),
        row(
            "pbseh_17_min_kernel_gram",
            "exact_identity",
            "available_exact",
            (
                "sum_i w_i u_i A_(i+1)+beta sum_i w_i u_i"
                "=(1/2)[E_B-sum_i w_i u_i^2-c beta^2]."
            ),
            "This is finite suffix Abel summation for the min(w_i,w_j) kernel.",
        ),
        row(
            "pbseh_18_joined_remainder",
            "exact_identity",
            "available_exact",
            (
                "R_B:=B_CC-A_CF+beta gamma+(c/2)(beta^2-Q)"
                "=(1/2)(E_B-D_B)+C_delta."
            ),
            "This is the complete joined boundary/future remainder.",
        ),
        row(
            "pbseh_19_diagonal_bound",
            "exact_bound",
            "available_exact",
            "0<=D_B<5pi^2/12.",
            (
                "The current part is below pi^2/4 and the future "
                "part cQ is below pi^2/6."
            ),
        ),
        row(
            "pbseh_20_logarithmic_error",
            "exact_bound",
            "available_exact",
            (
                "|R_B-(1/2)E_B|<="
                "pi H_R^odd+17pi^2/24."
            ),
            "The right side is O(1+log R) unconditionally.",
        ),
        row(
            "pbseh_21_subpower_equivalence",
            "exact_reduction",
            "available_exact",
            (
                "R_B=O_epsilon(K^epsilon) for every epsilon>0 "
                "if and only if E_B=O_epsilon(K^epsilon) for every "
                "epsilon>0."
            ),
            "Logarithms are absorbed at the all-epsilon subpower scale.",
        ),
        row(
            "pbseh_22_block_mean_square_sufficient",
            "exact_reduction",
            "available_exact",
            (
                "The earlier affine-tent energy is exactly "
                "H_K=sum_(i=1)^(K-1)(A_i+beta)^2+beta^2, and "
                "E_B<=pi(1/2+2pi)K^(-2)H_K."
            ),
            (
                "Thus its dyadic criterion dominates E_B, but no "
                "Mobius estimate for H_K is proved here."
            ),
        ),
        row(
            "pbseh_23_direct_route",
            "exact_reduction",
            "available_exact",
            (
                "Together with the inherited Y_(K,r) estimate for "
                "O_(CC,int), the E_B subpower estimate controls the "
                "complete off-diagonal O_(alpha,K); the prior "
                "affine-tent criterion is sufficient for its "
                "boundary/future part."
            ),
            "Both arithmetic estimates remain open.",
        ),
        row(
            "pbseh_24_mobius_finite_audit",
            "finite_validation",
            "validated_finite",
            (
                "At alpha=1/2 and K=16,...,1024, direct Mobius "
                "audits reproduce R_B=(E_B-D_B)/2+C_delta; at "
                "K=1024, E_B=0.097890613, D_B=1.325896381, and "
                "C_delta=-0.000050764."
            ),
            "Finite values prove no asymptotic estimate.",
        ),
        row(
            "pbseh_25_open_energy_gate",
            "theorem_target",
            "open",
            (
                "Prove E_B=O_epsilon(K^epsilon) on every dyadic "
                "scale for each member of one fixed cofinal "
                "positive-alpha sequence."
            ),
            "This is an anchored weighted Mertens-energy obligation.",
        ),
        row(
            "pbseh_26_proof_boundary",
            "proof_guard",
            "guard_active",
            (
                "The suffix reindexing, positive energy identity, "
                "and logarithmic error bound are exact."
            ),
            (
                "No E_B estimate, Y projection estimate, reciprocal-tail "
                "theorem, RH, PF-infinity, Lambda<=0, or Clay-prize result "
                "is proved."
            ),
        ),
    ]


MOBIUS_AUDIT_ROWS = [
    {
        "K": 16,
        "R": 8,
        "E_B": 0.110319464,
        "D_B": 1.292932999,
        "C_delta": -0.001824749,
        "R_B": -0.593131517,
    },
    {
        "K": 32,
        "R": 14,
        "E_B": 0.083247421,
        "D_B": 1.279015568,
        "C_delta": 0.001842047,
        "R_B": -0.596042026,
    },
    {
        "K": 64,
        "R": 23,
        "E_B": 0.135323930,
        "D_B": 1.285158679,
        "C_delta": 0.000036395,
        "R_B": -0.574880980,
    },
    {
        "K": 128,
        "R": 39,
        "E_B": 0.147238007,
        "D_B": 1.319393967,
        "C_delta": 0.000116640,
        "R_B": -0.585961340,
    },
    {
        "K": 256,
        "R": 64,
        "E_B": 0.115185299,
        "D_B": 1.326790649,
        "C_delta": 0.000251649,
        "R_B": -0.605551026,
    },
    {
        "K": 512,
        "R": 108,
        "E_B": 0.111068291,
        "D_B": 1.324491974,
        "C_delta": 0.000006107,
        "R_B": -0.606705735,
    },
    {
        "K": 1024,
        "R": 182,
        "E_B": 0.097890613,
        "D_B": 1.325896381,
        "C_delta": -0.000050764,
        "R_B": -0.614053648,
    },
]


NOTE_HEAD = r"""# Jensen-Window PF Mertens Planar Boundary Suffix-Energy Handoff

Date: 2026-07-24

Status: exact positive suffix-energy reduction for the complete joined
boundary/future remainder, with one open anchored Mobius-energy gate.
This is not a proof of RH, PF-infinity, `Lambda<=0`, or a Clay-prize
result.

```text
work/rh_compute/results/
jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.json
python work/rh_compute/scripts/
jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.py
python work/rh_compute/scripts/
check_jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.py
```

## Boundary And Shoulder As One Edge Kernel

Retain all notation from Corollary 11.22Z.8. In particular,

```text
u_i=mu(K+i), 1<=i<K,
p_i=P_(K,R)(i,K), c=P_(K,R)(K,K),

B_CC=sum_(i=1)^(K-2)S_i^bd P_(K,R)(i,K-1),
A_CF=sum_(i=1)^(K-3)S_i^sh p_(i+1).
```

Put

```text
A_i:=sum_(j=i)^(K-1)u_j,
w_0:=0,
w_i:=P_(K,R)(i,K-1), 1<=i<=K-1,
w_K:=c,
delta_i:=w_i-p_i,
H_R^odd:=sum_(r=1)^R 1/(2r-1).
```

Direct reindexing gives

```text
S_i^bd=S_(i-1)^sh+u_i A_(i+1)                 (PBSE.1)
```

with `S_0^sh=0`. Hence the mandatory joined boundary is

```text
B_CC-A_CF
 =sum_(i=1)^(K-2)w_i u_i A_(i+1)
  +sum_(i=2)^(K-2)delta_i S_(i-1)^sh.         (PBSE.2)
```

Equivalently, an edge `(a,a+g)` has kernel

```text
P_(K,R)(a,K-1)
 +sum_(i=a+1)^(K-1-g)
   [P_(K,R)(i,K-1)-P_(K,R)(i,K)].             (PBSE.3)
```

Modewise its numerator is

```text
phi_r(K-1)phi_r(a)
 +[phi_r(K-1)-e_r]
   sum_(i=a+1)^(K-1-g)phi_r(i).                (PBSE.4)
```

This is the cancellation-preserving arithmetic representation of
`phi_r(K-1)U_r-e_rV_r`.

## The Endpoint Correction Is Harmless

The exact endpoint jump from Corollary 11.22Z.8 gives

```text
delta_i
 =-(pi^2/N^2)sum_(r=1)^R
   phi_r(i)e_r sinc^2(theta_r/2).              (PBSE.5)
```

For every consecutive interval `I`, the finite sine-sum formula yields

```text
|sum_(i in I)delta_i|
 <=8pi H_R^odd/N^2.                            (PBSE.6)
```

Swap the first correction in (PBSE.2) back to its current edges before
taking absolute values. There are fewer than `K^2/2` edges, so

```text
|sum_(i=2)^(K-2)delta_i S_(i-1)^sh|
 <pi H_R^odd.                                  (PBSE.7)
```

Also,

```text
|delta_i|<=4pi^2 R/N^3,
|beta|<=sum_(m=2K)^(4K-1)(4K-m)/(2K)=N/2,

|beta sum_(i=1)^(K-1)u_i delta_i|<pi^2/2.      (PBSE.8)
```

Thus neither endpoint-smoothed correction carries a power of `K`.

## Positive Anchored Suffix Energy

The positive sine polynomial representation of `P_(K,R)` proves

```text
w_i-w_(i-1)>0, 1<=i<=K.                       (PBSE.9)
```

For `i<K`, direct spectral differencing gives

```text
w_i-w_(i-1)<=4pi H_R^odd/N^2,
w_K-w_(K-1)<4pi^2/N^2.                        (PBSE.10)
```

There is also a rank-uniform bound. By symmetry reduce to
`0<u<=pi/2`. If `u<1/2`, split the odd sine series where
`(2r-1)u<=1` and apply Abel summation to the tail; consecutive odd-sine
partial sums are bounded by `1/sin u`. If `u>=1/2`, apply the same
Dirichlet bound from the first term. This gives the explicit estimate

```text
0<S_R(u)<=C_sin:=1/2+2pi,             0<u<pi.
```

Using the positive integral representation of `P_(K,R)` on the two
new endpoint slivers then gives

```text
0<w_i-w_(i-1)<=4pi C_sin/N^2,          1<=i<=K.
                                                        (PBSE.10a)
```

Define

```text
E_B
 :=sum_(i=1)^(K-1)(w_i-w_(i-1))(A_i+beta)^2
   +(c-w_(K-1))beta^2,                        (PBSE.11)

D_B
 :=sum_(i=1)^(K-1)w_i u_i^2+cQ,

C_delta
 :=sum_(i=2)^(K-2)delta_i S_(i-1)^sh
   -beta sum_(i=1)^(K-1)u_i delta_i.
```

Finite suffix Abel summation for the min kernel generated by `w_i`
gives

```text
sum_i w_i u_i A_(i+1)+beta sum_i w_i u_i
 =(1/2)[E_B-sum_i w_i u_i^2-c beta^2].        (PBSE.12)
```

Since `gamma=sum_i u_i p_i` and the future-future term is
`(c/2)(beta^2-Q)`, the complete joined remainder is exactly

```text
R_B
 :=B_CC-A_CF+beta gamma+(c/2)(beta^2-Q)
 =(1/2)(E_B-D_B)+C_delta.                     (PBSE.13)
```

No mixed cell remains outside this formula.

## Lossless Subpower Gate

Monotonicity and the completed-kernel diagonal bound imply

```text
sum_i w_i u_i^2<pi^2/4,
cQ<pi^2/6,
0<=D_B<5pi^2/12.                              (PBSE.14)
```

Combining (PBSE.7), (PBSE.8), and (PBSE.14) gives the unconditional
two-sided approximation

```text
|R_B-(1/2)E_B|
 <=pi H_R^odd+17pi^2/24
 =O(1+log R).                                  (PBSE.15)
```

Consequently

```text
R_B=O_epsilon(K^epsilon) for every epsilon>0

if and only if

E_B=O_epsilon(K^epsilon) for every epsilon>0. (PBSE.16)
```

This is a genuine reduction because `E_B>=0`. It also meets an earlier
exact coordinate. The ordinary-Mobius affine-tent handoff defines

```text
H_K
 :=sum_(i=1)^(K-1)(A_i+beta)^2+beta^2.
```

This is exactly (MATBH.40), including the `i=K` term `beta^2`.
Equation (PBSE.10a) gives the unconditional domination

```text
E_B
 <=4pi C_sin H_K/N^2
 <pi C_sin K^(-2)H_K.                         (PBSE.17)
```

Consequently

```text
sum_(K dyadic)K^(-alpha)E_B
 <=pi C_sin
   sum_(K dyadic)K^(-2-alpha)H_K.              (PBSE.18)
```

The right side is the exact affine-tent criterion (MATBH.39). Thus that
already identified RH-equivalent open energy controls the complete
boundary/future part without an additional logarithmic weight. It is
still not proved. Together with the inherited `Y_(K,r)` square-function
estimate for `O_(CC,int)`, (PBSE.16) controls the complete signed
off-diagonal.

## Finite Mobius Diagnostic

At `alpha=1/2`, the exact identity has the following values:

| K | R | E_B | D_B | C_delta | R_B |
|---:|---:|---:|---:|---:|---:|
"""


NOTE_TAIL = r"""

The endpoint correction is much smaller than its unconditional
logarithmic envelope on this finite grid. That observation is
diagnostic only.

## Proof Boundary

The boundary/shoulder reindexing, endpoint correction bounds, positive
suffix-energy identity, diagonal bound, and logarithmic approximation
are exact and independently checked.

No estimate for `H_K` or `E_B`, no `Y_(K,r)` projection estimate, no
full reciprocal-tail theorem, no RH, no PF-infinity, no `Lambda<=0`,
and no Clay-prize result is proved.

## Claim Ledger

| ID | Role | Status | Statement | Boundary |
|---|---|---|---|---|
"""


def build_payload() -> dict:
    rows = build_rows()
    return {
        "kind": (
            "jensen_window_pf_mertens_planar_"
            "boundary_suffix_energy_handoff"
        ),
        "date": "2026-07-24",
        "status": (
            "exact positive suffix-energy reduction for the complete "
            "joined boundary/future remainder with one open anchored "
            "Mobius-energy gate"
        ),
        "source": {
            "parent_handoff": (
                "outputs/jensen_window_pf_mertens_planar_"
                "boundary_flux_anchor_handoff.md"
            ),
            "affine_tent_handoff": (
                "outputs/jensen_window_pf_mertens_"
                "affine_tent_bridge_handoff.md"
            ),
            "truncated_kernel": (
                "outputs/jensen_window_pf_mertens_"
                "truncated_half_odd_kernel_handoff.md"
            ),
            "formal_core": "outputs/formal_core.md",
        },
        "rows": rows,
        "mobius_audit": {
            "alpha": 0.5,
            "rows": MOBIUS_AUDIT_ROWS,
            "interpretation": (
                "deterministic finite identity checks only; no "
                "asymptotic estimate"
            ),
        },
        "summary": {
            "row_count": len(rows),
            "exact_reduction_count": 23,
            "conditional_calibration_count": 0,
            "proof_guard_count": 1,
            "finite_validation_count": 1,
            "open_energy_gate_count": 1,
            "boundary_suffix_split_proved": True,
            "endpoint_correction_logarithmic": True,
            "positive_suffix_energy_proved": True,
            "subpower_equivalence_proved": True,
            "anchored_mobius_energy_proved": False,
            "interior_projection_estimate_proved": False,
            "rh_proved": False,
            "pf_infinity_proved": False,
            "lambda_le_zero_proved": False,
        },
        "next_action": (
            "Seek a cancellation theorem for the anchored block suffixes "
            "A_i+beta under the positive increments w_i-w_(i-1), and "
            "couple it to the inherited Y_(K,r) projection without "
            "splitting the boundary/shoulder join."
        ),
        "proof_boundary": (
            "Exact finite algebra and unconditional logarithmic correction "
            "bound only; no anchored Mobius-energy estimate or global "
            "Riemann-hypothesis conclusion."
        ),
    }


def render_note(payload: dict) -> str:
    lines = [NOTE_HEAD.rstrip()]
    for item in payload["mobius_audit"]["rows"]:
        lines.append(
            "| {K} | {R} | {E_B:.9g} | {D_B:.9g} | "
            "{C_delta:.9g} | {R_B:.9g} |".format(**item)
        )
    lines.append(NOTE_TAIL.rstrip())
    for item in payload["rows"]:
        lines.append(
            "| `{id}` | `{role}` | `{status}` | {statement} | "
            "{proof_boundary} |".format(**item)
        )
    summary = payload["summary"]
    lines.extend(
        [
            "",
            "Summary:",
            "",
            f"- rows: {summary['row_count']}",
            f"- exact reductions: {summary['exact_reduction_count']}",
            (
                "- conditional calibrations: "
                f"{summary['conditional_calibration_count']}"
            ),
            f"- proof guards: {summary['proof_guard_count']}",
            f"- finite validations: {summary['finite_validation_count']}",
            (
                "- open anchored energy gates: "
                f"{summary['open_energy_gate_count']}"
            ),
            "",
        ]
    )
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote planar boundary suffix-energy handoff: "
        f"{payload['summary']['row_count']} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
