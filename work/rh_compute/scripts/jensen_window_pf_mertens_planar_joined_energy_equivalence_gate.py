#!/usr/bin/env python3
"""Build the planar joined-energy equivalence route guard."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.md"
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
            "pjeg_01_inherited_data",
            "exact_definition",
            "available_exact",
            (
                "Retain O_(CC,int), E_B, D_B, C_delta, x_r, e_r, "
                "beta, Q, c, and D_(C,r) from Corollaries 11.22Z.8-Z.10."
            ),
            "No arithmetic estimate is assumed.",
        ),
        row(
            "pjeg_02_lossless_energy",
            "exact_definition",
            "available_exact",
            (
                "Set L=sum_(r<=R)(x_r+beta e_r)^2/q_r^2="
                "E_perp+c(beta+gamma/c)^2."
            ),
            "L is nonnegative and is the inherited low-mode energy.",
        ),
        row(
            "pjeg_03_lossless_diagonal",
            "exact_definition",
            "available_exact",
            (
                "Set D_L=sum_(r<=R)D_(C,r)/q_r^2+cQ."
            ),
            "D_L is the current/future diagonal in the same mode coordinates.",
        ),
        row(
            "pjeg_04_block_compression",
            "exact_identity",
            "available_exact",
            "O_(alpha,K)=(L-D_L)/2.",
            "This is the exact three-block compression, not an estimate.",
        ),
        row(
            "pjeg_05_lossless_diagonal_bound",
            "exact_bound",
            "available_exact",
            "0<=D_L<7pi^2/24.",
            "The bound is uniform in alpha, K, and R.",
        ),
        row(
            "pjeg_06_offdiagonal_lossless_comparison",
            "exact_bound",
            "available_exact",
            "|O_(alpha,K)-L/2|=D_L/2<7pi^2/48.",
            "The full signed off-diagonal and positive energy differ boundedly.",
        ),
        row(
            "pjeg_07_suffix_remainder",
            "exact_identity",
            "available_exact",
            (
                "R_B=(E_B-D_B)/2+C_delta and "
                "O_(alpha,K)=O_(CC,int)+R_B."
            ),
            "The boundary shoulder remains included in R_B.",
        ),
        row(
            "pjeg_08_defect_bound",
            "exact_bound",
            "available_exact",
            (
                "|C_delta|<=pi H_R^odd+pi^2/2."
            ),
            "The two terms are the shoulder and anchor endpoint corrections.",
        ),
        row(
            "pjeg_09_joined_candidate",
            "exact_definition",
            "available_exact",
            "J_B:=O_(CC,int)+E_B/2.",
            "This is the proposed signed current-interior/boundary join.",
        ),
        row(
            "pjeg_10_joined_identity",
            "exact_identity",
            "available_exact",
            (
                "J_B=L/2+(D_B-D_L)/2-C_delta."
            ),
            "No Cauchy-Schwarz or absolute cell split is used.",
        ),
        row(
            "pjeg_11_joined_error_bound",
            "exact_bound",
            "available_exact",
            (
                "|J_B-L/2|<=pi H_R^odd+41pi^2/48."
            ),
            "Use D_B<5pi^2/12 and D_L<7pi^2/24.",
        ),
        row(
            "pjeg_12_joined_offdiagonal_bound",
            "exact_bound",
            "available_exact",
            (
                "|J_B-O_(alpha,K)|<=pi H_R^odd+17pi^2/24."
            ),
            "This is the same summable coordinate defect as in PBSE.15.",
        ),
        row(
            "pjeg_13_pointwise_subpower_equivalence",
            "exact_equivalence",
            "available_exact",
            (
                "For polynomially bounded R, J_B, O_(alpha,K), and L "
                "are all-epsilon subpower simultaneously."
            ),
            "The logarithmic defects are absorbed for every positive epsilon.",
        ),
        row(
            "pjeg_14_weighted_dyadic_equivalence",
            "exact_equivalence",
            "available_exact",
            (
                "For every fixed alpha>0, sum K^(-alpha)L is finite "
                "iff sum K^(-alpha)|J_B| is finite, iff "
                "sum K^(-alpha)|O_(alpha,K)| is finite."
            ),
            "The K^(-alpha) logarithmic and bounded defects are summable.",
        ),
        row(
            "pjeg_15_cofinal_composition",
            "theorem_composition",
            "available_exact",
            (
                "On one fixed cofinal positive-alpha sequence, the "
                "joined J_B criterion is the prior lossless RH criterion."
            ),
            "This is a route identification, not a proof of its antecedent.",
        ),
        row(
            "pjeg_16_route_collapse",
            "route_guard",
            "guard_active",
            (
                "Retaining O_(CC,int) with E_B/2 does not create a new "
                "weaker theorem interface; it returns L/2 up to a "
                "summable coordinate defect."
            ),
            "Separate current-interior/E_B work may not be advertised as easier.",
        ),
        row(
            "pjeg_17_finite_mobius_audit",
            "finite_validation",
            "validated_finite",
            (
                "At alpha=1/2 and K=16,...,1024, direct Mobius "
                "recomputation verifies the joined identity to float64 "
                "roundoff."
            ),
            "Finite values prove no asymptotic estimate.",
        ),
        row(
            "pjeg_18_open_lossless_gate",
            "theorem_target",
            "open",
            (
                "Prove dyadic summability of K^(-alpha)L for every "
                "member of one fixed cofinal positive-alpha sequence."
            ),
            "This remains the RH-equivalent low-mode arithmetic obligation.",
        ),
        row(
            "pjeg_19_next_route",
            "route_guard",
            "guard_active",
            (
                "Future work must attack L, E_B, or the strong signed "
                "edge form honestly at RH strength, or find structure "
                "outside this algebraically equivalent join."
            ),
            "Positivity and coordinate changes alone cannot supply cancellation.",
        ),
        row(
            "pjeg_20_proof_boundary",
            "proof_guard",
            "guard_active",
            (
                "The joined identity, diagonal bounds, subpower "
                "equivalence, and route collapse are exact."
            ),
            (
                "No L, E_B, or Y estimate, RH, PF-infinity, Lambda<=0, "
                "or Clay-prize result is proved."
            ),
        ),
    ]


DIAGNOSTICS = [
    {
        "K": 16,
        "R": 8,
        "L": 0.086714164,
        "D_L": 1.276405173,
        "O": -0.594845505,
        "E_B": 0.110319464,
        "O_CC_int": -0.001713988,
        "J_B": 0.053445744,
        "J_minus_L_over_2": 0.010088662,
    },
    {
        "K": 32,
        "R": 14,
        "L": 0.074555912,
        "D_L": 1.268291221,
        "O": -0.596867655,
        "E_B": 0.083247421,
        "O_CC_int": -0.000825629,
        "J_B": 0.040798082,
        "J_minus_L_over_2": 0.003520126,
    },
    {
        "K": 64,
        "R": 23,
        "L": 0.129088701,
        "D_L": 1.278746522,
        "O": -0.574828911,
        "E_B": 0.135323930,
        "O_CC_int": 0.000052069,
        "J_B": 0.067714034,
        "J_minus_L_over_2": 0.003169684,
    },
    {
        "K": 128,
        "R": 39,
        "L": 0.143665376,
        "D_L": 1.315481547,
        "O": -0.585908085,
        "E_B": 0.147238007,
        "O_CC_int": 0.000053254,
        "J_B": 0.073672258,
        "J_minus_L_over_2": 0.001839570,
    },
    {
        "K": 256,
        "R": 64,
        "L": 0.112963133,
        "D_L": 1.324398989,
        "O": -0.605717928,
        "E_B": 0.115185299,
        "O_CC_int": -0.000166902,
        "J_B": 0.057425747,
        "J_minus_L_over_2": 0.000944181,
    },
    {
        "K": 512,
        "R": 108,
        "L": 0.109596358,
        "D_L": 1.323092731,
        "O": -0.606748187,
        "E_B": 0.111068291,
        "O_CC_int": -0.000042452,
        "J_B": 0.055491694,
        "J_minus_L_over_2": 0.000693515,
    },
    {
        "K": 1024,
        "R": 182,
        "L": 0.096999647,
        "D_L": 1.325065222,
        "O": -0.614032788,
        "E_B": 0.097890613,
        "O_CC_int": 0.000020860,
        "J_B": 0.048966167,
        "J_minus_L_over_2": 0.000466344,
    },
]


NOTE_HEAD = r"""# Jensen-Window PF Mertens Planar Joined-Energy Equivalence Gate

Date: 2026-07-24

Status: exact equivalence between the proposed joined
current-interior/boundary coordinate and the inherited lossless low-mode
energy, with a summable coordinate defect. This is a route guard, not an
estimate and not a proof of RH, PF-infinity, `Lambda<=0`, or a Clay-prize
result.

```text
work/rh_compute/results/
jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.json
python work/rh_compute/scripts/
jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.py
python work/rh_compute/scripts/
check_jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.py
```

## The Full Off-Diagonal Already Knows The Positive Energy

Retain the notation of Corollaries 11.22Z.8-Z.10. Define

```text
L
 :=sum_(r=1)^R (x_r+beta e_r)^2/q_r^2
  =E_perp+c(beta+gamma/c)^2,

D_L
 :=sum_(r=1)^R D_(C,r)/q_r^2+cQ.                    (PJE.1)
```

The three-block compression from Corollary 11.22Z.8 is exactly

```text
O_(alpha,K)=(L-D_L)/2.                              (PJE.2)
```

The diagonal estimate from Corollary 11.22Z.3 gives

```text
0<=D_L<7pi^2/24,

|O_(alpha,K)-L/2|=D_L/2<7pi^2/48.                  (PJE.3)
```

Thus the complete signed off-diagonal is already the positive low-mode
energy up to a uniformly bounded diagonal.

## The Proposed Join Returns The Same Coordinate

Corollary 11.22Z.9 gives

```text
O_(alpha,K)=O_(CC,int)+R_B,

R_B=(E_B-D_B)/2+C_delta.                            (PJE.4)
```

Define the candidate joined quantity

```text
J_B:=O_(CC,int)+E_B/2.                              (PJE.5)
```

Eliminating `O_(alpha,K)` between (PJE.2) and (PJE.4) gives the exact
identity

```text
J_B
 =L/2+(D_B-D_L)/2-C_delta.                          (PJE.6)
```

The endpoint estimates already proved in Corollary 11.22Z.9 give

```text
|C_delta|<=pi H_R^odd+pi^2/2,
0<=D_B<5pi^2/12.
```

Together with (PJE.3),

```text
|J_B-L/2|
 <=pi H_R^odd+41pi^2/48,

|J_B-O_(alpha,K)|
 <=pi H_R^odd+17pi^2/24.                            (PJE.7)
```

For polynomially bounded `R`, every right side is `O(1+log K)`.
Consequently,

```text
J_B=O_epsilon(K^epsilon) for every epsilon>0
 iff
O_(alpha,K)=O_epsilon(K^epsilon) for every epsilon>0
 iff
L=O_epsilon(K^epsilon) for every epsilon>0.          (PJE.8)
```

For each fixed `alpha>0`, the stronger weighted statement is also exact:

```text
sum_(K dyadic)K^(-alpha)|J_B|<infinity
 iff
sum_(K dyadic)K^(-alpha)|O_(alpha,K)|<infinity
 iff
sum_(K dyadic)K^(-alpha)L<infinity.                  (PJE.9)
```

On one fixed cofinal positive-alpha sequence, the final member is the
prior lossless RH criterion.

## Route Decision

Keeping `O_(CC,int)` signed with `E_B/2` does not create a new weaker
interface. It returns `L/2` up to the same summable endpoint and diagonal
coordinate changes already isolated by the planar algebra. The options are
therefore honest:

```text
prove the RH-strength E_B anchored theorem;
prove the lossless L theorem;
prove the strong signed nested-edge theorem; or
find new structure outside this equivalent planar join.
```

The `Y_(K,r)` projection remains a useful exact coordinate, but pairing its
current-interior value with `E_B/2` cannot lower the theorem strength.

## Finite Mobius Diagnostic

At `alpha=1/2`, direct recomputation gives:

| K | R | L | D_L | O | E_B | O_CC,int | J_B | J_B-L/2 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
"""


NOTE_TAIL = r"""

The final column decreases on this finite grid because the two coordinate
systems are already close there. This is diagnostic only. The exact
identity (PJE.6), not the table, proves the route equivalence.

No estimate for `L`, `E_B`, or `Y_(K,r)`, no RH, no PF-infinity, no
`Lambda<=0`, and no Clay-prize result is proved.
"""


def note_text(rows: list[dict[str, str]]) -> str:
    lines = [NOTE_HEAD.rstrip()]
    for item in DIAGNOSTICS:
        lines.append(
            "| {K} | {R} | {L:.9f} | {D_L:.9f} | {O:.9f} | "
            "{E_B:.9f} | {O_CC_int:.9f} | {J_B:.9f} | "
            "{J_minus_L_over_2:.9f} |".format(**item)
        )
    lines.extend([NOTE_TAIL.rstrip(), "", "## Claim Ledger", ""])
    lines.append("| ID | Role | Status | Statement | Boundary |")
    lines.append("|---|---|---|---|---|")
    for item in rows:
        lines.append(
            "| `{id}` | `{role}` | `{status}` | {statement} | "
            "{proof_boundary} |".format(**item)
        )
    lines.extend(
        [
            "",
            "Summary:",
            "",
            f"- rows: {len(rows)}",
            "- exact/theorem-composed results: 16",
            "- finite validations: 1",
            "- active route/proof guards: 3",
            "- open lossless energy gates: 1",
            "",
        ]
    )
    return "\n".join(lines)


def payload(rows: list[dict[str, str]]) -> dict:
    return {
        "kind": (
            "jensen_window_pf_mertens_planar_"
            "joined_energy_equivalence_gate"
        ),
        "date": "2026-07-24",
        "status": (
            "exact joined current-interior/boundary equivalence "
            "to the lossless low-mode energy"
        ),
        "source": {
            "formal_core": "outputs/formal_core.md",
            "boundary_flux": (
                "outputs/"
                "jensen_window_pf_mertens_planar_"
                "boundary_flux_anchor_handoff.md"
            ),
            "suffix_energy": (
                "outputs/"
                "jensen_window_pf_mertens_planar_"
                "boundary_suffix_energy_handoff.md"
            ),
            "suffix_localization": (
                "outputs/"
                "jensen_window_pf_mertens_planar_"
                "boundary_suffix_localization_gate.md"
            ),
        },
        "rows": rows,
        "diagnostic": {
            "alpha": 0.5,
            "rows": DIAGNOSTICS,
        },
        "summary": {
            "row_count": len(rows),
            "exact_result_count": 16,
            "finite_validation_count": 1,
            "guard_count": 3,
            "open_energy_gate_count": 1,
            "joined_identity_proved": True,
            "weighted_dyadic_equivalence_proved": True,
            "joined_route_distinct": False,
            "lossless_energy_estimate_proved": False,
            "boundary_energy_estimate_proved": False,
            "interior_projection_estimate_proved": False,
            "rh_proved": False,
            "pf_infinity_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows = build_rows()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload(rows), indent=2) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(note_text(rows), encoding="utf-8")
    print(
        "built planar joined-energy equivalence gate: "
        f"{len(rows)} rows, {len(DIAGNOSTICS)} diagnostics, "
        "1 open lossless energy gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
