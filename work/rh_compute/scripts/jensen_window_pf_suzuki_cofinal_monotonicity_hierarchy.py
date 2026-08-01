#!/usr/bin/env python3
"""Build the Suzuki cofinal monotonicity-hierarchy theorem audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.md"
)


def row(
    row_id: str,
    role: str,
    status: str,
    statement: str,
    boundary: str,
) -> dict[str, str]:
    return {
        "id": row_id,
        "role": role,
        "status": status,
        "statement": statement,
        "proof_boundary": boundary,
    }


def build_payload() -> dict:
    rows = [
        row(
            "cmh_01_coefficient_identity",
            "published_exact",
            "available_exact",
            "c_omega(n)=n^omega*sum_(d|n)mu(d)d^(-2omega)"
            "=n^omega*product_(p|n)(1-p^(-2omega))>0.",
            "Exact arithmetic identity only.",
        ),
        row(
            "cmh_02_first_weight",
            "published_exact",
            "available_exact",
            "Suzuki's explicit beta-integral weight g_(omega,1) is "
            "supported on (0,1] for every 0<omega<1/2.",
            "Definition of the first smoothing weight only.",
        ),
        row(
            "cmh_03_weight_recursion",
            "published_exact",
            "available_exact",
            "g_(omega,k+1)(x)=integral_x^1 sqrt(y/x)"
            "*g_(omega,k)(y)dy/y.",
            "Exact logarithmic smoothing recursion only.",
        ),
        row(
            "cmh_04_summatory_definition",
            "published_exact",
            "available_exact",
            "H_(omega,k)(x)=x^(-1/2)*sum_(n<=x)c_omega(n)"
            "*g_(omega,k)(n/x), with H=0 on 0<x<1.",
            "Exact finite arithmetic sum only.",
        ),
        row(
            "cmh_05_log_antiderivative",
            "published_exact",
            "available_exact",
            "H_(omega,k+1)(x)=integral_1^x H_(omega,k)(y)dy/y.",
            "Exact smoothing identity; it supplies no eventual sign.",
        ),
        row(
            "cmh_06_mellin_identity",
            "published_exact",
            "available_exact",
            "integral_1^infinity H_(omega,k)(x)x^(1/2-s)dx/x"
            "=[xi(s-omega)/xi(s+omega)]/(s-1/2)^k, initially "
            "for Re(s)>1+omega.",
            "Published Mellin identity only; contour continuation may not "
            "cross uncancelled denominator zeros.",
        ),
        row(
            "cmh_07_k1_normalization",
            "exact_identification",
            "available_exact",
            "H_(omega,1)(x)=sqrt(x)*h_omega^<1>(x) in the canonical-"
            "system notation.",
            "Exact normalization bridge only.",
        ),
        row(
            "cmh_08_eventual_sign_landau",
            "published_implication",
            "available_exact",
            "If H_(omega,k) has one sign for all sufficiently large x, "
            "Landau's theorem and the absence of real xi zeros remove every "
            "uncancelled pole of xi(s-omega)/xi(s+omega) in Re(s)>1/2.",
            "The eventual-sign antecedent is not proved.",
        ),
        row(
            "cmh_09_pole_free_inner",
            "published_implication",
            "available_exact",
            "For the xi quotient, pole-freeness in Re(s)>1/2 plus Suzuki's "
            "high-strip estimate and Phragmen-Lindelof gives reduced "
            "meromorphic innerness.",
            "The xi-specific growth input is essential.",
        ),
        row(
            "cmh_10_fixed_shift_handoff",
            "exact_handoff",
            "available_exact",
            "Eventual one-sign behavior of one H_(omega,k) implies D(omega), "
            "the fixed-omega all-time determinant property.",
            "One fixed omega can still hide off-line zeros by exact "
            "horizontal cancellation.",
        ),
        row(
            "cmh_11_cofinal_handoff",
            "theorem_candidate",
            "available_exact",
            "If omega_j decreases to zero and every H_(omega_j,k_j) is "
            "eventually one-signed for some integers k_j>=1, then RH.",
            "Internally audited composition of published Landau machinery "
            "with the corpus fixed-omega phase theorem; external review is "
            "required.",
        ),
        row(
            "cmh_12_rh_k1_asymptotic",
            "published_conditional",
            "conditional_on_rh",
            "Under RH, H_(omega,1)(x)=1+O_omega(x^(-B_omega)) for some "
            "B_omega>0 and each fixed 0<omega<1/2.",
            "Conditional asymptotic; it cannot be used to prove its own "
            "antecedent.",
        ),
        row(
            "cmh_13_rh_higher_asymptotic",
            "published_conditional",
            "conditional_on_rh",
            "Under RH, H_(omega,k)(x) has positive leading term "
            "(log x)^(k-1)/(k-1)! for every fixed k>=2.",
            "Conditional contour shift only.",
        ),
        row(
            "cmh_14_cofinal_equivalence",
            "theorem_candidate",
            "available_exact",
            "RH is equivalent to the existence of omega_j decreasing to "
            "zero and integers k_j>=1 for which every H_(omega_j,k_j) has "
            "an eventual single sign.",
            "The forward implication uses Suzuki's RH-conditional "
            "asymptotics; the reverse implication is a corpus sharpening "
            "requiring external review.",
        ),
        row(
            "cmh_15_near_zero_first_weight",
            "exact_asymptotic",
            "available_exact",
            "g_(omega,1)(x)~-a_omega*x^(omega-1) as x->0+, where "
            "a_omega>0 for 0<omega<1/2.",
            "Exact endpoint sign, not a summatory sign.",
        ),
        row(
            "cmh_16_near_zero_all_weights",
            "exact_asymptotic",
            "available_exact",
            "g_(omega,k)(x)~-a_omega*(1/2-omega)^(-(k-1))"
            "*x^(omega-1) as x->0+.",
            "Derived recursively; it proves every smoothing weight remains "
            "negative near zero.",
        ),
        row(
            "cmh_17_near_one_all_weights",
            "exact_asymptotic",
            "available_exact",
            "For every k>=1, g_(omega,k)(x) is positive near x=1 and "
            "vanishes there like a positive constant times "
            "(1-x)^(omega+k-1).",
            "Exact endpoint sign only.",
        ),
        row(
            "cmh_18_signed_weight_guard",
            "countermodel_gate",
            "guard_validated",
            "Every g_(omega,k) is signed: negative near zero and positive "
            "near one, despite c_omega(n)>0.",
            "Coefficient positivity cannot prove any smoothing level "
            "termwise.",
        ),
        row(
            "cmh_19_finite_range_guard",
            "countermodel_gate",
            "guard_validated",
            "Finite positivity of H_(omega,k) cannot establish eventual sign "
            "on an unbounded half-line.",
            "No finite x grid closes the Landau antecedent.",
        ),
        row(
            "cmh_20_fixed_shift_guard",
            "countermodel_gate",
            "guard_validated",
            "Eventual sign at one fixed omega gives innerness only at that "
            "shift; exact horizontal zero cancellation prevents promotion "
            "to RH.",
            "Cofinality in omega is indispensable.",
        ),
        row(
            "cmh_21_open_smoothed_target",
            "open_arithmetic_gate",
            "open_target",
            "Prove eventual one-sign behavior for H_(omega_j,k_j) along one "
            "explicit sequence omega_j->0, with any chosen k_j>=1, by direct "
            "arithmetic inequalities that do not assume a zero-free region.",
            "This is the surviving RH-strength arithmetic obligation.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy",
        "date": "2026-07-23",
        "status": (
            "internally audited theorem-candidate reduction with one open "
            "arithmetic gate"
        ),
        "proof_boundary": (
            "This hierarchy is not a proof of eventual sign, RH, or "
            "Lambda<=0. The cofinal equivalence composes published Suzuki "
            "results with the internally audited fixed-omega phase theorem "
            "and requires independent expert review."
        ),
        "definitions": {
            "Q_omega": "xi(s-omega)/xi(s+omega)",
            "coefficient": (
                "c_omega(n)=n^omega*product_(p|n)(1-p^(-2omega))"
            ),
            "weight_recursion": (
                "g_(omega,k+1)(x)=integral_x^1 "
                "sqrt(y/x)g_(omega,k)(y)dy/y"
            ),
            "summatory_function": (
                "H_(omega,k)(x)=x^(-1/2)sum_(n<=x)"
                "c_omega(n)g_(omega,k)(n/x)"
            ),
        },
        "rows": rows,
        "sources": [
            {
                "title": (
                    "On monotonicity of certain weighted summatory "
                    "functions associated with L-functions"
                ),
                "url": "https://arxiv.org/abs/1204.1823",
                "use": (
                    "weight hierarchy, Mellin identity, Landau implication, "
                    "and RH-conditional eventual signs"
                ),
            },
            {
                "title": (
                    "A canonical system of differential equations arising "
                    "from the Riemann zeta-function"
                ),
                "url": "https://arxiv.org/abs/1204.1827",
                "use": "k=1 normalization, L2 criterion, and innerness",
            },
            {
                "title": "Corpus fixed-omega phase diagram",
                "url": (
                    "outputs/"
                    "jensen_window_pf_suzuki_fixed_omega_phase_diagram.md"
                ),
                "use": "cofinal innerness-to-RH sharpening",
            },
        ],
        "audit": {
            "row_count": len(rows),
            "published_or_exact_rows": sum(
                item["status"]
                in {"available_exact", "conditional_on_rh"}
                for item in rows
            ),
            "countermodel_gate_count": sum(
                item["role"] == "countermodel_gate" for item in rows
            ),
            "open_arithmetic_gate_count": sum(
                item["role"] == "open_arithmetic_gate" for item in rows
            ),
            "eventual_sign_proved": False,
            "rh_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note(payload: dict) -> str:
    return """# Jensen-Window PF Suzuki Cofinal Monotonicity Hierarchy

Date: 2026-07-23

Status: internally audited theorem-candidate reduction with one open
arithmetic gate. This is not a proof of eventual sign, RH, or
`Lambda <= 0`. Independent expert review is required.

```text
work/rh_compute/results/jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.json
python work/rh_compute/scripts/jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.py
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.py
```

## Published Hierarchy

For `0<omega<1/2`, Suzuki defines

```text
c_omega(n)
 =n^omega*sum_(d|n)mu(d)d^(-2omega)
 =n^omega*product_(p|n)(1-p^(-2omega))>0.
```

Starting from his explicit beta-integral weight `g_(omega,1)`, put

```text
g_(omega,k+1)(x)
 =integral_x^1 sqrt(y/x)*g_(omega,k)(y)dy/y,

H_(omega,k)(x)
 =x^(-1/2)*sum_(n<=x)c_omega(n)*g_(omega,k)(n/x).
```

The smoothing is a logarithmic antiderivative:

```text
H_(omega,k+1)(x)=integral_1^x H_(omega,k)(y)dy/y.
```

Suzuki's Mellin identity is

```text
integral_1^infinity H_(omega,k)(x)x^(1/2-s)dx/x
 =[xi(s-omega)/xi(s+omega)]/(s-1/2)^k,              (CMH.1)
```

initially for `Re(s)>1+omega`. At level `k=1`,

```text
H_(omega,1)(x)=sqrt(x)*h_omega^<1>(x)
```

in the canonical-system notation used by the fixed-omega phase diagram.

## Cofinal Equivalence

If one `H_(omega,k)` has a single sign for every sufficiently large `x`,
Landau's theorem applied to (CMH.1), together with the absence of real xi
zeros, removes every uncancelled pole of the reduced quotient in
`Re(s)>1/2`. Suzuki's high-strip estimate and Phragmen-Lindelof then give
meromorphic innerness. Therefore

```text
eventual one-sign of H_(omega,k) => D(omega).        (CMH.2)
```

The fixed-omega phase theorem makes cofinality decisive:

```text
If omega_j decreases to zero and H_(omega_j,k_j) is eventually
one-signed for every j, for any integers k_j>=1, then RH.       (CMH.3)
```

Conversely, Suzuki proves under RH that, for each fixed
`0<omega<1/2`,

```text
H_(omega,1)(x)=1+O_omega(x^(-B_omega))
```

for some `B_omega>0`, while every fixed `k>=2` has a positive leading
term `(log x)^(k-1)/(k-1)!`. Thus the sharpened equivalence is

```text
RH
iff there exist omega_j->0 and integers k_j>=1 such that every
    H_(omega_j,k_j) has one eventual sign.            (CMH.4)
```

Suzuki's published statement asks for all shifts in an interval. The
weakening to one cofinal sequence is the new corpus consequence and depends
on the internally audited fixed-omega phase theorem.

## Signed-Weight Guard

Higher smoothing does not make coefficient positivity a termwise proof.
For `0<omega<1/2`, the explicit first weight satisfies

```text
g_(omega,1)(x)~-a_omega*x^(omega-1),  a_omega>0,
```

as `x->0+`. The recursion gives, for every `k>=1`,

```text
g_(omega,k)(x)
 ~-a_omega*(1/2-omega)^(-(k-1))*x^(omega-1)
```

near zero. Near one, `g_(omega,k)` is positive and vanishes like a
positive constant times `(1-x)^(omega+k-1)`. Every smoothing weight is
therefore genuinely signed. The positive arithmetic coefficients cannot
close (CMH.4) term by term.

## Surviving Target

The broadened arithmetic obligation is

```text
Choose an explicit omega_j->0 and any convenient smoothing orders k_j>=1.
Prove eventual one-sign behavior of H_(omega_j,k_j) for every j by direct
arithmetic inequalities, without assuming a zero-free half-plane.
```

This is weaker than proving the unsmoothed `k=1` target at every shift.
The `k=2` level is the first natural focus because it is one logarithmic
antiderivative of the sampled function while retaining the same Landau
consequence.

Finite positivity remains nonpromotable, and one fixed shift remains
insufficient because exact horizontal zero cancellation can hide off-line
zeros.

## Sources

- Masatoshi Suzuki, `On monotonicity of certain weighted summatory functions associated with L-functions`: https://arxiv.org/abs/1204.1823
- Masatoshi Suzuki, `A canonical system of differential equations arising from the Riemann zeta-function`: https://arxiv.org/abs/1204.1827
- `outputs/jensen_window_pf_suzuki_fixed_omega_phase_diagram.md`

## Proof Boundary

The published hierarchy and Landau implication are exact source-backed
inputs. The cofinal weakening is an internally audited theorem candidate
obtained by composing them with the fixed-omega phase theorem. No eventual
sign is proved at any subcritical cofinal family, and this artifact does not
prove RH or Lambda<=0.
"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Suzuki cofinal monotonicity hierarchy: "
        f"{payload['audit']['row_count']} rows, "
        f"{payload['audit']['countermodel_gate_count']} guards, "
        f"{payload['audit']['open_arithmetic_gate_count']} open gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
