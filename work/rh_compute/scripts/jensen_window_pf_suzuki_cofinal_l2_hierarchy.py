#!/usr/bin/env python3
"""Build the Suzuki cofinal L2-hierarchy theorem audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_suzuki_cofinal_l2_hierarchy.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_suzuki_cofinal_l2_hierarchy.md"
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
            "sl2_01_shifted_quotient",
            "published_exact",
            "available_exact",
            "Q_omega(z)=xi(1/2+z-omega)/xi(1/2+z+omega), with "
            "Q_omega(0)=1 for 0<omega<1/2.",
            "Definition and the xi functional equation only.",
        ),
        row(
            "sl2_02_mellin_laplace_identity",
            "published_exact",
            "available_exact",
            "The Laplace transform of H_(omega,k)(exp(t)) is "
            "Q_omega(z)/z^k, initially for Re(z)>1/2+omega.",
            "Published Suzuki identity; no contour crosses a quotient pole.",
        ),
        row(
            "sl2_03_central_taylor_polynomial",
            "exact_definition",
            "available_exact",
            "Write Q_omega(z)=sum_(j>=0)q_(omega,j)z^j locally and "
            "T_(omega,k-1)(z)=sum_(j=0)^(k-1)q_(omega,j)z^j.",
            "A local Taylor definition at the zero-free real point z=0.",
        ),
        row(
            "sl2_04_log_polynomial",
            "exact_definition",
            "available_exact",
            "P_(omega,k)(t)=sum_(j=0)^(k-1)q_(omega,j)"
            "*t^(k-1-j)/(k-1-j)!.",
            "The inverse Laplace image of T_(omega,k-1)(z)/z^k.",
        ),
        row(
            "sl2_05_residual_transform",
            "exact_identity",
            "available_exact",
            "For r_(omega,k)(t)=H_(omega,k)(exp(t))-P_(omega,k)(t), "
            "Laplace(r)=[Q_omega-T_(omega,k-1)]/z^k in the initial "
            "right half-plane.",
            "Exact subtraction identity only.",
        ),
        row(
            "sl2_06_l2_to_hardy",
            "published_functional_analysis",
            "available_exact",
            "If r_(omega,k) belongs to L2(0,infinity), its Laplace "
            "transform belongs to H2 of the right half-plane.",
            "Standard L2 Paley-Wiener theorem.",
        ),
        row(
            "sl2_07_hardy_to_pole_free",
            "exact_implication",
            "available_exact",
            "The H2 continuation in sl2_06 gives "
            "Q_omega=T_(omega,k-1)+z^k Laplace(r), so Q_omega has no "
            "poles in Re(z)>0.",
            "Uses transform uniqueness from the common initial half-plane.",
        ),
        row(
            "sl2_08_pole_free_to_inner",
            "published_implication",
            "available_exact",
            "For the xi quotient, pole-freeness in Re(z)>0 plus boundary "
            "unimodularity, Suzuki's high-strip estimate, and "
            "Phragmen-Lindelof imply that Q_omega is inner.",
            "The xi-specific growth input is essential.",
        ),
        row(
            "sl2_09_inner_remainder_hardy",
            "exact_implication",
            "available_exact",
            "If Q_omega is inner, "
            "[Q_omega-T_(omega,k-1)]/z^k belongs to H2: Taylor "
            "cancellation controls z=0 and the boundary tail is O(1/|t|).",
            "Requires innerness, not boundary modulus one alone.",
        ),
        row(
            "sl2_10_hardy_to_l2_residual",
            "published_functional_analysis",
            "available_exact",
            "Paley-Wiener and Laplace uniqueness turn the H2 function in "
            "sl2_09 into the actual residual r_(omega,k) in L2.",
            "Standard half-plane Hardy/Laplace correspondence.",
        ),
        row(
            "sl2_11_fixed_shift_equivalence",
            "theorem_candidate",
            "available_exact",
            "For every fixed omega and k>=1, "
            "r_(omega,k) in L2(0,infinity) iff Q_omega is inner iff "
            "D(omega).",
            "Internally audited composition requiring independent review.",
        ),
        row(
            "sl2_12_cofinal_reverse",
            "theorem_candidate",
            "available_exact",
            "If omega_j tends to zero and one residual "
            "r_(omega_j,k_j) is L2 for every j, then RH.",
            "Uses the corpus fixed-omega phase theorem and cofinality.",
        ),
        row(
            "sl2_13_rh_forward",
            "published_conditional",
            "conditional_on_rh",
            "Under RH, Q_omega is inner and every residual "
            "r_(omega,k) is L2 for every fixed omega and k>=1.",
            "Conditional direction only; it cannot establish RH.",
        ),
        row(
            "sl2_14_cofinal_equivalence",
            "theorem_candidate",
            "available_exact",
            "RH iff there are omega_j->0 and integers k_j>=1 such that "
            "H_(omega_j,k_j)(exp(t))-P_(omega_j,k_j)(t) belongs to "
            "L2(0,infinity) for every j.",
            "New cofinal sharpening of Suzuki's published L2 criterion.",
        ),
        row(
            "sl2_15_jordan_error_form",
            "exact_handoff",
            "available_exact",
            "The same residual equals x^(-1/2) integral_0^1 "
            "E_omega(xu)W_(omega,k)(u)du/u-P_(omega,k)(log x).",
            "Uses the exact Jordan-error kernel reduction.",
        ),
        row(
            "sl2_16_level_one_target",
            "exact_specialization",
            "available_exact",
            "At k=1, P_(omega,1)=1, so the criterion is "
            "integral_1^infinity |H_(omega,1)(x)-1|^2 dx/x<infinity.",
            "Exact level-one specialization only.",
        ),
        row(
            "sl2_17_higher_level_target",
            "exact_specialization",
            "available_exact",
            "At k=2, q_(omega,1)=-2 xi'(1/2+omega)/xi(1/2+omega), "
            "so P_(omega,2)(t)=t+q_(omega,1) and "
            "r_(omega,2)(t)=-q_(omega,1)+integral_0^t "
            "(H_(omega,1)(exp(u))-1)du. Higher levels use the same "
            "central Taylor recursion.",
            "Exact residue-polynomial and antiderivative identities, not "
            "an unconditional decay estimate.",
        ),
        row(
            "sl2_18_boundary_trace_guard",
            "countermodel_gate",
            "guard_validated",
            "Q_bad(z)=(a+z)/(a-z) has modulus one on the imaginary axis "
            "and its regularized boundary trace is L2, but it has a pole "
            "at z=a. The boundary inverse Fourier transform is anti-causal, "
            "while the actual positive-time residual contains exp(a t) and "
            "is not L2.",
            "Boundary energy cannot replace half-plane analyticity or "
            "positive-time Hardy support.",
        ),
        row(
            "sl2_19_fixed_shift_guard",
            "countermodel_gate",
            "guard_validated",
            "An L2 residual at one fixed omega proves only D(omega); exact "
            "horizontal zero cancellation can still hide off-line zeros.",
            "Cofinality in omega remains indispensable.",
        ),
        row(
            "sl2_20_finite_energy_guard",
            "countermodel_gate",
            "guard_validated",
            "A finite numerical integral of |r_(omega,k)|^2 cannot certify "
            "integrability on the unbounded half-line.",
            "Finite energy samples remain nonpromotable evidence.",
        ),
        row(
            "sl2_21_riesz_comparison_guard",
            "literature_fit_guard",
            "guard_validated",
            "Classical smoothed-totient Riesz means retain a residue main "
            "term and obtain square-root bounds only under RH and extra "
            "zero-derivative hypotheses; they do not prove the Suzuki L2 "
            "premise.",
            "The similar smoothing vocabulary does not identify the two "
            "kernels or remove the open arithmetic gate.",
        ),
        row(
            "sl2_22_open_energy_target",
            "open_arithmetic_gate",
            "open_target",
            "For one explicit omega_j->0 and convenient k_j, prove the "
            "Jordan-error residual in sl2_15 is in L2(dx/x) by direct "
            "arithmetic estimates without assuming a zero-free region.",
            "This is an RH-strength obligation and remains open.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_suzuki_cofinal_l2_hierarchy",
        "date": "2026-07-23",
        "status": (
            "internally audited theorem-candidate L2 hierarchy with one "
            "open arithmetic energy gate"
        ),
        "proof_boundary": (
            "This hierarchy does not prove an L2 residual, RH, or "
            "Lambda<=0. It composes published Suzuki identities and "
            "Hardy-space machinery with the internally audited fixed-shift "
            "phase theorem; independent expert review is required."
        ),
        "definitions": {
            "shifted_quotient": (
                "Q_omega(z)=xi(1/2+z-omega)/xi(1/2+z+omega)"
            ),
            "central_polynomial": (
                "T_(omega,k-1)(z)=sum_(j=0)^(k-1)q_(omega,j)z^j"
            ),
            "log_polynomial": (
                "P_(omega,k)(t)=sum_(j=0)^(k-1)q_(omega,j)"
                "t^(k-1-j)/(k-1-j)!"
            ),
            "residual": (
                "r_(omega,k)(t)=H_(omega,k)(exp(t))-P_(omega,k)(t)"
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
                "use": "Mellin hierarchy, L2 criterion, and contour bounds",
            },
            {
                "title": (
                    "A canonical system of differential equations arising "
                    "from the Riemann zeta-function"
                ),
                "url": "https://arxiv.org/abs/1204.1827",
                "use": "level-one L2 criterion and innerness",
            },
            {
                "title": "Corpus fixed-omega phase diagram",
                "url": (
                    "outputs/"
                    "jensen_window_pf_suzuki_fixed_omega_phase_diagram.md"
                ),
                "use": "fixed-shift innerness and cofinal RH handoff",
            },
            {
                "title": "Corpus Jordan-error kernel reduction",
                "url": (
                    "outputs/"
                    "jensen_window_pf_suzuki_jordan_error_kernel_reduction.md"
                ),
                "use": "arithmetic error-convolution form",
            },
            {
                "title": (
                    "The Distribution of Error Terms of Smoothed Summatory "
                    "Totient Functions"
                ),
                "url": "https://arxiv.org/abs/2207.07722",
                "use": "Riesz-smoothing comparison and non-transfer guard",
            },
        ],
        "audit": {
            "row_count": len(rows),
            "countermodel_gate_count": sum(
                item["role"] == "countermodel_gate" for item in rows
            ),
            "literature_fit_guard_count": sum(
                item["role"] == "literature_fit_guard" for item in rows
            ),
            "open_arithmetic_gate_count": sum(
                item["role"] == "open_arithmetic_gate" for item in rows
            ),
            "fixed_shift_equivalence_candidate_count": sum(
                item["id"] == "sl2_11_fixed_shift_equivalence"
                for item in rows
            ),
            "l2_residual_proved": False,
            "rh_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note() -> str:
    return """# Jensen-Window PF Suzuki Cofinal L2 Hierarchy

Date: 2026-07-23

Status: internally audited theorem-candidate reduction with one open
arithmetic energy gate. This is not a proof of an `L2` residual, RH, or
`Lambda <= 0`. Independent expert review is required.

```text
work/rh_compute/results/jensen_window_pf_suzuki_cofinal_l2_hierarchy.json
python work/rh_compute/scripts/jensen_window_pf_suzuki_cofinal_l2_hierarchy.py
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_cofinal_l2_hierarchy.py
```

## Central Subtraction

For `0<omega<1/2`, put `z=s-1/2` and

```text
Q_omega(z)=xi(1/2+z-omega)/xi(1/2+z+omega).
```

The functional equation gives `Q_omega(0)=1`. Suzuki's smoothing hierarchy
has the Laplace form

```text
integral_0^infinity H_(omega,k)(exp(t))*exp(-z*t)dt
 =Q_omega(z)/z^k,                                    (SL2.1)
```

initially for `Re(z)>1/2+omega`.

Write the Taylor expansion at the central point as

```text
Q_omega(z)=sum_(j>=0)q_(omega,j)z^j,
T_(omega,k-1)(z)=sum_(j=0)^(k-1)q_(omega,j)z^j,

P_(omega,k)(t)
 =sum_(j=0)^(k-1)
  q_(omega,j)*t^(k-1-j)/(k-1-j)!.
```

Since the Laplace transform of `t^m/m!` is `z^(-m-1)`, the residual

```text
r_(omega,k)(t)
 =H_(omega,k)(exp(t))-P_(omega,k)(t)
```

satisfies the exact identity

```text
Laplace(r_(omega,k))(z)
 =[Q_omega(z)-T_(omega,k-1)(z)]/z^k.                (SL2.2)
```

At `k=1`, the subtraction is simply `P_(omega,1)=1`. At higher levels
the leading term is `t^(k-1)/(k-1)!`.

The first broadened target is completely explicit. The functional equation
gives

```text
q_(omega,1)
 =-2*xi'(1/2+omega)/xi(1/2+omega),

P_(omega,2)(t)=t+q_(omega,1),

r_(omega,2)(t)
 =-q_(omega,1)
  +integral_0^t (H_(omega,1)(exp(u))-1)du.            (SL2.2a)
```

Thus `k=2` asks for square integrability of a normalized cumulative
level-one discrepancy, not merely for ordinary smoothing.

## Fixed-Shift Equivalence

Suppose first that `r_(omega,k)` belongs to `L2(0,infinity)`. The
half-plane Paley-Wiener theorem makes its Laplace transform an `H2`
function in `Re(z)>0`. Equation (SL2.2), initially valid farther right,
then gives

```text
Q_omega(z)
 =T_(omega,k-1)(z)+z^k*Laplace(r_(omega,k))(z).
```

Thus the reduced quotient has no pole in the right half-plane. Boundary
unimodularity, Suzuki's high-strip estimate, and Phragmen-Lindelof promote
it to an inner function.

Conversely, suppose `Q_omega` is inner. Then

```text
F_(omega,k)(z)
 =[Q_omega(z)-T_(omega,k-1)(z)]/z^k
```

is bounded near `z=0` by Taylor cancellation. On the imaginary boundary
it is `O(1/|t|)` at infinity because `Q_omega` is bounded and the Taylor
polynomial has degree `k-1`. Hence `F_(omega,k)` belongs to `H2`.
Paley-Wiener and transform uniqueness recover the actual residual in
`L2(0,infinity)`. Therefore, for every fixed `omega` and every `k>=1`,

```text
r_(omega,k) in L2(0,infinity)
 iff Q_omega is inner
 iff D(omega).                                       (SL2.3)
```

This extends the published level-one `L2` criterion to the complete
logarithmic smoothing hierarchy.

## Cofinal Criterion

The fixed-omega phase theorem says that a sequence of shifts in `D`
accumulating at zero forces RH. Under RH every shifted quotient is inner.
Combining this with (SL2.3) gives the sharpened equivalence

```text
RH
iff there exist omega_j->0 and integers k_j>=1 such that
    H_(omega_j,k_j)(exp(t))-P_(omega_j,k_j)(t)
    belongs to L2(0,infinity) for every j.            (SL2.4)
```

The smoothing order may vary with the shift. This creates an energy target
parallel to, and weaker in shape than, eventual one-sign behavior.

## Arithmetic Form

The exact Jordan-error reduction turns the residual into

```text
x^(-1/2)*integral_0^1
 E_omega(x*u)W_(omega,k)(u)du/u
 -P_(omega,k)(log x).                                (SL2.5)
```

Thus the surviving arithmetic task is to prove, for one explicit
`omega_j->0` and convenient `k_j`, that (SL2.5) has finite squared norm
against `dx/x`, without assuming a zero-free half-plane. At level one this
is

```text
integral_1^infinity |H_(omega,1)(x)-1|^2 dx/x
 <infinity.                                          (SL2.6)
```

At level two it is the concrete condition

```text
integral_1^infinity |
 x^(-1/2)*integral_0^1 E_omega(x*u)W_(omega,2)(u)du/u
 -log(x)+2*xi'(1/2+omega)/xi(1/2+omega)
 |^2 dx/x < infinity.                                (SL2.7)
```

## Countermodel Guard

For `a>0`, the right-half-plane inner function

```text
Q_good(z)=(a-z)/(a+z)
```

has Taylor coefficients `q_0=1` and
`q_j=2*(-1)^j/a^j` for `j>=1`. Exact subtraction gives

```text
[Q_good(z)-T_(k-1)(z)]/z^k
 =2*(-1)^k/[a^(k-1)*(a+z)],
```

whose inverse Laplace transform is a decaying exponential in `L2`.

In contrast,

```text
Q_bad(z)=(a+z)/(a-z)
```

also has modulus one on the imaginary axis, but has a pole at `z=a`.
Its actual positive-time residual contains `exp(a*t)` and is not in `L2`.

There is a sharper warning. The regularized boundary trace

```text
B_(omega,k)(y)
 =[Q_omega(i*y)-T_(omega,k-1)(i*y)]/(i*y)^k
```

is already in `L2(R)` without RH: Taylor cancellation controls `y=0` and
boundary unimodularity gives an `O(1/|y|)` tail. For `Q_bad`, its inverse
Fourier transform is the anti-causal decaying function

```text
2*a^(-(k-1))*exp(a*t)*1_(t<0).
```

For `Q_good`, the inverse Fourier transform is supported on `t>=0`.
Paley-Wiener therefore makes positive-time support, not finite boundary
energy, the decisive condition. A Plancherel calculation of the boundary
norm alone is automatic and cannot prove (SL2.3).

One fixed shift is also insufficient because exact horizontal zero
cancellation can hide off-line zeros. A finite numerical energy integral
cannot certify the infinite tail.

## Riesz-Smoothing Comparison

Das, Lang, Wan, and Xu study positive logarithmic Riesz means of the
classical Euler-totient sum. Their square-root bound after two smoothing
operations assumes RH and a reciprocal-zeta-derivative moment bound. That
result is useful contour technology, but it neither identifies Suzuki's
signed archimedean kernel nor proves (SL2.5). In Suzuki's transform the
Jordan residue main term is annihilated exactly and the central Taylor
polynomial supplies the target profile.

## Sources

- Masatoshi Suzuki, `On monotonicity of certain weighted summatory functions associated with L-functions`: https://arxiv.org/abs/1204.1823
- Masatoshi Suzuki, `A canonical system of differential equations arising from the Riemann zeta-function`: https://arxiv.org/abs/1204.1827
- Sanjana Das, Hannah Lang, Hamilton Wan, and Nancy Xu, `The Distribution of Error Terms of Smoothed Summatory Totient Functions`: https://arxiv.org/abs/2207.07722
- `outputs/jensen_window_pf_suzuki_fixed_omega_phase_diagram.md`
- `outputs/jensen_window_pf_suzuki_jordan_error_kernel_reduction.md`

## Proof Boundary

Equations (SL2.1)-(SL2.3) use published Suzuki identities, standard
Paley-Wiener theory, and the internally audited fixed-shift innerness
handoff. Equation (SL2.4) is a new corpus theorem candidate requiring
independent review. No residual in (SL2.5) or (SL2.6) is proved
unconditionally, so this artifact does not prove RH or `Lambda <= 0`.
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
    args.note.write_text(render_note(), encoding="utf-8")
    print(
        "wrote Suzuki cofinal L2 hierarchy: "
        f"{len(payload['rows'])} rows -> {args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
