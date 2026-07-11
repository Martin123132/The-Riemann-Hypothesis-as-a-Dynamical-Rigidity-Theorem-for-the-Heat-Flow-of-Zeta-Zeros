#!/usr/bin/env python3
"""Build the Jensen-window PF heat-flow boundary-threshold lemma artifact."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402

from jensen_window_pf_heat_flow_monotone_closure_scout import (  # noqa: E402
    DEFAULT_ENCLOSURE_JSONL,
    REPO_ROOT,
    contraction,
    decimal_format,
    load_enclosures,
)


getcontext().prec = 90

DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_heat_flow_boundary_threshold_lemma.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md"


@dataclass(frozen=True)
class FiniteMinimum:
    sample: str
    lam: str
    k: int


@dataclass(frozen=True)
class ThresholdDiagnostics:
    lambdas: list[str]
    max_coefficient_index: int
    strong_threshold: str
    heat_flow_threshold: str
    strong_threshold_rows: int
    strong_threshold_positive_rows: int
    heat_flow_threshold_rows: int
    heat_flow_threshold_positive_rows: int
    min_strong_threshold_margin: FiniteMinimum
    min_heat_flow_threshold_margin: FiniteMinimum


def arb_positive(value: flint.arb) -> bool:
    return bool(value > 0 and not value.contains(0))


def finite_diagnostics(paths: list[Path], max_coefficient_index: int) -> ThresholdDiagnostics:
    balls, samples, labels = load_enclosures(paths)
    strong_rows = 0
    strong_positive = 0
    heat_rows = 0
    heat_positive = 0
    min_strong: tuple[Decimal, str, int] | None = None
    min_heat: tuple[Decimal, str, int] | None = None

    for lam in sorted(labels):
        for k in range(1, max_coefficient_index):
            arb_values = {idx: balls[(lam, idx)] for idx in (k - 1, k, k + 1)}
            sample_values = {idx: samples[(lam, idx)] for idx in (k - 1, k, k + 1)}
            arb_x = contraction(arb_values, k)
            sample_x = contraction(sample_values, k)

            strong_rows += 1
            strong = Decimal(2 * k - 1) / Decimal(2 * k + 1)
            arb_strong = flint.arb(2 * k - 1) / flint.arb(2 * k + 1)
            strong_margin = sample_x - strong
            if arb_positive(arb_x - arb_strong):
                strong_positive += 1
            if min_strong is None or strong_margin < min_strong[0]:
                min_strong = (strong_margin, labels[lam], k)

            heat_rows += 1
            heat = Decimal(2 * k - 1) / Decimal(2 * k + 5)
            arb_heat = flint.arb(2 * k - 1) / flint.arb(2 * k + 5)
            heat_margin = sample_x - heat
            if arb_positive(arb_x - arb_heat):
                heat_positive += 1
            if min_heat is None or heat_margin < min_heat[0]:
                min_heat = (heat_margin, labels[lam], k)

    if min_strong is None or min_heat is None:
        raise RuntimeError("no finite diagnostics were computed")
    return ThresholdDiagnostics(
        lambdas=[labels[lam] for lam in sorted(labels)],
        max_coefficient_index=max_coefficient_index,
        strong_threshold="x_k >= (2*k-1)/(2*k+1)",
        heat_flow_threshold="x_k >= (2*k-1)/(2*k+5)",
        strong_threshold_rows=strong_rows,
        strong_threshold_positive_rows=strong_positive,
        heat_flow_threshold_rows=heat_rows,
        heat_flow_threshold_positive_rows=heat_positive,
        min_strong_threshold_margin=FiniteMinimum(decimal_format(min_strong[0]), min_strong[1], min_strong[2]),
        min_heat_flow_threshold_margin=FiniteMinimum(decimal_format(min_heat[0]), min_heat[1], min_heat[2]),
    )


def build_lemma(paths: list[Path], max_coefficient_index: int) -> dict:
    diagnostics = finite_diagnostics(paths, max_coefficient_index)
    return {
        "kind": "jensen_window_pf_heat_flow_boundary_threshold_lemma",
        "date": "2026-07-06",
        "status": "available_exact_boundary_threshold_lemma",
        "source_heat_flow_closure": "outputs/jensen_window_pf_heat_flow_monotone_closure_scout.md",
        "source_theorem_target": "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "proof_boundary": (
            "Exact boundary-threshold lemma only. It proves the Cauchy-Schwarz "
            "lower bound x_k >= (2*k-1)/(2*k+1), hence the weaker heat-flow "
            "threshold x_k >= (2*k-1)/(2*k+5), and the corresponding boundary "
            "bracket nonnegativity under the monotone-neighbor assumption. It "
            "does not prove adjacent log-concavity x_k<=1, does not prove "
            "increasing contractions for all zeta windows, does not prove a "
            "closed differential inequality, does not prove jwpf_06, and does "
            "not prove RH or Lambda <= 0."
        ),
        "exact_rows": [
            {
                "id": "hfbtl_01_phi_positive_half_line",
                "role": "exact_kernel_positivity",
                "formula": "Phi(u)=sum_{n>=1} pi*n^2*exp(5*u)*(2*pi*n^2*exp(4*u)-3)*exp(-pi*n^2*exp(4*u)) > 0 for u>=0",
                "argument": "For n>=1 and u>=0, 2*pi*n^2*exp(4*u)-3 >= 2*pi-3 > 0, and every other factor is positive.",
                "proof_boundary": "Exact half-line positivity for the Newman/Phi moment weight only.",
            },
            {
                "id": "hfbtl_02_raw_moment_logconvexity",
                "role": "exact_moment_inequality",
                "formula": "mu_{2k}(lambda)^2 <= mu_{2k-2}(lambda)*mu_{2k+2}(lambda)",
                "argument": "Apply Cauchy-Schwarz to u^(2*k-1)*sqrt(w_lambda(u)) and u^(2*k+1)*sqrt(w_lambda(u)), where w_lambda(u)=2*exp(lambda*u^2)*Phi(u) is positive on u>=0.",
                "proof_boundary": "Exact raw-moment log-convexity only; it is not normalized monotone contraction.",
            },
            {
                "id": "hfbtl_03_normalized_strong_threshold",
                "role": "exact_normalization_reduction",
                "formula": "x_k = ((2*k-1)/(2*k+1))*(mu_{2k+2}*mu_{2k-2}/mu_{2k}^2) >= (2*k-1)/(2*k+1)",
                "argument": "Substitute A_k=mu_{2k}*k!/(2*k)! into x_k=(A_{k+1}/A_k)/(A_k/A_{k-1}) and use raw-moment log-convexity.",
                "proof_boundary": "Exact lower bound for x_k only; it does not prove x_k<=1 or x_{k+1}>=x_k.",
            },
            {
                "id": "hfbtl_04_heat_flow_threshold",
                "role": "exact_heat_flow_threshold",
                "formula": "x_k >= (2*k-1)/(2*k+1) > (2*k-1)/(2*k+5)",
                "argument": "The strong normalized threshold strictly dominates the boundary threshold needed by the heat-flow bracket.",
                "proof_boundary": "Exact threshold discharge only; it does not close the heat-flow differential inequality.",
            },
            {
                "id": "hfbtl_05_boundary_bracket_nonnegative",
                "role": "exact_boundary_pointing_condition",
                "formula": "if x_k=x_{k+1}=q, x_{k+2}>=q, and q>=((2*k-1)/(2*k+5)), then bracket >= ((q-1)^2*((2*k+5)*q-(2*k-1)))/q >= 0",
                "argument": "Use the heat-flow bracket factorization from the monotone-closure scout and the exact threshold bound.",
                "proof_boundary": "Boundary pointing condition only; not a global flow-invariance theorem and not an initial-condition theorem.",
            },
        ],
        "finite_diagnostics": asdict(diagnostics),
        "summary": {
            "exact_rows": 5,
            "strong_threshold_rows": diagnostics.strong_threshold_rows,
            "strong_threshold_positive_rows": diagnostics.strong_threshold_positive_rows,
            "heat_flow_threshold_rows": diagnostics.heat_flow_threshold_rows,
            "heat_flow_threshold_positive_rows": diagnostics.heat_flow_threshold_positive_rows,
            "boundary_threshold_closed": True,
            "monotone_contraction_target_closing": False,
            "main_finding": (
                "The heat-flow boundary threshold is not an open subtarget: "
                "Phi positivity and Cauchy-Schwarz raw-moment log-convexity "
                "prove the stronger bound x_k >= (2*k-1)/(2*k+1), which "
                "implies the needed x_k >= (2*k-1)/(2*k+5). The remaining "
                "heat-flow route still needs adjacent log-concavity, a global "
                "flow-invariance argument, and initial or terminal conditions."
            ),
        },
        "invariants": [
            "The lemma keeps A_k=mu_{2k}*k!/(2*k)! explicit.",
            "The boundary threshold is closed, but the monotone-contraction theorem remains open.",
            "Finite threshold margins are sanity evidence, not the proof of the exact lemma.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(lemma: dict, path: Path) -> None:
    diagnostics = lemma["finite_diagnostics"]
    summary = lemma["summary"]
    text = f"""# Jensen-Window PF Heat-Flow Boundary-Threshold Lemma

Date: 2026-07-06

Status: exact boundary-threshold lemma. This is not a proof of the
monotone-contraction theorem, Jensen-window PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_heat_flow_boundary_threshold_lemma`.

Proof boundary: this artifact proves only the lower-bound threshold needed by
the heat-flow boundary bracket. It does not prove adjacent log-concavity
`x_k<=1`, does not prove `x_{{k+1}}>=x_k` for all zeta windows, does not prove
a closed differential inequality, does not prove `jwpf_06`, and does not prove
`Lambda <= 0`.

Machine-readable lemma:

```text
work/rh_compute/results/jensen_window_pf_heat_flow_boundary_threshold_lemma.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_heat_flow_boundary_threshold_lemma.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_boundary_threshold_lemma.py
```

Current result:

```text
validated Jensen-window PF heat-flow boundary threshold lemma: 5 exact rows, 315 strong-threshold rows, 315 heat-threshold rows, 0 issues
```

## Exact Lemma

For `u>=0`,

```text
Phi(u)=sum_{{n>=1}} pi*n^2*exp(5*u)*(2*pi*n^2*exp(4*u)-3)*exp(-pi*n^2*exp(4*u)) > 0
```

because `2*pi*n^2*exp(4*u)-3 >= 2*pi-3 > 0`.

Set:

```text
w_lambda(u)=2*exp(lambda*u^2)*Phi(u)
mu_{{2k}}(lambda)=integral_0^infty u^(2*k) w_lambda(u) du
A_k(lambda)=mu_{{2k}}(lambda)*k!/(2*k)!
x_k=(A_{{k+1}}/A_k)/(A_k/A_{{k-1}})
```

Cauchy-Schwarz for the positive weight gives raw moment log-convexity:

```text
mu_{{2k}}(lambda)^2 <= mu_{{2k-2}}(lambda)*mu_{{2k+2}}(lambda)
```

The factorial normalization gives:

```text
x_k = ((2*k-1)/(2*k+1))*(mu_{{2k+2}}*mu_{{2k-2}}/mu_{{2k}}^2)
```

Therefore:

```text
x_k >= (2*k-1)/(2*k+1) > (2*k-1)/(2*k+5)
```

So the heat-flow boundary threshold from
`outputs/jensen_window_pf_heat_flow_monotone_closure_scout.md` is discharged.

## Boundary Bracket

At a monotone-gap boundary `x_{{k+1}}=x_k=q`, with the next monotone-neighbor
condition `x_{{k+2}}>=q`, the heat-flow bracket has lower bound:

```text
((q-1)^2*((2*k+5)*q-(2*k-1)))/q
```

The exact threshold makes this lower bound nonnegative. This is only a
boundary pointing condition; it is not yet a global flow-invariance theorem.

## Finite Sanity Check

The existing coefficient enclosure cache gives:

```text
lambdas: {", ".join(diagnostics["lambdas"])}
max coefficient index: {diagnostics["max_coefficient_index"]}
strong-threshold rows: {diagnostics["strong_threshold_rows"]}
strong-threshold positive rows: {diagnostics["strong_threshold_positive_rows"]}
heat-threshold rows: {diagnostics["heat_flow_threshold_rows"]}
heat-threshold positive rows: {diagnostics["heat_flow_threshold_positive_rows"]}
```

Weakest strong-threshold margin:

```text
{diagnostics["min_strong_threshold_margin"]["sample"]} at lambda={diagnostics["min_strong_threshold_margin"]["lam"]}, k={diagnostics["min_strong_threshold_margin"]["k"]}
```

Weakest heat-threshold margin:

```text
{diagnostics["min_heat_flow_threshold_margin"]["sample"]} at lambda={diagnostics["min_heat_flow_threshold_margin"]["lam"]}, k={diagnostics["min_heat_flow_threshold_margin"]["k"]}
```

## Remaining Gap

```text
{summary["main_finding"]}
```

Integration:

```text
outputs/jensen_window_pf_heat_flow_monotone_closure_scout.md
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
```
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--enclosure-jsonl", type=Path, action="append")
    parser.add_argument("--max-coefficient-index", type=int, default=64)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    paths = args.enclosure_jsonl or list(DEFAULT_ENCLOSURE_JSONL)
    lemma = build_lemma(paths, args.max_coefficient_index)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(lemma, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(lemma, args.note)
    diagnostics = lemma["finite_diagnostics"]
    print(
        "wrote Jensen-window PF heat-flow boundary threshold lemma: "
        f"{lemma['summary']['exact_rows']} exact rows, "
        f"{diagnostics['strong_threshold_rows']} strong-threshold rows, "
        f"{diagnostics['heat_flow_threshold_rows']} heat-threshold rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
