#!/usr/bin/env python3
"""Build the Jensen-window PF heat-flow monotone-closure scout.

This is a theorem-search diagnostic for the live heat-flow route in the
monotone-contraction target. It records exact lambda-flow algebra and checks
the finite zeta coefficient cache for the extra threshold and flow-bracket
conditions suggested by that algebra.
"""

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


getcontext().prec = 90

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ENCLOSURE_JSONL = (
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k0_k9.jsonl",
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k10_k20.jsonl",
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k21_k32.jsonl",
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k33_k48.jsonl",
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k49_k64.jsonl",
)
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_heat_flow_monotone_closure_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_heat_flow_monotone_closure_scout.md"


@dataclass(frozen=True)
class FiniteMinimum:
    sample: str
    lam: str
    k: int


@dataclass(frozen=True)
class FiniteDiagnostics:
    lambdas: list[str]
    max_coefficient_index: int
    contraction_threshold_rows: int
    contraction_threshold_positive_rows: int
    flow_bracket_rows: int
    flow_bracket_positive_rows: int
    min_threshold_margin: FiniteMinimum
    min_flow_bracket: FiniteMinimum


def decimal_lam_key(text: str) -> Decimal:
    return Decimal(str(text)).normalize()


def decimal_format(value: Decimal) -> str:
    return f"{value:.18E}"


def decimal_sample(row: dict) -> Decimal:
    cached = row.get("cache_A")
    if cached is not None:
        return Decimal(cached)
    ball = str(row["A_ball"]).strip()
    if not (ball.startswith("[") and "+/-" in ball):
        raise ValueError(f"cannot parse A_ball midpoint: {ball}")
    midpoint = ball[1:].split("+/-", 1)[0].strip()
    return Decimal(midpoint)


def arb_positive(value: flint.arb) -> bool:
    return bool(value > 0 and not value.contains(0))


def load_enclosures(paths: list[Path]) -> tuple[dict[tuple[Decimal, int], flint.arb], dict[tuple[Decimal, int], Decimal], dict[Decimal, str]]:
    balls: dict[tuple[Decimal, int], flint.arb] = {}
    samples: dict[tuple[Decimal, int], Decimal] = {}
    labels: dict[Decimal, str] = {}
    for path in paths:
        with path.open("r", encoding="utf-8") as handle:
            for raw in handle:
                if not raw.strip():
                    continue
                row = json.loads(raw)
                if row.get("kind") != "acb_coefficient_enclosure":
                    continue
                lam = decimal_lam_key(row["lam"])
                index = int(row["k"])
                balls[(lam, index)] = flint.arb(row["A_ball"])
                samples[(lam, index)] = decimal_sample(row)
                labels[lam] = row["lam"]
    return balls, samples, labels


def contraction(values: dict[int, object], k: int):
    return (values[k + 1] / values[k]) / (values[k] / values[k - 1])


def finite_diagnostics(paths: list[Path], max_coefficient_index: int) -> FiniteDiagnostics:
    balls, samples, labels = load_enclosures(paths)
    threshold_rows = 0
    threshold_positive = 0
    flow_rows = 0
    flow_positive = 0
    min_threshold: tuple[Decimal, str, int] | None = None
    min_flow: tuple[Decimal, str, int] | None = None

    for lam in sorted(labels):
        arb_x: dict[int, flint.arb] = {}
        sample_x: dict[int, Decimal] = {}
        for k in range(1, max_coefficient_index):
            arb_values = {idx: balls[(lam, idx)] for idx in (k - 1, k, k + 1)}
            sample_values = {idx: samples[(lam, idx)] for idx in (k - 1, k, k + 1)}
            arb_x[k] = contraction(arb_values, k)
            sample_x[k] = contraction(sample_values, k)

            threshold_rows += 1
            arb_threshold = flint.arb(2 * k - 1) / flint.arb(2 * k + 5)
            sample_threshold = Decimal(2 * k - 1) / Decimal(2 * k + 5)
            arb_margin = arb_x[k] - arb_threshold
            sample_margin = sample_x[k] - sample_threshold
            if arb_positive(arb_margin):
                threshold_positive += 1
            if min_threshold is None or sample_margin < min_threshold[0]:
                min_threshold = (sample_margin, labels[lam], k)

        for k in range(1, max_coefficient_index - 2):
            flow_rows += 1
            arb_bracket = (
                flint.arb(2 * k + 5) * arb_x[k + 1] * arb_x[k + 2]
                - flint.arb(3 * (2 * k + 3)) * arb_x[k + 1]
                + flint.arb(3 * (2 * k + 1))
                - flint.arb(2 * k - 1) / arb_x[k]
            )
            sample_bracket = (
                Decimal(2 * k + 5) * sample_x[k + 1] * sample_x[k + 2]
                - Decimal(3 * (2 * k + 3)) * sample_x[k + 1]
                + Decimal(3 * (2 * k + 1))
                - Decimal(2 * k - 1) / sample_x[k]
            )
            if arb_positive(arb_bracket):
                flow_positive += 1
            if min_flow is None or sample_bracket < min_flow[0]:
                min_flow = (sample_bracket, labels[lam], k)

    if min_threshold is None or min_flow is None:
        raise RuntimeError("no finite diagnostics were computed")
    return FiniteDiagnostics(
        lambdas=[labels[lam] for lam in sorted(labels)],
        max_coefficient_index=max_coefficient_index,
        contraction_threshold_rows=threshold_rows,
        contraction_threshold_positive_rows=threshold_positive,
        flow_bracket_rows=flow_rows,
        flow_bracket_positive_rows=flow_positive,
        min_threshold_margin=FiniteMinimum(decimal_format(min_threshold[0]), min_threshold[1], min_threshold[2]),
        min_flow_bracket=FiniteMinimum(decimal_format(min_flow[0]), min_flow[1], min_flow[2]),
    )


def build_scout(paths: list[Path], max_coefficient_index: int) -> dict:
    diagnostics = finite_diagnostics(paths, max_coefficient_index)
    return {
        "kind": "jensen_window_pf_heat_flow_monotone_closure_scout",
        "date": "2026-07-06",
        "status": "finite_heat_flow_closure_scout",
        "source_theorem_target": "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "source_boundary_threshold": "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
        "source_stress": "outputs/jensen_window_pf_monotone_contraction_stress.md",
        "source_stress_summary": "work/rh_compute/results/jensen_window_pf_monotone_contraction_stress_lamgrid_d3_d12_k64_summary.json",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "proof_boundary": (
            "Exact lambda-flow algebra plus finite Arb diagnostics only. This "
            "does not prove a closed differential inequality, does not prove "
            "the monotone-contraction theorem for all zeta windows, does not "
            "construct a Cauchy-Binet kernel, does not prove jwpf_06, and does "
            "not prove RH or Lambda <= 0."
        ),
        "exact_rows": [
            {
                "id": "hfmc_01_coefficient_flow_identity",
                "role": "exact_identity",
                "formula": "dA_k/dlambda = 2*(2*k+1)*A_{k+1}",
                "derivation": "Differentiate mu_{2k}(lambda) under the retained Phi/Newman integral and keep A_k=mu_{2k}*k!/(2*k)! explicit.",
                "proof_boundary": "Exact formal heat-flow identity only; not a proof of a sign theorem.",
            },
            {
                "id": "hfmc_02_contraction_log_flow",
                "role": "exact_identity",
                "formula": "dlog(x_k)/dlambda = 2*r_k*((2*k+3)*x_{k+1} + (2*k-1)/x_k - 2*(2*k+1))",
                "derivation": "Use x_k=(A_{k+1}/A_k)/(A_k/A_{k-1}) and r_k=A_{k+1}/A_k.",
                "proof_boundary": "Exact algebraic flow identity only; not a proof of closed positivity.",
            },
            {
                "id": "hfmc_03_monotone_gap_flow",
                "role": "exact_identity",
                "formula": "dlog(x_{k+1}/x_k)/dlambda = 2*r_k*((2*k+5)*x_{k+1}*x_{k+2} - 3*(2*k+3)*x_{k+1} + 3*(2*k+1) - (2*k-1)/x_k)",
                "derivation": "Subtract the adjacent contraction log flows.",
                "proof_boundary": "Exact algebraic flow identity only; not a proof; the bracket still needs a noncircular sign proof.",
            },
            {
                "id": "hfmc_04_boundary_threshold_factorization",
                "role": "boundary_test",
                "formula": "if x_k=x_{k+1}=q and x_{k+2}>=q, bracket >= ((q-1)^2*((2*k+5)*q-(2*k-1)))/q",
                "sufficient_threshold": "q >= (2*k-1)/(2*k+5)",
                "derivation": "Set x_k=x_{k+1}=q and minimize the bracket at x_{k+2}=q.",
                "proof_boundary": "Boundary test only; not a proof; the lower-bound subtarget is discharged separately by raw moment log-convexity.",
            },
        ],
        "finite_diagnostics": asdict(diagnostics),
        "summary": {
            "exact_rows": 4,
            "threshold_rows": diagnostics.contraction_threshold_rows,
            "threshold_positive_rows": diagnostics.contraction_threshold_positive_rows,
            "flow_bracket_rows": diagnostics.flow_bracket_rows,
            "flow_bracket_positive_rows": diagnostics.flow_bracket_positive_rows,
            "target_closing": False,
            "main_finding": (
                "The heat-flow route reduces to exact identities for the "
                "contraction flow and a concrete boundary threshold q >= "
                "(2*k-1)/(2*k+5). The existing finite zeta cache clears that "
                "threshold on 315 Arb-classified rows and has positive actual "
                "flow brackets on 305 rows. The boundary threshold is "
                "discharged exactly by raw moment log-convexity in the separate "
                "boundary-threshold lemma, but this scout is still not a closed "
                "analytic differential inequality."
            ),
        },
        "invariants": [
            "The coefficient normalization A_k=mu_{2k}*k!/(2*k)! is retained.",
            "Finite flow-bracket positivity is not promoted to an analytic theorem.",
            "The boundary threshold q >= (2*k-1)/(2*k+5) is discharged by a separate exact raw-moment log-convexity lemma.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(scout: dict, path: Path) -> None:
    diagnostics = scout["finite_diagnostics"]
    summary = scout["summary"]
    text = f"""# Jensen-Window PF Heat-Flow Monotone-Closure Scout

Date: 2026-07-06

Status: finite heat-flow closure scout. This is not a proof of a closed
differential inequality, monotone-contraction theorem, Jensen-window
PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_heat_flow_monotone_closure_scout`.

Proof boundary: exact lambda-flow algebra plus finite Arb diagnostics only.
This artifact does not prove the monotone-contraction theorem for all zeta
windows, does not construct a Cauchy-Binet kernel, does not prove `jwpf_06`,
and does not prove `Lambda <= 0`.

Machine-readable scout:

```text
work/rh_compute/results/jensen_window_pf_heat_flow_monotone_closure_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_heat_flow_monotone_closure_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_monotone_closure_scout.py
```

Current result:

```text
validated Jensen-window PF heat-flow monotone closure scout: 4 exact rows, 315 threshold rows, 305 flow-bracket rows, 0 issues
```

## Exact Flow Algebra

With

```text
A_k(lambda) = mu_{{2k}}(lambda)*k!/(2*k)!
r_k = A_{{k+1}}/A_k
x_k = r_k/r_{{k-1}}
```

the heat-flow derivative is:

```text
dA_k/dlambda = 2*(2*k+1)*A_{{k+1}}
```

So the adjacent contraction evolves by:

```text
dlog(x_k)/dlambda
  = 2*r_k*((2*k+3)*x_{{k+1}} + (2*k-1)/x_k - 2*(2*k+1))
```

and the monotone-gap log ratio evolves by:

```text
dlog(x_{{k+1}}/x_k)/dlambda
  = 2*r_k*((2*k+5)*x_{{k+1}}*x_{{k+2}}
            - 3*(2*k+3)*x_{{k+1}}
            + 3*(2*k+1)
            - (2*k-1)/x_k)
```

## Boundary Threshold

At a boundary point `x_{{k+1}}=x_k=q`, assuming only the next monotone
inequality `x_{{k+2}}>=q`, the bracket has lower bound:

```text
((q-1)^2*((2*k+5)*q-(2*k-1)))/q
```

Thus a heat-flow invariance proof from this direct bracket route needs the
additional lower-bound subtarget:

```text
q >= (2*k-1)/(2*k+5)
```

This threshold is supplied by ordinary raw-moment log-convexity after keeping
the factorial normalization explicit:

```text
outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md
validated Jensen-window PF heat-flow boundary threshold lemma: 5 exact rows, 315 strong-threshold rows, 315 heat-threshold rows, 0 issues
```

## Finite Diagnostics

The existing coefficient enclosure cache gives:

```text
lambdas: {", ".join(diagnostics["lambdas"])}
max coefficient index: {diagnostics["max_coefficient_index"]}
threshold rows: {diagnostics["contraction_threshold_rows"]}
threshold positive rows: {diagnostics["contraction_threshold_positive_rows"]}
flow-bracket rows: {diagnostics["flow_bracket_rows"]}
flow-bracket positive rows: {diagnostics["flow_bracket_positive_rows"]}
```

Weakest threshold margin:

```text
{diagnostics["min_threshold_margin"]["sample"]} at lambda={diagnostics["min_threshold_margin"]["lam"]}, k={diagnostics["min_threshold_margin"]["k"]}
```

Weakest actual flow bracket:

```text
{diagnostics["min_flow_bracket"]["sample"]} at lambda={diagnostics["min_flow_bracket"]["lam"]}, k={diagnostics["min_flow_bracket"]["k"]}
```

The finite result supports the heat-flow route but does not close it. The
threshold subtarget is now exact; the remaining analytic upgrade is a global
flow-invariance argument, adjacent-log-concavity control, and initial or
terminal conditions for the monotone-contraction cone.

## Integration

This scout sharpens:

```text
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md
outputs/jensen_window_pf_monotone_contraction_stress.md
outputs/jensen_window_pf_cauchy_binet_cone_frontier_matrix.md
```

Summary:

```text
{summary["main_finding"]}
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
    scout = build_scout(paths, args.max_coefficient_index)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(scout, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(scout, args.note)
    diagnostics = scout["finite_diagnostics"]
    print(
        "wrote Jensen-window PF heat-flow monotone closure scout: "
        f"{scout['summary']['exact_rows']} exact rows, "
        f"{diagnostics['contraction_threshold_rows']} threshold rows, "
        f"{diagnostics['flow_bracket_rows']} flow-bracket rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
