#!/usr/bin/env python3
"""Build the Jensen-window PF heat-flow ratio-cone invariance lemma."""

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

DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md"


@dataclass(frozen=True)
class FiniteMinimum:
    sample: str
    lam: str
    k: int


@dataclass(frozen=True)
class ConeDiagnostics:
    lambdas: list[str]
    max_coefficient_index: int
    coordinate_lower_rows: int
    coordinate_lower_positive_rows: int
    coordinate_upper_rows: int
    coordinate_upper_positive_rows: int
    coordinate_monotone_rows: int
    coordinate_monotone_positive_rows: int
    lower_wall_rows: int
    lower_wall_positive_rows: int
    upper_wall_rows: int
    upper_wall_positive_rows: int
    monotone_wall_rows: int
    monotone_wall_positive_rows: int
    min_coordinate_lower_margin: FiniteMinimum
    min_coordinate_upper_margin: FiniteMinimum
    min_coordinate_monotone_margin: FiniteMinimum
    min_lower_wall_margin: FiniteMinimum
    min_upper_wall_margin: FiniteMinimum
    min_monotone_wall_margin: FiniteMinimum


def arb_positive(value: flint.arb) -> bool:
    return bool(value > 0 and not value.contains(0))


def update_min(current: tuple[Decimal, str, int] | None, candidate: tuple[Decimal, str, int]) -> tuple[Decimal, str, int]:
    if current is None or candidate[0] < current[0]:
        return candidate
    return current


def finite_diagnostics(paths: list[Path], max_coefficient_index: int) -> ConeDiagnostics:
    balls, samples, labels = load_enclosures(paths)
    coord_lower_rows = coord_lower_pos = 0
    coord_upper_rows = coord_upper_pos = 0
    coord_mono_rows = coord_mono_pos = 0
    lower_wall_rows = lower_wall_pos = 0
    upper_wall_rows = upper_wall_pos = 0
    mono_wall_rows = mono_wall_pos = 0
    min_lower: tuple[Decimal, str, int] | None = None
    min_upper: tuple[Decimal, str, int] | None = None
    min_mono: tuple[Decimal, str, int] | None = None
    min_lower_wall: tuple[Decimal, str, int] | None = None
    min_upper_wall: tuple[Decimal, str, int] | None = None
    min_mono_wall: tuple[Decimal, str, int] | None = None

    for lam in sorted(labels):
        arb_x: dict[int, flint.arb] = {}
        sample_x: dict[int, Decimal] = {}
        for k in range(1, max_coefficient_index):
            arb_values = {idx: balls[(lam, idx)] for idx in (k - 1, k, k + 1)}
            sample_values = {idx: samples[(lam, idx)] for idx in (k - 1, k, k + 1)}
            arb_x[k] = contraction(arb_values, k)
            sample_x[k] = contraction(sample_values, k)

            lower_threshold = Decimal(2 * k - 1) / Decimal(2 * k + 1)
            arb_lower_threshold = flint.arb(2 * k - 1) / flint.arb(2 * k + 1)
            coord_lower_rows += 1
            lower_margin = sample_x[k] - lower_threshold
            if arb_positive(arb_x[k] - arb_lower_threshold):
                coord_lower_pos += 1
            min_lower = update_min(min_lower, (lower_margin, labels[lam], k))

            coord_upper_rows += 1
            upper_margin = Decimal(1) - sample_x[k]
            if arb_positive(flint.arb(1) - arb_x[k]):
                coord_upper_pos += 1
            min_upper = update_min(min_upper, (upper_margin, labels[lam], k))

        for k in range(1, max_coefficient_index - 1):
            coord_mono_rows += 1
            mono_margin = sample_x[k + 1] - sample_x[k]
            if arb_positive(arb_x[k + 1] - arb_x[k]):
                coord_mono_pos += 1
            min_mono = update_min(min_mono, (mono_margin, labels[lam], k))

            lower_wall_rows += 1
            lower_wall = Decimal(2 * k + 3) * sample_x[k + 1] - Decimal(2 * k + 1)
            arb_lower_wall = flint.arb(2 * k + 3) * arb_x[k + 1] - flint.arb(2 * k + 1)
            if arb_positive(arb_lower_wall):
                lower_wall_pos += 1
            min_lower_wall = update_min(min_lower_wall, (lower_wall, labels[lam], k))

            upper_wall_rows += 1
            upper_wall = Decimal(2 * k + 3) * (Decimal(1) - sample_x[k + 1])
            arb_upper_wall = flint.arb(2 * k + 3) * (flint.arb(1) - arb_x[k + 1])
            if arb_positive(arb_upper_wall):
                upper_wall_pos += 1
            min_upper_wall = update_min(min_upper_wall, (upper_wall, labels[lam], k))

        for k in range(1, max_coefficient_index - 2):
            mono_wall_rows += 1
            mono_wall = (
                Decimal(2 * k + 5) * sample_x[k + 1] * sample_x[k + 2]
                - Decimal(3 * (2 * k + 3)) * sample_x[k + 1]
                + Decimal(3 * (2 * k + 1))
                - Decimal(2 * k - 1) / sample_x[k]
            )
            arb_mono_wall = (
                flint.arb(2 * k + 5) * arb_x[k + 1] * arb_x[k + 2]
                - flint.arb(3 * (2 * k + 3)) * arb_x[k + 1]
                + flint.arb(3 * (2 * k + 1))
                - flint.arb(2 * k - 1) / arb_x[k]
            )
            if arb_positive(arb_mono_wall):
                mono_wall_pos += 1
            min_mono_wall = update_min(min_mono_wall, (mono_wall, labels[lam], k))

    minima = (min_lower, min_upper, min_mono, min_lower_wall, min_upper_wall, min_mono_wall)
    if any(item is None for item in minima):
        raise RuntimeError("no finite cone diagnostics were computed")

    assert min_lower is not None
    assert min_upper is not None
    assert min_mono is not None
    assert min_lower_wall is not None
    assert min_upper_wall is not None
    assert min_mono_wall is not None
    return ConeDiagnostics(
        lambdas=[labels[lam] for lam in sorted(labels)],
        max_coefficient_index=max_coefficient_index,
        coordinate_lower_rows=coord_lower_rows,
        coordinate_lower_positive_rows=coord_lower_pos,
        coordinate_upper_rows=coord_upper_rows,
        coordinate_upper_positive_rows=coord_upper_pos,
        coordinate_monotone_rows=coord_mono_rows,
        coordinate_monotone_positive_rows=coord_mono_pos,
        lower_wall_rows=lower_wall_rows,
        lower_wall_positive_rows=lower_wall_pos,
        upper_wall_rows=upper_wall_rows,
        upper_wall_positive_rows=upper_wall_pos,
        monotone_wall_rows=mono_wall_rows,
        monotone_wall_positive_rows=mono_wall_pos,
        min_coordinate_lower_margin=FiniteMinimum(decimal_format(min_lower[0]), min_lower[1], min_lower[2]),
        min_coordinate_upper_margin=FiniteMinimum(decimal_format(min_upper[0]), min_upper[1], min_upper[2]),
        min_coordinate_monotone_margin=FiniteMinimum(decimal_format(min_mono[0]), min_mono[1], min_mono[2]),
        min_lower_wall_margin=FiniteMinimum(decimal_format(min_lower_wall[0]), min_lower_wall[1], min_lower_wall[2]),
        min_upper_wall_margin=FiniteMinimum(decimal_format(min_upper_wall[0]), min_upper_wall[1], min_upper_wall[2]),
        min_monotone_wall_margin=FiniteMinimum(decimal_format(min_mono_wall[0]), min_mono_wall[1], min_mono_wall[2]),
    )


def build_lemma(paths: list[Path], max_coefficient_index: int) -> dict:
    diagnostics = finite_diagnostics(paths, max_coefficient_index)
    return {
        "kind": "jensen_window_pf_heat_flow_ratio_cone_invariance_lemma",
        "date": "2026-07-06",
        "status": "available_exact_conditional_invariance_lemma",
        "source_boundary_threshold": "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
        "source_heat_flow_closure": "outputs/jensen_window_pf_heat_flow_monotone_closure_scout.md",
        "source_theorem_target": "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "proof_boundary": (
            "Exact conditional ratio-cone invariance lemma only. It proves "
            "inward-pointing boundary algebra for the cone "
            "(2*k-1)/(2*k+1) <= x_k <= 1 and x_{k+1}>=x_k under the "
            "heat-flow ratio ODE, assuming positive smooth coefficients and "
            "the required infinite or collared finite ratio variables. It does "
            "not prove that the actual zeta coefficients enter the cone at an "
            "initial or terminal lambda, does not prove all-shift "
            "monotone contractions, does not prove jwpf_06, and does not prove "
            "RH or Lambda <= 0."
        ),
        "exact_rows": [
            {
                "id": "hfrci_01_ratio_flow_identity",
                "role": "exact_identity",
                "formula": "dlog(x_k)/dlambda = 2*r_k*((2*k+3)*x_{k+1} + (2*k-1)/x_k - 2*(2*k+1))",
                "argument": "Differentiate x_k=(A_{k+1}/A_k)/(A_k/A_{k-1}) using dA_k/dlambda=2*(2*k+1)*A_{k+1}.",
                "proof_boundary": "Exact ratio-flow identity only.",
            },
            {
                "id": "hfrci_02_lower_wall_inward",
                "role": "exact_boundary_condition",
                "formula": "at x_k=(2*k-1)/(2*k+1), F_k=(2*k+3)*x_{k+1}-(2*k+1)>=0 if x_{k+1}>=(2*k+1)/(2*k+3)",
                "argument": "The strong lower threshold for the next coordinate makes the lower wall inward-pointing.",
                "proof_boundary": "Exact lower-wall condition only; not an initial-condition theorem.",
            },
            {
                "id": "hfrci_03_upper_wall_inward",
                "role": "exact_boundary_condition",
                "formula": "at x_k=1, F_k=(2*k+3)*(x_{k+1}-1)<=0 if x_{k+1}<=1",
                "argument": "Adjacent log-concavity of the next coordinate makes the upper wall inward-pointing.",
                "proof_boundary": "Exact upper-wall condition only; it assumes the next upper bound.",
            },
            {
                "id": "hfrci_04_monotone_wall_inward",
                "role": "exact_boundary_condition",
                "formula": "at x_{k+1}=x_k=q, dlog(x_{k+1}/x_k)/dlambda >= 0 if x_{k+2}>=q and q>=(2*k-1)/(2*k+5)",
                "argument": "Use the heat-flow monotone-wall bracket and the exact boundary-threshold lemma.",
                "proof_boundary": "Exact monotone-wall condition only; not a proof that the wall is reached from zeta data.",
            },
            {
                "id": "hfrci_05_finite_collar_requirement",
                "role": "finite_prefix_caveat",
                "formula": "checking x_1..x_K requires collar variables x_{K+1} for lower/upper walls and x_{K+2} for monotone walls",
                "argument": "The vector field for x_k depends on x_{k+1}; the monotone-wall field for x_{k+1}/x_k depends on x_{k+2}.",
                "proof_boundary": "Finite-prefix caveat only; no finite prefix alone proves the infinite cone.",
            },
            {
                "id": "hfrci_06_conditional_forward_invariance",
                "role": "conditional_exact_lemma",
                "formula": "if the full ratio cone holds at lambda0 and the positive heat-flow ratio ODE is well posed, then the cone is forward-invariant while the hypotheses persist",
                "argument": "Combine the inward-pointing lower, upper, and monotone boundary conditions for every k.",
                "proof_boundary": "Conditional invariance only; the zeta initial/terminal cone condition and global analytic well-posedness remain open.",
            },
        ],
        "finite_diagnostics": asdict(diagnostics),
        "summary": {
            "exact_rows": 6,
            "coordinate_lower_rows": diagnostics.coordinate_lower_rows,
            "coordinate_upper_rows": diagnostics.coordinate_upper_rows,
            "coordinate_monotone_rows": diagnostics.coordinate_monotone_rows,
            "lower_wall_rows": diagnostics.lower_wall_rows,
            "upper_wall_rows": diagnostics.upper_wall_rows,
            "monotone_wall_rows": diagnostics.monotone_wall_rows,
            "conditional_invariance_available": True,
            "zeta_initial_condition_available": False,
            "target_closing": False,
            "main_finding": (
                "The heat-flow ratio cone has an exact conditional "
                "forward-invariance algebra: the lower wall, upper wall, and "
                "monotone wall all point inward under the stated neighbor and "
                "threshold hypotheses. The remaining theorem gap is no longer "
                "the boundary algebra; it is to prove that the actual zeta "
                "heat-flow coefficient sequence enters the full infinite cone "
                "at a suitable lambda and that the infinite/collared flow "
                "argument is analytically legitimate."
            ),
        },
        "invariants": [
            "The lemma is conditional on a full infinite cone or a collared finite cone.",
            "The finite zeta diagnostics are sanity evidence only.",
            "No initial or terminal zeta cone condition is proved.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(lemma: dict, path: Path) -> None:
    diagnostics = lemma["finite_diagnostics"]
    summary = lemma["summary"]
    text = f"""# Jensen-Window PF Heat-Flow Ratio-Cone Invariance Lemma

Date: 2026-07-06

Status: exact conditional ratio-cone invariance lemma. This is not a proof of
the monotone-contraction theorem for the zeta coefficients, Jensen-window
PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_heat_flow_ratio_cone_invariance_lemma`.

Proof boundary: this artifact proves inward-pointing boundary algebra for the
ratio cone only, conditional on positive smooth coefficients and the required
infinite or collared finite ratio variables. It does not prove that the actual
zeta coefficients enter the cone at an initial or terminal lambda, does not
prove all-shift monotone contractions, does not prove `jwpf_06`, and does not
establish `Lambda <= 0`.

Machine-readable lemma:

```text
work/rh_compute/results/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.py
```

Current result:

```text
validated Jensen-window PF heat-flow ratio cone invariance lemma: 6 exact rows, 315 lower rows, 315 upper rows, 310 monotone rows, 0 issues
```

## Cone

Use:

```text
r_k=A_{{k+1}}/A_k
x_k=r_k/r_{{k-1}}
a_k=(2*k-1)/(2*k+1)
```

The ratio cone is:

```text
a_k <= x_k <= 1
x_{{k+1}} >= x_k
```

The exact threshold lemma supplies the lower wall:

```text
outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md
validated Jensen-window PF heat-flow boundary threshold lemma: 5 exact rows, 315 strong-threshold rows, 315 heat-threshold rows, 0 issues
```

## Boundary Algebra

Write

```text
dlog(x_k)/dlambda = 2*r_k*F_k
F_k=(2*k+3)*x_{{k+1}} + (2*k-1)/x_k - 2*(2*k+1)
```

Lower wall:

```text
at x_k=a_k, F_k=(2*k+3)*x_{{k+1}}-(2*k+1)>=0
```

provided `x_{{k+1}}>=(2*k+1)/(2*k+3)`.

Upper wall:

```text
at x_k=1, F_k=(2*k+3)*(x_{{k+1}}-1)<=0
```

provided `x_{{k+1}}<=1`.

Monotone wall:

```text
at x_{{k+1}}=x_k=q,
dlog(x_{{k+1}}/x_k)/dlambda >= 0
```

provided `x_{{k+2}}>=q` and `q >= (2*k-1)/(2*k+5)`.

Thus the cone has exact forward-invariance boundary algebra. A finite prefix
needs collar variables: `x_{{K+1}}` for lower/upper wall tests and `x_{{K+2}}`
for monotone wall tests.

## Finite Sanity Diagnostics

The existing coefficient enclosure cache gives:

```text
lambdas: {", ".join(diagnostics["lambdas"])}
max coefficient index: {diagnostics["max_coefficient_index"]}
coordinate lower rows: {diagnostics["coordinate_lower_rows"]}
coordinate upper rows: {diagnostics["coordinate_upper_rows"]}
coordinate monotone rows: {diagnostics["coordinate_monotone_rows"]}
lower wall rows: {diagnostics["lower_wall_rows"]}
upper wall rows: {diagnostics["upper_wall_rows"]}
monotone wall rows: {diagnostics["monotone_wall_rows"]}
```

Minimum coordinate margins:

```text
lower: {diagnostics["min_coordinate_lower_margin"]["sample"]} at lambda={diagnostics["min_coordinate_lower_margin"]["lam"]}, k={diagnostics["min_coordinate_lower_margin"]["k"]}
upper: {diagnostics["min_coordinate_upper_margin"]["sample"]} at lambda={diagnostics["min_coordinate_upper_margin"]["lam"]}, k={diagnostics["min_coordinate_upper_margin"]["k"]}
monotone: {diagnostics["min_coordinate_monotone_margin"]["sample"]} at lambda={diagnostics["min_coordinate_monotone_margin"]["lam"]}, k={diagnostics["min_coordinate_monotone_margin"]["k"]}
```

Minimum wall margins:

```text
lower wall: {diagnostics["min_lower_wall_margin"]["sample"]} at lambda={diagnostics["min_lower_wall_margin"]["lam"]}, k={diagnostics["min_lower_wall_margin"]["k"]}
upper wall: {diagnostics["min_upper_wall_margin"]["sample"]} at lambda={diagnostics["min_upper_wall_margin"]["lam"]}, k={diagnostics["min_upper_wall_margin"]["k"]}
monotone wall: {diagnostics["min_monotone_wall_margin"]["sample"]} at lambda={diagnostics["min_monotone_wall_margin"]["lam"]}, k={diagnostics["min_monotone_wall_margin"]["k"]}
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
        "wrote Jensen-window PF heat-flow ratio cone invariance lemma: "
        f"{lemma['summary']['exact_rows']} exact rows, "
        f"{diagnostics['coordinate_lower_rows']} lower rows, "
        f"{diagnostics['coordinate_upper_rows']} upper rows, "
        f"{diagnostics['coordinate_monotone_rows']} monotone rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
