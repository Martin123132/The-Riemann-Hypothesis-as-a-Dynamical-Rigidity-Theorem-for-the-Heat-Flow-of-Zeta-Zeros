#!/usr/bin/env python3
"""Build a finite large-negative-lambda ratio-cone prefix scout."""

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
    REPO_ROOT,
    arb_positive,
    contraction,
    decimal_format,
    load_enclosures,
)


getcontext().prec = 100

DEFAULT_ENCLOSURE_JSONL = (
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k22.jsonl",
)
DEFAULT_ENCLOSURE_SUMMARY = (
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k22_summary.json"
)
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.md"


@dataclass(frozen=True)
class ConeMinimum:
    sample: str
    lam: str
    k: int


@dataclass(frozen=True)
class PrefixDiagnostics:
    lambdas: list[str]
    coefficient_rows: int
    coefficient_k_min: int
    coefficient_k_max: int
    lower_upper_k_max: int
    monotone_k_max: int
    lower_wall_rows: int
    lower_wall_positive_rows: int
    upper_wall_rows: int
    upper_wall_positive_rows: int
    monotone_gap_rows: int
    monotone_gap_positive_rows: int
    min_lower_wall_margin: ConeMinimum
    min_upper_wall_margin: ConeMinimum
    min_monotone_gap: ConeMinimum


def finite_diagnostics(
    paths: list[Path],
    lower_upper_k_max: int,
    monotone_k_max: int,
) -> PrefixDiagnostics:
    balls, samples, labels = load_enclosures(paths)
    lambdas = sorted(labels)
    if not lambdas:
        raise RuntimeError("no coefficient enclosures found")

    required_k_max = max(lower_upper_k_max + 1, monotone_k_max + 2)
    missing = [
        (labels[lam], k)
        for lam in lambdas
        for k in range(required_k_max + 1)
        if (lam, k) not in balls
    ]
    if missing:
        raise RuntimeError(f"missing coefficient enclosures: {missing[:5]}")

    lower_rows = 0
    lower_positive = 0
    upper_rows = 0
    upper_positive = 0
    monotone_rows = 0
    monotone_positive = 0
    min_lower: tuple[Decimal, str, int] | None = None
    min_upper: tuple[Decimal, str, int] | None = None
    min_mono: tuple[Decimal, str, int] | None = None

    for lam in lambdas:
        arb_x: dict[int, flint.arb] = {}
        sample_x: dict[int, Decimal] = {}
        for k in range(1, lower_upper_k_max + 1):
            arb_values = {idx: balls[(lam, idx)] for idx in (k - 1, k, k + 1)}
            sample_values = {idx: samples[(lam, idx)] for idx in (k - 1, k, k + 1)}
            arb_x[k] = contraction(arb_values, k)
            sample_x[k] = contraction(sample_values, k)

            lower_rows += 1
            arb_lower = flint.arb(2 * k - 1) / flint.arb(2 * k + 1)
            sample_lower = Decimal(2 * k - 1) / Decimal(2 * k + 1)
            arb_lower_margin = arb_x[k] - arb_lower
            sample_lower_margin = sample_x[k] - sample_lower
            if arb_positive(arb_lower_margin):
                lower_positive += 1
            if min_lower is None or sample_lower_margin < min_lower[0]:
                min_lower = (sample_lower_margin, labels[lam], k)

            upper_rows += 1
            arb_upper_margin = flint.arb(1) - arb_x[k]
            sample_upper_margin = Decimal(1) - sample_x[k]
            if arb_positive(arb_upper_margin):
                upper_positive += 1
            if min_upper is None or sample_upper_margin < min_upper[0]:
                min_upper = (sample_upper_margin, labels[lam], k)

        for k in range(1, monotone_k_max + 1):
            monotone_rows += 1
            arb_gap = arb_x[k + 1] - arb_x[k]
            sample_gap = sample_x[k + 1] - sample_x[k]
            if arb_positive(arb_gap):
                monotone_positive += 1
            if min_mono is None or sample_gap < min_mono[0]:
                min_mono = (sample_gap, labels[lam], k)

    if min_lower is None or min_upper is None or min_mono is None:
        raise RuntimeError("no finite diagnostics were computed")

    coefficient_indices = [index for (_lam, index) in balls]
    return PrefixDiagnostics(
        lambdas=[labels[lam] for lam in lambdas],
        coefficient_rows=len(balls),
        coefficient_k_min=min(coefficient_indices),
        coefficient_k_max=max(coefficient_indices),
        lower_upper_k_max=lower_upper_k_max,
        monotone_k_max=monotone_k_max,
        lower_wall_rows=lower_rows,
        lower_wall_positive_rows=lower_positive,
        upper_wall_rows=upper_rows,
        upper_wall_positive_rows=upper_positive,
        monotone_gap_rows=monotone_rows,
        monotone_gap_positive_rows=monotone_positive,
        min_lower_wall_margin=ConeMinimum(decimal_format(min_lower[0]), min_lower[1], min_lower[2]),
        min_upper_wall_margin=ConeMinimum(decimal_format(min_upper[0]), min_upper[1], min_upper[2]),
        min_monotone_gap=ConeMinimum(decimal_format(min_mono[0]), min_mono[1], min_mono[2]),
    )


def build_scout(
    paths: list[Path],
    enclosure_summary: Path,
    lower_upper_k_max: int,
    monotone_k_max: int,
) -> dict:
    diagnostics = finite_diagnostics(paths, lower_upper_k_max, monotone_k_max)
    run_id = paths[0].stem if len(paths) == 1 else "acb_enclosures_negative_lambda_cone_entry"
    return {
        "kind": "jensen_window_pf_negative_lambda_cone_entry_prefix_scout",
        "date": "2026-07-06",
        "status": "finite_negative_lambda_prefix_certificate",
        "source_target": "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "source_sign_scout": "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "enclosure_summary": enclosure_summary.relative_to(REPO_ROOT).as_posix(),
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_cone_entry_prefix_scout.py",
        "coefficient_reproduction_command": (
            "python work/rh_compute/scripts/acb_coefficient_enclosures.py "
            f"--lambdas=-25,-50,-100 --k-min 0 --k-max {diagnostics.coefficient_k_max} --n-sum 50 "
            "--cutoff 6 --dps 160 --digits 80 --abs-tol 1e-110 "
            "--constant-terms 30 --output-dir work/rh_compute/results "
            f"--run-id {run_id} --overwrite"
        ),
        "proof_boundary": (
            "Finite large-negative-lambda prefix certificate only. It certifies "
            "the ratio-cone inequalities on the listed lambda grid and finite "
            "prefix, but it does not prove an infinite cone-entry theorem, does "
            "not prove a collared finite flow theorem, does not control all k, "
            "does not prove monotone contractions for all zeta windows, does "
            "not prove jwpf_06, and leaves Lambda <= 0 unsettled."
        ),
        "prefix_statement": {
            "lambdas": diagnostics.lambdas,
            "coefficient_range": f"A_0..A_{diagnostics.coefficient_k_max}",
            "lower_upper_walls": f"1 <= k <= {diagnostics.lower_upper_k_max}",
            "monotone_wall": f"1 <= k <= {diagnostics.monotone_k_max}",
            "ratio": "x_k=(A_{k+1}/A_k)/(A_k/A_{k-1})",
            "lower_wall": "x_k >= (2*k-1)/(2*k+1)",
            "upper_wall": "x_k <= 1",
            "monotone_wall_formula": "x_{k+1} >= x_k",
        },
        "finite_rows": [
            {
                "id": "nlcep_01_acb_coefficient_enclosures",
                "role": "finite_coefficient_input",
                "claim": f"ACB/Arb enclosures cover A_0..A_{diagnostics.coefficient_k_max} for lambda=-25,-50,-100.",
                "rows": diagnostics.coefficient_rows,
                "proof_boundary": "Finite coefficient enclosures only.",
            },
            {
                "id": "nlcep_02_lower_wall_prefix",
                "role": "finite_lower_wall_certificate",
                "claim": "The checked contractions satisfy x_k >= (2*k-1)/(2*k+1).",
                "rows": diagnostics.lower_wall_rows,
                "positive_rows": diagnostics.lower_wall_positive_rows,
                "proof_boundary": "Finite prefix only; no all-k tail theorem.",
            },
            {
                "id": "nlcep_03_upper_wall_prefix",
                "role": "finite_upper_wall_certificate",
                "claim": "The checked contractions satisfy x_k <= 1.",
                "rows": diagnostics.upper_wall_rows,
                "positive_rows": diagnostics.upper_wall_positive_rows,
                "proof_boundary": "Finite prefix only; no all-k tail theorem.",
            },
            {
                "id": "nlcep_04_monotone_wall_prefix",
                "role": "finite_monotone_wall_certificate",
                "claim": "The checked contractions satisfy x_{k+1} >= x_k.",
                "rows": diagnostics.monotone_gap_rows,
                "positive_rows": diagnostics.monotone_gap_positive_rows,
                "proof_boundary": "Finite prefix only; no collared finite flow theorem.",
            },
        ],
        "finite_diagnostics": asdict(diagnostics),
        "summary": {
            "coefficient_rows": diagnostics.coefficient_rows,
            "lower_wall_rows": diagnostics.lower_wall_rows,
            "lower_wall_positive_rows": diagnostics.lower_wall_positive_rows,
            "upper_wall_rows": diagnostics.upper_wall_rows,
            "upper_wall_positive_rows": diagnostics.upper_wall_positive_rows,
            "monotone_gap_rows": diagnostics.monotone_gap_rows,
            "monotone_gap_positive_rows": diagnostics.monotone_gap_positive_rows,
            "ready_to_apply_rows": 0,
            "target_closing": False,
            "main_finding": (
                "The actual ACB-enclosed zeta coefficients satisfy the ratio "
                "cone on a finite large-negative-lambda prefix. This supports "
                "the cone-entry route beyond formal Taylor signs, but the open "
                "bridge remains an infinite or collared finite theorem with "
                "tail control."
            ),
        },
        "invariants": [
            "The scout uses actual A_k(lambda) enclosures, not only the Taylor expansion.",
            "The certificate is finite in lambda and k.",
            "No row is ready_to_apply for the cone-entry theorem.",
            "The infinite or collared finite cone-entry target remains open.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(scout: dict, path: Path) -> None:
    diagnostics = scout["finite_diagnostics"]
    summary = scout["summary"]
    result_json = scout.get("result_json", "work/rh_compute/results/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.json")
    note_ref = scout.get("note", path.relative_to(REPO_ROOT).as_posix())
    text = f"""# Jensen-Window PF Negative-Lambda Cone-Entry Prefix Scout

Date: 2026-07-06

Status: finite negative-lambda prefix certificate. This is not a proof of zeta
cone entry, monotone contractions, Jensen-window PF-infinity, RH, or the
Newman-direction goal.

Artifact kind: `jensen_window_pf_negative_lambda_cone_entry_prefix_scout`.

Proof boundary: this artifact certifies only a finite large-negative-lambda
prefix. It does not prove an infinite cone-entry theorem, does not prove a
collared finite flow theorem, does not control all `k`, does not prove
`jwpf_06`, and does not settle `Lambda <= 0`.

Machine-readable certificate:

```text
{result_json}
```

Coefficient enclosures:

```text
{scout['enclosure_jsonl'][0]}
{scout['enclosure_summary']}
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_cone_entry_prefix_scout.py --target {result_json} --note {note_ref}
```

Current result:

```text
validated Jensen-window PF negative-lambda cone-entry prefix scout: {summary['coefficient_rows']} coefficient rows, {summary['lower_wall_rows']} lower-wall rows, {summary['upper_wall_rows']} upper-wall rows, {summary['monotone_gap_rows']} monotone-gap rows, 0 issues
```

## Scope

The checked ratios are:

```text
x_k=(A_{{k+1}}/A_k)/(A_k/A_{{k-1}})
(2*k-1)/(2*k+1) <= x_k <= 1
x_{{k+1}} >= x_k
```

The finite prefix is:

```text
lambdas: {', '.join(diagnostics['lambdas'])}
coefficient range: A_{diagnostics['coefficient_k_min']}..A_{diagnostics['coefficient_k_max']}
lower and upper walls: 1 <= k <= {diagnostics['lower_upper_k_max']}
monotone wall: 1 <= k <= {diagnostics['monotone_k_max']}
```

The coefficient enclosures were generated by:

```text
{scout['coefficient_reproduction_command']}
```

## Certified Margins

```text
lower wall rows: {diagnostics['lower_wall_positive_rows']} / {diagnostics['lower_wall_rows']}
upper wall rows: {diagnostics['upper_wall_positive_rows']} / {diagnostics['upper_wall_rows']}
monotone gap rows: {diagnostics['monotone_gap_positive_rows']} / {diagnostics['monotone_gap_rows']}

minimum lower-wall margin:
  {diagnostics['min_lower_wall_margin']['sample']} at lambda={diagnostics['min_lower_wall_margin']['lam']}, k={diagnostics['min_lower_wall_margin']['k']}

minimum upper-wall margin:
  {diagnostics['min_upper_wall_margin']['sample']} at lambda={diagnostics['min_upper_wall_margin']['lam']}, k={diagnostics['min_upper_wall_margin']['k']}

minimum monotone gap:
  {diagnostics['min_monotone_gap']['sample']} at lambda={diagnostics['min_monotone_gap']['lam']}, k={diagnostics['min_monotone_gap']['k']}
```

## Consequence

This moves the cone-entry route from purely formal fixed-k asymptotics to a
finite interval-certified large-negative-lambda prefix. The remaining theorem
work is unchanged in kind: prove an all-`k` tail/uniform statement, or prove a
finite-prefix collar theorem strong enough to feed the conditional heat-flow
ratio-cone invariance lemma.

Related artifacts:

```text
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md
outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md
```
"""
    path.write_text(text, encoding="utf-8")


def write_json(scout: dict, path: Path) -> None:
    path.write_text(json.dumps(scout, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--enclosure-jsonl", type=Path, nargs="+", default=list(DEFAULT_ENCLOSURE_JSONL))
    parser.add_argument("--enclosure-summary", type=Path, default=DEFAULT_ENCLOSURE_SUMMARY)
    parser.add_argument("--lower-upper-k-max", type=int, default=21)
    parser.add_argument("--monotone-k-max", type=int, default=20)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    paths = [path if path.is_absolute() else REPO_ROOT / path for path in args.enclosure_jsonl]
    summary = args.enclosure_summary if args.enclosure_summary.is_absolute() else REPO_ROOT / args.enclosure_summary
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    scout = build_scout(paths, summary, args.lower_upper_k_max, args.monotone_k_max)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    scout["result_json"] = out_json.relative_to(REPO_ROOT).as_posix()
    scout["note"] = note.relative_to(REPO_ROOT).as_posix()
    write_json(scout, out_json)
    write_note(scout, note)
    print(
        "wrote Jensen-window PF negative-lambda cone-entry prefix scout: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
