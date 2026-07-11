#!/usr/bin/env python3
"""Build the finite-collar contract for the negative-lambda cone prefix."""

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
DEFAULT_PREFIX_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.json"
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_finite_collar_contract.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_finite_collar_contract.md"


@dataclass(frozen=True)
class CollarMinimum:
    sample: str
    lam: str
    k: int


@dataclass(frozen=True)
class CollarDiagnostics:
    lambdas: list[str]
    coefficient_k_max: int
    available_x_max: int
    lower_upper_k_max: int
    monotone_gap_k_max: int
    certified_active_depth: int
    first_collar_x: int
    second_collar_x: int
    active_lower_rows: int
    active_lower_positive_rows: int
    active_upper_rows: int
    active_upper_positive_rows: int
    active_monotone_rows: int
    active_monotone_positive_rows: int
    collar_lower_upper_rows: int
    collar_lower_upper_positive_rows: int
    collar_monotone_rows: int
    collar_monotone_positive_rows: int
    min_active_lower_margin: CollarMinimum
    min_active_upper_margin: CollarMinimum
    min_active_monotone_gap: CollarMinimum
    min_collar_lower_upper_margin: CollarMinimum
    min_collar_monotone_gap: CollarMinimum


def infer_available_ranges(paths: list[Path]) -> tuple[dict, dict, dict, int, int, int]:
    balls, samples, labels = load_enclosures(paths)
    if not labels:
        raise RuntimeError("no coefficient enclosures found")
    coefficient_k_max = max(index for (_lam, index) in balls)
    available_x_max = coefficient_k_max - 1
    lambdas = sorted(labels)

    lower_upper_k_max = 0
    for k in range(1, available_x_max + 1):
        ok = True
        for lam in lambdas:
            arb_values = {idx: balls[(lam, idx)] for idx in (k - 1, k, k + 1)}
            x = contraction(arb_values, k)
            lower = flint.arb(2 * k - 1) / flint.arb(2 * k + 1)
            if not (arb_positive(x - lower) and arb_positive(flint.arb(1) - x)):
                ok = False
                break
        if ok:
            lower_upper_k_max = k
        else:
            break

    monotone_gap_k_max = 0
    x_cache: dict[tuple[Decimal, int], flint.arb] = {}
    for lam in lambdas:
        for k in range(1, available_x_max + 1):
            arb_values = {idx: balls[(lam, idx)] for idx in (k - 1, k, k + 1)}
            x_cache[(lam, k)] = contraction(arb_values, k)
    for k in range(1, available_x_max):
        ok = True
        for lam in lambdas:
            if not arb_positive(x_cache[(lam, k + 1)] - x_cache[(lam, k)]):
                ok = False
                break
        if ok:
            monotone_gap_k_max = k
        else:
            break

    return balls, samples, labels, coefficient_k_max, lower_upper_k_max, monotone_gap_k_max


def build_diagnostics(paths: list[Path]) -> CollarDiagnostics:
    balls, samples, labels, coefficient_k_max, lower_upper_k_max, monotone_gap_k_max = infer_available_ranges(paths)
    lambdas = sorted(labels)
    available_x_max = coefficient_k_max - 1

    # For active depth K, lower/upper boundary tests need x_{K+1}, and
    # monotone-wall tests need x_{K+2}. Since the checked monotone gaps are
    # indexed by x_{k+1} - x_k, the available active depth is K <= gap_max-1.
    certified_active_depth = min(
        lower_upper_k_max - 1,
        monotone_gap_k_max - 1,
        coefficient_k_max - 3,
    )
    if certified_active_depth < 1:
        raise RuntimeError("no positive finite-collar active depth was certified")

    first_collar = certified_active_depth + 1
    second_collar = certified_active_depth + 2

    active_lower_rows = 0
    active_lower_pos = 0
    active_upper_rows = 0
    active_upper_pos = 0
    active_mono_rows = 0
    active_mono_pos = 0
    collar_lu_rows = 0
    collar_lu_pos = 0
    collar_mono_rows = 0
    collar_mono_pos = 0

    min_active_lower: tuple[Decimal, str, int] | None = None
    min_active_upper: tuple[Decimal, str, int] | None = None
    min_active_mono: tuple[Decimal, str, int] | None = None
    min_collar_lu: tuple[Decimal, str, int] | None = None
    min_collar_mono: tuple[Decimal, str, int] | None = None

    arb_x: dict[tuple[Decimal, int], flint.arb] = {}
    sample_x: dict[tuple[Decimal, int], Decimal] = {}
    for lam in lambdas:
        for k in range(1, second_collar + 1):
            arb_values = {idx: balls[(lam, idx)] for idx in (k - 1, k, k + 1)}
            sample_values = {idx: samples[(lam, idx)] for idx in (k - 1, k, k + 1)}
            arb_x[(lam, k)] = contraction(arb_values, k)
            sample_x[(lam, k)] = contraction(sample_values, k)

        for k in range(1, certified_active_depth + 1):
            active_lower_rows += 1
            lower = flint.arb(2 * k - 1) / flint.arb(2 * k + 1)
            sample_lower = Decimal(2 * k - 1) / Decimal(2 * k + 1)
            lower_margin = arb_x[(lam, k)] - lower
            sample_lower_margin = sample_x[(lam, k)] - sample_lower
            if arb_positive(lower_margin):
                active_lower_pos += 1
            if min_active_lower is None or sample_lower_margin < min_active_lower[0]:
                min_active_lower = (sample_lower_margin, labels[lam], k)

            active_upper_rows += 1
            upper_margin = flint.arb(1) - arb_x[(lam, k)]
            sample_upper_margin = Decimal(1) - sample_x[(lam, k)]
            if arb_positive(upper_margin):
                active_upper_pos += 1
            if min_active_upper is None or sample_upper_margin < min_active_upper[0]:
                min_active_upper = (sample_upper_margin, labels[lam], k)

            active_mono_rows += 1
            mono_gap = arb_x[(lam, k + 1)] - arb_x[(lam, k)]
            sample_mono_gap = sample_x[(lam, k + 1)] - sample_x[(lam, k)]
            if arb_positive(mono_gap):
                active_mono_pos += 1
            if min_active_mono is None or sample_mono_gap < min_active_mono[0]:
                min_active_mono = (sample_mono_gap, labels[lam], k)

        for k in (first_collar,):
            lower = flint.arb(2 * k - 1) / flint.arb(2 * k + 1)
            sample_lower = Decimal(2 * k - 1) / Decimal(2 * k + 1)
            lower_margin = arb_x[(lam, k)] - lower
            sample_lower_margin = sample_x[(lam, k)] - sample_lower
            collar_lu_rows += 1
            if arb_positive(lower_margin):
                collar_lu_pos += 1
            if min_collar_lu is None or sample_lower_margin < min_collar_lu[0]:
                min_collar_lu = (sample_lower_margin, labels[lam], k)

            upper_margin = flint.arb(1) - arb_x[(lam, k)]
            sample_upper_margin = Decimal(1) - sample_x[(lam, k)]
            collar_lu_rows += 1
            if arb_positive(upper_margin):
                collar_lu_pos += 1
            if min_collar_lu is None or sample_upper_margin < min_collar_lu[0]:
                min_collar_lu = (sample_upper_margin, labels[lam], k)

        collar_gap = arb_x[(lam, second_collar)] - arb_x[(lam, first_collar)]
        sample_collar_gap = sample_x[(lam, second_collar)] - sample_x[(lam, first_collar)]
        collar_mono_rows += 1
        if arb_positive(collar_gap):
            collar_mono_pos += 1
        if min_collar_mono is None or sample_collar_gap < min_collar_mono[0]:
            min_collar_mono = (sample_collar_gap, labels[lam], first_collar)

    required_mins = (min_active_lower, min_active_upper, min_active_mono, min_collar_lu, min_collar_mono)
    if any(item is None for item in required_mins):
        raise RuntimeError("missing finite-collar diagnostics")

    return CollarDiagnostics(
        lambdas=[labels[lam] for lam in lambdas],
        coefficient_k_max=coefficient_k_max,
        available_x_max=available_x_max,
        lower_upper_k_max=lower_upper_k_max,
        monotone_gap_k_max=monotone_gap_k_max,
        certified_active_depth=certified_active_depth,
        first_collar_x=first_collar,
        second_collar_x=second_collar,
        active_lower_rows=active_lower_rows,
        active_lower_positive_rows=active_lower_pos,
        active_upper_rows=active_upper_rows,
        active_upper_positive_rows=active_upper_pos,
        active_monotone_rows=active_mono_rows,
        active_monotone_positive_rows=active_mono_pos,
        collar_lower_upper_rows=collar_lu_rows,
        collar_lower_upper_positive_rows=collar_lu_pos,
        collar_monotone_rows=collar_mono_rows,
        collar_monotone_positive_rows=collar_mono_pos,
        min_active_lower_margin=CollarMinimum(decimal_format(min_active_lower[0]), min_active_lower[1], min_active_lower[2]),
        min_active_upper_margin=CollarMinimum(decimal_format(min_active_upper[0]), min_active_upper[1], min_active_upper[2]),
        min_active_monotone_gap=CollarMinimum(decimal_format(min_active_mono[0]), min_active_mono[1], min_active_mono[2]),
        min_collar_lower_upper_margin=CollarMinimum(decimal_format(min_collar_lu[0]), min_collar_lu[1], min_collar_lu[2]),
        min_collar_monotone_gap=CollarMinimum(decimal_format(min_collar_mono[0]), min_collar_mono[1], min_collar_mono[2]),
    )


def build_contract(paths: list[Path], prefix_json: Path) -> dict:
    diagnostics = build_diagnostics(paths)
    return {
        "kind": "jensen_window_pf_negative_lambda_finite_collar_contract",
        "date": "2026-07-06",
        "status": "finite_collar_contract",
        "source_prefix_scout": "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.md",
        "source_ratio_cone_lemma": "outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md",
        "source_cone_entry_target": "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "prefix_json": prefix_json.relative_to(REPO_ROOT).as_posix(),
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_finite_collar_contract.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_finite_collar_contract.py",
        "proof_boundary": (
            "Finite-collar accounting contract only. It extracts the active "
            "finite depth currently supported by the negative-lambda prefix "
            "and the ratio-cone collar requirements, but it does not prove an "
            "infinite cone-entry theorem, does not prove a tail theorem, does "
            "not prove a collared finite flow theorem beyond the stated depth, "
            "does not prove jwpf_06, and leaves Lambda <= 0 unsettled."
        ),
        "contract_rows": [
            {
                "id": "nlfcc_01_ratio_cone_collar_rule",
                "role": "exact_requirement_extraction",
                "claim": "Active finite depth K needs lower/upper control through x_{K+1} and monotone-gap control through x_{K+2}.",
                "proof_boundary": "Requirement extraction only; not a cone-entry theorem.",
            },
            {
                "id": "nlfcc_02_current_active_depth",
                "role": "finite_prefix_accounting",
                "claim": f"The current negative-lambda prefix supports active depth K={diagnostics.certified_active_depth}.",
                "proof_boundary": "Finite depth only; not all k.",
            },
            {
                "id": "nlfcc_03_collar_certificate",
                "role": "finite_collar_certificate",
                "claim": (
                    f"For K={diagnostics.certified_active_depth}, the first collar x_{diagnostics.first_collar_x} "
                    f"and second collar x_{diagnostics.second_collar_x} are available on the checked lambdas."
                ),
                "proof_boundary": "Finite lambda grid and finite depth only.",
            },
            {
                "id": "nlfcc_04_next_upgrade",
                "role": "open_upgrade",
                "claim": "Extending active depth requires either more ACB coefficients and collars or an analytic all-k tail theorem.",
                "proof_boundary": "Upgrade target only.",
            },
        ],
        "finite_diagnostics": asdict(diagnostics),
        "summary": {
            "certified_active_depth": diagnostics.certified_active_depth,
            "first_collar_x": diagnostics.first_collar_x,
            "second_collar_x": diagnostics.second_collar_x,
            "active_lower_rows": diagnostics.active_lower_rows,
            "active_upper_rows": diagnostics.active_upper_rows,
            "active_monotone_rows": diagnostics.active_monotone_rows,
            "collar_lower_upper_rows": diagnostics.collar_lower_upper_rows,
            "collar_monotone_rows": diagnostics.collar_monotone_rows,
            "ready_to_apply_rows": 0,
            "target_closing": False,
            "main_finding": (
                "The current negative-lambda prefix is collar-ready only for "
                f"active depth K={diagnostics.certified_active_depth}, with "
                f"x_{diagnostics.first_collar_x} and x_{diagnostics.second_collar_x} "
                "acting as collar variables. "
                "This makes the finite evidence usable as a precise local seed, "
                "but it still does not provide the all-k tail theorem or the "
                "finite-collar flow theorem needed for the cone-entry target."
            ),
        },
        "invariants": [
            "The contract uses the finite-collar requirement from the exact ratio-cone invariance lemma.",
            "The certified active depth is finite and tied to available ACB enclosures.",
            "No row is ready_to_apply for the cone-entry theorem.",
            "The all-k tail or finite-collar flow theorem remains open.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(contract: dict, path: Path) -> None:
    diagnostics = contract["finite_diagnostics"]
    summary = contract["summary"]
    result_json = contract.get("result_json", "work/rh_compute/results/jensen_window_pf_negative_lambda_finite_collar_contract.json")
    note_ref = contract.get("note", path.relative_to(REPO_ROOT).as_posix())
    text = f"""# Jensen-Window PF Negative-Lambda Finite-Collar Contract

Date: 2026-07-06

Status: finite-collar accounting contract. This is not a proof of zeta cone
entry, monotone contractions, Jensen-window PF-infinity, RH, or the
Newman-direction goal.

Artifact kind: `jensen_window_pf_negative_lambda_finite_collar_contract`.

Proof boundary: this artifact extracts the finite active depth currently
supported by the negative-lambda prefix and the exact ratio-cone collar
requirements. It does not prove an infinite cone-entry theorem, does not prove
an all-`k` tail theorem, does not prove a finite-collar flow theorem beyond the
stated depth, does not prove `jwpf_06`, and does not settle `Lambda <= 0`.

Machine-readable contract:

```text
{result_json}
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_finite_collar_contract.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_finite_collar_contract.py --target {result_json} --note {note_ref}
```

Current result:

```text
validated Jensen-window PF negative-lambda finite-collar contract: active depth K={summary['certified_active_depth']}, {summary['active_lower_rows']} active lower rows, {summary['active_upper_rows']} active upper rows, {summary['active_monotone_rows']} active monotone rows, {summary['collar_lower_upper_rows']} collar lower/upper rows, {summary['collar_monotone_rows']} collar monotone rows, 0 issues
```

## Collar Rule

The exact ratio-cone invariance lemma states:

```text
checking x_1..x_K requires collar variables x_{{K+1}} for lower/upper walls
and x_{{K+2}} for monotone wall tests
```

For the current negative-lambda prefix:

```text
lambdas: {', '.join(diagnostics['lambdas'])}
available coefficient range: A_0..A_{diagnostics['coefficient_k_max']}
available ratio range: x_1..x_{diagnostics['available_x_max']}
lower/upper range: x_1..x_{diagnostics['lower_upper_k_max']}
monotone-gap range: 1 <= k <= {diagnostics['monotone_gap_k_max']}

certified active depth: K={diagnostics['certified_active_depth']}
first collar: x_{diagnostics['first_collar_x']}
second collar: x_{diagnostics['second_collar_x']}
```

## Certified Rows

```text
active lower rows: {diagnostics['active_lower_positive_rows']} / {diagnostics['active_lower_rows']}
active upper rows: {diagnostics['active_upper_positive_rows']} / {diagnostics['active_upper_rows']}
active monotone rows: {diagnostics['active_monotone_positive_rows']} / {diagnostics['active_monotone_rows']}
collar lower/upper rows: {diagnostics['collar_lower_upper_positive_rows']} / {diagnostics['collar_lower_upper_rows']}
collar monotone rows: {diagnostics['collar_monotone_positive_rows']} / {diagnostics['collar_monotone_rows']}

minimum active lower margin:
  {diagnostics['min_active_lower_margin']['sample']} at lambda={diagnostics['min_active_lower_margin']['lam']}, k={diagnostics['min_active_lower_margin']['k']}

minimum active upper margin:
  {diagnostics['min_active_upper_margin']['sample']} at lambda={diagnostics['min_active_upper_margin']['lam']}, k={diagnostics['min_active_upper_margin']['k']}

minimum active monotone gap:
  {diagnostics['min_active_monotone_gap']['sample']} at lambda={diagnostics['min_active_monotone_gap']['lam']}, k={diagnostics['min_active_monotone_gap']['k']}

minimum collar lower/upper margin:
  {diagnostics['min_collar_lower_upper_margin']['sample']} at lambda={diagnostics['min_collar_lower_upper_margin']['lam']}, k={diagnostics['min_collar_lower_upper_margin']['k']}

minimum collar monotone gap:
  {diagnostics['min_collar_monotone_gap']['sample']} at lambda={diagnostics['min_collar_monotone_gap']['lam']}, k={diagnostics['min_collar_monotone_gap']['k']}
```

## Meaning

The negative-lambda prefix is now precisely collar-accounted: it can seed the
active finite depth `K={diagnostics['certified_active_depth']}` with collars
`x_{diagnostics['first_collar_x']}` and `x_{diagnostics['second_collar_x']}`.
It does not seed `K={diagnostics['certified_active_depth'] + 1}` because the
monotone-wall collar would require `x_{diagnostics['second_collar_x'] + 1}`,
hence additional coefficient enclosures through
`A_{diagnostics['second_collar_x'] + 2}`, or an analytic tail theorem.

Related artifacts:

```text
outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.md
outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
```
"""
    path.write_text(text, encoding="utf-8")


def write_json(contract: dict, path: Path) -> None:
    path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--enclosure-jsonl", type=Path, nargs="+", default=list(DEFAULT_ENCLOSURE_JSONL))
    parser.add_argument("--prefix-json", type=Path, default=DEFAULT_PREFIX_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    paths = [path if path.is_absolute() else REPO_ROOT / path for path in args.enclosure_jsonl]
    prefix_json = args.prefix_json if args.prefix_json.is_absolute() else REPO_ROOT / args.prefix_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    contract = build_contract(paths, prefix_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    contract["result_json"] = out_json.relative_to(REPO_ROOT).as_posix()
    contract["note"] = note.relative_to(REPO_ROOT).as_posix()
    write_json(contract, out_json)
    write_note(contract, note)
    print(
        "wrote Jensen-window PF negative-lambda finite-collar contract: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
