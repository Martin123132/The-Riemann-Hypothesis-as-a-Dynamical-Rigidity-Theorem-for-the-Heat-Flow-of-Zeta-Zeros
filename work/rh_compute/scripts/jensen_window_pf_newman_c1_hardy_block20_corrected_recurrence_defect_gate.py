#!/usr/bin/env python3
"""Enclose all recursive block-20 defects after the special-function correction."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import acb, arb, ctx


TELEMETRY = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10_block20_full/enabled/chain_telemetry.jsonl"
)
SPECIAL = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_special_function_residual_gate.json"
ATLAS = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_discrete_selector_atlas_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_corrected_recurrence_defect_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_corrected_recurrence_defect_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_corrected_recurrence_defect_gate.py"
PRECISIONS = (180, 260)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def fraction_decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def dyadic_fraction(payload: list[int]) -> Fraction:
    require(len(payload) == 2, "dyadic payload drift")
    mantissa, exponent = int(payload[0]), int(payload[1])
    return Fraction(mantissa * 2**exponent, 1) if exponent >= 0 else Fraction(mantissa, 2 ** (-exponent))


def interval_from_record(record: dict[str, Any]) -> arb:
    lower = dyadic_fraction(record["lower_dyadic"])
    upper = dyadic_fraction(record["upper_dyadic"])
    require(lower <= upper, "stored interval order drift")
    center = (lower + upper) / 2
    radius = (upper - lower) / 2
    return cells.fraction_ball(center) if radius == 0 else cells.interval_ball(center, radius)


def complex_from_record(record: dict[str, Any]) -> acb:
    return acb(interval_from_record(record["real"]), interval_from_record(record["imag"]))


def magnitude_upper(value: arb | acb) -> Fraction:
    return cells.bound_fraction(abs(value).upper())


def magnitude_lower(value: acb) -> Fraction:
    return cells.bound_fraction(abs(value).lower())


def complex_record(value: acb) -> dict[str, Any]:
    return point_balls.acb_record(value, 55)


def load_recursive_chains() -> dict[int, dict[str, Any]]:
    chains: dict[int, dict[str, Any]] = {}
    for line in TELEMETRY.read_text(encoding="utf-8").splitlines():
        record = json.loads(line)
        chain = int(record.get("chain", 0))
        if record["type"] == "chain":
            chains[chain] = {"header": record, "levels": {}, "steps": {}}
        elif record["type"] == "level":
            chains[chain]["levels"][int(record["level"])] = record
        elif record["type"] == "recurrence":
            chains[chain]["steps"][int(record["nit"])] = record
    recursive = {chain: data for chain, data in chains.items() if int(data["header"]["mit"]) == 2}
    require(len(recursive) == 374, "corrected-defect roster drift")
    return recursive


def corrected_step(step: dict[str, Any], child: acb, q_delta: acb) -> acb:
    value = point_balls.binary128_complex(step["multiplier_hex"]) * child
    value += point_balls.binary128_complex(step["qq_hex"]) + q_delta
    if step["conjugate"]:
        value = value.conjugate()
    if step["subtract_one"]:
        value -= 1
    return value


def evaluate(data: dict[str, Any], q_delta: acb, dps: int) -> dict[str, acb]:
    ctx.dps = dps
    ctx.threads = 1
    header = data["header"]
    levels = data["levels"]
    steps = data["steps"]
    require(set(levels) == {1, 2} and set(steps) == {1}, "recursive chain shape drift")
    tpm = point_balls.binary128_ball(header["tpm_hex"])
    parent = point_balls.adapt(levels[1], point_balls.direct_sum(levels[1], tpm))
    child = point_balls.adapt(levels[2], point_balls.direct_sum(levels[2], tpm))
    source_model = point_balls.apply_step(steps[1], child)
    corrected_model = corrected_step(steps[1], child, q_delta)
    source_defect = parent - source_model
    corrected_defect = parent - corrected_model
    logged_defect = point_balls.binary128_complex(steps[1]["local_defect_hex"])
    return {
        "source_defect": source_defect,
        "corrected_defect": corrected_defect,
        "logged_defect": logged_defect,
        "source_roundoff_gap": source_defect - logged_defect,
        "defect_change": corrected_defect - source_defect,
    }


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    classes = aggregate["rigorous_change_classification"]
    return f"""# Hardy block-20 corrected recurrence-defect gate

Date: 2026-08-06

Status: rigorous finite exact-point recurrence defects after the PSI/ERF replacement; not a proof or uniform recurrence theorem

## Comparison

For each of the 374 `MIT=2` calls, Arb independently sums the exact saved
length-104 parent and length-one-or-two child kernels.  It then compares

```text
D_source = parent - A(source child, source q),
D_corr   = parent - A(source child, source q + delta q_special),
```

where `A` includes the exact saved multiplier, conjugation, and subtract-one
branches.  The correction `delta q_special` is imported from the dyadic Arb
endpoints of Section 11.249, so no decimal re-rounding is used.  Evaluations at
180 and 260 decimal digits overlap for every source and corrected defect.

## Result

```text
maximum independent source/logged defect gap <= {aggregate['maximum_source_roundoff_gap_abs_upper']}
minimum source defect magnitude              >= {aggregate['minimum_source_defect_magnitude_lower']}
minimum corrected defect magnitude           >= {aggregate['minimum_corrected_defect_magnitude_lower']}
maximum source defect magnitude              <= {aggregate['maximum_source_defect_magnitude_upper']}
maximum corrected defect magnitude           <= {aggregate['maximum_corrected_defect_magnitude_upper']}
maximum special-function defect change       <= {aggregate['maximum_defect_change_abs_upper']}
maximum relative change budget               <= {aggregate['maximum_relative_change_upper']}
rigorously improved calls                      = {classes['improved']}
rigorously worsened calls                      = {classes['worsened']}
interval-indeterminate comparisons             = {classes['indeterminate']}
corrected defects excluding zero               = {aggregate['corrected_defects_excluding_zero']} / 374
```

The special-function replacement is therefore propagated through the actual
one-step recurrence, but it does not remove the finite discrepancy.  The
remaining nonzero defect is an explicit target for the intrinsic real `erfc`,
endpoint-saddle, Euler--Maclaurin, and other mathematical `t5/q` approximation
terms rather than an untracked numerical effect.

This is a finite low-height diagnostic and error certificate.  It does not
derive a sign or uniform bound for the remaining defect, transport it to large
heights, combine the 424 calls with an outer Hardy remainder, prove
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
"""


def main() -> int:
    for path in (TELEMETRY, SPECIAL, ATLAS, CHECKER):
        require(path.is_file(), f"missing corrected-defect dependency: {path}")
    special = json.loads(SPECIAL.read_text(encoding="utf-8"))
    delta_rows = {
        row["chain"]: complex_from_record(row["correlated_replacement"]["total_q_delta"])
        for row in special["rows"]
    }
    atlas = json.loads(ATLAS.read_text(encoding="utf-8"))
    expected = {row["chain"] for row in atlas["rows"] if row["mit"] == 2}
    recursive = load_recursive_chains()
    require(set(recursive) == set(delta_rows) == expected, "corrected-defect roster mismatch")

    rows: list[dict[str, Any]] = []
    classification: Counter[str] = Counter()
    maxima: dict[str, tuple[Fraction, int]] = {}
    minima: dict[str, tuple[Fraction, int]] = {}

    def observe_max(name: str, value: Fraction, chain: int) -> None:
        if name not in maxima or value > maxima[name][0]:
            maxima[name] = (value, chain)

    def observe_min(name: str, value: Fraction, chain: int) -> None:
        if name not in minima or value < minima[name][0]:
            minima[name] = (value, chain)

    for chain in sorted(recursive):
        low = evaluate(recursive[chain], delta_rows[chain], PRECISIONS[0])
        high = evaluate(recursive[chain], delta_rows[chain], PRECISIONS[1])
        for name in ("source_defect", "corrected_defect", "source_roundoff_gap", "defect_change"):
            require(low[name].overlaps(high[name]), f"chain {chain} {name} precision nonoverlap")

        source_abs = abs(high["source_defect"])
        corrected_abs = abs(high["corrected_defect"])
        source_lower = cells.bound_fraction(source_abs.lower())
        source_upper = cells.bound_fraction(source_abs.upper())
        corrected_lower = cells.bound_fraction(corrected_abs.lower())
        corrected_upper = cells.bound_fraction(corrected_abs.upper())
        roundoff_upper = magnitude_upper(high["source_roundoff_gap"])
        change_upper = magnitude_upper(high["defect_change"])
        relative_upper = change_upper / source_lower
        if corrected_upper < source_lower:
            change_class = "improved"
        elif corrected_lower > source_upper:
            change_class = "worsened"
        else:
            change_class = "indeterminate"
        classification[change_class] += 1

        observe_min("source", source_lower, chain)
        observe_min("corrected", corrected_lower, chain)
        observe_max("source", source_upper, chain)
        observe_max("corrected", corrected_upper, chain)
        observe_max("roundoff", roundoff_upper, chain)
        observe_max("change", change_upper, chain)
        observe_max("relative", relative_upper, chain)

        header = recursive[chain]["header"]
        rows.append(
            {
                "chain": chain,
                "sum_index": int(header["sum_index"]),
                "branch": int(header["branch"]),
                "child_length": int(header["kernel_length"]),
                "source_defect": complex_record(high["source_defect"]),
                "source_defect_magnitude_lower": fraction_decimal(source_lower),
                "source_defect_magnitude_upper": fraction_decimal(source_upper),
                "corrected_defect": complex_record(high["corrected_defect"]),
                "corrected_defect_magnitude_lower": fraction_decimal(corrected_lower),
                "corrected_defect_magnitude_upper": fraction_decimal(corrected_upper),
                "source_roundoff_gap_abs_upper": fraction_decimal(roundoff_upper),
                "defect_change_abs_upper": fraction_decimal(change_upper),
                "relative_change_upper": fraction_decimal(relative_upper),
                "change_classification": change_class,
                "corrected_defect_excludes_zero": not high["corrected_defect"].contains(0),
                "precision_overlap": True,
            }
        )

    exclude_zero = sum(row["corrected_defect_excludes_zero"] for row in rows)
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_corrected_recurrence_defect_gate",
        "status": "rigorous_374_call_special_corrected_recurrence_defects_enclosed",
        "scope": {
            "block": 20,
            "recursive_chain_count": 374,
            "precision_ladder_decimal_digits": list(PRECISIONS),
            "precision_overlap_count": 374,
            "parent_upper_index": 104,
            "child_upper_indices": [1, 2],
        },
        "correction_contract": {
            "source_defect": "adapted exact parent - adapted(multiplier*exact child + source q)",
            "corrected_defect": "adapted exact parent - adapted(multiplier*exact child + source q + rigorous delta q special)",
            "delta_source": "dyadic endpoints from the admitted block-20 special-function residual gate",
            "affine_magnitude_preservation": "Conjugation and subtract-one preserve the magnitude of the injected q change.",
        },
        "aggregate": {
            "minimum_source_defect_magnitude_lower": fraction_decimal(minima["source"][0]),
            "minimum_source_defect_witness": minima["source"][1],
            "minimum_corrected_defect_magnitude_lower": fraction_decimal(minima["corrected"][0]),
            "minimum_corrected_defect_witness": minima["corrected"][1],
            "maximum_source_defect_magnitude_upper": fraction_decimal(maxima["source"][0]),
            "maximum_source_defect_witness": maxima["source"][1],
            "maximum_corrected_defect_magnitude_upper": fraction_decimal(maxima["corrected"][0]),
            "maximum_corrected_defect_witness": maxima["corrected"][1],
            "maximum_source_roundoff_gap_abs_upper": fraction_decimal(maxima["roundoff"][0]),
            "maximum_source_roundoff_gap_witness": maxima["roundoff"][1],
            "maximum_defect_change_abs_upper": fraction_decimal(maxima["change"][0]),
            "maximum_defect_change_witness": maxima["change"][1],
            "maximum_relative_change_upper": fraction_decimal(maxima["relative"][0]),
            "maximum_relative_change_witness": maxima["relative"][1],
            "rigorous_change_classification": {
                name: classification[name] for name in ("improved", "worsened", "indeterminate")
            },
            "corrected_defects_excluding_zero": exclude_zero,
        },
        "rows": rows,
        "sources": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "special_residual": {"path": relative(SPECIAL), "sha256": file_hash(SPECIAL)},
            "selector_atlas": {"path": relative(ATLAS), "sha256": file_hash(ATLAS)},
            "point_ball_builder": {
                "path": relative(Path(point_balls.__file__).resolve()),
                "sha256": file_hash(Path(point_balls.__file__).resolve()),
            },
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "python_flint": {"version": point_balls.flint.__version__, "flint_version": point_balls.flint.__FLINT_VERSION__},
        },
        "next_handoff": {
            "analytic_target": "Bound or replace the remaining intrinsic-erfc, endpoint-saddle, Euler-Maclaurin, and t5/q truncation terms that make up the corrected nonzero defect.",
            "finite_accumulation": "Once a justified local-error contract is chosen, accumulate all 374 corrected recursive defects with the 50 direct-kernel certificates through the block weights.",
            "uniform_target": "Derive height-dependent bounds before extrapolating this low-height finite roster."
        },
        "proof_boundary": (
            "Rigorous finite exact-point source and special-corrected recurrence defects for 374 block-20 "
            "calls only. Nonzero defects diagnose a remaining approximation obligation; they do not prove "
            "its sign or a uniform bound, control the outer Hardy representation, establish height "
            "uniformity, Lambda <= 0, PF-infinity, RH, or a prize-level conclusion."
        ),
    }
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(
        "built Hardy block-20 corrected recurrence-defect gate: "
        f"374 defects, {exclude_zero} corrected nonzero, "
        f"{classification['improved']} improved/{classification['worsened']} worsened"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
