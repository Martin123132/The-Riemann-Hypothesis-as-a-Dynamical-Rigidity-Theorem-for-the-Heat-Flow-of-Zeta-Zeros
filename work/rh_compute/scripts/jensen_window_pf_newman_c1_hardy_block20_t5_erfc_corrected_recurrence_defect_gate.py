#!/usr/bin/env python3
"""Propagate the t5 intrinsic-erfc correction through block-20 recurrences."""

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
REAL_ERFC = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_t5_real_erfc_residual_gate.json"
PRIOR_CORRECTED = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_corrected_recurrence_defect_gate.json"
ATLAS = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_discrete_selector_atlas_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_t5_erfc_corrected_recurrence_defect_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_t5_erfc_corrected_recurrence_defect_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_t5_erfc_corrected_recurrence_defect_gate.py"
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
    require(len(recursive) == 374, "t5-erfc corrected-defect roster drift")
    return recursive


def corrected_step(step: dict[str, Any], child: acb, q_delta: acb) -> acb:
    value = point_balls.binary128_complex(step["multiplier_hex"]) * child
    value += point_balls.binary128_complex(step["qq_hex"]) + q_delta
    if step["conjugate"]:
        value = value.conjugate()
    if step["subtract_one"]:
        value -= 1
    return value


def evaluate(
    data: dict[str, Any],
    special_record: dict[str, Any],
    erfc_record: dict[str, Any],
    dps: int,
) -> dict[str, acb]:
    ctx.dps = dps
    ctx.threads = 1
    special_delta = complex_from_record(special_record)
    erfc_delta = complex_from_record(erfc_record)
    header = data["header"]
    levels = data["levels"]
    steps = data["steps"]
    require(set(levels) == {1, 2} and set(steps) == {1}, "recursive chain shape drift")
    tpm = point_balls.binary128_ball(header["tpm_hex"])
    parent = point_balls.adapt(levels[1], point_balls.direct_sum(levels[1], tpm))
    child = point_balls.adapt(levels[2], point_balls.direct_sum(levels[2], tpm))
    source_model = point_balls.apply_step(steps[1], child)
    special_model = corrected_step(steps[1], child, special_delta)
    full_model = corrected_step(steps[1], child, special_delta + erfc_delta)
    source_defect = parent - source_model
    special_defect = parent - special_model
    full_defect = parent - full_model
    return {
        "source_defect": source_defect,
        "special_defect": special_defect,
        "full_defect": full_defect,
        "special_change": special_defect - source_defect,
        "erfc_increment": full_defect - special_defect,
        "full_change": full_defect - source_defect,
        "erfc_q_delta": erfc_delta,
    }


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    classes = aggregate["intrinsic_erfc_change_classification"]
    return f"""# Hardy block-20 t5-erfc-corrected recurrence-defect gate

Date: 2026-08-06

Status: rigorous finite recurrence propagation of the admitted PSI/complex-ERF and intrinsic real-erfc replacements; not a proof or saddle remainder theorem

## Comparison

For each of the 374 recursive calls, Arb independently resums the exact saved
length-104 parent and length-one-or-two child.  The three models are

```text
D_source  = parent - A(child, q_source),
D_special = parent - A(child, q_source + delta_q_PSI/complex-ERF),
D_full    = parent - A(child, q_source + delta_q_PSI/complex-ERF
                                  + delta_q_intrinsic-real-erfc).
```

`A` includes the exact saved multiplier, conjugation, and subtract-one branch.
All imported correction balls use their stored dyadic endpoints.  Evaluations
at 180 and 260 decimal digits overlap, and every recomputed `D_special`
overlaps the previously admitted corrected-defect certificate.

## Result

```text
maximum intrinsic-real-erfc defect increment <= {aggregate['maximum_intrinsic_erfc_defect_increment_abs_upper']}
maximum increment / prior-defect ratio        <= {aggregate['maximum_intrinsic_erfc_relative_change_upper']}
minimum surviving full-defect magnitude       >= {aggregate['minimum_full_corrected_defect_magnitude_lower']}
maximum surviving full-defect magnitude       <= {aggregate['maximum_full_corrected_defect_magnitude_upper']}
minimum defect / max erfc-increment separation >= {aggregate['minimum_global_scale_separation_lower']}
rigorously improved by real-erfc replacement    = {classes['improved']}
rigorously worsened by real-erfc replacement    = {classes['worsened']}
interval-indeterminate magnitude comparisons    = {classes['indeterminate']}
full corrected defects excluding zero           = {aggregate['full_corrected_defects_excluding_zero']} / 374
```

Thus the source intrinsic real `erfc` accuracy is not the observed recurrence
wall: its largest possible propagated change is over thirty orders of magnitude
below even the smallest surviving defect.  The surviving finite discrepancy is
now quantitatively assigned to the still-unproved endpoint-saddle, finite-`ip`,
Euler--Maclaurin, and remaining `t5/q` approximation contract (plus separately
unbounded arithmetic outside the admitted replacements).

This finite diagnostic does not establish a formula or uniform bound for that
remainder, transport it in height, control block or outer Hardy accumulation,
establish `Lambda<=0`, PF-infinity, RH, or a prize-level theorem.
"""


def main() -> int:
    for path in (TELEMETRY, SPECIAL, REAL_ERFC, PRIOR_CORRECTED, ATLAS, CHECKER):
        require(path.is_file(), f"missing t5-erfc corrected-defect dependency: {path}")
    special = json.loads(SPECIAL.read_text(encoding="utf-8"))
    real_erfc = json.loads(REAL_ERFC.read_text(encoding="utf-8"))
    prior = json.loads(PRIOR_CORRECTED.read_text(encoding="utf-8"))
    special_records = {
        row["chain"]: row["correlated_replacement"]["total_q_delta"]
        for row in special["rows"]
    }
    erfc_records = {
        row["chain"]: row["correlated_replacement"]["total_q_delta"]
        for row in real_erfc["rows"]
    }
    prior_records = {row["chain"]: row["corrected_defect"] for row in prior["rows"]}
    atlas = json.loads(ATLAS.read_text(encoding="utf-8"))
    expected = {row["chain"] for row in atlas["rows"] if row["mit"] == 2}
    recursive = load_recursive_chains()
    require(
        set(recursive) == set(special_records) == set(erfc_records) == set(prior_records) == expected,
        "t5-erfc corrected-defect roster mismatch",
    )

    rows: list[dict[str, Any]] = []
    classification: Counter[str] = Counter()
    maxima: dict[str, tuple[Fraction, int]] = {}
    minima: dict[str, tuple[Fraction, int]] = {}
    prior_overlap_count = 0

    def observe_max(name: str, value: Fraction, chain: int) -> None:
        if name not in maxima or value > maxima[name][0]:
            maxima[name] = (value, chain)

    def observe_min(name: str, value: Fraction, chain: int) -> None:
        if name not in minima or value < minima[name][0]:
            minima[name] = (value, chain)

    for chain in sorted(recursive):
        low = evaluate(recursive[chain], special_records[chain], erfc_records[chain], PRECISIONS[0])
        high = evaluate(recursive[chain], special_records[chain], erfc_records[chain], PRECISIONS[1])
        for name in (
            "source_defect",
            "special_defect",
            "full_defect",
            "special_change",
            "erfc_increment",
            "full_change",
        ):
            require(low[name].overlaps(high[name]), f"chain {chain} {name} precision nonoverlap")

        prior_ball = complex_from_record(prior_records[chain])
        require(high["special_defect"].overlaps(prior_ball), f"chain {chain} prior special-defect nonoverlap")
        prior_overlap_count += 1

        special_abs = abs(high["special_defect"])
        full_abs = abs(high["full_defect"])
        special_lower = cells.bound_fraction(special_abs.lower())
        special_upper = cells.bound_fraction(special_abs.upper())
        full_lower = cells.bound_fraction(full_abs.lower())
        full_upper = cells.bound_fraction(full_abs.upper())
        erfc_increment_upper = magnitude_upper(high["erfc_increment"])
        erfc_q_upper = magnitude_upper(high["erfc_q_delta"])
        relative_upper = erfc_increment_upper / special_lower
        if full_upper < special_lower:
            change_class = "improved"
        elif full_lower > special_upper:
            change_class = "worsened"
        else:
            change_class = "indeterminate"
        classification[change_class] += 1

        observe_min("special", special_lower, chain)
        observe_min("full", full_lower, chain)
        observe_max("special", special_upper, chain)
        observe_max("full", full_upper, chain)
        observe_max("erfc_increment", erfc_increment_upper, chain)
        observe_max("erfc_q", erfc_q_upper, chain)
        observe_max("relative", relative_upper, chain)

        header = recursive[chain]["header"]
        rows.append(
            {
                "chain": chain,
                "sum_index": int(header["sum_index"]),
                "branch": int(header["branch"]),
                "child_length": int(header["kernel_length"]),
                "special_corrected_defect": complex_record(high["special_defect"]),
                "special_corrected_defect_magnitude_lower": fraction_decimal(special_lower),
                "special_corrected_defect_magnitude_upper": fraction_decimal(special_upper),
                "intrinsic_real_erfc_q_delta": complex_record(high["erfc_q_delta"]),
                "intrinsic_real_erfc_q_delta_abs_upper": fraction_decimal(erfc_q_upper),
                "intrinsic_real_erfc_defect_increment": complex_record(high["erfc_increment"]),
                "intrinsic_real_erfc_defect_increment_abs_upper": fraction_decimal(erfc_increment_upper),
                "full_corrected_defect": complex_record(high["full_defect"]),
                "full_corrected_defect_magnitude_lower": fraction_decimal(full_lower),
                "full_corrected_defect_magnitude_upper": fraction_decimal(full_upper),
                "intrinsic_real_erfc_relative_change_upper": fraction_decimal(relative_upper),
                "change_classification": change_class,
                "full_corrected_defect_excludes_zero": not high["full_defect"].contains(0),
                "precision_overlap": True,
                "prior_special_corrected_defect_overlap": True,
            }
        )

    exclude_zero = sum(row["full_corrected_defect_excludes_zero"] for row in rows)
    global_separation = minima["full"][0] / maxima["erfc_increment"][0]
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_t5_erfc_corrected_recurrence_defect_gate",
        "status": "rigorous_374_call_special_and_intrinsic_erfc_corrected_recurrence_defects_enclosed",
        "scope": {
            "block": 20,
            "recursive_chain_count": 374,
            "intrinsic_real_erfc_evaluation_count": 913,
            "precision_ladder_decimal_digits": list(PRECISIONS),
            "precision_overlap_count": 374,
            "prior_special_corrected_defect_overlap_count": prior_overlap_count,
            "parent_upper_index": 104,
            "child_upper_indices": [1, 2],
        },
        "correction_contract": {
            "source_defect": "adapted exact parent - adapted(multiplier*exact child + source q)",
            "special_corrected_defect": "source model plus rigorous PSI/complex-ERF q replacement",
            "full_corrected_defect": "special-corrected model plus rigorous intrinsic-real-erfc q replacement",
            "delta_sources": "stored dyadic endpoints from the admitted special-function and t5 real-erfc residual gates",
            "affine_magnitude_preservation": "Conjugation and subtract-one preserve the magnitude of the injected q change.",
        },
        "aggregate": {
            "minimum_special_corrected_defect_magnitude_lower": fraction_decimal(minima["special"][0]),
            "minimum_special_corrected_defect_witness": minima["special"][1],
            "minimum_full_corrected_defect_magnitude_lower": fraction_decimal(minima["full"][0]),
            "minimum_full_corrected_defect_witness": minima["full"][1],
            "maximum_special_corrected_defect_magnitude_upper": fraction_decimal(maxima["special"][0]),
            "maximum_special_corrected_defect_witness": maxima["special"][1],
            "maximum_full_corrected_defect_magnitude_upper": fraction_decimal(maxima["full"][0]),
            "maximum_full_corrected_defect_witness": maxima["full"][1],
            "maximum_intrinsic_erfc_q_delta_abs_upper": fraction_decimal(maxima["erfc_q"][0]),
            "maximum_intrinsic_erfc_q_delta_witness": maxima["erfc_q"][1],
            "maximum_intrinsic_erfc_defect_increment_abs_upper": fraction_decimal(maxima["erfc_increment"][0]),
            "maximum_intrinsic_erfc_defect_increment_witness": maxima["erfc_increment"][1],
            "maximum_intrinsic_erfc_relative_change_upper": fraction_decimal(maxima["relative"][0]),
            "maximum_intrinsic_erfc_relative_change_witness": maxima["relative"][1],
            "minimum_global_scale_separation_lower": fraction_decimal(global_separation),
            "intrinsic_erfc_change_classification": {
                name: classification[name] for name in ("improved", "worsened", "indeterminate")
            },
            "full_corrected_defects_excluding_zero": exclude_zero,
        },
        "rows": rows,
        "sources": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "special_residual": {"path": relative(SPECIAL), "sha256": file_hash(SPECIAL)},
            "t5_real_erfc_residual": {"path": relative(REAL_ERFC), "sha256": file_hash(REAL_ERFC)},
            "prior_corrected_defect": {"path": relative(PRIOR_CORRECTED), "sha256": file_hash(PRIOR_CORRECTED)},
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
            "analytic_target": "Derive a source-matched remainder formula for the endpoint-saddle, finite-ip, and Euler-Maclaurin terms remaining in t5/q.",
            "finite_target": "Test candidate remainder laws against all 374 surviving full-corrected complex defects without fitting and certifying on the same roster.",
            "uniform_target": "Prove height- and selector-uniform bounds before any extrapolation beyond this finite block."
        },
        "proof_boundary": (
            "Rigorous finite recurrence propagation of admitted PSI/complex-ERF and intrinsic-real-erfc "
            "replacements for 374 block-20 calls only. The surviving nonzero defects identify but do not "
            "derive or bound the endpoint-saddle, finite-ip, Euler-Maclaurin, or remaining q remainder. "
            "This does not control outer accumulation, establish height uniformity, Lambda <= 0, "
            "PF-infinity, RH, or a prize-level conclusion."
        ),
    }
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(
        "built Hardy block-20 t5-erfc-corrected recurrence-defect gate: "
        f"374 defects, {exclude_zero} nonzero, scale separation {artifact['aggregate']['minimum_global_scale_separation_lower']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
