#!/usr/bin/env python3
"""Enclose every intrinsic real-erfc residual inside block-20 t5."""

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

import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import acb, arb, ctx


SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
PROBE_RESULT = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/t5_component_probe/t1e10_block20/probe_result.json"
)
PROBE_OUTPUT = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/t5_component_probe/t1e10_block20/probe_output.txt"
)
ATLAS = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_discrete_selector_atlas_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_t5_real_erfc_residual_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_t5_real_erfc_residual_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_t5_real_erfc_residual_gate.py"
PRECISION_DPS = 220


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def exact_arb(payload: str) -> arb:
    return cells.fraction_ball(cells.binary128_fraction(payload))


def exact_acb(payload: list[str]) -> acb:
    require(len(payload) == 2, "complex binary128 payload drift")
    return acb(exact_arb(payload[0]), exact_arb(payload[1]))


def fraction_decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def interval_record(value: arb) -> dict[str, Any]:
    lower, upper = cells.bounds(value)
    return {
        "lower": fraction_decimal(lower),
        "upper": fraction_decimal(upper),
        "lower_dyadic": [lower.numerator, -(lower.denominator.bit_length() - 1)],
        "upper_dyadic": [upper.numerator, -(upper.denominator.bit_length() - 1)],
    }


def complex_interval_record(value: acb) -> dict[str, dict[str, Any]]:
    return {"real": interval_record(value.real), "imag": interval_record(value.imag)}


def magnitude_upper(value: arb | acb) -> Fraction:
    return cells.bound_fraction(abs(value).upper())


def argument_bucket(argument: Fraction) -> str:
    if argument < 0:
        return "negative"
    for upper, name in (
        (Fraction(1), "[0,1)"),
        (Fraction(2), "[1,2)"),
        (Fraction(4), "[2,4)"),
        (Fraction(8), "[4,8)"),
        (Fraction(16), "[8,16)"),
        (Fraction(32), "[16,32)"),
    ):
        if argument < upper:
            return name
    return "[32,infinity)"


def source_locations() -> dict[str, Any]:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    locations: dict[str, Any] = {}
    tokens = {
        "pi_definition": "p=4*ATAN(C)",
        "sqrt_pi_definition": "sp=sqrt(p)",
        "t5_start": "sum0=(0.0,0.0)",
        "t5_finish": "t5=t5+(SUM0+SUM1)",
        "q_combination": "qq=conjg(t1+t2+t4)+t3+t5",
    }
    for name, token in tokens.items():
        matches = [number for number, line in enumerate(lines, 1) if token in line]
        require(matches, f"source t5-erfc token missing: {token}")
        locations[name] = matches[0]
    erfc_matches = [number for number, line in enumerate(lines, 1) if "sav1=erfc(sav1)" in line]
    require(len(erfc_matches) == 2, "source t5 intrinsic-erfc call roster drift")
    locations["intrinsic_real_erfc_calls"] = erfc_matches
    return locations


def load_probe() -> tuple[dict[str, str], dict[int, dict[str, Any]]]:
    lines = [line.split() for line in PROBE_OUTPUT.read_text(encoding="ascii").splitlines() if line.strip()]
    require(lines and lines[0][0] == "#constants" and len(lines[0]) == 7, "t5 probe header drift")
    constants = {
        "pi": lines[0][1],
        "sqrt_pi": lines[0][2],
        "two_pi": lines[0][3],
        "minus_two_pi": lines[0][4],
        "epi4_real": lines[0][5],
        "epi4_imag": lines[0][6],
    }
    rows: dict[int, dict[str, Any]] = {}
    for fields in lines[1:]:
        require(len(fields) == 36, "t5 residual probe field count drift")
        chain = int(fields[0])
        calls: list[dict[str, Any]] = []
        for slot in range(3):
            offset = 9 + 9 * slot
            active = int(fields[offset])
            if not active:
                continue
            calls.append(
                {
                    "slot": slot + 1,
                    "region": "lower" if int(fields[offset + 1]) == 0 else "upper",
                    "index": int(fields[offset + 2]),
                    "argument_hex": fields[offset + 3],
                    "source_hex": fields[offset + 4],
                    "weight_hex": fields[offset + 5 : offset + 7],
                    "source_term_hex": fields[offset + 7 : offset + 9],
                }
            )
        rows[chain] = {
            "t5_psi_hex": fields[1:3],
            "lower_sum_hex": fields[3:5],
            "upper_sum_hex": fields[5:7],
            "t5_hex": fields[7:9],
            "calls": calls,
        }
    require(len(rows) == 374, "t5 residual probe row count drift")
    require(sum(len(row["calls"]) for row in rows.values()) == 913, "t5 residual call count drift")
    return constants, rows


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 t5 intrinsic real-erfc residual gate

Date: 2026-08-06

Status: rigorous finite exact-argument intrinsic-erfc residuals and correlated t5/q replacement; not a proof or saddle-error theorem

## Admitted source observation

The source-derived one-CPU probe reproduces all 374 saved binary128 `t5`
values bit-for-bit.  Its admitted call roster is

```text
374 recursive calls,
539 lower-endpoint intrinsic real-erfc evaluations,
374 upper-endpoint intrinsic real-erfc evaluations,
913 total intrinsic real-erfc evaluations.
```

For every call, this gate starts from the exact emitted binary128 argument `x`
and encloses `Arb.erfc(x) - source_erfc(x)`.  It does not reconstruct `x` from
decimal data.  In the source, pi enters the argument through `p=4*atan(1)` and
`sp=sqrt(p)`; the rigorous reference inserts no additional pi because the
already-generated binary128 argument is the exact comparison point.

## Linear correlated replacement

For each emitted source weight `W_k`, define

```text
delta_k = W_k * (Arb.erfc(x_k) - source_erfc(x_k)),
delta_q_real_erfc = sum_k delta_k.
```

The sum is taken jointly within each recursive call, preserving cancellation.
The source-term association gap `source_term - W_k*source_erfc(x_k)` is measured
separately and is not relabelled as intrinsic-erfc error.

## Result

```text
maximum intrinsic-erfc residual             <= {aggregate['maximum_erfc_residual_abs_upper']}
maximum emitted weight magnitude            <= {aggregate['maximum_weight_abs_upper']}
maximum one-call weighted replacement       <= {aggregate['maximum_weighted_call_delta_abs_upper']}
maximum correlated t5/q replacement         <= {aggregate['maximum_total_q_delta_abs_upper']}
median correlated t5/q replacement          <= {aggregate['median_total_q_delta_abs_upper']}
maximum source-term association gap         <= {aggregate['maximum_source_term_association_gap_abs_upper']}
source outputs enclosed by rigorous erfc ball = {aggregate['source_output_enclosed_count']} / 913
source outputs exactly zero                   = {aggregate['source_output_zero_count']}
```

This isolates the accuracy of the intrinsic real `erfc` calls at the 913 finite
source arguments.  It does not validate the saddle locations, endpoint formula,
finite `ip` endpoint sums, omitted Euler--Maclaurin terms, other binary128
arithmetic, recurrence accumulation, outer Hardy remainder, height uniformity,
RH, or a prize-level theorem.
"""


def main() -> int:
    for path in (SOURCE, PROBE_RESULT, PROBE_OUTPUT, ATLAS, CHECKER):
        require(path.is_file(), f"missing t5-erfc residual dependency: {path}")
    probe_artifact = json.loads(PROBE_RESULT.read_text(encoding="utf-8"))
    validation = probe_artifact["validation"]
    require(validation["saved_t5_bit_match_count"] == 374, "t5 probe not admitted")
    require(validation["intrinsic_erfc_call_count"] == 913, "t5 probe call roster drift")
    atlas = json.loads(ATLAS.read_text(encoding="utf-8"))
    expected_chains = {row["chain"] for row in atlas["rows"] if row["mit"] == 2}
    constants, probe = load_probe()
    require(set(probe) == expected_chains, "t5-erfc recursive roster mismatch")

    ctx.dps = PRECISION_DPS
    ctx.threads = 1
    argument_histogram: Counter[str] = Counter()
    region_histogram: Counter[str] = Counter()
    maxima: dict[str, tuple[Fraction, dict[str, Any]]] = {}
    total_delta_uppers: list[Fraction] = []
    source_output_enclosed = 0
    source_output_zero = 0
    rows: list[dict[str, Any]] = []

    def observe(name: str, value: Fraction, witness: dict[str, Any]) -> None:
        if name not in maxima or value > maxima[name][0]:
            maxima[name] = (value, witness)

    for chain in sorted(probe):
        call_rows: list[dict[str, Any]] = []
        total_delta = acb(0)
        total_association_gap = acb(0)
        for call in probe[chain]["calls"]:
            argument_fraction = cells.binary128_fraction(call["argument_hex"])
            argument = cells.fraction_ball(argument_fraction)
            source_value = exact_arb(call["source_hex"])
            true_value = argument.erfc()
            residual = true_value - source_value
            residual_upper = magnitude_upper(residual)
            weight = exact_acb(call["weight_hex"])
            source_term = exact_acb(call["source_term_hex"])
            weighted_delta = weight * residual
            association_gap = source_term - weight * source_value
            weight_upper = magnitude_upper(weight)
            weighted_upper = magnitude_upper(weighted_delta)
            association_upper = magnitude_upper(association_gap)
            total_delta += weighted_delta
            total_association_gap += association_gap

            bucket = argument_bucket(argument_fraction)
            argument_histogram[bucket] += 1
            region_histogram[call["region"]] += 1
            source_output_enclosed += true_value.contains(source_value)
            source_output_zero += cells.binary128_fraction(call["source_hex"]) == 0
            witness = {
                "chain": chain,
                "slot": call["slot"],
                "region": call["region"],
                "index": call["index"],
            }
            observe("residual", residual_upper, witness)
            observe("weight", weight_upper, witness)
            observe("weighted", weighted_upper, witness)
            observe("association", association_upper, witness)
            call_rows.append(
                {
                    **witness,
                    "argument_hex": call["argument_hex"],
                    "argument_bucket": bucket,
                    "source_erfc_hex": call["source_hex"],
                    "source_weight_hex": call["weight_hex"],
                    "source_term_hex": call["source_term_hex"],
                    "residual": interval_record(residual),
                    "residual_abs_upper": fraction_decimal(residual_upper),
                    "source_output_enclosed_by_reference": true_value.contains(source_value),
                    "weighted_delta": complex_interval_record(weighted_delta),
                    "weighted_delta_abs_upper": fraction_decimal(weighted_upper),
                    "source_term_association_gap": complex_interval_record(association_gap),
                    "source_term_association_gap_abs_upper": fraction_decimal(association_upper),
                }
            )

        total_upper = magnitude_upper(total_delta)
        total_association_upper = magnitude_upper(total_association_gap)
        total_delta_uppers.append(total_upper)
        observe("total", total_upper, {"chain": chain})
        observe("total_association", total_association_upper, {"chain": chain})
        rows.append(
            {
                "chain": chain,
                "t5_psi_hex": probe[chain]["t5_psi_hex"],
                "lower_sum_hex": probe[chain]["lower_sum_hex"],
                "upper_sum_hex": probe[chain]["upper_sum_hex"],
                "saved_t5_hex": probe[chain]["t5_hex"],
                "calls": call_rows,
                "correlated_replacement": {
                    "total_q_delta": complex_interval_record(total_delta),
                    "total_q_delta_abs_upper": fraction_decimal(total_upper),
                    "total_source_term_association_gap": complex_interval_record(total_association_gap),
                    "total_source_term_association_gap_abs_upper": fraction_decimal(total_association_upper),
                },
            }
        )

    require(sum(argument_histogram.values()) == 913, "t5-erfc argument histogram drift")
    require(region_histogram == Counter({"lower": 539, "upper": 374}), "t5-erfc region histogram drift")
    ordered = sorted(total_delta_uppers)
    thresholds = {
        "greater_than_1e-25": Fraction(1, 10**25),
        "greater_than_1e-30": Fraction(1, 10**30),
        "greater_than_1e-35": Fraction(1, 10**35),
        "greater_than_1e-40": Fraction(1, 10**40),
    }
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_t5_real_erfc_residual_gate",
        "status": "rigorous_913_exact_argument_intrinsic_erfc_residuals_and_correlated_q_replacement_enclosed",
        "scope": {
            "block": 20,
            "recursive_chain_count": 374,
            "lower_erfc_evaluation_count": 539,
            "upper_erfc_evaluation_count": 374,
            "intrinsic_real_erfc_evaluation_count": 913,
            "saved_t5_bit_replay_count": 374,
            "precision_decimal_digits": PRECISION_DPS,
        },
        "comparison_contract": {
            "reference": "Arb erfc(x) at the exact emitted binary128 source argument x",
            "residual": "Arb.erfc(x)-exact_binary128_source_erfc(x)",
            "correlated_q_replacement": "sum_k exact_binary128_emitted_weight_k * residual_k within each source q call",
            "association_gap": "exact_binary128_source_term_k - exact_binary128_emitted_weight_k*exact_binary128_source_erfc_k",
            "separation": "The association gap is reported separately and is not included in the intrinsic-erfc residual correction.",
        },
        "pi_provenance": {
            "source_argument_generation": "p=4*atan(1), sp=sqrt(p), followed by the source t5 argument expressions",
            "rigorous_reference": "The admitted binary128 argument is used directly; Arb introduces no independent pi.",
            "constants_hex": constants,
        },
        "aggregate": {
            "argument_bucket_histogram": dict(argument_histogram),
            "region_histogram": dict(region_histogram),
            "source_output_enclosed_count": source_output_enclosed,
            "source_output_zero_count": source_output_zero,
            "maximum_erfc_residual_abs_upper": fraction_decimal(maxima["residual"][0]),
            "maximum_erfc_residual_witness": maxima["residual"][1],
            "maximum_weight_abs_upper": fraction_decimal(maxima["weight"][0]),
            "maximum_weight_witness": maxima["weight"][1],
            "maximum_weighted_call_delta_abs_upper": fraction_decimal(maxima["weighted"][0]),
            "maximum_weighted_call_delta_witness": maxima["weighted"][1],
            "maximum_total_q_delta_abs_upper": fraction_decimal(maxima["total"][0]),
            "maximum_total_q_delta_witness": maxima["total"][1],
            "median_total_q_delta_abs_upper": fraction_decimal(ordered[len(ordered) // 2]),
            "p95_total_q_delta_abs_upper": fraction_decimal(ordered[(95 * (len(ordered) - 1)) // 100]),
            "q_delta_threshold_counts": {
                name: sum(value > threshold for value in total_delta_uppers)
                for name, threshold in thresholds.items()
            },
            "maximum_source_term_association_gap_abs_upper": fraction_decimal(maxima["association"][0]),
            "maximum_source_term_association_gap_witness": maxima["association"][1],
            "maximum_correlated_source_term_association_gap_abs_upper": fraction_decimal(maxima["total_association"][0]),
            "maximum_correlated_source_term_association_gap_witness": maxima["total_association"][1],
        },
        "rows": rows,
        "source": {"path": relative(SOURCE), "sha256": file_hash(SOURCE), "locations": source_locations()},
        "sources": {
            "probe_result": {"path": relative(PROBE_RESULT), "sha256": file_hash(PROBE_RESULT)},
            "probe_output": {"path": relative(PROBE_OUTPUT), "sha256": file_hash(PROBE_OUTPUT)},
            "selector_atlas": {"path": relative(ATLAS), "sha256": file_hash(ATLAS)},
            "local_cell_builder": {
                "path": relative(Path(cells.__file__).resolve()),
                "sha256": file_hash(Path(cells.__file__).resolve()),
            },
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "python_flint": {"version": cells.flint.__version__, "flint_version": cells.flint.__FLINT_VERSION__},
        },
        "next_handoff": {
            "first": "Inject the correlated real-erfc q correction alongside the admitted PSI/complex-ERF correction into all 374 recurrence defects.",
            "second": "Use the surviving corrected defect as the quantitative target for endpoint-saddle, finite-ip, and Euler-Maclaurin remainder analysis.",
        },
        "proof_boundary": (
            "Rigorous intrinsic real-erfc residuals at 913 exact finite source arguments and their explicitly "
            "weighted within-call replacement only. This does not bound saddle locations, endpoint formula "
            "error, finite-ip or Euler-Maclaurin remainders, all binary128 arithmetic, recurrence accumulation, "
            "the outer Hardy representation, height uniformity, RH, or a prize-level conclusion."
        ),
    }
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(
        "built Hardy block-20 t5 real-erfc residual gate: "
        f"913 calls, max residual {artifact['aggregate']['maximum_erfc_residual_abs_upper']}, "
        f"max correlated q shift {artifact['aggregate']['maximum_total_q_delta_abs_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
