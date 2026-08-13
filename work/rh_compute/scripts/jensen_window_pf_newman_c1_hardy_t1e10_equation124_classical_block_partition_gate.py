#!/usr/bin/env python3
"""Certify the equation-(124) classical partition of every saved ZP block."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_classical_upper_main_remainder_split_gate as classical_split


CHECKPOINT = classical_split.CHECKPOINT
CLASSICAL_SPLIT = classical_split.RESULT
SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/telemetry/zeta14cubicmult_telemetry.f90"
PAPER = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_Brereton_2607.15310.pdf"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation124_classical_block_partition_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation124_classical_block_partition_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISIONS = (70, 110)
EXPECTED_OUTPUTS = 15
EXPECTED_BLOCKS = 36
LOWER_START = 622
UPPER_END = 39894
REQUESTED_TOLERANCE = arb("0.005")
SOURCE_CUTOFF_MARKERS = (
    "raecutoff=RN1",
    "CE=raecutoff*(1-sqrt(1-(a/(raecutoff*C))**2))/CE",
    "NC=aint(CE,dp1)",
    "I_START = 1 + (iteration-1) * N",
    "I_END   =         iteration * N",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def load_stage_zero() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with CHECKPOINT.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                row = json.loads(line, parse_float=str)
                if int(row["stage"]) == 0:
                    rows.append(row)
    require(len(rows) == EXPECTED_BLOCKS, "stage-zero block roster drift")
    require([int(row["unit"]) for row in rows] == list(range(EXPECTED_BLOCKS)), "stage-zero unit order drift")
    require(all(len(row["zsum"]) == EXPECTED_OUTPUTS for row in rows), "stage-zero output roster drift")
    require(all(int(row["totblock"]) == EXPECTED_BLOCKS - 1 for row in rows), "total-block metadata drift")
    return rows


def require_source_cutoff_semantics() -> None:
    source_text = SOURCE.read_text(encoding="utf-8")
    for marker in SOURCE_CUTOFF_MARKERS:
        require(marker in source_text, f"source cutoff marker drift: {marker}")


def alpha_of_n(target_t: str, n: int, precision: int) -> arb:
    ctx.dps = precision
    return 2 * arb(n) + arb(target_t) / (arb.pi() * arb(n))


def cutoff_certificate(target_t: str, rn1: str, precision: int) -> dict[str, Any]:
    ctx.dps = precision
    t_float = float(Decimal(target_t))
    rn_float = float(Decimal(rn1))
    a_float = math.sqrt(8.0 * t_float / math.pi)
    lower_root_float = rn_float * (1.0 - math.sqrt(1.0 - (a_float / rn_float) ** 2)) / 4.0
    source_nc = math.floor(lower_root_float)
    cutoff = source_nc + 1
    require(LOWER_START <= cutoff <= UPPER_END, f"cutoff outside classical upper range: {cutoff}")
    rn = arb(rn1)
    t = arb(target_t)
    a = (8 * t / arb.pi()).sqrt()
    source_ce = rn * (1 - (1 - (a / rn) ** 2).sqrt()) / 4
    require(arb(source_nc) < source_ce, f"source CE lower inequality failed at n={source_nc}")
    require(source_ce < arb(cutoff), f"source CE upper inequality failed at n={cutoff}")
    alpha_at_cutoff = alpha_of_n(target_t, cutoff, precision)
    alpha_before = alpha_of_n(target_t, cutoff - 1, precision)
    require(alpha_at_cutoff < rn, f"lower cutoff inequality failed at n={cutoff}")
    require(rn < alpha_before, f"upper cutoff inequality failed at n={cutoff}")
    return {
        "precision_decimal_digits": precision,
        "cutoff_n": cutoff,
        "source_ce_ball": source_ce.str(precision, more=True),
        "source_nc_aint_positive": source_nc,
        "source_tail_start_n": cutoff,
        "source_ce_minus_nc_ball": (source_ce - source_nc).str(precision, more=True),
        "source_nc_plus_one_minus_ce_ball": (cutoff - source_ce).str(precision, more=True),
        "alpha_at_cutoff_ball": alpha_at_cutoff.str(precision, more=True),
        "alpha_before_cutoff_ball": alpha_before.str(precision, more=True),
        "rn1_minus_alpha_at_cutoff_ball": (rn - alpha_at_cutoff).str(precision, more=True),
        "alpha_before_cutoff_minus_rn1_ball": (alpha_before - rn).str(precision, more=True),
    }


def cumulative_targets(target_t: str, cutoffs: list[int], precision: int) -> dict[int, arb]:
    ctx.dps = precision
    ctx.threads = 1
    t = arb(target_t)
    theta = acb(arb("0.25"), t / 2).lgamma().imag - t * arb.pi().log() / 2
    wanted = set(cutoffs)
    targets: dict[int, arb] = {}
    total = arb(0)
    for n in range(UPPER_END, LOWER_START - 1, -1):
        nn = arb(n)
        total += 2 * (theta - t * nn.log()).cos() / nn.sqrt()
        if n in wanted:
            targets[n] = total
    require(set(targets) == wanted, "classical target cutoff roster drift")
    return targets


def sign_name(value: arb) -> str:
    if value > 0:
        return "positive"
    if value < 0:
        return "negative"
    return "contains_zero"


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Equation-(124) classical block-partition atlas

Date: 2026-08-09

Status: finite interval diagnostic validated; not a proof of RH

The source checkpoint journal preserves cumulative `ZP` values after block
zero and every block through 35.  Paper equation (124) supplies the coordinate

```text
alpha(n)=2n+t/(pi*n).
```

The executable fixes the discrete convention: it assigns `raecutoff=RN1`,
evaluates the lower root `CE` in equation (124), sets `NC=AINT(CE)`, and sums
the direct Riemann--Siegel part over `n=1,...,NC`.  Because `CE>0`, `AINT` is
the floor operation here, so the complementary classical tail starts at
`N_k=NC+1`.  No `RN1+1` cell-edge shift occurs in the source.

For each saved alpha endpoint `RN1_k`, the gate certifies the unique integer
`N_k` satisfying

```text
alpha(N_k)<RN1_k<alpha(N_k-1).
```

It then defines an exact classical comparison target, independently of the
non-exact hybrid formula,

```text
T_k(t)=2 sum_(n=N_k)^39894 cos(theta(t)-t log n)/sqrt(n).
```

The final cutoff is `N_35=622`, so `T_35` is exactly the full classical upper
main sum used in Section 11.298.  The 36 cutoffs partition all 39,273 terms
without overlap or omission.  Every cutoff inequality and all
{aggregate['row_count']} cumulative comparisons overlap at 70 and 110 decimal
digits; an independent checker repeats them at higher precision.

Block zero is already discrepant on all fifteen outputs.  Its signed source
minus exact-target residual is strictly negative and exceeds `0.005` in
absolute value on all fifteen:

```text
minimum |E_0| = {aggregate['minimum_block_zero_residual_absolute_ball']}
maximum |E_0| = {aggregate['maximum_block_zero_residual_absolute_ball']}
```

The net error added by blocks 1--35 is positive on all fifteen outputs, so it
partially cancels block zero rather than creating the final discrepancy.  The
saved final residual still matches the independently certified classical
upper-main residual on every output.

This atlas identifies the earliest block-aligned defect but not its analytic
cause.  Block zero combines the transition value, direct `ter` terms, and the
final `sqrt(8/a)` normalization.  Those columns must be separated against the
exact `n=N_0,...,39894` target.  Equation (124) is used only to define exact
integer boundaries; the paper explicitly says the hybrid representation is
not exact.  The atlas is finite at the saved heights and has no RH implication.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    require_source_cutoff_semantics()
    source_blocks = load_stage_zero()
    prior = json.loads(CLASSICAL_SPLIT.read_text(encoding="utf-8"))
    require(len(prior["output_rows"]) == EXPECTED_OUTPUTS, "classical-split output roster drift")

    output_rows: list[dict[str, Any]] = []
    cutoff_rosters: list[list[int]] = []
    all_cumulative_residuals: list[arb] = []
    all_increment_residuals: list[arb] = []
    block_zero_residuals: list[arb] = []
    final_residuals: list[arb] = []
    later_net_residuals: list[arb] = []

    for output_index in range(1, EXPECTED_OUTPUTS + 1):
        offset = Decimal(output_index - 8) / Decimal(100)
        target_t = format(Decimal("1e10") + offset, "f")
        low_certs = [cutoff_certificate(target_t, row["rn1"], PRECISIONS[0]) for row in source_blocks]
        high_certs = [cutoff_certificate(target_t, row["rn1"], PRECISIONS[1]) for row in source_blocks]
        cutoffs = [int(cert["cutoff_n"]) for cert in high_certs]
        require([int(cert["cutoff_n"]) for cert in low_certs] == cutoffs, "cutoff precision drift")
        require(all(cutoffs[k] < cutoffs[k - 1] for k in range(1, EXPECTED_BLOCKS)), "cutoff roster is not strictly decreasing")
        require(cutoffs[-1] == LOWER_START, "final cutoff does not close the classical upper sum")
        cutoff_rosters.append(cutoffs)

        low_targets = cumulative_targets(target_t, cutoffs, PRECISIONS[0])
        high_targets = cumulative_targets(target_t, cutoffs, PRECISIONS[1])
        ctx.dps = PRECISIONS[1]
        block_rows: list[dict[str, Any]] = []
        previous_source = arb(0)
        previous_target = arb(0)
        previous_residual = arb(0)
        for block, (source_row, low_cert, high_cert, cutoff) in enumerate(zip(source_blocks, low_certs, high_certs, cutoffs)):
            low_target = low_targets[cutoff]
            high_target = high_targets[cutoff]
            require(low_target.overlaps(high_target), f"target precision ladder missed at output {output_index}, block {block}")
            source_value = arb(source_row["zsum"][output_index - 1])
            cumulative_residual = source_value - high_target
            source_increment = source_value - previous_source
            target_increment = high_target - previous_target
            increment_residual = source_increment - target_increment
            increment_identity = increment_residual - (cumulative_residual - previous_residual)
            require(increment_identity.contains(0), f"increment identity excludes zero at output {output_index}, block {block}")
            block_rows.append(
                {
                    "block": block,
                    "source_rn1_binary128_decimal": source_row["rn1"],
                    "cutoff_n": cutoff,
                    "exact_target_term_count": UPPER_END - cutoff + 1,
                    "exact_block_term_count": UPPER_END - cutoff + 1 if block == 0 else cutoffs[block - 1] - cutoff,
                    "low_precision": {
                        "precision_decimal_digits": PRECISIONS[0],
                        "cutoff_certificate": low_cert,
                        "exact_cumulative_target_ball": low_target.str(PRECISIONS[0], more=True),
                    },
                    "high_precision": {
                        "precision_decimal_digits": PRECISIONS[1],
                        "cutoff_certificate": high_cert,
                        "exact_cumulative_target_ball": high_target.str(PRECISIONS[1], more=True),
                    },
                    "precision_enclosures_overlap": True,
                    "source_cumulative_zp_binary128_decimal": source_row["zsum"][output_index - 1],
                    "source_increment_ball": source_increment.str(PRECISIONS[1], more=True),
                    "exact_target_increment_ball": target_increment.str(PRECISIONS[1], more=True),
                    "cumulative_residual_ball": cumulative_residual.str(PRECISIONS[1], more=True),
                    "cumulative_residual_sign": sign_name(cumulative_residual),
                    "increment_residual_ball": increment_residual.str(PRECISIONS[1], more=True),
                    "increment_residual_sign": sign_name(increment_residual),
                    "increment_identity_residual_ball": increment_identity.str(PRECISIONS[1], more=True),
                }
            )
            all_cumulative_residuals.append(abs(cumulative_residual))
            all_increment_residuals.append(abs(increment_residual))
            previous_source = source_value
            previous_target = high_target
            previous_residual = cumulative_residual

        block_zero = arb(block_rows[0]["cumulative_residual_ball"])
        final_residual = arb(block_rows[-1]["cumulative_residual_ball"])
        later_net = final_residual - block_zero
        require(block_zero < -REQUESTED_TOLERANCE, f"block zero does not reject tolerance at output {output_index}")
        require(later_net > 0, f"later blocks do not partially cancel block zero at output {output_index}")
        prior_final = arb(prior["output_rows"][output_index - 1]["hybrid_minus_classical_main_ball"])
        require(final_residual.overlaps(prior_final), f"final residual misses prior split at output {output_index}")
        first_cumulative_nonzero = next(row["block"] for row in block_rows if row["cumulative_residual_sign"] != "contains_zero")
        first_increment_nonzero = next(row["block"] for row in block_rows if row["increment_residual_sign"] != "contains_zero")
        require(first_cumulative_nonzero == 0 and first_increment_nonzero == 0, "earliest discrepancy is not block zero")
        block_zero_residuals.append(abs(block_zero))
        final_residuals.append(abs(final_residual))
        later_net_residuals.append(abs(later_net))
        output_rows.append(
            {
                "output_index": output_index,
                "output_label": f"t{offset:+.2f}",
                "target_t": target_t,
                "first_nonzero_cumulative_residual_block": first_cumulative_nonzero,
                "first_nonzero_increment_residual_block": first_increment_nonzero,
                "block_zero_rejects_0p005": True,
                "later_blocks_net_partially_cancel_block_zero": True,
                "later_blocks_net_residual_ball": later_net.str(PRECISIONS[1], more=True),
                "later_net_to_block_zero_absolute_ratio_ball": (abs(later_net) / abs(block_zero)).str(PRECISIONS[1], more=True),
                "final_residual_matches_classical_split": True,
                "block_rows": block_rows,
            }
        )

    canonical_cutoffs = cutoff_rosters[0]
    require(all(roster == canonical_cutoffs for roster in cutoff_rosters), "cutoffs vary across shifted heights")
    block_aggregates: list[dict[str, Any]] = []
    for block in range(EXPECTED_BLOCKS):
        rows = [output["block_rows"][block] for output in output_rows]
        increment_balls = [abs(arb(row["increment_residual_ball"])) for row in rows]
        cumulative_balls = [abs(arb(row["cumulative_residual_ball"])) for row in rows]
        block_aggregates.append(
            {
                "block": block,
                "cutoff_n": canonical_cutoffs[block],
                "exact_block_term_count": rows[0]["exact_block_term_count"],
                "cumulative_negative_count": sum(row["cumulative_residual_sign"] == "negative" for row in rows),
                "cumulative_positive_count": sum(row["cumulative_residual_sign"] == "positive" for row in rows),
                "cumulative_contains_zero_count": sum(row["cumulative_residual_sign"] == "contains_zero" for row in rows),
                "increment_negative_count": sum(row["increment_residual_sign"] == "negative" for row in rows),
                "increment_positive_count": sum(row["increment_residual_sign"] == "positive" for row in rows),
                "increment_contains_zero_count": sum(row["increment_residual_sign"] == "contains_zero" for row in rows),
                "minimum_cumulative_residual_absolute_ball": min(cumulative_balls, key=lambda x: x.lower()).str(PRECISIONS[1], more=True),
                "maximum_cumulative_residual_absolute_ball": max(cumulative_balls, key=lambda x: x.upper()).str(PRECISIONS[1], more=True),
                "minimum_increment_residual_absolute_ball": min(increment_balls, key=lambda x: x.lower()).str(PRECISIONS[1], more=True),
                "maximum_increment_residual_absolute_ball": max(increment_balls, key=lambda x: x.upper()).str(PRECISIONS[1], more=True),
            }
        )

    later_ratios = [arb(output["later_net_to_block_zero_absolute_ratio_ball"]) for output in output_rows]
    aggregate = {
        "output_count": EXPECTED_OUTPUTS,
        "block_count": EXPECTED_BLOCKS,
        "row_count": EXPECTED_OUTPUTS * EXPECTED_BLOCKS,
        "precision_overlap_count": sum(row["precision_enclosures_overlap"] for output in output_rows for row in output["block_rows"]),
        "cutoff_inequality_certificate_count": EXPECTED_OUTPUTS * EXPECTED_BLOCKS,
        "source_ce_floor_certificate_count": EXPECTED_OUTPUTS * EXPECTED_BLOCKS,
        "source_cutoff_semantics_markers_verified": len(SOURCE_CUTOFF_MARKERS),
        "cutoff_roster_height_invariant": True,
        "cutoff_roster": canonical_cutoffs,
        "partition_total_terms": sum(row["exact_block_term_count"] for row in output_rows[0]["block_rows"]),
        "first_nonzero_cumulative_is_block_zero_count": sum(output["first_nonzero_cumulative_residual_block"] == 0 for output in output_rows),
        "first_nonzero_increment_is_block_zero_count": sum(output["first_nonzero_increment_residual_block"] == 0 for output in output_rows),
        "block_zero_strict_negative_count": sum(output["block_rows"][0]["cumulative_residual_sign"] == "negative" for output in output_rows),
        "block_zero_rejects_0p005_count": sum(output["block_zero_rejects_0p005"] for output in output_rows),
        "minimum_block_zero_residual_absolute_ball": min(block_zero_residuals, key=lambda x: x.lower()).str(PRECISIONS[1], more=True),
        "maximum_block_zero_residual_absolute_ball": max(block_zero_residuals, key=lambda x: x.upper()).str(PRECISIONS[1], more=True),
        "later_blocks_net_positive_count": sum(output["later_blocks_net_partially_cancel_block_zero"] for output in output_rows),
        "minimum_later_net_residual_absolute_ball": min(later_net_residuals, key=lambda x: x.lower()).str(PRECISIONS[1], more=True),
        "maximum_later_net_residual_absolute_ball": max(later_net_residuals, key=lambda x: x.upper()).str(PRECISIONS[1], more=True),
        "minimum_later_net_to_block_zero_ratio_ball": min(later_ratios, key=lambda x: x.lower()).str(PRECISIONS[1], more=True),
        "maximum_later_net_to_block_zero_ratio_ball": max(later_ratios, key=lambda x: x.upper()).str(PRECISIONS[1], more=True),
        "final_residual_matches_classical_split_count": sum(output["final_residual_matches_classical_split"] for output in output_rows),
        "minimum_final_residual_absolute_ball": min(final_residuals, key=lambda x: x.lower()).str(PRECISIONS[1], more=True),
        "maximum_final_residual_absolute_ball": max(final_residuals, key=lambda x: x.upper()).str(PRECISIONS[1], more=True),
    }
    require(aggregate["precision_overlap_count"] == aggregate["row_count"], "precision aggregate drift")
    require(aggregate["partition_total_terms"] == UPPER_END - LOWER_START + 1, "classical partition does not close")
    require(aggregate["block_zero_rejects_0p005_count"] == EXPECTED_OUTPUTS, "block-zero tolerance aggregate drift")
    require(aggregate["later_blocks_net_positive_count"] == EXPECTED_OUTPUTS, "later cancellation aggregate drift")
    require(aggregate["final_residual_matches_classical_split_count"] == EXPECTED_OUTPUTS, "final comparison aggregate drift")

    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation124_classical_block_partition_gate",
        "status": "rigorous_fifteen_output_thirty_six_block_equation124_induced_classical_partition_atlas",
        "passed": True,
        "scope": {
            "height_center": "1e10",
            "outputs": EXPECTED_OUTPUTS,
            "blocks": "0-35",
            "classical_upper_index_range": f"{LOWER_START}-{UPPER_END}",
            "precisions_decimal_digits": list(PRECISIONS),
        },
        "partition_definition": {
            "paper_equation": "124",
            "alpha_of_n": "alpha(n)=2*n+t/(pi*n)",
            "source_code_rule": "raecutoff=RN1; CE=equation124_lower_root(RN1); NC=AINT(CE); direct RS sum n=1,...,NC; N_k=NC+1",
            "source_rn1_endpoint_used_without_cell_edge_shift": True,
            "cutoff_rule": "alpha(N_k)<RN1_k<alpha(N_k-1)",
            "exact_cumulative_target": "T_k=2*sum_(n=N_k)^39894 cos(theta(t)-t*log(n))/sqrt(n)",
            "exact_increment_target": "Delta T_0=T_0; Delta T_k=T_k-T_(k-1)",
            "residual_identity": "E_k=Q_k-T_k; Delta E_k=E_k-E_(k-1)=Delta Q_k-Delta T_k",
        },
        "block_aggregates": block_aggregates,
        "output_rows": output_rows,
        "aggregate": aggregate,
        "decision": {
            "source_discrete_boundary_convention_resolved": True,
            "earliest_block_aligned_discrepancy": 0,
            "block_zero_alone_rejects_saved_tolerance_on_all_outputs": True,
            "blocks_1_through_35_net_partially_cancel_block_zero_on_all_outputs": True,
            "final_discrepancy_created_by_exact_remaining_rs_correction": False,
            "finite_atlas_is_height_uniform_theorem": False,
            "finite_atlas_has_rh_implication": False,
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "flint_threads": 1,
            "process_priority": priority,
        },
        "dependencies": {
            "source_checkpoint": {"path": relative(CHECKPOINT), "sha256": file_hash(CHECKPOINT)},
            "classical_upper_split": {"path": relative(CLASSICAL_SPLIT), "sha256": file_hash(CLASSICAL_SPLIT)},
            "source_file": {"path": relative(SOURCE), "sha256": file_hash(SOURCE)},
            "paper": {"path": relative(PAPER), "sha256": file_hash(PAPER)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": "Split source block zero into transition, direct ter-sum, and sqrt(8/a) normalization columns and compare each with the exact n=N_0,...,39894 classical target using an exact identity or explicit remainder.",
        "proof_boundary": "Rigorous finite equation-(124)-induced classical partition at fifteen saved heights only. It identifies block zero as the earliest block-aligned discrepancy but does not identify its internal analytic cause, prove the non-exact hybrid representation, establish height uniformity, prove Lambda<=0 or PF-infinity, prove RH, or reach a prize-level conclusion.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "validated equation-124 block partition: "
        f"{aggregate['row_count']} rows, block0-rejects-0.005={aggregate['block_zero_rejects_0p005_count']}, "
        f"later-net-cancels={aggregate['later_blocks_net_positive_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
