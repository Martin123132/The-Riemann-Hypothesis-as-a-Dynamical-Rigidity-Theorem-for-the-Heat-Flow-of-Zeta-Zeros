#!/usr/bin/env python3
"""Compare the source cutoff with a diagnostic nearest-root midpoint cutoff."""

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

from flint import arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_block_zero_source_formula_identity_gate as identity
import jensen_window_pf_newman_c1_hardy_t1e10_equation124_classical_block_partition_gate as partition


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_block_zero_equation62_cell_partition_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_block_zero_equation62_cell_partition_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISIONS = (70, 110)
EXPECTED_OUTPUTS = 15
M2 = identity.M2
N1 = identity.N1
ALPHAS = tuple(range(M2, N1 + 1, 2))
ROOT_ALPHAS = tuple(range(M2, N1 + 3, 2))
REQUESTED_TOLERANCE = arb("0.005")
EXPECTED_FINAL_DIAGNOSTIC_START = 37941
EXPECTED_FINAL_SOURCE_START = 37946


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


def continuous_lower_root(target_t: str, alpha: int, precision: int) -> arb:
    """Lower solution N of alpha=2N+t/(pi*N), paper equation (124)."""
    ctx.dps = precision
    t = arb(target_t)
    aa = arb(alpha)
    saddle = (8 * t / arb.pi()).sqrt()
    require(aa > saddle, f"alpha is not above the saddle: {alpha}")
    return aa * (1 - (1 - (saddle / aa) ** 2).sqrt()) / 4


def integer_certificates(target_t: str, precision: int) -> list[dict[str, Any]]:
    certificates: list[dict[str, Any]] = []
    previous: arb | None = None
    for alpha in ROOT_ALPHAS:
        root = continuous_lower_root(target_t, alpha, precision)
        nearest = math.floor(float(root.mid()) + 0.5)
        floor_n = math.floor(float(root.mid()))
        require(arb(nearest) - arb("0.5") < root < arb(nearest) + arb("0.5"), f"nearest-root ambiguity at alpha={alpha}")
        require(arb(floor_n) < root < arb(floor_n + 1), f"floor-root ambiguity at alpha={alpha}")
        if previous is not None:
            require(root < previous, f"root monotonicity lost at alpha={alpha}")
        previous = root
        certificates.append(
            {
                "alpha": alpha,
                "root_ball": root.str(precision, more=True),
                "nearest_integer": nearest,
                "distance_from_lower_half_integer_ball": (root - (arb(nearest) - arb("0.5"))).str(precision, more=True),
                "distance_from_upper_half_integer_ball": ((arb(nearest) + arb("0.5")) - root).str(precision, more=True),
                "floor_integer": floor_n,
                "distance_from_floor_ball": (root - floor_n).str(precision, more=True),
                "distance_to_next_integer_ball": (floor_n + 1 - root).str(precision, more=True),
            }
        )
    return certificates


def equation62_prefixes(offset: Decimal, precision: int) -> list[arb]:
    """Source-faithful normalized alpha prefixes, without the zero transition."""
    ctx.dps = precision
    ctx.threads = 1
    params = identity.central_parameters(precision)
    effective_offset = identity.default_real_shift(offset)
    delta = arb(format(effective_offset, "f"))
    pi = arb.pi()
    t = params["t"] + delta
    a = params["a"] + 4 * delta / (pi * params["a"])
    yphase = params["yphase"] + delta
    scale = (8 / a).sqrt()
    raw = arb(0)
    prefixes: list[arb] = []
    for alpha in ALPHAS:
        s = arb(alpha) / a
        b = (s * s - 1).sqrt()
        phase = t * (b * (b - s) + (s + b).log()) + yphase
        raw += phase.cos() / b.sqrt()
        prefixes.append(raw * scale)
    return prefixes


def sign_name(value: arb) -> str:
    if value > 0:
        return "positive"
    if value < 0:
        return "negative"
    return "contains_zero"


def output_atlas(offset: Decimal, precision: int) -> dict[str, Any]:
    target_t = format(Decimal("1e10") + offset, "f")
    roots = integer_certificates(target_t, precision)
    prefixes = equation62_prefixes(offset, precision)
    require(len(prefixes) == len(ALPHAS), "prefix roster drift")

    starts: set[int] = set()
    metadata: list[dict[str, int | bool]] = []
    for index, alpha in enumerate(ALPHAS):
        nearest = int(roots[index]["nearest_integer"])
        next_nearest = int(roots[index + 1]["nearest_integer"])
        midpoint_numerator = nearest + next_nearest
        midpoint_floor = midpoint_numerator // 2
        midpoint_ceil = (midpoint_numerator + 1) // 2
        source_start = int(roots[index]["floor_integer"]) + 1
        require(midpoint_floor <= midpoint_ceil <= source_start, f"cell ordering failed at alpha={alpha}")
        starts.update((midpoint_floor, midpoint_ceil, source_start))
        metadata.append(
            {
                "nearest_integer": nearest,
                "next_nearest_integer": next_nearest,
                "midpoint_numerator": midpoint_numerator,
                "midpoint_half_integral": bool(midpoint_numerator % 2),
                "midpoint_floor_start": midpoint_floor,
                "midpoint_ceil_start": midpoint_ceil,
                "source_endpoint_start": source_start,
            }
        )

    targets = partition.cumulative_targets(target_t, sorted(starts), precision)
    prefix_rows: list[dict[str, Any]] = []
    for index, alpha in enumerate(ALPHAS):
        meta = metadata[index]
        model = prefixes[index]
        floor_target = targets[int(meta["midpoint_floor_start"])]
        ceil_target = targets[int(meta["midpoint_ceil_start"])]
        source_target = targets[int(meta["source_endpoint_start"])]
        ceil_residual = model - ceil_target
        floor_residual = model - floor_target
        source_residual = model - source_target
        boundary_correction = ceil_target - source_target
        decomposition = ceil_residual + boundary_correction - source_residual
        require(decomposition.contains(0), f"boundary decomposition failed at prefix={index + 1}")
        prefix_rows.append(
            {
                "prefix_index": index + 1,
                "upper_alpha": alpha,
                **meta,
                "equation62_prefix_sum_ball": model.str(precision, more=True),
                "midpoint_floor_target_ball": floor_target.str(precision, more=True),
                "midpoint_ceil_target_ball": ceil_target.str(precision, more=True),
                "source_endpoint_target_ball": source_target.str(precision, more=True),
                "prefix_minus_midpoint_floor_target_ball": floor_residual.str(precision, more=True),
                "prefix_minus_midpoint_ceil_target_ball": ceil_residual.str(precision, more=True),
                "prefix_minus_source_endpoint_target_ball": source_residual.str(precision, more=True),
                "midpoint_ceil_target_minus_source_endpoint_target_ball": boundary_correction.str(precision, more=True),
                "boundary_decomposition_ball": decomposition.str(precision, more=True),
                "midpoint_ceil_residual_sign": sign_name(ceil_residual),
                "midpoint_ceil_residual_below_0p005": bool(abs(ceil_residual) < REQUESTED_TOLERANCE),
            }
        )
    return {"target_t": target_t, "root_certificates": roots, "prefix_rows": prefix_rows}


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    center = artifact["output_rows"][7]["final_prefix"]
    return f"""# Equation-(62) diagnostic midpoint-partition gate

Date: 2026-08-09

Status: finite diagnostic cutoff comparison validated; not a proof of RH

Paper equation (124) relates each alpha endpoint to a continuous lower root,
and equations (125)--(127) use an endpoint cutoff.  Neither the 2026 paper nor
the source prescribes the nearest-root Voronoi midpoint rule tested here.
This gate introduces that rule as a diagnostic alternative, extends it through
all {len(ALPHAS)} alpha terms, and certifies every root with Arb.

At `t=10^10`, the first nearest roots are `39852` and `39691`; their midpoint
is `39771.5`.  The introduced diagnostic therefore needs a lattice convention
at this height.  The artifact retains both floor and ceiling completions at
every half-integral boundary.  Its final boundary is not ambiguous: the roots for
alpha `159777` and `159779` are `37946` and `37936`, so their midpoint is the
integer `37941`.

The executable derives its complementary classical tail directly from the
terminal alpha endpoint, in agreement with equations (124)--(127).  It starts
at `37946`.
At the centre output the three relevant values are

```text
alpha prefix - source-endpoint target = {center['prefix_minus_source_endpoint_target_ball']}
alpha prefix - midpoint-cell target   = {center['prefix_minus_midpoint_ceil_target_ball']}
midpoint target - source target       = {center['midpoint_ceil_target_minus_source_endpoint_target_ball']}
```

The identity `source residual = midpoint residual + boundary correction` is
certified at every prefix and output.  Across all fifteen final prefixes, the
source-endpoint residual has magnitude above `0.005`, while the unambiguous
midpoint-cell residual has magnitude below `0.005`.  Its absolute range is

```text
{aggregate['minimum_final_midpoint_residual_absolute_ball']}
{aggregate['maximum_final_midpoint_residual_absolute_ball']}
```

Thus replacing the source cutoff by this diagnostic midpoint changes the
comparison by the five terms `37941..37945` and leaves a residual near
`0.00147`.  This is an exact finite comparison between two conventions.  It
does not show that the source omitted those terms or that the midpoint is the
published hybrid's correct boundary.

The smaller residual motivates a possible new cutoff theorem, but that theorem
must be derived independently and must transport the complete hybrid error.
The prefix atlas is diagnostic evidence only and has no RH implication.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    low_outputs: list[dict[str, Any]] = []
    high_outputs: list[dict[str, Any]] = []
    for output_index in range(1, EXPECTED_OUTPUTS + 1):
        offset = Decimal(output_index - 8) / Decimal(100)
        low_outputs.append(output_atlas(offset, PRECISIONS[0]))
        high_outputs.append(output_atlas(offset, PRECISIONS[1]))

    output_rows: list[dict[str, Any]] = []
    all_ceil_abs: list[arb] = []
    final_ceil_abs: list[arb] = []
    final_source_abs: list[arb] = []
    final_boundary_abs: list[arb] = []
    half_integral_count = 0
    ceil_below_count = 0
    for output_index, (low, high) in enumerate(zip(low_outputs, high_outputs), start=1):
        require(low["target_t"] == high["target_t"], "target precision roster drift")
        for low_root, high_root in zip(low["root_certificates"], high["root_certificates"]):
            require(low_root["alpha"] == high_root["alpha"], "root alpha precision drift")
            require(low_root["nearest_integer"] == high_root["nearest_integer"], "nearest integer precision drift")
            require(low_root["floor_integer"] == high_root["floor_integer"], "floor integer precision drift")
            require(arb(low_root["root_ball"]).overlaps(arb(high_root["root_ball"])), "root ball precision drift")
        for low_row, high_row in zip(low["prefix_rows"], high["prefix_rows"]):
            require(low_row["prefix_index"] == high_row["prefix_index"], "prefix precision roster drift")
            for key in (
                "nearest_integer",
                "next_nearest_integer",
                "midpoint_numerator",
                "midpoint_half_integral",
                "midpoint_floor_start",
                "midpoint_ceil_start",
                "source_endpoint_start",
            ):
                require(low_row[key] == high_row[key], f"discrete precision drift at prefix={low_row['prefix_index']}: {key}")
            for key in (
                "equation62_prefix_sum_ball",
                "midpoint_floor_target_ball",
                "midpoint_ceil_target_ball",
                "source_endpoint_target_ball",
                "prefix_minus_midpoint_floor_target_ball",
                "prefix_minus_midpoint_ceil_target_ball",
                "prefix_minus_source_endpoint_target_ball",
                "midpoint_ceil_target_minus_source_endpoint_target_ball",
            ):
                require(arb(low_row[key]).overlaps(arb(high_row[key])), f"ball precision drift at prefix={low_row['prefix_index']}: {key}")
            ceil_abs = abs(arb(high_row["prefix_minus_midpoint_ceil_target_ball"]))
            all_ceil_abs.append(ceil_abs)
            half_integral_count += int(high_row["midpoint_half_integral"])
            ceil_below_count += int(high_row["midpoint_ceil_residual_below_0p005"])

        final = high["prefix_rows"][-1]
        require(final["midpoint_ceil_start"] == EXPECTED_FINAL_DIAGNOSTIC_START, "final midpoint start drift")
        require(final["midpoint_floor_start"] == EXPECTED_FINAL_DIAGNOSTIC_START, "final midpoint unexpectedly ambiguous")
        require(final["source_endpoint_start"] == EXPECTED_FINAL_SOURCE_START, "final source start drift")
        midpoint_abs = abs(arb(final["prefix_minus_midpoint_ceil_target_ball"]))
        source_abs = abs(arb(final["prefix_minus_source_endpoint_target_ball"]))
        boundary_abs = abs(arb(final["midpoint_ceil_target_minus_source_endpoint_target_ball"]))
        require(midpoint_abs < REQUESTED_TOLERANCE, f"final midpoint target misses tolerance at output={output_index}")
        require(source_abs > REQUESTED_TOLERANCE, f"final source target no longer rejects tolerance at output={output_index}")
        final_ceil_abs.append(midpoint_abs)
        final_source_abs.append(source_abs)
        final_boundary_abs.append(boundary_abs)
        output_rows.append(
            {
                "output_index": output_index,
                "output_label": f"t{Decimal(output_index - 8) / Decimal(100):+.2f}",
                "target_t": high["target_t"],
                "root_certificates": high["root_certificates"],
                "prefix_rows": high["prefix_rows"],
                "final_prefix": final,
                "final_midpoint_residual_below_0p005": True,
                "final_source_endpoint_residual_above_0p005": True,
            }
        )

    aggregate = {
        "output_count": EXPECTED_OUTPUTS,
        "prefixes_per_output": len(ALPHAS),
        "total_prefix_comparisons": EXPECTED_OUTPUTS * len(ALPHAS),
        "total_root_certificates": EXPECTED_OUTPUTS * len(ROOT_ALPHAS),
        "half_integral_midpoint_count": half_integral_count,
        "midpoint_ceil_prefix_residual_below_0p005_count": ceil_below_count,
        "minimum_midpoint_ceil_prefix_residual_absolute_ball": min(all_ceil_abs, key=lambda x: x.lower()).str(PRECISIONS[1], more=True),
        "maximum_midpoint_ceil_prefix_residual_absolute_ball": max(all_ceil_abs, key=lambda x: x.upper()).str(PRECISIONS[1], more=True),
        "final_midpoint_residual_below_0p005_count": len(final_ceil_abs),
        "final_source_endpoint_residual_above_0p005_count": len(final_source_abs),
        "minimum_final_midpoint_residual_absolute_ball": min(final_ceil_abs, key=lambda x: x.lower()).str(PRECISIONS[1], more=True),
        "maximum_final_midpoint_residual_absolute_ball": max(final_ceil_abs, key=lambda x: x.upper()).str(PRECISIONS[1], more=True),
        "minimum_final_source_endpoint_residual_absolute_ball": min(final_source_abs, key=lambda x: x.lower()).str(PRECISIONS[1], more=True),
        "maximum_final_source_endpoint_residual_absolute_ball": max(final_source_abs, key=lambda x: x.upper()).str(PRECISIONS[1], more=True),
        "minimum_final_boundary_correction_absolute_ball": min(final_boundary_abs, key=lambda x: x.lower()).str(PRECISIONS[1], more=True),
        "maximum_final_boundary_correction_absolute_ball": max(final_boundary_abs, key=lambda x: x.upper()).str(PRECISIONS[1], more=True),
        "final_boundary_classical_terms": list(range(EXPECTED_FINAL_DIAGNOSTIC_START, EXPECTED_FINAL_SOURCE_START)),
    }

    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_block_zero_equation62_cell_partition_gate",
        "status": "finite_diagnostic_nearest_root_midpoint_comparison_validated",
        "passed": True,
        "scope": {
            "height_center": "1e10",
            "outputs": EXPECTED_OUTPUTS,
            "alpha_range": f"{M2}-{N1} odd",
            "root_alpha_range": f"{M2}-{N1 + 2} odd",
            "precisions_decimal_digits": list(PRECISIONS),
            "paper_equations": [62, 124, 125, 126, 127],
            "source_equations": [124, 125, 126, 127],
        },
        "conventions": {
            "diagnostic_midpoint": "Introduced here only: for alpha_j, take the nearest integer lower root of alpha=2N+t/(pi*N), then put the cumulative diagnostic boundary at (N_j+N_(j+1))/2.",
            "half_integral_midpoint_handling": "Both floor and ceiling lattice completions of the introduced diagnostic are retained; ceiling is the N>=midpoint set interpretation.",
            "source_endpoint": "Published/source convention: after terminal alpha A, CE is the continuous lower root, NC=floor(CE), and the complementary classical tail starts at NC+1.",
            "midpoint_rule_stated_in_paper_or_source": False,
            "final_midpoint_is_integral": True,
        },
        "output_rows": output_rows,
        "aggregate": aggregate,
        "decision": {
            "diagnostic_midpoint_rule_requires_lattice_completion_at_some_prefixes": half_integral_count > 0,
            "final_101_term_midpoint_is_ambiguous": False,
            "source_endpoint_and_diagnostic_midpoint_agree_at_final_prefix": False,
            "final_source_endpoint_residual_rejects_0p005_at_all_outputs": True,
            "final_diagnostic_midpoint_residual_satisfies_0p005_at_all_outputs": True,
            "five_term_cutoff_change_dominates_source_aligned_residual": True,
            "diagnostic_midpoint_is_a_published_cutoff_rule": False,
            "height_uniform_midpoint_remainder_proved": False,
            "finite_gate_has_rh_implication": False,
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "flint_threads": 1,
            "process_priority": priority,
        },
        "dependencies": {
            "source_identity_gate": {"path": relative(identity.RESULT), "sha256": file_hash(identity.RESULT)},
            "equation124_partition_gate": {"path": relative(partition.RESULT), "sha256": file_hash(partition.RESULT)},
            "paper": {"path": relative(partition.PAPER), "sha256": file_hash(partition.PAPER)},
            "source_file": {"path": relative(partition.SOURCE), "sha256": file_hash(partition.SOURCE)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": "Bound the source-aligned equations-(126)--(127) hybrid error under its endpoint convention. In parallel, derive an independent alternative-cutoff theorem before treating the diagnostic midpoint comparison as a replacement formula.",
        "proof_boundary": "Rigorous finite nearest-root diagnostic and prefix atlas at fifteen saved heights only. It compares a new midpoint convention with the published/source endpoint and observes a smaller residual; it does not repair the published hybrid, prove the midpoint convention valid, prove a height-uniform remainder, Lambda<=0, PF-infinity, RH, or a prize-level result.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "validated equation-(62) cell partition: "
        f"roots={aggregate['total_root_certificates']}, prefixes={aggregate['total_prefix_comparisons']}, "
        f"final-midpoint<0.005={aggregate['final_midpoint_residual_below_0p005_count']}/{EXPECTED_OUTPUTS}, "
        f"final-source>0.005={aggregate['final_source_endpoint_residual_above_0p005_count']}/{EXPECTED_OUTPUTS}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
