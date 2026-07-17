#!/usr/bin/env python3
"""Certify the exact-corridor localized curvature implication on 2<=u<=20."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
import json
import math
import os
from pathlib import Path
import sys
from time import perf_counter


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

from jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate import (  # noqa: E402
    sha256,
)
from jensen_window_pf_compound_order4_gaussian_cumulant_ray_target import (  # noqa: E402
    CORRIDOR_CAPS,
    CORRIDOR_K2,
    exact_interval,
)
from jensen_window_pf_compound_order4_localized_curvature_interval_core import (  # noqa: E402
    derivative_power,
    evaluate_localized_curvature_from_h_cover,
    signed_hurwitz_gamma_derivative,
)
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_interval,
    arb_lower_text,
    arb_rational,
)


DEFAULT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_chunks.jsonl"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_certificate.md"
)
SOURCE_CORRIDOR = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.json"
)
SOURCE_BRIDGE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_first_summand_curvature_bridge.json"
)
PRECISION_BITS = 256
CHUNK_BLOCKS = 20
DEFAULT_WORKERS = max(1, min(8, (os.cpu_count() or 4) - 1))


@dataclass(frozen=True)
class Segment:
    left: Fraction
    right: Fraction
    width: Fraction
    pad: Fraction


@dataclass(frozen=True)
class CurvatureRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


SEGMENTS = (
    Segment(Fraction(2), Fraction(23, 10), Fraction(1, 10_000), Fraction(1, 50_000)),
    Segment(Fraction(23, 10), Fraction(5, 2), Fraction(1, 1_000), Fraction(1, 50_000)),
    Segment(Fraction(5, 2), Fraction(3), Fraction(1, 1_000), Fraction(1, 100_000)),
    Segment(Fraction(3), Fraction(4), Fraction(1, 1_000), Fraction(1, 10**6)),
    Segment(Fraction(4), Fraction(5), Fraction(1, 1_000), Fraction(1, 10**8)),
    Segment(Fraction(5), Fraction(10), Fraction(1, 1_000), Fraction(1, 10**10)),
    Segment(Fraction(10), Fraction(20), Fraction(1, 1_000), Fraction(1, 10**19)),
)


def decimal_lower(value: Decimal, digits: int = 70) -> str:
    with localcontext() as context:
        context.prec = digits
        return format(value.next_minus(), "E")


def deterministic_blocks() -> list[tuple[Fraction, Fraction, Fraction, int]]:
    blocks = []
    for segment_index, segment in enumerate(SEGMENTS):
        count = (segment.right - segment.left) / segment.width
        if count.denominator != 1:
            raise RuntimeError(f"segment {segment_index} is not aligned")
        for offset in range(count.numerator):
            left = segment.left + offset * segment.width
            blocks.append((left, left + segment.width, segment.pad, segment_index))
    return blocks


BLOCKS = deterministic_blocks()


def deterministic_tasks() -> list[tuple[int, int, int]]:
    tasks = []
    first = 0
    chunk = 0
    while first < len(BLOCKS):
        count = min(CHUNK_BLOCKS, len(BLOCKS) - first)
        tasks.append((chunk, first, count))
        first += count
        chunk += 1
    return tasks


def initialize_worker() -> None:
    flint.ctx.prec = PRECISION_BITS


def corridor_h_derivatives(
    left: Fraction, right: Fraction, pad: Fraction
) -> tuple[dict[int, flint.arb], dict]:
    outer_left = left - pad
    outer_right = right + pad
    left_t = potential_jet_arb(arb_rational(left), 1)[1]
    right_t = potential_jet_arb(arb_rational(right), 1)[1]
    outer_left_t = potential_jet_arb(arb_rational(outer_left), 1)[1]
    outer_right_t = potential_jet_arb(arb_rational(outer_right), 1)[1]
    if not bool(outer_left_t < left_t - 2 and outer_right_t > right_t + 2):
        raise RuntimeError("outer mode pad does not cover the t+-2 collar")

    mode = arb_interval(outer_left, outer_right)
    jet = potential_jet_arb(mode, 8)
    t, curvature = jet[1], jet[2]
    if not bool(curvature > 0):
        raise RuntimeError("corridor mode curvature is not positive")
    q = flint.arb.pi() * (4 * mode).exp()
    kappa2 = exact_interval(
        1 + arb_rational(CORRIDOR_K2[0]) / q,
        1 + arb_rational(CORRIDOR_K2[1]) / q,
    )
    derivatives = {
        2: signed_hurwitz_gamma_derivative(2, t + flint.arb(1) / 2)
        - kappa2 / curvature
    }
    for order, cap in CORRIDOR_CAPS.items():
        baseline = math.factorial(order - 2) * q ** (
            1 - flint.arb(order) / 2
        )
        magnitude = exact_interval(baseline, arb_rational(cap) * baseline)
        cumulant = magnitude if order % 2 == 0 else -magnitude
        derivatives[order] = signed_hurwitz_gamma_derivative(
            order, t + flint.arb(1) / 2
        ) - cumulant / derivative_power(curvature, order)
    return derivatives, {
        "outer_mode": [str(outer_left), str(outer_right)],
        "left_t_pad": arb_lower_text(left_t - outer_left_t),
        "right_t_pad": arb_lower_text(outer_right_t - right_t),
    }


def certify_block(block_index: int) -> dict:
    left, right, pad, segment_index = BLOCKS[block_index]
    try:
        derivatives, diagnostics = corridor_h_derivatives(left, right, pad)
        result = evaluate_localized_curvature_from_h_cover(
            left,
            right,
            derivatives,
            cover_diagnostics=diagnostics,
        )
    except Exception as exc:
        return {
            "passed": False,
            "failure": "exception",
            "block_index": block_index,
            "mode": [str(left), str(right)],
            "detail": repr(exc),
        }
    if result.get("passed") is not True:
        return {
            "passed": False,
            "failure": result.get("failure"),
            "block_index": block_index,
            "mode": [str(left), str(right)],
            "result": result,
        }
    scaled_upper = Decimal(result["scaled_localized_upper"])
    scaled_margin = Decimal("3.5") - scaled_upper
    if scaled_margin <= 0:
        return {
            "passed": False,
            "failure": "nonpositive-scaled-margin",
            "block_index": block_index,
            "mode": [str(left), str(right)],
            "scaled_upper": str(scaled_upper),
        }
    return {
        "passed": True,
        "block_index": block_index,
        "segment_index": segment_index,
        "mode_left": str(left),
        "mode_right": str(right),
        "pad": str(pad),
        "scaled_localized_upper": result["scaled_localized_upper"],
        "scaled_margin_lower": decimal_lower(scaled_margin),
        "absolute_margin_lower": result["margin_lower"],
        "j_lower": result["j_lower"],
        "outer_mode": diagnostics["outer_mode"],
        "left_t_pad": diagnostics["left_t_pad"],
        "right_t_pad": diagnostics["right_t_pad"],
    }


def chunk_task(task: tuple[int, int, int]) -> dict:
    chunk_index, first_block, block_count = task
    selected = []
    maximum_scaled = None
    minimum_scaled_margin = None
    minimum_absolute_margin = None
    minimum_j = None
    minimum_t_pad = None
    for offset in range(block_count):
        block_index = first_block + offset
        row = certify_block(block_index)
        if row.get("passed") is not True:
            return {
                "kind": "exact_corridor_localized_curvature",
                "chunk_index": chunk_index,
                "first_block": first_block,
                "block_count": block_count,
                "passed": False,
                "failed_row": row,
            }
        if offset in {0, block_count // 2, block_count - 1}:
            selected.append(row)
        candidates = {
            "maximum_scaled": Decimal(row["scaled_localized_upper"]),
            "minimum_scaled_margin": Decimal(row["scaled_margin_lower"]),
            "minimum_absolute_margin": Decimal(row["absolute_margin_lower"]),
            "minimum_j": Decimal(row["j_lower"]),
            "minimum_t_pad": min(
                Decimal(row["left_t_pad"].replace("E", "e")),
                Decimal(row["right_t_pad"].replace("E", "e")),
            ),
        }
        if maximum_scaled is None or candidates["maximum_scaled"] > maximum_scaled[0]:
            maximum_scaled = (candidates["maximum_scaled"], block_index)
        for name, current in (
            ("minimum_scaled_margin", minimum_scaled_margin),
            ("minimum_absolute_margin", minimum_absolute_margin),
            ("minimum_j", minimum_j),
            ("minimum_t_pad", minimum_t_pad),
        ):
            candidate = candidates[name]
            if current is None or candidate < current[0]:
                value = (candidate, block_index)
                if name == "minimum_scaled_margin":
                    minimum_scaled_margin = value
                elif name == "minimum_absolute_margin":
                    minimum_absolute_margin = value
                elif name == "minimum_j":
                    minimum_j = value
                else:
                    minimum_t_pad = value
    assert all(
        value is not None
        for value in (
            maximum_scaled,
            minimum_scaled_margin,
            minimum_absolute_margin,
            minimum_j,
            minimum_t_pad,
        )
    )
    return {
        "kind": "exact_corridor_localized_curvature",
        "chunk_index": chunk_index,
        "first_block": first_block,
        "block_count": block_count,
        "mode_left": BLOCKS[first_block][0].__str__(),
        "mode_right": BLOCKS[first_block + block_count - 1][1].__str__(),
        "passed": True,
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "chunk_blocks": CHUNK_BLOCKS,
        },
        "extrema": {
            "maximum_scaled_localized_upper": format(maximum_scaled[0], "E"),
            "maximum_scaled_block": maximum_scaled[1],
            "minimum_scaled_margin_lower": format(minimum_scaled_margin[0], "E"),
            "minimum_scaled_margin_block": minimum_scaled_margin[1],
            "minimum_absolute_margin_lower": format(minimum_absolute_margin[0], "E"),
            "minimum_absolute_margin_block": minimum_absolute_margin[1],
            "minimum_j_lower": format(minimum_j[0], "E"),
            "minimum_j_block": minimum_j[1],
            "minimum_t_pad_lower": format(minimum_t_pad[0], "E"),
            "minimum_t_pad_block": minimum_t_pad[1],
        },
        "selected_rows": selected,
    }


def load_cache(path: Path, tasks: list[tuple[int, int, int]]) -> list[dict]:
    if not path.exists():
        return []
    records = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(records) > len(tasks):
        raise RuntimeError("localized-curvature cache has too many rows")
    for record, (chunk, first, count) in zip(records, tasks):
        if (
            record.get("kind") != "exact_corridor_localized_curvature"
            or record.get("chunk_index") != chunk
            or record.get("first_block") != first
            or record.get("block_count") != count
            or record.get("passed") is not True
            or record.get("parameters", {}).get("precision_bits") != PRECISION_BITS
        ):
            raise RuntimeError(f"localized-curvature cache mismatch at chunk {chunk}")
    return records


def build_cache(
    path: Path,
    tasks: list[tuple[int, int, int]],
    *,
    workers: int,
    overwrite: bool,
) -> list[dict]:
    if overwrite and path.exists():
        path.unlink()
    records = load_cache(path, tasks)
    remaining = tasks[len(records) :]
    if not remaining:
        return records
    path.parent.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    if workers == 1:
        initialize_worker()
        iterator = map(chunk_task, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=initialize_worker,
        )
        iterator = executor.map(chunk_task, remaining, chunksize=1)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for completed, record in enumerate(iterator, start=1):
                if record.get("passed") is not True:
                    raise RuntimeError(
                        "exact-corridor localized curvature failed: "
                        + json.dumps(record, sort_keys=True)
                    )
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 50 == 0:
                    handle.flush()
                    print(
                        f"exact-corridor curvature chunks: {len(records)}/{len(tasks)} "
                        f"({perf_counter()-started:.1f}s)"
                    )
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def aggregate_extremum(
    records: list[dict], value_key: str, block_key: str, *, maximum: bool
) -> dict:
    selected = (
        max(records, key=lambda row: Decimal(row["extrema"][value_key]))
        if maximum
        else min(records, key=lambda row: Decimal(row["extrema"][value_key]))
    )
    return {
        value_key: selected["extrema"][value_key],
        block_key: selected["extrema"][block_key],
    }


def build_artifact(records: list[dict], cache_path: Path) -> dict:
    tasks = deterministic_tasks()
    if len(records) != len(tasks):
        raise RuntimeError("localized-curvature artifact requires the complete cache")
    if records[0]["mode_left"] != "2" or records[-1]["mode_right"] != "20":
        raise RuntimeError("localized-curvature cache has the wrong endpoints")
    for previous, current in zip(records, records[1:]):
        if previous["mode_right"] != current["mode_left"]:
            raise RuntimeError("localized-curvature cache has a mode gap")
    corridor = json.loads(SOURCE_CORRIDOR.read_text(encoding="utf-8"))
    if corridor.get("summary", {}).get("global_exact_corridors_closed") is not True:
        raise RuntimeError("exact corridor source is not closed")
    bridge = json.loads(SOURCE_BRIDGE.read_text(encoding="utf-8"))
    if bridge.get("compact_curvature", {}).get("accepted_central_blocks") != 1073:
        raise RuntimeError("compact curvature source is not closed")

    extrema = {
        **aggregate_extremum(
            records,
            "maximum_scaled_localized_upper",
            "maximum_scaled_block",
            maximum=True,
        ),
        **aggregate_extremum(
            records,
            "minimum_scaled_margin_lower",
            "minimum_scaled_margin_block",
            maximum=False,
        ),
        **aggregate_extremum(
            records,
            "minimum_absolute_margin_lower",
            "minimum_absolute_margin_block",
            maximum=False,
        ),
        **aggregate_extremum(
            records,
            "minimum_j_lower",
            "minimum_j_block",
            maximum=False,
        ),
        **aggregate_extremum(
            records,
            "minimum_t_pad_lower",
            "minimum_t_pad_block",
            maximum=False,
        ),
    }
    segment_rows = [
        {
            "mode_interval": [str(segment.left), str(segment.right)],
            "block_width": str(segment.width),
            "outer_pad": str(segment.pad),
            "block_count": int((segment.right - segment.left) / segment.width),
        }
        for segment in SEGMENTS
    ]
    rows = [
        CurvatureRow(
            id="co4eclcfc_01_exact_corridor_input",
            role="exact_theorem_input",
            readiness="ready_to_apply",
            claim="All seven alternating-factorial exact cumulant corridors hold at every mode used by the cover.",
            formula="exact cumulant corridors hold globally on u>=2",
            proof_boundary="Inherited exact cumulant theorem only.",
        ),
        CurvatureRow(
            id="co4eclcfc_02_correlated_mode_cover",
            role="exact_interval_construction",
            readiness="ready_to_apply",
            claim="A nonuniform rational mode cover retains the correlation among q, t, curvature, and every corridor endpoint.",
            formula="width 10^-4 on [2,2.3]; width 10^-3 on [2.3,20]",
            proof_boundary="Finite mode-cover geometry only.",
            diagnostics={"segments": segment_rows},
        ),
        CurvatureRow(
            id="co4eclcfc_03_t_collar_cover",
            role="exact_interval_theorem",
            readiness="ready_to_apply",
            claim="Every rational outer pad covers the complete t+-2 derivative collar for its central mode block.",
            formula="t(u_left)-t(u_outer_left)>2 and t(u_outer_right)-t(u_right)>2",
            proof_boundary="Mode-to-t collar geometry only.",
            diagnostics={"minimum_t_pad_lower": extrema["minimum_t_pad_lower"]},
        ),
        CurvatureRow(
            id="co4eclcfc_04_localized_interval_evaluation",
            role="exact_interval_theorem",
            readiness="ready_to_apply",
            claim="Exact corridor boxes and Hurwitz-zeta derivatives keep the localized stable gap positive and the scaled localized curvature below 7/2 on every block.",
            formula="j_0>E_0 and t^2*U(t)<7/2",
            proof_boundary="Localized U inequality only.",
            diagnostics=extrema,
        ),
        CurvatureRow(
            id="co4eclcfc_05_finite_ray_curvature",
            role="exact_interval_theorem",
            readiness="ready_to_apply",
            claim="The first-summand stable-gap curvature ceiling holds throughout the complete finite ray segment.",
            formula="K_1(t)<=7/(2t^2) for every mode 2<=u<=20",
            proof_boundary="Finite mode segment only; u>=20 remains separate.",
        ),
        CurvatureRow(
            id="co4eclcfc_06_compact_composition",
            role="exact_theorem_composition",
            readiness="ready_to_apply",
            claim="Together with the inherited compact theorem, the first-summand curvature ceiling is proved from t=319 through the mode u=20.",
            formula="K_1(t)<=7/(2t^2), 319<=t<=V'(20)",
            proof_boundary="Finite real-parameter theorem only.",
        ),
        CurvatureRow(
            id="co4eclcfc_07_asymptotic_handoff",
            role="analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Prove the same exact-corridor localized inequality on the remaining asymptotic ray without extrapolating the finite cover.",
            formula="exact corridors => U(t)<=7/(2t^2), u>=20",
            proof_boundary="Open analytic ray; no global curvature or order-four entry is asserted.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_certificate",
        "date": "2026-07-13",
        "status": "rigorous exact-corridor localized-curvature theorem on 2<=u<=20",
        "proof_boundary": (
            "This artifact proves the exact cumulant corridors imply the localized "
            "first-summand curvature ceiling throughout 2<=u<=20 and composes that "
            "with the earlier compact theorem. It does not prove the u>=20 analytic "
            "ray, global order-four entry, PF-infinity, RH, or Lambda<=0."
        ),
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "mode_interval": ["2", "20"],
            "segments": segment_rows,
            "block_count": len(BLOCKS),
            "chunk_count": len(records),
            "chunk_blocks": CHUNK_BLOCKS,
        },
        "extrema": extrema,
        "cache": {
            "path": cache_path.relative_to(REPO_ROOT).as_posix(),
            "sha256": sha256(cache_path),
            "records": len(records),
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_rows": 6,
            "open_analytic_rows": 1,
            "mode_blocks": len(BLOCKS),
            "chunks": len(records),
            "corridors_used": 7,
            "t_collar_gates": 2 * len(BLOCKS),
            "positive_localized_blocks": len(BLOCKS),
            "finite_corridor_to_curvature_closed": True,
            "open_asymptotic_rays": 1,
        },
        "sources": [
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.md",
            "outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md",
            "outputs/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.md",
            "outputs/formal_core.md",
        ],
        "source_hashes": {
            SOURCE_CORRIDOR.relative_to(REPO_ROOT).as_posix(): sha256(SOURCE_CORRIDOR),
            SOURCE_BRIDGE.relative_to(REPO_ROOT).as_posix(): sha256(SOURCE_BRIDGE),
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_certificate.py"
        ),
        "remaining_target": (
            "Use q-leading mode geometry and the exact corridor constants to prove "
            "the localized U ceiling analytically for every u>=20."
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    extrema = artifact["extrema"]
    lines = [
        "# Jensen-Window PF Order-Four Exact-Corridor Localized-Curvature Finite Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous exact-corridor localized-curvature theorem on `2<=u<=20`.",
        "This is not a proof of the remaining asymptotic ray, order-four entry,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_certificate.json",
        "work/rh_compute/results/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_chunks.jsonl",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_certificate.py",
        "```",
        "",
        "## Correlated Cover",
        "",
        "The cover uses width `10^-4` on `2<=u<=2.3` and width `10^-3` on",
        "`2.3<=u<=20`. Each central block has a rational outer pad whose endpoint",
        "potential slopes certify the complete `t+-2` collar. The exact corridor",
        "boxes are evaluated with the same interval mode as `q`, `t`, curvature,",
        "and the Hurwitz-zeta derivatives, preserving their correlations.",
        "",
        "```text",
        f"mode blocks:                  {artifact['parameters']['block_count']}",
        f"t-collar endpoint gates:      {artifact['summary']['t_collar_gates']}",
        f"maximum t^2*U upper:          {extrema['maximum_scaled_localized_upper']}",
        f"minimum scaled margin:        {extrema['minimum_scaled_margin_lower']}",
        f"minimum localized J lower:    {extrema['minimum_j_lower']}",
        "```",
        "",
        "Every block proves `j_0>E_0` and `t^2 U(t)<7/2`. Therefore",
        "",
        "```text",
        "K_1(t)<=7/(2t^2) for every mode 2<=u<=20.",
        "```",
        "",
        "Composed with the earlier compact theorem, the ceiling now holds from",
        "`t=319` through `t=V'(20)`.",
        "",
        "## Remaining Boundary",
        "",
        "The only curvature segment left is `u>=20`. The finite cover is not",
        "extrapolated. The next theorem must use q-leading mode geometry and the",
        "exact corridor constants to prove the localized inequality analytically",
        "on that full ray.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.md",
        "outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite-cache", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    tasks = deterministic_tasks()
    records = build_cache(
        args.cache,
        tasks,
        workers=max(1, args.workers),
        overwrite=args.overwrite_cache,
    )
    artifact = build_artifact(records, args.cache)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    print(
        "certified order-four exact-corridor finite curvature: "
        "7 rows, 6 exact rows, 20700 mode blocks, 41400 t-collar gates, "
        "20700 positive localized blocks, 1 open asymptotic ray"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
