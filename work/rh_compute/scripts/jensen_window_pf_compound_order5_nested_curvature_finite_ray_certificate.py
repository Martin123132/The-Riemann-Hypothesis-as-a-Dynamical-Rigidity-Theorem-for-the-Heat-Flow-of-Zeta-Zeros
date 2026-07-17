#!/usr/bin/env python3
"""Certify the nested first-summand order-five curvature on 2<=u<=20."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
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

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as order4_compact  # noqa: E402
from jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_certificate import (  # noqa: E402
    corridor_h_derivatives,
)
from jensen_window_pf_compound_order5_nested_curvature_interval_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    evaluate_nested_curvature_from_h_cover,
)
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
)


DEFAULT_EXTENSION_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order5_nested_curvature_u2_extension_tiles.jsonl"
)
DEFAULT_RAY_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order5_nested_curvature_u2_u20_blocks.jsonl"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.md"
)
SOURCE_EXACT_CORRIDOR = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.json"
)
SOURCE_COMPACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order5_nested_curvature_compact_certificate.json"
)
SOURCE_BASE_CACHE = order4_compact.DEFAULT_CACHE
PRECISION_BITS = 256
MODE_TWO = Fraction(2)
COLLAR_END = Fraction(2001, 1000)
EXTENSION_START = order4_compact.OUTER_END
EXTENSION_END = Fraction(100051, 50000)
EXTENSION_ENVELOPE_END = Fraction(1001, 500)
BASE_TAIL_START = Fraction(19999, 10000)
CORRIDOR_PAD = Fraction(1, 1000)
COLLAR_T = 3
DEFAULT_WORKERS = max(1, min(8, (os.cpu_count() or 4) - 1))


@dataclass(frozen=True)
class Segment:
    left: Fraction
    right: Fraction
    width: Fraction


@dataclass(frozen=True)
class RayRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


SEGMENTS = (
    Segment(COLLAR_END, Fraction(2501, 1000), Fraction(1, 200)),
    Segment(Fraction(2501, 1000), Fraction(20), Fraction(1, 100)),
)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def source_contract() -> dict:
    corridor = load_json(SOURCE_EXACT_CORRIDOR)
    compact = load_json(SOURCE_COMPACT)
    if corridor.get("summary", {}).get("global_exact_corridors_closed") is not True:
        raise RuntimeError("exact cumulant corridor source is not closed")
    if compact.get("compact", {}).get("certified_t_range") != "320<=t<=V'(2)":
        raise RuntimeError("nested compact source contract changed")
    base_hash = sha256(SOURCE_BASE_CACHE)
    expected_hash = compact.get("source_contract", {}).get("cache_sha256")
    if base_hash != expected_hash:
        raise RuntimeError("base H-derivative cache hash changed")
    return {
        "exact_corridor_status": corridor.get("status"),
        "exact_corridor_formula": (
            "seven alternating-factorial exact cumulant corridors hold for every u>=2"
        ),
        "compact_range": compact["compact"]["certified_t_range"],
        "base_cache": SOURCE_BASE_CACHE.relative_to(REPO_ROOT).as_posix(),
        "base_cache_sha256": base_hash,
    }


def extension_tasks() -> list[tuple[int, Fraction, Fraction]]:
    start_index = order4_compact.mode_index(EXTENSION_START)
    tasks = []
    left = EXTENSION_START
    index = start_index
    while left < EXTENSION_END:
        right = min(left + order4_compact.TILE_WIDTH, EXTENSION_END)
        tasks.append((index, left, right))
        left = right
        index += 1
    return tasks


def ray_tasks() -> list[tuple[int, Fraction, Fraction]]:
    tasks: list[tuple[int, Fraction, Fraction]] = []
    index = 0
    for segment in SEGMENTS:
        left = segment.left
        while left < segment.right:
            right = min(left + segment.width, segment.right)
            tasks.append((index, left, right))
            left = right
            index += 1
    return tasks


def initialize_worker() -> None:
    flint.ctx.prec = PRECISION_BITS


def extension_task(task: tuple[int, Fraction, Fraction]) -> dict:
    return order4_compact.tile_task(task)


def ray_task(task: tuple[int, Fraction, Fraction]) -> dict:
    index, left, right = task
    flint.ctx.prec = PRECISION_BITS
    outer_left = left - CORRIDOR_PAD
    outer_right = right + CORRIDOR_PAD
    if outer_left < MODE_TWO:
        return {
            "index": index,
            "mode_left": str(left),
            "mode_right": str(right),
            "passed": False,
            "failure": "outer-mode-below-two",
        }
    try:
        derivatives, diagnostics = corridor_h_derivatives(
            left, right, CORRIDOR_PAD
        )
        left_t = potential_jet_arb(arb_rational(left), 1)[1]
        right_t = potential_jet_arb(arb_rational(right), 1)[1]
        outer_left_t = potential_jet_arb(arb_rational(outer_left), 1)[1]
        outer_right_t = potential_jet_arb(arb_rational(outer_right), 1)[1]
        if not bool(
            outer_left_t < left_t - COLLAR_T
            and outer_right_t > right_t + COLLAR_T
        ):
            raise RuntimeError("corridor mode pad does not cover t+-3")
        diagnostics.update(
            {
                "left_t_collar": (left_t - outer_left_t).str(35).replace(
                    "e", "E"
                ),
                "right_t_collar": (outer_right_t - right_t).str(35).replace(
                    "e", "E"
                ),
            }
        )
        result = evaluate_nested_curvature_from_h_cover(
            left,
            right,
            derivatives,
            cover_diagnostics=diagnostics,
        )
    except Exception as exc:
        return {
            "index": index,
            "mode_left": str(left),
            "mode_right": str(right),
            "passed": False,
            "failure": "exception",
            "detail": repr(exc),
        }
    return {
        "kind": "order5_nested_curvature_exact_corridor_block",
        "index": index,
        "mode_left": str(left),
        "mode_right": str(right),
        **result,
    }


def load_cache(
    path: Path,
    tasks: list[tuple[int, Fraction, Fraction]],
) -> list[dict]:
    if not path.exists():
        return []
    records = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(records) > len(tasks):
        raise RuntimeError(f"cache {path.name} has too many rows")
    for record, (index, left, right) in zip(records, tasks):
        if (
            record.get("index") != index
            or record.get("mode_left") != str(left)
            or record.get("mode_right") != str(right)
            or record.get("passed") is not True
        ):
            raise RuntimeError(f"cache prefix mismatch at {path.name}:{index}")
    return records


def build_cache(
    path: Path,
    tasks: list[tuple[int, Fraction, Fraction]],
    worker,
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
    start = perf_counter()
    if workers == 1:
        iterator = map(worker, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers, initializer=initialize_worker
        )
        iterator = executor.map(worker, remaining, chunksize=8)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for completed, record in enumerate(iterator, start=1):
                if record.get("passed") is not True:
                    raise RuntimeError(
                        f"cache row {record.get('index')} failed: "
                        f"{record.get('failure')} {record.get('detail', '')}"
                    )
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 250 == 0:
                    handle.flush()
                    print(
                        f"{path.stem}: {len(records)}/{len(tasks)} "
                        f"({perf_counter() - start:.1f}s)"
                    )
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def extension_envelope() -> dict:
    flint.ctx.prec = PRECISION_BITS
    mode = (
        (arb_rational(MODE_TWO) + arb_rational(EXTENSION_ENVELOPE_END)) / 2
        + flint.arb(
            0,
            (arb_rational(EXTENSION_ENVELOPE_END) - arb_rational(MODE_TWO))
            / 2,
        )
    )
    jet = potential_jet_arb(mode, 8)
    ratio = jet[8] / jet[2] ** 4
    cap = flint.arb(1) / 50_000
    if not bool(abs(ratio) < cap):
        raise RuntimeError("extension V8/a4 envelope failed")
    return {
        "mode_interval": [str(MODE_TWO), str(EXTENSION_ENVELOPE_END)],
        "ratio_ball": ratio.str(60).replace("e", "E"),
        "cap": "1/50000",
        "strict": True,
    }


def base_tail_records() -> list[dict]:
    rows = []
    with SOURCE_BASE_CACHE.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip():
                continue
            row = json.loads(raw)
            if Fraction(row["mode_left"]) >= BASE_TAIL_START:
                rows.append(row)
    if not rows:
        raise RuntimeError("base H cache has no u=2 tail rows")
    if rows[0]["mode_left"] != str(BASE_TAIL_START):
        raise RuntimeError("base H cache tail starts at the wrong mode")
    if rows[-1]["mode_right"] != str(EXTENSION_START):
        raise RuntimeError("base H cache tail ends at the wrong mode")
    return rows


def certify_initial_collar(
    extension: list[dict], envelope: dict
) -> dict:
    flint.ctx.prec = PRECISION_BITS
    records = base_tail_records() + extension
    derivatives, diagnostics = order4_compact.derivative_cover(
        records, 0, len(records) - 1
    )
    left_t = potential_jet_arb(arb_rational(MODE_TWO), 1)[1]
    right_t = potential_jet_arb(arb_rational(COLLAR_END), 1)[1]
    outer_left_t = order4_compact.interval_from_text(records[0]["t_left"])
    outer_right_t = order4_compact.interval_from_text(records[-1]["t_right"])
    if not bool(
        outer_left_t < left_t - COLLAR_T
        and outer_right_t > right_t + COLLAR_T
    ):
        raise RuntimeError("initial quadrature collar does not cover t+-3")
    diagnostics.update(
        {
            "left_t_collar": (left_t - outer_left_t).str(40).replace("e", "E"),
            "right_t_collar": (outer_right_t - right_t).str(40).replace(
                "e", "E"
            ),
            "extension_envelope": envelope,
        }
    )
    result = evaluate_nested_curvature_from_h_cover(
        MODE_TWO,
        COLLAR_END,
        derivatives,
        cover_diagnostics=diagnostics,
    )
    if result.get("passed") is not True:
        raise RuntimeError(f"initial quadrature collar failed: {result}")
    return result


def summarize_ray(records: list[dict]) -> dict:
    if len(records) != len(ray_tasks()):
        raise RuntimeError("finite ray cache is incomplete")
    for previous, current in zip(records, records[1:]):
        if previous["mode_right"] != current["mode_left"]:
            raise RuntimeError("finite ray cache has a mode gap")
    if records[0]["mode_left"] != str(COLLAR_END):
        raise RuntimeError("finite ray cache starts at the wrong mode")
    if records[-1]["mode_right"] != "20":
        raise RuntimeError("finite ray cache does not reach u=20")
    largest_scaled = max(
        range(len(records)),
        key=lambda index: Decimal(records[index]["scaled_curvature_upper"]),
    )
    weakest_margin = min(
        range(len(records)),
        key=lambda index: Decimal(records[index]["margin_lower"]),
    )
    weakest_J = min(
        range(len(records)), key=lambda index: Decimal(records[index]["J_lower"])
    )
    weakest_R = min(
        range(len(records)), key=lambda index: Decimal(records[index]["R_lower"])
    )
    selected_indices = sorted(
        {
            0,
            len(records) // 4,
            len(records) // 2,
            3 * len(records) // 4,
            len(records) - 1,
            largest_scaled,
            weakest_margin,
            weakest_J,
            weakest_R,
        }
    )

    def selected(index: int) -> dict:
        row = records[index]
        return {
            key: row[key]
            for key in (
                "index",
                "mode_left",
                "mode_right",
                "central_t_ball",
                "J_lower",
                "R_lower",
                "scaled_curvature_upper",
                "margin_lower",
            )
        }

    return {
        "block_count": len(records),
        "mode_range": [records[0]["mode_left"], records[-1]["mode_right"]],
        "largest_scaled_index": largest_scaled,
        "largest_scaled_curvature_upper": records[largest_scaled][
            "scaled_curvature_upper"
        ],
        "weakest_margin_index": weakest_margin,
        "weakest_margin_lower": records[weakest_margin]["margin_lower"],
        "weakest_J_index": weakest_J,
        "weakest_J_lower": records[weakest_J]["J_lower"],
        "weakest_R_index": weakest_R,
        "weakest_R_lower": records[weakest_R]["R_lower"],
        "selected": [selected(index) for index in selected_indices],
        "all_blocks_passed": True,
    }


def build_artifact(
    extension_cache: Path = DEFAULT_EXTENSION_CACHE,
    ray_cache: Path = DEFAULT_RAY_CACHE,
    *,
    workers: int = DEFAULT_WORKERS,
    overwrite: bool = False,
) -> dict:
    contract = source_contract()
    envelope = extension_envelope()
    extension = build_cache(
        extension_cache,
        extension_tasks(),
        extension_task,
        workers=workers,
        overwrite=overwrite,
    )
    collar = certify_initial_collar(extension, envelope)
    ray = build_cache(
        ray_cache,
        ray_tasks(),
        ray_task,
        workers=workers,
        overwrite=overwrite,
    )
    ray_summary = summarize_ray(ray)
    overall_scaled = max(
        Decimal(collar["scaled_curvature_upper"]),
        Decimal(ray_summary["largest_scaled_curvature_upper"]),
    )
    if overall_scaled >= CURVATURE_CONSTANT:
        raise RuntimeError("finite ray scaled curvature target failed")

    rows = [
        RayRow(
            id="co5ncfrc_01_extension_envelope",
            role="interval_analytic_lemma",
            readiness="ready_to_apply",
            claim="The existing eighth-derivative remainder envelope remains valid across the short quadrature extension above u=2.",
            formula="abs(V^(8)/V''^4)<1/50000 on 2<=u<=2.002",
            proof_boundary="Direct outward-rounded potential-jet interval.",
            diagnostics=envelope,
        ),
        RayRow(
            id="co5ncfrc_02_initial_collar",
            role="interval_theorem",
            readiness="ready_to_apply",
            claim="A compact quadrature collar closes the handoff immediately above mode two.",
            formula="q_1''(t)<=60/t^2 for 2<=u<=2.001",
            proof_boundary="Hashed base tiles plus 100 reproducible extension tiles.",
            diagnostics={
                "scaled_curvature_upper": collar["scaled_curvature_upper"],
                "margin_lower": collar["margin_lower"],
                "J_lower": collar["J_lower"],
                "R_lower": collar["R_lower"],
            },
        ),
        RayRow(
            id="co5ncfrc_03_exact_corridor_input",
            role="exact_theorem_input",
            readiness="ready_to_apply",
            claim="All seven exact cumulant corridors hold throughout every finite-ray derivative box.",
            formula=contract["exact_corridor_formula"],
            proof_boundary="Previously proved exact first-summand theorem.",
        ),
        RayRow(
            id="co5ncfrc_04_finite_ray",
            role="interval_theorem",
            readiness="ready_to_apply",
            claim="The nested first-summand curvature ceiling holds on the full finite mode ray through u=20.",
            formula="q_1''(t)<=60/t^2 for every 2<=u<=20",
            proof_boundary="Continuous interval cover; no point sampling is promoted.",
            diagnostics={
                "collar_blocks": 1,
                "exact_corridor_blocks": ray_summary["block_count"],
                "largest_scaled_curvature_upper": str(overall_scaled),
            },
        ),
        RayRow(
            id="co5ncfrc_05_open_asymptotic_ray",
            role="open_theorem_target",
            readiness="not_ready_to_apply",
            claim="Prove the same nested curvature ceiling on the remaining asymptotic ray.",
            formula="prove q_1''(t)<=60/t^2 for every mode u>=20",
            proof_boundary="Open asymptotic ray only.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate",
        "date": "2026-07-13",
        "status": (
            "rigorous first-summand nested-curvature theorem on 2<=u<=20 "
            "with one open asymptotic ray"
        ),
        "proof_boundary": (
            "This artifact proves q_1''(t)<=60/t^2 only on 2<=u<=20. "
            "It does not prove u>=20, complete order-five entry, PF-infinity, "
            "RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_compound_order5_nested_curvature_compact_certificate.md",
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.md",
            "outputs/formal_core.md",
        ],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.py"
        ),
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "workers": workers,
            "corridor_pad": str(CORRIDOR_PAD),
            "collar_t": COLLAR_T,
            "extension_mode_range": [str(EXTENSION_START), str(EXTENSION_END)],
            "initial_collar": [str(MODE_TWO), str(COLLAR_END)],
            "finite_ray": [str(COLLAR_END), "20"],
            "segments": [
                {
                    "left": str(segment.left),
                    "right": str(segment.right),
                    "width": str(segment.width),
                }
                for segment in SEGMENTS
            ],
        },
        "source_contract": contract,
        "extension": {
            "cache": extension_cache.relative_to(REPO_ROOT).as_posix(),
            "sha256": sha256(extension_cache),
            "rows": len(extension),
            "envelope": envelope,
        },
        "initial_collar": collar,
        "finite_ray": {
            "cache": ray_cache.relative_to(REPO_ROOT).as_posix(),
            "sha256": sha256(ray_cache),
            **ray_summary,
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": sum(row.readiness == "ready_to_apply" for row in rows),
            "open_rows": sum(
                row.readiness == "not_ready_to_apply" for row in rows
            ),
            "extension_tiles": len(extension),
            "initial_collar_blocks": 1,
            "exact_corridor_blocks": ray_summary["block_count"],
            "finite_ray_theorems": 1,
            "open_asymptotic_rays": 1,
            "largest_scaled_curvature_upper": str(overall_scaled),
        },
    }


def write_note(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    ray = artifact["finite_ray"]
    collar = artifact["initial_collar"]
    lines = [
        "# Jensen-Window PF Compound Order-Five Nested Curvature Finite-Ray Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous first-summand nested-curvature theorem on",
        "`2<=u<=20` with one open asymptotic ray. This is not a proof of",
        "complete order-five entry, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.py",
        "```",
        "",
        "## Mode-Two Collar",
        "",
        "The compact cache reaches slightly above `u=2`. One hundred additional",
        "quadrature tiles extend the common `t+-3` derivative collar through",
        "`u=2.001`. Direct potential-jet interval arithmetic proves",
        "",
        "```text",
        "abs(V^(8)/V''^4)<1/50000 on 2<=u<=2.002.",
        "```",
        "",
        "Thus the same paired Simpson remainder theorem applies. The collar gives",
        "",
        "```text",
        "q_1''(t)<=60/t^2 for every 2<=u<=2.001",
        f"largest scaled upper={collar['scaled_curvature_upper']}",
        f"margin lower={collar['margin_lower']}",
        "```",
        "",
        "The extension cache has SHA-256",
        "",
        "```text",
        artifact["extension"]["sha256"],
        "```",
        "",
        "## Exact-Corridor Cover",
        "",
        "The proved seven exact cumulant corridors generate common",
        "`H^(2),...,H^(8)` boxes on each central mode block. Every block has an",
        "explicit outer pad exceeding three units in the continuous `t` variable.",
        "The nested stable interval core then proves `J_1>0`, `R_1>0`, and",
        "the curvature ceiling on all",
        f"`{summary['exact_corridor_blocks']}` blocks from `u=2.001` to `u=20`.",
        "",
        "```text",
        f"largest t^2*q_1'' upper={summary['largest_scaled_curvature_upper']}",
        f"weakest J lower={ray['weakest_J_lower']}",
        f"weakest R lower={ray['weakest_R_lower']}",
        f"weakest absolute margin={ray['weakest_margin_lower']}",
        "```",
        "",
        "Combining the collar and exact-corridor blocks proves",
        "",
        "```text",
        "q_1''(t)<=60/t^2 for every mode 2<=u<=20.",
        "```",
        "",
        "## Remaining Ray",
        "",
        "The sole remaining part of the continuous first-summand theorem is",
        "",
        "```text",
        "q_1''(t)<=60/t^2 for every mode u>=20.",
        "```",
        "",
        "The existing normalized-H boxes `0<x_r<=1` for `2<=r<=8`, together",
        "with `x_2>=97/100`, `x_3>=24/25`, and `1/t<=10^-30`, leave a very",
        "large interval reserve for that asymptotic composition.",
        "",
        "| mode interval | t ball | t^2 q_1'' upper | margin lower | J lower | R lower |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in ray["selected"]:
        lines.append(
            f"| `{row['mode_left']}..{row['mode_right']}` | "
            f"`{row['central_t_ball']}` | `{row['scaled_curvature_upper']}` | "
            f"`{row['margin_lower']}` | `{row['J_lower']}` | `{row['R_lower']}` |"
        )
    lines.extend(
        [
            "",
            "```text",
            "outputs/jensen_window_pf_compound_order5_nested_curvature_compact_certificate.md",
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.md",
            "outputs/signed_hankel_jensen_bridge_target.md",
            "outputs/formal_core.md",
            "```",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extension-cache", type=Path, default=DEFAULT_EXTENSION_CACHE)
    parser.add_argument("--ray-cache", type=Path, default=DEFAULT_RAY_CACHE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact(
        args.extension_cache,
        args.ray_cache,
        workers=args.workers,
        overwrite=args.overwrite,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as handle:
        json.dump(artifact, handle, indent=2, sort_keys=True, default=str)
        handle.write("\n")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-five nested curvature finite-ray certificate: "
        f"{summary['extension_tiles']} extension tiles, "
        f"{summary['exact_corridor_blocks']} exact-corridor blocks, "
        f"{summary['finite_ray_theorems']} finite-ray theorem, "
        f"{summary['open_asymptotic_rays']} open asymptotic ray"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
