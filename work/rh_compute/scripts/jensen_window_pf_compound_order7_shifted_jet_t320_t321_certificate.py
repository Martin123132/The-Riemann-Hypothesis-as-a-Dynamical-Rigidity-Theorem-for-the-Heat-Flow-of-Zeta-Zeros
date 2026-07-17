#!/usr/bin/env python3
"""Certify r_1''(t)<=600/t^2 continuously on 320<=t<=321."""

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


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as compact  # noqa: E402
import jensen_window_pf_compound_order7_localized_tent_endpoint_certificate as endpoint_cache  # noqa: E402
from jensen_window_pf_compound_order7_shifted_jet_continuation_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    POINT_ORDER,
    PRECISION_BITS,
    QUADRATURE_PANELS,
    QUADRATURE_TAYLOR_ORDER,
    QUADRATURE_WINDOW_Y,
    REMAINDER_ORDER,
    continuation_row,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_shifted_jet_t320_t321_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order7_shifted_jet_t320_t321_certificate.md"
)
DEFAULT_WORKERS = max(1, min(12, (os.cpu_count() or 4) - 1))
BLOCK_COUNT = 12
START_T = Fraction(320)
END_T = Fraction(321)


@dataclass(frozen=True)
class CertificateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


WORKER_H_BOUNDS: dict[int, flint.arb] | None = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_h_bounds() -> dict[int, flint.arb]:
    flint.ctx.prec = PRECISION_BITS
    expected = endpoint_cache.tasks()
    records = endpoint_cache.load_cache(endpoint_cache.DEFAULT_CACHE, expected)
    if len(records) != len(expected):
        raise RuntimeError("complete the H2-H22 endpoint cache first")
    rows = endpoint_cache.arb_rows(records)
    return {
        order: compact.hull([row["H"][order] for row in rows])
        for order in range(2, 22)
    }


def initialize_worker() -> None:
    global WORKER_H_BOUNDS
    flint.ctx.prec = PRECISION_BITS
    WORKER_H_BOUNDS = load_h_bounds()


def deterministic_tasks() -> list[tuple[int, Fraction, Fraction]]:
    return [
        (
            index,
            START_T + Fraction(index, BLOCK_COUNT),
            START_T + Fraction(index + 1, BLOCK_COUNT),
        )
        for index in range(BLOCK_COUNT)
    ]


def block_task(task: tuple[int, Fraction, Fraction]) -> dict:
    index, anchor, right = task
    if WORKER_H_BOUNDS is None:
        h_bounds = load_h_bounds()
    else:
        h_bounds = WORKER_H_BOUNDS
    row = continuation_row(anchor, right, h_bounds)
    row["index"] = index
    return row


def source_contract() -> dict:
    expected = endpoint_cache.tasks()
    records = endpoint_cache.load_cache(endpoint_cache.DEFAULT_CACHE, expected)
    if len(records) != len(expected):
        raise RuntimeError("H2-H22 source cache is incomplete")
    return {
        "h_cache": endpoint_cache.DEFAULT_CACHE.relative_to(REPO_ROOT).as_posix(),
        "h_cache_sha256": sha256(endpoint_cache.DEFAULT_CACHE),
        "h_cache_rows": len(records),
        "h_derivative_orders": [2, 22],
        "h_mode_range": [str(endpoint_cache.MODE_START), str(endpoint_cache.MODE_END)],
        "tail_convexity_source": (
            "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.md"
        ),
        "tail_convexity_contract": (
            "V''>0 for u>=1/100 and V'<0 for 0<u<=1/100"
        ),
    }


def build_artifact(blocks: list[dict]) -> dict:
    if len(blocks) != BLOCK_COUNT:
        raise RuntimeError("continuation block count changed")
    for expected, block in zip(deterministic_tasks(), blocks):
        index, anchor, right = expected
        if (
            block.get("index") != index
            or block.get("anchor") != str(anchor)
            or block.get("right") != str(right)
            or block.get("passed") is not True
        ):
            raise RuntimeError(f"continuation block mismatch at index {index}")
    for previous, current in zip(blocks, blocks[1:]):
        if previous["right"] != current["anchor"]:
            raise RuntimeError("continuation blocks have a gap")

    largest = max(blocks, key=lambda row: Decimal(row["scaled_curvature_upper"]))
    weakest_margin = min(
        blocks, key=lambda row: Decimal(row["curvature_margin_lower"])
    )
    weakest_t = min(
        blocks,
        key=lambda row: Decimal(row["gaps"]["T"]["continued_lower"]),
    )
    maximum_panel_error = max(
        (row["quadrature"]["maximum_panel_error_upper"] for row in blocks),
        key=Decimal,
    )
    maximum_tail = max(
        (row["quadrature"]["maximum_tail_moment_upper"] for row in blocks),
        key=Decimal,
    )
    rows = [
        CertificateRow(
            "co7sj_01_exact_potential_quadrature",
            "interval_input",
            "ready_to_apply",
            "Exact-potential panel Taylor quadrature rigorously encloses H jets at every shifted anchor.",
            "H^(0),...,H^(10) at t_i+j, -5<=j<=5",
            "First Newman summand at lambda=-100 on the recorded point set.",
            {
                "panels": QUADRATURE_PANELS,
                "taylor_order": QUADRATURE_TAYLOR_ORDER,
                "window_y": QUADRATURE_WINDOW_Y,
                "maximum_panel_error_upper": maximum_panel_error,
                "maximum_tail_moment_upper": maximum_tail,
            },
        ),
        CertificateRow(
            "co7sj_02_shifted_stable_jets",
            "exact_interval_algebra",
            "ready_to_apply",
            "Exact centered-difference series algebra gives point jets through order ten for every stable layer.",
            "B,J,R,S,T,r jets through derivative order 10",
            "Outward-rounded point theorem only before continuation.",
        ),
        CertificateRow(
            "co7sj_03_bootstrap_majorant",
            "analytic_interval_bridge",
            "ready_to_apply",
            "Complete-monotonicity majorants and order-eleven Taylor remainders preserve every half-gap floor on each block.",
            "C(t)>=C(t_i)/2 for C in {B,J,R,S,T}",
            "Uses the hashed H2-H22 collar and a first-exit bootstrap argument.",
            {"blocks": BLOCK_COUNT, "remainder_order": REMAINDER_ORDER},
        ),
        CertificateRow(
            "co7sj_04_continuous_curvature",
            "interval_theorem",
            "ready_to_apply",
            "The fourth-nested first-summand curvature ceiling holds continuously from t=320 through t=321.",
            "r_1''(t)<=600/t^2 for 320<=t<=321",
            "Finite twelve-block endpoint theorem only.",
            {
                "largest_scaled_upper": largest["scaled_curvature_upper"],
                "weakest_margin_lower": weakest_margin["curvature_margin_lower"],
                "weakest_T_lower": weakest_t["gaps"]["T"]["continued_lower"],
            },
        ),
        CertificateRow(
            "co7sj_05_compact_handoff",
            "open_handoff",
            "not_ready_to_apply",
            "Continue the same shifted-jet or coarser interval theorem beyond t=321.",
            "r_1''(t)<=600/t^2 for t>=321",
            "The remaining compact, finite-ray, and asymptotic-ray ranges are separate tasks.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order7_shifted_jet_t320_t321_certificate",
        "date": "2026-07-13",
        "status": "rigorous continuous fourth-nested first-summand curvature theorem on 320<=t<=321",
        "proof_boundary": (
            "This artifact proves r_1''(t)<=600/t^2 only on 320<=t<=321. "
            "It does not prove the remaining t>=321 range, full-kernel order "
            "seven, PF-infinity, RH, or Lambda<=0."
        ),
        "source_contract": source_contract(),
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "point_order": POINT_ORDER,
            "remainder_order": REMAINDER_ORDER,
            "curvature_constant": CURVATURE_CONSTANT,
            "start_t": str(START_T),
            "end_t": str(END_T),
            "block_count": BLOCK_COUNT,
            "block_width": str(Fraction(1, BLOCK_COUNT)),
        },
        "blocks": blocks,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": 4,
            "open_rows": 1,
            "blocks": len(blocks),
            "all_blocks_passed": True,
            "continuous_curvature_theorems": 1,
            "largest_scaled_curvature_upper": largest["scaled_curvature_upper"],
            "weakest_curvature_margin_lower": weakest_margin[
                "curvature_margin_lower"
            ],
            "weakest_T_lower": weakest_t["gaps"]["T"]["continued_lower"],
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order7_shifted_jet_t320_t321_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order7_shifted_jet_t320_t321_certificate.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    source = artifact["source_contract"]
    lines = [
        "# Jensen-Window PF Compound Order-Seven Shifted-Jet Endpoint Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous continuous fourth-nested first-summand curvature",
        "theorem on `320<=t<=321`. This is not a proof of the remaining",
        "`t>=321` range, full-kernel order seven, PF-infinity, RH, or",
        "`Lambda <= 0`.",
        "",
        "## Method",
        "",
        "At each rational anchor `320+i/12`, exact-potential degree-30 panel",
        "Taylor quadrature encloses `H^(0),...,H^(10)` at all eleven shifted",
        "arguments. Exact centered-difference series algebra then produces the",
        "point jets for `B,J,R,S,T,r` through order ten.",
        "",
        "The hashed `H^(2),...,H^(22)` collar supplies order-eleven absolute",
        "majorants. Complete monotonicity of derivatives of",
        "`log(1-exp(-x))` and a first-exit argument preserve half of every",
        "point gap throughout each block.",
        "",
        "```text",
        source["h_cache"],
        f"sha256={source['h_cache_sha256']}",
        f"rows={source['h_cache_rows']}",
        "```",
        "",
        "## Theorem",
        "",
        "All twelve contiguous blocks pass:",
        "",
        "```text",
        "r_1''(t)<=600/t^2 for every real 320<=t<=321",
        f"largest scaled upper={summary['largest_scaled_curvature_upper']}",
        f"weakest curvature margin={summary['weakest_curvature_margin_lower']}",
        f"weakest T lower={summary['weakest_T_lower']}",
        "```",
        "",
        "## Remaining Range",
        "",
        "The next certificate begins at `t=321`. The compact continuation,",
        "finite saddle ray, and asymptotic ray remain separate proof tasks.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    args = parser.parse_args()

    tasks = deterministic_tasks()
    workers = max(1, args.workers)
    if workers == 1:
        initialize_worker()
        blocks = [block_task(task) for task in tasks]
    else:
        with ProcessPoolExecutor(
            max_workers=workers,
            initializer=initialize_worker,
        ) as executor:
            blocks = list(executor.map(block_task, tasks, chunksize=1))
    artifact = build_artifact(blocks)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-seven shifted-jet t=320..321 certificate: "
        f"{summary['blocks']} blocks, scaled upper "
        f"{summary['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
