#!/usr/bin/env python3
"""Replay retained-decay endpoint majorants on recursive t=1e10 blocks 20-28."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import jensen_window_pf_newman_c1_hardy_block20_retained_decay_endpoint_majorant_gate as retained
import jensen_window_pf_newman_c1_hardy_t1e10_crossblock_structural_atlas_gate as atlas
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import ctx


TELEMETRY = atlas.TELEMETRY
ATLAS = atlas.RESULT
BLOCK20 = retained.RESULT
CACHE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_retained_endpoint_gate_cache.jsonl"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_retained_endpoint_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_retained_endpoint_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_t1e10_crossblock_retained_endpoint_gate.py"
PRECISIONS = (60, 90)
EXPECTED_NEW_CALLS = 1040
EXPECTED_TOTAL_CALLS = 1414


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


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


def load_recursive() -> dict[int, dict[str, Any]]:
    chains: dict[int, dict[str, Any]] = {}
    for line in TELEMETRY.read_text(encoding="utf-8").splitlines():
        record = json.loads(line)
        chain = int(record.get("chain", 0))
        if record["type"] == "chain":
            chains[chain] = {"header": record, "levels": {}}
        elif record["type"] == "level":
            chains[chain]["levels"][int(record["level"])] = record
    return {
        chain: data
        for chain, data in chains.items()
        if int(data["header"]["mit"]) == 2 and 20 <= int(data["header"]["block"]) <= 28
    }


def cache_fingerprint() -> str:
    material = {
        "telemetry": file_hash(TELEMETRY),
        "atlas": file_hash(ATLAS),
        "block20": file_hash(BLOCK20),
        "retained_builder": file_hash(Path(retained.__file__).resolve()),
        "builder": file_hash(BUILDER),
        "checker": file_hash(CHECKER),
        "precisions": list(PRECISIONS),
    }
    return hashlib.sha256(json.dumps(material, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def load_cache(fingerprint: str) -> dict[int, dict[str, Any]]:
    if not CACHE.exists():
        CACHE.write_text(
            json.dumps({"kind": "t1e10_crossblock_retained_endpoint_cache", "fingerprint": fingerprint}, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return {}
    lines = [line for line in CACHE.read_text(encoding="utf-8").splitlines() if line.strip()]
    require(lines, "cross-block endpoint cache empty")
    header = json.loads(lines[0])
    require(header.get("fingerprint") == fingerprint, "cross-block endpoint cache fingerprint drift")
    rows: dict[int, dict[str, Any]] = {}
    for line in lines[1:]:
        row = json.loads(line)
        chain = int(row["chain"])
        require(chain not in rows, f"duplicate cached chain {chain}")
        rows[chain] = row
    return rows


def append_cache(row: dict[str, Any]) -> None:
    with CACHE.open("a", encoding="utf-8", buffering=1) as handle:
        handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def evaluate_row(chain: int, data: dict[str, Any]) -> dict[str, Any]:
    low = retained.evaluate_chain(data, PRECISIONS[0])
    high = retained.evaluate_chain(data, PRECISIONS[1])
    for name in ("upper_family", "lower_family", "complete"):
        require(low[name].overlaps(high[name]), f"chain {chain} retained endpoint {name} precision drift")
    require(low["selector"] == high["selector"], f"chain {chain} selector precision drift")
    header = data["header"]
    parent = data["levels"][1]
    return {
        "chain": chain,
        "block": int(header["block"]),
        "sum_index": int(header["sum_index"]),
        "branch": int(header["branch"]),
        "parent_length": int(parent["length"]),
        "child_length": int(data["levels"][2]["length"]),
        "selector": high["selector"],
        "delta": decimal(high["delta"]),
        "eta": decimal(high["eta"]),
        "upper_family_majorant_lower": decimal(retained.lower(high["upper_family"])),
        "upper_family_majorant_upper": decimal(retained.upper(high["upper_family"])),
        "lower_family_majorant_lower": decimal(retained.lower(high["lower_family"])),
        "lower_family_majorant_upper": decimal(retained.upper(high["lower_family"])),
        "complete_majorant_lower": decimal(retained.lower(high["complete"])),
        "complete_majorant_upper": decimal(retained.upper(high["complete"])),
        "precision_overlap": True,
    }


def block20_rows() -> dict[int, dict[str, Any]]:
    artifact = json.loads(BLOCK20.read_text(encoding="utf-8"))
    atlas_artifact = json.loads(ATLAS.read_text(encoding="utf-8"))
    child_lengths = {
        int(row["chain"]): int(row["child_length"])
        for row in atlas_artifact["rows"]
        if int(row["block"]) == 20 and int(row["mit"]) == 2
    }
    rows: dict[int, dict[str, Any]] = {}
    for source in artifact["rows"]:
        chain = int(source["chain"])
        rows[chain] = {
            "chain": chain,
            "block": 20,
            "sum_index": int(source["sum_index"]),
            "branch": int(source["branch"]),
            "parent_length": 104,
            "child_length": child_lengths[chain],
            "selector": source["selector"],
            "delta": source["delta"],
            "eta": source["eta"],
            "upper_family_majorant_lower": source["upper_family_majorant_lower"],
            "upper_family_majorant_upper": source["upper_family_majorant_upper"],
            "lower_family_majorant_lower": source["lower_family_majorant_lower"],
            "lower_family_majorant_upper": source["lower_family_majorant_upper"],
            "complete_majorant_lower": source["complete_majorant_lower"],
            "complete_majorant_upper": source["complete_majorant_upper"],
            "precision_overlap": True,
            "reused_admitted_block20_row": True,
        }
    require(len(rows) == 374, "admitted block-20 retained roster drift")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime-limit-seconds", type=float, default=10800.0)
    args = parser.parse_args()
    for path in (TELEMETRY, ATLAS, BLOCK20, CHECKER):
        require(path.is_file(), f"missing cross-block endpoint dependency: {path}")
    priority = set_low_priority()
    ctx.threads = 1
    started = time.monotonic()
    recursive = load_recursive()
    require(len(recursive) == EXPECTED_TOTAL_CALLS, "cross-block recursive roster drift")
    new_recursive = {
        chain: data for chain, data in recursive.items() if 21 <= int(data["header"]["block"]) <= 28
    }
    require(len(new_recursive) == EXPECTED_NEW_CALLS, "new recursive roster drift")
    fingerprint = cache_fingerprint()
    cached = load_cache(fingerprint)
    require(set(cached).issubset(new_recursive), "cache contains an unexpected chain")
    for chain in sorted(new_recursive):
        if chain in cached:
            continue
        row = evaluate_row(chain, new_recursive[chain])
        append_cache(row)
        cached[chain] = row
        if time.monotonic() - started >= args.runtime_limit_seconds:
            print(f"parked t=1e10 cross-block retained endpoint replay after {len(cached)}/{EXPECTED_NEW_CALLS} new calls")
            return 2

    require(len(cached) == EXPECTED_NEW_CALLS, "cross-block retained endpoint cache incomplete")
    all_rows = block20_rows()
    all_rows.update(cached)
    require(len(all_rows) == EXPECTED_TOTAL_CALLS, "cross-block retained endpoint combined roster drift")
    rows = [all_rows[chain] for chain in sorted(all_rows)]

    block_rows: list[dict[str, Any]] = []
    for block in range(20, 29):
        selected = [row for row in rows if int(row["block"]) == block]
        selectors = Counter(row["selector"] for row in selected)
        maximum = max((Fraction(row["complete_majorant_upper"]), int(row["chain"])) for row in selected)
        minimum = min((Fraction(row["complete_majorant_lower"]), int(row["chain"])) for row in selected)
        ordered = sorted(Fraction(row["complete_majorant_upper"]) for row in selected)
        block_rows.append(
            {
                "block": block,
                "recursive_call_count": len(selected),
                "selector_histogram": dict(sorted(selectors.items())),
                "minimum_complete_majorant_lower": decimal(minimum[0]),
                "minimum_complete_majorant_witness": minimum[1],
                "median_complete_majorant_upper": decimal(ordered[len(ordered) // 2]),
                "maximum_complete_majorant_upper": decimal(maximum[0]),
                "maximum_complete_majorant_witness": maximum[1],
            }
        )

    maximum = max((Fraction(row["complete_majorant_upper"]), int(row["chain"])) for row in rows)
    aggregate = {
        "recursive_call_count": len(rows),
        "newly_evaluated_call_count": len(cached),
        "reused_block20_call_count": len(rows) - len(cached),
        "precision_overlap_count": sum(bool(row["precision_overlap"]) for row in rows),
        "block_count": len(block_rows),
        "first_block": 20,
        "last_recursive_block": 28,
        "w3_count": sum(row["selector"] == "W3" for row in rows),
        "w4_count": sum(row["selector"] == "W4" for row in rows),
        "maximum_complete_majorant_upper": decimal(maximum[0]),
        "maximum_complete_majorant_witness": maximum[1],
        "maximum_complete_majorant_block": next(int(row["block"]) for row in rows if int(row["chain"]) == maximum[1]),
        "block20_maximum_complete_majorant_upper": max(
            (row["maximum_complete_majorant_upper"] for row in block_rows if int(row["block"]) == 20), key=Fraction
        ),
        "later_to_block20_maximum_ratio": decimal(
            max(Fraction(row["maximum_complete_majorant_upper"]) for row in block_rows[1:])
            / Fraction(block_rows[0]["maximum_complete_majorant_upper"])
        ),
    }
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_crossblock_retained_endpoint_gate",
        "status": "rigorous_retained_decay_endpoint_majorants_for_all_1414_recursive_calls_in_blocks_20_28",
        "scope": "Exact-point retained-decay endpoint majorants for every recursive cubic call in t=1e10 blocks 20 through 28.",
        "precisions_decimal_digits": list(PRECISIONS),
        "blocks": block_rows,
        "rows": rows,
        "aggregate": aggregate,
        "runtime": {
            "priority": priority,
            "active_workers": 1,
            "elapsed_seconds": time.monotonic() - started,
            "cache_rows_reused": EXPECTED_NEW_CALLS,
        },
        "dependencies": {
            "crossblock_atlas": {"path": relative(ATLAS), "sha256": file_hash(ATLAS)},
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "block20_retained_endpoint": {"path": relative(BLOCK20), "sha256": file_hash(BLOCK20)},
            "retained_endpoint_builder": {"path": relative(Path(retained.__file__).resolve()), "sha256": file_hash(Path(retained.__file__).resolve())},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "cache": {"path": relative(CACHE), "sha256": file_hash(CACHE)},
        },
        "next_obligation": (
            "Capture exact outer accumulation weights for blocks 21-35, transport these recursive majorants, and "
            "certify all 5,370 direct kernels before deriving a height-uniform block transition."
        ),
        "proof_boundary": (
            "Rigorous finite exact-point local endpoint majorants at t=1e10 only. No coefficient-cell transport is "
            "claimed outside block 20, no cross-block output weight accumulation or later direct-kernel error is yet "
            "included, and no height-uniform, outer Hardy, Lambda<=0, PF-infinity, RH, or prize-level conclusion follows."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built t=1e10 cross-block retained endpoint gate: "
        f"{aggregate['recursive_call_count']} calls, blocks 20-28, max {aggregate['maximum_complete_majorant_upper']}"
    )
    return 0


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    lines = [
        "# Hardy t=1e10 cross-block retained endpoint gate",
        "",
        "Date: 2026-08-09",
        f"Status: {artifact['status']}; finite exact-point theorem only, not a proof of RH",
        "",
        "The retained common-ray theorem from block 20 is replayed at 60 and 90",
        "decimal digits on every recursive call in blocks 21 through 28. The 374",
        "admitted block-20 rows are reused, giving the complete recursive tail.",
        "",
        "```text",
        "block  calls  W3  W4  median upper  maximum upper",
    ]
    for row in artifact["blocks"]:
        histogram = row["selector_histogram"]
        lines.append(
            f"{row['block']:>5} {row['recursive_call_count']:>6} {histogram.get('W3', 0):>3} "
            f"{histogram.get('W4', 0):>3}  {row['median_complete_majorant_upper']}  "
            f"{row['maximum_complete_majorant_upper']}"
        )
    lines.extend(
        [
            "```",
            "",
            f"The roster-wide maximum is `{aggregate['maximum_complete_majorant_upper']}` in block "
            f"`{aggregate['maximum_complete_majorant_block']}`. The maximum later-block to block-20 ratio is "
            f"`{aggregate['later_to_block20_maximum_ratio']}`.",
            "",
            "These are local endpoint majorants. Cross-block output weights, coefficient",
            "cells outside block 20, later direct kernels, height uniformity, the outer",
            "Hardy representation, `Lambda<=0`, PF-infinity, RH, and a prize theorem",
            "remain open.",
            "",
        ]
    )
    return "\n".join(lines)


if __name__ == "__main__":
    raise SystemExit(main())
