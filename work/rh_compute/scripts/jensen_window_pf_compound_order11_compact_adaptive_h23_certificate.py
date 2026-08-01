#!/usr/bin/env python3
"""Promote the complete adaptive-H23 compact cache to an order-eleven theorem."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import check_jensen_window_pf_compound_order11_compact_adaptive_h23_segments as segment_checker  # noqa: E402
import jensen_window_pf_compound_order11_compact_adaptive_h23_segments as source  # noqa: E402
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
    arb_upper_text,
)


DEFAULT_CACHE = source.DEFAULT_SEGMENT_CACHE
DEFAULT_RUN_CONTRACT = source.DEFAULT_RUN_CONTRACT
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_compact_adaptive_h23_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order11_compact_adaptive_h23_certificate.md"
)
MONOTONICITY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate.json"
)
POTENTIAL_SOURCE = (
    SCRIPT_DIR
    / "jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.py"
)
GENERATOR_PATH = (
    "work/rh_compute/scripts/"
    "jensen_window_pf_compound_order11_compact_adaptive_h23_certificate.py"
)
CHECKER_PATH = (
    "work/rh_compute/scripts/"
    "check_jensen_window_pf_compound_order11_compact_adaptive_h23_certificate.py"
)
SADDLE_RAY_START = Fraction(2001, 1000)
THEOREM = "y_1''(t)<=6000/t^2 for every real 5700<=t<=38020"
MONOTONICITY_FORMULA = "V''>0 for u>=1/100; V'<0 for 0<u<=1/100"


@dataclass(frozen=True)
class CompactRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT).as_posix()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def last_record(path: Path) -> dict:
    """Read only the final JSONL row so an incomplete cache fails cheaply."""
    if not path.exists() or path.stat().st_size == 0:
        raise RuntimeError(f"missing or empty compact segment cache: {path}")
    with path.open("rb") as handle:
        position = handle.seek(0, 2)
        data = b""
        while position > 0:
            width = min(1 << 20, position)
            position -= width
            handle.seek(position)
            data = handle.read(width) + data
            stripped = data.rstrip(b"\r\n")
            if b"\n" in stripped or position == 0:
                candidate = stripped.rsplit(b"\n", 1)[-1].strip()
                if not candidate:
                    break
                try:
                    record = json.loads(candidate)
                except json.JSONDecodeError as exc:
                    raise RuntimeError("invalid final compact JSONL row") from exc
                if not isinstance(record, dict):
                    raise RuntimeError("final compact JSONL row is not an object")
                return record
    raise RuntimeError("compact segment cache has no nonempty row")


def require_complete_tail(path: Path) -> None:
    tasks = source.deterministic_segments()
    final = last_record(path)
    expected_index, _, expected_right = tasks[-1]
    found_index = final.get("index")
    if (
        found_index != expected_index
        or final.get("segment_right") != str(expected_right)
        or final.get("passed") is not True
    ):
        found_count = found_index + 1 if isinstance(found_index, int) else "unknown"
        raise RuntimeError(
            "complete the compact segment cache before promotion: "
            f"expected {len(tasks)} segments through t={expected_right}, "
            f"found tail count {found_count}"
        )


def validate_complete_cache(cache: Path, run_contract: Path) -> dict:
    """Stream, validate, hash, and summarize the exact complete cache once."""
    contract_issues = segment_checker.validate_run_contract(run_contract)
    if contract_issues:
        raise RuntimeError("; ".join(contract_issues))
    require_complete_tail(cache)
    tasks = source.deterministic_segments()
    digest = hashlib.sha256()
    segment_count = 0
    block_count = 0
    largest: tuple[Decimal, int, int, str] | None = None
    smallest: tuple[Decimal, int, int, str] | None = None
    with cache.open("rb") as handle:
        for line_number, line in enumerate(handle, start=1):
            digest.update(line)
            if not line.strip():
                continue
            if segment_count >= len(tasks):
                raise RuntimeError("compact segment cache has too many rows")
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise RuntimeError(
                    f"invalid compact JSONL row {line_number}"
                ) from exc
            issues, scaled, margins = segment_checker.validate_segment_record(
                record,
                tasks[segment_count],
            )
            if issues:
                raise RuntimeError("; ".join(issues[:20]))
            blocks = record["blocks"]
            block_count += len(blocks)
            local_largest = max(scaled)
            local_largest_index = scaled.index(local_largest)
            largest_candidate = (
                local_largest,
                segment_count,
                local_largest_index,
                blocks[local_largest_index]["anchor"],
            )
            if largest is None or largest_candidate[0] > largest[0]:
                largest = largest_candidate
            local_smallest = min(margins)
            local_smallest_index = margins.index(local_smallest)
            smallest_candidate = (
                local_smallest,
                segment_count,
                local_smallest_index,
                blocks[local_smallest_index]["anchor"],
            )
            if smallest is None or smallest_candidate[0] < smallest[0]:
                smallest = smallest_candidate
            segment_count += 1
    if segment_count != source.EXPECTED_SEGMENTS:
        raise RuntimeError(
            f"expected {source.EXPECTED_SEGMENTS} segments, found {segment_count}"
        )
    if block_count != source.EXPECTED_QUARTER_BLOCKS:
        raise RuntimeError(
            "expected "
            f"{source.EXPECTED_QUARTER_BLOCKS} quarter blocks, found {block_count}"
        )
    if largest is None or smallest is None:
        raise RuntimeError("complete compact cache has no curvature summaries")
    if largest[0] >= Decimal(source.CURVATURE_CONSTANT) or smallest[0] <= 0:
        raise RuntimeError("complete compact cache does not fit the curvature cap")
    return {
        "segments": segment_count,
        "quarter_blocks": block_count,
        "sha256": digest.hexdigest(),
        "bytes": cache.stat().st_size,
        "largest_scaled_curvature_upper": str(largest[0]),
        "largest_location": {
            "segment": largest[1],
            "block": largest[2],
            "anchor": largest[3],
        },
        "smallest_margin_lower": str(smallest[0]),
        "smallest_location": {
            "segment": smallest[1],
            "block": smallest[2],
            "anchor": smallest[3],
        },
    }


def validate_monotonicity_source() -> dict:
    artifact = load_json(MONOTONICITY_SOURCE)
    rows = [
        row
        for row in artifact.get("rows", [])
        if row.get("id") == "co4ecetc_05_global_outward_monotonicity"
    ]
    if (
        artifact.get("kind")
        != "jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate"
        or len(rows) != 1
        or rows[0].get("readiness") != "ready_to_apply"
        or rows[0].get("formula") != MONOTONICITY_FORMULA
    ):
        raise RuntimeError("invalid saddle-monotonicity source")
    return {
        "path": relative(MONOTONICITY_SOURCE),
        "sha256": sha256(MONOTONICITY_SOURCE),
        "kind": artifact["kind"],
        "row": rows[0]["id"],
        "formula": rows[0]["formula"],
    }


def saddle_transition_upper() -> str:
    flint.ctx.prec = source.PRECISION_BITS
    transition = potential_jet_arb(arb_rational(SADDLE_RAY_START), 1)[1]
    if not bool(transition < arb_rational(source.END_T)):
        raise RuntimeError("u=2.001 saddle transition is not below t=38020")
    return arb_upper_text(transition)


def build_artifact(
    cache: Path = DEFAULT_CACHE,
    run_contract: Path = DEFAULT_RUN_CONTRACT,
) -> dict:
    summary = validate_complete_cache(cache, run_contract)
    monotonicity = validate_monotonicity_source()
    transition_upper = saddle_transition_upper()
    theorem_margin = Decimal(source.CURVATURE_CONSTANT) - Decimal(
        summary["largest_scaled_curvature_upper"]
    )
    rows = [
        CompactRow(
            "co11cah23c_01_contract",
            "immutable_source_contract",
            "ready_to_apply",
            "The complete segment cache is bound to the canonical adaptive-H23 run contract.",
            f"SHA256(cache)={summary['sha256']}",
            "Hash binding and deterministic live-source validation only.",
        ),
        CompactRow(
            "co11cah23c_02_partition",
            "interval_partition",
            "ready_to_apply",
            "The canonical quarter blocks form a gap-free partition of the compact real-t interval.",
            "[5700,38020]=union of 129280 certified quarter blocks",
            "First Newman summand at lambda=-100 on the displayed interval.",
            {
                "segments": summary["segments"],
                "quarter_blocks": summary["quarter_blocks"],
                "quarter_width": str(source.QUARTER_WIDTH),
            },
        ),
        CompactRow(
            "co11cah23c_03_curvature",
            "interval_theorem",
            "ready_to_apply",
            "Every compact quarter block has positive margin below the order-eleven curvature cap.",
            THEOREM,
            "Continuous compact first-summand theorem only.",
            {
                "largest_scaled_curvature_upper": summary[
                    "largest_scaled_curvature_upper"
                ],
                "largest_location": summary["largest_location"],
                "smallest_margin_lower": summary["smallest_margin_lower"],
                "smallest_location": summary["smallest_location"],
                "theorem_margin_lower": str(theorem_margin),
            },
        ),
        CompactRow(
            "co11cah23c_04_saddle_overlap",
            "exact_overlap",
            "ready_to_apply",
            "The compact endpoint lies strictly beyond the start of the complete monotone saddle ray.",
            "V'(2001/1000)<38020",
            "Joins first-summand curvature ranges only.",
            {
                "V_prime_2_001_upper": transition_upper,
                "monotonicity_formula": MONOTONICITY_FORMULA,
            },
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order11_compact_adaptive_h23_certificate",
        "date": "2026-07-22",
        "status": (
            "rigorous order-eleven first-summand curvature theorem on "
            "5700<=t<=38020"
        ),
        "proof_boundary": (
            "This artifact proves the compact first-summand theorem at "
            "lambda=-100 only. It does not by itself prove the global "
            "first-summand theorem, full-kernel transfer, PF-infinity, "
            "Lambda<=0, or RH."
        ),
        "theorem": THEOREM,
        "source_contract": {
            "run_contract": {
                "path": relative(run_contract),
                "sha256": sha256(run_contract),
            },
            "segment_cache": {
                "path": relative(cache),
                "sha256": summary["sha256"],
                "bytes": summary["bytes"],
                "segments": summary["segments"],
                "quarter_blocks": summary["quarter_blocks"],
            },
            "saddle_monotonicity": monotonicity,
            "potential_source": {
                "path": relative(POTENTIAL_SOURCE),
                "sha256": sha256(POTENTIAL_SOURCE),
            },
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            **summary,
            "theorem_margin_lower": str(theorem_margin),
            "saddle_transition_upper": transition_upper,
            "compact_first_summand_theorems": 1,
            "global_first_summand_theorems": 0,
            "full_kernel_theorems": 0,
            "rh_claims": 0,
        },
        "generator": GENERATOR_PATH,
        "checker": CHECKER_PATH,
    }


def write_note(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    lines = [
        "# Order-Eleven Compact Adaptive-H23 Certificate",
        "",
        "Date: 2026-07-22",
        "",
        f"Status: **{artifact['status']}**. This is not a proof of RH or `Lambda <= 0`.",
        "",
        "## Theorem",
        "",
        f"`{artifact['theorem']}`.",
        "",
        (
            f"All `{summary['segments']}` segments and "
            f"`{summary['quarter_blocks']}` quarter blocks pass."
        ),
        (
            "The largest scaled upper is "
            f"`{summary['largest_scaled_curvature_upper']}` and the smallest "
            f"margin is `{summary['smallest_margin_lower']}`."
        ),
        "",
        "## Saddle Overlap",
        "",
        (
            "The rigorous upper enclosure "
            f"`V'(2001/1000)<={summary['saddle_transition_upper']}<38020` "
            "joins the compact interval to the monotone saddle ray."
        ),
        "",
        "## Boundary",
        "",
        artifact["proof_boundary"],
        "",
        "This is not a proof of RH.",
        "",
        "## Reproduce",
        "",
        "```powershell",
        f"python {GENERATOR_PATH}",
        f"python {CHECKER_PATH}",
        "```",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument(
        "--run-contract",
        type=Path,
        default=DEFAULT_RUN_CONTRACT,
    )
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact(args.cache, args.run_contract)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    print(
        "wrote complete order-eleven compact adaptive-H23 theorem: "
        f"{artifact['summary']['segments']} segments, largest scaled upper "
        f"{artifact['summary']['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
