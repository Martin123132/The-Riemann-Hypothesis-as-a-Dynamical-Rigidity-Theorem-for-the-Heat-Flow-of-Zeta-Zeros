#!/usr/bin/env python3
"""Compose the complete order-twelve sparse-H23 physical lower theorem."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS = REPO_ROOT / "work/rh_compute/results"
OUTPUTS = REPO_ROOT / "outputs"
for candidate in (SCRIPT_DIR, REPO_ROOT / "work/rh_compute/vendor"):
    if str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import check_jensen_window_pf_compound_order12_sparse_h23_lower_bridge_segments as checker  # noqa: E402
import jensen_window_pf_compound_order12_sparse_h23_lower_bridge_segments as source  # noqa: E402


DEFAULT_OUT = RESULTS / "jensen_window_pf_compound_order12_sparse_h23_lower_bridge_certificate.json"
DEFAULT_NOTE = OUTPUTS / "jensen_window_pf_compound_order12_sparse_h23_lower_bridge_certificate.md"
GENERATOR = "work/rh_compute/scripts/jensen_window_pf_compound_order12_sparse_h23_lower_bridge_certificate.py"
CHECKER = "work/rh_compute/scripts/check_jensen_window_pf_compound_order12_sparse_h23_lower_bridge_certificate.py"
THEOREM = "v_1''(t)<=8000/t^2 for every real 1503<=t<=5700"


@dataclass(frozen=True)
class CertificateRow:
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


def build_artifact() -> dict:
    tasks = source.deterministic_segments()
    issues = checker.validate_run_contract(source.DEFAULT_RUN_CONTRACT)
    validated = checker.validate_segment_cache(source.DEFAULT_SEGMENT_CACHE, tasks)
    segment_issues, segment_count, block_count, maximum, minimum_margin = validated
    issues.extend(segment_issues)
    if segment_count != len(tasks):
        issues.append(f"expected {len(tasks)} segments, found {segment_count}")
    if block_count != source.EXPECTED_QUARTER_BLOCKS:
        issues.append(
            f"expected {source.EXPECTED_QUARTER_BLOCKS} quarter blocks, "
            f"found {block_count}"
        )
    if tasks[0][1] != source.START_T or tasks[-1][2] != source.END_T:
        issues.append("deterministic segment endpoints changed")
    if issues:
        raise RuntimeError("; ".join(issues[:20]))

    maximum_text = str(maximum)
    minimum_text = str(minimum_margin)
    rows = [
        CertificateRow(
            "co12sh23lbc_01_contract",
            "immutable_source_contract",
            "ready_to_apply",
            "The complete lower cache is bound to the inherited exact sparse source, mixed H24 wall, and order-twelve Taylor implementation.",
            f"SHA256(cache)={sha256(source.DEFAULT_SEGMENT_CACHE)}",
            "Hash, source-manifest, and deterministic geometry checks only.",
        ),
        CertificateRow(
            "co12sh23lbc_02_complete_cover",
            "interval_theorem",
            "ready_to_apply",
            "The canonical segments form a gap-free cover of the order-twelve physical lower interval.",
            "[1503,5700]=union of 16788 certified quarter blocks",
            "First Newman summand at lambda=-100 on the displayed finite interval.",
            {"segments": segment_count, "quarter_blocks": block_count},
        ),
        CertificateRow(
            "co12sh23lbc_03_curvature",
            "interval_theorem",
            "ready_to_apply",
            "Every lower quarter block has positive margin below the order-twelve curvature cap.",
            THEOREM,
            "This does not cover the compact interval t>5700 or promote the result to the full Newman kernel.",
            {
                "largest_scaled_curvature_upper": maximum_text,
                "smallest_margin_lower": minimum_text,
            },
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order12_sparse_h23_lower_bridge_certificate",
        "date": "2026-07-22",
        "status": "rigorous order-twelve first-summand curvature theorem on 1503<=t<=5700",
        "theorem": THEOREM,
        "proof_boundary": (
            "This artifact proves only the complete finite first-summand physical lower bridge. "
            "The compact interval, global composition, full-kernel transfer, fixed order twelve, "
            "PF-infinity, Lambda<=0, and RH are not proved here."
        ),
        "source_contract": {
            "segment_cache": relative(source.DEFAULT_SEGMENT_CACHE),
            "segment_cache_sha256": sha256(source.DEFAULT_SEGMENT_CACHE),
            "run_contract": relative(source.DEFAULT_RUN_CONTRACT),
            "run_contract_sha256": sha256(source.DEFAULT_RUN_CONTRACT),
        },
        "summary": {
            "segments": segment_count,
            "required_segments": len(tasks),
            "quarter_blocks": block_count,
            "required_quarter_blocks": source.EXPECTED_QUARTER_BLOCKS,
            "covered_interval": [str(source.START_T), str(source.END_T)],
            "largest_scaled_curvature_upper": maximum_text,
            "smallest_margin_lower": minimum_text,
            "lower_bridge_theorems": 1,
            "global_first_summand_theorems": 0,
            "rh_claims": 0,
        },
        "rows": [asdict(row) for row in rows],
        "generator": GENERATOR,
        "checker": CHECKER,
    }


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT).as_posix()


def write_note(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    lines = [
        "# Order-Twelve Sparse-H23 Physical Lower Certificate",
        "",
        "Date: 2026-07-22",
        "",
        f"Status: **{artifact['status']}**. This is not a proof of RH or `Lambda <= 0`.",
        "",
        "## Theorem",
        "",
        f"`{artifact['theorem']}`.",
        "",
        "```text",
        f"segments={summary['segments']}/{summary['required_segments']}",
        f"quarter blocks={summary['quarter_blocks']}/{summary['required_quarter_blocks']}",
        f"largest scaled upper={summary['largest_scaled_curvature_upper']}",
        f"smallest margin={summary['smallest_margin_lower']}",
        f"segment cache SHA-256={artifact['source_contract']['segment_cache_sha256']}",
        "```",
        "",
        "## Boundary",
        "",
        artifact["proof_boundary"],
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    artifact = build_artifact()
    DEFAULT_OUT.parent.mkdir(parents=True, exist_ok=True)
    DEFAULT_OUT.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(DEFAULT_NOTE, artifact)
    print(
        "wrote complete order-twelve lower-bridge certificate: "
        f"{artifact['summary']['segments']} segments, "
        f"{artifact['summary']['quarter_blocks']} quarter blocks"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
