#!/usr/bin/env python3
"""Compose the complete order-ten localized first-summand lower bridge."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order10_localized_lower_bridge_segments as segments  # noqa: E402


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_localized_lower_bridge_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order10_localized_lower_bridge_certificate.md"
)
SEGMENT_CACHE = segments.DEFAULT_SEGMENT_CACHE
RUN_CONTRACT = segments.DEFAULT_RUN_CONTRACT
THEOREM_CURVATURE_CONSTANT = 4200


@dataclass(frozen=True)
class BridgeRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_and_compose() -> dict:
    tasks = segments.deterministic_segments()
    records = segments.load_segment_cache(SEGMENT_CACHE, tasks)
    if len(records) != len(tasks):
        raise RuntimeError(
            f"lower-bridge segment cache is incomplete: {len(records)}/{len(tasks)}"
        )
    expected_t = segments.START_T
    block_count = 0
    regime_counts = {"near": 0, "middle": 0, "far": 0}
    regime_maxima = {}
    global_maximum = None
    global_minimum_margin = None
    for record, task in zip(records, tasks):
        index, regime, segment_left, segment_right = task
        if Fraction(record["segment_left"]) != expected_t:
            raise RuntimeError(f"segment coverage gap before row {index}")
        if Fraction(record["segment_right"]) != segment_right:
            raise RuntimeError(f"segment right edge changed at row {index}")
        for block in record.get("blocks", []):
            left = Fraction(block["anchor"])
            right = Fraction(block["right"])
            expansion = Fraction(block["expansion_anchor"])
            if left != expected_t:
                raise RuntimeError(f"block coverage gap before t={left}")
            expected_width = (
                segments.QUARTER_WIDTH
                if regime in {"near", "middle"}
                else segments.CELL_WIDTH
            )
            if right - left != expected_width:
                raise RuntimeError(
                    f"bad {regime} block width on {left}..{right}"
                )
            if regime in {"near", "middle"}:
                if expansion not in {left, right}:
                    raise RuntimeError(
                        f"quarter block lacks an endpoint expansion on {left}..{right}"
                    )
            elif expansion != left:
                raise RuntimeError(f"far block expansion changed at t={left}")
            if block.get("regime") != regime or block.get("passed") is not True:
                raise RuntimeError(f"failed or mislabelled block at t={left}")
            if "-chi(W)*(W')^2" not in str(block.get("proof_formula", "")):
                raise RuntimeError(f"proof formula changed at t={left}")
            scaled = Decimal(block["scaled_curvature_upper"])
            margin = Decimal(block["curvature_margin_lower"])
            if not (Decimal(0) <= scaled < Decimal(segments.CURVATURE_CONSTANT)):
                raise RuntimeError(f"bad scaled curvature at t={left}: {scaled}")
            if margin <= 0:
                raise RuntimeError(f"nonpositive curvature margin at t={left}")
            summary = {
                "anchor": str(left),
                "right": str(right),
                "expansion_anchor": str(expansion),
                "regime": regime,
                "scaled_curvature_upper": block["scaled_curvature_upper"],
                "curvature_margin_lower": block["curvature_margin_lower"],
            }
            if (
                global_maximum is None
                or scaled > Decimal(global_maximum["scaled_curvature_upper"])
            ):
                global_maximum = summary
            if (
                global_minimum_margin is None
                or margin
                < Decimal(global_minimum_margin["curvature_margin_lower"])
            ):
                global_minimum_margin = summary
            old_regime = regime_maxima.get(regime)
            if old_regime is None or scaled > Decimal(
                old_regime["scaled_curvature_upper"]
            ):
                regime_maxima[regime] = summary
            block_count += 1
            regime_counts[regime] += 1
            expected_t = right
        if expected_t != segment_right:
            raise RuntimeError(f"segment {index} block coverage is incomplete")
    if expected_t != segments.END_T:
        raise RuntimeError(f"lower bridge ends at {expected_t}, not {segments.END_T}")
    expected_counts = {"near": 996, "middle": 1200, "far": 7800}
    if regime_counts != expected_counts or block_count != 9996:
        raise RuntimeError(
            f"unexpected lower-bridge block counts: {regime_counts}, {block_count}"
        )
    return {
        "segment_rows": len(records),
        "block_count": block_count,
        "regime_counts": regime_counts,
        "regime_maxima": regime_maxima,
        "global_maximum": global_maximum,
        "global_minimum_margin": global_minimum_margin,
    }


def build_artifact() -> dict:
    canonical_contract = segments.canonical_run_contract()
    stored_contract = json.loads(RUN_CONTRACT.read_text(encoding="utf-8"))
    if stored_contract != canonical_contract:
        raise RuntimeError("stored lower-bridge run contract changed")
    composition = load_and_compose()
    if Decimal(composition["global_maximum"]["scaled_curvature_upper"]) >= Decimal(
        THEOREM_CURVATURE_CONSTANT
    ):
        raise RuntimeError("lower bridge does not fit below the 4200 transfer cap")
    theorem = "z_1''(t)<=4200/t^2 for every real 1251<=t<=5700"
    exact = {
        "seventh_gap": "W(t)=8*B(t)-w(t-1)+2*w(t)-w(t+1)",
        "order9_coordinate": "z(t)=2*w(t)-s(t)+log(1-exp(-W(t)))",
        "curvature_formula": (
            "z''=2*w''-s''+phi(W)*W''-chi(W)*(W')^2"
        ),
        "upper_formula": (
            "z''<=2*w''-s''+phi(W)*max(W'',0)"
        ),
        "theorem": theorem,
    }
    source_contract = {
        **canonical_contract,
        "segment_cache": {
            "path": SEGMENT_CACHE.relative_to(REPO_ROOT).as_posix(),
            "sha256": segments.sha256(SEGMENT_CACHE),
            "row_count": composition["segment_rows"],
        },
        "stored_run_contract": {
            "path": RUN_CONTRACT.relative_to(REPO_ROOT).as_posix(),
            "sha256": segments.sha256(RUN_CONTRACT),
        },
    }
    rows = [
        BridgeRow(
            "co10llbc_01_sources",
            "interval_input",
            "ready_to_apply",
            "SHA-bound H2-H24 caches and exact half-grid point jets cover every localized collar.",
            "26600 hundredth H rows; 42180 tenth H rows; 8931 point rows",
            "First Newman summand only.",
            source_contract,
        ),
        BridgeRow(
            "co10llbc_02_component_formula",
            "exact_inequality",
            "ready_to_apply",
            "The componentwise stable-log formula preserves cancellation and drops only a nonpositive term.",
            exact["upper_formula"],
            "Exact differential inequality inside the positive W cone.",
        ),
        BridgeRow(
            "co10llbc_03_near_middle",
            "interval_theorem",
            "ready_to_apply",
            "Bidirectional quarter blocks cover the near and middle regimes without gaps.",
            "z_1''(t)<=4200/t^2 for 1251<=t<=1800",
            "Hundredth tiles through 1500, then tenth tiles through 1800.",
            {
                "near_blocks": composition["regime_counts"]["near"],
                "middle_blocks": composition["regime_counts"]["middle"],
                "near_maximum": composition["regime_maxima"]["near"],
                "middle_maximum": composition["regime_maxima"]["middle"],
            },
        ),
        BridgeRow(
            "co10llbc_04_far",
            "interval_theorem",
            "ready_to_apply",
            "Half-unit blocks cover the far lower bridge through the saddle handoff.",
            "z_1''(t)<=4200/t^2 for 1800<=t<=5700",
            "Tenth-tile first-summand interval theorem.",
            {
                "far_blocks": composition["regime_counts"]["far"],
                "far_maximum": composition["regime_maxima"]["far"],
            },
        ),
        BridgeRow(
            "co10llbc_05_composition",
            "exact_theorem_composition",
            "ready_to_apply",
            "All 9,996 passing blocks compose into one continuous lower-bridge theorem.",
            theorem,
            "First summand only; upper saddle ranges and full-kernel transfer remain separate.",
            {
                "global_maximum": composition["global_maximum"],
                "global_minimum_margin": composition[
                    "global_minimum_margin"
                ],
            },
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order10_localized_lower_bridge_certificate",
        "date": "2026-07-16",
        "status": "rigorous order-ten first-summand curvature theorem on 1251<=t<=5700",
        "proof_boundary": (
            "This artifact proves only the displayed first-summand continuous "
            "lower bridge. It does not prove the upper saddle ranges, full-kernel "
            "curvature ceiling, endpoint tail, delayed entry, PF-infinity, RH, "
            "or Lambda<=0."
        ),
        "theorem": theorem,
        "exact": exact,
        "source_contract": source_contract,
        "composition": composition,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": len(rows),
            "open_rows": 0,
            "segment_rows": composition["segment_rows"],
            "accepted_blocks": composition["block_count"],
            "near_blocks": composition["regime_counts"]["near"],
            "middle_blocks": composition["regime_counts"]["middle"],
            "far_blocks": composition["regime_counts"]["far"],
            "lower_bridge_curvature_theorems": 1,
            "full_half_line_theorems": 0,
            "full_kernel_theorems": 0,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order10_localized_lower_bridge_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order10_localized_lower_bridge_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    composition = artifact["composition"]
    lines = [
        "# Order-Ten Localized Lower-Bridge Certificate",
        "",
        "Date: 2026-07-16",
        "",
        "Status: rigorous first-summand continuous curvature theorem on",
        "`1251<=t<=5700`. This is not yet a full-kernel or RH theorem.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order10_localized_lower_bridge_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order10_localized_lower_bridge_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order10_localized_lower_bridge_certificate.py",
        "```",
        "",
        "## Theorem",
        "",
        "```text",
        artifact["theorem"],
        artifact["exact"]["curvature_formula"],
        artifact["exact"]["upper_formula"],
        "```",
        "",
        "The bridge contains `9,996` contiguous real-interval blocks:",
        "",
        "```text",
        f"near quarter blocks: {composition['regime_counts']['near']}",
        f"middle quarter blocks: {composition['regime_counts']['middle']}",
        f"far half blocks: {composition['regime_counts']['far']}",
        "largest scaled upper: "
        + composition["global_maximum"]["scaled_curvature_upper"],
        "smallest cap margin: "
        + composition["global_minimum_margin"]["curvature_margin_lower"],
        "```",
        "",
        "Every serialized Arb ball was parsed at 512-bit precision. The complete",
        "segment cache and all source caches are SHA-bound by the run contract.",
        "",
        "The next analytic task is the upper first-summand saddle range from",
        "`t=5700` onward, followed by the first-summand-to-full-kernel transfer.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-ten localized lower bridge: "
        f"{summary['segment_rows']} segments, "
        f"{summary['accepted_blocks']} blocks, "
        f"largest scaled upper "
        f"{artifact['composition']['global_maximum']['scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
