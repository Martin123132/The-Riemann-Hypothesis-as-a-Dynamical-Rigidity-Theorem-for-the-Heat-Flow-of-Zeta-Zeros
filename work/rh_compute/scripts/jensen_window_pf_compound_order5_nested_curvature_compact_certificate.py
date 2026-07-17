#!/usr/bin/env python3
"""Certify the nested order-five curvature bound from t=320 to mode u=2."""

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
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as order4_compact  # noqa: E402
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


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order5_nested_curvature_compact_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order5_nested_curvature_compact_certificate.md"
)
SOURCE_COMPACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_localized_curvature_compact_certificate.json"
)
SOURCE_BRIDGE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order5_first_summand_curvature_bridge.json"
)
TARGET_T = 320
COLLAR_T = 3
INITIAL_CENTRAL_TILE_COUNT = 3000


@dataclass(frozen=True)
class CompactRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_contract() -> dict:
    compact = load_json(SOURCE_COMPACT)
    bridge = load_json(SOURCE_BRIDGE)
    cache = compact.get("cache", {})
    if compact.get("summary", {}).get("all_blocks_passed") is not True:
        raise RuntimeError("order-four compact H cache is not certified")
    if cache.get("path") != order4_compact.DEFAULT_CACHE.relative_to(
        REPO_ROOT
    ).as_posix():
        raise RuntimeError("order-four compact cache path changed")
    if bridge.get("exact", {}).get("continuous_target") != (
        "q_1''(t)<=60/t^2 for every real t>=320"
    ):
        raise RuntimeError("order-five continuous target changed")
    actual_hash = order4_compact.cache_sha256(order4_compact.DEFAULT_CACHE)
    if actual_hash != cache.get("sha256"):
        raise RuntimeError("order-four compact H cache hash changed")
    return {
        "cache_path": cache["path"],
        "cache_sha256": actual_hash,
        "cache_rows": cache["row_count"],
        "source_t_range": compact["mode_handoff"]["certified_t_range"],
        "continuous_target": bridge["exact"]["continuous_target"],
    }


def deterministic_tasks() -> list[tuple[int, Fraction, Fraction]]:
    return [
        (index, left, right)
        for index, (left, right) in enumerate(
            order4_compact.fraction_grid(
                order4_compact.OUTER_START,
                order4_compact.OUTER_END,
                order4_compact.TILE_WIDTH,
            )
        )
    ]


def first_target_tile(records: list[dict]) -> int:
    target = flint.arb(TARGET_T)
    for index, record in enumerate(records):
        left = order4_compact.interval_from_text(record["t_left"])
        right = order4_compact.interval_from_text(record["t_right"])
        if bool(right > target):
            if not bool(left < target):
                raise RuntimeError("compact target tile does not straddle t=320")
            return index
    raise RuntimeError("compact H cache does not reach t=320")


def collar_indices(
    records: list[dict], central_left: int, central_right: int
) -> tuple[int, int, dict]:
    lower_t = order4_compact.interval_from_text(records[central_left]["t_left"])
    upper_t = order4_compact.interval_from_text(
        records[central_right - 1]["t_right"]
    )
    left_index = central_left
    while left_index > 0:
        candidate = order4_compact.interval_from_text(
            records[left_index]["t_left"]
        )
        if bool(candidate < lower_t - COLLAR_T):
            break
        left_index -= 1
    right_index = central_right - 1
    while right_index + 1 < len(records):
        candidate = order4_compact.interval_from_text(
            records[right_index]["t_right"]
        )
        if bool(candidate > upper_t + COLLAR_T):
            break
        right_index += 1
    outer_lower = order4_compact.interval_from_text(
        records[left_index]["t_left"]
    )
    outer_upper = order4_compact.interval_from_text(
        records[right_index]["t_right"]
    )
    if not bool(
        outer_lower < lower_t - COLLAR_T
        and outer_upper > upper_t + COLLAR_T
    ):
        raise RuntimeError("cached H tiles do not contain the required t+-3 collar")
    return left_index, right_index, {
        "central_t_lower_ball": lower_t.str(40).replace("e", "E"),
        "central_t_upper_ball": upper_t.str(40).replace("e", "E"),
        "left_t_collar_ball": (lower_t - outer_lower).str(40).replace("e", "E"),
        "right_t_collar_ball": (outer_upper - upper_t).str(40).replace("e", "E"),
    }


def compact_block(
    records: list[dict], central_left: int, central_right: int
) -> dict:
    cover_left, cover_right, collar = collar_indices(
        records, central_left, central_right
    )
    derivatives, diagnostics = order4_compact.derivative_cover(
        records, cover_left, cover_right
    )
    diagnostics.update(collar)
    left = Fraction(records[central_left]["mode_left"])
    right = Fraction(records[central_right - 1]["mode_right"])
    result = evaluate_nested_curvature_from_h_cover(
        left,
        right,
        derivatives,
        cover_diagnostics=diagnostics,
    )
    result["central_tile_index_first"] = central_left
    result["central_tile_index_last"] = central_right - 1
    result["central_tile_count"] = central_right - central_left
    result["cover_tile_count"] = cover_right - cover_left + 1
    return result


def certify_adaptive_block(
    records: list[dict], left: int, right: int
) -> list[dict]:
    result = compact_block(records, left, right)
    if result.get("passed"):
        return [result]
    if right - left <= 1:
        raise RuntimeError(
            f"nested compact curvature failed on irreducible tile {left}: "
            f"{result.get('failure')} {result.get('detail', '')}"
        )
    midpoint = (left + right) // 2
    return certify_adaptive_block(records, left, midpoint) + certify_adaptive_block(
        records, midpoint, right
    )


def finite_certificate(records: list[dict]) -> dict:
    flint.ctx.prec = order4_compact.DEFAULT_PRECISION_BITS
    tasks = deterministic_tasks()
    if len(records) != len(tasks):
        raise RuntimeError(
            f"compact assembly needs {len(tasks)} cache rows, found {len(records)}"
        )
    central_left = first_target_tile(records)
    central_right = order4_compact.mode_index(Fraction(2))
    accepted: list[dict] = []
    cursor = central_left
    while cursor < central_right:
        endpoint = min(
            cursor + INITIAL_CENTRAL_TILE_COUNT,
            central_right,
        )
        accepted.extend(certify_adaptive_block(records, cursor, endpoint))
        cursor = endpoint

    for previous, current in zip(accepted, accepted[1:]):
        if previous["central_mode"][1] != current["central_mode"][0]:
            raise RuntimeError("nested compact cover has a mode gap")
    if accepted[-1]["central_mode"][1] != "2":
        raise RuntimeError("nested compact cover does not reach mode u=2")

    start_mode = Fraction(records[central_left]["mode_left"])
    start_t = potential_jet_arb(arb_rational(start_mode), 1)[1]
    end_t = potential_jet_arb(flint.arb(2), 1)[1]
    first_right_t = order4_compact.interval_from_text(
        records[central_left]["t_right"]
    )
    if not bool(start_t < TARGET_T < first_right_t):
        raise RuntimeError("first nested compact block does not cover t=320")

    worst_margin = min(
        range(len(accepted)),
        key=lambda index: Decimal(accepted[index]["margin_lower"]),
    )
    largest_scaled = max(
        range(len(accepted)),
        key=lambda index: Decimal(accepted[index]["scaled_curvature_upper"]),
    )
    weakest_J = min(
        range(len(accepted)),
        key=lambda index: Decimal(accepted[index]["J_lower"]),
    )
    weakest_R = min(
        range(len(accepted)),
        key=lambda index: Decimal(accepted[index]["R_lower"]),
    )
    selected_indices = sorted(
        {
            0,
            len(accepted) // 4,
            len(accepted) // 2,
            3 * len(accepted) // 4,
            len(accepted) - 1,
            worst_margin,
            largest_scaled,
            weakest_J,
            weakest_R,
        }
    )

    def summary_row(index: int) -> dict:
        row = accepted[index]
        return {
            "index": index,
            "central_mode": row["central_mode"],
            "central_t_ball": row["central_t_ball"],
            "central_tile_count": row["central_tile_count"],
            "cover_tile_count": row["cover_tile_count"],
            "J_lower": row["J_lower"],
            "R_lower": row["R_lower"],
            "q_second_upper": row["q_second_upper"],
            "target_lower": row["target_lower"],
            "margin_lower": row["margin_lower"],
            "scaled_curvature_upper": row["scaled_curvature_upper"],
        }

    return {
        "central_start_mode": str(start_mode),
        "central_end_mode": "2",
        "central_start_t_ball": start_t.str(50).replace("e", "E"),
        "central_end_t_ball": end_t.str(50).replace("e", "E"),
        "target_t": TARGET_T,
        "certified_t_range": "320<=t<=V'(2)",
        "accepted": accepted,
        "selected": [summary_row(index) for index in selected_indices],
        "summary": {
            "accepted_blocks": len(accepted),
            "minimum_central_tile_count": min(
                row["central_tile_count"] for row in accepted
            ),
            "maximum_central_tile_count": max(
                row["central_tile_count"] for row in accepted
            ),
            "worst_margin_index": worst_margin,
            "worst_margin_lower": accepted[worst_margin]["margin_lower"],
            "largest_scaled_index": largest_scaled,
            "largest_scaled_curvature_upper": accepted[largest_scaled][
                "scaled_curvature_upper"
            ],
            "weakest_J_index": weakest_J,
            "weakest_J_lower": accepted[weakest_J]["J_lower"],
            "weakest_R_index": weakest_R,
            "weakest_R_lower": accepted[weakest_R]["R_lower"],
            "all_blocks_passed": True,
        },
    }


def build_artifact() -> dict:
    contract = source_contract()
    records = order4_compact.load_cache(
        order4_compact.DEFAULT_CACHE, deterministic_tasks()
    )
    finite = finite_certificate(records)
    summary = finite["summary"]
    rows = [
        CompactRow(
            id="co5nccc_01_tent_hull_core",
            role="exact_interval_lemma",
            readiness="ready_to_apply",
            claim="A common H-derivative hull encloses every nested tent average needed by q_1''.",
            formula=(
                "B^(r)=tent(H^(r+2)); J^(r)=2B^(r)-tent(ell^(r+2)); "
                "R^(r)=3B^(r)-tent(h^(r+2))"
            ),
            proof_boundary="Inclusion-monotone interval algebra through derivative order eight.",
        ),
        CompactRow(
            id="co5nccc_02_three_unit_collar",
            role="exact_geometry_lemma",
            readiness="ready_to_apply",
            claim="The cached outer tiles cover every shift introduced by the two nested centered differences.",
            formula="central t plus or minus 3 lies in the common H^(2)..H^(8) cover",
            proof_boundary="Checked separately on every accepted central block.",
        ),
        CompactRow(
            id="co5nccc_03_positive_stable_coordinates",
            role="interval_theorem",
            readiness="ready_to_apply",
            claim="Both stable logarithmic coordinates remain strictly positive on every compact block.",
            formula="J_1(t)>0 and R_1(t)>0, 320<=t<=V'(2)",
            proof_boundary="Outward-rounded Arb hulls from the hashed quadrature cache.",
            diagnostics={
                "weakest_J_lower": summary["weakest_J_lower"],
                "weakest_R_lower": summary["weakest_R_lower"],
            },
        ),
        CompactRow(
            id="co5nccc_04_compact_curvature",
            role="interval_theorem",
            readiness="ready_to_apply",
            claim="The continuous first-summand nested curvature ceiling holds on the complete compact range.",
            formula="q_1''(t)<=60/t^2 for every 320<=t<=V'(2)",
            proof_boundary="Real-parameter interval cover, not point sampling.",
            diagnostics={
                "accepted_blocks": summary["accepted_blocks"],
                "largest_scaled_curvature_upper": summary[
                    "largest_scaled_curvature_upper"
                ],
                "worst_margin_lower": summary["worst_margin_lower"],
            },
        ),
        CompactRow(
            id="co5nccc_05_ray_handoff",
            role="open_theorem_target",
            readiness="not_ready_to_apply",
            claim="Complete the same curvature ceiling on the remaining analytic mode ray.",
            formula="prove q_1''(t)<=60/t^2 for every mode u>=2",
            proof_boundary="Open ray only; the compact t range is closed here.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order5_nested_curvature_compact_certificate",
        "date": "2026-07-13",
        "status": (
            "rigorous compact first-summand order-five curvature theorem with "
            "one open analytic ray"
        ),
        "proof_boundary": (
            "This artifact proves q_1''(t)<=60/t^2 only for "
            "320<=t<=V'(2). It does not prove the u>=2 ray, full order-five "
            "entry, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_compound_order5_first_summand_curvature_bridge.md",
            "outputs/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.md",
            "outputs/formal_core.md",
        ],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order5_nested_curvature_compact_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order5_nested_curvature_compact_certificate.py"
        ),
        "parameters": {
            "precision_bits": order4_compact.DEFAULT_PRECISION_BITS,
            "target_t": TARGET_T,
            "curvature_constant": CURVATURE_CONSTANT,
            "collar_t": COLLAR_T,
            "initial_central_tile_count": INITIAL_CENTRAL_TILE_COUNT,
        },
        "source_contract": contract,
        "compact": {
            key: value
            for key, value in finite.items()
            if key not in {"accepted"}
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": sum(row.readiness == "ready_to_apply" for row in rows),
            "open_rows": sum(
                row.readiness == "not_ready_to_apply" for row in rows
            ),
            "accepted_blocks": summary["accepted_blocks"],
            "positive_stable_layers": 2,
            "compact_curvature_theorems": 1,
            "open_rays": 1,
        },
    }


def write_note(path: Path, artifact: dict) -> None:
    compact = artifact["compact"]
    summary = compact["summary"]
    lines = [
        "# Jensen-Window PF Compound Order-Five Nested Curvature Compact Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous compact first-summand order-five curvature theorem",
        "with one open analytic ray. This is not a proof of full order-five",
        "entry, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order5_nested_curvature_compact_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order5_nested_curvature_compact_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order5_nested_curvature_compact_certificate.py",
        "```",
        "",
        "## Interval Core",
        "",
        "The two stable centered differences require a common three-unit",
        "collar. On each central block, a single outward-rounded hull for",
        "`H^(2),...,H^(8)` encloses",
        "",
        "```text",
        "B^(r)=integral_[-1,1](1-|s|)*H^(r+2)(t+s) ds",
        "J^(r)=2*B^(r)-integral_[-1,1](1-|s|)*ell^(r+2)(t+s) ds",
        "R^(r)=3*B^(r)-integral_[-1,1](1-|s|)*h^(r+2)(t+s) ds.",
        "```",
        "",
        "Truncated Arb Taylor-series arithmetic then evaluates",
        "`ell=log(1-exp(-B))`, `h=2ell+log(1-exp(-J))`, and",
        "`q=2h-ell+log(1-exp(-R))` through second order. These are finite",
        "jet identities, not asymptotic series.",
        "",
        "## Cached Cover",
        "",
        f"The source cache has `{artifact['source_contract']['cache_rows']}`",
        "paired interval-Simpson tiles and SHA-256",
        "",
        "```text",
        artifact["source_contract"]["cache_sha256"],
        "```",
        "",
        f"Adaptive assembly accepts `{summary['accepted_blocks']}` central blocks.",
        f"The first block straddles `t={compact['target_t']}` and the last reaches",
        "the exact mode `u=2`, so the certified range is",
        "",
        "```text",
        compact["certified_t_range"],
        "```",
        "",
        "Every block proves `J_1>0`, `R_1>0`, and a positive curvature",
        "margin. The global diagnostics are",
        "",
        "```text",
        f"weakest J lower: {summary['weakest_J_lower']}",
        f"weakest R lower: {summary['weakest_R_lower']}",
        f"largest t^2*q_1'' upper: {summary['largest_scaled_curvature_upper']}",
        f"worst absolute margin: {summary['worst_margin_lower']}",
        "```",
        "",
        "Thus the compact theorem is",
        "",
        "```text",
        "q_1''(t)<=60/t^2 for every 320<=t<=V'(2).",
        "```",
        "",
        "## Remaining Ray",
        "",
        "The only unproved part of the continuous first-summand target is",
        "",
        "```text",
        "q_1''(t)<=60/t^2 for every mode u>=2.",
        "```",
        "",
        "The exact-cumulant corridor cover on `2<=u<=20` and normalized-H",
        "boxes on `u>=20` are the natural inputs for that final ray theorem.",
        "",
        "| mode interval | t ball | t^2 q_1'' upper | margin lower | J lower | R lower |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in compact["selected"]:
        lines.append(
            f"| `{row['central_mode'][0]}..{row['central_mode'][1]}` | "
            f"`{row['central_t_ball']}` | `{row['scaled_curvature_upper']}` | "
            f"`{row['margin_lower']}` | `{row['J_lower']}` | `{row['R_lower']}` |"
        )
    lines.extend(
        [
            "",
            "```text",
            "outputs/jensen_window_pf_compound_order5_first_summand_curvature_bridge.md",
            "outputs/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.md",
            "outputs/signed_hankel_jensen_bridge_target.md",
            "outputs/formal_core.md",
            "```",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as handle:
        json.dump(artifact, handle, indent=2, sort_keys=True)
        handle.write("\n")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-five nested curvature compact certificate: "
        f"{summary['accepted_blocks']} blocks, "
        f"{summary['positive_stable_layers']} positive stable layers, "
        f"{summary['compact_curvature_theorems']} compact curvature theorem, "
        f"{summary['open_rays']} open ray"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
