#!/usr/bin/env python3
"""Check the order-seven shifted-jet t=320..1000 certificate."""

from __future__ import annotations

from decimal import Decimal
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
ARTIFACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate.md"
)
START_T = Fraction(320)
END_T = Fraction(1000)
STEP_RATIO = Fraction(161, 160)
BLOCK_GRID = Fraction(1, 10)
COLLAR_RADIUS = Fraction(5)
EXPECTED_IDS = [
    "co7sjb_01_exact_potential_point_jets",
    "co7sjb_02_dimensionless_remainder",
    "co7sjb_03_stable_coordinate_domain",
    "co7sjb_04_continuous_curvature",
    "co7sjb_05_ray_handoff",
]
CACHE_CONTRACTS = {
    "fine_hundredth": {
        "cache": (
            "work/rh_compute/results/"
            "jensen_window_pf_compound_order7_shifted_jet_compact_h2_h22_hundredth_tiles.jsonl"
        ),
        "manifest": (
            "work/rh_compute/results/"
            "jensen_window_pf_compound_order7_shifted_jet_compact_h2_h22_hundredth_cache.json"
        ),
        "start": Fraction(315),
        "end": Fraction(505),
        "width": Fraction(1, 100),
        "rows": 19000,
    },
    "coarse_tenth": {
        "cache": (
            "work/rh_compute/results/"
            "jensen_window_pf_compound_order7_shifted_jet_compact_h2_h22_tenth_tiles.jsonl"
        ),
        "manifest": (
            "work/rh_compute/results/"
            "jensen_window_pf_compound_order7_shifted_jet_compact_h2_h22_cache.json"
        ),
        "start": Fraction(315),
        "end": Fraction(1005),
        "width": Fraction(1, 10),
        "rows": 6900,
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def geometric_blocks() -> list[tuple[Fraction, Fraction]]:
    blocks = []
    anchor = START_T
    while anchor < END_T:
        candidate = anchor * STEP_RATIO
        right = min(
            math.floor(candidate / BLOCK_GRID) * BLOCK_GRID,
            END_T,
        )
        if right <= anchor:
            raise RuntimeError("checker block grid made no progress")
        blocks.append((anchor, right))
        anchor = right
    return blocks


def issue(issues: list[str], message: str) -> None:
    issues.append(message)


def check_sources(artifact: dict, issues: list[str]) -> dict[str, str]:
    source_contract = artifact.get("source_contract", {})
    hashes = {}
    for label, expected in CACHE_CONTRACTS.items():
        recorded = source_contract.get(label, {})
        if recorded.get("cache") != expected["cache"]:
            issue(issues, f"bad {label} cache path")
        if recorded.get("manifest") != expected["manifest"]:
            issue(issues, f"bad {label} manifest path")
        if recorded.get("rows") != expected["rows"]:
            issue(issues, f"bad {label} source row count")
        if recorded.get("t_range") != [
            str(expected["start"]),
            str(expected["end"]),
        ]:
            issue(issues, f"bad {label} source range")
        if recorded.get("tile_width") != str(expected["width"]):
            issue(issues, f"bad {label} tile width")
        if recorded.get("h_derivative_orders") != [2, 22]:
            issue(issues, f"bad {label} derivative-order contract")

        cache_path = REPO_ROOT / expected["cache"]
        manifest_path = REPO_ROOT / expected["manifest"]
        if not cache_path.exists() or not manifest_path.exists():
            issue(issues, f"missing {label} source files")
            continue
        digest = sha256(cache_path)
        hashes[label] = digest
        if digest != recorded.get("sha256"):
            issue(issues, f"{label} cache hash mismatch")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if (
            manifest.get("cache", {}).get("sha256") != digest
            or manifest.get("cache", {}).get("row_count") != expected["rows"]
            or manifest.get("cache", {}).get("all_rows_passed") is not True
        ):
            issue(issues, f"{label} manifest contract failed")
        line_count = sum(
            bool(line.strip())
            for line in cache_path.read_text(encoding="utf-8").splitlines()
        )
        if line_count != expected["rows"]:
            issue(issues, f"{label} physical row count mismatch")
    return hashes


def check_blocks(
    artifact: dict,
    hashes: dict[str, str],
    issues: list[str],
) -> None:
    expected = geometric_blocks()
    blocks = artifact.get("blocks", [])
    if len(blocks) != len(expected):
        issue(issues, f"bad continuation block count: {len(blocks)}")
    for index, (block, endpoints) in enumerate(zip(blocks, expected)):
        anchor, right = endpoints
        if block.get("index") != index:
            issue(issues, f"bad block index {index}")
        if block.get("anchor") != str(anchor) or block.get("right") != str(right):
            issue(issues, f"bad block endpoints {index}")
        if block.get("width") != str(right - anchor):
            issue(issues, f"bad block width {index}")
        if block.get("relative_width") != str((right - anchor) / anchor):
            issue(issues, f"bad relative width {index}")
        if (right - anchor) / anchor > Fraction(1, 160):
            issue(issues, f"relative step too wide at block {index}")
        if block.get("passed") is not True:
            issue(issues, f"block not passed {index}")

        collar_left = anchor - COLLAR_RADIUS
        collar_right = right + COLLAR_RADIUS
        expected_label = (
            "fine_hundredth" if collar_right <= Fraction(505) else "coarse_tenth"
        )
        contract = CACHE_CONTRACTS[expected_label]
        first = max(
            0,
            math.floor((collar_left - contract["start"]) / contract["width"]),
        )
        last = min(
            contract["rows"],
            math.ceil((collar_right - contract["start"]) / contract["width"]),
        )
        if block.get("source_label") != expected_label:
            issue(issues, f"bad source selection at block {index}")
        if block.get("source_sha256") != hashes.get(expected_label):
            issue(issues, f"bad source hash at block {index}")
        if (
            block.get("source_first_index") != first
            or block.get("source_last_index") != last - 1
            or block.get("source_rows") != last - first
            or block.get("collar_left") != str(collar_left)
            or block.get("collar_right") != str(collar_right)
        ):
            issue(issues, f"bad source collar at block {index}")

        point_scaled = Decimal(block.get("point_scaled_curvature", "600"))
        remainder = Decimal(
            block.get("scaled_curvature_remainder_upper", "600")
        )
        scaled = Decimal(block.get("scaled_curvature_upper", "600"))
        margin = Decimal(block.get("curvature_margin_lower", "0"))
        rho = Decimal(
            block.get("normalized_r_remainder_coefficient_upper", "0")
        )
        if not point_scaled < 600:
            issue(issues, f"point curvature failed at block {index}")
        if not remainder >= 0 or not rho > 0:
            issue(issues, f"bad normalized remainder at block {index}")
        if not scaled < 600 or not margin > 0:
            issue(issues, f"continued curvature failed at block {index}")
        if scaled + margin > Decimal(600):
            issue(issues, f"inconsistent curvature margin at block {index}")

        point_coordinates = block.get("point_coordinates", {})
        common_coordinates = block.get("common_coordinates", {})
        for name in ("B", "J", "R", "S", "T"):
            if Decimal(point_coordinates.get(name, "0")) <= 0:
                issue(issues, f"nonpositive point {name} at block {index}")
            coordinate = common_coordinates.get(name, {})
            lower = Decimal(coordinate.get("dimensionless_common_lower", "0"))
            upper = Decimal(coordinate.get("dimensionless_common_upper", "0"))
            if not 0 < lower <= upper:
                issue(issues, f"bad common {name} at block {index}")

        quadrature = block.get("quadrature", {})
        brackets = quadrature.get("mode_brackets", [])
        if quadrature.get("shift_count") != 11:
            issue(issues, f"bad quadrature shift count at block {index}")
        if Decimal(quadrature.get("maximum_panel_error_upper", "1")) >= Decimal(
            "1e-12"
        ):
            issue(issues, f"panel error too large at block {index}")
        if Decimal(quadrature.get("maximum_tail_moment_upper", "1")) >= Decimal(
            "1e-12"
        ):
            issue(issues, f"tail error too large at block {index}")
        if len(brackets) != 11 or any(
            len(pair) != 2 or Fraction(pair[0]) >= Fraction(pair[1])
            for pair in brackets
        ):
            issue(issues, f"bad mode brackets at block {index}")

    for previous, current in zip(blocks, blocks[1:]):
        if previous.get("right") != current.get("anchor"):
            issue(issues, "continuation block gap")
    if blocks and (
        blocks[0].get("anchor") != "320" or blocks[-1].get("right") != "1000"
    ):
        issue(issues, "continuation endpoints changed")

    if blocks:
        summary = artifact.get("summary", {})
        largest = max(
            (block["scaled_curvature_upper"] for block in blocks),
            key=Decimal,
        )
        weakest_margin = min(
            (block["curvature_margin_lower"] for block in blocks),
            key=Decimal,
        )
        weakest_t = min(
            (
                block["common_coordinates"]["T"][
                    "dimensionless_common_lower"
                ]
                for block in blocks
            ),
            key=Decimal,
        )
        source_counts = {
            label: sum(block.get("source_label") == label for block in blocks)
            for label in CACHE_CONTRACTS
        }
        if summary.get("blocks") != len(expected):
            issue(issues, "bad summary block count")
        if summary.get("source_blocks") != source_counts:
            issue(issues, "bad source-block summary")
        if summary.get("largest_scaled_curvature_upper") != largest:
            issue(issues, "bad largest-curvature summary")
        if summary.get("weakest_curvature_margin_lower") != weakest_margin:
            issue(issues, "bad weakest-margin summary")
        if summary.get("weakest_dimensionless_T_lower") != weakest_t:
            issue(issues, "bad weakest-T summary")
        if summary.get("all_blocks_passed") is not True:
            issue(issues, "summary does not mark all blocks passed")


def main() -> int:
    issues: list[str] = []
    if not ARTIFACT.exists():
        print(f"missing artifact: {ARTIFACT}")
        return 1
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate"
    ):
        issue(issues, "bad artifact kind")
    if artifact.get("date") != "2026-07-13":
        issue(issues, "bad artifact date")
    if artifact.get("status") != (
        "rigorous continuous fourth-nested first-summand curvature theorem on 320<=t<=1000"
    ):
        issue(issues, "bad artifact status")
    if "does not prove" not in artifact.get("proof_boundary", ""):
        issue(issues, "missing proof boundary")

    parameters = artifact.get("parameters", {})
    expected_parameters = {
        "precision_bits": 384,
        "point_order": 10,
        "remainder_order": 11,
        "curvature_constant": 600,
        "start_t": "320",
        "end_t": "1000",
        "step_ratio": "161/160",
        "block_grid": "1/10",
        "collar_radius": "5",
    }
    for key, value in expected_parameters.items():
        if parameters.get(key) != value:
            issue(issues, f"bad parameter {key}")

    hashes = check_sources(artifact, issues)
    check_blocks(artifact, hashes, issues)

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issue(issues, "bad certificate row ids")
    if [row.get("readiness") for row in rows] != [
        "ready_to_apply",
        "ready_to_apply",
        "ready_to_apply",
        "ready_to_apply",
        "not_ready_to_apply",
    ]:
        issue(issues, "bad certificate readiness sequence")
    if len(rows) == 5:
        if rows[3].get("formula") != (
            "r_1''(t)<=600/t^2 for 320<=t<=1000"
        ):
            issue(issues, "bad continuous theorem formula")
        if rows[4].get("formula") != "r_1''(t)<=600/t^2 for t>=1000":
            issue(issues, "bad ray handoff formula")

    if artifact.get("generator") != (
        "work/rh_compute/scripts/"
        "jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate.py"
    ):
        issue(issues, "bad generator path")
    if artifact.get("checker") != (
        "work/rh_compute/scripts/"
        "check_jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate.py"
    ):
        issue(issues, "bad checker path")
    if not NOTE.exists():
        issue(issues, "missing markdown note")
    else:
        note = NOTE.read_text(encoding="utf-8")
        for phrase in (
            "rigorous continuous fourth-nested first-summand curvature",
            "r_1''(t)<=600/t^2 for every real 320<=t<=1000",
            "not a proof of the remaining ray",
        ):
            if phrase not in note:
                issue(issues, f"missing note phrase: {phrase}")

    if issues:
        for item in issues:
            print(f"ISSUE: {item}")
        print(
            "failed order-seven shifted-jet t=320..1000 certificate: "
            f"{len(issues)} issues"
        )
        return 1
    print(
        "validated order-seven shifted-jet t=320..1000 certificate: "
        f"5 rows, {len(geometric_blocks())} contiguous blocks, "
        "11 shifts per anchor, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
