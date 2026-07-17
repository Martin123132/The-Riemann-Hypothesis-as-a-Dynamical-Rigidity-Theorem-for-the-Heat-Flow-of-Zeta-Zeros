#!/usr/bin/env python3
"""Scout the seventh-nested first-summand coordinate at rigorous exact points."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

from jensen_window_pf_compound_order4_localized_curvature_compact_certificate import (  # noqa: E402
    interval_from_text,
)
from jensen_window_pf_compound_order5_nested_curvature_interval_core import (  # noqa: E402
    series_add,
    series_scale,
    series_sub,
    stable_log_series,
)
from jensen_window_pf_compound_order7_shifted_jet_continuation_core import (  # noqa: E402
    centered_second_difference,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_upper_text,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_first_summand_point_scout.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order10_first_summand_point_scout.md"
)
POINT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_shifted_point_h0_h8_half_grid.jsonl"
)
POINT_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_shifted_point_h0_h8_cache.json"
)
REDUCTION_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_m100_tail_curvature_reduction.json"
)
SELECTED_ANCHORS = (
    Fraction(1251),
    Fraction(1300),
    Fraction(1500),
    Fraction(2000),
    Fraction(3000),
    Fraction(4000),
)
POINT_JET_ORDER = 4
POINT_CURVATURE_CAP = 250
FULL_TAIL_CURVATURE_TARGET = 5500
PRECISION_BITS = 512


@dataclass(frozen=True)
class ScoutRow:
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


def ball_text(value: flint.arb, digits: int = 70) -> str:
    return value.str(digits).replace("e", "E")


def load_source_contract() -> dict:
    manifest = json.loads(POINT_MANIFEST.read_text(encoding="utf-8"))
    reduction = json.loads(REDUCTION_SOURCE.read_text(encoding="utf-8"))
    cache = manifest.get("cache", {})
    parameters = manifest.get("parameters", {})
    if manifest.get("kind") != "jensen_window_pf_compound_order9_shifted_point_h0_h8_cache":
        raise RuntimeError("point-cache manifest kind changed")
    if cache.get("row_count") != 8929 or cache.get("all_rows_passed") is not True:
        raise RuntimeError("point-cache row contract changed")
    if cache.get("h_derivative_orders") != [0, 8]:
        raise RuntimeError("point-cache derivative range changed")
    if cache.get("sha256") != sha256(POINT_CACHE):
        raise RuntimeError("point-cache hash changed")
    if parameters.get("start_t") != "1243" or parameters.get("end_t") != "5707":
        raise RuntimeError("point-cache t range changed")
    if parameters.get("step_t") != "1/2" or parameters.get("precision_bits") != PRECISION_BITS:
        raise RuntimeError("point-cache grid contract changed")
    if reduction.get("summary", {}).get("open_curvature_targets") != 1:
        raise RuntimeError("order-ten curvature reduction source changed")
    if reduction.get("exact", {}).get("sufficient_ceiling") != (
        "Z_k<=5500/k^2 for every integer k>=1252"
    ):
        raise RuntimeError("order-ten curvature target changed")
    return {
        "point_manifest": POINT_MANIFEST.relative_to(REPO_ROOT).as_posix(),
        "point_manifest_sha256": sha256(POINT_MANIFEST),
        "point_cache": POINT_CACHE.relative_to(REPO_ROOT).as_posix(),
        "point_cache_sha256": cache["sha256"],
        "point_cache_rows": cache["row_count"],
        "point_cache_range": [parameters["start_t"], parameters["end_t"]],
        "point_cache_step": parameters["step_t"],
        "point_cache_precision_bits": parameters["precision_bits"],
        "reduction": REDUCTION_SOURCE.relative_to(REPO_ROOT).as_posix(),
        "reduction_sha256": sha256(REDUCTION_SOURCE),
    }


def needed_targets() -> set[Fraction]:
    return {
        anchor + shift
        for anchor in SELECTED_ANCHORS
        for shift in range(-8, 9)
    }


def load_h_jets() -> dict[Fraction, list[flint.arb]]:
    needed = needed_targets()
    rows: dict[Fraction, list[flint.arb]] = {}
    with POINT_CACHE.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            target = Fraction(record["target_t"])
            if target not in needed:
                continue
            derivatives = record.get("h_derivatives", {})
            if not all(str(degree) in derivatives for degree in range(POINT_JET_ORDER + 1)):
                raise RuntimeError(f"point cache misses derivatives at t={target}")
            rows[target] = [
                interval_from_text(derivatives[str(degree)]) / math.factorial(degree)
                for degree in range(POINT_JET_ORDER + 1)
            ]
    missing = sorted(needed - set(rows))
    if missing:
        raise RuntimeError(f"point cache misses selected targets: {missing[:5]}")
    return rows


def exact_point_hierarchy(
    anchor: Fraction,
    h_jets: dict[Fraction, list[flint.arb]],
) -> dict:
    order = POINT_JET_ORDER
    source = {shift: h_jets[anchor + shift] for shift in range(-8, 9)}
    b = {
        shift: centered_second_difference(source, shift)
        for shift in range(-7, 8)
    }
    ell = {shift: stable_log_series(value, order) for shift, value in b.items()}
    j_gap = {
        shift: series_sub(
            series_scale(b[shift], 2),
            centered_second_difference(ell, shift),
        )
        for shift in range(-6, 7)
    }
    h = {
        shift: series_add(
            series_scale(ell[shift], 2),
            stable_log_series(j_gap[shift], order),
        )
        for shift in range(-6, 7)
    }
    r_gap = {
        shift: series_sub(
            series_scale(b[shift], 3),
            centered_second_difference(h, shift),
        )
        for shift in range(-5, 6)
    }
    q = {
        shift: series_add(
            series_sub(series_scale(h[shift], 2), ell[shift]),
            stable_log_series(r_gap[shift], order),
        )
        for shift in range(-5, 6)
    }
    s_gap = {
        shift: series_sub(
            series_scale(b[shift], 4),
            centered_second_difference(q, shift),
        )
        for shift in range(-4, 5)
    }
    p = {
        shift: series_add(
            series_sub(series_scale(q[shift], 2), h[shift]),
            stable_log_series(s_gap[shift], order),
        )
        for shift in range(-4, 5)
    }
    t_gap = {
        shift: series_sub(
            series_scale(b[shift], 5),
            centered_second_difference(p, shift),
        )
        for shift in range(-3, 4)
    }
    r = {
        shift: series_add(
            series_sub(series_scale(p[shift], 2), q[shift]),
            stable_log_series(t_gap[shift], order),
        )
        for shift in range(-3, 4)
    }
    u_gap = {
        shift: series_sub(
            series_scale(b[shift], 6),
            centered_second_difference(r, shift),
        )
        for shift in range(-2, 3)
    }
    s = {
        shift: series_add(
            series_sub(series_scale(r[shift], 2), p[shift]),
            stable_log_series(u_gap[shift], order),
        )
        for shift in range(-2, 3)
    }
    v_gap = {
        shift: series_sub(
            series_scale(b[shift], 7),
            centered_second_difference(s, shift),
        )
        for shift in range(-1, 2)
    }
    w = {
        shift: series_add(
            series_sub(series_scale(s[shift], 2), r[shift]),
            stable_log_series(v_gap[shift], order),
        )
        for shift in range(-1, 2)
    }
    w_gap = series_sub(
        series_scale(b[0], 8),
        centered_second_difference(w, 0),
    )
    z = series_add(
        series_sub(series_scale(w[0], 2), s[0]),
        stable_log_series(w_gap, order),
    )
    coordinates = {
        "B": b[0][0],
        "J": j_gap[0][0],
        "R": r_gap[0][0],
        "S": s_gap[0][0],
        "T": t_gap[0][0],
        "U": u_gap[0][0],
        "V": v_gap[0][0],
        "W": w_gap[0],
    }
    if not all(bool(value > 0) for value in coordinates.values()):
        raise RuntimeError(f"selected point lost a positive coordinate at t={anchor}")
    anchor_arb = flint.arb(anchor.numerator) / anchor.denominator
    w_scaled_second = 2 * w[0][2] * anchor_arb**2
    z_scaled_second = 2 * z[2] * anchor_arb**2
    z_scaled_first = z[1] * anchor_arb
    cap_margin = flint.arb(POINT_CURVATURE_CAP) - z_scaled_second
    if not bool(cap_margin > 0):
        raise RuntimeError(f"selected point exceeds scout cap at t={anchor}")
    return {
        "t": str(anchor),
        "coordinate_balls": {
            name: ball_text(value) for name, value in coordinates.items()
        },
        "coordinate_lower_bounds": {
            name: arb_lower_text(value) for name, value in coordinates.items()
        },
        "w1_scaled_second_ball": ball_text(w_scaled_second),
        "z1_scaled_first_ball": ball_text(z_scaled_first),
        "z1_scaled_second_ball": ball_text(z_scaled_second),
        "z1_scaled_second_upper": arb_upper_text(z_scaled_second),
        "point_cap": POINT_CURVATURE_CAP,
        "point_cap_margin_lower": arb_lower_text(cap_margin),
        "passed": True,
        "_W": coordinates["W"],
        "_z_scaled_second": z_scaled_second,
    }


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    source_contract = load_source_contract()
    h_jets = load_h_jets()
    calculated = [exact_point_hierarchy(anchor, h_jets) for anchor in SELECTED_ANCHORS]
    maximum = max(
        calculated,
        key=lambda row: row["_z_scaled_second"].upper(),
    )
    minimum_w = min(
        calculated,
        key=lambda row: row["_W"].lower(),
    )
    selected = []
    for row in calculated:
        selected.append({key: value for key, value in row.items() if not key.startswith("_")})
    diagnostics = {
        "selected_points": selected,
        "maximum_scaled_curvature_anchor": maximum["t"],
        "maximum_scaled_curvature_upper": maximum["z1_scaled_second_upper"],
        "minimum_W_anchor": minimum_w["t"],
        "minimum_W_lower": arb_lower_text(minimum_w["_W"]),
        "point_curvature_cap": POINT_CURVATURE_CAP,
        "full_tail_curvature_target": FULL_TAIL_CURVATURE_TARGET,
    }
    rows = [
        ScoutRow(
            "co10fsp_01_exact_point_cache",
            "interval_input",
            "ready_to_apply",
            "The reusable exact-point cache supplies rigorous first-summand H jets through derivative order eight.",
            "H_1^(j)(t) enclosed for j=0,...,8 on the exact half-grid 1243<=t<=5707",
            "Computational first-summand input only.",
        ),
        ScoutRow(
            "co10fsp_02_seventh_nested_hierarchy",
            "exact_interval_algebra",
            "ready_to_apply",
            "The signed-condensation hierarchy extends through the new W gap and z logarithm at every selected point.",
            "W=8*B-Delta^2*w; z=2*w-s+log(1-exp(-W))",
            "Exact-point jets only; no interval between selected points is covered.",
        ),
        ScoutRow(
            "co10fsp_03_positive_coordinates",
            "rigorous_point_scout",
            "ready_to_apply",
            "Every selected point has strict positive B,J,R,S,T,U,V,W enclosures.",
            "B,J,R,S,T,U,V,W>0 at t in {1251,1300,1500,2000,3000,4000}",
            "Six exact points, not a continuous ray theorem.",
            diagnostics,
        ),
        ScoutRow(
            "co10fsp_04_curvature_headroom",
            "rigorous_point_scout",
            "ready_to_apply",
            "The seventh-nested first-summand scaled curvature stays below 250 at every selected point.",
            "t^2*z_1''(t)<250<5500 at every selected t",
            "The 5500 ceiling is a full-kernel discrete-tail target; this row is first-summand and pointwise only.",
            diagnostics,
        ),
        ScoutRow(
            "co10fsp_05_continuum_handoff",
            "analytic_handoff",
            "not_ready_to_apply",
            "Upgrade the selected-point headroom to a continuous first-summand theorem and then bound the higher-summand transfer.",
            "prove z_1''(t)<=C/t^2 on t>=1251 and |Z_k-Z_k^(1)|<D/k^2 for k>=1252 with C+D<5500",
            "Requires localized tent-convolution enclosures one stable layer beyond order nine.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order10_first_summand_point_scout",
        "date": "2026-07-16",
        "status": (
            "rigorous selected-point seventh-nested first-summand scout; "
            "continuous and full-kernel bounds remain open"
        ),
        "proof_boundary": (
            "This artifact rigorously evaluates the first-summand B through W "
            "coordinates and z_1 curvature at six exact points. It does not "
            "cover the intervals between those points, prove the full-kernel "
            "5500/k^2 ceiling, close the order-ten endpoint tail, prove delayed "
            "entry, PF-infinity, RH, or Lambda<=0."
        ),
        "source_contract": source_contract,
        "exact": {
            "seventh_gap": "W=8*B-Delta^2*w",
            "order9_log_coordinate": "z=2*w-s+log(1-exp(-W))",
            "scaled_curvature": "t^2*z_1''(t)",
            "selected_point_cap": "t^2*z_1''(t)<250",
            "full_tail_target": "Z_k<=5500/k^2 for every integer k>=1252",
        },
        "diagnostics": diagnostics,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": 4,
            "open_rows": 1,
            "selected_points": len(selected),
            "positive_coordinate_balls": len(selected) * 8,
            "point_curvature_cap_certificates": len(selected),
            "continuous_curvature_theorems": 0,
            "full_kernel_transfer_theorems": 0,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order10_first_summand_point_scout.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order10_first_summand_point_scout.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    diagnostics = artifact["diagnostics"]
    lines = [
        "# Order-Ten Seventh-Nested First-Summand Point Scout",
        "",
        "Date: 2026-07-16",
        "",
        "Status: rigorous selected-point first-summand scout. This is not a proof;",
        "it is not a",
        "continuous curvature theorem, a full-kernel tail theorem, or a proof",
        "of RH or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order10_first_summand_point_scout.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order10_first_summand_point_scout.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order10_first_summand_point_scout.py",
        "```",
        "",
        "## New Stable Layer",
        "",
        "```text",
        "W_1(t)=8*B_1(t)-w_1(t-1)+2*w_1(t)-w_1(t+1)",
        "z_1(t)=2*w_1(t)-s_1(t)+log(1-exp(-W_1(t)))",
        "```",
        "",
        "At each selected exact cache point",
        "`t in {1251,1300,1500,2000,3000,4000}`, all eight coordinates",
        "`B,J,R,S,T,U,V,W` are rigorously positive and",
        "",
        "```text",
        "t^2*z_1''(t)<250<5500.",
        "```",
        "",
        "The largest selected-point upper enclosure occurs at `t="
        + diagnostics["maximum_scaled_curvature_anchor"]
        + "` and is `"
        + diagnostics["maximum_scaled_curvature_upper"]
        + "`. The smallest selected `W` lower bound is `"
        + diagnostics["minimum_W_lower"]
        + "`.",
        "",
        "At `t=1251` the rigorous ball is",
        "",
        "```text",
        next(
            row["z1_scaled_second_ball"]
            for row in diagnostics["selected_points"]
            if row["t"] == "1251"
        ),
        "```",
        "",
        "## What Remains",
        "",
        "The headroom is substantial, but isolated points do not control a real",
        "interval. The next proof object must extend the localized tent hierarchy",
        "through `W,z`, cover `t>=1251`, and add a rigorous higher-summand",
        "transfer whose total stays below the full `5500/k^2` budget.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order10_m100_tail_curvature_reduction.md",
        "outputs/jensen_window_pf_compound_order9_first_summand_curvature_certificate.md",
        "outputs/jensen_window_pf_endpoint_order10_counterexample.md",
        "outputs/formal_core.md",
        "```",
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
        "wrote order-ten first-summand point scout: "
        f"{summary['selected_points']} selected points, "
        f"{summary['positive_coordinate_balls']} positive coordinates, "
        f"{summary['open_rows']} open continuum handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
