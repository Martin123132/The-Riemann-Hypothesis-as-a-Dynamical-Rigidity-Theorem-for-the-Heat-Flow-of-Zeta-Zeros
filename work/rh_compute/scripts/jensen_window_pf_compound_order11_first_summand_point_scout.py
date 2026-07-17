#!/usr/bin/env python3
"""Scout the eighth-nested first-summand coordinate at exact points."""

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


RESULTS = REPO_ROOT / "work/rh_compute/results"
DEFAULT_OUT = RESULTS / "jensen_window_pf_compound_order11_first_summand_point_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_compound_order11_first_summand_point_scout.md"
POINT_CACHE = RESULTS / "jensen_window_pf_compound_order9_shifted_point_h0_h8_half_grid.jsonl"
POINT_MANIFEST = RESULTS / "jensen_window_pf_compound_order9_shifted_point_h0_h8_cache.json"
ORDER10_REDUCTION = RESULTS / "jensen_window_pf_compound_order10_m100_tail_curvature_reduction.json"

SELECTED_ANCHORS = (
    Fraction(1252),
    Fraction(1300),
    Fraction(1500),
    Fraction(1700),
    Fraction(1800),
    Fraction(1900),
    Fraction(1950),
    Fraction(2000),
)
POINT_JET_ORDER = 4
POINT_CURVATURE_CAP = 1000
PRECISION_BITS = 512
GAP_NAMES = {2: "J", 3: "R", 4: "S", 5: "T", 6: "U", 7: "V", 8: "W", 9: "X"}


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
    reduction = json.loads(ORDER10_REDUCTION.read_text(encoding="utf-8"))
    cache = manifest.get("cache", {})
    parameters = manifest.get("parameters", {})
    if manifest.get("kind") != "jensen_window_pf_compound_order9_shifted_point_h0_h8_cache":
        raise RuntimeError("point-cache manifest kind changed")
    if cache.get("row_count") != 8929 or cache.get("all_rows_passed") is not True:
        raise RuntimeError("point-cache row contract changed")
    if cache.get("h_derivative_orders") != [0, 8] or cache.get("sha256") != sha256(POINT_CACHE):
        raise RuntimeError("point-cache derivative or hash contract changed")
    expected_parameters = {
        "start_t": "1243",
        "end_t": "5707",
        "step_t": "1/2",
        "precision_bits": PRECISION_BITS,
    }
    for key, expected in expected_parameters.items():
        if parameters.get(key) != expected:
            raise RuntimeError(f"point-cache parameter changed at {key}")
    exact = reduction.get("exact", {})
    if exact.get("canonical_factorization") != "Q_(9,n)=A_(n+8)^9*exp(z(n+8))":
        raise RuntimeError("order-ten source lost the z coordinate")
    if exact.get("seventh_gap") != "W(t)=8*B(t)-w(t-1)+2*w(t)-w(t+1)":
        raise RuntimeError("order-ten source W gap changed")
    return {
        "point_manifest": POINT_MANIFEST.relative_to(REPO_ROOT).as_posix(),
        "point_manifest_sha256": sha256(POINT_MANIFEST),
        "point_cache": POINT_CACHE.relative_to(REPO_ROOT).as_posix(),
        "point_cache_sha256": cache["sha256"],
        "point_cache_rows": cache["row_count"],
        "point_cache_range": [parameters["start_t"], parameters["end_t"]],
        "point_cache_step": parameters["step_t"],
        "point_cache_precision_bits": parameters["precision_bits"],
        "order10_reduction": ORDER10_REDUCTION.relative_to(REPO_ROOT).as_posix(),
        "order10_reduction_sha256": sha256(ORDER10_REDUCTION),
    }


def needed_targets() -> set[Fraction]:
    return {
        anchor + shift
        for anchor in SELECTED_ANCHORS
        for shift in range(-9, 10)
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
    source = {shift: h_jets[anchor + shift] for shift in range(-9, 10)}
    b = {
        shift: centered_second_difference(source, shift)
        for shift in range(-8, 9)
    }
    zero = [flint.arb(0) for _ in range(order + 1)]
    logs: dict[int, dict[int, list[flint.arb]]] = {
        0: {shift: zero[:] for shift in range(-8, 9)},
        1: {
            shift: stable_log_series(value, order)
            for shift, value in b.items()
        },
    }
    gaps: dict[int, dict[int, list[flint.arb]]] = {}
    for stage in range(2, 10):
        radius = 9 - stage
        stage_gaps = {
            shift: series_sub(
                series_scale(b[shift], stage),
                centered_second_difference(logs[stage - 1], shift),
            )
            for shift in range(-radius, radius + 1)
        }
        stage_logs = {}
        for shift in range(-radius, radius + 1):
            try:
                stable_log = stable_log_series(stage_gaps[shift], order)
            except ValueError as exc:
                raise RuntimeError(
                    f"stable log failed at t={anchor}, stage={stage}, "
                    f"shift={shift}, gap0={stage_gaps[shift][0]}"
                ) from exc
            stage_logs[shift] = series_add(
                series_sub(
                    series_scale(logs[stage - 1][shift], 2),
                    logs[stage - 2][shift],
                ),
                stable_log,
            )
        gaps[stage] = stage_gaps
        logs[stage] = stage_logs

    coordinates = {"B": b[0][0]}
    coordinates.update({GAP_NAMES[stage]: gaps[stage][0][0] for stage in range(2, 10)})
    if not all(bool(value > 0) for value in coordinates.values()):
        raise RuntimeError(f"selected point lost a positive coordinate at t={anchor}")

    anchor_arb = flint.arb(anchor.numerator) / anchor.denominator
    z = logs[8][0]
    y = logs[9][0]
    z_scaled_second = 2 * z[2] * anchor_arb**2
    y_scaled_first = y[1] * anchor_arb
    y_scaled_second = 2 * y[2] * anchor_arb**2
    cap_margin = flint.arb(POINT_CURVATURE_CAP) - y_scaled_second
    if not bool(cap_margin > 0):
        raise RuntimeError(f"selected point exceeds scout cap at t={anchor}")
    return {
        "t": str(anchor),
        "coordinate_balls": {name: ball_text(value) for name, value in coordinates.items()},
        "coordinate_lower_bounds": {name: arb_lower_text(value) for name, value in coordinates.items()},
        "z1_scaled_second_ball": ball_text(z_scaled_second),
        "y1_scaled_first_ball": ball_text(y_scaled_first),
        "y1_scaled_second_ball": ball_text(y_scaled_second),
        "y1_scaled_second_upper": arb_upper_text(y_scaled_second),
        "point_cap": POINT_CURVATURE_CAP,
        "point_cap_margin_lower": arb_lower_text(cap_margin),
        "passed": True,
        "_X": coordinates["X"],
        "_y_scaled_second": y_scaled_second,
    }


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    source_contract = load_source_contract()
    h_jets = load_h_jets()
    calculated = [exact_point_hierarchy(anchor, h_jets) for anchor in SELECTED_ANCHORS]
    maximum = max(calculated, key=lambda row: row["_y_scaled_second"].upper())
    minimum_x = min(calculated, key=lambda row: row["_X"].lower())
    selected = [
        {key: value for key, value in row.items() if not key.startswith("_")}
        for row in calculated
    ]
    diagnostics = {
        "selected_points": selected,
        "maximum_scaled_curvature_anchor": maximum["t"],
        "maximum_scaled_curvature_upper": maximum["y1_scaled_second_upper"],
        "minimum_X_anchor": minimum_x["t"],
        "minimum_X_lower": arb_lower_text(minimum_x["_X"]),
        "point_curvature_cap": POINT_CURVATURE_CAP,
    }
    rows = [
        ScoutRow(
            "co11fsp_01_exact_point_cache",
            "interval_input",
            "ready_to_apply",
            "The reusable exact-point cache supplies rigorous first-summand H jets through derivative order eight.",
            "H_1^(j)(t) enclosed for j=0,...,8 on the exact half-grid 1243<=t<=5707",
            "Computational first-summand input only.",
        ),
        ScoutRow(
            "co11fsp_02_generic_hierarchy",
            "exact_interval_algebra",
            "ready_to_apply",
            "The condensation recursion is rebuilt uniformly through stages two to nine.",
            "G_p=p*B-Delta^2*l_(p-1); l_p=2*l_(p-1)-l_(p-2)+log(1-exp(-G_p))",
            "Exact Taylor-series algebra at selected points.",
        ),
        ScoutRow(
            "co11fsp_03_eighth_nested_stage",
            "exact_interval_algebra",
            "ready_to_apply",
            "One stage beyond order ten defines the X gap and y coordinate.",
            "X=9*B-Delta^2*z; y=2*z-w+log(1-exp(-X))",
            "First Newman summand only.",
        ),
        ScoutRow(
            "co11fsp_04_positive_coordinates",
            "rigorous_point_scout",
            "ready_to_apply",
            "Every selected point has strict positive B,J,R,S,T,U,V,W,X enclosures.",
            "B,J,R,S,T,U,V,W,X>0 at all selected t",
            "Selected exact points, not intervals between them.",
            diagnostics,
        ),
        ScoutRow(
            "co11fsp_05_curvature_scale",
            "rigorous_point_scout",
            "ready_to_apply",
            "The selected-point enclosures measure the eighth-nested scaled curvature without finite differencing.",
            "t^2*y_1''(t) enclosed directly from fourth-order Arb Taylor jets",
            "This chooses a future target; it is not a continuous curvature theorem.",
            diagnostics,
        ),
        ScoutRow(
            "co11fsp_06_continuum_handoff",
            "analytic_handoff",
            "not_ready_to_apply",
            "Upgrade the point scale to a continuous first-summand theorem and control the full-kernel transfer.",
            "prove y_1''(t)<=C_11/t^2 and |Y_k-Y_k^(1)|<=D_11/k^2 on an explicit tail",
            "Requires a new localized eighth-stage enclosure, including a stable X floor beyond the dependency-limited point range, and a dominance transfer.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order11_first_summand_point_scout",
        "date": "2026-07-16",
        "status": "rigorous selected-point eighth-nested first-summand scout",
        "proof_boundary": (
            "This rigorously evaluates the first-summand hierarchy through X and "
            "y at eight exact points through t=2000. The reused point balls lose "
            "X sign separation farther out from dependency inflation; this does "
            "not assert that X changes sign. The scout does not cover intervals, "
            "prove a full-kernel order-eleven tail, PF-infinity, RH, or Lambda<=0."
        ),
        "source_contract": source_contract,
        "exact": {
            "generic_stage": "G_p=p*B-Delta^2*l_(p-1); l_p=2*l_(p-1)-l_(p-2)+log(1-exp(-G_p))",
            "eighth_gap": "X=9*B-Delta^2*z",
            "order10_log_coordinate": "y=2*z-w+log(1-exp(-X))",
            "scaled_curvature": "t^2*y_1''(t)",
        },
        "diagnostics": diagnostics,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": len(rows) - 1,
            "open_rows": 1,
            "selected_points": len(selected),
            "positive_coordinate_balls": len(selected) * 9,
            "point_curvature_enclosures": len(selected),
            "continuous_curvature_theorems": 0,
            "full_kernel_transfer_theorems": 0,
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order11_first_summand_point_scout.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order11_first_summand_point_scout.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    diagnostics = artifact["diagnostics"]
    lines = [
        "# Order-Eleven First-Summand Point Scout",
        "",
        "Date: 2026-07-16",
        "",
        "Status: rigorous selected-point eighth-nested first-summand scout. This is not a proof of RH or `Lambda <= 0`.",
        "",
        "```text",
        "X=9*B-Delta^2*z",
        "y=2*z-w+log(1-exp(-X))",
        "scaled curvature=t^2*y_1''(t)",
        "```",
        "",
        "| t | lower X | rigorous upper for t^2 y_1''(t) |",
        "|---:|---:|---:|",
    ]
    for row in diagnostics["selected_points"]:
        lines.append(
            f"| {row['t']} | `{row['coordinate_lower_bounds']['X']}` | "
            f"`{row['y1_scaled_second_upper']}` |"
        )
    lines.extend(
        [
            "",
            f"The largest selected upper occurs at `t={diagnostics['maximum_scaled_curvature_anchor']}`:",
            f"`{diagnostics['maximum_scaled_curvature_upper']}`.",
            "",
            "These are exact-point Arb enclosures, not a continuous ray theorem.",
            "With the reused 512-bit source balls, the next stage loses X sign",
            "separation beyond this range through interval dependency inflation;",
            "that is an enclosure limitation, not a certified sign change.",
            "The next analytic task is a localized eighth-stage continuum bound",
            "plus a first-summand-to-full-kernel transfer. This is not PF-infinity,",
            "RH, or `Lambda<=0`.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    maximum = artifact["diagnostics"]["maximum_scaled_curvature_upper"]
    print(f"wrote order-eleven first-summand point scout: 8 points, maximum scaled upper {maximum}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
