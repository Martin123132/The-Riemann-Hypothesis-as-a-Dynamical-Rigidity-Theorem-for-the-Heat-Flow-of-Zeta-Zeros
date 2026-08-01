#!/usr/bin/env python3
"""Check the Newman zero-slab degree-71 sector certificate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_zero_slab_degree71_sector_certificate.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_newman_zero_slab_degree71_sector_certificate.md"
)
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_newman_zero_slab_degree71_sector_certificate.py"
)


def load_builder():
    spec = importlib.util.spec_from_file_location(
        "newman_zero_slab_degree71",
        BUILDER,
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load zero-slab builder")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=RESULT)
    parser.add_argument("--note", type=Path, default=NOTE)
    args = parser.parse_args()
    issues: list[str] = []
    payload = json.loads(args.result.read_text(encoding="utf-8"))
    builder = load_builder()
    rebuilt = builder.build_payload()

    # Priority calls can report differently on non-Windows replays.
    rebuilt["priority_lowered"] = payload.get("priority_lowered")
    if payload != rebuilt:
        issues.append("stored payload differs from rigorous reconstruction")

    if payload.get("kind") != (
        "jensen_window_pf_newman_zero_slab_"
        "degree71_sector_certificate"
    ):
        issues.append("artifact kind changed")
    if payload.get("status") != (
        "rigorous uniform zero-slab and published-sector "
        "composition theorem"
    ):
        issues.append("artifact status changed")
    if len(payload.get("rows", [])) != 13:
        issues.append("expected thirteen theorem rows")

    exact = payload.get("exact", {})
    constants = exact.get("constants", {})
    if constants.get("precision_bits") != 192:
        issues.append("precision changed")
    if constants.get("time_interval") != "0<=t<=1/5":
        issues.append("heat interval changed")
    if constants.get("slab") != "|Re z|<=84/5, |Im z|<=1":
        issues.append("complex slab changed")
    if constants.get("split") != "5*pi/168":
        issues.append("directed split changed")

    bounds = exact.get("kernel_bounds", {})
    for label in ("higher_ratio_ball", "full_ratio_ball"):
        row = bounds.get(label, {})
        if not row.get("strictly_positive"):
            issues.append(f"{label} lost positivity")
        try:
            upper = builder.arb(row.get("upper", "2"))
        except Exception:
            issues.append(f"{label} upper endpoint is not parseable")
        else:
            if not upper < 1:
                issues.append(f"{label} is not below one")
    if not bounds.get("far_endpoint_gap", {}).get(
        "strictly_positive"
    ):
        issues.append("far exponent endpoint gap is not positive")

    integrals = exact.get("directed_integrals", {})
    for label in (
        "core",
        "first_tail",
        "higher_tail",
        "infinite_tail",
        "margin",
    ):
        if label not in integrals:
            issues.append(f"missing directed integral: {label}")
    if not integrals.get("margin", {}).get("strictly_positive"):
        issues.append("directed slab margin is not positive")
    try:
        core_lower = builder.arb(integrals["core"]["lower"])
        tail_upper = sum(
            (
                builder.arb(integrals[label]["upper"])
                for label in (
                    "first_tail",
                    "higher_tail",
                    "infinite_tail",
                )
            ),
            builder.arb(0),
        )
    except Exception as exc:
        issues.append(f"directed endpoint reconstruction failed: {exc}")
    else:
        if not core_lower - tail_upper > 0:
            issues.append("independent lower-minus-upper margin failed")

    zero_slab = exact.get("zero_slab", {})
    if zero_slab.get("zero_free") is not True:
        issues.append("zero-free slab conclusion changed")
    if "Re H_t(x+iy)>0" not in zero_slab.get("theorem", ""):
        issues.append("positive-real-part theorem missing")

    published = exact.get("published_inputs", {})
    if published.get("debruijn_strip", {}).get("theorem") != (
        "Theorem 3.2"
    ):
        issues.append("de Bruijn source theorem changed")
    if published.get("chasse_sector", {}).get("theorem") != (
        "Theorem 3.6 and derivative closure"
    ):
        issues.append("Chasse source theorem changed")
    if published.get("endpoint_comparison", {}).get("result") != (
        "Corollary 1.3"
    ):
        issues.append("endpoint comparison changed")

    sector = exact.get("sector_transfer", {})
    sine = Fraction(sector.get("sector_sine", "0"))
    inverse_square = Fraction(
        sector.get("inverse_sine_square", "0")
    )
    if sine != Fraction(840, 7081):
        issues.append("sector sine identity failed")
    if inverse_square != Fraction(50140561, 705600):
        issues.append("inverse sine-square identity failed")
    if inverse_square.numerator // inverse_square.denominator != 71:
        issues.append("degree floor is not 71")
    if sector.get("degree_cutoff") != 71:
        issues.append("stored degree cutoff changed")
    if "every 0<=d<=71" not in sector.get("theorem", ""):
        issues.append("degree-71 theorem statement missing")

    summary = payload.get("summary", {})
    expected_summary = {
        "rows": 13,
        "precision_bits": 192,
        "directed_integrals": 4,
        "positive_margin": 1,
        "zero_free_complex_slabs": 1,
        "published_inputs": 3,
        "uniform_heat_parameters": "continuum_0_to_1_over_5",
        "all_shifts": True,
        "maximum_proved_degree": 71,
        "quartic_layers_closed": 1,
        "quintic_layers_closed": 1,
        "all_degree_theorems": 0,
    }
    if summary != expected_summary:
        issues.append("summary counts changed")

    row_ids = {row.get("id") for row in payload.get("rows", [])}
    expected_ids = {
        f"nzsd71_{index:02d}_{suffix}"
        for index, suffix in (
            (1, "kernel_positivity"),
            (2, "near_lower_bound"),
            (3, "theta_tail_envelope"),
            (4, "far_tail_envelope"),
            (5, "directed_integrals"),
            (6, "zero_free_slab"),
            (7, "debruijn_strip"),
            (8, "squared_zero_sector"),
            (9, "chasse_transfer"),
            (10, "degree71_theorem"),
            (11, "quartic_consequence"),
            (12, "endpoint_scope"),
            (13, "nonpromotion"),
        )
    }
    if row_ids != expected_ids:
        issues.append("theorem row ids changed")

    note = args.note.read_text(encoding="utf-8")
    for marker in (
        "# Newman Zero-Slab Degree-71 Sector Certificate",
        "rigorous interval certificate and published-theorem composition",
        "`a=5*pi/168`",
        "Re H_t(x+iy)>0",
        "sin(delta)=840/7081",
        "|sin(delta)|^(-2)=50140561/705600",
        "for every 0<=d<=71, n>=0, and 0<=t<=1/5",
        "Degrees four and five are now closed",
        "d<=9.36*10^20",
        "The cutoff `71` is finite",
        "`Lambda<=0`",
        "RH",
    ):
        if marker not in note:
            issues.append(f"note marker missing: {marker}")

    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated Newman zero-slab degree-71 sector certificate: "
        "13 rows, 0 issues, 4 directed integrals, "
        "1 positive slab margin, 1 zero-free complex slab, "
        "3 published inputs, all shifts through degree 71, "
        "0 all-degree theorems"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
