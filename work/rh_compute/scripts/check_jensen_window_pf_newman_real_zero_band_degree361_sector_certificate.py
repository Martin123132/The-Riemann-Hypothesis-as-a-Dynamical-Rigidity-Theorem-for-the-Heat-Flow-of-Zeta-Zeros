#!/usr/bin/env python3
"""Check the Newman real-zero-band degree-361 sector certificate."""

from __future__ import annotations

from fractions import Fraction
import importlib.util
import json
import math
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_newman_real_zero_band_"
    "degree361_sector_certificate.py"
)
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_real_zero_band_"
    "degree361_sector_certificate.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_newman_real_zero_band_"
    "degree361_sector_certificate.md"
)


def load_builder():
    spec = importlib.util.spec_from_file_location(
        "newman_real_zero_band_degree361",
        BUILDER,
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load degree-361 builder")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def check_partition(payload: dict, issues: list[str]) -> None:
    certificate = payload["exact"]["boundary_certificate"]
    records = certificate.get("records", [])
    parsed: list[
        tuple[Fraction, Fraction, Fraction, Fraction]
    ] = []
    for index, row in enumerate(records):
        try:
            box = (
                Fraction(row["t_low"]),
                Fraction(row["t_high"]),
                Fraction(row["y_low"]),
                Fraction(row["y_high"]),
            )
        except (KeyError, ValueError, ZeroDivisionError) as exc:
            issues.append(f"invalid box {index}: {exc}")
            continue
        t_low, t_high, y_low, y_high = box
        if not (
            Fraction(0) <= t_low < t_high <= Fraction(1, 5)
        ):
            issues.append(f"box {index} leaves the time interval")
        if not (
            Fraction(0) <= y_low < y_high <= Fraction(1)
        ):
            issues.append(f"box {index} leaves the imaginary interval")
        if row.get("certified") is not True:
            issues.append(f"box {index} is not certified")
        parsed.append(box)
    area = sum(
        (t_high - t_low) * (y_high - y_low)
        for t_low, t_high, y_low, y_high in parsed
    )
    if area != Fraction(1, 5):
        issues.append(f"partition area mismatch: {area}")
    if len(records) != 200:
        issues.append(f"expected 200 boundary boxes, got {len(records)}")


def validate(payload: dict, note: str, builder) -> list[str]:
    issues: list[str] = []
    if payload.get("kind") != (
        "jensen_window_pf_newman_real_zero_band_"
        "degree361_sector_certificate"
    ):
        issues.append("artifact kind changed")
    if payload.get("status") != (
        "rigorous real-zero-band homotopy and "
        "published-sector composition theorem"
    ):
        issues.append("artifact status changed")
    rows = payload.get("rows", [])
    expected_ids = [
        f"nrzb361_{index:02d}_{suffix}"
        for index, suffix in (
            (1, "theta_enclosure"),
            (2, "taylor_enclosure"),
            (3, "vertical_boundary"),
            (4, "horizontal_boundary"),
            (5, "top_time_reality"),
            (6, "no_collision"),
            (7, "real_zero_band"),
            (8, "sector_transfer"),
            (9, "chasse_transfer"),
            (10, "degree361"),
            (11, "degree71_sharpening"),
            (12, "nonpromotion"),
        )
    ]
    if [row.get("id") for row in rows] != expected_ids:
        issues.append("row ids or order changed")
    if any(row.get("readiness") != "ready_to_apply" for row in rows[:-1]):
        issues.append("a proved row is not ready_to_apply")
    if not rows or rows[-1].get("readiness") != "open":
        issues.append("finite proof boundary was promoted")

    exact = payload.get("exact", {})
    certificate = exact.get("boundary_certificate", {})
    expected_certificate = {
        "precision_bits": 192,
        "time_taylor_order": 2,
        "imaginary_taylor_order": 4,
        "time_step": "1/50",
        "imaginary_step": "1/20",
        "integration_cutoff": "2",
        "retained_theta_terms": 4,
        "far_tail_radius": "1e-800",
        "boxes": 200,
        "unresolved_boxes": 0,
    }
    for key, expected in expected_certificate.items():
        if certificate.get(key) != expected:
            issues.append(
                f"certificate {key} changed: "
                f"{certificate.get(key)!r} != {expected!r}"
            )
    try:
        weak = builder.arb(
            certificate["weakest_negative_margin_lower"]
        )
        if not weak > 0:
            issues.append("weakest negative margin is not positive")
        for index, row in enumerate(certificate.get("records", [])):
            upper = builder.arb(row["upper"])
            margin = builder.arb(row["negative_margin_lower"])
            if not upper < 0 or not margin > 0:
                issues.append(f"box {index} lost strict negativity")
    except (KeyError, ValueError) as exc:
        issues.append(f"invalid interval endpoint: {exc}")
    check_partition(payload, issues)

    far_audit = exact.get("kernel", {}).get("far_tail_audit", {})
    exp8_lower = sum(
        Fraction(8**power, math.factorial(power))
        for power in range(8)
    )
    if exp8_lower != Fraction(47259, 35):
        issues.append("independent exp(8) partial sum changed")
    endpoint_upper = (
        Fraction(4, 5) + 29 - 3 * exp8_lower
    )
    if endpoint_upper != Fraction(-140734, 35):
        issues.append("far-tail endpoint exponent audit failed")
    derivative_upper = (
        Fraction(4, 5)
        + Fraction(29, 2)
        - 12 * exp8_lower
    )
    if not derivative_upper < -1:
        issues.append("far-tail exponent derivative is not below -1")
    if not Fraction(22, 7) ** 2 < 10:
        issues.append("independent pi-squared upper audit failed")
    if not 40 * 10**800 < 2**4020:
        issues.append("far-tail integer comparison failed")
    expected_far_audit = {
        "moment_range": "0<=m<=9",
        "exp8_partial_sum": "47259/35",
        "exponent": "u^2/5+(29/2)u-pi*exp(4u)",
        "endpoint_upper": "-140734/35",
        "derivative_rule": (
            "the exponent derivative is below -1 for u>=2"
        ),
        "prefactor_upper": "40",
        "integer_comparison": "40*10^800<2^4020",
        "conclusion": "40*exp(-140734/35)<10^-800",
    }
    if far_audit != expected_far_audit:
        issues.append("stored far-tail analytic audit changed")

    if "Re H_t(38+iy)<0" not in exact.get(
        "vertical_boundary_theorem",
        "",
    ):
        issues.append("vertical-boundary theorem missing")
    homotopy = exact.get("homotopy", {})
    for phrase in (
        "Boundary zero-freeness fixes",
        "leaving the real axis would require a real multiple zero",
    ):
        if phrase not in homotopy.get("argument", ""):
            issues.append(f"homotopy guard missing: {phrase}")
    if "imply Im z=0" not in homotopy.get("theorem", ""):
        issues.append("real-zero-band theorem missing")

    sector = exact.get("sector_transfer", {})
    sine = Fraction(sector.get("sector_sine", "0"))
    inverse = Fraction(sector.get("inverse_sine_square", "0"))
    if sine != Fraction(76, 1445):
        issues.append("sector sine identity failed")
    if inverse != Fraction(2088025, 5776):
        issues.append("inverse sine-square identity failed")
    if inverse.numerator // inverse.denominator != 361:
        issues.append("sector floor is not 361")
    if sector.get("degree_cutoff") != 361:
        issues.append("stored degree cutoff changed")
    if "every 0<=d<=361" not in sector.get("theorem", ""):
        issues.append("degree-361 theorem statement missing")

    sources = payload.get("sources", {})
    for source_name, source_path in (
        ("compact_no_collision", builder.CONTACT_SOURCE),
        ("degree71_sector", builder.DEGREE71_SOURCE),
    ):
        source = sources.get(source_name, {})
        if source.get("sha256") != builder.sha256(source_path):
            issues.append(f"{source_name} source hash changed")
    if sources.get("debruijn_strip", {}).get("theorem") != "Theorem 3.2":
        issues.append("de Bruijn theorem input changed")
    if sources.get("lambda_upper", {}).get("statement") != "Lambda<=1/5":
        issues.append("Lambda upper input changed")
    if sources.get("chasse_sector", {}).get("theorem") != (
        "Theorem 3.6 and derivative closure"
    ):
        issues.append("Chasse input changed")

    summary = payload.get("summary", {})
    expected_summary = {
        "rows": 12,
        "precision_bits": 192,
        "boundary_boxes": 200,
        "unresolved_boxes": 0,
        "vertical_zero_free_boundaries": 2,
        "real_zero_bands": 1,
        "all_shifts": True,
        "maximum_proved_degree": 361,
        "previous_maximum_degree": 71,
        "all_degree_theorems": 0,
    }
    if summary != expected_summary:
        issues.append("summary counts changed")

    for marker in (
        "# Newman Real-Zero-Band Degree-361 Sector Certificate",
        "Re H_t(38+iy)<0",
        "(H_t(x),H_t'(x))!=(0,0)",
        "H_t(z)=0 and |Re z|<=38 imply Im z=0",
        "sin(delta)=76/1445",
        "|sin(delta)|^(-2)=2088025/5776",
        "40*10^800<2^4020",
        "for every 0<=d<=361, n>=0, and 0<=t<=1/5",
        "The cutoff `361` is finite",
        "Degree `362`",
        "`Lambda<=0`",
        "RH",
    ):
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    return issues


def main() -> int:
    builder = load_builder()
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")
    issues = validate(payload, note, builder)
    rebuilt = builder.build_payload()
    rebuilt["priority_lowered"] = payload.get("priority_lowered")
    if rebuilt != payload:
        issues.append("stored payload differs from full Arb reconstruction")
    if builder.success_line(payload) not in note:
        issues.append("success line missing from note")
    if issues:
        print(
            "Newman real-zero-band degree-361 certificate: "
            f"{len(issues)} issues"
        )
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(builder.success_line(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
