#!/usr/bin/env python3
"""Independently validate the Q31 margin geometry audit."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR_DIR = REPO_ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR_DIR) not in sys.path:
    sys.path.insert(0, str(VENDOR_DIR))

from flint import arb, ctx


STEM = "jensen_window_pf_newman_theta_q31_margin_geometry_audit"
DEFAULT_ARTIFACT = (
    REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
)
Q31_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_modular_retained_q31_interval_certificate.json"
)
Q31_CACHE = Q31_SOURCE.with_suffix(".jsonl")
ARBITRARY_N_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_arbitrary_n_stable_remainder_gate.json"
)
SADDLE_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_modular_blend_adaptive_saddle_gate.json"
)
SCHEMA = "newman_theta_modular_retained_q31_interval_v1"
ZERO_HASH = "0" * 64
EXPECTED_IDS = [
    "ntq31mga_01_cache_provenance",
    "ntq31mga_02_weakest_leaf_transversality",
    "ntq31mga_03_saddle_coordinate_jet",
    "ntq31mga_04_finite_saddle_jet_floor",
    "ntq31mga_05_raw_ratio_extrapolation_guard",
    "ntq31mga_06_adaptive_count_noncoverage",
    "ntq31mga_07_saddle_normalized_candidate",
]


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def fraction_arb(text: str) -> arb:
    value = Fraction(text)
    return arb(value.numerator) / value.denominator


def endpoint_abs_lower(lower_text: str, upper_text: str) -> arb:
    lower = arb(lower_text)
    upper = arb(upper_text)
    if lower.lower() > 0:
        return lower.lower()
    if upper.upper() < 0:
        return -upper.upper()
    return arb(0)


def sign_pattern(lower_text: str, upper_text: str) -> str:
    lower = arb(lower_text)
    upper = arb(upper_text)
    if lower.lower() > 0:
        return "positive"
    if upper.upper() < 0:
        return "negative"
    return "contains_zero"


def read_cache(config_hash: str, issues: list[str]) -> list[dict]:
    raw = Q31_CACHE.read_bytes()
    if raw and not raw.endswith(b"\n"):
        issues.append("Q31 cache lacks final line terminator")
    records: list[dict] = []
    previous = ZERO_HASH
    for line_number, line in enumerate(
        raw.decode("utf-8").splitlines(), start=1
    ):
        try:
            record = json.loads(line)
        except Exception as exc:
            issues.append(f"cache JSON failed at line {line_number}: {exc}")
            continue
        stored_hash = record.get("row_sha256")
        payload = dict(record)
        payload.pop("row_sha256", None)
        actual_hash = sha256(
            canonical_json(payload).encode("utf-8")
        ).hexdigest()
        if stored_hash != actual_hash:
            issues.append(f"cache hash mismatch at line {line_number}")
        if record.get("sequence") != line_number:
            issues.append(f"cache sequence mismatch at line {line_number}")
        if record.get("schema") != SCHEMA:
            issues.append(f"cache schema mismatch at line {line_number}")
        if record.get("config_sha256") != config_hash:
            issues.append(f"cache config mismatch at line {line_number}")
        if record.get("previous_row_sha256") != previous:
            issues.append(f"cache chain mismatch at line {line_number}")
        previous = stored_hash or ""
        records.append(record)
    return records


def recompute_leaf(sequence: int, leaf: dict) -> dict:
    value_lower = endpoint_abs_lower(
        leaf["j_retained_lower"], leaf["j_retained_upper"]
    )
    derivative_lower = endpoint_abs_lower(
        leaf["j_retained_prime_lower"],
        leaf["j_retained_prime_upper"],
    )
    scale = 4 * (arb.pi() * fraction_arb(leaf["x_low"])).sqrt()
    scaled_derivative = scale * derivative_lower
    saddle_lower = (
        value_lower
        if value_lower.lower() >= scaled_derivative.lower()
        else scaled_derivative
    )
    unscaled_lower = (
        value_lower
        if value_lower.lower() >= derivative_lower.lower()
        else derivative_lower
    )
    return {
        "sequence": sequence,
        "leaf": leaf,
        "value_lower": value_lower,
        "derivative_lower": derivative_lower,
        "scale": scale,
        "scaled_derivative": scaled_derivative,
        "saddle_lower": saddle_lower,
        "unscaled_lower": unscaled_lower,
        "value_sign": sign_pattern(
            leaf["j_retained_lower"], leaf["j_retained_upper"]
        ),
        "derivative_sign": sign_pattern(
            leaf["j_retained_prime_lower"],
            leaf["j_retained_prime_upper"],
        ),
    }


def validate(path: Path) -> list[str]:
    issues: list[str] = []
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("kind") != STEM:
        issues.append("artifact kind mismatch")
    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    if [row.get("readiness") for row in rows] != (
        ["proved"] * 6 + ["not_ready_to_apply"]
    ):
        issues.append("row readiness mismatch")

    source = json.loads(Q31_SOURCE.read_text(encoding="utf-8"))
    if not source.get("summary", {}).get("theorem_ready"):
        issues.append("Q31 source no longer theorem-ready")
    arbitrary = json.loads(
        ARBITRARY_N_SOURCE.read_text(encoding="utf-8")
    )
    if arbitrary.get("kind") != (
        "jensen_window_pf_newman_theta_arbitrary_n_stable_remainder_gate"
    ):
        issues.append("arbitrary-N source kind mismatch")
    audit = artifact.get("source_audit", {})
    if audit.get("q31_artifact_sha256") != file_hash(Q31_SOURCE):
        issues.append("Q31 artifact hash mismatch")
    if audit.get("q31_cache_sha256") != file_hash(Q31_CACHE):
        issues.append("Q31 cache hash mismatch")
    if audit.get("arbitrary_n_artifact_sha256") != file_hash(
        ARBITRARY_N_SOURCE
    ):
        issues.append("arbitrary-N artifact hash mismatch")
    if audit.get("adaptive_saddle_artifact_sha256") != file_hash(
        SADDLE_SOURCE
    ):
        issues.append("adaptive-saddle artifact hash mismatch")

    records = read_cache(source["config_sha256"], issues)
    if audit.get("q31_config_sha256") != source["config_sha256"]:
        issues.append("Q31 config hash mismatch")
    if records and audit.get("q31_last_row_sha256") != (
        records[-1]["row_sha256"]
    ):
        issues.append("Q31 last-row hash mismatch")

    ctx.prec = 192
    recomputed: list[dict] = []
    branch_counts = {"value": 0, "derivative": 0}
    value_patterns = {"positive": 0, "negative": 0, "contains_zero": 0}
    derivative_patterns = {
        "positive": 0,
        "negative": 0,
        "contains_zero": 0,
    }
    for record in records:
        if record["result"].get("unresolved_boxes") != 0:
            issues.append(
                f"unresolved source box at sequence {record['sequence']}"
            )
        for leaf in record["result"]["evaluations"]:
            if not leaf.get("certified"):
                # Subdivided parent evaluations are retained in the cache.
                continue
            row = recompute_leaf(record["sequence"], leaf)
            recomputed.append(row)
            branch_counts[leaf["branch"]] += 1
            value_patterns[row["value_sign"]] += 1
            derivative_patterns[row["derivative_sign"]] += 1

    geometry = artifact.get("geometry", {})
    if geometry.get("records") != len(records):
        issues.append("record count mismatch")
    if geometry.get("leaves") != len(recomputed):
        issues.append("leaf count mismatch")
    if branch_counts != source["summary"]["branch_counts"]:
        issues.append("source branch count mismatch")
    if geometry.get("branch_counts") != branch_counts:
        issues.append("stored branch count mismatch")
    if geometry.get("value_sign_patterns") != value_patterns:
        issues.append("stored value sign patterns mismatch")
    if geometry.get("derivative_sign_patterns") != derivative_patterns:
        issues.append("stored derivative sign patterns mismatch")

    if recomputed:
        weakest_ratio = min(
            recomputed,
            key=lambda row: float(
                arb(row["leaf"]["certified_ratio_lower"]).lower()
            ),
        )
        weakest_unscaled = min(
            recomputed,
            key=lambda row: float(row["unscaled_lower"].lower()),
        )
        weakest_saddle = min(
            recomputed,
            key=lambda row: float(row["saddle_lower"].lower()),
        )
        stored_weakest = geometry["weakest_certified_ratio_leaves"][0]
        if stored_weakest["source_sequence"] != weakest_ratio["sequence"]:
            issues.append("weakest ratio source sequence mismatch")
        for key in ("t_low", "t_high", "x_low", "x_high", "branch"):
            if stored_weakest[key] != weakest_ratio["leaf"][key]:
                issues.append(f"weakest ratio {key} mismatch")
        if weakest_ratio["value_sign"] != "contains_zero":
            issues.append("weakest ratio no longer crosses retained zero")
        if weakest_ratio["derivative_sign"] != "positive":
            issues.append("weakest ratio derivative is not positive")

        stored_unscaled = arb(
            geometry["minimum_unscaled_first_jet_lower"]
        )
        if not stored_unscaled.overlaps(weakest_unscaled["unscaled_lower"]):
            issues.append("minimum unscaled jet mismatch")
        stored_saddle = arb(geometry["minimum_saddle_first_jet_lower"])
        if not stored_saddle.overlaps(weakest_saddle["saddle_lower"]):
            issues.append("minimum saddle jet mismatch")
        stored_saddle_leaf = geometry["weakest_saddle_first_jet_leaf"]
        if stored_saddle_leaf["source_sequence"] != weakest_saddle["sequence"]:
            issues.append("weakest saddle leaf sequence mismatch")

    coordinate = geometry.get("saddle_coordinate", {})
    if coordinate.get("definition") != "s=sqrt(x/(4*pi))":
        issues.append("saddle coordinate mismatch")
    if coordinate.get("derivative_conversion") != (
        "partial_s=4*sqrt(pi*x)*partial_x"
    ):
        issues.append("saddle derivative conversion mismatch")
    coordinate_marker = 16 * arb.pi()
    marker_rows = [
        row
        for row in recomputed
        if fraction_arb(row["leaf"]["x_low"]).lower() <= (
            coordinate_marker.upper()
        )
        and fraction_arb(row["leaf"]["x_high"]).upper() >= (
            coordinate_marker.lower()
        )
    ]
    if coordinate.get("integer_marker_leaf_count") != len(marker_rows):
        issues.append("integer-marker leaf count mismatch")
    if "only the leading asymptotic marker" not in (
        coordinate.get("exact_transition_surfaces", "")
    ):
        issues.append("exact saddle-transition guard missing")

    count_audit = geometry.get("adaptive_count_audit", {})
    if not (15**4 < 39**3 < 16**4):
        issues.append("x=38 adaptive-count arithmetic failed")
    if not (24**4 < 70**3 < 25**4):
        issues.append("x=69 adaptive-count arithmetic failed")
    expected_count = {
        "retained_count_in_Q31": 7,
        "kappa_one_count_at_x38": 16,
        "kappa_one_count_at_x69": 25,
    }
    for key, value in expected_count.items():
        if count_audit.get(key) != value:
            issues.append(f"adaptive count mismatch for {key}")

    candidate = artifact.get("candidate", "")
    for marker in (
        "s=sqrt(x/(4*pi))",
        "4*sqrt(pi*x)*|J_N'|",
        "A_N(t,x)",
        "E_(0,N)",
        "E_(1,N)",
        "n_*(a,t,x)=k",
        "a in {5,9}",
        "x=4*pi*k^2",
        "only the leading marker",
        "adaptive-count jump",
    ):
        if marker not in candidate:
            issues.append(f"candidate marker missing: {marker}")
    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "switch-defect",
        "amplitude",
        "x>69",
        "Lambda<=0",
        "RH",
        "Clay-prize",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    args = parser.parse_args()
    issues = validate(args.artifact)
    if issues:
        for issue in issues:
            print(f"ISSUE: {issue}")
        raise SystemExit(1)
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    geometry = artifact["geometry"]
    print(
        "validated Q31 margin geometry audit: "
        f"{geometry['records']} cache records, {geometry['leaves']} leaves, "
        f"{len(artifact['rows'])} rows, 0 issues, "
        "1 saddle-normalized cofinal handoff"
    )


if __name__ == "__main__":
    main()
