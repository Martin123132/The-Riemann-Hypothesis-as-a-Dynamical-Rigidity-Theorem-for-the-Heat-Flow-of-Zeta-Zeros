#!/usr/bin/env python3
"""Independently validate the order-eleven lambda-zero prefix certificate."""

from __future__ import annotations

import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order11_lambda0_prefix_certificate as generator  # noqa: E402
import jensen_window_pf_compound_order10_lambda0_prefix_certificate as order10  # noqa: E402


def main() -> int:
    artifact = json.loads(generator.DEFAULT_OUT.read_text(encoding="utf-8"))
    issues = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order11_lambda0_prefix_certificate"
    ):
        issues.append("bad artifact kind")
    expected_paths = (
        *generator.COEFFICIENT_SOURCES,
        generator.HANKEL_ROWS_SOURCE,
        generator.HANKEL_SUMMARY_SOURCE,
        generator.DEFINITION_SOURCE,
    )
    actual_sources = {
        row.get("path"): row.get("sha256") for row in artifact.get("sources", [])
    }
    for path in expected_paths:
        relative = order10.relative_path(path)
        if actual_sources.get(relative) != order10.sha256(path):
            issues.append(f"source hash mismatch: {relative}")
    if set(actual_sources) != {order10.relative_path(path) for path in expected_paths}:
        issues.append("source path set changed")

    expected_normalization = {
        "coefficient": "A_k(lambda)=mu_(2k)(lambda)*k!/(2k)!",
        "definition": (
            "H_(m,n)=det[A_(n+i+j)]_(0<=i,j<m), H_(0,n)=1, "
            "epsilon_m=(-1)^binom(m,2), Q_(m,n)=epsilon_m*H_(m,n)"
        ),
        "order": 11,
        "probe_parameter": 10,
        "probe_matrix_dimension": 11,
        "epsilon_11": -1,
        "probe_sigma": -1,
        "index_identification": "probe m=10 and shift=n equals Q_(11,n)",
    }
    if artifact.get("normalization") != expected_normalization:
        issues.append("normalization block changed")

    flint.ctx.dps = generator.PRECISION_DPS
    coefficients = generator.load_coefficients()
    source_rows = generator.load_probe_rows()
    rows = {int(row["n"]): row for row in artifact.get("finite", {}).get("rows", [])}
    if set(coefficients) != set(range(24)):
        issues.append("coefficient coverage is not A_0..A_23")
    if set(source_rows) != set(range(4)) or set(rows) != set(range(4)):
        issues.append("prefix row coverage is not n=0..3")
    if not issues:
        for n in range(4):
            raw = flint.arb_mat(
                [
                    [coefficients[n + i + j] for j in range(11)]
                    for i in range(11)
                ]
            ).det()
            signed = -raw
            source = source_rows[n]
            source_signed = flint.arb(source["signed_det"])
            stored = rows[n]
            checks = {
                "rebuild_positive": bool(signed > 0),
                "source_positive": bool(source_signed > 0),
                "source_overlap": bool(signed.overlaps(source_signed)),
                "stored_Q11_contains_rebuild": bool(
                    flint.arb(stored["Q11_ball"]).contains(signed)
                ),
                "stored_H11_contains_rebuild": bool(
                    flint.arb(stored["raw_H11_ball"]).contains(raw)
                ),
                "source_sigma": int(source.get("sigma", 0)) == -1,
                "source_certified": (
                    source.get("classification") == "positive"
                    and source.get("ok") is True
                    and source.get("contains_zero") is False
                    and source.get("sign_separated") is True
                ),
            }
            for name, passed in checks.items():
                if not passed:
                    issues.append(f"{name} failed at n={n}")
            expected = {
                "coefficient_range": [n, n + 20],
                "probe_m": 10,
                "matrix_dimension": 11,
                "probe_sigma": -1,
                "epsilon_11": -1,
                "source_Q11_ball": source["signed_det"],
                "source_overlap": True,
                "classification": "positive",
            }
            for key, value in expected.items():
                if stored.get(key) != value:
                    issues.append(f"stored {key} changed at n={n}")

    finite = artifact.get("finite", {})
    if finite.get("theorem") != (
        "Q_(11,n)(0)>0 for every integer 0<=n<=3"
    ):
        issues.append("finite theorem changed")
    if finite.get("coefficient_range") != [0, 23] or finite.get("precision_dps") != 520:
        issues.append("finite computation contract changed")
    expected_summary = {
        "coefficient_rows": 24,
        "prefix_rows": 4,
        "positive_Q11_rows": 4,
        "source_overlaps": 4,
        "inconclusive_rows": 0,
    }
    if artifact.get("summary") != expected_summary:
        issues.append("summary changed")
    note = generator.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "sigma=(-1)^(10*11/2)=-1=epsilon_11",
        "Q_(11,n)(0)>0 for every integer 0<=n<=3",
        "does not itself prove that ray",
        "PF-infinity, RH, or `Lambda<=0`",
    ):
        if marker not in note:
            issues.append(f"note misses marker: {marker}")

    if issues:
        for issue in issues:
            print(f"issue: {issue}")
        return 1
    print(
        "validated order-eleven lambda-zero prefix: "
        "24 coefficients, 4 positive Q11 rows, 4 source overlaps, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
