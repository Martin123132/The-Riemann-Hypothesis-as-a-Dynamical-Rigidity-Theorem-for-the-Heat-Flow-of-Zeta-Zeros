#!/usr/bin/env python3
"""Validate the legacy prime-curvature symmetry audit."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np

import legacy_prime_curvature_symmetry_audit as builder


EXPECTED_IDS = [
    f"lpc_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "provenance",
            "definition",
            "radial_shell",
            "transpose",
            "modular_residue",
            "no_exact_triangle",
            "boundary_axis",
            "censor_interpolate",
            "interpolation_symmetry",
            "skeleton",
            "flow",
            "zero_component",
            "zero_diagonal",
            "null_control",
            "spectral_sign",
            "laplacian",
            "pell_survivor",
            "prime_heat_bridge",
            "current_wall",
        ),
        start=1,
    )
]


def independent_exact_audit(issues: list[str]) -> None:
    primes = builder.primes_below(2_000)
    sampled = primes[primes >= 3][::7]
    if sampled.size != 44 or sampled[0] != 3 or sampled[-1] != 1999:
        issues.append("source prime sample drifted")

    data = builder.defect_data(sampled)
    defect = data["defect"]
    if defect.shape != (44, 44):
        issues.append("defect shape drifted")
    if np.max(np.abs(defect - defect.T)) != 0:
        issues.append("exact transpose symmetry failed")
    if int(np.count_nonzero(defect <= 500)) != 488:
        issues.append("threshold count drifted")
    if int(np.count_nonzero(defect == 0)) != 0:
        issues.append("unexpected exact odd-prime Pythagorean point")
    if set(int(value) for value in np.unique(defect % 8)) != {1, 2, 6}:
        issues.append("mod-8 defect classes drifted")

    p3 = int(np.flatnonzero(sampled == 3)[0])
    p3_nine = (
        np.count_nonzero(defect[p3, :] == 9)
        + np.count_nonzero(defect[:, p3] == 9)
        - int(defect[p3, p3] == 9)
    )
    if int(p3_nine) != 86:
        issues.append("p=3 boundary-axis count drifted")

    pell = []
    for prime in primes[(primes >= 3) & (primes < 2_000)]:
        root = math.isqrt(2 * int(prime) ** 2)
        if 2 * int(prime) ** 2 - root**2 == 1:
            pell.append((int(prime), root))
    if pell != [(5, 7), (29, 41)]:
        issues.append("negative-Pell prime hits drifted")


def independent_mellin_audit(issues: list[str]) -> None:
    mp.mp.dps = 40
    s = mp.mpf("1.3")
    for n in (2, 3, 11):
        numerical = mp.quad(
            lambda y: mp.e ** (-mp.pi * n * n * y) * y ** (s - 1),
            [0, mp.inf],
        )
        exact = mp.gamma(s) * mp.pi ** (-s) * n ** (-2 * s)
        if abs(numerical - exact) > mp.mpf("1e-32"):
            issues.append(f"Gaussian Mellin identity failed at n={n}")


def validate(path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []

    if artifact.get("kind") != builder.STEM:
        issues.append("artifact kind mismatch")
    if artifact.get("status") != (
        "exploratory_symmetry_and_artifact_audit_complete"
    ):
        issues.append("artifact status drifted")
    if artifact.get("sources") != builder.source_audit():
        issues.append("source hash or marker chain drifted")
    if artifact.get("attachments") != builder.attachment_audit():
        issues.append("attachment hash chain drifted")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    if len(rows) != 19:
        issues.append("row count drifted")
    if artifact.get("summary", {}).get("guard_rows") != 5:
        issues.append("guard-row count drifted")
    if artifact.get("summary", {}).get("open_hypothesis_rows") != 1:
        issues.append("open-hypothesis count drifted")

    arithmetic = artifact.get("arithmetic", {})
    expected_arithmetic = {
        "sampled_prime_count": 44,
        "pair_count": 1936,
        "selected_pair_count": 488,
        "pythagorean_equality_count": 0,
        "transpose_max_abs": 0,
        "allowed_defect_residues_mod_8": [1, 2, 6],
        "p_equals_3_axis_defect_9_count": 86,
        "pell_prime_hits_below_2000": [
            {"p": 5, "r": 7},
            {"p": 29, "r": 41},
        ],
        "sampled_pell_hits": [],
    }
    for key, expected in expected_arithmetic.items():
        if arithmetic.get(key) != expected:
            issues.append(f"arithmetic field drifted: {key}")

    interpolation = artifact.get("interpolation", {})
    cubic = interpolation.get("cubic_transpose", {})
    nearest = interpolation.get("nearest_transpose", {})
    if not (8.0 < cubic.get("rms", 0.0) < 11.0):
        issues.append("cubic transpose RMS outside reproduced range")
    if not (200.0 < cubic.get("max_abs", 0.0) < 300.0):
        issues.append("cubic transpose maximum outside reproduced range")
    if nearest.get("max_abs") != 0.0:
        issues.append("nearest interpolation lost transpose symmetry")
    if interpolation.get("flow", {}).get("curl_max_abs", 1.0) > 1e-12:
        issues.append("derived gradient curl is unexpectedly nonzero")
    if "not a mathematical skeletonization" not in (
        interpolation.get("skeleton", {}).get("operation", "")
    ):
        issues.append("skeleton artifact guard missing")

    spectral = artifact.get("spectral", {})
    raw = spectral.get("raw_symmetric_kernels", {})
    expected_signatures = {
        "defect": (22, 0, 22),
        "threshold_adjacency": (20, 0, 24),
        "closeness": (21, 0, 23),
        "exp_minus_defect_over_100": (20, 0, 24),
    }
    for key, expected in expected_signatures.items():
        record = raw.get(key, {})
        actual = (
            record.get("negative_count"),
            record.get("zero_count"),
            record.get("positive_count"),
        )
        if actual != expected:
            issues.append(f"spectral signature drifted: {key}")
    graph = spectral.get("graph_laplacian_of_closeness", {})
    if (
        graph.get("negative_count"),
        graph.get("zero_count"),
        graph.get("positive_count"),
    ) != (0, 1, 43):
        issues.append("graph-Laplacian PSD signature drifted")

    zero = artifact.get("zero_overlay", {})
    if not zero.get("diagonal_is_imposed"):
        issues.append("imposed zero diagonal guard missing")
    if not (14.0 < zero.get("correct_root_min", 0.0) < 15.0):
        issues.append("first zero height drifted")
    if not (77.0 < zero.get("correct_root_max", 0.0) < 78.0):
        issues.append("twentieth zero height drifted")
    for digits in ("30", "40", "50"):
        span = (
            zero.get("precision_records", {})
            .get(digits, {})
            .get("imaginary_residual_span", 1.0)
        )
        if not (0.0 < span < 1e-25):
            issues.append(f"root residual span guard failed at {digits} dps")
    if (
        zero.get("precision_instability", {}).get(
            "40_vs_50_max_abs", 0.0
        )
        < 0.5
    ):
        issues.append("precision-instability certificate weakened")
    if "precision-dependent solver residuals" not in zero.get(
        "verdict", ""
    ):
        issues.append("zero-overlay verdict drifted")

    null = artifact.get("null_controls", {})
    if (
        null.get("count_matched_odd_control", {}).get(
            "transpose_max_abs"
        )
        != 0
    ):
        issues.append("odd-control transpose symmetry failed")
    prime_null = null.get("prime_block_null", {})
    if prime_null.get("repetitions") != 2_000:
        issues.append("prime-block null repetition count drifted")
    if prime_null.get("thresholds") != [10, 25, 50, 100, 200, 300, 500]:
        issues.append("prime-block threshold grid drifted")

    figure = artifact.get("control_figure", {})
    figure_path = builder.REPO_ROOT / figure.get("path", "")
    if not figure_path.is_file():
        issues.append("control figure missing")
    elif builder.file_hash(figure_path) != figure.get("sha256"):
        issues.append("control figure hash drifted")

    note = builder.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Exact Geometry",
        "Rendering Controls",
        "Arithmetic Null",
        "Spectral Sign",
        "Zero Overlay",
        "Old Invariance Guard",
        "Zeta-Aware Replacement",
        "What It Says About The Wall",
        "Route Decision",
        "Boundary",
    ):
        if marker not in note:
            issues.append(f"note marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove",
        "zero/curvature relation",
        "Xi contact gap",
        "signed Mobius",
        "Lambda<=0",
        "RH",
        "prize-level",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    independent_exact_audit(issues)
    independent_mellin_audit(issues)
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=builder.DEFAULT_OUT)
    args = parser.parse_args()
    issues = validate(args.artifact)
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated legacy prime-curvature symmetry audit: "
        "19 rows, 44 sampled primes, 488 retained pairs, "
        "4 indefinite symmetric kernels, 1 imposed-zero-overlay bug, "
        "1 exact radial Mellin bridge"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
