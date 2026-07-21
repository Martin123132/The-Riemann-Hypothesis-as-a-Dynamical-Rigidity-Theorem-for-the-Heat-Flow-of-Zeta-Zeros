#!/usr/bin/env python3
"""Validate the promoted finite Toeplitz/PF certificate ledger.

This is a reproducibility gate for the proof-programme notes.  It checks the
summary JSON files and problem-row logs for the finite ranges that are
advertised as certified in the outputs/*.md notes.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


LAMBDAS: tuple[str, ...] = ("0", "1e-6", "1e-4", "1e-2", "1e-1")
SUFFIX_CANDIDATES: tuple[str, ...] = ("nonzero", "fast", "struct")


@dataclass(frozen=True)
class CertificateSpec:
    name: str
    matrix_size: int
    max_order: int
    tests: int
    positive: int
    structural_zero: int
    lambdas: tuple[str, ...] = LAMBDAS


@dataclass(frozen=True)
class ExplicitSpec:
    name: str
    summary_name: str
    matrix_size: int
    max_order: int
    tests: int
    positive: int
    structural_zero: int
    should_be_nonnegative: bool = True
    negative: int = 0


GRID_CERTS: tuple[CertificateSpec, ...] = (
    CertificateSpec("lambda grid, N=10, orders<=4", 10, 4, 60_625, 19_690, 40_935),
    CertificateSpec("lambda grid, N=10, orders<=5", 10, 5, 124_129, 39_094, 85_035),
    CertificateSpec("lambda grid, N=12, orders<=4", 12, 4, 297_925, 88_309, 209_616),
    CertificateSpec("lambda grid, N=12, orders<=5", 12, 5, 925_189, 258_193, 666_996),
    CertificateSpec("lambda grid, N=14, orders<=4", 14, 4, 1_142_974, 317_968, 825_006),
    CertificateSpec("lambda grid, N=14, orders<=5", 14, 5, 5_150_978, 1_319_969, 3_831_009),
    CertificateSpec("lambda grid, N=16, orders<=4", 16, 4, 3_640_656, 967_096, 2_673_560),
    CertificateSpec("lambda grid, N=16, orders<=5", 16, 5, 22_720_080, 5_471_960, 17_248_120),
    CertificateSpec("lambda grid, N=18, orders<=4", 18, 4, 10_053_189, 2_578_680, 7_474_509),
    CertificateSpec("lambda grid, N=18, orders<=5", 18, 5, 83_463_813, 19_183_464, 64_280_349),
    CertificateSpec("lambda grid, N=20, orders<=4", 20, 4, 24_810_125, 6_192_025, 18_618_100),
    CertificateSpec("lambda grid, N=20, orders<=5", 20, 5, 265_184_141, 58_773_841, 206_410_300),
    CertificateSpec("lambda grid, N=22, orders<=4", 22, 4, 55_934_670, 13_656_434, 42_278_236),
    CertificateSpec("lambda grid, N=22, orders<=5", 22, 5, 749_414_226, 161_341_895, 588_072_331),
    CertificateSpec("lambda grid, N=24, orders<=4", 24, 4, 117_085_204, 28_075_480, 89_009_724),
    CertificateSpec("lambda grid, N=26, orders<=4", 26, 4, 230_368_801, 54_414_126, 175_954_675),
    CertificateSpec("lambda grid, N=28, orders<=4", 28, 4, 430_101_469, 100_304_533, 329_796_936),
)


EXPLICIT_CERTS: tuple[ExplicitSpec, ...] = (
    ExplicitSpec(
        "lambda=0, N=24, orders<=5",
        "arb_toeplitz_taylor_lam0_N24_o5_r1e-80_nonzero_summary.json",
        24,
        5,
        1_923_675_220,
        404_448_400,
        1_519_226_820,
    ),
    ExplicitSpec(
        "lambda=1e-6, N=24, orders<=5",
        "arb_toeplitz_taylor_lam1e-6_N24_o5_r1e-80_nonzero_summary.json",
        24,
        5,
        1_923_675_220,
        404_448_400,
        1_519_226_820,
    ),
    ExplicitSpec(
        "lambda=1e-4, N=24, orders<=5",
        "arb_toeplitz_taylor_lam1e-4_N24_o5_r1e-80_nonzero_summary.json",
        24,
        5,
        1_923_675_220,
        404_448_400,
        1_519_226_820,
    ),
    ExplicitSpec(
        "lambda=1e-2, N=24, orders<=5",
        "arb_toeplitz_taylor_lam1em2_N24_o5_r1e-80_nonzero_summary.json",
        24,
        5,
        1_923_675_220,
        404_448_400,
        1_519_226_820,
    ),
    ExplicitSpec(
        "lambda=1e-1, N=24, orders<=5",
        "arb_toeplitz_taylor_lam1em1_N24_o5_r1e-80_nonzero_summary.json",
        24,
        5,
        1_923_675_220,
        404_448_400,
        1_519_226_820,
    ),
    ExplicitSpec(
        "lambda=0, N=30, orders<=4",
        "arb_toeplitz_taylor_lam0_N30_o4_r1e-80_nonzero_summary.json",
        30,
        4,
        767_707_750,
        177_089_980,
        590_617_770,
    ),
    ExplicitSpec(
        "lambda=1e-6, N=30, orders<=4",
        "arb_toeplitz_taylor_lam1e-6_N30_o4_r1e-80_nonzero_summary.json",
        30,
        4,
        767_707_750,
        177_089_980,
        590_617_770,
    ),
    ExplicitSpec(
        "lambda=1e-4, N=30, orders<=4",
        "arb_toeplitz_taylor_lam1e-4_N30_o4_r1e-80_nonzero_summary.json",
        30,
        4,
        767_707_750,
        177_089_980,
        590_617_770,
    ),
    ExplicitSpec(
        "lambda=1e-2, N=30, orders<=4",
        "arb_toeplitz_taylor_lam1e-2_N30_o4_r1e-80_nonzero_summary.json",
        30,
        4,
        767_707_750,
        177_089_980,
        590_617_770,
    ),
    ExplicitSpec(
        "lambda=1e-1, N=30, orders<=4",
        "arb_toeplitz_taylor_lam1e-1_N30_o4_r1e-80_nonzero_summary.json",
        30,
        4,
        767_707_750,
        177_089_980,
        590_617_770,
    ),
)


NEGATIVE_CONTROLS: tuple[ExplicitSpec, ...] = (
    ExplicitSpec(
        "beta negative control, lambda=0, N=10, orders<=3",
        "arb_toeplitz_beta_N10_o3_r1e-80_struct_summary.json",
        10,
        3,
        16_525,
        5_127,
        10_695,
        should_be_nonnegative=False,
        negative=703,
    ),
)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def lambda_prefixes(lam: str) -> tuple[str, ...]:
    if lam == "0":
        return ("", "lam0_")
    return (f"lam{lam}_",)


def resolve_grid_summary(results_dir: Path, spec: CertificateSpec, lam: str) -> Path:
    for prefix in lambda_prefixes(lam):
        for suffix in SUFFIX_CANDIDATES:
            if prefix:
                name = f"arb_toeplitz_taylor_{prefix}N{spec.matrix_size}_o{spec.max_order}_r1e-80_{suffix}_summary.json"
            else:
                name = f"arb_toeplitz_taylor_N{spec.matrix_size}_o{spec.max_order}_r1e-80_{suffix}_summary.json"
            path = results_dir / name
            if path.exists():
                return path
    raise FileNotFoundError(f"no summary found for {spec.name}, lambda={lam}")


def problem_path_for(summary_path: Path) -> Path:
    return summary_path.with_name(summary_path.name.replace("_summary.json", "_problem_rows.jsonl"))


def stderr_path_for(summary_path: Path) -> Path:
    return summary_path.with_name(summary_path.name.replace("_summary.json", ".err.log"))


def bad_count(summary: dict) -> int:
    counts = summary["counts"]
    return (
        int(counts.get("negative", 0))
        + int(counts.get("inconclusive_contains_zero", 0))
        + int(counts.get("unknown", 0))
        + int(counts.get("zero", 0))
    )


def validate_positive_summary(summary_path: Path, spec: CertificateSpec | ExplicitSpec) -> None:
    summary = load_json(summary_path)
    counts = summary["counts"]
    checks = {
        "matrix_size": spec.matrix_size,
        "max_order": spec.max_order,
        "tests": spec.tests,
        "positive": spec.positive,
        "structural_zero": spec.structural_zero,
    }
    actuals = {
        "matrix_size": int(summary["matrix_size"]),
        "max_order": int(summary["max_order"]),
        "tests": int(summary["tests"]),
        "positive": int(counts["positive"]),
        "structural_zero": int(counts["structural_zero"]),
    }
    for key, expected in checks.items():
        if actuals[key] != expected:
            raise AssertionError(f"{summary_path.name}: {key} {actuals[key]} != {expected}")

    evaluated = summary.get("evaluated_nonstructural_tests")
    if evaluated is not None and int(evaluated) != spec.positive:
        raise AssertionError(f"{summary_path.name}: evaluated_nonstructural_tests {evaluated} != {spec.positive}")
    if bad_count(summary) != 0:
        raise AssertionError(f"{summary_path.name}: expected zero bad counts, got {bad_count(summary)}")
    if not bool(summary.get("ok_for_nonnegative")):
        raise AssertionError(f"{summary_path.name}: ok_for_nonnegative is false")
    if int(summary.get("recorded_problem_rows", 0)) != 0:
        raise AssertionError(f"{summary_path.name}: recorded_problem_rows is nonzero")

    problem_path = problem_path_for(summary_path)
    if not problem_path.exists():
        raise FileNotFoundError(f"{summary_path.name}: missing problem-row file {problem_path.name}")
    if problem_path.stat().st_size != 0:
        raise AssertionError(f"{summary_path.name}: problem-row file is nonempty")

    stderr_path = stderr_path_for(summary_path)
    if stderr_path.exists() and stderr_path.stat().st_size != 0:
        raise AssertionError(f"{summary_path.name}: stderr log is nonempty")


def validate_negative_control(results_dir: Path, spec: ExplicitSpec) -> None:
    summary_path = results_dir / spec.summary_name
    summary = load_json(summary_path)
    counts = summary["counts"]
    if int(summary["tests"]) != spec.tests:
        raise AssertionError(f"{summary_path.name}: tests mismatch")
    if int(counts["positive"]) != spec.positive:
        raise AssertionError(f"{summary_path.name}: positive mismatch")
    if int(counts["structural_zero"]) != spec.structural_zero:
        raise AssertionError(f"{summary_path.name}: structural_zero mismatch")
    if int(counts["negative"]) != spec.negative:
        raise AssertionError(f"{summary_path.name}: negative count mismatch")
    if bool(summary.get("ok_for_nonnegative")):
        raise AssertionError(f"{summary_path.name}: negative control unexpectedly passed")
    if int(summary.get("recorded_problem_rows", 0)) <= 0:
        raise AssertionError(f"{summary_path.name}: negative control recorded no problem rows")
    problem_path = problem_path_for(summary_path)
    if not problem_path.exists() or problem_path.stat().st_size == 0:
        raise AssertionError(f"{summary_path.name}: negative-control problem-row file missing or empty")


def iter_positive_summaries(results_dir: Path) -> Iterable[tuple[str, Path, CertificateSpec | ExplicitSpec]]:
    for spec in GRID_CERTS:
        for lam in spec.lambdas:
            yield f"{spec.name}, lambda={lam}", resolve_grid_summary(results_dir, spec, lam), spec
    for spec in EXPLICIT_CERTS:
        yield spec.name, results_dir / spec.summary_name, spec


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", default="work/rh_compute/results", type=Path)
    args = parser.parse_args()
    toeplitz_dir = args.results_dir / "arb_toeplitz"

    checked_positive = 0
    for label, summary_path, spec in iter_positive_summaries(toeplitz_dir):
        validate_positive_summary(summary_path, spec)
        checked_positive += 1
        print(f"OK positive: {label} :: {summary_path.name}")

    for spec in NEGATIVE_CONTROLS:
        validate_negative_control(toeplitz_dir, spec)
        print(f"OK negative-control: {spec.name} :: {spec.summary_name}")

    print(
        "validated "
        f"{checked_positive} promoted positive certificate summaries "
        f"and {len(NEGATIVE_CONTROLS)} negative control"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
