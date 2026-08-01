#!/usr/bin/env python3
"""Validate the finite Suzuki k=2 residual scout."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_suzuki_k2_residual_scout.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_suzuki_k2_residual_scout.md"
)
EXPECTED_LABELS = ("1/4", "1/8", "1/16", "1/32")
REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Suzuki K2 Residual Scout",
    "Status: finite double-precision reconnaissance",
    "This is not a proof",
    "q_(omega,1)=-2*xi'(1/2+omega)/xi(1/2+omega)",
    "integrates only through `x=5000`",
    "positive-time Hardy support",
    "Boundary spectral energy is",
    "Finite double-precision time-domain quadrature only",
)


def xi(value):
    return (
        value
        * (value - 1)
        * mp.power(mp.pi, -value / 2)
        * mp.gamma(value / 2)
        * mp.zeta(value)
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = json.loads(args.result.read_text(encoding="utf-8"))
    note = args.note.read_text(encoding="utf-8")
    issues: list[str] = []
    if payload.get("kind") != "jensen_window_pf_suzuki_k2_residual_scout":
        issues.append("bad kind")
    summaries = payload.get("summaries", [])
    completed = bool(payload.get("runtime", {}).get("completed"))
    expected_count = 4 if completed else len(summaries)
    if len(summaries) != expected_count:
        issues.append("inconsistent completed summary count")
    labels = tuple(item.get("omega") for item in summaries)
    if labels != EXPECTED_LABELS[: len(labels)]:
        issues.append(f"bad omega order: {labels!r}")

    mp.mp.dps = 50
    for item in summaries:
        omega = mp.mpf(str(item["omega_float"]))
        point = mp.mpf("0.5") + omega
        expected_q1 = -2 * mp.diff(lambda s: mp.log(xi(s)), point)
        if abs(mp.mpf(str(item["q1"])) - expected_q1) > mp.mpf("1e-12"):
            issues.append(f"{item['omega']}: q1 mismatch")
        for key in (
            "initial_residual",
            "final_residual",
            "fine_finite_energy",
            "coarse_finite_energy",
            "tail_quarter_rms",
        ):
            if not math.isfinite(float(item[key])):
                issues.append(f"{item['omega']}: nonfinite {key}")
        if abs(float(item["initial_residual"]) + float(item["q1"])) > 1e-12:
            issues.append(f"{item['omega']}: bad initial residual")
        if float(item["recurrence_identity_error"]) > 1e-12:
            issues.append(f"{item['omega']}: recurrence mismatch")
        if float(item["fine_finite_energy"]) < 0:
            issues.append(f"{item['omega']}: negative finite energy")

    runtime = payload.get("runtime", {})
    if runtime.get("resource_mode") != "daytime_one_worker":
        issues.append("bad resource mode")
    if runtime.get("thread_limit") != 1:
        issues.append("thread limit is not one")
    audit = payload.get("audit", {})
    for flag in (
        "l2_tail_proved",
        "hardy_causality_proved",
        "rh_proved",
        "lambda_le_zero_proved",
    ):
        if audit.get(flag) is not False:
            issues.append(f"bad proof flag: {flag}")

    for needle in REQUIRED_NOTE_STRINGS:
        if needle not in note:
            issues.append(f"note missing: {needle}")
    lowered = note.lower()
    for forbidden in (
        "this proves rh",
        "therefore rh is true",
        "finite energy proves",
        "the l2 tail is proved",
    ):
        if forbidden in lowered:
            issues.append(f"forbidden promotion language: {forbidden}")

    if issues:
        for issue in issues:
            print(f"ISSUE {issue}")
        print(f"failed Suzuki k=2 residual scout: {len(issues)} issues")
        return 1
    print(
        "validated Suzuki k=2 residual scout: "
        f"{len(summaries)} shifts, "
        f"{payload['grid']['fine_count']} fine points, "
        "0 issues, finite evidence only"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
