#!/usr/bin/env python3
"""Independently check the actual-source complete-compensation ledger."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_complete_compensation_ledger_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(payload: dict) -> None:
    require(set(payload["source_sha256"]) == set(payload["source_audit"]), "source keys")
    require(len(payload["source_audit"]) == 8, "source count")
    for key, item in payload["source_audit"].items():
        path = REPO_ROOT / item["path"]
        require(path.is_file(), f"missing source {key}")
        digest = file_hash(path)
        require(digest == item["sha256"], f"source hash drift {key}")
        require(digest == payload["source_sha256"][key], f"source map drift {key}")


def check_two_tone() -> int:
    phi, theta, mismatch = sp.symbols("phi theta mismatch", real=True)
    raw = sp.cos(theta) * sp.sin(phi + theta) + mismatch * sp.sin(phi)
    tones = sp.sin(phi + 2 * theta) / 2 + (sp.Rational(1, 2) + mismatch) * sp.sin(phi)
    require(sp.simplify(sp.expand_trig(raw - tones)) == 0, "two-tone identity")

    audits = 0
    for p_num in range(-3, 4):
        for t_num in range(-3, 4):
            for e_num in (-2, 0, 3):
                value = {
                    phi: sp.Rational(p_num, 5),
                    theta: sp.Rational(t_num, 7),
                    mismatch: sp.Rational(e_num, 1000),
                }
                require(sp.simplify(raw.subs(value) - tones.subs(value)) == 0, "tone sample")
                audits += 1
    return audits


def check_compensation_arithmetic() -> int:
    reserve = Fraction(49, 100)
    threshold = 2 * reserve
    require(threshold == Fraction(98, 100), "threshold")

    audits = 0
    for delta in (Fraction(0), Fraction(1, 10), Fraction(49, 100), Fraction(97, 100)):
        require(delta < threshold, "delta range")
        finite_target = threshold - delta
        require(finite_target > 0, "positive finite target")
        # If E_F >= -(threshold-delta) and E_sec >= -delta, their sum
        # cannot be below -threshold.
        lower_sum = -finite_target - delta
        require(lower_sum == -threshold, "budget boundary")
        audits += 1

    for alpha, beta in (
        (Fraction(0), Fraction(0)),
        (Fraction(1, 5), Fraction(1, 4)),
        (Fraction(49, 100), Fraction(49, 100)),
    ):
        require(alpha + beta <= threshold, "retirement budget")
        require(-alpha - beta >= -threshold, "retirement lower bound")
        audits += 1
    return audits


def check_nonpromotion() -> int:
    audits = 0
    for epsilon in (Fraction(1, 2), Fraction(1, 10), Fraction(1, 10**6), Fraction(1, 10**18)):
        for target in (Fraction(1), Fraction(49, 100), Fraction(14300)):
            row = target / epsilon
            require(epsilon * row == target, "coefficient nonpromotion")
            audits += 1
    return audits


def check_rows(payload: dict) -> None:
    suffixes = (
        "normalization",
        "two_tone",
        "partition",
        "finite",
        "remote",
        "terminal",
        "correction",
        "affine",
        "pure_terminal",
        "moving",
        "quadratic",
        "ownership",
        "threshold",
        "budget",
        "nonpromotion",
        "mechanisms",
        "target",
    )
    expected = [f"ccl_{index:02d}_{suffix}" for index, suffix in enumerate(suffixes, start=1)]
    rows = payload["rows"]
    require([row["id"] for row in rows] == expected, "row ids")
    require(sum(row["readiness"] == "open" for row in rows) == 1, "open rows")
    require(rows[-1]["readiness"] == "open", "frontier row")


def check_partition(payload: dict) -> None:
    total = payload["exact"]["physical_partition"]["total"]
    labels = ("E_F", "E_R", "E_term", "E_Delta", "E_aff", "E_aa", "E_move", "E_quad")
    require(len(set(labels)) == 8, "unique package labels")
    for label in labels:
        require(label in total, f"missing package {label}")
    require(payload["exact"]["scale_audit"]["common_unit_secondary_bounds"] == 0, "secondary status")


def main() -> None:
    require(RESULT_PATH.is_file(), "missing result")
    require(NOTE_PATH.is_file(), "missing note")
    payload = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    require(payload["kind"] == KIND, "kind")
    require(payload["schema_version"] == 1, "schema")
    check_sources(payload)
    tone_audits = check_two_tone()
    budget_audits = check_compensation_arithmetic()
    nonpromotion_audits = check_nonpromotion()
    check_partition(payload)
    check_rows(payload)

    summary = payload["summary"]
    require(summary["rows"] == 17, "summary rows")
    require(summary["disjoint_scalar_packages"] == 8, "package count")
    require(summary["two_tone_identities"] == 1, "tone count")
    require(summary["sharp_adverse_thresholds"] == 1, "threshold count")
    require(summary["conditional_budget_lemmas"] == 1, "budget lemma")
    require(summary["common_unit_secondary_bounds"] == 0, "secondary bounds")
    require(summary["uniform_bias_mechanisms"] == 1, "bias mechanism")
    require(summary["phase_correlated_mechanisms"] == 1, "correlated mechanism")
    require(summary["live_signed_finite_package_targets"] == 1, "live target")

    note = NOTE_PATH.read_text(encoding="utf-8")
    require(payload["success"] in note, "note success")
    require("98/100" in note, "threshold note")
    require("not a proof of RH" in note, "proof boundary note")

    print(
        "validated Newman C1 actual-source complete-compensation ledger gate: "
        f"17 rows, 0 issues, 8 disjoint scalar packages, {tone_audits} tone audits, "
        f"{budget_audits} budget audits, {nonpromotion_audits} nonpromotion audits, "
        "1 sharp 98/100 adverse threshold, 0 common-unit secondary bounds, "
        "1 live signed finite-package target"
    )


if __name__ == "__main__":
    main()
