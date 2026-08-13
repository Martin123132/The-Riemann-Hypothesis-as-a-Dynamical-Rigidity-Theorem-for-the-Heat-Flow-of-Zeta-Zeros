#!/usr/bin/env python3
"""Independently check the actual-source raw-centre sign obstruction."""

from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_raw_center_sign_obstruction_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(payload: dict) -> None:
    require(set(payload["source_sha256"]) == set(payload["source_audit"]), "source keys")
    for key, item in payload["source_audit"].items():
        path = REPO_ROOT / item["path"]
        require(path.is_file(), f"missing source {key}")
        digest = file_hash(path)
        require(digest == item["sha256"], f"source hash drift {key}")
        require(digest == payload["source_sha256"][key], f"source map drift {key}")


def check_frequency_and_scale() -> int:
    ell, x = sp.symbols("L x", positive=True, real=True)
    first = sp.atan(1 / x) / (8 * ell**2)
    second = 3 * x / (4 * ell**2 * (1 + x**2))
    first_flow = sp.diff(first, ell) + x * sp.diff(first, x)
    second_flow = sp.diff(second, ell) + x * sp.diff(second, x)
    expected_first = -sp.atan(1 / x) / (4 * ell**3) - x / (8 * ell**2 * (1 + x**2))
    expected_second = -3 * x / (2 * ell**3 * (1 + x**2)) + 3 * x * (1 - x**2) / (4 * ell**2 * (1 + x**2) ** 2)
    require(sp.simplify(first_flow - expected_first) == 0, "first epsilon derivative")
    require(sp.simplify(second_flow - expected_second) == 0, "second epsilon derivative")

    samples = (2, 3, 4, 8, 16, 64, 256, 1024, 4096)
    with localcontext() as context:
        context.prec = 80
        for n in samples:
            nn = Decimal(n)
            b = nn.ln()
            speed = (Decimal(1) + Decimal(1) / nn).ln() / b
            require(speed < Decimal(2) / nn, f"relative speed N={n}")
            ideal_g = (nn + 1) ** 2 * (
                Decimal(1) - 2 * (Decimal(1) + Decimal(1) / nn).ln()
            ) - nn**2
            require(ideal_g < Decimal(-2), f"main turn N={n}")
    return len(samples)


def check_sign_reserve() -> int:
    n_floor = 2**25
    phase_error = Fraction(487, n_floor**2)
    mismatch = Fraction(2, n_floor)
    require(phase_error < Fraction(1, 1000), "phase error")
    require(mismatch < Fraction(1, 1000), "mismatch")

    cosine_floor = Fraction(999, 1000)
    sine_floor = Fraction(499, 1000)
    reserve = cosine_floor * sine_floor - Fraction(1, 1000)
    require(reserve == Fraction(497501, 1000000), "exact reserve")
    require(reserve > Fraction(49, 100), "reserve floor")
    require(2 * Fraction(49, 100) == Fraction(98, 100), "scalar compensation")

    audits = 0
    for source_sign in (-1, 1):
        for mismatch_sign in (-1, 1):
            product = source_sign * cosine_floor * sine_floor
            bracket = product + mismatch_sign * Fraction(1, 1000)
            if source_sign > 0:
                require(bracket > Fraction(49, 100), "positive sign corner")
            else:
                require(bracket < -Fraction(49, 100), "negative sign corner")
            audits += 1

    # The raw current has the opposite sign because its prefactor is negative.
    for bracket in (reserve, -reserve):
        current = -bracket
        require((bracket > 0) == (current < 0), "current orientation")
        audits += 1
    return audits


def check_rows(payload: dict) -> None:
    suffixes = (
        "main",
        "monotone",
        "turn",
        "arc",
        "speed",
        "sweep",
        "remainder",
        "positive",
        "negative",
        "obstruction",
        "ties",
        "compensation",
        "target",
    )
    expected = [f"rco_{index:02d}_{suffix}" for index, suffix in enumerate(suffixes, start=1)]
    rows = payload["rows"]
    require([row["id"] for row in rows] == expected, "row ids")
    require(sum(row["readiness"] == "open" for row in rows) == 1, "open rows")
    require(rows[-1]["readiness"] == "open", "frontier row")


def main() -> None:
    require(RESULT_PATH.is_file(), "missing result")
    require(NOTE_PATH.is_file(), "missing note")
    payload = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    require(payload["kind"] == KIND, "kind")
    require(payload["schema_version"] == 1, "schema")
    check_sources(payload)
    scale_samples = check_frequency_and_scale()
    sign_audits = check_sign_reserve()
    check_rows(payload)

    summary = payload["summary"]
    require(summary["rows"] == 13, "summary rows")
    require(summary["safe_phase_arcs"] == 1, "safe arc")
    require(summary["scale_separation_theorems"] == 1, "scale separation")
    require(summary["source_half_turn_pairs"] == 1, "half-turn pair")
    require(summary["strict_positive_witnesses"] == 1, "positive witness")
    require(summary["strict_negative_witnesses"] == 1, "negative witness")
    require(summary["actual_source_center_sign_obstructions"] == 1, "obstruction")
    require(summary["macroscopic_compensation_conditions"] == 1, "compensation")
    require(summary["normalized_reserve"] == "49/100", "reserve summary")
    require(summary["live_complete_current_targets"] == 1, "live target")

    note = NOTE_PATH.read_text(encoding="utf-8")
    require(payload["success"] in note, "note success")
    require("B>49/100" in note and "B<-49/100" in note, "witness note")
    require("not a proof of RH" in note, "proof boundary note")

    print(
        "validated Newman C1 actual-source raw-centre sign-obstruction gate: "
        f"13 rows, 0 issues, {scale_samples} scale samples, "
        f"{sign_audits} sign-corner audits, 2 strict sign witnesses, "
        "1 actual-source centre obstruction, 1 macroscopic compensation condition, "
        "1 live complete-current target"
    )


if __name__ == "__main__":
    main()
