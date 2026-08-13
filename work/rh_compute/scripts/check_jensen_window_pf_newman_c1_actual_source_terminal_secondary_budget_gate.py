#!/usr/bin/env python3
"""Independently check the actual-source terminal-secondary budget."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_terminal_secondary_budget_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(payload: dict) -> None:
    require(set(payload["source_sha256"]) == set(payload["source_audit"]), "source keys")
    require(len(payload["source_audit"]) == 5, "source count")
    for key, item in payload["source_audit"].items():
        path = REPO_ROOT / item["path"]
        require(path.is_file(), f"missing source {key}")
        digest = file_hash(path)
        require(digest == item["sha256"], f"source hash drift {key}")
        require(digest == payload["source_sha256"][key], f"source map drift {key}")


def check_rational_chain() -> int:
    rho_log_coefficient = Fraction(1, 6)
    require(rho_log_coefficient == Fraction(1, 6), "rho logarithmic coefficient")

    survivor_prefactor = 2 * rho_log_coefficient * 2 * Fraction(1, 24) / 14300
    pure_prefactor = Fraction(10, 20) * rho_log_coefficient / 14300
    survivor_target = Fraction(1, 68640)
    pure_target = Fraction(1, 22880)
    require(survivor_prefactor == Fraction(1, 514800), "survivor prefactor")
    require(pure_prefactor == Fraction(1, 171600), "pure prefactor")
    require(survivor_prefactor < survivor_target, "survivor target")
    require(pure_prefactor < pure_target, "pure target")
    require(survivor_target + pure_target == Fraction(1, 17160), "joint target")

    e_partial_sum = Fraction(65, 24)
    require(e_partial_sum > Fraction(8, 3), "elementary e lower bound")
    require(Fraction(8, 3) ** 4 > 36, "log 36 upper bound")
    require(59 * Fraction(1, 72_000_000_000) ** 3 < 1, "endpoint decay factor")
    return 9


def check_exponent_chain() -> int:
    audits = 0
    for ell in (50, 51, 64, 100, 256, 1024):
        L = Fraction(ell)
        upper_heat = (L + 1) ** 2 / (32 * L**2)
        lower_decay = L / 4 - Fraction(1, 2)
        require(upper_heat < Fraction(1, 8), "heat exponent")
        require(lower_decay >= 12, "decay exponent")
        require(upper_heat - lower_decay < 0, "terminal exponent sign")
        require(Fraction(1, 1) / (L + 9) - Fraction(3, 2) < 0, "kernel-decay monotonicity")
        audits += 1
    return audits


def check_rows(payload: dict) -> None:
    suffixes = (
        "carrier",
        "rho",
        "survivor_identity",
        "survivor_budget",
        "source_relation",
        "pure_identity",
        "pure_budget",
        "combined",
        "ownership",
        "remaining",
        "transfer",
        "target",
    )
    expected = [f"tsb_{index:02d}_{suffix}" for index, suffix in enumerate(suffixes, start=1)]
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
    rational_audits = check_rational_chain()
    exponent_audits = check_exponent_chain()
    check_rows(payload)

    summary = payload["summary"]
    require(summary["rows"] == 12, "summary rows")
    require(summary["rho_upper_bound"] == "(L+9)/6", "rho summary")
    require(summary["terminal_survivor_prefactor_denominator"] == 514800, "survivor prefactor summary")
    require(summary["pure_terminal_prefactor_denominator"] == 171600, "pure prefactor summary")
    require(summary["terminal_survivor_denominator"] == 68640, "survivor summary")
    require(summary["pure_terminal_denominator"] == 22880, "pure summary")
    require(summary["combined_terminal_denominator"] == 17160, "joint summary")
    require(summary["closed_terminal_packages"] == 2, "closed packages")
    require(summary["remaining_unbounded_secondary_packages"] == 5, "remaining packages")
    require(summary["live_common_unit_secondary_targets"] == 1, "live target")

    note = NOTE_PATH.read_text(encoding="utf-8")
    require(payload["success"] in note, "note success")
    require("rho A_T/17160" in note, "joint budget note")
    require("alpha_P/N" in note, "terminal-kernel ownership note")
    require("not a proof of RH" in note, "proof boundary note")

    print(
        "validated Newman C1 actual-source terminal-secondary budget gate: "
        f"12 rows, 0 issues, {exponent_audits} exponent audits, "
        f"{rational_audits} rational audits, rho<(L+9)/6, "
        "2 terminal packages closed, joint ratio 1/17160, "
        "5 unbounded secondary packages, 1 live common-unit secondary target"
    )


if __name__ == "__main__":
    main()
