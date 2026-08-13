#!/usr/bin/env python3
"""Independently check the actual-source full-phase remote budget."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_full_phase_remote_budget_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(payload: dict) -> None:
    require(set(payload["source_sha256"]) == set(payload["source_audit"]), "source keys")
    require(len(payload["source_audit"]) == 6, "source count")
    for key, item in payload["source_audit"].items():
        path = REPO_ROOT / item["path"]
        require(path.is_file(), f"missing source {key}")
        digest = file_hash(path)
        require(digest == item["sha256"], f"source hash drift {key}")
        require(digest == payload["source_sha256"][key], f"source map drift {key}")


def check_symbolic_identities() -> int:
    x, k = sp.symbols("x k", positive=True, real=True)
    pair = 1 / (x - k) + 1 / (x + k)
    require(sp.simplify(pair - 2 * x / (x**2 - k**2)) == 0, "pair identity")

    lam = sp.symbols("lambda", real=True)
    t, sigma = sp.symbols("t sigma", real=True)
    p = sp.Function("p")(lam)
    s = t * lam**2 / 4 - sigma * lam
    f = sp.exp(s) * p
    g = sp.diff(s, lam)
    transformed = sp.simplify(sp.diff(f, lam, 2) - sp.diff(f, lam))
    expected = sp.exp(s) * (
        sp.diff(p, lam, 2)
        + (2 * g - 1) * sp.diff(p, lam)
        + (g**2 - g + t / 2) * p
    )
    require(sp.simplify(transformed - expected) == 0, "amplitude second derivative")
    return 2


def check_rational_close() -> int:
    endpoint = Fraction(8, 27) * Fraction(6503, 3000) * Fraction(44, 175)
    remainder = Fraction(1, 14300)
    require(endpoint + remainder < Fraction(1, 6), "one-sixth close")
    require(Fraction(1, 6) + Fraction(1, 17160) == Fraction(2861, 17160), "known secondary")

    audits = 0
    for alpha in (Fraction(3), Fraction(10), Fraction(101, 2), Fraction(1000)):
        n_lower = 2 * alpha - 1
        require(n_lower > Fraction(3, 2) * alpha, "roster lower bound")
        require(alpha / n_lower < Fraction(2, 3), "weighted tail bound")
        audits += 1
    return audits + 2


def check_amplitude_constants() -> int:
    second_boundary = Fraction(29, 2)
    interior = Fraction(223, 8)
    require(second_boundary + interior == Fraction(339, 8), "IBP constant")
    require(Fraction(339, 72) < 5, "physical remainder constant")

    extra = Fraction(1, 3456) + Fraction(1, 21600)
    require(extra < Fraction(1, 1000), "terminal endpoint extras")
    require(Fraction(7, 6) + 1 + Fraction(1, 1000) == Fraction(6503, 3000), "endpoint amplitude")
    return 4


def check_rows(payload: dict) -> None:
    suffixes = (
        "phase",
        "denominator",
        "tail",
        "first",
        "amplitude",
        "second",
        "interior",
        "remainder",
        "anchor",
        "close",
        "remote",
        "cutoff",
        "reduced",
        "transfer",
        "target",
    )
    expected = [f"frb_{index:02d}_{suffix}" for index, suffix in enumerate(suffixes, start=1)]
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
    symbolic_audits = check_symbolic_identities()
    rational_audits = check_rational_close()
    constant_audits = check_amplitude_constants()
    check_rows(payload)

    summary = payload["summary"]
    require(summary["rows"] == 15, "summary rows")
    require(summary["paired_denominator_cancellations"] == 1, "pair cancellation")
    require(summary["closed_ideal_remote_channels"] == 2, "remote channels")
    require(summary["remote_budget_denominator"] == 6, "remote budget")
    require(summary["known_secondary_fraction"] == "2861/17160", "known secondary")
    require(summary["remaining_unbounded_secondary_packages"] == 4, "remaining packages")
    require(summary["live_common_unit_four_package_targets"] == 1, "live target")

    note = NOTE_PATH.read_text(encoding="utf-8")
    require(payload["success"] in note, "note success")
    require("rho A_T/6" in note, "remote budget note")
    require("not a proof of RH" in note, "proof boundary note")

    print(
        "validated Newman C1 actual-source full-phase remote budget gate: "
        f"15 rows, 0 issues, {symbolic_audits} symbolic audits, "
        f"{rational_audits} rational audits, {constant_audits} constant audits, "
        "1 paired denominator cancellation, 2 ideal remote channels closed, "
        "|E_R|<rho*A_T/6, 4 unbounded secondary packages, "
        "1 live common-unit four-package target"
    )


if __name__ == "__main__":
    main()
