#!/usr/bin/env python3
"""Independently check the actual-source physical-correction budget."""

from __future__ import annotations

import cmath
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_physical_correction_budget_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(payload: dict) -> None:
    require(set(payload["source_sha256"]) == set(payload["source_audit"]), "source keys")
    require(len(payload["source_audit"]) == 9, "source count")
    for key, item in payload["source_audit"].items():
        path = REPO_ROOT / item["path"]
        require(path.is_file(), f"missing source {key}")
        digest = file_hash(path)
        require(digest == item["sha256"], f"source hash drift {key}")
        require(digest == payload["source_sha256"][key], f"source map drift {key}")


def check_symbolic_identities() -> int:
    g, z, z0, w, w0, d, aux = sp.symbols("g z z0 w w0 d aux")
    physical = g * (z**2 + z * d + w) + z * aux
    ideal = z0**2 + w0
    separated = (g - 1) * ideal + g * (
        z**2 - z0**2 + z * d + w - w0
    ) + z * aux
    require(sp.expand(physical - ideal - separated) == 0, "channel subtraction")

    x, u = sp.symbols("x u", real=True)
    z_t = sp.I * (x + u) / 2
    z_h = sp.I * (x - u) / 2
    ux = sp.symbols("u_x", real=True)
    require(
        sp.simplify(z_t**2 - sp.I * ux + (x + u) ** 2 / 4 + sp.I * ux) == 0,
        "transpose ideal",
    )
    require(sp.simplify(z_h**2 + (x - u) ** 2 / 4) == 0, "Hermitian ideal")

    q = sp.symbols("q")
    require(sp.simplify(sum(q**j for j in range(7)) - (1 - q**7) / (1 - q)) == 0, "geometric kernel")
    return 4


def check_rational_constants() -> int:
    r0 = Fraction(53, 2)
    delta_b = (
        Fraction(2001, 1000)
        + Fraction(1003, 1000)
        * (
            Fraction(1, 1500) / r0
            + Fraction(1, 1000) / r0
            + Fraction(1, 500) / r0**2
        )
        + Fraction(3, 1000) / r0
    )
    require(delta_b == Fraction(8_431_968_331, 4_213_500_000), "delta-B coefficient")
    require(delta_b < 3, "delta-B close")

    h0 = Fraction(1, 72_000_000_000)
    nonterminal = 6 * h0 * 120**2 * r0**3 * 304
    require(nonterminal == Fraction(8_485_989, 1_250_000), "nonterminal value")
    require(nonterminal < Fraction(679, 100), "nonterminal threshold")

    terminal = 2 * 59 * h0**3 * 120**2 * r0**2
    require(terminal < Fraction(1, 1000), "terminal threshold")
    require(nonterminal + terminal < 7, "complete functional")
    require(Fraction(7, 14300) < Fraction(1, 2000), "anchor conversion")

    secondary = Fraction(1, 6) + Fraction(1, 17160) + Fraction(1, 2000)
    require(secondary == Fraction(143479, 858000), "secondary fraction")
    require(Fraction(49, 50) - secondary == Fraction(697361, 858000), "remaining reserve")

    derivative = -Fraction(1, 2) + Fraction(2, 120) + Fraction(3, 53) + Fraction(6, 304)
    terminal_derivative = (
        Fraction(1, 59)
        - Fraction(3, 2)
        + Fraction(2, 120)
        + Fraction(2, 53)
    )
    require(derivative < 0 and terminal_derivative < 0, "monotonicity")
    return 9


def check_kernel_samples() -> int:
    audits = 0
    v = 0.137
    for m in (1, 2, 3, 7, 31, 101):
        direct = sum(cmath.exp(2j * math.pi * j * v) for j in range(m))
        quotient = cmath.exp(1j * math.pi * (m - 1) * v) * math.sin(math.pi * m * v) / math.sin(math.pi * v)
        require(abs(direct - quotient) < 2e-12 * max(1, m), f"kernel identity m={m}")
        split_majorant = 2 * (
            m * (1 / (2 * m))
            + 0.5 * math.log(0.5 / (1 / (2 * m)))
        )
        require(abs(split_majorant - (1 + math.log(m))) < 1e-12, f"L1 split m={m}")
        audits += 1
    return audits


def check_rows(payload: dict) -> None:
    suffixes = (
        "fields",
        "factor",
        "bilinear",
        "polynomial",
        "recombine",
        "kernel",
        "l1",
        "scale",
        "carrier",
        "roster",
        "close",
        "terminal",
        "anchor",
        "fraction",
        "reduced",
        "target",
    )
    expected = [f"pcb_{index:02d}_{suffix}" for index, suffix in enumerate(suffixes, start=1)]
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
    rational_audits = check_rational_constants()
    kernel_audits = check_kernel_samples()
    check_rows(payload)

    exact = payload["exact"]
    require(exact["delta_b_coefficient_numerator"] == 8_431_968_331, "stored delta-B numerator")
    require(exact["delta_b_coefficient_denominator"] == 4_213_500_000, "stored delta-B denominator")
    require(exact["nonterminal_numerator"] == 8_485_989, "stored close numerator")
    require(exact["nonterminal_denominator"] == 1_250_000, "stored close denominator")

    summary = payload["summary"]
    require(summary["rows"] == 16, "summary rows")
    require(summary["closed_terminal_coefficient_channels"] == 2, "closed channels")
    require(summary["contiguous_kernel_l1_theorems"] == 1, "kernel theorem")
    require(summary["physical_correction_budget_denominator"] == 2000, "correction budget")
    require(summary["known_secondary_fraction"] == "143479/858000", "known secondary")
    require(summary["remaining_reserve_fraction"] == "697361/858000", "remaining reserve")
    require(summary["remaining_unbounded_secondary_packages"] == 3, "remaining packages")
    require(summary["live_common_unit_three_package_targets"] == 1, "live target")

    note = NOTE_PATH.read_text(encoding="utf-8")
    require(payload["success"] in note, "note success")
    require("rho A_T/2000" in note, "correction budget note")
    require("not a proof of RH" in note, "proof boundary note")

    print(
        "validated Newman C1 actual-source physical-correction budget gate: "
        f"16 rows, 0 issues, {symbolic_audits} symbolic audits, "
        f"{rational_audits} rational audits, {kernel_audits} kernel audits, "
        "2 terminal coefficient channels closed, 1 contiguous-kernel L1 theorem, "
        "|E_Delta|<rho*A_T/2000, 3 unbounded secondary packages, "
        "1 live common-unit three-package target"
    )


if __name__ == "__main__":
    main()
