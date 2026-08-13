#!/usr/bin/env python3
"""Independently check the actual-source phase-lock/centre-symmetry gate."""

from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_saddle_phase_lock_center_symmetry_gate"
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


def check_reverse_turn() -> tuple[int, int]:
    u = sp.symbols("u", positive=True, real=True)
    lower_log = u - u**2 / 2 + u**3 / 3 - u**4 / 4
    upper = sp.cancel(((1 + u) ** 2 * (1 - 2 * lower_log) - 1) / u**2)
    expected = -2 - 2 * u / 3 + u**2 / 6 + u**3 / 3 + u**4 / 2
    require(sp.expand(upper - expected) == 0, "log polynomial")
    require(
        Fraction(1, 12) + Fraction(1, 12) + Fraction(1, 16)
        == Fraction(11, 48),
        "tail split",
    )
    require(Fraction(2, 3) - Fraction(11, 48) == Fraction(7, 16), "turn gap")

    samples = (2, 3, 4, 7, 16, 64, 256, 1024, 4096)
    with localcontext() as context:
        context.prec = 90
        for n in samples:
            nn = Decimal(n)
            logarithm = (Decimal(1) + Decimal(1) / nn).ln()
            g_n = (nn + 1) ** 2 * (Decimal(1) - 2 * logarithm) - nn**2
            bound = -Decimal(2) - Decimal(7) / (Decimal(16) * nn)
            require(g_n < bound, f"sampled ideal decrement N={n}")

    correction = Fraction(1, 1) + Fraction(2, 9) + Fraction(17512, 9)
    require(correction == 1947, "physical correction coefficient")
    require(Fraction(1947, 2**25) < Fraction(21, 8), "physical turn survives")
    require(2**25 > 742, "cell floor")
    return len(samples), 4


def check_centre_symmetry() -> tuple[int, int]:
    a_t, b, u_n, phi, theta = sp.symbols(
        "A_T b u_N phi theta", positive=True, real=True
    )
    a_h = a_t * (1 + 2 * u_n / b)
    original = a_h * sp.sin(phi) + a_t * sp.sin(phi + 2 * theta)
    factored = 2 * a_t * (
        sp.cos(theta) * sp.sin(phi + theta) + u_n * sp.sin(phi) / b
    )
    require(sp.simplify(sp.expand_trig(original - factored)) == 0, "trig factorization")

    d_n, d_sum, pi = sp.symbols("D_N d pi", positive=True, real=True)
    a_h_direct = d_n * b * d_sum**2 / (4 * pi)
    a_t_direct = d_n * d_sum * b**2 / (4 * pi)
    require(
        sp.simplify(
            (a_h_direct - a_t_direct).subs(d_sum, b + 2 * u_n)
            - d_n * b * (b + 2 * u_n) * u_n / (2 * pi)
        )
        == 0,
        "anchor difference",
    )

    trig_samples = 0
    for phi_value in (0, sp.pi / 7, sp.pi / 2, 7 * sp.pi / 6):
        for theta_value in (-sp.pi / 3, 0, sp.pi / 5, sp.pi):
            for u_value, b_value in ((0, 2), (sp.Rational(1, 7), 3)):
                value = sp.simplify(
                    sp.expand_trig((original - factored).subs(
                        {
                            a_t: 5,
                            phi: phi_value,
                            theta: theta_value,
                            u_n: u_value,
                            b: b_value,
                        }
                    ))
                )
                numeric_value = complex(sp.N(value, 80))
                require(
                    abs(numeric_value) < 1e-70,
                    f"sampled centre identity phi={phi_value}, theta={theta_value}, "
                    f"u={u_value}, b={b_value}: {value} -> {numeric_value}",
                )
                trig_samples += 1

    mismatch_samples = 0
    for n in (2, 3, 8, 64, 1024):
        with localcontext() as context:
            context.prec = 60
            nn = Decimal(n)
            u_value = (Decimal(1) + Decimal(1) / nn).ln()
            b_value = nn.ln()
            require(u_value / b_value < Decimal(2) / nn, f"mismatch N={n}")
            mismatch_samples += 1
    return trig_samples, mismatch_samples


def check_rows(payload: dict) -> None:
    suffixes = (
        "normalizer",
        "lift",
        "lock",
        "identity",
        "endpoint",
        "decrement",
        "scale",
        "error",
        "reverse",
        "circle",
        "inventory",
        "ratio",
        "current",
        "symmetry",
        "branch",
        "mismatch",
        "ties",
        "joint",
        "complete",
    )
    expected = [f"asl_{index:02d}_{suffix}" for index, suffix in enumerate(suffixes, start=1)]
    rows = payload["rows"]
    require([row["id"] for row in rows] == expected, "row ids")
    require(sum(row["readiness"] == "open" for row in rows) == 2, "open rows")
    require(all(row["readiness"] == "open" for row in rows[-2:]), "frontier rows")


def main() -> None:
    require(RESULT_PATH.is_file(), "missing result")
    require(NOTE_PATH.is_file(), "missing note")
    payload = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    require(payload["kind"] == KIND, "kind")
    require(payload["schema_version"] == 1, "schema")
    check_sources(payload)
    decrement_samples, rational_audits = check_reverse_turn()
    trig_samples, mismatch_samples = check_centre_symmetry()
    check_rows(payload)

    summary = payload["summary"]
    require(summary["rows"] == 19, "summary rows")
    require(summary["reverse_full_turn_theorems"] == 1, "reverse theorem")
    require(summary["source_inclusive_circle_theorems"] == 1, "circle theorem")
    require(summary["phase_inventory_corollaries"] == 2, "phase inventory")
    require(summary["physical_error_coefficient"] == 1947, "error coefficient")
    require(summary["exact_center_factorizations"] == 2, "factorizations")
    require(summary["anchor_mismatch_bounds"] == 1, "mismatch bound")
    require(summary["live_joint_phase_targets"] == 1, "joint target")
    require(summary["open_rows"] == 2, "open row summary")

    note = NOTE_PATH.read_text(encoding="utf-8")
    require(payload["success"] in note, "note success")
    require("Delta chi_N<-2pi" in note, "reverse theorem note")
    require("cos(theta_eta)sin(phi_N+theta_eta)" in note, "symmetry note")
    require("not a proof of RH" in note, "proof boundary note")

    print(
        "validated Newman C1 actual-source saddle-phase lock/centre-symmetry gate: "
        f"19 rows, 0 issues, {decrement_samples} decrement samples, "
        f"{rational_audits} rational audits, {trig_samples} centre samples, "
        f"{mismatch_samples} mismatch samples, 1 reverse full-turn theorem, "
        "1 live joint-phase target"
    )


if __name__ == "__main__":
    main()
