#!/usr/bin/env python3
"""Independently check the non-A quadratic-theta reference evaluator gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
from mpmath import mp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_reference_evaluator_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
A = 159_577
L = 2_481_422


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def is_even_integer(value: Fraction) -> bool:
    return value.denominator == 1 and value.numerator % 2 == 0


def independent_minimal_period(x: Fraction, s: Fraction, limit: int) -> int:
    if x == 0:
        return 1
    for period in range(1, limit + 1):
        xp = x * period
        if xp.denominator == 1 and is_even_integer(xp * (period + A + 2 * s)):
            return period
    raise RuntimeError("recorded period is not witnessed")


def mp_cis_pi(value: Fraction):
    reduced = value % 2
    return mp.exp(mp.j * mp.pi * mp.mpf(reduced.numerator) / reduced.denominator)


def independent_periodic_oracle(x: Fraction, s: Fraction, period: int):
    mp.dps = 105
    cycles, remainder = divmod(L, period)
    z = [mp_cis_pi(x * (j * j + j * (A + 2 * s))) for j in range(period)]
    w = []
    alpha = []
    for j in range(period):
        alpha_j = A + 2 * j + 2 * s
        alpha.append(mp.mpf(alpha_j.numerator) / alpha_j.denominator)
        w.append(mp_cis_pi(x * alpha_j**2 / 4))
    z0 = mp.fsum(z)
    z1 = mp.fsum(j * z[j] for j in range(period))
    zr0 = mp.fsum(z[:remainder])
    zr1 = mp.fsum(j * z[j] for j in range(remainder))
    s0 = cycles * z0 + zr0
    s1 = period * cycles * (cycles - 1) * z0 / 2 + cycles * z1
    s1 += cycles * period * zr0 + zr1
    w0 = mp.fsum(w)
    wa = mp.fsum(alpha[j] * w[j] for j in range(period))
    wr0 = mp.fsum(w[:remainder])
    wra = mp.fsum(alpha[j] * w[j] for j in range(remainder))
    current = cycles * wa + period * cycles * (cycles - 1) * w0
    current += 2 * cycles * period * wr0 + wra
    return s0, s1, current


def recorded_mp(record: dict[str, str]):
    return mp.mpc(mp.mpf(record["real"]), mp.mpf(record["imag"]))


def main() -> int:
    require(RESULT.is_file() and NOTE.is_file() and BUILDER.is_file(), "reference-evaluator artifact missing")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("kind") == STEM and artifact.get("passed") is True, "artifact identity drift")
    require(artifact["scope"]["source_cell_count"] == L, "roster length drift")
    require(artifact["scope"]["workers"] == 1, "worker-count drift")
    require(artifact["algorithm"]["block_size"] == 256, "reseed block drift")
    require("math.fsum" in artifact["algorithm"]["summation"], "summation method drift")

    decisions = artifact["decision"]
    for key in (
        "short_high_precision_rows_passed",
        "full_roster_exact_periodic_rows_passed",
        "large_phase_trigonometric_calls_eliminated",
        "reference_evaluator_ready_as_fast_evaluator_oracle",
    ):
        require(decisions.get(key) is True, f"missing evaluator decision: {key}")
    for key in (
        "uniform_floating_error_theorem_proved",
        "interval_certified",
        "fast_incomplete_quadratic_Gauss_evaluator_built",
        "physical_x_s_quadrature_completed",
        "non_A_bound_proved",
        "R_after_A_bound_proved",
        "rh_implication",
    ):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    allowance = artifact["summary"]["mass_normalized_allowance"]
    require(allowance == 2.0e-10, "floating allowance drift")
    rows = artifact["short_high_precision_rows"] + artifact["full_roster_periodic_rows"]
    require(len(rows) == artifact["summary"]["row_count"] == 11, "row-count drift")
    for row in rows:
        require(row.get("passed") is True, "a stored reference row is not passed")
        require(max(row["mass_normalized_errors"].values()) <= allowance, "stored row exceeds allowance")

    periodic_rows = artifact["full_roster_periodic_rows"]
    require(len(periodic_rows) == 7, "periodic row-count drift")
    for row in periodic_rows:
        x = Fraction(row["x"])
        s = Fraction(row["s"])
        period = int(row["minimal_period"])
        require(independent_minimal_period(x, s, period) == period, "minimal period drift")
        congruences = row["period_congruences"]
        require(congruences["xP_is_integer"] is True, "integer congruence missing")
        require(congruences["endpoint_phase_is_even_integer"] is True, "even endpoint congruence missing")

    # Recompute the largest-period full-roster oracle independently.  This
    # checks both cycle formulas and the high-precision values stored by the builder.
    witness = max(periodic_rows, key=lambda row: row["minimal_period"])
    x = Fraction(witness["x"])
    s = Fraction(witness["s"])
    expected = independent_periodic_oracle(x, s, int(witness["minimal_period"]))
    stored = witness["periodic_high_precision"]
    for index, label in enumerate(("S0", "S1", "F")):
        serialization_scale = max(mp.mpf(1), abs(expected[index]))
        require(
            abs(expected[index] - recorded_mp(stored[label])) < mp.mpf("1e-72") * serialization_scale,
            f"independent periodic oracle drift: {label}",
        )

    for source in artifact["sources"].values():
        path = REPO_ROOT / source["path"]
        require(file_hash(path) == source["sha256"], f"source hash drift: {path}")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path}")

    note = NOTE.read_text(encoding="utf-8")
    for token in ("(RE1)", "(RE3)", "(RE4)", "(RE5)", "2,481,422", "No uniform floating-error theorem"):
        require(token in note, f"note token missing: {token}")
    require("RH" in note and "prize-level conclusion" in note, "proof boundary missing")
    print("independently checked the quadratic-theta reference evaluator and periodic oracle", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
