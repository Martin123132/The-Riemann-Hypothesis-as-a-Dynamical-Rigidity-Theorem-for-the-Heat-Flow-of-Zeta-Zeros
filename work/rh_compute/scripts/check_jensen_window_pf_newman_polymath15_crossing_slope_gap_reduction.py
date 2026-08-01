#!/usr/bin/env python3
"""Independently check the crossing Gram and slope-gap reduction."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import mpmath as mp
import sympy as sp

import jensen_window_pf_newman_polymath15_crossing_slope_gap_reduction as reduction


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_polymath15_crossing_slope_gap_reduction"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{STEM}.md"
BUILDER_PATH = REPO_ROOT / "work/rh_compute/scripts" / f"{STEM}.py"
COMPONENT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_critical_component_wronskian_gate.json"
)
PHASE_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction.json"
)
SINGLE_CARRIER_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_single_carrier_adiabatic_benchmark.json"
)


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def check_symbolic_identities(issues: list[str]) -> None:
    a = sp.Matrix(sp.symbols("a0:4", real=True))
    p = sp.Matrix(sp.symbols("p0:4", real=True))
    q = sp.Matrix(sp.symbols("q0:4", real=True))
    ell = sp.symbols("ell", real=True, nonzero=True)
    gram = p * p.T + q * q.T / ell**2
    direct = (p.dot(a)) ** 2 + (q.dot(a)) ** 2 / ell**2
    if sp.expand((a.T * gram * a)[0] - direct) != 0:
        issues.append("scaled Gram identity failed")
    if sp.simplify(gram[:3, :3].det()) != 0:
        issues.append("rank-at-most-two Gram check failed")

    m_plus, m_minus, h_plus, h_minus, d_zero = sp.symbols(
        "m_plus m_minus h_plus h_minus d_zero",
        real=True,
    )
    x_value = m_plus - m_minus
    mean_mass = (m_plus + m_minus) / 2
    direct_u = m_plus * h_plus - m_minus * h_minus + d_zero
    reconstructed_u = (
        mean_mass * (h_plus - h_minus)
        + x_value * (h_plus + h_minus) / 2
        + d_zero
    )
    if sp.expand(direct_u - reconstructed_u) != 0:
        issues.append("near-crossing slope identity failed")

    amplitudes = [sp.Integer(3), sp.Integer(3), sp.Integer(7), sp.Integer(7)]
    cosines = [
        sp.Rational(1, 2),
        -sp.Rational(1, 2),
        sp.Rational(1, 2),
        -sp.Rational(1, 2),
    ]
    sines = [
        sp.sqrt(3) / 2,
        sp.sqrt(3) / 2,
        -sp.sqrt(3) / 2,
        -sp.sqrt(3) / 2,
    ]
    speeds = [-4, -3, -2, -1]
    real_value = sum(
        amplitude * cosine
        for amplitude, cosine in zip(amplitudes, cosines, strict=True)
    )
    imag_value = sum(
        amplitude * sine
        for amplitude, sine in zip(amplitudes, sines, strict=True)
    )
    real_first = sum(
        -amplitude * speed * sine
        for amplitude, speed, sine in zip(
            amplitudes, speeds, sines, strict=True
        )
    )
    imag_first = sum(
        amplitude * speed * cosine
        for amplitude, speed, cosine in zip(
            amplitudes, speeds, cosines, strict=True
        )
    )
    real_second = sum(
        -amplitude * speed**2 * cosine
        for amplitude, speed, cosine in zip(
            amplitudes, speeds, cosines, strict=True
        )
    )
    expected = (0, -4 * sp.sqrt(3), 0, -5, -21)
    actual = (
        sp.simplify(real_value),
        sp.simplify(imag_value),
        sp.simplify(real_first),
        sp.simplify(imag_first),
        sp.simplify(real_second),
    )
    if actual != expected:
        issues.append(f"four-carrier jet drifted: {actual!r}")

    slopes = [
        sp.simplify(-speed * sine / cosine)
        for speed, sine, cosine in zip(
            speeds, sines, cosines, strict=True
        )
    ]
    mass = sp.Integer(5)
    positive_mean = sp.simplify(
        (
            amplitudes[0] * cosines[0] * slopes[0]
            + amplitudes[2] * cosines[2] * slopes[2]
        )
        / mass
    )
    negative_mean = sp.simplify(
        (
            -amplitudes[1] * cosines[1] * slopes[1]
            - amplitudes[3] * cosines[3] * slopes[3]
        )
        / mass
    )
    if positive_mean != -sp.sqrt(3) / 5:
        issues.append("four-carrier positive slope mean drifted")
    if negative_mean != -sp.sqrt(3) / 5:
        issues.append("four-carrier negative slope mean drifted")

    target_constant = sp.sqrt(sp.Integer(8_000_000))
    if sp.simplify(target_constant - 2000 * sp.sqrt(2)) != 0:
        issues.append("Xi target square-root constant drifted")


def main() -> int:
    issues: list[str] = []
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    note = NOTE_PATH.read_text(encoding="utf-8")

    if artifact.get("kind") != STEM:
        issues.append("kind drifted")
    if artifact.get("builder_sha256") != file_hash(BUILDER_PATH):
        issues.append("builder hash mismatch")
    expected_sources = {
        "critical_component_wronskian_gate": file_hash(COMPONENT_RESULT),
        "critical_wronskian_phase_reduction": file_hash(PHASE_RESULT),
        "single_carrier_adiabatic_benchmark": file_hash(
            SINGLE_CARRIER_RESULT
        ),
    }
    if artifact.get("source_sha256") != expected_sources:
        issues.append("source hash contract mismatch")

    expected_ids = [
        "crossing_real_component_jet",
        "scaled_gram_expansion",
        "crossing_rank_collapse",
        "near_crossing_slope_identity",
        "exact_crossing_slope_gap",
        "strong_slope_separation_floor",
        "rotating_half_plane_floor",
        "centered_log_moment_frame",
        "nonzero_cosine_four_carrier_obstruction",
        "xi_weighted_slope_gap_target",
    ]
    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != expected_ids:
        issues.append("row ids or order drifted")
    if not rows or rows[-1].get("status") != "open and not proved":
        issues.append("Xi weighted slope-gap target is not explicitly open")

    check_symbolic_identities(issues)

    fine = reduction.diagnostics(reduction.FINE_DPS)
    if artifact.get("diagnostics") != fine:
        issues.append("corrected-crossing diagnostics drifted")
    diagnostic_rows = fine.get("rows", [])
    if len(diagnostic_rows) != 4:
        issues.append("expected four corrected-crossing diagnostics")
    if any(
        mp.mpf(row.get("slope_identity_relative_error", "inf"))
        >= mp.mpf("1e-45")
        for row in diagnostic_rows
    ):
        issues.append("slope identity diagnostic failed")

    stability = artifact.get("precision_stability", {})
    if stability.get("status") != "passed":
        issues.append("precision stability did not pass")
    if mp.mpf(stability.get("maximum_abs_root_delta", "inf")) >= mp.mpf(
        stability.get("root_tolerance", "0")
    ):
        issues.append("root precision stability failed")
    if mp.mpf(
        stability.get("maximum_relative_quantity_delta", "inf")
    ) >= mp.mpf(stability.get("relative_quantity_tolerance", "0")):
        issues.append("quantity precision stability failed")

    expected_summary = {
        "row_count": 10,
        "exact_reductions": 8,
        "exact_countermodels": 1,
        "open_xi_targets": 1,
        "crossing_form_rank_upper_bound": 1,
        "generic_crossing_coercivity": False,
        "all_j_xi_theorem": False,
    }
    if artifact.get("exact_summary") != expected_summary:
        issues.append("exact summary drifted")

    required_note = (
        "rank at most two",
        "Slope-Gap Identity",
        "dimension at least `m-2`",
        "Centered Arithmetic Frame",
        "h_+=h_-=-sqrt(3)/5",
        "2000 sqrt(2) L exp(-3L/4)",
        "fixed `x`",
        "simplicity",
        "not its proof",
        "not logically necessary",
    )
    for marker in required_note:
        if marker not in note:
            issues.append(f"note missing marker {marker!r}")
    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "not proved",
        "does not prove an all-j",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof boundary missing {marker!r}")

    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated Newman crossing slope-gap reduction: "
        "10 rows, 8 exact reductions, 1 exact countermodel, "
        "4 corrected-crossing diagnostics, 1 open Xi target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
