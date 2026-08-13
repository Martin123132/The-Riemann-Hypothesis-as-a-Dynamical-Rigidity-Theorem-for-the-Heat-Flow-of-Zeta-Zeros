#!/usr/bin/env python3
"""Replay the C1 carrier-near odd-harmonic anchor gate."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_c1_carrier_near_odd_harmonic_anchor_gate as target  # noqa: E402


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def independent_symbolic_replay() -> None:
    d_odd, b, d, u_x, u_n, e_n, s_1 = sp.symbols(
        "D_N b d u_x u_N E_N S_1", real=True
    )
    p_h_1 = sp.I * b * d**2 / 4
    p_t_1 = sp.I * d * b**2 / 4 - u_x * d
    p_t_n = -2 * u_x * u_n
    kappa = 1 / (2 * sp.pi * sp.I)
    terminal = 2 * sp.I * u_n * e_n * kappa * (-sp.I * u_x) * s_1

    hermitian = sp.simplify(sp.I * d_odd * p_h_1 / sp.pi)
    transpose = sp.expand(
        sp.I * d_odd * (p_t_1 - e_n * p_t_n) / sp.pi + terminal
    )
    require(
        sp.simplify(hermitian + d_odd * b * d**2 / (4 * sp.pi)) == 0,
        "independent Hermitian replay failed",
    )
    expected_t = (
        -d_odd * d * b**2 / (4 * sp.pi)
        - sp.I * d_odd * u_x * d / sp.pi
        + sp.I * e_n * u_x * u_n * (2 * d_odd - s_1) / sp.pi
    )
    require(sp.simplify(transpose - expected_t) == 0, "independent transpose replay failed")


def main() -> int:
    require(target.RESULT.is_file(), f"missing result: {target.RESULT}")
    require(target.NOTE.is_file(), f"missing note: {target.NOTE}")

    stored = json.loads(target.RESULT.read_text(encoding="utf-8"))
    rebuilt = target.build_payload()
    require(stored == rebuilt, "stored payload differs from deterministic rebuild")
    require(target.NOTE.read_text(encoding="utf-8") == target.render_note(rebuilt), "note drifted")

    for audit in stored["source_audit"].values():
        path = target.REPO_ROOT / audit["path"]
        require(path.is_file(), f"missing pinned source: {path}")
        require(target.sha256_path(path) == audit["sha256"], f"source hash drifted: {path}")

    summary = stored["summary"]
    require(summary["rows"] == 18, "row count drifted")
    require(summary["half_cell_cases"] == 48, "half-cell case count drifted")
    require(summary["half_cell_identities"] == 96, "half-cell identity count drifted")
    require(summary["physical_roster_stress_cases"] == 30, "roster stress count drifted")
    require(summary["explicit_raw_anchors"] == 2, "anchor count drifted")
    require(summary["phase_rotation_guards"] == 1, "phase guard drifted")
    require(summary["cutoff_transfer_guards"] == 1, "cutoff guard drifted")
    require(summary["open_stationary_residual_targets"] == 1, "open target count drifted")

    exact = stored["exact"]
    require(exact["physical_bound"]["conclusion"] == "D_N>(1/2)log(a)>L/4", "D_N bound drifted")
    require(exact["hermitian_split"]["identity"] == "mathcal J_H^0=-A_H+R_H", "Hermitian split drifted")
    require(exact["transpose_split"]["identity"] == "mathcal J_T^0=-A_T+R_T+Q_T", "transpose split drifted")

    independent_symbolic_replay()
    print(target.success_line(stored))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
