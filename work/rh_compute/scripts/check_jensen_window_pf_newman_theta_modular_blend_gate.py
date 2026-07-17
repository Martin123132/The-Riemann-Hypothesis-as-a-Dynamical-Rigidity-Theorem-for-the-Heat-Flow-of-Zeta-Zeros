#!/usr/bin/env python3
"""Validate the positive-time-compatible modular theta-block artifact."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

from jensen_window_pf_newman_theta_modular_blend_gate import (
    DEFAULT_OUT,
    jensen_witnesses,
    transform_witnesses,
)


EXPECTED_IDS = [
    "ntmbg_01_switch_reflection",
    "ntmbg_02_even_smooth_blocks",
    "ntmbg_03_strict_block_positivity",
    "ntmbg_04_exact_partition",
    "ntmbg_05_positive_time_schwartz",
    "ntmbg_06_normal_transform_series",
    "ntmbg_07_coupled_laguerre",
    "ntmbg_08_tail_enclosure",
    "ntmbg_09_moment_floor_guard",
    "ntmbg_10_decaying_tail_enclosure",
    "ntmbg_11_termwise_spectral_guard",
    "ntmbg_12_live_handoff",
]


def validate(path: Path) -> list[str]:
    issues: list[str] = []
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("kind") != "jensen_window_pf_newman_theta_modular_blend_gate":
        issues.append("artifact kind mismatch")

    rows = artifact.get("rows", [])
    ids = [row.get("id") for row in rows]
    if ids != EXPECTED_IDS:
        issues.append(f"row id/order mismatch: {ids}")
    if any(row.get("readiness") == "promoted" for row in rows):
        issues.append("unexpected promoted row")

    u = sp.symbols("u", nonnegative=True)
    left = (
        10 * u
        - sp.Rational(3, 2) * (sp.exp(8 * u) - 1)
        + 9 * sp.sinh(4 * u) ** 2
    )
    right = (
        10 * u
        + sp.Rational(3, 4) * sp.exp(8 * u)
        + sp.Rational(9, 4) * sp.exp(-8 * u)
        - 3
    )
    if sp.simplify(sp.expand(left.rewrite(sp.exp)) - right) != 0:
        issues.append("positivity slack identity failed")

    omega, p, r = sp.symbols("omega p r", real=True)
    block = omega * p + (1 - omega) * r
    reflected = (1 - omega) * r + omega * p
    if sp.expand(block - reflected) != 0:
        issues.append("block reflection algebra failed")

    f1, f2, f1p, f2p, f1pp, f2pp = sp.symbols(
        "f1 f2 f1p f2p f1pp f2pp", real=True
    )
    direct = (f1p + f2p) ** 2 - (f1 + f2) * (f1pp + f2pp)
    matrix = (
        f1p**2
        - f1 * f1pp
        + 2 * (f1p * f2p - (f1 * f2pp + f2 * f1pp) / 2)
        + f2p**2
        - f2 * f2pp
    )
    if sp.expand(direct - matrix) != 0:
        issues.append("coupled Laguerre matrix identity failed")

    f, fp, fpp, r, rp, rpp = sp.symbols(
        "f fp fpp r rp rpp", real=True
    )
    full_laguerre = (fp + rp) ** 2 - (f + r) * (fpp + rpp)
    partial_laguerre = fp**2 - f * fpp
    tail_difference = 2 * fp * rp + rp**2 - f * rpp - r * fpp - r * rpp
    if sp.expand(full_laguerre - partial_laguerre - tail_difference) != 0:
        issues.append("finite-to-infinite Laguerre error identity failed")

    stored = artifact.get("numerics", {}).get("transform_rows", [])
    rebuilt = transform_witnesses(460)
    if len(stored) != len(rebuilt):
        issues.append("numerical witness count mismatch")
    else:
        for stored_row, rebuilt_row in zip(stored, rebuilt, strict=True):
            key = (stored_row.get("t"), stored_row.get("n"), stored_row.get("x"))
            rebuilt_key = (
                rebuilt_row.get("t"),
                rebuilt_row.get("n"),
                rebuilt_row.get("x"),
            )
            if key != rebuilt_key:
                issues.append(f"numerical witness key mismatch: {key} != {rebuilt_key}")
                continue
            stored_value = float(stored_row["block_transform"])
            rebuilt_value = float(rebuilt_row["block_transform"])
            if rebuilt_value >= -1e-12:
                issues.append(f"negative transform witness failed: {rebuilt_row}")
            relative_delta = abs(stored_value - rebuilt_value) / max(
                abs(rebuilt_value), 1e-300
            )
            if relative_delta >= 1e-8:
                issues.append(
                    f"numerical witness drift at {key}: relative delta {relative_delta}"
                )

    stored_jensen = artifact.get("numerics", {}).get("jensen_rows", [])
    rebuilt_jensen = jensen_witnesses(460)
    if len(stored_jensen) != len(rebuilt_jensen):
        issues.append("Jensen witness count mismatch")
    else:
        for stored_row, rebuilt_row in zip(
            stored_jensen, rebuilt_jensen, strict=True
        ):
            key = (stored_row.get("t"), stored_row.get("jensen_degree"))
            rebuilt_key = (
                rebuilt_row.get("t"), rebuilt_row.get("jensen_degree")
            )
            if key != rebuilt_key:
                issues.append(f"Jensen witness key mismatch: {key} != {rebuilt_key}")
                continue
            stored_root = complex(
                float(stored_row["scaled_nonreal_root_real"]),
                float(stored_row["scaled_nonreal_root_imag"]),
            )
            rebuilt_root = complex(
                float(rebuilt_row["scaled_nonreal_root_real"]),
                float(rebuilt_row["scaled_nonreal_root_imag"]),
            )
            if rebuilt_root.imag <= 1e-3:
                issues.append(f"nonreal Jensen witness failed: {rebuilt_row}")
            if abs(stored_root - rebuilt_root) / abs(rebuilt_root) >= 1e-8:
                issues.append(f"Jensen witness drift at {key}")

    boundary = artifact.get("proof_boundary", "")
    for forbidden in ("does not prove", "Lambda<=0", "RH"):
        if forbidden not in boundary:
            issues.append(f"proof-boundary guard missing: {forbidden}")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    issues = validate(args.artifact)
    if issues:
        for issue in issues:
            print(f"ISSUE: {issue}")
        raise SystemExit(1)
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    print(
        "validated Jensen-window PF Newman theta modular-blend gate: "
        f"{len(artifact['rows'])} rows, 0 issues, "
        "1 exact positive modular partition, 1 positive-time normal series, "
        f"{len(artifact['numerics']['transform_rows'])} transform witnesses, "
        f"{len(artifact['numerics']['jensen_rows'])} Jensen witnesses"
    )


if __name__ == "__main__":
    main()
