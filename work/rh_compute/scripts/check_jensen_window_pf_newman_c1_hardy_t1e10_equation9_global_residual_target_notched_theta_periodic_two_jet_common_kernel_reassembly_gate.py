#!/usr/bin/env python3
"""Independently check the two-jet common-kernel reassembly gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_periodic_two_jet_common_kernel_reassembly_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True and artifact.get("kind") == STEM, "artifact identity drift")
    decisions = artifact["decision"]
    for key in (
        "per_mode_two_jet_identity_proved",
        "target_harmonic_current_cancels_coefficientwise",
        "target_Phi_coefficients_restored_exactly",
        "common_kernel_reassembled_as_symmetric_pair_series",
        "zero_regulator_pair_series_absolute_and_uniform",
        "R_after_A_ownership_and_limit_order_preserved",
    ):
        require(decisions.get(key) is True, f"missing decision: {key}")
    for key in (
        "harmonic_current_is_independent_error_term",
        "raw_derivative_majorant_closes_target",
        "quantitative_R_after_A_bound_proved",
        "rh_implication",
    ):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    m, K, phi = sp.symbols("m K phi", nonzero=True)
    p_hat = 1 / (2 * sp.pi**2 * m**2)
    i_m = 2 * sp.pi * sp.I * m * (phi + K * p_hat / 2)
    require(sp.simplify(i_m - 2 * sp.pi * sp.I * m * phi - sp.I * K / (2 * sp.pi * m)) == 0, "independent per-mode identity failed")

    E, H_T, R_all, R_T, P_T = sp.symbols("E H_T R_all R_T P_T")
    source = E - sp.I * K * H_T / (2 * sp.pi) + R_all - R_T
    target = R_T + sp.I * K * H_T / (2 * sp.pi) - P_T
    require(sp.simplify(source + target - (E + R_all - P_T)) == 0, "independent harmonic cancellation failed")

    phi_plus, phi_minus, w = sp.symbols("phi_plus phi_minus w")
    pair = 2 * sp.pi * sp.I * m * w * phi_plus + 2 * sp.pi * sp.I * (-m) * w * phi_minus
    require(sp.expand(pair - 2 * sp.pi * sp.I * m * w * (phi_plus - phi_minus)) == 0, "independent pair reduction failed")

    certificate = artifact["common_kernel_certificate"]
    require("H_x+E_x" in certificate["finite_identity"], "half/zero current missing")
    require("hat Phi_x(m)-hat Phi_x(-m)" in certificate["finite_identity"], "symmetric pair missing")
    require("sum_(m=622)^39894" in certificate["finite_identity"], "target bulk missing")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(BUILDER.is_file() and NOTE.is_file(), "builder or note missing")
    note = NOTE.read_text(encoding="utf-8")
    require("cancel coefficientwise before any norm" in note, "harmonic cancellation missing")
    require("must not be counted after (TR5)" in note, "double-count guard missing")
    require("No useful numerical tail constant" in note, "proof boundary missing")
    print("independently checked two-jet common-kernel reassembly", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
