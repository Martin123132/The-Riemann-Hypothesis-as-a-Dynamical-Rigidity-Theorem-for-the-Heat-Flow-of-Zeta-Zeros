#!/usr/bin/env python3
"""Replay the inner/outer corrected-defect split on an alternate grid."""

from __future__ import annotations

import hashlib
import importlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = REPO_ROOT / "work/rh_compute/scripts"
sys.path.insert(0, str(SCRIPT_ROOT))
gate = importlib.import_module(
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_inner_outer_corrected_defect_projection_gate"
)

RESULT = REPO_ROOT / f"work/rh_compute/results/{gate.STEM}.json"
CHECK_PANELS = 32_768


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    gate.ctx.dps = gate.PRECISION
    require(RESULT.is_file(), "missing result artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "saved gate did not pass")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file() and file_hash(path) == record["sha256"], f"dependency hash drift: {path}")

    c = gate.arb(gate.C)
    tstar = gate.arb.pi() * c * c / 8
    beta = tstar ** (gate.arb(1) / 3)
    mu_radius = gate.arb.pi() / (16 * tstar)
    y_max = gate.arb.pi() * c / (2 * beta)
    scale = (gate.PANELS // CHECK_PANELS) ** 2
    saved_inner_error = gate.arb(artifact["certificate"]["inner"]["midpoint_error_ball"]).upper()
    saved_outer_error = gate.arb(artifact["certificate"]["outer"]["midpoint_error_ball"]).upper()

    # Three nodes, a different panel count, and independently reconstructed
    # DFT arrays test both branch assignments and midpoint convergence.
    for index in (0, 2, 4):
        mu = mu_radius * index / 4
        inner, outer = gate.grouped_components(mu, c, beta, y_max, panels=CHECK_PANELS)
        saved = artifact["node_rows"][index]
        saved_inner = gate.acb(
            gate.arb(saved["inner_beta4_value_ball"]["real_ball"]),
            gate.arb(saved["inner_beta4_value_ball"]["imag_ball"]),
        )
        saved_outer = gate.acb(
            gate.arb(saved["outer_beta4_value_ball"]["real_ball"]),
            gate.arb(saved["outer_beta4_value_ball"]["imag_ball"]),
        )
        require(abs(inner - saved_inner).upper() <= saved_inner_error * (scale + 1), f"inner alternate-grid mismatch at node {index}")
        require(abs(outer - saved_outer).upper() <= saved_outer_error * (scale + 1), f"outer alternate-grid mismatch at node {index}")

    decision = artifact["decision"]
    require(decision["inner_outer_beta4_grouped_integrals_certified"] is True, "component integral decision drift")
    require(decision["exact_finite_t_profile_correction_split_by_whole_kernel_norm"] is True, "profile split drift")
    require(decision["inner_plus_outer_reconstructs_prior_corrected_statistic"] is True, "sum reconstruction drift")
    require(decision["inner_corrected_projection_uniformly_negative"] is True, "inner negative sign drift")
    target = gate.arb(artifact["certificate"]["sufficient_paired_companion_real_upper_target"])
    require(target.lower() > gate.arb("6.84e-5"), "paired companion target drift")
    require(decision["total_corrected_sign_transferred_without_split"] is False, "total sign was transferred without split")
    require(decision["paired_companion_projection_bounded"] is False, "paired companion was overpromoted")
    require(decision["signed_fold_increment_proved"] is False, "fold increment was overpromoted")
    require(decision["complete_T_upper_proved"] is False, "T_upper was overpromoted")
    print("PASS: alternate-grid inner/outer corrected-defect projection replay; paired companion remains open")


if __name__ == "__main__":
    main()
