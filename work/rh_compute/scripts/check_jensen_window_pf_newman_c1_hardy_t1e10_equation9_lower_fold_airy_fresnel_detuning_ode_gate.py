#!/usr/bin/env python3
"""Independently validate the finite Airy-Fresnel detuning ODE gate."""

from __future__ import annotations

import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

from flint import arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_airy_fresnel_detuning_ode_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    priority = gate.set_low_priority()
    require(gate.RESULT.is_file(), "missing detuning-ODE result")
    ctx.dps = gate.PRECISION
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "detuning-ODE artifact does not pass")
    require(artifact["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash mismatch")
    require(artifact["sources"]["checker"]["sha256"] == gate.file_hash(gate.CHECKER), "checker hash mismatch")
    for name, path in (
        ("normal_form_gate", gate.NORMAL_FORM_GATE),
        ("core_gate", gate.CORE_GATE),
        ("lower_fold_gate", gate.LOWER_FOLD_GATE),
    ):
        require(artifact["dependencies"][name]["sha256"] == gate.file_hash(path), f"{name} hash mismatch")

    saved = artifact["certified_coefficients"]
    p = gate.parameters()
    require(arb(saved["beta_ball"]).overlaps(p["beta"]), "beta mismatch")
    require(arb(saved["lambda_ball"]).overlaps(p["lambda"]), "lambda mismatch")
    require(arb(saved["transport_coefficient_min_ball"]) > arb("0.998"), "transport lower margin failed")
    require(arb(saved["transport_coefficient_max_ball"]) < arb("1.002"), "transport upper margin failed")
    require(arb(saved["real_potential_min_ball"]) > arb("0.026"), "potential margin failed")

    central = gate.ode_residual(39894, p)
    saved_central = next(row for row in saved["rigorous_moment_residuals"] if row["mode"] == 39894)
    require(arb(central["detuning_ball"]).overlaps(arb(saved_central["detuning_ball"])), "central detuning mismatch")
    require(arb(central["residual"]["real_ball"]).contains(0), "central real residual excludes zero")
    require(arb(central["residual"]["imag_ball"]).contains(0), "central imaginary residual excludes zero")

    decisions = artifact["decision"]
    require(decisions["exact_finite_detuning_ODE_proved"] is True, "exact ODE decision missing")
    require(decisions["both_finite_endpoint_sources_retained"] is True, "endpoint-source decision missing")
    require(decisions["all_detuning_ODE_enclosure_proved"] is False, "full-interval boundary corrupted")
    require(decisions["rh_implication"] is False, "RH boundary corrupted")
    print(f"validated lower-fold detuning ODE independently: coefficients pass, central residual contains zero; priority={priority}")


if __name__ == "__main__":
    main()
