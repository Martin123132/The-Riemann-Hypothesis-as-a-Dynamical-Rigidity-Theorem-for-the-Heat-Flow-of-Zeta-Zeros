#!/usr/bin/env python3
"""Independently check the equation-(9)/source-hybrid nonidentification guard."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_endpoint_residual_source_hybrid_nonidentification_guard_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    ctx.dps = 135
    ctx.threads = 1
    require(gate.RESULT.is_file(), "missing nonidentification result")
    require(gate.NOTE.is_file(), "missing nonidentification note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == gate.STEM, "kind drift")
    require(artifact["passed"] is True, "gate is not passed")

    for name, record in artifact["dependencies"].items():
        path = REPO_ROOT / record["path"]
        require(path == gate.DEPENDENCIES[name], f"dependency path drift: {name}")
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {name}")
    require(file_hash(gate.BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(Path(__file__).resolve()) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    q_k, q_p, q_s, target, gamma = sp.symbols("q_k q_p q_s target gamma")
    source_error = q_s - target
    endpoint_error = q_k - target
    analytic_error = q_p - q_k
    evaluator_error = q_s - q_p
    require(
        sp.factor(source_error - endpoint_error - analytic_error - evaluator_error) == 0,
        "independent three-object decomposition failed",
    )
    require(
        sp.expand((q_k - gamma) - endpoint_error + (gamma - target)) == 0,
        "independent Gamma bridge failed",
    )

    dependencies = {
        name: json.loads((REPO_ROOT / record["path"]).read_text(encoding="utf-8"))
        for name, record in artifact["dependencies"].items()
    }
    source_row = next(row for row in dependencies["source_aligned_ledger"]["output_rows"] if int(row["output_index"]) == 8)
    classical_row = next(row for row in dependencies["classical_upper_split"]["output_rows"] if int(row["output_index"]) == 8)
    source_ball = arb(source_row["source_aligned_upper_block_sum_ball"])
    classical_ball = arb(classical_row["hybrid_minus_classical_main_ball"])
    stored_ball = arb(artifact["numerical_context"]["source_hybrid_minus_classical_target_ball"])
    require(source_ball.overlaps(classical_ball), "independent telemetry balls do not overlap")
    require(stored_ball.overlaps(source_ball), "stored telemetry ball drift")
    require(abs(stored_ball).lower() / arb("8.6e-6") > 800, "telemetry ratio replay failed")

    gamma_direct = arb(dependencies["Gamma_bulk"]["certificate"]["aggregate_two_real_absolute_bound_ball"])
    gamma_stored = arb(artifact["Gamma_bridge"]["aggregate_two_real_absolute_bound_ball"])
    require(gamma_direct.overlaps(gamma_stored), "Gamma bound drift")
    require(gamma_stored < arb("3e-19"), "Gamma bridge no longer tiny")

    paper_audit = dependencies["equation10_current"]["paper_audit"]
    require(paper_audit["exact_kummer_integral_precedes_asymptotic_grouping"] is True, "exact/asymptotic order drift")
    require(paper_audit["paper_drops_sine_term_before_equation51"] is True, "sine omission drift")
    require(paper_audit["equation62_and_hybrid_are_not_exact"] is True, "hybrid exactness drift")
    require(dependencies["double_Weber"]["decision"]["entire_finite_equation9_source_roster_reconstructed_exactly"] is True, "double-Weber object drift")

    decision = artifact["decision"]
    require(decision["source_hybrid_telemetry_equals_exact_endpoint_residual"] is False, "false identity promoted")
    require(decision["source_telemetry_lower_bound_transfer_permitted"] is False, "unsafe lower-bound transfer promoted")
    require(decision["exact_endpoint_residual_quantitatively_bounded_here"] is False, "endpoint residual overclaim")
    require(decision["false_rejection_shortcut_closed"] is True, "shortcut guard lost")
    require(decision["rh_implication"] is False, "RH overclaim")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in (
        "Q_K = exact finite equation-(9) Kummer roster",
        "cannot be transferred to the exact endpoint residual",
        "Partial evaluator certificates must not be",
        "no RH",
    ):
        require(token in note, f"note boundary token missing: {token}")

    print(
        "checked endpoint/source nonidentification independently: "
        "three-object identity, two telemetry certificates, Gamma bridge, and fail-closed decisions",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
