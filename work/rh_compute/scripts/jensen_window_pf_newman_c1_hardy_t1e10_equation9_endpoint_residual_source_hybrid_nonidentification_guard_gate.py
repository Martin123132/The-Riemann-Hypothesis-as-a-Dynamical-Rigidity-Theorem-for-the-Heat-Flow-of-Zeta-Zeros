#!/usr/bin/env python3
"""Reject identification of the exact equation-(9) residual with source telemetry."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
import sys

sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_endpoint_residual_source_hybrid_nonidentification_guard_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "source_aligned_ledger": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_source_aligned_hybrid_error_ledger_gate.json",
    "classical_upper_split": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_classical_upper_main_remainder_split_gate.json",
    "equation10_current": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation10_global_theta_current_reduction_gate.json",
    "symmetric_endpoint": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_symmetric_endpoint_tail_reassembly_gate.json",
    "paired_residual": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_paired_target_residual_normal_form_gate.json",
    "Gamma_bulk": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_global_morse_gamma_bulk_gate.json",
    "double_Weber": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_double_weber_exact_source_reconstruction_gate.json",
}

PRECISION = 110
REFERENCE_PHYSICAL_TARGET = arb("8.6e-6")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def symbolic_certificate() -> dict[str, str]:
    q_kummer, q_paper, q_source, target, gamma_bulk = sp.symbols(
        "Q_K Q_P Q_S T G"
    )
    r_kummer = q_kummer - target
    d_analytic = q_paper - q_kummer
    d_evaluator = q_source - q_paper
    e_source = q_source - target
    require(
        sp.expand(e_source - r_kummer - d_analytic - d_evaluator) == 0,
        "source-to-Kummer decomposition failed",
    )

    r_gamma = q_kummer - gamma_bulk
    d_gamma = gamma_bulk - target
    require(
        sp.expand(r_gamma - r_kummer + d_gamma) == 0,
        "Gamma-target decomposition failed",
    )

    free_residual = sp.symbols("r", real=True)
    witness_analytic = e_source - free_residual
    require(
        sp.expand(free_residual + witness_analytic - e_source) == 0,
        "nonidentification witness failed",
    )

    return {
        "objects": "Q_K=exact finite equation-(9) Kummer roster; Q_P=published asymptotic hybrid; Q_S=implemented source hybrid; T=exact classical upper main; G=Gamma bulk",
        "source_identity": "E_source=Q_S-T=(Q_K-T)+(Q_P-Q_K)+(Q_S-Q_P)",
        "endpoint_identity": "R_end^Gamma=Q_K-G=(Q_K-T)-(G-T)",
        "nonidentification_witness": "For any proposed value r of Q_K-T, choosing Q_P-Q_K=E_source-r and Q_S-Q_P=0 preserves the observed E_source.",
        "projection": "Q_K and G are compared only after the common equation-(9) normalization and 2 Re[exp(-i*pi/8) .] projection.",
    }


def central_row(rows: list[dict[str, Any]]) -> dict[str, Any]:
    matches = [row for row in rows if int(row["output_index"]) == 8]
    require(len(matches) == 1, "central output row drift")
    return matches[0]


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["numerical_context"]
    g = artifact["Gamma_bridge"]
    return f"""# Exact endpoint residual versus source telemetry: nonidentification guard

Date: 2026-08-13
Status: exact guard certificate; not a proof of the endpoint residual bound

Use three different upper objects:

```text
Q_K = exact finite equation-(9) Kummer roster,
Q_P = published asymptotic/generalized-Gauss hybrid,
Q_S = implemented source hybrid,
T   = exact classical upper main.
```

They obey the exact bookkeeping identity

```text
E_source=Q_S-T
        =(Q_K-T)+(Q_P-Q_K)+(Q_S-Q_P).                 (NI1)
```

The symmetric endpoint calculation and the double-Weber involution reconstruct
`Q_K`, not `Q_P` or `Q_S`.  The paper audit independently records that the
exact Kummer integral precedes the asymptotic grouping, that the sine member
of each equation-(10) pair is dropped before equation (51), and that the
published hybrid is not exact.  Consequently the saved source discrepancy

```text
Q_S-T = {c['source_hybrid_minus_classical_target_ball']}               (NI2)
```

cannot be transferred to the exact endpoint residual.  Its absolute size is
more than `{c['source_telemetry_to_reference_target_ratio_lower']}` times the
reference `8.6e-6` scale, but that fact belongs only to `Q_S-T`.

The Gamma calculation does close a different bridge.  With `G` denoting the
projected exact Gamma bulk,

```text
R_end^Gamma=Q_K-G=(Q_K-T)-(G-T),                      (NI3)
|G-T| <= {g['aggregate_two_real_absolute_bound_ball']}.
```

Thus the Gamma replacement changes the exact Kummer-to-classical residual by
less than `2.186e-19`; it does not identify that residual with source telemetry.
The missing quantities in (NI1) are the analytic Kummer-to-published defect
`Q_P-Q_K` and the complete published-to-implemented bridge `Q_S-Q_P` in the
same global normalization.  Partial evaluator certificates must not be
silently substituted for either complete quantity.

Pi provenance: `pi` comes from the equation-(9) Kummer/Fourier phase, the
Riemann-Siegel theta normalization, and the common `pi/8` reflection.  No
fitted or geometric parameter is introduced.

Proof boundary: this gate rejects a false equality and certifies the tiny
Gamma-to-classical carrier defect at `t=10^10`.  It proves no numerical bound
for `Q_K-T`, no complete B face estimate, no A-fold splice, no complete
source-hybrid error theorem, no height-uniform theorem, no `Lambda<=0`, no
PF-infinity, no RH, and no prize-level conclusion.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    ctx.dps = PRECISION
    ctx.threads = 1
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")

    equation10 = dependencies["equation10_current"]
    require(equation10["paper_audit"]["exact_kummer_integral_precedes_asymptotic_grouping"] is True, "Kummer ordering audit drift")
    require(equation10["paper_audit"]["paper_drops_sine_term_before_equation51"] is True, "sine omission audit drift")
    require(equation10["paper_audit"]["equation62_and_hybrid_are_not_exact"] is True, "hybrid exactness audit drift")
    require(equation10["route_decision"]["equation10_sine_term_is_part_of_exact_pair"] is True, "pair-current audit drift")

    require(dependencies["symmetric_endpoint"]["decision"]["bulk_extracted_joint_residual_identity_proved"] is True, "endpoint identity drift")
    require(dependencies["paired_residual"]["decision"]["finite_paired_half_domain_target_residual_identity_proved"] is True, "projection identity drift")
    require(dependencies["double_Weber"]["decision"]["entire_finite_equation9_source_roster_reconstructed_exactly"] is True, "source reconstruction drift")
    require(dependencies["double_Weber"]["decision"]["small_source_minus_target_bound_proved"] is False, "unexpected source-target closure")

    source_row = central_row(dependencies["source_aligned_ledger"]["output_rows"])
    classical_row = central_row(dependencies["classical_upper_split"]["output_rows"])
    source_gap = arb(source_row["source_aligned_upper_block_sum_ball"])
    independent_gap = arb(classical_row["hybrid_minus_classical_main_ball"])
    require(source_gap.overlaps(independent_gap), "source discrepancy certificates do not overlap")
    require(source_gap < 0, "central source discrepancy sign drift")
    ratio = abs(source_gap) / REFERENCE_PHYSICAL_TARGET
    require(ratio > arb(800), "source telemetry ratio drift")

    gamma_bound = arb(dependencies["Gamma_bulk"]["certificate"]["aggregate_two_real_absolute_bound_ball"])
    require(gamma_bound < arb("3e-19"), "Gamma-to-classical aggregate drift")

    artifact = {
        "kind": STEM,
        "status": "exact_Kummer_endpoint_residual_and_source_hybrid_telemetry_distinguished_false_rejection_shortcut_closed",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "source_roster": [159_577, 5_122_421],
            "classical_target_modes": [622, 39_894],
            "output_index": 8,
        },
        "symbolic_certificate": symbolic_certificate(),
        "numerical_context": {
            "source_hybrid_minus_classical_target_ball": source_gap.str(PRECISION, more=True),
            "independent_source_hybrid_minus_classical_target_ball": independent_gap.str(PRECISION, more=True),
            "reference_physical_target": str(REFERENCE_PHYSICAL_TARGET),
            "source_telemetry_to_reference_target_ratio_lower": ratio.lower().str(40, more=True),
            "transfer_to_exact_endpoint_residual_permitted": False,
        },
        "Gamma_bridge": {
            "identity": "R_end^Gamma=(Q_K-T)-(G-T)",
            "aggregate_two_real_absolute_bound_ball": gamma_bound.str(PRECISION, more=True),
            "same_common_projection_required": True,
        },
        "decision": {
            "exact_endpoint_object_is_finite_Kummer_roster": True,
            "double_Weber_reconstructs_implemented_source_hybrid": False,
            "source_hybrid_telemetry_equals_exact_endpoint_residual": False,
            "source_telemetry_lower_bound_transfer_permitted": False,
            "analytic_Kummer_to_published_hybrid_defect_bounded_here": False,
            "complete_published_to_implemented_hybrid_defect_bounded_here": False,
            "Gamma_bulk_to_classical_target_defect_bounded": True,
            "exact_endpoint_residual_quantitatively_bounded_here": False,
            "false_rejection_shortcut_closed": True,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "flint_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Bound the exact projected Q_K-T residual directly, or derive complete signed bounds for both Q_P-Q_K and Q_S-Q_P before using source telemetry. For the current endpoint route, keep the A, B, zero, negative, outer-positive, and endpoint-half sectors together.",
        "proof_boundary": "Exact object-separation and Gamma-carrier bridge at t=10^10 only. No quantitative exact endpoint residual, complete source-hybrid bridge, complete T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "certified endpoint/source nonidentification guard: "
        f"telemetry ratio>{ratio.lower()}, Gamma bridge<{gamma_bound.upper()}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
