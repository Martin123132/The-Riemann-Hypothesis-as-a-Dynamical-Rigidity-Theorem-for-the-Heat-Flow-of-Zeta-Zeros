#!/usr/bin/env python3
"""Certify the grouped complex A-endpoint lift for the transition block."""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "work" / "rh_compute" / "vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_transition_grouped_complex_A_endpoint_lift_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"

HEIGHT = 10_000_000_000
LOWER_START = 39_853
UPPER_END = 39_936
MODES = tuple(range(LOWER_START, UPPER_END + 1))

NATURAL_LIFT = ROOT / "work" / "rh_compute" / "results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_cell_"
    "natural_A_carrier_lift_gate.json"
)
LOCAL_AFFINE = ROOT / "work" / "rh_compute" / "results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_"
    "hyperbolic_wedge_affine_A_local_exact_domain_gate.json"
)
LOCAL_AFFINE_ROWS = ROOT / "work" / "rh_compute" / "results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_"
    "hyperbolic_wedge_affine_A_local_exact_domain_gate_rows.jsonl"
)
LOCAL_NONLINEAR = ROOT / "work" / "rh_compute" / "results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_"
    "hyperbolic_wedge_A_local_exact_nonlinear_amplitude_gate.json"
)
LOCAL_NONLINEAR_ROWS = ROOT / "work" / "rh_compute" / "results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_"
    "hyperbolic_wedge_A_local_exact_nonlinear_amplitude_gate_rows.jsonl"
)
EXTERIOR = ROOT / "work" / "rh_compute" / "results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_"
    "hyperbolic_wedge_A_exact_exterior_rational_common_phase_contour_gate.json"
)
EXTERIOR_REMAINDER = ROOT / "work" / "rh_compute" / "results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_"
    "hyperbolic_wedge_A_exact_exterior_endpoint_tail_remainder_gate.json"
)
A_ASSEMBLY = ROOT / "work" / "rh_compute" / "results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_"
    "hyperbolic_wedge_A_complete_endpoint_projector_assembly_gate.json"
)
EXACT_H = ROOT / "work" / "rh_compute" / "results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_exact_H_"
    "continuation_defect_target_gate.json"
)
GAMMA_BULK = ROOT / "work" / "rh_compute" / "results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_global_morse_"
    "gamma_bulk_gate.json"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
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


def complex_record(value: acb, digits: int = 70) -> dict[str, str]:
    return {
        "real_ball": value.real.str(digits, more=True),
        "imag_ball": value.imag.str(digits, more=True),
        "absolute_ball": abs(value).str(digits, more=True),
    }


def complex_from_record(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def add_complex_error(value: acb, radius: arb) -> acb:
    return acb(arb(value.real, radius), arb(value.imag, radius))


def load_raw_rows(path: Path, field: str) -> dict[int, acb]:
    rows: dict[int, acb] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        row = json.loads(line).get("row", {})
        mode = row.get("mode")
        if mode in MODES:
            require(mode not in rows, f"duplicate mode {mode} in {path.name}")
            rows[int(mode)] = complex_from_record(row[field])
    require(tuple(sorted(rows)) == MODES, f"incomplete mode roster in {path.name}")
    return rows


def ordered_sum(rows: dict[int, acb], reverse: bool = False) -> acb:
    modes = reversed(MODES) if reverse else MODES
    return sum((rows[mode] for mode in modes), acb(0))


def certificate(dependencies: dict[str, dict[str, Any]]) -> dict[str, Any]:
    pi, t, imaginary = arb.pi(), arb(HEIGHT), acb(0, 1)
    affine_rows = load_raw_rows(LOCAL_AFFINE_ROWS, "raw_A_endpoint_local_exact_affine_ball")
    nonlinear_rows = load_raw_rows(
        LOCAL_NONLINEAR_ROWS, "raw_A_endpoint_local_exact_minus_affine_ball"
    )
    local_affine = ordered_sum(affine_rows)
    local_nonlinear = ordered_sum(nonlinear_rows)
    reverse_affine = ordered_sum(affine_rows, reverse=True)
    reverse_nonlinear = ordered_sum(nonlinear_rows, reverse=True)
    require(local_affine.overlaps(reverse_affine), "affine summation-order drift")
    require(local_nonlinear.overlaps(reverse_nonlinear), "nonlinear summation-order drift")

    exterior_rational = complex_from_record(
        dependencies["exterior"]["certificate"]["canonical_rational_integral_ball"]
    )
    exterior_error = arb(
        dependencies["exterior_remainder"]["certificate"]
        ["canonical_complete_exterior_replacement_error_ball"]
    )
    exterior_exact = add_complex_error(exterior_rational, exterior_error)
    raw_endpoint = local_affine + local_nonlinear + exterior_exact

    c_t = (pi / (32 * t)) ** arb("0.25")
    paper_rotation = (-imaginary * pi / 8).exp()
    X_endpoint = c_t * paper_rotation * raw_endpoint
    physical_endpoint = 2 * X_endpoint.real
    saved_physical = arb(
        dependencies["A_assembly"]["certificate"]["components"]
        ["complete_exact_A_endpoint_carrier_ball"]
    )
    require(physical_endpoint.overlaps(saved_physical), "grouped complex lift misses saved A endpoint")

    H = arb(dependencies["exact_H"]["exact_H_certificate"]["high_precision_H_ball"])
    theta_zero = t * ((t / (2 * pi)).log() - 1) / 2 - pi / 8
    theta_correction = arb(
        dependencies["gamma_bulk"]["certificate"]["exact_theta_minus_theta_zero_ball"]
    )
    theta_phase = (imaginary * (theta_zero + theta_correction)).exp()
    natural_lift = H * theta_phase.conjugate() * X_endpoint
    hardy_reconstruction = 2 * (theta_phase * natural_lift).real
    require(hardy_reconstruction.overlaps(H * physical_endpoint), "natural grouped Hardy lift failed")

    canonical_real_lift = H * theta_phase.conjugate() * physical_endpoint / 2
    hardy_null = 2 * (theta_phase * (natural_lift - canonical_real_lift)).real
    require(hardy_null.contains(0), "grouped natural/canonical difference is not Hardy-null")
    require(abs(natural_lift - canonical_real_lift).lower() > arb("1e-7"), "grouped complex lift collapsed")

    return {
        "height": HEIGHT,
        "mode_range": [LOWER_START, UPPER_END],
        "mode_count": len(MODES),
        "raw_components": {
            "local_exact_affine_sum_ball": complex_record(local_affine),
            "local_exact_minus_affine_sum_ball": complex_record(local_nonlinear),
            "rational_exterior_sum_ball": complex_record(exterior_rational),
            "exterior_complex_replacement_error_radius_ball": exterior_error.str(70, more=True),
            "full_exact_exterior_sum_ball": complex_record(exterior_exact),
            "complete_raw_A_endpoint_sum_ball": complex_record(raw_endpoint),
        },
        "complex_preprojection": {
            "X_A_endpoint_ball": complex_record(X_endpoint),
            "physical_A_endpoint_reconstruction_ball": physical_endpoint.str(70, more=True),
            "saved_physical_A_endpoint_ball": saved_physical.str(70, more=True),
            "natural_Hardy_lift_sum_ball": complex_record(natural_lift),
            "natural_lift_Hardy_reconstruction_ball": hardy_reconstruction.str(70, more=True),
            "canonical_real_lift_sum_ball": complex_record(canonical_real_lift),
            "natural_minus_canonical_lift_absolute_ball": abs(
                natural_lift - canonical_real_lift
            ).str(70, more=True),
            "natural_minus_canonical_Hardy_projection_ball": hardy_null.str(70, more=True),
        },
        "transition_block_reduction": {
            "natural_endpoint_sum": "mathcal A_A^nat=sum_(m=39853)^39936 m^(-s)C_A,m",
            "joined_transition_packet": "sum_(m=39853)^39936 m^(-s)(beta_m-C_G)-mathcal A_A^nat",
            "equivalent_modewise_form": "sum_(m=39853)^39936 m^(-s)(beta_m-C_G-C_A,m)",
            "modewise_full_complex_rows_required_for_group_bound": False,
        },
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    p = c["complex_preprojection"]
    return f"""# Grouped complex A endpoint lift on the transition block

Date: 2026-08-27

Status: rigorous grouped complex endpoint lift; transition Fresnel/Gamma block
not yet evaluated.

The local affine and local exact-minus-affine caches retain the raw complex
endpoint contribution for every mode `39853..39936`.  The exact-exterior gate
retains their grouped complex rational contour value, while the endpoint-tail
gate supplies a rigorous complex modulus allowance for replacing the exact
tail.  Therefore, before physical projection,

```text
R_A^raw=R_A,loc^aff+R_A,loc^nonlinear+R_A,ext^exact,
X_A,endpoint=c_t exp(-i*pi/8)R_A^raw.
```

The certified complex preprojection is

```text
X_A,endpoint = {p['X_A_endpoint_ball']['real_ball']}
             + i {p['X_A_endpoint_ball']['imag_ball']}.
```

It reconstructs the previously certified physical endpoint through

```text
2 Re X_A,endpoint = {p['physical_A_endpoint_reconstruction_ball']}.
```

Using the exact `H(t)` and Riemann-Siegel phase gives the natural grouped lift

```text
mathcal A_A^nat=H(t)exp(-i theta(t))X_A,endpoint
               =sum_(m=39853)^39936 m^(-s)C_A,m,

Hardy_t[mathcal A_A^nat]=H(t)A_endpoint.
```

Thus the entire transition contribution can be bounded as one block:

```text
sum_(m=39853)^39936 m^(-s)(beta_m-C_G)
 -mathcal A_A^nat,
```

which is exactly the modewise sum of `m^(-s)(beta_m-C_G-C_A,m)`.  Separate
full complex exterior integration for 84 modes is not logically required.

Pi provenance: `pi` comes from the exact endpoint phase, bi-Morse
normalization, paper rotation, and Riemann-Siegel phase.  No fitted constant
is introduced.

Proof boundary: rigorous grouped complex A endpoint lift at `t=10^10` only.
This gate does not evaluate or bound the joined Fresnel/Gamma transition
packet, prove the `P_W` Mordell join, enclose `J_Z` or `D_K`, or prove a non-A,
all-height, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    ctx.dps = 100
    ctx.threads = 1
    paths = {
        "natural_lift": NATURAL_LIFT,
        "local_affine": LOCAL_AFFINE,
        "local_nonlinear": LOCAL_NONLINEAR,
        "exterior": EXTERIOR,
        "exterior_remainder": EXTERIOR_REMAINDER,
        "A_assembly": A_ASSEMBLY,
        "exact_H": EXACT_H,
        "gamma_bulk": GAMMA_BULK,
    }
    dependencies = {name: load_json(path) for name, path in paths.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency is not passed")
    require(
        dependencies["natural_lift"]["decision"]
        ["A_endpoint_aligned_on_common_m_minus_s_cell_carrier"] is True,
        "natural-lift dependency drift",
    )
    require(
        dependencies["exterior_remainder"]["decision"]["all_84_mode_remainders_certified"] is True,
        "exterior-remainder dependency drift",
    )
    certified = certificate(dependencies)

    builder = Path(__file__).resolve()
    checker = builder.with_name("check_" + builder.name)
    artifact = {
        "kind": STEM,
        "status": "grouped_full_complex_A_endpoint_natural_Hardy_lift_certified_transition_packet_open",
        "passed": True,
        "certificate": certified,
        "decision": {
            "grouped_full_complex_A_endpoint_preprojection_certified": True,
            "grouped_natural_A_Hardy_lift_certified": True,
            "physical_A_endpoint_reconstructed_from_complex_sum": True,
            "modewise_full_complex_exterior_rows_required_for_transition_bound": False,
            "joined_Fresnel_Gamma_A_transition_packet_bounded": False,
            "P_W_unowned_cell_Mordell_join_proved": False,
            "actual_height_J_Z_enclosed": False,
            "actual_height_D_K_enclosed": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in paths.items()
        },
        "row_sources": {
            "local_affine": {"path": relative(LOCAL_AFFINE_ROWS), "sha256": file_hash(LOCAL_AFFINE_ROWS)},
            "local_nonlinear": {"path": relative(LOCAL_NONLINEAR_ROWS), "sha256": file_hash(LOCAL_NONLINEAR_ROWS)},
        },
        "sources": {
            "builder": {"path": relative(builder), "sha256": file_hash(builder)},
            "checker": {"path": relative(checker), "sha256": file_hash(checker)},
        },
        "runtime": {
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.time() - started, 3),
        },
        "next_obligation": "Evaluate and rigorously enclose the grouped transition packet sum m^(-s)(beta_m-C_G)-mathcal A_A^nat without separating its carriers. Then derive the exact Mordell join of P_W with the lower and upper unowned cells.",
        "proof_boundary": "Rigorous grouped full complex A endpoint preprojection and natural Hardy lift at t=10^10 only. No joined transition packet bound, P_W/unowned-cell Mordell join, actual-height J_Z or D_K enclosure, non-A, all-height, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified grouped complex A endpoint lift", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
