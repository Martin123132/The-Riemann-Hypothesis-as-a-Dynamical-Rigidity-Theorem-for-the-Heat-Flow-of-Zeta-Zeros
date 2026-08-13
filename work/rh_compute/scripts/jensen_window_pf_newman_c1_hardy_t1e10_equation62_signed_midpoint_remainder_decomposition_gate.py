#!/usr/bin/env python3
"""Resolve the finite diagnostic-midpoint residual into signed channels."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_block_zero_equation62_cell_partition_gate as cells
import jensen_window_pf_newman_c1_hardy_t1e10_equation62_five_saddle_D_transport_gate as dgate


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation62_signed_midpoint_remainder_decomposition_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation62_signed_midpoint_remainder_decomposition_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 110
EXPECTED_OUTPUTS = 15
REQUESTED_TOLERANCE = arb("0.005")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


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


def load_artifact(path: Path, expected_kind: str) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    artifact = json.loads(path.read_text(encoding="utf-8"))
    require(artifact.get("kind") == expected_kind, f"dependency kind drift: {path}")
    require(bool(artifact.get("passed")), f"dependency gate failed: {path}")
    return artifact


def derive_row(cell_row: dict[str, Any], d_row: dict[str, Any]) -> dict[str, Any]:
    require(cell_row["output_index"] == d_row["output_index"], "output index drift")
    require(cell_row["output_label"] == d_row["output_label"], "output label drift")
    require(cell_row["target_t"] == d_row["target_t"], "target height drift")

    final = cell_row["final_prefix"]
    source_remainder = arb(final["prefix_minus_source_endpoint_target_ball"])
    midpoint_remainder = arb(final["prefix_minus_midpoint_ceil_target_ball"])
    exact_classical_collar = arb(final["midpoint_ceil_target_minus_source_endpoint_target_ball"])
    transported_D = arb(d_row["transformed_D_collar_real_ball"])
    d_exact_classical_collar = arb(d_row["exact_theta_classical_collar_ball"])
    require(exact_classical_collar.overlaps(d_exact_classical_collar), "five-term collar handoff drift")

    complementary_channel = source_remainder - transported_D
    D_boundary_channel = transported_D - exact_classical_collar
    signed_reconstruction = complementary_channel + D_boundary_channel
    direct_midpoint_identity = source_remainder - exact_classical_collar
    independent_triangle = abs(complementary_channel) + abs(D_boundary_channel)
    retained_fraction = abs(signed_reconstruction) / independent_triangle
    cancellation_fraction = 1 - retained_fraction

    require(complementary_channel < 0, "complementary channel sign drift")
    require(D_boundary_channel > 0, "D-boundary channel sign drift")
    require(midpoint_remainder > 0, "midpoint remainder sign drift")
    require(signed_reconstruction.overlaps(midpoint_remainder), "signed reconstruction misses midpoint remainder")
    require(direct_midpoint_identity.overlaps(midpoint_remainder), "direct midpoint identity drift")
    require(abs(midpoint_remainder) < REQUESTED_TOLERANCE, "midpoint target no longer closes")
    require(independent_triangle > REQUESTED_TOLERANCE, "independent triangle unexpectedly closes")

    return {
        "output_index": cell_row["output_index"],
        "output_label": cell_row["output_label"],
        "target_t": cell_row["target_t"],
        "source_endpoint_remainder_ball": source_remainder.str(PRECISION, more=True),
        "transported_D_real_ball": transported_D.str(PRECISION, more=True),
        "exact_five_term_classical_collar_ball": exact_classical_collar.str(PRECISION, more=True),
        "complementary_channel_source_minus_D_ball": complementary_channel.str(PRECISION, more=True),
        "D_minus_exact_classical_collar_ball": D_boundary_channel.str(PRECISION, more=True),
        "signed_channel_sum_ball": signed_reconstruction.str(PRECISION, more=True),
        "direct_midpoint_remainder_ball": midpoint_remainder.str(PRECISION, more=True),
        "signed_sum_minus_direct_midpoint_ball": (signed_reconstruction - midpoint_remainder).str(PRECISION, more=True),
        "independent_absolute_triangle_ball": independent_triangle.str(PRECISION, more=True),
        "signed_remainder_retained_fraction_ball": retained_fraction.str(PRECISION, more=True),
        "opposite_sign_cancellation_fraction_ball": cancellation_fraction.str(PRECISION, more=True),
        "complementary_channel_strict_negative": True,
        "D_boundary_channel_strict_positive": True,
        "midpoint_remainder_strict_positive_below_0p005": True,
        "independent_triangle_above_0p005": True,
    }


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    center = artifact["output_rows"][7]
    return f"""# Signed diagnostic-midpoint decomposition gate

Date: 2026-08-09

Status: finite certified signed decomposition; not a proof of RH

Let `E_src=Q-T_src` be the final equation-(62) residual at the published/source
endpoint, let `T_5=T_mid-T_src` be the exact-theta five-term change to the
introduced diagnostic midpoint, and let `Delta_D` be the certified real
five-saddle B9 transport.  Pure algebra gives

```text
E_mid = E_src-T_5
      = (E_src-Re Delta_D) + (Re Delta_D-T_5).
```

At the central output the two channels and their sum are

```text
E_src-Re Delta_D = {center['complementary_channel_source_minus_D_ball']}
Re Delta_D-T_5   = {center['D_minus_exact_classical_collar_ball']}
signed sum        = {center['signed_channel_sum_ball']}
```

Across all fifteen outputs, the complementary channel is strictly negative in

```text
{aggregate['minimum_complementary_channel_ball']}
{aggregate['maximum_complementary_channel_ball']}
```

while the D-to-boundary channel is strictly positive in

```text
{aggregate['minimum_D_boundary_channel_ball']}
{aggregate['maximum_D_boundary_channel_ball']}
```

Their signed sum reproduces the independently computed diagnostic midpoint residual in

```text
{aggregate['minimum_signed_midpoint_remainder_ball']}
{aggregate['maximum_signed_midpoint_remainder_ball']}
```

and is below `0.005` at every saved height.  In contrast, bounding the two
channels independently gives an absolute triangle between
`{aggregate['minimum_independent_triangle_ball']}` and
`{aggregate['maximum_independent_triangle_ball']}`, above `0.005` everywhere.
Only about 15 percent of that unsigned size survives the opposite-sign
cancellation.

This finite gate identifies correlation in the diagnostic alternative.
It does not prove that either sign persists for every admissible radius or
height.  More importantly, the `D` split is auxiliary: while the source and
diagnostic-midpoint residue rosters stay fixed, changing a contour selector gives
`dA=-d(Re Delta_D)`, `dB=d(Re Delta_D)`, and hence `d(A+B)=0`.  This does not
make `E_mid` the paper's intrinsic error: the midpoint cutoff itself still
needs an independent alternative-hybrid theorem.  The primary next theorem is
a direct signed bound for the source-aligned equations-(126)--(127) error.
Selector-uniform channel signs matter only if a separate midpoint theorem uses
this contour split.  There is no implication for `Lambda<=0`, PF-infinity, RH,
or the Clay prize.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    ctx.dps = PRECISION
    ctx.threads = 1
    cell_artifact = load_artifact(
        cells.RESULT,
        "jensen_window_pf_newman_c1_hardy_t1e10_block_zero_equation62_cell_partition_gate",
    )
    d_artifact = load_artifact(
        dgate.RESULT,
        "jensen_window_pf_newman_c1_hardy_t1e10_equation62_five_saddle_D_transport_gate",
    )
    require(len(cell_artifact["output_rows"]) == EXPECTED_OUTPUTS, "cell output roster drift")
    require(len(d_artifact["output_rows"]) == EXPECTED_OUTPUTS, "D output roster drift")

    rows = [derive_row(c, d) for c, d in zip(cell_artifact["output_rows"], d_artifact["output_rows"])]
    complementary = [arb(row["complementary_channel_source_minus_D_ball"]) for row in rows]
    D_boundary = [arb(row["D_minus_exact_classical_collar_ball"]) for row in rows]
    midpoint = [arb(row["signed_channel_sum_ball"]) for row in rows]
    triangles = [arb(row["independent_absolute_triangle_ball"]) for row in rows]
    retained = [arb(row["signed_remainder_retained_fraction_ball"]) for row in rows]
    cancelled = [arb(row["opposite_sign_cancellation_fraction_ball"]) for row in rows]

    aggregate = {
        "output_count": EXPECTED_OUTPUTS,
        "opposite_sign_channel_count": sum(a < 0 and b > 0 for a, b in zip(complementary, D_boundary)),
        "signed_midpoint_below_0p005_count": sum(abs(value) < REQUESTED_TOLERANCE for value in midpoint),
        "independent_triangle_above_0p005_count": sum(value > REQUESTED_TOLERANCE for value in triangles),
        "minimum_complementary_channel_ball": min(complementary, key=lambda x: x.lower()).str(PRECISION, more=True),
        "maximum_complementary_channel_ball": max(complementary, key=lambda x: x.upper()).str(PRECISION, more=True),
        "minimum_D_boundary_channel_ball": min(D_boundary, key=lambda x: x.lower()).str(PRECISION, more=True),
        "maximum_D_boundary_channel_ball": max(D_boundary, key=lambda x: x.upper()).str(PRECISION, more=True),
        "minimum_signed_midpoint_remainder_ball": min(midpoint, key=lambda x: x.lower()).str(PRECISION, more=True),
        "maximum_signed_midpoint_remainder_ball": max(midpoint, key=lambda x: x.upper()).str(PRECISION, more=True),
        "minimum_independent_triangle_ball": min(triangles, key=lambda x: x.lower()).str(PRECISION, more=True),
        "maximum_independent_triangle_ball": max(triangles, key=lambda x: x.upper()).str(PRECISION, more=True),
        "minimum_signed_remainder_retained_fraction_ball": min(retained, key=lambda x: x.lower()).str(PRECISION, more=True),
        "maximum_signed_remainder_retained_fraction_ball": max(retained, key=lambda x: x.upper()).str(PRECISION, more=True),
        "minimum_opposite_sign_cancellation_fraction_ball": min(cancelled, key=lambda x: x.lower()).str(PRECISION, more=True),
        "maximum_opposite_sign_cancellation_fraction_ball": max(cancelled, key=lambda x: x.upper()).str(PRECISION, more=True),
    }
    require(aggregate["opposite_sign_channel_count"] == EXPECTED_OUTPUTS, "opposite-sign aggregate drift")
    require(aggregate["signed_midpoint_below_0p005_count"] == EXPECTED_OUTPUTS, "signed target aggregate drift")
    require(aggregate["independent_triangle_above_0p005_count"] == EXPECTED_OUTPUTS, "triangle aggregate drift")

    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation62_signed_midpoint_remainder_decomposition_gate",
        "status": "finite_diagnostic_midpoint_residual_resolved_as_opposite_sign_channels",
        "passed": True,
        "scope": {
            "height_center": "1e10",
            "outputs": EXPECTED_OUTPUTS,
            "precision_decimal_digits": PRECISION,
            "classical_collar": list(dgate.COLLAR),
        },
        "exact_identity": {
            "definitions": "E_src=Q-T_src, T_5=T_mid-T_src, E_mid=Q-T_mid",
            "decomposition": "E_mid=(E_src-ReDelta_D)+(ReDelta_D-T_5)=E_src-T_5",
            "fixed_roster_selector_derivative": "dA=-d(ReDelta_D), dB=d(ReDelta_D), d(A+B)=0",
        },
        "output_rows": rows,
        "aggregate": aggregate,
        "decision": {
            "source_remainder_minus_D_strict_negative_all_outputs": True,
            "D_minus_exact_classical_collar_strict_positive_all_outputs": True,
            "signed_sum_equals_midpoint_remainder_all_outputs": True,
            "midpoint_remainder_below_0p005_all_outputs": True,
            "separate_absolute_triangle_bound_closes_target": False,
            "selector_split_changes_diagnostic_midpoint_residual": False,
            "channel_sign_uniformity_is_source_aligned_next_theorem": False,
            "diagnostic_midpoint_is_a_published_cutoff_rule": False,
            "height_uniform_sign_pairing_proved": False,
            "finite_gate_has_rh_implication": False,
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "flint_threads": 1,
            "process_priority": priority,
        },
        "dependencies": {
            "cell_partition_gate": {"path": relative(cells.RESULT), "sha256": file_hash(cells.RESULT)},
            "five_saddle_D_transport_gate": {"path": relative(dgate.RESULT), "sha256": file_hash(dgate.RESULT)},
            "lewis_2015_paper": {"path": relative(dgate.collar.LEGACY_PAPER), "sha256": file_hash(dgate.collar.LEGACY_PAPER)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": "Derive a direct signed, height-uniform error bound for the source-aligned equations-(126)--(127) hybrid. Pursue the diagnostic midpoint only through a separately derived alternative-cutoff theorem that preserves the observed cell cancellation.",
        "proof_boundary": "Rigorous finite signed decomposition of an introduced diagnostic comparison at fifteen saved heights only. It proves an algebraic cancellation identity but neither validates the midpoint cutoff nor bounds the source-aligned hybrid uniformly in height, proves Lambda<=0, PF-infinity, RH, or a prize-level result.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "validated signed midpoint remainder decomposition: "
        f"outputs={EXPECTED_OUTPUTS}, opposite-sign={aggregate['opposite_sign_channel_count']}, "
        f"signed<0.005={aggregate['signed_midpoint_below_0p005_count']}, "
        f"triangle>0.005={aggregate['independent_triangle_above_0p005_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
