#!/usr/bin/env python3
"""Independently check the grouped complex A-endpoint lift."""

from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "work" / "rh_compute" / "vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_transition_grouped_complex_A_endpoint_lift_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
MODES = tuple(range(39_853, 39_937))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def complex_from_record(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def load_rows(path: Path, field: str) -> dict[int, acb]:
    rows: dict[int, acb] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        row = json.loads(line)["row"]
        rows[int(row["mode"])] = complex_from_record(row[field])
    require(tuple(sorted(rows)) == MODES, "independent row roster drift")
    return rows


def set_low_priority() -> None:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
        else:
            process.nice(10)
    except Exception:
        pass


def main() -> int:
    set_low_priority()
    ctx.dps = 120
    ctx.threads = 1
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True and NOTE.is_file(), "saved gate or note missing")
    require(
        artifact["decision"]["grouped_full_complex_A_endpoint_preprojection_certified"] is True,
        "grouped complex decision drift",
    )
    for dependency in artifact["dependencies"].values():
        path = ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], "dependency hash drift")
    for row_source in artifact["row_sources"].values():
        path = ROOT / row_source["path"]
        require(path.is_file() and file_hash(path) == row_source["sha256"], "row hash drift")
    for source in artifact["sources"].values():
        path = ROOT / source["path"]
        require(path.is_file() and file_hash(path) == source["sha256"], "source hash drift")

    deps = {
        name: json.loads((ROOT / row["path"]).read_text(encoding="utf-8"))
        for name, row in artifact["dependencies"].items()
    }
    affine_path = ROOT / artifact["row_sources"]["local_affine"]["path"]
    nonlinear_path = ROOT / artifact["row_sources"]["local_nonlinear"]["path"]
    affine = load_rows(affine_path, "raw_A_endpoint_local_exact_affine_ball")
    nonlinear = load_rows(nonlinear_path, "raw_A_endpoint_local_exact_minus_affine_ball")
    local_affine = sum((affine[m] for m in reversed(MODES)), acb(0))
    local_nonlinear = sum((nonlinear[m] for m in reversed(MODES)), acb(0))
    exterior = complex_from_record(deps["exterior"]["certificate"]["canonical_rational_integral_ball"])
    error = arb(
        deps["exterior_remainder"]["certificate"]
        ["canonical_complete_exterior_replacement_error_ball"]
    )
    exterior = acb(arb(exterior.real, error), arb(exterior.imag, error))
    raw = local_affine + local_nonlinear + exterior

    pi, t, imaginary = arb.pi(), arb(10_000_000_000), acb(0, 1)
    c_t = (pi / (32 * t)) ** arb("0.25")
    X = c_t * (-imaginary * pi / 8).exp() * raw
    physical = 2 * X.real
    saved = artifact["certificate"]["complex_preprojection"]
    require(X.overlaps(complex_from_record(saved["X_A_endpoint_ball"])), "independent X endpoint drift")
    require(physical.overlaps(arb(saved["physical_A_endpoint_reconstruction_ball"])), "independent physical drift")

    H = arb(deps["exact_H"]["exact_H_certificate"]["high_precision_H_ball"])
    theta_zero = t * ((t / (2 * pi)).log() - 1) / 2 - pi / 8
    correction = arb(deps["gamma_bulk"]["certificate"]["exact_theta_minus_theta_zero_ball"])
    phase = (imaginary * (theta_zero + correction)).exp()
    natural = H * phase.conjugate() * X
    hardy = 2 * (phase * natural).real
    require(natural.overlaps(complex_from_record(saved["natural_Hardy_lift_sum_ball"])), "independent natural lift drift")
    require(hardy.overlaps(H * physical), "independent Hardy reconstruction failed")
    require(
        artifact["certificate"]["transition_block_reduction"]
        ["modewise_full_complex_rows_required_for_group_bound"] is False,
        "grouped-bound boundary drift",
    )
    print("independently checked grouped complex A endpoint lift", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
