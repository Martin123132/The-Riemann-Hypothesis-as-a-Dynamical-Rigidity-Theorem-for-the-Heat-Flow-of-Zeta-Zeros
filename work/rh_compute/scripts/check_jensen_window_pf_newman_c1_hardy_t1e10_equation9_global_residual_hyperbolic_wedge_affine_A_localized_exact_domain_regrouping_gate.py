#!/usr/bin/env python3
"""Independently check the localized exact A-domain regrouping gate."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_A_localized_exact_domain_regrouping_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing file: {path}")
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


def main() -> int:
    priority = set_low_priority()
    require(priority == "below_normal" or os.name != "nt", f"checker priority drift: {priority}")
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True and artifact.get("kind") == STEM, "artifact identity drift")
    decision = artifact["decision"]
    for key in (
        "global_tangent_plus_full_face_strip_equals_exact_affine_domain",
        "global_tangent_plus_local_face_retains_tangent_exterior",
        "tangent_exterior_must_be_removed_before_estimation",
        "eleven_artificial_tangent_exterior_stationary_modes_identified",
        "exact_face_phase_is_common_after_carrier_restoration",
        "ten_true_exact_face_stationary_transitions_remain_in_exterior",
        "uniform_exterior_inner_Fresnel_endpoint_below_minus_261p2",
    ):
        require(decision.get(key) is True, f"missing decision: {key}")
    require(decision.get("blanket_exterior_nonstationary_IBP_is_valid") is False, "IBP guard drift")
    for key in (
        "exterior_exact_affine_current_bound_proved",
        "transformed_amplitude_remainder_bound_proved",
        "complete_A_endpoint_block_proved",
        "R_Dir_bound_proved",
        "rh_implication",
    ):
        require(decision.get(key) is False, f"proof-boundary drift: {key}")

    ctx.dps = 120
    ctx.threads = 1
    pi, t, endpoint = arb.pi(), arb(10_000_000_000), arb(159_577)
    y0 = arb("0.0037")
    S0 = (t * (y0 - (1 + y0).log())).sqrt()
    kappa = (1 - 8 * t / (pi * endpoint**2)).sqrt()
    x_star = (1 - kappa) / 2
    fresh_tangent: list[int] = []
    fresh_exact: list[int] = []
    rows = {row["mode"]: row for row in artifact["interval_certificate"]["rows"]}
    require(set(rows) == set(range(39_853, 39_937)), "saved row roster drift")

    for mode_int in range(39_853, 39_937):
        mode = arb(mode_int)
        r = t / (2 * pi * mode**2)
        alpha = 2 * mode + t / (pi * mode)
        a = (pi * mode / alpha).sqrt() * (endpoint - alpha)
        denominator = 2 * pi * mode**2 + t
        rho = -(2 * t).sqrt() * (pi * endpoint * mode + denominator) / (2 * denominator ** arb("1.5"))
        c = rho**2 - 1
        root = None
        if c.lower() > 0:
            root = -(rho * a) / c
            if root.lower() > S0.upper():
                fresh_tangent.append(mode_int)
        y_star = (1 / x_star - 1) / r - 1
        if y_star.lower() > y0.upper():
            fresh_exact.append(mode_int)

        saved = rows[mode_int]
        for field, value in (
            ("a_ball", a),
            ("rho_ball", rho),
            ("c_ball", c),
            ("rho_a_ball", rho * a),
            ("exact_face_stationary_y_ball", y_star),
        ):
            require(arb(saved[field]).overlaps(value), f"fresh {field} misses at mode {mode_int}")
        if root is None:
            require(saved["tangent_stationary_S_ball"] is None, f"unexpected saved tangent root at {mode_int}")
        else:
            require(arb(saved["tangent_stationary_S_ball"]).overlaps(root), f"fresh tangent root misses at {mode_int}")

        g0 = 1 + r * (1 + y0)
        p0 = (pi / 2).sqrt() * (endpoint / g0.sqrt() - 2 * mode * g0.sqrt())
        require(arb(saved["P_exact_at_y0_ball"]).overlaps(p0), f"fresh P cutoff misses at mode {mode_int}")
        require(p0.upper() < arb("-261.2"), f"fresh P cutoff margin fails at mode {mode_int}")

    require(fresh_tangent == list(range(39_884, 39_895)), "fresh tangent exterior roster drift")
    require(fresh_exact == list(range(39_927, 39_937)), "fresh exact exterior roster drift")
    certificate = artifact["interval_certificate"]
    require(arb(certificate["S_cutoff_ball"]).overlaps(S0), "S cutoff overlap failed")
    require(arb(certificate["common_exact_face_kappa_ball"]).overlaps(kappa), "kappa overlap failed")
    require(arb(certificate["common_exact_face_stationary_x_ball"]).overlaps(x_star), "x-star overlap failed")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    require(NOTE.is_file(), "missing note")
    note = NOTE.read_text(encoding="utf-8")
    require("artificial exterior stationary point" in note, "tangent ghost guard missing")
    require("blanket exterior" in note and "not valid" in note, "nonstationary guard missing")
    require("No local exact-affine value" in note, "proof boundary missing")
    print("independently checked localized exact A domain and exterior saddle rosters", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
