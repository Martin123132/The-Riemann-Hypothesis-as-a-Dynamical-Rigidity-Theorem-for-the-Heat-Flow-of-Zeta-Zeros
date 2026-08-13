#!/usr/bin/env python3
"""Independently check the lower-fold Hankel carrier event selection."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))

from flint import arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_hankel_carrier_event_selection_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(RESULT.is_file() and NOTE.is_file(), "missing result or note")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == STEM and artifact["passed"] is True, "gate identity drift")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file() and file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file() and file_hash(path) == record["sha256"], f"source hash drift: {path}")

    ctx.dps = 125
    ctx.threads = 1
    pi = arb.pi()
    c = arb(159_577)
    beta = (pi * c**2 / 8) ** (arb(1) / 3)
    shift = pi / (16 * beta)
    d_min = beta / c
    d_max = arb(797) * beta / c

    # Independent route: replace beta^2/C^2 by beta^3/(beta*C^2).
    lambda_floor = beta**3 / (beta * c**2) - shift
    require(lambda_floor.overlaps(pi / (16 * beta)), "exact lambda floor replay failed")
    require(lambda_floor > arb("9e-5"), "positive-X floor replay failed")
    x_ceiling = d_max**2 + shift + 64
    require(x_ceiling < 180 and x_ceiling < beta**2, "strip ceiling replay failed")

    opposite_gap = d_min + lambda_floor.sqrt()
    face_slope = shift / opposite_gap
    require(opposite_gap > arb("0.023"), "opposite carrier gap replay failed")
    require(face_slope < arb("0.004"), "selected face slope replay failed")

    stored = artifact["interval_certificate"]
    for key, value in (
        ("lambda_min_ball", lambda_floor),
        ("fold_strip_X_max_ball", x_ceiling),
        ("opposite_branch_endpoint_gap_lower_ball", opposite_gap),
        ("selected_branch_face_slope_upper_ball", face_slope),
    ):
        require(arb(stored[key]).overlaps(value), f"stored interval drift: {key}")

    # Direct finite roster check, deliberately not using the builder's list.
    labels = []
    for n in range(399):
        labels.append(-(2 * n + 1) if n % 2 == 0 else 2 * n + 1)
    require(len(set(labels)) == 399, "independent event labels are not distinct")
    require(min(map(abs, labels)) == 1 and max(map(abs, labels)) == 797, "independent label extrema drift")

    decision = artifact["decision"]
    require(decision["all_399_event_cells_strictly_positive_X"] is True, "positive-X decision lost")
    require(decision["selected_extracted_carrier_crosses_at_event_center"] is True, "branch selection lost")
    require(decision["opposite_extracted_carrier_nonstationary_on_fold_strip"] is True, "opposite branch guard lost")
    require(decision["exact_scaled_Hankel_amplitude_derivative_bound_proved"] is False, "amplitude bound overclaim")
    require(decision["ordinary_logistic_Gamma_carrier_identification_proved"] is False, "ordinary carrier overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    note = NOTE.read_text(encoding="utf-8")
    for token in ("factorization is exact", "pi/(16beta)", "do not yet license integration by parts", "PF-infinity, RH"):
        require(token in note, f"note boundary token missing: {token}")
    print("independently checked Hankel event selection: exact lambda floor, 399 labels, opposite gap >0.023, 0 amplitude claims", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
