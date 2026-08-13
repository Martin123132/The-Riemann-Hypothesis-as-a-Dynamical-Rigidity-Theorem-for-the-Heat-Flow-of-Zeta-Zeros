#!/usr/bin/env python3
"""Independent replay of the exact bi-Morse face-fold chart."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_exact_bimorse_face_fold_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")

    # Replay the endpoint Morse map using the hyperbolic parameter eta.
    eta, m, endpoint = sp.symbols("eta m endpoint", real=True, positive=True)
    u = sp.exp(eta)
    z = 2 * m * u / endpoint
    phase = sp.pi * m**2 / z + sp.pi * endpoint**2 * z / 4
    p = 2 * sp.sinh(eta / 2)
    expected = sp.pi * m * endpoint + sp.pi * m * endpoint * p**2 / 2
    require(sp.simplify(sp.expand_trig(phase.rewrite(sp.exp)) - expected.rewrite(sp.exp)) == 0, "hyperbolic Morse replay failed")

    t = sp.symbols("t", positive=True, real=True)
    alpha = 2 * m + t / (sp.pi * m)
    event_t = sp.pi * m * (endpoint - 2 * m)
    require(sp.simplify(alpha.subs(t, event_t) - endpoint) == 0, "event replay failed")
    rho2_event = 2 * event_t / (sp.pi * m * endpoint)
    require(sp.simplify(rho2_event - 1 - (1 - 4 * m / endpoint)) == 0, "null defect replay failed")
    require(sp.simplify(event_t.subs(m, endpoint / 4) - sp.pi * endpoint**2 / 8) == 0, "corner replay failed")

    rows = {(row["endpoint"], row["mode"]): row for row in artifact["finite_height_ledger"]["rows"]}
    require(float(rows[(5_122_421, 621)]["event_face_fold_defect_1_minus_4m_over_D"]) > 0.999, "B defect replay failed")
    require(abs(float(rows[(159_577, 39_894)]["event_face_fold_defect_1_minus_4m_over_D"]) - 1 / 159_577) < 1e-18, "A defect replay failed")

    decision = artifact["decision"]
    require(decision["endpoint_phase_has_exact_global_Morse_coordinate"] is True, "endpoint Morse decision drift")
    require(decision["pair_triangle_phase_is_exactly_hyperbolic_quadratic"] is True, "bi-Morse decision drift")
    require(decision["A_fold_and_half_boundary_are_same_null_face"] is True, "null-face decision drift")
    require(decision["uniform_curved_face_remainder_proved"] is False, "uniform-bound overclaim")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    print("checked exact bi-Morse face-fold chart by hyperbolic replay", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
