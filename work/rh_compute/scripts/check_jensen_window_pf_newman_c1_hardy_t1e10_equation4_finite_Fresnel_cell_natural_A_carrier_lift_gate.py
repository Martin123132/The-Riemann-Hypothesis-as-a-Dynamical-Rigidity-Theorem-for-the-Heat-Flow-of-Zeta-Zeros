#!/usr/bin/env python3
"""Independently check the natural A endpoint-integral carrier lift."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Any

import mpmath as mp


ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_cell_natural_A_carrier_lift_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def independent_integral(t: mp.mpf, endpoint: int, mode: int) -> mp.mpc:
    m = mp.mpf(mode)
    left = -mp.mpf("0.51") + m / 120
    right = left + mp.mpf("0.58")
    rho = -mp.mpf("0.61") - m / 250
    face_zero = mp.mpf("0.27") + endpoint / 900
    lam = (3 * mp.pi * endpoint * m - 1j) / (
        2 * (mp.pi * endpoint * m - 1j) * mp.sqrt(mp.pi * endpoint * m)
    )
    mu = mp.mpf(5) / 12 * mp.sqrt(2 / t)

    def outer(S: mp.mpf) -> mp.mpc:
        face = face_zero + rho * S

        def inner(P: mp.mpf) -> mp.mpc:
            return (1 + lam * P + mu * S) * mp.exp(1j * (P**2 - S**2) / 2)

        return mp.quad(inner, [-mp.mpf("0.93"), face])

    return mp.quad(outer, [left, right]) / (2 * mp.pi)


def independent_row(t: mp.mpf, endpoint: int, mode: int, H: mp.mpf, theta: mp.mpf) -> dict[str, mp.mpf]:
    pi = mp.pi
    m = mp.mpf(mode)
    s = mp.mpf("0.5") + 1j * t
    c_t = (pi / (32 * t)) ** mp.mpf("0.25")
    r = t / (2 * pi * m**2)
    W0 = 4 * (pi * endpoint * m - 1j) * r ** mp.mpf("0.75") / (
        endpoint * mp.sqrt(pi * t)
    )
    theta_zero = t * (mp.log(t / (2 * pi)) - 1) / 2 - pi / 8
    full_phase = pi * m * endpoint - pi * m**2 - t / 2 + t * mp.log(t / (2 * pi * m**2)) / 2
    C = independent_integral(t, endpoint, mode)
    raw = -mp.exp(1j * full_phase) * W0 * C
    X = c_t * mp.exp(-1j * pi / 8) * raw
    eta = -c_t * mp.sqrt(m) * mp.exp(1j * theta_zero) * W0 * C
    coefficient = H * mp.exp(-1j * theta) * eta
    natural = mp.exp(-s * mp.log(m)) * coefficient
    physical = 2 * mp.re(X)
    canonical = H * mp.exp(-1j * theta) * physical / 2
    return {
        "phase": abs(mp.exp(1j * full_phase) - mp.exp(1j * (theta_zero + pi / 8 - t * mp.log(m)))),
        "factor": abs(X - mp.exp(-s * mp.log(m)) * eta),
        "hardy": abs(2 * mp.re(mp.exp(1j * theta) * natural) - H * physical),
        "null": abs(2 * mp.re(mp.exp(1j * theta) * (natural - canonical))),
        "nontrivial": abs(natural - canonical),
    }


def main() -> int:
    set_low_priority()
    mp.mp.dps = 105
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "saved gate is not passed")
    require(NOTE.is_file(), "companion note is missing")
    require(
        artifact["decision"]["natural_A_endpoint_integral_Hardy_lift_exact"] is True,
        "natural lift decision drift",
    )
    require(
        artifact["decision"]["modewise_full_endpoint_complex_enclosures_certified"] is False,
        "modewise enclosure boundary drift",
    )
    for dependency in artifact["dependencies"].values():
        path = ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], "dependency hash drift")
    rows_path = ROOT / artifact["witness_rows"]["path"]
    require(file_hash(rows_path) == artifact["witness_rows"]["sha256"], "witness-row hash drift")
    for source in artifact["sources"].values():
        path = ROOT / source["path"]
        require(path.is_file() and file_hash(path) == source["sha256"], "source hash drift")

    t = mp.mpf("31.625")
    endpoint = 23
    H = mp.mpf("0.9875")
    theta = mp.siegeltheta(t)
    tolerance = mp.mpf("1e-82")
    for mode in (1, 3, 6):
        row = independent_row(t, endpoint, mode, H, theta)
        require(row["phase"] < tolerance, "independent odd-endpoint phase parity failed")
        require(row["factor"] < tolerance, "independent m^-s factorization failed")
        require(row["hardy"] < tolerance, "independent natural Hardy lift failed")
        require(row["null"] < tolerance, "independent Hardy-null comparison failed")
        require(row["nontrivial"] > mp.mpf("1e-8"), "natural lift collapsed to canonical lift")

    expected = artifact["exact_reduction"]
    require(
        expected["transition_cell"] == "m^(-s)[beta_m-C_G-C_A,m] for 39853<=m<=39936",
        "transition coefficient formula drift",
    )
    print("independently checked natural A endpoint-integral carrier lift", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
