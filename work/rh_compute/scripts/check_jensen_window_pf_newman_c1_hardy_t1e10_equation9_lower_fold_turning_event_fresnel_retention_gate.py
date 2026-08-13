#!/usr/bin/env python3
"""Independently validate the turning-event Fresnel-retention gate."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, REPO_ROOT / "work/rh_compute/scripts"):
    sys.path.insert(0, str(candidate))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, acb, ctx


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_fresnel_retention_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_fresnel_retention_gate.md"
PRECISION = 120
C = 159_577
Q = 39_894
COUNT = 399


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def mode(index: int) -> int:
    return Q - index // 2 if index % 2 == 0 else Q + (index + 1) // 2


def event(index: int, pi: arb) -> arb:
    return pi * C**2 / 8 - pi * (arb(index * (index + 1) // 2) + arb(1) / 8)


def z_coordinate(height: arb, event_height: arb, event_mode: int, pi: arb) -> arb:
    return (event_height - height) / (pi * (height + 2 * pi * event_mode**2)).sqrt()


def fresnel(z: arb, pi: arb) -> acb:
    i = acb(0, 1)
    return (i * pi / 4).exp() / arb(2).sqrt() * ((-i * pi / 4).exp() * (pi / 2).sqrt() * z).erf()


def main() -> None:
    priority = set_low_priority()
    require(RESULT.is_file() and NOTE.is_file(), "missing Fresnel-retention artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "Fresnel-retention artifact is not passed")
    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file() and file_hash(path) == record["sha256"], f"hash drift: {path}")

    ctx.dps = PRECISION
    pi = arb.pi()
    top = pi * C**2 / 8
    lower = pi * (C - 2) ** 2 / 8
    beta = top ** (arb(1) / 3)

    isolation_max = arb(0)
    isolation_index = -1
    disjoint_max = arb(0)
    disjoint_index = -1
    for index in range(COUNT):
        m = mode(index)
        tau = event(index, pi)
        signed = 4 * m - C
        require(signed == (-1 if index % 2 == 0 else 1) * (2 * index + 1), f"event mode law failed at {index}")
        d = beta * signed / C
        lam = (top - tau) / beta
        require((lam - d**2).contains(0), f"lambda=d^2 failed at {index}")
        Delta = (2 * pi * m**2 - tau) / (2 * pi * m**2 + tau)
        require((Delta - arb(signed) / C).contains(0), f"logistic characteristic failed at {index}")

        if index == 0:
            isolation_radius = pi / 16
            disjoint_radius = pi / 16
        else:
            down = pi * (index + 1) if index + 1 < COUNT else tau - lower
            isolation_radius = min(pi * index, down)
            disjoint_radius = pi * index / 2
        for family, radius in (("isolation", isolation_radius), ("disjoint", disjoint_radius)):
            local = max(
                abs(z_coordinate(tau - radius, tau, m, pi)),
                abs(z_coordinate(tau + radius, tau, m, pi)),
            )
            if family == "isolation" and local > isolation_max:
                isolation_max, isolation_index = local, index
            if family == "disjoint" and local > disjoint_max:
                disjoint_max, disjoint_index = local, index

    saved_isolation = arb(artifact["atlas_certificate"]["isolation_envelope"]["maximum_abs_Fresnel_coordinate_ball"])
    saved_disjoint = arb(artifact["atlas_certificate"]["pairwise_disjoint_symmetric_buffers"]["maximum_abs_Fresnel_coordinate_ball"])
    require(isolation_max.overlaps(saved_isolation) and isolation_index == 396, "independent isolation maximum drift")
    require(disjoint_max.overlaps(saved_disjoint) and disjoint_index == 398, "independent disjoint maximum drift")
    require(isolation_max < arb("0.004976") and disjoint_max < arb("0.002501"), "atlas Fresnel envelope target failed")
    require(arb(1) / 2 - isolation_max / arb(2).sqrt() > arb("0.496"), "atlas full-Morse defect lower bound failed")

    tau = event(0, pi)
    radius = pi / 16
    lower_z = z_coordinate(tau - radius, tau, Q, pi)
    upper_z = z_coordinate(tau + radius, tau, Q, pi)
    require(lower_z > 0 and upper_z < 0, "prototype orientation drift")
    saved_lower_z = arb(artifact["prototype_event"]["lower_face_Fresnel_coordinate_ball"])
    saved_upper_z = arb(artifact["prototype_event"]["upper_face_Fresnel_coordinate_ball"])
    require(lower_z.overlaps(saved_lower_z) and upper_z.overlaps(saved_upper_z), "prototype coordinate misses saved enclosure")

    half = acb(arb(1) / 2, arb(1) / 2)
    full = acb(1, 1)
    for z in (lower_z, upper_z):
        transition = half - fresnel(z, pi)
        relative_defect = abs(transition - full) / arb(2).sqrt()
        require(relative_defect > arb("0.4999994"), "direct prototype full-Morse defect bound failed")
        require(abs(transition - half) <= abs(z), "Fresnel unit-speed bound failed")

    decisions = artifact["decision"]
    require(decisions["exact_Airy_event_Fresnel_coordinate_join_proved"] is True, "coordinate join decision missing")
    require(decisions["naive_full_interior_Fresnel_face_replacement_admissible"] is False, "inadmissibility decision corrupted")
    require(decisions["exact_Fresnel_retention_required"] is True, "Fresnel-retention decision missing")
    require(decisions["Airy_to_Fresnel_retaining_logistic_join_proved"] is False, "replacement join overpromoted")
    require(decisions["propagation_across_399_events_proved"] is False, "atlas propagation overpromoted")
    require(decisions["rh_implication"] is False, "RH boundary corrupted")

    note = NOTE.read_text(encoding="utf-8")
    for token in (
        "lambda(t)-d_m^2=(tau_m-t)/beta",
        "full-interior Fresnel replacement",
        "universal logistic Morse phase",
        "does not prove the Airy-to-Fresnel-retaining",
    ):
        require(token in note, f"note token missing: {token}")
    print(
        "validated turning-event Fresnel retention independently: "
        f"events={COUNT}, isolation<0.004976, full-Morse-defect>0.496; priority={priority}"
    )


if __name__ == "__main__":
    main()
