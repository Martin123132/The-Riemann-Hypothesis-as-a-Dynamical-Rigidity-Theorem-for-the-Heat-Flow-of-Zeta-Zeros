#!/usr/bin/env python3
"""Independently check the half-domain Kummer normalization repair."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

from mpmath import mp

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_domain_Kummer_normalization_repair_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_reflection_check() -> tuple[str, str, str]:
    mp.dps = 95
    t = mp.mpf("9.5")
    alpha = mp.mpf(7)
    a = mp.mpf(3) / 4 - mp.j * t / 2
    b = mp.mpf(3) / 4 + mp.j * t / 2
    z = mp.j * mp.pi * alpha**2 / 4
    integrand = lambda x: alpha * x ** (a - 1) * (1 - x) ** (b - 1) * mp.exp(z * x)
    half = mp.quad(integrand, [0, mp.mpf("0.1"), mp.mpf("0.3"), mp.mpf("0.5")])
    full = mp.quad(integrand, [0, mp.mpf("0.1"), mp.mpf("0.3"), mp.mpf("0.5"), mp.mpf("0.8"), 1])
    closed = alpha * mp.beta(a, b) * mp.hyp1f1(a, mp.mpf("1.5"), z)
    rotation = mp.exp(-mp.j * mp.pi / 8)
    c_t = (mp.pi / (32 * t)) ** mp.mpf("0.25")
    euler = abs(full - closed) / abs(closed)
    reflection = abs(rotation * full - 2 * mp.re(rotation * half)) / abs(full)
    projection = abs(c_t * mp.re(rotation * closed) - 2 * c_t * mp.re(rotation * half))
    projection /= abs(c_t * mp.re(rotation * closed))
    require(euler < mp.mpf("1e-65"), "independent Euler-Kummer check failed")
    require(reflection < mp.mpf("1e-65"), "independent half-reflection check failed")
    require(projection < mp.mpf("1e-65"), "independent projection check failed")
    return mp.nstr(euler, 35), mp.nstr(reflection, 35), mp.nstr(projection, 35)


def main() -> int:
    require(gate.RESULT.is_file(), "missing normalization-repair result")
    require(gate.NOTE.is_file(), "missing normalization-repair note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == gate.STEM, "kind drift")
    require(artifact["status"] == "exact_half_domain_physical_projection_and_full_Kummer_label_normalization_repaired", "status drift")
    require(artifact.get("passed") is True, "normalization repair failed")
    require(artifact["scope"]["physical_x_interval"] == ["0", "1/2"], "half-domain scope drift")
    corrected = artifact["corrected_certificate"]
    require("integral_0^(1/2)" in corrected["physical_functional"], "corrected half-domain formula missing")
    require(corrected["one_source_label"].startswith("K_t(alpha)=c_t "), "corrected Kummer prefactor drift")
    require("2c_t" not in corrected["one_source_label"], "factor two reintroduced")
    decision = artifact["decision"]
    require(decision["literal_full_domain_2Re_functional_rejected"] is True, "literal full-domain defect not rejected")
    require(decision["half_domain_physical_functional_certified"] is True, "half-domain functional not certified")
    require(decision["full_Kummer_label_prefactor_corrected_from_2c_to_c"] is True, "Kummer factor not repaired")
    require(decision["scalar_corridor_arithmetic_changed"] is False, "scalar corridor was spuriously rescaled")
    require(decision["Q_K_numerically_enclosed"] is False, "Q_K overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")

    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file {path}")
            require(file_hash(path) == record["sha256"], f"hash drift {path}")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("doubles every full Kummer label", "corrected physical functional", "does not rescale", "no numerical enclosure"):
        require(token in note, f"note boundary token missing: {token}")
    euler, reflection, projection = independent_reflection_check()
    print(
        "independently checked half-domain normalization repair; "
        f"Euler={euler}, reflection={reflection}, projection={projection}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
