#!/usr/bin/env python3
"""Independently check the finite Fresnel-cell common-carrier reduction."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_cell_common_carrier_subtraction_reduction_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finite_source(alpha_min: int, count: int, y: mp.mpf | mp.mpc) -> mp.mpc:
    return mp.fsum(mp.exp(1j * mp.pi * (alpha_min + 2 * j) * y) for j in range(count))


def main() -> int:
    require(gate.RESULT.is_file(), "missing common-carrier result")
    require(gate.NOTE.is_file(), "missing common-carrier note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == gate.STEM, "kind drift")
    require(
        artifact["status"]
        == "exact_Fresnel_cell_common_Gamma_and_A_carrier_subtraction_reduction_certified_no_amplitude_bound",
        "status drift",
    )
    require(artifact.get("passed") is True, "common-carrier gate failed")

    mp.mp.dps = 100
    t = mp.mpf("3.75")
    s = mp.mpf("0.5") + 1j * t
    alpha_min = 3
    count = 4
    for mode, u in ((1, mp.mpf("-0.27")), (3, mp.mpf("0.19")), (6, mp.mpf("0.37"))):
        direct = mp.exp(-1j * mp.pi * (mode + u) ** 2) * finite_source(alpha_min, count, mode + u)
        reduced = mp.exp(-1j * mp.pi * u * u - 2j * mp.pi * mode * u) * finite_source(alpha_min, count, u)
        require(abs(direct - reduced) < mp.mpf("1e-90"), "fresh chirp/source covariance failed")

        direct_cell = mp.quad(
            lambda y: y ** (-s) * mp.exp(-1j * mp.pi * y * y) * finite_source(alpha_min, count, y),
            [mp.mpf(mode) - mp.mpf("0.5"), mode, mp.mpf(mode) + mp.mpf("0.5")],
        )
        beta = mp.quad(
            lambda v: (1 + v / mode) ** (-s)
            * mp.exp(-1j * mp.pi * v * v - 2j * mp.pi * mode * v)
            * finite_source(alpha_min, count, v),
            [mp.mpf("-0.5"), 0, mp.mpf("0.5")],
        )
        require(abs(direct_cell - mode ** (-s) * beta) < mp.mpf("1e-85"), "fresh cell factorization failed")

    theta = mp.im(mp.loggamma(mp.mpf("0.25") + 0.5j * t)) - t * mp.log(mp.pi) / 2
    theta_zero = t * (mp.log(t / (2 * mp.pi)) - 1) / 2 - mp.pi / 8
    exponent = mp.mpf("0.75") + 0.5j * t
    bulk = (mp.exp(0.5j * t) * mp.gamma(exponent) / (0.5j * t) ** exponent) / (
        2 * mp.sqrt(mp.pi / t) * mp.exp(-1j * mp.pi / 4)
    )
    h_factor = (
        2 ** mp.mpf("1.25")
        * t ** mp.mpf("0.25")
        * mp.pi ** mp.mpf("1.5")
        * mp.exp(-3 * mp.pi * t / 4)
        / (
            (1 + mp.exp(-2 * mp.pi * t))
            * abs(mp.gamma(mp.mpf("0.25") + 0.5j * t))
            * abs(mp.gamma(mp.mpf("0.75") + 0.5j * t)) ** 2
        )
    )
    common = h_factor * mp.exp(-1j * theta) * bulk * mp.exp(1j * theta_zero)
    for mode in (3, 8, 13):
        lifted = 2 * mp.re(mp.exp(1j * theta) * common * mode ** (-s))
        physical = h_factor * 2 * mp.re(bulk * mp.exp(1j * theta_zero) * mode ** (-s))
        require(abs(lifted - physical) < mp.mpf("1e-90"), "fresh Gamma Hardy lift failed")

    cell = artifact["exact_cell_reduction"]
    require(cell["mode_cell"] == "B_m=m^(-s)beta_m", "cell factor drift")
    require("no infinite sum-integral interchange" in cell["construction"], "interchange guard lost")
    gamma = artifact["Gamma_common_carrier"]
    require(gamma["ordinary_cell_defect"] == "m^(-s)[beta_m-C_G]", "Gamma defect drift")
    a_lift = artifact["A_common_lift"]
    require(a_lift["natural_same_cell_integral_lift_derived"] is False, "A lift overclaim")
    rosters = artifact["mode_rosters"]
    require(rosters["ordinary_owned_cells"] == [622, 39_852], "ordinary roster drift")
    require(rosters["A_transition_cells"] == [39_853, 39_936], "A roster drift")
    require(rosters["A_transition_half_counts"] == [42, 42], "A half-count drift")
    require(rosters["symbolic_HJZ_reassembly_zero"] == "0", "H J_Z algebra drift")
    decision = artifact["decision"]
    require(decision["Gamma_carrier_subtracted_on_common_m_minus_s_coefficient"] is True, "Gamma subtraction lost")
    require(decision["P_W_joined_to_unowned_cell_packages"] is False, "P_W join overclaim")
    require(decision["owned_cell_amplitude_defects_bounded"] is False, "amplitude bound overclaim")
    require(decision["actual_height_J_Z_enclosed"] is False, "J_Z overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")

    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file {path}")
            require(file_hash(path) == record["sha256"], f"hash drift {path}")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in (
        "subtraction before norms",
        "39231 ordinary Gamma-subtracted cells",
        "natural integral lift",
        "Mordell/contour transformation",
        "No fitted constant",
        "No bound for any coefficient defect",
    ):
        require(token in note, f"note boundary token missing: {token}")

    print("independently checked Fresnel-cell common-carrier subtraction reduction", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
