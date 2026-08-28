#!/usr/bin/env python3
"""Independently check the branch/Fresnel connection and recovered mode roster."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_branch_Fresnel_Dirichlet_connection_and_mode_roster_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def half_lines(t: mp.mpf, alpha: int) -> tuple[mp.mpc, mp.mpc, mp.mpc, mp.mpc]:
    s = mp.mpf("0.5") + 1j * t
    c = mp.pi * alpha * (1 + 1j) / mp.sqrt(2)
    cuts = [0, mp.mpf("0.02"), mp.mpf("0.2"), 1, 3, 7, mp.inf]
    iminus = mp.quad(lambda x: x ** (-s) * mp.e ** (-mp.pi * x * x - c * x), cuts)
    iplus = mp.quad(lambda x: x ** (-s) * mp.e ** (-mp.pi * x * x + c * x), cuts)
    j = mp.quad(lambda x: x ** (-mp.conj(s)) * mp.e ** (-mp.pi * x * x + 1j * c * x), cuts)
    kappa = 1j * mp.sqrt(2 * mp.pi) * (2 * mp.pi) ** (1j * t) * mp.e ** (-mp.pi * t / 2) / mp.gamma(mp.mpf("0.5") + 1j * t)
    return iminus, iplus, j, kappa


def rotated_plus(t: mp.mpf, alpha: int, delta: mp.mpf) -> tuple[mp.mpc, mp.mpc]:
    s = mp.mpf("0.5") + 1j * t
    c = mp.pi * alpha * (1 + 1j) / mp.sqrt(2)
    phi = mp.pi / 4 - delta
    cuts = [0, mp.mpf("0.02"), mp.mpf("0.2"), 1, 3, 7, mp.inf]
    direct = mp.quad(lambda x: x ** (-s) * mp.e ** (-mp.pi * x * x + c * x), cuts)
    ray = mp.e ** (1j * phi * (1 - s)) * mp.quad(
        lambda y: y ** (-s) * mp.e ** (-mp.pi * mp.e ** (2j * phi) * y * y + c * mp.e ** (1j * phi) * y),
        cuts,
    )
    return direct, ray


def main() -> int:
    require(gate.RESULT.is_file(), "missing Fresnel-connection result")
    require(gate.NOTE.is_file(), "missing Fresnel-connection note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == gate.STEM, "kind drift")
    require(
        artifact["status"]
        == "exact_parabolic_cylinder_branch_connection_Fresnel_Dirichlet_sum_and_mode_roster_certified",
        "status drift",
    )
    require(artifact.get("passed") is True, "Fresnel-connection gate failed")

    mp.mp.dps = 110
    for t, alpha in ((mp.mpf("2.75"), 1), (mp.mpf("4.25"), 3)):
        iminus, iplus, j, kappa = half_lines(t, alpha)
        require(abs(iminus + 1j * mp.e ** (-mp.pi * t) * iplus - kappa * j) < mp.mpf("1e-40"), "fresh DLMF connection failed")
        require(abs(abs(kappa) - mp.sqrt(1 + mp.e ** (-2 * mp.pi * t))) < mp.mpf("1e-65"), "fresh kappa modulus failed")

    direct, ray = rotated_plus(mp.mpf("3.25"), 3, mp.mpf("0.16"))
    require(abs(direct - ray) < mp.mpf("1e-55"), "fresh in-sector contour rotation failed")

    t = mp.mpf(gate.HEIGHT)
    transition = mp.sqrt(8 * t / mp.pi)
    center = mp.sqrt(t / (2 * mp.pi))
    require(159575 < transition < 159577, "fresh alpha transition bracket failed")
    require(39894 < center < 39895, "fresh center bracket failed")

    def roots(alpha: int) -> tuple[mp.mpf, mp.mpf]:
        radical = mp.sqrt(alpha * alpha - 8 * t / mp.pi)
        return (alpha - radical) / 4, (alpha + radical) / 4

    lower_A, upper_A = roots(gate.SOURCE_ALPHA_MIN)
    lower_B, upper_B = roots(gate.SOURCE_ALPHA_MAX)
    require(621 < lower_B < 622, "fresh B lower saddle bracket failed")
    require(39852 < lower_A < 39853, "fresh A lower saddle bracket failed")
    require(39936 < upper_A < 39937, "fresh A upper saddle bracket failed")
    require(2560588 < upper_B < 2560589, "fresh B upper saddle bracket failed")

    connection = artifact["parabolic_cylinder_connection"]
    require(connection["source_split"].endswith("(S_alpha-T_alpha)+T_alpha"), "source split drift")
    fresnel = artifact["finite_Fresnel_branch_sum"]
    require(fresnel["integer_singularities_removable_by_finite_sum"] is True, "Dirichlet removability lost")
    require(fresnel["infinite_A21_interchange_used"] is False, "A21 reintroduced")
    roster = artifact["saddle_mode_roster"]
    require(roster["ordinary_integer_roster"] == [622, 39_852], "ordinary roster drift")
    require(roster["ordinary_integer_count"] == 39_231, "ordinary count drift")
    require(roster["transition_integer_roster"] == [39_853, 39_936], "transition roster drift")
    require(roster["transition_half_counts"] == [42, 42], "transition half counts drift")
    ownership = artifact["physical_ownership_match"]
    require(ownership["independent_saddle_roster_matches_existing_indices"] is True, "ownership match lost")
    require(ownership["carrier_amplitude_subtraction_on_common_Fresnel_integral_complete"] is False, "carrier subtraction overclaim")
    decision = artifact["decision"]
    require(decision["branch_Stokes_connection_complete"] is True, "Stokes connection lost")
    require(decision["finite_branch_sum_Fresnel_Dirichlet_representation_exact"] is True, "Fresnel sum lost")
    require(decision["common_integral_G_plus_A_transition_subtraction_complete"] is False, "joined subtraction overclaim")
    require(decision["actual_height_joined_remainder_enclosed"] is False, "joined remainder overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")

    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file {path}")
            require(file_hash(path) == record["sha256"], f"hash drift {path}")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in (
        "exactly the finite-label split",
        "Fresnel/Abel limit",
        "39231 modes",
        "42 modes on each side",
        "not a new amplitude estimate",
        "No fitted constant",
    ):
        require(token in note, f"note boundary token missing: {token}")

    print("independently checked branch/Stokes connection, Fresnel sum, and mode roster", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
