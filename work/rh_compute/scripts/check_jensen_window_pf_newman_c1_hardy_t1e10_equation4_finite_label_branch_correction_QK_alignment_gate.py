#!/usr/bin/env python3
"""Independently check the finite-label branch correction and Q_K alignment."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_label_branch_correction_QK_alignment_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finite_label_formula(t: mp.mpf, alpha: int) -> tuple[mp.mpc, mp.mpc, mp.mpc, mp.mpc, mp.mpc]:
    pi = mp.pi
    s = mp.mpf("0.5") + 1j * t
    a = mp.mpf("0.25") - 0.5j * t
    b = mp.mpf("0.75") - 0.5j * t
    xi = 1j * pi * alpha * alpha / 4
    c = pi * alpha * (1 + 1j) / mp.sqrt(2)
    d = mp.e ** (-pi * t / 4 - 1j * pi / 8)
    A = mp.gamma(a) * pi ** (-a) * mp.hyp1f1(a, mp.mpf("0.5"), xi) / 2
    B = c * mp.gamma(b) * pi ** (-b) * mp.hyp1f1(b, mp.mpf("1.5"), xi) / 2
    source = 1j * mp.e ** (pi * t) * d * (A - B)
    correction = d * (A + B)
    exact = source - correction
    return s, exact, source, correction, A + B


def direct_contour(t: mp.mpf, alpha: int) -> mp.mpc:
    pi = mp.pi
    s = mp.mpf("0.5") + 1j * t
    rotation = (1 - 1j) / mp.sqrt(2)

    def integrand(q: mp.mpf) -> mp.mpc:
        z = mp.mpf("0.5") + q * rotation
        return -mp.e ** (-1j * pi * z * z + 1j * pi * alpha * z) * mp.power(z, -s) * rotation

    saddle = mp.mpf(alpha) / (2 * mp.sqrt(2))
    return mp.quad(integrand, [-mp.inf, -5, 0, saddle - 2, saddle, saddle + 2, 8, mp.inf])


def hardy_phase(t: mp.mpf) -> mp.mpf:
    return mp.im(mp.loggamma(mp.mpf("0.25") + 0.5j * t)) - t * mp.log(mp.pi) / 2


def main() -> int:
    require(gate.RESULT.is_file(), "missing branch-correction result")
    require(gate.NOTE.is_file(), "missing branch-correction note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == gate.STEM, "kind drift")
    require(
        artifact["status"]
        == "exact_finite_label_branch_correction_and_augmented_QK_contour_alignment_certified",
        "status drift",
    )
    require(artifact.get("passed") is True, "branch-correction gate failed")

    mp.mp.dps = 75
    for t, alpha in ((mp.mpf("3.25"), 1), (mp.mpf("3.25"), 3), (mp.mpf("4.5"), 5)):
        _, exact, source, correction, half_line = finite_label_formula(t, alpha)
        quadrature = direct_contour(t, alpha)
        require(abs(exact - quadrature) < mp.mpf("1e-62"), "fresh finite-label contour check failed")
        require(abs(source - exact - correction) < mp.mpf("1e-68"), "fresh branch identity failed")

        theta = hardy_phase(t)
        projection_exact = 2 * mp.re(mp.e ** (1j * theta) * exact)
        a = mp.mpf("0.25") - 0.5j * t
        b = mp.mpf("0.75") - 0.5j * t
        xi = 1j * mp.pi * alpha * alpha / 4
        phi_one = mp.hyp1f1(a, mp.mpf("0.5"), xi)
        phi_two = mp.hyp1f1(b, mp.mpf("1.5"), xi)
        gamma_quarter = abs(mp.gamma(mp.mpf("0.25") + 0.5j * t))
        projection_exact_closed = (
            -mp.e ** (-mp.pi * t / 4)
            * gamma_quarter
            * mp.pi ** (-mp.mpf("0.25"))
            * mp.re(mp.e ** (-1j * mp.pi / 8) * phi_one)
        )
        projection_source_closed = (
            2
            * mp.pi ** mp.mpf("1.25")
            * alpha
            * mp.e ** (-3 * mp.pi * t / 4)
            * mp.re(mp.e ** (-1j * mp.pi / 8) * phi_two)
            / ((1 + mp.e ** (-2 * mp.pi * t)) * gamma_quarter)
        )
        require(abs(projection_exact - projection_exact_closed) < mp.mpf("1e-65"), "fresh Phi1 projection failed")
        require(
            abs(2 * mp.re(mp.e ** (1j * theta) * source) - projection_source_closed) < mp.mpf("1e-65"),
            "fresh Phi2 projection failed",
        )
        require(abs(half_line - (source - exact) / (mp.e ** (-mp.pi * t / 4 - 1j * mp.pi / 8))) < mp.mpf("1e-68"), "half-line correction drift")

    decision = artifact["decision"]
    require(decision["finite_RSI_label_equals_A33_source_label"] is False, "false label equality promoted")
    require(decision["uncorrected_finite_residual_permitted_as_HQK"] is False, "uncorrected Q_K route promoted")
    require(decision["branch_correction_rigorously_nonzero_at_surrogate"] is True, "branch witness lost")
    require(decision["A21_infinite_interchange_required_for_augmented_finite_roster"] is False, "A21 reintroduced")
    require(decision["physical_QK_prefactor_and_branch_alignment_complete_with_correction"] is True, "corrected alignment lost")
    require(decision["actual_window_augmented_common_contour_identity_exact"] is True, "actual-window identity lost")
    require(decision["branch_correction_sum_enclosed"] is False, "branch sum overclaim")
    require(decision["D_K_enclosed"] is False, "D_K overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")

    window = artifact["actual_window_augmented_identity"]
    require(window["alpha_min"] == 159_577 and window["alpha_max"] == 5_122_421, "source endpoints drift")
    require(window["window_label_count"] == 2_481_423, "source count drift")
    require(window["released_Dirichlet_block"] == [79_789, 2_561_211], "Dirichlet block drift")
    require("+B_W" in window["physical_identity"], "branch sum omitted from physical identity")

    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file {path}")
            require(file_hash(path) == record["sha256"], f"hash drift {path}")
    for record in artifact["references"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file(), f"missing reference {path}")
        require(file_hash(path) == record["sha256"], f"reference hash drift {path}")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in (
        "not an unresolved sign convention",
        "T_alpha=S_alpha-d I_alpha(c_alpha)",
        "uncorrected identification is false",
        "uses no infinite A21 interchange",
        "No fitted or unexplained circle constant",
        "branch sum `B_W` must",
    ):
        require(token in note, f"note boundary token missing: {token}")

    print("independently checked finite-label branch correction and augmented Q_K alignment", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
