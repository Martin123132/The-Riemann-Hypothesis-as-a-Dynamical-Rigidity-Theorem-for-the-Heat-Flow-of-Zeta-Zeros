#!/usr/bin/env python3
"""Validate the logarithmic-phase jet and leading Turan reduction."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = REPO_ROOT / "work/rh_compute/scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_logarithmic_phase_turan_reduction as gate  # noqa: E402


STEM = gate.STEM
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    f"ltp_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "domain",
            "polar",
            "fourier",
            "jets",
            "operator",
            "phases",
            "bochner",
            "boundary",
            "trace",
            "observations",
            "turan",
            "zero",
            "guard",
            "correction",
            "target",
        ),
        start=1,
    )
]
EXPECTED_COUNTS = {
    "rows": 15,
    "positive_amplitude_representations": 1,
    "fourier_jet_identities": 5,
    "differential_operator_identities": 1,
    "phase_difference_families": 1,
    "phase_sum_families": 1,
    "bochner_psd_kernels": 1,
    "leading_turan_identities": 1,
    "leading_zero_fibre_signs": 1,
    "integer_log_phase_guards": 1,
    "exact_correction_decompositions": 1,
    "signed_joint_targets": 1,
    "signed_joint_bounds": 0,
    "phi_b_bounds": 0,
    "xi_level_current_theorems": 0,
}
SUMMARY = (
    "validated logarithmic-phase Turan reduction: 15 rows, 5 Fourier "
    "jet identities, 1 Bochner kernel, 1 leading Turan identity, "
    "1 zero-fibre sign, 1 integer-log guard, 0 signed joint bounds"
)


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path, issues: list[str]) -> dict:
    if not path.is_file():
        issues.append(f"missing result: {path}")
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"invalid result: {exc}")
        return {}


def require_zero(expression: sp.Expr, label: str, issues: list[str]) -> None:
    if sp.simplify(sp.expand(expression)) != 0:
        issues.append(label)


def independent_four_carrier_jet_audit(issues: list[str]) -> None:
    i = sp.I
    z = sp.symbols("z", real=True)
    tau = sp.Rational(7, 11)
    omega = (sp.Integer(3) + 4 * i) / 5
    if sp.simplify(omega * sp.conjugate(omega) - 1) != 0:
        issues.append("independent unit phase failed")
    amplitudes = (
        sp.Rational(1, 2),
        sp.Rational(2, 3),
        sp.Rational(5, 7),
        sp.Rational(11, 13),
    )
    frequencies = (sp.Integer(0), sp.log(2), sp.log(3), sp.log(5))
    fourier = sum(
        amplitude * sp.exp(-i * z * frequency)
        for amplitude, frequency in zip(amplitudes, frequencies)
    )
    shifted = omega * fourier
    weights = tuple(
        omega * amplitude * sp.exp(-i * tau * frequency)
        for amplitude, frequency in zip(amplitudes, frequencies)
    )
    for order in range(5):
        moment = sum(
            frequency**order * weight
            for frequency, weight in zip(frequencies, weights)
        )
        jet = i**order * sp.diff(shifted, z, order).subs(z, tau)
        require_zero(
            jet - moment,
            f"independent order-{order} Fourier jet failed",
            issues,
        )

    coefficients = (
        sp.Rational(2, 5),
        -sp.Rational(3, 7),
        sp.Rational(5, 11),
        -sp.Rational(7, 13),
        sp.Rational(11, 17),
    )
    direct = sum(
        sum(coefficients[r] * frequency**r for r in range(5)) * weight
        for frequency, weight in zip(frequencies, weights)
    )
    differential = sum(
        coefficients[r] * i**r * sp.diff(shifted, z, r).subs(z, tau)
        for r in range(5)
    )
    require_zero(
        differential - direct,
        "independent polynomial differential functional failed",
        issues,
    )

    for n in range(4):
        for m in range(4):
            ratio_expected = (
                amplitudes[n]
                * amplitudes[m]
                * sp.exp(-i * tau * (frequencies[n] - frequencies[m]))
            )
            product_expected = (
                omega**2
                * amplitudes[n]
                * amplitudes[m]
                * sp.exp(-i * tau * (frequencies[n] + frequencies[m]))
            )
            require_zero(
                weights[n] * sp.conjugate(weights[m]) - ratio_expected,
                f"independent ratio phase failed at {n},{m}",
                issues,
            )
            require_zero(
                weights[n] * weights[m] - product_expected,
                f"independent product phase failed at {n},{m}",
                issues,
            )


def independent_bochner_audit(issues: list[str]) -> None:
    i = sp.I
    points = (
        sp.Integer(0),
        sp.Rational(1, 3),
        -sp.Rational(2, 5),
        sp.Rational(4, 7),
    )
    coefficients = (
        1 + i / 2,
        -sp.Rational(2, 3) + i / 5,
        sp.Rational(3, 7) - 2 * i / 9,
        -sp.Rational(5, 11) - i / 13,
    )
    amplitudes = (
        sp.Rational(1, 2),
        sp.Rational(3, 5),
        sp.Rational(7, 11),
        sp.Rational(13, 17),
    )
    frequencies = (sp.Integer(0), sp.log(2), sp.log(3), sp.log(7))

    def fourier(x: sp.Expr) -> sp.Expr:
        return sum(
            amplitude * sp.exp(-i * x * frequency)
            for amplitude, frequency in zip(amplitudes, frequencies)
        )

    direct = sum(
        coefficients[j]
        * sp.conjugate(coefficients[k])
        * fourier(points[j] - points[k])
        for j in range(4)
        for k in range(4)
    )
    factored = sum(
        amplitude
        * sum(
            coefficients[j] * sp.exp(-i * points[j] * frequency)
            for j in range(4)
        )
        * sp.conjugate(
            sum(
                coefficients[k] * sp.exp(-i * points[k] * frequency)
                for k in range(4)
            )
        )
        for amplitude, frequency in zip(amplitudes, frequencies)
    )
    require_zero(
        direct - factored,
        "independent Bochner Gram factorization failed",
        issues,
    )
    numerical = complex(sp.N(factored, 40))
    if abs(numerical.imag) > 1e-30 or numerical.real < -1e-30:
        issues.append("independent Bochner numerical PSD check failed")


def independent_leading_current_audit(issues: list[str]) -> None:
    i = sp.I
    y = sp.symbols("y", real=True)
    distances = (
        sp.Rational(1, 5),
        sp.Rational(2, 7),
        sp.Rational(3, 11),
        sp.Rational(5, 13),
    )
    weights = (
        1 + i / 3,
        -sp.Rational(2, 5) + i / 7,
        sp.Rational(4, 9) - 2 * i / 11,
        -sp.Rational(3, 8) - i / 13,
    )
    u_x = sp.Rational(2, 17)
    moments = tuple(
        sum(distance**order * weight for distance, weight in zip(distances, weights))
        for order in range(3)
    )
    complex_trace = sum(
        weight * sp.exp(-i * y * distance)
        for distance, weight in zip(distances, weights)
    )
    f = sp.expand_complex(sp.re(complex_trace))
    g = sp.expand_complex(sp.im(complex_trace))
    f0 = sp.simplify(f.subs(y, 0))
    f1 = sp.simplify(sp.diff(f, y).subs(y, 0))
    f2 = sp.simplify(sp.diff(f, y, 2).subs(y, 0))
    g0 = sp.simplify(g.subs(y, 0))
    require_zero(f0 - sp.re(moments[0]), "independent f(0) failed", issues)
    require_zero(f1 - sp.im(moments[1]), "independent f'(0) failed", issues)
    require_zero(f2 + sp.re(moments[2]), "independent f''(0) failed", issues)
    require_zero(g0 - sp.im(moments[0]), "independent g(0) failed", issues)

    direct_current = sp.expand(
        f0 * (f2 / 4 + u_x * g0 / 2) - (f1 / 2) ** 2
    )
    turan_current = sp.expand(
        (f0 * f2 - f1**2) / 4 + u_x * f0 * g0 / 2
    )
    require_zero(
        direct_current - turan_current,
        "independent leading Turan observation failed",
        issues,
    )

    pair_sum = 0
    for un, wn in zip(distances, weights):
        for um, wm in zip(distances, weights):
            h0 = -(un + um) ** 2 / 8
            t0 = -(un - um) ** 2 / 8 - i * u_x / 2
            pair_sum += h0 * wn * sp.conjugate(wm) + t0 * wn * wm
    pair_current = sp.expand_complex(sp.re(pair_sum) / 2)
    require_zero(
        pair_current - turan_current,
        "independent leading pair-current reduction failed",
        issues,
    )

    f_symbol, fp_symbol, fpp_symbol, g_symbol, ux_symbol = sp.symbols(
        "f fp fpp g ux", real=True
    )
    generic = (
        f_symbol * fpp_symbol - fp_symbol**2
    ) / 4 + ux_symbol * f_symbol * g_symbol / 2
    zero_fibre = sp.simplify(generic.subs(f_symbol, 0))
    if zero_fibre != -fp_symbol**2 / 4:
        issues.append("independent zero-fibre sign failed")


def independent_integer_log_guard(issues: list[str]) -> None:
    d = sp.log(2)
    tau = sp.pi / d
    omega = -sp.Integer(1)
    w1 = sp.simplify(omega * sp.Rational(1, 2))
    w2 = sp.simplify(omega * sp.exp(-sp.I * tau * d))
    if (w1, w2) != (-sp.Rational(1, 2), sp.Integer(1)):
        issues.append("independent integer-log phase weights failed")
        return
    s0 = w1 + w2
    s1 = 2 * d * w1 + d * w2
    s2 = (2 * d) ** 2 * w1 + d**2 * w2
    f0 = sp.re(s0)
    f1 = sp.im(s1)
    f2 = -sp.re(s2)
    g0 = sp.im(s0)
    ux = sp.symbols("u_x", real=True)
    current = sp.simplify(
        (f0 * f2 - f1**2) / 4 + ux * f0 * g0 / 2
    )
    if current != d**2 / 8:
        issues.append("independent integer-log positive-current guard failed")


def independent_source_audit(payload: dict, issues: list[str]) -> None:
    hashes = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in gate.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
        elif hashes.get(key) != file_hash(path):
            issues.append(f"source hash drifted: {key}")
    try:
        anchor = json.loads(
            gate.SOURCE_PATHS["phase_anchor"].read_text(encoding="utf-8")
        )
        mangoldt = json.loads(
            gate.SOURCE_PATHS["mangoldt"].read_text(encoding="utf-8")
        )
        two_carrier = json.loads(
            gate.SOURCE_PATHS["two_carrier"].read_text(encoding="utf-8")
        )
        growing = json.loads(
            gate.SOURCE_PATHS["growing_tail"].read_text(encoding="utf-8")
        )
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"independent source load failed: {exc}")
        return
    if "eta=f_1/|f_1|" not in anchor.get("exact", {}).get("unit_anchor", ""):
        issues.append("independent unit-anchor provenance failed")
    correction_free = mangoldt.get("symbolic_certificate", {}).get(
        "correction_free_moments", ""
    )
    if "0<=r<=4" not in correction_free:
        issues.append("independent five-moment provenance failed")
    if not two_carrier.get("symbolic_certificate", {}).get(
        "full_two_carrier_identity"
    ):
        issues.append("independent two-carrier provenance failed")
    if "u_x=h^2/(8*pi)" not in growing.get("exact", {}).get(
        "source_bounds", ""
    ):
        issues.append("independent physical-pi provenance failed")


def validate(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("unexpected artifact kind")
    if payload.get("date") != "2026-08-01":
        issues.append("artifact date drifted")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")
    rows = payload.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or order drifted")
    if sum(row.get("readiness") == "open" for row in rows) != 1:
        issues.append("open-row count drifted")

    status = payload.get("status", "")
    for token in (
        "positive-amplitude logarithmic Fourier-jet",
        "leading Turan-current",
        "signed joint bound open",
        "no Phi_B bound",
        "or RH",
    ):
        if token not in status:
            issues.append(f"status missing token: {token}")
    boundary = payload.get("proof_boundary", "")
    for token in (
        "no signed joint bound",
        "upper bound on Phi_B",
        "contact exclusion",
        "Lambda<=0",
        "RH",
    ):
        if token not in boundary:
            issues.append(f"proof boundary missing token: {token}")

    independent_source_audit(payload, issues)
    try:
        rebuilt = gate.build_payload()
        if payload != rebuilt:
            issues.append("full builder recomputation drifted")
    except Exception as exc:  # noqa: BLE001
        issues.append(f"builder recomputation failed: {exc}")

    counts = payload.get("counts", {})
    if counts.get("signed_joint_bounds") != 0:
        issues.append("unexpected signed joint bound")
    if counts.get("phi_b_bounds") != 0:
        issues.append("unexpected Phi_B bound")
    if counts.get("xi_level_current_theorems") != 0:
        issues.append("unexpected Xi-level theorem")

    independent_four_carrier_jet_audit(issues)
    independent_bochner_audit(issues)
    independent_leading_current_audit(issues)
    independent_integer_log_guard(issues)


def validate_note(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for token in (
        "Logarithmic-Phase Jet And Turan Reduction",
        "Date: 2026-08-01",
        "0 signed joint bounds",
        "G_r=i^r*Z^(r)(tau)",
        "Bochner factorization",
        "P_bulk^(0)=V*N-A*Q",
        "P_bulk^(0)=-(f'(0))^2/4<=0",
        "P_bulk^(0)=log(2)^2/8>0",
        "Pi Provenance",
        "h^2*Phi_B=E+P_bulk^(0)+R_corr",
        "RH, or prize-level conclusion",
    ):
        if token not in text:
            issues.append(f"note marker missing: {token}")


def main() -> int:
    issues: list[str] = []
    payload = load_json(RESULT, issues)
    if payload:
        validate(payload, issues)
    validate_note(issues)
    if issues:
        for issue in issues:
            print(f"ISSUE {issue}")
        return 1
    print(SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
