#!/usr/bin/env python3
"""Validate the six-moment endpoint-composition and retention gate."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = REPO_ROOT / "work/rh_compute/scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_endpoint_composition_retention_gate as gate  # noqa: E402


STEM = gate.STEM
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    f"mfcr_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "domain",
            "aggregates",
            "lower",
            "upper",
            "identity",
            "lower_expand",
            "upper_expand",
            "derivative",
            "real",
            "lower_floor",
            "upper_floor",
            "geometry",
            "cotangent",
            "geometric",
            "transfer",
            "invariance",
            "carrier",
            "normalize",
            "flow",
            "terminal",
            "tail",
            "interior",
            "pi",
            "boundary",
        ),
        start=1,
    )
]
EXPECTED_COUNTS = {
    "rows": 24,
    "endpoint_composed_value_identities": 6,
    "stable_endpoint_packages": 2,
    "roster_transfer_representations": 2,
    "retained_endpoint_floors": 7,
    "observation_expansions": 1,
    "enumerated_physical_modes": 0,
    "oscillatory_interior_bounds": 0,
    "endpoint_composed_h2_bounds": 0,
    "evaluated_signed_carriers": 0,
    "signed_flow_bounds": 0,
    "phi_b_bounds": 0,
}
SUMMARY = (
    "validated Morse-Fresnel endpoint-composition retention gate: 24 rows, "
    "6 composed identities, 2 roster representations, "
    "7 retained endpoint floors, 1 observation expansion, "
    "0 oscillatory interior bounds, 0 signed flow bounds"
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


def close(
    left: mp.mpf | mp.mpc,
    right: mp.mpf | mp.mpc,
    tolerance: str = "1e-45",
) -> bool:
    return abs(left - right) <= mp.mpf(tolerance) * max(1, abs(left), abs(right))


def independent_source_audit(payload: dict, issues: list[str]) -> None:
    saved = payload.get("source_audit", {})
    for key, path in gate.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
            continue
        relative = str(path.relative_to(REPO_ROOT)).replace("\\", "/")
        entry = saved.get(key, {})
        if entry.get("path") != relative:
            issues.append(f"source path drifted: {key}")
        if entry.get("sha256") != file_hash(path):
            issues.append(f"source hash drifted: {key}")

    expected = {
        "endpoint_coherence": {
            "endpoint_phase_collapses": 2,
            "endpoint_composed_bounds": 0,
        },
        "morse_fresnel_endpoint": {
            "exact_transition_integrals": 6,
            "stable_endpoint_sum_families": 3,
        },
        "flow_matrix": {
            "correction_free_moments": 6,
            "real_flow_observations": 8,
        },
    }
    for source, checks in expected.items():
        path = gate.SOURCE_PATHS[source]
        if not path.is_file():
            continue
        counts = json.loads(path.read_text(encoding="utf-8")).get("counts", {})
        for key, value in checks.items():
            if counts.get(key) != value:
                issues.append(f"source count drifted: {source}.{key}")


def independent_symbolic_audit(issues: list[str]) -> None:
    delta, e_b, a_b, kappa = sp.symbols("delta e_B A_B kappa")
    cal_1, cal_b, c_1, c_b = sp.symbols("cal_1 cal_B c_1 c_B")
    fresnel, interior, tail = sp.symbols("F J T")
    direct = (
        delta / 2
        + e_b * a_b / 2
        + fresnel
        + kappa * c_1
        - kappa * e_b * c_b
        + interior
        + e_b * cal_b
        - cal_1
        + tail
    )
    lower = delta / 2 - cal_1 + kappa * c_1
    upper = a_b / 2 + cal_b - kappa * c_b
    rebuilt = fresnel + lower + e_b * upper + interior + tail
    if sp.expand(direct - rebuilt) != 0:
        issues.append("independent endpoint collection failed")

    p0, p1, r0, r1 = sp.symbols("p0 p1 r0 r1", real=True)
    m00, m01, m11, l0, l1 = sp.symbols("m00 m01 m11 l0 l1", real=True)
    p = sp.Matrix([p0, p1])
    rho = sp.Matrix([r0, r1])
    matrix = sp.Matrix([[m00, m01], [m01, m11]])
    linear = sp.Matrix([l0, l1])
    full = ((p + rho).T * matrix * (p + rho))[0] + (linear.T * (p + rho))[0]
    split = (
        (p.T * matrix * p)[0]
        + (linear.T * p)[0]
        + 2 * (p.T * matrix * rho)[0]
        + (rho.T * matrix * rho)[0]
        + (linear.T * rho)[0]
    )
    if sp.expand(full - split) != 0:
        issues.append("independent observation expansion failed")


def independent_rational_audit(issues: list[str]) -> None:
    lower = Fraction(2603, 7200)
    if lower / 4096 >= Fraction(1, 10_000):
        issues.append("lower endpoint floor reserve failed")
    upper = Fraction(15_653, 1800)
    if upper / 4096 >= Fraction(1, 100):
        issues.append("upper endpoint floor reserve failed")

    sigma_bound = Fraction(201, 400)
    lower_numerator = 6 * sigma_bound + 10
    if lower_numerator != Fraction(2603, 200):
        issues.append("lower endpoint deviation constant failed")
    upper_braces = 6 * (50 + Fraction(51, 100)) + 10
    if upper_braces != Fraction(15_653, 50):
        issues.append("upper endpoint deviation constant failed")


def amplitude(j: int, u: mp.mpf, sigma: mp.mpf, t: mp.mpf) -> mp.mpf:
    lam = mp.log(u)
    return lam**j * mp.e ** (t * lam**2 / 4 - sigma * lam)


def amplitude_prime(j: int, u: mp.mpf, sigma: mp.mpf, t: mp.mpf) -> mp.mpf:
    lam = mp.log(u)
    exponential = mp.e ** (t * lam**2 / 4 - sigma * lam)
    polynomial = (j * lam ** (j - 1) if j else 0)
    polynomial += (t * lam / 2 - sigma) * lam**j
    return exponential * polynomial / u


def amplitude_second(j: int, u: mp.mpf, sigma: mp.mpf, t: mp.mpf) -> mp.mpf:
    lam = mp.log(u)
    g1 = t * lam / 2 - sigma
    g2 = t / 2
    exponential = mp.e ** (t * lam**2 / 4 - sigma * lam)
    polynomial = j * (j - 1) * lam ** (j - 2) if j >= 2 else mp.mpf(0)
    if j:
        polynomial += j * (2 * g1 - 1) * lam ** (j - 1)
    polynomial += (g2 + g1**2 - g1) * lam**j
    return exponential * polynomial / u**2


def endpoint_data(
    alpha: mp.mpf,
    mode: int,
    mu: mp.mpf,
    j: int,
    sigma: mp.mpf,
    t: mp.mpf,
) -> tuple[mp.mpf, mp.mpf, mp.mpf, mp.mpf]:
    u0 = alpha / mode
    v = mu / u0
    defect = v - 1 - mp.log(v)
    z = mp.sign(v - 1) * mp.sqrt(2 * defect)
    jacobian = v * z / (v - 1)
    y = mp.sqrt(alpha) * z
    b_endpoint = mp.sqrt(alpha) * amplitude(j, mu, sigma, t) * jacobian / mode
    b_zero = mp.sqrt(alpha) * amplitude(j, u0, sigma, t) / mode
    c_value = (b_endpoint - b_zero) / y
    return y, b_endpoint, b_zero, c_value


def independent_endpoint_geometry_audit(issues: list[str]) -> None:
    mp.mp.dps = 70
    alpha = mp.mpf("37.37")
    sigma = mp.mpf("0.501")
    t = mp.mpf(1) / 5000
    for mode in (2, 7, 19, 38, 71):
        for mu in (mp.mpf(1), mp.mpf(5)):
            q = alpha / mu - mode
            for j in (0, 1, 3, 5):
                y, b_endpoint, b_zero, c_value = endpoint_data(
                    alpha, mode, mu, j, sigma, t
                )
                if not close(b_endpoint / y, -amplitude(j, mu, sigma, t) / q, "1e-58"):
                    issues.append(f"endpoint b/y identity failed at r={mode}, mu={mu}, j={j}")
                rational = -amplitude(j, mu, sigma, t) / q - b_zero / y
                if not close(c_value, rational, "1e-58"):
                    issues.append(f"endpoint c identity failed at r={mode}, mu={mu}, j={j}")


def independent_stable_sum_audit(issues: list[str]) -> None:
    mp.mp.dps = 70
    m, n = 3, 18
    for x in (mp.mpf("7.37"), mp.mpf("13.21")):
        a_x = x - m + 1
        b_x = n + 1 - x
        stable_1 = mp.digamma(b_x) - mp.digamma(a_x)
        stable_2 = mp.zeta(2, a_x) + mp.zeta(2, b_x)
        stable_3 = mp.zeta(3, a_x) - mp.zeta(3, b_x)
        finite_1 = mp.fsum(1 / (x - mode) for mode in range(m, n + 1))
        finite_2 = mp.fsum(1 / (x - mode) ** 2 for mode in range(m, n + 1))
        finite_3 = mp.fsum(1 / (x - mode) ** 3 for mode in range(m, n + 1))
        cotangent = mp.cot(mp.pi * x)
        complete_1 = mp.pi * cotangent
        complete_2 = mp.pi**2 / mp.sin(mp.pi * x) ** 2
        complete_3 = mp.pi**3 * cotangent / mp.sin(mp.pi * x) ** 2
        if not close(stable_1, complete_1 - finite_1, "1e-60"):
            issues.append(f"stable S1 audit failed at x={x}")
        if not close(stable_2, complete_2 - finite_2, "1e-60"):
            issues.append(f"stable S2 audit failed at x={x}")
        if not close(stable_3, complete_3 - finite_3, "1e-60"):
            issues.append(f"stable S3 audit failed at x={x}")


def phase_value(alpha: mp.mpf, mode: int, u: mp.mpf) -> mp.mpc:
    return mp.e ** (2j * mp.pi * (alpha * mp.log(u) - mode * u))


def independent_twice_ibp_audit(issues: list[str]) -> None:
    mp.mp.dps = 60
    alpha = mp.mpf("7.37")
    upper = mp.mpf(4)
    sigma = mp.mpf("0.501")
    t = mp.mpf(1) / 5000
    kappa = 1 / (2 * mp.pi * 1j)

    for j in (0, 2):
        for mode in (-3, 0, 20):
            def q(u: mp.mpf) -> mp.mpf:
                return alpha / u - mode

            def q_prime(u: mp.mpf) -> mp.mpf:
                return -alpha / u**2

            def q_second(u: mp.mpf) -> mp.mpf:
                return 2 * alpha / u**3

            def c_operator(u: mp.mpf) -> mp.mpf:
                a0 = amplitude(j, u, sigma, t)
                a1 = amplitude_prime(j, u, sigma, t)
                a2 = amplitude_second(j, u, sigma, t)
                q0 = q(u)
                qp = q_prime(u)
                qpp = q_second(u)
                return (
                    a2 / q0**2
                    - 3 * a1 * qp / q0**3
                    - a0 * qpp / q0**3
                    + 3 * a0 * qp**2 / q0**4
                )

            direct = mp.quad(
                lambda u: amplitude(j, u, sigma, t) * phase_value(alpha, mode, u),
                [1, upper],
            )

            def second_boundary(u: mp.mpf) -> mp.mpc:
                q0 = q(u)
                d_over_q = (
                    amplitude_prime(j, u, sigma, t) / q0**2
                    - amplitude(j, u, sigma, t) * q_prime(u) / q0**3
                )
                return phase_value(alpha, mode, u) * d_over_q

            first_boundary = kappa * (
                phase_value(alpha, mode, upper)
                * amplitude(j, upper, sigma, t)
                / q(upper)
                - phase_value(alpha, mode, 1)
                * amplitude(j, 1, sigma, t)
                / q(1)
            )
            second_term = -kappa**2 * (second_boundary(upper) - second_boundary(1))
            integral_tail = kappa**2 * mp.quad(
                lambda u: phase_value(alpha, mode, u) * c_operator(u),
                [1, upper],
            )
            if not close(direct, first_boundary + second_term + integral_tail, "1e-42"):
                issues.append(f"twice-IBP transfer failed at j={j}, r={mode}")


def stable_sums(alpha: mp.mpf, mu: mp.mpf, m: int, n: int) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    x = alpha / mu
    a_x = x - m + 1
    b_x = n + 1 - x
    return (
        mp.digamma(b_x) - mp.digamma(a_x),
        mp.zeta(2, a_x) + mp.zeta(2, b_x),
        mp.zeta(3, a_x) - mp.zeta(3, b_x),
    )


def independent_package_real_part_audit(issues: list[str]) -> None:
    mp.mp.dps = 60
    alpha = mp.mpf("37.37")
    upper = mp.mpf(5)
    sigma = mp.mpf("0.501")
    t = mp.mpf(1) / 5000
    kappa = 1 / (2 * mp.pi * 1j)
    m = max(1, int(mp.floor(alpha / (2 * upper))))
    n = int(mp.ceil(2 * alpha))
    sums_1 = stable_sums(alpha, 1, m, n)
    sums_b = stable_sums(alpha, upper, m, n)

    for j in range(6):
        c_1 = mp.fsum(
            endpoint_data(alpha, mode, 1, j, sigma, t)[3]
            for mode in range(m, n + 1)
        )
        c_b = mp.fsum(
            endpoint_data(alpha, mode, upper, j, sigma, t)[3]
            for mode in range(m, n + 1)
        )
        a_1 = amplitude(j, 1, sigma, t)
        ap_1 = amplitude_prime(j, 1, sigma, t)
        a_b = amplitude(j, upper, sigma, t)
        ap_b = amplitude_prime(j, upper, sigma, t)
        cal_1 = kappa * a_1 * sums_1[0] - kappa**2 * (
            ap_1 * sums_1[1] + alpha * a_1 * sums_1[2]
        )
        cal_b_factored = kappa * a_b * sums_b[0] - kappa**2 * (
            ap_b * sums_b[1] + alpha * a_b * sums_b[2] / upper**2
        )
        lower = (mp.mpf(1) / 2 if j == 0 else 0) - cal_1 + kappa * c_1
        upper_package = a_b / 2 + cal_b_factored - kappa * c_b
        expected_lower_real = (mp.mpf(1) / 2 if j == 0 else 0) - (
            ap_1 * sums_1[1] + alpha * a_1 * sums_1[2]
        ) / (4 * mp.pi**2)
        expected_upper_real = a_b / 2 + (
            ap_b * sums_b[1] + alpha * a_b * sums_b[2] / upper**2
        ) / (4 * mp.pi**2)
        if not close(mp.re(lower), expected_lower_real, "1e-48"):
            issues.append(f"lower real-part identity failed at j={j}")
        if not close(mp.re(upper_package), expected_upper_real, "1e-48"):
            issues.append(f"upper real-part identity failed at j={j}")


def structural_audit(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("kind drifted")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")
    rows = payload.get("rows", [])
    ids = [row.get("id") for row in rows]
    if ids != EXPECTED_IDS:
        issues.append("row ids drifted")
    if len(ids) != len(set(ids)):
        issues.append("duplicate row ids")
    readiness = [row.get("readiness") for row in rows]
    if readiness.count("proved") != 23 or readiness.count("open") != 1:
        issues.append("row readiness boundary drifted")
    for key in (
        "enumerated_physical_modes",
        "oscillatory_interior_bounds",
        "endpoint_composed_h2_bounds",
        "evaluated_signed_carriers",
        "signed_flow_bounds",
        "phi_b_bounds",
    ):
        if payload.get("counts", {}).get(key) != 0:
            issues.append(f"{key} was overpromoted")
    boundary = payload.get("proof_boundary", "")
    for marker in (
        "does not bound the grouped oscillatory interior",
        "endpoint-composed h^2 remainder",
        "signed flow",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")
    try:
        rebuilt = gate.build_payload()
    except Exception as exc:  # pragma: no cover
        issues.append(f"builder replay failed: {exc}")
    else:
        if payload != rebuilt:
            issues.append("saved result differs from a fresh builder replay")


def note_audit(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for marker in (
        "# Six-Moment Morse-Fresnel Endpoint-Composition Retention Gate",
        "H_j=F_j+L_j+e(alpha log B)U_j+J_j+T_j",
        "|L_0|>4999/10000",
        "|U_j|>49A_j(B)/100",
        "one has exactly Psi'",
        "0 oscillatory interior bounds",
        "this is not a proof",
        f"python work/rh_compute/scripts/{STEM}.py",
        f"python work/rh_compute/scripts/check_{STEM}.py",
    ):
        if marker not in text:
            issues.append(f"note marker missing: {marker}")


def main() -> int:
    issues: list[str] = []
    payload = load_json(RESULT, issues)
    if payload:
        independent_source_audit(payload, issues)
        independent_symbolic_audit(issues)
        independent_rational_audit(issues)
        independent_endpoint_geometry_audit(issues)
        independent_stable_sum_audit(issues)
        independent_twice_ibp_audit(issues)
        independent_package_real_part_audit(issues)
        structural_audit(payload, issues)
    note_audit(issues)
    if issues:
        for issue in issues:
            print(f"FAIL: {issue}")
        return 1
    print(SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
