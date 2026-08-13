#!/usr/bin/env python3
"""Validate the six-moment logarithmic Morse-Fresnel endpoint reduction."""

from __future__ import annotations

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

import jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_endpoint_reduction as gate  # noqa: E402


STEM = gate.STEM
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    f"mfer_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "domain",
            "source",
            "phase",
            "coordinate",
            "jacobian",
            "series",
            "endpoints",
            "roster",
            "amplitude",
            "integral",
            "main",
            "limits",
            "residual",
            "removable_c",
            "smooth",
            "gap",
            "ibp",
            "operator",
            "sums",
            "boundary",
            "envelope",
            "convergence",
            "values",
            "route",
            "handoff",
            "pi",
            "boundary",
        ),
        start=1,
    )
]
EXPECTED_COUNTS = {
    "rows": 27,
    "bare_logarithmic_moments": 6,
    "fixed_transition_rosters": 1,
    "exact_morse_phase_identities": 1,
    "exact_transition_integrals": 6,
    "incomplete_fresnel_main_families": 6,
    "explicit_transition_residual_majorants": 6,
    "twofold_nonstationary_ibp_identities": 6,
    "stable_endpoint_sum_families": 3,
    "explicit_tail_integral_majorants": 6,
    "sharp_cutoff_jump_laws_required": 0,
    "hidden_big_o_terms_in_derived_decomposition": 0,
    "evaluated_physical_remainder_constants": 0,
    "signed_flow_bounds": 0,
    "phi_b_bounds": 0,
}
SUMMARY = (
    "validated six-moment Morse-Fresnel endpoint reduction: 27 rows, "
    "6 exact transition integrals, 6 transition majorants, "
    "3 stable endpoint sums, 6 tail majorants, "
    "0 evaluated physical constants, 0 signed flow bounds"
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
    if sp.simplify(sp.expand_complex(expression)) != 0:
        issues.append(label)


def close(
    left: mp.mpf | mp.mpc,
    right: mp.mpf | mp.mpc,
    tolerance: str = "1e-55",
) -> bool:
    return abs(left - right) <= mp.mpf(tolerance) * max(1, abs(left), abs(right))


def independent_source_audit(payload: dict, issues: list[str]) -> None:
    saved = payload.get("source_audit", {})
    for key, path in gate.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
            continue
        entry = saved.get(key, {})
        relative = str(path.relative_to(REPO_ROOT)).replace("\\", "/")
        if entry.get("path") != relative:
            issues.append(f"source path drifted: {key}")
        if entry.get("sha256") != file_hash(path):
            issues.append(f"source hash drifted: {key}")

    finite_path = gate.SOURCE_PATHS["finite_poisson_transport"]
    flow_path = gate.SOURCE_PATHS["flow_matrix_phase"]
    if finite_path.is_file():
        finite = json.loads(finite_path.read_text(encoding="utf-8"))
        counts = finite.get("counts", {})
        if counts.get("bare_logarithmic_moments") != 6:
            issues.append("finite-Poisson source lost six moments")
        if counts.get("dual_cutoff_jump_laws") != 2:
            issues.append("finite-Poisson source lost the two cutoff laws")
        if counts.get("explicit_uniform_remainder_constants") != 0:
            issues.append("finite-Poisson proof boundary drifted")
    if flow_path.is_file():
        flow = json.loads(flow_path.read_text(encoding="utf-8"))
        if flow.get("counts", {}).get("real_flow_observations") != 8:
            issues.append("flow source lost the eight-observation handoff")

    external = payload.get("external_source", {})
    if external.get("url") != gate.VANDEHEY_URL:
        issues.append("external source URL drifted")
    if external.get("usage") != "comparison and Fresnel-endpoint motivation only":
        issues.append("external source usage boundary drifted")


def independent_symbolic_audit(issues: list[str]) -> None:
    delta = sp.symbols("delta", positive=True)
    v = 1 + delta
    h = v - 1 - sp.log(v)
    z = sp.sqrt(2 * h)
    jacobian = v * z / (v - 1)
    require_zero(
        sp.series(z, delta, 0, 4).removeO()
        - (delta - delta**2 / 3 + 7 * delta**3 / 36),
        "independent Morse-coordinate series failed",
        issues,
    )
    require_zero(
        sp.series(jacobian, delta, 0, 3).removeO()
        - (1 + 2 * delta / 3 - 5 * delta**2 / 36),
        "independent Morse-Jacobian series failed",
        issues,
    )

    u, alpha, mode = sp.symbols("u alpha r", positive=True)
    scaled = sp.symbols("v", positive=True)
    saddle = alpha / mode
    phase = alpha * sp.log(saddle * scaled) - mode * saddle * scaled
    phase0 = alpha * (sp.log(saddle) - 1)
    require_zero(
        phase - phase0 + alpha * (scaled - 1 - sp.log(scaled)),
        "independent exact Morse phase failed",
        issues,
    )

    amplitude = sp.Function("A")(u)
    q = alpha / u - mode
    direct = sp.diff(sp.diff(amplitude / q, u) / q, u)
    q_prime = sp.diff(q, u)
    q_second = sp.diff(q, u, 2)
    expanded = (
        sp.diff(amplitude, u, 2) / q**2
        - 3 * sp.diff(amplitude, u) * q_prime / q**3
        - amplitude * q_second / q**3
        + 3 * amplitude * q_prime**2 / q**4
    )
    require_zero(
        direct - expanded,
        "independent twice-integrated tail operator failed",
        issues,
    )


def independent_morse_audit(issues: list[str]) -> None:
    mp.mp.dps = 80

    def defect(value: mp.mpf) -> mp.mpf:
        return value - 1 - mp.log(value)

    def coordinate(value: mp.mpf) -> mp.mpf:
        if value == 1:
            return mp.mpf(0)
        return mp.sign(value - 1) * mp.sqrt(2 * defect(value))

    samples = [
        mp.mpf("0.25"),
        mp.mpf("0.5"),
        mp.mpf("0.9"),
        mp.mpf("1.1"),
        mp.mpf("2"),
        mp.mpf("4"),
    ]
    values = [coordinate(value) for value in samples]
    if any(defect(value) <= 0 for value in samples):
        issues.append("global Morse defect lost positivity")
    if any(left >= right for left, right in zip(values, values[1:])):
        issues.append("global Morse coordinate lost monotonicity")

    for value in samples:
        z_value = coordinate(value)
        analytic_derivative = (value - 1) / (value * z_value)
        numerical_derivative = mp.diff(coordinate, value)
        jacobian = value * z_value / (value - 1)
        if not close(numerical_derivative, analytic_derivative, "1e-60"):
            issues.append(f"Morse derivative failed at v={value}")
        if not close(jacobian, 1 / analytic_derivative, "1e-70"):
            issues.append(f"Morse Jacobian failed at v={value}")

    alpha = mp.mpf(113) / 7
    mode = mp.mpf(5)
    saddle = alpha / mode
    phase0 = alpha * (mp.log(saddle) - 1)
    for value in samples:
        u = saddle * value
        y = mp.sqrt(alpha) * coordinate(value)
        phase = alpha * mp.log(u) - mode * u
        if not close(phase, phase0 - y**2 / 2, "1e-65"):
            issues.append(f"finite exact Morse phase failed at v={value}")


def independent_fresnel_audit(issues: list[str]) -> None:
    mp.mp.dps = 80
    root_two = mp.sqrt(2)
    imaginary = mp.mpc(0, 1)

    def primitive(value: mp.mpf) -> mp.mpc:
        scaled = root_two * value
        return (
            mp.fresnelc(scaled) - imaginary * mp.fresnels(scaled)
        ) / root_two

    half = (
        mp.fresnelc(mp.inf) - imaginary * mp.fresnels(mp.inf)
    ) / root_two
    expected_half = mp.e ** (-imaginary * mp.pi / 4) / 2
    if not close(half, expected_half, "1e-70"):
        issues.append("half-line Fresnel normalization failed")
    if not close(2 * half, mp.e ** (-imaginary * mp.pi / 4), "1e-70"):
        issues.append("full-line Fresnel normalization failed")

    lower = mp.mpf("-1.2")
    upper = mp.mpf("0.85")
    direct = mp.quad(
        lambda y: mp.e ** (-mp.pi * imaginary * y**2),
        [lower, 0, upper],
    )
    transformed = primitive(upper) - primitive(lower)
    if not close(direct, transformed, "1e-65"):
        issues.append("finite incomplete-Fresnel formula failed")


def independent_roster_audit(issues: list[str]) -> None:
    ranges = (
        (mp.mpf("0.7"), mp.mpf("1.3"), 8),
        (mp.mpf("11"), mp.mpf("17"), 4),
        (mp.mpf("2500"), mp.mpf("2600"), 47),
    )
    for alpha_minus, alpha_plus, endpoint in ranges:
        lower = max(1, int(mp.floor(alpha_minus / (2 * endpoint))))
        upper = int(mp.ceil(2 * alpha_plus))
        for step in range(101):
            alpha = alpha_minus + (alpha_plus - alpha_minus) * step / 100
            active_low = int(mp.ceil(alpha / endpoint))
            active_high = int(mp.floor(alpha))
            for mode in range(max(1, active_low), active_high + 1):
                if not lower <= mode <= upper:
                    issues.append("fixed transition roster missed a saddle mode")
                    return
        low_mode = lower - 1
        if low_mode >= 1:
            low_gap = alpha_minus / endpoint - low_mode
            if low_gap <= alpha_minus / (2 * endpoint):
                issues.append("low nonstationary roster gap failed")
        high_mode = upper + 1
        if high_mode - alpha_plus <= alpha_plus:
            issues.append("high nonstationary roster gap failed")


def independent_stable_sum_audit(issues: list[str]) -> None:
    mp.mp.dps = 90
    m = 3
    n = 19
    x = mp.mpf("7.25")
    a = x - m + 1
    b = n + 1 - x
    s1 = mp.digamma(b) - mp.digamma(a)
    s2 = mp.zeta(2, a) + mp.zeta(2, b)
    s3 = mp.zeta(3, a) - mp.zeta(3, b)
    cotangent = mp.pi * mp.cot(mp.pi * x) - mp.fsum(
        1 / (x - mode) for mode in range(m, n + 1)
    )
    if not close(s1, cotangent, "1e-75"):
        issues.append("digamma/cotangent endpoint identity failed")

    terms = 6000
    partial1 = mp.fsum(
        1 / (a + k) - 1 / (b + k) for k in range(terms)
    )
    partial2 = mp.fsum(
        1 / (a + k) ** 2 + 1 / (b + k) ** 2 for k in range(terms)
    )
    partial3 = mp.fsum(
        1 / (a + k) ** 3 - 1 / (b + k) ** 3 for k in range(terms)
    )
    c = min(a, b) + terms
    s1_tail_bound = abs(b - a) * (1 / c + 1 / c**2)
    s2_tail_bound = sum(
        1 / (base + terms) + 1 / (base + terms) ** 2
        for base in (a, b)
    )
    s3_tail_bound = sum(
        1 / (2 * (base + terms) ** 2) + 1 / (base + terms) ** 3
        for base in (a, b)
    )
    if abs(s1 - partial1) > s1_tail_bound:
        issues.append("direct symmetric S1 tail check failed")
    if abs(s2 - partial2) > s2_tail_bound:
        issues.append("direct absolute S2 tail check failed")
    if abs(s3 - partial3) > s3_tail_bound:
        issues.append("direct absolute S3 tail check failed")

    integer_x = mp.mpf(11)
    integer_s1 = mp.digamma(n + 1 - integer_x) - mp.digamma(
        integer_x - m + 1
    )
    if integer_s1 != 0:
        issues.append("removable integer S1 value failed")
    with mp.workdps(180):
        nearby = mp.mpf(11) + mp.mpf("1e-30")
        nearby_stable = mp.digamma(n + 1 - nearby) - mp.digamma(
            nearby - m + 1
        )
        nearby_cotangent = mp.pi * mp.cot(mp.pi * nearby) - mp.fsum(
            1 / (nearby - mode) for mode in range(m, n + 1)
        )
        if not close(nearby_stable, nearby_cotangent, "1e-100"):
            issues.append("near-integer removable S1 chart failed")


def independent_two_ibp_audit(issues: list[str]) -> None:
    mp.mp.dps = 75
    imaginary = mp.mpc(0, 1)
    kappa = 1 / (2 * mp.pi * imaginary)
    alpha = mp.mpf(10)
    lower = mp.mpf(1)
    upper = mp.mpf(3)

    def amplitude(u: mp.mpf) -> mp.mpf:
        return (1 + u) * mp.e ** (-u / 5)

    def amplitude_prime(u: mp.mpf) -> mp.mpf:
        return (4 - u) * mp.e ** (-u / 5) / 5

    def amplitude_second(u: mp.mpf) -> mp.mpf:
        return (u - 9) * mp.e ** (-u / 5) / 25

    def segmented_integral(function):
        points = [lower + (upper - lower) * k / 96 for k in range(97)]
        return mp.quad(function, points)

    for mode in (mp.mpf(-2), mp.mpf(30)):
        def phase(u: mp.mpf) -> mp.mpc:
            value = alpha * mp.log(u) - mode * u
            return mp.e ** (2 * mp.pi * imaginary * value)

        def q(u: mp.mpf) -> mp.mpf:
            return alpha / u - mode

        def q_prime(u: mp.mpf) -> mp.mpf:
            return -alpha / u**2

        def q_second(u: mp.mpf) -> mp.mpf:
            return 2 * alpha / u**3

        def first_operator(u: mp.mpf) -> mp.mpf:
            return (
                amplitude_prime(u) / q(u)
                - amplitude(u) * q_prime(u) / q(u) ** 2
            )

        def second_operator(u: mp.mpf) -> mp.mpf:
            return (
                amplitude_second(u) / q(u) ** 2
                - 3 * amplitude_prime(u) * q_prime(u) / q(u) ** 3
                - amplitude(u) * q_second(u) / q(u) ** 3
                + 3 * amplitude(u) * q_prime(u) ** 2 / q(u) ** 4
            )

        def bracket(function) -> mp.mpc:
            return phase(upper) * function(upper) - phase(lower) * function(lower)

        direct = segmented_integral(lambda u: phase(u) * amplitude(u))
        rebuilt = (
            kappa * bracket(lambda u: amplitude(u) / q(u))
            - kappa**2 * bracket(lambda u: first_operator(u) / q(u))
            + kappa**2
            * segmented_integral(lambda u: phase(u) * second_operator(u))
        )
        if not close(direct, rebuilt, "1e-50"):
            issues.append(f"finite twice-IBP identity failed for r={mode}")


def independent_transition_residual_audit(issues: list[str]) -> None:
    mp.mp.dps = 75
    imaginary = mp.mpc(0, 1)
    kappa = 1 / (2 * mp.pi * imaginary)
    lower = mp.mpf("-1.3")
    upper = mp.mpf("2.1")

    def amplitude(y: mp.mpf) -> mp.mpf:
        return 1 + y / 7 + y**2 / 5 + y**3 / 11

    def amplitude_prime(y: mp.mpf) -> mp.mpf:
        return mp.mpf(1) / 7 + 2 * y / 5 + 3 * y**2 / 11

    def amplitude_second(y: mp.mpf) -> mp.mpf:
        return mp.mpf(2) / 5 + 6 * y / 11

    b0 = amplitude(0)

    def oscillation(y: mp.mpf) -> mp.mpc:
        return mp.e ** (-mp.pi * imaginary * y**2)

    def c(y: mp.mpf) -> mp.mpf:
        return mp.mpf(1) / 7 + y / 5 + y**2 / 11

    def c_prime(y: mp.mpf) -> mp.mpf:
        return mp.mpf(1) / 5 + 2 * y / 11

    for sample in (lower, mp.mpf("-0.2"), 0, mp.mpf("0.4"), upper):
        if not close(sample * c(sample), amplitude(sample) - b0, "1e-70"):
            issues.append("removable c chart failed")
        if not close(c_prime(sample), mp.diff(c, sample), "1e-70"):
            issues.append("removable c-prime chart failed")

    direct = mp.quad(
        lambda y: (amplitude(y) - b0) * oscillation(y),
        [lower, 0, upper],
    )
    rebuilt = (
        -kappa * (c(upper) * oscillation(upper) - c(lower) * oscillation(lower))
        + kappa
        * mp.quad(lambda y: c_prime(y) * oscillation(y), [lower, 0, upper])
    )
    if not close(direct, rebuilt, "1e-55"):
        issues.append("finite transition residual identity failed")
    majorant = (
        abs(c(lower))
        + abs(c(upper))
        + mp.quad(lambda y: abs(c_prime(y)), [lower, 0, upper])
    ) / (2 * mp.pi)
    if abs(direct) > majorant + mp.mpf("1e-60"):
        issues.append("transition residual majorant failed")


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
    if readiness.count("open") != 1 or readiness.count("proved") != 26:
        issues.append("row readiness boundary drifted")
    for key in (
        "sharp_cutoff_jump_laws_required",
        "hidden_big_o_terms_in_derived_decomposition",
        "evaluated_physical_remainder_constants",
        "signed_flow_bounds",
        "phi_b_bounds",
    ):
        if payload.get("counts", {}).get(key) != 0:
            issues.append(f"{key} was overpromoted")
    boundary = payload.get("proof_boundary", "")
    for marker in (
        "exact fixed-roster Morse-Fresnel representation",
        "does not evaluate",
        "signed flow estimate",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")
    try:
        rebuilt = gate.build_payload()
    except Exception as exc:  # pragma: no cover - reported as gate failure
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
        "# Six-Moment Logarithmic Morse-Fresnel Endpoint Reduction",
        gate.VANDEHEY_URL,
        "exact endpoint",
        "psi(b_x)-psi(a_x)",
        "0 evaluated physical remainder constants",
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
        independent_morse_audit(issues)
        independent_fresnel_audit(issues)
        independent_roster_audit(issues)
        independent_stable_sum_audit(issues)
        independent_two_ibp_audit(issues)
        independent_transition_residual_audit(issues)
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
