#!/usr/bin/env python3
"""Validate the Morse-Fresnel aggregate-scaling route gate."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = REPO_ROOT / "work/rh_compute/scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_aggregate_scaling_gate as gate  # noqa: E402


STEM = gate.STEM
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    f"mfas_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "domain",
            "cutoff",
            "alpha",
            "geometry",
            "p_upper",
            "q_formula",
            "q_floor",
            "collar",
            "amplitude",
            "b_derivative",
            "c_endpoint",
            "mode_count",
            "majorant_floor",
            "h2_guard",
            "tail_coordinates",
            "hurwitz",
            "zeta",
            "tail_reduction",
            "decay",
            "derivatives",
            "constants",
            "route",
            "handoff",
            "boundary",
        ),
        start=1,
    )
]
EXPECTED_COUNTS = {
    "rows": 24,
    "bare_logarithmic_moments": 6,
    "certified_transition_collars": 1,
    "enumerated_transition_modes": 0,
    "termwise_absolute_transition_majorant_floors": 1,
    "majorant_h2_separation_factors": 1,
    "explicit_outside_roster_tail_bounds": 6,
    "explicit_tail_h2_scaling_bounds": 6,
    "bounded_endpoint_functionals": 0,
    "grouped_roster_bounds": 0,
    "evaluated_full_physical_remainders": 0,
    "signed_flow_bounds": 0,
    "phi_b_bounds": 0,
}
EXPECTED_C = (
    Fraction(234_803, 1_225),
    Fraction(1_145_600, 2_401),
    Fraction(232_001_200, 117_649),
    Fraction(69_600_360_000, 5_764_801),
    Fraction(27_840_144_000_000, 282_475_249),
    Fraction(13_920_072_000_000_000, 13_841_287_201),
)
EXPECTED_K = (6, 14, 55, 336, 2738, 27936)
SUMMARY = (
    "validated Morse-Fresnel aggregate-scaling gate: 24 rows, "
    "1 transition collar, 0 enumerated modes, majorant floor 91/96800, "
    "6 h^2 tail bounds, 0 grouped roster bounds, 0 signed flow bounds"
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
    if sp.simplify(expression) != 0:
        issues.append(label)


def close(left: mp.mpf, right: mp.mpf, tolerance: str = "1e-60") -> bool:
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
        "morse_fresnel_endpoint": {
            "exact_transition_integrals": 6,
            "explicit_tail_integral_majorants": 6,
        },
        "physical_q1": {"physical_q1_parameter_laws": 4},
        "growing_terminal_prefix": {"finite_height_current_theorems": 1},
    }
    for source, checks in expected.items():
        path = gate.SOURCE_PATHS[source]
        if not path.is_file():
            continue
        counts = json.loads(path.read_text(encoding="utf-8")).get("counts", {})
        for key, value in checks.items():
            if counts.get(key) != value:
                issues.append(f"source count drifted: {source}.{key}")


def defect_coefficient(order: int) -> int:
    return 16 * 3**order - 2 ** (order - 1) * (50 * order - 54) - 52


def independent_geometry_audit(issues: list[str]) -> None:
    s = sp.symbols("s", positive=True)
    exponential = sp.exp(s)
    p_geometry = 2 * (exponential * (s - 1) + 1) / (exponential - 1) ** 2
    q_geometry = sp.simplify((p_geometry - sp.diff(p_geometry, s)) / 2)
    q_expected = (
        2 * s * sp.exp(2 * s)
        - 3 * sp.exp(2 * s)
        + 4 * sp.exp(s)
        - 1
    ) / (sp.exp(s) - 1) ** 3
    require_zero(q_geometry - q_expected, "independent Q identity failed", issues)

    p_difference = (
        (exponential - 1) ** 2
        - 2 * (exponential * (s - 1) + 1)
    )
    require_zero(
        p_difference - (sp.exp(2 * s) - 2 * s * sp.exp(s) - 1),
        "independent P numerator identity failed",
        issues,
    )
    require_zero(
        p_difference
        - sp.expand((2 * sp.exp(s) * (sp.sinh(s) - s)).rewrite(sp.exp)),
        "independent P sinh factor failed",
        issues,
    )

    q_numerator = (
        2 * s * sp.exp(2 * s)
        - 3 * sp.exp(2 * s)
        + 4 * sp.exp(s)
        - 1
    )
    floor_numerator = sp.expand(
        25 * q_numerator - 16 * (sp.exp(s) - 1) ** 3
    )
    expected_numerator = (
        50 * s * sp.exp(2 * s)
        - 16 * sp.exp(3 * s)
        - 27 * sp.exp(2 * s)
        + 52 * sp.exp(s)
        - 9
    )
    require_zero(
        floor_numerator - expected_numerator,
        "independent Q-floor numerator failed",
        issues,
    )

    series = sp.series(expected_numerator, s, 0, 8).removeO().expand()
    expected_series = Fraction(2, 3) * s**3
    for order in range(4, 8):
        expected_series -= Fraction(defect_coefficient(order), math.factorial(order)) * s**order
    require_zero(series - expected_series, "Q-floor Taylor coefficients failed", issues)

    if defect_coefficient(4) != 76:
        issues.append("D_4 certificate failed")
    for order in range(4, 1001):
        difference = defect_coefficient(order + 1) - 2 * defect_coefficient(order)
        expected_difference = 16 * 3**order - 50 * 2**order + 52
        if difference != expected_difference or difference <= 0:
            issues.append(f"coefficient-defect induction failed at n={order}")
            break
    h4 = 16 * 3**4 - 50 * 2**4 + 52
    if h4 != 548:
        issues.append("coefficient-defect base positivity failed")

    tail = (
        Fraction(16_000)
        * Fraction(3, 10) ** 8
        / math.factorial(8)
        * Fraction(30, 29)
    )
    if tail != Fraction(2_187, 81_200_000):
        issues.append("Q-floor Taylor tail allowance failed")
    reserve = (
        Fraction(2, 3)
        - Fraction(19, 60)
        - Fraction(35, 600)
        - Fraction(187, 36_000)
        - Fraction(1_333, 4_200_000)
        - tail
    )
    if reserve != Fraction(8_364_091, 29_232_000) or reserve <= 0:
        issues.append("Q-floor exact rational reserve failed")


def independent_geometry_numeric_audit(issues: list[str]) -> None:
    mp.mp.dps = 90

    def p_geometry(value: mp.mpf) -> mp.mpf:
        if value == 0:
            return mp.mpf(1)
        exponential = mp.exp(value)
        return 2 * (exponential * (value - 1) + 1) / (exponential - 1) ** 2

    def q_geometry(value: mp.mpf) -> mp.mpf:
        if value == 0:
            return mp.mpf(2) / 3
        exponential = mp.exp(value)
        return (
            2 * value * exponential**2
            - 3 * exponential**2
            + 4 * exponential
            - 1
        ) / (exponential - 1) ** 3

    samples = ("0", "0.0001", "0.001", "0.01", "0.05", "0.1")
    for text in samples:
        value = mp.mpf(text)
        p_value = p_geometry(value)
        q_value = q_geometry(value)
        if p_value > 1 + mp.mpf("1e-75"):
            issues.append(f"P<=1 failed numerically at s={text}")
        if q_value < mp.mpf(16) / 25 - mp.mpf("1e-75"):
            issues.append(f"Q>=16/25 failed numerically at s={text}")
    if not close(p_geometry(mp.mpf(0)), mp.mpf(1)):
        issues.append("P removable value failed")
    if not close(q_geometry(mp.mpf(0)), mp.mpf(2) / 3):
        issues.append("Q removable value failed")


def independent_transition_audit(issues: list[str]) -> None:
    alpha, mode, v, t, sigma = sp.symbols(
        "alpha r v t sigma", positive=True
    )
    jacobian = sp.Function("J")(v)
    u = alpha * v / mode
    amplitude = sp.exp(t * sp.log(u) ** 2 / 4 - sigma * sp.log(u))
    transformed = sp.sqrt(alpha) * amplitude * jacobian / mode
    derivative_y = sp.diff(transformed, v) * jacobian / sp.sqrt(alpha)
    expected = amplitude / mode * (
        jacobian * sp.diff(jacobian, v)
        + jacobian**2 / v * (t * sp.log(u) / 2 - sigma)
    )
    require_zero(
        derivative_y - expected,
        "independent transformed-amplitude chain rule failed",
        issues,
    )

    amplitude_floor = Fraction(3_799, 4_000)
    bracket_floor = Fraction(16, 25) - Fraction(201, 400)
    derivative_floor = amplitude_floor * bracket_floor
    advertised = Fraction(13, 100)
    harmonic_floor = Fraction(1, 22)
    majorant_floor = advertised * harmonic_floor * Fraction(7, 44)
    if bracket_floor != Fraction(11, 80):
        issues.append("transition bracket floor failed")
    if derivative_floor != Fraction(41_789, 320_000):
        issues.append("transition derivative exact floor failed")
    if derivative_floor <= advertised:
        issues.append("transition derivative reserve failed")
    if majorant_floor != Fraction(91, 96_800):
        issues.append("termwise majorant floor failed")
    if majorant_floor * 2**50 <= 10**12:
        issues.append("h^2 scale-separation floor failed")

    mp.mp.dps = 80
    for alpha_value in (
        mp.mpf(22) + mp.mpf(1) / 3,
        mp.mpf(1000) + mp.mpf(1) / 3,
        mp.mpf(2) ** 50 + mp.mpf(1) / 3,
    ):
        lower = int(mp.ceil(alpha_value * mp.exp(-mp.mpf(1) / 10)))
        upper = int(mp.floor(alpha_value))
        count = upper - lower + 1
        harmonic = mp.digamma(upper + 1) - mp.digamma(lower)
        if count < alpha_value / 22:
            issues.append(f"representative collar count failed at alpha={alpha_value}")
        if harmonic < mp.mpf(1) / 22:
            issues.append(f"representative collar harmonic mass failed at alpha={alpha_value}")


def independent_hurwitz_audit(issues: list[str]) -> None:
    mp.mp.dps = 80
    for order in range(2, 7):
        for q_value in (mp.mpf(1), mp.mpf("1.25"), mp.mpf(3), mp.mpf(17)):
            zeta_value = mp.zeta(order, q_value)
            elementary = q_value ** (-order) + q_value ** (1 - order) / (order - 1)
            coarse = 2 * q_value ** (1 - order)
            if zeta_value > elementary + mp.mpf("1e-70"):
                issues.append(f"Hurwitz integral bound failed for k={order}, q={q_value}")
            if elementary > coarse + mp.mpf("1e-70"):
                issues.append(f"coarse Hurwitz bound failed for k={order}, q={q_value}")

    alpha = mp.mpf(100)
    lower_mode = 1
    upper_mode = 200
    for u in (mp.mpf(1), mp.mpf(2), mp.mpf(5), mp.mpf(10)):
        x = alpha / u
        a_x = x - lower_mode + 1
        b_x = upper_mode + 1 - x
        if a_x < alpha / (2 * u) or b_x < alpha:
            issues.append(f"representative roster distance failed at u={u}")
            continue
        z2 = mp.zeta(2, a_x) + mp.zeta(2, b_x)
        z3 = mp.zeta(3, a_x) + mp.zeta(3, b_x)
        z4 = mp.zeta(4, a_x) + mp.zeta(4, b_x)
        if z2 > 6 * u / alpha:
            issues.append(f"Z_2 envelope failed at u={u}")
        if z3 > 10 * u**2 / alpha**2:
            issues.append(f"Z_3 envelope failed at u={u}")
        if z4 > 18 * u**3 / alpha**3:
            issues.append(f"Z_4 envelope failed at u={u}")


def logarithmic_moment(order: int) -> Fraction:
    if order < 0:
        return Fraction(0)
    return Fraction(math.factorial(order)) / Fraction(49, 100) ** (order + 1)


def independent_tail_constant_audit(payload: dict, issues: list[str]) -> None:
    lam, t, sigma = sp.symbols("lambda t sigma", positive=True)
    order_symbol = sp.symbols("j", integer=True, nonnegative=True)
    g = t * lam**2 / 4 - sigma * lam
    amplitude = lam**order_symbol * sp.exp(g)
    first = sp.diff(amplitude, lam)
    second_minus_first = sp.diff(amplitude, lam, 2) - first
    expected_first = sp.exp(g) * (
        order_symbol * lam ** (order_symbol - 1)
        + (t * lam / 2 - sigma) * lam**order_symbol
    )
    expected_second_minus_first = sp.exp(g) * (
        order_symbol * (order_symbol - 1) * lam ** (order_symbol - 2)
        + (2 * order_symbol * (t * lam / 2 - sigma) - order_symbol)
        * lam ** (order_symbol - 1)
        + (
            t / 2
            + (t * lam / 2 - sigma) ** 2
            - (t * lam / 2 - sigma)
        )
        * lam**order_symbol
    )
    require_zero(first - expected_first, "amplitude first derivative failed", issues)
    require_zero(
        second_minus_first - expected_second_minus_first,
        "amplitude second-minus-first derivative failed",
        issues,
    )

    constants = payload.get("tail_certificate", {}).get("constants", [])
    if len(constants) != 6:
        issues.append("saved six-tail constant table is incomplete")
        constants = []
    derived_c: list[Fraction] = []
    derived_k: list[int] = []
    for order in range(6):
        d2_minus_d = (
            order * (order - 1) * logarithmic_moment(order - 2)
            + Fraction(101 * order, 50) * logarithmic_moment(order - 1)
            + Fraction(3_851, 5_000) * logarithmic_moment(order)
        )
        first_bound = (
            order * logarithmic_moment(order - 1)
            + Fraction(51, 100) * logarithmic_moment(order)
        )
        value_bound = logarithmic_moment(order)
        aggregate = 6 * d2_minus_d + 30 * first_bound + 74 * value_bound
        alpha_coefficient = aggregate / 36
        ceiling = math.ceil(alpha_coefficient)
        derived_c.append(aggregate)
        derived_k.append(ceiling)
        if constants:
            row = constants[order]
            if row.get("j") != order:
                issues.append(f"tail constant order drifted at j={order}")
            if row.get("C_j") != str(aggregate):
                issues.append(f"tail C_j drifted at j={order}")
            if row.get("K_j_ceiling") != ceiling:
                issues.append(f"tail K_j drifted at j={order}")
            expected_bound = f"<{ceiling}/alpha<={2 * ceiling}*h^2"
            if row.get("outside_roster_integral_tail_bound") != expected_bound:
                issues.append(f"tail h^2 statement drifted at j={order}")
    if tuple(derived_c) != EXPECTED_C:
        issues.append(f"independent C_j constants drifted: {derived_c}")
    if tuple(derived_k) != EXPECTED_K:
        issues.append(f"independent K_j ceilings drifted: {derived_k}")


def independent_finite_amplitude_audit(issues: list[str]) -> None:
    mp.mp.dps = 60
    t = mp.mpf(1) / 5000
    sigma = mp.mpf("0.501")
    upper = mp.mpf(49)
    points = [mp.mpf(0), mp.mpf(1), mp.mpf(5), mp.mpf(15), mp.mpf(30), upper]

    for order in range(6):
        def pieces(value: mp.mpf) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
            exponential = mp.exp(t * value**2 / 4 - sigma * value)
            power = value**order
            first_power = order * value ** (order - 1) if order >= 1 else mp.mpf(0)
            second_power = (
                order * (order - 1) * value ** (order - 2)
                if order >= 2
                else mp.mpf(0)
            )
            g_prime = t * value / 2 - sigma
            g_second = t / 2
            amplitude = exponential * power
            derivative = exponential * (first_power + g_prime * power)
            second = exponential * (
                second_power
                + 2 * g_prime * first_power
                + (g_second + g_prime**2) * power
            )
            return amplitude, derivative, second - derivative

        value_integral = mp.quad(lambda value: abs(pieces(value)[0]), points)
        first_integral = mp.quad(lambda value: abs(pieces(value)[1]), points)
        second_integral = mp.quad(lambda value: abs(pieces(value)[2]), points)
        value_bound = mp.mpf(logarithmic_moment(order).numerator) / logarithmic_moment(order).denominator
        first_fraction = (
            order * logarithmic_moment(order - 1)
            + Fraction(51, 100) * logarithmic_moment(order)
        )
        second_fraction = (
            order * (order - 1) * logarithmic_moment(order - 2)
            + Fraction(101 * order, 50) * logarithmic_moment(order - 1)
            + Fraction(3_851, 5_000) * logarithmic_moment(order)
        )
        first_bound = mp.mpf(first_fraction.numerator) / first_fraction.denominator
        second_bound = mp.mpf(second_fraction.numerator) / second_fraction.denominator
        if value_integral > value_bound:
            issues.append(f"finite amplitude envelope failed at j={order}")
        if first_integral > first_bound:
            issues.append(f"finite first-derivative envelope failed at j={order}")
        if second_integral > second_bound:
            issues.append(f"finite second-derivative envelope failed at j={order}")


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
        "enumerated_transition_modes",
        "bounded_endpoint_functionals",
        "grouped_roster_bounds",
        "evaluated_full_physical_remainders",
        "signed_flow_bounds",
        "phi_b_bounds",
    ):
        if payload.get("counts", {}).get(key) != 0:
            issues.append(f"{key} was overpromoted")

    boundary = payload.get("proof_boundary", "")
    for marker in (
        "fixed positive floor",
        "does not lower-bound the true signed transition residual",
        "grouped roster estimate",
        "signed flow bound",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")
    route = payload.get("exact", {}).get("route_decision", "")
    for marker in ("retire summing", "Group roster modes", "h^2 scaling"):
        if marker not in route:
            issues.append(f"route marker missing: {marker}")

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
        "# Six-Moment Morse-Fresnel Aggregate-Scaling Gate",
        "0 enumerated transition",
        "91/96800",
        "10^12*h^2",
        "Z_2<=6u/alpha",
        "j=5:",
        "K_j=27936",
        "majorant being proposed",
        "signed residual",
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
        independent_geometry_audit(issues)
        independent_geometry_numeric_audit(issues)
        independent_transition_audit(issues)
        independent_hurwitz_audit(issues)
        independent_tail_constant_audit(payload, issues)
        independent_finite_amplitude_audit(issues)
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
