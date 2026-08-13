#!/usr/bin/env python3
"""Validate the Morse-Fresnel endpoint-coherence and Abel route gate."""

from __future__ import annotations

import cmath
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

import jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_endpoint_coherence_abel_gate as gate  # noqa: E402


STEM = gate.STEM
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    f"mfec_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "domain",
            "gauge",
            "first",
            "second",
            "curvature",
            "resonance",
            "split",
            "abel",
            "vertical",
            "outer",
            "partial",
            "weighted",
            "scale",
            "core",
            "cfloor",
            "subblock",
            "residual",
            "lower",
            "upper",
            "trace",
            "route",
            "handoff",
            "pi",
            "boundary",
        ),
        start=1,
    )
]
EXPECTED_COUNTS = {
    "rows": 24,
    "exact_phase_difference_identities": 2,
    "curvature_envelopes": 1,
    "phase_partial_sum_bounds": 1,
    "weighted_abel_bounds": 1,
    "certified_terminal_subblocks": 1,
    "endpoint_phase_collapses": 2,
    "coherent_endpoint_trace_floors": 1,
    "enumerated_physical_modes": 0,
    "full_grouped_roster_bounds": 0,
    "endpoint_composed_bounds": 0,
    "evaluated_full_physical_remainders": 0,
    "signed_flow_bounds": 0,
    "phi_b_bounds": 0,
}
SUMMARY = (
    "validated Morse-Fresnel endpoint-coherence Abel gate: 24 rows, "
    "2 phase differences, phase sums <3sqrt(alpha), "
    "1 h-scale terminal subblock, 2 coherent endpoint phases, "
    "0 endpoint-composed bounds, 0 signed flow bounds"
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
        relative = str(path.relative_to(REPO_ROOT)).replace("\\", "/")
        entry = saved.get(key, {})
        if entry.get("path") != relative:
            issues.append(f"source path drifted: {key}")
        if entry.get("sha256") != file_hash(path):
            issues.append(f"source hash drifted: {key}")

    expected = {
        "aggregate_scaling": {
            "certified_transition_collars": 1,
            "explicit_tail_h2_scaling_bounds": 6,
            "grouped_roster_bounds": 0,
        },
        "morse_fresnel_endpoint": {"exact_transition_integrals": 6},
        "finite_poisson_transport": {
            "bare_logarithmic_moments": 6,
            "stationary_transport_identities": 5,
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
    alpha, r = sp.symbols("alpha r", positive=True)
    gauge = alpha * (sp.log(alpha / r) - 1) + r
    first = sp.expand_log(gauge.subs(r, r + 1) - gauge, force=True)
    first_expected = 1 - alpha * (sp.log(r + 1) - sp.log(r))
    require_zero(first - first_expected, "independent first difference failed", issues)
    second = sp.expand_log(
        gauge.subs(r, r + 2) - 2 * gauge.subs(r, r + 1) + gauge,
        force=True,
    )
    second_expected = alpha * (
        2 * sp.log(r + 1) - sp.log(r) - sp.log(r + 2)
    )
    require_zero(second - second_expected, "independent second difference failed", issues)

    require_zero(
        (r + 1) ** 2 / (r * (r + 2)) - (1 + 1 / (r * (r + 2))),
        "compact curvature identity failed",
        issues,
    )

    v1 = r / alpha
    h1 = v1 - 1 - sp.log(v1)
    saddle = alpha * (sp.log(alpha / r) - 1)
    require_zero(
        sp.expand_log(saddle - alpha * h1 + r, force=True),
        "independent lower endpoint phase failed",
        issues,
    )
    B = sp.symbols("B", positive=True)
    vB = B * r / alpha
    hB = vB - 1 - sp.log(vB)
    require_zero(
        sp.expand_log(
            saddle - alpha * hB - (alpha * sp.log(B) - r * B),
            force=True,
        ),
        "independent upper endpoint phase failed",
        issues,
    )


def independent_rational_audit(issues: list[str]) -> None:
    margin = Fraction(64**2, 4) - Fraction(301 * 64, 20) - 1
    if margin != Fraction(299, 5) or margin <= 0:
        issues.append("terminal phase-spread margin failed")
    derivative_at_64 = Fraction(64, 2) - Fraction(301, 20)
    if derivative_at_64 <= 0:
        issues.append("terminal phase-spread monotonicity failed")
    sector_floor = Fraction(99, 100)
    core_floor = sector_floor * Fraction(13, 100) * Fraction(1, 20)
    if core_floor != Fraction(1287, 200_000):
        issues.append("terminal core floor failed")
    if core_floor * 2**25 <= 200_000:
        issues.append("terminal core h^2 separation failed")

    A = Fraction(64)
    phase_upper = A + 1 + A * A / (A - 1)
    if phase_upper >= 3 * A:
        issues.append("3sqrt(alpha) phase bound base case failed")
    if Fraction(11, 525) >= Fraction(1, 40):
        issues.append("sector-angle rational comparison failed")
    if Fraction(3199, 3200) <= sector_floor:
        issues.append("sector cosine reserve failed")


def phase(alpha: float, mode: int) -> float:
    return alpha * (math.log(alpha / mode) - 1) + mode


def independent_phase_numeric_audit(issues: list[str]) -> None:
    for alpha in (4096.37, 10_000.37, 100_000.37):
        lower = math.ceil(alpha * math.exp(-0.1))
        upper = math.floor(alpha)
        root = math.sqrt(alpha)
        total = 0j
        maximum_prefix = 0.0
        for mode in range(lower, upper + 1):
            total += cmath.exp(2j * math.pi * phase(alpha, mode))
            maximum_prefix = max(maximum_prefix, abs(total))
        if abs(total) >= 3 * root:
            issues.append(f"full phase-sum bound failed at alpha={alpha}")
        if maximum_prefix >= 3 * root:
            issues.append(f"phase-prefix bound failed at alpha={alpha}")

        H = math.ceil(root)
        outer_upper = upper - H
        if outer_upper >= lower:
            direct = sum(
                cmath.exp(2j * math.pi * phase(alpha, mode))
                for mode in range(lower, outer_upper + 1)
            )

            def delta(mode: int) -> float:
                return 1 - alpha * math.log1p(1 / mode)

            def multiplier(mode: int) -> complex:
                return 1 / (cmath.exp(2j * math.pi * delta(mode)) - 1)

            rebuilt = (
                multiplier(outer_upper)
                * cmath.exp(2j * math.pi * phase(alpha, outer_upper + 1))
                - multiplier(lower)
                * cmath.exp(2j * math.pi * phase(alpha, lower))
            )
            rebuilt += sum(
                (multiplier(mode - 1) - multiplier(mode))
                * cmath.exp(2j * math.pi * phase(alpha, mode))
                for mode in range(lower + 1, outer_upper + 1)
            )
            if abs(direct - rebuilt) > 2e-8 * max(1, abs(direct), abs(rebuilt)):
                issues.append(f"finite Abel identity failed at alpha={alpha}")
            outer_bound = (alpha - H + 1) / (H - 1)
            if abs(direct) > outer_bound + 1e-8:
                issues.append(f"finite outer phase bound failed at alpha={alpha}")


def independent_weighted_abel_audit(issues: list[str]) -> None:
    alpha = 10_000.37
    lower = math.ceil(alpha * math.exp(-0.1))
    upper = math.floor(alpha)
    modes = list(range(lower, upper + 1))
    phases = [cmath.exp(2j * math.pi * phase(alpha, mode)) for mode in modes]
    weights = [
        (1 + 0.05j * math.sin(mode / 17)) / mode
        for mode in modes
    ]
    direct = sum(weight * value for weight, value in zip(weights, phases))
    prefixes: list[complex] = []
    running = 0j
    for value in phases:
        running += value
        prefixes.append(running)
    rebuilt = weights[-1] * prefixes[-1]
    rebuilt += sum(
        (weights[index] - weights[index + 1]) * prefixes[index]
        for index in range(len(modes) - 1)
    )
    if abs(direct - rebuilt) > 2e-11:
        issues.append("weighted Abel identity failed")
    variation = abs(weights[-1]) + sum(
        abs(weights[index + 1] - weights[index])
        for index in range(len(weights) - 1)
    )
    if abs(direct) >= 3 * math.sqrt(alpha) * variation:
        issues.append("weighted Abel bound failed")


def endpoint_c(
    alpha: mp.mpf,
    mode: int,
    sigma: mp.mpf,
    t: mp.mpf,
) -> tuple[mp.mpf, mp.mpf]:
    s = mp.log(alpha / mode)
    v = mp.e ** (-s)
    z_abs = mp.sqrt(2 * (v - 1 + s))
    jacobian = v * z_abs / (1 - v)
    amplitude = mp.e ** (t * s**2 / 4 - sigma * s)
    formula = (amplitude - jacobian) / (mode * z_abs)
    y1 = -mp.sqrt(alpha) * z_abs
    b1 = mp.sqrt(alpha) * jacobian / mode
    b0 = mp.sqrt(alpha) * amplitude / mode
    quotient = (b1 - b0) / y1
    return formula, quotient


def independent_terminal_core_audit(issues: list[str]) -> None:
    mp.mp.dps = 90
    sigma = mp.mpf("0.501")
    t = mp.mpf(1) / 5000
    for alpha_text in ("4096.37", "10000.37", "100000.37"):
        alpha = mp.mpf(alpha_text)
        root = mp.sqrt(alpha)
        upper = int(mp.floor(alpha))
        K = int(mp.floor(root / 20))
        endpoint_phase = mp.e ** (2j * mp.pi * (
            alpha * (mp.log(alpha / upper) - 1) + upper
        ))
        weighted_sum = mp.mpc(0)
        maximum_spread = mp.mpf(0)
        for k in range(K + 1):
            mode = upper - k
            value, quotient = endpoint_c(alpha, mode, sigma, t)
            if not close(value, quotient, "1e-70"):
                issues.append(f"endpoint c quotient failed at alpha={alpha_text}, k={k}")
            if value < mp.mpf(13) / (100 * mode):
                issues.append(f"endpoint c floor failed at alpha={alpha_text}, k={k}")
            gauge = alpha * (mp.log(alpha / mode) - 1) + mode
            upper_gauge = alpha * (mp.log(alpha / upper) - 1) + upper
            spread = gauge - upper_gauge
            maximum_spread = max(maximum_spread, spread)
            weighted_sum += value * mp.e ** (2j * mp.pi * gauge)
        if maximum_spread >= mp.mpf(1) / 300:
            issues.append(f"terminal phase sector failed at alpha={alpha_text}")
        projected = mp.re(weighted_sum * mp.conj(endpoint_phase))
        lower_bound = mp.mpf(1287) / (200_000 * root)
        if projected <= lower_bound:
            issues.append(f"terminal c subblock floor failed at alpha={alpha_text}")


def independent_residual_audit(issues: list[str]) -> None:
    mp.mp.dps = 75
    imaginary = mp.mpc(0, 1)
    kappa = 1 / (2 * mp.pi * imaginary)
    lower = mp.mpf("-0.8")
    upper = mp.mpf("1.3")
    phase0 = mp.mpf("0.271")
    ephase = mp.e ** (2 * mp.pi * imaginary * phase0)

    def c(value: mp.mpf) -> mp.mpf:
        return 1 + value / 5 + value**2 / 7

    def c_prime(value: mp.mpf) -> mp.mpf:
        return mp.mpf(1) / 5 + 2 * value / 7

    def oscillation(value: mp.mpf) -> mp.mpc:
        return mp.e ** (-mp.pi * imaginary * value**2)

    direct = ephase * mp.quad(
        lambda value: value * c(value) * oscillation(value),
        [lower, 0, upper],
    )
    rebuilt = (
        kappa * ephase * c(lower) * oscillation(lower)
        - kappa * ephase * c(upper) * oscillation(upper)
        + kappa
        * ephase
        * mp.quad(
            lambda value: c_prime(value) * oscillation(value),
            [lower, 0, upper],
        )
    )
    if not close(direct, rebuilt, "1e-60"):
        issues.append("finite residual integration-by-parts identity failed")

    alpha = mp.mpf("113.37")
    B = 17
    for mode in (103, 109, 113):
        saddle = alpha * (mp.log(alpha / mode) - 1)
        v1 = mp.mpf(mode) / alpha
        h1 = v1 - 1 - mp.log(v1)
        vB = B * mp.mpf(mode) / alpha
        hB = vB - 1 - mp.log(vB)
        lower_phase = saddle - alpha * h1
        upper_phase = saddle - alpha * hB
        if not close(lower_phase, -mode, "1e-65"):
            issues.append(f"finite lower phase collapse failed at r={mode}")
        if not close(upper_phase, alpha * mp.log(B) - mode * B, "1e-65"):
            issues.append(f"finite upper phase collapse failed at r={mode}")
        if not close(mp.e ** (2j * mp.pi * lower_phase), 1, "1e-65"):
            issues.append(f"finite lower phase coherence failed at r={mode}")


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
        "full_grouped_roster_bounds",
        "endpoint_composed_bounds",
        "evaluated_full_physical_remainders",
        "signed_flow_bounds",
        "phi_b_bounds",
    ):
        if payload.get("counts", {}).get(key) != 0:
            issues.append(f"{key} was overpromoted")
    boundary = payload.get("proof_boundary", "")
    for marker in (
        "does not lower-bound the complete signed residual",
        "endpoint-composed h^2 estimate",
        "signed flow",
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
        "# Six-Moment Morse-Fresnel Endpoint-Coherence Abel Gate",
        "phase sum strictly less than 3sqrt(alpha)",
        "1287/[200000sqrt(alpha)]",
        "200000h^2",
        "phi_r(alpha/r)-y_1^2/2=-r",
        "common value e(alpha log B)",
        "0 endpoint-composed bounds",
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
        independent_phase_numeric_audit(issues)
        independent_weighted_abel_audit(issues)
        independent_terminal_core_audit(issues)
        independent_residual_audit(issues)
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
