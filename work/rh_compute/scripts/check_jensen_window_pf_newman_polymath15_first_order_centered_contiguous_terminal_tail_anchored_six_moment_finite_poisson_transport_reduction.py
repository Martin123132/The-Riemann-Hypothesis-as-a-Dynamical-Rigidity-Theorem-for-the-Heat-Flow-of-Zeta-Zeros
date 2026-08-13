#!/usr/bin/env python3
"""Validate the six-moment finite-Poisson and transport reduction."""

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

import jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_finite_poisson_transport_reduction as gate  # noqa: E402


STEM = gate.STEM
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    f"fptr_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "domain",
            "source",
            "moments",
            "amplitude",
            "endpoint",
            "poisson",
            "shift",
            "phase",
            "saddle",
            "signature",
            "transport",
            "residual",
            "observations",
            "external",
            "dual",
            "upper",
            "lower",
            "c1",
            "dyadic",
            "global",
            "coarse",
            "handoff",
            "boundary",
        ),
        start=1,
    )
]
EXPECTED_COUNTS = {
    "rows": 23,
    "bare_logarithmic_moments": 6,
    "full_poisson_identities": 6,
    "explicit_primal_half_endpoints": 2,
    "stationary_main_families": 6,
    "stationary_transport_identities": 5,
    "residual_transport_recurrences": 5,
    "physical_flow_observations": 8,
    "required_value_transforms": 6,
    "dual_cutoff_jump_laws": 2,
    "dyadic_b_process_hypothesis_audits": 1,
    "imported_inexplicit_local_b_process_bounds": 1,
    "explicit_uniform_remainder_constants": 0,
    "smooth_endpoint_transition_bounds": 0,
    "signed_flow_bounds": 0,
    "phi_b_bounds": 0,
}
SUMMARY = (
    "validated six-moment finite-Poisson/transport reduction: 23 rows, "
    "6 moments, 5 transport recurrences, 2 cutoff jump laws, "
    "0 explicit uniform remainder constants, 0 signed flow bounds"
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


def independent_source_audit(payload: dict, issues: list[str]) -> None:
    saved = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in gate.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
        elif saved.get(key) != file_hash(path):
            issues.append(f"source hash drifted: {key}")
    source = payload.get("source_audit", {}).get("primary_analytic_source", {})
    if source.get("url") != gate.VANDEHEY_URL:
        issues.append("primary analytic source URL drifted")
    if "full starred Poisson formula (paper equation (10))" not in source.get(
        "specialized_items", []
    ):
        issues.append("primary Poisson source marker missing")


def independent_amplitude_audit(issues: list[str]) -> None:
    logarithm, t, sigma = sp.symbols("L t s", real=True)
    for order in range(6):
        amplitude = logarithm**order * sp.exp(
            t * logarithm**2 / 4 - sigma * logarithm
        )
        previous = (
            0
            if order == 0
            else logarithm ** (order - 1)
            * sp.exp(t * logarithm**2 / 4 - sigma * logarithm)
        )
        expected = order * previous + (
            t * logarithm / 2 - sigma
        ) * amplitude
        require_zero(
            sp.diff(amplitude, logarithm) - expected,
            f"independent amplitude recurrence failed at {order}",
            issues,
        )


def independent_endpoint_audit(issues: list[str]) -> None:
    alpha, b, t, sigma = sp.symbols("alpha B t s", positive=True)
    logarithm = sp.log(b)
    phase = sp.exp(2 * sp.pi * sp.I * alpha * logarithm)
    for order in range(5):
        lower = 1 if order == 0 else 0
        amplitude = logarithm**order * sp.exp(
            t * logarithm**2 / 4 - sigma * logarithm
        )
        next_amplitude = logarithm ** (order + 1) * sp.exp(
            t * logarithm**2 / 4 - sigma * logarithm
        )
        endpoint = (lower + amplitude * phase) / 2
        endpoint_next = next_amplitude * phase / 2
        require_zero(
            sp.diff(endpoint, alpha) / (2 * sp.pi)
            - sp.I * endpoint_next,
            f"independent endpoint shift failed at {order}",
            issues,
        )


def independent_stationary_audit(issues: list[str]) -> None:
    alpha, u, mode = sp.symbols("alpha u r", positive=True)
    phase = alpha * sp.log(u) - mode * u
    saddle = alpha / mode
    require_zero(
        sp.diff(phase, u).subs(u, saddle),
        "independent saddle failed",
        issues,
    )
    curvature = sp.diff(phase, u, 2).subs(u, saddle)
    prefactor = sp.sqrt(alpha) / mode
    require_zero(
        prefactor**2 * (-curvature) - 1,
        "independent stationary prefactor failed",
        issues,
    )
    expected_phase = alpha * (sp.log(alpha / mode) - 1)
    require_zero(
        sp.expand_log(phase.subs(u, saddle), force=True)
        - sp.expand_log(expected_phase, force=True),
        "independent stationary phase failed",
        issues,
    )


def independent_transport_audit(issues: list[str]) -> None:
    alpha, mode, t, sigma = sp.symbols(
        "alpha r t s", positive=True
    )
    logarithm = sp.log(alpha / mode)
    prefactor = sp.sqrt(alpha) / mode
    phase = sp.exp(
        2
        * sp.pi
        * sp.I
        * (alpha * (logarithm - 1) - sp.Rational(1, 8))
    )
    for order in range(5):
        amplitude = logarithm**order * sp.exp(
            t * logarithm**2 / 4 - sigma * logarithm
        )
        next_amplitude = logarithm ** (order + 1) * sp.exp(
            t * logarithm**2 / 4 - sigma * logarithm
        )
        previous = (
            0
            if order == 0
            else logarithm ** (order - 1)
            * sp.exp(t * logarithm**2 / 4 - sigma * logarithm)
        )
        main = prefactor * amplitude * phase
        next_main = prefactor * next_amplitude * phase
        correction = (
            prefactor
            * phase
            / (2 * sp.pi * alpha)
            * (
                order * previous
                + (sp.Rational(1, 2) - sigma + t * logarithm / 2)
                * amplitude
            )
        )
        require_zero(
            sp.diff(main, alpha) / (2 * sp.pi)
            - sp.I * next_main
            - correction,
            f"independent transport failed at {order}",
            issues,
        )


def finite_transport_audit(issues: list[str]) -> None:
    mp.mp.dps = 80
    alpha = mp.mpf(113) / 7
    mode = mp.mpf(5)
    t = mp.mpf(3) / 101
    sigma = mp.mpf(7) / 13

    def main(order: int, value: mp.mpf) -> mp.mpc:
        logarithm = mp.log(value / mode)
        amplitude = logarithm**order * mp.e ** (
            t * logarithm**2 / 4 - sigma * logarithm
        )
        return (
            mp.sqrt(value)
            / mode
            * amplitude
            * mp.e
            ** (
                2
                * mp.pi
                * mp.j
                * (value * (logarithm - 1) - mp.mpf(1) / 8)
            )
        )

    logarithm = mp.log(alpha / mode)
    for order in range(5):
        amplitude = logarithm**order * mp.e ** (
            t * logarithm**2 / 4 - sigma * logarithm
        )
        previous = (
            mp.mpf(0)
            if order == 0
            else logarithm ** (order - 1)
            * mp.e ** (t * logarithm**2 / 4 - sigma * logarithm)
        )
        phase = mp.e ** (
            2
            * mp.pi
            * mp.j
            * (alpha * (logarithm - 1) - mp.mpf(1) / 8)
        )
        correction = (
            mp.sqrt(alpha)
            / mode
            * phase
            / (2 * mp.pi * alpha)
            * (
                order * previous
                + (mp.mpf(1) / 2 - sigma + t * logarithm / 2)
                * amplitude
            )
        )
        numerical = mp.diff(lambda value: main(order, value), alpha) / (2 * mp.pi)
        expected = mp.j * main(order + 1, alpha) + correction
        if abs(numerical - expected) > mp.mpf("1e-65"):
            issues.append(f"finite transport audit failed at {order}")


def independent_dyadic_audit(issues: list[str]) -> None:
    alpha, u, block, b = sp.symbols("alpha u U B", positive=True)
    second = sp.diff(alpha * sp.log(u), u, 2)
    third = sp.diff(alpha * sp.log(u), u, 3)
    fourth = sp.diff(alpha * sp.log(u), u, 4)
    if second != -alpha / u**2:
        issues.append("independent second derivative failed")
    if third != 2 * alpha / u**3:
        issues.append("independent third derivative failed")
    if fourth != -6 * alpha / u**4:
        issues.append("independent fourth derivative failed")
    global_ratio = sp.simplify(
        sp.Abs(second.subs(u, 1)) / sp.Abs(second.subs(u, b))
    )
    dyadic_ratio = sp.simplify(
        sp.Abs(second.subs(u, block))
        / sp.Abs(second.subs(u, 2 * block))
    )
    if global_ratio != b**2:
        issues.append("independent global curvature ratio failed")
    if dyadic_ratio != 4:
        issues.append("independent dyadic curvature ratio failed")


def structural_audit(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("kind drifted")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")
    rows = payload.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids drifted")
    if len({row.get("id") for row in rows}) != len(rows):
        issues.append("duplicate row ids")
    for key in (
        "explicit_uniform_remainder_constants",
        "smooth_endpoint_transition_bounds",
        "signed_flow_bounds",
        "phi_b_bounds",
    ):
        if payload.get("counts", {}).get(key) != 0:
            issues.append(f"{key} was overpromoted")
    boundary = payload.get("proof_boundary", "")
    for marker in (
        "no explicit uniform Poisson/B-process remainder constant",
        "imaginary Hermitian coefficient bound",
        "signed flow estimate",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")


def note_audit(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for marker in (
        "# Six-Moment Finite-Poisson And Stationary-Transport Reduction",
        gate.VANDEHEY_URL,
        "H_j=E_j+lim_(R->infinity)",
        "e(-1/8)=exp(-i*pi/4)",
        "M_(j,r)'(xi)=iM_(j+1,r)+C_(j,r)",
        "upper crossings occur at xi=2*pi*q",
        "0 explicit uniform remainder constants",
        "this is not a proof",
        f"python work/rh_compute/scripts/check_{STEM}.py",
    ):
        if marker not in text:
            issues.append(f"note marker missing: {marker}")


def main() -> int:
    issues: list[str] = []
    payload = load_json(RESULT, issues)
    if payload:
        independent_source_audit(payload, issues)
        independent_amplitude_audit(issues)
        independent_endpoint_audit(issues)
        independent_stationary_audit(issues)
        independent_transport_audit(issues)
        finite_transport_audit(issues)
        independent_dyadic_audit(issues)
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
