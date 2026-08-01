#!/usr/bin/env python3
"""Independently validate the endpoint-relative phase-current gate."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_centered_endpoint_relative_phase_current_recurrence_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "c0_source_extraction": (
        REPO_ROOT
        / "work/rh_compute/scripts/"
        "jensen_window_pf_newman_polymath15_"
        "critical_RS_C1_endpoint_peeling_contract.py"
    ),
    "complex_endpoint_corrigendum": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_complex_endpoint_source_normalization_gate.json"
    ),
    "adjacent_saddle_recurrence": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_adjacent_saddle_recurrence.json"
    ),
    "adjacent_chart_stability": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_adjacent_chart_stability_certificate.json"
    ),
    "interior_projective_current": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_interior_projective_current_gate.json"
    ),
    "contact_transport": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_contact_signed_transport_reduction.json"
    ),
}
EXPECTED_IDS = [
    "erpc_01_endpoint_coordinates",
    "erpc_02_endpoint_current",
    "erpc_03_endpoint_phase_rate",
    "erpc_04_frame_current",
    "erpc_05_carrier_rate",
    "erpc_06_relative_current",
    "erpc_07_relative_rate",
    "erpc_08_near_terminal_limit",
    "erpc_09_midpoint_sign",
    "erpc_10_cutoff_sign",
    "erpc_11_terminal_sign_guard",
    "erpc_12_tail_recurrence",
    "erpc_13_current_defect",
    "erpc_14_projection_guard",
    "erpc_15_contact_minor_sum",
    "erpc_16_exceptional_fibre",
]
EXPECTED_COUNTS = {
    "rows": 16,
    "exact_endpoint_currents": 2,
    "division_free_relative_currents": 1,
    "asymptotic_sign_witnesses": 2,
    "terminal_tail_recurrences": 1,
    "contact_minor_sums": 1,
    "uniform_terminal_signs": 0,
    "abel_gaps": 0,
    "winding_bounds": 0,
}
SUMMARY = (
    "validated endpoint-relative phase-current recurrence gate: "
    "16 rows, 2 exact endpoint currents, 1 division-free relative current, "
    "2 asymptotic sign witnesses, 1 terminal-tail recurrence, "
    "1 contact-minor sum, 0 uniform terminal signs, 0 Abel gaps, "
    "0 winding bounds"
)


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path, label: str, issues: list[str]) -> dict:
    if not path.is_file():
        issues.append(f"missing {label}: {path}")
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"invalid {label}: {exc}")
        return {}


def validate_endpoint_current(issues: list[str]) -> None:
    k, kx, t0 = sp.symbols("k k_x T", real=True)
    hr, hi, xr, xi = sp.symbols("hr hi xr xi", real=True)
    h = hr + sp.I * hi
    hx = xr + sp.I * xi
    e = k * (t0 + sp.I) * h
    ex = kx * (t0 + sp.I) * h + k * h / 2 + k * (t0 + sp.I) * hx
    actual = sp.expand_complex(sp.im(ex * sp.conjugate(e)))
    expected = sp.expand_complex(
        k**2
        * (
            (t0**2 + 1) * sp.im(hx * sp.conjugate(h))
            - (hr**2 + hi**2) / 2
        )
    )
    if sp.simplify(actual - expected) != 0:
        issues.append("endpoint self-current identity failed")

    jr, ji, v, radial = sp.symbols("jr ji v radial", real=True)
    j = jr + sp.I * ji
    g = k * (t0 + sp.I) * j
    framed_ex = g + (radial + sp.I * v) * e
    framed = sp.expand_complex(
        sp.im(framed_ex * sp.conjugate(e))
    )
    framed_expected = sp.expand_complex(
        k**2
        * (t0**2 + 1)
        * (
            sp.im(j * sp.conjugate(h))
            + v * (hr**2 + hi**2)
        )
    )
    if sp.simplify(framed - framed_expected) != 0:
        issues.append("framed endpoint current failed")

    k_log = -sp.pi / 8 + 1 / (4 * t0) + 1 / (2 * (t0 + sp.I))
    if sp.simplify(
        sp.im(k_log) + 1 / (2 * (t0**2 + 1))
    ) != 0:
        issues.append("T_0+i phase-current term failed")


def validate_relative_current(issues: list[str]) -> None:
    k, t0, hr, hi, jr, ji = sp.symbols(
        "k T hr hi jr ji", real=True
    )
    v, b, u, delta, z2 = sp.symbols(
        "v b u delta z2", real=True
    )
    h = hr + sp.I * hi
    j = jr + sp.I * ji
    h2 = hr**2 + hi**2
    endpoint_mass = k**2 * (t0**2 + 1) * h2
    endpoint_current = (
        k**2
        * (t0**2 + 1)
        * (sp.im(j * sp.conjugate(h)) + v * h2)
    )
    carrier_current = (v + b * u + delta) * z2
    relative = sp.expand(
        endpoint_mass * carrier_current - z2 * endpoint_current
    )
    target = sp.expand_complex(
        k**2
        * (t0**2 + 1)
        * z2
        * (
            h2 * (b * u + delta)
            - sp.im(j * sp.conjugate(h))
        )
    )
    if sp.simplify(relative - target) != 0:
        issues.append("division-free relative current failed")
    if sp.simplify(target.subs({hr: 0, hi: 0})) != 0:
        issues.append("relative current did not vanish at H=0")
    if sp.simplify(target.subs(z2, 0)) != 0:
        issues.append("relative current did not vanish at z=0")


def validate_asymptotic_witnesses(issues: list[str]) -> None:
    p = sp.symbols("p", real=True)
    c0 = (
        sp.exp(sp.I * sp.pi * (p**2 / 2 + sp.Rational(3, 8)))
        - sp.I * sp.sqrt(2) * sp.cos(sp.pi * p / 2)
    ) / (2 * sp.cos(sp.pi * p))
    derivative = sp.diff(c0, p)
    if sp.simplify(derivative.subs(p, 0)) != 0:
        issues.append("C_0'(0) parity witness failed")

    phase_one = sp.simplify(
        sp.im(derivative.subs(p, 1) / c0.subs(p, 1))
        / sp.pi
    )
    radical = 1 - sp.sqrt(4 + 2 * sp.sqrt(2)) / 4
    if sp.simplify(phase_one - radical) != 0:
        issues.append("p=1 C_0 phase derivative failed")
    if radical.is_positive is not True:
        issues.append("p=1 terminal witness is not positive")

    midpoint_limit = -sp.Rational(1, 4)
    cutoff_limit = radical / 4
    if midpoint_limit.is_negative is not True:
        issues.append("p=0 terminal witness is not negative")
    if cutoff_limit.is_positive is not True:
        issues.append("p=1 terminal limit is not positive")

    m = sp.symbols("m", integer=True, nonnegative=True)
    theta = (1 - p) / 2
    generic = sp.Symbol("Phi", real=True) / (4 * sp.pi) - (m + theta) / 2
    stated = (
        sp.Symbol("Phi", real=True) / (4 * sp.pi)
        - m / 2
        - (1 - p) / 4
    )
    if sp.simplify(generic - stated) != 0:
        issues.append("near-terminal offset law failed")


def validate_recurrence(issues: list[str]) -> None:
    e0, e1, e2, z1, z2 = sp.symbols(
        "e0 e1 e2 z1 z2", complex=True
    )
    q0 = z1 + e1 - e0
    q1 = z2 + e2 - e1
    if sp.simplify(e2 + z2 + z1 - e0 - q0 - q1) != 0:
        issues.append("two-step tail recurrence failed")

    hr, hi, hxr, hxi = sp.symbols(
        "hr hi hxr hxi", real=True
    )
    rr, ri, rxr, rxi = sp.symbols(
        "rr ri rxr rxi", real=True
    )
    h = hr + sp.I * hi
    hx = hxr + sp.I * hxi
    r = rr + sp.I * ri
    rx = rxr + sp.I * rxi
    hp = r - h
    hpx = rx - hx
    actual = sp.expand_complex(sp.im(hpx * sp.conjugate(hp)))
    expected = sp.expand_complex(
        sp.im(rx * sp.conjugate(r))
        + sp.im(hx * sp.conjugate(h))
        - sp.im(rx * sp.conjugate(h) + hx * sp.conjugate(r))
    )
    if sp.simplify(actual - expected) != 0:
        issues.append("mixed current recurrence failed")

    y, sigma = sp.symbols("y sigma", real=True)
    q = sp.I * y
    qx = sp.I * sigma
    current = sp.expand_complex(
        sp.im(qx * sp.conjugate(1 + q))
    )
    if sp.simplify(current - sigma) != 0:
        issues.append("real-projection current witness failed")
    if sp.re(q) != 0 or sp.re(qx) != 0:
        issues.append("projection witness has nonzero real trace")


def validate_contact_sum(issues: list[str]) -> None:
    c0, d0, c1, c2, d1, d2 = sp.symbols(
        "c0 d0 c1 c2 d1 d2", real=True
    )
    minor_sum = c0 * d1 - d0 * c1 + c0 * d2 - d0 * c2
    value = c0 + c1 + c2
    slope = d0 + d1 + d2
    if sp.simplify(minor_sum - (c0 * slope - d0 * value)) != 0:
        issues.append("cumulative contact-minor identity failed")

    exceptional = {
        c0: 0,
        c1: 1,
        c2: -1,
        d0: 1,
        d1: 1,
        d2: 1,
    }
    if value.subs(exceptional) != 0:
        issues.append("exceptional witness is not on contact")
    if minor_sum.subs(exceptional) != 0:
        issues.append("exceptional minor sum did not vanish")
    if slope.subs(exceptional) == 0:
        issues.append("exceptional contact scalar unexpectedly vanished")


def validate_payload(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("unexpected artifact kind")
    status = payload.get("status", "")
    for token in (
        "phase-current",
        "terminal sign-reversal",
        "recurrence-defect",
    ):
        if token not in status:
            issues.append(f"status missing token: {token}")

    rows = payload.get("rows", [])
    row_ids = [row.get("id") for row in rows if isinstance(row, dict)]
    if row_ids != EXPECTED_IDS:
        issues.append("gate row ids or order drifted")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("gate counts drifted")

    proof_boundary = payload.get("proof_boundary", "")
    for token in (
        "no terminal-composed real-edge sign",
        "no signed cumulative-minor estimate",
        "no Abel-scalar gap",
        "no successor winding cap",
        "Lambda<=0",
        "RH",
    ):
        if token not in proof_boundary:
            issues.append(f"proof boundary missing token: {token}")

    source_audit = payload.get("source_audit", {})
    stored_hashes = source_audit.get("source_sha256", {})
    for key, path in SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
            continue
        if stored_hashes.get(key) != file_hash(path):
            issues.append(f"source hash drifted: {key}")

    exact = payload.get("exact", {})
    required_paths = (
        ("endpoint_self_current", "division_free"),
        ("endpoint_carrier_relative_current", "division_free"),
        ("near_terminal_asymptotic", "scaled_limit"),
        ("near_terminal_asymptotic", "midpoint_witness"),
        ("near_terminal_asymptotic", "cutoff_witness"),
        ("preprojection_tail_recurrence", "tail_value"),
        ("preprojection_tail_recurrence", "current_defect"),
        ("contact_minor_sum", "identity"),
        ("contact_minor_sum", "exceptional_blind_spot"),
        ("route_decision", "next_target"),
    )
    for section, field in required_paths:
        if not exact.get(section, {}).get(field):
            issues.append(f"missing exact field: {section}.{field}")


def validate_note(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for token in (
        "# Jensen-Window PF Newman Polymath-15 Endpoint-Relative",
        "## Endpoint Self-Current",
        "## Relative Carrier Current",
        "## Exact Sign-Reversal Guard",
        "## Tail Composition Before Projection",
        "## Contact-Minor Sum",
        "## Route Decision",
        "uniform terminal relative-current sign",
        "zero-real-projection",
        "This is not an Abel gap",
    ):
        if token not in text:
            issues.append(f"note missing token: {token}")
    fences = sum(1 for line in text.splitlines() if line.startswith("```"))
    if fences % 2:
        issues.append("note has unbalanced code fences")


def main() -> int:
    issues: list[str] = []
    payload = load_json(RESULT, "result", issues)
    if payload:
        validate_payload(payload, issues)
    validate_endpoint_current(issues)
    validate_relative_current(issues)
    validate_asymptotic_witnesses(issues)
    validate_recurrence(issues)
    validate_contact_sum(issues)
    validate_note(issues)

    if issues:
        print("endpoint-relative phase-current recurrence gate failed:")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
