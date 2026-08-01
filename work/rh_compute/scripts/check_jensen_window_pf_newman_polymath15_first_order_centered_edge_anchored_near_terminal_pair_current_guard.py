#!/usr/bin/env python3
"""Validate the edge-anchored near-terminal pair-current guard."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402

import jensen_window_pf_newman_polymath15_first_order_centered_edge_anchored_near_terminal_pair_current_guard as gate
import jensen_window_pf_newman_polymath15_first_order_centered_real_edge_projective_current_gate as edge_gate


STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "edge_anchored_near_terminal_pair_current_guard"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    "eant_01_path",
    "eant_02_phase_ratio",
    "eant_03_phase_derivation",
    "eant_04_carrier_source",
    "eant_05_carrier_jets",
    "eant_06_edge_jets",
    "eant_07_polarization",
    "eant_08_positive",
    "eant_09_negative",
    "eant_10_reversal",
    "eant_11_route",
    "eant_12_boundary",
]
EXPECTED_COUNTS = {
    "rows": 12,
    "near_terminal_phase_limits": 1,
    "carrier_jet_limits": 4,
    "edge_jet_limits": 4,
    "exact_pair_polarizations": 1,
    "arb_point_witnesses": 2,
    "positive_pair_symbols": 1,
    "negative_pair_symbols": 1,
    "uniform_termwise_edge_absorptions": 0,
    "complete_aggregate_counterexamples": 0,
    "abel_gaps": 0,
    "contact_exclusions": 0,
}
SUMMARY = (
    "validated edge-anchored near-terminal pair-current guard: "
    "12 rows, 1 phase limit, 4 carrier jets, 4 edge jets, "
    "2 Arb witnesses, 1 asymptotic sign reversal, "
    "0 uniform termwise absorptions, 0 aggregate counterexamples"
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


def independent_symbolic_audit(issues: list[str]) -> None:
    h = sp.symbols("h", positive=True)
    theta = sp.symbols("theta", real=True)
    m = sp.symbols("m", integer=True, positive=True)
    expression = (
        2
        * sp.pi
        * (
            sp.log(1 - h * (theta + m))
            - sp.log(1 - h * theta)
        )
        / h**2
        + 2 * sp.pi * m * (1 / h - theta)
    )
    target = -sp.pi * (4 * m * theta + m**2)
    if sp.simplify(sp.limit(expression, h, 0, dir="+") - target) != 0:
        issues.append("independent phase limit failed")

    a, b, mm, x, y, k = sp.symbols("A B M X Y k", real=True)
    d = k * y / 2
    n = y / (16 * sp.pi) - k**2 * x / 4
    pair = (a + x) * (mm + n) - (b + d) ** 2
    expanded = (
        a * mm
        - b**2
        + x * n
        - d**2
        + a * n
        + x * mm
        - 2 * b * d
    )
    if sp.simplify(pair - expanded) != 0:
        issues.append("independent pair polarization failed")


def independent_arb_audit(issues: list[str]) -> None:
    flint.ctx.prec = gate.PRECISION_BITS
    checks = (
        (Fraction(-5, 8), 3, "positive", Fraction(1, 10)),
        (Fraction(0), 3, "negative", Fraction(-2)),
    )
    for point, m, direction, margin in checks:
        p = edge_gate.aq(point)
        pi = flint.arb.pi()
        theta = (1 - p) / 2
        phase = p * p / 2 - p + edge_gate.aq(Fraction(3, 8))
        q = flint.acb(0, -pi * phase).exp()
        relative = (-1) ** m * flint.acb(
            0, -4 * pi * m * theta
        ).exp()
        z = -q * relative
        a, a_p, a_pp = edge_gate.direct_a_jet(p)
        b = -a_p / (4 * pi)
        mm = a_pp / (16 * pi**2)
        k = m + theta
        x, y = z.real, z.imag
        d = k * y / 2
        n = y / (16 * pi) - k**2 * x / 4
        pair = (a + x) * (mm + n) - (b + d) ** 2
        bound = edge_gate.aq(margin)
        if direction == "positive" and not pair > bound:
            issues.append("independent positive pair witness failed")
        if direction == "negative" and not pair < bound:
            issues.append("independent negative pair witness failed")


def validate(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("unexpected artifact kind")
    if payload.get("date") != "2026-07-31":
        issues.append("artifact date drifted")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")
    if [row.get("id") for row in payload.get("rows", [])] != EXPECTED_IDS:
        issues.append("row ids or order drifted")

    for token in ("near-terminal", "sign-reversal", "Arb"):
        if token not in payload.get("status", ""):
            issues.append(f"status missing token: {token}")
    for token in (
        "termwise edge absorption only",
        "no complete Xi aggregate",
        "Abel gap",
        "contact exclusion",
        "Lambda<=0",
        "RH",
    ):
        if token not in payload.get("proof_boundary", ""):
            issues.append(f"proof boundary missing token: {token}")

    hashes = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in gate.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
        elif hashes.get(key) != file_hash(path):
            issues.append(f"source hash drifted: {key}")

    try:
        if payload.get("symbolic_certificate") != gate.symbolic_certificate():
            issues.append("stored symbolic certificate drifted")
    except Exception as exc:  # noqa: BLE001
        issues.append(f"symbolic recomputation failed: {exc}")

    interval = payload.get("interval_certificate", {})
    if interval.get("positive_margin") != "P_3(-5/8)>1/10":
        issues.append("positive margin drifted")
    if interval.get("negative_margin") != "P_3(0)<-2":
        issues.append("negative margin drifted")

    independent_symbolic_audit(issues)
    independent_arb_audit(issues)


def validate_note(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for token in (
        "Near-Terminal Pair-Current Guard",
        "Four Jets",
        "Pair Current",
        "P_3(-5/8)>1/10",
        "P_3(0)<-2",
        "rejects termwise edge absorption",
        "It does not reject the complete",
        "does not prove an aggregate cross-current bound",
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
