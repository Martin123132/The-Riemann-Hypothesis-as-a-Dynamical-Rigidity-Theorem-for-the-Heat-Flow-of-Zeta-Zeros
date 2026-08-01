#!/usr/bin/env python3
"""Validate the adjacent-cutoff real-edge projective splice theorem."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_adjacent_cutoff_real_edge_projective_splice_gate as gate


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "adjacent_cutoff_real_edge_projective_splice_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    "acres_01_domain",
    "acres_02_alignment",
    "acres_03_endpoint_values",
    "acres_04_leading_wedge",
    "acres_05_finite_jets",
    "acres_06_wedge_expansion",
    "acres_07_error_budget",
    "acres_08_signed_wedge",
    "acres_09_nonvanishing",
    "acres_10_affine_current",
    "acres_11_q_ge_1",
    "acres_12_boundary",
]
EXPECTED_COUNTS = {
    "rows": 12,
    "exact_cutoff_endpoint_values": 4,
    "exact_leading_wedges": 1,
    "finite_connector_jet_families": 2,
    "rational_remainder_budgets": 1,
    "nonvanishing_affine_joins": 1,
    "signed_affine_splices": 1,
    "uniform_q_ge_1_model_splices": 1,
    "required_second_adjacent_x_jets": 0,
    "xi_level_edge_splices": 0,
    "abel_gaps": 0,
    "winding_bounds": 0,
}
SUMMARY = (
    "validated adjacent-cutoff real-edge projective splice gate: "
    "12 rows, 4 exact endpoint values, 1 exact leading wedge, "
    "1 rational remainder budget, 1 nonvanishing affine join, "
    "1 signed q>=1 model splice, 0 second adjacent x-jets, "
    "0 Xi-level splices"
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
    p = sp.symbols("p", real=True)
    pi = sp.pi
    theta = (1 - p) / 2
    f = sp.cos(pi * (p**2 / 2 + sp.Rational(3, 8))) / (
        2 * sp.cos(pi * p)
    )
    q_real = sp.cos(
        pi * (p**2 / 2 - p + sp.Rational(3, 8))
    )
    a_source = sp.simplify(f - q_real)
    b_source = sp.simplify(
        -sp.diff(f, p) / (4 * pi)
        + theta
        * sp.sin(pi * (p**2 / 2 - p + sp.Rational(3, 8)))
        / 2
    )
    am = sp.simplify(a_source.subs(p, -1))
    ap = sp.simplify(a_source.subs(p, 1))
    bm = sp.simplify(b_source.subs(p, -1))
    bp = sp.simplify(b_source.subs(p, 1))
    expected = (
        -sp.sqrt(2 + sp.sqrt(2)) / 4,
        -sp.sqrt(2 + sp.sqrt(2)) / 4,
        -3 * sp.sqrt(2 - sp.sqrt(2)) / 16,
        -sp.sqrt(2 - sp.sqrt(2)) / 16,
    )
    for name, value, target in zip(
        ("A_minus", "A_plus", "B_minus", "B_plus"),
        (am, ap, bm, bp),
        expected,
        strict=True,
    ):
        if sp.simplify(value - target) != 0:
            issues.append(f"independent endpoint identity failed: {name}")
    if sp.simplify(am * bp - bm * ap + sp.sqrt(2) / 32) != 0:
        issues.append("independent leading wedge failed")

    h, s = sp.symbols("h s", positive=True, real=True)
    e0m, e0p, e1m, e1p = sp.symbols(
        "e0m e0p e1m e1p", real=True
    )
    cm, cp = am + h * e0m, ap + h * e0p
    dm, dp = h * bm + h**2 * e1m, h * bp + h**2 * e1p
    wedge = sp.expand(cm * dp - dm * cp)
    cs, ds = (1 - s) * cm + s * cp, (1 - s) * dm + s * dp
    affine = sp.expand(cs * sp.diff(ds, s) - ds * sp.diff(cs, s))
    if sp.simplify(affine - wedge) != 0:
        issues.append("independent affine current failed")


def validate(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("unexpected artifact kind")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")
    if [row.get("id") for row in payload.get("rows", [])] != EXPECTED_IDS:
        issues.append("row ids or order drifted")

    for token in ("finite-height", "signed adjacent-cutoff", "first-order"):
        if token not in payload.get("status", ""):
            issues.append(f"status missing token: {token}")
    for token in (
        "higher-order Xi/source",
        "Xi-level edge sign",
        "cumulative-minor estimate",
        "Abel-scalar gap",
        "successor winding cap",
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

    independent_symbolic_audit(issues)
    try:
        if payload.get("symbolic_certificate") != gate.symbolic_certificate():
            issues.append("stored symbolic certificate drifted")
    except Exception as exc:  # noqa: BLE001
        issues.append(f"symbolic certificate recomputation failed: {exc}")

    arithmetic = payload.get("arithmetic_certificate", {})
    h0 = Fraction(1, gate.H_DENOMINATOR)
    error = 500 * h0 + 50_000 * h0**2
    if not error < Fraction(1, 600):
        issues.append("independent error cap failed")
    if not 100 * h0 < Fraction(1, 12):
        issues.append("independent nonvanishing cap failed")
    expected_arithmetic = {
        "h_upper": f"1/{gate.H_DENOMINATOR}",
        "endpoint_bounds": (
            "A<-1/3, |A|<1/2, |B_minus|<1, |B_plus|<1"
        ),
        "defect_bounds": "|e0_minus|,|e0_plus|<100; |e1_minus|,|e1_plus|<250",
        "second_order_coefficient": 500,
        "third_order_coefficient": 50_000,
        "error_at_h_upper": f"{error.numerator}/{error.denominator}",
        "error_cap": "1/600",
        "radical_margin": "sqrt(2)/32>1/24",
        "normalized_wedge": "Delta_cut/h<-1/25",
        "signed_wedge": "Delta_cut<-h/25<0",
        "nonvanishing": "c_minus<-1/4, c_plus<-1/4, c_s<-1/4",
    }
    if arithmetic != expected_arithmetic:
        issues.append("arithmetic certificate drifted")

    exact = payload.get("exact", {})
    for token in ("p_minus=-1", "p_plus=1", "S_plus", "-S_minus"):
        if token not in json.dumps(exact):
            issues.append(f"cutoff alignment token missing: {token}")


def validate_note(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for token in (
        "Cutoff Charts",
        "Exact Endpoint Wedge",
        "Finite-Height Budget",
        "Signed Join",
        "-sqrt(2)/32",
        "Delta_cut<-h/25<0",
        "c_s<-1/4",
        "second adjacent x-jet is not needed",
        "No Xi-level edge sign",
        "not an Xi-level",
    ):
        if token not in text:
            issues.append(f"note marker missing: {token}")


def main() -> int:
    issues: list[str] = []
    payload = load_json(RESULT, issues)
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
