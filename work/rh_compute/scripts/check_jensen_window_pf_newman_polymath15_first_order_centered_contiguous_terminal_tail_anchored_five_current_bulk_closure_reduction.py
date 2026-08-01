#!/usr/bin/env python3
"""Validate the terminal-tail-anchored five-current bulk closure."""

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

import jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_current_bulk_closure_reduction as closure  # noqa: E402


STEM = closure.STEM
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    f"ttfc_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "domain",
            "coefficient",
            "moments",
            "ladder",
            "value",
            "abel",
            "phi",
            "closure",
            "dyadic",
            "target",
            "guard",
            "handoff",
        ),
        start=1,
    )
]
EXPECTED_COUNTS = {
    "rows": 12,
    "normalized_coefficient_currents": 1,
    "complex_moments": 5,
    "moment_derivative_identities": 3,
    "centered_basis_transforms": 1,
    "centered_derivative_identities": 2,
    "five_current_algebraic_closures": 1,
    "multiplicative_rejoins": 1,
    "signed_bulk_targets": 1,
    "signed_bulk_bounds": 0,
    "bulk_aggregate_sign_closures": 0,
    "xi_level_current_theorems": 0,
}
SUMMARY = (
    "validated terminal-tail five-current bulk closure: 12 rows, "
    "5 complex currents, 3 derivative identities, "
    "1 algebraic closure, 0 signed bulk bounds"
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


def independent_moment_ladder(issues: list[str]) -> None:
    # Four atoms make this independent of the builder's three-atom audit.
    kappa, s_prime = sp.symbols("kappa s_prime")
    ell = sp.symbols("ell_1:5", real=True)
    z = sp.symbols("z_1:5")
    delta = sp.symbols("delta_1:5")
    zx = [
        (kappa - s_prime * ell[index] + delta[index]) * z[index]
        for index in range(4)
    ]

    def h(power: int) -> sp.Expr:
        return sum(ell[index] ** power * z[index] for index in range(4))

    def d(power: int) -> sp.Expr:
        return sum(
            ell[index] ** power * delta[index] * z[index]
            for index in range(4)
        )

    def hx(power: int) -> sp.Expr:
        return sum(ell[index] ** power * zx[index] for index in range(4))

    if sp.expand(hx(0) - (kappa * h(0) - s_prime * h(1) + d(0))) != 0:
        issues.append("independent H_0 ladder failed")
    if sp.expand(hx(1) - (kappa * h(1) - s_prime * h(2) + d(1))) != 0:
        issues.append("independent H_1 ladder failed")

    log_n = sp.symbols("log_N", real=True)
    r = log_n * h(0) - h(1)
    rx_direct = log_n * hx(0) - hx(1)
    rx_closed = (
        kappa * r
        - s_prime * (log_n * h(1) - h(2))
        + log_n * d(0)
        - d(1)
    )
    if sp.expand(rx_direct - rx_closed) != 0:
        issues.append("independent R_B closure failed")

    k0 = h(0)
    k1 = log_n * h(0) - h(1)
    k2 = log_n**2 * h(0) - 2 * log_n * h(1) + h(2)
    e0 = d(0)
    e1 = log_n * d(0) - d(1)
    kappa_n = kappa - s_prime * log_n
    if sp.expand(hx(0) - (kappa_n * k0 + s_prime * k1 + e0)) != 0:
        issues.append("independent centered K_0 ladder failed")
    if sp.expand(rx_direct - (kappa_n * k1 + s_prime * k2 + e1)) != 0:
        issues.append("independent centered K_1 ladder failed")


def independent_components(issues: list[str]) -> None:
    kr, ki, c, b = sp.symbols("kr ki c b", real=True)
    h0r, h0i, h1r, h1i, d0r, d0i = sp.symbols(
        "h0r h0i h1r h1i d0r d0i", real=True
    )
    value = (kr + sp.I * ki) * (h0r + sp.I * h0i)
    value -= (c + sp.I * b) * (h1r + sp.I * h1i)
    value += d0r + sp.I * d0i
    expected_r = kr * h0r - ki * h0i - c * h1r + b * h1i + d0r
    expected_i = kr * h0i + ki * h0r - c * h1i - b * h1r + d0i
    if sp.simplify(sp.re(sp.expand_complex(value)) - expected_r) != 0:
        issues.append("independent real component closure failed")
    if sp.simplify(sp.im(sp.expand_complex(value)) - expected_i) != 0:
        issues.append("independent imaginary component closure failed")


def validate(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("unexpected artifact kind")
    if payload.get("date") != "2026-07-31":
        issues.append("artifact date drifted")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")
    rows = payload.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or order drifted")
    if sum(row.get("readiness") == "open" for row in rows) != 2:
        issues.append("open-row count drifted")

    status = payload.get("status", "")
    for token in (
        "exact terminal-tail-anchored",
        "Phi_B bound open",
        "no bulk aggregate sign closure",
        "or RH",
    ):
        if token not in status:
            issues.append(f"status missing token: {token}")
    boundary = payload.get("proof_boundary", "")
    for token in (
        "no signed bound on Phi_B",
        "Xi-level current theorem",
        "Abel gap",
        "Lambda<=0",
        "RH",
    ):
        if token not in boundary:
            issues.append(f"proof boundary missing token: {token}")

    hashes = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in closure.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
        elif hashes.get(key) != file_hash(path):
            issues.append(f"source hash drifted: {key}")

    try:
        sources = closure.load_sources()
        if payload.get("source_audit") != closure.source_audit(sources):
            issues.append("source audit recomputation drifted")
        if payload.get("symbolic_certificate") != closure.symbolic_certificate():
            issues.append("symbolic certificate drifted")
        if payload.get("exact") != closure.exact_payload():
            issues.append("exact payload drifted")
    except Exception as exc:  # noqa: BLE001
        issues.append(f"certificate recomputation failed: {exc}")

    exact = payload.get("exact", {})
    if exact.get("sharp_target") != (
        "Prove Phi_B<=1/400, preferably Phi_B<=1/800, for the actual "
        "physical five-current vector and endpoint-terminal tail."
    ):
        issues.append("sharp target drifted")

    independent_moment_ladder(issues)
    independent_components(issues)


def validate_note(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for token in (
        "Terminal-Tail-Anchored Five-Current Bulk Closure",
        "Date: 2026-07-31",
        "0 signed bulk bounds",
        "H_j^(B)=sum_(n=1)^B",
        "H_(0,x)=kappa_0*H_0",
        "partial_x hat(R_B)=kappa_0*hat(R_B)",
        "lambda_n=log(N/n)",
        "partial_x K_1=chi_N*K_1+s_*'*K_2+E_1",
        "H_0^(B),H_1^(B),H_2^(B),D_0^(B),D_1^(B)",
        "Phi_B<=1/400",
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
