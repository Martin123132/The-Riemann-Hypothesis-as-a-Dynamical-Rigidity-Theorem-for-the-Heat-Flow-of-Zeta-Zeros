#!/usr/bin/env python3
"""Validate the terminal-tail-anchored Abel cross-current reduction."""

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

import jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_abel_cross_current_reduction as reduction  # noqa: E402


STEM = reduction.STEM
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    f"ttac_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "domain",
            "partition",
            "abel",
            "abel_x",
            "normalizer",
            "tail",
            "bulk_scalar",
            "bulk_x",
            "coordinates",
            "polarization",
            "remainder",
            "reserve",
            "guard",
            "handoff",
        ),
        start=1,
    )
]
EXPECTED_COUNTS = {
    "rows": 14,
    "exact_partitions": 1,
    "exact_abel_identities": 2,
    "normalized_coordinate_maps": 2,
    "current_polarization_identities": 2,
    "terminal_current_inputs": 1,
    "signed_bulk_targets": 1,
    "bulk_aggregate_closures": 0,
    "xi_level_current_theorems": 0,
}
SUMMARY = (
    "validated terminal-tail anchored Abel cross-current reduction: "
    "14 rows, 2 Abel identities, 2 current identities, "
    "1 signed Phi_B target, 0 bulk closures"
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


def independent_abel(issues: list[str]) -> None:
    # Use four generic atoms, independently of the builder's three-atom audit.
    kappa, h_1, h_2, h_3 = sp.symbols(
        "kappa h_1 h_2 h_3", real=True
    )
    z = sp.symbols("z_1:5")
    dz = sp.symbols("dz_1:5")
    weights = (
        kappa + h_1 + h_2 + h_3,
        kappa + h_2 + h_3,
        kappa + h_3,
        kappa,
    )
    direct = sum(weights[index] * z[index] for index in range(4))
    abel = kappa * sum(z)
    abel += h_1 * z[0]
    abel += h_2 * (z[0] + z[1])
    abel += h_3 * (z[0] + z[1] + z[2])
    if sp.expand(direct - abel) != 0:
        issues.append("independent four-atom Abel identity failed")

    direct_x = sum(weights[index] * dz[index] for index in range(4))
    abel_x = kappa * sum(dz)
    abel_x += h_1 * dz[0]
    abel_x += h_2 * (dz[0] + dz[1])
    abel_x += h_3 * (dz[0] + dz[1] + dz[2])
    if sp.expand(direct_x - abel_x) != 0:
        issues.append("independent differentiated Abel identity failed")


def independent_centering(issues: list[str]) -> None:
    c, b, u = sp.symbols("c b u", real=True)
    tr, ti, gr, fr, fi, rr, ri = sp.symbols(
        "tr ti gr fr fi rr ri", real=True
    )
    direct = gr + c * (u * fr + rr) - b * (u * fi + ri)
    direct -= c * u * (tr + fr)
    tail = gr - c * u * tr
    bulk = c * rr - b * ri - b * u * fi
    if sp.expand(direct - tail - bulk) != 0:
        issues.append("independent tail-anchored scalar failed")

    alternate = gr - c * u * tr + b * u * ti
    alternate -= b * u * (ti + fi)
    alternate += c * rr - b * ri
    if sp.expand(alternate - direct) != 0:
        issues.append("independent alternate re-anchoring failed")

    cx, bx, ux, fix, rrx, rix = sp.symbols(
        "cx bx ux fix rrx rix", real=True
    )
    derivative = (
        cx * rr
        + c * rrx
        - bx * ri
        - b * rix
        - (bx * u + b * ux) * fi
        - b * u * fix
    )
    expanded = (
        cx * rr
        + c * rrx
        - bx * ri
        - b * rix
        - bx * u * fi
        - b * ux * fi
        - b * u * fix
    )
    if sp.expand(derivative - expanded) != 0:
        issues.append("independent bulk scalar derivative failed")


def independent_current(issues: list[str]) -> None:
    ct, dt, et, nt, cb, db, eb, nb = sp.symbols(
        "ct dt et nt cb db eb nb", real=True
    )
    total = (ct + cb) * (nt + nb) - (dt + db) * (et + eb)
    tail = ct * nt - dt * et
    bulk = cb * nb - db * eb
    cross = ct * nb + cb * nt - dt * eb - db * et
    if sp.expand(total - tail - bulk - cross) != 0:
        issues.append("independent four-coordinate polarization failed")

    h = sp.symbols("h", positive=True, real=True)
    xt, at, xtx, atx, xb, ab, xbx, abx = sp.symbols(
        "xt at xtx atx xb ab xbx abx", real=True
    )
    p_total = (
        (xt + xb) * (atx + abx) - (at + ab) * (xtx + xbx)
    ) / h**2
    p_tail = (xt * atx - at * xtx) / h**2
    stated = (
        (xt + xb) * abx
        + xb * atx
        - (at + ab) * xbx
        - ab * xtx
    ) / h**2
    if sp.simplify(p_total - p_tail - stated) != 0:
        issues.append("independent signed Phi_B identity failed")

    c0, a0, c0x, a0x, sigma, sigmax = sp.symbols(
        "c0 a0 c0x a0x sigma sigmax", real=True, nonzero=True
    )
    cn = c0 / sigma
    an = a0 / sigma
    cnx = (c0x * sigma - c0 * sigmax) / sigma**2
    anx = (a0x * sigma - a0 * sigmax) / sigma**2
    lhs = (cn * anx - an * cnx) / h**2
    rhs = (c0 * a0x - a0 * c0x) / (h**2 * sigma**2)
    if sp.simplify(lhs - rhs) != 0:
        issues.append("independent normalizer cancellation failed")


def independent_guard(issues: list[str]) -> None:
    tail = (sp.Integer(1), sp.Integer(0), sp.Integer(0), sp.Integer(-1))
    bulk = (
        sp.Rational(-1, 2),
        sp.Integer(0),
        sp.Integer(0),
        sp.Integer(2),
    )
    p_tail = tail[0] * tail[3] - tail[1] * tail[2]
    p_bulk = bulk[0] * bulk[3] - bulk[1] * bulk[2]
    p_total = (tail[0] + bulk[0]) * (tail[3] + bulk[3])
    p_total -= (tail[1] + bulk[1]) * (tail[2] + bulk[2])
    p_cross = p_total - p_tail - p_bulk
    if (p_tail, p_bulk, p_cross, p_total) != (
        -1,
        -1,
        sp.Rational(5, 2),
        sp.Rational(1, 2),
    ):
        issues.append("independent negative-components guard failed")


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
    if sum(row.get("readiness") == "open" for row in rows) != 1:
        issues.append("open-row count drifted")

    status = payload.get("status", "")
    for token in ("exact q=1", "Phi_B target open", "no bulk closure", "or RH"):
        if token not in status:
            issues.append(f"status missing token: {token}")
    boundary = payload.get("proof_boundary", "")
    for token in (
        "no bulk aggregate closure",
        "Xi-level current theorem",
        "Abel gap",
        "Lambda<=0",
        "RH",
    ):
        if token not in boundary:
            issues.append(f"proof boundary missing token: {token}")

    hashes = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in reduction.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
        elif hashes.get(key) != file_hash(path):
            issues.append(f"source hash drifted: {key}")

    try:
        sources = reduction.load_sources()
        if payload.get("source_audit") != reduction.source_audit(sources):
            issues.append("source audit recomputation drifted")
        if payload.get("symbolic_certificate") != reduction.symbolic_certificate():
            issues.append("symbolic certificate drifted")
        if payload.get("exact") != reduction.exact_payload():
            issues.append("exact payload drifted")
    except Exception as exc:  # noqa: BLE001
        issues.append(f"certificate recomputation failed: {exc}")

    exact = payload.get("exact", {})
    if exact.get("sharp_join_target") != (
        "The signed inequality Phi_B<=1/400 implies P_ret<0. "
        "The stronger Phi_B<=1/800 preserves P_ret<-1/800."
    ):
        issues.append("sharp join target drifted")

    independent_abel(issues)
    independent_centering(issues)
    independent_current(issues)
    independent_guard(issues)


def validate_note(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for token in (
        "Terminal-Tail Anchored Abel Cross-Current Reduction",
        "Date: 2026-07-31",
        "B=N-M-1",
        "R_(B,x)=kappa_B*F_(B,x)",
        "P_ret=P_T+P_B+P_cross",
        "Phi_B=P_B+P_cross",
        "Phi_B<=1/400 implies P_ret<0",
        "O(B) time and O(1) streaming state",
        "0 bulk closures",
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
