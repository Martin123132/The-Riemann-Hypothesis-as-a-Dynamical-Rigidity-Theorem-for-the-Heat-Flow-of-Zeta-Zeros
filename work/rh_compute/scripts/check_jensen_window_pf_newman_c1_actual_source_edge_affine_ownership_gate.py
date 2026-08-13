#!/usr/bin/env python3
"""Independently check the actual-source edge-affine ownership gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_edge_affine_ownership_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(payload: dict) -> None:
    require(set(payload["source_sha256"]) == set(payload["source_audit"]), "source keys")
    require(len(payload["source_audit"]) == 6, "source count")
    for key, item in payload["source_audit"].items():
        path = REPO_ROOT / item["path"]
        require(path.is_file(), f"missing source {key}")
        digest = file_hash(path)
        require(digest == item["sha256"], f"source hash drift {key}")
        require(digest == payload["source_sha256"][key], f"source map drift {key}")


def check_matrix_ownership() -> int:
    J = sp.Rational(1, 2) * sp.Matrix(
        [
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, -1],
            [0, 0, -1, 0],
        ]
    )
    edge = sp.Matrix(sp.symbols("v_e n_e a_e q_e", real=True))
    carrier = sp.Matrix(sp.symbols("v_c n_c a_c q_c", real=True))
    error_lift = sp.Matrix(sp.symbols("dv dn da dq", real=True))

    t_edge = 2 * J * edge
    t_carrier = 2 * J * carrier
    t_terminal = 2 * J * (edge + carrier)
    require(t_edge == sp.Matrix([edge[1], edge[0], -edge[3], -edge[2]]), "edge row")
    require(sp.simplify(t_terminal - t_edge - t_carrier) == sp.zeros(4, 1), "row split")
    require(
        sp.simplify((t_terminal - t_edge).dot(error_lift) - 2 * (carrier.T * J * error_lift)[0]) == 0,
        "duplicate term",
    )
    require(sp.simplify(t_edge.dot(error_lift) - 2 * (edge.T * J * error_lift)[0]) == 0, "edge affine")
    return 4


def check_ideal_polynomial() -> int:
    x, u_x = sp.symbols("x u_x", real=True)
    v_e, n_e, a_e, q_e = sp.symbols("v_e n_e a_e q_e", real=True)
    c = 1
    a = sp.I * x / 2
    q = sp.I * x / 2
    n = -x**2 / 4 - sp.I * u_x / 2
    direct = sp.expand(sp.I * x * (n_e * c + v_e * n - q_e * a - a_e * q))
    expected = -sp.I * v_e * x**3 / 4 + (a_e + q_e) * x**2 / 2 + (v_e * u_x / 2 + sp.I * n_e) * x
    require(sp.simplify(direct - expected) == 0, "cubic identity")
    require(sp.Poly(direct, x).degree() == 3, "degree")
    require(sp.Poly(direct, x).coeff_monomial(x**3) == -sp.I * v_e / 4, "leading coefficient")
    return 3


def check_source_normalization() -> int:
    gamma, scale, o_r, o_i = sp.symbols("gamma scale o_r o_i", real=True, positive=True)
    omega = sp.cos(gamma) + sp.I * sp.sin(gamma)
    outer = o_r + sp.I * o_i
    nu = omega / scale
    direct = 2 * sp.re(nu * outer) / (sp.re(nu) ** 2 + sp.im(nu) ** 2)
    expected = 2 * scale * sp.re(omega * outer)
    require(sp.simplify(sp.expand_complex(direct - expected)) == 0, "normalization identity")
    require(sp.simplify(1 / sp.conjugate(nu) - scale * omega) == 0, "inverse source")
    return 2


def check_midpoint_and_budget() -> int:
    require(Fraction(2) < Fraction(23, 16) ** 2, "sqrt two rational upper bound")
    require(Fraction(2) < Fraction(6) ** 2, "midpoint imaginary sign guard")

    midpoint_real = sp.sqrt(2 - sp.sqrt(2)) / 4
    midpoint_imag = sp.sqrt(2 + sp.sqrt(2)) / 4 - sp.sqrt(2) / 2
    require(float(sp.N(midpoint_real, 80)) > 3 / 16, "midpoint real check")
    require(float(sp.N(midpoint_imag, 80)) < 0, "midpoint imaginary check")

    closed = Fraction(1, 6) + Fraction(1, 17160) + Fraction(1, 2000)
    reserve = Fraction(49, 50) - closed
    require(closed == Fraction(143479, 858000), "closed budget")
    require(reserve == Fraction(697361, 858000), "reserve")
    return 6


def check_rows(payload: dict) -> None:
    suffixes = (
        "split",
        "owner",
        "duplicate",
        "polynomial",
        "single",
        "source",
        "support",
        "terminal",
        "cubic",
        "midpoint",
        "class",
        "closed",
        "reclass",
        "target",
        "live",
    )
    expected = [f"eao_{index:02d}_{suffix}" for index, suffix in enumerate(suffixes, start=1)]
    rows = payload["rows"]
    require([row["id"] for row in rows] == expected, "row ids")
    require(sum(row["readiness"] == "open" for row in rows) == 1, "open rows")
    require(rows[-1]["readiness"] == "open", "frontier row")


def main() -> None:
    require(RESULT_PATH.is_file(), "missing result")
    require(NOTE_PATH.is_file(), "missing note")
    payload = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    require(payload["kind"] == KIND, "kind")
    require(payload["schema_version"] == 1, "schema")
    check_sources(payload)
    matrix_audits = check_matrix_ownership()
    polynomial_audits = check_ideal_polynomial()
    source_audits = check_source_normalization()
    rational_audits = check_midpoint_and_budget()
    check_rows(payload)

    summary = payload["summary"]
    require(summary["rows"] == 15, "summary rows")
    require(summary["edge_only_ownership_corrections"] == 1, "ownership summary")
    require(summary["exact_affine_polynomials"] == 1, "polynomial summary")
    require(summary["physical_cubic_nonsuppression_witnesses"] == 1, "witness summary")
    require(summary["absolute_affine_budgets"] == 0, "budget boundary")
    require(summary["signed_main_reclassifications"] == 1, "classification summary")
    require(summary["closed_secondary_fraction"] == "143479/858000", "closed fraction")
    require(summary["remaining_reserve_fraction"] == "697361/858000", "reserve fraction")
    require(summary["remaining_absolute_secondary_packages"] == 2, "remaining secondary")
    require(summary["live_signed_finite_affine_targets"] == 1, "live target")

    note = NOTE_PATH.read_text(encoding="utf-8")
    require(payload["success"] in note, "note success")
    require("P_E^0" in note, "cubic note")
    require("E_FA=E_F+E_aff" in note, "signed package note")
    require("not a proof of RH" in note, "proof boundary note")

    print(
        "validated Newman C1 actual-source edge-affine ownership gate: "
        f"15 rows, 0 issues, {matrix_audits} matrix audits, "
        f"{polynomial_audits} polynomial audits, {source_audits} source audits, "
        f"{rational_audits} rational audits, 1 edge-only ownership correction, "
        "0 absolute affine budgets, 1 signed-main reclassification, "
        "2 remaining absolute secondary packages, 1 live finite-affine target"
    )


if __name__ == "__main__":
    main()
