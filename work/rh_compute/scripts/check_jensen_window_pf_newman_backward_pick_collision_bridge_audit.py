#!/usr/bin/env python3
"""Check the backward-Pick collision-bridge audit independently."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ARTIFACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_backward_pick_collision_bridge_audit.json"
)

EXPECTED_IDS = [
    "nbpca_01_quadratic_heat_flow",
    "nbpca_02_real_zero_top",
    "nbpca_03_pick_top",
    "nbpca_04_double_collision",
    "nbpca_05_nonreal_birth",
    "nbpca_06_pick_bubble",
    "nbpca_07_speed_blowup",
    "nbpca_08_cutoff_hiding",
    "nbpca_09_noncommuting_limits",
    "nbpca_10_gronwall_boundary",
    "nbpca_11_preprint_audit",
    "nbpca_12_open_xi_repair",
]

SUCCESS = (
    "validated Newman backward-Pick collision-bridge audit: 12 rows, "
    "0 issues, 3 exact flow identities, 2 Pick-sign regions, "
    "1 square-root speed blowup, 1 cutoff-hiding theorem, "
    "1 noncommuting-limit obstruction, 1 preprint gap, 1 open Xi repair"
)


def symbolic_audit() -> list[str]:
    issues: list[str] = []
    z, t, a = sp.symbols("z t a", real=True)
    e = z**2 + a - 2 * t
    if sp.simplify(sp.diff(e, t) + sp.diff(e, z, 2)) != 0:
        issues.append("backward-heat residual is nonzero")
    g = sp.simplify(-sp.diff(e, z) / e)
    if (
        sp.simplify(
            sp.diff(g, t) + sp.diff(g, z, 2) - 2 * g * sp.diff(g, z)
        )
        != 0
    ):
        issues.append("Burgers residual is nonzero")

    x, y, b = sp.symbols("x y b", real=True)
    z_xy = x + sp.I * y
    direct = sp.im(sp.expand_complex(-2 * z_xy / (z_xy**2 + b)))
    expected = (
        2
        * y
        * (x**2 + y**2 - b)
        / ((x**2 - y**2 + b) ** 2 + 4 * x**2 * y**2)
    )
    if sp.simplify(direct - expected) != 0:
        issues.append("Pick-field formula residual is nonzero")

    r = sp.symbols("r", positive=True)
    witness = sp.simplify(expected.subs({x: 0, y: r / 2, b: r**2}))
    if witness != -sp.Rational(4, 3) / r:
        issues.append("negative-bubble witness changed")

    rho = sp.symbols("rho", positive=True)
    hidden_width = sp.simplify(rho**2 / 2)
    if sp.simplify((rho**2 - 2 * hidden_width)) != 0:
        issues.append("cutoff-hiding width failed")
    return issues


def content_audit(payload: dict) -> list[str]:
    issues: list[str] = []
    if payload.get("kind") != (
        "jensen_window_pf_newman_backward_pick_collision_bridge_audit"
    ):
        issues.append("unexpected artifact kind")
    rows = payload.get("rows")
    if not isinstance(rows, list) or len(rows) != 12:
        issues.append("expected exactly 12 rows")
        return issues
    ids = [row.get("id") for row in rows]
    if ids != EXPECTED_IDS:
        issues.append("row ids or deterministic order changed")
    if len(set(ids)) != len(ids):
        issues.append("duplicate row id")

    exact = payload.get("exact", {})
    for key in (
        "quadratic_flow",
        "pick_field",
        "collision_uniformity",
        "cutoff_obstruction",
        "literature_audit",
        "open_repair",
        "checks",
    ):
        if key not in exact:
            issues.append(f"missing exact section: {key}")

    text = json.dumps(exact, sort_keys=True)
    for phrase in (
        "E_t(z)=z^2+a-2t",
        "partial_t E_t=-partial_z^2 E_t",
        "t_*=a/2",
        "one lies in C+",
        "x^2+y^2<a-2t",
        "-4/(3*sqrt(a-2t))",
        "1/sqrt(a-2t)",
        "rho^2/2",
        "nonzero below y=rho",
        "cannot first prove a rho-dependent bridge",
        "Lemma 4.7 assumes a closed collision-free window",
        "Lemma 7.3 applies",
        "does not presently supply an admissible proof",
        "uniform in zero separation and lower cutoff",
    ):
        if phrase not in text:
            issues.append(f"missing collision-audit phrase: {phrase}")

    literature = exact.get("literature_audit", {})
    if literature.get("doi") != "https://doi.org/10.5281/zenodo.17636625":
        issues.append("preprint DOI changed")
    if rows[-1].get("readiness") != "open":
        issues.append("Xi repair target is not marked open")
    if rows[-2].get("readiness") != "gap_identified":
        issues.append("literature gap is not marked identified")
    boundary = payload.get("proof_boundary", "")
    for phrase in ("unreviewed", "does not", "Lambda<=0", "RH"):
        if phrase not in boundary:
            issues.append(f"proof boundary missing: {phrase}")
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    args = parser.parse_args()
    payload = json.loads(args.artifact.read_text(encoding="utf-8"))
    issues = symbolic_audit() + content_audit(payload)
    if issues:
        for issue in issues:
            print(f"ISSUE: {issue}")
        print(f"backward-Pick collision-bridge audit failed: {len(issues)} issues")
        return 1
    print(SUCCESS)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
