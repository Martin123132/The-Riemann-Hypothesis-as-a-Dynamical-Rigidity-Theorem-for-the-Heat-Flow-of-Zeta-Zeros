#!/usr/bin/env python3
"""Validate the Q208 closed-boundary winding certificate."""

from __future__ import annotations

from fractions import Fraction
import json

import check_jensen_window_pf_newman_q208_bottom_phase_cell_certificate as bottom_check
import check_jensen_window_pf_newman_q208_top_phase_cell_certificate as top_check
import jensen_window_pf_newman_q208_closed_boundary_winding_certificate as cert


def validate() -> list[str]:
    issues: list[str] = []
    cert.compact.flint.ctx.prec = cert.bottom.bridge.PRECISION_BITS
    if not cert.DEFAULT_OUT.exists():
        return ["closed-boundary result is missing"]
    if not cert.DEFAULT_NOTE.exists():
        issues.append("closed-boundary note is missing")

    issues.extend(
        f"bottom source: {issue}" for issue in bottom_check.validate()
    )
    issues.extend(f"top source: {issue}" for issue in top_check.validate())
    artifact = json.loads(cert.DEFAULT_OUT.read_text(encoding="utf-8"))
    if artifact.get("kind") != cert.STEM:
        issues.append("closed-boundary kind drifted")
    try:
        rebuilt = cert.build_certificate()
    except (RuntimeError, ValueError, KeyError) as exc:
        issues.append(f"independent composition failed: {exc}")
        return issues
    if artifact != rebuilt:
        issues.append("stored closed-boundary payload drifted")

    expected_hashes = {
        "bottom": cert.file_hash(cert.BOTTOM_RESULT),
        "top": cert.file_hash(cert.TOP_RESULT),
        "right": cert.file_hash(cert.RIGHT_RESULT),
        "phase": cert.file_hash(cert.PHASE_RESULT),
        "winding": cert.file_hash(cert.WINDING_RESULT),
        "positive_boundary": cert.file_hash(cert.BOUNDARY_RESULT),
    }
    if artifact.get("source_sha256") != expected_hashes:
        issues.append("closed-boundary source hashes drifted")

    raw_witnesses = artifact.get("witnesses", [])
    try:
        witnesses = [
            (Fraction(point[0]), Fraction(point[1]))
            for point in raw_witnesses
        ]
    except (IndexError, ValueError, ZeroDivisionError) as exc:
        issues.append(f"closed witnesses are not exact rationals: {exc}")
        witnesses = []
    if len(witnesses) != artifact.get("total_witnesses"):
        issues.append("closed witness count drifted")
    if artifact.get("total_cells") != artifact.get("total_witnesses"):
        issues.append("cyclic cell/witness counts differ")
    if len(artifact.get("witness_cell_labels", [])) != artifact.get(
        "total_cells"
    ):
        issues.append("witness-cell label count drifted")
    if witnesses:
        try:
            independent_winding = cert.phase.polygon_winding(witnesses)
        except ValueError as exc:
            issues.append(f"exact witness polygon is invalid: {exc}")
        else:
            if independent_winding != 0:
                issues.append(
                    f"independent exact winding is {independent_winding}"
                )
            if independent_winding != artifact.get("exact_winding"):
                issues.append("stored exact winding drifted")

    counts = artifact.get("cell_counts", {})
    bottom_payload = json.loads(
        cert.BOTTOM_RESULT.read_text(encoding="utf-8")
    )
    top_payload = json.loads(cert.TOP_RESULT.read_text(encoding="utf-8"))
    expected_bottom = bottom_payload["summary"]["open_phase_chain"][
        "cell_count"
    ]
    expected_top = top_payload["summary"]["open_phase_chain"]["cell_count"]
    if counts != {
        "bottom": expected_bottom,
        "right": 20,
        "top_reverse": expected_top,
        "axis": 1,
    }:
        issues.append("edge cell counts drifted")
    if artifact.get("right_edge_audit", {}).get("cell_count") != 20:
        issues.append("right-edge audit cell count drifted")
    if cert.compact.arb(
        artifact.get("right_edge_audit", {}).get(
            "minimum_proxy_derivative_lower",
            "0",
        )
    ).lower() <= 0:
        issues.append("right-edge proxy derivative is not positive")

    if artifact.get("proxy", {}).get(
        "contact_and_winding_transfer"
    ) is not True:
        issues.append("proxy winding transfer flag is false")
    if artifact.get("domain", {}).get("q208") != (
        "[1/1040,1/4]x[0,246]"
    ):
        issues.append("Q208 domain drifted")
    if artifact.get("conclusion") != (
        "Q_208 is contact-free and has zero first-jet winding; "
        "therefore Q_1 through Q_208 are certified by containment."
    ):
        issues.append("finite Q208 conclusion drifted")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "finite rectangle",
        "Lambda<=1/5",
        "Q208 is contact-free",
        "no stage j>=209",
        "no cofinal termination theorem",
        "no Lambda<=0",
        "no RH",
        "no Clay prize",
    ):
        if marker not in boundary:
            issues.append(f"proof boundary missing: {marker}")
    if cert.DEFAULT_NOTE.exists():
        note = cert.DEFAULT_NOTE.read_text(encoding="utf-8")
        for marker in (
            "Exact Boundary",
            "exact witness-polygon winding=0",
            "Q_1 through Q_208",
            "no parameter-uniform termination theorem",
            "not a proof of `Lambda<=0`",
        ):
            if marker not in note:
                issues.append(f"closed-boundary note missing: {marker}")
    return issues


def main() -> int:
    issues = validate()
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    artifact = json.loads(cert.DEFAULT_OUT.read_text(encoding="utf-8"))
    print(
        "validated Q208 closed-boundary winding certificate: "
        f"{artifact['total_cells']} cells, "
        f"{artifact['total_witnesses']} exact witnesses, "
        "winding=0, Q1..Q208 certified, 0 cofinal promotions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
