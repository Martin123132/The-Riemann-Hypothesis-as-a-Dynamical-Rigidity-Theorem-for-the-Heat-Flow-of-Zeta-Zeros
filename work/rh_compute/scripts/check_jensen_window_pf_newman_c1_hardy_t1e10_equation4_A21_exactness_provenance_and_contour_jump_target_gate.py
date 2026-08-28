#!/usr/bin/env python3
"""Independently check the A21 provenance and split-contour target gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
for candidate in (REPO_ROOT / "work/rh_compute/vendor", SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

from flint import acb, arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation4_A21_exactness_provenance_and_contour_jump_target_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fresh_partial(q_text: str, sign: int, terms: int) -> tuple[acb, acb, acb]:
    q_value = arb(q_text)
    w = acb(arb.pi() * q_value / arb(2).sqrt(), -arb.pi() * q_value / arb(2).sqrt())
    secant = 1 / w.cos()
    # Reverse summation order to alter the arithmetic path.
    partial = acb(0)
    for k in range(terms - 1, -1, -1):
        alpha = 2 * k + 1
        partial += 2 * (-1 if k % 2 else 1) * (acb(0, sign * alpha) * w).exp()
    ratio = -(acb(0, 2 * sign) * w).exp()
    return secant, secant - partial, secant * ratio**terms


def main() -> int:
    require(gate.RESULT.is_file(), "missing A21 provenance result")
    require(gate.NOTE.is_file(), "missing A21 provenance note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == gate.STEM, "kind drift")
    require(
        artifact["status"]
        == "A21_heuristic_continuation_quarantined_exact_split_contour_boundary_layer_and_DK_target_certified",
        "status drift",
    )
    require(artifact.get("passed") is True, "A21 provenance gate failed")

    ctx.dps = 150
    ctx.threads = 1
    rows = artifact["half_contour_identities"]["high_precision_rows"]
    require(len(rows) == len(gate.SAMPLE_Q), "saved q roster drift")
    for saved in rows:
        sign = 1 if saved["valid_expansion"] == "A13" else -1
        secant, residual, closed = fresh_partial(saved["q"], sign, saved["terms"])
        require(residual.overlaps(closed), f"fresh residual identity failed at q={saved['q']}")
        require(secant.real.overlaps(arb(saved["secant_real_ball"])), "saved secant real miss")
        require(secant.imag.overlaps(arb(saved["secant_imag_ball"])), "saved secant imag miss")
        require(residual.real.overlaps(arb(saved["residual_real_ball"])), "saved residual real miss")
        require(residual.imag.overlaps(arb(saved["residual_imag_ball"])), "saved residual imag miss")

    for q_text, sign in (("-0.2", 1), ("-0.7", 1), ("0.2", -1), ("0.7", -1)):
        _, residual, closed = fresh_partial(q_text, sign, 19)
        require(residual.overlaps(closed), f"altered q/N residual identity failed at q={q_text}")

    pi = arb.pi()
    root_two = arb(2).sqrt()
    expected = (-root_two * pi).exp()
    layer = artifact["nonuniform_boundary_layer_certificate"]
    require(expected.overlaps(arb(layer["expected_boundary_factor_ball"])), "boundary factor drift")
    for terms in (11, 23, 47):
        q_value = -arb(1) / terms
        w = acb(pi * q_value / root_two, -pi * q_value / root_two)
        ratio = -(acb(0, 2) * w).exp()
        require((abs(ratio) ** terms).overlaps(expected), "altered boundary-layer scaling failed")

    q_pos = arb("0.3")
    w_pos = acb(pi * q_pos / root_two, -pi * q_pos / root_two)
    q_neg = -q_pos
    w_neg = acb(pi * q_neg / root_two, -pi * q_neg / root_two)
    require(abs(-(acb(0, 2) * w_pos).exp()).lower() > 1, "A13 invalid-side term test lost")
    require(abs(-(acb(0, -2) * w_neg).exp()).lower() > 1, "A18 invalid-side term test lost")

    provenance = artifact["source_provenance_audit"]
    require("heuristically" in provenance["Lewis_2015_pdf_pages_1_based"]["38"], "A21 heuristic guard missing")
    require("suitable convergence" in provenance["Lewis_2015_pdf_pages_1_based"]["41"], "A33 conditional guard missing")
    require("delegated" in provenance["Lewis_Brereton_2026_pdf_pages_1_based"]["3"], "2026 provenance link missing")

    decision = artifact["decision"]
    require(decision["finite_QK_definition_exact"] is True, "finite Q_K was discarded")
    require(decision["finite_source_to_QK_transform_identity_retained"] is True, "finite transform identity was discarded")
    require(decision["equation4_exactness_imported_from_2015"] is False, "equation-(4) exactness overpromoted")
    require(decision["A21_permitted_as_exact_whole_contour_interchange"] is False, "A21 overpromoted")
    require(decision["A33_permitted_as_unconditional_identity"] is False, "A33 overpromoted")
    require(decision["exact_split_contour_route_selected"] is True, "split-contour route lost")
    require(decision["C_K_enclosed"] is False, "C_K overclaim")
    require(decision["D_K_enclosed"] is False, "D_K overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")

    target = artifact["contour_defect_target"]
    exact_h_path = REPO_ROOT / artifact["dependencies"]["exact_H_defect_target"]["path"]
    exact_h = json.loads(exact_h_path.read_text(encoding="utf-8"))["continuation_defect_target"]
    require(target["D_K_safe_lower_ball"] == exact_h["D_K_safe_lower_ball"], "D_K lower target drift")
    require(target["D_K_safe_upper_ball"] == exact_h["D_K_safe_upper_ball"], "D_K upper target drift")
    require(target["boundary_layer_alone_asserted_to_equal_D_K"] is False, "boundary layer overpromoted")

    for section in ("dependencies", "references", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file {path}")
            require(file_hash(path) == record["sha256"], f"hash drift {path}")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in (
        "provenance gap isolated",
        "not currently licensed",
        "not uniform",
        "q=O(1/N)",
        "not from an inserted geometric circle constant",
        "does not derive the contour defect",
    ):
        require(token in note, f"note boundary token missing: {token}")

    print(
        "independently checked A21 nonpromotion and split-contour boundary layer; "
        f"factor={expected.str(30, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
