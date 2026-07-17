#!/usr/bin/env python3
"""Validate the deep-Schur Toda and zero-boundary obstruction gate."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_deep_schur_toda_boundary_gate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    SOURCE_PATHS,
    SYMBOLIC_TODA_ORDERS,
    boundary_reset_audit,
    build_artifact,
    sha256,
    strict_schur_audit,
    strict_schur_jensen_counterexample,
    symbolic_toda_audit,
)


@dataclass(frozen=True)
class GateIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> GateIssue:
    return GateIssue(section, code, str(detail))


def parse_fraction(record: object) -> Fraction:
    if not isinstance(record, dict):
        raise TypeError("fraction record is not an object")
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def validate_artifact(path: Path) -> tuple[dict, list[GateIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[GateIssue] = []

    if artifact.get("kind") != "jensen_window_pf_deep_schur_toda_boundary_gate":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    status = str(artifact.get("status", ""))
    for marker in (
        "Toda",
        "zero-boundary",
        "strict-Schur/Jensen",
        "rejected",
        "open",
    ):
        if marker not in status:
            issues.append(issue("artifact", "bad-status", marker))
    boundary_text = str(artifact.get("proof_boundary", ""))
    for marker in (
        "m>=10",
        "reject",
        "Xi/Phi Jensen bridge",
        "RH",
        "Lambda<=0",
    ):
        if marker not in boundary_text:
            issues.append(issue("artifact", "bad-proof-boundary", marker))

    expected_summary = {
        "rows": 15,
        "ready_to_apply_rows": 13,
        "open_endpoint_rows": 0,
        "rejected_endpoint_rows": 1,
        "separate_bridge_rows": 1,
        "symbolic_toda_checks": 4,
        "boundary_translation_checks": 251,
        "deep_translation_checks": 31,
        "shallow_boundary_mismatches": 220,
        "strict_schur_checks": 138,
        "exact_strict_schur_jensen_counterexamples": 1,
        "rejected_moving_tail_pf_routes": 1,
        "rejected_generic_schur_jensen_routes": 1,
        "direct_literature_closing_theorems": 0,
        "negative_deep_rectangles": 4,
        "open_rectangle_hierarchies": 0,
        "rejected_rectangle_hierarchies": 1,
        "separate_xi_specific_bridges": 1,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    expected_ids = {
        f"jwpfdstb_{index:02d}_{suffix}"
        for index, suffix in (
            (1, "input_coordinate"),
            (2, "tau_coordinate"),
            (3, "toda_identity"),
            (4, "log_concavity_criterion"),
            (5, "no_cluster_promotion"),
            (6, "ordinary_tail"),
            (7, "conditional_translation"),
            (8, "boundary_counterexample"),
            (9, "no_lr_tail_lift"),
            (10, "strict_schur_model"),
            (11, "all_shape_strictness"),
            (12, "jensen_failure"),
            (13, "literature_fit"),
            (14, "rejected_rectangles"),
            (15, "xi_specific_bridge"),
        )
    }
    rows = artifact.get("rows", [])
    row_map = {row.get("id"): row for row in rows}
    if set(row_map) != expected_ids or len(rows) != len(expected_ids):
        issues.append(issue("rows", "bad-id-set", sorted(str(key) for key in row_map)))
    if row_map.get("jwpfdstb_14_rejected_rectangles", {}).get("readiness") != "rejected_by_counterexample":
        issues.append(issue("rows", "rectangle-target-not-rejected", row_map.get("jwpfdstb_14_rejected_rectangles")))
    if row_map.get("jwpfdstb_15_xi_specific_bridge", {}).get("readiness") != "separate_open_obligation":
        issues.append(issue("rows", "bridge-not-separated", row_map.get("jwpfdstb_15_xi_specific_bridge")))

    expected_exact = {
        "toda_identity": (
            "tau_(m+1,N)*tau_(m-1,N)=tau_(m,N)^2-"
            "tau_(m,N-1)*tau_(m,N+1)"
        ),
        "ordinary_tail": (
            "b_k^(s)=h_(s+k)/h_s for k>=0 and b_k^(s)=0 for k<0"
        ),
        "boundary_defect": (
            "s_(0,0)(b^(1))-h_1^(-2)*s_(1,1)(h)=h_0*h_2/h_1^2>0"
        ),
        "rectangle_target": (
            "s_((N^m))(h)>0 for every m>=10 and N>=m-1"
        ),
        "deep_failure": (
            "s_((N^10))(h)<0 for N=9,10,11,12"
        ),
        "rejected_rectangle_target": (
            "s_((N^m))(h)>0 for every m>=10 and N>=m-1 is false"
        ),
    }
    exact = artifact.get("exact", {})
    for key, expected in expected_exact.items():
        if exact.get(key) != expected:
            issues.append(issue("exact", f"bad-{key}", exact.get(key)))

    source_rows = artifact.get("sources", [])
    source_map = {row.get("path"): row for row in source_rows}
    if len(source_rows) != len(SOURCE_PATHS):
        issues.append(issue("sources", "bad-count", len(source_rows)))
    for source_path in SOURCE_PATHS:
        relative = source_path.relative_to(source_path.parents[3]).as_posix()
        row = source_map.get(relative)
        if row is None:
            issues.append(issue("sources", "missing-source", relative))
        elif row.get("sha256") != sha256(source_path):
            issues.append(issue("sources", "hash-mismatch", relative))

    toda = artifact.get("symbolic_toda_audit", {})
    if toda.get("orders") != list(SYMBOLIC_TODA_ORDERS):
        issues.append(issue("toda", "bad-orders", toda.get("orders")))
    if toda.get("all_residuals_zero") is not True:
        issues.append(issue("toda", "stored-residual-failure", toda))
    if any(row.get("residual") != "0" for row in toda.get("rows", [])):
        issues.append(issue("toda", "nonzero-stored-residual", toda.get("rows")))

    boundary = artifact.get("boundary_reset_audit", {})
    witness = boundary.get("order_two_witness", {})
    try:
        defect = parse_fraction(witness.get("defect"))
        tail = parse_fraction(witness.get("tail_schur"))
        translated = parse_fraction(witness.get("translated_deep_schur"))
    except Exception as exc:
        issues.append(issue("boundary", "bad-fraction-record", exc))
    else:
        if defect <= 0 or tail - translated != defect:
            issues.append(issue("boundary", "bad-reset-defect", witness))
    if witness.get("symbolic_defect") != "h_0*h_2/h_1^2":
        issues.append(issue("boundary", "bad-symbolic-defect", witness))

    strict = artifact.get("strict_schur_audit", {})
    if strict.get("all_values_at_least_plancherel_term") is not True:
        issues.append(issue("strict-schur", "stored-positivity-failure", strict))
    try:
        weakest_ratio = parse_fraction(strict.get("weakest_ratio"))
    except Exception as exc:
        issues.append(issue("strict-schur", "bad-ratio", exc))
    else:
        if weakest_ratio < 1:
            issues.append(issue("strict-schur", "ratio-below-one", weakest_ratio))

    counterexample = artifact.get("strict_schur_jensen_counterexample", {})
    if counterexample.get("generating_function") != "exp(z/100)/((1-z)*(1-2*z))":
        issues.append(issue("counterexample", "bad-generating-function", counterexample))
    if counterexample.get("primitive_polynomial_coefficients_ascending") != [
        6000000,
        54180000,
        126540900,
        90420901,
    ]:
        issues.append(issue("counterexample", "bad-polynomial", counterexample))
    try:
        discriminant = parse_fraction(counterexample.get("discriminant"))
    except Exception as exc:
        issues.append(issue("counterexample", "bad-discriminant", exc))
    else:
        expected_discriminant = Fraction(-222484532394597, 2000000000000)
        if discriminant != expected_discriminant or discriminant >= 0:
            issues.append(issue("counterexample", "wrong-discriminant", discriminant))
    if counterexample.get("strictly_negative_discriminant") is not True:
        issues.append(issue("counterexample", "sign-not-certified", counterexample))

    literature = artifact.get("primary_literature_audit", [])
    expected_literature = {
        "lam_postnikov_pylyavskyy_schur_log_concavity": "formal_identity_support_only",
        "wagner_hadamard_products": "hypotheses_not_met",
        "angarone_kim_oh_soskin_dual_jacobi_trudi": "restricted_shape_only",
        "craven_csordas_fox_wright": "fixed_multiplier_family_support",
    }
    literature_map = {row.get("id"): row for row in literature}
    if set(literature_map) != set(expected_literature):
        issues.append(issue("literature", "bad-id-set", sorted(literature_map)))
    for row_id, classification in expected_literature.items():
        row = literature_map.get(row_id, {})
        if row.get("classification") != classification:
            issues.append(issue("literature", "bad-classification", row_id))
        if not str(row.get("url", "")).startswith("https://"):
            issues.append(issue("literature", "bad-url", row_id))

    direct_audits = (
        ("toda", symbolic_toda_audit),
        ("boundary", boundary_reset_audit),
        ("strict-schur", strict_schur_audit),
        ("counterexample", strict_schur_jensen_counterexample),
    )
    stored_keys = (
        "symbolic_toda_audit",
        "boundary_reset_audit",
        "strict_schur_audit",
        "strict_schur_jensen_counterexample",
    )
    for (section, audit), key in zip(direct_audits, stored_keys):
        try:
            rebuilt = audit()
        except Exception as exc:
            issues.append(issue(section, "direct-audit-exception", exc))
        else:
            if rebuilt != artifact.get(key):
                issues.append(issue(section, "direct-audit-drift", key))

    try:
        canonical = build_artifact()
    except Exception as exc:
        issues.append(issue("rebuild", "exception", exc))
    else:
        if canonical != artifact:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[GateIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: exact rectangular Toda coordinate",
        "ordinary subtraction-free cluster",
        "moving-tail PF equivalence is false",
        "h_0*h_2/h_1^2>0",
        "strictly Schur-positive for every partition",
        "discriminant=-222484532394597/2000000000000<0",
        "general conjecture remains open",
        "additional Xi/Phi structure",
        "proposed all-order endpoint hierarchy",
        "width-log-concavity input is false at order ten",
        "s_((N^10))(h)<0 for N=9,10,11,12",
    )
    issues: list[GateIssue] = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "rh is proved",
        "the endpoint rectangles are proved",
        "the jensen bridge is proved",
        "all relevant literature",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-overclaim", forbidden))
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    if issues:
        for finding in issues:
            print(
                f"DEEP-SCHUR-TODA-BOUNDARY {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"deep-Schur Toda/boundary gate: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated deep-Schur Toda/boundary gate: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['symbolic_toda_checks']} symbolic Toda checks, "
        f"{summary['boundary_translation_checks']} boundary checks, "
        f"{summary['strict_schur_checks']} strict-Schur checks, "
        "1 exact full-PF Jensen counterexample, "
        "4 negative order-ten rectangles, "
        "1 rejected rectangle hierarchy, 1 Xi-specific bridge"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
