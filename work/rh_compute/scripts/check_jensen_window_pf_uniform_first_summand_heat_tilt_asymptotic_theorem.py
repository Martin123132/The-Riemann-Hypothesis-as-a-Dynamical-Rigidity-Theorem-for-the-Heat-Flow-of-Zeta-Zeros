#!/usr/bin/env python3
"""Validate the uniform first-summand heat-tilt asymptotic theorem."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "ufshta_01_integral_mapping",
    "ufshta_02_published_saddle_theorem",
    "ufshta_03_uniform_suitability",
    "ufshta_04_uniform_all_order_expansion",
    "ufshta_05_tilt_log_expansion",
    "ufshta_06_lambert_difference_bound",
    "ufshta_07_remainder_and_corrections",
    "ufshta_08_uniform_heat_tilt_theorem",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Uniform First-Summand Heat-Tilt Asymptotic Theorem",
    "uniform compact-heat first-summand heat-tilt asymptotic theorem",
    "Theorem 5.2",
    "f_T(t)=exp(-T*(log t)^2/16)",
    "w'(k)=w/(k*(1+w))",
    "w^(3R)/k^R=o(w/k^m)",
    "Delta_k^m log R_T^(1)(k)=O(log(k)/k^m)",
    "closes the sole",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class HeatTiltIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: object) -> HeatTiltIssue:
    return HeatTiltIssue(section=section, issue=name, detail=str(detail))


def validate_ref(ref: object) -> list[HeatTiltIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue("artifact", "bad-ref", repr(ref))]
    if ref.startswith("http"):
        return []
    if not (REPO_ROOT / ref).exists():
        return [issue("artifact", "missing-ref", ref)]
    return []


def validate_artifact(artifact: dict) -> list[HeatTiltIssue]:
    findings: list[HeatTiltIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem"
    ):
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("date") != "2026-07-13":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    for ref in artifact.get("sources", []):
        findings.extend(validate_ref(ref))
    for key in ("generator", "checker"):
        findings.extend(validate_ref(artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "first-summand",
        "does not by itself prove",
        "complete-kernel",
        "uniform order-four",
        "arbitrary-column",
        "rh",
        "lambda<=0",
    ):
        if required not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", required))

    expected_summary = {
        "rows": 8,
        "exact_or_published_rows": 8,
        "ready_to_apply_rows": 8,
        "open_analytic_rows": 0,
        "suitability_coefficients_checked": 8,
        "lambert_derivative_orders_checked": 7,
        "uniform_heat_tilt_theorems": 1,
    }
    summary = artifact.get("summary", {})
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            findings.append(issue("summary", f"bad-{key}", summary.get(key)))

    published = artifact.get("published_input", {})
    if published.get("url") != "https://arxiv.org/abs/2007.13582":
        findings.append(issue("published", "bad-url", published.get("url")))
    if published.get("theorem") != "Theorem 5.2":
        findings.append(issue("published", "bad-theorem", published.get("theorem")))

    suitability = artifact.get("uniform_suitability", {})
    coefficients = suitability.get("symbolic_coefficients_through_7", [])
    if len(coefficients) != 8:
        findings.append(issue("suitability", "bad-coefficient-count", len(coefficients)))
    for expected_degree, row in enumerate(coefficients):
        if row.get("degree") != expected_degree:
            findings.append(issue("suitability", "bad-degree", row))
        if row.get("T_degree", 99) > expected_degree or row.get("v_degree", 99) > expected_degree:
            findings.append(issue("suitability", "bad-polynomial-degree", row))

    lambert = artifact.get("lambert_derivatives", {})
    derivative_rows = lambert.get("rows", [])
    if [row.get("order") for row in derivative_rows] != list(range(1, 8)):
        findings.append(issue("lambert", "bad-orders", derivative_rows))
    if any(row.get("limit_w_to_infinity") in (None, "oo", "-oo") for row in derivative_rows):
        findings.append(issue("lambert", "bad-limit", derivative_rows))

    theorem = artifact.get("theorem", {})
    theorem_text = " ".join(str(value) for value in theorem.values())
    for required in (
        "I_pi(t^(5/4)*f_T(t);2k)",
        "-T*w^2/16",
        "O_R(w^(3R)/k^R)",
        "R>m",
        "Delta_k^m log R_T^(1)(k)=O(log(k)/k^m)",
        "backward-difference",
    ):
        if required not in theorem_text:
            findings.append(issue("theorem", "missing-formula", required))

    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        findings.append(issue("rows", "bad-row-ids", sorted(str(value) for value in ids)))
    if any(row.get("readiness") != "ready_to_apply" for row in rows):
        findings.append(issue("rows", "non-ready-row", rows))
    if any(not row.get("claim") or not row.get("proof_boundary") for row in rows):
        findings.append(issue("rows", "incomplete-row", rows))

    try:
        recomputed = build_artifact()
        for key in (
            "published_input",
            "uniform_suitability",
            "lambert_derivatives",
            "theorem",
            "rows",
            "summary",
        ):
            if artifact.get(key) != recomputed.get(key):
                findings.append(issue("recompute", f"mismatch-{key}", "stored differs"))
    except Exception as exc:
        findings.append(issue("recompute", "exception", exc))
    return findings


def validate_note(path: Path) -> list[HeatTiltIssue]:
    if not path.exists():
        return [issue("note", "missing", path)]
    text = path.read_text(encoding="utf-8")
    return [
        issue("note", "missing-text", required)
        for required in REQUIRED_NOTE_STRINGS
        if required not in text
    ]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    findings: list[HeatTiltIssue] = []
    if not args.input.exists():
        findings.append(issue("artifact", "missing", args.input))
    else:
        try:
            artifact = json.loads(args.input.read_text(encoding="utf-8"))
            findings.extend(validate_artifact(artifact))
        except Exception as exc:
            findings.append(issue("artifact", "exception", exc))
    findings.extend(validate_note(args.note))

    if args.json:
        print(json.dumps([asdict(value) for value in findings], indent=2))
    elif findings:
        print(f"uniform first-summand heat tilt: {len(findings)} issues")
        for finding in findings:
            print(f"- [{finding.section}] {finding.issue}: {finding.detail}")
    else:
        print(
            "validated uniform first-summand heat tilt: "
            "8 rows, 0 issues, 8 exact/published rows, 8 suitability coefficients, "
            "7 Lambert derivative orders, 0 open rows"
        )
    raise SystemExit(1 if findings else 0)


if __name__ == "__main__":
    main()
