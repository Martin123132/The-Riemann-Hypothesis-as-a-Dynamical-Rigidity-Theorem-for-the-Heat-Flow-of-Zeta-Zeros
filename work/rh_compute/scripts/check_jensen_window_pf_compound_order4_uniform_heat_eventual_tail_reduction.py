#!/usr/bin/env python3
"""Validate the order-four uniform-heat eventual-tail reduction."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order4_uniform_heat_eventual_tail_reduction import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "co4uhet_01_uniform_ratio_contract",
    "co4uhet_02_symbolic_cancellation",
    "co4uhet_03_uniform_determinant_transfer",
    "co4uhet_04_finite_confinement_composition",
    "co4uhet_05_heat_tilt_reduction",
    "co4uhet_06_first_summand_input",
    "co4uhet_07_uniform_heat_tilt_target",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Compound Order-Four Uniform-Heat Eventual-Tail Reduction",
    "exact uniform-asymptotic order-four reduction",
    "[h^0,...,h^6] det K=[0,0,0,0,0,0,768*G_2^6]",
    "Delta_k^m log R_T^(1)(k)=O(log(k)/k^m)",
    "Since `G_2->1` uniformly",
    "finite-confinement theorem",
    "only remaining",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class UniformIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: object) -> UniformIssue:
    return UniformIssue(section=section, issue=name, detail=str(detail))


def validate_ref(ref: object) -> list[UniformIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue("artifact", "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue("artifact", "missing-ref", ref)]
    return []


def validate_artifact(artifact: dict) -> list[UniformIssue]:
    findings: list[UniformIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order4_uniform_heat_eventual_tail_reduction"
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
        "does not prove",
        "xi-specific",
        "uniform eventual tail",
        "unconditional",
        "arbitrary-column",
        "rh",
        "lambda<=0",
    ):
        if required not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", required))

    expected_summary = {
        "rows": 7,
        "exact_or_input_rows": 6,
        "ready_to_apply_rows": 4,
        "conditional_rows": 2,
        "open_analytic_rows": 1,
        "symbolic_coefficients_checked": 7,
        "uniform_determinant_transfers": 1,
        "finite_confinement_compositions": 1,
        "open_heat_tilt_theorems": 1,
        "higher_theta_handoffs_closed": 1,
    }
    summary = artifact.get("summary", {})
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            findings.append(issue("summary", f"bad-{key}", summary.get(key)))

    cancellation = artifact.get("symbolic_cancellation", {})
    if cancellation.get("coefficients_h0_through_h6") != [
        "0",
        "0",
        "0",
        "0",
        "0",
        "0",
        "768*G2**6",
    ]:
        findings.append(issue("symbolic", "bad-coefficients", cancellation))
    if cancellation.get("permutations_checked") != 24:
        findings.append(issue("symbolic", "bad-permutation-count", cancellation))

    exact = artifact.get("exact", {})
    exact_text = " ".join(str(value) for value in exact.values())
    for required in (
        "r_(lambda,M,j)",
        "o(h(lambda,M)^6)",
        "768*G_2(lambda,M)^6",
        "H_(4,n)(lambda)>0",
        "Delta_k^m log R_T^(1)(k)=O(log(k)/k^m)",
        "closed:",
    ):
        if required not in exact_text:
            findings.append(issue("exact", "missing-formula", required))

    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        findings.append(issue("rows", "bad-row-ids", sorted(str(value) for value in ids)))
    open_rows = [row for row in rows if row.get("readiness") == "not_ready_to_apply"]
    if len(open_rows) != 1 or open_rows[0].get("id") != "co4uhet_07_uniform_heat_tilt_target":
        findings.append(issue("rows", "bad-open-row", open_rows))
    conditional = [row for row in rows if row.get("readiness") == "conditional"]
    if len(conditional) != 2:
        findings.append(issue("rows", "bad-conditional-count", len(conditional)))
    if any(not row.get("claim") or not row.get("proof_boundary") for row in rows):
        findings.append(issue("rows", "incomplete-row", rows))

    try:
        recomputed = build_artifact()
        for key in (
            "source_hashes",
            "source_diagnostics",
            "symbolic_cancellation",
            "exact",
            "rows",
            "summary",
        ):
            if artifact.get(key) != recomputed.get(key):
                findings.append(issue("recompute", f"mismatch-{key}", "stored differs"))
    except Exception as exc:
        findings.append(issue("recompute", "exception", exc))
    return findings


def validate_note(path: Path) -> list[UniformIssue]:
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
    findings: list[UniformIssue] = []
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
        print(f"order-four uniform-heat eventual tail: {len(findings)} issues")
        for finding in findings:
            print(f"- [{finding.section}] {finding.issue}: {finding.detail}")
    else:
        print(
            "validated order-four uniform-heat eventual tail: "
            "7 rows, 0 issues, 6 exact/input rows, 7 symbolic coefficients, "
            "1 conditional uniform transfer, 1 open heat-tilt theorem"
        )
    raise SystemExit(1 if findings else 0)


if __name__ == "__main__":
    main()
