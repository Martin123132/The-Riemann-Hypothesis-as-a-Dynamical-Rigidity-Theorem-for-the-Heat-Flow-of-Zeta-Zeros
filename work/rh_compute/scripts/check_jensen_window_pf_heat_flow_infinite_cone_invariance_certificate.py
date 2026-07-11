#!/usr/bin/env python3
"""Validate the infinite ratio-cone maximum-principle certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_heat_flow_infinite_cone_invariance_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "hfic_01_actual_analytic_trajectory",
    "hfic_02_defect_coordinates",
    "hfic_03_exact_adjacent_ode",
    "hfic_04_uniform_coefficient_bound",
    "hfic_05_infinite_minimum_principle",
    "hfic_06_full_cone_propagation",
    "hfic_07_endpoint_cone",
    "hfic_08_all_order_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Heat-Flow Infinite Cone-Invariance Certificate",
    "Status: exact infinite ratio-cone propagation theorem and open all-order handoff",
    "Artifact kind: `jensen_window_pf_heat_flow_infinite_cone_invariance_certificate`",
    "work/rh_compute/results/jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.json",
    "python work/rh_compute/scripts/jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.py",
    "d_k=1-x_k",
    "h_k=d_k-d_(k+1)=x_(k+1)-x_k",
    "h_k'=a_k*(h_(k+1)-h_k)+q_k*h_k+c_k",
    "a_k>=0, c_k>=0, q_k<=176*R",
    "spatial infimum is attained at a finite index",
    "m(t)=min(0,inf_(k>=1) z_k(t))",
    "D_+m(t)=min{z_k'(t):z_k(t)=m(t)}>=0",
    "connected component `(alpha,beta)` of `{t:m(t)<0}`",
    "negative component cannot form",
    "full infinite ratio cone at `lambda=0`",
    "ratio-cone propagation alone is not PF-infinity",
    "outputs/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.md",
    "outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md",
    "outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class CertificateIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: object) -> CertificateIssue:
    return CertificateIssue(section=section, issue=name, detail=str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(ref: object) -> list[CertificateIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue("artifact", "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue("artifact", "missing-ref", ref)]
    return []


def symbolic_issues() -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    k, r, d, h, j = sp.symbols("k r d h j", positive=True)
    d1 = d - h
    d2 = d - h - j
    d_prime = 2 * r * ((2 * k + 3) * (1 - d) * d1 - (2 * k - 1) * d)
    d1_prime = (
        2
        * r
        * (1 - d1)
        * ((2 * k + 5) * (1 - d1) * d2 - (2 * k + 1) * d1)
    )
    h_prime = sp.expand(d_prime - d1_prime)
    A = (2 * k + 5) * (1 - d + h) ** 2
    Q = (
        (8 * k + 20) * d**2
        - (10 * k + 25) * d * h
        - (6 * k + 25) * d
        + (4 * k + 10) * h**2
        + (6 * k + 19) * h
        + 6
    )
    B = sp.expand(Q - A)
    source = d**2 * (6 - (2 * k + 5) * d)
    expected = 2 * r * (source + B * h + A * j)
    if sp.expand(h_prime - expected) != 0:
        findings.append(issue("symbolic", "bad-adjacent-ode", sp.factor(h_prime - expected)))
    transport = 2 * r * (A * (j - h) + Q * h + source)
    if sp.expand(h_prime - transport) != 0:
        findings.append(issue("symbolic", "bad-transport-form", sp.factor(h_prime - transport)))
    boundary_margin = sp.factor(6 - (2 * k + 5) * 2 / (2 * k + 1))
    if sp.factor(boundary_margin - (8 * k - 4) / (2 * k + 1)) != 0:
        findings.append(issue("symbolic", "bad-source-margin", boundary_margin))
    return findings


def exact_bound_issues(diagnostics: dict) -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    terms = [
        Fraction(16, 3),
        Fraction(80, 9),
        Fraction(20, 3),
        Fraction(100, 9),
        Fraction(6, 1),
        Fraction(50, 3),
        Fraction(8, 3),
        Fraction(40, 9),
        Fraction(6, 1),
        Fraction(38, 3),
        Fraction(6, 1),
    ]
    total = sum(terms, Fraction(0, 1))
    if total != Fraction(778, 9) or not total < 88:
        findings.append(issue("bounds", "bad-Q-bound", total))
    stored = diagnostics.get("uniform_coefficient_bound", {})
    if stored.get("Q_absolute_sum") != "778/9" or stored.get("Q_cap") != 88:
        findings.append(issue("bounds", "bad-stored-Q-bound", stored))
    pointwise = diagnostics.get("pointwise_cone", {})
    if pointwise.get("uniform_tail") != "sup_lambda |h_k(lambda)|<=1/k on compact forward intervals":
        findings.append(issue("bounds", "bad-tail-bound", pointwise))
    maximum = diagnostics.get("maximum_principle", {})
    for text in (
        "negative spatial minimum",
        "uniformly",
        "finite active set",
        "D_+m(t)",
        "connected component",
        "C1",
        "h_k(lambda)>=0",
    ):
        if text not in " ".join(str(value) for value in maximum.values()):
            findings.append(issue("maximum", "missing-argument", text))
    return findings


def source_issues() -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    paths = {
        "entry": REPO_ROOT
        / "work/rh_compute/results/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.json",
        "upper": REPO_ROOT
        / "work/rh_compute/results/jensen_window_pf_kernel_mellin_upper_wall_certificate.json",
        "boundary": REPO_ROOT
        / "work/rh_compute/results/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.json",
    }
    try:
        sources = {key: load_json(path) for key, path in paths.items()}
    except Exception as exc:
        return [issue("sources", "load-failed", exc)]
    entry = sources["entry"].get("summary", {})
    if entry.get("full_cone_entry_rows") != 1 or entry.get("ready_to_apply_rows") != 3:
        findings.append(issue("sources", "bad-entry-source", entry))
    upper = sources["upper"].get("summary", {})
    if upper.get("all_k_upper_wall_certified") is not True:
        findings.append(issue("sources", "bad-pointwise-source", upper))
    boundary = sources["boundary"].get("summary", {})
    if boundary.get("conditional_invariance_available") is not True:
        findings.append(issue("sources", "bad-boundary-source", boundary))
    return findings


def validate_artifact(artifact: dict) -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    if artifact.get("kind") != "jensen_window_pf_heat_flow_infinite_cone_invariance_certificate":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != "exact infinite ratio-cone propagation theorem and open all-order handoff":
        findings.append(issue("artifact", "bad-status", artifact.get("status")))
    if artifact.get("date") != "2026-07-10":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    for key in (
        "source_cone_entry",
        "source_pointwise_walls",
        "source_boundary_algebra",
        "source_heat_target",
        "generator",
        "checker",
    ):
        findings.extend(validate_ref(artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for text in ("proves", "does not prove", "full infinite", "lambda=-100", "lambda=0", "rh", "lambda <= 0"):
        if text not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", text))

    summary = artifact.get("summary", {})
    expected_summary = {
        "certificate_rows": 8,
        "exact_ode_rows": 3,
        "infinite_maximum_principle_rows": 1,
        "full_cone_propagation_rows": 1,
        "endpoint_cone_rows": 1,
        "open_all_order_handoffs": 1,
        "ready_to_apply_rows": 3,
    }
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            findings.append(issue("summary", f"bad-{key}", summary.get(key)))

    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        findings.append(issue("rows", "bad-row-ids", sorted(str(value) for value in ids)))
    if len(rows) != 8:
        findings.append(issue("rows", "bad-row-count", len(rows)))
    ready = [row.get("id") for row in rows if row.get("readiness") == "ready_to_apply"]
    if ready != ["hfic_05_infinite_minimum_principle", "hfic_06_full_cone_propagation", "hfic_07_endpoint_cone"]:
        findings.append(issue("rows", "bad-ready-rows", ready))

    findings.extend(symbolic_issues())
    findings.extend(exact_bound_issues(artifact.get("diagnostics", {})))
    findings.extend(source_issues())
    try:
        recomputed = build_artifact()
        for key in ("diagnostics", "summary", "rows"):
            if artifact.get(key) != recomputed.get(key):
                findings.append(issue("recompute", f"mismatch-{key}", "stored differs"))
    except Exception as exc:
        findings.append(issue("recompute", "exception", exc))
    return findings


def validate_note(path: Path) -> list[CertificateIssue]:
    if not path.exists():
        return [issue("note", "missing", path)]
    text = path.read_text(encoding="utf-8")
    return [issue("note", "missing-text", value) for value in REQUIRED_NOTE_STRINGS if value not in text]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    findings: list[CertificateIssue] = []
    try:
        artifact = load_json(args.artifact)
    except Exception as exc:
        artifact = {}
        findings.append(issue("artifact", "load-failed", exc))
    if artifact:
        findings.extend(validate_artifact(artifact))
    findings.extend(validate_note(args.note))
    ok = not findings
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in findings]}, indent=2, sort_keys=True))
    else:
        for finding in findings:
            print(f"ISSUE {finding.section} [{finding.issue}] {finding.detail}")
        print(
            "validated Jensen-window PF heat-flow infinite cone-invariance certificate: "
            f"8 rows, {len(findings)} issues, 1 infinite maximum principle, "
            "1 full cone propagation, 1 endpoint cone theorem, "
            "1 open all-order handoff, 3 ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
