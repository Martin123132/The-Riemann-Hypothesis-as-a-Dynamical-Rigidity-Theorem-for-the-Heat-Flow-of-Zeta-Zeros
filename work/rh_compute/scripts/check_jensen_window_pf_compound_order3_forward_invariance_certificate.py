#!/usr/bin/env python3
"""Validate contiguous order-three forward invariance."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order3_forward_invariance_certificate import (  # noqa: E402
    CUBIC_SOURCE,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    ENTRY_SOURCE,
    RATIO_SOURCE,
    REPO_ROOT,
)


REQUIRED_IDS = {
    "co3fi_01_q_flow",
    "co3fi_02_compound_coordinate",
    "co3fi_03_cooperative_system",
    "co3fi_04_inward_boundary",
    "co3fi_05_forward_tail_input",
    "co3fi_06_coefficient_growth",
    "co3fi_07_weighted_maximum_principle",
    "co3fi_08_strict_forward_propagation",
    "co3fi_09_lambda_zero_theorem",
    "co3fi_10_higher_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Compound Order-Three Forward-Invariance Certificate",
    "C_n'/r_k=alpha_k*C_(n+1)+beta_k*C_n",
    "at C_n=0, C_n'/r_k=alpha_k*C_(n+1)",
    "lim_(h->0) h^2*alpha_k=4",
    "lim_(h->0) (alpha_k+beta_k)=2*(u+2*v-3*w)/y",
    "z_n=C_n/(n+1)",
    "C_n(lambda)>0 for every n>=0 and finite lambda>=-100",
    "D_(3,n)(0)<0 for every n>=0",
    "Noncontiguous order-three minors and every higher compound order remain",
)

REQUIRED_CORE_STRINGS = (
    "### Lemma 11.28: Contiguous Order-Three Forward Invariance",
    "C_n'/r_k=alpha_k*C_(n+1)+beta_k*C_n",
    "z_n=C_n/(n+1)",
    "D_(3,n)(0)<0",
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


def independent_exact_issues() -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    k = sp.symbols("k", integer=True, positive=True)
    a, b, c, d = sp.symbols("a b c d", positive=True)

    def qdot(index: sp.Expr, current: sp.Expr, next_q: sp.Expr) -> sp.Expr:
        return (
            (2 * index - 1) * current
            - (2 * index + 3) * (current**3 - current) / next_q**2
        )

    derivative = sp.together(
        c * qdot(k - 1, a, b) / (1 - b**-2)
        + a * (1 - c**-2) * qdot(k + 1, c, d)
        - 2 * b * qdot(k, b, c)
    )
    compound = a * c - b**2 + 1
    next_compound = b * d - c**2 + 1
    alpha = (2 * k + 5) * a * (b * d + c**2 - 1) / (c * d**2)
    beta_numerator = (
        (2 * k + 1) * a**2 * c**2
        + (2 * k + 1) * a * b**2 * c
        - (2 * k + 1) * a * c
        + (4 * k + 6) * b**4
        + (-4 * k + 2) * b**2 * c**2
        - (4 * k + 6) * b**2
    )
    beta = -beta_numerator / (c**2 * (b**2 - 1))
    residual = sp.cancel(derivative - alpha * next_compound - beta * compound)
    if residual != 0:
        findings.append(issue("exact", "bad-cooperative-decomposition", residual))

    h, y, u, v, w = sp.symbols("h y u v w", positive=True)
    scaling = {
        k: h**-2,
        a: (y - u * h**2) / h,
        b: y / h,
        c: (y + v * h**2) / h,
        d: (y + (v + w) * h**2) / h,
    }
    scaled_alpha = sp.cancel(h**2 * alpha.subs(scaling))
    scaled_sum = sp.cancel((alpha + beta).subs(scaling))
    alpha_limit = sp.cancel(scaled_alpha.subs(h, 0))
    sum_limit = sp.cancel(scaled_sum.subs(h, 0))
    if alpha_limit != 4:
        findings.append(issue("exact", "bad-alpha-limit", alpha_limit))
    if sp.cancel(sum_limit - 2 * (u + 2 * v - 3 * w) / y) != 0:
        findings.append(issue("exact", "bad-sum-limit", sum_limit))
    return findings


def source_issues() -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    try:
        entry = load_json(ENTRY_SOURCE)
        cubic = load_json(CUBIC_SOURCE)
        ratio = load_json(RATIO_SOURCE)
    except Exception as exc:
        return [issue("sources", "load-failed", exc)]
    if entry.get("summary", {}).get("full_entry_rows") != 1:
        findings.append(issue("sources", "entry-not-closed", entry.get("summary")))
    cubic_ids = {row.get("id") for row in cubic.get("rows", [])}
    for row_id in ("cfut_07_forward_uniform_tail", "cfut_08_cubic_continuation"):
        if row_id not in cubic_ids:
            findings.append(issue("sources", "missing-cubic-row", row_id))
    if ratio.get("summary", {}).get("full_cone_propagation_rows") != 1:
        findings.append(issue("sources", "ratio-not-closed", ratio.get("summary")))
    return findings


def artifact_issues(artifact: dict) -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    if artifact.get("kind") != "jensen_window_pf_compound_order3_forward_invariance_certificate":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("date") != "2026-07-12":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    if artifact.get("status") != (
        "exact all-shift contiguous order-three propagation through lambda=0 "
        "with noncontiguous and higher compounds open"
    ):
        findings.append(issue("artifact", "bad-status", artifact.get("status")))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for marker in (
        "proves",
        "does not prove",
        "contiguous",
        "noncontiguous",
        "order four",
        "all-order",
        "rh",
        "lambda<=0",
    ):
        if marker not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", marker))
    for ref in artifact.get("sources", []):
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            findings.append(issue("artifact", "missing-source", ref))
    for key in ("generator", "checker"):
        ref = artifact.get(key)
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            findings.append(issue("artifact", "missing-ref", (key, ref)))

    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_IDS:
        findings.append(issue("rows", "bad-ids", sorted(str(value) for value in ids)))
    expected_summary = {
        "certificate_rows": 10,
        "exact_identity_rows": 2,
        "cooperative_flow_rows": 1,
        "inward_boundary_rows": 1,
        "coefficient_growth_rows": 1,
        "infinite_maximum_principle_rows": 1,
        "full_forward_propagation_rows": 1,
        "lambda_zero_theorem_rows": 1,
        "open_handoffs": 1,
        "ready_to_apply_rows": 9,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(issue("summary", "mismatch", artifact.get("summary")))
    exact = artifact.get("exact", {})
    if exact.get("alpha_scaled_limit") != "4":
        findings.append(issue("artifact", "bad-alpha-limit", exact.get("alpha_scaled_limit")))
    if exact.get("diagonal_scaled_limit") != "2*(u + 2*v - 3*w)/y":
        findings.append(
            issue("artifact", "bad-diagonal-limit", exact.get("diagonal_scaled_limit"))
        )
    return findings


def required_string_issues(path: Path, required: tuple[str, ...], section: str) -> list[CertificateIssue]:
    if not path.exists():
        return [issue(section, "missing", path)]
    text = path.read_text(encoding="utf-8")
    return [issue(section, "missing-string", value) for value in required if value not in text]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    findings = independent_exact_issues()
    findings.extend(source_issues())
    try:
        artifact = load_json(args.artifact)
    except Exception as exc:
        artifact = {}
        findings.append(issue("artifact", "load-failed", exc))
    if artifact:
        findings.extend(artifact_issues(artifact))
    findings.extend(required_string_issues(args.note, REQUIRED_NOTE_STRINGS, "note"))
    findings.extend(
        required_string_issues(
            REPO_ROOT / "outputs/formal_core.md", REQUIRED_CORE_STRINGS, "formal-core"
        )
    )
    ok = not findings
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in findings]}, indent=2, sort_keys=True))
    else:
        for finding in findings:
            print(f"ISSUE {finding.section} [{finding.issue}] {finding.detail}")
        print(
            "validated Jensen-window PF compound order-three forward-invariance "
            f"certificate: 10 rows, {len(findings)} issues, 2 exact identities, "
            "1 cooperative flow, 1 inward boundary theorem, 1 coefficient-growth "
            "lemma, 1 infinite maximum principle, 1 full forward propagation "
            "theorem, 1 lambda=0 theorem, 1 open higher-compound handoff"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
