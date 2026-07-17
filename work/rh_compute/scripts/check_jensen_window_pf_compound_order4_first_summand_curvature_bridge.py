#!/usr/bin/env python3
"""Validate the order-four first-summand curvature bridge."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
import json
from pathlib import Path

import sympy as sp

from jensen_window_pf_compound_order4_first_summand_curvature_bridge import (
    DEFAULT_NOTE,
    DEFAULT_OUT,
    DEFAULT_SAMPLE_T,
    REPO_ROOT,
    exact_diagnostics,
    finite_j319_diagnostics,
    compact_curvature_diagnostics,
    lower_mode_collar_diagnostics,
)


REQUIRED_IDS = {
    "co4fcb_01_stable_gap_coordinate",
    "co4fcb_02_curvature_identity",
    "co4fcb_03_tent_bridge",
    "co4fcb_04_gap_floor",
    "co4fcb_04b_localized_curvature",
    "co4fcb_04c_compact_curvature_certificate",
    "co4fcb_05_continuous_curvature_target",
    "co4fcb_06_first_summand_transfer",
    "co4fcb_07_full_kernel_transfer",
    "co4fcb_08_all_tail_implication",
    "co4fcb_09_finite_scout",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Compound Order-Four First-Summand Curvature Bridge",
    "g(t)=d(t)^2*(1-exp(-J(t)))",
    "K(t)=(log g)''=2*ell''+phi(J)*J''-chi(J)*(J')^2",
    "J_1(t)>=2/(2*t+3)-2/(2*t-1)+(t-3)/(6*t^2)>=1/(7*t), t>=319",
    "K_1(t)<=7/(2*t^2), 319<=t<=V'(2)",
    "K_1(t)<=7/(2*t^2) on the mode ray u>=2",
    "This is a real-parameter interval theorem, not finite sampling.",
    "E_r=(1/12)*sup_|s|<=1 |ell^(r+4)(t+s)|",
    "j_0>E_0 and U<=7/(2*t^2) imply J>0 and K<=7/(2*t^2)",
    "|P_n-P_n^(1)|<=10112*(k+1)^2/(k-3)^6<=2/(5*k^2), k>=320",
    "P_n<=18/(5*k^2)+2/(5*k^2)=4/k^2, k=n+3>=320",
    "not interval enclosures",
)

REQUIRED_CORE_STRINGS = (
    "### Lemma 11.31: Order-Four First-Summand Curvature Reduction",
    "K_1(t)<=7/(2t^2)",
    "10112*(k+1)^2/(k-3)^6",
    "E_r=(1/12)*sup",
    "319<=t<=V'(2)",
)


@dataclass(frozen=True)
class GateIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: object) -> GateIssue:
    return GateIssue(section=section, issue=name, detail=str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def exact_issues() -> list[GateIssue]:
    findings: list[GateIssue] = []
    bm, b, bp = sp.symbols("bm b bp", positive=True)
    dm, d, dp = 1 - sp.exp(-bm), 1 - sp.exp(-b), 1 - sp.exp(-bp)
    j_value = 2 * b - sp.log(dm) + 2 * sp.log(d) - sp.log(dp)
    direct = d**2 - sp.exp(-2 * b) * dm * dp
    stable = d**2 * (1 - sp.exp(-j_value))
    if sp.simplify(direct - stable) != 0:
        findings.append(issue("exact", "bad-stable-factorization", direct - stable))

    diagnostics = exact_diagnostics()
    if diagnostics["j_floor_shifted_coefficients_descending"] != [
        4,
        3412,
        955637,
        87486770,
    ]:
        findings.append(
            issue(
                "exact",
                "bad-j-floor-polynomial",
                diagnostics["j_floor_shifted_coefficients_descending"],
            )
        )
    if diagnostics["perturbation_shifted_coefficients_descending"] != [
        1,
        1902,
        1482055,
        604691300,
        135889991935,
        15877422036942,
        748002501678169,
    ]:
        findings.append(
            issue(
                "exact",
                "bad-transfer-polynomial",
                diagnostics["perturbation_shifted_coefficients_descending"],
            )
        )
    k = sp.symbols("k", positive=True)
    margin = sp.factor(sp.Rational(18, 5) / k**2 - sp.Rational(7, 2) / (k**2 - 1))
    expected = (k**2 - 36) / (10 * k**2 * (k**2 - 1))
    if sp.factor(margin - expected) != 0:
        findings.append(issue("exact", "bad-first-summand-margin", margin))
    return findings


def source_issues() -> list[GateIssue]:
    findings: list[GateIssue] = []
    requirements = {
        "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.json": (
            "kappa_3,t(2*log(U))>=-37/(50*t^2), t>=318",
            "L_k^(1)>=1/(k+1/2)^2-37/(50*(k-1)^2)>=1/(4*k^2), k>=319",
        ),
        "work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.json": (
            "C_(k+1)>=C_k for every integer k>=1 at lambda=-100",
        ),
        "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.json": (
            "M_k=M_k^(1)*(1+delta_k), 0<=delta_k<=2/k^6",
        ),
    }
    for relative, formulas in requirements.items():
        path = REPO_ROOT / relative
        try:
            artifact = load_json(path)
        except Exception as exc:
            findings.append(issue("source", "load-failed", (relative, exc)))
            continue
        present = {row.get("formula") for row in artifact.get("rows", [])}
        for formula in formulas:
            if formula not in present:
                findings.append(issue("source", "missing-formula", (relative, formula)))
    return findings


def artifact_issues(artifact: dict) -> list[GateIssue]:
    findings: list[GateIssue] = []
    if artifact.get("kind") != "jensen_window_pf_compound_order4_first_summand_curvature_bridge":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    expected_status = (
        "exact order-four first-summand curvature reduction with proved gap floor, "
        "a compact curvature certificate, closed full-kernel perturbation budget, "
        "and one open analytic ray"
    )
    if artifact.get("status") != expected_status:
        findings.append(issue("artifact", "bad-status", artifact.get("status")))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for marker in (
        "proves",
        "does not prove",
        "analytic curvature",
        "order-four entry",
        "forward",
        "arbitrary-column",
        "rh",
        "lambda<=0",
    ):
        if marker not in boundary:
            findings.append(issue("artifact", "weak-boundary", marker))
    ids = {row.get("id") for row in artifact.get("rows", [])}
    if ids != REQUIRED_IDS:
        findings.append(issue("rows", "bad-ids", sorted(str(value) for value in ids)))
    expected_summary = {
        "bridge_rows": 11,
        "exact_identity_rows": 2,
        "exact_reduction_rows": 2,
        "localized_curvature_rows": 1,
        "compact_interval_certificate_rows": 1,
        "proved_gap_floor_rows": 1,
        "lower_mode_collar_blocks": 2,
        "proved_perturbation_rows": 1,
        "conditional_bridge_rows": 2,
        "open_analytic_handoffs": 1,
        "finite_scout_rows": len(DEFAULT_SAMPLE_T),
        "ready_to_apply_rows": 7,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(issue("summary", "mismatch", artifact.get("summary")))
    open_rows = [
        row
        for row in artifact.get("rows", [])
        if row.get("role") == "open_analytic_handoff"
    ]
    if len(open_rows) != 1 or open_rows[0].get("id") != "co4fcb_05_continuous_curvature_target":
        findings.append(issue("rows", "bad-open-frontier", open_rows))
    scouts = artifact.get("scout", {}).get("rows", [])
    if [row.get("t") for row in scouts] != list(DEFAULT_SAMPLE_T):
        findings.append(issue("scout", "bad-sample", [row.get("t") for row in scouts]))
    for row in scouts:
        try:
            if not row.get("positive_gap") or not row.get("below_target"):
                findings.append(issue("scout", "failed-row", row.get("t")))
            if Decimal(row["margin_below_7_over_2"]) <= 0:
                findings.append(issue("scout", "nonpositive-margin", row.get("t")))
        except Exception as exc:
            findings.append(issue("scout", "bad-decimal", exc))
    try:
        rebuilt_finite = finite_j319_diagnostics()
    except Exception as exc:
        findings.append(issue("finite", "recompute-failed", exc))
    else:
        if artifact.get("finite_j319") != rebuilt_finite:
            findings.append(issue("finite", "recompute-mismatch", "J_319 diagnostics changed"))
    try:
        rebuilt_collar = lower_mode_collar_diagnostics()
    except Exception as exc:
        findings.append(issue("collar", "recompute-failed", exc))
    else:
        if artifact.get("lower_mode_collar") != rebuilt_collar:
            findings.append(issue("collar", "recompute-mismatch", "lower-mode collar changed"))
    try:
        rebuilt_compact = compact_curvature_diagnostics()
    except Exception as exc:
        findings.append(issue("compact", "recompute-failed", exc))
    else:
        if artifact.get("compact_curvature") != rebuilt_compact:
            findings.append(
                issue("compact", "recompute-mismatch", "compact curvature diagnostics changed")
            )
    for ref in artifact.get("sources", []):
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            findings.append(issue("artifact", "missing-source", ref))
    for key in ("generator", "checker"):
        ref = artifact.get(key)
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            findings.append(issue("artifact", "missing-ref", (key, ref)))
    return findings


def required_string_issues(
    path: Path, required: tuple[str, ...], section: str
) -> list[GateIssue]:
    if not path.exists():
        return [issue(section, "missing", path)]
    text = path.read_text(encoding="utf-8")
    return [
        issue(section, "missing-string", value) for value in required if value not in text
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    findings = exact_issues()
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
        print(
            json.dumps(
                {"ok": ok, "issues": [asdict(item) for item in findings]},
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for finding in findings:
            print(f"ISSUE {finding.section} [{finding.issue}] {finding.detail}")
        print(
            "validated Jensen-window PF compound order-four first-summand curvature bridge: "
            f"11 rows, {len(findings)} issues, 2 exact identities, 2 exact "
            "reductions, 1 proved gap floor, 1 compact interval theorem, 1 proved "
            "full-kernel perturbation theorem, 1 open analytic ray, "
            f"{len(DEFAULT_SAMPLE_T)} positive finite scouts"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
