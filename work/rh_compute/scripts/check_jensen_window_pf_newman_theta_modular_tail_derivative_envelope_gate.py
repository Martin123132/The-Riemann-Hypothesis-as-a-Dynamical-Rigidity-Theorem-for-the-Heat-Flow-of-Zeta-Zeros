#!/usr/bin/env python3
"""Validate the exact modular-tail derivative envelope gate."""

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

import jensen_window_pf_newman_theta_modular_tail_derivative_envelope_gate as target  # noqa: E402


EXPECTED_IDS = [
    "ntmtdeg_01_heat_hermite_polynomials",
    "ntmtdeg_02_value_derivative_envelope",
    "ntmtdeg_03_first_jet_derivative_envelope",
    "ntmtdeg_04_reflected_phase_minimum",
    "ntmtdeg_05_stretched_exponential_tail",
    "ntmtdeg_06_fixed_collar_scale_guard",
    "ntmtdeg_07_three_quarter_scale",
    "ntmtdeg_08_revised_adaptive_candidate",
    "ntmtdeg_09_open_separation_handoff",
]

NOTE_MARKERS = (
    "# Newman Theta Modular-Tail Derivative Envelope Gate",
    "## Explicit Heat-Derivative Envelope",
    "## Reflected-Blend Phase",
    "## Correct Absolute-Error Scale",
    "## Open Separation Target",
    "P_b(T,v)=b!*sum_(ell=0)^floor(b/2)",
    "exp(-c_*N^(4/3))",
    "N_K(x)=ceil(K*(1+x)^(3/4))",
    "N(x,t)=max(N_sad(x,t),ceil(K*(1+x)^(3/4)))",
    "No new `Q_j`",
)

CORE_MARKERS = (
    "modular-tail derivative envelope",
    "P_b(T,v)",
    "exp(-c_*N^(4/3))",
    "N_K(x)=ceil(K*(1+x)^(3/4))",
    "fixed additive saddle collar",
)


@dataclass(frozen=True)
class GateIssue:
    section: str
    issue: str
    detail: str


def finding(section: str, issue: str, detail: object) -> GateIssue:
    return GateIssue(section, issue, str(detail))


def artifact_issues(artifact: dict) -> list[GateIssue]:
    issues: list[GateIssue] = []
    if artifact.get("kind") != target.STEM:
        issues.append(finding("artifact", "kind", artifact.get("kind")))
    if artifact.get("date") != "2026-07-24":
        issues.append(finding("artifact", "date", artifact.get("date")))
    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append(
            finding("rows", "id-order", [row.get("id") for row in rows])
        )
    if len(rows) != 9:
        issues.append(finding("rows", "count", len(rows)))
    if rows and rows[-1].get("readiness") != "not_ready_to_apply":
        issues.append(
            finding("rows", "open-target-promoted", rows[-1].get("readiness"))
        )

    exact = artifact.get("exact", {})
    hermite = exact.get("hermite_envelope", {})
    if hermite.get("verified_orders") != list(range(10)):
        issues.append(
            finding("exact", "verified-orders", hermite.get("verified_orders"))
        )
    try:
        c0 = sp.Float(exact["reflected_phase"]["c0_decimal"], 80)
        c_safe = sp.Float(
            exact["derivative_tail_theorem"]["c_safe_decimal"], 80
        )
        k_gamma = sp.Float(
            exact["frequency_scales"]["K_gamma_decimal"], 80
        )
        if abs(c_safe - c0 / 4) > sp.Float("1e-28"):
            issues.append(finding("exact", "safe-exponent", c_safe - c0 / 4))
        expected_k = (sp.pi / (8 * c_safe)) ** sp.Rational(3, 4)
        if abs(k_gamma - expected_k) > sp.Float("1e-28"):
            issues.append(
                finding("exact", "gamma-threshold", k_gamma - expected_k)
            )
    except (KeyError, TypeError, ValueError) as exc:
        issues.append(finding("exact", "constants", exc))

    audit = artifact.get("source_audit", {})
    if audit.get("modular_kind") != (
        "jensen_window_pf_newman_theta_modular_blend_gate"
    ):
        issues.append(
            finding("source-audit", "modular-kind", audit.get("modular_kind"))
        )
    if audit.get("c1_kind") != (
        "jensen_window_pf_newman_theta_adaptive_modular_c1_remainder_contract"
    ):
        issues.append(finding("source-audit", "c1-kind", audit.get("c1_kind")))
    if audit.get("modular_rows") != 12:
        issues.append(
            finding("source-audit", "modular-rows", audit.get("modular_rows"))
        )
    if audit.get("c1_rows") != 10:
        issues.append(
            finding("source-audit", "c1-rows", audit.get("c1_rows"))
        )

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not instantiate",
        "first-jet separation",
        "Lambda<=0",
        "RH",
        "Clay-prize",
    ):
        if marker not in boundary:
            issues.append(finding("artifact", "proof-boundary", marker))
    if artifact != target.build_artifact():
        issues.append(
            finding(
                "recompute",
                "artifact-mismatch",
                "stored artifact differs from exact rebuild",
            )
        )
    return issues


def text_issues(
    path: Path, markers: tuple[str, ...], section: str
) -> list[GateIssue]:
    text = path.read_text(encoding="utf-8")
    return [
        finding(section, "missing-marker", marker)
        for marker in markers
        if marker not in text
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=target.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=target.DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    issues = artifact_issues(artifact)
    issues.extend(text_issues(args.note, NOTE_MARKERS, "note"))
    issues.extend(
        text_issues(
            target.REPO_ROOT / "outputs/formal_core.md",
            CORE_MARKERS,
            "formal-core",
        )
    )
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"ISSUE {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Newman theta modular-tail derivative envelope gate: "
            f"9 rows, {len(issues)} issues, 3 exact inequalities/identities, "
            "3 exact tail-scale theorems/compositions, "
            "1 fixed-collar guard, 1 open separation handoff"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
