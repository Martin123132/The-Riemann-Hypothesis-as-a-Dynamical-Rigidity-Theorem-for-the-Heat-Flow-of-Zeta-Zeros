#!/usr/bin/env python3
"""Validate the fixed-block theta cofinal obstruction gate."""

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

import jensen_window_pf_newman_theta_fixed_block_cofinal_obstruction_gate as target  # noqa: E402


EXPECTED_IDS = [
    "ntfbco_01_uniform_regular_decay",
    "ntfbco_02_value_ibp_bound",
    "ntfbco_03_derivative_ibp_bound",
    "ntfbco_04_normalized_static_bars",
    "ntfbco_05_cofinal_obstruction",
    "ntfbco_06_two_block_endpoint_cancellation",
    "ntfbco_07_finite_frontier_diagnostic",
    "ntfbco_08_cancellation_aware_handoff",
    "ntfbco_09_proof_boundary",
]

REQUIRED_NOTE_MARKERS = (
    "# Newman Theta Fixed-Block Cofinal Obstruction Gate",
    "## Two Integrations By Parts",
    "## Exact Cofinal Obstruction",
    "## What The Bars Erase",
    "## Finite Scout",
    "## Live Handoff",
    "A_N/(epsilon_0*x^2)",
    "X_N=max(1,sqrt(2A_N/epsilon_0),B_N/(2epsilon_0))",
    "a_2=-sum_(n>=3)phi_n'(0)>0",
    "sum_(n>=3)phi_n'(0)=-a_2<0",
    "shells=5..30",
    "No common zero, new `Q_j`, `Lambda<=0`, RH",
)

REQUIRED_CORE_MARKERS = (
    "f_(N,t)(u)=exp(tu^2) sum_(1<=n<=N) Phi_n(u)",
    "A_N=sup_(0<=t<=1/5)(|f_(N,t)'(0)|+||f_(N,t)''||_1)",
    "|J_(N,t)|/(16x^4 epsilon_0)",
    "No fixed finite theta truncation with fixed positive static",
    "a_2=f_(2,t)'(0)=sum_(n=1)^2 Phi_n'(0)",
    "all 26 rows certify their rational box covers",
    "cancellation-aware first-jet remainder theorem",
    "Failure of the static sufficient bars is not evidence",
)


@dataclass(frozen=True)
class GateIssue:
    section: str
    issue: str
    detail: str


def finding(
    section: str, issue: str, detail: object
) -> GateIssue:
    return GateIssue(section, issue, str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def symbolic_issues() -> list[GateIssue]:
    issues: list[GateIssue] = []
    u = sp.symbols("u", real=True)
    n = sp.symbols("n", integer=True, positive=True)
    theta = (
        2 * sp.pi**2 * n**4 * sp.exp(9 * u)
        - 3 * sp.pi * n**2 * sp.exp(5 * u)
    ) * sp.exp(-sp.pi * n**2 * sp.exp(4 * u))
    expected = (
        sp.pi
        * n**2
        * sp.exp(-sp.pi * n**2)
        * (
            -8 * sp.pi**2 * n**4
            + 30 * sp.pi * n**2
            - 15
        )
    )
    if sp.simplify(sp.diff(theta, u).subs(u, 0) - expected) != 0:
        issues.append(
            finding(
                "symbolic",
                "endpoint-derivative",
                "theta derivative identity failed",
            )
        )

    a_two = sp.simplify(expected.subs(n, 1) + expected.subs(n, 2))
    if not sp.N(a_two, 90) > 0:
        issues.append(
            finding("symbolic", "two-block-sign", sp.N(a_two, 90))
        )
    if not (
        sp.N(a_two, 90) > sp.Float("8.26e-8", 90)
        and sp.N(a_two, 90) < sp.Float("8.27e-8", 90)
    ):
        issues.append(
            finding(
                "symbolic",
                "two-block-size",
                sp.N(a_two, 90),
            )
        )

    x, h, hp, e0, e1 = sp.symbols(
        "x h hp e0 e1", positive=True
    )
    value_ratio = (16 * x**4 * h) / (16 * x**4 * e0)
    derivative_ratio = (
        64 * x**3 * h + 16 * x**4 * hp
    ) / (64 * x**3 * e0 + 16 * x**4 * e1)
    if sp.simplify(value_ratio - h / e0) != 0:
        issues.append(
            finding("symbolic", "value-ratio", value_ratio)
        )
    expected_derivative = (4 * h + x * hp) / (4 * e0 + x * e1)
    if sp.simplify(derivative_ratio - expected_derivative) != 0:
        issues.append(
            finding(
                "symbolic",
                "derivative-ratio",
                derivative_ratio,
            )
        )

    return issues


def artifact_issues(artifact: dict) -> list[GateIssue]:
    issues: list[GateIssue] = []
    if artifact.get("kind") != target.STEM:
        issues.append(
            finding("artifact", "kind", artifact.get("kind"))
        )
    if artifact.get("date") != "2026-07-24":
        issues.append(
            finding("artifact", "date", artifact.get("date"))
        )
    if artifact.get("status") != (
        "exact fixed-block cofinal obstruction with "
        "cancellation-aware handoff"
    ):
        issues.append(
            finding("artifact", "status", artifact.get("status"))
        )

    boundary = str(artifact.get("proof_boundary", ""))
    for marker in (
        "fixed finite theta truncation",
        "cannot certify the cofinal",
        "does not prove a common zero",
        "certify any new Q_j",
        "Lambda<=0",
        "prove RH",
        "Clay prize",
    ):
        if marker not in boundary:
            issues.append(
                finding("artifact", "proof-boundary", marker)
            )

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append(
            finding(
                "rows",
                "id-order",
                [row.get("id") for row in rows],
            )
        )
    role_counts = {
        role: sum(row.get("role") == role for row in rows)
        for role in (
            "exact_lemma",
            "exact_identity",
            "exact_inequality",
            "route_guard",
            "exact_theorem",
            "finite_diagnostic",
            "open_handoff",
            "proof_guard",
        )
    }
    expected_counts = {
        "exact_lemma": 1,
        "exact_identity": 2,
        "exact_inequality": 1,
        "route_guard": 1,
        "exact_theorem": 1,
        "finite_diagnostic": 1,
        "open_handoff": 1,
        "proof_guard": 1,
    }
    if role_counts != expected_counts:
        issues.append(
            finding("rows", "role-counts", role_counts)
        )
    if rows and rows[7].get("readiness") != "not_ready_to_apply":
        issues.append(
            finding(
                "rows",
                "open-handoff-promoted",
                rows[7].get("readiness"),
            )
        )

    audit = artifact.get("source_audit", {}).get(
        "frontier_prefix", {}
    )
    expected_audit = {
        "start_j": 5,
        "end_j": 30,
        "shell_rows": 26,
        "certified_shell_rows": 26,
        "minimum_ratio_shell": 23,
        "terminal_prefix_row_sha256": (
            "a2e340abc6f317d0c8f6e794e3f9f4516e61db5e5ff23e159d82de080fde732d"
        ),
    }
    for key, expected in expected_audit.items():
        if audit.get(key) != expected:
            issues.append(
                finding(
                    "frontier",
                    key,
                    f"{audit.get(key)!r} != {expected!r}",
                )
            )
    if audit.get("total_evaluated_boxes", 0) <= 0:
        issues.append(
            finding(
                "frontier",
                "evaluated-boxes",
                audit.get("total_evaluated_boxes"),
            )
        )
    if audit.get("total_subdivisions") != 1292:
        issues.append(
            finding(
                "frontier",
                "subdivisions",
                audit.get("total_subdivisions"),
            )
        )

    exact = artifact.get("exact", {})
    if exact.get("static_tail_certificate", {}).get(
        "explicit_threshold"
    ) != (
        "X_N=max(1,sqrt(2A_N/epsilon_0),"
        "B_N/(2epsilon_0))"
    ):
        issues.append(
            finding(
                "exact",
                "threshold",
                exact.get("static_tail_certificate", {}).get(
                    "explicit_threshold"
                ),
            )
        )
    if "not evidence of a common zero" not in exact.get(
        "cofinal_obstruction", {}
    ).get("scope", ""):
        issues.append(
            finding(
                "exact",
                "scope",
                exact.get("cofinal_obstruction", {}).get("scope"),
            )
        )

    rebuilt = target.build_artifact()
    if artifact != rebuilt:
        issues.append(
            finding(
                "recompute",
                "artifact-mismatch",
                "stored artifact differs from exact rebuild",
            )
        )
    return issues


def note_issues(path: Path) -> list[GateIssue]:
    issues: list[GateIssue] = []
    text = path.read_text(encoding="utf-8")
    for marker in REQUIRED_NOTE_MARKERS:
        if marker not in text:
            issues.append(finding("note", "missing-marker", marker))
    return issues


def core_issues() -> list[GateIssue]:
    issues: list[GateIssue] = []
    path = target.REPO_ROOT / "outputs/formal_core.md"
    text = path.read_text(encoding="utf-8")
    for marker in REQUIRED_CORE_MARKERS:
        if marker not in text:
            issues.append(
                finding("formal-core", "missing-marker", marker)
            )
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--artifact", type=Path, default=target.DEFAULT_OUT
    )
    parser.add_argument(
        "--note", type=Path, default=target.DEFAULT_NOTE
    )
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    issues = symbolic_issues()
    issues.extend(artifact_issues(load_json(args.artifact)))
    issues.extend(note_issues(args.note))
    issues.extend(core_issues())
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
            print(
                f"ISSUE {item.section} "
                f"[{item.issue}] {item.detail}"
            )
        print(
            "validated Newman theta fixed-block cofinal obstruction "
            f"gate: 9 rows, {len(issues)} issues, "
            "3 exact integration-by-parts bounds, "
            "1 cofinal method obstruction, "
            "1 endpoint-cancellation theorem, "
            "26 finite diagnostic shells, "
            "1 cancellation-aware open handoff"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
