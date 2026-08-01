#!/usr/bin/env python3
"""Validate the adaptive modular theta C1 remainder contract."""

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

import jensen_window_pf_newman_theta_adaptive_modular_c1_remainder_contract as target  # noqa: E402


EXPECTED_IDS = [
    "ntamc1_01_positive_even_partition",
    "ntamc1_02_normal_transform_series",
    "ntamc1_03_decaying_tail_budget",
    "ntamc1_04_direct_first_jet_error",
    "ntamc1_05_contact_exclusion_disjunction",
    "ntamc1_06_fixed_block_guard",
    "ntamc1_07_adaptive_saddle_scale",
    "ntamc1_08_cancellation_stress",
    "ntamc1_09_open_adaptive_c1_target",
    "ntamc1_10_proof_boundary",
]

NOTE_MARKERS = (
    "# Newman Theta Adaptive Modular C1 Remainder Contract",
    "## Cancellation-Preserving Blocks",
    "## Decaying Remainder",
    "## Direct C1 Contract",
    "## Adaptive Scale",
    "## Open Quantitative Target",
    "sum_(n>=1)b_n(u)=Phi(u)",
    "|J_t-J_(N,t)|<=16*x^4*epsilon_(N,0,m)",
    "N_kappa(x,t)=ceil(max(n_*(5,t,x),n_*(9,t,x)))+kappa",
    "m>=5 and kappa>=0",
    "No new `Q_j`, `Lambda<=0`, RH",
)

CORE_MARKERS = (
    "omega(u)=(1+erf(3sinh(4u)))/2",
    "sum_(n>=1)b_n(u)=Phi(u)",
    "d_(N,j,m)(T)",
    "|J_t-J_(N,t)|",
    "N_kappa(x,t)",
    "modular C1 finite-to-infinite contract",
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


def symbolic_issues() -> list[GateIssue]:
    issues: list[GateIssue] = []
    x, s, sx, r, rx = sp.symbols(
        "x S S_x R R_x", real=True
    )
    full = 16 * x**4 * (s + r)
    partial = 16 * x**4 * s
    if sp.expand(full - partial - 16 * x**4 * r) != 0:
        issues.append(
            finding("symbolic", "value-error", full - partial)
        )
    full_prime = (
        64 * x**3 * (s + r)
        + 16 * x**4 * (sx + rx)
    )
    partial_prime = 64 * x**3 * s + 16 * x**4 * sx
    error = 64 * x**3 * r + 16 * x**4 * rx
    if sp.expand(full_prime - partial_prime - error) != 0:
        issues.append(
            finding(
                "symbolic",
                "derivative-error",
                full_prime - partial_prime - error,
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
    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append(
            finding(
                "rows",
                "id-order",
                [row.get("id") for row in rows],
            )
        )
    if len(rows) != 10:
        issues.append(finding("rows", "count", len(rows)))
    if rows and rows[8].get("readiness") != "not_ready_to_apply":
        issues.append(
            finding(
                "rows",
                "target-promoted",
                rows[8].get("readiness"),
            )
        )
    audit = artifact.get("source_audit", {})
    expected_kinds = {
        "modular_kind": (
            "jensen_window_pf_newman_theta_modular_blend_gate"
        ),
        "high_frequency_kind": (
            "jensen_window_pf_newman_theta_"
            "modular_blend_high_frequency_scout"
        ),
        "adaptive_kind": (
            "jensen_window_pf_newman_theta_"
            "modular_blend_adaptive_saddle_gate"
        ),
        "fixed_block_kind": (
            "jensen_window_pf_newman_theta_"
            "fixed_block_cofinal_obstruction_gate"
        ),
        "compact_kind": (
            "jensen_window_pf_newman_theta_"
            "compact_transversality_interval_certificate"
        ),
    }
    for key, expected in expected_kinds.items():
        if audit.get(key) != expected:
            issues.append(
                finding(
                    "source-audit",
                    key,
                    f"{audit.get(key)!r} != {expected!r}",
                )
            )
    if audit.get("high_frequency_rows") != 10:
        issues.append(
            finding(
                "source-audit",
                "high-frequency-rows",
                audit.get("high_frequency_rows"),
            )
        )
    try:
        cancellation = sp.Float(
            audit["minimum_x200_cancellation_digits"], 80
        )
        if not cancellation > 20:
            issues.append(
                finding(
                    "source-audit",
                    "cancellation",
                    cancellation,
                )
            )
    except (KeyError, TypeError, ValueError) as exc:
        issues.append(
            finding("source-audit", "cancellation", exc)
        )
    expected_counts = [3, 4, 4, 5, 6, 3, 4, 4, 5, 6]
    if audit.get("adaptive_first_counts") != expected_counts:
        issues.append(
            finding(
                "source-audit",
                "adaptive-counts",
                audit.get("adaptive_first_counts"),
            )
        )
    if audit.get("adaptive_collar_covers_all") is not True:
        issues.append(
            finding(
                "source-audit",
                "adaptive-collar",
                audit.get("adaptive_collar_covers_all"),
            )
        )

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not supply",
        "finite cancellation",
        "diagnostic only",
        "no new Q_j",
        "Lambda<=0",
        "RH",
        "Clay-prize",
    ):
        if marker not in boundary:
            issues.append(
                finding("artifact", "proof-boundary", marker)
            )
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
    parser.add_argument(
        "--artifact", type=Path, default=target.DEFAULT_OUT
    )
    parser.add_argument(
        "--note", type=Path, default=target.DEFAULT_NOTE
    )
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    issues = symbolic_issues()
    issues.extend(artifact_issues(artifact))
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
            print(
                f"ISSUE {item.section} "
                f"[{item.issue}] {item.detail}"
            )
        print(
            "validated Newman theta adaptive modular C1 remainder "
            f"contract: 10 rows, {len(issues)} issues, "
            "2 exact modular theorems, "
            "3 exact C1 inequalities/compositions, "
            "1 fixed-block guard, 1 adaptive saddle theorem, "
            "10 cancellation diagnostics, "
            "1 open quantitative target"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
