#!/usr/bin/env python3
"""Validate the modular-tail derivative-budget stress scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_theta_modular_tail_derivative_budget_scout"
DEFAULT_ARTIFACT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
GENERATOR = REPO_ROOT / "work/rh_compute/scripts" / f"{STEM}.py"

EXPECTED_IDS = [
    "ntmtdbs_01_m5_to_m9_envelopes",
    "ntmtdbs_02_envelope_convergence",
    "ntmtdbs_03_transition_jet_ladder",
    "ntmtdbs_04_direct_c1_composition",
    "ntmtdbs_05_fixed_collar_stress",
    "ntmtdbs_06_open_interval_handoff",
]

NOTE_MARKERS = (
    "# Newman Theta Modular-Tail Derivative-Budget Scout",
    "## Derivative Envelopes",
    "## Direct-C1 Stress",
    "## Revised Handoff",
    "kappa=2,3,4",
    "kappa=5",
    "N=max(N_sad,ceil(K*(1+x)^(3/4)))",
    "not a proof of `Lambda<=0` or RH",
)


@dataclass(frozen=True)
class GateIssue:
    section: str
    issue: str
    detail: str


def finding(section: str, issue: str, detail: object) -> GateIssue:
    return GateIssue(section, issue, str(detail))


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def artifact_issues(artifact: dict) -> list[GateIssue]:
    issues: list[GateIssue] = []
    if artifact.get("kind") != STEM:
        issues.append(finding("artifact", "kind", artifact.get("kind")))
    if artifact.get("date") != "2026-07-24":
        issues.append(finding("artifact", "date", artifact.get("date")))
    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append(
            finding("rows", "id-order", [row.get("id") for row in rows])
        )
    if len(rows) != 6:
        issues.append(finding("rows", "count", len(rows)))
    if rows and rows[-1].get("readiness") != "not_ready_to_apply":
        issues.append(
            finding("rows", "handoff-promoted", rows[-1].get("readiness"))
        )

    parameters = artifact.get("parameters", {})
    if parameters.get("orders") != [5, 6, 7, 8, 9]:
        issues.append(
            finding("parameters", "orders", parameters.get("orders"))
        )
    if parameters.get("N_values") != list(range(4, 11)):
        issues.append(
            finding("parameters", "N-values", parameters.get("N_values"))
        )
    if parameters.get("times") != ["0", "0.2"]:
        issues.append(
            finding("parameters", "times", parameters.get("times"))
        )
    if parameters.get("kappas") != [2, 3, 4, 5]:
        issues.append(
            finding("parameters", "kappas", parameters.get("kappas"))
        )

    derivative_rows = artifact.get("derivative_rows", [])
    if len(derivative_rows) != 35:
        issues.append(
            finding("derivatives", "row-count", len(derivative_rows))
        )
    if len(
        {(row.get("N"), row.get("m")) for row in derivative_rows}
    ) != 35:
        issues.append(finding("derivatives", "duplicate-key", "N,m"))
    for row in derivative_rows:
        for field in ("mu0", "mu1", "d0_envelope", "d1_envelope"):
            try:
                if not mp.mpf(row[field]) > 0:
                    issues.append(
                        finding("derivatives", "nonpositive", (row, field))
                    )
            except (KeyError, TypeError, ValueError) as exc:
                issues.append(finding("derivatives", "parse", exc))

    try:
        derivative_delta = mp.mpf(
            artifact["derivative_convergence"][
                "max_relative_coarse_fine_delta"
            ]
        )
        cap_delta = mp.mpf(
            artifact["derivative_convergence"][
                "max_relative_ncap24_ncap28_delta"
            ]
        )
        jet_full = mp.mpf(
            artifact["jet_convergence"]["max_relative_full_jet_delta"]
        )
        jet_partial = mp.mpf(
            artifact["jet_convergence"]["max_relative_partial_jet_delta"]
        )
        if not derivative_delta < mp.mpf("0.03"):
            issues.append(
                finding("convergence", "derivative-node-ladder", derivative_delta)
            )
        if not cap_delta < mp.mpf("1e-12"):
            issues.append(
                finding("convergence", "arithmetic-cap", cap_delta)
            )
        if not jet_full < mp.mpf("1e-20"):
            issues.append(finding("convergence", "full-jet", jet_full))
        if not jet_partial < mp.mpf("1e-20"):
            issues.append(
                finding("convergence", "partial-jet", jet_partial)
            )
    except (KeyError, TypeError, ValueError) as exc:
        issues.append(finding("convergence", "parse", exc))

    jet_rows = artifact.get("jet_rows", [])
    c1_rows = artifact.get("c1_rows", [])
    if len(jet_rows) != 64:
        issues.append(finding("jets", "row-count", len(jet_rows)))
    if len(c1_rows) != 64:
        issues.append(finding("C1", "row-count", len(c1_rows)))
    for row in c1_rows:
        order_rows = row.get("orders", [])
        if [item.get("m") for item in order_rows] != [5, 6, 7, 8, 9]:
            issues.append(
                finding("C1", "order-ladder", (row.get("t"), row.get("x")))
            )
            continue
        passing = [
            item["m"]
            for item in order_rows
            if item.get("diagnostic_disjunction_passes") is True
        ]
        expected_first = passing[0] if passing else None
        if row.get("first_passing_m_5_to_9") != expected_first:
            issues.append(
                finding(
                    "C1",
                    "first-passing-order",
                    (
                        row.get("t"),
                        row.get("x"),
                        row.get("kappa"),
                        row.get("first_passing_m_5_to_9"),
                        expected_first,
                    ),
                )
            )
        for item in order_rows:
            try:
                recomputed = (
                    mp.mpf(item["value_ratio"]) > 1
                    or mp.mpf(item["derivative_ratio"]) > 1
                )
                if recomputed != item["diagnostic_disjunction_passes"]:
                    issues.append(
                        finding("C1", "ratio-boolean", (row, item))
                    )
            except (KeyError, TypeError, ValueError) as exc:
                issues.append(finding("C1", "parse", exc))

    summary = artifact.get("summary", {})
    if summary.get("x200_kappa2_first_passing_orders") != [9, 8]:
        issues.append(
            finding(
                "summary",
                "x200-kappa2",
                summary.get("x200_kappa2_first_passing_orders"),
            )
        )
    if summary.get("x300_kappa2_to_4_all_fail_m5_to_m9") is not True:
        issues.append(
            finding("summary", "x300-failure-pattern", summary)
        )
    if summary.get("x300_kappa5_all_pass") is not True:
        issues.append(finding("summary", "x300-kappa5", summary))

    provenance = artifact.get("provenance", {})
    if provenance.get("generator_sha256") != file_hash(GENERATOR):
        issues.append(
            finding(
                "provenance",
                "generator-hash",
                provenance.get("generator_sha256"),
            )
        )
    source_path = REPO_ROOT / provenance.get("exact_envelope_result", "")
    if not source_path.is_file():
        issues.append(finding("provenance", "source-missing", source_path))
    elif provenance.get("exact_envelope_result_sha256") != file_hash(source_path):
        issues.append(
            finding(
                "provenance",
                "source-hash",
                provenance.get("exact_envelope_result_sha256"),
            )
        )

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "diagnostic only",
        "does not supply",
        "fixed-collar obstruction",
        "first-jet separation",
        "Lambda<=0",
        "RH",
        "Clay-prize",
    ):
        if marker not in boundary:
            issues.append(finding("artifact", "proof-boundary", marker))
    return issues


def note_issues(path: Path) -> list[GateIssue]:
    text = path.read_text(encoding="utf-8")
    return [
        finding("note", "missing-marker", marker)
        for marker in NOTE_MARKERS
        if marker not in text
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--artifact", type=Path, default=DEFAULT_ARTIFACT
    )
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    issues = artifact_issues(artifact)
    issues.extend(note_issues(args.note))
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
            "validated Newman theta modular-tail derivative-budget scout: "
            f"35 derivative rows, 64 C1 rows, {len(issues)} issues, "
            "2 node ladders, 1 arithmetic-cap stress, "
            "1 fixed-collar non-promotion pattern"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
