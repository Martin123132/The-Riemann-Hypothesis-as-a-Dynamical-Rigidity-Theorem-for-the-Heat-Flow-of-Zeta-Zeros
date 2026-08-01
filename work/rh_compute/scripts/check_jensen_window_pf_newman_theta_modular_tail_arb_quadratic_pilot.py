#!/usr/bin/env python3
"""Validate the rigorous Arb modular-tail quadratic pilot."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_theta_modular_tail_arb_quadratic_pilot as target  # noqa: E402

from flint import arb  # noqa: E402


EXPECTED_IDS = [
    "ntmtaqp_01_weighted_cauchy_schwarz",
    "ntmtaqp_02_heat_envelope_application",
    "ntmtaqp_03_arb_finite_compact_pilot",
    "ntmtaqp_04_omitted_tail_guard",
    "ntmtaqp_05_matrix_handoff",
]
NOTE_MARKERS = (
    "# Newman Theta Modular-Tail Arb Quadratic Pilot",
    "## Exact Weighted Reduction",
    "## Rigorous Finite Compact Pilot",
    "## Non-Promotion Guard",
    "## Next Certificate",
    "n=7..16",
    "n>=17",
    "No strict Laguerre conclusion",
)
CORE_MARKERS = (
    "weighted Cauchy-Schwarz reduction",
    "Arb quadratic pilot",
    "n=7,...,16",
    "outer-u tail",
)


@dataclass(frozen=True)
class PilotIssue:
    section: str
    issue: str
    detail: str


def finding(section: str, issue: str, detail: object) -> PilotIssue:
    return PilotIssue(section, issue, str(detail))


def artifact_issues(artifact: dict, recomputed: dict) -> list[PilotIssue]:
    issues: list[PilotIssue] = []
    if artifact.get("kind") != target.STEM:
        issues.append(finding("artifact", "kind", artifact.get("kind")))
    if artifact.get("date") != "2026-07-24":
        issues.append(finding("artifact", "date", artifact.get("date")))
    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append(
            finding("rows", "id-order", [row.get("id") for row in rows])
        )
    if len(rows) != 5:
        issues.append(finding("rows", "count", len(rows)))
    if rows and rows[-1].get("readiness") != "not_ready_to_apply":
        issues.append(
            finding("rows", "handoff-promoted", rows[-1].get("readiness"))
        )

    try:
        stored_ball = arb(artifact["pilot"]["real_enclosure"])
        recomputed_ball = arb(recomputed["pilot"]["real_enclosure"])
        if not stored_ball.overlaps(recomputed_ball):
            issues.append(
                finding(
                    "pilot",
                    "recompute-disjoint",
                    f"{stored_ball} versus {recomputed_ball}",
                )
            )
        if stored_ball.lower() <= 0:
            issues.append(finding("pilot", "not-positive", stored_ball))
        if stored_ball.rel_accuracy_bits() < 20:
            issues.append(
                finding(
                    "pilot",
                    "low-accuracy",
                    stored_ball.rel_accuracy_bits(),
                )
            )
    except (KeyError, TypeError, ValueError) as exc:
        issues.append(finding("pilot", "enclosure", exc))

    pilot = artifact.get("pilot", {})
    expected = {
        "precision_bits": 96,
        "derivative_order": 9,
        "n_start": 7,
        "n_stop": 16,
        "t_cap_exact": "1/5",
        "u_interval_exact": ["0", "11/5"],
        "integrand_entire": True,
    }
    for key, value in expected.items():
        if pilot.get(key) != value:
            issues.append(finding("pilot", key, pilot.get(key)))

    source = artifact.get("source", {})
    expected_hash = sha256(target.ENVELOPE_SOURCE.read_bytes()).hexdigest()
    if source.get("derivative_envelope_sha256") != expected_hash:
        issues.append(
            finding(
                "source",
                "derivative-envelope-hash",
                source.get("derivative_envelope_sha256"),
            )
        )

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "infinite arithmetic tail",
        "outer-u tail",
        "first-jet lower separation",
        "Lambda<=0",
        "RH",
        "Clay-prize",
    ):
        if marker not in boundary:
            issues.append(finding("artifact", "proof-boundary", marker))
    return issues


def text_issues(
    path: Path, markers: tuple[str, ...], section: str
) -> list[PilotIssue]:
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
    recomputed = target.build_artifact()
    issues = artifact_issues(artifact, recomputed)
    issues.extend(text_issues(args.note, NOTE_MARKERS, "note"))
    issues.extend(
        text_issues(
            target.REPO_ROOT / "outputs" / "formal_core.md",
            CORE_MARKERS,
            "formal-core",
        )
    )
    ok = not issues
    payload = {
        "ok": ok,
        "artifact": str(args.artifact),
        "rows": len(artifact.get("rows", [])),
        "issues": [asdict(issue) for issue in issues],
    }
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(
            "Newman theta modular-tail Arb quadratic pilot checker: "
            f"ok={ok} rows={payload['rows']} issues={len(issues)}"
        )
        for issue in issues:
            print(f"- [{issue.section}] {issue.issue}: {issue.detail}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
