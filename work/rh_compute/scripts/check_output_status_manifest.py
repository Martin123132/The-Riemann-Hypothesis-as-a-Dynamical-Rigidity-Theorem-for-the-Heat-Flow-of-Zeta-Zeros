#!/usr/bin/env python3
"""Validate output-note status classification and proof-boundary metadata.

Every referee-facing markdown artifact should identify what kind of object it
is near the top of the file: certificate, diagnostic, theorem target, roadmap,
countermodel library, progress report, and so on.  This scanner enforces that
discipline without judging the mathematics inside the note.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTPUTS = REPO_ROOT / "outputs"


ARTIFACT_KEYWORDS: tuple[str, ...] = (
    "artifact",
    "audit",
    "certificate",
    "correction",
    "diagnostic",
    "gate",
    "inventory",
    "ledger",
    "map",
    "note",
    "plan",
    "report",
    "roadmap",
    "scout",
    "status",
    "target",
)


BOUNDARY_MARKERS: tuple[str, ...] = (
    "not a proof",
    "not proof",
    "not evidence",
    "not yet",
    "not currently proved",
    "not a claim",
    "does not promote",
    "not promoted",
    "not a replacement",
    "not an all-order",
    "not a rigorous",
    "not a clay-ready",
)


@dataclass(frozen=True)
class StatusFinding:
    path: str
    issue: str
    detail: str


def markdown_files(root: Path) -> list[Path]:
    return sorted(path for path in root.glob("*.md") if path.is_file())


def status_block(lines: list[str], start: int) -> str:
    parts = [lines[start].strip()]
    for line in lines[start + 1 :]:
        stripped = line.strip()
        if not stripped:
            break
        if stripped.startswith("#") or stripped.startswith("```"):
            break
        parts.append(stripped)
    return " ".join(parts)


def validate_file(path: Path, root: Path) -> list[StatusFinding]:
    rel = str(path.relative_to(root))
    lines = path.read_text(encoding="utf-8").splitlines()
    findings: list[StatusFinding] = []

    if not lines or not lines[0].startswith("# "):
        findings.append(StatusFinding(rel, "missing-h1", "first line must be a markdown H1"))

    top = lines[:14]
    date_lines = [line for line in top if line.startswith("Date:")]
    if not date_lines:
        findings.append(StatusFinding(rel, "missing-date", "Date: line must appear near the top"))

    status_indices = [i for i, line in enumerate(top) if line.startswith("Status:")]
    if not status_indices:
        findings.append(StatusFinding(rel, "missing-status", "Status: line must appear near the top"))
        return findings

    block = status_block(lines, status_indices[0])
    lowered = block.lower()
    if not any(keyword in lowered for keyword in ARTIFACT_KEYWORDS):
        findings.append(
            StatusFinding(
                rel,
                "unclear-artifact-kind",
                "Status must identify an artifact kind such as certificate, diagnostic, ledger, note, report, or roadmap",
            )
        )
    if not any(marker in lowered for marker in BOUNDARY_MARKERS):
        findings.append(
            StatusFinding(
                rel,
                "missing-proof-boundary",
                "Status must contain proof-boundary language such as 'not a proof' or 'not a proof artifact'",
            )
        )
    return findings


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outputs-dir", type=Path, default=DEFAULT_OUTPUTS)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    files = markdown_files(args.outputs_dir)
    findings: list[StatusFinding] = []
    for path in files:
        findings.extend(validate_file(path, args.outputs_dir))

    ok = not findings
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "files_scanned": len(files),
                    "findings": [asdict(finding) for finding in findings],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for finding in findings:
            print(f"STATUS {finding.path} [{finding.issue}] {finding.detail}")
        print(
            "validated output artifact statuses: "
            f"scanned {len(files)} markdown files, {len(findings)} status issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
