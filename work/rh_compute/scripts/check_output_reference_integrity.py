#!/usr/bin/env python3
"""Validate path references in output markdown notes.

The proof-programme notes are only useful if their reproducibility links point
to real scripts, notes, and promoted result artifacts.  This scanner checks
concrete references to scripts, output notes, and result files.  Planned
deliverables in the compute plan are allowed to be absent but are reported.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import re


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTPUTS = REPO_ROOT / "outputs"


REF_PATTERN = re.compile(
    r"(?P<path>"
    r"work/rh_compute/scripts/[A-Za-z0-9_./-]+\.py"
    r"|work/rh_compute/results/[A-Za-z0-9_./-]+\.(?:jsonl(?:\.gz)?|json|csv|md|log)"
    r"|outputs/[A-Za-z0-9_./-]+\.md"
    r")"
)


@dataclass(frozen=True)
class Reference:
    source: str
    line: int
    path: str
    planned: bool


@dataclass(frozen=True)
class ReferenceIssue:
    source: str
    line: int
    path: str
    issue: str


def markdown_files(root: Path) -> list[Path]:
    return sorted(path for path in root.glob("*.md") if path.is_file())


def is_planned_deliverable(lines: list[str], index: int, source: Path) -> bool:
    if source.name != "Rigorous_Numerics_Compute_Plan.md":
        return False
    window = lines[max(0, index - 5) : index + 1]
    return any(line.strip().lower() == "deliverable:" for line in window)


def extract_references(path: Path, root: Path) -> list[Reference]:
    lines = path.read_text(encoding="utf-8").splitlines()
    refs: list[Reference] = []
    for index, line in enumerate(lines):
        for match in REF_PATTERN.finditer(line):
            ref = match.group("path").rstrip(".,:;)")
            refs.append(
                Reference(
                    source=str(path.relative_to(root)),
                    line=index + 1,
                    path=ref,
                    planned=is_planned_deliverable(lines, index, path),
                )
            )
    return refs


def validate_references(refs: list[Reference]) -> tuple[list[ReferenceIssue], list[Reference]]:
    issues: list[ReferenceIssue] = []
    planned_missing: list[Reference] = []
    seen: set[tuple[str, str, bool]] = set()
    for ref in refs:
        key = (ref.source, ref.path, ref.planned)
        if key in seen:
            continue
        seen.add(key)
        target = REPO_ROOT / ref.path
        if target.exists():
            continue
        if ref.planned:
            planned_missing.append(ref)
            continue
        issues.append(
            ReferenceIssue(
                source=ref.source,
                line=ref.line,
                path=ref.path,
                issue="missing-referenced-path",
            )
        )
    return issues, planned_missing


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outputs-dir", type=Path, default=DEFAULT_OUTPUTS)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    files = markdown_files(args.outputs_dir)
    refs: list[Reference] = []
    for path in files:
        refs.extend(extract_references(path, args.outputs_dir))
    issues, planned_missing = validate_references(refs)
    ok = not issues

    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "files_scanned": len(files),
                    "references_scanned": len(refs),
                    "planned_missing": [asdict(ref) for ref in planned_missing],
                    "issues": [asdict(issue) for issue in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for issue in issues:
            print(f"REFERENCE {issue.source}:{issue.line} [{issue.issue}] {issue.path}")
        for ref in planned_missing:
            print(f"PLANNED-MISSING {ref.source}:{ref.line} {ref.path}")
        print(
            "validated output references: "
            f"scanned {len(files)} markdown files, {len(refs)} path references, "
            f"{len(issues)} missing required paths, {len(planned_missing)} planned missing deliverables"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
