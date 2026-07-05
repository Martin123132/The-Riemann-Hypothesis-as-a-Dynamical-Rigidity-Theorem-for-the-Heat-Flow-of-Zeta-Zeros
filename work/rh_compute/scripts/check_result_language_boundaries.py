#!/usr/bin/env python3
"""Scan output notes for obvious result-language overclaims.

This is a conservative proof-safety linter.  It is not a mathematical
verifier.  It catches unqualified global-proof language in prose while
ignoring fenced examples such as "Forbidden statement" and "Blocked proof
step" snippets.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, asdict
import json
from pathlib import Path
import re


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTPUTS = REPO_ROOT / "outputs"


@dataclass(frozen=True)
class LanguageRule:
    name: str
    pattern: re.Pattern[str]


@dataclass(frozen=True)
class Finding:
    path: str
    line: int
    rule: str
    text: str


RULES: tuple[LanguageRule, ...] = (
    LanguageRule(
        "completed-rh-proof",
        re.compile(r"\b(?:we\s+)?(?:have\s+)?(?:proved|prove|proves|proven)\s+(?:the\s+)?(?:RH|Riemann Hypothesis)\b", re.I),
    ),
    LanguageRule(
        "proof-of-rh",
        re.compile(r"\b(?:proof|proved|proves|proven)\s+of\s+(?:the\s+)?(?:RH|Riemann Hypothesis)\b", re.I),
    ),
    LanguageRule(
        "completed-newman-direction",
        re.compile(r"\b(?:proved|prove|proves|proven)\s+.*Lambda\s*<=\s*0\b", re.I),
    ),
    LanguageRule(
        "finite-data-proves-global",
        re.compile(
            r"\b(?:finite|many|numerical|certified)\b.*\b(?:proves|prove|proved|implies|establishes)\b.*\b(?:RH|Lambda\s*<=\s*0|PF-infinity|Laguerre-Polya|all-order)\b",
            re.I,
        ),
    ),
    LanguageRule(
        "signed-hankel-implies-rh",
        re.compile(r"\bsigned[- ]Hankel\b.*\bimplies\b.*\b(?:RH|Lambda\s*<=\s*0)\b", re.I),
    ),
    LanguageRule(
        "jensen-passes-imply-rh",
        re.compile(r"\bJensen\b.*\b(?:passes|tests|asymptotics)\b.*\b(?:imply|prove|proves)\b.*\b(?:RH|Lambda\s*<=\s*0)\b", re.I),
    ),
    LanguageRule(
        "clay-ready-proof",
        re.compile(r"\bClay-ready proof\b", re.I),
    ),
    LanguageRule(
        "pf-infinity-proved",
        re.compile(r"\bPF-infinity\b.*\b(?:is|has been)\s+(?:proved|proven|established)\b", re.I),
    ),
)


ALLOW_MARKERS: tuple[str, ...] = (
    "not ",
    "not a ",
    "not yet",
    "no ",
    "does not",
    "do not",
    "cannot",
    "can't",
    "could ",
    "would ",
    "must ",
    "need ",
    "needs ",
    "needed ",
    "if ",
    "unless ",
    "equivalent",
    "lambda >= 0",
    "complementary newman",
    "conditional",
    "target",
    "route",
    "path",
    "toward",
    "missing",
    "blocked proof step",
    "forbidden",
    "do not use",
    "no manuscript claim",
    "remaining",
    "not enough",
    "not claim",
    "not claimed",
)


def markdown_files(root: Path) -> list[Path]:
    return sorted(path for path in root.glob("*.md") if path.is_file())


def is_allowed_context(line: str, previous_heading: str) -> bool:
    lowered = line.lower()
    heading = previous_heading.lower()
    if any(marker in lowered for marker in ALLOW_MARKERS):
        return True
    if any(marker in heading for marker in ("forbidden", "blocked proof step", "do not use", "boundary", "target")):
        return True
    return False


def scan_file(path: Path, root: Path) -> list[Finding]:
    findings: list[Finding] = []
    in_fence = False
    previous_heading = ""
    with path.open("r", encoding="utf-8") as handle:
        for lineno, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            stripped = line.strip()
            if stripped.startswith("```"):
                in_fence = not in_fence
                continue
            if in_fence:
                continue
            if stripped.startswith("#") or stripped.endswith(":"):
                previous_heading = stripped
            if stripped.startswith(">"):
                continue
            if not stripped:
                continue
            for rule in RULES:
                if not rule.pattern.search(line):
                    continue
                if is_allowed_context(line, previous_heading):
                    continue
                findings.append(
                    Finding(
                        path=str(path.relative_to(root)),
                        line=lineno,
                        rule=rule.name,
                        text=stripped,
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
    findings: list[Finding] = []
    for path in files:
        findings.extend(scan_file(path, args.outputs_dir))

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
            print(f"OVERCLAIM {finding.path}:{finding.line} [{finding.rule}] {finding.text}")
        print(
            "validated result-language boundaries: "
            f"scanned {len(files)} markdown files, {len(findings)} overclaims"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
