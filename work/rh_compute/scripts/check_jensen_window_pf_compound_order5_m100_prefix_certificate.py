#!/usr/bin/env python3
"""Validate the rigorous lambda=-100 order-five prefix certificate."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = SCRIPT_DIR.parent / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

from jensen_window_pf_compound_order5_m100_prefix_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    MAX_COEFFICIENT_INDEX,
    ORDER5_REDUCTION,
    PREFIX_LAST_N,
    PRECISION_BITS,
    REPO_ROOT,
    SOURCE_PATHS,
    build_artifact,
    sha256,
)


@dataclass(frozen=True)
class PrefixIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> PrefixIssue:
    return PrefixIssue(section, code, str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_artifact(path: Path) -> tuple[dict, list[PrefixIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = load_json(path)
    issues: list[PrefixIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order5_m100_prefix_certificate"
    ):
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != (
        "rigorous finite lambda=-100 contiguous order-five prefix through n=316"
    ):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))

    summary = artifact.get("summary", {})
    expected_summary = {
        "rows": 7,
        "ready_to_apply_rows": 6,
        "coefficients": 325,
        "positive_J_rows": 317,
        "positive_relative_rows": 317,
        "positive_H5_rows": 317,
        "inconclusive_rows": 0,
        "open_analytic_tails": 1,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    exact = artifact.get("exact", {})
    if exact.get("symbolic_residual") != "0":
        issues.append(issue("exact", "nonzero-residual", exact.get("symbolic_residual")))
    if exact.get("factorization") != "H_(5,n)=W_n*J_n":
        issues.append(issue("exact", "bad-factorization", exact.get("factorization")))
    if exact.get("scale_positive_in_lower_cone") is not True:
        issues.append(issue("exact", "scale-not-positive", exact))

    source_rows = artifact.get("source_diagnostics", [])
    if len(source_rows) != len(SOURCE_PATHS):
        issues.append(issue("sources", "bad-count", len(source_rows)))
    for source, recorded in zip(SOURCE_PATHS, source_rows):
        if not source.exists():
            issues.append(issue("sources", "missing-source", source))
            continue
        relative = source.relative_to(REPO_ROOT).as_posix()
        if recorded.get("source") != relative:
            issues.append(issue("sources", "bad-path", recorded))
        if recorded.get("sha256") != sha256(source):
            issues.append(issue("sources", "hash-mismatch", relative))
    if not ORDER5_REDUCTION.exists():
        issues.append(issue("sources", "missing-reduction", ORDER5_REDUCTION))

    finite = artifact.get("finite", {})
    if finite.get("n_range") != [0, PREFIX_LAST_N]:
        issues.append(issue("finite", "bad-n-range", finite.get("n_range")))
    if finite.get("coefficient_range") != [0, MAX_COEFFICIENT_INDEX]:
        issues.append(
            issue("finite", "bad-coefficient-range", finite.get("coefficient_range"))
        )
    if finite.get("precision_bits") != PRECISION_BITS:
        issues.append(issue("finite", "bad-precision", finite.get("precision_bits")))
    if finite.get("minimum_J_n") != PREFIX_LAST_N:
        issues.append(issue("finite", "bad-min-J-index", finite.get("minimum_J_n")))
    if finite.get("minimum_relative_n") != PREFIX_LAST_N:
        issues.append(
            issue(
                "finite",
                "bad-min-relative-index",
                finite.get("minimum_relative_n"),
            )
        )

    flint.ctx.prec = PRECISION_BITS
    finite_rows = finite.get("rows", [])
    if len(finite_rows) != PREFIX_LAST_N + 1:
        issues.append(issue("finite", "bad-row-count", len(finite_rows)))
    for expected_n, row in enumerate(finite_rows):
        if row.get("n") != expected_n:
            issues.append(issue("finite", "bad-index", row))
            continue
        try:
            stable = flint.arb(row.get("J_ball", "0"))
            relative = flint.arb(row.get("relative_ball", "0"))
            stable_lower = Decimal(row.get("J_lower", "0"))
            relative_lower = Decimal(row.get("relative_lower", "0"))
        except Exception as exc:
            issues.append(issue("finite", "ball-parse", (expected_n, exc)))
            continue
        if not stable.is_finite() or stable_lower <= 0:
            issues.append(issue("finite", "nonpositive-J-lower", expected_n))
        if not relative.is_finite() or relative_lower <= 0:
            issues.append(issue("finite", "nonpositive-relative-lower", expected_n))
        if row.get("H5_sign") != "positive_by_exact_scale":
            issues.append(issue("finite", "bad-H5-sign", row))

    try:
        minimum_relative = flint.arb(finite.get("minimum_relative_ball", "0"))
    except Exception as exc:
        issues.append(issue("finite", "minimum-relative-parse", exc))
    else:
        if not bool(minimum_relative > flint.arb("0.0062")):
            issues.append(
                issue("finite", "weak-minimum-relative", minimum_relative.str(30))
            )

    rows = artifact.get("rows", [])
    if len(rows) != 7:
        issues.append(issue("rows", "bad-count", len(rows)))
    ids = [row.get("id") for row in rows]
    if len(ids) != len(set(ids)):
        issues.append(issue("rows", "duplicate-id", ids))
    open_rows = [row for row in rows if row.get("readiness") == "not_ready_to_apply"]
    if len(open_rows) != 1 or "n>=317" not in open_rows[0].get("formula", ""):
        issues.append(issue("rows", "bad-open-tail", open_rows))

    try:
        canonical = build_artifact()
    except Exception as exc:  # pragma: no cover - diagnostic path
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[PrefixIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous finite `lambda=-100` contiguous order-five prefix",
        "This is not a proof of the analytic tail",
        "H_(5,n)=W_n*J_n",
        "A_k(-100)>0 for every 0<=k<=324",
        "J_n(-100)>0 for every 0<=n<=316",
        "H_(5,n)(-100)>0 for every 0<=n<=316",
        "minimum J row: n=316",
        "minimum relative row: n=316",
        "0.006269",
        "J_n(-100)>0 for every n>=317",
        "analytic tail remains open",
        "outputs/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.md",
        "outputs/signed_hankel_jensen_bridge_target.md",
        "outputs/formal_core.md",
    )
    issues: list[PrefixIssue] = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "all-shift order five is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-overclaim", forbidden))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    if issues:
        for finding in issues:
            print(
                f"ORDER5-M100-PREFIX {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-five lambda=-100 prefix: {len(issues)} issues")
        return 1

    summary = artifact["summary"]
    print(
        "validated order-five lambda=-100 prefix: "
        f"{summary['coefficients']} coefficients, "
        f"{summary['positive_J_rows']} positive J5 rows, "
        f"{summary['positive_relative_rows']} positive relative margins, "
        f"{summary['positive_H5_rows']} positive H5 signs, "
        f"{summary['inconclusive_rows']} inconclusive, "
        f"{summary['open_analytic_tails']} open analytic tail"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
