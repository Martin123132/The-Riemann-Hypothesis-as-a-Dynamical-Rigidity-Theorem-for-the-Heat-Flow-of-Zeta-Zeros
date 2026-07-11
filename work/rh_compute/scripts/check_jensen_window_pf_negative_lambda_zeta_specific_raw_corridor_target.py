#!/usr/bin/env python3
"""Validate the zeta-specific raw-corridor theorem target for negative lambda."""

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

from jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlzrct_01_raw_corridor_statement",
    "nlzrct_02_exact_equivalence_to_adaptive_cone",
    "nlzrct_03_k200_finite_anchor",
    "nlzrct_04_generic_stieltjes_shortcut_blocked",
    "nlzrct_05_positive_gaussian_mixture_shortcut_blocked",
    "nlzrct_06_signed_gaussian_local_route",
    "nlzrct_07_two_scale_remainder_requirement",
    "nlzrct_08_ratio_recurrence_route",
    "nlzrct_09_acceptance_gate",
}

ALLOWED_ROLES = {
    "open_statement",
    "exact_reduction",
    "finite_anchor",
    "rejected_shortcut",
    "live_route",
    "open_requirement",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Zeta-Specific Raw-Corridor Target",
    "Status: open theorem target",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.py",
    "validated Jensen-window PF negative-lambda zeta-specific raw-corridor target: 9 rows, 0 issues, 2 live routes, 2 rejected shortcuts, 0 ready-to-apply rows",
    "M_k = mu_{2k}",
    "R_k = M_(k+1)*M_(k-1)/M_k^2",
    "1 <= R_k <= (2*k+1)/(2*k-1)",
    "((2*k-1)*(2*k+3)/(2*k+1)^2)*R_k <= R_(k+1) <= (2+(2*k-1)*R_k)/(2*k+1)",
    "k200 raw-cone rows: 597 / 597",
    "k200 corridor rows: 594 / 594",
    "generic Stieltjes/raw-log-convexity proof",
    "positive Gaussian scale-mixture proof of the upper wall",
    "signed Gaussian perturbation plus two-scale uniform remainders",
    "direct zeta-specific ratio recurrence compatible with increasing scaled defect",
    "outputs/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.md",
    "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
    "outputs/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
)


@dataclass(frozen=True)
class RawCorridorIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> RawCorridorIssue:
    return RawCorridorIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[RawCorridorIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_symbolic_identities() -> list[RawCorridorIssue]:
    k, r, rn = sp.symbols("k r rn")
    x = (2 * k - 1) * r / (2 * k + 1)
    xn = (2 * k + 1) * rn / (2 * k + 3)
    s = (2 * k + 1 - (2 * k - 1) * r) / 2
    sn = (2 * k + 3 - (2 * k + 1) * rn) / 2
    lower = (2 * k - 1) * (2 * k + 3) * r / (2 * k + 1) ** 2
    upper = (2 + (2 * k - 1) * r) / (2 * k + 1)
    checks = {
        "upper-wall-to-defect-nonnegative": sp.simplify((1 - x) * (2 * k + 1) / 2 - s),
        "lower-wall-to-scaled-upper": sp.simplify((1 - s) - ((2 * k - 1) * (r - 1) / 2)),
        "monotone-bridge-lower-corridor": sp.simplify(
            (xn - x)
            - (((2 * k + 1) ** 2 * rn - (2 * k - 1) * (2 * k + 3) * r) / ((2 * k + 1) * (2 * k + 3)))
        ),
        "scaled-increase-upper-corridor": sp.simplify((sn - s) - ((2 + (2 * k - 1) * r - (2 * k + 1) * rn) / 2)),
        "corridor-width": sp.simplify((upper - lower) - (2 * (2 * k + 1 - (2 * k - 1) * r) / (2 * k + 1) ** 2)),
    }
    issues: list[RawCorridorIssue] = []
    for name, value in checks.items():
        if sp.simplify(value) != 0:
            issues.append(issue("symbolic", f"bad-{name}", str(value)))
    return issues


def validate_top_level(artifact: dict) -> list[RawCorridorIssue]:
    issues: list[RawCorridorIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "open_theorem_target":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    if artifact.get("target_id") != "target_negative_lambda_zeta_specific_raw_corridor":
        issues.append(issue("<artifact>", "bad-target-id", repr(artifact.get("target_id"))))
    for key in (
        "source_raw_bridge",
        "source_raw_obstruction",
        "source_adaptive_target",
        "source_uniform_remainder",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in ("open theorem target", "zeta-specific", "does not prove", "cone entry", "jwpf_06", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[RawCorridorIssue]:
    recomputed = build_artifact()
    issues: list[RawCorridorIssue] = []
    for key in ("target_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[RawCorridorIssue], int, int, int, int]:
    rows = artifact.get("target_rows", [])
    issues: list[RawCorridorIssue] = []
    if not isinstance(rows, list):
        return [issue("target_rows", "bad-rows", repr(type(rows)))], 0, 0, 0, 0
    by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(by_id)):
        issues.append(issue(missing, "missing-row", missing))
    live_routes = 0
    rejected_shortcuts = 0
    ready_to_apply = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("target_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "claim_if_proved", "source_artifacts", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if row.get("role") not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(row.get("role"))))
        readiness = row.get("readiness")
        if row.get("role") == "exact_reduction":
            if readiness != "available_exact":
                issues.append(issue(row_id, "bad-exact-readiness", repr(readiness)))
        elif readiness != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-open-readiness", repr(readiness)))
        if readiness == "ready_to_apply":
            ready_to_apply += 1
        if row.get("role") == "live_route":
            live_routes += 1
        if row.get("role") == "rejected_shortcut":
            rejected_shortcuts += 1
            text = f"{row.get('gap', '')} {row.get('proof_boundary', '')}".lower()
            if not any(marker in text for marker in ("reject", "blocked", "wrong direction", "witness")):
                issues.append(issue(row_id, "weak-rejection-language", text))
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("open", "not", "only", "finite", "rejected", "live")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    return issues, len(rows), live_routes, rejected_shortcuts, ready_to_apply


def validate_summary(
    artifact: dict,
    row_count: int,
    live_routes: int,
    rejected_shortcuts: int,
    ready_to_apply: int,
) -> list[RawCorridorIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "target_rows": 9,
        "exact_available_rows": 1,
        "finite_anchor_rows": 1,
        "live_routes": 2,
        "open_requirement_rows": 1,
        "rejected_shortcut_rows": 2,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "finite_raw_cone_rows": 597,
        "finite_corridor_rows": 594,
        "open_theorem_target": True,
    }
    issues: list[RawCorridorIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 9:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if live_routes != 2:
        issues.append(issue("summary", "bad-live-route-count", str(live_routes)))
    if rejected_shortcuts != 2:
        issues.append(issue("summary", "bad-rejected-shortcut-count", str(rejected_shortcuts)))
    if ready_to_apply != 0:
        issues.append(issue("summary", "bad-ready-to-apply-count", str(ready_to_apply)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "zeta-specific raw-corridor",
        "upper raw wall",
        "adaptive corridor",
        "597 raw-cone",
        "594 corridor",
        "generic stieltjes",
        "positive gaussian",
        "signed gaussian",
        "two-scale",
        "ratio recurrence",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("no row", "open_theorem_target", "stieltjes", "gaussian", "finite", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[RawCorridorIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[RawCorridorIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "raw-corridor theorem is proved",
        "cone entry is proved",
        "jwpf_06 is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[RawCorridorIssue], dict]:
    artifact = load_json(target_path)
    issues: list[RawCorridorIssue] = []
    issues.extend(validate_symbolic_identities())
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, live_routes, rejected_shortcuts, ready_to_apply = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, live_routes, rejected_shortcuts, ready_to_apply))
    issues.extend(validate_note(note_path))
    return issues, artifact.get("summary", {})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    target = args.target if args.target.is_absolute() else REPO_ROOT / args.target
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    issues, summary = validate(target, note)
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "summary": summary, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-ZETA-RAW-CORRIDOR {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda zeta-specific raw-corridor target: "
            f"{summary.get('target_rows')} rows, {len(issues)} issues, "
            f"{summary.get('live_routes')} live routes, "
            f"{summary.get('rejected_shortcut_rows')} rejected shortcuts, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
