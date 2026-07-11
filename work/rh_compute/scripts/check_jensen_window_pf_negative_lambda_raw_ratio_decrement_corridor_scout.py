#!/usr/bin/env python3
"""Validate the negative-lambda raw-ratio decrement-corridor scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
    exact_counterexample,
)


REQUIRED_ROW_IDS = {
    "nlrdc_01_exact_decrement_corridor",
    "nlrdc_02_exact_width_equivalence",
    "nlrdc_03_k200_decrement_corridor",
    "nlrdc_04_theta_coordinate_shape",
    "nlrdc_05_monotone_decrease_shortcut_blocked",
    "nlrdc_06_overfast_decrease_shortcut_blocked",
    "nlrdc_07_live_theta_recurrence_route",
    "nlrdc_08_open_requirements",
    "nlrdc_09_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_reduction",
    "finite_pattern",
    "exact_counterexample",
    "live_route",
    "open_requirement",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Raw-Ratio Decrement-Corridor Scout",
    "Status: exact finite theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.json",
    "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k200.jsonl",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.py",
    "validated Jensen-window PF negative-lambda raw-ratio decrement-corridor scout: 9 rows, 0 issues, 594 decrement-corridor rows, 591 theta-k-monotone rows, 2 exact counterexamples, 0 ready-to-apply rows",
    "M_k = mu_{2k}",
    "R_k = M_(k+1)*M_(k-1)/M_k^2",
    "D_k = R_k - R_(k+1)",
    "2*(R_k-1)/(2*k+1) <= D_k <= 4*R_k/(2*k+1)^2",
    "4*R_k/(2*k+1)^2 - 2*(R_k-1)/(2*k+1)",
    "= 2*(2*k+1-(2*k-1)*R_k)/(2*k+1)^2",
    "raw decrease rows: 594 / 594",
    "lower decrement rows: 594 / 594",
    "upper decrement rows: 594 / 594",
    "decrement-corridor rows: 594 / 594",
    "theta unit rows: 594 / 594",
    "theta-k monotone rows: 591 / 591",
    "theta lambda-order rows: 396 / 396",
    "min lower decrement margin: 2.241640098067743212E-6 at lambda=-25.0, k=198",
    "min upper decrement margin: 7.465690906027845625E-6 at lambda=-100.0, k=198",
    "theta range: 1.644126617921189501E-1 at lambda=-25.0, k=198 to 9.087309076587821269E-1 at lambda=-100.0, k=1",
    "R_1=2, R_2=3/2: raw cone and decrease hold, but the lower decrement wall fails",
    "R_1=2, R_2=1: raw cone and decrease hold, but the upper decrement wall fails",
    "outputs/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.md",
    "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
    "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
    "outputs/jensen_window_pf_negative_lambda_adaptive_envelope_obligations.md",
    "outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md",
)


@dataclass(frozen=True)
class DecrementIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> DecrementIssue:
    return DecrementIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[DecrementIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_symbolic_identities() -> list[DecrementIssue]:
    k, r, rn = sp.symbols("k r rn")
    decrement = r - rn
    dec_lower = 2 * (r - 1) / (2 * k + 1)
    dec_upper = 4 * r / (2 * k + 1) ** 2
    corridor_lower = (2 * k - 1) * (2 * k + 3) * r / (2 * k + 1) ** 2
    corridor_upper = (2 + (2 * k - 1) * r) / (2 * k + 1)
    theta = (decrement - dec_lower) / (dec_upper - dec_lower)
    checks = {
        "upper-corridor-to-lower-decrement": sp.simplify((corridor_upper - rn) - (decrement - dec_lower)),
        "lower-corridor-to-upper-decrement": sp.simplify((rn - corridor_lower) - (dec_upper - decrement)),
        "decrement-width": sp.simplify(
            (dec_upper - dec_lower) - (2 * (2 * k + 1 - (2 * k - 1) * r) / (2 * k + 1) ** 2)
        ),
        "theta-lower-margin": sp.simplify(theta * (dec_upper - dec_lower) - (decrement - dec_lower)),
        "theta-upper-margin": sp.simplify((1 - theta) * (dec_upper - dec_lower) - (dec_upper - decrement)),
    }
    issues: list[DecrementIssue] = []
    for name, value in checks.items():
        if sp.simplify(value) != 0:
            issues.append(issue("symbolic", f"bad-{name}", str(value)))
    return issues


def validate_top_level(artifact: dict) -> list[DecrementIssue]:
    issues: list[DecrementIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "exact finite theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_raw_bridge",
        "source_raw_corridor_target",
        "source_raw_obstruction",
        "source_adaptive_obligations",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    for ref in artifact.get("enclosure_jsonl", []):
        issues.extend(validate_ref("<artifact>", ref))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in ("exact", "finite", "does not prove", "raw-corridor", "cone entry", "jwpf_06", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[DecrementIssue]:
    refs = artifact.get("enclosure_jsonl", [])
    if not isinstance(refs, list) or not refs:
        return [issue("recompute", "missing-enclosures", repr(refs))]
    try:
        recomputed = build_artifact([REPO_ROOT / ref for ref in refs])
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[DecrementIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_witnesses(rows_by_id: dict[str, dict]) -> list[DecrementIssue]:
    issues: list[DecrementIssue] = []
    expected = {
        "nlrdc_05_monotone_decrease_shortcut_blocked": exact_counterexample(Fraction(2), Fraction(3, 2)),
        "nlrdc_06_overfast_decrease_shortcut_blocked": exact_counterexample(Fraction(2), Fraction(1)),
    }
    for row_id, witness in expected.items():
        row = rows_by_id.get(row_id, {})
        if row.get("witness") != witness:
            issues.append(issue(row_id, "bad-witness", repr(row.get("witness"))))
        if not witness["raw_cone_holds_at_k_and_next"] or not witness["raw_decrease_holds"]:
            issues.append(issue(row_id, "bad-witness-basics", repr(witness)))
    first = expected["nlrdc_05_monotone_decrease_shortcut_blocked"]
    second = expected["nlrdc_06_overfast_decrease_shortcut_blocked"]
    if not str(first["lower_margin_D_minus_lower"]).startswith("-"):
        issues.append(issue("witness", "first-does-not-fail-lower-decrement", repr(first)))
    if not str(second["upper_margin_upper_minus_D"]).startswith("-"):
        issues.append(issue("witness", "second-does-not-fail-upper-decrement", repr(second)))
    return issues


def validate_rows(artifact: dict) -> tuple[list[DecrementIssue], int, int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[DecrementIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0, 0, 0
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))

    exact_available = 0
    exact_counterexamples = 0
    live_routes = 0
    ready_to_apply = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("matrix_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "claim", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if row.get("role") not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(row.get("role"))))
        readiness = row.get("readiness")
        if row.get("role") == "exact_reduction":
            exact_available += 1
            if readiness != "available_exact":
                issues.append(issue(row_id, "bad-exact-readiness", repr(readiness)))
        elif readiness != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-open-readiness", repr(readiness)))
        if readiness == "ready_to_apply":
            ready_to_apply += 1
        if row.get("role") == "exact_counterexample":
            exact_counterexamples += 1
        if row.get("role") == "live_route":
            live_routes += 1
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("exact", "finite", "not", "only", "open", "live", "hygiene")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    issues.extend(validate_witnesses(rows_by_id))
    return issues, len(rows), exact_available, exact_counterexamples, live_routes if ready_to_apply == 0 else -live_routes


def validate_summary(
    artifact: dict,
    row_count: int,
    exact_available: int,
    exact_counterexamples: int,
    live_routes: int,
) -> list[DecrementIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 9,
        "exact_available_rows": 2,
        "finite_pattern_rows": 2,
        "exact_counterexample_rows": 2,
        "live_routes": 1,
        "open_requirement_rows": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "lambdas": ["-25.0", "-50.0", "-100.0"],
        "coefficient_k_max": 200,
        "checked_ratio_max": 199,
        "adjacent_rows": 594,
        "raw_decrease_rows": 594,
        "lower_decrement_rows": 594,
        "upper_decrement_rows": 594,
        "decrement_corridor_rows": 594,
        "theta_unit_rows": 594,
        "theta_k_monotone_rows": 591,
        "theta_k_monotone_total_rows": 591,
        "theta_lambda_order_rows": 396,
        "theta_lambda_order_total_rows": 396,
    }
    issues: list[DecrementIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 9:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if exact_available != 2:
        issues.append(issue("summary", "bad-exact-row-count", str(exact_available)))
    if exact_counterexamples != 2:
        issues.append(issue("summary", "bad-counterexample-count", str(exact_counterexamples)))
    if live_routes != 1:
        issues.append(issue("summary", "bad-live-route-count", str(live_routes)))
    extrema = {
        "min_lower_decrement_margin": ("-25.0", 198, "2.241640098067743212E-6"),
        "min_upper_decrement_margin": ("-100.0", 198, "7.465690906027845625E-6"),
        "min_decrement_corridor_width": ("-100.0", 198, "1.112322974466921331E-5"),
        "min_theta": ("-25.0", 198, "1.644126617921189501E-1"),
        "max_theta": ("-100.0", 1, "9.087309076587821269E-1"),
        "min_theta_k_drop": ("-25.0", 197, "3.722452214840226847E-4"),
        "min_theta_lambda_gap": ("-25.0->-50.0", 198, "5.048218901591782104E-2"),
    }
    for key, (lam, k, sample) in extrema.items():
        row = summary.get(key, {})
        if (row.get("lam"), row.get("k"), row.get("sample")) != (lam, k, sample):
            issues.append(issue("summary", f"bad-{key}", repr(row)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "decrement recurrence",
        "594 adjacent",
        "theta coordinate",
        "k-monotone",
        "lambda-ordered",
        "exact raw-cone monotone counterexamples",
        "zeta-specific structure",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("available_exact", "ready_to_apply", "zeta-specific", "blocked", "finite shape", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[DecrementIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[DecrementIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "raw-corridor theorem is proved",
        "decrement recurrence is proved",
        "cone entry is proved",
        "jwpf_06 is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[DecrementIssue], dict]:
    artifact = load_json(target_path)
    issues: list[DecrementIssue] = []
    issues.extend(validate_symbolic_identities())
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, exact_available, exact_counterexamples, live_routes = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, exact_available, exact_counterexamples, live_routes))
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
            print(f"JWPF-NEG-LAMBDA-RAW-RATIO-DECREMENT {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda raw-ratio decrement-corridor scout: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('decrement_corridor_rows')} decrement-corridor rows, "
            f"{summary.get('theta_k_monotone_rows')} theta-k-monotone rows, "
            f"{summary.get('exact_counterexample_rows')} exact counterexamples, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
