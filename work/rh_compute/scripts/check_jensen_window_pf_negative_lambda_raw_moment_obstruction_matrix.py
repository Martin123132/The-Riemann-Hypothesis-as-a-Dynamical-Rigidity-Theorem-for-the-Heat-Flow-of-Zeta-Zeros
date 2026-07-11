#!/usr/bin/env python3
"""Validate exact raw-moment obstruction gates for the adaptive negative-lambda route."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
    corridor_lower,
    corridor_upper,
    corridor_width,
    exact_values,
    raw_ratio,
    raw_upper_wall,
    scaled_defect_from_raw,
    two_atom_moments,
    x_from_raw,
)


REQUIRED_ROW_IDS = {
    "nlrmo_01_two_atom_family",
    "nlrmo_02_upper_raw_wall_failure",
    "nlrmo_03_scaled_upper_corridor_failure",
    "nlrmo_04_monotone_lower_corridor_failure",
    "nlrmo_05_corridor_width_identity",
    "nlrmo_06_zeta_prefix_contrast",
    "nlrmo_07_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_countermodel_family",
    "exact_counterexample",
    "exact_identity",
    "finite_contrast",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Raw-Moment Obstruction Matrix",
    "Status: exact countermodel gate",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.py",
    "validated Jensen-window PF negative-lambda raw-moment obstruction matrix: 7 matrix rows, 0 issues, 3 exact counterexamples, 0 ready-to-apply rows",
    "M_k = 1 + w*a^k",
    "R_k = M_(k+1)*M_(k-1)/M_k^2",
    "measure: delta_1 + (1/16)*delta_16",
    "R_1 = 289/64",
    "x_1 = 289/192",
    "s_1 = -97/128",
    "measure: delta_1 + (1)*delta_2",
    "R_2 = 27/25",
    "corridor upper = 28/27",
    "upper margin = -29/675",
    "s_2-s_1 = -29/450",
    "measure: delta_1 + (1/3)*delta_9",
    "R_1 = 7/3",
    "R_2 = 61/49",
    "corridor lower = 35/27",
    "lower margin = -68/1323",
    "x_2-x_1 = -68/2205",
    "corridor width = 2*(2*k+1-(2*k-1)*R_k)/(2*k+1)^2",
    "outputs/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.md",
    "validated 597/597 raw-cone rows and 594/594 corridor rows",
)


@dataclass(frozen=True)
class ObstructionIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> ObstructionIssue:
    return ObstructionIssue(section=section, issue=name, detail=detail)


def parse_fraction(text: object) -> Fraction:
    if not isinstance(text, str):
        raise ValueError(f"not a string fraction: {text!r}")
    return Fraction(text)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[ObstructionIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_exact_helpers() -> list[ObstructionIssue]:
    issues: list[ObstructionIssue] = []
    cases = [
        (Fraction(1, 16), Fraction(16), "289/64", "514/289"),
        (Fraction(1), Fraction(2), "10/9", "27/25"),
        (Fraction(1, 3), Fraction(9), "7/3", "61/49"),
    ]
    for weight, support, expected_r1, expected_r2 in cases:
        moments = two_atom_moments(weight, support, 3)
        if raw_ratio(moments, 1) != Fraction(expected_r1):
            issues.append(issue("helpers", "bad-r1", f"{weight=} {support=}"))
        if raw_ratio(moments, 2) != Fraction(expected_r2):
            issues.append(issue("helpers", "bad-r2", f"{weight=} {support=}"))
    k = 1
    raw = Fraction(10, 9)
    if corridor_lower(k, raw) != Fraction(50, 81):
        issues.append(issue("helpers", "bad-corridor-lower", str(corridor_lower(k, raw))))
    if corridor_upper(k, raw) != Fraction(28, 27):
        issues.append(issue("helpers", "bad-corridor-upper", str(corridor_upper(k, raw))))
    if corridor_width(k, raw) != Fraction(34, 81):
        issues.append(issue("helpers", "bad-corridor-width", str(corridor_width(k, raw))))
    if raw_upper_wall(1) != Fraction(3):
        issues.append(issue("helpers", "bad-upper-wall", str(raw_upper_wall(1))))
    if x_from_raw(1, Fraction(289, 64)) != Fraction(289, 192):
        issues.append(issue("helpers", "bad-x-from-raw", str(x_from_raw(1, Fraction(289, 64)))))
    if scaled_defect_from_raw(1, Fraction(289, 64)) != Fraction(-97, 128):
        issues.append(issue("helpers", "bad-scaled-defect", str(scaled_defect_from_raw(1, Fraction(289, 64)))))
    recomputed = exact_values(Fraction(1, 3), Fraction(9))
    if recomputed.get("lower_margin_R_next_minus_lower") != "-68/1323":
        issues.append(issue("helpers", "bad-exact-values", repr(recomputed)))
    return issues


def validate_top_level(artifact: dict) -> list[ObstructionIssue]:
    issues: list[ObstructionIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "exact countermodel gate":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in ("source_raw_moment_bridge", "source_adaptive_target", "source_cone_prefix", "generator", "checker"):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in ("countermodel", "blocks", "not a claim about actual zeta", "not an all-k theorem", "cone entry", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[ObstructionIssue]:
    recomputed = build_artifact()
    issues: list[ObstructionIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_witnesses(by_id: dict[str, dict]) -> list[ObstructionIssue]:
    issues: list[ObstructionIssue] = []
    upper = by_id["nlrmo_02_upper_raw_wall_failure"]["witness"]
    expected_upper = {
        "measure": "delta_1 + (1/16)*delta_16",
        "R_k": "289/64",
        "raw_upper_slack": "-97/64",
        "x_k": "289/192",
        "s_k": "-97/128",
        "corridor_width": "-97/288",
    }
    scaled = by_id["nlrmo_03_scaled_upper_corridor_failure"]["witness"]
    expected_scaled = {
        "measure": "delta_1 + (1)*delta_2",
        "R_k": "10/9",
        "R_next": "27/25",
        "raw_upper_slack": "17/9",
        "raw_next_upper_slack": "44/75",
        "corridor_upper": "28/27",
        "upper_margin_upper_minus_R_next": "-29/675",
        "s_next_minus_s_k": "-29/450",
    }
    lower = by_id["nlrmo_04_monotone_lower_corridor_failure"]["witness"]
    expected_lower = {
        "measure": "delta_1 + (1/3)*delta_9",
        "R_k": "7/3",
        "R_next": "61/49",
        "raw_upper_slack": "2/3",
        "raw_next_upper_slack": "62/147",
        "corridor_lower": "35/27",
        "lower_margin_R_next_minus_lower": "-68/1323",
        "x_next_minus_x_k": "-68/2205",
    }
    for section, witness, expected in (
        ("upper-witness", upper, expected_upper),
        ("scaled-upper-witness", scaled, expected_scaled),
        ("monotone-lower-witness", lower, expected_lower),
    ):
        for key, value in expected.items():
            if witness.get(key) != value:
                issues.append(issue(section, f"bad-{key}", f"{witness.get(key)!r} != {value!r}"))
    if parse_fraction(upper["R_k"]) <= parse_fraction(upper["raw_upper_wall"]):
        issues.append(issue("upper-witness", "does-not-violate-upper-wall", repr(upper)))
    if parse_fraction(scaled["upper_margin_upper_minus_R_next"]) >= 0:
        issues.append(issue("scaled-upper-witness", "does-not-violate-upper-corridor", repr(scaled)))
    if parse_fraction(lower["lower_margin_R_next_minus_lower"]) >= 0:
        issues.append(issue("monotone-lower-witness", "does-not-violate-lower-corridor", repr(lower)))
    if parse_fraction(scaled["raw_upper_slack"]) <= 0 or parse_fraction(scaled["raw_next_upper_slack"]) <= 0:
        issues.append(issue("scaled-upper-witness", "raw-upper-wall-not-held", repr(scaled)))
    if parse_fraction(lower["raw_upper_slack"]) <= 0 or parse_fraction(lower["raw_next_upper_slack"]) <= 0:
        issues.append(issue("monotone-lower-witness", "raw-upper-wall-not-held", repr(lower)))
    return issues


def validate_rows(artifact: dict) -> tuple[list[ObstructionIssue], int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[ObstructionIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0
    by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(by_id)):
        issues.append(issue(missing, "missing-row", missing))
    counterexamples = 0
    ready_to_apply = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("matrix_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "status", "claim", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if row.get("role") not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(row.get("role"))))
        if row.get("status") not in {"guard_validated", "available_exact", "finite_validated"}:
            issues.append(issue(row_id, "bad-status", repr(row.get("status"))))
        if row.get("role") == "exact_counterexample":
            counterexamples += 1
            if "witness" not in row:
                issues.append(issue(row_id, "missing-witness", row_id))
        if row.get("status") == "ready_to_apply":
            ready_to_apply += 1
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("not", "blocks", "only", "finite", "gate")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    if not issues:
        issues.extend(validate_witnesses(by_id))
    return issues, len(rows), counterexamples


def validate_summary(artifact: dict, row_count: int, counterexamples: int) -> list[ObstructionIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 7,
        "exact_counterexample_rows": 3,
        "exact_identity_rows": 2,
        "finite_contrast_rows": 1,
        "acceptance_gate_rows": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    issues: list[ObstructionIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 7:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if counterexamples != 3:
        issues.append(issue("summary", "bad-counterexample-count", str(counterexamples)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("generic positive stieltjes", "upper raw wall", "scaled-upper", "monotone-bridge", "zeta-specific", "all-k"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("positive two-atom", "not evidence against", "upper raw wall", "corridor occupancy", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[ObstructionIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ObstructionIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "adaptive scaled-defect target is proved",
        "raw upper wall is proved",
        "corridor theorem is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[ObstructionIssue], dict]:
    artifact = load_json(target_path)
    issues: list[ObstructionIssue] = []
    issues.extend(validate_exact_helpers())
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, counterexamples = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, counterexamples))
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
            print(f"JWPF-NEG-LAMBDA-RAW-MOMENT-OBSTRUCTION {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda raw-moment obstruction matrix: "
            f"{summary.get('matrix_rows')} matrix rows, {len(issues)} issues, "
            f"{summary.get('exact_counterexample_rows')} exact counterexamples, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
