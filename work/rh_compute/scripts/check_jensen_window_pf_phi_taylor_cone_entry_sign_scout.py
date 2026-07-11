#!/usr/bin/env python3
"""Validate the Phi Taylor sign scout for the heat-flow cone-entry route."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_TARGET = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_phi_taylor_cone_entry_sign_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md"

REQUIRED_COEFFICIENT_SIGNS = {
    "phi_taylor_c0": "positive",
    "phi_taylor_c2": "negative",
    "phi_taylor_c4": "positive",
    "phi_taylor_c6": "negative",
}

REQUIRED_COMBINATION_SIGNS = {
    "ptces_01_upper_wall_sign": "negative",
    "ptces_02_monotone_wall_sign": "positive",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Phi Taylor Cone-Entry Sign Scout",
    "Status: finite Taylor sign certificate",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_phi_taylor_cone_entry_sign_scout`",
    "work/rh_compute/results/jensen_window_pf_phi_taylor_cone_entry_sign_scout.json",
    "python work/rh_compute/scripts/jensen_window_pf_phi_taylor_cone_entry_sign_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_phi_taylor_cone_entry_sign_scout.py",
    "validated Jensen-window PF Phi Taylor cone-entry sign scout: 4 coefficient balls, 2 certified signs, 0 ready-to-apply rows, 0 issues",
    "2*b-a^2 < 0",
    "2*(a^3-3*a*b+3*c) > 0",
    "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
    "uniform-in-`k`",
    "collared finite",
)


@dataclass(frozen=True)
class TaylorScoutIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> TaylorScoutIssue:
    return TaylorScoutIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def arb_positive(value: flint.arb) -> bool:
    return bool(value > 0 and not value.contains(0))


def arb_negative(value: flint.arb) -> bool:
    return bool(value < 0 and not value.contains(0))


def validate_ref(section: str, ref: object) -> list[TaylorScoutIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(target: dict) -> list[TaylorScoutIssue]:
    issues: list[TaylorScoutIssue] = []
    if target.get("kind") != "jensen_window_pf_phi_taylor_cone_entry_sign_scout":
        issues.append(issue("<target>", "bad-kind", repr(target.get("kind"))))
    if target.get("status") != "finite_taylor_sign_certificate":
        issues.append(issue("<target>", "bad-status", repr(target.get("status"))))
    for key in ("source_target", "generator", "checker"):
        issues.extend(validate_ref("<target>", target.get(key)))
    boundary = str(target.get("proof_boundary", "")).lower()
    for required in (
        "local taylor",
        "does not prove zeta cone entry",
        "uniformly in k",
        "collared finite",
        "jwpf_06",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<target>", "weak-proof-boundary", required))
    params = target.get("parameters", {})
    if params.get("precision_bits", 0) < 128:
        issues.append(issue("parameters", "low-precision", repr(params.get("precision_bits"))))
    if params.get("tail_cutoff_n", 0) < 8:
        issues.append(issue("parameters", "weak-tail-cutoff", repr(params.get("tail_cutoff_n"))))
    return issues


def validate_coefficients(target: dict) -> tuple[list[TaylorScoutIssue], int]:
    rows = target.get("coefficient_enclosures", [])
    issues: list[TaylorScoutIssue] = []
    if not isinstance(rows, list):
        return [issue("coefficient_enclosures", "bad-rows", repr(type(rows)))], 0
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != set(REQUIRED_COEFFICIENT_SIGNS):
        issues.append(issue("coefficient_enclosures", "bad-row-ids", repr(sorted(ids))))
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("coefficient_enclosures", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "coefficient", "taylor_degree", "finite_sum", "tail_radius", "enclosure", "sign"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        expected_sign = REQUIRED_COEFFICIENT_SIGNS.get(row_id)
        if row.get("sign") != expected_sign:
            issues.append(issue(row_id, "bad-recorded-sign", f"{row.get('sign')!r} != {expected_sign!r}"))
        try:
            ball = flint.arb(str(row.get("enclosure")))
        except Exception as exc:  # pragma: no cover - defensive parser guard
            issues.append(issue(row_id, "bad-arb-ball", f"{type(exc).__name__}: {exc}"))
            continue
        if expected_sign == "positive" and not arb_positive(ball):
            issues.append(issue(row_id, "uncertified-positive", str(ball)))
        if expected_sign == "negative" and not arb_negative(ball):
            issues.append(issue(row_id, "uncertified-negative", str(ball)))
    return issues, len(rows)


def validate_combinations(target: dict) -> tuple[list[TaylorScoutIssue], int, int]:
    rows = target.get("sign_combinations", [])
    issues: list[TaylorScoutIssue] = []
    if not isinstance(rows, list):
        return [issue("sign_combinations", "bad-rows", repr(type(rows)))], 0, 0
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != set(REQUIRED_COMBINATION_SIGNS):
        issues.append(issue("sign_combinations", "bad-row-ids", repr(sorted(ids))))
    ready = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("sign_combinations", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "formula", "enclosure", "certified_sign", "consequence", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        expected_sign = REQUIRED_COMBINATION_SIGNS.get(row_id)
        if row.get("certified_sign") != expected_sign:
            issues.append(issue(row_id, "bad-recorded-sign", f"{row.get('certified_sign')!r} != {expected_sign!r}"))
        try:
            ball = flint.arb(str(row.get("enclosure")))
        except Exception as exc:  # pragma: no cover - defensive parser guard
            issues.append(issue(row_id, "bad-arb-ball", f"{type(exc).__name__}: {exc}"))
            continue
        if expected_sign == "positive" and not arb_positive(ball):
            issues.append(issue(row_id, "uncertified-positive", str(ball)))
        if expected_sign == "negative" and not arb_negative(ball):
            issues.append(issue(row_id, "uncertified-negative", str(ball)))
        boundary = str(row.get("proof_boundary", "")).lower()
        for required in ("local taylor", "does not give", "uniform-in-k"):
            if required not in boundary:
                issues.append(issue(row_id, "weak-row-boundary", required))
        if row.get("readiness") == "ready_to_apply":
            ready += 1
            issues.append(issue(row_id, "ready-to-apply-overclaim", row_id))
    return issues, len(rows), ready


def validate_summary(target: dict) -> list[TaylorScoutIssue]:
    summary = target.get("summary", {})
    issues: list[TaylorScoutIssue] = []
    expected = {
        "coefficient_balls": 4,
        "certified_signs": 2,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("local phi taylor", "desired certified signs", "uniform-in-k", "collared finite"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in target.get("invariants", [])).lower()
    for required in ("u=0", "fixed-k", "cone-entry theorem target remains open", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[TaylorScoutIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[TaylorScoutIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "cone entry is proved",
        "monotone contractions are proved for all zeta windows",
        "analytic monotone-contraction theorem is proved",
        "closed differential inequality is proved",
        "jwpf_06 is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[TaylorScoutIssue], int, int, int]:
    target = load_json(target_path)
    issues: list[TaylorScoutIssue] = []
    issues.extend(validate_top_level(target))
    coeff_issues, coeffs = validate_coefficients(target)
    combo_issues, combos, ready = validate_combinations(target)
    issues.extend(coeff_issues)
    issues.extend(combo_issues)
    issues.extend(validate_summary(target))
    issues.extend(validate_note(note_path))
    return issues, coeffs, combos, ready


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_TARGET)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, coeffs, combos, ready = validate(args.target, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "coefficient_balls": coeffs,
                    "certified_signs": combos,
                    "ready_to_apply_rows": ready,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-PHI-TAYLOR-SCOUT {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF Phi Taylor cone-entry sign scout: "
            f"{coeffs} coefficient balls, {combos} certified signs, "
            f"{ready} ready-to-apply rows, {len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
