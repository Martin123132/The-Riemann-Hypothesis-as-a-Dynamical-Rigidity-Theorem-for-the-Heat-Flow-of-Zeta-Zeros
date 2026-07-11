#!/usr/bin/env python3
"""Validate the Jensen-window PF heat-flow cone-entry asymptotic target."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_TARGET = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md"

REQUIRED_ROW_IDS = {
    "hfcet_01_exact_cone_entry_statement",
    "hfcet_02_fixed_k_asymptotic_formal_reduction",
    "hfcet_03_gaussian_leading_degeneracy",
    "hfcet_04_phi_taylor_sign_route",
    "hfcet_05_uniform_k_or_tail_collar_gap",
    "hfcet_06_finite_grid_support",
    "hfcet_07_forbidden_endpoint_shortcuts",
    "hfcet_08_conditional_application",
}

ALLOWED_ROLES = {
    "entry_theorem",
    "open_statement",
    "formal_reduction",
    "diagnostic_obstruction",
    "live_route",
    "hard_gap",
    "finite_evidence",
    "circular_route",
    "conditional_application",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Heat-Flow Cone-Entry Asymptotic Target",
    "Status: open forward-flow theorem target",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_heat_flow_cone_entry_asymptotic_target`",
    "work/rh_compute/results/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_cone_entry_asymptotic_target.py",
    "validated Jensen-window PF heat-flow cone-entry asymptotic target: 8 rows, 0 issues, 1 ready-to-apply rows, 0 live routes",
    "(2*k-1)/(2*k+1) <= x_k <= 1",
    "x_{k+1} >= x_k",
    "outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md",
    "validated Jensen-window PF heat-flow ratio cone invariance lemma: 6 exact rows, 315 lower rows, 315 upper rows, 310 monotone rows, 0 issues",
    "A_k(-T)=C*4^(-k)*T^(-k-1/2)",
    "log x_k",
    "log(x_{k+1}/x_k)",
    "2*b-a^2 < 0",
    "2*(a^3-3*a*b+3*c) > 0",
    "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md",
    "validated Jensen-window PF Phi Taylor cone-entry sign scout: 4 coefficient balls, 2 certified signs, 0 ready-to-apply rows, 0 issues",
    "leading Gaussian term gives `x_k -> 1`",
    "outputs/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.md",
    "explicit control preventing a first exit from escaping to k=infinity",
    "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
    "validated Jensen-window PF heat-flow boundary threshold lemma: 5 exact rows, 315 strong-threshold rows, 315 heat-threshold rows, 0 issues",
    "outputs/jensen_window_pf_monotone_contraction_stress.md",
    "validated Jensen-window PF monotone contraction stress: 2875 rows, 2875 positive rows, 0 issues",
    "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.md",
    "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md",
    "validated Jensen-window PF negative-lambda cone-entry prefix scout: 69 coefficient rows, 63 lower-wall rows, 63 upper-wall rows, 60 monotone-gap rows, 0 issues",
    "outputs/jensen_window_pf_negative_lambda_finite_collar_contract.md",
    "validated Jensen-window PF negative-lambda finite-collar contract: active depth K=19, 57 active lower rows, 57 active upper rows, 57 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues",
    "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
    "validated Jensen-window PF negative-lambda defect-tail theorem target: 8 rows, 0 issues, 0 ready-to-apply rows, 2 live routes",
    "0 <= d_k <= 2/(2*k+1)",
    "d_(k+1)<=d_k",
    "lambda=-25,-50,-100",
    "`K=19`, with first collar `x_20` and second collar `x_21`",
    "finite-collar flow theorem",
)


@dataclass(frozen=True)
class ConeEntryIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> ConeEntryIssue:
    return ConeEntryIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[ConeEntryIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(target: dict) -> list[ConeEntryIssue]:
    issues: list[ConeEntryIssue] = []
    if target.get("kind") != "jensen_window_pf_heat_flow_cone_entry_asymptotic_target":
        issues.append(issue("<target>", "bad-kind", repr(target.get("kind"))))
    if target.get("status") != "open_theorem_target":
        issues.append(issue("<target>", "bad-status", repr(target.get("status"))))
    if target.get("target_id") != "target_zeta_ratio_cone_entry":
        issues.append(issue("<target>", "bad-target-id", repr(target.get("target_id"))))
    boundary = str(target.get("proof_boundary", "")).lower()
    for required in ("open theorem target", "proves full ratio-cone entry", "collared finite", "jwpf_06", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<target>", "weak-proof-boundary", required))
    cone = target.get("cone", {})
    for key in ("ratio", "lower_wall", "upper_wall", "monotone_wall", "entry_target"):
        if key not in cone:
            issues.append(issue("cone", "missing-cone-field", key))
    formal = target.get("formal_asymptotic", {})
    for key in ("weight_taylor", "moment_expansion", "log_x_expansion", "monotone_log_gap_expansion", "entry_sign_requirements"):
        if key not in formal:
            issues.append(issue("formal_asymptotic", "missing-field", key))
    return issues


def validate_target_rows(target: dict) -> tuple[list[ConeEntryIssue], int, int, int]:
    rows = target.get("target_rows", [])
    issues: list[ConeEntryIssue] = []
    if not isinstance(rows, list):
        return [issue("target_rows", "bad-rows", repr(type(rows)))], 0, 0, 0
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        issues.append(issue("target_rows", "bad-row-ids", repr(sorted(ids))))
    ready = 0
    live = 0
    hard = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("target_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "source_artifacts", "claim_if_proved", "gap", "acceptance_test", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if row.get("role") not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(row.get("role"))))
        if row.get("readiness") == "ready_to_apply":
            ready += 1
            if row_id != "hfcet_01_exact_cone_entry_statement":
                issues.append(issue(row_id, "ready-to-apply-overclaim", row_id))
        if row.get("role") == "live_route":
            live += 1
        if row.get("role") == "hard_gap":
            hard += 1
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("open", "only", "gap", "circular", "finite", "formal", "theorem")):
            issues.append(issue(row_id, "weak-row-boundary", boundary))
    return issues, len(rows), ready, live + hard


def validate_summary(target: dict) -> list[ConeEntryIssue]:
    summary = target.get("summary", {})
    issues: list[ConeEntryIssue] = []
    expected = {
        "target_rows": 8,
        "ready_to_apply_rows": 1,
        "live_routes": 0,
        "formal_reduction_rows": 1,
        "finite_evidence_rows": 1,
        "hard_gap_rows": 1,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "full infinite ratio-cone entry is certified",
        "fixed-k",
        "local phi taylor signs",
        "repaired finite prefix",
        "analytic adjacent tail",
        "forward infinite-dimensional or collared finite-flow legitimacy",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in target.get("invariants", [])).lower()
    for required in (
        "cone-entry theorem row is ready_to_apply",
        "open_target",
        "fixed-k",
        "phi taylor sign scout",
        "negative-lambda prefix scout",
        "finite-collar contract",
        "finite grid",
        "lambda <= 0",
    ):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[ConeEntryIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ConeEntryIssue] = []
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


def validate(target_path: Path, note_path: Path) -> tuple[list[ConeEntryIssue], int, int, int]:
    target = load_json(target_path)
    issues: list[ConeEntryIssue] = []
    issues.extend(validate_top_level(target))
    row_issues, rows, ready, live_plus_hard = validate_target_rows(target)
    issues.extend(row_issues)
    issues.extend(validate_summary(target))
    issues.extend(validate_note(note_path))
    live_routes = int(target.get("summary", {}).get("live_routes", 0))
    return issues, rows, ready, live_routes


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_TARGET)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, rows, ready, live = validate(args.target, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "rows": rows,
                    "ready_to_apply_rows": ready,
                    "live_routes": live,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-HEAT-FLOW-CONE-ENTRY {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF heat-flow cone-entry asymptotic target: "
            f"{rows} rows, {len(issues)} issues, {ready} ready-to-apply rows, {live} live routes"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
