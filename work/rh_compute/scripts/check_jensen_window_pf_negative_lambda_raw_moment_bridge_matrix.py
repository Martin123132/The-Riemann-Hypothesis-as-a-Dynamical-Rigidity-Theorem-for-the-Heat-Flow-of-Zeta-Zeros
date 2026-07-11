#!/usr/bin/env python3
"""Validate the negative-lambda raw-moment bridge matrix."""

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

from jensen_window_pf_negative_lambda_raw_moment_bridge_matrix import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_matrix,
)


REQUIRED_ROW_IDS = {
    "nlrmb_01_raw_ratio_reindexing",
    "nlrmb_02_exact_cone_translation",
    "nlrmb_03_adaptive_corridor_identity",
    "nlrmb_04_raw_cone_prefix",
    "nlrmb_05_raw_corridor_prefix",
    "nlrmb_06_fixed_buffer_thresholds",
    "nlrmb_07_live_theorem_burden",
    "nlrmb_08_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_reindexing",
    "exact_reduction",
    "finite_certificate",
    "finite_pattern",
    "rejected_route",
    "live_route",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Raw-Moment Bridge Matrix",
    "Status: exact finite theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_raw_moment_bridge_matrix`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.json",
    "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k200.jsonl",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.py",
    "validated Jensen-window PF negative-lambda raw-moment bridge matrix: 8 matrix rows, 0 issues, 597 raw-cone rows, 594 corridor rows, 76 half-width failures",
    "M_k = mu_{2k}",
    "A_k = M_k*k!/(2*k)!",
    "R_k = M_(k+1)*M_(k-1)/M_k^2",
    "x_k = ((2*k-1)/(2*k+1))*R_k",
    "s_k = ((2*k+1)-(2*k-1)*R_k)/2",
    "0 <= s_k <= 1 iff 1 <= R_k <= (2*k+1)/(2*k-1)",
    "((2*k-1)*(2*k+3)/(2*k+1)^2)*R_k <= R_(k+1) <= (2+(2*k-1)*R_k)/(2*k+1)",
    "corridor width = 2*(2*k+1-(2*k-1)*R_k)/(2*k+1)^2",
    "lambdas: -25.0, -50.0, -100.0",
    "coefficient range: A_0..A_200",
    "checked raw ratios: R_1..R_199",
    "raw lower wall rows: 597 / 597",
    "raw upper wall rows: 597 / 597",
    "corridor occupancy rows: 594 / 594",
    "corridor width rows: 594 / 594",
    "half-width failure rows: 76",
    "one-third failure rows: 418",
    "max raw ratio: 2.973078641741057661E+0 at lambda=-100.0, k=1",
    "min raw lower margin: 2.329147017095538537E-3 at lambda=-25.0, k=199",
    "min raw upper slack: 2.211618643155480210E-3 at lambda=-100.0, k=199",
    "min bridge lower margin: 7.465690906027845625E-6 at lambda=-100.0, k=198",
    "min scaled upper margin: 2.241640098067743212E-6 at lambda=-25.0, k=198",
    "min corridor width: 1.112322974466921331E-5 at lambda=-100.0, k=198",
    "s_k <= 1/2 iff R_k >= 2*k/(2*k-1)",
    "s_k <= 1/3 iff R_k >= (6*k+1)/(3*(2*k-1))",
    "outputs/jensen_window_pf_negative_lambda_adaptive_envelope_obligations.md",
    "outputs/jensen_window_pf_negative_lambda_adaptive_envelope_matrix.md",
    "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_k200_scout.md",
    "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
    "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
    "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md",
    "outputs/jensen_window_pf_negative_lambda_half_width_tail_target.md",
)


@dataclass(frozen=True)
class BridgeIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> BridgeIssue:
    return BridgeIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[BridgeIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_symbolic_identities() -> list[BridgeIssue]:
    k, r, r_next = sp.symbols("k r r_next")
    x = (2 * k - 1) * r / (2 * k + 1)
    x_next = (2 * k + 1) * r_next / (2 * k + 3)
    s = (2 * k + 1) * (1 - x) / 2
    s_next = (2 * k + 3) * (1 - x_next) / 2
    lower_factor = (2 * k - 1) * (2 * k + 3) / (2 * k + 1) ** 2
    upper_bound = (2 + (2 * k - 1) * r) / (2 * k + 1)
    checks = {
        "scaled-defect-raw": sp.simplify(s - ((2 * k + 1 - (2 * k - 1) * r) / 2)),
        "lower-wall": sp.simplify((x - (2 * k - 1) / (2 * k + 1)) - ((2 * k - 1) * (r - 1) / (2 * k + 1))),
        "upper-wall": sp.simplify((1 - x) - ((2 * k - 1) * ((2 * k + 1) / (2 * k - 1) - r) / (2 * k + 1))),
        "half-width": sp.simplify((sp.Rational(1, 2) - s) - (((2 * k - 1) * r - 2 * k) / 2)),
        "one-third": sp.simplify((sp.Rational(1, 3) - s) - ((3 * (2 * k - 1) * r - (6 * k + 1)) / 6)),
        "bridge-gap": sp.simplify(
            (x_next - x)
            - (((2 * k + 1) ** 2 * r_next - (2 * k - 1) * (2 * k + 3) * r) / ((2 * k + 1) * (2 * k + 3)))
        ),
        "bridge-factor": sp.simplify(lower_factor - (1 - sp.Rational(4, 1) / (2 * k + 1) ** 2)),
        "scaled-k-increase": sp.simplify((s_next - s) - ((2 + (2 * k - 1) * r - (2 * k + 1) * r_next) / 2)),
        "corridor-width": sp.simplify(
            (upper_bound - lower_factor * r) - (2 * (2 * k + 1 - (2 * k - 1) * r) / (2 * k + 1) ** 2)
        ),
    }
    issues: list[BridgeIssue] = []
    for name, value in checks.items():
        if sp.simplify(value) != 0:
            issues.append(issue("symbolic", f"bad-{name}", str(value)))
    return issues


def validate_top_level(artifact: dict) -> list[BridgeIssue]:
    issues: list[BridgeIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_raw_moment_bridge_matrix":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "exact finite theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_adaptive_obligations",
        "source_adaptive_matrix",
        "source_cone_prefix",
        "source_boundary_threshold",
        "source_monotone_contraction_target",
        "source_bounded_log_curvature_target",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    for ref in artifact.get("enclosure_jsonl", []):
        issues.extend(validate_ref("<artifact>", ref))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in ("exact", "finite", "does not prove", "raw moment", "corridor", "cone entry", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[BridgeIssue]:
    refs = artifact.get("enclosure_jsonl", [])
    if not isinstance(refs, list) or not refs:
        return [issue("recompute", "missing-enclosures", repr(refs))]
    try:
        recomputed = build_matrix([REPO_ROOT / ref for ref in refs])
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[BridgeIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[BridgeIssue], int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[BridgeIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0, 0
    by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(by_id)):
        issues.append(issue(missing, "missing-row", missing))
    exact_available = 0
    live_routes = 0
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
        if row.get("role") in {"exact_reindexing", "exact_reduction"}:
            if readiness != "available_exact":
                issues.append(issue(row_id, "bad-exact-readiness", repr(readiness)))
            exact_available += 1
        elif readiness != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-finite-readiness", repr(readiness)))
        if row.get("role") == "live_route":
            live_routes += 1
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("exact", "finite", "not", "only", "theorem-search")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    return issues, len(rows), exact_available, live_routes


def validate_summary(artifact: dict, row_count: int, exact_available: int, live_routes: int) -> list[BridgeIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 8,
        "ready_to_apply_rows": 0,
        "coefficient_k_max": 200,
        "checked_ratio_max": 199,
        "lambdas": ["-25.0", "-50.0", "-100.0"],
        "raw_ratio_rows": 597,
        "raw_lower_positive_rows": 597,
        "raw_upper_positive_rows": 597,
        "raw_cone_positive_rows": 597,
        "half_width_positive_rows": 521,
        "half_width_failure_rows": 76,
        "one_third_positive_rows": 179,
        "one_third_failure_rows": 418,
        "corridor_rows": 594,
        "bridge_lower_positive_rows": 594,
        "scaled_upper_positive_rows": 594,
        "corridor_occupancy_positive_rows": 594,
        "corridor_width_positive_rows": 594,
        "target_closing": False,
    }
    issues: list[BridgeIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 8:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if exact_available != 3:
        issues.append(issue("summary", "bad-exact-row-count", str(exact_available)))
    if live_routes != 1:
        issues.append(issue("summary", "bad-live-route-count", str(live_routes)))
    extrema = {
        "max_raw_ratio": ("-100.0", 1, "2.973078641741057661E+0"),
        "min_raw_lower_margin": ("-25.0", 199, "2.329147017095538537E-3"),
        "min_raw_upper_slack": ("-100.0", 199, "2.211618643155480210E-3"),
        "min_raw_cone_margin": ("-100.0", 199, "2.211618643155480210E-3"),
        "min_bridge_lower_margin": ("-100.0", 198, "7.465690906027845625E-6"),
        "min_scaled_upper_margin": ("-25.0", 198, "2.241640098067743212E-6"),
        "min_corridor_occupancy_margin": ("-25.0", 198, "2.241640098067743212E-6"),
        "min_corridor_width": ("-100.0", 198, "1.112322974466921331E-5"),
    }
    for key, (lam, k, sample) in extrema.items():
        row = summary.get(key, {})
        if (row.get("lam"), row.get("k"), row.get("sample")) != (lam, k, sample):
            issues.append(issue("summary", f"bad-{key}", repr(row)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "raw-moment ratio corridor",
        "k200",
        "597/597 lower",
        "597/597 upper",
        "594/594 adjacent",
        "nonempty precisely under the upper raw wall",
        "all-k burden",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("available_exact", "not the bottleneck", "upper raw wall", "corridor occupancy", "half-width", "one-third", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[BridgeIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[BridgeIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "adaptive scaled-defect target is proved",
        "cone entry is proved",
        "jwpf_06 is proved",
        "all-k raw moment-growth theorem is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[BridgeIssue], dict]:
    artifact = load_json(target_path)
    issues: list[BridgeIssue] = []
    issues.extend(validate_symbolic_identities())
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, exact_available, live_routes = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, exact_available, live_routes))
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
            print(f"JWPF-NEG-LAMBDA-RAW-MOMENT-BRIDGE {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda raw-moment bridge matrix: "
            f"{summary.get('matrix_rows')} matrix rows, {len(issues)} issues, "
            f"{summary.get('raw_cone_positive_rows')} raw-cone rows, "
            f"{summary.get('corridor_occupancy_positive_rows')} corridor rows, "
            f"{summary.get('half_width_failure_rows')} half-width failures"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
