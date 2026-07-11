#!/usr/bin/env python3
"""Validate the Jensen-window PF positive spectral moment obstruction."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OBSTRUCTION = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_positive_spectral_moment_obstruction.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_positive_spectral_moment_obstruction.md"

REQUIRED_SYMBOLIC_IDS = {
    "psm_01_degree2_hankel_obstruction",
    "psm_02_signed_hankel_signature_conflict",
    "psm_03_transform_escape_hatch",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Positive Spectral Moment Obstruction",
    "Status: ordinary positive spectral moment obstruction",
    "This is not a proof",
    "Delta_2 = det [[mu_0, mu_1], [mu_1, mu_2]] = -g_2",
    "work/rh_compute/results/jensen_window_pf_positive_spectral_moment_obstruction.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_positive_spectral_moment_obstruction.py",
    "validated Jensen-window PF positive spectral moment obstruction: 3 symbolic rows, 735 finite Delta2 obstruction rows, 0 issues",
    "psm_01_degree2_hankel_obstruction",
    "psm_02_signed_hankel_signature_conflict",
    "psm_03_transform_escape_hatch",
    "mu_m=int x^m dnu(x)",
    "nonordinary positive transform",
    "positive kernel identity",
    "outputs/jensen_window_pf_reciprocal_signed_j_fraction_scout.md",
    "validated Jensen-window PF reciprocal signed J-fraction scout: 2 symbolic rows, 3675 signed Hankel rows, 2940 signed-lambda rows, 0 issues",
)


@dataclass(frozen=True)
class SpectralMomentIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> SpectralMomentIssue:
    return SpectralMomentIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: str) -> list[SpectralMomentIssue]:
    if ref.startswith(("http://", "https://")):
        return []
    if (REPO_ROOT / ref).exists():
        return []
    return [issue(section, "missing-ref", ref)]


def validate_top_level(obstruction: dict) -> list[SpectralMomentIssue]:
    issues: list[SpectralMomentIssue] = []
    if obstruction.get("kind") != "jensen_window_pf_positive_spectral_moment_obstruction":
        issues.append(issue("<obstruction>", "bad-kind", repr(obstruction.get("kind"))))
    expected = {
        "parent_positive_readout_row": "pr_01_positive_spectral_transform",
        "parent_fit_row": "or_06_positive_spectral_measure_after_transform",
        "parent_route_row": "rp_09_signed_or_modified_continued_fraction",
        "parent_target": "jwpf_06_sign_regular_to_jensen_pf_conversion",
    }
    for key, value in expected.items():
        if obstruction.get(key) != value:
            issues.append(issue("<obstruction>", f"bad-{key}", f"{obstruction.get(key)!r} != {value!r}"))
    boundary = str(obstruction.get("proof_boundary", "")).lower()
    for required in ("ordinary positive spectral moment", "not a proof against nonordinary", "not lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<obstruction>", "weak-proof-boundary", required))
    objects = obstruction.get("objects", {})
    for required in ("H(t)", "E(t)", "mu_m", "ordinary_moment_hankel", "positive_measure_requirement"):
        if required not in objects or not str(objects.get(required, "")).strip():
            issues.append(issue("objects", "missing-object", required))
    return issues


def validate_symbolic_rows(obstruction: dict) -> tuple[list[SpectralMomentIssue], int]:
    rows = obstruction.get("symbolic_rows", [])
    issues: list[SpectralMomentIssue] = []
    if not isinstance(rows, list) or not rows:
        return [issue("symbolic_rows", "missing-symbolic-rows", repr(rows))], 0
    seen: set[str] = set()
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("symbolic_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        seen.add(row_id)
        for key in ("statement", "consequence", "proof_boundary"):
            if not str(row.get(key, "")).strip():
                issues.append(issue(row_id, f"missing-{key}", key))
        text = " ".join(str(row.get(key, "")) for key in ("statement", "consequence", "proof_boundary")).lower()
        if row_id == "psm_01_degree2_hankel_obstruction":
            for required in ("delta_2", "-g_2", "positive semidefinite"):
                if required not in text:
                    issues.append(issue(row_id, "missing-symbolic-text", required))
        if row_id == "psm_03_transform_escape_hatch":
            for required in ("not the ordinary moment", "non-power", "positive kernel"):
                if required not in text:
                    issues.append(issue(row_id, "missing-escape-text", required))
        for forbidden in ("therefore rh", "we have proved lambda <= 0", "jwpf_06 is proved"):
            if forbidden in text:
                issues.append(issue(row_id, "forbidden-overclaim", forbidden))
    for missing in sorted(REQUIRED_SYMBOLIC_IDS - seen):
        issues.append(issue(missing, "missing-symbolic-row", missing))
    return issues, len(rows)


def validate_finite_obstruction(obstruction: dict) -> tuple[list[SpectralMomentIssue], int]:
    finite = obstruction.get("finite_obstruction", {})
    issues: list[SpectralMomentIssue] = []
    for key in ("source_summary", "source_rows"):
        value = finite.get(key)
        if not isinstance(value, str):
            issues.append(issue("finite_obstruction", f"missing-{key}", repr(value)))
        else:
            issues.extend(validate_ref("finite_obstruction", value))
    if finite.get("determinant_order_r") != 2:
        issues.append(issue("finite_obstruction", "bad-determinant-order", repr(finite.get("determinant_order_r"))))
    if finite.get("expected_delta_sign") != -1:
        issues.append(issue("finite_obstruction", "bad-expected-delta-sign", repr(finite.get("expected_delta_sign"))))
    if finite.get("finite_delta2_obstruction_rows") != 735:
        issues.append(issue("finite_obstruction", "bad-row-count", repr(finite.get("finite_delta2_obstruction_rows"))))
    if finite.get("all_checked_delta2_rows_negative") is not True:
        issues.append(issue("finite_obstruction", "not-all-negative", repr(finite.get("all_checked_delta2_rows_negative"))))

    row_path_value = finite.get("source_rows")
    if not isinstance(row_path_value, str):
        return issues, 0
    row_path = REPO_ROOT / row_path_value
    if not row_path.exists():
        return issues, 0

    count = 0
    bad_rows: list[str] = []
    with row_path.open("r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("kind") != "jensen_window_pf_reciprocal_signed_hankel_row":
                continue
            if row.get("determinant_order_r") != 2:
                continue
            count += 1
            if row.get("expected_delta_sign") != -1 or row.get("ok") is not True or row.get("signed_delta_classification") != "positive":
                bad_rows.append(f"line {line_no}: {row!r}")
                if len(bad_rows) >= 5:
                    break
    if count != 735:
        issues.append(issue("finite_rows", "bad-delta2-row-count", str(count)))
    for bad in bad_rows:
        issues.append(issue("finite_rows", "bad-delta2-row", bad))
    return issues, count


def validate_summary(obstruction: dict, symbolic_rows: int, finite_rows: int) -> list[SpectralMomentIssue]:
    issues: list[SpectralMomentIssue] = []
    summary = obstruction.get("summary", {})
    expected = {
        "symbolic_rows": 3,
        "finite_delta2_obstruction_rows": 735,
        "ready_to_apply_rows": 0,
        "ordinary_positive_moment_route_rejected": True,
        "nonordinary_positive_transform_route_still_live": True,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if symbolic_rows != 3:
        issues.append(issue("summary", "bad-symbolic-count", str(symbolic_rows)))
    if finite_rows != 735:
        issues.append(issue("summary", "bad-finite-count", str(finite_rows)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("ordinary positive spectral moment", "delta_2=-g_2", "735", "nonordinary"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in obstruction.get("invariants", [])).lower()
    for required in ("ordinary positive moment", "does not reject nonordinary", "does not prove"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[SpectralMomentIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[SpectralMomentIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "ordinary moment route proves",
        "jwpf_06 is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(obstruction_path: Path, note_path: Path) -> tuple[list[SpectralMomentIssue], int, int]:
    obstruction = load_json(obstruction_path)
    issues: list[SpectralMomentIssue] = []
    issues.extend(validate_top_level(obstruction))
    symbolic_issues, symbolic_rows = validate_symbolic_rows(obstruction)
    finite_issues, finite_rows = validate_finite_obstruction(obstruction)
    issues.extend(symbolic_issues)
    issues.extend(finite_issues)
    issues.extend(validate_summary(obstruction, symbolic_rows, finite_rows))
    issues.extend(validate_note(note_path))
    return issues, symbolic_rows, finite_rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--obstruction", type=Path, default=DEFAULT_OBSTRUCTION)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, symbolic_rows, finite_rows = validate(args.obstruction, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "symbolic_rows": symbolic_rows,
                    "finite_delta2_obstruction_rows": finite_rows,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-SPECTRAL-MOMENT {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF positive spectral moment obstruction: "
            f"{symbolic_rows} symbolic rows, {finite_rows} finite Delta2 obstruction rows, "
            f"{len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
