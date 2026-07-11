#!/usr/bin/env python3
"""Validate the Jensen-window PF nonpower-functional low-degree scout."""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import asdict, dataclass
import json
from math import comb
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SCOUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_nonpower_functional_low_degree_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_nonpower_functional_low_degree_scout.md"
MAX_FORMULA_M = 6
MAX_COMPOSITION_M = 8

REQUIRED_IDS = {
    "npf_01_exact_reciprocal_polynomials",
    "npf_02_signed_composition_expansion",
    "npf_03_positive_g_cone_obstruction",
    "npf_04_degree2_log_concavity_entry",
    "npf_05_degree3_cancellation_entry",
    "npf_06_tautological_basis_trap",
    "npf_07_positive_cone_contract",
}

ALLOWED_VERDICTS = {
    "exact_identity",
    "rejected_as_generic_positive_cone",
    "first_required_inequality",
    "higher_cancellation_required",
    "rejected_as_tautological",
    "live_if_cone_and_functional_constructed",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Nonpower Functional Low-Degree Scout",
    "Status: exact low-degree nonpower-functional scout",
    "This is not a proof",
    "npt_04_nonpower_positive_functional",
    "work/rh_compute/results/jensen_window_pf_nonpower_functional_low_degree_scout.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_nonpower_functional_low_degree_scout.py",
    "validated Jensen-window PF nonpower functional low-degree scout: 7 scout rows, 0 issues, 0 ready-to-apply rows, 1 live contract rows",
    "mu_2 = g1^2 - g2",
    "mu_3 = g1^3 - 2*g1*g2 + g3",
    "mu_6 = g1^6 - 5*g1^4*g2",
    "m=8: 64 positive, 64 negative",
    "npf_01_exact_reciprocal_polynomials",
    "npf_07_positive_cone_contract",
    "outputs/jensen_window_pf_nonordinary_positive_transform_ansatz_matrix.md",
    "outputs/jensen_window_pf_positive_spectral_moment_obstruction.md",
    "outputs/jensen_window_pf_nonpower_functional_cone_candidate_matrix.md",
    "validated Jensen-window PF nonpower functional cone candidate matrix: 8 cone rows, 0 issues, 0 ready-to-apply rows, 2 live cone rows",
    "Kill Gates",
)


@dataclass(frozen=True)
class NonpowerFunctionalIssue:
    section: str
    issue: str
    detail: str


Poly = dict[tuple[int, ...], int]


def issue(section: str, name: str, detail: str) -> NonpowerFunctionalIssue:
    return NonpowerFunctionalIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def add_poly(left: Poly, right: Poly) -> Poly:
    out: defaultdict[tuple[int, ...], int] = defaultdict(int)
    for exp, coeff in left.items():
        out[exp] += coeff
    for exp, coeff in right.items():
        out[exp] += coeff
    return {exp: coeff for exp, coeff in out.items() if coeff}


def mul_g(index: int, poly: Poly, max_m: int) -> Poly:
    out: defaultdict[tuple[int, ...], int] = defaultdict(int)
    for exp, coeff in poly.items():
        new_exp = list(exp)
        new_exp[index - 1] += 1
        out[tuple(new_exp)] += coeff
    return {exp: coeff for exp, coeff in out.items() if coeff}


def reciprocal_polynomials(max_m: int) -> list[Poly]:
    zero = (0,) * max_m
    polys: list[Poly] = [{zero: 1}]
    for m in range(1, max_m + 1):
        poly: Poly = {}
        for j in range(1, m + 1):
            sign = 1 if j % 2 == 1 else -1
            term = {exp: sign * coeff for exp, coeff in mul_g(j, polys[m - j], max_m).items()}
            poly = add_poly(poly, term)
        polys.append(poly)
    return polys


def monomial_string(exp: tuple[int, ...]) -> str:
    parts: list[str] = []
    for index, power in enumerate(exp, start=1):
        if power == 1:
            parts.append(f"g{index}")
        elif power > 1:
            parts.append(f"g{index}^{power}")
    return "*".join(parts) or "1"


def poly_string(poly: Poly) -> str:
    items = sorted(poly.items(), key=lambda item: (-sum(item[0]), item[0]))
    chunks: list[str] = []
    for exp, coeff in items:
        mono = monomial_string(exp)
        abs_coeff = abs(coeff)
        if mono == "1":
            body = str(abs_coeff)
        elif abs_coeff == 1:
            body = mono
        else:
            body = f"{abs_coeff}*{mono}"
        if not chunks:
            chunks.append(body if coeff > 0 else f"-{body}")
        else:
            chunks.append(f" + {body}" if coeff > 0 else f" - {body}")
    return "".join(chunks) if chunks else "0"


def composition_count(m: int) -> dict:
    if m == 1:
        return {
            "m": 1,
            "total_compositions": 1,
            "positive_sign_compositions": 1,
            "negative_sign_compositions": 0,
            "mixed_sign": False,
        }
    positive = 0
    negative = 0
    for length in range(1, m + 1):
        count = comb(m - 1, length - 1)
        if (m - length) % 2 == 0:
            positive += count
        else:
            negative += count
    return {
        "m": m,
        "total_compositions": 2 ** (m - 1),
        "positive_sign_compositions": positive,
        "negative_sign_compositions": negative,
        "mixed_sign": positive > 0 and negative > 0,
    }


def validate_ref(section: str, ref: str) -> list[NonpowerFunctionalIssue]:
    if ref.startswith(("http://", "https://")):
        return []
    if (REPO_ROOT / ref).exists():
        return []
    return [issue(section, "missing-ref", ref)]


def validate_top_level(scout: dict) -> list[NonpowerFunctionalIssue]:
    issues: list[NonpowerFunctionalIssue] = []
    if scout.get("kind") != "jensen_window_pf_nonpower_functional_low_degree_scout":
        issues.append(issue("<scout>", "bad-kind", repr(scout.get("kind"))))
    expected_parents = {
        "parent_ansatz_matrix": "jensen_window_pf_nonordinary_positive_transform_ansatz_matrix",
        "parent_ansatz_row": "npt_04_nonpower_positive_functional",
        "parent_positive_readout_target": "jensen_window_pf_positive_readout_theorem_target",
        "parent_target": "jwpf_06_sign_regular_to_jensen_pf_conversion",
    }
    for key, value in expected_parents.items():
        if scout.get(key) != value:
            issues.append(issue("<scout>", f"bad-{key}", f"{scout.get(key)!r} != {value!r}"))
    boundary = str(scout.get("proof_boundary", "")).lower()
    for required in ("exact low-degree", "does not construct", "does not prove coefficient positivity", "does not prove lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<scout>", "weak-proof-boundary", required))
    objects = scout.get("objects", {})
    for required in ("H(t)", "E(t)", "mu_m", "recurrence", "composition_expansion", "normalized_window_coefficients"):
        if required not in objects or not str(objects.get(required, "")).strip():
            issues.append(issue("objects", "missing-object", required))
    return issues


def validate_formulas(scout: dict) -> list[NonpowerFunctionalIssue]:
    issues: list[NonpowerFunctionalIssue] = []
    stored = scout.get("reciprocal_polynomials", {})
    polys = reciprocal_polynomials(MAX_FORMULA_M)
    for m, poly in enumerate(polys):
        key = f"mu_{m}"
        computed = poly_string(poly)
        if stored.get(key) != computed:
            issues.append(issue("reciprocal_polynomials", f"bad-{key}", f"{stored.get(key)!r} != {computed!r}"))
    return issues


def validate_composition_counts(scout: dict) -> list[NonpowerFunctionalIssue]:
    issues: list[NonpowerFunctionalIssue] = []
    rows = scout.get("composition_counts", [])
    if not isinstance(rows, list):
        return [issue("composition_counts", "bad-type", repr(type(rows)))]
    by_m = {row.get("m"): row for row in rows if isinstance(row, dict)}
    for m in range(1, MAX_COMPOSITION_M + 1):
        expected = composition_count(m)
        actual = by_m.get(m)
        if actual != expected:
            issues.append(issue("composition_counts", f"bad-m={m}", f"{actual!r} != {expected!r}"))
    return issues


def validate_rows(scout: dict) -> tuple[list[NonpowerFunctionalIssue], int, int, int]:
    rows = scout.get("scout_rows", [])
    issues: list[NonpowerFunctionalIssue] = []
    if not isinstance(rows, list) or not rows:
        return [issue("scout_rows", "missing-scout-rows", repr(rows))], 0, 0, 0
    seen: set[str] = set()
    ready_count = 0
    live_count = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("scout_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        if row_id in seen:
            issues.append(issue(row_id, "duplicate-row", row_id))
        seen.add(row_id)
        for key in ("id", "verdict", "readiness", "source_artifacts", "finding", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        verdict = row.get("verdict")
        if verdict not in ALLOWED_VERDICTS:
            issues.append(issue(row_id, "bad-verdict", repr(verdict)))
        if row.get("readiness") != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        else:
            pass
        if row.get("readiness") == "ready_to_apply":
            ready_count += 1
        if verdict == "live_if_cone_and_functional_constructed":
            live_count += 1
        refs = row.get("source_artifacts", [])
        if not isinstance(refs, list) or not refs:
            issues.append(issue(row_id, "missing-source-artifacts", repr(refs)))
        else:
            for ref in refs:
                if isinstance(ref, str):
                    issues.extend(validate_ref(row_id, ref))
                else:
                    issues.append(issue(row_id, "bad-ref", repr(ref)))
        text = " ".join(str(row.get(key, "")) for key in ("finding", "proof_boundary")).lower()
        if str(verdict).startswith("rejected") and not any(word in text for word in ("reject", "insufficient", "not a proof", "tautological")):
            issues.append(issue(row_id, "rejected-row-lacks-rejection", text))
        if verdict == "live_if_cone_and_functional_constructed":
            for required in ("cone", "basis", "functional", "not constructed"):
                if required not in text:
                    issues.append(issue(row_id, "live-contract-missing-text", required))
        for forbidden in ("therefore rh", "we have proved lambda <= 0", "jwpf_06 is proved", "positive functional is constructed"):
            if forbidden in text:
                issues.append(issue(row_id, "forbidden-overclaim", forbidden))
    for missing in sorted(REQUIRED_IDS - seen):
        issues.append(issue(missing, "missing-scout-row", missing))
    return issues, len(rows), ready_count, live_count


def validate_kill_gates(scout: dict) -> list[NonpowerFunctionalIssue]:
    gates = scout.get("kill_gates", [])
    issues: list[NonpowerFunctionalIssue] = []
    if not isinstance(gates, list) or len(gates) < 6:
        return [issue("kill_gates", "too-few-kill-gates", repr(gates))]
    text = " ".join(str(item) for item in gates).lower()
    for required in ("independent nonnegativity", "signed composition", "unknown target", "l(k_m)=mu_m", "finite", "lambda <= 0"):
        if required not in text:
            issues.append(issue("kill_gates", "missing-kill-gate-text", required))
    return issues


def validate_summary(scout: dict, rows: int, ready_count: int, live_count: int) -> list[NonpowerFunctionalIssue]:
    issues: list[NonpowerFunctionalIssue] = []
    summary = scout.get("summary", {})
    expected = {
        "scout_rows": 7,
        "exact_identity_rows": 2,
        "rejected_rows": 2,
        "low_degree_requirement_rows": 2,
        "live_contract_rows": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if rows != 7:
        issues.append(issue("summary", "bad-row-count", str(rows)))
    if ready_count != 0:
        issues.append(issue("summary", "ready-row-present", str(ready_count)))
    if live_count != 1:
        issues.append(issue("summary", "bad-live-count", str(live_count)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("signed reciprocal", "mu_2=", "mu_3=", "no positive cone", "basis", "functional"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in scout.get("invariants", [])).lower()
    for required in ("exact algebra only", "no scout row", "positive cone", "delta_2=-g_2"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[NonpowerFunctionalIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[NonpowerFunctionalIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "jwpf_06 is proved",
        "positive functional is constructed",
        "finite rows prove",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(scout_path: Path, note_path: Path) -> tuple[list[NonpowerFunctionalIssue], int, int, int]:
    scout = load_json(scout_path)
    issues: list[NonpowerFunctionalIssue] = []
    issues.extend(validate_top_level(scout))
    issues.extend(validate_formulas(scout))
    issues.extend(validate_composition_counts(scout))
    row_issues, rows, ready_count, live_count = validate_rows(scout)
    issues.extend(row_issues)
    issues.extend(validate_kill_gates(scout))
    issues.extend(validate_summary(scout, rows, ready_count, live_count))
    issues.extend(validate_note(note_path))
    return issues, rows, ready_count, live_count


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scout", type=Path, default=DEFAULT_SCOUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, rows, ready_count, live_count = validate(args.scout, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "scout_rows": rows,
                    "ready_to_apply_rows": ready_count,
                    "live_contract_rows": live_count,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-NONPOWER-FUNCTIONAL {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF nonpower functional low-degree scout: "
            f"{rows} scout rows, {len(issues)} issues, {ready_count} ready-to-apply rows, "
            f"{live_count} live contract rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
