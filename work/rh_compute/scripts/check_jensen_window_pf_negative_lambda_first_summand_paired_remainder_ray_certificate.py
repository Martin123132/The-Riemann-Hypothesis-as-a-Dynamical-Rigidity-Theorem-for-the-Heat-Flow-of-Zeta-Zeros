#!/usr/bin/env python3
"""Validate the analytic paired-remainder ray and global closure certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
from fractions import Fraction
import json
import math
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    Q_FLOOR,
    REPO_ROOT,
    build_artifact,
    tail_polynomial_ratio,
)


REQUIRED_ROW_IDS = {
    "fsprrc_01_mode_and_derivative_geometry",
    "fsprrc_02_adaptive_window",
    "fsprrc_03_first_order_weight_bound",
    "fsprrc_04_two_tail_bound",
    "fsprrc_05_raw_moment_comparison",
    "fsprrc_06_cumulant_composition",
    "fsprrc_07_ray_remainder_theorem",
    "fsprrc_08_global_remainder_closure",
    "fsprrc_09_cumulant_wall_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda First-Summand Paired-Remainder Ray Certificate",
    "Status: analytic asymptotic-ray theorem and global paired-remainder closure",
    "Artifact kind: `jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.py",
    "Y=sqrt(8*log(q))",
    "|alpha|<=2/sqrt(q)",
    "|kappa_3(Y)+alpha|<=120/q",
    "H_t>=-299/25000>-3/250>-79/1000, u>=5",
    "H_t>=-79/1000 for every real t>=318",
    "kappa_3,t(2*log(U))>=-37/(50*t^2) for every real t>=318",
    "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md",
    "outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class CertificateIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: object) -> CertificateIssue:
    return CertificateIssue(section=section, issue=name, detail=str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(ref: object) -> list[CertificateIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue("artifact", "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue("artifact", "missing-ref", ref)]
    return []


def gaussian_expectation(expr: sp.Expr, y: sp.Symbol) -> sp.Expr:
    total = sp.Integer(0)
    for (power,), coefficient in sp.Poly(sp.expand(expr), y).terms():
        if power % 2 == 0:
            total += coefficient * (sp.factorial2(power - 1) if power else 1)
    return sp.expand(total)


def symbolic_issues() -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    u, v, y, alpha = sp.symbols("u v y alpha", positive=True)
    q_polynomial = sp.Integer(1)
    q_bounds = [2, 5, 10, 22]
    for order, bound in enumerate(q_bounds, start=1):
        q_polynomial = sp.expand(u * sp.diff(q_polynomial, u) / 2 + 2 * u * q_polynomial)
        inverse = sp.Poly(sp.expand((q_polynomial / u**order).subs(u, 1 / v)), v)
        if any(coefficient < 0 for coefficient in inverse.all_coeffs()):
            findings.append(issue("symbolic", f"nonmonotone-q{order}", inverse))
        if sp.ceiling(q_polynomial.subs(u, 5) / 5**order) != bound:
            findings.append(issue("symbolic", f"bad-q{order}-bound", q_polynomial))
    bell_v4 = sum(
        2**parts
        * sp.factorial(parts - 1)
        * sp.bell(4, parts, tuple(q_bounds[: 5 - parts]))
        for parts in range(1, 5)
    )
    if bell_v4 != 4_120:
        findings.append(issue("symbolic", "bad-V4-bell-constant", bell_v4))

    linear_weight = 1 - alpha * y**3 / 6
    moments = [gaussian_expectation(y**order * linear_weight, y) for order in range(4)]
    expected_moments = [1, -alpha / 2, 1, -5 * alpha / 2]
    if any(sp.expand(actual - target) != 0 for actual, target in zip(moments, expected_moments)):
        findings.append(issue("symbolic", "bad-linear-model-moments", moments))
    base_cumulant = moments[3] - 3 * moments[1] * moments[2] + 2 * moments[1] ** 3
    if sp.expand(base_cumulant + alpha + alpha**3 / 4) != 0:
        findings.append(issue("symbolic", "bad-base-cumulant", base_cumulant))

    if not (
        Fraction(4, 1) * Fraction(39, 10) ** 3 > 14**2
        and Fraction(81, 1) / Fraction(39, 10) ** 3 < Fraction(36, 25)
        and Fraction(31, 1) / Fraction(39, 10) ** 2 < 3
    ):
        findings.append(issue("symbolic", "bad-derivative-ratio-gate", "exact fraction failed"))

    actual_tail = [tail_polynomial_ratio(order, Fraction(9, 10)) for order in range(4)]
    model_tail = [tail_polynomial_ratio(order, Fraction(1, 1)) for order in range(7)]
    if max(actual_tail) >= 2 or max(model_tail) >= 2:
        findings.append(issue("symbolic", "bad-tail-polynomial-gate", (actual_tail, model_tail)))
    return findings


def source_composition_issues() -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    paths = {
        "compact": REPO_ROOT
        / "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.json",
        "leading": REPO_ROOT
        / "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.json",
        "cumulant": REPO_ROOT
        / "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.json",
    }
    try:
        sources = {key: load_json(path) for key, path in paths.items()}
    except Exception as exc:
        return [issue("sources", "load-failed", exc)]
    compact = sources["compact"]
    leading = sources["leading"]
    cumulant = sources["cumulant"]
    if compact.get("status") != "interval compact paired-remainder and negative-skewness theorem with open asymptotic ray":
        findings.append(issue("sources", "bad-compact-status", compact.get("status")))
    if compact.get("summary", {}).get("compact_accepted_blocks") != 4_074:
        findings.append(issue("sources", "bad-compact-blocks", compact.get("summary")))
    if compact.get("summary", {}).get("open_ray_rows") != 1:
        findings.append(issue("sources", "bad-compact-open-ray", compact.get("summary")))
    if leading.get("summary", {}).get("positive_analytic_ray_gates") != 3:
        findings.append(issue("sources", "bad-leading-rays", leading.get("summary")))
    ray = leading.get("diagnostics", {}).get("ray", {})
    for key in (
        "cap_margin_fraction",
        "correction_cap_margin_fraction",
        "fifth_correction_cap_margin_fraction",
    ):
        try:
            if Fraction(ray[key]) <= 0:
                findings.append(issue("sources", f"bad-leading-{key}", ray[key]))
        except Exception as exc:
            findings.append(issue("sources", f"missing-leading-{key}", exc))
    if not (
        cumulant.get("summary", {}).get("exact_identity_rows") == 4
        and cumulant.get("summary", {}).get("conditional_bridge_rows") == 1
    ):
        findings.append(issue("sources", "bad-cumulant-source", cumulant.get("summary")))
    if Fraction(13, 20) + Fraction(1, 100) + Fraction(1, 1000) + Fraction(79, 1000) != Fraction(37, 50):
        findings.append(issue("sources", "bad-global-cumulant-composition", "fraction identity failed"))
    return findings


def exact_budget_issues(diagnostics: dict) -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    sqrt_bound = 30_000
    denominator = 1 - Fraction(6, Q_FLOOR)
    raw_coefficients = {
        "raw1": (Fraction(12, 1) + Fraction(6, sqrt_bound)) / denominator,
        "raw2": Fraction(35, 1) / denominator,
        "raw3": (Fraction(79, 1) + Fraction(30, sqrt_bound)) / denominator,
    }
    raw_caps = {"raw1": 13, "raw2": 36, "raw3": 80}
    if not all(raw_coefficients[key] < cap for key, cap in raw_caps.items()):
        findings.append(issue("budget", "bad-raw-errors", raw_coefficients))
    cumulant = (
        Fraction(2, sqrt_bound)
        + 80
        + 3 * (13 * (1 + Fraction(36, Q_FLOOR)) + Fraction(36, sqrt_bound))
        + 6 * 13 * (Fraction(1, sqrt_bound) + Fraction(13, Q_FLOOR)) ** 2
    )
    if not cumulant < 120:
        findings.append(issue("budget", "bad-cumulant-coefficient", cumulant))
    scaled = Fraction(6, 25) * 120 / sqrt_bound
    ray_cap = scaled + Fraction(1, 100) + Fraction(1, 1000)
    if scaled != Fraction(12, 12_500) or ray_cap != Fraction(299, 25_000):
        findings.append(issue("budget", "bad-exact-ray-cap", (scaled, ray_cap)))
    if not ray_cap < Fraction(3, 250) < Fraction(79, 1000):
        findings.append(issue("budget", "weak-ray-cap", ray_cap))

    stored = diagnostics.get("cumulant_composition", {})
    theorem = diagnostics.get("ray_theorem", {})
    expected = {
        "coefficient_upper": str(cumulant),
        "scaled_cap_fraction": str(scaled),
    }
    for key, value in expected.items():
        if stored.get(key) != value:
            findings.append(issue("budget", f"bad-stored-{key}", stored.get(key)))
    if theorem.get("ray_remainder_absolute_cap") != str(ray_cap):
        findings.append(issue("budget", "bad-stored-ray-cap", theorem))
    if theorem.get("floor_margin") != str(Fraction(79, 1000) - ray_cap):
        findings.append(issue("budget", "bad-stored-floor-margin", theorem))
    return findings


def validate_artifact(artifact: dict) -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != "analytic asymptotic-ray theorem and global paired-remainder closure":
        findings.append(issue("artifact", "bad-status", artifact.get("status")))
    if artifact.get("date") != "2026-07-10":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    for key in (
        "source_compact_certificate",
        "source_leading_certificate",
        "source_cumulant_bridge",
        "source_saddle_target",
        "generator",
        "checker",
    ):
        findings.extend(validate_ref(artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for text in ("proves", "does not prove", "u>=5", "global", "rh", "lambda <= 0"):
        if text not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", text))

    summary = artifact.get("summary", {})
    expected_summary = {
        "certificate_rows": 9,
        "analytic_ray_rows": 1,
        "global_remainder_closure_rows": 1,
        "open_ray_rows": 0,
        "ready_to_apply_rows": 2,
        "ray_remainder_lower": "-299/25000",
        "required_remainder_floor": "-79/1000",
        "floor_margin": "419/6250",
    }
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            findings.append(issue("summary", f"bad-{key}", summary.get(key)))

    diagnostics = artifact.get("diagnostics", {})
    try:
        handoff = diagnostics["mode_handoff"]
        central = diagnostics["central_window"]
        if not Decimal(handoff["potential_slope_at_compact_start_upper"]) < 318:
            findings.append(issue("handoff", "bad-compact-start", handoff))
        if not Decimal(handoff["q_at_499_lower"]) > Decimal("1e9"):
            findings.append(issue("handoff", "bad-window-q", handoff))
        if not Decimal("0.022") < Decimal(central["central_remainder_upper"]) < Decimal("0.023"):
            findings.append(issue("central", "bad-remainder", central))
        if not Decimal(central["endpoint_slope_loss_upper"]) < Decimal("0.001"):
            findings.append(issue("central", "bad-slope-loss", central))
        if not Decimal(central["standardized_log_shift_upper"]) < Decimal("0.000025"):
            findings.append(issue("central", "bad-log-shift", central))
        if not Decimal(central["q_log_drift_upper"]) < Decimal("0.0005"):
            findings.append(issue("central", "bad-q-drift", central))
        if not Decimal(central["tail_log_margin_lower"]) > 0:
            findings.append(issue("central", "bad-tail-log", central))
    except Exception as exc:
        findings.append(issue("diagnostics", "bad-decimal", exc))

    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        findings.append(issue("rows", "bad-row-ids", sorted(str(value) for value in ids)))
    if len(rows) != 9:
        findings.append(issue("rows", "bad-row-count", len(rows)))
    ready = [row.get("id") for row in rows if row.get("readiness") == "ready_to_apply"]
    if ready != ["fsprrc_08_global_remainder_closure", "fsprrc_09_cumulant_wall_handoff"]:
        findings.append(issue("rows", "bad-ready-rows", ready))

    findings.extend(symbolic_issues())
    findings.extend(source_composition_issues())
    findings.extend(exact_budget_issues(diagnostics))
    try:
        recomputed = json.loads(json.dumps(build_artifact()))
        for key in ("diagnostics", "summary", "rows"):
            if artifact.get(key) != recomputed.get(key):
                findings.append(issue("recompute", f"mismatch-{key}", "stored differs"))
    except Exception as exc:
        findings.append(issue("recompute", "exception", exc))
    return findings


def validate_note(path: Path) -> list[CertificateIssue]:
    if not path.exists():
        return [issue("note", "missing", path)]
    text = path.read_text(encoding="utf-8")
    return [issue("note", "missing-text", value) for value in REQUIRED_NOTE_STRINGS if value not in text]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    findings: list[CertificateIssue] = []
    try:
        artifact = load_json(args.artifact)
    except Exception as exc:
        artifact = {}
        findings.append(issue("artifact", "load-failed", exc))
    if artifact:
        findings.extend(validate_artifact(artifact))
    findings.extend(validate_note(args.note))
    ok = not findings
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in findings]}, indent=2, sort_keys=True))
    else:
        for finding in findings:
            print(f"ISSUE {finding.section} [{finding.issue}] {finding.detail}")
        print(
            "validated Jensen-window PF negative-lambda first-summand paired-remainder ray certificate: "
            f"9 rows, {len(findings)} issues, 1 analytic ray theorem, "
            "1 global remainder closure, 0 open rays, 2 ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
