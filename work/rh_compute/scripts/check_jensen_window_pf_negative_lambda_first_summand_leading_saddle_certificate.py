#!/usr/bin/env python3
"""Validate the first-summand leading-saddle interval certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_SOURCE_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "fslsc_01_log_variable_potential",
    "fslsc_02_mode_localization",
    "fslsc_03_compact_leading_cap",
    "fslsc_04_ray_leading_cap",
    "fslsc_05_global_leading_saddle_theorem",
    "fslsc_06_remainder_reduction",
    "fslsc_07_remainder_scout",
    "fslsc_08_full_wall_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda First-Summand Leading-Saddle Certificate",
    "Status: interval leading-saddle theorem and open remainder target",
    "Artifact kind: `jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.py",
    "validated Jensen-window PF negative-lambda first-summand leading-saddle certificate: 8 rows, 0 issues, 40740 positive leading intervals, 40740 positive cubic-correction intervals, 40740 positive fifth-correction intervals, 3 positive analytic ray gates, 9 positive seventh-remainder samples, 1 open remainder, 0 ready-to-apply rows",
    "V(x)=100*u^2+q-5*u-log(2*q-3)-log(u)",
    "t(u)^2*V'''(u)/V''(u)^3<=13/20",
    "V''>=(39/10)*u^2*q",
    ">=-79/1000, t>=318",
    "outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md",
    "outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md",
    "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md",
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


def symbolic_derivative_issues() -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    u, q = sp.symbols("u q", positive=True)

    def derivative_x(value: sp.Expr) -> sp.Expr:
        return sp.factor(u * sp.diff(value, u) / 2 + 2 * u * q * sp.diff(value, q))

    potential = 100 * u**2 + q - 5 * u - sp.log(2 * q - 3) - sp.log(u)
    derived = []
    current = potential
    for _ in range(7):
        current = derivative_x(current)
        derived.append(current)
    expected = [
        (
            8 * q**2 * u
            + 400 * q * u**2
            - 30 * q * u
            - 2 * q
            - 600 * u**2
            + 15 * u
            + 3
        )
        / (2 * (2 * q - 3)),
        u
        * (
            64 * q**3 * u
            + 16 * q**3
            + 1408 * q**2 * u
            - 84 * q**2
            - 4560 * q * u
            + 120 * q
            + 3600 * u
            - 45
        )
        / (4 * (2 * q - 3) ** 2),
        u
        * (
            512 * q**4 * u**2
            + 384 * q**4 * u
            + 32 * q**4
            - 2304 * q**3 * u**2
            + 4672 * q**3 * u
            - 216 * q**3
            + 2688 * q**2 * u**2
            - 25632 * q**2 * u
            + 492 * q**2
            - 2880 * q * u**2
            + 41040 * q * u
            - 450 * q
            - 21600 * u
            + 135
        )
        / (8 * (2 * q - 3) ** 3),
    ]
    r = 2 / (2 * q - 3)
    q1 = q * (2 * u)
    q2 = q * u * (4 * u + 1)
    q3 = q * u * (16 * u**2 + 12 * u + 1) / 2
    q4 = q * u * (64 * u**3 + 96 * u**2 + 28 * u + 1) / 4
    q5 = q * u * (256 * u**4 + 640 * u**3 + 400 * u**2 + 60 * u + 1) / 8
    q6 = q * u * (1024 * u**5 + 3840 * u**4 + 4160 * u**3 + 1440 * u**2 + 124 * u + 1) / 16
    q7 = q * u * (4096 * u**6 + 21504 * u**5 + 35840 * u**4 + 22400 * u**3 + 4816 * u**2 + 252 * u + 1) / 32
    g1, g2, g3, g4, g5, g6, g7 = (
        1 - r,
        r**2,
        -2 * r**3,
        6 * r**4,
        -24 * r**5,
        120 * r**6,
        -720 * r**7,
    )
    expected.extend(
        [
            100 * u**2
            + g1 * q4
            + g2 * (4 * q1 * q3 + 3 * q2**2)
            + 6 * g3 * q1**2 * q2
            + g4 * q1**4
            - 5 * u / 16,
            100 * u**2
            + g1 * q5
            + g2 * (5 * q1 * q4 + 10 * q2 * q3)
            + g3 * (10 * q1**2 * q3 + 15 * q1 * q2**2)
            + 10 * g4 * q1**3 * q2
            + g5 * q1**5
            - 5 * u / 32,
        ]
    )
    q_derivatives = [q1, q2, q3, q4, q5, q6, q7]
    g_derivatives = [g1, g2, g3, g4, g5, g6, g7]
    for n in (6, 7):
        composition = sum(
            g_derivatives[m - 1]
            * sp.bell(n, m, tuple(q_derivatives[: n - m + 1]))
            for m in range(1, n + 1)
        )
        expected.append(100 * u**2 + composition - 5 * u / (2**n))
    for order, (actual, target) in enumerate(zip(derived, expected), start=1):
        if sp.cancel(actual - target) != 0:
            findings.append(issue("symbolic", f"bad-V{order}", sp.factor(actual - target)))
    return findings


def validate_artifact(artifact: dict, source_path: Path) -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != "interval leading-saddle theorem and open remainder target":
        findings.append(issue("artifact", "bad-status", artifact.get("status")))
    if artifact.get("date") != "2026-07-10":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    for key in (
        "source_cumulant_bridge",
        "source_saddle_target",
        "source_json",
        "generator",
        "checker",
    ):
        findings.extend(validate_ref(artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for text in ("proves", "does not prove", "remainder", "cumulant", "rh", "lambda <= 0"):
        if text not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", text))

    summary = artifact.get("summary", {})
    expected_summary = {
        "certificate_rows": 8,
        "compact_subintervals": 40_740,
        "positive_compact_subintervals": 40_740,
        "positive_correction_compact_subintervals": 40_740,
        "positive_fifth_correction_compact_subintervals": 40_740,
        "positive_curvature_subintervals": 40_740,
        "positive_analytic_ray_gates": 3,
        "sample_rows": 9,
        "positive_remainder_sample_rows": 9,
        "minimum_remainder_sample_margin_at_t": 100_000,
        "open_remainder_rows": 1,
        "ready_to_apply_rows": 0,
    }
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            findings.append(issue("summary", f"bad-{key}", summary.get(key)))
    try:
        if not Decimal("0.079") < Decimal(summary["minimum_remainder_sample_margin"]) < Decimal("0.081"):
            findings.append(issue("summary", "bad-remainder-margin", summary))
    except Exception as exc:
        findings.append(issue("summary", "bad-remainder-decimal", exc))

    diagnostics = artifact.get("diagnostics", {})
    compact = diagnostics.get("compact", {})
    localization = diagnostics.get("mode_localization", {})
    ray = diagnostics.get("ray", {})
    try:
        if not Decimal(compact["minimum_curvature_lower"]) > Decimal(600):
            findings.append(issue("compact", "weak-curvature", compact))
        if not Decimal(compact["minimum_leading_cap_margin_lower"]) > Decimal("0.02"):
            findings.append(issue("compact", "weak-cap-margin", compact))
        if not Decimal(compact["minimum_correction_cap_margin_lower"]) > Decimal("0.005"):
            findings.append(issue("compact", "weak-correction-margin", compact))
        if not Decimal(compact["minimum_fifth_correction_cap_margin_lower"]) > Decimal("0.0008"):
            findings.append(issue("compact", "weak-fifth-margin", compact))
        if not Decimal(localization["original_u_log_integrand_slope_at_t318_lower"]) > Decimal("0.3"):
            findings.append(issue("localization", "weak-slope", localization))
        if not Decimal(localization["potential_slope_at_u_start_upper"]) < Decimal(318):
            findings.append(issue("localization", "bad-potential-slope", localization))
        if ray.get("cap_margin_fraction") != "12561/43940":
            findings.append(issue("ray", "bad-cap-margin", ray))
        if not Decimal(ray["q_at_5_minus_500_lower"]) > 0:
            findings.append(issue("ray", "weak-q-endpoint", ray))
        if ray.get("correction_cap_margin_fraction") != "390971184173/39097152900000":
            findings.append(issue("ray", "bad-correction-margin", ray))
        if not Decimal(ray["q_at_5_minus_billion_lower"]) > 0:
            findings.append(issue("ray", "weak-billion-endpoint", ray))
        if ray.get("fifth_correction_cap_margin_fraction") != "9276816051497935508833/9276816051500400000000000":
            findings.append(issue("ray", "bad-fifth-margin", ray))
    except Exception as exc:
        findings.append(issue("diagnostics", "bad-decimal", exc))

    scout_rows = diagnostics.get("remainder_sample_rows", [])
    if len(scout_rows) != 9 or not all(row.get("above_remainder_floor") for row in scout_rows):
        findings.append(issue("scout", "bad-remainder-rows", scout_rows))
    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        findings.append(issue("rows", "bad-row-ids", sorted(str(value) for value in ids)))
    if len(rows) != 8:
        findings.append(issue("rows", "bad-row-count", len(rows)))

    findings.extend(symbolic_derivative_issues())
    try:
        recomputed = build_artifact(source_path)
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
    parser.add_argument("--source-json", type=Path, default=DEFAULT_SOURCE_JSON)
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
        findings.extend(validate_artifact(artifact, args.source_json))
    findings.extend(validate_note(args.note))
    ok = not findings
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in findings]}, indent=2, sort_keys=True))
    else:
        for finding in findings:
            print(f"ISSUE {finding.section} [{finding.issue}] {finding.detail}")
        print(
            "validated Jensen-window PF negative-lambda first-summand leading-saddle certificate: "
            f"8 rows, {len(findings)} issues, 40740 positive leading intervals, "
            "40740 positive cubic-correction intervals, 40740 positive fifth-correction "
            "intervals, 3 positive analytic ray gates, 9 positive seventh-remainder "
            "samples, 1 open remainder, 0 ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
