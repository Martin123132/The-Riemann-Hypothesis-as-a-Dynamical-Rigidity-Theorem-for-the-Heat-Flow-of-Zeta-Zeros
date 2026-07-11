#!/usr/bin/env python3
"""Validate the first-summand paired seventh-order remainder certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "fsprc_01_standardized_pairing",
    "fsprc_02_left_tail_monotonicity",
    "fsprc_03_eighth_derivative_envelope",
    "fsprc_04_paired_midpoint_rule",
    "fsprc_05_two_tail_budget",
    "fsprc_06_compact_remainder_theorem",
    "fsprc_07_far_ray_target",
    "fsprc_08_wall_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda First-Summand Paired-Remainder Certificate",
    "Status: interval compact paired-remainder and negative-skewness theorem with open asymptotic ray",
    "Artifact kind: `jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.py",
    "H_t=(t^2/a_t^(3/2))*(kappa_3(Y_t)+alpha_t+C3_t+C5_t)",
    "sup_|y|<=6 V^(8)(x_t+y/sqrt(a_t))/a_t^4<=1/50000",
    "H_t>=-79/1000 for every real mode parameter 0.9264<=u_t<=5",
    "kappa_3,t(2*log(U))<0 for every real mode parameter 0.9264<=u_t<=5",
    "Minimum outward-rounded margin",
    "H_t>=-79/1000 for u_t>=5",
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


def series_multiply(left: list[sp.Expr], right: list[sp.Expr]) -> list[sp.Expr]:
    degree = len(left) - 1
    return [
        sp.expand(sum(left[index] * right[order - index] for index in range(order + 1)))
        for order in range(degree + 1)
    ]


def series_divide(numerator: list[sp.Expr], denominator: list[sp.Expr]) -> list[sp.Expr]:
    quotient: list[sp.Expr] = []
    for order, coefficient in enumerate(numerator):
        correction = sum(
            denominator[index] * quotient[order - index]
            for index in range(1, order + 1)
        )
        quotient.append(sp.expand(coefficient - correction))
    return quotient


def symbolic_edgeworth_issues() -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    y = sp.symbols("y")
    alpha, beta, gamma, delta, epsilon = sp.symbols("alpha beta gamma delta epsilon")
    degree = 5
    perturbation = [
        sp.Integer(0),
        alpha * y**3 / sp.factorial(3),
        beta * y**4 / sp.factorial(4),
        gamma * y**5 / sp.factorial(5),
        delta * y**6 / sp.factorial(6),
        epsilon * y**7 / sp.factorial(7),
    ]
    weight = [sp.Integer(1)]
    for order in range(1, degree + 1):
        coefficient = -sum(
            index * perturbation[index] * weight[order - index]
            for index in range(1, order + 1)
        ) / order
        weight.append(sp.expand(coefficient))
    moments = [
        [gaussian_expectation(y**power * weight[order], y) for order in range(degree + 1)]
        for power in range(4)
    ]
    raw = [series_divide(moment, moments[0]) for moment in moments]
    raw1_raw2 = series_multiply(raw[1], raw[2])
    raw1_cubed = series_multiply(series_multiply(raw[1], raw[1]), raw[1])
    cumulant = [
        sp.expand(raw[3][order] - 3 * raw1_raw2[order] + 2 * raw1_cubed[order])
        for order in range(degree + 1)
    ]
    expected = [
        0,
        -alpha,
        0,
        -(8 * alpha**3 - 7 * alpha * beta + gamma) / 2,
        0,
        -(
            525 * alpha**5
            - 954 * alpha**3 * beta
            + 234 * alpha**2 * gamma
            + 298 * alpha * beta**2
            - 37 * alpha * delta
            - 57 * beta * gamma
            + 3 * epsilon
        )
        / 24,
    ]
    for order, (actual, target) in enumerate(zip(cumulant, expected)):
        if sp.expand(actual - target) != 0:
            findings.append(issue("symbolic", f"bad-cumulant-order-{order}", actual - target))
    return findings


def symbolic_ray_issues() -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    u, v = sp.symbols("u v", positive=True)
    q_polynomial = sp.Integer(1)
    q_bounds = [2, 5, 10, 22, 51, 124, 314, 822, 2218]
    for order, expected_bound in enumerate(q_bounds, start=1):
        q_polynomial = sp.expand(u * sp.diff(q_polynomial, u) / 2 + 2 * u * q_polynomial)
        inverse_form = sp.Poly(sp.expand((q_polynomial / u**order).subs(u, 1 / v)), v)
        if any(coefficient < 0 for coefficient in inverse_form.all_coeffs()):
            findings.append(issue("symbolic", f"nonmonotone-q{order}-bound", inverse_form))
        endpoint = sp.cancel(q_polynomial.subs(u, 5) / 5**order)
        if sp.ceiling(endpoint) != expected_bound:
            findings.append(issue("symbolic", f"bad-q{order}-bound", endpoint))
    bell_constants = {}
    for order in (8, 9):
        bell_constants[order] = sum(
            2**parts
            * sp.factorial(parts - 1)
            * sp.bell(order, parts, tuple(q_bounds[: order - parts + 1]))
            for parts in range(1, order + 1)
        )
    expected_constants = {8: 2_376_267_836, 9: 97_312_467_060}
    if bell_constants != expected_constants:
        findings.append(issue("symbolic", "bad-bell-constants", bell_constants))
    if 822 * 1_000_000_000 + bell_constants[8] >= 1_000 * 1_000_000_000:
        findings.append(issue("symbolic", "bad-V8-ray-majorant", bell_constants[8]))
    if 512 * 1_000_000_000 - bell_constants[9] <= 0:
        findings.append(issue("symbolic", "bad-V9-ray-margin", bell_constants[9]))
    return findings


def normalized_recompute(value: dict) -> dict:
    copy = json.loads(json.dumps(value))

    def remove_elapsed(item: object) -> None:
        if isinstance(item, dict):
            item.pop("elapsed_seconds", None)
            for child in item.values():
                remove_elapsed(child)
        elif isinstance(item, list):
            for child in item:
                remove_elapsed(child)

    remove_elapsed(copy)
    return copy


def validate_artifact(artifact: dict) -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != "interval compact paired-remainder and negative-skewness theorem with open asymptotic ray":
        findings.append(issue("artifact", "bad-status", artifact.get("status")))
    if artifact.get("date") != "2026-07-10":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    for key in ("source_saddle_expansion", "source_cumulant_bridge", "generator", "checker"):
        findings.extend(validate_ref(artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for text in ("proves", "does not prove", "u>=5", "cumulant", "cone entry", "rh", "lambda <= 0"):
        if text not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", text))

    summary = artifact.get("summary", {})
    expected_summary = {
        "certificate_rows": 8,
        "positive_curvature_intervals": 9_164,
        "positive_v9_intervals": 42_800,
        "positive_eighth_envelope_intervals": 40_736,
        "compact_base_blocks": 4_074,
        "compact_accepted_blocks": 4_074,
        "compact_negative_cumulant_blocks": 4_074,
        "open_ray_rows": 1,
        "ready_to_apply_rows": 0,
    }
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            findings.append(issue("summary", f"bad-{key}", summary.get(key)))
    try:
        margin = Decimal(summary["compact_minimum_margin_lower"])
        if not Decimal("0.00006") < margin < Decimal("0.00007"):
            findings.append(issue("summary", "bad-compact-margin", margin))
    except Exception as exc:
        findings.append(issue("summary", "bad-compact-decimal", exc))

    diagnostics = artifact.get("diagnostics", {})
    parameters = diagnostics.get("parameters", {})
    expected_parameters = {
        "precision_bits": 192,
        "mode_start": "579/625",
        "mode_end": "5",
        "base_width": "1/1000",
        "maximum_refinement_depth": 4,
        "primary_panels": 300,
        "mid_panels": 600,
        "retry_panels": 900,
        "high_panel_start": "19/5",
        "window_y": 6,
        "eighth_envelope": "1/50000",
        "remainder_floor": "-79/1000",
    }
    if parameters != expected_parameters:
        findings.append(issue("parameters", "mismatch", parameters))

    geometry = diagnostics.get("geometry", {})
    eighth = diagnostics.get("eighth_derivative_envelope", {})
    compact = diagnostics.get("compact_remainder", {})
    try:
        if not Decimal(geometry["minimum_curvature_lower"]) > Decimal("0.008"):
            findings.append(issue("geometry", "weak-curvature", geometry))
        if not Decimal(geometry["tiny_u_q_upper"]) < Decimal("3.3"):
            findings.append(issue("geometry", "weak-tiny-q", geometry))
        for key in ("tiny_u_crude_slope_upper", "tiny_u_slope_upper"):
            if not Decimal(geometry[key]) < 0:
                findings.append(issue("geometry", f"bad-{key}", geometry[key]))
        if not Decimal(geometry["left_reference_slope_lower"]) > 199:
            findings.append(issue("geometry", "weak-reference-slope", geometry))
        if not Decimal(geometry["minimum_v9_lower"]) > Decimal("1e6"):
            findings.append(issue("geometry", "weak-V9", geometry))
        if not Decimal(geometry["v9_ray_q_lower"]) > Decimal("1e9"):
            findings.append(issue("geometry", "weak-ray-q", geometry))
        if not Decimal(geometry["v9_ray_margin_lower"]) > Decimal("4e11"):
            findings.append(issue("geometry", "weak-ray-V9", geometry))
        if "97312467060" not in geometry.get("v9_ray_argument", ""):
            findings.append(issue("geometry", "stale-ray-constant", geometry.get("v9_ray_argument")))

        if not Decimal(eighth["minimum_margin_lower"]) > Decimal("9e-7"):
            findings.append(issue("eighth", "weak-margin", eighth))
        if not Decimal(eighth["minimum_window_u_lower"]) > Decimal("0.82"):
            findings.append(issue("eighth", "bad-left-window", eighth))
        if not Decimal(eighth["maximum_window_u_upper"]) < Decimal("5.001"):
            findings.append(issue("eighth", "bad-right-window", eighth))
        ray_margin = Fraction(1, 50_000) - Fraction(1_000 * 2**9, 1) / (
            Fraction(39, 10) ** 4 * 1_000_000_000**3
        )
        if eighth.get("ray_margin_fraction") != str(ray_margin):
            findings.append(issue("eighth", "bad-ray-margin", eighth.get("ray_margin_fraction")))

        if compact.get("attempt_count") != 4_401:
            findings.append(issue("compact", "bad-attempt-count", compact.get("attempt_count")))
        if compact.get("maximum_refinement_depth") != 0:
            findings.append(issue("compact", "unexpected-refinement", compact))
        if compact.get("negative_cumulant_blocks") != compact.get("accepted_block_count"):
            findings.append(issue("compact", "nonnegative-cumulant-block", compact))
        if not Decimal(compact["maximum_cumulant3_upper"]) < 0:
            findings.append(issue("compact", "weak-cumulant-sign", compact))
        if not Decimal(compact["end_t_lower"]) > Decimal("1.5e10"):
            findings.append(issue("compact", "weak-end-t", compact))
        worst = compact.get("minimum_margin_block", {})
        if worst != {"left": "11681/2500", "right": "23367/5000", "depth": 0, "panels": 600}:
            findings.append(issue("compact", "bad-worst-block", worst))
        selected = compact.get("selected_blocks", [])
        if len(selected) != 5 or not all(row.get("passed") and row.get("cumulant3_negative") for row in selected):
            findings.append(issue("compact", "bad-selected-blocks", selected))
    except Exception as exc:
        findings.append(issue("diagnostics", "bad-decimal", exc))

    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        findings.append(issue("rows", "bad-row-ids", sorted(str(value) for value in ids)))
    if len(rows) != 8:
        findings.append(issue("rows", "bad-row-count", len(rows)))
    open_rows = [row for row in rows if row.get("readiness") == "open_target"]
    if len(open_rows) != 1 or open_rows[0].get("id") != "fsprc_07_far_ray_target":
        findings.append(issue("rows", "bad-open-row", open_rows))

    findings.extend(symbolic_edgeworth_issues())
    findings.extend(symbolic_ray_issues())
    try:
        recomputed = build_artifact()
        stored = normalized_recompute(artifact)
        regenerated = normalized_recompute(recomputed)
        for key in ("diagnostics", "summary", "rows"):
            if stored.get(key) != regenerated.get(key):
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
            "validated Jensen-window PF negative-lambda first-summand paired-remainder certificate: "
            f"8 rows, {len(findings)} issues, 40736 eighth-envelope intervals, "
            "4074 compact remainder blocks, 1 open asymptotic ray, 0 ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
