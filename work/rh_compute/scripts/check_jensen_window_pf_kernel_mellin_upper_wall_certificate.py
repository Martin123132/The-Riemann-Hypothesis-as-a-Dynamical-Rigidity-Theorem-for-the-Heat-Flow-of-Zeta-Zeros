#!/usr/bin/env python3
"""Validate the kernel Mellin upper-wall interval theorem certificate."""

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

from jensen_window_pf_kernel_mellin_upper_wall_certificate import (  # noqa: E402
    DEFAULT_EVENNESS_JSON,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "kmuw_01_exact_y_mellin_reindexing",
    "kmuw_02_full_kernel_evenness_input",
    "kmuw_03_compact_interval_log_concavity",
    "kmuw_04_n1_analytic_ray",
    "kmuw_05_full_n_tail_ray",
    "kmuw_06_global_kernel_sqrt_log_concavity",
    "kmuw_07_berwald_borell_upper_wall",
    "kmuw_08_remaining_monotone_wall",
}

ALLOWED_ROLES = {
    "exact_reduction",
    "exact_input",
    "interval_certificate",
    "exact_analytic_bound",
    "interval_analytic_bound",
    "interval_theorem",
    "literature_theorem_application",
    "non_promotion_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Kernel Mellin Upper-Wall Certificate",
    "Date: 2026-07-10",
    "Status: interval theorem certificate. This is not a proof",
    "Artifact kind: `jensen_window_pf_kernel_mellin_upper_wall_certificate`",
    "work/rh_compute/results/jensen_window_pf_kernel_mellin_upper_wall_certificate.json",
    "python work/rh_compute/scripts/jensen_window_pf_kernel_mellin_upper_wall_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_kernel_mellin_upper_wall_certificate.py",
    (
        "validated Jensen-window PF kernel Mellin upper-wall certificate: 8 rows, 0 issues, "
        "200 positive compact intervals, 1 positive analytic ray, 1 remaining open cone clause, "
        "0 ready-to-apply rows"
    ),
    "M_k(lambda) = integral_0^infinity y^(k-1/2)*exp(lambda*y)*g(y) dy",
    "A_k(lambda) = sqrt(pi)*4^(-k)*M_k(lambda)/Gamma(k+1/2)",
    "Q(y)=g'(y)^2-g(y)*g''(y) > 0",
    "D(u)=L'(u)-u*L''(u) > 0",
    "A_k(lambda)^2 >= A_(k-1)(lambda)*A_(k+1)(lambda)",
    "x_k(lambda) <= 1",
    "`x_(k+1)>=x_k`",
    "https://doi.org/10.1007/BF01362702",
    "https://www.math.tau.ac.il/~klartagb/papers/log_concave_bernstein.pdf",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma.md",
    "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
    "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
    "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
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


def validate_ref(section: str, ref: object) -> list[CertificateIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if ref.startswith("https://"):
        return []
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_symbolic_identities() -> list[CertificateIssue]:
    u, phi, phi1, phi2 = sp.symbols("u phi phi1 phi2", nonzero=True)
    log_first = phi1 / phi
    log_second = phi2 / phi - (phi1 / phi) ** 2
    g_first = phi1 / (2 * u)
    g_second = (u * phi2 - phi1) / (4 * u**3)
    D = log_first - u * log_second
    Q = g_first**2 - phi * g_second
    checks = {
        "sqrt-log-curvature": sp.simplify(Q / phi**2 - D / (4 * u**3)),
        "sqrt-second-derivative": sp.simplify(
            (u * log_second - log_first) / (4 * u**3) + D / (4 * u**3)
        ),
    }
    k = sp.symbols("k", integer=True, positive=True)
    raw_ratio = sp.symbols("R", positive=True)
    x = (2 * k - 1) * raw_ratio / (2 * k + 1)
    checks["raw-upper-wall"] = sp.simplify(
        (1 - x)
        - (2 * k - 1)
        / (2 * k + 1)
        * ((2 * k + 1) / (2 * k - 1) - raw_ratio)
    )
    findings = []
    for name, value in checks.items():
        if value != 0:
            findings.append(issue("symbolic", f"bad-{name}", value))
    for index in range(1, 16):
        left = sp.factorial(2 * index)
        right = (
            4**index
            * sp.factorial(index)
            * sp.gamma(sp.Rational(2 * index + 1, 2))
            / sp.sqrt(sp.pi)
        )
        if sp.simplify(left - right) != 0:
            findings.append(issue("symbolic", "bad-duplication", index))
    return findings


def validate_top_level(artifact: dict) -> list[CertificateIssue]:
    findings = []
    if artifact.get("kind") != "jensen_window_pf_kernel_mellin_upper_wall_certificate":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != "interval theorem certificate":
        findings.append(issue("artifact", "bad-status", artifact.get("status")))
    if artifact.get("date") != "2026-07-10":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    for key in (
        "source_evenness_json",
        "source_evenness_note",
        "source_boundary_threshold_note",
        "source_cone_entry_target",
        "source_defect_tail_target",
        "generator",
        "checker",
    ):
        findings.extend(validate_ref("artifact", artifact.get(key)))
    for ref in artifact.get("literature_sources", []):
        findings.extend(validate_ref("literature", ref))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for text in (
        "full kernel",
        "all-k upper cone wall",
        "does not prove",
        "adjacent-k",
        "rh",
        "lambda <= 0",
    ):
        if text not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", text))
    return findings


def validate_summary(artifact: dict) -> list[CertificateIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "certificate_rows": 8,
        "ready_to_apply_rows": 0,
        "compact_subintervals": 200,
        "positive_compact_subintervals": 200,
        "global_kernel_sqrt_log_concavity_certified": True,
        "all_k_upper_wall_certified": True,
        "remaining_open_cone_clauses": 1,
    }
    findings = []
    for key, value in expected.items():
        if summary.get(key) != value:
            findings.append(issue("summary", f"bad-{key}", summary.get(key)))
    for key in ("compact_minimum_Q_lower", "ray_full_D_lower"):
        try:
            if Decimal(str(summary.get(key))) <= 0:
                findings.append(issue("summary", f"nonpositive-{key}", summary.get(key)))
        except Exception as exc:
            findings.append(issue("summary", f"bad-decimal-{key}", exc))
    finding = str(summary.get("main_finding", "")).lower()
    for text in ("strictly log-concave", "upper", "every real lambda", "remains open", "not cone entry"):
        if text not in finding:
            findings.append(issue("summary", "weak-main-finding", text))
    return findings


def validate_diagnostics(artifact: dict) -> list[CertificateIssue]:
    diagnostics = artifact.get("diagnostics", {})
    compact = diagnostics.get("compact", {})
    ray = diagnostics.get("ray", {})
    findings = []
    flags = {
        "source_full_kernel_evenness_certified": diagnostics.get(
            "source_full_kernel_evenness_certified"
        ),
        "source_order42_residual_zero_certified": diagnostics.get(
            "source_order42_residual_zero_certified"
        ),
        "compact_log_concavity_certified": compact.get(
            "compact_log_concavity_certified"
        ),
        "ray_log_concavity_certified": ray.get("ray_log_concavity_certified"),
        "global_kernel_sqrt_log_concavity_certified": diagnostics.get(
            "global_kernel_sqrt_log_concavity_certified"
        ),
        "all_k_upper_wall_certified": diagnostics.get("all_k_upper_wall_certified"),
    }
    for name, value in flags.items():
        if value is not True:
            findings.append(issue("diagnostics", f"false-{name}", value))
    if compact.get("subinterval_count") != 200:
        findings.append(issue("compact", "bad-subinterval-count", compact.get("subinterval_count")))
    if compact.get("positive_Q_subintervals") != 200:
        findings.append(issue("compact", "bad-positive-Q-count", compact.get("positive_Q_subintervals")))
    if compact.get("positive_g_subintervals") != 200:
        findings.append(issue("compact", "bad-positive-g-count", compact.get("positive_g_subintervals")))
    if compact.get("coefficient_ball_count") != 21:
        findings.append(issue("compact", "bad-coefficient-count", compact.get("coefficient_ball_count")))
    positive_fields = (
        (compact, "minimum_g_lower"),
        (compact, "minimum_Q_lower"),
        (ray, "n1_D_lower"),
        (ray, "n1_third_derivative_margin_lower"),
        (ray, "endpoint_decay_margin_lower"),
        (ray, "full_D_lower"),
    )
    for section, key in positive_fields:
        try:
            if Decimal(str(section.get(key))) <= 0:
                findings.append(issue("diagnostics", f"nonpositive-{key}", section.get(key)))
        except Exception as exc:
            findings.append(issue("diagnostics", f"bad-decimal-{key}", exc))
    try:
        perturbation = Decimal(str(ray.get("perturbation_upper")))
        n1_margin = Decimal(str(ray.get("n1_D_lower")))
        if not (Decimal(0) < perturbation < n1_margin):
            findings.append(issue("ray", "bad-perturbation-order", f"{perturbation}, {n1_margin}"))
    except Exception as exc:
        findings.append(issue("ray", "bad-perturbation-decimal", exc))
    if diagnostics.get("remaining_open_cone_clause") != "x_(k+1)(lambda)>=x_k(lambda)":
        findings.append(
            issue(
                "diagnostics",
                "bad-open-clause",
                diagnostics.get("remaining_open_cone_clause"),
            )
        )
    return findings


def validate_rows(artifact: dict) -> list[CertificateIssue]:
    rows = artifact.get("rows", [])
    findings = []
    if not isinstance(rows, list):
        return [issue("rows", "not-list", type(rows).__name__)]
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        findings.append(issue("rows", "bad-row-ids", sorted(str(item) for item in ids)))
    for row in rows:
        if not isinstance(row, dict):
            findings.append(issue("rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id"))
        if row.get("role") not in ALLOWED_ROLES:
            findings.append(issue(row_id, "bad-role", row.get("role")))
        if row.get("readiness") not in {"proved", "available_exact", "interval_validated", "open"}:
            findings.append(issue(row_id, "bad-readiness", row.get("readiness")))
        boundary = str(row.get("proof_boundary", ""))
        if not boundary:
            findings.append(issue(row_id, "missing-proof-boundary", ""))
        for ref in row.get("source_artifacts", []):
            findings.extend(validate_ref(row_id, ref))
    open_row = next((row for row in rows if row.get("id") == "kmuw_08_remaining_monotone_wall"), {})
    if open_row.get("readiness") != "open" or not open_row.get("gap"):
        findings.append(issue("rows", "missing-open-gate", open_row))
    return findings


def validate_recomputed(artifact: dict, evenness_path: Path) -> list[CertificateIssue]:
    try:
        recomputed = build_artifact(evenness_path)
    except Exception as exc:
        return [issue("recompute", "exception", exc)]
    findings = []
    for key in ("summary", "diagnostics", "rows"):
        if artifact.get(key) != recomputed.get(key):
            findings.append(issue("recompute", f"mismatch-{key}", "stored artifact differs"))
    return findings


def validate_note(path: Path) -> list[CertificateIssue]:
    if not path.exists():
        return [issue("note", "missing-note", path)]
    text = path.read_text(encoding="utf-8")
    findings = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            findings.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in ("proves rh", "proves lambda <= 0", "clay-ready proof"):
        if forbidden in lowered:
            findings.append(issue("note", "overclaim", forbidden))
    return findings


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--evenness-json", type=Path, default=DEFAULT_EVENNESS_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    findings = []
    try:
        artifact = load_json(args.artifact)
    except Exception as exc:
        artifact = {}
        findings.append(issue("artifact", "load-failed", exc))
    if artifact:
        findings.extend(validate_top_level(artifact))
        findings.extend(validate_summary(artifact))
        findings.extend(validate_diagnostics(artifact))
        findings.extend(validate_rows(artifact))
        findings.extend(validate_recomputed(artifact, args.evenness_json))
    findings.extend(validate_symbolic_identities())
    findings.extend(validate_note(args.note))
    ok = not findings
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in findings]}, indent=2, sort_keys=True))
    else:
        for finding in findings:
            print(f"ISSUE {finding.section} [{finding.issue}] {finding.detail}")
        print(
            "validated Jensen-window PF kernel Mellin upper-wall certificate: "
            f"8 rows, {len(findings)} issues, 200 positive compact intervals, "
            "1 positive analytic ray, 1 remaining open cone clause, 0 ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
