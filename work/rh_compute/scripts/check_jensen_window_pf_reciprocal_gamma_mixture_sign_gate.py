#!/usr/bin/env python3
"""Validate the reciprocal-gamma mixture sign gate."""

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

from jensen_window_pf_reciprocal_gamma_mixture_sign_gate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
    gamma_value,
)


REQUIRED_IDS = {
    "rgmsg_01_coefficient_split",
    "rgmsg_02_karlin_theorem",
    "rgmsg_03_fixed_scale",
    "rgmsg_04_contiguous_formula",
    "rgmsg_05_rowwise_integral",
    "rgmsg_06_integrand_sign_change",
    "rgmsg_07_mixture_countermodel",
    "rgmsg_08_tilted_cv",
    "rgmsg_09_xi_order_two",
    "rgmsg_10_higher_order_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Reciprocal-Gamma Mixture Sign Gate",
    "Karlin's Lemma 9.2",
    "gamma_k=k!/(2k)!=sqrt(pi)/(4^k*Gamma(k+1/2))",
    "D_sym(t_0,t_1)=(t_0^2-6*t_0*t_1+t_1^2)/24",
    "rho=10*delta_(1/100)+delta_1+delta_2",
    "A_0*A_2-A_1^2=5197/2000>0",
    "CV_n^2<=2/(2*n+1)",
    "higher compound",
    "https://doi.org/10.1090/S0002-9947-1964-0168010-2",
)

REQUIRED_CORE_STRINGS = (
    "### Lemma 11.25: Reciprocal-Gamma Mixture Sign Gate",
    "gamma_k=k!/(2k)!",
    "A_0*A_2-A_1^2=5197/2000>0",
    "CV_n^2<=2/(2n+1)",
)


@dataclass(frozen=True)
class GateIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: object) -> GateIssue:
    return GateIssue(section=section, issue=name, detail=str(detail))


def symbolic_issues() -> list[GateIssue]:
    findings: list[GateIssue] = []
    for size in range(1, 7):
        matrix = sp.Matrix(
            [[gamma_value(i + j) for j in range(size)] for i in range(size)]
        )
        formula = (
            (-1) ** (size * (size - 1) // 2)
            * sp.pi ** sp.Rational(size, 2)
            * sp.Rational(1, 4) ** (size * (size - 1))
            * sp.prod(sp.factorial(r) for r in range(1, size))
            / sp.prod(
                sp.gamma(sp.Rational(1, 2) + size - 1 + i)
                for i in range(size)
            )
        )
        if sp.simplify(matrix.det() - formula) != 0:
            findings.append(issue("symbolic", f"bad-contiguous-size-{size}", matrix.det()))

    t0, t1 = sp.symbols("t0 t1", positive=True)
    sym_integrand = (
        gamma_value(0) * gamma_value(2) * (t0**2 + t1**2) / 2
        - gamma_value(1) ** 2 * t0 * t1
    )
    target = (t0**2 - 6 * t0 * t1 + t1**2) / 24
    if sp.expand(sym_integrand - target) != 0:
        findings.append(issue("symbolic", "bad-order-two-integrand", sym_integrand))
    if sp.simplify(target.subs({t0: 1, t1: 1}) + sp.Rational(1, 6)) != 0:
        findings.append(issue("symbolic", "bad-negative-witness", target.subs({t0: 1, t1: 1})))
    if sp.simplify(target.subs({t0: 1, t1: 10}) - sp.Rational(41, 24)) != 0:
        findings.append(issue("symbolic", "bad-positive-witness", target.subs({t0: 1, t1: 10})))

    weights = (10, 1, 1)
    scales = (sp.Rational(1, 100), sp.Integer(1), sp.Integer(2))
    moments = [sum(w * t**order for w, t in zip(weights, scales)) for order in range(3)]
    coefficients = [gamma_value(order) * moments[order] for order in range(3)]
    mixture_minor = sp.factor(coefficients[0] * coefficients[2] - coefficients[1] ** 2)
    if mixture_minor != sp.Rational(5197, 2000):
        findings.append(issue("symbolic", "bad-mixture-countermodel", mixture_minor))

    n = sp.symbols("n", integer=True, nonnegative=True)
    theta = (2 * n + 1) / (2 * n + 3)
    if sp.simplify(1 / theta - 1 - 2 / (2 * n + 1)) != 0:
        findings.append(issue("symbolic", "bad-cv-threshold", 1 / theta - 1))
    return findings


def validate_artifact(artifact: dict) -> list[GateIssue]:
    findings: list[GateIssue] = []
    if artifact.get("kind") != "jensen_window_pf_reciprocal_gamma_mixture_sign_gate":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("date") != "2026-07-12":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    if artifact.get("status") != "exact reciprocal-gamma sign regularity with positive-mixture obstruction":
        findings.append(issue("artifact", "bad-status", artifact.get("status")))

    boundary = str(artifact.get("proof_boundary", "")).lower()
    for phrase in (
        "published",
        "proves",
        "does not prove",
        "mixture countermodel",
        "higher-order",
        "signed-hankel/jensen",
        "pf-infinity",
        "rh",
        "lambda<=0",
    ):
        if phrase not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", phrase))

    for ref in artifact.get("sources", []):
        if ref.startswith("https://"):
            continue
        if not (REPO_ROOT / ref).exists():
            findings.append(issue("artifact", "missing-source", ref))

    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows}
    if len(rows) != 10:
        findings.append(issue("rows", "bad-count", len(rows)))
    if ids != REQUIRED_IDS:
        findings.append(issue("rows", "bad-ids", sorted(str(value) for value in ids)))
    role_counts = {
        role: sum(row.get("role") == role for row in rows)
        for role in {
            "exact_identity",
            "published_theorem",
            "exact_theorem",
            "exact_countermodel",
            "exact_equivalence",
            "exact_composition",
            "open_handoff",
        }
    }
    expected_counts = {
        "exact_identity": 3,
        "published_theorem": 1,
        "exact_theorem": 1,
        "exact_countermodel": 2,
        "exact_equivalence": 1,
        "exact_composition": 1,
        "open_handoff": 1,
    }
    if role_counts != expected_counts:
        findings.append(issue("rows", "bad-role-counts", role_counts))

    source = artifact.get("exact", {}).get("karlin_sign_regularity", {}).get("source")
    if source != "https://doi.org/10.1090/S0002-9947-1964-0168010-2":
        findings.append(issue("source", "bad-primary-source", source))
    if artifact != build_artifact():
        findings.append(issue("recompute", "mismatch", "stored artifact differs from exact rebuild"))
    return findings


def required_string_issues(path: Path, strings: tuple[str, ...], section: str) -> list[GateIssue]:
    text = path.read_text(encoding="utf-8")
    return [issue(section, "missing-string", value) for value in strings if value not in text]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    findings = symbolic_issues()
    findings.extend(validate_artifact(artifact))
    findings.extend(required_string_issues(args.note, REQUIRED_NOTE_STRINGS, "note"))
    findings.extend(
        required_string_issues(
            REPO_ROOT / "outputs/formal_core.md", REQUIRED_CORE_STRINGS, "formal-core"
        )
    )
    ok = not findings
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in findings]}, indent=2, sort_keys=True))
    else:
        for finding in findings:
            print(f"ISSUE {finding.section} [{finding.issue}] {finding.detail}")
        print(
            "validated Jensen-window PF reciprocal-gamma mixture sign gate: "
            f"10 rows, {len(findings)} issues, 3 exact kernel identities, "
            "1 published all-order theorem, 1 fixed-scale theorem, "
            "2 exact mixture countermodels, 1 tilted-variance equivalence, "
            "1 Xi order-two composition, 1 higher-order handoff"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
