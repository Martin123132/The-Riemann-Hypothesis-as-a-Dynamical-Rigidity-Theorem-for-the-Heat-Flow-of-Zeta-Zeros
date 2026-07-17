#!/usr/bin/env python3
"""Validate the endpoint theta cell-renormalization gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys

import mpmath as mp
import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_newman_theta_cell_renormalization_gate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_IDS = {
    "ntcrg_01_continuous_index",
    "ntcrg_02_continuum_annihilation",
    "ntcrg_03_cell_kernel_sum",
    "ntcrg_04_weighted_l1",
    "ntcrg_05_uniform_transform",
    "ntcrg_06_euler_zeta",
    "ntcrg_07_origin_sign",
    "ntcrg_08_coupled_laguerre",
    "ntcrg_09_positive_time_obstruction",
    "ntcrg_10_live_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Newman Theta Cell-Renormalization Gate",
    "integral_0^infinity phi_a(u)da=0",
    "r_n(u)=g_n(u)-integral_(n-1)^n g_a(u)da",
    "sum_(n>=1)integral_R |u|^k*|r_n(u)|du<infinity",
    "sum_(n>=1)e_n(s)=lim_(N->infinity)",
    "J_n(0)>0 for every n>=1",
    "L[H_0](x)=sum_(m,n>=1)M_(m,n)(x)",
    "r_n(u)=-(pi/2)*(3*n-1)*exp(-5u)+O_n(exp(-9u))",
    "the endpoint matrix cannot be deformed block by block to `t>0`",
)

REQUIRED_CORE_STRINGS = (
    "### Lemma 11.24: Theta Cell Renormalization",
    "r_n(u)=g_n(u)-integral_(n-1)^n g_a(u)da",
    "J_n(0)>0",
    "r_n(u)=-(pi/2)*(3n-1)exp(-5u)+O_n(exp(-9u))",
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
    a, u, n = sp.symbols("a u n", positive=True)
    q = sp.pi * sp.exp(4 * u)
    moment_2 = sp.sqrt(sp.pi) / (4 * q ** sp.Rational(3, 2))
    moment_4 = 3 * sp.sqrt(sp.pi) / (8 * q ** sp.Rational(5, 2))
    cancellation = sp.simplify(
        2 * sp.pi**2 * sp.exp(9 * u) * moment_4
        - 3 * sp.pi * sp.exp(5 * u) * moment_2
    )
    if cancellation != 0:
        findings.append(issue("symbolic", "bad-continuum-cancellation", cancellation))

    cell_gap = sp.expand(n**2 - sp.integrate(a**2, (a, n - 1, n)))
    if cell_gap != n - sp.Rational(1, 3):
        findings.append(issue("symbolic", "bad-tail-cell-gap", cell_gap))

    s = sp.symbols("s")
    antiderivative = a ** (1 - s) / (1 - s)
    antiderivative_error = sp.simplify(sp.diff(antiderivative, a) - a ** (-s))
    if antiderivative_error != 0:
        findings.append(issue("symbolic", "bad-cell-power", antiderivative_error))

    mp.mp.dps = 60
    qhat_zero = -mp.gamma(mp.mpf(1) / 4) / (32 * mp.pi ** (mp.mpf(1) / 4))
    for n_value in (1, 2, 10, 100):
        bracket = n_value ** (-mp.mpf(1) / 2) - 2 * (
            mp.sqrt(n_value) - mp.sqrt(n_value - 1)
        )
        if not (bracket < 0 and qhat_zero * bracket / 2 > 0):
            findings.append(issue("numeric", "bad-origin-sign", (n_value, bracket)))

    s_value = mp.mpf("0.5") + mp.mpf("0.7") * 1j
    errors = []
    for count in (200, 800, 3200):
        value = mp.fsum(k ** (-s_value) for k in range(1, count + 1))
        value -= count ** (1 - s_value) / (1 - s_value)
        errors.append(abs(value - mp.zeta(s_value)))
    if not (errors[2] < errors[1] < errors[0]):
        findings.append(issue("numeric", "euler-limit-not-converging", errors))
    return findings


def validate_artifact(artifact: dict) -> list[GateIssue]:
    findings: list[GateIssue] = []
    if artifact.get("kind") != "jensen_window_pf_newman_theta_cell_renormalization_gate":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("date") != "2026-07-11":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    if artifact.get("status") != "exact endpoint cell renormalization with positive-time obstruction":
        findings.append(issue("artifact", "bad-status", artifact.get("status")))

    boundary = str(artifact.get("proof_boundary", "")).lower()
    for phrase in (
        "proves",
        "does not prove",
        "t=0",
        "positive newman",
        "matrix positivity",
        "wiener density",
        "rh",
        "lambda<=0",
    ):
        if phrase not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", phrase))

    for ref in artifact.get("sources", []):
        if not (REPO_ROOT / ref).exists():
            findings.append(issue("artifact", "missing-source", ref))

    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows}
    if len(rows) != 10:
        findings.append(issue("rows", "bad-count", len(rows)))
    if ids != REQUIRED_IDS:
        findings.append(issue("rows", "bad-ids", sorted(str(value) for value in ids)))
    expected_roles = {
        "exact_identity": 5,
        "exact_theorem": 3,
        "nonpromotion_gate": 1,
        "open_handoff": 1,
    }
    for role, count in expected_roles.items():
        actual = sum(row.get("role") == role for row in rows)
        if actual != count:
            findings.append(issue("rows", f"bad-{role}-count", actual))

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
            "validated Jensen-window PF Newman theta cell-renormalization gate: "
            f"10 rows, {len(findings)} issues, 4 exact kernel/transform identities, "
            "3 exact convergence/sign theorems, 1 coupled Laguerre identity, "
            "1 positive-time obstruction, 1 modular handoff"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
