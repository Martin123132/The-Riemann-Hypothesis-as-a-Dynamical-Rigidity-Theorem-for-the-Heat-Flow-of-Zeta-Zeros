#!/usr/bin/env python3
"""Validate the theta/Bessel higher-shift regularization gate."""

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

from jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "ntbhsr_01_theta_symmetrization",
    "ntbhsr_02_higher_shift_expansion",
    "ntbhsr_03_coefficient_formula",
    "ntbhsr_04_coefficient_sign",
    "ntbhsr_05_fixed_block_bessel_transform",
    "ntbhsr_06_zero_frequency_divergence",
    "ntbhsr_07_fubini_boundary",
    "ntbhsr_08_continuation_guard",
    "ntbhsr_09_coupled_matrix_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Newman Theta/Bessel Higher-Shift Regularization Gate",
    "Status: exact higher-shift expansion with spectral summation obstruction.",
    "c_(n,m)=pi^m*n^(2m)*(2*pi^2*n^4-3*m)/m!",
    "g_n(u)=exp(-2*pi*n^2*cosh(4u))*sum_(m>=0)",
    "c_(1,m)>0 for 0<=m<=6 and c_(1,m)<0 for m>=7",
    "sum_(n>=1)I_n(0)=-infinity",
    "Qhat(x)*zeta((1+i*x)/2)=xi((1+i*x)/2)/4",
    "M_(i,j)=F_i'*F_j'-(F_i*F_j''+F_j*F_i'')/2",
    "Naive ordinary higher-shift summation is now closed.",
    "outputs/jensen_window_pf_newman_theta_summand_spectral_square_gate.md",
)

REQUIRED_CORE_STRINGS = (
    "### Lemma 11.23: Theta/Bessel Higher-Shift Regularization Gate",
    "c_(n,m)=pi^m*n^(2m)*(2*pi^2*n^4-3m)/m!.",
    "sum_(n>=1)I_n(0)=-infinity",
    "M_(i,j)=F_i'*F_j'-(F_i*F_j''+F_j*F_i'')/2.",
    "it proves neither strict Laguerre positivity, Wiener",
)


@dataclass(frozen=True)
class GateIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: object) -> GateIssue:
    return GateIssue(section=section, issue=name, detail=str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(ref: object) -> list[GateIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue("artifact", "bad-ref", repr(ref))]
    if ref.startswith("https://"):
        return []
    if not (REPO_ROOT / ref).exists():
        return [issue("artifact", "missing-ref", ref)]
    return []


def symbolic_issues() -> list[GateIssue]:
    findings: list[GateIssue] = []
    n, m = sp.symbols("n m", integer=True, positive=True)
    p = sp.pi * n**2
    combined = 2 * sp.pi**2 * n**4 * p**m - 3 * sp.pi * n**2 * m * p ** (m - 1)
    target = p**m * (2 * sp.pi**2 * n**4 - 3 * m)
    if sp.simplify(combined - target) != 0:
        findings.append(issue("symbolic", "bad-coefficient", combined - target))

    if not (sp.N(2 * sp.pi**2 - 18, 80) > 0 and sp.N(2 * sp.pi**2 - 21, 80) < 0):
        findings.append(issue("symbolic", "bad-first-sign-threshold", sp.N(2 * sp.pi**2 / 3, 80)))

    for n_value, expected_last in ((1, 6), (2, 105), (3, 532)):
        threshold = sp.N(2 * sp.pi**2 * n_value**4 / 3, 80)
        if int(sp.floor(threshold)) != expected_last:
            findings.append(issue("symbolic", f"bad-threshold-n{n_value}", threshold))

    f1, f2, f1p, f2p, f1pp, f2pp, c1, c2 = sp.symbols(
        "f1 f2 f1p f2p f1pp f2pp c1 c2", real=True
    )
    direct = (c1 * f1p + c2 * f2p) ** 2 - (c1 * f1 + c2 * f2) * (
        c1 * f1pp + c2 * f2pp
    )
    matrix = (
        c1**2 * (f1p**2 - f1 * f1pp)
        + 2 * c1 * c2 * (f1p * f2p - (f1 * f2pp + f2 * f1pp) / 2)
        + c2**2 * (f2p**2 - f2 * f2pp)
    )
    if sp.expand(direct - matrix) != 0:
        findings.append(issue("symbolic", "bad-matrix-identity", direct - matrix))

    qhat_zero = -sp.gamma(sp.Rational(1, 4)) / (32 * sp.pi ** sp.Rational(1, 4))
    if not sp.N(qhat_zero, 80) < 0:
        findings.append(issue("symbolic", "bad-qhat-zero-sign", qhat_zero))

    mp.mp.dps = 80
    for n_value, u_value in (
        (1, mp.mpf("0")),
        (1, mp.mpf("0.35")),
        (2, mp.mpf("0.15")),
    ):
        def phi(value: mp.mpf) -> mp.mpf:
            return (
                2 * mp.pi**2 * n_value**4 * mp.exp(9 * value)
                - 3 * mp.pi * n_value**2 * mp.exp(5 * value)
            ) * mp.exp(-mp.pi * n_value**2 * mp.exp(4 * value))

        direct_value = (phi(u_value) + phi(-u_value)) / 2
        series_value = mp.fsum(
            mp.pi**m_value
            * n_value ** (2 * m_value)
            * (2 * mp.pi**2 * n_value**4 - 3 * m_value)
            / mp.factorial(m_value)
            * mp.cosh((9 - 4 * m_value) * u_value)
            for m_value in range(220)
        ) * mp.exp(-2 * mp.pi * n_value**2 * mp.cosh(4 * u_value))
        relative_error = abs(series_value - direct_value) / max(
            abs(direct_value), mp.mpf("1e-70")
        )
        if relative_error >= mp.mpf("1e-60"):
            findings.append(
                issue(
                    "numeric",
                    "higher-shift-normalization",
                    f"n={n_value}, u={u_value}, relative_error={relative_error}",
                )
            )
    return findings


def validate_artifact(artifact: dict) -> list[GateIssue]:
    findings: list[GateIssue] = []
    if artifact.get("kind") != "jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("date") != "2026-07-11":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    if artifact.get("status") != "exact higher-shift expansion with spectral summation obstruction":
        findings.append(issue("artifact", "bad-status", artifact.get("status")))

    boundary = str(artifact.get("proof_boundary", "")).lower()
    for phrase in (
        "proves",
        "does not reject",
        "does not prove",
        "zero frequency",
        "strict laguerre",
        "wiener density",
        "rh",
        "lambda<=0",
    ):
        if phrase not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", phrase))

    for ref in artifact.get("sources", []):
        findings.extend(validate_ref(ref))

    rows = artifact.get("rows", [])
    row_ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if len(rows) != 9:
        findings.append(issue("rows", "bad-count", len(rows)))
    if row_ids != REQUIRED_ROW_IDS:
        findings.append(issue("rows", "bad-ids", sorted(str(value) for value in row_ids)))
    if sum(row.get("role") == "exact_identity" for row in rows) != 3:
        findings.append(issue("rows", "bad-exact-identity-count", rows))
    if sum(row.get("role") == "exact_theorem" for row in rows) != 2:
        findings.append(issue("rows", "bad-exact-theorem-count", rows))
    if sum(row.get("role") == "nonpromotion_gate" for row in rows) != 3:
        findings.append(issue("rows", "bad-nonpromotion-count", rows))
    if sum(row.get("role") == "open_handoff" for row in rows) != 1:
        findings.append(issue("rows", "bad-handoff-count", rows))

    exact = artifact.get("exact", {})
    if exact.get("spectral_summation_obstruction", {}).get("divergence") != "sum_(n>=1)I_n(0)=-infinity":
        findings.append(issue("exact", "missing-divergence", exact.get("spectral_summation_obstruction")))
    sign_rows = exact.get("coefficient_sign", {}).get("rows", [])
    expected_thresholds = [(1, 6, 7), (2, 105, 106), (3, 532, 533)]
    actual_thresholds = [
        (row.get("n"), row.get("last_positive_m"), row.get("first_negative_m"))
        for row in sign_rows
    ]
    if actual_thresholds != expected_thresholds:
        findings.append(issue("exact", "bad-sign-rows", actual_thresholds))

    if artifact != build_artifact():
        findings.append(issue("recompute", "mismatch", "stored artifact differs from exact rebuild"))
    return findings


def validate_note(path: Path) -> list[GateIssue]:
    findings: list[GateIssue] = []
    text = path.read_text(encoding="utf-8")
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            findings.append(issue("note", "missing-string", required))
    return findings


def validate_formal_core() -> list[GateIssue]:
    findings: list[GateIssue] = []
    text = (REPO_ROOT / "outputs/formal_core.md").read_text(encoding="utf-8")
    for required in REQUIRED_CORE_STRINGS:
        if required not in text:
            findings.append(issue("formal-core", "missing-string", required))
    return findings


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    findings = symbolic_issues()
    findings.extend(validate_artifact(load_json(args.artifact)))
    findings.extend(validate_note(args.note))
    findings.extend(validate_formal_core())
    ok = not findings
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in findings]}, indent=2, sort_keys=True))
    else:
        for finding in findings:
            print(f"ISSUE {finding.section} [{finding.issue}] {finding.detail}")
        print(
            "validated Jensen-window PF Newman theta/Bessel higher-shift regularization gate: "
            f"9 rows, {len(findings)} issues, 3 exact expansion identities, "
            "1 coefficient sign theorem, 1 fixed-block Bessel theorem, "
            "3 spectral non-promotion gates, 1 coupled modular handoff"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
