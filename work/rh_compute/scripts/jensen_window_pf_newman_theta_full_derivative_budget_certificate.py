#!/usr/bin/env python3
"""Compose compact and outer bounds into full m=9 derivative budgets."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path
import sys

try:
    import psutil

    psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
except Exception:
    pass


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR_DIR = REPO_ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR_DIR) not in sys.path:
    sys.path.insert(0, str(VENDOR_DIR))

import flint
from flint import arb


STEM = "jensen_window_pf_newman_theta_full_derivative_budget_certificate"
DEFAULT_OUT = (
    REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
COMPACT_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_stable_remainder_arb_quadratic_matrix.json"
)
OUTER_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_stable_remainder_outer_tail_gate.json"
)
CONTRACT_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_adaptive_modular_c1_remainder_contract.json"
)
DATE = "2026-07-24"
PRECISION_BITS = 256
N_VALUES = tuple(range(4, 11))
X_CALIBRATION = 38


@dataclass(frozen=True)
class CertificateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def serialize(value: arb) -> dict:
    if not value.is_finite() or value.lower() <= 0:
        raise RuntimeError(f"invalid positive budget: {value}")
    log10_value = value.log() / arb(10).log()
    return {
        "enclosure": value.str(65, more=True),
        "upper": value.upper().str(65),
        "log10_enclosure": log10_value.str(50, more=True),
        "relative_accuracy_bits": value.rel_accuracy_bits(),
    }


def source_audit() -> tuple[dict, dict, dict, dict]:
    compact = json.loads(COMPACT_SOURCE.read_text(encoding="utf-8"))
    outer = json.loads(OUTER_SOURCE.read_text(encoding="utf-8"))
    contract = json.loads(CONTRACT_SOURCE.read_text(encoding="utf-8"))
    cache = compact.get("cache", {})
    if not cache.get("complete") or cache.get("records") != 223:
        raise RuntimeError("stable compact matrix is not complete")
    if len(compact.get("compact_rows", [])) != len(N_VALUES):
        raise RuntimeError("stable compact retained rows are incomplete")
    if len(outer.get("outer_aggregates", [])) != len(N_VALUES):
        raise RuntimeError("outer retained rows are incomplete")
    expected_contract = (
        "d_(N,j,m)(T)=(1/2)*sup_(0<=t<=T)"
        "||partial_u^m[(iu)^j*exp(tu^2)r_N(u)]||_1"
    )
    exact = contract.get("exact", {})
    stored_contract = exact.get("tail_budgets", {}).get("derivative")
    if stored_contract != expected_contract:
        raise RuntimeError(
            f"direct-C1 derivative contract drifted: {stored_contract!r}"
        )
    hashes = {
        "compact_source_sha256": file_hash(COMPACT_SOURCE),
        "outer_source_sha256": file_hash(OUTER_SOURCE),
        "contract_source_sha256": file_hash(CONTRACT_SOURCE),
        "contract_formula": expected_contract,
    }
    return compact, outer, contract, hashes


def build_budgets(compact: dict, outer: dict) -> list[dict]:
    compact_by_n = {
        int(row["N"]): row for row in compact["compact_rows"]
    }
    outer_by_n = {
        int(row["N"]): row for row in outer["outer_aggregates"]
    }
    x = arb(X_CALIBRATION)
    rows: list[dict] = []
    previous_d0: arb | None = None
    previous_d1: arb | None = None
    for retained in N_VALUES:
        compact_row = compact_by_n[retained]
        outer_row = outer_by_n[retained]
        compact_d0 = arb(compact_row["compact_full_d0"]["upper"])
        compact_d1 = arb(compact_row["compact_full_d1"]["upper"])
        outer_d0 = arb(outer_row["outer_d0"]["upper"])
        outer_d1 = arb(outer_row["outer_d1"]["upper"])
        full_d0 = compact_d0 + outer_d0
        full_d1 = compact_d1 + outer_d1
        if (
            previous_d0 is not None
            and full_d0.upper() >= previous_d0.lower()
        ):
            raise RuntimeError(f"d0 monotonic separation failed at N={retained}")
        if (
            previous_d1 is not None
            and full_d1.upper() >= previous_d1.lower()
        ):
            raise RuntimeError(f"d1 monotonic separation failed at N={retained}")
        previous_d0 = full_d0
        previous_d1 = full_d1

        j_error = 16 * full_d0 / x**5
        j_prime_error = 64 * full_d0 / x**6 + 16 * full_d1 / x**5
        rows.append(
            {
                "N": retained,
                "compact_d0": serialize(compact_d0),
                "outer_d0": serialize(outer_d0),
                "full_d0_m9_t1_5": serialize(full_d0),
                "compact_d1": serialize(compact_d1),
                "outer_d1": serialize(outer_d1),
                "full_d1_m9_t1_5": serialize(full_d1),
                "x38_J_error_upper": serialize(j_error),
                "x38_J_prime_error_upper": serialize(j_prime_error),
                "uniform_frequency_statement": (
                    "|J-J_N|<=16*d0/x^5; "
                    "|J'-J_N'|<=64*d0/x^6+16*d1/x^5, x>=38"
                ),
            }
        )
    return rows


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    compact, outer, _contract, audit = source_audit()
    budgets = build_budgets(compact, outer)
    rows = [
        CertificateRow(
            id="ntfdbc_01_exact_domain_split",
            role="exact_identity",
            readiness="proved",
            claim=(
                "Each weighted derivative norm splits at u=11/5 into "
                "the certified compact and outer pieces."
            ),
            formula="d_(N,j,9)<=d_compact_(N,j)+d_outer_(N,j)",
            proof_boundary=(
                "This is an upper bound because both pieces were enclosed "
                "after the same positive derivative envelope."
            ),
        ),
        CertificateRow(
            id="ntfdbc_02_compact_full_arithmetic",
            role="directed_rounding_input",
            readiness="proved",
            claim=(
                "The 223-row stable compact matrix covers the complete "
                "arithmetic remainder on 0<=u<=11/5."
            ),
            formula=(
                "int_0^(11/5)W|D^q r_N|"
                "<=sqrt(MQ_f)+B_q*M"
            ),
            proof_boundary="Certified only for 4<=N<=10.",
        ),
        CertificateRow(
            id="ntfdbc_03_outer_full_arithmetic",
            role="directed_rounding_input",
            readiness="proved",
            claim=(
                "The explicit phase envelopes cover the complete arithmetic "
                "remainder on u>11/5."
            ),
            formula=(
                "outer_entry<=forward_initial/lambda_4"
                "+sum defect_initial/lambda_8"
            ),
            proof_boundary="Certified only for 4<=N<=10.",
        ),
        CertificateRow(
            id="ntfdbc_04_full_m9_budgets",
            role="rigorous_derivative_budget",
            readiness="proved",
            claim=(
                "Full d_(N,0,9)(1/5) and d_(N,1,9)(1/5) upper balls are "
                "available for every N from 4 through 10."
            ),
            formula=(
                "d0_N=compact_d0_N+outer_d0_N; "
                "d1_N=compact_d1_N+outer_d1_N"
            ),
            proof_boundary=(
                "Finite retained range; no extrapolation to adaptive N>10."
            ),
        ),
        CertificateRow(
            id="ntfdbc_05_direct_c1_errors",
            role="exact_reduction",
            readiness="proved",
            claim=(
                "The full budgets give explicit J and J' remainder bars for "
                "every x>=38."
            ),
            formula=(
                "|J-J_N|<=16*d0_N/x^5; "
                "|J'-J_N'|<=64*d0_N/x^6+16*d1_N/x^5"
            ),
            proof_boundary=(
                "These are upper error bars, not retained first-jet lower "
                "separation."
            ),
        ),
        CertificateRow(
            id="ntfdbc_06_finite_monotone_calibration",
            role="certified_diagnostic",
            readiness="proved",
            claim=(
                "The seven certified d0 and d1 upper budgets are strictly "
                "separated downward as N increases."
            ),
            formula="d_(N+1,j,9)<d_(N,j,9), 4<=N<10, j in {0,1}",
            proof_boundary=(
                "A finite checked trend is not an all-N monotonicity theorem."
            ),
        ),
        CertificateRow(
            id="ntfdbc_07_cofinal_separation_handoff",
            role="proof_search_target",
            readiness="not_ready_to_apply",
            claim=(
                "Prove retained J/J' value-or-derivative separation on all "
                "transition cells and extend budgets to the required "
                "unbounded adaptive N range."
            ),
            formula=(
                "|J_N|>16*d0_N/x^5 OR "
                "|J_N'|>64*d0_N/x^6+16*d1_N/x^5"
            ),
            proof_boundary=(
                "This open disjunction is the noncircular Newman handoff."
            ),
        ),
    ]
    return {
        "kind": STEM,
        "date": DATE,
        "status": (
            "complete finite-N full derivative budgets with open cofinal "
            "first-jet separation"
        ),
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "m_order": 9,
            "t_cap_exact": "1/5",
            "u_split_exact": "11/5",
            "n_values": list(N_VALUES),
            "frequency_calibration": X_CALIBRATION,
        },
        "source_audit": audit,
        "budgets": budgets,
        "rows": [asdict(row) for row in rows],
        "proof_boundary": (
            "This certificate completes rigorous d0/d1 derivative-error "
            "budgets only for N=4..10. It does not prove retained first-jet "
            "lower separation, cover unbounded adaptive N or all transition "
            "cells, prove Lambda<=0 or RH, or supply a Clay-prize conclusion."
        ),
        "versions": {
            "python_flint": getattr(flint, "__version__", "unknown"),
        },
    }


def render_note(artifact: dict) -> str:
    lines = [
        "# Newman Theta Full Derivative-Budget Certificate",
        "",
        f"Date: {DATE}",
        "",
        "Status: rigorous full derivative-error budgets for `N=4..10`.",
        "This is not the retained first-jet separation theorem and not a",
        "proof of `Lambda<=0`, RH, or a Clay-prize result.",
        "",
        "## Complete Bounds",
        "",
        "For `m=9`, `T=1/5`, and every `x>=38`,",
        "",
        "```text",
        "|J-J_N|<=16*d0_N/x^5",
        "|J'-J_N'|<=64*d0_N/x^6+16*d1_N/x^5",
        "```",
        "",
        "| N | full d0 upper | full d1 upper | J error at 38 | J' error at 38 |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in artifact["budgets"]:
        lines.append(
            f"| {row['N']} | "
            f"`{row['full_d0_m9_t1_5']['upper']}` | "
            f"`{row['full_d1_m9_t1_5']['upper']}` | "
            f"`{row['x38_J_error_upper']['upper']}` | "
            f"`{row['x38_J_prime_error_upper']['upper']}` |"
        )
    lines.extend(
        [
            "",
            "Each full bound is compact plus outer, and each piece covers",
            "the complete arithmetic remainder rather than a finite `n` sum.",
            "",
            "## Open Handoff",
            "",
            "The remaining proof obligation is a retained `J/J'` lower",
            "separation on every transition cell, together with a cofinal",
            "extension beyond the finite `N<=10` calibration.",
            "",
        ]
    )
    return "\n".join(lines)


def write_artifact(artifact: dict, out: Path, note: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    note.write_text(render_note(artifact), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    write_artifact(artifact, args.out, args.note)
    print(
        "wrote Newman theta full derivative-budget certificate: "
        f"{len(artifact['budgets'])} retained counts, "
        f"{len(artifact['rows'])} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
