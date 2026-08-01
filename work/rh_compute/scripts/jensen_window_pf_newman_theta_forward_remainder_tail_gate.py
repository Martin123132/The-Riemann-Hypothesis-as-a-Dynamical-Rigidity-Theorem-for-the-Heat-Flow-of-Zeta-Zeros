#!/usr/bin/env python3
"""Build the exact stable forward-remainder arithmetic-tail gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_theta_forward_remainder_tail_gate"
DEFAULT_OUT = REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
MODULAR_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_modular_blend_gate.json"
)
DATE = "2026-07-24"
MAX_ORDER = 9
FORWARD_CAP = 12
FIRST_OMITTED = FORWARD_CAP + 1


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def polynomial_rows() -> list[dict]:
    x = sp.symbols("X", nonnegative=True)
    current = 2 * x - 3
    rows: list[dict] = []
    for order in range(MAX_ORDER + 1):
        polynomial = sp.Poly(sp.expand(current), x)
        if polynomial.degree() != order + 1:
            raise RuntimeError(f"unexpected degree at order {order}")
        coefficient_norm = sum(abs(value) for value in polynomial.all_coeffs())
        rows.append(
            {
                "q": order,
                "P_q": str(polynomial.as_expr()),
                "degree": polynomial.degree(),
                "A_q": int(coefficient_norm),
            }
        )
        expected_next = sp.expand(
            (5 - 4 * x) * polynomial.as_expr()
            + 4 * x * sp.diff(polynomial.as_expr(), x)
        )
        current = expected_next
    return rows


def tail_rows(polynomials: list[dict]) -> list[dict]:
    rows: list[dict] = []
    k = sp.Integer(FIRST_OMITTED)
    for polynomial in polynomials:
        q = polynomial["q"]
        a = 4 * q + 9
        d = 2 * q + 4
        monotonic_margin = sp.simplify(4 * sp.pi * k**2 - a)
        if not bool(sp.N(monotonic_margin, 50) > 0):
            raise RuntimeError(f"u-monotonicity failed at q={q}")
        rho = sp.exp(sp.Rational(d, FIRST_OMITTED) - sp.pi * (2 * k + 1))
        if not bool(sp.N(rho, 50) < 1):
            raise RuntimeError(f"geometric ratio failed at q={q}")
        bound = (
            polynomial["A_q"]
            * sp.pi ** (q + 2)
            * k ** (2 * q + 4)
            * sp.exp(-sp.pi * k**2)
            / (1 - rho)
        )
        rows.append(
            {
                "q": q,
                "a_q": a,
                "d_q": d,
                "u_monotonicity_margin": str(monotonic_margin),
                "rho_exact": str(rho),
                "rho_decimal": str(sp.N(rho, 30)),
                "B_q_exact": str(bound),
                "B_q_decimal": str(sp.N(bound, 30)),
                "log10_B_q": str(sp.N(sp.log(bound, 10), 30)),
            }
        )
    return rows


def source_audit() -> dict:
    source = json.loads(MODULAR_SOURCE.read_text(encoding="utf-8"))
    exact = source.get("exact", {})
    identity = exact.get("decaying_tail_enclosure", {}).get("tail_kernel")
    required = (
        "r_N(u)=sum_(n>N)b_n(u)=Phi(u)-sum_(n<=N)b_n(u)"
    )
    if identity != required:
        raise RuntimeError(f"modular tail identity drifted: {identity!r}")
    series = exact.get("theta_summand", {}).get("kernel_series", "")
    if "Phi(u)=sum_(n>=1)phi_n(u)" not in series:
        raise RuntimeError("forward theta series marker missing")
    return {
        "kind": source["kind"],
        "sha256": file_hash(MODULAR_SOURCE),
        "tail_identity": identity,
        "forward_series": series,
    }


def build_artifact() -> dict:
    polynomials = polynomial_rows()
    tails = tail_rows(polynomials)
    stable_remainder = (
        "delta_n(u)=(erfc(3*sinh(4u))/2)"
        "*(phi_n(u)-phi_n(-u)); "
        "r_N^(q)(u)=D^q[sum_(n=N+1)^12 phi_n(u)"
        "+sum_(n=1)^N delta_n(u)]+tau_q(u), "
        "tau_q(u)=sum_(n>=13)D^q phi_n(u), 1<=N<=12"
    )
    rows = [
        GateRow(
            id="ntfrtg_01_stable_modular_remainder",
            role="exact_identity",
            readiness="proved",
            claim=(
                "The blended tail can be evaluated as the rapidly convergent "
                "forward theta series minus the retained blended blocks."
            ),
            formula=(
                "r_N=Phi-sum_(n<=N)b_n="
                "sum_(n=N+1)^12 phi_n+sum_(n<=N)delta_n+"
                "sum_(n>=13)phi_n; "
                "delta_n=phi_n-b_n=(erfc(3*sinh(4u))/2)"
                "*(phi_n(u)-phi_n(-u))"
            ),
            proof_boundary=(
                "This is the stored theta modular partition identity and uses "
                "no zeta-zero information; the displayed finite split applies "
                "for 1<=N<=12."
            ),
        ),
        GateRow(
            id="ntfrtg_02_forward_derivative_recurrence",
            role="exact_identity",
            readiness="proved",
            claim=(
                "Every forward derivative is one Gaussian times an explicit "
                "integer polynomial generated by a first-order recurrence."
            ),
            formula=(
                "D^q phi_n=pi*n^2*exp(5u)*P_q(X)*exp(-X), "
                "X=pi*n^2*exp(4u); P_0=2X-3; "
                "P_(q+1)=(5-4X)P_q+4X P_q'"
            ),
            proof_boundary=(
                "Exact differentiation identity through q=9; it does not "
                "bound the retained finite remainder."
            ),
            diagnostics=polynomials,
        ),
        GateRow(
            id="ntfrtg_03_uniform_forward_tail",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "For n>=13 and u>=0, each forward derivative is maximized "
                "at u=0 after replacing P_q by its coefficient 1-norm."
            ),
            formula=(
                "|D^q phi_n(u)|<=A_q*pi^(q+2)*"
                "n^(2q+4)*exp(-pi*n^2)"
            ),
            proof_boundary=(
                "Uses X>=1 and 4*pi*13^2>4q+9 for 0<=q<=9."
            ),
        ),
        GateRow(
            id="ntfrtg_04_geometric_arithmetic_sum",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "The n>=13 derivative tail has a closed explicit geometric "
                "majorant B_q for every q<=9."
            ),
            formula=(
                "sup_(u>=0)|tau_q(u)|<=B_q="
                "A_q*pi^(q+2)*13^(2q+4)*exp(-169pi)/(1-rho_q)"
            ),
            proof_boundary=(
                "Uses log(1+1/n)<=1/n and "
                "rho_q=exp((2q+4)/13-27pi)<1."
            ),
            diagnostics=tails,
        ),
        GateRow(
            id="ntfrtg_05_compact_quadratic_correction",
            role="exact_reduction",
            readiness="proved",
            claim=(
                "The exact forward arithmetic tail enters each compact "
                "weighted L1 term as an explicit additive B_q M_(p,b)."
            ),
            formula=(
                "int_I W|f+tau_q|<="
                "sqrt(M_(p,b)*Q_f)+B_q*M_(p,b), I=[0,11/5]"
            ),
            proof_boundary=(
                "This closes the n>=13 compact arithmetic tail only; it does "
                "not cover u>11/5 or retained first-jet separation."
            ),
        ),
        GateRow(
            id="ntfrtg_06_stable_matrix_handoff",
            role="proof_search_target",
            readiness="not_ready_to_apply",
            claim=(
                "Build the directed-rounding matrix from the stable finite "
                "remainder, then append an independent outer-u bound."
            ),
            formula=stable_remainder,
            proof_boundary=(
                "No complete d0/d1 tail budget, transition-cell theorem, "
                "Lambda<=0, RH, or Clay-prize conclusion is supplied."
            ),
        ),
    ]
    return {
        "kind": STEM,
        "date": DATE,
        "status": "exact stable remainder and compact arithmetic-tail theorem",
        "source_audit": source_audit(),
        "parameters": {
            "max_derivative_order": MAX_ORDER,
            "forward_cap": FORWARD_CAP,
            "first_omitted": FIRST_OMITTED,
            "retained_range": ["1", "12"],
            "compact_interval": ["0", "11/5"],
        },
        "polynomials": polynomials,
        "tail_bounds": tails,
        "stable_remainder": stable_remainder,
        "compact_correction": (
            "sqrt(M_(p,b)Q_f)+B_q M_(p,b)"
        ),
        "rows": [asdict(row) for row in rows],
        "proof_boundary": (
            "This artifact proves a cancellation-preserving representation "
            "of r_N and explicit uniform bounds for the forward arithmetic "
            "tail n>=13 through derivative order nine. It does not bound the "
            "outer interval u>11/5, certify retained first-jet J/J' lower "
            "separation, "
            "cover transition cells, prove Lambda<=0 or RH, or supply a "
            "Clay-prize proof."
        ),
    }


def render_note(artifact: dict) -> str:
    lines = [
        "# Newman Theta Forward-Remainder Tail Gate",
        "",
        f"Date: {DATE}",
        "",
        "Status: exact stable-remainder and compact arithmetic-tail theorem.",
        "This is not a complete derivative budget and not a proof of",
        "`Lambda<=0`, RH, or a Clay-prize result.",
        "",
        "## Stable Remainder",
        "",
        "The exact modular partition gives",
        "",
        "```text",
        "r_N=Phi-sum_(n<=N)b_n",
        "   =sum_(n>=1)phi_n-sum_(n<=N)b_n.",
        "delta_n=phi_n-b_n",
        "       =(erfc(3*sinh(4u))/2)*(phi_n(u)-phi_n(-u)).",
        artifact["stable_remainder"],
        "```",
        "",
        "This representation preserves modular cancellation while leaving",
        "only the rapidly convergent forward theta tail to bound.",
        "",
        "## Derivative Recurrence",
        "",
        "With `X=pi*n^2*exp(4u)`,",
        "",
        "```text",
        "D^q phi_n=pi*n^2*exp(5u)*P_q(X)*exp(-X)",
        "P_0=2X-3",
        "P_(q+1)=(5-4X)P_q+4X P_q'",
        "```",
        "",
        "For `A_q` equal to the coefficient 1-norm of `P_q`, `n>=13`,",
        "and `u>=0`,",
        "",
        "```text",
        "|D^q phi_n(u)|",
        " <=A_q*pi^(q+2)*n^(2q+4)*exp(-pi*n^2).",
        "```",
        "",
        "## Explicit Tail Constants",
        "",
        "| q | A_q | log10(B_q) |",
        "|---:|---:|---:|",
    ]
    for polynomial, tail in zip(
        artifact["polynomials"], artifact["tail_bounds"], strict=True
    ):
        lines.append(
            f"| {polynomial['q']} | {polynomial['A_q']} | "
            f"`{tail['log10_B_q']}` |"
        )
    lines.extend(
        [
            "",
            "Here `sup_(u>=0)|tau_q(u)|<=B_q` for the omitted forward",
            "series `n>=13`.",
            "",
            "## Compact Matrix Correction",
            "",
            "For every matrix weight on `I=[0,11/5]`,",
            "",
            "```text",
            "M_(p,b)=integral_I W_(p,b)",
            "Q_f=integral_I W_(p,b)f^2",
            artifact["compact_correction"],
            "```",
            "",
            "The remaining obligation is an Arb matrix for the stable finite",
            "remainder plus an explicit `u>11/5` tail. No retained first-jet",
            "lower separation or Newman conclusion is claimed.",
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
        "wrote Newman theta forward-remainder tail gate: "
        f"{len(artifact['rows'])} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
