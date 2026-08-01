#!/usr/bin/env python3
"""Build the exact arbitrary-N stable forward-remainder theorem."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_theta_arbitrary_n_stable_remainder_gate"
DEFAULT_OUT = REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
MODULAR_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_modular_blend_gate.json"
)
FINITE_PREDECESSOR = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_forward_remainder_tail_gate.json"
)
DATE = "2026-07-24"
MAX_ORDER = 9
MIN_RETAINED = 1
MIN_FIRST_OMITTED = 2
WITNESS_K = (2, 3, 4, 5, 8, 11, 13, 21, 33, 65)


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
        coefficient_norm = int(
            sum(abs(value) for value in polynomial.all_coeffs())
        )
        rows.append(
            {
                "q": order,
                "P_q": str(polynomial.as_expr()),
                "degree": polynomial.degree(),
                "A_q": coefficient_norm,
            }
        )
        current = sp.expand(
            (5 - 4 * x) * polynomial.as_expr()
            + 4 * x * sp.diff(polynomial.as_expr(), x)
        )
    return rows


def exact_tail_bound(q: int, coefficient_norm: int, k: int) -> tuple:
    d = 2 * q + 4
    rho = sp.exp(sp.Rational(d, k) - sp.pi * (2 * k + 1))
    bound = (
        coefficient_norm
        * sp.pi ** (q + 2)
        * k**d
        * sp.exp(-sp.pi * k**2)
        / (1 - rho)
    )
    return rho, bound


def witness_rows(polynomials: list[dict]) -> list[dict]:
    rows: list[dict] = []
    for k in WITNESS_K:
        for polynomial in polynomials:
            q = polynomial["q"]
            rho, bound = exact_tail_bound(q, polynomial["A_q"], k)
            rows.append(
                {
                    "K": k,
                    "q": q,
                    "rho_exact": str(rho),
                    "rho_decimal": str(sp.N(rho, 30)),
                    "B_qK_exact": str(bound),
                    "B_qK_decimal": str(sp.N(bound, 30)),
                    "log10_B_qK": str(sp.N(sp.log(bound, 10), 30)),
                }
            )
    return rows


def source_audit() -> dict:
    modular = json.loads(MODULAR_SOURCE.read_text(encoding="utf-8"))
    exact = modular.get("exact", {})
    identity = exact.get("decaying_tail_enclosure", {}).get("tail_kernel")
    required = (
        "r_N(u)=sum_(n>N)b_n(u)=Phi(u)-sum_(n<=N)b_n(u)"
    )
    if identity != required:
        raise RuntimeError(f"modular tail identity drifted: {identity!r}")
    series = exact.get("theta_summand", {}).get("kernel_series", "")
    if "Phi(u)=sum_(n>=1)phi_n(u)" not in series:
        raise RuntimeError("forward theta series marker missing")

    predecessor = json.loads(
        FINITE_PREDECESSOR.read_text(encoding="utf-8")
    )
    if predecessor.get("kind") != (
        "jensen_window_pf_newman_theta_forward_remainder_tail_gate"
    ):
        raise RuntimeError("finite predecessor kind drifted")
    return {
        "modular_source": {
            "kind": modular["kind"],
            "sha256": file_hash(MODULAR_SOURCE),
            "tail_identity": identity,
            "forward_series": series,
        },
        "finite_predecessor": {
            "kind": predecessor["kind"],
            "sha256": file_hash(FINITE_PREDECESSOR),
            "first_omitted": predecessor["parameters"]["first_omitted"],
        },
    }


def build_artifact() -> dict:
    polynomials = polynomial_rows()
    witnesses = witness_rows(polynomials)
    exact_split = (
        "delta_n(u)=phi_n(u)-b_n(u)="
        "(erfc(3*sinh(4u))/2)*(phi_n(u)-phi_n(-u)); "
        "r_N^(q)(u)=sum_(n=1)^N D^q delta_n(u)"
        "+sum_(n=N+1)^M D^q phi_n(u)+T_(M+1,q)(u), "
        "T_(K,q)=sum_(n>=K)D^q phi_n, "
        "N>=1, M>=N, K=M+1"
    )
    cap_free_split = (
        "r_N^(q)(u)=sum_(n=1)^N D^q delta_n(u)"
        "+T_(N+1,q)(u)"
    )
    tail_formula = (
        "B_(q,K)=A_q*pi^(q+2)*K^(2q+4)*exp(-pi*K^2)"
        "/(1-exp((2q+4)/K-pi*(2K+1)))"
    )
    rows = [
        GateRow(
            id="ntansrg_01_all_n_modular_identity",
            role="exact_identity",
            readiness="proved",
            claim=(
                "The cancellation-stable modular remainder identity holds "
                "for every integer retained count N>=1."
            ),
            formula=(
                "r_N=sum_(n>N)b_n=Phi-sum_(n<=N)b_n="
                "sum_(n>N)phi_n+sum_(n<=N)delta_n"
            ),
            proof_boundary=(
                "This is an algebraic consequence of the exact positive "
                "modular partition; it does not estimate delta_n."
            ),
        ),
        GateRow(
            id="ntansrg_02_arbitrary_cap_split",
            role="exact_identity",
            readiness="proved",
            claim=(
                "Any cap M>=N gives an exact finite stable block and a "
                "forward tail beginning at K=M+1."
            ),
            formula=exact_split,
            proof_boundary=(
                "Termwise differentiation is justified by locally uniform "
                "superexponential convergence of the forward theta series."
            ),
        ),
        GateRow(
            id="ntansrg_03_cap_free_split",
            role="exact_identity",
            readiness="proved",
            claim=(
                "Taking M=N removes the artificial fixed cap completely."
            ),
            formula=cap_free_split,
            proof_boundary=(
                "The finite switch-defect derivatives remain to be bounded "
                "inside the full d0/d1 construction."
            ),
        ),
        GateRow(
            id="ntansrg_04_forward_derivative_recurrence",
            role="exact_identity",
            readiness="proved",
            claim=(
                "Every forward derivative is one Gaussian times an integer "
                "polynomial generated by a first-order recurrence."
            ),
            formula=(
                "D^q phi_n=pi*n^2*exp(5u)*P_q(X)*exp(-X), "
                "X=pi*n^2*exp(4u); P_0=2X-3; "
                "P_(q+1)=(5-4X)P_q+4X P_q'"
            ),
            proof_boundary="Exact through q=9.",
            diagnostics=polynomials,
        ),
        GateRow(
            id="ntansrg_05_uniform_endpoint_domination",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "For every K>=2, n>=K, u>=0, and q<=9, the coefficient-"
                "norm forward derivative majorant is maximized at u=0."
            ),
            formula=(
                "|D^q phi_n(u)|<=A_q*pi^(q+2)*"
                "n^(2q+4)*exp(-pi*n^2)"
            ),
            proof_boundary=(
                "The worst exponent derivative is "
                "45-16*pi<45-48<0, using pi>3."
            ),
        ),
        GateRow(
            id="ntansrg_06_uniform_geometric_ratio",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "The endpoint majorant sequence has a uniform geometric "
                "ratio below one for every K>=2 and q<=9."
            ),
            formula=(
                "a_(n+1)/a_n<=rho_(q,K)="
                "exp((2q+4)/K-pi*(2K+1))"
                "<=exp(11-5*pi)<exp(-4)<1"
            ),
            proof_boundary=(
                "Uses log(1+1/n)<=1/n and monotonicity in n, K, and q."
            ),
        ),
        GateRow(
            id="ntansrg_07_explicit_arbitrary_tail",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "Every arbitrary forward tail has a closed explicit uniform "
                "bound B_(q,K), with no finite upper bound on K or N."
            ),
            formula=(
                "sup_(u>=0)|T_(K,q)(u)|<=B_(q,K); " + tail_formula
            ),
            proof_boundary=(
                "Valid for every integer K>=2 and 0<=q<=9. The witness "
                "table is diagnostic; the theorem itself is symbolic."
            ),
            diagnostics=witnesses,
        ),
        GateRow(
            id="ntansrg_08_adaptive_cofinal_handoff",
            role="proof_search_target",
            readiness="not_ready_to_apply",
            claim=(
                "For an adaptive N(x,t), the forward arithmetic tail is "
                "controlled explicitly by B_(q,N+1); the remaining cofinal "
                "work is the switch defect and retained first-jet margin."
            ),
            formula=(
                "N(x,t)>=max(1,ceil(kappa*(1+x)^(3/4))) "
                "implies sup_(u>=0)|T_(N+1,q)|<=B_(q,N+1)"
            ),
            proof_boundary=(
                "This does not bound sum_(n<=N)D^q delta_n, give full "
                "d0/d1 budgets, prove a retained J/J' lower separation, "
                "cover all x>38, prove Lambda<=0 or RH, or supply a "
                "Clay-prize proof."
            ),
        ),
    ]
    return {
        "kind": STEM,
        "date": DATE,
        "status": "exact arbitrary-N stable forward-remainder theorem",
        "source_audit": source_audit(),
        "parameters": {
            "max_derivative_order": MAX_ORDER,
            "minimum_retained_count": MIN_RETAINED,
            "minimum_first_omitted": MIN_FIRST_OMITTED,
            "retained_count_upper_bound": None,
            "cap_upper_bound": None,
            "witness_K": list(WITNESS_K),
        },
        "stable_remainder": exact_split,
        "cap_free_remainder": cap_free_split,
        "tail_bound_formula": tail_formula,
        "endpoint_extrema": {
            "u_monotonicity_worst_case": "16*pi-45>48-45>0",
            "ratio_exponent_worst_case": "11-5*pi<11-15=-4",
        },
        "polynomials": polynomials,
        "witness_bounds": witnesses,
        "rows": [asdict(row) for row in rows],
        "proof_boundary": (
            "This artifact removes the artificial N<=12 forward-tail cap "
            "and proves explicit uniform q<=9 forward arithmetic-tail "
            "constants for every N>=1. It does not yet bound the finite "
            "switch-defect sum, assemble arbitrary-N full d0/d1 budgets, "
            "prove retained first-jet separation, establish a terminating "
            "cofinal cover for x>38, prove Lambda<=0 or RH, or supply a "
            "Clay-prize proof."
        ),
    }


def render_note(artifact: dict) -> str:
    lines = [
        "# Newman Theta Arbitrary-N Stable Remainder Gate",
        "",
        f"Date: {DATE}",
        "",
        "Status: exact arbitrary-`N` forward-remainder theorem.",
        "This is not the cofinal retained-separation theorem and not a proof",
        "of `Lambda<=0`, RH, or a Clay-prize result.",
        "",
        "## Exact Split",
        "",
        "For every integer `N>=1`,",
        "",
        "```text",
        "r_N=sum_(n>N)b_n",
        "   =sum_(n>N)phi_n+sum_(n<=N)delta_n",
        "delta_n=(erfc(3*sinh(4u))/2)*(phi_n(u)-phi_n(-u)).",
        "```",
        "",
        "For every `M>=N`, differentiation through order nine gives",
        "",
        "```text",
        artifact["stable_remainder"],
        "```",
        "",
        "Taking `M=N` removes the former fixed cap:",
        "",
        "```text",
        artifact["cap_free_remainder"],
        "```",
        "",
        "## Uniform Tail",
        "",
        "With `X=pi*n^2*exp(4u)`,",
        "",
        "```text",
        "D^q phi_n=pi*n^2*exp(5u)*P_q(X)*exp(-X)",
        "P_0=2X-3",
        "P_(q+1)=(5-4X)P_q+4X P_q'.",
        "```",
        "",
        "If `A_q` is the coefficient 1-norm of `P_q`, then for every",
        "`K>=2`, `n>=K`, `u>=0`, and `0<=q<=9`,",
        "",
        "```text",
        "|D^q phi_n(u)|",
        " <=A_q*pi^(q+2)*n^(2q+4)*exp(-pi*n^2).",
        "rho_(q,K)=exp((2q+4)/K-pi*(2K+1))<1.",
        artifact["tail_bound_formula"],
        "sup_(u>=0)|T_(K,q)(u)|<=B_(q,K).",
        "```",
        "",
        "The two uniform endpoint checks reduce exactly to",
        "`16*pi-45>48-45>0` and `11-5*pi<11-15=-4`.",
        "",
        "## Witness Scale",
        "",
        "| K | q | log10(B_(q,K)) |",
        "|---:|---:|---:|",
    ]
    selected = {(2, 9), (3, 9), (4, 9), (8, 9), (13, 9), (33, 9), (65, 9)}
    for row in artifact["witness_bounds"]:
        if (row["K"], row["q"]) in selected:
            lines.append(
                f"| {row['K']} | {row['q']} | "
                f"`{row['log10_B_qK']}` |"
            )
    lines.extend(
        [
            "",
            "The table is only a scale check. The theorem is the symbolic",
            "all-`K` inequality above.",
            "",
            "## Cofinal Boundary",
            "",
            "For an adaptive retained count, the forward arithmetic tail is",
            "now explicit as `B_(q,N+1)`. The unresolved terms are the full",
            "finite switch-defect derivative sum and, more importantly, a",
            "uniform retained value-or-derivative lower margin. Those are",
            "separate gates and are not inferred from Q31 or this theorem.",
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
        "wrote Newman theta arbitrary-N stable remainder gate: "
        f"{len(artifact['rows'])} rows, "
        f"{len(artifact['witness_bounds'])} witness bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
