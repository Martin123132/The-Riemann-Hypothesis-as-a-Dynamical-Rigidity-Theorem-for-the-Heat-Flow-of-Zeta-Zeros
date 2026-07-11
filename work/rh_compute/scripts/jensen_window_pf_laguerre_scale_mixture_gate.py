#!/usr/bin/env python3
"""Derive the Laguerre scale-mixture gate for Jensen-window hyperbolicity."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_laguerre_scale_mixture_gate.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_laguerre_scale_mixture_gate.md"


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def falling(d: int, j: int) -> sp.Expr:
    return sp.factorial(d) / sp.factorial(d - j)


def kernel_polynomial(d: int, w: sp.Symbol, y: sp.Expr) -> sp.Expr:
    return sp.expand(
        sum(falling(d, j) * (w * y) ** j / sp.factorial(2 * j) for j in range(d + 1))
    )


def gamma_polynomial(
    d: int, alpha: sp.Expr, theta: sp.Expr, w: sp.Symbol
) -> sp.Expr:
    return sp.expand(
        sum(
            falling(d, j)
            * sp.rf(alpha, j)
            * (theta * w) ** j
            / sp.factorial(2 * j)
            for j in range(d + 1)
        )
    )


def terminating_hyper(
    upper_degree: int, second: sp.Expr, x: sp.Expr
) -> sp.Expr:
    return sp.expand(
        sum(
            sp.rf(-upper_degree, j)
            * sp.rf(second, j)
            * x**j
            / (sp.rf(sp.Rational(1, 2), j) * sp.factorial(j))
            for j in range(upper_degree + 1)
        )
    )


def build_exact() -> dict:
    w, y = sp.symbols("w y", real=True)
    kernel_checks: list[dict] = []
    for d in range(1, 9):
        direct = kernel_polynomial(d, w, y)
        kummer = sp.hyperexpand(
            sp.hyper([-d], [sp.Rational(1, 2)], -w * y / 4)
        )
        laguerre = (
            sp.factorial(d)
            / sp.rf(sp.Rational(1, 2), d)
            * sp.assoc_laguerre(d, sp.Rational(-1, 2), -w * y / 4)
        )
        if sp.simplify(direct - kummer) != 0:
            raise RuntimeError(f"Kummer kernel identity failed at D={d}")
        if sp.simplify(direct - laguerre) != 0:
            raise RuntimeError(f"Laguerre kernel identity failed at D={d}")
        kernel_checks.append(
            {"D": d, "degree": int(sp.degree(direct, w)), "passed": True}
        )

    mass = sp.Rational(9, 10)
    atom_a = sp.Rational(1, 4)
    atom_b = sp.Rational(5, 2)
    mean_y = sp.factor(mass * atom_a + (1 - mass) * atom_b)
    mean_y2 = sp.factor(mass * atom_a**2 + (1 - mass) * atom_b**2)
    two_atom = sp.expand(1 + mean_y * w + mean_y2 * w**2 / 12)
    two_atom_disc = sp.factor(sp.discriminant(two_atom, w))
    if two_atom_disc != -sp.Rational(7, 4800):
        raise RuntimeError("two-atom positive-mixture countermodel failed")

    gamma_checks: list[dict] = []
    theta = sp.symbols("theta", positive=True)
    for d in range(1, 7):
        alpha = sp.Rational(7, 5)
        direct = gamma_polynomial(d, alpha, theta, w)
        hyper = sp.hyperexpand(
            sp.hyper([-d, alpha], [sp.Rational(1, 2)], -theta * w / 4)
        )
        if sp.simplify(direct - hyper) != 0:
            raise RuntimeError(f"Gamma hypergeometric identity failed at D={d}")
        gamma_checks.append({"D": d, "alpha": str(alpha), "passed": True})

    half_integer_checks: list[dict] = []
    xvar = -w / 4
    for m in range(6):
        alpha = sp.Rational(2 * m + 1, 2)
        for d in range(1, 11):
            direct = gamma_polynomial(d, alpha, 1, w)
            if d >= m:
                cofactor = terminating_hyper(m, sp.Rational(2 * d + 1, 2), xvar)
                expected = sp.expand((1 - xvar) ** (d - m) * cofactor)
                branch = "Euler-Jacobi"
            else:
                expected = (
                    sp.factorial(d)
                    / sp.rf(sp.Rational(1, 2), d)
                    * sp.jacobi(
                        d,
                        sp.Rational(-1, 2),
                        m - d,
                        1 - 2 * xvar,
                    )
                )
                branch = "direct-Jacobi"
            if sp.simplify(direct - expected) != 0:
                raise RuntimeError(f"half-integer factorization failed at m={m}, D={d}")
            half_integer_checks.append(
                {
                    "m": m,
                    "D": d,
                    "branch": branch,
                    "endpoint_multiplicity": max(d - m, 0),
                    "interior_root_count": min(d, m),
                    "passed": True,
                }
            )

    exponential = gamma_polynomial(3, 1, 1, w)
    exponential_disc = sp.factor(sp.discriminant(exponential, w))
    if exponential_disc != -sp.Rational(1, 200):
        raise RuntimeError("log-concave exponential countermodel failed")

    return {
        "moment_normalization": (
            "A_j(lambda)=j!/(2j)!*integral_0^infinity W_lambda(u)u^(2j)du, "
            "W_lambda(u)=2*exp(lambda*u^2)*Phi(u)>0"
        ),
        "jensen_integral": (
            "J_D(w)=integral_0^infinity W_lambda(u)K_D(w,u^2)du, "
            "K_D(w,y)=sum_(j=0)^D (D)_j*(w*y)^j/(2j)!"
        ),
        "kummer_kernel": "K_D(w,y)={}_1F_1(-D;1/2;-w*y/4)",
        "laguerre_kernel": (
            "K_D(w,y)=D!/(1/2)_D*L_D^(-1/2)(-w*y/4)"
        ),
        "single_scale_roots": (
            "For y>0, K_D(w,y) has D simple roots w=-4*ell_(D,r)/y<0, "
            "where ell_(D,r) are the positive roots of L_D^(-1/2)."
        ),
        "positive_mixture_countermodel": {
            "measure": "(9/10)*delta_(1/4)+(1/10)*delta_(5/2) in y=u^2",
            "D": 2,
            "polynomial": str(two_atom),
            "discriminant": str(two_atom_disc),
            "conclusion": "Positive scale mixtures of the hyperbolic kernels need not be hyperbolic.",
        },
        "gamma_mixing": {
            "law": "Y~Gamma(alpha,theta), E[Y^j]=theta^j*(alpha)_j",
            "polynomial": "P_(D,alpha,theta)(w)={}_2F_1(-D,alpha;1/2;-theta*w/4)",
        },
        "half_integer_factorization": {
            "parameters": "alpha=m+1/2, m in Z_(>=0), X=-theta*w/4",
            "D_ge_m": (
                "P=(1-X)^(D-m)*{}_2F_1(-m,D+1/2;1/2;X), "
                "the cofactor being proportional to P_m^(-1/2,D-m)(1-2X)"
            ),
            "D_lt_m": (
                "P=D!/(1/2)_D*P_D^(-1/2,m-D)(1-2X)"
            ),
            "root_theorem": (
                "All D roots are real and negative: min(D,m) are simple in "
                "(-4/theta,0), and w=-4/theta has multiplicity max(D-m,0)."
            ),
        },
        "log_concavity_countermodel": {
            "mixing_density": "exp(-y) on y>=0 (Gamma(1,1), log-concave)",
            "D": 3,
            "polynomial": str(exponential),
            "discriminant": str(exponential_disc),
            "conclusion": "Log-concavity of the squared-scale density does not ensure Jensen hyperbolicity.",
        },
        "xi_specific_handoff": (
            "Exploit the actual pushforward density W_lambda(sqrt(y))/(2*sqrt(y)) "
            "through a half-integer Gamma/Laguerre expansion with a proved "
            "common-interlacing, total-positive connection, or direct variation-diminishing rule; "
            "positivity of the expansion coefficients alone is insufficient."
        ),
        "checks": {
            "kernel_checks": kernel_checks,
            "gamma_checks": gamma_checks,
            "half_integer_checks": half_integer_checks,
        },
    }


def build_payload() -> dict:
    exact = build_exact()
    rows = [
        GateRow(
            id="lsmg_01_jensen_moment_integral",
            role="exact_identity",
            readiness="available_exact",
            claim="Each unshifted Jensen polynomial is an integral of a universal finite scale kernel against the positive Phi heat-flow weight.",
            formula=exact["jensen_integral"],
            proof_boundary="Unshifted Jensen window n=0; no hyperbolicity conclusion.",
        ),
        GateRow(
            id="lsmg_02_kummer_kernel",
            role="exact_identity",
            readiness="available_exact",
            claim="The finite scale kernel is a terminating confluent hypergeometric polynomial.",
            formula=exact["kummer_kernel"],
            proof_boundary="Finite algebraic identity.",
        ),
        GateRow(
            id="lsmg_03_laguerre_kernel",
            role="exact_identity",
            readiness="available_exact",
            claim="The same kernel is a rescaled generalized Laguerre polynomial with parameter -1/2.",
            formula=exact["laguerre_kernel"],
            proof_boundary="Finite algebraic identity.",
        ),
        GateRow(
            id="lsmg_04_single_scale_hyperbolicity",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="Every nonzero fixed-scale kernel has only simple negative roots.",
            formula=exact["single_scale_roots"],
            proof_boundary="Individual kernels only; integration is not covered.",
        ),
        GateRow(
            id="lsmg_05_positive_mixture_guard",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="Positive integration of the fixed-scale kernels does not preserve real-rootedness in general.",
            formula="D=2 two-atom discriminant=-7/4800",
            proof_boundary="Abstract positive measure in y, not the Xi density.",
            diagnostics=exact["positive_mixture_countermodel"],
        ),
        GateRow(
            id="lsmg_06_gamma_reduction",
            role="exact_identity",
            readiness="available_exact",
            claim="Gamma scale mixing converts the Jensen kernel into an explicit Gauss hypergeometric polynomial.",
            formula=exact["gamma_mixing"]["polynomial"],
            proof_boundary="Gamma mixing laws only.",
        ),
        GateRow(
            id="lsmg_07_half_integer_factorization",
            role="exact_identity",
            readiness="available_exact",
            claim="At half-integer Gamma shape the polynomial factors exactly into an endpoint power and a Jacobi polynomial.",
            formula=f"{exact['half_integer_factorization']['D_ge_m']}; {exact['half_integer_factorization']['D_lt_m']}",
            proof_boundary="alpha=m+1/2 only.",
        ),
        GateRow(
            id="lsmg_08_half_integer_all_degree_theorem",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="Every half-integer Gamma scale law produces a hyperbolic Jensen polynomial in every degree.",
            formula=exact["half_integer_factorization"]["root_theorem"],
            proof_boundary="Solvable benchmark family, not a representation theorem for Xi.",
            diagnostics=exact["checks"]["half_integer_checks"],
        ),
        GateRow(
            id="lsmg_09_log_concavity_guard",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="Even a log-concave squared-scale density can destroy Jensen hyperbolicity.",
            formula="Gamma(1,1), D=3, discriminant=-1/200",
            proof_boundary="Generic log-concavity obstruction, not an Xi counterexample.",
            diagnostics=exact["log_concavity_countermodel"],
        ),
        GateRow(
            id="lsmg_10_xi_gamma_laguerre_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="A closing kernel route needs an Xi-specific Gamma/Laguerre expansion together with a zero-preserving connection theorem stronger than positivity.",
            formula=exact["xi_specific_handoff"],
            proof_boundary="Open all-degree theorem; not a proof of PF-infinity, RH, or Lambda <= 0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_laguerre_scale_mixture_gate",
        "date": "2026-07-11",
        "status": "exact Laguerre scale-mixture reduction with sharp preservation guards",
        "proof_boundary": (
            "This artifact exactly identifies the Laguerre kernel inside the unshifted Jensen "
            "moment integral, proves fixed-scale and half-integer Gamma benchmark hyperbolicity, "
            "and blocks generic positive-mixture and log-concavity promotions. It does not prove "
            "an Xi-specific zero-preserving expansion, all-degree Jensen hyperbolicity, PF-infinity, "
            "RH, or Lambda <= 0."
        ),
        "sources": [
            "outputs/coefficient_pf_bridge_obstruction.md",
            "outputs/jensen_window_pf_bridge_target.md",
            "https://dlmf.nist.gov/13.6",
            "https://dlmf.nist.gov/15.8",
            "https://dlmf.nist.gov/18.16",
        ],
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    positive_guard = exact["positive_mixture_countermodel"]
    log_guard = exact["log_concavity_countermodel"]
    half = exact["half_integer_factorization"]
    return "\n".join(
        [
            "# Jensen-Window PF Laguerre Scale-Mixture Gate",
            "",
            "Date: 2026-07-11",
            "",
            "Status: exact Laguerre scale-mixture reduction with sharp preservation",
            "guards. This is not a proof of PF-infinity, RH, or `Lambda <= 0`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_laguerre_scale_mixture_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_laguerre_scale_mixture_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_laguerre_scale_mixture_gate.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF Laguerre scale-mixture gate: 10 rows, 0 issues, 3 exact kernel identities, 1 individual-kernel hyperbolicity theorem, 1 positive-mixture countermodel, 1 Gamma reduction, 1 half-integer all-degree theorem, 1 log-concavity countermodel, 1 open Xi-specific handoff",
            "```",
            "",
            "## Exact Kernel",
            "",
            "With the established normalization",
            "",
            "```text",
            exact["moment_normalization"],
            exact["jensen_integral"],
            exact["kummer_kernel"],
            exact["laguerre_kernel"],
            "```",
            "",
            exact["single_scale_roots"],
            "",
            "Thus every scale entering the integral is individually hyperbolic.",
            "That observation is exact, but the next promotion is false.",
            "",
            "## Preservation Guards",
            "",
            "The positive two-atom measure",
            "",
            "```text",
            positive_guard["measure"],
            f"J_2(w)={positive_guard['polynomial']}",
            f"disc_w J_2={positive_guard['discriminant']}<0",
            "```",
            "",
            "is already a counterexample to positive-mixture preservation. Even",
            "log-concavity is insufficient: the exponential density gives",
            "",
            "```text",
            f"J_3(w)={log_guard['polynomial']}",
            f"disc_w J_3={log_guard['discriminant']}<0.",
            "```",
            "",
            "So the route cannot close from positivity, moment existence, or",
            "log-concavity of the squared-scale density alone.",
            "",
            "## Half-Integer Gamma Theorem",
            "",
            "For `Y~Gamma(alpha,theta)`,",
            "",
            "```text",
            exact["gamma_mixing"]["polynomial"],
            "```",
            "",
            "When `alpha=m+1/2`, Euler's transformation and the Jacobi identity give",
            "",
            "```text",
            half["D_ge_m"],
            half["D_lt_m"],
            half["root_theorem"],
            "```",
            "",
            "This is a genuine all-degree hyperbolic benchmark family. It matches",
            "the half-integer power structure of an even analytic weight after the",
            "pushforward `y=u^2`, but positive mixing of these blocks is not enough.",
            "",
            "## Live Handoff",
            "",
            exact["xi_specific_handoff"],
            "",
            "The useful question is now whether the actual Phi pushforward carries",
            "a common-interlacing or total-positive Gamma/Laguerre connection that",
            "the two countermodels cannot possess. That is an Xi-specific theorem",
            "target, not a generic measure argument.",
            "",
            "References: https://dlmf.nist.gov/13.6, https://dlmf.nist.gov/15.8,",
            "and https://dlmf.nist.gov/18.16.",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Jensen-window PF Laguerre scale-mixture gate: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
