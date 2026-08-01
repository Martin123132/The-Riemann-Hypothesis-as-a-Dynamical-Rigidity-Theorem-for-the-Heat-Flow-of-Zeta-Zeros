#!/usr/bin/env python3
"""Build the exact derivative envelope for the modular theta tail."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_theta_modular_tail_derivative_envelope_gate"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
RESULT_ROOT = REPO_ROOT / "work/rh_compute/results"
MODULAR_SOURCE = (
    RESULT_ROOT / "jensen_window_pf_newman_theta_modular_blend_gate.json"
)
C1_SOURCE = (
    RESULT_ROOT
    / "jensen_window_pf_newman_theta_adaptive_modular_c1_remainder_contract.json"
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def source_audit() -> dict:
    modular = load(MODULAR_SOURCE)
    c1 = load(C1_SOURCE)
    modular_text = json.dumps(modular, sort_keys=True)
    c1_text = json.dumps(c1, sort_keys=True)
    for marker in (
        "omega(u)=(1+erf(3*sinh(4u)))/2",
        "|R_(N,t)^(j)(x)|<=d_(N,j,m)(T)/|x|^m",
        "b_n(u)>0 for every n>=1",
    ):
        if marker not in modular_text:
            raise RuntimeError(f"modular source marker missing: {marker}")
    for marker in (
        "m>=5",
        "N_kappa(x,t)",
        "direct C1 sufficient disjunction",
    ):
        if marker not in c1_text:
            raise RuntimeError(f"C1 source marker missing: {marker}")
    return {
        "modular_kind": modular["kind"],
        "c1_kind": c1["kind"],
        "modular_rows": len(modular["rows"]),
        "c1_rows": len(c1["rows"]),
    }


def build_exact() -> dict:
    u, t, cap_t = sp.symbols("u t T", nonnegative=True, real=True)
    q = sp.symbols("q", positive=True, real=True)
    hermite: list[sp.Expr] = []
    for order in range(10):
        polynomial = sp.factor(
            sp.diff(sp.exp(t * u**2), u, order) / sp.exp(t * u**2)
        )
        formula = sp.factor(
            sp.factorial(order)
            * sum(
                (2 * u) ** (order - 2 * ell)
                * t ** (order - ell)
                / (
                    sp.factorial(ell)
                    * sp.factorial(order - 2 * ell)
                )
                for ell in range(order // 2 + 1)
            )
        )
        if sp.simplify(polynomial - formula) != 0:
            raise RuntimeError(f"Hermite envelope identity failed at {order}")
        hermite.append(formula)

    r = sp.Function("r")
    product_five = sp.diff(sp.exp(t * u**2) * r(u), u, 5)
    expected_five = sp.exp(t * u**2) * sum(
        sp.binomial(5, b) * hermite[b] * sp.diff(r(u), u, 5 - b)
        for b in range(6)
    )
    if sp.simplify(product_five - expected_five) != 0:
        raise RuntimeError("fifth derivative product identity failed")
    first_jet = sp.diff(u * sp.exp(t * u**2) * r(u), u, 5)
    expected_first_jet = (
        u * product_five
        + 5 * sp.diff(sp.exp(t * u**2) * r(u), u, 4)
    )
    if sp.simplify(first_jet - expected_first_jet) != 0:
        raise RuntimeError("first-jet fifth derivative identity failed")

    y_star = (8 * q / 9) ** sp.Rational(1, 3)
    phase = q / y_star + sp.Rational(9, 16) * y_star**2
    phase_expected = (
        sp.Rational(3, 2)
        * (sp.Rational(9, 8)) ** sp.Rational(1, 3)
        * q ** sp.Rational(2, 3)
    )
    if sp.simplify(phase - phase_expected) != 0:
        raise RuntimeError("reflected phase minimization failed")
    c0 = sp.factor(phase_expected.subs(q, sp.pi))
    threshold_slack = sp.simplify(
        sp.pi * 2**2 / sp.sqrt(2) - c0 * 2 ** sp.Rational(4, 3)
    )
    if not bool(sp.N(threshold_slack, 50) > 0):
        raise RuntimeError("small-y reflected phase threshold failed")
    c_safe = sp.factor(c0 / 4)
    gamma_k = sp.factor((sp.pi / (8 * c_safe)) ** sp.Rational(3, 4))

    return {
        "hermite_envelope": {
            "definition": (
                "P_b(T,v)=b!*sum_(ell=0)^floor(b/2) "
                "(2v)^(b-2ell)*T^(b-ell)/(ell!*(b-2ell)!)"
            ),
            "derivative_identity": (
                "partial_u^b exp(tu^2)=exp(tu^2)P_b(t,u)"
            ),
            "uniform_bound": (
                "|partial_u^b exp(tu^2)|<="
                "exp(Tu^2)P_b(T,|u|), 0<=t<=T"
            ),
            "verified_orders": list(range(10)),
        },
        "tail_envelopes": {
            "tail": "r_N(u)=sum_(n>N)b_n(u)",
            "E": (
                "E_(N,k)(u;T)=exp(Tu^2)*sum_(b=0)^k "
                "binom(k,b)P_b(T,u)|r_N^(k-b)(u)|, u>=0"
            ),
            "value_budget": (
                "d_(N,0,m)(T)<=integral_0^infinity E_(N,m)(u;T)du"
            ),
            "first_jet_budget": (
                "d_(N,1,m)(T)<=integral_0^infinity "
                "[u*E_(N,m)(u;T)+m*E_(N,m-1)(u;T)]du"
            ),
            "m5_value_integrand": (
                "exp(Tu^2)*sum_(b=0)^5 "
                "binom(5,b)P_b(T,u)|r_N^(5-b)(u)|"
            ),
            "m5_first_jet_integrand": (
                "u*E_(N,5)(u;T)+5*E_(N,4)(u;T)"
            ),
            "parity_reason": (
                "r_N is even, so one half of the real-line L1 norm is "
                "the positive-half-line integral."
            ),
        },
        "reflected_phase": {
            "variables": "y=exp(4u)>=1, q=pi*n^2",
            "switch_bound": (
                "1-omega(u)<=exp(-9*sinh(4u)^2)/2"
            ),
            "large_y": (
                "9*sinh(4u)^2>=9y^2/16 for y>=sqrt(2)"
            ),
            "phase": "q/y+9y^2/16",
            "minimizer": "y_*=(8q/9)^(1/3)",
            "minimum": (
                "(3/2)*(9/8)^(1/3)*q^(2/3)"
            ),
            "global_n_bound": (
                "pi*n^2/y+9*sinh(4u)^2>=c_0*n^(4/3), "
                "u>=0, n>=2"
            ),
            "c0_exact": str(c0),
            "c0_decimal": str(sp.N(c0, 30)),
        },
        "derivative_tail_theorem": {
            "statement": (
                "For fixed T,p,k there are effective C_(T,p,k)>0 and "
                "A_(p,k)>=0 such that integral_R |u|^p*exp(Tu^2)*"
                "|r_N^(k)(u)|du<=C_(T,p,k)*N^A*"
                "exp(-c_*N^(4/3)) for N>=2."
            ),
            "safe_exponent": "c_*=c_0/4",
            "c_safe_exact": str(c_safe),
            "c_safe_decimal": str(sp.N(c_safe, 30)),
            "proof": (
                "Differentiate the two blended summands. Forward terms are "
                "polynomials times exp(-pi*n^2*exp(4u)). Reflected terms are "
                "polynomials times exp(-pi*n^2*exp(-4u)-"
                "9*sinh(4u)^2). The phase bound supplies exp(-c_0*n^(4/3)); "
                "two successive halvings absorb all fixed derivative, "
                "exp(Tu^2), integration, and arithmetic-summation factors."
            ),
            "budget_consequence": (
                "For every fixed j,m,T, d_(N,j,m)(T)<="
                "C_(j,m,T)*N^A*exp(-c_*N^(4/3))."
            ),
        },
        "frequency_scales": {
            "saddle_count": (
                "N_sad(x,t)=ceil(max(n_*(5,t,x),n_*(9,t,x)))="
                "Theta(sqrt(x))"
            ),
            "saddle_bound_scale": (
                "The proved absolute derivative estimate at "
                "N_sad+O(1) is only polynomial(x)*exp(-c*x^(2/3))."
            ),
            "warning": (
                "This upper-bound scale does not prove failure of the exact "
                "Fourier remainder, but a fixed additive saddle collar is not "
                "by itself a gamma-scale absolute-error theorem."
            ),
            "gamma_compatible_count": (
                "N_K(x)=ceil(K*(1+x)^(3/4))"
            ),
            "gamma_compatible_bound": (
                "d_(N_K,j,m)(T)/x^m<=poly(x)*"
                "exp(-c_*K^(4/3)*x)"
            ),
            "safe_threshold": (
                "K>K_gamma=(pi/(8*c_*))^(3/4)"
            ),
            "K_gamma_exact": str(gamma_k),
            "K_gamma_decimal": str(sp.N(gamma_k, 30)),
            "combined_candidate": (
                "N(x,t)=max(N_sad(x,t),ceil(K*(1+x)^(3/4)))"
            ),
        },
        "remaining_target": {
            "statement": (
                "Instantiate effective constants for T=1/5 and m=5 or a "
                "higher optimized order, then prove the retained direct-C1 "
                "value-or-derivative separation on all transition cells."
            ),
            "not_supplied": (
                "The derivative envelope does not give a lower separation "
                "bound for (S_(N,t),S_(N,t)') and does not intervalize the "
                "adaptive transition cells."
            ),
        },
    }


def build_artifact() -> dict:
    exact = build_exact()
    audit = source_audit()
    rows = [
        GateRow(
            "ntmtdeg_01_heat_hermite_polynomials",
            "exact_identity",
            "available_exact",
            "All bounded-time heat-factor derivatives have an explicit positive polynomial envelope.",
            exact["hermite_envelope"]["definition"],
            "Exact scalar calculus, verified through the orders used by the scout.",
        ),
        GateRow(
            "ntmtdeg_02_value_derivative_envelope",
            "exact_inequality",
            "ready_to_apply",
            "The modular value remainder has an explicit positive-half-line derivative envelope.",
            exact["tail_envelopes"]["value_budget"],
            "The displayed integral still requires numerical or analytic enclosure.",
        ),
        GateRow(
            "ntmtdeg_03_first_jet_derivative_envelope",
            "exact_inequality",
            "ready_to_apply",
            "The modular first-derivative remainder has an explicit companion envelope.",
            exact["tail_envelopes"]["first_jet_budget"],
            "Exact Leibniz reduction; no retained first-jet separation follows.",
        ),
        GateRow(
            "ntmtdeg_04_reflected_phase_minimum",
            "exact_theorem",
            "ready_to_apply",
            "The reflected modular blend has a uniform n^(4/3) decay phase.",
            exact["reflected_phase"]["global_n_bound"],
            "Kernel-derivative decay only; not a Fourier sign theorem.",
            exact["reflected_phase"],
        ),
        GateRow(
            "ntmtdeg_05_stretched_exponential_tail",
            "exact_theorem",
            "ready_to_apply",
            "Every fixed weighted derivative norm of the omitted modular tail has an effective stretched-exponential bound.",
            exact["derivative_tail_theorem"]["statement"],
            "Constants are effective but not numerically instantiated here.",
            exact["derivative_tail_theorem"],
        ),
        GateRow(
            "ntmtdeg_06_fixed_collar_scale_guard",
            "nonpromotion_gate",
            "guard_validated",
            "A fixed additive saddle collar is not yet a gamma-scale absolute-error theorem.",
            exact["frequency_scales"]["saddle_bound_scale"],
            exact["frequency_scales"]["warning"],
        ),
        GateRow(
            "ntmtdeg_07_three_quarter_scale",
            "exact_theorem",
            "ready_to_apply",
            "An x^(3/4) block count converts the proved tail estimate to an exponential-in-x bound.",
            exact["frequency_scales"]["gamma_compatible_bound"],
            "Tail control only; it does not prove the direct-C1 disjunction.",
            exact["frequency_scales"],
        ),
        GateRow(
            "ntmtdeg_08_revised_adaptive_candidate",
            "exact_composition",
            "ready_to_apply",
            "The absolute-error architecture must retain both the saddle transition and the larger derivative-tail scale.",
            exact["frequency_scales"]["combined_candidate"],
            "A sufficient truncation architecture, not a proof that its retained jet stays separated.",
        ),
        GateRow(
            "ntmtdeg_09_open_separation_handoff",
            "open_handoff",
            "not_ready_to_apply",
            exact["remaining_target"]["statement"],
            exact["remaining_target"]["not_supplied"],
            "No new Q_j, Lambda<=0, RH, or Clay-prize result is proved.",
        ),
    ]
    return {
        "kind": STEM,
        "date": "2026-07-24",
        "status": (
            "exact modular-tail derivative envelope and effective "
            "three-quarter-scale absolute-error theorem"
        ),
        "proof_boundary": (
            "This artifact proves explicit integration-by-parts envelopes and "
            "an effective N^A*exp(-c*N^(4/3)) bound for every fixed modular-tail "
            "derivative norm. It derives an x^(3/4) sufficient scale for "
            "exponential absolute tail control and guards against promoting a "
            "fixed additive saddle collar. It does not instantiate referee-ready "
            "constants, prove retained first-jet separation, strict Laguerre "
            "positivity, Lambda<=0, RH, or a Clay-prize result."
        ),
        "source_audit": audit,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "sources": [
            "outputs/jensen_window_pf_newman_theta_modular_blend_gate.md",
            "outputs/jensen_window_pf_newman_theta_adaptive_modular_c1_remainder_contract.md",
            "outputs/formal_core.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return "\n".join(
        [
            "# Newman Theta Modular-Tail Derivative Envelope Gate",
            "",
            "Date: 2026-07-24",
            "",
            "Status: exact derivative envelope and effective tail-scale theorem.",
            "This is not a proof of `Lambda<=0`, RH, or a Clay-prize result.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Explicit Heat-Derivative Envelope",
            "",
            "For `v>=0`, define",
            "",
            "```text",
            exact["hermite_envelope"]["definition"],
            exact["hermite_envelope"]["uniform_bound"],
            "```",
            "",
            "With `r_N=sum_(n>N)b_n`, put",
            "",
            "```text",
            exact["tail_envelopes"]["E"],
            exact["tail_envelopes"]["value_budget"],
            exact["tail_envelopes"]["first_jet_budget"],
            "```",
            "",
            "The second formula is the exact identity",
            "`D^m[u*f]=u*D^m[f]+m*D^(m-1)[f]`, followed by the",
            "triangle inequality. Evenness turns half of the real-line norm",
            "into the displayed positive-half-line integral.",
            "",
            "## Reflected-Blend Phase",
            "",
            "Set `y=exp(4u)` and `q=pi*n^2`. For `y>=sqrt(2)`,",
            "",
            "```text",
            exact["reflected_phase"]["large_y"],
            exact["reflected_phase"]["phase"],
            exact["reflected_phase"]["minimizer"],
            exact["reflected_phase"]["minimum"],
            "```",
            "",
            "On `1<=y<=sqrt(2)`, the forward `q/y` term is stronger.",
            "Consequently, for every `u>=0` and `n>=2`,",
            "",
            "```text",
            exact["reflected_phase"]["global_n_bound"],
            f"c_0={exact['reflected_phase']['c0_decimal']}",
            "```",
            "",
            "Differentiating either blended summand only introduces fixed",
            "polynomial factors. Splitting the exponential phase twice absorbs",
            "those factors, the bounded Newman heat weight, integration, and",
            "the arithmetic tail sum. Thus",
            "",
            "```text",
            exact["derivative_tail_theorem"]["statement"],
            exact["derivative_tail_theorem"]["budget_consequence"],
            f"c_*={exact['derivative_tail_theorem']['c_safe_decimal']}",
            "```",
            "",
            "## Correct Absolute-Error Scale",
            "",
            "The saddle transition still identifies which arithmetic blocks",
            "matter spectrally, but the proved absolute derivative envelope has",
            "a different scale:",
            "",
            "```text",
            exact["frequency_scales"]["saddle_count"],
            exact["frequency_scales"]["saddle_bound_scale"],
            exact["frequency_scales"]["gamma_compatible_count"],
            exact["frequency_scales"]["gamma_compatible_bound"],
            exact["frequency_scales"]["safe_threshold"],
            f"K_gamma={exact['frequency_scales']['K_gamma_decimal']}",
            exact["frequency_scales"]["combined_candidate"],
            "```",
            "",
            "The saddle-scale statement is a non-promotion guard: a slower",
            "available upper bound does not prove that the exact oscillatory",
            "remainder is large. It does show that a fixed additive collar",
            "does not itself supply the required absolute-error theorem.",
            "",
            "## Open Separation Target",
            "",
            exact["remaining_target"]["statement"],
            "",
            exact["remaining_target"]["not_supplied"],
            "",
            "No new `Q_j`, strict Laguerre theorem, `Lambda<=0`, RH, or",
            "Clay-prize conclusion is claimed.",
            "",
        ]
    )


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
        "wrote Newman theta modular-tail derivative envelope gate: "
        f"{len(artifact['rows'])} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
