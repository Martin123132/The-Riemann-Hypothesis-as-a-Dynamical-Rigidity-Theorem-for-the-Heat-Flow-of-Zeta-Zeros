#!/usr/bin/env python3
"""Build the fixed-omega Suzuki determinant phase diagram."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_suzuki_fixed_omega_phase_diagram.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_suzuki_fixed_omega_phase_diagram.md"
)


@dataclass(frozen=True)
class PhaseRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def exact_statements() -> dict[str, str]:
    return {
        "determinant_property": (
            "D(omega): for one (equivalently every) integer nu with "
            "nu*omega>1, det(I+/-K_(omega,nu)[t])!=0 for every t>=0"
        ),
        "fixed_pair_equivalence": (
            "D(omega) iff Theta_omega is meromorphic inner in C+ iff "
            "m_xi(rho-2*omega)>=m_xi(rho) for every xi-zero rho with "
            "Re(rho)>1/2+omega"
        ),
        "power_independence": (
            "Theta_omega^nu is inner for one positive integer nu iff "
            "Theta_omega is inner iff Theta_omega^mu is inner for every "
            "positive integer mu"
        ),
        "boundary_unimodularity": (
            "Theta_omega(u)=conj(xi(1/2+omega-i*u))/"
            "xi(1/2+omega-i*u), so |Theta_omega(u)|=1 wherever defined"
        ),
        "pole_coordinate": (
            "rho=beta+i*gamma maps to "
            "z_rho=-gamma+i*(beta-1/2-omega); "
            "the numerator there is xi(rho-2*omega)"
        ),
        "multiplicity_rule": (
            "ord_pole(z_rho)=max(0,m_xi(rho)-m_xi(rho-2*omega))"
        ),
        "no_poles_to_inner": (
            "For the xi quotient, no poles in C+ plus Suzuki's high-strip "
            "bound and Phragmen-Lindelof imply Theta_omega in H-infinity(C+)"
        ),
        "inner_to_determinants": (
            "Theta_omega inner gives a full L2 isometry; Suzuki's "
            "noncompact-support lemma and compact truncations give "
            "||K_(omega,nu)[t]||<1 for every finite t"
        ),
        "determinants_to_inner": (
            "All finite strict contractions extend to a bounded causal "
            "convolution multiplier equal to Theta_omega^nu, hence "
            "Theta_omega is meromorphic inner"
        ),
        "cancellation_set": (
            "C_xi={(Re(rho)-Re(rho'))/2>0: xi(rho)=xi(rho')=0 and "
            "Im(rho)=Im(rho')}; C_xi is countable"
        ),
        "zero_width": (
            "delta_xi=sup_{xi(rho)=0}(Re(rho)-1/2) in [0,1/2]"
        ),
        "phase_diagram": (
            "[delta_xi,infinity) subset D subset "
            "[delta_xi,infinity) union (C_xi intersect (0,delta_xi))"
        ),
        "local_exclusion": (
            "If RH is false, one off-line zero rho and zero isolation give "
            "epsilon_rho>0 such that D(omega) fails for every "
            "0<omega<epsilon_rho"
        ),
        "closure_criterion": (
            "RH iff 0 belongs to closure(D) iff D contains a sequence "
            "omega_n->0"
        ),
        "summatory_coefficients": (
            "c_omega(n)=n^omega*product_(p|n)(1-p^(-2*omega))>0"
        ),
        "summatory_function": (
            "h_omega^<1>(x)=x^(-1)*sum_(n<=x)c_omega(n)*"
            "g_omega^<1>(n/x)"
        ),
        "summatory_l2": (
            "Theta_omega inner iff "
            "x^(-1/2)*1_(1,infinity)(x)-h_omega^<1>(x) "
            "belongs to L2(1,infinity)"
        ),
        "summatory_sign": (
            "If h_omega^<1>(x) has one sign for all sufficiently large x, "
            "then Theta_omega is inner"
        ),
        "open_scalar_target": (
            "Prove the L2 residual, or the stronger eventual one-sign "
            "condition, for a sequence omega_n->0 by direct arithmetic "
            "bounds without assuming a zero-free half-plane"
        ),
        "growth_guard": (
            "No poles and real-boundary modulus one do not imply innerness "
            "for an arbitrary meromorphic function: exp(-i*z) is unbounded "
            "in C+; the xi high-strip bound is essential"
        ),
    }


def polynomial_checks() -> dict[str, str]:
    z = sp.symbols("z")
    i = sp.I

    f_chain = lambda w: (w**2 + 1) * (w**2 + 9)
    chain_ratio = sp.factor(sp.cancel(f_chain(z - i) / f_chain(z + i)))
    expected_chain = (z - 4 * i) / (z + 4 * i)
    if sp.simplify(chain_ratio - expected_chain) != 0:
        raise RuntimeError("fixed-shift cancellation chain did not factor")

    f_defect = lambda w: (w**2 + 1) * (w**2 + 9) ** 2
    defect_ratio = sp.factor(
        sp.cancel(f_defect(z - i) / f_defect(z + i))
    )
    expected_defect = (
        (z - 4 * i) ** 2
        * (z + 2 * i)
        / ((z - 2 * i) * (z + 4 * i) ** 2)
    )
    if sp.simplify(defect_ratio - expected_defect) != 0:
        raise RuntimeError("multiplicity-defect ratio did not factor")

    return {
        "chain_function": "F(z)=(z^2+1)*(z^2+9)",
        "chain_shift": "omega=1",
        "chain_ratio": "F(z-i)/F(z+i)=(z-4i)/(z+4i)",
        "chain_conclusion": (
            "The reduced ratio is inner although F has nonreal zeros; "
            "one fixed shift cannot force all zeros real."
        ),
        "defect_function": "F_def(z)=(z^2+1)*(z^2+9)^2",
        "defect_ratio": (
            "F_def(z-i)/F_def(z+i)="
            "(z-4i)^2*(z+2i)/((z-2i)*(z+4i)^2)"
        ),
        "defect_conclusion": (
            "The unmatched multiplicity leaves a pole at z=2i."
        ),
    }


def build_rows(exact: dict[str, str]) -> list[PhaseRow]:
    return [
        PhaseRow(
            "sfp_01_determinant_property",
            "definition",
            "available_exact",
            "The admissible fixed-shift determinant property is independent of the auxiliary power.",
            exact["determinant_property"],
            "Independence is proved only after the inner-function equivalence.",
        ),
        PhaseRow(
            "sfp_02_all_history_to_multiplier",
            "exact_lemma",
            "internally_audited",
            "All-time determinant nonvanishing produces a contractive causal multiplier.",
            exact["determinants_to_inner"],
            "Uses the separately audited dense-form extension and causal multiplier theorem.",
        ),
        PhaseRow(
            "sfp_03_boundary_symmetry",
            "exact_identity",
            "available_exact",
            "The xi functional equations make the quotient unimodular on the real boundary.",
            exact["boundary_unimodularity"],
            "Boundary values alone do not supply boundedness in C+.",
        ),
        PhaseRow(
            "sfp_04_power_independence",
            "exact_lemma",
            "available_exact",
            "Innerness of one positive integral power is equivalent to innerness of the base quotient.",
            exact["power_independence"],
            "A pole or modulus excess survives every positive power.",
        ),
        PhaseRow(
            "sfp_05_pole_coordinate",
            "exact_identity",
            "available_exact",
            "Every upper-half-plane pole is an off-line xi zero in explicit coordinates.",
            exact["pole_coordinate"],
            "Only zeros strictly to the right of 1/2+omega lie in C+.",
        ),
        PhaseRow(
            "sfp_06_multiplicity_cancellation",
            "exact_identity",
            "available_exact",
            "Pole removal is exactly a horizontal shifted-zero multiplicity inequality.",
            exact["multiplicity_rule"],
            "A common zero cancels only to the smaller multiplicity.",
        ),
        PhaseRow(
            "sfp_07_no_poles_to_inner",
            "published_analytic_step",
            "ready_to_apply",
            "For this xi quotient, absence of poles is sufficient for meromorphic innerness.",
            exact["no_poles_to_inner"],
            "This is Suzuki's fixed-omega Phragmen-Lindelof argument and is not a generic boundary principle.",
        ),
        PhaseRow(
            "sfp_08_inner_to_strict_compressions",
            "published_operator_step",
            "ready_to_apply",
            "Meromorphic innerness gives strict norm less than one for every compact truncation.",
            exact["inner_to_determinants"],
            "The proof uses Suzuki's noncompact-support lemma and nu*omega>1 compactness.",
        ),
        PhaseRow(
            "sfp_09_fixed_pair_equivalence",
            "theorem_candidate",
            "internally_audited",
            "The fixed-omega Fredholm problem is exactly a zero-pair cancellation problem.",
            exact["fixed_pair_equivalence"],
            "The combined equivalence should receive independent expert review.",
        ),
        PhaseRow(
            "sfp_10_cancellation_set",
            "exact_definition",
            "available_exact",
            "All exceptional shifts that can hide a right-hand zero form a countable set.",
            exact["cancellation_set"],
            "Countability does not make an explicitly chosen omega provably generic.",
        ),
        PhaseRow(
            "sfp_11_generic_shift",
            "exact_corollary",
            "internally_audited",
            "Outside the cancellation set, D(omega) is equivalent to the zero-free half-plane Re(s)>1/2+omega.",
            "omega not in C_xi => [D(omega) iff xi(s)!=0 for Re(s)>1/2+omega]",
            "This is a generic-parameter equivalence, not a zero-free proof.",
        ),
        PhaseRow(
            "sfp_12_zero_width",
            "exact_definition",
            "available_exact",
            "The maximal horizontal zero displacement is the determinant phase threshold.",
            exact["zero_width"],
            "The value delta_xi=0 is exactly RH.",
        ),
        PhaseRow(
            "sfp_13_phase_diagram",
            "exact_corollary",
            "internally_audited",
            "Apart from countably many possible cancellation shifts below threshold, D is a sharp half-line.",
            exact["phase_diagram"],
            "The exceptional points need not themselves satisfy D.",
        ),
        PhaseRow(
            "sfp_14_local_failure_if_not_rh",
            "exact_corollary",
            "internally_audited",
            "If RH fails, all sufficiently small shifts fail the determinant property.",
            exact["local_exclusion"],
            "The exclusion radius depends on an unknown off-line zero.",
        ),
        PhaseRow(
            "sfp_15_closure_criterion",
            "exact_corollary",
            "internally_audited",
            "Cofinality can be weakened to the topological statement that zero is a limit point of successful shifts.",
            exact["closure_criterion"],
            "This remains an RH-equivalent criterion, not evidence that such shifts exist.",
        ),
        PhaseRow(
            "sfp_16_chain_countermodel",
            "exact_countermodel",
            "guard_validated",
            "A finite vertical zero chain makes one fixed quotient inner despite nonreal zeros.",
            "F(z)=(z^2+1)*(z^2+9), omega=1 => Theta_F(z)=(z-4i)/(z+4i)",
            "This tests fixed-shift logic only; its rational quotient is not Suzuki's decaying xi kernel.",
        ),
        PhaseRow(
            "sfp_17_multiplicity_countermodel",
            "exact_countermodel",
            "guard_validated",
            "Insufficient lower multiplicity leaves an uncancelled upper-half-plane pole.",
            "F_def(z)=(z^2+1)*(z^2+9)^2 => pole at z=2i",
            "This validates the direction and multiplicity in the zero-shift rule.",
        ),
        PhaseRow(
            "sfp_18_growth_guard",
            "exact_countermodel",
            "guard_validated",
            "Boundary unimodularity and pole-freeness are not enough without the xi growth input.",
            exact["growth_guard"],
            "Suzuki's high-strip estimate and Phragmen-Lindelof step must remain explicit.",
        ),
        PhaseRow(
            "sfp_19_summatory_l2",
            "published_scalar_equivalence",
            "ready_to_apply",
            "Suzuki's Jordan-totient weighted residual is an exact scalar test for the same innerness gate.",
            exact["summatory_l2"],
            "The residual is global on (1,infinity); finite samples do not certify it.",
        ),
        PhaseRow(
            "sfp_20_summatory_sign",
            "published_sufficient_condition",
            "ready_to_apply",
            "Eventual one-sign behavior of the weighted summatory function is a stronger scalar route.",
            exact["summatory_sign"],
            "Positive coefficients alone do not decide the signed weight transform.",
        ),
        PhaseRow(
            "sfp_21_open_scalar_target",
            "open_arithmetic_gate",
            "not_ready_to_apply",
            "The operator gate can be attacked through one explicit arithmetic summatory function.",
            exact["open_scalar_target"],
            "This target is RH-strength on a cofinal sequence and must not be inferred from finite positivity.",
        ),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    checks = payload["polynomial_checks"]
    return "\n".join(
        [
            "# Jensen-Window PF Suzuki Fixed-Omega Phase Diagram",
            "",
            "Date: 2026-07-23",
            "",
            "Status: internally audited theorem-candidate note with exact",
            "countermodel guards and one open arithmetic target. It is not a proof",
            "of the determinant premise, RH, or `Lambda <= 0`.",
            "Independent expert review is required before publication.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_suzuki_fixed_omega_phase_diagram.json",
            "python work/rh_compute/scripts/jensen_window_pf_suzuki_fixed_omega_phase_diagram.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_suzuki_fixed_omega_phase_diagram.py",
            "```",
            "",
            "## Fixed-Omega Equivalence",
            "",
            "For `omega>0`, let",
            "",
            "```text",
            exact["determinant_property"],
            "Theta_omega(z)=xi(1/2-omega-i*z)/xi(1/2+omega-i*z).",
            "```",
            "",
            "Combining the audited causal-multiplier direction with Suzuki's",
            "fixed-parameter Phragmen-Lindelof and strict-compression arguments gives",
            "",
            "```text",
            exact["fixed_pair_equivalence"],
            "```",
            "",
            "This combined equivalence is a corpus theorem candidate. Its directions",
            "are kept explicit below.",
            "",
            "## Determinants To Innerness",
            "",
            "All-time determinant nonvanishing is equivalent to strict contraction",
            "of every compact self-adjoint truncation. Compatible compactly supported",
            "forms therefore extend to a bounded full-line Hankel form. Reflection",
            "turns it into a bounded causal translation-invariant convolution",
            "operator, whose multiplier is `Theta_omega^nu` in the high half-plane.",
            "Hence",
            "",
            "```text",
            exact["determinants_to_inner"],
            exact["power_independence"],
            "```",
            "",
            "Real-boundary unimodularity is exact:",
            "",
            "```text",
            exact["boundary_unimodularity"],
            "```",
            "",
            "## Innerness To Determinants",
            "",
            "Suzuki's earlier fixed-omega argument observes that the quotient is",
            "uniformly bounded in a sufficiently high strip. If it has no poles in",
            "`C+`, Phragmen-Lindelof propagates that bound down to the real boundary.",
            "Thus",
            "",
            "```text",
            exact["no_poles_to_inner"],
            "```",
            "",
            "Once the quotient is inner, its boundary multiplier is an `L2` isometry.",
            "Suzuki's noncompact-support lemma excludes equality after finite",
            "compression, while `nu*omega>1` makes each truncation compact. Therefore",
            "",
            "```text",
            exact["inner_to_determinants"],
            "```",
            "",
            "The growth input is indispensable: `exp(-i*z)` has no poles and has",
            "modulus one on the real line, but grows like `exp(Im(z))` in `C+`.",
            "",
            "## Exact Zero-Pair Rule",
            "",
            "For `rho=beta+i*gamma`,",
            "",
            "```text",
            exact["pole_coordinate"],
            exact["multiplicity_rule"],
            "```",
            "",
            "Consequently a fixed shift can hide an off-line zero only through an",
            "equal-height horizontal zero pair separated by `2*omega`.",
            "",
            "## Phase Diagram",
            "",
            "Define",
            "",
            "```text",
            exact["cancellation_set"],
            exact["zero_width"],
            "```",
            "",
            "Then",
            "",
            "```text",
            exact["phase_diagram"],
            "omega not in C_xi => [D(omega) iff xi(s)!=0 for Re(s)>1/2+omega].",
            "```",
            "",
            "The possible subthreshold successes are therefore countable and",
            "entirely attributable to exact zero cancellation. More strongly,",
            "isolation of any one off-line zero gives",
            "",
            "```text",
            exact["local_exclusion"],
            exact["closure_criterion"],
            "```",
            "",
            "## Fixed-Shift Countermodels",
            "",
            "The cancellation guard is realized exactly by",
            "",
            "```text",
            checks["chain_function"],
            checks["chain_ratio"],
            "```",
            "",
            "The quotient is a Blaschke factor although `F` has zeros at",
            "`+/-i` and `+/-3i`. Raising the upper pair's multiplicity breaks the",
            "condition:",
            "",
            "```text",
            checks["defect_function"],
            checks["defect_ratio"],
            checks["defect_conclusion"],
            "```",
            "",
            "These polynomial examples test the cancellation logic. They are not",
            "surrogates for Suzuki's decaying arithmetic kernel.",
            "",
            "## Jordan-Totient Scalar Target",
            "",
            "Suzuki also gives a scalar arithmetic formulation. Put",
            "",
            "```text",
            exact["summatory_coefficients"],
            "g_omega^<1>(x)=integral_x^1 sqrt(y/x)*g_omega(y)dy/y,",
            exact["summatory_function"],
            "```",
            "",
            "where `g_omega` is Suzuki's explicit beta-integral weight. Then",
            "",
            "```text",
            exact["summatory_l2"],
            exact["summatory_sign"],
            "```",
            "",
            "Finite samples do not certify it: both conditions concern the full",
            "unbounded half-line.",
            "",
            "The sharp scalar obligation is therefore",
            "",
            "```text",
            exact["open_scalar_target"],
            "```",
            "",
            "Although every `c_omega(n)` is positive, the weight is signed. For",
            "`omega=1/2`, its explicit primitive is negative near zero and positive",
            "elsewhere, so coefficient positivity alone cannot establish the target.",
            "",
            "## Sources",
            "",
            "- Masatoshi Suzuki, `Hamiltonians arising from L-functions in the Selberg class`: https://arxiv.org/abs/1606.05726",
            "- Masatoshi Suzuki, `A canonical system of differential equations arising from the Riemann zeta-function`: https://arxiv.org/abs/1204.1827",
            "- Masatoshi Suzuki, `On monotonicity of certain weighted summatory functions associated with L-functions`: https://arxiv.org/abs/1204.1823",
            "- Chris Guiver, Hartmut Logemann, and Mark R. Opmeer, causal `L2` multiplier theorem: https://link.springer.com/article/10.1007/s00498-024-00387-4",
            "",
            "## Proof Boundary",
            "",
            payload["proof_boundary"],
            "",
        ]
    )


def build_payload() -> dict:
    exact = exact_statements()
    checks = polynomial_checks()
    rows = build_rows(exact)
    return {
        "kind": "jensen_window_pf_suzuki_fixed_omega_phase_diagram",
        "date": "2026-07-23",
        "status": (
            "internally audited fixed-omega equivalence and phase-diagram "
            "theorem candidate"
        ),
        "proof_boundary": (
            "This artifact identifies the fixed-omega Suzuki determinant "
            "property with meromorphic innerness and an exact shifted-zero "
            "multiplicity condition, subject to independent review of the "
            "combined operator argument. It does not prove that the condition "
            "holds at any subcritical omega, does not prove eventual sign or "
            "the L2 summatory residual, and does not prove RH or Lambda<=0."
        ),
        "exact": exact,
        "polynomial_checks": checks,
        "rows": [asdict(row) for row in rows],
        "sources": [
            {
                "title": "Hamiltonians arising from L-functions in the Selberg class",
                "url": "https://arxiv.org/abs/1606.05726",
                "use": "kernel construction, noncompact-support lemma, strict finite compression",
            },
            {
                "title": "A canonical system of differential equations arising from the Riemann zeta-function",
                "url": "https://arxiv.org/abs/1204.1827",
                "use": "fixed-omega no-pole Phragmen-Lindelof step and summatory criterion",
            },
            {
                "title": "On monotonicity of certain weighted summatory functions associated with L-functions",
                "url": "https://arxiv.org/abs/1204.1823",
                "use": "Jordan-totient weighted sign and L2 criteria",
            },
            {
                "title": "Operator-valued multiplier theorems for causal translation-invariant operators",
                "url": "https://link.springer.com/article/10.1007/s00498-024-00387-4",
                "use": "causal L2 multiplier equivalence",
            },
        ],
        "audit": {
            "row_count": len(rows),
            "exact_countermodel_count": 3,
            "published_scalar_target_count": 2,
            "open_arithmetic_gate_count": 1,
            "rh_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Suzuki fixed-omega phase diagram: "
        f"{len(payload['rows'])} rows, "
        "3 exact countermodel guards, 2 published scalar targets, "
        "1 open arithmetic gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
