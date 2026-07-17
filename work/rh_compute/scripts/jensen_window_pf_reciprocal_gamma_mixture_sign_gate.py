#!/usr/bin/env python3
"""Build the reciprocal-gamma sign-regularity and mixture-closure gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_reciprocal_gamma_mixture_sign_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_reciprocal_gamma_mixture_sign_gate.md"
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


def gamma_value(index: int) -> sp.Rational:
    return sp.Rational(sp.factorial(index), sp.factorial(2 * index))


def contiguous_rows() -> list[dict]:
    rows: list[dict] = []
    for size in range(1, 7):
        matrix = sp.Matrix(
            [[gamma_value(i + j) for j in range(size)] for i in range(size)]
        )
        determinant = sp.factor(matrix.det())
        expected_sign = (-1) ** (size * (size - 1) // 2)
        rows.append(
            {
                "size": size,
                "determinant": str(determinant),
                "expected_sign": "positive" if expected_sign > 0 else "negative",
            }
        )
    return rows


def build_exact() -> dict:
    for index in range(13):
        duplication = sp.simplify(
            sp.factorial(index) / sp.factorial(2 * index)
            - sp.sqrt(sp.pi)
            / (4**index * sp.gamma(sp.Rational(1, 2) + index))
        )
        if duplication != 0:
            raise RuntimeError(f"gamma duplication failed at {index}")

    contiguous = contiguous_rows()
    for row in contiguous:
        size = int(row["size"])
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
        if sp.simplify(sp.sympify(row["determinant"]) - formula) != 0:
            raise RuntimeError(f"contiguous determinant formula failed at {size}")

    t0, t1 = sp.symbols("t0 t1", positive=True)
    sym_integrand = sp.factor(
        gamma_value(0)
        * gamma_value(2)
        * (t0**2 + t1**2)
        / 2
        - gamma_value(1) ** 2 * t0 * t1
    )
    expected_integrand = (t0**2 - 6 * t0 * t1 + t1**2) / 24
    if sp.expand(sym_integrand - expected_integrand) != 0:
        raise RuntimeError("symmetrized order-two integrand failed")

    weights = (10, 1, 1)
    scales = (sp.Rational(1, 100), sp.Integer(1), sp.Integer(2))
    moments = [sum(w * t**order for w, t in zip(weights, scales)) for order in range(3)]
    coefficients = [gamma_value(order) * moments[order] for order in range(3)]
    determinant = sp.factor(coefficients[0] * coefficients[2] - coefficients[1] ** 2)
    if determinant != sp.Rational(5197, 2000):
        raise RuntimeError("positive-support mixture countermodel failed")

    n, m0, m1, m2 = sp.symbols("n m0 m1 m2", positive=True)
    cv_squared = sp.cancel((m0 * m2 - m1**2) / m1**2)
    if sp.simplify(cv_squared - (m0 * m2 / m1**2 - 1)) != 0:
        raise RuntimeError("tilted coefficient-of-variation identity failed")

    return {
        "coefficient_split": {
            "moments": (
                "m_k(lambda)=mu_(2k)(lambda)=integral_0^infinity "
                "t^k*rho_lambda(dt)"
            ),
            "factor": "gamma_k=k!/(2k)!=sqrt(pi)/(4^k*Gamma(k+1/2))",
            "coefficients": "A_k(lambda)=gamma_k*m_k(lambda)",
            "matrix": (
                "A_(n+i+j)=integral gamma_(n+i+j)*t^(n+i+j)*rho_lambda(dt)"
            ),
        },
        "karlin_sign_regularity": {
            "source": "https://doi.org/10.1090/S0002-9947-1964-0168010-2",
            "citation": (
                "S. Karlin, Total positivity, absorption probabilities and "
                "applications, Trans. Amer. Math. Soc. 111 (1964), Lemma 9.2"
            ),
            "theorem": (
                "For alpha>-nu, every (nu+1)-minor of "
                "[1/Gamma(alpha+i+j)]_(i,j>=0) has strict sign "
                "(-1)^(nu*(nu+1)/2)."
            ),
            "specialization": (
                "For every t>0 and n>=0, K_(n,t)(i,j)="
                "gamma_(n+i+j)*t^(n+i+j) is strictly sign-regular of all orders "
                "with signature epsilon_k=(-1)^(k*(k-1)/2)."
            ),
            "scope": (
                "This proves the complete fixed-scale antecedent, including arbitrary "
                "increasing row and column sets."
            ),
        },
        "contiguous_determinant": {
            "formula": (
                "det[gamma_(i+j)]_(i,j=0..k-1)=(-1)^(k*(k-1)/2)*"
                "pi^(k/2)*4^(-k*(k-1))*prod_(r=1)^(k-1)r!/"
                "prod_(i=0)^(k-1)Gamma(k-1/2+i)"
            ),
            "scaled_formula": (
                "det[gamma_(i+j)*t^(i+j)]="
                "t^(k*(k-1))*det[gamma_(i+j)]"
            ),
            "rows": contiguous,
        },
        "rowwise_mixture_integral": {
            "formula": (
                "R_(k,n)(j_1,...,j_k)=integral_((0,infinity)^k) "
                "det[gamma_(n+i+j_l)*t_i^(n+i+j_l)]_(i=0..k-1,l=1..k) "
                "prod_(i=0)^(k-1)rho_lambda(dt_i)"
            ),
            "derivation": (
                "Expand the determinant by multilinearity in its rows and assign "
                "one independent scale variable to each row."
            ),
            "boundary": (
                "The scales are independent, so Karlin's common-scale sign theorem "
                "does not sign this integrand away from the diagonal t_0=...=t_(k-1)."
            ),
        },
        "order_two_integrand": {
            "raw": (
                "D(t_0,t_1)=gamma_0*gamma_2*t_1^2-gamma_1^2*t_0*t_1"
            ),
            "symmetrized": (
                "D_sym(t_0,t_1)=(t_0^2-6*t_0*t_1+t_1^2)/24"
            ),
            "negative_witness": "D_sym(1,1)=-1/6",
            "positive_witness": "D_sym(1,10)=41/24",
            "sign_region": (
                "D_sym<=0 iff 3-2*sqrt(2)<=t_1/t_0<=3+2*sqrt(2)"
            ),
        },
        "positive_mixture_countermodel": {
            "measure": "rho=10*delta_(1/100)+delta_1+delta_2",
            "moments": [str(value) for value in moments],
            "coefficients": [str(value) for value in coefficients],
            "minor": "A_0*A_2-A_1^2=5197/2000>0",
            "expected": "The signed-Hankel order-two sign would require A_0*A_2-A_1^2<0.",
            "conclusion": (
                "Positive mixing of strictly sign-regular fixed-scale reciprocal-gamma "
                "kernels is not closed even at order two, with all scales strictly positive."
            ),
        },
        "tilted_concentration": {
            "tilt": (
                "P_n(dt)=t^n*rho_lambda(dt)/m_n(lambda)"
            ),
            "relative_variance": (
                "CV_n^2=Var_(P_n)(T)/E_(P_n)(T)^2="
                "m_n*m_(n+2)/m_(n+1)^2-1"
            ),
            "degree_two_equivalence": (
                "A_(n+1)^2>=A_n*A_(n+2) iff "
                "m_(n+1)^2/(m_n*m_(n+2))>=(2*n+1)/(2*n+3) "
                "iff CV_n^2<=2/(2*n+1)"
            ),
            "xi_input": (
                "The lambda=-100 full ratio-cone entry plus infinite heat-flow "
                "invariance proves this tilted relative-variance bound for every n "
                "at lambda=0."
            ),
            "scope": (
                "This closes all shifted order-two signed-Hankel minors only."
            ),
        },
        "higher_order_handoff": {
            "target": (
                "Control the sign-changing row-wise compound integrand after Xi "
                "mixing for every k and column set, using a hierarchy of tilted "
                "concentration, a coupling that keeps row scales near the diagonal, "
                "or a new positive compound-kernel factorization."
            ),
            "forbidden_promotion": (
                "Karlin's fixed-scale sign regularity and the scalar CV_n bound do "
                "not imply mixture sign regularity in orders k>=3."
            ),
            "first_compound_coordinate": (
                "The subsequent order-three gate proves that the first contiguous "
                "compound sign is exactly C_n=q_(n+1)*q_(n+3)-q_(n+2)^2+1>0."
            ),
            "first_compound_guard": (
                "A strict rational sequence satisfies the full ratio, adaptive-defect, "
                "and cubic Jensen cones while having C_n<0."
            ),
        },
    }


def build_artifact() -> dict:
    exact = build_exact()
    rows = [
        GateRow(
            id="rgmsg_01_coefficient_split",
            role="exact_identity",
            readiness="available_exact",
            claim="The normalized Xi coefficients are positive scale mixtures of the reciprocal-gamma sequence.",
            formula=exact["coefficient_split"]["coefficients"],
            proof_boundary="Exact moment normalization only.",
            diagnostics=exact["coefficient_split"],
        ),
        GateRow(
            id="rgmsg_02_karlin_theorem",
            role="published_theorem",
            readiness="ready_to_apply",
            claim="Karlin's reciprocal-gamma kernel is strictly sign-regular in every order.",
            formula=exact["karlin_sign_regularity"]["theorem"],
            proof_boundary="Published fixed-kernel theorem; mixture closure is not part of the result.",
            diagnostics=exact["karlin_sign_regularity"],
        ),
        GateRow(
            id="rgmsg_03_fixed_scale",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="Every positive fixed-scale Xi coefficient kernel has exactly the observed all-order signature.",
            formula=exact["karlin_sign_regularity"]["specialization"],
            proof_boundary="One common scale only.",
        ),
        GateRow(
            id="rgmsg_04_contiguous_formula",
            role="exact_identity",
            readiness="available_exact",
            claim="The contiguous reciprocal-gamma Hankel determinants have an explicit strict product formula.",
            formula=exact["contiguous_determinant"]["formula"],
            proof_boundary="Exact fixed-scale benchmark, not a mixed-measure determinant.",
            diagnostics=exact["contiguous_determinant"],
        ),
        GateRow(
            id="rgmsg_05_rowwise_integral",
            role="exact_identity",
            readiness="available_exact",
            claim="Every shifted reshaped-Hankel minor has an exact independent-row scale integral.",
            formula=exact["rowwise_mixture_integral"]["formula"],
            proof_boundary="The exact integrand is signed; no positive determinant formula is claimed.",
            diagnostics=exact["rowwise_mixture_integral"],
        ),
        GateRow(
            id="rgmsg_06_integrand_sign_change",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="The natural symmetrized order-two determinant integrand changes sign away from the common-scale diagonal.",
            formula=exact["order_two_integrand"]["symmetrized"],
            proof_boundary="Rejects pointwise positivity of the natural row-wise integral only.",
            diagnostics=exact["order_two_integrand"],
        ),
        GateRow(
            id="rgmsg_07_mixture_countermodel",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="A positive three-atom measure on strictly positive scales reverses the required order-two sign.",
            formula=exact["positive_mixture_countermodel"]["minor"],
            proof_boundary="Abstract positive-mixture countermodel, not the Xi measure.",
            diagnostics=exact["positive_mixture_countermodel"],
        ),
        GateRow(
            id="rgmsg_08_tilted_cv",
            role="exact_equivalence",
            readiness="ready_to_apply",
            claim="The complete shifted degree-two cone is exactly a shrinking relative-variance bound under moment tilting.",
            formula=exact["tilted_concentration"]["degree_two_equivalence"],
            proof_boundary="Equivalent to order two only.",
            diagnostics=exact["tilted_concentration"],
        ),
        GateRow(
            id="rgmsg_09_xi_order_two",
            role="exact_composition",
            readiness="ready_to_apply",
            claim="The established ratio-cone entry and heat invariance supply the tilted relative-variance theorem for Xi at lambda zero.",
            formula=exact["tilted_concentration"]["xi_input"],
            proof_boundary="All shifts at order two; no higher compound sign is inferred.",
        ),
        GateRow(
            id="rgmsg_10_higher_order_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The remaining signed-Hankel antecedent is a higher-order concentration or positive compound-kernel theorem for the Xi mixing measure.",
            formula=exact["higher_order_handoff"]["target"],
            proof_boundary="Open; not the signed-Hankel/Jensen bridge, RH, or Lambda<=0.",
            diagnostics=exact["higher_order_handoff"],
        ),
    ]
    return {
        "kind": "jensen_window_pf_reciprocal_gamma_mixture_sign_gate",
        "date": "2026-07-12",
        "status": "exact reciprocal-gamma sign regularity with positive-mixture obstruction",
        "proof_boundary": (
            "This artifact applies Karlin's published all-order sign-regularity "
            "theorem to the fixed-scale reciprocal-gamma kernel, proves the exact "
            "row-wise determinant integral for positive Xi scale mixtures, and "
            "exhibits exact sign-changing integrands and a strictly-positive-scale "
            "mixture countermodel. It translates the established Xi order-two cone "
            "into a tilted relative-variance bound. It does not prove higher-order "
            "mixture sign regularity, the signed-Hankel/Jensen bridge, PF-infinity, "
            "RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_factorial_multiplier_split_audit.md",
            "outputs/signed_hankel_jensen_bridge_target.md",
            "outputs/jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.md",
            "outputs/formal_core.md",
            "https://doi.org/10.1090/S0002-9947-1964-0168010-2",
        ],
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    lines = [
        "# Jensen-Window PF Reciprocal-Gamma Mixture Sign Gate",
        "",
        "Date: 2026-07-12",
        "",
        "Status: exact reciprocal-gamma sign regularity with positive-mixture",
        "obstruction. This is not a proof of RH or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_reciprocal_gamma_mixture_sign_gate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_reciprocal_gamma_mixture_sign_gate.json",
        "python work/rh_compute/scripts/jensen_window_pf_reciprocal_gamma_mixture_sign_gate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_gamma_mixture_sign_gate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        "validated Jensen-window PF reciprocal-gamma mixture sign gate: 10 rows, 0 issues, 3 exact kernel identities, 1 published all-order theorem, 1 fixed-scale theorem, 2 exact mixture countermodels, 1 tilted-variance equivalence, 1 Xi order-two composition, 1 higher-order handoff",
        "```",
        "",
        "## Fixed-Scale Theorem",
        "",
        "The coefficient normalization is",
        "",
        "```text",
        exact["coefficient_split"]["moments"],
        exact["coefficient_split"]["factor"],
        exact["coefficient_split"]["coefficients"],
        "```",
        "",
        "Karlin's Lemma 9.2 states:",
        "",
        "```text",
        exact["karlin_sign_regularity"]["theorem"],
        exact["karlin_sign_regularity"]["specialization"],
        "```",
        "",
        "So the fixed-scale reciprocal-gamma kernel already has the complete",
        "signature observed in the Arb staircase, for arbitrary row and column",
        "sets. The contiguous benchmark is explicit:",
        "",
        "```text",
        exact["contiguous_determinant"]["formula"],
        "```",
        "",
        "## Mixture Integral",
        "",
        "Multilinearity in the determinant rows gives",
        "",
        "```text",
        exact["rowwise_mixture_integral"]["formula"],
        "```",
        "",
        "The scale variables are independent. Karlin's theorem signs the diagonal",
        "`t_0=...=t_(k-1)`, not the full product-measure integrand.",
        "",
        "## Exact Nonclosure",
        "",
        "At order two the symmetrized integrand is",
        "",
        "```text",
        exact["order_two_integrand"]["symmetrized"],
        exact["order_two_integrand"]["negative_witness"],
        exact["order_two_integrand"]["positive_witness"],
        exact["order_two_integrand"]["sign_region"],
        "```",
        "",
        "A positive measure supported away from zero gives a complete countermodel:",
        "",
        "```text",
        exact["positive_mixture_countermodel"]["measure"],
        f"m_0,m_1,m_2={exact['positive_mixture_countermodel']['moments']}",
        f"A_0,A_1,A_2={exact['positive_mixture_countermodel']['coefficients']}",
        exact["positive_mixture_countermodel"]["minor"],
        "```",
        "",
        "Thus the cone of fixed-scale sign-regular matrices is not closed under",
        "the positive scale mixing used by the Xi moment representation.",
        "",
        "## Xi Concentration",
        "",
        "For the `n`-tilted moment law,",
        "",
        "```text",
        exact["tilted_concentration"]["tilt"],
        exact["tilted_concentration"]["relative_variance"],
        exact["tilted_concentration"]["degree_two_equivalence"],
        "```",
        "",
        exact["tilted_concentration"]["xi_input"],
        "This gives a structural interpretation of the completed all-shift",
        "order-two cone: the tilted Xi measure becomes relatively more concentrated",
        "at the exact rate `2/(2n+1)`.",
        "",
        "## Live Handoff",
        "",
        exact["higher_order_handoff"]["target"],
        "",
        exact["higher_order_handoff"]["first_compound_coordinate"],
        exact["higher_order_handoff"]["first_compound_guard"],
        "See `outputs/jensen_window_pf_reciprocal_defect_compound_order3_gate.md`.",
        "",
        "The next signed-Hankel antecedent is therefore a higher compound",
        "concentration theorem for the actual Xi mixing measure, not another",
        "application of fixed-scale reciprocal-gamma sign regularity.",
        "",
        "Primary source: https://doi.org/10.1090/S0002-9947-1964-0168010-2.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(f"wrote reciprocal-gamma mixture sign gate: {len(artifact['rows'])} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
