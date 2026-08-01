#!/usr/bin/env python3
"""Build the exact score-Abel radial dimension-lift gate."""

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
    "jensen_window_pf_newman_score_abel_radial_dimension_lift_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_score_abel_radial_dimension_lift_gate.md"
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


def build_exact() -> dict:
    coefficient_checks: list[str] = []
    normalization_checks: list[str] = []
    characteristic_checks: list[str] = []
    gaussian_moment_checks: list[str] = []
    gaussian_discriminant_checks: list[str] = []
    xi_quadratic_scaling_checks: list[str] = []

    # The nth derivative of I_0(2*sqrt(v*z)) is the 0F1 kernel.
    for shift in range(7):
        for power in range(8):
            direct = 1 / (
                sp.factorial(shift + power) * sp.factorial(power)
            )
            hyper = 1 / (
                sp.factorial(shift)
                * sp.rf(shift + 1, power)
                * sp.factorial(power)
            )
            if sp.simplify(direct - hyper) != 0:
                raise RuntimeError(
                    "0F1 derivative coefficient identity failed "
                    f"at n={shift}, k={power}"
                )
            coefficient_checks.append(f"n={shift},k={power}")

    # Polar integration in dimension 2n+2 normalizes the proposed density.
    for shift in range(9):
        sphere_area = (
            2 * sp.pi ** (shift + 1) / sp.factorial(shift)
        )
        radial_jacobian = 2 ** (2 * shift + 1)
        density_constant = (4 * sp.pi) ** (shift + 1)
        coefficient = sp.simplify(
            sphere_area * radial_jacobian / density_constant
        )
        if coefficient != 1 / sp.factorial(shift):
            raise RuntimeError(
                f"radial normalization failed at n={shift}"
            )
        normalization_checks.append(f"n={shift}")

    # Coefficients of the radial characteristic function match F^(n)/A_n.
    for shift in range(7):
        for power in range(7):
            a_n = sp.symbols(f"A_{shift}", positive=True)
            a_nk = sp.symbols(
                f"A_{shift + power}", positive=True
            )
            moment_ratio = (
                sp.factorial(shift + power)
                * a_nk
                / (sp.factorial(shift) * a_n)
            )
            radial_coefficient = sp.simplify(
                moment_ratio
                / (
                    sp.rf(shift + 1, power)
                    * sp.factorial(power)
                )
            )
            derivative_coefficient = (
                a_nk / (a_n * sp.factorial(power))
            )
            if (
                sp.simplify(
                    radial_coefficient - derivative_coefficient
                )
                != 0
            ):
                raise RuntimeError(
                    "radial characteristic coefficient failed "
                    f"at n={shift}, k={power}"
                )
            characteristic_checks.append(
                f"n={shift},k={power}"
            )

    x = sp.symbols("x", nonzero=True, real=True)
    a_n, a_np1, f_np1 = sp.symbols(
        "A_n A_np1 F_np1", positive=True
    )
    phi_prime = -2 * x * f_np1 / a_n
    dimension_walk_rhs = sp.simplify(
        -a_n * phi_prime / (2 * x * a_np1)
    )
    if sp.simplify(dimension_walk_rhs - f_np1 / a_np1) != 0:
        raise RuntimeError("radial dimension-walk identity failed")

    abel_inner = sp.expand_func(
        sp.beta(sp.Rational(1, 2), sp.Rational(1, 2))
    ) / 2
    if sp.simplify(abel_inner - sp.pi / 2) != 0:
        raise RuntimeError("Abel marginal beta integral failed")

    # The two-Gaussian model is an exact Newman flow and gives a strict
    # quadratic Jensen failure for every shift.
    p, q, alpha, beta = sp.symbols(
        "p q alpha beta", positive=True
    )
    for shift in range(9):
        current = p * alpha**shift + q * beta**shift
        first = p * alpha ** (shift + 1) + q * beta ** (shift + 1)
        second = (
            p * alpha ** (shift + 2)
            + q * beta ** (shift + 2)
        )
        determinant = sp.factor(first**2 - current * second)
        expected = -p * q * (alpha * beta) ** shift * (
            alpha - beta
        ) ** 2
        if sp.simplify(determinant - expected) != 0:
            raise RuntimeError(
                "two-Gaussian Jensen discriminant failed "
                f"at n={shift}"
            )
        gaussian_discriminant_checks.append(f"n={shift}")

    a = sp.symbols("a", positive=True)
    for power in range(9):
        abel_moment = sp.simplify(
            4
            * sp.sqrt(sp.pi * a)
            * sp.factorial(power)
            / (4 * a) ** (power + 1)
        )
        generating_moment = (
            sp.factorial(power)
            * sp.sqrt(sp.pi / a)
            / (4 * a) ** power
        )
        if sp.simplify(abel_moment - generating_moment) != 0:
            raise RuntimeError(
                f"Gaussian Abel moment failed at k={power}"
            )
        gaussian_moment_checks.append(f"k={power}")

    c, t = sp.symbols("c t", positive=True)
    a_ct = c**2 - t
    b_ct = 2 * c**2 - t
    delta = sp.simplify(1 / (4 * a_ct) - 1 / (4 * b_ct))
    if sp.simplify(delta - c**2 / (4 * a_ct * b_ct)) != 0:
        raise RuntimeError("two-Gaussian exponential gap failed")
    k = sp.symbols("k", integer=True)
    zero = (
        2 * a_ct * b_ct * sp.log(a_ct / b_ct) / c**2
        + 4
        * sp.pi
        * sp.I
        * a_ct
        * b_ct
        * (2 * k + 1)
        / c**2
    )
    logarithmic_equation = (
        sp.log(a_ct / b_ct) / 2
        + sp.I * sp.pi * (2 * k + 1)
    )
    if sp.simplify(delta * zero - logarithmic_equation) != 0:
        raise RuntimeError("two-Gaussian zero lattice failed")

    strong_constant = (2 - 4 / sp.E) * c**2 - sp.Rational(2, 5)
    if not bool(sp.N(strong_constant.subs(c, 1), 30) > 0):
        raise RuntimeError("two-Gaussian strong-concavity bound failed")

    # Translate the already-proved log-concavity of A_k into the Abel
    # moment concentration inequality without losing the factorial scale.
    for shift in range(12):
        m0, m1, m2 = sp.symbols(
            f"M_{shift} M_{shift + 1} M_{shift + 2}",
            positive=True,
        )
        a0 = m0 / sp.factorial(shift)
        a1 = m1 / sp.factorial(shift + 1)
        a2 = m2 / sp.factorial(shift + 2)
        scaled = sp.simplify(
            sp.factorial(shift + 1) ** 2
            * (a1**2 - a0 * a2)
        )
        expected = sp.simplify(
            m1**2
            - sp.Rational(shift + 1, shift + 2) * m0 * m2
        )
        if sp.simplify(scaled - expected) != 0:
            raise RuntimeError(
                "Xi quadratic Abel scaling failed "
                f"at n={shift}"
            )
        xi_quadratic_scaling_checks.append(f"n={shift}")

    return {
        "generating_and_moments": (
            "Let mathcal_F_t(z)=2*H_t(sqrt(-z))="
            "sum_(k>=0) A_k(t)*z^k/k!. For the established Abel "
            "density r_t, M_k=integral_0^infinity v^k*r_t(v)dv="
            "k!*A_k."
        ),
        "radial_density": (
            "For n>=0 and y in R^(2n+2), define "
            "g_(n,t)(y)=r_t(|y|^2/4)/((4*pi)^(n+1)*A_n(t)). "
            "Then g_(n,t)>=0 and integral_R^(2n+2) g_(n,t)(y)dy=1."
        ),
        "radial_mixing_law": (
            "Under g_(n,t), V=|Y|^2/4 has probability law "
            "dmu_(n,t)(v)=v^n*r_t(v)dv/M_n(t). Thus each shift n "
            "uses the n-th size bias of the same Abel measure."
        ),
        "abel_marginal": (
            "For n=0, g_(0,t) is the isotropic planar Abel/Radon "
            "lift of f_t/A_0: integral_R g_(0,t)(u,y)dy=f_t(u)/A_0. "
            "A direct proof swaps the Abel integrals and uses "
            "integral_u^s r*dr/(sqrt(r^2-u^2)*sqrt(s^2-r^2))=pi/2."
        ),
        "derivative_bessel": (
            "mathcal_F_t^(n)(z)=1/n!*integral_0^infinity "
            "v^n*0F1(;n+1;v*z)*r_t(v)dv, and "
            "0F1(;n+1;-v*x^2)=n!*(x*sqrt(v))^(-n)"
            "*J_n(2*x*sqrt(v)), with the continuous value at x=0."
        ),
        "radial_characteristic": (
            "With the Fourier convention exp(i*xi dot y), "
            "hat(g_(n,t))(xi)=mathcal_F_t^(n)(-|xi|^2)/A_n(t). "
            "Consequently every normalized derivative is radial "
            "positive definite in its own dimension 2n+2."
        ),
        "dimension_walk": (
            "If phi_(n,t)(x)=mathcal_F_t^(n)(-x^2)/A_n(t), then "
            "phi_(n+1,t)(x)=-A_n(t)*phi_(n,t)'(x)/"
            "(2*x*A_(n+1)(t)), x>0, with continuous extension at 0."
        ),
        "jensen_window": (
            "P_(D,n)(w)=sum_(k=0)^D binom(D,k)A_(n+k)w^k "
            "is exactly the degree-D Jensen polynomial of "
            "mathcal_F_t^(n), equivalently of the dimension-(2n+2) "
            "radial characteristic profile after normalization."
        ),
        "schoenberg_scope_guard": (
            "The all-dimensions Schoenberg theorem concerns one fixed "
            "radial profile in every dimension. Here both phi_n and its "
            "size-biased radial law change with n. It therefore gives no "
            "complete-monotonicity or Laguerre-Polya promotion; indeed "
            "phi_(0,0) has real Xi zeros and cannot be a nonzero "
            "Gaussian mixture valid in every dimension."
        ),
        "gaussian_newman_flow": (
            "For c>=1, 0<=t<=1/5, let a=c^2-t, b=2*c^2-t and "
            "f^G_(t,c)(u)=exp(-a*u^2)+exp(-b*u^2)="
            "exp(t*u^2)*(exp(-c^2*u^2)+exp(-2*c^2*u^2)). "
            "Then partial_t f^G=u^2*f^G."
        ),
        "gaussian_shape": (
            "Writing q=exp(-c^2*u^2), "
            "(log f^G)''=-2*a-2*c^2*q/(1+q)"
            "+4*c^4*u^2*q/(1+q)^2 <= "
            "-((2-4/e)*c^2-2/5)<0. Hence the family is uniformly "
            "strongly log-concave, and the certified curvature can be "
            "made arbitrarily large by increasing c."
        ),
        "gaussian_abel": (
            "r^G_(t,c)(v)=4*sqrt(pi*a)*exp(-4*a*v)"
            "+4*sqrt(pi*b)*exp(-4*b*v). It is positive on v>0 and "
            "therefore generates the same normalized radial lifts and "
            "dimension walk in every shift."
        ),
        "gaussian_generating": (
            "mathcal_G_(t,c)(z)=sqrt(pi/a)*exp(z/(4*a))"
            "+sqrt(pi/b)*exp(z/(4*b))."
        ),
        "gaussian_zero_lattice": (
            "The complete zero set of mathcal_G_(t,c) is "
            "z_k=2*a*b*log(a/b)/c^2"
            "+4*pi*i*a*b*(2*k+1)/c^2, k in Z. Every zero is nonreal."
        ),
        "gaussian_jensen_failure": (
            "Put p=sqrt(pi/a), q=sqrt(pi/b), alpha=1/(4*a), "
            "beta=1/(4*b), so A_k=p*alpha^k+q*beta^k. Then "
            "disc P_(2,n)=4*(A_(n+1)^2-A_n*A_(n+2))="
            "-4*p*q*(alpha*beta)^n*(alpha-beta)^2<0 for every n>=0."
        ),
        "gaussian_squared_variable_guard": (
            "The Gaussian model does not satisfy the Xi-specific input "
            "used below: y maps to f^G_(t,c)(sqrt(y))="
            "exp(-a*y)+exp(-b*y), whose logarithmic second derivative is "
            "(a-b)^2*exp(-(a+b)*y)/(exp(-a*y)+exp(-b*y))^2>0. "
            "It is strictly log-convex, not log-concave, in y."
        ),
        "xi_quadratic_closure": (
            "The established Xi theorem says y maps to Phi(sqrt(y)) is "
            "strictly log-concave on [0,infinity). Multiplication by "
            "exp(t*y) preserves this for every real t. If "
            "E_k(t)=integral_0^infinity y^(k-1/2)*exp(t*y)"
            "*Phi(sqrt(y))dy, then "
            "A_k(t)=sqrt(pi)*4^(-k)*E_k(t)/Gamma(k+1/2). "
            "Berwald-Borell therefore gives "
            "A_(n+1)^2>=A_n*A_(n+2) for every n>=0 and real t. "
            "Since M_n=n!*A_n, equivalently "
            "M_(n+1)^2>=(n+1)/(n+2)*M_n*M_(n+2), so every shifted "
            "degree-two Jensen polynomial is hyperbolic."
        ),
        "xi_cubic_closure": (
            "The validated reciprocal-defect entry and forward-uniform "
            "heat theorem separately prove that every shifted degree-three "
            "Jensen polynomial is hyperbolic at every finite t>=-100, in "
            "particular throughout 0<=t<=1/5. This conclusion uses the Xi "
            "coefficient heat trajectory and is not a consequence of the "
            "radial ladder or of quadratic log-concavity."
        ),
        "nonpromotion": (
            "Own-dimension radial positive definiteness, the exact "
            "dimension walk, smooth full-support Abel density, Newman "
            "evolution, and even arbitrarily strong log-concavity do not "
            "imply one quadratic Jensen window, let alone PF-infinity."
        ),
        "live_target": (
            "A separate Xi real-zero-band and sector theorem now closes every "
            "shifted Jensen layer through degree 361 uniformly on "
            "0<=t<=1/5. Thus the outer-root quartic threshold u<=U(a,p) "
            "and the adjacent quintic obstruction are no longer open for "
            "the actual Xi heat flow on this interval. The radial ladder "
            "remains an exact coordinate, not the source of that theorem. "
            "The first unproved direct degree is 362, but one more bounded "
            "degree is not the terminal need: the live target is an "
            "unbounded cofinal degree sequence or an all-degree "
            "zero-preserving theorem. A separate "
            "strong-log-concave local Mellin countermodel satisfies all "
            "neighboring ratio and cubic signs while failing its quartic, "
            "and an exact length-ten positive rational segment satisfies the strengthened "
            "scaled-defect and reciprocal-increment corridors together with "
            "all 120, 126, and 56 supported signed-Hankel conditions of "
            "orders two, three, and four across multiple shifts while still "
            "having u>U and a nonhyperbolic adjacent quintic. "
            "A downstream exact Bernstein certificate proves that every "
            "signed order-three/order-four continuation of this fixed prefix "
            "through x_12 fails the necessary increasing-x_13 compatibility "
            "margin, so this witness cannot become an infinite countermodel. "
            "A distinct exact rational tail from the same outer contact "
            "nevertheless clears every stated scalar gate and all 364, 715, "
            "and 792 finite signed minors of orders two through four through "
            "x_13, where its compatibility margin is positive. Its own "
            "x_14 compatibility is negative, but only for that selected tail. "
            "A later length-14 survivor passes 4,043 supported signed minors "
            "through order eight while its adjacent quintic is nonhyperbolic. "
            "These remain valid generic nonpromotion guards, while the Xi "
            "sector theorem excludes them through degree 361 by theta-kernel "
            "geometry. The surviving condition must control unbounded degree."
        ),
        "checks": {
            "derivative_0f1_coefficients": coefficient_checks,
            "radial_normalizations": normalization_checks,
            "radial_characteristic_coefficients": characteristic_checks,
            "gaussian_abel_moments": gaussian_moment_checks,
            "gaussian_quadratic_discriminants": (
                gaussian_discriminant_checks
            ),
            "xi_quadratic_abel_scalings": (
                xi_quadratic_scaling_checks
            ),
        },
    }


def build_payload() -> dict:
    exact = build_exact()
    rows = [
        GateRow(
            id="sardl_01_imported_abel_coordinate",
            role="imported_exact_bridge",
            readiness="ready_to_apply",
            claim="The score/Beta Abel density supplies the moments of the Jensen coefficient sequence.",
            formula=exact["generating_and_moments"],
            proof_boundary="Imports the validated one-sided phase/moment bridge.",
        ),
        GateRow(
            id="sardl_02_radial_probability",
            role="exact_probability_lift",
            readiness="ready_to_apply",
            claim="Every coefficient shift has a normalized radial probability density in an explicit even dimension.",
            formula=(
                exact["radial_density"]
                + " "
                + exact["radial_mixing_law"]
            ),
            proof_boundary="Polar-coordinate identity; the dimension depends on n.",
        ),
        GateRow(
            id="sardl_03_planar_abel_marginal",
            role="exact_radon_identity",
            readiness="ready_to_apply",
            claim="The unshifted radial law is exactly a planar lift of the normalized Newman half-kernel.",
            formula=exact["abel_marginal"],
            proof_boundary="The displayed identity is a one-coordinate marginal, not an all-dimensional extension of one fixed profile.",
        ),
        GateRow(
            id="sardl_04_bessel_derivative",
            role="exact_bessel_identity",
            readiness="ready_to_apply",
            claim="Every derivative of the Jensen generating function is an order-n Bessel mixture of a size-biased Abel law.",
            formula=exact["derivative_bessel"],
            proof_boundary="Exact entire-function identity under the already established moment bounds.",
        ),
        GateRow(
            id="sardl_05_radial_fourier",
            role="exact_positive_definite_theorem",
            readiness="ready_to_apply",
            claim="The nth normalized derivative is a radial characteristic function in dimension 2n+2.",
            formula=exact["radial_characteristic"],
            proof_boundary="Positive definiteness is asserted only in the matching dimension.",
        ),
        GateRow(
            id="sardl_06_dimension_walk",
            role="exact_dimension_walk",
            readiness="ready_to_apply",
            claim="Successive coefficient shifts obey the classical radial derivative dimension walk.",
            formula=exact["dimension_walk"],
            proof_boundary="The normalization and radial law also change with n.",
        ),
        GateRow(
            id="sardl_07_jensen_identification",
            role="exact_window_identification",
            readiness="ready_to_apply",
            claim="Every shifted Jensen polynomial is attached to one member of the radial dimension ladder.",
            formula=exact["jensen_window"],
            proof_boundary="A characteristic function need not have hyperbolic Jensen polynomials.",
        ),
        GateRow(
            id="sardl_08_schoenberg_scope",
            role="exact_scope_guard",
            readiness="guard_validated",
            claim="Schoenberg's all-dimension classification cannot be applied across this moving family.",
            formula=exact["schoenberg_scope_guard"],
            proof_boundary="Blocks a tempting theorem mismatch; it does not weaken the exact fixed-dimension lift.",
        ),
        GateRow(
            id="sardl_09_gaussian_flow_shape",
            role="exact_newman_flow_countermodel",
            readiness="guard_validated",
            claim="A two-Gaussian family shares the Newman deformation and arbitrarily strong log-concavity.",
            formula=(
                exact["gaussian_newman_flow"]
                + " "
                + exact["gaussian_shape"]
            ),
            proof_boundary="The model lacks the Xi theta-series arithmetic and double-exponential tail.",
        ),
        GateRow(
            id="sardl_10_gaussian_radial_lift",
            role="exact_radial_countermodel",
            readiness="guard_validated",
            claim="The countermodel has a positive explicit Abel density and every radial dimension lift.",
            formula=exact["gaussian_abel"],
            proof_boundary="Establishes all generic radial hypotheses used above.",
        ),
        GateRow(
            id="sardl_11_gaussian_jensen_failure",
            role="exact_jensen_countermodel",
            readiness="guard_validated",
            claim="Despite those properties, the countermodel has nonreal entire zeros and fails every shifted quadratic Jensen test.",
            formula=(
                exact["gaussian_generating"]
                + " "
                + exact["gaussian_zero_lattice"]
                + " "
                + exact["gaussian_jensen_failure"]
            ),
            proof_boundary="A generic obstruction, not a counterexample to an Xi-specific arithmetic theorem.",
        ),
        GateRow(
            id="sardl_12_xi_quadratic_closure",
            role="imported_exact_quadratic_closure",
            readiness="ready_to_apply",
            claim="Xi-specific squared-variable log-concavity already closes every shifted quadratic Jensen window.",
            formula=(
                exact["xi_quadratic_closure"]
                + " "
                + exact["gaussian_squared_variable_guard"]
            ),
            proof_boundary="Imports the validated global log-concavity and Berwald-Borell theorem; generic strong log-concavity in u is not the required hypothesis.",
        ),
        GateRow(
            id="sardl_13_xi_cubic_closure",
            role="imported_exact_cubic_closure",
            readiness="ready_to_apply",
            claim="A separate Xi-specific reciprocal-defect flow theorem already closes every shifted cubic Jensen window throughout the target heat interval.",
            formula=exact["xi_cubic_closure"],
            proof_boundary="Imports the validated cubic entry and forward-uniform propagation gates; it is not derived from radial positive definiteness.",
        ),
        GateRow(
            id="sardl_14_imported_degree361_closure",
            role="imported_exact_bounded_degree_closure",
            readiness="ready_to_apply",
            claim="A separate Xi real-zero-band and sector theorem closes every shifted Jensen window through degree 361 on the full target heat interval.",
            formula=exact["live_target"],
            proof_boundary="The cutoff is finite; degree 362, every unbounded cofinal sequence, and every all-degree replacement remain open.",
        ),
        GateRow(
            id="sardl_15_nonpromotion",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="The radial theorem is exact but cannot itself be promoted to PF-infinity or the Newman conclusion.",
            formula=exact["nonpromotion"],
            proof_boundary="This artifact does not prove PF-infinity, Lambda<=0, or RH.",
        ),
    ]
    return {
        "kind": (
            "jensen_window_pf_newman_score_abel_"
            "radial_dimension_lift_gate"
        ),
        "date": "2026-07-25",
        "status": (
            "exact score-Abel radial dimension ladder for every Jensen "
            "shift, reconciled with the completed Xi degree-361 sector "
            "layers and a sharp two-Gaussian Newman-flow nonpromotion gate"
        ),
        "proof_boundary": (
            "This artifact proves the explicit dimension-(2n+2) radial "
            "probability lift, its Fourier/Bessel and dimension-walk "
            "identities, and the exact attachment of every shifted Jensen "
            "window. The two-Gaussian model proves that these facts, even "
            "with arbitrarily strong log-concavity, do not imply a "
            "quadratic Jensen inequality. Separate Xi-specific theorems "
            "and the independent Xi real-zero-band/sector theorem closes every "
            "shifted window through degree 361 on 0<=t<=1/5. This artifact "
            "does not prove degree 362, an unbounded cofinal degree "
            "sequence, an all-degree zero-preserving theorem, PF-infinity, "
            "Lambda<=0, or RH."
        ),
        "sources": [
            "outputs/jensen_window_pf_newman_one_sided_phase_moment_bridge_gate.md",
            "outputs/jensen_window_pf_laguerre_scale_mixture_gate.md",
            "outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md",
            "outputs/jensen_window_pf_cubic_forward_uniform_tail_certificate.md",
            "outputs/jensen_window_pf_strong_logconcave_local_quartic_countermodel.md",
            "outputs/jensen_window_pf_quartic_double_root_threshold_lemma.md",
            "outputs/jensen_window_pf_quartic_signed_hankel_branch_exclusion_lemma.md",
            "outputs/jensen_window_pf_quartic_outer_threshold_order4_nonpromotion_gate.md",
            "outputs/jensen_window_pf_quartic_outer_branch_length13_obstruction.md",
            "outputs/jensen_window_pf_quartic_outer_branch_alternate_length13_survivor_gate.md",
            "outputs/jensen_window_pf_quartic_outer_contact_length14_survivor_gate.md",
            "outputs/jensen_window_pf_newman_zero_slab_degree71_sector_certificate.md",
            "outputs/jensen_window_pf_newman_real_zero_band_degree361_sector_certificate.md",
            "outputs/jensen_window_pf_newman_gasper_fake_xi_remainder_gate.md",
            "https://dlmf.nist.gov/10.39",
            "https://doi.org/10.2307/1968466",
            "https://arxiv.org/abs/2408.11612",
        ],
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    return "\n".join(
        [
            "# Score-Abel Radial Dimension-Lift Gate",
            "",
            "Date: 2026-07-25",
            "",
            "Status: exact radial dimension ladder for every shifted Jensen",
            "window, with a sharp two-Gaussian nonpromotion gate. This is",
            "not a proof of PF-infinity, RH, or `Lambda <= 0`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_score_abel_radial_dimension_lift_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_score_abel_radial_dimension_lift_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_score_abel_radial_dimension_lift_gate.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated score-Abel radial dimension-lift gate: 15 rows, 0 issues, 1 imported Abel coordinate, 1 radial probability ladder, 1 planar marginal, 1 Bessel derivative identity, 1 own-dimension positive-definiteness theorem, 1 dimension walk, 1 Jensen identification, 1 Schoenberg scope guard, 3 exact Gaussian countermodel rows, 2 imported Xi low-degree closures, 1 imported degree-361 closure, 1 nonpromotion gate",
            "```",
            "",
            "## Radial Probability Ladder",
            "",
            "Start from the exact Abel representation already obtained:",
            "",
            "```text",
            exact["generating_and_moments"],
            "```",
            "",
            "For every shift `n>=0`, define",
            "",
            "```text",
            exact["radial_density"],
            exact["radial_mixing_law"],
            "```",
            "",
            "Indeed, the sphere area in dimension `2n+2` and the change",
            "`v=|y|^2/4` leave the factor `1/n!`; the remaining radial",
            "moment is `M_n=n!*A_n`. Thus the density has mass one.",
            "",
            "At the bottom of the ladder there is a direct geometric",
            "interpretation:",
            "",
            "```text",
            exact["abel_marginal"],
            "```",
            "",
            "So the Abel transform was not merely a moment device: it is the",
            "isotropic planar density whose one-coordinate marginal is the",
            "normalized Newman kernel.",
            "",
            "## Bessel-Fourier Identity",
            "",
            "Termwise differentiation gives",
            "",
            "```text",
            exact["derivative_bessel"],
            "```",
            "",
            "The right side is precisely the spherical characteristic",
            "kernel in dimension `2n+2`. Therefore",
            "",
            "```text",
            exact["radial_characteristic"],
            exact["dimension_walk"],
            "```",
            "",
            "This makes every derivative a characteristic function, but in",
            "its own increasing dimension and under its own size-biased law.",
            "",
            "## Jensen Windows",
            "",
            "```text",
            exact["jensen_window"],
            "```",
            "",
            "The radial ladder therefore captures every shift, not only the",
            "unshifted window. It also tells us exactly why positive",
            "definiteness alone is not the missing theorem.",
            "",
            "## Schoenberg Scope Guard",
            "",
            "```text",
            exact["schoenberg_scope_guard"],
            "```",
            "",
            "The fixed-profile all-dimension hypothesis is absent. Treating",
            "the sequence `phi_n` as one Schoenberg profile would silently",
            "replace the actual size-biased family by a stronger false",
            "premise.",
            "",
            "## Two-Gaussian Newman Countermodel",
            "",
            "The failure can be made exact while retaining the same heat-flow",
            "deformation:",
            "",
            "```text",
            exact["gaussian_newman_flow"],
            exact["gaussian_shape"],
            exact["gaussian_abel"],
            "```",
            "",
            "Thus the countermodel has every radial lift above, full Abel",
            "support, and an arbitrarily large certified strong-concavity",
            "constant. Nevertheless,",
            "",
            "```text",
            exact["gaussian_generating"],
            exact["gaussian_zero_lattice"],
            exact["gaussian_jensen_failure"],
            "```",
            "",
            "So it fails the degree-two Jensen test for every shift. The",
            "obstruction is stronger than a single accidental polynomial:",
            "the entire generating function has a complete nonreal zero",
            "lattice.",
            "",
            "## Quadratic And Cubic Reconciliation",
            "",
            "The Gaussian countermodel does not conflict with the existing",
            "Xi low-degree theorems:",
            "",
            "```text",
            exact["gaussian_squared_variable_guard"],
            exact["xi_quadratic_closure"],
            exact["xi_cubic_closure"],
            "```",
            "",
            "Thus the Abel concentration inequality in degree two is not",
            "open, and the separate reciprocal-defect heat theorem also",
            "closes degree three throughout the target heat interval.",
            "",
            "## Degree-361 Reconciliation",
            "",
            "```text",
            exact["live_target"],
            "```",
            "",
            "The radial theorem is therefore useful as a coordinate and a",
            "theorem-search filter. The Xi sector theorem now closes the",
            "earlier quartic and quintic frontier, and indeed every shifted",
            "degree through 361, without promoting the generic radial facts.",
            "The direct finite-degree frontier begins at degree 362. The",
            "mathematically relevant next step is an unbounded cofinal",
            "sequence or a genuine all-degree theorem.",
            "",
            "References: https://dlmf.nist.gov/10.39 and",
            "https://doi.org/10.2307/1968466; for modern dimension-walk",
            "context see https://arxiv.org/abs/2408.11612.",
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
    args.out.write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote score-Abel radial dimension-lift gate: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
