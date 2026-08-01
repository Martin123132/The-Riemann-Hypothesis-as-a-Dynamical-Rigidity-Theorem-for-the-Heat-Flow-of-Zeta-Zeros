#!/usr/bin/env python3
"""Derive the parabolic-frequency contact normal form and hierarchy guards."""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_parabolic_frequency_contact_normal_hierarchy_gate"
DATE = "2026-07-26"
DEFAULT_OUT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / f"outputs/{STEM}.md"
SOURCE_FILES = {
    "ray_aligned_reduction": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_ray_aligned_parabolic_frequency_reduction.json"
    ),
    "critical_transversality": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_transversality_target.json"
    ),
    "critical_c1_global_remainder": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate.json"
    ),
    "correlation_hierarchy": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate.json"
    ),
    "crossing_slope_gap": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_crossing_slope_gap_reduction.json"
    ),
    "scaled_successor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_time_dependent_scaled_successor_lemma.json"
    ),
    "lehmer_margin": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_lehmer_margin_gate.json"
    ),
}


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_hashes() -> dict[str, str]:
    missing = [str(path) for path in SOURCE_FILES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing source artifacts: {missing}")
    return {name: file_hash(path) for name, path in SOURCE_FILES.items()}


def laguerre_coefficients() -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    y = sp.symbols("y", real=True)
    h = sp.symbols("h0:7", real=True)
    plus = sum(h[k] * (sp.I * y) ** k / sp.factorial(k) for k in range(7))
    minus = sum(h[k] * (-sp.I * y) ** k / sp.factorial(k) for k in range(7))
    product = sp.expand(plus * minus)
    return tuple(sp.simplify(product.coeff(y, 2 * n)) for n in range(1, 4))


def build_exact() -> dict:
    t, ell = sp.symbols("t ell", positive=True)
    q = 2 * t * ell**2
    scale = sp.sqrt(2 * t / (1 + q))
    beta = sp.simplify(ell * scale)
    if sp.simplify(beta**2 - q / (1 + q)) != 0:
        raise RuntimeError("frequency-chart weight identity failed")
    if sp.simplify(scale**2 * ell**2 - beta**2) != 0:
        raise RuntimeError("scaled derivative remainder identity failed")
    slope_sum = sp.symbols("S", nonnegative=True)
    decay = sp.exp(-3 * ell / 4)
    epsilon_0 = 2500 * decay
    epsilon_1 = 5000 * decay
    explicit_slope_budget = sp.simplify(
        ell * epsilon_1 / 2 + epsilon_0 * slope_sum / 4
    )
    expected_slope_budget = 625 * decay * (4 * ell + slope_sum)
    if sp.simplify(explicit_slope_budget - expected_slope_budget) != 0:
        raise RuntimeError("explicit robust slope budget identity failed")
    mass_plus, mass_minus, h_plus, h_minus = sp.symbols(
        "M_plus M_minus h_plus h_minus", real=True
    )
    mean_mass = (mass_plus + mass_minus) / 2
    real_value = mass_plus - mass_minus
    slope_decomposition = sp.expand(
        mean_mass * (h_plus - h_minus)
        + real_value * (h_plus + h_minus) / 2
    )
    if sp.simplify(
        slope_decomposition - (mass_plus * h_plus - mass_minus * h_minus)
    ) != 0:
        raise RuntimeError("zero-mass slope convention identity failed")

    a0, a1, a2, a3, z0, z1, z2, z3 = sp.symbols(
        "a0 a1 a2 a3 z0 z1 z2 z3", real=True
    )
    heat_xx = (a2 * z0 + 2 * a1 * z1 + a0 * z2).subs({z0: 0, z1: 0})
    heat_xxx = (
        a3 * z0 + 3 * a2 * z1 + 3 * a1 * z2 + a0 * z3
    ).subs({z0: 0, z1: 0})
    expected_xx = a0 * z2
    expected_xxx = a0 * (z3 + 3 * (a1 / a0) * z2)
    if sp.simplify(heat_xx - expected_xx) != 0:
        raise RuntimeError("normalized contact second-jet identity failed")
    if sp.simplify(heat_xxx - expected_xxx) != 0:
        raise RuntimeError("normalized contact third-jet identity failed")

    h1, h2, h3, s, s_t = sp.symbols(
        "h1 h2 h3 s s_t", real=True, finite=True
    )
    scaled_heat_norm = h2**2 + (s_t * h1 - s * h3) ** 2
    contact_heat_norm = sp.simplify(scaled_heat_norm.subs(h1, 0))
    if contact_heat_norm != h2**2 + s**2 * h3**2:
        raise RuntimeError("contact scaled-heat norm identity failed")

    laguerre_1, laguerre_2, laguerre_3 = laguerre_coefficients()
    h0, h1c, h2c, h3c, h4c, h5c, h6c = sp.symbols(
        "h0 h1 h2 h3 h4 h5 h6", real=True
    )
    expected_l1 = h1c**2 - h0 * h2c
    expected_l2 = h2c**2 / 4 - h1c * h3c / 3 + h0 * h4c / 12
    expected_l3 = (
        h3c**2 / 36
        - h2c * h4c / 24
        + h1c * h5c / 60
        - h0 * h6c / 360
    )
    for actual, expected, label in (
        (laguerre_1, expected_l1, "L1"),
        (laguerre_2, expected_l2, "L2"),
        (laguerre_3, expected_l3, "L3"),
    ):
        if sp.simplify(actual - expected) != 0:
            raise RuntimeError(f"{label} coefficient identity failed")

    f1 = laguerre_1
    f2 = 3 * laguerre_2
    f3 = sp.Rational(45, 2) * laguerre_3
    f1_xixi = (h2c**2 - h0 * h4c) / 4
    hierarchy_defect = sp.simplify(f2 - 3 * f1_xixi)
    if sp.simplify(hierarchy_defect - (h0 * h4c - h1c * h3c)) != 0:
        raise RuntimeError("correlation hierarchy defect identity failed")

    double_contact = {h0: 0, h1c: 0}
    if sp.simplify(f1.subs(double_contact)) != 0:
        raise RuntimeError("F1 double-contact value failed")
    if sp.simplify(f1_xixi.subs(double_contact) - h2c**2 / 4) != 0:
        raise RuntimeError("F1 double-contact curvature failed")
    if sp.simplify(f2.subs(double_contact) - 3 * h2c**2 / 4) != 0:
        raise RuntimeError("F2 double-contact value failed")
    order_hankel = sp.expand(f1 * f3 - f2**2)
    if sp.simplify(
        order_hankel.subs(double_contact) + sp.Rational(9, 16) * h2c**4
    ) != 0:
        raise RuntimeError("correlation-order Hankel contact signature failed")

    y, tau = sp.symbols("y tau", real=True)
    quadratic_jet = {
        h0: y**2 - 2 * tau,
        h1c: 2 * y,
        h2c: 2,
        h3c: 0,
        h4c: 0,
        h5c: 0,
        h6c: 0,
    }
    quadratic_f1 = sp.simplify(f1.subs(quadratic_jet))
    quadratic_f2 = sp.simplify(f2.subs(quadratic_jet))
    quadratic_f3 = sp.simplify(f3.subs(quadratic_jet))
    quadratic_hankel = sp.simplify(order_hankel.subs(quadratic_jet))
    if quadratic_f1 != 2 * y**2 + 4 * tau:
        raise RuntimeError("quadratic F1 identity failed")
    if quadratic_f2 != 3 or quadratic_f3 != 0 or quadratic_hankel != -9:
        raise RuntimeError("quadratic hierarchy countermodel failed")
    if sp.simplify(sp.diff(y**2 - 2 * tau, tau) + sp.diff(y**2 - 2 * tau, y, 2)) != 0:
        raise RuntimeError("quadratic Newman heat equation failed")

    return {
        "scale": "s_pf=sqrt(2t)/sqrt(1+q), q=2tL^2",
        "chart_weight": "beta=L*s_pf=sqrt(q/(1+q))",
        "corrected_main_jet": (
            "J=2X-Q, J_x=2(U-BY)-Q_x, "
            "T_pf[J]=J^2+s_pf^2*J_x^2"
        ),
        "contact_remainder_transfer": (
            "if Z=J+r=Z_x=J_x+r_x=0, |r|<=epsilon_0, "
            "|r_x|<=L*epsilon_1, then "
            "T_pf[J]<=epsilon_0^2+beta^2*epsilon_1^2"
        ),
        "critical_c1_budget": (
            "for L>=50 and 0<tL<=25: "
            "epsilon_0=2500*exp(-3L/4), "
            "epsilon_1=5000*exp(-3L/4)"
        ),
        "sharp_contact_box": (
            "a contact requires |J|<=epsilon_0 and "
            "|J_x|<=L*epsilon_1"
        ),
        "explicit_sharp_contact_box": (
            "a contact requires |J|<=2500*exp(-3L/4) and "
            "|J_x|<=5000*L*exp(-3L/4)"
        ),
        "sharp_c1_target": (
            "|2X-Q|>epsilon_0 or "
            "|2(U-BY)-Q_x|>L*epsilon_1"
        ),
        "crossing_band_target": (
            "on |2X-Q|<=epsilon_0 prove "
            "|2(U-BY)-Q_x|>L*epsilon_1"
        ),
        "corrected_component_box": (
            "for E_c with J=2Re(E_c), put X_c=Re(E_c), "
            "U_c=Re(E_c,x); a contact requires "
            "|X_c|<=epsilon_0/2 and |U_c|<=L*epsilon_1/2"
        ),
        "near_crossing_slope_identity": (
            "U_c=M(h_+-h_-)+(X_c/2)(h_++h_-)+D_0"
        ),
        "zero_mass_slope_convention": (
            "set h_+=0 when M_+=0 and h_-=0 when M_-=0; "
            "the near-crossing identity then remains valid, including "
            "the all-zero-real-part case U_c=D_0"
        ),
        "robust_band_slope_target": (
            "|M(h_+-h_-)+D_0|>L*epsilon_1/2"
            "+epsilon_0*|h_++h_-|/4"
        ),
        "explicit_robust_band_slope_target": (
            "|M(h_+-h_-)+D_0|>"
            "625*exp(-3L/4)*(4L+|h_++h_-|)"
        ),
        "phase_matched_contact_equation": (
            "at a contact, M(h_+-h_-)+D_0="
            "-r_x/2+(r/4)(h_++h_-)"
        ),
        "signed_remainder_peeling": (
            "if r=r_0+delta_0, r_x=r_(0,x)+delta_1, "
            "|delta_0|<=eta_0, and |delta_1|<=L*eta_1, then contact "
            "requires |M(h_+-h_-)+D_0+r_(0,x)/2"
            "-(r_0/4)(h_++h_-)|<="
            "L*eta_1/2+(eta_0/4)|h_++h_-|"
        ),
        "signed_remainder_target": (
            "|M(h_+-h_-)+D_0+r_(0,x)/2"
            "-(r_0/4)(h_++h_-)|>"
            "L*eta_1/2+(eta_0/4)|h_++h_-|"
        ),
        "uniform_wedge_multiplicity_guard": (
            "if at t=0 a contact has |J(0)|<epsilon_0 and "
            "|J_x(0)|<L*epsilon_1, continuity keeps J,J_x inside the "
            "same open box for all sufficiently small t>0; therefore a "
            "full-wedge pointwise C1-box exclusion rules out endpoint "
            "multiplicity and is stronger than Lambda<=0 or RH"
        ),
        "strict_c1_target": (
            "(2X-Q)^2+s_pf^2*(2(U-BY)-Q_x)^2 "
            ">epsilon_0^2+beta^2*epsilon_1^2"
        ),
        "domain_split": (
            "q<=1 iff t<=1/(2L^2) iff c=tL<=1/(2L); "
            "q>=1 is the frequency chart and contains every fixed "
            "0<c<=c_*+epsilon as L->infinity"
        ),
        "frequency_affine_reduction": (
            "put A0=J^2-epsilon_0^2 and "
            "A1=(J_x/L)^2-epsilon_1^2; for q>=1, "
            "D_pf=A0+beta^2*A1>=min(A0+A1/2,A0+A1), "
            "so positivity of both endpoint deficits suffices"
        ),
        "normalizer_contact_jet": (
            "for H=A Z and a=(log A)_x, at Z=Z_x=0: "
            "H_xx=A Z_xx, H_xxx=A(Z_xxx+3a Z_xx)"
        ),
        "scaled_heat_contact": (
            "partial_t(H,s_pf H_x)=-A("
            "Z_xx,s_pf*(Z_xxx+3a Z_xx)) at a contact"
        ),
        "scale_derivative_guard": (
            "the s_pf,t H_x term vanishes at H_x=0; the scale cannot "
            "regularize a genuine collision"
        ),
        "route_order_guard": (
            "direct corrected-main exclusion uses a C1 remainder; "
            "relative heat-jet closure exposes Z_xx and Z_xxx and therefore "
            "needs C3 control or an independent closed factorization"
        ),
        "correlation_normalization": (
            "F_1(2x)=H_x^2-H H_xx; "
            "F_2(2x)=3H_xx^2/4-H_x H_xxx+H H_xxxx/4"
        ),
        "double_contact_signature": (
            "H=H_x=0, H_xx!=0 implies F_1=F_1'=0, "
            "F_1''=H_xx^2/4, F_2=3H_xx^2/4, "
            "partial_t F_1=H_xx^2"
        ),
        "hierarchy_defect": (
            "F_2-3*partial_xi^2 F_1=H H_xxxx-H_x H_xxx"
        ),
        "order_hankel_contact": (
            "at a nondegenerate double contact, "
            "F_1*F_3-F_2^2=-9*H_xx^4/16<0"
        ),
        "quadratic_countermodel": (
            "H_tau(y)=y^2-2tau solves H_tau=-H_yy; "
            "F_1=2y^2+4tau, F_2=3, F_3=0, "
            "F_1F_3-F_2^2=-9, and at y=0 "
            "||partial_tau(H,sH_y)||/||(H,sH_y)||=1/|tau|"
        ),
    }


def gate_rows(exact: dict) -> list[GateRow]:
    return [
        GateRow(
            id="pfcn_01_positive_normalizer_contact_transfer",
            role="exact_identity",
            readiness="available_exact",
            claim="A positive normalizer preserves first-jet contacts exactly.",
            formula="H=A Z, A>0: H=H_x=0 iff Z=Z_x=0",
            proof_boundary="Algebraic equivalence only.",
        ),
        GateRow(
            id="pfcn_02_corrected_main_scaled_jet",
            role="exact_identity",
            readiness="available_exact",
            claim="The corrected finite main has an exact parabolic-frequency first-jet norm.",
            formula=exact["corrected_main_jet"],
            proof_boundary="No lower bound for this norm is asserted.",
        ),
        GateRow(
            id="pfcn_03_c1_remainder_transfer",
            role="exact_reduction",
            readiness="available_exact",
            claim="A refined C1 remainder puts every possible contact in one sharp coordinate box.",
            formula=(
                exact["critical_c1_budget"]
                + "; "
                + exact["explicit_sharp_contact_box"]
                + "; sufficient exclusion: "
                + exact["sharp_c1_target"]
            ),
            proof_boundary="The joint Xi small-box exclusion remains open.",
        ),
        GateRow(
            id="pfcn_04_two_chart_domain_split",
            role="exact_reduction",
            readiness="available_exact",
            claim="The parabolic chart is ultra-small, while two endpoint inequalities control the fixed-c frequency layer.",
            formula=(
                exact["domain_split"]
                + "; "
                + exact["frequency_affine_reduction"]
            ),
            proof_boundary="A domain decomposition, not contact exclusion.",
        ),
        GateRow(
            id="pfcn_05_normalized_contact_higher_jet",
            role="exact_identity",
            readiness="available_exact",
            claim="At contact, the physical second and third jets have a simple normalized form.",
            formula=exact["normalizer_contact_jet"],
            proof_boundary="Does not bound either higher jet.",
        ),
        GateRow(
            id="pfcn_06_scaled_heat_contact_vector",
            role="exact_identity",
            readiness="available_exact",
            claim="The scaled heat derivative at contact exposes the normalized second and third jets.",
            formula=exact["scaled_heat_contact"],
            proof_boundary="This is the numerator of a relative estimate, not a bound.",
        ),
        GateRow(
            id="pfcn_07_scale_derivative_nonregularization",
            role="exact_route_guard",
            readiness="closed_route",
            claim="Time dependence of the positive scale cannot regularize a first-jet collision.",
            formula=exact["scale_derivative_guard"],
            proof_boundary="Rejects a scale-only closure; Xi structure may still bound the higher jets.",
        ),
        GateRow(
            id="pfcn_08_derivative_order_route_guard",
            role="exact_route_guard",
            readiness="closed_route",
            claim="Direct corrected-main exclusion has a strictly lower remainder-derivative burden than bare relative transport.",
            formula=exact["route_order_guard"],
            proof_boundary="A closed Xi factorization could still make the relative route competitive.",
        ),
        GateRow(
            id="pfcn_09_correlation_first_two_layers",
            role="exact_identity",
            readiness="available_exact",
            claim="The first two correlation transforms are explicit differential polynomials in H.",
            formula=exact["correlation_normalization"],
            proof_boundary="Correlation identities do not provide a favorable pointwise sign.",
        ),
        GateRow(
            id="pfcn_10_double_contact_hierarchy_signature",
            role="exact_identity",
            readiness="available_exact",
            claim="A nondegenerate contact has a rigid quadratic-touch signature in the correlation hierarchy.",
            formula=exact["double_contact_signature"],
            proof_boundary="The hierarchy is compatible with this signature.",
        ),
        GateRow(
            id="pfcn_11_hierarchy_defect",
            role="exact_identity",
            readiness="available_exact",
            claim="The proposed F2 versus F1 curvature comparison reduces to a fourth-order Wronskian defect.",
            formula=exact["hierarchy_defect"],
            proof_boundary="No Xi sign for the defect is asserted.",
        ),
        GateRow(
            id="pfcn_12_order_hankel_countermodel",
            role="exact_countermodel",
            readiness="closed_counterexample",
            claim="Correlation-order log-convexity is not a generic or Laguerre-Polya closure.",
            formula=exact["quadratic_countermodel"],
            proof_boundary="The quadratic flow is not Xi; it rejects only a generic promotion.",
        ),
        GateRow(
            id="pfcn_13_positive_definite_nonpromotion",
            role="literature_guard",
            readiness="literature_nonpromotion",
            claim="Positive definiteness of the first Wronskian is not pointwise strict positivity.",
            formula=(
                "Dimitrov-Xu make F_1 positive definite because its inverse "
                "Fourier correlation is nonnegative; positive-definite "
                "functions may have real zeros, so this does not exclude contact"
            ),
            proof_boundary="The cited density criterion for Laguerre-Polya membership is itself RH-equivalent here.",
        ),
        GateRow(
            id="pfcn_14_phase_matched_remainder_equation",
            role="exact_reduction",
            readiness="available_exact",
            claim="At a hypothetical contact the signed remainder obeys an exact phase-matched slope equation.",
            formula=(
                exact["zero_mass_slope_convention"]
                + "; "
                + exact["phase_matched_contact_equation"]
            ),
            proof_boundary="The equation retains remainder correlation but does not contradict it.",
        ),
        GateRow(
            id="pfcn_15_signed_remainder_peeling",
            role="exact_reduction",
            readiness="available_exact",
            claim="A signed leading remainder can be moved to the signal side so only its smaller residual is paid as an absolute error.",
            formula=exact["signed_remainder_peeling"],
            proof_boundary="No Xi signed leading remainder or residual improvement is supplied here.",
        ),
        GateRow(
            id="pfcn_16_endpoint_multiplicity_quantifier_guard",
            role="exact_route_guard",
            readiness="closed_route",
            claim="A pointwise sharp-box theorem uniform down to t=0 is stronger than the multiplicity-compatible Newman conclusion.",
            formula=exact["uniform_wedge_multiplicity_guard"],
            proof_boundary="The stronger theorem remains sufficient; it must not be presented as logically necessary.",
        ),
        GateRow(
            id="pfcn_17_live_xi_c1_target",
            role="theorem_target",
            readiness="open",
            claim="The direct frequency-layer route is a signed Riemann-Siegel remainder peeling theorem; the unpeeled magnitude gap on the corrected-component value band is a stronger fallback.",
            formula=(
                exact["phase_matched_contact_equation"]
                + "; preferred signed refinement: "
                + exact["signed_remainder_target"]
                + "; unpeeled sufficient magnitude relaxation: "
                + exact["robust_band_slope_target"]
                + "; equivalently under the certified critical C1 budget: "
                + exact["explicit_robust_band_slope_target"]
            ),
            proof_boundary="The Xi phase-matched contradiction, stronger arithmetic gap, ultra-small theorem, and bounded-L shoulder are unproved.",
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    rows = gate_rows(exact)
    return {
        "kind": STEM,
        "date": DATE,
        "status": (
            "exact parabolic-frequency contact normal form, hierarchy route "
            "guards, and one open Xi C1 target"
        ),
        "proof_boundary": (
            "The artifact proves contact-normal, corrected-main, normalizer, "
            "heat-jet, and correlation identities. It rejects scale-only "
            "regularization and generic correlation-order log-convexity. It "
            "does not prove the Xi phase-matched contact contradiction, "
            "strict C1 small-ball theorem, contact "
            "exclusion in the wedge, Lambda<=0, RH, or the Clay prize."
        ),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": 17,
            "exact_identities": 12,
            "exact_route_guards": 3,
            "literature_guards": 1,
            "open_xi_targets": 1,
        },
        "sources": {
            "Dimitrov_Xu": "https://arxiv.org/abs/1606.05011",
            "Rodgers_Tao": "https://arxiv.org/abs/1801.05914",
            "Polymath15": "https://arxiv.org/abs/1904.12438",
        },
        "source_sha256": source_hashes(),
        "builder_sha256": file_hash(Path(__file__)),
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return f"""# Newman Parabolic-Frequency Contact-Normal Hierarchy Gate

Date: {DATE}

Status: exact contact-normal reduction, exact route guards, and one open
Xi `C1` theorem target. This is not a proof of `Lambda<=0` or RH.

## Corrected Contact Normal Form

Write `H=A Z`, where the Polymath-15 normalizer `A` is positive, and use

```text
Z=J+r,
{exact["corrected_main_jet"]},
{exact["scale"]},
{exact["chart_weight"]}.
```

At a first-jet contact, `Z=Z_x=0`, so `J=-r` and `J_x=-r_x`. Therefore

```text
{exact["critical_c1_budget"]},
{exact["contact_remainder_transfer"]}.             (1)
{exact["sharp_contact_box"]}.                      (1a)
{exact["explicit_sharp_contact_box"]}.             (1b)
```

The sharp noncircular sufficient theorem is

```text
{exact["sharp_c1_target"]},                         (2)
```

or, equivalently,

```text
{exact["crossing_band_target"]}.                    (2a)
```

This requires derivative separation only on the thin value band. The
smooth scaled inequality

```text
{exact["strict_c1_target"]}                         (2b)
```

is sufficient but stronger than the sharp box exclusion. All three use only
the refined value and first-derivative remainder.

For the endpoint-corrected complex component sum `E_c`, write

```text
{exact["corrected_component_box"]},
{exact["near_crossing_slope_identity"]},
{exact["zero_mass_slope_convention"]}.
```

The triangle inequality reduces the frequency-layer band theorem to

```text
{exact["robust_band_slope_target"]},                (2c)
{exact["explicit_robust_band_slope_target"]}.       (2d)
```

This is a weighted arithmetic slope gap only where the corrected value is
inside its remainder band. It is weaker than the previous strip-wide cone
and retains the finite remainder terms omitted by an exact-crossing floor.

At an actual contact the signed remainder gives the sharper phase-matched equation

```text
{exact["phase_matched_contact_equation"]}.          (2e)
```

The magnitude gap (2c)-(2d) is obtained by applying the triangle inequality
to (2e). It is a sufficient fallback, but it discards the phase and derivative
correlation of the actual Riemann-Siegel remainder. The primary
contact-conditioned search should retain (2e) and insert a signed asymptotic
model for `r,r_x`.

More precisely, peeling a signed leading remainder gives the exact reduction

```text
{exact["signed_remainder_peeling"]}.                (2f)
```

Thus the next analytic deliverable is an explicit Riemann-Siegel
`r_0,r_(0,x)` and smaller certified residual `eta_0,eta_1`, followed by

```text
{exact["signed_remainder_target"]}.                 (2g)
```

The unpeeled magnitude target is the special case `r_0=r_(0,x)=0`.

There is also an essential endpoint quantifier guard:

```text
{exact["uniform_wedge_multiplicity_guard"]}.        (2h)
```

Thus a pointwise box exclusion uniform down to `t=0` may still prove the
desired result, but it is stronger than the multiplicity-compatible theorem.
The cofinal boundary-degree and delta-localized transport routes remain live.

## Exact Chart Split

```text
{exact["domain_split"]}.                            (3)
```

Thus `q<=1` is the ultra-small layer `t<=1/(2L^2)`. Every asymptotic ray
with fixed positive `c=tL` eventually lies in the frequency chart. The
proof search should not treat the whole critical layer as parabolic.

For the frequency layer the scale parameter can be removed exactly:

```text
{exact["frequency_affine_reduction"]}.              (3a)
```

It is enough to prove the two endpoint deficits at derivative weights
`1/2` and `1`; no continuum of beta values needs a separate theorem.

## Why Scale Alone Cannot Close Theorem

At `Z=Z_x=0`,

```text
{exact["normalizer_contact_jet"]},
{exact["scaled_heat_contact"]}.                     (4)
```

The term `s_(pf,t) H_x` vanishes exactly at a contact. Consequently a
time-dependent scale conditions the approach to a collision but does not
regularize the collision itself. A direct proof of (2) uses a `C1`
remainder, whereas a bare relative heat-jet proof exposes `Z_xx,Z_xxx`
and needs `C3` control or an independently closed factorization.

The exact Newman heat flow

```text
H_tau(y)=y^2-2tau
```

has, at `y=0`,

```text
||partial_tau(H_tau,sH_(tau,y))||/
||(H_tau,sH_(tau,y))||=1/|tau|.
```

The pole is nonintegrable through the double contact for every smooth
positive scale. This is an exact generic guard, not an Xi counterexample.

## Correlation-Hierarchy Translation

The first layers satisfy

```text
{exact["correlation_normalization"]},
{exact["hierarchy_defect"]}.                        (5)
```

At a nondegenerate double contact,

```text
{exact["double_contact_signature"]},
{exact["order_hankel_contact"]}.                    (6)
```

A tempting order-Hankel condition `F_1 F_3-F_2^2>=0` would therefore
exclude contact. It is not a viable generic closure: for the real-rooted
quadratic heat flow,

```text
{exact["quadratic_countermodel"]}.                  (7)
```

Hence correlation-order log-convexity is not even necessary for a
Laguerre-Polya polynomial. Any surviving hierarchy inequality must use
additional Xi/theta structure.

Dimitrov and Xu prove that the first Wronskian is positive definite as a
function because its inverse Fourier correlation is nonnegative. Positive
definiteness is not pointwise strict positivity and permits real zeros; their
density characterization becomes an RH-equivalent condition for the Xi
kernel rather than the missing implication.

## Route Decision

The direct contact route is the signed peeling target (2g), divided into:

```text
q<=1: an ultra-small-time Hermite/parabolic contact theorem;
q>=1 and 0<tL<=c_*+epsilon: prove (2g), with the robust band slope gap
as a stronger sufficient fallback;
bounded L: a terminating finite shoulder theorem.
```

Weighted crossing slopes remain an admissible way to prove the second line
at the contact equations themselves. A strip-wide cone is no longer needed.
An endpoint-uniform floor is not logically required, and no proposed
majorant may contain the inverse unknown first-jet norm.

Primary sources:

```text
https://arxiv.org/abs/1606.05011
https://arxiv.org/abs/1801.05914
https://arxiv.org/abs/1904.12438
```

## Proof Boundary

The exact identities and countermodels above do not prove (2a), the
ultra-small-time Xi theorem, the critical asymptotic Xi theorem, the
bounded-L shoulder, `Lambda<=0`, RH, PF-infinity, or the Clay prize.

Machine-audited files:

```text
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    summary = artifact["summary"]
    print(
        "built Newman parabolic-frequency contact-normal hierarchy gate: "
        f"{summary['rows']} rows, "
        f"{summary['exact_identities']} exact identities, "
        f"{summary['exact_route_guards']} exact route guards, "
        f"{summary['literature_guards']} literature guard, "
        f"{summary['open_xi_targets']} open Xi target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
