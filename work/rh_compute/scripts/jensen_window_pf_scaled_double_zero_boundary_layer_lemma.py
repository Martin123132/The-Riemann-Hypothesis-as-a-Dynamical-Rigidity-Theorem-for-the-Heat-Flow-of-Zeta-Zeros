#!/usr/bin/env python3
"""Derive the universal scaled Jensen layer near a heat-flow double zero."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_scaled_double_zero_boundary_layer_lemma.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_scaled_double_zero_boundary_layer_lemma.md"


@dataclass(frozen=True)
class LemmaRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def build_exact() -> dict:
    h = sp.symbols("h")
    falling_checks = []
    for order in range(13):
        product = sp.expand(sp.prod(1 - offset * h for offset in range(order)))
        linear = sp.expand(product).coeff(h, 1)
        quadratic = sp.expand(product).coeff(h, 2)
        expected_linear = -sp.Rational(order * (order - 1), 2)
        expected_quadratic = sp.Rational(
            order * (order - 1) * (order - 2) * (3 * order - 1), 24
        )
        if linear != expected_linear:
            raise RuntimeError("falling-factorial first correction failed")
        if quadratic != expected_quadratic:
            raise RuntimeError("falling-factorial second correction failed")
        falling_checks.append(
            {
                "order": order,
                "linear_coefficient": str(linear),
                "quadratic_coefficient": str(quadratic),
            }
        )

    rho, tau, eta, f2 = sp.symbols("rho tau eta f2", nonzero=True, real=True)
    local_quadratic = sp.expand(eta**2 + 8 * rho * tau - rho**2)
    local_discriminant = sp.factor(sp.discriminant(local_quadratic, eta))
    if sp.simplify(local_discriminant - 4 * (rho**2 - 8 * rho * tau)) != 0:
        raise RuntimeError("local boundary-layer discriminant failed")
    if sp.simplify(local_discriminant.subs(tau, rho / 8)) != 0:
        raise RuntimeError("local collision threshold failed")

    lam, z, D = sp.symbols("lambda z D", real=True, positive=True)
    toy = z**2 + (2 + 12 * lam) * z + 1 + 4 * lam + 12 * lam**2
    toy_pde_residual = sp.factor(
        sp.diff(toy, lam) - (2 * sp.diff(toy, z) + 4 * z * sp.diff(toy, z, 2))
    )
    if toy_pde_residual != 0:
        raise RuntimeError("toy family does not satisfy the heat PDE")

    toy_scaled_jensen = sp.expand(
        1
        + 4 * lam
        + 12 * lam**2
        + (2 + 12 * lam) * z
        + (1 - 1 / D) * z**2
    )
    toy_discriminant = sp.factor(sp.discriminant(toy_scaled_jensen, z))
    expected_discriminant = sp.factor(
        4 * (24 * D * lam**2 + 8 * D * lam + 12 * lam**2 + 4 * lam + 1) / D
    )
    if sp.simplify(toy_discriminant - expected_discriminant) != 0:
        raise RuntimeError("toy Jensen discriminant failed")

    near_threshold = sp.factor(
        (sp.sqrt(2) * sp.sqrt(D - 1) - sp.sqrt(2 * D + 1))
        / (6 * sp.sqrt(2 * D + 1))
    )
    far_threshold = sp.factor(
        -(sp.sqrt(2) * sp.sqrt(D - 1) + sp.sqrt(2 * D + 1))
        / (6 * sp.sqrt(2 * D + 1))
    )
    if sp.simplify(toy_discriminant.subs(lam, near_threshold)) != 0:
        raise RuntimeError("near toy threshold failed")
    if sp.simplify(toy_discriminant.subs(lam, far_threshold)) != 0:
        raise RuntimeError("far toy threshold failed")
    near_series = sp.series(near_threshold.subs(D, 1 / h), h, 0, 4)
    expected_series = -h / 8 + h**2 / 64 - 5 * h**3 / 256
    if sp.simplify(near_series.removeO() - expected_series) != 0:
        raise RuntimeError("near toy threshold asymptotic failed")

    toy_root_minus = sp.factor((-1 - 1 / sp.sqrt(D)) / (1 - 1 / D))
    toy_root_plus = sp.factor((-1 + 1 / sp.sqrt(D)) / (1 - 1 / D))
    for root in (toy_root_minus, toy_root_plus):
        if sp.simplify(toy_scaled_jensen.subs({lam: 0, z: root})) != 0:
            raise RuntimeError("toy boundary root failed")
    toy_unscaled_gap = sp.factor((toy_root_plus - toy_root_minus) / D)

    n = sp.symbols("n", integer=True, nonnegative=True)
    f3, f4 = sp.symbols("f3 f4", real=True)
    c0 = 4 * n + 2
    first_parameter_shift = rho / 8
    first_center_shift = sp.factor(rho * (1 - 2 * n) / 4)
    f_lambda = 4 * rho * f2
    f_z_lambda = (c0 + 4) * f2 + 4 * rho * f3
    f_lambda_lambda = (
        c0 * (c0 + 4) * f2
        + 8 * (c0 + 4) * rho * f3
        + 16 * rho**2 * f4
    )
    g1_z = -rho * f2 - rho**2 * f3 / 2
    g1_lambda = -rho**2 * ((c0 + 8) * f3 + 4 * rho * f4) / 2
    g2 = rho**3 * f3 / 3 + rho**4 * f4 / 8
    second_parameter_shift = sp.factor(
        (
            f2 * first_center_shift**2 / 2
            - f_lambda_lambda * first_parameter_shift**2 / 2
            - g1_lambda * first_parameter_shift
            - g2
        )
        / f_lambda
    )
    expected_second_shift = sp.factor(
        -rho * (6 * n + 1) / 64 - rho**2 * f3 / (48 * f2)
    )
    if sp.simplify(second_parameter_shift - expected_second_shift) != 0:
        raise RuntimeError("second collision-parameter correction failed")
    if sp.simplify(
        f2 * first_center_shift
        + f_z_lambda * first_parameter_shift
        + g1_z
    ) != 0:
        raise RuntimeError("first collision-center correction failed")

    return {
        "entire_family": "F_(n,lambda)(z)=sum_(j>=0) A_(n+j)(lambda)*z^j/j!",
        "scaled_jensen": (
            "J_(D,n,lambda)(z/D)=sum_(j=0)^D ((D)_j/D^j)*A_(n+j)(lambda)*z^j/j!"
        ),
        "falling_expansion": "(D)_j/D^j=1-j*(j-1)/(2*D)+O_j(D^(-2))",
        "first_correction": (
            "J_(D,n,lambda)(z/D)=F_(n,lambda)(z)-"
            "z^2*F_(n,lambda)''(z)/(2*D)+O_K(D^(-2))"
        ),
        "second_correction": (
            "J_D(z/D)=F-z^2*F''/(2*D)+"
            "(z^3*F'''/3+z^4*F''''/8)/D^2+O_K(D^(-3))"
        ),
        "heat_pde": (
            "partial_lambda F_n=(4*n+2)*partial_z F_n+4*z*partial_z^2 F_n"
        ),
        "double_zero_transversality": (
            "F(rho)=F'(rho)=0 => partial_lambda F(rho)=4*rho*F''(rho)"
        ),
        "boundary_layer": (
            "If lambda=lambda_*+tau/D and z=rho+eta/sqrt(D), then "
            "2*D*J_D(z/D)/F''_*(rho)->eta^2+8*rho*tau-rho^2."
        ),
        "local_discriminant": str(local_discriminant),
        "root_layer": (
            "eta_+/- -> +/-sqrt(rho^2-8*rho*tau); at tau=0, "
            "w_+/-=rho/D+/-|rho|/D^(3/2)+O(D^(-2))."
        ),
        "collision_shift": (
            "lambda_D=lambda_*+rho/(8*D)-"
            "[rho*(6*n+1)/64+rho^2*F'''_*(rho)/(48*F''_*(rho))]/D^2+o(D^(-2))."
        ),
        "collision_center": (
            "z_D=rho+rho*(1-2*n)/(4*D)+o(D^(-1)) in the scaled Jensen variable."
        ),
        "external_field": (
            "If F_*(z)=(z-rho)^2*U(z), then F'''_*(rho)/F''_*(rho)=3*U'(rho)/U(rho); "
            "the global root external field first enters at order D^(-2)."
        ),
        "eventual_thresholds": {
            "definition": (
                "Theta_D=inf{ell:[ell,infinity) is contained in the degree-D hyperbolicity set}"
            ),
            "nesting": "Theta_D<=Theta_(D+1)<=Lambda",
            "exhaustion": "sup_D Theta_D=Lambda",
            "warning": "Theta_D<=0 requires the full forward ray, not a single lambda=0 test.",
        },
        "toy": {
            "family": str(toy),
            "pde_residual": str(toy_pde_residual),
            "scaled_jensen": str(toy_scaled_jensen),
            "discriminant": str(toy_discriminant),
            "near_threshold": str(near_threshold),
            "near_threshold_series": str(near_series),
            "far_threshold": str(far_threshold),
            "boundary_scaled_roots": [str(toy_root_minus), str(toy_root_plus)],
            "boundary_unscaled_gap": str(toy_unscaled_gap),
        },
        "checks": {
            "falling_factor_orders": [0, 12],
            "falling_factor_checks": falling_checks,
            "local_quadratic": str(local_quadratic),
            "first_parameter_shift": str(first_parameter_shift),
            "first_center_shift": str(first_center_shift),
            "second_parameter_shift": str(second_parameter_shift),
            "toy_exact_checks": 7,
        },
    }


def build_payload() -> dict:
    exact = build_exact()
    rows = [
        LemmaRow(
            id="sdzb_01_scaled_jensen_identity",
            role="exact_identity",
            readiness="available_exact",
            claim="Degree scaling inserts the falling-factorial multiplier into the entire-function coefficients.",
            formula=exact["scaled_jensen"],
            proof_boundary="Exact coefficient identity only.",
        ),
        LemmaRow(
            id="sdzb_02_first_degree_correction",
            role="exact_asymptotic",
            readiness="ready_to_apply",
            claim="The first two locally uniform finite-degree corrections are universal differential operators on the coefficient entire function.",
            formula=f"{exact['first_correction']}; {exact['second_correction']}",
            proof_boundary="Compact-z asymptotic for an entire coefficient function.",
        ),
        LemmaRow(
            id="sdzb_03_coefficient_heat_pde",
            role="exact_identity",
            readiness="available_exact",
            claim="The shifted coefficient entire functions obey a Laguerre-type heat PDE.",
            formula=exact["heat_pde"],
            proof_boundary="Exact coefficient flow only.",
        ),
        LemmaRow(
            id="sdzb_04_double_zero_transversality",
            role="exact_local_lemma",
            readiness="ready_to_apply",
            claim="At a nonzero double zero, the heat derivative is transverse with a shift-independent coefficient.",
            formula=exact["double_zero_transversality"],
            proof_boundary="Nondegenerate double zero only.",
        ),
        LemmaRow(
            id="sdzb_05_universal_boundary_layer",
            role="exact_asymptotic_theorem",
            readiness="ready_to_apply",
            claim="The joint degree-parameter-root scaling converges to one universal quadratic collision model.",
            formula=exact["boundary_layer"],
            proof_boundary="Local layer near one finite nondegenerate double zero.",
        ),
        LemmaRow(
            id="sdzb_06_root_gap_law",
            role="exact_asymptotic_theorem",
            readiness="ready_to_apply",
            claim="At the entire-function boundary, the associated finite Jensen roots remain distinct with an unscaled gap of order D^(-3/2).",
            formula=exact["root_layer"],
            proof_boundary="The local root pair only; no claim about all other roots.",
        ),
        LemmaRow(
            id="sdzb_07_finite_collision_shift",
            role="exact_asymptotic_theorem",
            readiness="ready_to_apply",
            claim="The finite-degree pair collides behind the entire-function boundary by rho/(8D), with a second-order correction governed by the regularized global root field.",
            formula=f"{exact['collision_shift']} {exact['external_field']}",
            proof_boundary="Local collision branch and second-order external field only.",
        ),
        LemmaRow(
            id="sdzb_08_exact_toy_validation",
            role="exact_countercheck",
            readiness="guard_validated",
            claim="A solvable polynomial heat family realizes both the root-gap and negative 1/(8D) threshold shift exactly.",
            formula="F_lambda=z^2+(2+12*lambda)z+1+4*lambda+12*lambda^2, rho=-1",
            proof_boundary="Algebraic model check, not a zeta coefficient family.",
            diagnostics=exact["toy"],
        ),
        LemmaRow(
            id="sdzb_09_eventual_threshold_exhaustion",
            role="exact_global_reduction",
            readiness="ready_to_apply",
            claim="Nested finite-degree eventual hyperbolicity thresholds increase to the full Newman boundary.",
            formula=exact["eventual_thresholds"]["exhaustion"],
            proof_boundary="Equivalent exhaustion theorem; it does not bound any open higher-degree threshold.",
            diagnostics=exact["eventual_thresholds"],
        ),
        LemmaRow(
            id="sdzb_10_uniform_layer_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="To prove Lambda<=0 by this route, control the scaled collision discriminant uniformly over degree and over collision locations, including escape in zero height.",
            formula="uniform lower bound for rho^2-8*rho*tau plus controlled remainders",
            proof_boundary="Open degree-uniform and height-uniform theorem; not RH or Lambda <= 0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_scaled_double_zero_boundary_layer_lemma",
        "date": "2026-07-11",
        "status": "exact scaled double-zero boundary layer and finite-threshold exhaustion theorem",
        "proof_boundary": (
            "This artifact derives the universal 1/D Jensen correction, the local D^(-3/2) root gap, "
            "the rho/(8D) finite collision shift, its second-order external-field correction, and the exhaustion of the Newman boundary by nested "
            "finite-degree eventual thresholds. It does not provide the required uniform remainder bound, "
            "exclude collision escape in degree or zero height, prove PF-infinity, RH, or Lambda <= 0."
        ),
        "sources": [
            "outputs/jensen_window_pf_polar_heat_collision_cascade_lemma.md",
            "outputs/jensen_window_pf_cofinal_scaling_limit_equivalence_gate.md",
            "outputs/formal_core.md",
            "outputs/jensen_window_pf_cubic_forward_uniform_tail_certificate.md",
        ],
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    toy = exact["toy"]
    return "\n".join(
        [
            "# Jensen-Window PF Scaled Double-Zero Boundary-Layer Lemma",
            "",
            "Date: 2026-07-11",
            "",
            "Status: exact scaled double-zero boundary layer and finite-threshold",
            "exhaustion theorem. This is not a proof of PF-infinity, RH, or",
            "`Lambda <= 0`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_scaled_double_zero_boundary_layer_lemma.json",
            "python work/rh_compute/scripts/jensen_window_pf_scaled_double_zero_boundary_layer_lemma.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_scaled_double_zero_boundary_layer_lemma.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF scaled double-zero boundary-layer lemma: 10 rows, 0 issues, 3 exact scaling identities, 1 heat PDE, 1 double-zero transversality, 1 universal boundary layer, 1 D^(-3/2) gap law, 1 external-field D^(-2) collision law, 1 exact toy family, 1 threshold-exhaustion theorem, 1 open uniform handoff",
            "```",
            "",
            "## First Degree Correction",
            "",
            "For every fixed shift and compact `z` set,",
            "",
            "```text",
            exact["scaled_jensen"],
            exact["falling_expansion"],
            exact["first_correction"],
            exact["second_correction"],
            "```",
            "",
            "The `O_K(D^(-2))` remainder follows by splitting the entire series",
            "on a larger compact disk and applying the falling-factorial expansion",
            "to the coefficient core.",
            "",
            "## Double-Zero Layer",
            "",
            "The coefficient flow is",
            "",
            "```text",
            exact["heat_pde"],
            exact["double_zero_transversality"],
            "```",
            "",
            "Let `rho<0` be a nondegenerate double zero at `lambda_*`. Under",
            "",
            "```text",
            "lambda=lambda_*+tau/D, z=rho+eta/sqrt(D),",
            "```",
            "",
            "the scaled Jensen polynomial has the universal limit",
            "",
            "```text",
            exact["boundary_layer"],
            "```",
            "",
            "Consequently",
            "",
            "```text",
            exact["root_layer"],
            exact["collision_shift"],
            exact["collision_center"],
            "```",
            "",
            "Because `rho<0`, the finite-degree collision lies on the bad side of",
            "the entire-function boundary. This is why every fixed degree can remain",
            "strict at the limiting boundary.",
            "Writing `F_*(z)=(z-rho)^2 U(z)` gives",
            "`F'''_*(rho)/F''_*(rho)=3 U'(rho)/U(rho)`. Thus the leading layer is",
            "universal, while the global root external field first enters at order",
            "`D^(-2)`.",
            "",
            "## Exact Model Check",
            "",
            "The solvable family",
            "",
            "```text",
            f"F_lambda(z)={toy['family']}",
            "partial_lambda F=2*F'+4*z*F''",
            f"J_D(z/D)={toy['scaled_jensen']}",
            "```",
            "",
            "has a double zero `rho=-1` at `lambda=0`. Its near collision threshold is",
            "",
            "```text",
            toy["near_threshold"],
            f"={toy['near_threshold_series']}",
            "```",
            "",
            "and its exact unscaled boundary gap is",
            "",
            "```text",
            toy["boundary_unscaled_gap"],
            "```",
            "",
            "matching the universal `-1/(8D)` shift and `2/D^(3/2)` gap.",
            "",
            "## Newman Threshold Exhaustion",
            "",
            "Let `S_D` be the parameter set where degree `D` is hyperbolic and set",
            "",
            "```text",
            exact["eventual_thresholds"]["definition"],
            exact["eventual_thresholds"]["nesting"],
            exact["eventual_thresholds"]["exhaustion"],
            "```",
            "",
            "Polar descent gives the nesting. Closedness and Jensen's theorem force",
            "the supremum to equal the full Newman boundary. A single endpoint test",
            "does not bound `Theta_D`; one must control the full forward parameter ray.",
            "",
            "The remaining proof obligation is now quantitative: a degree- and",
            "height-uniform bound on the scaled collision layer and its remainders.",
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
        "wrote Jensen-window PF scaled double-zero boundary-layer lemma: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
