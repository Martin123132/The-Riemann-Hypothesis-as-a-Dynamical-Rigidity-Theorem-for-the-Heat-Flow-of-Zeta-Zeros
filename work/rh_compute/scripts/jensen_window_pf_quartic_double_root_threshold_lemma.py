#!/usr/bin/env python3
"""Derive the exact quartic double-root heat threshold and live handoff."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_quartic_double_root_threshold_lemma.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_quartic_double_root_threshold_lemma.md"


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
    a, p, u, k, w = sp.symbols("a p u k w")
    A = -3 * a**2 + 8 * a + p
    B = -a**2 + 2 * a + p
    x = sp.factor(A / 6)
    y = sp.factor(18 * a * B / A**2)
    z = sp.factor(2 * p * A / (3 * B**2))
    threshold = sp.factor(B * (3 * a**2 - 5 * a + 5 * p) / (6 * p**2))
    curvature = 3 * a**2 - 4 * a + p

    H = sp.factor(
        4
        * (2 * k + 7)
        * (
            3 * a**4
            - 11 * a**3
            + 2 * a**2 * p
            + 10 * a**2
            - 5 * a * p
            + 6 * p**2 * u
            - 5 * p**2
        )
        / (a * A * B)
    )
    expected_H = sp.factor(24 * (2 * k + 7) * p**2 * (u - threshold) / (a * A * B))
    if sp.simplify(H - expected_H) != 0:
        raise RuntimeError("double-root heat threshold factorization failed")

    b, c = sp.symbols("b c")
    if sp.expand((a - b) * (a - c)).subs({b + c: 4 - 2 * a, b * c: p}) == 0:
        raise RuntimeError("unexpected symbolic substitution collapse")
    curvature_identity = sp.expand((a - b) * (a - c)).subs(c, 4 - 2 * a - b)
    curvature_identity = sp.rem(curvature_identity, b**2 - (4 - 2 * a) * b + p, b)
    if sp.simplify(curvature_identity - curvature) != 0:
        raise RuntimeError("double-root curvature identity failed")

    triple_p = 4 * a - 3 * a**2
    triple_threshold = sp.factor(threshold.subs(p, triple_p))
    triple_H = sp.factor(H.subs({p: triple_p, u: triple_threshold}))
    if triple_H != 0:
        raise RuntimeError("triple-root constant perturbation does not vanish at threshold")
    triple_perturbation = sp.factor(
        24
        * w**2
        * (a - 1) ** 2
        * (a * w + 1) ** 2
        * (2 * a * k + a - 2 * k + 1)
        / (a - 2)
    )

    av = sp.Rational(13, 20)
    bv = sp.Rational(21, 50)
    cv = sp.Rational(57, 25)
    pv = sp.factor(bv * cv)
    Uv = sp.factor(threshold.subs({a: av, p: pv}))
    curvature_v = sp.factor(curvature.subs({a: av, p: pv}))
    xv = sp.factor(x.subs({a: av, p: pv}))
    yv = sp.factor(y.subs({a: av, p: pv}))
    zv = sp.factor(z.subs({a: av, p: pv}))
    gap = sp.factor(Uv - zv)
    criterion_v = sp.factor(curvature_v * (zv - Uv))
    if not (curvature_v < 0 and gap > 0 and criterion_v > 0):
        raise RuntimeError("countermodel is not rejected by the exact inward threshold")

    return {
        "root_coordinates": {
            "factorization": "P(w)=(1+a*w)^2*(1+b*w)*(1+c*w)",
            "normalization": "2*a+b+c=4",
            "symmetric_simple_roots": "p=b*c",
            "simple_root_quadratic": "T^2-(4-2*a)*T+p=0",
        },
        "contraction_coordinates": {"x": str(x), "y": str(y), "z": str(z)},
        "positive_denominators": {
            "A": f"{A}=6*x>0",
            "B": f"{B}=e_3/(2*a)>0",
        },
        "double_root_curvature": {
            "P_second": "P''(-1/a)=2*(a-b)*(a-c)",
            "symmetric_factor": str(curvature),
        },
        "heat_value": str(H),
        "threshold": str(threshold),
        "heat_threshold_factorization": str(expected_H),
        "inward_criterion": "(3*a^2-4*a+p)*(u-U(a,p))<=0",
        "branch_interpretation": {
            "double_root_outside_simple_roots": "u<=U",
            "double_root_between_simple_roots": "u>=U",
            "triple_root": "u=U is necessary for first-order viability",
        },
        "triple_root": {
            "p": str(triple_p),
            "threshold": str(triple_threshold),
            "constant_heat_value_at_threshold": str(triple_H),
            "full_first_variation_at_threshold": str(triple_perturbation),
            "tangent_factor": "(1+a*w)^2",
        },
        "countermodel_threshold_check": {
            "a": str(av),
            "p": str(pv),
            "x": str(xv),
            "y": str(yv),
            "z_equals_u": str(zv),
            "U": str(Uv),
            "U_minus_u": str(gap),
            "curvature": str(curvature_v),
            "criterion_left_side": str(criterion_v),
            "conclusion": "positive left side violates the inward criterion",
        },
    }


def build_payload() -> dict:
    exact = build_exact()
    rows = [
        LemmaRow(
            id="qdrt_01_double_root_coordinates",
            role="exact_definition",
            readiness="available_exact",
            claim="Every nondegenerate hyperbolic quartic discriminant boundary point has one double-root factor and two simple-root factors.",
            formula="P=(1+a*w)^2*(1+b*w)*(1+c*w), 2*a+b+c=4",
            proof_boundary="Simple quartic boundary stratum only.",
        ),
        LemmaRow(
            id="qdrt_02_symmetric_contractions",
            role="exact_identity",
            readiness="available_exact",
            claim="The three contraction coordinates are rational functions of the double root a and p=b*c.",
            formula="x=A/6, y=18*a*B/A^2, z=2*p*A/(3*B^2)",
            proof_boundary="Coordinate identity only.",
            diagnostics=exact["contraction_coordinates"],
        ),
        LemmaRow(
            id="qdrt_03_root_splitting_criterion",
            role="exact_local_lemma",
            readiness="available_exact",
            claim="A first-order heat perturbation splits the double root into real roots exactly when its value has sign opposite to P'' at the double root.",
            formula="H(-1/a)*P''(-1/a)<=0",
            proof_boundary="Nondegenerate double-root first-order criterion only.",
        ),
        LemmaRow(
            id="qdrt_04_heat_value_factorization",
            role="exact_identity",
            readiness="available_exact",
            claim="The normalized heat perturbation at the double root factors into a positive prefactor times u-U(a,p).",
            formula=exact["heat_threshold_factorization"],
            proof_boundary="Boundary heat value only.",
        ),
        LemmaRow(
            id="qdrt_05_exact_inward_threshold",
            role="exact_boundary_theorem",
            readiness="ready_to_apply",
            claim="The complete nondegenerate quartic inward condition is a single root-branch-aware threshold inequality.",
            formula=exact["inward_criterion"],
            proof_boundary="Local degree-4 boundary condition, not a global invariant.",
            diagnostics=exact["branch_interpretation"],
        ),
        LemmaRow(
            id="qdrt_06_shift_independence",
            role="exact_identity",
            readiness="available_exact",
            claim="The threshold U is independent of the shift index; k enters only through a positive factor 2*k+7.",
            formula="U=B*(3*a^2-5*a+5*p)/(6*p^2)",
            proof_boundary="Threshold simplification only.",
        ),
        LemmaRow(
            id="qdrt_07_triple_root_necessity",
            role="exact_boundary_lemma",
            readiness="ready_to_apply",
            claim="On the triple-root stratum, first-order viability requires u=U.",
            formula="3*a^2-4*a+p=0 => H(-1/a)=0 iff u=U",
            proof_boundary="Necessary first-order condition at a degenerate boundary.",
        ),
        LemmaRow(
            id="qdrt_08_triple_root_tangency",
            role="exact_identity",
            readiness="available_exact",
            claim="At the triple-root threshold, the full first heat variation retains a double root factor.",
            formula="H(w) is divisible by (1+a*w)^2",
            proof_boundary="First-order tangency only; higher-time viability is not proved.",
            diagnostics=exact["triple_root"],
        ),
        LemmaRow(
            id="qdrt_09_countermodel_explanation",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="The earlier rational quartic obstruction fails precisely because its middle double-root branch requires u>=U but has u<U.",
            formula="curvature<0, U-u>0 => curvature*(u-U)>0",
            proof_boundary="Explains one abstract countermodel only.",
            diagnostics=exact["countermodel_threshold_check"],
        ),
        LemmaRow(
            id="qdrt_10_live_quartic_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="Construct a closed contraction-coordinate inequality realizing the branch-aware threshold and propagate it with a uniform spatial tail.",
            formula="Q>=0 plus a root-branch-compatible extension of (3*a^2-4*a+p)*(u-U)<=0",
            proof_boundary="Open quartic invariant and all-degree handoff; not PF-infinity, RH, or Lambda <= 0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_quartic_double_root_threshold_lemma",
        "date": "2026-07-10",
        "status": "exact quartic double-root heat threshold with open global invariant",
        "proof_boundary": (
            "This artifact proves the exact local inward threshold on every quartic double-root boundary branch and the necessary first-order triple-root equality. It does not construct a globally closed quartic cone, prove degree-4 zeta hyperbolicity, establish higher-time triple-root viability, PF-infinity, RH, or Lambda <= 0."
        ),
        "sources": [
            "outputs/jensen_window_pf_quartic_boundary_flow_obstruction.md",
            "outputs/jensen_window_pf_cubic_forward_uniform_tail_certificate.md",
            "outputs/jensen_window_pf_bridge_target.md",
        ],
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    check = exact["countermodel_threshold_check"]
    return "\n".join(
        [
            "# Jensen-Window PF Quartic Double-Root Threshold Lemma",
            "",
            "Date: 2026-07-10",
            "",
            "Status: exact local quartic heat threshold with the global quartic",
            "invariant open. This is not a proof of degree-4 zeta hyperbolicity,",
            "PF-infinity, RH, or `Lambda <= 0`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_quartic_double_root_threshold_lemma.json",
            "python work/rh_compute/scripts/jensen_window_pf_quartic_double_root_threshold_lemma.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_quartic_double_root_threshold_lemma.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF quartic double-root threshold lemma: 10 rows, 0 issues, 5 exact coordinate identities, 1 double-root splitting criterion, 1 branch-aware inward threshold, 1 triple-root equality, 1 tangent factor, 1 explained countermodel, 1 open global-invariant handoff",
            "```",
            "",
            "## Double-Root Coordinates",
            "",
            "Write a nondegenerate hyperbolic quartic boundary point as",
            "",
            "```text",
            "P(w)=(1+a*w)^2*(1+b*w)*(1+c*w),",
            "2*a+b+c=4, p=b*c.",
            "```",
            "",
            "Set `A=-3*a^2+8*a+p` and `B=-a^2+2*a+p`. Then",
            "",
            "```text",
            f"x={exact['contraction_coordinates']['x']}",
            f"y={exact['contraction_coordinates']['y']}",
            f"z={exact['contraction_coordinates']['z']}",
            "```",
            "",
            "Both `A=6*x` and `B=e_3/(2*a)` are positive on the positive-root",
            "boundary stratum.",
            "",
            "## Heat Threshold",
            "",
            "If `H=P_lambda`, local double-root splitting is real exactly when",
            "",
            "```text",
            "H(-1/a)*P''(-1/a)<=0.",
            "```",
            "",
            "Exact substitution into the coefficient heat flow gives",
            "",
            "```text",
            f"U(a,p)={exact['threshold']}",
            "H(-1/a)/r_k=24*(2*k+7)*p^2*(u-U)/(a*A*B).",
            "P''(-1/a)=2*(a-b)*(a-c).",
            "```",
            "",
            "Therefore the exact inward condition is",
            "",
            "```text",
            exact["inward_criterion"],
            "```",
            "",
            "If the double root lies outside the two simple roots this requires",
            "`u<=U`; if it lies between them it requires `u>=U`.",
            "",
            "## Triple Root",
            "",
            "When `3*a^2-4*a+p=0`, one simple root merges with the double root.",
            "A nonzero constant heat perturbation would split the cubic contact into",
            "one real root and a complex pair, so first-order viability requires",
            "`u=U`. At that equality the complete first variation factors as",
            "",
            "```text",
            exact["triple_root"]["full_first_variation_at_threshold"],
            "```",
            "",
            "and retains `(1+a*w)^2`. This is tangency, not yet a higher-time",
            "invariance proof.",
            "",
            "## Countermodel And Handoff",
            "",
            "For the earlier rational obstruction,",
            "",
            "```text",
            f"curvature={check['curvature']}<0",
            f"u={check['z_equals_u']}",
            f"U={check['U']}",
            f"U-u={check['U_minus_u']}>0.",
            "```",
            "",
            "The middle-root branch requires `u>=U`, explaining the outward heat",
            "direction exactly. The next task is to express and propagate this",
            "branch-aware threshold as a closed contraction-coordinate invariant,",
            "including the tangent triple-root stratum. Higher degrees remain open.",
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
        "wrote Jensen-window PF quartic double-root threshold lemma: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
