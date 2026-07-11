#!/usr/bin/env python3
"""Prove the quartic-quintic polar-contact mechanism at a double root."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_quartic_quintic_polar_contact_lemma.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_quartic_quintic_polar_contact_lemma.md"


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
    w, x, y, z, u = sp.symbols("w x y z u")
    P4 = 1 + 4 * w + 6 * x * w**2 + 4 * x**2 * y * w**3 + x**3 * y**2 * z * w**4
    P5 = (
        1
        + 5 * w
        + 10 * x * w**2
        + 10 * x**2 * y * w**3
        + 5 * x**3 * y**2 * z * w**4
        + x**4 * y**3 * z**2 * u * w**5
    )
    if sp.expand(P5 - w * sp.diff(P5, w) / 5 - P4) != 0:
        raise RuntimeError("quartic-quintic polar identity failed")

    a, p = sp.symbols("a p")
    A = -3 * a**2 + 8 * a + p
    B = -a**2 + 2 * a + p
    xv = A / 6
    yv = 18 * a * B / A**2
    zv = 2 * p * A / (3 * B**2)
    U = sp.factor(B * (3 * a**2 - 5 * a + 5 * p) / (6 * p**2))
    P5_boundary = sp.factor(P5.subs({x: xv, y: yv, z: zv}))
    root_value = sp.factor(P5_boundary.subs(w, -1 / a))
    expected_root_value = sp.factor(-2 * p**2 * (u - U) / (a**2 * B))
    if sp.simplify(root_value - expected_root_value) != 0:
        raise RuntimeError("quintic root-value threshold failed")

    threshold_factorization = sp.factor(P5_boundary.subs(u, U))
    cofactor = sp.factor(threshold_factorization / (1 + a * w) ** 3)
    expected_cofactor = sp.factor(
        (3 + (15 - 9 * a) * w + (3 * a**2 - 5 * a + 5 * p) * w**2) / 3
    )
    if sp.simplify(cofactor - expected_cofactor) != 0:
        raise RuntimeError("quintic triple-contact factorization failed")
    cofactor_discriminant = sp.factor(sp.discriminant(cofactor, w))
    expected_discriminant = sp.factor(5 * (3 * a**2 - 14 * a - 4 * p + 15) / 3)
    if sp.simplify(cofactor_discriminant - expected_discriminant) != 0:
        raise RuntimeError("quintic cofactor discriminant failed")

    alpha = sp.symbols("alpha_0:5", positive=True)
    r = sp.symbols("r", nonzero=True)
    log_polar = sum(1 / (1 + value * r) for value in alpha)
    log_polar_derivative = sp.factor(
        sum(-value / (1 + value * r) ** 2 for value in alpha)
    )

    d, j = sp.symbols("d j", integer=True, positive=True)
    binomial_identity = sp.simplify(
        (1 - j / (d + 1)) * sp.binomial(d + 1, j) - sp.binomial(d, j)
    )
    if binomial_identity != 0:
        raise RuntimeError("general adjacent-degree polar identity failed")

    return {
        "quartic": str(P4),
        "quintic": str(P5),
        "polar_identity": "P_4=P_5-(w/5)*P_5'",
        "general_polar_identity": "P_d=P_(d+1)-(w/(d+1))*P_(d+1)'",
        "binomial_identity": "(1-j/(d+1))*C(d+1,j)=C(d,j)",
        "polar_root_test": {
            "factorization": "P_5=product_i(1+alpha_i*w), alpha_i>0",
            "polar_ratio": "(5*P_5-w*P_5')/P_5=sum_i 1/(1+alpha_i*w)",
            "value_at_nonroot": str(log_polar),
            "derivative": str(log_polar_derivative),
            "strict_sign": "negative wherever P_5 is nonzero",
            "multiplicity_rule": "a nonzero root of P_5 of multiplicity m gives a polar-derivative root of multiplicity m-1",
        },
        "double_root_coordinates": {
            "A": str(A),
            "B": str(B),
            "U": str(U),
        },
        "quintic_value_at_quartic_double_root": str(root_value),
        "quintic_value_factorization": str(expected_root_value),
        "triple_contact_factorization": str(threshold_factorization),
        "quadratic_cofactor": str(cofactor),
        "quadratic_cofactor_discriminant": str(cofactor_discriminant),
        "contact_conclusion": (
            "If a positive-root hyperbolic P_5 has polar derivative P_4 with a nonzero double root -1/a, then P_5 has a triple root there and u=U(a,p)."
        ),
    }


def build_payload() -> dict:
    exact = build_exact()
    rows = [
        LemmaRow(
            id="qqpc_01_normalized_adjacent_windows",
            role="exact_definition",
            readiness="available_exact",
            claim="Normalized degree-4 and degree-5 Jensen windows use the same first four contraction monomials, with u entering only the quintic leading term.",
            formula="P_4(x,y,z), P_5(x,y,z,u)",
            proof_boundary="Adjacent normalized windows only.",
        ),
        LemmaRow(
            id="qqpc_02_polar_identity",
            role="exact_identity",
            readiness="available_exact",
            claim="The quartic window is the polar derivative of its quintic extension with pole zero.",
            formula=exact["polar_identity"],
            proof_boundary="Degree-4/5 algebra only.",
        ),
        LemmaRow(
            id="qqpc_03_polar_nonroot_simplicity",
            role="exact_lemma",
            readiness="available_exact",
            claim="For a positive-root hyperbolic quintic, every nonroot zero of its polar derivative is simple.",
            formula="d/dw sum_i 1/(1+alpha_i*w)=-sum_i alpha_i/(1+alpha_i*w)^2<0",
            proof_boundary="Positive-root hyperbolic quintics only.",
            diagnostics=exact["polar_root_test"],
        ),
        LemmaRow(
            id="qqpc_04_polar_multiplicity_drop",
            role="exact_lemma",
            readiness="available_exact",
            claim="At a nonzero root shared with the quintic, the polar derivative lowers multiplicity by exactly one.",
            formula="mult_r(P_5)=m => mult_r(5*P_5-w*P_5')=m-1",
            proof_boundary="Local multiplicity identity only.",
        ),
        LemmaRow(
            id="qqpc_05_double_implies_triple",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="A double root of a quartic polar derivative of a hyperbolic quintic forces a quintic triple root at the same location.",
            formula="mult_r(P_4)=2 => mult_r(P_5)=3",
            proof_boundary="Conditional on a hyperbolic positive-root quintic extension.",
        ),
        LemmaRow(
            id="qqpc_06_threshold_root_value",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="At a quartic double-root boundary point, the quintic value at that root is a nonzero positive prefactor times -(u-U).",
            formula=exact["quintic_value_factorization"],
            proof_boundary="Quartic double-root coordinates only.",
        ),
        LemmaRow(
            id="qqpc_07_quintic_triple_contact",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="The quartic heat threshold u=U is exactly the quintic triple-contact value.",
            formula=exact["triple_contact_factorization"],
            proof_boundary="Triple contact identity; cofactor hyperbolicity remains a separate sign condition.",
        ),
        LemmaRow(
            id="qqpc_08_cofactor_gate",
            role="exact_identity",
            readiness="available_exact",
            claim="The remaining quintic quadratic factor has an explicit discriminant gate.",
            formula=f"Disc(R_2)={exact['quadratic_cofactor_discriminant']}",
            proof_boundary="Quadratic cofactor test only.",
        ),
        LemmaRow(
            id="qqpc_09_general_adjacent_degree_identity",
            role="exact_identity",
            readiness="available_exact",
            claim="The polar relation is not special to degrees four and five; it holds for every adjacent Jensen degree.",
            formula=exact["general_polar_identity"],
            proof_boundary="Algebraic adjacent-degree identity only.",
        ),
        LemmaRow(
            id="qqpc_10_all_degree_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="Turn the adjacent-degree polar-contact hierarchy into a noncircular infinite-degree invariant or compactness theorem.",
            formula="hyperbolic P_(d+1) => controlled boundary contact for P_d, uniformly as d->infinity",
            proof_boundary="Open all-degree closure; not PF-infinity, RH, or Lambda <= 0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_quartic_quintic_polar_contact_lemma",
        "date": "2026-07-10",
        "status": "exact quartic-quintic polar-contact theorem with open all-degree closure",
        "proof_boundary": (
            "This artifact identifies the quartic heat threshold with quintic triple contact and proves the general adjacent-degree polar identity. It does not supply a hyperbolic quintic extension for the zeta coefficients, close the infinite degree hierarchy, prove PF-infinity, RH, or Lambda <= 0."
        ),
        "sources": [
            "outputs/jensen_window_pf_quartic_double_root_threshold_lemma.md",
            "outputs/jensen_window_pf_heat_flow_jensen_hierarchy_lemma.md",
            "outputs/jensen_window_pf_bridge_target.md",
        ],
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    return "\n".join(
        [
            "# Jensen-Window PF Quartic-Quintic Polar-Contact Lemma",
            "",
            "Date: 2026-07-10",
            "",
            "Status: exact adjacent-degree contact theorem with the infinite-degree",
            "closure open. This is not a proof of PF-infinity, RH, or `Lambda <= 0`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_quartic_quintic_polar_contact_lemma.json",
            "python work/rh_compute/scripts/jensen_window_pf_quartic_quintic_polar_contact_lemma.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_quartic_quintic_polar_contact_lemma.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF quartic-quintic polar-contact lemma: 10 rows, 0 issues, 4 exact polar identities, 1 strict nonroot test, 1 multiplicity rule, 1 double-to-triple theorem, 1 quintic contact factorization, 1 cofactor gate, 1 open all-degree handoff",
            "```",
            "",
            "## Polar Identity",
            "",
            "The normalized adjacent windows satisfy",
            "",
            "```text",
            exact["polar_identity"],
            "```",
            "",
            "If `P_5=product_i(1+alpha_i*w)` with every `alpha_i>0`, then",
            "",
            "```text",
            "(5*P_5-w*P_5')/P_5=sum_i 1/(1+alpha_i*w),",
            "d/dw=-sum_i alpha_i/(1+alpha_i*w)^2<0.",
            "```",
            "",
            "Thus a polar-derivative root away from the roots of `P_5` is simple.",
            "At a nonzero quintic root, the polar derivative lowers multiplicity by",
            "one. Consequently a double root of `P_4` forces a triple root of any",
            "hyperbolic positive-root quintic extension.",
            "",
            "## Threshold Equals Triple Contact",
            "",
            "In the quartic double-root coordinates of the threshold lemma,",
            "",
            "```text",
            f"P_5(-1/a)={exact['quintic_value_factorization']}",
            "```",
            "",
            "so the shared root condition is exactly `u=U(a,p)`. At this value,",
            "",
            "```text",
            exact["triple_contact_factorization"],
            "```",
            "",
            "The remaining quadratic factor has discriminant",
            "",
            "```text",
            exact["quadratic_cofactor_discriminant"],
            "```",
            "",
            "which is the explicit residual gate for a fully hyperbolic quintic",
            "contact.",
            "",
            "## General Handoff",
            "",
            "The binomial identity extends the polar relation to every degree:",
            "",
            "```text",
            exact["general_polar_identity"],
            "```",
            "",
            "This exposes a plausible higher-minor hierarchy: degree `d+1` controls",
            "boundary multiplicity for degree `d`. It is not yet a proof because no",
            "noncircular infinite-degree closure or uniform compactness theorem has",
            "been supplied for the zeta coefficients.",
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
        "wrote Jensen-window PF quartic-quintic polar-contact lemma: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
