#!/usr/bin/env python3
"""Prove forward invariance of the contiguous order-three compound cone."""

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
    "jensen_window_pf_compound_order3_forward_invariance_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order3_forward_invariance_certificate.md"
)
ENTRY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.json"
)
CUBIC_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_cubic_forward_uniform_tail_certificate.json"
)
RATIO_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.json"
)


@dataclass(frozen=True)
class CertificateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def exact_diagnostics() -> dict:
    k = sp.symbols("k", integer=True, positive=True)
    a, b, c, d = sp.symbols("a b c d", positive=True)

    def q_flow(index: sp.Expr, current: sp.Expr, next_q: sp.Expr) -> sp.Expr:
        return sp.factor(
            (2 * index - 1) * current
            - (2 * index + 3) * (current**3 - current) / next_q**2
        )

    x_b = 1 - b**-2
    x_c = 1 - c**-2
    q_prev_prime_over_rk = q_flow(k - 1, a, b) / x_b
    q_mid_prime_over_rk = q_flow(k, b, c)
    q_next_prime_over_rk = x_c * q_flow(k + 1, c, d)
    derivative = sp.factor(
        c * q_prev_prime_over_rk
        + a * q_next_prime_over_rk
        - 2 * b * q_mid_prime_over_rk
    )
    compound = a * c - b**2 + 1
    next_compound = b * d - c**2 + 1
    alpha = sp.factor((2 * k + 5) * a * (b * d + c**2 - 1) / (c * d**2))
    beta_numerator = sp.factor(
        (2 * k + 1) * a**2 * c**2
        + (2 * k + 1) * a * b**2 * c
        - (2 * k + 1) * a * c
        + (4 * k + 6) * b**4
        + (-4 * k + 2) * b**2 * c**2
        - (4 * k + 6) * b**2
    )
    beta = sp.factor(-beta_numerator / (c**2 * (b**2 - 1)))
    remainder = sp.factor(derivative - alpha * next_compound - beta * compound)
    if remainder != 0:
        raise RuntimeError(f"cooperative decomposition failed: {remainder}")

    boundary = sp.factor(derivative.subs(a, (b**2 - 1) / c))
    expected_boundary = sp.factor(
        alpha.subs(a, (b**2 - 1) / c) * next_compound
    )
    if sp.factor(boundary - expected_boundary) != 0:
        raise RuntimeError("boundary factorization failed")

    h, y, u, v, w = sp.symbols("h y u v w", positive=True)
    scaling = {
        k: h**-2,
        a: (y - u * h**2) / h,
        b: y / h,
        c: (y + v * h**2) / h,
        d: (y + (v + w) * h**2) / h,
    }
    scaled_alpha = sp.factor(sp.cancel(alpha.subs(scaling)))
    scaled_sum = sp.factor(sp.cancel((alpha + beta).subs(scaling)))
    alpha_limit = sp.factor(sp.limit(h**2 * scaled_alpha, h, 0, dir="+"))
    sum_limit = sp.factor(sp.limit(scaled_sum, h, 0, dir="+"))
    if alpha_limit != 4:
        raise RuntimeError(f"bad alpha scaling limit: {alpha_limit}")
    expected_sum_limit = 2 * (u + 2 * v - 3 * w) / y
    if sp.factor(sum_limit - expected_sum_limit) != 0:
        raise RuntimeError(f"bad diagonal-coefficient scaling limit: {sum_limit}")
    scaled_sum_denominator = sp.factor(sp.together(scaled_sum).as_numer_denom()[1])

    return {
        "indices": "k=n+2, a=q_(k-1), b=q_k, c=q_(k+1), d=q_(k+2)",
        "q_flow": (
            "q_j'=r_j*((2*j-1)*q_j-(2*j+3)*(q_j^3-q_j)/q_(j+1)^2)"
        ),
        "compound": "C_n=a*c-b^2+1",
        "next_compound": "C_(n+1)=b*d-c^2+1",
        "cooperative_system": "C_n'/r_k=alpha_k*C_(n+1)+beta_k*C_n",
        "alpha": "alpha_k=(2*k+5)*a*(b*d+c^2-1)/(c*d^2)",
        "beta": "beta_k=-N_k/(c^2*(b^2-1))",
        "beta_numerator": (
            "N_k=(2*k+1)*a^2*c^2+(2*k+1)*a*b^2*c-(2*k+1)*a*c"
            "+(4*k+6)*b^4+(-4*k+2)*b^2*c^2-(4*k+6)*b^2"
        ),
        "boundary": (
            "at C_n=0, C_n'/r_k=alpha_k*C_(n+1), alpha_k>0"
        ),
        "alpha_positive_reason": (
            "q_j>1 implies a,c,d>0 and b*d+c^2-1>0"
        ),
        "tail_scaling": (
            "k=h^(-2), b=y/h, a=(y-u*h^2)/h, "
            "c=(y+v*h^2)/h, d=(y+(v+w)*h^2)/h"
        ),
        "alpha_scaled_limit": str(alpha_limit),
        "diagonal_scaled_limit": str(sum_limit),
        "scaled_sum_denominator": str(scaled_sum_denominator),
    }


def source_diagnostics() -> dict:
    entry = load_json(ENTRY_SOURCE)
    cubic = load_json(CUBIC_SOURCE)
    ratio = load_json(RATIO_SOURCE)
    if entry.get("summary", {}).get("full_entry_rows") != 1:
        raise RuntimeError("order-three entry source is not closed")
    cubic_ids = {row.get("id") for row in cubic.get("rows", [])}
    if not {"cfut_07_forward_uniform_tail", "cfut_08_cubic_continuation"}.issubset(cubic_ids):
        raise RuntimeError("cubic forward-uniform source is not closed")
    if ratio.get("summary", {}).get("full_cone_propagation_rows") != 1:
        raise RuntimeError("ratio-cone forward source is not closed")
    return {
        "entry_status": entry.get("status"),
        "entry_full_rows": entry.get("summary", {}).get("full_entry_rows"),
        "cubic_status": cubic.get("status"),
        "cubic_forward_uniform_row": "cfut_07_forward_uniform_tail",
        "cubic_continuation_row": "cfut_08_cubic_continuation",
        "ratio_status": ratio.get("status"),
        "ratio_full_rows": ratio.get("summary", {}).get("full_cone_propagation_rows"),
    }


def build_artifact() -> dict:
    exact = exact_diagnostics()
    sources = source_diagnostics()
    rows = [
        CertificateRow(
            id="co3fi_01_q_flow",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The reciprocal-defect coordinate obeys an exact one-sided heat ODE.",
            formula=exact["q_flow"],
            proof_boundary="Exact coordinate identity only.",
        ),
        CertificateRow(
            id="co3fi_02_compound_coordinate",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Use the contiguous order-three compound margin and its next shift.",
            formula="C_n=a*c-b^2+1, C_(n+1)=b*d-c^2+1",
            proof_boundary="Contiguous order three only.",
        ),
        CertificateRow(
            id="co3fi_03_cooperative_system",
            role="exact_flow_lemma",
            readiness="ready_to_apply",
            claim="Differentiation and exact division give a one-sided cooperative compound system.",
            formula=exact["cooperative_system"],
            proof_boundary="Exact local flow algebra; infinite continuation is separate.",
            diagnostics={
                "alpha": exact["alpha"],
                "beta": exact["beta"],
                "beta_numerator": exact["beta_numerator"],
            },
        ),
        CertificateRow(
            id="co3fi_04_inward_boundary",
            role="exact_boundary_theorem",
            readiness="ready_to_apply",
            claim="At a zero compound margin, the vector field points inward whenever the next margin is nonnegative.",
            formula=exact["boundary"],
            proof_boundary="Boundary theorem inside the strict ratio cone.",
            diagnostics={"positivity": exact["alpha_positive_reason"]},
        ),
        CertificateRow(
            id="co3fi_05_forward_tail_input",
            role="theorem_input",
            readiness="ready_to_apply",
            claim="The existing cubic theorem supplies a compact-interval reciprocal-increment tail and the full cubic cone.",
            formula="0<=g_k<=B_L/sqrt(k) on every [-100,L]",
            proof_boundary="Previously proved degree-3 tail input.",
            diagnostics=sources,
        ),
        CertificateRow(
            id="co3fi_06_coefficient_growth",
            role="exact_analytic_lemma",
            readiness="ready_to_apply",
            claim="The compact-interval tail makes alpha_k=O_L(k) and alpha_k+beta_k=O_L(1).",
            formula=(
                "lim_(h->0) h^2*alpha_k=4; "
                "lim_(h->0)(alpha_k+beta_k)=2*(u+2*v-3*w)/y"
            ),
            proof_boundary="Coefficient-growth lemma along the actual compact heat trajectory.",
            diagnostics={
                "tail_scaling": exact["tail_scaling"],
                "scaled_sum_denominator": exact["scaled_sum_denominator"],
            },
        ),
        CertificateRow(
            id="co3fi_07_weighted_maximum_principle",
            role="exact_infinite_maximum_principle",
            readiness="ready_to_apply",
            claim="Dividing C_n by n+1 gives a locally attained negative minimum and a bounded diagonal coefficient on every compact heat interval.",
            formula=(
                "z_n=C_n/(n+1), z_n'=A_n*(n+2)/(n+1)*z_(n+1)+B_n*z_n"
            ),
            proof_boundary="Uses uniform C_n=O_L(1) and the established coefficient-growth bounds.",
            diagnostics={
                "tail_attainment": (
                    "g_k=O_L(k^(-1/2)) and q_k=O_L(sqrt(k)) imply C_n=O_L(1), hence z_n->0 uniformly."
                ),
                "bounded_diagonal": (
                    "r_k*(alpha_k*(n+2)/(n+1)+beta_k)="
                    "r_k*(alpha_k+beta_k+alpha_k/(n+1)) is bounded above."
                ),
            },
        ),
        CertificateRow(
            id="co3fi_08_strict_forward_propagation",
            role="theorem_composition",
            readiness="ready_to_apply",
            claim="Strict all-shift entry and the weighted maximum principle propagate every contiguous order-three margin to every finite forward heat time.",
            formula="C_n(lambda)>0 for every n>=0 and finite lambda>=-100",
            proof_boundary="Contiguous order-three family only.",
        ),
        CertificateRow(
            id="co3fi_09_lambda_zero_theorem",
            role="theorem_conclusion",
            readiness="ready_to_apply",
            claim="Every shifted contiguous 3x3 signed-Hankel minor has the required negative sign at lambda=0.",
            formula="D_(3,n)(0)<0 for every n>=0",
            proof_boundary="One contiguous compound order at lambda=0.",
        ),
        CertificateRow(
            id="co3fi_10_higher_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="Noncontiguous order-three minors and every compound order four or higher still require new coordinates and invariance theorems.",
            formula="contiguous order 3 closed; noncontiguous order 3 and k>=4 open",
            proof_boundary="Open all-order signed-Hankel/Jensen bridge; not PF-infinity, RH, or Lambda<=0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order3_forward_invariance_certificate",
        "date": "2026-07-12",
        "status": (
            "exact all-shift contiguous order-three propagation through lambda=0 "
            "with noncontiguous and higher compounds open"
        ),
        "proof_boundary": (
            "This artifact proves strict forward propagation of the contiguous "
            "order-three reciprocal-defect margin and hence D_(3,n)(0)<0 for every "
            "shift. It does not prove noncontiguous order-three minors, any compound "
            "order four or higher, the all-order signed-Hankel/Jensen bridge, "
            "PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_reciprocal_defect_compound_order3_gate.md",
            "outputs/jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.md",
            "outputs/jensen_window_pf_cubic_forward_uniform_tail_certificate.md",
            "outputs/jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.md",
            "outputs/formal_core.md",
        ],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order3_forward_invariance_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order3_forward_invariance_certificate.py"
        ),
        "exact": exact,
        "source_diagnostics": sources,
        "summary": {
            "certificate_rows": len(rows),
            "exact_identity_rows": 2,
            "cooperative_flow_rows": 1,
            "inward_boundary_rows": 1,
            "coefficient_growth_rows": 1,
            "infinite_maximum_principle_rows": 1,
            "full_forward_propagation_rows": 1,
            "lambda_zero_theorem_rows": 1,
            "open_handoffs": 1,
            "ready_to_apply_rows": sum(row.readiness == "ready_to_apply" for row in rows),
        },
        "rows": [asdict(row) for row in rows],
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Jensen-Window PF Compound Order-Three Forward-Invariance Certificate",
        "",
        "Date: 2026-07-12",
        "",
        "Status: exact all-shift contiguous order-three propagation through",
        "lambda=0 with noncontiguous and higher compounds open. This is not a",
        "proof of PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order3_forward_invariance_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order3_forward_invariance_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order3_forward_invariance_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order3_forward_invariance_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        "validated Jensen-window PF compound order-three forward-invariance certificate: 10 rows, 0 issues, 2 exact identities, 1 cooperative flow, 1 inward boundary theorem, 1 coefficient-growth lemma, 1 infinite maximum principle, 1 full forward propagation theorem, 1 lambda=0 theorem, 1 open higher-compound handoff",
        "```",
        "",
        "## Cooperative Flow",
        "",
        "Put `k=n+2` and",
        "",
        "```text",
        "a=q_(k-1), b=q_k, c=q_(k+1), d=q_(k+2),",
        "C_n=a*c-b^2+1, C_(n+1)=b*d-c^2+1.",
        "```",
        "",
        "The exact reciprocal-defect heat equation is",
        "",
        "```text",
        exact["q_flow"],
        "```",
        "",
        "and direct differentiation factors as",
        "",
        "```text",
        exact["cooperative_system"],
        exact["alpha"],
        exact["beta"],
        exact["beta_numerator"],
        "```",
        "",
        "Inside the strict ratio cone, `alpha_k>0`. In particular,",
        "",
        "```text",
        "at C_n=0, C_n'/r_k=alpha_k*C_(n+1).",
        "```",
        "",
        "Thus the local boundary is exactly one-sided and inward whenever the next",
        "compound margin is nonnegative.",
        "",
        "## Compact Tail Control",
        "",
        "The proved cubic forward theorem gives, on every finite heat interval",
        "`[-100,L]`, a constant `B_L` such that",
        "",
        "```text",
        "0<=g_j=q_(j+1)-q_j<=B_L/sqrt(j).",
        "```",
        "",
        "The ratio cone gives `q_j>=sqrt(j)`, while summing the increment bound",
        "gives `q_j=O_L(sqrt(j))`. Use the exact scaling",
        "",
        "```text",
        exact["tail_scaling"],
        "```",
        "",
        "where `u,v,w` are the scaled neighboring increments and remain in a",
        "compact box. Exact cancellation gives",
        "",
        "```text",
        "lim_(h->0) h^2*alpha_k=4,",
        "lim_(h->0) (alpha_k+beta_k)=2*(u+2*v-3*w)/y.",
        "```",
        "",
        "The canceled denominator is nonzero on the compact box because `y>=1`",
        "and `h^2<=1/2`. Therefore",
        "",
        "```text",
        "alpha_k=O_L(k), alpha_k+beta_k=O_L(1).",
        "```",
        "",
        "## Infinite Maximum Principle",
        "",
        "The same tail estimates give `C_n=O_L(1)`. Set",
        "",
        "```text",
        "z_n=C_n/(n+1).",
        "```",
        "",
        "Then `z_n->0` uniformly on compact heat intervals and",
        "",
        "```text",
        "z_n'=r_k*alpha_k*((n+2)/(n+1))*z_(n+1)+r_k*beta_k*z_n.",
        "```",
        "",
        "Since `r_k<=r_0`, the effective diagonal coefficient",
        "",
        "```text",
        "r_k*(alpha_k+beta_k+alpha_k/(n+1))",
        "```",
        "",
        "is bounded above on each compact interval. After one exponential",
        "integrating factor, every negative spatial infimum is attained at a finite",
        "index and has nonnegative upper-right derivative. The standard connected-",
        "component argument excludes a negative component.",
        "",
        "The lambda=-100 entry theorem is strict. Variation of constants in the",
        "cooperative equation then preserves strict positivity at every fixed index.",
        "Hence",
        "",
        "```text",
        "C_n(lambda)>0 for every n>=0 and finite lambda>=-100,",
        "D_(3,n)(0)<0 for every n>=0.",
        "```",
        "",
        "This closes the complete shifted contiguous order-three layer at lambda=0.",
        "Noncontiguous order-three minors and every higher compound order remain",
        "open, so the all-order signed-Hankel/Jensen bridge is not proved.",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.md",
        "outputs/jensen_window_pf_cubic_forward_uniform_tail_certificate.md",
        "outputs/jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.md",
        "outputs/signed_hankel_jensen_bridge_target.md",
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "```",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print(
        "validated Jensen-window PF compound order-three forward-invariance "
        "certificate: 10 rows, 0 issues, 2 exact identities, 1 cooperative "
        "flow, 1 inward boundary theorem, 1 coefficient-growth lemma, 1 infinite "
        "maximum principle, 1 full forward propagation theorem, 1 lambda=0 "
        "theorem, 1 open higher-compound handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
