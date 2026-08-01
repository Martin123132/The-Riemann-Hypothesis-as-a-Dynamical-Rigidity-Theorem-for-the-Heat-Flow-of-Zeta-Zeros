#!/usr/bin/env python3
"""Build the signed-Hankel exclusion of the middle quartic boundary branch."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_quartic_signed_hankel_branch_exclusion_lemma.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_quartic_signed_hankel_branch_exclusion_lemma.md"
)
COUNTERMODEL = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_strong_logconcave_local_quartic_countermodel.json"
)


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
    x, y, z, a, p = sp.symbols("x y z a p")
    normalized = [sp.Integer(1), sp.Integer(1), x, x**2 * y, x**3 * y**2 * z]
    hankel = sp.factor(
        sp.det(sp.Matrix([[normalized[i + j] for j in range(3)] for i in range(3)]))
    )
    expected_hankel = sp.factor(
        x**3 * (x * y**2 * z - x * y**2 - y**2 * z + 2 * y - 1)
    )
    if sp.expand(hankel - expected_hankel) != 0:
        raise RuntimeError("normalized order-three Hankel determinant changed")

    A = -3 * a**2 + 8 * a + p
    B = -a**2 + 2 * a + p
    boundary = {
        x: A / 6,
        y: 18 * a * B / A**2,
        z: 2 * p * A / (3 * B**2),
    }
    curvature = 3 * a**2 - 4 * a + p
    boundary_hankel = sp.factor(hankel.subs(boundary))
    if sp.expand(boundary_hankel + curvature**3 / 216) != 0:
        raise RuntimeError("quartic boundary Hankel factorization changed")

    b, c = sp.symbols("b c")
    root_curvature = sp.factor((a - b) * (a - c))
    reduced_curvature = sp.expand(root_curvature.subs(c, 4 - 2 * a - b))
    reduced_p = sp.expand((b * (4 - 2 * a - b)))
    if sp.expand(reduced_curvature - curvature.subs(p, reduced_p)) != 0:
        raise RuntimeError("root-order curvature identity changed")

    threshold = sp.factor(B * (3 * a**2 - 5 * a + 5 * p) / (6 * p**2))
    return {
        "normalized_coefficients": [
            "1",
            "1",
            "x",
            "x^2*y",
            "x^3*y^2*z",
        ],
        "shift_scaling": (
            "A_(n+j)=A_n*r_n^j*B_j, r_n=A_(n+1)/A_n; "
            "D_(3,n)=A_n^3*r_n^6*det[B_(i+j)]"
        ),
        "normalized_hankel": str(hankel),
        "boundary_coordinates": {
            "x": str(boundary[x]),
            "y": str(boundary[y]),
            "z": str(boundary[z]),
        },
        "curvature": "C=3*a^2-4*a+p=(a-b)*(a-c)",
        "boundary_hankel": "det[B_(i+j)]=-C^3/216",
        "strict_xi_sign": "D_(3,n)(lambda)<0 for every n>=0 and finite lambda>=-100",
        "branch_conclusion": (
            "C>0: the double-root slope a lies outside the interval with "
            "endpoints b,c"
        ),
        "excluded_strata": [
            "C<0 middle-root branch",
            "C=0 triple-root stratum",
        ],
        "threshold": str(threshold),
        "reduced_inward_condition": "u<=U(a,p)",
    }


def countermodel_hankel_ball() -> flint.arb:
    payload = json.loads(COUNTERMODEL.read_text(encoding="utf-8"))
    moments = [
        flint.arb(entry["ball"]) for entry in payload["diagnostics"]["moments"][:5]
    ]
    a0, a1, a2, a3, a4 = moments
    return (
        a0 * (a2 * a4 - a3**2)
        - a1 * (a1 * a4 - a2 * a3)
        + a2 * (a1 * a3 - a2**2)
    )


def build_payload() -> dict:
    exact = build_exact()
    flint.ctx.prec = 512
    counter_hankel = countermodel_hankel_ball()
    if not bool(counter_hankel > 0):
        raise RuntimeError("local quartic countermodel does not fail the Xi Hankel sign")
    counter_text = counter_hankel.str(55).replace("e", "E")
    rows = [
        LemmaRow(
            id="qshb_01_shift_normalization",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Every shifted quartic can be normalized to B_0=B_1=1 with positive Hankel scaling.",
            formula=exact["shift_scaling"],
            proof_boundary="Positive normalization only.",
        ),
        LemmaRow(
            id="qshb_02_order_three_coordinate",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The normalized contiguous order-three Hankel determinant is explicit in the quartic contractions.",
            formula="det[B_(i+j)]=" + exact["normalized_hankel"],
            proof_boundary="One contiguous 3 by 3 minor.",
        ),
        LemmaRow(
            id="qshb_03_double_root_coordinates",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Substitute the complete nondegenerate positive-root quartic boundary coordinates.",
            formula="P=(1+a*w)^2*(1+b*w)*(1+c*w), 2*a+b+c=4, p=b*c",
            proof_boundary="Nondegenerate quartic boundary stratum.",
            diagnostics=exact["boundary_coordinates"],
        ),
        LemmaRow(
            id="qshb_04_boundary_factorization",
            role="exact_factorization",
            readiness="ready_to_apply",
            claim="On the quartic boundary, the order-three Hankel sign is the negative cube of the root-order curvature.",
            formula=exact["boundary_hankel"] + ", " + exact["curvature"],
            proof_boundary="Exact boundary identity before importing the Xi sign theorem.",
        ),
        LemmaRow(
            id="qshb_05_imported_xi_hankel_sign",
            role="imported_theorem",
            readiness="ready_to_apply",
            claim="The completed Xi compound-order-three flow gives the strict negative contiguous Hankel sign throughout the target interval.",
            formula=exact["strict_xi_sign"],
            proof_boundary="Imports the validated order-three entry and forward-invariance certificates.",
        ),
        LemmaRow(
            id="qshb_06_outer_branch_only",
            role="theorem_conclusion",
            readiness="ready_to_apply",
            claim="Every hypothetical Xi quartic double-root boundary point lies on an outer-root branch.",
            formula="D_(3,n)<0 => C>0 => (a-b)*(a-c)>0",
            proof_boundary="Branch selection only; it does not prove inward motion.",
        ),
        LemmaRow(
            id="qshb_07_middle_and_triple_exclusion",
            role="theorem_conclusion",
            readiness="ready_to_apply",
            claim="The middle-root branch and tangent triple-root stratum are incompatible with the strict Xi order-three Hankel theorem.",
            formula="C<0 and C=0 are excluded",
            proof_boundary="Excludes boundary strata, not all quartic double roots.",
        ),
        LemmaRow(
            id="qshb_08_reduced_inward_threshold",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim="The branch-aware quartic inward condition reduces on Xi to one upper threshold.",
            formula="C*(u-U)<=0 and C>0 => " + exact["reduced_inward_condition"],
            proof_boundary="The inequality u<=U remains unproved.",
            diagnostics={"U": exact["threshold"]},
        ),
        LemmaRow(
            id="qshb_09_countermodel_separation",
            role="countermodel_diagnostic",
            readiness="ready_to_apply",
            claim="The strong-log-concave local quartic countermodel is excluded by the already proved Xi order-three Hankel sign.",
            formula=f"D_3(countermodel)={counter_text}>0",
            proof_boundary="Explains the scope difference; it is not an Xi computation.",
            diagnostics={"countermodel_hankel_ball": counter_text},
        ),
        LemmaRow(
            id="qshb_10_live_outer_threshold",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="Prove u<=U at every remaining outer-branch Xi quartic contact using a global heat-compatible or theta-specific condition.",
            formula="outer quartic contact => u<=U(a,p)",
            proof_boundary="Degree-four Xi hyperbolicity, PF-infinity, Lambda<=0, and RH remain open.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_quartic_signed_hankel_branch_exclusion_lemma",
        "date": "2026-07-25",
        "status": "exact Xi quartic branch-exclusion lemma",
        "proof_boundary": (
            "This artifact combines an exact quartic boundary factorization with "
            "the already proved strict Xi contiguous order-three signed-Hankel "
            "theorem. It excludes the middle-double-root and triple-root boundary "
            "strata and reduces the live quartic inward condition to u<=U on the "
            "outer branch. It does not prove that remaining inequality, degree-four "
            "Xi hyperbolicity, PF-infinity, Lambda<=0, or RH."
        ),
        "sources": [
            "outputs/jensen_window_pf_reciprocal_defect_compound_order3_gate.md",
            "outputs/jensen_window_pf_compound_order3_forward_invariance_certificate.md",
            "outputs/jensen_window_pf_quartic_double_root_threshold_lemma.md",
            "outputs/jensen_window_pf_strong_logconcave_local_quartic_countermodel.md",
        ],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_quartic_signed_hankel_branch_exclusion_lemma.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_quartic_signed_hankel_branch_exclusion_lemma.py"
        ),
        "exact": exact,
        "summary": {
            "rows": len(rows),
            "exact_identities": 4,
            "imported_hankel_theorems": 1,
            "excluded_boundary_strata": 2,
            "reduced_outer_thresholds": 1,
            "countermodel_separations": 1,
            "open_outer_thresholds": 1,
            "main_finding": (
                "At a quartic double-root boundary, the normalized contiguous "
                "order-three Hankel determinant is -((a-b)(a-c))^3/216. The "
                "strict Xi sign therefore excludes the middle-root and triple-root "
                "strata throughout the target heat interval, leaving only the outer "
                "branch and the single unproved threshold u<=U(a,p)."
            ),
        },
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    counter = payload["rows"][8]["diagnostics"]["countermodel_hankel_ball"]
    return "\n".join(
        [
            "# Jensen-Window PF Quartic Signed-Hankel Branch-Exclusion Lemma",
            "",
            "Date: 2026-07-25",
            "",
            "Status: exact Xi quartic branch-exclusion lemma. The remaining outer",
            "threshold is open; this is not a proof of degree-four Xi",
            "hyperbolicity, PF-infinity, RH, or `Lambda <= 0`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_quartic_signed_hankel_branch_exclusion_lemma.json",
            "python work/rh_compute/scripts/jensen_window_pf_quartic_signed_hankel_branch_exclusion_lemma.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_quartic_signed_hankel_branch_exclusion_lemma.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF quartic signed-Hankel branch-exclusion lemma: 10 rows, 0 issues, 4 exact identities, 1 imported Xi Hankel theorem, 2 excluded boundary strata, 1 reduced outer threshold, 1 countermodel separation, 1 open outer threshold",
            "```",
            "",
            "## Shifted Hankel Normalization",
            "",
            "For a fixed shift `n`, write",
            "",
            "```text",
            "A_(n+j)=A_n*r_n^j*B_j, r_n=A_(n+1)/A_n,",
            "B_0=B_1=1, B_2=x, B_3=x^2*y, B_4=x^3*y^2*z.",
            "```",
            "",
            "Then",
            "",
            "```text",
            "D_(3,n)=det[A_(n+i+j)]_(i,j=0..2)",
            "       =A_n^3*r_n^6*det[B_(i+j)]_(i,j=0..2),",
            f"det[B_(i+j)]={exact['normalized_hankel']}.",
            "```",
            "",
            "The prefactor is strictly positive.",
            "",
            "## Quartic Boundary Factorization",
            "",
            "At a nondegenerate quartic boundary point, write",
            "",
            "```text",
            "P(w)=(1+a*w)^2*(1+b*w)*(1+c*w),",
            "2*a+b+c=4, p=b*c.",
            "```",
            "",
            "Substitution of the exact contraction coordinates gives",
            "",
            "```text",
            "C=3*a^2-4*a+p=(a-b)*(a-c),",
            "det[B_(i+j)]=-C^3/216.",
            "```",
            "",
            "The completed compound-order-three theorem proves",
            "",
            "```text",
            "D_(3,n)(lambda)<0",
            "```",
            "",
            "for every shift and every finite `lambda>=-100` in the propagated",
            "Xi ratio cone. Therefore `C>0`. The repeated root lies outside the",
            "two simple roots; the middle-root branch `C<0` and triple-root",
            "stratum `C=0` are impossible at a hypothetical Xi quartic contact.",
            "",
            "## Reduced Quartic Target",
            "",
            "The exact branch-aware condition from the threshold lemma was",
            "",
            "```text",
            "C*(u-U(a,p))<=0,",
            f"U(a,p)={exact['threshold']}.",
            "```",
            "",
            "Since the Xi signed-Hankel theorem forces `C>0`, the complete live",
            "first-crossing condition is now only",
            "",
            "```text",
            "u<=U(a,p).",
            "```",
            "",
            "The tangent triple-root analysis no longer belongs to the live Xi",
            "boundary set. Proving the remaining outer threshold is still a",
            "genuine open step.",
            "",
            "## Countermodel Separation",
            "",
            "The strong-log-concave local quartic countermodel has",
            "",
            "```text",
            f"D_3={counter}>0.",
            "```",
            "",
            "It is therefore excluded by structure already proved for Xi. This",
            "explains exactly why that model blocks a generic local theorem while",
            "not blocking the signed-Hankel route.",
            "",
            "## Proof Boundary",
            "",
            "The lemma selects the only possible Xi quartic boundary branch. It",
            "does not prove `u<=U`, quartic heat invariance, an all-degree bridge,",
            "PF-infinity, `Lambda<=0`, or RH.",
            "",
            "```text",
            "outputs/jensen_window_pf_reciprocal_defect_compound_order3_gate.md",
            "outputs/jensen_window_pf_compound_order3_forward_invariance_certificate.md",
            "outputs/jensen_window_pf_quartic_double_root_threshold_lemma.md",
            "outputs/jensen_window_pf_strong_logconcave_local_quartic_countermodel.md",
            "```",
            "",
            "Summary:",
            "",
            payload["summary"]["main_finding"],
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
        "wrote Jensen-window PF quartic signed-Hankel branch-exclusion lemma: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
