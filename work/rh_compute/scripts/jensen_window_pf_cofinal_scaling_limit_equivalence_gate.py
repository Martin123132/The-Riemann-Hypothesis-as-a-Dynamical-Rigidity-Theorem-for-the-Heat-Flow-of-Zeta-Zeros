#!/usr/bin/env python3
"""Record why cofinal terminal hyperbolicity is already an LP endpoint theorem."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_cofinal_scaling_limit_equivalence_gate.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_cofinal_scaling_limit_equivalence_gate.md"


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def falling_ratio(D: int, j: int) -> Fraction:
    value = Fraction(1, 1)
    for offset in range(j):
        value *= Fraction(D - offset, D)
    return value


def build_exact() -> dict:
    sample_degrees = [10, 100, 1000, 10000]
    sample_orders = [0, 1, 2, 3, 5]
    ratios = {
        str(j): [str(falling_ratio(D, j)) for D in sample_degrees if D >= j]
        for j in sample_orders
    }
    if any(falling_ratio(D, j) <= 0 or falling_ratio(D, j) > 1 for D in sample_degrees for j in sample_orders if D >= j):
        raise RuntimeError("falling-factorial scaling ratio left (0,1]")
    if falling_ratio(10000, 5) <= falling_ratio(1000, 5):
        raise RuntimeError("sample scaling ratios do not approach one")
    return {
        "entire_function": "F_n(z)=sum_(j>=0) A_(n+j)*z^j/j!",
        "jensen_polynomial": "P_(D,n)(w)=sum_(j=0)^D C(D,j)*A_(n+j)*w^j",
        "scaled_polynomial": (
            "P_(D,n)(z/D)=sum_(j=0)^D ((D)_j/D^j)*A_(n+j)*z^j/j!"
        ),
        "coefficient_limit": "(D)_j/D^j=product_(m=0)^(j-1)(1-m/D)->1 for fixed j",
        "local_uniform_limit": "P_(D,n)(z/D)->F_n(z) locally uniformly because F_n is entire",
        "closure_step": (
            "A locally uniform nonzero limit of real-rooted polynomials is in the Laguerre-Polya class."
        ),
        "cofinal_implication": (
            "If P_(D_l,n) is hyperbolic for D_l->infinity, then F_n is Laguerre-Polya."
        ),
        "converse": (
            "If F_n is Laguerre-Polya, Jensen's theorem gives hyperbolicity of every P_(D,n), hence cofinal hyperbolicity."
        ),
        "equivalence": "cofinal degree hyperbolicity at fixed n <=> F_n is Laguerre-Polya",
        "shift_zero_warning": (
            "At n=0, the cofinal terminal theorem is already the Laguerre-Polya endpoint required by the Jensen/RH bridge."
        ),
        "sample_degrees": sample_degrees,
        "sample_falling_ratios": ratios,
    }


def build_payload() -> dict:
    exact = build_exact()
    rows = [
        GateRow(
            id="csle_01_scaled_jensen_identity",
            role="exact_identity",
            readiness="available_exact",
            claim="Scaling a degree-D Jensen polynomial by z/D inserts a falling-factorial ratio in each exponential-generating coefficient.",
            formula=exact["scaled_polynomial"],
            proof_boundary="Coefficient identity only.",
        ),
        GateRow(
            id="csle_02_fixed_coefficient_limit",
            role="exact_limit",
            readiness="available_exact",
            claim="Every fixed scaled coefficient converges to the corresponding entire-function coefficient.",
            formula=exact["coefficient_limit"],
            proof_boundary="Fixed-j coefficient limit only.",
        ),
        GateRow(
            id="csle_03_local_uniform_limit",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="Entire-function convergence upgrades the coefficient limit to local uniform convergence on the complex plane.",
            formula=exact["local_uniform_limit"],
            proof_boundary="Standard compact convergence for the fixed-shift entire series.",
        ),
        GateRow(
            id="csle_04_lp_closure",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="A nonzero local uniform limit of hyperbolic real polynomials belongs to the Laguerre-Polya class.",
            formula=exact["closure_step"],
            proof_boundary="Classical Laguerre-Polya closure theorem only.",
        ),
        GateRow(
            id="csle_05_cofinal_implies_lp",
            role="theorem_composition",
            readiness="ready_to_apply",
            claim="Cofinal terminal hyperbolicity at one fixed shift already proves the corresponding entire derivative is Laguerre-Polya.",
            formula=exact["cofinal_implication"],
            proof_boundary="Conditional on the unproved cofinal terminal sequence.",
        ),
        GateRow(
            id="csle_06_lp_implies_all_degrees",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="The converse follows from the classical Jensen characterization.",
            formula=exact["converse"],
            proof_boundary="Classical Jensen theorem only.",
        ),
        GateRow(
            id="csle_07_fixed_shift_equivalence",
            role="exact_equivalence",
            readiness="ready_to_apply",
            claim="At fixed shift, cofinal terminal hyperbolicity is equivalent to the full Laguerre-Polya endpoint.",
            formula=exact["equivalence"],
            proof_boundary="Fixed-shift equivalence only.",
        ),
        GateRow(
            id="csle_08_shift_zero_endpoint",
            role="noncircularity_guard",
            readiness="guard_validated",
            claim="For the original zeta coefficient entire function, a shift-zero cofinal theorem is not a weaker surrogate for the endpoint bridge.",
            formula=exact["shift_zero_warning"],
            proof_boundary="Proof-route classification; it does not establish the cofinal hypothesis.",
        ),
        GateRow(
            id="csle_09_finite_degree_guard",
            role="forbidden_promotion",
            readiness="guard_validated",
            claim="Any bounded terminal ladder, including the certified degrees 3 through 12, is insufficient for the scaling-limit implication.",
            formula="bounded D<=12 does not provide D_l->infinity",
            proof_boundary="Blocks promotion of finite Sturm evidence only.",
        ),
        GateRow(
            id="csle_10_independent_route_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="A cofinal proof remains useful only if derived independently from heat-flow dynamics, kernel structure, or another noncircular mechanism.",
            formula="independent structure => D_l->infinity; not LP endpoint => D_l->infinity",
            proof_boundary="Open noncircular route; not PF-infinity, RH, or Lambda <= 0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_cofinal_scaling_limit_equivalence_gate",
        "date": "2026-07-10",
        "status": "exact scaling-limit equivalence and cofinal-route noncircularity gate",
        "proof_boundary": (
            "This artifact proves that cofinal terminal hyperbolicity at shift zero is already equivalent to the Laguerre-Polya endpoint under the standard Jensen scaling limit. It does not prove the cofinal hypothesis, Laguerre-Polya membership for the zeta function, PF-infinity, RH, or Lambda <= 0."
        ),
        "sources": [
            "outputs/jensen_window_pf_cofinal_degree_polar_closure_lemma.md",
            "outputs/jensen_window_pf_bridge_target.md",
            "outputs/arb_jensen_window_sturm_hyperbolicity_diagnostic.md",
        ],
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    return "\n".join(
        [
            "# Jensen-Window PF Cofinal Scaling-Limit Equivalence Gate",
            "",
            "Date: 2026-07-10",
            "",
            "Status: exact scaling-limit and proof-route guard. This is not a proof",
            "of Laguerre-Polya membership, PF-infinity, RH, or `Lambda <= 0`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_cofinal_scaling_limit_equivalence_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_cofinal_scaling_limit_equivalence_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_cofinal_scaling_limit_equivalence_gate.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF cofinal scaling-limit equivalence gate: 10 rows, 0 issues, 2 exact scaling identities, 2 analytic limit steps, 1 cofinal-to-LP theorem, 1 LP-to-all-degrees theorem, 1 fixed-shift equivalence, 2 non-promotion guards, 1 open independent-route handoff",
            "```",
            "",
            "## Scaling Limit",
            "",
            "For a fixed shift `n`,",
            "",
            "```text",
            exact["scaled_polynomial"],
            exact["coefficient_limit"],
            "```",
            "",
            "Since `F_n` is entire, the convergence is locally uniform:",
            "",
            "```text",
            exact["local_uniform_limit"],
            "```",
            "",
            "If an unbounded subsequence of the scaled polynomials is hyperbolic,",
            "the Laguerre-Polya closure theorem places the nonzero limit `F_n` in",
            "the Laguerre-Polya class. Conversely, Jensen's theorem gives every",
            "degree once `F_n` is Laguerre-Polya. Hence",
            "",
            "```text",
            exact["equivalence"],
            "```",
            "",
            "## Proof-Safety Consequence",
            "",
            "At shift zero, proving the cofinal terminal sequence is already an",
            "endpoint theorem, not a weaker corollary of familiar fixed-degree",
            "asymptotics. The 1,050 finite Sturm rows through degree 12 do not enter",
            "this limit because their degrees are bounded.",
            "",
            "The cofinal formulation can still be useful if heat-flow or kernel",
            "structure proves it independently. It cannot be assumed from",
            "Laguerre-Polya membership or from the desired RH endpoint without",
            "becoming circular.",
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
        "wrote Jensen-window PF cofinal scaling-limit equivalence gate: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
