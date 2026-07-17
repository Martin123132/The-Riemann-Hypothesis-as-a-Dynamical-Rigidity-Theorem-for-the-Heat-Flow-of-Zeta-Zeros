#!/usr/bin/env python3
"""Reduce the positive order-ten endpoint tail to one curvature ceiling."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_m100_tail_curvature_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order10_m100_tail_curvature_reduction.md"
)
FINITE_SPLICE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_m100_finite_splice_certificate.json"
)
ORDER9_ENTRY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_m100_entry_certificate.json"
)
ORDER9_TAIL_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_m100_tail_curvature_reduction.json"
)
ORDER9_CURVATURE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_first_summand_curvature_certificate.json"
)
DEFECT_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order4_m100_entry_certificate.json"
)
HEAT_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_all_order_endpoint_heat_reduction.json"
)

FINITE_FIRST_N = 4
FINITE_LAST_N = 1242
TAIL_FIRST_N = 1243
TAIL_FIRST_K = TAIL_FIRST_N + 9
CURVATURE_CONSTANT = 5500
DEFECT_NUMERATOR = 2259
DEFECT_DENOMINATOR = 250


@dataclass(frozen=True)
class ReductionRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_sources() -> dict:
    finite_splice = load_json(FINITE_SPLICE_SOURCE)
    order9_entry = load_json(ORDER9_ENTRY_SOURCE)
    order9_tail = load_json(ORDER9_TAIL_SOURCE)
    order9_curvature = load_json(ORDER9_CURVATURE_SOURCE)
    defect = load_json(DEFECT_SOURCE)
    heat = load_json(HEAT_SOURCE)

    finite = finite_splice.get("finite", {})
    if finite.get("preserved_negative_indices") != [0, 1, 2, 3]:
        raise RuntimeError("order-ten negative endpoint prefix changed")
    if finite.get("combined_positive_range") != [FINITE_FIRST_N, FINITE_LAST_N]:
        raise RuntimeError("order-ten positive endpoint range changed")
    if finite_splice.get("summary", {}).get("direct_Q10_checks") != 2:
        raise RuntimeError("order-ten finite splice is not closed")
    if order9_entry.get("summary", {}).get("all_shift_m100_entry_theorems") != 1:
        raise RuntimeError("all-shift order-nine endpoint input is not closed")
    if "Q_(8,n)=A_(n+7)^8*exp(w(n+7))" not in str(
        order9_tail.get("exact", {}).get("canonical_factorization", "")
    ):
        raise RuntimeError("order-nine nested coordinate changed")
    if (
        order9_curvature.get("summary", {}).get(
            "global_first_summand_curvature_theorems"
        )
        != 1
    ):
        raise RuntimeError("order-nine first-summand curvature input is not closed")
    defect_anchor = defect.get("tail_arithmetic", {}).get("defect_buffer")
    if defect_anchor != "-3*log(x_k)>3*d_k>=753/(250*(2*k+1))":
        raise RuntimeError(f"scaled defect anchor changed: {defect_anchor!r}")
    if heat.get("summary", {}).get("all_fixed_order_tail_theorems") != 1:
        raise RuntimeError("all-fixed-order heat tail input is not closed")
    if heat.get("summary", {}).get("all_order_cooperative_recursions") != 1:
        raise RuntimeError("cooperative heat-flow input is not closed")

    paths = (
        FINITE_SPLICE_SOURCE,
        ORDER9_ENTRY_SOURCE,
        ORDER9_TAIL_SOURCE,
        ORDER9_CURVATURE_SOURCE,
        DEFECT_SOURCE,
        HEAT_SOURCE,
    )
    return {
        "finite_endpoint": (
            "Q_(10,n)(-100)>0 for 4<=n<=1242 and "
            "Q_(10,n)(-100)<0 for 0<=n<=3"
        ),
        "order9_endpoint": (
            "Q_(9,n)(-100)>0 for every integer n>=0"
        ),
        "previous_continuous_curvature": order9_curvature.get("theorem"),
        "defect_anchor": "d_k>=251/(250*(2*k+1)), k>=320",
        "heat_tail": heat.get("exact", {}).get("fixed_order_tail"),
        "heat_flow": heat.get("exact", {}).get("cooperative_flow"),
        "files": [
            {
                "path": path.relative_to(REPO_ROOT).as_posix(),
                "sha256": sha256(path),
            }
            for path in paths
        ],
    }


def exact_reduction() -> dict:
    a, w, s, gap = sp.symbols("A w s W", positive=True)
    q8_center = a**8 * sp.exp(w)
    q7_denominator = a**7 * sp.exp(s)
    q9 = sp.factor(q8_center**2 * (1 - sp.exp(-gap)) / q7_denominator)
    target = a**9 * sp.exp(2 * w - s) * (1 - sp.exp(-gap))
    if sp.simplify(q9 - target) != 0:
        raise RuntimeError("canonical order-nine factorization failed")

    k, m = sp.symbols("k m", integer=True, nonnegative=True)
    comparison = sp.expand(
        DEFECT_NUMERATOR * k**2
        - DEFECT_DENOMINATOR * CURVATURE_CONSTANT * (2 * k + 1)
    )
    shifted = sp.expand(comparison.subs(k, TAIL_FIRST_K + m))
    coefficients = sp.Poly(shifted, m).all_coeffs()
    if any(value <= 0 for value in coefficients):
        raise RuntimeError("order-ten tail comparison is not coefficient-positive")

    return {
        "completed_coordinates": (
            "V(t)=7*B(t)-s(t-1)+2*s(t)-s(t+1); "
            "w(t)=2*s(t)-r(t)+log(1-exp(-V(t)))"
        ),
        "seventh_gap": "W(t)=8*B(t)-w(t-1)+2*w(t)-w(t+1)",
        "order9_coordinate": (
            "z(t)=2*w(t)-s(t)+log(1-exp(-W(t)))"
        ),
        "canonical_factorization": "Q_(9,n)=A_(n+8)^9*exp(z(n+8))",
        "canonical_factorization_residual": "0",
        "curvature_identity": (
            "E_n=log(Q_(9,n)*Q_(9,n+2)/Q_(9,n+1)^2)="
            "9*log(x_k)+Z_k, Z_k=z(k-1)-2*z(k)+z(k+1), k=n+9"
        ),
        "stable_margin": "L_n=exp(-E_n)-1",
        "sign_equivalence": "Q_(10,n)>0 iff L_n>0 iff E_n<0",
        "tail_index": "n>=1243 iff k=n+9>=1252",
        "sufficient_ceiling": "Z_k<=5500/k^2 for every integer k>=1252",
        "log_buffer": (
            "-9*log(x_k)>=9*d_k>=2259/(250*(2*k+1)), k>=320"
        ),
        "rational_comparison": (
            "5500/k^2<2259/(250*(2*k+1)), k>=1252"
        ),
        "cleared_polynomial": str(comparison),
        "shifted_polynomial_k_1252_plus_m": str(shifted),
        "shifted_coefficients": [str(value) for value in coefficients],
        "finite_splice": (
            "[Q_(10,n)(-100)>0 for 4<=n<=1242] plus "
            "[Q_(10,n)(-100)>0 for n>=1243] implies "
            "Q_(10,n)(-100)>0 for every n>=4"
        ),
        "conditional_heat_handoff": (
            "[Q_(10,n)(-100)>0 for every integer n>=4] implies "
            "Q_(10,n)(lambda)>0 for every n>=4 and -100<=lambda<=0"
        ),
        "negative_prefix_boundary": (
            "Q_(10,n)(-100)<0 for n=0,1,2,3 remains unchanged"
        ),
    }


def build_artifact() -> dict:
    sources = validate_sources()
    exact = exact_reduction()
    rows = [
        ReductionRow(
            "co10m100tcr_01_order9_input",
            "theorem_input",
            "ready_to_apply",
            "The completed signed order-nine endpoint layer supplies the positive lower cone for order ten.",
            sources["order9_endpoint"],
            "Endpoint lambda=-100 only; no order-ten sign is imported.",
        ),
        ReductionRow(
            "co10m100tcr_02_finite_sign_chart",
            "interval_theorem_input",
            "ready_to_apply",
            "The rigorous order-ten endpoint certificates have four negative shifts followed by a positive block through n=1242.",
            sources["finite_endpoint"],
            "Finite endpoint theorem; it neither fills the tail nor removes the four counterexamples.",
        ),
        ReductionRow(
            "co10m100tcr_03_nested_coordinate",
            "exact_identity",
            "ready_to_apply",
            "One more stable gap gives a cancellation-free coordinate for Q9.",
            exact["seventh_gap"] + "; " + exact["order9_coordinate"] + "; " + exact["canonical_factorization"],
            "Exact signed-condensation factorization only.",
        ),
        ReductionRow(
            "co10m100tcr_04_curvature",
            "exact_reduction",
            "ready_to_apply",
            "Signed order ten is equivalent to the seventh-nested curvature lying below the ninth defect buffer.",
            exact["curvature_identity"] + "; " + exact["sign_equivalence"],
            "Exact logarithmic reduction in the positive Q8,Q9 cone.",
        ),
        ReductionRow(
            "co10m100tcr_05_defect",
            "theorem_input",
            "ready_to_apply",
            "The completed coefficient-defect theorem supplies an inverse-linear order-ten buffer.",
            exact["log_buffer"],
            "Previously proved lambda=-100 defect anchor, scaled by nine.",
        ),
        ReductionRow(
            "co10m100tcr_06_arithmetic",
            "exact_inequality",
            "ready_to_apply",
            "The proposed inverse-square curvature ceiling lies strictly inside the defect buffer on the complete tail.",
            exact["rational_comparison"],
            "Coefficient-positive rational comparison.",
            {"shifted_coefficients": exact["shifted_coefficients"]},
        ),
        ReductionRow(
            "co10m100tcr_07_conditional_tail",
            "conditional_theorem",
            "ready_to_apply",
            "The scalar curvature ceiling would prove the complete positive order-ten endpoint tail.",
            "[Z_k<=5500/k^2 for k>=1252] => [Q_(10,n)(-100)>0 for n>=1243]",
            "Conditional only on the displayed full-kernel scalar ceiling.",
        ),
        ReductionRow(
            "co10m100tcr_08_conditional_splice",
            "conditional_theorem",
            "ready_to_apply",
            "The analytic tail would splice exactly to the known positive finite block.",
            exact["finite_splice"],
            "The negative shifts n=0,1,2,3 are deliberately excluded.",
        ),
        ReductionRow(
            "co10m100tcr_09_heat_handoff",
            "conditional_forward_theorem",
            "ready_to_apply",
            "The endpoint splice would make the complete n>=4 order-ten tail forward invariant under the cooperative heat flow.",
            exact["conditional_heat_handoff"],
            "Uses the compact-uniform eventual tail and positive Q9 lower layer; it says nothing yet about shifts 0 through 3.",
        ),
        ReductionRow(
            "co10m100tcr_10_open_curvature",
            "analytic_theorem_target",
            "not_ready_to_apply",
            "Prove the seventh-nested full-kernel curvature ceiling on the endpoint tail.",
            exact["sufficient_ceiling"],
            "Requires a new continuum enclosure and a controlled first-summand-to-full-kernel transfer.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order10_m100_tail_curvature_reduction",
        "date": "2026-07-16",
        "status": (
            "exact order-ten positive endpoint-tail reduction with one open "
            "seventh-nested curvature ceiling"
        ),
        "proof_boundary": (
            "This artifact proves the order-nine canonical normalization, "
            "order-ten sign reduction, exact tail arithmetic, finite-splice "
            "coordinate, and conditional n>=4 heat handoff. It does not prove "
            "the 5500/k^2 curvature ceiling, the endpoint tail, delayed entry "
            "for shifts 0 through 3, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": sources,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_to_apply_rows": 9,
            "open_rows": 1,
            "exact_factorizations": 1,
            "exact_curvature_reductions": 1,
            "coefficient_positive_comparisons": 1,
            "conditional_tail_theorems": 1,
            "conditional_finite_splices": 1,
            "conditional_heat_handoffs": 1,
            "open_curvature_targets": 1,
            "preserved_negative_endpoint_shifts": 4,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order10_m100_tail_curvature_reduction.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order10_m100_tail_curvature_reduction.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Jensen-Window PF Compound Order-Ten Lambda=-100 Tail Curvature Reduction",
        "",
        "Date: 2026-07-16",
        "",
        "Status: exact positive-tail reduction with one open seventh-nested",
        "curvature ceiling. This does not remove the four rigorous negative",
        "endpoint shifts and is not a proof of RH or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order10_m100_tail_curvature_reduction.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order10_m100_tail_curvature_reduction.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order10_m100_tail_curvature_reduction.py",
        "```",
        "",
        "## Canonical Coordinate",
        "",
        "```text",
        exact["seventh_gap"],
        exact["order9_coordinate"],
        exact["canonical_factorization"],
        "```",
        "",
        "## Order-Ten Sign",
        "",
        "```text",
        exact["curvature_identity"],
        exact["stable_margin"],
        exact["sign_equivalence"],
        exact["tail_index"],
        "```",
        "",
        "## Tail Budget",
        "",
        "```text",
        exact["log_buffer"],
        exact["sufficient_ceiling"],
        exact["rational_comparison"],
        exact["shifted_polynomial_k_1252_plus_m"] + ">0 for m>=0,",
        "coefficients=" + str(exact["shifted_coefficients"]),
        "```",
        "",
        "Thus the full-kernel ceiling would prove `Q_(10,n)(-100)>0` for",
        "every `n>=1243`. It would meet the rigorous finite block",
        "`4<=n<=1242` exactly and give positivity for every `n>=4`.",
        "",
        "## Boundary And Handoff",
        "",
        "The endpoint signs at `n=0,1,2,3` remain negative. Conditional on",
        "the positive `n>=4` splice, the existing cooperative system and",
        "uniform eventual tail propagate `Q_(10,n)(lambda)>0` for every",
        "`n>=4` and `-100<=lambda<=0`. The four delayed-entry crossings",
        "remain a separate finite-shift theorem.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order10_m100_finite_splice_certificate.md",
        "outputs/jensen_window_pf_compound_order9_m100_entry_certificate.md",
        "outputs/jensen_window_pf_compound_order9_first_summand_curvature_certificate.md",
        "outputs/jensen_window_pf_all_order_endpoint_heat_reduction.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-ten endpoint-tail curvature reduction: "
        f"{summary['rows']} rows, "
        f"{summary['coefficient_positive_comparisons']} positive comparison, "
        f"{summary['open_curvature_targets']} open curvature target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
