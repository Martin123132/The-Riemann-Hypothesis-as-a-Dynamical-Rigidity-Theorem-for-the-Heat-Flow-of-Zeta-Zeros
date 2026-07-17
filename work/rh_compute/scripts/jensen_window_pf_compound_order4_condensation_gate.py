#!/usr/bin/env python3
"""Build the contiguous order-four condensation target and countermodel gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import itertools
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate as m100  # noqa: E402
import flint  # noqa: E402


LOCAL_REPAIR = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "acb_enclosures_negative_lambda_m100_order3_k208_k217_dps220.jsonl"
)
COLLAR_321 = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "acb_enclosures_negative_lambda_m100_order3_k321_dps220.jsonl"
)
COLLAR_322 = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "negative_lambda_m100_order4_k322_dps220.jsonl"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order4_condensation_gate.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_compound_order4_condensation_gate.md"
DEFAULT_PRECISION_BITS = 1024


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def arb_positive(value: flint.arb) -> bool:
    return bool(value > 0 and not value.contains(0))


def arb_text(value: flint.arb, digits: int = 60) -> str:
    return value.str(digits).replace("e", "E")


def lower_text(value: flint.arb, digits: int = 60) -> str:
    return m100.arb_lower_text(value, digits)


def merged_coefficients() -> tuple[dict[int, flint.arb], list[dict]]:
    values, sources = m100.merged_coefficients()
    diagnostics = list(sources)
    for precedence, source in enumerate(
        (LOCAL_REPAIR, COLLAR_321, COLLAR_322), start=len(sources)
    ):
        loaded = m100.load_source(source)
        overwritten = len(set(values).intersection(loaded))
        values.update(loaded)
        diagnostics.append(
            {
                "precedence": precedence,
                "source": source.relative_to(REPO_ROOT).as_posix(),
                "lambda_minus_100_rows": len(loaded),
                "index_range": [min(loaded), max(loaded)],
                "overwritten_rows": overwritten,
            }
        )
    if set(values) != set(range(323)):
        raise RuntimeError("order-four source merge does not cover A_0..A_322")
    return values, diagnostics


def finite_diagnostics() -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    values, sources = merged_coefficients()
    ratios = {k: values[k + 1] / values[k] for k in range(322)}
    contractions = {k: ratios[k] / ratios[k - 1] for k in range(1, 322)}
    defects = {k: 1 - contractions[k] for k in contractions}
    gaps = {
        n: defects[n + 2] ** 2
        - contractions[n + 2] ** 2 * defects[n + 1] * defects[n + 3]
        for n in range(319)
    }
    if not all(arb_positive(value) for value in gaps.values()):
        raise RuntimeError("order-three gap input failed")
    margins = {
        n: gaps[n + 1] ** 2
        / (contractions[n + 3] ** 3 * gaps[n] * gaps[n + 2])
        - 1
        for n in range(317)
    }
    if not all(arb_positive(value) for value in margins.values()):
        failed = [n for n, value in margins.items() if not arb_positive(value)]
        raise RuntimeError(f"order-four prefix margins failed: {failed[:10]}")
    minimum_n = min(margins, key=lambda n: margins[n].lower())
    maximum_n = max(margins, key=lambda n: margins[n].upper())
    penalties = {
        n: (gaps[n] * gaps[n + 2] / gaps[n + 1] ** 2).log() for n in range(317)
    }
    scaled_caps = {
        n: flint.arb(2) / (n + 3) ** 2 - penalties[n] for n in range(317)
    }
    if not all(arb_positive(value) for value in scaled_caps.values()):
        failed = [n for n, value in scaled_caps.items() if not arb_positive(value)]
        raise RuntimeError(f"order-four scaled penalty caps failed: {failed[:10]}")
    minimum_scaled_cap_n = min(scaled_caps, key=lambda n: scaled_caps[n].lower())
    maximum_scaled_penalty_n = max(
        penalties, key=lambda n: (penalties[n] * (n + 3) ** 2).upper()
    )
    selected_n = (0, 1, 2, 10, 50, 100, 200, 300, 315, 316)
    return {
        "parameters": {"lambda": "-100", "precision_bits": DEFAULT_PRECISION_BITS},
        "merged_sources": sources,
        "coefficient_range": [0, 322],
        "positive_coefficients": len(values),
        "order_three_gap_range": [0, 318],
        "positive_order_three_gaps": len(gaps),
        "order_four_margin_range": [0, 316],
        "positive_order_four_margins": len(margins),
        "scaled_penalty_cap": "P_n<2/(n+3)^2",
        "positive_scaled_penalty_caps": len(scaled_caps),
        "minimum_scaled_penalty_cap_at_n": minimum_scaled_cap_n,
        "minimum_scaled_penalty_cap_lower": lower_text(scaled_caps[minimum_scaled_cap_n]),
        "maximum_scaled_penalty_at_n": maximum_scaled_penalty_n,
        "maximum_scaled_penalty_ball": arb_text(
            penalties[maximum_scaled_penalty_n] * (maximum_scaled_penalty_n + 3) ** 2
        ),
        "minimum_margin_at_n": minimum_n,
        "minimum_margin_ball": arb_text(margins[minimum_n]),
        "minimum_margin_lower": lower_text(margins[minimum_n]),
        "maximum_margin_at_n": maximum_n,
        "maximum_margin_ball": arb_text(margins[maximum_n]),
        "selected_margins": [
            {"n": n, "margin_ball": arb_text(margins[n]), "margin_lower": lower_text(margins[n])}
            for n in selected_n
        ],
    }


def exact_countermodel() -> dict:
    scaled_defects = (
        sp.Rational(3, 10),
        sp.Rational(2, 5),
        sp.Rational(41, 100),
        sp.Rational(49, 100),
        sp.Rational(11, 20),
    )
    defects = tuple(
        sp.factor(2 * scaled_defects[k - 1] / (2 * k + 1)) for k in range(1, 6)
    )
    contractions = tuple(sp.factor(1 - value) for value in defects)
    q = tuple(sp.sqrt(1 / value) for value in defects)
    increments = tuple(sp.factor(q[j + 1] - q[j]) for j in range(4))
    compound = tuple(
        sp.factor(q[n] * q[n + 2] - q[n + 1] ** 2 + 1) for n in range(3)
    )
    gaps = tuple(
        sp.factor(
            defects[n + 1] ** 2
            - contractions[n + 1] ** 2 * defects[n] * defects[n + 2]
        )
        for n in range(3)
    )
    frontier = sp.factor(gaps[1] ** 2 - contractions[2] ** 3 * gaps[0] * gaps[2])
    ratio = sp.factor(gaps[1] ** 2 / (contractions[2] ** 3 * gaps[0] * gaps[2]))
    cubic_frontiers = tuple(
        sp.factor(
            contractions[k] ** 2 * contractions[k + 1] ** 2
            - 6 * contractions[k] * contractions[k + 1]
            + 4 * contractions[k]
            + 4 * contractions[k + 1]
            - 3
        )
        for k in range(4)
    )

    ratios = [sp.Integer(1)]
    for contraction in contractions:
        ratios.append(sp.factor(ratios[-1] * contraction))
    coefficients = [sp.Integer(1)]
    for ratio_item in ratios:
        coefficients.append(sp.factor(coefficients[-1] * ratio_item))

    h2 = [
        sp.factor(sp.Matrix([[coefficients[n + i + j] for j in range(2)] for i in range(2)]).det())
        for n in range(5)
    ]
    h3 = [
        sp.factor(sp.Matrix([[coefficients[n + i + j] for j in range(3)] for i in range(3)]).det())
        for n in range(3)
    ]
    h4 = sp.factor(sp.Matrix([[coefficients[i + j] for j in range(4)] for i in range(4)]).det())
    arbitrary_h3 = []
    for n in range(3):
        max_column = 4 - n
        for columns in itertools.combinations(range(max_column + 1), 3):
            determinant = sp.factor(
                sp.Matrix(
                    [[coefficients[n + i + column] for column in columns] for i in range(3)]
                ).det()
            )
            arbitrary_h3.append({"n": n, "columns": list(columns), "determinant": str(determinant)})

    lower_walls = tuple(sp.Rational(2 * k - 1, 2 * k + 1) for k in range(1, 6))
    checks = {
        "pointwise_ratio_cone": all(contractions[k - 1] > lower_walls[k - 1] for k in range(1, 6)),
        "increasing_contractions": all(contractions[k + 1] > contractions[k] for k in range(4)),
        "decreasing_defects": all(defects[k + 1] < defects[k] for k in range(4)),
        "increasing_scaled_defects": all(
            scaled_defects[k + 1] > scaled_defects[k] for k in range(4)
        ),
        "strict_cubic_frontiers": all(value < 0 for value in cubic_frontiers),
        "positive_compound_gaps": all(value > 0 for value in gaps),
        "negative_h2": all(value < 0 for value in h2),
        "negative_contiguous_h3": all(value < 0 for value in h3),
        "negative_arbitrary_h3": all(sp.Rational(row["determinant"]) < 0 for row in arbitrary_h3),
        "wrong_negative_h4": h4 < 0,
    }
    checks = {key: bool(value) for key, value in checks.items()}
    if not all(checks.values()):
        raise RuntimeError(f"order-four countermodel failed: {checks}")
    return {
        "reciprocal_defects": [str(value) for value in q],
        "scaled_defects": [str(value) for value in scaled_defects],
        "defects": [str(value) for value in defects],
        "contractions": [str(value) for value in contractions],
        "lower_walls": [str(value) for value in lower_walls],
        "increments": [str(value) for value in increments],
        "compound_margins": [str(value) for value in compound],
        "compound_gaps": [str(value) for value in gaps],
        "cubic_frontiers": [str(value) for value in cubic_frontiers],
        "order_four_gap_frontier": str(frontier),
        "order_four_gap_ratio": str(ratio),
        "coefficient_ratios": [str(value) for value in ratios],
        "coefficients": [str(value) for value in coefficients],
        "h2_determinants": [str(value) for value in h2],
        "h3_determinants": [str(value) for value in h3],
        "arbitrary_h3_determinants": arbitrary_h3,
        "h4_determinant": str(h4),
        "checks": checks,
    }


def exact_diagnostics() -> dict:
    h2, h3a, h3b, h3c, h4 = sp.symbols("h2 h3a h3b h3c h4")
    condensation = sp.Eq(h4 * h2, h3a * h3c - h3b**2)
    return {
        "condensation": (
            "H_(4,n)*H_(2,n+2)=H_(3,n)*H_(3,n+2)-H_(3,n+1)^2"
        ),
        "sign_equivalence": (
            "H_(2,n+2)<0 and H_(3,j)<0 => H_(4,n)>0 iff "
            "T_(n+1)^2>T_n*T_(n+2), T_j=-H_(3,j)>0"
        ),
        "gap_definition": (
            "G_n=d_(n+2)^2-x_(n+2)^2*d_(n+1)*d_(n+3)>0"
        ),
        "h3_factor": "T_n=A_n^3*r_n^6*x_(n+1)^3*G_n",
        "ratio_identity": (
            "T_(n+1)^2/(T_n*T_(n+2))="
            "G_(n+1)^2/(x_(n+3)^3*G_n*G_(n+2))"
        ),
        "order_four_target": (
            "G_(n+1)^2>x_(n+3)^3*G_n*G_(n+2) for every n>=0"
        ),
        "log_penalty": (
            "P_n=log(G_n*G_(n+2)/G_(n+1)^2), k=n+3"
        ),
        "anchor_buffer": (
            "lambda=-100, k>=319: -3*log(x_k)>3*d_k>="
            "753/(250*(2*k+1))"
        ),
        "sufficient_tail": "P_n<=4/(n+3)^2 for every n>=317",
        "rational_comparison": (
            "753/(250*(2*k+1))-4/k^2>0 for k>=320; "
            "cleared numerator at k=320 is 76466200"
        ),
        "symbolic_condensation": str(condensation),
    }


def build_artifact() -> dict:
    exact = exact_diagnostics()
    finite = finite_diagnostics()
    countermodel = exact_countermodel()
    rows = [
        GateRow(
            id="co4cg_01_condensation",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Desnanot-Jacobi condensation reduces the contiguous order-four determinant to three neighboring lower compounds.",
            formula=exact["condensation"],
            proof_boundary="Contiguous Hankel determinants only.",
        ),
        GateRow(
            id="co4cg_02_log_concavity_equivalence",
            role="exact_equivalence",
            readiness="ready_to_apply",
            claim="Inside the completed order-two and order-three cones, the order-four sign is exactly log-concavity of the order-three determinant magnitudes.",
            formula=exact["sign_equivalence"],
            proof_boundary="Exact order-four coordinate, not its proof for Xi.",
        ),
        GateRow(
            id="co4cg_03_gap_coordinate",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The order-four condition is buffered log-concavity of the positive order-three defect gaps.",
            formula=exact["order_four_target"],
            proof_boundary="Exact scalar coordinate for contiguous order four.",
            diagnostics={"gap": exact["gap_definition"], "ratio": exact["ratio_identity"]},
        ),
        GateRow(
            id="co4cg_04_m100_prefix",
            role="interval_certificate",
            readiness="finite_validated",
            claim="Every available repaired lambda=-100 prefix satisfies the strict order-four gap coordinate.",
            formula=(
                "order-four margin>0 and P_n<2/(n+3)^2, 0<=n<=316"
            ),
            proof_boundary="Finite prefix only; no analytic tail.",
            diagnostics=finite,
        ),
        GateRow(
            id="co4cg_05_lower_cone_countermodel",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="A rational positive sequence lies in the strict ratio and cubic cones and has all available order-three signs, yet fails order four.",
            formula="s=(3/10,2/5,41/100,49/100,11/20)",
            proof_boundary="Abstract finite countermodel, not the Xi sequence.",
            diagnostics=countermodel,
        ),
        GateRow(
            id="co4cg_06_wrong_h4",
            role="forbidden_promotion",
            readiness="guard_validated",
            claim="Completed lower compound layers cannot be promoted to order four without the new log-concavity inequality.",
            formula=f"H_(4,0)={countermodel['h4_determinant']}<0",
            proof_boundary="Rejects a structural implication only.",
        ),
        GateRow(
            id="co4cg_07_tail_target",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="Prove the buffered G-gap log-concavity on the lambda=-100 tail and splice it to the interval prefix.",
            formula=(
                "co4fcb_05 continuous ceiling => P_n<=4/(n+3)^2, n>=317"
            ),
            proof_boundary="Open all-index entry theorem.",
        ),
        GateRow(
            id="co4cg_08_flow_and_columns",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="After entry, derive a cooperative order-four flow and an arbitrary-column transfer theorem.",
            formula="order-four entry + forward invariance + all-column transfer",
            proof_boundary="Open; not all-order sign regularity, PF-infinity, RH, or Lambda<=0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order4_condensation_gate",
        "date": "2026-07-12",
        "status": (
            "exact contiguous order-four condensation coordinate with 317-row "
            "lambda=-100 prefix and strict lower-order countermodel"
        ),
        "proof_boundary": (
            "This artifact derives the exact contiguous order-four coordinate, "
            "certifies a finite lambda=-100 prefix, and supplies a strict rational "
            "countermodel to promotion from the completed lower cones. It does not "
            "prove the all-index order-four tail, forward invariance, arbitrary "
            "order-four columns, the all-order bridge, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.md",
            "outputs/jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.md",
            "outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md",
            "outputs/signed_hankel_jensen_bridge_target.md",
            "outputs/formal_core.md",
            LOCAL_REPAIR.relative_to(REPO_ROOT).as_posix(),
            COLLAR_321.relative_to(REPO_ROOT).as_posix(),
            COLLAR_322.relative_to(REPO_ROOT).as_posix(),
        ],
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order4_condensation_gate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order4_condensation_gate.py",
        "exact": exact,
        "finite": finite,
        "countermodel": countermodel,
        "summary": {
            "gate_rows": len(rows),
            "exact_identity_rows": 2,
            "exact_equivalence_rows": 1,
            "prefix_order_four_margins": finite["positive_order_four_margins"],
            "prefix_scaled_penalty_caps": finite["positive_scaled_penalty_caps"],
            "strict_lower_order_countermodels": 1,
            "forbidden_promotions": 1,
            "open_handoffs": 2,
            "ready_to_apply_rows": sum(row.readiness == "ready_to_apply" for row in rows),
        },
        "rows": [asdict(row) for row in rows],
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    finite = artifact["finite"]
    countermodel = artifact["countermodel"]
    lines = [
        "# Jensen-Window PF Compound Order-Four Condensation Gate",
        "",
        "Date: 2026-07-12",
        "",
        "Status: exact contiguous order-four condensation coordinate with 317-row",
        "lambda=-100 prefix and strict lower-order countermodel. This is not a",
        "proof of PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order4_condensation_gate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_condensation_gate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_condensation_gate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_condensation_gate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        "validated Jensen-window PF compound order-four condensation gate: 8 rows, 0 issues, 2 exact identities, 1 exact sign equivalence, 317 positive lambda=-100 prefix margins, 1 strict lower-order countermodel, 1 forbidden promotion, 2 open handoffs",
        "```",
        "",
        "## Exact Condensation",
        "",
        "Write `H_(m,n)=det[A_(n+i+j)]_(i,j=0..m-1)`. Desnanot-Jacobi gives",
        "",
        "```text",
        exact["condensation"],
        "```",
        "",
        "The completed lower layers have `H_(2,n)<0` and `H_(3,n)<0`. Put",
        "`T_n=-H_(3,n)>0`. Then",
        "",
        "```text",
        "H_(4,n)>0 iff T_(n+1)^2>T_n*T_(n+2).",
        "```",
        "",
        "In the order-three defect coordinate",
        "",
        "```text",
        exact["gap_definition"],
        exact["h3_factor"],
        exact["ratio_identity"],
        "```",
        "",
        "so the exact new target is",
        "",
        "```text",
        exact["order_four_target"],
        "```",
        "",
        "## Repaired Prefix",
        "",
        "At `lambda=-100`, 1024-bit Arb arithmetic certifies all `317` available",
        "margins on `0<=n<=316`. The smallest lower bound is",
        "",
        "```text",
        f"{finite['minimum_margin_lower']} at n={finite['minimum_margin_at_n']}.",
        "```",
        "",
        "The largest checked margin is approximately `0.0263397201892780` at",
        "`n=0`; the margins decrease to approximately `0.00471789538689370` at",
        "the current collar. This is rigorous finite evidence, not a tail theorem.",
        "",
        "The same balls prove the stronger finite curvature cap",
        "",
        "```text",
        "P_n=log(G_n*G_(n+2)/G_(n+1)^2)<2/(n+3)^2, 0<=n<=316.",
        "```",
        "",
        "For the analytic tail put `k=n+3>=320`. The completed scaled-defect",
        "theorem and `s_319>251/500` give",
        "",
        "```text",
        "-3*log(x_k)>3*d_k>=753/(250*(2*k+1)).",
        "```",
        "",
        "Therefore the single sufficient tail estimate",
        "",
        "```text",
        "P_n<=4/(n+3)^2, n>=317,",
        "```",
        "",
        "would prove the order-four inequality. The rational comparison is exact:",
        "after clearing denominators its numerator is `76466200>0` at `k=320`",
        "and increases thereafter.",
        "",
        "## Strict Countermodel",
        "",
        "Take scaled defects",
        "",
        "```text",
        f"s_1,...,s_5={countermodel['scaled_defects']}",
        f"d_1,...,d_5={countermodel['defects']}",
        f"x_1,...,x_5={countermodel['contractions']}",
        f"cubic frontiers={countermodel['cubic_frontiers']}",
        f"G_0,G_1,G_2={countermodel['compound_gaps']}",
        "```",
        "",
        "The contractions satisfy every strict pointwise ratio wall and increase;",
        "the defects decrease; the scaled defects increase; all four cubic",
        "frontiers are strictly negative; and every available order-two,",
        "contiguous order-three, and arbitrary-column order-three determinant has",
        "the required sign. Nevertheless",
        "",
        "```text",
        f"G_1^2-x_3^3*G_0*G_2={countermodel['order_four_gap_frontier']}<0,",
        f"H_(4,0)={countermodel['h4_determinant']}<0.",
        "```",
        "",
        "Thus complete lower compound layers do not imply order four. The next",
        "analytic task is the all-index lambda=-100 G-gap inequality, followed by",
        "a cooperative flow law and arbitrary-column transfer.",
        "",
        "```text",
        "outputs/jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.md",
        "outputs/jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.md",
        "outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md",
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
        "validated Jensen-window PF compound order-four condensation gate: "
        "8 rows, 0 issues, 2 exact identities, 1 exact sign equivalence, "
        "317 positive lambda=-100 prefix margins, 1 strict lower-order "
        "countermodel, 1 forbidden promotion, 2 open handoffs"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
