#!/usr/bin/env python3
"""Certify all-shift contiguous order-three entry at lambda=-100."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys


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
ADAPTIVE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate.json"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.md"
)
DEFAULT_PRECISION_BITS = 1024
ANCHOR = Fraction(251, 500)
TAIL_M_MIN = 641


@dataclass(frozen=True)
class CertificateRow:
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


def upper_text(value: flint.arb, digits: int = 60) -> str:
    return value.upper().str(digits, radius=False).replace("e", "E")


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def merged_order3_coefficients() -> tuple[dict[int, flint.arb], list[dict]]:
    values, sources = m100.merged_coefficients()
    diagnostics = list(sources)
    for precedence, source in enumerate((LOCAL_REPAIR, COLLAR_321), start=len(sources)):
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
    if set(values) != set(range(322)):
        raise RuntimeError("order-three source merge does not cover A_0..A_321")
    return values, diagnostics


def exact_tail_data() -> dict:
    m = TAIL_M_MIN
    bracket = Fraction(1) + Fraction(2, m) + Fraction(2, m**2) + Fraction(1, m**3)
    lower = Fraction(1) - Fraction(250, 251) * bracket
    if lower != Fraction(57_613_471, 66_107_054_971) or lower <= 0:
        raise RuntimeError("exact tail endpoint comparison failed")
    return {
        "anchor": f"{ANCHOR.numerator}/{ANCHOR.denominator}",
        "center_index": "k=n+2",
        "tail_center_range": "k>=320",
        "tail_shift_range": "n>=318",
        "m_definition": "m=2*k+1>=641",
        "sqrt_bound": (
            "r_m:=1-sqrt(1-2/m)<1/m+1/m^2; after squaring this is "
            "m^2-2*m-1>0"
        ),
        "increment_definitions": "u_k=q_k-q_(k-1), v_k=q_(k+1)-q_k",
        "increment_bounds": (
            "q_k*u_k<(250/251)*(1+1/m), "
            "u_k*v_k<(250/251)*(1/m)*(1+1/m)^2"
        ),
        "compound_lower_bound": (
            "C_(k-2)>1-(250/251)*(1+2/m+2/m^2+1/m^3)"
        ),
        "endpoint_bracket": str(bracket),
        "uniform_tail_lower": str(lower),
        "uniform_tail_lower_decimal": f"{float(lower):.18g}",
    }


def finite_diagnostics() -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    values, sources = merged_order3_coefficients()
    if not all(arb_positive(value) for value in values.values()):
        raise RuntimeError("order-three coefficient positivity failed")

    contractions = {
        k: values[k + 1] * values[k - 1] / values[k] ** 2 for k in range(1, 321)
    }
    defects = {k: 1 - contractions[k] for k in contractions}
    if not all(arb_positive(value) for value in defects.values()):
        raise RuntimeError("order-three defect positivity failed")
    q = {k: flint.arb(1) / defects[k].sqrt() for k in defects}
    increments = {k: q[k + 1] - q[k] for k in range(1, 320)}
    increment_decreases = {
        k: increments[k] - increments[k + 1] for k in range(1, 319)
    }
    if not all(arb_positive(value) and bool(value < 1) for value in increments.values()):
        raise RuntimeError("prefix reciprocal-defect increment cone failed")
    if not all(arb_positive(value) for value in increment_decreases.values()):
        raise RuntimeError("prefix reciprocal-defect increment decrease failed")

    compound = {
        n: q[n + 1] * q[n + 3] - q[n + 2] ** 2 + 1 for n in range(318)
    }
    defect_gaps = {
        n: defects[n + 2] ** 2
        - contractions[n + 2] ** 2 * defects[n + 1] * defects[n + 3]
        for n in compound
    }
    compound_increases = {n: compound[n + 1] - compound[n] for n in range(317)}
    if not all(arb_positive(value) and bool(value < 1) for value in compound.values()):
        raise RuntimeError("prefix compound margin failed")
    if not all(arb_positive(value) for value in defect_gaps.values()):
        raise RuntimeError("prefix defect compound gap failed")
    if not all(arb_positive(value) for value in compound_increases.values()):
        raise RuntimeError("prefix compound-margin growth failed")

    s319 = flint.arb(639) * defects[319] / 2
    anchor_margin = s319 - flint.arb(ANCHOR.numerator) / ANCHOR.denominator
    if not arb_positive(anchor_margin):
        raise RuntimeError("scaled-defect tail anchor failed")

    minimum_compound_n = min(compound, key=lambda n: compound[n].lower())
    minimum_defect_gap_n = min(defect_gaps, key=lambda n: defect_gaps[n].lower())
    minimum_increment_decrease_k = min(
        increment_decreases, key=lambda k: increment_decreases[k].lower()
    )
    minimum_compound_increase_n = min(
        compound_increases, key=lambda n: compound_increases[n].lower()
    )
    selected_n = (0, 1, 2, 10, 50, 100, 200, 300, 316, 317)
    return {
        "parameters": {"lambda": "-100", "precision_bits": DEFAULT_PRECISION_BITS},
        "merged_sources": sources,
        "coefficient_range": [0, 321],
        "positive_coefficients": len(values),
        "contraction_range": [1, 320],
        "positive_defects": len(defects),
        "increment_range": [1, 319],
        "positive_subunit_increments": len(increments),
        "decreasing_increment_range": [1, 318],
        "positive_increment_decreases": len(increment_decreases),
        "minimum_increment_decrease_at_k": minimum_increment_decrease_k,
        "minimum_increment_decrease_lower": lower_text(
            increment_decreases[minimum_increment_decrease_k]
        ),
        "compound_shift_range": [0, 317],
        "positive_compound_margins": len(compound),
        "minimum_compound_at_n": minimum_compound_n,
        "minimum_compound_ball": arb_text(compound[minimum_compound_n]),
        "minimum_compound_lower": lower_text(compound[minimum_compound_n]),
        "positive_defect_gaps": len(defect_gaps),
        "minimum_defect_gap_at_n": minimum_defect_gap_n,
        "minimum_defect_gap_lower": lower_text(defect_gaps[minimum_defect_gap_n]),
        "increasing_compound_range": [0, 316],
        "positive_compound_increases": len(compound_increases),
        "minimum_compound_increase_at_n": minimum_compound_increase_n,
        "minimum_compound_increase_lower": lower_text(
            compound_increases[minimum_compound_increase_n]
        ),
        "scaled_defect_anchor": {
            "k": 319,
            "s_ball": arb_text(s319),
            "rational_floor": f"{ANCHOR.numerator}/{ANCHOR.denominator}",
            "margin_lower": lower_text(anchor_margin),
        },
        "selected_compound_rows": [
            {
                "n": n,
                "C_ball": arb_text(compound[n]),
                "C_lower": lower_text(compound[n]),
                "defect_gap_lower": lower_text(defect_gaps[n]),
            }
            for n in selected_n
        ],
        "q_320_upper": upper_text(q[320]),
    }


def build_artifact() -> dict:
    finite = finite_diagnostics()
    tail = exact_tail_data()
    adaptive = load_json(ADAPTIVE_SOURCE)
    summary = adaptive.get("summary", {})
    if summary.get("open_requirements") != 0 or summary.get("defect_conclusion_rows") != 4:
        raise RuntimeError("adaptive-defect source is not closed")

    rows = [
        CertificateRow(
            id="m100co3_01_repaired_source_merge",
            role="interval_input",
            readiness="finite_validated",
            claim="Merge the full lambda=-100 source with the local order-three precision repair and A_321 collar.",
            formula="A_k(-100), 0<=k<=321",
            proof_boundary="Finite interval input only.",
            diagnostics={"merged_sources": finite["merged_sources"]},
        ),
        CertificateRow(
            id="m100co3_02_compound_coordinate",
            role="exact_theorem_input",
            readiness="ready_to_apply",
            claim="The contiguous order-three signed-Hankel sign is exactly the reciprocal-defect margin sign.",
            formula="D_(3,n)<0 iff C_n=q_(n+1)*q_(n+3)-q_(n+2)^2+1>0",
            proof_boundary="Contiguous 3x3 minors only.",
        ),
        CertificateRow(
            id="m100co3_03_prefix_compound",
            role="interval_certificate",
            readiness="interval_validated",
            claim="Every repaired-prefix reciprocal-defect compound margin is strictly positive.",
            formula="C_n(-100)>0, 0<=n<=317",
            proof_boundary="Finite prefix only.",
            diagnostics={
                "positive_compound_margins": finite["positive_compound_margins"],
                "minimum_compound_at_n": finite["minimum_compound_at_n"],
                "minimum_compound_lower": finite["minimum_compound_lower"],
            },
        ),
        CertificateRow(
            id="m100co3_04_prefix_defect_gap",
            role="interval_crosscheck",
            readiness="interval_validated",
            claim="The equivalent defect compound gap is independently positive on the same prefix.",
            formula="d_(n+2)^2-x_(n+2)^2*d_(n+1)*d_(n+3)>0, 0<=n<=317",
            proof_boundary="Equivalent interval cross-check, not an independent all-order theorem.",
            diagnostics={
                "positive_defect_gaps": finite["positive_defect_gaps"],
                "minimum_defect_gap_at_n": finite["minimum_defect_gap_at_n"],
                "minimum_defect_gap_lower": finite["minimum_defect_gap_lower"],
            },
        ),
        CertificateRow(
            id="m100co3_05_prefix_shape",
            role="interval_diagnostic",
            readiness="diagnostic_validated",
            claim="On the repaired prefix the q increments decrease while the compound margins increase.",
            formula="Delta*q_k>Delta*q_(k+1)>0 and C_(n+1)>C_n",
            proof_boundary="Finite shape evidence; neither monotonicity is promoted beyond the prefix.",
            diagnostics={
                "positive_subunit_increments": finite["positive_subunit_increments"],
                "positive_increment_decreases": finite["positive_increment_decreases"],
                "positive_compound_increases": finite["positive_compound_increases"],
            },
        ),
        CertificateRow(
            id="m100co3_06_scaled_defect_anchor",
            role="interval_anchor",
            readiness="interval_validated",
            claim="The scaled defect at the splice exceeds a simple rational floor.",
            formula="s_319=(639/2)*d_319>251/500",
            proof_boundary="One certified anchor at lambda=-100.",
            diagnostics=finite["scaled_defect_anchor"],
        ),
        CertificateRow(
            id="m100co3_07_adaptive_defect_input",
            role="exact_theorem_input",
            readiness="ready_to_apply",
            claim="The completed adaptive-defect theorem propagates decreasing d_k and increasing s_k for every index.",
            formula="d_(k+1)<=d_k and s_(k+1)>=s_k, k>=1",
            proof_boundary="All-k theorem at lambda=-100 only.",
        ),
        CertificateRow(
            id="m100co3_08_tail_increment_bound",
            role="exact_tail_lemma",
            readiness="ready_to_apply",
            claim="Scaled-defect growth and the rational anchor bound both neighboring reciprocal-defect increments.",
            formula=tail["increment_bounds"],
            proof_boundary="Tail estimate for centers k>=320 at lambda=-100.",
            diagnostics=tail,
        ),
        CertificateRow(
            id="m100co3_09_tail_compound",
            role="exact_tail_theorem",
            readiness="ready_to_apply",
            claim="The exact increment bounds leave a uniform positive compound margin on the whole tail.",
            formula=(
                "C_n>57613471/66107054971 for every n>=318 at lambda=-100"
            ),
            proof_boundary="Contiguous order-three tail only.",
            diagnostics=tail,
        ),
        CertificateRow(
            id="m100co3_10_full_entry",
            role="interval_analytic_theorem",
            readiness="ready_to_apply",
            claim="Every shifted contiguous order-three signed-Hankel minor has the required negative sign at lambda=-100.",
            formula="C_n(-100)>0 and D_(3,n)(-100)<0 for every n>=0",
            proof_boundary="All shifts, but only contiguous order-three minors at one heat parameter.",
        ),
        CertificateRow(
            id="m100co3_11_forward_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="Establish a forward-invariance theorem for the compound margin or its stronger scaled-defect sufficient cone.",
            formula="order-three entry at -100 + legitimate compound invariance => order-three sign at 0",
            proof_boundary="Open; no noncontiguous order-three theorem, all-order bridge, RH, or Lambda<=0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate",
        "date": "2026-07-12",
        "status": (
            "all-shift contiguous order-three entry theorem at lambda=-100 with "
            "open forward propagation"
        ),
        "proof_boundary": (
            "This artifact proves C_n(-100)>0 and hence the required negative "
            "contiguous 3x3 signed-Hankel sign for every shift n>=0. It does not "
            "prove forward compound-cone invariance, the sign of noncontiguous "
            "order-three minors, higher orders, the all-order signed-Hankel/Jensen "
            "bridge, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_reciprocal_defect_compound_order3_gate.md",
            "outputs/jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate.md",
            LOCAL_REPAIR.relative_to(REPO_ROOT).as_posix(),
            COLLAR_321.relative_to(REPO_ROOT).as_posix(),
        ],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.py"
        ),
        "finite": finite,
        "tail": tail,
        "summary": {
            "certificate_rows": len(rows),
            "positive_coefficients": finite["positive_coefficients"],
            "prefix_compound_rows": finite["positive_compound_margins"],
            "prefix_defect_gap_rows": finite["positive_defect_gaps"],
            "prefix_decreasing_increment_rows": finite["positive_increment_decreases"],
            "prefix_increasing_compound_rows": finite["positive_compound_increases"],
            "tail_start_n": 318,
            "tail_theorem_rows": 1,
            "full_entry_rows": 1,
            "open_forward_handoffs": 1,
            "ready_to_apply_rows": sum(row.readiness == "ready_to_apply" for row in rows),
        },
        "rows": [asdict(row) for row in rows],
    }


def write_note(path: Path, artifact: dict) -> None:
    finite = artifact["finite"]
    tail = artifact["tail"]
    lines = [
        "# Jensen-Window PF Negative-Lambda -100 Compound Order-Three Entry Certificate",
        "",
        "Date: 2026-07-12",
        "",
        "Status: all-shift contiguous order-three entry theorem at lambda=-100",
        "with open forward propagation. This is not a proof of RH or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        "validated Jensen-window PF negative-lambda -100 compound order-three entry certificate: 11 rows, 0 issues, 322 positive coefficients, 318 prefix compound margins, 318 prefix defect gaps, 1 exact tail theorem, 1 all-shift entry theorem, 1 open forward handoff",
        "```",
        "",
        "## Repaired Prefix",
        "",
        "The full repaired source is supplemented by a ten-row precision repair",
        "around `k=210` and a one-row collar at `k=321`. Direct 1024-bit Arb",
        "arithmetic proves",
        "",
        "```text",
        "C_n=q_(n+1)*q_(n+3)-q_(n+2)^2+1>0, 0<=n<=317,",
        "d_(n+2)^2-x_(n+2)^2*d_(n+1)*d_(n+3)>0, 0<=n<=317.",
        "```",
        "",
        f"The weakest compound margin is `{finite['minimum_compound_lower']}` at `n={finite['minimum_compound_at_n']}`.",
        f"The weakest defect-gap enclosure is `{finite['minimum_defect_gap_lower']}` at `n={finite['minimum_defect_gap_at_n']}`.",
        "All 319 prefix reciprocal-defect increments are positive and below one;",
        "the 318 adjacent increment differences are positive, and the 317",
        "adjacent compound-margin differences are positive. These shape facts are",
        "finite diagnostics, not assumptions in the tail proof.",
        "",
        "The splice anchor is",
        "",
        "```text",
        f"s_319=(639/2)*d_319={finite['scaled_defect_anchor']['s_ball']}",
        f"s_319-251/500>{finite['scaled_defect_anchor']['margin_lower']}.",
        "```",
        "",
        "## Exact Tail",
        "",
        "Put `k=n+2`, `m=2*k+1`, and",
        "",
        "```text",
        "u_k=q_k-q_(k-1),",
        "v_k=q_(k+1)-q_k.",
        "```",
        "",
        "The completed adaptive-defect theorem gives decreasing `d_k` and",
        "increasing `s_k=(2*k+1)*d_k/2`. Hence for every `k>=320`,",
        "",
        "```text",
        "s_j>251/500 for j>=319,",
        "q_j^2<250*(2*j+1)/251.",
        "```",
        "",
        "Scaled-defect growth also gives",
        "",
        "```text",
        "u_k<=q_k*(1-sqrt((2*k-1)/(2*k+1))).",
        "```",
        "",
        "For `m>=641`, exact squaring proves",
        "",
        "```text",
        "1-sqrt(1-2/m)<1/m+1/m^2",
        "```",
        "",
        "because the cleared condition is `m^2-2*m-1>0`. The same estimate at",
        "the next index yields",
        "",
        "```text",
        tail["increment_bounds"],
        "```",
        "",
        "Since",
        "",
        "```text",
        "C_(k-2)=1-u_k*v_k+q_k*(v_k-u_k)",
        "        >=1-u_k*v_k-q_k*u_k,",
        "```",
        "",
        "we obtain",
        "",
        "```text",
        tail["compound_lower_bound"],
        "C_n>57613471/66107054971>0, n>=318.",
        "```",
        "",
        "The endpoint comparison is exact; after setting `m=641` the rational",
        "lower margin is approximately `0.0008715177377856874`.",
        "",
        "## Theorem And Boundary",
        "",
        "Combining the Arb prefix with the analytic tail proves",
        "",
        "```text",
        "C_n(-100)>0,",
        "D_(3,n)(-100)=det[A_(n+i+j)]_(i,j=0..2)<0,",
        "for every integer n>=0.",
        "```",
        "",
        "This closes all-shift contiguous order-three entry at one heat parameter.",
        "Forward propagation to `lambda=0`, noncontiguous order-three minors, and",
        "every higher compound order remain open.",
        "",
        "```text",
        "outputs/jensen_window_pf_reciprocal_defect_compound_order3_gate.md",
        "outputs/jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate.md",
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
        "validated Jensen-window PF negative-lambda -100 compound order-three "
        "entry certificate: 11 rows, 0 issues, 322 positive coefficients, "
        "318 prefix compound margins, 318 prefix defect gaps, 1 exact tail "
        "theorem, 1 all-shift entry theorem, 1 open forward handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
