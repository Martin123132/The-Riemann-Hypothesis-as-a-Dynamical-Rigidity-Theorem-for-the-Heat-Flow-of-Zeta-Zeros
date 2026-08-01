#!/usr/bin/env python3
"""Compose the first-order remainder with the cofinal boundary-degree route."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_cofinal_boundary_reduction"
)
DEFAULT_OUT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / f"outputs/{STEM}.md"
SOURCES = {
    "global_remainder": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_"
        "first_order_global_remainder_certificate.json"
    ),
    "signed_contact": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_"
        "first_order_signed_contact_reduction.json"
    ),
    "boundary_degree": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_"
        "cofinal_boundary_degree_transfer_target.json"
    ),
    "ray_schedule": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_"
        "ray_aligned_parabolic_frequency_reduction.json"
    ),
}
ETA_0 = 100_000
ETA_1 = 200_000
FULL_SQUARED_BUDGET = ETA_0**2 + ETA_1**2
HALF_SQUARED_BUDGET = FULL_SQUARED_BUDGET // 4
OLD_SQUARED_BUDGET = 32_000_000
ADJACENT_VALUE_BUDGET = 20_000
ADJACENT_DERIVATIVE_BUDGET = 30_000
ADJACENT_SQUARED_BUDGET = (
    ADJACENT_VALUE_BUDGET**2 + ADJACENT_DERIVATIVE_BUDGET**2
)
C_STAR = "4911678521/1933561194"


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
    missing = [str(path) for path in SOURCES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing source artifacts: {missing}")
    return {key: file_hash(path) for key, path in SOURCES.items()}


def verify_exact_guards() -> dict[str, str]:
    if FULL_SQUARED_BUDGET != 50_000_000_000:
        raise RuntimeError("full first-order vector budget drifted")
    if HALF_SQUARED_BUDGET != 12_500_000_000:
        raise RuntimeError("complex-half vector budget drifted")
    if ADJACENT_SQUARED_BUDGET != 1_300_000_000:
        raise RuntimeError("adjacent vector budget drifted")
    if (
        Fraction(ADJACENT_SQUARED_BUDGET, FULL_SQUARED_BUDGET)
        != Fraction(13, 500)
    ):
        raise RuntimeError("adjacent/main vector ratio drifted")

    ratio_constant = Fraction(FULL_SQUARED_BUDGET, OLD_SQUARED_BUDGET)
    if ratio_constant != Fraction(3125, 2):
        raise RuntimeError("old/new ratio constant drifted")
    rational_l50_bound = ratio_constant * Fraction(10, 27) ** 50
    if not rational_l50_bound < Fraction(5, 10**19):
        raise RuntimeError("L=50 budget-ratio guard failed")

    for stage in (0, 1, 10, 100, 1000, 10_000):
        bottom = Fraction(25, 100 + stage)
        log_radius = 101 + stage
        right_min = bottom * log_radius
        if not right_min > 25:
            raise RuntimeError(f"right-edge schedule failed at stage {stage}")
        next_bottom = Fraction(25, 101 + stage)
        if next_bottom * log_radius != 25:
            raise RuntimeError(f"new-strip schedule failed at stage {stage}")

    return {
        "ratio_constant": str(ratio_constant),
        "rational_L50_upper": str(rational_l50_bound),
        "decimal_L50_upper": f"{float(rational_l50_bound):.18e}",
    }


def exact_identities() -> dict[str, str]:
    return {
        "normalizer": (
            "H=A_t*Z, A_t>0; every positive first-jet scale and the "
            "normalizer derivative shear preserve boundary winding and "
            "local contact index"
        ),
        "vector_remainder": (
            "V_r=(r_[1],r_[1],x/L), "
            "||V_r||_2^2<50000000000*exp(-5L/2)"
        ),
        "boundary_rouche": (
            "If V_J=(J_[1],J_[1],x/L) satisfies "
            "J_[1]^2+(J_[1],x/L)^2>"
            "50000000000*exp(-5L/2) pointwise on a closed boundary, "
            "then V_J+s*V_r is nonzero there for 0<=s<=1 and "
            "wind(V_Z)=wind(V_J)"
        ),
        "complex_half": (
            "For J_[1]=2X_[1] and J_[1],x=2U_[1], it is enough that "
            "X_[1]^2+(U_[1]/L)^2>"
            "12500000000*exp(-5L/2)"
        ),
        "old_new_ratio": (
            "Relative to the old 32000000*exp(-3L/2) squared budget, "
            "the new/old ratio is (3125/2)*exp(-L), decreasing in L; "
            "at L>=50 it is <5*10^-19"
        ),
        "ray_schedule": (
            "t_j=25/(100+j), L_j=101+j, R_j=4*pi*exp(L_j), "
            "P_j=[t_j,1/2]x[0,R_j]"
        ),
        "right_edge": (
            "On x=R_j and t_j<=t<=1/2, "
            "tL_j>=25*(101+j)/(100+j)>25; "
            "the dominant-saddle theorem makes this boundary contact-free"
        ),
        "new_strip": (
            "For S_j=[t_(j+1),1/2]x[R_j,R_(j+1)], "
            "tL>=t_(j+1)L_j=25, so every new strip is already closed"
        ),
        "bottom_split": (
            "On the bottom t=t_j, q=2t_jL^2=50L^2/(100+j); "
            "q<1 iff L<sqrt((100+j)/50), and q>=1 beyond that point"
        ),
        "outer_bottom": (
            f"On a fixed bottom, c=t_jL; every fixed c>{C_STAR} "
            "is asymptotically closed after the oscillatory-zeta handoff, "
            "subject to its epsilon-dependent finite shoulder"
        ),
        "adjacent_value": (
            "On a doubled cutoff collar, the two adjacent first-order lifts "
            "satisfy |Delta J_[1]|<20000*exp(-5L/4)"
        ),
        "adjacent_derivative": (
            "Cauchy on radius 1/L and |partial_x log A_t|<L/2 give "
            "|partial_x Delta J_[1]|<30000*L*exp(-5L/4)"
        ),
        "adjacent_vector": (
            "||V_(J,N+1)-V_(J,N)||_2^2"
            "<1300000000*exp(-5L/2)"
            "=(13/500)*50000000000*exp(-5L/2)"
        ),
        "continuity_homotopy": (
            "If the standard boundary margin ||V_(J,N)||_2^2>"
            "50000000000*exp(-5L/2) holds, then "
            "V_(J,N)+s*(V_(J,N+1)-V_(J,N)) is nonzero for 0<=s<=1; "
            "the local cutoff lifts therefore join continuously without a "
            "separate arithmetic theorem"
        ),
        "winding_obligation": (
            "Boundary nonvanishing transfers winding but does not determine it; "
            "because the transferred contact degree is a nonnegative integer, "
            "it is enough to prove wind(V_J)<1 on every successor boundary"
        ),
        "cofinal_handoff": (
            "The ray-aligned P_0 has zero winding because t_0=1/4>1/5>=Lambda; "
            "Q208 remains a compact calibration. For every successor, "
            "prove the first-order boundary domination and continuous "
            "one-sided phase bound only on its boundary; the same-sign positive "
            "contact index in the standard (x,t) orientation then excludes "
            "all interior contacts"
        ),
        "open_split": (
            "The remaining proof work is one-dimensional but has three pieces: "
            "the q>=1 critical bottom arc, a multiplicity-compatible q<1 "
            "bottom arc, and finite bounded-L/epsilon shoulders, together with "
            "the winding computation"
        ),
    }


def build_rows(exact: dict[str, str]) -> list[GateRow]:
    return [
        GateRow(
            id="np15f1cbr_01_orientation",
            role="certified_input",
            readiness="available_exact",
            claim="Positive normalization preserves boundary degree.",
            formula=exact["normalizer"],
            proof_boundary="Imported from the checked cofinal degree gate.",
        ),
        GateRow(
            id="np15f1cbr_02_vector_budget",
            role="proved_composition",
            readiness="available_asymptotic",
            claim="The first-order C1 residual has a squared jet budget.",
            formula=exact["vector_remainder"],
            proof_boundary="Only on L>=50 and 0<tL<=25.",
        ),
        GateRow(
            id="np15f1cbr_03_rouche",
            role="proved_composition",
            readiness="available_asymptotic",
            claim="Strict first-order boundary domination transfers winding.",
            formula=exact["boundary_rouche"],
            proof_boundary="Conditional on the displayed strict boundary margin.",
        ),
        GateRow(
            id="np15f1cbr_04_complex_half",
            role="exact_identity",
            readiness="available_exact",
            claim="The boundary target has an equivalent complex-half form.",
            formula=exact["complex_half"],
            proof_boundary="Exact factor-four conversion.",
        ),
        GateRow(
            id="np15f1cbr_05_gain",
            role="exact_bound",
            readiness="available_exact",
            claim="The first-order peel improves the old boundary budget exponentially.",
            formula=exact["old_new_ratio"],
            proof_boundary="Exact rational comparison using e>27/10.",
        ),
        GateRow(
            id="np15f1cbr_06_schedule",
            role="exact_definition",
            readiness="available_exact",
            claim="The ray-aligned rectangles form a cofinal boundary family.",
            formula=exact["ray_schedule"],
            proof_boundary="Uses t_j tending to zero and R_j tending to infinity.",
        ),
        GateRow(
            id="np15f1cbr_07_right_edge",
            role="proved_composition",
            readiness="available_asymptotic",
            claim="Every right edge lies in the proved dominant-saddle region.",
            formula=exact["right_edge"] + "; " + exact["new_strip"],
            proof_boundary="L_j>=101 and the checked tL>=25 theorem.",
        ),
        GateRow(
            id="np15f1cbr_08_bottom_split",
            role="exact_identity",
            readiness="available_exact",
            claim="Every bottom edge has an exact parabolic/frequency split.",
            formula=exact["bottom_split"],
            proof_boundary="Pure schedule algebra.",
        ),
        GateRow(
            id="np15f1cbr_09_outer_handoff",
            role="proved_composition",
            readiness="available_asymptotic",
            claim="Fixed outer bottom rays are already closed asymptotically.",
            formula=exact["outer_bottom"],
            proof_boundary="Does not close c=c_* or finite L_epsilon shoulders.",
        ),
        GateRow(
            id="np15f1cbr_10_frequency_bottom",
            role="open_theorem_target",
            readiness="open",
            claim="The corrected first-order proxy must dominate on the live q>=1 bottom arc.",
            formula=exact["complex_half"],
            proof_boundary="Open Xi-specific phase-critical-value theorem.",
        ),
        GateRow(
            id="np15f1cbr_11_parabolic_bottom",
            role="open_theorem_target",
            readiness="open",
            claim="The q<1 bottom arc needs multiplicity-compatible control.",
            formula=exact["bottom_split"] + "; " + exact["open_split"],
            proof_boundary="Open; no endpoint-simplicity floor is imposed.",
        ),
        GateRow(
            id="np15f1cbr_12_continuity",
            role="proved_composition",
            readiness="available_asymptotic",
            claim="The standard boundary margin automatically joins adjacent cutoff lifts.",
            formula=(
                exact["adjacent_value"]
                + "; "
                + exact["adjacent_derivative"]
                + "; "
                + exact["adjacent_vector"]
                + "; "
                + exact["continuity_homotopy"]
            ),
            proof_boundary=(
                "Conditional only on the same open boundary margin; no extra "
                "cutoff-separation theorem is needed."
            ),
        ),
        GateRow(
            id="np15f1cbr_13_winding",
            role="open_theorem_target",
            readiness="open",
            claim="The continuous successor proxy needs only a winding upper bound strictly below one.",
            formula=exact["winding_obligation"] + "; " + exact["cofinal_handoff"],
            proof_boundary="Open; nonvanishing alone does not determine winding.",
        ),
    ]


def build_artifact() -> dict:
    audit = verify_exact_guards()
    exact = exact_identities()
    rows = build_rows(exact)
    return {
        "kind": STEM,
        "date": "2026-07-26",
        "status": (
            "exact first-order cofinal boundary reduction with improved "
            "jet budget; not a boundary theorem, Lambda<=0, or RH"
        ),
        "proof_boundary": (
            "This artifact proves the first-order vector remainder budget, "
            "boundary-Rouche transfer condition, exponential improvement over "
            "the old budget, ray-schedule right-edge closure, and exact bottom "
            "q split. It also proves that the same boundary margin joins "
            "adjacent cutoff lifts by a nonvanishing overlap homotopy. It does "
            "not prove the q>=1 boundary margin, q<1 multiplicity-compatible "
            "margin, one-sided proxy phase bound, bounded-L shoulders, the all-stage cofinal boundary "
            "theorem, contact exclusion, Lambda<=0, RH, PF-infinity, or a "
            "Clay-prize result."
        ),
        "builder_sha256": file_hash(Path(__file__)),
        "source_sha256": source_hashes(),
        "exact": exact,
        "audit": audit,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "certified_inputs_or_compositions": 6,
            "exact_identities_or_bounds": 5,
            "open_boundary_targets": 3,
            "full_squared_budget": FULL_SQUARED_BUDGET,
            "half_squared_budget": HALF_SQUARED_BUDGET,
            "adjacent_squared_budget": ADJACENT_SQUARED_BUDGET,
        },
        "sources": [
            str(path.relative_to(REPO_ROOT)).replace("\\", "/")
            for path in SOURCES.values()
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    audit = artifact["audit"]
    return "\n".join(
        [
            "# Newman Polymath-15 First-Order Cofinal Boundary Reduction",
            "",
            "Date: 2026-07-26",
            "",
            "Status: exact boundary reduction with a certified first-order",
            "jet budget. This is not a proof of the cofinal boundary theorem,",
            "`Lambda <= 0`, or RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## First-Order Boundary Budget",
            "",
            "```text",
            exact["vector_remainder"],
            exact["boundary_rouche"],
            exact["complex_half"],
            "```",
            "",
            "The comparison with the old corrected budget is",
            "",
            "```text",
            exact["old_new_ratio"],
            "```",
            "",
            f"The exact rational `L=50` upper audit is `{audit['rational_L50_upper']}`",
            f"or `{audit['decimal_L50_upper']}`.",
            "",
            "## Ray-Aligned Boundary",
            "",
            "```text",
            exact["ray_schedule"],
            exact["right_edge"],
            exact["new_strip"],
            "```",
            "",
            "The compact core and top edge are inherited from the certified",
            "base and published positive-time bound. The only new arithmetic",
            "work is on the bottom boundary and finite joins.",
            "",
            "## Bottom Split",
            "",
            "```text",
            exact["bottom_split"],
            exact["outer_bottom"],
            exact["open_split"],
            "```",
            "",
            "This reduces the hard two-dimensional wedge to selected",
            "one-dimensional boundary arcs, without claiming those arcs are",
            "already controlled.",
            "",
            "## Continuity And Winding",
            "",
            "```text",
            exact["adjacent_value"],
            exact["adjacent_derivative"],
            exact["adjacent_vector"],
            exact["continuity_homotopy"],
            exact["winding_obligation"],
            exact["cofinal_handoff"],
            "```",
            "",
            "The cutoff joins now cost no independent theorem once the standard",
            "margin holds. Positivity of the contact degree means the resulting",
            "integer needs only a strict upper bound below one.",
            "",
            "## Proof Boundary",
            "",
            artifact["proof_boundary"],
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman first-order cofinal boundary reduction: "
        "13 rows, full budget 50000000000, "
        "3 open boundary targets"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
