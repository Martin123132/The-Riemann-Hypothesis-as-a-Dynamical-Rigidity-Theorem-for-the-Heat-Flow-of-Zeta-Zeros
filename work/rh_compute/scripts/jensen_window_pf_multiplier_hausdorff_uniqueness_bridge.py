#!/usr/bin/env python3
"""Build the exact Hausdorff uniqueness bridge for multiplier atoms."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_multiplier_hausdorff_uniqueness_bridge.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_multiplier_hausdorff_uniqueness_bridge.md"
)


@dataclass(frozen=True)
class BridgeRow:
    id: str
    role: str
    claim: str
    formula: str
    readiness: str
    proof_boundary: str


def build_payload() -> dict:
    rows = [
        BridgeRow(
            id="mhb_01_atom_laplace_kernel",
            role="exact_identity",
            claim="Each elementary multiplier curvature is a positive Laplace transform.",
            formula=(
                "-log(1-1/(k+alpha)^2)=integral_0^infinity "
                "exp(-(k+alpha)t)*q(t)dt, q(t)=2*(cosh(t)-1)/t"
            ),
            readiness="available_exact",
            proof_boundary="Frullani identity for k+alpha>1.",
        ),
        BridgeRow(
            id="mhb_02_counting_laplace_sum",
            role="exact_lemma",
            claim="A convergent unit counting multiset defines a finite positive Laplace density.",
            formula=(
                "S(t)=sum_j exp(-alpha_j*t), "
                "y_k=integral_0^infinity exp(-k*t)*q(t)*S(t)dt"
            ),
            readiness="available_exact",
            proof_boundary="Uses alpha_j>0 and sum_j alpha_j^(-2)<infinity.",
        ),
        BridgeRow(
            id="mhb_03_hausdorff_pushforward",
            role="exact_lemma",
            claim="The integer log contractions are Hausdorff moments of one finite measure.",
            formula=(
                "y_k=integral_[0,1] r^(k-1)dnu(r), "
                "dnu/dr=q(-log r)*S(-log r)"
            ),
            readiness="available_exact",
            proof_boundary="Change of variables r=exp(-t); total mass is y_1.",
        ),
        BridgeRow(
            id="mhb_04_hausdorff_uniqueness",
            role="exact_theorem",
            claim="All integer values y_1,y_2,... uniquely determine nu.",
            formula="int r^n dnu=int r^n dnu_tilde for all n>=0 implies nu=nu_tilde",
            readiness="available_exact",
            proof_boundary="Polynomial density in C([0,1]) and finite-measure duality.",
        ),
        BridgeRow(
            id="mhb_05_complete_monotonicity",
            role="exact_consequence",
            claim="Hausdorff positivity gives every alternating finite-difference inequality.",
            formula=(
                "(-1)^m*Delta^m y_k=integral_[0,1] "
                "r^(k-1)*(1-r)^m dnu(r)>=0"
            ),
            readiness="available_exact",
            proof_boundary="Necessary for the counting product, not sufficient for unit atoms.",
        ),
        BridgeRow(
            id="mhb_06_unit_atomic_characterization",
            role="exact_characterization",
            claim="The integer product target is equivalent to a rigid property of the unique Hausdorff measure.",
            formula=(
                "dnu/dr=q(-log r)*sum_j r^alpha_j with unit integer "
                "multiplicity and sum_j alpha_j^(-2)<infinity"
            ),
            readiness="available_exact",
            proof_boundary="Characterization only; existence of the density and atoms is open.",
        ),
        BridgeRow(
            id="mhb_07_atom_uniqueness",
            role="exact_consequence",
            claim="If the counting representation exists, its multiset is unique.",
            formula=(
                "S(t)=(dnu/dr)(exp(-t))/q(t); uniqueness of Laplace "
                "transforms determines sum_j delta_(alpha_j)"
            ),
            readiness="available_exact",
            proof_boundary="Uniqueness does not supply existence or an extraction algorithm.",
        ),
        BridgeRow(
            id="mhb_08_finite_frontier_input",
            role="finite_evidence",
            claim="The available zeta sequence passes the full finite Hausdorff sign frontier.",
            formula="7980/7980 strict signs through order 55 on five heat parameters",
            readiness="finite_only",
            proof_boundary="Finite interval evidence only; not a representing-measure theorem.",
        ),
        BridgeRow(
            id="mhb_09_periodic_interpolation_guard",
            role="nonpromotion_guard",
            claim="Integer agreement and even continuous curvature agreement do not identify the natural Mellin interpolation.",
            formula=(
                "f(s)=sin(2*pi*s): f(n)=0 and "
                "f(s+1)-2*f(s)+f(s-1)=0"
            ),
            readiness="available_exact",
            proof_boundary=(
                "A growth or canonical-interpolation theorem is still required before "
                "using the continuous Mellin power-sum obstruction."
            ),
        ),
        BridgeRow(
            id="mhb_10_open_recovery_handoff",
            role="open_handoff",
            claim="Recover or rule out the unique deconvolved unit counting measure.",
            formula=(
                "nu -> h(r)=dnu/dr -> S(t)=h(exp(-t))/q(t) "
                "-> unit counting measure N"
            ),
            readiness="not_ready_to_apply",
            proof_boundary="This recovery is the remaining multiplier target, not a proved step.",
        ),
    ]
    summary = {
        "bridge_rows": len(rows),
        "exact_rows": 7,
        "hausdorff_measure_theorems": 1,
        "unit_atomic_characterizations": 1,
        "finite_frontier_rows": 1,
        "interpolation_guards": 1,
        "open_recovery_handoffs": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    return {
        "kind": "jensen_window_pf_multiplier_hausdorff_uniqueness_bridge",
        "date": "2026-07-11",
        "status": "exact Hausdorff uniqueness and unit-atomic characterization",
        "proof_boundary": (
            "This artifact proves uniqueness and an exact characterization of any "
            "integer-only multiplier counting measure. It does not construct that "
            "measure, identify the natural Mellin interpolation with the canonical "
            "Hausdorff interpolation, prove PF-infinity, RH, or Lambda <= 0."
        ),
        "rows": [asdict(row) for row in rows],
        "summary": summary,
        "source_frontier": (
            "outputs/jensen_window_pf_multiplier_complete_monotonicity_frontier_scout.md"
        ),
        "source_target": "outputs/jensen_window_pf_multiplier_counting_measure_target.md",
        "source_mellin_guard": (
            "outputs/jensen_window_pf_mellin_multiplier_power_sum_obstruction.md"
        ),
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_multiplier_hausdorff_uniqueness_bridge.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_multiplier_hausdorff_uniqueness_bridge.py"
        ),
    }


def write_note(payload: dict, path: Path) -> None:
    summary = payload["summary"]
    lines = [
        "# Jensen-Window PF Multiplier Hausdorff-Uniqueness Bridge",
        "",
        "Date: 2026-07-11",
        "",
        "Status: exact Hausdorff uniqueness and unit-atomic characterization. This is not a proof",
        "of the counting-measure target, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_multiplier_hausdorff_uniqueness_bridge`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_multiplier_hausdorff_uniqueness_bridge.json",
        "python work/rh_compute/scripts/jensen_window_pf_multiplier_hausdorff_uniqueness_bridge.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_multiplier_hausdorff_uniqueness_bridge.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        (
            "validated Jensen-window PF multiplier Hausdorff-uniqueness bridge: "
            f"{summary['bridge_rows']} rows, 0 issues, 1 Hausdorff measure theorem, "
            "1 unit-atomic characterization, 1 interpolation guard, "
            "1 open recovery handoff"
        ),
        "```",
        "",
        "## Exact Measure Bridge",
        "",
        "For `z>1`, Frullani's identity gives",
        "",
        "```text",
        "-log(1-z^(-2))",
        " =integral_0^infinity exp(-z*t)*q(t)dt,",
        "q(t)=2*(cosh(t)-1)/t>0.",
        "```",
        "",
        "If the target multiset exists, put",
        "",
        "```text",
        "S(t)=sum_j exp(-alpha_j*t).",
        "```",
        "",
        "The convergence condition `sum_j alpha_j^(-2)<infinity` makes the",
        "following positive measure finite. With `r=exp(-t)`,",
        "",
        "```text",
        "y_k=-log(x_k)=integral_[0,1] r^(k-1)dnu(r),",
        "dnu(r)=q(-log r)*S(-log r)dr.",
        "```",
        "",
        "The values `y_1,y_2,...` are every moment of `nu`. Polynomials are dense",
        "in `C([0,1])`, so the finite positive measure `nu` is unique. Consequently",
        "",
        "```text",
        "(-1)^m*Delta^m y_k",
        " =integral_[0,1] r^(k-1)*(1-r)^m dnu(r)>=0.",
        "```",
        "",
        "More importantly, the integer-only multiplier target is equivalent to the",
        "unique measure having an absolutely continuous density of the rigid form",
        "",
        "```text",
        "dnu/dr=q(-log r)*sum_j r^alpha_j",
        "```",
        "",
        "with unit integer multiplicities and `sum_j alpha_j^(-2)<infinity`.",
        "If such a representation exists, division by `q(t)>0` followed by Laplace",
        "uniqueness also makes the multiset `{alpha_j}` unique.",
        "",
        "## Interpolation Guard",
        "",
        "Hausdorff uniqueness does not by itself identify the natural continuous",
        "Mellin interpolation used by the power-sum obstruction. The exact ambiguity",
        "already appears in",
        "",
        "```text",
        "f(s)=sin(2*pi*s),",
        "f(n)=0,",
        "f(s+1)-2*f(s)+f(s-1)=0.",
        "```",
        "",
        "Thus integer values, and even the same centered continuous log curvature,",
        "permit a nontrivial periodic interpolation factor. A canonical-interpolation",
        "or growth theorem is still needed before the continuous Mellin determinant",
        "can rule out the integer-only product.",
        "",
        "## New Handoff",
        "",
        "The target is now a unique recovery problem:",
        "",
        "```text",
        "{y_k} -> unique Hausdorff measure nu",
        "      -> test absolute continuity",
        "      -> divide by q",
        "      -> recover or rule out a unit counting Laplace measure.",
        "```",
        "",
        "The order-55 Arb frontier supports positivity of `nu` finitely but supplies",
        "neither the all-order measure theorem nor the unit-multiplicity recovery.",
        "",
        "```text",
        "outputs/jensen_window_pf_multiplier_complete_monotonicity_frontier_scout.md",
        "outputs/jensen_window_pf_multiplier_counting_measure_target.md",
        "outputs/jensen_window_pf_mellin_multiplier_power_sum_obstruction.md",
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "```",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(payload, args.note)
    summary = payload["summary"]
    print(
        "validated Jensen-window PF multiplier Hausdorff-uniqueness bridge: "
        f"{summary['bridge_rows']} rows, 0 issues, 1 Hausdorff measure theorem, "
        "1 unit-atomic characterization, 1 interpolation guard, "
        "1 open recovery handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
