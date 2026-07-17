#!/usr/bin/env python3
"""Rule out the unit-atomic multiplier product using two Arb inequalities."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import jensen_window_pf_multiplier_leading_atom_bound_certificate as leading


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_multiplier_unit_atomic_obstruction_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_multiplier_unit_atomic_obstruction_certificate.md"
)
CUTOFF = leading.ROOT_LOWER
ORDER = 6


@dataclass(frozen=True)
class ObstructionRow:
    id: str
    role: str
    claim: str
    formula: str
    readiness: str
    proof_boundary: str
    diagnostics: dict | None = None


def build_payload() -> dict:
    leading.stress.flint.ctx.dps = leading.ARB_DPS
    differences = leading.first_alternating_differences(
        leading.load_log_contractions()
    )
    lower_gap = leading.atom_difference(CUTOFF, ORDER) - differences[ORDER]
    observed_ratio = differences[ORDER + 1] / differences[ORDER]
    atom_ratio_cap = (
        leading.atom_difference(CUTOFF, ORDER + 1)
        / leading.atom_difference(CUTOFF, ORDER)
    )
    ratio_gap = observed_ratio - atom_ratio_cap
    if not (lower_gap > 0 and not lower_gap.contains(0)):
        raise RuntimeError("unit-atom cutoff exclusion is not certified")
    if not (ratio_gap > 0 and not ratio_gap.contains(0)):
        raise RuntimeError("unit-atom ratio contradiction is not certified")

    rows = [
        ObstructionRow(
            id="muao_01_atom_kernel",
            role="exact_identity",
            claim="A unit atom contributes the positive Hausdorff difference kernel.",
            formula=(
                "g_m(alpha)=integral_0^1 r^(alpha-1)*(1-r)^(m+2)/(-log r)dr"
            ),
            readiness="available_exact",
            proof_boundary="Exact kernel identity only.",
        ),
        ObstructionRow(
            id="muao_02_unit_sum",
            role="conditional_reduction",
            claim="The target product would express every first difference moment as a unit sum.",
            formula="a_m=(-1)^m*Delta^m y_1=sum_j g_m(alpha_j)",
            readiness="conditional_on_rejected_target",
            proof_boundary="This is precisely the unit-atomic hypothesis being tested.",
        ),
        ObstructionRow(
            id="muao_03_cutoff_exclusion",
            role="interval_obstruction",
            claim="No unit atom can lie at or below the certified cutoff.",
            formula=f"g_6({CUTOFF})>a_6 implies alpha_j>{CUTOFF} for every j",
            readiness="ready_to_apply",
            proof_boundary="Uses unit multiplicity; arbitrary subunit weights are not excluded.",
            diagnostics={"gap_ball": str(lower_gap)},
        ),
        ObstructionRow(
            id="muao_04_ratio_monotonicity",
            role="exact_lemma",
            claim="The consecutive one-atom ratio is strictly decreasing in alpha.",
            formula=(
                "R_m(alpha)=g_(m+1)(alpha)/g_m(alpha)=E_alpha[1-r], "
                "R_m'(alpha)=Cov_alpha(1-r,log r)<0"
            ),
            readiness="available_exact",
            proof_boundary="Strict opposite-monotonicity covariance on 0<r<1.",
        ),
        ObstructionRow(
            id="muao_05_weighted_ratio_cap",
            role="exact_conditional_consequence",
            claim="If all target atoms exceed the cutoff, the total ratio is capped by the cutoff atom ratio.",
            formula=(
                "a_7/a_6=sum_j (g_6(alpha_j)/a_6)*R_6(alpha_j)"
                f"<R_6({CUTOFF})"
            ),
            readiness="conditional_on_rejected_target",
            proof_boundary="Positive weighted-average consequence of the target hypothesis.",
        ),
        ObstructionRow(
            id="muao_06_interval_ratio_contradiction",
            role="interval_obstruction",
            claim="The actual coefficient ratio is strictly above the required cap.",
            formula=f"a_7/a_6-R_6({CUTOFF})>0",
            readiness="ready_to_apply",
            proof_boundary="Arb contradiction at lambda=0 using A_0..A_57 sources.",
            diagnostics={
                "observed_ratio_ball": str(observed_ratio),
                "atom_ratio_cap_ball": str(atom_ratio_cap),
                "gap_ball": str(ratio_gap),
            },
        ),
        ObstructionRow(
            id="muao_07_unit_product_rejected",
            role="proved_obstruction",
            claim="The normalized zeta coefficient sequence has no admissible unit-atomic elementary multiplier product.",
            formula=(
                "not exists {alpha_j}: y_k=sum_j -log(1-1/(k+alpha_j)^2), "
                "alpha_j>0, sum_j alpha_j^(-2)<infinity"
            ),
            readiness="ready_to_apply",
            proof_boundary="Rejects this sufficient subclass, not all multiplier sequences.",
        ),
        ObstructionRow(
            id="muao_08_route_boundary",
            role="route_closure",
            claim="The multiplier counting-measure route is retired; other all-order Jensen/PF routes remain open.",
            formula="unit-atomic route closed as false; general bridge targets unchanged",
            readiness="ready_to_apply",
            proof_boundary="Not PF-infinity, the all-order Jensen bridge, RH, or Lambda <= 0.",
        ),
    ]
    summary = {
        "certificate_rows": len(rows),
        "source_difference_orders": len(differences),
        "cutoff_order": ORDER,
        "atom_cutoff": CUTOFF,
        "ratio_order": ORDER,
        "unit_atomic_routes_ruled_out": 1,
        "arbitrary_positive_weight_routes_ruled_out": 0,
        "open_requirements": 0,
        "target_closing": True,
    }
    return {
        "kind": "jensen_window_pf_multiplier_unit_atomic_obstruction_certificate",
        "date": "2026-07-11",
        "status": (
            "interval-certified unit-atomic obstruction; not a proof of PF-infinity, "
            "RH, or Lambda <= 0"
        ),
        "proof_boundary": (
            "This artifact rules out the specific convergent unit-atomic elementary "
            "multiplier product for the normalized lambda-zero zeta coefficients. It "
            "does not rule out every multiplier sequence or other all-order Jensen/PF "
            "bridges, and does not prove PF-infinity, RH, or Lambda <= 0."
        ),
        "rows": [asdict(row) for row in rows],
        "summary": summary,
        "arb_dps": leading.ARB_DPS,
        "source_enclosures": [
            path.relative_to(REPO_ROOT).as_posix() for path in leading.SOURCES
        ],
        "source_leading_bound": (
            "outputs/jensen_window_pf_multiplier_leading_atom_bound_certificate.md"
        ),
        "source_hausdorff_bridge": (
            "outputs/jensen_window_pf_multiplier_hausdorff_uniqueness_bridge.md"
        ),
        "source_target": "outputs/jensen_window_pf_multiplier_counting_measure_target.md",
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_multiplier_unit_atomic_obstruction_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_multiplier_unit_atomic_obstruction_certificate.py"
        ),
    }


def write_note(payload: dict, path: Path) -> None:
    rows = {row["id"]: row for row in payload["rows"]}
    lines = [
        "# Jensen-Window PF Multiplier Unit-Atomic Obstruction Certificate",
        "",
        "Date: 2026-07-11",
        "",
        "Status: interval-certified unit-atomic obstruction; not a proof of",
        "PF-infinity, the all-order Jensen bridge, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_multiplier_unit_atomic_obstruction_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_multiplier_unit_atomic_obstruction_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_multiplier_unit_atomic_obstruction_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_multiplier_unit_atomic_obstruction_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        "validated Jensen-window PF multiplier unit-atomic obstruction certificate: 8 rows, 0 issues, 1 atom cutoff, 1 ratio cap, 1 unit-atomic route ruled out, 0 open requirements",
        "```",
        "",
        "## Target Consequences",
        "",
        "If the proposed unit counting multiset existed, then",
        "",
        "```text",
        "a_m=(-1)^m*Delta^m y_1=sum_j g_m(alpha_j),",
        "g_m(alpha)=integral_0^1 r^(alpha-1)*(1-r)^(m+2)/(-log r)dr.",
        "```",
        "",
        "The kernel is positive and strictly decreasing in `alpha`. Arb proves",
        "",
        "```text",
        f"g_6({CUTOFF})-a_6={rows['muao_03_cutoff_exclusion']['diagnostics']['gap_ball']}>0.",
        "```",
        "",
        f"Because every target atom has unit weight, this forces `alpha_j>{CUTOFF}`",
        "for every `j`.",
        "",
        "## Ratio Contradiction",
        "",
        "Put `R_m(alpha)=g_(m+1)(alpha)/g_m(alpha)`. Under the probability",
        "density proportional to the integrand defining `g_m`,",
        "",
        "```text",
        "R_m(alpha)=E_alpha[1-r],",
        "R_m'(alpha)=Cov_alpha(1-r,log r)<0.",
        "```",
        "",
        "The covariance is strictly negative because `1-r` decreases and `log r`",
        "increases. Therefore the target would require",
        "",
        "```text",
        f"a_7/a_6<R_6({CUTOFF}).",
        "```",
        "",
        "The actual intervals give",
        "",
        "```text",
        f"a_7/a_6={rows['muao_06_interval_ratio_contradiction']['diagnostics']['observed_ratio_ball']},",
        f"R_6({CUTOFF})={rows['muao_06_interval_ratio_contradiction']['diagnostics']['atom_ratio_cap_ball']},",
        f"gap={rows['muao_06_interval_ratio_contradiction']['diagnostics']['gap_ball']}>0.",
        "```",
        "",
        "This is incompatible with the required weighted-average cap. Hence the",
        "normalized zeta coefficients do not admit the proposed convergent",
        "unit-atomic elementary multiplier product.",
        "",
        "## Boundary",
        "",
        "The contradiction uses unit multiplicity in the cutoff step. It does not rule",
        "out arbitrary positive subunit weights, but such fractional weights were",
        "already insufficient for multiplier preservation. More general multiplier",
        "sequences and the other all-order Jensen/PF routes remain open.",
        "",
        "```text",
        "outputs/jensen_window_pf_multiplier_leading_atom_bound_certificate.md",
        "outputs/jensen_window_pf_multiplier_hausdorff_uniqueness_bridge.md",
        "outputs/jensen_window_pf_multiplier_counting_measure_target.md",
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
    print(
        "validated Jensen-window PF multiplier unit-atomic obstruction certificate: "
        "8 rows, 0 issues, 1 atom cutoff, 1 ratio cap, "
        "1 unit-atomic route ruled out, 0 open requirements"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
