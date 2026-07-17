#!/usr/bin/env python3
"""Certify conditional leading-atom bounds for the multiplier product."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path

import jensen_window_pf_monotone_contraction_stress as stress


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULTS = REPO_ROOT / "work/rh_compute/results"
DEFAULT_OUT = RESULTS / "jensen_window_pf_multiplier_leading_atom_bound_certificate.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_multiplier_leading_atom_bound_certificate.md"
SOURCES = tuple(
    RESULTS / f"acb_enclosures_edrei_boundary_lam0_{chunk}_dps220_tol1e-140.jsonl"
    for chunk in ("k0_k53", "k54_k55", "k56_k57")
)
ARB_DPS = 250
ROOT_ORDER = 6
ROOT_LOWER = "4.863538496"
ROOT_UPPER = "4.863538497"
COUNT_CUTOFF = "5.5"
COUNT_ORDER = 3


@dataclass(frozen=True)
class CertificateRow:
    id: str
    role: str
    claim: str
    formula: str
    readiness: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_log_contractions():
    for path in SOURCES:
        if not path.exists():
            raise FileNotFoundError(path)
    balls, _samples, labels = stress.load_enclosures(list(SOURCES))
    if len(labels) != 1:
        raise RuntimeError(f"expected one lambda, found {len(labels)}")
    lam = next(iter(labels))
    coefficients = [balls[(lam, k)] for k in range(58)]
    return [-value.log() for value in stress.contractions(coefficients)]


def first_alternating_differences(values) -> list:
    current = list(values)
    out = []
    for order in range(len(current)):
        out.append((1 if order % 2 == 0 else -1) * current[0])
        current = [
            current[index + 1] - current[index]
            for index in range(len(current) - 1)
        ]
    return out


def atom_difference(alpha: str, order: int):
    arb = stress.flint.arb
    value = arb(0)
    a = arb(alpha)
    for shift in range(order + 1):
        z = arb(1 + shift) + a
        kernel = -(1 - 1 / (z * z)).log()
        value += ((-1) ** shift) * math.comb(order, shift) * kernel
    return value


def build_payload() -> dict:
    stress.flint.ctx.dps = ARB_DPS
    differences = first_alternating_differences(load_log_contractions())
    if len(differences) != 56:
        raise RuntimeError(f"expected 56 difference orders, found {len(differences)}")
    if not all(value > 0 and not value.contains(0) for value in differences):
        raise RuntimeError("source alternating differences are not all strictly positive")

    lower_gap = atom_difference(ROOT_LOWER, ROOT_ORDER) - differences[ROOT_ORDER]
    upper_gap = atom_difference(ROOT_UPPER, ROOT_ORDER) - differences[ROOT_ORDER]
    if not (lower_gap > 0 and not lower_gap.contains(0)):
        raise RuntimeError("root lower endpoint is not certified below beta6")
    if not (upper_gap < 0 and not upper_gap.contains(0)):
        raise RuntimeError("root upper endpoint is not certified above beta6")

    other_gaps = {
        order: atom_difference(ROOT_LOWER, order) - differences[order]
        for order in range(len(differences))
        if order != ROOT_ORDER
    }
    if not all(gap < 0 and not gap.contains(0) for gap in other_gaps.values()):
        raise RuntimeError("another finite order is not certified weaker than order 6")
    weakest_other_order = max(other_gaps, key=lambda order: other_gaps[order].upper())

    count_kernel = atom_difference(COUNT_CUTOFF, COUNT_ORDER)
    count_ratio = differences[COUNT_ORDER] / count_kernel
    count_exclusion_gap = 2 * count_kernel - differences[COUNT_ORDER]
    if not (count_exclusion_gap > 0 and not count_exclusion_gap.contains(0)):
        raise RuntimeError("two-atom cutoff exclusion is not certified")

    rows = [
        CertificateRow(
            id="mlab_01_atom_difference_kernel",
            role="exact_identity",
            claim="Each atom contributes a positive finite-difference kernel.",
            formula=(
                "g_m(alpha)=integral_0^1 r^(alpha-1)*(1-r)^(m+2)/(-log r)dr"
            ),
            readiness="available_exact",
            proof_boundary="Exact consequence of the Hausdorff density formula.",
        ),
        CertificateRow(
            id="mlab_02_kernel_monotonicity",
            role="exact_lemma",
            claim="The atom kernel is strictly decreasing in alpha.",
            formula="partial_alpha g_m(alpha)=-Beta(alpha,m+3)<0",
            readiness="available_exact",
            proof_boundary="Exact monotonicity only.",
        ),
        CertificateRow(
            id="mlab_03_conditional_sum",
            role="exact_conditional_reduction",
            claim="Under the multiplier target, every first difference moment is the sum of atom kernels.",
            formula="a_m=(-1)^m*Delta^m y_1=sum_j g_m(alpha_j)",
            readiness="conditional_on_target",
            proof_boundary="Uses the still-open unit-atomic factorization hypothesis.",
        ),
        CertificateRow(
            id="mlab_04_root_lower_endpoint",
            role="interval_certificate",
            claim="Order 6 places the root beta6 above the rational decimal lower endpoint.",
            formula=f"g_6({ROOT_LOWER})-a_6>0",
            readiness="ready_to_apply_conditionally",
            proof_boundary="Conditional lower bound for any target atom multiset.",
            diagnostics={"gap_ball": str(lower_gap)},
        ),
        CertificateRow(
            id="mlab_05_root_upper_endpoint",
            role="interval_certificate",
            claim="Order 6 places beta6 below the rational decimal upper endpoint.",
            formula=f"g_6({ROOT_UPPER})-a_6<0",
            readiness="ready_to_apply_conditionally",
            proof_boundary="Brackets beta6, not the unknown smallest atom itself.",
            diagnostics={"gap_ball": str(upper_gap)},
        ),
        CertificateRow(
            id="mlab_06_strongest_finite_order",
            role="interval_certificate",
            claim="Every other available order gives a strictly weaker root lower bound.",
            formula=f"g_m({ROOT_LOWER})-a_m<0 for 0<=m<=55, m!=6",
            readiness="ready_to_apply_conditionally",
            proof_boundary="Strongest only within the certified finite order-55 frontier.",
            diagnostics={
                "checked_other_orders": len(other_gaps),
                "weakest_other_order": weakest_other_order,
                "weakest_other_gap_ball": str(other_gaps[weakest_other_order]),
            },
        ),
        CertificateRow(
            id="mlab_07_low_atom_count",
            role="interval_certificate",
            claim="At most one target atom can lie at or below 11/2.",
            formula="N(11/2)*g_3(11/2)<=a_3 and a_3/g_3(11/2)<2",
            readiness="ready_to_apply_conditionally",
            proof_boundary="Does not prove that an atom actually lies below 11/2.",
            diagnostics={
                "ratio_ball": str(count_ratio),
                "two_atom_exclusion_gap_ball": str(count_exclusion_gap),
            },
        ),
        CertificateRow(
            id="mlab_08_open_existence_handoff",
            role="open_handoff",
            claim="Existence and location of a first atom remain open.",
            formula=(
                f"if target holds: alpha_min>{ROOT_LOWER} and N(11/2)<=1; "
                "no atom below 11/2 is proved"
            ),
            readiness="not_ready_to_apply_unconditionally",
            proof_boundary="This certificate constrains but does not construct the target measure.",
        ),
    ]
    summary = {
        "certificate_rows": len(rows),
        "difference_orders": len(differences),
        "root_order": ROOT_ORDER,
        "root_bracket": [ROOT_LOWER, ROOT_UPPER],
        "conditional_alpha_min_lower_bound": ROOT_LOWER,
        "other_orders_certified_weaker": len(other_gaps),
        "count_cutoff": "11/2",
        "count_order": COUNT_ORDER,
        "max_atoms_at_or_below_cutoff": 1,
        "open_existence_handoffs": 1,
        "target_closing": False,
    }
    return {
        "kind": "jensen_window_pf_multiplier_leading_atom_bound_certificate",
        "date": "2026-07-11",
        "status": (
            "conditional interval leading-atom bounds; not a proof of the "
            "multiplier factorization"
        ),
        "proof_boundary": (
            "Assuming the still-open unit-atomic multiplier representation, this "
            "artifact proves alpha_min>4.863538496 and at most one atom at or below "
            "11/2. It does not prove that any atom lies below 11/2, construct the "
            "counting measure, prove PF-infinity, RH, or Lambda <= 0."
        ),
        "rows": [asdict(row) for row in rows],
        "summary": summary,
        "arb_dps": ARB_DPS,
        "source_enclosures": [path.relative_to(REPO_ROOT).as_posix() for path in SOURCES],
        "source_hausdorff_bridge": (
            "outputs/jensen_window_pf_multiplier_hausdorff_uniqueness_bridge.md"
        ),
        "source_frontier": (
            "outputs/jensen_window_pf_multiplier_complete_monotonicity_frontier_scout.md"
        ),
        "source_target": "outputs/jensen_window_pf_multiplier_counting_measure_target.md",
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_multiplier_leading_atom_bound_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_multiplier_leading_atom_bound_certificate.py"
        ),
    }


def write_note(payload: dict, path: Path) -> None:
    summary = payload["summary"]
    rows = {row["id"]: row for row in payload["rows"]}
    lines = [
        "# Jensen-Window PF Multiplier Leading-Atom Bound Certificate",
        "",
        "Date: 2026-07-11",
        "",
        "Status: conditional interval leading-atom bounds; not a proof of the",
        "multiplier factorization, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_multiplier_leading_atom_bound_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_multiplier_leading_atom_bound_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_multiplier_leading_atom_bound_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_multiplier_leading_atom_bound_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        (
            "validated Jensen-window PF multiplier leading-atom bound certificate: "
            f"{summary['certificate_rows']} rows, 0 issues, "
            f"{summary['difference_orders']} difference orders, "
            f"beta6 in ({summary['root_bracket'][0]},{summary['root_bracket'][1]}), "
            f"alpha_min>{summary['conditional_alpha_min_lower_bound']}, "
            "N(11/2)<=1, 1 open existence handoff"
        ),
        "```",
        "",
        "## Conditional Kernel",
        "",
        "Under the open unit-atomic multiplier target, define",
        "",
        "```text",
        "a_m=(-1)^m*Delta^m y_1,",
        "g_m(alpha)=integral_0^1 r^(alpha-1)*(1-r)^(m+2)/(-log r)dr.",
        "```",
        "",
        "Then",
        "",
        "```text",
        "a_m=sum_j g_m(alpha_j),",
        "partial_alpha g_m(alpha)=-Beta(alpha,m+3)<0.",
        "```",
        "",
        "If `beta_m` is the unique solution of `g_m(beta_m)=a_m`, monotonicity",
        "gives `alpha_min>=beta_m`.",
        "",
        "## Interval Bound",
        "",
        "At order 6, 250-digit Arb arithmetic proves",
        "",
        "```text",
        f"g_6({ROOT_LOWER})-a_6={rows['mlab_04_root_lower_endpoint']['diagnostics']['gap_ball']}>0,",
        f"g_6({ROOT_UPPER})-a_6={rows['mlab_05_root_upper_endpoint']['diagnostics']['gap_ball']}<0.",
        "```",
        "",
        "Therefore",
        "",
        "```text",
        f"{ROOT_LOWER}<beta_6<{ROOT_UPPER},",
        f"alpha_min>{ROOT_LOWER}.",
        "```",
        "",
        "All other orders `0..55` are interval-certified weaker at the lower",
        "endpoint, so order 6 gives the strongest bound in the available triangle.",
        "",
        "## Conditional Count",
        "",
        "For any cutoff `A`, positivity and monotonicity give",
        "",
        "```text",
        "N(A)*g_m(A)<=a_m.",
        "```",
        "",
        "At `A=11/2`, order 3 gives",
        "",
        "```text",
        f"a_3/g_3(11/2)={rows['mlab_07_low_atom_count']['diagnostics']['ratio_ball']}<2,",
        "N(11/2)<=1.",
        "```",
        "",
        "This does not prove that an atom lies below `11/2`. It says that if one",
        "does, it is the unique atom in that interval. The subsequent certificate",
        "`outputs/jensen_window_pf_multiplier_unit_atomic_obstruction_certificate.md`",
        "composes the cutoff with a consecutive-moment ratio and rules out the measure",
        "entirely; the conditional bounds here are retained as its first step.",
        "",
        "```text",
        "outputs/jensen_window_pf_multiplier_hausdorff_uniqueness_bridge.md",
        "outputs/jensen_window_pf_multiplier_complete_monotonicity_frontier_scout.md",
        "outputs/jensen_window_pf_multiplier_unit_atomic_obstruction_certificate.md",
        "outputs/jensen_window_pf_multiplier_counting_measure_target.md",
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
        "validated Jensen-window PF multiplier leading-atom bound certificate: "
        f"{summary['certificate_rows']} rows, 0 issues, "
        f"{summary['difference_orders']} difference orders, "
        f"beta6 in ({summary['root_bracket'][0]},{summary['root_bracket'][1]}), "
        f"alpha_min>{summary['conditional_alpha_min_lower_bound']}, "
        "N(11/2)<=1, 1 open existence handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
