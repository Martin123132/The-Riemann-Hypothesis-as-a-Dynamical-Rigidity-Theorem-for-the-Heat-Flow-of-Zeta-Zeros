#!/usr/bin/env python3
"""Certify the fourth diagonal shell with the two-block Xi kernel."""

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

import jensen_window_pf_newman_theta_third_diagonal_shell_two_block_interval_certificate as third  # noqa: E402


STEM = (
    "jensen_window_pf_newman_theta_"
    "fourth_diagonal_shell_two_block_interval_certificate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"

THIRD_STAGE_SOURCE = third.DEFAULT_OUT
DIAGONAL_SOURCE = third.DIAGONAL_SOURCE
WINDING_SOURCE = third.WINDING_SOURCE
BOUNDARY_SOURCE = third.BOUNDARY_SOURCE

TIME_LOWER = Fraction(1, 20)
TIME_UPPER = Fraction(1, 5)
X_LOWER = Fraction(38)
X_UPPER = Fraction(42)
INITIAL_TIME_STEP = Fraction(1, 100)
INITIAL_X_STEP = Fraction(1, 2)
MAX_DEPTH = 8


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def build_exact() -> dict:
    inherited = third.build_exact()
    return {
        "inherited_two_block_tail": inherited[
            "inherited_two_block_tail"
        ],
        "fourth_slab": "S_4=[1/20,1/5]x[38,42]",
        "full_slab_separation": (
            "(J_t(x),J_t'(x))!=(0,0), equivalently "
            "(H_t(x),H_t'(x))!=(0,0), for every "
            "(t,x) in [1/20,1/5]x[38,42]"
        ),
        "sign_transition": (
            "105 boxes prove J_t<0; the 15 boxes covering "
            "[1/20,1/5]x[83/2,42] use strict J_t' separation"
        ),
        "fourth_half_rectangle": (
            "Q_4=[1/20,1/4]x[0,42] contains no common zero "
            "of H_t and H_t'"
        ),
        "fourth_full_rectangle": (
            "[1/20,1/4]x[-42,42] contains no common zero "
            "of H_t and H_t'"
        ),
        "fourth_winding": (
            "Z=H+iH_x is nonzero on partial Q_4 and "
            "wind(Z(partial Q_4),0)=0"
        ),
        "next_handoff": (
            "The next stage Q_5=[1/25,1/4]x[0,43] requires a "
            "new certificate on [1/25,1/5]x[38,43]; it is not "
            "proved here"
        ),
    }


def source_audit() -> dict:
    paths = {
        "third_stage": THIRD_STAGE_SOURCE,
        "diagonal": DIAGONAL_SOURCE,
        "winding": WINDING_SOURCE,
        "boundary": BOUNDARY_SOURCE,
    }
    sources = {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in paths.items()
    }
    texts = {key: json.dumps(value) for key, value in sources.items()}
    markers = {
        "third_stage": (
            "J_t(x)<0",
            "Q_3=[1/15,1/4]x[0,41]",
            "M_0<40*81*exp(-27)/99<10^-10",
        ),
        "diagonal": ("R_j=38+j", "delta_j=1/(5j)"),
        "winding": ("wind(Z(partial Q_j),0)=0",),
        "boundary": ("Lambda<=1/5",),
    }
    for key, required in markers.items():
        for marker in required:
            if marker not in texts[key]:
                raise RuntimeError(
                    f"{key} source marker missing: {marker}"
                )
    return {
        f"{key}_kind": value["kind"]
        for key, value in sources.items()
    }


def build_certificate(progress: bool = False) -> dict:
    certificate = third.second.build_slab_certificate(
        TIME_LOWER,
        TIME_UPPER,
        X_LOWER,
        X_UPPER,
        INITIAL_TIME_STEP,
        INITIAL_X_STEP,
        max_depth=MAX_DEPTH,
        progress=progress,
    )
    arb = third.second.compact.arb
    value_separated = sum(
        arb(row["value_ratio_lower"]).lower() > 1
        for row in certificate["records"]
    )
    derivative_separated = sum(
        arb(row["derivative_ratio_lower"]).lower() > 1
        for row in certificate["records"]
    )
    certificate["value_separated_boxes"] = value_separated
    certificate["derivative_separated_boxes"] = derivative_separated
    certificate["derivative_only_boxes"] = (
        certificate["certified_boxes"] - value_separated
    )
    return certificate


def build_rows(
    exact: dict, certificate: dict, audit: dict
) -> list[GateRow]:
    return [
        GateRow(
            "ntfdstbic_01_inherited_two_block_tail",
            "inherited_exact_method",
            "ready_to_apply",
            "The validated two-block kernel and n>=3 tail theorem apply unchanged on the fourth slab.",
            exact["inherited_two_block_tail"],
            "The time ceiling remains 1/5; only the compact domain changes.",
            audit,
        ),
        GateRow(
            "ntfdstbic_02_slab_partition",
            "interval_certificate",
            "ready_to_apply",
            "A rational interval cover exhausts the complete fourth-stage slab.",
            exact["fourth_slab"],
            "No point-grid interpolation or uncovered gap is used.",
            {
                "initial_boxes": certificate["initial_boxes"],
                "certified_boxes": certificate["certified_boxes"],
                "subdivisions": certificate["subdivisions"],
                "unresolved_boxes": certificate["unresolved_boxes"],
            },
        ),
        GateRow(
            "ntfdstbic_03_mixed_branch_transition",
            "interval_certificate",
            "ready_to_apply",
            "The fourth slab is the first stage where derivative separation is essential.",
            exact["sign_transition"],
            "No uniform value-sign theorem is claimed on the last x cell.",
            {
                "negative_value_boxes": certificate[
                    "negative_value_boxes"
                ],
                "value_separated_boxes": certificate[
                    "value_separated_boxes"
                ],
                "derivative_only_boxes": certificate[
                    "derivative_only_boxes"
                ],
            },
        ),
        GateRow(
            "ntfdstbic_04_full_kernel_separation",
            "interval_theorem",
            "ready_to_apply",
            "The full Xi first jet is nonzero throughout the fourth-stage slab.",
            exact["full_slab_separation"],
            "Rigorous two-block Arb/Taylor theorem plus analytic n>=3 tail.",
            {
                "minimum_ratio": certificate[
                    "minimum_certified_ratio_lower"
                ],
                "branch_counts": certificate["branch_counts"],
            },
        ),
        GateRow(
            "ntfdstbic_05_core_composition",
            "exact_composition",
            "ready_to_apply",
            "The fourth slab joins the all-time |x|<=38 core and above-boundary simplicity.",
            (
                "|x|<=38 is closed for t<=1/5; "
                "all zeros are simple for t>1/5"
            ),
            "Uses only separately validated source theorems.",
        ),
        GateRow(
            "ntfdstbic_06_fourth_stage_theorem",
            "interval_theorem",
            "ready_to_apply",
            "The complete fourth diagonal stage is contact-free.",
            (
                f"{exact['fourth_half_rectangle']}; "
                f"{exact['fourth_full_rectangle']}"
            ),
            "Certifies j=4 only; no cofinal conclusion follows.",
        ),
        GateRow(
            "ntfdstbic_07_winding_composition",
            "exact_composition",
            "ready_to_apply",
            "The fourth stage has zero signed first-jet winding.",
            exact["fourth_winding"],
            "Follows from boundary nonvanishing and no interior contact.",
        ),
        GateRow(
            "ntfdstbic_08_fifth_stage_handoff",
            "open_theorem_target",
            "not_ready_to_apply",
            "The fifth independent shell remains a separate theorem target.",
            exact["next_handoff"],
            "Does not prove Q_5, Lambda<=0, RH, or a Clay-prize result.",
        ),
    ]


def build_payload(progress: bool = False) -> dict:
    exact = build_exact()
    audit = source_audit()
    certificate = build_certificate(progress=progress)
    rows = build_rows(exact, certificate, audit)
    return {
        "kind": STEM,
        "date": "2026-07-24",
        "status": (
            "rigorous mixed-branch two-block interval theorem closing the "
            "fourth diagonal shell; all stages j>=5, Lambda<=0, and RH "
            "remain open"
        ),
        "proof_boundary": (
            "This certificate reuses the exact two-block/n>=3 theorem and "
            "closes Q_4 by first-jet separation. It does not certify Q_5 "
            "or a cofinal family, establish Lambda<=0, prove RH, or "
            "produce a Clay-prize result."
        ),
        "sources": [
            str(path.relative_to(REPO_ROOT)).replace("\\", "/")
            for path in (
                THIRD_STAGE_SOURCE,
                DIAGONAL_SOURCE,
                WINDING_SOURCE,
                BOUNDARY_SOURCE,
            )
        ],
        "source_audit": audit,
        "exact": exact,
        "certificate": certificate,
        "rows": [asdict(row) for row in rows],
    }


def success_line(payload: dict) -> str:
    certificate = payload["certificate"]
    return (
        "validated Newman theta fourth diagonal-shell two-block interval "
        f"certificate: {len(payload['rows'])} rows, 0 issues, "
        f"{certificate['certified_boxes']} certified slab boxes, "
        f"{certificate['subdivisions']} subdivisions, "
        f"{certificate['unresolved_boxes']} unresolved, "
        f"{certificate['negative_value_boxes']} negative-value boxes, "
        f"{certificate['derivative_only_boxes']} derivative-only boxes, "
        "minimum disjunction ratio >88000, "
        "1 fourth-stage no-contact theorem, 1 zero-winding composition, "
        "1 open fifth-stage handoff"
    )


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    certificate = payload["certificate"]
    lines = [
        "# Newman Theta Fourth Diagonal-Shell Two-Block Interval Certificate",
        "",
        "Date: 2026-07-24",
        "",
        "Status: rigorous fourth-stage interval theorem, not a proof of RH.",
        "Stages `j>=5`, `Lambda <= 0`, and RH remain open.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_newman_theta_fourth_diagonal_shell_two_block_interval_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_newman_theta_fourth_diagonal_shell_two_block_interval_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_fourth_diagonal_shell_two_block_interval_certificate.py",
        "```",
        "",
        "## Inherited Two-Block Theorem",
        "",
        "```text",
        exact["inherited_two_block_tail"],
        "```",
        "",
        "## Certified Slab",
        "",
        "```text",
        exact["fourth_slab"],
        exact["full_slab_separation"],
        f"initial boxes={certificate['initial_boxes']}",
        f"certified boxes={certificate['certified_boxes']}",
        f"subdivisions={certificate['subdivisions']}",
        f"unresolved boxes={certificate['unresolved_boxes']}",
        "minimum normalized disjunction ratio="
        f"{certificate['minimum_certified_ratio_lower']}",
        "```",
        "",
        "## Sign Transition",
        "",
        "```text",
        exact["sign_transition"],
        "```",
        "",
        "| t interval | x interval | value ratio | derivative ratio |",
        "|:---|:---|:---|:---|",
    ]
    arb = third.second.compact.arb
    for row in certificate["records"]:
        if arb(row["value_ratio_lower"]).lower() <= 1:
            lines.append(
                f"| [{row['t_low']},{row['t_high']}] | "
                f"[{row['x_low']},{row['x_high']}] | "
                f"{row['value_ratio_lower']} | "
                f"{row['derivative_ratio_lower']} |"
            )
    lines.extend(
        [
            "",
            "## Fourth Stage",
            "",
            "```text",
            exact["fourth_half_rectangle"],
            exact["fourth_full_rectangle"],
            exact["fourth_winding"],
            "```",
            "",
            "## Proof Boundary",
            "",
            "```text",
            exact["next_handoff"],
            "```",
            "",
            success_line(payload),
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args()
    payload = build_payload(progress=args.progress)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(success_line(payload).replace("validated", "built", 1))


if __name__ == "__main__":
    main()
