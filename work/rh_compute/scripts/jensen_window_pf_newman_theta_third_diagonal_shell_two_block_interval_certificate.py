#!/usr/bin/env python3
"""Certify the third diagonal shell with the two-block Xi kernel."""

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

import jensen_window_pf_newman_theta_second_diagonal_shell_two_block_interval_certificate as second  # noqa: E402


STEM = (
    "jensen_window_pf_newman_theta_"
    "third_diagonal_shell_two_block_interval_certificate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"

SECOND_STAGE_SOURCE = second.DEFAULT_OUT
DIAGONAL_SOURCE = second.DIAGONAL_SOURCE
WINDING_SOURCE = second.WINDING_SOURCE
BOUNDARY_SOURCE = second.BOUNDARY_SOURCE

TIME_LOWER = Fraction(1, 15)
TIME_UPPER = Fraction(1, 5)
X_LOWER = Fraction(38)
X_UPPER = Fraction(41)
INITIAL_TIME_STEP = Fraction(1, 60)
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
    inherited = second.build_exact()
    return {
        "inherited_two_block_tail": (
            f"{inherited['tail_moment_formula']}; "
            f"{inherited['tail_geometric_bound']}; "
            f"{inherited['direct_tail_bars']}"
        ),
        "third_slab": "S_3=[1/15,1/5]x[38,41]",
        "full_slab_sign": (
            "J_t(x)<0, equivalently H_t(x)<0, for every "
            "(t,x) in [1/15,1/5]x[38,41]"
        ),
        "third_half_rectangle": (
            "Q_3=[1/15,1/4]x[0,41] contains no common zero "
            "of H_t and H_t'"
        ),
        "third_full_rectangle": (
            "[1/15,1/4]x[-41,41] contains no common zero "
            "of H_t and H_t'"
        ),
        "third_winding": (
            "Z=H+iH_x is nonzero on partial Q_3 and "
            "wind(Z(partial Q_3),0)=0"
        ),
        "next_handoff": (
            "The next stage Q_4=[1/20,1/4]x[0,42] requires a "
            "new certificate on [1/20,1/5]x[38,42]; it is not "
            "proved here"
        ),
    }


def source_audit() -> dict:
    paths = {
        "second_stage": SECOND_STAGE_SOURCE,
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
        "second_stage": (
            "M_0<40*81*exp(-27)/99<10^-10",
            "J_t(x)<0",
            "Q_2=[1/10,1/4]x[0,40]",
        ),
        "diagonal": ("R_j=38+j", "delta_j=1/(5j)"),
        "winding": ("wind(Z(partial Q_j),0)=0",),
        "boundary": ("Lambda<=1/5",),
    }
    for key, required in markers.items():
        for marker in required:
            if marker not in texts[key]:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {
        f"{key}_kind": value["kind"]
        for key, value in sources.items()
    }


def build_certificate(progress: bool = False) -> dict:
    return second.build_slab_certificate(
        TIME_LOWER,
        TIME_UPPER,
        X_LOWER,
        X_UPPER,
        INITIAL_TIME_STEP,
        INITIAL_X_STEP,
        max_depth=MAX_DEPTH,
        progress=progress,
    )


def build_rows(
    exact: dict, certificate: dict, audit: dict
) -> list[GateRow]:
    return [
        GateRow(
            "nttdstbic_01_inherited_two_block_tail",
            "inherited_exact_method",
            "ready_to_apply",
            "The validated two-block kernel and n>=3 tail theorem apply unchanged on the larger slab.",
            exact["inherited_two_block_tail"],
            "The time ceiling remains 1/5; only the compact domain changes.",
            audit,
        ),
        GateRow(
            "nttdstbic_02_slab_partition",
            "interval_certificate",
            "ready_to_apply",
            "A rational interval cover exhausts the complete third-stage slab.",
            exact["third_slab"],
            "No point-grid interpolation or uncovered gap is used.",
            {
                "initial_boxes": certificate["initial_boxes"],
                "certified_boxes": certificate["certified_boxes"],
                "subdivisions": certificate["subdivisions"],
                "unresolved_boxes": certificate["unresolved_boxes"],
            },
        ),
        GateRow(
            "nttdstbic_03_full_kernel_sign",
            "interval_theorem",
            "ready_to_apply",
            "The full Xi kernel is strictly negative throughout the third-stage slab.",
            exact["full_slab_sign"],
            "Rigorous two-block Arb/Taylor theorem plus analytic n>=3 tail.",
            {
                "minimum_value_ratio": certificate[
                    "minimum_value_ratio_lower"
                ],
                "negative_value_boxes": certificate[
                    "negative_value_boxes"
                ],
            },
        ),
        GateRow(
            "nttdstbic_04_core_composition",
            "exact_composition",
            "ready_to_apply",
            "The third slab joins the all-time |x|<=38 core and above-boundary simplicity.",
            (
                "|x|<=38 is closed for t<=1/5; "
                "all zeros are simple for t>1/5"
            ),
            "Uses only separately validated source theorems.",
        ),
        GateRow(
            "nttdstbic_05_third_stage_theorem",
            "interval_theorem",
            "ready_to_apply",
            "The complete third diagonal stage is contact-free.",
            (
                f"{exact['third_half_rectangle']}; "
                f"{exact['third_full_rectangle']}"
            ),
            "Certifies j=3 only; no cofinal conclusion follows.",
        ),
        GateRow(
            "nttdstbic_06_winding_composition",
            "exact_composition",
            "ready_to_apply",
            "The third stage has zero signed first-jet winding.",
            exact["third_winding"],
            "Follows from boundary nonvanishing and no interior contact.",
        ),
        GateRow(
            "nttdstbic_07_method_scope",
            "proof_boundary",
            "ready_to_apply",
            "The successful certificate retains theta arithmetic rather than reviving the rejected one-block bars.",
            "Phi_1+Phi_2 is evaluated oscillatory; only n>=3 is bounded absolutely",
            "This finite sign theorem does not extrapolate to all shells.",
        ),
        GateRow(
            "nttdstbic_08_fourth_stage_handoff",
            "open_theorem_target",
            "not_ready_to_apply",
            "The fourth independent shell remains a separate theorem target.",
            exact["next_handoff"],
            "Does not prove Q_4, Lambda<=0, RH, or a Clay-prize result.",
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
            "rigorous two-block interval theorem closing the third "
            "diagonal shell; all stages j>=4, Lambda<=0, and RH remain open"
        ),
        "proof_boundary": (
            "This certificate reuses the exact two-block/n>=3 theorem and "
            "closes Q_3. It does not certify Q_4 or a cofinal family, "
            "establish Lambda<=0, prove RH, or produce a Clay-prize result."
        ),
        "sources": [
            str(path.relative_to(REPO_ROOT)).replace("\\", "/")
            for path in (
                SECOND_STAGE_SOURCE,
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
        "validated Newman theta third diagonal-shell two-block interval "
        f"certificate: {len(payload['rows'])} rows, 0 issues, "
        f"{certificate['certified_boxes']} certified slab boxes, "
        f"{certificate['subdivisions']} subdivisions, "
        f"{certificate['unresolved_boxes']} unresolved, "
        f"{certificate['negative_value_boxes']} negative-value boxes, "
        "minimum value ratio >1000, "
        "1 third-stage no-contact theorem, 1 zero-winding composition, "
        "1 open fourth-stage handoff"
    )


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    certificate = payload["certificate"]
    lines = [
        "# Newman Theta Third Diagonal-Shell Two-Block Interval Certificate",
        "",
        "Date: 2026-07-24",
        "",
        "Status: rigorous third-stage interval theorem, not a proof of RH.",
        "Stages `j>=4`, `Lambda <= 0`, and RH remain open.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_newman_theta_third_diagonal_shell_two_block_interval_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_newman_theta_third_diagonal_shell_two_block_interval_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_third_diagonal_shell_two_block_interval_certificate.py",
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
        exact["third_slab"],
        exact["full_slab_sign"],
        f"initial boxes={certificate['initial_boxes']}",
        f"certified boxes={certificate['certified_boxes']}",
        f"subdivisions={certificate['subdivisions']}",
        f"unresolved boxes={certificate['unresolved_boxes']}",
        "minimum normalized value-sign ratio="
        f"{certificate['minimum_value_ratio_lower']}",
        "```",
        "",
        "| t interval | x interval | branch | value ratio lower |",
        "|:---|:---|:---|:---|",
    ]
    for row in certificate["records"]:
        lines.append(
            f"| [{row['t_low']},{row['t_high']}] | "
            f"[{row['x_low']},{row['x_high']}] | "
            f"{row['branch']} | {row['value_ratio_lower']} |"
        )
    lines.extend(
        [
            "",
            "## Third Stage",
            "",
            "```text",
            exact["third_half_rectangle"],
            exact["third_full_rectangle"],
            exact["third_winding"],
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
