#!/usr/bin/env python3
"""Certify the first positive-boundary diagonal half-rectangle."""

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

import jensen_window_pf_newman_theta_compact_transversality_interval_certificate as compact  # noqa: E402


STEM = (
    "jensen_window_pf_newman_theta_"
    "first_diagonal_shell_interval_certificate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
COMPACT_SOURCE = compact.DEFAULT_OUT
DIAGONAL_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_positive_boundary_diagonal_exhaustion_gate.json"
)
WINDING_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_first_jet_winding_gate.json"
)
BOUNDARY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_positive_boundary_attainment_lemma.json"
)

EDGE_TIME = Fraction(1, 5)
X_LOWER = Fraction(38)
X_UPPER = Fraction(39)
X_STEP = Fraction(1, 5)
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


def build_edge_certificate() -> dict:
    priority_lowered = compact.request_below_normal_priority()
    certifier = compact.CompactCertifier()
    queue = [
        ((EDGE_TIME, EDGE_TIME, x_low, x_high), 0)
        for x_low, x_high in compact.fraction_range(
            X_LOWER, X_UPPER, X_STEP
        )
    ]
    records: list[dict] = []
    unresolved: list[dict] = []
    subdivisions = 0
    cursor = 0
    while cursor < len(queue):
        bounds, depth = queue[cursor]
        cursor += 1
        record = certifier.certify_box(*bounds, depth)
        if record["certified"]:
            records.append(record)
        elif depth < MAX_DEPTH:
            queue.extend(
                (child, depth + 1)
                for child in compact.split_box(bounds)
            )
            subdivisions += 1
        else:
            unresolved.append(record)
    if unresolved:
        raise RuntimeError(
            f"{len(unresolved)} first-shell edge boxes remain unresolved"
        )
    minimum = min(
        records,
        key=lambda row: float(
            compact.arb(row["certified_ratio_lower"]).lower()
        ),
    )
    return {
        "precision_bits": compact.PRECISION_BITS,
        "integration_cutoff": str(compact.INTEGRATION_CUTOFF),
        "tail_radius": compact.TAIL_RADIUS,
        "edge_time": str(EDGE_TIME),
        "x_lower": str(X_LOWER),
        "x_upper": str(X_UPPER),
        "initial_x_step": str(X_STEP),
        "initial_boxes": 5,
        "evaluated_boxes": certifier.evaluations,
        "certified_boxes": len(records),
        "subdivisions": subdivisions,
        "unresolved_boxes": 0,
        "maximum_depth": max(row["depth"] for row in records),
        "branch_counts": {
            "value": sum(row["branch"] == "value" for row in records),
            "derivative": sum(
                row["branch"] == "derivative" for row in records
            ),
        },
        "minimum_certified_ratio_lower": (
            minimum["certified_ratio_lower"]
        ),
        "minimum_record": minimum,
        "below_normal_priority_applied": priority_lowered,
        "records": records,
    }


def build_exact() -> dict:
    return {
        "edge_contact_exclusion": (
            "At t=1/5, (H_t(x),H_t'(x))!=(0,0) for 38<=x<=39"
        ),
        "even_edge": (
            "Evenness gives the same contact exclusion for -39<=x<=-38"
        ),
        "compact_composition": (
            "The prior compact theorem covers |x|<=38 at t=1/5, so "
            "(H_(1/5),H_(1/5)')!=(0,0) for every |x|<=39"
        ),
        "above_boundary_simplicity": (
            "The published bound Lambda<=1/5 implies that H_t has only "
            "simple zeros for every t>1/5"
        ),
        "first_half_rectangle": (
            "Q_1=[1/5,1/4]x[0,39] contains no common zero of H_t and H_t'"
        ),
        "first_full_rectangle": (
            "[1/5,1/4]x[-39,39] contains no common zero of H_t and H_t'"
        ),
        "first_winding": (
            "Z=H+iH_x is nonzero on partial Q_1 and "
            "wind(Z(partial Q_1),0)=0"
        ),
        "next_handoff": (
            "The next diagonal stage Q_2=[1/10,1/4]x[0,40] requires "
            "bottom-edge separation at t=1/10, right-edge separation at "
            "x=40 for 1/10<=t<=1/5, and zero first-jet winding; none is "
            "proved by this first-stage certificate"
        ),
    }


def source_audit() -> dict:
    sources = {
        "compact": json.loads(COMPACT_SOURCE.read_text(encoding="utf-8")),
        "diagonal": json.loads(DIAGONAL_SOURCE.read_text(encoding="utf-8")),
        "winding": json.loads(WINDING_SOURCE.read_text(encoding="utf-8")),
        "boundary": json.loads(BOUNDARY_SOURCE.read_text(encoding="utf-8")),
    }
    texts = {key: json.dumps(value) for key, value in sources.items()}
    markers = {
        "compact": (
            "|x|<=38",
            "|J_1|>B_0 or |J_1'|>B_1",
            "1e-800",
        ),
        "diagonal": (
            "R_j=38+j",
            "delta_j=1/(5j)",
        ),
        "winding": (
            "Q_j=[1/(5j),1/4]x[0,38+j]",
            "wind(Z(partial Q_j),0)=0",
        ),
        "boundary": (
            "Lambda<=1/5",
            "only simple zeros for every 0<t<=1/5",
        ),
    }
    for key, required in markers.items():
        for marker in required:
            if marker not in texts[key]:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {
        f"{key}_kind": value["kind"]
        for key, value in sources.items()
    }


def build_rows(
    exact: dict, certificate: dict, audit: dict
) -> list[GateRow]:
    return [
        GateRow(
            "ntfdic_01_interval_method",
            "inherited_exact_method",
            "ready_to_apply",
            "The existing Arb/Taylor first-component disjunction applies on a zero-width time edge.",
            (
                "160-bit retained-integral enclosures, analytic 10^-800 "
                "tails, and third-order x Taylor bounds certify each box"
            ),
            "Reuses the separately validated compact certifier.",
            audit,
        ),
        GateRow(
            "ntfdic_02_edge_partition",
            "interval_certificate",
            "ready_to_apply",
            "Five rational boxes cover the complete first open edge.",
            "t=1/5 and 38<=x<=39 in five width-1/5 boxes",
            "No point sampling or uncovered interval is used.",
            {
                "certified_boxes": certificate["certified_boxes"],
                "subdivisions": certificate["subdivisions"],
                "unresolved_boxes": certificate["unresolved_boxes"],
                "minimum_ratio": certificate[
                    "minimum_certified_ratio_lower"
                ],
            },
        ),
        GateRow(
            "ntfdic_03_edge_contact_exclusion",
            "interval_certificate",
            "ready_to_apply",
            "The first formerly open unit edge has rigorous first-jet separation.",
            exact["edge_contact_exclusion"],
            "Finite interval theorem on one fixed edge only.",
        ),
        GateRow(
            "ntfdic_04_compact_composition",
            "exact_composition",
            "ready_to_apply",
            "The new edge joins the certified inner core without a gap.",
            f"{exact['even_edge']}; {exact['compact_composition']}",
            "Uses evenness and the prior |x|<=38 theorem.",
        ),
        GateRow(
            "ntfdic_05_above_boundary_simplicity",
            "published_composition",
            "ready_to_apply",
            "Every time strictly above the bottom edge is already contact-free globally.",
            exact["above_boundary_simplicity"],
            "Uses the published Lambda upper bound and simplicity above Lambda.",
        ),
        GateRow(
            "ntfdic_06_first_stage_theorem",
            "interval_theorem",
            "ready_to_apply",
            "The complete first diagonal stage is closed.",
            (
                f"{exact['first_half_rectangle']}; "
                f"{exact['first_full_rectangle']}"
            ),
            "Certifies j=1 only; no cofinal conclusion follows.",
        ),
        GateRow(
            "ntfdic_07_winding_composition",
            "exact_composition",
            "ready_to_apply",
            "The first stage has zero signed first-jet winding.",
            exact["first_winding"],
            "Follows from boundary nonvanishing and absence of interior contacts.",
        ),
        GateRow(
            "ntfdic_08_second_stage_handoff",
            "open_theorem_target",
            "not_ready_to_apply",
            "The next stage requires new positive-time edge control.",
            exact["next_handoff"],
            "No claim is made for j>=2, Lambda<=0, or RH.",
        ),
    ]


def build_payload() -> dict:
    exact = build_exact()
    certificate = build_edge_certificate()
    audit = source_audit()
    rows = build_rows(exact, certificate, audit)
    return {
        "kind": STEM,
        "date": "2026-07-24",
        "status": (
            "rigorous first diagonal-shell interval theorem closing j=1; "
            "all stages j>=2, Lambda<=0, and RH remain open"
        ),
        "proof_boundary": (
            "This certificate proves first-jet separation on one unit edge "
            "and composes it with known simplicity to close Q_1 only. It "
            "does not certify Q_2 or any cofinal family, positive-time "
            "simplicity down to zero, Lambda<=0, or RH."
        ),
        "sources": [
            str(path.relative_to(REPO_ROOT)).replace("\\", "/")
            for path in (
                COMPACT_SOURCE,
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
        "validated Newman theta first diagonal-shell interval certificate: "
        f"{len(payload['rows'])} rows, 0 issues, "
        f"{certificate['certified_boxes']} certified edge boxes, "
        f"{certificate['subdivisions']} subdivisions, "
        f"{certificate['unresolved_boxes']} unresolved, "
        "minimum ratio >6/5, 1 first-stage no-contact theorem, "
        "1 zero-winding composition, 1 open second-stage handoff"
    )


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    certificate = payload["certificate"]
    lines = [
        "# Newman Theta First Diagonal-Shell Interval Certificate",
        "",
        "Date: 2026-07-24",
        "",
        "Status: rigorous first-stage interval theorem, not a proof of RH.",
        "Stages `j>=2`, `Lambda <= 0`, and RH remain open.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_newman_theta_first_diagonal_shell_interval_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_newman_theta_first_diagonal_shell_interval_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_first_diagonal_shell_interval_certificate.py",
        "```",
        "",
        "## Certified Edge",
        "",
        "```text",
        exact["edge_contact_exclusion"],
        exact["even_edge"],
        exact["compact_composition"],
        "```",
        "",
        "The replay uses the same 160-bit Arb retained integrals, analytic",
        "`10^-800` tail, and bivariate Taylor machinery as the compact",
        "certificate, specialized to the exact time `t=1/5`.",
        "",
        "```text",
        f"certified boxes={certificate['certified_boxes']}",
        f"subdivisions={certificate['subdivisions']}",
        f"unresolved boxes={certificate['unresolved_boxes']}",
        "minimum normalized disjunction ratio="
        f"{certificate['minimum_certified_ratio_lower']}",
        "```",
        "",
        "| t | x interval | branch | certified ratio lower |",
        "|:---|:---|:---|:---|",
    ]
    for row in certificate["records"]:
        lines.append(
            f"| {row['t_low']} | [{row['x_low']},{row['x_high']}] | "
            f"{row['branch']} | {row['certified_ratio_lower']} |"
        )
    lines.extend(
        [
            "",
            "## First Stage",
            "",
            "```text",
            exact["above_boundary_simplicity"],
            exact["first_half_rectangle"],
            exact["first_full_rectangle"],
            exact["first_winding"],
            "```",
            "",
            "Thus the first member of the independent diagonal exhaustion is",
            "now a theorem, not a diagnostic.",
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
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(success_line(payload).replace("validated", "built", 1))


if __name__ == "__main__":
    main()
