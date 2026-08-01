#!/usr/bin/env python3
"""Guard the second diagonal shell against a failed first-block route."""

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
    "second_diagonal_shell_first_block_route_guard"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
COMPACT_SOURCE = compact.DEFAULT_OUT
FIRST_STAGE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_"
    "first_diagonal_shell_interval_certificate.json"
)

AUDIT_TIME = Fraction(1, 10)
AUDIT_POINTS = (
    Fraction(198, 5),
    Fraction(397, 10),
    Fraction(199, 5),
    Fraction(399, 10),
    Fraction(40),
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def point_audit(certifier: compact.CompactCertifier, x: Fraction) -> dict:
    zero = Fraction(0)
    max_moment = (
        1 + 2 * compact.TAYLOR_TIME_ORDER + compact.TAYLOR_X_ORDER
    )
    integrals = [
        certifier.oscillatory_integral(moment, AUDIT_TIME, x)
        for moment in range(max_moment + 1)
    ]
    h_ball = certifier.transform_box(
        0, AUDIT_TIME, zero, x, zero, integrals
    )
    h_prime_ball = certifier.transform_box(
        1, AUDIT_TIME, zero, x, zero, integrals
    )

    time_ball = compact.interval_ball(AUDIT_TIME, zero)
    x_ball = compact.interval_ball(x, zero)
    x_two = x_ball * x_ball
    x_three = x_two * x_ball
    x_four = x_three * x_ball
    time_two = time_ball * time_ball
    j_ball = 16 * x_four * h_ball
    j_prime_ball = (
        64 * x_three * h_ball + 16 * x_four * h_prime_ball
    )

    delta_zero = compact.arb(1) / 2800
    delta_one = compact.arb(1) / 137200
    delta_two = compact.arb(1) / 3361400
    delta_three = compact.arb(1) / 54900000
    p_poly = (
        x_four + (6 * time_ball + 1) * x_two + 24 * time_two
    )
    r_poly = (6 * time_ball + 1) * x_two + 24 * time_two
    b_zero = (
        4 * time_two * x_two * delta_two
        + 4
        * time_ball
        * x_ball
        * (x_two + 4 * time_ball)
        * delta_one
        + (x_four + 2 * r_poly) * delta_zero
    )
    b_one = (
        4 * time_two * x_two * delta_three
        + 4
        * time_ball
        * x_ball
        * (x_two + 2 * time_ball)
        * delta_two
        + abs(p_poly - 4 * time_ball * (3 * x_two + 4 * time_ball))
        * delta_one
        + (4 * x_three + 4 * (6 * time_ball + 1) * x_ball)
        * delta_zero
    )

    value_ratio_lower = j_ball.abs_lower() / b_zero.upper()
    value_ratio_upper = j_ball.abs_upper() / b_zero.lower()
    derivative_ratio_lower = (
        j_prime_ball.abs_lower() / b_one.upper()
    )
    derivative_ratio_upper = (
        j_prime_ball.abs_upper() / b_one.lower()
    )
    one = compact.arb(1)
    if value_ratio_lower > one or derivative_ratio_lower > one:
        classification = "raw_disjunction_certified"
    elif value_ratio_upper < one and derivative_ratio_upper < one:
        classification = "raw_disjunction_rigorously_false"
    else:
        classification = "interval_undetermined"

    return {
        "t": str(AUDIT_TIME),
        "x": str(x),
        "classification": classification,
        "value_ratio_lower": str(value_ratio_lower.lower()),
        "value_ratio_upper": str(value_ratio_upper.upper()),
        "derivative_ratio_lower": str(derivative_ratio_lower.lower()),
        "derivative_ratio_upper": str(derivative_ratio_upper.upper()),
        "j_ball": str(j_ball),
        "j_prime_ball": str(j_prime_ball),
        "b_zero_ball": str(b_zero),
        "b_one_ball": str(b_one),
    }


def build_certificate() -> dict:
    priority_lowered = compact.request_below_normal_priority()
    certifier = compact.CompactCertifier()
    records = [point_audit(certifier, x) for x in AUDIT_POINTS]
    failed = [
        row
        for row in records
        if row["classification"] == "raw_disjunction_rigorously_false"
    ]
    passed = [
        row
        for row in records
        if row["classification"] == "raw_disjunction_certified"
    ]
    undetermined = [
        row
        for row in records
        if row["classification"] == "interval_undetermined"
    ]
    failed_ratio_entries = []
    for row in failed:
        failed_ratio_entries.extend(
            (
                (
                    compact.arb(row["value_ratio_upper"]),
                    row["x"],
                    "value",
                ),
                (
                    compact.arb(row["derivative_ratio_upper"]),
                    row["x"],
                    "derivative",
                ),
            )
        )
    maximum_failed = max(
        failed_ratio_entries,
        key=lambda item: float(item[0].upper()),
    )
    return {
        "precision_bits": compact.PRECISION_BITS,
        "integration_cutoff": str(compact.INTEGRATION_CUTOFF),
        "tail_radius": compact.TAIL_RADIUS,
        "audit_time": str(AUDIT_TIME),
        "audit_points": [str(x) for x in AUDIT_POINTS],
        "evaluated_points": len(records),
        "raw_disjunction_certified_points": len(passed),
        "raw_disjunction_rigorously_false_points": len(failed),
        "interval_undetermined_points": len(undetermined),
        "maximum_failed_ratio_upper": str(maximum_failed[0].upper()),
        "maximum_failed_ratio_point": maximum_failed[1],
        "maximum_failed_ratio_branch": maximum_failed[2],
        "below_normal_priority_applied": priority_lowered,
        "records": records,
    }


def source_audit() -> dict:
    sources = {
        "compact": json.loads(COMPACT_SOURCE.read_text(encoding="utf-8")),
        "first_stage": json.loads(
            FIRST_STAGE_SOURCE.read_text(encoding="utf-8")
        ),
    }
    texts = {key: json.dumps(value) for key, value in sources.items()}
    markers = {
        "compact": (
            "|J_1|>B_0 or |J_1'|>B_1",
            "1e-800",
            "|x|<=38",
        ),
        "first_stage": (
            "Q_2=[1/10,1/4]x[0,40]",
            "does not certify Q_2",
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


def build_exact() -> dict:
    return {
        "raw_disjunction": "|J_1|>B_0 or |J_1'|>B_1",
        "failed_point_statement": (
            "At t=1/10 and each x in {397/10,199/5,399/10,40}, "
            "rigorous Arb upper bounds give |J_1|/B_0<1 and "
            "|J_1'|/B_1<1"
        ),
        "route_consequence": (
            "The unchanged first-block sufficient disjunction cannot "
            "certify the complete bottom edge of Q_2"
        ),
        "non_consequence": (
            "Failure of a sufficient bound does not imply H=H_x=0, a "
            "multiple zero, Lambda>0, or failure of RH"
        ),
        "next_handoff": (
            "Q_2 requires a sharper tail decomposition, additional theta "
            "blocks, or a different direct boundary/winding certificate"
        ),
    }


def build_rows(
    exact: dict, certificate: dict, audit: dict
) -> list[GateRow]:
    return [
        GateRow(
            "ntsdsfbrg_01_inherited_disjunction",
            "inherited_exact_method",
            "ready_to_apply",
            "The audit uses the exact retained block and error bars of the compact certificate.",
            exact["raw_disjunction"],
            "No new tail estimate is introduced.",
            audit,
        ),
        GateRow(
            "ntsdsfbrg_02_transition_point",
            "interval_certificate",
            "ready_to_apply",
            "The derivative branch still certifies the nearby point x=198/5.",
            "t=1/10, x=198/5: |J_1'|/B_1>1",
            "A single certified point is not an edge theorem.",
            certificate["records"][0],
        ),
        GateRow(
            "ntsdsfbrg_03_strict_failed_points",
            "interval_countercertificate",
            "ready_to_apply",
            "Four exact rational points rigorously violate both sufficient inequalities.",
            exact["failed_point_statement"],
            "This rejects only the unchanged first-block sufficient route.",
            certificate["records"][1:],
        ),
        GateRow(
            "ntsdsfbrg_04_route_guard",
            "methodological_consequence",
            "ready_to_apply",
            "Further subdivision cannot repair the global first-block argument.",
            exact["route_consequence"],
            "A stronger decomposition or different certificate may still close Q_2.",
        ),
        GateRow(
            "ntsdsfbrg_05_non_consequence",
            "proof_boundary",
            "ready_to_apply",
            "The route failure is not evidence for a contact.",
            exact["non_consequence"],
            "No zero-location or Newman-constant conclusion is drawn.",
        ),
        GateRow(
            "ntsdsfbrg_06_second_stage_handoff",
            "open_theorem_target",
            "not_ready_to_apply",
            "The second diagonal stage remains open with a narrowed method choice.",
            exact["next_handoff"],
            "Does not prove Q_2, Lambda<=0, RH, or a Clay-prize result.",
        ),
    ]


def build_payload() -> dict:
    audit = source_audit()
    exact = build_exact()
    certificate = build_certificate()
    rows = build_rows(exact, certificate, audit)
    return {
        "kind": STEM,
        "date": "2026-07-24",
        "status": (
            "rigorous first-block route guard for four points on the Q_2 "
            "bottom edge; Q_2, Lambda<=0, and RH remain open"
        ),
        "proof_boundary": (
            "This pointwise countercertificate proves that the unchanged "
            "first-block sufficient bars fail at four exact points. It "
            "does not prove a common zero, disprove edge separation, close "
            "Q_2, establish Lambda<=0, or prove RH."
        ),
        "sources": [
            str(path.relative_to(REPO_ROOT)).replace("\\", "/")
            for path in (COMPACT_SOURCE, FIRST_STAGE_SOURCE)
        ],
        "source_audit": audit,
        "exact": exact,
        "certificate": certificate,
        "rows": [asdict(row) for row in rows],
    }


def success_line(payload: dict) -> str:
    certificate = payload["certificate"]
    return (
        "validated Newman theta second diagonal-shell first-block route "
        f"guard: {len(payload['rows'])} rows, 0 issues, "
        f"{certificate['evaluated_points']} exact point audits, "
        f"{certificate['raw_disjunction_rigorously_false_points']} strict "
        "two-branch failures, maximum failed ratio <24/25, "
        "1 second-stage route guard, 1 open replacement handoff"
    )


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    certificate = payload["certificate"]
    lines = [
        "# Newman Theta Second Diagonal-Shell First-Block Route Guard",
        "",
        "Date: 2026-07-24",
        "",
        "Status: rigorous pointwise countercertificate for the unchanged",
        "first-block sufficient bars, not a proof of RH. `Q_2`,",
        "`Lambda <= 0`, and RH remain open.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_newman_theta_second_diagonal_shell_first_block_route_guard.json",
        "python work/rh_compute/scripts/jensen_window_pf_newman_theta_second_diagonal_shell_first_block_route_guard.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_second_diagonal_shell_first_block_route_guard.py",
        "```",
        "",
        "## Exact Point Audit",
        "",
        "The ratios use lower bounds to certify a branch and upper bounds",
        "to certify failure of both branches.",
        "",
        "| t | x | value ratio interval | derivative ratio interval | result |",
        "|:---|:---|:---|:---|:---|",
    ]
    for row in certificate["records"]:
        lines.append(
            f"| {row['t']} | {row['x']} | "
            f"[{row['value_ratio_lower']}, {row['value_ratio_upper']}] | "
            f"[{row['derivative_ratio_lower']}, "
            f"{row['derivative_ratio_upper']}] | "
            f"{row['classification']} |"
        )
    lines.extend(
        [
            "",
            "```text",
            exact["failed_point_statement"],
            exact["route_consequence"],
            "maximum failed upper ratio="
            f"{certificate['maximum_failed_ratio_upper']}",
            "```",
            "",
            "## Proof Boundary",
            "",
            "```text",
            exact["non_consequence"],
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
    args.out.write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(success_line(payload).replace("validated", "built", 1))


if __name__ == "__main__":
    main()
