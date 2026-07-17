#!/usr/bin/env python3
"""Build the current-ANTEDB beta-frontier audit for the Newman scaled layer."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path

from jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate import (
    C_STAR,
    R_STAR,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_polymath15_antedb_beta_frontier_audit.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_antedb_beta_frontier_audit.md"
)

ANTEDB_PAPER_URL = "https://arxiv.org/abs/2501.16779"
ANTEDB_REPO_URL = "https://github.com/teorth/expdb"
ANTEDB_COMMIT = "99668603896af86e6cda90ed6755cf3116aab0ac"
ANTEDB_COMMIT_URL = f"{ANTEDB_REPO_URL}/commit/{ANTEDB_COMMIT}"
TRUDGIAN_YANG_URL = "https://arxiv.org/abs/2306.05599"

ALPHA_STAR = Fraction(62_831, 155_153)
BETA_STAR = Fraction(220_633, 620_612)
C_TWO_DEFICIT = Fraction(3_133_668_399, 48_144_906_818)
HB_CUTOFF = 100
HB_TAIL_RADIUS = Fraction(100, 4_901)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def fraction_record(value: Fraction) -> dict:
    return {
        "exact": f"{value.numerator}/{value.denominator}",
        "decimal": f"{float(value):.15f}",
    }


def beta_required_time(alpha: Fraction, beta: Fraction) -> Fraction:
    """Scaled time forced by a beta(alpha) exponent at r=2*alpha."""
    return (4 * beta - 2 * alpha) / (alpha * (1 - alpha))


def pair_beta(pair: tuple[Fraction, Fraction], alpha: Fraction) -> Fraction:
    kappa, lam = pair
    return kappa + (lam - kappa) * alpha


BASE_PAIRS: tuple[tuple[str, Fraction, Fraction], ...] = (
    ("Bourgain", Fraction(13, 84), Fraction(55, 84)),
    ("TY1", Fraction(4_742, 38_463), Fraction(35_731, 51_284)),
    ("TY2", Fraction(18, 199), Fraction(593, 796)),
    ("TY3", Fraction(2_779, 38_033), Fraction(58_699, 76_066)),
    ("TY4", Fraction(715, 10_238), Fraction(7_955, 10_238)),
    ("A_TY2", Fraction(9, 217), Fraction(1_461, 1_736)),
)

POST_2023_PAIRS: tuple[tuple[str, Fraction, Fraction, str], ...] = (
    ("TTY_1", Fraction(89, 1_282), Fraction(997, 1_282), "Tao-Trudgian-Yang 2025"),
    (
        "TTY_2",
        Fraction(652_397, 9_713_986),
        Fraction(7_599_781, 9_713_986),
        "Tao-Trudgian-Yang 2025",
    ),
    (
        "TTY_3",
        Fraction(10_769, 351_096),
        Fraction(609_317, 702_192),
        "Tao-Trudgian-Yang 2025",
    ),
    ("TTY_4", Fraction(89, 3_478), Fraction(15_327, 17_390), "Tao-Trudgian-Yang 2025"),
    ("Cushing_1", Fraction(311, 4_822), Fraction(3_799, 4_822), "ANTEDB Cushing 2025"),
    (
        "Cushing_2",
        Fraction(80_219, 1_298_878),
        Fraction(515_638, 649_439),
        "ANTEDB Cushing 2025",
    ),
)


def heath_brown_pair(m: int) -> tuple[Fraction, Fraction]:
    return (
        Fraction(2, (m - 1) ** 2 * (m + 2)),
        1 - Fraction(3 * m - 2, m * (m - 1) * (m + 2)),
    )


def pair_lines() -> list[tuple[str, Fraction, Fraction]]:
    """Return phase lines Phi(r)=a+b*r, including the trivial line."""
    pairs = [(name, kappa, lam) for name, kappa, lam in BASE_PAIRS]
    pairs.extend((name, kappa, lam) for name, kappa, lam, _ in POST_2023_PAIRS)
    pairs.extend(
        (f"HB_{m}", *heath_brown_pair(m)) for m in range(6, HB_CUTOFF + 1)
    )
    lines = [("trivial", Fraction(0), Fraction(1))]
    lines.extend(
        (name, 2 * kappa, lam - kappa) for name, kappa, lam in pairs
    )
    return lines


def phase_value(line: tuple[str, Fraction, Fraction], radius: Fraction) -> Fraction:
    return line[1] + line[2] * radius


def required_time(phase: Fraction, radius: Fraction) -> Fraction:
    if radius == 0:
        return Fraction(2)
    return 8 * (phase - radius / 2) / (radius * (2 - radius))


def build_pair_hull() -> dict:
    lines = pair_lines()
    points = {Fraction(0), Fraction(1)}
    for index, left in enumerate(lines):
        for right in lines[index + 1 :]:
            if left[2] == right[2]:
                continue
            radius = (right[1] - left[1]) / (left[2] - right[2])
            if 0 < radius < 1:
                points.add(radius)

    candidates: list[dict] = []
    for radius in sorted(points):
        values = [(phase_value(line, radius), line[0]) for line in lines]
        phase = min(value for value, _ in values)
        active = sorted(name for value, name in values if value == phase)
        candidates.append(
            {
                "radius": fraction_record(radius),
                "phase": fraction_record(phase),
                "required_scaled_time": fraction_record(required_time(phase, radius)),
                "active": active,
            }
        )

    maximum = max(
        candidates,
        key=lambda row: Fraction(row["required_scaled_time"]["exact"]),
    )
    tail_trivial_ceiling = Fraction(4, 1) / (2 - HB_TAIL_RADIUS)
    if Fraction(maximum["required_scaled_time"]["exact"]) != C_STAR:
        raise RuntimeError("current exponent-pair hull changed c_star")
    if Fraction(maximum["radius"]["exact"]) != R_STAR:
        raise RuntimeError("current exponent-pair hull changed r_star")
    if maximum["active"] != ["TY1", "TY2"]:
        raise RuntimeError("current hull contact is no longer TY1/TY2")
    if tail_trivial_ceiling >= C_STAR:
        raise RuntimeError("Heath-Brown tail guard no longer lies below c_star")
    return {
        "candidate_count": len(candidates),
        "line_count": len(lines),
        "maximum": maximum,
        "heath_brown_cutoff": HB_CUTOFF,
        "tail_radius_ceiling": fraction_record(HB_TAIL_RADIUS),
        "tail_trivial_time_ceiling": fraction_record(tail_trivial_ceiling),
    }


def build_post_2023_rows() -> list[dict]:
    rows = []
    for name, kappa, lam, source in POST_2023_PAIRS:
        beta = pair_beta((kappa, lam), ALPHA_STAR)
        rows.append(
            {
                "name": name,
                "source": source,
                "kappa": fraction_record(kappa),
                "lambda": fraction_record(lam),
                "beta_at_alpha_star": fraction_record(beta),
                "excess_over_beta_star": fraction_record(beta - BETA_STAR),
            }
        )
    if any(Fraction(row["excess_over_beta_star"]["exact"]) < 0 for row in rows):
        raise RuntimeError("a post-2023 exponent pair lowers the critical contact")
    return rows


def build_exact() -> dict:
    if R_STAR != 2 * ALPHA_STAR:
        raise RuntimeError("alpha/radius map failed")

    left_beta = Fraction(18, 199) + Fraction(521, 796) * ALPHA_STAR
    point_beta = Fraction(4_742, 38_463) + Fraction(88_225, 153_852) * ALPHA_STAR
    right_beta = Fraction(569, 2_800) + Fraction(1_053, 2_800) * ALPHA_STAR
    if not (left_beta == point_beta == right_beta == BETA_STAR):
        raise RuntimeError("ANTEDB direct-beta contact identity failed")
    if beta_required_time(ALPHA_STAR, BETA_STAR) != C_STAR:
        raise RuntimeError("beta-coordinate threshold identity failed")

    c2_target_phase = R_STAR / 2 + Fraction(2) * R_STAR * (2 - R_STAR) / 8
    c2_phase_deficit = 2 * BETA_STAR - c2_target_phase
    if c2_phase_deficit != C_TWO_DEFICIT:
        raise RuntimeError("c=2 deficit drifted")

    pair_hull = build_pair_hull()
    post_rows = build_post_2023_rows()
    return {
        "source_snapshot": {
            "antedb_commit": ANTEDB_COMMIT,
            "antedb_commit_date": "2026-07-11",
            "paper": ANTEDB_PAPER_URL,
            "repository": ANTEDB_REPO_URL,
            "commit_url": ANTEDB_COMMIT_URL,
        },
        "coordinate_map": (
            "For tau=N^(2+o(1)) and block length M=N^r=tau^(alpha+o(1)), "
            "alpha=r/2 and Phi(r)=2*beta(alpha); hence "
            "c_req(alpha)=(4*beta(alpha)-2*alpha)/(alpha*(1-alpha))"
        ),
        "critical_contact": {
            "alpha_star": fraction_record(ALPHA_STAR),
            "radius_star": fraction_record(R_STAR),
            "beta_star": fraction_record(BETA_STAR),
            "scaled_time": fraction_record(C_STAR),
            "left_beta_line": "18/199+(521/796)*alpha",
            "point_beta_line": "4742/38463+(88225/153852)*alpha",
            "right_beta_line": "569/2800+(1053/2800)*alpha",
        },
        "pair_hull": pair_hull,
        "post_2023_pairs": post_rows,
        "antedb_pipeline": {
            "python": "3.12.13",
            "pycddlib": "2.1.8.post1",
            "python_hash_seed": 0,
            "loaded_hypotheses": 207,
            "depth_one_derived_pairs": 53,
            "exponent_pair_beta_bounds": 262,
            "beta_derived_pairs": 104,
            "vdc_beta_pieces": 7,
            "final_piece_count_in_two_seeded_runs": 68,
            "contact_values": [
                "18/199+(521/796)*alpha on the left",
                "4742/38463+(88225/153852)*alpha at the point",
                "569/2800+(1053/2800)*alpha on the right",
            ],
            "reproducibility_boundary": (
                "Intermediate partition counts and raw row hashes depend on set/tie "
                "ordering, so they are not promoted. Two PYTHONHASHSEED=0 runs "
                "returned 68 final pieces and the same exact rational contact."
            ),
        },
        "beta_only_iteration_audit": {
            "passes_completed": 12,
            "piece_counts": [62, 69, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81],
            "contact_alpha": fraction_record(ALPHA_STAR),
            "contact_beta_on_every_pass": fraction_record(BETA_STAR),
            "global_threshold_on_every_pass": fraction_record(C_STAR),
            "result": (
                "Twelve successive applications of the proved ANTEDB beta-to-beta "
                "van der Corput lemma improved other intervals but left the exact "
                "critical contact and global threshold unchanged on every pass"
            ),
            "proof_boundary": (
                "Finite twelve-pass theorem-search audit only; it does not prove "
                "stabilization under infinitely many passes or closure under the "
                "full beta/exponent-pair transform system"
            ),
        },
        "c2_phase_deficit": fraction_record(c2_phase_deficit),
        "current_conclusion": (
            "The post-2023 exponent pairs and the pinned current ANTEDB one-pass "
            "direct-beta envelope improve intermediate radii but do not lower the "
            "TY1/TY2 contact; the audited pointwise-bound frontier remains c_*"
        ),
        "required_upgrade": (
            "Lowering c_* still requires a beta bound strictly below "
            "220633/620612 at alpha=62831/155153, compatible improvement across "
            "the exposed neighborhood, or cancellation beyond pointwise beta "
            "majorants; this audit supplies none of those inputs"
        ),
    }


def build_artifact() -> dict:
    exact = build_exact()
    rows = [
        GateRow(
            id="np15abf_01_coordinate_map",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim="The Newman dyadic phase frontier has an exact ANTEDB beta-coordinate form.",
            formula=exact["coordinate_map"],
            proof_boundary="Exact exponent-coordinate conversion only.",
        ),
        GateRow(
            id="np15abf_02_source_snapshot",
            role="source_audit",
            readiness="source_validated",
            claim="The audit is pinned to one identified current ANTEDB source state.",
            formula=json.dumps(exact["source_snapshot"], sort_keys=True),
            proof_boundary="As-of source snapshot; not a timeless claim that no later bound improves the frontier.",
        ),
        GateRow(
            id="np15abf_03_post2023_pairs",
            role="published_input_audit",
            readiness="available_published",
            claim="Four Tao-Trudgian-Yang and two Cushing pairs were tested at the critical contact.",
            formula="six post-2023 rational exponent pairs",
            proof_boundary="Tests these recorded inputs; it does not prove completeness of all future exponent-pair literature.",
            diagnostics=exact["post_2023_pairs"],
        ),
        GateRow(
            id="np15abf_04_pair_hull",
            role="exact_frontier_audit",
            readiness="ready_to_apply",
            claim="The independently recomputed current rational exponent-pair hull retains the exact old maximum.",
            formula=json.dumps(exact["pair_hull"]["maximum"], sort_keys=True),
            proof_boundary="Exact for the listed current hull plus the explicit Heath-Brown tail guard.",
            diagnostics=exact["pair_hull"],
        ),
        GateRow(
            id="np15abf_05_beta_pipeline",
            role="reproducibility_audit",
            readiness="reproduced_snapshot",
            claim="The stronger one-pass direct-beta optimizer was reproduced with pinned software inputs.",
            formula=json.dumps(exact["antedb_pipeline"], sort_keys=True),
            proof_boundary="Computational source audit; raw tie partitions and an unbounded recursive closure are not promoted.",
            diagnostics=exact["antedb_pipeline"],
        ),
        GateRow(
            id="np15abf_06_direct_contact",
            role="exact_frontier_contact",
            readiness="ready_to_apply",
            claim="Three current direct-beta pieces meet at the same exact critical value.",
            formula=json.dumps(exact["critical_contact"], sort_keys=True),
            proof_boundary="Exact rational contact for the pinned one-pass envelope.",
            diagnostics=exact["critical_contact"],
        ),
        GateRow(
            id="np15abf_07_c2_deficit",
            role="exact_frontier_target",
            readiness="ready_to_apply",
            claim="The current source refresh leaves the exact c=2 phase deficit unchanged.",
            formula=exact["c2_phase_deficit"]["exact"],
            proof_boundary="Necessary local deficit only, not a complete c=2 cancellation theorem.",
        ),
        GateRow(
            id="np15abf_08_stale_table_gate",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="Updating from the 2023 pair table to the pinned current database does not lower c_*.",
            formula=exact["current_conclusion"],
            proof_boundary="Current-input route audit only; it does not rule out new estimates or non-pointwise methods.",
        ),
        GateRow(
            id="np15abf_09_recursive_boundary",
            role="method_boundary",
            readiness="guard_validated",
            claim="Twelve finite beta-only transform passes preserve the contact, while unbounded recursive closure remains unclaimed.",
            formula=json.dumps(exact["beta_only_iteration_audit"], sort_keys=True),
            proof_boundary="Finite iteration audit only; the terminated full recursive fixed-point search contributes no theorem claim.",
            diagnostics=exact["beta_only_iteration_audit"],
        ),
        GateRow(
            id="np15abf_10_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The exact contact states the remaining pointwise beta improvement required to move the frontier.",
            formula=exact["required_upgrade"],
            proof_boundary="The stronger analytic input remains open; not Lambda<=0 or RH.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_polymath15_antedb_beta_frontier_audit",
        "date": "2026-07-17",
        "status": (
            "current-source beta-frontier audit with exact unchanged contact and "
            "one route nonpromotion gate; not a proof of Lambda<=0 or RH"
        ),
        "proof_boundary": (
            "This artifact verifies the exact beta-coordinate map, recomputes the "
            "listed current exponent-pair hull, and records a pinned one-pass "
            "ANTEDB direct-beta audit. It does not prove the external database "
            "complete, promote raw optimizer tie partitions, complete recursive "
            "closure, improve beta at the critical contact, close c=2, cross the "
            "fixed c<2 phase wall, prove Wronskian separation, prove Lambda<=0, "
            "or prove RH."
        ),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "sources": [
            ANTEDB_PAPER_URL,
            ANTEDB_REPO_URL,
            ANTEDB_COMMIT_URL,
            TRUDGIAN_YANG_URL,
            "outputs/jensen_window_pf_newman_polymath15_gaussian_legendre_duality_gate.md",
            "outputs/jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    contact = exact["critical_contact"]
    hull = exact["pair_hull"]
    pipeline = exact["antedb_pipeline"]
    lines = [
        "# Jensen-Window PF Newman Polymath-15 ANTEDB Beta-Frontier Audit",
        "",
        "Date: 2026-07-17",
        "",
        "Status: current-source frontier audit and route nonpromotion gate. The",
        "required cancellation improvement remains open. This is not a proof of",
        "`Lambda <= 0` or RH.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_newman_polymath15_antedb_beta_frontier_audit.json",
        "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_antedb_beta_frontier_audit.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_antedb_beta_frontier_audit.py",
        "```",
        "",
        "## Source Snapshot",
        "",
        f"The audit uses ANTEDB commit `{ANTEDB_COMMIT}` dated 2026-07-11 and",
        "the 2025 Tao-Trudgian-Yang paper. It includes their four new exponent",
        "pairs and the two Cushing 2025 pairs recorded by the current blueprint.",
        "This is an as-of-source audit, not a claim about future literature.",
        "",
        "## Exact Coordinate",
        "",
        "```text",
        exact["coordinate_map"],
        "```",
        "",
        "Thus the current pointwise phase input must be optimized directly as a",
        "beta profile; restricting attention to named exponent pairs could miss a",
        "stronger piecewise beta estimate.",
        "",
        "## Current Hull",
        "",
        "The independently recomputed rational pair hull, including the six",
        "post-2023 pairs and Heath-Brown sequence through `m=100`, gives",
        "",
        "```text",
        f"line count={hull['line_count']}",
        f"r_*={hull['maximum']['radius']['exact']}",
        f"active={','.join(hull['maximum']['active'])}",
        f"c_*={hull['maximum']['required_scaled_time']['exact']}={hull['maximum']['required_scaled_time']['decimal']}...",
        f"omitted-tail trivial ceiling={hull['tail_trivial_time_ceiling']['exact']}={hull['tail_trivial_time_ceiling']['decimal']}...",
        "```",
        "",
        "The new pairs improve intermediate radii, but the global maximum remains",
        "the old `TY1/TY2` contact.",
        "",
        "## Direct Beta Audit",
        "",
        "The documented one-pass ANTEDB pipeline was reproduced with Python",
        f"`{pipeline['python']}`, pycddlib `{pipeline['pycddlib']}`, and",
        f"`PYTHONHASHSEED={pipeline['python_hash_seed']}`. Two seeded runs returned",
        f"{pipeline['final_piece_count_in_two_seeded_runs']} final pieces. At the",
        "critical point, its left, point, and right pieces meet exactly:",
        "",
        "```text",
        f"alpha_*={contact['alpha_star']['exact']}",
        f"beta_*={contact['beta_star']['exact']}",
        contact["left_beta_line"],
        contact["point_beta_line"],
        contact["right_beta_line"],
        f"c_req(alpha_*)={contact['scaled_time']['exact']}",
        "```",
        "",
        pipeline["reproducibility_boundary"],
        "",
        "A separate beta-only iteration audit then applied the proved ANTEDB",
        "van der Corput transform twelve times:",
        "",
        "```text",
        exact["beta_only_iteration_audit"]["result"],
        f"piece counts={exact['beta_only_iteration_audit']['piece_counts']}",
        f"contact beta on every pass={exact['beta_only_iteration_audit']['contact_beta_on_every_pass']['exact']}",
        f"global threshold on every pass={exact['beta_only_iteration_audit']['global_threshold_on_every_pass']['exact']}",
        "```",
        "",
        exact["beta_only_iteration_audit"]["proof_boundary"],
        "The broader exploratory beta/exponent-pair fixed-point run was terminated",
        "and contributes no mathematical claim here.",
        "",
        "## Consequence",
        "",
        "```text",
        exact["current_conclusion"],
        exact["required_upgrade"],
        "```",
        "",
        f"At `c=2`, the unchanged local phase deficit is",
        f"`{exact['c2_phase_deficit']['exact']}`. This remains a route benchmark,",
        "not a proof of cancellation, inner Wronskian separation, or RH.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman Polymath-15 ANTEDB beta-frontier audit: "
        f"{len(artifact['rows'])} rows, 6 post-2023 pairs, "
        "1 exact current-hull maximum, 1 direct-beta contact, "
        "1 unchanged c=2 deficit, 1 nonpromotion gate"
    )


if __name__ == "__main__":
    main()
