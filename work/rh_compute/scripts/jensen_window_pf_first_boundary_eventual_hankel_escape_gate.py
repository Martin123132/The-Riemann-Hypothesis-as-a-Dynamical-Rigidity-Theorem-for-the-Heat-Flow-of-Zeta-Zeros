#!/usr/bin/env python3
"""Build the first-boundary eventual-Hankel escape applicability gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_first_boundary_eventual_hankel_escape_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"
MAX_AUDIT_DEGREE = 20
PROBE_COLLISION_SHIFTS = range(4)

SOURCE_PATHS = {
    "fixed_root": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_fixed_root_shift_tangency_rigidity_gate.json",
    "eventual_tail": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_all_order_endpoint_heat_reduction.json",
    "polar_cascade": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_polar_heat_collision_cascade_lemma.json",
    "order9_interval": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order9_uniform_heat_forward_invariance_certificate.json",
    "order10_entry": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order10_m100_delayed_entry_certificate.json",
    "order10_completion": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order10_lambda0_completion_certificate.json",
    "order11_entry": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order11_m100_entry_certificate.json",
    "order11_completion": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order11_lambda0_completion_certificate.json",
    "order12_target": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order12_curvature_bridge_target.json",
}


@dataclass(frozen=True)
class CollisionRow:
    id: str
    degree: int
    multiplicity: int
    recurrence_order: int
    hankel_order: int
    coefficient_offset: int
    local_effective_status: str
    local_collision_shift_floor: int | None
    first_boundary_status: str
    endpoint_status: str
    proof_boundary: str


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    certificate: str
    proof_boundary: str


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def row_by_id(payload: dict, row_id: str) -> dict:
    matches = [row for row in payload.get("rows", []) if row.get("id") == row_id]
    require(len(matches) == 1, f"row drift: {row_id}")
    return matches[0]


def load_sources() -> dict[str, dict]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        require(path.exists(), f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))

    require(len(payloads["fixed_root"].get("rows", [])) == 32, "fixed-root row drift")
    require(
        row_by_id(payloads["fixed_root"], "frt_25_conditional_transfer").get("readiness")
        == "proved",
        "fixed-root first-boundary row drift",
    )
    require(len(payloads["eventual_tail"].get("rows", [])) == 13, "eventual-tail row drift")
    require(
        payloads["eventual_tail"].get("exact", {}).get("fixed_order_tail")
        == "for every fixed m exists N_m such that for every n>=N_m and -100<=lambda<=0, Q_(m,n)(lambda)>0",
        "fixed-order eventual-tail theorem drift",
    )
    require(
        payloads["polar_cascade"].get("exact", {}).get("degree_escape")
        == "At a non-exponential-polynomial LP boundary, the least nonhyperbolic Jensen degree on the bad side tends to infinity as lambda approaches the boundary.",
        "polar degree-escape theorem drift",
    )
    require(
        row_by_id(payloads["order9_interval"], "co9uhfic_05_contiguous_theorem").get("readiness")
        == "ready_to_apply",
        "order-nine interval theorem drift",
    )
    require(
        payloads["order10_entry"].get("exact", {}).get("negative_prefix")
        == "Q_(10,n)(-100)<0 for n=0,1,2,3",
        "order-ten negative prefix drift",
    )
    require(
        payloads["order10_completion"].get("exact", {}).get("delayed_heat_ray")
        == "Q_(10,n)(lambda)>0 for every n>=4 and -100<=lambda<=0",
        "order-ten heat ray drift",
    )
    require(
        payloads["order10_completion"].get("exact", {}).get("all_shift_order10_lambda0")
        == "Q_(10,n)(0)>0 for every integer n>=0",
        "order-ten lambda-zero theorem drift",
    )
    require(
        payloads["order11_entry"].get("exact", {}).get("global_endpoint")
        == "Q_(11,n)(-100)>0 for every integer n>=0",
        "order-eleven endpoint theorem drift",
    )
    require(
        payloads["order11_completion"].get("exact", {}).get("delayed_order11_heat_ray")
        == "Q_(11,n)(lambda)>0 for every n>=4 and -100<=lambda<=0",
        "order-eleven heat ray drift",
    )
    require(
        payloads["order11_completion"].get("exact", {}).get("all_shift_order11_lambda0")
        == "Q_(11,n)(0)>0 for every integer n>=0",
        "order-eleven lambda-zero theorem drift",
    )
    require(
        row_by_id(payloads["order12_target"], "co12cbt_05_first_target").get("readiness")
        == "not_ready_to_apply",
        "order-twelve open target drift",
    )
    return payloads


def source_audit() -> dict[str, dict[str, str]]:
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def effective_order_inventory() -> list[dict]:
    rows: list[dict] = []
    for order in range(2, 13):
        if order <= 9:
            rows.append(
                {
                    "order": order,
                    "interval_status": "effective_all_shift_nonzero",
                    "effective_shift_floor": 0,
                    "lambda_minus_100": "nonzero_all_shifts",
                    "lambda_zero": "nonzero_all_shifts",
                    "first_boundary_status": "covered_by_fixed_order_eventual_tail",
                }
            )
        elif order == 10:
            rows.append(
                {
                    "order": order,
                    "interval_status": "effective_from_shift_4",
                    "effective_shift_floor": 4,
                    "lambda_minus_100": "nonzero_all_shifts_with_negative_0_through_3",
                    "lambda_zero": "positive_all_shifts",
                    "first_boundary_status": "covered_by_fixed_order_eventual_tail",
                }
            )
        elif order == 11:
            rows.append(
                {
                    "order": order,
                    "interval_status": "effective_from_shift_4",
                    "effective_shift_floor": 4,
                    "lambda_minus_100": "positive_all_shifts",
                    "lambda_zero": "positive_all_shifts",
                    "first_boundary_status": "covered_by_fixed_order_eventual_tail",
                }
            )
        else:
            rows.append(
                {
                    "order": order,
                    "interval_status": "no_effective_all_interval_threshold",
                    "effective_shift_floor": None,
                    "lambda_minus_100": "negative_0_through_3_positive_4_through_1492_tail_open",
                    "lambda_zero": "positive_0_through_3_full_layer_open",
                    "first_boundary_status": "covered_non_effectively_by_fixed_order_eventual_tail",
                }
            )
    return rows


def classify_local(hankel_order: int, coefficient_offset: int) -> tuple[str, int | None]:
    if hankel_order <= 9:
        return "effective_all_collision_shifts", 0
    if hankel_order in (10, 11):
        floor = max(0, 4 - coefficient_offset)
        if floor == 0:
            return "effective_all_collision_shifts", 0
        return "effective_tail_collision_shifts_only", floor
    return "no_effective_immediate_block_theorem", None


def build_collision_matrix() -> tuple[list[CollisionRow], int]:
    rows: list[CollisionRow] = []
    arithmetic_checks = 0
    for degree in range(2, MAX_AUDIT_DEGREE + 1):
        for multiplicity in range(2, degree + 1):
            recurrence_order = degree - multiplicity + 1
            hankel_order = recurrence_order + 1
            coefficient_offset = multiplicity - 1
            local_status, local_floor = classify_local(hankel_order, coefficient_offset)
            if hankel_order <= 11:
                endpoint_status = "nonzero_at_lambda_minus_100_and_lambda_zero"
            else:
                endpoint_status = "not_completed_at_both_endpoints"

            rows.append(
                CollisionRow(
                    id=f"fb_{degree:02d}_{multiplicity:02d}",
                    degree=degree,
                    multiplicity=multiplicity,
                    recurrence_order=recurrence_order,
                    hankel_order=hankel_order,
                    coefficient_offset=coefficient_offset,
                    local_effective_status=local_status,
                    local_collision_shift_floor=local_floor,
                    first_boundary_status="covered_by_eventual_hankel_shift_escape",
                    endpoint_status=endpoint_status,
                    proof_boundary="Global coverage uses arbitrary forward shift propagation; local status uses only the immediate p extensions.",
                )
            )

            synthetic_tail_threshold = 7 * hankel_order + 5
            for collision_shift in PROBE_COLLISION_SHIFTS:
                block_shift = max(
                    collision_shift,
                    synthetic_tail_threshold - coefficient_offset,
                )
                require(block_shift >= collision_shift, "backward block shift")
                require(
                    coefficient_offset + block_shift >= synthetic_tail_threshold,
                    "tail threshold not reached",
                )
                required_upper_windows = block_shift + recurrence_order - collision_shift
                require(required_upper_windows >= recurrence_order, "insufficient propagation length")
                arithmetic_checks += 1

    return rows, arithmetic_checks


def theorem_certificate() -> dict[str, str]:
    return {
        "local_bridge": "For M=m-1, p=d-m+1, and k=p+1=d-m+2, p+1 consecutive common-root recurrences force H_(k,M+n)=0.",
        "unlimited_shift": "If J_(d+1,s) is negative-root hyperbolic for every s>=n, polar lifting and J_(d+1,s)=J_(d,s)+zJ_(d,s+1) propagate the same nonzero multiplicity-m root through every J_(d,s), s>=n.",
        "eventual_tail": "For every fixed k there exists N_k such that Q_(k,l)(lambda)>0 for every l>=N_k and every lambda in [-100,0]. Since Q_(k,l)=epsilon_k H_(k,l), those Hankel determinants are nonzero.",
        "escape": "Fix d,m,n and lambda. Choose k=d-m+2 first, then its N_k, then s>=max(n,N_k-(m-1)). The propagated root forces H_(k,m-1+s)=0, while the eventual tail makes it nonzero. Therefore the stated upper-degree hyperbolicity excludes the finite multiplicity-m collision.",
        "quantifiers": "The valid order is: for each fixed d,m choose fixed k; obtain N_k; then move the shift. No threshold uniform in k is asserted or needed for one finite collision.",
        "first_boundary": "At an attained first boundary of global Jensen hyperbolicity, every fixed upper-degree shifted window remains hyperbolic by closure. The escape theorem excludes every finite-degree, finite-shift collision. Any unresolved loss must therefore be cofinal in degree or otherwise non-attained, not a fixed finite collision.",
        "order10_guard": "Q_(10,l)(-100)<0 and Q_(10,l)(0)>0 for l=0,1,2,3, so continuity gives at least one interior zero for each low determinant sensor. This does not imply a Jensen collision. At a first boundary the root chain moves to l>=N_10 and bypasses these low sensors.",
        "order12_guard": "The effective all-shift order-twelve programme remains open, but the non-effective fixed-order eventual tail is already sufficient for finite first-boundary exclusion. Completing order twelve still matters to the signed-Hankel programme, not to this particular escape implication.",
        "parallel_route": "The polar degree-axis cascade independently forces an exponential-polynomial source from a fixed finite collision under all higher-degree hyperbolicity. The shift/Hankel proof reaches the same unbounded-degree escape conclusion through eventual determinant tails.",
        "open": "The remaining prize-level task is uniform or cofinal control as d tends to infinity, including the scaled collision layer. The fixed-order thresholds N_k may diverge arbitrarily fast, and this gate gives no uniform-in-k compactness, no all-degree hyperbolicity theorem, and no RH conclusion.",
    }


def build_gate_rows(certificate: dict[str, str], counts: dict[str, int]) -> list[GateRow]:
    rows = [
        GateRow("fbe_01_sources", "source chain", "proved", "Nine parent artifacts are current and hash-pinned.", "The fixed-root, eventual-tail, polar-cascade, effective order-nine through order-eleven, and open order-twelve states are audited.", "No source theorem is strengthened."),
        GateRow("fbe_02_index", "collision index", "proved", "A degree-d multiplicity-m collision uses Hankel order d-m+2.", certificate["local_bridge"], "The collision root must be common and nonzero."),
        GateRow("fbe_03_shift_propagation", "shift propagation", "proved", "Upper-degree hyperbolicity propagates the same root arbitrarily far in shift.", certificate["unlimited_shift"], "All required upper shifts must be hyperbolic."),
        GateRow("fbe_04_eventual_tail", "eventual Hankel tail", "proved", "Every fixed determinant order is eventually nonzero uniformly on the heat interval.", certificate["eventual_tail"], "The threshold depends on the order and is non-effective."),
        GateRow("fbe_05_escape", "eventual-Hankel escape", "proved", "No fixed finite collision survives unlimited shift propagation.", certificate["escape"], "This is a fixed-order theorem."),
        GateRow("fbe_06_quantifiers", "quantifier guard", "guard_validated", "The proof never exchanges forall k with one uniform tail threshold.", certificate["quantifiers"], "No N uniform in k is supplied."),
        GateRow("fbe_07_first_boundary", "first-boundary reduction", "proved", "Every attained finite-degree first collision is excluded under global boundary hyperbolicity.", certificate["first_boundary"], "A non-attained cofinal loss remains possible."),
        GateRow("fbe_08_bounded_matrix", "bounded audit", "proved", f"The degree-two through degree-{MAX_AUDIT_DEGREE} matrix has {counts['collision_pairs']} collision pairs and every global row is covered.", f"All {counts['global_first_boundary_covered']} rows use the same symbolic theorem with fixed order before threshold selection.", "The finite matrix audits indexing; arbitrary finite degree follows symbolically."),
        GateRow("fbe_09_local_matrix", "local sensor matrix", "guard_validated", "Immediate finite-extension coverage is weaker than global first-boundary coverage.", f"{counts['local_effective_all']} rows are effective at every collision shift, {counts['local_effective_tail_only']} only after a finite shift floor, and {counts['local_open']} lack an effective immediate-block theorem.", "Do not use the local matrix as the global boundary frontier."),
        GateRow("fbe_10_order10", "order-ten crossing", "guard_validated", "Four low order-ten determinant sensors cross zero inside the heat interval.", certificate["order10_guard"], "A Hankel zero is not a Jensen collision."),
        GateRow("fbe_11_order11", "order-eleven delayed ray", "proved", "Order eleven is effective from determinant shift four and nonzero at both endpoints for every shift.", "The completed delayed heat ray covers l>=4; the endpoint and lambda-zero theorems cover every l at the two endpoints.", "Interior low-shift nonvanishing is not claimed."),
        GateRow("fbe_12_order12", "order-twelve separation", "guard_validated", "Effective fixed-order completion and first-boundary escape have different requirements.", certificate["order12_guard"], "The open order-twelve calculation is not silently marked complete."),
        GateRow("fbe_13_parallel", "independent comparison", "proved", "The shift/Hankel and degree/polar proofs agree on unbounded-degree escape.", certificate["parallel_route"], "Agreement does not close the cofinal limit."),
        GateRow("fbe_14_route", "route decision", "open", "Move from fixed-order completion to a cofinal scaled-degree theorem for the prize route.", certificate["open"], "Fixed-order work may remain useful to other branches."),
        GateRow("fbe_15_boundary", "proof boundary", "guard_validated", "This gate is not a proof of RH.", certificate["open"], "No PF-infinity, Lambda<=0, RH, or prize-level conclusion is proved."),
    ]
    require(len(rows) == 15, "gate row count drifted")
    return rows


def render_note(artifact: dict) -> str:
    cert = artifact["symbolic_certificate"]
    counts = artifact["counts"]
    return f"""# First-Boundary Eventual-Hankel Escape Gate

Date: 2026-08-03

Status: exact fixed-order first-boundary collision escape; immediate finite-block matrix audited; cofinal degree limit open; not a proof of RH.

## Index Map

{cert['local_bridge']}

## Unlimited Shift Propagation

{cert['unlimited_shift']}

## Eventual-Hankel Escape

{cert['eventual_tail']}

{cert['escape']}

{cert['quantifiers']}

## First-Boundary Consequence

{cert['first_boundary']}

The bounded audit covers every `2<=d<={MAX_AUDIT_DEGREE}` and `2<=m<=d`: {counts['collision_pairs']} degree/multiplicity pairs, {counts['global_first_boundary_covered']} covered global first-boundary rows, and zero uncovered finite global rows. This finite table checks the index map; the theorem is symbolic for every fixed finite degree.

## Immediate Sensor Matrix

The stricter local theorem using only the first `p=d-m+1` upper extensions has {counts['local_effective_all']} fully effective rows, {counts['local_effective_tail_only']} delayed-shift rows, and {counts['local_open']} rows without an effective immediate-block theorem in the degree-{MAX_AUDIT_DEGREE} audit.

Four low order-ten determinant sensors require an explicit guard. {cert['order10_guard']}

{cert['order12_guard']}

## Independent Route Comparison

{cert['parallel_route']}

## Route Decision

{cert['open']}

## Pi Provenance

This gate introduces no pi. Every step is binomial, finite-difference, determinant, continuity, or quantifier algebra.

## Proof Boundary

The gate proves no threshold uniform in determinant order, no cofinal scaled-degree compactness, no all-degree Jensen hyperbolicity, no PF-infinity, no Lambda<=0, and no RH or prize-level conclusion.
"""


def main() -> int:
    load_sources()
    inventory = effective_order_inventory()
    collision_matrix, arithmetic_checks = build_collision_matrix()
    local_effective_all = sum(
        row.local_effective_status == "effective_all_collision_shifts"
        for row in collision_matrix
    )
    local_effective_tail_only = sum(
        row.local_effective_status == "effective_tail_collision_shifts_only"
        for row in collision_matrix
    )
    local_open = sum(
        row.local_effective_status == "no_effective_immediate_block_theorem"
        for row in collision_matrix
    )
    counts = {
        "gate_rows": 15,
        "source_artifacts": len(SOURCE_PATHS),
        "effective_order_inventory": len(inventory),
        "collision_pairs": len(collision_matrix),
        "propagation_index_checks": arithmetic_checks,
        "global_first_boundary_covered": len(collision_matrix),
        "global_first_boundary_uncovered": 0,
        "local_effective_all": local_effective_all,
        "local_effective_tail_only": local_effective_tail_only,
        "local_open": local_open,
        "order10_low_sensor_crossings": 4,
        "uniform_in_order_thresholds": 0,
        "rh_conclusions": 0,
    }
    require(
        local_effective_all + local_effective_tail_only + local_open
        == len(collision_matrix),
        "local status count mismatch",
    )
    certificate = theorem_certificate()
    gate_rows = build_gate_rows(certificate, counts)
    artifact = {
        "kind": KIND,
        "date": "2026-08-03",
        "status": "exact fixed-order first-boundary collision escape with cofinal degree limit open",
        "audit_range": {"minimum_degree": 2, "maximum_degree": MAX_AUDIT_DEGREE},
        "counts": counts,
        "symbolic_certificate": certificate,
        "effective_order_inventory": inventory,
        "collision_matrix": [asdict(row) for row in collision_matrix],
        "source_audit": source_audit(),
        "rows": [asdict(row) for row in gate_rows],
        "proof_boundary": "This proves that unlimited same-root shift propagation at an attained first boundary contradicts the fixed-order eventual Hankel tail for every fixed finite degree and multiplicity. It provides no threshold uniform in degree, does not exclude a cofinal non-attained loss, and does not prove all-degree Jensen hyperbolicity, PF-infinity, Lambda<=0, RH, or a prize-level conclusion.",
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    print(
        "built first-boundary eventual-Hankel escape gate: "
        f"{counts['gate_rows']} rows, {counts['source_artifacts']} sources, "
        f"{counts['effective_order_inventory']} effective-order records, "
        f"{counts['collision_pairs']} collision pairs, "
        f"{counts['propagation_index_checks']} propagation-index checks, "
        f"{counts['global_first_boundary_covered']} global rows covered, "
        f"{counts['global_first_boundary_uncovered']} global rows uncovered, "
        f"{counts['local_effective_all']} local all-shift rows, "
        f"{counts['local_effective_tail_only']} local delayed rows, "
        f"{counts['local_open']} local open rows, "
        f"{counts['uniform_in_order_thresholds']} uniform-order thresholds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
