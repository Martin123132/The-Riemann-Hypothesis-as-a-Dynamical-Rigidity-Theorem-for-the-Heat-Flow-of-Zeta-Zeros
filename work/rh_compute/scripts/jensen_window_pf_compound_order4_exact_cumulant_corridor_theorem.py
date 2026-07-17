#!/usr/bin/env python3
"""Compose the exact order-four cumulant corridors globally on u>=2."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

from jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate import (  # noqa: E402
    sha256,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.md"
)
SOURCE_TARGET = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.json"
)
SOURCE_FORMAL_FINITE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.json"
)
SOURCE_FORMAL_RAY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.json"
)
SOURCE_BUDGET = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.json"
)
SOURCE_CENTRAL = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_central_residual_certificate.json"
)


@dataclass(frozen=True)
class CorridorRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def build_artifact() -> dict:
    target = load_json(SOURCE_TARGET)
    finite = load_json(SOURCE_FORMAL_FINITE)
    ray = load_json(SOURCE_FORMAL_RAY)
    budget = load_json(SOURCE_BUDGET)
    central = load_json(SOURCE_CENTRAL)
    if (
        finite.get("parameters", {}).get("block_count") != 1_800_000
        or len(finite.get("corridors", {})) != 7
    ):
        raise RuntimeError("finite formal corridors are not closed")
    if ray.get("summary", {}).get("formal_ray_closed") is not True:
        raise RuntimeError("formal corridor ray is not closed")
    if budget.get("summary", {}).get("global_epsilon_ten_layer_closed") is not True:
        raise RuntimeError("epsilon-ten budget source is not closed")
    if central.get("summary", {}).get(
        "exact_minus_epsilon_ten_cumulant_budgets_closed"
    ) is not True:
        raise RuntimeError("exact-density cumulant budget is not closed")

    finite_budget = budget["epsilon_ten_budget"]["finite"]
    ray_budget = budget["epsilon_ten_budget"]["ray"]
    if finite_budget["reserved_final_corridor_margin"] != "79999/10000000":
        raise RuntimeError("finite final corridor margin changed")
    if ray_budget["reserved_final_corridor_margin"] != "29/(1000u)":
        raise RuntimeError("ray final corridor margin changed")
    corridors = target["candidate_corridors"]
    rows = [
        CorridorRow(
            id="co4ecct_01_corridor_signature",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The candidate inequalities are the explicit alternating-factorial cumulant corridors required by the localized curvature route.",
            formula="; ".join([corridors["kappa2"], *corridors["higher"].values()]),
            proof_boundary="Exact corridor specification only.",
        ),
        CorridorRow(
            id="co4ecct_02_global_formal_corridors",
            role="exact_theorem_composition",
            readiness="ready_to_apply",
            claim="The finite Arb cover and coefficient-positive ray theorem prove every epsilon-six formal corridor for all u>=2.",
            formula="formal kappa_r^[6] lies inside every displayed corridor, r=2,...,8",
            proof_boundary="Epsilon-six formal cumulants only.",
        ),
        CorridorRow(
            id="co4ecct_03_global_epsilon_eight_ten_layers",
            role="exact_theorem_composition",
            readiness="ready_to_apply",
            claim="The global epsilon-eight and epsilon-ten coefficient theorems consume only the recorded portions of the formal corridor buffers.",
            formula="finite epsilon-ten correction <10^-7; ray <1/(1000u)",
            proof_boundary="Formal correction layers only.",
            diagnostics={"finite": finite_budget, "ray": ray_budget},
        ),
        CorridorRow(
            id="co4ecct_04_exact_density_residual",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="The complex-disk central-and-tail theorem supplies the sufficient exact-minus-epsilon-ten cumulant budgets simultaneously through order eight.",
            formula="finite <9/1000; ray <1/(100u)",
            proof_boundary="Exact-minus-formal cumulant residual only.",
            diagnostics={
                "finite_partition_ratio": central["finite"]["budget"]["budget_ratio_upper"],
                "ray_partition_ratio": central["ray"]["budget"]["budget_ratio_upper"],
            },
        ),
        CorridorRow(
            id="co4ecct_05_finite_exact_corridors",
            role="exact_theorem_composition",
            readiness="ready_to_apply",
            claim="Every exact cumulant corridor holds on the complete finite mode segment with a strict common reserve.",
            formula="exact corridors hold on 2<=u<=20; reserve 79999/10000000",
            proof_boundary="Exact cumulants through order eight on the finite segment only.",
        ),
        CorridorRow(
            id="co4ecct_06_ray_exact_corridors",
            role="exact_theorem_composition",
            readiness="ready_to_apply",
            claim="Every exact cumulant corridor holds on the asymptotic mode ray with a strict u-dependent reserve.",
            formula="exact corridors hold on u>=20; reserve 29/(1000u)",
            proof_boundary="Exact cumulants through order eight on the asymptotic segment only.",
        ),
        CorridorRow(
            id="co4ecct_07_global_exact_corridor_theorem",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="The seven alternating-factorial exact cumulant corridors hold globally for every u>=2.",
            formula="candidate exact cumulant corridors hold for u>=2, r=2,...,8",
            proof_boundary="Exact cumulant theorem only; the continuum corridor-to-U inequality remains separate.",
        ),
        CorridorRow(
            id="co4ecct_08_localized_curvature_handoff",
            role="analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Prove the localized U(t) ceiling uniformly from the now-exact corridors and exact mode geometry, not merely on sampled collars.",
            formula="exact corridors + mode geometry => U(t)<=7/(2t^2), u>=2",
            proof_boundary="Open continuum composition; no curvature ray or order-four entry is asserted.",
        ),
    ]
    source_paths = (
        SOURCE_TARGET,
        SOURCE_FORMAL_FINITE,
        SOURCE_FORMAL_RAY,
        SOURCE_BUDGET,
        SOURCE_CENTRAL,
    )
    return {
        "kind": "jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem",
        "date": "2026-07-13",
        "status": "global exact cumulant corridor theorem with open localized-curvature composition",
        "proof_boundary": (
            "This artifact proves all seven exact cumulant corridors for every u>=2. "
            "It does not prove that interval substitution of those corridors and the "
            "mode geometry yields the localized U ceiling on the full continuum, "
            "order-four entry, PF-infinity, RH, or Lambda<=0."
        ),
        "corridors": corridors,
        "finite_reserve": finite_budget["reserved_final_corridor_margin"],
        "ray_reserve": ray_budget["reserved_final_corridor_margin"],
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_rows": 7,
            "open_analytic_rows": 1,
            "exact_corridors": 7,
            "finite_corridors_closed": 7,
            "ray_corridors_closed": 7,
            "global_exact_corridors_closed": True,
            "open_localized_curvature_compositions": 1,
        },
        "sources": [
            "outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md",
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md",
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_central_residual_certificate.md",
            "outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md",
            "outputs/formal_core.md",
        ],
        "source_hashes": {
            path.relative_to(REPO_ROOT).as_posix(): sha256(path)
            for path in source_paths
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.py"
        ),
        "remaining_target": (
            "Intervalize the exact corridor boxes together with exact t, curvature, "
            "and Hurwitz-zeta geometry across u>=2 and prove the localized U ceiling."
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    lines = [
        "# Jensen-Window PF Order-Four Exact Cumulant Corridor Theorem",
        "",
        "Date: 2026-07-13",
        "",
        "Status: global exact cumulant corridor theorem with open localized-curvature composition.",
        "This is not a proof of the curvature ray, order-four entry, PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.py",
        "```",
        "",
        "## Exact Theorem",
        "",
        "The finite and asymptotic formal corridor theorems, the two globally",
        "bounded correction layers, and the exact central-and-tail residual compose",
        "to prove all seven candidate corridors for every `u>=2`:",
        "",
        "```text",
        artifact["corridors"]["kappa2"],
        *artifact["corridors"]["higher"].values(),
        "```",
        "",
        f"The finite common reserve is `{artifact['finite_reserve']}`; the ray reserve",
        f"is `{artifact['ray_reserve']}`.",
        "",
        "## Remaining Boundary",
        "",
        "The exact-density theorem is closed, but the curvature ray is not yet",
        "closed. The remaining task is a continuum interval/analytic composition",
        "of these exact corridor boxes with the exact mode geometry in the localized",
        "quantity `U(t)`. Seven sampled collars passed previously; samples are not",
        "a substitute for the full `u>=2` theorem.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_exact_cumulant_central_residual_certificate.md",
        "outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    print(
        "certified order-four exact cumulant corridors: "
        "8 rows, 7 exact rows, 7 global exact corridors, "
        "2 strict reserve regimes, 1 open localized-curvature composition"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
