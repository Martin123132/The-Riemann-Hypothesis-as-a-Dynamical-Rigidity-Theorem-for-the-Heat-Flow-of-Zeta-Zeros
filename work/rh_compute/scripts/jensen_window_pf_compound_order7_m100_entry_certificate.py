#!/usr/bin/env python3
"""Compose the all-shift signed order-seven entry theorem at lambda=-100."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order7_m100_entry_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_compound_order7_m100_entry_certificate.md"
)
PREFIX_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order7_m100_prefix_certificate.json"
)
TAIL_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order7_m100_tail_curvature_reduction.json"
)
BRIDGE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order7_first_summand_curvature_bridge.json"
)
SHIFTED_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate.json"
)
COMPACT_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order7_nested_curvature_compact_certificate.json"
)
FINITE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order7_nested_curvature_finite_ray_certificate.json"
)
ASYMPTOTIC_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order7_nested_curvature_asymptotic_ray_certificate.json"
)


@dataclass(frozen=True)
class EntryRow:
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


def row_by_id(artifact: dict, row_id: str) -> dict:
    row = next(
        (candidate for candidate in artifact.get("rows", []) if candidate.get("id") == row_id),
        None,
    )
    if row is None or row.get("readiness") != "ready_to_apply":
        raise RuntimeError(f"required source row {row_id} is unavailable")
    return row


def validate_sources() -> dict:
    prefix = load_json(PREFIX_SOURCE)
    tail = load_json(TAIL_SOURCE)
    bridge = load_json(BRIDGE_SOURCE)
    shifted = load_json(SHIFTED_SOURCE)
    compact = load_json(COMPACT_SOURCE)
    finite = load_json(FINITE_SOURCE)
    asymptotic = load_json(ASYMPTOTIC_SOURCE)

    if prefix.get("summary", {}).get("positive_Q7_rows") != 315:
        raise RuntimeError("order-seven prefix source changed")
    if prefix.get("exact", {}).get("prefix") != (
        "Q_(7,n)(-100)>0 for every 0<=n<=314"
    ):
        raise RuntimeError("order-seven prefix theorem changed")
    if tail.get("exact", {}).get("sufficient_ceiling") != (
        "R_k<=900/k^2 for every real/integer k>=321"
    ):
        raise RuntimeError("order-seven scalar tail reduction changed")
    if tail.get("exact", {}).get("tail_index") != (
        "k=n+6, so n>=315 iff k>=321"
    ):
        raise RuntimeError("order-seven tail indexing changed")
    full_transfer = bridge.get("exact", {}).get("full_transfer")
    if full_transfer != (
        "|R_k-R_k^(1)|<=N_(k-1)+2*N_k+N_(k+1)<262/k^2, k>=321"
    ):
        raise RuntimeError("order-seven full-kernel transfer changed")
    if bridge.get("exact", {}).get("continuous_target") != (
        "r_1''(t)<=600/t^2 for every real t>=320"
    ):
        raise RuntimeError("order-seven continuous bridge target changed")

    shifted_row = row_by_id(shifted, "co7sjb_04_continuous_curvature")
    compact_row = row_by_id(compact, "co7ncc_03_compact_theorem")
    finite_row = row_by_id(finite, "co7ncfr_03_dimensionless_cover")
    asymptotic_row = row_by_id(asymptotic, "co7ncarc_06_asymptotic_consequence")
    coverage_row = row_by_id(asymptotic, "co7ncarc_07_global_handoff")
    expected_ranges = {
        "shifted": "r_1''(t)<=600/t^2 for 320<=t<=1000",
        "compact": "r_1''(t)<=600/t^2 for every real 1000<=t<=V'(2)",
        "finite_ray": "r_1''(t)<=600/t^2 for every saddle mode 2<=u<=20",
        "asymptotic_ray": "t^2*r_1''(t)<100<600 for u>=20",
        "coverage": (
            "[320<=t<=1000] union [1000<=t<=V'(2)] union "
            "[2<=u<=20] union [u>=20]"
        ),
    }
    observed_ranges = {
        "shifted": shifted_row.get("formula"),
        "compact": compact_row.get("formula"),
        "finite_ray": finite_row.get("formula"),
        "asymptotic_ray": asymptotic_row.get("formula"),
        "coverage": coverage_row.get("formula"),
    }
    for key, expected in expected_ranges.items():
        if observed_ranges[key] != expected:
            raise RuntimeError(f"order-seven {key} source changed")
    if shifted.get("summary", {}).get("continuous_curvature_theorems") != 1:
        raise RuntimeError("order-seven shifted curvature source is not closed")
    if compact.get("summary", {}).get("compact_curvature_theorems") != 1:
        raise RuntimeError("order-seven compact curvature source is not closed")
    if finite.get("summary", {}).get("finite_ray_theorems") != 1:
        raise RuntimeError("order-seven finite-ray source is not closed")
    if asymptotic.get("summary", {}).get("asymptotic_curvature_theorems") != 1:
        raise RuntimeError("order-seven asymptotic source is not closed")

    return {
        "prefix": prefix["exact"]["prefix"],
        **observed_ranges,
        "asymptotic_scaled_upper": asymptotic["dimensionless_interval"][
            "scaled_curvature_upper"
        ],
        "full_transfer": full_transfer,
        "complete_ceiling": tail["exact"]["sufficient_ceiling"],
        "conditional_tail": bridge["exact"]["conditional_endpoint_tail"],
    }


def exact_composition(sources: dict) -> dict:
    return {
        "global_first_curvature": (
            "r_1''(t)<=600/t^2 for every real t>=320"
        ),
        "coverage": sources["coverage"] + " covers every real t>=320",
        "tent_transfer": (
            "R_k^(1)<600*[-log(1-1/k^2)]<601/k^2, k>=321"
        ),
        "full_ceiling": (
            "R_k<R_k^(1)+|R_k-R_k^(1)|<601/k^2+262/k^2="
            "863/k^2<900/k^2, k>=321"
        ),
        "tail_sign": (
            "R_k<900/k^2 => Q_(7,n)(-100)>0, k=n+6, n>=315"
        ),
        "all_shift_entry": (
            "Q_(7,n)(-100)=-H_(7,n)(-100)>0 for every integer n>=0"
        ),
        "next_handoff": (
            "compose all-shift lambda=-100 entry with the completed uniform "
            "order-seven tail and cooperative flow"
        ),
        "asymptotic_reserve": (
            f"u>=20 scaled upper {sources['asymptotic_scaled_upper']}<100<600"
        ),
    }


def build_artifact() -> dict:
    sources = validate_sources()
    exact = exact_composition(sources)
    rows = [
        EntryRow(
            "co7m100ec_01_stable_coordinate",
            "exact_identity_input",
            "ready_to_apply",
            "The signed order-seven sign is controlled by one centered stable curvature coordinate.",
            "Q_(7,n)>0 iff L_n>0 iff E_n<0, k=n+6",
            "Previously proved exact factorization inside the positive order-six cone.",
        ),
        EntryRow(
            "co7m100ec_02_shifted_curvature",
            "interval_theorem_input",
            "ready_to_apply",
            "The shifted-jet certificate closes the continuous endpoint and intermediate compact range.",
            sources["shifted"],
            "Outward-rounded shifted H-jet blocks for the first Newman summand.",
        ),
        EntryRow(
            "co7m100ec_03_compact_curvature",
            "interval_theorem_input",
            "ready_to_apply",
            "The common-collar nested certificate continues the theorem through saddle mode two.",
            sources["compact"],
            "Outward-rounded rational mode blocks from hashed H2-H12 caches.",
        ),
        EntryRow(
            "co7m100ec_04_finite_ray",
            "interval_theorem_input",
            "ready_to_apply",
            "The exact-corridor dimensionless cover closes the finite mode ray.",
            sources["finite_ray"],
            "One exact collar and 17999 rational mode blocks.",
        ),
        EntryRow(
            "co7m100ec_05_asymptotic_ray",
            "interval_analytic_theorem_input",
            "ready_to_apply",
            "The normalized-H interval closes the complete asymptotic mode ray.",
            sources["asymptotic_ray"] + "; " + exact["asymptotic_reserve"],
            "Uniform z in [0,10^-30] interval theorem.",
        ),
        EntryRow(
            "co7m100ec_06_global_curvature",
            "exact_theorem_composition",
            "ready_to_apply",
            "The four continuous certificates cover the complete first-summand tail.",
            exact["global_first_curvature"] + "; " + exact["coverage"],
            "Continuous first summand only.",
        ),
        EntryRow(
            "co7m100ec_07_first_discrete",
            "exact_analytic_theorem",
            "ready_to_apply",
            "The tent identity transfers the real-parameter theorem to discrete first-summand curvature.",
            exact["tent_transfer"],
            "Exact tent integration and an elementary logarithm bound.",
        ),
        EntryRow(
            "co7m100ec_08_full_kernel",
            "exact_perturbation_theorem",
            "ready_to_apply",
            "The proved higher-theta transfer closes the complete-kernel scalar ceiling.",
            sources["full_transfer"] + "; " + exact["full_ceiling"],
            "Complete Newman kernel at lambda=-100, k>=321.",
        ),
        EntryRow(
            "co7m100ec_09_tail",
            "exact_tail_theorem",
            "ready_to_apply",
            "The scalar ceiling proves every missing signed order-seven tail sign.",
            exact["tail_sign"],
            "Analytic tail n>=315 only.",
        ),
        EntryRow(
            "co7m100ec_10_all_shift_entry",
            "exact_entry_theorem",
            "ready_to_apply",
            "The rigorous prefix and analytic tail prove signed contiguous order-seven entry at every shift.",
            sources["prefix"] + "; " + exact["all_shift_entry"],
            "One heat parameter only; not forward propagation or PF-infinity.",
        ),
        EntryRow(
            "co7m100ec_11_forward_handoff",
            "theorem_handoff",
            "ready_to_apply",
            "The endpoint input required by the proved conditional cooperative-flow theorem is now available.",
            exact["next_handoff"],
            "Composition handoff; the forward theorem is recorded separately.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order7_m100_entry_certificate",
        "date": "2026-07-13",
        "status": "all-shift signed contiguous order-seven entry theorem at lambda=-100",
        "proof_boundary": (
            "This artifact proves Q_(7,n)(-100)=-H_(7,n)(-100)>0 for every "
            "n>=0. It does not by itself prove forward order-seven invariance, "
            "arbitrary-column order seven, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_compound_order7_m100_prefix_certificate.md",
            "outputs/jensen_window_pf_compound_order7_m100_tail_curvature_reduction.md",
            "outputs/jensen_window_pf_compound_order7_first_summand_curvature_bridge.md",
            "outputs/jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate.md",
            "outputs/jensen_window_pf_compound_order7_nested_curvature_compact_certificate.md",
            "outputs/jensen_window_pf_compound_order7_nested_curvature_finite_ray_certificate.md",
            "outputs/jensen_window_pf_compound_order7_nested_curvature_asymptotic_ray_certificate.md",
            "outputs/formal_core.md",
        ],
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order7_m100_entry_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order7_m100_entry_certificate.py",
        "source_contract": sources,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": len(rows),
            "continuous_curvature_theorems": 1,
            "complete_scalar_ceiling_theorems": 1,
            "analytic_tail_theorems": 1,
            "all_shift_m100_entry_theorems": 1,
            "open_rows": 0,
        },
    }


def write_note(path: Path, artifact: dict) -> None:
    sources = artifact["source_contract"]
    exact = artifact["exact"]
    lines = [
        "# Jensen-Window PF Compound Order-Seven Lambda=-100 Entry Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: all-shift signed contiguous order-seven entry theorem at",
        "`lambda=-100`. This certificate is not a proof of forward order-seven",
        "invariance, arbitrary-column order seven, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order7_m100_entry_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order7_m100_entry_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order7_m100_entry_certificate.py",
        "```",
        "",
        "## Continuous First-Summand Theorem",
        "",
        "The fourth-nested stable curvature is covered on four closed ranges:",
        "",
        "```text",
        sources["shifted"],
        sources["compact"],
        sources["finite_ray"],
        sources["asymptotic_ray"],
        "```",
        "",
        f"The asymptotic scaled upper is `{sources['asymptotic_scaled_upper']}`.",
        "The certified mode-to-t coverage gives",
        "",
        "```text",
        exact["global_first_curvature"],
        exact["coverage"],
        "```",
        "",
        "## Discrete And Full-Kernel Transfer",
        "",
        "```text",
        exact["tent_transfer"],
        sources["full_transfer"],
        exact["full_ceiling"],
        "```",
        "",
        "The last line is strictly below the scalar ceiling isolated by the",
        "tail reduction, so",
        "",
        "```text",
        exact["tail_sign"],
        "```",
        "",
        "## Entry Theorem",
        "",
        "The 1024-bit prefix proves",
        "",
        "```text",
        sources["prefix"],
        "```",
        "",
        "Combining prefix and tail gives",
        "",
        "```text",
        exact["all_shift_entry"],
        "```",
        "",
        "The sole endpoint input in the conditional order-seven heat-flow",
        "reduction is therefore available for the forward composition.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order7_m100_prefix_certificate.md",
        "outputs/jensen_window_pf_compound_order7_first_summand_curvature_bridge.md",
        "outputs/jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate.md",
        "outputs/jensen_window_pf_compound_order7_nested_curvature_compact_certificate.md",
        "outputs/jensen_window_pf_compound_order7_nested_curvature_finite_ray_certificate.md",
        "outputs/jensen_window_pf_compound_order7_nested_curvature_asymptotic_ray_certificate.md",
        "outputs/jensen_window_pf_compound_order7_uniform_tail_flow_reduction.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-seven lambda=-100 entry certificate: "
        f"{summary['rows']} rows, "
        f"{summary['continuous_curvature_theorems']} continuous curvature theorem, "
        f"{summary['analytic_tail_theorems']} analytic tail theorem, "
        f"{summary['all_shift_m100_entry_theorems']} all-shift entry theorem, "
        f"{summary['open_rows']} open rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
