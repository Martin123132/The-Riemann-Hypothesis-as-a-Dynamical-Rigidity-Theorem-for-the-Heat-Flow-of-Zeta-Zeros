#!/usr/bin/env python3
"""Compose the all-shift signed order-six entry theorem at lambda=-100."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order6_m100_entry_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_compound_order6_m100_entry_certificate.md"
)
PREFIX_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order6_m100_prefix_certificate.json"
)
TAIL_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order6_m100_tail_curvature_reduction.json"
)
BRIDGE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order6_first_summand_curvature_bridge.json"
)
COMPACT_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order6_nested_curvature_compact_certificate.json"
)
FINITE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.json"
)
ASYMPTOTIC_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate.json"
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


def validate_sources() -> dict:
    prefix = load_json(PREFIX_SOURCE)
    tail = load_json(TAIL_SOURCE)
    bridge = load_json(BRIDGE_SOURCE)
    compact = load_json(COMPACT_SOURCE)
    finite = load_json(FINITE_SOURCE)
    asymptotic = load_json(ASYMPTOTIC_SOURCE)
    if prefix.get("summary", {}).get("positive_Q6_rows") != 317:
        raise RuntimeError("order-six prefix source changed")
    if tail.get("exact", {}).get("sufficient_ceiling") != (
        "P_k<=320/k^2 for every real/integer k>=322"
    ):
        raise RuntimeError("order-six scalar tail reduction changed")
    full_transfer = bridge.get("exact", {}).get("full_transfer")
    if full_transfer != (
        "|P_k-P_k^(1)|<=Y_(k-1)+2*Y_k+Y_(k+1)<100/k^2, k>=322"
    ):
        raise RuntimeError("order-six full-kernel transfer changed")
    if compact.get("compact", {}).get("certified_t_range") != (
        "321<=t<=V'(2)"
    ):
        raise RuntimeError("order-six compact curvature source changed")
    if finite.get("summary", {}).get("finite_ray_theorems") != 1:
        raise RuntimeError("order-six finite-ray source is not closed")
    if asymptotic.get("summary", {}).get("asymptotic_curvature_theorems") != 1:
        raise RuntimeError("order-six asymptotic source is not closed")
    return {
        "prefix": "Q_(6,n)(-100)>0 for every 0<=n<=316",
        "compact": "p_1''(t)<=200/t^2 for 321<=t<=V'(2)",
        "finite_ray": "p_1''(t)<=200/t^2 for every saddle mode 2<=u<=20",
        "asymptotic_ray": "p_1''(t)<=200/t^2 for every mode u>=20",
        "asymptotic_scaled_upper": asymptotic["dimensionless_interval"][
            "scaled_curvature_upper"
        ],
        "full_transfer": full_transfer,
        "complete_ceiling": tail["exact"]["sufficient_ceiling"],
        "conditional_tail": (
            "[P_k<=320/k^2 for k>=322] => "
            "[Q_(6,n)(-100)>0 for n>=317]"
        ),
    }


def exact_composition(sources: dict) -> dict:
    return {
        "global_first_curvature": (
            "p_1''(t)<=200/t^2 for every real t>=321"
        ),
        "coverage": (
            "[321<=t<=V'(2)] union [2<=u<=20] union [u>=20] "
            "covers every t>=321"
        ),
        "tent_transfer": (
            "P_k^(1)<=200*[-log(1-1/k^2)]<201/k^2, k>=322"
        ),
        "full_ceiling": (
            "P_k<=P_k^(1)+|P_k-P_k^(1)|<201/k^2+100/k^2="
            "301/k^2<320/k^2, k>=322"
        ),
        "tail_sign": (
            "P_k<=320/k^2 => Q_(6,n)(-100)>0, k=n+5, n>=317"
        ),
        "all_shift_entry": (
            "Q_(6,n)(-100)=-H_(6,n)(-100)>0 for every integer n>=0"
        ),
        "target_discharge": "target_compound_order6_m100_entry is discharged",
        "next_handoff": (
            "compose all-shift lambda=-100 entry with the completed uniform "
            "tail and cooperative flow"
        ),
        "asymptotic_reserve": (
            f"u>=20 scaled upper {sources['asymptotic_scaled_upper']}<100<200"
        ),
    }


def build_artifact() -> dict:
    sources = validate_sources()
    exact = exact_composition(sources)
    rows = [
        EntryRow(
            "co6m100ec_01_stable_coordinate",
            "exact_identity_input",
            "ready_to_apply",
            "The signed order-six sign is controlled by one centered stable curvature coordinate.",
            "Q_(6,n)>0 iff D_n=5*log(x_k)+P_k<0, k=n+5",
            "Previously proved exact factorization inside the positive order-five cone.",
        ),
        EntryRow(
            "co6m100ec_02_compact_curvature",
            "interval_theorem_input",
            "ready_to_apply",
            "The continuous first-summand curvature ceiling holds from the splice through mode two.",
            sources["compact"],
            "Outward-rounded real-parameter Arb blocks from hashed H2-H10 caches.",
        ),
        EntryRow(
            "co6m100ec_03_finite_ray",
            "interval_theorem_input",
            "ready_to_apply",
            "The exact-corridor dimensionless cover closes the finite mode ray.",
            sources["finite_ray"],
            "One collar block and 17999 rational mode blocks.",
        ),
        EntryRow(
            "co6m100ec_04_asymptotic_ray",
            "interval_analytic_theorem_input",
            "ready_to_apply",
            "The dimensionless normalized-H box closes the asymptotic mode ray.",
            sources["asymptotic_ray"] + "; " + exact["asymptotic_reserve"],
            "Uniform z in [0,10^-30] interval theorem.",
        ),
        EntryRow(
            "co6m100ec_05_global_curvature",
            "exact_theorem_composition",
            "ready_to_apply",
            "The three continuous certificates cover the complete first-summand tail.",
            exact["global_first_curvature"] + "; " + exact["coverage"],
            "Continuous first summand only.",
        ),
        EntryRow(
            "co6m100ec_06_first_discrete",
            "exact_analytic_theorem",
            "ready_to_apply",
            "The tent identity transfers the real-parameter theorem to the discrete first-summand curvature.",
            exact["tent_transfer"],
            "Exact tent integration and elementary logarithm bound.",
        ),
        EntryRow(
            "co6m100ec_07_full_kernel",
            "exact_perturbation_theorem",
            "ready_to_apply",
            "The proved higher-theta transfer closes the complete-kernel scalar ceiling.",
            sources["full_transfer"] + "; " + exact["full_ceiling"],
            "Complete Newman kernel at lambda=-100, k>=322.",
        ),
        EntryRow(
            "co6m100ec_08_tail",
            "exact_tail_theorem",
            "ready_to_apply",
            "The scalar ceiling proves every missing signed order-six tail sign.",
            exact["tail_sign"],
            "Analytic tail n>=317 only.",
        ),
        EntryRow(
            "co6m100ec_09_all_shift_entry",
            "exact_entry_theorem",
            "ready_to_apply",
            "The rigorous prefix and analytic tail prove signed contiguous order-six entry at every shift.",
            exact["all_shift_entry"] + "; " + exact["target_discharge"],
            "One heat parameter only; not forward propagation or PF-infinity.",
        ),
        EntryRow(
            "co6m100ec_10_forward_handoff",
            "theorem_handoff",
            "ready_to_apply",
            "The endpoint input required by the proved conditional cooperative-flow theorem is now available.",
            exact["next_handoff"],
            "Composition handoff; the forward theorem is recorded separately.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order6_m100_entry_certificate",
        "date": "2026-07-13",
        "status": "all-shift signed contiguous order-six entry theorem at lambda=-100",
        "proof_boundary": (
            "This artifact proves Q_(6,n)(-100)=-H_(6,n)(-100)>0 for every "
            "n>=0. It does not by itself prove forward order-six invariance, "
            "arbitrary-column order six, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_compound_order6_m100_prefix_certificate.md",
            "outputs/jensen_window_pf_compound_order6_m100_tail_curvature_reduction.md",
            "outputs/jensen_window_pf_compound_order6_first_summand_curvature_bridge.md",
            "outputs/jensen_window_pf_compound_order6_nested_curvature_compact_certificate.md",
            "outputs/jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.md",
            "outputs/jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate.md",
            "outputs/formal_core.md",
        ],
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order6_m100_entry_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order6_m100_entry_certificate.py",
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
            "discharged_local_targets": 1,
            "open_rows": 0,
        },
    }


def write_note(path: Path, artifact: dict) -> None:
    sources = artifact["source_contract"]
    exact = artifact["exact"]
    lines = [
        "# Jensen-Window PF Compound Order-Six Lambda=-100 Entry Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: all-shift signed contiguous order-six entry theorem at",
        "`lambda=-100`. This certificate is not a proof of forward order-six",
        "invariance, arbitrary-column order six, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order6_m100_entry_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order6_m100_entry_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order6_m100_entry_certificate.py",
        "```",
        "",
        "## Continuous First-Summand Theorem",
        "",
        "The third-nested stable curvature is covered on three overlapping ranges:",
        "",
        "```text",
        sources["compact"],
        sources["finite_ray"],
        sources["asymptotic_ray"],
        "```",
        "",
        f"The asymptotic scaled upper is `{sources['asymptotic_scaled_upper']}`.",
        "These ranges cover every real `t>=321`, hence",
        "",
        "```text",
        exact["global_first_curvature"],
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
        "prior tail reduction, so",
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
        "Thus `target_compound_order6_m100_entry` is discharged. The next",
        "composition is the already-derived uniform tail and cooperative flow",
        "from `lambda=-100` through `lambda=0`.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order6_m100_prefix_certificate.md",
        "outputs/jensen_window_pf_compound_order6_first_summand_curvature_bridge.md",
        "outputs/jensen_window_pf_compound_order6_nested_curvature_compact_certificate.md",
        "outputs/jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.md",
        "outputs/jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate.md",
        "outputs/jensen_window_pf_compound_order6_uniform_tail_flow_reduction.md",
        "outputs/signed_hankel_jensen_bridge_target.md",
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
        "wrote order-six lambda=-100 entry certificate: "
        f"{summary['rows']} rows, "
        f"{summary['continuous_curvature_theorems']} continuous curvature theorem, "
        f"{summary['analytic_tail_theorems']} analytic tail theorem, "
        f"{summary['all_shift_m100_entry_theorems']} all-shift entry theorem, "
        f"{summary['open_rows']} open rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
