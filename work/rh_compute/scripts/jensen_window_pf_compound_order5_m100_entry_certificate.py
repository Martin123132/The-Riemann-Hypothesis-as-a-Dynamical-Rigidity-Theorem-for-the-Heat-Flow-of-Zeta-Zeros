#!/usr/bin/env python3
"""Compose the all-shift contiguous order-five entry theorem at lambda=-100."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_compound_order5_m100_entry_certificate.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_compound_order5_m100_entry_certificate.md"
PREFIX_SOURCE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_compound_order5_m100_prefix_certificate.json"
TAIL_SOURCE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_compound_order5_m100_tail_curvature_reduction.json"
BRIDGE_SOURCE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_compound_order5_first_summand_curvature_bridge.json"
COMPACT_SOURCE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_compound_order5_nested_curvature_compact_certificate.json"
FINITE_SOURCE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.json"
ASYMPTOTIC_SOURCE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate.json"


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
    if prefix.get("summary", {}).get("positive_H5_rows") != 317:
        raise RuntimeError("order-five prefix source changed")
    if tail.get("exact", {}).get("sufficient_ceiling") != "C_n<=100/k^2 for every k=n+4>=321":
        raise RuntimeError("order-five scalar tail reduction changed")
    full_transfer = bridge.get("exact", {}).get("full_transfer")
    if full_transfer != "|C_n-C_n^(1)|<=E_(k-1)+2*E_k+E_(k+1)<=37/k^2, k=n+4>=321":
        raise RuntimeError("order-five full-kernel transfer changed")
    if compact.get("compact", {}).get("certified_t_range") != "320<=t<=V'(2)":
        raise RuntimeError("order-five compact curvature source changed")
    if finite.get("summary", {}).get("finite_ray_theorems") != 1:
        raise RuntimeError("order-five finite-ray source is not closed")
    if asymptotic.get("summary", {}).get("asymptotic_curvature_theorems") != 1:
        raise RuntimeError("order-five asymptotic source is not closed")
    return {
        "prefix": "H_(5,n)(-100)>0 for every 0<=n<=316",
        "compact": "q_1''(t)<=60/t^2 for 320<=t<=V'(2)",
        "finite_ray": "q_1''(t)<=60/t^2 for every mode 2<=u<=20",
        "asymptotic_ray": "q_1''(t)<=60/t^2 for every mode u>=20",
        "asymptotic_scaled_upper": asymptotic["dimensionless_interval"]["scaled_curvature_upper"],
        "full_transfer": full_transfer,
        "complete_ceiling": tail["exact"]["sufficient_ceiling"],
        "conditional_tail": tail["conclusions"]["conditional_tail"],
    }


def exact_composition(sources: dict) -> dict:
    return {
        "global_first_curvature": "q_1''(t)<=60/t^2 for every real t>=320",
        "coverage": "[320<=t<=V'(2)] union [2<=u<=20] union [u>=20] covers every t>=320",
        "tent_transfer": "C_n^(1)<=60*[-log(1-1/k^2)]<=60/(k^2-1)<63/k^2, k=n+4>=321",
        "full_ceiling": "C_n<=C_n^(1)+|C_n-C_n^(1)|<63/k^2+37/k^2=100/k^2, k>=321",
        "tail_sign": "C_n<=100/(n+4)^2 => J_n(-100)>0 and H_(5,n)(-100)>0, n>=317",
        "all_shift_entry": "H_(5,n)(-100)>0 for every integer n>=0",
        "target_discharge": "target_compound_order5_m100_entry is discharged",
        "next_handoff": "compose all-shift lambda=-100 entry with the completed uniform tail and cooperative flow",
        "asymptotic_reserve": f"u>=20 scaled upper {sources['asymptotic_scaled_upper']}<10<60",
    }


def build_artifact() -> dict:
    sources = validate_sources()
    exact = exact_composition(sources)
    rows = [
        EntryRow("co5m100ec_01_stable_coordinate", "exact_identity_input", "ready_to_apply", "The order-five sign is controlled by the stable scalar curvature coordinate.", "J_n>0 iff C_n<-4*log(x_(n+4))", "Previously proved exact factorization inside the positive lower cone."),
        EntryRow("co5m100ec_02_compact_curvature", "interval_theorem_input", "ready_to_apply", "The continuous first-summand curvature ceiling holds from the splice through mode two.", sources["compact"], "Thirty-six real-parameter Arb blocks from a hashed quadrature cache."),
        EntryRow("co5m100ec_03_finite_ray", "interval_theorem_input", "ready_to_apply", "The exact-corridor interval cover closes the finite mode ray.", sources["finite_ray"], "One collar block and 1850 exact-corridor blocks."),
        EntryRow("co5m100ec_04_asymptotic_ray", "interval_analytic_theorem_input", "ready_to_apply", "The dimensionless normalized-H box closes the asymptotic mode ray.", sources["asymptotic_ray"] + "; " + exact["asymptotic_reserve"], "Uniform z in [0,10^-30] interval theorem."),
        EntryRow("co5m100ec_05_global_curvature", "exact_theorem_composition", "ready_to_apply", "The three continuous certificates cover the complete first-summand tail.", exact["global_first_curvature"] + "; " + exact["coverage"], "Continuous first summand only."),
        EntryRow("co5m100ec_06_first_discrete", "exact_analytic_theorem", "ready_to_apply", "The tent identity transfers the real-parameter theorem to the discrete first-summand curvature.", exact["tent_transfer"], "Exact tent integration and elementary logarithm bound."),
        EntryRow("co5m100ec_07_full_kernel", "exact_perturbation_theorem", "ready_to_apply", "The proved higher-theta transfer closes the complete-kernel scalar ceiling.", sources["full_transfer"] + "; " + exact["full_ceiling"], "Complete Newman kernel at lambda=-100, k>=321."),
        EntryRow("co5m100ec_08_tail", "exact_tail_theorem", "ready_to_apply", "The scalar ceiling proves every missing order-five tail sign.", exact["tail_sign"], "Analytic tail n>=317 only."),
        EntryRow("co5m100ec_09_all_shift_entry", "exact_entry_theorem", "ready_to_apply", "The rigorous prefix and analytic tail prove contiguous order-five entry at every shift.", exact["all_shift_entry"] + "; " + exact["target_discharge"], "One heat parameter only; not forward propagation or PF-infinity."),
        EntryRow("co5m100ec_10_forward_handoff", "theorem_handoff", "ready_to_apply", "The endpoint input required by the proved conditional cooperative-flow theorem is now available.", exact["next_handoff"], "Composition handoff; the forward theorem is recorded separately."),
    ]
    return {
        "kind": "jensen_window_pf_compound_order5_m100_entry_certificate",
        "date": "2026-07-13",
        "status": "all-shift contiguous order-five entry theorem at lambda=-100",
        "proof_boundary": "This artifact proves H_(5,n)(-100)>0 for every n>=0. It does not by itself prove forward order-five invariance, arbitrary-column order five, PF-infinity, RH, or Lambda<=0.",
        "sources": [
            "outputs/jensen_window_pf_compound_order5_m100_prefix_certificate.md",
            "outputs/jensen_window_pf_compound_order5_m100_tail_curvature_reduction.md",
            "outputs/jensen_window_pf_compound_order5_first_summand_curvature_bridge.md",
            "outputs/jensen_window_pf_compound_order5_nested_curvature_compact_certificate.md",
            "outputs/jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.md",
            "outputs/jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate.md",
            "outputs/formal_core.md",
        ],
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order5_m100_entry_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order5_m100_entry_certificate.py",
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
        "# Jensen-Window PF Compound Order-Five Lambda=-100 Entry Certificate",
        "", "Date: 2026-07-13", "",
        "Status: all-shift contiguous order-five entry theorem at `lambda=-100`.",
        "This is not a proof of forward order-five invariance by itself, arbitrary-column order five, PF-infinity, RH, or `Lambda <= 0`.",
        "This is not by itself a proof of the forward or arbitrary-column theorem.",
        "", "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order5_m100_entry_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order5_m100_entry_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order5_m100_entry_certificate.py",
        "```", "", "## Continuous First-Summand Theorem", "",
        "The nested stable curvature is now covered on three overlapping ranges:", "", "```text",
        sources["compact"], sources["finite_ray"], sources["asymptotic_ray"], "```", "",
        f"The asymptotic scaled upper is `{sources['asymptotic_scaled_upper']}`. These ranges cover every real `t>=320`, hence", "", "```text",
        exact["global_first_curvature"], "```", "", "## Discrete And Full-Kernel Transfer", "", "```text",
        exact["tent_transfer"], sources["full_transfer"], exact["full_ceiling"], "```", "",
        "The last line is exactly the scalar ceiling isolated by the prior tail reduction, so", "", "```text",
        exact["tail_sign"], "```", "", "## Entry Theorem", "",
        "The 1024-bit prefix proves", "", "```text", sources["prefix"], "```", "",
        "Combining prefix and tail gives", "", "```text", exact["all_shift_entry"], "```", "",
        "Thus `target_compound_order5_m100_entry` is discharged. The next composition is the already-derived uniform tail and cooperative flow from `lambda=-100` through `lambda=0`.",
        "", "```text",
        "outputs/jensen_window_pf_compound_order5_m100_prefix_certificate.md",
        "outputs/jensen_window_pf_compound_order5_first_summand_curvature_bridge.md",
        "outputs/jensen_window_pf_compound_order5_nested_curvature_compact_certificate.md",
        "outputs/jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.md",
        "outputs/jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate.md",
        "outputs/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.md",
        "outputs/signed_hankel_jensen_bridge_target.md", "outputs/formal_core.md", "```",
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
    with args.out.open("w", encoding="utf-8") as handle:
        json.dump(artifact, handle, indent=2, sort_keys=True)
        handle.write("\n")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print("wrote order-five lambda=-100 entry certificate: "
          f"{summary['rows']} rows, {summary['continuous_curvature_theorems']} continuous curvature theorem, "
          f"{summary['analytic_tail_theorems']} analytic tail theorem, "
          f"{summary['all_shift_m100_entry_theorems']} all-shift entry theorem, {summary['open_rows']} open rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
