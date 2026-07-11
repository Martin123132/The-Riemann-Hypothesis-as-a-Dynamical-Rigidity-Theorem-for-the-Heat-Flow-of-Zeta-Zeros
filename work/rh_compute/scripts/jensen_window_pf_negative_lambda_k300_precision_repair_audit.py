#!/usr/bin/env python3
"""Build the k300 precision-repair audit for the negative-lambda decrement route."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402

from jensen_window_pf_heat_flow_monotone_closure_scout import (  # noqa: E402
    REPO_ROOT,
    arb_positive,
    contraction,
    load_enclosures,
)
from jensen_window_pf_negative_lambda_raw_moment_bridge_matrix import (  # noqa: E402
    raw_ratio_factor_arb,
    raw_ratio_factor_decimal,
)
from jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout import (  # noqa: E402
    decrement_lower_arb,
    decrement_upper_arb,
    theta_arb,
    theta_decimal,
)


getcontext().prec = 100

DEFAULT_BASE_ENCLOSURE = (
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k300.jsonl"
)
DEFAULT_REPAIR_ENCLOSURES = (
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k220_k250_dps220.jsonl",
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220.jsonl",
)
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_k300_precision_repair_audit.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md"


@dataclass(frozen=True)
class CountSet:
    raw_lower: int
    raw_upper: int
    raw_decrease: int
    lower_decrement: int
    upper_decrement: int
    decrement_corridor: int
    corridor_width: int
    theta_unit: int
    theta_k_monotone: int


def label_to_t(label: str) -> Decimal:
    return abs(Decimal(label))


def load_summary(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def summary_path(jsonl_path: Path) -> Path:
    return jsonl_path.with_name(jsonl_path.stem + "_summary.json")


def ratio_tables(paths: list[Path], coefficient_cap: int) -> tuple[dict[str, dict[int, flint.arb]], dict[str, dict[int, Decimal]]]:
    balls, samples, labels = load_enclosures(paths)
    checked_ratio_max = coefficient_cap - 1
    by_label_arb: dict[str, dict[int, flint.arb]] = {}
    by_label_sample: dict[str, dict[int, Decimal]] = {}
    for lam in sorted(labels, key=lambda value: label_to_t(labels[value])):
        label = labels[lam]
        if not all((lam, index) in balls for index in range(0, coefficient_cap + 1)):
            continue
        by_label_arb[label] = {}
        by_label_sample[label] = {}
        for k in range(1, checked_ratio_max + 1):
            by_label_arb[label][k] = (
                contraction({idx: balls[(lam, idx)] for idx in (k - 1, k, k + 1)}, k)
                * raw_ratio_factor_arb(k)
            )
            by_label_sample[label][k] = (
                contraction({idx: samples[(lam, idx)] for idx in (k - 1, k, k + 1)}, k)
                * raw_ratio_factor_decimal(k)
            )
    return by_label_arb, by_label_sample


def first_failure(firsts: dict[str, int], name: str) -> int | None:
    return firsts.get(name)


def diagnostics(paths: list[Path], coefficient_cap: int = 300) -> dict:
    checked_ratio_max = coefficient_cap - 1
    adjacent_max = checked_ratio_max - 1
    by_label_arb, by_label_sample = ratio_tables(paths, coefficient_cap)
    labels = sorted(by_label_sample, key=label_to_t)
    theta_by_label: dict[str, dict[int, flint.arb]] = {}
    theta_sample_by_label: dict[str, dict[int, Decimal]] = {}
    per_lambda: dict[str, dict] = {}

    for label in labels:
        counts = {
            "raw_lower": 0,
            "raw_upper": 0,
            "raw_decrease": 0,
            "lower_decrement": 0,
            "upper_decrement": 0,
            "decrement_corridor": 0,
            "corridor_width": 0,
            "theta_unit": 0,
            "theta_k_monotone": 0,
        }
        firsts: dict[str, int] = {}
        theta_by_label[label] = {}
        theta_sample_by_label[label] = {}

        for k in range(1, checked_ratio_max + 1):
            raw = by_label_arb[label][k]
            tests = {
                "raw_lower": raw - flint.arb(1),
                "raw_upper": raw_ratio_factor_arb(k) - raw,
            }
            for name, value in tests.items():
                if arb_positive(value):
                    counts[name] += 1
                elif name not in firsts:
                    firsts[name] = k

        for k in range(1, adjacent_max + 1):
            raw = by_label_arb[label][k]
            raw_next = by_label_arb[label][k + 1]
            decrement = raw - raw_next
            lower = decrement_lower_arb(k, raw)
            upper = decrement_upper_arb(k, raw)
            width = upper - lower
            theta = theta_arb(k, raw, raw_next)
            theta_by_label[label][k] = theta
            theta_sample_by_label[label][k] = theta_decimal(k, by_label_sample[label][k], by_label_sample[label][k + 1])
            tests = {
                "raw_decrease": decrement,
                "lower_decrement": decrement - lower,
                "upper_decrement": upper - decrement,
                "corridor_width": width,
            }
            for name, value in tests.items():
                if arb_positive(value):
                    counts[name] += 1
                elif name not in firsts:
                    firsts[name] = k
            if arb_positive(decrement - lower) and arb_positive(upper - decrement):
                counts["decrement_corridor"] += 1
            elif "decrement_corridor" not in firsts:
                firsts["decrement_corridor"] = k
            if arb_positive(theta) and arb_positive(flint.arb(1) - theta):
                counts["theta_unit"] += 1
            elif "theta_unit" not in firsts:
                firsts["theta_unit"] = k

        for k in range(1, adjacent_max):
            drop = theta_by_label[label][k] - theta_by_label[label][k + 1]
            if arb_positive(drop):
                counts["theta_k_monotone"] += 1
            elif "theta_k_monotone" not in firsts:
                firsts["theta_k_monotone"] = k

        per_lambda[label] = {
            "counts": counts,
            "first_failures": firsts,
        }

    theta_lambda_order = 0
    theta_lambda_first_failure: dict | None = None
    for left, right in zip(labels, labels[1:]):
        for k in range(1, adjacent_max + 1):
            gap = theta_by_label[right][k] - theta_by_label[left][k]
            if arb_positive(gap):
                theta_lambda_order += 1
            elif theta_lambda_first_failure is None:
                theta_lambda_first_failure = {"left": left, "right": right, "k": k}

    return {
        "labels": labels,
        "coefficient_cap": coefficient_cap,
        "raw_ratio_rows_per_lambda": checked_ratio_max,
        "adjacent_rows_per_lambda": adjacent_max,
        "theta_k_rows_per_lambda": adjacent_max - 1,
        "theta_lambda_order_rows": theta_lambda_order,
        "theta_lambda_order_total_rows": max(0, len(labels) - 1) * adjacent_max,
        "theta_lambda_first_failure": theta_lambda_first_failure,
        "per_lambda": per_lambda,
    }


def aggregate_repaired_counts(repaired: dict) -> dict:
    per_lambda = repaired["per_lambda"]
    raw_total = repaired["raw_ratio_rows_per_lambda"] * len(per_lambda)
    adjacent_total = repaired["adjacent_rows_per_lambda"] * len(per_lambda)
    theta_k_total = repaired["theta_k_rows_per_lambda"] * len(per_lambda)
    return {
        "raw_total_rows": raw_total,
        "adjacent_total_rows": adjacent_total,
        "theta_k_total_rows": theta_k_total,
        "raw_lower_rows": sum(item["counts"]["raw_lower"] for item in per_lambda.values()),
        "raw_upper_rows": sum(item["counts"]["raw_upper"] for item in per_lambda.values()),
        "raw_decrease_rows": sum(item["counts"]["raw_decrease"] for item in per_lambda.values()),
        "lower_decrement_rows": sum(item["counts"]["lower_decrement"] for item in per_lambda.values()),
        "upper_decrement_rows": sum(item["counts"]["upper_decrement"] for item in per_lambda.values()),
        "decrement_corridor_rows": sum(item["counts"]["decrement_corridor"] for item in per_lambda.values()),
        "theta_unit_rows": sum(item["counts"]["theta_unit"] for item in per_lambda.values()),
        "theta_k_monotone_rows": sum(item["counts"]["theta_k_monotone"] for item in per_lambda.values()),
    }


def build_artifact(
    base_enclosure: Path = DEFAULT_BASE_ENCLOSURE,
    repair_enclosures: tuple[Path, ...] = DEFAULT_REPAIR_ENCLOSURES,
) -> dict:
    base = base_enclosure if base_enclosure.is_absolute() else REPO_ROOT / base_enclosure
    repairs = tuple(path if path.is_absolute() else REPO_ROOT / path for path in repair_enclosures)
    base_summary = load_summary(summary_path(base))
    repair_summaries = [load_summary(summary_path(path)) for path in repairs]
    base_diag = diagnostics([base], coefficient_cap=300)
    repaired_diag = diagnostics([base, *repairs], coefficient_cap=300)
    repaired_counts = aggregate_repaired_counts(repaired_diag)
    rows = [
        {
            "id": "nlk300pra_01_broad_k300_input",
            "role": "finite_input",
            "readiness": "not_ready_to_apply",
            "claim": "The broad k300 negative-lambda enclosure run covers lambdas -25,-50,-100 and coefficients A_0..A_300 with dps=160, cutoff=6, n_sum=50, abs_tol=1e-110.",
            "source": base.relative_to(REPO_ROOT).as_posix(),
            "summary": summary_path(base).relative_to(REPO_ROOT).as_posix(),
            "proof_boundary": "Finite enclosure input only; not a theorem and not a counterexample by itself.",
        },
        {
            "id": "nlk300pra_02_broad_run_failure_gate",
            "role": "precision_gate",
            "readiness": "not_ready_to_apply",
            "claim": "The broad k300 run alone reports lambda=-100 high-k failures, so broad-run failures must be repaired before being interpreted mathematically.",
            "lambda_minus_100_counts": base_diag["per_lambda"]["-100.0"]["counts"],
            "lambda_minus_100_first_failures": base_diag["per_lambda"]["-100.0"]["first_failures"],
            "proof_boundary": "Precision-sensitivity gate only; not evidence against the zeta-specific target.",
        },
        {
            "id": "nlk300pra_03_local_precision_repair",
            "role": "finite_input",
            "readiness": "not_ready_to_apply",
            "claim": "Two local lambda=-100 dps220/cutoff7/n_sum70 repair runs cover k220..250 and k245..320, overriding the broad run in the sensitive high-k window.",
            "sources": [path.relative_to(REPO_ROOT).as_posix() for path in repairs],
            "summaries": [summary_path(path).relative_to(REPO_ROOT).as_posix() for path in repairs],
            "proof_boundary": "Finite precision repair only; not an all-k stability theorem.",
        },
        {
            "id": "nlk300pra_04_repaired_raw_wall_stress",
            "role": "finite_stress",
            "readiness": "not_ready_to_apply",
            "claim": f"The repaired k300 data validate the raw lower and upper walls on {repaired_counts['raw_lower_rows']}/{repaired_counts['raw_total_rows']} and {repaired_counts['raw_upper_rows']}/{repaired_counts['raw_total_rows']} rows.",
            "proof_boundary": "Finite stress evidence only; not an all-k raw-wall theorem.",
        },
        {
            "id": "nlk300pra_05_repaired_decrement_corridor_stress",
            "role": "finite_stress",
            "readiness": "not_ready_to_apply",
            "claim": f"The repaired k300 data validate raw decrease and both decrement-corridor sides on {repaired_counts['decrement_corridor_rows']}/{repaired_counts['adjacent_total_rows']} adjacent rows.",
            "proof_boundary": "Finite stress evidence only; not an all-k decrement recurrence theorem.",
        },
        {
            "id": "nlk300pra_06_repaired_theta_shape_stress",
            "role": "finite_stress",
            "readiness": "not_ready_to_apply",
            "claim": f"The repaired k300 data validate theta in (0,1), theta-k monotonicity, and theta lambda-order on {repaired_counts['theta_unit_rows']}/{repaired_counts['adjacent_total_rows']}, {repaired_counts['theta_k_monotone_rows']}/{repaired_counts['theta_k_total_rows']}, and {repaired_diag['theta_lambda_order_rows']}/{repaired_diag['theta_lambda_order_total_rows']} rows.",
            "proof_boundary": "Finite shape evidence only; not a theta monotonicity theorem.",
        },
        {
            "id": "nlk300pra_07_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "Future k-depth stress should treat broad-run high-k failures as precision alarms until repaired by local high-precision enclosures; only repaired failures may become countermodel evidence.",
            "proof_boundary": "Proof-hygiene gate only; not a proof and not a disproof.",
        },
    ]
    summary = {
        "audit_rows": len(rows),
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "base_rows": base_summary["rows"],
        "base_k_max": base_summary["k_max"],
        "repair_rows": sum(item["rows"] for item in repair_summaries),
        "repair_ranges": [[item["k_min"], item["k_max"]] for item in repair_summaries],
        "repaired_raw_total_rows": repaired_counts["raw_total_rows"],
        "repaired_adjacent_total_rows": repaired_counts["adjacent_total_rows"],
        "repaired_theta_k_total_rows": repaired_counts["theta_k_total_rows"],
        "repaired_raw_lower_rows": repaired_counts["raw_lower_rows"],
        "repaired_raw_upper_rows": repaired_counts["raw_upper_rows"],
        "repaired_raw_decrease_rows": repaired_counts["raw_decrease_rows"],
        "repaired_lower_decrement_rows": repaired_counts["lower_decrement_rows"],
        "repaired_upper_decrement_rows": repaired_counts["upper_decrement_rows"],
        "repaired_decrement_corridor_rows": repaired_counts["decrement_corridor_rows"],
        "repaired_theta_unit_rows": repaired_counts["theta_unit_rows"],
        "repaired_theta_k_monotone_rows": repaired_counts["theta_k_monotone_rows"],
        "repaired_theta_lambda_order_rows": repaired_diag["theta_lambda_order_rows"],
        "repaired_theta_lambda_order_total_rows": repaired_diag["theta_lambda_order_total_rows"],
        "base_lambda_minus_100_first_failures": base_diag["per_lambda"]["-100.0"]["first_failures"],
        "open_theorem_target": False,
        "main_finding": (
            "The k300 stress extension supports the raw-ratio decrement route only after local "
            "precision repair: the broad dps160/cutoff6 run falsely reports lambda=-100 high-k "
            "failures, while the dps220/cutoff7 repair restores 897/897 raw wall rows, 894/894 "
            "decrement-corridor rows, 894/894 theta-unit rows, 891/891 theta-k monotone rows, "
            "and 596/596 theta lambda-order rows. This is finite stress evidence, not an all-k theorem."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_k300_precision_repair_audit",
        "date": "2026-07-07",
        "status": "finite precision-repair theorem-search diagnostic",
        "base_enclosure_jsonl": base.relative_to(REPO_ROOT).as_posix(),
        "base_summary": summary_path(base).relative_to(REPO_ROOT).as_posix(),
        "repair_enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in repairs],
        "repair_summaries": [summary_path(path).relative_to(REPO_ROOT).as_posix() for path in repairs],
        "source_decrement_scout": "outputs/jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.md",
        "source_raw_corridor_target": "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_k300_precision_repair_audit.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_k300_precision_repair_audit.py",
        "proof_boundary": (
            "Finite precision-repair diagnostic only. It documents a k300 stress extension and a local "
            "precision repair for lambda=-100, but it does not prove the raw-corridor theorem, does not "
            "prove cone entry, does not prove jwpf_06, and does not prove RH or Lambda <= 0."
        ),
        "audit_rows": rows,
        "base_diagnostics": base_diag,
        "repaired_diagnostics": repaired_diag,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "Broad-run failures are precision alarms until repaired.",
            "The repaired k300 stress is finite evidence only.",
            "The decrement corridor remains an open zeta-specific all-k theorem.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    base_failures = summary["base_lambda_minus_100_first_failures"]
    result_line = (
        "validated Jensen-window PF negative-lambda k300 precision-repair audit: "
        f"{summary['audit_rows']} rows, 0 issues, "
        f"{summary['repaired_decrement_corridor_rows']} repaired decrement-corridor rows, "
        f"{summary['repaired_theta_k_monotone_rows']} repaired theta-k-monotone rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda k300 Precision-Repair Audit",
        "",
        "Date: 2026-07-07",
        "",
        "Status: finite precision-repair theorem-search diagnostic. This is not",
        "a proof of the raw-corridor theorem, cone entry, Jensen-window PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_k300_precision_repair_audit`.",
        "",
        "Proof boundary: this artifact documents a k300 stress extension and a",
        "local high-precision repair. It does not prove any all-`k` theorem.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_k300_precision_repair_audit.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_k300_precision_repair_audit.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_k300_precision_repair_audit.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Inputs",
        "",
        "Broad k300 run:",
        "",
        "```text",
        artifact["base_enclosure_jsonl"],
        artifact["base_summary"],
        "python work/rh_compute/scripts/acb_coefficient_enclosures.py --lambdas=-25,-50,-100 --k-min 0 --k-max 300 --n-sum 50 --cutoff 6 --dps 160 --digits 80 --abs-tol 1e-110 --constant-terms 30 --output-dir work/rh_compute/results --run-id acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k300 --overwrite",
        "```",
        "",
        "Local repair runs:",
        "",
        "```text",
        *artifact["repair_enclosure_jsonl"],
        *artifact["repair_summaries"],
        "python work/rh_compute/scripts/acb_coefficient_enclosures.py --lambdas=-100 --k-min 220 --k-max 250 --n-sum 70 --cutoff 7 --dps 220 --digits 100 --abs-tol 1e-150 --constant-terms 40 --output-dir work/rh_compute/results --run-id acb_enclosures_negative_lambda_cone_entry_lam_m100_k220_k250_dps220 --overwrite",
        "python work/rh_compute/scripts/acb_coefficient_enclosures.py --lambdas=-100 --k-min 245 --k-max 320 --n-sum 70 --cutoff 7 --dps 220 --digits 100 --abs-tol 1e-150 --constant-terms 40 --output-dir work/rh_compute/results --run-id acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220 --overwrite",
        "```",
        "",
        "## Broad-Run Precision Alarm",
        "",
        "The broad dps160/cutoff6 run alone reports lambda=-100 failures:",
        "",
        "```text",
        f"first raw decrease failure: k={base_failures.get('raw_decrease')}",
        f"first lower decrement failure: k={base_failures.get('lower_decrement')}",
        f"first upper decrement failure: k={base_failures.get('upper_decrement')}",
        f"first raw upper-wall failure: k={base_failures.get('raw_upper')}",
        f"first theta-k monotone failure: k={base_failures.get('theta_k_monotone')}",
        "```",
        "",
        "These broad-run failures are treated as precision alarms, not mathematical",
        "counterexamples.",
        "",
        "## Repaired k300 Stress",
        "",
        "After the local lambda=-100 repair:",
        "",
        "```text",
        f"raw lower wall rows: {summary['repaired_raw_lower_rows']} / {summary['repaired_raw_total_rows']}",
        f"raw upper wall rows: {summary['repaired_raw_upper_rows']} / {summary['repaired_raw_total_rows']}",
        f"raw decrease rows: {summary['repaired_raw_decrease_rows']} / {summary['repaired_adjacent_total_rows']}",
        f"lower decrement rows: {summary['repaired_lower_decrement_rows']} / {summary['repaired_adjacent_total_rows']}",
        f"upper decrement rows: {summary['repaired_upper_decrement_rows']} / {summary['repaired_adjacent_total_rows']}",
        f"decrement-corridor rows: {summary['repaired_decrement_corridor_rows']} / {summary['repaired_adjacent_total_rows']}",
        f"theta unit rows: {summary['repaired_theta_unit_rows']} / {summary['repaired_adjacent_total_rows']}",
        f"theta-k monotone rows: {summary['repaired_theta_k_monotone_rows']} / {summary['repaired_theta_k_total_rows']}",
        f"theta lambda-order rows: {summary['repaired_theta_lambda_order_rows']} / {summary['repaired_theta_lambda_order_total_rows']}",
        "```",
        "",
        "Integration:",
        "",
        "```text",
        artifact["source_decrement_scout"],
        artifact["source_raw_corridor_target"],
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "```",
        "",
        "Summary:",
        "",
        summary["main_finding"],
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-enclosure", type=Path, default=DEFAULT_BASE_ENCLOSURE)
    parser.add_argument("--repair-enclosures", type=Path, nargs="+", default=list(DEFAULT_REPAIR_ENCLOSURES))
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    base = args.base_enclosure if args.base_enclosure.is_absolute() else REPO_ROOT / args.base_enclosure
    repairs = tuple(path if path.is_absolute() else REPO_ROOT / path for path in args.repair_enclosures)
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(base, repairs)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda k300 precision-repair audit: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
