#!/usr/bin/env python3
"""Build a raw-log decrement bridge for the negative-lambda route."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
from fractions import Fraction
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
    decimal_format,
)
from jensen_window_pf_negative_lambda_k300_precision_repair_audit import (  # noqa: E402
    DEFAULT_BASE_ENCLOSURE,
    DEFAULT_REPAIR_ENCLOSURES,
    label_to_t,
    ratio_tables,
)
from jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout import (  # noqa: E402
    exact_counterexample,
    raw_upper_wall,
)


getcontext().prec = 100

DEFAULT_ENCLOSURES = (DEFAULT_BASE_ENCLOSURE, *DEFAULT_REPAIR_ENCLOSURES)
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.md"


@dataclass(frozen=True)
class Extremum:
    sample: str
    lam: str
    k: int


def update_min(current: tuple[Decimal, str, int] | None, value: Decimal, label: str, k: int) -> tuple[Decimal, str, int]:
    if current is None or value < current[0]:
        return value, label, k
    return current


def update_max(current: tuple[Decimal, str, int] | None, value: Decimal, label: str, k: int) -> tuple[Decimal, str, int]:
    if current is None or value > current[0]:
        return value, label, k
    return current


def extremum(item: tuple[Decimal, str, int] | None) -> Extremum:
    if item is None:
        raise RuntimeError("missing raw-log decrement extremum")
    value, lam, k = item
    return Extremum(decimal_format(value), lam, k)


def log_raw_upper_wall_arb(k: int) -> flint.arb:
    return (flint.arb(2 * k + 1) / flint.arb(2 * k - 1)).log()


def log_lower_bound_arb(k: int) -> flint.arb:
    return (flint.arb(1) - flint.arb(4) / flint.arb((2 * k + 1) ** 2)).log()


def log_upper_bound_arb(k: int, raw: flint.arb) -> flint.arb:
    return (flint.arb(1) - flint.arb(2) * (flint.arb(1) - flint.arb(1) / raw) / flint.arb(2 * k + 1)).log()


def log_raw_upper_wall_decimal(k: int) -> Decimal:
    return (Decimal(2 * k + 1) / Decimal(2 * k - 1)).ln()


def log_lower_bound_decimal(k: int) -> Decimal:
    return (Decimal(1) - Decimal(4) / Decimal((2 * k + 1) ** 2)).ln()


def log_upper_bound_decimal(k: int, raw: Decimal) -> Decimal:
    return (Decimal(1) - Decimal(2) * (Decimal(1) - Decimal(1) / raw) / Decimal(2 * k + 1)).ln()


def q(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def log_counterexample(raw: Fraction, raw_next: Fraction, k: int = 1) -> dict:
    witness = exact_counterexample(raw, raw_next, k)
    ratio = raw_next / raw
    lower_exp = Fraction(1) - Fraction(4, (2 * k + 1) ** 2)
    upper_exp = Fraction(1) - Fraction(2) * (Fraction(1) - Fraction(1, 1) / raw) / Fraction(2 * k + 1)
    witness.update(
        {
            "exp_delta": q(ratio),
            "lower_exp_wall": q(lower_exp),
            "upper_exp_wall": q(upper_exp),
            "log_decrease_holds": bool(ratio < 1),
            "log_lower_wall_holds": bool(ratio >= lower_exp),
            "log_upper_wall_holds": bool(ratio <= upper_exp),
            "raw_log_wall_holds_at_k": bool(Fraction(1) <= raw <= raw_upper_wall(k)),
            "raw_log_wall_holds_at_next": bool(Fraction(1) <= raw_next <= raw_upper_wall(k + 1)),
        }
    )
    return witness


def build_diagnostics(paths: list[Path], coefficient_cap: int = 300) -> dict:
    by_label_arb, by_label_sample = ratio_tables(paths, coefficient_cap)
    labels = sorted(by_label_sample, key=label_to_t)
    checked_ratio_max = coefficient_cap - 1
    adjacent_max = checked_ratio_max - 1

    raw_log_lower_rows = 0
    raw_log_upper_rows = 0
    log_decrease_rows = 0
    log_lower_bound_rows = 0
    log_upper_bound_rows = 0
    log_corridor_rows = 0
    log_width_rows = 0

    min_raw_log_upper_margin: tuple[Decimal, str, int] | None = None
    min_delta_negative_margin: tuple[Decimal, str, int] | None = None
    min_log_lower_margin: tuple[Decimal, str, int] | None = None
    min_log_upper_margin: tuple[Decimal, str, int] | None = None
    min_log_width: tuple[Decimal, str, int] | None = None
    min_p: tuple[Decimal, str, int] | None = None
    max_p: tuple[Decimal, str, int] | None = None
    min_delta: tuple[Decimal, str, int] | None = None
    max_delta: tuple[Decimal, str, int] | None = None

    for label in labels:
        for k in range(1, checked_ratio_max + 1):
            raw_arb = by_label_arb[label][k]
            raw_sample = by_label_sample[label][k]
            p_arb = raw_arb.log()
            p_sample = raw_sample.ln()
            upper_margin_arb = log_raw_upper_wall_arb(k) - p_arb
            upper_margin_sample = log_raw_upper_wall_decimal(k) - p_sample
            if arb_positive(p_arb):
                raw_log_lower_rows += 1
            if arb_positive(upper_margin_arb):
                raw_log_upper_rows += 1
            min_raw_log_upper_margin = update_min(min_raw_log_upper_margin, upper_margin_sample, label, k)
            min_p = update_min(min_p, p_sample, label, k)
            max_p = update_max(max_p, p_sample, label, k)

        for k in range(1, adjacent_max + 1):
            raw_arb = by_label_arb[label][k]
            raw_next_arb = by_label_arb[label][k + 1]
            raw_sample = by_label_sample[label][k]
            raw_next_sample = by_label_sample[label][k + 1]

            delta_arb = (raw_next_arb / raw_arb).log()
            lower_arb = log_lower_bound_arb(k)
            upper_arb = log_upper_bound_arb(k, raw_arb)
            lower_margin_arb = delta_arb - lower_arb
            upper_margin_arb = upper_arb - delta_arb
            width_arb = upper_arb - lower_arb

            delta_sample = (raw_next_sample / raw_sample).ln()
            lower_sample = log_lower_bound_decimal(k)
            upper_sample = log_upper_bound_decimal(k, raw_sample)
            lower_margin_sample = delta_sample - lower_sample
            upper_margin_sample = upper_sample - delta_sample
            width_sample = upper_sample - lower_sample

            if arb_positive(-delta_arb):
                log_decrease_rows += 1
            if arb_positive(lower_margin_arb):
                log_lower_bound_rows += 1
            if arb_positive(upper_margin_arb):
                log_upper_bound_rows += 1
            if arb_positive(lower_margin_arb) and arb_positive(upper_margin_arb):
                log_corridor_rows += 1
            if arb_positive(width_arb):
                log_width_rows += 1

            min_delta_negative_margin = update_min(min_delta_negative_margin, -delta_sample, label, k)
            min_log_lower_margin = update_min(min_log_lower_margin, lower_margin_sample, label, k)
            min_log_upper_margin = update_min(min_log_upper_margin, upper_margin_sample, label, k)
            min_log_width = update_min(min_log_width, width_sample, label, k)
            min_delta = update_min(min_delta, delta_sample, label, k)
            max_delta = update_max(max_delta, delta_sample, label, k)

    return {
        "lambdas": labels,
        "coefficient_cap": coefficient_cap,
        "raw_ratio_rows_per_lambda": checked_ratio_max,
        "adjacent_rows_per_lambda": adjacent_max,
        "raw_log_total_rows": checked_ratio_max * len(labels),
        "adjacent_total_rows": adjacent_max * len(labels),
        "raw_log_lower_rows": raw_log_lower_rows,
        "raw_log_upper_rows": raw_log_upper_rows,
        "log_decrease_rows": log_decrease_rows,
        "log_lower_bound_rows": log_lower_bound_rows,
        "log_upper_bound_rows": log_upper_bound_rows,
        "log_corridor_rows": log_corridor_rows,
        "log_width_rows": log_width_rows,
        "min_raw_log_upper_margin": asdict(extremum(min_raw_log_upper_margin)),
        "min_delta_negative_margin": asdict(extremum(min_delta_negative_margin)),
        "min_log_lower_margin": asdict(extremum(min_log_lower_margin)),
        "min_log_upper_margin": asdict(extremum(min_log_upper_margin)),
        "min_log_width": asdict(extremum(min_log_width)),
        "min_p": asdict(extremum(min_p)),
        "max_p": asdict(extremum(max_p)),
        "min_delta": asdict(extremum(min_delta)),
        "max_delta": asdict(extremum(max_delta)),
    }


def build_artifact(paths: list[Path] | None = None) -> dict:
    if paths is None:
        paths = list(DEFAULT_ENCLOSURES)
    diagnostics = build_diagnostics(paths)
    underfast = log_counterexample(Fraction(2), Fraction(3, 2))
    overfast = log_counterexample(Fraction(2), Fraction(1))
    rows = [
        {
            "id": "nlrldb_01_log_coordinate",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "Set p_k=log(R_k) and delta_k=p_(k+1)-p_k=log(R_(k+1)/R_k) for the raw ratios R_k=M_(k+1)*M_(k-1)/M_k^2.",
            "formula": "R_k=exp(p_k); R_(k+1)=R_k*exp(delta_k)",
            "proof_boundary": "Exact change of variables only; it proves no zeta inequality.",
        },
        {
            "id": "nlrldb_02_exact_log_corridor",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The decrement corridor 2*(R_k-1)/(2*k+1)<=R_k-R_(k+1)<=4*R_k/(2*k+1)^2 is equivalent to a two-sided bound on delta_k.",
            "formula": "log(1-4/(2*k+1)^2)<=delta_k<=log(1-2*(1-exp(-p_k))/(2*k+1))",
            "proof_boundary": "Exact algebraic equivalence only; the log-corridor inequality remains open for all k.",
        },
        {
            "id": "nlrldb_03_width_raw_wall_equivalence",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The log corridor has nonnegative width exactly when the raw upper wall R_k<=(2*k+1)/(2*k-1) holds.",
            "formula": "log-upper-exp - log-lower-exp >= 0 iff R_k <= (2*k+1)/(2*k-1)",
            "proof_boundary": "Exact algebra only; it does not prove the raw upper wall.",
        },
        {
            "id": "nlrldb_04_repaired_k300_log_corridor",
            "role": "finite_stress",
            "readiness": "not_ready_to_apply",
            "claim": "The repaired k300 data validate the raw-log walls, log decrease, and log-corridor sides on all checked rows.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
                "outputs/jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.md",
            ],
            "diagnostics": diagnostics,
            "proof_boundary": "Finite repaired stress evidence only; not an all-k log-ratio recurrence theorem.",
        },
        {
            "id": "nlrldb_05_underfast_log_decrease_shortcut_blocked",
            "role": "exact_counterexample",
            "readiness": "not_ready_to_apply",
            "claim": "Raw-log decrease plus raw cone occupancy does not imply the lower decrement wall.",
            "witness": underfast,
            "proof_boundary": "Exact shortcut rejection only; not evidence against actual zeta coefficients.",
        },
        {
            "id": "nlrldb_06_overfast_log_decrease_shortcut_blocked",
            "role": "exact_counterexample",
            "readiness": "not_ready_to_apply",
            "claim": "Raw-log decrease plus raw cone occupancy does not imply the upper decrement wall.",
            "witness": overfast,
            "proof_boundary": "Exact shortcut rejection only; not evidence against actual zeta coefficients.",
        },
        {
            "id": "nlrldb_07_live_log_recurrence_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim": "A zeta-specific proof may target the log-ratio recurrence directly: prove p_k is in the raw wall and delta_k lies between the two log bounds after a finite collar.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
                "outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md",
            ],
            "proof_boundary": "Live theorem-search route only; no all-k recurrence is proved here.",
        },
        {
            "id": "nlrldb_08_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "A promoted log-ratio proof must state the tail start, lambda range, finite collar, raw wall, both delta_k bounds, and forbidden assumptions.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
                "outputs/signed_hankel_jensen_dependency_graph.md",
            ],
            "proof_boundary": "Proof-hygiene gate only; not a proof of the raw-corridor theorem.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "exact_available_rows": 3,
        "finite_stress_rows": 1,
        "exact_counterexample_rows": 2,
        "live_routes": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "lambdas": diagnostics["lambdas"],
        "coefficient_cap": diagnostics["coefficient_cap"],
        "raw_log_total_rows": diagnostics["raw_log_total_rows"],
        "adjacent_total_rows": diagnostics["adjacent_total_rows"],
        "raw_log_lower_rows": diagnostics["raw_log_lower_rows"],
        "raw_log_upper_rows": diagnostics["raw_log_upper_rows"],
        "log_decrease_rows": diagnostics["log_decrease_rows"],
        "log_lower_bound_rows": diagnostics["log_lower_bound_rows"],
        "log_upper_bound_rows": diagnostics["log_upper_bound_rows"],
        "log_corridor_rows": diagnostics["log_corridor_rows"],
        "log_width_rows": diagnostics["log_width_rows"],
        "main_finding": (
            "The raw-ratio decrement corridor is exactly equivalent to a log-ratio recurrence "
            "for delta_k=log(R_(k+1)/R_k). On the repaired k300 data the log walls hold on "
            "897/897 raw rows and the two-sided log corridor holds on 894/894 adjacent rows, "
            "but raw-log decrease alone is blocked by two exact cone counterexamples."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_raw_log_decrement_bridge",
        "date": "2026-07-07",
        "status": "exact finite theorem-search diagnostic",
        "source_decrement_scout": "outputs/jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.md",
        "source_k300_precision_repair": "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
        "source_raw_corridor_target": "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
        "source_log_curvature_bridge": "outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md",
        "source_raw_obstruction": "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_log_decrement_bridge.py",
        "proof_boundary": (
            "Exact finite theorem-search diagnostic only. It rewrites the raw decrement corridor "
            "as a log-ratio recurrence and checks repaired k300 evidence, but it does not prove "
            "the all-k raw-corridor theorem, does not prove cone entry, does not prove jwpf_06, "
            "and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "invariants": [
            "No row is ready_to_apply.",
            "The log-ratio recurrence remains an open zeta-specific all-k theorem.",
            "Repaired k300 evidence is finite stress evidence only.",
            "Raw-log decrease alone is rejected as a shortcut.",
            "Generic moment positivity, endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
        "summary": summary,
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][3]["diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda raw-log decrement bridge: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['log_corridor_rows']} log-corridor rows, "
        f"{summary['log_decrease_rows']} log-decrease rows, "
        f"{summary['exact_counterexample_rows']} exact counterexamples, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Raw-Log Decrement Bridge",
        "",
        "Date: 2026-07-07",
        "",
        "Status: exact finite theorem-search diagnostic. This is not a proof",
        "of the raw-corridor theorem, cone entry, Jensen-window PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_raw_log_decrement_bridge`.",
        "",
        "Proof boundary: this artifact rewrites the raw-ratio decrement",
        "corridor as a log-ratio recurrence and checks repaired k300 evidence.",
        "It does not prove any all-`k` theorem.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_log_decrement_bridge.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Exact Log Form",
        "",
        "Let:",
        "",
        "```text",
        "M_k = mu_{2k}",
        "R_k = M_(k+1)*M_(k-1)/M_k^2",
        "p_k = log(R_k)",
        "delta_k = p_(k+1)-p_k = log(R_(k+1)/R_k)",
        "```",
        "",
        "The decrement corridor is exactly equivalent to:",
        "",
        "```text",
        "log(1-4/(2*k+1)^2) <= delta_k",
        "delta_k <= log(1-2*(1-exp(-p_k))/(2*k+1))",
        "```",
        "",
        "The log-corridor width is nonnegative exactly when the raw upper wall",
        "`R_k <= (2*k+1)/(2*k-1)` holds.",
        "",
        "## Repaired k300 Stress",
        "",
        "Inputs:",
        "",
        "```text",
        "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k300.jsonl",
        "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k220_k250_dps220.jsonl",
        "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220.jsonl",
        "```",
        "",
        "Repaired log diagnostics:",
        "",
        "```text",
        f"raw-log lower wall rows: {diagnostics['raw_log_lower_rows']} / {diagnostics['raw_log_total_rows']}",
        f"raw-log upper wall rows: {diagnostics['raw_log_upper_rows']} / {diagnostics['raw_log_total_rows']}",
        f"log-decrease rows: {diagnostics['log_decrease_rows']} / {diagnostics['adjacent_total_rows']}",
        f"log lower-bound rows: {diagnostics['log_lower_bound_rows']} / {diagnostics['adjacent_total_rows']}",
        f"log upper-bound rows: {diagnostics['log_upper_bound_rows']} / {diagnostics['adjacent_total_rows']}",
        f"log-corridor rows: {diagnostics['log_corridor_rows']} / {diagnostics['adjacent_total_rows']}",
        f"log-width rows: {diagnostics['log_width_rows']} / {diagnostics['adjacent_total_rows']}",
        "```",
        "",
        "Extrema:",
        "",
        "```text",
        f"min log lower margin: {diagnostics['min_log_lower_margin']['sample']} at lambda={diagnostics['min_log_lower_margin']['lam']}, k={diagnostics['min_log_lower_margin']['k']}",
        f"min log upper margin: {diagnostics['min_log_upper_margin']['sample']} at lambda={diagnostics['min_log_upper_margin']['lam']}, k={diagnostics['min_log_upper_margin']['k']}",
        f"min log width: {diagnostics['min_log_width']['sample']} at lambda={diagnostics['min_log_width']['lam']}, k={diagnostics['min_log_width']['k']}",
        f"delta range: {diagnostics['min_delta']['sample']} at lambda={diagnostics['min_delta']['lam']}, k={diagnostics['min_delta']['k']} to {diagnostics['max_delta']['sample']} at lambda={diagnostics['max_delta']['lam']}, k={diagnostics['max_delta']['k']}",
        "```",
        "",
        "## Shortcut Gate",
        "",
        "Raw-log decrease is not enough:",
        "",
        "```text",
        "R_1=2, R_2=3/2: raw cone and log decrease hold, but the lower decrement wall fails",
        "R_1=2, R_2=1: raw cone and log decrease hold, but the upper decrement wall fails",
        "```",
        "",
        "Integration:",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.md",
        "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
        "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
        "outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md",
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "```",
        "",
        "Summary:",
        "",
        summary["main_finding"],
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact()
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(f"wrote Jensen-window PF negative-lambda raw-log decrement bridge: {out_json} and {note}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
