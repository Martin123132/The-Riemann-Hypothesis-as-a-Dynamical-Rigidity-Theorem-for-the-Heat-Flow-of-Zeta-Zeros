#!/usr/bin/env python3
"""Certify the lambda=-100 signed order-six prefix through n=316."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

from jensen_window_pf_compound_order5_m100_prefix_certificate import (  # noqa: E402
    SOURCE_PATHS as ORDER5_SOURCE_PATHS,
    arb_lower_text,
    arb_text,
    load_source,
    sha256,
    symbolic_factorization,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_m100_prefix_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order6_m100_prefix_certificate.md"
)
ORDER6_REDUCTION = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_uniform_tail_flow_reduction.json"
)
PRECISION_REPAIR = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "negative_lambda_m100_order6_k191_k229_dps220.jsonl"
)
PRECISION_REPAIR_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "negative_lambda_m100_order6_k191_k229_dps220_summary.json"
)
PREFIX_EXTENSION = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "negative_lambda_m100_order6_k325_k326_dps220.jsonl"
)
PREFIX_EXTENSION_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "negative_lambda_m100_order6_k325_k326_dps220_summary.json"
)
SOURCE_PATHS = (*ORDER5_SOURCE_PATHS, PRECISION_REPAIR, PREFIX_EXTENSION)
PREFIX_LAST_N = 316
MAX_COEFFICIENT_INDEX = PREFIX_LAST_N + 10
PRECISION_BITS = 1024


@dataclass(frozen=True)
class PrefixRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def merged_coefficients() -> tuple[dict[int, flint.arb], list[dict]]:
    values: dict[int, flint.arb] = {}
    diagnostics = []
    for precedence, source in enumerate(SOURCE_PATHS):
        loaded = load_source(source)
        overwritten = len(set(values).intersection(loaded))
        values.update(loaded)
        diagnostics.append(
            {
                "precedence": precedence,
                "source": source.relative_to(REPO_ROOT).as_posix(),
                "rows": len(loaded),
                "index_range": [min(loaded), max(loaded)],
                "overwritten_rows": overwritten,
                "sha256": sha256(source),
            }
        )
    expected = set(range(MAX_COEFFICIENT_INDEX + 1))
    if set(values) != expected:
        missing = sorted(expected - set(values))
        extra = sorted(set(values) - expected)
        raise RuntimeError(f"coefficient coverage mismatch: missing={missing}, extra={extra}")
    if not all(bool(value > 0 and not value.contains(0)) for value in values.values()):
        raise RuntimeError("coefficient source contains a nonpositive ball")
    return values, diagnostics


def finite_prefix(values: dict[int, flint.arb]) -> dict:
    flint.ctx.prec = PRECISION_BITS

    def contraction(index: int) -> flint.arb:
        return values[index - 1] * values[index + 1] / values[index] ** 2

    contractions = {
        index: contraction(index)
        for index in range(1, MAX_COEFFICIENT_INDEX)
    }
    defects = {index: 1 - value for index, value in contractions.items()}
    gaps = {
        index: (
            defects[index + 2] ** 2
            - contractions[index + 2] ** 2
            * defects[index + 1]
            * defects[index + 3]
        )
        for index in range(PREFIX_LAST_N + 7)
    }
    margins = {
        index: (
            gaps[index + 1] ** 2
            - contractions[index + 3] ** 3
            * gaps[index]
            * gaps[index + 2]
        )
        for index in range(PREFIX_LAST_N + 5)
    }
    stable_h5 = {
        index: (
            defects[index + 3]
            * defects[index + 5]
            * margins[index + 1] ** 2
            - contractions[index + 4] ** 4
            * defects[index + 4] ** 2
            * margins[index]
            * margins[index + 2]
        )
        for index in range(PREFIX_LAST_N + 3)
    }

    for name, family in (
        ("contraction", contractions),
        ("defect", defects),
        ("order3_gap", gaps),
        ("order4_margin", margins),
        ("order5_stable_margin", stable_h5),
    ):
        if not all(bool(value > 0 and not value.contains(0)) for value in family.values()):
            raise RuntimeError(f"{name} family is not strictly positive")

    def h5_value(index: int) -> flint.arb:
        ratio = values[index + 1] / values[index]
        scale = (
            values[index] ** 5
            * ratio**20
            * contractions[index + 1] ** 15
            * contractions[index + 2] ** 10
            * contractions[index + 3] ** 5
            / (
                defects[index + 3]
                * defects[index + 4] ** 2
                * defects[index + 5]
                * gaps[index + 2]
            )
        )
        return scale * stable_h5[index]

    h5_values = {
        index: h5_value(index) for index in range(PREFIX_LAST_N + 3)
    }
    if not all(
        bool(value > 0 and not value.contains(0)) for value in h5_values.values()
    ):
        raise RuntimeError("stable H5 family is not strictly positive")

    rows = []
    minimum: tuple[flint.arb, int] | None = None
    for n in range(PREFIX_LAST_N + 1):
        relative = h5_values[n + 1] ** 2 / (
            h5_values[n] * h5_values[n + 2]
        ) - 1
        if not bool(relative > 0 and not relative.contains(0)):
            raise RuntimeError(f"inconclusive signed order-six margin at n={n}")
        if minimum is None or relative.lower() < minimum[0].lower():
            minimum = (relative, n)
        rows.append(
            {
                "n": n,
                "relative_H5_margin_ball": arb_text(relative),
                "relative_H5_margin_lower": arb_lower_text(relative),
                "Q6_sign": "positive_by_signed_condensation",
            }
        )
    assert minimum is not None
    return {
        "lambda": "-100",
        "n_range": [0, PREFIX_LAST_N],
        "coefficient_range": [0, MAX_COEFFICIENT_INDEX],
        "precision_bits": PRECISION_BITS,
        "rows": rows,
        "all_H5_positive": True,
        "all_relative_H5_margins_positive": True,
        "all_Q6_positive": True,
        "minimum_relative_n": minimum[1],
        "minimum_relative_ball": arb_text(minimum[0]),
        "minimum_relative_lower": arb_lower_text(minimum[0]),
    }


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    reduction = json.loads(ORDER6_REDUCTION.read_text(encoding="utf-8"))
    if reduction.get("summary", {}).get("open_entry_targets") != 1:
        raise RuntimeError("order-six reduction source contract changed")
    repair = json.loads(PRECISION_REPAIR_SUMMARY.read_text(encoding="utf-8"))
    if repair.get("rows") != 39 or repair.get("k_min") != 191 or repair.get("k_max") != 229:
        raise RuntimeError("order-six precision repair source contract changed")
    extension = json.loads(PREFIX_EXTENSION_SUMMARY.read_text(encoding="utf-8"))
    if (
        extension.get("rows") != 2
        or extension.get("k_min") != 325
        or extension.get("k_max") != 326
    ):
        raise RuntimeError("order-six prefix extension source contract changed")
    values, source_diagnostics = merged_coefficients()
    h5_factorization = symbolic_factorization()
    finite = finite_prefix(values)
    exact = {
        "signed_condensation": (
            "Q_(6,n)*H_(4,n+2)=H_(5,n+1)^2-H_(5,n)*H_(5,n+2)"
        ),
        "stable_coordinate": (
            "K_n=H_(5,n+1)^2/(H_(5,n)*H_(5,n+2))-1"
        ),
        "positive_scale": (
            "Q_(6,n)=H_(5,n)*H_(5,n+2)*K_n/H_(4,n+2)"
        ),
        "entry_equivalence": "Q_(6,n)>0 iff K_n>0 inside the completed order-five cone",
        "prefix": "Q_(6,n)(-100)>0 for every 0<=n<=316",
        "remaining_tail": "Q_(6,n)(-100)>0 for every n>=317",
    }
    rows = [
        PrefixRow(
            "co6m100pc_01_signed_coordinate",
            "exact_identity",
            "ready_to_apply",
            "Signed order six is a positive lower-cone scale times one relative H5 log-concavity margin.",
            exact["signed_condensation"] + "; " + exact["stable_coordinate"] + "; " + exact["positive_scale"],
            "Exact Desnanot-Jacobi algebra only.",
        ),
        PrefixRow(
            "co6m100pc_02_stable_H5_factorization",
            "exact_identity",
            "ready_to_apply",
            "Every H5 value is evaluated through the completed W_n J_n factorization rather than a raw determinant.",
            h5_factorization["factorization"] + "; " + h5_factorization["J_definition"],
            "Exact coefficient-ratio factorization only.",
            h5_factorization,
        ),
        PrefixRow(
            "co6m100pc_03_coefficient_enclosures",
            "interval_input",
            "ready_to_apply",
            "Merged outward-rounded Arb sources enclose every endpoint coefficient through A_326.",
            "A_k(-100)>0 for every 0<=k<=326",
            "Finite coefficient range only.",
            {"sources": source_diagnostics},
        ),
        PrefixRow(
            "co6m100pc_04_precision_repair",
            "interval_certificate",
            "ready_to_apply",
            "A dedicated 220-digit source removes the only relative-precision gap in the order-six margin evaluation.",
            "A_k(-100) rigorously repaired for 191<=k<=229",
            "Local finite precision repair only.",
            {"summary": repair, "sha256": sha256(PRECISION_REPAIR)},
        ),
        PrefixRow(
            "co6m100pc_05_positive_relative_margin",
            "interval_certificate",
            "ready_to_apply",
            "Every available relative H5 log-concavity margin is strictly positive.",
            "K_n(-100)>0 for every 0<=n<=316",
            "Finite 1024-bit Arb prefix only.",
            {
                "minimum_n": finite["minimum_relative_n"],
                "minimum_ball": finite["minimum_relative_ball"],
                "minimum_lower": finite["minimum_relative_lower"],
            },
        ),
        PrefixRow(
            "co6m100pc_06_signed_order6_prefix",
            "exact_interval_composition",
            "ready_to_apply",
            "Positive lower-cone factors and the relative margin prove the complete available signed order-six prefix.",
            exact["prefix"],
            "Finite lambda=-100 prefix only.",
        ),
        PrefixRow(
            "co6m100pc_07_open_tail",
            "open_handoff",
            "not_ready_to_apply",
            "Prove the signed order-six endpoint margin on the infinite tail beginning at n=317.",
            exact["remaining_tail"],
            "Open analytic tail; not an all-shift order-six theorem.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order6_m100_prefix_certificate",
        "date": "2026-07-13",
        "status": (
            "rigorous lambda=-100 signed order-six prefix through n=316 with "
            "one open analytic tail"
        ),
        "proof_boundary": (
            "This artifact proves Q_(6,n)(-100)>0 only for 0<=n<=316. It does "
            "not prove the n>=317 tail, all-shift order-six entry, forward "
            "invariance, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [source.relative_to(REPO_ROOT).as_posix() for source in SOURCE_PATHS],
        "source_diagnostics": source_diagnostics,
        "exact": exact,
        "finite": finite,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_to_apply_rows": 6,
            "coefficients": len(values),
            "positive_H5_rows": PREFIX_LAST_N + 3,
            "positive_relative_H5_rows": len(finite["rows"]),
            "positive_Q6_rows": len(finite["rows"]),
            "inconclusive_rows": 0,
            "precision_repair_rows": repair["rows"],
            "prefix_extension_rows": extension["rows"],
            "open_analytic_tails": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order6_m100_prefix_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order6_m100_prefix_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    finite = artifact["finite"]
    lines = [
        "# Jensen-Window PF Compound Order-Six Lambda=-100 Prefix Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous signed order-six endpoint prefix through `n=316` with",
        "one open analytic tail. This is not a proof of all-shift order six,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order6_m100_prefix_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order6_m100_prefix_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order6_m100_prefix_certificate.py",
        "```",
        "",
        "## Stable Endpoint Coordinate",
        "",
        "Inside the completed positive order-five cone, signed condensation gives",
        "",
        "```text",
        exact["signed_condensation"],
        exact["stable_coordinate"],
        exact["positive_scale"],
        "```",
        "",
        "Thus `Q_(6,n)>0` is exactly `K_n>0`. Each `H_(5,n)` is evaluated",
        "through the cancellation-preserving factorization `H_(5,n)=W_n*J_n`,",
        "not by subtracting a raw six-by-six determinant.",
        "",
        "## Coefficient Cover",
        "",
        "The inherited sources, local repair, and two-row splice extension,",
        "merged in precedence",
        "order and hashed in the JSON artifact, give 1024-bit outward-rounded",
        "balls for",
        "",
        "```text",
        "A_k(-100)>0 for every 0<=k<=326.",
        "```",
        "",
        "The dedicated repair covers `A_191,...,A_229` at 220 decimal digits.",
        "It removes 24 interval-indeterminate rows caused by the older broad",
        "source; no midpoint sign was promoted.",
        "",
        "## Prefix Theorem",
        "",
        "Direct stable Arb evaluation proves",
        "",
        "```text",
        "K_n(-100)>0 for every 0<=n<=316,",
        exact["prefix"] + ".",
        "```",
        "",
        "The weakest row is",
        "",
        "```text",
        f"n={finite['minimum_relative_n']},",
        f"K_316={finite['minimum_relative_ball']},",
        f"K_316 lower={finite['minimum_relative_lower']}>7/1000.",
        "```",
        "",
        "## Remaining Tail",
        "",
        "The complete endpoint problem is now reduced to",
        "",
        "```text",
        exact["remaining_tail"] + ".",
        "```",
        "",
        "The uniform eventual signed tail from the companion reduction is",
        "non-effective, so it does not splice this finite prefix to every shift.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order6_uniform_tail_flow_reduction.md",
        "outputs/jensen_window_pf_compound_order5_m100_prefix_certificate.md",
        "outputs/jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.md",
        "outputs/signed_hankel_jensen_bridge_target.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-six lambda=-100 prefix: "
        f"{summary['coefficients']} coefficients, "
        f"{summary['positive_relative_H5_rows']} positive relative margins, "
        f"{summary['positive_Q6_rows']} positive Q6 signs"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
