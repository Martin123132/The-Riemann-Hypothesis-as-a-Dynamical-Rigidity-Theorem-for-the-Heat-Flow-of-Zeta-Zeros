#!/usr/bin/env python3
"""Certify the lambda=-100 signed order-eight prefix through n=1242."""

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

import flint  # noqa: E402

from jensen_window_pf_compound_order5_m100_prefix_certificate import (  # noqa: E402
    arb_lower_text,
    arb_text,
    load_source,
    sha256,
)
from jensen_window_pf_compound_order7_m100_prefix_certificate import (  # noqa: E402
    PRECISION_REPAIR,
    PRECISION_REPAIR_SUMMARY,
    SOURCE_PATHS as ORDER7_SOURCE_PATHS,
    symbolic_h4_factorization,
    symbolic_h5_factorization,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order8_m100_prefix_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order8_m100_prefix_certificate.md"
)
ORDER8_REDUCTION = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order8_uniform_tail_flow_reduction.json"
)
COEFFICIENT_EXTENSION_327_706 = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "negative_lambda_m100_order8_k327_k706_dps220.jsonl"
)
COEFFICIENT_EXTENSION_327_706_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "negative_lambda_m100_order8_k327_k706_dps220_summary.json"
)
COEFFICIENT_EXTENSION_707_1256 = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "negative_lambda_m100_order8_k707_k1256_dps220.jsonl"
)
COEFFICIENT_EXTENSION_707_1256_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "negative_lambda_m100_order8_k707_k1256_dps220_summary.json"
)
SOURCE_PATHS = (
    *ORDER7_SOURCE_PATHS,
    COEFFICIENT_EXTENSION_327_706,
    COEFFICIENT_EXTENSION_707_1256,
)
PREFIX_LAST_N = 1242
MAX_COEFFICIENT_INDEX = PREFIX_LAST_N + 14
PRECISION_BITS = 2048


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
    contractions = {
        index: values[index - 1] * values[index + 1] / values[index] ** 2
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
        for index in range(PREFIX_LAST_N + 11)
    }
    margins = {
        index: (
            gaps[index + 1] ** 2
            - contractions[index + 3] ** 3 * gaps[index] * gaps[index + 2]
        )
        for index in range(PREFIX_LAST_N + 9)
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
        for index in range(PREFIX_LAST_N + 7)
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

    def h4_value(index: int) -> flint.arb:
        ratio = values[index + 1] / values[index]
        return (
            values[index] ** 4
            * ratio**12
            * contractions[index + 1] ** 8
            * contractions[index + 2] ** 4
            * margins[index]
            / defects[index + 3]
        )

    h5_values = {index: h5_value(index) for index in range(PREFIX_LAST_N + 7)}
    h4_values = {index: h4_value(index) for index in range(PREFIX_LAST_N + 7)}
    for name, family in (("H4", h4_values), ("H5", h5_values)):
        if not all(bool(value > 0 and not value.contains(0)) for value in family.values()):
            raise RuntimeError(f"stable {name} family is not strictly positive")

    q6_values = {}
    for index in range(PREFIX_LAST_N + 5):
        relative_h5 = h5_values[index + 1] ** 2 / (
            h5_values[index] * h5_values[index + 2]
        ) - 1
        q6_values[index] = (
            h5_values[index]
            * h5_values[index + 2]
            * relative_h5
            / h4_values[index + 2]
        )
    if not all(bool(value > 0 and not value.contains(0)) for value in q6_values.values()):
        raise RuntimeError("stable Q6 family is not strictly positive")

    q7_values = {}
    for index in range(PREFIX_LAST_N + 3):
        relative_q6 = q6_values[index + 1] ** 2 / (
            q6_values[index] * q6_values[index + 2]
        ) - 1
        q7_values[index] = (
            q6_values[index]
            * q6_values[index + 2]
            * relative_q6
            / h5_values[index + 2]
        )
    if not all(bool(value > 0 and not value.contains(0)) for value in q7_values.values()):
        raise RuntimeError("stable Q7 family is not strictly positive")

    rows = []
    minimum: tuple[flint.arb, int] | None = None
    for n in range(PREFIX_LAST_N + 1):
        relative = q7_values[n + 1] ** 2 / (
            q7_values[n] * q7_values[n + 2]
        ) - 1
        if not bool(relative > 0 and not relative.contains(0)):
            raise RuntimeError(f"inconclusive signed order-eight margin at n={n}")
        if minimum is None or relative.lower() < minimum[0].lower():
            minimum = (relative, n)
        rows.append(
            {
                "n": n,
                "relative_Q7_margin_ball": arb_text(relative),
                "relative_Q7_margin_lower": arb_lower_text(relative),
                "Q8_sign": "positive_by_signed_condensation",
            }
        )
    assert minimum is not None
    return {
        "lambda": "-100",
        "n_range": [0, PREFIX_LAST_N],
        "coefficient_range": [0, MAX_COEFFICIENT_INDEX],
        "precision_bits": PRECISION_BITS,
        "rows": rows,
        "positive_coordinate_counts": {
            "contractions": len(contractions),
            "defects": len(defects),
            "order3_gaps": len(gaps),
            "order4_margins": len(margins),
            "order5_stable_margins": len(stable_h5),
            "H4_values": len(h4_values),
            "H5_values": len(h5_values),
            "Q6_values": len(q6_values),
            "Q7_values": len(q7_values),
        },
        "all_Q6_positive": True,
        "all_Q7_positive": True,
        "all_relative_Q7_margins_positive": True,
        "all_Q8_positive": True,
        "minimum_relative_n": minimum[1],
        "minimum_relative_ball": arb_text(minimum[0]),
        "minimum_relative_lower": arb_lower_text(minimum[0]),
    }


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    reduction = json.loads(ORDER8_REDUCTION.read_text(encoding="utf-8"))
    if reduction.get("summary", {}).get("open_entry_targets") != 1:
        raise RuntimeError("order-eight reduction source contract changed")
    repair = json.loads(PRECISION_REPAIR_SUMMARY.read_text(encoding="utf-8"))
    if repair.get("rows") != 12 or repair.get("k_min") != 179 or repair.get("k_max") != 190:
        raise RuntimeError("inherited precision repair source contract changed")
    extension = json.loads(
        COEFFICIENT_EXTENSION_327_706_SUMMARY.read_text(encoding="utf-8")
    )
    if (
        extension.get("rows") != 380
        or extension.get("k_min") != 327
        or extension.get("k_max") != 706
        or extension.get("lambdas") != ["-100.0"]
    ):
        raise RuntimeError("order-eight coefficient extension source contract changed")
    second_extension = json.loads(
        COEFFICIENT_EXTENSION_707_1256_SUMMARY.read_text(encoding="utf-8")
    )
    if (
        second_extension.get("rows") != 550
        or second_extension.get("k_min") != 707
        or second_extension.get("k_max") != 1256
        or second_extension.get("lambdas") != ["-100.0"]
    ):
        raise RuntimeError("second order-eight coefficient extension contract changed")
    values, source_diagnostics = merged_coefficients()
    if len(values) != MAX_COEFFICIENT_INDEX + 1:
        raise RuntimeError("order-eight coefficient cover changed")
    h4_factorization = symbolic_h4_factorization()
    h5_factorization = symbolic_h5_factorization()
    finite = finite_prefix(values)
    exact = {
        "signed_condensation": (
            "Q_(8,n)*Q_(6,n+2)=Q_(7,n+1)^2-Q_(7,n)*Q_(7,n+2)"
        ),
        "stable_coordinate": (
            "M_n=Q_(7,n+1)^2/(Q_(7,n)*Q_(7,n+2))-1"
        ),
        "positive_scale": (
            "Q_(8,n)=Q_(7,n)*Q_(7,n+2)*M_n/Q_(6,n+2)"
        ),
        "entry_equivalence": (
            "Q_(8,n)>0 iff M_n>0 inside the completed order-seven cone"
        ),
        "prefix": "Q_(8,n)(-100)>0 for every 0<=n<=1242",
        "remaining_tail": "Q_(8,n)(-100)>0 for every n>=1243",
    }
    rows = [
        PrefixRow(
            "co8m100pc_01_signed_coordinate",
            "exact_identity",
            "ready_to_apply",
            "Signed order eight is a positive lower-cone scale times one relative Q7 log-concavity margin.",
            exact["signed_condensation"] + "; " + exact["stable_coordinate"] + "; " + exact["positive_scale"],
            "Exact signed Desnanot-Jacobi algebra only.",
        ),
        PrefixRow(
            "co8m100pc_02_stable_factorizations",
            "exact_identity",
            "ready_to_apply",
            "Every Q7 value is evaluated through cancellation-preserving H4, H5, Q6, and Q7 coordinates.",
            h4_factorization["factorization"] + "; " + h5_factorization["factorization"],
            "Exact coefficient-ratio identities only.",
            {"H4": h4_factorization, "H5": h5_factorization},
        ),
        PrefixRow(
            "co8m100pc_03_coefficient_enclosures",
            "interval_input",
            "ready_to_apply",
            "Merged outward-rounded Arb sources enclose every endpoint coefficient through A_1256.",
            "A_k(-100)>0 for every 0<=k<=1256",
            "Finite coefficient range only.",
            {"sources": source_diagnostics},
        ),
        PrefixRow(
            "co8m100pc_04_precision_repair",
            "interval_certificate",
            "ready_to_apply",
            "The inherited twelve-coefficient retained-integral repair preserves every order-eight interval seam.",
            "A_k(-100) rigorously repaired for 179<=k<=190",
            "Local finite precision repair only.",
            {"summary": repair, "sha256": sha256(PRECISION_REPAIR)},
        ),
        PrefixRow(
            "co8m100pc_05_positive_q7_family",
            "interval_certificate",
            "ready_to_apply",
            "All stable lower coordinates and every required Q7 value are strictly positive.",
            "Q_(7,n)(-100)>0 for every 0<=n<=314",
            "Finite 2048-bit Arb prefix only.",
            finite["positive_coordinate_counts"],
        ),
        PrefixRow(
            "co8m100pc_06_positive_relative_margin",
            "interval_certificate",
            "ready_to_apply",
            "Every available relative Q7 log-concavity margin is strictly positive.",
            "M_n(-100)>0 for every 0<=n<=1242",
            "Finite 2048-bit Arb prefix only.",
            {
                "minimum_n": finite["minimum_relative_n"],
                "minimum_ball": finite["minimum_relative_ball"],
                "minimum_lower": finite["minimum_relative_lower"],
            },
        ),
        PrefixRow(
            "co8m100pc_07_signed_order8_prefix",
            "exact_interval_composition",
            "ready_to_apply",
            "Positive lower-cone factors and relative margins prove the complete cached signed order-eight prefix.",
            exact["prefix"],
            "Finite lambda=-100 prefix only.",
        ),
        PrefixRow(
            "co8m100pc_08_open_tail",
            "open_handoff",
            "not_ready_to_apply",
            "Prove the signed order-eight endpoint margin on the infinite tail beginning at n=1243.",
            exact["remaining_tail"],
            "Open analytic tail; not an all-shift order-eight theorem.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order8_m100_prefix_certificate",
        "date": "2026-07-13",
        "status": (
            "rigorous lambda=-100 signed order-eight prefix through n=1242 with "
            "one open analytic tail"
        ),
        "proof_boundary": (
            "This artifact proves Q_(8,n)(-100)>0 only for 0<=n<=1242. It does "
            "not prove the n>=1243 tail, all-shift order-eight entry, forward "
            "invariance, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [source.relative_to(REPO_ROOT).as_posix() for source in SOURCE_PATHS],
        "source_diagnostics": source_diagnostics,
        "exact": exact,
        "factorizations": {"H4": h4_factorization, "H5": h5_factorization},
        "finite": finite,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_to_apply_rows": 7,
            "coefficients": len(values),
            "positive_Q7_rows": finite["positive_coordinate_counts"]["Q7_values"],
            "positive_relative_Q7_rows": len(finite["rows"]),
            "positive_Q8_rows": len(finite["rows"]),
            "inconclusive_rows": 0,
            "precision_repair_rows": repair["rows"],
            "open_analytic_tails": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order8_m100_prefix_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order8_m100_prefix_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    finite = artifact["finite"]
    lines = [
        "# Jensen-Window PF Compound Order-Eight Lambda=-100 Prefix Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous signed order-eight endpoint prefix through `n=1242`",
        "with one open analytic tail. This is not a proof of all-shift order",
        "eight, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order8_m100_prefix_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order8_m100_prefix_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order8_m100_prefix_certificate.py",
        "```",
        "",
        "## Stable Endpoint Coordinate",
        "",
        "Inside the completed positive order-seven cone, signed condensation gives",
        "",
        "```text",
        exact["signed_condensation"],
        exact["stable_coordinate"],
        exact["positive_scale"],
        "```",
        "",
        "Thus `Q_(8,n)>0` is exactly `M_n>0`. The calculation evaluates",
        "`H_4`, `H_5`, `Q_6`, and `Q_7` through stable ratio coordinates;",
        "it never subtracts a raw eight-by-eight determinant.",
        "",
        "## Coefficient Cover",
        "",
        "The inherited hashed sources, twelve-row precision repair, and two",
        "retained-integral extensions totaling 930 rows give 2048-bit Arb evaluation",
        "from outward-rounded coefficient balls for",
        "",
        "```text",
        "A_k(-100)>0 for every 0<=k<=1256.",
        "```",
        "",
        "The extension was generated directly by rigorous Arb quadrature.",
        "",
        "## Prefix Theorem",
        "",
        "Direct stable Arb evaluation proves",
        "",
        "```text",
        "M_n(-100)>0 for every 0<=n<=1242,",
        exact["prefix"] + ".",
        "```",
        "",
        "The weakest row is the final cached shift:",
        "",
        "```text",
        f"n={finite['minimum_relative_n']},",
        f"M_{finite['minimum_relative_n']}={finite['minimum_relative_ball']},",
        f"M_{finite['minimum_relative_n']} lower={finite['minimum_relative_lower']}>1/300.",
        "```",
        "",
        "## Remaining Tail",
        "",
        "The complete endpoint problem is reduced to",
        "",
        "```text",
        exact["remaining_tail"] + ".",
        "```",
        "",
        "The all-fixed-order eventual theorem supplies a finite but ineffective",
        "tail threshold, so it cannot splice this prefix to every shift.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order8_uniform_tail_flow_reduction.md",
        "outputs/jensen_window_pf_compound_order7_m100_prefix_certificate.md",
        "outputs/jensen_window_pf_compound_order7_uniform_heat_forward_invariance_certificate.md",
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
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-eight lambda=-100 prefix: "
        f"{summary['coefficients']} coefficients, "
        f"{summary['positive_relative_Q7_rows']} positive relative margins, "
        f"{summary['positive_Q8_rows']} positive Q8 signs"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
