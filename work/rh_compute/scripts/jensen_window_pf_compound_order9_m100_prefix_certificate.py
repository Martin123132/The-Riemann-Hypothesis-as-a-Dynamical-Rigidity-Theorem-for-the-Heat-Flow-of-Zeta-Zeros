#!/usr/bin/env python3
"""Certify the lambda=-100 signed order-nine prefix through n=1240."""

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
    symbolic_h4_factorization,
    symbolic_h5_factorization,
)
from jensen_window_pf_compound_order8_m100_prefix_certificate import (  # noqa: E402
    SOURCE_PATHS as ORDER8_SOURCE_PATHS,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_m100_prefix_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order9_m100_prefix_certificate.md"
)
ORDER9_REDUCTION = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_uniform_tail_flow_reduction.json"
)
ORDER9_PRECISION_REPAIR = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "negative_lambda_m100_order9_k153_k178_dps220.jsonl"
)
ORDER9_PRECISION_REPAIR_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "negative_lambda_m100_order9_k153_k178_dps220_summary.json"
)
SOURCE_PATHS = (*ORDER8_SOURCE_PATHS, ORDER9_PRECISION_REPAIR)
PREFIX_LAST_N = 1240
MAX_COEFFICIENT_INDEX = PREFIX_LAST_N + 16
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


def merged_coefficients(
    source_paths: tuple[Path, ...] = SOURCE_PATHS,
    max_coefficient_index: int = MAX_COEFFICIENT_INDEX,
) -> tuple[dict[int, flint.arb], list[dict]]:
    values: dict[int, flint.arb] = {}
    diagnostics = []
    for precedence, source in enumerate(source_paths):
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
    expected = set(range(max_coefficient_index + 1))
    if set(values) != expected:
        missing = sorted(expected - set(values))
        extra = sorted(set(values) - expected)
        raise RuntimeError(f"coefficient coverage mismatch: missing={missing}, extra={extra}")
    if not all(bool(value > 0 and not value.contains(0)) for value in values.values()):
        raise RuntimeError("coefficient source contains a nonpositive ball")
    return values, diagnostics


def finite_prefix(
    values: dict[int, flint.arb],
    prefix_last_n: int = PREFIX_LAST_N,
    max_coefficient_index: int = MAX_COEFFICIENT_INDEX,
) -> dict:
    flint.ctx.prec = PRECISION_BITS
    contractions = {
        index: values[index - 1] * values[index + 1] / values[index] ** 2
        for index in range(1, max_coefficient_index)
    }
    defects = {index: 1 - value for index, value in contractions.items()}
    gaps = {
        index: (
            defects[index + 2] ** 2
            - contractions[index + 2] ** 2
            * defects[index + 1]
            * defects[index + 3]
        )
        for index in range(prefix_last_n + 13)
    }
    margins = {
        index: (
            gaps[index + 1] ** 2
            - contractions[index + 3] ** 3 * gaps[index] * gaps[index + 2]
        )
        for index in range(prefix_last_n + 11)
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
        for index in range(prefix_last_n + 9)
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

    h5_values = {index: h5_value(index) for index in range(prefix_last_n + 9)}
    h4_values = {index: h4_value(index) for index in range(prefix_last_n + 9)}
    for name, family in (("H4", h4_values), ("H5", h5_values)):
        if not all(bool(value > 0 and not value.contains(0)) for value in family.values()):
            raise RuntimeError(f"stable {name} family is not strictly positive")

    q6_values = {}
    for index in range(prefix_last_n + 7):
        relative_h5 = h5_values[index + 1] ** 2 / (
            h5_values[index] * h5_values[index + 2]
        ) - 1
        q6_values[index] = (
            h5_values[index]
            * h5_values[index + 2]
            * relative_h5
            / h4_values[index + 2]
        )
    q7_values = {}
    for index in range(prefix_last_n + 5):
        relative_q6 = q6_values[index + 1] ** 2 / (
            q6_values[index] * q6_values[index + 2]
        ) - 1
        q7_values[index] = (
            q6_values[index]
            * q6_values[index + 2]
            * relative_q6
            / h5_values[index + 2]
        )
    q8_values = {}
    for index in range(prefix_last_n + 3):
        relative_q7 = q7_values[index + 1] ** 2 / (
            q7_values[index] * q7_values[index + 2]
        ) - 1
        q8_values[index] = (
            q7_values[index]
            * q7_values[index + 2]
            * relative_q7
            / q6_values[index + 2]
        )
    for name, family in (
        ("Q6", q6_values),
        ("Q7", q7_values),
        ("Q8", q8_values),
    ):
        if not all(bool(value > 0 and not value.contains(0)) for value in family.values()):
            raise RuntimeError(f"stable {name} family is not strictly positive")

    rows = []
    minimum: tuple[flint.arb, int] | None = None
    for n in range(prefix_last_n + 1):
        relative = q8_values[n + 1] ** 2 / (
            q8_values[n] * q8_values[n + 2]
        ) - 1
        if not bool(relative > 0 and not relative.contains(0)):
            raise RuntimeError(f"inconclusive signed order-nine margin at n={n}")
        if minimum is None or relative.lower() < minimum[0].lower():
            minimum = (relative, n)
        rows.append(
            {
                "n": n,
                "relative_Q8_margin_ball": arb_text(relative),
                "relative_Q8_margin_lower": arb_lower_text(relative),
                "Q9_sign": "positive_by_signed_condensation",
            }
        )
    assert minimum is not None
    return {
        "lambda": "-100",
        "n_range": [0, prefix_last_n],
        "coefficient_range": [0, max_coefficient_index],
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
            "Q8_values": len(q8_values),
        },
        "all_Q6_positive": True,
        "all_Q7_positive": True,
        "all_Q8_positive": True,
        "all_relative_Q8_margins_positive": True,
        "all_Q9_positive": True,
        "minimum_relative_n": minimum[1],
        "minimum_relative_ball": arb_text(minimum[0]),
        "minimum_relative_lower": arb_lower_text(minimum[0]),
    }


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    reduction = json.loads(ORDER9_REDUCTION.read_text(encoding="utf-8"))
    if reduction.get("summary", {}).get("open_entry_targets") != 1:
        raise RuntimeError("order-nine reduction source contract changed")
    inherited_repair = json.loads(PRECISION_REPAIR_SUMMARY.read_text(encoding="utf-8"))
    if (
        inherited_repair.get("rows") != 12
        or inherited_repair.get("k_min") != 179
        or inherited_repair.get("k_max") != 190
    ):
        raise RuntimeError("inherited precision repair source contract changed")
    order9_repair = json.loads(
        ORDER9_PRECISION_REPAIR_SUMMARY.read_text(encoding="utf-8")
    )
    if (
        order9_repair.get("rows") != 26
        or order9_repair.get("k_min") != 153
        or order9_repair.get("k_max") != 178
        or order9_repair.get("lambdas") != ["-100.0"]
    ):
        raise RuntimeError("order-nine precision repair source contract changed")
    values, source_diagnostics = merged_coefficients()
    h4_factorization = symbolic_h4_factorization()
    h5_factorization = symbolic_h5_factorization()
    finite = finite_prefix(values)
    exact = {
        "signed_condensation": (
            "Q_(9,n)*Q_(7,n+2)=Q_(8,n+1)^2-Q_(8,n)*Q_(8,n+2)"
        ),
        "stable_coordinate": (
            "M_n=Q_(8,n+1)^2/(Q_(8,n)*Q_(8,n+2))-1"
        ),
        "positive_scale": (
            "Q_(9,n)=Q_(8,n)*Q_(8,n+2)*M_n/Q_(7,n+2)"
        ),
        "entry_equivalence": (
            "Q_(9,n)>0 iff M_n>0 inside the completed order-eight cone"
        ),
        "prefix": "Q_(9,n)(-100)>0 for every 0<=n<=1240",
        "remaining_tail": "Q_(9,n)(-100)>0 for every n>=1241",
    }
    rows = [
        PrefixRow(
            "co9m100pc_01_signed_coordinate",
            "exact_identity",
            "ready_to_apply",
            "Signed order nine is a positive lower-cone scale times one relative Q8 log-concavity margin.",
            exact["signed_condensation"] + "; " + exact["stable_coordinate"] + "; " + exact["positive_scale"],
            "Exact signed Desnanot-Jacobi algebra only.",
        ),
        PrefixRow(
            "co9m100pc_02_stable_factorizations",
            "exact_identity",
            "ready_to_apply",
            "Every Q8 value is evaluated through cancellation-preserving H4, H5, Q6, Q7, and Q8 coordinates.",
            h4_factorization["factorization"] + "; " + h5_factorization["factorization"],
            "Exact coefficient-ratio identities only.",
            {"H4": h4_factorization, "H5": h5_factorization},
        ),
        PrefixRow(
            "co9m100pc_03_coefficient_enclosures",
            "interval_input",
            "ready_to_apply",
            "Merged outward-rounded Arb sources enclose every endpoint coefficient through A_1256.",
            "A_k(-100)>0 for every 0<=k<=1256",
            "Finite coefficient range only.",
            {"sources": source_diagnostics},
        ),
        PrefixRow(
            "co9m100pc_04_precision_repairs",
            "interval_certificate",
            "ready_to_apply",
            "Two retained-integral repairs preserve all sixth-nested interval seams.",
            "A_k(-100) rigorously repaired for 153<=k<=190",
            "Two contiguous finite high-precision repairs only.",
            {
                "order9": order9_repair,
                "order9_sha256": sha256(ORDER9_PRECISION_REPAIR),
                "inherited": inherited_repair,
                "inherited_sha256": sha256(PRECISION_REPAIR),
            },
        ),
        PrefixRow(
            "co9m100pc_05_positive_q8_family",
            "interval_certificate",
            "ready_to_apply",
            "All stable lower coordinates and every required Q8 value are strictly positive.",
            "Q_(8,n)(-100)>0 for every 0<=n<=1242",
            "Finite 2048-bit Arb prefix only.",
            finite["positive_coordinate_counts"],
        ),
        PrefixRow(
            "co9m100pc_06_positive_relative_margin",
            "interval_certificate",
            "ready_to_apply",
            "Every available relative Q8 log-concavity margin is strictly positive.",
            "M_n(-100)>0 for every 0<=n<=1240",
            "Finite 2048-bit Arb prefix only.",
            {
                "minimum_n": finite["minimum_relative_n"],
                "minimum_ball": finite["minimum_relative_ball"],
                "minimum_lower": finite["minimum_relative_lower"],
            },
        ),
        PrefixRow(
            "co9m100pc_07_signed_order9_prefix",
            "exact_interval_composition",
            "ready_to_apply",
            "Positive lower-cone factors and relative margins prove the complete cached signed order-nine prefix.",
            exact["prefix"],
            "Finite lambda=-100 prefix only.",
        ),
        PrefixRow(
            "co9m100pc_08_open_tail",
            "open_handoff",
            "not_ready_to_apply",
            "Prove the signed order-nine endpoint margin on the infinite tail beginning at n=1241.",
            exact["remaining_tail"],
            "Open analytic tail; not an all-shift order-nine theorem.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order9_m100_prefix_certificate",
        "date": "2026-07-13",
        "status": (
            "rigorous lambda=-100 signed order-nine prefix through n=1240 with "
            "one open analytic tail"
        ),
        "proof_boundary": (
            "This artifact proves Q_(9,n)(-100)>0 only for 0<=n<=1240. It does "
            "not prove the n>=1241 tail, all-shift order-nine entry, forward "
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
            "positive_Q8_rows": finite["positive_coordinate_counts"]["Q8_values"],
            "positive_relative_Q8_rows": len(finite["rows"]),
            "positive_Q9_rows": len(finite["rows"]),
            "inconclusive_rows": 0,
            "precision_repair_rows": (
                inherited_repair["rows"] + order9_repair["rows"]
            ),
            "open_analytic_tails": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order9_m100_prefix_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order9_m100_prefix_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    finite = artifact["finite"]
    lines = [
        "# Jensen-Window PF Compound Order-Nine Lambda=-100 Prefix Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous signed order-nine endpoint prefix through `n=1240`",
        "with one open analytic tail. This is not a proof of all-shift order",
        "nine, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order9_m100_prefix_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order9_m100_prefix_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order9_m100_prefix_certificate.py",
        "```",
        "",
        "## Stable Endpoint Coordinate",
        "",
        "```text",
        exact["signed_condensation"],
        exact["stable_coordinate"],
        exact["positive_scale"],
        "```",
        "",
        "Thus `Q_(9,n)>0` is exactly `M_n>0`. No raw nine-by-nine",
        "near-cancelling determinant is evaluated.",
        "",
        "## Coefficient Cover And Repair",
        "",
        "The inherited sources and two rigorous retained-integral repairs give",
        "2048-bit Arb coefficient balls through `A_1256`. The new repair covers",
        "`A_153,...,A_178`; the inherited repair covers `A_179,...,A_190`.",
        "All 38 repaired rows use `n_sum=70`, cutoff 7, and 220 decimal digits.",
        "",
        "## Prefix Theorem",
        "",
        "```text",
        "Q_(8,n)(-100)>0 for every 0<=n<=1242,",
        "M_n(-100)>0 for every 0<=n<=1240,",
        exact["prefix"] + ".",
        "```",
        "",
        "The weakest row is the final cached shift:",
        "",
        "```text",
        f"n={finite['minimum_relative_n']},",
        f"M_{finite['minimum_relative_n']}={finite['minimum_relative_ball']},",
        f"M_{finite['minimum_relative_n']} lower={finite['minimum_relative_lower']}>1/250.",
        "```",
        "",
        "## Remaining Tail",
        "",
        "```text",
        exact["remaining_tail"] + ".",
        "```",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order9_uniform_tail_flow_reduction.md",
        "outputs/jensen_window_pf_compound_order8_m100_prefix_certificate.md",
        "outputs/jensen_window_pf_compound_order8_uniform_heat_forward_invariance_certificate.md",
        "outputs/signed_hankel_jensen_bridge_target.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-nine lambda=-100 prefix: "
        f"{summary['coefficients']} coefficients, "
        f"{summary['positive_relative_Q8_rows']} positive relative margins, "
        f"{summary['positive_Q9_rows']} positive Q9 signs"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
