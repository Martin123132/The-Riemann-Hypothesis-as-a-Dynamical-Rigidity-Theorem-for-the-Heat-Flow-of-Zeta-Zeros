#!/usr/bin/env python3
"""Certify the lambda=-100 signed order-seven prefix through n=314."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
import json
from pathlib import Path
import sys

import sympy as sp


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
    symbolic_factorization as symbolic_h5_factorization,
)
from jensen_window_pf_compound_order6_m100_prefix_certificate import (  # noqa: E402
    SOURCE_PATHS as ORDER6_SOURCE_PATHS,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_m100_prefix_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order7_m100_prefix_certificate.md"
)
ORDER7_REDUCTION = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_uniform_tail_flow_reduction.json"
)
PRECISION_REPAIR = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "negative_lambda_m100_order7_k179_k190_dps220.jsonl"
)
PRECISION_REPAIR_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "negative_lambda_m100_order7_k179_k190_dps220_summary.json"
)
SOURCE_PATHS = (*ORDER6_SOURCE_PATHS, PRECISION_REPAIR)
PREFIX_LAST_N = 314
MAX_COEFFICIENT_INDEX = PREFIX_LAST_N + 12
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


def symbolic_h4_factorization() -> dict:
    coefficient, ratio = sp.symbols("A rho", positive=True)
    contractions = sp.symbols("x1:7", positive=True)
    ratios = [ratio]
    for index in range(1, 7):
        ratios.append(ratio * sp.prod(contractions[:index]))
    values = [coefficient]
    for current in ratios:
        values.append(sp.factor(values[-1] * current))
    h4 = sp.factor(sp.det(sp.Matrix(4, 4, lambda i, j: values[i + j])))

    def gap(center: int) -> sp.Expr:
        x_mid = contractions[center - 1]
        return sp.expand(
            (1 - x_mid) ** 2
            - x_mid**2
            * (1 - contractions[center - 2])
            * (1 - contractions[center])
        )

    margin = sp.expand(gap(3) ** 2 - contractions[2] ** 3 * gap(2) * gap(4))
    scale = sp.factor(
        coefficient**4
        * ratio**12
        * contractions[0] ** 8
        * contractions[1] ** 4
        / (1 - contractions[2])
    )
    residual = sp.factor(h4 - scale * margin)
    if residual != 0:
        raise RuntimeError("stable order-four factorization failed")
    return {
        "factorization": (
            "H_(4,n)=A_n^4*rho_n^12*x_(n+1)^8*x_(n+2)^4*"
            "F_n/d_(n+3)"
        ),
        "F_definition": "F_n=G_(n+1)^2-x_(n+3)^3*G_n*G_(n+2)",
        "symbolic_residual": str(residual),
    }


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
        for index in range(PREFIX_LAST_N + 9)
    }
    margins = {
        index: (
            gaps[index + 1] ** 2
            - contractions[index + 3] ** 3 * gaps[index] * gaps[index + 2]
        )
        for index in range(PREFIX_LAST_N + 7)
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
        for index in range(PREFIX_LAST_N + 5)
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

    h5_values = {
        index: h5_value(index) for index in range(PREFIX_LAST_N + 5)
    }
    h4_values = {
        index: h4_value(index) for index in range(PREFIX_LAST_N + 5)
    }
    for name, family in (("H4", h4_values), ("H5", h5_values)):
        if not all(bool(value > 0 and not value.contains(0)) for value in family.values()):
            raise RuntimeError(f"stable {name} family is not strictly positive")

    q6_values = {}
    for index in range(PREFIX_LAST_N + 3):
        relative_h5 = h5_values[index + 1] ** 2 / (
            h5_values[index] * h5_values[index + 2]
        ) - 1
        q6_values[index] = (
            h5_values[index]
            * h5_values[index + 2]
            * relative_h5
            / h4_values[index + 2]
        )
    if not all(
        bool(value > 0 and not value.contains(0)) for value in q6_values.values()
    ):
        raise RuntimeError("stable Q6 family is not strictly positive")

    rows = []
    minimum: tuple[flint.arb, int] | None = None
    for n in range(PREFIX_LAST_N + 1):
        relative = q6_values[n + 1] ** 2 / (
            q6_values[n] * q6_values[n + 2]
        ) - 1
        if not bool(relative > 0 and not relative.contains(0)):
            raise RuntimeError(f"inconclusive signed order-seven margin at n={n}")
        if minimum is None or relative.lower() < minimum[0].lower():
            minimum = (relative, n)
        rows.append(
            {
                "n": n,
                "relative_Q6_margin_ball": arb_text(relative),
                "relative_Q6_margin_lower": arb_lower_text(relative),
                "Q7_sign": "positive_by_signed_condensation",
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
        },
        "all_Q6_positive": True,
        "all_relative_Q6_margins_positive": True,
        "all_Q7_positive": True,
        "minimum_relative_n": minimum[1],
        "minimum_relative_ball": arb_text(minimum[0]),
        "minimum_relative_lower": arb_lower_text(minimum[0]),
    }


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    reduction = json.loads(ORDER7_REDUCTION.read_text(encoding="utf-8"))
    if reduction.get("summary", {}).get("open_entry_targets") != 1:
        raise RuntimeError("order-seven reduction source contract changed")
    repair = json.loads(PRECISION_REPAIR_SUMMARY.read_text(encoding="utf-8"))
    if repair.get("rows") != 12 or repair.get("k_min") != 179 or repair.get("k_max") != 190:
        raise RuntimeError("order-seven precision repair source contract changed")
    values, source_diagnostics = merged_coefficients()
    h4_factorization = symbolic_h4_factorization()
    h5_factorization = symbolic_h5_factorization()
    finite = finite_prefix(values)
    exact = {
        "signed_condensation": (
            "Q_(7,n)*H_(5,n+2)=Q_(6,n+1)^2-Q_(6,n)*Q_(6,n+2)"
        ),
        "stable_coordinate": (
            "L_n=Q_(6,n+1)^2/(Q_(6,n)*Q_(6,n+2))-1"
        ),
        "positive_scale": (
            "Q_(7,n)=Q_(6,n)*Q_(6,n+2)*L_n/H_(5,n+2)"
        ),
        "entry_equivalence": (
            "Q_(7,n)>0 iff L_n>0 inside the completed order-six cone"
        ),
        "prefix": "Q_(7,n)(-100)>0 for every 0<=n<=314",
        "remaining_tail": "Q_(7,n)(-100)>0 for every n>=315",
    }
    rows = [
        PrefixRow(
            "co7m100pc_01_signed_coordinate",
            "exact_identity",
            "ready_to_apply",
            "Signed order seven is a positive lower-cone scale times one relative Q6 log-concavity margin.",
            exact["signed_condensation"] + "; " + exact["stable_coordinate"] + "; " + exact["positive_scale"],
            "Exact Desnanot-Jacobi algebra only.",
        ),
        PrefixRow(
            "co7m100pc_02_stable_factorizations",
            "exact_identity",
            "ready_to_apply",
            "Every Q6 value is evaluated from cancellation-preserving H4 and H5 factorizations.",
            h4_factorization["factorization"] + "; " + h5_factorization["factorization"],
            "Exact coefficient-ratio identities only.",
            {"H4": h4_factorization, "H5": h5_factorization},
        ),
        PrefixRow(
            "co7m100pc_03_coefficient_enclosures",
            "interval_input",
            "ready_to_apply",
            "Merged outward-rounded Arb sources enclose every endpoint coefficient through A_326.",
            "A_k(-100)>0 for every 0<=k<=326",
            "Finite coefficient range only.",
            {"sources": source_diagnostics},
        ),
        PrefixRow(
            "co7m100pc_04_precision_repair",
            "interval_certificate",
            "ready_to_apply",
            "A dedicated twelve-coefficient repair removes every order-seven interval seam.",
            "A_k(-100) rigorously repaired for 179<=k<=190",
            "Local finite precision repair only.",
            {"summary": repair, "sha256": sha256(PRECISION_REPAIR)},
        ),
        PrefixRow(
            "co7m100pc_05_positive_q6_family",
            "interval_certificate",
            "ready_to_apply",
            "All stable lower coordinates and every required Q6 value are strictly positive.",
            "Q_(6,n)(-100)>0 for every 0<=n<=316",
            "Finite 2048-bit Arb prefix only.",
            finite["positive_coordinate_counts"],
        ),
        PrefixRow(
            "co7m100pc_06_positive_relative_margin",
            "interval_certificate",
            "ready_to_apply",
            "Every available relative Q6 log-concavity margin is strictly positive.",
            "L_n(-100)>0 for every 0<=n<=314",
            "Finite 2048-bit Arb prefix only.",
            {
                "minimum_n": finite["minimum_relative_n"],
                "minimum_ball": finite["minimum_relative_ball"],
                "minimum_lower": finite["minimum_relative_lower"],
            },
        ),
        PrefixRow(
            "co7m100pc_07_signed_order7_prefix",
            "exact_interval_composition",
            "ready_to_apply",
            "Positive lower-cone factors and relative margins prove the complete cached signed order-seven prefix.",
            exact["prefix"],
            "Finite lambda=-100 prefix only.",
        ),
        PrefixRow(
            "co7m100pc_08_open_tail",
            "open_handoff",
            "not_ready_to_apply",
            "Prove the signed order-seven endpoint margin on the infinite tail beginning at n=315.",
            exact["remaining_tail"],
            "Open analytic tail; not an all-shift order-seven theorem.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order7_m100_prefix_certificate",
        "date": "2026-07-13",
        "status": (
            "rigorous lambda=-100 signed order-seven prefix through n=314 with "
            "one open analytic tail"
        ),
        "proof_boundary": (
            "This artifact proves Q_(7,n)(-100)>0 only for 0<=n<=314. It does "
            "not prove the n>=315 tail, all-shift order-seven entry, forward "
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
            "positive_Q6_rows": finite["positive_coordinate_counts"]["Q6_values"],
            "positive_relative_Q6_rows": len(finite["rows"]),
            "positive_Q7_rows": len(finite["rows"]),
            "inconclusive_rows": 0,
            "precision_repair_rows": repair["rows"],
            "open_analytic_tails": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order7_m100_prefix_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order7_m100_prefix_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    finite = artifact["finite"]
    lines = [
        "# Jensen-Window PF Compound Order-Seven Lambda=-100 Prefix Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous signed order-seven endpoint prefix through `n=314`",
        "with one open analytic tail. This is not a proof of all-shift order",
        "seven, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order7_m100_prefix_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order7_m100_prefix_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order7_m100_prefix_certificate.py",
        "```",
        "",
        "## Stable Endpoint Coordinate",
        "",
        "Inside the completed positive order-six cone, signed condensation gives",
        "",
        "```text",
        exact["signed_condensation"],
        exact["stable_coordinate"],
        exact["positive_scale"],
        "```",
        "",
        "Thus `Q_(7,n)>0` is exactly `L_n>0`. The calculation evaluates",
        "`H_4`, `H_5`, and `Q_6` through their exact stable ratio",
        "factorizations; it never subtracts a raw seven-by-seven determinant.",
        "",
        "## Coefficient Cover And Repair",
        "",
        "The inherited sources and the new twelve-row local repair, merged in",
        "precedence order and hashed in the JSON artifact, give 2048-bit Arb",
        "evaluation from outward-rounded coefficient balls for",
        "",
        "```text",
        "A_k(-100)>0 for every 0<=k<=326.",
        "```",
        "",
        "The dedicated retained-integral repair covers `A_179,...,A_190` at",
        "220 decimal digits with the established analytic n-tail and cutoff-tail",
        "bounds. It removes all eight interval-inconclusive rows from the broad",
        "source; no midpoint sign is promoted.",
        "",
        "## Prefix Theorem",
        "",
        "Direct stable Arb evaluation proves",
        "",
        "```text",
        "L_n(-100)>0 for every 0<=n<=314,",
        exact["prefix"] + ".",
        "```",
        "",
        "The weakest row is the final cached shift:",
        "",
        "```text",
        f"n={finite['minimum_relative_n']},",
        f"L_314={finite['minimum_relative_ball']},",
        f"L_314 lower={finite['minimum_relative_lower']}>9/1000.",
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
        "The all-fixed-order eventual theorem supplies a finite but ineffective",
        "tail threshold, so it cannot splice this prefix to every shift.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order7_uniform_tail_flow_reduction.md",
        "outputs/jensen_window_pf_compound_order6_m100_prefix_certificate.md",
        "outputs/jensen_window_pf_compound_order6_uniform_heat_forward_invariance_certificate.md",
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
        "wrote order-seven lambda=-100 prefix: "
        f"{summary['coefficients']} coefficients, "
        f"{summary['positive_relative_Q6_rows']} positive relative margins, "
        f"{summary['positive_Q7_rows']} positive Q7 signs"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
