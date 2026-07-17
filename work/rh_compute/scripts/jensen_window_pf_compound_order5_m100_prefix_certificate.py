#!/usr/bin/env python3
"""Certify the lambda=-100 contiguous order-five prefix through n=316."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists() and str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order5_m100_prefix_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order5_m100_prefix_certificate.md"
)
ORDER5_REDUCTION = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order5_uniform_tail_flow_reduction.json"
)
SOURCE_PATHS = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k300.jsonl",
    REPO_ROOT
    / "work/rh_compute/results/"
    "acb_enclosures_negative_lambda_m100_order3_k208_k217_dps220.jsonl",
    REPO_ROOT
    / "work/rh_compute/results/"
    "acb_enclosures_negative_lambda_cone_entry_lam_m100_k220_k250_dps220.jsonl",
    REPO_ROOT
    / "work/rh_compute/results/"
    "acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220.jsonl",
    REPO_ROOT
    / "work/rh_compute/results/"
    "acb_enclosures_negative_lambda_m100_order3_k321_dps220.jsonl",
    REPO_ROOT
    / "work/rh_compute/results/negative_lambda_m100_order4_k322_dps220.jsonl",
    REPO_ROOT
    / "work/rh_compute/results/negative_lambda_m100_order5_k323_k324_dps220.jsonl",
)
PREFIX_LAST_N = 316
MAX_COEFFICIENT_INDEX = PREFIX_LAST_N + 8
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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def arb_text(value: flint.arb, digits: int = 60) -> str:
    return value.str(digits).replace("e", "E")


def arb_lower_text(value: flint.arb, digits: int = 60) -> str:
    rounded = value.lower().str(digits, radius=False)
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(rounded).next_minus(), "E")


def load_source(path: Path) -> dict[int, flint.arb]:
    values: dict[int, flint.arb] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip():
                continue
            row = json.loads(raw)
            if row.get("kind") != "acb_coefficient_enclosure":
                continue
            if row.get("lam") not in {"-100", "-100.0"}:
                continue
            values[int(row["k"])] = flint.arb(row["A_ball"])
    return values


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


def symbolic_factorization() -> dict:
    coefficient, ratio = sp.symbols("A rho", positive=True)
    contractions = sp.symbols("x1:8", positive=True)
    ratios = [ratio]
    for index in range(1, 8):
        ratios.append(ratio * sp.prod(contractions[:index]))
    values = [coefficient]
    for current in ratios:
        values.append(sp.factor(values[-1] * current))
    h5 = sp.factor(sp.det(sp.Matrix(5, 5, lambda i, j: values[i + j])))

    def gap(center: int) -> sp.Expr:
        x_mid = contractions[center - 1]
        return sp.expand(
            (1 - x_mid) ** 2
            - x_mid**2
            * (1 - contractions[center - 2])
            * (1 - contractions[center])
        )

    def margin(shift: int) -> sp.Expr:
        return sp.expand(
            gap(shift + 3) ** 2
            - contractions[shift + 2] ** 3
            * gap(shift + 2)
            * gap(shift + 4)
        )

    stable = sp.factor(
        (1 - contractions[2])
        * (1 - contractions[4])
        * margin(1) ** 2
        - contractions[3] ** 4
        * (1 - contractions[3]) ** 2
        * margin(0)
        * margin(2)
    )
    positive_scale = sp.factor(
        coefficient**5
        * ratio**20
        * contractions[0] ** 15
        * contractions[1] ** 10
        * contractions[2] ** 5
        / (
            (1 - contractions[2])
            * (1 - contractions[3]) ** 2
            * (1 - contractions[4])
            * gap(4)
        )
    )
    residual = sp.factor(h5 - positive_scale * stable)
    if residual != 0:
        raise RuntimeError("stable order-five factorization failed")
    return {
        "x_definition": "x_k=A_(k-1)*A_(k+1)/A_k^2, d_k=1-x_k",
        "G_definition": "G_n=d_(n+2)^2-x_(n+2)^2*d_(n+1)*d_(n+3)",
        "F_definition": "F_n=G_(n+1)^2-x_(n+3)^3*G_n*G_(n+2)",
        "J_definition": (
            "J_n=d_(n+3)*d_(n+5)*F_(n+1)^2-"
            "x_(n+4)^4*d_(n+4)^2*F_n*F_(n+2)"
        ),
        "positive_scale": (
            "W_n=A_n^5*rho_n^20*x_(n+1)^15*x_(n+2)^10*x_(n+3)^5/"
            "(d_(n+3)*d_(n+4)^2*d_(n+5)*G_(n+2))"
        ),
        "factorization": "H_(5,n)=W_n*J_n",
        "scale_positive_in_lower_cone": True,
        "symbolic_residual": str(residual),
    }


def finite_prefix(values: dict[int, flint.arb]) -> dict:
    flint.ctx.prec = PRECISION_BITS

    def contraction(index: int) -> flint.arb:
        return values[index - 1] * values[index + 1] / values[index] ** 2

    contractions = {
        index: contraction(index) for index in range(1, MAX_COEFFICIENT_INDEX)
    }
    defects = {index: 1 - value for index, value in contractions.items()}

    def gap(index: int) -> flint.arb:
        return (
            defects[index + 2] ** 2
            - contractions[index + 2] ** 2
            * defects[index + 1]
            * defects[index + 3]
        )

    gaps = {index: gap(index) for index in range(PREFIX_LAST_N + 5)}

    def margin(index: int) -> flint.arb:
        return (
            gaps[index + 1] ** 2
            - contractions[index + 3] ** 3
            * gaps[index]
            * gaps[index + 2]
        )

    margins = {index: margin(index) for index in range(PREFIX_LAST_N + 3)}

    def stable(index: int) -> flint.arb:
        return (
            defects[index + 3]
            * defects[index + 5]
            * margins[index + 1] ** 2
            - contractions[index + 4] ** 4
            * defects[index + 4] ** 2
            * margins[index]
            * margins[index + 2]
        )

    def relative(index: int) -> flint.arb:
        denominator = (
            contractions[index + 4] ** 4
            * defects[index + 4] ** 2
            * margins[index]
            * margins[index + 2]
        )
        return (
            defects[index + 3]
            * defects[index + 5]
            * margins[index + 1] ** 2
            / denominator
            - 1
        )

    for name, family in (
        ("contraction", contractions),
        ("defect", defects),
        ("order3_gap", gaps),
        ("order4_margin", margins),
    ):
        if not all(bool(value > 0 and not value.contains(0)) for value in family.values()):
            raise RuntimeError(f"{name} family is not strictly positive")

    rows = []
    minimum_stable: tuple[flint.arb, int] | None = None
    minimum_relative: tuple[flint.arb, int] | None = None
    for n in range(PREFIX_LAST_N + 1):
        stable_value = stable(n)
        relative_value = relative(n)
        if not bool(stable_value > 0 and not stable_value.contains(0)):
            raise RuntimeError(f"inconclusive stable order-five margin at n={n}")
        if not bool(relative_value > 0 and not relative_value.contains(0)):
            raise RuntimeError(f"inconclusive relative order-five margin at n={n}")
        if (
            minimum_stable is None
            or stable_value.lower() < minimum_stable[0].lower()
        ):
            minimum_stable = (stable_value, n)
        if (
            minimum_relative is None
            or relative_value.lower() < minimum_relative[0].lower()
        ):
            minimum_relative = (relative_value, n)
        rows.append(
            {
                "n": n,
                "J_ball": arb_text(stable_value),
                "J_lower": arb_lower_text(stable_value),
                "relative_ball": arb_text(relative_value),
                "relative_lower": arb_lower_text(relative_value),
                "H5_sign": "positive_by_exact_scale",
            }
        )

    assert minimum_stable is not None
    assert minimum_relative is not None
    return {
        "lambda": "-100",
        "n_range": [0, PREFIX_LAST_N],
        "coefficient_range": [0, MAX_COEFFICIENT_INDEX],
        "precision_bits": PRECISION_BITS,
        "rows": rows,
        "all_J_positive": True,
        "all_relative_positive": True,
        "all_H5_positive": True,
        "minimum_J_n": minimum_stable[1],
        "minimum_J_ball": arb_text(minimum_stable[0]),
        "minimum_J_lower": arb_lower_text(minimum_stable[0]),
        "minimum_relative_n": minimum_relative[1],
        "minimum_relative_ball": arb_text(minimum_relative[0]),
        "minimum_relative_lower": arb_lower_text(minimum_relative[0]),
    }


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    reduction = json.loads(ORDER5_REDUCTION.read_text(encoding="utf-8"))
    if reduction.get("summary", {}).get("open_entry_targets") != 1:
        raise RuntimeError("order-five reduction source contract changed")
    values, source_diagnostics = merged_coefficients()
    exact = symbolic_factorization()
    finite = finite_prefix(values)
    rows = [
        PrefixRow(
            id="co5m100pc_01_stable_factorization",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The raw order-five determinant is a positive lower-cone scale times one bounded stable margin.",
            formula=exact["factorization"] + "; " + exact["J_definition"],
            proof_boundary="Exact coefficient-ratio algebra only.",
            diagnostics=exact,
        ),
        PrefixRow(
            id="co5m100pc_02_coefficient_enclosures",
            role="interval_input",
            readiness="ready_to_apply",
            claim="Merged outward-rounded Arb sources enclose every A_k(-100) needed by the prefix.",
            formula="A_k(-100)>0 for every 0<=k<=324",
            proof_boundary="Finite coefficient range only.",
            diagnostics={"sources": source_diagnostics},
        ),
        PrefixRow(
            id="co5m100pc_03_positive_lower_cone",
            role="interval_certificate",
            readiness="ready_to_apply",
            claim="All denominator factors in the exact stable scale are strictly positive on the complete prefix.",
            formula="d_k>0, G_n>0, F_n>0 on every index used for 0<=n<=316",
            proof_boundary="Finite Arb prefix only.",
        ),
        PrefixRow(
            id="co5m100pc_04_positive_stable_margin",
            role="interval_certificate",
            readiness="ready_to_apply",
            claim="Every stable order-five margin J_n is strictly positive through n=316.",
            formula="J_n>0 for every 0<=n<=316",
            proof_boundary="Finite Arb prefix only.",
            diagnostics={
                "minimum_n": finite["minimum_J_n"],
                "minimum_ball": finite["minimum_J_ball"],
                "minimum_lower": finite["minimum_J_lower"],
            },
        ),
        PrefixRow(
            id="co5m100pc_05_positive_relative_margin",
            role="interval_certificate",
            readiness="ready_to_apply",
            claim="The scale-free relative order-five margin remains strictly positive throughout the prefix.",
            formula="relative_margin_n>0 for every 0<=n<=316",
            proof_boundary="Finite Arb prefix only.",
            diagnostics={
                "minimum_n": finite["minimum_relative_n"],
                "minimum_ball": finite["minimum_relative_ball"],
                "minimum_lower": finite["minimum_relative_lower"],
            },
        ),
        PrefixRow(
            id="co5m100pc_06_prefix_theorem",
            role="finite_theorem",
            readiness="ready_to_apply",
            claim="The actual contiguous order-five determinant is strictly positive through n=316 at lambda=-100.",
            formula="H_(5,n)(-100)>0 for every 0<=n<=316",
            proof_boundary="Rigorous finite prefix; no all-shift promotion.",
        ),
        PrefixRow(
            id="co5m100pc_07_open_tail",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="Prove the stable order-five margin for the infinite tail beginning at n=317.",
            formula="J_n(-100)>0 for every n>=317",
            proof_boundary="Open analytic tail; not an all-shift order-five theorem.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order5_m100_prefix_certificate",
        "date": "2026-07-13",
        "status": (
            "rigorous finite lambda=-100 contiguous order-five prefix through n=316"
        ),
        "proof_boundary": (
            "This artifact proves H_(5,n)(-100)>0 only for 0<=n<=316. "
            "It does not prove the analytic tail, all-shift order-five entry, "
            "forward invariance, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [
            source.relative_to(REPO_ROOT).as_posix() for source in SOURCE_PATHS
        ]
        + [
            "outputs/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.md",
            "outputs/formal_core.md",
        ],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order5_m100_prefix_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order5_m100_prefix_certificate.py"
        ),
        "exact": exact,
        "source_diagnostics": source_diagnostics,
        "finite": finite,
        "summary": {
            "rows": len(rows),
            "ready_to_apply_rows": sum(
                row.readiness == "ready_to_apply" for row in rows
            ),
            "coefficients": len(values),
            "positive_J_rows": len(finite["rows"]),
            "positive_relative_rows": len(finite["rows"]),
            "positive_H5_rows": len(finite["rows"]),
            "inconclusive_rows": 0,
            "open_analytic_tails": 1,
        },
        "rows": [asdict(row) for row in rows],
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    finite = artifact["finite"]
    lines = [
        "# Jensen-Window PF Compound Order-Five Lambda=-100 Prefix Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous finite `lambda=-100` contiguous order-five prefix",
        "through `n=316`. This is not a proof of the analytic tail, all-shift",
        "order five, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order5_m100_prefix_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order5_m100_prefix_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order5_m100_prefix_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order5_m100_prefix_certificate.py",
        "```",
        "",
        "## Stable Coordinate",
        "",
        "Use",
        "",
        "```text",
        exact["x_definition"],
        exact["G_definition"],
        exact["F_definition"],
        exact["J_definition"],
        "```",
        "",
        "Exact symbolic determinant algebra gives",
        "",
        "```text",
        exact["factorization"],
        exact["positive_scale"],
        "```",
        "",
        "Every factor in `W_n` is positive inside the already completed lower",
        "cones. Thus `H_(5,n)` has exactly the sign of `J_n`.",
        "",
        "## Arb Prefix",
        "",
        "Seven source files, merged in recorded precedence order, give 1024-bit",
        "outward-rounded balls for",
        "",
        "```text",
        "A_k(-100)>0 for every 0<=k<=324.",
        "```",
        "",
        "Direct Arb evaluation of the stable coordinate proves",
        "",
        "```text",
        "J_n(-100)>0 for every 0<=n<=316,",
        "relative_margin_n(-100)>0 for every 0<=n<=316,",
        "H_(5,n)(-100)>0 for every 0<=n<=316.",
        "```",
        "",
        "The weakest bounds occur at the final row:",
        "",
        "```text",
        f"minimum J row: n={finite['minimum_J_n']}",
        f"J_316={finite['minimum_J_ball']}",
        f"minimum relative row: n={finite['minimum_relative_n']}",
        f"relative_316={finite['minimum_relative_ball']}",
        "```",
        "",
        "The relative lower bound is about `0.006269`, leaving a visible strict",
        "buffer rather than a zero-containing endpoint enclosure.",
        "",
        "## Remaining Tail",
        "",
        "The sole endpoint-entry target is now",
        "",
        "```text",
        "J_n(-100)>0 for every n>=317,",
        "equivalently H_(5,n)(-100)>0 for every n>=317.",
        "```",
        "",
        "The compact-uniform eventual theorem proves positivity after some",
        "non-effective threshold, but it does not splice this explicit prefix to",
        "all shifts. That analytic tail remains open.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.md",
        "outputs/signed_hankel_jensen_bridge_target.md",
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
    with args.out.open("w", encoding="utf-8") as handle:
        json.dump(artifact, handle, indent=2, sort_keys=True)
        handle.write("\n")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-five lambda=-100 prefix: "
        f"{summary['coefficients']} coefficients, "
        f"{summary['positive_J_rows']} positive rows, "
        f"{summary['inconclusive_rows']} inconclusive"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
