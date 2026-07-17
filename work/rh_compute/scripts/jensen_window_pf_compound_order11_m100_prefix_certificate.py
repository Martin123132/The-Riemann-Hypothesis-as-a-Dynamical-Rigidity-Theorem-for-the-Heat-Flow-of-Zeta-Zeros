#!/usr/bin/env python3
"""Certify the lambda=-100 signed order-eleven prefix through n=1242."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
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

import jensen_window_pf_compound_order9_m100_prefix_certificate as prefix  # noqa: E402
import jensen_window_pf_compound_order10_m100_finite_splice_certificate as order10  # noqa: E402
import jensen_window_pf_endpoint_order10_counterexample as endpoint  # noqa: E402


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_m100_prefix_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order11_m100_prefix_certificate.md"
)
ORDER10_SOURCE = order10.DEFAULT_OUT
PRECISION_REPAIR = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "negative_lambda_m100_order11_k134_k152_dps220.jsonl"
)
PRECISION_REPAIR_SUMMARY = PRECISION_REPAIR.with_name(
    "negative_lambda_m100_order11_k134_k152_dps220_summary.json"
)
ENDPOINT_EXTENSION = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "negative_lambda_m100_order11_k1261_k1262_dps220.jsonl"
)
ENDPOINT_EXTENSION_SUMMARY = ENDPOINT_EXTENSION.with_name(
    "negative_lambda_m100_order11_k1261_k1262_dps220_summary.json"
)
SOURCE_PATHS = (*order10.SOURCE_PATHS, PRECISION_REPAIR, ENDPOINT_EXTENSION)
PREFIX_LAST_N = 1242
MAX_COEFFICIENT_INDEX = 1262
DIRECT_LAST_N = 4
PRECISION_BITS = 4096


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


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT).as_posix()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def validate_enclosure_summary(
    path: Path,
    *,
    rows: int,
    k_min: int,
    k_max: int,
) -> dict:
    summary = load_json(path)
    expected = {
        "kind": "acb_coefficient_enclosure_summary",
        "rows": rows,
        "k_min": k_min,
        "k_max": k_max,
        "lambdas": ["-100.0"],
        "n_sum": 70,
        "cutoff": "7",
        "dps": 220,
        "abs_tol": "1e-150",
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            raise RuntimeError(
                f"coefficient summary {path.name} changed at {key}: "
                f"{summary.get(key)!r}"
            )
    return summary


def validate_sources() -> dict:
    old = load_json(ORDER10_SOURCE)
    finite = old.get("finite", {})
    if (
        old.get("kind")
        != "jensen_window_pf_compound_order10_m100_finite_splice_certificate"
        or finite.get("combined_positive_range") != [4, 1242]
        or finite.get("preserved_negative_indices") != [0, 1, 2, 3]
    ):
        raise RuntimeError("order-ten endpoint source contract changed")
    repair = validate_enclosure_summary(
        PRECISION_REPAIR_SUMMARY,
        rows=19,
        k_min=134,
        k_max=152,
    )
    extension = validate_enclosure_summary(
        ENDPOINT_EXTENSION_SUMMARY,
        rows=2,
        k_min=1261,
        k_max=1262,
    )
    return {
        "order10_endpoint": {
            "path": relative(ORDER10_SOURCE),
            "sha256": sha256(ORDER10_SOURCE),
            "claim": (
                "Q_(10,n)(-100)<0 for 0<=n<=3 and "
                "Q_(10,n)(-100)>0 for 4<=n<=1242"
            ),
        },
        "precision_repair": {
            "path": relative(PRECISION_REPAIR),
            "sha256": sha256(PRECISION_REPAIR),
            "summary_path": relative(PRECISION_REPAIR_SUMMARY),
            "summary_sha256": sha256(PRECISION_REPAIR_SUMMARY),
            "summary": repair,
        },
        "endpoint_extension": {
            "path": relative(ENDPOINT_EXTENSION),
            "sha256": sha256(ENDPOINT_EXTENSION),
            "summary_path": relative(ENDPOINT_EXTENSION_SUMMARY),
            "summary_sha256": sha256(ENDPOINT_EXTENSION_SUMMARY),
            "summary": extension,
        },
    }


def load_coefficients() -> tuple[dict[int, flint.arb], list[dict]]:
    values, diagnostics = prefix.merged_coefficients(
        SOURCE_PATHS,
        MAX_COEFFICIENT_INDEX,
    )
    return values, diagnostics


def _require_positive(name: str, family: dict[int, flint.arb]) -> None:
    failed = [
        index
        for index, value in family.items()
        if endpoint.sign_class(value) != "positive"
    ]
    if failed:
        raise RuntimeError(f"{name} lost strict positivity at {failed[:10]}")


def _next_signed_layer(
    denominator: dict[int, flint.arb],
    previous: dict[int, flint.arb],
    count: int,
    *,
    require_positive_margin: bool,
) -> tuple[dict[int, flint.arb], dict[int, flint.arb]]:
    values = {}
    margins = {}
    for index in range(count):
        margin = (
            previous[index + 1] ** 2
            / (previous[index] * previous[index + 2])
            - 1
        )
        if require_positive_margin and endpoint.sign_class(margin) != "positive":
            raise RuntimeError(f"signed layer margin failed at n={index}")
        margins[index] = margin
        values[index] = (
            previous[index]
            * previous[index + 2]
            * margin
            / denominator[index + 2]
        )
    return values, margins


def stable_order11_prefix(values: dict[int, flint.arb]) -> dict:
    """Rebuild Q6 through Q11 without subtracting raw Hankel determinants."""
    maximum = MAX_COEFFICIENT_INDEX
    contractions = {
        index: values[index - 1] * values[index + 1] / values[index] ** 2
        for index in range(1, maximum)
    }
    defects = {index: 1 - value for index, value in contractions.items()}
    gaps = {
        index: (
            defects[index + 2] ** 2
            - contractions[index + 2] ** 2
            * defects[index + 1]
            * defects[index + 3]
        )
        for index in range(maximum - 3)
    }
    order4_margins = {
        index: (
            gaps[index + 1] ** 2
            - contractions[index + 3] ** 3 * gaps[index] * gaps[index + 2]
        )
        for index in range(maximum - 5)
    }
    stable_h5 = {
        index: (
            defects[index + 3]
            * defects[index + 5]
            * order4_margins[index + 1] ** 2
            - contractions[index + 4] ** 4
            * defects[index + 4] ** 2
            * order4_margins[index]
            * order4_margins[index + 2]
        )
        for index in range(maximum - 7)
    }
    for name, family in (
        ("contraction", contractions),
        ("defect", defects),
        ("order3_gap", gaps),
        ("order4_margin", order4_margins),
        ("order5_stable_margin", stable_h5),
    ):
        _require_positive(name, family)

    def h5_value(index: int) -> flint.arb:
        ratio = values[index + 1] / values[index]
        return (
            values[index] ** 5
            * ratio**20
            * contractions[index + 1] ** 15
            * contractions[index + 2] ** 10
            * contractions[index + 3] ** 5
            * stable_h5[index]
            / (
                defects[index + 3]
                * defects[index + 4] ** 2
                * defects[index + 5]
                * gaps[index + 2]
            )
        )

    def h4_value(index: int) -> flint.arb:
        ratio = values[index + 1] / values[index]
        return (
            values[index] ** 4
            * ratio**12
            * contractions[index + 1] ** 8
            * contractions[index + 2] ** 4
            * order4_margins[index]
            / defects[index + 3]
        )

    h5_values = {
        index: h5_value(index) for index in range(maximum - 7)
    }
    h4_values = {
        index: h4_value(index) for index in range(maximum - 7)
    }
    _require_positive("H4", h4_values)
    _require_positive("H5", h5_values)

    q6, q6_margins = _next_signed_layer(
        h4_values, h5_values, maximum - 9, require_positive_margin=True
    )
    q7, q7_margins = _next_signed_layer(
        h5_values, q6, maximum - 11, require_positive_margin=True
    )
    q8, q8_margins = _next_signed_layer(
        q6, q7, maximum - 13, require_positive_margin=True
    )
    q9, q9_margins = _next_signed_layer(
        q7, q8, maximum - 15, require_positive_margin=True
    )
    q10, q10_margins = _next_signed_layer(
        q8, q9, maximum - 17, require_positive_margin=False
    )
    for name, family in (("Q6", q6), ("Q7", q7), ("Q8", q8), ("Q9", q9)):
        _require_positive(name, family)

    q10_classes = {
        kind: [
            index
            for index, value in q10.items()
            if endpoint.sign_class(value) == kind
        ]
        for kind in ("positive", "negative", "inconclusive")
    }
    if q10_classes["negative"] != [0, 1, 2, 3]:
        raise RuntimeError(f"unexpected Q10 negative set: {q10_classes['negative']}")
    if q10_classes["inconclusive"]:
        raise RuntimeError("Q10 reconstruction contains inconclusive rows")
    if q10_classes["positive"] != list(range(4, maximum - 17)):
        raise RuntimeError("Q10 positive range is not contiguous from n=4")

    q11 = {}
    numerators = {}
    relative_q10 = {}
    rows = []
    minimum: tuple[flint.arb, int] | None = None
    for index in range(PREFIX_LAST_N + 1):
        numerator = q10[index + 1] ** 2 - q10[index] * q10[index + 2]
        if endpoint.sign_class(numerator) != "positive":
            raise RuntimeError(f"order-eleven numerator failed at n={index}")
        value = numerator / q9[index + 2]
        if endpoint.sign_class(value) != "positive":
            raise RuntimeError(f"Q11 failed at n={index}")
        numerators[index] = numerator
        q11[index] = value
        relative_text = None
        relative_lower = None
        if index >= 4:
            relative = (
                q10[index + 1] ** 2 / (q10[index] * q10[index + 2]) - 1
            )
            if endpoint.sign_class(relative) != "positive":
                raise RuntimeError(f"relative Q10 margin failed at n={index}")
            relative_q10[index] = relative
            relative_text = prefix.arb_text(relative, 80)
            relative_lower = prefix.arb_lower_text(relative, 80)
            if minimum is None or relative.lower() < minimum[0].lower():
                minimum = (relative, index)
        rows.append(
            {
                "n": index,
                "Q10_signs": [
                    endpoint.sign_class(q10[index + shift])
                    for shift in range(3)
                ],
                "condensation_numerator_ball": prefix.arb_text(numerator, 80),
                "condensation_numerator_lower": prefix.arb_lower_text(numerator, 80),
                "relative_Q10_margin_ball": relative_text,
                "relative_Q10_margin_lower": relative_lower,
                "Q11_ball": prefix.arb_text(value, 80),
                "Q11_lower": prefix.arb_lower_text(value, 80),
                "Q11_sign": "positive_by_signed_condensation",
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
            "order4_margins": len(order4_margins),
            "order5_stable_margins": len(stable_h5),
            "H4_values": len(h4_values),
            "H5_values": len(h5_values),
            "Q6_values": len(q6),
            "Q7_values": len(q7),
            "Q8_values": len(q8),
            "Q9_values": len(q9),
            "Q11_numerators": len(numerators),
            "Q11_values": len(q11),
        },
        "Q10_sign_map": {
            "negative_indices": q10_classes["negative"],
            "positive_range": [4, max(q10)],
            "inconclusive_indices": q10_classes["inconclusive"],
        },
        "exceptional_Q11_indices": [0, 1, 2, 3],
        "positive_cone_Q11_range": [4, PREFIX_LAST_N],
        "minimum_relative_n": minimum[1],
        "minimum_relative_ball": prefix.arb_text(minimum[0], 80),
        "minimum_relative_lower": prefix.arb_lower_text(minimum[0], 80),
        "internal": {
            "q9": q9,
            "q10": q10,
            "q11": q11,
            "q10_margins": q10_margins,
            "q6_margins": q6_margins,
            "q7_margins": q7_margins,
            "q8_margins": q8_margins,
            "q9_margins": q9_margins,
        },
    }


def direct_audit(
    values: dict[int, flint.arb],
    stable: dict,
) -> list[dict]:
    rows = []
    internal = stable["internal"]
    for shift in range(DIRECT_LAST_N + 1):
        q9 = endpoint.direct_signed_hankel(values, 9, shift + 2)
        q10_left = endpoint.direct_signed_hankel(values, 10, shift)
        q10_center = endpoint.direct_signed_hankel(values, 10, shift + 1)
        q10_right = endpoint.direct_signed_hankel(values, 10, shift + 2)
        q11 = endpoint.direct_signed_hankel(values, 11, shift)
        schur = endpoint.direct_deep_schur(values, 11, shift + 10)
        numerator = q10_center**2 - q10_left * q10_right
        toda_residual = q11 * q9 - numerator
        schur_residual = schur - q11 / values[0] ** 11
        checks = {
            "q9_positive": endpoint.sign_class(q9) == "positive",
            "numerator_positive": endpoint.sign_class(numerator) == "positive",
            "q11_positive": endpoint.sign_class(q11) == "positive",
            "deep_schur_positive": endpoint.sign_class(schur) == "positive",
            "toda_residual_contains_zero": bool(toda_residual.contains(0)),
            "schur_residual_contains_zero": bool(schur_residual.contains(0)),
            "stable_q11_overlaps_direct": bool(internal["q11"][shift].overlaps(q11)),
        }
        if not all(checks.values()):
            raise RuntimeError(f"direct Q11 audit failed at n={shift}: {checks}")
        rows.append(
            {
                "n": shift,
                "N": shift + 10,
                "shape": f"(({shift + 10}^11))",
                "Q9_ball": endpoint.ball_text(q9),
                "Q10_left_ball": endpoint.ball_text(q10_left),
                "Q10_center_ball": endpoint.ball_text(q10_center),
                "Q10_right_ball": endpoint.ball_text(q10_right),
                "condensation_numerator_ball": endpoint.ball_text(numerator),
                "Q11_ball": endpoint.ball_text(q11),
                "deep_schur_ball": endpoint.ball_text(schur),
                "toda_residual_ball": endpoint.ball_text(toda_residual),
                "schur_residual_ball": endpoint.ball_text(schur_residual),
                "checks": checks,
            }
        )
    return rows


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    sources = validate_sources()
    values, source_diagnostics = load_coefficients()
    stable = stable_order11_prefix(values)
    direct = direct_audit(values, stable)
    stable.pop("internal")
    exact = {
        "orientation": "epsilon_11=(-1)^55=-1, Q_(11,n)=-H_(11,n)",
        "condensation": (
            "Q_(11,n)*Q_(9,n+2)=Q_(10,n+1)^2-"
            "Q_(10,n)*Q_(10,n+2)"
        ),
        "relative_coordinate": (
            "M_n=Q_(10,n+1)^2/(Q_(10,n)*Q_(10,n+2))-1"
        ),
        "exceptional_prefix": (
            "Q_(10,n+1)^2-Q_(10,n)*Q_(10,n+2)>0 for 0<=n<=3"
        ),
        "positive_cone_block": "M_n>0 for every 4<=n<=1242",
        "finite_prefix": "Q_(11,n)(-100)>0 for every 0<=n<=1242",
        "analytic_handoff": "prove Q_(11,n)(-100)>0 for every n>=1243",
    }
    rows = [
        PrefixRow(
            "co11m100pc_01_condensation",
            "exact_identity",
            "ready_to_apply",
            "Signed order eleven is the order-ten condensation numerator over a positive order-nine denominator.",
            exact["condensation"],
            "Exact signed Desnanot-Jacobi identity only.",
        ),
        PrefixRow(
            "co11m100pc_02_coefficients",
            "interval_input",
            "ready_to_apply",
            "Merged retained-integral Arb enclosures cover every coefficient through A_1262.",
            "A_k(-100)>0 for every 0<=k<=1262",
            "Finite endpoint coefficient range only.",
            {"sources": source_diagnostics},
        ),
        PrefixRow(
            "co11m100pc_03_precision_repair",
            "interval_certificate",
            "ready_to_apply",
            "A contiguous nineteen-row repair removes the sole inherited precision seam.",
            "A_k(-100) re-enclosed at 220 dps for 134<=k<=152",
            "Precision repair only; it does not infer signs by numerical approximation.",
        ),
        PrefixRow(
            "co11m100pc_04_stable_chain",
            "interval_certificate",
            "ready_to_apply",
            "Cancellation-preserving H4/H5/Q6/Q7/Q8/Q9/Q10 coordinates rebuild every finite numerator.",
            exact["exceptional_prefix"] + "; " + exact["positive_cone_block"],
            "Outward-rounded finite Arb arithmetic only.",
            stable["positive_coordinate_counts"],
        ),
        PrefixRow(
            "co11m100pc_05_direct_audit",
            "interval_certificate",
            "ready_to_apply",
            "Independent 11x11 Hankel and Jacobi-Trudi determinants verify the exceptional prefix and first positive-cone row.",
            "Q_(11,n)(-100)>0 for 0<=n<=4",
            "Five direct low-index audits only.",
            {"rows": direct},
        ),
        PrefixRow(
            "co11m100pc_06_finite_prefix",
            "interval_theorem",
            "ready_to_apply",
            "Every finite order-eleven condensation numerator is strictly positive with no unresolved row.",
            exact["finite_prefix"],
            "Finite prefix only; no analytic tail or heat propagation is claimed.",
        ),
        PrefixRow(
            "co11m100pc_07_handoff",
            "open_handoff",
            "not_ready_to_apply",
            "The finite endpoint theorem now meets the proposed analytic tail at the next integer.",
            exact["analytic_handoff"],
            "The order-eleven continuous curvature theorem remains open.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order11_m100_prefix_certificate",
        "date": "2026-07-16",
        "status": "rigorous lambda=-100 signed order-eleven prefix through n=1242",
        "proof_boundary": (
            "This artifact proves Q_(11,n)(-100)>0 only for 0<=n<=1242. "
            "It does not prove the analytic order-eleven tail, heat propagation, "
            "PF-infinity, RH, or Lambda<=0."
        ),
        "sources": sources,
        "source_diagnostics": source_diagnostics,
        "exact": exact,
        "finite": stable,
        "direct_audit": direct,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": sum(row.readiness == "ready_to_apply" for row in rows),
            "open_rows": sum(row.readiness != "ready_to_apply" for row in rows),
            "coefficients": len(values),
            "precision_repair_rows": 19,
            "endpoint_extension_rows": 2,
            "Q11_numerators": PREFIX_LAST_N + 1,
            "positive_Q11_rows": PREFIX_LAST_N + 1,
            "negative_Q11_rows": 0,
            "inconclusive_Q11_rows": 0,
            "exceptional_Q11_rows": 4,
            "positive_cone_Q11_rows": PREFIX_LAST_N - 3,
            "direct_Q11_checks": len(direct),
            "finite_prefix_theorems": 1,
            "open_analytic_tail_targets": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order11_m100_prefix_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order11_m100_prefix_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    finite = artifact["finite"]
    exact = artifact["exact"]
    lines = [
        "# Order-Eleven Lambda=-100 Prefix Certificate",
        "",
        "Date: 2026-07-16",
        "",
        "Status: rigorous finite endpoint theorem. This is not a proof of the",
        "analytic order-eleven tail, heat propagation, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order11_m100_prefix_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order11_m100_prefix_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order11_m100_prefix_certificate.py",
        "```",
        "",
        "## Exact Reduction",
        "",
        "```text",
        exact["orientation"],
        exact["condensation"],
        exact["relative_coordinate"],
        "```",
        "",
        "The four rows `n=0,1,2,3` lie outside the positive `Q10` cone, so",
        "their condensation numerators are enclosed directly through the stable",
        "chain and audited with independent `11x11` determinants. For",
        "`4<=n<=1242`, the relative `Q10` margin is rigorously positive.",
        "",
        "## Finite Theorem",
        "",
        "```text",
        exact["exceptional_prefix"],
        exact["positive_cone_block"],
        exact["finite_prefix"],
        "minimum relative margin at n="
        + str(finite["minimum_relative_n"])
        + ": "
        + finite["minimum_relative_ball"],
        "```",
        "",
        "There are `1243` positive rows, zero negative rows, and zero",
        "inconclusive rows. The remaining endpoint handoff begins exactly at",
        "`n=1243` and requires the continuous order-eleven curvature theorem.",
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
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-eleven finite prefix: "
        f"{summary['coefficients']} coefficients, "
        f"{summary['positive_Q11_rows']} positive Q11 rows, "
        f"{summary['inconclusive_Q11_rows']} inconclusive"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
