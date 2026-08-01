#!/usr/bin/env python3
"""Certify the finite lambda=-100 order-twelve endpoint collar through n=1492."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
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

import jensen_window_pf_compound_order11_m100_prefix_certificate as order11  # noqa: E402
import jensen_window_pf_compound_order9_m100_prefix_certificate as prefix  # noqa: E402
import jensen_window_pf_endpoint_order10_counterexample as endpoint  # noqa: E402


RESULTS = REPO_ROOT / "work/rh_compute/results"
ORDER11_PREFIX = RESULTS / "jensen_window_pf_compound_order11_m100_prefix_certificate.json"
ORDER12_PARTIAL = RESULTS / "jensen_window_pf_compound_order12_m100_partial_prefix_certificate.json"
DEFAULT_OUT = RESULTS / "jensen_window_pf_compound_order12_m100_endpoint_completion_certificate.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_compound_order12_m100_endpoint_completion_certificate.md"

CHUNK_SPECS = (
    (1263, 1325, "negative_lambda_m100_order12_k1263_k1325_dps220"),
    (1326, 1388, "negative_lambda_m100_order12_k1326_k1388_dps220"),
    (1389, 1451, "negative_lambda_m100_order12_k1389_k1451_dps220"),
    (1452, 1514, "negative_lambda_m100_order12_k1452_k1514_dps220"),
)
CHUNK_PATHS = tuple(RESULTS / f"{run_id}.jsonl" for _, _, run_id in CHUNK_SPECS)
SOURCE_PATHS = (*order11.SOURCE_PATHS, *CHUNK_PATHS)
MAX_COEFFICIENT_INDEX = 1514
Q11_LAST_N = 1494
COLLAR_FIRST_N = 1241
ENDPOINT_FIRST_N = 1241
ENDPOINT_LAST_N = 1492


@dataclass(frozen=True)
class CompletionRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT).as_posix()


def source_record(path: Path, artifact: dict) -> dict:
    return {
        "path": relative(path),
        "sha256": order11.sha256(path),
        "kind": artifact.get("kind"),
        "status": artifact.get("status"),
    }


def validate_chunk(k_min: int, k_max: int, run_id: str) -> dict:
    path = RESULTS / f"{run_id}.jsonl"
    summary_path = RESULTS / f"{run_id}_summary.json"
    summary = load_json(summary_path)
    expected_summary = {
        "kind": "acb_coefficient_enclosure_summary",
        "rows": k_max - k_min + 1,
        "k_min": k_min,
        "k_max": k_max,
        "lambdas": ["-100.0"],
        "n_sum": 70,
        "cutoff": "7",
        "dps": 220,
        "abs_tol": "1e-150",
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            raise RuntimeError(f"chunk {run_id} summary changed at {key}")

    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    if [row.get("k") for row in rows] != list(range(k_min, k_max + 1)):
        raise RuntimeError(f"chunk {run_id} rows are not contiguous")
    for row in rows:
        expected = {
            "kind": "acb_coefficient_enclosure",
            "lam": "-100.0",
            "n_sum": 70,
            "cutoff": "7",
            "dps": 220,
            "abs_tol": "1e-150",
        }
        for key, value in expected.items():
            if row.get(key) != value:
                raise RuntimeError(f"chunk {run_id} row {row.get('k')} changed at {key}")
        if endpoint.sign_class(flint.arb(row["A_ball"])) != "positive":
            raise RuntimeError(f"chunk {run_id} has nonpositive A_{row.get('k')}")
    return {
        "path": relative(path),
        "sha256": order11.sha256(path),
        "summary_path": relative(summary_path),
        "summary_sha256": order11.sha256(summary_path),
        "summary": summary,
    }


def validate_sources() -> dict:
    inherited = load_json(ORDER11_PREFIX)
    partial = load_json(ORDER12_PARTIAL)
    inherited_finite = inherited.get("finite", {})
    if (
        inherited.get("status")
        != "rigorous lambda=-100 signed order-eleven prefix through n=1242"
        or inherited_finite.get("n_range") != [0, 1242]
        or inherited.get("summary", {}).get("positive_Q11_rows") != 1243
    ):
        raise RuntimeError("inherited order-eleven prefix contract changed")
    partial_finite = partial.get("finite", {})
    if (
        partial.get("status")
        != "rigorous delayed order-twelve endpoint sign chart through n=1240 at lambda=-100"
        or partial_finite.get("negative_indices") != [0, 1, 2, 3]
        or partial_finite.get("positive_range") != [4, 1240]
        or partial.get("summary", {}).get("inconclusive_Q12_rows") != 0
    ):
        raise RuntimeError("inherited order-twelve partial contract changed")
    chunks = [validate_chunk(*spec) for spec in CHUNK_SPECS]
    return {
        "inherited_order11_prefix": source_record(ORDER11_PREFIX, inherited),
        "inherited_order12_partial": source_record(ORDER12_PARTIAL, partial),
        "coefficient_chunks": chunks,
    }


def build_artifact() -> dict:
    flint.ctx.prec = order11.PRECISION_BITS
    sources = validate_sources()
    values, source_diagnostics = prefix.merged_coefficients(
        SOURCE_PATHS,
        MAX_COEFFICIENT_INDEX,
    )
    stable = order11.stable_order11_prefix(
        values,
        maximum=MAX_COEFFICIENT_INDEX,
        prefix_last_n=Q11_LAST_N,
    )
    internal = stable["internal"]
    q10 = internal["q10"]
    q11 = internal["q11"]

    inherited = load_json(ORDER11_PREFIX)
    inherited_rows = inherited["finite"]["rows"]
    if [row.get("n") for row in inherited_rows] != list(range(1243)):
        raise RuntimeError("inherited Q11 rows lost contiguity")
    overlap_failures = []
    for row in inherited_rows:
        n = int(row["n"])
        if not q11[n].overlaps(flint.arb(row["Q11_ball"])):
            overlap_failures.append(n)
    if overlap_failures:
        raise RuntimeError(f"enlarged Q11 rebuild lost overlap at {overlap_failures[:10]}")

    collar_rows = [
        {
            "n": n,
            "Q11_ball": prefix.arb_text(q11[n], 80),
            "Q11_lower": prefix.arb_lower_text(q11[n], 80),
            "Q11_sign": endpoint.sign_class(q11[n]),
        }
        for n in range(COLLAR_FIRST_N, Q11_LAST_N + 1)
    ]
    if any(row["Q11_sign"] != "positive" for row in collar_rows):
        raise RuntimeError("order-eleven collar contains a nonpositive row")

    endpoint_rows = []
    minimum_relative: tuple[flint.arb, int] | None = None
    minimum_q12: tuple[flint.arb, int] | None = None
    for n in range(ENDPOINT_FIRST_N, ENDPOINT_LAST_N + 1):
        denominator = q10[n + 2]
        relative_margin = q11[n + 1] ** 2 / (q11[n] * q11[n + 2]) - 1
        raw_numerator = q11[n + 1] ** 2 - q11[n] * q11[n + 2]
        factored_numerator = q11[n] * q11[n + 2] * relative_margin
        q12 = factored_numerator / denominator
        checks = {
            "Q10_denominator_positive": endpoint.sign_class(denominator) == "positive",
            "relative_Q11_margin_positive": endpoint.sign_class(relative_margin) == "positive",
            "raw_numerator_positive": endpoint.sign_class(raw_numerator) == "positive",
            "factored_numerator_positive": endpoint.sign_class(factored_numerator) == "positive",
            "raw_and_factored_overlap": bool(raw_numerator.overlaps(factored_numerator)),
            "Q12_positive": endpoint.sign_class(q12) == "positive",
        }
        if not all(checks.values()):
            failed = [name for name, passed in checks.items() if not passed]
            raise RuntimeError(f"order-twelve endpoint failed at n={n}: {failed}")
        if minimum_relative is None or relative_margin.lower() < minimum_relative[0].lower():
            minimum_relative = (relative_margin, n)
        if minimum_q12 is None or q12.lower() < minimum_q12[0].lower():
            minimum_q12 = (q12, n)
        endpoint_rows.append(
            {
                "n": n,
                "Q10_denominator_index": n + 2,
                "Q10_denominator_sign": "positive",
                "relative_Q11_margin_ball": prefix.arb_text(relative_margin, 80),
                "relative_Q11_margin_lower": prefix.arb_lower_text(relative_margin, 80),
                "raw_numerator_ball": prefix.arb_text(raw_numerator, 80),
                "factored_numerator_ball": prefix.arb_text(factored_numerator, 80),
                "Q12_ball": prefix.arb_text(q12, 80),
                "Q12_lower": prefix.arb_lower_text(q12, 80),
                "Q12_sign": "positive",
                "checks": checks,
            }
        )
    assert minimum_relative is not None and minimum_q12 is not None

    exact = {
        "condensation": "Q_(12,n)*Q_(10,n+2)=Q_(11,n+1)^2-Q_(11,n)*Q_(11,n+2)",
        "relative_coordinate": "R_n=Q_(11,n+1)^2/(Q_(11,n)*Q_(11,n+2))-1",
        "inherited_signs": "Q_(12,n)(-100)<0 for n=0,1,2,3 and Q_(12,n)(-100)>0 for every 4<=n<=1240",
        "collar_theorem": "Q_(12,n)(-100)>0 for every 1241<=n<=1492",
        "finite_sign_chart": "Q_(12,n)(-100)<0 for n=0,1,2,3 and Q_(12,n)(-100)>0 for every 4<=n<=1492",
        "remaining_tail_target": "Q_(12,n)(-100)>0 for every n>=1493",
    }
    rows = [
        CompletionRow(
            "co12m100ecc_01_condensation",
            "exact_identity",
            "ready_to_apply",
            "Signed order twelve is the order-eleven log-concavity numerator over the signed order-ten denominator.",
            exact["condensation"],
            "Exact signed Desnanot-Jacobi identity only.",
        ),
        CompletionRow(
            "co12m100ecc_02_coefficients",
            "interval_input",
            "ready_to_apply",
            "Four hash-bound retained-integral chunks extend rigorous coefficient coverage through A_1514.",
            "A_k(-100)>0 for every 0<=k<=1514",
            "Finite lambda=-100 coefficient range only.",
            {"chunks": len(CHUNK_SPECS), "new_coefficients": 252},
        ),
        CompletionRow(
            "co12m100ecc_03_q11_collar",
            "interval_certificate",
            "ready_to_apply",
            "The stable H4/H5/Q6/.../Q11 chain rebuilds the enlarged collar with strict signs.",
            "Q_(11,n)(-100)>0 for every 0<=n<=1494",
            "Finite order-eleven reconstruction used only as input to the order-twelve collar.",
        ),
        CompletionRow(
            "co12m100ecc_04_overlap",
            "reproduction_gate",
            "ready_to_apply",
            "Every inherited order-eleven ball overlaps the enlarged reconstruction.",
            "overlap(Q11_old(n),Q11_new(n)) for every 0<=n<=1242",
            "Agreement gate only; it does not infer an analytic tail.",
            {"overlap_rows": len(inherited_rows), "failures": overlap_failures},
        ),
        CompletionRow(
            "co12m100ecc_05_endpoint_collar",
            "interval_theorem",
            "ready_to_apply",
            "All 252 remaining finite order-twelve endpoint rows have strictly positive raw and factored condensation numerators.",
            exact["collar_theorem"],
            "Finite endpoint collar only.",
        ),
        CompletionRow(
            "co12m100ecc_06_finite_composition",
            "exact_theorem_composition",
            "ready_to_apply",
            "The inherited sign chart and new collar meet without an index gap.",
            exact["finite_sign_chart"],
            "The first four endpoint signs remain negative and must not be promoted.",
        ),
        CompletionRow(
            "co12m100ecc_07_tail",
            "theorem_target",
            "not_ready_to_apply",
            "The conditional curvature bridge begins at the next shift.",
            exact["remaining_tail_target"],
            "Requires the separate v_1 curvature theorem; no tail sign is inferred here.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order12_m100_endpoint_completion_certificate",
        "date": "2026-07-22",
        "status": "rigorous lambda=-100 order-twelve finite endpoint sign chart through n=1492",
        "proof_boundary": (
            "This artifact proves the delayed finite order-twelve sign chart at lambda=-100: "
            "four negative rows followed by positivity through n=1492. It does not prove the "
            "analytic tail, all-shift endpoint positivity, PF-infinity, Lambda<=0, or RH."
        ),
        "sources": sources,
        "source_diagnostics": source_diagnostics,
        "exact": exact,
        "order11_collar": {
            "n_range": [COLLAR_FIRST_N, Q11_LAST_N],
            "coefficient_range": [0, MAX_COEFFICIENT_INDEX],
            "precision_bits": order11.PRECISION_BITS,
            "rows": collar_rows,
            "inherited_overlap_rows": len(inherited_rows),
            "inherited_overlap_failures": overlap_failures,
            "Q10_sign_map": stable["Q10_sign_map"],
        },
        "finite": {
            "lambda": "-100",
            "n_range": [ENDPOINT_FIRST_N, ENDPOINT_LAST_N],
            "rows": endpoint_rows,
            "minimum_relative_n": minimum_relative[1],
            "minimum_relative_ball": prefix.arb_text(minimum_relative[0], 80),
            "minimum_relative_lower": prefix.arb_lower_text(minimum_relative[0], 80),
            "minimum_Q12_n": minimum_q12[1],
            "minimum_Q12_ball": prefix.arb_text(minimum_q12[0], 80),
            "minimum_Q12_lower": prefix.arb_lower_text(minimum_q12[0], 80),
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": 6,
            "open_rows": 1,
            "coefficient_rows": MAX_COEFFICIENT_INDEX + 1,
            "new_coefficient_rows": 252,
            "coefficient_chunks": len(CHUNK_SPECS),
            "rebuilt_Q11_rows": Q11_LAST_N + 1,
            "inherited_Q11_overlap_rows": len(inherited_rows),
            "Q11_collar_rows": len(collar_rows),
            "endpoint_collar_rows": len(endpoint_rows),
            "positive_endpoint_collar_rows": len(endpoint_rows),
            "negative_endpoint_collar_rows": 0,
            "inconclusive_endpoint_collar_rows": 0,
            "finite_sign_chart_theorems": 1,
            "open_analytic_tail_targets": 1,
            "all_shift_endpoint_positivity_theorems": 0,
            "rh_claims": 0,
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order12_m100_endpoint_completion_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order12_m100_endpoint_completion_certificate.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    finite = artifact["finite"]
    lines = [
        "# Order-Twelve Lambda=-100 Endpoint Completion",
        "",
        "Date: 2026-07-22",
        "",
        "Status: rigorous finite endpoint handoff through `n=1492`. This is not",
        "all-shift endpoint positivity or PF-infinity, and it is not a proof of RH or `Lambda<=0`.",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order12_m100_endpoint_completion_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order12_m100_endpoint_completion_certificate.py",
        "```",
        "",
        "## Exact Reduction",
        "",
        "```text",
        exact["condensation"],
        exact["relative_coordinate"],
        "```",
        "",
        "Four retained-integral chunks add exactly 252 coefficients, extending",
        "coverage from `A_1262` through `A_1514`. The enlarged stable chain",
        "overlaps every one of the 1,243 inherited order-eleven balls.",
        "",
        "## Finite Theorem",
        "",
        "```text",
        exact["collar_theorem"],
        exact["finite_sign_chart"],
        "minimum relative Q11 margin at n=" + str(finite["minimum_relative_n"]) + ": " + finite["minimum_relative_ball"],
        "```",
        "",
        "All 252 new endpoint rows have positive raw and factored numerators,",
        "positive denominators, positive Q12 balls, and zero inconclusive rows.",
        "The first four endpoint rows remain rigorously negative. The sole",
        "remaining endpoint obligation is the analytic tail `n>=1493`, supplied",
        "conditionally by the separate ninth-coordinate curvature target.",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-twelve endpoint completion: "
        f"{summary['endpoint_collar_rows']} positive collar rows, "
        f"{summary['inconclusive_endpoint_collar_rows']} inconclusive"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
