#!/usr/bin/env python3
"""Certify the available lambda=-100 order-twelve partial prefix."""

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

import jensen_window_pf_endpoint_order10_counterexample as endpoint  # noqa: E402


SOURCE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_compound_order11_m100_prefix_certificate.json"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_compound_order12_m100_partial_prefix_certificate.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_compound_order12_m100_partial_prefix_certificate.md"
LAST_N = 1240


@dataclass(frozen=True)
class CertificateRow:
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


def quotient_sign(numerator_sign: str, denominator_sign: str) -> str:
    if numerator_sign not in {"positive", "negative"}:
        return "inconclusive"
    return "positive" if numerator_sign == denominator_sign else "negative"


def build_finite(source: dict) -> dict:
    finite = source.get("finite", {})
    source_rows = finite.get("rows", [])
    if [row.get("n") for row in source_rows] != list(range(1243)):
        raise RuntimeError("order-eleven source rows are not contiguous through n=1242")
    q11 = {int(row["n"]): flint.arb(row["Q11_ball"]) for row in source_rows}
    rows = []
    classes = {"positive": [], "negative": [], "inconclusive": []}
    for n in range(LAST_N + 1):
        numerator = q11[n + 1] ** 2 - q11[n] * q11[n + 2]
        numerator_sign = endpoint.sign_class(numerator)
        denominator_sign = "negative" if n + 2 <= 3 else "positive"
        q12_sign = quotient_sign(numerator_sign, denominator_sign)
        classes[q12_sign].append(n)
        rows.append(
            {
                "n": n,
                "Q10_denominator_index": n + 2,
                "Q10_denominator_sign": denominator_sign,
                "condensation_numerator_ball": str(numerator),
                "condensation_numerator_sign": numerator_sign,
                "Q12_sign": q12_sign,
            }
        )
    if classes["negative"] != [0, 1, 2, 3]:
        raise RuntimeError(f"unexpected Q12 negative set: {classes['negative']}")
    if classes["positive"] != list(range(4, LAST_N + 1)):
        raise RuntimeError("Q12 positive range is not contiguous from n=4")
    if classes["inconclusive"]:
        raise RuntimeError("Q12 partial prefix contains inconclusive rows")
    return {
        "lambda": "-100",
        "n_range": [0, LAST_N],
        "rows": rows,
        "negative_indices": classes["negative"],
        "positive_range": [4, LAST_N],
        "inconclusive_indices": classes["inconclusive"],
    }


def build_artifact() -> dict:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    if source.get("status") != "rigorous lambda=-100 signed order-eleven prefix through n=1242":
        raise RuntimeError("order-eleven prefix source changed")
    if source.get("summary", {}).get("positive_Q11_rows") != 1243:
        raise RuntimeError("order-eleven prefix source is not complete")
    if source.get("finite", {}).get("Q10_sign_map") != {
        "negative_indices": [0, 1, 2, 3],
        "positive_range": [4, 1244],
        "inconclusive_indices": [],
    }:
        raise RuntimeError("order-ten denominator sign map changed")
    finite = build_finite(source)
    exact = {
        "orientation": "epsilon_12=(-1)^66=1, Q_(12,n)=H_(12,n)",
        "condensation": "Q_(12,n)*Q_(10,n+2)=Q_(11,n+1)^2-Q_(11,n)*Q_(11,n+2)",
        "negative_prefix": "Q_(12,n)(-100)<0 for n=0,1,2,3",
        "positive_block": "Q_(12,n)(-100)>0 for every 4<=n<=1240",
        "remaining_finite_target": "Q_(12,n)(-100)>0 for every 1241<=n<=1492",
    }
    rows = [
        CertificateRow("co12m100ppc_01_condensation", "exact_identity", "ready_to_apply", "Signed order twelve is the order-eleven condensation numerator over the signed order-ten denominator.", exact["condensation"], "Exact signed Desnanot-Jacobi identity only."),
        CertificateRow("co12m100ppc_02_negative_prefix", "interval_theorem", "ready_to_apply", "The first four endpoint shifts are rigorously negative.", exact["negative_prefix"], "This forbids all-shift endpoint positivity at lambda=-100."),
        CertificateRow("co12m100ppc_03_positive_block", "interval_theorem", "ready_to_apply", "Every following available endpoint shift is rigorously positive.", exact["positive_block"], "Finite block only; the source ends before the analytic-tail handoff."),
        CertificateRow("co12m100ppc_04_remaining", "finite_theorem_target", "not_ready_to_apply", "Extend the positive block to meet the conditional analytic tail at n=1493.", exact["remaining_finite_target"], "Requires fresh coefficient enclosures through the enlarged collar."),
    ]
    return {
        "kind": "jensen_window_pf_compound_order12_m100_partial_prefix_certificate",
        "date": "2026-07-22",
        "status": "rigorous delayed order-twelve endpoint sign chart through n=1240 at lambda=-100",
        "proof_boundary": "This proves four negative and 1237 positive order-twelve endpoint signs only through n=1240. It does not prove the remaining finite block, analytic tail, lambda-zero completion, PF-infinity, RH, or Lambda<=0.",
        "source": {
            "path": SOURCE.relative_to(REPO_ROOT).as_posix(),
            "sha256": sha256(SOURCE),
            "kind": source.get("kind"),
            "status": source.get("status"),
        },
        "exact": exact,
        "finite": finite,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": 3,
            "open_rows": 1,
            "checked_shifts": LAST_N + 1,
            "negative_Q12_rows": 4,
            "positive_Q12_rows": LAST_N - 3,
            "inconclusive_Q12_rows": 0,
            "remaining_positive_shifts": 252,
            "all_shift_endpoint_theorems": 0,
            "rh_claims": 0,
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order12_m100_partial_prefix_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order12_m100_partial_prefix_certificate.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Order-Twelve Partial Endpoint Prefix",
        "",
        "Date: 2026-07-22",
        "",
        "Status: rigorous delayed endpoint sign chart through `n=1240` at",
        "`lambda=-100`. This is not all-shift endpoint positivity and is not a proof of RH or `Lambda<=0`.",
        "",
        "```text",
        exact["condensation"],
        "CERTIFIED:",
        exact["negative_prefix"],
        exact["positive_block"],
        "OPEN TARGET:",
        exact["remaining_finite_target"],
        "```",
        "",
        "All 1,241 derived signs are outward-rounded consequences of the",
        "hash-bound order-eleven prefix. There are zero inconclusive rows.",
        "The certified split is `4 negative, 1237 positive`.",
        "The four negative endpoint shifts are compatible with delayed heat",
        "recovery at lambda zero; they must not be promoted to positivity at",
        "lambda=-100. The remaining 252 positive shifts are not proved by this",
        "partial artifact; the separate endpoint-completion certificate now",
        "supplies them from fresh source coefficients.",
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
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print("wrote order-twelve partial endpoint prefix: 4 negative, 1237 positive, 0 inconclusive")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
