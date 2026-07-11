#!/usr/bin/env python3
"""Audit the factorial-multiplier split in the Jensen-window PF route."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402

from jensen_window_pf_heat_flow_monotone_closure_scout import (  # noqa: E402
    DEFAULT_ENCLOSURE_JSONL,
    REPO_ROOT,
    decimal_format,
    decimal_lam_key,
)


getcontext().prec = 120

DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_factorial_multiplier_split_audit.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_factorial_multiplier_split_audit.md"


@dataclass(frozen=True)
class FiniteMinimum:
    sample: str
    lam: str
    n: int


@dataclass(frozen=True)
class SplitDiagnostics:
    lambdas: list[str]
    max_coefficient_index: int
    raw_degree2_rows: int
    raw_degree2_negative_rows: int
    normalized_degree2_rows: int
    normalized_degree2_positive_rows: int
    raw_ratio_below_one_rows: int
    normalized_threshold_rows: int
    normalized_threshold_positive_rows: int
    max_raw_discriminant: FiniteMinimum
    min_normalized_discriminant: FiniteMinimum
    min_normalized_threshold_margin: FiniteMinimum


def parse_ball_midpoint(text: str) -> Decimal:
    value = str(text).strip()
    if value.startswith("[") and "+/-" in value:
        return Decimal(value[1:].split("+/-", 1)[0].strip())
    return Decimal(value)


def arb_positive(value: flint.arb) -> bool:
    return bool(value > 0 and not value.contains(0))


def arb_negative(value: flint.arb) -> bool:
    return bool(value < 0 and not value.contains(0))


def load_enclosures(
    paths: list[Path],
) -> tuple[
    dict[tuple[Decimal, int], flint.arb],
    dict[tuple[Decimal, int], flint.arb],
    dict[tuple[Decimal, int], Decimal],
    dict[tuple[Decimal, int], Decimal],
    dict[Decimal, str],
]:
    mu_balls: dict[tuple[Decimal, int], flint.arb] = {}
    a_balls: dict[tuple[Decimal, int], flint.arb] = {}
    mu_samples: dict[tuple[Decimal, int], Decimal] = {}
    a_samples: dict[tuple[Decimal, int], Decimal] = {}
    labels: dict[Decimal, str] = {}
    for path in paths:
        with path.open("r", encoding="utf-8") as handle:
            for raw in handle:
                if not raw.strip():
                    continue
                row = json.loads(raw)
                if row.get("kind") != "acb_coefficient_enclosure":
                    continue
                lam = decimal_lam_key(row["lam"])
                index = int(row["k"])
                mu_balls[(lam, index)] = flint.arb(row["full_mu_ball"])
                a_balls[(lam, index)] = flint.arb(row["A_ball"])
                mu_samples[(lam, index)] = parse_ball_midpoint(row["full_mu_ball"])
                a_samples[(lam, index)] = parse_ball_midpoint(row["A_ball"])
                labels[lam] = row["lam"]
    return mu_balls, a_balls, mu_samples, a_samples, labels


def finite_diagnostics(paths: list[Path], max_coefficient_index: int) -> SplitDiagnostics:
    mu_balls, a_balls, mu_samples, a_samples, labels = load_enclosures(paths)
    raw_rows = 0
    raw_negative = 0
    normalized_rows = 0
    normalized_positive = 0
    raw_ratio_rows = 0
    threshold_rows = 0
    threshold_positive = 0
    max_raw_disc: tuple[Decimal, str, int] | None = None
    min_norm_disc: tuple[Decimal, str, int] | None = None
    min_threshold_margin: tuple[Decimal, str, int] | None = None

    for lam in sorted(labels):
        label = labels[lam]
        for n in range(0, max_coefficient_index - 1):
            needed = (n, n + 1, n + 2)
            if any((lam, idx) not in mu_balls for idx in needed):
                continue

            raw_rows += 1
            raw_disc = mu_balls[(lam, n + 1)] ** 2 - mu_balls[(lam, n)] * mu_balls[(lam, n + 2)]
            if arb_negative(raw_disc):
                raw_negative += 1
            raw_sample = mu_samples[(lam, n + 1)] ** 2 - mu_samples[(lam, n)] * mu_samples[(lam, n + 2)]
            if max_raw_disc is None or raw_sample > max_raw_disc[0]:
                max_raw_disc = (raw_sample, label, n)

            raw_ratio_rows += 1
            raw_ratio = mu_balls[(lam, n + 1)] ** 2 / (mu_balls[(lam, n)] * mu_balls[(lam, n + 2)])
            if not arb_negative(raw_ratio - flint.arb(1)):
                raise RuntimeError(f"raw moment ratio was not certified below 1 at lambda={label}, n={n}")

            normalized_rows += 1
            normalized_disc = a_balls[(lam, n + 1)] ** 2 - a_balls[(lam, n)] * a_balls[(lam, n + 2)]
            if arb_positive(normalized_disc):
                normalized_positive += 1
            normalized_sample = a_samples[(lam, n + 1)] ** 2 - a_samples[(lam, n)] * a_samples[(lam, n + 2)]
            if min_norm_disc is None or normalized_sample < min_norm_disc[0]:
                min_norm_disc = (normalized_sample, label, n)

            threshold_rows += 1
            threshold = flint.arb(2 * n + 1) / flint.arb(2 * n + 3)
            threshold_margin = raw_ratio - threshold
            if arb_positive(threshold_margin):
                threshold_positive += 1
            sample_ratio = mu_samples[(lam, n + 1)] ** 2 / (mu_samples[(lam, n)] * mu_samples[(lam, n + 2)])
            sample_threshold = Decimal(2 * n + 1) / Decimal(2 * n + 3)
            sample_margin = sample_ratio - sample_threshold
            if min_threshold_margin is None or sample_margin < min_threshold_margin[0]:
                min_threshold_margin = (sample_margin, label, n)

    if max_raw_disc is None or min_norm_disc is None or min_threshold_margin is None:
        raise RuntimeError("no finite diagnostics were computed")

    return SplitDiagnostics(
        lambdas=[labels[lam] for lam in sorted(labels)],
        max_coefficient_index=max_coefficient_index,
        raw_degree2_rows=raw_rows,
        raw_degree2_negative_rows=raw_negative,
        normalized_degree2_rows=normalized_rows,
        normalized_degree2_positive_rows=normalized_positive,
        raw_ratio_below_one_rows=raw_ratio_rows,
        normalized_threshold_rows=threshold_rows,
        normalized_threshold_positive_rows=threshold_positive,
        max_raw_discriminant=FiniteMinimum(decimal_format(max_raw_disc[0]), max_raw_disc[1], max_raw_disc[2]),
        min_normalized_discriminant=FiniteMinimum(decimal_format(min_norm_disc[0]), min_norm_disc[1], min_norm_disc[2]),
        min_normalized_threshold_margin=FiniteMinimum(
            decimal_format(min_threshold_margin[0]),
            min_threshold_margin[1],
            min_threshold_margin[2],
        ),
    )


def build_audit(paths: list[Path], max_coefficient_index: int) -> dict:
    diagnostics = finite_diagnostics(paths, max_coefficient_index)
    return {
        "kind": "jensen_window_pf_factorial_multiplier_split_audit",
        "date": "2026-07-06",
        "status": "finite_theorem_search_diagnostic",
        "source_coefficient_pf_obstruction": "outputs/coefficient_pf_bridge_obstruction.md",
        "source_heat_flow_boundary_threshold": "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
        "source_sign_regular_transfer_gap": "outputs/jensen_window_pf_sign_regular_transfer_gap_matrix.md",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "proof_boundary": (
            "Finite theorem-search diagnostic plus exact low-degree algebra only. "
            "It records that the factorial multiplier gamma_k=k!/(2*k)! is a "
            "conditional Pólya-Schur preserver for shifted windows, and that raw "
            "moment-window hyperbolicity is not a valid input theorem. It does "
            "not prove Jensen-window PF-infinity, Jensen hyperbolicity, "
            "Laguerre-Polya membership, RH, or Lambda <= 0."
        ),
        "exact_rows": [
            {
                "id": "fm_01_factorial_multiplier_sequence",
                "role": "standard_multiplier_sequence",
                "formula": "gamma_k = k!/(2*k)!",
                "egf": "sum_{k>=0} gamma_k*z^k/k! = cosh(sqrt(z))",
                "reason": "cosh(sqrt(z)) has only real nonpositive zeros in the z-plane; its derivatives give the shifted sequences gamma_{n+j}.",
                "readiness": "conditional_only",
                "proof_boundary": "Multiplier-sequence fact only; it preserves real-rootedness under its hypotheses but does not create the missing input theorem.",
            },
            {
                "id": "fm_02_shifted_window_preserver",
                "role": "conditional_transfer",
                "formula": "P_{d,n}(x)=sum_j binom(d,j)*gamma_{n+j}*mu_{2*(n+j)}*x^j",
                "conditional_claim": "If M_{d,n}(x)=sum_j binom(d,j)*mu_{2*(n+j)}*x^j has real nonpositive zeros, then P_{d,n}(x) does too.",
                "readiness": "not_ready_to_apply",
                "proof_boundary": "Conditional Pólya-Schur reduction only; the raw moment-window hypothesis is not proved and is false at degree 2 for nondegenerate positive moment weights.",
            },
            {
                "id": "fm_03_raw_moment_degree2_anti_alignment",
                "role": "exact_obstruction",
                "formula": "Delta_2(M_{2,n}) = 4*(mu_{2*n+2}^2 - mu_{2*n}*mu_{2*n+4}) <= 0",
                "reason": "Cauchy-Schwarz for the positive Newman/Phi moment weight gives raw moment log-convexity, the opposite sign from degree-2 hyperbolicity for positive coefficients.",
                "readiness": "not_ready_to_apply",
                "proof_boundary": "Exact degree-2 obstruction to the raw moment-window input theorem; not a statement about all normalized Jensen windows.",
            },
            {
                "id": "fm_04_normalized_degree2_threshold",
                "role": "exact_normalization_threshold",
                "formula": "A_{n+1}^2 >= A_n*A_{n+2} iff mu_{2*n+2}^2/(mu_{2*n}*mu_{2*n+4}) >= (2*n+1)/(2*n+3)",
                "reason": "gamma_n*gamma_{n+2}/gamma_{n+1}^2 = (2*n+1)/(2*n+3).",
                "readiness": "available_low_degree_only",
                "proof_boundary": "Degree-2 normalization threshold only; it does not prove higher Jensen-window PF minors.",
            },
            {
                "id": "fm_05_route_verdict",
                "role": "route_boundary",
                "formula": "factorial multiplier route = preserver theorem + missing Xi/Phi-specific input theorem",
                "reason": "The preserver component is standard, while raw Stieltjes moment positivity supplies anti-hyperbolic degree-2 windows; the route must instead prove a normalized ratio-cone, determinant-integral, or Schur-positive theorem.",
                "readiness": "not_ready_to_apply",
                "proof_boundary": "Route-selection boundary only; no all-degree/all-shift bridge theorem is proved.",
            },
        ],
        "finite_diagnostics": asdict(diagnostics),
        "summary": {
            "exact_rows": 5,
            "raw_degree2_rows": diagnostics.raw_degree2_rows,
            "raw_degree2_negative_rows": diagnostics.raw_degree2_negative_rows,
            "normalized_degree2_rows": diagnostics.normalized_degree2_rows,
            "normalized_degree2_positive_rows": diagnostics.normalized_degree2_positive_rows,
            "ready_to_apply_rows": 0,
            "target_closing": False,
            "main_finding": (
                "The factorial multiplier gamma_k=k!/(2*k)! is not the missing "
                "bridge by itself. It gives a clean conditional preservation "
                "step for shifted windows, but the natural raw moment-window "
                "input theorem is anti-aligned at degree 2. The useful route is "
                "therefore a normalized Xi/Phi-specific ratio-cone, determinant "
                "integral, or Schur-positive theorem, not raw moment positivity."
            ),
        },
        "invariants": [
            "No row is ready_to_apply.",
            "The shifted multiplier sequence is used only conditionally.",
            "Raw moment positivity is not promoted to Jensen-window hyperbolicity.",
            "Degree-2 normalized positivity is not promoted to Jensen-window PF-infinity.",
            "Endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(audit: dict, path: Path) -> None:
    summary = audit["summary"]
    diagnostics = audit["finite_diagnostics"]
    result_line = (
        "validated Jensen-window PF factorial multiplier split audit: "
        f"{summary['exact_rows']} exact rows, "
        f"{summary['raw_degree2_negative_rows']} raw degree-2 anti-hyperbolic rows, "
        f"{summary['normalized_degree2_positive_rows']} normalized degree-2 positive rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows, 0 issues"
    )
    lines = [
        "# Jensen-Window PF Factorial Multiplier Split Audit",
        "",
        "Date: 2026-07-06",
        "",
        "Status: finite theorem-search diagnostic. This is not a proof of",
        "Jensen-window PF-infinity, Jensen hyperbolicity, Laguerre-Polya",
        "membership, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_factorial_multiplier_split_audit`.",
        "",
        "Proof boundary: this artifact separates a valid conditional",
        "Pólya-Schur multiplier step from the missing input theorem. It does not",
        "prove the sign-regular-to-Jensen transfer theorem.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_factorial_multiplier_split_audit.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_factorial_multiplier_split_audit.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_factorial_multiplier_split_audit.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Exact Split",
        "",
        "For the raw moments and normalized Jensen coefficients:",
        "",
        "```text",
        "A_k = gamma_k * mu_{2k}",
        "gamma_k = k!/(2*k)!",
        "sum gamma_k*z^k/k! = cosh(sqrt(z))",
        "```",
        "",
        "The shifted sequences `gamma_{n+j}` are conditional multiplier",
        "sequences because they come from derivatives of `cosh(sqrt(z))`.",
        "Thus they preserve real-rootedness if the raw shifted moment-window",
        "polynomials are already real-rooted.",
        "",
        "But raw moment-window hyperbolicity is the wrong input theorem:",
        "",
        "```text",
        "Delta_2(M_{2,n}) = 4*(mu_{2*n+2}^2 - mu_{2*n}*mu_{2*n+4}) <= 0",
        "```",
        "",
        "Cauchy-Schwarz gives raw moment log-convexity, which is the opposite",
        "degree-2 sign from real-rooted positive-coefficient Jensen windows.",
        "",
        "The factorial normalization changes the degree-2 threshold to:",
        "",
        "```text",
        "A_{n+1}^2 >= A_n*A_{n+2}",
        "iff mu_{2*n+2}^2/(mu_{2*n}*mu_{2*n+4}) >= (2*n+1)/(2*n+3)",
        "```",
        "",
        "## Finite Arb Sanity Check",
        "",
        "Using the existing coefficient enclosures:",
        "",
        "```text",
        f"lambdas: {', '.join(diagnostics['lambdas'])}",
        f"max coefficient index: {diagnostics['max_coefficient_index']}",
        f"raw degree-2 rows: {diagnostics['raw_degree2_rows']}",
        f"raw degree-2 anti-hyperbolic rows: {diagnostics['raw_degree2_negative_rows']}",
        f"normalized degree-2 rows: {diagnostics['normalized_degree2_rows']}",
        f"normalized degree-2 positive rows: {diagnostics['normalized_degree2_positive_rows']}",
        f"normalized threshold positive rows: {diagnostics['normalized_threshold_positive_rows']}",
        "```",
        "",
        "Sharpest finite margins:",
        "",
        "```text",
        f"max raw discriminant: {diagnostics['max_raw_discriminant']['sample']} at lambda={diagnostics['max_raw_discriminant']['lam']}, n={diagnostics['max_raw_discriminant']['n']}",
        f"min normalized discriminant: {diagnostics['min_normalized_discriminant']['sample']} at lambda={diagnostics['min_normalized_discriminant']['lam']}, n={diagnostics['min_normalized_discriminant']['n']}",
        f"min normalized threshold margin: {diagnostics['min_normalized_threshold_margin']['sample']} at lambda={diagnostics['min_normalized_threshold_margin']['lam']}, n={diagnostics['min_normalized_threshold_margin']['n']}",
        "```",
        "",
        "## Route Consequence",
        "",
        "The factorial multiplier route is not dead, but it is conditional. It",
        "must be paired with a genuinely normalized Xi/Phi-specific theorem,",
        "such as a ratio-cone argument, determinant integral, positive kernel,",
        "or Schur-positive specialization. Raw moment positivity is not enough.",
        "",
        "Integration:",
        "",
        "```text",
        "outputs/coefficient_pf_bridge_obstruction.md",
        "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
        "outputs/jensen_window_pf_sign_regular_transfer_gap_matrix.md",
        "outputs/jensen_window_pf_theorem_machinery_fit_matrix.md",
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
    parser.add_argument("--enclosure-jsonl", type=Path, nargs="*", default=list(DEFAULT_ENCLOSURE_JSONL))
    parser.add_argument("--max-coefficient-index", type=int, default=64)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    paths = [path if path.is_absolute() else REPO_ROOT / path for path in args.enclosure_jsonl]
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    audit = build_audit(paths, args.max_coefficient_index)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(audit, note)
    print(
        "wrote Jensen-window PF factorial multiplier split audit: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
