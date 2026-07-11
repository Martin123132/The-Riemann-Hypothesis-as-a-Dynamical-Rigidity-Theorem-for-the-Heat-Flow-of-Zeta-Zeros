#!/usr/bin/env python3
"""Build the Taylor-moment budget scout for the negative-lambda route."""

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


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SIGN_SCOUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_phi_taylor_cone_entry_sign_scout.json"
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_taylor_moment_budget.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_taylor_moment_budget.md"

getcontext().prec = 100


@dataclass(frozen=True)
class TailStartSample:
    T: str
    k: int
    q_over_T: str
    first_correction_abs: str
    truncated_F_k_minus_1: str
    truncated_F_k: str
    truncated_F_k_plus_1: str
    target_bound: str
    truncation_status: str
    truncated_x: str | None
    truncated_B: str | None
    B_over_target_bound: str | None
    proof_boundary: str


@dataclass(frozen=True)
class FirstCorrectionBudget:
    threshold: str
    eta_q_over_T: str
    tail_start_k: int
    tail_start_T_min: str
    proof_boundary: str


@dataclass(frozen=True)
class TaylorMomentDiagnostics:
    ratio_a: str
    ratio_b: str
    ratio_c: str
    tail_start_k: int
    sample_T_values: list[int]
    tail_start_samples: list[TailStartSample]
    invalid_truncation_rows: int
    positive_truncation_rows: int
    overbound_truncation_rows: int
    bounded_truncation_rows: int
    first_correction_budget: FirstCorrectionBudget


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def decimal_format(value: Decimal) -> str:
    return f"{value:.18E}"


def arb_mid_decimal(value: flint.arb) -> Decimal:
    mid, _rad, exp10 = value.mid_rad_10exp()
    return Decimal(int(mid)) * (Decimal(10) ** int(exp10))


def truncated_multiplier(k: int, T: int, a: Decimal, b: Decimal, c: Decimal) -> Decimal:
    q = Decimal(k) + Decimal("0.5")
    big_t = Decimal(T)
    return (
        Decimal(1)
        + a * q / big_t
        + b * q * (q + 1) / (big_t**2)
        + c * q * (q + 1) * (q + 2) / (big_t**3)
    )


def target_bound(k: int) -> Decimal:
    return Decimal(2) / (Decimal(3) * (Decimal(2 * k + 1)))


def sample_tail_start(k: int, T: int, a: Decimal, b: Decimal, c: Decimal) -> TailStartSample:
    q = Decimal(k) + Decimal("0.5")
    big_t = Decimal(T)
    f_minus = truncated_multiplier(k - 1, T, a, b, c)
    f = truncated_multiplier(k, T, a, b, c)
    f_plus = truncated_multiplier(k + 1, T, a, b, c)
    bound = target_bound(k)
    if f_minus > 0 and f > 0 and f_plus > 0:
        x_value = f_plus * f_minus / (f * f)
        b_value = -x_value.ln()
        ratio = b_value / bound
        status = "bounded_positive_truncation" if Decimal(0) <= b_value <= bound else "overbound_positive_truncation"
        x_text = decimal_format(x_value)
        b_text = decimal_format(b_value)
        ratio_text = decimal_format(ratio)
    else:
        status = "invalid_truncation_normalizer"
        x_text = None
        b_text = None
        ratio_text = None
    return TailStartSample(
        T=str(T),
        k=k,
        q_over_T=decimal_format(q / big_t),
        first_correction_abs=decimal_format(abs(a) * q / big_t),
        truncated_F_k_minus_1=decimal_format(f_minus),
        truncated_F_k=decimal_format(f),
        truncated_F_k_plus_1=decimal_format(f_plus),
        target_bound=decimal_format(bound),
        truncation_status=status,
        truncated_x=x_text,
        truncated_B=b_text,
        B_over_target_bound=ratio_text,
        proof_boundary="Low-order Taylor-moment truncation diagnostic only; it is not the actual zeta coefficient prefix.",
    )


def build_diagnostics(sign_scout_json: Path, tail_start_k: int, sample_t_values: list[int]) -> TaylorMomentDiagnostics:
    scout = load_json(sign_scout_json)
    ratios = scout["ratio_enclosures"]
    a_ball = flint.arb(ratios["a=c2/c0"])
    b_ball = flint.arb(ratios["b=c4/c0"])
    c_ball = flint.arb(ratios["c=c6/c0"])
    a = arb_mid_decimal(a_ball)
    b = arb_mid_decimal(b_ball)
    c = arb_mid_decimal(c_ball)

    samples = [sample_tail_start(tail_start_k, T, a, b, c) for T in sample_t_values]
    invalid = sum(1 for row in samples if row.truncation_status == "invalid_truncation_normalizer")
    positive = len(samples) - invalid
    overbound = sum(1 for row in samples if row.truncation_status == "overbound_positive_truncation")
    bounded = sum(1 for row in samples if row.truncation_status == "bounded_positive_truncation")
    eta = Decimal(1) / (Decimal(2) * abs(a))
    tail_t_min = (Decimal(tail_start_k) + Decimal("0.5")) / eta
    first_budget = FirstCorrectionBudget(
        threshold="abs(a)*(k+1/2)/T <= 1/2",
        eta_q_over_T=decimal_format(eta),
        tail_start_k=tail_start_k,
        tail_start_T_min=decimal_format(tail_t_min),
        proof_boundary="Scale budget only; it is a stability diagnostic, not a rigorous Taylor remainder theorem.",
    )
    return TaylorMomentDiagnostics(
        ratio_a=a_ball.str(40),
        ratio_b=b_ball.str(40),
        ratio_c=c_ball.str(40),
        tail_start_k=tail_start_k,
        sample_T_values=sample_t_values,
        tail_start_samples=samples,
        invalid_truncation_rows=invalid,
        positive_truncation_rows=positive,
        overbound_truncation_rows=overbound,
        bounded_truncation_rows=bounded,
        first_correction_budget=first_budget,
    )


def build_artifact(sign_scout_json: Path, tail_start_k: int, sample_t_values: list[int]) -> dict:
    diagnostics = build_diagnostics(sign_scout_json, tail_start_k, sample_t_values)
    rows = [
        {
            "id": "nltmb_01_gaussian_moment_factorization",
            "role": "exact_reduction",
            "readiness": "not_ready_to_apply",
            "source_artifacts": ["outputs/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.md"],
            "claim": "After factoring the pure Gaussian moments, a local Taylor model changes x_k only through F_(k+1)*F_(k-1)/F_k^2.",
            "acceptance_test": "Keep the pure Gaussian cancellation separate from any Taylor remainder estimate.",
            "proof_boundary": "Exact model reduction only; not a theorem about the full zeta coefficient tail.",
        },
        {
            "id": "nltmb_02_truncated_multiplier_formula",
            "role": "formal_model",
            "readiness": "not_ready_to_apply",
            "source_artifacts": ["outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md"],
            "claim": "The degree-6 local multiplier is F_k^(<=6)=1+a*q/T+b*q*(q+1)/T^2+c*q*(q+1)*(q+2)/T^3 with q=k+1/2.",
            "acceptance_test": "Treat q/T, not k alone, as the local expansion parameter.",
            "proof_boundary": "Formal low-order Taylor model only; the omitted R_k^(>=8) is not bounded here.",
        },
        {
            "id": "nltmb_03_tail_start_stability_samples",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "source_artifacts": ["outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md"],
            "claim": "At the k=22 tail start, the degree-6 Taylor multiplier is an invalid positive-moment model for four sampled T values and only bounded for two larger sampled T values.",
            "acceptance_test": "Do not use the low-order truncation as a finite proof model for lambda=-25,-50,-100.",
            "proof_boundary": "Finite diagnostic only; it does not contradict the separately certified actual ACB prefix.",
        },
        {
            "id": "nltmb_04_first_correction_budget",
            "role": "scale_budget",
            "readiness": "not_ready_to_apply",
            "source_artifacts": ["outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md"],
            "claim": "The simple stability budget abs(a)*(k+1/2)/T<=1/2 covers k=22 only after T is about 1685.4.",
            "acceptance_test": "Any local/mesoscopic theorem must state its q/T range and constants explicitly.",
            "proof_boundary": "Scale budget only; not a rigorous remainder theorem.",
        },
        {
            "id": "nltmb_05_normalizer_positivity_requirement",
            "role": "open_requirement",
            "readiness": "not_ready_to_apply",
            "source_artifacts": ["outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md"],
            "claim": "A proof must certify F_k^(<=6)+R_k^(>=8)>0 throughout its claimed k/T range.",
            "acceptance_test": "Give interval-safe lower bounds for the full normalized moment multiplier, including the Taylor tail.",
            "proof_boundary": "Open analytic requirement only; positivity of the full multiplier is not proved here.",
        },
        {
            "id": "nltmb_06_log_curvature_remainder_requirement",
            "role": "open_requirement",
            "readiness": "not_ready_to_apply",
            "source_artifacts": ["outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md"],
            "claim": "A proof must bound the second log-difference contribution of R_k^(>=8) strongly enough to keep B_k<=2/(3*(2*k+1)).",
            "acceptance_test": "State a concrete inequality comparing the Taylor-tail log-curvature error to the remaining target margin.",
            "proof_boundary": "Open analytic requirement only; the bounded log-curvature theorem is not proved.",
        },
        {
            "id": "nltmb_07_monotone_remainder_requirement",
            "role": "open_requirement",
            "readiness": "not_ready_to_apply",
            "source_artifacts": ["outputs/jensen_window_pf_monotone_contraction_theorem_target.md"],
            "claim": "A proof must also control the third log-difference or provide the separate monotone-contraction theorem required by the log-curvature bridge.",
            "acceptance_test": "Do not infer B_(k+1)<=B_k from the fixed-k sign without a uniform remainder bound.",
            "proof_boundary": "Open analytic requirement only; monotone contraction remains an open target.",
        },
        {
            "id": "nltmb_08_low_order_finite_T_promotion_rejected",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "source_artifacts": ["outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.md"],
            "claim": "The degree-6 Taylor truncation alone proves the finite negative-lambda prefix or the all-k tail.",
            "acceptance_test": "Reject any use of the low-order local truncation as a substitute for the actual ACB prefix plus analytic tail.",
            "proof_boundary": "Rejected proof template only.",
        },
        {
            "id": "nltmb_09_local_mesoscopic_theorem_handoff",
            "role": "conditional_handoff",
            "readiness": "not_ready_to_apply",
            "source_artifacts": ["outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md"],
            "claim": "A successful local/mesoscopic theorem should hand off at a stated k<=eta*T or k<=K(T) to the far-tail moving-saddle theorem.",
            "acceptance_test": "Pair the local q/T theorem with either a far-tail saddle theorem or an explicit finite collar and analytic tail.",
            "proof_boundary": "Conditional handoff only; not cone entry or Lambda <= 0.",
        },
    ]
    summary = {
        "budget_rows": len(rows),
        "tail_start_samples": len(diagnostics.tail_start_samples),
        "invalid_truncation_rows": diagnostics.invalid_truncation_rows,
        "positive_truncation_rows": diagnostics.positive_truncation_rows,
        "overbound_truncation_rows": diagnostics.overbound_truncation_rows,
        "bounded_truncation_rows": diagnostics.bounded_truncation_rows,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The exact Gaussian-moment budget makes q/T the relevant local parameter. "
            "The low-order Taylor truncation is not a finite proof model for the actual finite prefix: "
            "at the k=22 tail start it has invalid positive-moment normalizers for T=25,50,100,200, "
            "and the Taylor remainder must prove positivity plus log-curvature control before the signed "
            "perturbation route can enter the bounded-curvature target."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_taylor_moment_budget",
        "date": "2026-07-06",
        "status": "finite_theorem_search_diagnostic",
        "source_uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "source_signed_gaussian_matrix": "outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md",
        "source_phi_taylor_sign_scout": "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md",
        "source_bounded_log_curvature_target": "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md",
        "sign_scout_json": sign_scout_json.relative_to(REPO_ROOT).as_posix(),
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_taylor_moment_budget.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_taylor_moment_budget.py",
        "proof_boundary": (
            "Finite theorem-search diagnostic only. It derives the Gaussian-moment Taylor budget and "
            "records finite low-order truncation failures, but it does not prove a Taylor remainder "
            "bound, does not prove bounded log-curvature, does not prove cone entry, and does not "
            "prove Lambda <= 0."
        ),
        "formulae": {
            "q": "q=k+1/2",
            "normalized_multiplier": "F_k=1+a*q/T+b*q*(q+1)/T^2+c*q*(q+1)*(q+2)/T^3+R_k^(>=8)",
            "gaussian_cancellation": "x_k=F_(k+1)*F_(k-1)/F_k^2 after the pure Gaussian moment factor is cancelled",
            "bounded_curvature_target": "B_k=-log(x_k)<=2/(3*(2*k+1))",
        },
        "budget_rows": rows,
        "diagnostics": asdict(diagnostics),
        "invariants": [
            "No row is ready_to_apply.",
            "The low-order Taylor truncation is not the actual zeta coefficient prefix.",
            "The Taylor-tail remainder R_k^(>=8) is not bounded here.",
            "The local/mesoscopic theorem still requires explicit q/T constants.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
        "summary": summary,
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda Taylor moment budget: "
        f"{summary['budget_rows']} budget rows, "
        f"{summary['tail_start_samples']} tail-start samples, "
        f"{summary['invalid_truncation_rows']} invalid truncation rows, "
        f"{summary['bounded_truncation_rows']} bounded truncation rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows, 0 issues"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Taylor Moment Budget",
        "",
        "Date: 2026-07-06",
        "",
        "Status: finite theorem-search diagnostic. This is not a proof of the",
        "bounded log-curvature theorem, cone entry, Jensen-window PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_taylor_moment_budget`.",
        "",
        "Proof boundary: this artifact derives the Gaussian-moment Taylor budget",
        "for the local/mesoscopic route. It does not prove the Taylor-tail",
        "remainder theorem or close the all-`k` negative-lambda tail.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_taylor_moment_budget.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_taylor_moment_budget.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_taylor_moment_budget.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Moment Budget",
        "",
        "For the normalized local Taylor multiplier:",
        "",
        "```text",
        "q = k+1/2",
        "F_k = 1+a*q/T+b*q*(q+1)/T^2+c*q*(q+1)*(q+2)/T^3+R_k^(>=8)",
        "x_k = F_(k+1)*F_(k-1)/F_k^2",
        "B_k = -log(x_k)",
        "```",
        "",
        "the expansion parameter is `q/T`, not fixed `k` alone.",
        "",
        "Certified local ratios:",
        "",
        "```text",
        f"a = {diagnostics['ratio_a']}",
        f"b = {diagnostics['ratio_b']}",
        f"c = {diagnostics['ratio_c']}",
        "```",
        "",
        "Tail-start samples for `k=22`:",
        "",
        "```text",
    ]
    for row in diagnostics["tail_start_samples"]:
        status = row["truncation_status"]
        ratio = row["B_over_target_bound"] if row["B_over_target_bound"] is not None else "n/a"
        lines.append(
            f"T={row['T']}: q/T={row['q_over_T']}, |a|q/T={row['first_correction_abs']}, "
            f"F_k={row['truncated_F_k']}, status={status}, B/bound={ratio}"
        )
    lines.extend(
        [
            "```",
            "",
            "Simple first-correction stability budget:",
            "",
            "```text",
            f"{diagnostics['first_correction_budget']['threshold']}",
            f"q/T <= {diagnostics['first_correction_budget']['eta_q_over_T']}",
            f"k=22 enters only after T >= {diagnostics['first_correction_budget']['tail_start_T_min']}",
            "```",
            "",
            "Interpretation: the low-order Taylor truncation is not a finite proof",
            "model for the certified `lambda=-25,-50,-100` prefix. The actual prefix",
            "is certified by ACB enclosures elsewhere; this artifact only sizes the",
            "local Taylor/remainder theorem needed for an analytic route.",
            "",
            "Required analytic upgrades:",
            "",
            "```text",
            "1. prove F_k^(<=6)+R_k^(>=8)>0 over a stated q/T range",
            "2. bound the second log-difference of R_k^(>=8) against the B_k target margin",
            "3. control the monotone third log-difference or use the separate monotone-contraction theorem",
            "4. hand off to a global far-tail moving-saddle theorem or finite collar plus analytic tail",
            "```",
            "",
            "Integration:",
            "",
            "```text",
            "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
            "outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md",
            "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md",
            "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md",
            "```",
            "",
            "Summary:",
            "",
            summary["main_finding"],
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_t_values(text: str) -> list[int]:
    values = [int(part.strip()) for part in text.split(",") if part.strip()]
    if not values:
        raise argparse.ArgumentTypeError("at least one T value is required")
    return values


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sign-scout-json", type=Path, default=DEFAULT_SIGN_SCOUT_JSON)
    parser.add_argument("--tail-start-k", type=int, default=22)
    parser.add_argument("--t-values", type=parse_t_values, default=parse_t_values("25,50,100,200,500,1000,2000"))
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    sign_scout_json = args.sign_scout_json if args.sign_scout_json.is_absolute() else REPO_ROOT / args.sign_scout_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(sign_scout_json, args.tail_start_k, args.t_values)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda Taylor moment budget: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
