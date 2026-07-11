#!/usr/bin/env python3
"""Build the uniform remainder target for the negative-lambda perturbation route."""

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
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_uniform_remainder_target.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md"

getcontext().prec = 100


@dataclass(frozen=True)
class LeadingScaleRow:
    T: str
    leading_only_k_max_floor: int
    leading_only_k_max_sample: str
    target_tail_start_k: int
    proof_boundary: str


@dataclass(frozen=True)
class UniformRemainderDiagnostics:
    deficit_coefficient: str
    monotone_coefficient: str
    activation_depth: int
    activation_t_lower: str
    target_tail_start_k: int
    leading_scale_rows: list[LeadingScaleRow]
    rows_below_tail_start: int


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def decimal_format(value: Decimal) -> str:
    return f"{value:.18E}"


def arb_mid_decimal(value: flint.arb) -> Decimal:
    mid, _rad, exp10 = value.mid_rad_10exp()
    return Decimal(int(mid)) * (Decimal(10) ** int(exp10))


def build_diagnostics(sign_scout_json: Path, activation_depth: int, t_values: list[int]) -> UniformRemainderDiagnostics:
    scout = load_json(sign_scout_json)
    ratios = scout["ratio_enclosures"]
    a = flint.arb(ratios["a=c2/c0"])
    b = flint.arb(ratios["b=c4/c0"])
    c = flint.arb(ratios["c=c6/c0"])
    deficit = a * a - 2 * b
    monotone = 2 * (a**3 - 3 * a * b + 3 * c)
    if not bool(deficit > 0 and not deficit.contains(0)):
        raise RuntimeError(f"deficit coefficient is not positive: {deficit}")
    if not bool(monotone > 0 and not monotone.contains(0)):
        raise RuntimeError(f"monotone coefficient is not positive: {monotone}")

    activation = (deficit * flint.arb(3 * (2 * activation_depth + 1)) / flint.arb(2)).sqrt()
    d2 = arb_mid_decimal(deficit)
    rows: list[LeadingScaleRow] = []
    for t_value in t_values:
        T = Decimal(t_value)
        # Leading-only inequality: D2/T^2 <= 2/(3*(2*k+1)).
        k_max = (Decimal(2) * T * T / (Decimal(3) * d2) - Decimal(1)) / Decimal(2)
        rows.append(
            LeadingScaleRow(
                T=str(t_value),
                leading_only_k_max_floor=max(-1, int(k_max)),
                leading_only_k_max_sample=decimal_format(k_max),
                target_tail_start_k=activation_depth + 1,
                proof_boundary="Leading fixed-k estimate only; it ignores all remainders and cannot prove the tail.",
            )
        )
    rows_below = sum(1 for row in rows if row.leading_only_k_max_floor < row.target_tail_start_k)
    return UniformRemainderDiagnostics(
        deficit_coefficient=deficit.str(40),
        monotone_coefficient=monotone.str(40),
        activation_depth=activation_depth,
        activation_t_lower=activation.str(40),
        target_tail_start_k=activation_depth + 1,
        leading_scale_rows=rows,
        rows_below_tail_start=rows_below,
    )


def build_target(sign_scout_json: Path, activation_depth: int, t_values: list[int]) -> dict:
    diagnostics = build_diagnostics(sign_scout_json, activation_depth, t_values)
    summary = {
        "target_rows": 8,
        "ready_to_apply_rows": 0,
        "open_requirement_rows": 2,
        "live_routes": 2,
        "rejected_routes": 1,
        "leading_scale_rows": len(diagnostics.leading_scale_rows),
        "leading_scale_rows_below_tail_start": diagnostics.rows_below_tail_start,
        "target_closing": False,
        "main_finding": (
            "The signed Gaussian perturbation route now needs a two-scale uniform remainder theorem: "
            "the fixed-k Taylor signs give the correct local direction, but the bounded-curvature target "
            "shrinks like 1/k, so a leading fixed-k D2/T^2 estimate cannot be promoted to the all-k tail. "
            "A proof must either control the moving saddle globally in k/T or combine a finite collar with "
            "a separate far-tail theorem."
        ),
    }
    rows = [
        {
            "id": "nlurt_01_fixed_k_scope",
            "role": "formal_scope",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md",
                "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md",
            ],
            "claim_if_proved": "The signed perturbation expansion gives the correct fixed-k signs as T tends to infinity.",
            "gap": "Fixed-k signs do not control the infinite tail where k is unbounded.",
            "acceptance_test": "State the scale in k/T where each remainder bound is valid.",
            "proof_boundary": "Formal fixed-k scope only; not a uniform theorem.",
        },
        {
            "id": "nlurt_02_shrinking_target_obstruction",
            "role": "exact_obstruction",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md",
            ],
            "claim_if_proved": "The bounded-curvature target is O(1/k), while the leading fixed-k deficit D2/T^2 is independent of k.",
            "gap": "For any fixed T, a leading-only fixed-k estimate eventually exceeds the all-k target bound.",
            "acceptance_test": "Do not promote fixed-k asymptotics to an all-k tail without changing scale.",
            "proof_boundary": "Exact scale obstruction only; not a statement about the full zeta coefficient tail.",
        },
        {
            "id": "nlurt_03_leading_scale_diagnostic",
            "role": "finite_scale_diagnostic",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md",
            ],
            "claim_if_proved": "The leading-only cutoff estimates quantify why the fixed-k model is not the finite-prefix proof at T=25,50,100.",
            "gap": "These estimates ignore remainders and do not decide the actual finite prefix.",
            "acceptance_test": "Use them only to size the asymptotic route.",
            "proof_boundary": "Scale diagnostic only.",
        },
        {
            "id": "nlurt_04_local_uniform_remainder_requirement",
            "role": "open_requirement",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md",
            ],
            "claim_if_proved": "For a controlled local or mesoscopic range, bound the signed Taylor remainder by a fixed fraction of the positive deficit and monotone correction.",
            "gap": "No explicit k/T-dependent remainder bound is proved.",
            "acceptance_test": "Give constants and a range such as k<=alpha*T or k<=K(T), with interval-safe inequalities.",
            "proof_boundary": "Open analytic requirement only; the required k/T-dependent constants are not proved.",
        },
        {
            "id": "nlurt_05_far_tail_saddle_requirement",
            "role": "open_requirement",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md",
            ],
            "claim_if_proved": "For k beyond the local range, prove the bounded log-curvature estimate using a moving-saddle or global Phi-kernel argument.",
            "gap": "The saddle location moves away from u=0 as k/T grows, so the local Phi Taylor expansion is not enough.",
            "acceptance_test": "Produce a global estimate valid to infinity, or a finite cutoff plus a separate tail theorem.",
            "proof_boundary": "Open analytic requirement only.",
        },
        {
            "id": "nlurt_06_collared_finite_tail_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_finite_collar_contract.md",
                "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
            ],
            "claim_if_proved": "Extend the finite collar to a larger K, then prove a tail theorem starting at K+1.",
            "gap": "Finite extension alone never proves the infinite tail.",
            "acceptance_test": "Pair any new A_k enclosures with an analytic tail theorem and explicit starting index.",
            "proof_boundary": "Live route only; not a closed proof.",
        },
        {
            "id": "nlurt_07_fixed_k_promotion_rejected",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md",
            ],
            "claim_if_proved": "The fixed-k Taylor signs alone prove the all-k negative-lambda tail.",
            "gap": "This promotion ignores the shrinking target and moving saddle.",
            "acceptance_test": "Reject any manuscript step that uses only fixed-k signs for the all-k tail.",
            "proof_boundary": "Rejected proof template only.",
        },
        {
            "id": "nlurt_08_conditional_bounded_curvature_application",
            "role": "conditional_application",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md",
                "outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md",
            ],
            "claim_if_proved": "A two-scale uniform remainder theorem would supply the bounded-curvature half of the buffered defect-tail route.",
            "gap": "The monotone-contraction theorem and the uniform remainder theorem are both still open.",
            "acceptance_test": "After proof, re-run the bounded-curvature and defect-tail targets without changing forbidden endpoint assumptions.",
            "proof_boundary": "Conditional application only; not cone entry or Lambda <= 0.",
        },
    ]
    return {
        "kind": "jensen_window_pf_negative_lambda_uniform_remainder_target",
        "date": "2026-07-06",
        "status": "open_theorem_target",
        "target_id": "target_negative_lambda_uniform_remainder",
        "source_signed_gaussian_matrix": "outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md",
        "source_bounded_log_curvature_target": "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md",
        "source_defect_tail_target": "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
        "sign_scout_json": sign_scout_json.relative_to(REPO_ROOT).as_posix(),
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_uniform_remainder_target.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_uniform_remainder_target.py",
        "proof_boundary": (
            "Open theorem target only. It states the two-scale uniform remainder theorem needed by the "
            "signed Gaussian perturbation route, but it does not prove that theorem, does not prove the "
            "bounded log-curvature target, does not prove cone entry, and does not prove Lambda <= 0."
        ),
        "target_rows": rows,
        "finite_diagnostics": asdict(diagnostics),
        "invariants": [
            "No row is ready_to_apply.",
            "The target remains open_theorem_target.",
            "Fixed-k Taylor signs are not promoted to an all-k tail theorem.",
            "The moving-saddle or far-tail requirement remains open.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
        "summary": summary,
    }


def write_note(target: dict, path: Path) -> None:
    summary = target["summary"]
    diagnostics = target["finite_diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda uniform remainder target: "
        f"{summary['target_rows']} rows, 0 issues, {summary['ready_to_apply_rows']} ready-to-apply rows, "
        f"{summary['open_requirement_rows']} open requirements, {summary['leading_scale_rows']} leading-scale rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Uniform Remainder Target",
        "",
        "Date: 2026-07-06",
        "",
        "Status: open theorem target. This is not a proof of the bounded",
        "log-curvature theorem, cone entry, Jensen-window PF-infinity, RH, or",
        "`Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_uniform_remainder_target`.",
        "",
        "Proof boundary: this artifact states the two-scale uniform remainder",
        "theorem needed by the signed perturbation route. It does not prove that",
        "theorem or close the negative-lambda tail.",
        "",
        "Machine-readable target:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_uniform_remainder_target.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_uniform_remainder_target.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_uniform_remainder_target.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Scale Obstruction",
        "",
        "The fixed-`k` signed perturbation gives a leading deficit:",
        "",
        "```text",
        "B_k ~ D2/T^2",
        f"D2 = {diagnostics['deficit_coefficient']}",
        "```",
        "",
        "but the target bound shrinks with `k`:",
        "",
        "```text",
        "B_k <= 2/(3*(2*k+1))",
        "```",
        "",
        "Therefore a fixed-`k` expansion cannot be the all-tail theorem by itself.",
        "",
        "Leading-only scale diagnostics:",
        "",
        "```text",
    ]
    for row in diagnostics["leading_scale_rows"]:
        lines.append(
            f"T={row['T']}: k_max~{row['leading_only_k_max_sample']} "
            f"(floor {row['leading_only_k_max_floor']}), tail starts at k={row['target_tail_start_k']}"
        )
    lines.extend(
        [
            "```",
            "",
            "These are scale diagnostics only. They ignore all remainders and do not",
            "contradict the certified finite prefix.",
            "",
            "## Required Theorem",
            "",
            "A proof must provide one of the following:",
            "",
            "```text",
            "1. a uniform local/mesoscopic remainder theorem, plus a global far-tail saddle theorem",
            "2. a finite collar to K, plus an analytic tail theorem from K+1 onward",
            "```",
            "",
            "Integration:",
            "",
            "```text",
            "outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md",
            "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md",
            "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
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
    parser.add_argument("--activation-depth", type=int, default=21)
    parser.add_argument("--t-values", type=parse_t_values, default=parse_t_values("25,50,100"))
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    sign_scout_json = args.sign_scout_json if args.sign_scout_json.is_absolute() else REPO_ROOT / args.sign_scout_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    target = build_target(sign_scout_json, args.activation_depth, args.t_values)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(target, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(target, note)
    print(
        "wrote Jensen-window PF negative-lambda uniform remainder target: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
