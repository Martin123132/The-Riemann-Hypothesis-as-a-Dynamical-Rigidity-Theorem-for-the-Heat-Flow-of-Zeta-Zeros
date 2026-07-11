#!/usr/bin/env python3
"""Build the signed Gaussian perturbation matrix for the negative-lambda route."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from pathlib import Path
import json
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SIGN_SCOUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_phi_taylor_cone_entry_sign_scout.json"
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md"


@dataclass(frozen=True)
class PerturbationDiagnostics:
    ratio_a: str
    ratio_b: str
    ratio_c: str
    upper_wall_correction: str
    deficit_coefficient: str
    monotone_wall_correction: str
    certified_taylor_signs: int
    fixed_k_activation_depth: int
    fixed_k_activation_t_lower: str


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def arb_positive(value: flint.arb) -> bool:
    return bool(value > 0 and not value.contains(0))


def arb_negative(value: flint.arb) -> bool:
    return bool(value < 0 and not value.contains(0))


def build_diagnostics(sign_scout_json: Path, activation_depth: int) -> PerturbationDiagnostics:
    scout = load_json(sign_scout_json)
    ratios = scout["ratio_enclosures"]
    a = flint.arb(ratios["a=c2/c0"])
    b = flint.arb(ratios["b=c4/c0"])
    c = flint.arb(ratios["c=c6/c0"])
    upper_wall = 2 * b - a * a
    deficit = -upper_wall
    monotone_wall = 2 * (a**3 - 3 * a * b + 3 * c)
    if not arb_negative(upper_wall):
        raise RuntimeError(f"upper-wall correction is not certified negative: {upper_wall}")
    if not arb_positive(deficit):
        raise RuntimeError(f"deficit coefficient is not certified positive: {deficit}")
    if not arb_positive(monotone_wall):
        raise RuntimeError(f"monotone-wall correction is not certified positive: {monotone_wall}")

    sign_rows = scout.get("sign_combinations", [])
    certified = sum(1 for row in sign_rows if row.get("certified_sign") in {"positive", "negative"})
    # Leading fixed-k bound D2/T^2 <= 2/(3*(2*K+1)).
    activation = (deficit * flint.arb(3 * (2 * activation_depth + 1)) / flint.arb(2)).sqrt()
    return PerturbationDiagnostics(
        ratio_a=a.str(40),
        ratio_b=b.str(40),
        ratio_c=c.str(40),
        upper_wall_correction=upper_wall.str(40),
        deficit_coefficient=deficit.str(40),
        monotone_wall_correction=monotone_wall.str(40),
        certified_taylor_signs=certified,
        fixed_k_activation_depth=activation_depth,
        fixed_k_activation_t_lower=activation.str(40),
    )


def build_artifact(sign_scout_json: Path, activation_depth: int) -> dict:
    diagnostics = build_diagnostics(sign_scout_json, activation_depth)
    summary = {
        "matrix_rows": 8,
        "ready_to_apply_rows": 0,
        "certified_taylor_signs": diagnostics.certified_taylor_signs,
        "fixed_k_activation_estimates": 1,
        "live_routes": 2,
        "rejected_templates": 1,
        "target_closing": False,
        "main_finding": (
            "The surviving Gaussian-baseline route is a signed or tilted perturbation, not a positive "
            "Gaussian mixture: the certified Phi Taylor signs give a positive leading deficit "
            "D2=a^2-2*b and positive monotone correction 2*(a^3-3*a*b+3*c), but the fixed-k "
            "expansion still needs uniform remainder control or a finite-collar tail theorem before it "
            "can close the bounded log-curvature target."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix",
        "date": "2026-07-06",
        "status": "finite theorem-search diagnostic",
        "source_phi_taylor_sign_scout": "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md",
        "source_gaussian_curvature_matrix": "outputs/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.md",
        "source_bounded_log_curvature_target": "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md",
        "sign_scout_json": sign_scout_json.relative_to(REPO_ROOT).as_posix(),
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.py",
        "proof_boundary": (
            "Finite theorem-search diagnostic only. It packages the fixed-k signed Gaussian perturbation "
            "route and its certified local Taylor signs, but it does not prove a uniform asymptotic "
            "remainder, does not prove the bounded log-curvature theorem, does not prove cone entry, "
            "and does not prove Lambda <= 0."
        ),
        "matrix_rows": [
            {
                "id": "nlsgp_01_gaussian_baseline",
                "role": "exact_baseline",
                "claim": "The pure Gaussian raw-moment baseline has x_k=1 and B_k=0.",
                "proof_boundary": "Exact baseline only; it is degenerate for strict cone entry.",
            },
            {
                "id": "nlsgp_02_signed_taylor_perturbation_model",
                "role": "formal_model",
                "claim": "For Phi(u)=c0*(1+a*u^2+b*u^4+c*u^6+...), the fixed-k perturbation has log x_k=(2*b-a^2)/T^2+O_k(T^-3).",
                "proof_boundary": "Formal fixed-k model only; not a uniform tail theorem.",
            },
            {
                "id": "nlsgp_03_positive_deficit_sign",
                "role": "certified_local_sign",
                "claim": "The certified sign 2*b-a^2<0 gives positive leading Gaussian deficit D2=a^2-2*b>0.",
                "proof_boundary": "Certified local Taylor sign only.",
            },
            {
                "id": "nlsgp_04_monotone_deficit_sign",
                "role": "certified_local_sign",
                "claim": "The certified sign 2*(a^3-3*a*b+3*c)>0 gives fixed-k increasing contractions, equivalently decreasing B_k to the next order.",
                "proof_boundary": "Certified local Taylor sign only.",
            },
            {
                "id": "nlsgp_05_fixed_k_activation_scale",
                "role": "finite_depth_estimate",
                "claim": "Ignoring remainders, D2/T^2<=2/(3*(2*K+1)) gives a finite-depth activation scale for a chosen K.",
                "proof_boundary": "Leading-term estimate only; not a rigorous remainder bound.",
            },
            {
                "id": "nlsgp_06_uniform_tail_gap",
                "role": "open_gap",
                "claim": "The fixed-k expansion cannot by itself prove the all-k tail because the target bound shrinks like 1/k.",
                "proof_boundary": "Gap statement only; the needed uniform theorem remains open.",
            },
            {
                "id": "nlsgp_07_positive_mixture_template_rejected",
                "role": "rejected_template",
                "claim": "Positive Gaussian scale mixtures are unsuitable for the upper-wall route because they push the deficit nonpositive.",
                "proof_boundary": "Rejects this proof template only; not a statement about the actual signed Phi perturbation.",
            },
            {
                "id": "nlsgp_08_live_signed_perturbation_route",
                "role": "live_route",
                "claim": "A signed or tilted Gaussian perturbation theorem with explicit uniform remainders could prove the bounded log-curvature target.",
                "proof_boundary": "Live route only; no such theorem is proved.",
            },
        ],
        "finite_diagnostics": asdict(diagnostics),
        "invariants": [
            "No row is ready_to_apply.",
            "The route is fixed-k/formal until uniform remainders are proved.",
            "Positive Gaussian scale mixtures remain rejected for the upper-wall proof template.",
            "The bounded log-curvature theorem remains open.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
        "summary": summary,
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["finite_diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda signed Gaussian perturbation matrix: "
        f"{summary['matrix_rows']} matrix rows, "
        f"{summary['certified_taylor_signs']} certified Taylor signs, "
        f"{summary['fixed_k_activation_estimates']} fixed-k activation estimates, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows, 0 issues"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Signed Gaussian Perturbation Matrix",
        "",
        "Date: 2026-07-06",
        "",
        "Status: finite theorem-search diagnostic. This is not a proof of the",
        "bounded log-curvature theorem, cone entry, Jensen-window PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix`.",
        "",
        "Proof boundary: this artifact packages the fixed-k signed Gaussian",
        "perturbation route. It does not prove a uniform asymptotic remainder",
        "or an all-`k` tail theorem.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Perturbation Route",
        "",
        "For:",
        "",
        "```text",
        "Phi(u)=c0*(1+a*u^2+b*u^4+c*u^6+...)",
        "lambda=-T",
        "```",
        "",
        "the fixed-k expansion gives:",
        "",
        "```text",
        "log x_k = (2*b-a^2)/T^2 + O_k(T^-3)",
        "B_k = -log(x_k) = (a^2-2*b)/T^2 + O_k(T^-3)",
        "log(x_(k+1)/x_k) = 2*(a^3-3*a*b+3*c)/T^3 + O_k(T^-4)",
        "```",
        "",
        "Certified signs:",
        "",
        "```text",
        f"a = {diagnostics['ratio_a']}",
        f"b = {diagnostics['ratio_b']}",
        f"c = {diagnostics['ratio_c']}",
        f"2*b-a^2 = {diagnostics['upper_wall_correction']}",
        f"D2=a^2-2*b = {diagnostics['deficit_coefficient']}",
        f"2*(a^3-3*a*b+3*c) = {diagnostics['monotone_wall_correction']}",
        "```",
        "",
        "Leading fixed-k activation estimate:",
        "",
        "```text",
        f"K = {diagnostics['fixed_k_activation_depth']}",
        f"T >= {diagnostics['fixed_k_activation_t_lower']}",
        "```",
        "",
        "This estimate ignores remainders and therefore cannot be used as a proof.",
        "It is a scale diagnostic for the signed perturbation route only.",
        "",
        "Integration:",
        "",
        "```text",
        "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md",
        "outputs/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.md",
        "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md",
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
    parser.add_argument("--sign-scout-json", type=Path, default=DEFAULT_SIGN_SCOUT_JSON)
    parser.add_argument("--activation-depth", type=int, default=21)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    sign_scout_json = args.sign_scout_json if args.sign_scout_json.is_absolute() else REPO_ROOT / args.sign_scout_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(sign_scout_json, args.activation_depth)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda signed Gaussian perturbation matrix: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
