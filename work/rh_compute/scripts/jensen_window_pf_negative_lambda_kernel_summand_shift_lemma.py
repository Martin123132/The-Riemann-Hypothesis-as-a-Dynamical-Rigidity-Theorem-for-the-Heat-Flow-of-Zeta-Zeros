#!/usr/bin/env python3
"""Build the exact Newman-kernel summand-shift lemma and a compact tail bound."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.md"
DEFAULT_PRECISION_BITS = 512
DEFAULT_T = 100
DEFAULT_K = 300
DEFAULT_V = "1.5"
DEFAULT_COMPACT_N_MAX = 20
DEFAULT_COMPACT_RATIO_CAP = "2.122e-29"


@dataclass(frozen=True)
class LemmaRow:
    id: str
    role: str
    readiness: str
    claim: str
    proof_boundary: str
    formula: str | None = None
    diagnostics: dict | None = None
    gap: str | None = None


def arb_text(value: flint.arb, digits: int = 70) -> str:
    return value.str(digits).replace("e", "E")


def arb_upper_text(value: flint.arb, digits: int = 70) -> str:
    rounded = value.upper().str(digits, radius=False)
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(rounded).next_plus(), "E")


def compact_ratio(n_int: int, k: int, T: int, V: flint.arb) -> tuple[flint.arb, flint.arb]:
    n = flint.arb(n_int)
    a = n.log() / 2
    base = (V - a) / V
    if not bool(base > 0 and base < 1):
        raise ValueError(f"n={n_int} does not meet 0<a_n<V")
    ratio = base ** (2 * k) * (flint.arb(T) * (2 * a * V - a**2)).exp() / n.sqrt()
    return a, ratio


def build_artifact() -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    V = flint.arb(DEFAULT_V)
    cap = flint.arb(DEFAULT_COMPACT_RATIO_CAP)
    compact_rows = []
    compact_sum = flint.arb(0)
    for n in range(2, DEFAULT_COMPACT_N_MAX + 1):
        a, ratio = compact_ratio(n, DEFAULT_K, DEFAULT_T, V)
        compact_sum += ratio
        compact_rows.append(
            {
                "n": n,
                "a_n_ball": arb_text(a, 50),
                "ratio_bound_ball": arb_text(ratio, 60),
                "ratio_bound_upper": arb_upper_text(ratio, 60),
            }
        )
    a21 = flint.arb(21).log() / 2
    if not bool(compact_sum < cap and a21 > V):
        raise RuntimeError("compact summand-shift bound did not certify")

    diagnostics = {
        "T": DEFAULT_T,
        "K": DEFAULT_K,
        "V": DEFAULT_V,
        "compact_n_min": 2,
        "compact_n_max": DEFAULT_COMPACT_N_MAX,
        "compact_summand_count": DEFAULT_COMPACT_N_MAX - 1,
        "compact_ratio_rows": compact_rows,
        "compact_ratio_sum_ball": arb_text(compact_sum, 70),
        "compact_ratio_sum_upper": arb_upper_text(compact_sum, 70),
        "compact_ratio_cap": DEFAULT_COMPACT_RATIO_CAP,
        "compact_ratio_below_cap": True,
        "a_21_ball": arb_text(a21, 50),
        "a_21_above_V": True,
        "first_far_only_summand": 21,
    }
    rows = [
        LemmaRow(
            id="ksl_01_kernel_summand_definition",
            role="exact_input",
            readiness="available_exact",
            claim="Write Phi(u)=sum_(n>=1) phi_n(u) with the standard positive Newman-kernel summands.",
            formula="phi_n(u)=(2*pi^2*n^4*exp(9u)-3*pi*n^2*exp(5u))*exp(-pi*n^2*exp(4u))",
            proof_boundary="Exact kernel decomposition only.",
        ),
        LemmaRow(
            id="ksl_02_exact_summand_shift",
            role="exact_lemma",
            readiness="available_exact",
            claim="Every summand is a translated copy of the first: phi_n(u)=n^(-1/2)*phi_1(u+a_n), a_n=(log n)/2.",
            formula="exp(4*a_n)=n^2, n^(-1/2)*exp(9*a_n)=n^4, n^(-1/2)*exp(5*a_n)=n^2",
            proof_boundary="Exact algebra for every integer n>=1 and real u.",
        ),
        LemmaRow(
            id="ksl_03_exact_moment_shift",
            role="exact_reduction",
            readiness="available_exact",
            claim="After v=u+a_n, the nth negative-lambda moment is an explicitly shifted first-summand integral.",
            formula="M_(k,n)(-T)=2*n^(-1/2)*integral_(a_n)^infinity (v-a_n)^(2k)*exp(-T*(v-a_n)^2)*phi_1(v)dv",
            proof_boundary="Exact change of variables only; no tail estimate.",
        ),
        LemmaRow(
            id="ksl_04_relative_integrand_ratio",
            role="exact_reduction",
            readiness="available_exact",
            claim="Relative to the n=1 integrand at the same v, the shifted ratio has a closed positive form.",
            formula="rho_(n,k,T)(v)=n^(-1/2)*(1-a_n/v)^(2k)*exp(T*(2*a_n*v-a_n^2))",
            proof_boundary="Exact for v>a_n; phi_1 cancels pointwise.",
        ),
        LemmaRow(
            id="ksl_05_compact_monotonicity",
            role="exact_lemma",
            readiness="available_exact",
            claim="For fixed n,k,T, rho is increasing in v, and for fixed n,T,v it decreases in k.",
            formula="d_v log rho=2*k*a_n/(v*(v-a_n))+2*T*a_n>0; rho_(k+1)/rho_k=(1-a_n/v)^2<1",
            proof_boundary="Exact compact comparison for n>=2, T>0, and v>a_n.",
        ),
        LemmaRow(
            id="ksl_06_t100_k300_compact_certificate",
            role="interval_certificate",
            readiness="interval_validated",
            claim="For T=100, k>=300, and v<=3/2, the sum of all shifted n=2..20 integrands is below 2.122e-29 times the n=1 integrand.",
            formula="sum_(n=2)^20 rho_(n,k,100)(v)<=sum_(n=2)^20 rho_(n,300,100)(3/2)<2.122e-29",
            diagnostics={
                "compact_ratio_sum_upper": diagnostics["compact_ratio_sum_upper"],
                "compact_ratio_cap": diagnostics["compact_ratio_cap"],
                "summand_count": diagnostics["compact_summand_count"],
            },
            proof_boundary="Complete compact shifted-v contribution only; the v>=3/2 far tail is excluded.",
        ),
        LemmaRow(
            id="ksl_07_far_tail_partition",
            role="exact_reduction",
            readiness="available_exact",
            claim="Since a_21>3/2, every n>=21 shifted integral lies wholly in the v>=3/2 far-tail region.",
            formula="a_n=(log n)/2>=a_21>3/2 for n>=21",
            diagnostics={"a_21_ball": diagnostics["a_21_ball"]},
            proof_boundary="Exact support partition only; no far-tail magnitude bound.",
        ),
        LemmaRow(
            id="ksl_08_fixed_t_saddle_handoff",
            role="open_requirement",
            readiness="not_ready_to_apply",
            claim="Close the lambda=-100, k>=300 tail by proving the dominant n=1 adjacent-wall inequality and bounding the complete v>=3/2 shifted far tail below its stability margin.",
            gap="A uniform dominant-saddle inequality and a far-tail perturbation budget are still missing.",
            proof_boundary="Open theorem handoff; not cone entry, Jensen-window PF-infinity, RH, or Lambda <= 0.",
        ),
    ]
    summary = {
        "lemma_rows": len(rows),
        "available_exact_rows": sum(row.readiness == "available_exact" for row in rows),
        "interval_validated_rows": sum(row.readiness == "interval_validated" for row in rows),
        "open_requirement_rows": sum(row.readiness == "not_ready_to_apply" for row in rows),
        "compact_summand_count": diagnostics["compact_summand_count"],
        "compact_ratio_sum_upper": diagnostics["compact_ratio_sum_upper"],
        "first_far_only_summand": diagnostics["first_far_only_summand"],
        "ready_to_apply_rows": 0,
        "main_finding": (
            "The exact identity phi_n(u)=n^(-1/2)phi_1(u+(log n)/2) converts the full "
            "Newman kernel into shifted copies of one dominant summand. For lambda=-100, "
            "k>=300, and shifted variable v<=3/2, Arb bounds the complete n=2..20 "
            "relative contribution below 2.122e-29; n>=21 begins entirely beyond the split. "
            "The remaining tail theorem is now the dominant n=1 adjacent-wall inequality plus "
            "a v>=3/2 far-tail perturbation bound."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_kernel_summand_shift_lemma",
        "status": "exact kernel-shift lemma with compact interval certificate",
        "date": "2026-07-10",
        "proof_boundary": (
            "This artifact proves an exact all-n summand-shift identity and a complete compact "
            "relative bound for T=100, k>=300. It does not bound the v>=3/2 far tail, does not "
            "prove the dominant n=1 adjacent wall, does not prove cone entry, and does not prove "
            "RH or Lambda <= 0."
        ),
        "source_kernel_log_concavity": "outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md",
        "source_t1156_counterexample": "outputs/jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.md",
        "source_k300_audit": "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
        "source_cone_entry_target": "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.py",
        "parameters": {
            "T": DEFAULT_T,
            "K": DEFAULT_K,
            "V": DEFAULT_V,
            "compact_n_range": [2, DEFAULT_COMPACT_N_MAX],
            "precision_bits": DEFAULT_PRECISION_BITS,
        },
        "diagnostics": diagnostics,
        "summary": summary,
        "rows": [asdict(row) for row in rows],
    }


def write_note(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["diagnostics"]
    lines = [
        "# Jensen-Window PF Negative-Lambda Kernel Summand-Shift Lemma",
        "",
        "Date: 2026-07-10",
        "",
        "Status: exact kernel-shift lemma with compact interval certificate. This",
        "is not a proof of cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_kernel_summand_shift_lemma`.",
        "",
        "Machine-readable result:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.json",
        "```",
        "",
        "Generator and checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        "validated Jensen-window PF negative-lambda kernel summand-shift lemma: 8 rows, 0 issues, 6 exact rows, 1 compact interval row, 1 open far-tail row, 0 ready-to-apply rows",
        "```",
        "",
        "## Exact Shift",
        "",
        "Write the positive Newman-kernel summands as",
        "",
        "```text",
        "phi_n(u)=(2*pi^2*n^4*exp(9u)-3*pi*n^2*exp(5u))*exp(-pi*n^2*exp(4u)).",
        "```",
        "",
        "For `a_n=(log n)/2`, direct substitution gives",
        "",
        "```text",
        "phi_n(u)=n^(-1/2)*phi_1(u+a_n).",
        "```",
        "",
        "Therefore the nth contribution to the raw moment at `lambda=-T` is",
        "",
        "```text",
        "M_(k,n)(-T)=2*n^(-1/2)*integral_(a_n)^infinity",
        "  (v-a_n)^(2k)*exp(-T*(v-a_n)^2)*phi_1(v)dv.",
        "```",
        "",
        "Relative to the first-summand integrand at the same `v`,",
        "",
        "```text",
        "rho_(n,k,T)(v)=n^(-1/2)*(1-a_n/v)^(2k)",
        "                   *exp(T*(2*a_n*v-a_n^2)).",
        "d_v log rho=2*k*a_n/(v*(v-a_n))+2*T*a_n>0.",
        "```",
        "",
        "Thus `rho` increases in `v` and decreases in `k`.",
        "",
        "## Compact Certificate",
        "",
        "At `T=100`, `k>=300`, and `v<=3/2`, only `n=2..20` can occur before",
        "the split because `a_21>3/2`. Arb certifies",
        "",
        "```text",
        f"sum_(n=2)^20 rho_(n,300,100)(3/2) = {diagnostics['compact_ratio_sum_ball']}",
        f"upper bound < {diagnostics['compact_ratio_cap']}",
        "```",
        "",
        "Hence the complete shifted `n=2..20` compact contribution is below",
        "`2.122e-29` times the first-summand moment for every `k>=300`.",
        "",
        "## Remaining Tail",
        "",
        "The revised `lambda=-100` tail theorem now has two explicit pieces:",
        "",
        "```text",
        "1. prove the adjacent-k wall for the dominant n=1 moment when k>=300;",
        "2. bound all shifted contributions on v>=3/2 below the resulting margin.",
        "```",
        "",
        "The k300 Arb cache supplies the finite collar through `k=300`; this lemma",
        "does not yet supply either infinite-tail piece.",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
        "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "outputs/jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.md",
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "```",
        "",
        "Summary:",
        "",
        summary["main_finding"],
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print(
        "validated Jensen-window PF negative-lambda kernel summand-shift lemma: "
        "8 rows, 0 issues, 6 exact rows, 1 compact interval row, "
        "1 open far-tail row, 0 ready-to-apply rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
