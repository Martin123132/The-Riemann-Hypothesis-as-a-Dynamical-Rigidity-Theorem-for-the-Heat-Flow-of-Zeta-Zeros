#!/usr/bin/env python3
"""Prove eventual contiguous order-four positivity at lambda zero."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import itertools
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.md"
)
PREFIX_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/arb_xi_lambda0_order4_prefix_certificate.json"
)
PREFIX_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "arb_xi_lambda0_order4_prefix_coefficients_n0_n506_bits24576.jsonl"
)
PREFIX_LAST_N = 500
NEEDED_MAX_K = PREFIX_LAST_N + 6
ASYMPTOTIC_ORDER = 6


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


def series_multiply(left: list[sp.Expr], right: list[sp.Expr]) -> list[sp.Expr]:
    return [
        sp.expand(sum(left[k] * right[n - k] for k in range(n + 1)))
        for n in range(ASYMPTOTIC_ORDER + 1)
    ]


def exp_log_series(backward_shift: int, symbols: tuple[sp.Symbol, ...]) -> list[sp.Expr]:
    log_coefficients = [sp.Integer(0)] * (ASYMPTOTIC_ORDER + 1)
    for m in range(2, ASYMPTOTIC_ORDER + 2):
        log_coefficients[m - 1] = -(
            symbols[m - 2] * sp.Integer(backward_shift) ** m
        )
    exponential = [sp.Integer(0)] * (ASYMPTOTIC_ORDER + 1)
    exponential[0] = sp.Integer(1)
    for n in range(1, ASYMPTOTIC_ORDER + 1):
        exponential[n] = sp.expand(
            sum(
                k * log_coefficients[k] * exponential[n - k]
                for k in range(1, n + 1)
            )
            / n
        )
    return exponential


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def determinant_cancellation() -> dict:
    symbols = sp.symbols("G2:8")
    entry_series = {
        backward_shift: exp_log_series(backward_shift, symbols)
        for backward_shift in range(7)
    }
    determinant = [sp.Integer(0)] * (ASYMPTOTIC_ORDER + 1)
    for permutation in itertools.permutations(range(4)):
        product = [sp.Integer(1)] + [sp.Integer(0)] * ASYMPTOTIC_ORDER
        for row in range(4):
            backward_shift = 6 - row - permutation[row]
            product = series_multiply(product, entry_series[backward_shift])
        sign = permutation_sign(permutation)
        determinant = [
            sp.expand(current + sign * contribution)
            for current, contribution in zip(determinant, product)
        ]
    factored = [sp.factor(value) for value in determinant]
    expected = [sp.Integer(0)] * ASYMPTOTIC_ORDER + [768 * symbols[0] ** 6]
    if factored != expected:
        raise RuntimeError(f"unexpected determinant expansion: {factored}")
    return {
        "matrix_entry": (
            "K_(i,j)(h)=exp(-sum_(m>=2) G_m*h^(m-1)*(6-i-j)^m)"
        ),
        "truncation": "m=2,...,7 and powers h^0,...,h^6",
        "coefficients_h0_through_h6": [str(value) for value in factored],
        "first_nonzero_term": "768*G2^6*h^6",
        "independence": "G3,...,G7 cancel from the first nonzero coefficient",
        "permutations_checked": 24,
    }


def finite_prefix() -> dict:
    source = json.loads(PREFIX_SOURCE.read_text(encoding="utf-8"))
    if source.get("kind") != "arb_xi_lambda0_order4_prefix_certificate":
        raise RuntimeError("lambda-zero prefix source has the wrong kind")
    summary = source.get("summary", {})
    if summary.get("positive_H4_rows") != PREFIX_LAST_N + 1:
        raise RuntimeError("lambda-zero prefix source has incomplete H4 coverage")
    if summary.get("positive_stable_margin_rows") != PREFIX_LAST_N + 1:
        raise RuntimeError("lambda-zero prefix source has incomplete stable margins")
    finite = source.get("finite", {})
    if finite.get("n_range") != [0, PREFIX_LAST_N]:
        raise RuntimeError("lambda-zero prefix source has the wrong range")
    return {
        "lambda": "0",
        "n_range": [0, PREFIX_LAST_N],
        "coefficient_range": [0, NEEDED_MAX_K],
        "all_H4_positive": True,
        "all_stable_margins_positive": True,
        "rows": finite.get("rows", []),
        "minimum_margin_n": finite.get("minimum_margin_n"),
        "minimum_margin_ball": finite.get("minimum_margin_ball"),
        "minimum_normalized_n": finite.get("minimum_normalized_n"),
        "minimum_normalized_ball": finite.get("minimum_normalized_ball"),
        "source_sha256": {
            str(PREFIX_SOURCE.relative_to(REPO_ROOT)).replace("\\", "/"): sha256(
                PREFIX_SOURCE
            ),
            str(PREFIX_CACHE.relative_to(REPO_ROOT)).replace("\\", "/"): sha256(
                PREFIX_CACHE
            ),
        },
    }


def xi_normalization() -> dict:
    return {
        "newman_identity": "H_0(x)=xi((1+i*x)/2)/8",
        "xi_taylor": "xi(1/2+z)=sum_(k>=0) gamma(k)*z^(2k)/k!",
        "moment_definition": (
            "A_k(0)=2*k!/(2k)!*integral_0^infinity Phi(u)*u^(2k)du"
        ),
        "coefficient_identity": "A_k(0)=gamma(k)/4^(k+1)",
        "determinant_identity": (
            "H_(4,n)[A(0)]=4^(-4*n-16)*H_(4,n)[gamma]"
        ),
        "sign_consequence": "the A and gamma contiguous order-four determinants have the same sign",
    }


def published_ratio_input() -> dict:
    return {
        "source": (
            "Griffin-Ono-Rolen-Thorner-Tripp-Wagner, "
            "Jensen Polynomials for the Riemann Xi Function"
        ),
        "arxiv": "https://arxiv.org/abs/1910.01227",
        "journal": "Advances in Mathematics 397 (2022), 108186",
        "theorem": "Theorem 2.1, especially parts (1) and (2)",
        "uniformizer": (
            "Delta(M)^2=(1/2)*(1-gamma(M-2)*gamma(M)/gamma(M-1)^2)"
        ),
        "exact_ratio": (
            "log(gamma(M-j)/gamma(M))="
            "-sum_(m>=1) G_m(M)*Delta(M)^(2*m-2)*j^m"
        ),
        "range": "integers 1<=j<M",
        "limits": "Delta(M)~1/sqrt(2*M) and G_m(M)->2^(m-1)/(m*(m-1))",
        "needed_limit": "G_2(M)->1",
        "tail_control": "|G_m(M)|<<_C(2*C)^m uniformly in m and M",
    }


def eventual_theorem(cancellation: dict) -> dict:
    if cancellation["first_nonzero_term"] != "768*G2^6*h^6":
        raise RuntimeError("symbolic cancellation source changed")
    return {
        "index_choice": "M=n+6, j=6-i-j_column, h=Delta(M)^2",
        "affine_factor": (
            "the m=1 ratio term factors from rows and columns as exp(-12*G_1(M))"
        ),
        "gamma_asymptotic": (
            "H_(4,n)[gamma]=gamma(M)^4*exp(-12*G_1(M))*"
            "(768*G_2(M)^6*Delta(M)^12+o(Delta(M)^12))"
        ),
        "A_asymptotic": (
            "H_(4,n)[A(0)]=4^(-4*n-16)*gamma(M)^4*exp(-12*G_1(M))*"
            "(768*G_2(M)^6*Delta(M)^12+o(Delta(M)^12))"
        ),
        "strict_sign": "there exists N_H4 such that H_(4,n)(0)>0 for every n>=N_H4",
        "threshold_effective_here": False,
    }


def build_artifact() -> dict:
    normalization = xi_normalization()
    published = published_ratio_input()
    cancellation = determinant_cancellation()
    prefix = finite_prefix()
    eventual = eventual_theorem(cancellation)
    rows = [
        CertificateRow(
            id="co4l0epc_01_xi_normalization",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The lambda-zero Newman moments are a positive geometric rescaling of the standard Xi Taylor coefficients.",
            formula=normalization["coefficient_identity"],
            proof_boundary="Exact Fourier/Taylor normalization only.",
            diagnostics=normalization,
        ),
        CertificateRow(
            id="co4l0epc_02_hankel_scale",
            role="exact_sign_equivalence",
            readiness="ready_to_apply",
            claim="The geometric coefficient rescaling preserves the contiguous order-four sign.",
            formula=normalization["determinant_identity"],
            proof_boundary="Contiguous four-by-four determinant only.",
        ),
        CertificateRow(
            id="co4l0epc_03_published_ratio_input",
            role="published_theorem_input",
            readiness="ready_to_apply",
            claim="The published Xi ratio theorem supplies an exact convergent backward-ratio expansion with the required graded powers of Delta.",
            formula=published["exact_ratio"],
            proof_boundary="Uses the cited published theorem; no Jensen-hyperbolicity-to-Hankel implication is assumed.",
            diagnostics=published,
        ),
        CertificateRow(
            id="co4l0epc_04_symbolic_cancellation",
            role="exact_symbolic_lemma",
            readiness="ready_to_apply",
            claim="The normalized four-by-four determinant vanishes through order h^5 and has a universal positive first coefficient.",
            formula=cancellation["first_nonzero_term"],
            proof_boundary="Exact formal determinant expansion through the first nonzero order.",
            diagnostics=cancellation,
        ),
        CertificateRow(
            id="co4l0epc_05_eventual_positivity",
            role="exact_theorem_composition",
            readiness="ready_to_apply",
            claim="The actual lambda-zero contiguous order-four determinant is strictly positive for every sufficiently large shift.",
            formula=eventual["strict_sign"],
            proof_boundary="Eventual theorem only; the threshold N_H4 is not made explicit here.",
            diagnostics=eventual,
        ),
        CertificateRow(
            id="co4l0epc_06_arb_prefix",
            role="interval_certificate",
            readiness="ready_to_apply",
            claim="Rigorous coefficient balls prove the lambda-zero contiguous order-four sign on the complete available prefix.",
            formula="H_(4,n)(0)>0 for 0<=n<=500",
            proof_boundary="Finite prefix only.",
            diagnostics={
                "n_range": prefix["n_range"],
                "minimum_margin_n": prefix["minimum_margin_n"],
                "minimum_margin_ball": prefix["minimum_margin_ball"],
                "minimum_normalized_n": prefix["minimum_normalized_n"],
                "minimum_normalized_ball": prefix["minimum_normalized_ball"],
            },
        ),
        CertificateRow(
            id="co4l0epc_07_finite_gap_reduction",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim="All-shift lambda-zero order-four positivity is reduced to an effective finite splice between the certified prefix and the eventual theorem.",
            formula="it remains to prove 501<=n<N_H4, or prove N_H4<=501",
            proof_boundary="The finite set is not identified until N_H4 is effectivized.",
        ),
        CertificateRow(
            id="co4l0epc_08_effective_splice_target",
            role="analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Extract explicit constants from the Xi ratio expansion that force the determinant remainder below its 768*G_2^6 main term.",
            formula="find explicit N_H4 and certify every remaining n>=501 below it",
            proof_boundary="Open effective remainder bound; no all-shift lambda-zero theorem is asserted.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate",
        "date": "2026-07-13",
        "status": (
            "exact eventual lambda-zero order-four positivity theorem plus a rigorous "
            "501-shift prefix, with the effective splice still open"
        ),
        "proof_boundary": (
            "This artifact proves eventual H_(4,n)(0)>0 and the finite prefix "
            "0<=n<=500. It is not a proof of all-shift order four at lambda zero, "
            "forward order-four invariance, PF-infinity, RH, or Lambda<=0."
        ),
        "normalization": normalization,
        "published_ratio_input": published,
        "symbolic_cancellation": cancellation,
        "eventual_theorem": eventual,
        "finite_prefix": prefix,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_rows": len(rows) - 1,
            "symbolic_coefficients_checked": ASYMPTOTIC_ORDER + 1,
            "finite_prefix_rows": PREFIX_LAST_N + 1,
            "eventual_positivity_theorems": 1,
            "open_effective_splices": 1,
            "all_shift_lambda0_order4_proved": False,
        },
        "sources": [
            "https://arxiv.org/abs/1910.01227",
            "outputs/arb_xi_lambda0_order4_prefix_certificate.md",
            "work/rh_compute/results/arb_xi_lambda0_order4_prefix_certificate.json",
            "work/rh_compute/results/arb_xi_lambda0_order4_prefix_coefficients_n0_n506_bits24576.jsonl",
            "outputs/formal_core.md",
        ],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.py"
        ),
        "remaining_target": "effectivize N_H4 and close the finite splice above n=500",
    }


def write_note(path: Path, artifact: dict) -> None:
    normalization = artifact["normalization"]
    published = artifact["published_ratio_input"]
    cancellation = artifact["symbolic_cancellation"]
    eventual = artifact["eventual_theorem"]
    prefix = artifact["finite_prefix"]
    lines = [
        "# Jensen-Window PF Compound Order-Four Lambda-Zero Eventual Positivity",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact eventual lambda-zero order-four positivity theorem plus a",
        "rigorous `0<=n<=500` prefix, with the effective splice still open. This is",
        "not a proof of all-shift order four at lambda zero, PF-infinity, RH, or",
        "`Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.py",
        "```",
        "",
        "## Exact Xi Normalization",
        "",
        "The standard identities",
        "",
        "```text",
        normalization["newman_identity"],
        normalization["xi_taylor"],
        normalization["moment_definition"],
        "```",
        "",
        "give",
        "",
        "```text",
        normalization["coefficient_identity"],
        normalization["determinant_identity"],
        "```",
        "",
        "so the lambda-zero `A_k` and standard `gamma(k)` determinants have the",
        "same sign.",
        "",
        "## Published Ratio Input",
        "",
        "Theorem 2.1 of Griffin, Ono, Rolen, Thorner, Tripp, and Wagner,",
        "[Jensen Polynomials for the Riemann Xi Function](https://arxiv.org/abs/1910.01227),",
        "gives, for `1<=j<M`,",
        "",
        "```text",
        published["exact_ratio"],
        published["limits"],
        published["tail_control"],
        "```",
        "",
        "This use is directly through the Xi coefficient asymptotics. No claim that",
        "Jensen hyperbolicity by itself implies the signed Hankel sign is used.",
        "",
        "## Universal Determinant Term",
        "",
        "Set `M=n+6`, `h=Delta(M)^2`, and remove the affine `G_1` factor from",
        "rows and columns. The remaining matrix has entries",
        "",
        "```text",
        cancellation["matrix_entry"],
        "```",
        "",
        "Exact 24-permutation truncated-series algebra gives",
        "",
        "```text",
        "[h^0,...,h^6] det(K) =",
        str(cancellation["coefficients_h0_through_h6"]),
        cancellation["first_nonzero_term"],
        "```",
        "",
        "All `G_3,...,G_7` contributions cancel from the first nonzero term. The",
        "published tail bound controls `m>=8`, while `G_2(M)->1`. Consequently",
        "",
        "```text",
        eventual["gamma_asymptotic"],
        eventual["strict_sign"],
        "```",
        "",
        "This is an actual eventual positivity theorem for the zeta/Xi sequence,",
        "not merely numerical evidence.",
        "",
        "## Rigorous Prefix",
        "",
        "A direct 24,576-bit Arb expansion of `xi(1/2+z)` through `A_506(0)`",
        "give",
        "",
        "```text",
        "H_(4,n)(0)>0 for every 0<=n<=500",
        f"smallest recorded stable margin occurs at n={prefix['minimum_margin_n']}",
        f"F_500={prefix['minimum_margin_ball']}",
        f"F_500/G_501^2={prefix['minimum_normalized_ball']}",
        "```",
        "",
        "All 501 raw `H_4` balls and all 501 stable `F_n` balls are strictly",
        "positive. The full interval rows and source hashes are stored in the JSON",
        "artifact.",
        "",
        "## Remaining Splice",
        "",
        "The current result is",
        "",
        "```text",
        "finite theorem:   H_(4,n)(0)>0 for 0<=n<=500",
        "eventual theorem: H_(4,n)(0)>0 for n>=N_H4 for some finite N_H4",
        "open splice:      make N_H4 explicit and cover 501<=n<N_H4",
        "```",
        "",
        "The next analytic job is to extract explicit constants from the convergent",
        "ratio expansion, bound the normalized determinant remainder against",
        "`768*G_2^6`, and then run Arb only over the resulting finite gap.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_forward_flow_reduction.md",
        "outputs/jensen_window_pf_compound_order4_m100_entry_certificate.md",
        "outputs/arb_xi_lambda0_order4_prefix_certificate.md",
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
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print(
        "built lambda-zero order-four eventual positivity certificate: "
        f"{artifact['summary']['finite_prefix_rows']} prefix rows, "
        "1 eventual theorem, 1 open effective splice"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
