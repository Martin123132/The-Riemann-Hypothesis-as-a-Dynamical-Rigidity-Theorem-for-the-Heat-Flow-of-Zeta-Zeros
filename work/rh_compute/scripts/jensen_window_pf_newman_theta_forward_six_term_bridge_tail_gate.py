#!/usr/bin/env python3
"""Prove a six-term ordinary-theta tail theorem on the finite Newman bridge."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
from math import factorial
from pathlib import Path
import sys

try:
    import psutil

    psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
except Exception:
    pass


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR_DIR = REPO_ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR_DIR) not in sys.path:
    sys.path.insert(0, str(VENDOR_DIR))

import flint
from flint import arb
import sympy as sp

from jensen_window_pf_newman_theta_switch_defect_explicit_constant_gate import (
    forward_tail,
    serialize_positive,
)


STEM = "jensen_window_pf_newman_theta_forward_six_term_bridge_tail_gate"
DEFAULT_OUT = (
    REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
RESULT_ROOT = REPO_ROOT / "work" / "rh_compute" / "results"
SCRIPT_ROOT = REPO_ROOT / "work" / "rh_compute" / "scripts"
MODULAR_SOURCE = RESULT_ROOT / "jensen_window_pf_newman_theta_modular_blend_gate.json"
BASE_STEM = "jensen_window_pf_newman_theta_switch_defect_explicit_constant_gate"
BASE_RESULT = RESULT_ROOT / f"{BASE_STEM}.json"
BASE_BUILDER = SCRIPT_ROOT / f"{BASE_STEM}.py"
BASE_CHECKER = SCRIPT_ROOT / f"check_{BASE_STEM}.py"
CONTRACT_RESULT = (
    RESULT_ROOT
    / "jensen_window_pf_newman_theta_adaptive_modular_c1_remainder_contract.json"
)
DATE = "2026-07-25"
PRECISION_BITS = 256
T_CAP = Fraction(1, 5)
X_MIN = 38
X_MAX = 245
RETAINED_TERMS = 6
FIRST_OMITTED = RETAINED_TERMS + 1
WITNESS_FIRST_OMITTED = (6, 7, 8)
RATIO_TARGET = Fraction(1, 100_000_000_000)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def serialize_ball(value: arb) -> dict:
    if not value.is_finite():
        raise RuntimeError(f"nonfinite enclosure: {value}")
    return {
        "enclosure": value.str(70, more=True),
        "lower": value.lower().str(70),
        "upper": value.upper().str(70),
        "relative_accuracy_bits": value.rel_accuracy_bits(),
    }


def source_audit() -> dict:
    modular = json.loads(MODULAR_SOURCE.read_text(encoding="utf-8"))
    theta = modular.get("exact", {}).get("theta_summand", {})
    expected_definition = (
        "phi_n(u)=(2*pi^2*n^4*exp(9u)-3*pi*n^2*exp(5u))"
        "*exp(-pi*n^2*exp(4u))"
    )
    if theta.get("definition") != expected_definition:
        raise RuntimeError("ordinary theta summand definition drifted")
    if theta.get("kernel_series") != "Phi(u)=sum_(n>=1)phi_n(u)=Phi(-u)":
        raise RuntimeError("ordinary theta series marker drifted")
    positivity = (
        modular.get("exact", {})
        .get("strict_positivity", {})
        .get("positive_side")
    )
    if positivity != "phi_n(u)>0 for n>=1 and u>=0":
        raise RuntimeError("positive-half-line theta marker drifted")

    contract = json.loads(CONTRACT_RESULT.read_text(encoding="utf-8"))
    direct = contract.get("exact", {}).get("direct_c1_contract", {})
    expected_contract = {
        "partial_value": "J_(N,t)(x)=16*x^4*S_(N,t)(x)",
        "partial_derivative": (
            "J_(N,t)'(x)=64*x^3*S_(N,t)(x)+16*x^4*S_(N,t)'(x)"
        ),
    }
    for key, value in expected_contract.items():
        if direct.get(key) != value:
            raise RuntimeError(f"direct C1 contract drifted: {key}")
    return {
        "modular_source_sha256": file_hash(MODULAR_SOURCE),
        "base_result_sha256": file_hash(BASE_RESULT),
        "base_builder_sha256": file_hash(BASE_BUILDER),
        "base_checker_sha256": file_hash(BASE_CHECKER),
        "contract_result_sha256": file_hash(CONTRACT_RESULT),
        "theta_definition": expected_definition,
        "theta_series": theta["kernel_series"],
        "positive_half_line": positivity,
        "direct_c1_contract": expected_contract,
    }


def raw_forward_tail(moment: int, first: int) -> arb:
    if moment not in (0, 1):
        raise ValueError("this gate needs only moments zero and one")
    # P_0(X)=2X-3 has coefficient 1-norm five.
    return 5 * forward_tail(
        moment,
        kernel_power=2,
        beta=sp.Rational(9, 4),
        first=first,
    )


def witness_row(first: int) -> dict:
    e0 = raw_forward_tail(0, first)
    e1 = raw_forward_tail(1, first)
    x = arb(X_MAX)
    gamma = (-arb.pi() * x / 8).exp()
    value_error = 16 * x**4 * e0
    derivative_error = 64 * x**3 * e0 + 16 * x**4 * e1
    return {
        "retained_terms": first - 1,
        "first_omitted": first,
        "E0": serialize_positive(e0),
        "E1": serialize_positive(e1),
        "x_worst": X_MAX,
        "gamma_scale": serialize_positive(gamma),
        "J_error": serialize_positive(value_error),
        "J_prime_error": serialize_positive(derivative_error),
        "J_error_over_gamma": serialize_positive(value_error / gamma),
        "J_prime_error_over_gamma": serialize_positive(
            derivative_error / gamma
        ),
    }


def monotonicity_diagnostics() -> dict:
    rate = arb.pi() / 8
    slopes = {
        "value_x4_at_x38": rate + arb(4) / X_MIN,
        "derivative_x3_at_x38": rate + arb(3) / X_MIN,
        "derivative_x4_at_x38": rate + arb(4) / X_MIN,
    }
    if min(value.lower() for value in slopes.values()) <= 0:
        raise RuntimeError("gamma-normalized bridge multiplier is not increasing")
    theta = 1 - (arb(121) / 80) / (
        arb.pi() * FIRST_OMITTED**2
    )
    rho = (
        arb(2) / FIRST_OMITTED
        - arb.pi() * (2 * FIRST_OMITTED + 1)
    ).exp()
    if theta.lower() <= 0 or rho.upper() >= 1:
        raise RuntimeError("forward-tail denominator diagnostics failed")
    return {
        "spectral_rate_pi_over_8": serialize_positive(rate),
        "value_x4_log_slope_at_x38": serialize_positive(
            slopes["value_x4_at_x38"]
        ),
        "derivative_x3_log_slope_at_x38": serialize_positive(
            slopes["derivative_x3_at_x38"]
        ),
        "derivative_x4_log_slope_at_x38": serialize_positive(
            slopes["derivative_x4_at_x38"]
        ),
        "theta_K7": serialize_positive(theta),
        "rho_K7": serialize_positive(rho),
        "monotonicity_formula": (
            "d/dx log(exp(pi*x/8)*x^s)=pi/8+s/x>0, "
            "x>0, s in {3,4}"
        ),
    }


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    audit = source_audit()
    witnesses = [witness_row(first) for first in WITNESS_FIRST_OMITTED]
    by_first = {row["first_omitted"]: row for row in witnesses}
    target = arb(RATIO_TARGET.numerator) / RATIO_TARGET.denominator
    for field in ("J_error_over_gamma", "J_prime_error_over_gamma"):
        if arb(by_first[FIRST_OMITTED][field]["enclosure"]).upper() >= target:
            raise RuntimeError(f"six-term ratio target failed: {field}")
        if arb(by_first[FIRST_OMITTED - 1][field]["enclosure"]).lower() <= 1:
            raise RuntimeError(f"five-term failure guard drifted: {field}")
        if (
            arb(by_first[FIRST_OMITTED + 1][field]["enclosure"]).upper()
            >= arb(by_first[FIRST_OMITTED][field]["enclosure"]).lower()
        ):
            raise RuntimeError(f"tail decrease guard failed: {field}")
    monotonicity = monotonicity_diagnostics()

    tail_formula = (
        "E_p(K)<=5*p!*exp(1/80)*pi*K^2*exp(-pi*K^2)"
        "/[4*theta_K*(1-rho_K)]; "
        "theta_K=1-121/(80*pi*K^2), "
        "rho_K=exp(2/K-pi*(2K+1)), p in {0,1}"
    )
    rows = [
        GateRow(
            id="ntfstbtg_01_forward_half_line_partition",
            role="exact_identity",
            readiness="proved",
            claim=(
                "On the positive half-line, the Xi theta kernel is the "
                "ordinary positive summand series."
            ),
            formula=(
                "Phi(u)=sum_(n>=1)phi_n(u), phi_n(u)>0 for u>=0"
            ),
            proof_boundary=(
                "This is the stored theta identity; individual ordinary "
                "summands are not even."
            ),
        ),
        GateRow(
            id="ntfstbtg_02_six_term_transform",
            role="exact_definition",
            readiness="proved",
            claim=(
                "The first six ordinary summands define a legitimate "
                "finite value and first-derivative transform."
            ),
            formula=(
                "S_(6,t)^F(x)=sum_(n=1)^6 integral_0^infinity "
                "exp(tu^2)phi_n(u)cos(xu)du"
            ),
            proof_boundary=(
                "The residual is treated directly in L1; no evenness or "
                "integration-by-parts boundary cancellation is assumed."
            ),
        ),
        GateRow(
            id="ntfstbtg_03_raw_forward_tail_majorant",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "The omitted value and first-moment tails have explicit "
                "Gaussian arithmetic majorants uniform in 0<=t<=1/5."
            ),
            formula=tail_formula,
            proof_boundary=(
                "Uses P_0(X)=2X-3, coefficient norm five, "
                "u^p<=p!*exp(u), u^2<=exp(4u)/16, and a geometric "
                "sum from K onward."
            ),
        ),
        GateRow(
            id="ntfstbtg_04_direct_c1_error",
            role="exact_composition",
            readiness="proved",
            claim=(
                "Raw L1 and first-moment tails give direct J/J' errors "
                "without frequency integration by parts."
            ),
            formula=(
                "|J-J_6^F|<=16*x^4*E0(7); "
                "|J'-(J_6^F)'|<=64*x^3*E0(7)+16*x^4*E1(7)"
            ),
            proof_boundary=(
                "Triangle inequalities only; this is an omitted-tail upper "
                "bound, not retained first-jet separation."
            ),
        ),
        GateRow(
            id="ntfstbtg_05_bridge_endpoint_monotonicity",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "Both gamma-normalized direct error bounds increase with x "
                "on the finite bridge, so x=245 is the unique worst endpoint."
            ),
            formula=monotonicity["monotonicity_formula"],
            proof_boundary="Applies to the positive terms with powers x^3 and x^4.",
        ),
        GateRow(
            id="ntfstbtg_06_six_term_gamma_scale_theorem",
            role="exact_theorem",
            readiness="proved",
            claim=(
                "Six ordinary forward terms leave less than 10^-11 of the "
                "natural gamma scale in both J and J' on 38<=x<=245."
            ),
            formula=(
                "|J-J_6^F|<10^-11*exp(-pi*x/8); "
                "|J'-(J_6^F)'|<10^-11*exp(-pi*x/8)"
            ),
            proof_boundary=(
                "Uniform for 0<=t<=1/5; five retained terms explicitly fail "
                "this gamma-scale endpoint test."
            ),
            diagnostics=by_first[FIRST_OMITTED],
        ),
        GateRow(
            id="ntfstbtg_07_finite_bridge_reduction",
            role="exact_reduction",
            readiness="proved",
            claim=(
                "A retained interval certificate for 38<=x<=245 needs only "
                "six ordinary oscillatory integrals plus the displayed "
                "direct tail bars."
            ),
            formula=(
                "max(|J_6^F|,|(J_6^F)'|)"
                ">10^-11*exp(-pi*x/8) is a sufficient pointwise target"
            ),
            proof_boundary=(
                "This reduction does not assert that the retained "
                "disjunction holds."
            ),
        ),
        GateRow(
            id="ntfstbtg_08_retained_bridge_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim=(
                "Certify retained six-term first-jet separation on the "
                "remaining low-time and right-hand finite bridge rectangles."
            ),
            formula=(
                "[1/1035,1/155]x[38,69] union "
                "[1/1035,1/5]x[69,245]"
            ),
            proof_boundary=(
                "No retained interval cover, cofinal retained theorem, "
                "strict Laguerre positivity, Lambda<=0, RH, or Clay-prize "
                "proof is supplied."
            ),
        ),
    ]
    return {
        "kind": STEM,
        "date": DATE,
        "status": (
            "six-term ordinary-theta gamma-scale omitted-tail theorem "
            "on the finite Newman bridge"
        ),
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "t_cap_exact": "1/5",
            "x_interval": [X_MIN, X_MAX],
            "retained_terms": RETAINED_TERMS,
            "first_omitted": FIRST_OMITTED,
            "witness_first_omitted": list(WITNESS_FIRST_OMITTED),
            "ratio_target_exact": "1/100000000000",
        },
        "source_audit": audit,
        "tail_majorant": {
            "kernel_polynomial": "P_0(X)=2X-3",
            "kernel_coefficient_norm": 5,
            "beta": "9/4",
            "kernel_power": 2,
            "formula": tail_formula,
            "uniformity": "0<=t<=1/5",
        },
        "monotonicity": monotonicity,
        "witnesses": witnesses,
        "theorem": {
            "domain": "0<=t<=1/5, 38<=x<=245",
            "retained_transform": (
                "S_(6,t)^F=sum_(n=1)^6 integral_0^infinity "
                "exp(tu^2)phi_n(u)cos(xu)du"
            ),
            "value_error": (
                "|J_t-J_(6,t)^F|<10^-11*exp(-pi*x/8)"
            ),
            "first_derivative_error": (
                "|J_t'-(J_(6,t)^F)'|<10^-11*exp(-pi*x/8)"
            ),
        },
        "rows": [asdict(row) for row in rows],
        "proof_boundary": (
            "This artifact proves only direct omitted-tail upper bounds for "
            "the six-term ordinary forward theta transform on "
            "0<=t<=1/5 and 38<=x<=245. It does not certify retained "
            "first-jet separation on the finite bridge, prove the cofinal "
            "retained transition theorem, establish strict Laguerre "
            "positivity, prove Lambda<=0 or RH, or supply a Clay-prize proof."
        ),
    }


def render_note(artifact: dict) -> str:
    lines = [
        "# Newman Theta Forward Six-Term Bridge Tail Gate",
        "",
        f"Date: {DATE}",
        "",
        "Status: exact finite-bridge omitted-tail theorem. This is not a",
        "retained-separation theorem and not a proof of `Lambda<=0` or RH.",
        "",
        "## Why Ordinary Terms Work Here",
        "",
        "For `u>=0`,",
        "",
        "```text",
        "Phi(u)=sum_(n>=1)phi_n(u),",
        "phi_n(u)>0.",
        "```",
        "",
        "The ordinary summands are not even. Accordingly, this gate does not",
        "reuse the ninefold modular integration-by-parts estimate. It bounds",
        "the raw positive-half-line `L1` tail and its first moment directly.",
        "",
        "Define",
        "",
        "```text",
        "S_(6,t)^F(x)",
        " =sum_(n=1)^6 integral_0^infinity",
        "   exp(tu^2)phi_n(u)cos(xu)du.",
        "```",
        "",
        "## Explicit Tail",
        "",
        "For `K>=6`, `p in {0,1}`, and",
        "",
        "```text",
        "E_p(K)=sum_(n>=K) integral_0^infinity",
        "       u^p exp(u^2/5)phi_n(u)du,",
        "```",
        "",
        "the directed majorant is",
        "",
        "```text",
        artifact["tail_majorant"]["formula"],
        "```",
        "",
        "This uses the coefficient norm five of `P_0(X)=2X-3`. No",
        "cancellation or zero information enters the bound.",
        "",
        "## Direct C1 Contract",
        "",
        "With `K=7`,",
        "",
        "```text",
        "|J-J_6^F|<=16*x^4*E_0(7),",
        "|J'-(J_6^F)'|<=64*x^3*E_0(7)+16*x^4*E_1(7).",
        "```",
        "",
        "After division by `exp(-pi*x/8)`, every positive multiplier is",
        "increasing because",
        "",
        "```text",
        "d/dx log(exp(pi*x/8)*x^s)=pi/8+s/x>0.",
        "```",
        "",
        "Thus `x=245` is the worst endpoint throughout `38<=x<=245`.",
        "",
        "## Directed Witnesses",
        "",
        "| retained | first omitted | J error / gamma at 245 | J' error / gamma at 245 |",
        "|---:|---:|---:|---:|",
    ]
    for row in artifact["witnesses"]:
        lines.append(
            f"| {row['retained_terms']} | {row['first_omitted']} | "
            f"`{row['J_error_over_gamma']['enclosure']}` | "
            f"`{row['J_prime_error_over_gamma']['enclosure']}` |"
        )
    lines.extend(
        [
            "",
            "Five retained terms are an explicit failure guard. Six terms",
            "put both ratios below `10^-11`; seven terms give a much smaller",
            "cross-check.",
            "",
            "## Theorem",
            "",
            "For every `0<=t<=1/5` and `38<=x<=245`,",
            "",
            "```text",
            "|J_t-J_(6,t)^F|<10^-11*exp(-pi*x/8),",
            "|J_t'-(J_(6,t)^F)'|<10^-11*exp(-pi*x/8).",
            "```",
            "",
            "The next finite task is therefore a retained six-term interval",
            "cover on",
            "",
            "```text",
            "[1/1035,1/155] x [38,69]",
            "union [1/1035,1/5] x [69,245].",
            "```",
            "",
            "That retained cover and the cofinal retained theorem remain open.",
            "",
        ]
    )
    return "\n".join(lines)


def write_artifact(artifact: dict, out: Path, note: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    note.write_text(render_note(artifact), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    write_artifact(artifact, args.out, args.note)
    print(
        "wrote six-term ordinary-theta bridge tail gate: "
        f"{len(artifact['rows'])} rows, "
        f"{len(artifact['witnesses'])} witnesses"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
