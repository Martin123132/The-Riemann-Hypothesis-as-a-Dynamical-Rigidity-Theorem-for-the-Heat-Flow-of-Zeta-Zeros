#!/usr/bin/env python3
"""Prove a cofinal ordinary-theta tail theorem with a square-root count."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
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


STEM = "jensen_window_pf_newman_theta_forward_adaptive_sqrt_tail_gate"
DEFAULT_OUT = (
    REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
RESULT_ROOT = REPO_ROOT / "work" / "rh_compute" / "results"
SCRIPT_ROOT = REPO_ROOT / "work" / "rh_compute" / "scripts"
FINITE_STEM = "jensen_window_pf_newman_theta_forward_six_term_bridge_tail_gate"
FINITE_RESULT = RESULT_ROOT / f"{FINITE_STEM}.json"
FINITE_BUILDER = SCRIPT_ROOT / f"{FINITE_STEM}.py"
FINITE_CHECKER = SCRIPT_ROOT / f"check_{FINITE_STEM}.py"
MODULAR_RESULT = (
    RESULT_ROOT / "jensen_window_pf_newman_theta_adaptive_gamma_scale_tail_gate.json"
)
DATE = "2026-07-25"
PRECISION_BITS = 256
X_THRESHOLD = 245
FIRST_OMITTED_THRESHOLD = 7
RATIO_TARGET = Fraction(1, 50)


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


def serialize_positive(value: arb) -> dict:
    if value.lower() <= 0:
        raise RuntimeError(f"expected a positive enclosure: {value}")
    return serialize_ball(value)


def source_audit() -> dict:
    finite = json.loads(FINITE_RESULT.read_text(encoding="utf-8"))
    if finite.get("kind") != FINITE_STEM:
        raise RuntimeError("finite direct-tail source kind drifted")
    majorant = finite.get("tail_majorant", {})
    expected_formula = (
        "E_p(K)<=5*p!*exp(1/80)*pi*K^2*exp(-pi*K^2)"
        "/[4*theta_K*(1-rho_K)]; "
        "theta_K=1-121/(80*pi*K^2), "
        "rho_K=exp(2/K-pi*(2K+1)), p in {0,1}"
    )
    if majorant.get("formula") != expected_formula:
        raise RuntimeError("finite direct-tail formula drifted")
    if majorant.get("uniformity") != "0<=t<=1/5":
        raise RuntimeError("finite direct-tail time cap drifted")
    modular = json.loads(MODULAR_RESULT.read_text(encoding="utf-8"))
    if modular.get("kind") != (
        "jensen_window_pf_newman_theta_adaptive_gamma_scale_tail_gate"
    ):
        raise RuntimeError("modular comparison source kind drifted")
    return {
        "finite_result_sha256": file_hash(FINITE_RESULT),
        "finite_builder_sha256": file_hash(FINITE_BUILDER),
        "finite_checker_sha256": file_hash(FINITE_CHECKER),
        "modular_comparison_result_sha256": file_hash(MODULAR_RESULT),
        "ordinary_tail_formula": expected_formula,
        "ordinary_tail_uniformity": majorant["uniformity"],
    }


def geometry_diagnostics() -> dict:
    x = arb(X_THRESHOLD)
    pi = arb.pi()
    g = x / 8 + 2 * (1 + x).log()
    root_g = g.sqrt()
    k = arb(FIRST_OMITTED_THRESHOLD)
    theta = 1 - arb(121) / (80 * pi * k**2)
    rho = (arb(2) / k - pi * (2 * k + 1)).exp()
    denominator = theta * (1 - rho)
    constant = 20 * (arb(1) / 80).exp() * pi / denominator

    g_prime = arb(1) / 8 + 2 / (1 + x)
    log_slope_margin = (
        2 * (1 + x).log() + 1 - 2 * x / (1 + x)
    )
    common_decay_slope = 2 * pi / (1 + x) - 5 / x
    tuning_rate_margin = pi - arb(1) / 42 - 3
    ceiling_factor = 2 * g + 2
    value_envelope = (
        constant * x**4 * ceiling_factor / (1 + x) ** (2 * pi)
    )
    derivative_envelope = (
        constant
        * x**3
        * (x + 4)
        * ceiling_factor
        / (1 + x) ** (2 * pi)
    )
    saddle_ratio = (pi / 2).sqrt()

    if not (g.lower() > 36 and g.upper() < 49):
        raise RuntimeError("x=245 adaptive count is not K=7")
    for label, value in (
        ("theta", theta),
        ("denominator", denominator),
        ("log-slope margin", log_slope_margin),
        ("common decay slope", common_decay_slope),
        ("tuning rate margin", tuning_rate_margin),
    ):
        if value.lower() <= 0:
            raise RuntimeError(f"{label} is not positive")
    target = arb(RATIO_TARGET.numerator) / RATIO_TARGET.denominator
    if max(value_envelope.upper(), derivative_envelope.upper()) >= target:
        raise RuntimeError("square-root tail envelope misses 1/32")

    return {
        "g_at_x245": serialize_positive(g),
        "sqrt_g_at_x245": serialize_positive(root_g),
        "K_at_x245": FIRST_OMITTED_THRESHOLD,
        "theta_K7": serialize_positive(theta),
        "rho_K7": serialize_positive(rho),
        "denominator_K7": serialize_positive(denominator),
        "coefficient_C": serialize_positive(constant),
        "log_slope_margin_at_x245": serialize_positive(log_slope_margin),
        "ceiling_factor_at_x245": serialize_positive(ceiling_factor),
        "envelope_decay_slope_at_x245": serialize_positive(
            common_decay_slope
        ),
        "tuning_rate_margin": serialize_positive(tuning_rate_margin),
        "value_ratio_envelope_at_x245": serialize_positive(value_envelope),
        "derivative_ratio_envelope_at_x245": serialize_positive(
            derivative_envelope
        ),
        "ratio_target_exact": "1/50",
        "count_to_saddle_limit": serialize_positive(saddle_ratio),
    }


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    audit = source_audit()
    diagnostics = geometry_diagnostics()
    count = (
        "K_h(x)=ceil(sqrt(x/8+2*log(1+x)+h)), "
        "N_(F,h)(x)=K_h(x)-1, h>=0"
    )
    rows = [
        GateRow(
            id="ntfastg_01_ordinary_positive_partition",
            role="exact_identity",
            readiness="proved",
            claim=(
                "The ordinary positive-half-line theta series may be "
                "truncated directly without the modular switch."
            ),
            formula=(
                "Phi(u)=sum_(n>=1)phi_n(u), phi_n(u)>0 for u>=0"
            ),
            proof_boundary=(
                "Individual ordinary summands are not even; the proof uses "
                "a direct L1 tail and no endpoint integration by parts."
            ),
        ),
        GateRow(
            id="ntfastg_02_all_k_gaussian_tail",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "For every first omitted K>=7, the value and first-moment "
                "ordinary tails obey the stored explicit Gaussian bound."
            ),
            formula=audit["ordinary_tail_formula"],
            proof_boundary="Uniform for 0<=t<=1/5 and p in {0,1}.",
        ),
        GateRow(
            id="ntfastg_03_adaptive_sqrt_count",
            role="exact_definition",
            readiness="proved",
            claim=(
                "A tunable logarithmically enlarged square-root count "
                "supplies the cofinal ordinary truncation."
            ),
            formula=count,
            proof_boundary=(
                "For x>=245 and h>=0, K_h>=7 and "
                "K_h^2>=x/8+2*log(1+x)+h."
            ),
        ),
        GateRow(
            id="ntfastg_04_gamma_conversion",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "The adaptive Gaussian exponent leaves a power stronger "
                "than every direct C1 polynomial, together with a tunable "
                "exponential factor."
            ),
            formula=(
                "exp(-pi*K_h^2+pi*x/8)"
                "<=exp(-pi*h)*(1+x)^(-2*pi)"
            ),
            proof_boundary="This is immediate from the defining ceiling.",
        ),
        GateRow(
            id="ntfastg_05_ceiling_geometry",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "The ceiling loss is polynomial in h and is absorbed by "
                "the tunable exponential."
            ),
            formula=(
                "K_h^2<=2(g+h)+2; "
                "[1+h/(g+1)]exp(-pi*h)<=exp(-3h), "
                "g=x/8+2log(1+x)"
            ),
            proof_boundary=(
                "Uses g+1>42 and pi-1/42>3 on x>=245."
            ),
            diagnostics={
                "ceiling_factor": diagnostics["ceiling_factor_at_x245"],
                "tuning_rate_margin": diagnostics["tuning_rate_margin"],
            },
        ),
        GateRow(
            id="ntfastg_06_envelope_monotonicity",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "Both gamma-normalized C1 envelopes decrease for x>=245."
            ),
            formula=(
                "V(x)=C*x^4*(2g+2)/(1+x)^(2pi); "
                "D(x)=C*x^3*(x+4)*(2g+2)/(1+x)^(2pi); "
                "g'/(g+1)<1/x and both log derivatives are "
                "<5/x-2pi/(1+x)<0"
            ),
            proof_boundary=(
                "Uses 2log(1+x)+1-2x/(1+x)>0 and x>=245."
            ),
        ),
        GateRow(
            id="ntfastg_07_cofinal_sqrt_tail_theorem",
            role="exact_theorem",
            readiness="proved",
            claim=(
                "The tunable square-root ordinary truncation leaves less "
                "than exp(-3h)/50 of the gamma scale in both first-jet "
                "coordinates for every x>=245 and h>=0."
            ),
            formula=(
                "|J-J_(N_(F,h))^F|<exp(-3h-pi*x/8)/50; "
                "|J'-(J_(N_(F,h))^F)'|<exp(-3h-pi*x/8)/50"
            ),
            proof_boundary=(
                "This is an omitted-tail theorem only; it does not lower "
                "bound the retained ordinary first jet."
            ),
            diagnostics={
                "value_endpoint": diagnostics[
                    "value_ratio_envelope_at_x245"
                ],
                "derivative_endpoint": diagnostics[
                    "derivative_ratio_envelope_at_x245"
                ],
            },
        ),
        GateRow(
            id="ntfastg_08_retained_sqrt_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim=(
                "Use the saddle-scale retained ordinary sum as the preferred "
                "cofinal finite main, or transfer it to the corrected "
                "Riemann-Siegel main before proving first-jet separation."
            ),
            formula=(
                "For fixed h, N_(F,h)(x)/sqrt(x/(4*pi))->sqrt(pi/2); "
                "for any eta>0 choose "
                "h>=max(0,log(1/(50eta))/3)"
            ),
            proof_boundary=(
                "The margin dial removes any need for a t-independent "
                "approximation error near t=0. It does not supply the "
                "matching retained separation, strict Laguerre theorem, "
                "Lambda<=0, RH, or Clay-prize proof."
            ),
        ),
    ]
    return {
        "kind": STEM,
        "date": DATE,
        "status": (
            "cofinal ordinary-theta gamma-scale omitted-tail theorem with "
            "a logarithmically enlarged square-root retained count"
        ),
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "t_cap_exact": "1/5",
            "x_threshold": X_THRESHOLD,
            "adaptive_first_omitted": (
                "ceil(sqrt(x/8+2*log(1+x)+h))"
            ),
            "adaptive_retained": (
                "ceil(sqrt(x/8+2*log(1+x)+h))-1"
            ),
            "tuning_parameter": "h>=0",
            "first_omitted_threshold": FIRST_OMITTED_THRESHOLD,
            "ratio_target_exact": "exp(-3h)/50",
        },
        "source_audit": audit,
        "diagnostics": diagnostics,
        "theorem": {
            "domain": "0<=t<=1/5, x>=245",
            "count": count,
            "retained_transform": (
                "S_(N_F,t)^F=sum_(n=1)^N_F integral_0^infinity "
                "exp(tu^2)phi_n(u)cos(xu)du"
            ),
            "value_error": (
                "|J_t-J_(N_(F,h),t)^F|<exp(-3h-pi*x/8)/50"
            ),
            "first_derivative_error": (
                "|J_t'-(J_(N_(F,h),t)^F)'|"
                "<exp(-3h-pi*x/8)/50"
            ),
            "saddle_count_comparison": (
                "For fixed h, N_(F,h)(x)/sqrt(x/(4*pi))->sqrt(pi/2)"
            ),
            "arbitrary_tolerance_corollary": (
                "For eta>0, h>=max(0,log(1/(50eta))/3) gives both "
                "errors <eta*exp(-pi*x/8)"
            ),
        },
        "rows": [asdict(row) for row in rows],
        "proof_boundary": (
            "This artifact proves only cofinal omitted-tail upper bounds for "
            "an adaptive ordinary forward-theta truncation on "
            "0<=t<=1/5 and x>=245. It replaces the x^(3/4) modular count on "
            "the error side by a logarithmically enlarged square-root count, "
            "but it does not prove the retained first-jet disjunction. The "
            "tunable h removes the need for a t-independent approximation "
            "floor near t=0. No retained separation or strict Laguerre "
            "positivity, Lambda<=0, RH, PF-infinity, or Clay-prize result "
            "follows."
        ),
    }


def render_note(artifact: dict) -> str:
    d = artifact["diagnostics"]
    return "\n".join(
        [
            "# Newman Theta Forward Adaptive Square-Root Tail Gate",
            "",
            f"Date: {DATE}",
            "",
            "Status: exact cofinal omitted-tail theorem. This is not a",
            "retained-separation theorem and not a proof of `Lambda<=0` or RH.",
            "",
            "## Adaptive Ordinary Count",
            "",
            "For `x>=245`, set",
            "",
            "```text",
            "g(x)=x/8+2*log(1+x),",
            "K_h(x)=ceil(sqrt(g(x)+h)), h>=0,",
            "N_(F,h)(x)=K_h(x)-1.",
            "```",
            "",
            "The ordinary positive-half-line theta series is retained through",
            "`n=N_F(x)`. Its summands are not individually even, so this proof",
            "uses the direct `L1` and first-moment tail; it does not integrate",
            "individual summands by parts.",
            "",
            "## Gaussian Tail",
            "",
            "For `p=0,1` and `K>=7`,",
            "",
            "```text",
            "E_p(K)<=5*p!*exp(1/80)*pi*K^2*exp(-pi*K^2)",
            "         /[4*theta_K*(1-rho_K)],",
            "theta_K=1-121/(80*pi*K^2),",
            "rho_K=exp(2/K-pi*(2K+1)).",
            "```",
            "",
            "The defining ceiling gives",
            "",
            "```text",
            "exp(-pi*K_h^2+pi*x/8)",
            "  <=exp(-pi*h)*(1+x)^(-2*pi).",
            "```",
            "",
            "The ceiling bound `K_h^2<=2(g+h)+2`, together with",
            "`g+1>42` and `pi-1/42>3`, gives",
            "",
            "```text",
            "[1+h/(g+1)]*exp(-pi*h)<=exp(-3h).",
            "```",
            "",
            "## Cofinal Envelopes",
            "",
            "With",
            "",
            "```text",
            "C=20*exp(1/80)*pi/[theta_7*(1-rho_7)],",
            "V(x)=C*x^4*(2g+2)/(1+x)^(2*pi),",
            "D(x)=C*x^3*(x+4)*(2g+2)/(1+x)^(2*pi),",
            "```",
            "",
            "both envelopes decrease on `x>=245`. At the left endpoint,",
            "",
            "| envelope | upper ratio to `exp(-pi*x/8)` |",
            "|---|---:|",
            (
                "| value | "
                f"`{d['value_ratio_envelope_at_x245']['enclosure']}` |"
            ),
            (
                "| derivative | "
                f"`{d['derivative_ratio_envelope_at_x245']['enclosure']}` |"
            ),
            "",
            "Both are strictly below `1/50`.",
            "",
            "## Theorem",
            "",
            "For every `0<=t<=1/5`, `x>=245`, and `h>=0`,",
            "",
            "```text",
            "|J_t-J_(N_(F,h),t)^F|<exp(-3h-pi*x/8)/50,",
            "|J_t'-(J_(N_(F,h),t)^F)'|<exp(-3h-pi*x/8)/50.",
            "```",
            "",
            "The count is on the saddle scale:",
            "",
            "```text",
            "For fixed h, N_(F,h)(x)/sqrt(x/(4*pi))->sqrt(pi/2).",
            "```",
            "",
            "This replaces the `x^(3/4)` modular count on the omitted-error",
            "side with only about `sqrt(pi/2)` times the Riemann-Siegel count.",
            "For any `eta>0`, taking",
            "`h>=max(0,log(1/(50*eta))/3)` makes both errors smaller than",
            "`eta*exp(-pi*x/8)`. Thus the approximation can follow a margin",
            "that vanishes as `t->0`; retained first-jet separation remains open.",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "wrote Newman theta forward adaptive square-root tail gate: "
        f"{len(artifact['rows'])} rows, x>={X_THRESHOLD}, "
        "value/derivative ratio target exp(-3h)/50"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
