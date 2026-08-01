#!/usr/bin/env python3
"""Transfer the adaptive ordinary-theta C1 tail to the corrected RS scale."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR_DIR = REPO_ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR_DIR) not in sys.path:
    sys.path.insert(0, str(VENDOR_DIR))

import flint
from flint import arb
import sympy as sp


STEM = (
    "jensen_window_pf_newman_theta_forward_"
    "sqrt_to_corrected_rs_C1_transfer_gate"
)
DEFAULT_OUT = (
    REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
RESULT_ROOT = REPO_ROOT / "work" / "rh_compute" / "results"
SCRIPT_ROOT = REPO_ROOT / "work" / "rh_compute" / "scripts"
ADAPTIVE_STEM = (
    "jensen_window_pf_newman_theta_forward_adaptive_sqrt_tail_gate"
)
GLOBAL_STEM = (
    "jensen_window_pf_newman_polymath15_"
    "critical_C1_global_remainder_certificate"
)
CELL_STEM = (
    "jensen_window_pf_newman_polymath15_"
    "critical_C1_cell_remainder_certificate"
)
NORMALIZER_STEM = (
    "jensen_window_pf_newman_polymath15_normalized_laguerre_bridge"
)
DATE = "2026-07-25"
PRECISION_BITS = 256
L_MIN = 50
X_COARSE = 245


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def source_paths(stem: str) -> tuple[Path, Path, Path]:
    return (
        RESULT_ROOT / f"{stem}.json",
        SCRIPT_ROOT / f"{stem}.py",
        SCRIPT_ROOT / f"check_{stem}.py",
    )


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
    adaptive_paths = source_paths(ADAPTIVE_STEM)
    global_paths = source_paths(GLOBAL_STEM)
    cell_paths = source_paths(CELL_STEM)
    normalizer_paths = source_paths(NORMALIZER_STEM)

    adaptive = json.loads(adaptive_paths[0].read_text(encoding="utf-8"))
    if adaptive.get("kind") != ADAPTIVE_STEM:
        raise RuntimeError("adaptive ordinary-tail source kind drifted")
    theorem = adaptive.get("theorem", {})
    expected_value = (
        "|J_t-J_(N_(F,h),t)^F|<exp(-3h-pi*x/8)/50"
    )
    expected_first = (
        "|J_t'-(J_(N_(F,h),t)^F)'|"
        "<exp(-3h-pi*x/8)/50"
    )
    if theorem.get("value_error") != expected_value:
        raise RuntimeError("adaptive value-tail theorem drifted")
    if theorem.get("first_derivative_error") != expected_first:
        raise RuntimeError("adaptive derivative-tail theorem drifted")

    global_gate = json.loads(
        global_paths[0].read_text(encoding="utf-8")
    )
    if global_gate.get("kind") != GLOBAL_STEM:
        raise RuntimeError("global corrected C1 source kind drifted")
    expected_global = (
        "|r|<2500*exp(-3L/4), |r'|<5000*L*exp(-3L/4), "
        "r^2+(r'/L)^2<32000000*exp(-3L/2)"
    )
    if global_gate.get("exact", {}).get("global_c1") != expected_global:
        raise RuntimeError("global corrected C1 theorem drifted")

    cell = json.loads(cell_paths[0].read_text(encoding="utf-8"))
    if cell.get("kind") != CELL_STEM:
        raise RuntimeError("fixed-cell corrected C1 source kind drifted")
    expected_split = (
        "Z_t=J_hat_(N,t)+r, where J_hat agrees on the real axis "
        "with the corrected main J=P-Q"
    )
    if cell.get("exact", {}).get("normalized_split") != expected_split:
        raise RuntimeError("corrected normalized split drifted")
    expected_eab = "e_A+e_B<1000*exp(-3L/4)"
    if cell.get("exact", {}).get("eAB") != expected_eab:
        raise RuntimeError("bulk A+B remainder envelope drifted")

    normalizer = json.loads(
        normalizer_paths[0].read_text(encoding="utf-8")
    )
    if normalizer.get("kind") != NORMALIZER_STEM:
        raise RuntimeError("normalizer source kind drifted")
    published = normalizer.get("exact", {}).get("published_input", {})
    expected_mt = "M_t(s)=exp((t/4)*alpha(s)^2)*M_0(s)"
    expected_alpha = (
        "alpha(s)=1/(2s)+1/(s-1)+(1/2)Log(s/(2*pi))"
    )
    if published.get("M_t") != expected_mt:
        raise RuntimeError("published normalizer formula drifted")
    if published.get("alpha") != expected_alpha:
        raise RuntimeError("published alpha formula drifted")

    audit = {
        "adaptive_value_error": expected_value,
        "adaptive_first_derivative_error": expected_first,
        "global_corrected_C1": expected_global,
        "corrected_normalized_split": expected_split,
        "bulk_A_plus_B_envelope": expected_eab,
        "published_M_t": expected_mt,
        "published_alpha": expected_alpha,
    }
    for label, paths in (
        ("adaptive", adaptive_paths),
        ("global", global_paths),
        ("cell", cell_paths),
        ("normalizer", normalizer_paths),
    ):
        audit[f"{label}_result_sha256"] = digest(paths[0])
        audit[f"{label}_builder_sha256"] = digest(paths[1])
        audit[f"{label}_checker_sha256"] = digest(paths[2])
    return audit


def symbolic_audit() -> None:
    x = sp.symbols("x", positive=True)
    delta_j = sp.Function("delta_j")(x)
    delta_h = delta_j / (16 * x**4)
    expected_delta_h_prime = (
        sp.diff(delta_j, x) / (16 * x**4)
        - delta_j / (4 * x**5)
    )
    if sp.simplify(sp.diff(delta_h, x) - expected_delta_h_prime) != 0:
        raise RuntimeError("characteristic-to-heat derivative identity failed")

    a = sp.Function("a")(x)
    e = sp.Function("e")(x)
    normalized = e / a
    expected = sp.diff(e, x) / a - sp.diff(a, x) * e / a**2
    if sp.simplify(sp.diff(normalized, x) - expected) != 0:
        raise RuntimeError("normalized derivative identity failed")

    ordinary, corrected, r_rs, e_theta = sp.symbols(
        "ordinary corrected r_rs e_theta", real=True
    )
    if sp.expand(
        (ordinary + e_theta) - (corrected + r_rs)
        - (ordinary - corrected - r_rs + e_theta)
    ) != 0:
        raise RuntimeError("dual-approximation transfer identity failed")


def interval_diagnostics() -> dict:
    flint.ctx.prec = PRECISION_BITS
    pi = arb.pi()
    ell = arb(L_MIN)
    x = 4 * pi * ell.exp()
    coarse = arb(X_COARSE)

    prefactor_margin = (
        pi ** (arb(1) / 4)
        * (-arb(1) / (12 * coarse**2)).exp()
        - 1
    )
    coarse_ell = (coarse / (4 * pi)).log()
    alpha_square_floor_coarse = (
        coarse_ell / 2 - 1 / coarse**2
    ) ** 2 - (pi / 4) ** 2
    alpha_square_floor = (
        ell / 2 - 1 / x**2
    ) ** 2 - (pi / 4) ** 2
    log_a_ratio = (
        (ell / 4 + arb(1) / 2)
        * (1 + arb(1) / (2 * x))
        / (ell / 2)
    )
    scaled_derivative_factor = (
        arb(1) / 2 + (1 + 4 / x) / ell
    )
    x_power = ((arb(23) / 4) * x.log()).exp()
    theta_value_over_rs_envelope = (
        (3 * ell / 4).exp() / (25 * x_power)
    )
    theta_derivative_over_rs_envelope = (
        scaled_derivative_factor * theta_value_over_rs_envelope
    )
    composite_constant = (
        arb(2501) ** 2 + arb(5001) ** 2
    ).sqrt()
    ordinary_to_corrected_squared_target_ratio = (
        (arb(2) / 625)
        / arb(32_000_000)
        / (((arb(23) / 2) * x.log()).exp())
        * (3 * ell / 2).exp()
    )

    guards = {
        "normalizer_prefactor_margin_at_x245": prefactor_margin,
        "Re_alpha_squared_floor_at_x245": alpha_square_floor_coarse,
        "Re_alpha_squared_floor_at_L50": alpha_square_floor,
        "log_A_derivative_bound_ratio_at_L50": log_a_ratio,
        "scaled_theta_derivative_factor_at_L50": (
            scaled_derivative_factor
        ),
        "theta_value_envelope_over_exp_minus_3L_over_4_at_L50": (
            theta_value_over_rs_envelope
        ),
        "theta_scaled_derivative_envelope_over_exp_minus_3L_over_4_at_L50": (
            theta_derivative_over_rs_envelope
        ),
        "composite_transfer_vector_constant": composite_constant,
        "ordinary_to_corrected_squared_target_ratio_at_L50_h0": (
            ordinary_to_corrected_squared_target_ratio
        ),
    }
    if prefactor_margin.lower() <= 0:
        raise RuntimeError("normalizer lower-bound prefactor failed")
    if alpha_square_floor_coarse.lower() <= 0:
        raise RuntimeError("coarse Re(alpha^2) lower bound failed")
    if alpha_square_floor.lower() <= 0:
        raise RuntimeError("Re(alpha^2) lower bound failed")
    if log_a_ratio.upper() >= 1:
        raise RuntimeError("log-normalizer derivative bound failed")
    if scaled_derivative_factor.upper() >= arb(53) / 100:
        raise RuntimeError("scaled theta derivative factor failed")
    if theta_value_over_rs_envelope.upper() >= 1:
        raise RuntimeError("theta value envelope is not below RS scale")
    if theta_derivative_over_rs_envelope.upper() >= 1:
        raise RuntimeError("theta derivative envelope is not below RS scale")
    if composite_constant.upper() >= 5600:
        raise RuntimeError("composite transfer constant failed")
    if ordinary_to_corrected_squared_target_ratio.upper() >= 1:
        raise RuntimeError("target-scale comparison failed")
    return {
        "precision_bits": PRECISION_BITS,
        "L_min": L_MIN,
        "x_at_L50": serialize_positive(x),
        **{
            key: serialize_positive(value)
            for key, value in guards.items()
        },
    }


def build_artifact() -> dict:
    symbolic_audit()
    audit = source_audit()
    diagnostics = interval_diagnostics()
    domain = (
        "L=log(x/(4*pi))>=50, 0<t<=1/5, t*L<=25, h>=0"
    )
    ordinary = (
        "O_(h,t)=J_(N_(F,h),t)^F/(16*x^4*A_t)"
    )
    theta_value = "E_F0=exp(-3h)/(25*x^(23/4))"
    theta_scaled_first = (
        "E_F1=E_F0*(1/2+(1+4/x)/L)"
    )
    rows = [
        GateRow(
            id="ntfsrsg_01_exact_normalizer_amplitude",
            role="exact_identity",
            readiness="proved",
            claim=(
                "The published real-axis normalizer has an elementary "
                "closed amplitude relative to the gamma scale."
            ),
            formula=(
                "A_0/exp(-pi*x/8)="
                "pi^(1/4)*(1+x^2)^(7/8)"
                "*exp((x*atan(1/x)-1)/4)/32"
            ),
            proof_boundary=(
                "Principal Log, s=(1-i*x)/2, and x>0; uses "
                "atan(x)+atan(1/x)=pi/2."
            ),
        ),
        GateRow(
            id="ntfsrsg_02_uniform_normalizer_lower_bound",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "Positive heat time does not reduce the normalizer "
                "amplitude on the transfer region."
            ),
            formula=(
                "A_t=A_0*exp(t*Re(alpha(s)^2)/4)>="
                "A_0>=exp(-pi*x/8)*x^(7/4)/32"
            ),
            proof_boundary=(
                "For x>=245, Re(alpha)^2>0; atan(y)>=y-y^3/3 "
                "for 0<=y<=1."
            ),
            diagnostics={
                "prefactor_margin": diagnostics[
                    "normalizer_prefactor_margin_at_x245"
                ],
                "alpha_square_floor": diagnostics[
                    "Re_alpha_squared_floor_at_x245"
                ],
            },
        ),
        GateRow(
            id="ntfsrsg_03_normalizer_first_derivative",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "The logarithmic normalizer derivative is harmless at "
                "the scale-adapted first-jet norm."
            ),
            formula=(
                "|(log A_t)'|<L/2; |alpha|<=L/2+1, "
                "|alpha'|<=2/x"
            ),
            proof_boundary=(
                "Uses t<=25/L, L>=50, and ds/dx=-i/2."
            ),
            diagnostics=diagnostics[
                "log_A_derivative_bound_ratio_at_L50"
            ],
        ),
        GateRow(
            id="ntfsrsg_04_normalized_ordinary_tail",
            role="exact_theorem",
            readiness="proved",
            claim=(
                "After conversion from the characteristic J_t=16*x^4*H_t, "
                "the adaptive ordinary main approximates Z_t=H_t/A_t "
                "with a power-23/4 first-jet error."
            ),
            formula=(
                f"{ordinary}; |Z_t-O_(h,t)|<{theta_value}; "
                f"|(Z_t-O_(h,t))'|/L<{theta_scaled_first}"
            ),
            proof_boundary=(
                "Combines the adaptive ordinary C1 tail with the exact "
                "normalizer lower and derivative bounds."
            ),
            diagnostics=diagnostics[
                "scaled_theta_derivative_factor_at_L50"
            ],
        ),
        GateRow(
            id="ntfsrsg_05_direct_ordinary_transversality",
            role="exact_sufficient_target",
            readiness="proved",
            claim=(
                "A tiny explicit lower bound for the retained ordinary "
                "first jet excludes a double zero directly."
            ),
            formula=(
                "T_L[O]=O^2+(O'/L)^2>"
                "2*exp(-6h)/(625*x^(23/2)) "
                "implies (H_t,H_t')!=(0,0)"
            ),
            proof_boundary=(
                "At a double zero O=-(Z-O) and O'=-(Z-O)'; "
                "the scaled derivative factor is <53/100<1."
            ),
        ),
        GateRow(
            id="ntfsrsg_06_dual_approximation_identity",
            role="exact_identity",
            readiness="proved",
            claim=(
                "The ordinary and corrected Riemann-Siegel mains are "
                "related through their independently certified remainders."
            ),
            formula=(
                "Z_t=O_(h,t)+e_F=J_hat_(N,t)+r_RS; "
                "O_(h,t)-J_hat_(N,t)=r_RS-e_F"
            ),
            proof_boundary=(
                "This is a relation between two approximations of Z_t, "
                "not a termwise identity between their finite sums."
            ),
        ),
        GateRow(
            id="ntfsrsg_07_corrected_rs_transfer_envelope",
            role="exact_theorem",
            readiness="proved",
            claim=(
                "The dual identity gives a global C1 comparison across "
                "all corrected cutoff transitions."
            ),
            formula=(
                "|O_(h,t)-J_hat_(N,t)|<2501*exp(-3L/4); "
                "|O_(h,t)'-J_hat_(N,t)'|/L<5001*exp(-3L/4); "
                "Euclidean scaled-jet distance<5600*exp(-3L/4)"
            ),
            proof_boundary=(
                "The added ordinary error is below exp(-3L/4); the "
                "constants deliberately round outward."
            ),
            diagnostics={
                "theta_value_ratio": diagnostics[
                    "theta_value_envelope_over_exp_minus_3L_over_4_at_L50"
                ],
                "theta_derivative_ratio": diagnostics[
                    "theta_scaled_derivative_envelope_over_exp_minus_3L_over_4_at_L50"
                ],
                "vector_constant": diagnostics[
                    "composite_transfer_vector_constant"
                ],
            },
        ),
        GateRow(
            id="ntfsrsg_08_scale_mismatch",
            role="rigorous_architecture_guard",
            readiness="proved",
            claim=(
                "The currently certified corrected-RS envelope does not "
                "preserve the much smaller ordinary-tail target."
            ),
            formula=(
                "ordinary amplitude envelope=O(exp(-3h-23L/4)); "
                "corrected-RS amplitude envelope=O(exp(-3L/4)); "
                "their certified-envelope ratio costs O(exp(5L+3h)); "
                "already e_A+e_B<1000*exp(-3L/4)"
            ),
            proof_boundary=(
                "This compares proved upper envelopes, not the unknown "
                "actual errors. Transfer is exact, but the two sufficient "
                "separation thresholds are not quantitatively equivalent. "
                "Higher endpoint C_k terms alone cannot improve the stored "
                "bulk A+B envelope."
            ),
            diagnostics=diagnostics[
                "ordinary_to_corrected_squared_target_ratio_at_L50_h0"
            ],
        ),
        GateRow(
            id="ntfsrsg_09_retained_arithmetic_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim=(
                "Prove a noncircular first-jet lower bound either for the "
                "adaptive ordinary main at its tiny direct threshold or "
                "for the corrected Riemann-Siegel main at its published "
                "remainder threshold."
            ),
            formula=(
                "ordinary target: T_L[O_h]>"
                "2*exp(-6h)/(625*x^(23/2)); "
                "corrected target: T_L[J_hat]>"
                "32000000*exp(-3L/2)"
            ),
            proof_boundary=(
                "Neither lower bound is proved here. No retained "
                "separation, endpoint simplicity, Lambda<=0, RH, "
                "PF-infinity, or Clay-prize proof follows."
            ),
        ),
    ]
    return {
        "kind": STEM,
        "date": DATE,
        "status": (
            "exact C1 transfer between the adaptive ordinary-theta main "
            "and the corrected Polymath-15 Riemann-Siegel main, with a "
            "quantified certified-envelope mismatch"
        ),
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "L_min": L_MIN,
            "t_cap_exact": "1/5",
            "scaled_time_cap": "t*L<=25",
            "tuning_parameter": "h>=0",
            "domain": domain,
        },
        "source_audit": audit,
        "diagnostics": diagnostics,
        "theorem": {
            "domain": domain,
            "ordinary_normalized_main": ordinary,
            "ordinary_value_error": theta_value,
            "ordinary_scaled_first_derivative_error": theta_scaled_first,
            "ordinary_direct_sufficient_target": (
                "T_L[O_h]>2*exp(-6h)/(625*x^(23/2))"
            ),
            "dual_approximation_identity": (
                "Z_t=O_(h,t)+e_F=J_hat_(N,t)+r_RS"
            ),
            "corrected_value_transfer": (
                "|O_(h,t)-J_hat_(N,t)|<2501*exp(-3L/4)"
            ),
            "corrected_scaled_derivative_transfer": (
                "|O_(h,t)'-J_hat_(N,t)'|/L"
                "<5001*exp(-3L/4)"
            ),
            "corrected_vector_transfer": (
                "||(O-J_hat),(O'-J_hat')/L||_2"
                "<5600*exp(-3L/4)"
            ),
        },
        "rows": [asdict(row) for row in rows],
        "proof_boundary": (
            "This gate proves an exact first-jet transfer only on the "
            "critical overlap L>=50, 0<t<=1/5, tL<=25. It shows that the "
            "adaptive square-root ordinary truncation has a tunable "
            "power-23/4 normalized tail and gives a direct sufficient "
            "retained first-jet target. It also proves that the currently "
            "certified corrected-RS remainder envelope is exponentially "
            "coarser, so 'equivalent after transfer' cannot mean equal "
            "quantitative separation thresholds. Higher endpoint C_k terms alone "
            "leave the stored bulk A+B envelope at the corrected scale. It "
            "does not prove either "
            "retained lower bound, endpoint simplicity at t=0, Lambda<=0, "
            "RH, PF-infinity, or a Clay-prize result."
        ),
    }


def render_note(artifact: dict) -> str:
    d = artifact["diagnostics"]
    return "\n".join(
        [
            "# Theta Square-Root to Corrected Riemann-Siegel C1 Transfer",
            "",
            f"Date: {DATE}",
            "",
            "Status: exact transfer and quantitative architecture guard.",
            "This is not retained separation and not a proof of `Lambda<=0`",
            "or RH.",
            "",
            "## Common Normalization",
            "",
            "On",
            "",
            "```text",
            "L=log(x/(4*pi))>=50, 0<t<=1/5, tL<=25, h>=0,",
            "Z_t=H_t/A_t,",
            "O_(h,t)=J_(N_(F,h),t)^F/(16*x^4*A_t),",
            "```",
            "",
            "the exact real-axis normalizer satisfies",
            "",
            "```text",
            "A_0/exp(-pi*x/8)",
            " =pi^(1/4)*(1+x^2)^(7/8)",
            "  *exp((x*atan(1/x)-1)/4)/32.",
            "```",
            "",
            "Elementary bounds on `alpha` give",
            "",
            "```text",
            "A_t>=A_0>=exp(-pi*x/8)*x^(7/4)/32,",
            "|(log A_t)'|<L/2.",
            "```",
            "",
            "The certified prefactor margin at the coarse threshold is",
            f"`{d['normalizer_prefactor_margin_at_x245']['enclosure']}`.",
            "",
            "## Normalized Ordinary Tail",
            "",
            "Converting `J_t=16*x^4*H_t` and differentiating exactly gives",
            "",
            "```text",
            "|Z_t-O_(h,t)|<E_F0,",
            "|(Z_t-O_(h,t))'|/L<E_F1,",
            "E_F0=exp(-3h)/(25*x^(23/4)),",
            "E_F1=E_F0*(1/2+(1+4/x)/L).",
            "```",
            "",
            "The final factor is below `53/100`. Therefore the direct",
            "ordinary sufficient target is",
            "",
            "```text",
            "T_L[O_h]=O_h^2+(O_h'/L)^2",
            "  >2*exp(-6h)/(625*x^(23/2)).",
            "```",
            "",
            "Any point satisfying this inequality cannot be a double zero",
            "of `H_t`.",
            "",
            "## Corrected-Main Transfer",
            "",
            "The two finite mains are not termwise identical. They satisfy",
            "",
            "```text",
            "Z_t=O_(h,t)+e_F=J_hat_(N,t)+r_RS,",
            "O_(h,t)-J_hat_(N,t)=r_RS-e_F.",
            "```",
            "",
            "Combining the two independently certified remainders yields",
            "",
            "```text",
            "|O_(h,t)-J_hat_(N,t)|<2501*exp(-3L/4),",
            "|O_(h,t)'-J_hat_(N,t)'|/L<5001*exp(-3L/4),",
            "scaled first-jet distance<5600*exp(-3L/4).",
            "```",
            "",
            "## Architecture Consequence",
            "",
            "The transfer is exact, but the certified envelopes have very",
            "different scales:",
            "",
            "```text",
            "ordinary amplitude: O(exp(-3h-23L/4)),",
            "corrected-RS amplitude: O(exp(-3L/4)).",
            "```",
            "",
            "Thus the current corrected-RS envelope costs about `exp(5L+3h)`",
            "relative to the ordinary one. This is a comparison of proved",
            "upper bounds, not a claim that the actual RS error is that large.",
            "Moreover, the current bulk estimate already has",
            "",
            "```text",
            "e_A+e_B<1000*exp(-3L/4).",
            "```",
            "",
            "Therefore adding higher classical endpoint `C_k` corrections alone",
            "cannot close the certified five-exponent gap. Such a route would",
            "have to upgrade the bulk heat-flow saddle approximation and the",
            "endpoint expansion together.",
            "A proof must therefore establish either the tiny direct ordinary",
            "first-jet target or the existing corrected-main target",
            "`T_L[J_hat]>32000000*exp(-3L/2)`. The phrase 'equivalent after",
            "transfer' is valid algebraically but not as equality of current",
            "quantitative proof burdens.",
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
        "wrote theta square-root to corrected RS C1 transfer gate: "
        f"{len(artifact['rows'])} rows, L>={L_MIN}, "
        "exact dual approximation and quantified envelope mismatch"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
