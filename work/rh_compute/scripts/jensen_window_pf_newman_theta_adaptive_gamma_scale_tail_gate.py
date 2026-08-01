#!/usr/bin/env python3
"""Prove an explicit cofinal gamma-scale bound for the modular theta tail."""

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

from jensen_window_pf_newman_theta_switch_defect_explicit_constant_gate import (
    arb_power_rational,
    compose_budget,
    polynomial_rows,
    safe_phase_constant,
    serialize_positive,
    switch_derivative_rows,
)


STEM = "jensen_window_pf_newman_theta_adaptive_gamma_scale_tail_gate"
DEFAULT_OUT = (
    REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
RESULT_ROOT = REPO_ROOT / "work" / "rh_compute" / "results"
BASE_STEM = "jensen_window_pf_newman_theta_switch_defect_explicit_constant_gate"
BASE_RESULT = RESULT_ROOT / f"{BASE_STEM}.json"
BASE_BUILDER = (
    REPO_ROOT / "work" / "rh_compute" / "scripts" / f"{BASE_STEM}.py"
)
BASE_CHECKER = (
    REPO_ROOT / "work" / "rh_compute" / "scripts" / f"check_{BASE_STEM}.py"
)
CONTRACT_RESULT = (
    RESULT_ROOT
    / "jensen_window_pf_newman_theta_adaptive_modular_c1_remainder_contract.json"
)
DATE = "2026-07-25"
PRECISION_BITS = 256
N_THRESHOLD = 63
X_THRESHOLD = 245
BOUNDARY_N = (62, 63, 64)
RATIO_TARGET = Fraction(1, 4)
MAX_FORWARD_DEGREE = 20
MAX_COMPACT_DEGREE = 22
MAX_STRETCHED_DEGREE = 22


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
    base = json.loads(BASE_RESULT.read_text(encoding="utf-8"))
    if base.get("kind") != BASE_STEM:
        raise RuntimeError("base explicit-constant artifact kind drifted")
    parameters = base.get("parameters", {})
    if (
        parameters.get("m_order") != 9
        or parameters.get("t_cap_exact") != "1/5"
        or parameters.get("retained_count_upper_bound") is not None
    ):
        raise RuntimeError("base explicit-constant parameters drifted")
    contract = json.loads(CONTRACT_RESULT.read_text(encoding="utf-8"))
    contract_text = json.dumps(contract, sort_keys=True)
    for marker in (
        "x^(4-m)*d_(N,0,m)",
        "4*x^(3-m)*d_(N,0,m)+x^(4-m)*d_(N,1,m)",
    ):
        if marker not in contract_text:
            raise RuntimeError(f"adaptive error-contract marker missing: {marker}")
    return {
        "base_result_sha256": file_hash(BASE_RESULT),
        "base_builder_sha256": file_hash(BASE_BUILDER),
        "base_checker_sha256": file_hash(BASE_CHECKER),
        "contract_result_sha256": file_hash(CONTRACT_RESULT),
        "base_kind": base["kind"],
        "m_order": 9,
        "t_cap_exact": "1/5",
        "retained_count_upper_bound": None,
    }


def endpoint_x(retained: int) -> arb:
    return arb_power_rational(arb(retained), Fraction(4, 3)) - 1


def boundary_row(
    retained: int,
    polynomials: list[dict],
    switches: list[dict],
) -> dict:
    budget = compose_budget(retained + 1, polynomials, switches)
    d0 = arb(budget["d0_m9_t1_5"]["total"]["enclosure"])
    d1 = arb(budget["d1_m9_t1_5"]["total"]["enclosure"])
    x_right = endpoint_x(retained)
    gamma_scale = (-arb.pi() * x_right / 8).exp()
    value_error = 16 * d0 / x_right**5
    derivative_error = 64 * d0 / x_right**6 + 16 * d1 / x_right**5
    return {
        "N": retained,
        "K": retained + 1,
        "cell": (
            "((N-1)^(4/3)-1,N^(4/3)-1], "
            "N(x)=ceil((1+x)^(3/4))"
        ),
        "x_right": serialize_positive(x_right),
        "d0": serialize_positive(d0),
        "d1": serialize_positive(d1),
        "gamma_scale": serialize_positive(gamma_scale),
        "J_error": serialize_positive(value_error),
        "J_prime_error": serialize_positive(derivative_error),
        "J_error_over_gamma": serialize_positive(value_error / gamma_scale),
        "J_prime_error_over_gamma": serialize_positive(
            derivative_error / gamma_scale
        ),
    }


def monotonicity_diagnostics() -> dict:
    pi = arb.pi()
    spectral_rate = pi / 8
    c_safe = safe_phase_constant()
    n = arb(N_THRESHOLD)
    n_next = arb(N_THRESHOLD + 1)
    delta = (
        arb_power_rational(n_next, Fraction(4, 3))
        - arb_power_rational(n, Fraction(4, 3))
    )
    forward_log_ratio = (
        arb(MAX_FORWARD_DEGREE) / (N_THRESHOLD + 1)
        - pi * (2 * (N_THRESHOLD + 1) + 1)
        + spectral_rate * delta
    )
    compact_log_ratio = (
        arb(MAX_COMPACT_DEGREE) / (N_THRESHOLD + 1)
        - pi / arb(2).sqrt() * (2 * (N_THRESHOLD + 1) + 1)
        + spectral_rate * delta
    )
    large_derivative_margin = (
        arb(4)
        * (c_safe - spectral_rate)
        / 3
        * arb_power_rational(n, Fraction(1, 3))
        - arb(MAX_STRETCHED_DEGREE) / n
    )
    value_cell_slope = spectral_rate - arb(5) / X_THRESHOLD
    derivative_cell_slope = spectral_rate - arb(6) / X_THRESHOLD
    if c_safe.lower() <= spectral_rate.upper():
        raise RuntimeError("safe reflected exponent does not beat gamma rate")
    if forward_log_ratio.upper() >= 0:
        raise RuntimeError("forward normalized ratio is not decreasing")
    if compact_log_ratio.upper() >= 0:
        raise RuntimeError("compact normalized ratio is not decreasing")
    if large_derivative_margin.lower() <= 0:
        raise RuntimeError("large reflected normalized tail is not decreasing")
    if min(value_cell_slope.lower(), derivative_cell_slope.lower()) <= 0:
        raise RuntimeError("within-cell normalized error is not increasing")
    return {
        "spectral_rate_pi_over_8": serialize_positive(spectral_rate),
        "reflected_safe_rate": serialize_positive(c_safe),
        "reflected_rate_margin": serialize_positive(c_safe - spectral_rate),
        "worst_forward_log_ratio_upper_at_N63": serialize_ball(
            forward_log_ratio
        ),
        "worst_compact_log_ratio_upper_at_N63": serialize_ball(
            compact_log_ratio
        ),
        "worst_large_derivative_margin_at_N63": serialize_positive(
            large_derivative_margin
        ),
        "value_within_cell_log_slope_at_x245": serialize_positive(
            value_cell_slope
        ),
        "derivative_within_cell_log_slope_at_x245": serialize_positive(
            derivative_cell_slope
        ),
        "large_tail_condition": (
            "(4/3)*(c_safe-pi/8)*N^(1/3)-A/N>0 "
            "for N>=63 and A<=22"
        ),
    }


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    audit = source_audit()
    polynomials = polynomial_rows()
    switches = switch_derivative_rows()
    boundaries = [
        boundary_row(retained, polynomials, switches)
        for retained in BOUNDARY_N
    ]
    by_n = {row["N"]: row for row in boundaries}
    passing = by_n[N_THRESHOLD]
    target = arb(RATIO_TARGET.numerator) / RATIO_TARGET.denominator
    for field in ("J_error_over_gamma", "J_prime_error_over_gamma"):
        if arb(passing[field]["enclosure"]).upper() >= target:
            raise RuntimeError(f"threshold ratio does not beat 1/4: {field}")
        if arb(by_n[N_THRESHOLD - 1][field]["enclosure"]).lower() <= 1:
            raise RuntimeError(f"N=62 failure guard drifted: {field}")
    diagnostics = monotonicity_diagnostics()

    rows = [
        GateRow(
            id="ntagstg_01_adaptive_cell_partition",
            role="exact_identity",
            readiness="proved",
            claim=(
                "The adaptive count N(x)=ceil((1+x)^(3/4)) partitions the "
                "positive axis into exact half-open cells."
            ),
            formula=(
                "C_N=((N-1)^(4/3)-1,N^(4/3)-1], "
                "N(x)=N on C_N"
            ),
            proof_boundary="Exact ceiling arithmetic.",
        ),
        GateRow(
            id="ntagstg_02_within_cell_monotonicity",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "For x>=245 and fixed N, each gamma-normalized m=9 error "
                "multiplier increases with x, so the cell maximum is its "
                "right endpoint."
            ),
            formula=(
                "d/dx log(exp(pi*x/8)/x^s)=pi/8-s/x>0, "
                "s in {5,6}"
            ),
            proof_boundary="The d0 and d1 budgets are constant on each cell.",
        ),
        GateRow(
            id="ntagstg_03_forward_endpoint_monotonicity",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "Every forward template decreases after gamma normalization "
                "from N=63 onward."
            ),
            formula=(
                "log(F_(N+1)/F_N)+pi*(x_(N+1)-x_N)/8"
                "<=D/(N+1)-pi*(2N+3)+pi*Delta_N/8<0, D<=20"
            ),
            proof_boundary=(
                "The denominator theta increases and the geometric ratio "
                "rho_F decreases with N."
            ),
        ),
        GateRow(
            id="ntagstg_04_compact_endpoint_monotonicity",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "Every compact-reflected template decreases after gamma "
                "normalization from N=63 onward."
            ),
            formula=(
                "log(C_(N+1)/C_N)+pi*(x_(N+1)-x_N)/8"
                "<=D/(N+1)-pi*(2N+3)/sqrt(2)+pi*Delta_N/8<0, D<=22"
            ),
            proof_boundary="The compact geometric ratio decreases with N.",
        ),
        GateRow(
            id="ntagstg_05_large_endpoint_monotonicity",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "Every large-reflected unimodal tail decreases after gamma "
                "normalization from N=63 onward."
            ),
            formula=(
                "G_A(y)=exp(pi*y^(4/3)/8)"
                "[f_A(y)+integral_y^infinity f_A]; "
                "G_A'<0 if (4/3)(c_safe-pi/8)y^(1/3)-A/y>0"
            ),
            proof_boundary=(
                "For f_A=y^A exp(-c_safe*y^(4/3)), "
                "integral_y^infinity f_A<=f_A/h_A and h_A is increasing."
            ),
        ),
        GateRow(
            id="ntagstg_06_positive_composition",
            role="exact_composition",
            readiness="proved",
            claim=(
                "Positive Leibniz and Hermite coefficients preserve the "
                "three template monotonicity statements for full d0 and d1."
            ),
            formula=(
                "d0_N,d1_N are fixed positive linear combinations of the "
                "forward, compact-reflected, and large-reflected templates"
            ),
            proof_boundary="No cancellation is used in the error theorem.",
        ),
        GateRow(
            id="ntagstg_07_cofinal_gamma_scale_theorem",
            role="exact_theorem",
            readiness="proved",
            claim=(
                "For every real x>=245, the adaptive modular tail contributes "
                "less than one quarter of the natural exp(-pi*x/8) scale to "
                "both the value and first derivative."
            ),
            formula=(
                "N(x)=ceil((1+x)^(3/4)); "
                "|J-J_N|<exp(-pi*x/8)/4; "
                "|J'-J_N'|<exp(-pi*x/8)/4"
            ),
            proof_boundary=(
                "This is an omitted-tail upper theorem, not a lower bound "
                "for the retained first jet."
            ),
            diagnostics=passing,
        ),
        GateRow(
            id="ntagstg_08_retained_margin_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim=(
                "The cofinal Newman task is reduced on x>=245 to proving "
                "that the retained adaptive first jet cannot enter the "
                "quarter-gamma error square."
            ),
            formula=(
                "max(|J_N|,|J_N'|)>exp(-pi*x/8)/4 "
                "is a sufficient retained-margin target"
            ),
            proof_boundary=(
                "No retained lower margin, transition-cell termination, "
                "Lambda<=0, RH, or Clay-prize proof is supplied."
            ),
        ),
    ]
    return {
        "kind": STEM,
        "date": DATE,
        "status": "explicit cofinal adaptive gamma-scale omitted-tail theorem",
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "m_order": 9,
            "t_cap_exact": "1/5",
            "adaptive_count": "ceil((1+x)^(3/4))",
            "x_threshold": X_THRESHOLD,
            "n_threshold": N_THRESHOLD,
            "ratio_target_exact": "1/4",
            "boundary_n": list(BOUNDARY_N),
        },
        "source_audit": audit,
        "cell_geometry": {
            "cell_N": "((N-1)^(4/3)-1,N^(4/3)-1]",
            "right_endpoint": "x_N=N^(4/3)-1",
            "first_complete_cell_above_x245": N_THRESHOLD,
        },
        "monotonicity": diagnostics,
        "boundary_witnesses": boundaries,
        "theorem": {
            "domain": "0<=t<=1/5, x>=245",
            "retained_count": "N(x)=ceil((1+x)^(3/4))",
            "value_error": "|J_t-J_(N,t)|<exp(-pi*x/8)/4",
            "first_derivative_error": (
                "|J_t'-J_(N,t)'|<exp(-pi*x/8)/4"
            ),
            "uniformity": (
                "The d0/d1 budgets use T=1/5 and hence are uniform in "
                "0<=t<=1/5."
            ),
        },
        "rows": [asdict(row) for row in rows],
        "proof_boundary": (
            "This artifact proves an explicit cofinal gamma-scale upper bound "
            "for the omitted modular tail with N(x)=ceil((1+x)^(3/4)) and "
            "x>=245. It does not lower-bound the retained J_N or J_N' first "
            "jet, certify the finite bridge 69<x<245, prove a terminating "
            "retained transition cover, establish strict Laguerre positivity, "
            "prove Lambda<=0 or RH, or supply a Clay-prize result."
        ),
    }


def render_note(artifact: dict) -> str:
    lines = [
        "# Newman Theta Adaptive Gamma-Scale Tail Gate",
        "",
        f"Date: {DATE}",
        "",
        "Status: explicit cofinal omitted-tail theorem. This is not a",
        "retained-separation theorem and not a proof of `Lambda<=0` or RH.",
        "",
        "## Adaptive Cells",
        "",
        "Set",
        "",
        "```text",
        "N(x)=ceil((1+x)^(3/4)).",
        "C_N=((N-1)^(4/3)-1,N^(4/3)-1].",
        "```",
        "",
        "On `C_N`, the all-`N` `d0/d1` budgets are constant. For `x>=245`,",
        "`exp(pi*x/8)/x^5` and `exp(pi*x/8)/x^6` increase, so each",
        "gamma-normalized error is largest at the right endpoint.",
        "",
        "## Cofinal Monotonicity",
        "",
        "The forward and compact-reflected endpoint bounds have decreasing",
        "geometric ratios, and their Gaussian exponents dominate the",
        "increment of `exp(pi*x_N/8)`. For the large reflected term, write",
        "",
        "```text",
        "f_A(y)=y^A*exp(-c_safe*y^(4/3)),",
        "c_safe=3*3^(2/3)*pi^(2/3)/16.",
        "```",
        "",
        "The normalized supremum-plus-integral tail decreases whenever",
        "",
        "```text",
        "(4/3)*(c_safe-pi/8)*y^(1/3)-A/y>0.",
        "```",
        "",
        "The worst case `A=22,y=63` is strictly positive by directed Arb.",
        "All Hermite and Leibniz coefficients are positive, so the complete",
        "`d0/d1` budgets inherit this decrease.",
        "",
        "## Boundary Witnesses",
        "",
        "| N | right x | J error / gamma | J' error / gamma |",
        "|---:|---:|---:|---:|",
    ]
    for row in artifact["boundary_witnesses"]:
        lines.append(
            f"| {row['N']} | `{row['x_right']['enclosure']}` | "
            f"`{row['J_error_over_gamma']['enclosure']}` | "
            f"`{row['J_prime_error_over_gamma']['enclosure']}` |"
        )
    lines.extend(
        [
            "",
            "`N=62` is an explicit non-promotion guard. At `N=63`, both",
            "ratios are below `1/4`, and the monotonicity theorem propagates",
            "that bound through every later adaptive cell.",
            "",
            "## Theorem",
            "",
            "For every `0<=t<=1/5` and every real `x>=245`,",
            "",
            "```text",
            "N=N(x),",
            "|J_t-J_(N,t)|<exp(-pi*x/8)/4,",
            "|J_t'-J_(N,t)'|<exp(-pi*x/8)/4.",
            "```",
            "",
            "The remaining cofinal obligation is a retained first-jet lower",
            "margin strong enough to dominate this quarter-gamma square.",
            "The finite bridge `69<x<245` is also not certified here.",
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
        "wrote Newman theta adaptive gamma-scale tail gate: "
        f"{len(artifact['rows'])} rows, "
        f"{len(artifact['boundary_witnesses'])} boundary witnesses"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
