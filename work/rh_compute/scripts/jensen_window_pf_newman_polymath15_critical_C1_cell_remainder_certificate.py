#!/usr/bin/env python3
"""Certify a corrected C1 remainder on critical fixed-cutoff cells."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_C1_cell_remainder_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_critical_C1_cell_remainder_certificate.md"
)
SOURCE_URL = "https://arxiv.org/abs/1904.12438"
PRECISION_BITS = 256
L_MIN = 50
AB_CONSTANT = 1000
EC_CONSTANT = 100
LIFT_CONSTANT = 100
RAW_CONSTANT = 2500
FIRST_CONSTANT = 5000
NORM_SQUARED_CONSTANT = 32_000_000


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def arb_text(value: flint.arb, digits: int = 70) -> str:
    return value.str(digits)


def interval_budget() -> dict:
    flint.ctx.prec = PRECISION_BITS
    ell = flint.arb(L_MIN)
    pi = flint.arb.pi()
    exp_one = flint.arb(1).exp()

    coefficient_constant = (
        6 * 4 * flint.arb(2).sqrt() * (flint.arb(1) / 4).exp()
    )
    if not coefficient_constant < 50:
        raise RuntimeError("critical coefficient-sum constant failed")

    x_lower = 4 * pi * (ell - 1).exp()
    u_scaled = 42 * ell.exp() / (x_lower - flint.arb(7))
    if not u_scaled < 10:
        raise RuntimeError("published eAB exponent constant failed")
    eab_constant = 50 * 2 * u_scaled
    if not eab_constant < AB_CONSTANT:
        raise RuntimeError("published eAB collar constant failed")

    ec_positive_exponent = (
        3 * ((ell + 1) ** 2 + pi**2 / 4).sqrt() + flint.arb(179) / 50
    ) / (x_lower - flint.arb(9))
    if not ec_positive_exponent < 1:
        raise RuntimeError("published eC positive exponent failed")
    ec_bracket_constant = 5 * (flint.arb(1) / 2).exp() + (
        flint.arb(173) / 25
    ) * (ell / 2).exp() / (x_lower - 12)
    if not ec_bracket_constant < 10:
        raise RuntimeError("published eC saddle bracket failed")
    ec_constant = (
        (flint.arb(1) / 4).exp()
        * ec_positive_exponent.exp()
        * ec_bracket_constant
    )
    if not ec_constant < EC_CONSTANT:
        raise RuntimeError("published eC collar constant failed")

    normalizer_ratio = (flint.arb(1) / 2).exp()
    if not normalizer_ratio < 2:
        raise RuntimeError("normalizer collar ratio failed")

    # |S'| <= 58*K*|T| follows from |C0|<5 and |C0'|<60. The
    # two sharp-conjugate components, |T|/T0<2, endpoint prefactor <3,
    # and sqrt(T0)>2*exp(L/2) leave the following normalized source cap.
    lift_source_constant = (
        flint.arb(58)
        * 2
        * 3
        / 2
        * (pi / (8 * ell)).exp()
    )
    if not lift_source_constant < 300:
        raise RuntimeError("holomorphic-lift source constant failed")
    lift_constant = (
        lift_source_constant * (pi / (8 * ell)).exp() / ell
    )
    if not lift_constant < LIFT_CONSTANT:
        raise RuntimeError("holomorphic-lift defect constant failed")

    combined_raw_constant = (
        2 * (flint.arb(AB_CONSTANT) + EC_CONSTANT) + LIFT_CONSTANT
    )
    if not combined_raw_constant < RAW_CONSTANT:
        raise RuntimeError("combined raw collar constant failed")
    first_constant = 2 * flint.arb(RAW_CONSTANT)
    if not first_constant <= FIRST_CONSTANT:
        raise RuntimeError("normalized first-jet constant failed")
    norm_squared_constant = (
        flint.arb(RAW_CONSTANT) ** 2 + flint.arb(FIRST_CONSTANT) ** 2
    )
    if not norm_squared_constant < NORM_SQUARED_CONSTANT:
        raise RuntimeError("transversality error constant failed")

    endpoint_error = (
        flint.arb(NORM_SQUARED_CONSTANT) * (-3 * ell / 2).exp()
    )
    if not endpoint_error < flint.arb(1) / 10**24:
        raise RuntimeError("L=50 transversality error endpoint failed")
    return {
        "precision_bits": PRECISION_BITS,
        "L_min": str(L_MIN),
        "coefficient_sum_constant_ball": arb_text(coefficient_constant),
        "coefficient_sum_constant_lt": "50",
        "eAB_exponent_scaled_ball": arb_text(u_scaled),
        "eAB_exponent_scaled_lt": "10",
        "eAB_constant_ball": arb_text(eab_constant),
        "eAB_constant_lt": str(AB_CONSTANT),
        "eC_positive_exponent_ball": arb_text(ec_positive_exponent),
        "eC_positive_exponent_lt": "1",
        "eC_bracket_constant_ball": arb_text(ec_bracket_constant),
        "eC_bracket_constant_lt": "10",
        "eC_constant_ball": arb_text(ec_constant),
        "eC_constant_lt": str(EC_CONSTANT),
        "normalizer_ratio_ball": arb_text(normalizer_ratio),
        "normalizer_ratio_lt": "2",
        "lift_source_constant_ball": arb_text(lift_source_constant),
        "lift_source_constant_lt": "300",
        "lift_constant_ball": arb_text(lift_constant),
        "lift_constant_lt": str(LIFT_CONSTANT),
        "combined_raw_constant_ball": arb_text(combined_raw_constant),
        "combined_raw_constant_lt": str(RAW_CONSTANT),
        "first_constant": str(FIRST_CONSTANT),
        "norm_squared_constant_ball": arb_text(norm_squared_constant),
        "norm_squared_constant_lt": str(NORM_SQUARED_CONSTANT),
        "L50_error_norm_squared_ball": arb_text(endpoint_error),
        "L50_error_norm_squared_lt": "10^-24",
    }


def build_exact() -> dict:
    return {
        "region": (
            "L=log(x/(4*pi))>=50, 0<t<=1/2, 0<t*L<=25, "
            "rho=1/L, and the radius-rho disk stays in one prescribed-N cell"
        ),
        "collar_geometry": (
            "On |z-x|<=1/L: X=Re z satisfies L-1<=log(X/(4*pi))<=L+1 "
            "and 0<=|Im z|<=1/L"
        ),
        "coefficient_sum": (
            "sum_(n<=N)(1+|gamma|N^|kappa|n^y)b_n^t*n^(-Re s_*) "
            "<=50*exp(L/4)"
        ),
        "eAB": "e_A+e_B<1000*exp(-3L/4)",
        "eC": "e_C<100*exp(-3L/4)",
        "normalizer": "sup_collar |B_t(z)|/A_t(x)<2",
        "lift": (
            "|C_hat_(N,t)(z)-C_paper(z)|/A_t(x)<"
            "100*exp(-3L/4)"
        ),
        "analytic_remainder": (
            "R_hat=H_t-(A_(t,N)+B_(t,N)-C_hat_(N,t)) is holomorphic and "
            "sup_collar |R_hat|/A_t(x)<2500*exp(-3L/4)"
        ),
        "normalized_split": (
            "Z_t=J_hat_(N,t)+r, where J_hat agrees on the real axis with "
            "the corrected main J=P-Q"
        ),
        "c1_caps": (
            "|r|<2500*exp(-3L/4), "
            "|r'|<5000*L*exp(-3L/4)"
        ),
        "error_norm": (
            "r^2+(r'/L)^2<32000000*exp(-3L/2)"
        ),
        "cell_target": (
            "T_L[J]>32000000*exp(-3L/2) excludes a double zero on a fixed-N cell"
        ),
        "transition_gap": (
            "Disks crossing a cutoff still require an adjacent-N corrected-lift comparison"
        ),
    }


def build_artifact() -> dict:
    exact = build_exact()
    interval = interval_budget()
    rows = [
        GateRow(
            id="np15c1crc_01_region",
            role="exact_scope",
            readiness="ready_to_apply",
            claim="The critical asymptotic cell region has a scale-adapted holomorphic collar.",
            formula=f"{exact['region']}; {exact['collar_geometry']}",
            proof_boundary="Fixed prescribed cutoff N throughout the disk.",
        ),
        GateRow(
            id="np15c1crc_02_coefficient_sum",
            role="analytic_bound",
            readiness="ready_to_apply",
            claim="The complete two-saddle Dirichlet coefficient mass has a uniform square-root cutoff bound.",
            formula=exact["coefficient_sum"],
            proof_boundary="Uses the published Re(s_*) and gamma/kappa estimates.",
            diagnostics={
                "constant": interval["coefficient_sum_constant_ball"],
            },
        ),
        GateRow(
            id="np15c1crc_03_published_eAB",
            role="interval_certificate",
            readiness="certified",
            claim="The published finite-sum approximation errors gain three quarters of the logarithmic height.",
            formula=exact["eAB"],
            proof_boundary="Uniform for bounded c=tL on the stated collar.",
            diagnostics={
                "exponent": interval["eAB_exponent_scaled_ball"],
                "constant": interval["eAB_constant_ball"],
            },
        ),
        GateRow(
            id="np15c1crc_04_published_eC",
            role="interval_certificate",
            readiness="certified",
            claim="After extracting the endpoint term, the published contour error gains one saddle factor.",
            formula=exact["eC"],
            proof_boundary="Uses the refined e_C, not the coarser e_C0.",
            diagnostics={
                "positive_exponent": interval["eC_positive_exponent_ball"],
                "constant": interval["eC_constant_ball"],
            },
        ),
        GateRow(
            id="np15c1crc_05_holomorphic_lift",
            role="analytic_transfer",
            readiness="ready_to_apply",
            claim="The certified C0 strip bound makes the endpoint holomorphic-lift defect lower order on the collar.",
            formula=exact["lift"],
            proof_boundary="Composes the exact lift defect equation with |C0|<5 and |C0'|<60.",
            diagnostics={"constant": interval["lift_constant_ball"]},
        ),
        GateRow(
            id="np15c1crc_06_raw_collar",
            role="proved_bound",
            readiness="ready_to_apply",
            claim="The corrected fixed-cutoff approximation has one holomorphic scalar remainder bound.",
            formula=f"{exact['normalizer']}; {exact['analytic_remainder']}",
            proof_boundary="Fixed-cutoff cell interiors only.",
            diagnostics={
                "normalizer": interval["normalizer_ratio_ball"],
                "raw_constant": interval["combined_raw_constant_ball"],
            },
        ),
        GateRow(
            id="np15c1crc_07_C1_transfer",
            role="proved_bound",
            readiness="ready_to_apply",
            claim="Cauchy's estimate and normalizer differentiation give an explicit corrected C1 remainder.",
            formula=f"{exact['normalized_split']}; {exact['c1_caps']}",
            proof_boundary="Uses |(log A_t)'|<L/2 and collar radius 1/L.",
        ),
        GateRow(
            id="np15c1crc_08_transversality_budget",
            role="proved_bound",
            readiness="ready_to_apply",
            claim="The exact double-zero exclusion target now has a fully explicit fixed-cell error side.",
            formula=f"{exact['error_norm']}; {exact['cell_target']}",
            proof_boundary="The corrected-main lower bound remains open.",
            diagnostics={
                "constant": interval["norm_squared_constant_ball"],
                "L50_value": interval["L50_error_norm_squared_ball"],
            },
        ),
        GateRow(
            id="np15c1crc_09_transition_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="Cutoff-crossing Cauchy disks remain outside this fixed-cell certificate.",
            formula=exact["transition_gap"],
            proof_boundary="Needs one adjacent corrected-lift comparison.",
        ),
        GateRow(
            id="np15c1crc_10_nonpromotion",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="Closing the remainder side does not establish corrected finite transversality.",
            formula="explicit error budget != T_L[J] lower bound",
            proof_boundary="Not Lambda<=0 or RH.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_polymath15_critical_C1_cell_remainder_certificate",
        "date": "2026-07-17",
        "status": (
            "explicit corrected C1 remainder certificate on critical fixed-cutoff "
            "cells with L>=50; not corrected transversality, Lambda<=0, or RH"
        ),
        "proof_boundary": (
            "This artifact proves the stated corrected normalized C1 remainder "
            "bounds on radius-1/L disks contained in one prescribed-N cell. It "
            "does not cover cutoff-crossing disks, prove the corrected-main "
            "first-jet lower bound, close lower frequencies, exclude a positive "
            "Newman boundary, prove Lambda<=0, or prove RH."
        ),
        "exact": exact,
        "interval": interval,
        "rows": [asdict(row) for row in rows],
        "sources": [
            SOURCE_URL,
            "outputs/jensen_window_pf_newman_polymath15_endpoint_holomorphic_lift.md",
            "outputs/jensen_window_pf_newman_polymath15_endpoint_C0_strip_certificate.md",
            "outputs/jensen_window_pf_newman_polymath15_critical_transversality_target.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    interval = artifact["interval"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Polymath-15 Critical C1 Cell Remainder Certificate",
            "",
            "Date: 2026-07-17",
            "",
            "Status: explicit corrected `C1` remainder on critical fixed-cutoff",
            "cells. This is not a proof of `Lambda <= 0` or RH; corrected",
            "transversality remains open.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_C1_cell_remainder_certificate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_C1_cell_remainder_certificate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_C1_cell_remainder_certificate.py",
            "```",
            "",
            "The scalar approximation comes from the A+B-C theorem in",
            f"[Polymath 15]({SOURCE_URL}). The endpoint lift and its `C0` strip",
            "bounds are supplied by the two preceding local certificates.",
            "",
            "## Region",
            "",
            "```text",
            exact["region"],
            exact["collar_geometry"],
            "```",
            "",
            "## Published Errors",
            "",
            "The coefficient envelope and explicit exponential bracket give",
            "",
            "```text",
            exact["coefficient_sum"],
            exact["eAB"],
            f"saved eAB constant = {interval['eAB_constant_ball']} < 1000",
            "```",
            "",
            "The refined endpoint remainder gives",
            "",
            "```text",
            exact["eC"],
            f"saved eC constant = {interval['eC_constant_ball']} < 100",
            "```",
            "",
            "Using `e_C` is essential: `e_C0` retains a leading factor one and",
            "would lose the extra saddle decay.",
            "",
            "## Analytic Collar",
            "",
            "The exact endpoint lift satisfies",
            "",
            "```text",
            exact["lift"],
            f"saved lift source constant = {interval['lift_source_constant_ball']} < 300",
            exact["normalizer"],
            exact["analytic_remainder"],
            "```",
            "",
            "The saved combined constant is",
            f"`{interval['combined_raw_constant_ball']} < 2500`.",
            "",
            "## C1 Transfer",
            "",
            "Cauchy's estimate on radius `1/L`, followed by differentiation of",
            "the positive real normalizer, yields",
            "",
            "```text",
            exact["normalized_split"],
            exact["c1_caps"],
            exact["error_norm"],
            "```",
            "",
            "At `L=50` the saved error norm is already",
            f"`{interval['L50_error_norm_squared_ball']} < 10^-24`.",
            "Thus the fixed-cell collision condition is excluded whenever",
            "",
            "```text",
            exact["cell_target"],
            "```",
            "",
            "## Remaining Work",
            "",
            "The error side is now explicit on fixed-cutoff cells. Two obligations",
            "remain separate: compare adjacent corrected lifts on cutoff-crossing",
            "disks, and establish the arithmetic lower bound for `T_L[J]`, which",
            "is the remaining RH-level obligation.",
            "Neither is supplied by this certificate.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman Polymath-15 critical C1 cell remainder certificate: "
        f"{len(artifact['rows'])} rows, 5 interval constants, "
        "1 explicit C1 budget, 1 transition handoff"
    )


if __name__ == "__main__":
    main()
