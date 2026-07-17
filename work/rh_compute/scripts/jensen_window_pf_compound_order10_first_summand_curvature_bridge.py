#!/usr/bin/env python3
"""Transfer order-ten seventh-nested curvature from first summand to full kernel."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_first_summand_curvature_bridge.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order10_first_summand_curvature_bridge.md"
)
TAIL_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_m100_tail_curvature_reduction.json"
)
POWER12_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_negative_lambda_first_summand_power12_rebalanced_dominance_extension.json"
)
ORDER9_ENTRY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_m100_entry_certificate.json"
)
ORDER9_BRIDGE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_first_summand_curvature_bridge.json"
)
ORDER4_BRIDGE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_first_summand_curvature_bridge.json"
)
ORDER4_ENTRY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_m100_entry_certificate.json"
)
FINITE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_m100_finite_splice_certificate.json"
)

POWER_START = 320
B_ERROR_START = POWER_START + 1
FIRST_COORDINATE_START = 1251
TAIL_FIRST_K = 1252
FIRST_CONTINUOUS_CONSTANT = 4200
FIRST_DISCRETE_CONSTANT = 4201
FULL_TRANSFER_CONSTANT = 10
FULL_CURVATURE_CONSTANT = 4211
TARGET_CURVATURE_CONSTANT = 5500


@dataclass(frozen=True)
class BridgeRow:
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
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def source_record(path: Path, artifact: dict) -> dict:
    return {
        "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
        "sha256": sha256(path),
        "kind": artifact.get("kind"),
        "status": artifact.get("status"),
    }


def validate_sources() -> dict:
    tail = load_json(TAIL_SOURCE)
    power12 = load_json(POWER12_SOURCE)
    order9 = load_json(ORDER9_ENTRY_SOURCE)
    order9_bridge = load_json(ORDER9_BRIDGE_SOURCE)
    order4_bridge = load_json(ORDER4_BRIDGE_SOURCE)
    order4_entry = load_json(ORDER4_ENTRY_SOURCE)
    finite = load_json(FINITE_SOURCE)

    if tail.get("exact", {}).get("sufficient_ceiling") != (
        "Z_k<=5500/k^2 for every integer k>=1252"
    ):
        raise RuntimeError("order-ten tail reduction contract changed")
    if power12.get("summary", {}).get("full_tail_power") != 12:
        raise RuntimeError("power-twelve dominance exponent changed")
    if power12.get("summary", {}).get("tail_start_k") != POWER_START:
        raise RuntimeError("power-twelve dominance start changed")
    power12_bound = power12.get("diagnostics", {}).get(
        "full_tail_relative_bound"
    )
    if power12_bound != (
        "0<=delta_k=(M_k-M_k^(1))/M_k^(1)<2/k^12 "
        "for every integer k>=320"
    ):
        raise RuntimeError("power-twelve dominance contract changed")
    if order9.get("exact", {}).get("first_discrete_ceiling") != (
        "w_1''(t)<=4200/t^2 => Y_k^(1)<=4200*[-log(1-1/k^2)]"
        "<4201/k^2, k>=1251"
    ):
        raise RuntimeError("order-nine first curvature ceiling changed")
    if order9.get("exact", {}).get("full_curvature_ceiling") != (
        "w_1''(t)<=4200/t^2 on t>=1250 => "
        "Y_k<4201/k^2+550/k^2=4751/k^2<4900/k^2, k>=1251"
    ):
        raise RuntimeError("order-nine full curvature ceiling changed")
    if order9_bridge.get("exact", {}).get("sixth_gap_floor") != (
        "min(V_j,V_j^(1))>=4/(3*j), j>=1250"
    ):
        raise RuntimeError("order-nine sixth-gap floor changed")
    if order9_bridge.get("exact", {}).get("order8_coordinate_error") != (
        "C_j=2*N_j+Y1_j+(2*j/3)*P_j; |s_j-s_j^(1)|<=C_j"
    ):
        raise RuntimeError("order-eight coordinate transfer changed")
    if order9_bridge.get("exact", {}).get("order9_coordinate_error") != (
        "F_j=2*C_j+N_j+(3*j/4)*D_j; |w_j-w_j^(1)|<=F_j"
    ):
        raise RuntimeError("order-nine coordinate transfer changed")
    if order4_bridge.get("exact", {}).get("integer_B_bounds") != (
        "1/(2*m+1)<=B_1(m)<=3/(2*m-1), integer m>=319"
    ):
        raise RuntimeError("first-summand B floor changed")
    if order4_entry.get("tail_arithmetic", {}).get("defect_buffer") != (
        "-3*log(x_k)>3*d_k>=753/(250*(2*k+1))"
    ):
        raise RuntimeError("full-kernel B floor changed")
    if finite.get("exact", {}).get("combined_positive_block") != (
        "Q_(10,n)(-100)>0 for every 4<=n<=1242"
    ):
        raise RuntimeError("order-ten finite positive block changed")
    if finite.get("exact", {}).get("preserved_negative_prefix") != (
        "Q_(10,n)(-100)<0 for n=0,1,2,3"
    ):
        raise RuntimeError("order-ten negative prefix changed")

    artifacts = (
        (TAIL_SOURCE, tail),
        (POWER12_SOURCE, power12),
        (ORDER9_ENTRY_SOURCE, order9),
        (ORDER9_BRIDGE_SOURCE, order9_bridge),
        (ORDER4_BRIDGE_SOURCE, order4_bridge),
        (ORDER4_ENTRY_SOURCE, order4_entry),
        (FINITE_SOURCE, finite),
    )
    return {
        "sources": [source_record(path, artifact) for path, artifact in artifacts],
        "tail_ceiling": tail["exact"]["sufficient_ceiling"],
        "power12_bound": power12_bound,
        "power12_B_error": power12["diagnostics"]["adjacent_B_error"],
        "order9_first_ceiling": order9["exact"]["first_discrete_ceiling"],
        "order9_full_ceiling": order9["exact"]["full_curvature_ceiling"],
        "first_B_floor": order4_bridge["exact"]["integer_B_bounds"],
        "full_B_floor": "B(j)>=d_j>=251/(250*(2*j+1)), j>=320",
        "finite_positive": finite["exact"]["combined_positive_block"],
        "finite_negative": finite["exact"]["preserved_negative_prefix"],
    }


def shifted_positive_polynomial(
    expression: sp.Expr, variable: sp.Symbol, start: int
) -> dict:
    numerator, denominator = sp.fraction(sp.factor(expression))
    shift = sp.symbols("m", integer=True, nonnegative=True)
    polynomial = sp.Poly(sp.expand(numerator.subs(variable, start + shift)), shift)
    coefficients = polynomial.all_coeffs()
    if any(coefficient <= 0 for coefficient in coefficients):
        raise RuntimeError(
            f"shifted positivity failed at {variable}>={start}: {polynomial.as_expr()}"
        )
    return {
        "start": start,
        "degree": polynomial.degree(),
        "coefficient_count": len(coefficients),
        "leading_coefficient": str(coefficients[0]),
        "constant_coefficient": str(coefficients[-1]),
        "minimum_coefficient": str(min(coefficients)),
        "denominator_factorization": str(sp.factor(denominator)),
        "shifted_coefficients": [str(value) for value in coefficients],
    }


def stencil_factor(power: int, start: int) -> Fraction:
    return Fraction(start, start - 1) ** power + 3


def decimal_text(value: Fraction, digits: int = 50) -> str:
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(value.numerator) / Decimal(value.denominator), "E")


def envelope_row(
    name: str,
    symbol: str,
    start: int,
    power: int,
    constant: Fraction,
    recurrence: str,
) -> dict:
    return {
        "name": name,
        "symbol": symbol,
        "start": start,
        "power": power,
        "constant": str(constant),
        "constant_decimal": decimal_text(constant),
        "bound": f"{symbol}_j<={constant}/j^{power}, j>={start}",
        "recurrence": recurrence,
    }


def rational_power_envelope() -> dict:
    # If E_m<=C/m^p for m>=J-1, its positive centered stencil is at most
    # C*((J/(J-1))^p+3)/j^p for every j>=J.
    a = 2 * (Fraction(321, 320) ** 12 + 3)
    ell = 4 * a

    j = 322
    first_gap = 2 * a / j + ell * stencil_factor(11, j)
    log_first_gap = 2 * ell / j + 8 * first_gap

    j = 323
    second_gap = 3 * a / j**2 + log_first_gap * stencil_factor(10, j)
    order5 = (
        2 * log_first_gap / j
        + Fraction(5, 7) * second_gap
        + ell / j**2
    )

    j = 324
    third_gap = 4 * a / j**3 + order5 * stencil_factor(9, j)
    order6 = (
        2 * order5 / j
        + log_first_gap / j**2
        + Fraction(2, 3) * third_gap
    )

    j = 325
    fourth_gap = 5 * a / j**4 + order6 * stencil_factor(8, j)
    order7 = (
        2 * order6 / j
        + order5 / j**2
        + Fraction(2, 3) * fourth_gap
    )

    j = 326
    fifth_gap = 6 * a / j**5 + order7 * stencil_factor(7, j)

    # The inherited fifth-gap floor begins at j=1249, so the order-eight
    # stable-log envelope starts there even though its inputs hold earlier.
    j = 1249
    order8 = (
        2 * order7 / j
        + order6 / j**2
        + Fraction(2, 3) * fifth_gap
    )

    j = 1250
    sixth_gap = 7 * a / j**6 + order8 * stencil_factor(6, j)
    order9 = (
        2 * order8 / j
        + order7 / j**2
        + Fraction(3, 4) * sixth_gap
    )

    j = 1251
    seventh_gap = 8 * a / j**7 + order9 * stencil_factor(5, j)

    j = FIRST_COORDINATE_START
    order10 = 2 * order9 / j + order8 / j**2 + 5 * seventh_gap

    k = TAIL_FIRST_K
    transfer_scaled = order10 * stencil_factor(4, k) / k**2
    if transfer_scaled >= FULL_TRANSFER_CONSTANT:
        raise RuntimeError("order-ten rational transfer envelope exhausted")

    rows = [
        envelope_row(
            "moment wall",
            "a",
            321,
            12,
            a,
            "a_j=2*((j-1)^(-12)+2*j^(-12)+(j+1)^(-12))",
        ),
        envelope_row("log defect", "L", 321, 11, ell, "L_j=4*j*a_j"),
        envelope_row(
            "first gap", "U1", 322, 11, first_gap, "U1_j=2*a_j+stencil(L)_j"
        ),
        envelope_row(
            "log first gap",
            "V1",
            322,
            10,
            log_first_gap,
            "V1_j=2*L_j+8*j*U1_j",
        ),
        envelope_row(
            "second gap",
            "W1",
            323,
            10,
            second_gap,
            "W1_j=3*a_j+stencil(V1)_j",
        ),
        envelope_row(
            "order-five coordinate",
            "E",
            323,
            9,
            order5,
            "E_j=2*V1_j+(5*j/7)*W1_j+L_j",
        ),
        envelope_row(
            "third gap", "Z", 324, 9, third_gap, "Z_j=4*a_j+stencil(E)_j"
        ),
        envelope_row(
            "order-six coordinate",
            "Y1",
            324,
            8,
            order6,
            "Y1_j=2*E_j+V1_j+(2*j/3)*Z_j",
        ),
        envelope_row(
            "fourth gap", "O", 325, 8, fourth_gap, "O_j=5*a_j+stencil(Y1)_j"
        ),
        envelope_row(
            "order-seven coordinate",
            "N",
            325,
            7,
            order7,
            "N_j=2*Y1_j+E_j+(2*j/3)*O_j",
        ),
        envelope_row(
            "fifth gap", "P", 326, 7, fifth_gap, "P_j=6*a_j+stencil(N)_j"
        ),
        envelope_row(
            "order-eight coordinate",
            "C",
            1249,
            6,
            order8,
            "C_j=2*N_j+Y1_j+(2*j/3)*P_j",
        ),
        envelope_row(
            "sixth gap", "D", 1250, 6, sixth_gap, "D_j=7*a_j+stencil(C)_j"
        ),
        envelope_row(
            "order-nine coordinate",
            "F",
            1250,
            5,
            order9,
            "F_j=2*C_j+N_j+(3*j/4)*D_j",
        ),
        envelope_row(
            "seventh gap",
            "D7",
            1251,
            5,
            seventh_gap,
            "D7_j=8*a_j+stencil(F)_j",
        ),
        envelope_row(
            "order-ten coordinate",
            "G",
            FIRST_COORDINATE_START,
            4,
            order10,
            "G_j=2*F_j+C_j+5*j*D7_j",
        ),
    ]
    return {
        "stencil_lemma": (
            "E_m<=C/m^p for m>=J-1 implies "
            "E_(j-1)+2*E_j+E_(j+1)<=C*((J/(J-1))^p+3)/j^p, j>=J"
        ),
        "rows": rows,
        "transfer_scaled_exact": str(transfer_scaled),
        "transfer_scaled_decimal": decimal_text(transfer_scaled),
        "transfer_reserve_exact": str(Fraction(FULL_TRANSFER_CONSTANT) - transfer_scaled),
        "transfer_reserve_decimal": decimal_text(
            Fraction(FULL_TRANSFER_CONSTANT) - transfer_scaled
        ),
        "transfer_bound": (
            "|Z_k-Z_k^(1)|<=G_(k-1)+2*G_k+G_(k+1)<10/k^2, k>=1252"
        ),
    }


def exact_diagnostics() -> dict:
    j = sp.symbols("j", integer=True, positive=True)
    first_floor = shifted_positive_polynomial(
        8 / (2 * j + 1) - sp.Rational(4201, 1) / j**2 - 1 / (5 * j),
        j,
        FIRST_COORDINATE_START,
    )
    full_floor = shifted_positive_polynomial(
        sp.Rational(1004, 125) / (2 * j + 1)
        - sp.Rational(4751, 1) / j**2
        - 1 / (5 * j),
        j,
        FIRST_COORDINATE_START,
    )

    t = sp.symbols("t", positive=True)
    phi = 1 / (sp.exp(t) - 1)
    derivative_residual = sp.simplify(
        sp.diff(sp.log(1 - sp.exp(-t)), t) - phi
    )
    if derivative_residual != 0:
        raise RuntimeError("stable logarithm derivative identity failed")

    return {
        "continuous_coordinates": (
            "W(t)=8*B(t)-w(t-1)+2*w(t)-w(t+1); "
            "z(t)=2*w(t)-s(t)+log(1-exp(-W(t)))"
        ),
        "discrete_identity": (
            "Z_k=z(k-1)-2*z(k)+z(k+1)="
            "integral_[-1,1](1-|v|)*z''(k+v) dv"
        ),
        "first_seventh_gap_floor": (
            "W_j^(1)>=8/(2*j+1)-4201/j^2>1/(5*j), j>=1251"
        ),
        "full_seventh_gap_floor": (
            "W_j>=1004/(125*(2*j+1))-4751/j^2>1/(5*j), j>=1251"
        ),
        "seventh_gap_floor": "min(W_j,W_j^(1))>1/(5*j), j>=1251",
        "stable_log_lipschitz": (
            "phi(min(W_j,W_j^(1)))<=1/min(W_j,W_j^(1))<5*j"
        ),
        "first_floor_polynomial": first_floor,
        "full_floor_polynomial": full_floor,
        "power_envelope": rational_power_envelope(),
        "continuous_target": "z_1''(t)<=4200/t^2 for every real t>=1251",
        "tent_transfer": (
            "z_1''(t)<=4200/t^2 => Z_k^(1)<=4200*[-log(1-1/k^2)]"
            "<4201/k^2, k>=1252"
        ),
        "full_transfer": (
            "|Z_k-Z_k^(1)|<10/k^2 for every integer k>=1252"
        ),
        "conditional_full_ceiling": (
            "z_1''(t)<=4200/t^2 on t>=1251 => "
            "Z_k<4201/k^2+10/k^2=4211/k^2<5500/k^2, k>=1252"
        ),
        "conditional_endpoint_tail": (
            "z_1''(t)<=4200/t^2 on t>=1251 => "
            "Q_(10,n)(-100)>0 for every n>=1243"
        ),
        "conditional_delayed_entry": (
            "z_1''(t)<=4200/t^2 on t>=1251 => "
            "Q_(10,n)(-100)>0 for every n>=4"
        ),
        "preserved_negative_prefix": "Q_(10,n)(-100)<0 for n=0,1,2,3",
        "stable_derivative_residual": str(derivative_residual),
    }


def build_artifact() -> dict:
    sources = validate_sources()
    exact = exact_diagnostics()
    envelope = exact["power_envelope"]
    rows = [
        BridgeRow(
            "co10fscb_01_continuous_coordinate",
            "exact_identity",
            "ready_to_apply",
            "The seventh stable gap gives a cancellation-free coordinate for Q9 and its order-ten centered curvature.",
            exact["continuous_coordinates"] + "; " + exact["discrete_identity"],
            "Exact stable-coordinate algebra only.",
        ),
        BridgeRow(
            "co10fscb_02_power12_input",
            "theorem_input",
            "ready_to_apply",
            "Effective inverse-twelfth-power first-summand dominance supplies the moment-wall perturbation.",
            sources["power12_bound"] + "; " + sources["power12_B_error"],
            "Complete versus first Newman summand at lambda=-100.",
        ),
        BridgeRow(
            "co10fscb_03_seventh_gap_floor",
            "exact_rational_theorem",
            "ready_to_apply",
            "Completed order-nine curvature and the coefficient-defect floor keep both seventh gaps uniformly positive.",
            exact["first_seventh_gap_floor"] + "; " + exact["full_seventh_gap_floor"],
            "Integer j>=1251 only.",
            {
                "first_polynomial": exact["first_floor_polynomial"],
                "full_polynomial": exact["full_floor_polynomial"],
            },
        ),
        BridgeRow(
            "co10fscb_04_power_envelope",
            "exact_rational_theorem",
            "ready_to_apply",
            "A staggered rational power envelope propagates the twelfth-power wall through all seven stable logarithms.",
            envelope["stencil_lemma"] + "; " + exact["stable_log_lipschitz"],
            "Triangle inequalities with exact positive rational constants.",
            {
                "envelope_rows": len(envelope["rows"]),
                "transfer_scaled": envelope["transfer_scaled_decimal"],
                "reserve_below_10": envelope["transfer_reserve_decimal"],
            },
        ),
        BridgeRow(
            "co10fscb_05_full_transfer",
            "exact_transfer_theorem",
            "ready_to_apply",
            "The complete theta kernel changes the order-ten centered curvature by fewer than ten inverse squares.",
            exact["full_transfer"],
            "Complete versus first Newman summand at lambda=-100, k>=1252.",
            {
                "scaled_envelope_exact": envelope["transfer_scaled_exact"],
                "scaled_envelope_decimal": envelope["transfer_scaled_decimal"],
            },
        ),
        BridgeRow(
            "co10fscb_06_tent_transfer",
            "conditional_exact_transfer",
            "ready_to_apply",
            "A 4200 continuous first-summand ceiling fits inside the endpoint target after tent integration and full-kernel transfer.",
            exact["tent_transfer"] + "; " + exact["conditional_full_ceiling"],
            "Conditional only on the displayed continuous first-summand theorem.",
        ),
        BridgeRow(
            "co10fscb_07_endpoint_splice",
            "conditional_theorem_composition",
            "ready_to_apply",
            "The conditional analytic tail meets the completed positive block while preserving the four proved negative endpoint rows.",
            exact["conditional_endpoint_tail"]
            + "; "
            + sources["finite_positive"]
            + "; "
            + exact["preserved_negative_prefix"],
            "No claim for n=0,1,2,3 and no RH/Jensen/Lambda conclusion.",
        ),
        BridgeRow(
            "co10fscb_08_open_continuous_theorem",
            "open_analytic_target",
            "open",
            "Complete the global continuous first-summand ceiling with a strict 4200 constant.",
            exact["continuous_target"],
            "Compact interval computation remains the only open premise in this bridge.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order10_first_summand_curvature_bridge",
        "date": "2026-07-16",
        "status": (
            "exact order-ten first/full curvature transfer with one open global "
            "first-summand curvature theorem"
        ),
        "proof_boundary": (
            "This artifact proves the seventh-gap floor and full-kernel transfer, "
            "and reduces the endpoint tail to the displayed continuous first-summand "
            "theorem. It does not yet prove that theorem, order-ten delayed entry, "
            "PF-infinity, RH, or Lambda<=0."
        ),
        "source_contract": sources,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": sum(row.readiness == "ready_to_apply" for row in rows),
            "open_rows": sum(row.readiness == "open" for row in rows),
            "seventh_gap_floor_theorems": 1,
            "full_kernel_transfer_theorems": 1,
            "power_envelope_rows": len(envelope["rows"]),
            "transfer_constant": FULL_TRANSFER_CONSTANT,
            "conditional_full_constant": FULL_CURVATURE_CONSTANT,
            "target_constant": TARGET_CURVATURE_CONSTANT,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order10_first_summand_curvature_bridge.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order10_first_summand_curvature_bridge.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    envelope = exact["power_envelope"]
    lines = [
        "# Jensen-Window PF Order-Ten First/Full Curvature Bridge",
        "",
        "Date: 2026-07-16",
        "",
        "Status: exact full-kernel transfer with one open continuous",
        "first-summand theorem. This is not a proof of PF-infinity, RH, or",
        "`Lambda <= 0`, and it is not an order-ten entry theorem.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order10_first_summand_curvature_bridge.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order10_first_summand_curvature_bridge.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order10_first_summand_curvature_bridge.py",
        "```",
        "",
        "## Stable Floor",
        "",
        "```text",
        exact["first_seventh_gap_floor"],
        exact["full_seventh_gap_floor"],
        exact["seventh_gap_floor"],
        exact["stable_log_lipschitz"],
        "```",
        "",
        "Both floor inequalities are certified by shifted polynomials with",
        "strictly positive coefficients.",
        "",
        "## Rational Error Envelope",
        "",
        "Starting from the proved inverse-twelfth-power moment defect, the exact",
        "positive-stencil lemma propagates `C/j^p` bounds through the seven",
        "stable logarithms. The final scaled transfer envelope is",
        "",
        "```text",
        "k^2*|Z_k-Z_k^(1)| <= " + envelope["transfer_scaled_decimal"],
        exact["full_transfer"],
        "```",
        "",
        "The power-ten estimate used at order nine is not asserted to be",
        "impossible here; the power-twelve input is the effective estimate that",
        "makes this particular rigorous envelope close cleanly.",
        "",
        "## Remaining Premise",
        "",
        "```text",
        exact["continuous_target"],
        exact["tent_transfer"],
        exact["conditional_full_ceiling"],
        exact["conditional_endpoint_tail"],
        "```",
        "",
        "Together with the finite block, that premise would prove positivity only",
        "for `n>=4`; the four negative endpoint rows remain negative.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-ten first/full curvature bridge: "
        f"{summary['power_envelope_rows']} rational envelope rows, "
        f"transfer <{summary['transfer_constant']}/k^2, "
        f"{summary['open_rows']} open continuous theorem"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
