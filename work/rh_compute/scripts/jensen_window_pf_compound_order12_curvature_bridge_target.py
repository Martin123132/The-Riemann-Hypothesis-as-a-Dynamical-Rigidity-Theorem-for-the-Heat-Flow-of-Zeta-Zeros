#!/usr/bin/env python3
"""Derive the conditional order-twelve curvature and transfer target."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order10_first_summand_curvature_bridge as bridge_core  # noqa: E402


RESULTS = REPO_ROOT / "work/rh_compute/results"
POWER14_SOURCE = RESULTS / "jensen_window_pf_negative_lambda_first_summand_power14_rebalanced_dominance_extension.json"
ORDER11_CURVATURE_SOURCE = RESULTS / "jensen_window_pf_compound_order11_first_summand_curvature_certificate.json"
ORDER11_ENTRY_SOURCE = RESULTS / "jensen_window_pf_compound_order11_m100_entry_certificate.json"
ORDER10_ENTRY_SOURCE = RESULTS / "jensen_window_pf_compound_order10_m100_delayed_entry_certificate.json"
ORDER11_COMPLETION_SOURCE = RESULTS / "jensen_window_pf_compound_order11_lambda0_completion_certificate.json"
DELAYED_HEAT_SOURCE = RESULTS / "jensen_window_pf_delayed_cooperative_heat_tail_lemma.json"
DEFECT_SOURCE = RESULTS / "jensen_window_pf_compound_order4_m100_entry_certificate.json"
PARTIAL_PREFIX_SOURCE = RESULTS / "jensen_window_pf_compound_order12_m100_partial_prefix_certificate.json"
ENDPOINT_COMPLETION_SOURCE = RESULTS / "jensen_window_pf_compound_order12_m100_endpoint_completion_certificate.json"
LAMBDA0_PREFIX_SOURCE = RESULTS / "jensen_window_pf_compound_order12_lambda0_prefix_certificate.json"
DEFAULT_OUT = RESULTS / "jensen_window_pf_compound_order12_curvature_bridge_target.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_compound_order12_curvature_bridge_target.md"

POWER_START = 380
WALL_START = 381
GAP_FLOOR_START = 1503
TAIL_FIRST_K = 1504
TAIL_FIRST_N = 1493
CONTINUUM_CONSTANT = 8000
FIRST_DISCRETE_CONSTANT = 8001
TRANSFER_CONSTANT = 100
FULL_CEILING = 8101
SAFE_CEILING = 8200

Y_CONTINUOUS = "y_1''(t)<=6000/t^2 for every real t>=1252"
V_CONTINUOUS_TARGET = "v_1''(t)<=8000/t^2 for every real t>=1503"
V_FIRST_DISCRETE = "v_1''(t)<=8000/t^2 => V_k^(1)<=8000*[-log(1-1/k^2)]<8001/k^2, k>=1504"
V_FULL_TRANSFER = "|V_k-V_k^(1)|<100/k^2 for every integer k>=1504"
V_FULL_CEILING = "v_1''(t)<=8000/t^2 on t>=1503 => V_k<8001/k^2+100/k^2=8101/k^2<8200/k^2, k>=1504"
ORDER12_FINITE_PREFIX_THEOREM = "Q_(12,n)(-100)>0 for every 1241<=n<=1492"
ORDER12_LAMBDA0_PREFIX_TARGET = "Q_(12,n)(0)>0 for every integer 0<=n<=3"
ORDER11_HEAT_RAY = "Q_(11,n)(lambda)>0 for every n>=4 and -100<=lambda<=0"
ORDER12_DELAYED_HANDOFF = (
    "[Q_(11,n)(lambda)>0 for every n>=4 and -100<=lambda<=0 and "
    "Q_(12,n)(-100)>0 for every n>=4] implies "
    "Q_(12,n)(lambda)>0 for every n>=4 and -100<=lambda<=0"
)


@dataclass(frozen=True)
class TargetRow:
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


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def source_record(path: Path, artifact: dict) -> dict:
    return {
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "sha256": sha256(path),
        "kind": artifact.get("kind"),
        "status": artifact.get("status"),
    }


def validate_sources() -> list[dict]:
    power14 = load_json(POWER14_SOURCE)
    curvature = load_json(ORDER11_CURVATURE_SOURCE)
    entry11 = load_json(ORDER11_ENTRY_SOURCE)
    entry10 = load_json(ORDER10_ENTRY_SOURCE)
    completion11 = load_json(ORDER11_COMPLETION_SOURCE)
    delayed = load_json(DELAYED_HEAT_SOURCE)
    defect = load_json(DEFECT_SOURCE)
    partial = load_json(PARTIAL_PREFIX_SOURCE)
    endpoint_completion = load_json(ENDPOINT_COMPLETION_SOURCE)
    lambda0_prefix = load_json(LAMBDA0_PREFIX_SOURCE)
    if power14.get("summary", {}).get("full_tail_power") != 14:
        raise RuntimeError("power-fourteen dominance source changed")
    if power14.get("summary", {}).get("tail_start_k") != POWER_START:
        raise RuntimeError("power-fourteen dominance start changed")
    if power14.get("summary", {}).get("positive_analytic_gates") != 14:
        raise RuntimeError("power-fourteen dominance source is not closed")
    if power14.get("diagnostics", {}).get("full_tail_relative_bound") != (
        "0<=delta_k=(M_k-M_k^(1))/M_k^(1)<2/k^14 for every integer k>=380"
    ):
        raise RuntimeError("power-fourteen dominance contract changed")
    if curvature.get("theorem") != Y_CONTINUOUS:
        raise RuntimeError("order-eleven curvature theorem changed")
    if curvature.get("summary", {}).get("global_first_summand_curvature_theorems") != 1:
        raise RuntimeError("order-eleven curvature theorem is not closed")
    exact11 = entry11.get("exact", {})
    if exact11.get("global_endpoint") != "Q_(11,n)(-100)>0 for every integer n>=0":
        raise RuntimeError("order-eleven endpoint theorem changed")
    if exact11.get("first_discrete_ceiling") != (
        "y_1''(t)<=6000/t^2 => Y_k^(1)<=6000*[-log(1-1/k^2)]<6001/k^2, k>=1253"
    ):
        raise RuntimeError("order-eleven first discrete ceiling changed")
    if exact11.get("full_curvature_ceiling") != (
        "[z_1''(t)<=4200/t^2 on t>=1251 and y_1''(t)<=6000/t^2 on t>=1252] => "
        "Y_k<6001/k^2+37/k^2=6038/k^2<6100/k^2, k>=1253"
    ):
        raise RuntimeError("order-eleven full curvature ceiling changed")
    if entry10.get("exact", {}).get("delayed_entry") != (
        "Q_(10,n)(-100)>0 for every integer n>=4"
    ):
        raise RuntimeError("order-ten positive endpoint ray changed")
    if completion11.get("exact", {}).get("delayed_order11_heat_ray") != ORDER11_HEAT_RAY:
        raise RuntimeError("order-eleven heat ray changed")
    if delayed.get("exact", {}).get("shifted_single_layer_implication") != (
        "[Q_(m-1,n)(lambda)>0 for every n>=n0 on -100<=lambda<=0, the fixed-order m eventual tail holds, and Q_(m,n)(-100)>0 for every n>=n0] => [Q_(m,n)(lambda)>0 for every n>=n0 and -100<=lambda<=0]"
    ):
        raise RuntimeError("generic delayed heat theorem changed")
    if defect.get("tail_arithmetic", {}).get("defect_buffer") != (
        "-3*log(x_k)>3*d_k>=753/(250*(2*k+1))"
    ):
        raise RuntimeError("coefficient defect source changed")
    if partial.get("exact", {}).get("negative_prefix") != (
        "Q_(12,n)(-100)<0 for n=0,1,2,3"
    ):
        raise RuntimeError("order-twelve negative endpoint prefix changed")
    if partial.get("exact", {}).get("positive_block") != (
        "Q_(12,n)(-100)>0 for every 4<=n<=1240"
    ):
        raise RuntimeError("order-twelve positive endpoint block changed")
    if partial.get("summary", {}).get("inconclusive_Q12_rows") != 0:
        raise RuntimeError("order-twelve partial prefix is inconclusive")
    if endpoint_completion.get("exact", {}).get("collar_theorem") != (
        ORDER12_FINITE_PREFIX_THEOREM
    ):
        raise RuntimeError("order-twelve endpoint collar theorem changed")
    if endpoint_completion.get("exact", {}).get("finite_sign_chart") != (
        "Q_(12,n)(-100)<0 for n=0,1,2,3 and "
        "Q_(12,n)(-100)>0 for every 4<=n<=1492"
    ):
        raise RuntimeError("order-twelve finite endpoint sign chart changed")
    if (
        endpoint_completion.get("summary", {}).get("positive_endpoint_collar_rows")
        != 252
        or endpoint_completion.get("summary", {}).get(
            "inconclusive_endpoint_collar_rows"
        )
        != 0
    ):
        raise RuntimeError("order-twelve endpoint collar is not complete")
    if lambda0_prefix.get("finite", {}).get("theorem") != ORDER12_LAMBDA0_PREFIX_TARGET:
        raise RuntimeError("order-twelve lambda-zero prefix changed")
    if lambda0_prefix.get("summary", {}).get("positive_Q12_rows") != 4:
        raise RuntimeError("order-twelve lambda-zero prefix is not complete")
    payloads = (
        power14,
        curvature,
        entry11,
        entry10,
        completion11,
        delayed,
        defect,
        partial,
        endpoint_completion,
        lambda0_prefix,
    )
    paths = (
        POWER14_SOURCE,
        ORDER11_CURVATURE_SOURCE,
        ORDER11_ENTRY_SOURCE,
        ORDER10_ENTRY_SOURCE,
        ORDER11_COMPLETION_SOURCE,
        DELAYED_HEAT_SOURCE,
        DEFECT_SOURCE,
        PARTIAL_PREFIX_SOURCE,
        ENDPOINT_COMPLETION_SOURCE,
        LAMBDA0_PREFIX_SOURCE,
    )
    return [source_record(path, payload) for path, payload in zip(paths, payloads)]


def rational_power_envelope() -> dict:
    sf = bridge_core.stencil_factor
    row = bridge_core.envelope_row
    a = 2 * (Fraction(WALL_START, POWER_START) ** 14 + 3)
    ell = 4 * a

    j = 382
    first_gap = 2 * a / j + ell * sf(13, j)
    log_first_gap = 2 * ell / j + 8 * first_gap

    j = 383
    second_gap = 3 * a / j**2 + log_first_gap * sf(12, j)
    order5 = 2 * log_first_gap / j + Fraction(5, 7) * second_gap + ell / j**2

    j = 384
    third_gap = 4 * a / j**3 + order5 * sf(11, j)
    order6 = 2 * order5 / j + log_first_gap / j**2 + Fraction(2, 3) * third_gap

    j = 385
    fourth_gap = 5 * a / j**4 + order6 * sf(10, j)
    order7 = 2 * order6 / j + order5 / j**2 + Fraction(2, 3) * fourth_gap

    j = 386
    fifth_gap = 6 * a / j**5 + order7 * sf(9, j)

    j = 1249
    order8 = 2 * order7 / j + order6 / j**2 + Fraction(2, 3) * fifth_gap

    j = 1250
    sixth_gap = 7 * a / j**6 + order8 * sf(8, j)
    order9 = 2 * order8 / j + order7 / j**2 + Fraction(3, 4) * sixth_gap

    j = 1251
    seventh_gap = 8 * a / j**7 + order9 * sf(7, j)
    order10 = 2 * order9 / j + order8 / j**2 + 5 * seventh_gap

    j = 1252
    eighth_gap = 9 * a / j**8 + order10 * sf(6, j)
    order11 = 2 * order10 / j + order9 / j**2 + eighth_gap

    j = GAP_FLOOR_START
    ninth_gap = 10 * a / j**9 + order11 * sf(5, j)
    order12 = 2 * order11 / j + order10 / j**2 + ninth_gap

    k = TAIL_FIRST_K
    transfer_scaled = order12 * sf(4, k) / k**2
    if transfer_scaled >= TRANSFER_CONSTANT:
        raise RuntimeError("order-twelve rational transfer envelope exhausted")

    rows = [
        row("moment wall", "a", 381, 14, a, "a_j=2*((j-1)^(-14)+2*j^(-14)+(j+1)^(-14))"),
        row("log defect", "L", 381, 13, ell, "L_j=4*j*a_j"),
        row("first gap", "U1", 382, 13, first_gap, "U1_j=2*a_j+stencil(L)_j"),
        row("log first gap", "V1", 382, 12, log_first_gap, "V1_j=2*L_j+8*j*U1_j"),
        row("second gap", "W1", 383, 12, second_gap, "W1_j=3*a_j+stencil(V1)_j"),
        row("order-five coordinate", "E", 383, 11, order5, "E_j=2*V1_j+(5*j/7)*W1_j+L_j"),
        row("third gap", "Z", 384, 11, third_gap, "Z_j=4*a_j+stencil(E)_j"),
        row("order-six coordinate", "Y1", 384, 10, order6, "Y1_j=2*E_j+V1_j+(2*j/3)*Z_j"),
        row("fourth gap", "O", 385, 10, fourth_gap, "O_j=5*a_j+stencil(Y1)_j"),
        row("order-seven coordinate", "N", 385, 9, order7, "N_j=2*Y1_j+E_j+(2*j/3)*O_j"),
        row("fifth gap", "P", 386, 9, fifth_gap, "P_j=6*a_j+stencil(N)_j"),
        row("order-eight coordinate", "C", 1249, 8, order8, "C_j=2*N_j+Y1_j+(2*j/3)*P_j"),
        row("sixth gap", "D", 1250, 8, sixth_gap, "D_j=7*a_j+stencil(C)_j"),
        row("order-nine coordinate", "F", 1250, 7, order9, "F_j=2*C_j+N_j+(3*j/4)*D_j"),
        row("seventh gap", "D7", 1251, 7, seventh_gap, "D7_j=8*a_j+stencil(F)_j"),
        row("order-ten coordinate", "G", 1251, 6, order10, "G_j=2*F_j+C_j+5*j*D7_j"),
        row("eighth gap", "D8", 1252, 6, eighth_gap, "D8_j=9*a_j+stencil(G)_j"),
        row("order-eleven coordinate", "H", 1252, 5, order11, "H_j=2*G_j+F_j+j*D8_j"),
        row("ninth gap", "D9", GAP_FLOOR_START, 5, ninth_gap, "D9_j=10*a_j+stencil(H)_j"),
        row("order-twelve coordinate", "J", GAP_FLOOR_START, 4, order12, "J_j=2*H_j+G_j+j*D9_j"),
    ]
    return {
        "stencil_lemma": "E_m<=C/m^p for m>=J-1 implies E_(j-1)+2*E_j+E_(j+1)<=C*((J/(J-1))^p+3)/j^p, j>=J",
        "rows": rows,
        "transfer_scaled_exact": str(transfer_scaled),
        "transfer_scaled_decimal": bridge_core.decimal_text(transfer_scaled),
        "transfer_reserve_exact": str(Fraction(TRANSFER_CONSTANT) - transfer_scaled),
        "transfer_reserve_decimal": bridge_core.decimal_text(Fraction(TRANSFER_CONSTANT) - transfer_scaled),
        "transfer_bound": V_FULL_TRANSFER,
    }


def exact_diagnostics() -> dict:
    j = sp.symbols("j", integer=True, positive=True)
    first_floor = bridge_core.shifted_positive_polynomial(
        10 / (2 * j + 1) - sp.Rational(6001) / j**2 - 1 / j,
        j,
        GAP_FLOOR_START,
    )
    full_floor = bridge_core.shifted_positive_polynomial(
        sp.Rational(2510, 250) / (2 * j + 1) - sp.Rational(6038) / j**2 - 1 / j,
        j,
        GAP_FLOOR_START,
    )
    k, m = sp.symbols("k m", integer=True, nonnegative=True)
    comparison = sp.expand(2761 * k**2 - 250 * SAFE_CEILING * (2 * k + 1))
    shifted = sp.expand(comparison.subs(k, TAIL_FIRST_K + m))
    coefficients = sp.Poly(shifted, m).all_coeffs()
    if any(value <= 0 for value in coefficients):
        raise RuntimeError("order-twelve endpoint comparison is not coefficient-positive")

    a, y, z, gap = sp.symbols("A y z U", positive=True)
    q10_center = a**10 * sp.exp(y)
    q9_denominator = a**9 * sp.exp(z)
    q11 = sp.factor(q10_center**2 * (1 - sp.exp(-gap)) / q9_denominator)
    target = a**11 * sp.exp(2 * y - z) * (1 - sp.exp(-gap))
    if sp.simplify(q11 - target) != 0:
        raise RuntimeError("canonical order-eleven factorization failed")

    return {
        "ninth_gap": "U(t)=10*B(t)-y(t-1)+2*y(t)-y(t+1)",
        "order11_coordinate": "v(t)=2*y(t)-z(t)+log(1-exp(-U(t)))",
        "canonical_factorization": "Q_(11,n)=A_(n+10)^11*exp(v(n+10))",
        "canonical_factorization_residual": "0",
        "curvature_identity": "F_n=log(Q_(11,n)*Q_(11,n+2)/Q_(11,n+1)^2)=11*log(x_k)+V_k, V_k=v(k-1)-2*v(k)+v(k+1), k=n+11",
        "sign_equivalence": "Q_(12,n)>0 iff F_n<0, provided Q_(10,n+2)>0 and Q_(11,n),Q_(11,n+1),Q_(11,n+2)>0",
        "first_U_floor": "U_j^(1)>=10/(2*j+1)-6001/j^2>1/j, j>=1503",
        "full_U_floor": "U_j>=2510/(250*(2*j+1))-6038/j^2>1/j, j>=1503",
        "U_floor": "min(U_j,U_j^(1))>1/j, j>=1503",
        "stable_log_lipschitz": "phi(min(U_j,U_j^(1)))<j",
        "first_floor_polynomial": first_floor,
        "full_floor_polynomial": full_floor,
        "power_envelope": rational_power_envelope(),
        "continuous_target": V_CONTINUOUS_TARGET,
        "tent_transfer": V_FIRST_DISCRETE,
        "full_transfer": V_FULL_TRANSFER,
        "conditional_full_ceiling": V_FULL_CEILING,
        "defect_buffer": "-11*log(x_k)>=11*d_k>=2761/(250*(2*k+1)), k>=320",
        "sufficient_ceiling": "V_k<=8200/k^2 for every integer k>=1504",
        "rational_comparison": "8200/k^2<2761/(250*(2*k+1)), k>=1504",
        "cleared_polynomial": str(comparison),
        "shifted_polynomial_k_1504_plus_m": str(shifted),
        "shifted_coefficients": [str(value) for value in coefficients],
        "conditional_endpoint_tail": "v_1''(t)<=8000/t^2 on t>=1503 => Q_(12,n)(-100)>0 for every n>=1493",
        "certified_endpoint_sign_chart": "Q_(12,n)(-100)<0 for n=0,1,2,3 and Q_(12,n)(-100)>0 for every 4<=n<=1240",
        "finite_prefix_theorem": ORDER12_FINITE_PREFIX_THEOREM,
        "delayed_heat_theorem": ORDER12_DELAYED_HANDOFF,
        "lambda0_prefix_target": ORDER12_LAMBDA0_PREFIX_TARGET,
    }


def build_artifact() -> dict:
    sources = validate_sources()
    exact = exact_diagnostics()
    envelope = exact["power_envelope"]
    rows = [
        TargetRow("co12cbt_01_power14", "theorem_input", "ready_to_apply", "Power-fourteen dominance supplies the perturbative reserve for one more stable logarithm.", "delta_k<2/k^14, k>=380", "Full versus first Newman moment wall at lambda=-100."),
        TargetRow("co12cbt_02_coordinate", "exact_identity", "ready_to_apply", "The next signed-condensation stage gives a canonical positive coordinate for Q11 on its endpoint ray.", exact["ninth_gap"] + "; " + exact["order11_coordinate"] + "; " + exact["canonical_factorization"], "Conditional only on the inherited positive Q10 and Q11 tail."),
        TargetRow("co12cbt_03_floor", "exact_bound", "ready_to_apply", "The completed order-eleven curvature theorem gives a common inverse-linear U floor.", exact["first_U_floor"] + "; " + exact["full_U_floor"], "Uses the completed first and full order-eleven discrete curvature ceilings."),
        TargetRow("co12cbt_04_transfer", "conditional_transfer_theorem", "ready_to_apply", "The exact twenty-row rational envelope keeps full and first ninth-nested curvature within 100/k^2.", exact["full_transfer"], "Uses the common U floor and power-fourteen wall.", {"scaled_transfer": envelope["transfer_scaled_decimal"]}),
        TargetRow("co12cbt_05_first_target", "analytic_theorem_target", "not_ready_to_apply", "Prove the continuous ninth-nested first-summand curvature ceiling.", exact["continuous_target"], "New continuous half-line theorem; no interval certificate exists yet."),
        TargetRow("co12cbt_06_full_ceiling", "conditional_theorem", "ready_to_apply", "The first-summand target and transfer fit below the explicit full-kernel budget.", exact["conditional_full_ceiling"], "Conditional on the displayed continuous premise."),
        TargetRow("co12cbt_07_endpoint_tail", "conditional_theorem", "ready_to_apply", "The 8200/k^2 budget lies strictly inside the eleventh coefficient-defect buffer.", exact["rational_comparison"] + "; " + exact["conditional_endpoint_tail"], "Coefficient-positive rational comparison on k>=1504."),
        TargetRow("co12cbt_08_partial_prefix", "interval_theorem", "ready_to_apply", "The inherited endpoint source proves the delayed sign chart through n=1240.", exact["certified_endpoint_sign_chart"], "The four negative endpoint shifts are preserved."),
        TargetRow("co12cbt_09_finite_prefix", "interval_theorem", "ready_to_apply", "The enlarged coefficient collar closes the positive endpoint block up to the conditional tail.", exact["finite_prefix_theorem"], "The 252 shifts n=1241,...,1492 are rigorously positive; no analytic tail is inferred."),
        TargetRow("co12cbt_10_heat_handoff", "conditional_forward_theorem", "ready_to_apply", "Generic shifted cooperative descent propagates any completed n>=4 order-twelve endpoint ray.", exact["delayed_heat_theorem"], "The theorem is ready; its endpoint premise now depends only on the continuum target and its conditional tail."),
        TargetRow("co12cbt_11_lambda0_prefix", "finite_theorem_input", "ready_to_apply", "Fresh 520-digit determinants supply the four lambda-zero shifts outside delayed heat descent.", exact["lambda0_prefix_target"], "Lambda zero and shifts zero through three only."),
    ]
    return {
        "kind": "jensen_window_pf_compound_order12_curvature_bridge_target",
        "date": "2026-07-22",
        "status": "exact conditional order-twelve curvature bridge with one open continuum target",
        "proof_boundary": (
            "This derives the order-twelve canonical coordinate, common gap floor, "
            "power-fourteen transfer envelope, conditional endpoint-tail arithmetic, "
            "completed finite endpoint collar, and delayed heat handoff. It does not prove "
            "the 8000/t^2 continuum target, analytic endpoint tail, all-shift order twelve, "
            "PF-infinity, RH, or Lambda<=0."
        ),
        "sources": sources,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": 10,
            "open_rows": 1,
            "exact_factorizations": 1,
            "power_envelope_rows": len(envelope["rows"]),
            "conditional_transfer_theorems": 1,
            "conditional_endpoint_tail_theorems": 1,
            "open_continuum_targets": 1,
            "open_finite_endpoint_targets": 0,
            "lambda0_prefix_theorems": 1,
            "open_lambda0_prefix_targets": 0,
            "conditional_heat_handoffs": 1,
            "orders_above_12": 0,
            "rh_claims": 0,
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order12_curvature_bridge_target.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order12_curvature_bridge_target.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    envelope = exact["power_envelope"]
    lines = [
        "# Order-Twelve Curvature Bridge Target",
        "",
        "Date: 2026-07-22",
        "",
        "Status: exact conditional reduction with one open continuum target.",
        "This is not a proof of order twelve, PF-infinity, RH, or `Lambda<=0`.",
        "",
        "## Coordinate",
        "",
        "```text",
        exact["ninth_gap"],
        exact["order11_coordinate"],
        exact["canonical_factorization"],
        exact["curvature_identity"],
        "```",
        "",
        "## Transfer",
        "",
        "```text",
        exact["U_floor"],
        exact["stable_log_lipschitz"],
        exact["full_transfer"],
        f"exact scaled transfer={envelope['transfer_scaled_decimal']}<100",
        "```",
        "",
        "The power envelope has twenty exact rational rows.",
        "",
        "## Endpoint Budget",
        "",
        "```text",
        exact["continuous_target"],
        exact["tent_transfer"],
        exact["conditional_full_ceiling"],
        exact["defect_buffer"],
        exact["rational_comparison"],
        exact["conditional_endpoint_tail"],
        exact["certified_endpoint_sign_chart"],
        exact["finite_prefix_theorem"],
        "```",
        "",
        "## Open Targets",
        "",
        "```text",
        exact["continuous_target"],
        "```",
        "",
        "The power-fourteen input, ninth-gap floor, exact transfer, endpoint-tail",
        "implication, completed finite endpoint collar, lambda-zero prefix, and delayed",
        "heat handoff are ready. The single target displayed under `Open Targets` remains",
        "open. Separately, " + exact["lambda0_prefix_target"] + " is proved and is",
        "not an open target. No order above twelve, PF-infinity, RH, or `Lambda<=0`",
        "conclusion is claimed.",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    envelope = artifact["exact"]["power_envelope"]
    print(
        "wrote order-twelve curvature bridge target: 20 envelope rows, "
        f"scaled transfer {envelope['transfer_scaled_decimal']}<100"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
