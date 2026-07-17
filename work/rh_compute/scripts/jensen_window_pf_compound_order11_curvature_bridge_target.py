#!/usr/bin/env python3
"""Derive the conditional order-eleven curvature and transfer target."""

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

import jensen_window_pf_compound_order10_first_summand_curvature_bridge as order10_bridge  # noqa: E402


RESULTS = REPO_ROOT / "work/rh_compute/results"
POWER13_SOURCE = RESULTS / "jensen_window_pf_negative_lambda_first_summand_power13_rebalanced_dominance_extension.json"
ORDER10_BRIDGE_SOURCE = RESULTS / "jensen_window_pf_compound_order10_first_summand_curvature_bridge.json"
ORDER10_TAIL_SOURCE = RESULTS / "jensen_window_pf_compound_order10_m100_tail_curvature_reduction.json"
DEFECT_SOURCE = RESULTS / "jensen_window_pf_compound_order4_m100_entry_certificate.json"
POINT_SCOUT_SOURCE = RESULTS / "jensen_window_pf_compound_order11_first_summand_point_scout.json"
ORDER11_PREFIX_SOURCE = RESULTS / "jensen_window_pf_compound_order11_m100_prefix_certificate.json"
DELAYED_HEAT_SOURCE = RESULTS / "jensen_window_pf_delayed_cooperative_heat_tail_lemma.json"
ORDER10_COMPLETION_SOURCE = RESULTS / "jensen_window_pf_compound_order10_lambda0_completion_certificate.json"
DEFAULT_OUT = RESULTS / "jensen_window_pf_compound_order11_curvature_bridge_target.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_compound_order11_curvature_bridge_target.md"

POWER_START = 340
B_ERROR_START = 341
X_FLOOR_START = 1252
TAIL_FIRST_K = 1253
TAIL_FIRST_N = 1243
INHERITED_Z_CONTINUOUS = "z_1''(t)<=4200/t^2 for every real t>=1251"
Y_CONTINUOUS_TARGET = "y_1''(t)<=6000/t^2 for every real t>=1252"
Y_FIRST_DISCRETE = "y_1''(t)<=6000/t^2 => Y_k^(1)<=6000*[-log(1-1/k^2)]<6001/k^2, k>=1253"
Y_FULL_TRANSFER = "|Y_k-Y_k^(1)|<37/k^2 for every integer k>=1253"
Y_FULL_CEILING = "[z_1''(t)<=4200/t^2 on t>=1251 and y_1''(t)<=6000/t^2 on t>=1252] => Y_k<6001/k^2+37/k^2=6038/k^2<6100/k^2, k>=1253"
ORDER11_FINITE_PREFIX = "Q_(11,n)(-100)>0 for every 0<=n<=1242"
ORDER10_DELAYED_HEAT_RAY = "Q_(10,n)(lambda)>0 for every n>=4 and -100<=lambda<=0"
ORDER11_DELAYED_HANDOFF = (
    "[Q_(10,n)(lambda)>0 for every n>=4 and -100<=lambda<=0 and "
    "Q_(11,n)(-100)>0 for every n>=4] implies "
    "Q_(11,n)(lambda)>0 for every n>=4 and -100<=lambda<=0"
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
    power13 = load_json(POWER13_SOURCE)
    bridge = load_json(ORDER10_BRIDGE_SOURCE)
    tail = load_json(ORDER10_TAIL_SOURCE)
    defect = load_json(DEFECT_SOURCE)
    scout = load_json(POINT_SCOUT_SOURCE)
    prefix = load_json(ORDER11_PREFIX_SOURCE)
    delayed = load_json(DELAYED_HEAT_SOURCE)
    completion = load_json(ORDER10_COMPLETION_SOURCE)
    if power13.get("diagnostics", {}).get("full_tail_relative_bound") != (
        "0<=delta_k=(M_k-M_k^(1))/M_k^(1)<2/k^13 for every integer k>=340"
    ):
        raise RuntimeError("power-thirteen dominance source changed")
    if power13.get("summary", {}).get("positive_analytic_gates") != 14:
        raise RuntimeError("power-thirteen dominance source is not closed")
    exact_bridge = bridge.get("exact", {})
    if exact_bridge.get("continuous_target") != INHERITED_Z_CONTINUOUS:
        raise RuntimeError("order-ten continuous target changed")
    if exact_bridge.get("full_transfer") != "|Z_k-Z_k^(1)|<10/k^2 for every integer k>=1252":
        raise RuntimeError("order-ten transfer source changed")
    if exact_bridge.get("conditional_full_ceiling") != (
        "z_1''(t)<=4200/t^2 on t>=1251 => "
        "Z_k<4201/k^2+10/k^2=4211/k^2<5500/k^2, k>=1252"
    ):
        raise RuntimeError("order-ten full ceiling changed")
    if tail.get("exact", {}).get("canonical_factorization") != "Q_(9,n)=A_(n+8)^9*exp(z(n+8))":
        raise RuntimeError("order-ten canonical coordinate changed")
    if defect.get("tail_arithmetic", {}).get("defect_buffer") != "-3*log(x_k)>3*d_k>=753/(250*(2*k+1))":
        raise RuntimeError("coefficient defect source changed")
    if scout.get("exact", {}).get("eighth_gap") != "X=9*B-Delta^2*z":
        raise RuntimeError("order-eleven point scout coordinate changed")
    if scout.get("summary", {}).get("point_curvature_enclosures") != 8:
        raise RuntimeError("order-eleven point scout is not closed")
    if prefix.get("exact", {}).get("finite_prefix") != ORDER11_FINITE_PREFIX:
        raise RuntimeError("order-eleven finite endpoint prefix changed")
    if prefix.get("summary", {}).get("finite_prefix_theorems") != 1:
        raise RuntimeError("order-eleven finite endpoint prefix is not closed")
    if delayed.get("exact", {}).get("order11_delayed_handoff") != ORDER11_DELAYED_HANDOFF:
        raise RuntimeError("order-eleven delayed heat theorem changed")
    if delayed.get("summary", {}).get("order11_n4_specializations") != 1:
        raise RuntimeError("order-eleven delayed heat theorem is not closed")
    if completion.get("exact", {}).get("delayed_heat_ray") != ORDER10_DELAYED_HEAT_RAY:
        raise RuntimeError("order-ten delayed heat ray changed")
    payloads = (power13, bridge, tail, defect, scout, prefix, delayed, completion)
    paths = (
        POWER13_SOURCE,
        ORDER10_BRIDGE_SOURCE,
        ORDER10_TAIL_SOURCE,
        DEFECT_SOURCE,
        POINT_SCOUT_SOURCE,
        ORDER11_PREFIX_SOURCE,
        DELAYED_HEAT_SOURCE,
        ORDER10_COMPLETION_SOURCE,
    )
    return [source_record(path, payload) for path, payload in zip(paths, payloads)]


def rational_power_envelope() -> dict:
    sf = order10_bridge.stencil_factor
    row = order10_bridge.envelope_row
    a = 2 * (Fraction(341, 340) ** 13 + 3)
    ell = 4 * a

    j = 342
    first_gap = 2 * a / j + ell * sf(12, j)
    log_first_gap = 2 * ell / j + 8 * first_gap

    j = 343
    second_gap = 3 * a / j**2 + log_first_gap * sf(11, j)
    order5 = 2 * log_first_gap / j + Fraction(5, 7) * second_gap + ell / j**2

    j = 344
    third_gap = 4 * a / j**3 + order5 * sf(10, j)
    order6 = 2 * order5 / j + log_first_gap / j**2 + Fraction(2, 3) * third_gap

    j = 345
    fourth_gap = 5 * a / j**4 + order6 * sf(9, j)
    order7 = 2 * order6 / j + order5 / j**2 + Fraction(2, 3) * fourth_gap

    j = 346
    fifth_gap = 6 * a / j**5 + order7 * sf(8, j)

    j = 1249
    order8 = 2 * order7 / j + order6 / j**2 + Fraction(2, 3) * fifth_gap

    j = 1250
    sixth_gap = 7 * a / j**6 + order8 * sf(7, j)
    order9 = 2 * order8 / j + order7 / j**2 + Fraction(3, 4) * sixth_gap

    j = 1251
    seventh_gap = 8 * a / j**7 + order9 * sf(6, j)
    order10 = 2 * order9 / j + order8 / j**2 + 5 * seventh_gap

    j = X_FLOOR_START
    eighth_gap = 9 * a / j**8 + order10 * sf(5, j)
    order11 = 2 * order10 / j + order9 / j**2 + eighth_gap

    k = TAIL_FIRST_K
    transfer_scaled = order11 * sf(4, k) / k**2
    if transfer_scaled >= 37:
        raise RuntimeError("order-eleven rational transfer envelope exhausted")

    rows = [
        row("moment wall", "a", 341, 13, a, "a_j=2*((j-1)^(-13)+2*j^(-13)+(j+1)^(-13))"),
        row("log defect", "L", 341, 12, ell, "L_j=4*j*a_j"),
        row("first gap", "U1", 342, 12, first_gap, "U1_j=2*a_j+stencil(L)_j"),
        row("log first gap", "V1", 342, 11, log_first_gap, "V1_j=2*L_j+8*j*U1_j"),
        row("second gap", "W1", 343, 11, second_gap, "W1_j=3*a_j+stencil(V1)_j"),
        row("order-five coordinate", "E", 343, 10, order5, "E_j=2*V1_j+(5*j/7)*W1_j+L_j"),
        row("third gap", "Z", 344, 10, third_gap, "Z_j=4*a_j+stencil(E)_j"),
        row("order-six coordinate", "Y1", 344, 9, order6, "Y1_j=2*E_j+V1_j+(2*j/3)*Z_j"),
        row("fourth gap", "O", 345, 9, fourth_gap, "O_j=5*a_j+stencil(Y1)_j"),
        row("order-seven coordinate", "N", 345, 8, order7, "N_j=2*Y1_j+E_j+(2*j/3)*O_j"),
        row("fifth gap", "P", 346, 8, fifth_gap, "P_j=6*a_j+stencil(N)_j"),
        row("order-eight coordinate", "C", 1249, 7, order8, "C_j=2*N_j+Y1_j+(2*j/3)*P_j"),
        row("sixth gap", "D", 1250, 7, sixth_gap, "D_j=7*a_j+stencil(C)_j"),
        row("order-nine coordinate", "F", 1250, 6, order9, "F_j=2*C_j+N_j+(3*j/4)*D_j"),
        row("seventh gap", "D7", 1251, 6, seventh_gap, "D7_j=8*a_j+stencil(F)_j"),
        row("order-ten coordinate", "G", 1251, 5, order10, "G_j=2*F_j+C_j+5*j*D7_j"),
        row("eighth gap", "D8", 1252, 5, eighth_gap, "D8_j=9*a_j+stencil(G)_j"),
        row("order-eleven coordinate", "H", 1252, 4, order11, "H_j=2*G_j+F_j+j*D8_j"),
    ]
    return {
        "stencil_lemma": "E_m<=C/m^p for m>=J-1 implies E_(j-1)+2*E_j+E_(j+1)<=C*((J/(J-1))^p+3)/j^p, j>=J",
        "rows": rows,
        "transfer_scaled_exact": str(transfer_scaled),
        "transfer_scaled_decimal": order10_bridge.decimal_text(transfer_scaled),
        "transfer_reserve_exact": str(Fraction(37) - transfer_scaled),
        "transfer_reserve_decimal": order10_bridge.decimal_text(Fraction(37) - transfer_scaled),
        "transfer_bound": Y_FULL_TRANSFER,
    }


def exact_diagnostics() -> dict:
    j = sp.symbols("j", integer=True, positive=True)
    first_floor = order10_bridge.shifted_positive_polynomial(
        9 / (2 * j + 1) - sp.Rational(4201) / j**2 - 1 / j,
        j,
        X_FLOOR_START,
    )
    full_floor = order10_bridge.shifted_positive_polynomial(
        sp.Rational(2259, 250) / (2 * j + 1) - sp.Rational(4211) / j**2 - 1 / j,
        j,
        X_FLOOR_START,
    )
    k, m = sp.symbols("k m", integer=True, nonnegative=True)
    comparison = sp.expand(2510 * k**2 - 250 * 6100 * (2 * k + 1))
    shifted = sp.expand(comparison.subs(k, TAIL_FIRST_K + m))
    coefficients = sp.Poly(shifted, m).all_coeffs()
    if any(value <= 0 for value in coefficients):
        raise RuntimeError("order-eleven endpoint comparison is not coefficient-positive")

    a, z, w, gap = sp.symbols("A z w X", positive=True)
    q9_center = a**9 * sp.exp(z)
    q8_denominator = a**8 * sp.exp(w)
    q10 = sp.factor(q9_center**2 * (1 - sp.exp(-gap)) / q8_denominator)
    target = a**10 * sp.exp(2 * z - w) * (1 - sp.exp(-gap))
    if sp.simplify(q10 - target) != 0:
        raise RuntimeError("canonical order-ten factorization failed")

    return {
        "eighth_gap": "X(t)=9*B(t)-z(t-1)+2*z(t)-z(t+1)",
        "order10_coordinate": "y(t)=2*z(t)-w(t)+log(1-exp(-X(t)))",
        "canonical_factorization": "Q_(10,n)=A_(n+9)^10*exp(y(n+9))",
        "canonical_factorization_residual": "0",
        "curvature_identity": "E_n=log(Q_(10,n)*Q_(10,n+2)/Q_(10,n+1)^2)=10*log(x_k)+Y_k, Y_k=y(k-1)-2*y(k)+y(k+1), k=n+10",
        "sign_equivalence": "Q_(11,n)>0 iff E_n<0, provided Q_(9,n+2)>0 and Q_(10,n),Q_(10,n+1),Q_(10,n+2)>0",
        "first_X_floor": "X_j^(1)>=9/(2*j+1)-4201/j^2>1/j, j>=1252",
        "full_X_floor": "X_j>=2259/(250*(2*j+1))-4211/j^2>1/j, j>=1252",
        "X_floor": "min(X_j,X_j^(1))>1/j, j>=1252",
        "stable_log_lipschitz": "phi(min(X_j,X_j^(1)))<j",
        "first_floor_polynomial": first_floor,
        "full_floor_polynomial": full_floor,
        "power_envelope": rational_power_envelope(),
        "continuous_target": Y_CONTINUOUS_TARGET,
        "tent_transfer": Y_FIRST_DISCRETE,
        "full_transfer": Y_FULL_TRANSFER,
        "conditional_full_ceiling": Y_FULL_CEILING,
        "defect_buffer": "-10*log(x_k)>=10*d_k>=2510/(250*(2*k+1)), k>=320",
        "sufficient_ceiling": "Y_k<=6100/k^2 for every integer k>=1253",
        "rational_comparison": "6100/k^2<2510/(250*(2*k+1)), k>=1253",
        "cleared_polynomial": str(comparison),
        "shifted_polynomial_k_1253_plus_m": str(shifted),
        "shifted_coefficients": [str(value) for value in coefficients],
        "conditional_endpoint_tail": "[z_1''(t)<=4200/t^2 on t>=1251 and y_1''(t)<=6000/t^2 on t>=1252] => Q_(11,n)(-100)>0 for every n>=1243",
        "finite_prefix_theorem": ORDER11_FINITE_PREFIX,
        "delayed_heat_theorem": ORDER11_DELAYED_HANDOFF,
    }


def build_artifact() -> dict:
    sources = validate_sources()
    exact = exact_diagnostics()
    envelope = exact["power_envelope"]
    rows = [
        TargetRow("co11cbt_01_power13", "theorem_input", "ready_to_apply", "Power-thirteen dominance supplies the perturbative reserve for one more stable logarithm.", "delta_k<2/k^13, k>=340", "Full versus first Newman moment wall at lambda=-100."),
        TargetRow("co11cbt_02_coordinate", "exact_identity", "ready_to_apply", "The next signed-condensation stage gives a canonical positive coordinate for Q10 on its positive ray.", exact["eighth_gap"] + "; " + exact["order10_coordinate"] + "; " + exact["canonical_factorization"], "Conditional on the inherited positive order-ten ray."),
        TargetRow("co11cbt_03_floor", "conditional_exact_bound", "ready_to_apply", "The inherited order-ten curvature premise gives a common inverse-linear X floor.", exact["first_X_floor"] + "; " + exact["full_X_floor"], "Conditional on the inherited z curvature theorem."),
        TargetRow("co11cbt_04_transfer", "conditional_transfer_theorem", "ready_to_apply", "The exact eighteen-row rational envelope keeps full and first eighth-nested curvature within 37/k^2.", exact["full_transfer"], "Uses the common X floor and power-thirteen wall.", {"scaled_transfer": envelope["transfer_scaled_decimal"]}),
        TargetRow("co11cbt_05_first_target", "analytic_theorem_target", "not_ready_to_apply", "Prove the continuous eighth-nested first-summand curvature ceiling.", exact["continuous_target"], "New localized continuum theorem; selected-point evidence is not enough."),
        TargetRow("co11cbt_06_full_ceiling", "conditional_theorem", "ready_to_apply", "The first-summand target and transfer fit below the explicit full-kernel budget.", exact["conditional_full_ceiling"], "Conditional on both displayed continuous premises."),
        TargetRow("co11cbt_07_endpoint_tail", "conditional_theorem", "ready_to_apply", "The 6100/k^2 budget lies strictly inside the tenth coefficient-defect buffer.", exact["rational_comparison"] + "; " + exact["conditional_endpoint_tail"], "Coefficient-positive rational comparison on k>=1253."),
        TargetRow("co11cbt_08_finite_prefix", "finite_theorem_input", "ready_to_apply", "The rigorous endpoint prefix meets the conditional analytic tail exactly at n=1243.", exact["finite_prefix_theorem"], "Finite endpoint theorem only; the analytic tail remains conditional on the continuum target."),
        TargetRow("co11cbt_09_heat_handoff", "conditional_forward_theorem", "ready_to_apply", "Shifted cooperative descent propagates any completed n>=4 order-eleven endpoint ray.", exact["delayed_heat_theorem"], "The theorem is ready; its order-eleven endpoint-ray premise still depends on the continuum target."),
    ]
    return {
        "kind": "jensen_window_pf_compound_order11_curvature_bridge_target",
        "date": "2026-07-16",
        "status": "exact conditional order-eleven curvature bridge with one open continuum target",
        "proof_boundary": (
            "This derives the order-eleven canonical coordinate, exact power-thirteen "
            "transfer envelope, finite endpoint prefix, delayed heat theorem, and "
            "conditional endpoint-tail arithmetic. It does not prove the 6000/t^2 "
            "continuum target, all-shift order eleven, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": sources,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": 8,
            "open_rows": 1,
            "exact_factorizations": 1,
            "power_envelope_rows": len(envelope["rows"]),
            "conditional_transfer_theorems": 1,
            "conditional_endpoint_tail_theorems": 1,
            "open_continuum_targets": 1,
            "finite_prefix_theorems": 1,
            "conditional_heat_handoffs": 1,
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order11_curvature_bridge_target.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order11_curvature_bridge_target.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    envelope = exact["power_envelope"]
    lines = [
        "# Order-Eleven Curvature Bridge Target",
        "",
        "Date: 2026-07-16",
        "",
        "Status: exact conditional reduction with one continuum target still open. This is not a proof of RH or `Lambda <= 0`.",
        "",
        "## Coordinate",
        "",
        "```text",
        exact["eighth_gap"],
        exact["order10_coordinate"],
        exact["canonical_factorization"],
        exact["curvature_identity"],
        "```",
        "",
        "## Transfer",
        "",
        "Under the inherited order-ten `z` curvature theorem,",
        "",
        "```text",
        exact["X_floor"],
        exact["stable_log_lipschitz"],
        exact["full_transfer"],
        f"exact scaled transfer={envelope['transfer_scaled_decimal']}<37",
        "```",
        "",
        "The power envelope has eighteen exact rational rows.",
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
        exact["finite_prefix_theorem"],
        exact["delayed_heat_theorem"],
        "```",
        "",
        "## Open Work",
        "",
        "The finite endpoint prefix through `n=1242`, the delayed order-ten heat",
        "ray, and the shifted order-eleven heat theorem are now rigorous inputs.",
        "Only the continuum target `y_1''(t)<=6000/t^2` remains open in this",
        "composition. This is not all-shift order eleven,",
        "PF-infinity, RH, or `Lambda<=0`.",
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
    print(f"wrote order-eleven curvature bridge target: 18 envelope rows, scaled transfer {envelope['transfer_scaled_decimal']}<37")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
