#!/usr/bin/env python3
"""Build the exact Gaussian-cumulant target for the remaining order-four ray."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
import math
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402
import sympy as sp  # noqa: E402

from jensen_window_pf_compound_order4_localized_curvature_interval_core import (  # noqa: E402
    DEFAULT_PRECISION_BITS,
    evaluate_localized_curvature_from_h_cover,
    signed_hurwitz_gamma_derivative,
)
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_interval,
    arb_lower_text,
    arb_rational,
    arb_upper_text,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md"
)
FORMAL_FINITE_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.json"
)
FORMAL_RAY_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.json"
)
MAX_EPSILON_ORDER = 6
CORRIDOR_K2 = (Fraction(2, 5), Fraction(4, 5))
CORRIDOR_CAPS = {
    3: Fraction(6, 5),
    4: Fraction(27, 20),
    5: Fraction(3, 2),
    6: Fraction(17, 10),
    7: Fraction(2, 1),
    8: Fraction(5, 2),
}
SAMPLE_COLLARS = (
    (Fraction(2), Fraction(1, 50_000)),
    (Fraction(5, 2), Fraction(1, 100_000)),
    (Fraction(3), Fraction(1, 1_000_000)),
    (Fraction(4), Fraction(1, 100_000_000)),
    (Fraction(5), Fraction(1, 10_000_000_000)),
    (Fraction(10), Fraction(1, 10**19)),
    (Fraction(20), Fraction(1, 10**37)),
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


def exact_interval(lower: flint.arb, upper: flint.arb) -> flint.arb:
    lower_endpoint = flint.arb(lower.lower())
    upper_endpoint = flint.arb(upper.upper())
    return (lower_endpoint + upper_endpoint) / 2 + flint.arb(
        0, (upper_endpoint - lower_endpoint) / 2
    )


def formal_cumulant_expansion() -> dict:
    """Expand log E exp(zY-R_epsilon(Y)) exactly through epsilon^6."""
    epsilon, z, y = sp.symbols("epsilon z y")
    symbols = {order: sp.symbols(f"L_{order}") for order in range(3, 9)}
    perturbation = [sp.Integer(0) for _ in range(MAX_EPSILON_ORDER + 1)]
    for order in range(3, 9):
        perturbation[order - 2] = symbols[order] * y**order / sp.factorial(order)

    exponential = [sp.Integer(0) for _ in range(MAX_EPSILON_ORDER + 1)]
    exponential[0] = sp.Integer(1)
    for degree in range(1, MAX_EPSILON_ORDER + 1):
        exponential[degree] = sp.expand(
            -sum(
                index * perturbation[index] * exponential[degree - index]
                for index in range(1, degree + 1)
            )
            / degree
        )

    maximum_y_degree = max(sp.Poly(value, y).degree() for value in exponential)
    tilted_gaussian_moments = [sp.Integer(1)]
    for _degree in range(maximum_y_degree):
        previous = tilted_gaussian_moments[-1]
        tilted_gaussian_moments.append(sp.expand(z * previous + sp.diff(previous, z)))

    def tilted_expectation(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sum(
                coefficient * tilted_gaussian_moments[monomial[0]]
                for monomial, coefficient in sp.Poly(polynomial, y).terms()
            )
        )

    partition = [tilted_expectation(value) for value in exponential]
    logarithm = [sp.Integer(0) for _ in range(MAX_EPSILON_ORDER + 1)]
    for degree in range(1, MAX_EPSILON_ORDER + 1):
        logarithm[degree] = sp.expand(
            partition[degree]
            - sum(
                index * logarithm[index] * partition[degree - index]
                for index in range(1, degree)
            )
            / degree
        )

    cumulants: dict[str, dict] = {}
    leading_signature: dict[str, dict] = {}
    for order in range(2, 9):
        terms = []
        for degree in range(1, MAX_EPSILON_ORDER + 1):
            coefficient = sp.factor(sp.diff(logarithm[degree], z, order).subs(z, 0))
            if coefficient != 0:
                terms.append(
                    {
                        "epsilon_power": degree,
                        "coefficient": sp.sstr(coefficient),
                    }
                )
        cumulants[str(order)] = {
            "gaussian_base": 1 if order == 2 else 0,
            "terms": terms,
        }
        expected_power = 2 if order == 2 else order - 2
        leading = next(
            term for term in terms if term["epsilon_power"] == expected_power
        )
        coefficient = sp.sympify(leading["coefficient"])
        at_unity = sp.expand(coefficient.subs({value: 1 for value in symbols.values()}))
        expected = sp.Rational(1, 2) if order == 2 else (-1) ** order * math.factorial(order - 2)
        if at_unity != expected:
            raise RuntimeError(f"bad leading cumulant signature at order {order}")
        leading_signature[str(order)] = {
            "epsilon_power": expected_power,
            "coefficient_at_L_equals_1": str(at_unity),
            "expected": str(expected),
        }
    return {
        "grading": "lambda_r=L_r*epsilon^(r-2), epsilon=q^(-1/2)",
        "standardized_potential": (
            "R_epsilon(y)=sum_(r=3)^8 L_r*epsilon^(r-2)*y^r/r!"
        ),
        "partition_recurrence": (
            "n*E_n=-sum_(j=1)^n j*R_j*E_(n-j); "
            "n*S_n=n*Z_n-sum_(j=1)^(n-1) j*S_j*Z_(n-j)"
        ),
        "cumulants": cumulants,
        "leading_signature": leading_signature,
    }


def candidate_h_cover(center: Fraction, collar: Fraction) -> tuple[dict[int, flint.arb], dict]:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    central_u = arb_rational(center)
    central_t = potential_jet_arb(central_u, 1)[1]
    left = center - collar
    right = center + collar
    left_t = potential_jet_arb(arb_rational(left), 1)[1]
    right_t = potential_jet_arb(arb_rational(right), 1)[1]
    if not bool(left_t < central_t - 2 and right_t > central_t + 2):
        raise RuntimeError(f"candidate collar does not cover t+-2 at u={center}")

    mode = arb_interval(left, right)
    jet = potential_jet_arb(mode, 8)
    t, curvature = jet[1], jet[2]
    q = flint.arb.pi() * (4 * mode).exp()

    kappa2 = exact_interval(
        1 + arb_rational(CORRIDOR_K2[0]) / q,
        1 + arb_rational(CORRIDOR_K2[1]) / q,
    )
    derivatives = {
        2: signed_hurwitz_gamma_derivative(2, t + flint.arb(1) / 2)
        - kappa2 / curvature
    }
    for order, cap in CORRIDOR_CAPS.items():
        baseline = math.factorial(order - 2) * q ** (1 - flint.arb(order) / 2)
        magnitude = exact_interval(baseline, arb_rational(cap) * baseline)
        cumulant = magnitude if order % 2 == 0 else -magnitude
        derivatives[order] = signed_hurwitz_gamma_derivative(
            order, t + flint.arb(1) / 2
        ) - cumulant / curvature ** (flint.arb(order) / 2)

    normalized_potential = {}
    for order in range(3, 9):
        epsilon = q ** (-flint.arb(1) / 2)
        value = jet[order] / curvature ** (flint.arb(order) / 2)
        normalized_potential[str(order)] = (
            value / epsilon ** (order - 2)
        ).str(30).replace("e", "E")
    return derivatives, {
        "mode": str(center),
        "mode_collar": [str(left), str(right)],
        "central_t": central_t.str(35).replace("e", "E"),
        "collar_t_left_upper": arb_upper_text(left_t),
        "collar_t_right_lower": arb_lower_text(right_t),
        "normalized_potential_L3_to_L8": normalized_potential,
    }


def conditional_collar_scout() -> list[dict]:
    rows = []
    for center, collar in SAMPLE_COLLARS:
        derivatives, diagnostics = candidate_h_cover(center, collar)
        result = evaluate_localized_curvature_from_h_cover(
            center,
            center,
            derivatives,
            cover_diagnostics=diagnostics,
        )
        if not result.get("passed"):
            raise RuntimeError(f"candidate cumulant corridor failed at u={center}")
        rows.append(
            {
                **diagnostics,
                "j_lower": result["j_lower"],
                "scaled_localized_upper": result["scaled_localized_upper"],
                "target_lower": result["target_lower"],
                "margin_lower": result["margin_lower"],
                "passed_conditionally": True,
            }
        )
    return rows


def corridor_formulas() -> dict:
    return {
        "kappa2": "2/5<=q*(kappa_2-1)<=4/5",
        "higher": {
            str(order): (
                f"1<=(-1)^{order}*kappa_{order}*q^({order}/2-1)/"
                f"{math.factorial(order-2)}<={cap}"
            )
            for order, cap in CORRIDOR_CAPS.items()
        },
    }


def formal_theorem_diagnostics() -> dict:
    finite = json.loads(FORMAL_FINITE_JSON.read_text(encoding="utf-8"))
    ray = json.loads(FORMAL_RAY_JSON.read_text(encoding="utf-8"))
    if finite.get("kind") != "jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate":
        raise RuntimeError("formal finite source has the wrong kind")
    if ray.get("kind") != "jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate":
        raise RuntimeError("formal ray source has the wrong kind")
    if ray.get("summary", {}).get("formal_ray_closed") is not True:
        raise RuntimeError("formal asymptotic ray is not closed")
    return {
        "finite": {
            "source": FORMAL_FINITE_JSON.relative_to(REPO_ROOT).as_posix(),
            "mode_interval": finite["parameters"]["mode_interval"],
            "block_count": finite["parameters"]["block_count"],
            "weakest_margin": finite["weakest_margin"],
        },
        "ray": {
            "source": FORMAL_RAY_JSON.relative_to(REPO_ROOT).as_posix(),
            "mode_interval": "u>=20",
            "formal_transfer": ray["scalar_geometry"]["formal_transfer"],
            "leading_corridor_gates": ray["summary"]["leading_corridor_gates"],
            "jet_remainder_sign_gates": ray["summary"]["jet_remainder_sign_gates"],
        },
        "global_formal_theorem": "candidate corridors hold for the epsilon^6 formal cumulants for every u>=2",
    }


def build_artifact() -> dict:
    exact = formal_cumulant_expansion()
    corridors = corridor_formulas()
    scouts = conditional_collar_scout()
    formal_theorems = formal_theorem_diagnostics()
    rows = [
        TargetRow(
            id="co4gcrt_01_epsilon_grading",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Grade each standardized potential derivative by its natural Gaussian ray power.",
            formula=exact["grading"],
            proof_boundary="Exact notation and scaling only.",
        ),
        TargetRow(
            id="co4gcrt_02_partition_recurrence",
            role="exact_formal_algebra",
            readiness="ready_to_apply",
            claim="Compute the tilted partition and its logarithm by exact truncated power-series recurrences.",
            formula=exact["partition_recurrence"],
            proof_boundary="Exact formal coefficients through epsilon^6; no analytic remainder bound.",
        ),
        TargetRow(
            id="co4gcrt_03_cumulant_polynomials",
            role="exact_formal_algebra",
            readiness="ready_to_apply",
            claim="Extract the complete formal cumulant polynomials needed for kappa_2 through kappa_8.",
            formula="kappa_r=partial_z^r[z^2/2+log Z_epsilon(z)] at z=0, through epsilon^6",
            proof_boundary="Exact finite formal expansion only.",
            diagnostics=exact["cumulants"],
        ),
        TargetRow(
            id="co4gcrt_04_factorial_signature",
            role="exact_formal_theorem",
            readiness="ready_to_apply",
            claim="At the exponential-potential limit L_r=1, the leading non-Gaussian cumulants have the alternating factorial signature.",
            formula="kappa_r~(-1)^r*(r-2)!*q^(1-r/2), r=3,...,8",
            proof_boundary="Leading formal signature only; not a bound for the exact density.",
            diagnostics=exact["leading_signature"],
        ),
        TargetRow(
            id="co4gcrt_05_candidate_corridor",
            role="analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Prove explicit two-sided exact cumulant corridors preserving the alternating factorial signature on u>=2.",
            formula="; ".join([corridors["kappa2"], *corridors["higher"].values()]),
            proof_boundary="Open analytic corridor theorem; the displayed constants are targets, not proved bounds.",
        ),
        TargetRow(
            id="co4gcrt_06_conditional_collar_scout",
            role="conditional_interval_scout",
            readiness="conditional",
            claim="Assuming the candidate corridors throughout each displayed t+-2 collar, Arb interval boxes clear the localized curvature ceiling at representative ray modes.",
            formula="candidate corridors => U(t)<7/(2*t^2) on the recorded point collars",
            proof_boundary="Conditional finite collar scout only; not a continuum ray certificate.",
            diagnostics={"rows": scouts},
        ),
        TargetRow(
            id="co4gcrt_06a_formal_finite_corridor",
            role="interval_certificate",
            readiness="ready_to_apply",
            claim="A deterministic Arb cover proves every candidate corridor for the epsilon-six formal cumulants on 2<=u<=20.",
            formula="formal cumulant corridors hold on 2<=u<=20",
            proof_boundary="Formal polynomial theorem only; exact-minus-formal density error remains open.",
            diagnostics=formal_theorems["finite"],
        ),
        TargetRow(
            id="co4gcrt_06b_formal_asymptotic_corridor",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="Coefficient-positive leading-model and jet-remainder gates prove every formal corridor on u>=20.",
            formula="formal cumulant corridors hold on u>=20",
            proof_boundary="Formal polynomial theorem only; exact-minus-formal density error remains open.",
            diagnostics=formal_theorems["ray"],
        ),
        TargetRow(
            id="co4gcrt_07_central_remainder_target",
            role="analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Bound the central Gaussian perturbation remainder after subtracting the exact epsilon^6 cumulant polynomial.",
            formula="|kappa_r-kappa_r^[6]|<=explicit C_r*epsilon^(next parity order), r=2,...,8",
            proof_boundary="Open central-window theorem.",
        ),
        TargetRow(
            id="co4gcrt_08_tail_and_composition_target",
            role="analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Control both adaptive tails and compose the remainder with exact L_r geometry to prove the candidate corridors for every u>=2.",
            formula="central remainder + two tails + exact mode geometry => cumulant corridors => curvature ray",
            proof_boundary="Sole analytic route proposed here; no ray theorem is asserted.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order4_gaussian_cumulant_ray_target",
        "date": "2026-07-12",
        "status": "exact formal cumulant algebra, global formal-corridor theorem, and conditional exact-density handoff",
        "proof_boundary": (
            "This artifact proves the formal epsilon^6 cumulant algebra, composes a "
            "global formal-corridor theorem, and checks that the same corridor is "
            "compatible with the localized curvature ceiling on finite t-collars. It "
            "does not prove the exact-density cumulant corridor, the continuum "
            "u>=2 ray, complete order-four entry, PF-infinity, RH, or Lambda<=0."
        ),
        "exact": exact,
        "candidate_corridors": corridors,
        "conditional_collar_scout": scouts,
        "formal_theorems": formal_theorems,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_rows": 4,
            "formal_theorem_rows": 2,
            "open_analytic_rows": 3,
            "conditional_scout_rows": 1,
            "positive_conditional_collars": len(scouts),
            "ready_to_apply_rows": sum(row.readiness == "ready_to_apply" for row in rows),
        },
        "sources": [
            "outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md",
            "outputs/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.md",
            "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.md",
            "outputs/formal_core.md",
        ],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.py"
        ),
        "remaining_target": (
            "Prove the exact-minus-formal cumulant errors fit inside the certified "
            "corridor margins by a central remainder and adaptive two-tail theorem."
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    corridors = artifact["candidate_corridors"]
    lines = [
        "# Jensen-Window PF Order-Four Gaussian Cumulant Ray Target",
        "",
        "Date: 2026-07-12",
        "",
        "Status: exact formal cumulant algebra, global formal-corridor theorem, and",
        "conditional exact-density handoff. This is not a proof of the exact-density",
        "cumulant corridor or continuum ray,",
        "order-four entry, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order4_gaussian_cumulant_ray_target`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.py",
        "```",
        "",
        "## Exact Expansion",
        "",
        "Put `epsilon=q^(-1/2)` and",
        "",
        "```text",
        exact["standardized_potential"],
        "```",
        "",
        "Exact tilted-Gaussian power-series algebra through `epsilon^6` gives:",
        "",
    ]
    for order in range(2, 9):
        row = exact["cumulants"][str(order)]
        terms = [str(row["gaussian_base"])] if row["gaussian_base"] else []
        terms.extend(
            f"epsilon^{term['epsilon_power']}*({term['coefficient']})"
            for term in row["terms"]
        )
        lines.extend(["```text", f"kappa_{order}=" + "+".join(terms), "```", ""])
    lines.extend(
        [
            "At `L_3=...=L_8=1`, the leading terms are",
            "",
            "```text",
            "kappa_3~-epsilon, kappa_4~2*epsilon^2, kappa_5~-6*epsilon^3,",
            "kappa_6~24*epsilon^4, kappa_7~-120*epsilon^5, kappa_8~720*epsilon^6.",
            "```",
            "",
            "## Candidate Corridor",
            "",
            "The exact remainder theorem should prove",
            "",
            "```text",
            corridors["kappa2"],
            *corridors["higher"].values(),
            "```",
            "",
            "These are theorem targets, not inferred from the formal series.",
            "",
            "## Proved Formal Corridor",
            "",
            f"A `{artifact['formal_theorems']['finite']['block_count']}`-block Arb certificate proves",
            "all seven corridors for the exact epsilon-six formal polynomial on",
            "`2<=u<=20`. Coefficient-positive analytic gates prove the same formal",
            "corridors on every `u>=20`. Thus the formal model is closed for all",
            "`u>=2`; only exact-minus-formal density errors remain.",
            "",
            "## Conditional Collar Test",
            "",
            "Assuming the corridor throughout each displayed mode collar, Arb boxes",
            "cover the full `t+-2` derivative interval and give:",
            "",
            "| u | central t | t^2 U upper | margin lower |",
            "|---:|---:|---:|---:|",
        ]
    )
    for row in artifact["conditional_collar_scout"]:
        lines.append(
            f"| `{row['mode']}` | `{row['central_t']}` | "
            f"`{row['scaled_localized_upper']}` | `{row['margin_lower']}` |"
        )
    lines.extend(
        [
            "",
            "This is a conditional finite compatibility check, not a continuum ray",
            "certificate. It shows that the proposed corridor constants retain enough",
            "margin for the localized order-four inequality.",
            "",
            "## Remaining Theorem",
            "",
            "Prove an explicit central-window remainder after subtracting the exact",
            "`epsilon^6` polynomial, control both adaptive tails, and compose those",
            "bounds with the exact normalized-potential geometry for every `u>=2`.",
            "That theorem would close the sole remaining curvature ray.",
            "",
            "```text",
            "outputs/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.md",
            "outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.md",
            "outputs/formal_core.md",
            "```",
        ]
    )
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
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print(
        "built order-four Gaussian cumulant ray target: "
        "10 rows, 4 exact formal rows, 2 proved formal-corridor rows, "
        "7 positive conditional collars, "
        "3 open analytic rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
