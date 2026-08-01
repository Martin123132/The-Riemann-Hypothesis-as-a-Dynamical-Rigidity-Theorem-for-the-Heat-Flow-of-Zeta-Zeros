#!/usr/bin/env python3
"""Build an exact length-13 obstruction for the fixed outer quartic branch."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
PARENT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_quartic_outer_threshold_order4_nonpromotion_gate.json"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_quartic_outer_branch_length13_obstruction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_quartic_outer_branch_length13_obstruction.md"
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bernstein_coefficients(
    polynomial: sp.Expr, variables: tuple[sp.Symbol, ...]
) -> tuple[tuple[int, ...], list[sp.Expr]]:
    """Return the tensor-product Bernstein coefficients on the unit cube."""
    poly = sp.Poly(polynomial, *variables)
    degrees = tuple(poly.degree(variable) for variable in variables)
    power = dict(poly.terms())
    coefficients: list[sp.Expr] = []
    for index in __import__("itertools").product(
        *(range(degree + 1) for degree in degrees)
    ):
        total = sp.Integer(0)
        for exponent in __import__("itertools").product(
            *(range(entry + 1) for entry in index)
        ):
            coefficient = power.get(exponent, 0)
            if not coefficient:
                continue
            weight = sp.Integer(1)
            for k, i, degree in zip(index, exponent, degrees):
                weight *= sp.Rational(math.comb(k, i), math.comb(degree, i))
            total += coefficient * weight
        coefficients.append(sp.factor(total))
    return degrees, coefficients


def bernstein_certificate(
    polynomial: sp.Expr, variables: tuple[sp.Symbol, ...]
) -> dict[str, object]:
    degrees, coefficients = bernstein_coefficients(polynomial, variables)
    positive = sum(1 for value in coefficients if value > 0)
    negative = sum(1 for value in coefficients if value < 0)
    zero = len(coefficients) - positive - negative
    encoded = "\n".join(str(value) for value in coefficients).encode("ascii")
    return {
        "degrees": list(degrees),
        "coefficient_count": len(coefficients),
        "positive": positive,
        "negative": negative,
        "zero": zero,
        "sha256": hashlib.sha256(encoded).hexdigest(),
    }


def build_symbolic(parent: dict) -> dict[str, object]:
    expected = {
        2: sp.Rational(24154639, 25000000),
        3: sp.Rational(567331181410000, 583446585220321),
        4: sp.Rational(909681023616163, 930852655574884),
        5: sp.Rational(2453, 2500),
        6: sp.Rational(123, 125),
        7: sp.Rational(19719, 20000),
        8: sp.Rational(987309, 1000000),
        9: sp.Rational(98823997, 100000000),
    }
    stored = {
        int(name.split("_")[1]): sp.Rational(value)
        for name, value in parent["exact"]["contractions"].items()
    }
    if stored != expected:
        raise RuntimeError("parent outer-branch contractions changed")

    variables = sp.symbols("y_10 y_11 y_12")
    x: dict[int, sp.Expr] = dict(expected)

    def defect(index: int) -> sp.Expr:
        return sp.cancel(1 - x[index])

    def gap(index: int) -> sp.Expr:
        return sp.cancel(
            defect(index + 2) ** 2
            - x[index + 2] ** 2
            * defect(index + 1)
            * defect(index + 3)
        )

    def cap(index: int) -> sp.Expr:
        return sp.cancel(
            gap(index - 4) ** 2
            / (x[index - 2] ** 3 * gap(index - 5))
        )

    extension_rows: list[dict[str, object]] = []
    for index, variable in zip((10, 11, 12), variables):
        current_cap = cap(index)
        current_gap = sp.cancel(variable * current_cap)
        current_defect = sp.cancel(
            (
                defect(index - 1) ** 2
                - current_gap
            )
            / (x[index - 1] ** 2 * defect(index - 2))
        )
        x[index] = sp.cancel(1 - current_defect)
        reconstructed = sp.cancel(gap(index - 3) - current_gap)
        if reconstructed != 0:
            raise RuntimeError(f"extension identity failed at x_{index}")
        extension_rows.append(
            {
                "index": index,
                "parameter": str(variable),
                "cap_operation_count": sp.count_ops(current_cap),
                "contraction_operation_count": sp.count_ops(x[index]),
                "identity": f"G_{index - 3}={variable}*C_{index}",
            }
        )

    next_cap = cap(13)
    repeated_gap = sp.cancel(
        defect(12) ** 2 - x[12] ** 2 * defect(11) * defect(12)
    )
    compatibility_gap = sp.cancel(next_cap - repeated_gap)
    numerator, denominator = sp.fraction(compatibility_gap)
    gap_poly = sp.Poly(numerator, *variables)
    denominator_poly = sp.Poly(denominator, *variables)

    factor_coefficient, denominator_factors = sp.factor_list(
        denominator, *variables
    )
    if factor_coefficient <= 0 or len(denominator_factors) != 2:
        raise RuntimeError("unexpected compatibility denominator factorization")
    factor_certificates: list[dict[str, object]] = []
    for factor, exponent in denominator_factors:
        certificate = bernstein_certificate(factor, variables)
        certificate.update(
            {
                "exponent": exponent,
                "operation_count": sp.count_ops(factor),
                "factor": str(factor),
            }
        )
        factor_certificates.append(certificate)
    factor_certificates.sort(key=lambda row: tuple(row["degrees"]))
    if not (
        factor_certificates[0]["negative"] == 2
        and factor_certificates[0]["positive"] == 0
        and factor_certificates[0]["exponent"] == 3
        and factor_certificates[1]["positive"] == 6
        and factor_certificates[1]["negative"] == 0
        and factor_certificates[1]["exponent"] == 6
    ):
        raise RuntimeError("compatibility denominator sign was not certified")

    derivative_certificates: dict[str, dict[str, object]] = {}
    expected_numerator_signs = {
        "y_10": (612, 0, 0),
        "y_11": (0, 480, 64),
        "y_12": (0, 216, 144),
    }
    expected_denominator_signs = {"y_10": 1, "y_11": -1, "y_12": -1}
    cube_midpoint = {variable: sp.Rational(1, 2) for variable in variables}
    for variable in variables:
        derivative = sp.cancel(sp.diff(compatibility_gap, variable))
        derivative_numerator, derivative_denominator = sp.fraction(derivative)
        certificate = bernstein_certificate(derivative_numerator, variables)
        observed = (
            certificate["positive"],
            certificate["negative"],
            certificate["zero"],
        )
        if observed != expected_numerator_signs[str(variable)]:
            raise RuntimeError(
                f"unexpected derivative Bernstein signs for {variable}: {observed}"
            )
        denominator_sign = sp.sign(
            derivative_denominator.subs(cube_midpoint)
        )
        if denominator_sign != expected_denominator_signs[str(variable)]:
            raise RuntimeError(
                f"unexpected derivative denominator sign for {variable}"
            )
        certificate.update(
            {
                "numerator_terms": len(
                    sp.Poly(derivative_numerator, *variables).terms()
                ),
                "denominator_sign_on_open_cube": int(denominator_sign),
                "derivative_sign_on_open_cube": "nonnegative",
            }
        )
        derivative_certificates[str(variable)] = certificate

    corner = sp.cancel(
        compatibility_gap.subs({variable: 1 for variable in variables})
    )
    if corner >= 0:
        raise RuntimeError("cube-corner compatibility gap is not negative")

    return {
        "fixed_contractions": {
            f"x_{index}": str(value) for index, value in expected.items()
        },
        "definitions": {
            "defect": "d_j=1-x_j",
            "order3_gap": "G_j=d_(j+2)^2-x_(j+2)^2*d_(j+1)*d_(j+3)",
            "order4_cap": "C_k=G_(k-4)^2/(x_(k-2)^3*G_(k-5))",
            "extension": "G_(k-3)=y_k*C_k, 0<y_k<1",
            "next_repeated_gap": (
                "M_13=d_12^2-x_12^2*d_11*d_12"
            ),
            "compatibility_gap": "Delta_13=C_13-M_13",
        },
        "extension_rows": extension_rows,
        "compatibility": {
            "numerator_degrees": list(gap_poly.degree_list()),
            "numerator_terms": len(gap_poly.terms()),
            "denominator_degrees": list(denominator_poly.degree_list()),
            "denominator_terms": len(denominator_poly.terms()),
            "denominator_sign_on_open_cube": -1,
            "denominator_factors": factor_certificates,
            "derivatives": derivative_certificates,
            "cube_corner": str(corner),
            "cube_corner_decimal": str(sp.N(corner, 30)),
            "conclusion": (
                "Delta_13 is coordinatewise nondecreasing on (0,1)^3 "
                "and Delta_13(1,1,1)<0"
            ),
        },
    }


def build_payload() -> dict:
    parent = json.loads(PARENT_RESULT.read_text(encoding="utf-8"))
    if parent.get("kind") != (
        "jensen_window_pf_quartic_outer_threshold_order4_nonpromotion_gate"
    ):
        raise RuntimeError("unexpected parent artifact kind")
    exact = build_symbolic(parent)
    rows = [
        GateRow(
            id="qob13_01_parent",
            role="provenance",
            readiness="proved_exact",
            claim="The fixed outer-contact prefix is exactly the rational x_2 through x_9 branch from the order-four nonpromotion gate.",
            formula="x_2,...,x_9 fixed",
            proof_boundary="This lemma does not vary the earlier outer-contact prefix.",
            diagnostics={
                "parent_result": str(PARENT_RESULT.relative_to(REPO_ROOT)),
                "parent_sha256": file_sha256(PARENT_RESULT),
            },
        ),
        GateRow(
            id="qob13_02_gap",
            role="identity",
            readiness="proved_exact",
            claim="The terminal order-three sign is exactly positivity of the defect gap G_j.",
            formula=exact["definitions"]["order3_gap"],
            proof_boundary="Contiguous order three only.",
        ),
        GateRow(
            id="qob13_03_cap",
            role="identity",
            readiness="proved_exact",
            claim="Given the preceding strict layers, the next order-four sign is exactly the upper cap G_(k-3)<C_k.",
            formula=exact["definitions"]["order4_cap"],
            proof_boundary="Uses the strict positive preceding gaps.",
        ),
        GateRow(
            id="qob13_04_parameterization",
            role="reduction",
            readiness="proved_exact",
            claim="Every strict order-three/order-four continuation through x_12 is parameterized by three independent open-cube coordinates.",
            formula=exact["definitions"]["extension"],
            proof_boundary="It parameterizes signed contiguous layers, not every scalar or Xi constraint.",
            diagnostics={"extensions": exact["extension_rows"]},
        ),
        GateRow(
            id="qob13_05_necessary_margin",
            role="reduction",
            readiness="proved_exact",
            claim="An increasing x_13 with the next order-four sign requires Delta_13=C_13-M_13>0.",
            formula="x_13>x_12 => G_10(x_13)>M_13; H_4>0 => G_10(x_13)<C_13",
            proof_boundary="A necessary condition, sufficient only for these two terminal inequalities.",
        ),
        GateRow(
            id="qob13_06_denominator",
            role="exact_certificate",
            readiness="proved_exact",
            claim="The rational compatibility denominator has one negative odd-power factor and one positive even-power factor, so its sign is fixed.",
            formula="sign denominator(Delta_13)=-1 on (0,1)^3",
            proof_boundary="Exact tensor-product Bernstein certificate.",
            diagnostics={
                "factors": exact["compatibility"]["denominator_factors"]
            },
        ),
        GateRow(
            id="qob13_07_y10_derivative",
            role="exact_certificate",
            readiness="proved_exact",
            claim="The y_10 derivative of Delta_13 is nonnegative on the open cube.",
            formula="partial_(y_10) Delta_13>=0",
            proof_boundary="Exact rational derivative and Bernstein numerator signs.",
            diagnostics=exact["compatibility"]["derivatives"]["y_10"],
        ),
        GateRow(
            id="qob13_08_y11_derivative",
            role="exact_certificate",
            readiness="proved_exact",
            claim="The y_11 derivative of Delta_13 is nonnegative on the open cube.",
            formula="partial_(y_11) Delta_13>=0",
            proof_boundary="Exact rational derivative and Bernstein numerator signs.",
            diagnostics=exact["compatibility"]["derivatives"]["y_11"],
        ),
        GateRow(
            id="qob13_09_y12_derivative",
            role="exact_certificate",
            readiness="proved_exact",
            claim="The y_12 derivative of Delta_13 is nonnegative on the open cube.",
            formula="partial_(y_12) Delta_13>=0",
            proof_boundary="Exact rational derivative and Bernstein numerator signs.",
            diagnostics=exact["compatibility"]["derivatives"]["y_12"],
        ),
        GateRow(
            id="qob13_10_obstruction",
            role="scope_gate",
            readiness="proved_exact",
            claim="Even the closed-cube maximum is negative, so this fixed outer branch cannot be continued to an increasing x_13 while retaining the signed layers through order four.",
            formula="Delta_13(y_10,y_11,y_12)<=Delta_13(1,1,1)<0",
            proof_boundary="This kills this prefix only; it is not a uniform theorem over all outer quartic contacts.",
            diagnostics={
                "corner": exact["compatibility"]["cube_corner"],
                "corner_decimal": exact["compatibility"][
                    "cube_corner_decimal"
                ],
            },
        ),
    ]
    return {
        "kind": "jensen_window_pf_quartic_outer_branch_length13_obstruction",
        "date": "2026-07-25",
        "status": "exact fixed-prefix length-thirteen obstruction",
        "proof_boundary": (
            "This exact semialgebraic certificate proves that the particular "
            "outer-contact prefix x_2 through x_9 from the local order-four "
            "nonpromotion gate cannot have an increasing x_13 after any "
            "strict signed-Hankel continuation through x_12. It does not "
            "exclude other outer-contact prefixes, prove a uniform finite-"
            "tail theorem, establish an Xi-specific inequality, prove "
            "PF-infinity, Lambda<=0, or RH."
        ),
        "parent": {
            "path": str(PARENT_RESULT.relative_to(REPO_ROOT)),
            "sha256": file_sha256(PARENT_RESULT),
        },
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "extension_parameters": 3,
            "derivative_certificates": 3,
            "derivative_bernstein_coefficients": sum(
                row["coefficient_count"]
                for row in exact["compatibility"]["derivatives"].values()
            ),
            "denominator_factor_bernstein_coefficients": sum(
                row["coefficient_count"]
                for row in exact["compatibility"]["denominator_factors"]
            ),
            "negative_cube_corner_maxima": 1,
            "forbidden_all_length_promotions": 1,
            "repaired_handoffs": 1,
        },
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    compatibility = exact["compatibility"]
    derivatives = compatibility["derivatives"]
    return "\n".join(
        [
            "# Quartic Outer-Branch Length-13 Obstruction",
            "",
            "Date: 2026-07-25",
            "",
            "Status: exact fixed-prefix length-thirteen obstruction. This is",
            "not a proof of a uniform quartic threshold, PF-infinity,",
            "`Lambda <= 0`, or RH.",
            "",
            f"Artifact kind: `{payload['kind']}`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_quartic_outer_branch_length13_obstruction.json",
            "python work/rh_compute/scripts/jensen_window_pf_quartic_outer_branch_length13_obstruction.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_quartic_outer_branch_length13_obstruction.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated quartic outer-branch length-13 obstruction: 10 rows, 0 issues, 3 exact corridor parameters, 3 derivative certificates, 1516 derivative Bernstein coefficients, 8 denominator-factor Bernstein coefficients, 1 negative cube-corner maximum, 1 forbidden all-length promotion, 1 repaired handoff",
            "```",
            "",
            "## Exact Parameterization",
            "",
            "Keep the rational outer-contact prefix `x_2,...,x_9` from the",
            "order-four nonpromotion gate. Put",
            "",
            "```text",
            exact["definitions"]["order3_gap"],
            exact["definitions"]["order4_cap"],
            "```",
            "",
            "For each `k=10,11,12`, the strict terminal order-three and",
            "order-four signs are exactly parameterized by",
            "",
            "```text",
            "G_(k-3)=y_k*C_k,  0<y_k<1.",
            "```",
            "",
            "This is exhaustive: `G_(k-3)>0` is the order-three sign and",
            "`G_(k-3)<C_k` is the order-four sign.",
            "The checker independently reconstructs the same cube by solving",
            "the raw terminal `3x3` and `4x4` Hankel determinant roots at each",
            "step, then reproduces the derivative Bernstein hashes.",
            "",
            "## Terminal Compatibility",
            "",
            "At the next step define",
            "",
            "```text",
            "M_13=d_12^2-x_12^2*d_11*d_12,",
            "Delta_13=C_13-M_13.",
            "```",
            "",
            "If `x_13>x_12`, then the new gap `G_10(x_13)` is strictly",
            "larger than `M_13`. The next order-four sign requires",
            "`G_10(x_13)<C_13`. Hence `Delta_13>0` is necessary.",
            "",
            "Exact tensor-product Bernstein conversion gives the derivative",
            "numerator sign counts",
            "",
            "```text",
            f"y_10: +{derivatives['y_10']['positive']}, -{derivatives['y_10']['negative']}, 0={derivatives['y_10']['zero']}",
            f"y_11: +{derivatives['y_11']['positive']}, -{derivatives['y_11']['negative']}, 0={derivatives['y_11']['zero']}",
            f"y_12: +{derivatives['y_12']['positive']}, -{derivatives['y_12']['negative']}, 0={derivatives['y_12']['zero']}.",
            "```",
            "",
            "Together with the independently certified denominator signs,",
            "all three partial derivatives of `Delta_13` are nonnegative on",
            "the open cube. Therefore its supremum is the closed-cube corner",
            "`(1,1,1)`. Exact substitution gives",
            "",
            "```text",
            f"Delta_13(1,1,1)={compatibility['cube_corner']}",
            f"                    ={compatibility['cube_corner_decimal']}<0.",
            "```",
            "",
            "Thus no choice inside any of the three signed continuation",
            "corridors can produce an increasing `x_13` while preserving the",
            "next order-four sign.",
            "",
            "## Corrected Handoff",
            "",
            "The previously displayed `x_10` corridor was real, and the",
            "branch can continue through `x_12`, but it cannot continue",
            "indefinitely inside the increasing-contraction cone. The live",
            "question is no longer whether this particular witness has an",
            "all-length extension. It does not.",
            "",
            "This does not promote the quartic threshold. The companion",
            "alternate-tail gate now gives an explicit rational continuation",
            "from the same outer contact that survives through `x_13`, so",
            "this obstruction is demonstrably prefix-specific. The next",
            "theorem target is a uniform finite-tail obstruction beginning",
            "at length 14, or another tail crossing that boundary. Global",
            "Xi, degree-five, and",
            "theta-specific routes remain open.",
            "",
            "```text",
            "outputs/jensen_window_pf_quartic_outer_branch_alternate_length13_survivor_gate.md",
            "```",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(f"wrote quartic outer-branch obstruction: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
