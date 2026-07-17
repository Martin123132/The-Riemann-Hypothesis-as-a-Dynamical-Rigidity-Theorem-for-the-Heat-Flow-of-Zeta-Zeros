#!/usr/bin/env python3
"""Audit Toda propagation, the shifted-tail boundary, and the Jensen gap."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_deep_schur_toda_boundary_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_deep_schur_toda_boundary_gate.md"
)
DEEP_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_endpoint_deep_schur_coordinate.json"
)
HEAT_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_all_order_endpoint_heat_reduction.json"
)
SOURCE_PATHS = (DEEP_SOURCE, HEAT_SOURCE)

SYMBOLIC_TODA_ORDERS = tuple(range(1, 5))
BOUNDARY_MAX_ORDER = 5
BOUNDARY_MAX_PART = 4
STRICT_SCHUR_MAX_SIZE = 10
COUNTEREXAMPLE_EPSILON = Fraction(1, 100)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_record(path: Path, payload: dict) -> dict:
    return {
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "sha256": sha256(path),
        "kind": payload.get("kind"),
        "status": payload.get("status"),
    }


def validate_sources() -> tuple[list[dict], dict]:
    deep = load_json(DEEP_SOURCE)
    heat = load_json(HEAT_SOURCE)
    if deep.get("summary", {}).get("endpoint_deep_equivalences") != 1:
        raise RuntimeError("deep-Schur source does not contain the endpoint equivalence")
    if deep.get("summary", {}).get("rigorous_endpoint_pf_counterexamples") != 1:
        raise RuntimeError("deep-Schur source does not contain the endpoint PF boundary")
    if deep.get("summary", {}).get("negative_deep_rectangles") != 4:
        raise RuntimeError("order-ten deep-rectangle counterexample source changed")
    if deep.get("summary", {}).get("rejected_rectangle_hierarchies") != 1:
        raise RuntimeError("deep rectangle hierarchy is not marked rejected")
    if heat.get("summary", {}).get("endpoint_interval_equivalences") != 1:
        raise RuntimeError("heat source does not contain the endpoint/interval equivalence")
    expected_rectangle = (
        "Q_(m,n)(-100)=A_0(-100)^m*s_((n+m-1)^m)(h)"
    )
    if deep.get("exact", {}).get("rectangle_identity") != expected_rectangle:
        raise RuntimeError("deep-Schur rectangle identity changed")
    return (
        [source_record(DEEP_SOURCE, deep), source_record(HEAT_SOURCE, heat)],
        {
            "rectangle_identity": expected_rectangle,
            "rectangle_target": deep["exact"]["rectangle_target"],
            "deep_failure": deep["exact"]["deep_failure"],
            "rejected_rectangle_target": deep["exact"][
                "rejected_rectangle_target"
            ],
            "negative_deep_rectangles": deep["summary"][
                "negative_deep_rectangles"
            ],
            "rejected_rectangle_hierarchies": deep["summary"][
                "rejected_rectangle_hierarchies"
            ],
            "deep_support_boundary": deep["exact"]["support_boundary"],
            "endpoint_heat_equivalence": heat["exact"][
                "endpoint_interval_equivalence"
            ],
        },
    )


def determinant(matrix: list[list[Fraction]]) -> Fraction:
    size = len(matrix)
    if size == 0:
        return Fraction(1)
    work = [row[:] for row in matrix]
    value = Fraction(1)
    for column in range(size):
        pivot_row = next(
            (row for row in range(column, size) if work[row][column] != 0),
            None,
        )
        if pivot_row is None:
            return Fraction(0)
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            value = -value
        pivot = work[column][column]
        value *= pivot
        for row in range(column + 1, size):
            quotient = work[row][column] / pivot
            for index in range(column + 1, size):
                work[row][index] -= quotient * work[column][index]
    return value


def fraction_record(value: Fraction) -> dict:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "text": str(value),
    }


def hankel_symbolic(values: tuple[sp.Symbol, ...], order: int, shift: int) -> sp.Expr:
    if order == 0:
        return sp.Integer(1)
    return sp.Matrix(
        [[values[shift + i + j] for j in range(order)] for i in range(order)]
    ).det(method="domain-ge")


def signed_hankel_symbolic(
    values: tuple[sp.Symbol, ...], order: int, shift: int
) -> sp.Expr:
    return (-1) ** (order * (order - 1) // 2) * hankel_symbolic(
        values, order, shift
    )


def symbolic_toda_audit() -> dict:
    rows = []
    for order in SYMBOLIC_TODA_ORDERS:
        values = sp.symbols(f"a0:{2 * order + 1}")
        residual = sp.factor(
            signed_hankel_symbolic(values, order + 1, 0)
            * signed_hankel_symbolic(values, order - 1, 2)
            - signed_hankel_symbolic(values, order, 1) ** 2
            + signed_hankel_symbolic(values, order, 0)
            * signed_hankel_symbolic(values, order, 2)
        )
        if residual != 0:
            raise RuntimeError(f"Toda identity failed at central order {order}")
        rows.append({"central_order": order, "residual": "0"})
    return {
        "orders": list(SYMBOLIC_TODA_ORDERS),
        "rows": rows,
        "checks": len(rows),
        "all_residuals_zero": True,
        "scope": (
            "Finite symbolic specializations audit the index and sign convention; "
            "arbitrary-order validity is Desnanot-Jacobi."
        ),
    }


def partitions_in_box(length: int, max_part: int):
    def visit(prefix: tuple[int, ...], ceiling: int):
        if len(prefix) == length:
            yield prefix
            return
        for part in range(ceiling, -1, -1):
            yield from visit(prefix + (part,), part)

    yield from visit((), max_part)


def partitions_of_size(total: int, ceiling: int | None = None):
    if total == 0:
        yield ()
        return
    if ceiling is None or ceiling > total:
        ceiling = total
    for first in range(ceiling, 0, -1):
        for rest in partitions_of_size(total - first, first):
            yield (first,) + rest


def counterexample_h(index: int) -> Fraction:
    if index < 0:
        return Fraction(0)
    epsilon = COUNTEREXAMPLE_EPSILON
    return sum(
        epsilon**power
        / math.factorial(power)
        * (2 ** (index - power + 1) - 1)
        for power in range(index + 1)
    )


def boundary_reset_audit() -> dict:
    checks = 0
    deep_checks = 0
    shallow_checks = 0
    shallow_mismatches = 0
    for order in range(1, BOUNDARY_MAX_ORDER + 1):
        shift = order - 1
        normalizer = counterexample_h(shift)
        for partition in partitions_in_box(order, BOUNDARY_MAX_PART):
            checks += 1
            tail_matrix = []
            translated_matrix = []
            for row in range(order):
                tail_row = []
                translated_row = []
                for column in range(order):
                    index = partition[row] - row + column
                    tail_row.append(
                        counterexample_h(shift + index) / normalizer
                        if index >= 0
                        else Fraction(0)
                    )
                    translated_row.append(
                        counterexample_h(shift + index) / normalizer
                    )
                tail_matrix.append(tail_row)
                translated_matrix.append(translated_row)
            tail_value = determinant(tail_matrix)
            translated_value = determinant(translated_matrix)
            if partition[-1] >= order - 1:
                deep_checks += 1
                if tail_value != translated_value:
                    raise RuntimeError(
                        f"deep shifted-tail identity failed at {partition}"
                    )
            else:
                shallow_checks += 1
                if tail_value != translated_value:
                    shallow_mismatches += 1

    h0 = counterexample_h(0)
    h1 = counterexample_h(1)
    h2 = counterexample_h(2)
    tail_empty = Fraction(1)
    translated_empty = (h1 * h1 - h0 * h2) / (h1 * h1)
    defect = tail_empty - translated_empty
    if defect != h0 * h2 / (h1 * h1) or defect <= 0:
        raise RuntimeError("order-two boundary-reset defect failed")
    if shallow_mismatches != shallow_checks:
        raise RuntimeError("bounded shallow boundary audit changed")
    return {
        "checks": checks,
        "deep_translation_checks": deep_checks,
        "shallow_checks": shallow_checks,
        "shallow_mismatches": shallow_mismatches,
        "order_two_witness": {
            "order": 2,
            "shift": 1,
            "partition": [0, 0],
            "tail_schur": fraction_record(tail_empty),
            "translated_deep_schur": fraction_record(translated_empty),
            "defect": fraction_record(defect),
            "symbolic_defect": "h_0*h_2/h_1^2",
        },
        "exact_fit_condition": (
            "det[b_(mu_i-i+j)] = h_s^(-r)*s_(mu+(s^r))(h) only when "
            "mu_r>=r-1, so every tail Jacobi-Trudi index is nonnegative"
        ),
    }


def hook_product(partition: tuple[int, ...]) -> int:
    product = 1
    for row, row_length in enumerate(partition):
        for column in range(row_length):
            below = sum(
                1 for later_length in partition[row + 1 :] if later_length > column
            )
            product *= row_length - column + below
    return product


def schur_value(partition: tuple[int, ...]) -> Fraction:
    order = len(partition)
    return determinant(
        [
            [counterexample_h(partition[row] - row + column) for column in range(order)]
            for row in range(order)
        ]
    )


def strict_schur_audit() -> dict:
    checks = 0
    weakest_ratio: Fraction | None = None
    weakest_partition: tuple[int, ...] | None = None
    for size in range(1, STRICT_SCHUR_MAX_SIZE + 1):
        for partition in partitions_of_size(size):
            value = schur_value(partition)
            tableaux = math.factorial(size) // hook_product(partition)
            plancherel_term = (
                Fraction(tableaux, math.factorial(size))
                * COUNTEREXAMPLE_EPSILON**size
            )
            if value < plancherel_term or plancherel_term <= 0:
                raise RuntimeError(f"strict Schur audit failed at {partition}")
            ratio = value / plancherel_term
            if weakest_ratio is None or ratio < weakest_ratio:
                weakest_ratio = ratio
                weakest_partition = partition
            checks += 1
    return {
        "checks": checks,
        "max_partition_size": STRICT_SCHUR_MAX_SIZE,
        "all_values_at_least_plancherel_term": True,
        "weakest_ratio": fraction_record(weakest_ratio or Fraction(0)),
        "weakest_partition": list(weakest_partition or ()),
        "all_order_strictness_formula": (
            "s_lambda[X+Pl_epsilon]>=s_lambda[Pl_epsilon]="
            "f^lambda*epsilon^|lambda|/|lambda|!>0"
        ),
    }


def cubic_discriminant(a: Fraction, b: Fraction, c: Fraction, d: Fraction) -> Fraction:
    return (
        18 * a * b * c * d
        - 4 * b**3 * d
        + b**2 * c**2
        - 4 * a * c**3
        - 27 * a**2 * d**2
    )


def strict_schur_jensen_counterexample() -> dict:
    coefficients = [
        Fraction(math.comb(3, index)) * counterexample_h(index)
        for index in range(4)
    ]
    discriminant = cubic_discriminant(
        coefficients[3], coefficients[2], coefficients[1], coefficients[0]
    )
    if discriminant >= 0:
        raise RuntimeError("strict-Schur Jensen counterexample lost its sign")
    common_denominator = math.lcm(*(value.denominator for value in coefficients))
    primitive = [int(value * common_denominator) for value in coefficients]
    divisor = math.gcd(*primitive)
    primitive = [value // divisor for value in primitive]
    return {
        "epsilon": fraction_record(COUNTEREXAMPLE_EPSILON),
        "generating_function": "exp(z/100)/((1-z)*(1-2*z))",
        "edrei_parameters": {
            "exponential": "1/100",
            "denominator_poles": ["1", "2"],
            "numerator_zeros": [],
        },
        "strict_schur_reason": (
            "the Plancherel summand alone is strictly positive for every partition"
        ),
        "h_0_to_h_3": [fraction_record(counterexample_h(index)) for index in range(4)],
        "jensen_degree": 3,
        "jensen_shift": 0,
        "jensen_coefficients_ascending": [
            fraction_record(value) for value in coefficients
        ],
        "primitive_polynomial_coefficients_ascending": primitive,
        "primitive_polynomial": (
            "6000000+54180000*x+126540900*x^2+90420901*x^3"
        ),
        "discriminant": fraction_record(discriminant),
        "strictly_negative_discriminant": True,
        "conclusion": (
            "one real zero and one nonreal conjugate pair; the cubic Jensen "
            "polynomial is not hyperbolic"
        ),
    }


def literature_audit() -> list[dict]:
    return [
        {
            "id": "lam_postnikov_pylyavskyy_schur_log_concavity",
            "classification": "formal_identity_support_only",
            "finding": (
                "Schur log-concavity supports formal rectangular difference identities, "
                "but numerical positivity after specialization still needs the relevant "
                "Schur values to be positive."
            ),
            "does_not_supply": "the zeta endpoint rectangle signs",
            "url": "https://arxiv.org/abs/math/0502446",
        },
        {
            "id": "wagner_hadamard_products",
            "classification": "hypotheses_not_met",
            "finding": (
                "Hadamard products of totally positive Toeplitz matrices are not closed "
                "in general; Wagner proves narrower sufficient classes."
            ),
            "does_not_supply": (
                "a conversion from deep endpoint minors to binomially weighted Jensen minors"
            ),
            "url": "https://doi.org/10.1016/0022-247X(92)90261-B",
        },
        {
            "id": "angarone_kim_oh_soskin_dual_jacobi_trudi",
            "classification": "restricted_shape_only",
            "finding": (
                "The newer preprint leaves the general Hadamard Jacobi-Trudi conjecture "
                "open and proves a 3x2-avoiding ribbon-like case."
            ),
            "does_not_supply": (
                "the open rectangles (N^m), which contain 3x2 blocks for m>=3,N>=2"
            ),
            "url": "https://arxiv.org/abs/2511.08969",
        },
        {
            "id": "craven_csordas_fox_wright",
            "classification": "fixed_multiplier_family_support",
            "finding": (
                "Reciprocal-Gamma arithmetic progressions provide genuine multiplier "
                "sequences and explain the fixed-scale model."
            ),
            "does_not_supply": (
                "closure under the Xi scale mixture or an endpoint-to-Jensen theorem"
            ),
            "url": "https://doi.org/10.1016/j.jmaa.2005.03.058",
        },
    ]


def build_artifact() -> dict:
    sources, source_contract = validate_sources()
    toda = symbolic_toda_audit()
    boundary = boundary_reset_audit()
    strict_audit = strict_schur_audit()
    counterexample = strict_schur_jensen_counterexample()
    literature = literature_audit()
    exact = {
        "tau_coordinate": (
            "tau_(m,N)=s_((N^m))(h)=Q_(m,N-m+1)(-100)/A_0(-100)^m"
        ),
        "toda_identity": (
            "tau_(m+1,N)*tau_(m-1,N)=tau_(m,N)^2-"
            "tau_(m,N-1)*tau_(m,N+1)"
        ),
        "next_order_criterion": (
            "tau_(m+1,N)>0 iff tau_(m,N)^2>"
            "tau_(m,N-1)*tau_(m,N+1), assuming tau_(m-1,N)>0"
        ),
        "ordinary_tail": (
            "b_k^(s)=h_(s+k)/h_s for k>=0 and b_k^(s)=0 for k<0"
        ),
        "conditional_tail_translation": boundary["exact_fit_condition"],
        "boundary_defect": (
            "s_(0,0)(b^(1))-h_1^(-2)*s_(1,1)(h)=h_0*h_2/h_1^2>0"
        ),
        "generic_bridge_counterexample": (
            "H(z)=exp(z/100)/((1-z)(1-2z)) is strictly Schur-positive, "
            "but J_(3,0)(x) has negative discriminant"
        ),
        "rectangle_target": source_contract["rectangle_target"],
        "deep_failure": source_contract["deep_failure"],
        "rejected_rectangle_target": source_contract[
            "rejected_rectangle_target"
        ],
        "bridge_target": (
            "prove every binomially weighted Jensen window for the actual Xi/Phi "
            "sequence using hypotheses stronger than unweighted Schur positivity "
            "and weaker than all-shift signed-Hankel positivity"
        ),
    }
    rows = [
        GateRow(
            "jwpfdstb_01_input_coordinate",
            "source_contract",
            "ready_to_apply",
            "Use the normalized endpoint specialization and exact rectangular deep-Schur map.",
            source_contract["rectangle_identity"],
            "Imports the exact coordinate and its separately certified order-ten failure.",
        ),
        GateRow(
            "jwpfdstb_02_tau_coordinate",
            "exact_reindexing",
            "ready_to_apply",
            "Rectangular endpoint Schur values form a two-index tau array.",
            exact["tau_coordinate"],
            "A notation change only.",
        ),
        GateRow(
            "jwpfdstb_03_toda_identity",
            "exact_determinant_identity",
            "ready_to_apply",
            "Desnanot-Jacobi gives the rectangular discrete Toda relation.",
            exact["toda_identity"],
            "The right side is a difference, so the identity alone has no positivity direction.",
            toda,
        ),
        GateRow(
            "jwpfdstb_04_log_concavity_criterion",
            "exact_sign_equivalence",
            "ready_to_apply",
            "The next rectangle order is exactly strict width-log-concavity of the current row; the first required row fails at four widths.",
            exact["next_order_criterion"],
            "This restates rather than proves the next-order inequality.",
        ),
        GateRow(
            "jwpfdstb_05_no_cluster_promotion",
            "route_rejection_guard",
            "ready_to_apply",
            "Ordinary subtraction-free cluster positivity cannot be read off from this signed Toda form.",
            "positive square minus positive neighbor product has undetermined sign",
            "A different zeta-specific invariant could still control the difference.",
        ),
        GateRow(
            "jwpfdstb_06_ordinary_tail",
            "candidate_coordinate",
            "ready_to_apply",
            "Resetting a shifted coefficient tail creates a genuine one-sided Schur specialization.",
            exact["ordinary_tail"],
            "The reset changes negative-index entries and must be audited before any PF claim.",
        ),
        GateRow(
            "jwpfdstb_07_conditional_translation",
            "exact_boundary_theorem",
            "ready_to_apply",
            "Tail translation agrees with the shifted deep Schur determinant only when the tail shape is itself deep.",
            exact["conditional_tail_translation"],
            "It does not extend positivity to shallow tail shapes.",
            boundary,
        ),
        GateRow(
            "jwpfdstb_08_boundary_counterexample",
            "exact_counterexample",
            "ready_to_apply",
            "The smallest shallow tail shape already has a nonzero reset defect.",
            exact["boundary_defect"],
            "This rejects the proposed moving-tail PF equivalence independently of the direct order-ten deep failures.",
            boundary["order_two_witness"],
        ),
        GateRow(
            "jwpfdstb_09_no_lr_tail_lift",
            "route_rejection_guard",
            "ready_to_apply",
            "Littlewood-Richardson positivity cannot fill the missing tail boundary because its required straight tail Schurs are not translated deep Schurs.",
            "moving-tail PF_m does not follow from the deep cone by zero-reset translation",
            "A strip/interior Toeplitz theorem with the original boundary remains conceivable.",
        ),
        GateRow(
            "jwpfdstb_10_strict_schur_model",
            "exact_countermodel_construction",
            "ready_to_apply",
            "An Edrei specialization with a positive exponential component is strictly Schur-positive for every partition.",
            counterexample["generating_function"],
            "This is a generic proof-safety model, not the zeta endpoint sequence.",
        ),
        GateRow(
            "jwpfdstb_11_all_shape_strictness",
            "exact_all_order_theorem",
            "ready_to_apply",
            "The Plancherel summand proves strict positivity at every Schur shape in the countermodel.",
            strict_audit["all_order_strictness_formula"],
            "The bounded determinant audit checks implementation; the displayed formula is all-order.",
            strict_audit,
        ),
        GateRow(
            "jwpfdstb_12_jensen_failure",
            "exact_counterexample",
            "ready_to_apply",
            "The same strictly Schur-positive sequence has a nonhyperbolic cubic Jensen polynomial.",
            counterexample["primitive_polynomial"]
            + "; discriminant="
            + counterexample["discriminant"]["text"],
            "This proves unweighted Schur/PF positivity alone cannot be the Jensen bridge.",
            counterexample,
        ),
        GateRow(
            "jwpfdstb_13_literature_fit",
            "primary_literature_gate",
            "ready_to_apply",
            "The audited Toda, Hadamard, Jacobi-Trudi, and multiplier results support formal pieces but do not supply the surviving Xi/Phi Jensen theorem.",
            "direct_closing_theorems_in_audited_set=0",
            "A bounded route audit, not an exhaustive literature claim.",
            {"rows": literature},
        ),
        GateRow(
            "jwpfdstb_14_rejected_rectangles",
            "countermodel_gate",
            "rejected_by_counterexample",
            "The zeta-specific endpoint rectangular hierarchy fails at its first uncompleted order.",
            exact["deep_failure"] + "; " + exact["rejected_rectangle_target"],
            "Rejects this all-order signed-Hankel route, not RH or Jensen hyperbolicity.",
        ),
        GateRow(
            "jwpfdstb_15_xi_specific_bridge",
            "separate_open_obligation",
            "separate_open_obligation",
            "Use Xi/Phi analytic structure weaker than the rejected signed/deep hierarchy to prove every Jensen-window hyperbolicity statement.",
            exact["bridge_target"],
            "The strict-Schur counterexample proves that a generic Schur/PF implication is false.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_deep_schur_toda_boundary_gate",
        "date": "2026-07-16",
        "status": (
            "exact rectangular Toda coordinate, moving-tail zero-boundary "
            "obstruction, strict-Schur/Jensen separation, and rejected "
            "all-order endpoint hierarchy; the Xi/Phi Jensen bridge remains open"
        ),
        "proof_boundary": (
            "This artifact proves the rectangular Toda identity, rejects the "
            "ordinary shifted-tail PF shortcut, and gives an exact strictly "
            "Schur-positive sequence whose cubic Jensen polynomial is not "
            "hyperbolic. Rigorous order-ten data reject the proposed m>=10 "
            "endpoint hierarchy at four initial widths. This rejects that "
            "route, not the Xi/Phi Jensen bridge, RH, or Lambda<=0."
        ),
        "sources": sources,
        "source_contract": source_contract,
        "exact": exact,
        "symbolic_toda_audit": toda,
        "boundary_reset_audit": boundary,
        "strict_schur_audit": strict_audit,
        "strict_schur_jensen_counterexample": counterexample,
        "primary_literature_audit": literature,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_to_apply_rows": sum(
                row.readiness == "ready_to_apply" for row in rows
            ),
            "open_endpoint_rows": sum(
                row.readiness == "not_ready_to_apply" for row in rows
            ),
            "rejected_endpoint_rows": sum(
                row.readiness == "rejected_by_counterexample" for row in rows
            ),
            "separate_bridge_rows": sum(
                row.readiness == "separate_open_obligation" for row in rows
            ),
            "symbolic_toda_checks": toda["checks"],
            "boundary_translation_checks": boundary["checks"],
            "deep_translation_checks": boundary["deep_translation_checks"],
            "shallow_boundary_mismatches": boundary["shallow_mismatches"],
            "strict_schur_checks": strict_audit["checks"],
            "exact_strict_schur_jensen_counterexamples": 1,
            "rejected_moving_tail_pf_routes": 1,
            "rejected_generic_schur_jensen_routes": 1,
            "direct_literature_closing_theorems": 0,
            "negative_deep_rectangles": source_contract[
                "negative_deep_rectangles"
            ],
            "open_rectangle_hierarchies": 0,
            "rejected_rectangle_hierarchies": source_contract[
                "rejected_rectangle_hierarchies"
            ],
            "separate_xi_specific_bridges": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_deep_schur_toda_boundary_gate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_deep_schur_toda_boundary_gate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    summary = artifact["summary"]
    boundary = artifact["boundary_reset_audit"]
    counterexample = artifact["strict_schur_jensen_counterexample"]
    lines = [
        "# Jensen-Window PF Deep-Schur Toda And Boundary Gate",
        "",
        "Date: 2026-07-16",
        "",
        "Status: exact rectangular Toda coordinate, moving-tail zero-boundary",
        "obstruction, strict-Schur/Jensen separation, and a rigorous rejection",
        "of the proposed all-order endpoint hierarchy. The order-ten failures",
        "reject that route, not the Xi/Phi Jensen bridge. This is not a proof",
        "of RH or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_deep_schur_toda_boundary_gate.json",
        "python work/rh_compute/scripts/jensen_window_pf_deep_schur_toda_boundary_gate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_deep_schur_toda_boundary_gate.py",
        "```",
        "",
        "## Rectangular Toda Coordinate",
        "",
        "Put",
        "",
        "```text",
        exact["tau_coordinate"],
        "```",
        "",
        "Desnanot-Jacobi becomes the discrete Toda identity",
        "",
        "```text",
        exact["toda_identity"],
        "```",
        "",
        "Consequently, when the lower row is positive,",
        "",
        "```text",
        exact["next_order_criterion"],
        "```",
        "",
        "This is exact but not a propagation theorem. The right side is a",
        "difference of positive terms, so ordinary subtraction-free cluster",
        "positivity does not determine its sign. The proposed strict",
        "width-log-concavity input is false at order ten and widths 9 through 12.",
        "",
        "## Moving-Tail Boundary Check",
        "",
        "A tempting translation is the normalized ordinary tail",
        "",
        "```text",
        exact["ordinary_tail"],
        "```",
        "",
        "But the ordinary tail resets every negative index to zero. The",
        "shifted deep determinant retains `h_(s+k)` whenever `s+k>=0`.",
        "The two Jacobi-Trudi matrices therefore agree only under the same",
        "deep condition one was trying to escape:",
        "",
        "```text",
        exact["conditional_tail_translation"],
        "```",
        "",
        "The smallest exact witness is",
        "",
        "```text",
        exact["boundary_defect"],
        "```",
        "",
        f"The bounded audit checks `{summary['boundary_translation_checks']}` shapes:",
        f"`{summary['deep_translation_checks']}` deep translations agree, while all",
        f"`{summary['shallow_boundary_mismatches']}` shallow shapes have a reset",
        "mismatch in the exact rational test specialization. Thus the proposed",
        "moving-tail PF equivalence is false, and a Littlewood-Richardson lift",
        "from those tails is unavailable.",
        "",
        "## Strict-Schur Jensen Counterexample",
        "",
        "The generic bridge can be rejected even under a much stronger",
        "antecedent. Take",
        "",
        "```text",
        f"H(z)={counterexample['generating_function']}",
        "```",
        "",
        "This is an Edrei specialization with positive exponential parameter",
        "`1/100` and positive denominator parameters `1,2`. More strongly,",
        "the Plancherel component gives, for every partition `lambda`,",
        "",
        "```text",
        artifact["strict_schur_audit"]["all_order_strictness_formula"],
        "```",
        "",
        "so the sequence is strictly Schur-positive for every partition.",
        "Nevertheless its cubic",
        "Jensen polynomial at shift zero is",
        "",
        "```text",
        counterexample["primitive_polynomial"],
        f"discriminant={counterexample['discriminant']['text']}<0",
        "```",
        "",
        "A real cubic with negative discriminant has one real zero and one",
        "nonreal conjugate pair. Hence even strict full Schur positivity of",
        "the unweighted sequence does not imply the binomially weighted Jensen",
        "window is hyperbolic. The bridge must use additional Xi/Phi structure.",
        "",
        "## Primary-Literature Fit",
        "",
    ]
    for row in artifact["primary_literature_audit"]:
        lines.extend(
            [
                f"### {row['id']}",
                "",
                f"Classification: `{row['classification']}`.",
                "",
                row["finding"],
                "",
                f"It does not supply: {row['does_not_supply']}.",
                "",
                f"Source: {row['url']}",
                "",
            ]
        )
    lines.extend(
        [
            "No direct closing theorem was found in this bounded audited set.",
            "The newer Jacobi-Trudi Hadamard result is especially important to",
            "state accurately: its general conjecture remains open, and the",
            "proved 3x2-avoiding case does not supply a replacement for the",
            "rejected rectangle hierarchy.",
            "",
            "## Live Handoffs",
            "",
            "The endpoint hierarchy is rejected:",
            "",
            "```text",
            exact["deep_failure"],
            exact["rejected_rectangle_target"],
            "```",
            "",
            "The logically separate Xi/Phi bridge remains open:",
            "",
            "```text",
            exact["bridge_target"],
            "```",
            "",
            "The next route cannot use bare Toda positivity, ordinary shifted-tail",
            "PF, generic Schur positivity, or the false all-order rectangle",
            "antecedent. It must identify a weaker Xi/Phi-specific mechanism for",
            "the binomially weighted Jensen windows.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote deep-Schur Toda/boundary gate: "
        f"{summary['rows']} rows, "
        f"{summary['symbolic_toda_checks']} symbolic Toda checks, "
        f"{summary['boundary_translation_checks']} boundary checks, "
        f"{summary['strict_schur_checks']} strict-Schur checks, "
        "1 exact full-PF Jensen counterexample, "
        f"{summary['negative_deep_rectangles']} negative order-ten rectangles, "
        "1 rejected rectangle hierarchy"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
