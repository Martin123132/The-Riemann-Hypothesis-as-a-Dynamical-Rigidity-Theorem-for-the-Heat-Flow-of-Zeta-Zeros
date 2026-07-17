#!/usr/bin/env python3
"""Build the normalized deep-Schur coordinate for the static endpoint hierarchy."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path
import re

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_endpoint_deep_schur_coordinate.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_endpoint_deep_schur_coordinate.md"

ALL_ORDER_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_all_order_endpoint_heat_reduction.json"
)
TRANSFER_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_order4_noncontiguous_total_positivity_transfer.json"
)
ENDPOINT_BALL_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k22.jsonl"
)
ENDPOINT_BALL_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k22_summary.json"
)
ORDER10_COUNTEREXAMPLE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_endpoint_order10_counterexample.json"
)
SOURCE_PATHS = (
    ALL_ORDER_SOURCE,
    TRANSFER_SOURCE,
    ENDPOINT_BALL_SOURCE,
    ENDPOINT_BALL_SUMMARY,
    ORDER10_COUNTEREXAMPLE_SOURCE,
)

SYMBOLIC_ORDERS = tuple(range(1, 6))
ARBITRARY_MAX_ORDER = 6
ARBITRARY_COLUMN_BOUND = 8
ARBITRARY_SHIFTS = tuple(range(4))
INVERSE_MAX_ORDER = 5
INVERSE_MAX_PART = 10
RECTANGLE_MAX_ORDER = 64
RECTANGLE_MAX_SHIFT = 63


@dataclass(frozen=True)
class CoordinateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


@dataclass(frozen=True)
class RationalInterval:
    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        if self.lower > self.upper:
            raise ValueError("reversed interval")

    @classmethod
    def ball(cls, center: Fraction, radius: Fraction) -> "RationalInterval":
        if radius < 0:
            raise ValueError("negative radius")
        return cls(center - radius, center + radius)

    def __add__(self, other: "RationalInterval") -> "RationalInterval":
        return RationalInterval(self.lower + other.lower, self.upper + other.upper)

    def __neg__(self) -> "RationalInterval":
        return RationalInterval(-self.upper, -self.lower)

    def __sub__(self, other: "RationalInterval") -> "RationalInterval":
        return self + (-other)

    def __mul__(self, other: "RationalInterval") -> "RationalInterval":
        products = (
            self.lower * other.lower,
            self.lower * other.upper,
            self.upper * other.lower,
            self.upper * other.upper,
        )
        return RationalInterval(min(products), max(products))

    def __truediv__(self, other: "RationalInterval") -> "RationalInterval":
        if other.lower <= 0 <= other.upper:
            raise ZeroDivisionError("interval divisor contains zero")
        reciprocal = RationalInterval(1 / other.upper, 1 / other.lower)
        return self * reciprocal

    def power(self, exponent: int) -> "RationalInterval":
        if exponent < 0:
            raise ValueError("negative interval power")
        result = RationalInterval(Fraction(1), Fraction(1))
        for _ in range(exponent):
            result = result * self
        return result


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_record(path: Path, payload: dict | None = None) -> dict:
    row = {
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "sha256": sha256(path),
    }
    if payload is not None:
        row["kind"] = payload.get("kind")
        row["status"] = payload.get("status")
    return row


def validate_sources() -> tuple[list[dict], dict]:
    all_order = load_json(ALL_ORDER_SOURCE)
    transfer = load_json(TRANSFER_SOURCE)
    endpoint_summary = load_json(ENDPOINT_BALL_SUMMARY)
    order10 = load_json(ORDER10_COUNTEREXAMPLE_SOURCE)

    all_order_summary = all_order.get("summary", {})
    if all_order_summary.get("endpoint_interval_equivalences") != 1:
        raise RuntimeError("all-order endpoint equivalence source is not closed")
    if all_order_summary.get("completed_base_order") != 9:
        raise RuntimeError("completed endpoint base is no longer order nine")
    if transfer.get("summary", {}).get("fixed_order_transfer_theorems") != 1:
        raise RuntimeError("fixed-order arbitrary-column transfer is not closed")
    if endpoint_summary.get("kind") != "acb_coefficient_enclosure_summary":
        raise RuntimeError("endpoint coefficient source has the wrong kind")
    if "-100.0" not in endpoint_summary.get("lambdas", []):
        raise RuntimeError("endpoint coefficient source does not include lambda=-100")
    if endpoint_summary.get("k_min") != 0 or endpoint_summary.get("k_max", -1) < 3:
        raise RuntimeError("endpoint coefficient source does not cover A_0,...,A_3")
    if order10.get("summary", {}).get("negative_deep_rectangles") != 4:
        raise RuntimeError("order-ten deep-rectangle counterexample source changed")
    if order10.get("summary", {}).get("rejected_all_order_endpoint_hierarchies") != 1:
        raise RuntimeError("deep endpoint hierarchy is not marked rejected")

    records = [
        source_record(ALL_ORDER_SOURCE, all_order),
        source_record(TRANSFER_SOURCE, transfer),
        source_record(ENDPOINT_BALL_SOURCE),
        source_record(ENDPOINT_BALL_SUMMARY, endpoint_summary),
        source_record(ORDER10_COUNTEREXAMPLE_SOURCE, order10),
    ]
    contract = {
        "known_endpoint_base": (
            "Q_(m,n)(-100)>0 for 1<=m<=9 and every n>=0"
        ),
        "candidate_endpoint_hierarchy": (
            "Q_(m,n)(-100)>0 for every m>=10 and n>=0"
        ),
        "order10_counterexample": (
            "Q_(10,n)(-100)<0 for n=0,1,2,3"
        ),
        "fixed_order_transfer": (
            "strict signed contiguous layers through m imply all consecutive-row "
            "arbitrary-column signed minors through m"
        ),
        "coefficient_enclosure": (
            "rigorous acb enclosures at lambda=-100 include A_0,...,A_3"
        ),
    }
    return records, contract


def decimal_fraction(text: str) -> Fraction:
    return Fraction(Decimal(text))


def fraction_scientific(value: Fraction, digits: int = 32) -> str:
    with localcontext() as context:
        context.prec = digits + 12
        decimal = Decimal(value.numerator) / Decimal(value.denominator)
        return f"{decimal:.{digits}E}"


BALL_PATTERN = re.compile(r"^\[([^ ]+) \+/- ([^\]]+)\]$")


def endpoint_intervals() -> tuple[dict[int, RationalInterval], list[dict]]:
    intervals: dict[int, RationalInterval] = {}
    rows: list[dict] = []
    with ENDPOINT_BALL_SOURCE.open("r", encoding="utf-8") as handle:
        for line in handle:
            payload = json.loads(line)
            if payload.get("lam") != "-100.0" or payload.get("k") not in range(4):
                continue
            match = BALL_PATTERN.fullmatch(str(payload.get("A_ball", "")))
            if match is None:
                raise RuntimeError("could not parse serialized acb A-ball")
            center = decimal_fraction(match.group(1))
            display_radius = decimal_fraction(match.group(2))
            stored_radius = decimal_fraction(str(payload["A_rad"]))
            radius = max(display_radius, stored_radius)
            interval = RationalInterval.ball(center, radius)
            if interval.lower <= 0:
                raise RuntimeError("endpoint coefficient interval is not strictly positive")
            k = int(payload["k"])
            intervals[k] = interval
            rows.append(
                {
                    "k": k,
                    "center": match.group(1),
                    "conservative_radius": fraction_scientific(radius),
                    "lower": fraction_scientific(interval.lower),
                    "upper": fraction_scientific(interval.upper),
                }
            )
    if set(intervals) != set(range(4)):
        raise RuntimeError("did not recover exactly A_0,...,A_3 at lambda=-100")
    rows.sort(key=lambda row: row["k"])
    return intervals, rows


def endpoint_pf_counterexample() -> dict:
    values, source_rows = endpoint_intervals()
    a0, a1, a2, a3 = (values[index] for index in range(4))
    two = RationalInterval(Fraction(2), Fraction(2))
    determinant = a1.power(3) - two * a0 * a1 * a2 + a0.power(2) * a3
    normalized = determinant / a0.power(3)
    if determinant.upper >= 0 or normalized.upper >= 0:
        raise RuntimeError("endpoint s_(1,1,1) interval is not strictly negative")
    return {
        "lambda": -100,
        "toeplitz_rows": [0, 1, 2],
        "toeplitz_columns": [1, 2, 3],
        "partition": [1, 1, 1],
        "formula_A": "A_1^3-2*A_0*A_1*A_2+A_0^2*A_3",
        "formula_h": "h_1^3-2*h_1*h_2+h_3",
        "unnormalized_interval": {
            "lower": fraction_scientific(determinant.lower),
            "upper": fraction_scientific(determinant.upper),
        },
        "normalized_interval": {
            "lower": fraction_scientific(normalized.lower),
            "upper": fraction_scientific(normalized.upper),
        },
        "strictly_negative": True,
        "excluded_from_deep_order_three": True,
        "exclusion_reason": "lambda_3=1<3-1=2",
        "coefficient_intervals": source_rows,
        "consequence": (
            "the normalized endpoint sequence is not PF_3 and hence not PF-infinity"
        ),
    }


def epsilon(order: int) -> int:
    return (-1) ** (order * (order - 1) // 2)


def arbitrary_partition(n: int, columns: tuple[int, ...]) -> tuple[int, ...]:
    order = len(columns)
    return tuple(
        n + columns[order - 1 - q] + q for q in range(order)
    )


def inverse_columns(partition: tuple[int, ...]) -> tuple[int, tuple[int, ...]]:
    order = len(partition)
    n = partition[-1] - (order - 1)
    columns = tuple(
        partition[order - 1 - ell] - partition[-1] + ell
        for ell in range(order)
    )
    return n, columns


def rectangle_index_audit() -> dict:
    checks = 0
    for order in range(1, RECTANGLE_MAX_ORDER + 1):
        for n in range(RECTANGLE_MAX_SHIFT + 1):
            columns = tuple(range(order))
            partition = arbitrary_partition(n, columns)
            expected = (n + order - 1,) * order
            if partition != expected:
                raise RuntimeError("contiguous columns did not map to a rectangle")
            n_back, columns_back = inverse_columns(partition)
            if n_back != n or columns_back != columns:
                raise RuntimeError("rectangle inverse map failed")
            checks += 1
    return {
        "orders": [1, RECTANGLE_MAX_ORDER],
        "shifts": [0, RECTANGLE_MAX_SHIFT],
        "checks": checks,
        "all_passed": True,
    }


def symbolic_rectangle_audit() -> dict:
    rows = []
    for order in SYMBOLIC_ORDERS:
        n = 2
        values = sp.symbols(f"a0:{n + 2 * order}")
        hankel = sp.Matrix(
            [[values[n + i + j] for j in range(order)] for i in range(order)]
        )
        reversed_transpose = hankel[:, ::-1].T
        rectangle = n + order - 1
        jacobi_trudi = sp.Matrix(
            [
                [values[rectangle - (i + 1) + (j + 1)] for j in range(order)]
                for i in range(order)
            ]
        )
        if reversed_transpose != jacobi_trudi:
            raise RuntimeError(f"symbolic rectangle matrix map failed at order {order}")
        residual = sp.factor(
            jacobi_trudi.det(method="domain-ge")
            - epsilon(order) * hankel.det(method="domain-ge")
        )
        if residual != 0:
            raise RuntimeError(f"symbolic rectangle determinant failed at order {order}")
        rows.append(
            {
                "order": order,
                "shift": n,
                "rectangle_part": rectangle,
                "matrix_entries_identical": True,
                "determinant_residual": str(residual),
            }
        )
    return {
        "orders": list(SYMBOLIC_ORDERS),
        "rows": rows,
        "all_residuals_zero": True,
    }


def arbitrary_column_audit() -> dict:
    checks = 0
    minimum_gaps: set[int] = set()
    for order in range(1, ARBITRARY_MAX_ORDER + 1):
        for n in ARBITRARY_SHIFTS:
            for columns in combinations(range(ARBITRARY_COLUMN_BOUND), order):
                partition = arbitrary_partition(n, columns)
                if any(
                    partition[index] < partition[index + 1]
                    for index in range(order - 1)
                ):
                    raise RuntimeError("column map did not produce a partition")
                if partition[-1] < order - 1:
                    raise RuntimeError("column map escaped the deep cone")
                for row in range(order):
                    for column in range(order):
                        reversed_hankel_index = (
                            n + column + columns[order - 1 - row]
                        )
                        jacobi_trudi_index = partition[row] - row + column
                        if reversed_hankel_index != jacobi_trudi_index:
                            raise RuntimeError("arbitrary-column Jacobi-Trudi map failed")
                n_back, columns_back = inverse_columns(partition)
                if n_back != n + columns[0]:
                    raise RuntimeError("canonical inverse shift was incorrect")
                translated = tuple(value - columns[0] for value in columns)
                if columns_back != translated:
                    raise RuntimeError("canonical inverse columns were incorrect")
                if order > 1:
                    minimum_gaps.add(
                        min(
                            columns[index + 1] - columns[index]
                            for index in range(order - 1)
                        )
                    )
                checks += 1
    return {
        "max_order": ARBITRARY_MAX_ORDER,
        "column_offsets": [0, ARBITRARY_COLUMN_BOUND - 1],
        "shifts": list(ARBITRARY_SHIFTS),
        "checks": checks,
        "observed_minimum_column_gaps": sorted(minimum_gaps),
        "all_passed": True,
        "canonicalization": (
            "the inverse absorbs j_0 into n and returns columns j_l-j_0"
        ),
    }


def bounded_partitions(
    order: int, minimum: int, maximum: int
) -> list[tuple[int, ...]]:
    rows: list[tuple[int, ...]] = []

    def extend(prefix: tuple[int, ...], ceiling: int) -> None:
        if len(prefix) == order:
            rows.append(prefix)
            return
        for value in range(ceiling, minimum - 1, -1):
            extend(prefix + (value,), value)

    extend((), maximum)
    return rows


def inverse_partition_audit() -> dict:
    checks = 0
    by_order = []
    for order in range(1, INVERSE_MAX_ORDER + 1):
        count = 0
        for partition in bounded_partitions(
            order, order - 1, INVERSE_MAX_PART
        ):
            n, columns = inverse_columns(partition)
            if n < 0 or columns[0] != 0:
                raise RuntimeError("deep partition inverse left the canonical domain")
            if any(
                columns[index] >= columns[index + 1]
                for index in range(order - 1)
            ):
                raise RuntimeError("deep partition inverse columns are not increasing")
            if arbitrary_partition(n, columns) != partition:
                raise RuntimeError("deep partition round trip failed")
            count += 1
            checks += 1
        by_order.append({"order": order, "partitions": count})
    return {
        "max_order": INVERSE_MAX_ORDER,
        "max_part": INVERSE_MAX_PART,
        "checks": checks,
        "rows": by_order,
        "all_passed": True,
    }


def exact_statements() -> dict:
    return {
        "normalization": (
            "h_k=A_k(-100)/A_0(-100) for k>=0, h_k=0 for k<0, h_0=1"
        ),
        "hankel_coordinates": (
            "H_(m,n)=det[A_(n+i+j)]_(0<=i,j<m), "
            "epsilon_m=(-1)^binom(m,2), Q_(m,n)=epsilon_m*H_(m,n)"
        ),
        "column_reversal": (
            "Q_(m,n)(-100)=det[A_(n+m-1+i-j)(-100)]_(0<=i,j<m)"
        ),
        "rectangle_identity": (
            "Q_(m,n)(-100)=A_0(-100)^m*s_((n+m-1)^m)(h)"
        ),
        "arbitrary_minor": (
            "R_(m,n)(j_0,...,j_(m-1))="
            "det[A_(n+i+j_l)]_(0<=i,l<m)"
        ),
        "partition_map": (
            "lambda_(q+1)=n+j_(m-1-q)+q for 0<=q<m"
        ),
        "arbitrary_schur_identity": (
            "epsilon_m*R_(m,n)(j_0,...,j_(m-1))="
            "A_0(-100)^m*s_lambda(h)"
        ),
        "inverse_map": (
            "n=lambda_m-(m-1), j_l=lambda_(m-l)-lambda_m+l for 0<=l<m"
        ),
        "deep_cone": (
            "D_m={lambda_1>=...>=lambda_m>=m-1}; "
            "s_lambda(h)=det[h_(lambda_i-i+j)]_(1<=i,j<=m)"
        ),
        "support_boundary": (
            "min_(i,j)(lambda_i-i+j)=lambda_m-m+1>=0 exactly on D_m"
        ),
        "rectangle_target": (
            "s_((N^m))(h)>0 for every m>=10 and N>=m-1"
        ),
        "deep_target": (
            "s_lambda(h)>0 for every m>=10 and lambda in D_m"
        ),
        "endpoint_deep_equivalence": (
            "[Q_(m,n)(-100)>0 for every m>=10,n>=0] iff "
            "[s_lambda(h)>0 for every m>=10,lambda in D_m]"
        ),
        "pf_failure": "s_(1,1,1)(h)<0 at lambda=-100",
        "deep_failure": (
            "s_((N^10))(h)<0 for N=9,10,11,12"
        ),
        "rejected_rectangle_target": (
            "s_((N^m))(h)>0 for every m>=10 and N>=m-1 is false"
        ),
        "dynamic_composition": (
            "deep endpoint positivity => all-order endpoint hierarchy => "
            "all-order signed-Hankel positivity on -100<=lambda<=0"
        ),
        "bridge_boundary": (
            "deep Schur positivity of h_k=A_k(-100)/A_0(-100) is not "
            "Schur positivity of the binomially weighted Jensen-window sequences"
        ),
    }


def literature_audit() -> list[dict]:
    return [
        {
            "id": "gasca_pena_initial_minors",
            "citation": (
                "M. Gasca and J. M. Pena, Total positivity and Neville "
                "elimination, Linear Algebra Appl. 165 (1992), 25-44"
            ),
            "url": "https://doi.org/10.1016/0024-3795(92)90226-Z",
            "classification": "applicable_transfer_only",
            "finding": (
                "Positive initial minors characterize strict total positivity of "
                "each finite reversed block. This is exactly the fixed-order "
                "rectangle-to-arbitrary-shape transfer already used here."
            ),
            "does_not_supply": "new positivity for any open rectangle",
        },
        {
            "id": "edrei_pf_classification",
            "citation": (
                "A. Edrei, Proof of a conjecture of Schoenberg on the generating "
                "function of a totally positive sequence, Canad. J. Math. 5 "
                "(1953), 86-94"
            ),
            "url": "https://doi.org/10.4153/CJM-1953-010-3",
            "classification": "strictly_too_strong",
            "finding": (
                "The theorem classifies full Toeplitz total positivity, whereas "
                "the endpoint specialization has a rigorously negative order-three "
                "Toeplitz minor outside the deep cone."
            ),
            "does_not_supply": "restricted deep-cone positivity",
        },
        {
            "id": "pena_hankel_toeplitz_orientation",
            "citation": (
                "J. M. Pena, Positive Hankel Matrices, Eigenvalues and Total "
                "Positivity, Mathematics 13 (2025), 2278, Proposition 6"
            ),
            "url": "https://doi.org/10.3390/math13142278",
            "classification": "structural_orientation_support",
            "finding": (
                "Finite Hankel column reversal produces the signature "
                "(-1)^binom(k,2) in the associated Toeplitz matrix. It supports "
                "the orientation coordinate but assumes total positivity on the "
                "Hankel side and does not prove the Xi endpoint signs."
            ),
            "does_not_supply": "the open all-order endpoint hierarchy",
        },
        {
            "id": "kushel_eventual_total_positivity",
            "citation": (
                "O. Y. Kushel, Matrices with totally positive powers and their "
                "generalizations, Operators and Matrices 9 (2015), 943-964"
            ),
            "url": "https://doi.org/10.7153/oam-09-56",
            "classification": "coordinate_mismatch",
            "finding": (
                "Eventual total positivity there means total positivity of all "
                "sufficiently high matrix powers, not positivity of Toeplitz minors "
                "whose Jacobi-Trudi shapes lie beyond a moving depth boundary."
            ),
            "does_not_supply": "an eventual-depth Toeplitz theorem",
        },
    ]


def build_artifact() -> dict:
    sources, source_contract = validate_sources()
    exact = exact_statements()
    rectangle_audit = rectangle_index_audit()
    symbolic_audit = symbolic_rectangle_audit()
    arbitrary_audit = arbitrary_column_audit()
    inverse_audit = inverse_partition_audit()
    counterexample = endpoint_pf_counterexample()
    literature = literature_audit()

    rows = [
        CoordinateRow(
            "jwpfedsc_01_normalization",
            "exact_definition",
            "ready_to_apply",
            "Normalize the endpoint sequence so it is a genuine complete-homogeneous specialization with h_0=1.",
            exact["normalization"],
            "A_0(-100)>0 is rigorously enclosed away from zero.",
            counterexample["coefficient_intervals"][0],
        ),
        CoordinateRow(
            "jwpfedsc_02_column_reversal",
            "exact_identity",
            "ready_to_apply",
            "Reversing the contiguous Hankel columns absorbs the signed orientation and produces a solid Toeplitz determinant.",
            exact["column_reversal"],
            "One finite column reversal and one transpose; no sign hypothesis is used.",
            symbolic_audit,
        ),
        CoordinateRow(
            "jwpfedsc_03_rectangles",
            "exact_jacobi_trudi_identity",
            "ready_to_apply",
            "Every signed contiguous endpoint minor is exactly a normalized rectangular Schur value.",
            exact["rectangle_identity"],
            "The rectangle has m parts N=n+m-1, hence N>=m-1.",
            rectangle_audit,
        ),
        CoordinateRow(
            "jwpfedsc_04_arbitrary_columns",
            "exact_jacobi_trudi_identity",
            "ready_to_apply",
            "Every consecutive-row arbitrary-column signed minor maps to one deep partition.",
            exact["partition_map"] + "; " + exact["arbitrary_schur_identity"],
            "Increasing column gaps are exactly the adjacent partition inequalities.",
            arbitrary_audit,
        ),
        CoordinateRow(
            "jwpfedsc_05_inverse_map",
            "exact_bijection",
            "ready_to_apply",
            "Every partition with m parts and smallest part at least m-1 comes from a unique canonical shift and column set with j_0=0.",
            exact["inverse_map"],
            "Noncanonical representations differ only by moving j_0 into n.",
            inverse_audit,
        ),
        CoordinateRow(
            "jwpfedsc_06_deep_support_boundary",
            "exact_coordinate_boundary",
            "ready_to_apply",
            "The depth inequality is exactly the condition that the whole Jacobi-Trudi matrix avoids negative-index structural zeros.",
            exact["deep_cone"] + "; " + exact["support_boundary"],
            "This defines the term deep; it does not assert positivity.",
        ),
        CoordinateRow(
            "jwpfedsc_07_endpoint_equivalence",
            "exact_all_order_reduction",
            "ready_to_apply",
            "With the completed lower base and fixed-order transfer, the candidate endpoint hierarchy is equivalent to positivity on the entire deep Schur cone.",
            exact["endpoint_deep_equivalence"],
            "The equivalence is exact, but the order-ten counterexample makes both universal positivity statements false.",
            {"source_contract": source_contract},
        ),
        CoordinateRow(
            "jwpfedsc_08_minimal_rectangle_target",
            "exact_target_coordinate",
            "ready_to_apply",
            "Within the candidate hierarchy, only rectangular deep shapes are independent; the remaining deep shapes would follow structurally from lower rectangles.",
            exact["rectangle_target"],
            "An exact coordinate for a candidate statement later rejected at order ten.",
        ),
        CoordinateRow(
            "jwpfedsc_09_endpoint_pf_counterexample",
            "rigorous_exclusion_certificate",
            "ready_to_apply",
            "The normalized endpoint sequence fails ordinary PF already at the first excluded order-three solid block.",
            exact["pf_failure"],
            "Exact rational interval propagation from rigorous acb coefficient enclosures.",
            counterexample,
        ),
        CoordinateRow(
            "jwpfedsc_10_strictly_weaker_than_pf",
            "route_separation_theorem",
            "ready_to_apply",
            "The deep cone is a proper subfamily of the full Schur/PF inequality family, and an additional PF inequality already fails for the actual endpoint sequence.",
            "s_(1,1,1)(h)<0 but (1,1,1) is not in D_3 because 1<2",
            "This blocks Edrei/ASW as a direct endpoint theorem; the separate order-ten certificate now also rejects the deep cone itself.",
        ),
        CoordinateRow(
            "jwpfedsc_11_literature_fit",
            "primary_literature_gate",
            "ready_to_apply",
            "The audited primary theorems support the coordinate and transfer or impose stronger incompatible hypotheses; none supplies the candidate rectangle signs, which are now directly false at order ten.",
            "direct_closing_theorems_in_audited_set=0",
            "A bounded theorem search, not a claim that no relevant theorem exists anywhere.",
            literature,
        ),
        CoordinateRow(
            "jwpfedsc_12_dynamic_composition",
            "exact_conditional_consequence",
            "ready_to_apply",
            "Abstract static deep endpoint positivity would feed into the completed cooperative heat induction.",
            exact["dynamic_composition"],
            "The implication remains exact, but its antecedent is false for the actual endpoint sequence.",
            {"source": sources[0]},
        ),
        CoordinateRow(
            "jwpfedsc_13_rejected_deep_target",
            "countermodel_gate",
            "rejected_by_counterexample",
            "The first four order-ten rectangles are negative, so the proposed deep endpoint hierarchy is false.",
            exact["deep_failure"] + "; " + exact["rejected_rectangle_target"],
            "Rejects the all-order deep-cone antecedent, not RH or Jensen hyperbolicity.",
            {"source": sources[4]},
        ),
        CoordinateRow(
            "jwpfedsc_14_jensen_bridge_boundary",
            "route_separation_guard",
            "separate_open_obligation",
            "The binomially weighted Jensen-window problem survives, but its Xi/Phi-specific antecedent must be weaker than the rejected deep cone.",
            exact["bridge_boundary"],
            "The endpoint specialization and every finite Jensen-window specialization are different Toeplitz sequences.",
        ),
    ]

    return {
        "kind": "jensen_window_pf_endpoint_deep_schur_coordinate",
        "date": "2026-07-16",
        "status": (
            "exact normalized deep-Schur coordinate and rigorous endpoint PF "
            "separation, with the all-order rectangle hierarchy rejected at "
            "order ten"
        ),
        "proof_boundary": (
            "This artifact proves the Jacobi-Trudi index maps, their inverse, the "
            "equivalence with the existing endpoint hierarchy, and a rigorous "
            "negative PF_3 minor outside the deep cone. The order-ten source "
            "rejects the deep endpoint hierarchy itself. This does not settle "
            "the separate Jensen-window PF bridge, Jensen hyperbolicity, RH, or "
            "Lambda<=0."
        ),
        "sources": sources,
        "source_contract": source_contract,
        "exact": exact,
        "rectangle_index_audit": rectangle_audit,
        "symbolic_rectangle_audit": symbolic_audit,
        "arbitrary_column_audit": arbitrary_audit,
        "inverse_partition_audit": inverse_audit,
        "endpoint_pf_counterexample": counterexample,
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
            "symbolic_orders": len(SYMBOLIC_ORDERS),
            "rectangle_index_checks": rectangle_audit["checks"],
            "arbitrary_column_checks": arbitrary_audit["checks"],
            "inverse_partition_checks": inverse_audit["checks"],
            "deep_cone_bijections": 1,
            "endpoint_deep_equivalences": 1,
            "rigorous_endpoint_pf_counterexamples": 1,
            "direct_literature_closing_theorems": 0,
            "completed_base_order": 9,
            "first_failed_order": 10,
            "negative_deep_rectangles": 4,
            "open_rectangle_hierarchies": 0,
            "rejected_rectangle_hierarchies": 1,
            "separate_jensen_pf_bridges": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_endpoint_deep_schur_coordinate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_endpoint_deep_schur_coordinate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    summary = artifact["summary"]
    counterexample = artifact["endpoint_pf_counterexample"]
    literature = artifact["primary_literature_audit"]
    lines = [
        "# Jensen-Window PF Endpoint Deep-Schur Coordinate",
        "",
        "Date: 2026-07-16",
        "",
        "Status: exact normalized deep-Schur coordinate and rigorous endpoint",
        "PF separation. The proposed all-order deep rectangle hierarchy is",
        "rejected by an order-ten endpoint counterexample. This does not settle",
        "the separate Jensen-window PF bridge. This is not a proof of RH or",
        "`Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_endpoint_deep_schur_coordinate.json",
        "python work/rh_compute/scripts/jensen_window_pf_endpoint_deep_schur_coordinate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_endpoint_deep_schur_coordinate.py",
        "```",
        "",
        "## Normalized Endpoint Specialization",
        "",
        "Use the genuine complete-homogeneous specialization",
        "",
        "```text",
        exact["normalization"],
        "```",
        "",
        "The normalization is essential: ordinary Schur theory has `h_0=1`.",
        "Every size-`m` determinant in the unnormalized `A_k` acquires only the",
        "strictly positive common factor `A_0(-100)^m`.",
        "",
        "## Rectangular Identity",
        "",
        "Reverse the `m` Hankel columns. Since the reversal has sign",
        "`epsilon_m=(-1)^binom(m,2)`, transpose the resulting matrix, and put",
        "`N=n+m-1`. The Jacobi-Trudi matrix is obtained entry by entry:",
        "",
        "```text",
        exact["column_reversal"],
        exact["rectangle_identity"],
        "```",
        "",
        "Thus the candidate static endpoint statement was precisely",
        "",
        "```text",
        exact["rectangle_target"],
        "```",
        "",
        "This is an exact change of coordinates; the counterexample below",
        "shows that the universal positivity statement is false.",
        "",
        "## Arbitrary Columns And The Deep Cone",
        "",
        "For increasing offsets `0<=j_0<...<j_(m-1)`, define `R_(m,n)` as in",
        "the existing arbitrary-column transfer theorem. Reversal and transpose",
        "give",
        "",
        "```text",
        exact["partition_map"],
        exact["arbitrary_schur_identity"],
        "```",
        "",
        "The adjacent partition differences are",
        "`lambda_r-lambda_(r+1)=j_(m-r)-j_(m-r-1)-1>=0`, and",
        "`lambda_m=n+j_0+m-1>=m-1`. Conversely, for every",
        "`lambda_1>=...>=lambda_m>=m-1`, the canonical inverse is",
        "",
        "```text",
        exact["inverse_map"],
        "```",
        "",
        "It has `j_0=0` and strictly increasing columns. Hence the arbitrary",
        "signed minors are in bijection with",
        "",
        "```text",
        exact["deep_cone"],
        "```",
        "",
        "The threshold is not decorative:",
        "",
        "```text",
        exact["support_boundary"],
        "```",
        "",
        "So the deep cone is exactly the Jacobi-Trudi region that never touches",
        "the artificial values `h_k=0` for `k<0`.",
        "",
        "## Exact Endpoint Equivalence",
        "",
        "The completed orders through nine provide every lower initial minor.",
        "At each fixed higher order, the Gasca-Pena initial-minor theorem then",
        "transfers rectangular positivity to every arbitrary column set. Since",
        "rectangles are themselves deep shapes, this proves",
        "",
        "```text",
        exact["endpoint_deep_equivalence"],
        "```",
        "",
        "Within this candidate hierarchy, only the rectangles are independent;",
        "the rest of the deep cone would be a structural consequence once all",
        "lower rectangles were available.",
        "",
        "## Rigorous PF Separation",
        "",
        "The rigorous acb endpoint enclosures for `A_0,...,A_3` were propagated",
        "with exact rational interval arithmetic. They give",
        "",
        "```text",
        "s_(1,1,1)(h)=h_1^3-2*h_1*h_2+h_3",
        f"interval = [{counterexample['normalized_interval']['lower']},",
        f"            {counterexample['normalized_interval']['upper']} ]",
        "```",
        "",
        "The upper endpoint is strictly negative. Therefore the actual normalized",
        "endpoint sequence is not `PF_3`, hence not `PF-infinity`. But",
        "`(1,1,1)` lies outside `D_3` because its smallest part is `1<2`.",
        "Thus full Schur positivity imposes a genuinely additional inequality",
        "that is already false here. The ordinary Edrei route is unavailable,",
        "while this failed inequality is not part of the rectangular endpoint",
        "family.",
        "",
        "## Deep Rectangle Counterexample",
        "",
        "The first-open-order certificate now gives required deep failures:",
        "",
        "```text",
        exact["deep_failure"],
        exact["rejected_rectangle_target"],
        "```",
        "",
        "These are not shallow boundary shapes. At order ten the depth condition",
        "is `N>=9`, so `(9^10)` through `(12^10)` belong to the required cone.",
        "",
        "## Primary-Literature Fit Gate",
        "",
    ]
    for row in literature:
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
            "No direct closing theorem was found in this audited primary-source",
            "set. That is a bounded route audit, not a claim about all literature.",
            "The phrase `eventual total positivity` is especially unsafe here:",
            "in the audited matrix literature it refers to high matrix powers,",
            "whereas our coordinate concerns minors beyond a moving Toeplitz depth.",
            "",
            "## Conditional Heat Composition",
            "",
            "The abstract implication into the heat theorem remains exact:",
            "",
            "```text",
            exact["dynamic_composition"],
            "```",
            "",
            "For the actual endpoint sequence its antecedent is false, so this",
            "composition is unavailable as an RH route.",
            "",
            "The separate Jensen-window obligation remains:",
            "",
            "```text",
            exact["bridge_boundary"],
            "```",
            "",
            "## Machine Audit",
            "",
            f"The generator records `{summary['rectangle_index_checks']}` rectangle-map checks,",
            f"`{summary['arbitrary_column_checks']}` arbitrary-column checks,",
            f"`{summary['inverse_partition_checks']}` bounded inverse-partition checks, and",
            f"symbolic determinant residuals through order `{max(SYMBOLIC_ORDERS)}`.",
            "The finite checks audit implementation. General validity comes from",
            "the explicit entrywise maps and finite determinant identities above.",
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
        "wrote endpoint deep-Schur coordinate: "
        f"{summary['rows']} rows, "
        f"{summary['arbitrary_column_checks']} arbitrary-column checks, "
        f"{summary['inverse_partition_checks']} inverse checks, "
        "1 rigorous PF counterexample, 4 negative deep rectangles, "
        "1 rejected rectangle hierarchy"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
