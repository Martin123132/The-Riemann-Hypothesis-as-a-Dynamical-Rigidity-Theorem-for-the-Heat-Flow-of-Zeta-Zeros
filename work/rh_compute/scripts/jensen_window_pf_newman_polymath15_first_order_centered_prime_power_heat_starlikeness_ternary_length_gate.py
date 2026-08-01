#!/usr/bin/env python3
"""Build the ternary all-length prime-power heat-starlikeness gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction
from hashlib import sha256
import json
from math import comb
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_ternary_length_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
LENGTH_STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_length_propagation_gate"
)
BASE_STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_base_certificate"
)
PHASE_STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_logarithmic_phase_flow_gate"
)
SOURCES = {
    "length_propagation": (
        REPO_ROOT / "work/rh_compute/results" / f"{LENGTH_STEM}.json"
    ),
    "base_certificate": (
        REPO_ROOT / "work/rh_compute/results" / f"{BASE_STEM}.json"
    ),
    "phase_flow": (
        REPO_ROOT / "work/rh_compute/results" / f"{PHASE_STEM}.json"
    ),
}

Q_LOWER = Fraction(17, 20)
SQRT_THREE_RECIPROCAL_UPPER = Fraction(577351, 1_000_000)
WIDE_A_UPPER = Fraction(29, 50)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def decimal_text(value: Fraction, digits: int = 30) -> str:
    with localcontext() as context:
        context.prec = digits
        return format(
            Decimal(value.numerator) / Decimal(value.denominator), ".20e"
        )


def add_term(
    terms: dict[tuple[int, int], Fraction],
    q_degree: int,
    a_degree: int,
    coefficient: int | Fraction,
) -> None:
    key = (q_degree, a_degree)
    terms[key] = terms.get(key, Fraction(0)) + Fraction(coefficient)
    if not terms[key]:
        del terms[key]


def q_monomial_bernstein(
    monomial_degree: int, elevated_degree: int
) -> list[Fraction]:
    denominator = comb(elevated_degree, monomial_degree)
    powers = [Q_LOWER**k for k in range(monomial_degree + 1)]
    values: list[Fraction] = []
    for index in range(elevated_degree + 1):
        lower = max(0, monomial_degree - (elevated_degree - index))
        upper = min(monomial_degree, index)
        numerator = sum(
            Fraction(
                comb(index, j)
                * comb(elevated_degree - index, monomial_degree - j)
            )
            * powers[monomial_degree - j]
            for j in range(lower, upper + 1)
        )
        values.append(numerator / denominator)
    return values


def second_monomial_bernstein(
    monomial_degree: int,
    elevated_degree: int,
    upper: Fraction,
) -> list[Fraction]:
    denominator = comb(elevated_degree, monomial_degree)
    scale = upper**monomial_degree
    return [
        (
            scale * Fraction(comb(index, monomial_degree), denominator)
            if index >= monomial_degree
            else Fraction(0)
        )
        for index in range(elevated_degree + 1)
    ]


def bernstein_certificate(
    source: dict[tuple[int, int], Fraction],
    second_upper: Fraction,
    claimed_lower: Fraction,
) -> dict:
    terms = dict(source)
    add_term(terms, 0, 0, -claimed_lower)
    q_degree = max(q_degree for q_degree, _ in terms)
    second_degree = max(second_degree for _, second_degree in terms)
    q_vectors = {
        degree: q_monomial_bernstein(degree, q_degree)
        for degree in {q_degree for q_degree, _ in terms}
    }
    second_vectors = {
        degree: second_monomial_bernstein(
            degree, second_degree, second_upper
        )
        for degree in {second_degree for _, second_degree in terms}
    }
    minimum: Fraction | None = None
    argminimum = (-1, -1)
    for q_index in range(q_degree + 1):
        for second_index in range(second_degree + 1):
            value = sum(
                coefficient
                * q_vectors[q_power][q_index]
                * second_vectors[second_power][second_index]
                for (q_power, second_power), coefficient in terms.items()
            )
            if minimum is None or value < minimum:
                minimum = value
                argminimum = (q_index, second_index)
    if minimum is None:
        raise RuntimeError("empty Bernstein polynomial")
    return {
        "q_interval": ["17/20", "1"],
        "second_interval": ["0", fraction_text(second_upper)],
        "power_degrees": [q_degree, second_degree],
        "bernstein_shape": [q_degree + 1, second_degree + 1],
        "claimed_lower": fraction_text(claimed_lower),
        "minimum_shifted_coefficient": fraction_text(minimum),
        "minimum_shifted_coefficient_decimal": decimal_text(minimum),
        "argminimum": list(argminimum),
        "positive": minimum > 0,
    }


def nonterminal_cluster_terms(d: int) -> dict[tuple[int, int], Fraction]:
    # This is C(d+2,d) in the variables (q,a).
    B = d + 2
    terms: dict[tuple[int, int], Fraction] = {}
    add_term(terms, 0, 0, B)
    add_term(terms, 3, 1, -2 * (B + 1))
    add_term(terms, 4, 2, B + 2)
    add_term(terms, 2 * d + 6, 2, B + 2)
    add_term(terms, 2 * d + 7, 3, -2 * (B + 3))
    add_term(terms, 4 * d + 8, 4, B + 4)
    return terms


def zero_reserve_terms() -> dict[tuple[int, int], Fraction]:
    # The second variable is v and b=alpha*q^3*v.
    terms: dict[tuple[int, int], Fraction] = {}
    alpha = SQRT_THREE_RECIPROCAL_UPPER
    for k in range(4):
        A = 2 * k + 2
        c_q_degree = 14 * k - 2 * k * k
        c_b_degree = 2 * k
        for coefficient, extra_b, extra_q in (
            (A, 0, 0),
            (-2 * (A + 1), 1, 6 - 2 * k),
            (A + 2, 2, 10 - 4 * k),
        ):
            b_degree = c_b_degree + extra_b
            q_degree = c_q_degree + extra_q + 3 * b_degree
            add_term(
                terms,
                q_degree,
                b_degree,
                Fraction(coefficient) * alpha**b_degree,
            )
    return terms


def d0_terms(M: int) -> dict[tuple[int, int], Fraction]:
    exponents = [k * (2 * M - 2 - k) for k in range(M)]
    terms: dict[tuple[int, int], Fraction] = {}
    for k in range(M):
        add_term(terms, 2 * exponents[k], 2 * k, 2 * (k + 1))
    for k in range(M - 1):
        add_term(
            terms,
            exponents[k] + exponents[k + 1],
            2 * k + 1,
            -2 * (2 * k + 3),
        )
    for k in range(M - 2):
        add_term(
            terms,
            exponents[k] + exponents[k + 2],
            2 * k + 2,
            2 * k + 4,
        )
    return terms


def terminal_terms(M: int) -> dict[tuple[int, int], Fraction]:
    terms: dict[tuple[int, int], Fraction] = {}
    base_q = M * (M - 2)
    base_a = M - 2
    add_term(terms, base_q, base_a, M)
    add_term(terms, base_q + 1, base_a + 1, -2 * (M + 1))
    add_term(
        terms,
        base_q + 2 * M - 2,
        base_a + 2,
        M + 2,
    )
    return terms


def combine_terms(
    left: dict[tuple[int, int], Fraction],
    right: dict[tuple[int, int], Fraction],
    scale: int,
) -> dict[tuple[int, int], Fraction]:
    result = dict(left)
    for (q_degree, a_degree), coefficient in right.items():
        add_term(result, q_degree, a_degree, scale * coefficient)
    return result


def build_certificates() -> dict:
    nonterminal = []
    for d in range(1, 6):
        certificate = bernstein_certificate(
            nonterminal_cluster_terms(d),
            WIDE_A_UPPER,
            Fraction(17, 100),
        )
        certificate["d"] = d
        nonterminal.append(certificate)

    zero_reserve = bernstein_certificate(
        zero_reserve_terms(),
        Fraction(1),
        Fraction(11, 250),
    )
    m5_d0 = bernstein_certificate(
        d0_terms(5),
        SQRT_THREE_RECIPROCAL_UPPER,
        Fraction(1, 25),
    )

    finite_energy = []
    for M in range(5, 17):
        energy = combine_terms(
            d0_terms(M), terminal_terms(M), (M - 1) ** 2
        )
        certificate = bernstein_certificate(
            energy,
            SQRT_THREE_RECIPROCAL_UPPER,
            Fraction(1, 25),
        )
        certificate["M"] = M
        finite_energy.append(certificate)

    certificates = {
        "nonterminal_clusters_d1_to_d5": nonterminal,
        "four_group_zero_reserve": zero_reserve,
        "m5_zero_reserve": m5_d0,
        "finite_energy_m5_to_m16": finite_energy,
    }
    if not all(
        row["positive"]
        for rows in (
            nonterminal,
            [zero_reserve],
            [m5_d0],
            finite_energy,
        )
        for row in rows
    ):
        raise RuntimeError("a ternary Bernstein certificate is not positive")
    return certificates


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCES.items()
    }


def audit_sources(payloads: dict[str, dict]) -> dict:
    length_summary = payloads["length_propagation"].get("summary", {})
    if length_summary.get("analytic_all_length_large_prime_families") != 1:
        raise RuntimeError("length-propagation source drifted")
    if payloads["base_certificate"].get("summary", {}).get(
        "exact_bernstein_base_certificates"
    ) != 3:
        raise RuntimeError("base-certificate source drifted")
    actual = payloads["base_certificate"].get("exact", {}).get(
        "actual_transfer", {}
    )
    if actual.get("coefficient_error") != "|R_k-1|<r_x=18000/x":
        raise RuntimeError("actual-transfer source drifted")
    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCES.items()
        },
        "source_kinds": {
            key: payload.get("kind") for key, payload in payloads.items()
        },
    }


def build_rows() -> list[GateRow]:
    return [
        GateRow(
            "pphstl_01_domain",
            "ternary complete-chain domain",
            "proved",
            "Fix p=3 and every complete-chain length M>=5.",
            "17/20<=q<=1; 73/100<=s<=1; a=s/sqrt(3)<=1/sqrt(3).",
            "Lengths M=2,3 and all singleton or joined chains remain separate.",
        ),
        GateRow(
            "pphstl_02_fourier",
            "Fourier current",
            "proved",
            "Use the exact Fourier coefficients and Fejer decomposition from the preceding length gate.",
            "J_0=A_0/2+sum_(d=1)^(M-1)A_d cos(d theta); D_d=A_d-2A_(d+1)+A_(d+2).",
            "Convexity is proved below rather than assumed.",
        ),
        GateRow(
            "pphstl_03_common_groups",
            "common-summand positivity",
            "proved",
            "For d>=1 every common interior group is positive.",
            "A-2(A+1)r+(A+2)rr_+ >=(1-r)[A-(A+2)r]>0.",
            "This uses A>=3 and 1/sqrt(3)<3/5<=A/(A+2).",
        ),
        GateRow(
            "pphstl_04_terminal_cluster",
            "terminal pairing",
            "proved",
            "For 1<=d<=M-3, pair the last common group with the two boundary terms.",
            "C(B,d)=I_B+a^2 q^(2d+6)G_(B+2,d); B=2M-d-4>=d+2.",
            "The unpaired index d=M-2 is retained explicitly.",
        ),
        GateRow(
            "pphstl_05_cluster_slope",
            "cluster monotonicity",
            "proved",
            "C(B,d) is strictly increasing in B.",
            "partial_B C >=(5-8/sqrt(3))/3>0.",
            "The last comparison is exact because 75>64.",
        ),
        GateRow(
            "pphstl_06_small_d",
            "finite nonterminal boundary",
            "proved",
            "Exact ordinary Bernstein coefficients prove C(d+2,d)>17/100 for d=1,...,5.",
            "q in [17/20,1], a in [0,29/50].",
            "The certificate covers a larger a-box than the physical ternary domain.",
        ),
        GateRow(
            "pphstl_07_large_d",
            "analytic nonterminal tail",
            "proved",
            "The paired terminal cluster is positive for every d>=6.",
            "3C(d+2,d)>=5d+14-(8d+28)/sqrt(3)>0.",
            "The d=6 endpoint is 44-76/sqrt(3)>0; the slope uses 75>64.",
        ),
        GateRow(
            "pphstl_08_unique_defect",
            "defect localization",
            "proved",
            "All Fourier second differences except possibly D_(M-2) are positive.",
            "D_d>0 for 1<=d<=M-3; D_(M-1)=A_(M-1)>0.",
            "No sign is assigned to D_(M-2).",
        ),
        GateRow(
            "pphstl_09_zero_reserve",
            "uniform first reserve",
            "proved",
            "The first four common groups give D_0>11/250 for M>=6; M=5 has D_0>1/25.",
            "b=a q^(2M-9)<=alpha q^3; alpha=577351/10^6>1/sqrt(3).",
            "Both lower bounds are exact rational Bernstein certificates.",
        ),
        GateRow(
            "pphstl_10_finite_energy",
            "finite defect absorption",
            "proved",
            "For 5<=M<=16 exact Bernstein arithmetic absorbs the possible terminal defect.",
            "D_0+(M-1)^2 D_(M-2)>1/25.",
            "The certificate uses the enlarged box 0<=a<=alpha.",
        ),
        GateRow(
            "pphstl_11_tail_implication",
            "negative-defect localization",
            "proved",
            "If D_(M-2)<0, its q-power is forced below a strict scalar threshold.",
            "q^(2M-2)<R_M=3[2(M+1)/sqrt(3)-M]/(M+2).",
            "This follows from the retained positive terminal square.",
        ),
        GateRow(
            "pphstl_12_tail_decay",
            "analytic defect tail",
            "proved",
            "For every M>=17 the weighted negative defect has magnitude below 1/100.",
            "(M-1)^2(-D_(M-2))<((M-1)^2 M/3)5^(-(M-2)/2)<1/100.",
            "R_M<3/5 and the displayed majorant decreases from M=17.",
        ),
        GateRow(
            "pphstl_13_fejer_absorption",
            "Fejer lower bound",
            "proved",
            "Use 0<=F_(M-2)<=M-1 and discard only proved nonnegative terms.",
            "J_0>=1/2[D_0+(M-1)^2D_(M-2)] when D_(M-2)<0.",
            "When the defect is nonnegative, J_0>=D_0/2.",
        ),
        GateRow(
            "pphstl_14_ideal_theorem",
            "all-length ternary theorem",
            "proved",
            "Every correction-free ternary complete chain of length M>=5 is strictly heat-starlike.",
            "p=3, M>=5: J_0>17/1000.",
            "This is an all-length theorem, not a finite scan.",
        ),
        GateRow(
            "pphstl_15_actual_sums",
            "length-uniform actual transfer",
            "proved",
            "Geometric coefficient sums make the actual perturbation uniform in M.",
            "sum c_k<5/2; sum k c_k<15/4; |J_ang-J_0|<864000/x.",
            "The relative coefficient estimate remains |R_k-1|<18000/x.",
        ),
        GateRow(
            "pphstl_16_actual_current",
            "actual fixed-ray theorem",
            "proved",
            "The rho and epsilon corrections cost less than the retained ideal margin.",
            "|J_ray-J_0|<864020/x+633450/x^2<1/1000; J_ray>2/125.",
            "This is local to one complete ternary prime-power ray.",
        ),
        GateRow(
            "pphstl_17_pi_provenance",
            "constant provenance",
            "proved",
            "No new circle constant enters the ternary decomposition.",
            "Any pi in x=4*pi*exp(L) is inherited from the completed-zeta saddle normalization.",
            "The Fejer kernels use the ordinary 2*pi period of exp(i theta).",
        ),
        GateRow(
            "pphstl_18_boundary",
            "proof boundary",
            "open",
            "The next all-length complete-chain family is dyadic M>=10.",
            "Open: p=2,M>=10; short/singleton chains; p-free, endpoint, cutoff, adjacent, Abel, and winding joins.",
            "No contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion is claimed.",
        ),
    ]


def build_payload() -> dict:
    source_payloads = load_sources()
    certificates = build_certificates()
    rows = build_rows()
    return {
        "schema_version": 1,
        "kind": "prime_power_heat_starlikeness_ternary_length_gate",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "title": (
            "Ternary all-length prime-power heat-starlikeness gate"
        ),
        "source_audit": audit_sources(source_payloads),
        "assumptions": {
            "prime": 3,
            "block_lengths": "M>=5",
            "q_interval": ["17/20", "1"],
            "s_interval": ["73/100", "1"],
            "a_definition": "a=s/sqrt(3)",
            "a_upper_rational": "577351/1000000",
            "a_upper_check": (
                "3*(577351)^2-10^12=2531603>0"
            ),
        },
        "rows": [asdict(row) for row in rows],
        "exact": {
            "fourier": {
                "coefficients": (
                    "A_0=2sum_(k=0)^(M-1)(k+1)c_k^2; "
                    "A_d=sum_(k=0)^(M-1-d)(2k+d+2)c_kc_(k+d)"
                ),
                "difference": "D_d=A_d-2A_(d+1)+A_(d+2)",
                "fejer": (
                    "J_0=(1/2)sum_(d=0)^(M-1)(d+1)D_dF_d; "
                    "0<=F_d<=d+1"
                ),
            },
            "nonterminal": {
                "interior_group": (
                    "A-2(A+1)r+(A+2)rr_+"
                    ">=(1-r)[A-(A+2)r]"
                ),
                "cluster": (
                    "C(B,d)=B-2(B+1)aq^3+(B+2)a^2q^4"
                    "+a^2q^(2d+6)[B+2-2(B+3)aq"
                    "+(B+4)a^2q^(2d+2)]"
                ),
                "slope_lower": "(5-8/sqrt(3))/3>0",
                "large_d_lower": (
                    "3C(d+2,d)>=5d+14-(8d+28)/sqrt(3)>0, d>=6"
                ),
                "certificates": certificates[
                    "nonterminal_clusters_d1_to_d5"
                ],
            },
            "zero_reserve": {
                "alpha": "577351/1000000",
                "alpha_check": "3alpha^2-1=2531603/10^12>0",
                "four_group_certificate": certificates[
                    "four_group_zero_reserve"
                ],
                "m5_certificate": certificates["m5_zero_reserve"],
                "theorem": (
                    "D_0>11/250 for M>=6; D_0>1/25 for M=5"
                ),
            },
            "terminal_defect": {
                "formula": (
                    "D_(M-2)=a^(M-2)q^[M(M-2)]"
                    "[M-2(M+1)aq+(M+2)a^2q^(2M-2)]"
                ),
                "finite_energy": (
                    "D_0+(M-1)^2D_(M-2)>1/25, 5<=M<=16"
                ),
                "finite_certificates": certificates[
                    "finite_energy_m5_to_m16"
                ],
                "negative_implication": (
                    "D_(M-2)<0 implies "
                    "q^(2M-2)<R_M=3[2(M+1)/sqrt(3)-M]/(M+2)"
                ),
                "tail_bound": (
                    "M>=17: (M-1)^2(-D_(M-2))<1/100"
                ),
            },
            "ideal_theorem": {
                "finite": "5<=M<=16: J_0>1/50",
                "tail": "M>=17: J_0>17/1000",
                "uniform": "p=3,M>=5: J_0>17/1000",
            },
            "actual_transfer": {
                "coefficient_error": "|R_k-1|<r_x=18000/x",
                "geometric_sums": (
                    "sum c_k<5/2; sum k c_k<15/4; "
                    "absolute current mass<16"
                ),
                "frozen_error": "|J_ang-J_0|<48r_x=864000/x",
                "ray_error": (
                    "|J_ray-J_0|<864020/x+633450/x^2<1/1000"
                ),
                "actual_theorem": "p=3,M>=5: J_ray>2/125",
            },
        },
        "proof_boundary": (
            "The exact Fourier grouping, nonterminal cluster theorem, "
            "finite rational Bernstein absorption, analytic terminal tail, "
            "and length-uniform actual transfer prove J_0>17/1000 and "
            "J_ray>2/125 for every p=3,M>=5 complete chain. Open: "
            "p=2,M>=10, p=2 lengths 2..7, p=3 lengths 2..3, singleton "
            "chains, p-free bases, recurrent endpoints, cutoff equality, "
            "adjacent charts, the Xi Abel gap, joined winding, contact "
            "exclusion, Lambda<=0, PF-infinity, RH, and a prize-level proof."
        ),
        "summary": {
            "rows": len(rows),
            "exact_small_d_cluster_certificates": 5,
            "exact_zero_reserve_certificates": 2,
            "exact_finite_energy_certificates": 12,
            "analytic_terminal_tail_families": 1,
            "ideal_all_length_ternary_families": 1,
            "actual_all_length_ternary_families": 1,
            "open_small_prime_length_families": 1,
            "joined_abel_gaps": 0,
            "successor_winding_bounds": 0,
        },
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    finite_rows = exact["terminal_defect"]["finite_certificates"]
    cluster_rows = exact["nonterminal"]["certificates"]
    return "\n".join(
        [
            "# Ternary All-Length Prime-Power Heat-Starlikeness Gate",
            "",
            "Date: 2026-07-30",
            "",
            "Status: exact all-length complete-chain theorem for p=3,M>=5.",
            "This does not promote RH or Lambda<=0.",
            "",
            "```text",
            "work/rh_compute/results/"
            "jensen_window_pf_newman_polymath15_first_order_centered_"
            "prime_power_heat_starlikeness_ternary_length_gate.json",
            "python work/rh_compute/scripts/"
            "jensen_window_pf_newman_polymath15_first_order_centered_"
            "prime_power_heat_starlikeness_ternary_length_gate.py",
            "python work/rh_compute/scripts/"
            "check_jensen_window_pf_newman_polymath15_first_order_centered_"
            "prime_power_heat_starlikeness_ternary_length_gate.py",
            "```",
            "",
            "## Outcome",
            "",
            "```text",
            exact["ideal_theorem"]["uniform"],
            exact["actual_transfer"]["actual_theorem"],
            "```",
            "",
            "This closes the complete-chain ternary family from length five",
            "onward. It does not close short chains or any global join.",
            "",
            "## Unique Defect",
            "",
            "With `D_d=A_d-2A_(d+1)+A_(d+2)`, every common group for",
            "`d>=1` is positive because `1/sqrt(3)<3/5`. Pairing the last",
            "common group with the boundary pair gives",
            "",
            "```text",
            exact["nonterminal"]["cluster"],
            exact["nonterminal"]["slope_lower"],
            exact["nonterminal"]["large_d_lower"],
            "```",
            "",
            "The exact small-`d` certificates are:",
            "",
            "| d | claimed lower | shifted Bernstein minimum | degree |",
            "|---:|---:|---:|---:|",
            *[
                "| {d} | {claimed_lower} | {minimum_shifted_coefficient_decimal} | "
                "{power_degrees[0]} x {power_degrees[1]} |".format(**row)
                for row in cluster_rows
            ],
            "",
            "Therefore only `D_(M-2)` may be negative.",
            "",
            "## Positive Reserve",
            "",
            "```text",
            exact["zero_reserve"]["theorem"],
            exact["terminal_defect"]["finite_energy"],
            exact["terminal_defect"]["negative_implication"],
            exact["terminal_defect"]["tail_bound"],
            "```",
            "",
            "The exact finite absorption certificates are:",
            "",
            "| M | claimed lower | shifted Bernstein minimum | degree |",
            "|---:|---:|---:|---:|",
            *[
                "| {M} | {claimed_lower} | {minimum_shifted_coefficient_decimal} | "
                "{power_degrees[0]} x {power_degrees[1]} |".format(**row)
                for row in finite_rows
            ],
            "",
            "Since `0<=F_(M-2)<=M-1`, the finite certificates give",
            "`J_0>1/50`; the analytic tail gives `J_0>17/1000`.",
            "",
            "## Actual Fixed-Ray Transfer",
            "",
            "```text",
            exact["actual_transfer"]["coefficient_error"],
            exact["actual_transfer"]["geometric_sums"],
            exact["actual_transfer"]["frozen_error"],
            exact["actual_transfer"]["ray_error"],
            exact["actual_transfer"]["actual_theorem"],
            "```",
            "",
            "No new `pi` is inserted here. Any `pi` in `x=4*pi*exp(L)`",
            "is inherited from the completed-zeta saddle normalization;",
            "the Fejer kernels use the ordinary `2*pi` angular period.",
            "",
            "## Boundary",
            "",
            payload["proof_boundary"],
            "",
        ]
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    arguments = build_parser().parse_args()
    payload = build_payload()
    arguments.out.parent.mkdir(parents=True, exist_ok=True)
    arguments.note.parent.mkdir(parents=True, exist_ok=True)
    arguments.out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    arguments.note.write_text(render_note(payload), encoding="utf-8")
    summary = payload["summary"]
    print(
        "built ternary all-length prime-power heat-starlikeness gate: "
        f"{summary['rows']} rows, "
        f"{summary['exact_small_d_cluster_certificates']} small-d clusters, "
        f"{summary['exact_finite_energy_certificates']} finite energies, "
        "1 analytic terminal tail, 1 ideal all-length ternary family, "
        "1 actual all-length ternary family, "
        "1 open small-prime length family, "
        "0 joined Abel gaps, 0 successor winding bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
