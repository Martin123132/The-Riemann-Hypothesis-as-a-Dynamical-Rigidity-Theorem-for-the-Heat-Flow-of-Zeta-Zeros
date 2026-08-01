#!/usr/bin/env python3
"""Build the dyadic all-length Fourier-defect localization gate."""

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
    "prime_power_heat_starlikeness_dyadic_defect_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
LENGTH_STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_length_propagation_gate"
)
TERNARY_STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_ternary_length_gate"
)
SOURCES = {
    "length_propagation": (
        REPO_ROOT / "work/rh_compute/results" / f"{LENGTH_STEM}.json"
    ),
    "ternary_length": (
        REPO_ROOT / "work/rh_compute/results" / f"{TERNARY_STEM}.json"
    ),
}

Q_LOWER = Fraction(47, 50)
A_UPPER = Fraction(707107, 1_000_000)


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


def decimal_text(value: Fraction, digits: int = 40) -> str:
    with localcontext() as context:
        context.prec = digits
        return format(
            Decimal(value.numerator) / Decimal(value.denominator), ".20e"
        )


def add_term(
    terms: dict[tuple[int, int], Fraction],
    q_degree: int,
    second_degree: int,
    coefficient: int | Fraction,
) -> None:
    if q_degree < 0 or second_degree < 0:
        raise ValueError("negative polynomial degree")
    key = (q_degree, second_degree)
    terms[key] = terms.get(key, Fraction()) + Fraction(coefficient)
    if not terms[key]:
        del terms[key]


def interval_monomial_bernstein(
    monomial_degree: int,
    elevated_degree: int,
    lower: Fraction,
    upper: Fraction,
) -> list[Fraction]:
    width = upper - lower
    values = []
    for index in range(elevated_degree + 1):
        value = Fraction()
        for power in range(min(monomial_degree, index) + 1):
            value += (
                Fraction(comb(monomial_degree, power))
                * lower ** (monomial_degree - power)
                * width**power
                * Fraction(
                    comb(index, power), comb(elevated_degree, power)
                )
            )
        values.append(value)
    return values


def bernstein_certificate(
    source: dict[tuple[int, int], Fraction],
    second_interval: tuple[Fraction, Fraction],
    claimed_lower: Fraction,
) -> dict:
    terms = dict(source)
    add_term(terms, 0, 0, -claimed_lower)
    q_degree = max(q_degree for q_degree, _ in terms)
    second_degree = max(second_degree for _, second_degree in terms)
    q_vectors = {
        degree: interval_monomial_bernstein(
            degree, q_degree, Q_LOWER, Fraction(1)
        )
        for degree in {q_degree for q_degree, _ in terms}
    }
    second_vectors = {
        degree: interval_monomial_bernstein(
            degree,
            second_degree,
            second_interval[0],
            second_interval[1],
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
        "q_interval": ["47/50", "1"],
        "second_interval": [
            fraction_text(second_interval[0]),
            fraction_text(second_interval[1]),
        ],
        "power_degrees": [q_degree, second_degree],
        "bernstein_shape": [q_degree + 1, second_degree + 1],
        "claimed_lower": fraction_text(claimed_lower),
        "minimum_shifted_coefficient": fraction_text(minimum),
        "minimum_shifted_coefficient_decimal": decimal_text(minimum),
        "argminimum": list(argminimum),
        "positive": minimum > 0,
    }


def terminal_cluster_terms(
    B: int, d: int
) -> dict[tuple[int, int], Fraction]:
    terms: dict[tuple[int, int], Fraction] = {}
    for q_degree, a_degree, coefficient in (
        (0, 0, B),
        (3, 1, -2 * (B + 1)),
        (4, 2, B + 2),
        (2 * d + 6, 2, B + 2),
        (2 * d + 7, 3, -2 * (B + 3)),
        (4 * d + 8, 4, B + 4),
    ):
        add_term(terms, q_degree, a_degree, coefficient)
    return terms


def initial_cluster_terms(
    d: int, count: int
) -> tuple[dict[tuple[int, int], Fraction], int]:
    raw: list[tuple[int, int, int]] = []
    minimum_q_degree = 0
    for k in range(count):
        A = 2 * k + d + 2
        entries = (
            (2 * k * (d + 1 - k), 2 * k, A),
            (2 * k * (d - k), 2 * k + 1, -2 * (A + 1)),
            (
                2 * k * (d - 1 - k) - 2,
                2 * k + 2,
                A + 2,
            ),
        )
        raw.extend(entries)
        minimum_q_degree = min(
            minimum_q_degree, *(q_degree for q_degree, _, _ in entries)
        )
    shift = -minimum_q_degree
    terms: dict[tuple[int, int], Fraction] = {}
    for q_degree, ratio_degree, coefficient in raw:
        add_term(
            terms, q_degree + shift, ratio_degree, coefficient
        )
    return terms, shift


def normalized_difference_terms(
    M: int, d: int
) -> dict[tuple[int, int], Fraction]:
    exponents = [k * (2 * M - 2 - k) for k in range(M)]
    terms: dict[tuple[int, int], Fraction] = {}
    for increment, scale in ((0, 1), (1, -2), (2, 1)):
        offset = d + increment
        for k in range(M - offset):
            add_term(
                terms,
                exponents[k] + exponents[k + offset] - exponents[d],
                2 * k + increment,
                scale * (2 * k + offset + 2),
            )
    return terms


def build_certificates() -> dict:
    d1_terms, d1_shift = initial_cluster_terms(1, 7)
    d2_terms, d2_shift = initial_cluster_terms(2, 3)
    certificates = {
        "d1_initial_seven": {
            **bernstein_certificate(
                d1_terms,
                (Fraction(3, 5), A_UPPER),
                Fraction(1, 200),
            ),
            "d": 1,
            "common_group_count": 7,
            "q_shift": d1_shift,
        },
        "d2_initial_three": {
            **bernstein_certificate(
                d2_terms,
                (Fraction(2, 3), A_UPPER),
                Fraction(1, 25),
            ),
            "d": 2,
            "common_group_count": 3,
            "q_shift": d2_shift,
        },
        "d1_m10_full": bernstein_certificate(
            normalized_difference_terms(10, 1),
            (Fraction(0), A_UPPER),
            Fraction(1, 100),
        ),
        "d1_terminal_B15": bernstein_certificate(
            terminal_cluster_terms(15, 1),
            (Fraction(0), A_UPPER),
            Fraction(1),
        ),
        "d2_terminal_B14": bernstein_certificate(
            terminal_cluster_terms(14, 2),
            (Fraction(0), A_UPPER),
            Fraction(1),
        ),
        "d3_to_d14_terminal": [],
    }
    for d in range(3, 15):
        certificate = bernstein_certificate(
            terminal_cluster_terms(d + 2, d),
            (Fraction(0), A_UPPER),
            Fraction(1, 25),
        )
        certificate["d"] = d
        certificates["d3_to_d14_terminal"].append(certificate)
    flattened = [
        certificates["d1_initial_seven"],
        certificates["d2_initial_three"],
        certificates["d1_m10_full"],
        certificates["d1_terminal_B15"],
        certificates["d2_terminal_B14"],
        *certificates["d3_to_d14_terminal"],
    ]
    if not all(certificate["positive"] for certificate in flattened):
        failed = [
            certificate
            for certificate in flattened
            if not certificate["positive"]
        ]
        raise RuntimeError(f"nonpositive dyadic certificate: {failed}")
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
    length = payloads["length_propagation"].get("summary", {})
    ternary = payloads["ternary_length"].get("summary", {})
    if length.get("exact_dyadic_next_length_certificates") != 1:
        raise RuntimeError("length-propagation source drifted")
    if ternary.get("ideal_all_length_ternary_families") != 1:
        raise RuntimeError("ternary source drifted")
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
            "pphsdd_01_domain",
            "dyadic complete-chain domain",
            "proved",
            "Fix p=2 and every complete-chain length M>=10.",
            "47/50<=q<=1; 22/25<=s<=1; a=s/sqrt(2)<=1/sqrt(2).",
            "The already certified M=9 family and shorter chains remain separate.",
        ),
        GateRow(
            "pphsdd_02_fourier_groups",
            "exact second-difference grouping",
            "proved",
            "Expand every D_d into common interior groups and one terminal cluster.",
            "D_d=A_d-2A_(d+1)+A_(d+2); C(B,d)=I_B+a^2q^(2d+6)G_(B+2,d).",
            "No coefficient sign is inferred from a numerical scan.",
        ),
        GateRow(
            "pphsdd_03_common_threshold",
            "ordinary common groups",
            "proved",
            "Every common group with coefficient A>=5 is strictly positive.",
            "T_A>=(1-r)[A-(A+2)r]>0; r<=1/sqrt(2)<5/7.",
            "Only A=3 and A=4 require low-offset clusters.",
        ),
        GateRow(
            "pphsdd_04_cluster_slope",
            "terminal-cluster monotonicity",
            "proved",
            "The terminal cluster is strictly increasing in B throughout the dyadic box.",
            "partial_B C=E_d>(3-2sqrt(2))/4>1/25.",
            "The square-completion retains the q-dependent terminal term.",
        ),
        GateRow(
            "pphsdd_05_d1_initial",
            "first low-offset cluster",
            "proved",
            "For M>=11 the first seven d=1 common groups absorb the sole A=3 group.",
            "r_1<=3/5 gives termwise positivity; r_1>=3/5 gives q^74 S_1>1/200.",
            "M=10 is handled by a separate full normalized certificate.",
        ),
        GateRow(
            "pphsdd_06_d1_m10",
            "shortest dyadic exception",
            "proved",
            "The exact normalized M=10 first-offset difference is positive.",
            "D_1/c_1>1/100 for M=10.",
            "This is a bivariate exact-rational certificate, not floating point evidence.",
        ),
        GateRow(
            "pphsdd_07_d2_initial",
            "second low-offset cluster",
            "proved",
            "The first three d=2 common groups absorb the sole A=4 group.",
            "r_2<=2/3 gives termwise positivity; r_2>=2/3 gives q^6 S_2>1/25.",
            "The cluster exists for every M>=10.",
        ),
        GateRow(
            "pphsdd_08_low_terminal",
            "low-offset terminal anchors",
            "proved",
            "The actual d=1,2 terminal clusters are positive at their shortest lengths.",
            "C(15,1)>1; C(14,2)>1; larger B follows from the slope theorem.",
            "The negative formal anchors C(3,1),C(4,2) are never used.",
        ),
        GateRow(
            "pphsdd_09_middle_terminal",
            "finite terminal offsets",
            "proved",
            "Exact Bernstein coefficients prove every minimal terminal cluster for 3<=d<=14.",
            "C(d+2,d)>1/25 for d=3,...,14.",
            "The certificates use the enlarged box 0<=a<=707107/10^6.",
        ),
        GateRow(
            "pphsdd_10_terminal_tail",
            "analytic terminal tail",
            "proved",
            "Every minimal terminal cluster with d>=15 is positive.",
            "C(d+2,d)=dE_d+K_d; E_d>delta; K_d>=kappa; 15delta+kappa=(75-53sqrt(2))/6>0.",
            "The final radical comparison is exact because 75^2-2*53^2=7.",
        ),
        GateRow(
            "pphsdd_11_localization",
            "all-length defect localization",
            "proved",
            "Every nonterminal positive offset has positive Fourier second difference.",
            "p=2,M>=10: D_d>0 for 1<=d<=M-3, and D_(M-1)=A_(M-1)>0.",
            "No sign is assigned to D_0 or D_(M-2).",
        ),
        GateRow(
            "pphsdd_12_two_defects",
            "sharp remaining obstruction",
            "open",
            "The dyadic Fejer problem is reduced to exactly two possible defects.",
            "Possible negative indices: d=0 and d=M-2.",
            "D_0 is genuinely negative at q=s=1, so ternary reserve absorption cannot be copied.",
        ),
        GateRow(
            "pphsdd_13_boundary",
            "proof boundary",
            "open",
            "A positive-kernel or equivalent phase-aware absorption of the two defects is still required.",
            "Open: prove J_0>0 uniformly for p=2,M>=10, then transfer to J_ray.",
            "No global join, Lambda<=0, PF-infinity, RH, or prize-level conclusion is claimed.",
        ),
    ]


def build_payload() -> dict:
    source_payloads = load_sources()
    certificates = build_certificates()
    rows = build_rows()
    return {
        "schema_version": 1,
        "kind": "prime_power_heat_starlikeness_dyadic_defect_gate",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "title": "Dyadic all-length Fourier-defect localization gate",
        "source_audit": audit_sources(source_payloads),
        "assumptions": {
            "prime": 2,
            "block_lengths": "M>=10",
            "q_interval": ["47/50", "1"],
            "s_interval": ["22/25", "1"],
            "a_definition": "a=s/sqrt(2)",
            "a_upper_rational": "707107/1000000",
            "a_upper_check": (
                "2*(707107)^2-10^12=618898>0"
            ),
            "common_threshold_check": (
                "7*707107=4949749<5000000"
            ),
        },
        "rows": [asdict(row) for row in rows],
        "exact": {
            "fourier": {
                "coefficients": (
                    "A_d=sum_(k=0)^(M-1-d)(2k+d+2)c_kc_(k+d)"
                ),
                "difference": "D_d=A_d-2A_(d+1)+A_(d+2)",
                "common_group": (
                    "c_kc_(k+d)[A-2(A+1)r+(A+2)rr_+]"
                ),
                "terminal_cluster": (
                    "C(B,d)=B-2(B+1)aq^3+(B+2)a^2q^4"
                    "+a^2q^(2d+6)[B+2-2(B+3)aq"
                    "+(B+4)a^2q^(2d+2)]"
                ),
                "actual_B": "B=2M-d-4",
            },
            "slope": {
                "formula": (
                    "E_d=1-2aq^3+a^2q^4+a^2q^(2d+6)"
                    "[1-2aq+a^2q^(2d+2)]"
                ),
                "completion": (
                    "E_d>=(1-a)^2-(2a-1)_+^2/4"
                ),
                "lower": (
                    "E_d>(3-2sqrt(2))/4>1/25"
                ),
                "rational_check": (
                    "(71/25)^2-8=41/625>0"
                ),
            },
            "low_offsets": {
                "normalized_weight": (
                    "c_kc_(k+d)/c_d=R^(2k)q^[2k(d+1-k)], "
                    "R=r_d"
                ),
                "d1_initial_certificate": certificates[
                    "d1_initial_seven"
                ],
                "d1_m10_certificate": certificates["d1_m10_full"],
                "d2_initial_certificate": certificates[
                    "d2_initial_three"
                ],
                "d1_terminal_certificate": certificates[
                    "d1_terminal_B15"
                ],
                "d2_terminal_certificate": certificates[
                    "d2_terminal_B14"
                ],
            },
            "middle_offsets": {
                "certificates": certificates["d3_to_d14_terminal"],
                "theorem": "C(d+2,d)>1/25 for 3<=d<=14",
            },
            "tail": {
                "decomposition": "C(d+2,d)=dE_d+K_d",
                "E_lower": "delta=(3-2sqrt(2))/4",
                "K_lower": "kappa=5/4-(4/3)sqrt(2)",
                "endpoint": (
                    "15delta+kappa=(75-53sqrt(2))/6>0"
                ),
                "endpoint_check": "75^2-2*53^2=7>0",
                "theorem": "C(d+2,d)>0 for d>=15",
            },
            "localization": {
                "theorem": (
                    "p=2,M>=10: D_d>0 for 1<=d<=M-3; "
                    "D_(M-1)>0"
                ),
                "possible_negative_indices": ["0", "M-2"],
                "genuine_zero_defect": (
                    "D_0<0 occurs at q=s=1"
                ),
            },
        },
        "proof_boundary": (
            "The exact common-group decomposition, dyadic cluster-slope "
            "bound, low-offset rational certificates, and analytic "
            "terminal tail prove D_d>0 for every 1<=d<=M-3 and "
            "D_(M-1)>0 when p=2,M>=10. Exactly D_0 and D_(M-2) "
            "remain as possible Fourier defects. This gate does not yet "
            "prove uniform dyadic heat-starlikeness, an actual-ray "
            "theorem, a global join, Lambda<=0, PF-infinity, RH, or a "
            "prize-level result."
        ),
        "summary": {
            "rows": len(rows),
            "exact_initial_cluster_certificates": 2,
            "exact_shortest_length_certificates": 1,
            "exact_low_terminal_certificates": 2,
            "exact_middle_terminal_certificates": 12,
            "analytic_terminal_tail_families": 1,
            "proved_positive_offset_families": 1,
            "remaining_possible_defects": 2,
            "ideal_all_length_dyadic_families": 0,
            "actual_all_length_dyadic_families": 0,
            "joined_abel_gaps": 0,
            "successor_winding_bounds": 0,
        },
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    middle = exact["middle_offsets"]["certificates"]
    low = exact["low_offsets"]
    return "\n".join(
        [
            "# Dyadic All-Length Fourier-Defect Localization Gate",
            "",
            "Date: 2026-07-30",
            "",
            "Status: exact defect-localization theorem for p=2,M>=10.",
            "This does not yet prove dyadic heat-starlikeness or RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Outcome",
            "",
            "```text",
            exact["localization"]["theorem"],
            "possible negative indices: D_0 and D_(M-2)",
            "```",
            "",
            "The old numerical pattern is now an all-length theorem. The",
            "zero-mode defect is genuine, so the remaining step needs a",
            "two-defect positive-kernel or equivalent phase-aware argument.",
            "",
            "## Terminal Slope",
            "",
            "For the last common group and its two boundary terms,",
            "",
            "```text",
            exact["fourier"]["terminal_cluster"],
            exact["slope"]["formula"],
            exact["slope"]["completion"],
            exact["slope"]["lower"],
            "```",
            "",
            "Thus `C(B,d)` increases with the actual length parameter `B`.",
            "",
            "## Low Offsets",
            "",
            "| certificate | lower | shifted Bernstein minimum | degree |",
            "|---|---:|---:|---:|",
            *[
                "| {label} | {claimed_lower} | "
                "{minimum_shifted_coefficient_decimal} | "
                "{power_degrees[0]} x {power_degrees[1]} |".format(
                    label=label, **certificate
                )
                for label, certificate in (
                    ("d=1 first seven", low["d1_initial_certificate"]),
                    ("d=1, M=10 full", low["d1_m10_certificate"]),
                    ("d=2 first three", low["d2_initial_certificate"]),
                    ("C(15,1)", low["d1_terminal_certificate"]),
                    ("C(14,2)", low["d2_terminal_certificate"]),
                )
            ],
            "",
            "The `d=1` first-seven certificate is used for `M>=11`; the",
            "separate normalized certificate handles `M=10`. The `d=2`",
            "three-group cluster is available for every `M>=10`.",
            "",
            "## Middle And Tail",
            "",
            "| d | lower | shifted Bernstein minimum | degree |",
            "|---:|---:|---:|---:|",
            *[
                "| {d} | {claimed_lower} | "
                "{minimum_shifted_coefficient_decimal} | "
                "{power_degrees[0]} x {power_degrees[1]} |".format(
                    **certificate
                )
                for certificate in middle
            ],
            "",
            "For `d>=15`, no finite scan is used:",
            "",
            "```text",
            exact["tail"]["decomposition"],
            exact["tail"]["E_lower"],
            exact["tail"]["K_lower"],
            exact["tail"]["endpoint"],
            exact["tail"]["endpoint_check"],
            "```",
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
        "built dyadic all-length Fourier-defect localization gate: "
        f"{summary['rows']} rows, "
        f"{summary['exact_initial_cluster_certificates']} initial clusters, "
        f"{summary['exact_shortest_length_certificates']} shortest-length "
        "certificate, "
        f"{summary['exact_low_terminal_certificates']} low terminal anchors, "
        f"{summary['exact_middle_terminal_certificates']} middle terminal "
        "clusters, 1 analytic terminal tail, "
        "1 proved positive-offset family, 2 possible defects, "
        "0 ideal all-length dyadic families, "
        "0 actual all-length dyadic families, "
        "0 joined Abel gaps, 0 successor winding bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
