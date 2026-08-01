#!/usr/bin/env python3
"""Build the dyadic all-length prime-power heat-starlikeness gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from dyadic_heat_starlikeness_bernstein import (
    RATIO_UPPER,
    certify_finite_core,
    certify_initial_core,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_dyadic_length_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
LENGTH_STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_length_propagation_gate"
)
DEFECT_STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_dyadic_defect_gate"
)
SOURCES = {
    "length_propagation": (
        REPO_ROOT / "work/rh_compute/results" / f"{LENGTH_STEM}.json"
    ),
    "dyadic_defect": (
        REPO_ROOT / "work/rh_compute/results" / f"{DEFECT_STEM}.json"
    ),
}

X_LEAVES = tuple(
    "".join(directions)
    for directions in product("LR", repeat=4)
)


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


def terminal_square_guards(M: int, target: Fraction) -> dict:
    n = M - 1
    A = 2 * (M + 1) * RATIO_UPPER - M
    denominator = M - (M + 1) * RATIO_UPPER
    low_coefficient = Fraction(M + 1) * A / denominator
    low_base = A / (2 * (M + 2))
    low_margin = target**2 - low_coefficient**2 * low_base**n

    high_coefficient = Fraction(4 * (M + 1), M) * A
    high_radicand = (
        (A / (M + 2)) ** n
        * Fraction(n**n, (n + 2) ** (n + 2))
    )
    high_margin = (
        target**2 - high_coefficient**2 * high_radicand
    )
    return {
        "M": M,
        "target": fraction_text(target),
        "A_upper": fraction_text(A),
        "denominator_lower": fraction_text(denominator),
        "low_square_margin": fraction_text(low_margin),
        "low_square_margin_decimal": decimal_text(low_margin),
        "high_square_margin": fraction_text(high_margin),
        "high_square_margin_decimal": decimal_text(high_margin),
        "positive": low_margin > 0 and high_margin > 0,
    }


def tail_sequence_guards() -> dict:
    M = 16
    alpha_numerator = M - 1
    low_coefficient = Fraction(21, 10) * (M + 1)
    low_base = Fraction(9, 40)
    low_margin = (
        Fraction(1, 250) ** 2
        - low_coefficient**2 * low_base**alpha_numerator
    )
    high_coefficient = (
        Fraction(9, 10)
        * Fraction((M + 1) * (M + 2), 21 * M)
    )
    high_base = Fraction(9, 20)
    high_margin = (
        Fraction(1, 250) ** 2
        - high_coefficient**2 * high_base**alpha_numerator
    )
    fejer_max_margin = (
        17**17 - 4 * 21**2 * 15**15
    )
    return {
        "base_ratio": "A_M/(M+2)<9/20 for M>=15",
        "denominator_ratio": (
            "A_M/[M-(M+1)u]<21/10 for M>=16"
        ),
        "half_power_max": (
            "max t^(15/2)(1-t)<1/21"
        ),
        "half_power_integer_margin": fejer_max_margin,
        "low_M16_square_margin": fraction_text(low_margin),
        "low_M16_square_margin_decimal": decimal_text(low_margin),
        "high_M16_square_margin": fraction_text(high_margin),
        "high_M16_square_margin_decimal": decimal_text(high_margin),
        "low_successor_ratio": "<9/17",
        "high_successor_ratio": "<51/64",
        "positive": (
            fejer_max_margin > 0
            and low_margin > 0
            and high_margin > 0
        ),
    }


def build_certificates() -> dict:
    finite = [
        certify_finite_core(
            M,
            x_leaf_paths=X_LEAVES,
        ).to_dict()
        for M in range(10, 15)
    ]
    initial = certify_initial_core(
        x_leaf_paths=X_LEAVES
    ).to_dict()
    if not all(
        certificate["positive"]
        for certificate in [*finite, initial]
    ):
        raise RuntimeError("a dyadic core certificate is nonpositive")

    finite_terminal = [
        terminal_square_guards(M, Fraction(1, 25))
        for M in range(10, 15)
    ]
    m15_terminal = terminal_square_guards(15, Fraction(1, 250))
    tail = tail_sequence_guards()
    if not all(row["positive"] for row in finite_terminal):
        raise RuntimeError("a finite terminal guard is nonpositive")
    if not m15_terminal["positive"] or not tail["positive"]:
        raise RuntimeError("an all-length terminal guard is nonpositive")
    return {
        "finite_seven_kernel": finite,
        "all_length_initial_four_kernel": initial,
        "finite_terminal": finite_terminal,
        "m15_terminal": m15_terminal,
        "tail_terminal": tail,
    }


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
    defect_summary = payloads["dyadic_defect"].get("summary", {})
    if length_summary.get("exact_dyadic_next_length_certificates") != 1:
        raise RuntimeError("length-propagation source drifted")
    if defect_summary.get("remaining_possible_defects") != 2:
        raise RuntimeError("dyadic-defect source drifted")
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
            "pphsdl_01_domain",
            "dyadic complete-chain domain",
            "proved",
            "Fix p=2 and every complete-chain length M>=10.",
            "47/50<=q<=1; 22/25<=s<=1; a=s/sqrt(2).",
            "The certified M=8 and M=9 families are inherited separately.",
        ),
        GateRow(
            "pphsdl_02_localization",
            "two-defect input",
            "proved",
            "All Fourier second differences except D_0 and D_(M-2) are positive.",
            "D_d>0 for 1<=d<=M-3; D_(M-1)>0.",
            "This imports the independently checked all-length defect gate.",
        ),
        GateRow(
            "pphsdl_03_fejer",
            "positive-kernel reconstruction",
            "proved",
            "Use unnormalized Fejer numerators K_d=|1+z+...+z^d|^2.",
            "2J_0=sum_(d=0)^(M-1)D_dK_d; K_d>=0.",
            "Only the two localized defects need absorption.",
        ),
        GateRow(
            "pphsdl_04_terminal_pair",
            "completed-square terminal pair",
            "proved",
            "The final positive kernel controls a negative penultimate defect by one exact square completion.",
            "P|Y+u|^2-N|Y|^2>=[-PN/(P-N)]; P=D_(M-1), N=-D_(M-2).",
            "The denominator P-N is proved positive, not assumed.",
        ),
        GateRow(
            "pphsdl_05_terminal_parameters",
            "one-variable terminal loss",
            "proved",
            "A negative terminal defect reduces to a scalar t in (0,1).",
            "L_M=(M+1)A^(alpha+1)/(M+2)^alpha * t^alpha(1-t)/[M-(M+1)r+At].",
            "Here alpha=(M-1)/2, A=2(M+1)r-M, and r<=1/sqrt(2).",
        ),
        GateRow(
            "pphsdl_06_finite_core",
            "finite seven-kernel core",
            "proved",
            "For each 10<=M<=14, sixteen exact angular leaves certify the first eight Fejer terms.",
            "D_0+sum_(d=1)^7D_dK_d>1/10.",
            "The certificate covers 0<=r=aq<=707107/10^6.",
        ),
        GateRow(
            "pphsdl_07_finite_terminal",
            "finite terminal absorption",
            "proved",
            "Exact rational square guards bound the terminal loss for 10<=M<=14.",
            "PN/(P-N)<1/25.",
            "Both t<=1/2 and t>=1/2 regimes are checked exactly.",
        ),
        GateRow(
            "pphsdl_08_finite_theorem",
            "finite ideal theorem",
            "proved",
            "The seven-kernel core absorbs the complete terminal pair.",
            "10<=M<=14: 2J_0>1/10-1/25=3/50; J_0>3/100.",
            "All omitted middle kernels are nonnegative by defect localization.",
        ),
        GateRow(
            "pphsdl_09_initial_groups",
            "all-length initial grouping",
            "proved",
            "For M>=15, retain k=0,...,7 in D_0,...,D_4 and discard only positive groups.",
            "c_k=R^kq^[k(25-k)]; R=r_12; 0<=R<=1/sqrt(2).",
            "The D_0 terminal cluster is positive from M=15 onward.",
        ),
        GateRow(
            "pphsdl_10_initial_core",
            "all-length four-kernel core",
            "proved",
            "Sixteen exact angular leaves certify the length-free retained core.",
            "sum_(d=0)^4 K_d sum_(k=0)^7 T_(k,d)>1/100.",
            "The certificate uses y=q^2 and the enlarged full ratio interval.",
        ),
        GateRow(
            "pphsdl_11_m15_terminal",
            "first analytic terminal row",
            "proved",
            "The split scalar loss is below 1/250 at M=15.",
            "max(low-t bound, high-t bound)<1/250.",
            "The comparisons are exact rational square inequalities.",
        ),
        GateRow(
            "pphsdl_12_terminal_tail",
            "all-length terminal tail",
            "proved",
            "The same loss stays below 1/250 for every M>=16.",
            "low successor ratio<9/17; high successor ratio<51/64.",
            "Both majorants decrease geometrically from their exact M=16 checks.",
        ),
        GateRow(
            "pphsdl_13_tail_theorem",
            "all-length ideal theorem",
            "proved",
            "The initial four-kernel core absorbs the completed-square terminal loss.",
            "M>=15: 2J_0>1/100-1/250=3/500; J_0>3/1000.",
            "This is an analytic all-length conclusion after one fixed certificate.",
        ),
        GateRow(
            "pphsdl_14_ideal_theorem",
            "uniform dyadic ideal theorem",
            "proved",
            "Every correction-free dyadic complete chain of length M>=10 is strictly heat-starlike.",
            "p=2,M>=10: J_0>3/1000.",
            "Together with earlier gates, complete dyadic chains are covered from M=8.",
        ),
        GateRow(
            "pphsdl_15_actual_sums",
            "uniform coefficient transfer",
            "proved",
            "The dyadic geometric coefficient sums are independent of M.",
            "sum c_k<4; sum k c_k<12; absolute current mass<64.",
            "These are the same length-uniform bounds used at M=9.",
        ),
        GateRow(
            "pphsdl_16_actual_current",
            "actual fixed-ray theorem",
            "proved",
            "The rho and epsilon corrections cost less than 1/1000 uniformly in M.",
            "|J_ray-J_0|<3456096/x+1621632/x^2<1/1000; J_ray>1/500.",
            "This remains local to one complete dyadic prime-power ray.",
        ),
        GateRow(
            "pphsdl_17_pi_provenance",
            "constant provenance",
            "proved",
            "No new circle constant enters the dyadic kernel argument.",
            "z=exp(i theta), x_ang=cos(theta), and K_d have the ordinary 2*pi angular period.",
            "Any pi in the saddle scale is inherited from completed-zeta normalization.",
        ),
        GateRow(
            "pphsdl_18_boundary",
            "proof boundary",
            "open",
            "The all-length dyadic complete-chain obstruction is closed, but global joins remain.",
            "Open: short p=2 lengths 2..7, short p=3 lengths 2..3, singleton and p-free chains, endpoint, adjacent, Abel, and winding joins.",
            "No contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level proof is claimed.",
        ),
    ]


def build_payload() -> dict:
    source_payloads = load_sources()
    certificates = build_certificates()
    rows = build_rows()
    return {
        "schema_version": 1,
        "kind": "prime_power_heat_starlikeness_dyadic_length_gate",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "title": "Dyadic all-length prime-power heat-starlikeness gate",
        "source_audit": audit_sources(source_payloads),
        "assumptions": {
            "prime": 2,
            "block_lengths": "M>=10",
            "q_interval": ["47/50", "1"],
            "s_interval": ["22/25", "1"],
            "a_definition": "a=s/sqrt(2)",
            "ratio_upper": "707107/1000000",
            "ratio_upper_check": (
                "2*(707107)^2-10^12=618898>0"
            ),
            "angular_variable": "x_ang=cos(theta)",
        },
        "rows": [asdict(row) for row in rows],
        "exact": {
            "fejer": {
                "numerator": "K_d=|sum_(j=0)^d z^j|^2",
                "reconstruction": (
                    "2J_0=sum_(d=0)^(M-1)D_dK_d"
                ),
                "localization": (
                    "D_d>0 for 1<=d<=M-3 and D_(M-1)>0"
                ),
            },
            "terminal_pair": {
                "completion": (
                    "P|Y+u|^2-N|Y|^2=(P-N)"
                    "|Y+P*u/(P-N)|^2-PN/(P-N)"
                ),
                "positive_denominator": (
                    "P-N=c_(M-2)[M-(M+1)r+(M+2)h]>0"
                ),
                "parameters": (
                    "r=aq; v=q^(2M-4); h=r^2v; "
                    "A=2(M+1)r-M; alpha=(M-1)/2; "
                    "t=v(M+2)r^2/A"
                ),
                "loss": (
                    "L_M=(M+1)A^(alpha+1)/(M+2)^alpha"
                    "*t^alpha(1-t)/[M-(M+1)r+At]"
                ),
            },
            "finite": {
                "core_certificates": certificates[
                    "finite_seven_kernel"
                ],
                "terminal_guards": certificates["finite_terminal"],
                "core_theorem": (
                    "10<=M<=14: D_0+sum_(d=1)^7D_dK_d>1/10"
                ),
                "terminal_theorem": (
                    "10<=M<=14: terminal pair>-1/25"
                ),
                "ideal_theorem": "10<=M<=14: J_0>3/100",
            },
            "all_length": {
                "retained_groups": (
                    "T_(k,d), 0<=k<=7, 0<=d<=4"
                ),
                "d0_terminal_guard": (
                    "24delta+kappa=77/4-(40/3)sqrt(2)>0; "
                    "231^2-2*160^2=2161"
                ),
                "coefficient_anchor": (
                    "R=r_12; c_k=R^kq^[k(25-k)]"
                ),
                "initial_core_certificate": certificates[
                    "all_length_initial_four_kernel"
                ],
                "initial_core_theorem": (
                    "sum_(d=0)^4K_d sum_(k=0)^7T_(k,d)>1/100"
                ),
                "m15_terminal_guard": certificates["m15_terminal"],
                "tail_terminal_guard": certificates["tail_terminal"],
                "terminal_theorem": "M>=15: terminal pair>-1/250",
                "ideal_theorem": "M>=15: J_0>3/1000",
            },
            "ideal_theorem": {
                "uniform": "p=2,M>=10: J_0>3/1000",
                "complete_dyadic_coverage": (
                    "with preceding M=8 and M=9 gates: p=2,M>=8"
                ),
            },
            "actual_transfer": {
                "coefficient_error": "|R_k-1|<r_x=18000/x",
                "geometric_sums": (
                    "sum c_k<4; sum k c_k<12; "
                    "absolute current mass<64"
                ),
                "frozen_error": "|J_ang-J_0|<192r_x=3456000/x",
                "ray_error": (
                    "|J_ray-J_0|<3456096/x+1621632/x^2<1/1000"
                ),
                "actual_theorem": "p=2,M>=10: J_ray>1/500",
            },
        },
        "proof_boundary": (
            "The exact seven-kernel finite certificates, all-length "
            "four-kernel initial-group certificate, completed-square "
            "terminal estimate, and length-uniform actual transfer prove "
            "J_0>3/1000 and J_ray>1/500 for every p=2,M>=10 "
            "complete chain. Together with preceding gates, complete "
            "dyadic chains are covered for M>=8. Open: p=2 lengths 2..7, "
            "p=3 lengths 2..3, singleton chains, p-free bases, recurrent "
            "endpoints, cutoff equality, adjacent charts, the Xi Abel "
            "gap, joined winding, contact exclusion, Lambda<=0, "
            "PF-infinity, RH, and a prize-level proof."
        ),
        "summary": {
            "rows": len(rows),
            "exact_finite_core_certificates": 5,
            "exact_initial_core_certificates": 1,
            "exact_finite_terminal_guards": 5,
            "analytic_terminal_tail_families": 1,
            "ideal_all_length_dyadic_families": 1,
            "actual_all_length_dyadic_families": 1,
            "complete_dyadic_minimum_length": 8,
            "remaining_possible_defects": 0,
            "joined_abel_gaps": 0,
            "successor_winding_bounds": 0,
        },
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    finite = exact["finite"]["core_certificates"]
    initial = exact["all_length"]["initial_core_certificate"]
    return "\n".join(
        [
            "# Dyadic All-Length Prime-Power Heat-Starlikeness Gate",
            "",
            "Date: 2026-07-30",
            "",
            "Status: exact all-length complete-chain theorem for p=2,M>=10.",
            "This does not promote RH or Lambda<=0.",
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
            exact["ideal_theorem"]["uniform"],
            exact["actual_transfer"]["actual_theorem"],
            exact["ideal_theorem"]["complete_dyadic_coverage"],
            "```",
            "",
            "## Two-Defect Absorption",
            "",
            "With `K_d=|1+z+...+z^d|^2`,",
            "",
            "```text",
            exact["fejer"]["reconstruction"],
            exact["terminal_pair"]["completion"],
            exact["terminal_pair"]["loss"],
            "```",
            "",
            "The last two kernels are kept together. Completing the square",
            "replaces the old worst-case quadratic Fejer loss by one scalar",
            "`t^alpha(1-t)` loss.",
            "",
            "## Finite Rows",
            "",
            "| M | target | shifted Bernstein minimum | degree |",
            "|---:|---:|---:|---:|",
            *[
                "| {block_length} | {claimed_rational_lower} | "
                "{global_minimum_lower} | "
                "{power_degrees[0]} x {power_degrees[1]} x "
                "{power_degrees[2]} |".format(**certificate)
                for certificate in finite
            ],
            "",
            "```text",
            exact["finite"]["core_theorem"],
            exact["finite"]["terminal_theorem"],
            exact["finite"]["ideal_theorem"],
            "```",
            "",
            "## All-Length Tail",
            "",
            "For `M>=15`, retain seven initial groups in each of",
            "`D_0,...,D_4`. With `R=r_12`,",
            "",
            "```text",
            exact["all_length"]["coefficient_anchor"],
            exact["all_length"]["initial_core_theorem"],
            exact["all_length"]["terminal_theorem"],
            exact["all_length"]["ideal_theorem"],
            "```",
            "",
            "The length-free core certificate has:",
            "",
            "| target | shifted Bernstein minimum | degree |",
            "|---:|---:|---:|",
            "| {claimed_rational_lower} | {global_minimum_lower} | "
            "{power_degrees[0]} x {power_degrees[1]} x "
            "{power_degrees[2]} |".format(**initial),
            "",
            "## Actual Fixed-Ray Transfer",
            "",
            "```text",
            exact["actual_transfer"]["coefficient_error"],
            exact["actual_transfer"]["geometric_sums"],
            exact["actual_transfer"]["ray_error"],
            exact["actual_transfer"]["actual_theorem"],
            "```",
            "",
            "No new `pi` is inserted here. The Fejer numerators use the",
            "ordinary `2*pi` period of `z=exp(i theta)`; any `pi` in the",
            "saddle scale is inherited from completed-zeta normalization.",
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
        "built dyadic all-length prime-power heat-starlikeness gate: "
        f"{summary['rows']} rows, "
        f"{summary['exact_finite_core_certificates']} finite cores, "
        f"{summary['exact_initial_core_certificates']} initial core, "
        f"{summary['exact_finite_terminal_guards']} finite terminal guards, "
        "1 analytic terminal tail, 1 ideal all-length dyadic family, "
        "1 actual all-length dyadic family, "
        "complete dyadic minimum length 8, 0 remaining defects, "
        "0 joined Abel gaps, 0 successor winding bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
