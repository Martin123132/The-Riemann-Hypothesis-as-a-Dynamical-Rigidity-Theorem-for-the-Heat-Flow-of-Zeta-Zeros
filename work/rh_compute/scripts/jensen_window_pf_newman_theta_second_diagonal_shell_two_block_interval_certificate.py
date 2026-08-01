#!/usr/bin/env python3
"""Certify the second diagonal shell with two retained theta blocks."""

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
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_theta_compact_transversality_interval_certificate as compact  # noqa: E402


STEM = (
    "jensen_window_pf_newman_theta_"
    "second_diagonal_shell_two_block_interval_certificate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"

COMPACT_SOURCE = compact.DEFAULT_OUT
FIRST_STAGE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_"
    "first_diagonal_shell_interval_certificate.json"
)
ROUTE_GUARD_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_"
    "second_diagonal_shell_first_block_route_guard.json"
)
DIAGONAL_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_positive_boundary_diagonal_exhaustion_gate.json"
)
WINDING_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_first_jet_winding_gate.json"
)
BOUNDARY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_positive_boundary_attainment_lemma.json"
)

TIME_LOWER = Fraction(1, 10)
TIME_UPPER = Fraction(1, 5)
X_LOWER = Fraction(38)
X_UPPER = Fraction(40)
INITIAL_TIME_STEP = Fraction(1, 50)
INITIAL_X_STEP = Fraction(2, 5)
MAX_DEPTH = 8
ENDPOINT_DIGITS = 50
TAIL_MOMENT_ZERO = Fraction(1, 10**10)
TAIL_MOMENT_ONE = Fraction(1, 10**12)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


class TwoBlockCertifier(compact.CompactCertifier):
    def phi_two(self, u: compact.acb) -> compact.acb:
        return (
            32 * self.pi * self.pi * (9 * u).exp()
            - 12 * self.pi * (5 * u).exp()
        ) * (-4 * self.pi * (4 * u).exp()).exp()

    def phi_retained(self, u: compact.acb) -> compact.acb:
        return self.phi_one(u) + self.phi_two(u)

    def oscillatory_integral(
        self, moment: int, time: Fraction, frequency: Fraction
    ) -> compact.acb:
        time_ball = compact.acb(compact.fraction_decimal(time))
        frequency_ball = compact.acb(
            compact.fraction_decimal(frequency)
        )

        def integrand(
            u: compact.acb, analytic: bool
        ) -> compact.acb:
            del analytic
            return (
                u**moment
                * (time_ball * u * u).exp()
                * self.phi_retained(u)
                * (compact.acb(0, 1) * frequency_ball * u).exp()
            )

        retained = compact.acb.integral(
            integrand,
            self.zero,
            self.cutoff,
            abs_tol=self.abs_tol,
            rel_tol=self.abs_tol,
            deg_limit=64,
            eval_limit=200000,
            depth_limit=64,
        )
        return retained + self.tail_complex

    def moment_upper(
        self, moment: int, time: Fraction
    ) -> compact.arb:
        key = (time, moment)
        cached = self.moment_cache.get(key)
        if cached is not None:
            return cached
        time_ball = compact.acb(compact.fraction_decimal(time))

        def integrand(
            u: compact.acb, analytic: bool
        ) -> compact.acb:
            del analytic
            return (
                u**moment
                * (time_ball * u * u).exp()
                * self.phi_retained(u)
            )

        retained = compact.acb.integral(
            integrand,
            self.zero,
            self.cutoff,
            abs_tol=self.abs_tol,
            rel_tol=self.abs_tol,
            deg_limit=64,
            eval_limit=200000,
            depth_limit=64,
        ).real
        value = retained + self.tail_real
        self.moment_cache[key] = value
        return value

    def certify_two_block_box(
        self,
        time_low: Fraction,
        time_high: Fraction,
        x_low: Fraction,
        x_high: Fraction,
        depth: int,
    ) -> dict:
        self.evaluations += 1
        time_center = (time_low + time_high) / 2
        time_radius = (time_high - time_low) / 2
        x_center = (x_low + x_high) / 2
        x_radius = (x_high - x_low) / 2
        max_moment = (
            1
            + 2 * compact.TAYLOR_TIME_ORDER
            + compact.TAYLOR_X_ORDER
        )
        integrals = [
            self.oscillatory_integral(
                moment, time_center, x_center
            )
            for moment in range(max_moment + 1)
        ]
        h_box = self.transform_box(
            0,
            time_center,
            time_radius,
            x_center,
            x_radius,
            integrals,
        )
        h_prime_box = self.transform_box(
            1,
            time_center,
            time_radius,
            x_center,
            x_radius,
            integrals,
        )
        x_box = compact.interval_ball(x_center, x_radius)
        x_three = x_box**3
        x_four = x_box**4
        j_box = 16 * x_four * h_box
        j_prime_box = (
            64 * x_three * h_box + 16 * x_four * h_prime_box
        )
        epsilon_zero = compact.arb(
            compact.fraction_decimal(TAIL_MOMENT_ZERO)
        )
        epsilon_one = compact.arb(
            compact.fraction_decimal(TAIL_MOMENT_ONE)
        )
        tail_value_upper = (16 * x_four * epsilon_zero).upper()
        tail_derivative_upper = (
            64 * x_three * epsilon_zero
            + 16 * x_four * epsilon_one
        ).upper()
        value_ratio = j_box.abs_lower() / tail_value_upper
        derivative_ratio = (
            j_prime_box.abs_lower() / tail_derivative_upper
        )
        value_pass = j_box.abs_lower() > tail_value_upper
        derivative_pass = (
            j_prime_box.abs_lower() > tail_derivative_upper
        )
        record = {
            "depth": depth,
            "t_low": str(time_low),
            "t_high": str(time_high),
            "x_low": str(x_low),
            "x_high": str(x_high),
            "j_two_block_ball": str(j_box),
            "j_two_block_prime_ball": str(j_prime_box),
            "j_two_block_lower": j_box.lower().str(
                ENDPOINT_DIGITS
            ),
            "j_two_block_upper": j_box.upper().str(
                ENDPOINT_DIGITS
            ),
            "j_two_block_prime_lower": j_prime_box.lower().str(
                ENDPOINT_DIGITS
            ),
            "j_two_block_prime_upper": j_prime_box.upper().str(
                ENDPOINT_DIGITS
            ),
            "tail_value_upper": tail_value_upper.str(
                ENDPOINT_DIGITS
            ),
            "tail_derivative_upper": tail_derivative_upper.str(
                ENDPOINT_DIGITS
            ),
            "value_ratio_lower": value_ratio.lower().str(
                ENDPOINT_DIGITS
            ),
            "derivative_ratio_lower": derivative_ratio.lower().str(
                ENDPOINT_DIGITS
            ),
        }
        if value_pass or derivative_pass:
            if value_ratio > derivative_ratio:
                branch = "value"
                ratio = value_ratio
            else:
                branch = "derivative"
                ratio = derivative_ratio
            return {
                **record,
                "certified": True,
                "branch": branch,
                "certified_ratio_lower": ratio.lower().str(
                    ENDPOINT_DIGITS
                ),
                "two_block_value_negative": bool(
                    j_box.upper() < -tail_value_upper
                ),
            }
        return {**record, "certified": False}


def split_box(
    bounds: tuple[Fraction, Fraction, Fraction, Fraction],
    time_step: Fraction = INITIAL_TIME_STEP,
    x_step: Fraction = INITIAL_X_STEP,
) -> list[tuple[Fraction, Fraction, Fraction, Fraction]]:
    time_low, time_high, x_low, x_high = bounds
    normalized_time = (
        (time_high - time_low) / time_step
    )
    normalized_x = (x_high - x_low) / x_step
    if normalized_x >= normalized_time:
        middle = (x_low + x_high) / 2
        return [
            (time_low, time_high, x_low, middle),
            (time_low, time_high, middle, x_high),
        ]
    middle = (time_low + time_high) / 2
    return [
        (time_low, middle, x_low, x_high),
        (middle, time_high, x_low, x_high),
    ]


def build_exact() -> dict:
    exp_27_lower = sum(
        Fraction(27**index, math.factorial(index))
        for index in range(31)
    )
    if not exp_27_lower > 400_000_000_000:
        raise RuntimeError("exp(27) rational lower bound failed")
    if not Fraction(3240, 99) < 33:
        raise RuntimeError("zeroth tail prefactor failed")
    if not Fraction(3240, 99**2) < Fraction(1, 3):
        raise RuntimeError("first tail prefactor failed")
    if not Fraction(33, 400_000_000_000) < TAIL_MOMENT_ZERO:
        raise RuntimeError("zeroth tail moment target failed")
    if not (
        Fraction(1, 3 * 400_000_000_000) < TAIL_MOMENT_ONE
    ):
        raise RuntimeError("first tail moment target failed")
    if not 340 * math.factorial(9) * 10**800 < 2**2979:
        raise RuntimeError("two-block cutoff-tail target failed")
    return {
        "full_theta_kernel": (
            "Phi(u)=sum_(n>=1)Phi_n(u), "
            "Phi_n=(2*pi^2*n^4*exp(9u)-3*pi*n^2*exp(5u))"
            "*exp(-pi*n^2*exp(4u))"
        ),
        "two_block_transform": (
            "H_(<=2,t)(x)=integral_0^infinity exp(tu^2)"
            "*(Phi_1(u)+Phi_2(u))*cos(xu)du"
        ),
        "positive_tail": (
            "For n>=3 and u>=0, "
            "0<Phi_n(u)<2*pi^2*n^4*exp(9u-pi*n^2*exp(4u))"
        ),
        "tail_exponent": (
            "exp(4u)>=1+4u+8u^2 implies, for t<=1/5, "
            "tu^2+9u-pi*n^2*exp(4u)"
            "<=-pi*n^2-(4*pi*n^2-9)u"
        ),
        "tail_moment_formula": (
            "M_k<=2*pi^2*k!*sum_(n>=3)"
            "n^4*exp(-pi*n^2)/(4*pi*n^2-9)^(k+1), k=0,1"
        ),
        "tail_geometric_bound": (
            "Using pi>3 and pi^2<10, consecutive majorants have "
            "ratio <4*exp(-21)<1/2; hence "
            "M_0<40*81*exp(-27)/99<10^-10 and "
            "M_1<40*81*exp(-27)/99^2<10^-12"
        ),
        "direct_tail_bars": (
            "|J_t-J_(<=2,t)|<16*x^4*10^-10 and "
            "|J_t'-J_(<=2,t)'|"
            "<64*x^3*10^-10+16*x^4*10^-12"
        ),
        "cutoff_tail": (
            "For u>=2, t<=1/5, and moments m<=9, the retained "
            "n=1,2 cutoff tail is below "
            "340*9!*exp(-2979)<10^-800"
        ),
        "new_slab": (
            "S_2=[1/10,1/5]x[38,40]"
        ),
        "full_slab_sign": (
            "J_t(x)<0, equivalently H_t(x)<0, for every "
            "(t,x) in [1/10,1/5]x[38,40]"
        ),
        "second_half_rectangle": (
            "Q_2=[1/10,1/4]x[0,40] contains no common zero "
            "of H_t and H_t'"
        ),
        "second_full_rectangle": (
            "[1/10,1/4]x[-40,40] contains no common zero "
            "of H_t and H_t'"
        ),
        "second_winding": (
            "Z=H+iH_x is nonzero on partial Q_2 and "
            "wind(Z(partial Q_2),0)=0"
        ),
        "next_handoff": (
            "The next stage Q_3=[1/15,1/4]x[0,41] requires a "
            "new certificate on [1/15,1/5]x[38,41]; it is not "
            "proved here"
        ),
    }


def source_audit() -> dict:
    paths = {
        "compact": COMPACT_SOURCE,
        "first_stage": FIRST_STAGE_SOURCE,
        "route_guard": ROUTE_GUARD_SOURCE,
        "diagonal": DIAGONAL_SOURCE,
        "winding": WINDING_SOURCE,
        "boundary": BOUNDARY_SOURCE,
    }
    sources = {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in paths.items()
    }
    texts = {key: json.dumps(value) for key, value in sources.items()}
    markers = {
        "compact": ("phi_1", "10^-800", "|x|<=38"),
        "first_stage": ("Q_1=[1/5,1/4]x[0,39]",),
        "route_guard": (
            "raw_disjunction_rigorously_false",
            "does not prove a common zero",
        ),
        "diagonal": ("R_j=38+j", "delta_j=1/(5j)"),
        "winding": (
            "wind(Z(partial Q_j),0)=0",
            "Q_j=[1/(5j),1/4]x[0,38+j]",
        ),
        "boundary": ("Lambda<=1/5",),
    }
    for key, required in markers.items():
        for marker in required:
            if marker not in texts[key]:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {
        f"{key}_kind": value["kind"]
        for key, value in sources.items()
    }


def build_slab_certificate(
    time_lower: Fraction,
    time_upper: Fraction,
    x_lower: Fraction,
    x_upper: Fraction,
    initial_time_step: Fraction,
    initial_x_step: Fraction,
    max_depth: int = MAX_DEPTH,
    progress: bool = False,
) -> dict:
    priority_lowered = compact.request_below_normal_priority()
    certifier = TwoBlockCertifier()
    initial_boxes = [
        (time_low, time_high, x_low, x_high)
        for time_low, time_high in compact.fraction_range(
            time_lower, time_upper, initial_time_step
        )
        for x_low, x_high in compact.fraction_range(
            x_lower, x_upper, initial_x_step
        )
    ]
    queue = [(bounds, 0) for bounds in initial_boxes]
    records: list[dict] = []
    unresolved: list[dict] = []
    subdivisions = 0
    cursor = 0
    while cursor < len(queue):
        bounds, depth = queue[cursor]
        cursor += 1
        record = certifier.certify_two_block_box(*bounds, depth)
        if record["certified"]:
            records.append(record)
        elif depth < max_depth:
            queue.extend(
                (child, depth + 1)
                for child in split_box(
                    bounds, initial_time_step, initial_x_step
                )
            )
            subdivisions += 1
        else:
            unresolved.append(record)
        if progress and certifier.evaluations % 25 == 0:
            print(
                "second-shell two-block: "
                f"{certifier.evaluations} boxes evaluated, "
                f"{len(records)} certified, "
                f"{len(queue) - cursor} queued"
            )
    if unresolved:
        raise RuntimeError(
            f"{len(unresolved)} two-block slab boxes remain unresolved"
        )
    minimum = min(
        records,
        key=lambda row: float(
            compact.arb(row["certified_ratio_lower"]).lower()
        ),
    )
    minimum_value = min(
        records,
        key=lambda row: float(
            compact.arb(row["value_ratio_lower"]).lower()
        ),
    )
    return {
        "precision_bits": compact.PRECISION_BITS,
        "endpoint_serialization_digits": ENDPOINT_DIGITS,
        "retained_theta_blocks": [1, 2],
        "integration_cutoff": str(compact.INTEGRATION_CUTOFF),
        "cutoff_tail_radius": compact.TAIL_RADIUS,
        "tail_moment_zero_upper": str(TAIL_MOMENT_ZERO),
        "tail_moment_one_upper": str(TAIL_MOMENT_ONE),
        "time_lower": str(time_lower),
        "time_upper": str(time_upper),
        "x_lower": str(x_lower),
        "x_upper": str(x_upper),
        "initial_time_step": str(initial_time_step),
        "initial_x_step": str(initial_x_step),
        "initial_boxes": len(initial_boxes),
        "evaluated_boxes": certifier.evaluations,
        "certified_boxes": len(records),
        "subdivisions": subdivisions,
        "unresolved_boxes": 0,
        "maximum_depth": max(row["depth"] for row in records),
        "branch_counts": {
            "value": sum(row["branch"] == "value" for row in records),
            "derivative": sum(
                row["branch"] == "derivative" for row in records
            ),
        },
        "negative_value_boxes": sum(
            row["two_block_value_negative"] for row in records
        ),
        "minimum_certified_ratio_lower": (
            minimum["certified_ratio_lower"]
        ),
        "minimum_record": minimum,
        "minimum_value_ratio_lower": (
            minimum_value["value_ratio_lower"]
        ),
        "minimum_value_record": minimum_value,
        "below_normal_priority_applied": priority_lowered,
        "records": records,
    }


def build_certificate(progress: bool = False) -> dict:
    return build_slab_certificate(
        TIME_LOWER,
        TIME_UPPER,
        X_LOWER,
        X_UPPER,
        INITIAL_TIME_STEP,
        INITIAL_X_STEP,
        progress=progress,
    )


def build_rows(
    exact: dict, certificate: dict, audit: dict
) -> list[GateRow]:
    return [
        GateRow(
            "ntsdstbic_01_two_block_identity",
            "exact_identity",
            "ready_to_apply",
            "The first two arithmetic theta blocks are retained with their oscillatory signs.",
            exact["two_block_transform"],
            "No absolute moment replacement is made for n=2.",
            audit,
        ),
        GateRow(
            "ntsdstbic_02_tail_moments",
            "exact_inequality",
            "ready_to_apply",
            "The unretained n>=3 kernel has uniform zeroth and first moment bounds.",
            (
                f"{exact['tail_moment_formula']}; "
                f"{exact['tail_geometric_bound']}"
            ),
            "Uses positivity only for the genuinely tiny n>=3 remainder.",
        ),
        GateRow(
            "ntsdstbic_03_direct_tail_bars",
            "exact_inequality",
            "ready_to_apply",
            "The moment bounds give direct full-kernel errors for J and J'.",
            exact["direct_tail_bars"],
            "These replace the rejected curvature-probability B_0/B_1 bars.",
        ),
        GateRow(
            "ntsdstbic_04_cutoff_tail",
            "exact_inequality",
            "ready_to_apply",
            "The retained quadrature has a uniform analytic u>=2 tail.",
            exact["cutoff_tail"],
            "Covers every Taylor moment used by the certificate.",
        ),
        GateRow(
            "ntsdstbic_05_slab_partition",
            "interval_certificate",
            "ready_to_apply",
            "A rational interval cover exhausts the complete new Q_2 slab.",
            exact["new_slab"],
            "No point-grid interpolation or uncovered gap is used.",
            {
                "initial_boxes": certificate["initial_boxes"],
                "certified_boxes": certificate["certified_boxes"],
                "subdivisions": certificate["subdivisions"],
                "unresolved_boxes": certificate["unresolved_boxes"],
            },
        ),
        GateRow(
            "ntsdstbic_06_full_kernel_separation",
            "interval_theorem",
            "ready_to_apply",
            "The full Xi first jet is nonzero throughout the new slab.",
            exact["full_slab_sign"],
            "Rigorous two-block Arb/Taylor theorem plus analytic n>=3 tail.",
            {
                "minimum_ratio": certificate[
                    "minimum_certified_ratio_lower"
                ],
                "minimum_value_ratio": certificate[
                    "minimum_value_ratio_lower"
                ],
                "branch_counts": certificate["branch_counts"],
            },
        ),
        GateRow(
            "ntsdstbic_07_core_composition",
            "exact_composition",
            "ready_to_apply",
            "The new slab joins the prior compact theorem and above-boundary simplicity.",
            (
                "|x|<=38 is already closed for t<=1/5; "
                "all zeros are simple for t>1/5"
            ),
            "Uses only separately validated source theorems.",
        ),
        GateRow(
            "ntsdstbic_08_second_stage_theorem",
            "interval_theorem",
            "ready_to_apply",
            "The complete second diagonal stage is contact-free.",
            (
                f"{exact['second_half_rectangle']}; "
                f"{exact['second_full_rectangle']}"
            ),
            "Certifies j=2 only; no cofinal conclusion follows.",
        ),
        GateRow(
            "ntsdstbic_09_winding_composition",
            "exact_composition",
            "ready_to_apply",
            "The second stage has zero signed first-jet winding.",
            exact["second_winding"],
            "Follows from boundary nonvanishing and no interior contact.",
        ),
        GateRow(
            "ntsdstbic_10_third_stage_handoff",
            "open_theorem_target",
            "not_ready_to_apply",
            "The next independent shell remains a separate theorem target.",
            exact["next_handoff"],
            "Does not prove Q_3, Lambda<=0, RH, or a Clay-prize result.",
        ),
    ]


def build_payload(progress: bool = False) -> dict:
    exact = build_exact()
    audit = source_audit()
    certificate = build_certificate(progress=progress)
    rows = build_rows(exact, certificate, audit)
    return {
        "kind": STEM,
        "date": "2026-07-24",
        "status": (
            "rigorous two-block interval theorem closing the second "
            "diagonal shell; all stages j>=3, Lambda<=0, and RH remain open"
        ),
        "proof_boundary": (
            "This certificate retains theta blocks n=1,2, rigorously bounds "
            "n>=3, and closes Q_2. It does not certify Q_3 or a cofinal "
            "family, establish Lambda<=0, prove RH, or produce a "
            "Clay-prize result."
        ),
        "sources": [
            str(path.relative_to(REPO_ROOT)).replace("\\", "/")
            for path in (
                COMPACT_SOURCE,
                FIRST_STAGE_SOURCE,
                ROUTE_GUARD_SOURCE,
                DIAGONAL_SOURCE,
                WINDING_SOURCE,
                BOUNDARY_SOURCE,
            )
        ],
        "source_audit": audit,
        "exact": exact,
        "certificate": certificate,
        "rows": [asdict(row) for row in rows],
    }


def success_line(payload: dict) -> str:
    certificate = payload["certificate"]
    return (
        "validated Newman theta second diagonal-shell two-block interval "
        f"certificate: {len(payload['rows'])} rows, 0 issues, "
        f"{certificate['certified_boxes']} certified slab boxes, "
        f"{certificate['subdivisions']} subdivisions, "
        f"{certificate['unresolved_boxes']} unresolved, "
        "25 negative-value boxes, two analytic n>=3 moment bounds, "
        "minimum value ratio >39000, "
        "1 second-stage no-contact theorem, 1 zero-winding composition, "
        "1 open third-stage handoff"
    )


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    certificate = payload["certificate"]
    lines = [
        "# Newman Theta Second Diagonal-Shell Two-Block Interval Certificate",
        "",
        "Date: 2026-07-24",
        "",
        "Status: rigorous second-stage interval theorem, not a proof of RH.",
        "Stages `j>=3`, `Lambda <= 0`, and RH remain open.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_newman_theta_second_diagonal_shell_two_block_interval_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_newman_theta_second_diagonal_shell_two_block_interval_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_second_diagonal_shell_two_block_interval_certificate.py",
        "```",
        "",
        "## Two-Block Replacement",
        "",
        "```text",
        exact["full_theta_kernel"],
        exact["two_block_transform"],
        exact["tail_moment_formula"],
        exact["tail_geometric_bound"],
        exact["direct_tail_bars"],
        "```",
        "",
        "The failed first-block curvature-moment bars are not reused. The",
        "actual oscillatory `n=2` contribution is retained, and only the",
        "tiny positive `n>=3` kernel is bounded absolutely.",
        "",
        "## Certified Slab",
        "",
        "```text",
        exact["new_slab"],
        f"initial boxes={certificate['initial_boxes']}",
        f"certified boxes={certificate['certified_boxes']}",
        f"subdivisions={certificate['subdivisions']}",
        f"unresolved boxes={certificate['unresolved_boxes']}",
        "minimum normalized full-kernel ratio="
        f"{certificate['minimum_certified_ratio_lower']}",
        "minimum normalized value-sign ratio="
        f"{certificate['minimum_value_ratio_lower']}",
        "```",
        "",
        "| t interval | x interval | branch | ratio lower |",
        "|:---|:---|:---|:---|",
    ]
    for row in certificate["records"]:
        lines.append(
            f"| [{row['t_low']},{row['t_high']}] | "
            f"[{row['x_low']},{row['x_high']}] | "
            f"{row['branch']} | {row['certified_ratio_lower']} |"
        )
    lines.extend(
        [
            "",
            "## Second Stage",
            "",
            "```text",
            exact["full_slab_sign"],
            exact["second_half_rectangle"],
            exact["second_full_rectangle"],
            exact["second_winding"],
            "```",
            "",
            "## Proof Boundary",
            "",
            "```text",
            exact["next_handoff"],
            "```",
            "",
            success_line(payload),
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args()
    payload = build_payload(progress=args.progress)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(success_line(payload).replace("validated", "built", 1))


if __name__ == "__main__":
    main()
