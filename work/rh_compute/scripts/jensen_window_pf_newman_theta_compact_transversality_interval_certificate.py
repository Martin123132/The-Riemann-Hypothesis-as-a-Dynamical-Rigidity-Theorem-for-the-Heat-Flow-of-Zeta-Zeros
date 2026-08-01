#!/usr/bin/env python3
"""Certify theta contact exclusion on the first compact frequency window."""

from __future__ import annotations

import argparse
import ctypes
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
from fractions import Fraction
import json
import math
import os
from pathlib import Path
import sys
from time import perf_counter


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402
from flint import acb, arb  # noqa: E402


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_"
    "compact_transversality_interval_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_theta_"
    "compact_transversality_interval_certificate.md"
)

PRECISION_BITS = 160
TAYLOR_TIME_ORDER = 2
TAYLOR_X_ORDER = 3
INTEGRATION_CUTOFF = Fraction(2)
TAIL_RADIUS = "1e-800"
INITIAL_TIME_STEP = Fraction(1, 50)
INITIAL_X_STEP = Fraction(1, 5)
TIME_LOWER = Fraction(0)
TIME_UPPER = Fraction(1, 5)
X_LOWER = Fraction(1, 4)
X_UPPER = Fraction(38)
MAX_DEPTH = 8


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def request_below_normal_priority() -> bool:
    if os.name != "nt":
        return False
    try:
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.GetCurrentProcess.restype = ctypes.c_void_p
        kernel32.SetPriorityClass.argtypes = [
            ctypes.c_void_p,
            ctypes.c_uint32,
        ]
        kernel32.SetPriorityClass.restype = ctypes.c_int
        handle = kernel32.GetCurrentProcess()
        return bool(kernel32.SetPriorityClass(handle, 0x00004000))
    except (AttributeError, OSError):
        return False


def fraction_decimal(value: Fraction) -> str:
    getcontext().prec = 90
    decimal = Decimal(value.numerator) / Decimal(value.denominator)
    return format(decimal, "f")


def interval_ball(center: Fraction, radius: Fraction) -> arb:
    return arb(
        f"[{fraction_decimal(center)} +/- {fraction_decimal(radius)}]"
    )


def symmetric_ball(radius: Fraction) -> arb:
    return arb(f"[0 +/- {fraction_decimal(radius)}]")


def arb_power(value: arb, exponent: int) -> arb:
    result = arb(1)
    for _ in range(exponent):
        result *= value
    return result


def fraction_range(
    lower: Fraction, upper: Fraction, step: Fraction
) -> list[tuple[Fraction, Fraction]]:
    boxes: list[tuple[Fraction, Fraction]] = []
    left = lower
    while left < upper:
        right = min(left + step, upper)
        boxes.append((left, right))
        left = right
    return boxes


class CompactCertifier:
    def __init__(self) -> None:
        flint.ctx.prec = PRECISION_BITS
        self.pi = acb.pi()
        self.zero = acb(0)
        self.cutoff = acb(fraction_decimal(INTEGRATION_CUTOFF))
        self.abs_tol = arb("1e-35")
        self.tail_real = arb(f"[0 +/- {TAIL_RADIUS}]")
        self.tail_complex = acb(self.tail_real, self.tail_real)
        self.moment_cache: dict[tuple[Fraction, int], arb] = {}
        self.evaluations = 0

    def phi_one(self, u: acb) -> acb:
        return (
            2 * self.pi * self.pi * (9 * u).exp()
            - 3 * self.pi * (5 * u).exp()
        ) * (-self.pi * (4 * u).exp()).exp()

    def oscillatory_integral(
        self, moment: int, time: Fraction, frequency: Fraction
    ) -> acb:
        time_ball = acb(fraction_decimal(time))
        frequency_ball = acb(fraction_decimal(frequency))

        def integrand(u: acb, analytic: bool) -> acb:
            del analytic
            return (
                u**moment
                * (time_ball * u * u).exp()
                * self.phi_one(u)
                * (acb(0, 1) * frequency_ball * u).exp()
            )

        retained = acb.integral(
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

    def moment_upper(self, moment: int, time: Fraction) -> arb:
        key = (time, moment)
        cached = self.moment_cache.get(key)
        if cached is not None:
            return cached
        time_ball = acb(fraction_decimal(time))

        def integrand(u: acb, analytic: bool) -> acb:
            del analytic
            return (
                u**moment
                * (time_ball * u * u).exp()
                * self.phi_one(u)
            )

        retained = acb.integral(
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

    def transform_box(
        self,
        derivative_order: int,
        time_center: Fraction,
        time_radius: Fraction,
        x_center: Fraction,
        x_radius: Fraction,
        integrals: list[acb],
    ) -> arb:
        tau = symmetric_ball(time_radius)
        xi = symmetric_ball(x_radius)
        polynomial = arb(0)
        for time_order in range(TAYLOR_TIME_ORDER + 1):
            for x_order in range(TAYLOR_X_ORDER + 1):
                moment = (
                    derivative_order + 2 * time_order + x_order
                )
                coefficient = (
                    acb(0, 1) ** (derivative_order + x_order)
                    * integrals[moment]
                ).real
                polynomial += (
                    coefficient
                    * arb_power(tau, time_order)
                    * arb_power(xi, x_order)
                    / (
                        math.factorial(time_order)
                        * math.factorial(x_order)
                    )
                )

        time_high = time_center + time_radius
        time_remainder_moment = (
            derivative_order + 2 * (TAYLOR_TIME_ORDER + 1)
        )
        time_remainder = (
            arb_power(arb(fraction_decimal(time_radius)), 3)
            / math.factorial(3)
            * self.moment_upper(
                time_remainder_moment, time_high
            ).upper()
        )

        x_remainder = arb(0)
        for time_order in range(TAYLOR_TIME_ORDER + 1):
            moment = (
                derivative_order
                + 2 * time_order
                + TAYLOR_X_ORDER
                + 1
            )
            x_remainder += (
                arb_power(arb(fraction_decimal(time_radius)), time_order)
                / math.factorial(time_order)
                * arb_power(
                    arb(fraction_decimal(x_radius)),
                    TAYLOR_X_ORDER + 1,
                )
                / math.factorial(TAYLOR_X_ORDER + 1)
                * self.moment_upper(moment, time_center).upper()
            )

        total_remainder = (time_remainder + x_remainder).upper()
        return polynomial + arb(0, str(total_remainder))

    def certify_box(
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

        max_integral_moment = (
            1
            + 2 * TAYLOR_TIME_ORDER
            + TAYLOR_X_ORDER
        )
        integrals = [
            self.oscillatory_integral(moment, time_center, x_center)
            for moment in range(max_integral_moment + 1)
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

        time_box = interval_ball(time_center, time_radius)
        x_box = interval_ball(x_center, x_radius)
        x_two = x_box * x_box
        x_three = x_two * x_box
        x_four = x_three * x_box
        time_two = time_box * time_box
        j_box = 16 * x_four * h_box
        j_prime_box = (
            64 * x_three * h_box + 16 * x_four * h_prime_box
        )

        delta_zero = arb(1) / 2800
        delta_one = arb(1) / 137200
        delta_two = arb(1) / 3361400
        delta_three = arb(1) / 54900000
        p_poly = (
            x_four
            + (6 * time_box + 1) * x_two
            + 24 * time_two
        )
        r_poly = (6 * time_box + 1) * x_two + 24 * time_two
        b_zero = (
            4 * time_two * x_two * delta_two
            + 4
            * time_box
            * x_box
            * (x_two + 4 * time_box)
            * delta_one
            + (x_four + 2 * r_poly) * delta_zero
        )
        b_one = (
            4 * time_two * x_two * delta_three
            + 4
            * time_box
            * x_box
            * (x_two + 2 * time_box)
            * delta_two
            + abs(
                p_poly
                - 4 * time_box * (3 * x_two + 4 * time_box)
            )
            * delta_one
            + (
                4 * x_three
                + 4 * (6 * time_box + 1) * x_box
            )
            * delta_zero
        )

        value_lower = j_box.abs_lower()
        derivative_lower = j_prime_box.abs_lower()
        b_zero_upper = b_zero.upper()
        b_one_upper = b_one.upper()
        value_ratio = value_lower / b_zero_upper
        derivative_ratio = derivative_lower / b_one_upper
        value_pass = value_lower > b_zero_upper
        derivative_pass = derivative_lower > b_one_upper

        if value_pass or derivative_pass:
            if value_ratio > derivative_ratio:
                branch = "value"
                ratio = value_ratio
            else:
                branch = "derivative"
                ratio = derivative_ratio
            return {
                "certified": True,
                "depth": depth,
                "t_low": str(time_low),
                "t_high": str(time_high),
                "x_low": str(x_low),
                "x_high": str(x_high),
                "branch": branch,
                "certified_ratio_lower": str(ratio.lower()),
                "value_ratio_lower": str(value_ratio.lower()),
                "derivative_ratio_lower": str(
                    derivative_ratio.lower()
                ),
            }

        return {
            "certified": False,
            "depth": depth,
            "t_low": str(time_low),
            "t_high": str(time_high),
            "x_low": str(x_low),
            "x_high": str(x_high),
            "value_ratio_lower": str(value_ratio.lower()),
            "derivative_ratio_lower": str(
                derivative_ratio.lower()
            ),
        }


def split_box(
    bounds: tuple[Fraction, Fraction, Fraction, Fraction],
) -> list[tuple[Fraction, Fraction, Fraction, Fraction]]:
    time_low, time_high, x_low, x_high = bounds
    normalized_time_width = (
        (time_high - time_low) / INITIAL_TIME_STEP
    )
    normalized_x_width = (x_high - x_low) / INITIAL_X_STEP
    if normalized_x_width >= normalized_time_width:
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
    if not Fraction(248, 371925) == (
        Fraction(9, 2900) - Fraction(5, 2052)
    ):
        raise RuntimeError("near-origin overlap margin changed")
    exp_eight_lower = sum(
        Fraction(8**k, math.factorial(k)) for k in range(8)
    )
    if not exp_eight_lower > 1000:
        raise RuntimeError("exp(8) rational lower bound failed")
    if not Fraction(22, 7) ** 2 < 10:
        raise RuntimeError("pi-squared rational upper bound failed")
    if not 20 * math.factorial(9) < 8_000_000:
        raise RuntimeError("tail prefactor bound failed")
    if not 2**10 > 10**3:
        raise RuntimeError("binary-decimal tail conversion failed")
    if not (
        20 * math.factorial(9) * 10**800 < 2**2979
    ):
        raise RuntimeError("tail radius integer comparison failed")
    return {
        "contact_implication": (
            "|J_1|>B_0 or |J_1'|>B_1 implies "
            "(J_t,J_t')!=(0,0)"
        ),
        "first_component_transform": (
            "H_(1,t)(x)=integral_0^infinity "
            "exp(tu^2)*phi_1(u)*cos(xu)du"
        ),
        "contact_formulas": (
            "J_1=16x^4H_(1,t), "
            "J_1'=64x^3H_(1,t)+16x^4H_(1,t)'"
        ),
        "tail_proof": {
            "range": "u>=2, 0<=t<=1/5, 0<=m<=9",
            "rational_audit": (
                "sum_(k=0)^7 8^k/k!>1000, pi>3, pi^2<10, "
                "20*9!*10^800<2^2979"
            ),
            "exponent": (
                "t*u^2+9u-pi*exp(4u)<-2979-u"
            ),
            "moment_tail": (
                "integral_2^infinity u^m*exp(tu^2)*phi_1(u)du"
                "<20*9!*exp(-2979)<10^-800"
            ),
            "radius_used": TAIL_RADIUS,
        },
        "taylor_enclosure": {
            "time_order": TAYLOR_TIME_ORDER,
            "x_order": TAYLOR_X_ORDER,
            "coefficient": (
                "partial_t^a partial_x^b H="
                "Re[i^b integral u^(2a+b)exp(tu^2)"
                "*phi_1(u)exp(ixu)du]"
            ),
            "time_remainder": (
                "h_t^(A+1)/(A+1)! times "
                "M_(q+2A+2)(t_high)"
            ),
            "x_remainder": (
                "sum_(a=0)^A h_t^a/a! * h_x^(B+1)/(B+1)! "
                "*M_(q+2a+B+1)(t_center)"
            ),
        },
        "proved_rectangle": (
            "0<=t<=1/5 and 1/4<=x<=38"
        ),
        "origin_overlap": (
            "The separate exact moment lemma proves H_t(x)>0 "
            "for |x|<=1/4 with margin 248/371925."
        ),
    }


def build_certificate(progress: bool = False) -> dict:
    priority_lowered = request_below_normal_priority()
    certifier = CompactCertifier()
    initial_boxes = [
        (time_low, time_high, x_low, x_high)
        for time_low, time_high in fraction_range(
            TIME_LOWER, TIME_UPPER, INITIAL_TIME_STEP
        )
        for x_low, x_high in fraction_range(
            X_LOWER, X_UPPER, INITIAL_X_STEP
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
        record = certifier.certify_box(*bounds, depth)
        if record["certified"]:
            records.append(record)
        elif depth < MAX_DEPTH:
            children = split_box(bounds)
            queue.extend((child, depth + 1) for child in children)
            subdivisions += 1
        else:
            unresolved.append(record)
        if progress and certifier.evaluations % 250 == 0:
            print(
                "compact transversality: "
                f"{certifier.evaluations} boxes evaluated, "
                f"{len(records)} certified, "
                f"{len(queue) - cursor} queued"
            )

    if unresolved:
        raise RuntimeError(
            f"{len(unresolved)} compact boxes remain uncertified"
        )

    minimum_record = min(
        records,
        key=lambda row: float(
            arb(row["certified_ratio_lower"]).lower()
        ),
    )
    branch_counts = {
        "value": sum(row["branch"] == "value" for row in records),
        "derivative": sum(
            row["branch"] == "derivative" for row in records
        ),
    }
    return {
        "precision_bits": PRECISION_BITS,
        "integration_cutoff": str(INTEGRATION_CUTOFF),
        "tail_radius": TAIL_RADIUS,
        "initial_time_step": str(INITIAL_TIME_STEP),
        "initial_x_step": str(INITIAL_X_STEP),
        "initial_boxes": len(initial_boxes),
        "evaluated_boxes": certifier.evaluations,
        "certified_boxes": len(records),
        "subdivisions": subdivisions,
        "unresolved_boxes": 0,
        "maximum_depth": max(row["depth"] for row in records),
        "branch_counts": branch_counts,
        "minimum_certified_ratio_lower": (
            minimum_record["certified_ratio_lower"]
        ),
        "minimum_record": minimum_record,
        "below_normal_priority_requested": True,
        "below_normal_priority_applied": priority_lowered,
        "records": records,
    }


def build_rows(exact: dict, certificate: dict) -> list[GateRow]:
    return [
        GateRow(
            "ntctic_01_contact_implication",
            "exact_implication",
            "ready_to_apply",
            "The componentwise bars convert either strict first-block inequality into contact exclusion.",
            exact["contact_implication"],
            "Imported from the independently checked theta/operator gate.",
        ),
        GateRow(
            "ntctic_02_first_component_transform",
            "exact_identity",
            "ready_to_apply",
            "The first contact block is the first positive Xi-kernel transform.",
            exact["contact_formulas"],
            "Follows from the theta primitive identity term by term.",
        ),
        GateRow(
            "ntctic_03_retained_integral_tail",
            "exact_inequality",
            "ready_to_apply",
            "The omitted u>2 tail is rigorously below the attached Arb radius.",
            exact["tail_proof"]["moment_tail"],
            exact["tail_proof"]["range"],
        ),
        GateRow(
            "ntctic_04_bivariate_taylor_enclosure",
            "exact_enclosure",
            "ready_to_apply",
            "Mixed transform jets and positive moments enclose every parameter box.",
            exact["taylor_enclosure"]["coefficient"],
            "Second order in time and third order in frequency, with explicit remainders.",
        ),
        GateRow(
            "ntctic_05_partition_integrity",
            "interval_certificate",
            "ready_to_apply",
            "The adaptive rational box partition covers the full compact rectangle.",
            (
                f"{certificate['certified_boxes']} certified boxes, "
                f"{certificate['unresolved_boxes']} unresolved"
            ),
            exact["proved_rectangle"],
        ),
        GateRow(
            "ntctic_06_compact_contact_exclusion",
            "interval_theorem",
            "ready_to_apply",
            "The Xi heat flow has no multiple real zero on the certified compact rectangle.",
            (
                "|J_1|>B_0 or |J_1'|>B_1 on every box; "
                "(H_t,H_t')!=(0,0)"
            ),
            exact["proved_rectangle"],
        ),
        GateRow(
            "ntctic_07_origin_composition",
            "exact_composition",
            "ready_to_apply",
            "The compact interval theorem overlaps the exact origin-positivity theorem.",
            exact["origin_overlap"],
            "Together they prove no contact for |x|<=38.",
        ),
        GateRow(
            "ntctic_08_high_frequency_handoff",
            "open_theorem_target",
            "open",
            "Only frequencies beyond the certified compact window remain in this route.",
            (
                "Prove (H_t,H_t')!=(0,0) for x>38 and "
                "0<t<=1/5 using the corrected Riemann-Siegel partition."
            ),
            "No claim beyond x=38 is made by this interval certificate.",
        ),
    ]


def build_payload(progress: bool = False) -> dict:
    exact = build_exact()
    certificate = build_certificate(progress=progress)
    rows = build_rows(exact, certificate)
    return {
        "kind": (
            "jensen_window_pf_newman_theta_"
            "compact_transversality_interval_certificate"
        ),
        "date": "2026-07-24",
        "status": (
            "rigorous Arb/Taylor contact exclusion for |x|<=38 and "
            "0<=t<=1/5; the high-frequency transversality target, "
            "Lambda<=0, and RH remain open"
        ),
        "proof_boundary": (
            "The exact origin theorem covers |x|<=1/4 and this certificate "
            "covers 1/4<=|x|<=38 by evenness. It does not certify x>38 or "
            "the residual corrected Riemann-Siegel phase layer."
        ),
        "exact": exact,
        "certificate": certificate,
        "rows": [asdict(row) for row in rows],
    }


def success_line(payload: dict) -> str:
    certificate = payload["certificate"]
    return (
        "validated Newman theta compact-transversality interval certificate: "
        "8 rows, 0 issues, 4 exact identities/inequalities, "
        f"{certificate['certified_boxes']} certified boxes, "
        f"{certificate['subdivisions']} adaptive subdivisions, "
        f"{certificate['unresolved_boxes']} unresolved boxes, "
        "1 compact no-contact theorem, 1 origin composition, "
        "1 open high-frequency handoff"
    )


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    certificate = payload["certificate"]
    minimum = certificate["minimum_record"]
    return "\n".join(
        [
            "# Newman Theta Compact-Transversality Interval Certificate",
            "",
            "Date: 2026-07-24",
            "",
            "Status: rigorous compact contact exclusion. This is not a proof",
            "of `Lambda <= 0` or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_theta_compact_transversality_interval_certificate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_theta_compact_transversality_interval_certificate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_compact_transversality_interval_certificate.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            success_line(payload),
            "```",
            "",
            "## Certified Theorem",
            "",
            "For every",
            "",
            "```text",
            "0<=t<=1/5 and |x|<=38",
            "```",
            "",
            "the exact Newman heat-flow transform satisfies",
            "",
            "```text",
            "(H_t(x),H_t'(x))!=(0,0).",
            "```",
            "",
            "For `|x|<=1/4` this follows from the independent rational",
            "moment margin `248/371925`. Evenness reduces the remaining",
            "region to the certified rectangle `1/4<=x<=38`.",
            "",
            "## Interval Method",
            "",
            exact["first_component_transform"],
            "",
            exact["contact_formulas"],
            "",
            "The retained integral is evaluated by Arb on `0<=u<=2`.",
            "For every transform moment used in the Taylor model,",
            "",
            "```text",
            exact["tail_proof"]["exponent"],
            exact["tail_proof"]["moment_tail"],
            "```",
            "",
            "The omitted tail is therefore included with radius `1e-800`.",
            "A second-order time and third-order frequency Taylor model",
            "uses exact mixed-jet integrals and positive moment remainders.",
            "",
            "## Partition",
            "",
            "```text",
            f"precision bits={certificate['precision_bits']}",
            f"initial boxes={certificate['initial_boxes']}",
            f"evaluated boxes={certificate['evaluated_boxes']}",
            f"certified leaf boxes={certificate['certified_boxes']}",
            f"adaptive subdivisions={certificate['subdivisions']}",
            f"maximum depth={certificate['maximum_depth']}",
            f"unresolved boxes={certificate['unresolved_boxes']}",
            (
                "value/derivative branches="
                f"{certificate['branch_counts']['value']}/"
                f"{certificate['branch_counts']['derivative']}"
            ),
            (
                "minimum certified ratio lower="
                f"{certificate['minimum_certified_ratio_lower']}"
            ),
            (
                "minimum box="
                f"t:[{minimum['t_low']},{minimum['t_high']}], "
                f"x:[{minimum['x_low']},{minimum['x_high']}], "
                f"branch={minimum['branch']}"
            ),
            "```",
            "",
            "Every leaf proves `|J_1|>B_0` or `|J_1'|>B_1`; the exact",
            "componentwise tail theorem then excludes full contact.",
            "",
            "## Proof Boundary",
            "",
            "No value beyond `|x|=38` is certified here. The next theorem",
            "must use the existing corrected Riemann-Siegel partition:",
            "the dominant and oscillatory-zeta regions are already closed,",
            "while the scaled critical phase-transversality layer remains",
            "open. The compact theorem does not imply `Lambda<=0` or RH.",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args()

    started = perf_counter()
    payload = build_payload(progress=args.progress)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Newman theta compact-transversality interval certificate: "
        f"{payload['certificate']['certified_boxes']} boxes in "
        f"{perf_counter() - started:.1f}s"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
