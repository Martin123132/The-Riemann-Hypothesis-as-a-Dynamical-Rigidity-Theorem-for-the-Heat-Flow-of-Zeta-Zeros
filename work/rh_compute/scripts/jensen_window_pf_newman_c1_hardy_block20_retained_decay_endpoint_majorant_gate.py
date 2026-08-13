#!/usr/bin/env python3
"""Sharpen endpoint majorants by retaining common-ray homotopy decay."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Callable


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_block20_transported_endpoint_budget_gate as transport
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import acb, arb, ctx


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_retained_decay_endpoint_majorant_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_retained_decay_endpoint_majorant_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_block20_retained_decay_endpoint_majorant_gate.py"
PRIOR = transport.MAJORANT
WEIGHT_FIXTURE = transport.WEIGHT_FIXTURE
WEIGHTS = transport.WEIGHTS
PRECISIONS = (60, 90)
ORIGIN = Fraction(1, 10**12)
CUTOFF = Fraction(64)
REQUESTED_ERROR_SCALE = Fraction(5, 1000)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def upper(value: arb) -> Fraction:
    return cells.bound_fraction(value.upper())


def lower(value: arb) -> Fraction:
    return cells.bound_fraction(value.lower())


def upper_abs(value: acb) -> Fraction:
    return upper(abs(value))


def error_ball(radius: Fraction) -> arb:
    return arb(0) if radius == 0 else cells.interval_ball(Fraction(0), radius)


def integrate_real(
    function: Callable[[acb, bool], acb], lower_limit: Fraction, upper_limit: Fraction, dps: int
) -> arb:
    ctx.dps = dps
    ctx.threads = 1
    tolerance = arb(f"1e-{dps - 20}")
    value = acb.integral(
        function,
        acb(cells.fraction_ball(lower_limit)),
        acb(cells.fraction_ball(upper_limit)),
        abs_tol=tolerance,
        rel_tol=tolerance,
        deg_limit=64,
        eval_limit=200000,
        depth_limit=64,
    )
    require(value.imag.contains(0), "real majorant quadrature acquired an imaginary part")
    return value.real


def ray_sines(family: str, phi3_sign: int) -> tuple[arb, arb, arb, str]:
    root_two = arb(2).sqrt()
    root_three = arb(3).sqrt()
    root_six = arb(6).sqrt()
    if family == "b" and phi3_sign >= 0:
        return (root_six + root_two) / 4, arb(1) / 2, root_two / 2, "7*pi/12"
    if family == "b":
        return root_two / 2, arb(1), root_two / 2, "3*pi/4"
    if family == "c" and phi3_sign >= 0:
        return arb(1) / 2, root_three / 2, arb(1), "pi/6"
    return (root_six + root_two) / 4, arb(1) / 2, root_two / 2, "5*pi/12"


def exponential_moment_two(rate: arb, cutoff: arb) -> arb:
    exponential = (-rate * cutoff).exp()
    return exponential * (cutoff**2 / rate + 2 * cutoff / rate**2 + 2 / rate**3)


def exponential_moment_three(rate: arb, cutoff: arb) -> arb:
    exponential = (-rate * cutoff).exp()
    return exponential * (
        cutoff**3 / rate + 3 * cutoff**2 / rate**2 + 6 * cutoff / rate**3 + 6 / rate**4
    )


def generic_retained_decay_signed(
    quadratic_fraction: Fraction,
    cubic_signed_fraction: Fraction,
    gap_fraction: Fraction,
    family: str,
    dps: int,
) -> tuple[arb, Fraction, Fraction, str]:
    require(quadratic_fraction > 0 and gap_fraction > 0, "generic retained-decay domain failure")
    ctx.dps = dps
    quadratic = cells.fraction_ball(quadratic_fraction)
    cubic = cells.fraction_ball(abs(cubic_signed_fraction))
    gap = cells.fraction_ball(gap_fraction)
    sine_one, sine_two, sine_three, angle = ray_sines(
        family, 1 if cubic_signed_fraction >= 0 else -1
    )
    pi = arb.pi()
    alpha = 2 * pi * sine_one
    beta = 2 * pi * quadratic * sine_two
    gamma = 2 * pi * cubic * sine_three
    epsilon = cells.fraction_ball(ORIGIN)
    cutoff = cells.fraction_ball(CUTOFF)

    def integrand(radius: acb, _analytic: bool) -> acb:
        homotopy_decay = beta * radius**2 + gamma * radius**3
        homotopy_factor = (1 - (-homotopy_decay).exp()) / homotopy_decay
        geometric_denominator = 1 - (-alpha * radius).exp()
        return (
            2
            * pi
            * (quadratic * radius**2 + cubic * radius**3)
            * (-alpha * gap * radius).exp()
            * homotopy_factor
            / geometric_denominator
        )

    central = integrate_real(integrand, ORIGIN, CUTOFF, dps)
    origin_radius = upper(
        2
        * pi
        * (alpha * epsilon).exp()
        / alpha
        * (quadratic * epsilon**2 / 2 + cubic * epsilon**3 / 3)
    )
    denominator_floor = 1 - (-alpha * cutoff).exp()
    rate = alpha * gap
    tail_radius = upper(
        2
        * pi
        / denominator_floor
        * (
            quadratic * exponential_moment_two(rate, cutoff)
            + cubic * exponential_moment_three(rate, cutoff)
        )
    )
    return central + error_ball(origin_radius + tail_radius), origin_radius, tail_radius, angle


def gaussian_moment_two(rate: arb, cutoff: arb) -> arb:
    root_rate = rate.sqrt()
    exponential = (-rate * cutoff**2).exp()
    return (
        cutoff * exponential / (2 * rate)
        + arb.pi().sqrt() * (root_rate * cutoff).erfc() / (4 * root_rate**3)
    )


def gaussian_moment_three(rate: arb, cutoff: arb) -> arb:
    exponential = (-rate * cutoff**2).exp()
    return exponential * (1 + rate * cutoff**2) / (2 * rate**2)


def exceptional_b_retained_linear(
    phi2_fraction: Fraction,
    phi3_fraction: Fraction,
    length: int,
    delta_fraction: Fraction,
    dps: int,
) -> tuple[arb, Fraction, str]:
    require(phi2_fraction > 0 and delta_fraction >= 0, "W2 retained-linear domain failure")
    local_quadratic = phi2_fraction + 3 * phi3_fraction * length
    require(local_quadratic > 0, "W2 local quadratic is nonpositive")
    ctx.dps = dps
    sine_one, _, _, angle = ray_sines("b", 1 if phi3_fraction >= 0 else -1)
    sigma_two = cells.fraction_ball(Fraction(1, 2) if phi3_fraction >= 0 else Fraction(1))
    floor = phi2_fraction if phi3_fraction >= 0 else local_quadratic
    pi = arb.pi()
    cubic = cells.fraction_ball(abs(phi3_fraction))
    linear_rate = 2 * pi * sine_one * cells.fraction_ball(delta_fraction)
    quadratic_rate = 2 * pi * cells.fraction_ball(floor) * sigma_two

    def integrand(radius: acb, _analytic: bool) -> acb:
        return (
            2
            * pi
            * cubic
            * (3 * length * radius**2 + radius**3)
            * (-linear_rate * radius - quadratic_rate * radius**2).exp()
        )

    central = integrate_real(integrand, Fraction(0), CUTOFF, dps)
    cutoff = cells.fraction_ball(CUTOFF)
    tail_radius = upper(
        2
        * pi
        * cubic
        * (
            3 * length * gaussian_moment_two(quadratic_rate, cutoff)
            + gaussian_moment_three(quadratic_rate, cutoff)
        )
    )
    return central + error_ball(tail_radius), tail_radius, angle


def exceptional_c_retained_linear(
    phi2_fraction: Fraction,
    phi3_fraction: Fraction,
    eta_fraction: Fraction,
    dps: int,
) -> tuple[arb, Fraction, str]:
    require(phi2_fraction > 0 and eta_fraction >= 0, "W3 retained-linear domain failure")
    ctx.dps = dps
    sine_one, _, _, angle = ray_sines("c", 1 if phi3_fraction >= 0 else -1)
    pi = arb.pi()
    cubic = cells.fraction_ball(abs(phi3_fraction))
    linear_rate = 2 * pi * sine_one * cells.fraction_ball(eta_fraction)
    quadratic_rate = pi * cells.fraction_ball(phi2_fraction)

    def integrand(radius: acb, _analytic: bool) -> acb:
        return 2 * pi * cubic * radius**3 * (-linear_rate * radius - quadratic_rate * radius**2).exp()

    central = integrate_real(integrand, Fraction(0), CUTOFF, dps)
    cutoff = cells.fraction_ball(CUTOFF)
    tail_radius = upper(2 * pi * cubic * gaussian_moment_three(quadratic_rate, cutoff))
    return central + error_ball(tail_radius), tail_radius, angle


def evaluate_chain(data: dict[str, Any], dps: int) -> dict[str, Any]:
    parent = data["levels"][1]
    phi1, phi2, phi3 = (
        cells.binary128_fraction(value) for value in parent["coefficients_hex"]
    )
    length = int(parent["length"])
    child_length = int(data["levels"][2]["length"])
    xi = phi1 + 2 * phi2 * length + 3 * phi3 * length**2
    endpoint_quadratic = phi2 + 3 * phi3 * length
    mode_b = child_length + 1
    delta = Fraction(mode_b) - xi
    eta = 1 + phi1
    require(delta > 0 and eta > 0, "exceptional selector gap is nonpositive")

    b_exceptional = exceptional_b_retained_linear(phi2, phi3, length, delta, dps)
    b_zero = generic_retained_decay_signed(phi2, phi3, Fraction(mode_b) - phi1, "b", dps)
    b_length = generic_retained_decay_signed(
        endpoint_quadratic, phi3, Fraction(mode_b + 1) - xi, "b", dps
    )
    if phi1 < 0:
        c_exceptional_gap = eta
        c_zero_gap = 2 + phi1
        c_length_gap = 1 + xi
        selector = "W3"
    else:
        c_exceptional_gap = phi1
        c_zero_gap = 1 + phi1
        c_length_gap = xi
        selector = "W4"
    require(c_exceptional_gap > 0, "selected lower exceptional gap is nonpositive")
    c_exceptional = exceptional_c_retained_linear(
        phi2, phi3, c_exceptional_gap, dps
    )
    c_zero = generic_retained_decay_signed(phi2, phi3, c_zero_gap, "c", dps)
    c_length = generic_retained_decay_signed(endpoint_quadratic, phi3, c_length_gap, "c", dps)
    upper_family = b_exceptional[0] + b_zero[0] + b_length[0]
    lower_family = c_exceptional[0] + c_zero[0] + c_length[0]
    return {
        "selector": selector,
        "delta": delta,
        "eta": eta,
        "c_exceptional_gap": c_exceptional_gap,
        "b_exceptional": b_exceptional,
        "b_zero": b_zero,
        "b_length": b_length,
        "c_exceptional": c_exceptional,
        "c_zero": c_zero,
        "c_length": c_length,
        "upper_family": upper_family,
        "lower_family": lower_family,
        "complete": upper_family + lower_family,
    }


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 retained-decay endpoint majorant gate

Date: 2026-08-07

Status: rigorous finite-roster application of a selector-gap-finite retained-decay theorem; not a proof of the complete evaluator or RH

## Summed generic modes

For one endpoint family choose the same sign-aware ray as in Sections
11.270--11.272.  Let `alpha=2*pi*sin(theta)`,
`beta=2*pi*B*|sin(2 theta)|`, `gamma=2*pi*C*|sin(3 theta)|`, and

```text
phi(x)=(1-exp(-x))/x,  phi(0)=1.
```

Integrating the linear-to-cubic homotopy parameter before applying the
triangle inequality and then summing the modes geometrically gives

```text
S(B,C,q)=2*pi integral_0^infinity
 (B*r^2+C*r^3) exp(-alpha*q*r)
 phi(beta*r^2+gamma*r^3)/(1-exp(-alpha*r)) dr.       (1)
```

This is no larger than the Hurwitz bound in Section 11.272, but it retains
all common-ray quadratic and cubic decay.  Its apparent origin singularity is
removable.  The gate integrates from `10^-12` to `64`, bounds the omitted
origin explicitly, and bounds infinity by exact exponential moments.

For W2 and W3 the reciprocal boundary has already cancelled.  Retaining the
nonnegative selector-gap decay gives

```text
E_b=2*pi*C integral (3*N*r^2+r^3) exp(-alpha*delta*r-a_b*r^2) dr,
E_c=2*pi*C integral r^3 exp(-alpha*eta*r-a_c*r^2) dr. (2)
```

Both remain finite at `delta=0` or `eta=0`; no reciprocal selector gap has
been reintroduced.

## Finite block result

All 374 calls overlap across the 60/90-digit Arb ladder and dominate their
independently enclosed upper and lower endpoint corrections.  Transport
through the exact source-observed block weights gives

```text
maximum prior transported bound       {aggregate['maximum_prior_transported_majorant_upper']}
maximum retained-decay bound           {aggregate['maximum_retained_decay_transported_majorant_upper']}
requested input scale                  {aggregate['requested_error_scale']}
maximum retained bound/scale ratio     {aggregate['maximum_retained_decay_to_requested_scale_ratio']}
outputs below requested scale          {aggregate['outputs_within_requested_scale']} / 15
minimum improvement factor             {aggregate['minimum_transport_improvement_factor']}
```

This closes the **local endpoint-model portion** of the finite block-20
nominal budget without using observed signed cancellation.  It does not close
coefficient transport, recurrence arithmetic, Legendre truncation, other
blocks, the outer Hardy representation, or a height-uniform theorem.
"""


def main() -> int:
    for path in (PRIOR, WEIGHT_FIXTURE, WEIGHTS, CHECKER):
        require(path.is_file(), f"missing retained-decay dependency: {path}")
    recursive = native_q.load_recursive_chains()
    prior = json.loads(PRIOR.read_text(encoding="utf-8"))
    prior_rows = {int(row["chain"]): row for row in prior["rows"]}
    require(set(prior_rows) == set(recursive), "prior/recursive roster mismatch")

    rows: list[dict[str, Any]] = []
    evaluated_high: dict[int, dict[str, Any]] = {}
    maximum_origin = Fraction(0)
    maximum_tail = Fraction(0)
    for chain in sorted(recursive):
        low = evaluate_chain(recursive[chain], PRECISIONS[0])
        high = evaluate_chain(recursive[chain], PRECISIONS[1])
        for name in ("upper_family", "lower_family", "complete"):
            require(low[name].overlaps(high[name]), f"chain {chain} {name} precision nonoverlap")
        upper_lower, upper_upper = lower(high["upper_family"]), upper(high["upper_family"])
        lower_lower, lower_upper = lower(high["lower_family"]), upper(high["lower_family"])
        complete_lower, complete_upper = lower(high["complete"]), upper(high["complete"])
        prior_row = prior_rows[chain]
        require(
            upper_lower > Fraction(prior_row["upper_family_actual_upper"]),
            f"chain {chain} sharpened upper majorant failure",
        )
        require(
            lower_lower > Fraction(prior_row["lower_family_actual_upper"]),
            f"chain {chain} sharpened lower majorant failure",
        )
        require(
            complete_lower > Fraction(prior_row["correction_actual_upper"]),
            f"chain {chain} sharpened complete majorant failure",
        )
        origins = [high[name][1] for name in ("b_zero", "b_length", "c_zero", "c_length")]
        tails = [
            high["b_exceptional"][1],
            high["c_exceptional"][1],
            high["b_zero"][2],
            high["b_length"][2],
            high["c_zero"][2],
            high["c_length"][2],
        ]
        maximum_origin = max(maximum_origin, *origins)
        maximum_tail = max(maximum_tail, *tails)
        evaluated_high[chain] = high
        header = recursive[chain]["header"]
        rows.append(
            {
                "chain": chain,
                "sum_index": int(header["sum_index"]),
                "branch": int(header["branch"]),
                "selector": high["selector"],
                "delta": decimal(high["delta"]),
                "eta": decimal(high["eta"]),
                "c_exceptional_gap": decimal(high["c_exceptional_gap"]),
                "b_exceptional_upper": decimal(upper(high["b_exceptional"][0])),
                "b_zero_generic_upper": decimal(upper(high["b_zero"][0])),
                "b_length_generic_upper": decimal(upper(high["b_length"][0])),
                "c_exceptional_upper": decimal(upper(high["c_exceptional"][0])),
                "c_zero_generic_upper": decimal(upper(high["c_zero"][0])),
                "c_length_generic_upper": decimal(upper(high["c_length"][0])),
                "upper_family_majorant_lower": decimal(upper_lower),
                "upper_family_majorant_upper": decimal(upper_upper),
                "lower_family_majorant_lower": decimal(lower_lower),
                "lower_family_majorant_upper": decimal(lower_upper),
                "complete_majorant_lower": decimal(complete_lower),
                "complete_majorant_upper": decimal(complete_upper),
                "prior_complete_majorant_upper": prior_row["complete_majorant_upper"],
                "complete_improvement_factor_lower": decimal(
                    Fraction(prior_row["complete_majorant_upper"]) / complete_upper
                ),
                "precision_overlap": True,
            }
        )

    weights = transport.load_weights()
    output_rows: list[dict[str, Any]] = []
    for output_index in range(1, 16):
        prior_total = Fraction(0)
        sharp_total = Fraction(0)
        for row in rows:
            chain = int(row["chain"])
            weight = weights[(int(row["sum_index"]), output_index, int(row["branch"]))]
            factor = weight["amplitude_fraction"] * weight["phase_abs_upper"]
            prior_total += factor * Fraction(row["prior_complete_majorant_upper"])
            sharp_total += factor * Fraction(row["complete_majorant_upper"])
        output_rows.append(
            {
                "output_index": output_index,
                "output_label": transport.output_label(output_index),
                "prior_transported_majorant_upper": decimal(prior_total),
                "retained_decay_transported_majorant_upper": decimal(sharp_total),
                "transport_improvement_factor": decimal(prior_total / sharp_total),
                "retained_decay_to_requested_scale_ratio": decimal(
                    sharp_total / REQUESTED_ERROR_SCALE
                ),
                "within_requested_scale": sharp_total < REQUESTED_ERROR_SCALE,
            }
        )

    max_prior = max((Fraction(row["prior_transported_majorant_upper"]), row["output_index"]) for row in output_rows)
    max_sharp = max((Fraction(row["retained_decay_transported_majorant_upper"]), row["output_index"]) for row in output_rows)
    max_ratio = max((Fraction(row["retained_decay_to_requested_scale_ratio"]), row["output_index"]) for row in output_rows)
    min_improvement = min((Fraction(row["transport_improvement_factor"]), row["output_index"]) for row in output_rows)
    aggregate = {
        "call_count": len(rows),
        "precision_overlap_count": sum(bool(row["precision_overlap"]) for row in rows),
        "w3_count": sum(row["selector"] == "W3" for row in rows),
        "w4_count": sum(row["selector"] == "W4" for row in rows),
        "maximum_origin_radius_upper": decimal(maximum_origin),
        "maximum_tail_radius_upper": decimal(maximum_tail),
        "requested_error_scale": decimal(REQUESTED_ERROR_SCALE),
        "maximum_prior_transported_majorant_upper": decimal(max_prior[0]),
        "maximum_prior_transported_majorant_witness": max_prior[1],
        "maximum_retained_decay_transported_majorant_upper": decimal(max_sharp[0]),
        "maximum_retained_decay_transported_majorant_witness": max_sharp[1],
        "maximum_retained_decay_to_requested_scale_ratio": decimal(max_ratio[0]),
        "maximum_retained_decay_to_requested_scale_ratio_witness": max_ratio[1],
        "minimum_transport_improvement_factor": decimal(min_improvement[0]),
        "minimum_transport_improvement_factor_witness": min_improvement[1],
        "outputs_within_requested_scale": sum(bool(row["within_requested_scale"]) for row in output_rows),
    }
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_retained_decay_endpoint_majorant_gate",
        "status": "rigorous_retained_decay_endpoint_majorant_closes_finite_block20_nominal_local_scale",
        "theorem": {
            "generic": "2*pi*integral (B*r^2+C*r^3)*exp(-alpha*q*r)*phi(beta*r^2+gamma*r^3)/(1-exp(-alpha*r)) dr",
            "homotopy_factor": "phi(x)=(1-exp(-x))/x with phi(0)=1",
            "exceptional_b": "2*pi*C*integral (3*N*r^2+r^3)*exp(-alpha*delta*r-a_b*r^2) dr",
            "exceptional_c": "2*pi*C*integral r^3*exp(-alpha*eta*r-a_c*r^2) dr",
            "gap_boundary": "The exceptional bounds remain finite at delta=0 and eta=0.",
        },
        "quadrature": {
            "precision_ladder_decimal_digits": list(PRECISIONS),
            "origin_cutoff": decimal(ORIGIN),
            "far_cutoff": decimal(CUTOFF),
            "origin_bound": "2*pi*exp(alpha*epsilon)/alpha*(B*epsilon^2/2+C*epsilon^3/3)",
            "generic_tail": "phi<=1, geometric denominator frozen at R, exact exponential moments I2 and I3",
            "exceptional_tail": "linear decay dropped, exact Gaussian moments I2 and I3",
        },
        "rows": rows,
        "transported_outputs": output_rows,
        "aggregate": aggregate,
        "dependencies": {
            "prior_majorant": {"path": relative(PRIOR), "sha256": file_hash(PRIOR)},
            "weight_fixture": {"path": relative(WEIGHT_FIXTURE), "sha256": file_hash(WEIGHT_FIXTURE)},
            "weights": {"path": relative(WEIGHTS), "sha256": file_hash(WEIGHTS)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": (
            "Promote the retained-decay theorem, then keep coefficient-neighborhood transport, recurrence "
            "arithmetic, Legendre truncation, cross-block accumulation, and the outer Hardy remainder as "
            "separate budgets."
        ),
        "proof_boundary": (
            "Rigorous finite block-20 application of selector-gap-finite retained-decay endpoint majorants. "
            "Closing the 0.005 local endpoint-model scale is not a complete evaluator error theorem and does "
            "not control coefficient neighborhoods, other blocks, the outer Hardy representation, Lambda <= 0, "
            "PF-infinity, RH, or a prize-level conclusion."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print(
        "built retained-decay endpoint majorant: "
        f"max transported {aggregate['maximum_retained_decay_transported_majorant_upper']}, "
        f"{aggregate['outputs_within_requested_scale']}/15 below 0.005"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
