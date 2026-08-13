#!/usr/bin/env python3
"""Rigorous four-branch Arb pilot for the logarithmic endpoint rays."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import acb, arb, ctx


SECTOR = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_cubic_sector_contour_repair_gate.json"
)
LOG_REDUCTION = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_reduction_gate.json"
)
EQ69 = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_exact_eq69_saddle_family_gate.json"
)
W2_W5 = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_w2_w5_correspondence_gate.json"
)
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_pilot_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_pilot_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = (
    REPO_ROOT
    / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_pilot_gate.py"
)
WITNESS_CHAINS = (1, 2, 3, 36)
PRECISIONS = (70, 110)
ORIGIN_CUTOFF = Fraction(1, 10**16)
RAY_CUTOFF = Fraction(64)
CENTRAL_BREAKS = (
    Fraction(1, 10**16),
    Fraction(1, 10**12),
    Fraction(1, 10**8),
    Fraction(1, 10**4),
    Fraction(1, 100),
    Fraction(1, 10),
    Fraction(1, 2),
    Fraction(1),
    Fraction(4),
    Fraction(16),
    Fraction(64),
)
I = acb(0, 1)


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


def upper_abs(value: acb) -> Fraction:
    return upper(abs(value))


def fraction_from_record(record: dict[str, Any]) -> Fraction:
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def dyadic(pair: list[int]) -> Fraction:
    mantissa, exponent = (int(value) for value in pair)
    if exponent >= 0:
        return Fraction(mantissa * 2**exponent)
    return Fraction(mantissa, 2 ** (-exponent))


def error_box(radius: Fraction) -> acb:
    text = decimal(radius)
    component = arb(f"[0 +/- {text}]")
    return acb(component, component)


def harmonic(number: int) -> Fraction:
    return sum((Fraction(1, index) for index in range(1, number + 1)), Fraction(0))


def direction(angle: str) -> tuple[acb, Fraction]:
    if angle == "pi/2":
        return I, Fraction(1)
    root_three_over_two = arb(3).sqrt() / 2
    if angle == "pi/6":
        return acb(root_three_over_two, cells.fraction_ball(Fraction(1, 2))), Fraction(1, 2)
    require(angle == "5*pi/6", f"unknown ray angle: {angle}")
    return acb(-root_three_over_two, cells.fraction_ball(Fraction(1, 2))), Fraction(1, 2)


def exact_ball(value: str) -> arb:
    return cells.fraction_ball(cells.binary128_fraction(value))


def kernel(exponential: acb, start: int, hypergeometric: bool) -> acb:
    if hypergeometric:
        return exponential**start * exponential.hypgeom_2f1(1, start, start + 1) / start
    value = -(-exponential).log1p()
    for index in range(1, start):
        value -= exponential**index / index
    return value


def ray_integral(
    data: dict[str, Any],
    sector_row: dict[str, Any],
    endpoint: int,
    family: str,
    dps: int,
    ray_cutoff: Fraction = RAY_CUTOFF,
    central_breaks: tuple[Fraction, ...] = CENTRAL_BREAKS,
) -> dict[str, Any]:
    ctx.dps = dps
    ctx.threads = 1
    parent = data["levels"][1]
    phi1, phi2, phi3 = (exact_ball(value) for value in parent["coefficients_hex"])
    contract = sector_row[f"{family}_contract"]
    angle = contract["angle"]
    ray_direction, sine = direction(angle)
    phase_sign = -1 if family == "b" else 1
    start = int(data["levels"][2]["length"]) + 1 if family == "b" else 1
    endpoint_ball = arb(endpoint)
    pi = arb.pi()

    def phase(z: acb) -> acb:
        return phi1 * z + phi2 * z**2 + phi3 * z**3

    def phase_derivative(z: acb) -> acb:
        return phi1 + 2 * phi2 * z + 3 * phi3 * z**2

    def make_integrand(hypergeometric: bool):
        def integrand(radius: acb, _analytic: bool) -> acb:
            u = ray_direction * radius
            z = endpoint_ball + u
            exponential = (2 * pi * I * u).exp()
            summed_tail = kernel(exponential, start, hypergeometric)
            phase_factor = (phase_sign * 2 * pi * I * phase(z)).exp()
            return phase_derivative(z) * phase_factor * summed_tail * ray_direction

        return integrand

    tolerance = arb(f"1e-{dps - 25}")
    central = acb(0)
    require(central_breaks[0] == ORIGIN_CUTOFF, "ray central origin cutoff drift")
    require(central_breaks[-1] == ray_cutoff, "ray central far cutoff drift")
    for lower_break, upper_break in zip(central_breaks[:-1], central_breaks[1:]):
        use_hypergeometric = lower_break >= Fraction(1, 2)
        central += acb.integral(
            make_integrand(use_hypergeometric),
            acb(cells.fraction_ball(lower_break)),
            acb(cells.fraction_ball(upper_break)),
            abs_tol=tolerance,
            rel_tol=tolerance,
            deg_limit=64,
            eval_limit=200000,
            depth_limit=64,
        )

    epsilon = cells.fraction_ball(ORIGIN_CUTOFF)
    far = cells.fraction_ball(ray_cutoff)
    angular_decay = 2 * pi * cells.fraction_ball(sine)
    endpoint_derivative = phi1 + 2 * phi2 * endpoint_ball + 3 * phi3 * endpoint_ball**2
    endpoint_second = 2 * phi2 + 6 * phi3 * endpoint_ball
    amplitude_origin = (
        abs(endpoint_derivative)
        + abs(endpoint_second) * epsilon
        + 3 * abs(phi3) * epsilon**2
    )
    kernel_origin_integral = (angular_decay * epsilon).exp() * (
        epsilon * (1 - (angular_decay * epsilon).log())
        + angular_decay * epsilon**2 / 2
    )
    origin_radius = abs(ray_direction) * amplitude_origin * kernel_origin_integral

    linear_decay = cells.fraction_ball(fraction_from_record(contract["linear_decay"]))
    quadratic_decay = cells.fraction_ball(dyadic(contract["quadratic_decay"]["lower_dyadic"]))
    cubic_decay = cells.fraction_ball(fraction_from_record(contract["cubic_decay"]))
    decay_rate = 2 * pi * (
        linear_decay + quadratic_decay * far + cubic_decay * far**2
    )
    p0 = abs(endpoint_derivative)
    p1 = abs(endpoint_second)
    p2 = 3 * abs(phi3)
    inverse = 1 / decay_rate
    polynomial_tail = (-decay_rate * far).exp() * (
        p0 * inverse
        + p1 * (far * inverse + inverse**2)
        + p2 * (far**2 * inverse + 2 * far * inverse**2 + 2 * inverse**3)
    )
    geometric_denominator = start * (1 - (-angular_decay * far).exp())
    tail_radius = abs(ray_direction) * polynomial_tail / geometric_denominator

    origin_upper = upper(origin_radius)
    tail_upper = upper(tail_radius)
    total = central + error_box(origin_upper + tail_upper)
    return {
        "family": family,
        "endpoint": endpoint,
        "angle": angle,
        "start": start,
        "central": central,
        "origin_upper": origin_upper,
        "tail_upper": tail_upper,
        "total": total,
    }


def finite_zero_mode(data: dict[str, Any], dps: int) -> acb:
    ctx.dps = dps
    parent = data["levels"][1]
    phi1, phi2, phi3 = (exact_ball(value) for value in parent["coefficients_hex"])
    pi = arb.pi()

    def integrand(y: acb, _analytic: bool) -> acb:
        phase = phi1 * y + phi2 * y**2 + phi3 * y**3
        return (-2 * pi * I * phase).exp()

    tolerance = arb(f"1e-{dps - 25}")
    total = acb(0)
    for lower_break in range(0, 104, 13):
        total += acb.integral(
            integrand,
            acb(lower_break),
            acb(lower_break + 13),
            abs_tol=tolerance,
            rel_tol=tolerance,
            deg_limit=64,
            eval_limit=200000,
            depth_limit=64,
        )
    return total


def evaluate(
    data: dict[str, Any],
    sector_row: dict[str, Any],
    eq69_row: dict[str, Any],
    w2_row: dict[str, Any],
    dps: int,
    ray_cutoff: Fraction = RAY_CUTOFF,
    central_breaks: tuple[Fraction, ...] = CENTRAL_BREAKS,
) -> dict[str, Any]:
    ctx.dps = dps
    parent = data["levels"][1]
    phi1, phi2, phi3 = (exact_ball(value) for value in parent["coefficients_hex"])
    length = int(parent["length"])
    transformed_length = int(data["levels"][2]["length"])
    b_zero = ray_integral(data, sector_row, 0, "b", dps, ray_cutoff, central_breaks)
    b_length = ray_integral(data, sector_row, length, "b", dps, ray_cutoff, central_breaks)
    c_zero = ray_integral(data, sector_row, 0, "c", dps, ray_cutoff, central_breaks)
    c_length = ray_integral(data, sector_row, length, "c", dps, ray_cutoff, central_breaks)
    b_source = b_zero["total"] - b_length["total"]
    c_source = (c_length["total"] - c_zero["total"]).conjugate()

    pi = arb.pi()
    endpoint_phase = phi1 * length + phi2 * length**2 + phi3 * length**3
    endpoint = (-2 * pi * I * endpoint_phase).exp()
    endpoint_half = (1 + endpoint) / 2
    harmonic_ball = cells.fraction_ball(harmonic(transformed_length))
    a_endpoint = I * (endpoint - 1) * harmonic_ball / (2 * pi)
    zero_mode = finite_zero_mode(data, dps) if cells.binary128_fraction(parent["coefficients_hex"][0]) > 0 else acb(0)
    exact_nonsaddle = endpoint_half + a_endpoint + b_source + c_source + zero_mode

    parent_sum = point_balls.direct_sum(parent, -2 * pi)
    saddle_sum = native_q.complex_from_record(eq69_row["paper_phase_integral_sum"])
    independent_target = parent_sum - saddle_sum
    formula_target_gap = exact_nonsaddle - independent_target
    paper_nonsaddle = native_q.complex_from_record(w2_row["paper_nonsaddle_aggregate"])
    paper_residual = native_q.complex_from_record(w2_row["paper_residual"])
    exact_minus_paper = exact_nonsaddle - paper_nonsaddle
    correction_residual_gap = exact_minus_paper - paper_residual
    rays = (b_zero, b_length, c_zero, c_length)
    return {
        "rays": rays,
        "b_source": b_source,
        "c_source": c_source,
        "endpoint_half": endpoint_half,
        "a_endpoint": a_endpoint,
        "zero_mode": zero_mode,
        "exact_nonsaddle": exact_nonsaddle,
        "independent_target": independent_target,
        "formula_target_gap": formula_target_gap,
        "paper_nonsaddle": paper_nonsaddle,
        "paper_residual": paper_residual,
        "exact_minus_paper": exact_minus_paper,
        "correction_residual_gap": correction_residual_gap,
        "maximum_origin_upper": max(ray["origin_upper"] for ray in rays),
        "maximum_tail_upper": max(ray["tail_upper"] for ray in rays),
    }


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    rows = "\n".join(
        f"chain {row['chain']:>3}  Phi1 {row['phi1_sign']:<8} Phi3 {row['phi3_sign']:<8} "
        f"|formula-target| <= {row['formula_target_gap_magnitude_upper']}"
        for row in artifact["rows"]
    )
    return f"""# Hardy block-20 logarithmic endpoint-ray Arb pilot gate

Date: 2026-08-07

Status: rigorous four-branch finite-input Arb pilot; not a proof of all-call or height-uniform W2--W4 control or RH

## Certified decomposition

This pilot evaluates the exact logarithmic-ray identity from Formal Core
Section 11.264 on one witness from each `(sign(Phi1),sign(Phi3))` branch.
For each endpoint ray it integrates on

```text
epsilon <= r <= R,   epsilon={artifact['quadrature']['origin_cutoff']}, R={artifact['quadrature']['ray_cutoff']}.
```

Near the origin the exact absolute majorant is

```text
|F'(x+r*e^(i theta))| * exp(a*epsilon)
* [epsilon*(1-log(a*epsilon))+a*epsilon^2/2],
a=2*pi*sin(theta).                                      (1)
```

It follows from
`sum_(k>=0) exp(-a*k*r)/(m+k) <= exp(a*r)*[-log(1-exp(-a*r))]`
and `1-exp(-a*r)>=a*r*exp(-a*r)`.  Thus the logarithmic endpoint is
enclosed rather than sampled.

For `r>=R`, the certified phase bound
`Im phase >= ell*r+q*r^2+c*r^3` is replaced by the linear lower bound
`(ell+q*R+c*R^2)r`.  The resulting polynomial-times-exponential integral
is evaluated explicitly, including the geometric tail of the summed
`1/n` kernel.

On the compact pieces the kernel is evaluated as a principal logarithm near
the origin and as

```text
K_m(w)=w^m*2F1(1,m;m+1;w)/m                            (2)
```

for small `|w|`, avoiding cancellation of the first `m-1` terms.

## Four-branch result

The 70- and 110-digit constructions overlap for every ray and complete
nonsaddle value.  Both independent identities contain zero:

```text
Q_exact^- - (P^- - I69^-) = 0,
(Q_exact^- - Q_paper^-) - R_paper = 0.
```

```text
{rows}
```

Aggregate bounds:

```text
maximum omitted-origin radius       <= {aggregate['maximum_origin_radius_upper']}
maximum far-tail radius             <= {aggregate['maximum_tail_radius_upper']}
maximum |formula-target gap|        <= {aggregate['maximum_formula_target_gap_magnitude_upper']}
maximum |correction-residual gap|    <= {aggregate['maximum_correction_residual_gap_magnitude_upper']}
```

This independently confirms the source orientation, endpoint harmonic term,
positive-`Phi1` zero mode, and all four ray signs with rigorous balls.

## Pi provenance and boundary

Pi is the mathematical constant inherited from `exp(2*pi*i*phase)`.  It
enters the exact kernel, phase decay, and endpoint harmonic coefficient; no
fitted normalization is used.

This is a rigorous pilot on four fixed low-height calls.  It does not yet
enclose the other 370 calls, derive a symbolic height-uniform W2--W4 error,
control recursive accumulation or the outer Hardy remainder, or prove
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
"""


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def main() -> int:
    for path in (SECTOR, LOG_REDUCTION, EQ69, W2_W5, CHECKER):
        require(path.is_file(), f"missing log-ray Arb pilot dependency: {path}")
    priority = set_low_priority()
    sector = json.loads(SECTOR.read_text(encoding="utf-8"))
    reduction = json.loads(LOG_REDUCTION.read_text(encoding="utf-8"))
    eq69 = json.loads(EQ69.read_text(encoding="utf-8"))
    w2 = json.loads(W2_W5.read_text(encoding="utf-8"))
    require(
        reduction["status"] == "exact_logarithmic_endpoint_ray_reduction_with_nonrigorous_orientation_diagnostics",
        "log-ray Arb pilot reduction gate not admitted",
    )
    sector_rows = {int(row["chain"]): row for row in sector["rows"]}
    eq69_rows = {int(row["chain"]): row for row in eq69["rows"]}
    w2_rows = {int(row["chain"]): row for row in w2["rows"]}
    recursive = native_q.load_recursive_chains()

    rows: list[dict[str, Any]] = []
    maximum_origin = Fraction(0)
    maximum_tail = Fraction(0)
    maximum_formula_gap = Fraction(0)
    maximum_correction_gap = Fraction(0)
    for chain in WITNESS_CHAINS:
        low = evaluate(recursive[chain], sector_rows[chain], eq69_rows[chain], w2_rows[chain], PRECISIONS[0])
        high = evaluate(recursive[chain], sector_rows[chain], eq69_rows[chain], w2_rows[chain], PRECISIONS[1])
        for index, (low_ray, high_ray) in enumerate(zip(low["rays"], high["rays"])):
            require(low_ray["total"].overlaps(high_ray["total"]), f"chain {chain} ray {index} precision nonoverlap")
        for name in (
            "b_source",
            "c_source",
            "endpoint_half",
            "a_endpoint",
            "zero_mode",
            "exact_nonsaddle",
            "independent_target",
            "formula_target_gap",
            "exact_minus_paper",
            "correction_residual_gap",
        ):
            require(low[name].overlaps(high[name]), f"chain {chain} {name} precision nonoverlap")
        require(high["formula_target_gap"].contains(0), f"chain {chain} formula/target identity excludes zero")
        require(high["correction_residual_gap"].contains(0), f"chain {chain} correction/residual identity excludes zero")
        formula_gap = upper_abs(high["formula_target_gap"])
        correction_gap = upper_abs(high["correction_residual_gap"])
        maximum_origin = max(maximum_origin, high["maximum_origin_upper"])
        maximum_tail = max(maximum_tail, high["maximum_tail_upper"])
        maximum_formula_gap = max(maximum_formula_gap, formula_gap)
        maximum_correction_gap = max(maximum_correction_gap, correction_gap)
        parent = recursive[chain]["levels"][1]
        phi1 = cells.binary128_fraction(parent["coefficients_hex"][0])
        phi3 = cells.binary128_fraction(parent["coefficients_hex"][2])
        rows.append(
            {
                "chain": chain,
                "phi1_sign": "positive" if phi1 > 0 else "negative",
                "phi3_sign": "positive" if phi3 > 0 else "negative",
                "b_angle": sector_rows[chain]["b_contract"]["angle"],
                "c_angle": sector_rows[chain]["c_contract"]["angle"],
                "precision_overlap": True,
                "origin_radius_upper": decimal(high["maximum_origin_upper"]),
                "tail_radius_upper": decimal(high["maximum_tail_upper"]),
                "exact_nonsaddle": point_balls.acb_record(high["exact_nonsaddle"], 50),
                "independent_target": point_balls.acb_record(high["independent_target"], 50),
                "formula_target_gap": point_balls.acb_record(high["formula_target_gap"], 30),
                "formula_target_gap_magnitude_upper": decimal(formula_gap),
                "exact_minus_paper": point_balls.acb_record(high["exact_minus_paper"], 50),
                "paper_residual": point_balls.acb_record(high["paper_residual"], 50),
                "correction_residual_gap": point_balls.acb_record(high["correction_residual_gap"], 30),
                "correction_residual_gap_magnitude_upper": decimal(correction_gap),
                "formula_target_identity_contains_zero": True,
                "correction_residual_identity_contains_zero": True,
            }
        )

    require(len({(row["phi1_sign"], row["phi3_sign"]) for row in rows}) == 4, "log-ray Arb pilot sign coverage drift")
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_pilot_gate",
        "status": "rigorous_four_branch_logarithmic_endpoint_ray_arb_pilot_validated",
        "scope": {
            "block": 20,
            "witness_chains": list(WITNESS_CHAINS),
            "precision_ladder_decimal_digits": list(PRECISIONS),
            "precision_overlap_count": len(rows),
            "worker_count": 1,
            "process_priority": priority,
        },
        "quadrature": {
            "origin_cutoff": decimal(ORIGIN_CUTOFF),
            "ray_cutoff": decimal(RAY_CUTOFF),
            "central_breaks": [decimal(value) for value in CENTRAL_BREAKS],
            "near_kernel": "principal-log identity (11.264.1)",
            "far_kernel": "K_m(w)=w^m*2F1(1,m;m+1;w)/m",
            "origin_bound": "exp(a*eps)*M(eps)*[eps*(1-log(a*eps))+a*eps^2/2]",
            "tail_bound": "linearize ell*r+q*r^2+c*r^3 below by (ell+q*R+c*R^2)r and integrate the quadratic amplitude exactly",
        },
        "aggregate": {
            "witness_count": len(rows),
            "maximum_origin_radius_upper": decimal(maximum_origin),
            "maximum_tail_radius_upper": decimal(maximum_tail),
            "maximum_formula_target_gap_magnitude_upper": decimal(maximum_formula_gap),
            "maximum_correction_residual_gap_magnitude_upper": decimal(maximum_correction_gap),
            "formula_target_identity_count": len(rows),
            "correction_residual_identity_count": len(rows),
        },
        "rows": rows,
        "sources": {
            "sector_gate": {"path": relative(SECTOR), "sha256": file_hash(SECTOR)},
            "log_reduction_gate": {"path": relative(LOG_REDUCTION), "sha256": file_hash(LOG_REDUCTION)},
            "equation_69_gate": {"path": relative(EQ69), "sha256": file_hash(EQ69)},
            "w2_w5_gate": {"path": relative(W2_W5), "sha256": file_hash(W2_W5)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": (
            "Replace the fixed pilot roster by deterministic resumable evaluation of all 374 calls, selecting each "
            "finite ray cutoff from its exact decay contract and preserving two-precision overlap."
        ),
        "proof_boundary": (
            "Rigorous Arb quadrature, origin enclosure, and far-tail enclosure for four fixed low-height calls only. "
            "This does not certify the other 370 calls, prove a height-uniform W2--W4 remainder or recurrence, control "
            "the outer Hardy representation, or prove Lambda<=0, PF-infinity, RH, or a prize-level conclusion."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print("built logarithmic endpoint-ray Arb pilot: 4/4 sign branches, both exact identities enclose zero")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
