#!/usr/bin/env python3
"""Certify W4 contour recombination and complete endpoint Hurwitz majorants."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_full_gate as full
import jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_pilot_gate as pilot
import jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_pilot_gate as rays
import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_block20_shared_contour_exceptional_majorant_gate as exceptional
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import acb, arb, ctx


RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_w4_contour_endpoint_hurwitz_majorant_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_w4_contour_endpoint_hurwitz_majorant_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = (
    REPO_ROOT
    / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_w4_contour_endpoint_hurwitz_majorant_gate.py"
)
PRECISIONS = pilot.PRECISIONS
HALF_RAY_CUTOFF = Fraction(64)
HALF_RAY_BREAKS = (
    Fraction(0),
    Fraction(1, 100),
    Fraction(1, 10),
    Fraction(1, 2),
    Fraction(1),
    Fraction(4),
    Fraction(16),
    HALF_RAY_CUTOFF,
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


def lower(value: arb) -> Fraction:
    return cells.bound_fraction(value.lower())


def upper_abs(value: acb) -> Fraction:
    return upper(abs(value))


def lower_abs(value: acb) -> Fraction:
    return lower(abs(value))


def zero_half_ray(data: dict[str, Any], endpoint: int, dps: int) -> tuple[acb, Fraction, str]:
    """Lower half-ray integral of exp(-2*pi*i*F) with an explicit tail."""
    ctx.dps = dps
    ctx.threads = 1
    parent = data["levels"][1]
    phi1_fraction, phi2_fraction, phi3_fraction = (
        cells.binary128_fraction(value) for value in parent["coefficients_hex"]
    )
    phi1, phi2, phi3 = (cells.fraction_ball(value) for value in (phi1_fraction, phi2_fraction, phi3_fraction))
    root_two = arb(2).sqrt()
    root_six = arb(6).sqrt()
    if phi3_fraction >= 0:
        angle = "-pi/6"
        direction = acb(arb(3).sqrt() / 2, -cells.fraction_ball(Fraction(1, 2)))
        sine_one = cells.fraction_ball(Fraction(1, 2))
        sine_two = arb(3).sqrt() / 2
        cubic_sine = arb(1)
    else:
        angle = "-5*pi/12"
        direction = acb((root_six - root_two) / 4, -(root_six + root_two) / 4)
        sine_one = (root_six + root_two) / 4
        sine_two = cells.fraction_ball(Fraction(1, 2))
        cubic_sine = root_two / 2
    endpoint_ball = arb(endpoint)
    pi = arb.pi()

    def phase(z: acb) -> acb:
        return phi1 * z + phi2 * z**2 + phi3 * z**3

    def integrand(radius: acb, _analytic: bool) -> acb:
        u = direction * radius
        return (-2 * pi * I * phase(endpoint_ball + u)).exp() * direction

    tolerance = arb(f"1e-{dps - 25}")
    total = acb(0)
    for lower_break, upper_break in zip(HALF_RAY_BREAKS[:-1], HALF_RAY_BREAKS[1:]):
        total += acb.integral(
            integrand,
            acb(cells.fraction_ball(lower_break)),
            acb(cells.fraction_ball(upper_break)),
            abs_tol=tolerance,
            rel_tol=tolerance,
            deg_limit=64,
            eval_limit=200000,
            depth_limit=64,
        )

    slope = phi1 + 2 * phi2 * endpoint_ball + 3 * phi3 * endpoint_ball**2
    quadratic = phi2 + 3 * phi3 * endpoint_ball
    require(slope > 0 and quadratic > 0, "W4 half-ray decay coefficient not positive")
    far = cells.fraction_ball(HALF_RAY_CUTOFF)
    decay_rate = 2 * pi * (
        slope * sine_one + quadratic * sine_two * far + abs(phi3) * cubic_sine * far**2
    )
    tail_radius = (-decay_rate * far).exp() / decay_rate
    tail_upper = upper(tail_radius)
    return total + rays.error_box(tail_upper), tail_upper, angle


def generic_linear_majorant(quadratic: arb, cubic_abs: arb, gap_start: arb) -> arb:
    """Sum exact-minus-linear mode bounds from the positive gap_start onward."""
    require(quadratic > 0 and gap_start > 0, "generic Hurwitz majorant domain failure")
    pi = arb.pi()
    zeta_three = arb(3).zeta(gap_start)
    zeta_four = arb(4).zeta(gap_start)
    return 4 * quadratic * zeta_three / pi**2 + 12 * cubic_abs * zeta_four / pi**3


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# W4 contour recombination and complete endpoint Hurwitz majorants

Date: 2026-08-07

Status: rigorous symbolic selector-gap-uniform local endpoint theorem with an all-374 finite-roster application; not a proof of the recurrence budget or RH

## W4 is an endpoint half-ray difference

For `Phi1>0`, deform the finite zero-mode integral to the two lower
half-rays on which `exp(-2*pi*i*F)` decays.  Cauchy's theorem gives

```text
zero_mode=H_exact(0,0)-H_exact(N,0).                  (1)
```

The quadratic half-ray in W4 is `H_quad(0,0)`, while

```text
H_linear(N,0)=-i*exp(-2*pi*i*F(N))/(2*pi*F'(N)).     (2)
```

Hence

```text
zero_mode-W4
 =[H_exact(0,0)-H_quad(0,0)]
  -[H_exact(N,0)-H_linear(N,0)].                     (3)
```

The explicit W4 endpoint term is therefore the missing linear `n=0` ray,
not an independent error.  Arb verifies (1)--(3) at 70 and 110 digits on all
{aggregate['w4_contour_identity_count']} positive-`Phi1` calls.  The largest
identity gap is below `{aggregate['maximum_w4_contour_identity_gap_upper']}`.

## Complete endpoint majorant

For an exact-minus-linear mode with starting gap `q>0`, local quadratic
coefficient `B>0`, and `C=|Phi3|`, a common-ray homotopy and the lower bound
`sin(theta)>=1/2` give

```text
sum_(k>=0)|D_(q+k)|
 <= 4*B*zeta(3,q)/pi^2 + 12*C*zeta(4,q)/pi^3.        (4)
```

The upper family uses (4) at `q=m-Phi1` at zero and `q=m+1-xi` at `N`,
plus the W2 bound of Section 11.271.  The lower W3 family starts at
`q=2+Phi1` and `q=1+xi`; the W4 family starts at `q=1+Phi1` and `q=xi`,
plus the common quadratic zero-mode bound.  Thus every Hurwitz argument is
bounded away from the exceptional selector gap that was removed before
estimation.

Arb certifies the resulting upper, lower, and complete local bounds on all
374 calls:

```text
maximum upper-family bound   <= {aggregate['maximum_upper_family_majorant_upper']}
maximum lower-family bound   <= {aggregate['maximum_lower_family_majorant_upper']}
maximum complete local bound <= {aggregate['maximum_complete_majorant_upper']}
worst complete bound/actual ratio <= {aggregate['maximum_complete_majorant_to_correction_ratio_upper']}
```

Pi in (1)--(4) comes from the Fourier phase `exp(2*pi*i*n*u)` and the exact
Gaussian/exponential moments; it is not fitted.  The theorem is uniform in
the selector gaps but conservative.  It does not yet prove that its bound
fits the complete recurrence budget, control coefficient transport, Legendre
tails, cross-block scaling, the outer Hardy remainder, `Lambda<=0`,
PF-infinity, RH, or a prize-level conclusion.
"""


def main() -> int:
    rays.set_low_priority()
    recursive = native_q.load_recursive_chains()
    source = json.loads(full.RESULT.read_text(encoding="utf-8"))
    source_rows = {int(row["chain"]): row for row in source["rows"]}
    decomposition = json.loads(pilot.SOURCE.read_text(encoding="utf-8"))
    decomposition_rows = {int(row["chain"]): row for row in decomposition["rows"]}
    require(set(recursive) == set(source_rows) == set(decomposition_rows), "W4/Hurwitz roster drift")

    rows: list[dict[str, Any]] = []
    max_w4_gap = (Fraction(0), 0)
    max_half_tail = (Fraction(0), 0)
    min_upper_slack: tuple[Fraction, int] | None = None
    min_lower_slack: tuple[Fraction, int] | None = None
    min_complete_slack: tuple[Fraction, int] | None = None
    max_upper_bound = (Fraction(0), 0)
    max_lower_bound = (Fraction(0), 0)
    max_complete_bound = (Fraction(0), 0)
    max_upper_ratio = (Fraction(0), 0)
    max_lower_ratio = (Fraction(0), 0)
    max_complete_ratio = (Fraction(0), 0)
    w4_count = 0
    for chain in sorted(recursive):
        data = recursive[chain]
        parent = data["levels"][1]
        phi1_fraction, phi2_fraction, phi3_fraction = (
            cells.binary128_fraction(value) for value in parent["coefficients_hex"]
        )
        ctx.dps = 110
        phi1, phi2, phi3 = (
            cells.fraction_ball(value) for value in (phi1_fraction, phi2_fraction, phi3_fraction)
        )
        cubic_abs = abs(phi3)
        length = int(parent["length"])
        length_ball = arb(length)
        xi = phi1 + 2 * phi2 * length_ball + 3 * phi3 * length_ball**2
        endpoint_quadratic = phi2 + 3 * phi3 * length_ball
        mode_b = int(data["levels"][2]["length"]) + 1
        b_exceptional = exceptional.b_majorant(phi2_fraction, phi3_fraction, length)[0]
        upper_bound_ball = b_exceptional + generic_linear_majorant(
            phi2, cubic_abs, arb(mode_b) - phi1
        ) + generic_linear_majorant(endpoint_quadratic, cubic_abs, arb(mode_b + 1) - xi)

        c_exceptional = exceptional.c_majorant(phi2_fraction, phi3_fraction)[0]
        if phi1_fraction < 0:
            lower_bound_ball = c_exceptional + generic_linear_majorant(
                phi2, cubic_abs, arb(2) + phi1
            ) + generic_linear_majorant(endpoint_quadratic, cubic_abs, arb(1) + xi)
            w4_gap = Fraction(0)
            half_tail = Fraction(0)
            w4_angles: list[str] = []
        else:
            require(phi1_fraction > 0, "W4 selector on zero")
            w4_count += 1
            low_h0, _, _ = zero_half_ray(data, 0, PRECISIONS[0])
            low_hn, _, _ = zero_half_ray(data, length, PRECISIONS[0])
            high_h0, tail_h0, angle_h0 = zero_half_ray(data, 0, PRECISIONS[1])
            high_hn, tail_hn, angle_hn = zero_half_ray(data, length, PRECISIONS[1])
            require(low_h0.overlaps(high_h0) and low_hn.overlaps(high_hn), f"chain {chain} W4 half-ray precision nonoverlap")
            pi = arb.pi()
            phase_length = phi1 * length_ball + phi2 * length_ball**2 + phi3 * length_ball**3
            endpoint_minus = (-2 * pi * I * phase_length).exp()
            quadratic_zero = pilot.paper_fresnel(phi2, phi1).conjugate()
            endpoint_linear_zero = -I * endpoint_minus / (2 * pi * xi)
            zero_mode = native_q.complex_from_record(decomposition_rows[chain]["zero_mode"])
            paper_w4 = native_q.complex_from_record(decomposition_rows[chain]["paper_w3_or_w4"])
            zero_contour_gap = (high_h0 - high_hn) - zero_mode
            w4_model_gap = (quadratic_zero - endpoint_linear_zero) - paper_w4
            w4_recombination_gap = (
                (high_h0 - quadratic_zero) - (high_hn - endpoint_linear_zero)
            ) - (zero_mode - paper_w4)
            for name, value in (
                ("zero contour", zero_contour_gap),
                ("W4 model", w4_model_gap),
                ("W4 recombination", w4_recombination_gap),
            ):
                require(value.contains(0), f"chain {chain} {name} identity excludes zero")
            w4_gap = max(
                upper_abs(zero_contour_gap), upper_abs(w4_model_gap), upper_abs(w4_recombination_gap)
            )
            half_tail = max(tail_h0, tail_hn)
            w4_angles = [angle_h0, angle_hn]
            max_w4_gap = max(max_w4_gap, (w4_gap, chain))
            max_half_tail = max(max_half_tail, (half_tail, chain))
            lower_bound_ball = c_exceptional + generic_linear_majorant(
                phi2, cubic_abs, arb(1) + phi1
            ) + generic_linear_majorant(endpoint_quadratic, cubic_abs, xi)

        upper_actual = native_q.complex_from_record(source_rows[chain]["upper_recombined"])
        lower_actual = native_q.complex_from_record(source_rows[chain]["lower_recombined"])
        correction = native_q.complex_from_record(decomposition_rows[chain]["exact_correction"])
        upper_bound_lower, upper_bound_upper = lower(upper_bound_ball), upper(upper_bound_ball)
        lower_bound_lower, lower_bound_upper = lower(lower_bound_ball), upper(lower_bound_ball)
        complete_bound_ball = upper_bound_ball + lower_bound_ball
        complete_bound_lower, complete_bound_upper = lower(complete_bound_ball), upper(complete_bound_ball)
        upper_actual_lower, upper_actual_upper = lower_abs(upper_actual), upper_abs(upper_actual)
        lower_actual_lower, lower_actual_upper = lower_abs(lower_actual), upper_abs(lower_actual)
        correction_lower, correction_upper = lower_abs(correction), upper_abs(correction)
        require(upper_actual_lower > 0 and lower_actual_lower > 0 and correction_lower > 0, f"chain {chain} endpoint value contains zero")
        require(upper_bound_lower > upper_actual_upper, f"chain {chain} upper Hurwitz majorant failure")
        require(lower_bound_lower > lower_actual_upper, f"chain {chain} lower Hurwitz majorant failure")
        require(complete_bound_lower > correction_upper, f"chain {chain} complete Hurwitz majorant failure")
        upper_slack = upper_bound_lower - upper_actual_upper
        lower_slack = lower_bound_lower - lower_actual_upper
        complete_slack = complete_bound_lower - correction_upper
        upper_ratio = upper_bound_upper / upper_actual_lower
        lower_ratio = lower_bound_upper / lower_actual_lower
        complete_ratio = complete_bound_upper / correction_lower
        min_upper_slack = (upper_slack, chain) if min_upper_slack is None else min(min_upper_slack, (upper_slack, chain))
        min_lower_slack = (lower_slack, chain) if min_lower_slack is None else min(min_lower_slack, (lower_slack, chain))
        min_complete_slack = (complete_slack, chain) if min_complete_slack is None else min(min_complete_slack, (complete_slack, chain))
        max_upper_bound = max(max_upper_bound, (upper_bound_upper, chain))
        max_lower_bound = max(max_lower_bound, (lower_bound_upper, chain))
        max_complete_bound = max(max_complete_bound, (complete_bound_upper, chain))
        max_upper_ratio = max(max_upper_ratio, (upper_ratio, chain))
        max_lower_ratio = max(max_lower_ratio, (lower_ratio, chain))
        max_complete_ratio = max(max_complete_ratio, (complete_ratio, chain))
        rows.append(
            {
                "chain": chain,
                "sum_index": int(data["header"]["sum_index"]),
                "branch": int(data["header"]["branch"]),
                "selector": source_rows[chain]["selector"],
                "w4_half_ray_angles": w4_angles,
                "w4_contour_identity_gap_upper": decimal(w4_gap),
                "w4_half_ray_tail_upper": decimal(half_tail),
                "upper_family_majorant_lower": decimal(upper_bound_lower),
                "upper_family_majorant_upper": decimal(upper_bound_upper),
                "upper_family_actual_upper": decimal(upper_actual_upper),
                "lower_family_majorant_lower": decimal(lower_bound_lower),
                "lower_family_majorant_upper": decimal(lower_bound_upper),
                "lower_family_actual_upper": decimal(lower_actual_upper),
                "complete_majorant_lower": decimal(complete_bound_lower),
                "complete_majorant_upper": decimal(complete_bound_upper),
                "correction_actual_upper": decimal(correction_upper),
                "upper_majorant_slack_lower": decimal(upper_slack),
                "lower_majorant_slack_lower": decimal(lower_slack),
                "complete_majorant_slack_lower": decimal(complete_slack),
                "upper_majorant_to_actual_ratio_upper": decimal(upper_ratio),
                "lower_majorant_to_actual_ratio_upper": decimal(lower_ratio),
                "complete_majorant_to_correction_ratio_upper": decimal(complete_ratio),
            }
        )

    require(min_upper_slack is not None and min_lower_slack is not None and min_complete_slack is not None, "Hurwitz aggregate missing")
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_w4_contour_endpoint_hurwitz_majorant_gate",
        "status": "rigorous_W4_contour_and_complete_selector_gap_uniform_endpoint_majorants_validated_all_374",
        "theorem": {
            "W4_contour": "zero_mode-W4=(H_exact(0,0)-H_quad(0,0))-(H_exact(N,0)-H_linear(N,0))",
            "generic_tail": "sum_(k>=0)|D_(q+k)|<=4*B*zeta(3,q)/pi^2+12*|Phi3|*zeta(4,q)/pi^3",
            "upper_arguments": "q=m-Phi1 at 0 and q=m+1-xi at N, plus the W2 exceptional bound",
            "lower_W3_arguments": "q=2+Phi1 at 0 and q=1+xi at N, plus the W3 exceptional bound",
            "lower_W4_arguments": "q=1+Phi1 at 0 and q=xi at N, plus the quadratic zero-mode bound",
        },
        "rows": rows,
        "aggregate": {
            "call_count": len(rows),
            "w4_contour_identity_count": w4_count,
            "upper_majorant_count": sum(Fraction(row["upper_majorant_slack_lower"]) > 0 for row in rows),
            "lower_majorant_count": sum(Fraction(row["lower_majorant_slack_lower"]) > 0 for row in rows),
            "complete_majorant_count": sum(Fraction(row["complete_majorant_slack_lower"]) > 0 for row in rows),
            "maximum_w4_contour_identity_gap_upper": decimal(max_w4_gap[0]),
            "maximum_w4_contour_identity_gap_witness": max_w4_gap[1],
            "maximum_w4_half_ray_tail_upper": decimal(max_half_tail[0]),
            "maximum_w4_half_ray_tail_witness": max_half_tail[1],
            "minimum_upper_majorant_slack_lower": decimal(min_upper_slack[0]),
            "minimum_upper_majorant_slack_witness": min_upper_slack[1],
            "minimum_lower_majorant_slack_lower": decimal(min_lower_slack[0]),
            "minimum_lower_majorant_slack_witness": min_lower_slack[1],
            "minimum_complete_majorant_slack_lower": decimal(min_complete_slack[0]),
            "minimum_complete_majorant_slack_witness": min_complete_slack[1],
            "maximum_upper_family_majorant_upper": decimal(max_upper_bound[0]),
            "maximum_upper_family_majorant_witness": max_upper_bound[1],
            "maximum_lower_family_majorant_upper": decimal(max_lower_bound[0]),
            "maximum_lower_family_majorant_witness": max_lower_bound[1],
            "maximum_complete_majorant_upper": decimal(max_complete_bound[0]),
            "maximum_complete_majorant_witness": max_complete_bound[1],
            "maximum_upper_majorant_to_actual_ratio_upper": decimal(max_upper_ratio[0]),
            "maximum_upper_majorant_to_actual_ratio_witness": max_upper_ratio[1],
            "maximum_lower_majorant_to_actual_ratio_upper": decimal(max_lower_ratio[0]),
            "maximum_lower_majorant_to_actual_ratio_witness": max_lower_ratio[1],
            "maximum_complete_majorant_to_correction_ratio_upper": decimal(max_complete_ratio[0]),
            "maximum_complete_majorant_to_correction_ratio_witness": max_complete_ratio[1],
        },
        "dependencies": {
            "exceptional_full": {"path": relative(full.RESULT), "sha256": file_hash(full.RESULT)},
            "decomposition": {"path": relative(pilot.SOURCE), "sha256": file_hash(pilot.SOURCE)},
            "exceptional_majorant": {"path": relative(exceptional.RESULT), "sha256": file_hash(exceptional.RESULT)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "proof_boundary": (
            "Exact W4 contour identity and symbolic selector-gap-uniform local endpoint bounds, with rigorous "
            "application to all 374 saved calls. No proof that the conservative bound fits the complete recurrence, "
            "coefficient transport, Legendre and outer Hardy remainders, Lambda<=0, PF-infinity, RH, or a prize-level conclusion."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print("built W4/Hurwitz endpoint majorant gate: 209/209 W4 contours and 374/374 complete local bounds")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
