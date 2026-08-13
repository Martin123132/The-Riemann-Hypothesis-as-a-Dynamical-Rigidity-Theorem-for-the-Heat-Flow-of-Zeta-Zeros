#!/usr/bin/env python3
"""Certify singularity-free exceptional-mode recombinations on four branches."""

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

import jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_pilot_gate as decomposition
import jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_pilot_gate as rays
import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import acb, arb, ctx


SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_full_gate.json"
)
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_pilot_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_pilot_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = (
    REPO_ROOT
    / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_pilot_gate.py"
)
WITNESS_CHAINS = decomposition.WITNESS_CHAINS
PRECISIONS = decomposition.PRECISIONS
MODE_BREAKS = (Fraction(0), Fraction(1, 100), Fraction(1, 10), Fraction(1, 2), Fraction(1), Fraction(4), Fraction(16))
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


def lower_abs(value: acb) -> Fraction:
    return cells.bound_fraction(abs(value).lower())


def paper_fresnel(phi2: arb, argument: arb) -> acb:
    """Paper's quadratic half-ray integral with global coefficient Phi2."""
    pi = arb.pi()
    root_two = arb(2).sqrt()
    epi_minus = acb(1 / root_two, -1 / root_two)
    xr = 2 * phi2
    z = epi_minus * (pi / xr).sqrt() * argument
    return I * epi_minus * acb(0, -pi * argument**2 / xr).exp() * (1 - z.erf()) / (2 * xr.sqrt())


def exceptional_mode_integral(
    data: dict[str, Any],
    sector_row: dict[str, Any],
    endpoint: int,
    family: str,
    mode: int,
    dps: int,
    ray_cutoff: Fraction,
) -> tuple[acb, Fraction]:
    """Direct Arb integral of one exact logarithmic-ray Fourier mode."""
    ctx.dps = dps
    ctx.threads = 1
    parent = data["levels"][1]
    phi1, phi2, phi3 = (rays.exact_ball(value) for value in parent["coefficients_hex"])
    contract = sector_row[f"{family}_contract"]
    ray_direction, _ = rays.direction(contract["angle"])
    phase_sign = -1 if family == "b" else 1
    endpoint_ball = arb(endpoint)
    pi = arb.pi()

    def phase(z: acb) -> acb:
        return phi1 * z + phi2 * z**2 + phi3 * z**3

    def phase_derivative(z: acb) -> acb:
        return phi1 + 2 * phi2 * z + 3 * phi3 * z**2

    def integrand(radius: acb, _analytic: bool) -> acb:
        u = ray_direction * radius
        z = endpoint_ball + u
        exponent = phase_sign * 2 * pi * I * phase(z) + 2 * pi * I * mode * u
        return phase_derivative(z) * exponent.exp() * ray_direction / mode

    tolerance = arb(f"1e-{dps - 25}")
    total = acb(0)
    breaks = MODE_BREAKS + (ray_cutoff,)
    for lower_break, upper_break in zip(breaks[:-1], breaks[1:]):
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

    far = cells.fraction_ball(ray_cutoff)
    linear_decay = cells.fraction_ball(rays.fraction_from_record(contract["linear_decay"]))
    quadratic_decay = cells.fraction_ball(rays.dyadic(contract["quadratic_decay"]["lower_dyadic"]))
    cubic_decay = cells.fraction_ball(rays.fraction_from_record(contract["cubic_decay"]))
    decay_rate = 2 * pi * (linear_decay + quadratic_decay * far + cubic_decay * far**2)
    endpoint_derivative = phi1 + 2 * phi2 * endpoint_ball + 3 * phi3 * endpoint_ball**2
    endpoint_second = 2 * phi2 + 6 * phi3 * endpoint_ball
    p0 = abs(endpoint_derivative)
    p1 = abs(endpoint_second)
    p2 = 3 * abs(phi3)
    inverse = 1 / decay_rate
    polynomial_tail = (-decay_rate * far).exp() * (
        p0 * inverse
        + p1 * (far * inverse + inverse**2)
        + p2 * (far**2 * inverse + 2 * far * inverse**2 + 2 * inverse**3)
    )
    tail_radius = polynomial_tail / mode
    tail_upper = upper(tail_radius)
    return total + rays.error_box(tail_upper), tail_upper


def recombine(
    data: dict[str, Any],
    sector_row: dict[str, Any],
    source_row: dict[str, Any],
    dps: int,
) -> dict[str, Any]:
    ctx.dps = dps
    ctx.threads = 1
    parent = data["levels"][1]
    phi1_fraction, phi2_fraction, phi3_fraction = (
        cells.binary128_fraction(value) for value in parent["coefficients_hex"]
    )
    phi1, phi2, phi3 = (cells.fraction_ball(value) for value in (phi1_fraction, phi2_fraction, phi3_fraction))
    length = int(parent["length"])
    transformed_length = int(data["levels"][2]["length"])
    mode_b = transformed_length + 1
    length_ball = arb(length)
    xi = phi1 + 2 * phi2 * length_ball + 3 * phi3 * length_ball**2
    delta = arb(mode_b) - xi
    pi = arb.pi()
    phase_length = phi1 * length_ball + phi2 * length_ball**2 + phi3 * length_ball**3
    endpoint_minus = (-2 * pi * I * phase_length).exp()
    fresnel_b = paper_fresnel(phi2, delta)
    b_quadratic_phase_integral = -fresnel_b.conjugate()
    b_linear_mode = endpoint_minus * (I / (2 * pi * delta) - I / (2 * pi * mode_b))
    b_quadratic_mode = endpoint_minus * (b_quadratic_phase_integral - I / (2 * pi * mode_b))
    paper_w2 = native_q.complex_from_record(source_row["paper_w2"])
    w2_model_gap = (b_linear_mode - b_quadratic_mode) - paper_w2

    ray_cutoff = Fraction(int(source_row["ray_cutoff"]))
    b_exact_mode, b_tail = exceptional_mode_integral(
        data, sector_row, length, "b", mode_b, dps, ray_cutoff
    )
    b_remainder = native_q.complex_from_record(source_row["b_remainder"])
    upper_generic = b_remainder + (b_exact_mode - b_linear_mode)
    upper_exceptional = b_exact_mode - b_quadratic_mode
    upper_recombined = upper_generic - upper_exceptional
    upper_existing = b_remainder - paper_w2
    upper_recombination_gap = upper_recombined - upper_existing

    c_remainder = native_q.complex_from_record(source_row["c_remainder"])
    zero_mode = native_q.complex_from_record(source_row["zero_mode"])
    paper_w34 = native_q.complex_from_record(source_row["paper_w3_or_w4"])
    if phi1_fraction < 0:
        selector = "W3"
        eta = 1 + phi1
        c_exact_mode, c_tail = exceptional_mode_integral(data, sector_row, 0, "c", 1, dps, ray_cutoff)
        c_linear_mode = I / (2 * pi) - I / (2 * pi * eta)
        c_quadratic_mode = I / (2 * pi) - paper_fresnel(phi2, eta)
        w34_model_gap = (c_linear_mode - c_quadratic_mode).conjugate() - paper_w34
        lower_generic = c_remainder + (c_exact_mode - c_linear_mode).conjugate()
        lower_exceptional = (c_exact_mode - c_quadratic_mode).conjugate()
        lower_recombined = lower_generic - lower_exceptional
        lower_components = (lower_generic, -lower_exceptional)
    else:
        require(phi1_fraction > 0, "Phi1 selector on zero")
        selector = "W4"
        c_tail = Fraction(0)
        quadratic_zero_mode = paper_fresnel(phi2, phi1).conjugate()
        endpoint_boundary = I * endpoint_minus / (2 * pi * xi)
        w34_model_gap = quadratic_zero_mode + endpoint_boundary - paper_w34
        lower_zero_remainder = zero_mode - quadratic_zero_mode
        lower_c_boundary_remainder = c_remainder - endpoint_boundary
        lower_recombined = lower_zero_remainder + lower_c_boundary_remainder
        lower_components = (lower_zero_remainder, lower_c_boundary_remainder)

    lower_existing = c_remainder + zero_mode - paper_w34
    lower_recombination_gap = lower_recombined - lower_existing
    correction = native_q.complex_from_record(source_row["exact_correction"])
    total_recombined = upper_recombined + lower_recombined
    total_gap = total_recombined - correction
    component_triangle = upper_abs(upper_generic) + upper_abs(upper_exceptional)
    component_triangle += sum(upper_abs(value) for value in lower_components)
    return {
        "selector": selector,
        "b_exact_mode": b_exact_mode,
        "b_linear_mode": b_linear_mode,
        "b_quadratic_mode": b_quadratic_mode,
        "w2_model_gap": w2_model_gap,
        "upper_generic": upper_generic,
        "upper_exceptional": upper_exceptional,
        "upper_recombined": upper_recombined,
        "upper_recombination_gap": upper_recombination_gap,
        "w34_model_gap": w34_model_gap,
        "lower_recombined": lower_recombined,
        "lower_recombination_gap": lower_recombination_gap,
        "total_recombined": total_recombined,
        "total_gap": total_gap,
        "correction": correction,
        "component_triangle": component_triangle,
        "b_tail_upper": b_tail,
        "c_tail_upper": c_tail,
    }


def build_note(artifact: dict[str, Any]) -> str:
    rows = "\n".join(
        f"chain {row['chain']:>3} ({row['selector']}): model gap <= {row['maximum_model_identity_gap_upper']}, "
        f"recombination gap <= {row['total_recombination_gap_upper']}"
        for row in artifact["rows"]
    )
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 exceptional-mode recombination pilot

Date: 2026-08-07

Status: rigorous four-branch finite-input identity and quadrature certificate; not a height-uniform bound and not a proof of RH

## Why the reciprocal poles disappear

The factor `2*pi` comes from the Fourier phase `exp(2*pi*i*n*u)` in the
Poisson endpoint rays.  The paper's exceptional quadratic models retain its
global coefficient `Phi2`; they do not use the exact local Taylor coefficient
`Phi2+3*Phi3*x`.  Let `M` denote a single endpoint-ray mode and let `Z_quad`
denote the positive-`Phi1` quadratic zero-mode half-ray.  Direct evaluation of
the paper formulas gives

```text
W2 = M_b,lin(N)-M_b,quad(N),
W3 = conjugate(M_c,lin(0)-M_c,quad(0)),
W4 = Z_quad(0)+i*exp(-2*pi*i*F(N))/(2*pi*F'(N)).       (1)
```

Therefore the upper endpoint pair is exactly

```text
R_b-W2
 = G_b-[M_b,exact(N)-M_b,quad(N)],                    (2)
```

where `G_b` contains only full-minus-linear generic modes.  In the W3 branch,

```text
R_c-W3
 = G_c-conjugate[M_c,exact(0)-M_c,quad(0)].           (3)
```

Equations (2)--(3) contain neither `1/delta`, for
`delta=ceil(F'(N))-F'(N)`, nor `1/eta`, for `eta=1+Phi1`.  In the W4 branch the
lower pair instead splits into the finite differences
`Z_exact-Z_quad` and `R_c-i*exp(-2*pi*i*F(N))/(2*pi*F'(N))`.

## Arb certificate

The exceptional exact modes are integrated directly on the already certified
cubic rays, with explicit polynomial-times-exponential tail balls.  The 70-
and 110-digit enclosures overlap, and (1)--(3) plus the complete correction
identity enclose zero on every pilot branch:

```text
{rows}
```

Aggregate bounds:

```text
maximum model identity gap       <= {aggregate['maximum_model_identity_gap_upper']}
maximum recombination gap        <= {aggregate['maximum_total_recombination_gap_upper']}
maximum exceptional-mode tail    <= {aggregate['maximum_exceptional_mode_tail_upper']}
maximum four-piece triangle/correction ratio <= {aggregate['maximum_component_triangle_to_correction_ratio_upper']}
```

This removes the apparent small-gap singularities at the identity level.  It
does not yet supply a height-uniform majorant for the finite differences, an
all-374 replay, recursive accumulation, outer Hardy control, `Lambda<=0`,
PF-infinity, RH, or a prize-level conclusion.
"""


def main() -> int:
    recursive = native_q.load_recursive_chains()
    sector = json.loads(rays.SECTOR.read_text(encoding="utf-8"))
    sector_rows = {int(row["chain"]): row for row in sector["rows"]}
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    source_rows = {int(row["chain"]): row for row in source["rows"]}
    rows: list[dict[str, Any]] = []
    max_model = Fraction(0)
    max_recombination = Fraction(0)
    max_tail = Fraction(0)
    max_ratio = Fraction(0)
    for chain in WITNESS_CHAINS:
        low = recombine(recursive[chain], sector_rows[chain], source_rows[chain], PRECISIONS[0])
        high = recombine(recursive[chain], sector_rows[chain], source_rows[chain], PRECISIONS[1])
        for name in (
            "b_exact_mode",
            "b_linear_mode",
            "b_quadratic_mode",
            "w2_model_gap",
            "upper_generic",
            "upper_exceptional",
            "upper_recombined",
            "upper_recombination_gap",
            "w34_model_gap",
            "lower_recombined",
            "lower_recombination_gap",
            "total_recombined",
            "total_gap",
        ):
            require(low[name].overlaps(high[name]), f"chain {chain} {name} precision nonoverlap")
        for name in ("w2_model_gap", "w34_model_gap", "upper_recombination_gap", "lower_recombination_gap", "total_gap"):
            require(high[name].contains(0), f"chain {chain} {name} excludes zero")
        model_gap = max(upper_abs(high["w2_model_gap"]), upper_abs(high["w34_model_gap"]))
        recombination_gap = upper_abs(high["total_gap"])
        tail = max(high["b_tail_upper"], high["c_tail_upper"])
        correction_lower = lower_abs(high["correction"])
        require(correction_lower > 0, f"chain {chain} correction contains zero")
        ratio = high["component_triangle"] / correction_lower
        max_model = max(max_model, model_gap)
        max_recombination = max(max_recombination, recombination_gap)
        max_tail = max(max_tail, tail)
        max_ratio = max(max_ratio, ratio)
        rows.append(
            {
                "chain": chain,
                "selector": high["selector"],
                "b_exact_mode": point_balls.acb_record(high["b_exact_mode"], 45),
                "b_exact_minus_quadratic": point_balls.acb_record(high["upper_exceptional"], 45),
                "upper_recombined": point_balls.acb_record(high["upper_recombined"], 45),
                "lower_recombined": point_balls.acb_record(high["lower_recombined"], 45),
                "maximum_model_identity_gap_upper": decimal(model_gap),
                "total_recombination_gap_upper": decimal(recombination_gap),
                "exceptional_mode_tail_upper": decimal(tail),
                "component_triangle_to_correction_ratio_upper": decimal(ratio),
            }
        )

    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_pilot_gate",
        "status": "rigorous_four_branch_exceptional_mode_recombination_validated",
        "precisions_decimal_digits": list(PRECISIONS),
        "witness_chains": list(WITNESS_CHAINS),
        "identities": {
            "W2": "W2=M_b_linear(N)-M_b_quadratic_global_Phi2(N)",
            "W3": "W3=conjugate(M_c_linear(0)-M_c_quadratic_global_Phi2(0))",
            "W4": "W4=Z_quadratic_global_Phi2(0)+i*exp(-2*pi*i*F(N))/(2*pi*F'(N))",
            "upper": "R_b-W2=G_b-(M_b_exact(N)-M_b_quadratic_global_Phi2(N))",
            "lower_W3": "R_c-W3=G_c-conjugate(M_c_exact(0)-M_c_quadratic_global_Phi2(0))",
        },
        "rows": rows,
        "aggregate": {
            "branch_count": len(rows),
            "maximum_model_identity_gap_upper": decimal(max_model),
            "maximum_total_recombination_gap_upper": decimal(max_recombination),
            "maximum_exceptional_mode_tail_upper": decimal(max_tail),
            "maximum_component_triangle_to_correction_ratio_upper": decimal(max_ratio),
        },
        "dependencies": {
            "source": {"path": relative(SOURCE), "sha256": file_hash(SOURCE)},
            "sector": {"path": relative(rays.SECTOR), "sha256": file_hash(rays.SECTOR)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "proof_boundary": (
            "Rigorous four-branch finite-input identity and mode-quadrature certificate only. No all-call or "
            "height-uniform finite-difference majorant, recursive accumulation, outer Hardy control, Lambda<=0, "
            "PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print("built exceptional-mode recombination pilot: 4/4 branches, singular reciprocal terms removed algebraically")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
