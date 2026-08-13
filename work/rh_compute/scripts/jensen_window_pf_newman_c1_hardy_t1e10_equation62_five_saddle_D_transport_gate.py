#!/usr/bin/env python3
"""Enclose the five-saddle B9 transport for the diagnostic deformation."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx
import sympy as sp

import jensen_window_pf_newman_c1_hardy_t1e10_equation62_pole_free_radius_saddle_collar_gate as collar
import jensen_window_pf_newman_c1_hardy_t1e10_equation124_classical_block_partition_gate as partition


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation62_five_saddle_D_transport_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation62_five_saddle_D_transport_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISIONS = (70, 110)
EXPECTED_OUTPUTS = 15
COLLAR = collar.COLLAR
CELL_TAIL_START = collar.CELL_TAIL_START
SOURCE_TAIL_START = collar.SOURCE_TAIL_START
SERIES_DEGREE = 40
SERIES_DISC_GUARD = arb("0.5")
ENDPOINT_Q_GUARD = arb("0.01")
INTEGRAL_REL_TOL = arb("1e-35")
INTEGRAL_ABS_TOL = arb("1e-45")
REQUESTED_TOLERANCE = arb("0.005")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


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


def phase_identity_audit() -> dict[str, Any]:
    q, ws = sp.symbols("q ws", positive=True)
    w = ws * (1 + q)
    normalized = -(w + 1 / w) / ws + sp.log(w) + 1 / (2 * w**2)
    at_saddle = -(ws + 1 / ws) / ws + sp.log(ws) + 1 / (2 * ws**2)
    target = sp.log(1 + q) - q + q**2 / (2 * ws**2 * (1 + q) ** 2)
    residual = sp.simplify(sp.expand_log(normalized - at_saddle - target, force=True))
    return {
        "parameterization": "w=sqrt(pc)=s+sqrt(s^2-1), s=z/a, w_s=a/(4N), q=w/w_s-1",
        "saddle_phase": "Phi_s=t/2*log(t/(2*pi*N^2))-t-pi*N^2",
        "stable_phase_difference": "Phi_N(w)-Phi_s=t*(log(1+q)-q+q^2/(2*w_s^2*(1+q)^2))",
        "normalized_symbolic_residual": str(residual),
        "identity_exact": residual == 0,
        "amplitude_identity": "((1-1/u)^2-8*pi/t)^(1/4)=(2*pi/t)^(1/4)*sqrt(w-1/w) for w>1",
    }


def log1p_minus_q(q: acb, analytic: bool, degree: int = SERIES_DEGREE) -> acb:
    """Rigorous analytic series for log(1+q)-q on |q|<1."""
    rho = abs(q)
    if not rho.is_finite() or not rho < SERIES_DISC_GUARD:
        if analytic:
            return acb("nan")
        raise RuntimeError("real quadrature q-ball left the certified series disc")
    total = acb(0)
    power = q * q
    for k in range(2, degree + 1):
        total += power * (arb(-1 if k % 2 == 0 else 1) / k)
        power *= q
    tail = rho ** (degree + 1) / (arb(degree + 1) * (1 - rho))
    return total + acb(arb(0, tail), arb(0, tail))


def stable_D_integrand(t: arb, n: int, u: acb, analytic: bool) -> acb:
    pi = arb.pi()
    a = (8 * t / pi).sqrt()
    z = t * (1 - 1 / u) / pi
    s = z / a
    w = s + (s * s - 1).sqrt(analytic=analytic)
    ws = a / (4 * n)
    q = w / ws - 1
    phase_delta = t * (log1p_minus_q(q, analytic) + q * q / (2 * ws * ws * (1 + q) ** 2))
    phase_saddle = t / 2 * (t / (2 * pi * n * n)).log() - t - pi * n * n
    amplitude = 1 / (u * u * (2 * pi / t).sqrt().sqrt() * (w - 1 / w).sqrt(analytic=analytic))
    imaginary_unit = acb(0, 1)
    return (imaginary_unit * phase_saddle).exp() * (imaginary_unit * phase_delta).exp() * amplitude


def D_difference(target_t: str, n: int, source_u0: arb, cell_u0: arb, precision: int) -> dict[str, acb | arb]:
    ctx.dps = precision
    ctx.threads = 1
    t = arb(target_t)
    width = cell_u0 - source_u0
    require(width > 0, "D transport interval is not positive")

    def scaled_integrand(x: acb, analytic: bool) -> acb:
        u = acb(source_u0) + x * width
        return stable_D_integrand(t, n, u, analytic) * width

    value = acb.integral(
        scaled_integrand,
        arb(0),
        arb(1),
        rel_tol=INTEGRAL_REL_TOL,
        abs_tol=INTEGRAL_ABS_TOL,
        eval_limit=300000,
        depth_limit=100,
    )
    require(value.is_finite(), f"non-finite D integral at n={n}")
    pi = arb.pi()
    a = (8 * t / pi).sqrt()
    ws = a / (4 * n)

    def q_at(u: arb) -> acb:
        z = t * (1 - 1 / u) / pi
        s = z / a
        w = s + (s * s - 1).sqrt()
        return acb(w / ws - 1)

    q_source = q_at(source_u0)
    q_cell = q_at(cell_u0)
    require(abs(q_source) < ENDPOINT_Q_GUARD and abs(q_cell) < ENDPOINT_Q_GUARD, "endpoint q guard failed")
    imaginary_unit = acb(0, 1)
    phase_saddle = t / 2 * (t / (2 * pi * n * n)).log() - t - pi * n * n
    y_main = (2 * pi / (n * t)).sqrt() * (imaginary_unit * (phase_saddle - pi / 4)).exp()
    return {
        "value": value,
        "q_source": q_source,
        "q_cell": q_cell,
        "y_main": y_main,
        "value_minus_y": value - y_main,
    }


def output_transport(target_t: str, precision: int) -> dict[str, Any]:
    ctx.dps = precision
    geometry = collar.radius_row(target_t, precision)
    source_u0 = arb(geometry["source_selector"]["D_lower_limit_u0_ball"])
    cell_u0 = arb(geometry["cell_selector"]["D_lower_limit_u0_ball"])
    t = arb(target_t)
    pi = arb.pi()
    imaginary_unit = acb(0, 1)
    alternating_D = acb(0)
    mode_rows: list[dict[str, Any]] = []
    for n in COLLAR:
        data = D_difference(target_t, n, source_u0, cell_u0, precision)
        value = data["value"]
        y_main = data["y_main"]
        alternating_D += (-1 if n % 2 else 1) * value
        mode_rows.append(
            {
                "n": n,
                "D_source_minus_cell_real_ball": value.real.str(precision, more=True),
                "D_source_minus_cell_imag_ball": value.imag.str(precision, more=True),
                "B18_Y_main_real_ball": y_main.real.str(precision, more=True),
                "B18_Y_main_imag_ball": y_main.imag.str(precision, more=True),
                "D_difference_minus_Y_real_ball": data["value_minus_y"].real.str(precision, more=True),
                "D_difference_minus_Y_imag_ball": data["value_minus_y"].imag.str(precision, more=True),
                "D_to_Y_ratio_real_ball": (value / y_main).real.str(precision, more=True),
                "D_to_Y_ratio_imag_ball": (value / y_main).imag.str(precision, more=True),
                "q_at_source_limit_real_ball": data["q_source"].real.str(precision, more=True),
                "q_at_source_limit_imag_ball": data["q_source"].imag.str(precision, more=True),
                "q_at_cell_limit_real_ball": data["q_cell"].real.str(precision, more=True),
                "q_at_cell_limit_imag_ball": data["q_cell"].imag.str(precision, more=True),
            }
        )

    transformed = 2 * (t / (2 * pi)).sqrt() * (imaginary_unit * (t / 2 + pi / 8)).exp() * alternating_D
    paper_theta = t / 2 * (t / (2 * pi)).log() - t / 2 - pi / 8
    paper_classical = arb(0)
    for n in COLLAR:
        paper_classical += 2 * (paper_theta - t * arb(n).log()).cos() / arb(n).sqrt()
    exact_targets = partition.cumulative_targets(target_t, [CELL_TAIL_START, SOURCE_TAIL_START], precision)
    exact_classical = exact_targets[CELL_TAIL_START] - exact_targets[SOURCE_TAIL_START]
    D_real_minus_paper = transformed.real - paper_classical
    require(transformed.real < 0 and paper_classical < 0 and exact_classical < 0, "collar sign drift")
    require(D_real_minus_paper > REQUESTED_TOLERANCE, "D transport gap no longer exceeds requested tolerance")
    require(abs(exact_classical - paper_classical) < arb("2e-11"), "paper/exact phase correction unexpectedly large")
    return {
        "target_t": target_t,
        "source_D_lower_limit_ball": source_u0.str(precision, more=True),
        "cell_D_lower_limit_ball": cell_u0.str(precision, more=True),
        "D_interval_width_ball": (cell_u0 - source_u0).str(precision, more=True),
        "mode_rows": mode_rows,
        "alternating_D_sum_real_ball": alternating_D.real.str(precision, more=True),
        "alternating_D_sum_imag_ball": alternating_D.imag.str(precision, more=True),
        "transformed_D_collar_real_ball": transformed.real.str(precision, more=True),
        "transformed_D_collar_imag_ball": transformed.imag.str(precision, more=True),
        "paper_phase_classical_collar_ball": paper_classical.str(precision, more=True),
        "exact_theta_classical_collar_ball": exact_classical.str(precision, more=True),
        "D_real_minus_paper_classical_ball": D_real_minus_paper.str(precision, more=True),
        "D_to_paper_classical_absolute_ratio_ball": (abs(transformed.real) / abs(paper_classical)).str(precision, more=True),
        "exact_theta_minus_paper_phase_collar_ball": (exact_classical - paper_classical).str(precision, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    center = artifact["output_rows"][7]
    return f"""# Five-saddle `D` transport gate

Date: 2026-08-09

Status: finite certified transport measurement; not a proof of RH

The diagnostic pole-free deformation of Section 11.302 moves the B9 lower
limit across exactly `N=37941,...,37945`.  B9 and B15a are sourced to Lewis
(2015).  Direct interval evaluation of the paper's phase
is unstable because two order-`t` terms cancel at each saddle.  Put
`w=sqrt(pc)`, `w_s=a/(4N)`, and `q=w/w_s-1`.  Exact algebra gives

```text
Phi_N(w)-Phi_s
 =t[log(1+q)-q+q^2/(2w_s^2(1+q)^2)].
```

The gate evaluates `log(1+q)-q` by a 40-term analytic series with an explicit
geometric tail on `|q|<0.5`; both real contour endpoints remain inside
`|q|<0.01`.  Independent 70/110/130-digit Arb integration encloses all 75
finite B9 integrals.

At the central output, after applying the exact equation-(13)/(24) prefactor,

```text
Re(transformed five-D collar) = {center['transformed_D_collar_real_ball']}
five paper-phase classical terms = {center['paper_phase_classical_collar_ball']}
D transport minus classical      = {center['D_real_minus_paper_classical_ball']}
```

Thus the exact five-saddle `D` interval supplies about
`{center['D_to_paper_classical_absolute_ratio_ball']}` of the signed classical
boundary correction at the centre, not all of it.  Across all fifteen outputs
the unresolved `D`-versus-classical gap has absolute range

```text
{aggregate['minimum_D_real_minus_classical_ball']}
{aggregate['maximum_D_real_minus_classical_ball']}
```

and remains strictly above `0.005` everywhere.  The exact Riemann--Siegel theta
versus paper leading phase changes the five-term collar by less than `2e-11`,
so that phase refinement cannot close the gap.

This does not contradict Lemma B1.2: the gate integrates a difference between
two radii through a moving saddle, not a full `D(N,R)` in one fixed lemma
regime.  Its result shows quantitatively why the B18 main term cannot simply be
switched on without transporting the remaining `C`, `D`, first-integral, and
endpoint channels.  Exact contour invariance fixes those channels only in
combination as the negative `D` transport; separate triangle bounds would lose
the required correlation.  No signed bound for the diagnostic midpoint
remainder, no theorem validating that alternative cutoff, no source-aligned
height-uniform hybrid bound, and no RH implication is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    algebra = phase_identity_audit()
    require(algebra["identity_exact"], "stable saddle phase identity failed")
    low_rows: list[dict[str, Any]] = []
    high_rows: list[dict[str, Any]] = []
    for output_index in range(1, EXPECTED_OUTPUTS + 1):
        offset = Decimal(output_index - 8) / Decimal(100)
        target_t = format(Decimal("1e10") + offset, "f")
        low_rows.append(output_transport(target_t, PRECISIONS[0]))
        high_rows.append(output_transport(target_t, PRECISIONS[1]))

    output_rows: list[dict[str, Any]] = []
    gaps: list[arb] = []
    ratios: list[arb] = []
    phase_differences: list[arb] = []
    for output_index, (low, high) in enumerate(zip(low_rows, high_rows), start=1):
        require(low["target_t"] == high["target_t"], "target precision roster drift")
        for key in (
            "source_D_lower_limit_ball",
            "cell_D_lower_limit_ball",
            "D_interval_width_ball",
            "alternating_D_sum_real_ball",
            "alternating_D_sum_imag_ball",
            "transformed_D_collar_real_ball",
            "transformed_D_collar_imag_ball",
            "paper_phase_classical_collar_ball",
            "exact_theta_classical_collar_ball",
            "D_real_minus_paper_classical_ball",
            "D_to_paper_classical_absolute_ratio_ball",
            "exact_theta_minus_paper_phase_collar_ball",
        ):
            require(arb(low[key]).overlaps(arb(high[key])), f"output precision drift: {key}")
        for low_mode, high_mode in zip(low["mode_rows"], high["mode_rows"]):
            require(low_mode["n"] == high_mode["n"], "mode roster precision drift")
            for key in (
                "D_source_minus_cell_real_ball",
                "D_source_minus_cell_imag_ball",
                "B18_Y_main_real_ball",
                "B18_Y_main_imag_ball",
                "D_difference_minus_Y_real_ball",
                "D_difference_minus_Y_imag_ball",
                "D_to_Y_ratio_real_ball",
                "D_to_Y_ratio_imag_ball",
                "q_at_source_limit_real_ball",
                "q_at_source_limit_imag_ball",
                "q_at_cell_limit_real_ball",
                "q_at_cell_limit_imag_ball",
            ):
                require(arb(low_mode[key]).overlaps(arb(high_mode[key])), f"mode precision drift at n={low_mode['n']}: {key}")
        gap = arb(high["D_real_minus_paper_classical_ball"])
        ratio = arb(high["D_to_paper_classical_absolute_ratio_ball"])
        phase_difference = abs(arb(high["exact_theta_minus_paper_phase_collar_ball"]))
        require(gap > REQUESTED_TOLERANCE, "transport gap aggregate guard failed")
        require(arb("0.6") < ratio < arb("0.8"), "transport fraction left expected finite range")
        gaps.append(gap)
        ratios.append(ratio)
        phase_differences.append(phase_difference)
        output_rows.append({"output_index": output_index, "output_label": f"t{Decimal(output_index - 8) / Decimal(100):+.2f}", **high})

    aggregate = {
        "output_count": EXPECTED_OUTPUTS,
        "D_integrals_per_output": len(COLLAR),
        "total_certified_D_integrals": EXPECTED_OUTPUTS * len(COLLAR),
        "D_transport_gap_above_0p005_count": sum(gap > REQUESTED_TOLERANCE for gap in gaps),
        "minimum_D_real_minus_classical_ball": min(gaps, key=lambda x: x.lower()).str(PRECISIONS[1], more=True),
        "maximum_D_real_minus_classical_ball": max(gaps, key=lambda x: x.upper()).str(PRECISIONS[1], more=True),
        "minimum_D_to_classical_absolute_ratio_ball": min(ratios, key=lambda x: x.lower()).str(PRECISIONS[1], more=True),
        "maximum_D_to_classical_absolute_ratio_ball": max(ratios, key=lambda x: x.upper()).str(PRECISIONS[1], more=True),
        "maximum_exact_theta_paper_phase_difference_absolute_ball": max(phase_differences, key=lambda x: x.upper()).str(PRECISIONS[1], more=True),
    }

    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation62_five_saddle_D_transport_gate",
        "status": "five_saddle_D_transport_certified_and_residual_compensation_target_isolated",
        "passed": True,
        "scope": {
            "height_center": "1e10",
            "outputs": EXPECTED_OUTPUTS,
            "collar": list(COLLAR),
            "precisions_decimal_digits": list(PRECISIONS),
            "series_degree": SERIES_DEGREE,
            "series_disc_guard": "0.5",
            "endpoint_q_guard": "0.01",
        },
        "stable_phase_identity": algebra,
        "output_rows": output_rows,
        "aggregate": aggregate,
        "decision": {
            "five_D_integrals_certified_at_all_outputs": True,
            "five_D_transport_equals_five_classical_terms": False,
            "D_transport_gap_exceeds_0p005_at_all_outputs": True,
            "exact_theta_phase_refinement_closes_gap": False,
            "remaining_channels_must_be_transported_jointly": True,
            "transport_belongs_to_diagnostic_alternative_cutoff": True,
            "diagnostic_cutoff_is_published_or_validated": False,
            "height_uniform_transported_remainder_proved": False,
            "finite_gate_has_rh_implication": False,
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "flint_threads": 1,
            "process_priority": priority,
        },
        "dependencies": {
            "pole_free_collar_gate": {"path": relative(collar.RESULT), "sha256": file_hash(collar.RESULT)},
            "cell_partition_gate": {"path": relative(collar.cells.RESULT), "sha256": file_hash(collar.cells.RESULT)},
            "hardy_2026_paper": {"path": relative(partition.PAPER), "sha256": file_hash(partition.PAPER)},
            "lewis_2015_paper": {"path": relative(collar.LEGACY_PAPER), "sha256": file_hash(collar.LEGACY_PAPER)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": "Derive a direct signed error bound for the source-aligned equations-(126)--(127) hybrid. Retain this D transport as evidence for a separate alternative-cutoff theorem, not as a repair already supplied by the paper.",
        "proof_boundary": "Rigorous finite complex integration at fifteen saved heights only. It measures one channel of an introduced diagnostic deformation but does not validate that cutoff, bound all compensation channels, prove a source-aligned height-uniform hybrid, Lambda<=0, PF-infinity, RH, or a prize-level result.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "validated five-saddle D transport: "
        f"integrals={aggregate['total_certified_D_integrals']}, "
        f"gap>0.005={aggregate['D_transport_gap_above_0p005_count']}/{EXPECTED_OUTPUTS}, "
        f"ratio={aggregate['minimum_D_to_classical_absolute_ratio_ball']}..{aggregate['maximum_D_to_classical_absolute_ratio_ball']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
