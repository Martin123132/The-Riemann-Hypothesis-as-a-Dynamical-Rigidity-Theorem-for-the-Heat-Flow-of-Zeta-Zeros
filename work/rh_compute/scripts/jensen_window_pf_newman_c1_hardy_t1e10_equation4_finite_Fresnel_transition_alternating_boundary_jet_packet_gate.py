#!/usr/bin/env python3
"""Certify the joined transition packet by an alternating boundary-jet collapse."""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "work" / "rh_compute" / "vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_transition_alternating_boundary_jet_packet_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"

HEIGHT = 10_000_000_000
SOURCE_ENDPOINT = 159_577
SOURCE_COUNT = 2_481_423
LATER_LABEL_COUNT = SOURCE_COUNT - 1
TRANSITION_START = 39_853
TRANSITION_END = 39_936
CELL_LEFT = arb("39852.5")
CELL_RIGHT = arb("39936.5")
JET_ORDER = 12
DERIVATIVE_SLABS = 128
INTEGRAL_PANELS = 32
SERIES_TERMS = 30
LOG_SERIES_DISC_GUARD = arb("0.2")

GROUPED_A = ROOT / "work" / "rh_compute" / "results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_"
    "transition_grouped_complex_A_endpoint_lift_gate.json"
)
COMMON_CELL = ROOT / "work" / "rh_compute" / "results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_cell_"
    "common_carrier_subtraction_reduction_gate.json"
)
FRESNEL_CONNECTION = ROOT / "work" / "rh_compute" / "results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_branch_"
    "Fresnel_Dirichlet_connection_and_mode_roster_gate.json"
)
EXACT_H = ROOT / "work" / "rh_compute" / "results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_exact_H_"
    "continuation_defect_target_gate.json"
)
GAMMA_BULK = ROOT / "work" / "rh_compute" / "results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_global_morse_"
    "gamma_bulk_gate.json"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


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


def complex_record(value: acb, digits: int = 70) -> dict[str, str]:
    return {
        "real_ball": value.real.str(digits, more=True),
        "imag_ball": value.imag.str(digits, more=True),
        "absolute_ball": abs(value).str(digits, more=True),
    }


def complex_from_record(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def add_complex_error(value: acb, radius: arb) -> acb:
    return acb(arb(value.real, radius), arb(value.imag, radius))


def real_interval(left: arb, right: arb) -> arb:
    require(left.upper() <= right.lower(), "invalid real interval")
    midpoint = (left + right) / 2
    radius = (right - left) / 2
    return arb(midpoint, radius)


def centered_parameters(t: arb, endpoint: int) -> dict[str, Any]:
    pi = arb.pi()
    nu = (t / (2 * pi)).sqrt()
    delta = pi * (endpoint - 4 * nu)
    phase_nu = pi * endpoint * nu - pi * nu**2 - t * nu.log()
    F_nu = nu ** arb("-0.5") * (acb(0, 1) * phase_nu).exp()
    return {"pi": pi, "t": t, "endpoint": endpoint, "nu": nu, "delta": delta, "F_nu": F_nu}


def centered_first_label_value(y: acb, p: dict[str, Any], terms: int, analytic: bool) -> acb:
    z = y - p["nu"]
    w = z / p["nu"]
    rho = abs(w)
    if not rho.is_finite() or not rho < LOG_SERIES_DISC_GUARD:
        if analytic:
            return acb("nan")
        raise RuntimeError("centered logarithm series escaped its disk")
    radius = rho.upper()
    remainder = w * 0
    power = w**3
    for degree in range(3, terms + 1):
        remainder += power * ((-1) ** (degree + 1)) / degree
        power *= w
    tail = radius ** (terms + 1) / ((terms + 1) * (1 - radius))
    remainder += acb(arb(0, tail), arb(0, tail))
    root = (1 + w).sqrt(analytic=analytic)
    relative_phase = p["delta"] * z - p["t"] * remainder
    return p["F_nu"] / root * (acb(0, 1) * relative_phase).exp()


def first_label_integral(
    p: dict[str, Any], left: arb, right: arb, panels: int, terms: int, tolerance: str
) -> tuple[acb, list[dict[str, str]]]:
    total = acb(0)
    rows: list[dict[str, str]] = []

    def integrand(y: acb, analytic: bool) -> acb:
        return centered_first_label_value(y, p, terms, analytic)

    width = right - left
    for index in range(panels):
        panel_left = left + width * index / panels
        panel_right = left + width * (index + 1) / panels
        value = acb.integral(
            integrand,
            panel_left,
            panel_right,
            abs_tol=arb(tolerance),
            rel_tol=arb(tolerance),
            eval_limit=300_000,
            depth_limit=40,
        )
        total += value
        rows.append(
            {
                "left": panel_left.str(35, more=True),
                "right": panel_right.str(35, more=True),
                "absolute_ball": abs(value).str(35, more=True),
            }
        )
    return total, rows


def log_derivatives(y: arb, p: dict[str, Any], order: int) -> list[acb]:
    w = (y - p["nu"]) / p["nu"]
    q = p["delta"] - (p["t"] / p["nu"]) * w**2 / (1 + w)
    q_prime = -(p["t"] / p["nu"]**2) * w * (2 + w) / (1 + w) ** 2
    derivatives: list[acb] = []
    for n in range(1, order + 1):
        real_part = ((-1) ** n) * math.factorial(n - 1) / (2 * y**n)
        if n == 1:
            imag_part = q
        elif n == 2:
            imag_part = q_prime
        else:
            imag_part = ((-1) ** n) * math.factorial(n - 1) * p["t"] / y**n
        derivatives.append(acb(real_part, imag_part))
    return derivatives


def bell_jet(log_jets: list[acb]) -> list[acb]:
    complete = [acb(1)]
    for n in range(len(log_jets)):
        value = acb(0)
        for k in range(n + 1):
            value += math.comb(n, k) * complete[n - k] * log_jets[k]
        complete.append(value)
    return complete


def endpoint_derivatives(y: arb, p: dict[str, Any], order: int, terms: int) -> list[acb]:
    F = centered_first_label_value(acb(y), p, terms, False)
    jets = bell_jet(log_derivatives(y, p, order - 1))
    return [F * jets[r] for r in range(order)]


def alternating_power_sum(power: int, count: int) -> arb:
    require(count > 0 and count % 2 == 0, "alternating sum expects a positive even count")
    half = count // 2
    if power == 1:
        return arb(half + 1).digamma() - arb(count + 1).digamma()
    exponent = arb(power)
    H_half = exponent.zeta() - exponent.zeta(arb(half + 1))
    H_full = exponent.zeta() - exponent.zeta(arb(count + 1))
    return arb(2) ** (1 - power) * H_half - H_full


def derivative_integral_bound(
    p: dict[str, Any], left: arb, right: arb, order: int, slabs: int
) -> tuple[arb, arb]:
    width = right - left
    integral_bound = arb(0)
    maximum_jet = arb(0)
    for index in range(slabs):
        slab_left = left + width * index / slabs
        slab_right = left + width * (index + 1) / slabs
        y = real_interval(slab_left, slab_right)
        jet = bell_jet(log_derivatives(y, p, order))[order]
        jet_bound = abs(jet).upper()
        F_bound = slab_left ** arb("-0.5")
        integral_bound += (slab_right - slab_left) * F_bound * jet_bound
        maximum_jet = max(maximum_jet, jet_bound)
    return integral_bound, maximum_jet


def transition_source_certificate(
    t: arb,
    endpoint: int,
    source_count: int,
    left: arb,
    right: arb,
    order: int,
    derivative_slabs: int,
    integral_panels: int,
    series_terms: int,
    tolerance: str,
) -> dict[str, Any]:
    require(source_count >= 3 and source_count % 2 == 1, "source count must be odd")
    later_count = source_count - 1
    require(later_count % 2 == 0, "later-label count parity drift")
    require((right - left).contains(int(float((right - left).mid()))), "cell length is not integral")
    p = centered_parameters(t, endpoint)
    first, panel_rows = first_label_integral(
        p, left, right, integral_panels, series_terms, tolerance
    )
    left_derivatives = endpoint_derivatives(left, p, order, series_terms)
    right_derivatives = endpoint_derivatives(right, p, order, series_terms)
    boundary = acb(0)
    alternating_rows: list[dict[str, Any]] = []
    imaginary = acb(0, 1)
    for r in range(order):
        power = r + 1
        alternating = alternating_power_sum(power, later_count)
        derivative_jump = right_derivatives[r] - left_derivatives[r]
        term = ((-1) ** r) * derivative_jump * alternating / (imaginary * 2 * p["pi"]) ** power
        boundary += term
        alternating_rows.append(
            {
                "power": power,
                "alternating_sum_ball": alternating.str(55, more=True),
                "boundary_term": complex_record(term, 55),
            }
        )
    derivative_bound, maximum_jet = derivative_integral_bound(
        p, left, right, order, derivative_slabs
    )
    harmonic_bound = arb(order).zeta()
    remainder = derivative_bound * harmonic_bound / (2 * p["pi"]) ** order
    complete = add_complex_error(first + boundary, remainder)
    return {
        "height_ball": t.str(55, more=True),
        "source_endpoint": endpoint,
        "source_count": source_count,
        "later_label_count": later_count,
        "cell_interval": [left.str(35, more=True), right.str(35, more=True)],
        "cell_length_ball": (right - left).str(35, more=True),
        "jet_order": order,
        "derivative_slabs": derivative_slabs,
        "integral_panels": integral_panels,
        "logarithm_series_terms": series_terms,
        "center_nu_ball": p["nu"].str(55, more=True),
        "center_detuning_ball": p["delta"].str(55, more=True),
        "first_label_integral_ball": complex_record(first),
        "boundary_jet_sum_ball": complex_record(boundary),
        "derivative_integral_bound_ball": derivative_bound.str(55, more=True),
        "maximum_normalized_order_jet_ball": maximum_jet.str(55, more=True),
        "later_label_remainder_bound_ball": remainder.str(55, more=True),
        "complete_transition_source_integral_ball": complex_record(complete),
        "alternating_rows": alternating_rows,
        "panel_rows": panel_rows,
    }


def joined_packet_certificate(
    source: dict[str, Any], dependencies: dict[str, dict[str, Any]]
) -> dict[str, Any]:
    pi, t, imaginary = arb.pi(), arb(HEIGHT), acb(0, 1)
    H = arb(dependencies["exact_H"]["exact_H_certificate"]["high_precision_H_ball"])
    B = complex_from_record(dependencies["gamma_bulk"]["certificate"]["Gamma_bulk_factor_ball"])
    theta_correction = arb(
        dependencies["gamma_bulk"]["certificate"]["exact_theta_minus_theta_zero_ball"]
    )
    theta_zero = t * ((t / (2 * pi)).log() - 1) / 2 - pi / 8
    theta_phase = (imaginary * (theta_zero + theta_correction)).exp()
    C_G = H * theta_phase.conjugate() * B * (imaginary * theta_zero).exp()
    s = acb(arb("0.5"), t)
    mode_sum = acb(0)
    for mode in range(TRANSITION_START, TRANSITION_END + 1):
        mode_sum += (-s * arb(mode).log()).exp()
    gamma_lift = C_G * mode_sum
    source_integral = complex_from_record(source["complete_transition_source_integral_ball"])
    natural_A = complex_from_record(
        dependencies["grouped_A"]["certificate"]["complex_preprojection"]
        ["natural_Hardy_lift_sum_ball"]
    )
    joined = source_integral - gamma_lift - natural_A
    hardy_joined = 2 * (theta_phase * joined).real
    return {
        "common_Gamma_coefficient_ball": complex_record(C_G),
        "transition_m_minus_s_sum_ball": complex_record(mode_sum),
        "transition_Gamma_lift_ball": complex_record(gamma_lift),
        "grouped_natural_A_lift_ball": complex_record(natural_A),
        "joined_transition_packet_ball": complex_record(joined),
        "joined_transition_Hardy_projection_ball": hardy_joined.str(70, more=True),
        "identity": "T_A^join=integral_(39852.5)^39936.5 y^(-s)e^(-i*pi*y^2)D_W(y)dy-C_G*sum_(39853..39936)m^(-s)-mathcal_A_A^nat",
    }


def render_note(artifact: dict[str, Any]) -> str:
    s = artifact["source_certificate"]
    j = artifact["joined_packet"]
    return f"""# Alternating boundary-jet transition packet

Date: 2026-08-27

Status: rigorous grouped Fresnel/Gamma/A transition packet at `t=10^10`;
the full `P_W`/unowned-cell join remains open.

On the contiguous transition interval

```text
[a,b]=[39852.5,39936.5],       b-a=84,
F(y)=y^(-s)exp(i*pi*(159577*y-y^2)),
```

the finite source is

```text
sum_(j=0)^(M-1) integral_a^b F(y)exp(2*pi*i*j*y)dy,
M=2481423.
```

Because both endpoints are half integers and their difference is integral,
`exp(2*pi*i*j*a)=exp(2*pi*i*j*b)=(-1)^j`.  Repeated integration by parts gives
the exact finite expansion

```text
I_transition=I_0
 +sum_(r=0)^(R-1) (-1)^r [F^(r)(b)-F^(r)(a)]
    A_(r+1)(M-1)/(2*pi*i)^(r+1)
 +R_R,

A_p(N)=sum_(j=1)^N (-1)^j/j^p.
```

The remainder is bounded without summing the source roster:

```text
|R_R| <= zeta(R)/(2*pi)^R integral_a^b |F^(R)(y)|dy.
```

With `R={s['jet_order']}` and {s['derivative_slabs']} derivative slabs,

```text
|R_R| <= {s['later_label_remainder_bound_ball']}.
```

The complete transition source integral is

```text
{s['complete_transition_source_integral_ball']['real_ball']}
+ i {s['complete_transition_source_integral_ball']['imag_ball']}.
```

Subtracting the exact extended Gamma lift and the grouped natural A lift before
projection gives

```text
T_A^join = {j['joined_transition_packet_ball']['real_ball']}
          +i {j['joined_transition_packet_ball']['imag_ball']},

Hardy_t[T_A^join] = {j['joined_transition_Hardy_projection_ball']}.
```

Pi provenance: every `pi` comes from the finite Fresnel source, its exact
integer Fourier spacing, the Gamma normalization, or the Riemann-Siegel
phase.  No fitted constant is introduced.

Proof boundary: rigorous transition block at the one saved height only.  This
gate does not join `P_W` to cells `1..621` and `39937..infinity`, bound the
ordinary packet, enclose complete `J_Z` or `D_K`, or prove a non-A, all-height,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    ctx.dps = 105
    ctx.threads = 1
    paths = {
        "grouped_A": GROUPED_A,
        "common_cell": COMMON_CELL,
        "Fresnel_connection": FRESNEL_CONNECTION,
        "exact_H": EXACT_H,
        "gamma_bulk": GAMMA_BULK,
    }
    dependencies = {name: load_json(path) for name, path in paths.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency is not passed")
    require(
        dependencies["grouped_A"]["decision"]["grouped_natural_A_Hardy_lift_certified"] is True,
        "grouped A dependency drift",
    )
    require(
        dependencies["common_cell"]["decision"]["finite_Fresnel_cell_partition_exact"] is True,
        "cell-partition dependency drift",
    )
    source = transition_source_certificate(
        arb(HEIGHT),
        SOURCE_ENDPOINT,
        SOURCE_COUNT,
        CELL_LEFT,
        CELL_RIGHT,
        JET_ORDER,
        DERIVATIVE_SLABS,
        INTEGRAL_PANELS,
        SERIES_TERMS,
        "1e-34",
    )
    require(
        arb(source["later_label_remainder_bound_ball"]).upper() < arb("1e-12"),
        "boundary-jet remainder misses transition cap",
    )
    joined = joined_packet_certificate(source, dependencies)

    builder = Path(__file__).resolve()
    checker = builder.with_name("check_" + builder.name)
    artifact = {
        "kind": STEM,
        "status": "actual_height_grouped_Fresnel_Gamma_A_transition_packet_interval_certified",
        "passed": True,
        "source_certificate": source,
        "joined_packet": joined,
        "decision": {
            "all_2481423_source_labels_collapsed_exactly_on_transition_interval": True,
            "alternating_boundary_jet_remainder_rigorously_bounded": True,
            "grouped_Fresnel_Gamma_A_transition_packet_enclosed": True,
            "P_W_unowned_cell_Mordell_join_proved": False,
            "ordinary_owned_packet_bounded": False,
            "actual_height_J_Z_enclosed": False,
            "actual_height_D_K_enclosed": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in paths.items()
        },
        "sources": {
            "builder": {"path": relative(builder), "sha256": file_hash(builder)},
            "checker": {"path": relative(checker), "sha256": file_hash(checker)},
        },
        "runtime": {
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.time() - started, 3),
        },
        "next_obligation": "Extend the same endpoint-jet/Mordell bookkeeping to join P_W with cells 1..621 and 39937..infinity, and enclose the ordinary Gamma-subtracted packet before assembling J_Z.",
        "proof_boundary": "Rigorous grouped Fresnel/Gamma/A transition packet at t=10^10 only. No P_W/unowned-cell Mordell join, ordinary packet bound, complete J_Z or D_K enclosure, non-A, all-height, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified alternating boundary-jet transition packet", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
