#!/usr/bin/env python3
"""Build the centered first-order adjacent-saddle recurrence."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_centered_adjacent_saddle_recurrence"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "centered_reduction": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_real_residual_reduction.json"
    ),
    "endpoint_recurrence": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "endpoint_C0_recurrence_transition_gate.json"
    ),
    "global_first_order": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_"
        "first_order_global_remainder_certificate.json"
    ),
}
LOW_DPS = 110
HIGH_DPS = 155
N_VALUES = (100_000_000_000, 10_000_000_000_000)
THETA_VALUES = ("0", "0.5", "0.9")
C_VALUES = ("0", "1", "25")


@dataclass(frozen=True)
class RecurrenceRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def source_hashes() -> dict[str, str]:
    missing = [str(path) for path in SOURCES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {key: file_hash(path) for key, path in SOURCES.items()}


def source_audit() -> dict[str, str]:
    payloads = {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCES.items()
    }
    markers = {
        "centered_reduction": (
            "A_a=-(1/2)Im(M_(1,a))",
            "|U-A_a|<1e-6*exp(-5L/4)",
        ),
        "endpoint_recurrence": (
            "C0(p+2)+C0(p)=exp(pi*i*(p^2/2+p+3/8))",
            "one added two-saddle Dirichlet block",
        ),
        "global_first_order": (
            "|Delta lift_[1]|/A_t(x)<20000*exp(-5L/4)",
            "|r_[1](x)|<100000*exp(-5L/4)",
        ),
    }
    for key, required in markers.items():
        text = json.dumps(payloads[key], sort_keys=True)
        for marker in required:
            if marker not in text:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {key: str(payload.get("kind", "")) for key, payload in payloads.items()}


def symbolic_audit() -> dict[str, str]:
    p, a = sp.symbols("p a", positive=True, real=True)
    y = p + 1
    phase = sp.pi * sp.I * (
        p**2 / 2 + p + sp.Rational(3, 8)
    )
    r_phase = sp.exp(phase)
    r3 = sp.diff(r_phase, p, 3)
    r4 = sp.diff(r_phase, p, 4)
    expected_r3 = (
        -3 * sp.pi**2 * y - sp.I * sp.pi**3 * y**3
    ) * r_phase
    expected_r4 = (
        sp.pi**4 * y**4
        - 6 * sp.I * sp.pi**3 * y**2
        - 3 * sp.pi**2
    ) * r_phase
    if sp.simplify(r3 - expected_r3) != 0:
        raise RuntimeError("R third derivative identity failed")
    if sp.simplify(r4 - expected_r4) != 0:
        raise RuntimeError("R fourth derivative identity failed")

    h_sum = sp.simplify(r_phase + r3 / (12 * sp.pi**2 * a))
    displacement = sp.symbols("r", real=True)
    h_ratio = sp.simplify(
        (h_sum / r_phase).subs(p, 2 * a * displacement - 1)
    )
    expected_h_ratio = (
        1
        - displacement / 2
        - 2 * sp.I * sp.pi * a**2 * displacement**3 / 3
    )
    if sp.simplify(h_ratio - expected_h_ratio) != 0:
        raise RuntimeError("corrected endpoint recurrence failed")

    p_x = -1 / (4 * sp.pi * a)
    a_x = 1 / (8 * sp.pi * a)
    h_x = sp.simplify(
        p_x
        * (
            sp.diff(r_phase, p)
            + r4 / (12 * sp.pi**2 * a)
        )
        - a_x * r3 / (12 * sp.pi**2 * a**2)
    )
    h_x_ratio = sp.simplify(
        (h_x / r_phase).subs(p, 2 * a * displacement - 1)
    )
    expected_h_x_ratio = (
        (1 + displacement) / (16 * sp.pi * a**2)
        - sp.pi * a**2 * displacement**4 / 3
        + sp.I
        * (
            -displacement / 2
            + displacement**2 / 2
            + displacement**3 / 12
        )
    )
    if sp.simplify(h_x_ratio - expected_h_x_ratio) != 0:
        raise RuntimeError("corrected endpoint derivative recurrence failed")

    lam, mu, saddle_x, ell = sp.symbols(
        "lambda mu s_x ell"
    )
    f, e, d_x, kappa, h_block, h_block_x = sp.symbols(
        "f e d_x kappa H H_x"
    )
    f_x = lam * f - saddle_x * ell * f + e * d_x
    kappa_x = (lam + mu) * kappa
    q = f + kappa * h_block
    q_x = f_x + kappa_x * h_block + kappa * h_block_x
    scalar_complex = (
        -saddle_x * ell * f
        + kappa * (h_block_x + mu * h_block)
    )
    if sp.simplify(
        scalar_complex - (q_x - lam * q - e * d_x)
    ) != 0:
        raise RuntimeError("centered adjacent scalar identity failed")

    epsilon, w, time = sp.symbols(
        "epsilon w time", positive=True, real=True
    )
    shift = sp.Rational(1, 2) - sp.I * sp.pi * time / 8
    log_displacement_first = w
    saddle_defect_first = 2 * w**3 / 3
    omega_first = sp.simplify(
        -shift * log_displacement_first
        - sp.I * sp.pi * time * log_displacement_first / 8
        - sp.I * sp.pi * saddle_defect_first
    )
    endpoint_first = -w / 2 - 2 * sp.I * sp.pi * w**3 / 3
    if sp.simplify(omega_first - endpoint_first) != 0:
        raise RuntimeError("first adjacent ratio coefficient did not cancel")

    return {
        "R_derivatives": (
            "R'''=(-3*pi^2*(p+1)-i*pi^3*(p+1)^3)*R; "
            "R''''=(pi^4*(p+1)^4-6i*pi^3*(p+1)^2-3*pi^2)*R"
        ),
        "endpoint_sum": (
            "J_a=H_a(p+2)+H_a(p)=R(p)*"
            "[1-r/2-(2*pi*i/3)*a^2*r^3], r=(N+1-a)/a"
        ),
        "endpoint_derivative": (
            "J_(a,x)/R=(1+r)/(16*pi*a^2)-pi*a^2*r^4/3"
            "+i*(-r/2+r^2/2+r^3/12)"
        ),
        "scalar_jump": (
            "Delta A_a=Re[-s_*'*ell*f_(N+1)+"
            "kappa_N*(J_(a,x)+mu_a*J_a)]"
        ),
        "mismatch_identity": (
            "Q_N=f_(N+1)+kappa_N*J_a=E_[1],N+1-E_[1],N; "
            "Delta A_a=Re(Q_(N,x)-lambda_a*Q_N-e_(N+1)*d_(N+1,x))"
        ),
        "real_projection": (
            "Delta A_a=Delta U-u_a*Delta X+v_a*Delta Y"
            "-Re(e_(N+1)*d_(N+1,x))"
        ),
        "first_coefficient": (
            "For epsilon=1/a and w=N+1-a, "
            "[epsilon]log(main_sharp/endpoint_C1)=0"
        ),
        "formal_scale": (
            "main_sharp/endpoint_C1=1+O(a^-2); after centered "
            "x-differentiation Delta A_a=O(a^-7/2)=O(exp(-7L/4))"
        ),
    }


def alpha(s: mp.mpc) -> mp.mpc:
    return 1 / (2 * s) + 1 / (s - 1) + mp.log(s / (2 * mp.pi)) / 2


def alpha_prime(s: mp.mpc) -> mp.mpc:
    return -1 / (2 * s**2) - 1 / (s - 1) ** 2 + 1 / (2 * s)


def alpha_second(s: mp.mpc) -> mp.mpc:
    return 1 / s**3 + 2 / (s - 1) ** 3 - 1 / (2 * s**2)


def m_zero(s: mp.mpc) -> mp.mpc:
    return (
        mp.mpf(1)
        / 8
        * (s * (s - 1) / 2)
        * mp.pi ** (-s / 2)
        * mp.sqrt(2 * mp.pi)
        * mp.exp(
            (s / 2 - mp.mpf("0.5")) * mp.log(s / 2) - s / 2
        )
    )


def m_time(s: mp.mpc, time: mp.mpf) -> mp.mpc:
    alpha_value = alpha(s)
    return mp.exp(time * alpha_value**2 / 4) * m_zero(s)


def boundary_time(a: mp.mpf, c_value: mp.mpf) -> mp.mpf:
    if c_value == 0:
        return mp.mpf(0)
    initial = c_value / (2 * mp.log(a))
    return mp.findroot(
        lambda value: value * mp.log(a**2 - value / 16) - c_value,
        initial,
    )


def diagnostic_row(
    n_cutoff: int,
    theta_text: str,
    c_text: str,
) -> dict[str, str | int]:
    theta = mp.mpf(theta_text)
    c_value = mp.mpf(c_text)
    a = mp.mpf(n_cutoff) + theta
    time = boundary_time(a, c_value)
    x = 4 * mp.pi * a**2 - mp.pi * time / 4
    ell_height = mp.log(x / (4 * mp.pi))
    temperature = 2 * mp.pi * a**2

    s = (1 - mp.j * x) / 2
    alpha_value = alpha(s)
    alpha_one = alpha_prime(s)
    alpha_two = alpha_second(s)
    s_star = s + time * alpha_value / 2
    s_star_x = -mp.j / 2 - mp.j * time * alpha_one / 4
    normalizer = m_time(s, time)
    phase = normalizer / abs(normalizer)

    entering_n = n_cutoff + 1
    log_n = mp.log(entering_n)
    centered_log = mp.log(entering_n / a)
    alpha_n = alpha_value - log_n
    e_n = phase * mp.exp(
        time * log_n**2 / 4 - s_star * log_n
    )
    d_n = (
        1 / (6 * s)
        + alpha_one
        * (time / 4 + time**2 * alpha_n**2 / 8)
    )
    d_n_x = (
        -mp.j
        / 2
        * (
            -1 / (6 * s**2)
            + alpha_two
            * (time / 4 + time**2 * alpha_n**2 / 8)
            + (time**2 / 4) * alpha_n * alpha_one**2
        )
    )
    f_n = e_n * (1 + d_n)

    p = 1 - 2 * (a - n_cutoff)
    y = p + 1
    relative = (entering_n - a) / a
    recurrence_phase = mp.exp(
        mp.pi * mp.j * (p**2 / 2 + p + mp.mpf(3) / 8)
    )
    endpoint_sum = recurrence_phase * (
        1
        - y / (4 * a)
        - mp.j * mp.pi * y**3 / (12 * a)
    )
    endpoint_sum_x = recurrence_phase * (
        (1 + relative) / (16 * mp.pi * a**2)
        - mp.pi * a**2 * relative**4 / 3
        + mp.j
        * (
            -relative / 2
            + relative**2 / 2
            + relative**3 / 12
        )
    )

    u_rs = mp.exp(
        -mp.j
        * (
            temperature
            / 2
            * mp.log(temperature / (2 * mp.pi))
            - temperature / 2
            - mp.pi / 8
        )
    )
    kappa = (
        (-1) ** n_cutoff
        * mp.exp(time * mp.pi**2 / 64)
        * m_zero(mp.j * temperature)
        * u_rs
        * mp.exp(mp.j * mp.pi / 8)
        / abs(normalizer)
    )
    beta_prime = -mp.re(
        alpha_value + time * alpha_value * alpha_one / 2
    ) / 2
    lambda_a = mp.j * beta_prime - s_star_x * mp.log(a)
    kappa_log_x = (
        -mp.pi / 8
        + 1 / (4 * temperature)
        + 1 / (2 * (temperature + mp.j))
    )
    m_time_log_x = (
        -mp.j
        / 2
        * (
            alpha_value
            + time * alpha_value * alpha_one / 2
        )
    )
    mu_a = (
        kappa_log_x - m_time_log_x + s_star_x * mp.log(a)
    )

    finite_piece = -s_star_x * centered_log * f_n
    endpoint_piece = kappa * (
        endpoint_sum_x + mu_a * endpoint_sum
    )
    scalar_jump = mp.re(finite_piece + endpoint_piece)
    scale_three = mp.exp(
        3 * ell_height / 4 + c_value * ell_height / 16
    )
    scale_seven = mp.exp(
        7 * ell_height / 4 + c_value * ell_height / 16
    )
    leading_denominator = (
        abs(mp.re(finite_piece)) + abs(mp.re(endpoint_piece))
    )
    cancellation_ratio = (
        abs(scalar_jump) / leading_denominator
        if leading_denominator
        else mp.mpf(0)
    )
    return {
        "N": n_cutoff,
        "theta": theta_text,
        "c": c_text,
        "L": mp.nstr(ell_height, 50),
        "t": mp.nstr(time, 50),
        "finite_real_scaled_e3": mp.nstr(
            mp.re(finite_piece) * scale_three, 50
        ),
        "endpoint_real_scaled_e3": mp.nstr(
            mp.re(endpoint_piece) * scale_three, 50
        ),
        "scalar_jump_scaled_e7": mp.nstr(
            scalar_jump * scale_seven, 50
        ),
        "explicit_d_scaled_e7": mp.nstr(
            abs(mp.re(e_n * d_n_x)) * scale_seven, 50
        ),
        "cancellation_ratio": mp.nstr(cancellation_ratio, 50),
        "lambda_abs": mp.nstr(abs(lambda_a), 50),
    }


def selected_diagnostics(
    low_dps: int = LOW_DPS,
    high_dps: int = HIGH_DPS,
) -> dict:
    def run(dps: int) -> list[dict]:
        mp.mp.dps = dps
        return [
            diagnostic_row(n_cutoff, theta, c_value)
            for n_cutoff in N_VALUES
            for theta in THETA_VALUES
            for c_value in C_VALUES
        ]

    low_rows = run(low_dps)
    high_rows = run(high_dps)
    keys = ("N", "theta", "c")
    if [
        tuple(row[key] for key in keys) for row in low_rows
    ] != [
        tuple(row[key] for key in keys) for row in high_rows
    ]:
        raise RuntimeError("diagnostic row keys drifted")

    numeric_keys = (
        "L",
        "t",
        "finite_real_scaled_e3",
        "endpoint_real_scaled_e3",
        "scalar_jump_scaled_e7",
        "explicit_d_scaled_e7",
        "cancellation_ratio",
        "lambda_abs",
    )
    max_relative_delta = mp.mpf(0)
    for low, high in zip(low_rows, high_rows, strict=True):
        for key in numeric_keys:
            low_value = mp.mpf(low[key])
            high_value = mp.mpf(high[key])
            denominator = max(abs(high_value), mp.mpf("1e-80"))
            max_relative_delta = max(
                max_relative_delta,
                abs(low_value - high_value) / denominator,
            )

    max_leading = max(
        max(
            abs(mp.mpf(row["finite_real_scaled_e3"])),
            abs(mp.mpf(row["endpoint_real_scaled_e3"])),
        )
        for row in high_rows
    )
    max_scalar = max(
        abs(mp.mpf(row["scalar_jump_scaled_e7"]))
        for row in high_rows
    )
    max_explicit_d = max(
        mp.mpf(row["explicit_d_scaled_e7"])
        for row in high_rows
    )
    max_cancellation_ratio = max(
        mp.mpf(row["cancellation_ratio"])
        for row in high_rows
    )
    if max_relative_delta >= mp.mpf("1e-45"):
        raise RuntimeError("diagnostic precision drift exceeded 1e-45")
    if max_leading >= mp.mpf("0.24"):
        raise RuntimeError("selected leading scale exceeded 0.24")
    if max_scalar >= 2:
        raise RuntimeError("selected scalar residual exceeded 2")
    return {
        "role": "selected_high_precision_diagnostic_only",
        "proof_boundary": (
            "Eighteen selected rows at two L>=50 heights and three saddle "
            "positions are not a uniform theta or parameter theorem."
        ),
        "low_dps": low_dps,
        "high_dps": high_dps,
        "row_count": len(high_rows),
        "max_relative_delta": mp.nstr(max_relative_delta, 30),
        "max_leading_scaled_e3": mp.nstr(max_leading, 30),
        "max_scalar_jump_scaled_e7": mp.nstr(max_scalar, 30),
        "max_explicit_d_scaled_e7": mp.nstr(max_explicit_d, 30),
        "max_cancellation_ratio": mp.nstr(
            max_cancellation_ratio, 30
        ),
        "rows": high_rows,
    }


def build_artifact() -> dict:
    exact = symbolic_audit()
    diagnostics = selected_diagnostics()
    rows = [
        RecurrenceRow(
            id="nfocasr_01_recurrence_derivatives",
            role="exact_identity",
            readiness="ready_to_apply",
            claim=(
                "The adjacent C0 recurrence differentiates exactly through "
                "the orders needed by the retained C1 endpoint."
            ),
            formula=exact["R_derivatives"],
            proof_boundary="Algebraic identity for the entire recurrence phase.",
        ),
        RecurrenceRow(
            id="nfocasr_02_endpoint_sum",
            role="exact_identity",
            readiness="ready_to_apply",
            claim=(
                "The two adjacent corrected endpoints collapse to one "
                "explicit recurrence block."
            ),
            formula=exact["endpoint_sum"],
            proof_boundary="Fixed a and adjacent integer lifts N,N+1.",
        ),
        RecurrenceRow(
            id="nfocasr_03_endpoint_derivative",
            role="exact_identity",
            readiness="ready_to_apply",
            claim=(
                "The x derivative of the recurrence block is an explicit "
                "quartic polynomial in the saddle displacement."
            ),
            formula=exact["endpoint_derivative"],
            proof_boundary="Uses a_x=1/(8*pi*a) and p_x=-1/(4*pi*a).",
        ),
        RecurrenceRow(
            id="nfocasr_04_scalar_jump",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim=(
                "The adjacent change in A_a is one entering centered moment "
                "plus one differentiated recurrence endpoint."
            ),
            formula=exact["scalar_jump"],
            proof_boundary="No division by E_[1], a Wronskian, or a jet norm.",
        ),
        RecurrenceRow(
            id="nfocasr_05_mismatch_identity",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim=(
                "The two leading pieces are exactly the centered derivative "
                "of the adjacent complex-lift mismatch."
            ),
            formula=exact["mismatch_identity"],
            proof_boundary=(
                "The complex mismatch need not be small; only its real trace "
                "and centered real derivative enter A_a."
            ),
        ),
        RecurrenceRow(
            id="nfocasr_06_real_projection",
            role="exact_identity",
            readiness="ready_to_apply",
            claim=(
                "The chart jump separates into the adjacent real value/jet "
                "and already-controlled frame nuisances."
            ),
            formula=exact["real_projection"],
            proof_boundary="Real critical line only.",
        ),
        RecurrenceRow(
            id="nfocasr_07_first_coefficient",
            role="asymptotic_identity",
            readiness="ready_to_apply",
            claim=(
                "The order-a^-1 saddle phase, heat shift, logarithmic weight, "
                "and C1 endpoint coefficient cancel before absolute values."
            ),
            formula=exact["first_coefficient"],
            proof_boundary=(
                "Formal coefficient identity from the exact positive-sharp "
                "ratio; a uniform remainder constant is separate."
            ),
        ),
        RecurrenceRow(
            id="nfocasr_08_formal_scale",
            role="asymptotic_reduction",
            readiness="conditional_ready",
            claim=(
                "After the vanished first coefficient, the natural scalar "
                "chart-jump scale is two powers of a below the edge block."
            ),
            formula=exact["formal_scale"],
            proof_boundary=(
                "Conditional on differentiating the exact ratio with a "
                "uniform O(a^-2) remainder."
            ),
        ),
        RecurrenceRow(
            id="nfocasr_09_selected_audit",
            role="high_precision_diagnostic",
            readiness="diagnostic_only",
            claim=(
                "Selected L>=50 rows show equal-and-opposite e^-3L/4 real "
                "pieces and an e^-7L/4 centered scalar residual."
            ),
            formula=(
                "max edge-scaled piece <0.24; "
                "max |Delta A_a|*exp(7L/4+cL/16)<2"
            ),
            proof_boundary=diagnostics["proof_boundary"],
            diagnostics={
                "rows": diagnostics["row_count"],
                "max_edge": diagnostics["max_leading_scaled_e3"],
                "max_scalar": diagnostics[
                    "max_scalar_jump_scaled_e7"
                ],
                "precision_delta": diagnostics["max_relative_delta"],
            },
        ),
        RecurrenceRow(
            id="nfocasr_10_route_audit",
            role="strategy_reduction",
            readiness="ready_to_apply",
            claim=(
                "The last saddle and endpoint cannot be bounded separately "
                "to obtain a useful sign; their chart-invariant combination "
                "must be retained."
            ),
            formula=(
                "absolute values at e^-3L/4 destroy the exact adjacent "
                "cancellation"
            ),
            proof_boundary=(
                "This removes an invalid dominance route; it does not supply "
                "the bulk scalar sign."
            ),
        ),
        RecurrenceRow(
            id="nfocasr_11_uniform_target",
            role="open_theorem_target",
            readiness="not_ready_to_apply",
            claim=(
                "A uniform differentiated positive-sharp ratio bound would "
                "make the centered scalar independent of the adjacent chart "
                "below the contact scale."
            ),
            formula=(
                "Prove |rho-1|<=C0/a^2 and |rho_x|<=C1/a^3 for "
                "0<=N+1-a<=1, then |Delta A_a|<=C_A*exp(-7L/4) "
                "with an explicit C_A"
            ),
            proof_boundary=(
                "Open uniform estimate on L>=50 and 0<tL<=25; selected "
                "rows do not prove it."
            ),
        ),
        RecurrenceRow(
            id="nfocasr_12_bulk_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim=(
                "Once chart stability is certified, the RH-level work is a "
                "signed small-value/crossing theorem for the chart-invariant "
                "bulk centered moment."
            ),
            formula=(
                "|A_a|>(100000L+1)*exp(-5L/4) and an A_a-positive "
                "successor count below one turn"
            ),
            proof_boundary=(
                "The adjacent recurrence supplies no lower bound or signed "
                "crossing count."
            ),
        ),
        RecurrenceRow(
            id="nfocasr_13_nonpromotion",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim=(
                "An exact edge cancellation and finite high-height rows are "
                "not a transversality theorem."
            ),
            formula=(
                "adjacent-saddle recurrence != contact exclusion"
            ),
            proof_boundary=(
                "Not the q>=1 scalar theorem, q<1 theorem, Lambda<=0, RH, "
                "PF-infinity, or a Clay-prize conclusion."
            ),
        ),
    ]
    return {
        "kind": STEM,
        "date": "2026-07-26",
        "status": (
            "exact centered adjacent-saddle recurrence, first-coefficient "
            "cancellation, and selected L>=50 diagnostics; not a uniform "
            "chart-stability theorem, contact exclusion, Lambda<=0, or RH"
        ),
        "proof_boundary": (
            "This artifact proves the differentiated C0 recurrence, corrected "
            "endpoint sum and derivative, exact adjacent A_a jump identities, "
            "and the vanishing order-a^-1 coefficient. The selected diagnostics "
            "do not prove the uniform differentiated ratio bound. It does not "
            "prove the q>=1 scalar lower bound or signed crossing budget, the "
            "q<1 multiplicity-compatible theorem, finite phase cells, one-sided "
            "successor winding, contact exclusion, Lambda<=0, RH, PF-infinity, "
            "or a Clay-prize conclusion."
        ),
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "exact": exact,
        "diagnostics": diagnostics,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    diagnostics = artifact["diagnostics"]
    lines = [
        "# Jensen-Window PF Newman Polymath-15 First-Order Centered Adjacent-Saddle Recurrence",
        "",
        "Date: 2026-07-26",
        "",
        "Status: exact adjacent-saddle identities and selected high-height",
        "diagnostics. The uniform differentiated-ratio theorem remains open;",
        "this is not a proof of `Lambda<=0` or RH.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"python work/rh_compute/scripts/{STEM}.py",
        f"python work/rh_compute/scripts/check_{STEM}.py",
        "```",
        "",
        "## Exact Adjacent Recurrence",
        "",
        "Put `n=N+1`, `p=1-2(a-N)`,",
        "`r=(n-a)/a=(p+1)/(2a)`, and",
        "`R(p)=exp(pi*i*(p^2/2+p+3/8))`. Differentiating the exact",
        "`C_0` recurrence gives",
        "",
        "```text",
        exact["R_derivatives"],
        "",
        exact["endpoint_sum"],
        "",
        exact["endpoint_derivative"],
        "```",
        "",
        "If `H_a=F+F'''/(12pi^2a)` and `g_N=-kappa_NH_a(p)`,",
        "then `kappa_(N+1)=-kappa_N` and",
        "",
        "```text",
        "g_(N+1)-g_N=kappa_N*J_a.",
        "```",
        "",
        "No subtraction of two long finite sums occurs.",
        "",
        "## Centered Scalar Jump",
        "",
        "The moment changes by",
        "`M_(1,a;N+1)-M_(1,a;N)=log((N+1)/a)f_(N+1)`.",
        "Consequently",
        "",
        "```text",
        exact["scalar_jump"],
        "",
        exact["mismatch_identity"],
        "",
        exact["real_projection"],
        "```",
        "",
        "The complex half itself is chart dependent: its adjacent mismatch",
        "can have a large imaginary part. The exact real trace and its real",
        "derivative are the relevant quantities; replacing them by a complex",
        "absolute-value estimate would throw away the cancellation.",
        "",
        "## Leading Cancellation",
        "",
        "For the reflected positive sharp block, set `epsilon=1/a` and",
        "`w=N+1-a`. The exact adjacent log-ratio has first coefficient",
        "",
        "```text",
        "-(1/2-i*pi*t/8)w-i*pi*t*w/8-i*(2*pi/3)w^3",
        " =-w/2-i*(2*pi/3)w^3,",
        "```",
        "",
        "which is exactly the first coefficient of",
        "`log(1-w/(2a)-i*(2pi/3)w^3/a)`. The heat shift cancels too.",
        "The finite correction starts only at `a^-2`. Therefore",
        "",
        "```text",
        exact["first_coefficient"],
        exact["formal_scale"],
        "```",
        "",
        "This explains why the separate finite and endpoint pieces live at",
        "`exp(-3L/4)` while their centered real combination is expected at",
        "`exp(-7L/4)`. A uniform remainder constant is still required.",
        "",
        "## Selected High-Height Audit",
        "",
        (
            f"Independent {diagnostics['low_dps']}- and "
            f"{diagnostics['high_dps']}-digit runs agree with maximum "
            f"relative drift `{diagnostics['max_relative_delta']}`."
        ),
        "",
        "| N | theta | c=tL | L | finite e3 | endpoint e3 | Delta A e7 |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in diagnostics["rows"]:
        lines.append(
            "| {N} | {theta} | {c} | {L} | {finite} | {endpoint} | {scalar} |".format(
                N=row["N"],
                theta=row["theta"],
                c=row["c"],
                L=mp.nstr(mp.mpf(row["L"]), 10),
                finite=mp.nstr(
                    mp.mpf(row["finite_real_scaled_e3"]), 9
                ),
                endpoint=mp.nstr(
                    mp.mpf(row["endpoint_real_scaled_e3"]), 9
                ),
                scalar=mp.nstr(
                    mp.mpf(row["scalar_jump_scaled_e7"]), 10
                ),
            )
        )
    lines.extend(
        [
            "",
            "Here `finite e3` and `endpoint e3` multiply the real pieces by",
            "`exp(3L/4+cL/16)`, while `Delta A e7` multiplies their sum by",
            "`exp(7L/4+cL/16)`. Across these selected rows,",
            "",
            "```text",
            (
                "max edge-scaled magnitude = "
                f"{diagnostics['max_leading_scaled_e3']} < 0.24"
            ),
            (
                "max scalar-jump magnitude = "
                f"{diagnostics['max_scalar_jump_scaled_e7']} < 2"
            ),
            "```",
            "",
            "These are diagnostics at eighteen prescribed points, not an",
            "interval proof over `theta`, `L`, or `t`.",
            "",
            "## Route Decision",
            "",
            "The edge contribution cannot be assigned a sign by bounding the",
            "last centered finite term and endpoint separately. Their",
            "`exp(-3L/4)` pieces are equal-and-opposite in the proof-facing",
            "real coordinate. The adjacent recurrence must be differentiated",
            "before absolute values are taken.",
            "",
            "The next technical theorem is explicit:",
            "",
            "```text",
            "Prove |rho-1|<=C0/a^2 and |rho_x|<=C1/a^3",
            "uniformly for L>=50, 0<tL<=25, 0<=N+1-a<=1.",
            "Deduce |A_(N+1)-A_N|<=C_A*exp(-7L/4)",
            "with any explicit constant small at the exp(-5L/4) contact scale.",
            "```",
            "",
            "Once that chart-stability lemma is certified, the remaining",
            "RH-level object is the chart-invariant bulk scalar lower bound and",
            "its oriented positive-crossing count. This artifact supplies",
            "neither theorem and does not prove contact exclusion, `Lambda<=0`,",
            "RH, PF-infinity, or a Clay-prize conclusion.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman first-order centered adjacent-saddle recurrence: "
        "13 rows, exact Delta A identity, vanished a^-1 coefficient, "
        "18 selected L>=50 diagnostics, 2 open scalar obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
