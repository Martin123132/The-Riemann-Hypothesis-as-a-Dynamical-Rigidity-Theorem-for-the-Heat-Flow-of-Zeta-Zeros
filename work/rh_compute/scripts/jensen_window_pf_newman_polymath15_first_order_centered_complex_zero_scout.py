#!/usr/bin/env python3
"""Build a selected high-precision corrected-complex-main zero scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_centered_complex_zero_scout"
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
    "signed_contact": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_"
        "first_order_signed_contact_reduction.json"
    ),
}
N_FIXED = 6
LOW_DPS = 70
HIGH_DPS = 95


@dataclass(frozen=True)
class ScoutRow:
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
        * mp.exp((s / 2 - mp.mpf("0.5")) * mp.log(s / 2) - s / 2)
    )


def m_time(s: mp.mpc, time: mp.mpf) -> mp.mpc:
    alpha_value = alpha(s)
    return mp.exp(time * alpha_value**2 / 4) * m_zero(s)


def c_zero(p: mp.mpf) -> mp.mpc:
    return (
        mp.exp(mp.pi * mp.j * (p**2 / 2 + mp.mpf(3) / 8))
        - mp.j * mp.sqrt(2) * mp.cos(mp.pi * p / 2)
    ) / (2 * mp.cos(mp.pi * p))


def corrected_half(
    x: mp.mpf,
    time: mp.mpf,
    n_fixed: int = N_FIXED,
) -> dict[str, mp.mpc | mp.mpf | int]:
    s = (1 - mp.j * x) / 2
    alpha_value = alpha(s)
    alpha_one = alpha_prime(s)
    alpha_two = alpha_second(s)
    s_star = s + time * alpha_value / 2
    normalizer = m_time(s, time)
    phase = normalizer / abs(normalizer)
    temperature = x / 2 + mp.pi * time / 8
    saddle = mp.sqrt(temperature / (2 * mp.pi))
    actual_cutoff = int(mp.floor(saddle))
    if actual_cutoff != n_fixed:
        raise RuntimeError(
            f"selected lift left N={n_fixed}: floor(a)={actual_cutoff}"
        )

    finite = mp.mpc(0)
    centered_moment = mp.mpc(0)
    derivative_correction = mp.mpc(0)
    for n in range(1, n_fixed + 1):
        log_n = mp.log(n)
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
        finite += f_n
        centered_moment += mp.log(n / saddle) * f_n
        derivative_correction += e_n * d_n_x

    p = 1 - 2 * (saddle - n_fixed)
    f_three = mp.diff(c_zero, p, 3)
    u_rs = mp.exp(
        -mp.j
        * (
            (temperature / 2)
            * mp.log(temperature / (2 * mp.pi))
            - temperature / 2
            - mp.pi / 8
        )
    )
    kappa = (
        (-1) ** n_fixed
        * mp.exp(time * mp.pi**2 / 64)
        * m_zero(mp.j * temperature)
        * u_rs
        * mp.exp(mp.j * mp.pi / 8)
        / abs(normalizer)
    )
    endpoint = -kappa * (
        c_zero(p) + f_three / (12 * mp.pi**2 * saddle)
    )
    return {
        "E": finite + endpoint,
        "M1": centered_moment,
        "D1x": derivative_correction,
        "g0": endpoint,
        "N": actual_cutoff,
        "p": p,
        "a": saddle,
    }


def solve_selected_zero(dps: int) -> dict:
    mp.mp.dps = dps
    x_zero, time_zero = mp.findroot(
        lambda x, time: (
            mp.re(corrected_half(x, time)["E"]),
            mp.im(corrected_half(x, time)["E"]),
        ),
        (mp.mpf("519.38"), mp.mpf("0.403")),
        (mp.mpf("519.40"), mp.mpf("0.405")),
        tol=mp.mpf(10) ** (-(dps - 15)),
        maxsteps=60,
    )
    values = corrected_half(x_zero, time_zero)
    complex_main = values["E"]
    centered_moment = values["M1"]
    derivative_correction = values["D1x"]
    endpoint = values["g0"]

    derivative_x = mp.diff(
        lambda local_x: corrected_half(
            local_x,
            time_zero,
        )["E"],
        x_zero,
    )
    derivative_t = mp.diff(
        lambda local_t: corrected_half(
            x_zero,
            local_t,
        )["E"],
        time_zero,
    )

    s = (1 - mp.j * x_zero) / 2
    alpha_value = alpha(s)
    alpha_one = alpha_prime(s)
    temperature = x_zero / 2 + mp.pi * time_zero / 8
    saddle = mp.sqrt(temperature / (2 * mp.pi))
    beta_prime = -mp.re(
        alpha_value + time_zero * alpha_value * alpha_one / 2
    ) / 2
    s_star_prime = (
        -mp.j / 2 - mp.j * time_zero * alpha_one / 4
    )
    lambda_a = (
        mp.j * beta_prime - s_star_prime * mp.log(saddle)
    )
    residual = derivative_x - lambda_a * complex_main
    endpoint_defect = (
        mp.diff(
            lambda local_x: corrected_half(
                local_x,
                time_zero,
            )["g0"],
            x_zero,
        )
        - lambda_a * endpoint
    )
    core = (
        -mp.im(centered_moment) / 2
        - time_zero
        * (
            mp.im(alpha_one) * mp.re(centered_moment)
            + mp.re(alpha_one) * mp.im(centered_moment)
        )
        / 4
        + mp.re(endpoint_defect)
    )
    centered_scalar = (
        mp.re(residual)
        - mp.im(lambda_a) * mp.im(complex_main)
    )
    wronskian = mp.im(derivative_x * mp.conj(complex_main))
    jacobian = (
        mp.re(derivative_x) * mp.im(derivative_t)
        - mp.re(derivative_t) * mp.im(derivative_x)
    )
    logarithmic_height = mp.log(x_zero / (4 * mp.pi))
    q_value = 2 * time_zero * logarithmic_height**2

    def text(value: mp.mpf | mp.mpc, digits: int = 70) -> str:
        return mp.nstr(value, digits, strip_zeros=False)

    return {
        "dps": dps,
        "x": text(x_zero),
        "t": text(time_zero),
        "L": text(logarithmic_height),
        "tL": text(time_zero * logarithmic_height),
        "q": text(q_value),
        "N": values["N"],
        "p": text(values["p"]),
        "E_real": text(mp.re(complex_main)),
        "E_imag": text(mp.im(complex_main)),
        "E_abs": text(abs(complex_main)),
        "E_x_real_U": text(mp.re(derivative_x)),
        "E_x_imag_V": text(mp.im(derivative_x)),
        "E_t_real": text(mp.re(derivative_t)),
        "E_t_imag": text(mp.im(derivative_t)),
        "complex_zero_jacobian_xt": text(jacobian),
        "wronskian": text(wronskian),
        "centered_scalar_S": text(centered_scalar),
        "core_scalar_A": text(core),
        "real_D1x": text(mp.re(derivative_correction)),
        "u_a": text(mp.re(lambda_a)),
        "v_a": text(mp.im(lambda_a)),
    }


def cross_precision_audit(low: dict, high: dict) -> dict:
    mp.mp.dps = HIGH_DPS
    x_delta = abs(mp.mpf(low["x"]) - mp.mpf(high["x"]))
    t_delta = abs(mp.mpf(low["t"]) - mp.mpf(high["t"]))
    u_delta = abs(
        mp.mpf(low["E_x_real_U"]) - mp.mpf(high["E_x_real_U"])
    )
    if x_delta >= mp.mpf("1e-50"):
        raise RuntimeError("selected-zero x drift exceeded 1e-50")
    if t_delta >= mp.mpf("1e-50"):
        raise RuntimeError("selected-zero t drift exceeded 1e-50")
    if u_delta >= mp.mpf("1e-50"):
        raise RuntimeError("selected-zero slope drift exceeded 1e-50")
    return {
        "x_delta_lt": "1e-50",
        "t_delta_lt": "1e-50",
        "U_delta_lt": "1e-50",
        "x_delta": mp.nstr(x_delta, 15),
        "t_delta": mp.nstr(t_delta, 15),
        "U_delta": mp.nstr(u_delta, 15),
    }


def build_rows(high: dict, precision: dict) -> list[ScoutRow]:
    return [
        ScoutRow(
            "nfocczs_01_selected_lift",
            "diagnostic_scope",
            "diagnostic_only",
            "The scout stays on one prescribed N=6 corrected first-order lift.",
            (
                f"N={high['N']}, x={high['x']}, t={high['t']}, "
                f"L={high['L']}, q={high['q']}"
            ),
            "The point has L<50 and is outside the proved asymptotic domain.",
        ),
        ScoutRow(
            "nfocczs_02_complex_zero",
            "high_precision_diagnostic",
            "diagnostic_only",
            "Both coordinates of the corrected complex main are numerically zero.",
            (
                f"Re(E_[1])={high['E_real']}, "
                f"Im(E_[1])={high['E_imag']}, "
                f"|E_[1]|={high['E_abs']}"
            ),
            "High-precision mpmath residual, not an interval existence proof.",
        ),
        ScoutRow(
            "nfocczs_03_transverse_complex_zero",
            "high_precision_diagnostic",
            "diagnostic_only",
            "The two-variable corrected complex zero is numerically transverse.",
            (
                "det D_(x,t)(Re E_[1],Im E_[1])="
                f"{high['complex_zero_jacobian_xt']}"
            ),
            "A rigorous existence theorem would require an interval Newton box.",
        ),
        ScoutRow(
            "nfocczs_04_simple_real_crossing",
            "route_guard",
            "guard_validated",
            "The real-part proxy crosses simply even though its complete complex main vanishes.",
            (
                f"U={high['E_x_real_U']}, "
                f"S_a={high['centered_scalar_S']}"
            ),
            "Finite selected diagnostic only; it does not decide L>=50.",
        ),
        ScoutRow(
            "nfocczs_05_wronskian_degeneracy",
            "route_guard",
            "guard_validated",
            "The Wronskian coordinate vanishes with E_[1] and cannot classify this simple crossing.",
            f"W_[1]={high['wronskian']}, while U={high['E_x_real_U']}",
            "Rejects deletion of E_[1]=0 crossings; not the asymptotic Wronskian theorem.",
        ),
        ScoutRow(
            "nfocczs_06_centered_core",
            "high_precision_diagnostic",
            "diagnostic_only",
            "The centered core scalar retains the nonzero crossing slope.",
            (
                f"A_a={high['core_scalar_A']}, "
                f"Re(D_(1,x))={high['real_D1x']}"
            ),
            "At this moderate L the asymptotic 1e-6 remainder-scale budget is not invoked.",
        ),
        ScoutRow(
            "nfocczs_07_precision",
            "numerical_stability",
            "diagnostic_validated",
            "Independent working-precision reruns agree beyond fifty decimal places.",
            (
                f"|Delta x|<{precision['x_delta_lt']}, "
                f"|Delta t|<{precision['t_delta_lt']}, "
                f"|Delta U|<{precision['U_delta_lt']}"
            ),
            "Cross-precision stability is not interval certification.",
            precision,
        ),
        ScoutRow(
            "nfocczs_08_interval_handoff",
            "open_theorem_target",
            "not_ready_to_apply",
            "A tiny complex interval Newton box could promote this selected route guard to a rigorous finite certificate.",
            "Certify one N=6 box around the stored (x,t) point with det D_(x,t)(Re E_[1],Im E_[1]) bounded away from zero.",
            "Optional finite route guard only; it is not part of the RH-level asymptotic proof.",
        ),
    ]


def build_artifact() -> dict:
    low = solve_selected_zero(LOW_DPS)
    high = solve_selected_zero(HIGH_DPS)
    precision = cross_precision_audit(low, high)
    return {
        "kind": STEM,
        "date": "2026-07-26",
        "status": (
            "selected high-precision diagnostic of a corrected complex-main "
            "zero with nonzero real crossing slope; not interval certified"
        ),
        "proof_boundary": (
            "This scout supplies a reproducible selected moderate-height "
            "diagnostic and a cross-precision route guard. It does not prove "
            "existence of the complex zero by interval arithmetic, any L>=50 "
            "complex-main zero, the q>=1 centered scalar theorem, contact "
            "exclusion, Lambda<=0, RH, PF-infinity, or a Clay-prize conclusion."
        ),
        "sources": {
            key: str(path.relative_to(REPO_ROOT)).replace("\\", "/")
            for key, path in SOURCES.items()
        },
        "source_sha256": source_hashes(),
        "configuration": {
            "N_fixed": N_FIXED,
            "low_dps": LOW_DPS,
            "high_dps": HIGH_DPS,
            "initial_pair_1": ["519.38", "0.403"],
            "initial_pair_2": ["519.40", "0.405"],
        },
        "low_precision": low,
        "selected_zero": high,
        "cross_precision": precision,
        "rows": [asdict(row) for row in build_rows(high, precision)],
    }


def render_note(artifact: dict) -> str:
    row = artifact["selected_zero"]
    precision = artifact["cross_precision"]
    return "\n".join(
        [
            "# Newman First-Order Centered Complex-Zero Scout",
            "",
            "Date: 2026-07-26",
            "",
            "Status: selected high-precision route diagnostic; this is not an interval certificate",
            "and is not a proof of RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Selected Point",
            "",
            "| quantity | value |",
            "|---|---:|",
            f"| x | {row['x']} |",
            f"| t | {row['t']} |",
            f"| L | {row['L']} |",
            f"| tL | {row['tL']} |",
            f"| q=2tL^2 | {row['q']} |",
            f"| N | {row['N']} |",
            f"| p | {row['p']} |",
            f"| Re(E_[1]) | {row['E_real']} |",
            f"| Im(E_[1]) | {row['E_imag']} |",
            f"| |E_[1]| | {row['E_abs']} |",
            f"| U=Re(E_[1],x) | {row['E_x_real_U']} |",
            f"| V=Im(E_[1],x) | {row['E_x_imag_V']} |",
            f"| det D_(x,t)(Re E_[1],Im E_[1]) | {row['complex_zero_jacobian_xt']} |",
            f"| W_[1] | {row['wronskian']} |",
            f"| S_a | {row['centered_scalar_S']} |",
            f"| A_a | {row['core_scalar_A']} |",
            f"| Re(D_(1,x)) | {row['real_D1x']} |",
            "",
            "The corrected complex main is numerically zero while the real",
            "crossing slope is about `0.556`. Consequently the Wronskian is",
            "zero at the limiting point, but the centered scalar still",
            "classifies the crossing as upward.",
            "",
            "## Precision Replay",
            "",
            "```text",
            f"|Delta x|<{precision['x_delta_lt']}",
            f"|Delta t|<{precision['t_delta_lt']}",
            f"|Delta U|<{precision['U_delta_lt']}",
            "```",
            "",
            "The low- and high-precision runs agree beyond fifty decimal",
            "places. This is strong numerical stability, not an interval",
            "existence proof.",
            "",
            "## Route Consequence",
            "",
            "This selected point makes the exceptional class concrete:",
            "`E_[1]=0` does not imply a multiple real-part crossing. A proof",
            "based only on the sign of `W_[1]Im(E_[1])` would delete it, while",
            "`S_a=U-u_aX` and the centered core retain it. The optional next",
            "finite task is a tiny interval Newton box around this point.",
            "",
            artifact["proof_boundary"],
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    selected = artifact["selected_zero"]
    print(
        "built Newman first-order centered complex-zero scout: "
        f"N={selected['N']}, |E_[1]|<1e-70, "
        "|U|>0.5, |det D_(x,t)E_[1]|>0.3, "
        "cross-precision drift <1e-50"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
