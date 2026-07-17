#!/usr/bin/env python3
"""High-precision scout for the corrected critical finite coercivity."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_scaled_coercivity_scout.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_critical_scaled_coercivity_scout.md"
)
X_VALUES = (1000, 10_000, 100_000, 1_000_000)
C_VALUES = ("0", "0.5", "1", "2", "3", "4")


def alpha(s: mp.mpc) -> mp.mpc:
    return 1 / (2 * s) + 1 / (s - 1) + mp.log(s / (2 * mp.pi)) / 2


def alpha_prime(s: mp.mpc) -> mp.mpc:
    return -1 / (2 * s**2) - 1 / (s - 1) ** 2 + 1 / (2 * s)


def alpha_second(s: mp.mpc) -> mp.mpc:
    return 1 / s**3 + 2 / (s - 1) ** 3 - 1 / (2 * s**2)


def m_zero(s: mp.mpc) -> mp.mpc:
    return (
        mp.mpf(1) / 8
        * (s * (s - 1) / 2)
        * mp.pi ** (-s / 2)
        * mp.sqrt(2 * mp.pi)
        * mp.exp((s / 2 - mp.mpf("0.5")) * mp.log(s / 2) - s / 2)
    )


def m_time(s: mp.mpc, t: mp.mpf) -> mp.mpc:
    a = alpha(s)
    return mp.exp(t * a**2 / 4) * m_zero(s)


def c_zero(p: mp.mpf) -> mp.mpc:
    return (
        mp.exp(mp.pi * mp.j * (p**2 / 2 + mp.mpf(3) / 8))
        - mp.j * mp.sqrt(2) * mp.cos(mp.pi * p / 2)
    ) / (2 * mp.cos(mp.pi * p))


def xi(s: mp.mpc) -> mp.mpc:
    return (
        mp.mpf("0.5")
        * s
        * (s - 1)
        * mp.pi ** (-s / 2)
        * mp.gamma(s / 2)
        * mp.zeta(s)
    )


def compute_point(x_int: int, c_text: str, dps: int) -> dict:
    mp.mp.dps = dps
    x = mp.mpf(x_int)
    c = mp.mpf(c_text)
    L = mp.log(x / (4 * mp.pi))
    t = c / L if c else mp.mpf(0)
    t_prime = x / 2 + mp.pi * t / 8
    saddle = mp.sqrt(t_prime / (2 * mp.pi))
    n_cutoff = int(mp.floor(saddle))
    p = 1 - 2 * (saddle - n_cutoff)

    s = (1 - mp.j * x) / 2
    a = alpha(s)
    a1 = alpha_prime(s)
    a2 = alpha_second(s)
    s_star = s + t * a / 2
    s_star_first = -mp.j / 2 - mp.j * t * a1 / 4
    s_star_second = -t * a2 / 8
    log_m_first = a + t * a * a1 / 2
    log_m_second = a1 + t * (a1**2 + a * a2) / 2
    beta_first = -mp.re(log_m_first) / 2
    beta_second = -mp.im(log_m_second) / 4
    potential = mp.re(log_m_second) / 4

    d0 = mp.mpc(0)
    d1 = mp.mpc(0)
    d2 = mp.mpc(0)
    for n in range(1, n_cutoff + 1):
        log_n = mp.log(n)
        term = mp.exp(t * log_n**2 / 4 - s_star * log_n)
        d0 += term
        d1 += -s_star_first * log_n * term
        d2 += (
            s_star_first**2 * log_n**2 - s_star_second * log_n
        ) * term

    normalizer = m_time(s, t)
    amplitude = abs(normalizer)
    phase = normalizer / amplitude
    p0 = 2 * mp.re(phase * d0)
    p1 = 2 * mp.re(phase * (d1 + mp.j * beta_first * d0))
    p2 = 2 * mp.re(
        phase
        * (
            d2
            + 2 * mp.j * beta_first * d1
            + (mp.j * beta_second - beta_first**2) * d0
        )
    )
    uncorrected = p1**2 - p0 * p2 + potential * p0**2

    def normalized_correction(xx: mp.mpf) -> mp.mpf:
        ss = (1 - mp.j * xx) / 2
        local_amplitude = abs(m_time(ss, t))
        local_t_prime = xx / 2 + mp.pi * t / 8
        local_saddle = mp.sqrt(local_t_prime / (2 * mp.pi))
        local_p = 1 - 2 * (local_saddle - n_cutoff)
        u_rs = mp.exp(
            -mp.j
            * (
                (local_t_prime / 2)
                * mp.log(local_t_prime / (2 * mp.pi))
                - local_t_prime / 2
                - mp.pi / 8
            )
        )
        raw = (
            2
            * ((-1) ** n_cutoff)
            * mp.exp(t * mp.pi**2 / 64)
            * mp.re(
                m_zero(mp.j * local_t_prime)
                * c_zero(local_p)
                * u_rs
                * mp.exp(mp.j * mp.pi / 8)
            )
        )
        return raw / local_amplitude

    q0 = normalized_correction(x)
    q1 = mp.diff(normalized_correction, x, 1)
    q2 = mp.diff(normalized_correction, x, 2)
    j0 = p0 - q0
    j1 = p1 - q1
    j2 = p2 - q2
    corrected = j1**2 - j0 * j2 + potential * j0**2

    row = {
        "x": x_int,
        "c": c_text,
        "L": mp.nstr(L, 35, strip_zeros=False),
        "t": mp.nstr(t, 35, strip_zeros=False),
        "N": n_cutoff,
        "p": mp.nstr(p, 35, strip_zeros=False),
        "uncorrected_curvature": mp.nstr(
            uncorrected, 35, strip_zeros=False
        ),
        "normalized_correction": mp.nstr(q0, 35, strip_zeros=False),
        "corrected_curvature": mp.nstr(corrected, 35, strip_zeros=False),
        "corrected_main_jet": [
            mp.nstr(value, 35, strip_zeros=False)
            for value in (j0, j1, j2)
        ],
    }
    if c == 0:
        def exact_normalized(xx: mp.mpf) -> mp.mpf:
            return mp.re(xi((1 + mp.j * xx) / 2) / 8) / abs(
                m_time((1 - mp.j * xx) / 2, t)
            )

        z0 = exact_normalized(x)
        z1 = mp.diff(exact_normalized, x, 1)
        z2 = mp.diff(exact_normalized, x, 2)
        exact_curvature = z1**2 - z0 * z2 + potential * z0**2
        row["exact_xi_curvature"] = mp.nstr(
            exact_curvature, 35, strip_zeros=False
        )
        row["relative_corrected_to_exact_delta"] = mp.nstr(
            abs(corrected - exact_curvature) / abs(exact_curvature),
            25,
            strip_zeros=False,
        )
    return row


def compute_rows(dps: int) -> list[dict]:
    return [
        compute_point(x, c, dps)
        for x in X_VALUES
        for c in C_VALUES
    ]


def compare_rows(coarse: list[dict], fine: list[dict]) -> dict:
    maximum = mp.mpf(0)
    stable_signs = True
    for left, right in zip(coarse, fine, strict=True):
        if (left["x"], left["c"]) != (right["x"], right["c"]):
            raise RuntimeError("critical scout coarse/fine key mismatch")
        lv = mp.mpf(left["corrected_curvature"])
        rv = mp.mpf(right["corrected_curvature"])
        stable_signs &= (lv > 0) == (rv > 0)
        maximum = max(maximum, abs(lv - rv) / abs(rv))
    return {
        "stable_corrected_signs": bool(stable_signs),
        "max_relative_corrected_curvature_delta": mp.nstr(maximum, 25),
    }


def build_artifact(coarse_dps: int, fine_dps: int) -> dict:
    coarse = compute_rows(coarse_dps)
    fine = compute_rows(fine_dps)
    convergence = compare_rows(coarse, fine)
    if not convergence["stable_corrected_signs"]:
        raise RuntimeError("corrected critical signs are precision-unstable")
    if mp.mpf(convergence["max_relative_corrected_curvature_delta"]) >= mp.mpf("1e-30"):
        raise RuntimeError("corrected critical curvature precision failed")
    if any(mp.mpf(row["corrected_curvature"]) <= 0 for row in fine):
        raise RuntimeError("sampled corrected critical curvature is nonpositive")
    exact_rows = [row for row in fine if row["c"] == "0"]
    if any(
        mp.mpf(row["relative_corrected_to_exact_delta"]) >= mp.mpf("1e-3")
        for row in exact_rows
    ):
        raise RuntimeError("corrected t=0 xi cross-check failed")
    return {
        "kind": "jensen_window_pf_newman_polymath15_critical_scaled_coercivity_scout",
        "date": "2026-07-17",
        "status": (
            "high-precision corrected critical-coercivity diagnostics with exact "
            "xi cross-checks at t=0; not a proof of a uniform sign, Lambda<=0, or RH"
        ),
        "proof_boundary": (
            "These are 24 point diagnostics with independent precision reruns. "
            "They do not certify intervals, all scaled times, all frequencies, "
            "the corrected coercivity target, Lambda<=0, or RH."
        ),
        "parameters": {
            "coarse_dps": coarse_dps,
            "fine_dps": fine_dps,
            "x_values": list(X_VALUES),
            "c_values": list(C_VALUES),
        },
        "convergence": convergence,
        "rows": fine,
        "decision": (
            "The endpoint correction materially improves some rows and tracks the "
            "exact t=0 xi curvature. Continue with a corrected-main interval or "
            "asymptotic coercivity theorem; do not promote point positivity."
        ),
        "sources": [
            "outputs/jensen_window_pf_newman_polymath15_critical_scaled_coercivity_target.md",
            "https://arxiv.org/abs/1904.12438",
        ],
    }


def render_note(artifact: dict) -> str:
    table = [
        "| x | c=tL | N | p | C[P] | Q=Ct/A | C[P-Q] | exact xi C (c=0) |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in artifact["rows"]:
        table.append(
            "| {x} | {c} | {n} | {p} | {cp} | {q} | {cj} | {exact} |".format(
                x=row["x"],
                c=row["c"],
                n=row["N"],
                p=mp.nstr(mp.mpf(row["p"]), 6),
                cp=mp.nstr(mp.mpf(row["uncorrected_curvature"]), 8),
                q=mp.nstr(mp.mpf(row["normalized_correction"]), 6),
                cj=mp.nstr(mp.mpf(row["corrected_curvature"]), 8),
                exact=(
                    mp.nstr(mp.mpf(row["exact_xi_curvature"]), 8)
                    if "exact_xi_curvature" in row
                    else "-"
                ),
            )
        )
    exact_deltas = [
        mp.mpf(row["relative_corrected_to_exact_delta"])
        for row in artifact["rows"]
        if "relative_corrected_to_exact_delta" in row
    ]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Polymath-15 Critical Scaled Coercivity Scout",
            "",
            "Date: 2026-07-17",
            "",
            "Status: high-precision corrected-main diagnostics. This is not a proof",
            "of a uniform sign, `Lambda <= 0`, or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_scaled_coercivity_scout.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_scaled_coercivity_scout.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_scaled_coercivity_scout.py",
            "```",
            "",
            "The table evaluates the uncorrected finite curvature, the explicit",
            "Polymath-15 endpoint correction, and the corrected curvature for six",
            "scaled times at four frequencies:",
            "",
            *table,
            "",
            "All 24 corrected point values are positive. Coarse/fine precision",
            "reruns have maximum relative corrected-curvature delta",
            f"`{artifact['convergence']['max_relative_corrected_curvature_delta']}`.",
            "",
            "At `c=0`, the corrected finite expression is independently compared",
            "with `H_0(x)=xi((1+i*x)/2)/8`, including derivatives of the real",
            "normalization. The largest relative curvature discrepancy is",
            f"`{mp.nstr(max(exact_deltas), 10)}` and the smallest is",
            f"`{mp.nstr(min(exact_deltas), 10)}`.",
            "",
            "The correction is not cosmetic: it changes the curvature visibly and",
            "systematically improves agreement with exact xi. The sampled margins",
            "are substantial, but point positivity is not an interval or asymptotic",
            "theorem. The live obligation remains a uniform corrected coercivity",
            "bound throughout the critical scaled layer.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--coarse-dps", type=int, default=50)
    parser.add_argument("--fine-dps", type=int, default=70)
    args = parser.parse_args()
    artifact = build_artifact(args.coarse_dps, args.fine_dps)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman Polymath-15 critical scaled coercivity scout: "
        f"{len(artifact['rows'])} rows, 24 positive diagnostics, 4 exact xi checks"
    )


if __name__ == "__main__":
    main()
