#!/usr/bin/env python3
"""High-precision scout for cancellation in finite modular-blend sums."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_theta_modular_blend_high_frequency_scout.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_theta_modular_blend_high_frequency_scout.md"
)
DEFAULT_X = (80, 100, 120, 150, 200)
DEFAULT_TIMES = (0, mp.mpf("0.5"))
PANELS = (
    ("0", ".04"),
    (".04", ".08"),
    (".08", ".14"),
    (".14", ".22"),
    (".22", ".34"),
    (".34", ".5"),
    (".5", ".72"),
    (".72", "1"),
    ("1", "1.4"),
    ("1.4", "2"),
    ("2", "2.8"),
)


def theta_summand(u: mp.mpf, n: int) -> mp.mpf:
    return (
        2 * mp.pi**2 * n**4 * mp.exp(9 * u)
        - 3 * mp.pi * n**2 * mp.exp(5 * u)
    ) * mp.exp(-mp.pi * n**2 * mp.exp(4 * u))


def blend_weight(u: mp.mpf) -> mp.mpf:
    return (1 + mp.erf(3 * mp.sinh(4 * u))) / 2


def xi(s: mp.mpc) -> mp.mpc:
    return (
        mp.mpf("0.5")
        * s
        * (s - 1)
        * mp.power(mp.pi, -s / 2)
        * mp.gamma(s / 2)
        * mp.zeta(s)
    )


def xi_jets(x: mp.mpf) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    s = mp.mpf("0.5") + mp.j * x / 2
    return (
        mp.re(xi(s) / 8),
        mp.re(mp.j * mp.diff(xi, s, 1) / 16),
        mp.re(-mp.diff(xi, s, 2) / 32),
    )


def quadrature_rows(node_count: int, theta_terms: int) -> list[tuple]:
    nodes, weights = mp.gauss_quadrature(node_count, "legendre")
    rows: list[tuple] = []
    for left_text, right_text in PANELS:
        left = mp.mpf(left_text)
        right = mp.mpf(right_text)
        for node, weight in zip(nodes, weights):
            u = (right - left) * node / 2 + (right + left) / 2
            qweight = (right - left) * weight / 2
            omega = blend_weight(u)
            forward = [theta_summand(u, n) for n in range(1, theta_terms + 1)]
            reflected = [theta_summand(-u, n) for n in range(1, 4)]
            partial = mp.fsum(
                omega * forward[n - 1] + (1 - omega) * reflected[n - 1]
                for n in range(1, 4)
            )
            full = mp.fsum(forward)
            rows.append((u, qweight, partial, full))
    return rows


def transform_jet(
    rows: list[tuple], t: mp.mpf, x: mp.mpf, column: int
) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    weighted = [
        (u, qweight * mp.exp(t * u * u) * row[column])
        for row in rows
        for u, qweight in [(row[0], row[1])]
    ]
    return (
        mp.fsum(weight * mp.cos(x * u) for u, weight in weighted),
        mp.fsum(-weight * u * mp.sin(x * u) for u, weight in weighted),
        mp.fsum(-weight * u * u * mp.cos(x * u) for u, weight in weighted),
    )


def laguerre(jet: tuple[mp.mpf, mp.mpf, mp.mpf]) -> mp.mpf:
    value, first, second = jet
    return first * first - value * second


def format_mpf(value: mp.mpf, digits: int) -> str:
    return mp.nstr(value, digits, strip_zeros=False)


def compute_rows(
    node_count: int,
    dps: int,
    theta_terms: int,
    x_values: tuple[int, ...] = DEFAULT_X,
) -> list[dict]:
    mp.mp.dps = dps
    qrows = quadrature_rows(node_count, theta_terms)
    output: list[dict] = []
    for t in DEFAULT_TIMES:
        t = mp.mpf(t)
        for x_int in x_values:
            x = mp.mpf(x_int)
            partial_jet = transform_jet(qrows, t, x, 2)
            full_jet = transform_jet(qrows, t, x, 3)
            partial_laguerre = laguerre(partial_jet)
            full_laguerre = laguerre(full_jet)
            ratio = full_laguerre / partial_laguerre
            row = {
                "t": format_mpf(t, 8),
                "x": x_int,
                "partial_terms": 3,
                "theta_terms": theta_terms,
                "partial_jet": [format_mpf(value, dps - 8) for value in partial_jet],
                "full_jet": [format_mpf(value, dps - 8) for value in full_jet],
                "partial_laguerre": format_mpf(partial_laguerre, dps - 8),
                "full_laguerre": format_mpf(full_laguerre, dps - 8),
                "tail_correction": format_mpf(
                    full_laguerre - partial_laguerre, dps - 8
                ),
                "full_to_partial_ratio": format_mpf(ratio, dps - 8),
                "cancellation_digits": format_mpf(-mp.log10(abs(ratio)), 20),
            }
            if t == 0:
                exact_jet = xi_jets(x)
                exact_laguerre = laguerre(exact_jet)
                jet_relative = max(
                    abs(left - right) / max(abs(right), mp.mpf(10) ** (-(dps - 10)))
                    for left, right in zip(full_jet, exact_jet, strict=True)
                )
                laguerre_relative = abs(full_laguerre - exact_laguerre) / abs(
                    exact_laguerre
                )
                row["exact_xi_jet"] = [
                    format_mpf(value, dps - 8) for value in exact_jet
                ]
                row["exact_xi_laguerre"] = format_mpf(
                    exact_laguerre, dps - 8
                )
                row["max_relative_xi_jet_delta"] = format_mpf(
                    jet_relative, 20
                )
                row["relative_xi_laguerre_delta"] = format_mpf(
                    laguerre_relative, 20
                )
            output.append(row)
    return output


def compare_rows(coarse: list[dict], fine: list[dict]) -> dict:
    max_partial_delta = mp.mpf(0)
    max_full_delta = mp.mpf(0)
    for left, right in zip(coarse, fine, strict=True):
        if (left["t"], left["x"]) != (right["t"], right["x"]):
            raise RuntimeError("coarse/fine high-frequency row mismatch")
        for field, accumulator in (
            ("partial_laguerre", "partial"),
            ("full_laguerre", "full"),
        ):
            left_value = mp.mpf(left[field])
            right_value = mp.mpf(right[field])
            delta = abs(left_value - right_value) / abs(right_value)
            if accumulator == "partial":
                max_partial_delta = max(max_partial_delta, delta)
            else:
                max_full_delta = max(max_full_delta, delta)
    return {
        "max_relative_partial_laguerre_delta": format_mpf(max_partial_delta, 20),
        "max_relative_full_laguerre_delta": format_mpf(max_full_delta, 20),
    }


def build_artifact(
    dps: int, coarse_nodes: int, fine_nodes: int, theta_terms: int
) -> dict:
    coarse = compute_rows(coarse_nodes, dps, theta_terms)
    fine = compute_rows(fine_nodes, dps, theta_terms)
    convergence = compare_rows(coarse, fine)
    if mp.mpf(convergence["max_relative_partial_laguerre_delta"]) >= mp.mpf("1e-30"):
        raise RuntimeError("partial Laguerre quadrature convergence failed")
    if mp.mpf(convergence["max_relative_full_laguerre_delta"]) >= mp.mpf("1e-25"):
        raise RuntimeError("full Laguerre quadrature convergence failed")
    if any(mp.mpf(row["partial_laguerre"]) <= 0 for row in fine):
        raise RuntimeError("partial Laguerre scout positivity failed")
    if any(mp.mpf(row["full_laguerre"]) <= 0 for row in fine):
        raise RuntimeError("full Laguerre scout positivity failed")
    return {
        "kind": "jensen_window_pf_newman_theta_modular_blend_high_frequency_scout",
        "date": "2026-07-16",
        "status": (
            "high-precision high-frequency cancellation diagnostic; not a proof "
            "of a global Laguerre inequality, Lambda<=0, or RH"
        ),
        "proof_boundary": (
            "The rows are high-precision point diagnostics with independent "
            "coarse/fine quadrature and an exact xi cross-check at t=0. They show "
            "that a fixed three-block Laguerre expression is almost completely "
            "cancelled by the arithmetic tail at high frequency. They do not "
            "certify an interval, all frequencies, all times, Lambda<=0, or RH."
        ),
        "parameters": {
            "dps": dps,
            "coarse_nodes_per_panel": coarse_nodes,
            "fine_nodes_per_panel": fine_nodes,
            "panels": [list(panel) for panel in PANELS],
            "theta_terms": theta_terms,
            "partial_terms": 3,
            "x_values": list(DEFAULT_X),
            "times": ["0", "0.5"],
        },
        "convergence": convergence,
        "rows": fine,
        "decision": (
            "Do not promote a fixed finite modular-blend sum by absolute tail "
            "dominance. Preserve the coupled infinite cancellation or let the "
            "arithmetic truncation grow with frequency."
        ),
        "sources": [
            "outputs/jensen_window_pf_newman_theta_modular_blend_gate.md",
            "outputs/jensen_window_pf_newman_strict_laguerre_correlation_target.md",
        ],
    }


def render_note(artifact: dict) -> str:
    rows = artifact["rows"]
    table = [
        "| t | x | L[three blocks] | L[full] | full/partial | cancelled digits |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        table.append(
            "| {t} | {x} | {partial} | {full} | {ratio} | {digits} |".format(
                t=row["t"],
                x=row["x"],
                partial=mp.nstr(mp.mpf(row["partial_laguerre"]), 12),
                full=mp.nstr(mp.mpf(row["full_laguerre"]), 12),
                ratio=mp.nstr(mp.mpf(row["full_to_partial_ratio"]), 10),
                digits=mp.nstr(mp.mpf(row["cancellation_digits"]), 8),
            )
        )
    return "\n".join(
        [
            "# Jensen-Window PF Newman Modular-Blend High-Frequency Scout",
            "",
            "Date: 2026-07-16",
            "",
            "Status: high-precision cancellation diagnostic. This is not a proof",
            "of a global Laguerre inequality, `Lambda <= 0`, or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_theta_modular_blend_high_frequency_scout.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_theta_modular_blend_high_frequency_scout.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_modular_blend_high_frequency_scout.py",
            "```",
            "",
            "The entire modular blend makes every finite block smooth and",
            "positive, but that does not make a fixed finite sum dominant at high",
            "frequency. The first three blended blocks and the full theta kernel give",
            "",
            *table,
            "",
            "At `t=0`, the full quadrature jets are independently checked against",
            "`H_0(x)=xi((1+i*x)/2)/8` and its first two derivatives. Coarse/fine",
            "Gauss-Legendre quadrature has maximum relative deltas",
            "",
            "```text",
            f"partial L: {artifact['convergence']['max_relative_partial_laguerre_delta']}",
            f"full L:    {artifact['convergence']['max_relative_full_laguerre_delta']}",
            "```",
            "",
            "By `x=150`, about eleven to twelve decimal digits of the positive",
            "three-block Laguerre expression are cancelled. By `x=200`, the",
            "cancellation is about twenty-one digits. The full sampled expression",
            "remains positive, but its sign is the small coupled remainder, not a",
            "finite-base margin.",
            "",
            "Decision: a fixed finite modular-blend truncation plus an absolute",
            "tail estimate is not the global architecture. A viable continuation",
            "must preserve the infinite arithmetic cancellation, or use a",
            "frequency-adaptive arithmetic truncation with a sign-aware remainder.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--dps", type=int, default=75)
    parser.add_argument("--coarse-nodes", type=int, default=120)
    parser.add_argument("--fine-nodes", type=int, default=160)
    parser.add_argument("--theta-terms", type=int, default=13)
    args = parser.parse_args()
    artifact = build_artifact(
        args.dps, args.coarse_nodes, args.fine_nodes, args.theta_terms
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman modular-blend high-frequency scout: "
        f"{len(artifact['rows'])} rows, "
        f"max cancellation {max(float(row['cancellation_digits']) for row in artifact['rows']):.2f} digits"
    )


if __name__ == "__main__":
    main()
