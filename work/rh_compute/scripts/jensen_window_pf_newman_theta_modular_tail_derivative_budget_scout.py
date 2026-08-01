#!/usr/bin/env python3
"""Scout explicit modular-tail derivative budgets and direct-C1 ratios."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from math import comb, factorial
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import erf
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_theta_modular_blend_adaptive_saddle_gate as saddle  # noqa: E402


STEM = "jensen_window_pf_newman_theta_modular_tail_derivative_budget_scout"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
ENVELOPE_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_modular_tail_derivative_envelope_gate.json"
)
T_CAP = 0.2
ORDERS = tuple(range(5, 10))
N_VALUES = tuple(range(4, 11))
X_VALUES = (80, 150, 200, 240, 260, 280, 300, 312)
TIMES = ("0", "0.2")
KAPPAS = (2, 3, 4, 5)
DERIVATIVE_PANELS = (
    (0.0, 0.02),
    (0.02, 0.05),
    (0.05, 0.1),
    (0.1, 0.2),
    (0.2, 0.35),
    (0.35, 0.5),
    (0.5, 0.7),
    (0.7, 0.9),
    (0.9, 1.15),
    (1.15, 1.45),
    (1.45, 1.8),
    (1.8, 2.2),
)


@dataclass(frozen=True)
class ScoutRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def quadrature(
    node_count: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nodes, weights = leggauss(node_count)
    u_rows: list[np.ndarray] = []
    weight_rows: list[np.ndarray] = []
    panel_rows: list[np.ndarray] = []
    for panel_index, (left, right) in enumerate(DERIVATIVE_PANELS):
        u_rows.append((right - left) * nodes / 2 + (right + left) / 2)
        weight_rows.append((right - left) * weights / 2)
        panel_rows.append(np.full(node_count, panel_index, dtype=int))
    return (
        np.concatenate(u_rows),
        np.concatenate(weight_rows),
        np.concatenate(panel_rows),
    )


def block_derivative_functions(max_order: int) -> list:
    u, n = sp.symbols("u n", positive=True, real=True)
    pi = sp.pi

    def phi(z: sp.Expr) -> sp.Expr:
        return (
            2 * pi**2 * n**4 * sp.exp(9 * z)
            - 3 * pi * n**2 * sp.exp(5 * z)
        ) * sp.exp(-pi * n**2 * sp.exp(4 * z))

    omega = (1 + sp.erf(3 * sp.sinh(4 * u))) / 2
    block = omega * phi(u) + (1 - omega) * phi(-u)
    derivatives: list[sp.Expr] = []
    current = block
    for _ in range(max_order + 1):
        derivatives.append(current)
        current = sp.diff(current, u)
    return [
        sp.lambdify((u, n), expr, modules=["numpy", {"erf": erf}])
        for expr in derivatives
    ]


def heat_polynomial(order: int, u: np.ndarray) -> np.ndarray:
    value = np.zeros_like(u)
    for ell in range(order // 2 + 1):
        value += (
            factorial(order)
            * (2 * u) ** (order - 2 * ell)
            * T_CAP ** (order - ell)
            / (factorial(ell) * factorial(order - 2 * ell))
        )
    return value


def evaluate_derivative_table(
    functions: list,
    node_count: int,
    n_cap: int,
) -> tuple[list[dict], dict]:
    u, weights, panel_index = quadrature(node_count)
    values = np.empty((len(functions), n_cap, len(u)))
    for derivative_order, function in enumerate(functions):
        for n in range(1, n_cap + 1):
            values[derivative_order, n - 1] = function(u, n)
    if not np.isfinite(values).all():
        raise RuntimeError("nonfinite blended-block derivative sample")

    polynomials = [heat_polynomial(order, u) for order in range(len(functions))]
    heat = np.exp(T_CAP * u * u)
    rows: list[dict] = []
    max_last_panel_fraction = 0.0
    for retained in N_VALUES:
        derivatives = [
            np.sum(values[order, retained:n_cap], axis=0)
            for order in range(len(functions))
        ]
        moment0_integrand = heat * np.abs(derivatives[0])
        moment1_integrand = heat * u * np.abs(derivatives[0])
        moment0 = float(np.dot(weights, moment0_integrand))
        moment1 = float(np.dot(weights, moment1_integrand))
        for order in ORDERS:
            raw_value = sum(
                comb(order, b)
                * polynomials[b]
                * np.abs(derivatives[order - b])
                for b in range(order + 1)
            )
            raw_lower = sum(
                comb(order - 1, b)
                * polynomials[b]
                * np.abs(derivatives[order - 1 - b])
                for b in range(order)
            )
            value_integrand = heat * raw_value
            first_integrand = heat * (u * raw_value + order * raw_lower)
            d0 = float(np.dot(weights, value_integrand))
            d1 = float(np.dot(weights, first_integrand))
            tail_mask = panel_index == len(DERIVATIVE_PANELS) - 1
            last_fraction = max(
                float(np.dot(weights[tail_mask], value_integrand[tail_mask]))
                / max(d0, 1e-300),
                float(np.dot(weights[tail_mask], first_integrand[tail_mask]))
                / max(d1, 1e-300),
            )
            max_last_panel_fraction = max(
                max_last_panel_fraction, last_fraction
            )
            rows.append(
                {
                    "N": retained,
                    "m": order,
                    "mu0": format(moment0, ".17e"),
                    "mu1": format(moment1, ".17e"),
                    "d0_envelope": format(d0, ".17e"),
                    "d1_envelope": format(d1, ".17e"),
                    "last_panel_fraction": format(last_fraction, ".17e"),
                }
            )
    return rows, {
        "nodes_per_panel": node_count,
        "n_cap": n_cap,
        "panels": [list(panel) for panel in DERIVATIVE_PANELS],
        "max_last_panel_fraction": format(max_last_panel_fraction, ".17e"),
    }


def compare_derivative_tables(
    coarse: list[dict],
    fine: list[dict],
    cap_check: list[dict],
) -> dict:
    max_coarse_fine = 0.0
    max_cap_delta = 0.0
    for left, right, capped in zip(coarse, fine, cap_check, strict=True):
        if (left["N"], left["m"]) != (right["N"], right["m"]):
            raise RuntimeError("derivative table key mismatch")
        if (capped["N"], capped["m"]) != (right["N"], right["m"]):
            raise RuntimeError("arithmetic-cap table key mismatch")
        for field in ("mu0", "mu1", "d0_envelope", "d1_envelope"):
            left_value = float(left[field])
            right_value = float(right[field])
            cap_value = float(capped[field])
            max_coarse_fine = max(
                max_coarse_fine,
                abs(left_value - right_value) / max(abs(right_value), 1e-300),
            )
            max_cap_delta = max(
                max_cap_delta,
                abs(cap_value - right_value) / max(abs(right_value), 1e-300),
            )
    return {
        "max_relative_coarse_fine_delta": format(
            max_coarse_fine, ".17e"
        ),
        "max_relative_ncap24_ncap28_delta": format(
            max_cap_delta, ".17e"
        ),
    }


def derivative_lookup(rows: list[dict]) -> dict[tuple[int, int], dict]:
    return {(row["N"], row["m"]): row for row in rows}


def jet_rows(node_count: int, dps: int, theta_terms: int) -> list[dict]:
    mp.mp.dps = dps
    quadrature_rows = saddle.quadrature_rows(node_count, theta_terms)
    full_column = [row[3] for row in quadrature_rows]
    output: list[dict] = []
    for t_text in TIMES:
        t = mp.mpf(t_text)
        for x_int in X_VALUES:
            x = mp.mpf(x_int)
            _, n5 = saddle.transition_index(5, t, x)
            _, n9 = saddle.transition_index(9, t, x)
            base = int(mp.ceil(max(n5, n9)))
            full_jet = saddle.jet_from_values(
                quadrature_rows, full_column, t, x
            )
            for kappa in KAPPAS:
                retained = base + kappa
                partial_column = [
                    mp.fsum(row[2][:retained]) for row in quadrature_rows
                ]
                partial_jet = saddle.jet_from_values(
                    quadrature_rows, partial_column, t, x
                )
                output.append(
                    {
                        "t": t_text,
                        "x": x_int,
                        "n_star": mp.nstr(max(n5, n9), 30),
                        "base_count": base,
                        "kappa": kappa,
                        "N": retained,
                        "full_H": mp.nstr(full_jet[0], dps - 10),
                        "full_H1": mp.nstr(full_jet[1], dps - 10),
                        "partial_H": mp.nstr(partial_jet[0], dps - 10),
                        "partial_H1": mp.nstr(partial_jet[1], dps - 10),
                    }
                )
    return output


def compare_jet_rows(coarse: list[dict], fine: list[dict]) -> dict:
    max_full_delta = mp.mpf(0)
    max_partial_delta = mp.mpf(0)
    for left, right in zip(coarse, fine, strict=True):
        key = ("t", "x", "kappa", "N")
        if tuple(left[item] for item in key) != tuple(
            right[item] for item in key
        ):
            raise RuntimeError("jet table key mismatch")
        for field in ("full_H", "full_H1"):
            lv = mp.mpf(left[field])
            rv = mp.mpf(right[field])
            max_full_delta = max(
                max_full_delta,
                abs(lv - rv) / max(abs(rv), mp.mpf("1e-100")),
            )
        for field in ("partial_H", "partial_H1"):
            lv = mp.mpf(left[field])
            rv = mp.mpf(right[field])
            max_partial_delta = max(
                max_partial_delta,
                abs(lv - rv) / max(abs(rv), mp.mpf("1e-100")),
            )
    return {
        "max_relative_full_jet_delta": mp.nstr(max_full_delta, 25),
        "max_relative_partial_jet_delta": mp.nstr(max_partial_delta, 25),
    }


def compose_c1_rows(
    jets: list[dict], derivatives: list[dict]
) -> list[dict]:
    lookup = derivative_lookup(derivatives)
    rows: list[dict] = []
    for jet in jets:
        retained = jet["N"]
        x = mp.mpf(jet["x"])
        partial = mp.mpf(jet["partial_H"])
        partial_prime = mp.mpf(jet["partial_H1"])
        order_rows: list[dict] = []
        first_passing_order: int | None = None
        for order in ORDERS:
            budget = lookup[(retained, order)]
            mu0 = mp.mpf(budget["mu0"])
            mu1 = mp.mpf(budget["mu1"])
            epsilon0 = min(
                mu0, mp.mpf(budget["d0_envelope"]) / x**order
            )
            epsilon1 = min(
                mu1, mp.mpf(budget["d1_envelope"]) / x**order
            )
            value_ratio = abs(partial) / epsilon0
            derivative_ratio = abs(4 * partial + x * partial_prime) / (
                4 * epsilon0 + x * epsilon1
            )
            passes = value_ratio > 1 or derivative_ratio > 1
            if passes and first_passing_order is None:
                first_passing_order = order
            order_rows.append(
                {
                    "m": order,
                    "epsilon0": mp.nstr(epsilon0, 25),
                    "epsilon1": mp.nstr(epsilon1, 25),
                    "value_ratio": mp.nstr(value_ratio, 25),
                    "derivative_ratio": mp.nstr(derivative_ratio, 25),
                    "diagnostic_disjunction_passes": bool(passes),
                }
            )
        rows.append(
            {
                **jet,
                "first_passing_m_5_to_9": first_passing_order,
                "orders": order_rows,
            }
        )
    return rows


def build_artifact(
    derivative_coarse_nodes: int,
    derivative_fine_nodes: int,
    jet_coarse_nodes: int,
    jet_fine_nodes: int,
    dps: int,
    theta_terms: int,
) -> dict:
    envelope = json.loads(ENVELOPE_RESULT.read_text(encoding="utf-8"))
    if envelope.get("kind") != (
        "jensen_window_pf_newman_theta_modular_tail_derivative_envelope_gate"
    ):
        raise RuntimeError("exact envelope source kind mismatch")
    functions = block_derivative_functions(max(ORDERS))
    coarse, coarse_meta = evaluate_derivative_table(
        functions, derivative_coarse_nodes, 28
    )
    fine, fine_meta = evaluate_derivative_table(
        functions, derivative_fine_nodes, 28
    )
    cap_check, cap_meta = evaluate_derivative_table(
        functions, derivative_fine_nodes, 24
    )
    derivative_convergence = compare_derivative_tables(
        coarse, fine, cap_check
    )
    coarse_jets = jet_rows(jet_coarse_nodes, dps, theta_terms)
    fine_jets = jet_rows(jet_fine_nodes, dps, theta_terms)
    jet_convergence = compare_jet_rows(coarse_jets, fine_jets)
    c1_rows = compose_c1_rows(fine_jets, fine)

    x200_k2 = [
        row
        for row in c1_rows
        if row["x"] == 200 and row["kappa"] == 2
    ]
    x300 = [row for row in c1_rows if row["x"] == 300]
    if not x200_k2 or not x300:
        raise RuntimeError("required stress rows missing")
    x300_by_kappa = {
        kappa: [
            row["first_passing_m_5_to_9"]
            for row in x300
            if row["kappa"] == kappa
        ]
        for kappa in KAPPAS
    }
    summary = {
        "x200_kappa2_first_passing_orders": [
            row["first_passing_m_5_to_9"] for row in x200_k2
        ],
        "x300_first_passing_orders_by_kappa": x300_by_kappa,
        "x300_kappa2_to_4_all_fail_m5_to_m9": all(
            value is None
            for kappa in (2, 3, 4)
            for value in x300_by_kappa[kappa]
        ),
        "x300_kappa5_all_pass": all(
            value is not None for value in x300_by_kappa[5]
        ),
        "scope": (
            "Quadrature and finite arithmetic truncation diagnostics only. "
            "The derivative envelopes are not interval enclosures, and these "
            "rows do not prove or disprove any all-frequency collar theorem."
        ),
    }
    rows = [
        ScoutRow(
            "ntmtdbs_01_m5_to_m9_envelopes",
            "finite_diagnostic",
            "diagnostic_only",
            "Positive-integrand modular-tail envelopes were evaluated for m=5 through 9 and N=4 through 10.",
            "d0_hat=integral E_(N,m); d1_hat=integral[uE_(N,m)+mE_(N,m-1)]",
            "Floating-point quadrature with an arithmetic cap; not a rigorous upper enclosure.",
            fine,
        ),
        ScoutRow(
            "ntmtdbs_02_envelope_convergence",
            "numerical_guard",
            "diagnostic_only",
            "Independent node ladders and arithmetic caps stress the derivative-envelope table.",
            "coarse/fine and n_cap=24/28 comparisons",
            "Convergence evidence only; no directed-rounding claim.",
            derivative_convergence,
        ),
        ScoutRow(
            "ntmtdbs_03_transition_jet_ladder",
            "finite_diagnostic",
            "diagnostic_only",
            "High-precision retained first jets were sampled across two Newman times and frequencies through x=312.",
            "N=ceil(max(n_*5,n_*9))+kappa, kappa=2,3,4,5",
            "Point rows only; no transition cell is intervalized.",
            fine_jets,
        ),
        ScoutRow(
            "ntmtdbs_04_direct_c1_composition",
            "finite_diagnostic",
            "diagnostic_only",
            "The approximate derivative envelopes were composed with the exact direct-C1 disjunction.",
            "max(|S|/epsilon0,|4S+xS'|/(4epsilon0+xepsilon1))>1",
            "Approximate sufficient-ratio scout, not a certificate.",
            c1_rows,
        ),
        ScoutRow(
            "ntmtdbs_05_fixed_collar_stress",
            "nonpromotion_gate",
            "guard_validated",
            "The tested fixed collars deteriorate immediately below an arithmetic saddle-count jump.",
            "At x=300, kappa=2,3,4 fail every sampled m=5..9; kappa=5 passes at both sampled times.",
            "Finite diagnostic guard; it rejects promotion of the sampled collars only.",
            summary,
        ),
        ScoutRow(
            "ntmtdbs_06_open_interval_handoff",
            "open_handoff",
            "not_ready_to_apply",
            "Rigorous work should intervalize the x^(3/4) count and its retained first jet, not a fixed additive collar.",
            "N=max(N_sad,ceil(K*(1+x)^(3/4)))",
            "No strict Laguerre theorem, Lambda<=0, RH, or Clay-prize result is proved.",
        ),
    ]
    return {
        "kind": STEM,
        "date": "2026-07-24",
        "status": (
            "finite modular-tail derivative-budget and direct-C1 "
            "transition stress scout"
        ),
        "proof_boundary": (
            "This artifact evaluates the exact positive derivative envelopes "
            "with convergent but non-interval quadrature and composes them with "
            "high-precision retained jets. It is diagnostic only. It does not "
            "supply directed-rounding tail bounds, intervalize a transition "
            "cell, prove a fixed-collar obstruction, prove retained first-jet "
            "separation, Lambda<=0, RH, or a Clay-prize result."
        ),
        "parameters": {
            "T_cap": T_CAP,
            "orders": list(ORDERS),
            "N_values": list(N_VALUES),
            "x_values": list(X_VALUES),
            "times": list(TIMES),
            "kappas": list(KAPPAS),
            "derivative_coarse": coarse_meta,
            "derivative_fine": fine_meta,
            "derivative_cap_check": cap_meta,
            "jet_coarse_nodes_per_panel": jet_coarse_nodes,
            "jet_fine_nodes_per_panel": jet_fine_nodes,
            "dps": dps,
            "theta_terms": theta_terms,
        },
        "provenance": {
            "generator": str(Path(__file__).relative_to(REPO_ROOT)).replace(
                "\\", "/"
            ),
            "generator_sha256": file_hash(Path(__file__)),
            "exact_envelope_result": str(
                ENVELOPE_RESULT.relative_to(REPO_ROOT)
            ).replace("\\", "/"),
            "exact_envelope_result_sha256": file_hash(ENVELOPE_RESULT),
        },
        "derivative_convergence": derivative_convergence,
        "jet_convergence": jet_convergence,
        "derivative_rows": fine,
        "jet_rows": fine_jets,
        "c1_rows": c1_rows,
        "summary": summary,
        "rows": [asdict(row) for row in rows],
        "sources": [
            "outputs/jensen_window_pf_newman_theta_modular_tail_derivative_envelope_gate.md",
            "outputs/jensen_window_pf_newman_theta_adaptive_modular_c1_remainder_contract.md",
            "outputs/formal_core.md",
        ],
    }


def render_note(artifact: dict) -> str:
    summary = artifact["summary"]
    derivative_rows = artifact["derivative_rows"]
    selected = [
        row
        for row in derivative_rows
        if row["m"] in (5, 9)
    ]
    table_lines = [
        (
            f"| {row['N']} | {row['m']} | {row['d0_envelope']} | "
            f"{row['d1_envelope']} |"
        )
        for row in selected
    ]
    return "\n".join(
        [
            "# Newman Theta Modular-Tail Derivative-Budget Scout",
            "",
            "Date: 2026-07-24",
            "",
            "Status: finite floating-point stress scout. This is not an",
            "interval certificate and not a proof of `Lambda<=0` or RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Derivative Envelopes",
            "",
            "At `T=1/5`, the exact positive-integrand formulas were evaluated",
            "for `m=5,...,9`. Selected rows are:",
            "",
            "| N | m | d0 envelope | d1 envelope |",
            "|---:|---:|---:|---:|",
            *table_lines,
            "",
            "The node ladder, arithmetic-cap comparison, and final-panel mass",
            "are recorded in the JSON artifact. They are convergence checks,",
            "not directed-rounding enclosures.",
            "",
            "## Direct-C1 Stress",
            "",
            "For each sampled row the scout composes",
            "",
            "```text",
            "max(|S|/epsilon0,|4S+xS'|/(4epsilon0+xepsilon1))>1",
            "N=ceil(max(n_*5,n_*9))+kappa.",
            "```",
            "",
            f"At `x=200`, `kappa=2` first passes at orders "
            f"`{summary['x200_kappa2_first_passing_orders']}` for the two",
            "sampled times. Immediately below the next arithmetic-count jump,",
            "`x=300`, every `kappa=2,3,4` row fails all sampled orders",
            "`m=5,...,9`, whereas `kappa=5` passes at both sampled times.",
            "",
            "This is a finite non-promotion guard. It does not prove that every",
            "fixed collar fails, and it does not certify the passing rows.",
            "",
            "## Revised Handoff",
            "",
            "The interval programme should test",
            "",
            "```text",
            "N=max(N_sad,ceil(K*(1+x)^(3/4)))",
            "```",
            "",
            "with rigorous derivative integrals and retained first-jet boxes.",
            "No strict Laguerre theorem, `Lambda<=0`, RH, or Clay-prize result",
            "is claimed.",
            "",
        ]
    )


def write_artifact(artifact: dict, out: Path, note: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    note.write_text(render_note(artifact), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--derivative-coarse-nodes", type=int, default=180)
    parser.add_argument("--derivative-fine-nodes", type=int, default=260)
    parser.add_argument("--jet-coarse-nodes", type=int, default=180)
    parser.add_argument("--jet-fine-nodes", type=int, default=220)
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--theta-terms", type=int, default=15)
    args = parser.parse_args()
    artifact = build_artifact(
        args.derivative_coarse_nodes,
        args.derivative_fine_nodes,
        args.jet_coarse_nodes,
        args.jet_fine_nodes,
        args.dps,
        args.theta_terms,
    )
    write_artifact(artifact, args.out, args.note)
    print(
        "wrote Newman theta modular-tail derivative-budget scout: "
        f"{len(artifact['derivative_rows'])} derivative rows, "
        f"{len(artifact['c1_rows'])} C1 rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
