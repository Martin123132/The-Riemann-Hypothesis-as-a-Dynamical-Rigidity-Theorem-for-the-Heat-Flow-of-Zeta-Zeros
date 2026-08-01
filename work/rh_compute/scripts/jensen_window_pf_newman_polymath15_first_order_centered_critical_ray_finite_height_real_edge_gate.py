#!/usr/bin/env python3
"""Build the full critical-ray finite-height real-edge sign gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = REPO_ROOT / "work/rh_compute/scripts"
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_polymath15_first_order_centered_q1_finite_height_real_edge_remainder_gate as q1  # noqa: E402


STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "critical_ray_finite_height_real_edge_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "q1_finite_height_gate": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "q1_finite_height_real_edge_remainder_gate.json"
    ),
    "cofinal_real_edge_gate": q1.SOURCE_PATHS["cofinal_real_edge"],
    "complex_endpoint_source": q1.SOURCE_PATHS["complex_endpoint_source"],
    "critical_c1_source": q1.SOURCE_PATHS["critical_c1_source"],
    "first_dirichlet_correction": q1.SOURCE_PATHS[
        "first_dirichlet_correction"
    ],
}

H_DENOMINATOR = q1.H_DENOMINATOR
T_MAX = Fraction(1, 2)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def source_audit() -> dict:
    missing = [str(path) for path in SOURCE_PATHS.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCE_PATHS.items()
        },
        "published_source": {
            "citation": (
                "D.H.J. Polymath, Effective approximation of heat flow "
                "evolution of the Riemann xi function, and a new upper "
                "bound for the de Bruijn-Newman constant, arXiv:1904.12438"
            ),
            "location": "equation (53), C_0, first correction, and Proposition 6.2",
            "url": "https://arxiv.org/abs/1904.12438",
        },
    }


def exact_payload() -> dict:
    inherited = q1.exact_payload()
    return {
        "domain": {
            "critical_ray": (
                "L>=50, 0<tL<=25, theta=(1-p)/2 in [0,1], "
                "a=N+theta, h=1/a, T_0=2pi/h^2"
            ),
            "effective_box": (
                "0<t<=1/2, 0<h<exp(-25)<1/72000000000, "
                "0<=theta<=1"
            ),
            "q_consequence": (
                "q=2tL^2; therefore the theorem contains every q>=1 "
                "point in the critical range 0<tL<=25"
            ),
        },
        "endpoint": inherited["endpoint"],
        "terminal": inherited["terminal"],
        "stable_log": inherited["stable_log"],
        "endpoint_rate": inherited["endpoint_rate"],
        "normalized_edge": inherited["normalized_edge"],
        "extension": {
            "unchanged_identities": (
                "The q=1 proof used no q identity after replacing t by its "
                "upper bound. Equations for H, d, W, W_x, mu, mu_x, and "
                "the four normalized edge jets are exact for every "
                "0<tL<=25 point on a prescribed-N critical chart."
            ),
            "only_changed_caps": (
                "rho/h^2: 1/80000 -> 1/32; "
                "t/4+t^2|w|^2/8: <1/10000 -> <1/2; "
                "the stable heat terms are then re-audited before exponentiation"
            ),
        },
    }


def majorant_certificate() -> dict:
    h0 = Fraction(1, H_DENOMINATOR)
    t = T_MAX
    rho_coefficient = t / 16
    rho_max = rho_coefficient * h0**2
    u_coefficient = 1 / (1 - h0)
    u_remainder_coefficient = 1 / (2 * (1 - h0))
    log_r_coefficient = rho_coefficient / (1 - rho_max)
    y_coefficient = Fraction(1, 12) / (1 - rho_max)
    eta_coefficient = Fraction(7, 2) * y_coefficient
    chi_r_coefficient = (
        log_r_coefficient / 2
        + Fraction(5, 4) * y_coefficient**2 * h0**2
    )
    u_r_coefficient = u_coefficient + chi_r_coefficient * h0
    chi_absolute = (
        chi_r_coefficient * h0**2 + 1 + eta_coefficient * h0**2
    )
    w_absolute = u_r_coefficient * h0 + 1 + eta_coefficient * h0**2
    alpha_prime_coefficient = (
        6 * y_coefficient + 7 * y_coefficient**2 * h0**2
    )
    alpha_second_coefficient = (
        2 * y_coefficient**2
        + 24 * y_coefficient**3 * h0**2
    )
    correction_block = t / 4 + t**2 * 9 / 8
    d_coefficient = (
        y_coefficient / 3
        + alpha_prime_coefficient * correction_block
    )
    d_x_coefficient = Fraction(1, 2) * (
        Fraction(2, 3) * y_coefficient**2
        + alpha_second_coefficient * correction_block
        + (t**2 / 4) * 3 * alpha_prime_coefficient**2
    )
    log_one_plus_d_coefficient = d_coefficient / (
        1 - d_coefficient * h0**2
    )

    real_w_coefficient = (
        u_coefficient / 2
        + Fraction(7, 4) * log_r_coefficient * h0
        + Fraction(7, 8) * y_coefficient**2 * h0**3
        + y_coefficient**2 * h0**3 / 12
        + (t / 4)
        * (
            u_r_coefficient**2 * h0
            + 2 * eta_coefficient * h0
            + eta_coefficient**2 * h0**3
        )
        + log_one_plus_d_coefficient * h0
    )
    imaginary_w_coefficient = (
        Fraction(8, 3) / (1 - h0)
        + (t / 2) * u_coefficient
        + t**2 * h0 / (128 * (1 - rho_max))
        + Fraction(3, 8) * y_coefficient * h0
        + (t / 4)
        * 2
        * u_r_coefficient
        * (1 + eta_coefficient * h0**2)
        + log_one_plus_d_coefficient * h0
    )
    w_coefficient = real_w_coefficient + imaginary_w_coefficient
    w_x_coefficient = (
        eta_coefficient / 2
        + Fraction(1, 8)
        + (u_remainder_coefficient + chi_r_coefficient) / 2
        + (t / 4) * alpha_prime_coefficient * 3
        + d_x_coefficient
        * h0**2
        / (1 - d_coefficient * h0**2)
    )

    mu_coefficient = (
        Fraction(1, 24)
        + Fraction(1, 12)
        + eta_coefficient / 2
        + chi_r_coefficient / 2
        + (t / 4) * alpha_prime_coefficient * 3
    )
    chi_x_real_coefficient = y_coefficient / 2 * (
        3 * y_coefficient**2 * h0**2 + rho_coefficient
    )
    chi_x_imag_coefficient = Fraction(7, 2) * y_coefficient**2
    chi_x_coefficient = chi_x_real_coefficient + chi_x_imag_coefficient
    mu_x_coefficient = (
        Fraction(1, 96)
        + chi_x_coefficient / 2
        + (t / 4)
        * (
            alpha_second_coefficient * 3 / 2
            + alpha_prime_coefficient * chi_x_coefficient * h0**2
        )
    )
    s1_defect_coefficient = t * alpha_prime_coefficient / 4
    real_s1_coefficient = t * 6 * y_coefficient / 4
    s2_coefficient = t * alpha_second_coefficient / 8

    checks = {
        "u_over_h": (u_coefficient, Fraction(2)),
        "u_minus_h_theta_over_h2": (u_remainder_coefficient, Fraction(1)),
        "rho_over_h2": (rho_coefficient, Fraction(1, 16)),
        "log_r_over_h2": (log_r_coefficient, Fraction(1, 16)),
        "y_over_h2": (y_coefficient, Fraction(1, 6)),
        "eta_over_h2": (eta_coefficient, Fraction(1)),
        "chi_R_over_h2": (chi_r_coefficient, Fraction(1, 16)),
        "chi_absolute": (chi_absolute, Fraction(2)),
        "w_absolute": (w_absolute, Fraction(2)),
        "alpha_prime_over_h2": (alpha_prime_coefficient, Fraction(1)),
        "alpha_second_over_h4": (alpha_second_coefficient, Fraction(1)),
        "correction_block": (correction_block, Fraction(1, 2)),
        "d_over_h2": (d_coefficient, Fraction(1, 4)),
        "d_x_over_h4": (d_x_coefficient, Fraction(1, 32)),
        "Re_W_over_h": (real_w_coefficient, Fraction(1)),
        "Im_W_over_h": (imaginary_w_coefficient, Fraction(13, 4)),
        "W_over_h": (w_coefficient, Fraction(4)),
        "W_x_over_h2": (w_x_coefficient, Fraction(1)),
        "mu_over_h2": (mu_coefficient, Fraction(1)),
        "mu_x_over_h4": (mu_x_coefficient, Fraction(1)),
        "s1_defect_over_h2": (s1_defect_coefficient, Fraction(1, 15)),
        "real_s1_over_h2": (real_s1_coefficient, Fraction(1, 15)),
        "s2_over_h4": (s2_coefficient, Fraction(1, 1000)),
    }
    failed = [name for name, (value, cap) in checks.items() if not value < cap]
    if failed:
        raise RuntimeError("failed critical-ray majorants: " + ", ".join(failed))
    if not (4 * q1.aq(h0)).exp() < 2:
        raise RuntimeError("exp(4h)<2 audit failed")

    # The endpoint estimates are unchanged.  The enlarged terminal caps are
    # |exp(W)-1|<8h and |Z_x-i*h*theta*Q/2|<6h^2.
    c_defect = Fraction(86 + 8) + Fraction(2, 3) * h0
    c_x_defect = Fraction(232 + 6) + h0 / 2 + h0**2 / 18
    d_edge_defect = (
        Fraction(236 + 1 + 4)
        + h0 / 2
        + 4 * h0 / 15
        + Fraction(2, 3) * h0
    )
    d_edge_x_defect = (
        Fraction(810)
        + Fraction(2, 3) * h0
        + h0**2 / 24
        + 4 * h0**2 / 1000
        + Fraction(1, 24) * (4 + 2 * h0 / 15)
        + 2 * h0 / 15
        + Fraction(1, 2)
        + 3
        + 5 * (h0**2 / 500 + h0 / 360)
        + 8 * h0 / 15
    )
    jet_checks = {
        "c_defect_over_h": (c_defect, Fraction(100)),
        "c_x_defect_over_h2": (c_x_defect, Fraction(250)),
        "d_edge_defect_over_h2": (d_edge_defect, Fraction(250)),
        "d_edge_x_defect_over_h3": (d_edge_x_defect, Fraction(900)),
    }
    failed_jets = [
        name for name, (value, cap) in jet_checks.items() if not value < cap
    ]
    if failed_jets:
        raise RuntimeError("failed critical-ray edge jets: " + ", ".join(failed_jets))

    all_checks = {**checks, **jet_checks}
    chain = {
        name: {
            "upper": fraction_text(value),
            "cap": fraction_text(cap),
        }
        for name, (value, cap) in all_checks.items()
    }

    e0, e1, e2, e3 = 100, 250, 250, 900
    linear = e0 * 3 + 4 * e3 + 3 * (e1 + e2)
    quadratic = e0 * e3 + e1 * e2
    if linear != 5400 or quadratic != 152500:
        raise RuntimeError("critical-ray current arithmetic drifted")
    combined = Fraction(linear) + h0 * quadratic
    if not combined < 5500:
        raise RuntimeError("critical-ray 5500h budget failed")
    if not Fraction(5500, H_DENOMINATOR) < Fraction(1, 10_000_000):
        raise RuntimeError("critical-ray finite-height target failed")

    return {
        "L_min": 50,
        "critical_range": "0<tL<=25",
        "t_max": "1/2",
        "h_upper": f"1/{H_DENOMINATOR}",
        "normalized_majorant_chain": chain,
        "source_bounds": {
            "rho": "rho<h^2/16",
            "log_r": "|log r|<h^2/16",
            "alpha_prime": "|alpha'|<h^2",
            "alpha_second": "|alpha''|<h^4",
            "correction": "|d|<h^2/4, |d_x|<h^4/32",
            "stable_log": "|W|<4h, |W_x|<h^2",
            "terminal": "|Z+Q|<8h, |Z_x-i*h*theta*Q/2|<6h^2",
            "endpoint_rate": "|mu|<h^2, |mu_x|<h^4",
        },
        "edge_jet_defects": {
            "c": "c=A+h*e_0, |e_0|<100",
            "d": "d_edge=h*B+h^2*e_1, |e_1|<250",
            "c_x": "c_x=h*B+h^2*e_2, |e_2|<250",
            "d_x": "d_(edge,x)=h^2*M+h^3*e_3, |e_3|<900",
        },
        "linear_error_constant": linear,
        "quadratic_error_constant": quadratic,
        "uniform_remainder": (
            "|a^2 J_edge/S_a^2-K_edge(p)|<5500h<1/10000000"
        ),
        "finite_height_margin": (
            "a^2 J_edge/S_a^2<-3749/10000000<0"
        ),
    }


def build_rows(exact: dict, derivatives: dict, budget: dict) -> list[GateRow]:
    return [
        GateRow(
            "crfher_01_domain",
            "effective_domain",
            "ready_to_apply",
            "The complete critical range lies in one explicit parameter box.",
            exact["domain"]["critical_ray"] + "; " + exact["domain"]["effective_box"],
            "A prescribed-N chart is used before adjacent-cutoff splicing.",
        ),
        GateRow(
            "crfher_02_exact_source",
            "source_inheritance",
            "ready_to_apply",
            "The stable endpoint-terminal identities do not require q=1.",
            exact["extension"]["unchanged_identities"],
            "Only the parameter majorants are enlarged.",
        ),
        GateRow(
            "crfher_03_changed_caps",
            "analytic_majorant",
            "certified",
            "The two q=1-specific numerical caps are replaced before bounding.",
            exact["extension"]["only_changed_caps"],
            "No q=1 constant is reused without recomputation.",
        ),
        GateRow(
            "crfher_04_stable_heat",
            "exact_desingularization",
            "ready_to_apply",
            "The terminal heat constants still cancel before absolute values.",
            exact["stable_log"]["removable_rewrites"],
            "The continuous W->0 branch is retained on the enlarged box.",
        ),
        GateRow(
            "crfher_05_source_bounds",
            "analytic_majorant",
            "certified",
            "The full critical range preserves the terminal and endpoint rates.",
            "|W|<4h, |W_x|<h^2, |mu|<h^2, |mu_x|<h^4",
            "All normalized rational caps are independently stored.",
            diagnostics=budget["source_bounds"],
        ),
        GateRow(
            "crfher_06_source_derivatives",
            "interval_inheritance",
            "certified",
            "The C_0 derivative theorem is independent of t and transfers unchanged.",
            "4096 direct Arb boxes plus two Cauchy charts bound C_0 through order five",
            "The interval certificate is recomputed by this gate.",
            diagnostics={
                "boxes": derivatives["direct_boxes"],
                "bounds": derivatives["global_C0_derivative_bounds"],
            },
        ),
        GateRow(
            "crfher_07_terminal_jets",
            "finite_height_jet",
            "certified",
            "The enlarged heat range retains one-order terminal defects.",
            "|Z+Q|<8h and |Z_x-i*h*theta*Q/2|<6h^2",
            "No phase interpolation is used.",
        ),
        GateRow(
            "crfher_08_four_jets",
            "finite_height_jet",
            "certified",
            "The same four edge-jet caps hold on the complete critical range.",
            "; ".join(budget["edge_jet_defects"].values()),
            "All H_xx, mu_x, s_*'', u_x, and Z_x terms are retained.",
        ),
        GateRow(
            "crfher_09_current_budget",
            "division_free_perturbation",
            "certified",
            "The critical-ray finite-height remainder stays below the q=1 target.",
            budget["uniform_remainder"],
            "Polynomial at A=0 and at every zero real projection.",
            diagnostics={
                "linear": budget["linear_error_constant"],
                "quadratic": budget["quadratic_error_constant"],
            },
        ),
        GateRow(
            "crfher_10_uniform_sign",
            "effective_interval_theorem",
            "certified",
            "The retained first-order real-edge current is clockwise throughout the critical range.",
            budget["finite_height_margin"],
            "This is a model-edge theorem on one prescribed-N chart.",
        ),
        GateRow(
            "crfher_11_q_ge_1",
            "domain_corollary",
            "ready_to_apply",
            "Every q>=1 critical point is included.",
            exact["domain"]["q_consequence"],
            "No adjacent-cutoff signed splice follows from pointwise chart signs alone.",
        ),
        GateRow(
            "crfher_12_nonpromotion",
            "nonpromotion_guard",
            "guard_validated",
            "The model sign is not silently promoted to the complete Xi edge.",
            "The omitted higher-order Xi/source remainder and global contact insertion remain separate.",
            "No Xi-level edge sign, Abel gap, winding cap, Lambda<=0, or RH theorem is claimed.",
        ),
        GateRow(
            "crfher_13_handoff",
            "route_decision",
            "open_quantitative_handoff",
            "The pointwise q>=1 retained-model edge wall is closed.",
            "Prove the signed adjacent-cutoff splice, then insert the edge block into the cumulative contact scalar with the existing global remainder budgets.",
            "The splice and cumulative Abel-scalar gap remain open.",
        ),
    ]


def render_note(artifact: dict) -> str:
    budget = artifact["majorant_certificate"]
    return "\n".join(
        [
            "# Newman Polymath-15 Critical-Ray Finite-Height Real-Edge Gate",
            "",
            "Date: 2026-07-31",
            "",
            "Status: effective sign for the retained first-order real-edge model on the full critical range. This is not a proof of an Xi-level edge sign, Lambda <= 0, or RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Domain Extension",
            "",
            "For L>=50 and 0<tL<=25,",
            "",
            "```text",
            "0<t<=1/2,          h=1/a<exp(-25)<1/72000000000.",
            "```",
            "",
            "The exact stable formulas from the q=1 gate depend on t but do not use q=1. This gate redoes every t-dependent majorant on the enlarged box rather than extrapolating the q=1 constants.",
            "",
            "## Enlarged Source Bounds",
            "",
            "The two visibly changed inputs are",
            "",
            "```text",
            "rho/h^2 <= 1/32,",
            "t/4+t^2|w|^2/8 < 1/2.",
            "```",
            "",
            "After the same constant and removable-quotient cancellations, the full normalized chain gives",
            "",
            "```text",
            "|alpha'|<h^2,            |alpha''|<h^4,",
            "|d|<h^2/4,               |d_x|<h^4/32,",
            "|W|<4h,                  |W_x|<h^2,",
            "|Z+Q|<8h,                |Z_x-i h theta Q/2|<6h^2,",
            "|mu|<h^2,                |mu_x|<h^4.",
            "```",
            "",
            "The checker stores and re-evaluates all 27 normalized rational inequalities. The 4096-box C_0 derivative cover and both factored Cauchy charts are recomputed as well.",
            "",
            "## Four-Jet Transfer",
            "",
            "The complete normalized edge still satisfies",
            "",
            "```text",
            "c=A+h e_0,                       |e_0|<100,",
            "d_edge=hB+h^2 e_1,               |e_1|<250,",
            "c_x=hB+h^2 e_2,                  |e_2|<250,",
            "d_(edge,x)=h^2M+h^3 e_3,         |e_3|<900.",
            "```",
            "",
            "Therefore the same division-free perturbation budget holds:",
            "",
            "```text",
            budget["uniform_remainder"],
            budget["finite_height_margin"],
            "```",
            "",
            "uniformly for L>=50, 0<tL<=25, and -1<=p<=1.",
            "",
            "## q>=1 Consequence",
            "",
            "Since q=2tL^2, every q>=1 point in the critical range is included. Thus the pointwise prescribed-chart q>=1 extension requested after the q=1 gate is closed for the retained first-order model.",
            "",
            "## Boundary and Handoff",
            "",
            "A pointwise sign on each prescribed-N chart does not itself sign the transition when the cutoff changes. The next gate is the adjacent-cutoff signed splice. The omitted higher-order Xi/source remainder and its insertion into the cumulative contact scalar remain separate. No Xi-level edge sign, Abel-scalar gap, winding cap, contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
            "",
        ]
    )


def build_artifact() -> dict:
    q1.verify_symbolics()
    derivatives = q1.c0_derivative_certificate()
    exact = exact_payload()
    budget = majorant_certificate()
    rows = build_rows(exact, derivatives, budget)
    return {
        "kind": STEM,
        "date": "2026-07-31",
        "status": (
            "effective finite-height negative critical-ray real-edge "
            "current for the retained first-order model"
        ),
        "source_audit": source_audit(),
        "exact": exact,
        "derivative_certificate": derivatives,
        "majorant_certificate": budget,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "direct_arb_boxes": derivatives["direct_boxes"],
            "removable_cauchy_charts": 2,
            "normalized_majorant_inequalities": len(
                budget["normalized_majorant_chain"]
            ),
            "finite_edge_jet_defects": 4,
            "uniform_critical_ray_model_signs": 1,
            "uniform_q_ge_1_model_signs": 1,
            "adjacent_cutoff_signed_splices": 0,
            "xi_level_edge_signs": 0,
            "abel_gaps": 0,
            "winding_bounds": 0,
        },
        "proof_boundary": (
            "This artifact proves an effective negative real-edge current "
            "for the retained first-order model on L>=50, 0<tL<=25, and "
            "-1<=p<=1, including every q>=1 point in that critical range. "
            "It does not prove an adjacent-cutoff signed splice, bound the "
            "omitted higher-order Xi/source contribution to this current, "
            "prove an Xi-level edge sign, cumulative-minor estimate, "
            "Abel-scalar gap, successor winding cap, contact exclusion, "
            "Q209, the cofinal descendant theorem, Lambda<=0, PF-infinity, "
            "RH, or a prize-level conclusion."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built critical-ray finite-height real-edge gate: 13 rows, "
        "4096 Arb boxes, 27 normalized majorants, 4 finite edge jets, "
        "1 full critical-ray model sign, 1 q>=1 model sign, "
        "0 cutoff splices, 0 Xi-level signs"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
