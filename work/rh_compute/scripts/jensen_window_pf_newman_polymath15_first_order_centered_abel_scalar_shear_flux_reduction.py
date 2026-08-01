#!/usr/bin/env python3
"""Build the endpoint-complete Abel-scalar shear-flux reduction."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
import math
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "abel_scalar_shear_flux_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "carrier_kernel": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "carrier_kernel_abel_prefix_reduction.json"
    ),
    "normalized_flux": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "normalized_prefix_phase_flux_reduction.json"
    ),
    "phase_cylinders": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "phase_cylinder_jacobian_transport_guard.json"
    ),
    "oriented_successor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_oriented_successor_winding_reduction.json"
    ),
    "real_residual": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "real_residual_reduction.json"
    ),
    "degree_transfer": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_cofinal_boundary_degree_transfer_target.json"
    ),
    "odd_fibre_guard": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "endpoint_odd_fibre_correlation_feasibility_guard.json"
    ),
}


@dataclass(frozen=True)
class FluxRow:
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
        "carrier_kernel": (
            "mathsf_A=Re(W_A)=mathcal_C_N+c*u_N*mathsf_X",
            "even when W_0=0",
        ),
        "normalized_flux": (
            "Gamma_ell=mathsf_X+i*mathsf_A/ell",
            "0<=kappa_j",
        ),
        "phase_cylinders": (
            "J_(theta,tau)",
            "H_2-log(a)H_1",
        ),
        "oriented_successor": (
            "kappa_j=sum_(p in D_j)",
            "one_sided_trap",
        ),
        "real_residual": (
            "|E_[1]|<101*exp(L/4)",
            "exact_crossing_sign",
        ),
        "degree_transfer": (
            "derivative shear",
            "determinant",
        ),
        "odd_fibre_guard": (
            "falsification coordinate",
            "|mathcal_C_N|>A_L+epsilon_term",
        ),
    }
    for key, required in markers.items():
        text = json.dumps(payloads[key], sort_keys=True)
        for marker in required:
            if marker not in text:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {
        key: str(payload.get("kind", ""))
        for key, payload in payloads.items()
    }


def shear_flux_audit() -> dict[str, str]:
    x_value, scalar, alpha, ell, homotopy = sp.symbols(
        "X C alpha ell s", real=True
    )
    x_rate, scalar_rate, ell_rate = sp.symbols(
        "X_s C_s ell_s", real=True
    )
    shear = sp.Matrix(
        [
            [1, 0],
            [homotopy * alpha / ell, 1],
        ]
    )
    if sp.simplify(shear.det() - 1) != 0:
        raise RuntimeError("Abel-scalar shear determinant failed")

    proxy = x_value + sp.I * scalar / ell
    proxy_rate = x_rate + sp.I * (
        scalar_rate / ell - scalar * ell_rate / ell**2
    )
    numerator = sp.simplify(
        sp.im(sp.conjugate(proxy) * proxy_rate)
    )
    denominator = sp.simplify(
        sp.re(sp.conjugate(proxy) * proxy)
    )
    expected_numerator = (
        ell * x_value * scalar_rate
        - x_value * scalar * ell_rate
        - ell * scalar * x_rate
    )
    expected_denominator = ell**2 * x_value**2 + scalar**2
    if sp.simplify(
        numerator / denominator
        - expected_numerator / expected_denominator
    ) != 0:
        raise RuntimeError("Abel-scalar argument flux failed")
    crossing_speed = sp.simplify(
        (numerator / denominator).subs(x_value, 0)
    )
    if sp.simplify(crossing_speed + ell * x_rate / scalar) != 0:
        raise RuntimeError("Abel-scalar crossing speed failed")

    return {
        "shear_matrix": (
            "[[1,0],[s*alpha/ell,1]] with determinant 1."
        ),
        "flux": (
            "partial_s arg(Psi_ell)="
            "{ell*mathsf_X*(mathcal_C_N)_s"
            "-mathsf_X*mathcal_C_N*ell_s"
            "-ell*mathcal_C_N*(mathsf_X)_s}/"
            "{ell^2*mathsf_X^2+mathcal_C_N^2}."
        ),
        "crossing_speed": (
            "At mathsf_X=0 and mathcal_C_N!=0, "
            "partial_s arg(Psi_ell)="
            "-ell*(mathsf_X)_s/mathcal_C_N."
        ),
    }


def orientation_audit() -> dict[str, str | float | list[float]]:
    level = 50.0
    ratios = [
        math.exp(-1.5 * level) / (4.0 * level),
        math.exp(-level) / level,
        101.0
        * math.exp(-0.5 * level)
        / (100000.0 * level + 1.0),
        1.0e-7 / (100000.0 * level + 1.0),
    ]
    total = sum(ratios)
    if not total < 2.03e-14:
        raise RuntimeError("crossing-orientation budget drifted")
    return {
        "error_bound": (
            "On |mathsf_X|<=delta_L, "
            "|partial_x mathsf_X-mathcal_C_N| is at most "
            "epsilon_term+2*exp(-L)*delta_L"
            "+[101*exp(-7L/4)+1e-7*exp(-5L/4)]/|f_1|."
        ),
        "ratio_bound": (
            "After division by A_L, the four terms are below "
            "exp(-3L/2)/(4L), exp(-L)/L, "
            "101*exp(-L/2)/(100000L+1), and "
            "1e-7/(100000L+1). They decrease for L>=50 and "
            "sum to less than 2.03e-14 at L=50."
        ),
        "at_L_50_terms": ratios,
        "at_L_50_total": total,
    }


def heat_recrossing_audit() -> dict[str, str | dict[str, int]]:
    x_value, time, frequency, phase = sp.symbols(
        "x t m phi", real=True, positive=True
    )
    field = sp.exp(frequency**2 * time) * sp.sin(
        frequency * x_value + phase
    )
    slope = sp.diff(field, x_value)
    if sp.simplify(sp.diff(field, time) + sp.diff(field, x_value, 2)) != 0:
        raise RuntimeError("backward-heat Fourier audit failed")
    flux_numerator = sp.trigsimp(
        field * sp.diff(slope, x_value)
        - slope * sp.diff(field, x_value)
    )
    expected = -frequency**2 * sp.exp(2 * frequency**2 * time)
    if sp.simplify(flux_numerator - expected) != 0:
        raise RuntimeError("backward-heat winding audit failed")

    counts: dict[str, int] = {}
    for integer_frequency in (1, 2, 5, 9):
        zero_count = 0
        upward_count = 0
        for index in range(-2, 2 * integer_frequency + 3):
            location = (
                index * math.pi - math.pi / 4
            ) / integer_frequency
            if 0.0 <= location < 2.0 * math.pi:
                zero_count += 1
                if index % 2 == 0:
                    upward_count += 1
        if zero_count != 2 * integer_frequency:
            raise RuntimeError("Fourier zero count failed")
        if upward_count != integer_frequency:
            raise RuntimeError("Fourier upward count failed")
        counts[str(integer_frequency)] = upward_count

    return {
        "family": (
            "For integer m>=1, "
            "X_m(x,t)=exp(m^2*t)*sin(m*x+pi/4) satisfies "
            "partial_t X_m=-partial_x^2 X_m. Put "
            "mathcal_C_m=partial_x X_m and alpha=0."
        ),
        "gap": (
            "At t=0, on |X_m|<=delta<1, "
            "|mathcal_C_m|>=m*sqrt(1-delta^2). "
            "The pointwise crossing gap can be made arbitrarily large."
        ),
        "winding": (
            "On the half-open interval 0<=x<2*pi there are exactly "
            "m upward crossings, and "
            "wind[X_m+i*mathcal_C_m]=-m. The argument-flux numerator "
            "is -m^2*exp(2m^2*t)<0."
        ),
        "checked_upward_counts": counts,
    }


def build_exact() -> dict:
    shear = shear_flux_audit()
    orientation = orientation_audit()
    heat = heat_recrossing_audit()
    return {
        "pi_provenance": (
            "The pi in a^2=x/(4pi)+t/16 and the Xi endpoint remains "
            "inherited from the completed-zeta normalization and the "
            "Riemann-Siegel saddle. The 2*pi in a winding normalization "
            "is the period of exp(i*theta). The pi used in the generic "
            "sine guard is only the standard sine period and introduces "
            "no new constant into the Xi formulas."
        ),
        "slope_shear_identity": (
            "On every fixed-N q>=1 chart, put alpha=c*u_N. The exact "
            "all-fiber identity is "
            "mathsf_A=mathcal_C_N+alpha*mathsf_X."
        ),
        "pointwise_gap_nonvanishing": (
            "If |mathsf_X|<=delta_L implies "
            "|mathcal_C_N|>A_L+epsilon_term, then "
            "(mathsf_X,mathcal_C_N) is nonzero everywhere on that "
            "q>=1 arc: inside the band the second coordinate is nonzero, "
            "and outside it the first coordinate is nonzero."
        ),
        "shear_homotopy": (
            "For 0<=s<=1 and ell>0 define "
            "Gamma_(s,ell)=mathsf_X+i*"
            "(mathcal_C_N+s*alpha*mathsf_X)/ell. "
            "A zero would force mathsf_X=mathcal_C_N=0. Thus the "
            "pointwise gap gives a nonvanishing homotopy from "
            "Psi_ell=mathsf_X+i*mathcal_C_N/ell to "
            "Gamma_ell=mathsf_X+i*mathsf_A/ell."
        ),
        "shear_matrix": shear["shear_matrix"],
        "positive_scale_homotopy": (
            "Changing ell through positive values is the diagonal "
            "orientation-preserving map diag(1,1/ell). Once "
            "(mathsf_X,mathcal_C_N) is nonzero, it preserves winding "
            "and positive-imaginary-ray intersection signs."
        ),
        "reduced_proxy": (
            "The q>=1 signed integer may therefore be computed from the "
            "division-free Abel proxy "
            "Psi_ell=mathsf_X+i*mathcal_C_N/ell. This reduction includes "
            "W_0=0 and removes the terminal shear c*u_N*mathsf_X from "
            "the topological coordinate."
        ),
        "reduced_derivative": (
            "Inside a fixed-N chart, alpha_x=Re(s_*'')*u_N"
            "+c/(4*T_0) and "
            "(mathcal_C_N)_x=(mathsf_A)_x-alpha_x*mathsf_X"
            "-alpha*(mathsf_X)_x. The existing endpoint-complete formulas "
            "for (mathsf_X)_x and (mathsf_A)_x retain "
            "H_0,H_1,H_2,D_0,D_1 and all endpoint derivatives, so the "
            "reduced flux remains an O(N) prefix computation."
        ),
        "general_flux": shear["flux"],
        "crossing_speed": shear["crossing_speed"],
        "orientation_error": orientation["error_bound"],
        "orientation_ratio": orientation["ratio_bound"],
        "orientation_consequence": (
            "Under the pointwise gap, "
            "|mathcal_C_N|>A_L+epsilon_term>A_L exceeds the certified "
            "orientation error. Hence sign(partial_x mathsf_X)="
            "sign(mathcal_C_N) throughout the contact band. At "
            "mathsf_X=0, mathcal_C_N>0 is exactly an upward crossing "
            "and contributes -1 on an increasing-x edge."
        ),
        "zero_fibre": (
            "At W_0=0 one has mathsf_X=mathsf_Y=0 and "
            "Psi_ell=i*mathcal_C_N/ell. The pointwise gap keeps this "
            "value nonzero; no ratio, Wronskian, or Schur coordinate "
            "is used."
        ),
        "chart_boundary": (
            "The shear is exact inside each canonical fixed-N chart. "
            "The already proved adjacent real-projection homotopies "
            "must still be used at cutoff changes; no complex "
            "chart-invariance of mathcal_C_N is inferred."
        ),
        "heat_family": heat["family"],
        "heat_gap": heat["gap"],
        "heat_winding": heat["winding"],
        "heat_nonpromotion": (
            "The backward-heat Fourier family proves that an arbitrarily "
            "strong pointwise crossing gap, exact identity "
            "mathcal_C=partial_x mathsf_X, and heat evolution do not "
            "bound the one-sided crossing count. It is a generic route "
            "guard, not an Xi counterexample."
        ),
        "q_ge_1_pointwise_target": (
            "Pointwise Xi target: prove "
            "|mathsf_X|<=delta_L implies "
            "|mathcal_C_N|>A_L+epsilon_term in each prescribed "
            "canonical q>=1 chart, retaining W_0=0 and the endpoint."
        ),
        "q_ge_1_flux_target": (
            "Independent signed Xi target: compose Psi_ell through the "
            "q>=1 horizontal arcs, adjacent cutoff homotopies, vertical "
            "connector, finite shoulders, and chart joins, and prove "
            "Delta_(partial D_j)arg(Psi_j)<2*pi. Since the transferred "
            "integer kappa_j is already nonnegative, this gives "
            "0<=kappa_j<1 and hence kappa_j=0."
        ),
        "minimal_route_verdict": (
            "The smallest surviving linked-point route consists of two "
            "independent Xi theorems: the Abel-scalar contact-band gap "
            "and the endpoint-complete Abel-proxy phase budget below "
            "one turn. Full Schur disk stability and full-circle "
            "three-cylinder zero-freeness are stronger than required."
        ),
        "q_lt_1_target": (
            "The q=2tL^2<1 multiplicity-compatible parabolic/Hermite "
            "first-jet chart and its finite connectors remain separate."
        ),
        "proof_boundary": (
            "The all-fiber shear identity, determinant-one homotopy, "
            "positive-scale invariance, reduced argument flux, crossing "
            "orientation handoff, zero-fibre inclusion, and generic "
            "backward-heat many-crossing guard are exact or certified. "
            "The pointwise Xi Abel-scalar gap and the Xi phase budget "
            "below one turn are independent open theorems. No q<1 "
            "closure, complete boundary composition, contact exclusion, "
            "Lambda<=0, PF-infinity, RH proof, or Clay-prize conclusion "
            "is asserted."
        ),
        "diagnostics": {
            "shear_flux": shear,
            "orientation": orientation,
            "heat_recrossing": heat,
        },
    }


def build_rows(exact: dict) -> list[FluxRow]:
    return [
        FluxRow(
            "assfr_00_pi_provenance",
            "source_provenance",
            "available_exact",
            "No unexplained pi is introduced.",
            exact["pi_provenance"],
            "The generic sine guard is not an Xi normalization.",
        ),
        FluxRow(
            "assfr_01_slope_shear",
            "exact_identity",
            "available_exact",
            "The physical slope is a real shear of the Abel scalar.",
            exact["slope_shear_identity"],
            "Valid inside one fixed-N q>=1 chart.",
        ),
        FluxRow(
            "assfr_02_gap_nonvanishing",
            "conditional_exact_reduction",
            "available_exact",
            "The proposed gap makes the Abel pair nonzero.",
            exact["pointwise_gap_nonvanishing"],
            "Conditional on the still-open Xi pointwise gap.",
        ),
        FluxRow(
            "assfr_03_shear_homotopy",
            "exact_homotopy",
            "available_exact",
            "The terminal slope shear can be removed topologically.",
            exact["shear_homotopy"],
            "Requires nonvanishing of the Abel pair.",
            exact["diagnostics"]["shear_flux"],
        ),
        FluxRow(
            "assfr_04_positive_scale",
            "exact_homotopy",
            "available_exact",
            "A positive jet scale does not change winding.",
            exact["positive_scale_homotopy"],
            "The scale must remain strictly positive.",
        ),
        FluxRow(
            "assfr_05_reduced_proxy",
            "exact_reduction",
            "available_exact",
            "The q>=1 integer reduces to the Abel proxy.",
            exact["reduced_proxy"],
            "Cutoff and finite connector homotopies remain required.",
        ),
        FluxRow(
            "assfr_06_general_flux",
            "exact_identity",
            "available_exact",
            "The Abel proxy has an exact division-free phase current.",
            exact["reduced_derivative"] + " " + exact["general_flux"],
            "Applied only where the proxy is nonzero.",
        ),
        FluxRow(
            "assfr_07_crossing_speed",
            "exact_identity",
            "available_exact",
            "The local ray-intersection speed is explicit.",
            exact["crossing_speed"],
            "The scalar must be nonzero at the crossing.",
        ),
        FluxRow(
            "assfr_08_orientation_handoff",
            "certified_bound",
            "finite_validated",
            "The pointwise gap fixes the true crossing orientation.",
            exact["orientation_error"]
            + " "
            + exact["orientation_ratio"]
            + " "
            + exact["orientation_consequence"],
            "This does not bound the number of crossings.",
            exact["diagnostics"]["orientation"],
        ),
        FluxRow(
            "assfr_09_zero_fibre",
            "exact_identity",
            "available_exact",
            "The reduced proxy includes W_0=0 without division.",
            exact["zero_fibre"],
            "No Wronskian or ordinary ratio is substituted.",
        ),
        FluxRow(
            "assfr_10_chart_boundary",
            "nonpromotion_guard",
            "guard_validated",
            "The Abel scalar is not promoted to a complex chart invariant.",
            exact["chart_boundary"],
            "Adjacent real-projection homotopies remain explicit.",
        ),
        FluxRow(
            "assfr_11_heat_family",
            "countermodel",
            "guard_validated",
            "A backward-heat family has arbitrarily many crossings.",
            exact["heat_family"],
            "Generic heat-flow family, not Xi.",
            exact["diagnostics"]["heat_recrossing"],
        ),
        FluxRow(
            "assfr_12_heat_gap",
            "countermodel_property",
            "guard_validated",
            "The same family satisfies an arbitrarily strong local gap.",
            exact["heat_gap"],
            "The frequency is unrestricted.",
        ),
        FluxRow(
            "assfr_13_heat_winding",
            "countermodel_property",
            "guard_validated",
            "Its first-jet winding is unbounded.",
            exact["heat_winding"],
            "This rejects implication from pointwise gap to one-turn flux.",
        ),
        FluxRow(
            "assfr_14_nonpromotion",
            "countermodel_consequence",
            "guard_validated",
            "Pointwise transversality and signed flux are independent.",
            exact["heat_nonpromotion"],
            "No claim about the actual Xi crossing count is made.",
        ),
        FluxRow(
            "assfr_15_pointwise_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The Xi Abel-scalar contact-band gap remains open.",
            exact["q_ge_1_pointwise_target"],
            "No generic coefficient family can prove it.",
        ),
        FluxRow(
            "assfr_16_flux_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The Xi Abel-proxy one-turn phase budget remains open.",
            exact["q_ge_1_flux_target"],
            exact["minimal_route_verdict"],
        ),
        FluxRow(
            "assfr_17_q_lt_1_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The multiplicity-compatible small-q route remains separate.",
            exact["q_lt_1_target"],
            "No small-q or connector closure follows.",
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    rows = build_rows(exact)
    return {
        "kind": STEM,
        "date": "2026-07-28",
        "status": (
            "exact endpoint-complete Abel-scalar shear/scale homotopy, "
            "reduced physical phase flux, certified crossing-orientation "
            "handoff, and backward-heat many-crossing guard; the Xi "
            "pointwise gap and one-turn flux remain independent open theorems"
        ),
        "proof_boundary": exact["proof_boundary"],
        "sources": {
            key: str(path.relative_to(REPO_ROOT))
            for key, path in SOURCES.items()
        },
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return "\n".join(
        [
            "# Newman Abel-Scalar Shear-Flux Reduction",
            "",
            "Date: 2026-07-28",
            "",
            "Status: exact linked-point topological reduction and route guard;",
            "not a proof of the Xi pointwise gap, one-turn phase budget,",
            "`Lambda<=0`, PF-infinity, RH, or a Clay-prize result.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Pi Provenance",
            "",
            exact["pi_provenance"],
            "",
            "## Exact Shear",
            "",
            "```text",
            exact["slope_shear_identity"],
            exact["pointwise_gap_nonvanishing"],
            exact["shear_homotopy"],
            exact["shear_matrix"],
            exact["positive_scale_homotopy"],
            "```",
            "",
            "## Reduced Proxy",
            "",
            "```text",
            exact["reduced_proxy"],
            exact["reduced_derivative"],
            exact["general_flux"],
            exact["crossing_speed"],
            "```",
            "",
            "## Crossing Orientation",
            "",
            exact["orientation_error"],
            "",
            exact["orientation_ratio"],
            "",
            exact["orientation_consequence"],
            "",
            "## Zero Fibre",
            "",
            exact["zero_fibre"],
            "",
            exact["chart_boundary"],
            "",
            "## Backward-Heat Guard",
            "",
            "```text",
            exact["heat_family"],
            exact["heat_gap"],
            exact["heat_winding"],
            "```",
            "",
            exact["heat_nonpromotion"],
            "",
            "## Minimal Linked-Point Route",
            "",
            "```text",
            exact["q_ge_1_pointwise_target"],
            exact["q_ge_1_flux_target"],
            exact["q_lt_1_target"],
            "```",
            "",
            exact["minimal_route_verdict"],
            "",
            "## Boundary",
            "",
            exact["proof_boundary"],
            "",
        ]
    )


def write_outputs(artifact: dict, output: Path, note: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    note.write_text(render_note(artifact), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    write_outputs(artifact, args.output, args.note)
    print(
        "built Newman Abel-scalar shear-flux reduction: "
        "18 rows, 1 exact shear/scale homotopy, "
        "1 reduced O(N) physical flux, "
        "1 certified crossing-orientation handoff, "
        "1 backward-heat many-crossing guard, 3 open obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
