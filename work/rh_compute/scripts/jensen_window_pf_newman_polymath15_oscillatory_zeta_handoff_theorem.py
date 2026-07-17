#!/usr/bin/env python3
"""Build the oscillatory zeta handoff for the Newman scaled boundary layer."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem.md"
)
POLYMATH_SOURCE_URL = "https://arxiv.org/abs/1904.12438"
EXPONENT_PAIR_SOURCE_URL = "https://arxiv.org/abs/2306.05599"
BOURGAIN_SOURCE_URL = "https://arxiv.org/abs/1408.5794"


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def a_process(pair: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    kappa, lam = pair
    denominator = 2 * kappa + 2
    return kappa / denominator, (kappa + lam + 1) / denominator


def exponent_pair(index: int) -> tuple[Fraction, Fraction]:
    fixed = {
        0: (Fraction(13, 84), Fraction(55, 84)),
        1: (Fraction(4742, 38463), Fraction(35731, 51284)),
        2: (Fraction(18, 199), Fraction(593, 796)),
        3: (Fraction(2779, 38033), Fraction(58699, 76066)),
        4: (Fraction(715, 10238), Fraction(7955, 10238)),
    }
    if index in fixed:
        return fixed[index]
    if 5 <= index <= 8:
        return a_process(fixed[index - 4])
    if index >= 9:
        m = index - 4
        return (
            Fraction(2, (m - 1) ** 2 * (m + 2)),
            1 - Fraction(3 * m - 2, m * (m - 1) * (m + 2)),
        )
    raise ValueError(f"unsupported exponent-pair index: {index}")


def phase_exponent(index: int, radius: Fraction) -> Fraction:
    kappa, lam = exponent_pair(index)
    return 2 * kappa + (lam - kappa) * radius


def weighted_exponent(index: int, radius: Fraction, scaled_time: Fraction) -> Fraction:
    return (
        phase_exponent(index, radius)
        - (Fraction(1, 2) + scaled_time / 4) * radius
        + scaled_time * radius**2 / 8
    )


def transition_radius(upper_index: int, lower_index: int) -> Fraction:
    upper_kappa, upper_lam = exponent_pair(upper_index)
    lower_kappa, lower_lam = exponent_pair(lower_index)
    numerator = 2 * (lower_kappa - upper_kappa)
    denominator = (
        upper_lam
        - upper_kappa
        - lower_lam
        + lower_kappa
    )
    return numerator / denominator


def required_scaled_time(index: int, radius: Fraction) -> Fraction:
    numerator = 8 * (phase_exponent(index, radius) - radius / 2)
    return numerator / (radius * (2 - radius))


def fraction_record(value: Fraction) -> dict:
    return {
        "exact": f"{value.numerator}/{value.denominator}",
        "decimal": f"{float(value):.15f}",
    }


def build_envelope() -> dict:
    transition_rows: list[dict] = []
    for upper_index in range(10, 0, -1):
        lower_index = upper_index - 1
        radius = transition_radius(upper_index, lower_index)
        upper_phase = phase_exponent(upper_index, radius)
        lower_phase = phase_exponent(lower_index, radius)
        if upper_phase != lower_phase:
            raise RuntimeError("exponent-pair transition identity failed")
        threshold = required_scaled_time(upper_index, radius)
        transition_rows.append(
            {
                "upper_index": upper_index,
                "lower_index": lower_index,
                "radius": fraction_record(radius),
                "phase_exponent": fraction_record(upper_phase),
                "required_scaled_time": fraction_record(threshold),
            }
        )

    critical_transition = next(
        row
        for row in transition_rows
        if row["upper_index"] == 2 and row["lower_index"] == 1
    )
    critical_radius = transition_radius(2, 1)
    critical_scaled_time = required_scaled_time(2, critical_radius)
    if critical_transition["required_scaled_time"]["exact"] != (
        f"{critical_scaled_time.numerator}/{critical_scaled_time.denominator}"
    ):
        raise RuntimeError("critical transition serialization failed")

    absolute_radius = 2 - 4 / critical_scaled_time
    boundaries = [absolute_radius]
    boundaries.extend(
        transition_radius(index, index - 1) for index in range(10, 0, -1)
    )
    boundaries.append(Fraction(1, 1))
    if boundaries != sorted(boundaries):
        raise RuntimeError("exponent-pair envelope boundaries are not ordered")

    interval_rows: list[dict] = []
    finite_indices = tuple(range(11))
    for position, index in enumerate(range(10, -1, -1)):
        left = boundaries[position]
        right = boundaries[position + 1]
        midpoint = (left + right) / 2
        for test_radius in (left, midpoint, right):
            selected = phase_exponent(index, test_radius)
            if selected != min(
                phase_exponent(other, test_radius) for other in finite_indices
            ):
                raise RuntimeError(
                    f"pair {index} is not active at radius {test_radius}"
                )
        kappa, lam = exponent_pair(index)
        derivative_offset = lam - kappa - Fraction(1, 2)
        if derivative_offset < 0:
            raise RuntimeError("required-time endpoint argument lost nonnegativity")
        left_required = required_scaled_time(index, left)
        right_required = required_scaled_time(index, right)
        interval_rows.append(
            {
                "index": index,
                "kappa": fraction_record(kappa),
                "lambda": fraction_record(lam),
                "left_radius": fraction_record(left),
                "right_radius": fraction_record(right),
                "left_required_scaled_time": fraction_record(left_required),
                "right_required_scaled_time": fraction_record(right_required),
                "derivative_offset": fraction_record(derivative_offset),
            }
        )

    endpoint_thresholds = [
        Fraction(row["left_required_scaled_time"]["exact"])
        for row in interval_rows
    ] + [
        Fraction(row["right_required_scaled_time"]["exact"])
        for row in interval_rows
    ]
    if max(endpoint_thresholds) != critical_scaled_time:
        raise RuntimeError("critical scaled threshold is not the envelope maximum")
    if sum(value == critical_scaled_time for value in endpoint_thresholds) != 2:
        raise RuntimeError("critical threshold should occur on both sides of one transition")

    for row in transition_rows:
        row["is_critical_maximum"] = (
            Fraction(row["required_scaled_time"]["exact"])
            == critical_scaled_time
        )

    return {
        "pair_indices": list(range(10, -1, -1)),
        "transition_rows": transition_rows,
        "interval_rows": interval_rows,
        "critical_radius": fraction_record(critical_radius),
        "critical_scaled_time": fraction_record(critical_scaled_time),
        "absolute_split_radius_at_threshold": fraction_record(absolute_radius),
        "endpoint_maximum_argument": (
            "For G(r)=a+b*r, d=b-1/2>=0, the derivative numerator of "
            "8*(G(r)-r/2)/(r*(2-r)) is d*r^2+2*a*r-2*a. It is "
            "strictly increasing, so any interior critical point is a minimum; "
            "the maximum on each active interval is at an endpoint."
        ),
    }


def build_exact() -> dict:
    envelope = build_envelope()
    amplitude, u, u1, v, v1, theta, potential = sp.symbols(
        "amplitude u u1 v v1 theta potential", real=True
    )
    value = 2 * amplitude * sp.cos(theta)
    first = 2 * amplitude * (u * sp.cos(theta) - v * sp.sin(theta))
    second = 2 * amplitude * (
        (u**2 + u1 - v**2) * sp.cos(theta)
        - (2 * u * v + v1) * sp.sin(theta)
    )
    curvature = sp.trigsimp(first**2 - value * second + potential * value**2)
    expected = 4 * amplitude**2 * (
        v**2
        + (potential - u1) * sp.cos(theta) ** 2
        + v1 * sp.sin(2 * theta) / 2
    )
    if sp.trigsimp(sp.expand_trig(curvature - expected)) != 0:
        raise RuntimeError("phase-amplitude curvature identity failed")

    c_star = envelope["critical_scaled_time"]["exact"]
    return {
        "scale": (
            "L=log(x/(4*pi)), c=t*L, N=floor(sqrt(exp(L)+t/16)), "
            "L=2*log(N)+o(1), sigma=Re(s_*)=1/2+c/4+o(1), "
            "tau=|Im(s_*)|=2*pi*N^2*(1+o(1))"
        ),
        "moments": (
            "D_k=sum_(n<=N)exp((t/4)log(n)^2)log(n)^k*n^(-s_*), "
            "Z_k=sum_(n>=1)log(n)^k*n^(-s_*), k=0,1,2"
        ),
        "phase_bound": (
            "For M=N^r and exponent pair (kappa,lambda), "
            "sup_I|sum_(n in I)n^(-i*tau)| "
            "<<_eta N^(2*kappa+(lambda-kappa)*r+eta)"
        ),
        "weighted_exponent": (
            "B_(kappa,lambda,c)(r)=2*kappa+(lambda-kappa)*r-"
            "(1/2+c/4)*r+c*r^2/8"
        ),
        "critical_threshold": (
            f"c_*={c_star}=2.540223984760008...; it occurs at "
            "r_*=125662/155153 where exponent-pair rows 2 and 1 cross"
        ),
        "low_range": (
            "For c>=c_*+epsilon choose "
            "2-4/c_*<rho<2-4/(c_*+epsilon). On n<=N^rho, "
            "exp((t/4)log(n)^2)n^(-sigma)<=n^(-1-delta_epsilon), "
            "so sum (w_n-1)log(n)^k*n^(-s_*)=O_epsilon(1/L)."
        ),
        "high_range": (
            "On N^rho<n<=N, dyadic Abel summation and the 11-pair envelope "
            "give sum (w_n-1)log(n)^k*n^(-s_*)=O_epsilon(N^(-delta_epsilon))."
        ),
        "abel_variation": (
            "A_k(u)=(exp(q)-1)u^(-sigma)log(u)^k, q=t*log(u)^2/4; "
            "|A_k'(u)|<=u^(-sigma-1){exp(q)*(t/2)*log(u)^(k+1)+"
            "sigma*(exp(q)-1)*log(u)^k+k*(exp(q)-1)*log(u)^(k-1)}. "
            "Since t*log(u)=O(1) on u<=N, one dyadic total variation is "
            "O_epsilon(log(N)^(k+1)*max(exp(q)u^(-sigma)))."
        ),
        "tail": (
            "sum_(n>N)log(n)^k*n^(-s_*)=O_epsilon(N^(1-sigma)log(N)^k)=o_epsilon(1)"
        ),
        "jet_handoff": (
            "D_k-Z_k=o_epsilon(1), k=0,1,2, uniformly for "
            "c_*+epsilon<=c<=25; hence D and its first two x-jets converge "
            "to zeta(s_*) and its corresponding x-jets"
        ),
        "zeta_floor": (
            "sigma>1 and the Euler product give |zeta(s_*)|>=1/zeta(sigma_0)>0 "
            "for a fixed sigma_0>1 depending only on epsilon"
        ),
        "phase_amplitude": (
            "For D=a*exp(i*psi), theta=beta+psi, u=(log(a))', v=theta': "
            "C[P]=4*a^2*(v^2+(V-u')cos(theta)^2+(v'/2)sin(2*theta))"
        ),
        "phase_geometry": (
            "a>=m_epsilon>0, u'=O_epsilon(1), "
            "v=-L/4+O_epsilon(1), v'=O_epsilon(1), V=o(1)"
        ),
        "main_curvature": (
            "C_t[P_(N,t)]>=c_epsilon*L^2 for all sufficiently large L "
            "when c_*+epsilon<=t*L<=25"
        ),
        "remainder_transfer": (
            "On radius-1/L cutoff collars, the normalized Polymath-15 A+B "
            "remainder plus one adjacent-N block is "
            "E=O_epsilon(N^(-1/2-c_*/8+o(1))); its j-th jet is "
            "O_epsilon(L^j*E), j=0,1,2, and its curvature cost is o_epsilon(L^2)."
        ),
        "remainder_powers": (
            "sum_(n<=N)w_n*n^(-sigma)<=N^(1/2-c/8+o(1)); "
            "the e_A+e_B bracket is O(N^-2), so e_A+e_B="
            "O(N^(-3/2-c/8+o(1))); e_C0 and one adjacent-N block are "
            "O(N^(-1/2-c/8+o(1)))."
        ),
        "theorem": (
            f"For every epsilon>0 there exists L_epsilon such that 0<t<=1/2, "
            f"L>=L_epsilon, and t*L>={c_star}+epsilon imply L_t(x)>0"
        ),
        "critical_layer": f"0<t*log(x/(4*pi))<={c_star}+o(1)",
        "envelope": envelope,
    }


def build_artifact() -> dict:
    exact = build_exact()
    envelope = exact["envelope"]
    rows = [
        GateRow(
            id="np15ozht_01_scaled_geometry",
            role="exact_asymptotic_reduction",
            readiness="ready_to_apply",
            claim="The finite saddle sum has length N while its logarithmic phase has frequency asymptotic to 2*pi*N^2.",
            formula=exact["scale"],
            proof_boundary="Uniform for bounded c=tL and the real-axis Polymath-15 exponent s_*.",
        ),
        GateRow(
            id="np15ozht_02_exponent_pair_input",
            role="published_analytic_input",
            readiness="available_published",
            claim="Eleven published exponent pairs provide scale-dependent cancellation for n^(-i*tau).",
            formula=exact["phase_bound"],
            proof_boundary="Uses the exponent-pair definition and convex-hull vertices in Trudgian-Yang; epsilon losses are absorbed by a strict c margin.",
            diagnostics=envelope["interval_rows"],
        ),
        GateRow(
            id="np15ozht_03_exact_envelope",
            role="exact_rational_certificate",
            readiness="ready_to_apply",
            claim="The finite scale-dependent exponent-pair envelope has one exact worst transition.",
            formula=f"{exact['weighted_exponent']}; {exact['critical_threshold']}",
            proof_boundary="All transition radii and endpoint maxima are checked with exact rational arithmetic.",
            diagnostics=envelope["transition_rows"],
        ),
        GateRow(
            id="np15ozht_04_low_absolute_range",
            role="proved_asymptotic_estimate",
            readiness="ready_to_apply",
            claim="A short initial range still admits the original absolute O(1/L) comparison.",
            formula=exact["low_range"],
            proof_boundary="The split rho depends on epsilon but is fixed as L tends to infinity.",
        ),
        GateRow(
            id="np15ozht_05_high_oscillatory_range",
            role="proved_asymptotic_estimate",
            readiness="ready_to_apply",
            claim="Dyadic partial summation makes the weighted perturbation vanish throughout the complementary range.",
            formula=exact["high_range"],
            proof_boundary="Uses bounded total variation of the smooth weight on each dyadic block and O(log N) blocks.",
        ),
        GateRow(
            id="np15ozht_06_zeta_jet_handoff",
            role="proved_asymptotic_theorem",
            readiness="ready_to_apply",
            claim="The weighted finite Dirichlet moments converge through order two to the ordinary zeta moments.",
            formula=f"{exact['jet_handoff']}; {exact['tail']}",
            proof_boundary="Uniform only with a fixed positive margin above c_* and bounded c.",
        ),
        GateRow(
            id="np15ozht_07_euler_floor",
            role="analytic_theorem",
            readiness="ready_to_apply",
            claim="The limiting zeta amplitude is uniformly separated from zero without assuming RH.",
            formula=exact["zeta_floor"],
            proof_boundary="Uses only absolute convergence of the Euler product in Re(s)>1.",
        ),
        GateRow(
            id="np15ozht_08_phase_curvature_identity",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The nonvanishing complex amplitude converts finite curvature into a dominant phase-speed square and bounded corrections.",
            formula=exact["phase_amplitude"],
            proof_boundary="Exact symbolic identity; no division is used before the Euler floor is established.",
        ),
        GateRow(
            id="np15ozht_09_main_curvature",
            role="proved_asymptotic_theorem",
            readiness="ready_to_apply",
            claim="The finite real Riemann-Siegel main sum has positive normalized curvature above c_*.",
            formula=f"{exact['phase_geometry']}; {exact['main_curvature']}",
            proof_boundary="The threshold L_epsilon is existential and is not numerically optimized.",
        ),
        GateRow(
            id="np15ozht_10_remainder_transfer",
            role="proved_asymptotic_transfer",
            readiness="ready_to_apply",
            claim="The published analytic remainder and cutoff jump are lower order in normalized C2 curvature.",
            formula=exact["remainder_transfer"],
            proof_boundary="Uses the Polymath-15 A+B estimate, radius-1/L Cauchy transfer, and the adjacent-cutoff lift already audited locally.",
        ),
        GateRow(
            id="np15ozht_11_exact_h_theorem",
            role="proved_theorem",
            readiness="ready_to_apply",
            claim="Exact Newman first-Laguerre positivity holds asymptotically above the oscillatory threshold.",
            formula=exact["theorem"],
            proof_boundary="The bounded branch c<=25 is combined with the existing exact c>=25 global ray.",
        ),
        GateRow(
            id="np15ozht_12_live_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The asymptotically unresolved RH-level layer is reduced to the scaled strip at and below c_*.",
            formula=exact["critical_layer"],
            proof_boundary="Open; this artifact is not Lambda<=0 or RH.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem",
        "date": "2026-07-17",
        "status": (
            "quantified oscillatory asymptotic strict Laguerre theorem for exact H_t "
            "above c_*=4911678521/1933561194; not a proof of Lambda<=0 or RH"
        ),
        "proof_boundary": (
            "For every fixed epsilon>0 this artifact proves the existence of an "
            "L_epsilon beyond which L_t(x)>0 whenever 0<t<=1/2 and "
            "t*log(x/(4*pi))>=4911678521/1933561194+epsilon. It does not provide "
            "a practical L_epsilon, cover the scaled layer at or below that threshold, "
            "prove positivity for all x and t, prove Lambda<=0, or prove RH."
        ),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "sources": [
            POLYMATH_SOURCE_URL,
            EXPONENT_PAIR_SOURCE_URL,
            BOURGAIN_SOURCE_URL,
            "outputs/jensen_window_pf_newman_polymath15_normalized_laguerre_bridge.md",
            "outputs/jensen_window_pf_newman_polymath15_scaled_boundary_layer_asymptotic_theorem.md",
            "outputs/jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    envelope = exact["envelope"]
    lines = [
        "# Jensen-Window PF Newman Polymath-15 Oscillatory Zeta Handoff Theorem",
        "",
        "Date: 2026-07-17",
        "",
        "Status: quantified asymptotic strict first-Laguerre positivity for the",
        "exact Newman heat flow above an explicit oscillatory threshold. This is",
        "not a proof of `Lambda <= 0` or RH.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem.json",
        "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem.py",
        "```",
        "",
        "The finite heat-flow approximation is imported from D. H. J. Polymath,",
        f"[Effective approximation of heat flow evolution of the Riemann xi function]({POLYMATH_SOURCE_URL}).",
        "The cancellation input is the exponent-pair definition and explicit",
        "known-pair hull in Trudgian and Yang,",
        f"[Toward optimal exponent pairs]({EXPONENT_PAIR_SOURCE_URL}); the central",
        f"pair is originally due to [Bourgain]({BOURGAIN_SOURCE_URL}).",
        "",
        "## Scaled Geometry",
        "",
        "```text",
        exact["scale"],
        exact["moments"],
        "```",
        "",
        "The ordinary zeta tail is already absolutely small for `c>2`. The new",
        "task is to show that the heat weight minus one contributes `o(1)`.",
        "",
        "## Exponent-Pair Envelope",
        "",
        "On a dyadic block `M=N^r`, the published exponent-pair definition gives",
        "",
        "```text",
        exact["phase_bound"],
        exact["weighted_exponent"],
        "```",
        "",
        "The exponent pair may change with `r`. Exact rational optimization of",
        "eleven available pairs gives",
        "",
        "```text",
        exact["critical_threshold"],
        "```",
        "",
        "| pair change | r | required c | maximum |",
        "|---:|---:|---:|:---:|",
    ]
    for row in envelope["transition_rows"]:
        lines.append(
            "| {upper}->{lower} | {radius} | {threshold} | {maximum} |".format(
                upper=row["upper_index"],
                lower=row["lower_index"],
                radius=row["radius"]["decimal"],
                threshold=row["required_scaled_time"]["decimal"],
                maximum="yes" if row["is_critical_maximum"] else "",
            )
        )
    lines.extend(
        [
            "",
            envelope["endpoint_maximum_argument"],
            "Thus the table is a proof certificate for the full continuous scale",
            "interval, not a sampled optimization.",
            "",
            "## Zeta Handoff",
            "",
            "Fix `epsilon>0`. Choose the low/high split just below the absolute",
            "summability endpoint. Then",
            "",
            "```text",
            exact["low_range"],
            exact["high_range"],
            exact["abel_variation"],
            exact["tail"],
            "```",
            "",
            "Partial summation costs only the total variation of the smooth heat",
            "weight on each dyadic block. Its logarithmic derivative is bounded",
            "for bounded `c`, logarithmic moment factors are absorbed in the strict",
            "power margin, and there are only `O(log N)` blocks. Consequently",
            "",
            "```text",
            exact["jet_handoff"],
            exact["zeta_floor"],
            "```",
            "",
            "## Curvature Transfer",
            "",
            "Writing `D=a*exp(i*psi)` and `theta=beta+psi`, exact differentiation",
            "gives",
            "",
            "```text",
            exact["phase_amplitude"],
            exact["phase_geometry"],
            exact["main_curvature"],
            "```",
            "",
            "The Euler-product floor makes division by `D` legitimate for all",
            "sufficiently large `L`. The `v^2` term is of order `L^2`; all possible",
            "negative logarithmic corrections remain bounded.",
            "",
            "On radius-`1/L` cutoff collars the published scalar remainder and one",
            "adjacent cutoff block satisfy",
            "",
            "```text",
            exact["remainder_powers"],
            exact["remainder_transfer"],
            "```",
            "",
            "The exact normalized-curvature perturbation identity therefore",
            "preserves the main sign. Combining with the existing `c>=25` ray gives",
            "",
            "```text",
            exact["theorem"],
            "```",
            "",
            "## Remaining Gap",
            "",
            "The asymptotically unresolved layer is now",
            "",
            "```text",
            exact["critical_layer"],
            "```",
            "",
            "The improvement uses established cancellation estimates; it is not a",
            "finite-height verification and does not prove `Lambda <= 0` or RH.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman Polymath-15 oscillatory zeta handoff theorem: "
        f"{len(artifact['rows'])} rows, 11 exponent pairs, "
        "critical scaled threshold 4911678521/1933561194"
    )


if __name__ == "__main__":
    main()
