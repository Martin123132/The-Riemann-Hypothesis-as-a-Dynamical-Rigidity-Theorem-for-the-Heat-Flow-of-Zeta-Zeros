#!/usr/bin/env python3
"""Build the Hermitian reciprocal swap and near-diagonal reduction."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "contiguous_terminal_tail_anchored_six_moment_hermitian_"
    "reciprocal_pairing_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "flow_matrix_phase": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_flow_matrix_"
        "phase_reduction.json"
    ),
    "saddle_flow": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_saddle_flow_"
        "reduction.json"
    ),
    "two_carrier": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_moment_two_carrier_"
        "kernel_reduction.json"
    ),
    "reciprocal_saddle": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "reciprocal_saddle_self_duality_gate.json"
    ),
}


@dataclass(frozen=True)
class PairingRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def require_zero(expression: sp.Expr, label: str) -> None:
    if sp.simplify(sp.expand_trig(sp.expand_complex(expression))) != 0:
        raise RuntimeError(label)


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCE_PATHS.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCE_PATHS.items()
    }


def source_audit(payloads: dict[str, dict]) -> dict:
    phase = payloads["flow_matrix_phase"]
    if phase.get("counts", {}).get("hermitian_reciprocal_diagonal_nulls") != 1:
        raise RuntimeError("Hermitian reciprocal source drifted")
    if phase.get("counts", {}).get("signed_flow_bounds") != 0:
        raise RuntimeError("flow proof boundary drifted")

    flow = payloads["saddle_flow"]
    leading = flow.get("leading_pair_flow_certificate", {})
    if "H_0^+=-Sigma_nm^2/8" not in leading.get("kernel_origin", ""):
        raise RuntimeError("leading Hermitian flow source drifted")

    carrier = payloads["two_carrier"].get("symbolic_certificate", {})
    if "H^+_(n,m)" not in carrier.get("hermitian_symmetrization", ""):
        raise RuntimeError("two-carrier Hermitian source drifted")
    if "H^+_(n,m)=-(u_n+u_m)^2/8" not in carrier.get(
        "leading_q1_kernels", ""
    ):
        raise RuntimeError("leading real-kernel source drifted")

    reciprocal = payloads["reciprocal_saddle"].get("exact", {})
    if "a_omega^2=omega/(2*pi)" not in reciprocal.get("pi_provenance", ""):
        raise RuntimeError("reciprocal pi provenance drifted")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCE_PATHS.items()
        },
        "source_kinds": {
            key: payload.get("kind", "") for key, payload in payloads.items()
        },
        "imported_poisson_remainder_bounds": 0,
        "imported_signed_flow_bounds": 0,
    }


def swap_pair_certificate() -> dict:
    alpha, mu, real_part, imag_part = sp.symbols(
        "alpha mu R I", real=True
    )
    angle = 2 * sp.pi * alpha * mu
    carrier = (real_part + sp.I * imag_part) * (
        sp.cos(angle) + sp.I * sp.sin(angle)
    )
    reverse_carrier = sp.conjugate(carrier)
    ordered_flow = sp.I * mu * carrier
    reverse_flow = -sp.I * mu * reverse_carrier
    paired_flow = sp.expand_complex(ordered_flow + reverse_flow)
    expected_flow = -2 * mu * (
        real_part * sp.sin(angle) + imag_part * sp.cos(angle)
    )
    require_zero(paired_flow - expected_flow, "Hermitian swap pairing failed")

    undifferentiated_pair = sp.expand_complex(carrier + reverse_carrier)
    phase_derivative = sp.diff(undifferentiated_pair, alpha) / (2 * sp.pi)
    require_zero(
        phase_derivative - expected_flow,
        "Hermitian phase derivative identity failed",
    )
    require_zero(expected_flow.subs(mu, 0), "Hermitian diagonal null failed")

    outer_half_pair = sp.simplify(expected_flow / 2)
    leading_real_pair = sp.simplify(outer_half_pair.subs(imag_part, 0))
    correction_pair = sp.simplify(outer_half_pair - leading_real_pair)

    negative_witness = sp.simplify(
        outer_half_pair.subs(
            {alpha: 1, mu: sp.Rational(1, 4), real_part: 1, imag_part: 0}
        )
    )
    positive_witness = sp.simplify(
        outer_half_pair.subs(
            {alpha: 3, mu: sp.Rational(1, 4), real_part: 1, imag_part: 0}
        )
    )
    resonance_witness = sp.simplify(
        outer_half_pair.subs(
            {alpha: 2, mu: sp.Rational(1, 2), real_part: 0, imag_part: 1}
        )
    )
    if (
        negative_witness != -sp.Rational(1, 4)
        or positive_witness != sp.Rational(1, 4)
        or resonance_witness != -sp.Rational(1, 2)
    ):
        raise RuntimeError("Hermitian pairing sign guards failed")

    return {
        "swap_covariance": (
            "Let Z_(k,l)=A_(k,l)e(alpha*mu_(k,l)), where "
            "e(t)=exp(2*pi*i*t), mu_(k,l)=log(l/k), and the interior "
            "Hermitian stationary coefficients obey "
            "A_(l,k)=conj(A_(k,l)). Then Z_(l,k)=conj(Z_(k,l)) and "
            "mu_(l,k)=-mu_(k,l)."
        ),
        "ordered_pair": (
            "Writing A_(k,l)=R_(k,l)+iI_(k,l), the two ordered flow "
            "terms satisfy i*mu*Z_(k,l)-i*mu*conj(Z_(k,l))="
            "-2mu{R sin(2*pi*alpha*mu)+I cos(2*pi*alpha*mu)}."
        ),
        "outer_half_pair": (
            "After the existing outer factor 1/2, one unordered pair "
            "k<l contributes P_H(k,l)=-mu{R sin(2*pi*alpha*mu)+"
            "I cos(2*pi*alpha*mu)}."
        ),
        "phase_derivative": (
            "The ordered flow pair is exactly (1/(2*pi))*partial_alpha "
            "[Z_(k,l)+conj(Z_(k,l))] when the stationary coefficient is "
            "held fixed. Equivalently, since xi=2*pi*alpha, it is the "
            "xi phase derivative."
        ),
        "real_defect_split": (
            "The real channel is -mu*R*sin(2*pi*alpha*mu), while the "
            "imaginary correction is -mu*I*cos(2*pi*alpha*mu). Thus the "
            "real leading kernel gains a sine zero at integral phase, but "
            "an imaginary correction generally gains only the single "
            "mu factor."
        ),
        "diagonal_null": (
            "At k=l one has mu=0, so the full Hermitian differentiated "
            "stationary pair vanishes, including every complex correction."
        ),
        "sign_and_resonance_guards": (
            "For R=1,I=0,mu=1/4 the outer-half pair is -1/4 at alpha=1 "
            "and +1/4 at alpha=3, so swap symmetry has no sign. For "
            "R=0,I=1,mu=1/2,alpha=2 the phase is integral but the pair is "
            "-1/2, so phase resonance does not kill the imaginary defect."
        ),
        "symbolic": {
            "paired_ordered_flow": str(sp.simplify(paired_flow)),
            "outer_half_pair": str(outer_half_pair),
            "real_channel": str(leading_real_pair),
            "imaginary_channel": str(correction_pair),
            "phase_derivative_difference": "0",
            "negative_real_witness": str(negative_witness),
            "positive_real_witness": str(positive_witness),
            "imaginary_resonance_witness": str(resonance_witness),
        },
    }


def reciprocal_leading_certificate() -> dict:
    a, k, ell, weight_k, weight_l = sp.symbols(
        "a k l W_k W_l", positive=True, real=True
    )
    alpha = a**2
    n_star = alpha / k
    m_star = alpha / ell
    mu = sp.expand_log(sp.log(ell / k), force=True)
    tau = sp.expand_log(sp.log(k * ell / alpha), force=True)
    u_n = sp.expand_log(sp.log(a / n_star), force=True)
    u_m = sp.expand_log(sp.log(a / m_star), force=True)
    require_zero(u_n - u_m + mu, "reciprocal difference coordinate failed")
    require_zero(u_n + u_m - tau, "reciprocal sum coordinate failed")

    hessian = sp.diag(-k**2 / alpha, ell**2 / alpha)
    inverse_sqrt_determinant = alpha / (k * ell)
    require_zero(
        inverse_sqrt_determinant**2 * (-hessian.det()) - 1,
        "reciprocal stationary factor failed",
    )
    leading_kernel = -tau**2 / 8
    coefficient = sp.simplify(
        leading_kernel * inverse_sqrt_determinant * weight_k * weight_l
    )
    positive_scale = sp.simplify(-coefficient)
    leading_pair = sp.simplify(
        positive_scale * mu * sp.sin(2 * sp.pi * alpha * mu)
    )

    return {
        "reciprocal_coordinates": (
            "At n_*=a^2/k and m_*=a^2/l, put "
            "mu=log(l/k) and tau=log(kl/a^2). Then "
            "u_n-u_m=-mu and u_n+u_m=tau."
        ),
        "stationary_factor": (
            "The Hermitian Hessian is diag(-k^2/a^2,+l^2/a^2), has "
            "signature zero, and contributes the positive leading factor "
            "a^2/(kl); no Maslov phase is introduced by this mixed "
            "signature."
        ),
        "leading_coefficient": (
            "For the ideal q=1 kernel H_0^+=-tau^2/8 and positive real "
            "normalized saddle amplitudes W_k,W_l, the interior stationary "
            "coefficient is A^0_(k,l)=-C_(k,l), where "
            "C_(k,l)=a^2*tau^2*W_kW_l/(8kl)>=0."
        ),
        "leading_unordered_pair": (
            "The corresponding outer-half unordered Hermitian flow main is "
            "P_H^0(k,l)=C_(k,l)*mu*sin(2*pi*a^2*mu). It has the exact "
            "diagonal zero and a second zero whenever a^2*mu is integral, "
            "but it has no fixed sign between those phases."
        ),
        "symbolic": {
            "u_difference": str(sp.simplify(u_n - u_m)),
            "u_sum": str(sp.simplify(u_n + u_m)),
            "hessian": str(hessian),
            "inverse_sqrt_abs_det": str(inverse_sqrt_determinant),
            "leading_stationary_coefficient": str(coefficient),
            "leading_unordered_pair": str(leading_pair),
        },
    }


def offset_and_defect_certificate() -> dict:
    alpha, k, d = sp.symbols("alpha k d", positive=True, real=True)
    x = d / k
    mu = sp.log(1 + x)
    phase_cycles = alpha * mu
    t = sp.symbols("t", nonnegative=True, real=True)
    taylor_defect = x - sp.log(1 + x)
    integral_defect = sp.integrate(t / (1 + t), (t, 0, x))
    require_zero(
        taylor_defect - integral_defect,
        "logarithmic Taylor defect identity failed",
    )
    require_zero(
        sp.diff(phase_cycles, d) - alpha / (k + d),
        "offset phase first derivative failed",
    )
    require_zero(
        sp.diff(phase_cycles, d, 2) + alpha / (k + d) ** 2,
        "offset phase second derivative failed",
    )

    return {
        "offset_coordinates": (
            "For k<l write l=k+d with integer d>=1. Then "
            "mu_(k,d)=log(1+d/k), F_k(d)=alpha*mu_(k,d), "
            "F_k'(d)=alpha/(k+d), and F_k''(d)=-alpha/(k+d)^2. "
            "The two-dimensional ordered main is therefore reduced to "
            "one triangular family of one-dimensional logarithmic phases."
        ),
        "log_bounds": (
            "For x=d/k>=0, x-log(1+x)=integral_0^x t/(1+t)dt. Hence "
            "0<=x-log(1+x)<=x^2/2 and 0<=mu<=d/k."
        ),
        "lattice_defect": (
            "Let rho_k=dist(alpha/k,Z). Choosing an integer q with "
            "|alpha/k-q|=rho_k gives "
            "dist(alpha*log(1+d/k),Z) <= "
            "d*rho_k+alpha*d^2/(2k^2)."
        ),
        "proof_of_lattice_defect": (
            "The integer dq may be subtracted from the phase. The triangle "
            "inequality leaves d|alpha/k-q| plus "
            "alpha[d/k-log(1+d/k)], and the displayed logarithmic defect "
            "bounds the latter term."
        ),
        "sine_bound": (
            "Because |sin(2*pi*y)|<=2*pi*dist(y,Z), the real leading "
            "channel satisfies |P_H^0(k,k+d)| <= "
            "2*pi*C_(k,k+d)*(d/k)*"
            "[d*rho_k+alpha*d^2/(2k^2)]."
        ),
        "full_mixed_bound": (
            "For a full stationary coefficient R+iI, the exact swap formula "
            "gives |P_H(k,k+d)| <= "
            "2*pi*|R|*(d/k)*[d*rho_k+alpha*d^2/(2k^2)]"
            "+|I|*d/k. The first term exposes reciprocal-lattice "
            "suppression; the second is the correction obligation that "
            "cannot be discarded."
        ),
        "resonant_nonresonant_handoff": (
            "A proof may split by rho_k. Small rho_k is controlled by the "
            "displayed sine defect if R and I are bounded sharply. On the "
            "complement, the monotone derivative alpha/(k+d) identifies "
            "the remaining one-dimensional exponential-sum problem. This "
            "is a route decomposition, not a completed estimate."
        ),
        "symbolic": {
            "mu": str(mu),
            "phase_cycles": str(phase_cycles),
            "phase_first_derivative": str(sp.diff(phase_cycles, d)),
            "phase_second_derivative": str(sp.diff(phase_cycles, d, 2)),
            "taylor_defect_integral_difference": "0",
        },
    }


def cutoff_certificate() -> dict:
    mode_count = sp.symbols("M", integer=True, nonnegative=True)
    unordered = mode_count * (mode_count - 1) / 2
    require_zero(
        mode_count**2 - mode_count - 2 * unordered,
        "dual square pairing count failed",
    )
    return {
        "dual_square": (
            "Let K={ceil(a^2/B),...,floor(a^2)}. The interior two-variable "
            "stationary main is indexed by KxK, which is invariant under "
            "(k,l)<->(l,k). Its |K| diagonal modes vanish and its "
            "|K|(|K|-1) off-diagonal modes form exactly "
            "|K|(|K|-1)/2 unordered swap pairs."
        ),
        "boundary_guard": (
            "There is no unpaired mode in the square interior stationary "
            "main. This combinatorial fact does not pair the one-variable "
            "endpoint family, hard-cutoff half-weights, nonstationary modes, "
            "Poisson remainders, or the adjacent terminal recurrence."
        ),
        "symbolic_pairing_residual": "0",
    }


def exact_payload() -> dict:
    return {
        "domain": (
            "Use the physical q=1, L>=50, fixed-N terminal/bulk chart and "
            "the auxiliary segment from Sections 11.164-11.165. The present "
            "reduction concerns only the interior Hermitian reciprocal "
            "stationary main at xi=T_0; all endpoint and remainder terms "
            "remain explicit obligations."
        ),
        "pi_provenance": (
            "The sine is sin(2*pi*a^2*mu) because Poisson summation uses "
            "e(t)=exp(2*pi*i*t) and T_0=2*pi*a^2. The flow multiplier is "
            "i*mu because partial_xi e[xi*mu/(2*pi)]=i*mu e[xi*mu/(2*pi)]. "
            "Thus the factor 1/(2*pi) in the alpha-derivative identity is "
            "forced by the same two inherited normalizations; no new pi is "
            "inserted."
        ),
        "next_target": (
            "Derive the endpoint-complete discrete Poisson/B-process formula "
            "with uniform value and first-x derivative remainders. In its "
            "Hermitian main, bound the real leading offset sums using the "
            "reciprocal-lattice split, prove a source-specific bound for the "
            "imaginary coefficient I, and retain the endpoint, transpose, "
            "hard-cutoff, and adjacent-recurrence terms before taking "
            "absolute values."
        ),
        "exact_target": (
            "Prove Psi(T_0)-epsilon*integral_0^1 "
            "Psi'(T_0-theta*epsilon)dtheta<=h^2/400, preferably h^2/800."
        ),
        "proof_boundary": (
            "This proves the exact Hermitian reciprocal swap pairing, its "
            "outer-half sine/cosine formula, the real-leading reciprocal "
            "coordinates and stationary factor, the offset phase derivatives, "
            "the logarithmic and reciprocal-lattice defect bounds, complete "
            "pairing of the square interior dual main, and three algebraic "
            "nonpromotion guards. It proves no discrete Poisson or B-process "
            "formula, stationary-phase remainder, hard-cutoff boundary "
            "estimate, imaginary-correction bound, signed offset-sum gain, "
            "signed flow estimate, Phi_B bound, contact exclusion, retained "
            "aggregate or Xi theorem, Q209, cofinal descendant theorem, "
            "Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }


def build_rows(
    exact: dict,
    pair: dict,
    leading: dict,
    offset: dict,
    cutoff: dict,
) -> list[PairingRow]:
    return [
        PairingRow("hrpr_01_domain", "fixed-chart domain", "proved", "The reduction is restricted to the interior Hermitian stationary main.", exact["domain"], "Endpoint and remainder families are not absorbed."),
        PairingRow("hrpr_02_covariance", "swap covariance", "proved", "Hermitian dual coefficients and phases are conjugate under mode swap.", pair["swap_covariance"], "This assumes the interior stationary coefficient inherited from the Hermitian kernel."),
        PairingRow("hrpr_03_pair", "unordered pair formula", "proved", "Two ordered dual modes collapse to one real sine/cosine expression.", pair["ordered_pair"] + " " + pair["outer_half_pair"], "No absolute value has been taken."),
        PairingRow("hrpr_04_derivative", "phase derivative", "proved", "The paired flow is the exact frequency derivative of the paired stationary phase.", pair["phase_derivative"], exact["pi_provenance"]),
        PairingRow("hrpr_05_diagonal", "diagonal null", "proved", "The entire complex Hermitian reciprocal diagonal vanishes.", pair["diagonal_null"], "This says nothing by itself about adjacent modes."),
        PairingRow("hrpr_06_split", "real/imaginary split", "proved", "The leading real and corrective imaginary channels have different phase zeros.", pair["real_defect_split"], "The imaginary correction must be bounded independently."),
        PairingRow("hrpr_07_coordinates", "reciprocal coordinates", "proved", "The saddle difference and sum coordinates are mu and tau.", leading["reciprocal_coordinates"], "Continuous saddle algebra only."),
        PairingRow("hrpr_08_stationary", "mixed-signature factor", "proved", "The Hermitian stationary factor is positive and carries no Maslov phase.", leading["stationary_factor"], "No uniform stationary-phase remainder is supplied."),
        PairingRow("hrpr_09_leading", "leading sine channel", "proved", "The ideal q=1 Hermitian main is one explicit unsigned sine channel.", leading["leading_coefficient"] + " " + leading["leading_unordered_pair"], "The sine oscillation has no pointwise sign."),
        PairingRow("hrpr_10_offset", "offset reduction", "proved", "Swap pairing rewrites the off-diagonal square as triangular logarithmic offset sums.", offset["offset_coordinates"], "This is an exact phase organization, not a sum estimate."),
        PairingRow("hrpr_11_log", "logarithmic defect", "proved", "The nonlinear log phase has a uniform quadratic Taylor defect.", offset["log_bounds"], offset["proof_of_lattice_defect"]),
        PairingRow("hrpr_12_lattice", "reciprocal-lattice defect", "proved", "Near an integral reciprocal saddle, the leading sine gains an explicit second defect factor.", offset["lattice_defect"] + " " + offset["sine_bound"], "This bound can be weak when rho_k is not small."),
        PairingRow("hrpr_13_mixed", "full mixed bound", "proved", "The exact bound isolates the surviving imaginary-correction obligation.", offset["full_mixed_bound"], "No source-specific bound for I is proved here."),
        PairingRow("hrpr_14_cutoff", "dual-square pairing", "proved", "Every off-diagonal interior main mode has exactly one swap partner.", cutoff["dual_square"], cutoff["boundary_guard"]),
        PairingRow("hrpr_15_guards", "sign and resonance guards", "proved", "Swap symmetry and integral phase cannot be promoted to a sign theorem.", pair["sign_and_resonance_guards"], "The witnesses are algebraic, not asserted Xi configurations."),
        PairingRow("hrpr_16_handoff", "analytic handoff", "open", "The paired main is ready for a resonant/nonresonant estimate inside the endpoint-complete formula.", offset["resonant_nonresonant_handoff"] + " " + exact["next_target"], exact["exact_target"]),
        PairingRow("hrpr_17_boundary", "proof boundary", "proved", "The exact reduction is separated from the missing analytic theorem.", exact["proof_boundary"], "This is not a proof of RH."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    pair = payload["swap_pair_certificate"]
    leading = payload["reciprocal_leading_certificate"]
    offset = payload["offset_and_defect_certificate"]
    cutoff = payload["cutoff_certificate"]
    return f"""# Hermitian Reciprocal Pairing Reduction

Date: 2026-08-01

Status: exact reduction artifact for the interior Hermitian reciprocal main;
`0 Poisson remainder bounds`, `0 signed flow bounds`, and this is not a proof
of RH.

## Domain

{exact['domain']}

## Swap Pairing

{pair['swap_covariance']}

{pair['ordered_pair']}

{pair['outer_half_pair']}

{pair['phase_derivative']}

{pair['diagonal_null']}

{pair['real_defect_split']}

## Leading Reciprocal Main

{leading['reciprocal_coordinates']}

{leading['stationary_factor']}

{leading['leading_coefficient']}

{leading['leading_unordered_pair']}

## Offset And Lattice Defect

{offset['offset_coordinates']}

{offset['log_bounds']}

{offset['lattice_defect']}

{offset['proof_of_lattice_defect']}

{offset['sine_bound']}

{offset['full_mixed_bound']}

{offset['resonant_nonresonant_handoff']}

## Cutoff And Guards

{cutoff['dual_square']}

{cutoff['boundary_guard']}

{pair['sign_and_resonance_guards']}

## Pi Provenance

{exact['pi_provenance']}

## Handoff

{exact['next_target']}

{exact['exact_target']}

## Proof Boundary

{exact['proof_boundary']}

## Reproduce

```text
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```
"""


def build_payload() -> dict:
    sources = load_sources()
    pair = swap_pair_certificate()
    leading = reciprocal_leading_certificate()
    offset = offset_and_defect_certificate()
    cutoff = cutoff_certificate()
    exact = exact_payload()
    rows = build_rows(exact, pair, leading, offset, cutoff)
    return {
        "kind": STEM,
        "date": "2026-08-01",
        "status": (
            "exact Hermitian reciprocal swap, offset, and lattice-defect "
            "reduction complete; endpoint-complete signed estimate open"
        ),
        "source_audit": source_audit(sources),
        "swap_pair_certificate": pair,
        "reciprocal_leading_certificate": leading,
        "offset_and_defect_certificate": offset,
        "cutoff_certificate": cutoff,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "hermitian_swap_pair_identities": 1,
            "exact_full_diagonal_nulls": 1,
            "leading_real_sine_channels": 1,
            "reciprocal_lattice_defect_bounds": 1,
            "unpaired_interior_dual_main_modes": 0,
            "sign_guards": 2,
            "imaginary_resonance_guards": 1,
            "poisson_remainder_bounds": 0,
            "hard_cutoff_boundary_bounds": 0,
            "imaginary_correction_bounds": 0,
            "signed_offset_sum_bounds": 0,
            "signed_flow_bounds": 0,
            "phi_b_bounds": 0,
        },
        "proof_boundary": exact["proof_boundary"],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    atomic_write(args.out, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    atomic_write(args.note, render_note(payload))
    print(
        "built Hermitian reciprocal pairing reduction: 17 rows, "
        "1 swap identity, 1 lattice-defect bound, "
        "0 Poisson remainder bounds, 0 signed flow bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
