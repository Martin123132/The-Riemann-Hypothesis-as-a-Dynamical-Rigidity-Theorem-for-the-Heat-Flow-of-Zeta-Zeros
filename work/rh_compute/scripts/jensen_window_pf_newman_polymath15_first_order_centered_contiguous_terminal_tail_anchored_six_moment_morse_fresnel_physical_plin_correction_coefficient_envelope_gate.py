#!/usr/bin/env python3
"""Build the physical correction-coefficient envelope gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
    "physical_plin_correction_coefficient_envelope_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "critical_ray": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "critical_ray_finite_height_real_edge_gate.json"
    ),
    "absolute_anchor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "absolute_phase_anchor_reduction.json"
    ),
    "real_residual": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "real_residual_reduction.json"
    ),
    "two_carrier": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_moment_two_carrier_"
        "kernel_reduction.json"
    ),
    "physical_plin": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
        "physical_plin_amplitude_gate.json"
    ),
}

H_DENOMINATOR = 72_000_000_000
SMALLNESS_DENOMINATOR = 16_892


@dataclass(frozen=True)
class EnvelopeRow:
    id: str
    role: str
    readiness: str
    claim: str
    certificate: str
    proof_boundary: str


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCE_PATHS.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCE_PATHS.items()
    }


def source_audit(payloads: dict[str, dict]) -> dict:
    critical = payloads["critical_ray"]
    critical_bounds = critical.get("majorant_certificate", {}).get(
        "source_bounds", {}
    )
    required_critical = {
        "alpha_prime": "|alpha'|<h^2",
        "alpha_second": "|alpha''|<h^4",
        "correction": "|d|<h^2/4, |d_x|<h^4/32",
    }
    for key, value in required_critical.items():
        if critical_bounds.get(key) != value:
            raise RuntimeError(f"critical-ray source drifted at {key}")
    domain = critical.get("exact", {}).get("domain", {})
    if "h<exp(-25)<1/72000000000" not in domain.get(
        "effective_box", ""
    ):
        raise RuntimeError("critical-ray h cap drifted")
    chain = critical.get("majorant_certificate", {}).get(
        "normalized_majorant_chain", {}
    )
    for key, cap in (
        ("chi_absolute", "2/1"),
        ("alpha_prime_over_h2", "1/1"),
        ("alpha_second_over_h4", "1/1"),
        ("d_over_h2", "1/4"),
        ("d_x_over_h4", "1/32"),
    ):
        if chain.get(key, {}).get("cap") != cap:
            raise RuntimeError(f"critical-ray normalized cap drifted at {key}")

    anchor = payloads["absolute_anchor"].get("exact", {}).get(
        "first_coefficient", ""
    )
    if "|d_1|<2189/x<1/2" not in anchor:
        raise RuntimeError("absolute-anchor denominator bound drifted")

    derivative = (
        payloads["real_residual"]
        .get("exact", {})
        .get("main_and_derivative_bounds", {})
        .get("correction_derivative", "")
    )
    if "|d_(n,x)|<4223/x^2" not in derivative:
        raise RuntimeError("all-carrier derivative bound drifted")

    physical = (
        payloads["two_carrier"]
        .get("finite_sum_certificate", {})
        .get("physical_chi_bound", "")
    )
    for token in ("x>a^2=h^-2", "16892h^2<1"):
        if token not in physical:
            raise RuntimeError(f"two-carrier scale relation drifted: {token}")

    coefficient_counts = payloads["physical_plin"].get("counts", {})
    if coefficient_counts.get("exact_centered_complex_coefficients") != 6:
        raise RuntimeError("physical P_lin source lost six coefficients")
    if coefficient_counts.get("maximum_polynomial_degree") != 5:
        raise RuntimeError("physical P_lin source lost degree-five closure")

    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def require_zero(expression: sp.Expr, label: str) -> None:
    if sp.cancel(sp.expand(expression)) != 0:
        raise RuntimeError(label)


def symbolic_certificate() -> dict:
    A, alpha, alpha_x = sp.symbols("A alpha alpha_x")
    B, B_x, a_0, a_0_x = sp.symbols("B B_x a_0 a_0_x", nonzero=True)
    d_A, d_A_x = sp.symbols("d_A d_A_x")
    chi = alpha - A

    rho_2 = B / a_0
    rho_1 = -2 * alpha * rho_2
    rho_2_x = B_x / a_0 - B * a_0_x / a_0**2
    rho_1_x = -2 * alpha_x * rho_2 - 2 * alpha * rho_2_x

    c_0 = 1 + rho_1 * A + rho_2 * A**2
    c_1 = rho_1 + 2 * A * rho_2
    c_2 = rho_2
    e_0 = rho_1_x * A + rho_2_x * A**2
    e_1 = rho_1_x + 2 * A * rho_2_x
    e_2 = rho_2_x

    one_plus_d_A = a_0 - 2 * B * alpha * A + B * A**2
    d_A_x_identity = (
        a_0_x
        - 2 * B_x * alpha * A
        - 2 * B * alpha_x * A
        + B_x * A**2
    )
    ratio_c_0 = (1 + d_A) / a_0
    ratio_e_0 = d_A_x / a_0 - (1 + d_A) * a_0_x / a_0**2

    require_zero(c_1 + 2 * chi * rho_2, "centered c_1 identity failed")
    require_zero(
        e_1 + 2 * alpha_x * rho_2 + 2 * chi * rho_2_x,
        "centered d_1 identity failed",
    )
    require_zero(
        (c_0 - ratio_c_0).subs(d_A, one_plus_d_A - 1),
        "centered c_0 ratio failed",
    )
    require_zero(
        (e_0 - ratio_e_0)
        .subs(d_A, one_plus_d_A - 1)
        .subs(d_A_x, d_A_x_identity),
        "centered d_0 ratio failed",
    )

    return {
        "source_correction": (
            "Write the original first Dirichlet correction as frak_d(ell)="
            "1/(6s)+alpha'(s){t/4+t^2[alpha(s)-ell]^2/8}. Put "
            "A=log(a), chi=alpha-A, a_0=1+frak_d(0), and "
            "B=alpha'(s)t^2/8."
        ),
        "normalized_coefficients": (
            "rho_2=B/a_0 and rho_1=-2alpha*rho_2. At fixed t and with "
            "s_x=-i/2, B_x=-i*alpha''(s)t^2/16 and "
            "rho_(2,x)=B_x/a_0-B*frak_d_(1,x)/a_0^2."
        ),
        "centered_value_coefficients": (
            "For C(A+y)=c_0+c_1*y+c_2*y^2, c_0=[1+frak_d(A)]/a_0, "
            "c_1=-2chi*rho_2, and c_2=rho_2."
        ),
        "centered_derivative_coefficients": (
            "For D(A+y)=d_0+d_1*y+d_2*y^2 with D=partial_x C at fixed "
            "ell, d_0=frak_d_x(A)/a_0-[1+frak_d(A)]frak_d_(1,x)/a_0^2, "
            "d_1=-2alpha_x*rho_2-2chi*rho_(2,x), and d_2=rho_(2,x), "
            "where alpha_x=-i*alpha'/2."
        ),
        "symbolic_differences": {
            "c_0_ratio": "0",
            "c_1_centering": "0",
            "d_0_ratio": "0",
            "d_1_centering": "0",
        },
    }


def bound_certificate() -> dict:
    h_0 = Fraction(1, H_DENOMINATOR)
    smallness = Fraction(SMALLNESS_DENOMINATOR) * h_0**2
    if not smallness < 1:
        raise RuntimeError("16892 h^2 smallness failed")

    b_bound = Fraction(1, 32)
    b_x_bound = Fraction(1, 64)
    rho_2_bound = 2 * b_bound
    rho_2_x_bound = (
        2 * b_x_bound
        + 4
        * b_bound
        * Fraction(4_223, SMALLNESS_DENOMINATOR)
    )
    if rho_2_bound != Fraction(1, 16):
        raise RuntimeError("rho_2 coefficient drifted")
    if rho_2_x_bound != Fraction(1, 16):
        raise RuntimeError("rho_(2,x) coefficient drifted")

    c_0_raw = 2 * (Fraction(1, 4) + 2_189)
    c_0_round = Fraction(4_379)
    if not c_0_raw < c_0_round:
        raise RuntimeError("c_0 rounding failed")
    if not 3 * c_0_round < SMALLNESS_DENOMINATOR:
        raise RuntimeError("c_0 nonvanishing guard failed")

    c_1_bound = 2 * 2 * rho_2_bound
    c_2_bound = rho_2_bound
    if c_1_bound != Fraction(1, 4):
        raise RuntimeError("c_1 coefficient drifted")

    d_0_raw = (
        2 * Fraction(1, 32)
        + 4
        * (1 + Fraction(1, 4 * SMALLNESS_DENOMINATOR))
        * 4_223
    )
    d_0_round = Fraction(16_893)
    if not d_0_raw < d_0_round:
        raise RuntimeError("d_0 rounding failed")
    d_1_bound = (
        2 * Fraction(1, 2) * rho_2_bound
        + 2 * 2 * rho_2_x_bound
    )
    d_2_bound = rho_2_x_bound
    if d_1_bound != Fraction(5, 16):
        raise RuntimeError("d_1 coefficient drifted")

    return {
        "domain": (
            "L>=50, 0<tL<=25, 0<t<=1/2, h=1/a<1/72000000000, "
            "x>a^2=h^-2, and 16892h^2<1."
        ),
        "separated_source_bounds": {
            "critical_geometry": (
                "|alpha'|<h^2, |alpha''|<h^4, |chi|<2. The same "
                "normalized correction majorant applies at ell=A because "
                "alpha-A=chi, giving |frak_d(A)|<h^2/4 and "
                "|frak_d_x(A)|<h^4/32."
            ),
            "all_carrier_anchor": (
                "|frak_d(0)|<2189/x<2189h^2<1/2 and "
                "|frak_d_x(0)|<4223/x^2<4223h^4."
            ),
        },
        "denominator": (
            "a_0=1+frak_d(0), so |a_0|>1/2, |a_0|^-1<2, and "
            "|a_0|^-2<4."
        ),
        "quotient_bounds": {
            "B": "|B|<h^2/32",
            "B_x": "|B_x|<h^4/64",
            "rho_2": "|rho_2|<h^2/16",
            "rho_2_x": (
                "|rho_(2,x)|<h^4/32+(4223/8)h^6<h^4/16"
            ),
        },
        "centered_value_bounds": {
            "c_0": "|c_0-1|<4379h^2<1/3; hence 2/3<|c_0|<4/3",
            "c_1": "|c_1|<h^2/4",
            "c_2": "|c_2|<h^2/16",
        },
        "centered_derivative_bounds": {
            "d_0": "|d_0|<16893h^4",
            "d_1": "|d_1|<5h^4/16",
            "d_2": "|d_2|<h^4/16",
        },
        "rational_diagnostics": {
            "h_upper": fraction_text(h_0),
            "16892_h_upper_squared": fraction_text(smallness),
            "rho_2_over_h2": fraction_text(rho_2_bound),
            "rho_2_x_over_h4_using_16892h2_lt_1": fraction_text(
                rho_2_x_bound
            ),
            "c_0_minus_1_raw_over_h2": fraction_text(c_0_raw),
            "c_0_minus_1_rounded_over_h2": fraction_text(c_0_round),
            "d_0_raw_over_h4_using_16892h2_lt_1": fraction_text(d_0_raw),
            "d_0_rounded_over_h4": fraction_text(d_0_round),
            "d_1_over_h4": fraction_text(d_1_bound),
        },
    }


def build_rows(symbolic: dict, bounds: dict) -> list[EnvelopeRow]:
    return [
        EnvelopeRow(
            "pccenv_01_domain",
            "effective_domain",
            "ready_to_apply",
            "The coefficient envelope uses the complete critical chart.",
            bounds["domain"],
            "No q<1 parabolic theorem or Xi-level sign is inferred.",
        ),
        EnvelopeRow(
            "pccenv_02_source_separation",
            "source_audit",
            "certified",
            "Terminal-centered and n=1 denominator estimates retain distinct roles.",
            (
                bounds["separated_source_bounds"]["critical_geometry"]
                + " "
                + bounds["separated_source_bounds"]["all_carrier_anchor"]
            ),
            "A terminal-only d bound is not applied to n=1.",
        ),
        EnvelopeRow(
            "pccenv_03_exact_correction",
            "exact_identity",
            "certified",
            "The first correction is an exact quadratic in log n.",
            symbolic["source_correction"],
            "No asymptotic correction remainder enters this coefficient identity.",
        ),
        EnvelopeRow(
            "pccenv_04_exact_quotient",
            "exact_identity",
            "certified",
            "rho_2 and its physical x derivative are exact quotients.",
            symbolic["normalized_coefficients"],
            "The derivative is taken at fixed t.",
        ),
        EnvelopeRow(
            "pccenv_05_denominator",
            "nonvanishing_guard",
            "certified",
            "The first-carrier normalizer never vanishes.",
            bounds["denominator"],
            "This is a modulus bound and uses no choice of phase branch.",
        ),
        EnvelopeRow(
            "pccenv_06_B_bounds",
            "analytic_majorant",
            "certified",
            "The quadratic correction coefficient and its derivative have h-scale bounds.",
            "|B|<h^2/32 and |B_x|<h^4/64.",
            "Strictness comes from the certified alpha bounds.",
        ),
        EnvelopeRow(
            "pccenv_07_rho_bounds",
            "quotient_majorant",
            "certified",
            "Both normalized quadratic coefficients meet the candidate targets.",
            (
                bounds["quotient_bounds"]["rho_2"]
                + "; "
                + bounds["quotient_bounds"]["rho_2_x"]
            ),
            "No P_lin or mode-sum norm follows yet.",
        ),
        EnvelopeRow(
            "pccenv_08_value_centering",
            "exact_identity",
            "certified",
            "The value channel centers on chi rather than separate large logarithms.",
            symbolic["centered_value_coefficients"],
            "A=log(a) is a parameter in the polynomial shift.",
        ),
        EnvelopeRow(
            "pccenv_09_derivative_centering",
            "exact_identity",
            "certified",
            "The derivative channel also centers without a free log(a) loss.",
            symbolic["centered_derivative_coefficients"],
            "D is the x derivative at fixed original ell before shifting.",
        ),
        EnvelopeRow(
            "pccenv_10_c0_ratio",
            "normalized_ratio",
            "certified",
            "The constant value channel is controlled as one normalized correction ratio.",
            bounds["centered_value_bounds"]["c_0"],
            "No expansion of log(a)^2 is bounded termwise.",
        ),
        EnvelopeRow(
            "pccenv_11_value_envelopes",
            "coefficient_envelope",
            "certified",
            "All three centered C coefficients are uniformly bounded.",
            "; ".join(bounds["centered_value_bounds"].values()),
            "These are coefficient bounds, not a bound for C on an unbounded y-line.",
        ),
        EnvelopeRow(
            "pccenv_12_d0_ratio",
            "normalized_ratio",
            "certified",
            "The constant derivative channel is controlled before large logs appear.",
            bounds["centered_derivative_bounds"]["d_0"],
            "The all-carrier d_(1,x) estimate is retained in the denominator term.",
        ),
        EnvelopeRow(
            "pccenv_13_derivative_envelopes",
            "coefficient_envelope",
            "certified",
            "All three centered D coefficients are uniformly bounded.",
            "; ".join(bounds["centered_derivative_bounds"].values()),
            "The large constant 16893 remains visible for subsequent propagation.",
        ),
        EnvelopeRow(
            "pccenv_14_constant_channel",
            "nonvanishing_guard",
            "certified",
            "The centered value channel has a nonvanishing constant coefficient.",
            "2/3<|c_0|<4/3.",
            "This does not make C zero-free away from y=0.",
        ),
        EnvelopeRow(
            "pccenv_15_handoff",
            "route_decision",
            "open_quantitative_handoff",
            "The six correction-channel envelopes are ready for P_lin propagation.",
            (
                "Insert them into the exact gamma/eta convolution for p_0,...,p_5, "
                "retaining p_5=i*rho_2*(X_T+V_p)*(s_*')^2."
            ),
            "Physical observation and saddle-operator bounds remain open.",
        ),
        EnvelopeRow(
            "pccenv_16_boundary",
            "proof_boundary",
            "guard_validated",
            "The coefficient theorem is not promoted to a signed flow theorem.",
            (
                "There are zero P_lin norm bounds, grouped Morse bounds, quadratic "
                "pairing bounds, Phi_B bounds, or Xi-level signs in this gate."
            ),
            "No contact exclusion, Lambda<=0, PF-infinity, RH, or prize theorem is claimed.",
        ),
    ]


def render_note(artifact: dict) -> str:
    symbolic = artifact["symbolic_certificate"]
    bounds = artifact["bound_certificate"]
    return "\n".join(
        [
            "# Physical P_lin Correction-Coefficient Envelope Gate",
            "",
            "Date: 2026-08-02",
            "",
            "Status: exact quotient identities and certified physical envelopes for the six terminal-centered C/D coefficients. This is not a proof of a physical P_lin norm, signed Morse sum, or RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Source Separation",
            "",
            bounds["domain"],
            "",
            "To avoid the collision with the centered coefficient d_1, write the original first Dirichlet correction as frak_d(ell). The critical-ray gate supplies the alpha, chi, and ell=A correction bounds; the all-carrier gates separately supply the ell=0 denominator and derivative bounds.",
            "",
            "```text",
            bounds["separated_source_bounds"]["critical_geometry"],
            bounds["separated_source_bounds"]["all_carrier_anchor"],
            "```",
            "",
            "## Exact Quotients",
            "",
            symbolic["source_correction"],
            "",
            "```text",
            "rho_2=B/a_0,                 rho_1=-2*alpha*rho_2,",
            "B_x=-i*alpha''*t^2/16,",
            "rho_(2,x)=B_x/a_0-B*frak_d_(1,x)/a_0^2.",
            "```",
            "",
            bounds["denominator"],
            "",
            "Therefore",
            "",
            "```text",
            bounds["quotient_bounds"]["B"],
            bounds["quotient_bounds"]["B_x"],
            bounds["quotient_bounds"]["rho_2"],
            bounds["quotient_bounds"]["rho_2_x"],
            "```",
            "",
            "The last step uses 16892=4*4223 and 16892h^2<1. Thus both candidate quotient bounds are now certified rather than merely suggested.",
            "",
            "## Centered Coefficients",
            "",
            "Put A=log(a), chi=alpha-A, and shift ell=A+y. The exact rows are",
            "",
            "```text",
            "c_0=[1+frak_d(A)]/a_0,",
            "c_1=-2*chi*rho_2,             c_2=rho_2,",
            "",
            "d_0=frak_d_x(A)/a_0-[1+frak_d(A)]frak_d_(1,x)/a_0^2,",
            "d_1=-2*alpha_x*rho_2-2*chi*rho_(2,x),",
            "d_2=rho_(2,x),                alpha_x=-i*alpha'/2.",
            "```",
            "",
            "Consequently, without separately majorizing log(a) or log(a)^2,",
            "",
            "```text",
            bounds["centered_value_bounds"]["c_0"],
            bounds["centered_value_bounds"]["c_1"],
            bounds["centered_value_bounds"]["c_2"],
            "",
            bounds["centered_derivative_bounds"]["d_0"],
            bounds["centered_derivative_bounds"]["d_1"],
            bounds["centered_derivative_bounds"]["d_2"],
            "```",
            "",
            "The 16893 constant is intentionally retained. It is small after multiplication by h^4, but it must remain visible when the p_j convolution is bounded.",
            "",
            "## Next Target",
            "",
            "Insert these six envelopes into the exact gamma/eta recurrence of the physical P_lin gate. Bound p_0,...,p_5 only after the physical observation coefficients, R, R_x, and delta are joined; retain the exact top factor p_5=i*rho_2*(X_T+V_p)*(s_*')^2 and the total-value fibre.",
            "",
            "## Proof Boundary",
            "",
            "This gate proves no P_lin norm, grouped c-prime sum, reciprocal cancellation estimate, quadratic residual bound, h^2 flow error, Phi_B upper bound, contact exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
            "",
            "## Reproduce",
            "",
            "```powershell",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
        ]
    )


def build_artifact() -> dict:
    payloads = load_sources()
    audit = source_audit(payloads)
    symbolic = symbolic_certificate()
    bounds = bound_certificate()
    rows = build_rows(symbolic, bounds)
    return {
        "kind": STEM,
        "date": "2026-08-02",
        "status": (
            "exact rho_2 quotient and derivative bounds plus six certified "
            "terminal-centered correction-coefficient envelopes"
        ),
        "source_audit": audit,
        "symbolic_certificate": symbolic,
        "bound_certificate": bounds,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "sources": len(audit),
            "quotient_bounds": 2,
            "centered_coefficients": 6,
            "coefficient_envelopes": 6,
            "nonvanishing_constant_channels": 1,
            "physical_plin_bounds": 0,
            "grouped_interior_bounds": 0,
            "signed_flow_bounds": 0,
        },
        "proof_boundary": (
            "Exact correction quotient and six coefficient envelopes only; no "
            "P_lin norm, grouped Morse bound, quadratic residual estimate, signed "
            "flow theorem, Phi_B bound, contact exclusion, Xi theorem, Lambda<=0, "
            "PF-infinity, RH, or prize-level conclusion."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()

    artifact = build_artifact()
    atomic_write(args.output, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(args.note, render_note(artifact))
    counts = artifact["counts"]
    print(
        "built physical correction-coefficient envelope gate: "
        f"{counts['rows']} rows, {counts['quotient_bounds']} quotient bounds, "
        f"{counts['centered_coefficients']} centered coefficients, "
        f"{counts['coefficient_envelopes']} coefficient envelopes, "
        f"{counts['nonvanishing_constant_channels']} nonvanishing constant channel, "
        f"{counts['physical_plin_bounds']} P_lin bounds, "
        f"{counts['signed_flow_bounds']} signed flow bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
