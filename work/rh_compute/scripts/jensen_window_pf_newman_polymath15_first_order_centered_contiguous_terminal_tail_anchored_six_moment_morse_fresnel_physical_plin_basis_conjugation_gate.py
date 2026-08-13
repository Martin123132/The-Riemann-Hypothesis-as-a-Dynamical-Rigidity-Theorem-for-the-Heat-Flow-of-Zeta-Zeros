#!/usr/bin/env python3
"""Build the centered-basis Morse conjugation and tail-separation gate."""

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
    "physical_plin_basis_conjugation_gate"
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
    "q1_saddle": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_moment_q1_saddle_"
        "phase_variance_reduction.json"
    ),
    "observation_image": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
        "observation_image_compression_gate.json"
    ),
    "physical_plin": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
        "physical_plin_amplitude_gate.json"
    ),
    "joined_envelope": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
        "physical_plin_joined_coefficient_envelope_gate.json"
    ),
}

K_VALUES = (6, 14, 55, 336, 2738, 27936)
H_DENOMINATOR = 72_000_000_000
PI_CAPS = (8764, 8765, 8762, 8760)


@dataclass(frozen=True)
class GateRow:
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


def require_zero(expression: sp.Expr, label: str) -> None:
    if sp.simplify(sp.expand(expression)) != 0:
        raise RuntimeError(label)


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


def source_audit(payloads: dict[str, dict]) -> dict[str, dict[str, str]]:
    domain = payloads["critical_ray"].get("exact", {}).get("domain", {})
    effective = domain.get("effective_box", "")
    if "h<exp(-25)<1/72000000000" not in effective:
        raise RuntimeError("critical-ray exponential height cap drifted")

    q1 = payloads["q1_saddle"].get("physical_q1_certificate", {})
    for key, marker in (
        ("domain", "a^2=exp(L)+1/(32L^2)"),
        ("amplitude_defect", "delta_a=log(a)-Re(alpha)"),
        ("b_a", "B_a=1/2-delta_a/(4L^2)"),
        ("profile", "exp[t*log(a)^2/4-sigma*log(a)]"),
    ):
        if marker not in q1.get(key, ""):
            raise RuntimeError(f"q=1 saddle source drifted: {key}")

    observation = payloads["observation_image"]
    tail = observation.get("tail_certificate", {})
    if tuple(tail.get("constants", ())) != K_VALUES:
        raise RuntimeError("outside-tail constants drifted")
    if "P(lambda)=sum_(j=0)^5p_jlambda^j" not in tail.get(
        "projection_template", ""
    ):
        raise RuntimeError("outside-tail basis drifted")

    plin = payloads["physical_plin"]
    morse = plin.get("morse_amplitude_certificate", {})
    for key, marker in (
        ("cprime_formula", "zJ^2(P'+gP)/v"),
        ("saddle_formula", "P''+(2g_0+1)P'"),
    ):
        if marker not in morse.get(key, ""):
            raise RuntimeError(f"Morse P_lin source drifted: {key}")
    rows = plin.get("coefficient_certificate", {}).get(
        "observation_rows", ""
    )
    if "P_lin=F_alpha+i*x*F_beta" not in rows:
        raise RuntimeError("physical joined-row factorization drifted")

    joined = payloads["joined_envelope"]
    if joined.get("counts", {}).get("centered_plin_coefficients_propagated") != 6:
        raise RuntimeError("joined centered-coefficient count drifted")
    if joined.get("counts", {}).get("numerical_observation_row_bounds") != 0:
        raise RuntimeError("joined proof boundary drifted")
    p_bounds = joined.get("bound_certificate", {}).get(
        "p_deviation_envelopes", {}
    )
    if "|p_5|<h^2|beta_N|/63" not in p_bounds.get("Pi_5", ""):
        raise RuntimeError("joined p_5 bound drifted")

    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def symbolic_certificate() -> dict:
    lam, ell, x = sp.symbols("lambda ell x", real=True)
    p = sp.symbols("p_0:6")
    q_symbols = sp.symbols("q_0:6")

    centered = sum(p[k] * x**k for k in range(6))
    q = [
        sp.expand(
            sum(
                sp.binomial(k, j) * (-ell) ** (k - j) * p[k]
                for k in range(j, 6)
            )
        )
        for j in range(6)
    ]
    uncentered = sum(q[j] * lam**j for j in range(6))
    require_zero(
        uncentered.subs(lam, ell + x) - centered,
        "forward centered-to-lambda map failed",
    )

    p_from_q = [
        sp.expand(
            sum(
                sp.binomial(k, j) * ell ** (k - j) * q_symbols[k]
                for k in range(j, 6)
            )
        )
        for j in range(6)
    ]
    q_roundtrip = [
        sp.expand(
            sum(
                sp.binomial(k, j) * (-ell) ** (k - j) * p_from_q[k]
                for k in range(j, 6)
            )
        )
        for j in range(6)
    ]
    for j in range(6):
        require_zero(q_roundtrip[j] - q_symbols[j], f"inverse basis row {j}")

    derivative_q = [
        sp.expand(sp.diff(uncentered, lam, j).subs(lam, 0) / sp.factorial(j))
        for j in range(6)
    ]
    for j in range(6):
        require_zero(derivative_q[j] - q[j], f"lambda jet row {j}")

    hat_k = [
        sp.expand(
            sum(
                sp.binomial(k, j) * K_VALUES[j] * ell ** (k - j)
                for j in range(k + 1)
            )
        )
        for k in range(6)
    ]

    t, sigma, b_a = sp.symbols("t sigma B_a", real=True)
    S = t * lam**2 / 4 - sigma * lam
    S_tilde = t * x**2 / 4 - b_a * x
    b_sub = sigma - t * ell / 2
    alpha_real, delta_a = sp.symbols("alpha_real delta_a", real=True)
    physical_sigma = sp.Rational(1, 2) + t * alpha_real / 2
    physical_b = sp.Rational(1, 2) - t * delta_a / 2
    require_zero(
        (
            physical_b - (physical_sigma - t * ell / 2)
        ).subs(delta_a, ell - alpha_real),
        "physical B_a identity failed",
    )
    require_zero(
        (S.subs(lam, ell + x) - S.subs(lam, ell) - S_tilde).subs(
            b_a, b_sub
        ),
        "centered exponent conjugation failed",
    )
    g = t * lam / 2 - sigma
    g_tilde = t * x / 2 - b_a
    require_zero(
        (g.subs(lam, ell + x) - g_tilde).subs(b_a, b_sub),
        "centered logarithmic derivative failed",
    )

    P = centered
    first = sp.diff(P, x) + g_tilde * P
    saddle = (
        sp.diff(P, x, 2)
        + (2 * g_tilde + 1) * sp.diff(P, x)
        + (g_tilde**2 + g_tilde + t / 2 + sp.Rational(1, 6)) * P
    )
    conjugated = sp.exp(-S_tilde) * (
        sp.diff(sp.exp(S_tilde) * P, x, 2)
        + sp.diff(sp.exp(S_tilde) * P, x)
        + sp.exp(S_tilde) * P / 6
    )
    require_zero(conjugated - saddle, "saddle conjugation identity failed")

    alpha_row = sp.symbols("alpha_0:4", real=True)
    beta_row = sp.symbols("beta_0:4", real=True)
    fields = [sp.Function(f"F_{j}")(x) for j in range(4)]
    f_alpha = sum(alpha_row[j] * fields[j] for j in range(4))
    f_beta = sum(beta_row[j] * fields[j] for j in range(4))
    p_joined = f_alpha + sp.I * x * f_beta

    def d_g(expression: sp.Expr) -> sp.Expr:
        return sp.diff(expression, x) + g_tilde * expression

    def saddle_op(expression: sp.Expr) -> sp.Expr:
        return (
            sp.diff(expression, x, 2)
            + (2 * g_tilde + 1) * sp.diff(expression, x)
            + (g_tilde**2 + g_tilde + t / 2 + sp.Rational(1, 6))
            * expression
        )

    expected_first = sum(
        alpha_row[j] * d_g(fields[j])
        + sp.I
        * beta_row[j]
        * (x * d_g(fields[j]) + fields[j])
        for j in range(4)
    )
    require_zero(d_g(p_joined) - expected_first, "joined first order failed")

    expected_saddle = sum(
        alpha_row[j] * saddle_op(fields[j])
        + sp.I
        * beta_row[j]
        * (
            x * saddle_op(fields[j])
            + 2 * sp.diff(fields[j], x)
            + (2 * g_tilde + 1) * fields[j]
        )
        for j in range(4)
    )
    require_zero(
        saddle_op(p_joined) - expected_saddle,
        "joined saddle operator failed",
    )

    a_coef, b_coef, e_x, e_0, x_0 = sp.symbols(
        "A_z B_z E_x E_0 x_0"
    )
    p_at_0 = p_joined.subs(x, x_0)
    off_saddle = e_x * (
        a_coef * d_g(p_joined) + b_coef * p_joined
    ) + e_0 * p_at_0
    expected_off = 0
    for j, field in enumerate(fields):
        field_0 = field.subs(x, x_0)
        alpha_kernel = e_x * (
            a_coef * d_g(field) + b_coef * field
        ) + e_0 * field_0
        beta_kernel = e_x * (
            a_coef * (x * d_g(field) + field)
            + b_coef * x * field
        ) + e_0 * x_0 * field_0
        expected_off += (
            alpha_row[j] * alpha_kernel
            + sp.I * beta_row[j] * beta_kernel
        )
    require_zero(off_saddle - expected_off, "joined off-saddle grouping failed")

    return {
        "basis_translation": {
            "forward": (
                "q_j=sum_(k=j)^5 binom(k,j)(-ell)^(k-j)p_k, "
                "where ell=log(a)."
            ),
            "inverse": (
                "p_j=sum_(k=j)^5 binom(k,j)ell^(k-j)q_k."
            ),
            "jet_interpretation": (
                "q_j=P^(j)(0)/j! and p_j=P^(j)(ell)/j! for "
                "P(lambda)=sum p_k(lambda-ell)^k."
            ),
            "forward_rows": [f"q_{j}={sp.sstr(q[j])}" for j in range(6)],
            "inverse_rows": [
                f"p_{j}={sp.sstr(p_from_q[j])}" for j in range(6)
            ],
            "symbolic_zero_checks": 13,
        },
        "translated_tail_weights": {
            "definition": (
                "Khat_k(ell)=sum_(j=0)^k binom(k,j)K_j ell^(k-j), "
                "K=(6,14,55,336,2738,27936)."
            ),
            "rows": [
                f"Khat_{k}(ell)={sp.sstr(hat_k[k])}" for k in range(6)
            ],
            "sharpness": (
                "For a centered monomial P=p_k(lambda-ell)^k, the old "
                "coefficientwise majorant is exactly "
                "2h^2 Khat_k(ell)|p_k|. Thus no smaller universal "
                "coefficientwise weight follows from that template."
            ),
        },
        "physical_conjugation": {
            "parameter": (
                "B_a=sigma-(t/2)ell=1/2-(t/2)delta_a, "
                "ell=log(a), delta_a=ell-Re(alpha)."
            ),
            "exponent": (
                "S(ell+x)=S(ell)+S_tilde(x), "
                "S_tilde(x)=t*x^2/4-B_a*x."
            ),
            "log_derivative": (
                "g(ell+x)=t*x/2-B_a=:g_tilde(x)."
            ),
            "first_order": (
                "P'(lambda)+g(lambda)P(lambda)="
                "Ptilde'(x)+g_tilde(x)Ptilde(x)."
            ),
            "saddle_operator": (
                "L_B[P]=P''+(2g_tilde+1)P'+"
                "(g_tilde^2+g_tilde+t/2+1/6)P."
            ),
            "conjugate_form": (
                "L_B[P]=exp(-S_tilde)*(D_x^2+D_x+1/6)"
                "[exp(S_tilde)P]."
            ),
            "mode_coordinate": (
                "x_0=lambda_0-ell=log(alpha_P/(a*r)); log(a) enters "
                "the sampled interval through x_0, but not as an extra "
                "coefficient in L_B."
            ),
        },
        "joined_operator": {
            "field_vector": (
                "F=(C,mathcal N,A,Q), F_f=f dot F, and "
                "P_lin=alpha dot F+i*x*beta dot F."
            ),
            "first_order": (
                "(D_x+g_tilde)P_lin=alpha dot[(D_x+g_tilde)F]+"
                "i beta dot{x(D_x+g_tilde)F+F}."
            ),
            "saddle": (
                "L_B[P_lin]=alpha dot L_B[F]+i beta dot"
                "{x L_B[F]+2F'+(2g_tilde+1)F}."
            ),
            "off_saddle": (
                "With A_z=zJ^2/v, B_z=zJJ'-J, E_x=exp(S_tilde(x)), "
                "and E_0=exp(S_tilde(x_0)), factor exp(S(ell)) from "
                "the c-prime numerator. The remaining centered bracket is "
                "alpha dot K_alpha+i beta dot K_beta, "
                "where K_alpha=E_x[A_z(D+g_tilde)F+B_zF]+E_0F(x_0) "
                "and K_beta=E_x[A_z{x(D+g_tilde)F+F}+B_zxF]+"
                "E_0x_0F(x_0). Thus c-prime equals exp(S(ell)) times "
                "this bracket divided by r*sqrt(alpha_P)*z^2."
            ),
            "symbolic_zero_checks": 3,
        },
    }


def bound_certificate(symbolic: dict) -> dict:
    ell = sp.symbols("ell", nonnegative=True)
    hat_values = [
        sp.Poly(
            sum(
                sp.binomial(k, j) * K_VALUES[j] * ell ** (k - j)
                for j in range(k + 1)
            ),
            ell,
        )
        for k in range(6)
    ]
    at_50 = [int(poly.eval(50)) for poly in hat_values]
    expected_at_50 = [6, 314, 16455, 863586, 45394938, 2390362436]
    if at_50 != expected_at_50:
        raise RuntimeError("translated tail weights at L=50 drifted")

    denominator = H_DENOMINATOR**2
    h4_remainders = {
        "Pi_0": Fraction(2 * 16893, denominator),
        "Pi_1": Fraction(2 * 16893, denominator)
        + Fraction(10, 16 * denominator),
        "Pi_2": Fraction(2 * 16893, denominator)
        + Fraction(10, 16 * denominator)
        + Fraction(2, 16 * denominator),
        "Pi_3": Fraction(10, 16 * denominator)
        + Fraction(2, 16 * denominator),
    }
    if any(value >= 1 for value in h4_remainders.values()):
        raise RuntimeError("h^4 correction absorption failed")

    c_50 = (
        PI_CAPS[0] * at_50[0]
        + PI_CAPS[1] * at_50[1]
        + PI_CAPS[2] * at_50[2]
        + PI_CAPS[3] * at_50[3]
        + Fraction(2, 3) * at_50[4]
        + Fraction(1, 63) * at_50[5]
    )
    perturbation_ratio = 2 * c_50 / denominator
    target_ratio = Fraction(1, 300_000_000_000)
    if not perturbation_ratio < target_ratio:
        raise RuntimeError("translated perturbation tail rounding failed")

    return {
        "domain": (
            "On the physical q=1 chart, L>=50, "
            "a^2=exp(L)+1/(32L^2), ell=log(a), h=1/a, and "
            "h<exp(-25)<1/72000000000. Hence 0<ell<L and h^2<exp(-L)."
        ),
        "row_scale": (
            "R=||alpha||_1+||beta||_1, with alpha and beta the actual "
            "joined retained-observation rows. No numerical bound on R is assumed."
        ),
        "coarse_joined_bounds": {
            "ideal": (
                "|bar(gamma)_j|<=R and |bar(eta)_j|<=R; the Section "
                "11.175 errors give Gamma_j<2R and H_j<2R."
            ),
            "p_0": "|p_0-bar(gamma)_0|<8764h^2R",
            "p_1": "|p_1-bar(gamma)_1|<8765h^2R",
            "p_2": "|p_2-bar(gamma)_2|<8762h^2R",
            "p_3": "|p_3-bar(gamma)_3|<8760h^2R",
            "p_4": "|p_4|<(2/3)h^2R",
            "p_5": "|p_5|<h^2R/63",
        },
        "monotone_weight_bound": (
            "For every positive-coefficient polynomial W of degree at most "
            "five, L W'(L)<=5W(L), so exp(-L)W(L) decreases for L>=50. "
            "Therefore h^2 Khat_k(ell)<Khat_k(50)/(72000000000)^2."
        ),
        "weights_at_50": {
            f"Khat_{j}": value for j, value in enumerate(at_50)
        },
        "perturbation_polynomial_at_50": fraction_text(c_50),
        "perturbation_ratio": fraction_text(perturbation_ratio),
        "tail_separation": (
            "Let Pbar=sum_(j=0)^3bar(gamma)_j(lambda-ell)^j. Then "
            "|T[P_lin-Pbar]|<h^2 R/300000000000, while "
            "|T[Pbar]|<2h^2 sum_(j=0)^3 Khat_j(ell)|bar(gamma)_j|."
        ),
        "interpretation": (
            "All correction-induced coefficient departures, including the "
            "degree-four and degree-five channels, are negligible in the raw "
            "outside tail relative to the retained row scale. The unsuppressed "
            "ideal cubic is the remaining coefficient-aware tail theorem."
        ),
    }


def build_rows(symbolic: dict, bounds: dict) -> list[GateRow]:
    basis = symbolic["basis_translation"]
    weights = symbolic["translated_tail_weights"]
    conjugation = symbolic["physical_conjugation"]
    joined = symbolic["joined_operator"]
    coarse = bounds["coarse_joined_bounds"]
    return [
        GateRow(
            "pbc_01_domain",
            "physical_domain",
            "ready_to_apply",
            "The basis and conjugation gate stays on the complete q=1 chart.",
            bounds["domain"],
            "No extension outside the source chart is asserted.",
        ),
        GateRow(
            "pbc_02_forward",
            "exact_basis_map",
            "certified",
            "The centered coefficients have one exact lambda-basis image.",
            basis["forward"],
            "The signs and binomial factors are retained exactly.",
        ),
        GateRow(
            "pbc_03_inverse",
            "exact_basis_map",
            "certified",
            "The basis map is exactly invertible.",
            basis["inverse"],
            "No coefficient information is discarded by centering.",
        ),
        GateRow(
            "pbc_04_jet",
            "jet_interpretation",
            "certified",
            "Both coefficient vectors are Taylor jets of one polynomial.",
            basis["jet_interpretation"],
            "Large translated coefficients can be coordinate effects.",
        ),
        GateRow(
            "pbc_05_tail_weights",
            "translated_majorant",
            "certified",
            "The old tail template induces six explicit centered weights.",
            weights["definition"] + " " + "; ".join(weights["rows"]),
            "This is still a coefficientwise absolute majorant.",
        ),
        GateRow(
            "pbc_06_sharpness",
            "nonpromotion_guard",
            "guard_validated",
            "The powers of log(a) cannot be removed inside the old six-term majorant.",
            weights["sharpness"],
            "Grouped cross-coefficient cancellation may still improve the true tail.",
        ),
        GateRow(
            "pbc_07_parameter",
            "physical_parameter_identity",
            "certified",
            "The physical amplitude rate absorbs the coordinate origin.",
            conjugation["parameter"],
            "This uses the q=1 source identity, not a fitted parameter.",
        ),
        GateRow(
            "pbc_08_exponent",
            "exact_conjugation",
            "certified",
            "The Morse exponent has an exact centered form.",
            conjugation["exponent"] + " " + conjugation["log_derivative"],
            "The sampled coordinate x can still range over a long interval.",
        ),
        GateRow(
            "pbc_09_first_order",
            "exact_conjugation",
            "certified",
            "The first-order Morse amplitude has no explicit log(a) coefficient.",
            conjugation["first_order"],
            "This is an identity, not a mode-sum estimate.",
        ),
        GateRow(
            "pbc_10_saddle",
            "exact_conjugation",
            "certified",
            "The removable saddle operator is a fixed conjugated differential operator.",
            conjugation["saddle_operator"] + " " + conjugation["conjugate_form"],
            "No sign follows from this factorization alone.",
        ),
        GateRow(
            "pbc_11_mode",
            "coordinate_guard",
            "guard_validated",
            "Explicit log(a) coefficients cancel without erasing the physical mode range.",
            conjugation["mode_coordinate"],
            "The distinction prevents both fake losses and fake uniformity.",
        ),
        GateRow(
            "pbc_12_first_joined",
            "joined_row_identity",
            "certified",
            "The first-order amplitude remains two joined row contractions.",
            joined["field_vector"] + " " + joined["first_order"],
            "No alpha or beta component is bounded separately.",
        ),
        GateRow(
            "pbc_13_saddle_joined",
            "joined_row_identity",
            "certified",
            "The saddle operator preserves the joined alpha/beta structure.",
            joined["saddle"],
            "This identifies the exact vector kernels still requiring cancellation.",
        ),
        GateRow(
            "pbc_14_off_saddle_joined",
            "joined_row_identity",
            "certified",
            "The complete off-saddle c-prime numerator, after its common physical factor, has only two joined contractions.",
            joined["off_saddle"],
            "The roster sum and its phases remain unevaluated.",
        ),
        GateRow(
            "pbc_15_row_scale",
            "relative_envelope",
            "certified",
            "The Section 11.175 envelopes admit one coarse joined row scale.",
            bounds["row_scale"] + " " + " ".join(coarse.values()),
            "R is a bookkeeping scale, not a proved numerical observation bound.",
        ),
        GateRow(
            "pbc_16_tail_perturbation",
            "quantitative_tail_separation",
            "certified",
            "Every departure from the ideal cubic has a uniformly tiny raw outside-tail contribution.",
            bounds["tail_separation"] + " " + bounds["monotone_weight_bound"],
            "The fixed source normalizer and real projection are not absorbed here.",
        ),
        GateRow(
            "pbc_17_handoff",
            "open_quantitative_handoff",
            "open",
            "The remaining linear theorem is the ideal cubic inside the grouped physical functional.",
            bounds["interpretation"],
            "It still requires retained-row or direct phase-correlation control.",
        ),
        GateRow(
            "pbc_18_boundary",
            "proof_boundary",
            "guard_validated",
            "Basis cancellation is not promoted to a signed flow theorem.",
            (
                "There are zero numerical row bounds, grouped roster bounds, "
                "quadratic-pairing bounds, Phi_B bounds, or Xi signs here."
            ),
            "No contact exclusion, Lambda<=0, PF-infinity, RH, or prize theorem is claimed.",
        ),
    ]


def render_note(artifact: dict) -> str:
    symbolic = artifact["symbolic_certificate"]
    basis = symbolic["basis_translation"]
    weights = symbolic["translated_tail_weights"]
    conjugation = symbolic["physical_conjugation"]
    joined = symbolic["joined_operator"]
    bounds = artifact["bound_certificate"]
    coarse = bounds["coarse_joined_bounds"]
    return "\n".join(
        [
            "# Physical P_lin Centered-Basis Morse Conjugation Gate",
            "",
            "Date: 2026-08-02",
            "",
            "Status: exact centered-to-lambda translation, physical Morse conjugation, joined alpha/beta c-prime formulas, and quantitative outside-tail perturbation separation. This is not a proof of the grouped ideal-cubic remainder estimate, signed flow, or RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Exact Basis Map",
            "",
            "Put ell=log(a), x=lambda-ell, and",
            "",
            "```text",
            "P(lambda)=sum_(k=0)^5 p_k x^k=sum_(j=0)^5 q_j lambda^j.",
            basis["forward"],
            basis["inverse"],
            "```",
            "",
            basis["jet_interpretation"],
            "",
            "The six forward rows are",
            "",
            "```text",
            *basis["forward_rows"],
            "```",
            "",
            "## Tail Translation",
            "",
            weights["definition"],
            "",
            "```text",
            *weights["rows"],
            "```",
            "",
            "Therefore the old lambda-basis estimate gives",
            "",
            "```text",
            "|mathcal T[P]|<2h^2 sum_(k=0)^5 Khat_k(ell)|p_k|.",
            "```",
            "",
            weights["sharpness"],
            "",
            "## Physical Morse Conjugation",
            "",
            bounds["domain"],
            "",
            "```text",
            conjugation["parameter"],
            conjugation["exponent"],
            conjugation["log_derivative"],
            "```",
            "",
            "Thus every explicit ell coefficient cancels from the differential operator:",
            "",
            "```text",
            conjugation["first_order"],
            conjugation["saddle_operator"],
            conjugation["conjugate_form"],
            "```",
            "",
            conjugation["mode_coordinate"],
            "",
            "## Joined Rows",
            "",
            joined["field_vector"],
            "",
            "```text",
            joined["first_order"],
            joined["saddle"],
            "```",
            "",
            "For the full off-saddle numerator:",
            "",
            "```text",
            joined["off_saddle"],
            "```",
            "",
            "This is the direct grouped formulation promised in Section 11.175: no scalar current has been converted into an unproved row norm.",
            "",
            "## Quantitative Tail Separation",
            "",
            bounds["row_scale"],
            "",
            "The detailed Section 11.175 envelopes imply the safe coarse bounds",
            "",
            "```text",
            coarse["ideal"],
            coarse["p_0"],
            coarse["p_1"],
            coarse["p_2"],
            coarse["p_3"],
            coarse["p_4"],
            coarse["p_5"],
            "```",
            "",
            bounds["monotone_weight_bound"],
            "",
            "At L=50 the translated weights are",
            "",
            "```text",
            *[
                f"{key}={value}"
                for key, value in bounds["weights_at_50"].items()
            ],
            "```",
            "",
            bounds["tail_separation"],
            "",
            "The certified rational coefficient multiplying h^2 R is",
            "",
            "```text",
            bounds["perturbation_ratio"],
            "<1/300000000000.",
            "```",
            "",
            bounds["interpretation"],
            "",
            "## Next Target",
            "",
            "Estimate the ideal cubic with bar(gamma)_0,...,bar(gamma)_3 left inside the joined Morse/endpoint phase sum. The exact alpha/beta vector kernels above are the admissible starting point. A separate numerical norm on either retained row is optional, not assumed.",
            "",
            "The four quadratic observation pairings remain a separate obligation, and the fixed source normalizer nu and common phase c_xi must be restored before comparison with the signed flow benchmark.",
            "",
            "## Proof Boundary",
            "",
            "This gate proves no numerical retained-observation bound, grouped ideal-cubic c-prime sum, reciprocal cancellation theorem, quadratic residual bound, h^2 flow error, Phi_B upper bound, contact exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
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
    source = source_audit(payloads)
    symbolic = symbolic_certificate()
    bounds = bound_certificate(symbolic)
    rows = build_rows(symbolic, bounds)
    return {
        "kind": STEM,
        "date": "2026-08-02",
        "status": (
            "exact centered-to-lambda basis map, physical Morse conjugation, "
            "joined alpha/beta c-prime identities, and raw outside-tail "
            "perturbation separation complete; grouped ideal-cubic estimate open; "
            "no signed flow, Phi_B bound, Xi theorem, or RH"
        ),
        "proof_boundary": (
            "This gate proves exact basis and Morse-operator identities and a "
            "row-relative outside-tail perturbation bound. It proves no numerical "
            "retained-observation bound, grouped ideal-cubic c-prime sum, quadratic "
            "residual bound, signed flow estimate, Phi_B upper bound, contact "
            "exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, "
            "RH, or prize-level conclusion."
        ),
        "source_audit": source,
        "symbolic_certificate": symbolic,
        "bound_certificate": bounds,
        "counts": {
            "rows": len(rows),
            "basis_coefficients": 6,
            "inverse_basis_coefficients": 6,
            "translated_tail_weights": 6,
            "exact_morse_conjugations": 3,
            "joined_operator_identities": 3,
            "tail_perturbation_bounds": 1,
            "numerical_observation_row_bounds": 0,
            "grouped_ideal_cubic_bounds": 0,
            "quadratic_remainder_bounds": 0,
            "signed_flow_bounds": 0,
            "phi_b_bounds": 0,
        },
        "rows": [asdict(row) for row in rows],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()

    artifact = build_artifact()
    atomic_write(args.out, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(args.note, render_note(artifact))
    print(
        "built physical P_lin basis-conjugation gate: "
        f"{artifact['counts']['rows']} rows, "
        f"{artifact['counts']['basis_coefficients']} basis coefficients, "
        f"{artifact['counts']['translated_tail_weights']} translated weights, "
        f"{artifact['counts']['joined_operator_identities']} joined identities, "
        f"{artifact['counts']['grouped_ideal_cubic_bounds']} ideal-cubic bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
