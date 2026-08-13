#!/usr/bin/env python3
"""Build the ideal-cubic nilpotent jet and mode-transport gate."""

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
    "physical_plin_ideal_cubic_nilpotent_transport_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "basis_conjugation": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
        "physical_plin_basis_conjugation_gate.json"
    ),
    "endpoint_abel": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
        "endpoint_coherence_abel_gate.json"
    ),
    "q1_saddle": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_moment_q1_saddle_"
        "phase_variance_reduction.json"
    ),
    "physical_plin": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
        "physical_plin_amplitude_gate.json"
    ),
}


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


def require_zero(expression: sp.Expr | sp.MatrixBase, label: str) -> None:
    if isinstance(expression, sp.MatrixBase):
        if any(sp.simplify(sp.expand(item)) != 0 for item in expression):
            raise RuntimeError(label)
        return
    if sp.simplify(sp.expand(expression)) != 0:
        raise RuntimeError(label)


def matrix_one_norm(matrix: sp.MatrixBase) -> sp.Expr:
    return max(
        sp.simplify(sum(abs(matrix[i, j]) for i in range(matrix.rows)))
        for j in range(matrix.cols)
    )


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
    basis = payloads["basis_conjugation"]
    if basis.get("counts", {}).get("joined_operator_identities") != 3:
        raise RuntimeError("basis-conjugation parent drifted")
    joined = basis.get("symbolic_certificate", {}).get("joined_operator", {})
    if "F=(C,mathcal N,A,Q)" not in joined.get("field_vector", ""):
        raise RuntimeError("joined field vector drifted")
    if "L_B[P_lin]" not in joined.get("saddle", ""):
        raise RuntimeError("joined saddle parent drifted")

    abel = payloads["endpoint_abel"]
    phase = abel.get("phase_certificate", {})
    for key, marker in (
        ("first_difference", "log(1+1/r)"),
        ("second_difference", "alpha log((r+1)^2/[r(r+2)])"),
        ("curvature_envelope", "5/(4alpha)"),
    ):
        if marker not in phase.get(key, ""):
            raise RuntimeError(f"endpoint Abel phase drifted: {key}")
    if "3sqrt(alpha)" not in abel.get("abel_certificate", {}).get(
        "partial_sum_bound", ""
    ):
        raise RuntimeError("endpoint Abel partial-sum bound drifted")

    q1 = payloads["q1_saddle"].get("physical_q1_certificate", {})
    if "t=1/(2L^2)" not in q1.get("domain", ""):
        raise RuntimeError("q=1 t law drifted")
    if "B_a>49/100" not in q1.get("b_a_bounds", ""):
        raise RuntimeError("q=1 B_a bound drifted")

    ideal = payloads["physical_plin"].get("ideal_certificate", {})
    if "R=Q=i*x/2" not in ideal.get("specialization", ""):
        raise RuntimeError("physical ideal specialization drifted")
    if "x^3/4" not in ideal.get("ideal_polynomial", ""):
        raise RuntimeError("physical ideal cubic source drifted")

    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def symbolic_certificate() -> dict:
    i = sp.I
    x, u_x = sp.symbols("x u_x", real=True)
    identity_4 = sp.eye(4)
    zero_4 = sp.zeros(4)
    d_0 = sp.Matrix(
        [
            [0, 0, 0, 0],
            [0, 0, i / 2, i / 2],
            [i / 2, 0, 0, 0],
            [i / 2, 0, 0, 0],
        ]
    )
    f_at_zero = sp.Matrix([1, -i * u_x / 2, 0, 0])
    field = sp.Matrix([1, -x**2 / 4 - i * u_x / 2, i * x / 2, i * x / 2])
    require_zero(d_0**3, "four-vector generator is not cubic nilpotent")
    require_zero(sp.diff(field, x) - d_0 * field, "four-vector jet failed")
    require_zero(
        field - (identity_4 + x * d_0 + x**2 * d_0**2 / 2) * f_at_zero,
        "four-vector exponential failed",
    )

    d_hat = d_0.row_join(zero_4).col_join(identity_4.row_join(d_0))
    identity_8 = sp.eye(8)
    y_at_zero = f_at_zero.col_join(sp.zeros(4, 1))
    augmented = field.col_join(x * field)
    require_zero(d_hat**4, "augmented generator is not order-four nilpotent")
    if d_hat**3 == sp.zeros(8):
        raise RuntimeError("augmented nilpotence order collapsed")
    require_zero(
        sp.diff(augmented, x) - d_hat * augmented,
        "augmented jet failed",
    )
    augmented_exp = sum(
        ((x**k / sp.factorial(k)) * d_hat**k for k in range(4)),
        sp.zeros(8),
    )
    require_zero(
        augmented - augmented_exp * y_at_zero,
        "augmented exponential failed",
    )

    alpha = sp.symbols("alpha_V alpha_N alpha_A alpha_Q", real=True)
    beta = sp.symbols("beta_V beta_N beta_A beta_Q", real=True)
    row = sp.Matrix([[*alpha, *[i * item for item in beta]]])
    ideal_polynomial = sp.expand((row * augmented)[0])
    ideal_coefficients = [
        sp.Poly(ideal_polynomial, x).coeff_monomial(x**j) for j in range(4)
    ]
    expected_coefficients = [
        alpha[0] - i * u_x * alpha[1] / 2,
        i * (alpha[2] + alpha[3]) / 2
        + i * beta[0]
        + u_x * beta[1] / 2,
        -alpha[1] / 4 - (beta[2] + beta[3]) / 2,
        -i * beta[1] / 4,
    ]
    for index, (actual, expected) in enumerate(
        zip(ideal_coefficients, expected_coefficients)
    ):
        require_zero(actual - expected, f"ideal cubic coefficient {index}")
    anti = sp.Matrix([[0, 0, 1, -1]])
    require_zero((anti * field)[0], "A-Q null direction failed")

    t, delta_a = sp.symbols("t delta_a", real=True)
    epsilon = t * (x + delta_a) / 2
    g = -sp.Rational(1, 2) + epsilon
    h_scalar = g**2 + g + t / 2 + sp.Rational(1, 6)
    require_zero(
        h_scalar - (epsilon**2 + t / 2 - sp.Rational(1, 12)),
        "physical saddle scalar failed",
    )
    m_hat = (
        d_hat**2
        - identity_8 / 12
        + 2 * epsilon * d_hat
        + (epsilon**2 + t / 2) * identity_8
    )
    direct_operator = (
        sp.diff(augmented, x, 2)
        + (2 * g + 1) * sp.diff(augmented, x)
        + h_scalar * augmented
    )
    require_zero(direct_operator - m_hat * augmented, "augmented saddle matrix")
    m_core = d_hat**2 - identity_8 / 12
    core_determinant = sp.factor(m_core.det())
    if core_determinant != sp.Rational(1, 12**8):
        raise RuntimeError("ideal saddle core determinant drifted")

    d = sp.symbols("d", nonnegative=True, real=True)
    transport = sum(
        (((-d) ** k / sp.factorial(k)) * d_hat**k for k in range(4)),
        sp.zeros(8),
    )
    require_zero(
        augmented.subs(x, x - d) - transport * augmented,
        "augmented mode transport failed",
    )
    d_1, d_2 = sp.symbols("d_1 d_2", nonnegative=True, real=True)

    def transport_at(step: sp.Expr) -> sp.Matrix:
        return sum(
            (
                ((-step) ** k / sp.factorial(k)) * d_hat**k
                for k in range(4)
            ),
            sp.zeros(8),
        )

    require_zero(
        transport_at(d_2) * transport_at(d_1) - transport_at(d_1 + d_2),
        "transport semigroup failed",
    )
    second_transport = sp.expand(
        transport_at(d_1 + d_2) - 2 * transport_at(d_1) + identity_8
    )
    expected_second = (
        (d_1 - d_2) * d_hat
        + ((d_1 + d_2) ** 2 / 2 - d_1**2) * d_hat**2
        + (2 * d_1**3 - (d_1 + d_2) ** 3) * d_hat**3 / 6
    )
    require_zero(second_transport - expected_second, "second transport failed")

    epsilon_symbol = sp.symbols("epsilon", real=True)

    def mode_matrix(e_value: sp.Expr) -> sp.Matrix:
        return (
            d_hat**2
            - identity_8 / 12
            + 2 * e_value * d_hat
            + (e_value**2 + t / 2) * identity_8
        )

    next_epsilon = epsilon_symbol - t * d / 2
    saddle_first_transport = sp.expand(
        mode_matrix(next_epsilon) * transport - mode_matrix(epsilon_symbol)
    )
    second_next_epsilon = epsilon_symbol - t * (d_1 + d_2) / 2
    saddle_second_transport = sp.expand(
        mode_matrix(second_next_epsilon) * transport_at(d_1 + d_2)
        - 2
        * mode_matrix(epsilon_symbol - t * d_1 / 2)
        * transport_at(d_1)
        + mode_matrix(epsilon_symbol)
    )

    norms = {
        "D_0": matrix_one_norm(d_0),
        "D_0_squared": matrix_one_norm(d_0**2),
        "D_hat": matrix_one_norm(d_hat),
        "D_hat_squared": matrix_one_norm(d_hat**2),
        "D_hat_cubed": matrix_one_norm(d_hat**3),
    }

    return {
        "four_vector_jet": {
            "field": (
                "F_0(x)=(1,-x^2/4-i*u_(N,x)/2,i*x/2,i*x/2)^T."
            ),
            "generator": (
                "D_0=[[0,0,0,0],[0,0,i/2,i/2],"
                "[i/2,0,0,0],[i/2,0,0,0]]."
            ),
            "jet": (
                "F_0'=D_0F_0, D_0^3=0, and "
                "F_0(x)=exp(xD_0)F_0(0)."
            ),
            "null_direction": (
                "The row direction (0,0,1,-1) annihilates F_0(x) "
                "identically, so only alpha_A+alpha_Q and beta_A+beta_Q survive."
            ),
        },
        "augmented_jet": {
            "definition": "Y(x)=(F_0(x),xF_0(x))^T.",
            "generator": (
                "Dhat=[[D_0,0],[I_4,D_0]], with Dhat^4=0 but Dhat^3!=0."
            ),
            "jet": (
                "Y'=Dhat Y and Y(x)=exp(xDhat)Y(0), where the exponential "
                "terminates after its cubic term."
            ),
            "row": (
                "bar(P)=(alpha^T,i beta^T)Y; its four coefficients are "
                "the Section 11.175 bar(gamma)_0,...,bar(gamma)_3."
            ),
            "ranks": [int((d_hat**k).rank()) for k in range(1, 4)],
        },
        "saddle_matrix": {
            "epsilon": (
                "epsilon_g=g_tilde+1/2=t(x+delta_a)/2. On the physical "
                "amplitude interval |epsilon_g|<1/200 and t<=1/5000."
            ),
            "matrix": (
                "M(epsilon_g)=Dhat^2-I_8/12+2epsilon_g Dhat+"
                "(epsilon_g^2+t/2)I_8, and L_B[Y]=M(epsilon_g)Y."
            ),
            "core": (
                "At t=epsilon_g=0, M_0=Dhat^2-I_8/12 and "
                "det(M_0)=12^(-8), so nilpotence supplies transport but no "
                "universal saddle zero."
            ),
            "core_determinant": str(core_determinant),
        },
        "mode_transport": {
            "step": (
                "d_r=log(1+1/r), x_(r+1)=x_r-d_r, and "
                "epsilon_(r+1)=epsilon_r-t*d_r/2."
            ),
            "transport": (
                "E(d)=exp(-dDhat)=I-dDhat+d^2Dhat^2/2-d^3Dhat^3/6, "
                "so Y_(r+1)=E(d_r)Y_r."
            ),
            "semigroup": "E(d_2)E(d_1)=E(d_1+d_2) exactly.",
            "second_difference": (
                "Delta^2Y_r=[(d_r-d_(r+1))Dhat+"
                "(((d_r+d_(r+1))^2)/2-d_r^2)Dhat^2+"
                "(2d_r^3-(d_r+d_(r+1))^3)Dhat^3/6]Y_r."
            ),
            "saddle_first": (
                "Delta[M_rY_r]=[M(epsilon_r-td_r/2)E(d_r)-"
                "M(epsilon_r)]Y_r."
            ),
            "saddle_second": (
                "Delta^2[M_rY_r]=[M(epsilon_r-t(d_r+d_(r+1))/2)"
                "E(d_r+d_(r+1))-2M(epsilon_r-td_r/2)E(d_r)+"
                "M(epsilon_r)]Y_r."
            ),
            "symbolic_matrix_zero_checks": 9,
            "first_transport_nonzero_entries": sum(
                1 for item in saddle_first_transport if item != 0
            ),
            "second_transport_nonzero_entries": sum(
                1 for item in saddle_second_transport if item != 0
            ),
        },
        "matrix_norms": {key: str(value) for key, value in norms.items()},
    }


def bound_certificate(symbolic: dict) -> dict:
    alpha_min = 4096
    first_at_min = (
        Fraction(20, 9)
        + Fraction(5, 4) * Fraction(100, 81 * alpha_min)
        + Fraction(1, 4) * Fraction(1000, 729 * alpha_min**2)
    )
    if not first_at_min < Fraction(9, 4):
        raise RuntimeError("first transport collar rounding failed")

    second_at_min = (
        Fraction(5, 2)
        + Fraction(500, 81)
        + Fraction(5, 2) * Fraction(1000, 729 * alpha_min)
    )
    if not second_at_min < 9:
        raise RuntimeError("second transport collar rounding failed")

    perturbation_norm = Fraction(4, 200) + Fraction(1, 40000) + Fraction(
        1, 10000
    )
    if perturbation_norm != Fraction(161, 8000):
        raise RuntimeError("physical matrix perturbation rounding failed")

    return {
        "physical_interval": (
            "On 0<=lambda<=log(B)<ell with ell=log(a)<L, "
            "delta_a=ell-Re(alpha_s), and t=1/(2L^2), one has "
            "|x+delta_a|=|lambda-Re(alpha_s)|<L. Hence "
            "|epsilon_g|<tL/2=1/(4L)<=1/200 and t<=1/5000."
        ),
        "operator_perturbation": (
            "In the induced matrix one-norm, ||Dhat||_1=2. Therefore "
            "||M(epsilon_g)-(Dhat^2-I_8/12)||_1<"
            "4/200+1/40000+1/10000=161/8000."
        ),
        "phase_step_link": (
            "On the terminal collar, d_r=log(1+1/r), "
            "Delta G=1-alpha*d_r, and "
            "d_r-d_(r+1)=Delta^2G/alpha. Thus "
            "0<d_r-d_(r+1)<=5/(4alpha^2), while "
            "d_r<=exp(1/10)/alpha<10/(9alpha)."
        ),
        "first_difference_bound": (
            "Using ||Dhat||_1=2, ||Dhat^2||_1=5/2, and "
            "||Dhat^3||_1=3/2, the exact transport gives "
            "||Y_(r+1)-Y_r||_1<(9/(4alpha))||Y_r||_1."
        ),
        "second_difference_bound": (
            "The exact second-transport polynomial and the phase-curvature "
            "bounds give ||Y_(r+2)-2Y_(r+1)+Y_r||_1<"
            "(9/alpha^2)||Y_r||_1."
        ),
        "rational_diagnostics": {
            "first_difference_scaled_at_alpha_4096": fraction_text(
                first_at_min
            ),
            "first_difference_target": "9/4",
            "second_difference_scaled_at_alpha_4096": fraction_text(
                second_at_min
            ),
            "second_difference_target": "9",
            "operator_perturbation": fraction_text(perturbation_norm),
        },
        "abel_handoff": (
            "The ideal cubic field now has first variation O(alpha^-1) and "
            "second variation O(alpha^-2) on exactly the collar where the "
            "phase curvature is O(alpha^-1). Insert the exact transported "
            "saddle matrix into the weighted Abel identity before bounding "
            "the Fresnel integral or endpoint packages."
        ),
        "nonpromotion": (
            "The invertible core M_0 and arbitrary joined rows prevent a "
            "universal zero or sign. Slow mode transport is an Abel input, "
            "not a grouped roster bound or signed flow theorem."
        ),
    }


def build_rows(symbolic: dict, bounds: dict) -> list[GateRow]:
    field = symbolic["four_vector_jet"]
    augmented = symbolic["augmented_jet"]
    saddle = symbolic["saddle_matrix"]
    transport = symbolic["mode_transport"]
    return [
        GateRow("icnt_01_field", "ideal_field", "certified", "The surviving ideal cubic is generated by one four-vector polynomial.", field["field"], "The correction-driven difference remains governed by Section 11.176."),
        GateRow("icnt_02_generator", "nilpotent_jet", "certified", "The four-vector field has one constant nilpotent generator.", field["generator"] + " " + field["jet"], "This is an algebraic transport identity, not a phase estimate."),
        GateRow("icnt_03_null", "row_quotient", "certified", "The ideal channel identifies the A and Q field directions.", field["null_direction"], "The physical corrections need not preserve this exact null direction."),
        GateRow("icnt_04_augmented", "joined_nilpotent_jet", "certified", "The base and lifted fields form one order-four nilpotent jet.", augmented["definition"] + " " + augmented["generator"] + " " + augmented["jet"], "Both retained rows remain inside one augmented row."),
        GateRow("icnt_05_cubic", "joined_cubic", "certified", "The augmented row reproduces all four ideal coefficients.", augmented["row"], "No numerical retained-row bound is inserted."),
        GateRow("icnt_06_epsilon", "physical_operator_parameter", "certified", "The physical Morse operator is a small perturbation of one constant nilpotent core.", saddle["epsilon"] + " " + bounds["physical_interval"], "Small operator perturbation does not make the row contraction small."),
        GateRow("icnt_07_matrix", "saddle_matrix", "certified", "The removable saddle operator is one explicit matrix polynomial.", saddle["matrix"], "The full c-prime integral still contains its Fresnel and endpoint geometry."),
        GateRow("icnt_08_core", "nonpromotion_guard", "guard_validated", "The ideal saddle core is invertible rather than identically canceling.", saddle["core"], "Nilpotent transport supplies no universal saddle zero or sign."),
        GateRow("icnt_09_perturbation", "matrix_envelope", "certified", "The physical saddle matrix stays uniformly close to its constant core.", bounds["operator_perturbation"], "This is a matrix coefficient bound, not a roster sum."),
        GateRow("icnt_10_step", "mode_lattice", "certified", "Adjacent logarithmic modes have one exact step law.", transport["step"], "The fixed roster and endpoint crossings remain unchanged."),
        GateRow("icnt_11_transport", "exact_mode_transport", "certified", "Every ideal joined field shifts by a terminating cubic exponential.", transport["transport"] + " " + transport["semigroup"], "No numerical mode enumeration is required."),
        GateRow("icnt_12_phase_link", "phase_transport_bridge", "certified", "The variation of the logarithmic step is exactly the normalized phase curvature.", bounds["phase_step_link"], "This bridge is local to the certified terminal collar."),
        GateRow("icnt_13_first", "variation_bound", "certified", "The ideal joined field has O(alpha^-1) first variation on the collar.", bounds["first_difference_bound"], "The field norm remains explicit."),
        GateRow("icnt_14_second", "variation_bound", "certified", "The ideal joined field has O(alpha^-2) second variation on the collar.", transport["second_difference"] + " " + bounds["second_difference_bound"], "This does not by itself improve the existing phase partial-sum theorem."),
        GateRow("icnt_15_saddle_first", "transported_saddle", "certified", "The physical saddle field has an exact first-difference matrix.", transport["saddle_first"], "Its joined row contraction and Fresnel integral are not bounded here."),
        GateRow("icnt_16_saddle_second", "transported_saddle", "certified", "The physical saddle field has an exact second-difference matrix.", transport["saddle_second"], "Endpoint packages must be recomposed before an Abel estimate."),
        GateRow("icnt_17_handoff", "abel_handoff", "open", "The exact transport is ready for a coefficient-aware Abel estimate.", bounds["abel_handoff"], "No grouped ideal-cubic roster bound has yet been proved."),
        GateRow("icnt_18_boundary", "proof_boundary", "guard_validated", "Nilpotent smoothness is not promoted to the RH bridge.", bounds["nonpromotion"], "No quadratic residual bound, signed flow, Phi_B theorem, contact exclusion, Lambda<=0, PF-infinity, RH, or prize theorem is claimed."),
    ]


def render_note(artifact: dict) -> str:
    symbolic = artifact["symbolic_certificate"]
    field = symbolic["four_vector_jet"]
    augmented = symbolic["augmented_jet"]
    saddle = symbolic["saddle_matrix"]
    transport = symbolic["mode_transport"]
    bounds = artifact["bound_certificate"]
    return "\n".join(
        [
            "# Physical P_lin Ideal-Cubic Nilpotent Transport Gate",
            "",
            "Date: 2026-08-02",
            "",
            "Status: exact ideal four-vector and augmented nilpotent jets, physical saddle matrix, logarithmic mode transport, phase-curvature bridge, and first/second collar variation bounds. This is not a proof of the grouped roster estimate, signed flow, or RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Four-Vector Jet",
            "",
            "```text",
            field["field"],
            field["generator"],
            field["jet"],
            "```",
            "",
            field["null_direction"],
            "",
            "## Joined Augmentation",
            "",
            "```text",
            augmented["definition"],
            augmented["generator"],
            augmented["jet"],
            augmented["row"],
            "```",
            "",
            "The ranks of Dhat, Dhat^2, and Dhat^3 are respectively "
            + ", ".join(str(value) for value in augmented["ranks"])
            + ".",
            "",
            "## Physical Saddle Matrix",
            "",
            bounds["physical_interval"],
            "",
            "```text",
            saddle["epsilon"],
            saddle["matrix"],
            saddle["core"],
            "```",
            "",
            bounds["operator_perturbation"],
            "",
            "## Exact Mode Transport",
            "",
            "```text",
            transport["step"],
            transport["transport"],
            transport["semigroup"],
            "```",
            "",
            bounds["phase_step_link"],
            "",
            "The exact second field difference is",
            "",
            "```text",
            transport["second_difference"],
            "```",
            "",
            "and the certified collar bounds are",
            "",
            "```text",
            bounds["first_difference_bound"],
            bounds["second_difference_bound"],
            "```",
            "",
            "For the physical saddle field itself, no finite differencing or mode enumeration is needed:",
            "",
            "```text",
            transport["saddle_first"],
            transport["saddle_second"],
            "```",
            "",
            "## Handoff",
            "",
            bounds["abel_handoff"],
            "",
            "The exact identity d_r-d_(r+1)=Delta^2G/alpha ties the second variation of the coefficient jet to the already certified phase curvature. This is the structural input needed for a second-order or blockwise Abel calculation; it does not itself perform that calculation.",
            "",
            bounds["nonpromotion"],
            "",
            "## Proof Boundary",
            "",
            "This gate proves no numerical retained-observation bound, grouped ideal-cubic c-prime or roster estimate, endpoint cancellation, quadratic residual bound, h^2 flow error, Phi_B upper bound, contact exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
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
            "exact ideal four-vector and augmented nilpotent jets, physical "
            "saddle matrix, logarithmic mode transport, phase-curvature bridge, "
            "and first/second collar variation bounds complete; grouped roster "
            "estimate open; no signed flow, Phi_B theorem, Xi theorem, or RH"
        ),
        "proof_boundary": (
            "This gate proves nilpotent ideal-cubic transport and collar "
            "variation estimates. It proves no numerical retained-observation "
            "bound, grouped ideal-cubic c-prime or roster estimate, endpoint "
            "cancellation, quadratic residual bound, signed flow estimate, "
            "Phi_B upper bound, contact exclusion, retained aggregate or Xi "
            "theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
        "source_audit": source,
        "symbolic_certificate": symbolic,
        "bound_certificate": bounds,
        "counts": {
            "rows": len(rows),
            "four_vector_nilpotence_order": 3,
            "augmented_nilpotence_order": 4,
            "ideal_cubic_coefficients": 4,
            "invisible_ideal_row_directions": 2,
            "exact_mode_transport_identities": 5,
            "collar_variation_bounds": 2,
            "grouped_ideal_cubic_bounds": 0,
            "numerical_observation_row_bounds": 0,
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
    counts = artifact["counts"]
    print(
        "built ideal-cubic nilpotent transport gate: "
        f"{counts['rows']} rows, nilpotence orders "
        f"{counts['four_vector_nilpotence_order']}/"
        f"{counts['augmented_nilpotence_order']}, "
        f"{counts['ideal_cubic_coefficients']} cubic coefficients, "
        f"{counts['exact_mode_transport_identities']} transport identities, "
        f"{counts['collar_variation_bounds']} collar bounds, "
        f"{counts['grouped_ideal_cubic_bounds']} grouped bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
