#!/usr/bin/env python3
"""Build the centered bulk adjacent-pair transfer gate."""

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
    "jensen_window_pf_newman_polymath15_"
    "first_order_centered_bulk_pair_transfer_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "first_correction": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_"
        "dirichlet_first_correction_gate.json"
    ),
    "centered_reduction": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_real_residual_reduction.json"
    ),
    "adjacent_stability": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_adjacent_chart_stability_certificate.json"
    ),
}

L_MIN = 50
TL_MAX = 25
D_ABSOLUTE_CONSTANT = 2_189
D_DIFFERENCE_CONSTANT = 44
AMPLITUDE_DECAY_NUMERATOR = 2
AMPLITUDE_DECAY_DENOMINATOR = 5
PAIR_DETERMINANT_DENOMINATOR = 16
PAIR_SEPARATION_DENOMINATOR = 48


@dataclass(frozen=True)
class TransferRow:
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
        "first_correction": (
            "d_(t,n)(s)=1/(6s)+alpha'(s)*(t/4+t^2*alpha_n^2/8)",
            "m_[1](t,n;s)=m_(t,n)(s)*(1+d_(t,n)(s))",
        ),
        "centered_reduction": (
            "M_(1,a)=sum_(n=1)^N log(n/a)f_n",
            "|U-A_a|<1e-6*exp(-5L/4)",
        ),
        "adjacent_stability": (
            "|A_(a,N+1)-A_(a,N)|<5000exp(-7L/4)",
            "chart-invariant bulk",
        ),
    }
    for key, required in markers.items():
        text = json.dumps(payloads[key], sort_keys=True)
        for marker in required:
            if marker not in text:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {key: str(value.get("kind", "")) for key, value in payloads.items()}


def symbolic_audit() -> dict[str, str]:
    c, b, ell = sp.symbols("c b ell", real=True)
    x_part, y_part = sp.symbols("x_n y_n", real=True)
    component = sp.Matrix([x_part, y_part])
    transfer = sp.Matrix([[1, 0], [-c * ell, b * ell]])
    if sp.factor(transfer.det() - b * ell) != 0:
        raise RuntimeError("single-component transfer determinant failed")

    slope = c + sp.I * b
    complex_component = x_part + sp.I * y_part
    centered_real = sp.expand_complex(
        sp.re(-slope * ell * complex_component)
    )
    if sp.simplify(
        centered_real - ell * (-c * x_part + b * y_part)
    ) != 0:
        raise RuntimeError("single-component centered projection failed")

    rho, cosine, sine = sp.symbols(
        "rho cosine sine", positive=True, real=True
    )
    ell_next = sp.symbols("ell_next", real=True)
    rotation = sp.Matrix([[cosine, -sine], [sine, cosine]])
    transfer_next = sp.Matrix(
        [[1, 0], [-c * ell_next, b * ell_next]]
    )
    pair_matrix = transfer + rho * transfer_next * rotation
    pair_det = sp.expand(pair_matrix.det())
    expected_pair_det = (
        b * (ell + rho**2 * ell_next)
        + rho
        * (
            b * (ell + ell_next) * cosine
            + c * (ell_next - ell) * sine
        )
    )
    unit_circle_residual = sp.factor(
        pair_det
        - expected_pair_det
        - b * rho**2 * ell_next * (cosine**2 + sine**2 - 1)
    )
    if unit_circle_residual != 0:
        raise RuntimeError("locked-pair transfer determinant failed")

    big_l, big_m = sp.symbols("L_0 M_0", positive=True)
    factor_identity = sp.factor(
        (big_l + rho**2 * big_m) ** 2
        - rho**2 * (big_l + big_m) ** 2
    )
    expected_factor = (
        (1 - rho**2)
        * (big_l - rho * big_m)
        * (big_l + rho * big_m)
    )
    if sp.simplify(factor_identity - expected_factor) != 0:
        raise RuntimeError("determinant margin factorization failed")

    alpha_n, alpha_prime, time, h = sp.symbols(
        "alpha_n alpha_prime time h"
    )
    d_n = alpha_prime * time**2 * alpha_n**2 / 8
    d_next = alpha_prime * time**2 * (alpha_n - h) ** 2 / 8
    d_difference = sp.factor(d_next - d_n)
    expected_d_difference = (
        alpha_prime * time**2 * (h**2 - 2 * h * alpha_n) / 8
    )
    if sp.simplify(d_difference - expected_d_difference) != 0:
        raise RuntimeError("correction-ratio difference failed")

    f_1, f_2, f_3, f_4 = sp.symbols("f_1 f_2 f_3 f_4")
    e_1, e_2, e_3, e_4 = sp.symbols(
        "ell_1 ell_2 ell_3 ell_4"
    )
    partials = [
        f_1,
        f_1 + f_2,
        f_1 + f_2 + f_3,
        f_1 + f_2 + f_3 + f_4,
    ]
    abel = (
        e_4 * partials[3]
        - (e_2 - e_1) * partials[0]
        - (e_3 - e_2) * partials[1]
        - (e_4 - e_3) * partials[2]
    )
    direct = e_1 * f_1 + e_2 * f_2 + e_3 * f_3 + e_4 * f_4
    if sp.expand(abel - direct) != 0:
        raise RuntimeError("finite Abel basis identity failed")

    targets = (
        sp.Matrix([1, 0]),
        sp.Matrix([0, 1]),
        sp.Matrix([-1, -1]),
    )
    ell_1, ell_2, ell_3 = sp.symbols(
        "ell_1 ell_2 ell_3", nonzero=True, real=True
    )
    coefficients = []
    for target, ell_value in zip(targets, (ell_1, ell_2, ell_3)):
        coefficient = sp.Matrix(
            [
                target[0],
                (target[1] / ell_value + c * target[0]) / b,
            ]
        )
        coefficients.append(coefficient)
    counter_sum = sum(
        (
            sp.Matrix([[1, 0], [-c * ell_value, b * ell_value]])
            * coefficient
            for ell_value, coefficient in zip(
                (ell_1, ell_2, ell_3), coefficients
            )
        ),
        sp.zeros(2, 1),
    )
    if sp.simplify(counter_sum) != sp.zeros(2, 1):
        raise RuntimeError("three-component cancellation family failed")

    p_1, q_1, p_2, q_2 = sp.symbols(
        "p_1 q_1 p_2 q_2", positive=True, real=True
    )
    phase_x, phase_y = sp.symbols("phase_x phase_y", real=True)
    phase_vector = sp.Matrix([phase_x, phase_y])
    block_1 = sp.diag(p_1, q_1) * phase_vector
    block_2 = sp.diag(p_2, q_2) * phase_vector
    turning = sp.det(sp.Matrix.hstack(block_1, block_2))
    expected_turning = (
        (p_1 * q_2 - q_1 * p_2) * phase_x * phase_y
    )
    if sp.factor(turning - expected_turning) != 0:
        raise RuntimeError("common-phase turning countermodel failed")

    return {
        "reflected_basis": (
            "2*(X,A_a)^T=(g_0,G_a)^T+(conj(g_0),conj(G_a))^T"
            "+sum_n{(1,-s_*'*ell_n)^T*f_n"
            "+(1,-conj(s_*')*ell_n)^T*conj(f_n)}, "
            "ell_n=log(n/a)"
        ),
        "real_transfer": (
            "For s_*'=c+i*b and f_n=x_n+i*y_n, "
            "V_n=(X_n,A_n)^T=T_(ell_n)(x_n,y_n)^T, "
            "T_ell=[[1,0],[-c*ell,b*ell]], det(T_ell)=b*ell"
        ),
        "abel_basis": (
            "F_k=sum_(n<=k)f_n, h_k=log((k+1)/k): "
            "M_(1,a)=ell_N*F_N-sum_(k=1)^(N-1)h_k*F_k; "
            "A_a=Re(G_a-s_*'*ell_N*F_N"
            "+s_*'*sum_(k=1)^(N-1)h_k*F_k)"
        ),
        "correction_difference": (
            "d_(n+1)-d_n=alpha'*t^2*(h_n^2-2*h_n*alpha_n)/8, "
            "alpha_n=alpha-log(n)"
        ),
        "ratio_lock": (
            "R_n=f_(n+1)/f_n=rho_n*exp(i*delta_n)="
            "exp[-s_*h_n+(t/4)*(2log(n)h_n+h_n^2)]"
            "*(1+d_(n+1))/(1+d_n)"
        ),
        "pair_matrix": (
            "V_n+V_(n+1)=B_n*(Re f_n,Im f_n)^T, "
            "B_n=T_(ell_n)+rho_n*T_(ell_(n+1))*R(delta_n)"
        ),
        "pair_determinant": (
            "det(B_n)=b*(ell_n+rho_n^2*ell_(n+1))"
            "+rho_n*[b*(ell_n+ell_(n+1))*cos(delta_n)"
            "+c*(ell_(n+1)-ell_n)*sin(delta_n)]"
        ),
        "positive_coordinates": (
            "With b=-B, ell_n=-L_0, ell_(n+1)=-M_0, "
            "h_n=L_0-M_0>0: det(B_n)="
            "B*(L_0+rho_n^2*M_0)"
            "+rho_n*[B*(L_0+M_0)*cos(delta_n)"
            "+c*h_n*sin(delta_n)]"
        ),
        "margin_factorization": (
            "(L_0+rho^2*M_0)^2-rho^2*(L_0+M_0)^2"
            "=(1-rho^2)*(L_0-rho*M_0)*(L_0+rho*M_0)"
        ),
        "three_component_countermodel": (
            "For arbitrary c, b!=0, and ell_j!=0, any target "
            "w_j=(P_j,Q_j) is realized by "
            "f_j=P_j+i*(Q_j/ell_j+c*P_j)/b. "
            "The nonzero targets (1,0),(0,1),(-1,-1) sum to zero."
        ),
        "turning_form": (
            "For W_k=B_(2k-1)u_k and "
            "u_(k+1)=rho_(2k-1)rho_(2k)"
            "*R(delta_(2k-1)+delta_(2k))*u_k, "
            "det(W_k,W_(k+1)) is the quadratic form "
            "rho_(2k-1)rho_(2k)*u_k^T"
            "*B_(2k-1)^T*J*B_(2k+1)"
            "*R(delta_(2k-1)+delta_(2k))*u_k"
        ),
        "turning_countermodel": (
            "If the frame drift and all three ratio phases are zero, "
            "positive pair matrices may be "
            "diag(p_1,q_1),diag(p_2,q_2). Then "
            "det(W_1,W_2)=(p_1*q_2-q_1*p_2)*x*y, "
            "which takes both signs under the common phase whenever "
            "p_1*q_2!=q_1*p_2."
        ),
    }


def numeric_audit() -> dict[str, str]:
    x_min = 4 * math.pi * math.exp(L_MIN)
    d_absolute = D_ABSOLUTE_CONSTANT / x_min
    d_log_difference = (
        2 * D_DIFFERENCE_CONSTANT / x_min
    )
    uncorrected_per_h = -33 / 80 + 1 / (16 * x_min)
    corrected_per_h = uncorrected_per_h + d_log_difference
    if not d_absolute < 0.5:
        raise RuntimeError("absolute first-correction disk failed")
    if not corrected_per_h < -2 / 5:
        raise RuntimeError("corrected amplitude decay failed")

    h_max = math.log(2)
    one_minus_ratio_per_h = (
        1 - math.exp(-2 * h_max / 5)
    ) / h_max
    if not one_minus_ratio_per_h > 0.25:
        raise RuntimeError("one-minus-ratio margin failed")

    time_max = TL_MAX / L_MIN
    saddle = math.sqrt(
        x_min / (4 * math.pi) + time_max / 16
    )
    saddle_to_x_quarter = saddle / (x_min / 4)
    if not saddle_to_x_quarter < 1:
        raise RuntimeError("saddle/frame comparison failed")

    return {
        "x_min": format(x_min, ".17e"),
        "d_absolute_upper": format(d_absolute, ".17e"),
        "d_log_ratio_per_h_upper": format(
            d_log_difference, ".17e"
        ),
        "uncorrected_log_ratio_per_h_upper": format(
            uncorrected_per_h, ".17e"
        ),
        "corrected_log_ratio_per_h_upper": format(
            corrected_per_h, ".17e"
        ),
        "one_minus_ratio_over_h_lower": format(
            one_minus_ratio_per_h, ".17e"
        ),
        "saddle_over_x_quarter": format(
            saddle_to_x_quarter, ".17e"
        ),
    }


def build_exact() -> dict:
    symbolic = symbolic_audit()
    numeric = numeric_audit()
    return {
        "domain": (
            "L=log(x/(4*pi))>=50, 0<=tL<=25, "
            "a^2=x/(4*pi)+t/16, N=floor(a), 1<=n<=N-1"
        ),
        "reflected_sharp_basis": symbolic["reflected_basis"],
        "real_component_transfer": symbolic["real_transfer"],
        "abel_chart_basis": symbolic["abel_basis"],
        "ratio_lock": symbolic["ratio_lock"],
        "correction_difference": symbolic["correction_difference"],
        "amplitude_monotonicity": (
            "|d_n|<2189/x<1/2, "
            "|d_(n+1)-d_n|<44*h_n/x, "
            "log(rho_n)<-(2/5)h_n, hence "
            "0<rho_n<=exp(-2h_n/5)<1 and "
            "1-rho_n>h_n/4"
        ),
        "pair_matrix": symbolic["pair_matrix"],
        "pair_determinant": symbolic["pair_determinant"],
        "positive_coordinates": symbolic["positive_coordinates"],
        "determinant_bound": (
            "C_x>0, D_x>0, B=1/2+t*C_x/4>=1/2, "
            "c=t*D_x/4<=1/(4x)<=h_n/16. "
            "Therefore det(B_n)>=h_n*[B*(1-rho_n)-c*rho_n]"
            ">=h_n^2/16."
        ),
        "two_term_separation": (
            "||B_n||_F<3L, so "
            "||(X_n+X_(n+1),A_n+A_(n+1))||_2"
            ">h_n^2*|f_n|/(48L)."
        ),
        "chart_edge": (
            "The disjoint adjacent blocks stop before the unpaired edge. "
            "The certified recurrence "
            "|A_(a,N+1)-A_(a,N)|<5000exp(-7L/4) "
            "must be retained for that last-saddle/endpoint block."
        ),
        "three_component_countermodel": symbolic[
            "three_component_countermodel"
        ],
        "locked_block_countermodel": (
            "For any three invertible locked-pair matrices B_1,B_2,B_3, "
            "set u_j=B_j^(-1)w_j for "
            "w_1=(1,0), w_2=(0,1), w_3=(-1,-1). "
            "Every pair obeys its internal ratio lock and is nonzero, "
            "yet sum_j B_j u_j=0. Cross-block Xi ratios are essential."
        ),
        "turning_form": symbolic["turning_form"],
        "turning_countermodel": symbolic["turning_countermodel"],
        "surviving_xi_target": (
            "Use the full consecutive ratio chain together with the exact "
            "normalizer/endpoint phase anchor, not isolated pair "
            "invertibility or relative phases alone, to prove a quantitative "
            "half-plane, signed-area, or one-sided winding bound for the "
            "paired block polygon plus the recurrent edge block on the "
            "contact value band."
        ),
        "numeric_audit": numeric,
        "factorization": symbolic["margin_factorization"],
    }


def build_rows(exact: dict) -> list[TransferRow]:
    return [
        TransferRow(
            "nfocbptg_01_reflected_basis",
            "exact_reduction",
            "available_exact",
            "The centered pair has one division-free reflected sharp basis.",
            exact["reflected_sharp_basis"],
            "This is a real projection identity, not a lower bound.",
        ),
        TransferRow(
            "nfocbptg_02_abel_basis",
            "exact_reduction",
            "available_exact",
            "Finite Abel summation isolates chart-independent positive logarithmic increments.",
            exact["abel_chart_basis"],
            "No sign is asserted for the complex partial sums.",
        ),
        TransferRow(
            "nfocbptg_03_component_transfer",
            "exact_lemma",
            "available_exact",
            "Each interior complex coefficient maps locally and orientation-preservingly to value/slope coordinates.",
            exact["real_component_transfer"],
            "The determinant degenerates at ell=0 and says nothing about sums of components.",
        ),
        TransferRow(
            "nfocbptg_04_ratio_lock",
            "exact_lemma",
            "available_exact",
            "Consecutive corrected Xi coefficients obey an exact phase-amplitude ratio law.",
            exact["ratio_lock"],
            "The law is local; its global oscillatory consequences remain to be proved.",
        ),
        TransferRow(
            "nfocbptg_05_amplitude_decay",
            "asymptotic_theorem",
            "certified",
            "Corrected coefficient amplitudes strictly decrease across every interior adjacent pair.",
            exact["amplitude_monotonicity"],
            "Valid only on the displayed L>=50 first-order domain.",
            exact["numeric_audit"],
        ),
        TransferRow(
            "nfocbptg_06_pair_determinant",
            "exact_lemma",
            "available_exact",
            "The locked two-frequency pair has a common-phase-independent transfer determinant.",
            exact["pair_determinant"],
            "Nonzero pair determinant is a local two-term statement.",
        ),
        TransferRow(
            "nfocbptg_07_positive_pair",
            "asymptotic_theorem",
            "certified",
            "Every interior locked adjacent pair has a uniformly positive transfer determinant.",
            exact["determinant_bound"],
            "This does not imply that different pair outputs have a common direction.",
        ),
        TransferRow(
            "nfocbptg_08_pair_separation",
            "asymptotic_theorem",
            "certified",
            "Each locked adjacent two-term block has an explicit first-jet lower separation.",
            exact["two_term_separation"],
            "The bound is before summing other blocks and the endpoint.",
        ),
        TransferRow(
            "nfocbptg_09_edge_composition",
            "composition_rule",
            "available_exact",
            "The unpaired last saddle is composed through the certified adjacent recurrence.",
            exact["chart_edge"],
            "It must not be bounded as an independent exp(-3L/4) component.",
        ),
        TransferRow(
            "nfocbptg_10_component_countermodel",
            "countermodel",
            "guard_validated",
            "Orientation-preserving component transfers alone do not prevent exact aggregate cancellation.",
            exact["three_component_countermodel"],
            "The coefficients are generic, not the Xi coefficient chain.",
        ),
        TransferRow(
            "nfocbptg_11_pair_countermodel",
            "countermodel",
            "guard_validated",
            "Positive locked-pair invertibility and relative ratio phases do not give aggregate nonvanishing or a fixed turning sign.",
            (
                exact["locked_block_countermodel"]
                + " "
                + exact["turning_countermodel"]
            ),
            "The surviving theorem must use cross-block compatibility and the absolute Xi phase anchor.",
        ),
        TransferRow(
            "nfocbptg_12_xi_handoff",
            "open_theorem_target",
            "not_ready_to_apply",
            "The next arithmetic input must use the full phase-locked Xi chain to control cancellation between blocks.",
            exact["surviving_xi_target"],
            "No anchored half-plane, signed-area, crossing-count, or winding theorem is proved here.",
        ),
        TransferRow(
            "nfocbptg_13_scalar_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The q>=1 contact-band scalar lower bound and one-sided successor count remain open.",
            (
                "|X|<=50000exp(-5L/4) => "
                "|A_a|>(100000L+1)exp(-5L/4), followed by an "
                "A_a-positive successor count below one turn"
            ),
            "The local pair theorem is not this global implication.",
        ),
        TransferRow(
            "nfocbptg_14_nonpromotion",
            "nonpromotion_guard",
            "guard_validated",
            "The local determinant theorem is kept separate from contact exclusion and RH.",
            "adjacent locked-pair separation != aggregate Xi separation",
            (
                "No q<1 theorem, finite phase closure, contact exclusion, "
                "Lambda<=0, PF-infinity, RH, or Clay-prize conclusion."
            ),
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    rows = build_rows(exact)
    return {
        "kind": STEM,
        "date": "2026-07-26",
        "status": (
            "exact reflected bulk basis and uniform adjacent locked-pair "
            "transfer theorem on L>=50, 0<=tL<=25; aggregate Xi "
            "cancellation remains open"
        ),
        "proof_boundary": (
            "This artifact proves the reflected and Abel coordinate "
            "identities, the corrected consecutive ratio law, strict "
            "amplitude decrease, the exact locked-pair determinant, and "
            "det(B_n)>=h_n^2/16 with a two-term separation bound. It also "
            "gives exact countermodels showing that component or isolated-"
            "pair invertibility does not prove aggregate nonvanishing. It "
            "does not prove the q>=1 bulk scalar lower bound or signed "
            "crossing count, the q<1 multiplicity-compatible theorem, "
            "finite phase cells, one-sided successor winding, contact "
            "exclusion, Lambda<=0, RH, PF-infinity, or a Clay-prize "
            "conclusion."
        ),
        "sources": {
            key: str(path.relative_to(REPO_ROOT)).replace("\\", "/")
            for key, path in SOURCES.items()
        },
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "constants": {
            "L_min": L_MIN,
            "tL_max": TL_MAX,
            "d_absolute": D_ABSOLUTE_CONSTANT,
            "d_difference": D_DIFFERENCE_CONSTANT,
            "amplitude_decay": "2/5",
            "pair_determinant_denominator": PAIR_DETERMINANT_DENOMINATOR,
            "pair_separation_denominator": PAIR_SEPARATION_DENOMINATOR,
        },
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    numeric = exact["numeric_audit"]
    return "\n".join(
        [
            "# Newman First-Order Centered Bulk Pair-Transfer Gate",
            "",
            "Date: 2026-07-26",
            "",
            "Status: exact reflected/Abel bulk coordinates and a uniform",
            "adjacent locked-pair determinant theorem. The aggregate Xi",
            "cancellation problem remains open; this is not a proof of RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Reflected Bulk Coordinates",
            "",
            "On",
            "",
            "```text",
            exact["domain"],
            "```",
            "",
            "put `ell_n=log(n/a)`. The centered pair has the exact",
            "division-free reflected-sharp form",
            "",
            "```text",
            exact["reflected_sharp_basis"],
            "```",
            "",
            "For `s_*'=c+i*b` and `f_n=x_n+i*y_n`, the real",
            "component transfer is",
            "",
            "```text",
            exact["real_component_transfer"],
            "```",
            "",
            "Here `b=-1/2-t*C_x/4<0`; for every strict interior",
            "`ell_n<0`, the determinant `b*ell_n` is positive. This is",
            "only a local orientation statement.",
            "",
            "Finite Abel summation gives a second exact bulk basis:",
            "",
            "```text",
            exact["abel_chart_basis"],
            "```",
            "",
            "The weights `h_k=log(1+1/k)` are positive and independent",
            "of the saddle chart. The chart dependence is confined to",
            "the terminal/endpoint block already controlled by the adjacent",
            "stability theorem.",
            "",
            "## Exact Ratio Lock",
            "",
            "Consecutive corrected coefficients are not arbitrary:",
            "",
            "```text",
            exact["ratio_lock"],
            exact["correction_difference"],
            "```",
            "",
            "The imported critical-frame bounds give",
            "",
            "```text",
            exact["amplitude_monotonicity"],
            "```",
            "",
            f"At `L=50`, `x={numeric['x_min']}` and the checked upper",
            f"bound for `|d_n|` is `{numeric['d_absolute_upper']}`.",
            "The complete corrected logarithmic ratio per `h_n` is at most",
            f"`{numeric['corrected_log_ratio_per_h_upper']}<-2/5`.",
            "",
            "## Locked-Pair Determinant",
            "",
            "Let `R(delta)` be the real rotation matrix. Then",
            "",
            "```text",
            exact["pair_matrix"],
            exact["pair_determinant"],
            "```",
            "",
            "Writing `b=-B`, `ell_n=-L_0`,",
            "`ell_(n+1)=-M_0`, and `h_n=L_0-M_0`, this becomes",
            "",
            "```text",
            exact["positive_coordinates"],
            "```",
            "",
            "The phase can be eliminated without estimating it:",
            "",
            "```text",
            exact["determinant_bound"],
            "```",
            "",
            "Thus every adjacent corrected pair is orientation-positive",
            "for every common phase. The determinant also yields",
            "",
            "```text",
            exact["two_term_separation"],
            "```",
            "",
            "This is a genuine two-frequency first-jet theorem. It holds",
            "on the full displayed first-order domain, hence also on its",
            "`q=2tL^2>=1` sublayer.",
            "",
            "## Falsification Gate",
            "",
            "Local orientation is not aggregate separation. For generic",
            "component coefficients one has the exact family",
            "",
            "```text",
            exact["three_component_countermodel"],
            "```",
            "",
            "Even preserving each internal adjacent ratio is insufficient:",
            "",
            "```text",
            exact["locked_block_countermodel"],
            "```",
            "",
            "The latter construction breaks only the ratios *between*",
            "successive pairs. Keeping those relative ratios still leaves a",
            "common absolute phase. In exact matrix notation,",
            "",
            "```text",
            exact["turning_form"],
            "```",
            "",
            "and the generic common-phase guard is",
            "",
            "```text",
            exact["turning_countermodel"],
            "```",
            "",
            "Thus neither another isolated determinant estimate nor",
            "relative ratio phases alone can supply the orientation. The",
            "normalizer/endpoint phase anchor is part of the indispensable",
            "Xi input.",
            "",
            "## Xi-Specific Handoff",
            "",
            "```text",
            exact["surviving_xi_target"],
            "```",
            "",
            "A concrete next calculation is to substitute the actual",
            "normalizer phase and recurrent endpoint block into the displayed",
            "turning form, then reduce it in the saddle variable. Its purpose",
            "is to decide whether the anchored stationary-phase block polygon",
            "admits a quantitative half-plane or one-sided winding theorem.",
            "The unanchored sign is now rigorously rejected, so it must not be",
            "rediscovered by a large numerical sweep.",
            "",
            "The unpaired last saddle remains composed through",
            "",
            "```text",
            exact["chart_edge"],
            "```",
            "",
            "This artifact does not prove the aggregate contact-band scalar",
            "bound, its signed crossing count, the `q<1` chart, finite phase",
            "cells, contact exclusion, `Lambda<=0`, PF-infinity, RH, or a",
            "Clay-prize conclusion.",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2) + "\n", encoding="utf-8"
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman centered bulk pair-transfer gate: "
        "14 rows, rho_n<exp(-2h_n/5), det(B_n)>=h_n^2/16, "
        "1 exact aggregate countermodel family, 2 open Xi obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
