#!/usr/bin/env python3
"""Build the branch-free absolute-phase anchor reduction."""

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
    "first_order_centered_absolute_phase_anchor_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "endpoint_phase": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "endpoint_holomorphic_lift.json"
    ),
    "centered_reduction": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_real_residual_reduction.json"
    ),
    "wronskian_reduction": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_wronskian_crossing_reduction.json"
    ),
    "adjacent_stability": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_adjacent_chart_stability_certificate.json"
    ),
    "bulk_pair_gate": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_bulk_pair_transfer_gate.json"
    ),
}

L_MIN = 50
D_ABSOLUTE_CONSTANT = 2_189
VALUE_CHART_CONSTANT = 1_100
SCALAR_CHART_CONSTANT = 5_000
NORMALIZED_VALUE_CHART_CONSTANT = 2_200
NORMALIZED_SCALAR_CHART_CONSTANT = 10_000


@dataclass(frozen=True)
class AnchorRow:
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
        "endpoint_phase": (
            "M_0(i*T)*U(T)*exp(pi*i/8)",
            "T^(3/2)+i*T^(1/2)",
        ),
        "centered_reduction": (
            "M_(1,a)=sum_(n=1)^N log(n/a)f_n",
            "G_a=-kappa_N*(H_(a,x)+mu_a*H_a)",
            "|U-A_a|<1e-6*exp(-5L/4)",
        ),
        "wronskian_reduction": (
            "W_[1]=Im(E_[1],x*conj(E_[1]))",
            "If E_[1]=0 then X=Y=W_[1]=0",
        ),
        "adjacent_stability": (
            "|Delta X|<1100*exp(-5L/4)",
            "|A_(a,N+1)-A_(a,N)|<5000*exp(-7L/4)",
        ),
        "bulk_pair_gate": (
            "|d_n|<2189/x<1/2",
            "normalizer/endpoint phase anchor",
        ),
    }
    for key, required in markers.items():
        text = json.dumps(payloads[key], sort_keys=True)
        for marker in required:
            if marker not in text:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {key: str(value.get("kind", "")) for key, value in payloads.items()}


def symbolic_audit() -> dict[str, str]:
    f_r, f_i = sp.symbols("f_r f_i", real=True)
    z_0r, z_0i, z_ar, z_ai = sp.symbols(
        "z_0r z_0i z_ar z_ai", real=True
    )
    f_vector = sp.Matrix([f_r, f_i])
    shape_matrix = sp.Matrix(
        [[z_0r, -z_0i], [z_ar, -z_ai]]
    )
    projected = shape_matrix * f_vector
    expected_projected = sp.Matrix(
        [
            sp.re(
                (f_r + sp.I * f_i) * (z_0r + sp.I * z_0i)
            ),
            sp.re(
                (f_r + sp.I * f_i) * (z_ar + sp.I * z_ai)
            ),
        ]
    )
    if sp.simplify(projected - expected_projected) != sp.zeros(2, 1):
        raise RuntimeError("anchored real-projection matrix failed")

    determinant = sp.factor(shape_matrix.det())
    expected_determinant = -sp.im(
        (z_0r - sp.I * z_0i) * (z_ar + sp.I * z_ai)
    )
    if sp.simplify(determinant - expected_determinant) != 0:
        raise RuntimeError("anchored shape determinant failed")

    x, y, a_part, b_part = sp.symbols(
        "X Y A B", real=True
    )
    u, v, d_r, d_i = sp.symbols(
        "u v D_r D_i", real=True
    )
    e_value = x + sp.I * y
    c_value = a_part + sp.I * b_part
    d_value = d_r + sp.I * d_i
    derivative = (u + sp.I * v) * e_value + c_value + d_value
    wronskian = sp.expand_complex(
        sp.im(derivative * sp.conjugate(e_value))
    )
    unscaled_determinant = sp.expand_complex(
        -sp.im(c_value * sp.conjugate(e_value))
    )
    collapse = (
        v * (x**2 + y**2)
        + sp.im(d_value * sp.conjugate(e_value))
        - wronskian
    )
    if sp.simplify(unscaled_determinant - collapse) != 0:
        raise RuntimeError("Wronskian-collapse identity failed")
    if sp.simplify(
        unscaled_determinant.subs(x, 0) - a_part * y
    ) != 0:
        raise RuntimeError("crossing determinant identity failed")

    parity, beta, time_0, h_value = sp.symbols(
        "parity beta T_0 H_a", real=True
    )
    endpoint_core = -beta * (time_0 + sp.I)
    kappa = parity * endpoint_core
    g_value = -kappa * h_value
    expected_g = parity * beta * (time_0 + sp.I) * h_value
    if sp.simplify(g_value - expected_g) != 0:
        raise RuntimeError("endpoint phase simplification failed")

    m_r, m_i, norm, d_1 = sp.symbols(
        "M_r M_i norm d_1", nonzero=True
    )
    normalizer = m_r + sp.I * m_i
    first_coefficient = normalizer * (1 + d_1) / norm
    relative_endpoint = sp.cancel(g_value / norm / first_coefficient)
    expected_relative = sp.cancel(
        expected_g / (normalizer * (1 + d_1))
    )
    if sp.simplify(relative_endpoint - expected_relative) != 0:
        raise RuntimeError("relative endpoint normalization failed")

    return {
        "first_coefficient": (
            "f_1=phi*(1+d_1), |f_1|=|1+d_1|, "
            "|d_1|<2189/x<1/2"
        ),
        "unit_anchor": (
            "eta=f_1/|f_1|=phi*(1+d_1)/|1+d_1|"
        ),
        "relative_coefficients": (
            "q_n=f_n/f_1="
            "exp[t*log(n)^2/4-s_*log(n)]"
            "*(1+d_n)/(1+d_1), q_1=1"
        ),
        "endpoint_phase": (
            "beta=(sqrt(pi)/8)*exp(t*pi^2/64-pi*T_0/4)"
            "*T_0^(1/2)>0; "
            "kappa_N=(-1)^(N+1)*beta*(T_0+i)/|M_t(s)|; "
            "g_0=(-1)^N*beta*(T_0+i)*H_a/|M_t(s)|; "
            "G_a=(-1)^N*beta*(T_0+i)"
            "*(H_(a,x)+mu_a*H_a)/|M_t(s)|"
        ),
        "relative_endpoint": (
            "r_0=g_0/f_1=(-1)^N*beta*(T_0+i)*H_a"
            "/[M_t(s)*(1+d_1)]; "
            "r_A=G_a/f_1=(-1)^N*beta*(T_0+i)"
            "*(H_(a,x)+mu_a*H_a)/[M_t(s)*(1+d_1)]"
        ),
        "aggregate_shapes": (
            "P_0=sum_(n=1)^N q_n, "
            "P_1=sum_(n=1)^N ell_n*q_n, ell_n=log(n/a); "
            "Z_0=P_0+r_0=E_[1]/f_1; "
            "Z_A=-s_*'*P_1+r_A=C_a/f_1, "
            "C_a=-s_*'*M_(1,a)+G_a, A_a=Re(C_a)"
        ),
        "real_projection": (
            "(X,A_a)^T=S(Z_0,Z_A)*(Re(f_1),Im(f_1))^T, "
            "S=[[Re(Z_0),-Im(Z_0)],"
            "[Re(Z_A),-Im(Z_A)]]; "
            "X/|f_1|=Re(eta*Z_0), "
            "A_a/|f_1|=Re(eta*Z_A)"
        ),
        "contact_split": (
            "If Z_0!=0, then X=A_a=0 iff "
            "Re(eta*Z_0)=0 and Im(conj(Z_0)*Z_A)=0. "
            "If Z_0=0, then X=A_a=0 iff Re(eta*Z_A)=0."
        ),
        "shape_determinant": (
            "Delta_anchor=det(S)=-Im(conj(Z_0)*Z_A); "
            "|f_1|^2*Delta_anchor=-Im(C_a*conj(E_[1]))"
        ),
        "singular_value_bound": (
            "||(X,A_a)||_2>=|f_1|*|Delta_anchor|"
            "/sqrt(|Z_0|^2+|Z_A|^2) "
            "when |Z_0|^2+|Z_A|^2>0"
        ),
        "wronskian_collapse": (
            "E_[1],x=lambda_a*E_[1]+C_a+D_(1,x), "
            "lambda_a=u_a+i*v_a; hence "
            "|f_1|^2*Delta_anchor="
            "v_a*|E_[1]|^2+Im(D_(1,x)*conj(E_[1]))-W_[1]"
        ),
        "crossing_sign": (
            "At X=0, |f_1|^2*Delta_anchor=A_a*Y. "
            "If E_[1]=iY!=0, then "
            "A_a>0 iff Delta_anchor*Y>0."
        ),
        "complex_zero_guard": (
            "At E_[1]=0 one has Z_0=Delta_anchor=W_[1]=0, "
            "while A_a=|f_1|*Re(eta*Z_A) can be nonzero. "
            "The determinant cannot classify this crossing."
        ),
    }


def numeric_audit() -> dict[str, str]:
    x_min = 4 * math.pi * math.exp(L_MIN)
    d_upper = D_ABSOLUTE_CONSTANT / x_min
    first_lower = 1 - d_upper
    value_normalized = VALUE_CHART_CONSTANT / first_lower
    scalar_normalized = SCALAR_CHART_CONSTANT / first_lower
    if not d_upper < 0.5:
        raise RuntimeError("first-coefficient nonzero disk failed")
    if not value_normalized < NORMALIZED_VALUE_CHART_CONSTANT:
        raise RuntimeError("normalized value chart budget failed")
    if not scalar_normalized < NORMALIZED_SCALAR_CHART_CONSTANT:
        raise RuntimeError("normalized scalar chart budget failed")
    return {
        "x_min": format(x_min, ".17e"),
        "d_1_upper": format(d_upper, ".17e"),
        "f_1_lower": format(first_lower, ".17e"),
        "normalized_value_chart_constant_raw": format(
            value_normalized, ".17e"
        ),
        "normalized_scalar_chart_constant_raw": format(
            scalar_normalized, ".17e"
        ),
    }


def build_exact() -> dict:
    symbolic = symbolic_audit()
    numeric = numeric_audit()
    return {
        "domain": (
            "L=log(x/(4*pi))>=50, 0<=tL<=25, "
            "a^2=x/(4*pi)+t/16, N=floor(a)"
        ),
        "first_coefficient": symbolic["first_coefficient"],
        "unit_anchor": symbolic["unit_anchor"],
        "relative_coefficients": symbolic["relative_coefficients"],
        "endpoint_phase": symbolic["endpoint_phase"],
        "relative_endpoint": symbolic["relative_endpoint"],
        "aggregate_shapes": symbolic["aggregate_shapes"],
        "real_projection": symbolic["real_projection"],
        "contact_split": symbolic["contact_split"],
        "shape_determinant": symbolic["shape_determinant"],
        "singular_value_bound": symbolic["singular_value_bound"],
        "wronskian_collapse": symbolic["wronskian_collapse"],
        "crossing_sign": symbolic["crossing_sign"],
        "complex_zero_guard": symbolic["complex_zero_guard"],
        "chart_covariance": (
            "The anchor f_1 is independent of the prescribed cutoff. "
            "|Re(eta*Delta Z_0)|<2200*exp(-5L/4) and "
            "|Re(eta*Delta Z_A)|<10000*exp(-7L/4). "
            "The adjacent theorem controls these real projections, not "
            "the complex shape increments; Q_N can have a large imaginary "
            "part. Therefore Delta_anchor has no certified adjacent-chart "
            "absolute bound."
        ),
        "live_q_ge_1_target": (
            "On q=2tL^2>=1 in the canonical N=floor(a) chart, prove "
            "|Re(eta*Z_0)|<=50000*exp(-5L/4)/|f_1| implies "
            "|Re(eta*Z_A)|>(100000L+1)*exp(-5L/4)/|f_1|, "
            "then bound N_(Re(eta*Z_0)=0,Re(eta*Z_A)>0) "
            "below one successor turn. This formulation retains Z_0=0."
        ),
        "ordinary_exceptional_split": (
            "A possible proof split is: on |Z_0|>=tau use anchored "
            "phase-difference information together with the exact "
            "determinant collapse; on |Z_0|<tau prove the direct real "
            "projection Re(eta*Z_A) bound. The second branch is essential "
            "and cannot be inferred from Delta_anchor."
        ),
        "q_lt_1_target": (
            "Construct the separate q=2tL^2<1 "
            "multiplicity-compatible local chart and finite phase-cell "
            "theorem; the present first-order anchor does not provide it."
        ),
        "numeric_audit": numeric,
    }


def build_rows(exact: dict) -> list[AnchorRow]:
    return [
        AnchorRow(
            "nfocaapr_01_first_coefficient",
            "asymptotic_lemma",
            "certified",
            "The first corrected Xi coefficient is uniformly nonzero.",
            exact["first_coefficient"],
            "Valid on the displayed first-order L>=50 domain.",
            exact["numeric_audit"],
        ),
        AnchorRow(
            "nfocaapr_02_branch_free_anchor",
            "exact_definition",
            "ready_to_apply",
            "The coefficient chain has a canonical unit phase without choosing an argument branch.",
            exact["unit_anchor"] + "; " + exact["relative_coefficients"],
            "The anchor uses the fixed n=1 coefficient.",
        ),
        AnchorRow(
            "nfocaapr_03_endpoint_phase",
            "exact_identity",
            "ready_to_apply",
            "The endpoint normalizer phase reduces to the explicit algebraic direction T_0+i.",
            exact["endpoint_phase"],
            "Positive T_0 and the established principal branches.",
        ),
        AnchorRow(
            "nfocaapr_04_relative_endpoint",
            "exact_reduction",
            "ready_to_apply",
            "Division by f_1 puts the finite chain and endpoint in one absolute phase frame.",
            exact["relative_endpoint"],
            "M_t(s) and 1+d_1 are nonzero on the working domain.",
        ),
        AnchorRow(
            "nfocaapr_05_aggregate_shapes",
            "exact_reduction",
            "ready_to_apply",
            "Two relative complex shapes encode the retained value and centered scalar.",
            exact["aggregate_shapes"],
            "No lower bound is asserted.",
        ),
        AnchorRow(
            "nfocaapr_06_real_projection",
            "exact_identity",
            "ready_to_apply",
            "The original real pair is exactly the anchored projection of the two shapes.",
            exact["real_projection"],
            "This is branch-free but not phase-independent.",
        ),
        AnchorRow(
            "nfocaapr_07_contact_split",
            "exact_equivalence",
            "ready_to_apply",
            "Ordinary contacts split into anchor orthogonality and shape collinearity, while complex-main zeros form a separate branch.",
            exact["contact_split"],
            "The Z_0=0 clause must not be divided away.",
        ),
        AnchorRow(
            "nfocaapr_08_shape_determinant",
            "exact_identity",
            "ready_to_apply",
            "The aggregate shape determinant measures relative value/slope orientation.",
            exact["shape_determinant"],
            "A nonzero determinant is only a sufficient contact exclusion.",
        ),
        AnchorRow(
            "nfocaapr_09_singular_value",
            "exact_inequality",
            "ready_to_apply",
            "A determinant lower bound would imply quantitative first-jet separation away from the zero denominator case.",
            exact["singular_value_bound"],
            "This stronger sufficient route is not yet proved arithmetically.",
        ),
        AnchorRow(
            "nfocaapr_10_wronskian_collapse",
            "exact_obstruction",
            "guard_validated",
            "The aggregate determinant is the centered Wronskian plus already isolated frame and derivative-correction terms.",
            exact["wronskian_collapse"],
            "It is not an independent replacement for the direct scalar theorem.",
        ),
        AnchorRow(
            "nfocaapr_11_crossing_sign",
            "exact_equivalence",
            "ready_to_apply",
            "On an ordinary real-part crossing the determinant carries the centered crossing sign.",
            exact["crossing_sign"],
            "Requires E_[1] nonzero for the sign equivalence.",
        ),
        AnchorRow(
            "nfocaapr_12_complex_zero_guard",
            "nondivision_guard",
            "guard_validated",
            "The determinant and Wronskian both lose the complex-main-zero crossing class.",
            exact["complex_zero_guard"],
            "Direct control of Re(eta*Z_A) remains necessary.",
        ),
        AnchorRow(
            "nfocaapr_13_chart_guard",
            "covariance_guard",
            "guard_validated",
            "Only the anchored real projections inherit the adjacent-chart certificate.",
            exact["chart_covariance"],
            "The complex determinant is a canonical-chart surrogate, not a certified chart-invariant scalar.",
        ),
        AnchorRow(
            "nfocaapr_14_q_ge_1_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The live q>=1 obligation is the direct anchored real-projection lower bound and its one-sided crossing count.",
            (
                exact["live_q_ge_1_target"]
                + " "
                + exact["ordinary_exceptional_split"]
            ),
            "No such arithmetic lower bound or count is proved here.",
        ),
        AnchorRow(
            "nfocaapr_15_q_lt_1_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The q<1 multiplicity-compatible chart remains a separate obligation.",
            exact["q_lt_1_target"],
            "It is not implied by the q>=1 first-order reduction.",
        ),
        AnchorRow(
            "nfocaapr_16_nonpromotion",
            "nonpromotion_guard",
            "guard_validated",
            "The exact phase reduction is kept separate from the missing global arithmetic theorem and RH.",
            (
                "absolute phase anchor and determinant collapse "
                "!= anchored projection lower bound"
            ),
            (
                "No q>=1 scalar theorem, q<1 theorem, finite phase closure, "
                "contact exclusion, Lambda<=0, PF-infinity, RH, or "
                "Clay-prize conclusion."
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
            "exact branch-free first-coefficient phase anchor and endpoint "
            "normalization; the aggregate determinant collapses to the "
            "Wronskian and does not remove the direct real-projection "
            "obligation"
        ),
        "proof_boundary": (
            "This artifact proves the nonvanishing first-coefficient "
            "anchor, explicit endpoint phase, relative aggregate shapes, "
            "contact split, determinant identity, Wronskian collapse, "
            "ordinary crossing sign, and real-projection chart covariance. "
            "It proves that the determinant route is not independent and "
            "does not classify E_[1]=0. It does not prove the q>=1 anchored "
            "projection lower bound or crossing count, the q<1 "
            "multiplicity-compatible theorem, finite phase cells, "
            "one-sided successor winding, contact exclusion, Lambda<=0, "
            "RH, PF-infinity, or a Clay-prize conclusion."
        ),
        "sources": {
            key: str(path.relative_to(REPO_ROOT)).replace("\\", "/")
            for key, path in SOURCES.items()
        },
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "constants": {
            "L_min": L_MIN,
            "d_absolute": D_ABSOLUTE_CONSTANT,
            "value_chart": VALUE_CHART_CONSTANT,
            "scalar_chart": SCALAR_CHART_CONSTANT,
            "normalized_value_chart": NORMALIZED_VALUE_CHART_CONSTANT,
            "normalized_scalar_chart": NORMALIZED_SCALAR_CHART_CONSTANT,
        },
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    numeric = exact["numeric_audit"]
    return "\n".join(
        [
            "# Newman First-Order Centered Absolute-Phase Anchor Reduction",
            "",
            "Date: 2026-07-26",
            "",
            "Status: exact branch-free phase normalization and a decisive",
            "determinant-collapse guard. The direct anchored projection",
            "theorem remains open; this is not a proof of RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## First-Coefficient Anchor",
            "",
            "The n=1 term removes every logarithmic phase:",
            "",
            "```text",
            exact["first_coefficient"],
            exact["unit_anchor"],
            "```",
            "",
            f"At `L=50`, the certified `|d_1|` upper bound is",
            f"`{numeric['d_1_upper']}`, so the anchor is uniformly",
            "nonzero without choosing a branch of `arg`.",
            "",
            "Relative coefficients are therefore exact:",
            "",
            "```text",
            exact["relative_coefficients"],
            "```",
            "",
            "## Endpoint In The Same Frame",
            "",
            "The established endpoint phase cancellation gives",
            "",
            "```text",
            exact["endpoint_phase"],
            "```",
            "",
            "After division by the first coefficient, the absolute",
            "normalizer phase cancels algebraically:",
            "",
            "```text",
            exact["relative_endpoint"],
            "```",
            "",
            "Define the two aggregate shapes by",
            "",
            "```text",
            exact["aggregate_shapes"],
            "```",
            "",
            "Then the retained real value and centered scalar are exactly",
            "",
            "```text",
            exact["real_projection"],
            "```",
            "",
            "## Exact Contact Split",
            "",
            "The contact system now has a branch-free geometric form:",
            "",
            "```text",
            exact["contact_split"],
            "```",
            "",
            "The second clause is essential. It retains the complex-main-",
            "zero class that any division by `Z_0` would delete.",
            "",
            "For the ordinary branch, the natural shape determinant is",
            "",
            "```text",
            exact["shape_determinant"],
            exact["singular_value_bound"],
            "```",
            "",
            "## Determinant Collapse",
            "",
            "The determinant is not a new independent arithmetic route.",
            "Substituting the exact centered derivative decomposition gives",
            "",
            "```text",
            exact["wronskian_collapse"],
            "```",
            "",
            "At a real-part crossing this reduces to",
            "",
            "```text",
            exact["crossing_sign"],
            "```",
            "",
            "but at a complex-main zero,",
            "",
            "```text",
            exact["complex_zero_guard"],
            "```",
            "",
            "Thus a lower bound for `Delta_anchor` is a stronger sufficient",
            "Wronskian theorem with the same blind spot. It cannot replace",
            "the direct centered scalar.",
            "",
            "## Chart-Covariance Guard",
            "",
            "```text",
            exact["chart_covariance"],
            "```",
            "",
            "This matters because the adjacent cancellation is real, not",
            "complex-absolute. Treating `Delta_anchor` as chart invariant",
            "would discard the cancellation already proved in the endpoint",
            "recurrence.",
            "",
            "## Surviving Theorem",
            "",
            "The exact q>=1 target is now",
            "",
            "```text",
            exact["live_q_ge_1_target"],
            "```",
            "",
            "A logically complete proof may split the domain as follows:",
            "",
            "```text",
            exact["ordinary_exceptional_split"],
            "```",
            "",
            "The independent q<1 obligation is",
            "",
            "```text",
            exact["q_lt_1_target"],
            "```",
            "",
            "This reduction proves no anchored projection lower bound,",
            "signed crossing budget, q<1 chart, finite phase closure,",
            "contact exclusion, `Lambda<=0`, PF-infinity, RH, or",
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
        "built Newman centered absolute-phase anchor reduction: "
        "16 rows, 1 branch-free anchor, 1 determinant collapse, "
        "2 open arithmetic obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
