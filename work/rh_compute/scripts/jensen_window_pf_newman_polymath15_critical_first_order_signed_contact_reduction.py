#!/usr/bin/env python3
"""Build the first-order signed-contact reduction for the critical layer."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_critical_"
    "first_order_signed_contact_reduction"
)
DEFAULT_OUT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / f"outputs/{STEM}.md"
SOURCES = {
    "global_remainder": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_"
        "first_order_global_remainder_certificate.json"
    ),
    "first_correction": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_"
        "dirichlet_first_correction_gate.json"
    ),
    "contact_normal": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_"
        "parabolic_frequency_contact_normal_hierarchy_gate.json"
    ),
    "oscillatory_handoff": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "oscillatory_zeta_handoff_theorem.json"
    ),
}
C_STAR = "4911678521/1933561194"
ETA_0 = 100_000
ETA_1 = 200_000


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_hashes() -> dict[str, str]:
    missing = [str(path) for path in SOURCES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing source artifacts: {missing}")
    return {key: file_hash(path) for key, path in SOURCES.items()}


def verify_symbolic_identities() -> None:
    s, time = sp.symbols("s t", nonzero=True)
    alpha = sp.Function("alpha")
    alpha_value = alpha(s)
    alpha_first = sp.symbols("alpha_first")
    alpha_second = sp.symbols("alpha_second")
    log_n = sp.symbols("ell", real=True)
    alpha_n = alpha_value - log_n
    moment = time / 4 + time**2 * alpha_n**2 / 8
    correction = 1 / (6 * s) + alpha_first * moment
    correction_s = sp.diff(
        1 / (6 * s)
        + sp.diff(alpha_value, s)
        * (time / 4 + time**2 * alpha_n**2 / 8),
        s,
    ).subs(
        {
            sp.diff(alpha_value, s): alpha_first,
            sp.diff(alpha_value, s, 2): alpha_second,
        }
    )
    correction_x = sp.expand(-sp.I * correction_s / 2)
    expected_correction_x = sp.simplify(
        -sp.I
        / 2
        * (
            -1 / (6 * s**2)
            + alpha_second * moment
            + time**2 * alpha_n * alpha_first**2 / 4
        )
    )
    if sp.simplify(correction_x - expected_correction_x) != 0:
        raise RuntimeError("first-correction derivative identity failed")

    j_zero, delta_j, remainder = sp.symbols("J_0 DeltaJ r_1")
    old_remainder = delta_j + remainder
    j_one = j_zero + delta_j
    if sp.expand(j_zero + old_remainder - (j_one + remainder)) != 0:
        raise RuntimeError("old/new signed repartition failed")

    mass_plus, mass_minus = sp.symbols("M_plus M_minus", real=True)
    h_plus, h_minus, defect = sp.symbols(
        "h_plus h_minus D_0", real=True
    )
    mass = (mass_plus + mass_minus) / 2
    value = mass_plus - mass_minus
    slope = mass_plus * h_plus - mass_minus * h_minus + defect
    expected_slope = (
        mass * (h_plus - h_minus)
        + value * (h_plus + h_minus) / 2
        + defect
    )
    if sp.expand(slope - expected_slope) != 0:
        raise RuntimeError("rate-free near-crossing identity failed")

    rem_value, rem_slope = sp.symbols("r r_x", real=True)
    abstract_value, abstract_slope = sp.symbols("X U", real=True)
    signed_signal = abstract_slope - abstract_value * (
        h_plus + h_minus
    ) / 2
    contact_signal = signed_signal.subs(
        {
            abstract_slope: -rem_slope / 2,
            abstract_value: -rem_value / 2,
        }
    )
    expected_contact = -rem_slope / 2 + rem_value * (
        h_plus + h_minus
    ) / 4
    if sp.simplify(contact_signal - expected_contact) != 0:
        raise RuntimeError("signed contact identity failed")

    u_ref, v_ref, x_real, y_imag, r_real = sp.symbols(
        "u_ref v_ref X Y R", real=True
    )
    reference_real_slope = u_ref * x_real - v_ref * y_imag + r_real
    centered_contact = sp.simplify(
        (r_real - v_ref * y_imag).subs(
            r_real,
            -rem_slope / 2
            - u_ref * (-rem_value / 2)
            + v_ref * y_imag,
        )
    )
    if (
        sp.simplify(
            centered_contact - (-rem_slope + u_ref * rem_value) / 2
        )
        != 0
    ):
        raise RuntimeError("centered contact identity failed")
    if (
        sp.simplify(
            reference_real_slope.subs(
                {
                    x_real: -rem_value / 2,
                    r_real: (
                        -rem_slope + u_ref * rem_value
                    )
                    / 2
                    + v_ref * y_imag,
                }
            )
            + rem_slope / 2
        )
        != 0
    ):
        raise RuntimeError("centered reference reconstruction failed")

    tau = sp.symbols("tau", positive=True)
    y = sp.symbols("y")
    root = sp.sqrt(2 * tau)
    toy = y**2 - 2 * tau
    if sp.simplify(sp.diff(toy, y).subs(y, root) - 2 * root) != 0:
        raise RuntimeError("quadratic boundary-splitting guard failed")


def exact_identities() -> dict[str, str]:
    return {
        "region": (
            "L=log(x/(4*pi))>=50, 0<tL<=25, "
            "q=2*t*L^2; the live frequency layer has q>=1 and "
            f"0<tL<={C_STAR}+o(1)"
        ),
        "local_lift": (
            "On each critical radius-1/L disk choose one prescribed-N "
            "analytic lift; an adjacent-lift difference is retained in r_[1]"
        ),
        "base_components": (
            "s=(1-i*x)/2, s_*=s+t*alpha(s)/2, "
            "phi=M_t(s)/|M_t(s)|, "
            "e_n=phi*exp(t*log(n)^2/4-s_*log(n))"
        ),
        "first_correction": (
            "alpha_n=alpha(s)-log(n), "
            "d_n=1/(6s)+alpha'(s)*(t/4+t^2*alpha_n^2/8), "
            "f_n=e_n*(1+d_n)"
        ),
        "correction_derivative": (
            "d_(n,x)=(-i/2)*[-1/(6s^2)+alpha''(s)"
            "*(t/4+t^2*alpha_n^2/8)"
            "+(t^2/4)*alpha_n*alpha'(s)^2]"
        ),
        "endpoint_component": (
            "T_0=x/2+pi*t/8, a=sqrt(T_0/(2*pi)), "
            "p=1-2*(a-N), F=C_0, "
            "kappa_N=(-1)^N*exp(t*pi^2/64)*M_0(iT_0)"
            "*U*exp(pi*i/8)/|M_t(s)|, "
            "g_0=-kappa_N*(F(p)+F'''(p)/(12*pi^2*a))"
        ),
        "corrected_complex_main": (
            "E_[1]=g_0+sum_(n=1)^N f_n, J_[1]=2*Re(E_[1])"
        ),
        "repartition": (
            "J_[1]=J_[0]+DeltaJ, r_[0]=DeltaJ+r_[1], "
            "so J_[0]+r_[0]=J_[1]+r_[1] and the value/derivative "
            "contact systems are identical"
        ),
        "rate_free_partition": (
            "For components z_j in {g_0,f_1,...,f_N}, put "
            "c_j=Re(z_j), d_j=Re(z_(j,x)); P={c_j>0}, "
            "N_-={c_j<0}, Z={c_j=0}; "
            "M_+=sum_P c_j, M_-=sum_N_(-c_j), "
            "h_+=sum_P d_j/M_+, h_-="
            "sum_N_(-d_j)/M_-, D_0=sum_Z d_j, with an empty mean set to 0"
        ),
        "near_crossing": (
            "X_[1]=Re(E_[1])=M_+-M_-, "
            "U_[1]=Re(E_[1],x)=M*(h_+-h_-)"
            "+(X_[1]/2)*(h_++h_-)+D_0, M=(M_++M_-)/2"
        ),
        "contact_equation": (
            "At H=H_x=0, X_[1]=-r_[1]/2 and "
            "U_[1]=-r_[1],x/2; hence "
            "M*(h_+-h_-)+D_0=-r_[1],x/2"
            "+(r_[1]/4)*(h_++h_-)"
        ),
        "certified_contact_box": (
            "|X_[1]|<50000*exp(-5L/4), "
            "|U_[1]|<100000*L*exp(-5L/4)"
        ),
        "signed_band_target": (
            "|X_[1]|<=50000*exp(-5L/4) implies "
            "|U_[1]|>100000*L*exp(-5L/4)"
        ),
        "partition_target": (
            "Whenever |X_[1]|<=50000*exp(-5L/4), prove "
            "|M*(h_+-h_-)+D_0|>"
            "25000*exp(-5L/4)*(4L+|h_++h_-|)"
        ),
        "saddle_center": (
            "lambda_a=phi'/phi-s_*'*log(a)=u_a+i*v_a; "
            "R_a=-s_*'*sum_(n=1)^N log(n/a)f_n"
            "+sum_(n=1)^N e_n*d_(n,x)+(g_(0,x)-lambda_a*g_0); "
            "E_[1],x=lambda_a*E_[1]+R_a"
        ),
        "centered_contact": (
            "At contact, Re(R_a)-v_a*Im(E_[1])="
            "-(r_[1],x-u_a*r_[1])/2, so "
            "|Re(R_a)-v_a*Im(E_[1])|<"
            "(100000L+50000|u_a|)*exp(-5L/4)"
        ),
        "outer_handoff": (
            f"For every fixed epsilon>0, tL>={C_STAR}+epsilon "
            "is already closed for sufficiently large L by the "
            "oscillatory-zeta theorem"
        ),
        "multiplicity_guard": (
            "H_tau(y)=y^2-2tau has no positive-time contact but its "
            "root slopes are 2*sqrt(2tau)->0; a fixed positive slope "
            "floor uniform to tau=0 would impose endpoint simplicity"
        ),
        "three_open_pieces": (
            "Prove the signed band only on q>=1 in the live asymptotic "
            "frequency layer; use a multiplicity-compatible parabolic/Hermite "
            "or boundary-degree theorem on q<1; terminate the bounded-L shoulder"
        ),
    }


def build_rows(exact: dict[str, str]) -> list[GateRow]:
    return [
        GateRow(
            id="np15f1scr_01_remainder",
            role="certified_input",
            readiness="available_asymptotic",
            claim="The first-order residual has explicit global C1 bounds.",
            formula=(
                "|r_[1]|<100000*exp(-5L/4), "
                "|r_[1],x|<200000*L*exp(-5L/4)"
            ),
            proof_boundary="Only on the stated L>=50, 0<tL<=25 layer.",
        ),
        GateRow(
            id="np15f1scr_02_local_lift",
            role="exact_scope",
            readiness="available_exact",
            claim="Every critical disk admits one fixed-N analytic main.",
            formula=exact["region"] + "; " + exact["local_lift"],
            proof_boundary=(
                "The adjacent-lift error is already included in the certified "
                "r_[1] constants."
            ),
        ),
        GateRow(
            id="np15f1scr_03_finite_components",
            role="exact_definition",
            readiness="available_exact",
            claim="The retained finite blocks are explicit complex components.",
            formula=exact["base_components"] + "; " + exact["first_correction"],
            proof_boundary="A definition using the checked first correction.",
        ),
        GateRow(
            id="np15f1scr_04_correction_derivative",
            role="exact_identity",
            readiness="available_exact",
            claim="The derivative of every signed finite correction is explicit.",
            formula=exact["correction_derivative"],
            proof_boundary="Symbolically differentiated from d_n.",
        ),
        GateRow(
            id="np15f1scr_05_endpoint",
            role="exact_definition",
            readiness="available_exact",
            claim="C0 and C1 form one endpoint pseudo-component.",
            formula=exact["endpoint_component"],
            proof_boundary="On a chosen prescribed-N analytic lift.",
        ),
        GateRow(
            id="np15f1scr_06_complex_main",
            role="exact_identity",
            readiness="available_exact",
            claim="The full first-order real main is one complex component sum.",
            formula=exact["corrected_complex_main"],
            proof_boundary="Exact on the real axis of the selected local lift.",
        ),
        GateRow(
            id="np15f1scr_07_repartition",
            role="exact_identity",
            readiness="available_exact",
            claim="Absorbing DeltaJ preserves the exact contact system.",
            formula=exact["repartition"],
            proof_boundary="No asymptotic estimate is used.",
        ),
        GateRow(
            id="np15f1scr_08_partition",
            role="exact_definition",
            readiness="available_exact",
            claim="A rate-free sign partition remains valid at component zeros.",
            formula=exact["rate_free_partition"],
            proof_boundary="No division by a zero-real-part component.",
        ),
        GateRow(
            id="np15f1scr_09_near_crossing",
            role="exact_identity",
            readiness="available_exact",
            claim="The corrected slope has an exact weighted-mean decomposition.",
            formula=exact["near_crossing"],
            proof_boundary="Pure finite-sum algebra.",
        ),
        GateRow(
            id="np15f1scr_10_contact",
            role="exact_identity",
            readiness="available_exact",
            claim="A true Xi contact forces a signed corrected-main equation.",
            formula=exact["contact_equation"],
            proof_boundary="Uses H=A(J_[1]+r_[1]) with A>0.",
        ),
        GateRow(
            id="np15f1scr_11_box",
            role="proved_composition",
            readiness="available_asymptotic",
            claim="The global residual constants give a sharp first-order box.",
            formula=exact["certified_contact_box"],
            proof_boundary="Necessary at contact, not an exclusion theorem.",
        ),
        GateRow(
            id="np15f1scr_12_frequency_target",
            role="open_theorem_target",
            readiness="open",
            claim="A strict crossing-band theorem would close the frequency layer.",
            formula=(
                exact["signed_band_target"]
                + "; "
                + exact["partition_target"]
            ),
            proof_boundary=(
                "Open Xi-specific arithmetic inequality; required only on q>=1, "
                "not uniformly down to t=0."
            ),
        ),
        GateRow(
            id="np15f1scr_13_saddle_center",
            role="exact_identity",
            readiness="available_exact",
            claim="Saddle centering isolates one logarithmic moment and endpoint defect.",
            formula=exact["saddle_center"] + "; " + exact["centered_contact"],
            proof_boundary=(
                "Exact observable reduction; no lower bound for it is proved."
            ),
        ),
        GateRow(
            id="np15f1scr_14_outer_handoff",
            role="proved_composition",
            readiness="available_asymptotic",
            claim="The existing oscillatory theorem removes the outer fixed rays.",
            formula=exact["outer_handoff"],
            proof_boundary=(
                "Existential in L for each epsilon and does not cover the endpoint "
                "c_* or the inner layer."
            ),
        ),
        GateRow(
            id="np15f1scr_15_multiplicity_guard",
            role="route_guard",
            readiness="available_exact",
            claim="The ultra-small layer cannot inherit a uniform slope floor.",
            formula=(
                exact["multiplicity_guard"]
                + "; "
                + exact["three_open_pieces"]
            ),
            proof_boundary=(
                "The quadratic model is a nonpromotion guard, not an Xi theorem; "
                "all three listed pieces remain open."
            ),
        ),
    ]


def build_artifact() -> dict:
    verify_symbolic_identities()
    exact = exact_identities()
    rows = build_rows(exact)
    return {
        "kind": STEM,
        "date": "2026-07-26",
        "status": (
            "exact first-order signed-contact reduction with certified "
            "frequency-layer box; not contact exclusion, Lambda<=0, or RH"
        ),
        "proof_boundary": (
            "This artifact proves the corrected component representation, "
            "old/new repartition equivalence, rate-free near-crossing identity, "
            "explicit contact box, and saddle-centered observable. It does not "
            "prove the strict Xi crossing-band inequality, the q<1 "
            "multiplicity-compatible theorem, the bounded-L shoulder, contact "
            "exclusion, Lambda<=0, RH, PF-infinity, or a Clay-prize result."
        ),
        "builder_sha256": file_hash(Path(__file__)),
        "source_sha256": source_hashes(),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "certified_inputs_or_compositions": 3,
            "exact_definitions_or_identities": 10,
            "route_guards": 1,
            "open_frequency_targets": 1,
            "eta_0": ETA_0,
            "eta_1": ETA_1,
        },
        "sources": [
            str(path.relative_to(REPO_ROOT)).replace("\\", "/")
            for path in SOURCES.values()
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return "\n".join(
        [
            "# Newman Polymath-15 Critical First-Order Signed-Contact Reduction",
            "",
            "Date: 2026-07-26",
            "",
            "Status: exact first-order contact reduction with certified error",
            "constants. This is not a proof of contact exclusion, `Lambda <= 0`,",
            "or RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Corrected Complex Main",
            "",
            "On each radius-`1/L` proof disk, select one prescribed-`N`",
            "analytic lift. The adjacent-lift discrepancy is already part of",
            "the certified residual.",
            "",
            "```text",
            exact["base_components"],
            exact["first_correction"],
            exact["correction_derivative"],
            exact["endpoint_component"],
            exact["corrected_complex_main"],
            "```",
            "",
            "The previous peel can now be absorbed without changing the zero",
            "system:",
            "",
            "```text",
            exact["repartition"],
            "```",
            "",
            "## Rate-Free Contact Identity",
            "",
            "Use real component values and derivatives, so isolated component",
            "zeros do not create singular logarithmic rates.",
            "",
            "```text",
            exact["rate_free_partition"],
            exact["near_crossing"],
            exact["contact_equation"],
            "```",
            "",
            "## Certified Contact Box",
            "",
            "The global first-order theorem supplies",
            "",
            "```text",
            exact["certified_contact_box"],
            "```",
            "",
            "Therefore the exact open frequency-layer theorem is",
            "",
            "```text",
            exact["signed_band_target"],
            "```",
            "",
            "A stronger partition form is",
            "",
            "```text",
            exact["partition_target"],
            "```",
            "",
            "Neither inequality is proved here.",
            "",
            "## Saddle-Centered Observable",
            "",
            "Centering at the moving Riemann-Siegel saddle gives",
            "",
            "```text",
            exact["saddle_center"],
            exact["centered_contact"],
            "```",
            "",
            "This is the proof-facing arithmetic object: a centered logarithmic",
            "moment, the explicit `d_n` derivative, and one endpoint defect.",
            "",
            "## Domain Split",
            "",
            "```text",
            exact["region"],
            exact["outer_handoff"],
            exact["three_open_pieces"],
            "```",
            "",
            "The restriction to `q>=1` is logical, not cosmetic:",
            "",
            "```text",
            exact["multiplicity_guard"],
            "```",
            "",
            "A slope floor uniform to `t=0` would demand endpoint simplicity,",
            "which RH does not require. The ultra-small parabolic layer therefore",
            "needs a multiplicity-compatible Hermite or boundary-degree argument.",
            "",
            "## Proof Boundary",
            "",
            artifact["proof_boundary"],
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman critical first-order signed-contact reduction: "
        "15 rows, eta_0=100000, eta_1=200000, "
        "1 open frequency target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
