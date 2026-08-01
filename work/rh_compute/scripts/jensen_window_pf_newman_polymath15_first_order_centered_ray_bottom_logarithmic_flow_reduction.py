#!/usr/bin/env python3
"""Build the centered ray-bottom logarithmic-flow reduction."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_centered_ray_bottom_logarithmic_flow_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "rectangular_boundary": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_rectangular_boundary_degree_reduction.json"
    ),
    "normalized_phase_flux": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_normalized_prefix_phase_flux_reduction.json"
    ),
    "abel_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_carrier_kernel_abel_prefix_reduction.json"
    ),
    "abel_shear": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_abel_scalar_shear_flux_reduction.json"
    ),
    "target_reconciliation": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_target_reconciliation_gate.json"
    ),
    "connector_cap": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_dominant_ray_connector_phase_cap.json"
    ),
}


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCES.items()
    }


def audit_sources(payloads: dict[str, dict]) -> dict:
    rectangular = payloads["rectangular_boundary"].get("exact", {})
    if rectangular.get("ray_bottom", {}).get("hard_interval") != (
        "max(B_epsilon,sqrt((100+j)/50))<=L"
        "<=min(L_j,(c_*+epsilon)*(100+j)/25)"
    ):
        raise RuntimeError("hard-bottom source interval drifted")

    normalized = payloads["normalized_phase_flux"].get("exact", {})
    if "gamma_n=q_(n,x)/q_n" not in normalized.get(
        "coefficient_current", ""
    ):
        raise RuntimeError("normalized carrier current drifted")
    if "Z_(A,x)=" not in normalized.get("normalized_flux_derivative", ""):
        raise RuntimeError("normalized endpoint derivative drifted")

    prefix = payloads["abel_prefix"].get("exact", {})
    if "mathsf_A=Re(W_A)=mathcal_C_N+c*u_N*mathsf_X" not in (
        prefix.get("contact_scalar", "")
    ):
        raise RuntimeError("Abel contact scalar drifted")

    shear = payloads["abel_shear"].get("exact", {})
    if "partial_s arg(Psi_ell)" not in shear.get("general_flux", ""):
        raise RuntimeError("Abel proxy flux drifted")
    if "2.03e-14" not in shear.get("orientation_ratio", ""):
        raise RuntimeError("crossing orientation ratio drifted")

    target = payloads["target_reconciliation"].get("exact", {})
    if target.get("abel_target", {}).get("pointwise") != (
        "|mathsf_X|<=delta_L => "
        "|mathcal_C_N|>A_L+epsilon_term"
    ):
        raise RuntimeError("pointwise Abel target drifted")

    connector = payloads["connector_cap"].get("exact", {})
    if "H_j<3*pi/2" not in connector.get(
        "reduced_horizontal_target", ""
    ):
        raise RuntimeError("joined horizontal target drifted")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCES.items()
        },
        "source_kinds": {
            key: payload.get("kind") for key, payload in payloads.items()
        },
    }


def build_rows() -> list[GateRow]:
    return [
        GateRow(
            "rlf_01_ray_schedule",
            "one-dimensional domain",
            "proved",
            "Restrict a generic horizontal ray to the cofinal schedule.",
            "t=t_j=25/(100+j), x(L)=4*pi*exp(L), "
            "a_j(L)^2=exp(L)+t_j/16.",
            "The result is local to one horizontal ray.",
        ),
        GateRow(
            "rlf_02_hard_interval",
            "outer arithmetic domain",
            "proved",
            "The unresolved low-c q>=1 ray piece is one closed interval.",
            "L_- =max(B_epsilon,sqrt((100+j)/50)), "
            "L_+=min(101+j,(c_*+epsilon)*(100+j)/25).",
            "If L_->L_+, the hard interval is empty.",
        ),
        GateRow(
            "rlf_03_cutoff_cells",
            "cutoff geometry",
            "proved",
            "The hard interval has an exact canonical fixed-N partition.",
            "lambda_(j,N)=log(N^2-t_j/16), "
            "I_(j,N)=[lambda_(j,N),lambda_(j,N+1)).",
            "Each cell is intersected with [L_-,L_+].",
        ),
        GateRow(
            "rlf_04_equality_ownership",
            "cutoff convention",
            "proved",
            "At a cutoff equality the canonical chart is the new N chart.",
            "L=lambda_(j,N) => a_j(L)=N and floor(a_j(L))=N.",
            "The adjacent old chart is retained only through its overlap homotopy.",
        ),
        GateRow(
            "rlf_05_log_chain_rule",
            "exact reparametrization",
            "proved",
            "Fixed-ray logarithmic differentiation is positive x differentiation.",
            "D_j:=d/dL|_(t=t_j,N fixed)=x*partial_x, x=4*pi*exp(L)>0.",
            "The reparametrization preserves crossing orientation and winding.",
        ),
        GateRow(
            "rlf_06_terminal_distance",
            "exact coefficient",
            "proved",
            "The terminal logarithmic distance has an explicit L rate.",
            "u_(N,L)=x/(4*T_0)=exp(L)/"
            "[2*(exp(L)+t_j/16)], T_0=2*pi*a^2.",
            "N is held fixed inside the cell.",
        ),
        GateRow(
            "rlf_07_carrier_rates",
            "exact carrier flow",
            "proved",
            "Every normalized carrier has an exact logarithmic rate.",
            "widehat_gamma_n=(D_j q_n)/q_n=x*gamma_n="
            "x[-s_*'*log(n)+epsilon_n-epsilon_1].",
            "The terminal and correction factors are not dropped.",
        ),
        GateRow(
            "rlf_08_inward_current",
            "certified q>=1 sign",
            "proved",
            "The fixed-x inward current becomes a uniform L-current floor.",
            "Re(widehat_gamma_n)<-log(2)/(16L^2)<0 for n>=2.",
            "This is an individual-carrier statement, not an aggregate sign.",
        ),
        GateRow(
            "rlf_09_ordered_rotation",
            "certified relative phase",
            "proved",
            "Relative angular ordering is preserved under the positive scaling.",
            "Im(widehat_gamma_m-widehat_gamma_n)>="
            "(x/2)log(m/n)-16892/x>0 for n<m.",
            "Ordered carrier rotation does not sign a prefix current.",
        ),
        GateRow(
            "rlf_10_prefix_flow",
            "exact O(N) derivative",
            "proved",
            "Every normalized prefix differentiates by the scaled carrier rates.",
            "D_j G_k=sum_(n=1)^k widehat_gamma_n*q_n.",
            "The formula includes k=N and is valid only inside a fixed-N cell.",
        ),
        GateRow(
            "rlf_11_endpoint_value_flow",
            "endpoint-complete derivative",
            "proved",
            "The normalized value retains the recurrent endpoint.",
            "D_j Z_0=x*r_(0,x)+sum_(n=1)^N widehat_gamma_n*q_n.",
            "No bulk-only replacement is made.",
        ),
        GateRow(
            "rlf_12_endpoint_slope_flow",
            "endpoint-complete derivative",
            "proved",
            "The normalized slope has an exact O(N) logarithmic derivative.",
            "D_j Z_A=(x*s_*''*u_N+s_*'*u_(N,L))*Z_0"
            "+s_*'*u_N*D_j Z_0+D_j B_N.",
            "D_j B_N retains r_A, r_0, and all their derivatives.",
        ),
        GateRow(
            "rlf_13_projector_flow",
            "absolute phase anchor",
            "proved",
            "The moving unit projector contributes its exact frame current.",
            "eta_L=i*Omega_eta*eta; D_j mathsf_X=Pi_eta(D_j Z_0)"
            "-Omega_eta*Im(eta*Z_0), with the analogous mathsf_A formula.",
            "The absolute projector cannot be frozen.",
        ),
        GateRow(
            "rlf_14_abel_scalar_flow",
            "division-free derivative",
            "proved",
            "The terminal shear gives an endpoint-complete Abel derivative.",
            "alpha=c*u_N; D_j mathcal_C_N=D_j mathsf_A"
            "-(D_j alpha)*mathsf_X-alpha*D_j mathsf_X; equivalently "
            "differentiate Pi_eta(B_N)-b*u_N*mathsf_Y directly.",
            "D_j alpha=x*Re(s_*'')*u_N+c*u_(N,L).",
        ),
        GateRow(
            "rlf_15_ray_proxy",
            "topological coordinate",
            "conditional",
            "Under the pointwise gap use the physical L-scaled Abel proxy.",
            "Psi_(j,N)(L)=mathsf_X+i*mathcal_C_N/L.",
            "The pointwise Xi Abel gap remains open.",
        ),
        GateRow(
            "rlf_16_exact_phase_flux",
            "one-dimensional flux",
            "proved",
            "The ray proxy has one exact signed phase numerator.",
            "D_j arg(Psi)=[L*mathsf_X*D_j mathcal_C_N"
            "-mathsf_X*mathcal_C_N-L*mathcal_C_N*D_j mathsf_X]"
            "/[L^2*mathsf_X^2+mathcal_C_N^2].",
            "No absolute-value sum over carriers is licensed.",
        ),
        GateRow(
            "rlf_17_orientation_defect",
            "source-specific decomposition",
            "proved",
            "Separate the negative crossing core from the endpoint-complete defect.",
            "E_j=D_j mathsf_X-x*mathcal_C_N; numerator="
            "-L*x*mathcal_C_N^2+mathsf_X*(L*D_j mathcal_C_N"
            "-mathcal_C_N)-L*mathcal_C_N*E_j.",
            "The middle term has no proved Xi sign.",
        ),
        GateRow(
            "rlf_18_scaled_orientation",
            "crossing orientation",
            "conditional",
            "The certified x-orientation estimate is unchanged in L scale.",
            "|E_j|/(x*A_L)<2.03e-14 on |mathsf_X|<=delta_L.",
            "It becomes effective only after the open Abel gap is supplied.",
        ),
        GateRow(
            "rlf_19_ray_intersections",
            "signed crossings",
            "conditional",
            "The pointwise gap makes every real crossing simple and oriented.",
            "mathsf_X=0 => D_j arg(Psi)=-L*D_j mathsf_X/mathcal_C_N<0.",
            "For increasing L, mathcal_C_N>0 contributes -1; reversal gives +1.",
        ),
        GateRow(
            "rlf_20_cutoff_join",
            "chart topology",
            "conditional",
            "Join Abel proxies through the physical-slope chart at each cutoff.",
            "Psi_N ~ (mathsf_X+i*mathsf_A/L)_N "
            "~ (mathsf_X+i*mathsf_A/L)_(N+1) ~ Psi_(N+1).",
            "No complex adjacent-chart bound for mathcal_C_N is asserted.",
        ),
        GateRow(
            "rlf_21_local_flux_target",
            "sharpened sufficient condition",
            "open",
            "One scalar inequality would force clockwise flux inside the contact band.",
            "|mathsf_X|*|L*D_j mathcal_C_N-mathcal_C_N|"
            "+L*|mathcal_C_N|*|E_j|"
            "<L*x*mathcal_C_N^2.",
            "This condition is sufficient, not known and not necessary.",
        ),
        GateRow(
            "rlf_22_recrossing_guard",
            "nonpromotion guard",
            "proved",
            "A strong gap and correctly oriented crossings do not bound their number.",
            "The exact backward-heat and two-frequency guards have arbitrarily "
            "many same-orientation crossings.",
            "Individual inward motion and ordered rotation cannot prove one turn.",
        ),
        GateRow(
            "rlf_23_joined_ledger",
            "successor composition",
            "open",
            "The hard rays must be compared only after every open-path seam is joined.",
            "2*pi*kappa_j=H_j+C_j, |C_j|<pi/2; prove H_j<3*pi/2.",
            "Endpoint tracks, q<1 arcs, shoulders, and cutoff joins are counted once.",
        ),
        GateRow(
            "rlf_24_handoff",
            "route decision",
            "open",
            "Attack the endpoint-complete signed ray defect before broader searches.",
            "First bound mathsf_X*(L*D_j mathcal_C_N-mathcal_C_N) "
            "on each hard cell; then compare the two adjacent rays in H_j.",
            "No Abel gap, winding bound, contact exclusion, or RH conclusion follows.",
        ),
    ]


def build_payload() -> dict:
    payloads = load_sources()
    source_audit = audit_sources(payloads)
    rows = build_rows()
    proof_boundary = (
        "This artifact proves the exact hard-ray cutoff partition, logarithmic "
        "chain rule, scaled carrier and prefix currents, endpoint-complete "
        "normalized value and slope derivatives, Abel-scalar derivative, ray "
        "phase-flux identity, defect decomposition, and orientation transfer. "
        "It does not prove the Xi Abel gap, a sign for the aggregate defect, "
        "the joined H_j<3*pi/2 bound, q<1 or finite-shoulder closure, contact "
        "exclusion, Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion."
    )
    return {
        "kind": STEM,
        "date": "2026-07-29",
        "status": (
            "exact_ray_bottom_logarithmic_flow_with_open_"
            "abel_gap_and_joined_phase_target"
        ),
        "proof_boundary": proof_boundary,
        "source_audit": source_audit,
        "exact": {
            "ray": {
                "schedule": (
                    "t_j=25/(100+j), x(L)=4*pi*exp(L), "
                    "a_j(L)^2=exp(L)+t_j/16"
                ),
                "q": "q_j(L)=2*t_j*L^2=50*L^2/(100+j)",
                "hard_interval": (
                    "L_- =max(B_epsilon,sqrt((100+j)/50)), "
                    "L_+=min(101+j,(c_*+epsilon)*(100+j)/25)"
                ),
            },
            "cutoff_partition": {
                "cutoff": "lambda_(j,N)=log(N^2-t_j/16)",
                "cell": (
                    "I_(j,N)=[lambda_(j,N),lambda_(j,N+1)), "
                    "intersected with [L_-,L_+]"
                ),
                "equality": (
                    "At L=lambda_(j,N), a_j(L)=N and the canonical chart "
                    "has cutoff N"
                ),
            },
            "logarithmic_derivative": {
                "operator": (
                    "D_j=d/dL|_(t=t_j,N fixed)=x*partial_x, "
                    "x=4*pi*exp(L)>0"
                ),
                "terminal_distance": (
                    "u_(N,L)=x/(4*T_0)="
                    "exp(L)/[2*(exp(L)+t_j/16)]"
                ),
                "phase_current": "eta_L=i*Omega_eta*eta, Omega_eta=x*omega_eta",
            },
            "carrier_flow": {
                "rate": (
                    "widehat_gamma_n=x*gamma_n="
                    "x[-s_*'*log(n)+epsilon_n-epsilon_1]"
                ),
                "prefix": (
                    "D_j G_k=sum_(n=1)^k widehat_gamma_n*q_n"
                ),
                "inward": (
                    "Re(widehat_gamma_n)<-log(2)/(16L^2)<0 for n>=2"
                ),
                "ordered_rotation": (
                    "Im(widehat_gamma_m-widehat_gamma_n)>="
                    "(x/2)*log(m/n)-16892/x>0 for n<m"
                ),
            },
            "endpoint_flow": {
                "value": (
                    "D_j Z_0=x*r_(0,x)+"
                    "sum_(n=1)^N widehat_gamma_n*q_n"
                ),
                "v_tilde": (
                    "D_j V_tilde_N=x*r_(A,x)-x*s_*''*u_N*r_0"
                    "-s_*'*u_(N,L)*r_0-s_*'*u_N*x*r_(0,x)"
                ),
                "b": (
                    "D_j B_N=D_j V_tilde_N+x*s_*''*"
                    "sum_(k=1)^(N-1)h_k*G_k+s_*'*"
                    "sum_(k=1)^(N-1)h_k*D_j G_k"
                ),
                "slope": (
                    "D_j Z_A=(x*s_*''*u_N+s_*'*u_(N,L))*Z_0"
                    "+s_*'*u_N*D_j Z_0+D_j B_N"
                ),
            },
            "projected_flow": {
                "value": (
                    "D_j mathsf_X=Pi_eta(D_j Z_0)"
                    "-Omega_eta*Im(eta*Z_0)"
                ),
                "slope": (
                    "D_j mathsf_A=Pi_eta(D_j Z_A)"
                    "-Omega_eta*Im(eta*Z_A)"
                ),
                "abel": (
                    "alpha=c*u_N; D_j mathcal_C_N=D_j mathsf_A"
                    "-(D_j alpha)*mathsf_X-alpha*D_j mathsf_X"
                ),
                "alpha": (
                    "D_j alpha=x*Re(s_*'')*u_N+c*u_(N,L)"
                ),
                "direct_abel": (
                    "mathcal_C_N=Pi_eta(B_N)-b*u_N*mathsf_Y, "
                    "mathsf_Y=Im(eta*Z_0)"
                ),
                "direct_abel_derivative": (
                    "D_j mathcal_C_N=Pi_eta(D_j B_N)"
                    "-Omega_eta*Im(eta*B_N)"
                    "-[x*Im(s_*'')*u_N+b*u_(N,L)]*mathsf_Y"
                    "-b*u_N*[Im(eta*D_j Z_0)+Omega_eta*mathsf_X]"
                ),
            },
            "ray_proxy": {
                "definition": (
                    "Psi_(j,N)(L)=mathsf_X+i*mathcal_C_N/L"
                ),
                "flux": (
                    "D_j arg(Psi)=[L*mathsf_X*D_j mathcal_C_N"
                    "-mathsf_X*mathcal_C_N-L*mathcal_C_N*D_j mathsf_X]"
                    "/[L^2*mathsf_X^2+mathcal_C_N^2]"
                ),
                "orientation_defect": (
                    "E_j=D_j mathsf_X-x*mathcal_C_N"
                ),
                "defect_flux": (
                    "numerator=-L*x*mathcal_C_N^2"
                    "+mathsf_X*(L*D_j mathcal_C_N-mathcal_C_N)"
                    "-L*mathcal_C_N*E_j"
                ),
                "scaled_error": (
                    "|E_j|/(x*A_L)<2.03e-14 "
                    "on |mathsf_X|<=delta_L"
                ),
                "crossing": (
                    "At mathsf_X=0, D_j arg(Psi)="
                    "-L*D_j mathsf_X/mathcal_C_N"
                ),
            },
            "cutoff_join": {
                "path": (
                    "Psi_N ~ (mathsf_X+i*mathsf_A/L)_N "
                    "~ (mathsf_X+i*mathsf_A/L)_(N+1) ~ Psi_(N+1)"
                ),
                "scope": (
                    "Use the terminal-shear homotopy inside each chart and "
                    "the certified real first-jet homotopy between charts"
                ),
            },
            "conditional_targets": {
                "pointwise_gap": (
                    "|mathsf_X|<=delta_L => "
                    "|mathcal_C_N|>A_L+epsilon_term"
                ),
                "local_clockwise_sufficient": (
                    "|mathsf_X|*|L*D_j mathcal_C_N-mathcal_C_N|"
                    "+L*|mathcal_C_N|*|E_j|"
                    "<L*x*mathcal_C_N^2"
                ),
                "joined_horizontal": (
                    "After top, bottom, q<1, shoulder, cutoff, chart, and "
                    "connector-endpoint pieces are joined, prove H_j<3*pi/2"
                ),
            },
            "nonpromotion": {
                "aggregate_current": (
                    "Individual inward and ordered carrier currents do not "
                    "assign a sign to D_j mathcal_C_N or the proxy flux"
                ),
                "recrossing": (
                    "A pointwise crossing gap and correct local orientation "
                    "do not bound the number of same-orientation crossings"
                ),
            },
            "open_handoff": (
                "On every nonempty hard interval, work cell by cell with the "
                "exact endpoint-complete defect "
                "mathsf_X*(L*D_j mathcal_C_N-mathcal_C_N). Test whether Xi "
                "multiplicative completion and the recurrent endpoint give a "
                "signed or integrable upper bound before taking absolute "
                "values. Then compare the two oppositely oriented adjacent "
                "rays only inside the fully joined H_j ledger. Preserve q=1, "
                "W_0=0, cutoff equalities, and both adjacent charts."
            ),
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_cutoff_partitions": 1,
            "exact_logarithmic_flow_identities": 8,
            "exact_phase_flux_identities": 2,
            "orientation_transfers": 1,
            "cutoff_join_contracts": 1,
            "nonpromotion_guards": 2,
            "open_Xi_targets": 2,
            "proved_abel_gaps": 0,
            "proved_horizontal_phase_bounds": 0,
            "pointwise_contact_exclusions": 0,
        },
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    lines = [
        "# First-Order Ray-Bottom Logarithmic-Flow Reduction",
        "",
        "Date: 2026-07-29",
        "",
        "Status: exact fixed-cutoff L-flow and signed phase reduction with",
        "open Xi Abel and joined-ray targets. This is not a proof of contact",
        "exclusion, Lambda<=0, or RH.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"python work/rh_compute/scripts/{STEM}.py",
        f"python work/rh_compute/scripts/check_{STEM}.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        (
            "validated Newman centered ray-bottom logarithmic-flow reduction: "
            "24 rows, 1 exact cutoff partition, 8 logarithmic-flow identities, "
            "2 phase-flux identities, 1 orientation transfer, 1 cutoff-join "
            "contract, 2 nonpromotion guards, 2 open Xi targets, "
            "0 Abel gaps, 0 horizontal phase bounds"
        ),
        "```",
        "",
        "## Ray And Cutoff Cells",
        "",
        "On a generic cofinal ray,",
        "",
        "```text",
        exact["ray"]["schedule"],
        exact["ray"]["q"],
        exact["ray"]["hard_interval"],
        "```",
        "",
        "The canonical cutoff changes at",
        "",
        "```text",
        exact["cutoff_partition"]["cutoff"],
        exact["cutoff_partition"]["cell"],
        exact["cutoff_partition"]["equality"],
        "```",
        "",
        "Thus every equality boundary is owned explicitly. The old and new",
        "charts are compared by the certified adjacent homotopy rather than",
        "by identifying their complex Abel coordinates.",
        "",
        "## Exact Logarithmic Flow",
        "",
        "Inside one fixed-N cell,",
        "",
        "```text",
        exact["logarithmic_derivative"]["operator"],
        exact["logarithmic_derivative"]["terminal_distance"],
        exact["carrier_flow"]["rate"],
        exact["carrier_flow"]["prefix"],
        "```",
        "",
        "The positive factor `x` preserves the certified carrier signs:",
        "",
        "```text",
        exact["carrier_flow"]["inward"],
        exact["carrier_flow"]["ordered_rotation"],
        "```",
        "",
        "These remain individual-carrier facts. The pair kernels from the",
        "normalized-prefix gate still contain unrestricted sine and cosine",
        "terms, so no aggregate sign is inferred.",
        "",
        "## Endpoint-Complete Derivatives",
        "",
        "Write `V_tilde_N=r_A-s_*'*u_N*r_0` and",
        "`B_N=V_tilde_N+s_*' sum_k h_k G_k`. Then",
        "",
        "```text",
        exact["endpoint_flow"]["value"],
        exact["endpoint_flow"]["v_tilde"],
        exact["endpoint_flow"]["b"],
        exact["endpoint_flow"]["slope"],
        "```",
        "",
        "The moving unit phase is retained:",
        "",
        "```text",
        exact["logarithmic_derivative"]["phase_current"],
        exact["projected_flow"]["value"],
        exact["projected_flow"]["slope"],
        exact["projected_flow"]["abel"],
        exact["projected_flow"]["alpha"],
        exact["projected_flow"]["direct_abel"],
        exact["projected_flow"]["direct_abel_derivative"],
        "```",
        "",
        "This is an exact O(N) derivative with the recurrent endpoint and",
        "all endpoint derivatives still present.",
        "",
        "## Signed Ray Flux",
        "",
        "Under the still-open pointwise Abel gap, define",
        "",
        "```text",
        exact["ray_proxy"]["definition"],
        exact["ray_proxy"]["flux"],
        "```",
        "",
        "Put",
        "",
        "```text",
        exact["ray_proxy"]["orientation_defect"],
        exact["ray_proxy"]["defect_flux"],
        exact["ray_proxy"]["scaled_error"],
        "```",
        "",
        "The negative square is the source-specific crossing core. The",
        "remaining unsigned obstruction is the endpoint-complete term",
        "`mathsf_X*(L*D_j mathcal_C_N-mathcal_C_N)`. At a crossing,",
        "",
        "```text",
        exact["ray_proxy"]["crossing"],
        "```",
        "",
        "so the open Abel gap and the certified defect ratio imply the same",
        "orientation in `L` as in `x`. Increasing-L and decreasing-L edges",
        "carry opposite intersection signs.",
        "",
        "A sufficient, but unproved, pointwise clockwise condition is",
        "",
        "```text",
        exact["conditional_targets"]["local_clockwise_sufficient"],
        "```",
        "",
        "It is recorded as a testable scalar target, not as a theorem.",
        "",
        "## Cutoff And Successor Composition",
        "",
        "At each cutoff use the joined path",
        "",
        "```text",
        exact["cutoff_join"]["path"],
        "```",
        "",
        exact["cutoff_join"]["scope"] + ".",
        "",
        "A strong pointwise gap does not bound the number of recrossings, and",
        "ordered individual carriers do not sign the aggregate flux. The",
        "closed successor target therefore remains",
        "",
        "```text",
        exact["conditional_targets"]["joined_horizontal"],
        "```",
        "",
        "Open handoff:",
        "",
        exact["open_handoff"],
        "",
        "The `pi` in `x=4*pi*exp(L)`, the saddle, and the cutoff cells is",
        "the completed-zeta/Riemann-Siegel constant already recorded in the",
        "source chain. The `2*pi` phase normalization is the period of",
        "`exp(i*theta)`; no circle or fitted polygon is introduced here.",
        "",
        "## Boundary",
        "",
        payload["proof_boundary"],
        "",
    ]
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "built Newman centered ray-bottom logarithmic-flow reduction: "
        "24 rows, 1 exact cutoff partition, 8 logarithmic-flow identities, "
        "2 phase-flux identities, 1 orientation transfer, 1 cutoff-join "
        "contract, 2 nonpromotion guards, 2 open Xi targets, "
        "0 Abel gaps, 0 horizontal phase bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
