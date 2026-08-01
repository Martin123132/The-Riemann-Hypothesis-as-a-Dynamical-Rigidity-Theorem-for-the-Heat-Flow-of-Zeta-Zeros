#!/usr/bin/env python3
"""Build the ray-aligned energy endpoint-margin and conditioning audit."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
import re

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_ray_aligned_energy_endpoint_margin_gate"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "energy_current": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_parabolic_frequency_energy_current_gate.json"
    ),
    "ray_schedule": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_ray_aligned_parabolic_frequency_reduction.json"
    ),
    "compact": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_theta_compact_transversality_interval_certificate.json"
    ),
    "theta_operator": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_theta_curvature_probability_operator_gate.json"
    ),
    "dominant_global": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_dominant_saddle_global_ray_certificate.json"
    ),
    "dominant_cell": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_dominant_saddle_cell_interior_certificate.json"
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


BALL_RE = re.compile(
    r"^\[(?P<center>[0-9.eE+-]+) \+/- (?P<radius>[0-9.eE+-]+)\]$"
)


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def ball_lower(value: str) -> Decimal:
    match = BALL_RE.match(value)
    if not match:
        raise ValueError(f"unsupported ball string: {value}")
    return Decimal(match.group("center")) - Decimal(match.group("radius"))


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCES.items()
    }


def audit_sources(payloads: dict[str, dict]) -> dict:
    energy = payloads["energy_current"]
    if energy.get("summary", {}).get("pointwise_bridge_rows") != 2:
        raise RuntimeError("energy pointwise bridge drifted")
    if energy.get("summary", {}).get("pointwise_contact_exclusion"):
        raise RuntimeError("energy gate was unexpectedly promoted")

    schedule = payloads["ray_schedule"]
    if schedule.get("summary", {}).get("new_strip_open_antecedents") != 0:
        raise RuntimeError("ray-aligned new-strip closure drifted")
    if schedule.get("summary", {}).get("old_collar_open_antecedents") != 1:
        raise RuntimeError("ray-aligned old-collar handoff drifted")

    compact = payloads["compact"]
    certificate = compact.get("certificate", {})
    records = certificate.get("records", [])
    edge_records = [
        row for row in records if row.get("x_high") == "38"
    ]
    if len(edge_records) != 10:
        raise RuntimeError("x=38 compact edge record count drifted")
    if any(row.get("x_low") != "757/20" for row in edge_records):
        raise RuntimeError("x=38 compact edge cell drifted")
    if any(row.get("branch") != "derivative" for row in edge_records):
        raise RuntimeError("x=38 compact edge is no longer derivative-only")
    if any(
        ball_lower(row["certified_ratio_lower"]) <= Decimal(7) / Decimal(4)
        for row in edge_records
    ):
        raise RuntimeError("x=38 derivative ratio no longer exceeds 7/4")
    intervals = sorted(
        (Fraction(row["t_low"]), Fraction(row["t_high"]))
        for row in edge_records
    )
    if intervals[0][0] != 0 or intervals[-1][1] != Fraction(1, 5):
        raise RuntimeError("x=38 compact time coverage drifted")
    if any(right != next_left for (_, right), (next_left, _) in zip(
        intervals, intervals[1:]
    )):
        raise RuntimeError("x=38 compact time intervals no longer tile")

    operator = payloads["theta_operator"].get("exact", {}).get(
        "dominant_block_disjunction", {}
    )
    if "delta_0=1/2800" not in operator.get("tail_constants", ""):
        raise RuntimeError("theta tail delta_0 drifted")
    if "4*x^3" not in operator.get("tail_derivative_bound", ""):
        raise RuntimeError("theta B_1 cubic term drifted")
    if "|J_t'-J_(1,t)'|<B_1" not in operator.get("bounds", ""):
        raise RuntimeError("theta derivative tail bound drifted")

    global_ray = payloads["dominant_global"]
    theorem = global_ray.get("exact", {}).get("theorem", "")
    if "C_t[H_t/A_t](x)>3/40*L^2" not in theorem:
        raise RuntimeError("dominant global curvature margin drifted")
    if global_ray.get("parameters", {}).get("tL_min") != "25/1":
        raise RuntimeError("dominant global tL threshold drifted")

    cell = payloads["dominant_cell"]
    cell_exact = cell.get("exact", {})
    if cell_exact.get("normalized_jets") != (
        "eps0<=delta, eps1<=0.35*L*delta, "
        "eps2<=0.13*L^2*delta, delta=1/8000"
    ):
        raise RuntimeError("dominant normalized jet caps drifted")
    if cell_exact.get("main_jets") != (
        "|P|<=2.256, |P'|<=0.59L, |P''|<=0.154L^2"
    ):
        raise RuntimeError("dominant main jet caps drifted")
    cell_rows = {row["id"]: row for row in cell.get("rows", [])}
    normalized_row = cell_rows.get("np15dscic_05_normalized_jets", {})
    boundary = normalized_row.get("proof_boundary", "")
    if "|(log A)'|<=0.27L" not in boundary or "|V|<=10^-21" not in boundary:
        raise RuntimeError("dominant normalizer derivative caps drifted")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCES.items()
        },
        "compact_edge_records": len(edge_records),
        "compact_edge_time_intervals": [
            f"[{left},{right}]" for left, right in intervals
        ],
        "compact_edge_branches": sorted(
            {row["branch"] for row in edge_records}
        ),
        "compact_edge_ratio_lower_min": min(
            row["certified_ratio_lower"] for row in edge_records
        ),
        "inherited_boundary": (
            "The endpoint sources are quantitative, but neither source "
            "supplies an upper bound for the energy bulk through the full "
            "old descendant collar."
        ),
    }


def build_exact() -> dict:
    x_left = 38
    active_stage = 25
    ratio_floor = Fraction(7, 4)
    delta_0 = Fraction(1, 2800)
    j_prime_margin = (
        (ratio_floor - 1) * 4 * x_left**3 * delta_0
    )
    if j_prime_margin != Fraction(20577, 350):
        raise RuntimeError("left J-prime margin failed")

    z_0 = Fraction(282, 125) + Fraction(1, 8000)
    z_1 = Fraction(59, 100) + Fraction(7, 160000)
    z_2 = Fraction(77, 500) + Fraction(13, 800000)
    ell_cap = Fraction(27, 100)
    ell_prime_over_l2 = Fraction(1, 10**23)
    h_xx_cap = (
        z_2
        + 2 * ell_cap * z_1
        + (ell_cap**2 + ell_prime_over_l2) * z_0
    )
    if not h_xx_cap < Fraction(13, 20):
        raise RuntimeError("normalized H_xx cap failed")

    active_l_min = 101 + active_stage
    q_min = 50 * active_l_min
    if q_min != 6300:
        raise RuntimeError("active-stage q floor failed")
    beta_squared_floor = Fraction(q_min, q_min + 1)
    y_floor = Fraction(99, 1000)
    curvature_floor = Fraction(3, 40)
    test_upper = (
        y_floor**2 / beta_squared_floor
        + Fraction(13, 20) * y_floor
    )
    curvature_gap = curvature_floor - test_upper
    if curvature_gap != Fraction(593211, 700000000):
        raise RuntimeError("right normalized margin gap failed")
    if curvature_gap <= 0:
        raise RuntimeError("right normalized margin is not strict")

    asymptotic_prefactor = (
        25 * j_prime_margin**2
        / (512 * x_left**8)
    )
    if asymptotic_prefactor != Fraction(9, 231853260800):
        raise RuntimeError("left floor asymptotic prefactor failed")

    return {
        "active_schedule": {
            "range": "j>=25",
            "reason": (
                "t_j=25/(100+j)<=1/5 exactly when j>=25; earlier "
                "finite stages remain separate base/transition bookkeeping "
                "and are not used in this cofinal endpoint estimate"
            ),
            "log_radius": "L_j=101+j>=126",
            "times": (
                "t_(j+1)=25/L_j and t_j=25/(L_j-1)"
            ),
            "radius": "R_j=4*pi*exp(L_j)",
            "width": "W_j=R_j-38=4*pi*exp(L_j)-38",
        },
        "left_source": {
            "endpoint": "x=38, 0<t<=1/5",
            "edge_records": 10,
            "edge_cell": "757/20<=x<=38",
            "branches": "all ten edge records use the derivative branch",
            "ratio": "|J_(1,t)'|/B_1>7/4",
            "tail": "|J_t'-J_(1,t)'|<B_1",
            "b1_floor": "B_1(t,38)>=4*38^3/2800",
            "full_margin": (
                "|J_t'(38)|>3*38^3/2800="
                + fraction_text(j_prime_margin)
            ),
        },
        "left_energy_margin": {
            "characteristic_derivative": (
                "J_t'(x)=64x^3H_t(x)+16x^4H_t'(x)"
            ),
            "cauchy": (
                "|J_t'(38)|<=sqrt((64*38^3)^2"
                "+(16*38^4/s_pf(t,38))^2)*sqrt(E_pf(t,38))"
            ),
            "pointwise": (
                "E_pf(t,38)>M_38^2/"
                "[4096*38^6+256*38^8/s_pf(t,38)^2], "
                "M_38=20577/350"
            ),
            "stage_uniform": (
                "E_pf(t,38)>="
                "M_38^2/[4096*38^6+256*38^8"
                "*(log(38/(4*pi))^2+L_j/50)]"
            ),
            "notation": "call the square root of this lower bound m_(L,j)",
        },
        "right_source": {
            "endpoint": "x=R_j, j>=25",
            "ray": "tL_j>=25 and L_j>=126",
            "curvature": (
                "C=(H_x^2-H H_xx)/A_t^2>3L_j^2/40"
            ),
            "normalization": (
                "Z=H/A_t, g=H_x/A_t, h_2=H_xx/A_t, "
                "C=g^2-Zh_2"
            ),
            "jet_caps": {
                "|Z|": f"{fraction_text(z_0)}",
                "|Z'|/L": f"{fraction_text(z_1)}",
                "|Z''|/L^2": f"{fraction_text(z_2)}",
                "|ell|/L": "27/100",
                "|ell'|/L^2": "<=1/10^23",
                "|h_2|/L^2": (
                    f"<{fraction_text(Fraction(13, 20))}"
                ),
                "exact_h2_cap_before_rounding": fraction_text(h_xx_cap),
            },
        },
        "right_energy_margin": {
            "scaled_norm": (
                "y=sqrt(E_pf(t,R_j))/A_t(R_j)"
                "=sqrt(Z^2+s_pf^2 g^2)"
            ),
            "beta_floor": (
                "s_pf^2 L_j^2=q/(1+q)>="
                f"{fraction_text(beta_squared_floor)}, q>=6300"
            ),
            "curvature_upper": (
                "C/L_j^2<=y^2/(s_pf^2L_j^2)+(13/20)y"
            ),
            "test": (
                "at y=99/1000 the upper side is "
                f"{fraction_text(test_upper)}"
                f"=3/40-{fraction_text(curvature_gap)}"
            ),
            "conclusion": (
                "sqrt(E_pf(t,R_j))>(99/1000)A_t(R_j)"
            ),
        },
        "certified_contact_floor": {
            "two_endpoint": (
                "a contact in [38,R_j] implies integral_38^R_j B_pf dx"
                ">=[sqrt(E_pf(t,38))+sqrt(E_pf(t,R_j))]^2/W_j"
            ),
            "source_lower": (
                "integral_38^R_j B_pf dx"
                ">=[m_(L,j)+(99/1000)A_t(R_j)]^2/W_j"
            ),
            "left_only": (
                "integral_38^R_j B_pf dx>=m_(L,j)^2/W_j"
            ),
            "strict_budget_target": (
                "U_j(t)<2[m_(L,j)+(99/1000)A_t(R_j)]^2/W_j"
            ),
            "upper_budget_requirement": (
                "D_j(t)<=U_j(t), with U_j derived independently from "
                "Xi source bounds and D_j=d/dt integral E_pf"
                "+2[J_pf(R_j)-J_pf(38)]"
            ),
        },
        "conditioning": {
            "left_floor_asymptotic": (
                "m_(L,j)^2/W_j~"
                "[9/(231853260800*pi)]*exp(-L_j)/L_j"
            ),
            "prefactor_without_pi": fraction_text(asymptotic_prefactor),
            "meaning": (
                "The currently certified raw floor loses one exponential "
                "radius factor. This is a conditioning statement about the "
                "available sufficient bound, not a claim that the actual "
                "endpoint energy has this exact size."
            ),
        },
        "dominant_coverage": {
            "threshold": "x_dom(t)=4*pi*exp(25/t)",
            "bottom": "x_dom(t_(j+1))=R_j",
            "top": "x_dom(t_j)=R_(j-1)",
            "covered_tail": "[x_dom(t),R_j]",
            "uncovered_for_bulk": "[38,x_dom(t))",
            "consequence": (
                "The dominant-ray C2 bounds control only the right tail; "
                "at the collar bottom they meet the interval only at R_j."
            ),
        },
        "numerator_gap": (
            "The compact certificate stops at x=38 and the dominant-ray "
            "certificate starts at x_dom(t). Neither supplies a full-collar "
            "upper bound U_j for D_j. Pointwise noncontact in other bands "
            "does not itself bound the integrated H_x and H_xx bulk."
        ),
        "nonpromotion": {
            "tautology": (
                "Substituting D_j=2 integral B_pf is an identity, not the "
                "strict source-level upper estimate required by the test."
            ),
            "normalizer": (
                "The right margin is relative to A_t(R_j); A_t cannot be "
                "dropped or replaced by one without a proved comparison."
            ),
            "finite_stage": (
                "A finite evaluation of rho_j cannot prove the all-stage "
                "bulk budget."
            ),
            "conditioning": (
                "Failure of the coarse certified floor to beat a candidate "
                "upper bound would reject that bound, not disprove the "
                "underlying contact-free theorem."
            ),
        },
        "route_decision": (
            "Retain both endpoint margins, but do not launch a whole-collar "
            "bulk computation until a source-level interior majorant exists. "
            "The preferred refinement is a localized energy criterion on "
            "bounded logarithmic or phase cells with independently certified "
            "endpoint margins; otherwise return to the contact-normal C1 "
            "arithmetic target."
        ),
        "open_handoff": (
            "Derive an exact localization of the energy/current identity on "
            "bounded cells and inventory which cell endpoints already have "
            "quantitative first-jet margins from compact, oscillatory-zeta, "
            "dominant-saddle, or phase-cell sources. Determine whether the "
            "cell boundary fluxes telescope without requiring the same open "
            "interior theorem."
        ),
    }


def diagnostic_rows() -> list[dict[str, str | int]]:
    x = 38.0
    l_0 = math.log(x / (4 * math.pi))
    margin = 3 * x**3 / 2800
    rows: list[dict[str, str | int]] = []
    for stage in (25, 100, 1000, 10_000):
        ell = 101 + stage
        time = 25 / ell
        inverse_s_squared = l_0**2 + 1 / (2 * time)
        left_energy = margin**2 / (
            4096 * x**6 + 256 * x**8 * inverse_s_squared
        )
        log10_width = (math.log(4 * math.pi) + ell) / math.log(10)
        log10_floor = math.log10(left_energy) - log10_width
        rows.append(
            {
                "stage": stage,
                "L_j": ell,
                "bottom_time": f"{time:.17g}",
                "certified_left_energy_lower": f"{left_energy:.17g}",
                "certified_left_norm_lower": (
                    f"{math.sqrt(left_energy):.17g}"
                ),
                "left_only_floor_log10": f"{log10_floor:.17g}",
                "status": "conditioning_diagnostic_not_proof_input",
            }
        )
    return rows


def build_payload() -> dict:
    payloads = load_sources()
    source_audit = audit_sources(payloads)
    exact = build_exact()
    rows = [
        GateRow(
            id="raem_01_active_schedule",
            role="exact_geometry",
            readiness="proved_exact",
            claim="The compact quantitative endpoint source is available on every nontrivial ray-aligned stage j>=25.",
            formula=f"{exact['active_schedule']['range']}; {exact['active_schedule']['times']}",
            proof_boundary="No new assertion is made about the finite initial-stage bookkeeping.",
        ),
        GateRow(
            id="raem_02_compact_edge_inventory",
            role="source_audit",
            readiness="certified",
            claim="The x=38 boundary is covered by ten contiguous derivative-branch interval boxes.",
            formula=f"{exact['left_source']['edge_records']} records; {exact['left_source']['branches']}",
            proof_boundary="Uses only the stored rigorous compact certificate.",
        ),
        GateRow(
            id="raem_03_compact_ratio",
            role="exact_composition",
            readiness="available_exact",
            claim="Every compact edge box clears the derivative tail budget by a ratio greater than seven quarters.",
            formula=f"{exact['left_source']['ratio']}; {exact['left_source']['tail']}",
            proof_boundary="The ratio is an Arb lower enclosure, not a sampled midpoint.",
        ),
        GateRow(
            id="raem_04_left_full_margin",
            role="exact_inequality",
            readiness="proved_exact",
            claim="The full characteristic derivative has one explicit uniform raw margin at x=38.",
            formula=exact["left_source"]["full_margin"],
            proof_boundary="Uses the positive cubic term in B_1 and delta_0=1/2800.",
        ),
        GateRow(
            id="raem_05_left_energy_conversion",
            role="exact_inequality",
            readiness="proved_exact",
            claim="Weighted Cauchy-Schwarz converts the characteristic derivative margin into an E_pf endpoint margin.",
            formula=exact["left_energy_margin"]["pointwise"],
            proof_boundary="Valid for t>0; no division by the unknown jet occurs.",
        ),
        GateRow(
            id="raem_06_left_stage_margin",
            role="exact_inequality",
            readiness="proved_exact",
            claim="The bottom time gives one explicit margin uniform over each descendant collar.",
            formula=exact["left_energy_margin"]["stage_uniform"],
            proof_boundary="Uses monotonicity of s_pf^2 in t.",
        ),
        GateRow(
            id="raem_07_right_curvature",
            role="source_audit",
            readiness="certified",
            claim="The exact dominant-ray theorem supplies a strict normalizer-relative curvature margin at R_j.",
            formula=exact["right_source"]["curvature"],
            proof_boundary="Applies because tL_j>=25, including cutoff transitions.",
        ),
        GateRow(
            id="raem_08_right_hxx_cap",
            role="exact_composition",
            readiness="proved_exact",
            claim="The dominant main and remainder C2 caps give a uniform normalized second-jet upper bound.",
            formula="|H_xx|/(A_t L_j^2)<13/20",
            proof_boundary="Uses the certified normalizer log-derivative bounds.",
        ),
        GateRow(
            id="raem_09_right_beta_floor",
            role="exact_geometry",
            readiness="proved_exact",
            claim="The parabolic-frequency scale is uniformly frequency-conditioned at every active right endpoint.",
            formula=exact["right_energy_margin"]["beta_floor"],
            proof_boundary="Uses j>=25 and tL_j>=25.",
        ),
        GateRow(
            id="raem_10_right_energy_margin",
            role="exact_inequality",
            readiness="proved_exact",
            claim="Curvature positivity and the C2 cap force a quantitative normalized first-jet margin.",
            formula=f"{exact['right_energy_margin']['test']}; {exact['right_energy_margin']['conclusion']}",
            proof_boundary="The raw margin retains the explicit normalizer A_t(R_j).",
        ),
        GateRow(
            id="raem_11_two_endpoint_floor",
            role="exact_composition",
            readiness="available_exact",
            claim="The two source margins instantiate the energy gate's exact contact floor.",
            formula=exact["certified_contact_floor"]["source_lower"],
            proof_boundary="Conditional consequence if a contact exists; not a contact proof by itself.",
        ),
        GateRow(
            id="raem_12_strict_budget_target",
            role="conditional_criterion",
            readiness="not_ready_to_apply",
            claim="One independently derived full-collar bulk upper bound below the certified floor would exclude contact.",
            formula=exact["certified_contact_floor"]["strict_budget_target"],
            proof_boundary=exact["certified_contact_floor"]["upper_budget_requirement"],
        ),
        GateRow(
            id="raem_13_floor_conditioning",
            role="exact_asymptotic",
            readiness="proved_exact",
            claim="The left-only certified whole-collar floor loses one exponential radius factor.",
            formula=exact["conditioning"]["left_floor_asymptotic"],
            proof_boundary=exact["conditioning"]["meaning"],
        ),
        GateRow(
            id="raem_14_dominant_coverage",
            role="exact_geometry",
            readiness="proved_exact",
            claim="The dominant-ray C2 theorem covers only a moving right tail of each old collar.",
            formula=f"{exact['dominant_coverage']['covered_tail']}; {exact['dominant_coverage']['uncovered_for_bulk']}",
            proof_boundary=exact["dominant_coverage"]["consequence"],
        ),
        GateRow(
            id="raem_15_interior_numerator_gap",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The audited sources do not provide the full-collar numerator upper bound.",
            formula=exact["numerator_gap"],
            proof_boundary="This is the remaining Xi source obligation.",
        ),
        GateRow(
            id="raem_16_tautology_guard",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="The exact energy identity cannot be recycled as its own strict upper estimate.",
            formula=exact["nonpromotion"]["tautology"],
            proof_boundary="Blocks a circular bulk-budget proof.",
        ),
        GateRow(
            id="raem_17_normalizer_and_finite_guards",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="Neither the right normalizer nor all-stage quantification may be discarded.",
            formula=f"{exact['nonpromotion']['normalizer']} {exact['nonpromotion']['finite_stage']}",
            proof_boundary=exact["nonpromotion"]["conditioning"],
        ),
        GateRow(
            id="raem_18_localization_handoff",
            role="route_decision",
            readiness="available_exact",
            claim="The next energy test should be localized to bounded cells with quantitative endpoint ownership.",
            formula=exact["open_handoff"],
            proof_boundary="No localized telescoping or cellwise Xi bulk theorem is yet proved.",
        ),
    ]
    return {
        "kind": STEM,
        "date": "2026-07-28",
        "status": (
            "exact ray-aligned endpoint margins and whole-collar "
            "bulk-budget conditioning audit"
        ),
        "proof_boundary": (
            "This artifact proves the compact left derivative margin, its "
            "parabolic-frequency energy conversion, the dominant-ray "
            "normalizer-relative right energy margin, the instantiated "
            "two-endpoint contact floor, its cofinal conditioning, and the "
            "moving dominant-tail coverage geometry. It does not prove a "
            "full-collar Xi bulk upper bound, the strict no-contact budget, "
            "a localized telescoping theorem, Q209, the cofinal descendant "
            "theorem, Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion."
        ),
        "source_audit": source_audit,
        "exact": exact,
        "diagnostics": diagnostic_rows(),
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "compact_edge_records": source_audit["compact_edge_records"],
            "left_raw_margins": 1,
            "right_normalized_margins": 1,
            "exact_contact_floors": 1,
            "open_full_collar_bulk_budgets": 1,
            "nonpromotion_guards": sum(
                row.role == "nonpromotion_gate" for row in rows
            ),
            "pointwise_contact_exclusions": 0,
            "cofinal_descendant_theorem": False,
        },
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    left = exact["left_source"]
    left_energy = exact["left_energy_margin"]
    right = exact["right_source"]
    right_energy = exact["right_energy_margin"]
    floor = exact["certified_contact_floor"]
    diagnostics = payload["diagnostics"]
    diagnostic_lines = [
        (
            f"| {row['stage']} | {row['L_j']} | {row['bottom_time']} | "
            f"{row['certified_left_norm_lower']} | "
            f"{row['left_only_floor_log10']} |"
        )
        for row in diagnostics
    ]
    return "\n".join(
        [
            "# Ray-Aligned Energy Endpoint-Margin Gate",
            "",
            "Date: 2026-07-28",
            "",
            "Status: exact endpoint margins and whole-collar conditioning",
            "audit. This is not a proof of the descendant theorem,",
            "Lambda<=0, or RH.",
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
            "validated ray-aligned energy endpoint-margin gate: 18 rows, "
            "10 compact edge records, 1 left raw margin, "
            "1 right normalized margin, 1 exact contact floor, "
            "1 open full-collar bulk budget, 2 nonpromotion guards, "
            "0 pointwise contact exclusions",
            "```",
            "",
            "## Active Stages",
            "",
            "The quantitative compact endpoint source applies to",
            "",
            "```text",
            exact["active_schedule"]["range"],
            exact["active_schedule"]["log_radius"],
            exact["active_schedule"]["times"],
            exact["active_schedule"]["width"],
            "```",
            "",
            exact["active_schedule"]["reason"],
            "",
            "Pi provenance remains unchanged: `pi` is the ordinary circle",
            "constant inherited from `L=log(x/(4*pi))`.",
            "",
            "## Left Endpoint",
            "",
            "At `x=38`, all ten stored compact edge boxes use the",
            "derivative branch and tile `0<=t<=1/5`. They prove",
            "",
            "```text",
            left["ratio"],
            left["tail"],
            left["b1_floor"],
            left["full_margin"],
            "```",
            "",
            "The last line is a margin for the full exact heat-flow",
            "characteristic derivative, not only its first theta block.",
            "",
            "Because",
            "",
            "```text",
            left_energy["characteristic_derivative"],
            left_energy["cauchy"],
            "```",
            "",
            "we obtain",
            "",
            "```text",
            left_energy["pointwise"],
            left_energy["stage_uniform"],
            "```",
            "",
            "No unknown jet norm was used as a divisor.",
            "",
            "## Right Endpoint",
            "",
            "At `x=R_j`, the complete dominant-ray theorem gives",
            "",
            "```text",
            right["ray"],
            right["curvature"],
            right["normalization"],
            "```",
            "",
            "The main-sum, remainder, and normalizer derivative bounds",
            "compose to",
            "",
            "```text",
            "|H_xx|/(A_t L_j^2)<13/20",
            right_energy["beta_floor"],
            right_energy["curvature_upper"],
            right_energy["test"],
            right_energy["conclusion"],
            "```",
            "",
            "The factor `A_t(R_j)` is retained. It is not normalized away.",
            "",
            "## Exact Contact Floor",
            "",
            "The prior energy gate now gives the fully sourced implication",
            "",
            "```text",
            floor["source_lower"],
            "```",
            "",
            "A sufficient strict theorem would be",
            "",
            "```text",
            floor["upper_budget_requirement"],
            floor["strict_budget_target"],
            "```",
            "",
            "That upper budget is not currently proved.",
            "",
            "## Conditioning",
            "",
            "The left-only certified floor has the exact asymptotic",
            "",
            "```text",
            exact["conditioning"]["left_floor_asymptotic"],
            "```",
            "",
            exact["conditioning"]["meaning"],
            "",
            "| stage j | L_j | bottom t | certified left norm | "
            "log10(left-only floor) |",
            "|---:|---:|---:|---:|---:|",
            *diagnostic_lines,
            "",
            "The table is a deterministic conditioning diagnostic, not a",
            "finite-stage proof or counterexample.",
            "",
            "## Interior Numerator Gap",
            "",
            "For fixed `t`, the dominant theorem begins at",
            "",
            "```text",
            exact["dominant_coverage"]["threshold"],
            exact["dominant_coverage"]["bottom"],
            exact["dominant_coverage"]["top"],
            exact["dominant_coverage"]["covered_tail"],
            exact["dominant_coverage"]["uncovered_for_bulk"],
            "```",
            "",
            exact["numerator_gap"],
            "",
            "## Route Decision",
            "",
            exact["nonpromotion"]["tautology"],
            "",
            exact["nonpromotion"]["normalizer"],
            "",
            exact["nonpromotion"]["finite_stage"],
            "",
            exact["route_decision"],
            "",
            "Open handoff:",
            "",
            exact["open_handoff"],
            "",
            "## Boundary",
            "",
            payload["proof_boundary"],
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "built ray-aligned energy endpoint-margin gate: 18 rows, "
        "10 compact edge records, 1 left raw margin, "
        "1 right normalized margin, 1 exact contact floor, "
        "1 open full-collar bulk budget, 2 nonpromotion guards, "
        "0 pointwise contact exclusions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
