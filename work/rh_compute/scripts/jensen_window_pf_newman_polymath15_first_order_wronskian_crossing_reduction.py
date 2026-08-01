#!/usr/bin/env python3
"""Build the first-order Wronskian and oriented-crossing reduction."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_wronskian_crossing_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "global_remainder": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_"
        "first_order_global_remainder_certificate.json"
    ),
    "signed_contact": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_"
        "first_order_signed_contact_reduction.json"
    ),
    "oriented_successor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_oriented_successor_winding_reduction.json"
    ),
    "old_wronskian": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_wronskian_phase_reduction.json"
    ),
    "component_wronskian": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_component_wronskian_gate.json"
    ),
    "slope_gap": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "crossing_slope_gap_reduction.json"
    ),
}

ETA_0 = 100_000
ETA_1 = 200_000
HALF_SQUARED_BUDGET = 12_500_000_000


@dataclass(frozen=True)
class ReductionRow:
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
        "global_remainder": (
            '"eta_0": 100000',
            '"eta_1": 200000',
            "exp(-5L/4)",
        ),
        "signed_contact": (
            "E_[1]=g_0+sum_(n=1)^N f_n",
            "E_[1],x=lambda_a*E_[1]+R_a",
        ),
        "oriented_successor": (
            "positive imaginary ray",
            "wind(proxy_j)<1",
        ),
        "old_wronskian": (
            "W_E",
            "phase",
        ),
        "component_wronskian": (
            "inertia (1,1",
            "ordered",
        ),
        "slope_gap": (
            "h_+-h_-",
            "common half-plane",
        ),
    }
    for key, required in markers.items():
        text = json.dumps(payloads[key], sort_keys=True)
        for marker in required:
            if marker not in text:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {key: str(payload.get("kind", "")) for key, payload in payloads.items()}


def symbolic_audit() -> dict:
    x_part, y_part, u_part, v_part = sp.symbols(
        "X Y U V", real=True
    )
    complex_main = x_part + sp.I * y_part
    complex_derivative = u_part + sp.I * v_part
    wronskian = sp.expand_complex(
        sp.im(complex_derivative * sp.conjugate(complex_main))
    )
    if sp.simplify(wronskian - (v_part * x_part - u_part * y_part)) != 0:
        raise RuntimeError("Wronskian coordinate identity failed")
    crossing_wronskian = sp.simplify(wronskian.subs(x_part, 0))
    if crossing_wronskian != -u_part * y_part:
        raise RuntimeError("crossing Wronskian identity failed")

    scale = sp.symbols("L", positive=True)
    proxy_norm = x_part**2 + (u_part / scale) ** 2
    crossing_norm = sp.simplify(proxy_norm.subs(x_part, 0))
    if crossing_norm != u_part**2 / scale**2:
        raise RuntimeError("crossing norm identity failed")

    threshold = sp.sqrt(HALF_SQUARED_BUDGET)
    if sp.simplify(threshold - 50_000 * sp.sqrt(5)) != 0:
        raise RuntimeError("crossing threshold simplification failed")
    return {
        "coordinates": (
            "E_[1]=X+iY, E_[1],x=U+iV, "
            "W_[1]=Im(E_[1],x*conj(E_[1]))"
        ),
        "wronskian": "W_[1]=V*X-U*Y",
        "proxy": "J_[1]=2X, J_[1],x=2U",
        "crossing": "At X=0, W_[1]=-U*Y.",
        "proxy_crossing_norm": (
            "At X=0, X^2+(U/L)^2=U^2/L^2."
        ),
        "threshold": "sqrt(12500000000)=50000*sqrt(5)",
        "sympy_wronskian": str(wronskian),
        "sympy_crossing": str(crossing_wronskian),
    }


def build_exact() -> dict:
    symbolic = symbolic_audit()
    return {
        "symbolic": symbolic,
        "contact": {
            "equations": (
                "At a full contact, X=-r_[1]/2 and U=-r_[1],x/2."
            ),
            "identity": "W_[1]=V*X-U*Y",
            "bound": (
                "2*|W_[1]|<exp(-5L/4)*"
                "(100000*|E_[1],x|+200000*L*|E_[1]|)"
            ),
            "sufficient_disjunction": (
                "For every point in the live layer, prove either "
                "|X|>50000*exp(-5L/4) or "
                "2*|W_[1]|>exp(-5L/4)*"
                "(100000*|E_[1],x|+200000*L*|E_[1]|)."
            ),
        },
        "boundary_crossing": {
            "margin_equivalence": (
                "At X=0 and E_[1]!=0, "
                "X^2+(U/L)^2>12500000000*exp(-5L/2) iff "
                "|W_[1]|>50000*sqrt(5)*L*|E_[1]|*exp(-5L/4)."
            ),
            "phase_speed": (
                "Where E_[1]!=0, theta_[1],x=W_[1]/|E_[1]|^2."
            ),
            "upward_sign": (
                "At X=0, E_[1]=iY!=0, an upward zero U>0 is "
                "equivalent to W_[1]*Y<0."
            ),
            "exceptional_zero": (
                "If E_[1]=0 then X=Y=W_[1]=0; a proxy crossing with U>0 "
                "is not classified by W_[1]*Y and must be counted separately."
            ),
            "count_decomposition": (
                "N_up=N_(X=0,E_[1]!=0,W_[1]*Y<0)"
                "+N_(E_[1]=0,U>0), with half-open endpoint conventions."
            ),
        },
        "successor_link": {
            "crossing_formula": (
                "Insert the decomposed N_up on both horizontal edges into "
                "kappa_j=N_up(t_j;R_j)-N_up(t_(j+1);R_j)+I_i(V_j), "
                "then add the finite shoulder/join intersection cells."
            ),
            "one_sided_target": (
                "It is enough to prove the resulting oriented count is <1; "
                "positive contact degree then forces kappa_j=0."
            ),
            "distinction": (
                "The absolute Wronskian margin proves boundary nonvanishing. "
                "The signed W_[1]*Y counts and exceptional complex-main zeros "
                "are additional data needed for the one-sided winding bound."
            ),
        },
        "route_guards": {
            "phase_sign": (
                "The aggregate phase derivative is not known to have one sign; "
                "stored corrected crossings exhibit both signs below the "
                "asymptotic L>=50 layer."
            ),
            "component_form": (
                "The instantaneous component-rate Hermitian form has inertia "
                "(1,1,m-2) when rates differ, so it is not positive semidefinite."
            ),
            "ordered_speeds": (
                "Positive amplitudes and ordered same-sign component speeds "
                "admit an exact double-crossing countermodel."
            ),
            "finite_diagnostics": (
                "Finite sign diagnostics shape the route but neither prove nor "
                "disprove the L>=50 heat-coupled theorem."
            ),
        },
        "open_target": {
            "frequency_layer": (
                "On q=2tL^2>=1 in the live 0<tL<c_*+o(1) layer, prove the "
                "first-order Wronskian disjunction and a signed crossing/phase "
                "budget including E_[1]=0 exceptional crossings."
            ),
            "parabolic_layer": (
                "Treat q<1 by a multiplicity-compatible Hermite or degree chart; "
                "do not demand a fixed phase-speed floor at t=0."
            ),
            "shoulders": (
                "Close bounded L, L_epsilon, vertical connector, and chart-join "
                "phase cells separately."
            ),
        },
        "proof_handoff": (
            "Use the Wronskian magnitude to certify first-jet separation at "
            "real-part crossings, and use the signs of W_[1]*Im(E_[1]) plus "
            "the explicit E_[1]=0 count to upper-bound the successor winding. "
            "Do not replace this signed arithmetic theorem by phase "
            "monotonicity or component-rate ordering."
        ),
    }


def build_rows(exact: dict) -> list[ReductionRow]:
    return [
        ReductionRow(
            "nfowc_01_coordinates",
            "exact_identity",
            "available_exact",
            "The first-order complex main has an exact real-coordinate Wronskian.",
            exact["symbolic"]["coordinates"] + "; " + exact["symbolic"]["wronskian"],
            "Algebraic identity only.",
            exact["symbolic"],
        ),
        ReductionRow(
            "nfowc_02_contact_wronskian",
            "exact_contact_reduction",
            "available_exact",
            "A full contact forces an explicit exponentially small first-order Wronskian.",
            exact["contact"]["bound"],
            "Uses the certified eta_0 and eta_1 remainder bounds.",
            exact["contact"],
        ),
        ReductionRow(
            "nfowc_03_contact_disjunction",
            "exact_sufficient_target",
            "available_exact",
            "A value-or-Wronskian disjunction excludes every full contact in the live frequency layer.",
            exact["contact"]["sufficient_disjunction"],
            "The Xi-specific strict disjunction remains unproved.",
        ),
        ReductionRow(
            "nfowc_04_crossing_margin",
            "exact_equivalence",
            "available_exact",
            "At a nonzero complex-main crossing, the first-order vector margin is a normalized Wronskian margin.",
            exact["boundary_crossing"]["margin_equivalence"],
            "The E_[1]=0 case is excluded by the strict Wronskian inequality.",
        ),
        ReductionRow(
            "nfowc_05_upward_sign",
            "exact_crossing_reduction",
            "available_exact",
            "The sign of W_[1]*Im(E_[1]) identifies upward real-part crossings when E_[1] is nonzero.",
            exact["boundary_crossing"]["upward_sign"],
            "A signed count, not a global phase-monotonicity claim.",
        ),
        ReductionRow(
            "nfowc_06_exceptional_zero",
            "nonpromotion_guard",
            "guard_validated",
            "Zeros of the complete complex main form a separate crossing class.",
            exact["boundary_crossing"]["exceptional_zero"],
            "They cannot be divided out or hidden in the Wronskian sign.",
        ),
        ReductionRow(
            "nfowc_07_count_decomposition",
            "exact_crossing_reduction",
            "available_exact",
            "Every upward crossing splits into a Wronskian-sign crossing or an exceptional complex-main zero.",
            exact["boundary_crossing"]["count_decomposition"],
            "Assumes the standard half-open endpoint convention.",
        ),
        ReductionRow(
            "nfowc_08_successor_link",
            "exact_composition",
            "available_exact",
            "The decomposed crossing count plugs directly into the oriented successor recurrence.",
            exact["successor_link"]["crossing_formula"],
            "Finite connector and shoulder cells remain explicit.",
            exact["successor_link"],
        ),
        ReductionRow(
            "nfowc_09_phase_sign_guard",
            "nonpromotion_guard",
            "guard_validated",
            "A one-signed aggregate phase derivative is not an available shortcut.",
            exact["route_guards"]["phase_sign"],
            "The stored examples are below L=50 and are route guards only.",
            exact["route_guards"],
        ),
        ReductionRow(
            "nfowc_10_frequency_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "Prove the first-order Wronskian magnitude and signed crossing budget in the q>=1 layer.",
            exact["open_target"]["frequency_layer"],
            "This is the live Xi-specific arithmetic target.",
        ),
        ReductionRow(
            "nfowc_11_parabolic_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "Supply the multiplicity-compatible q<1 crossing theorem.",
            exact["open_target"]["parabolic_layer"],
            "No endpoint-simplicity floor is permitted.",
        ),
        ReductionRow(
            "nfowc_12_shoulder_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "Terminate all finite phase shoulders and connector cells.",
            exact["open_target"]["shoulders"],
            "Finite but not currently certified.",
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    return {
        "kind": STEM,
        "date": "2026-07-26",
        "status": (
            "exact first-order Wronskian, upward-crossing, and exceptional-zero "
            "reduction; three Xi crossing obligations remain open"
        ),
        "proof_boundary": (
            "This artifact proves the first-order contact Wronskian bound, "
            "crossing-margin equivalence, upward-sign identity, exceptional-zero "
            "split, and successor-count composition. It does not prove the "
            "q>=1 Wronskian theorem, q<1 multiplicity-compatible theorem, "
            "finite shoulders, one-sided successor winding, contact exclusion, "
            "Lambda<=0, RH, PF-infinity, or a Clay-prize conclusion."
        ),
        "constants": {
            "eta_0": ETA_0,
            "eta_1": ETA_1,
            "half_squared_budget": HALF_SQUARED_BUDGET,
            "crossing_threshold": "50000*sqrt(5)",
        },
        "sources": {
            key: str(path.relative_to(REPO_ROOT)).replace("\\", "/")
            for key, path in SOURCES.items()
        },
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "exact": exact,
        "rows": [asdict(row) for row in build_rows(exact)],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return "\n".join(
        [
            "# Newman First-Order Wronskian Crossing Reduction",
            "",
            "Date: 2026-07-26",
            "",
            "Status: exact Wronskian and oriented-crossing reduction. The",
            "Xi arithmetic theorem remains open; this is not a proof of RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## First-Order Wronskian",
            "",
            "```text",
            exact["symbolic"]["coordinates"],
            exact["symbolic"]["wronskian"],
            exact["symbolic"]["proxy"],
            exact["contact"]["equations"],
            exact["contact"]["bound"],
            exact["contact"]["sufficient_disjunction"],
            "```",
            "",
            "## Boundary Crossings",
            "",
            "```text",
            exact["symbolic"]["crossing"],
            exact["boundary_crossing"]["margin_equivalence"],
            exact["boundary_crossing"]["phase_speed"],
            exact["boundary_crossing"]["upward_sign"],
            exact["boundary_crossing"]["exceptional_zero"],
            exact["boundary_crossing"]["count_decomposition"],
            "```",
            "",
            "A Wronskian magnitude bound supplies separation, but the winding",
            "also needs signed crossings and the explicit `E_[1]=0` class.",
            "",
            "## Successor Composition",
            "",
            "```text",
            exact["successor_link"]["crossing_formula"],
            exact["successor_link"]["one_sided_target"],
            exact["successor_link"]["distinction"],
            "```",
            "",
            "## Route Guards",
            "",
            "```text",
            exact["route_guards"]["phase_sign"],
            exact["route_guards"]["component_form"],
            exact["route_guards"]["ordered_speeds"],
            exact["route_guards"]["finite_diagnostics"],
            "```",
            "",
            "## Remaining Xi Input",
            "",
            "```text",
            exact["open_target"]["frequency_layer"],
            exact["open_target"]["parabolic_layer"],
            exact["open_target"]["shoulders"],
            "```",
            "",
            "## Live Handoff",
            "",
            "```text",
            exact["proof_handoff"],
            "```",
            "",
            artifact["proof_boundary"],
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman first-order Wronskian crossing reduction: "
        "12 rows, eta_0=100000, eta_1=200000, "
        "2 exact crossing classes, 3 open Xi obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
