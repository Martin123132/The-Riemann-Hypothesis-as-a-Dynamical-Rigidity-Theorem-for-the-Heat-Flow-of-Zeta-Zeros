#!/usr/bin/env python3
"""Build the dominant-ray connector phase cap and reduced phase ledger."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "dominant_ray_connector_phase_cap"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "normalized_bridge": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_normalized_laguerre_bridge.json"
    ),
    "dominant_arithmetic": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "dominant_saddle_arithmetic_ray_certificate.json"
    ),
    "dominant_global": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "dominant_saddle_global_ray_certificate.json"
    ),
    "cofinal_boundary": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_cofinal_boundary_reduction.json"
    ),
    "oriented_successor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_oriented_successor_winding_reduction.json"
    ),
    "abel_shear": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "abel_scalar_shear_flux_reduction.json"
    ),
}


@dataclass(frozen=True)
class ConnectorRow:
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
        "normalized_bridge": (
            "M_t(s)=exp((t/4)*alpha(s)^2)*M_0(s)",
            "P_0=2*cos(beta)",
        ),
        "dominant_arithmetic": (
            "0.24L<=b=|beta'|<=0.26L",
            "Q0<=0.256, Q1<=0.06901L",
        ),
        "dominant_global": (
            '"eps0": "1/8000"',
            '"eps1_over_L": "7/160000"',
        ),
        "cofinal_boundary": (
            "t_j=25/(100+j), L_j=101+j",
            "ray_schedule",
        ),
        "oriented_successor": (
            "kappa_j=sum_(p in D_j) floor(m_p/2)",
            "positive_integer",
        ),
        "abel_shear": (
            "Delta_(partial D_j)arg(Psi_j)<2*pi",
            "complete boundary composition",
        ),
    }
    for key, required in markers.items():
        text = json.dumps(payloads[key], sort_keys=True)
        for marker in required:
            if marker not in text:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {
        key: str(payload.get("kind", ""))
        for key, payload in payloads.items()
    }


def rational_audit() -> dict[str, str]:
    value_error = Fraction(32, 125) + Fraction(1, 8000)
    slope_error = Fraction(6901, 100000) + Fraction(7, 160000)
    error_square = value_error**2 + slope_error**2
    error_cap = Fraction(4, 15)
    if not error_square < error_cap**2:
        raise RuntimeError("dominant connector vector error cap failed")

    leading_floor = Fraction(12, 25)
    relative_cap = error_cap / leading_floor
    if relative_cap != Fraction(5, 9):
        raise RuntimeError("dominant connector relative cap failed")

    sine_lower = Fraction(3, 5) - Fraction(3, 5) ** 3 / 6
    if sine_lower != Fraction(141, 250):
        raise RuntimeError("sine Taylor fraction drifted")
    if not sine_lower > relative_cap:
        raise RuntimeError("connector cone-angle cap failed")

    pi_upper = Fraction(22, 7)
    leading_phase_cap = Fraction(25, 384) * pi_upper + Fraction(1, 24)
    if leading_phase_cap != Fraction(331, 1344):
        raise RuntimeError("leading connector phase cap drifted")

    full_phase_cap = Fraction(6, 5) + leading_phase_cap
    if full_phase_cap != Fraction(9719, 6720):
        raise RuntimeError("full connector phase cap drifted")
    if not full_phase_cap < Fraction(3, 2):
        raise RuntimeError("strict connector phase cap failed")

    return {
        "value_error": str(value_error),
        "slope_error": str(slope_error),
        "error_square": str(error_square),
        "euclidean_error_cap": str(error_cap),
        "leading_norm_floor": str(leading_floor),
        "relative_error_cap": str(relative_cap),
        "sine_lower_at_3_over_5": str(sine_lower),
        "leading_phase_cap": str(leading_phase_cap),
        "full_phase_cap": str(full_phase_cap),
    }


def build_exact() -> dict:
    audit = rational_audit()
    return {
        "pi_provenance": (
            "The pi in x=4*pi*exp(L), the completed-zeta normalizer, and "
            "the Riemann-Siegel saddle is the same universal Euclidean "
            "constant already present in xi. The 2*pi phase turn is the "
            "period of exp(i*theta). The classical rational bounds "
            "3<pi<22/7 are used only to turn the source-derived angular "
            "bounds into strict rational inequalities; no circle or "
            "polygon is selected by the connector argument."
        ),
        "connector_geometry": (
            "For L=L_j=101+j and x=R_j=4*pi*exp(L), the positively "
            "oriented right edge of D_j traverses "
            "I_L=[25/L,25/(L-1)]. Its length is "
            "Delta t=25/[L*(L-1)], and t*L>=25 throughout."
        ),
        "normalized_vector": (
            "Put Z_t=H_t/A_t and V(t)=(Z_t(x),partial_x Z_t(x)/L). "
            "The positive normalizer and its derivative act on the full "
            "first jet by a positive determinant shear, so V may be used "
            "on the connector without changing the closed-boundary degree."
        ),
        "normalizer_phase": (
            "At fixed x, M_t(s)=exp[t*alpha(s)^2/4]*M_0(s), so "
            "beta_t=arg M_t(s) obeys "
            "partial_t beta=Im(alpha^2)/4="
            "Re(alpha)*Im(alpha)/2. For x=4*pi*exp(L), "
            "Re(alpha)<L/2 and |Im(alpha)|<pi/4; hence "
            "|Delta beta|<25*pi/[16*(L-1)]<=pi/64."
        ),
        "leading_ellipse": (
            "With b=-partial_x beta and lambda=b/L, the certified bounds "
            "give 6/25<=lambda<=13/50. The one-saddle normalized jet is "
            "V_0(t)=(2*cos(beta),2*lambda*sin(beta)), so "
            "|V_0|>=12/25."
        ),
        "uniform_vector_error": (
            "The finite arithmetic tail and exact-H collar give "
            "|V_1-V_(0,1)|<=32/125+1/8000=2049/8000 and "
            "|V_2-V_(0,2)|<=6901/100000+7/160000="
            "55243/800000. Their Euclidean norm is strictly below 4/15, "
            "so |V-V_0|/|V_0|<5/9 uniformly, including cutoff transitions."
        ),
        "uniform_cone": (
            "Because sin(3/5)>3/5-(3/5)^3/6=141/250>5/9, "
            "the relative angle delta(t)=arg V(t)-arg V_0(t) has one "
            "continuous branch with |delta(t)|<3/5 on the whole connector. "
            "Equivalently V/V_0 remains in a disk contained in one open "
            "half-plane, so no hidden relative turn is possible."
        ),
        "leading_phase_variation": (
            "For theta_0=arg(cos(beta)+i*lambda*sin(beta)), "
            "|d theta_0|<=(1/lambda)|d beta|"
            "+[1/(2*lambda)]|d lambda|. At fixed x, lambda is affine in t; "
            "its certified range has width 1/50. Therefore "
            "|Delta theta_0|<25*pi/384+1/24"
            "<331/1344, using pi<22/7."
        ),
        "connector_phase_cap": (
            "A continuous argument lift of V on I_L has total angular "
            "range below 6/5+331/1344=9719/6720<3/2<pi/2. In particular "
            "|Delta_(right D_j)arg V|<pi/2 for every j>=0. The exact "
            "dominant-ray connector therefore costs strictly less than "
            "one quarter turn."
        ),
        "axis_phase": (
            "On the symmetry axis H_t(0)>0 and partial_x H_t(0)=0. "
            "Its normalized first-jet path lies on the positive real ray "
            "and contributes zero phase."
        ),
        "composite_ledger": (
            "Insert the exact normalized connector into the boundary proxy "
            "through the already required nonvanishing homotopies, assigning "
            "the two endpoint tracks once to the complementary path. Write "
            "C_j for the connector phase and H_j for every other oriented "
            "piece: top and bottom q>=1 Abel arcs, q<1 arcs, finite/core "
            "shoulders, adjacent-cutoff and chart joins, and the two "
            "connector endpoint tracks. Then "
            "2*pi*kappa_j=H_j+C_j, |C_j|<pi/2, and kappa_j is a "
            "nonnegative integer."
        ),
        "reduced_horizontal_target": (
            "The strict signed inequality H_j<3*pi/2 is sufficient for "
            "H_j+C_j<2*pi, hence 0<=kappa_j<1 and kappa_j=0. The vertical "
            "connector is no longer an open arithmetic phase theorem; the "
            "remaining phase budget is the completely joined horizontal "
            "composite, including its endpoint tracks."
        ),
        "seam_boundary": (
            "The cap is proved for the exact normalized Xi first jet on the "
            "right edge, not for an unjoined Abel chart. Replacing the "
            "horizontal exact jet by the Abel proxy requires the contact-band "
            "gap and explicit nonvanishing endpoint tracks. Those tracks are "
            "part of H_j and may not be discarded or charged twice."
        ),
        "endpoint_only_guard": (
            "Endpoint agreement alone cannot cap an open-path phase: "
            "G_m(u)=exp(2*pi*i*m*u), 0<=u<=1, has G_m(0)=G_m(1)=1 but "
            "Delta arg G_m=2*pi*m. The connector proof instead uses the "
            "uniform pointwise cone |arg(V/V_0)|<3/5, which forbids this "
            "hidden winding."
        ),
        "remaining_target": (
            "Prove the Abel-scalar pointwise gap on every q>=1 horizontal "
            "chart and prove H_j<3*pi/2 after the top, bottom, q<1, finite "
            "shoulder, cutoff/chart, and connector-seam pieces are joined. "
            "No absolute-value sum over local carrier fluxes is licensed."
        ),
        "proof_boundary": (
            "The ray geometry, normalizer time-phase identity, one-saddle "
            "ellipse, uniform exact-H cone, strict connector phase cap, "
            "zero axis phase, reduced horizontal ledger, and endpoint-only "
            "winding guard are exact or inherited from certified source "
            "bounds. This does not prove the Xi Abel-scalar gap, the "
            "3*pi/2 horizontal-composite inequality, q<1 closure, finite "
            "shoulders or endpoint tracks, complete boundary composition, "
            "contact exclusion, Lambda<=0, PF-infinity, RH, or a Clay-prize "
            "conclusion."
        ),
        "diagnostics": audit,
    }


def build_rows(exact: dict) -> list[ConnectorRow]:
    return [
        ConnectorRow(
            "drcpc_00_pi_provenance",
            "source_provenance",
            "available_exact",
            "No new pi is introduced by the connector cap.",
            exact["pi_provenance"],
            "The rational pi bounds normalize an angle already present.",
        ),
        ConnectorRow(
            "drcpc_01_connector_geometry",
            "exact_geometry",
            "available_exact",
            "Every successor right edge lies on one dominant ray.",
            exact["connector_geometry"],
            "Valid for the ray-aligned successor schedule j>=0.",
        ),
        ConnectorRow(
            "drcpc_02_normalized_vector",
            "exact_homotopy",
            "available_exact",
            "The exact normalized first jet may represent the connector.",
            exact["normalized_vector"],
            "Closed-boundary degree is preserved; endpoint tracks remain.",
        ),
        ConnectorRow(
            "drcpc_03_normalizer_phase",
            "exact_identity",
            "available_exact",
            "The one-saddle phase moves by less than pi/64.",
            exact["normalizer_phase"],
            "Uses L>=101 and the exact Polymath normalizer.",
        ),
        ConnectorRow(
            "drcpc_04_leading_ellipse",
            "certified_bound",
            "ready_to_apply",
            "The leading connector jet has a uniform radial floor.",
            exact["leading_ellipse"],
            "Uses the certified dominant-saddle beta-prime bounds.",
        ),
        ConnectorRow(
            "drcpc_05_uniform_error",
            "certified_bound",
            "ready_to_apply",
            "The exact normalized jet stays within 4/15 of the ellipse.",
            exact["uniform_vector_error"],
            "Includes every cutoff transition on the dominant ray.",
            exact["diagnostics"],
        ),
        ConnectorRow(
            "drcpc_06_uniform_cone",
            "exact_geometry",
            "available_exact",
            "The exact connector has no hidden relative winding.",
            exact["uniform_cone"],
            "Uniform path control, not endpoint control, is essential.",
        ),
        ConnectorRow(
            "drcpc_07_leading_phase",
            "certified_bound",
            "ready_to_apply",
            "The leading ellipse changes phase by less than 331/1344.",
            exact["leading_phase_variation"],
            "The lambda range is used as an affine total-variation cap.",
        ),
        ConnectorRow(
            "drcpc_08_connector_cap",
            "proved_boundary_lemma",
            "ready_to_apply",
            "The exact right connector costs less than one quarter turn.",
            exact["connector_phase_cap"],
            "This is a phase cap, not merely boundary nonvanishing.",
            exact["diagnostics"],
        ),
        ConnectorRow(
            "drcpc_09_axis_phase",
            "proved_boundary_lemma",
            "ready_to_apply",
            "The symmetry-axis phase contribution is zero.",
            exact["axis_phase"],
            "Uses evenness and H_t(0)>0.",
        ),
        ConnectorRow(
            "drcpc_10_composite_ledger",
            "exact_reduction",
            "available_exact",
            "The successor winding splits into connector and complement.",
            exact["composite_ledger"],
            "Every seam is assigned exactly once.",
        ),
        ConnectorRow(
            "drcpc_11_horizontal_target",
            "conditional_reduction",
            "available_exact",
            "A 3*pi/2 horizontal bound now suffices.",
            exact["reduced_horizontal_target"],
            "Conditional on a compatible nonvanishing composite.",
        ),
        ConnectorRow(
            "drcpc_12_seam_boundary",
            "nonpromotion_guard",
            "guard_validated",
            "The connector cap cannot be copied to an unjoined Abel chart.",
            exact["seam_boundary"],
            "Endpoint tracks remain part of the open theorem.",
        ),
        ConnectorRow(
            "drcpc_13_endpoint_guard",
            "countermodel",
            "guard_validated",
            "Matching endpoint phases do not bound path winding.",
            exact["endpoint_only_guard"],
            "Generic path guard, not an Xi counterexample.",
        ),
        ConnectorRow(
            "drcpc_14_pointwise_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The Xi Abel-scalar gap remains open.",
            exact["remaining_target"],
            "The connector theorem supplies no horizontal gap.",
        ),
        ConnectorRow(
            "drcpc_15_horizontal_budget",
            "open_theorem_target",
            "not_ready_to_apply",
            "The joined horizontal 3*pi/2 budget remains open.",
            exact["remaining_target"],
            "It includes both connector endpoint tracks.",
        ),
        ConnectorRow(
            "drcpc_16_small_q_and_finite",
            "open_theorem_target",
            "not_ready_to_apply",
            "The q<1, finite-shoulder, and join pieces remain open.",
            exact["remaining_target"],
            "No complete boundary composition is inferred.",
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    rows = [asdict(row) for row in build_rows(exact)]
    return {
        "kind": STEM,
        "date": "2026-07-28",
        "status": (
            "exact dominant-ray right-connector phase cap below one quarter "
            "turn and reduced 3*pi/2 horizontal-composite phase ledger; "
            "the horizontal Xi inequalities remain open"
        ),
        "sources": {key: str(path.relative_to(REPO_ROOT)) for key, path in SOURCES.items()},
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "exact": exact,
        "rows": rows,
        "proof_boundary": exact["proof_boundary"],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return f"""# Newman Dominant-Ray Connector Phase Cap

Date: 2026-07-28

Status: exact right-connector phase cap and reduced horizontal phase
ledger; not a proof of the remaining Xi boundary inequalities,
`Lambda<=0`, PF-infinity, RH, or a Clay-prize result.

```text
work/rh_compute/results/{STEM}.json
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```

## Pi Provenance

{exact["pi_provenance"]}

## Ray Geometry

```text
{exact["connector_geometry"]}
{exact["normalized_vector"]}
```

## One-Saddle Ellipse

```text
{exact["normalizer_phase"]}
{exact["leading_ellipse"]}
```

## Uniform Cone

```text
{exact["uniform_vector_error"]}
{exact["uniform_cone"]}
```

## Connector Cap

```text
{exact["leading_phase_variation"]}
{exact["connector_phase_cap"]}
```

The bound is uniform across every cutoff transition. It controls the
continuous argument lift of the complete connector, not just its endpoints.

## Phase Ledger

```text
{exact["axis_phase"]}
{exact["composite_ledger"]}
{exact["reduced_horizontal_target"]}
```

The vertical arithmetic connector has therefore been discharged. The two
connector endpoint tracks remain in the horizontal composite so that no
open-path phase is silently lost.

## Seam Guard

{exact["seam_boundary"]}

## Endpoint-Only Guard

```text
{exact["endpoint_only_guard"]}
```

This is a generic path countermodel, not an Xi counterexample.

## Remaining Target

```text
{exact["remaining_target"]}
```

## Boundary

{exact["proof_boundary"]}
"""


def write_artifact(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    write_artifact(args.out, artifact)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman dominant-ray connector phase cap: "
        "17 rows, 1 uniform cone, 1 strict quarter-turn connector cap, "
        "1 reduced 3*pi/2 horizontal ledger, "
        "1 endpoint-only winding guard, 3 open obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
