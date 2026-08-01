#!/usr/bin/env python3
"""Build the positive-boundary diagonal compact-exhaustion gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_positive_boundary_diagonal_exhaustion_gate"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
BOUNDARY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_positive_boundary_attainment_lemma.json"
)
COMPACT_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_compact_transversality_"
    "interval_certificate.json"
)
DELTA_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_positive_boundary_delta_localization_gate.json"
)
FIELD_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_local_odd_count_reduction_lemma.json"
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def symbolic_countermodel() -> dict:
    z, c, b, tau, eps = sp.symbols(
        "z c b tau eps", real=True, nonzero=True
    )
    a2 = sp.cancel(c**2 * (3 - b * c) / (1 - b * c))
    polynomial = sp.expand((z**2 - c**2) ** 2 * (z**2 - a2))
    regular = sp.cancel(polynomial / (z - c) ** 2)
    field = sp.factor(
        sp.diff(regular, z).subs(z, c) / regular.subs(z, c)
    )
    if field != b:
        raise RuntimeError(f"arbitrary-center field identity failed: {field}")

    stiffness = sp.factor(
        sp.Rational(1, 2) / c**2
        + 2 * (c**2 + a2) / (c**2 - a2) ** 2
    )
    expected_stiffness = b**2 - 3 * b / c + sp.Rational(5, 2) / c**2
    if sp.simplify(stiffness - expected_stiffness) != 0:
        raise RuntimeError("arbitrary-center stiffness identity failed")

    heat = sp.expand(
        sum(
            (-tau) ** k
            * sp.diff(polynomial, z, 2 * k)
            / sp.factorial(k)
            for k in range(4)
        )
    )
    if sp.simplify(sp.diff(heat, tau) + sp.diff(heat, z, 2)) != 0:
        raise RuntimeError("polynomial backward-heat identity failed")

    for sign in (-1, 1):
        displacement = sign * sp.sqrt(2) * eps + 2 * b * eps**2
        residual = sp.series(
            heat.subs({z: c + displacement, tau: eps**2}),
            eps,
            0,
            4,
        ).removeO()
        if sp.simplify(residual) != 0:
            raise RuntimeError(
                f"Puiseux collision expansion failed for sign {sign}"
            )

    return {
        "a_squared": "a_c^2=c^2*(3-b*c)/(1-b*c)",
        "polynomial": "P_c(z)=(z^2-c^2)^2*(z^2-a_c^2)",
        "field": (
            "B_c=1/c+2*c/(c^2-a_c^2)=b at the positive double zero c"
        ),
        "stiffness": (
            "K_c=1/(2*c^2)+1/(c-a_c)^2+1/(c+a_c)^2"
            "=b^2-3*b/c+5/(2*c^2)"
        ),
        "heat_flow": (
            "F_(lambda+tau)=exp(-tau*d_z^2)P_c solves "
            "partial_tau F=-partial_z^2 F"
        ),
        "local_split": (
            "z_+/-=c+/-sqrt(2*tau)+2*b*tau+O(tau^(3/2))"
        ),
        "hyperbolicity": (
            "For tau>=0, exp(-tau*D^2) preserves real-rootedness because "
            "it is the coefficientwise limit of "
            "[(1-sqrt(tau/n)D)(1+sqrt(tau/n)D)]^n, and each first-order "
            "factor preserves real-rootedness by derivative interlacing"
        ),
        "boundary": (
            "For tau<0 sufficiently close to zero the two roots near c are "
            "nonreal, while every root is real for tau>=0; after the time "
            "shift, this polynomial flow has Newman-style boundary lambda"
        ),
        "classical_specialization": (
            "For b=-pi/8, a_c^2=c^2*(24+pi*c)/(8+pi*c), "
            "B_c=-pi/8, and "
            "K_c=pi^2/64+3*pi/(8*c)+5/(2*c^2)->pi^2/64"
        ),
    }


def build_exact() -> dict:
    return {
        "boundary_attainment": (
            "If Lambda>0, then 0<Lambda<=1/5 and H_Lambda has a finite "
            "real multiple zero c"
        ),
        "abstract_exhaustion": (
            "For any sequences 0<delta_j<=1/5 with delta_j->0 and "
            "R_j>0 with R_j->infinity, Lambda<=0 iff "
            "(H_t,H_t')!=(0,0) on "
            "[delta_j,1/5]x[-R_j,R_j] for every j"
        ),
        "eventual_witness": (
            "If Lambda>0 with boundary collision c, then every sufficiently "
            "large j satisfies delta_j<=Lambda and R_j>=|c|, so the same "
            "fixed point (Lambda,c) lies in every later rectangle"
        ),
        "linear_exhaustion": (
            "With delta_j=1/(5j) and R_j=38+j, Lambda<=0 iff "
            "(H_t,H_t')!=(0,0) on "
            "[1/(5j),1/5]x[-(38+j),38+j] for every j>=1"
        ),
        "explicit_witness_index": (
            "Under Lambda>0, J=ceil(max(1,1/(5*Lambda),|c|-38)) "
            "places (Lambda,c) in every linear-exhaustion rectangle j>=J"
        ),
        "compact_shell": (
            "Because (H_t,H_t')!=(0,0) is already certified for "
            "0<=t<=1/5 and |x|<=38, the linear exhaustion only leaves "
            "1/(5j)<=t<=1/5 and 38<|x|<=38+j"
        ),
        "compact_minimum": (
            "For each fixed j, no common zero on its compact rectangle is "
            "equivalent to min[H_t(x)^2+H_t'(x)^2]>0 there"
        ),
        "independent_rates": (
            "No coupling between delta_j and R_j is logically required: "
            "R_j may diverge arbitrarily slowly, for example "
            "R_j=38+log(1+j), independently of delta_j->0"
        ),
        "old_radius_reclassified": (
            "The radius 4*pi*exp(25/delta) remains a valid stronger "
            "per-delta localization that contains every possible contact "
            "above that time floor, but it is not required by the cofinal "
            "positive-boundary contradiction"
        ),
        "countermodel": symbolic_countermodel(),
        "proof_handoff": (
            "Use the independent diagonal exhaustion as the minimal compact "
            "endgame: prove delta-dependent C1 separation on its growing "
            "finite shells, or derive an Xi-specific nodal-flux/resolvent "
            "continuation across them. Do not impose an exponential "
            "radius/time-floor coupling, and do not infer any radius bound "
            "from local field, drift, or stiffness data alone"
        ),
    }


def source_audit() -> dict:
    sources = {
        "boundary": json.loads(BOUNDARY_SOURCE.read_text(encoding="utf-8")),
        "compact": json.loads(COMPACT_SOURCE.read_text(encoding="utf-8")),
        "delta": json.loads(DELTA_SOURCE.read_text(encoding="utf-8")),
        "field": json.loads(FIELD_SOURCE.read_text(encoding="utf-8")),
    }
    texts = {key: json.dumps(value) for key, value in sources.items()}
    markers = {
        "boundary": (
            "H_Lambda has a finite real multiple zero",
            "Lambda<=0 if and only if H_t has only simple zeros",
        ),
        "compact": (
            "|x|<=38",
            "0<=t<=1/5",
            "no contact",
        ),
        "delta": (
            "R_delta=4*pi*exp(25/delta)",
            "delta_j=1/(5j)",
        ),
        "field": (
            "B(c)=-pi/8",
            "positive square-root birth",
        ),
    }
    for key, required in markers.items():
        for marker in required:
            if marker not in texts[key]:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {
        f"{key}_kind": value["kind"]
        for key, value in sources.items()
    }


def linear_examples() -> list[dict]:
    return [
        {
            "j": j,
            "delta_j": f"1/{5 * j}",
            "radius": str(38 + j),
            "new_shell": f"38<|x|<={38 + j}",
        }
        for j in (1, 2, 5, 10, 100)
    ]


def field_examples() -> list[dict]:
    b = -math.pi / 8
    rows: list[dict] = []
    for c in (1.0, 10.0, 100.0, 1000.0):
        a = c * math.sqrt((3 - b * c) / (1 - b * c))
        stiffness = b * b - 3 * b / c + 2.5 / (c * c)
        rows.append(
            {
                "c": c,
                "a_c": a,
                "a_c_minus_c": a - c,
                "field": b,
                "stiffness": stiffness,
            }
        )
    return rows


def build_rows(exact: dict, audit: dict) -> list[GateRow]:
    countermodel = exact["countermodel"]
    return [
        GateRow(
            "npbdeg_01_boundary_attainment",
            "published_input",
            "ready_to_apply",
            "A positive Newman boundary supplies one fixed finite collision.",
            exact["boundary_attainment"],
            "Uses the separately validated positive-boundary attainment lemma.",
            audit,
        ),
        GateRow(
            "npbdeg_02_abstract_exhaustion",
            "exact_equivalence",
            "ready_to_apply",
            "Any independent compact exhaustion detects the fixed boundary collision.",
            f"{exact['abstract_exhaustion']}; {exact['eventual_witness']}",
            "This is an equivalence of the full cofinal family, not a certificate for any member.",
        ),
        GateRow(
            "npbdeg_03_linear_exhaustion",
            "exact_equivalence",
            "ready_to_apply",
            "A linear radius and reciprocal time floor are an explicit sufficient exhaustion.",
            f"{exact['linear_exhaustion']}; {exact['explicit_witness_index']}",
            "Every rectangle remains an open Xi transversality problem outside the certified core.",
            linear_examples(),
        ),
        GateRow(
            "npbdeg_04_compact_shell",
            "exact_composition",
            "ready_to_apply",
            "The existing compact theorem removes the fixed inner core from every stage.",
            exact["compact_shell"],
            "Uses the 1,900-box interval certificate and exact origin collar.",
        ),
        GateRow(
            "npbdeg_05_compact_minimum",
            "exact_reduction",
            "ready_to_apply",
            "Each diagonal stage is an ordinary positive first-jet minimum problem.",
            exact["compact_minimum"],
            "No new positive lower bound is proved.",
        ),
        GateRow(
            "npbdeg_06_independent_rates",
            "nonpromotion_gate",
            "guard_validated",
            "The exponential radius/time-floor coupling is sufficient but not logically necessary.",
            f"{exact['independent_rates']}; {exact['old_radius_reclassified']}",
            "This does not make the infinitely many compact stages automatic.",
        ),
        GateRow(
            "npbdeg_07_arbitrary_center_field",
            "exact_countermodel",
            "guard_validated",
            "The classical collision field can occur at an arbitrary prescribed height.",
            (
                f"{countermodel['a_squared']}; {countermodel['polynomial']}; "
                f"{countermodel['field']}; {countermodel['stiffness']}"
            ),
            "Generic even polynomial model, not an Xi representation.",
            field_examples(),
        ),
        GateRow(
            "npbdeg_08_shifted_boundary_flow",
            "exact_countermodel",
            "guard_validated",
            "The same arbitrary-height model can be placed at any positive backward-heat boundary time.",
            (
                f"{countermodel['heat_flow']}; {countermodel['local_split']}; "
                f"{countermodel['hyperbolicity']}; {countermodel['boundary']}"
            ),
            "Uses only real-rooted polynomial interlacing and local double-root splitting.",
        ),
        GateRow(
            "npbdeg_09_local_field_guard",
            "nonpromotion_gate",
            "guard_validated",
            "Classical field, drift, and bounded positive stiffness cannot recover a height/time coupling.",
            countermodel["classical_specialization"],
            "An Xi-specific global invariant or arithmetic estimate is still required.",
        ),
        GateRow(
            "npbdeg_10_proof_handoff",
            "open_theorem_target",
            "not_ready_to_apply",
            "Close the independent diagonal shells by Xi-specific transversality or continuation.",
            exact["proof_handoff"],
            "Open for the unbounded family; this gate proves no new Xi shell beyond |x|<=38.",
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    audit = source_audit()
    rows = build_rows(exact, audit)
    return {
        "kind": STEM,
        "date": "2026-07-24",
        "status": (
            "exact independent diagonal compact-exhaustion theorem and "
            "arbitrary-height local-field route guard; the unbounded Xi "
            "shell family remains open and this is not a proof of "
            "Lambda<=0 or RH"
        ),
        "proof_boundary": (
            "The gate removes an unnecessary exponential coupling between "
            "the lower time cutoff and spatial radius. It does not certify "
            "any new Xi shell, positive-time simplicity, Lambda<=0, or RH."
        ),
        "sources": [
            str(path.relative_to(REPO_ROOT)).replace("\\", "/")
            for path in (
                BOUNDARY_SOURCE,
                COMPACT_SOURCE,
                DELTA_SOURCE,
                FIELD_SOURCE,
            )
        ],
        "source_audit": audit,
        "exact": exact,
        "linear_examples": linear_examples(),
        "field_examples": field_examples(),
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    countermodel = exact["countermodel"]
    lines = [
        "# Newman Positive-Boundary Diagonal Exhaustion Gate",
        "",
        "Date: 2026-07-24",
        "",
        "Status: exact independent compact-exhaustion theorem and route guard.",
        "The growing Xi shell family remains open; this is not a proof of",
        "`Lambda <= 0` or RH.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_newman_positive_boundary_diagonal_exhaustion_gate.json",
        "python work/rh_compute/scripts/jensen_window_pf_newman_positive_boundary_diagonal_exhaustion_gate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_newman_positive_boundary_diagonal_exhaustion_gate.py",
        "```",
        "",
        "## Diagonal Exhaustion",
        "",
        "Positive-boundary attainment gives",
        "",
        "```text",
        exact["boundary_attainment"],
        "```",
        "",
        "The fixed finite witness changes the quantifier geometry:",
        "",
        "```text",
        exact["abstract_exhaustion"],
        exact["eventual_witness"],
        "```",
        "",
        "A particularly simple choice, composed with the certified",
        "`|x|<=38` core, is",
        "",
        "```text",
        exact["linear_exhaustion"],
        exact["explicit_witness_index"],
        exact["compact_shell"],
        "```",
        "",
        "| j | delta_j | radius | remaining shell |",
        "|---:|---:|---:|:---|",
    ]
    for row in artifact["linear_examples"]:
        lines.append(
            f"| {row['j']} | {row['delta_j']} | {row['radius']} | "
            f"{row['new_shell']} |"
        )
    lines.extend(
        [
            "",
            exact["compact_minimum"],
            "",
            "## Quantifier Correction",
            "",
            "```text",
            exact["independent_rates"],
            exact["old_radius_reclassified"],
            "```",
            "",
            "The earlier `4*pi*exp(25/delta)` radius remains useful when one",
            "wants a single fixed-delta rectangle containing every possible",
            "contact. It is not the minimal cofinal contradiction. One fixed",
            "finite boundary contact is eventually captured by any spatial",
            "radius tending to infinity.",
            "",
            "## Arbitrary-Height Field Guard",
            "",
            "Fix `c>0` and a target field `b<0`. Put",
            "",
            "```text",
            countermodel["a_squared"],
            countermodel["polynomial"],
            countermodel["field"],
            countermodel["stiffness"],
            "```",
            "",
            "The shifted backward-heat flow satisfies",
            "",
            "```text",
            countermodel["heat_flow"],
            countermodel["local_split"],
            countermodel["hyperbolicity"],
            countermodel["boundary"],
            "```",
            "",
            "At the classical field value:",
            "",
            "```text",
            countermodel["classical_specialization"],
            "```",
            "",
            "| c | a_c | a_c-c | B_c | K_c |",
            "|---:|---:|---:|---:|---:|",
        ]
    )
    for row in artifact["field_examples"]:
        lines.append(
            f"| {row['c']:.0f} | {row['a_c']:.12g} | "
            f"{row['a_c_minus_c']:.12g} | {row['field']:.12g} | "
            f"{row['stiffness']:.12g} |"
        )
    lines.extend(
        [
            "",
            "Thus local field, center drift, positive stiffness, evenness,",
            "and backward-heat evolution do not bind collision height to",
            "boundary time. Any such bound must use global Xi structure.",
            "",
            "## Live Handoff",
            "",
            "```text",
            exact["proof_handoff"],
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman positive-boundary diagonal-exhaustion gate: "
        f"{len(artifact['rows'])} rows, 5 exact exhaustion reductions, "
        "1 compact-shell composition, 1 arbitrary-height classical-field "
        "boundary countermodel, 1 independent-rate open target"
    )


if __name__ == "__main__":
    main()
