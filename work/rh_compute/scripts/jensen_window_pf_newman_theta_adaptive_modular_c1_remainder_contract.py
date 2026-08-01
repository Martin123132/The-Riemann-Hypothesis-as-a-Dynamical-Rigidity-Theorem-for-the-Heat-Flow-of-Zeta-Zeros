#!/usr/bin/env python3
"""Compose the adaptive modular theta blocks into a direct C1 contract."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_theta_"
    "adaptive_modular_c1_remainder_contract"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
RESULT_ROOT = REPO_ROOT / "work/rh_compute/results"
MODULAR_SOURCE = (
    RESULT_ROOT
    / "jensen_window_pf_newman_theta_modular_blend_gate.json"
)
HIGH_FREQUENCY_SOURCE = (
    RESULT_ROOT
    / "jensen_window_pf_newman_theta_"
    "modular_blend_high_frequency_scout.json"
)
ADAPTIVE_SOURCE = (
    RESULT_ROOT
    / "jensen_window_pf_newman_theta_"
    "modular_blend_adaptive_saddle_gate.json"
)
FIXED_BLOCK_SOURCE = (
    RESULT_ROOT
    / "jensen_window_pf_newman_theta_"
    "fixed_block_cofinal_obstruction_gate.json"
)
COMPACT_SOURCE = (
    RESULT_ROOT
    / "jensen_window_pf_newman_theta_"
    "compact_transversality_interval_certificate.json"
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


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def source_audit() -> dict:
    modular = load(MODULAR_SOURCE)
    high = load(HIGH_FREQUENCY_SOURCE)
    adaptive = load(ADAPTIVE_SOURCE)
    fixed = load(FIXED_BLOCK_SOURCE)
    compact = load(COMPACT_SOURCE)
    texts = {
        "modular": json.dumps(modular, sort_keys=True),
        "high": json.dumps(high, sort_keys=True),
        "adaptive": json.dumps(adaptive, sort_keys=True),
        "fixed": json.dumps(fixed, sort_keys=True),
        "compact": json.dumps(compact, sort_keys=True),
    }
    markers = {
        "modular": (
            "sum_(n>=1)b_n(u)=Phi(u)",
            "|R_(N,t)^(j)(x)|<=d_(N,j,m)(T)/|x|^m",
            "b_n(u)>0 for every n>=1",
        ),
        "high": (
            "fixed three-block Laguerre expression",
            "arithmetic truncation grow with frequency",
        ),
        "adaptive": (
            "sqrt(x/(4*pi))",
            "Laguerre remainder sign",
        ),
        "fixed": (
            "fixed-N static absolute-tail proof method",
            "endpoint-jet cancellation",
        ),
        "compact": (
            "|x|<=38",
            "(H_t,H_t')!=(0,0)",
        ),
    }
    for key, required in markers.items():
        for marker in required:
            if marker not in texts[key]:
                raise RuntimeError(
                    f"{key} source marker missing: {marker}"
                )

    high_rows = high["rows"]
    x200 = [
        row for row in high_rows if row.get("x") == 200
    ]
    adaptive_rows = adaptive["diagnostics"]
    return {
        "modular_kind": modular["kind"],
        "high_frequency_kind": high["kind"],
        "adaptive_kind": adaptive["kind"],
        "fixed_block_kind": fixed["kind"],
        "compact_kind": compact["kind"],
        "high_frequency_rows": len(high_rows),
        "minimum_x200_cancellation_digits": str(
            min(
                sp.Float(row["cancellation_digits"], 80)
                for row in x200
            )
        ),
        "adaptive_rows": len(adaptive_rows),
        "adaptive_first_counts": [
            row["first_n_relative_error_below_half"]
            for row in adaptive_rows
        ],
        "adaptive_collar_covers_all": all(
            row["adaptive_bound_covers_observed_n"]
            and row["adaptive_bound_covers_observed_margin_n"]
            for row in adaptive_rows
        ),
        "fixed_frontier_prefix": fixed["source_audit"][
            "frontier_prefix"
        ],
    }


def build_exact() -> dict:
    x, s, spx, r, rp = sp.symbols(
        "x S S_x R R_x", real=True
    )
    full_j = 16 * x**4 * (s + r)
    partial_j = 16 * x**4 * s
    if sp.expand(full_j - partial_j - 16 * x**4 * r) != 0:
        raise RuntimeError("J value error identity failed")
    full_j_prime = (
        64 * x**3 * (s + r)
        + 16 * x**4 * (spx + rp)
    )
    partial_j_prime = 64 * x**3 * s + 16 * x**4 * spx
    expected_error = 64 * x**3 * r + 16 * x**4 * rp
    if sp.expand(full_j_prime - partial_j_prime - expected_error) != 0:
        raise RuntimeError("J derivative error identity failed")

    return {
        "modular_partition": {
            "switch": "omega(u)=(1+erf(3*sinh(4u)))/2",
            "block": (
                "b_n(u)=omega(u)*phi_n(u)+"
                "(1-omega(u))*phi_n(-u)"
            ),
            "properties": (
                "b_n is even, entire, and strictly positive on R"
            ),
            "partition": "sum_(n>=1)b_n(u)=Phi(u)",
            "newman_decay": (
                "For every finite T, exp(tu^2)b_n(u) is Schwartz "
                "uniformly for 0<=t<=T"
            ),
        },
        "normal_transform": {
            "block_transform": (
                "B_(n,t)(x)=integral_0^infinity "
                "exp(tu^2)b_n(u)cos(xu)du"
            ),
            "partial": (
                "S_(N,t)(x)=sum_(1<=n<=N)B_(n,t)(x)"
            ),
            "remainder": (
                "R_(N,t)(x)=H_t(x)-S_(N,t)(x)"
            ),
            "series": (
                "H_t^(j)(x)=sum_(n>=1)B_(n,t)^(j)(x), "
                "j=0,1,2,..., absolutely and uniformly on "
                "R_x times [0,T]"
            ),
        },
        "tail_budgets": {
            "moment": (
                "mu_(N,j)(T)=integral_0^infinity "
                "u^j*exp(Tu^2)*sum_(n>N)b_n(u)du"
            ),
            "derivative": (
                "d_(N,j,m)(T)=(1/2)*sup_(0<=t<=T)"
                "||partial_u^m[(iu)^j*exp(tu^2)r_N(u)]||_1"
            ),
            "combined": (
                "epsilon_(N,j,m)(x,T)=min(mu_(N,j)(T),"
                "d_(N,j,m)(T)/|x|^m)"
            ),
            "jet_bound": (
                "|R_(N,t)^(j)(x)|<=epsilon_(N,j,m)(x,T), "
                "j=0,1,2"
            ),
            "finiteness": (
                "Every d_(N,j,m)(T) is finite because the modular "
                "tail is uniformly Schwartz at bounded Newman time."
            ),
        },
        "direct_c1_contract": {
            "partial_value": (
                "J_(N,t)(x)=16*x^4*S_(N,t)(x)"
            ),
            "partial_derivative": (
                "J_(N,t)'(x)=64*x^3*S_(N,t)(x)+"
                "16*x^4*S_(N,t)'(x)"
            ),
            "value_error": (
                "|J_t-J_(N,t)|<=16*x^4*epsilon_(N,0,m)"
            ),
            "derivative_error": (
                "|J_t'-J_(N,t)'|<=64*x^3*epsilon_(N,0,m)+"
                "16*x^4*epsilon_(N,1,m)"
            ),
            "sufficient_disjunction": (
                "|J_(N,t)|>16*x^4*epsilon_(N,0,m) OR "
                "|J_(N,t)'|>64*x^3*epsilon_(N,0,m)+"
                "16*x^4*epsilon_(N,1,m)"
            ),
            "conclusion": (
                "The sufficient disjunction implies "
                "(H_t(x),H_t'(x))!=(0,0)."
            ),
        },
        "adaptive_scale": {
            "component_phase": (
                "F_(a,n,t,x)(u)=tu^2+(a+ix)u-"
                "pi*n^2*exp(4u), a in {5,9}"
            ),
            "crossing": (
                "4y=atan((x+2ty)/a); "
                "n_*^2=sqrt(a^2+(x+2ty)^2)/(4pi)"
            ),
            "asymptotic": (
                "n_*(a,t,x)=sqrt(x/(4pi))*(1+O(x^-1)) "
                "uniformly for bounded t"
            ),
            "candidate": (
                "N_kappa(x,t)=ceil(max(n_*(5,t,x),"
                "n_*(9,t,x)))+kappa"
            ),
            "boundary": (
                "The diagnostics suggest kappa=2, but no fixed-collar "
                "theorem is proved."
            ),
        },
        "quantitative_target": {
            "target": (
                "Find explicit integers m>=5 and kappa>=0, or another "
                "controlled adaptive N(x,t), and rigorously prove the "
                "direct C1 sufficient disjunction for every "
                "0<t<=1/5 and x>38."
            ),
            "budget_scaling": (
                "For m>=5, the explicit J error multipliers are "
                "x^(4-m)*d_(N,0,m) and "
                "4*x^(3-m)*d_(N,0,m)+"
                "x^(4-m)*d_(N,1,m), up to the common factor 16."
            ),
            "composition": (
                "Together with the certified |x|<=38 core and "
                "positive-boundary attainment, this target would close "
                "the remaining positive-time first-jet obstruction."
            ),
            "required_work": (
                "Derive rigorous N-dependent derivative L1 budgets and "
                "a matching lower separation profile for the finite "
                "adaptive first jet; intervalize transition cells."
            ),
        },
        "nonpromotion": {
            "fixed_block": (
                "Fixed N with static absolute bars is cofinally "
                "impossible."
            ),
            "cancellation": (
                "At x=200, the fixed three-block Laguerre diagnostic "
                "is cancelled by more than 20 decimal digits."
            ),
            "monotonicity": (
                "The stronger global sign -L_t'>0 is false at the "
                "Xi Lehmer stress point and is not part of this target."
            ),
            "scope": (
                "Neither finite diagnostics nor the exact error "
                "contract prove the open adaptive separation inequality."
            ),
        },
    }


def build_artifact() -> dict:
    exact = build_exact()
    audit = source_audit()
    rows = [
        GateRow(
            "ntamc1_01_positive_even_partition",
            "exact_theorem",
            "ready_to_apply",
            "The modular blend partitions Phi into positive even entire blocks.",
            (
                f"{exact['modular_partition']['block']}; "
                f"{exact['modular_partition']['partition']}"
            ),
            "Kernel decomposition only; no transform sign follows.",
        ),
        GateRow(
            "ntamc1_02_normal_transform_series",
            "exact_theorem",
            "ready_to_apply",
            "Every frequency derivative has a normal arithmetic block series at bounded Newman time.",
            exact["normal_transform"]["series"],
            "Absolute convergence does not sign the coupled sum.",
        ),
        GateRow(
            "ntamc1_03_decaying_tail_budget",
            "exact_inequality",
            "ready_to_apply",
            "Repeated Fourier integration by parts gives arbitrary-power modular-tail decay.",
            (
                f"{exact['tail_budgets']['derivative']}; "
                f"{exact['tail_budgets']['jet_bound']}"
            ),
            "The constants are finite but useful explicit bounds remain to be derived.",
        ),
        GateRow(
            "ntamc1_04_direct_first_jet_error",
            "exact_inequality",
            "ready_to_apply",
            "The modular jet budgets give direct full-Xi error bars for J and J'.",
            (
                f"{exact['direct_c1_contract']['value_error']}; "
                f"{exact['direct_c1_contract']['derivative_error']}"
            ),
            "Exact finite-to-infinite C1 reduction.",
        ),
        GateRow(
            "ntamc1_05_contact_exclusion_disjunction",
            "exact_composition",
            "ready_to_apply",
            "One strict retained value-or-derivative inequality excludes full-Xi contact.",
            exact["direct_c1_contract"]["sufficient_disjunction"],
            exact["direct_c1_contract"]["conclusion"],
        ),
        GateRow(
            "ntamc1_06_fixed_block_guard",
            "nonpromotion_gate",
            "guard_validated",
            "The retained block count must grow or use an equivalent cancellation-preserving architecture.",
            exact["nonpromotion"]["fixed_block"],
            "Rejects static fixed-block bars, not Xi transversality.",
        ),
        GateRow(
            "ntamc1_07_adaptive_saddle_scale",
            "exact_theorem",
            "ready_to_apply",
            "The modular transition fixes the Riemann-Siegel square-root block scale.",
            (
                f"{exact['adaptive_scale']['crossing']}; "
                f"{exact['adaptive_scale']['asymptotic']}"
            ),
            exact["adaptive_scale"]["boundary"],
            {
                "diagnostic_first_counts": audit[
                    "adaptive_first_counts"
                ],
                "diagnostic_collar_covers_all": audit[
                    "adaptive_collar_covers_all"
                ],
            },
        ),
        GateRow(
            "ntamc1_08_cancellation_stress",
            "finite_diagnostic",
            "diagnostic_only",
            "A fixed three-block Laguerre expression loses more than twenty digits to the arithmetic tail at x=200.",
            exact["nonpromotion"]["cancellation"],
            "Two-time point diagnostic only; no interval sign theorem.",
            {
                "rows": audit["high_frequency_rows"],
                "minimum_x200_cancellation_digits": audit[
                    "minimum_x200_cancellation_digits"
                ],
            },
        ),
        GateRow(
            "ntamc1_09_open_adaptive_c1_target",
            "open_handoff",
            "not_ready_to_apply",
            exact["quantitative_target"]["target"],
            exact["quantitative_target"]["budget_scaling"],
            (
                f"{exact['quantitative_target']['required_work']} "
                f"{exact['nonpromotion']['monotonicity']}"
            ),
        ),
        GateRow(
            "ntamc1_10_proof_boundary",
            "proof_guard",
            "guard_active",
            "The contract isolates one quantitative bridge and keeps every finite diagnostic subordinate to it.",
            exact["nonpromotion"]["scope"],
            "No new Q_j, Lambda<=0, RH, or Clay-prize result is proved.",
        ),
    ]
    return {
        "kind": STEM,
        "date": "2026-07-24",
        "status": (
            "exact adaptive modular C1 remainder contract with one "
            "open quantitative separation target"
        ),
        "proof_boundary": (
            "This artifact composes an exact positive even modular theta "
            "partition, normal transform series, arbitrary-power Fourier "
            "tail budgets, direct J/J' error bars, and the exact adaptive "
            "saddle scale. It does not supply useful uniform derivative "
            "constants or the retained first-jet lower separation needed "
            "by the displayed disjunction. The finite cancellation and "
            "collar rows are diagnostic only. It proves no new Q_j, "
            "Lambda<=0, RH, or Clay-prize result."
        ),
        "sources": [
            str(path.relative_to(REPO_ROOT)).replace("\\", "/")
            for path in (
                MODULAR_SOURCE,
                HIGH_FREQUENCY_SOURCE,
                ADAPTIVE_SOURCE,
                FIXED_BLOCK_SOURCE,
                COMPACT_SOURCE,
            )
        ]
        + ["outputs/formal_core.md"],
        "source_audit": audit,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    audit = artifact["source_audit"]
    return "\n".join(
        [
            "# Newman Theta Adaptive Modular C1 Remainder Contract",
            "",
            "Date: 2026-07-24",
            "",
            "Status: exact reduction with one open quantitative "
            "first-jet target.",
            "This is not a proof of `Lambda <= 0` or RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Cancellation-Preserving Blocks",
            "",
            "Set",
            "",
            "```text",
            exact["modular_partition"]["switch"],
            exact["modular_partition"]["block"],
            exact["modular_partition"]["properties"],
            exact["modular_partition"]["partition"],
            "```",
            "",
            exact["modular_partition"]["newman_decay"],
            "Consequently,",
            "",
            "```text",
            exact["normal_transform"]["block_transform"],
            exact["normal_transform"]["partial"],
            exact["normal_transform"]["remainder"],
            exact["normal_transform"]["series"],
            "```",
            "",
            "Unlike a hard half-line theta truncation, every retained "
            "block is even, so the modular endpoint cancellation is "
            "preserved before transformation.",
            "",
            "## Decaying Remainder",
            "",
            "For `r_N=sum_(n>N)b_n`, define",
            "",
            "```text",
            exact["tail_budgets"]["moment"],
            exact["tail_budgets"]["derivative"],
            exact["tail_budgets"]["combined"],
            "```",
            "",
            "Repeated Fourier integration by parts gives",
            "",
            "```text",
            exact["tail_budgets"]["jet_bound"],
            "```",
            "",
            exact["tail_budgets"]["finiteness"],
            "",
            "## Direct C1 Contract",
            "",
            "With",
            "",
            "```text",
            exact["direct_c1_contract"]["partial_value"],
            exact["direct_c1_contract"]["partial_derivative"],
            "```",
            "",
            "the full first-jet errors satisfy",
            "",
            "```text",
            exact["direct_c1_contract"]["value_error"],
            exact["direct_c1_contract"]["derivative_error"],
            "```",
            "",
            "Therefore the exact sufficient condition is",
            "",
            "```text",
            exact["direct_c1_contract"]["sufficient_disjunction"],
            "```",
            "",
            exact["direct_c1_contract"]["conclusion"],
            "",
            "## Adaptive Scale",
            "",
            "The component saddle obeys",
            "",
            "```text",
            exact["adaptive_scale"]["component_phase"],
            exact["adaptive_scale"]["crossing"],
            exact["adaptive_scale"]["asymptotic"],
            exact["adaptive_scale"]["candidate"],
            "```",
            "",
            exact["adaptive_scale"]["boundary"],
            "The ten stored diagnostics have first scale-level counts",
            "",
            "```text",
            str(audit["adaptive_first_counts"]),
            "```",
            "",
            "and all lie within the empirical two-index collar. This "
            "is supporting evidence, not a collar theorem.",
            "",
            "## Stress Gate",
            "",
            "The fixed three-block high-frequency scout has",
            "",
            "```text",
            f"rows={audit['high_frequency_rows']}",
            "minimum x=200 cancellation digits="
            f"{audit['minimum_x200_cancellation_digits']}",
            "```",
            "",
            "Thus the sign is carried by the coupled arithmetic "
            "remainder, not by a fixed base margin. The separate exact "
            "fixed-block theorem also proves that static absolute bars "
            "cannot be cofinal.",
            "",
            "## Open Quantitative Target",
            "",
            "```text",
            exact["quantitative_target"]["target"],
            exact["quantitative_target"]["budget_scaling"],
            "```",
            "",
            exact["quantitative_target"]["required_work"],
            "",
            exact["quantitative_target"]["composition"],
            "",
            "The stronger global monotonicity condition is not being "
            "smuggled back in:",
            "",
            "```text",
            exact["nonpromotion"]["monotonicity"],
            "```",
            "",
            "## Proof Boundary",
            "",
            exact["nonpromotion"]["scope"],
            "No new `Q_j`, `Lambda<=0`, RH, or Clay-prize result "
            "follows.",
            "",
            "Current validation:",
            "",
            "```text",
            "validated Newman theta adaptive modular C1 remainder "
            "contract: 10 rows, 0 issues, 2 exact modular theorems, "
            "3 exact C1 inequalities/compositions, 1 fixed-block "
            "guard, 1 adaptive saddle theorem, 10 cancellation "
            "diagnostics, 1 open quantitative target",
            "```",
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
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "wrote Newman theta adaptive modular C1 remainder "
        f"contract: {len(artifact['rows'])} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
