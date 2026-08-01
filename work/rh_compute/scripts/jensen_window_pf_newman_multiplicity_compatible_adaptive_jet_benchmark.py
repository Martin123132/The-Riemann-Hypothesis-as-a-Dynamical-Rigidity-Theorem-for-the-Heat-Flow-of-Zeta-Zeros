#!/usr/bin/env python3
"""Build the multiplicity-compatible adaptive-jet benchmark."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_multiplicity_compatible_adaptive_jet_benchmark"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{STEM}.md"
SUCCESSOR_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma.json"
)
SINGLE_CARRIER_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_single_carrier_adiabatic_benchmark.json"
)
CROSSING_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_crossing_slope_gap_reduction.json"
)
ATTAINMENT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_positive_boundary_attainment_lemma.json"
)
DELTA_LOCALIZATION_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_positive_boundary_delta_localization_gate.json"
)
DATE = "2026-07-25"


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def rows() -> list[dict]:
    return [
        {
            "id": "fixed_x_uniform_floor_guard",
            "kind": "exact_implication",
            "claim": (
                "A positive first-jet floor uniform as t decreases to zero "
                "at fixed x persists at t=0 and therefore excludes a "
                "multiple endpoint zero."
            ),
            "formula": (
                "If V(t,x)->V(0,x) and ||V(t,x)||>=b(x)>0 for every "
                "0<t<=t_0, then ||V(0,x)||>=b(x)"
            ),
        },
        {
            "id": "backward_heat_monomial",
            "kind": "exact_identity",
            "claim": (
                "The backward-heat evolution of an m-fold local zero is "
                "a finite polynomial with an exact coefficient formula."
            ),
            "formula": (
                "P_m(t,y)=exp(-t*D_y^2)y^m"
                "=m!*sum_(k=0)^floor(m/2)"
                "((-t)^k*y^(m-2k))/(k!*(m-2k)!)"
            ),
        },
        {
            "id": "hermite_rescaling",
            "kind": "exact_identity",
            "claim": (
                "The local heat polynomial is exactly a probabilists' "
                "Hermite polynomial on the parabolic scale."
            ),
            "formula": (
                "With s=sqrt(2t), z=y/s, "
                "P_m(t,y)=s^m*He_m(z)"
            ),
        },
        {
            "id": "real_simple_splitting",
            "kind": "exact_lemma",
            "claim": (
                "For every t>0 the m-fold endpoint zero splits into m "
                "simple real zeros, and P_m and its first derivative "
                "have no common zero."
            ),
            "formula": (
                "zeros(P_m(t,.))=sqrt(2t)*zeros(He_m); "
                "He_m'=m*He_(m-1)"
            ),
        },
        {
            "id": "adaptive_first_jet",
            "kind": "exact_identity",
            "claim": (
                "The time-dependent positive scaling balances both local "
                "jet components at the same parabolic order."
            ),
            "formula": (
                "W_m=(P_m,s*partial_y P_m)"
                "=s^m*(He_m(z),m*He_(m-1)(z)); s=sqrt(2t)"
            ),
        },
        {
            "id": "adaptive_heat_jet",
            "kind": "exact_identity",
            "claim": (
                "At fixed y, the adaptive jet derivative has an exact "
                "Hermite representation two parabolic orders lower."
            ),
            "formula": (
                "partial_t W_m=s^(m-2)*(A_m(z),B_m(z)); "
                "A_m=-m*(m-1)*He_(m-2); "
                "B_m=m*He_(m-1)-m*(m-1)*(m-2)*He_(m-3)"
            ),
        },
        {
            "id": "multiplicity_condition_number",
            "kind": "exact_inequality",
            "claim": (
                "For each fixed multiplicity, the adaptive heat-jet "
                "condition number is O(1/t), with a finite explicit "
                "one-variable Hermite supremum."
            ),
            "formula": (
                "C_m=sup_(z in R) sqrt(A_m(z)^2+B_m(z)^2)"
                "/sqrt(He_m(z)^2+m^2*He_(m-1)(z)^2)<infinity; "
                "||partial_t W_m||/||W_m||<=C_m/(2t)"
            ),
        },
        {
            "id": "cofinal_relative_transport",
            "kind": "exact_inequality",
            "claim": (
                "The multiplicity-compatible relative transport cost on "
                "one cofinal successor tends to zero and preserves every "
                "positive-time jet even if its endpoint clearance vanishes."
            ),
            "formula": (
                "For K_m=C_m/2 and t_j=1/(5j), "
                "||W_m(t_(j+1),y)||"
                ">=(j/(j+1))^K_m*||W_m(t_j,y)||; "
                "K_m*log(1+1/j)->0"
            ),
        },
        {
            "id": "time_dependent_scaling_transfer",
            "kind": "exact_identity",
            "claim": (
                "For any heat solution H_t=-H_xx, the adaptive jet "
                "has an exact transport derivative and is related to "
                "the ordinary first jet by an orientation-preserving map."
            ),
            "formula": (
                "W=(H,sH_x), s=sqrt(2t); "
                "partial_t W=(-H_xx,H_x/s-s*H_xxx); "
                "diag(1,s) is in GL+(2,R) for t>0"
            ),
        },
        {
            "id": "xi_domain_split_target",
            "kind": "open_xi_target",
            "claim": (
                "Replace the unnecessarily fixed-x-uniform small-ball "
                "shortcut by a two-regime successor theorem: arithmetic "
                "entry/cone control for new high-frequency strips and "
                "multiplicity-compatible relative transport for their "
                "older descendants."
            ),
            "formula": (
                "entry strip: use a corrected Xi half-plane or slope-gap "
                "certificate; descendant collars: prove "
                "||partial_t(H,sH_x)||<=(K/t)||(H,sH_x)|| with s=sqrt(2t), "
                "or a rigorously comparable local bound"
            ),
            "status": "open and not proved",
        },
    ]


def render_note() -> str:
    return f"""# Newman Multiplicity-Compatible Adaptive-Jet Benchmark

Date: {DATE}

Status: exact hypothetical heat-flow benchmark, theorem-domain correction, and open Xi handoff; not a proof of an all-stage successor, `Lambda<=0`, or RH.

## Why The Uniform Floor Is Too Strong

Suppose a continuous first jet `V(t,x)` at one fixed `x` obeys

```text
||V(t,x)||>=b(x)>0 for every 0<t<=t_0.
```

Taking `t` down to zero gives `||V(0,x)||>=b(x)`. Thus a positive
fixed-`x` floor uniform through `0<tL<=c_*` excludes a multiple endpoint
zero. RH itself does not assert simplicity, and the cofinal successor does
not need such a floor. Its bottom times are `t_j=1/(5j)`, and an old cell
may have clearance tending to zero while remaining nonzero at every
positive stage.

The existing positive-boundary attainment and delta-localization gates
already establish this quantifier warning and the arbitrary-multiplicity
Hermite split. The new contribution below is to put that split into the
adiabatic successor coordinates and derive its relative transport cost.

## Multiple-Zero Thought Experiment

Let an endpoint profile have an `m`-fold local zero and evolve by the same
backward heat equation:

```text
P_m(t,y)=exp(-t D_y^2)y^m
        =m! sum_(k=0)^floor(m/2)
          (-t)^k y^(m-2k)/(k!(m-2k)!).             (1)
```

For `s=sqrt(2t)` and `z=y/s`, (1) is exactly

```text
P_m(t,y)=s^m He_m(z),                               (2)
```

where `He_m` is the probabilists' Hermite polynomial. Its zeros are real
and simple. Since `He_m'=m He_(m-1)`, consecutive Hermite polynomials have
no common zero. Therefore an `m`-fold zero at `t=0` becomes `m` distinct
real zeros for every `t>0`; it does not create a positive-time contact.

This is precisely the kind of mathematically legitimate limiting
thought experiment that a simplicity-assuming proof would miss.

## Adaptive Jet

Use the time-dependent positive scaling

```text
W_m(t,y)=(P_m(t,y),s partial_y P_m(t,y)),
s=sqrt(2t).
```

Then both components have the same parabolic order:

```text
W_m=s^m(He_m(z),m He_(m-1)(z)).                    (3)
```

At fixed `y`,

```text
partial_t W_m=s^(m-2)(A_m(z),B_m(z)),              (4)
A_m=-m(m-1)He_(m-2),
B_m=m He_(m-1)-m(m-1)(m-2)He_(m-3),
```

with negative-index Hermite terms interpreted as zero. Define

```text
C_m=sup_(z in R)
 sqrt(A_m(z)^2+B_m(z)^2)
 /sqrt(He_m(z)^2+m^2 He_(m-1)(z)^2).              (5)
```

The denominator in (5) never vanishes, and the quotient tends to zero as
`|z|` tends to infinity, so `C_m<infinity`. Equations (3)-(5) give

```text
||partial_t W_m||/||W_m||<=C_m/(2t).               (6)
```

The singularity in (6) is the correct one: it permits the endpoint
clearance to vanish as a power of `t` without permitting a contact at
positive `t`.

## Cofinal Cost

Put `K_m=C_m/2`. Gronwall between consecutive bottom times gives

```text
||W_m(t_(j+1),y)||
 >=(t_(j+1)/t_j)^K_m ||W_m(t_j,y)||
 =(j/(j+1))^K_m ||W_m(t_j,y)||.                    (7)
```

The logarithmic cost is

```text
K_m log(1+1/j)->0.                                 (8)
```

Every factor in (7) is positive, although their infinite product may tend
to zero. That is enough for the successor induction, whose stages only
need positive-time origin exclusion.

For a general heat solution `H_t=-H_xx`, the same adaptive jet
`W=(H,sH_x)` obeys

```text
partial_t W=(-H_xx,H_x/s-sH_xxx),  s=sqrt(2t).     (9)
```

For every `t>0`, `diag(1,s)` lies in `GL+(2,R)`, so this scaling preserves
the contact set, orientation, and first-jet degree.

## Corrected Xi Strategy

The exact benchmark suggests a two-regime theorem.

```text
entry regime:
  certify each newly added high-frequency right strip by a corrected
  Riemann-Siegel half-plane, slope-gap, or equivalent arithmetic theorem;

descendant regime:
  transport previously certified cells down their later collars using
  the adaptive jet and a relative K/t estimate, or a rigorously comparable
  local factorization estimate.
```

This would use the history carried by the successor induction instead of
reproving one absolute lower bound on the entire old edge at every stage.
It is compatible with finite-multiplicity zeros at `t=0`.

## Proof Boundary

The finite heat-polynomial formula, Hermite rescaling, simple real
splitting, adaptive-jet identities, `O(1/t)` condition number, cofinal
relative transport cost, and fixed-`x` uniform-floor guard are exact. They
are a benchmark, not a local factorization theorem for Xi. No uniform
adaptive relative bound, entry-strip cone, all-`j` composition,
`Lambda<=0`, RH, PF-infinity, or Clay-prize conclusion is proved.

Machine-audited files:

```text
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```
"""


def main() -> int:
    artifact = {
        "kind": STEM,
        "date": DATE,
        "status": (
            "exact multiplicity-compatible adaptive-jet thought "
            "experiment and theorem-domain correction with one open "
            "two-regime Xi successor handoff"
        ),
        "proof_boundary": (
            "The monomial heat evolution, Hermite rescaling, adaptive "
            "jet and heat derivative, finite O(1/t) condition number "
            "for each fixed multiplicity, cofinal relative transport, "
            "and fixed-x uniform-floor implication are exact. No Xi "
            "adaptive relative estimate or entry-strip theorem is proved. "
            "This does not prove an all-j successor, Lambda<=0, RH, "
            "PF-infinity, or the Clay prize."
        ),
        "source_sha256": {
            "adiabatic_successor_lemma": file_hash(SUCCESSOR_RESULT),
            "single_carrier_adiabatic_benchmark": file_hash(
                SINGLE_CARRIER_RESULT
            ),
            "crossing_slope_gap_reduction": file_hash(CROSSING_RESULT),
            "positive_boundary_attainment_lemma": file_hash(
                ATTAINMENT_RESULT
            ),
            "positive_boundary_delta_localization_gate": file_hash(
                DELTA_LOCALIZATION_RESULT
            ),
        },
        "builder_sha256": file_hash(Path(__file__)),
        "rows": rows(),
        "exact_summary": {
            "row_count": 10,
            "exact_benchmark_reductions": 9,
            "open_xi_targets": 1,
            "adaptive_condition_number": "C_m/(2t)",
            "cofinal_step_cost": "K_m*log(1+1/j)->0",
            "requires_endpoint_simplicity": False,
            "all_j_xi_theorem": False,
        },
        "next_exact_obligation": (
            "Generalize the phase-cell successor lemma to a continuously "
            "time-dependent positive jet scaling, then prove an Xi-specific "
            "relative adaptive-jet estimate on descendant cells and a "
            "separate arithmetic half-plane or slope-gap theorem on newly "
            "entering right strips."
        ),
    }
    RESULT_PATH.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    NOTE_PATH.write_text(render_note(), encoding="utf-8")
    print(
        "built Newman multiplicity-compatible adaptive-jet benchmark: "
        "10 rows, 9 exact reductions, 1 open two-regime Xi handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
