#!/usr/bin/env python3
"""Build the time-dependent scaled-jet successor lemma."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_time_dependent_scaled_successor_lemma"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{STEM}.md"
SUCCESSOR_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma.json"
)
CROSSING_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_crossing_slope_gap_reduction.json"
)
ADAPTIVE_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_multiplicity_compatible_adaptive_jet_benchmark.json"
)
DATE = "2026-07-25"


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def rows() -> list[dict]:
    return [
        {
            "id": "positive_scaled_jet",
            "kind": "exact_definition",
            "claim": (
                "For any continuous positive scale s(t,x), the scaled "
                "first jet has exactly the same contact set as the "
                "ordinary first jet."
            ),
            "formula": (
                "V_s=(H,s*H_x), s>0; "
                "V_s=0 iff (H,H_x)=0"
            ),
        },
        {
            "id": "scaled_degree_homotopy",
            "kind": "exact_lemma",
            "claim": (
                "Positive scaling preserves boundary degree and local "
                "contact index through an explicit GL+(2,R) homotopy."
            ),
            "formula": (
                "A_r=diag(1,(1-r)+r*s), 0<=r<=1; "
                "det(A_r)=(1-r)+r*s>0"
            ),
        },
        {
            "id": "time_dependent_heat_jet",
            "kind": "exact_identity",
            "claim": (
                "The time derivative of a time-dependent scaled jet has "
                "one explicit scale-variation term."
            ),
            "formula": (
                "If H_t=-H_xx, then "
                "partial_t V_s=(-H_xx,s_t*H_x-s*H_xxx)"
            ),
        },
        {
            "id": "additive_scaled_transport",
            "kind": "exact_inequality",
            "claim": (
                "The original convex-cell additive transport theorem "
                "remains valid with the full time-dependent derivative."
            ),
            "formula": (
                "If V_s(t_j,I) subset C, d=dist(0,C)>0, and "
                "M>=sup_collar||partial_t V_s||, then "
                "||V_s(t,x)||>=d-|t-t_j|*M"
            ),
        },
        {
            "id": "relative_scaled_transport",
            "kind": "exact_inequality",
            "claim": (
                "An integrable relative heat-jet bound transports every "
                "nonzero old jet multiplicatively and cannot create a "
                "positive-time contact."
            ),
            "formula": (
                "If ||partial_t V_s(t,x)||<=kappa(t)||V_s(t,x)||, then "
                "||V_s(t,x)||>=exp(-integral_t^t_j kappa(r)dr)"
                "*||V_s(t_j,x)||"
            ),
        },
        {
            "id": "relative_cell_clearance",
            "kind": "exact_inequality",
            "claim": (
                "A whole old phase cell inherits a positive descendant "
                "clearance without requiring an additive ratio below one."
            ),
            "formula": (
                "If dist(0,V_s(t_j,I))>=d_j>0, then "
                "dist(0,V_s(t,I))>=d_j*exp(-integral_t^t_j kappa)"
            ),
        },
        {
            "id": "adaptive_hermite_specialization",
            "kind": "exact_composition",
            "claim": (
                "The multiplicity-compatible scale s=sqrt(2t) inserts "
                "the exact Hermite condition number into the generalized "
                "successor."
            ),
            "formula": (
                "s=sqrt(2t), s_t=1/s; "
                "kappa_m(t)=C_m/(2t)"
            ),
        },
        {
            "id": "cofinal_multiplicative_factor",
            "kind": "exact_composition",
            "claim": (
                "On t_j=1/(5j), the relative descendant factor is "
                "positive at every stage and its one-step logarithmic "
                "cost tends to zero."
            ),
            "formula": (
                "exp(-integral_(t_(j+1))^t_j C_m/(2t)dt)"
                "=(j/(j+1))^(C_m/2)>0; "
                "(C_m/2)*log(1+1/j)->0"
            ),
        },
        {
            "id": "positive_scale_switch",
            "kind": "exact_lemma",
            "claim": (
                "Any two positive jet scales can be blended without "
                "changing contacts or degree, allowing separate entry "
                "and descendant coordinates."
            ),
            "formula": (
                "s_r=(1-r)*s_entry+r*s_descendant>0; "
                "diag(1,s_r) is in GL+(2,R)"
            ),
        },
        {
            "id": "two_regime_successor_target",
            "kind": "open_xi_target",
            "claim": (
                "A finite base, a new-strip entry cone, and a relative "
                "scaled-jet bound on every descendant collar would give "
                "the cofinal successor without an endpoint-simplicity "
                "assumption."
            ),
            "formula": (
                "entry strips: q_j.V_(s_entry)>=eta_j>0; "
                "descendant collars: "
                "||partial_t V_(s_desc)||<=kappa_j(t)||V_(s_desc)|| "
                "with finite collar integrals; compose through a positive "
                "scale homotopy"
            ),
            "status": "open and not proved",
        },
    ]


def render_note() -> str:
    return f"""# Newman Time-Dependent Scaled-Jet Successor Lemma

Date: {DATE}

Status: exact generalized successor and relative-transport lemma with an open two-regime Xi antecedent; not a proof of an all-stage successor, `Lambda<=0`, or RH.

## Positive Scaled Jet

Let `H_t=-H_xx` and choose any continuous scale

```text
s(t,x)>0,                    V_s=(H,sH_x).          (1)
```

Then `V_s=0` exactly when `(H,H_x)=0`. On any contact-free boundary,

```text
A_r=diag(1,(1-r)+r s),       0<=r<=1,               (2)
det(A_r)=(1-r)+r s>0
```

is an explicit homotopy from the ordinary first jet to (1). Hence positive
scaling preserves the contact set, local index, boundary winding, and degree.

## Exact Transport Derivative

Unlike the earlier scale `1/ell(x)`, the new scale may depend on time.
Direct differentiation gives

```text
partial_t V_s=(-H_xx,s_t H_x-sH_xxx).              (3)
```

Thus the original additive phase-cell theorem remains valid verbatim after
its derivative bound is replaced by the norm of (3). If

```text
V_s(t_j,I) subset C,          d=dist(0,C)>0,
M>=sup_collar ||partial_t V_s||,
```

then

```text
||V_s(t,x)||>=d-|t-t_j|M.                           (4)
```

## Relative Transport

The time-dependent coordinate also exposes a stronger interface. Suppose
on one descendant collar

```text
||partial_t V_s(t,x)||<=kappa(t)||V_s(t,x)||,       (5)
```

where `kappa` is integrable over that positive-time collar. Gronwall in
reverse time gives

```text
||V_s(t,x)||
 >=exp(-integral_t^t_j kappa(r)dr)||V_s(t_j,x)||.   (6)
```

Consequently an old cell with clearance `d_j` has descendant clearance

```text
d_j exp(-integral_t^t_j kappa)>0.                  (7)
```

No additive ratio below one and no floor uniform as `t` tends to zero is
needed. A finite relative integral on each individual collar is enough.

## Hermite Specialization

For the multiplicity-compatible scale

```text
s=sqrt(2t),                   s_t=1/s,              (8)
```

the exact Hermite benchmark supplies

```text
kappa_m(t)=C_m/(2t).                                 (9)
```

On the cofinal times `t_j=1/(5j)`, (6) becomes

```text
||V_s(t_(j+1),x)||
 >=(j/(j+1))^(C_m/2)||V_s(t_j,x)||,                (10)
```

and the one-step logarithmic cost
`(C_m/2)log(1+1/j)` tends to zero. Every factor in (10) is positive,
although the endpoint limit may have zero clearance.

## Scale Interface

Entry-strip arithmetic may be cleaner in a frequency scale such as `1/L`,
whereas old-cell multiplicity transport is cleaner in `sqrt(2t)`. Any two
positive scales have the interpolation

```text
s_r=(1-r)s_entry+r s_descendant>0.                 (11)
```

Equation (2) applied to (11) proves that switching coordinates cannot create
a contact or alter degree. A rigorous implementation still has to bound the
transport derivative through a chosen overlap; the topology itself is exact.

## Conditional Successor

For

```text
P_j=[t_j,1/5]x[0,R_j],
C_j=[t_(j+1),t_j]x[0,R_j],
S_j=[t_(j+1),1/5]x[R_j,R_(j+1)],
```

a multiplicity-compatible all-stage theorem would follow from:

```text
1. one certified finite base P_J;
2. a strict first-jet half-plane cone on each new strip S_j in a positive
   entry scale;
3. an integrable relative bound (5) on every old-cell descendant collar,
   possibly after a positive scale switch;
4. explicit finite interfaces and source/remainder bounds.
```

The ordinary successor proof then applies to the scaled boundary paths,
and (2) transfers the resulting zero degree back to the ordinary first jet.

## Proof Boundary

The contact equivalence, positive-scaling degree homotopy, derivative (3),
additive and relative transport inequalities, cell-clearance factor,
Hermite specialization, cofinal factor, and positive scale switch are
exact. The Xi entry cone, relative descendant estimate, overlap bounds,
finite transition, and all-`j` composition are open. This lemma proves no
Q209 theorem, `Lambda<=0`, RH, PF-infinity, or Clay-prize conclusion.

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
            "exact time-dependent scaled-jet successor and relative "
            "transport lemma with one open two-regime Xi antecedent"
        ),
        "proof_boundary": (
            "Positive-scale contact/degree invariance, the full scaled "
            "heat derivative, additive and relative transport, Hermite "
            "specialization, cofinal factor, and scale-switch homotopy "
            "are exact. The Xi entry cone, relative descendant bound, "
            "interface estimates, and all-j composition are not proved. "
            "This does not prove Q209, Lambda<=0, RH, PF-infinity, or "
            "the Clay prize."
        ),
        "source_sha256": {
            "adiabatic_successor_lemma": file_hash(SUCCESSOR_RESULT),
            "crossing_slope_gap_reduction": file_hash(CROSSING_RESULT),
            "adaptive_jet_benchmark": file_hash(ADAPTIVE_RESULT),
        },
        "builder_sha256": file_hash(Path(__file__)),
        "rows": rows(),
        "exact_summary": {
            "row_count": 10,
            "exact_reductions": 9,
            "open_xi_targets": 1,
            "time_dependent_positive_scaling": True,
            "relative_transport_available": True,
            "requires_endpoint_simplicity": False,
            "all_j_xi_theorem": False,
        },
        "next_exact_obligation": (
            "Choose an explicit positive entry/descendant scale profile "
            "and prove Xi-specific relative adaptive-jet bounds on old "
            "cells together with a corrected Riemann-Siegel half-plane "
            "or weighted slope-gap theorem on each new right strip."
        ),
    }
    RESULT_PATH.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    NOTE_PATH.write_text(render_note(), encoding="utf-8")
    print(
        "built Newman time-dependent scaled-jet successor lemma: "
        "10 rows, 9 exact reductions, 1 open two-regime Xi antecedent"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
