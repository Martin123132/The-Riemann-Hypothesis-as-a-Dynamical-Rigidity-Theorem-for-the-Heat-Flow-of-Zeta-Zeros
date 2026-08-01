#!/usr/bin/env python3
"""Build the ray-aligned exhaustion and parabolic-frequency scale reduction."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_ray_aligned_parabolic_frequency_reduction"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
DATE = "2026-07-26"

C_STAR = Fraction(4_911_678_521, 1_933_561_194)
EXPLICIT_RAY = Fraction(25)
SCHEDULE_OFFSET = 100
TOP_TIME = Fraction(1, 2)
BASE_TIME = Fraction(1, 4)

SOURCE_FILES = {
    "dominant_saddle_global_ray": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "dominant_saddle_global_ray_certificate.json"
    ),
    "oscillatory_zeta_handoff": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "oscillatory_zeta_handoff_theorem.json"
    ),
    "positive_boundary_attainment": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_positive_boundary_attainment_lemma.json"
    ),
    "diagonal_exhaustion": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_positive_boundary_"
        "diagonal_exhaustion_gate.json"
    ),
    "compact_transversality": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_theta_"
        "compact_transversality_interval_certificate.json"
    ),
    "adaptive_jet": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_multiplicity_"
        "compatible_adaptive_jet_benchmark.json"
    ),
    "scaled_successor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_time_dependent_"
        "scaled_successor_lemma.json"
    ),
}


@dataclass(frozen=True)
class ReductionRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def source_hashes() -> dict[str, str]:
    missing = [str(path) for path in SOURCE_FILES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing source artifacts: {missing}")
    return {name: file_hash(path) for name, path in SOURCE_FILES.items()}


def build_exact() -> dict[str, str]:
    t, ell = sp.symbols("t ell", positive=True)
    q = 2 * t * ell**2
    scale = sp.sqrt(2 * t / (1 + q))
    parabolic = sp.sqrt(2 * t)
    frequency = 1 / ell
    alpha = sp.simplify(scale / parabolic)
    beta = sp.simplify(scale / frequency)

    if sp.simplify(scale ** -2 - (1 / (2 * t) + ell**2)) != 0:
        raise RuntimeError("reciprocal-square scale identity failed")
    if sp.simplify(alpha**2 + beta**2 - 1) != 0:
        raise RuntimeError("chart partition identity failed")
    if sp.simplify(sp.diff(sp.log(scale), t) - 1 / (2 * t * (1 + q))) != 0:
        raise RuntimeError("scale logarithmic time derivative failed")
    if sp.simplify(sp.diff(sp.log(alpha), t) + ell**2 / (1 + q)) != 0:
        raise RuntimeError("parabolic chart derivative failed")
    if sp.simplify(sp.diff(sp.log(beta), t) - 1 / (2 * t * (1 + q))) != 0:
        raise RuntimeError("frequency chart derivative failed")

    x = sp.symbols("x", positive=True)
    ell_x = sp.log(x / (4 * sp.pi))
    scale_x = sp.sqrt(2 * t / (1 + 2 * t * ell_x**2))
    expected_x = -2 * t * ell_x / (x * (1 + 2 * t * ell_x**2))
    if sp.simplify(sp.diff(sp.log(scale_x), x) - expected_x) != 0:
        raise RuntimeError("scale logarithmic spatial derivative failed")

    return {
        "scale": "s_pf(t,x)=(L(x)^2+(2t)^(-1))^(-1/2)",
        "equivalent_scale": "s_pf=sqrt(2t)/sqrt(1+2tL^2)",
        "dimensionless_parameter": "q=2tL^2",
        "parabolic_chart": "alpha=s_pf/sqrt(2t)=1/sqrt(1+q)",
        "frequency_chart": "beta=L*s_pf=sqrt(q/(1+q))",
        "partition": "alpha^2+beta^2=1",
        "time_derivative": "partial_t log(s_pf)=1/(2t(1+q))",
        "parabolic_derivative": "partial_t log(alpha)=-L^2/(1+q)",
        "frequency_derivative": "partial_t log(beta)=1/(2t(1+q))",
        "space_derivative": (
            "partial_x log(s_pf)=-2tL/(x(1+2tL^2)), "
            "L=log(x/(4pi))"
        ),
    }


def schedule_samples() -> list[dict[str, str | int]]:
    rows: list[dict[str, str | int]] = []
    c = EXPLICIT_RAY
    a = SCHEDULE_OFFSET
    for stage in (0, 1, 10, 100, 1000, 10_000):
        time = c / (a + stage)
        entry_log = Fraction(a + stage + 1)
        next_time = c / (a + stage + 1)
        entry_q = 2 * next_time * entry_log**2
        beta_squared = entry_q / (1 + entry_q)
        one_step_log_cost = math.log((a + stage + 1) / (a + stage))
        rows.append(
            {
                "stage": stage,
                "time": str(time),
                "entry_log_radius": str(entry_log),
                "next_time": str(next_time),
                "new_strip_min_tL": str(next_time * entry_log),
                "entry_q": str(entry_q),
                "entry_beta_squared": (
                    f"{float(beta_squared):.17g}"
                ),
                "unit_K_log_cost": f"{one_step_log_cost:.17g}",
            }
        )
    return rows


def reduction_rows() -> list[ReductionRow]:
    return [
        ReductionRow(
            id="rapf_01_explicit_schedule",
            role="exact_exhaustion",
            readiness="proved_exact",
            claim=(
                "An exponential-radius exhaustion keeps every newly entering "
                "strip on the already proved explicit tL=25 ray."
            ),
            formula=(
                "c=25, a=100, t_j=c/(a+j), "
                "L_j=log(R_j/(4pi))=a+j+1, "
                "R_j=4pi*exp(L_j)"
            ),
            proof_boundary=(
                "This changes the exhaustion geometry; it does not control "
                "the old descendant collars."
            ),
        ),
        ReductionRow(
            id="rapf_02_new_strip_identity",
            role="exact_exhaustion",
            readiness="proved_exact",
            claim=(
                "The minimum scaled time on the new strip from R_j to "
                "R_(j+1) is exactly 25."
            ),
            formula=(
                "t_(j+1)*L(x)>=t_(j+1)*L_j="
                "[25/(a+j+1)]*(a+j+1)=25"
            ),
            proof_boundary=(
                "Uses x>=R_j and monotonicity of L(x)."
            ),
        ),
        ReductionRow(
            id="rapf_03_outer_ray_composition",
            role="theorem_composition",
            readiness="available_exact",
            claim=(
                "The existing dominant-saddle theorem excludes every "
                "first-jet contact on every new strip of the explicit "
                "schedule."
            ),
            formula=(
                "L>=101>50 and tL>=25 imply "
                "H_x^2-H*H_xx>0"
            ),
            proof_boundary=(
                "Imports the checked dominant-saddle global-ray theorem. "
                "No claim is made below tL=25."
            ),
        ),
        ReductionRow(
            id="rapf_04_cofinal_geometry",
            role="exact_exhaustion",
            readiness="proved_exact",
            claim=(
                "The ray-aligned rectangles are cofinal in positive time "
                "and finite frequency."
            ),
            formula=(
                "t_j=25/(100+j)->0, "
                "R_j=4pi*exp(101+j)->infinity"
            ),
            proof_boundary=(
                "Cofinality alone does not certify the rectangles."
            ),
        ),
        ReductionRow(
            id="rapf_05_descendant_localization",
            role="exact_reduction",
            readiness="proved_exact",
            claim=(
                "After the schedule change, the only all-stage successor "
                "antecedent is transport on the old collar."
            ),
            formula=(
                "C_j^out=[t_(j+1),t_j]x[38,R_j]; "
                "S_j=[t_(j+1),1/2]x[R_j,R_(j+1)] "
                "is already contact-free; x<=38 is independently closed"
            ),
            proof_boundary=(
                "A relative old-collar estimate is still open and is the "
                "remaining infinite-dimensional Xi input."
            ),
        ),
        ReductionRow(
            id="rapf_06_blended_scale",
            role="exact_scale",
            readiness="proved_exact",
            claim=(
                "The reciprocal-square blend is a smooth positive scale "
                "joining the parabolic and frequency coordinates."
            ),
            formula=(
                "s_pf=(L^2+(2t)^(-1))^(-1/2)"
                "=sqrt(2t)/sqrt(1+2tL^2)"
            ),
            proof_boundary=(
                "Positive scaling preserves contacts and degree but does "
                "not itself prove a relative Xi estimate."
            ),
        ),
        ReductionRow(
            id="rapf_07_chart_partition",
            role="exact_scale",
            readiness="proved_exact",
            claim=(
                "The two chart weights form an exact unit partition and "
                "the blended scale stays within sqrt(2) of the smaller "
                "natural scale."
            ),
            formula=(
                "alpha=s_pf/sqrt(2t), beta=L*s_pf, "
                "alpha^2+beta^2=1; "
                "min(sqrt(2t),1/L)/sqrt(2)<=s_pf"
                "<=min(sqrt(2t),1/L)"
            ),
            proof_boundary=(
                "The comparison is conditioning control, not a lower bound "
                "for the Xi first jet."
            ),
        ),
        ReductionRow(
            id="rapf_08_scale_derivatives",
            role="exact_scale",
            readiness="proved_exact",
            claim=(
                "The scale-variation terms are explicit and locally "
                "integrable at every positive time."
            ),
            formula=(
                "partial_t log(s_pf)=1/[2t(1+2tL^2)]; "
                "partial_x log(s_pf)=-2tL/[x(1+2tL^2)]"
            ),
            proof_boundary=(
                "The time derivative is singular at t=0, as required by "
                "the arbitrary-multiplicity Hermite model."
            ),
        ),
        ReductionRow(
            id="rapf_09_two_chart_relative_transfer",
            role="exact_reduction",
            readiness="proved_exact",
            claim=(
                "A relative bound in either natural chart transfers to the "
                "blended chart with at most sqrt(2) conditioning loss in "
                "that chart's regime."
            ),
            formula=(
                "q<=1: kappa_pf<=sqrt(2)kappa_par+L^2/(1+q); "
                "q>=1: kappa_pf<=sqrt(2)kappa_freq"
                "+1/[2t(1+q)]"
            ),
            proof_boundary=(
                "The inequalities are conditional on Xi relative bounds "
                "kappa_par or kappa_freq; neither bound is asserted here."
            ),
        ),
        ReductionRow(
            id="rapf_10_entry_conditioning",
            role="exact_scale",
            readiness="proved_exact",
            claim=(
                "Every new point enters deeply in the frequency chart, "
                "then moves continuously toward the parabolic chart as "
                "its fixed spatial coordinate descends toward t=0."
            ),
            formula=(
                "on entry q=2t_(j+1)L_j^2=50L_j>=5050; "
                "for fixed x, q=2tL(x)^2->0 as t->0"
            ),
            proof_boundary=(
                "Chart conditioning is exact; descendant noncontact still "
                "requires the relative Xi inequality."
            ),
        ),
        ReductionRow(
            id="rapf_11_one_step_cost",
            role="exact_asymptotic",
            readiness="proved_exact",
            claim=(
                "The ray-aligned old-collar logarithmic time ratio tends "
                "to zero, so a K/t descendant estimate has vanishing "
                "one-step cost."
            ),
            formula=(
                "log(t_j/t_(j+1))="
                "log((100+j+1)/(100+j))->0"
            ),
            proof_boundary=(
                "An infinite product may tend to zero; the induction needs "
                "only a positive factor at each finite stage."
            ),
        ),
        ReductionRow(
            id="rapf_12_asymptotic_ray_variant",
            role="theorem_composition",
            readiness="available_asymptotic",
            claim=(
                "The same schedule can be placed on any fixed "
                "c=c_*+epsilon ray above its existential threshold."
            ),
            formula=(
                "c_*=4911678521/1933561194; "
                "t_(j+1)L_j=c_*+epsilon"
            ),
            proof_boundary=(
                "The oscillatory-zeta theorem supplies an existential "
                "L_epsilon; the explicit c=25 schedule avoids hiding it."
            ),
        ),
        ReductionRow(
            id="rapf_13_conditional_rh_composition",
            role="conditional_theorem",
            readiness="open_antecedent",
            claim=(
                "The published upper-time input, explicit new-strip ray, "
                "and relative blended-jet transport on every old collar "
                "would certify a cofinal family and imply Lambda<=0. "
                "Equivalently, one relative theorem on the inner wedge "
                "0<t<=min(1/4,25/L(x)) is sufficient."
            ),
            formula=(
                "if ||partial_t(H,s_pf H_x)||"
                "<=kappa_j(t)||(H,s_pf H_x)|| on C_j^out "
                "with integral_(t_(j+1))^(t_j)kappa_j<infinity "
                "for every j, then every P_j is contact-free; equivalently "
                "tau(x)=min(1/4,25/L(x)) and an integrable relative bound "
                "on 0<t<=tau(x), x>=38, suffices"
            ),
            proof_boundary=(
                "The relative old-collar antecedent is unproved. Therefore "
                "this row is not RH or a Clay-prize proof."
            ),
        ),
        ReductionRow(
            id="rapf_14_open_xi_descendant_theorem",
            role="theorem_target",
            readiness="open",
            claim=(
                "Construct an a priori two-chart relative majorant for the "
                "exact Xi heat flow on the ray-aligned descendant collars, "
                "with a finite bounded-L transition and without dividing by "
                "the unknown first-jet norm."
            ),
            formula=(
                "q<=1: establish kappa_par=O(1/t) in the Hermite chart; "
                "q>=1 below the entry ray: establish an integrable "
                "frequency-chart bound; kappa must not be defined as "
                "||partial_t V_pf||/||V_pf||"
            ),
            proof_boundary=(
                "On each compact positive-time interval, mere existence of "
                "an integrable relative bound is equivalent to nonvanishing: "
                "the reverse implication takes "
                "kappa=||partial_t V_pf||/||V_pf||. The required a priori "
                "majorant is the remaining RH-level analytic obligation."
            ),
        ),
    ]


def render_note(artifact: dict) -> str:
    samples = artifact["schedule_samples"]
    sample_lines = [
        "| stage | t_j | L_j | t_(j+1) | min tL | entry q | beta^2 | unit-K log cost |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in samples:
        sample_lines.append(
            "| {stage} | {time} | {entry_log_radius} | {next_time} | "
            "{new_strip_min_tL} | {entry_q} | {entry_beta_squared} | "
            "{unit_K_log_cost} |".format(**row)
        )
    table = "\n".join(sample_lines)
    return f"""# Newman Ray-Aligned Parabolic-Frequency Reduction

Date: {DATE}

Status: exact exhaustion redesign, exact positive scale, and conditional
one-antecedent cofinal reduction. The Xi descendant estimate remains open;
this is not a proof of `Lambda<=0` or RH.

## Why Change The Exhaustion

The linear-radius schedule used previously,

```text
t_j=1/(5j), R_j=j+38,
```

forces each new right strip into `tL->0`, even though the existing
high-frequency theorem already proves strict first-Laguerre positivity on
`L>=50`, `tL>=25`. The diagonal-exhaustion theorem permits any time floors
tending to zero and any radii tending to infinity. They need not be coupled
linearly.

Choose instead

```text
c=25, a=100,
t_j=c/(a+j),
L_j=log(R_j/(4pi))=a+j+1,
R_j=4pi exp(L_j).                                  (1)
```

The base time is `t_0=1/4>1/5`, so the published `Lambda<=1/5` bound makes
the base rectangle `[1/4,1/2]x[0,R_0]` contact-free.

## Every New Strip Is Already Closed

Let

```text
P_j=[t_j,1/2]x[0,R_j],
C_j^out=[t_(j+1),t_j]x[38,R_j],
S_j=[t_(j+1),1/2]x[R_j,R_(j+1)].
```

For every point of `S_j`,

```text
tL(x)>=t_(j+1)L_j
      =[25/(a+j+1)](a+j+1)=25.                    (2)
```

Also `L(x)>=L_j>=101>50`. The checked dominant-saddle global-ray theorem
therefore gives

```text
H_x^2-H H_xx>0                                    (3)
```

throughout every new strip. Hence a new strip cannot contain a first-jet
contact. This removes the formerly open all-stage right-strip cone or
weighted slope-gap theorem from the induction antecedents.

The schedule is cofinal:

```text
t_j->0, R_j->infinity.                             (4)
```

The independently certified compact theorem handles `0<=x<=38` for all
relevant times. What remains is exactly the outer old descendant collar
`C_j^out`.

## Parabolic-Frequency Scale

For `L=log(x/(4pi))>0`, define

```text
s_pf(t,x)=(L^2+(2t)^(-1))^(-1/2)
         =sqrt(2t)/sqrt(1+2tL^2),                  (5)
V_pf=(H,s_pf H_x).
```

Put `q=2tL^2` and compare (5) with the parabolic scale
`s_par=sqrt(2t)` and frequency scale `s_freq=1/L`:

```text
alpha=s_pf/s_par=1/sqrt(1+q),
beta=s_pf/s_freq=sqrt(q/(1+q)),
alpha^2+beta^2=1.                                 (6)
```

Consequently

```text
min(s_par,s_freq)/sqrt(2)<=s_pf<=min(s_par,s_freq). (7)
```

The transition `q=1` has `alpha=beta=1/sqrt(2)`. There is no badly
conditioned chart interface.

The exact scale derivatives are

```text
partial_t log(s_pf)=1/[2t(1+q)],
partial_t log(alpha)=-L^2/(1+q),
partial_t log(beta)=1/[2t(1+q)],
partial_x log(s_pf)=-2tL/[x(1+q)].                 (8)
```

At entry, (1)-(2) give

```text
q_entry=2t_(j+1)L_j^2=50L_j>=5050,                (9)
```

so the scale is almost exactly `1/L`. For any fixed descendant coordinate
`x`, `q=2tL(x)^2` tends to zero with `t`, and the same scale becomes almost
exactly `sqrt(2t)`, the multiplicity-compatible Hermite coordinate.

## Relative Two-Chart Handoff

Let `V_par=(H,s_par H_x)` and `V_freq=(H,s_freq H_x)`. Since
`V_pf=diag(1,alpha)V_par=diag(1,beta)V_freq`, the exact rescaling inequality
gives

```text
q<=1:
  kappa_pf<=sqrt(2) kappa_par+L^2/(1+q);

q>=1:
  kappa_pf<=sqrt(2) kappa_freq+1/[2t(1+q)].        (10)
```

Here `kappa_par` or `kappa_freq` denotes a proved relative bound for the
corresponding Xi jet. Equation (10) does not assert either bound. It shows
that once such bounds are established in their natural regimes, the smooth
interface costs at most `sqrt(2)`.

The schedule also has

```text
log(t_j/t_(j+1))
 =log((a+j+1)/(a+j))->0,                           (11)
```

so a Hermite-type `K/t` estimate has vanishing one-step logarithmic cost.

## Conditional Cofinal Theorem

Suppose that on every outer old collar `C_j^out` there is an integrable
`kappa_j(t)` such that

```text
||partial_t V_pf(t,x)||<=kappa_j(t)||V_pf(t,x)||.  (12)
```

Gronwall transports every old nonzero jet through `C_j^out`. Equation (3)
certifies `S_j`, and the positive scale preserves contact sets and degree.
Induction then certifies every `P_j`; cofinality (4) gives positive-time
simplicity and therefore `Lambda<=0`.

This composition is exact, but (12) for Xi is not proved.

The rectangular condition has an equivalent direct wedge formulation. Put

```text
tau(x)=min(1/4,25/L(x)),             x>=38.         (13)
```

At `t=tau(x)`, the first jet is nonzero: when `L(x)<=100`, this follows
from `1/4>Lambda`; when `L(x)>=100`, it follows from `tL=25` and the
dominant-saddle theorem. Therefore it is sufficient to prove, for every
fixed `x>=38`,

```text
||partial_t V_pf(t,x)||<=kappa_x(t)||V_pf(t,x)||,
0<t<=tau(x),    integral_delta^tau(x) kappa_x<infinity
for every delta>0.                                  (14)
```

Gronwall then excludes a contact at every positive `t`. The exponential
rectangles (1) are a finite-stage exhaustion of this one wedge target.

There is an essential nonpromotion guard. For any `C^1` vector path `V` on
`[delta,tau]`, a nonzero endpoint together with an integrable relative bound
implies nonvanishing by Gronwall. Conversely, if `V` is already nonvanishing,
then

```text
kappa(t)=||partial_t V(t)||/||V(t)||
```

is continuous and supplies such a bound. Thus bare existence of `kappa_x`
is equivalent to the desired positive-time noncontact statement. A useful
Xi theorem must construct an a priori majorant from independently bounded
arithmetic or analytic quantities, without dividing by the unknown jet norm.
The scale `s_pf` conditions the two asymptotic charts; it does not itself
provide that majorant.

## Explicit Schedule Diagnostics

{table}

## Sharper Asymptotic Variant

For any fixed `epsilon>0`, the same construction may replace `25` by
`c_*+epsilon`, where

```text
c_*={C_STAR.numerator}/{C_STAR.denominator}.
```

The oscillatory-zeta theorem then closes every new strip once its
existential `L_epsilon` threshold is passed. The `c=25`, `L>=50`
construction above is preferred as the fully explicit reduction.

## Proof Boundary

The exhaustion identities, new-strip composition, blended scale,
two-chart conditioning, derivative formulas, and conditional induction are
rigorous. The result does not prove the relative Xi bound (12), contact-free
old collars at every stage, `Lambda<=0`, RH, PF-infinity, or the Clay prize.

Machine-audited files:

```text
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```
"""


def build_artifact() -> dict:
    rows = reduction_rows()
    return {
        "kind": STEM,
        "date": DATE,
        "status": (
            "exact ray-aligned exhaustion and parabolic-frequency scale "
            "reduction with one open Xi descendant antecedent"
        ),
        "proof_boundary": (
            "The artifact proves exact exhaustion, scale, chart, and "
            "conditional composition identities. It does not prove the Xi "
            "relative old-collar estimate, Lambda<=0, RH, PF-infinity, or "
            "the Clay prize."
        ),
        "constants": {
            "c_star": str(C_STAR),
            "explicit_ray": str(EXPLICIT_RAY),
            "schedule_offset": SCHEDULE_OFFSET,
            "top_time": str(TOP_TIME),
            "base_time": str(BASE_TIME),
            "minimum_log_radius": SCHEDULE_OFFSET + 1,
        },
        "exact": build_exact(),
        "schedule_samples": schedule_samples(),
        "source_sha256": source_hashes(),
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_reductions": sum(
                row.readiness in {"proved_exact", "available_exact"}
                for row in rows
            ),
            "asymptotic_compositions": sum(
                row.readiness == "available_asymptotic" for row in rows
            ),
            "open_antecedents": sum(
                row.readiness in {"open", "open_antecedent"} for row in rows
            ),
            "new_strip_open_antecedents": 0,
            "old_collar_open_antecedents": 1,
        },
        "builder_sha256": file_hash(Path(__file__)),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    summary = artifact["summary"]
    print(
        "built Newman ray-aligned parabolic-frequency reduction: "
        f"{summary['rows']} rows, "
        f"{summary['exact_reductions']} exact reductions, "
        f"{summary['asymptotic_compositions']} asymptotic composition, "
        f"{summary['open_antecedents']} open bookkeeping rows for "
        "1 Xi descendant theorem"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
