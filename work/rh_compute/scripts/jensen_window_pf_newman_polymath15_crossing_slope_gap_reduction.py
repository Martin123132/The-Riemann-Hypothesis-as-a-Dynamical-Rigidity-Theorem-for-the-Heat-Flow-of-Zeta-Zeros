#!/usr/bin/env python3
"""Build the crossing-restricted Gram and slope-gap reduction."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import mpmath as mp

import jensen_window_pf_newman_polymath15_critical_component_wronskian_gate as component_gate
import jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction as phase_gate


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_polymath15_crossing_slope_gap_reduction"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{STEM}.md"
COMPONENT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_critical_component_wronskian_gate.json"
)
PHASE_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction.json"
)
SINGLE_CARRIER_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_single_carrier_adiabatic_benchmark.json"
)
DATE = "2026-07-25"
COARSE_DPS = 55
FINE_DPS = 75


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def rows() -> list[dict]:
    return [
        {
            "id": "crossing_real_component_jet",
            "kind": "exact_identity",
            "claim": (
                "For e_j=a_j*exp(i*theta_j) with "
                "e_j'=(u_j+i*v_j)e_j, the two real observables are "
                "linear forms in the positive amplitude vector."
            ),
            "formula": (
                "p_j=cos(theta_j); "
                "q_j=u_j*cos(theta_j)-v_j*sin(theta_j); "
                "X=sum_j a_j*p_j; U=sum_j a_j*q_j"
            ),
        },
        {
            "id": "scaled_gram_expansion",
            "kind": "exact_identity",
            "claim": (
                "The corrected scaled first-jet square is a rank-at-most-two "
                "Gram form with an explicit diagonal/off-diagonal expansion."
            ),
            "formula": (
                "Q_L=X^2+(U/L)^2=a^T*(p*p^T+q*q^T/L^2)*a"
                "=sum_j a_j^2*(p_j^2+q_j^2/L^2)"
                "+2*sum_(j<k)a_j*a_k*(p_j*p_k+q_j*q_k/L^2)"
            ),
        },
        {
            "id": "crossing_rank_collapse",
            "kind": "exact_impossibility",
            "claim": (
                "On the crossing hyperplane p.a=0, Q_L=(q.a)^2/L^2 "
                "has rank at most one. For m>=3 its kernel inside that "
                "hyperplane has dimension at least m-2."
            ),
            "formula": (
                "ker(Q_L|_(p_perp)) contains p_perp intersect q_perp; "
                "there is no positive uniform lower eigenvalue on all "
                "crossing amplitude vectors"
            ),
        },
        {
            "id": "near_crossing_slope_identity",
            "kind": "exact_identity",
            "claim": (
                "Splitting c_j=a_j*p_j by sign converts U into two "
                "weighted slope means plus an explicit zero-real-part term."
            ),
            "formula": (
                "h_j=q_j/p_j=d_j/c_j; "
                "M_+=sum_(c_j>0)c_j, M_-=sum_(c_j<0)(-c_j), "
                "M=(M_++M_-)/2, X=M_+-M_-; "
                "U=M*(h_+-h_-)+(X/2)*(h_++h_-)+D_0"
            ),
        },
        {
            "id": "exact_crossing_slope_gap",
            "kind": "exact_identity",
            "claim": (
                "At X=0 the positive and negative real masses agree, and "
                "the complete derivative observable is their common mass "
                "times a weighted slope gap, plus D_0."
            ),
            "formula": (
                "M_+=M_-=M; U=M*(h_+-h_-)+D_0; "
                "Q_L=(M*(h_+-h_-)+D_0)^2/L^2"
            ),
        },
        {
            "id": "strong_slope_separation_floor",
            "kind": "exact_conditional_inequality",
            "claim": (
                "A componentwise ordering of the crossing slopes is a "
                "sufficient, but deliberately stronger than necessary, "
                "route to a quantitative first-jet floor."
            ),
            "formula": (
                "If D_0=0 and either inf_P h_j-sup_N h_j>=gamma>0 "
                "or inf_N h_j-sup_P h_j>=gamma>0, then "
                "|U|>=M*gamma and Q_L>=M^2*gamma^2/L^2"
            ),
        },
        {
            "id": "rotating_half_plane_floor",
            "kind": "exact_conditional_inequality",
            "claim": (
                "A strict common half-plane for the component first-jet "
                "vectors gives a global floor and excludes a positive-weight "
                "double crossing."
            ),
            "formula": (
                "g_j=(p_j,q_j/L); if a unit eta and kappa>0 obey "
                "eta.g_j>=kappa for every j, then "
                "sqrt(Q_L)>=kappa*sum_j a_j"
            ),
        },
        {
            "id": "centered_log_moment_frame",
            "kind": "exact_identity",
            "claim": (
                "At t=0 the Dirichlet carriers admit an arbitrary centered "
                "log-frequency frame, isolating the endpoint defect."
            ),
            "formula": (
                "v_mu=beta_0'+mu/2; "
                "E'=i*v_mu*E+(i/2)*sum_(n=1)^N(log(n)-mu)e_n"
                "+R_(0,mu), R_(0,mu)=e_0'-i*v_mu*e_0; "
                "at X=0, U=-v_mu*Y-(1/2)Im(sum(log(n)-mu)e_n)"
                "+Re(R_(0,mu))"
            ),
        },
        {
            "id": "nonzero_cosine_four_carrier_obstruction",
            "kind": "exact_countermodel",
            "claim": (
                "Positive amplitudes, strictly ordered negative speeds, "
                "and nonzero real part in every component still allow the "
                "two crossing slope means to agree exactly."
            ),
            "formula": (
                "E_*(x)=3*exp(i*pi/3-4ix)+3*exp(2i*pi/3-3ix)"
                "+7*exp(-i*pi/3-2ix)+7*exp(4i*pi/3-ix); "
                "E_*(0)=-4i*sqrt(3), E_*'(0)=-5i, "
                "Re(E_*''(0))=-21, M=5, h_+=h_-=-sqrt(3)/5"
            ),
        },
        {
            "id": "xi_weighted_slope_gap_target",
            "kind": "open_xi_target",
            "claim": (
                "At an exact crossing, the currently recorded uniform Xi "
                "small-ball shortcut is exactly a weighted-slope "
                "separation, not a generic spectral-coercivity assertion. "
                "This fixed-x-uniform shortcut is stronger than the "
                "multiplicity-compatible cofinal theorem logically needed."
            ),
            "formula": (
                "|M*(h_+-h_-)+D_0|"
                ">2000*sqrt(2)*L*exp(-3*L/4) for L>=50; "
                "equivalently Q_L>8000000*exp(-3*L/2) when X=0"
            ),
            "status": "open and not proved",
        },
    ]


def component_row(
    root: mp.mpf,
    time: mp.mpf,
    cutoff: int,
    label: str,
) -> dict:
    components = component_gate.corrected_components(root, time, cutoff)
    derivatives = [
        mp.diff(
            lambda value, index=index: component_gate.corrected_components(
                value, time, cutoff
            )[index],
            root,
        )
        for index in range(len(components))
    ]
    amplitudes = [abs(component) for component in components]
    real_parts = [mp.re(component) for component in components]
    real_derivatives = [mp.re(derivative) for derivative in derivatives]
    x_value = mp.fsum(real_parts)
    u_value = mp.fsum(real_derivatives)
    ell = mp.log(root / (4 * mp.pi))

    positive = [index for index, value in enumerate(real_parts) if value > 0]
    negative = [index for index, value in enumerate(real_parts) if value < 0]
    zero = [index for index, value in enumerate(real_parts) if value == 0]
    if not positive or not negative:
        raise RuntimeError(f"{label}: crossing sign partition degenerated")

    mass_positive = mp.fsum(real_parts[index] for index in positive)
    mass_negative = mp.fsum(-real_parts[index] for index in negative)
    mean_mass = (mass_positive + mass_negative) / 2
    positive_slopes = [
        real_derivatives[index] / real_parts[index] for index in positive
    ]
    negative_slopes = [
        real_derivatives[index] / real_parts[index] for index in negative
    ]
    mean_positive = (
        mp.fsum(real_derivatives[index] for index in positive)
        / mass_positive
    )
    mean_negative = (
        -mp.fsum(real_derivatives[index] for index in negative)
        / mass_negative
    )
    zero_derivative = mp.fsum(real_derivatives[index] for index in zero)
    reconstructed = (
        mean_mass * (mean_positive - mean_negative)
        + x_value * (mean_positive + mean_negative) / 2
        + zero_derivative
    )
    identity_error = abs(reconstructed - u_value) / max(1, abs(u_value))

    forward_separation = min(positive_slopes) - max(negative_slopes)
    reverse_separation = min(negative_slopes) - max(positive_slopes)

    angles = sorted(
        (mp.atan2(
            real_derivatives[index] / (amplitudes[index] * ell),
            real_parts[index] / amplitudes[index],
        ) % (2 * mp.pi))
        for index in range(len(components))
    )
    circular_gaps = [
        angles[index + 1] - angles[index]
        for index in range(len(angles) - 1)
    ]
    circular_gaps.append(angles[0] + 2 * mp.pi - angles[-1])
    covering_arc = 2 * mp.pi - max(circular_gaps)

    component_norm_sum = mp.fsum(
        mp.sqrt(
            real_parts[index] ** 2
            + (real_derivatives[index] / ell) ** 2
        )
        for index in range(len(components))
    )
    aggregate_norm = mp.sqrt(x_value**2 + (u_value / ell) ** 2)
    cancellation_factor = component_norm_sum / aggregate_norm

    return {
        "label": label,
        "time": mp.nstr(time, 20, strip_zeros=False),
        "cutoff": cutoff,
        "component_count": len(components),
        "root": mp.nstr(root, 40, strip_zeros=False),
        "ell": mp.nstr(ell, 30, strip_zeros=False),
        "abs_X": mp.nstr(abs(x_value), 25, strip_zeros=False),
        "U": mp.nstr(u_value, 35, strip_zeros=False),
        "positive_component_count": len(positive),
        "negative_component_count": len(negative),
        "zero_real_component_count": len(zero),
        "positive_real_mass": mp.nstr(
            mass_positive, 30, strip_zeros=False
        ),
        "negative_real_mass": mp.nstr(
            mass_negative, 30, strip_zeros=False
        ),
        "weighted_positive_slope_mean": mp.nstr(
            mean_positive, 35, strip_zeros=False
        ),
        "weighted_negative_slope_mean": mp.nstr(
            mean_negative, 35, strip_zeros=False
        ),
        "weighted_slope_gap": mp.nstr(
            mean_positive - mean_negative, 35, strip_zeros=False
        ),
        "slope_identity_relative_error": mp.nstr(
            identity_error, 15, strip_zeros=False
        ),
        "forward_strong_separation_margin": mp.nstr(
            forward_separation, 30, strip_zeros=False
        ),
        "reverse_strong_separation_margin": mp.nstr(
            reverse_separation, 30, strip_zeros=False
        ),
        "strong_componentwise_separation_holds": (
            forward_separation > 0 or reverse_separation > 0
        ),
        "minimum_component_covering_arc_radians": mp.nstr(
            covering_arc, 30, strip_zeros=False
        ),
        "strict_common_half_plane_holds": covering_arc < mp.pi,
        "scaled_first_jet_cancellation_factor": mp.nstr(
            cancellation_factor, 30, strip_zeros=False
        ),
    }


def diagnostics(dps: int) -> dict:
    mp.mp.dps = dps
    time = mp.mpf(0)
    result_rows = []
    for bracket, cutoff, label in component_gate.ROOT_SPECS:
        crossing = lambda value: 2 * mp.re(
            phase_gate.corrected_complex_main(value, time, cutoff)
        )
        root = mp.findroot(
            crossing,
            tuple(mp.mpf(endpoint) for endpoint in bracket),
        )
        result_rows.append(component_row(root, time, cutoff, label))
    if any(
        mp.mpf(row["slope_identity_relative_error"]) >= mp.mpf("1e-45")
        for row in result_rows
    ):
        raise RuntimeError("crossing slope identity diagnostic failed")
    return {
        "role": "finite_route_shaping_diagnostics_only",
        "proof_boundary": (
            "These four t=0 corrected crossings validate the numerical "
            "evaluation of the exact identity and test stronger sufficient "
            "criteria. They prove no L>=50 or positive-time separation."
        ),
        "dps": dps,
        "rows": result_rows,
    }


def compare_diagnostics(coarse: dict, fine: dict) -> dict:
    maximum_root_delta = mp.mpf(0)
    maximum_relative_quantity_delta = mp.mpf(0)
    fields = (
        "U",
        "weighted_positive_slope_mean",
        "weighted_negative_slope_mean",
        "weighted_slope_gap",
        "minimum_component_covering_arc_radians",
        "scaled_first_jet_cancellation_factor",
    )
    for left, right in zip(
        coarse["rows"], fine["rows"], strict=True
    ):
        if left["label"] != right["label"]:
            raise RuntimeError("diagnostic row ordering drifted")
        maximum_root_delta = max(
            maximum_root_delta,
            abs(mp.mpf(left["root"]) - mp.mpf(right["root"])),
        )
        for field in fields:
            coarse_value = mp.mpf(left[field])
            fine_value = mp.mpf(right[field])
            maximum_relative_quantity_delta = max(
                maximum_relative_quantity_delta,
                abs(coarse_value - fine_value) / max(1, abs(fine_value)),
            )
        for field in (
            "strong_componentwise_separation_holds",
            "strict_common_half_plane_holds",
        ):
            if left[field] != right[field]:
                raise RuntimeError(f"diagnostic boolean drifted: {field}")
    if maximum_root_delta >= mp.mpf("1e-35"):
        raise RuntimeError("diagnostic roots are not precision stable")
    if maximum_relative_quantity_delta >= mp.mpf("1e-30"):
        raise RuntimeError("diagnostic quantities are not precision stable")
    return {
        "coarse_dps": coarse["dps"],
        "fine_dps": fine["dps"],
        "maximum_abs_root_delta": mp.nstr(
            maximum_root_delta, 25, strip_zeros=False
        ),
        "maximum_relative_quantity_delta": mp.nstr(
            maximum_relative_quantity_delta, 25, strip_zeros=False
        ),
        "root_tolerance": "1e-35",
        "relative_quantity_tolerance": "1e-30",
        "status": "passed",
    }


def render_note(artifact: dict) -> str:
    diagnostic_lines = []
    for row in artifact["diagnostics"]["rows"]:
        diagnostic_lines.append(
            "| {label} | {cutoff} | {root} | {weighted_slope_gap} | "
            "{strong_componentwise_separation_holds} | "
            "{strict_common_half_plane_holds} | "
            "{scaled_first_jet_cancellation_factor} |".format(**row)
        )
    table = "\n".join(diagnostic_lines)
    return f"""# Newman Crossing Slope-Gap Reduction

Date: {DATE}

Status: exact reduction, impossibility gate, and finite diagnostic note; not a proof of the Xi separation theorem, `Lambda<=0`, or RH.

## Real Component Jet

Write the corrected complex main as

```text
E=sum_j e_j,                    e_j=a_j exp(i theta_j),
e_j'=(u_j+i v_j)e_j,           a_j>0,
X=Re(E),                        U=Re(E').
```

Set

```text
p_j=cos(theta_j),
q_j=u_j cos(theta_j)-v_j sin(theta_j).
```

Then `X=sum a_j p_j` and `U=sum a_j q_j`. For `L>0`,

```text
Q_L=X^2+(U/L)^2
   =a^T (p p^T+q q^T/L^2) a                         (1)
   =sum_j a_j^2(p_j^2+q_j^2/L^2)
    +2 sum_(j<k) a_j a_k(p_j p_k+q_j q_k/L^2).
```

This is the requested exact diagonal/off-diagonal expansion. Its matrix has
rank at most two. On the crossing hyperplane `p.a=0`, it collapses to

```text
Q_L=(q.a)^2/L^2.                                    (2)
```

For at least three components, the kernel of (2) inside `p_perp` has
dimension at least `m-2`. Therefore no argument using only a positive lower
eigenvalue of this unrestricted crossing form can work.

## Slope-Gap Identity

Put `c_j=a_j p_j`, `d_j=a_j q_j`, and, when `c_j!=0`,
`h_j=d_j/c_j=u_j-v_j tan(theta_j)`. Split the indices into
`P={{j:c_j>0}}`, `N={{j:c_j<0}}`, and `Z={{j:c_j=0}}`, and define

```text
M_+=sum_P c_j,                 M_-=sum_N (-c_j),
M=(M_++M_-)/2,                X=M_+-M_-,
h_+=(sum_P c_j h_j)/M_+,
h_-=(sum_N (-c_j) h_j)/M_-,
D_0=sum_Z d_j.
```

Direct substitution gives the robust near-crossing identity

```text
U=M(h_+-h_-)+(X/2)(h_++h_-)+D_0.                  (3)
```

At an exact crossing, `M_+=M_-=M`, so

```text
U=M(h_+-h_-)+D_0,
Q_L=(M(h_+-h_-)+D_0)^2/L^2.                       (4)
```

Equation (4) is the surviving reduction. It uses the crossing constraint
instead of asking an indefinite full form to become positive.

A strong sufficient condition is immediate. If `D_0=0` and all positive-side
slopes exceed all negative-side slopes by at least `gamma>0`, or conversely,
then

```text
|U|>=M gamma,                 Q_L>=M^2 gamma^2/L^2. (5)
```

The actual weighted-mean gap in (4) is weaker than this pairwise ordering and
is the better theorem target.

## Half-Plane Form

Let `g_j=(p_j,q_j/L)`. If some unit vector `eta` and `kappa>0` obey

```text
eta.g_j>=kappa  for every j,
```

then positive amplitudes give

```text
sqrt(Q_L)=||sum_j a_j g_j||>=kappa sum_j a_j.      (6)
```

This is exactly the component-level version of the rotating half-plane cone
needed on a successor strip. It is sufficient, not asserted for Xi.

## Centered Arithmetic Frame

At `t=0`, the Dirichlet carriers have
`v_n=beta_0'+log(n)/2`. For any real center `mu`, put
`v_mu=beta_0'+mu/2` and
`R_(0,mu)=e_0'-i v_mu e_0` for the endpoint pseudo-component. Then

```text
E'=i v_mu E
   +(i/2) sum_(n=1)^N (log(n)-mu)e_n+R_(0,mu).    (7)
```

At `X=0`, with `E=iY`,

```text
U=-v_mu Y
  -(1/2) Im sum_(n=1)^N (log(n)-mu)e_n
  +Re R_(0,mu).                                   (8)
```

This identifies the missing input as a centered logarithmic trigonometric
moment separation with the endpoint retained, rather than as a generic
frequency-ordering claim.

## Exact Generic Obstruction

The four-carrier thought experiment

```text
E_*(x)=3 exp(i*pi/3-4ix)+3 exp(2i*pi/3-3ix)
      +7 exp(-i*pi/3-2ix)+7 exp(4i*pi/3-ix)       (9)
```

has positive amplitudes, speeds `-4<-3<-2<-1`, and no component with zero
real part at `x=0`. Nevertheless,

```text
E_*(0)=-4i sqrt(3),       E_*'(0)=-5i,
Re(E_*''(0))=-21,         M=5,
h_+=h_-=-sqrt(3)/5.                              (10)
```

Thus the slope gap can vanish exactly under all of those generic assumptions.
The new arithmetic theorem must use the actual zeta phases, amplitudes, and
endpoint relation.

## Corrected-Crossing Diagnostics

These are reproducible `t=0` route-shaping diagnostics, not interval
certificates and not evidence for the `L>=50` theorem.

| label | N | crossing | h_+-h_- | pairwise separation | common half-plane | scaled cancellation |
|---|---:|---:|---:|---|---|---:|
{table}

The coarse/fine precision comparison passed with maximum root delta
`{artifact["precision_stability"]["maximum_abs_root_delta"]}` and maximum
relative quantity delta
`{artifact["precision_stability"]["maximum_relative_quantity_delta"]}`.

## Xi Handoff

At an exact corrected-main crossing, the live target

```text
Q_L>8000000 exp(-3L/2)
```

is equivalent to

```text
|M(h_+-h_-)+D_0|
  >2000 sqrt(2) L exp(-3L/4),       L>=50.         (11)
```

This is a sharper statement of the missing theorem, not its proof. A viable
proof may bound the centered log moment in (8), prove the weighted slope gap
(11), or establish an equivalent Xi-specific half-plane cone on newly
entering high-frequency strips.

There is an important quantifier guard. If (11) is imposed uniformly at one
fixed `x` for every positive `t` tending to zero, continuity promotes it to
a nonzero endpoint first jet and therefore adds a high-zero simplicity
burden. The existing positive-boundary delta-localization gate already
proves that such a uniform floor is sufficient but not logically necessary
for `Lambda<=0`. The multiplicity-compatible alternative is to use (11), or
a weaker entry certificate, only where new cells enter and then transport
their descendants with a time-adaptive relative estimate. The finite
diagnostics, exact countermodel, Q208 theorem, and one finite successor do
not prove either route, an all-`j` successor, `Lambda<=0`, RH, PF-infinity,
or the Clay prize.

Machine-audited files:

```text
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```
"""


def main() -> int:
    coarse = diagnostics(COARSE_DPS)
    fine = diagnostics(FINE_DPS)
    artifact = {
        "kind": STEM,
        "date": DATE,
        "status": (
            "exact crossing Gram/slope-gap reduction with a generic "
            "impossibility guard, finite corrected-crossing diagnostics, "
            "and one open Xi arithmetic target"
        ),
        "proof_boundary": (
            "The Gram expansion, rank collapse, near/exact crossing "
            "slope identities, sufficient half-plane and slope-ordering "
            "bounds, centered t=0 log-moment frame, and four-carrier "
            "countermodel are exact. The diagnostics are finite only. "
            "The Xi weighted slope-gap estimate is not proved, and this "
            "does not prove an all-j successor, Lambda<=0, RH, "
            "PF-infinity, or the Clay prize."
        ),
        "source_sha256": {
            "critical_component_wronskian_gate": file_hash(COMPONENT_RESULT),
            "critical_wronskian_phase_reduction": file_hash(PHASE_RESULT),
            "single_carrier_adiabatic_benchmark": file_hash(
                SINGLE_CARRIER_RESULT
            ),
        },
        "builder_sha256": file_hash(Path(__file__)),
        "rows": rows(),
        "diagnostics": fine,
        "precision_stability": compare_diagnostics(coarse, fine),
        "exact_summary": {
            "row_count": 10,
            "exact_reductions": 8,
            "exact_countermodels": 1,
            "open_xi_targets": 1,
            "crossing_form_rank_upper_bound": 1,
            "generic_crossing_coercivity": False,
            "all_j_xi_theorem": False,
        },
        "next_exact_obligation": (
            "Use the actual corrected Riemann-Siegel phase/amplitude "
            "relations and endpoint recurrence to prove a weighted "
            "slope-gap, centered logarithmic-moment, or rotating-half-plane "
            "entry theorem on new high-frequency strips. Do not require "
            "one fixed-x floor uniformly to t=0 unless endpoint simplicity "
            "is explicitly accepted as an extra burden; combine the entry "
            "theorem with multiplicity-compatible descendant transport."
        ),
    }
    RESULT_PATH.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    NOTE_PATH.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman crossing slope-gap reduction: "
        "10 rows, 8 exact reductions, 1 exact countermodel, "
        "4 corrected-crossing diagnostics, 1 open Xi target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
