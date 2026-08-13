#!/usr/bin/env python3
"""Build the Morse-Fresnel aggregate-scaling route gate."""

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
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
    "aggregate_scaling_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "morse_fresnel_endpoint": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
        "endpoint_reduction.json"
    ),
    "physical_q1": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_moment_q1_saddle_"
        "phase_variance_reduction.json"
    ),
    "growing_terminal_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_growing_prefix_finite_height_gate.json"
    ),
}


@dataclass(frozen=True)
class ScalingRow:
    id: str
    role: str
    readiness: str
    claim: str
    certificate: str
    proof_boundary: str


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCE_PATHS.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCE_PATHS.items()
    }


def source_audit(payloads: dict[str, dict]) -> dict:
    morse = payloads["morse_fresnel_endpoint"].get("counts", {})
    if morse.get("exact_transition_integrals") != 6:
        raise RuntimeError("Morse-Fresnel source lost six transition integrals")
    if morse.get("explicit_tail_integral_majorants") != 6:
        raise RuntimeError("Morse-Fresnel source lost six tail majorants")
    physical = payloads["physical_q1"].get("counts", {})
    if physical.get("physical_q1_parameter_laws") != 4:
        raise RuntimeError("physical q=1 source drifted")
    prefix = payloads["growing_terminal_prefix"].get("counts", {})
    if prefix.get("finite_height_current_theorems") != 1:
        raise RuntimeError("growing-prefix source drifted")
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def geometry_certificate() -> dict:
    s = sp.symbols("s", positive=True)
    exponential = sp.exp(s)
    p_geometry = sp.simplify(
        2 * (exponential * (s - 1) + 1) / (exponential - 1) ** 2
    )
    q_geometry = sp.simplify((p_geometry - sp.diff(p_geometry, s)) / 2)
    q_expected = (
        2 * s * sp.exp(2 * s)
        - 3 * sp.exp(2 * s)
        + 4 * sp.exp(s)
        - 1
    ) / (sp.exp(s) - 1) ** 3
    if sp.simplify(q_geometry - q_expected) != 0:
        raise RuntimeError("Morse geometry Q identity failed")
    p_difference = sp.simplify(
        (sp.exp(s) - 1) ** 2
        - 2 * (sp.exp(s) * (s - 1) + 1)
    )
    p_expected = sp.exp(2 * s) - 2 * s * sp.exp(s) - 1
    if sp.simplify(p_difference - p_expected) != 0:
        raise RuntimeError("Morse geometry P inequality identity failed")
    p_factor = 2 * sp.exp(s) * (sp.sinh(s) - s)
    if sp.simplify(p_expected - sp.expand(p_factor.rewrite(sp.exp))) != 0:
        raise RuntimeError("Morse geometry P sinh factor failed")

    r_expression = sp.expand(
        25
        * (
            2 * s * sp.exp(2 * s)
            - 3 * sp.exp(2 * s)
            + 4 * sp.exp(s)
            - 1
        )
        - 16 * (sp.exp(s) - 1) ** 3
    )
    r_expected = (
        50 * s * sp.exp(2 * s)
        - 16 * sp.exp(3 * s)
        - 27 * sp.exp(2 * s)
        + 52 * sp.exp(s)
        - 9
    )
    if sp.simplify(r_expression - r_expected) != 0:
        raise RuntimeError("Q-floor numerator identity failed")

    tail_bound = (
        Fraction(16000)
        * Fraction(3, 10) ** 8
        / math.factorial(8)
        * Fraction(30, 29)
    )
    lower = (
        Fraction(2, 3)
        - Fraction(19, 60)
        - Fraction(35, 600)
        - Fraction(187, 36000)
        - Fraction(1333, 4_200_000)
        - tail_bound
    )
    if lower != Fraction(8_364_091, 29_232_000) or lower <= 0:
        raise RuntimeError("Q-floor rational reserve failed")
    return {
        "coordinate_parameter": "s=-log(v) on 0<v<=1",
        "P": str(p_geometry),
        "Q": str(q_geometry),
        "P_upper_identity": str(p_factor),
        "Q_floor_numerator": str(r_expected),
        "Q_floor": "16/25",
        "Q_floor_tail_bound": str(tail_bound),
        "Q_floor_reserve": str(lower),
        "coefficient_defect": (
            "D_n=16*3^n-2^(n-1)(50n-54)-52 for n>=4; "
            "D_4=76 and D_(n+1)-2D_n=16*3^n-50*2^n+52>0"
        ),
    }


def transition_certificate() -> dict:
    amplitude_floor = Fraction(3799, 4000)
    bracket_floor = Fraction(11, 80)
    derivative_floor = amplitude_floor * bracket_floor
    advertised_floor = Fraction(13, 100)
    if derivative_floor != Fraction(41_789, 320_000):
        raise RuntimeError("transition derivative floor drifted")
    if derivative_floor <= advertised_floor:
        raise RuntimeError("transition derivative reserve failed")
    harmonic_floor = Fraction(1, 22)
    pi_free_majorant_floor = (
        advertised_floor * harmonic_floor * Fraction(7, 44)
    )
    if pi_free_majorant_floor != Fraction(91, 96_800):
        raise RuntimeError("transition majorant floor drifted")
    h2_ratio_floor = pi_free_majorant_floor * 2**50
    if h2_ratio_floor <= 10**12:
        raise RuntimeError("transition h^2 separation failed")
    return {
        "collar": (
            "R_alpha={r integer: ceil(alpha*exp(-1/10))<=r<=floor(alpha)}; "
            "then u_0=alpha/r and 0<=log(u_0)<=1/10"
        ),
        "b_derivative": (
            "partial_y b_(0,r)=(A_0(u)/r){Q(v)+P(v)[t log(u)/2-sigma]}, "
            "P=J^2/v, Q=J partial_v J"
        ),
        "physical_bounds": (
            "P<=1, Q>=16/25, sigma<=201/400, and "
            "A_0(u)>=1-sigma log(u)>=3799/4000 on the collar"
        ),
        "bracket_floor": str(bracket_floor),
        "amplitude_floor": str(amplitude_floor),
        "derivative_floor_exact": str(derivative_floor),
        "derivative_floor": str(advertised_floor),
        "endpoint_c_floor": (
            "c_(0,r)(y_1)>=13/(100r), because c(y_1) is the average of "
            "partial_y b on [y_1,0]"
        ),
        "mode_count": (
            "exp(-1/10)<10/11 and alpha>=22 imply #R_alpha>=alpha/22, "
            "hence sum_(r in R_alpha)1/r>=1/22"
        ),
        "harmonic_floor": str(harmonic_floor),
        "termwise_majorant_floor": str(pi_free_majorant_floor),
        "h2_ratio_floor": str(h2_ratio_floor),
        "h2_ratio_statement": (
            "At L>=50, a^2>2^50, so the sum of the j=0 termwise absolute "
            "transition majorants exceeds 10^12*h^2. This lower-bounds the "
            "chosen majorant, not the true signed residual."
        ),
    }


def logarithmic_moment(order: int) -> Fraction:
    if order < 0:
        return Fraction(0)
    return Fraction(math.factorial(order)) / Fraction(49, 100) ** (order + 1)


def tail_constants() -> tuple[list[dict], list[int]]:
    q_bound = Fraction(51, 100)
    second_bound = Fraction(3851, 5000)
    records: list[dict] = []
    ceilings: list[int] = []
    for order in range(6):
        d2_minus_d = (
            order * (order - 1) * logarithmic_moment(order - 2)
            + Fraction(101 * order, 50) * logarithmic_moment(order - 1)
            + second_bound * logarithmic_moment(order)
        )
        first = (
            order * logarithmic_moment(order - 1)
            + q_bound * logarithmic_moment(order)
        )
        value = logarithmic_moment(order)
        aggregate = 6 * d2_minus_d + 30 * first + 74 * value
        alpha_coefficient = aggregate / 36
        ceiling = (alpha_coefficient.numerator + alpha_coefficient.denominator - 1) // alpha_coefficient.denominator
        ceilings.append(ceiling)
        records.append(
            {
                "j": order,
                "D2_minus_D_integral_bound": str(d2_minus_d),
                "D_integral_bound": str(first),
                "value_integral_bound": str(value),
                "C_j": str(aggregate),
                "C_j_decimal": f"{float(aggregate):.12f}",
                "K_j_ceiling": ceiling,
                "outside_roster_integral_tail_bound": f"<{ceiling}/alpha<={2 * ceiling}*h^2",
            }
        )
    expected = [6, 14, 55, 336, 2738, 27936]
    if ceilings != expected:
        raise RuntimeError(f"tail constants drifted: {ceilings}")
    return records, ceilings


def tail_certificate() -> dict:
    records, ceilings = tail_constants()
    return {
        "hurwitz_lemma": (
            "For k>=2 and q>=1, zeta(k,q)<=q^(-k)+q^(1-k)/(k-1)"
            "<=2q^(1-k)."
        ),
        "roster_distances": (
            "For x=alpha/u, a_x=x-m+1>=alpha/(2u) and "
            "b_x=n+1-x>=alpha_+>=alpha."
        ),
        "zeta_envelopes": (
            "Z_2<=6u/alpha, Z_3<=10u^2/alpha^2, "
            "Z_4<=18u^3/alpha^3."
        ),
        "tail_reduction": (
            "The outside-roster integral remainder is at most "
            "[4pi^2 alpha]^(-1) integral_0^(log B){"
            "6|(D^2-D)A_j|+30|DA_j|+74|A_j|}d lambda."
        ),
        "amplitude_decay": (
            "On q=1, L>=50 and 0<=lambda<=log B<log a<L, "
            "g'=t lambda/2-sigma<=-49/100, so |A_j|<="
            "lambda^j exp(-49lambda/100)."
        ),
        "derivative_bounds": (
            "With |g'|<=51/100 and g''=t/2<=1/10000, "
            "|DA_j|<=[j lambda^(j-1)+(51/100)lambda^j]e^(-49lambda/100) "
            "and |(D^2-D)A_j|<=[j(j-1)lambda^(j-2)+(101j/50)"
            "lambda^(j-1)+(3851/5000)lambda^j]e^(-49lambda/100)."
        ),
        "moment_integral": (
            "integral_0^infinity lambda^k exp(-49lambda/100)d lambda="
            "k!(100/49)^(k+1)."
        ),
        "alpha_to_h2": (
            "alpha=a^2-epsilon/(2pi)>a^2/2, so 1/alpha<2h^2. "
            "Also pi>3, hence 1/(4pi^2)<1/36."
        ),
        "constants": records,
        "K_j": ceilings,
        "endpoint_boundary": (
            "These bounds cover only the absolutely convergent integral "
            "remainder after calB_j(B)-calB_j(1) is extracted. The two "
            "endpoint functionals remain exact and unbounded here."
        ),
    }


def exact_payload() -> dict:
    return {
        "domain": (
            "Use q=1, L>=50, t=1/(2L^2), a^2=exp(L)+1/(32L^2), "
            "a=N+theta with 0<=theta<=1, B=N-M-1<a, h=1/a, and "
            "alpha=xi/(2pi) on T_0-epsilon<=xi<=T_0. Then alpha>a^2/2>22."
        ),
        "route_decision": (
            "Keep the exact Morse-Fresnel chart, but retire summing its "
            "per-mode transition majorants by absolute values as a bare "
            "h^2-scale theorem. Group roster modes through Hermitian pairing, "
            "Poisson/Abel summation, or a smooth aggregate transform before "
            "absolute values. The outside-roster integral tail may be bounded "
            "absolutely because its six explicit bounds have h^2 scaling."
        ),
        "next_target": (
            "Construct one grouped roster theorem that keeps the incomplete-"
            "Fresnel transition exact while summing its residuals with phase "
            "cancellation. Compose calB_j(B)-calB_j(1), the physical endpoint "
            "E_j, c_xi, transpose terms, and the terminal recurrence before "
            "inserting the six values into the eight observations."
        ),
        "proof_boundary": (
            "This proves a fixed positive floor for the sum of one particular "
            "termwise absolute j=0 transition majorant, and explicit O(h^2) "
            "bounds for the six outside-roster integral tails. It does not "
            "lower-bound the true signed transition residual, bound either "
            "endpoint functional, include the physical observation/current "
            "coefficients, prove a grouped roster estimate, signed flow bound, "
            "Phi_B bound, contact exclusion, retained aggregate or Xi theorem, "
            "Q209, Lambda<=0, PF-infinity, RH, or a prize-level conclusion."
        ),
    }


def build_rows(
    exact: dict,
    geometry: dict,
    transition: dict,
    tail: dict,
) -> list[ScalingRow]:
    return [
        ScalingRow("mfas_01_domain", "physical domain", "proved", "The scaling chart is the certified q=1 terminal-tail chart.", exact["domain"], "No aggregate sign follows."),
        ScalingRow("mfas_02_cutoff", "cutoff scale", "proved", "The bulk cutoff obeys B<a.", "a=N+theta and B=N-M-1 with M>=1.", "Only the inequality B<a is used."),
        ScalingRow("mfas_03_alpha", "frequency scale", "proved", "The auxiliary alpha remains larger than a^2/2.", tail["alpha_to_h2"], "No replacement alpha=a^2 is made."),
        ScalingRow("mfas_04_geometry", "Morse geometry", "proved", "The derivative geometry closes on P and Q.", f"P={geometry['P']}; Q={geometry['Q']}", "Lower endpoint collar only."),
        ScalingRow("mfas_05_p_upper", "Jacobian upper", "proved", "P is at most one on the lower collar.", geometry["P_upper_identity"], "Uses sinh(s)>=s."),
        ScalingRow("mfas_06_q_formula", "Jacobian slope", "proved", "Q=(P-P')/2 exactly.", geometry["Q"], "Removable at s=0."),
        ScalingRow("mfas_07_q_floor", "rational geometry floor", "proved", "Q>=16/25 on 0<=s<=1/10.", f"reserve={geometry['Q_floor_reserve']}; tail={geometry['Q_floor_tail_bound']}", "Exact Taylor-tail certificate."),
        ScalingRow("mfas_08_collar", "transition collar", "proved", "A macroscopic lower-endpoint saddle collar lies in the fixed roster.", transition["collar"], "No modes are enumerated."),
        ScalingRow("mfas_09_amplitude", "amplitude floor", "proved", "A_0 is uniformly positive on the collar.", transition["physical_bounds"], "Order zero only."),
        ScalingRow("mfas_10_b_derivative", "transition derivative", "proved", "The transformed derivative has an explicit 1/r factor.", transition["b_derivative"], "Exact chain rule."),
        ScalingRow("mfas_11_c_endpoint", "residual endpoint", "proved", "Each collar mode contributes a positive c(y_1) floor.", transition["endpoint_c_floor"], "Floor concerns the majorant, not Q itself."),
        ScalingRow("mfas_12_mode_count", "analytic mode count", "proved", "The collar harmonic mass is at least 1/22.", transition["mode_count"], "Floor/ceiling loss included."),
        ScalingRow("mfas_13_majorant_floor", "absolute aggregation guard", "proved", "The summed j=0 termwise majorant has a fixed positive floor.", f">{transition['termwise_majorant_floor']}", "Not a true-residual lower bound."),
        ScalingRow("mfas_14_h2_guard", "scale separation", "proved", "That chosen majorant exceeds 10^12 h^2 at L>=50.", transition["h2_ratio_statement"], "Bare value-transform comparison only."),
        ScalingRow("mfas_15_tail_coordinates", "tail distances", "proved", "Outside-roster distances grow with alpha/u.", tail["roster_distances"], "Uses the enlarged fixed roster."),
        ScalingRow("mfas_16_hurwitz", "Hurwitz compression", "proved", "All k>=2 outside sums have an elementary integral envelope.", tail["hurwitz_lemma"], "S_1 remains in the exact endpoint functional."),
        ScalingRow("mfas_17_zeta", "zeta envelopes", "proved", "Z_2,Z_3,Z_4 have explicit alpha powers.", tail["zeta_envelopes"], "Uniform in u."),
        ScalingRow("mfas_18_tail_reduction", "tail normal form", "proved", "The integral tail reduces to three logarithmic amplitude norms.", tail["tail_reduction"], "Endpoint functional excluded."),
        ScalingRow("mfas_19_decay", "amplitude decay", "proved", "Every amplitude has a fixed exponential logarithmic envelope.", tail["amplitude_decay"], "Physical q=1 chart."),
        ScalingRow("mfas_20_derivatives", "amplitude derivatives", "proved", "The two needed logarithmic derivatives have polynomial envelopes.", tail["derivative_bounds"], "Constants are deliberately coarse."),
        ScalingRow("mfas_21_constants", "six tail constants", "proved", "All six outside-roster integral tails scale as h^2.", f"K_j={tail['K_j']}", tail["endpoint_boundary"]),
        ScalingRow("mfas_22_route", "route decision", "proved", "Absolute tail, grouped transition is the surviving architecture.", exact["route_decision"], "Architectural theorem, not signed closure."),
        ScalingRow("mfas_23_handoff", "next analytic target", "open", "A grouped roster estimate is now the live scalar theorem.", exact["next_target"], "Endpoint and current coefficients remain."),
        ScalingRow("mfas_24_boundary", "proof boundary", "proved", "The scaling gate is not a signed-flow theorem.", exact["proof_boundary"], "This is not a proof of RH."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    geometry = payload["geometry_certificate"]
    transition = payload["transition_certificate"]
    tail = payload["tail_certificate"]
    constants = "\n".join(
        f"j={row['j']}: C_j={row['C_j']}, K_j={row['K_j_ceiling']}, "
        f"tail {row['outside_roster_integral_tail_bound']}"
        for row in tail["constants"]
    )
    return f"""# Six-Moment Morse-Fresnel Aggregate-Scaling Gate

Date: 2026-08-02

Status: exact scaling gate with one absolute-transition route obstruction,
six explicit outside-roster integral-tail bounds, `0 enumerated transition
modes`, `0 signed flow bounds`, and this is not a proof of RH.

## Physical Domain

{exact['domain']}

## Lower-Saddle Collar Geometry

Write `s=-log(v)`. The exact geometry is

```text
P(s)={geometry['P']},

Q(s)={geometry['Q']}.
```

The identity controlling `P<=1` is

```text
{geometry['P_upper_identity']}=2 exp(s)[sinh(s)-s]>=0.
```

For `0<=s<=1/10`, exact Taylor-tail arithmetic gives

```text
Q(s)>=16/25,

Q-floor reserve={geometry['Q_floor_reserve']},

tail allowance={geometry['Q_floor_tail_bound']}.
```

## Absolute Transition Aggregation Guard

{transition['collar']}

{transition['b_derivative']}

{transition['physical_bounds']}

Therefore

```text
partial_y b_(0,r)>=13/(100r),

c_(0,r)(y_1)>=13/(100r).
```

{transition['mode_count']}

After the exact factor `1/(2pi)` in the transition majorant and `pi<22/7`,

```text
sum_(r in R_alpha) transition_majorant_(0,r)
 >{transition['termwise_majorant_floor']}.
```

{transition['h2_ratio_statement']}

This is a lower bound on the nonnegative *majorant being proposed*, not on
the signed residual `sum_r Q_(0,r)`. It retires termwise absolute aggregation
at bare-value `h^2` scale and positively points to grouped cancellation.

## Outside-Roster Integral Tail

{tail['hurwitz_lemma']}

{tail['roster_distances']}

Hence

```text
{tail['zeta_envelopes']}
```

and the exact Section 11.168 tail envelope reduces to

```text
{tail['tail_reduction']}
```

{tail['amplitude_decay']}

{tail['derivative_bounds']}

{tail['moment_integral']}

{tail['alpha_to_h2']}

The resulting explicit constants are

```text
{constants}
```

{tail['endpoint_boundary']}

## Route Decision

{exact['route_decision']}

## Next Target

{exact['next_target']}

## Proof Boundary

{exact['proof_boundary']}

## Reproduce

```text
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```
"""


def build_payload() -> dict:
    sources = load_sources()
    exact = exact_payload()
    geometry = geometry_certificate()
    transition = transition_certificate()
    tail = tail_certificate()
    rows = build_rows(exact, geometry, transition, tail)
    return {
        "kind": STEM,
        "date": "2026-08-02",
        "status": (
            "termwise absolute transition aggregation retired at bare h^2 "
            "scale; six outside-roster integral tails proved O(h^2); grouped "
            "roster theorem open"
        ),
        "source_audit": source_audit(sources),
        "exact": exact,
        "geometry_certificate": geometry,
        "transition_certificate": transition,
        "tail_certificate": tail,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "bare_logarithmic_moments": 6,
            "certified_transition_collars": 1,
            "enumerated_transition_modes": 0,
            "termwise_absolute_transition_majorant_floors": 1,
            "majorant_h2_separation_factors": 1,
            "explicit_outside_roster_tail_bounds": 6,
            "explicit_tail_h2_scaling_bounds": 6,
            "bounded_endpoint_functionals": 0,
            "grouped_roster_bounds": 0,
            "evaluated_full_physical_remainders": 0,
            "signed_flow_bounds": 0,
            "phi_b_bounds": 0,
        },
        "proof_boundary": exact["proof_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    atomic_write(args.output, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    atomic_write(args.note, render_note(payload))
    counts = payload["counts"]
    print(
        "built Morse-Fresnel aggregate-scaling gate: "
        f"{counts['rows']} rows, "
        f"{counts['certified_transition_collars']} transition collar, "
        f"{counts['enumerated_transition_modes']} enumerated modes, "
        f"{counts['termwise_absolute_transition_majorant_floors']} absolute-majorant floor, "
        f"{counts['explicit_tail_h2_scaling_bounds']} h^2 tail bounds, "
        f"{counts['grouped_roster_bounds']} grouped roster bounds, "
        f"{counts['signed_flow_bounds']} signed flow bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
