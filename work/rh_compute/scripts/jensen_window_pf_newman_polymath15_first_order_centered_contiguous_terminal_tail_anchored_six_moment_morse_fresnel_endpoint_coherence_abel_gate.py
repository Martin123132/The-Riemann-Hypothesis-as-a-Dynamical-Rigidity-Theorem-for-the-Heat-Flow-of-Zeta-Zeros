#!/usr/bin/env python3
"""Build the Morse-Fresnel endpoint-coherence and Abel route gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
    "endpoint_coherence_abel_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "aggregate_scaling": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
        "aggregate_scaling_gate.json"
    ),
    "morse_fresnel_endpoint": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
        "endpoint_reduction.json"
    ),
    "finite_poisson_transport": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_finite_poisson_"
        "transport_reduction.json"
    ),
}


@dataclass(frozen=True)
class CoherenceRow:
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
    aggregate = payloads["aggregate_scaling"].get("counts", {})
    if aggregate.get("certified_transition_collars") != 1:
        raise RuntimeError("aggregate-scaling source lost its collar")
    if aggregate.get("explicit_tail_h2_scaling_bounds") != 6:
        raise RuntimeError("aggregate-scaling source lost six tail bounds")
    if aggregate.get("grouped_roster_bounds") != 0:
        raise RuntimeError("aggregate-scaling source overpromoted the grouped roster")
    morse = payloads["morse_fresnel_endpoint"].get("counts", {})
    if morse.get("exact_transition_integrals") != 6:
        raise RuntimeError("Morse-Fresnel source lost six transition integrals")
    poisson = payloads["finite_poisson_transport"].get("counts", {})
    if poisson.get("bare_logarithmic_moments") != 6:
        raise RuntimeError("finite-Poisson source lost six moments")
    if poisson.get("stationary_transport_identities") != 5:
        raise RuntimeError("finite-Poisson source lost five transport identities")
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def phase_certificate() -> dict:
    alpha, r = sp.symbols("alpha r", positive=True)
    gauge = alpha * (sp.log(alpha / r) - 1) + r
    first = sp.expand_log(gauge.subs(r, r + 1) - gauge, force=True)
    first_expected = 1 - alpha * (sp.log(r + 1) - sp.log(r))
    if sp.simplify(first - first_expected) != 0:
        raise RuntimeError("first phase difference failed")
    second = sp.expand_log(
        gauge.subs(r, r + 2) - 2 * gauge.subs(r, r + 1) + gauge,
        force=True,
    )
    second_expected = alpha * (
        2 * sp.log(r + 1) - sp.log(r) - sp.log(r + 2)
    )
    if sp.simplify(second - second_expected) != 0:
        raise RuntimeError("second phase difference failed")

    v = r / alpha
    defect = v - 1 - sp.log(v)
    saddle_phase = alpha * (sp.log(alpha / r) - 1)
    lower = sp.expand_log(saddle_phase - alpha * defect + r, force=True)
    if sp.simplify(lower) != 0:
        raise RuntimeError("lower endpoint phase collapse failed")
    B = sp.symbols("B", positive=True)
    upper_v = B * r / alpha
    upper_defect = upper_v - 1 - sp.log(upper_v)
    upper = sp.expand_log(
        saddle_phase
        - alpha * upper_defect
        - (alpha * sp.log(B) - r * B),
        force=True,
    )
    if sp.simplify(upper) != 0:
        raise RuntimeError("upper endpoint phase collapse failed")

    return {
        "phase_gauge": (
            "G_alpha(r)=alpha[log(alpha/r)-1]+r; because r is integral, "
            "e(G_alpha(r))=e(phi_r(alpha/r))."
        ),
        "first_difference": "Delta G=1-alpha log(1+1/r)",
        "second_difference": (
            "Delta^2 G=alpha log((r+1)^2/[r(r+2)])="
            "alpha log(1+1/[r(r+2)])"
        ),
        "curvature_envelope": (
            "On the collar and alpha>=4096, "
            "1/(2alpha)<=Delta^2G<=5/(4alpha)."
        ),
        "endpoint_resonance": (
            "Delta G is strictly increasing and approaches the integer "
            "frequency zero at r=floor(alpha); no uniform first-difference "
            "gap exists."
        ),
        "lower_endpoint_phase": "phi_r(alpha/r)-y_1^2/2=-r",
        "upper_endpoint_phase": "phi_r(alpha/r)-y_B^2/2=alpha log(B)-rB",
    }


def abel_certificate() -> dict:
    spread_margin_at_64 = Fraction(64**2, 4) - Fraction(301 * 64, 20) - 1
    if spread_margin_at_64 != Fraction(299, 5):
        raise RuntimeError("terminal spread margin failed")
    sector_floor = Fraction(99, 100)
    core_floor = sector_floor * Fraction(13, 100) * Fraction(1, 20)
    if core_floor != Fraction(1287, 200_000):
        raise RuntimeError("terminal core floor failed")
    h2_ratio = core_floor * 2**25
    if h2_ratio <= 200_000:
        raise RuntimeError("terminal core h^2 separation failed")
    return {
        "outer_split": (
            "Let R=floor(alpha), H=ceil(sqrt(alpha)), and split before "
            "R-H. On the outer part delta_r=Delta G lies in (-1/9,0)."
        ),
        "summation_identity": (
            "For z_r=e(G(r)) and A_r=[e(delta_r)-1]^(-1), "
            "sum_(m<=r<=n)z_r=A_nz_(n+1)-A_mz_m+"
            "sum_(m<r<=n)(A_(r-1)-A_r)z_r."
        ),
        "vertical_monotonicity": (
            "A(delta)=-1/2-(i/2)cot(pi delta). Since delta_r increases "
            "inside (-1/2,0), its total variation telescopes vertically."
        ),
        "outer_bound": (
            "Every outer interval ending at n<=R-H has phase sum at most "
            "1/|delta_n|<=(alpha-H+1)/(H-1)."
        ),
        "partial_sum_bound": (
            "Every contiguous interval in the collar has phase sum "
            "strictly less than 3sqrt(alpha)."
        ),
        "weighted_abel": (
            "For arbitrary complex weights w_r, "
            "|sum w_r e(G(r))|<3sqrt(alpha){|w_n|+sum|w_(r+1)-w_r|}."
        ),
        "generic_scale": (
            "A bounded-variation norm C/alpha therefore gives only "
            "3C/sqrt(alpha), an h-scale estimate rather than h^2."
        ),
        "terminal_core": (
            "Put K=floor(sqrt(alpha)/20). For 0<=k<=K, "
            "0<=G(R-k)-G(R)<1/300."
        ),
        "sector_floor": str(sector_floor),
        "actual_c_core_floor": (
            "Using c_(0,r)(y_1)>=13/(100r), "
            "|sum_(k=0)^K c_(0,R-k)(y_1)e(G(R-k))|>"
            "1287/[200000sqrt(alpha)]>=1287h/200000."
        ),
        "core_floor": str(core_floor),
        "h2_ratio_floor": str(h2_ratio),
        "h2_separation": (
            "Since a>2^25, the certified terminal subblock exceeds "
            "200000h^2. This is a subblock/Abel nonpromotion guard, not a "
            "lower bound for the complete signed transition residual."
        ),
    }


def endpoint_certificate() -> dict:
    return {
        "endpoint_c_formula": (
            "For s=log(alpha/r), Z(s)=sqrt(2[exp(-s)-1+s]), "
            "c_(j,r)(y_1)={A_j(exp(s))-A_j(1)J(exp(-s))}/[rZ(s)]."
        ),
        "residual_identity": (
            "With kappa=1/(2pi i), exact integration by parts gives "
            "Q_(j,r)=kappa c_(j,r)(y_1)-kappa e(alpha log B)"
            "c_(j,r)(y_B)+kappa e(phi_r(alpha/r))"
            "integral_(y_1)^(y_B)c'_(j,r)(y)e(-y^2/2)dy."
        ),
        "lower_coherence": (
            "The lower boundary phase is e(-r)=1 for every integral mode."
        ),
        "upper_coherence": (
            "Because B and r are integral, the upper boundary phase is the "
            "common value e(alpha log B)."
        ),
        "coherent_trace_floor": (
            "On R_alpha, |kappa sum c_(0,r)(y_1)|>91/96800. This is one "
            "exact coherent component; the upper trace and c' integral may "
            "cancel it in the full residual."
        ),
        "composition_requirement": (
            "Neither saddle-phase Abel cancellation nor Hermitian pairing "
            "acts on the two separated coherent endpoint traces. They must "
            "remain joined to the physical endpoint E_j, the extracted "
            "calB_j endpoint functionals, and the interior c' integral."
        ),
    }


def exact_payload() -> dict:
    return {
        "domain": (
            "Use q=1, L>=50, a^2=exp(L)+1/(32L^2), h=1/a, "
            "alpha=xi/(2pi) on T_0-epsilon<=xi<=T_0, and the lower collar "
            "ceil(alpha exp(-1/10))<=r<=floor(alpha). Then alpha>2^49>4096."
        ),
        "route_decision": (
            "Do not seek the missing h^2 theorem from a coefficient-blind "
            "first- or second-difference estimate for the saddle phase. Its "
            "endpoint stationary core naturally has h scale, and residual "
            "integration by parts removes the saddle phase entirely from both "
            "physical boundary traces. Reassemble the complete lower and "
            "upper endpoint packages before estimating; apply phase "
            "cancellation only to the still-oscillatory interior term."
        ),
        "next_target": (
            "Derive a single endpoint-composed identity for E_j, the roster "
            "boundary traces, calB_j(B)-calB_j(1), and the terminal recurrence. "
            "Test exact cancellation of its coherent pieces before bounding "
            "the c' interior by grouped saddle-phase methods. Then insert all "
            "six composed values into the eight current observations."
        ),
        "pi_provenance": (
            "The phase e(x)=exp(2pi i x) supplies every discrete phase and "
            "kappa=1/(2pi i); A(delta) uses the equivalent identity "
            "1/[e(delta)-1]=-1/2-(i/2)cot(pi delta). The rational sector "
            "comparison uses only pi<22/7. No geometric pi is fitted."
        ),
        "proof_boundary": (
            "This proves exact discrete saddle-phase differences and curvature, "
            "an explicit O(sqrt(alpha)) phase partial-sum bound, a certified "
            "h-scale terminal subblock guard for the actual lower c coefficient, "
            "and exact coherence of both residual endpoint traces. It does not "
            "lower-bound the complete signed residual, prove cancellation of "
            "the coherent endpoint package, bound either calB_j functional, "
            "prove a grouped roster or endpoint-composed h^2 estimate, evaluate "
            "the full physical remainder, prove a signed flow or Phi_B bound, "
            "exclude contact, establish a retained aggregate or Xi theorem, "
            "prove Q209, Lambda<=0, PF-infinity, RH, or a prize-level conclusion."
        ),
    }


def build_rows(
    exact: dict,
    phase: dict,
    abel: dict,
    endpoint: dict,
) -> list[CoherenceRow]:
    return [
        CoherenceRow("mfec_01_domain", "physical domain", "proved", "The endpoint-coherence audit uses the certified q=1 collar.", exact["domain"], "No sign conclusion follows."),
        CoherenceRow("mfec_02_gauge", "integer phase gauge", "proved", "Adding r preserves every mode phase.", phase["phase_gauge"], "Uses integral Poisson modes."),
        CoherenceRow("mfec_03_first", "first difference", "proved", "The gauged saddle phase has an exact first difference.", phase["first_difference"], "No gap from the integers is asserted."),
        CoherenceRow("mfec_04_second", "second difference", "proved", "The discrete curvature is exact and positive.", phase["second_difference"], "Lower collar only."),
        CoherenceRow("mfec_05_curvature", "curvature envelope", "proved", "The collar curvature is comparable to 1/alpha.", phase["curvature_envelope"], "Constants are explicit and coarse."),
        CoherenceRow("mfec_06_resonance", "endpoint resonance", "proved", "The first difference approaches zero at the physical endpoint.", phase["endpoint_resonance"], "First-derivative cancellation is nonuniform."),
        CoherenceRow("mfec_07_split", "Fresnel split", "proved", "A sqrt(alpha) terminal core isolates the resonance.", abel["outer_split"], "No modes are physically enumerated."),
        CoherenceRow("mfec_08_abel", "outer Abel identity", "proved", "The nonresonant outer phase sum has an exact telescoping form.", abel["summation_identity"], "Uses no imported exponential-sum constant."),
        CoherenceRow("mfec_09_vertical", "variation collapse", "proved", "The Abel multiplier has one-dimensional monotone variation.", abel["vertical_monotonicity"], "Valid for delta in (-1/2,0)."),
        CoherenceRow("mfec_10_outer", "outer phase bound", "proved", "Every outer phase interval has an explicit first-difference bound.", abel["outer_bound"], "Endpoint core excluded."),
        CoherenceRow("mfec_11_partial", "full phase bound", "proved", "Every collar subinterval has phase sum below 3sqrt(alpha).", abel["partial_sum_bound"], "This is square-root scale."),
        CoherenceRow("mfec_12_weighted", "weighted Abel law", "proved", "Arbitrary weights enter only through their discrete BV norm.", abel["weighted_abel"], "Coefficient-blind theorem."),
        CoherenceRow("mfec_13_scale", "Abel scale guard", "proved", "A natural C/alpha BV input yields only h scale.", abel["generic_scale"], "Does not preclude coefficient-specific cancellation."),
        CoherenceRow("mfec_14_core", "terminal core", "proved", "A physical terminal block lies in one narrow phase sector.", abel["terminal_core"], "K=floor(sqrt(alpha)/20)."),
        CoherenceRow("mfec_15_cfloor", "actual c coefficient", "proved", "The lower residual coefficient remains positive on that block.", "c_(0,r)(y_1)>=13/(100r)", "Imported from the checked scaling gate."),
        CoherenceRow("mfec_16_subblock", "partial-sum obstruction", "proved", "The actual c-weighted terminal subblock has h-scale mass.", abel["actual_c_core_floor"], abel["h2_separation"]),
        CoherenceRow("mfec_17_residual", "exact residual IBP", "proved", "Residual integration by parts exposes both physical traces exactly.", endpoint["residual_identity"], "No absolute value has been taken."),
        CoherenceRow("mfec_18_lower", "lower trace coherence", "proved", "Every lower boundary trace has the same phase.", endpoint["lower_coherence"], phase["lower_endpoint_phase"]),
        CoherenceRow("mfec_19_upper", "upper trace coherence", "proved", "Every upper boundary trace has the same phase.", endpoint["upper_coherence"], phase["upper_endpoint_phase"]),
        CoherenceRow("mfec_20_trace", "coherent trace floor", "proved", "The separated lower trace is macroscopically coherent.", endpoint["coherent_trace_floor"], "Other exact residual terms may cancel it."),
        CoherenceRow("mfec_21_route", "route decision", "proved", "Endpoint composition must precede absolute estimation.", exact["route_decision"], endpoint["composition_requirement"]),
        CoherenceRow("mfec_22_handoff", "next theorem", "open", "The endpoint-composed six-value identity is now the live target.", exact["next_target"], "No endpoint cancellation theorem is claimed."),
        CoherenceRow("mfec_23_pi", "pi provenance", "proved", "Every pi has analytic phase provenance.", exact["pi_provenance"], "No fitted geometry."),
        CoherenceRow("mfec_24_boundary", "proof boundary", "proved", "The route gate is not a grouped h^2 theorem.", exact["proof_boundary"], "This is not a proof of RH."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    phase = payload["phase_certificate"]
    abel = payload["abel_certificate"]
    endpoint = payload["endpoint_certificate"]
    return f"""# Six-Moment Morse-Fresnel Endpoint-Coherence Abel Gate

Date: 2026-08-02

Status: exact phase and endpoint-coherence route gate with `0 enumerated
physical modes`, `0 endpoint-composed bounds`, `0 signed flow bounds`, and
this is not a proof of RH.

## Physical Collar

{exact['domain']}

## Discrete Saddle Phase

{phase['phase_gauge']}

```text
{phase['first_difference']},

{phase['second_difference']}.
```

{phase['curvature_envelope']}

{phase['endpoint_resonance']}

## Abel Bound

{abel['outer_split']}

{abel['summation_identity']}

{abel['vertical_monotonicity']}

{abel['outer_bound']}

Combining the outer interval with the terminal core proves

```text
{abel['partial_sum_bound']}
```

and weighted summation by parts gives

```text
{abel['weighted_abel']}
```

{abel['generic_scale']}

## Terminal-Core Guard

{abel['terminal_core']}

The angular width is below `pi/150<1/40`, so projection onto the endpoint
phase has cosine greater than `99/100`. Together with the checked lower
coefficient floor this gives

```text
{abel['actual_c_core_floor']}
```

{abel['h2_separation']}

## Exact Endpoint Coherence

{endpoint['endpoint_c_formula']}

{endpoint['residual_identity']}

The two phase collapses are

```text
{phase['lower_endpoint_phase']},

{phase['upper_endpoint_phase']}.
```

Thus {endpoint['lower_coherence']} {endpoint['upper_coherence']}

{endpoint['coherent_trace_floor']}

{endpoint['composition_requirement']}

## Route Decision

{exact['route_decision']}

## Next Target

{exact['next_target']}

## Pi Provenance

{exact['pi_provenance']}

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
    phase = phase_certificate()
    abel = abel_certificate()
    endpoint = endpoint_certificate()
    rows = build_rows(exact, phase, abel, endpoint)
    return {
        "kind": STEM,
        "date": "2026-08-02",
        "status": (
            "endpoint-stationary Abel scale and coherent residual traces "
            "proved; endpoint-composed cancellation theorem open"
        ),
        "source_audit": source_audit(sources),
        "exact": exact,
        "phase_certificate": phase,
        "abel_certificate": abel,
        "endpoint_certificate": endpoint,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "exact_phase_difference_identities": 2,
            "curvature_envelopes": 1,
            "phase_partial_sum_bounds": 1,
            "weighted_abel_bounds": 1,
            "certified_terminal_subblocks": 1,
            "endpoint_phase_collapses": 2,
            "coherent_endpoint_trace_floors": 1,
            "enumerated_physical_modes": 0,
            "full_grouped_roster_bounds": 0,
            "endpoint_composed_bounds": 0,
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
        "built Morse-Fresnel endpoint-coherence Abel gate: "
        f"{counts['rows']} rows, "
        f"{counts['exact_phase_difference_identities']} phase differences, "
        f"{counts['phase_partial_sum_bounds']} sqrt(alpha) phase bound, "
        f"{counts['certified_terminal_subblocks']} h-scale subblock, "
        f"{counts['endpoint_phase_collapses']} coherent endpoint phases, "
        f"{counts['endpoint_composed_bounds']} endpoint-composed bounds, "
        f"{counts['signed_flow_bounds']} signed flow bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
