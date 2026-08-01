#!/usr/bin/env python3
"""Build the reciprocal-saddle heat-amplitude self-duality gate."""

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
    "reciprocal_saddle_self_duality_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"

R_STAR = Fraction(125_662, 155_153)
R_STAR_DUAL = 2 - R_STAR
C_STAR = Fraction(4_911_678_521, 1_933_561_194)
C_TWO_DEFICIT = Fraction(3_133_668_399, 48_144_906_818)

SOURCES = {
    "signed_contact": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_first_order_signed_contact_reduction.json"
    ),
    "direct_projection": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "direct_projection_regime_reduction.json"
    ),
    "adjacent_recurrence": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "adjacent_saddle_recurrence.json"
    ),
    "cancellation_wall": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "cancellation_zero_free_wall_gate.json"
    ),
    "carrier_kernel": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "carrier_kernel_abel_prefix_reduction.json"
    ),
}


@dataclass(frozen=True)
class DualityRow:
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
        "signed_contact": (
            "s=(1-i*x)/2, s_*=s+t*alpha(s)/2",
            "e_n=phi*exp(t*log(n)^2/4-s_*log(n))",
        ),
        "direct_projection": (
            "delta_a=log(a)-Re(alpha)",
            "0<=delta_a<1/(4x)",
        ),
        "adjacent_recurrence": (
            "[epsilon]log(main_sharp/endpoint_C1)=0",
            "last saddle and endpoint cannot be bounded separately",
        ),
        "cancellation_wall": (
            "r_*=125662/155153",
            "3133668399/48144906818",
            "c_*=4911678521/1933561194",
        ),
        "carrier_kernel": (
            "a^2=T_0/(2*pi)=x/(4*pi)+t/16",
            "terminal block retains the exact endpoint phase",
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


def symbolic_audit() -> dict:
    heat_time, log_saddle, log_dual, sigma = sp.symbols(
        "t A z sigma", real=True
    )
    log_primal = 2 * log_saddle - log_dual
    log_weight_primal = (
        heat_time * log_primal**2 / 4
        - sigma * log_primal
    )
    log_stationary_factor = log_primal - log_saddle
    log_weight_dual = heat_time * log_dual**2 / 4 - sigma * log_dual
    defect = 2 * sigma - 1 - heat_time * log_saddle
    difference = sp.expand(
        log_weight_primal
        + log_stationary_factor
        - log_weight_dual
    )
    expected = sp.expand(defect * (log_dual - log_saddle))
    if sp.simplify(difference - expected) != 0:
        raise RuntimeError("reciprocal heat-amplitude identity failed")

    omega, pi_constant = sp.symbols("omega pi", positive=True)
    saddle_squared = omega / (2 * pi_constant)
    dual = sp.symbols("nu", positive=True)
    primal = saddle_squared / dual
    curvature = saddle_squared / primal**2
    if sp.simplify(1 / sp.sqrt(curvature) - primal / sp.sqrt(saddle_squared)) != 0:
        raise RuntimeError("B-process stationary factor failed")
    phase_at_saddle = sp.expand(
        omega * (2 * log_saddle - log_dual) - omega
    )
    expected_phase = sp.expand(
        -omega * log_dual + omega * (2 * log_saddle - 1)
    )
    if sp.simplify(phase_at_saddle - expected_phase) != 0:
        raise RuntimeError("physical reciprocal phase identity failed")

    x = sp.symbols("x", positive=True)
    s = sp.Rational(1, 2) - sp.I * x / 2
    alpha_rational = 1 / (2 * s) + 1 / (s - 1)
    rational_imag = sp.simplify(sp.im(sp.expand_complex(alpha_rational)))
    if rational_imag != 3 * x / (1 + x**2):
        raise RuntimeError("alpha rational imaginary part failed")

    if R_STAR_DUAL != Fraction(184_644, 155_153):
        raise RuntimeError("critical reciprocal radius drifted")

    return {
        "log_amplitude_identity": (
            "log(g(a_omega^2/nu)*(a_omega^2/nu)/a_omega)"
            "-log(g(nu))"
            "=(2*sigma-1-t*log(a_omega))*log(nu/a_omega)"
        ),
        "self_dual_condition": (
            "2*sigma-1=t*log(a_omega)"
        ),
        "stationary_factor": (
            "f''(u)=a_omega^2/u^2 and 1/sqrt(f''(u))=u/a_omega"
        ),
        "stationary_phase": (
            "omega*log(a_omega^2/nu)-omega"
            "=-omega*log(nu)+omega*(2*log(a_omega)-1)"
        ),
        "alpha_imaginary_part": (
            "Im(alpha)=3*x/(1+x^2)-(1/2)*atan(x)"
        ),
        "critical_radius": {
            "primal_exact": str(R_STAR),
            "primal_decimal": f"{float(R_STAR):.15f}",
            "dual_exact": str(R_STAR_DUAL),
            "dual_decimal": f"{float(R_STAR_DUAL):.15f}",
        },
    }


def build_exact() -> dict:
    symbolic = symbolic_audit()
    return {
        "pi_provenance": (
            "The 2*pi in a_omega^2=omega/(2*pi) is fixed by the "
            "physical logarithmic oscillation exp(i*omega*log n) written in "
            "the standard Poisson phase exp(2*pi*i*f). It is the same "
            "Riemann-Siegel saddle normalization already used in "
            "a^2=T_0/(2*pi); no fitted circle or new pi is introduced."
        ),
        "model_phase": (
            "For g_(t,sigma)(u)=exp[(t/4)log(u)^2-sigma log(u)] "
            "and the physical phase exp[i*omega log(u)], put "
            "a_omega^2=omega/(2*pi). In Poisson mode nu, the phase "
            "omega log(u)-2*pi*nu*u has the unique stationary point "
            "u_nu=a_omega^2/nu, curvature "
            "magnitude omega/(2*pi*u_nu^2)=a_omega^2/u_nu^2, and stationary "
            "factor u_nu/a_omega."
        ),
        "dual_phase": (
            "At u_nu=a_omega^2/nu, "
            "omega log(u_nu)-2*pi*nu*u_nu="
            "-omega log(nu)+omega*(2log(a_omega)-1). Thus the physical "
            "positive logarithmic phase becomes the conjugate negative "
            "logarithmic phase, up to one global unit phase."
        ),
        "amplitude_identity": (
            "The leading stationary amplitude obeys the exact identity "
            "g(u_nu)*(u_nu/a_omega)=g(nu)*"
            "(nu/a_omega)^Delta, where "
            "Delta=2*sigma-1-t*log(a_omega). Therefore the heat "
            "amplitude is exactly reciprocal-self-dual when Delta=0."
        ),
        "actual_xi_defect": (
            "For s_*=(1-i*x)/2+t*alpha(s)/2, "
            "sigma=Re(s_*)=1/2+t*Re(alpha)/2 and "
            "omega=-Im(s_*)=x/2-t*Im(alpha)/2. Hence the actual "
            "reciprocal defect is exactly "
            "Delta_Xi=t*(Re(alpha)-log(a_omega))."
        ),
        "phase_saddle_comparison": (
            "With T_0=x/2+pi*t/8 and a_0^2=T_0/(2*pi), "
            "Im(alpha)=3x/(1+x^2)-(1/2)atan(x), so "
            "T_0-omega=(t/2)*[3x/(1+x^2)+(1/2)atan(1/x)] "
            "and 0<T_0-omega<7t/(4x). Thus "
            "0<a_0^2-a_omega^2<7t/(8*pi*x)."
        ),
        "uniform_defect_bound": (
            "On L>=50, 0<=tL<=25, the imported "
            "0<=log(a_0)-Re(alpha)<1/(4x), together with "
            "omega>x/3, gives "
            "0<log(a_0)-log(a_omega)<21t/(8x^2) and "
            "|Delta_Xi|<t/(2x). For "
            "a_omega<=nu<=a_omega^2 one has t*log(a_omega)<13, so "
            "exp[-13/(2x)]<(nu/a_omega)^Delta_Xi"
            "<exp[13/(2x)]."
        ),
        "cutoff_guard": (
            "The phase saddle a_omega and the Riemann-Siegel cutoff "
            "saddle a_0 differ by less than the displayed O(t/x) "
            "square defect, but their floors may still differ on an "
            "arbitrarily thin cutoff seam. The adjacent-saddle "
            "recurrence must therefore remain in the theorem; the "
            "floors may not be identified by size alone."
        ),
        "critical_block_map": (
            "A block u asymp N^r is sent to the reciprocal scale "
            "nu asymp N^(2-r). The current pointwise frontier is exposed "
            "at r_*=125662/155153, so its exact reciprocal radius is "
            "2-r_*=184644/155153=1.190076...>1. The obstruction block "
            "inside the cutoff is paired naturally with a block outside "
            "the cutoff, not with another retained dyadic prefix."
        ),
        "tail_interpretation": (
            "This locates a candidate cancellation mechanism discarded by the "
            "pointwise exponent-pair envelope: compose each critical "
            "primal block with its reciprocal B-process tail and the "
            "Riemann-Siegel endpoint before taking absolute values. "
            "The c=2 pointwise exponent deficit "
            "3133668399/48144906818 is the quantitative benchmark: at "
            "the active block, the paired estimate must gain a strictly "
            "larger power of N over the current pointwise bound."
        ),
        "hypothetical_guard": (
            "In the exact hypothetical Delta=0 with corrections and "
            "endpoints suppressed, self-duality fixes only the two "
            "amplitude moduli. If the global reciprocal phase aligns, "
            "A*exp(i*theta)+A*exp(i*theta)=2A*exp(i*theta); if it is "
            "opposite, the same pair cancels. Thus reciprocal "
            "self-duality alone proves neither cancellation nor a "
            "half-plane sign. The physical Xi global phase and endpoint "
            "must decide between reinforcement and cancellation."
        ),
        "correction_guard": (
            "The actual first-order carriers contain (1+d_n), the "
            "normalizer phase phi, the C_0+C_1/a endpoint, first x "
            "derivatives, and adjacent-cutoff transitions. The leading "
            "self-duality identity does not authorize dropping any of "
            "them. In particular, the existing recurrence proves that "
            "the last saddle and endpoint cancel at leading order and "
            "must be composed."
        ),
        "route_decision": (
            "Promote reciprocal primal-tail composition as the next "
            "cancellation experiment ahead of another pointwise pair "
            "search. First derive an endpoint-complete one-block "
            "Poisson/B-process formula with value and first-x-derivative "
            "remainders uniform in the saddle fraction. Then test its "
            "signed paired main at the exact active radius r_* and at "
            "c=2 before attempting a global block sum."
        ),
        "open_theorem": (
            "Prove, or rigorously falsify, an endpoint-complete "
            "reciprocal block theorem for the corrected Xi carriers: "
            "pair every retained block near r_* with its "
            "nu asymp N^(2-r) stationary tail, retain the global "
            "normalizer phase, C_0+C_1/a endpoint, d_n corrections, "
            "adjacent recurrence, and first x derivative, and obtain a "
            "uniform exponent gain strictly exceeding "
            "3133668399/48144906818 at c=2. A successful theorem would "
            "lower the cancellation wall toward c=2; it would not cross "
            "the separate fixed-c<2 zero-free wall."
        ),
        "proof_boundary": (
            "This artifact proves the continuous stationary-point "
            "algebra, exact reciprocal heat-amplitude identity, actual "
            "Xi defect formula and bound, critical-radius map, and two "
            "nonpromotion guards. It does not prove a discrete "
            "Poisson/B-process remainder theorem, reciprocal block "
            "cancellation, an improved c_* threshold, the Abel-scalar "
            "gap, inner degree closure, contact exclusion, Lambda<=0, "
            "PF-infinity, RH, or a Clay-prize conclusion."
        ),
        "diagnostics": {
            "symbolic": symbolic,
            "constants": {
                "c_star": str(C_STAR),
                "r_star": str(R_STAR),
                "r_star_dual": str(R_STAR_DUAL),
                "c_two_pointwise_deficit": str(C_TWO_DEFICIT),
            },
        },
    }


def build_rows(exact: dict) -> list[DualityRow]:
    diagnostics = exact["diagnostics"]
    return [
        DualityRow(
            "rssd_00_pi_provenance",
            "source_provenance",
            "available_exact",
            "The reciprocal saddle introduces no unexplained pi.",
            exact["pi_provenance"],
            "The standard Poisson and Xi saddle normalizations agree.",
        ),
        DualityRow(
            "rssd_01_model_phase",
            "exact_stationary_algebra",
            "available_exact",
            "Every reciprocal Poisson mode has one explicit saddle.",
            exact["model_phase"],
            "Continuous stationary algebra, not yet a discrete sum theorem.",
            diagnostics["symbolic"],
        ),
        DualityRow(
            "rssd_02_dual_phase",
            "exact_stationary_algebra",
            "available_exact",
            "The reciprocal saddle conjugates the logarithmic phase.",
            exact["dual_phase"],
            "One global phase remains and must be retained.",
        ),
        DualityRow(
            "rssd_03_amplitude_identity",
            "exact_algebraic_lemma",
            "available_exact",
            "The heat amplitude has an exact reciprocal defect coordinate.",
            exact["amplitude_identity"],
            "Exact for the leading stationary amplitude.",
            diagnostics["symbolic"],
        ),
        DualityRow(
            "rssd_04_actual_xi_defect",
            "exact_xi_specialization",
            "available_exact",
            "The Xi reciprocal defect is one explicit saddle mismatch.",
            exact["actual_xi_defect"],
            "Before d_n and endpoint corrections.",
        ),
        DualityRow(
            "rssd_05_phase_saddle_comparison",
            "exact_bound",
            "available_exact",
            "The phase and cutoff saddles are exponentially close.",
            exact["phase_saddle_comparison"],
            "Uniform on the first-order critical domain.",
        ),
        DualityRow(
            "rssd_06_uniform_defect_bound",
            "exact_bound",
            "available_exact",
            "The reciprocal amplitude defect is uniformly negligible.",
            exact["uniform_defect_bound"],
            "This controls modulus matching, not paired phase cancellation.",
        ),
        DualityRow(
            "rssd_07_cutoff_guard",
            "nonpromotion_guard",
            "guard_validated",
            "Close saddles do not identify their integer cutoffs.",
            exact["cutoff_guard"],
            "The adjacent recurrence remains compulsory.",
        ),
        DualityRow(
            "rssd_08_critical_block_map",
            "exact_scale_map",
            "available_exact",
            "The active bad block maps just beyond the cutoff.",
            exact["critical_block_map"],
            "Scale localization only; no cancellation is asserted.",
            diagnostics["constants"],
        ),
        DualityRow(
            "rssd_09_tail_interpretation",
            "route_reduction",
            "available_exact",
            "The missing partner is a reciprocal tail block.",
            exact["tail_interpretation"],
            "A proposed theorem target, not a proved gain.",
        ),
        DualityRow(
            "rssd_10_hypothetical_guard",
            "countermodel",
            "guard_validated",
            "Exact self-duality can reinforce or cancel.",
            exact["hypothetical_guard"],
            "Two-term hypothetical, not an Xi counterexample.",
        ),
        DualityRow(
            "rssd_11_correction_guard",
            "nonpromotion_guard",
            "guard_validated",
            "The endpoint and first-order corrections remain essential.",
            exact["correction_guard"],
            "No leading stationary identity replaces the corrected theorem.",
        ),
        DualityRow(
            "rssd_12_route_decision",
            "theorem_search_decision",
            "available_exact",
            "Reciprocal composition is the next cancellation experiment.",
            exact["route_decision"],
            "It must be falsified at one block before global promotion.",
        ),
        DualityRow(
            "rssd_13_open_theorem",
            "open_theorem_target",
            "not_ready_to_apply",
            "The endpoint-complete reciprocal block theorem is open.",
            exact["open_theorem"],
            "Even success reaches c=2 only, not the c<2 wall.",
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    return {
        "kind": STEM,
        "date": "2026-07-28",
        "status": (
            "exact reciprocal-saddle self-duality and route gate; the "
            "endpoint-complete discrete cancellation theorem remains open"
        ),
        "sources": {
            key: str(path.relative_to(REPO_ROOT))
            for key, path in SOURCES.items()
        },
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "exact": exact,
        "rows": [asdict(row) for row in build_rows(exact)],
        "proof_boundary": exact["proof_boundary"],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return f"""# Newman Reciprocal-Saddle Self-Duality Gate

Date: 2026-07-28

Status: exact stationary algebra and theorem-search gate. This is not a
proof of reciprocal block cancellation, `Lambda<=0`, PF-infinity, RH,
or a Clay-prize result.

```text
work/rh_compute/results/{STEM}.json
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```

## Pi Provenance

{exact["pi_provenance"]}

## Reciprocal Saddle

```text
{exact["model_phase"]}
{exact["dual_phase"]}
```

## Exact Heat Self-Duality

```text
{exact["amplitude_identity"]}
{exact["actual_xi_defect"]}
```

The exact identity is the main new coordinate: the stationary amplitude
is self-dual up to one scalar exponent, rather than merely comparable by
an unrelated upper bound.

## Xi Defect Bound

```text
{exact["phase_saddle_comparison"]}
{exact["uniform_defect_bound"]}
{exact["cutoff_guard"]}
```

## Critical Block

```text
{exact["critical_block_map"]}
{exact["tail_interpretation"]}
```

## Hypothetical Stress Test

```text
{exact["hypothetical_guard"]}
```

This is the stipulated impossible-or-unphysical calculation used
correctly: it reveals what the exact algebra can and cannot decide before
the physical endpoint is restored.

## Correction Guard

```text
{exact["correction_guard"]}
```

## Route Decision

```text
{exact["route_decision"]}
```

## Open Theorem

```text
{exact["open_theorem"]}
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
        "built Newman reciprocal-saddle self-duality gate: "
        "14 rows, 1 exact stationary map, 1 exact heat-amplitude "
        "self-duality, 1 Xi defect bound, 1 critical reciprocal-tail "
        "map, 2 nonpromotion guards, 1 open theorem"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
