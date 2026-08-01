#!/usr/bin/env python3
"""Build the logarithmic-phase jet and leading Turan-current reduction."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "contiguous_terminal_tail_anchored_five_moment_logarithmic_"
    "phase_turan_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "phase_anchor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "absolute_phase_anchor_reduction.json"
    ),
    "mangoldt": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_current_mangoldt_"
        "normal_form_gate.json"
    ),
    "two_carrier": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_moment_two_carrier_"
        "kernel_reduction.json"
    ),
    "growing_tail": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_growing_prefix_finite_height_gate.json"
    ),
}


@dataclass(frozen=True)
class PhaseTuranRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def require_zero(expression: sp.Expr, label: str) -> None:
    if sp.simplify(sp.expand(expression)) != 0:
        raise RuntimeError(label)


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCE_PATHS.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCE_PATHS.items()
    }


def source_audit(payloads: dict[str, dict]) -> dict:
    anchor = payloads["phase_anchor"].get("exact", {})
    unit_anchor = anchor.get("unit_anchor", "")
    relative = anchor.get("relative_coefficients", "")
    if "eta=f_1/|f_1|" not in unit_anchor:
        raise RuntimeError("phase-anchor unit source drifted")
    for marker in (
        "q_n=f_n/f_1",
        "exp[t*log(n)^2/4-s_*log(n)]",
        "q_1=1",
    ):
        if marker not in relative:
            raise RuntimeError(f"phase-anchor carrier source drifted: {marker}")

    mangoldt = payloads["mangoldt"]
    correction_free = mangoldt.get("symbolic_certificate", {}).get(
        "correction_free_moments", ""
    )
    for marker in ("(eta/S_a)sum", "log(n)^r", "0<=r<=4"):
        if marker not in correction_free:
            raise RuntimeError(f"five-moment source drifted: {marker}")
    if mangoldt.get("counts", {}).get("signed_type_ii_bounds") != 0:
        raise RuntimeError("Mangoldt proof boundary drifted")

    two_carrier = payloads["two_carrier"]
    symbolic = two_carrier.get("symbolic_certificate", {})
    for key in ("full_two_carrier_identity", "leading_q1_kernels"):
        if not symbolic.get(key):
            raise RuntimeError(f"two-carrier source missing: {key}")
    if two_carrier.get("counts", {}).get("phi_b_bounds") != 0:
        raise RuntimeError("two-carrier proof boundary drifted")

    source_bounds = payloads["growing_tail"].get("exact", {}).get(
        "source_bounds", ""
    )
    if "u_x=h^2/(8*pi)" not in source_bounds:
        raise RuntimeError("physical u_x provenance drifted")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCE_PATHS.items()
        },
        "source_kinds": {
            key: payload.get("kind", "") for key, payload in payloads.items()
        },
        "imported_signed_joint_bounds": 0,
        "imported_phi_b_bounds": 0,
    }


def symbolic_certificate() -> dict:
    i = sp.I
    z, tau, phi = sp.symbols("z tau phi", real=True)
    ell = sp.symbols("ell_1:4", real=True)
    amplitudes = sp.symbols("A_1:4", positive=True, real=True)
    omega = sp.exp(i * phi)
    fourier = sum(
        amplitude * sp.exp(-i * z * frequency)
        for amplitude, frequency in zip(amplitudes, ell)
    )
    shifted = omega * fourier
    weights = tuple(
        omega * amplitude * sp.exp(-i * tau * frequency)
        for amplitude, frequency in zip(amplitudes, ell)
    )

    jet_differences: list[str] = []
    for order in range(5):
        moment = sum(
            frequency**order * weight
            for frequency, weight in zip(ell, weights)
        )
        jet = i**order * sp.diff(shifted, z, order).subs(z, tau)
        difference = sp.simplify(jet - moment)
        require_zero(difference, f"order-{order} logarithmic jet failed")
        jet_differences.append(str(difference))

    coefficients = sp.symbols("p_0:5")
    polynomial_moment = sum(
        sum(coefficients[r] * frequency**r for r in range(5)) * weight
        for frequency, weight in zip(ell, weights)
    )
    differential_polynomial = sum(
        coefficients[r] * i**r * sp.diff(shifted, z, r).subs(z, tau)
        for r in range(5)
    )
    polynomial_difference = sp.simplify(
        differential_polynomial - polynomial_moment
    )
    require_zero(polynomial_difference, "polynomial jet functional failed")

    difference_phase = sp.simplify(
        weights[0] * sp.conjugate(weights[1])
        - amplitudes[0]
        * amplitudes[1]
        * sp.exp(-i * tau * (ell[0] - ell[1]))
    )
    sum_phase = sp.simplify(
        weights[0] * weights[1]
        - omega**2
        * amplitudes[0]
        * amplitudes[1]
        * sp.exp(-i * tau * (ell[0] + ell[1]))
    )
    require_zero(difference_phase, "phase-difference identity failed")
    require_zero(sum_phase, "phase-sum identity failed")

    x_points = sp.symbols("x_1:4", real=True)
    v_real = sp.symbols("v_1_R v_2_R v_3_R", real=True)
    v_imag = sp.symbols("v_1_I v_2_I v_3_I", real=True)
    vectors = tuple(
        real + i * imag for real, imag in zip(v_real, v_imag)
    )
    bochner_left = sum(
        vectors[j]
        * sp.conjugate(vectors[k])
        * sum(
            amplitude
            * sp.exp(-i * (x_points[j] - x_points[k]) * frequency)
            for amplitude, frequency in zip(amplitudes, ell)
        )
        for j in range(3)
        for k in range(3)
    )
    bochner_right = sum(
        amplitude
        * sum(
            vectors[j] * sp.exp(-i * x_points[j] * frequency)
            for j in range(3)
        )
        * sp.conjugate(
            sum(
                vectors[k] * sp.exp(-i * x_points[k] * frequency)
                for k in range(3)
            )
        )
        for amplitude, frequency in zip(amplitudes, ell)
    )
    bochner_difference = sp.simplify(sp.expand(bochner_left - bochner_right))
    require_zero(bochner_difference, "Bochner Gram factorization failed")

    return {
        "polar_carrier": (
            "Write s_*=sigma+i*tau and polarize the nonzero common factor "
            "eta/S_a=omega*|eta/S_a|, |omega|=1. Then "
            "w_n=omega*A_n*exp(-i*tau*ell_n), where ell_n=log(n) and "
            "A_n=|eta/S_a|*exp[t*ell_n^2/4-sigma*ell_n]>0."
        ),
        "fourier_polynomial": (
            "F(z)=sum_(n<=B)A_n*exp(-i*z*ell_n), Z(z)=omega*F(z), "
            "and G_r=i^r*Z^(r)(tau) for 0<=r<=4."
        ),
        "jet_differences": jet_differences,
        "polynomial_functional": (
            "For every polynomial P of degree at most four, "
            "sum_(n<=B)P(ell_n)w_n=[P(i*partial_z)Z(z)]_(z=tau)."
        ),
        "polynomial_difference": str(polynomial_difference),
        "pair_phases": (
            "w_n*conj(w_m)=A_nA_m*exp[-i*tau*log(n/m)], while "
            "w_nw_m=omega^2*A_nA_m*exp[-i*tau*log(nm)]."
        ),
        "bochner_factorization": (
            "For real x_j and complex v_j, "
            "sum_(j,k)v_j*conj(v_k)F(x_j-x_k)="
            "sum_(n<=B)A_n*|sum_j v_j*exp(-i*x_j*ell_n)|^2>=0."
        ),
        "bochner_difference": str(bochner_difference),
    }


def leading_current_certificate() -> dict:
    s0r, s0i, s1r, s1i, s2r, s2i, u_x = sp.symbols(
        "S0_R S0_I S1_R S1_I S2_R S2_I u_x", real=True
    )
    f0 = s0r
    f1 = s1i
    f2 = -s2r
    g0 = s0i
    leading_current = sp.expand(
        f0 * (f2 / 4 + u_x * g0 / 2) - (f1 / 2) ** 2
    )
    turan_form = sp.expand(
        (f0 * f2 - f1**2) / 4 + u_x * f0 * g0 / 2
    )
    require_zero(leading_current - turan_form, "leading Turan form failed")

    s0 = s0r + sp.I * s0i
    s1 = s1r + sp.I * s1i
    s2 = s2r + sp.I * s2i
    hermitian_sum = -(
        s2 * sp.conjugate(s0)
        + 2 * s1 * sp.conjugate(s1)
        + s0 * sp.conjugate(s2)
    ) / 8
    transpose_sum = -(
        2 * s2 * s0 - 2 * s1**2
    ) / 8 - sp.I * u_x * s0**2 / 2
    pair_current = sp.expand_complex(sp.re(hermitian_sum + transpose_sum) / 2)
    require_zero(
        pair_current - leading_current,
        "leading pair-kernel to Turan reduction failed",
    )

    y = sp.symbols("y", real=True)
    trace = sp.Function("f", real=True)(y)
    log_identity = sp.simplify(
        trace**2 * sp.diff(sp.log(trace), y, 2)
        - (trace * sp.diff(trace, y, 2) - sp.diff(trace, y) ** 2)
    )
    require_zero(log_identity, "logarithmic Turan identity failed")

    zero_fibre = sp.simplify(turan_form.subs(s0r, 0))
    if zero_fibre != -s1i**2 / 4:
        raise RuntimeError("zero-fibre sign failed")

    return {
        "shifted_trace": (
            "Put u_n=log(a/n)=L_a-ell_n and "
            "f(y)=Re[e^(-i*y*L_a)Z(tau-y)]="
            "Re sum_(n<=B)w_n*exp(-i*y*u_n); define g(y) as the "
            "imaginary part of the same sum."
        ),
        "trace_derivatives": (
            "For S_r=sum_(n<=B)u_n^r*w_n, f(0)=Re(S_0), "
            "f'(0)=Im(S_1), f''(0)=-Re(S_2), and g(0)=Im(S_0)."
        ),
        "leading_observations": (
            "At C=1, D=0, s_*'=-i/2, s_*''=0, b=-1/2, b_x=0, "
            "and chi_N=-i*u_N/2, the real bulk observations are "
            "V=f, Q=A=f'/2, and N=f''/4+(u_x/2)g at y=0."
        ),
        "leading_turan": (
            "P_bulk^(0)=V*N-A*Q="
            "[f*f''-(f')^2]/4+(u_x/2)*f*g at y=0."
        ),
        "pair_kernel_difference": str(sp.simplify(pair_current - leading_current)),
        "logarithmic_identity": (
            "On every interval where f is nonzero, "
            "f*f''-(f')^2=f^2*(log|f|)''. The division-free left side "
            "is the primary identity."
        ),
        "logarithmic_identity_difference": str(log_identity),
        "zero_fibre": (
            "If f(0)=0, then P_bulk^(0)=-(f'(0))^2/4<=0; the "
            "u_x*f*g term vanishes without division."
        ),
        "zero_fibre_expression": str(zero_fibre),
    }


def finite_guard_certificate() -> dict:
    d = sp.log(2)
    tau = sp.pi / d
    omega = -sp.Integer(1)
    amplitudes = (sp.Rational(1, 2), sp.Integer(1))
    frequencies = (sp.Integer(0), d)
    weights = tuple(
        sp.simplify(
            omega * amplitude * sp.exp(-sp.I * tau * frequency)
        )
        for amplitude, frequency in zip(amplitudes, frequencies)
    )
    if weights != (-sp.Rational(1, 2), sp.Integer(1)):
        raise RuntimeError("integer-log phase guard weights failed")

    y = sp.symbols("y", real=True)
    u = (2 * d, d)
    complex_trace = sum(
        weight * sp.exp(-sp.I * y * distance)
        for weight, distance in zip(weights, u)
    )
    f = sp.expand_complex(sp.re(complex_trace))
    g = sp.expand_complex(sp.im(complex_trace))
    f0 = sp.simplify(f.subs(y, 0))
    f1 = sp.simplify(sp.diff(f, y).subs(y, 0))
    f2 = sp.simplify(sp.diff(f, y, 2).subs(y, 0))
    g0 = sp.simplify(g.subs(y, 0))
    u_x = sp.symbols("u_x", real=True)
    current = sp.simplify(
        (f0 * f2 - f1**2) / 4 + u_x * f0 * g0 / 2
    )
    if current != d**2 / 8:
        raise RuntimeError("integer-log positive-current guard failed")

    return {
        "support": "n=1,2; a=4; ell_1=0, ell_2=log(2); u_1=2log(2), u_2=log(2)",
        "phase_choice": (
            "Take tau=pi/log(2), omega=-1, A_1=1/2, A_2=1. Then "
            "exp[-i*tau*log(2)]=exp(-i*pi)=-1, so w_1=-1/2 and w_2=1."
        ),
        "trace": "f(y)=cos(log(2)*y)-(1/2)cos(2log(2)*y)",
        "jet": (
            "f(0)=1/2, f'(0)=0, f''(0)=log(2)^2, g(0)=0"
        ),
        "leading_current": "P_bulk^(0)=log(2)^2/8>0 for every real u_x",
        "exact_current": str(current),
        "pi_provenance": (
            "The pi in this guard is chosen solely to make a half-turn on "
            "the log(2) frequency: tau*log(2)=pi and exp(-i*pi)=-1. "
            "It is not an inserted physical constant and it proves only a "
            "nonimplication. The separate physical u_x=h^2/(8*pi) is "
            "inherited from the certified growing-tail source."
        ),
        "scope": (
            "This exact two-atom example has positive amplitudes, integer "
            "logarithmic support, and one common logarithmic phase, but it "
            "is not an attained large-q=1 Xi chart. It proves that those "
            "abstract properties and Bochner positivity alone cannot sign "
            "the leading bulk current."
        ),
    }


def exact_payload() -> dict:
    return {
        "domain": (
            "Use the fixed physical q=1, L>=50 terminal/bulk chart with "
            "B=N-M-1 and the correction-free carrier "
            "w_n=(eta/S_a)exp[t*log(n)^2/4-s_*log(n)]. The source gives "
            "eta/S_a nonzero; no assumption that S_a is positive or real "
            "is made."
        ),
        "bochner_boundary": (
            "Bochner positivity controls the unweighted phase-difference "
            "matrix F(x_j-x_k). It does not control the phase-sum family "
            "w_nw_m, the polynomially weighted Hermitian kernel, the "
            "endpoint-linear sum, or their joint real projection."
        ),
        "exact_correction_decomposition": (
            "Let E=Re sum_(n<=B)L_nw_n and let H_0^+,T_0^+ be the ideal "
            "kernels H_0^+(n,m)=-(u_n+u_m)^2/8 and "
            "T_0^+(n,m)=-(lambda_n-lambda_m)^2/8-i*u_x/2. Then exactly "
            "h^2*Phi_B=E+P_bulk^(0)+R_corr, where "
            "R_corr=(1/2)Re sum_(n,m<=B){[H^+-H_0^+]w_nconj(w_m)+"
            "[T^+-T_0^+]w_nw_m}."
        ),
        "correction_boundary": (
            "Every kernel difference in R_corr vanishes at the ideal "
            "specialization, but no quantitative bound on E or R_corr is "
            "proved here. They may not be discarded or absorbed into an "
            "unnamed error."
        ),
        "sharp_target": (
            "Use the actual physical tau, amplitudes A_n, terminal "
            "aggregate, and correction coefficients to prove "
            "E+P_bulk^(0)+R_corr<=h^2/400, preferably h^2/800. A viable "
            "argument must control the phase-sum family jointly with the "
            "weighted phase-difference family and endpoint."
        ),
        "proof_boundary": (
            "This proves the positive-amplitude logarithmic Fourier-jet "
            "representation, Bochner positivity of its difference kernel, "
            "the exact leading Turan-current identity, zero-fibre sign, an "
            "integer-log common-phase nonpromotion guard, and an exact "
            "correction decomposition. It proves no signed joint bound, "
            "upper bound on Phi_B, contact exclusion, retained aggregate "
            "sign, Xi residual transfer, Q209, cofinal descendant theorem, "
            "Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }


def build_rows(
    exact: dict, symbolic: dict, leading: dict, guard: dict
) -> list[PhaseTuranRow]:
    return [
        PhaseTuranRow("ltp_01_domain", "fixed-chart domain", "proved", "The reduction uses the certified correction-free physical carrier.", exact["domain"], "No sign is imported from the source."),
        PhaseTuranRow("ltp_02_polar", "positive-amplitude polarization", "proved", "The common carrier factor can be polarized without assuming S_a is real.", symbolic["polar_carrier"], "Only nonvanishing of eta/S_a is used."),
        PhaseTuranRow("ltp_03_fourier", "logarithmic Fourier polynomial", "proved", "All five moments belong to one logarithmic-frequency exponential polynomial.", symbolic["fourier_polynomial"], "The frequencies are the actual ell_n=log(n)."),
        PhaseTuranRow("ltp_04_jets", "order-four jet", "proved", "The five complex moments are exactly the first five tau jets.", "G_r=i^r*Z^(r)(tau), 0<=r<=4", "This is an identity, not a derivative estimate.", {"differences": symbolic["jet_differences"]}),
        PhaseTuranRow("ltp_05_operator", "polynomial functional", "proved", "Every degree-at-most-four carrier polynomial is a differential functional of Z.", symbolic["polynomial_functional"], "No moment is treated as independent.", {"difference": symbolic["polynomial_difference"]}),
        PhaseTuranRow("ltp_06_phases", "ratio/product phases", "proved", "The two-carrier families have distinct logarithmic phase geometries.", symbolic["pair_phases"], "The global unit phase cancels only in the ratio family."),
        PhaseTuranRow("ltp_07_bochner", "Bochner kernel", "proved", "The unweighted phase-difference kernel is positive semidefinite.", symbolic["bochner_factorization"], "Positivity follows from A_n>0.", {"difference": symbolic["bochner_difference"]}),
        PhaseTuranRow("ltp_08_boundary", "positivity boundary", "proved", "Bochner positivity does not sign the complete current.", exact["bochner_boundary"], "No transpose or endpoint sign is inferred."),
        PhaseTuranRow("ltp_09_trace", "shifted real trace", "proved", "The leading distance moments are derivatives of one real trace.", leading["shifted_trace"] + " " + leading["trace_derivatives"], "The imaginary quadrature g is retained."),
        PhaseTuranRow("ltp_10_observations", "leading observations", "proved", "The four ideal observations reduce to f, f', f'', and g.", leading["leading_observations"], "This is only the exact ideal specialization."),
        PhaseTuranRow("ltp_11_turan", "division-free Turan current", "proved", "The ideal bulk current is a Turan numerator plus one quadrature term.", leading["leading_turan"] + " " + leading["logarithmic_identity"], "The division-free form remains primary.", {"pair_difference": leading["pair_kernel_difference"], "log_difference": leading["logarithmic_identity_difference"]}),
        PhaseTuranRow("ltp_12_zero", "zero-fibre sign", "proved", "On the physical zero fibre f(0)=0, the leading bulk current is nonpositive.", leading["zero_fibre"], "This does not control points where f is nonzero.", {"expression": leading["zero_fibre_expression"]}),
        PhaseTuranRow("ltp_13_guard", "integer-log phase guard", "proved", "Positive amplitudes and common logarithmic phase alone do not give a global sign.", guard["support"] + "; " + guard["phase_choice"] + " " + guard["leading_current"], guard["scope"], {"exact_current": guard["exact_current"], "pi_provenance": guard["pi_provenance"]}),
        PhaseTuranRow("ltp_14_correction", "exact correction decomposition", "proved", "The full current is the leading Turan current plus explicit endpoint and kernel corrections.", exact["exact_correction_decomposition"], exact["correction_boundary"]),
        PhaseTuranRow("ltp_15_target", "physical joint target", "open", "Only the actual Xi phase-amplitude-endpoint coupling can now close this route.", exact["sharp_target"], exact["proof_boundary"]),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    symbolic = payload["symbolic_certificate"]
    leading = payload["leading_current_certificate"]
    guard = payload["finite_guard_certificate"]
    return f"""# Logarithmic-Phase Jet And Turan Reduction

Date: 2026-08-01

Status: exact structural reduction with `0 signed joint bounds` and `0 Phi_B
bounds`. This is not a proof of RH; the physical endpoint/correction estimate
is open.

## Physical Carrier

{exact['domain']}

```text
{symbolic['polar_carrier']}

{symbolic['fourier_polynomial']}

{symbolic['polynomial_functional']}
```

Thus the ten real bulk coordinates are the order-four real/imaginary jet of
one positive-amplitude logarithmic-frequency polynomial, not ten freely
chosen parameters.

## Pair Phases

```text
{symbolic['pair_phases']}
```

The ratio family has the exact Bochner factorization

```text
{symbolic['bochner_factorization']}
```

{exact['bochner_boundary']}

## Shifted Trace

```text
{leading['shifted_trace']}

{leading['trace_derivatives']}

{leading['leading_observations']}
```

Consequently,

```text
{leading['leading_turan']}
```

{leading['logarithmic_identity']}

Most importantly,

```text
{leading['zero_fibre']}
```

## Exact Nonpromotion Guard

```text
{guard['support']}
{guard['phase_choice']}
{guard['trace']}
{guard['jet']}
{guard['leading_current']}
```

{guard['scope']}

### Pi Provenance

{guard['pi_provenance']}

## Full-Current Decomposition

```text
{exact['exact_correction_decomposition']}
```

{exact['correction_boundary']}

## Open Physical Target

```text
{exact['sharp_target']}
```

## Boundary

{exact['proof_boundary']}

## Reproduce

```text
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```
"""


def build_payload() -> dict:
    payloads = load_sources()
    source = source_audit(payloads)
    symbolic = symbolic_certificate()
    leading = leading_current_certificate()
    guard = finite_guard_certificate()
    exact = exact_payload()
    rows = build_rows(exact, symbolic, leading, guard)
    return {
        "kind": STEM,
        "date": "2026-08-01",
        "status": (
            "exact positive-amplitude logarithmic Fourier-jet, Bochner "
            "difference-kernel, leading Turan-current, zero-fibre sign, "
            "and correction decomposition; signed joint bound open; no "
            "Phi_B bound, Xi-level theorem, Lambda<=0, or RH"
        ),
        "proof_boundary": exact["proof_boundary"],
        "source_audit": source,
        "symbolic_certificate": symbolic,
        "leading_current_certificate": leading,
        "finite_guard_certificate": guard,
        "exact": exact,
        "counts": {
            "rows": len(rows),
            "positive_amplitude_representations": 1,
            "fourier_jet_identities": 5,
            "differential_operator_identities": 1,
            "phase_difference_families": 1,
            "phase_sum_families": 1,
            "bochner_psd_kernels": 1,
            "leading_turan_identities": 1,
            "leading_zero_fibre_signs": 1,
            "integer_log_phase_guards": 1,
            "exact_correction_decompositions": 1,
            "signed_joint_targets": 1,
            "signed_joint_bounds": 0,
            "phi_b_bounds": 0,
            "xi_level_current_theorems": 0,
        },
        "rows": [asdict(row) for row in rows],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    atomic_write(args.out, json.dumps(payload, indent=2) + "\n")
    atomic_write(args.note, render_note(payload))
    print(
        "wrote logarithmic-phase Turan reduction: "
        f"{payload['counts']['rows']} rows, 5 jet identities, "
        "1 Bochner kernel, 1 zero-fibre sign, 0 signed joint bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
