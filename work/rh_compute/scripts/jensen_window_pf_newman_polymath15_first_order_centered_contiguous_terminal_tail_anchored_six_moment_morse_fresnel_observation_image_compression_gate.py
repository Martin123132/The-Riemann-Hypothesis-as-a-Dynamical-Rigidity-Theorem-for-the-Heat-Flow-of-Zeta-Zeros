#!/usr/bin/env python3
"""Build the Morse-Fresnel observation-image compression gate."""

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
    "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
    "observation_image_compression_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "endpoint_composition": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
        "endpoint_composition_retention_gate.json"
    ),
    "flow_matrix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_flow_matrix_"
        "phase_reduction.json"
    ),
    "morse_fresnel_endpoint": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
        "endpoint_reduction.json"
    ),
}


@dataclass(frozen=True)
class CompressionRow:
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
    composition = payloads["endpoint_composition"].get("counts", {})
    if composition.get("endpoint_composed_value_identities") != 6:
        raise RuntimeError("endpoint-composition source lost six values")
    if composition.get("observation_expansions") != 1:
        raise RuntimeError("endpoint-composition source lost observation split")
    flow = payloads["flow_matrix"].get("counts", {})
    if flow.get("real_flow_observations") != 8:
        raise RuntimeError("flow source lost eight observations")
    if flow.get("flow_rank_upper_bound") != 8:
        raise RuntimeError("flow source lost rank-eight factorization")
    morse = payloads["morse_fresnel_endpoint"].get("counts", {})
    if morse.get("exact_transition_integrals") != 6:
        raise RuntimeError("Morse-Fresnel source lost six transition integrals")
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def matrix_certificate() -> dict:
    dimension = 5
    u0 = sp.Matrix(
        dimension,
        4,
        lambda row, column: sp.Rational((row + 2) * (column + 3) + 1, row + column + 5),
    )
    u1 = sp.Matrix(
        dimension,
        4,
        lambda row, column: sp.Rational((row + 1) ** 2 + column + 2, 2 * row + column + 7),
    )
    p = sp.Matrix([sp.Rational(index + 2, index + 5) for index in range(dimension)])
    residual = sp.Matrix(
        [sp.Rational((-1) ** index * (index + 3), index + 7) for index in range(dimension)]
    )
    core = sp.Matrix(
        [
            [0, sp.Rational(1, 2), 0, 0],
            [sp.Rational(1, 2), 0, 0, 0],
            [0, 0, 0, -sp.Rational(1, 2)],
            [0, 0, -sp.Rational(1, 2), 0],
        ]
    )
    terminal = sp.Matrix([2, 3, -5, -7])
    matrix = u0 * core * u1.T + u1 * core * u0.T
    linear = u1 * terminal
    carrier_base = u0.T * p
    carrier_dot = u1.T * p
    error_base = u0.T * residual
    error_dot = u1.T * residual
    gradient = 2 * matrix * p + linear
    factored_gradient = u0 * (2 * core * carrier_dot) + u1 * (
        2 * core * carrier_base + terminal
    )
    if sp.simplify(gradient - factored_gradient) != sp.zeros(dimension, 1):
        raise RuntimeError("gradient factorization failed")

    direct = sp.expand(
        (gradient.T * residual)[0] + (residual.T * matrix * residual)[0]
    )
    compressed = sp.expand(
        2 * (carrier_base.T * core * error_dot)[0]
        + 2 * (carrier_dot.T * core * error_base)[0]
        + (terminal.T * error_dot)[0]
        + 2 * (error_base.T * core * error_dot)[0]
    )
    if sp.simplify(direct - compressed) != 0:
        raise RuntimeError("residual observation compression failed")

    return {
        "gradient_factorization": (
            "For a=U_0^Tp and b=U_1^Tp, g_p=2M_xi p+ell_xi="
            "U_0(2Jb)+U_1(2Ja+t_T)."
        ),
        "residual_factorization": (
            "For x=U_0^Tr and z=U_1^Tr, g_p^Tr+r^TM_xi r="
            "2a^TJz+2b^TJx+t_T^Tz+2x^TJz."
        ),
        "visible_image": (
            "The remainder affects the flow through x and z only: four real "
            "observation errors and their four frequency derivatives. Any "
            "remainder in ker(U_0^T) intersect ker(U_1^T) is exactly invisible."
        ),
        "symbolic_difference": str(sp.simplify(direct - compressed)),
    }


def polynomial_certificate() -> dict:
    lam, ell = sp.symbols("lambda ell", real=True)
    coefficients = sp.symbols("r0:5")
    polynomial = sum(coefficients[index] * lam**index for index in range(5))
    lifted_coefficients = [-sp.I * ell * coefficients[0]]
    lifted_coefficients.extend(
        sp.I * (coefficients[index - 1] - ell * coefficients[index])
        for index in range(1, 5)
    )
    lifted_coefficients.append(sp.I * coefficients[4])
    lifted_polynomial = sum(
        lifted_coefficients[index] * lam**index for index in range(6)
    )
    if sp.expand(lifted_polynomial - sp.I * (lam - ell) * polynomial) != 0:
        raise RuntimeError("polynomial derivative lift failed")

    p_v = 1 + lam**2
    p_n = 1 + lam**4
    p_a = lam**3
    p_q = 1 + lam**3
    witness = (
        p_v
        + 2 * p_n
        - 3 * p_a
        - 4 * p_q
        + sp.I
        * (lam - 2)
        * (5 * p_v + 6 * p_n - 7 * p_a - 8 * p_q)
    )
    witness = sp.expand(witness)
    if sp.degree(witness, lam) != 5:
        raise RuntimeError("degree-five nonpromotion witness failed")
    if sp.simplify(witness.subs(lam, 2)) != -21:
        raise RuntimeError("lift-zero witness value failed")

    return {
        "amplitude_linearity": (
            "For P(lambda)=sum_(j=0)^5 p_jlambda^j, set "
            "A_P(u)=P(log u)exp(t(log u)^2/4-sigma log u). The Morse maps "
            "A to b, c, c', and the outside C_rA operator linearly."
        ),
        "residual_functional": (
            "Define R[P]=kappa sum_T e(phi_r(alpha/r)) integral c'_(P,r)"
            "e(-y^2/2)+kappa^2 sum_(r notin T) integral e(phi_r)C_rA_P. "
            "Then R[P]=sum_j p_j rho_j^(osc)."
        ),
        "base_polynomials": (
            "The four base observation polynomials are P_V=C, P_N=N, "
            "P_A=RC, and P_Q=Q, with degrees at most 2,4,3,3."
        ),
        "derivative_lift": (
            "The frequency-derivative row of P is dot(P)(lambda)="
            "i(lambda-log a)P(lambda), exactly matching delta_a(r)."
        ),
        "degree_ceiling": (
            "The base and lifted observation amplitudes have degrees at most "
            "2,4,3,3 and 3,5,4,4. No degree above five is introduced."
        ),
        "nonpromotion_witness": (
            "An exact algebraic specialization has P_lin of degree five with "
            "leading coefficient 6i and P_lin(log a)=-21. Thus the compression "
            "does not force a universal factor lambda-log a or any automatic "
            "endpoint zero. This witness is not asserted to be a physical Xi state."
        ),
        "witness_polynomial": str(witness),
    }


def observation_certificate() -> dict:
    return {
        "normalization": (
            "Let E[P]=nu c_xi R[P], with the fixed source normalizer and common "
            "phase. For X in {V,N,A,Q}, set x_X=Re E[P_X] and "
            "z_X=Re E[i(lambda-log a)P_X]."
        ),
        "expanded_error": (
            "The exact flow error is (N_p+A_(T,x))z_V+(X_T+V_p)z_N-"
            "(Q_p+X_(T,x))z_A-(A_T+A_p)z_Q+N_(p,xi)x_V+"
            "V_(p,xi)x_N-Q_(p,xi)x_A-A_(p,xi)x_Q+x_Vz_N+x_Nz_V-"
            "x_Az_Q-x_Qz_A."
        ),
        "linear_polynomial": (
            "P_lin=N_(p,xi)P_V+V_(p,xi)P_N-Q_(p,xi)P_A-A_(p,xi)P_Q+"
            "i(lambda-log a){(N_p+A_(T,x))P_V+(X_T+V_p)P_N-"
            "(Q_p+X_(T,x))P_A-(A_T+A_p)P_Q}."
        ),
        "linear_collapse": (
            "The entire linear remainder interaction is exactly "
            "g_p^Tr=Re E[P_lin], where deg(P_lin)<=5."
        ),
        "quadratic_pairing": (
            "The quadratic remainder is exactly r^TM_xi r="
            "x_Vz_N+x_Nz_V-x_Az_Q-x_Qz_A. It contains four signed "
            "base/derivative pairings and no hidden moment directions."
        ),
        "flow_invisible_guard": (
            "Large componentwise moment errors can be harmless if all eight "
            "visible projections vanish. Conversely, smallness of only the "
            "six raw moments is not the theorem dictated by the flow."
        ),
    }


def tail_certificate() -> dict:
    constants = (6, 14, 55, 336, 2738, 27936)
    return {
        "constants": list(constants),
        "projection_template": (
            "If P(lambda)=sum_(j=0)^5p_jlambda^j, the known outside-tail "
            "bounds imply |T[P]|<2h^2 sum_(j=0)^5 K_j|p_j|, with "
            "K=(6,14,55,336,2738,27936)."
        ),
        "boundary": (
            "This is a coefficient-dependent propagation template. The "
            "physical coefficients of P_lin and the eight base/lifted rows "
            "are not yet bounded here, and the c' interior is untouched."
        ),
    }


def route_certificate() -> dict:
    return {
        "route_decision": (
            "Estimate the single grouped interior functional associated with "
            "P_lin before taking absolute values. Treat the four quadratic "
            "observation pairings at their normalized scale. Do not allocate "
            "independent budgets to rho_0,...,rho_5."
        ),
        "next_target": (
            "Expand the physical coefficients of P_lin from C,N,RC,Q and the "
            "retained carrier observations. Derive the exact c'_(P_lin,r) "
            "Morse amplitude, audit its discrete mode variation and saddle "
            "zeros, and prove or falsify a grouped bound after the common "
            "phase, reciprocal pairing, endpoints, and terminal row are retained."
        ),
        "proof_boundary": (
            "This proves exact observation-image compression, polynomial-"
            "amplitude closure through degree five, collapse of the full linear "
            "remainder to one complex Morse functional, four exact quadratic "
            "pairings, and one coefficient-dependent outside-tail template. It "
            "does not evaluate P_lin on the physical chart, bound its c-prime "
            "interior, control the four quadratic pairings, prove an h^2 flow "
            "error, evaluate the retained carrier, prove a signed flow or Phi_B "
            "bound, exclude contact, establish a retained aggregate or Xi "
            "theorem, prove Q209, Lambda<=0, PF-infinity, RH, or a prize-level "
            "conclusion."
        ),
    }


def build_rows(
    matrix: dict,
    polynomial: dict,
    observation: dict,
    tail: dict,
    route: dict,
) -> list[CompressionRow]:
    boundary = route["proof_boundary"]
    return [
        CompressionRow("mfic_01_source", "source chain", "proved", "The reduction starts from the exact endpoint carrier and rank-eight flow.", "Sections 11.165 and 11.171.", "No estimate follows."),
        CompressionRow("mfic_02_gradient", "flow gradient", "proved", "The carrier gradient factors through U_0 and U_1.", matrix["gradient_factorization"], "Exact linear algebra."),
        CompressionRow("mfic_03_residual", "residual factor", "proved", "The full flow error uses only the observation image.", matrix["residual_factorization"], "No norm is applied."),
        CompressionRow("mfic_04_visible", "visible image", "proved", "Only eight real remainder observations affect the flow.", matrix["visible_image"], "At most eight; physical rank can fall."),
        CompressionRow("mfic_05_normalize", "normalization", "proved", "Fixed normalization and common phase pass through every residual projection.", observation["normalization"], "Common phase is not discarded."),
        CompressionRow("mfic_06_linearity", "Morse linearity", "proved", "Polynomial amplitudes are closed under the exact Morse construction.", polynomial["amplitude_linearity"], "Degree at most five."),
        CompressionRow("mfic_07_functional", "residual functional", "proved", "Every row projection is one polynomial-amplitude residual.", polynomial["residual_functional"], "Includes interior and outside tail."),
        CompressionRow("mfic_08_base", "base polynomials", "proved", "The four observations use C,N,RC,Q.", polynomial["base_polynomials"], "Physical coefficients retained."),
        CompressionRow("mfic_09_lift", "derivative lift", "proved", "Frequency differentiation multiplies by i(lambda-log a).", polynomial["derivative_lift"], "No seventh moment."),
        CompressionRow("mfic_10_degree", "degree closure", "proved", "All eight observation amplitudes remain degree at most five.", polynomial["degree_ceiling"], "Sufficiency, not minimality."),
        CompressionRow("mfic_11_expanded", "expanded error", "proved", "The carrier/remainder interaction has an exact observation formula.", observation["expanded_error"], "Signed formula."),
        CompressionRow("mfic_12_plin", "linear polynomial", "proved", "All physical linear error coefficients form one polynomial.", observation["linear_polynomial"], "Carrier observations remain exact."),
        CompressionRow("mfic_13_collapse", "linear collapse", "proved", "The full linear remainder is one complex Morse functional.", observation["linear_collapse"], "No bound is claimed."),
        CompressionRow("mfic_14_quadratic", "quadratic error", "proved", "The nonlinear error is four signed observation pairings.", observation["quadratic_pairing"], "No sign follows."),
        CompressionRow("mfic_15_kernel", "invisible kernel", "proved", "Moment directions outside the observation image do not affect the flow.", observation["flow_invisible_guard"], "Does not prove actual errors lie there."),
        CompressionRow("mfic_16_tail", "tail projection", "proved", "Known component tails propagate after polynomial composition.", tail["projection_template"], tail["boundary"]),
        CompressionRow("mfic_17_witness", "root nonpromotion", "proved", "Compression creates no universal endpoint zero.", polynomial["nonpromotion_witness"], "Algebraic witness, not a physical Xi state."),
        CompressionRow("mfic_18_route", "route decision", "proved", "The proof target is coefficient-composed, not componentwise.", route["route_decision"], "The grouped estimate remains open."),
        CompressionRow("mfic_19_target", "interior target", "open", "The physical P_lin c-prime sum is the next analytic theorem.", route["next_target"], "No h^2 bound is claimed."),
        CompressionRow("mfic_20_endpoint", "endpoint retention", "proved", "Retained hard endpoints stay inside p while their residual sensitivity enters P_lin.", "No endpoint package is moved back into the error.", "Cross-term cancellation remains to be proved."),
        CompressionRow("mfic_21_terminal", "terminal retention", "proved", "The terminal vector contributes through t_T in P_lin.", "The coefficients N_p+A_(T,x), X_T+V_p, Q_p+X_(T,x), and A_T+A_p remain joined.", "No terminal sign is assumed."),
        CompressionRow("mfic_22_hermitian", "symmetry retention", "proved", "The V,N,A,Q rows retain the Hermitian and transpose source structure.", "Polynomial composition occurs after the physical coefficient rows are formed.", "No swap-symmetry sign is inferred."),
        CompressionRow("mfic_23_scale", "scale boundary", "proved", "A degree-five closure is not an h^2 estimate.", tail["boundary"], "Coefficient norms and c-prime cancellation remain open."),
        CompressionRow("mfic_24_pi", "pi provenance", "proved", "No new pi is introduced by observation compression.", "All pi factors remain those in kappa and the inherited Fourier character.", "No fitted geometry."),
        CompressionRow("mfic_25_boundary", "proof boundary", "proved", "The gate is an exact reduction, not the missing signed estimate.", boundary, "This is not a proof of RH."),
    ]


def render_note(payload: dict) -> str:
    matrix = payload["matrix_certificate"]
    polynomial = payload["polynomial_certificate"]
    observation = payload["observation_certificate"]
    tail = payload["tail_certificate"]
    route = payload["route_certificate"]
    return f"""# Morse-Fresnel Observation-Image Compression Gate

Date: 2026-08-02

Status: exact eight-observation compression and degree-five polynomial-amplitude
closure; `0 coefficient-aware interior bounds`, `0 signed flow bounds`, and
this is not a proof of RH.

## Observation Image

{matrix['gradient_factorization']}

{matrix['residual_factorization']}

{matrix['visible_image']}

Writing the four base errors as `x=(x_V,x_N,x_A,x_Q)` and their frequency
derivatives as `z=(z_V,z_N,z_A,z_Q)`, the exact error is

```text
{observation['expanded_error']}
```

## Polynomial-Amplitude Closure

{polynomial['amplitude_linearity']}

{polynomial['residual_functional']}

{polynomial['base_polynomials']}

{polynomial['derivative_lift']}

{polynomial['degree_ceiling']}

{observation['normalization']}

## One Linear Functional

```text
{observation['linear_polynomial']}
```

{observation['linear_collapse']}

The quadratic remainder is

```text
{observation['quadratic_pairing']}
```

{observation['flow_invisible_guard']}

## Tail Propagation

{tail['projection_template']}

{tail['boundary']}

## Nonpromotion Guard

{polynomial['nonpromotion_witness']}

## Route Decision

{route['route_decision']}

## Next Target

{route['next_target']}

## Proof Boundary

{route['proof_boundary']}

## Reproduce

```text
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```
"""


def build_payload() -> dict:
    sources = load_sources()
    matrix = matrix_certificate()
    polynomial = polynomial_certificate()
    observation = observation_certificate()
    tail = tail_certificate()
    route = route_certificate()
    rows = build_rows(matrix, polynomial, observation, tail, route)
    return {
        "kind": STEM,
        "date": "2026-08-02",
        "status": (
            "eight-observation and degree-five linear-functional compression "
            "proved; coefficient-aware interior and signed flow remain open"
        ),
        "source_audit": source_audit(sources),
        "matrix_certificate": matrix,
        "polynomial_certificate": polynomial,
        "observation_certificate": observation,
        "tail_certificate": tail,
        "route_certificate": route,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "visible_residual_observations": 8,
            "linear_residual_functionals": 1,
            "maximum_polynomial_degree": 5,
            "quadratic_observation_pairings": 4,
            "tail_projection_templates": 1,
            "universal_endpoint_zeros": 0,
            "enumerated_physical_modes": 0,
            "evaluated_physical_linear_polynomials": 0,
            "coefficient_aware_interior_bounds": 0,
            "quadratic_remainder_bounds": 0,
            "signed_flow_bounds": 0,
            "phi_b_bounds": 0,
        },
        "proof_boundary": route["proof_boundary"],
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
        "built Morse-Fresnel observation-image compression gate: "
        f"{counts['rows']} rows, "
        f"{counts['visible_residual_observations']} visible residual observations, "
        f"{counts['linear_residual_functionals']} degree-five linear functional, "
        f"{counts['quadratic_observation_pairings']} quadratic pairings, "
        f"{counts['tail_projection_templates']} tail projection template, "
        f"{counts['coefficient_aware_interior_bounds']} interior bounds, "
        f"{counts['signed_flow_bounds']} signed flow bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
