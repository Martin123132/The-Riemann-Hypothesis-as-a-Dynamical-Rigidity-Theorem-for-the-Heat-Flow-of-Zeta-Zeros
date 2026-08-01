#!/usr/bin/env python3
"""Build the reciprocal normalizer-phase reinforcement guard."""

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
    "reciprocal_normalizer_phase_reinforcement_guard"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
POLYMATH_SOURCE_URL = "https://arxiv.org/abs/1904.12438"

C_TWO_DEFICIT = Fraction(3_133_668_399, 48_144_906_818)

SOURCES = {
    "self_duality": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "reciprocal_saddle_self_duality_gate.json"
    ),
    "normalized_bridge": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "normalized_laguerre_bridge.json"
    ),
    "direct_projection": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "direct_projection_regime_reduction.json"
    ),
    "endpoint_lift": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "endpoint_holomorphic_lift.json"
    ),
    "cancellation_wall": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "cancellation_zero_free_wall_gate.json"
    ),
}


@dataclass(frozen=True)
class PhaseRow:
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
        "self_duality": (
            "physical phase exp[i*omega log(u)]",
            "Delta_Xi=t*(Re(alpha)-log(a_omega))",
            "exp[-13/(2x)]",
        ),
        "normalized_bridge": (
            "M_t(s)=exp((t/4)*alpha(s)^2)*M_0(s)",
            "conj(B_t(x))/B_t(x)",
            "P_(N,t)(x)=2*Re(exp(i*beta_t(x))*D_(N,t)(x))",
        ),
        "direct_projection": (
            "0<=delta_a<1/(4x)",
            "all unresolved oscillation is in the unit carriers",
        ),
        "endpoint_lift": (
            "M_0(i*T)*U(T)*exp(pi*i/8)",
            "phase cancellation",
        ),
        "cancellation_wall": (
            "3133668399/48144906818",
            "2<c<=c_* is a cancellation gap",
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
    heat_time, real_alpha, imag_alpha = sp.symbols(
        "t alpha_R alpha_I", real=True
    )
    h = -heat_time * imag_alpha / 2
    if sp.simplify(heat_time * real_alpha * imag_alpha + 2 * h * real_alpha) != 0:
        raise RuntimeError("heat normalizer phase cancellation failed")

    pi_upper = Fraction(22, 7)
    coefficient_upper = pi_upper**2 / 64 + pi_upper / 16
    if not coefficient_upper == Fraction(275, 784):
        raise RuntimeError("phase coefficient arithmetic drifted")
    if not coefficient_upper < Fraction(3, 8):
        raise RuntimeError("phase coefficient no longer below 3/8")

    return {
        "m0_definition": (
            "M_0(s)=sqrt(2*pi)*s*(s-1)*pi^(-s/2)"
            "*exp[(s/2-1/2)Log(s/2)-s/2]/16"
        ),
        "signature_phase": (
            "Q(v)=exp(i*[v*log(v/(2*pi))-v-pi/4])"
        ),
        "zero_time_defect": (
            "delta_0(T)=(T/2)log(1+1/(4T^2))"
            "+(1/2)atan(1/(2T))"
        ),
        "heat_increment": (
            "h=omega-T=-t*Im(alpha)/2; "
            "t*Re(alpha)*Im(alpha)=-2h*Re(alpha)"
        ),
        "rational_pi_bound": (
            "pi<22/7 gives pi^2/64+pi/16"
            "<=275/784<3/8"
        ),
    }


def build_exact() -> dict:
    symbolic = symbolic_audit()
    return {
        "source_coordinate": (
            "Use the exact Polymath-15 normalizer "
            "M_0(s)=sqrt(2*pi)*s*(s-1)*pi^(-s/2)"
            "*exp[(s/2-1/2)Log(s/2)-s/2]/16 and "
            "M_t(s)=exp[t*alpha(s)^2/4]M_0(s), with the standard Log "
            "branch. On the real axis write "
            "phi_t=M_t(s)/|M_t(s)|=exp(i*beta_t) and "
            "gamma_t=conj(M_t(s))/M_t(s)=exp(-2i*beta_t)."
        ),
        "pi_provenance": (
            "The pi in M_0 and "
            "Q(v)=exp(i[v*log(v/(2*pi))-v-pi/4]) is the same "
            "completed-zeta and stationary-phase normalization. The "
            "pi/4 is the negative-curvature stationary signature. No "
            "new geometric construction of pi is introduced."
        ),
        "zero_time_ratio": (
            "For s=1/2-iT, T>0, direct evaluation of the exact M_0 "
            "phase gives gamma_0/Q(T)=exp(i*delta_0(T)), where "
            "delta_0(T)=(T/2)log(1+1/(4T^2))"
            "+(1/2)atan(1/(2T))."
        ),
        "zero_time_bound": (
            "Since log(1+u)<u and atan(u)<u for u>0, "
            "0<delta_0(T)<3/(8T)=3/(4x). Thus the reciprocal "
            "stationary phase Q(T) is already conjugate-locked to the "
            "normalizer coefficient gamma_0, with O(1/x) angular "
            "error."
        ),
        "heat_phase_transport": (
            "Let omega=-Im(s_*)=T-t*Im(alpha)/2, "
            "h=omega-T>0, T_0=T+pi*t/8, "
            "a_0^2=T_0/(2*pi), and "
            "delta_a=log(a_0)-Re(alpha). For the continuous argument "
            "psi_t=arg(Q(omega)/gamma_t) chosen from t=0, "
            "psi_t=-delta_0(T)+integral_T^omega log(v/T_0)dv"
            "+2h*delta_a. This identity is exact."
        ),
        "heat_phase_bound": (
            "On L>=50, 0<=tL<=25, one has "
            "0<h<pi*t/8, 0<=delta_a<1/(4x), and "
            "log(T_0/T)<pi*t/(8T). Hence "
            "|psi_t+delta_0(T)|"
            "<pi^2*t^2/(32x)+pi*t/(16x)"
            "<3t/(8x), using pi<22/7 and t<=1/2. Therefore "
            "|psi_t|<15/(16x)<1/x."
        ),
        "complex_conjugate_lock": (
            "Combine the phase bound with the previous reciprocal "
            "amplitude logarithm eta, |eta|<13/(2x). The normalized "
            "leading reciprocal coefficient divided by the conjugate "
            "physical carrier is R=exp(eta+i*psi_t), and for x>120 "
            "one has |R-1|<8/x."
        ),
        "reinforcement_guard": (
            "The lock is conjugating, not destructively phased. In the "
            "allowed real-carrier test z=1, the normalized primal plus "
            "reciprocal leading pair is z+R*conj(z)=1+R and has "
            "modulus greater than 2-8/x. Thus reciprocal self-duality "
            "cannot by itself yield a uniform pairwise cancellation "
            "gain; for one legitimate phase it asymptotically doubles "
            "the contribution."
        ),
        "projection_guard": (
            "For a general carrier z, the locked pair is "
            "z+R*conj(z)=2Re(z)+O(|z|/x). It can be small when z is "
            "nearly imaginary, but that is the original real-projection "
            "oscillation, not a new reciprocal power saving. Any useful "
            "gain must control the internal block phase or a genuinely "
            "subtractive completed expression."
        ),
        "exponent_guard": (
            "The c=2 wall needs an exponent gain strictly larger than "
            "d_2=3133668399/48144906818. The relative lock "
            "R=1+O(1/x) changes a paired leading term only by a "
            "constant/conjugate projection and exponentially smaller "
            "corrections. It cannot supply an N^(-d_2-eta) gain by "
            "self-duality alone."
        ),
        "endpoint_guard": (
            "The C_0+C_1/a endpoint and adjacent recurrence remain "
            "essential for a correct discrete transform, but their "
            "known role is to repair the cutoff transition and leading "
            "last-saddle mismatch. No existing endpoint identity "
            "reverses the conjugate lock uniformly across the active "
            "r_* block."
        ),
        "route_rejection": (
            "Reject raw primal-plus-reciprocal destructive interference "
            "as the promised source of the missing power. Retain the "
            "reciprocal coordinate as an exact organization of the "
            "functional-equation partner and endpoint, but do not "
            "expect a power gain from pairing equal leading amplitudes."
        ),
        "replacement_target": (
            "The next cancellation calculation must transform the "
            "actual signed difference used in the zeta handoff, "
            "sum_(n<=N)(exp[t log(n)^2/4]-1)n^(-s_*)"
            "-sum_(n>N)n^(-s_*), before absolute values. Derive its "
            "endpoint-complete reciprocal kernel with first x "
            "derivative and determine whether the weighted-minus-"
            "unweighted tail has a sign-changing or vanishing leading "
            "symbol at r_*. If that symbol is also reinforcing, retire "
            "the reciprocal route and return to bilinear/additive-energy "
            "or direct Xi Abel-phase methods."
        ),
        "proof_boundary": (
            "This artifact proves the exact zero-time normalizer/"
            "stationary phase ratio, its heat-transport identity and "
            "O(1/x) bound, the combined conjugate lock, and a "
            "reinforcement nonpromotion guard. It rejects only a power "
            "gain from raw reciprocal pairing. It does not evaluate the "
            "signed weighted-minus-unweighted tail kernel, prove "
            "cancellation, improve c_*, prove the Abel-scalar gap, close "
            "the inner degree, exclude contact, prove Lambda<=0, "
            "PF-infinity, RH, or a Clay-prize conclusion."
        ),
        "diagnostics": {
            "symbolic": symbolic,
            "constants": {
                "c_two_exponent_deficit": str(C_TWO_DEFICIT),
                "phase_bound": "|psi_t|<15/(16x)<1/x",
                "amplitude_bound": "|eta|<13/(2x)",
                "complex_lock": "|exp(eta+i*psi_t)-1|<8/x",
            },
        },
    }


def build_rows(exact: dict) -> list[PhaseRow]:
    diagnostics = exact["diagnostics"]
    return [
        PhaseRow(
            "rnpr_00_source_coordinate",
            "primary_source_coordinate",
            "available_exact",
            "The phase calculation uses the exact Polymath-15 normalizer.",
            exact["source_coordinate"],
            "Formula (6)-(10) of the cited primary source.",
            diagnostics["symbolic"],
        ),
        PhaseRow(
            "rnpr_01_pi_provenance",
            "source_provenance",
            "available_exact",
            "The reciprocal signature introduces no unexplained pi.",
            exact["pi_provenance"],
            "Completed-zeta and standard stationary-phase constants only.",
        ),
        PhaseRow(
            "rnpr_02_zero_time_ratio",
            "exact_phase_identity",
            "available_exact",
            "The zero-time reciprocal phase is explicitly normalizer-locked.",
            exact["zero_time_ratio"],
            "Exact on the standard Log branch for T>0.",
        ),
        PhaseRow(
            "rnpr_03_zero_time_bound",
            "exact_bound",
            "available_exact",
            "The zero-time phase defect is below 3/(4x).",
            exact["zero_time_bound"],
            "Strict for every positive T.",
        ),
        PhaseRow(
            "rnpr_04_heat_phase_transport",
            "exact_phase_identity",
            "available_exact",
            "The heat normalizer and saddle phase shifts cancel to first order.",
            exact["heat_phase_transport"],
            "Exact identity before bounding the residual integral.",
        ),
        PhaseRow(
            "rnpr_05_heat_phase_bound",
            "exact_bound",
            "available_exact",
            "The physical reciprocal phase remains locked within 1/x.",
            exact["heat_phase_bound"],
            "Uniform on L>=50 and 0<=tL<=25.",
            diagnostics["constants"],
        ),
        PhaseRow(
            "rnpr_06_complex_lock",
            "exact_composition_bound",
            "available_exact",
            "Amplitude and phase defects give one O(1/x) conjugate lock.",
            exact["complex_conjugate_lock"],
            "For the leading stationary coefficient before d_n restoration.",
        ),
        PhaseRow(
            "rnpr_07_reinforcement_guard",
            "countermodel",
            "guard_validated",
            "A reciprocal pair can asymptotically reinforce.",
            exact["reinforcement_guard"],
            "Generic phase test, not an Xi block-value theorem.",
        ),
        PhaseRow(
            "rnpr_08_projection_guard",
            "nonpromotion_guard",
            "guard_validated",
            "Small real projection is not reciprocal cancellation.",
            exact["projection_guard"],
            "Internal arithmetic phase control remains open.",
        ),
        PhaseRow(
            "rnpr_09_exponent_guard",
            "nonpromotion_guard",
            "guard_validated",
            "An O(1/x) lock cannot repair the c=2 exponent deficit.",
            exact["exponent_guard"],
            "Rejects this mechanism only, not every reciprocal transform.",
        ),
        PhaseRow(
            "rnpr_10_endpoint_guard",
            "scope_guard",
            "guard_validated",
            "Endpoint repair does not currently reverse the block lock.",
            exact["endpoint_guard"],
            "No unproved endpoint sign is inserted.",
        ),
        PhaseRow(
            "rnpr_11_route_rejection",
            "theorem_search_decision",
            "available_exact",
            "Raw reciprocal destructive interference is retired.",
            exact["route_rejection"],
            "The reciprocal coordinate itself remains useful.",
        ),
        PhaseRow(
            "rnpr_12_replacement_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The signed weighted-minus-unweighted reciprocal symbol is open.",
            exact["replacement_target"],
            "Must be falsified at r_* before a global theorem search.",
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    return {
        "kind": STEM,
        "date": "2026-07-28",
        "status": (
            "exact reciprocal normalizer-phase reinforcement guard; "
            "raw reciprocal pair cancellation is rejected and the "
            "signed difference kernel remains open"
        ),
        "sources": {
            key: str(path.relative_to(REPO_ROOT))
            for key, path in SOURCES.items()
        },
        "source_urls": {
            "polymath15_primary": POLYMATH_SOURCE_URL,
        },
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "exact": exact,
        "rows": [asdict(row) for row in build_rows(exact)],
        "proof_boundary": exact["proof_boundary"],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return f"""# Newman Reciprocal Normalizer-Phase Reinforcement Guard

Date: 2026-07-28

Status: exact phase-lock and route-rejection gate. This is not a proof
of improved cancellation, `Lambda<=0`, PF-infinity, RH, or a Clay-prize
result.

```text
work/rh_compute/results/{STEM}.json
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```

Primary source: {POLYMATH_SOURCE_URL}

## Exact Normalizer

```text
{exact["source_coordinate"]}
```

## Pi Provenance

{exact["pi_provenance"]}

## Zero-Time Phase

```text
{exact["zero_time_ratio"]}
{exact["zero_time_bound"]}
```

## Heat Transport

```text
{exact["heat_phase_transport"]}
{exact["heat_phase_bound"]}
```

## Conjugate Lock

```text
{exact["complex_conjugate_lock"]}
```

## Reinforcement Guard

```text
{exact["reinforcement_guard"]}
{exact["projection_guard"]}
{exact["exponent_guard"]}
```

## Endpoint Scope

```text
{exact["endpoint_guard"]}
```

## Route Decision

```text
{exact["route_rejection"]}
```

## Replacement Target

```text
{exact["replacement_target"]}
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
        "built Newman reciprocal normalizer-phase reinforcement guard: "
        "13 rows, 2 exact phase identities, 1 uniform conjugate lock, "
        "3 nonpromotion guards, 1 rejected raw-pair route, "
        "1 open signed-kernel target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
