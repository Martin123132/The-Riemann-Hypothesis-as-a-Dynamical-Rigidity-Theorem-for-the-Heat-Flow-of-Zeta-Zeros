#!/usr/bin/env python3
"""Certify contiguous Airy phase packets and their grouped Green bound."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_detuning_airy_lattice_grouped_green_packet_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "cubic_resonance": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_detuning_airy_lattice_cubic_resonance_gate.json",
    "airy_green_kernel": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_beta4_detuning_airy_green_kernel_gate.json",
}

C = 159_577
Q = (C - 1) // 4
J_MIN = -198
J_MAX = 200
JS = tuple(range(J_MIN, J_MAX + 1))
N = len(JS)
PRECISION = 100
COEFFICIENTS = (
    (3, 2),
    (3, 8),
    (-1, 16),
    (3, 128),
    (-3, 256),
)
DLMF_URL = "https://dlmf.nist.gov/9.8"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def binomial_coefficients() -> tuple[arb, ...]:
    return tuple(arb(p) / q for p, q in COEFFICIENTS)


def quintic_phase(j: int, u: arb, c: arb, pi: arb) -> arb:
    xj = (arb(8) * j - 2) / c + u
    x0 = -arb(2) / c + u
    truncated = arb(0)
    for degree, coefficient in enumerate(binomial_coefficients(), start=1):
        truncated += coefficient * (xj**degree - x0**degree)
    return pi * c**2 * truncated / 12 - pi * j * (c - 1 + 2 * j)


def quintic_phase_u_derivative(j: int, u: arb, c: arb, pi: arb) -> arb:
    xj = (arb(8) * j - 2) / c + u
    x0 = -arb(2) / c + u
    derivative = arb(0)
    for degree, coefficient in enumerate(binomial_coefficients(), start=1):
        if degree >= 2:
            derivative += degree * coefficient * (xj ** (degree - 1) - x0 ** (degree - 1))
    return pi * c**2 * derivative / 12


def phase_packet_certificate(cubic: dict[str, Any], green: dict[str, Any]) -> dict[str, Any]:
    pi = arb.pi()
    c = arb(C)
    u_max = arb(1) / (2 * c**2)
    u_interval = arb(u_max / 2, u_max / 2)
    imaginary = acb(0, 1)

    center_terms = [(imaginary * quintic_phase(j, arb(0), c, pi)).exp() for j in JS]
    prefix = [acb(0)]
    for term in center_terms:
        prefix.append(prefix[-1] + term)

    max_prefix = arb(0)
    max_prefix_end = J_MIN
    for index in range(1, N + 1):
        value = abs(prefix[index]).upper()
        if value > max_prefix:
            max_prefix = value
            max_prefix_end = JS[index - 1]

    max_suffix = arb(0)
    max_suffix_start = J_MAX
    full = prefix[-1]
    for index in range(N):
        value = abs(full - prefix[index]).upper()
        if value > max_suffix:
            max_suffix = value
            max_suffix_start = JS[index]

    max_contiguous = arb(0)
    max_contiguous_indices = (J_MIN, J_MIN)
    for start in range(N):
        for stop in range(start + 1, N + 1):
            value = abs(prefix[stop] - prefix[start]).upper()
            if value > max_contiguous:
                max_contiguous = value
                max_contiguous_indices = (JS[start], JS[stop - 1])

    height_variation = arb(0)
    max_single_height_variation = arb(0)
    for j in JS:
        derivative = quintic_phase_u_derivative(j, u_interval, c, pi)
        variation = abs(derivative).upper() * u_max.upper()
        height_variation += variation
        max_single_height_variation = max(max_single_height_variation, variation)

    phase_error = arb(
        cubic["interval_certificate"]["uniform_exact_phase_normal_form_error_ball"]
    ).upper()
    normal_form_error_total = N * phase_error
    universal_perturbation = height_variation + normal_form_error_total

    exact_contiguous_upper = max_contiguous + universal_perturbation
    full_center_abs = abs(full)
    exact_full_lower = full_center_abs.lower() - universal_perturbation
    exact_full_upper = full_center_abs.upper() + universal_perturbation
    require(exact_full_lower > arb("31.4"), "full exact Airy phase packet lost its noncancellation lower bound")
    require(exact_contiguous_upper < arb("52.53"), "contiguous exact Airy phase packet exceeded 52.53")

    r_min = arb(green["interval_certificate"]["r_min_ball"])
    r_max = arb(green["interval_certificate"]["r_max_ball"])
    pointwise_kernel = arb(green["interval_certificate"]["green_kernel_absolute_envelope_ball"])
    grouped_kernel_upper = (exact_contiguous_upper / r_min.sqrt()).upper()
    termwise_kernel_upper = (N * pointwise_kernel).upper()
    gain = (termwise_kernel_upper / grouped_kernel_upper).lower()
    require(grouped_kernel_upper < arb("0.0246"), "grouped Green packet exceeded 0.0246")
    require(gain > arb("7.5"), "grouped Green packet did not improve the termwise bound by 7.5")

    # DLMF 9.8.21 has first correction +7/(32s^3) for N(-s)^2.
    # After that term the signed remainder is negative, so the displayed
    # two-term expression is an upper bound.  Combining it with the M bound
    # gives pi*M(-r)*N(-s) <= (s/r)^(1/4)*sqrt(1+7/(32s^3)).
    derivative_modulus_factor = (
        (r_max / r_min) ** (arb(1) / 4)
        * (1 + arb(7) / (32 * r_min**3)).sqrt()
    ).upper()
    grouped_derivative_kernel_upper = (exact_contiguous_upper * derivative_modulus_factor).upper()
    termwise_derivative_kernel_upper = (N * derivative_modulus_factor).upper()
    derivative_gain = (termwise_derivative_kernel_upper / grouped_derivative_kernel_upper).lower()
    require(grouped_derivative_kernel_upper < arb("52.8"), "grouped derivative Green packet exceeded 52.8")
    require(derivative_gain > arb("7.5"), "grouped derivative packet did not improve termwise bound by 7.5")

    return {
        "event_offsets": [J_MIN, J_MAX],
        "mode_roster": [Q + J_MIN, Q + J_MAX],
        "roster_size": N,
        "height_interval": "0<=u<=1/(2C^2)",
        "center_quintic_full_sum": {
            "real_ball": full.real.str(PRECISION, more=True),
            "imag_ball": full.imag.str(PRECISION, more=True),
            "absolute_ball": full_center_abs.str(PRECISION, more=True),
        },
        "center_quintic_max_prefix_upper_ball": max_prefix.str(PRECISION, more=True),
        "center_quintic_max_prefix_end_j": max_prefix_end,
        "center_quintic_max_suffix_upper_ball": max_suffix.str(PRECISION, more=True),
        "center_quintic_max_suffix_start_j": max_suffix_start,
        "center_quintic_max_contiguous_upper_ball": max_contiguous.str(PRECISION, more=True),
        "center_quintic_max_contiguous_indices": list(max_contiguous_indices),
        "uniform_sum_of_height_phase_variations_ball": height_variation.str(PRECISION, more=True),
        "uniform_max_single_height_phase_variation_ball": max_single_height_variation.str(PRECISION, more=True),
        "per_term_exact_Airy_normal_form_error_ball": phase_error.str(PRECISION, more=True),
        "full_roster_normal_form_error_sum_ball": normal_form_error_total.str(PRECISION, more=True),
        "universal_phase_packet_perturbation_ball": universal_perturbation.str(PRECISION, more=True),
        "exact_Airy_full_packet_absolute_lower_ball": exact_full_lower.str(PRECISION, more=True),
        "exact_Airy_full_packet_absolute_upper_ball": exact_full_upper.str(PRECISION, more=True),
        "exact_Airy_any_contiguous_packet_upper_ball": exact_contiguous_upper.str(PRECISION, more=True),
        "termwise_399_kernel_bound_ball": termwise_kernel_upper.str(PRECISION, more=True),
        "grouped_any_contiguous_Green_kernel_bound_ball": grouped_kernel_upper.str(PRECISION, more=True),
        "termwise_to_grouped_improvement_lower_ball": gain.str(PRECISION, more=True),
        "derivative_modulus_factor_upper_ball": derivative_modulus_factor.str(PRECISION, more=True),
        "termwise_399_derivative_kernel_bound_ball": termwise_derivative_kernel_upper.str(PRECISION, more=True),
        "grouped_any_contiguous_derivative_Green_kernel_bound_ball": grouped_derivative_kernel_upper.str(PRECISION, more=True),
        "derivative_termwise_to_grouped_improvement_lower_ball": derivative_gain.str(PRECISION, more=True),
        "abel_transfer": (
            "M(-r_j) is positive and nonincreasing in j. Abel summation bounds every weighted contiguous phasor "
            "by M at its left endpoint times the largest unweighted contiguous partial phasor."
        ),
        "green_transfer": (
            "K(r_j,s)=pi*M(-r_j)*M(-s)*sin(vartheta(r_j)-vartheta(s)); hence every contiguous lattice block "
            "is below B*(r_left*s)^(-1/4)<=B*r_min^(-1/2)."
        ),
        "exact_gauge_lattice_lock": (
            "For d_m=beta(4m/C-1), chi(d_m)=2beta^2*d_m+beta*d_m^2="
            "2pi*m^2-beta^3, so exp(-i chi(d_m))=exp(i beta^3) for every integer m."
        ),
        "derivative_green_transfer": (
            "partial_s K(r_j,s)=pi*M(-r_j)*N(-s)*sin(varphi(-s)-vartheta(r_j)); DLMF 9.8.21 "
            "and its signed remainder give N(-s)^2<=sqrt(s)/pi*(1+7/(32s^3))."
        ),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Grouped Airy phase and Green packets on the event lattice

Date: 2026-08-13

Status: rigorous finite all-height contiguous-packet and grouped Green-kernel
bound; no completed source or endpoint cancellation

Write the exact Airy event phase relative to event zero as

```text
vartheta(r_j)-vartheta(r_0)=2pi*j*(2q+j)+P_5(j,u)+delta_j(u),
|delta_j(u)|<{c['per_term_exact_Airy_normal_form_error_ball']}.       (GP1)
```

The integer phase disappears inside every exponential.  At `u=0`, direct
Arb summation of the 399 quintic phasors gives

```text
sum exp(iP_5(j,0))
 ={c['center_quintic_full_sum']['real_ball']}
 +i*{c['center_quintic_full_sum']['imag_ball']},
absolute value ={c['center_quintic_full_sum']['absolute_ball']}.       (GP2)
```

Enumeration of all `399*400/2` nonempty contiguous blocks proves

```text
max_(a<=b) |sum_(j=a)^b exp(iP_5(j,0))|
 <{c['center_quintic_max_contiguous_upper_ball']},
candidate block {c['center_quintic_max_contiguous_indices']}.         (GP3)
```

This is finite exact-roster enumeration with Arb balls, not floating-point
sampling.  Differentiating the explicit polynomial `P_5` in `u` and using
`0<=u<=1/(2C^2)` gives

```text
sum_j sup_u |P_5(j,u)-P_5(j,0)|
 <{c['uniform_sum_of_height_phase_variations_ball']}.                 (GP4)
```

The exponential Lipschitz inequality and (GP1)--(GP4) therefore prove,
uniformly on the complete top corridor,

```text
|sum_(j=-198)^200 exp(i[vartheta(r_j)-vartheta(r_0)])|
 >{c['exact_Airy_full_packet_absolute_lower_ball']},

sup_(a<=b) |sum_(j=a)^b exp(i[vartheta(r_j)-vartheta(r_0)])|
 <{c['exact_Airy_any_contiguous_packet_upper_ball']}<52.53.           (GP5)
```

The first line is a useful negative result: the cubic packet does not
collapse to an `O(1)` phasor.  The second line still gives substantial,
certified cancellation for every active contiguous roster block.

Using

```text
K(r_j,s)=pi M(-r_j)M(-s)
                 sin(vartheta(r_j)-vartheta(s)),                       (GP6)
```

DLMF 9.8(iii) makes `M(-r_j)` positive and nonincreasing in `j`.
Finite Abel summation, DLMF 9.8.20, and `r,s>=r_min` then give

```text
sup_(contiguous I) |sum_(j in I) K(r_j,s)|
 <{c['grouped_any_contiguous_Green_kernel_bound_ball']}<0.0246.       (GP7)
```

The corresponding termwise pointwise estimate is
`{c['termwise_399_kernel_bound_ball']}`; (GP7) improves it by a factor
greater than `{c['termwise_to_grouped_improvement_lower_ball']}` while
retaining phase cancellation.

Equation (GP7) is the kernel packet required after interchanging the finite
event sum with a variation-of-constants integral: at fixed `s`, the active
event indices are contiguous.

Returning from `H(r_m)` to the original transform introduces no modewise
gauge defect.  Indeed, on the exact event lattice,

```text
chi(d_m)=2beta^2 d_m+beta d_m^2
        =2pi m^2-beta^3,
exp[-i chi(d_m)]=exp(i beta^3).                      (GP8)
```

The last factor is common to every integer mode.  Thus the grouped packet
survives the inverse gauge exactly rather than approximately.

The initial-data channel uses

```text
partial_s K(r_j,s)=pi M(-r_j)N(-s)
                    sin(varphi(-s)-vartheta(r_j)).                    (GP9)
```

The event phasor is unchanged.  DLMF 9.8.21 and its signed remainder give

```text
N(-s)^2<=sqrt(s)/pi [1+7/(32s^3)].                  (GP10)
```

The same Abel argument therefore proves

```text
sup_(contiguous I) |sum_(j in I) partial_s K(r_j,s)|
 <{c['grouped_any_contiguous_derivative_Green_kernel_bound_ball']}<52.8. (GP11)
```

Its termwise counterpart is
`{c['termwise_399_derivative_kernel_bound_ball']}`, again improving by a
factor greater than
`{c['derivative_termwise_to_grouped_improvement_lower_ball']}`.

Equations (GP7) and (GP11) do not yet compose the initial data, coherent
endpoint sources, and interior forcing with the completed Poisson remainder.

Pi provenance: the phase integer and `P_5` inherit
`beta^3=pi C^2/8`; (GP6) uses the Airy Wronskian normalization.  The Airy
phase, modulus monotonicity, and modulus remainder are sourced from
`{DLMF_URL}`.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No completed endpoint cancellation, grouped 398-mode finite-integral splice,
all-corridor continuation, complete
`Q_K-T` or `T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    require(dependencies["cubic_resonance"].get("passed") is True, "cubic-resonance dependency failed")
    require(dependencies["airy_green_kernel"].get("passed") is True, "Airy Green-kernel dependency failed")

    ctx.dps = PRECISION
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "all_height_contiguous_Airy_phase_packets_and_grouped_Green_kernel_certified",
        "passed": True,
        "certificate": phase_packet_certificate(
            dependencies["cubic_resonance"], dependencies["airy_green_kernel"]
        ),
        "external_theorem": {
            "source": "NIST Digital Library of Mathematical Functions, Section 9.8",
            "url": DLMF_URL,
            "used_statements": [
                "Airy modulus and phase definitions 9.8.1--9.8.4",
                "Airy modulus monotonicity in Section 9.8(iii)",
                "Airy modulus expansion 9.8.20 and the signed remainder statement following 9.8.23",
                "Airy derivative-modulus expansion 9.8.21 and the signed remainder statement following 9.8.23",
            ],
        },
        "decision": {
            "uniform_any_contiguous_exact_Airy_phase_packet_below_52_53": True,
            "full_exact_Airy_phase_packet_above_31_4": True,
            "phase_only_O1_collapse_rejected": True,
            "grouped_any_contiguous_Green_kernel_below_0_0246": True,
            "inverse_gauge_phase_common_on_integer_event_lattice": True,
            "termwise_absolute_kernel_sum_used": False,
            "derivative_kernel_packet_proved": True,
            "completed_endpoint_cancellation_proved": False,
            "grouped_398_mode_splice_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)} for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
        "next_action": (
            "Insert both grouped packets into the exact variation-of-constants event sum. Keep the initial-data, "
            "endpoint-source, interior-forcing, and completed Poisson channels composed before taking a norm."
        ),
        "proof_boundary": (
            "Finite all-height exact Airy phase-packet and grouped K/partial_s K kernel bounds only. No completed "
            "endpoint cancellation, grouped finite-integral splice, complete Q_K-T or T_upper, Lambda<=0, "
            "PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("certified all-height Airy phase and grouped Green packets", flush=True)


if __name__ == "__main__":
    main()
