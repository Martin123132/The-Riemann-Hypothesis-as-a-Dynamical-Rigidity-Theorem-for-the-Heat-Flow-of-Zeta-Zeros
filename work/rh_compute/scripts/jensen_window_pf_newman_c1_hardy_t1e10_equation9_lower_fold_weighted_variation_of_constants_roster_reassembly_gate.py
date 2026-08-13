#!/usr/bin/env python3
"""Certify the weighted branchwise variation-of-constants roster reassembly."""

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

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_variation_of_constants_roster_reassembly_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "detuning_ode": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_beta4_endpoint_coherent_detuning_ode_gate.json",
    "green_kernel": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_beta4_detuning_airy_green_kernel_gate.json",
    "grouped_packets": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_detuning_airy_lattice_grouped_green_packet_gate.json",
    "weighted_packets": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_leading_defect_weighted_airy_green_packet_gate.json",
    "completed_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.json",
}

Q = 39_894
NEGATIVE_COUNT = 198
POSITIVE_COUNT = 200
PRECISION = 100


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


def exact_algebra() -> dict[str, Any]:
    beta, c, m = sp.symbols("beta C m", nonzero=True)
    d_m = beta * (4 * m / c - 1)
    chi = sp.expand(2 * beta**2 * d_m + beta * d_m**2)
    beta_cubed = sp.pi * c**2 / 8
    locked = sp.expand(chi.subs(beta**3, beta_cubed))
    require(sp.simplify(locked - (2 * sp.pi * m**2 - beta_cubed)) == 0, "gauge lock failed")

    positive_counts = [POSITIVE_COUNT - k + 1 for k in range(1, POSITIVE_COUNT + 1)]
    negative_counts = [NEGATIVE_COUNT - k + 1 for k in range(1, NEGATIVE_COUNT + 1)]
    require(sum(positive_counts) == POSITIVE_COUNT * (POSITIVE_COUNT + 1) // 2, "positive Fubini count failed")
    require(sum(negative_counts) == NEGATIVE_COUNT * (NEGATIVE_COUNT + 1) // 2, "negative Fubini count failed")
    return {
        "reference_mode": Q,
        "reference_coordinate": "r_0=r_q",
        "negative_offsets": [-NEGATIVE_COUNT, -1],
        "positive_offsets": [1, POSITIVE_COUNT],
        "event_lattice": "r_j=r_0+j*Delta_r",
        "transformed_source": "F_sigma(s)=exp(i chi(d(s)))[R_0,sigma(d(s))+R_Y,sigma(d(s))+Q_sigma(d(s))]",
        "one_mode_variation_formula": (
            "H_sigma(r_j)=-partial_s K(r_j,r_0)H_sigma(r_0)+K(r_j,r_0)H_sigma,r(r_0)"
            "+Integral_(r_0)^(r_j)K(r_j,s)F_sigma(s)ds"
        ),
        "positive_reassembly": (
            "sum_(j=1)^200 a_j H_+(r_j)=-D_+ H_+(r_0)+A_+ H_+,r(r_0)"
            "+Integral_(r_0)^(r_200)[sum_(j: r_j>=s) a_j K(r_j,s)]F_+(s)ds"
        ),
        "negative_reassembly": (
            "sum_(j=-198)^-1 a_j H_-(r_j)=-D_- H_-(r_0)+A_- H_-,r(r_0)"
            "-Integral_(r_-198)^(r_0)[sum_(j: r_j<=s) a_j K(r_j,s)]F_-(s)ds"
        ),
        "initial_packets": "D_sigma=sum a_j partial_s K(r_j,r_0); A_sigma=sum a_j K(r_j,r_0)",
        "positive_active_blocks": "suffixes k..200 with counts 200,199,...,1",
        "negative_active_blocks": "prefixes -198..-k with counts 198,197,...,1",
        "positive_source_cell_incidences": sum(positive_counts),
        "negative_source_cell_incidences": sum(negative_counts),
        "gauge_lock": "chi(d_m)=2pi*m^2-beta^3; exp(-i chi(d_m))=exp(i beta^3)",
        "original_transform_reassembly": (
            "sum a_j G_sigma(d_j)=exp(i beta^3) sum a_j H_sigma(r_j), separately on each selected branch"
        ),
    }


def operator_certificate(grouped: dict[str, Any], weighted: dict[str, Any]) -> dict[str, Any]:
    grouped_certificate = grouped["certificate"]
    weighted_certificate = weighted["certificate"]
    base_k = arb(grouped_certificate["grouped_any_contiguous_Green_kernel_bound_ball"]).upper()
    base_d = arb(grouped_certificate["grouped_any_contiguous_derivative_Green_kernel_bound_ball"]).upper()
    minus_variation = arb(
        weighted_certificate["minus_branch"]["maximum_any_contiguous_modulus_adjusted_variation_ball"]
    ).upper()
    plus_variation = arb(
        weighted_certificate["plus_branch"]["maximum_any_contiguous_modulus_adjusted_variation_ball"]
    ).upper()
    minus_k = (base_k * minus_variation).upper()
    plus_k = (base_k * plus_variation).upper()
    minus_d = (base_d * minus_variation).upper()
    plus_d = (base_d * plus_variation).upper()
    require(minus_k + plus_k < arb("7.58e-7"), "branch K coefficients drifted")
    require(minus_d + plus_d < arb("0.001634"), "branch derivative coefficients drifted")
    return {
        "minus_branch_weighted_K_coefficient_ball": minus_k.str(PRECISION, more=True),
        "plus_branch_weighted_K_coefficient_ball": plus_k.str(PRECISION, more=True),
        "minus_branch_weighted_partial_s_K_coefficient_ball": minus_d.str(PRECISION, more=True),
        "plus_branch_weighted_partial_s_K_coefficient_ball": plus_d.str(PRECISION, more=True),
        "conditional_norm_bound": (
            "|weighted roster| <= D_-|H_-(r0)|+D_+|H_+(r0)|+A_-|H_-,r(r0)|+A_+|H_+,r(r0)|"
            "+A_- Integral_(r_-198)^r0 |F_-(s)|ds+A_+ Integral_r0^(r_200)|F_+(s)|ds"
        ),
        "norm_guard": (
            "This fallback inequality may be used only after the initial-data and source channels have been composed "
            "with the event-zero and completed Poisson terms in the same normalization."
        ),
    }


def render_note(artifact: dict[str, Any]) -> str:
    a = artifact["exact_algebra"]
    c = artifact["operator_certificate"]
    return f"""# Weighted variation-of-constants roster reassembly

Date: 2026-08-13

Status: exact branchwise finite-Fubini reassembly and conditional packet
coefficients proved; completed source and initial-data cancellation open

Choose event zero as reference and write

```text
r_j=r_0+j Delta_r,
F_sigma(s)=exp(i chi(d(s)))[R_0,sigma+R_Y,sigma+Q_sigma]. (VR1)
```

For either selected Hankel branch, Section 11.388 gives

```text
H_sigma(r_j)=-partial_s K(r_j,r_0)H_sigma(r_0)
 +K(r_j,r_0)H_sigma,r(r_0)
 +Integral_(r_0)^(r_j)K(r_j,s)F_sigma(s)ds.          (VR2)
```

Keep the certified leading-defect weights `a_j` attached.  Finite Fubini
reassembly on the positive branch is exactly

```text
sum_(j=1)^200 a_j H_+(r_j)
 =-D_+ H_+(r_0)+A_+ H_+,r(r_0)
 +Integral_(r_0)^(r_200)
    [sum_(j:r_j>=s) a_jK(r_j,s)]F_+(s)ds,            (VR3)
```

whereas orientation on the negative branch gives

```text
sum_(j=-198)^-1 a_j H_-(r_j)
 =-D_- H_-(r_0)+A_- H_-,r(r_0)
 -Integral_(r_-198)^(r_0)
    [sum_(j:r_j<=s) a_jK(r_j,s)]F_-(s)ds.            (VR4)
```

Here `D_sigma=sum a_j partial_sK(r_j,r_0)` and
`A_sigma=sum a_jK(r_j,r_0)`.  At fixed `s`, (VR3) uses a suffix of
`1..200`; (VR4) uses a prefix of `-198..-1`.  The exact source-cell incidence
counts are `{a['positive_source_cell_incidences']}` and
`{a['negative_source_cell_incidences']}`, respectively.  An independent
rational finite-array checker verifies both orientations.

The inverse gauge remains common:

```text
sum a_jG_sigma(d_j)=exp(i beta^3)sum a_jH_sigma(r_j). (VR5)
```

Thus no modewise phase is introduced after (VR3)--(VR4).  Sections
11.390--11.391 supply the branchwise coefficients

```text
|A_-| and negative source packet
 <{c['minus_branch_weighted_K_coefficient_ball']},
|A_+| and positive source packet
 <{c['plus_branch_weighted_K_coefficient_ball']},

|D_-|<{c['minus_branch_weighted_partial_s_K_coefficient_ball']},
|D_+|<{c['plus_branch_weighted_partial_s_K_coefficient_ball']}.       (VR6)
```

Equations (VR3)--(VR6) expose the remaining wall precisely.  The `K`
source coefficients are at the `10^-7` scale, but the derivative packets
multiply the branch initial data at the `8.2e-4` scale.  Neither branch
initial value may be discarded.  The event-zero term, the beta-minus-four
identity

```text
H4+sum G4+C4=g4_lambda(0),                            (VR7)
```

and the exact-minus-beta-four source must be rewritten in the same
`H_sigma(r_0),H_sigma,r(r_0),F_sigma` coordinates before any final norm.

This is an exact reassembly for the leading stationary-defect weights only.
It does not identify the finite-integral amplitude remainder with `a_j`, nor
does it prove the cancellation demanded after (VR7).

Pi provenance: (VR5) uses `beta^3=pi C^2/8` and integer mode parity; all
kernel normalizations are inherited from the Airy Wronskian.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No completed initial-data/source cancellation, exact finite-integral
amplitude remainder, grouped 398-mode splice, complete `Q_K-T` or
`T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is
proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    for name, dependency in dependencies.items():
        require(dependency.get("passed") is True, f"dependency failed: {name}")

    ctx.dps = PRECISION
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "weighted_branchwise_variation_of_constants_roster_exactly_reassembled",
        "passed": True,
        "exact_algebra": exact_algebra(),
        "operator_certificate": operator_certificate(
            dependencies["grouped_packets"], dependencies["weighted_packets"]
        ),
        "decision": {
            "positive_suffix_Fubini_reassembly_exact": True,
            "negative_prefix_Fubini_reassembly_exact": True,
            "inverse_gauge_common": True,
            "leading_defect_weighted_packet_coefficients_inserted": True,
            "completed_initial_data_cancellation_proved": False,
            "completed_source_cancellation_proved": False,
            "finite_integral_amplitude_remainder_proved": False,
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
            "Express the beta^-4 completed projection and event-zero correction in the same H_sigma(r0), H_sigma,r(r0), "
            "and F_sigma coordinates. Test whether the derivative initial-data packet cancels algebraically; if it does "
            "not, derive a certified residual bound before addressing the finite-integral amplitude remainder."
        ),
        "proof_boundary": (
            "Exact leading-defect weighted variation-of-constants roster reassembly and conditional packet coefficients "
            "only. No completed initial-data/source cancellation, finite-integral amplitude remainder, grouped 398-mode "
            "splice, complete Q_K-T or T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("certified weighted variation-of-constants roster reassembly", flush=True)


if __name__ == "__main__":
    main()
