#!/usr/bin/env python3
"""Separate completed-point control from weighted exact-profile control."""

from __future__ import annotations

from decimal import Decimal, getcontext
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_completed_point_weighted_profile_obligation_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "fourier_completion": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_399_mode_fourier_completion_gate.json",
    "completed_height_cell": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_nonzero_height_completed_strip_cell_gate.json",
    "beta4_completed_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.json",
    "weighted_completion": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_completion_partial_sum_obligation_gate.json",
    "weighted_kernel": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_opposite_branch_dirichlet_gate.json",
    "grouped_selected": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_selected_branch_grouped_integral_gate.json",
}

CENTER = 39_895
LOWER_MODE = 40_093
UPPER_MODE = 40_094
LOWER_FREQUENCY = LOWER_MODE - CENTER
UPPER_FREQUENCY = UPPER_MODE - CENTER


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


def parse_ball(text: str) -> tuple[Decimal, Decimal]:
    body = text.strip().removeprefix("[").removesuffix("]").strip()
    if "+/-" not in body:
        return Decimal(body), Decimal(0)
    midpoint, radius = body.split("+/-", 1)
    return Decimal(midpoint.strip()), Decimal(radius.strip())


def positive_lower(ball: tuple[Decimal, Decimal]) -> Decimal:
    return ball[0] - ball[1]


def positive_upper(ball: tuple[Decimal, Decimal]) -> Decimal:
    return ball[0] + ball[1]


def exact_counterprofile_check() -> None:
    """Check Fourier/end-point algebra without analytic approximation."""

    z = Fraction(37, 41)
    coefficients = {
        LOWER_FREQUENCY: z,
        UPPER_FREQUENCY: -z,
    }
    require(LOWER_FREQUENCY == 198 and UPPER_FREQUENCY == 199, "counterprofile frequency drift")
    require(sum(coefficients.values(), Fraction()) == 0, "counterprofile changed s=0")
    # Every integer Fourier character is one at s=1, so the same exact sum
    # proves that the second endpoint is fixed as well.
    require(sum(coefficients.values(), Fraction()) == 0, "counterprofile changed s=1")

    weights = {
        LOWER_FREQUENCY: Fraction(11, 29),
        UPPER_FREQUENCY: Fraction(7, 31),
    }
    shift = sum((weights[n] * value for n, value in coefficients.items()), Fraction())
    require(shift == z * (weights[LOWER_FREQUENCY] - weights[UPPER_FREQUENCY]), "weighted shift drift")
    require(shift != 0, "exact counterprofile lost weighted distinction")


def numerical_certificate(dependencies: dict[str, Any]) -> dict[str, Any]:
    completion = dependencies["weighted_completion"]["coefficient_certificate"]
    kernel = dependencies["weighted_kernel"]["certificate"]
    height = dependencies["completed_height_cell"]["uniform_error"]
    geometry = dependencies["fourier_completion"]["geometry"]

    difference = parse_ball(completion["last_two_weight_difference_absolute_ball"])
    kernel_l1 = parse_ball(kernel["two_branch_kernel_L1_bound_ball"])
    completed_value = parse_ball(height["uniform_reduced_completed_strip_error_bound_ball"])
    spacing = parse_ball(geometry["detuning_spacing_h_ball"])

    require(positive_lower(difference) > Decimal("1.917e-7"), "last-weight distinction lost")
    require(positive_upper(kernel_l1) < Decimal("3.840283e-5"), "weighted-kernel L1 drift")
    require(positive_upper(completed_value) < Decimal("9.497e-12"), "completed scalar error drift")
    require(positive_lower(spacing) > Decimal("0.054"), "detuning spacing lost positivity")
    profile_scale_lower = Decimal(1) / positive_upper(spacing)
    profile_scale_upper = Decimal(1) / positive_lower(spacing)
    completed_profile_upper = positive_upper(completed_value) * profile_scale_upper
    raw_fold_transfer_upper = positive_upper(kernel_l1) * profile_scale_upper
    require(completed_profile_upper < Decimal("1.759e-10"), "canonical completed profile scale drift")
    require(raw_fold_transfer_upper < Decimal("7.112e-4"), "raw-fold transfer scale drift")

    return {
        "counterprofile_modes": [LOWER_MODE, UPPER_MODE],
        "counterprofile_frequencies": [LOWER_FREQUENCY, UPPER_FREQUENCY],
        "last_two_weight_difference_absolute_ball": completion["last_two_weight_difference_absolute_ball"],
        "last_two_weight_difference_lower_bound": str(positive_lower(difference)),
        "two_piece_weighted_kernel_L1_bound_ball": kernel["two_branch_kernel_L1_bound_ball"],
        "two_piece_weighted_kernel_L1_upper_bound": str(positive_upper(kernel_l1)),
        "uniform_completed_scalar_reduced_error_ball": height["uniform_reduced_completed_strip_error_bound_ball"],
        "uniform_completed_scalar_reduced_error_upper_bound": str(positive_upper(completed_value)),
        "detuning_spacing_h_ball": geometry["detuning_spacing_h_ball"],
        "canonical_profile_scale_identity": "Y/(2pi)=1/h",
        "canonical_profile_scale_interval": [str(profile_scale_lower), str(profile_scale_upper)],
        "uniform_completed_canonical_profile_value_upper_bound": str(completed_profile_upper),
        "raw_fold_profile_transfer_coefficient_upper_bound": str(raw_fold_transfer_upper),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Exact completed point versus weighted profile

Date: 2026-08-13

Status: exact logical separation and conditional whole-kernel transfer;
the physical exact-profile identity and uniform profile error remain open

Let `delta J_lambda(s)` denote the reduced exact-minus-beta-minus-four fold
integral before the all-mode evaluation, normalized so that the saved scalar
theorem is `delta J_lambda(0)`.  Since `hY=2pi`, the corresponding canonical
Fourier profile is

```text
delta g_lambda(s)=[Y/(2pi)]delta J_lambda(s)
                  =delta J_lambda(s)/h.                (EP0)
```

Use the selector-centred Fourier convention

```text
Delta G_(39895+n)(lambda)
 =Integral_0^1 delta g_lambda(s)e^(-2pi i n s)ds.       (EP1)
```

Here (EP1) is the common-profile representation that a finite-`t` remainder
theorem would have to establish; this gate does **not** assume that the still
missing physical amplitude/Fresnel remainder already has this form.  The
completed-strip theorems evaluate the all-mode object at `s=0`.  In the raw
reduced-fold and canonical Fourier normalizations they prove only

```text
|delta J_lambda(0)|
 <{c['uniform_completed_scalar_reduced_error_upper_bound']}<9.497e-12,

|delta g_lambda(0)|
 <{c['uniform_completed_canonical_profile_value_upper_bound']}<1.759e-10. (EP2)
```

Equation (EP2) does not control the weighted Fourier correction.  For an
arbitrary complex `z`, define the endpoint-preserving counterprofile

```text
phi_z(s)=z[e^(2pi i 198s)-e^(2pi i 199s)].             (EP3)
```

Both `phi_z(0)` and `phi_z(1)` vanish exactly.  Thus adding (EP3) leaves the
Fourier midpoint, endpoint half-current, their completed value, and the
event-zero coefficient unchanged.  Its only nonzero Fourier coefficients
in the selector roster are

```text
Delta G_40093=z,             Delta G_40094=-z.          (EP4)
```

The certified weights satisfy

```text
|w_40093-w_40094|
 ={c['last_two_weight_difference_absolute_ball']}
 >{c['last_two_weight_difference_lower_bound']}>1.917e-7.              (EP5)
```

Consequently (EP3) changes the weighted correction by
`z(w_40093-w_40094)` while preserving all data used by the completed-point
estimate.  This is a logical insufficiency result for the existing scalar
theorem, not a claim that the physical exact profile can be varied freely.

There is a cancellation-preserving sufficient replacement.  If the missing
physical correction is first proved to satisfy (EP1), define the whole
weighted kernel

```text
P_w(s)=sum_(n=-199)^199 w_(39895+n)e^(-2pi i n s),
w_39894=0.                                               (EP6)
```

The two independently enclosed pieces of this same kernel obey

```text
Integral_0^1|P_w(s)|ds
 <={c['two_piece_weighted_kernel_L1_upper_bound']}<3.840283e-5.          (EP7)
```

Therefore one uniform common-profile theorem

```text
sup_(lambda,s)|delta g_lambda(s)|<=epsilon_profile       (EP8)
```

would imply directly, with no coefficientwise estimate,

```text
|sum_m w_m Delta G_m|
 <=Integral_0^1|delta g_lambda(s)P_w(s)|ds
 <3.840283e-5 epsilon_profile.                           (EP9)
```

Equivalently, a raw-fold profile estimate
`sup|delta J_lambda|<=epsilon_J` would give the same weighted correction
below

```text
{c['raw_fold_profile_transfer_coefficient_upper_bound']} epsilon_J
 <7.112e-4 epsilon_J.                                   (EP10)
```

The next obligation is not 398 separate mode estimates.  It is to derive
the exact finite-`t` common profile before the all-mode collapse, prove that
the missing weighted physical remainder is its Fourier pairing with (EP6),
and bound that profile uniformly on the ordinary top corridor.  A completed
cutoff-family estimate is an equivalent admissible route.

Pi provenance: the `2pi` in (EP1), (EP3), and (EP6) is forced by the exact
selector relation `hY=2pi` and the integer Fourier--Poisson character.  This
gate introduces no fitted period, geometric circle, or polygon constant.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No exact physical common-profile identity, uniform exact-minus-beta-minus-four
profile bound, weighted finite-`t` amplitude/Fresnel remainder, completed
source/initial-data splice, complete `Q_K-T` or `T_upper`, all-corridor or
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion
is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    getcontext().prec = 90
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    for name, dependency in dependencies.items():
        require(dependency.get("passed") is True, f"dependency failed: {name}")

    exact_counterprofile_check()
    certificate = numerical_certificate(dependencies)
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "completed_point_insufficient_weighted_profile_transfer_obligation_certified",
        "passed": True,
        "exact_algebra": {
            "selector_center": CENTER,
            "profile_normalization": "delta g(s)=Y delta J(s)/(2pi)=delta J(s)/h because hY=2pi",
            "candidate_common_profile": "Delta G_(39895+n)=int_0^1 delta g(s) exp(-2pi i n s) ds",
            "completed_value": "The all-mode midpoint plus endpoint half-current evaluates delta J(0), hence delta g(0)=delta J(0)/h",
            "counterprofile": "phi_z(s)=z(exp(2pi i 198s)-exp(2pi i 199s))",
            "counterprofile_endpoints": "phi_z(0)=phi_z(1)=0",
            "counterprofile_coefficients": "Delta G_40093=z, Delta G_40094=-z",
            "counterprofile_weighted_shift": "z(w_40093-w_40094)",
            "conditional_transfer": "If ||delta g||_infinity<=epsilon_profile, then |sum w_m Delta G_m|<=epsilon_profile ||P_w||_1",
        },
        "certificate": certificate,
        "countermodel_scope": (
            "Logical non-identifiability from completed endpoint data alone; it does not assert that the physical exact profile "
            "admits arbitrary perturbations."
        ),
        "decision": {
            "completed_scalar_value_controls_weighted_profile_correction": False,
            "endpoint_half_current_and_midpoint_preserved_by_counterprofile": True,
            "event_zero_preserved_by_counterprofile": True,
            "conditional_whole_kernel_transfer_exact": True,
            "physical_exact_common_profile_operator_identity_proved": False,
            "uniform_exact_minus_beta4_profile_bound_proved": False,
            "weighted_exact_finite_integral_remainder_proved": False,
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
            "Derive the exact finite-t common profile before all-mode evaluation at s=0 and prove that the missing weighted "
            "amplitude/Fresnel correction is its Fourier pairing with P_w. Then certify a uniform profile or completed-cutoff "
            "bound and apply the whole-kernel L1 constant; do not estimate 398 modes separately."
        ),
        "proof_boundary": (
            "Exact completed-point insufficiency counterprofile and a conditional cancellation-preserving whole-kernel transfer only. "
            "No physical exact common-profile identity, uniform exact-minus-beta^-4 profile estimate, weighted finite-integral "
            "amplitude/Fresnel remainder, completed source/initial-data splice, complete Q_K-T or T_upper, all-corridor or "
            "height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("certified exact completed-point versus weighted-profile obligation", flush=True)


if __name__ == "__main__":
    main()
