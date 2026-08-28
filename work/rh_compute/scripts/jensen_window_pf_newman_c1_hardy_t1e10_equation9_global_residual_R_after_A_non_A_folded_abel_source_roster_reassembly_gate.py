#!/usr/bin/env python3
"""Certify the Abel-zero source-roster reassembly and transition schedule."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import time
from typing import Any


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_folded_abel_source_roster_reassembly_gate"
)
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

ABEL_STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_folded_source_kernel_abel_cubic_symmetry_gate"
)
COMMON_STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_common_kernel_theta_current_gate"
)
RECON_STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_"
    "double_weber_exact_source_reconstruction_gate"
)
PORTCULLIS_STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_"
    "finite_poisson_portcullis_saddle_reduction_gate"
)
COVERAGE_STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_"
    "lower_fold_ordinary_mode_coverage_ledger_gate"
)
REFLECTION_STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_"
    "half_kummer_reflection_branch_reduction_gate"
)
DEPENDENCIES = {
    "folded_Abel_limit": REPO_ROOT / f"work/rh_compute/results/{ABEL_STEM}.json",
    "non_A_common_kernel": REPO_ROOT / f"work/rh_compute/results/{COMMON_STEM}.json",
    "double_Weber_source_reconstruction": REPO_ROOT / f"work/rh_compute/results/{RECON_STEM}.json",
    "portcullis_saddle": REPO_ROOT / f"work/rh_compute/results/{PORTCULLIS_STEM}.json",
    "ordinary_fold_coverage": REPO_ROOT / f"work/rh_compute/results/{COVERAGE_STEM}.json",
    "half_reflection": REPO_ROOT / f"work/rh_compute/results/{REFLECTION_STEM}.json",
}

A = 159_577
B = 5_122_421
L = 2_481_422
FIRST = 622
ORDINARY_LAST = 39_852
A_WINDOW_FIRST = 39_853
TARGET_LAST = 39_894
A_WINDOW_LAST = 39_936
ORDINARY_CORE_LAST = 39_694
FOLD_FIRST = 39_695
HEIGHT = 10_000_000_000


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


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


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency {relative(path)}")
    return json.loads(path.read_text(encoding="utf-8"))


def fraction_record(value: Fraction) -> dict[str, Any]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "fraction": str(value),
        "decimal": format(float(value), ".17g"),
    }


def source_roster_certificate() -> dict[str, Any]:
    require(B - A == 2 * L, "source roster arithmetic drift")
    require(A % 2 == 1 and B % 2 == 1, "source endpoints are not odd")
    source_count = L + 1
    h_zero = (A + B) // 2
    f0_zero = L * (A + B - 2) // 2
    f1_zero = L * (A + B + 2) // 2
    abel_half_zero = (f0_zero + f1_zero) // 2
    source_zero = source_count * (A + B) // 2
    require(source_zero == h_zero + abel_half_zero, "x=0 endpoint reconstruction failed")

    # A finite coefficient replay: H contributes one half at each outer label,
    # while [F(0)+F(1)]/2 contributes the other halves and every interior label.
    surrogate_length = 7
    coefficients = [Fraction(0) for _ in range(surrogate_length + 1)]
    coefficients[0] += Fraction(1, 2)
    coefficients[-1] += Fraction(1, 2)
    for index in range(surrogate_length):
        coefficients[index] += Fraction(1, 2)
    for index in range(1, surrogate_length + 1):
        coefficients[index] += Fraction(1, 2)
    require(all(value == 1 for value in coefficients), "formal source coefficient replay failed")

    return {
        "source_odd_roster": [A, B],
        "source_coordinate": "f_x(n)=(A+2n)exp(i*pi*x*(A+2n)^2/4), 0<=n<=L",
        "source_count": source_count,
        "half_current": "H_x=[f_x(0)+f_x(L)]/2",
        "folded_endpoint_values": {
            "F_x(0)": "sum_(n=0)^(L-1) f_x(n)",
            "F_x(1)": "sum_(n=1)^L f_x(n)",
        },
        "endpoint_reassembly": (
            "H_x+[F_x(0)+F_x(1)]/2=sum_(n=0)^L f_x(n)"
        ),
        "surrogate_coefficient_replay": [str(value) for value in coefficients],
        "x_zero": {
            "H_0": h_zero,
            "F_0_0": f0_zero,
            "F_0_1": f1_zero,
            "Abel_endpoint_half_sum": abel_half_zero,
            "complete_source_sum": source_zero,
            "identity_checked": True,
        },
    }


def ownership_certificate() -> dict[str, Any]:
    low_count = ORDINARY_LAST - FIRST + 1
    high_count = A_WINDOW_LAST - A_WINDOW_FIRST + 1
    notch_count = A_WINDOW_LAST - FIRST + 1
    require((low_count, high_count, notch_count) == (39_231, 84, 39_315), "notch counts drift")

    # Coefficient order: P_m, A_m, A_-m, B_m.  The Abel source supplies -I_m,
    # where I_m=P_m+A_m+B_m.
    abel = (-1, -1, 0, -1)
    low_row = (0, 1, 0, 1)
    high_row = (0, 0, -1, 1)
    low_net = tuple(left + right for left, right in zip(abel, low_row))
    high_net = tuple(left + right for left, right in zip(abel, high_row))
    require(low_net == (-1, 0, 0, 0), "ordinary row reduction failed")
    require(high_net == (-1, -1, -1, 0), "A-window row reduction failed")

    partition = [
        {
            "sector": "ordinary_core",
            "mode_range": [FIRST, ORDINARY_CORE_LAST],
            "count": ORDINARY_CORE_LAST - FIRST + 1,
        },
        {
            "sector": "fold_owned_remaining_interior",
            "mode_range": [FOLD_FIRST, ORDINARY_LAST],
            "count": ORDINARY_LAST - FOLD_FIRST + 1,
        },
        {
            "sector": "extracted_A_transition_window",
            "mode_range": [A_WINDOW_FIRST, A_WINDOW_LAST],
            "count": high_count,
        },
    ]
    require(sum(row["count"] for row in partition) == notch_count, "route partition count failed")

    return {
        "extended_notch": [FIRST, A_WINDOW_LAST],
        "ordinary_endpoint_row": {
            "mode_range": [FIRST, ORDINARY_LAST],
            "count": low_count,
            "row": "A_m+B_m=I_m-P_m",
            "net_after_Abel_source": "-P_m",
        },
        "A_window_endpoint_row": {
            "mode_range": [A_WINDOW_FIRST, A_WINDOW_LAST],
            "count": high_count,
            "row": "B_m-A_-m=I_m-P_m-A_m-A_-m",
            "net_after_Abel_source": "-P_m-A_m-A_-m",
        },
        "coefficient_order": ["P_m", "A_m", "A_-m", "B_m"],
        "ordinary_net_vector": list(low_net),
        "A_window_net_vector": list(high_net),
        "exact_zero_regulator_source_block": (
            "sum_(n=0)^L f_x(n)-sum_(m=622)^39936 P_m(x)"
            "-sum_(m=39853)^39936[A_m(x)+A_-m(x)]"
        ),
        "non_A_reassembly": (
            "Inside the inherited ordered physical limit, replace H_x plus the folded Abel source plus "
            "its two endpoint rows by the exact zero-regulator source block; retain "
            "-chi_W*Btr-O-mathfrak_E_(A,42) unchanged."
        ),
        "post_A_positive_notch_route_partition": partition,
        "independent_norm_permission_created": False,
    }


def transition_certificate() -> dict[str, Any]:
    b_events = [(Fraction(2 * mode, B), "B_entry", mode) for mode in range(FIRST, A_WINDOW_LAST + 1)]
    a_events = [(Fraction(2 * mode, A), "A_exit", mode) for mode in range(FIRST, TARGET_LAST + 1)]
    events = sorted(b_events + a_events)
    require(len(events) == 78_588, "transition event count drift")
    require(len({value for value, _, _ in events}) == len(events), "coincident rational transitions found")

    minimum_gap: Fraction | None = None
    minimum_pair: tuple[tuple[Fraction, str, int], tuple[Fraction, str, int]] | None = None
    for left, right in zip(events, events[1:]):
        gap = right[0] - left[0]
        require(gap > 0, "transition ordering failed")
        if minimum_gap is None or gap < minimum_gap:
            minimum_gap = gap
            minimum_pair = (left, right)
    require(minimum_gap is not None and minimum_pair is not None, "empty transition roster")
    require(minimum_gap == Fraction(882, A * B), "minimum transition gap drift")
    require(
        minimum_pair[0][1:] == ("A_exit", 630)
        and minimum_pair[1][1:] == ("B_entry", 20_223),
        "minimum transition pair drift",
    )

    midpoint_exit_last = A // 4
    require(midpoint_exit_last == TARGET_LAST, "midpoint A-exit floor drift")
    require(4 * TARGET_LAST < A < 4 * (TARGET_LAST + 1), "midpoint boundary arithmetic failed")

    return {
        "finite_phase": "phi_m(y)=x*y^2/4-m*y",
        "finite_phase_derivative": "partial_y phi_m=x*y/2-m",
        "stationary_point": "y_*(m,x)=2m/x",
        "stationary_interval": "y_* in [A,B] iff x_B(m)=2m/B <= x <= x_A(m)=2m/A",
        "half_interval": ["0", "1/2"],
        "B_entry_events": {
            "mode_range": [FIRST, A_WINDOW_LAST],
            "count": len(b_events),
            "first": fraction_record(b_events[0][0]),
            "last": fraction_record(b_events[-1][0]),
        },
        "A_exit_events_on_half_interval": {
            "mode_range": [FIRST, TARGET_LAST],
            "count": len(a_events),
            "first": fraction_record(a_events[0][0]),
            "last": fraction_record(a_events[-1][0]),
        },
        "finite_y_stationary_modes_persisting_through_midpoint": {
            "mode_range": [TARGET_LAST + 1, A_WINDOW_LAST],
            "count": A_WINDOW_LAST - TARGET_LAST,
            "reason": "x_A(m)>1/2 iff 4m>A",
        },
        "total_distinct_events": len(events),
        "minimum_adjacent_gap": fraction_record(minimum_gap),
        "minimum_gap_pair": [
            {"kind": minimum_pair[0][1], "mode": minimum_pair[0][2], "x": str(minimum_pair[0][0])},
            {"kind": minimum_pair[1][1], "mode": minimum_pair[1][2], "x": str(minimum_pair[1][0])},
        ],
        "common_physical_phase": (
            "Psi_m(x)=-pi*m^2/x+(t/2)log((1-x)/x)"
        ),
        "physical_saddle": "x_m=2*pi*m^2/(t+2*pi*m^2)",
        "joint_source_coordinate": "alpha_m=2m+t/(pi*m)",
        "saddle_endpoint_equivalence": (
            "x_B(m)<=x_m<=x_A(m) iff A<=alpha_m<=B"
        ),
        "saved_height_classification": {
            "interior_joint_saddles": {"mode_range": [FIRST, ORDINARY_LAST], "count": 39_231},
            "lower_A_transition": {"mode_range": [A_WINDOW_FIRST, TARGET_LAST], "count": 42},
            "reflected_A_transition": {"mode_range": [TARGET_LAST + 1, A_WINDOW_LAST], "count": 42},
        },
    }


def render_note(artifact: dict[str, Any]) -> str:
    source = artifact["source_roster_certificate"]
    ownership = artifact["ownership_reassembly_certificate"]
    transitions = artifact["transition_schedule_certificate"]
    gap = transitions["minimum_adjacent_gap"]
    return f"""# Folded Abel source-roster reassembly

Date: 2026-08-26

Status: exact source/ownership reassembly and rational transition schedule;
the physical non-A bound remains open.

## Source reconstruction

Write

```text
f_x(n)=(A+2n)exp(i*pi*x*(A+2n)^2/4),       0<=n<=L,
A={A}, B={B}, L={L}.                                      (SR1)
```

The inherited endpoint half-current and the Abel boundary layer are

```text
H_x=[f_x(0)+f_x(L)]/2,
F_x(0)=sum_(n=0)^(L-1)f_x(n),
F_x(1)=sum_(n=1)^L f_x(n).                              (SR2)
```

Therefore, coefficient by coefficient,

```text
H_x+[F_x(0)+F_x(1)]/2=sum_(n=0)^L f_x(n).              (SR3)
```

This is the complete odd source roster `A,A+2,...,B`, containing
`{source['source_count']}` labels.  At `x=0`, (SR3) is the exact integer

```text
{source['x_zero']['complete_source_sum']}.             (SR4)
```

No split `1/x` formula is needed at that endpoint.

## Ownership reduction

The Abel-limit source contributes `-I_m` on the extended notch.  On
`622..39852`, the endpoint row is `A_m+B_m=I_m-P_m`; on
`39853..39936`, it is
`B_m-A_-m=I_m-P_m-A_m-A_-m`.  Combining each row before a norm gives

```text
H_x+AbelSource+EndpointRows
 =sum_(n=0)^L f_x(n)
  -sum_(m=622)^39936 P_m(x)
  -sum_(m=39853)^39936 [A_m(x)+A_-m(x)].              (SR5)
```

Thus the 39,315 raw finite Fresnel coefficients are a valid independent
cross-check, but they are not the primary physical quadrature: every one of
their owned contributions cancels algebraically to the already defined bulk
and A-face atoms.  The B-trace, B-outer, and translated A-face subtractions
remain unchanged in the common non-A limit.

The post-A positive notch has the exact route partition

```text
622..39694:   39073 ordinary-core modes,
39695..39852:   158 remaining fold-owned interior modes,
39853..39936:    84 extracted A-transition modes.      (SR6)
```

This partition selects evaluator coordinates; it does not license separate
absolute bounds.

## Transition schedule

For `phi_m(y)=x*y^2/4-m*y`, the finite coefficient has stationary point
`y_*=2m/x`.  It enters at `x_B(m)=2m/B` and exits at `x_A(m)=2m/A`.
On `0<=x<=1/2` there are `{transitions['total_distinct_events']}` distinct
rational events.  The closest pair is

```text
x_A(630) < x_B(20223),
x_B(20223)-x_A(630)={gap['fraction']}
                    ={gap['decimal']}... .             (SR7)
```

All modes `39895..39936` have `x_A(m)>1/2`, so their finite-`y` stationary
points persist through the midpoint.  Their joint physical saddles lie in
the reflected `x>1/2` half; exact Kummer reflection supplies the full source
without introducing a second half-domain stationary family.  After quadratic
completion, the common physical phase and saddle are

```text
Psi_m(x)=-pi*m^2/x+(t/2)log((1-x)/x),
x_m=2*pi*m^2/(t+2*pi*m^2),
alpha_m=2m+t/(pi*m),                                   (SR8)
```

with `x_B<=x_m<=x_A` exactly when `A<=alpha_m<=B`.  At the saved height this
reproduces 39,231 interior modes, 42 lower A-transition modes, and 42
reflected A-transition modes.

## Route decision

Use the exact double-Weber source reconstruction, the closed full-line bulk,
the remaining 158 fold modes, and the already owned A/B channels in one
signed physical assembly.  Use the endpoint-plus-39,315-Fresnel form only as
an altered-representation check at sparse points.  The next numerical gate is
a bounded source-minus-owned-carriers pilot split at the rational events in
(SR7), followed by interval `x` quadrature if the signed scale is viable.

## Pi provenance

Every `pi` in (SR1)--(SR8) is inherited from the Gaussian Abel/Fourier
kernel, the original Kummer quadratic phase, or the physical Kummer weight.
The rational event locations themselves contain no `pi`; the `pi` in the
physical saddle is forced by differentiating the declared phase.

## Boundary

{artifact['proof_boundary']}
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    require(CHECKER.is_file(), "missing independent checker")
    loaded = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in loaded.values()), "a dependency gate failed")
    require(
        loaded["folded_Abel_limit"]["decision"][
            "complete_cell_Abel_limit_reduced_to_endpoint_half_sum_and_finite_modes"
        ] is True,
        "Abel endpoint dependency drift",
    )
    require(
        loaded["non_A_common_kernel"]["decision"]["extended_notch_and_two_jet_forms_equivalent"] is True,
        "non-A ownership dependency drift",
    )
    require(
        loaded["double_Weber_source_reconstruction"]["decision"][
            "entire_finite_equation9_source_roster_reconstructed_exactly"
        ] is True,
        "source reconstruction dependency drift",
    )
    stationary = loaded["portcullis_saddle"]["stationary_classification"]
    require(stationary["lower_branch_interior_modes"] == "622..39852", "interior roster drift")
    require(stationary["turning_gap_modes"] == "39853..39936", "A-window roster drift")
    coverage = loaded["ordinary_fold_coverage"]["mode_coverage"]
    require(
        coverage["ordinary_atlas_validity_overlap"]["range"] == [39_695, 39_852],
        "fold overlap drift",
    )
    require(
        loaded["half_reflection"]["decision"]["odd_alpha_half_kummer_reflection_proved"] is True,
        "half reflection dependency drift",
    )

    source = source_roster_certificate()
    ownership = ownership_certificate()
    transitions = transition_certificate()
    proof_boundary = (
        "The exact Abel-zero reassembly of the inherited half-current, folded endpoint boundary layer, "
        "and two post-A endpoint rows into the complete 2481423-label odd source minus the contiguous "
        "P block and 84 paired A atoms; the exact x=0 source value; and the complete rational finite-phase "
        "transition schedule on 0<=x<=1/2 only. The route partition is not permission for independent "
        "norms. No numerical signed assembly, fast-evaluator error theorem, physical x quadrature, non-A "
        "bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or "
        "prize-level conclusion is proved."
    )
    artifact = {
        "kind": STEM,
        "status": (
            "exact_folded_Abel_source_roster_and_owned_carrier_reassembly_certified_"
            "source_minus_carriers_quadrature_open"
        ),
        "passed": True,
        "scope": {
            "height": HEIGHT,
            "source_odd_roster": [A, B],
            "source_coordinate_length": L,
            "extended_notch": [FIRST, A_WINDOW_LAST],
            "half_x_interval": ["0", "1/2"],
        },
        "source_roster_certificate": source,
        "ownership_reassembly_certificate": ownership,
        "transition_schedule_certificate": transitions,
        "decision": {
            "half_current_plus_Abel_boundary_layer_equals_complete_source_roster": True,
            "all_39315_owned_finite_mode_contributions_reassembled_before_norms": True,
            "zero_regulator_source_block_matches_post_A_ownership": True,
            "raw_39315_Fresnel_physical_quadrature_selected": False,
            "raw_39315_Fresnel_form_retained_as_cross_check": True,
            "source_roster_minus_owned_bulk_and_A_face_route_selected": True,
            "complete_rational_transition_schedule_certified": True,
            "physical_x_quadrature_completed": False,
            "non_A_bound_proved": False,
            "rh_implication": False,
        },
        "next_obligation": (
            "Build one bounded altered-representation pilot for the signed complete-source-minus-owned-"
            "carriers assembly on 0<=x<=1/2. Evaluate the complete source with the validated O(L) oracle "
            "or a separately checked recursive Gauss current, subtract the closed P block and joined A/B "
            "channels, and cross-check sparse rational x points against the endpoint-plus-39315-Fresnel "
            "form. Split at the exact rational transition roster and promote only an interval physical "
            "quadrature whose cumulative error is below the non-A allowance."
        ),
        "proof_boundary": proof_boundary,
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "threads": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.time() - started, 3),
        },
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("certified folded Abel source-roster reassembly and transition schedule", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
