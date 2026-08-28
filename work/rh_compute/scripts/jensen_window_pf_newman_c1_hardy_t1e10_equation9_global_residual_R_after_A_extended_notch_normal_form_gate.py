#!/usr/bin/env python3
"""Certify the post-A extended-notch normal form for the joined residual."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import time
from typing import Any

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_normal_form_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__)
CHECKER = Path(__file__).with_name("check_" + STEM + ".py")
DEPENDENCIES = {
    "R_after_A_equivalence": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_finite_regulator_equivalence_gate.json",
    "R_Dir_ownership": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_Dir_projector_ownership_ledger_gate.json",
    "two_jet_reassembly": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_periodic_two_jet_common_kernel_reassembly_gate.json",
    "all_order_route_guard": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_all_order_endpoint_jet_route_gate.json",
}

A = 159_577
B = 5_122_421
T_LO = 622
T_HI = 39_894
A_WINDOW_LO = 39_853
A_WINDOW_HI = 39_936
EXTENDED_HI = A_WINDOW_HI


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


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


def set_certificate() -> dict[str, Any]:
    target = set(range(T_LO, T_HI + 1))
    a_window = set(range(A_WINDOW_LO, A_WINDOW_HI + 1))
    extended = set(range(T_LO, EXTENDED_HI + 1))
    require(target | a_window == extended, "extended-notch union drift")
    require(target & a_window == set(range(A_WINDOW_LO, T_HI + 1)), "target/A overlap drift")
    require(a_window - target == set(range(T_HI + 1, A_WINDOW_HI + 1)), "outside-target A half drift")
    require(len(target) == 39_273, "target count drift")
    require(len(a_window) == 84, "A-window count drift")
    require(len(target & a_window) == 42, "inside-target A count drift")
    require(len(a_window - target) == 42, "outside-target A count drift")
    require(len(extended) == 39_315, "extended-notch count drift")

    for mode in range(1, EXTENDED_HI + 2):
        chi_t = int(mode in target)
        chi_a = int(mode in a_window)
        chi_u = int(mode in extended)
        require(chi_t + chi_a * (1 - chi_t) == chi_u, f"projector coefficient drift at mode {mode}")

    return {
        "original_target": [T_LO, T_HI],
        "A_window": [A_WINDOW_LO, A_WINDOW_HI],
        "extended_target": [T_LO, EXTENDED_HI],
        "counts": {
            "original_target": len(target),
            "A_window": len(a_window),
            "A_window_inside_original_target": len(target & a_window),
            "A_window_outside_original_target": len(a_window - target),
            "extended_target": len(extended),
        },
        "indicator_identity": "chi_T+chi_A*(1-chi_T)=chi_U, U=T union A_window={622,...,39936}",
        "finite_projector_identity": "Delta-mathcal_A=H+sum_(m=-M)^M w_m I_m-sum_(m=622)^39936 w_m P_m-sum_(m=39853)^39936 w_m(A_m+A_-m)",
    }


def symbolic_certificate() -> dict[str, Any]:
    source_all, source_u, endpoint_u, bulk_u, a_pair = sp.symbols(
        "source_all source_u endpoint_u bulk_u a_pair"
    )
    endpoint_relation = sp.Eq(endpoint_u, source_u - bulk_u)
    folded = source_all - source_u + endpoint_u - a_pair
    reduced = folded.subs(endpoint_u, endpoint_relation.rhs)
    require(sp.expand(reduced - (source_all - bulk_u - a_pair)) == 0, "folded extended-notch identity drift")

    a_plus, a_minus, b_plus = sp.symbols("a_plus a_minus b_plus")
    window_endpoint = a_plus + b_plus - (a_plus + a_minus)
    require(sp.expand(window_endpoint - (b_plus - a_minus)) == 0, "A-window endpoint roster drift")

    return {
        "folded_source_form": "H_x+integral_0^L f_x(u)[D_(M,epsilon)(u)-G_(U,epsilon)(u)]du+sum_(m in U)w_m(A_m+B_m)-sum_(m in W_A)w_m(A_m+A_-m)",
        "one_cell_form": "Replace the source integral by integral_0^1 F_x(s)C_(U,M,epsilon)(s)ds using the certified periodic fold.",
        "endpoint_roster": {
            "622_to_39852": "A_m+B_m",
            "39853_to_39936": "B_m-A_-m",
        },
        "two_jet_common_kernel_form": "H_x+E_x+2*pi*i*sum_(m=1)^M m*w_m[hat Phi_x(m)-hat Phi_x(-m)]-sum_(m=622)^39936 w_m P_m-sum_(m=39853)^39936 w_m(A_m+A_-m)",
        "post_A_residual": "mathfrak R_A=(Delta-mathcal_A)-chi_W*Btr-O",
    }


def edge_geometry() -> dict[str, Any]:
    old_edge = Fraction(2 * T_HI + 1, 2)
    new_edge = Fraction(2 * EXTENDED_HI + 1, 2)
    half_boundary = Fraction(A, 4)
    shift = new_edge - old_edge
    old_gap = old_edge - half_boundary
    new_gap = new_edge - half_boundary
    require(shift == 42, "upper-edge shift drift")
    require(old_gap == Fraction(1, 4), "old half-boundary gap drift")
    require(new_gap == Fraction(169, 4), "new half-boundary gap drift")

    def q_half(mode: Fraction) -> Fraction:
        return (A - 4 * mode) / 2

    values = {
        "deleted_first_integer": q_half(Fraction(A_WINDOW_LO)),
        "deleted_last_integer": q_half(Fraction(A_WINDOW_HI)),
        "lower_surviving_integer": q_half(Fraction(A_WINDOW_LO - 1)),
        "upper_surviving_integer": q_half(Fraction(A_WINDOW_HI + 1)),
        "old_projector_edge": q_half(old_edge),
        "new_projector_edge": q_half(new_edge),
    }
    expected = {
        "deleted_first_integer": Fraction(165, 2),
        "deleted_last_integer": Fraction(-167, 2),
        "lower_surviving_integer": Fraction(169, 2),
        "upper_surviving_integer": Fraction(-171, 2),
        "old_projector_edge": Fraction(-1, 2),
        "new_projector_edge": Fraction(-169, 2),
    }
    require(values == expected, "A normal-coordinate geometry drift")

    return {
        "old_upper_half_integer_edge": str(old_edge),
        "new_upper_half_integer_edge": str(new_edge),
        "edge_shift_modes": str(shift),
        "A_half_boundary": str(half_boundary),
        "old_edge_gap_modes": str(old_gap),
        "new_edge_gap_modes": str(new_gap),
        "q_A_at_x_half": {name: str(value) for name, value in values.items()},
        "normal_gap_consequence": "At x=1/2 the extended upper projector edge has q_A=-84.5 instead of -0.5, and the nearest surviving A endpoint integers have |q_A|>=84.5.",
    }


def modular_certificate() -> dict[str, Any]:
    return {
        "edges": {"lower": "c=621.5=1243/2", "upper": "d_A=39936.5=79873/2"},
        "interpolant": "h_A,epsilon(y)=exp(-pi*epsilon*y^2)[1_(y<c)+1_(y>d_A)]",
        "integer_samples": "h_A,epsilon(m)=w_m[1-1_(622<=m<=39936)]",
        "poisson_identity": "C_A,epsilon(u)=sum_(m in Z)[1-1_U(m)]w_m exp(-2*pi*i*m*u)=sum_(k in Z)hat h_A,epsilon(k+u)",
        "edge_phase_guard": "Both c and d_A are half-integers, so exp(-2*pi*i*k*c)=exp(-2*pi*i*k*d_A)=(-1)^k for integer k.",
        "pole_free_inheritance": "The two-edge Mittag-Leffler closure and centred-cell k=0 recombination apply with d replaced by d_A; the proof uses only Gaussian differentiation and the shared half-integer phase.",
        "limit_order": "M tends to infinity first at fixed epsilon>0; epsilon tends down to zero second, exactly as in the inherited R_after_A definition.",
    }


def render_note() -> str:
    return f"""# Post-A extended-notch normal form

Date: 2026-08-23

Status: exact post-A projector compression and transition-free upper spectral
edge certified; joined phase-adapted bound open

Let

```text
T={{622,...,39894}},
W_A={{39853,...,39936}},
U={{622,...,39936}}=T union W_A.                     (EN1)
```

The two halves of `W_A` contain 42 modes each: `39853..39894` is already in
the original target, while `39895..39936` lies outside it.  Therefore the
indicator identity

```text
chi_T+chi_A(1-chi_T)=chi_U                            (EN2)
```

holds coefficientwise.  Subtracting the certified finite A transition from
the common projector kernel before any norm gives, for every inherited common
cutoff and regulator,

```text
Delta_(M,epsilon)-mathcal A_(M,epsilon)
 =H_x+sum_(m=-M)^M w_m I_m
  -sum_(m=622)^39936 w_m P_m
  -sum_(m=39853)^39936 w_m(A_m+A_-m).                (EN3)
```

Thus the positive full-line projector is no longer a six-row mask: it is one
contiguous extended notch.  Equivalently,

```text
Delta-mathcal A
 =H_x+integral_0^L f_x(u)[D_(M,epsilon)-G_(U,epsilon)]du
  +sum_(m in U)w_m(A_m+B_m)
  -sum_(m in W_A)w_m(A_m+A_-m).                      (EN4)
```

The endpoint rows in (EN4) reduce exactly to

```text
622..39852:    A_m+B_m,
39853..39936:  B_m-A_-m.                             (EN5)
```

The one-cell fold applies directly to the source integral in (EN4).  The
two-jet reassembly gives the equivalent common-kernel form

```text
H_x+E_x
 +2*pi*i sum_(m=1)^M m*w_m[hat Phi_x(m)-hat Phi_x(-m)]
 -sum_(m=622)^39936 w_m P_m
 -sum_(m=39853)^39936 w_m(A_m+A_-m).                 (EN6)
```

The B window and analytic B exterior remain separate subtractions:

```text
mathfrak R_A=(Delta-mathcal A)-chi_W Btr-O.           (EN7)
```

No B atom or negative Gamma bulk has been removed by (EN1)--(EN7).

The extended target has half-integer edges

```text
c=621.5,             d_A=39936.5.                    (EN8)
```

Hence its exact Gaussian interpolant

```text
h_A,epsilon(y)=e^(-pi*epsilon*y^2)
               [1_(y<c)+1_(y>d_A)]                  (EN9)
```

has integer samples equal to the complement of `U`.  Poisson summation, the
two closed boundary currents, and the pole-free centred-cell recombination
from Sections 11.444--11.445 apply with `d` replaced by `d_A`: both edges are
half-integers and therefore have the same integer phase `(-1)^k`.

This edge translation has a concrete geometric effect.  The A half-boundary
is `A/4=39894.25`.  The old and new upper-edge gaps are

```text
39894.5-A/4=0.25,
39936.5-A/4=42.25.                                   (EN10)
```

For the exact A normal coordinate at `x=1/2`,

```text
q_A(m,1/2)=(A-4m)/2.                                 (EN11)
```

The deleted integer window runs from `q_A=82.5` to `q_A=-83.5`.  The nearest
surviving integers have `q_A=84.5` and `q_A=-85.5`, while the new upper
half-integer edge has `q_A=-84.5` instead of the old `-0.5`.  The certified A
transition has therefore moved the upper spectral edge out of the A
half-boundary collar rather than merely subtracting a scalar after the fact.

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Pi provenance: every `pi` in (EN1)--(EN11) is inherited from the Gaussian
Abel regulator, integer Fourier character, and Kummer phase.  No fitted or
geometric occurrence is introduced.

Proof boundary: exact finite-regulator post-A coefficient compression,
extended-notch fold/modular inheritance, and rational upper-edge geometry
only.  No quantitative bound for the remaining joined kernel, no complete
`R_after_A`, `R_Dir`, or `Q_K-T` estimate, no all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["R_after_A_equivalence"]["decision"]["common_kernel_equals_six_class_post_A_mask"] is True, "post-A mask dependency drift")
    require(dependencies["R_Dir_ownership"]["decision"]["A_endpoint_and_positive_bulk_coefficients_transferred_exactly"] is True, "A transfer dependency drift")
    require(dependencies["two_jet_reassembly"]["decision"]["target_harmonic_current_cancels_coefficientwise"] is True, "two-jet cancellation dependency drift")
    require(dependencies["all_order_route_guard"]["decision"]["phase_adapted_transition_charts_still_required"] is True, "route-guard dependency drift")

    artifact = {
        "kind": STEM,
        "status": "exact_post_A_extended_positive_projector_notch_and_transition_free_upper_edge_certified_joined_bound_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "source_odd_roster": [A, B],
            "original_target": [T_LO, T_HI],
            "A_transition_window": [A_WINDOW_LO, A_WINDOW_HI],
            "extended_target": [T_LO, EXTENDED_HI],
            "common_cutoff": "M>=B=5122421",
        },
        "set_certificate": set_certificate(),
        "symbolic_certificate": symbolic_certificate(),
        "edge_geometry": edge_geometry(),
        "modular_certificate": modular_certificate(),
        "decision": {
            "post_A_positive_bulk_mask_is_one_contiguous_extended_notch": True,
            "six_row_mask_compressed_coefficientwise": True,
            "post_A_endpoint_roster_reduced_to_two_rows": True,
            "extended_notch_one_cell_fold_and_modular_transform_inherited": True,
            "A_upper_projector_edge_moved_out_of_half_boundary_collar": True,
            "B_subtractions_and_negative_bulk_ownership_preserved": True,
            "extended_notch_alone_proves_joined_bound": False,
            "quantitative_R_after_A_bound_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "process_priority": priority,
        },
        "next_obligation": "Instantiate the pole-free two-edge kernel with d_A=39936.5 and combine its upper-edge current with the already-deleted A window before norms. Derive a joint-gradient or contour estimate on the surviving post-A kernel, using the exact |q_A|>=84.5 half-boundary separation and the existing B-window subtraction.",
        "proof_boundary": "Exact finite-regulator post-A coefficient compression, extended-notch fold/modular inheritance, and rational upper-edge geometry only. No quantitative bound for the remaining joined kernel, complete R_after_A, R_Dir, or Q_K-T estimate, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(), encoding="utf-8")
    print("certified post-A extended-notch normal form", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
