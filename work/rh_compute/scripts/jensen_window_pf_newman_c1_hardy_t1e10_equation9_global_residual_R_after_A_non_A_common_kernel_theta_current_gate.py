#!/usr/bin/env python3
"""Derive executable common-kernel forms for the non-A post-A remainder."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_common_kernel_theta_current_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "finite_regulator_equivalence": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_finite_regulator_equivalence_gate.json",
    "extended_notch_normal_form": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_normal_form_gate.json",
    "extended_notch_pole_free_kernel": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_pole_free_kernel_gate.json",
    "one_cell_fold": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_punctured_cell_fold_gate.json",
    "two_jet_reassembly": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_periodic_two_jet_common_kernel_reassembly_gate.json",
    "A_face_embedding": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_fixed_regulator_embedding_gate.json",
}

T_HEIGHT = 10_000_000_000
A = 159_577
B = 5_122_421
L = (B - A) // 2
FIRST_TARGET = 622
LAST_EXTENDED_TARGET = 39_936
R_AFTER_A_TARGET = Decimal("0.0368147039947")
A_FACE_THRESHOLD = Decimal("0.00364")
NON_A_TARGET = R_AFTER_A_TARGET - A_FACE_THRESHOLD


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


def symbolic_theta_certificate() -> dict[str, Any]:
    n, a, s, x = sp.symbols("n a s x", real=True)
    alpha = a + 2 * n + 2 * s
    base = (a + 2 * s) ** 2 / 4
    reduced = n**2 + n * (a + 2 * s)
    require(sp.expand(alpha**2 / 4 - base - reduced) == 0, "quadratic theta phase split drift")

    alpha_next = alpha + 2
    first_ratio_phase = sp.expand((alpha_next**2 - alpha**2) / 4)
    require(sp.simplify(first_ratio_phase - (alpha + 1)) == 0, "first phase ratio drift")
    alpha_next2 = alpha + 4
    second_ratio_phase = sp.expand((alpha_next2**2 - 2 * alpha_next**2 + alpha**2) / 4)
    require(sp.simplify(second_ratio_phase - 2) == 0, "second phase ratio drift")

    # F_0 is checked without summing a symbolic exponential.
    l = sp.symbols("L", positive=True, integer=True)
    zero_sum = l * (a + 2 * s) + 2 * l * (l - 1) / 2
    require(sp.simplify(zero_sum - l * (a + 2 * s + l - 1)) == 0, "x=0 theta-current extension drift")

    return {
        "phase_factorization": "(A+2n+2s)^2/4=(A+2s)^2/4+n^2+n(A+2s)",
        "quadratic_sum": "S_0(x,s)=sum_(n=0)^(L-1) exp(i*pi*x[n^2+n(A+2s)])",
        "first_moment": "S_1(x,s)=sum_(n=0)^(L-1) n exp(i*pi*x[n^2+n(A+2s)])",
        "theta_current": "F_x(s)=exp(i*pi*x(A+2s)^2/4)[(A+2s)S_0(x,s)+2S_1(x,s)]",
        "term_ratio": "z_(n+1)/z_n=exp(i*pi*x[A+2n+2s+1])",
        "ratio_of_ratios": "(z_(n+2)/z_(n+1))/(z_(n+1)/z_n)=exp(2*pi*i*x)",
        "zero_x_extension": "F_0(s)=L(A+2s+L-1)",
        "roster": {"A": A, "B": B, "L": L, "n_range": [0, L - 1]},
    }


def coefficient_roster(finite: dict[str, Any]) -> list[dict[str, Any]]:
    rows = finite["finite_equivalence_certificate"]["positive_class_rows"]
    expected = [
        ([1, 621], "P_m+P_-m+A_m+A_-m+B_m+B_-m"),
        ([622, 39_852], "P_-m+A_m+A_-m+B_m+B_-m"),
        ([39_853, 39_894], "P_-m+B_m+B_-m"),
        ([39_895, 39_936], "P_-m+B_m+B_-m"),
        ([39_937, "M"], "P_m+P_-m+A_m+A_-m+B_m+B_-m"),
    ]
    require(len(rows) == len(expected), "positive class row count drift")
    for row, (mode_range, summand) in zip(rows, expected):
        require(row["mode_range"] == mode_range, f"mode range drift: {mode_range}")
        require(row["joined_summand"] == summand, f"joined summand drift: {mode_range}")
    return rows


def render_note(artifact: dict[str, Any]) -> str:
    theta = artifact["theta_current_certificate"]
    return f"""# Non-A common kernel and quadratic-theta current

Date: 2026-08-24

Status: exact executable saved-height reduction; quantitative non-A bound open

At one common cutoff `M>=B=5122421` and regulator `epsilon>0`, Section
11.459 defines

```text
mathfrak K_nonA,(M,epsilon)
 =mathfrak R_A,(M,epsilon)-mathfrak E_(A,42),epsilon. (NK1)
```

The coefficient-mask evaluator is

```text
H_x+I_0
 +sum_(m=1)^621       w_m(P_m+P_-m+A_m+A_-m+B_m+B_-m)
 +sum_(m=622)^39852   w_m(P_-m+A_m+A_-m+B_m+B_-m)
 +sum_(m=39853)^39936 w_m(P_-m+B_m+B_-m)
 +sum_(m=39937)^M     w_m(P_m+P_-m+A_m+A_-m+B_m+B_-m)
 -chi_W Btr_(M,epsilon)-O_(M,epsilon)
 -mathfrak E_(A,42),epsilon.                         (NK2)
```

The two adjacent 42-mode rows have been combined in (NK2) because their
surviving coefficient vectors are identical.  No coefficient changed when
the analytic A-face channel was extracted.

Let `U={{622,...,39936}}`.  The exactly equivalent extended-notch form is

```text
H_x+integral_0^L f_x(u)[D_(M,epsilon)-G_(U,epsilon)]du
 +sum_(m=622)^39852 w_m(A_m+B_m)
 +sum_(m=39853)^39936 w_m(B_m-A_-m)
 -chi_W Btr_(M,epsilon)-O_(M,epsilon)
 -mathfrak E_(A,42),epsilon.                         (NK3)
```

Equivalently, after the exact two-jet reassembly,

```text
H_x+E_x+2*pi*i sum_(m=1)^M m w_m
       [hat Phi_x(m)-hat Phi_x(-m)]
 -sum_(m=622)^39936 w_m P_m
 -sum_(m=39853)^39936 w_m(A_m+A_-m)
 -chi_W Btr_(M,epsilon)-O_(M,epsilon)
 -mathfrak E_(A,42),epsilon.                         (NK4)
```

The pair series in (NK4) is absolutely and uniformly convergent after the
two endpoint jets are removed.  Its inherited raw derivative majorant is
still far too large for the target, so this identity is a convergence-safe
cross-check rather than the selected numerical enclosure.

The source term in (NK3) folds exactly to one cell:

```text
integral_0^L f_x(u)C_(U,M,epsilon)(u)du
 =integral_0^1 F_x(s)C_(U,M,epsilon)(s)ds,           (NK5)

F_x(s)=sum_(n=0)^(L-1)(A+2n+2s)
       exp(i*pi*x(A+2n+2s)^2/4),
L={L}.                                               (NK6)
```

After `M` tends to infinity at fixed `epsilon`, the certified pole-free
centred-cell formula evaluates `C_U,epsilon` with a cubic dual tail.  The
remaining large object in (NK5) is not arbitrary: put

```text
S_0=sum_(n=0)^(L-1) exp(i*pi*x[n^2+n(A+2s)]),
S_1=sum_(n=0)^(L-1) n exp(i*pi*x[n^2+n(A+2s)]).      (NK7)
```

Exact quadratic completion gives

```text
F_x(s)=exp(i*pi*x(A+2s)^2/4)
       [(A+2s)S_0+2S_1].                            (NK8)
```

For direct reference evaluation, if `z_n` is the summand in `S_0`, then

```text
z_(n+1)/z_n=exp(i*pi*x[A+2n+2s+1]),

(z_(n+2)/z_(n+1))/(z_(n+1)/z_n)=exp(2*pi*i*x).      (NK9)
```

Thus one initial phase and two multiplicative recurrences evaluate both
`S_0` and `S_1` without repeated large-argument trigonometric calls.  At
`x=0`, the removable value is exactly

```text
F_0(s)=L(A+2s+L-1).                                  (NK10)
```

Equations (NK7)--(NK10) define the reference oracle for a fast incomplete
quadratic-Gauss recursion.  The next numerical stage must cross-check that
fast evaluator against this recurrence on short and full rosters, then
compose it with the pole-free `C_U,epsilon` kernel on a bounded lattice.
Neither a finite-Gauss error theorem nor the physical `x,s` quadrature has
yet been proved.

The inherited physical limit is

```text
R_nonA=lim_(epsilon down 0)lim_(M to infinity)
       mathcal P_t[mathfrak K_nonA,(M,epsilon)],      (NK11)
```

and the sufficient target remains

```text
|R_nonA|<0.0331747039947.                            (NK12)
```

No row of (NK2) may be normed independently: its zero, half-current,
negative, B-completion, and remote-positive cancellations are precisely why
(NK3)--(NK5) are retained as common-kernel evaluators.

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Pi provenance: every `pi` in (NK3)--(NK10) is inherited from the original
Kummer quadratic phase, integer Fourier character, Gaussian Abel regulator,
and Poisson normalization.  The algebraic phase factorization introduces no
new geometric or fitted occurrence.

Proof boundary: (NK2)--(NK10) are exact saved-height common-kernel and finite
quadratic-theta reductions.  No fast-evaluator error theorem, physical
quadrature, non-A bound, joined `R_after_A`, `R_Dir`, `Q_K-T`, all-height
theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a common-kernel dependency is not passed")

    finite = dependencies["finite_regulator_equivalence"]
    normal = dependencies["extended_notch_normal_form"]
    pole = dependencies["extended_notch_pole_free_kernel"]
    fold = dependencies["one_cell_fold"]
    pair = dependencies["two_jet_reassembly"]
    embedding = dependencies["A_face_embedding"]

    require(finite["scope"]["height"] == T_HEIGHT, "height drift")
    require(normal["scope"]["extended_target"] == [FIRST_TARGET, LAST_EXTENDED_TARGET], "extended target drift")
    require(normal["decision"]["post_A_endpoint_roster_reduced_to_two_rows"] is True, "two-row endpoint roster missing")
    require(pole["decision"]["pole_free_closed_centred_cell_formula_certified"] is True, "pole-free kernel missing")
    require(fold["decision"]["long_u_integral_folded_exactly_to_one_periodic_cell"] is True, "one-cell fold missing")
    require(fold["decision"]["folded_source_is_exact_finite_quadratic_theta_current"] is True, "finite theta current missing")
    require(pair["decision"]["zero_regulator_pair_series_absolute_and_uniform"] is True, "two-jet convergence missing")
    require(pair["decision"]["raw_derivative_majorant_closes_target"] is False, "raw-majorant obstruction drift")
    require(embedding["decision"]["fixed_regulator_A_face_channel_embedded_exactly"] is True, "A-face embedding missing")
    require(embedding["decision"]["six_class_mask_unchanged"] is True, "mask guard missing")
    require(CHECKER.is_file(), "independent checker missing")

    roster = coefficient_roster(finite)
    theta = symbolic_theta_certificate()
    require(L == 2_481_422, "source roster length drift")
    require(NON_A_TARGET == Decimal("0.0331747039947"), "non-A target drift")

    artifact = {
        "kind": STEM,
        "status": "exact_non_A_common_kernel_and_quadratic_theta_current_reduction_certified_fast_evaluator_and_bound_open",
        "passed": True,
        "scope": {
            "height": T_HEIGHT,
            "source_odd_roster": [A, B],
            "source_cell_count": L,
            "extended_target": [FIRST_TARGET, LAST_EXTENDED_TARGET],
            "common_cutoff": f"M>=B={B}",
            "common_regulator": "epsilon>0, M to infinity first, epsilon down to zero second",
        },
        "coefficient_mask_form": {
            "zero_sector": "H_x+I_0",
            "positive_rows": roster,
            "subtractions": ["chi_W Btr_(M,epsilon)", "O_(M,epsilon)", "mathfrak E_(A,42),epsilon"],
            "mask_changed_by_A_face_extraction": False,
        },
        "extended_notch_form": {
            "source_kernel": "integral_0^L f_x(u)[D_(M,epsilon)-G_(U,epsilon)]du, U={622,...,39936}",
            "endpoint_rows": {
                "622..39852": "A_m+B_m",
                "39853..39936": "B_m-A_-m",
            },
            "common_subtractions": "chi_W Btr+O+mathfrak E_(A,42)",
        },
        "two_jet_form": {
            "kernel": "H_x+E_x+2*pi*i*sum_(m=1)^M m*w_m[hat Phi_x(m)-hat Phi_x(-m)]-sum_(m=622)^39936 w_m P_m-sum_(m=39853)^39936 w_m(A_m+A_-m)-chi_W Btr-O-mathfrak E_(A,42)",
            "absolute_uniform_pair_convergence": True,
            "raw_derivative_majorant_quantitatively_useful": False,
        },
        "one_cell_form": {
            "identity": "integral_0^L f_x(u)C_(U,M,epsilon)(u)du=integral_0^1 F_x(s)C_(U,M,epsilon)(s)ds",
            "target_kernel": "C_(U,M,epsilon)=D_(M,epsilon)-G_(U,epsilon)",
            "fixed_epsilon_limit": "C_U,epsilon has the certified pole-free centred-cell formula and cubic dual tail",
        },
        "theta_current_certificate": theta,
        "route_decision": {
            "direct_five_row_norm_route_admissible": False,
            "direct_mask_is_exact_reference_evaluator": True,
            "two_jet_form_is_convergence_safe_cross_check": True,
            "raw_two_jet_majorant_selected": False,
            "one_cell_pole_free_theta_current_route_selected_for_pilot": True,
            "reference_quadratic_recurrence_exact": True,
            "fast_incomplete_quadratic_Gauss_evaluator_built": False,
            "fast_evaluator_error_theorem_proved": False,
            "physical_x_s_quadrature_proved": False,
        },
        "budget_certificate": {
            "R_after_A_target": str(R_AFTER_A_TARGET),
            "A_face_absolute_threshold": str(A_FACE_THRESHOLD),
            "non_A_sufficient_absolute_target": str(NON_A_TARGET),
        },
        "decision": {
            "non_A_common_kernel_defined_on_one_regulator": True,
            "coefficient_mask_and_extended_notch_forms_equivalent": True,
            "extended_notch_and_two_jet_forms_equivalent": True,
            "long_source_integral_folded_to_one_cell": True,
            "folded_source_reduced_to_quadratic_sum_and_first_moment": True,
            "reference_term_recurrence_proved": True,
            "zero_x_theta_current_extension_proved": True,
            "separate_positive_row_norms_forbidden": True,
            "non_A_bound_proved": False,
            "R_after_A_bound_proved": False,
            "R_Dir_bound_proved": False,
            "QK_minus_T_bound_proved": False,
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
        "runtime": {"elapsed_seconds": round(time.time() - started, 3), "workers": 1},
        "next_obligation": "Implement a bounded reference evaluator for S_0 and S_1 using the exact two-ratio recurrence, then cross-check a fast incomplete quadratic-Gauss recursion on deterministic short rosters and selected full-roster points. Compose only a validated evaluator with the pole-free C_U,epsilon kernel; retain all B, zero, half-current, negative, and remote-positive terms in the common physical quadrature.",
        "proof_boundary": "Exact saved-height non-A coefficient-mask, extended-notch, two-jet, one-cell, and quadratic-theta-current forms only. No fast-evaluator error theorem, physical quadrature, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("certified exact non-A common kernel and quadratic-theta current reduction", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
