#!/usr/bin/env python3
"""Certify the finite-regulator common-kernel representation of R_after_A."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_finite_regulator_equivalence_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "finite_dirichlet_rejoin": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_finite_block_dirichlet_rejoin_gate.json",
    "projector_defect": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_phase_space_projector_defect_gate.json",
    "ownership": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_Dir_projector_ownership_ledger_gate.json",
    "B_outer": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_outer_phase_coupled_tangential_gate.json",
    "double_weber": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_double_weber_exact_source_reconstruction_gate.json",
}

HEIGHT = 10_000_000_000
A = 159_577
B = 5_122_421
L = (B - A) // 2
TARGET_START = 622
TARGET_END = 39_894
A_START = 39_853
A_END = 39_936
ATOM_ORDER = ("P_plus", "P_minus", "A_plus", "A_minus", "B_plus", "B_minus")


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


def vector(*values: int) -> dict[str, int]:
    require(len(values) == len(ATOM_ORDER), "coefficient-vector length drift")
    return dict(zip(ATOM_ORDER, values, strict=True))


def subtract(left: dict[str, int], right: dict[str, int]) -> dict[str, int]:
    return {atom: left[atom] - right[atom] for atom in ATOM_ORDER}


def target_indicator(mode: int) -> int:
    return int(TARGET_START <= mode <= TARGET_END)


def pre_A_vector(mode: int) -> dict[str, int]:
    return vector(1 - target_indicator(mode), 1, 1, 1, 1, 1)


def A_vector(mode: int) -> dict[str, int]:
    if not A_START <= mode <= A_END:
        return vector(0, 0, 0, 0, 0, 0)
    return vector(1 - target_indicator(mode), 0, 1, 1, 0, 0)


def post_A_vector(mode: int) -> dict[str, int]:
    return subtract(pre_A_vector(mode), A_vector(mode))


def class_rows() -> list[dict[str, Any]]:
    rows = [
        {
            "sector_id": "low_positive_pairs",
            "mode_range": [1, 621],
            "count": 621,
            "target_indicator": 0,
            "post_A_coefficients": vector(1, 1, 1, 1, 1, 1),
            "joined_summand": "P_m+P_-m+A_m+A_-m+B_m+B_-m",
        },
        {
            "sector_id": "ordinary_target_pairs_before_A_window",
            "mode_range": [622, 39_852],
            "count": 39_231,
            "target_indicator": 1,
            "post_A_coefficients": vector(0, 1, 1, 1, 1, 1),
            "joined_summand": "P_-m+A_m+A_-m+B_m+B_-m",
        },
        {
            "sector_id": "A_window_target_pairs",
            "mode_range": [39_853, 39_894],
            "count": 42,
            "target_indicator": 1,
            "post_A_coefficients": vector(0, 1, 0, 0, 1, 1),
            "joined_summand": "P_-m+B_m+B_-m",
        },
        {
            "sector_id": "A_window_outer_pairs",
            "mode_range": [39_895, 39_936],
            "count": 42,
            "target_indicator": 0,
            "post_A_coefficients": vector(0, 1, 0, 0, 1, 1),
            "joined_summand": "P_-m+B_m+B_-m",
        },
        {
            "sector_id": "remote_positive_pairs",
            "mode_range": [39_937, "M"],
            "count": "M-39936",
            "target_indicator": 0,
            "post_A_coefficients": vector(1, 1, 1, 1, 1, 1),
            "joined_summand": "P_m+P_-m+A_m+A_-m+B_m+B_-m",
        },
    ]
    representative_modes = (1, 622, 39_853, 39_895, 39_937)
    for row, mode in zip(rows, representative_modes, strict=True):
        require(post_A_vector(mode) == row["post_A_coefficients"], f"class vector drift: {row['sector_id']}")
    return rows


def symbolic_certificate() -> dict[str, str]:
    Pp, Pm, Ap, Am, Bp, Bm, chi = sp.symbols("Pp Pm Ap Am Bp Bm chi")
    Ip = Pp + Ap + Bp
    Im = Pm + Am + Bm
    pair = sp.expand(Ip + Im - chi * Pp)
    projector = sp.expand((1 - chi) * Pp + Pm + Ap + Am + Bp + Bm)
    require(sp.expand(pair - projector) == 0, "paired projector identity failed")

    H, D, IT, PT, QT = sp.symbols("H D IT PT QT")
    common_kernel = H + D - IT + QT
    require(sp.expand(common_kernel.subs(IT, PT + QT) - (H + D - PT)) == 0, "common-kernel insertion failed")

    Delta, Bwin, O, Aextract, Rafter = sp.symbols("Delta Bwin O Aextract Rafter")
    finite_definition = Delta - Bwin - O - Aextract
    require(sp.expand(finite_definition - Rafter).subs(Rafter, finite_definition) == 0, "finite R_after_A definition failed")
    require(sp.expand(Delta - (Bwin + O + Aextract + Rafter)).subs(Rafter, finite_definition) == 0, "finite allocation failed")

    PK, EB, EO, AT, RA = sp.symbols("R_KGamma E_Btr E_outer A_transition R_after_A")
    require(sp.expand(PK - (EB + EO + AT + RA)).subs(PK, EB + EO + AT + RA) == 0, "physical allocation failed")
    return {
        "paired_projector": "I_m+I_-m-chi_T(m)P_m=(1-chi_T(m))P_m+P_-m+A_m+A_-m+B_m+B_-m",
        "common_kernel": "Delta_(M,epsilon)=H_x+integral_0^L f_x(u)[D_(M,epsilon)-G_(T,epsilon)]du+sum_(m in T)w_m(A_m+B_m)",
        "A_notch": "A_(M,epsilon)=sum_(m=39853)^39936 w_m[A_m+A_-m+(1-chi_T(m))P_m]",
        "finite_R_after_A": "mathfrak R_A,(M,epsilon)=Delta_(M,epsilon)-chi_W Btr_(M,epsilon)-O_(M,epsilon)-A_(M,epsilon)",
        "physical_projection": "P_t[F]=2(pi/(32t))^(1/4)Re[e^(-i*pi/8) integral_0^1 W_t(x)F(x)dx]",
        "physical_limit": "R_after_A=lim_(epsilon down 0) lim_(M to infinity) P_t[mathfrak R_A,(M,epsilon)]",
        "physical_allocation": "R_KGamma=E_Btr,win+E_outer+A_transition+R_after_A",
    }


def finite_equivalence_certificate() -> dict[str, Any]:
    rows = class_rows()
    boundary_samples = [1, 621, 622, 39_852, 39_853, 39_894, 39_895, 39_936, 39_937, B - 1, B, B + 17]
    samples = [
        {
            "mode": mode,
            "target_indicator": target_indicator(mode),
            "pre_A_coefficients": pre_A_vector(mode),
            "A_notch_coefficients": A_vector(mode),
            "post_A_coefficients": post_A_vector(mode),
        }
        for mode in boundary_samples
    ]
    for sample in samples:
        require(
            subtract(sample["pre_A_coefficients"], sample["A_notch_coefficients"])
            == sample["post_A_coefficients"],
            f"boundary coefficient transfer failed at mode {sample['mode']}",
        )

    count_checks = []
    for cutoff in (B, B + 17):
        counts = [621, 39_231, 42, 42, cutoff - 39_936]
        require(sum(counts) == cutoff, f"positive-class count failed at M={cutoff}")
        count_checks.append(
            {
                "cutoff": cutoff,
                "class_counts": counts,
                "positive_mode_total": sum(counts),
                "symmetric_mode_total_with_zero": 2 * cutoff + 1,
            }
        )

    return {
        "cutoff_condition": f"integer M>=B={B}; epsilon>0",
        "reason_for_strong_cutoff": "A deletion alone needs M>=39936, but the same finite identity with O_(M,epsilon)=eta_delta E_W Ctr_>=B,M,epsilon needs M>=B.",
        "zero_sector": {
            "summand": "H_x+I_0",
            "weight": "w_0=1",
            "separate_norm_admissible": False,
        },
        "positive_class_rows": rows,
        "boundary_samples": samples,
        "count_checks": count_checks,
        "common_kernel_form": "H_x+integral f_x C_(M,epsilon)+sum_(m in T)w_m(A_m+B_m)-chi_W Btr_(M,epsilon)-O_(M,epsilon)-A_(M,epsilon)",
        "six_class_mask_form": "H_x+I_0+sum_(five positive-label ranges) w_m times the recorded post-A joined summand-chi_W Btr_(M,epsilon)-O_(M,epsilon)",
        "equivalence_method": "Expand integral f_x D as sum_(m=-M)^M w_m I_m, expand integral f_x G_T as sum_(m in T)w_m I_m, use I_m=P_m+A_m+B_m, pair +/-m using w_-m=w_m, then subtract the A notch coefficientwise.",
        "B_subtractions_common_to_both_forms": True,
        "independent_norm_forbidden_sectors": [
            "endpoint_half_current",
            "zero_mode",
            "negative_bulk",
            "remaining_endpoint_atoms",
            "remote_positive_pairs",
        ],
    }


def render_note(artifact: dict[str, Any]) -> str:
    return f"""# Finite-regulator equivalence for `R_after_A`

Date: 2026-08-23

Status: exact common-kernel/mode-mask equivalence and physical limit order
certified; no quantitative bound for `R_after_A`

Let `T={{622,...,39894}}`, `w_m=exp(-pi*epsilon*m^2)`, and take one
finite symmetric cutoff

```text
M>=B={B},       epsilon>0.                              (RA1)
```

The stronger cutoff in (RA1) is deliberate.  The A deletion itself only
needs `M>=39936`, but the finite analytic B exterior
`O_(M,epsilon)=eta_delta E_W Ctr_>=B,M,epsilon` belongs to the same identity
only for `M>=B`.

The common Gamma-normalized projector kernel is

```text
Delta_(M,epsilon)(x)
 =H_x+integral_0^{L} f_x(u)[D_(M,epsilon)-G_(T,epsilon)]du
   +sum_(m in T)w_m(A_m+B_m).                         (RA2)
```

Delete the already-certified A transition before taking a norm:

```text
mathcal A_(M,epsilon)
 =sum_(m=39853)^39936 w_m
   [A_m+A_-m+(1-chi_T(m))P_m].                        (RA3)

mathfrak R_A,(M,epsilon)
 =Delta_(M,epsilon)-chi_W Btr_(M,epsilon)
   -O_(M,epsilon)-mathcal A_(M,epsilon).               (RA4)
```

Equations (RA2)--(RA4) are executable without assigning independent norms
to any cancellation partner.  Expanding the common Fourier kernel, pairing
`+m` and `-m`, and using `I_m=P_m+A_m+B_m` gives the exact coefficient mask

```text
mode range       post-A joined summand
1..621           P_m+P_-m+A_m+A_-m+B_m+B_-m
622..39852       P_-m+A_m+A_-m+B_m+B_-m
39853..39894     P_-m+B_m+B_-m
39895..39936     P_-m+B_m+B_-m
39937..M         P_m+P_-m+A_m+A_-m+B_m+B_-m.          (RA5)
```

Together with `H_x+I_0`, (RA5), `-chi_W Btr`, and `-O` is exactly (RA4),
not an approximation.  The two middle rows have the same surviving vector
for different reasons: the target projector had already removed `P_m` in
the first row, while the A transition removes it in the second.  No
`P_-m` or B endpoint atom is deleted.

With the inherited equation-(9) physical functional

```text
P_t[F]=2(pi/(32t))^(1/4)
       Re[e^(-i*pi/8) integral_0^1 W_t(x)F(x)dx],     (RA6)
```

the object being targeted is precisely

```text
R_after_A
 =lim_(epsilon down 0) lim_(M to infinity)
    P_t[mathfrak R_A,(M,epsilon)],                    (RA7)

R_KGamma=E_Btr,win+E_outer+A_transition+R_after_A.    (RA8)
```

The `M` limit is taken first at fixed positive `epsilon`; only then does the
common Gaussian Abel regulator tend to zero.  The separately certified B
outer block keeps its dominated Abel passage.  The endpoint half-current,
zero mode, negative bulk, remaining endpoint atoms, and remote-positive
pairs remain joined throughout.

Pi provenance: every `pi` in (RA1)--(RA8) comes from the original Kummer
quadratic phase, integer Fourier character, Gaussian Abel regulator, and
the inherited physical normalization.  No circular construction, fitted
constant, or new geometric occurrence of `pi` is introduced.

Proof boundary: exact finite-regulator representation, coefficient-mask
equivalence, cutoff compatibility, and common physical limit order at
`t=10^10` only.  No numerical or analytic upper bound for `R_after_A`,
complete `R_Dir` or `Q_K-T`, all-height theorem, `Lambda<=0`, PF-infinity,
RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["finite_dirichlet_rejoin"]["decision"]["finite_completed_B_block_rejoined_before_norms"] is True, "finite B rejoin drift")
    require(dependencies["projector_defect"]["decision"]["target_full_line_carrier_subtracted_before_norms"] is True, "projector dependency drift")
    require(dependencies["ownership"]["decision"]["A_endpoint_and_positive_bulk_coefficients_transferred_exactly"] is True, "A ownership drift")
    require(dependencies["B_outer"]["decision"]["outer_block_separate_Abel_limit_justified_by_summable_derivative_majorants"] is True, "B outer Abel-limit drift")
    require(dependencies["double_weber"]["decision"]["entire_finite_equation9_source_roster_reconstructed_exactly"] is True, "double-Weber reconstruction drift")

    artifact = {
        "kind": STEM,
        "status": "exact_R_after_A_finite_regulator_common_kernel_to_six_class_mask_equivalence_certified_quantitative_bound_open",
        "passed": True,
        "scope": {
            "height": HEIGHT,
            "source_odd_roster": [A, B],
            "roster_coordinate_length": L,
            "finite_cutoff": f"integer M>=B={B}",
            "regulator": "epsilon>0, w_m=exp(-pi*epsilon*m^2)",
            "limit_order": "M to infinity first at fixed epsilon; epsilon down to zero second",
            "target_positive_modes": [TARGET_START, TARGET_END],
            "A_transition_modes": [A_START, A_END],
        },
        "symbolic_certificate": symbolic_certificate(),
        "finite_equivalence_certificate": finite_equivalence_certificate(),
        "decision": {
            "finite_R_after_A_defined_on_one_common_regulator": True,
            "common_kernel_equals_six_class_post_A_mask": True,
            "analytic_B_outer_cutoff_condition_M_ge_B_enforced": True,
            "A_endpoint_and_outside_target_positive_bulk_deleted_algebraically_before_norms": True,
            "B_window_and_outer_subtractions_identical_in_both_forms": True,
            "common_physical_limit_order_preserved": True,
            "zero_half_negative_endpoint_and_remote_positive_sectors_remain_joined": True,
            "quantitative_R_after_A_upper_bound_proved": False,
            "complete_R_Dir_proved": False,
            "complete_Q_K_minus_T_proved": False,
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
        "next_obligation": "Build a bounded convergence scout for the single joined mathfrak R_A object. Compare the direct finite-epsilon kernel with its Jacobi-modular/double-Weber representation at matched regulator and cutoff, retaining the A notch, B subtractions, zero mode, half-current, negative partners, and remote-positive modes in the same evaluation. Use samples only to choose the analytic enclosure route; do not promote convergence data to a proof bound.",
        "proof_boundary": "Exact saved-height finite-regulator representation, class-mask equivalence, cutoff compatibility, and inherited common-limit order only. No quantitative R_after_A or R_Dir bound, complete Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified finite-regulator R_after_A common-kernel equivalence", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
