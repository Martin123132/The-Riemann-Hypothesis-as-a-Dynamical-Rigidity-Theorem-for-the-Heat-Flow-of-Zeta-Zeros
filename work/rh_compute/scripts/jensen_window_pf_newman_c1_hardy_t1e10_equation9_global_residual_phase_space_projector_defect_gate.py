#!/usr/bin/env python3
"""Reduce the Gamma-normalized residual to one signed projector defect."""

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

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_phase_space_projector_defect_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "uniform_Gamma_insertion": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_uniform_gamma_target_insertion_gate.json",
    "finite_rejoin": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_finite_block_dirichlet_rejoin_gate.json",
    "symmetric_Poisson": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate.json",
    "pair_triangle": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_separable_triangle_geometry_gate.json",
    "bi_Morse": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_exact_bimorse_face_fold_gate.json",
    "ownership": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_ordinary_mode_coverage_ledger_gate.json",
}

HEIGHT = 10_000_000_000
A = 159_577
B = 5_122_421
L = 2_481_422
TARGET_START = 622
TARGET_END = 39_894
TARGET_COUNT = TARGET_END - TARGET_START + 1
OWNERSHIP_SPLIT = (39_694, 39_695)
WORKING_TARGET = arb("0.00014057919999999995")


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


def symbolic_certificate() -> dict[str, str]:
    half, outside_interval, target_interval, target_full, target_tail = sp.symbols(
        "H O_U I_T P_T Q_T"
    )
    all_interval = outside_interval + target_interval
    require(
        sp.expand(
            half + all_interval - target_full
            - (half + outside_interval + target_tail)
        ).subs(target_tail, target_interval - target_full)
        == 0,
        "finite projector split failed",
    )
    exterior_target = target_full - target_interval
    require(
        sp.expand(
            half + outside_interval - exterior_target
            - (half + all_interval - target_full)
        )
        == 0,
        "oriented rectangle identity failed",
    )

    mode = sp.symbols("m", integer=True, positive=True)
    t, endpoint = sp.symbols("t D", positive=True)
    alpha = 2 * mode + t / (sp.pi * mode)
    x_mode = 2 * sp.pi * mode**2 / (t + 2 * sp.pi * mode**2)
    z_endpoint = 2 * mode / endpoint
    expected_gap = 2 * sp.pi * mode**2 * (endpoint - alpha) / (
        endpoint * (t + 2 * sp.pi * mode**2)
    )
    require(sp.simplify(x_mode - z_endpoint - expected_gap) == 0, "face-gap identity failed")
    require(
        sp.simplify(
            x_mode
            - sp.Rational(1, 2)
            - (2 * sp.pi * mode**2 - t) / (2 * (t + 2 * sp.pi * mode**2))
        )
        == 0,
        "half-boundary identity failed",
    )

    return {
        "finite_source": "S_(M,epsilon)=H_x+sum_(m=-M)^M w_m I_m(x)",
        "target_full_line": "P_(T,epsilon)=sum_(m=622)^39894 w_m P_m(x)",
        "finite_defect": "Delta_(M,epsilon)=H_x+sum_(m=-M)^M w_m I_m-sum_(m=622)^39894 w_m P_m",
        "outside_interval": "O_(M,epsilon)=sum_(m in [-M,621] union [39895,M]) w_m I_m",
        "target_exterior": "X_(T,epsilon)=sum_(m=622)^39894 w_m[P_m-I_m]",
        "oriented_projector_defect": "Delta_(M,epsilon)=H_x+O_(M,epsilon)-X_(T,epsilon)",
        "integral_form": "I_m=int_[0,L] K_x(u,m)du, P_m=Abel-int_R K_x(u,m)du, K_x=f_x(u)exp(-2*pi*i*m*u)",
        "face_gap": "x_m-z_D=2*pi*m^2[D-alpha_m]/[D(t+2*pi*m^2)], alpha_m=2m+t/(pi*m)",
        "half_boundary_gap": "x_m-1/2=[2*pi*m^2-t]/[2(t+2*pi*m^2)]",
    }


def mode_partition_certificate() -> dict[str, Any]:
    rows: list[dict[str, int]] = []
    for cutoff in (TARGET_END, B, B + 17):
        low_count = cutoff + TARGET_START
        high_count = max(0, cutoff - TARGET_END)
        outside_count = 2 * cutoff + 1 - TARGET_COUNT
        require(low_count + high_count == outside_count, f"outside mode count failed at M={cutoff}")
        rows.append(
            {
                "cutoff": cutoff,
                "low_outside_count": low_count,
                "high_outside_count": high_count,
                "outside_count": outside_count,
                "target_count": TARGET_COUNT,
                "full_symmetric_count": 2 * cutoff + 1,
            }
        )
    return {
        "target_modes": [TARGET_START, TARGET_END],
        "target_count": TARGET_COUNT,
        "outside_modes_for_M_ge_39894": "[-M,621] union [39895,M] (empty upper interval when M=39894)",
        "count_identity": "|F_M\\T|=(M+622)+(M-39894)=2M+1-39273 for M>=39894",
        "rows": rows,
    }


def signed(value: arb) -> int:
    if value.lower() > 0:
        return 1
    if value.upper() < 0:
        return -1
    return 0


def transition_certificate() -> dict[str, Any]:
    ctx.dps = 90
    pi = arb.pi()
    t = arb(HEIGHT)

    def alpha(mode: int) -> arb:
        m = arb(mode)
        return 2 * m + t / (pi * m)

    def face_gap(endpoint: int, mode: int) -> arb:
        m = arb(mode)
        return 2 * pi * m**2 * (arb(endpoint) - alpha(mode)) / (
            arb(endpoint) * (t + 2 * pi * m**2)
        )

    def half_gap(mode: int) -> arb:
        m = arb(mode)
        return (2 * pi * m**2 - t) / (2 * (t + 2 * pi * m**2))

    values = {
        "B_face_m621": face_gap(B, 621),
        "B_face_m622": face_gap(B, 622),
        "A_face_m39852": face_gap(A, 39_852),
        "A_face_m39853": face_gap(A, 39_853),
        "half_boundary_m39894": half_gap(39_894),
        "half_boundary_m39895": half_gap(39_895),
        "ownership_A_face_m39694": face_gap(A, OWNERSHIP_SPLIT[0]),
        "ownership_A_face_m39695": face_gap(A, OWNERSHIP_SPLIT[1]),
        "ownership_half_m39694": half_gap(OWNERSHIP_SPLIT[0]),
        "ownership_half_m39695": half_gap(OWNERSHIP_SPLIT[1]),
    }
    require(signed(values["B_face_m621"]) == -1 and signed(values["B_face_m622"]) == 1, "B transition drift")
    require(signed(values["A_face_m39852"]) == -1 and signed(values["A_face_m39853"]) == 1, "A transition drift")
    require(signed(values["half_boundary_m39894"]) == -1 and signed(values["half_boundary_m39895"]) == 1, "half-boundary transition drift")
    require(signed(values["ownership_A_face_m39694"]) == signed(values["ownership_A_face_m39695"]) == -1, "ownership split crosses A face")
    require(signed(values["ownership_half_m39694"]) == signed(values["ownership_half_m39695"]) == -1, "ownership split crosses half boundary")
    require(4 * 39_894 == A - 1 and 4 * 39_895 == A + 3, "A null-face lattice bracket drift")
    return {
        "height": HEIGHT,
        "values": {key: value.str(70, more=True) for key, value in values.items()},
        "signs": {key: signed(value) for key, value in values.items()},
        "stationary_occupancy_interfaces": {
            "B_face": "621|622",
            "A_face": "39852|39853",
            "half_boundary_corner": "39894|39895",
        },
        "proof_ownership_interface": "39694|39695; no stationary sign changes there",
        "A_null_face_lattice_bracket": "4*39894=A-1 and 4*39895=A+3",
    }


def raw_norm_obstruction(symmetric_poisson: dict[str, Any]) -> dict[str, str]:
    ctx.dps = 90
    coefficient = arb(
        symmetric_poisson["certified_constants"]["tail_coefficient_C_AB_over_2pi2_ball"]
    )
    at_target_edge = coefficient / arb(TARGET_END)
    at_outer_threshold = coefficient / arb(B)
    required_cutoff = coefficient / WORKING_TARGET
    require(at_target_edge.lower() > WORKING_TARGET, "raw C2 bound unexpectedly closes at target edge")
    require(at_outer_threshold.lower() > WORKING_TARGET, "raw C2 bound unexpectedly closes at outer threshold")
    require(required_cutoff.lower() > arb("1e29"), "raw C2 obstruction scale drift")
    return {
        "inherited_tail_bound": "|K_t-K_(t,M)|<=C_AB/(2*pi^2*M)",
        "tail_coefficient_ball": coefficient.str(80, more=True),
        "bound_at_M_39894_ball": at_target_edge.str(80, more=True),
        "bound_at_M_5122421_ball": at_outer_threshold.str(80, more=True),
        "cutoff_ratio_needed_below_working_target_ball": required_cutoff.str(80, more=True),
        "working_target": WORKING_TARGET.str(40, more=True),
        "interpretation": "This proves only that the inherited raw C2 majorant cannot certify the required scale at either natural cutoff; it is not a lower bound for the actual defect.",
    }


def render_note(artifact: dict[str, Any]) -> str:
    obstruction = artifact["raw_norm_obstruction"]
    transitions = artifact["transition_certificate"]
    return f"""# Endpoint-completed phase-space projector defect

Date: 2026-08-13

Status: exact signed-projector reduction and transition audit; not a bound for
the compressed endpoint defect

Put `U=[0,L]`, `L={L}`, `T={{622,...,39894}}`, and

```text
K_x(u,m)=f_x(u)exp(-2*pi*i*m*u),
I_m(x)=integral_U K_x(u,m)du,
P_m(x)=Abel-integral_R K_x(u,m)du.                   (PD1)
```

For `M>=39894` and the common weight `w_m=exp(-pi*epsilon*m^2)`, the
Gamma-normalized finite residual integrand is exactly

```text
Delta_(M,epsilon)
 =H_x+sum_(m=-M)^M w_m I_m
     -sum_(m=622)^39894 w_m P_m.                    (PD2)
```

Let `F_M={{-M,...,M}}`.  Splitting the first sum by the target projector gives

```text
Delta_(M,epsilon)
 =H_x+sum_(m in F_M minus T) w_m I_m
     -sum_(m in T) w_m(P_m-I_m),                    (PD3)

F_M minus T=[-M,621] union [39895,M].                (PD4)
```

Thus the residual is one endpoint-completed, oriented phase-space rectangle:
the finite source interval over frequencies outside `T`, minus the two
full-line spatial tails over frequencies inside `T`, with the Poisson
half-current retained.  The full-line integrals and the limit in (PD1)--(PD3)
use the already-certified common Abel prescription.  No divergent piece is
split independently.

After the inherited physical projection and common limit,

```text
lim Delta_(M,epsilon)=R_KGamma=Q_K-G
 =E_Btr,win+E_outer+R_Dir.                           (PD5)
```

Therefore (PD3) is an exact representation of the whole Gamma-normalized
residual, not of `R_Dir` alone.  The latter is obtained only after subtracting
the certified B window and outer block.

The saved-height stationary audit identifies three occupancy interfaces:

```text
B face:                 {transitions['stationary_occupancy_interfaces']['B_face']},
A face:                 {transitions['stationary_occupancy_interfaces']['A_face']},
half-boundary corner:   {transitions['stationary_occupancy_interfaces']['half_boundary_corner']}. (PD6)
```

By contrast, `39694|39695` is only the certified ordinary/fold proof-ownership
boundary.  Both the A-face gap and the half-boundary gap retain the same sign
across it, so it is not an admissible place to split norms.  At the upper
corner, `4*39894=A-1` and `4*39895=A+3`; this is the lattice bracket of the
same A null face identified by the exact bi-Morse chart.

The old height-uniform `C^2` interchange majorant remains rigorous but is
quantitatively useless here:

```text
C_AB/(2*pi^2) = {obstruction['tail_coefficient_ball']},
[C_AB/(2*pi^2)]/39894 = {obstruction['bound_at_M_39894_ball']},
[C_AB/(2*pi^2)]/5122421 = {obstruction['bound_at_M_5122421_ball']},
M needed by that majorant for the working target
  > {obstruction['cutoff_ratio_needed_below_working_target_ball']}. (PD7)
```

This is not evidence that the actual projector defect is large.  It proves
that a raw derivative norm destroys too much cancellation.  The quantitative
theorem must estimate (PD3) as one signed object, use the exact hyperbolic
triangle coordinates away from (PD6), and install compatible B-face, A-face,
and half-corner charts before taking absolute values.

Pi provenance: `pi` in (PD1)--(PD7) comes from the equation-(9) Fourier
character, Gaussian Abel weight, Kummer phase, and the exact stationary
equations.  No fitted or geometric constant is introduced.

Proof boundary: exact finite/common-Abel projector algebra, mode counts,
saved-height transition classification, and failure of one inherited raw
majorant only.  No phase-adapted signed estimate, bound for `R_Dir`, complete
`Q_K-T` or `T_upper`, all-height theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["uniform_Gamma_insertion"]["decision"]["Dirichlet_target_kernel_is_exact_Gamma_bulk_insertion"] is True, "Gamma insertion dependency drift")
    require(dependencies["finite_rejoin"]["decision"]["finite_completed_B_block_rejoined_before_norms"] is True, "finite rejoin dependency drift")
    require(dependencies["symmetric_Poisson"]["decision"]["symmetric_poisson_interchange_proved"] is True, "Poisson dependency drift")
    require(dependencies["pair_triangle"]["decision"]["B_A_and_half_boundary_events_share_one_triangle_geometry"] is True, "triangle dependency drift")
    require(dependencies["bi_Morse"]["decision"]["A_fold_and_half_boundary_are_same_null_face"] is True, "bi-Morse dependency drift")
    require(dependencies["ownership"]["decision"]["mode_ownership_disjoint_and_exhaustive"] is True, "ownership dependency drift")

    artifact = {
        "kind": STEM,
        "status": "exact_endpoint_completed_phase_space_projector_defect_and_three_transition_interfaces_certified_quantitative_signed_estimate_open",
        "passed": True,
        "scope": {
            "height": HEIGHT,
            "source_interval": [0, L],
            "source_odd_roster": [A, B],
            "target_modes": [TARGET_START, TARGET_END],
            "target_count": TARGET_COUNT,
            "common_regulator": "w_m=exp(-pi*epsilon*m^2), then epsilon down to zero after the symmetric cutoff limit",
        },
        "symbolic_certificate": symbolic_certificate(),
        "mode_partition_certificate": mode_partition_certificate(),
        "transition_certificate": transition_certificate(),
        "raw_norm_obstruction": raw_norm_obstruction(dependencies["symmetric_Poisson"]),
        "decision": {
            "Gamma_normalized_residual_is_one_endpoint_completed_oriented_projector_defect": True,
            "target_full_line_carrier_subtracted_before_norms": True,
            "common_Abel_prescription_preserved": True,
            "three_saved_height_stationary_occupancy_interfaces_identified": True,
            "ordinary_fold_ownership_interface_is_not_stationary_transition": True,
            "A_null_face_and_half_boundary_corner_link_preserved": True,
            "raw_C2_interchange_majorant_closes_working_target": False,
            "phase_adapted_signed_projector_bound_proved": False,
            "compressed_R_Dir_target_proved": False,
            "complete_T_upper_proved": False,
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
            "flint_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Prove a cancellation-preserving quantitative estimate for the endpoint-completed signed projector defect (PD3). Use the exact hyperbolic pair-triangle chart on stationary-stable blocks, compatible local charts at B 621|622, A 39852|39853, and the half-boundary corner 39894|39895, and subtract the already-certified B window and outer block only after the common estimate. Do not norm at the bookkeeping split 39694|39695 or reuse the raw C2 tail majorant.",
        "proof_boundary": "Exact finite/common-Abel phase-space projector identity, outside-mode count, saved-height analytic-interface audit, and a quantitative obstruction to the inherited raw C2 majorant only. No phase-adapted signed estimate, R_Dir bound, complete Q_K-T or T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified endpoint-completed phase-space projector defect and transition audit", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
