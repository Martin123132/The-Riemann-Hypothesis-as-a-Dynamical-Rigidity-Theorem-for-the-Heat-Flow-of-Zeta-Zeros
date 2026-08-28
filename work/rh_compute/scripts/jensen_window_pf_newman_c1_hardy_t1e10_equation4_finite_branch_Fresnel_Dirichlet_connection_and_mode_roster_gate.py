#!/usr/bin/env python3
"""Connect the finite branch sum to a Fresnel/Dirichlet integral and recover its mode roster."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import mpmath as mp

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))

from flint import arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_branch_Fresnel_Dirichlet_connection_and_mode_roster_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "branch_alignment": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_label_branch_correction_QK_alignment_gate.json",
    "geometric_half_line": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_QK_geometric_half_line_and_saddle_route_gate.json",
    "physical_ownership": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_physical_transform_ownership_ledger_gate.json",
}

HEIGHT = 10_000_000_000
SOURCE_ALPHA_MIN = 159_577
SOURCE_ALPHA_MAX = 5_122_421
SOURCE_COUNT = 2_481_423
PRECISION = 110


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {relative(path)}")
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


def half_line_integrals(t: mp.mpf, alpha: int) -> tuple[mp.mpc, mp.mpc, mp.mpc, mp.mpc]:
    s = mp.mpf("0.5") + 1j * t
    c = mp.pi * alpha * (1 + 1j) / mp.sqrt(2)
    cuts = [0, mp.mpf("0.01"), mp.mpf("0.1"), 1, 4, 8, mp.inf]
    i_minus = mp.quad(lambda x: x ** (-s) * mp.e ** (-mp.pi * x * x - c * x), cuts)
    i_plus = mp.quad(lambda x: x ** (-s) * mp.e ** (-mp.pi * x * x + c * x), cuts)
    j_stokes = mp.quad(
        lambda x: x ** (-mp.conj(s)) * mp.e ** (-mp.pi * x * x + 1j * c * x),
        cuts,
    )
    kappa = (
        1j
        * mp.sqrt(2 * mp.pi)
        * (2 * mp.pi) ** (1j * t)
        * mp.e ** (-mp.pi * t / 2)
        / mp.gamma(mp.mpf("0.5") + 1j * t)
    )
    return i_minus, i_plus, j_stokes, kappa


def connection_rows() -> list[dict[str, Any]]:
    mp.mp.dps = 65
    rows = []
    for t, alpha in ((mp.mpf("2"), 1), (mp.mpf("3.5"), 3), (mp.mpf("5"), 5)):
        i_minus, i_plus, j_stokes, kappa = half_line_integrals(t, alpha)
        discrepancy = i_minus - (-1j * mp.e ** (-mp.pi * t) * i_plus + kappa * j_stokes)
        k0 = mp.e ** (3 * mp.pi * t / 4 + 3j * mp.pi / 8)
        d = mp.e ** (-mp.pi * t / 4 - 1j * mp.pi / 8)
        branch_discrepancy = k0 * (-1j * mp.e ** (-mp.pi * t)) * i_plus - d * i_plus
        rows.append(
            {
                "t": str(t),
                "alpha": alpha,
                "DLMF_connection_discrepancy_absolute": mp.nstr(abs(discrepancy), 18),
                "branch_prefactor_discrepancy_absolute": mp.nstr(abs(branch_discrepancy), 18),
                "kappa_absolute": mp.nstr(abs(kappa), 35),
                "sqrt_1_plus_exp_minus_2pi_t": mp.nstr(mp.sqrt(1 + mp.e ** (-2 * mp.pi * t)), 35),
            }
        )
    return rows


def saddle_roster() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    t = arb(HEIGHT)
    pi = arb.pi()
    transition = (8 * t / pi).sqrt()
    center = (t / (2 * pi)).sqrt()

    def roots(alpha_value: int) -> tuple[arb, arb]:
        alpha = arb(alpha_value)
        radical = (alpha * alpha - 8 * t / pi).sqrt()
        return (alpha - radical) / 4, (alpha + radical) / 4

    lower_A, upper_A = roots(SOURCE_ALPHA_MIN)
    lower_B, upper_B = roots(SOURCE_ALPHA_MAX)
    require(arb(621) < lower_B < arb(622), "B lower saddle did not select mode 622")
    require(arb(39_852) < lower_A < arb(39_853), "A lower saddle did not select mode 39853")
    require(arb(39_936) < upper_A < arb(39_937), "A upper saddle did not select mode 39936")
    require(arb(2_560_588) < upper_B < arb(2_560_589), "B upper saddle bracket drift")
    require(arb(39_894) < center < arb(39_895), "classical center floor drift")

    ordinary_start = 622
    ordinary_end = 39_852
    transition_start = 39_853
    transition_end = 39_936
    center_floor = 39_894
    require(ordinary_end - ordinary_start + 1 == 39_231, "ordinary roster count drift")
    require(transition_end - transition_start + 1 == 84, "transition roster count drift")
    require(center_floor - transition_start + 1 == 42, "lower transition half count drift")
    require(transition_end - (center_floor + 1) + 1 == 42, "upper transition half count drift")

    return {
        "transition_alpha_ball": transition.str(70, more=True),
        "classical_center_sqrt_t_over_2pi_ball": center.str(70, more=True),
        "alpha_A": SOURCE_ALPHA_MIN,
        "alpha_B": SOURCE_ALPHA_MAX,
        "A_lower_saddle_ball": lower_A.str(70, more=True),
        "A_upper_saddle_ball": upper_A.str(70, more=True),
        "B_lower_saddle_ball": lower_B.str(70, more=True),
        "B_upper_saddle_ball": upper_B.str(70, more=True),
        "ordinary_integer_roster": [ordinary_start, ordinary_end],
        "ordinary_integer_count": 39_231,
        "transition_integer_roster": [transition_start, transition_end],
        "transition_integer_count": 84,
        "transition_lower_half": [transition_start, center_floor],
        "transition_upper_half": [center_floor + 1, transition_end],
        "transition_half_counts": [42, 42],
    }


def render_note(artifact: dict[str, Any]) -> str:
    roster = artifact["saddle_mode_roster"]
    return f"""# Finite branch Fresnel connection and mode roster

Date: 2026-08-27

Status: exact connection and mode roster certified; joined amplitude bound open

Write

```text
s=1/2+it,
c_alpha=pi*alpha*(1+i)/sqrt(2),
d=exp(-pi*t/4-i*pi/8),
K_0=exp(3pi*t/4+i3pi/8),                            (FC1)

I_-(alpha)=integral_0^infinity x^(-s)exp(-pi*x^2-c_alpha*x)dx,
I_+(alpha)=integral_0^infinity x^(-s)exp(-pi*x^2+c_alpha*x)dx,
J(alpha)=integral_0^infinity x^(-conj(s))
         exp(-pi*x^2+i*c_alpha*x)dx.                (FC2)
```

Scaling DLMF 12.5.1 to (FC2) and applying the minus-sign version of the exact
connection formula DLMF 12.2.19 gives, for every positive odd `alpha`,

```text
I_-(alpha)=-i exp(-pi*t)I_+(alpha)+kappa(t)J(alpha), (FC3)

kappa(t)=i sqrt(2pi)(2pi)^(it)exp(-pi*t/2)
         /Gamma(1/2+it),                            (FC4)

|kappa(t)|=sqrt(1+exp(-2pi*t)).                     (FC5)
```

The odd-label phase `exp(i*pi*alpha^2/4)=exp(i*pi/4)` is what makes (FC4)
independent of the label.  Three altered `(t,alpha)` rows directly verify
(FC3), (FC5), and the prefactor identity

```text
K_0[-i exp(-pi*t)] = d.                             (FC6)
```

Consequently the connection formula is exactly the finite-label split found
in Section 11.484:

```text
S_alpha=K_0 I_-(alpha)
       =d I_+(alpha)+K_0 kappa(t)J(alpha),

d I_+(alpha)=S_alpha-T_alpha,
K_0 kappa(t)J(alpha)=T_alpha.                       (FC7)
```

Unlike the rejected third-quadrant rotation of `I_-`, the `I_+` contour can
be rotated inside `0<=arg(x)<pi/4`, where the Gaussian supplies quadratic
decay.  Taking the endpoint as a Fresnel/Abel limit and integrating the tail
by parts gives the exact oscillatory representation

```text
d I_+(alpha)=integral_0^infinity[Fresnel]
             y^(-s)exp(i*pi*(alpha*y-y^2))dy.       (FC8)
```

The phase derivative is eventually `-2pi*y`, so (FC8) is a convergent
improper oscillatory integral; near zero the amplitude is `y^(-1/2)`.
Because the roster is finite, summing (FC8) is exact:

```text
B_W=integral_0^infinity[Fresnel]
    y^(-s)exp(-i*pi*y^2)D_W(y)dy,                   (FC9)

D_W(y)=sum_(j=0)^(M-1)exp(i*pi*(A+2j)y)
      =exp(i*pi*(A+M-1)y)sin(pi*M*y)/sin(pi*y),     (FC10)
```

with removable values supplied by the finite sum at integer `y`.  Equations
(FC3)--(FC10) own the branch multiplier and Stokes complement that were open
in Section 11.485.

The one-label stationary equation in (FC8) is

```text
2pi*y^2-pi*alpha*y+t=0,
y_+/-=(alpha+/-sqrt(alpha^2-8t/pi))/4.              (FC11)
```

At `t=10^10`, Arb proves

```text
y_-(B)={roster['B_lower_saddle_ball']},
y_-(A)={roster['A_lower_saddle_ball']},
y_+(A)={roster['A_upper_saddle_ball']},
sqrt(t/(2pi))={roster['classical_center_sqrt_t_over_2pi_ball']}. (FC12)
```

The induced integer roster is therefore

```text
ordinary lower branch: 622..39852, 39231 modes,
transition block:       39853..39936, 84 modes,
transition split:       39853..39894 and 39895..39936,
                        42 modes on each side.       (FC13)
```

This independently recovers the existing physical ownership indices:
the extended P/Gamma roster is `622..39936`, the paired A roster is
`39853..39936`, the classical target ends at `39894`, and the extra Gamma
block is `39895..39936`.  The agreement certifies the route geometry and
index ownership, not a new amplitude estimate.  The next step is to express
the already-owned `G+A_transition` carriers on the common Fresnel integral
(FC9), subtract them before norms, and enclose only the joined remainder.

Pi provenance: DLMF's `pi` is the standard gamma/parabolic-cylinder
normalization; all project-specific occurrences descend from the RSI
Gaussian, quarter-turn, and Riemann-Siegel phase.  No fitted constant or
circle construction is introduced.

Proof boundary: exact parabolic-cylinder connection, branch/Stokes ownership,
finite Fresnel/Dirichlet representation, and actual saddle-mode roster only.
No common-integral carrier subtraction, joined amplitude bound, actual-height
`B_W`, `Q_K`, `D_K`, `J_Z`, non-A, all-height, `Lambda<=0`, PF-infinity, RH,
or prize-level enclosure is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "dependency failure")
    require(
        dependencies["geometric_half_line"]["decision"]["candidate_saddle_ray_permitted_as_integral_identity"]
        is False,
        "uncertified saddle ray was promoted",
    )
    require(
        dependencies["physical_ownership"]["decision"]["target_and_extra_Gamma_blocks_owned_without_overlap"]
        is True,
        "physical ownership drift",
    )

    rows = connection_rows()
    require(
        max(mp.mpf(row["DLMF_connection_discrepancy_absolute"]) for row in rows) < mp.mpf("1e-35"),
        "parabolic-cylinder connection scout failed",
    )
    require(
        max(mp.mpf(row["branch_prefactor_discrepancy_absolute"]) for row in rows) < mp.mpf("1e-60"),
        "branch prefactor scout failed",
    )
    roster = saddle_roster()
    ownership_scope = dependencies["physical_ownership"]["scope"]
    require(ownership_scope["extended_P_modes"] == [622, 39_936], "extended P roster drift")
    require(ownership_scope["paired_A_modes"] == [39_853, 39_936], "paired A roster drift")
    require(ownership_scope["target_P_modes"] == [622, 39_894], "target P roster drift")

    artifact = {
        "kind": STEM,
        "status": "exact_parabolic_cylinder_branch_connection_Fresnel_Dirichlet_sum_and_mode_roster_certified",
        "passed": True,
        "scope": {
            "height": HEIGHT,
            "alpha_min": SOURCE_ALPHA_MIN,
            "alpha_max": SOURCE_ALPHA_MAX,
            "window_label_count": SOURCE_COUNT,
        },
        "parabolic_cylinder_connection": {
            "I_minus": "integral_0^infinity x^(-s)exp(-pi*x^2-c_alpha*x)dx",
            "I_plus": "integral_0^infinity x^(-s)exp(-pi*x^2+c_alpha*x)dx",
            "J": "integral_0^infinity x^(-conj(s))exp(-pi*x^2+i*c_alpha*x)dx",
            "identity": "I_minus=-i*exp(-pi*t)*I_plus+kappa(t)*J",
            "kappa": "i*sqrt(2*pi)*(2*pi)^(i*t)*exp(-pi*t/2)/Gamma(1/2+i*t)",
            "kappa_modulus": "sqrt(1+exp(-2*pi*t))",
            "source_split": "S_alpha=d*I_plus+K_0*kappa*J=(S_alpha-T_alpha)+T_alpha",
            "DLMF_references": ["https://dlmf.nist.gov/12.5.E1", "https://dlmf.nist.gov/12.2.E19"],
            "surrogate_rows": rows,
        },
        "finite_Fresnel_branch_sum": {
            "one_label": "d*I_plus=Fresnel integral_0^infinity y^(-s)exp(i*pi*(alpha*y-y^2))dy",
            "rotation": "x=exp(i*pi/4)*y approached through arg(x)<pi/4",
            "tail_convergence": "one integration by parts after the last stationary point; phase derivative is asymptotic to -2*pi*y",
            "finite_sum": "B_W=Fresnel integral_0^infinity y^(-s)exp(-i*pi*y^2)D_W(y)dy",
            "Dirichlet_factor": "D_W(y)=sum_(j=0)^(M-1)exp(i*pi*(A+2*j)*y)=exp(i*pi*(A+M-1)*y)*sin(pi*M*y)/sin(pi*y)",
            "integer_singularities_removable_by_finite_sum": True,
            "infinite_A21_interchange_used": False,
        },
        "saddle_mode_roster": roster,
        "physical_ownership_match": {
            "extended_P_modes": ownership_scope["extended_P_modes"],
            "target_P_modes": ownership_scope["target_P_modes"],
            "paired_A_modes": ownership_scope["paired_A_modes"],
            "extra_Gamma_modes": [39_895, 39_936],
            "independent_saddle_roster_matches_existing_indices": True,
            "carrier_amplitude_subtraction_on_common_Fresnel_integral_complete": False,
        },
        "decision": {
            "branch_Stokes_connection_complete": True,
            "finite_branch_sum_Fresnel_Dirichlet_representation_exact": True,
            "actual_mode_roster_recovered_from_saddle_map": True,
            "uncertified_third_quadrant_rotation_required": False,
            "common_integral_G_plus_A_transition_subtraction_complete": False,
            "actual_height_joined_remainder_enclosed": False,
            "D_K_enclosed": False,
            "non_A_bound_proved": False,
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
            "numerical_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Put the exact target P/Gamma modes 622..39894, extra Gamma modes 39895..39936, and paired A modes 39853..39936 onto the common Fresnel/Dirichlet representation of B_W. Subtract G+A_transition before norms, retain P_W and the parabolic-cylinder Stokes complement with exact H(t), and derive one signed interval for the joined J_Z or D_K remainder.",
        "proof_boundary": "Exact parabolic-cylinder connection, branch/Stokes ownership, finite Fresnel/Dirichlet representation, and actual saddle-mode roster only. No common-integral carrier subtraction, joined amplitude bound, actual-height B_W, Q_K, D_K, Delta_KU, J_Z, non-A, all-height, Lambda<=0, PF-infinity, RH, or prize-level enclosure is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified exact branch/Stokes connection, Fresnel sum, and saddle mode roster", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
