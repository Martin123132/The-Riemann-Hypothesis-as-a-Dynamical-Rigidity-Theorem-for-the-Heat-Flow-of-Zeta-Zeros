#!/usr/bin/env python3
"""Restore the physical carrier and projector orientation at the A corner."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_corner_projector_orientation_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "corner": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_null_coordinate_affine_corner_evaluation_gate.json",
    "triangle": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_separable_triangle_geometry_gate.json",
    "projector": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_phase_space_projector_defect_gate.json",
    "paired_target": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_paired_target_residual_normal_form_gate.json",
    "universal_morse": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_universal_logistic_morse_characteristic_fold_reduction_gate.json",
}

HEIGHT = 10_000_000_000
A = 159_577
TARGET_END = 39_894
MODES = (39_894, 39_895)


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


def complex_record(value: acb, digits: int = 80) -> dict[str, str]:
    return {
        "real_ball": value.real.str(digits, more=True),
        "imag_ball": value.imag.str(digits, more=True),
        "absolute_ball": abs(value).str(digits, more=True),
    }


def symbolic_certificate() -> dict[str, str]:
    m, d = sp.symbols("m d", integer=True)
    endpoint = 2 * d + 1
    parity_exponent = sp.expand(m * endpoint - m**2)
    require(
        sp.expand(parity_exponent - (2 * m * d + m * (1 - m))) == 0,
        "odd-endpoint parity decomposition failed",
    )

    t, mode = sp.symbols("t m_pos", positive=True, real=True)
    D = sp.symbols("D", positive=True, real=True)
    r = t / (2 * sp.pi * mode**2)
    psi_star = -t / 2 - sp.pi * mode**2 + t * sp.log(r) / 2
    phi_star = sp.pi * mode * D
    base = t * (sp.log(t / (2 * sp.pi)) - 1) / 2 - t * sp.log(mode)
    require(
        sp.simplify(sp.expand_log(phi_star + psi_star - base, force=True) - sp.pi * mode * (D - mode)) == 0,
        "stationary carrier reduction failed",
    )

    chi = sp.symbols("chi")
    P_plus, P_minus, Q_plus, Q_minus = sp.symbols("P_plus P_minus Q_plus Q_minus")
    I_plus = P_plus + Q_plus
    I_minus = P_minus + Q_minus
    paired_residual = sp.expand(I_plus + I_minus - chi * P_plus)
    expected = sp.expand((1 - chi) * P_plus + P_minus + Q_plus + Q_minus)
    require(sp.expand(paired_residual - expected) == 0, "paired projector orientation failed")
    require(sp.diff(expected, Q_plus) == 1 and sp.diff(expected, Q_minus) == 1, "endpoint coefficient drift")

    return {
        "stationary_phase": "Phi_D,m=pi*m*D-pi*m^2-t/2+(t/2)log[t/(2pi*m^2)]",
        "odd_endpoint_parity": "For D=2d+1, m(D-m)=2md+m(1-m) is even, so exp(i*pi*m(D-m))=1.",
        "common_carrier": "exp(i Phi_D,m)=exp(i{t[log(t/(2pi))-1]/2-t log m}) for integer m and odd D",
        "paper_rotation": "exp(-i*pi/8)exp(i Phi_D,m)=exp(i[theta_0-t log m]), theta_0=t[log(t/(2pi))-1]/2-pi/8",
        "affine_endpoint_half_mode": "K_D,m^aff=epsilon_D exp(i Phi_D,m) W0_D,m C_aff,D,m, epsilon_A=-1",
        "physical_projection": "q_D,m^aff=2(pi/(32t))^(1/4) Re[exp(-i*pi/8)K_D,m^aff]",
        "paired_projector_mode": "I_m+I_-m-chi_T(m)P_m=(1-chi_T)P_m+P_-m+Q_m+Q_-m",
        "endpoint_orientation": "The paired endpoint term Q_m+Q_-m has coefficient +1 on both sides of the target edge; only the positive full-line bulk coefficient changes from 0 to 1.",
    }


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def interval_certificate(corner: dict[str, Any]) -> dict[str, Any]:
    ctx.dps = 110
    ctx.threads = 1
    pi, t, endpoint = arb.pi(), arb(HEIGHT), arb(A)
    imaginary = acb(0, 1)
    normalizer = (pi / (32 * t)) ** (arb(1) / 4)
    paper_rotation = (-imaginary * pi / 8).exp()
    source_rows = {row["mode"]: row for row in corner["rows"]}

    rows: list[dict[str, Any]] = []
    scalar_sum = arb(0)
    affine_sum = arb(0)
    correction_sum = arb(0)
    normalized_half_sum = acb(0)
    for mode_int in MODES:
        source = source_rows[mode_int]
        mode = arb(mode_int)
        r = t / (2 * pi * mode**2)
        K = pi * endpoint * mode
        W0 = acb(K, -1) * (4 * r ** arb("0.75") / (endpoint * (pi * t).sqrt()))
        phase = t * ((t / (2 * pi)).log() - 1) / 2 - t * mode.log()
        common_carrier = (imaginary * phase).exp()
        epsilon_A = arb(-1)

        scalar = parse_complex(source["canonical_wedge_C"])
        affine = parse_complex(source["affine_wedge_C_aff"])
        correction = affine - scalar
        raw_scalar = epsilon_A * common_carrier * W0 * scalar
        raw_affine = epsilon_A * common_carrier * W0 * affine
        raw_correction = epsilon_A * common_carrier * W0 * correction
        normalized_half = normalizer * raw_affine
        physical_scalar = 2 * normalizer * (paper_rotation * raw_scalar).real
        physical_affine = 2 * normalizer * (paper_rotation * raw_affine).real
        physical_correction = 2 * normalizer * (paper_rotation * raw_correction).real
        target_indicator = 1 if mode_int <= TARGET_END else 0
        scalar_sum += physical_scalar
        affine_sum += physical_affine
        correction_sum += physical_correction
        normalized_half_sum += normalized_half
        rows.append(
            {
                "endpoint": A,
                "mode": mode_int,
                "target_indicator": target_indicator,
                "positive_full_line_bulk_projector_coefficient": 1 - target_indicator,
                "paired_endpoint_projector_coefficient": 1,
                "common_stationary_carrier": complex_record(common_carrier),
                "W0": complex_record(W0),
                "normalized_affine_half_mode": complex_record(normalized_half),
                "physical_scalar_A_tangent_wedge_ball": physical_scalar.str(80, more=True),
                "physical_affine_A_tangent_wedge_ball": physical_affine.str(80, more=True),
                "physical_affine_correction_ball": physical_correction.str(80, more=True),
            }
        )

    require(rows[0]["target_indicator"] == 1 and rows[1]["target_indicator"] == 0, "target edge bracket drift")
    require(all(row["paired_endpoint_projector_coefficient"] == 1 for row in rows), "endpoint projector sign drift")
    require(all(arb(row["physical_affine_A_tangent_wedge_ball"]).upper() < 0 for row in rows), "corner physical signs drift")
    require(affine_sum.upper() < arb("-0.0077"), "two-mode affine A carrier lost its negative sign")
    require(correction_sum.upper() < arb("-6.3e-6"), "two-mode affine correction scale drift")
    return {
        "height": HEIGHT,
        "endpoint": A,
        "paper_normalizer_ball": normalizer.str(80, more=True),
        "rows": rows,
        "two_mode_normalized_half_sum": complex_record(normalized_half_sum),
        "two_mode_physical_scalar_A_tangent_wedge_ball": scalar_sum.str(80, more=True),
        "two_mode_physical_affine_A_tangent_wedge_ball": affine_sum.str(80, more=True),
        "two_mode_physical_affine_correction_ball": correction_sum.str(80, more=True),
        "interpretation": "The two affine tangent-wedge A terms reinforce after the exact carrier is restored. This is one signed endpoint chart, not the complete projector residual, and neither exact remainder has been included.",
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    rows = {row["mode"]: row for row in c["rows"]}
    return f"""# Affine A-corner carrier and projector orientation

Date: 2026-08-13

Status: exact phase/projector dictionary and rigorous two-mode tangent-wedge
projection; not an exact A-corner or global residual bound

The stationary phase of the endpoint-`D` triangle is

```text
Phi_D,m=pi*m*D-pi*m^2-t/2+(t/2)log[t/(2pi*m^2)].       (PO1)
```

For an odd endpoint `D=2d+1`,

```text
m(D-m)=2md+m(1-m) is even.                              (PO2)
```

Consequently both odd endpoints have the same exact carrier,

```text
exp(i Phi_D,m)
 =exp(i{{t[log(t/(2pi))-1]/2-t log m}}),

exp(-i*pi/8)exp(i Phi_D,m)
 =exp(i[theta_0-t log m]),
theta_0=t[log(t/(2pi))-1]/2-pi/8.                      (PO3)
```

Thus the affine A endpoint term and its equation-(9) projection are

```text
K_A,m^aff=-exp(i Phi_A,m)W0_A,m C_aff,A,m,
q_A,m^aff=2(pi/(32t))^(1/4)
 Re[exp(-i*pi/8)K_A,m^aff].                            (PO4)
```

The minus sign in (PO4) is `epsilon_A=-1`; it is not supplied by the target
projector.  Indeed, with `chi_T(m)` the positive target indicator,

```text
I_m+I_-m-chi_T(m)P_m
 =(1-chi_T)P_m+P_-m+Q_m+Q_-m.                         (PO5)
```

The paired endpoint term `Q_m+Q_-m` therefore has coefficient `+1` on both
sides of `39894|39895`.  Only the positive full-line bulk coefficient jumps
from zero at mode 39894 to one at mode 39895.

Restoring the common carrier and the paper normalizer gives

```text
mode 39894: q_A^aff={rows[39894]['physical_affine_A_tangent_wedge_ball']},
mode 39895: q_A^aff={rows[39895]['physical_affine_A_tangent_wedge_ball']},

two-mode affine sum={c['two_mode_physical_affine_A_tangent_wedge_ball']},
two-mode scalar sum={c['two_mode_physical_scalar_A_tangent_wedge_ball']},
affine-minus-scalar={c['two_mode_physical_affine_correction_ball']}. (PO6)
```

The two canonical affine A terms reinforce rather than cancel, and the
physical affine correction alone is more negative than `6.3e-6`.  That is
already comparable with the final `8.6e-6` target scale, so it cannot be
dropped.  The roughly `-0.00774` two-mode value is not an error estimate: it
must remain grouped with the rest of the A roster, the exact curved-face and
amplitude remainders, the full-line projector jump, and the other endpoint
sectors.

Pi provenance: `pi` in (PO1)--(PO6) comes from the exact equation-(9)
endpoint phase, Fourier parity, odd-square reflection, Gaussian
normalization, and paper prefactor.  No fitted constant is introduced.

Proof boundary: exact phase reduction, paired-projector orientation, and
rigorous normalized tangent-wedge values at modes 39894 and 39895 only.  No
exact A-corner remainder, complete A roster, `R_Dir` estimate, complete
`Q_K-T` or `T_upper`, all-height theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["triangle"]["decision"]["pair_Volterra_integral_reduced_to_separable_triangle"] is True, "triangle dependency drift")
    require(dependencies["projector"]["decision"]["Gamma_normalized_residual_is_one_endpoint_completed_oriented_projector_defect"] is True, "projector dependency drift")
    require(dependencies["paired_target"]["decision"]["target_carrier_subtracted_exactly_once"] is True, "target dependency drift")
    require(dependencies["universal_morse"]["decision"]["classical_inverse_sqrt_mode_carrier_recovered"] is True, "Morse carrier dependency drift")
    require(dependencies["corner"]["decision"]["canonical_affine_wedge_values_certified"] is True, "corner dependency drift")

    artifact = {
        "kind": STEM,
        "status": "exact_A_corner_common_phase_and_paired_projector_orientation_with_two_mode_physical_affine_values_certified",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(dependencies["corner"]),
        "decision": {
            "odd_endpoint_phase_reduced_to_common_classical_carrier": True,
            "A_endpoint_sign_identified_as_epsilon_A_minus_one": True,
            "paired_endpoint_coefficient_continuous_across_target_edge": True,
            "positive_full_line_bulk_projector_jump_identified": True,
            "equation9_normalization_and_half_projection_restored": True,
            "two_mode_affine_A_tangent_wedge_projection_certified": True,
            "two_mode_affine_correction_is_target_scale_significant": True,
            "exact_A_corner_remainder_bound_proved": False,
            "complete_A_roster_sum_proved": False,
            "R_Dir_bound_proved": False,
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
            "arb_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Evaluate the same carrier-oriented affine wedge in deterministic mode order over 39853..39936. Sum before taking norms and keep the positive full-line projector jump at 39894|39895 separate from the endpoint carrier. Then bound the exact curved-face and transformed-amplitude remainders in the null chart.",
        "proof_boundary": "Exact common carrier, paired-projector orientation, and normalized affine tangent-wedge values at modes 39894 and 39895 only. No exact A-corner remainder, complete A roster, R_Dir estimate, complete Q_K-T or T_upper, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified affine A-corner carrier and paired-projector orientation", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
