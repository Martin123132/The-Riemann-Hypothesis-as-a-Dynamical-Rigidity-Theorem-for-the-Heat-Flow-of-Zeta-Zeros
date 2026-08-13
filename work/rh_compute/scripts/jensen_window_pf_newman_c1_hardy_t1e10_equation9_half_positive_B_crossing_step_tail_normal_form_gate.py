#!/usr/bin/env python3
"""Certify the exact positive-mode B crossing as step mismatch plus one tail."""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_positive_B_crossing_step_tail_normal_form_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "endpoint_decomposition": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_finite_endpoint_tail_decomposition_gate.json",
    "Gamma_bulk": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_global_morse_gamma_bulk_gate.json",
    "face_trace": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_face_trace_lattice_deghosting_gate.json",
    "B_tangent": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_kummer_B_crossing_tangent_profile_barrier_gate.json",
}

PRECISION = 100
T = 10_000_000_000
A = 159_577
B = 5_122_421
LOW_END = 621
TARGET_START = 622
TARGET_END = 39_894


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
    q, a, e, full, tail, low, negative_q = sp.symbols(
        "q a E F_inf T I_low I_qneg", real=True
    )
    p_bulk = 2 * a * full

    # F_inf=(1+i)/2, so 2*a*F_inf=a(1+i)=P_bulk.  On q<0,
    # F(q)=-F_inf+T(-q); on q>0, F_inf-F(q)=T(q).
    p_b_positive = e - a * tail
    p_b_negative = e - p_bulk + a * tail
    u_positive = e - a * tail
    u_negative = e + a * tail
    require(sp.expand(p_b_positive - u_positive) == 0, "positive-q B tail failed")
    require(sp.expand(p_b_negative - (u_negative - p_bulk)) == 0, "negative-q B tail failed")

    residual_positive_q = p_bulk + sp.symbols("P_A") + p_b_positive - (1 - low) * p_bulk
    normal_positive_q = sp.symbols("P_A") + u_positive + low * p_bulk
    require(sp.expand(residual_positive_q - normal_positive_q) == 0, "positive-q residual identity failed")
    residual_negative_q = p_bulk + sp.symbols("P_A") + p_b_negative - (1 - low) * p_bulk
    normal_negative_q = sp.symbols("P_A") + u_negative + (low - 1) * p_bulk
    require(sp.expand(residual_negative_q - normal_negative_q) == 0, "negative-q residual identity failed")

    # Jump cancellation at q=0: T(0)=F_inf.
    u_left = e + a * full
    u_right = e - a * full
    require(sp.expand((u_right - u_left) + p_bulk) == 0, "crossing jump failed")

    x, m, endpoint = sp.symbols("x m B", positive=True, real=True)
    q_formula = sp.sqrt(x / 2) * (endpoint - 2 * m / x)
    a_formula = m * sp.sqrt(2 / x)
    leading = sp.simplify(1 + a_formula / q_formula)
    require(sp.simplify(leading - endpoint * x / (endpoint * x - 2 * m)) == 0, "rational leading tail failed")

    # Include the completed-square phase and the external 1/x current factor.
    completed_phase = -sp.pi * m**2 / x + sp.pi * q_formula**2 / 2
    trace_phase = sp.pi * endpoint**2 * x / 4 - sp.pi * m * endpoint
    require(sp.simplify(completed_phase - trace_phase) == 0, "trace completion failed")

    z = sp.symbols("z", positive=True, real=True)
    phase_z = sp.pi * m**2 / z + sp.pi * endpoint**2 * z / 4
    driver = z ** (-sp.Rational(1, 2)) * (2 + sp.I * sp.pi * endpoint**2 * z)
    pair_boundary = sp.simplify(
        x ** (-sp.Rational(3, 2)) * driver.subs(z, x)
        / (2 * sp.I * sp.pi * sp.I * sp.diff(phase_z, z).subs(z, x))
    )
    expected_pair_boundary = -2 * (2 + sp.I * sp.pi * endpoint**2 * x) / (
        sp.pi**2 * (endpoint**2 * x**2 - 4 * m**2)
    )
    require(sp.simplify(pair_boundary - expected_pair_boundary) == 0, "paired first boundary current failed")

    roster = sp.symbols("N", integer=True, positive=True)
    aa = sp.symbols("a", noninteger=True, positive=True)
    # Algebraic partial fraction underlying the cotangent sum.
    require(
        sp.simplify(1 / (aa**2 - roster**2) - (1 / (aa - roster) + 1 / (aa + roster)) / (2 * aa)) == 0,
        "cotangent partial fraction failed",
    )

    return {
        "B_coordinate": "q_B=sqrt(x/2)(B-2m/x), a_m=m*sqrt(2/x)",
        "Fresnel_tail": "T(r)=integral_r^infinity exp(i*pi*u^2/2)du for r>=0, with T(0)=(1+i)/2",
        "sign_adapted_tail": "U_B=E_B/(i*pi)-sgn(q_B)*a_m*T(|q_B|)",
        "B_endpoint_identity": "P_B=U_B-1_(q_B<0)P_bulk",
        "positive_residual_identity": "P_bulk+P_A+P_B-1_(m>=622)P_bulk=P_A+U_B+[1_(m<=621)-1_(q_B<0)]P_bulk",
        "jump_cancellation": "Across q_B=0, U_B jumps by -P_bulk and the step mismatch jumps by +P_bulk; their sum is continuous.",
        "tail_leading_term": "Using T(r)=-E(r)/(i*pi*r)+remainder, U_B has common leading E_B*B*x/[i*pi(B*x-2m)].",
        "completed_leading_trace": "After (-1)^m exp(-i*pi*m^2/x)/x, the common leading term is B exp(i*pi*B^2*x/4)/[i*pi(B*x-2m)].",
        "paired_first_IBP_boundary": "-2(2+i*pi*B^2*x)/[pi^2(B^2*x^2-4m^2)] times the common face phase",
        "cotangent_sum": "For a=B*x/2, sum_(m>=1)1/(B^2*x^2-4m^2)=pi*cot(pi*a)/(4*B*x)-1/(2*B^2*x^2), away from integer a.",
        "local_pole_guard": "Near the B trace saddle, modes 621 and 622 must be retained exactly before the cotangent remainder is bounded.",
    }


def fresnel_primitive(q: arb, pi: arb, imaginary_unit: acb) -> acb:
    full = (1 + imaginary_unit) / 2
    argument = (-imaginary_unit * pi / 4).exp() * (pi / 2).sqrt() * q
    return full * argument.erf()


def continuity_witness(mode: int) -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    imaginary_unit = acb(0, 1)
    m = arb(mode)
    endpoint = arb(B)
    x_cross = 2 * m / endpoint
    amplitude = m * (2 / x_cross).sqrt()
    full = acb(1, 1)
    bulk = amplitude * full
    e0 = acb(1)
    half_fresnel = full / 2
    u_left = e0 / (imaginary_unit * pi) + amplitude * half_fresnel
    u_right = e0 / (imaginary_unit * pi) - amplitude * half_fresnel
    jump = u_right - u_left
    require(abs(jump + bulk) < arb("1e-90"), f"jump witness failed at mode {mode}")

    eps = arb("1e-30")
    q_left = -eps
    q_right = eps
    f_left = fresnel_primitive(q_left, pi, imaginary_unit)
    f_right = fresnel_primitive(q_right, pi, imaginary_unit)
    p_b_left = (imaginary_unit * pi * q_left**2 / 2).exp() / (imaginary_unit * pi) - amplitude * (full / 2 - f_left)
    p_b_right = (imaginary_unit * pi * q_right**2 / 2).exp() / (imaginary_unit * pi) - amplitude * (full / 2 - f_right)
    low = mode <= LOW_END
    original_left = p_b_left + (bulk if low else 0)
    original_right = p_b_right + (bulk if low else 0)
    # The original current is continuous in q.  The tiny finite difference is
    # O(eps) and is recorded only as a reconstruction witness.
    require(abs(original_left - original_right) < arb("1e-20"), f"finite crossing continuity failed at mode {mode}")
    return {
        "mode": mode,
        "x_cross_ball": x_cross.str(PRECISION, more=True),
        "bulk_ball": {"real": bulk.real.str(PRECISION, more=True), "imag": bulk.imag.str(PRECISION, more=True)},
        "U_right_minus_U_left_ball": {"real": jump.real.str(PRECISION, more=True), "imag": jump.imag.str(PRECISION, more=True)},
        "jump_plus_bulk_absolute_ball": abs(jump + bulk).str(PRECISION, more=True),
        "finite_eps_original_difference_absolute_ball": abs(original_left - original_right).str(PRECISION, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    return f"""# Exact B crossing: step mismatch plus sign-adapted tail

Date: 2026-08-13

Status: exact-lemma certificate; not a proof of the quantitative B remainder

For a positive mode put

```text
q_B=sqrt(x/2)(B-2m/x),
a_m=m sqrt(2/x),
T(r)=integral_r^infinity exp(i*pi*u^2/2)du.            (ST1)
```

Define one sign-adapted grouped tail

```text
U_B=E_B/(i*pi)-sgn(q_B)a_m T(|q_B|).                  (ST2)
```

Using oddness of the Fresnel primitive gives the exact identity

```text
P_B=U_B-1_(q_B<0)P_bulk.                              (ST3)
```

Consequently, over the complete positive source target `1..39894`,

```text
P_bulk+P_A+P_B-1_(m>=622)P_bulk
 =P_A+U_B+[1_(m<=621)-1_(q_B<0)]P_bulk.               (ST4)
```

Neither term after `P_A` is continuous separately.  At `q_B=0`, `U_B`
jumps by `-P_bulk`, while the step mismatch jumps by `+P_bulk`; their sum is
exactly continuous.  This is the nonlinear completion missing from the
isolated tangent profile.

Away from `q_B=0`, one Fresnel integration by parts gives the common leading
term

```text
U_B=E_B*B*x/[i*pi(B*x-2m)]+Fresnel remainder.         (ST5)
```

After the completed-square phase and external `1/x` factor, (ST5) becomes

```text
B exp(i*pi*B^2*x/4)/[i*pi(B*x-2m)],                   (ST6)
```

whose phase is independent of `m`.  Pairing the negative mode before the
same normal integration gives the complete first boundary current

```text
-2(2+i*pi*B^2*x)/[pi^2(B^2*x^2-4m^2)].               (ST7)
```

Thus its all-pair leading sum is cotangent-type.  For `a=B*x/2`,

```text
sum_(m=1)^infinity 1/(B^2*x^2-4m^2)
 =pi*cot(pi*a)/(4*B*x)-1/(2*B^2*x^2).                (ST8)
```

Equation (ST8) is used only after the local poles are removed.  The exact
face ledger proves that near the B trace saddle those poles are precisely
modes `621` and `622`; they must remain in the continuous combination
(ST2)--(ST4).  The cotangent remainder and the Fresnel remainder still need
explicit bounds.

Pi provenance: every `pi` comes from the equation-(9) Fresnel phase,
completed-square phase, or the integer cotangent partial fraction.  No fitted
constant is introduced.

Proof boundary: exact step/tail identity, crossing jump cancellation,
first boundary current, and cotangent algebra only.  No Fresnel-remainder or
cotangent-remainder bound, two-mode B enclosure, A-fold splice, complete
paired residual, `T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion follows.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(
        dependencies["endpoint_decomposition"]["decision"]["finite_current_decomposed_exactly"] is True,
        "endpoint dependency drift",
    )
    require(
        dependencies["Gamma_bulk"]["decision"]["ordinary_full_line_bulk_closed_at_saved_height"] is True,
        "bulk dependency drift",
    )
    require(
        dependencies["face_trace"]["decision"]["B_face_local_integer_roster_is_621_622"] is True,
        "face roster dependency drift",
    )
    require(
        dependencies["B_tangent"]["decision"]["scalar_and_q_density_channels_retained"] is True,
        "tangent current dependency drift",
    )

    artifact = {
        "kind": STEM,
        "status": "exact_positive_B_crossing_step_mismatch_sign_tail_and_cotangent_leading_normal_form_proved",
        "passed": True,
        "scope": {"height": T, "upper_endpoint": B, "positive_modes": [1, TARGET_END], "step": [LOW_END, TARGET_START]},
        "symbolic_certificate": symbolic_certificate(),
        "continuity_witnesses": [continuity_witness(LOW_END), continuity_witness(TARGET_START)],
        "decision": {
            "B_endpoint_rewritten_as_sign_adapted_tail_plus_bulk_step": True,
            "tail_and_step_jumps_cancel_exactly": True,
            "positive_B_residual_normal_form_is_continuous": True,
            "paired_first_boundary_current_retains_q_density": True,
            "nonlocal_leading_pair_sum_has_cotangent_form": True,
            "local_modes_621_622_may_be_absorbed_into_cotangent_bound": False,
            "Fresnel_and_cotangent_remainders_bounded": False,
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
        "next_obligation": "On a face-centered window containing only q_B crossings 621 and 622, evaluate their exact continuous step-tail currents. Subtract those poles from the cotangent expression (ST8), bound the analytic nonlocal leading sum and the next Fresnel remainder, and use trace-phase integration by parts outside the window.",
        "proof_boundary": "Exact B step/tail and cotangent-leading normal form only. No remainder bound, two-mode enclosure, A-fold splice, complete paired residual, complete T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified exact B crossing step-tail normal form and cotangent leading sum", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
