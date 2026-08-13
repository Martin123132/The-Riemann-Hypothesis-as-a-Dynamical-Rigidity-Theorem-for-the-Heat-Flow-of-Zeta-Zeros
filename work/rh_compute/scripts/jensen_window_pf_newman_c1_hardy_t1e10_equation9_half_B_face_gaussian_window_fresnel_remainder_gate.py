#!/usr/bin/env python3
"""Certify the two-crossing B window and its summable Fresnel remainder."""

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

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_gaussian_window_fresnel_remainder_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "face_trace": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_face_trace_lattice_deghosting_gate.json",
    "step_tail": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_positive_B_crossing_step_tail_normal_form_gate.json",
}

PRECISION = 100
T = 10_000_000_000
B = 5_122_421
WINDOW_XI = 70
LOCAL_MODES = (621, 622)
FINITE_END = 39_894


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
    x, endpoint, mode = sp.symbols("x B m", positive=True, real=True)
    c = endpoint * x / 2
    q = sp.sqrt(2 / x) * (c - mode)
    amplitude = mode * sp.sqrt(2 / x)

    first = sp.simplify((1 + amplitude / q) / (sp.I * sp.pi * x))
    first_expected = endpoint / (sp.I * sp.pi * (endpoint * x - 2 * mode))
    require(sp.simplify(first - first_expected) == 0, "first completed tail term failed")

    second = sp.simplify(-amplitude / (sp.pi**2 * x * q**3))
    second_expected = -mode / (2 * sp.pi**2 * (c - mode) ** 3)
    require(sp.simplify(second - second_expected) == 0, "second completed tail term failed")

    third = sp.simplify(3 * sp.I * amplitude / (sp.pi**3 * x * q**5))
    third_expected = 3 * sp.I * mode * x / (4 * sp.pi**3 * (c - mode) ** 5)
    require(sp.simplify(third - third_expected) == 0, "third completed tail term failed")

    remainder = sp.simplify(30 * amplitude / (sp.pi**4 * x * abs(q) ** 7))
    remainder_expected = 15 * mode * x**2 / (4 * sp.pi**4 * abs(c - mode) ** 7)
    require(sp.simplify(remainder - remainder_expected) == 0, "completed remainder bound failed")

    return {
        "trace_scale": "xi=sqrt(H_B)(x-x_B^-), H_B=Phi_B''(x_B^-)",
        "crossing_scale": "xi_m=(2*sqrt(H_B)/B)(m-m_B^-)",
        "Fresnel_recurrence": "J_k(r)=-E(r)/(i*pi*r^(2k+1))+(2k+1)J_(k+1)(r)/(i*pi)",
        "tail_expansion": "T(r)=iE/(pi*r)+E/(pi^2*r^3)-3iE/(pi^3*r^5)+R_T, |R_T|<=30/(pi^4*r^7)",
        "signed_tail_expansion": "U_B=E/(i*pi)(1+a_m/q_B)-a_m*E/(pi^2*q_B^3)+3i*a_m*E/(pi^3*q_B^5)+R_U",
        "signed_tail_remainder": "|R_U|<=30*a_m/(pi^4*|q_B|^7)",
        "completed_terms": "After parity, exp(-i*pi*m^2/x), and 1/x: B/[i*pi(Bx-2m)]-m/[2*pi^2(c-m)^3]+3i*m*x/[4*pi^3(c-m)^5], c=Bx/2",
        "completed_remainder": "|R_ext,m|<=15*m*x^2/[4*pi^4*|c-m|^7]",
        "summation_guard": "Only R_ext is absolutely summed here. The three rational leading terms and local modes 621,622 remain cancellation-coupled.",
    }


def interval_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    t = arb(T)
    endpoint = arb(B)
    pi = arb.pi()
    x0 = (1 - (1 - 8 * t / (pi * endpoint**2)).sqrt()) / 2
    hessian = t * (1 - 2 * x0) / (2 * x0**2 * (1 - x0) ** 2)
    sqrt_hessian = hessian.sqrt()
    m0 = endpoint * x0 / 2
    lattice_spacing = 2 * sqrt_hessian / endpoint
    x_low = x0 - arb(WINDOW_XI) / sqrt_hessian
    x_high = x0 + arb(WINDOW_XI) / sqrt_hessian
    c_low = endpoint * x_low / 2
    c_high = endpoint * x_high / 2

    def crossing_xi(mode: int) -> arb:
        return lattice_spacing * (arb(mode) - m0)

    xi = {mode: crossing_xi(mode) for mode in (620, 621, 622, 623)}
    require(xi[620] < -arb(WINDOW_XI), "mode 620 enters the B window")
    require(-arb(WINDOW_XI) < xi[621] < arb(WINDOW_XI), "mode 621 leaves the B window")
    require(-arb(WINDOW_XI) < xi[622] < arb(WINDOW_XI), "mode 622 leaves the B window")
    require(xi[623] > arb(WINDOW_XI), "mode 623 enters the B window")

    q_620_min = (endpoint * x_low - 2 * arb(620)) / (2 * x_low).sqrt()
    q_623_min = abs((endpoint * x_high - 2 * arb(623)) / (2 * x_high).sqrt())
    require(q_620_min > arb(85), "left nonlocal q separation below 85")
    require(q_623_min > arb(75), "right nonlocal q separation below 75")

    weighted_sum = arb(0)
    for mode in range(1, FINITE_END + 1):
        if mode in LOCAL_MODES:
            continue
        distance = c_low - mode if mode <= 620 else arb(mode) - c_high
        require(distance > 0, f"nonlocal distance is not positive at mode {mode}")
        weighted_sum += arb(mode) / distance**7

    first_outer = FINITE_END + 1
    y0 = arb(first_outer) - c_high
    outer_tail = (
        arb(first_outer) / y0**7
        + 1 / (5 * y0**5)
        + c_high / (6 * y0**6)
    )
    infinite_weighted_sum = weighted_sum + outer_tail
    finite_remainder = 15 * x_high**2 * weighted_sum / (4 * pi**4)
    infinite_remainder = 15 * x_high**2 * infinite_weighted_sum / (4 * pi**4)
    require(finite_remainder < arb("7.48e-6"), "finite nonlocal remainder exceeds 7.48e-6")
    require(infinite_remainder < arb("7.48e-6"), "infinite nonlocal remainder exceeds 7.48e-6")

    return {
        "height": T,
        "endpoint": B,
        "window_xi": WINDOW_XI,
        "x_trace_ball": x0.str(PRECISION, more=True),
        "trace_hessian_ball": hessian.str(PRECISION, more=True),
        "sqrt_trace_hessian_ball": sqrt_hessian.str(PRECISION, more=True),
        "lattice_center_ball": m0.str(PRECISION, more=True),
        "standardized_lattice_spacing_ball": lattice_spacing.str(PRECISION, more=True),
        "x_window_low_ball": x_low.str(PRECISION, more=True),
        "x_window_high_ball": x_high.str(PRECISION, more=True),
        "c_window_low_ball": c_low.str(PRECISION, more=True),
        "c_window_high_ball": c_high.str(PRECISION, more=True),
        "crossing_xi_balls": {str(mode): value.str(PRECISION, more=True) for mode, value in xi.items()},
        "local_crossing_modes": list(LOCAL_MODES),
        "mode_620_q_min_ball": q_620_min.str(PRECISION, more=True),
        "mode_623_q_abs_min_ball": q_623_min.str(PRECISION, more=True),
        "all_nonlocal_positive_modes_q_abs_lower_bound": "75",
        "finite_nonlocal_weighted_sum_ball": weighted_sum.str(PRECISION, more=True),
        "finite_nonlocal_completed_remainder_ball": finite_remainder.str(PRECISION, more=True),
        "outer_tail_weighted_sum_upper_ball": outer_tail.str(PRECISION, more=True),
        "infinite_nonlocal_completed_remainder_ball": infinite_remainder.str(PRECISION, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    xi = c["crossing_xi_balls"]
    return f"""# Gaussian B-face window and summable Fresnel remainder

Date: 2026-08-13

Status: exact asymptotic-remainder lemma plus saved-height interval
certificate; not a proof of the complete B-face estimate

Use the exact lower B trace saddle and scale

```text
xi=sqrt(H_B)(x-x_B^-),
H_B=Phi_B''(x_B^-),
ell_B=2sqrt(H_B)/B.                                   (GW1)
```

At `t=10^10`, interval arithmetic gives

```text
H_B={c['trace_hessian_ball']},
ell_B={c['standardized_lattice_spacing_ball']}.        (GW2)
```

The crossing center of mode `m` is
`xi_m=ell_B(m-m_B^-)`.  The four nearest relevant centers are

```text
xi_620={xi['620']},
xi_621={xi['621']},
xi_622={xi['622']},
xi_623={xi['623']}.                                    (GW3)
```

Therefore the closed window `|xi|<=70` contains exactly the B crossings
`621` and `622`.  Throughout it, every other positive mode has

```text
|q_B|>75,                                             (GW4)
```

with the nearest left and right margins certified respectively by
`{c['mode_620_q_min_ball']}` and
`{c['mode_623_q_abs_min_ball']}`.

For `E(r)=exp(i*pi*r^2/2)`, repeated integration by parts gives

```text
T(r)=iE/(pi*r)+E/(pi^2*r^3)-3iE/(pi^3*r^5)+R_T,
|R_T|<=30/(pi^4*r^7).                                 (GW5)
```

Substitution into the sign-adapted current from Section 11.362, followed by
the exact completed phase, parity, and external `1/x` factor, gives

```text
B/[i*pi(Bx-2m)]
-m/[2*pi^2(c-m)^3]
+3i*m*x/[4*pi^3(c-m)^5] + R_ext,m,  c=Bx/2,          (GW6)

|R_ext,m|<=15*m*x^2/[4*pi^4*|c-m|^7].                (GW7)
```

Unlike the leading rational terms, (GW7) is absolutely summable.  Uniformly
on `|xi|<=70`, summing every positive nonlocal mode gives

```text
sum_(1<=m<=39894, m notin {{621,622}})|R_ext,m|
 <= {c['finite_nonlocal_completed_remainder_ball']},

sum_(m>=1, m notin {{621,622}})|R_ext,m|
 <= {c['infinite_nonlocal_completed_remainder_ball']} < 7.48e-6.   (GW8)
```

The infinite tail in (GW8) is bounded by monotone integral comparison; it is
not numerically truncated.

This closes the nonlocal third-order Fresnel remainder on the natural B
window.  It deliberately does not sum the three pole-bearing terms in
(GW6) separately.  Those must be pole-subtracted with modes 621 and 622 and
combined with the exact step current, negative/outer completion, and the
grouped endpoint difference.

Pi provenance: every `pi` is inherited from the equation-(9) Fresnel phase,
the exact face trace, or integration by parts.  No fitted constant is used.

Proof boundary: exact window roster, nonlocal `q_B` separation, three-term
Fresnel expansion, and its absolutely summed remainder only.  The three
rational lattice sums, exact two-mode current, outside-window tail, and
A-fold splice remain open.  This result does not establish the
complete paired residual, `T_upper`, `Lambda<=0`, PF-infinity, RH, or a
prize-level conclusion.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(
        dependencies["face_trace"]["decision"]["B_face_local_integer_roster_is_621_622"] is True,
        "face-roster dependency drift",
    )
    require(
        dependencies["step_tail"]["decision"]["B_endpoint_rewritten_as_sign_adapted_tail_plus_bulk_step"] is True,
        "step-tail dependency drift",
    )

    artifact = {
        "kind": STEM,
        "status": "B_face_two_crossing_window_and_nonlocal_three_term_Fresnel_remainder_certified",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(),
        "decision": {
            "xi_70_window_contains_exactly_B_crossings_621_622": True,
            "all_nonlocal_positive_modes_have_q_abs_above_75": True,
            "three_term_Fresnel_tail_with_explicit_remainder_proved": True,
            "nonlocal_completed_remainder_absolutely_summable": True,
            "infinite_nonlocal_completed_remainder_below_7p48e_minus_6": True,
            "leading_rational_terms_summed_without_pole_subtraction": False,
            "exact_local_modes_621_622_enclosed": False,
            "complete_B_face_estimate_proved": False,
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
        "next_obligation": "On |xi|<=70, remove modes 621 and 622 before summing each rational lattice term in GW6, combine their poles with the exact continuous step-tail currents, and enclose that two-mode-plus-analytic remainder. Treat |xi|>70 through the grouped Abel-theta endpoint difference, not isolated B currents.",
        "proof_boundary": "Exact B window and summed third-order Fresnel remainder only. No leading rational sum, local two-mode enclosure, outside-window estimate, A-fold splice, complete paired residual, complete T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified B Gaussian window and nonlocal Fresnel remainder", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
