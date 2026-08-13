#!/usr/bin/env python3
"""Certify the cubic resonance normal form on the selected Airy lattice."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_detuning_airy_lattice_cubic_resonance_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "airy_green_kernel": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_beta4_detuning_airy_green_kernel_gate.json",
    "event_pairing": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_event_ordered_398_pairing_multiplier_gate.json",
}

C = 159_577
Q = (C - 1) // 4
J_MIN = -198
J_MAX = 200
PRECISION = 100
DLMF_URL = "https://dlmf.nist.gov/9.8"
BINOMIAL = (sp.Rational(3, 2), sp.Rational(3, 8), sp.Rational(-1, 16), sp.Rational(3, 128), sp.Rational(-3, 256))


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


def symbolic_certificate() -> dict[str, str]:
    c, u = sp.symbols("C u", positive=True, real=True)
    j = sp.symbols("j", real=True)
    xj = (8 * j - 2) / c + u
    x0 = -2 / c + u
    truncated = sum(BINOMIAL[k - 1] * (xj**k - x0**k) for k in range(1, 6))
    resonant_integer = sp.pi * j * (c - 1 + 2 * j)
    p5 = sp.factor(sp.pi * c**2 * truncated / 12 - resonant_integer)
    at_center = sp.factor(p5.subs(u, 0))
    expected_center = -sp.pi * j * (16 * j**2 - 12 * j + 3) / (6 * c)
    expected_center += sp.pi * ((8 * j - 2) ** 4 - 16) / (512 * c**2)
    expected_center -= sp.pi * ((8 * j - 2) ** 5 + 32) / (1024 * c**3)
    require(sp.expand(at_center - expected_center) == 0, "center quintic phase polynomial failed")

    q = sp.symbols("q", integer=True)
    require(sp.expand((4 * q + 1 - 1 + 2 * j) * j - 2 * j * (2 * q + j)) == 0, "resonant integer parity failed")
    cubic = -sp.pi * j * (16 * j**2 - 12 * j + 3) / (6 * c)
    return {
        "scaled_height": "u=lambda/beta^2; 0<=u<=1/(2C^2)",
        "relative_coordinate": "r_(q+j)=beta^2[1+(8j-2)/C+u]",
        "exact_resonant_integer": "pi*j*(C-1+2j)=2pi*j*(2q+j)",
        "quintic_normal_form": str(p5),
        "center_cubic_term": str(sp.factor(cubic)),
        "center_quartic_term": "pi*((8j-2)^4-16)/(512C^2)",
        "center_quintic_term": "-pi*((8j-2)^5+32)/(1024C^3)",
        "interpretation": "The constant, linear, and quadratic lattice phase is an exact multiple of 2pi; the first nontrivial relative phase is cubic in j.",
    }


def arb_binomial() -> tuple[arb, ...]:
    return (arb(3) / 2, arb(3) / 8, -arb(1) / 16, arb(3) / 128, -arb(3) / 256)


def quintic_phase(j: int, u: arb, c: arb, pi: arb) -> arb:
    xj = (arb(8) * j - 2) / c + u
    x0 = -arb(2) / c + u
    truncated = arb(0)
    for degree, coefficient in enumerate(arb_binomial(), start=1):
        truncated += coefficient * (xj**degree - x0**degree)
    return pi * c**2 * truncated / 12 - pi * j * (c - 1 + 2 * j)


def exact_leading_phase(j: int, u: arb, c: arb, pi: arb) -> arb:
    xj = (arb(8) * j - 2) / c + u
    x0 = -arb(2) / c + u
    return pi * c**2 / 12 * ((1 + xj) ** (arb(3) / 2) - (1 + x0) ** (arb(3) / 2))


def interval_certificate() -> dict[str, Any]:
    pi = arb.pi()
    c = arb(C)
    beta = (pi * c**2 / 8) ** (arb(1) / 3)
    u_max = 1 / (2 * c**2)
    u = arb(u_max / 2, u_max / 2)
    x_max = (arb(8) * J_MAX - 2) / c + u_max
    x_min = (arb(8) * J_MIN - 2) / c
    z_max = max(abs(x_min), abs(x_max))
    x0_max = arb(2) / c
    taylor_remainder = (
        pi * c**2 / 12
        * arb(7) / 1024
        * (1 - z_max) ** (-arb(9) / 2)
        * (z_max**6 + x0_max**6)
    )
    r_min = beta**2 * (1 + x_min)
    airy_phase_remainder = arb(5) / (24 * r_min ** (arb(3) / 2))
    total_remainder = taylor_remainder + airy_phase_remainder
    require(total_remainder < arb("4.82e-5"), "quintic Airy phase remainder exceeds 4.82e-5")

    # For j=1 the leading-phase residual is strictly increasing in u:
    # its derivative is pi*C^2/8 times
    # sqrt(1+6/C+u)-sqrt(1-2/C+u), which is positive.  Evaluating the
    # two endpoints separately avoids the severe dependency inflation caused
    # by putting the same interval u into both square roots.
    j_one_low = exact_leading_phase(1, arb(0), c, pi) - 2 * pi * (2 * Q + 1)
    j_one_high = exact_leading_phase(1, u_max, c, pi) - 2 * pi * (2 * Q + 1)
    j_one_mid = (j_one_low.lower() + j_one_high.upper()) / 2
    j_one_radius = (j_one_high.upper() - j_one_low.lower()) / 2 + airy_phase_remainder
    j_one_exact = arb(j_one_mid, j_one_radius)
    require(j_one_exact.upper() < arb("-1.7e-5"), "central adjacent phase lost its negative separation")
    require(j_one_exact.lower() > arb("-2.4e-5"), "central adjacent phase escaped expected interval")

    rows = []
    maximum_quintic_error = arb(0)
    for index in (J_MIN, -150, -100, -50, -2, -1, 1, 2, 50, 100, 150, J_MAX):
        point_u = arb(0)
        exact = exact_leading_phase(index, point_u, c, pi) - pi * index * (c - 1 + 2 * index)
        normal = quintic_phase(index, point_u, c, pi)
        error = abs(exact - normal)
        maximum_quintic_error = max(maximum_quintic_error, error.upper())
        rows.append({
            "j": index,
            "leading_phase_mod_2pi_ball": exact.str(PRECISION, more=True),
            "quintic_normal_form_ball": normal.str(PRECISION, more=True),
            "difference_absolute_ball": error.str(PRECISION, more=True),
        })
    return {
        "j_interval": [J_MIN, J_MAX],
        "mode_interval": [Q + J_MIN, Q + J_MAX],
        "u_interval": "0<=u<=1/(2C^2)",
        "u_max_ball": u_max.str(PRECISION, more=True),
        "maximum_abs_scaled_argument_ball": z_max.str(PRECISION, more=True),
        "uniform_binomial_degree5_remainder_ball": taylor_remainder.str(PRECISION, more=True),
        "uniform_exact_Airy_phase_correction_ball": airy_phase_remainder.str(PRECISION, more=True),
        "uniform_exact_phase_normal_form_error_ball": total_remainder.str(PRECISION, more=True),
        "central_adjacent_exact_phase_residual_ball": j_one_exact.str(PRECISION, more=True),
        "sample_rows": rows,
        "sampled_maximum_quintic_error_ball": maximum_quintic_error.str(PRECISION, more=True),
        "pairwise_error_rule": "Subtracting two q-referenced phase normal forms incurs at most twice the uniform error.",
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    return f"""# Cubic resonance normal form on the Airy event lattice

Date: 2026-08-13

Status: exact resonance cancellation and uniform quintic Airy-phase normal
form proved; no cubic exponential-sum bound

Let `vartheta(r)` be the exact Airy phase defined by
`Ai(-r)=M(-r)sin(vartheta)` and `Bi(-r)=M(-r)cos(vartheta)`.  For
`m=q+j`, put

```text
u=lambda/beta^2,       0<=u<=1/(2C^2),
r_j=beta^2[1+(8j-2)/C+u].                            (CR1)
```

The upper bound on `u` is exact: `lambda<=pi/(16beta)` and
`beta^3=pi C^2/8` imply `lambda/beta^2<=1/(2C^2)`.

For the leading Airy phase `xi(r)=2r^(3/2)/3`, the constant, linear, and
quadratic binomial terms satisfy

```text
pi*j*(C-1+2j)=2pi*j*(2q+j).                          (CR2)
```

They are therefore an exact integer multiple of `2pi`, not a small
numerical accident.  If

```text
x_j=(8j-2)/C+u,       x_0=-2/C+u,

P_5(j,u)=pi C^2/12 sum_(k=1)^5 binom(3/2,k)
                         (x_j^k-x_0^k)
          -pi*j*(C-1+2j),                            (CR3)
```

then at the selector centre

```text
P_5(j,0)= -pi*j*(16j^2-12j+3)/(6C)
 +pi[(8j-2)^4-16]/(512C^2)
 -pi[(8j-2)^5+32]/(1024C^3).                         (CR4)
```

The first nontrivial term is cubic in `j`.  This is the same fold mechanism
seen in the physical Airy chart, now recovered inside the detuning Green
lattice.

For every `-198<=j<=200` and the full top-corridor height interval,

```text
|x_j|<={c['maximum_abs_scaled_argument_ball']},
|binomial remainder through degree 5|
 <{c['uniform_binomial_degree5_remainder_ball']},
|exact Airy phase correction|
 <{c['uniform_exact_Airy_phase_correction_ball']}.   (CR5)
```

DLMF 9.8.22 and its signed first-neglected-term remainder supply the last
line.  Combining both errors proves

```text
|vartheta(r_j)-vartheta(r_0)
 -2pi*j*(2q+j)-P_5(j,u)|
 <{c['uniform_exact_phase_normal_form_error_ball']}<4.82e-5. (CR6)
```

For the nearest adjacent lattice point, the exact residual throughout the
top corridor is

```text
{c['central_adjacent_exact_phase_residual_ball']}.    (CR7)
```

Thus a generic Dirichlet denominator is necessarily poor near `j=0`: the
phase increment is almost resonant because two complete asymptotic orders
cancel exactly.  The correct next object is the finite cubic-quintic sum
generated by (CR3), with the event-pair amplitudes and completed endpoint
sources still attached.

Pi provenance: (CR2) comes from the exact selector identity
`beta^3=pi C^2/8`; the Airy phase and its error are sourced from
`{DLMF_URL}`.  No fitted frequency or unrelated circle is used.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No cubic exponential-sum estimate, grouped Green-operator bound, completed
endpoint cancellation, finite-integral 398-mode splice, complete `Q_K-T` or
`T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    require(dependencies["airy_green_kernel"].get("passed") is True, "Airy Green-kernel dependency failed")
    require(dependencies["event_pairing"].get("passed") is True, "event-pairing dependency failed")

    ctx.dps = PRECISION
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "Airy_event_lattice_constant_linear_quadratic_resonance_removed_and_quintic_phase_certified",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(),
        "external_theorem": {
            "source": "NIST Digital Library of Mathematical Functions, Section 9.8",
            "url": DLMF_URL,
            "used_statements": [
                "Airy phase definitions 9.8.1--9.8.4",
                "Airy phase expansion 9.8.22",
                "signed first-neglected-term remainder statement following 9.8.23",
            ],
        },
        "decision": {
            "constant_linear_quadratic_lattice_phase_is_exact_2pi_integer": True,
            "first_nontrivial_relative_phase_is_cubic": True,
            "uniform_quintic_exact_Airy_phase_normal_form_proved": True,
            "uniform_phase_normal_form_error_below_4_82e_minus_5": True,
            "generic_adjacent_Dirichlet_bound_rejected_near_event_zero": True,
            "cubic_exponential_sum_bound_proved": False,
            "grouped_Green_operator_bound_proved": False,
            "grouped_398_mode_splice_proved": False,
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
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
        "next_action": (
            "Derive a finite van der Corput or Poisson bound for the cubic-quintic phase P5(j,u), split only at its "
            "certified curvature transitions, and retain the coherent endpoint-source and forcing amplitudes. Reject any "
            "bound that replaces the central near-resonant block by a first-difference Dirichlet denominator."
        ),
        "proof_boundary": (
            "Exact resonance cancellation and a uniform quintic normal form for the exact Airy phase only. No cubic "
            "exponential-sum estimate, grouped Green-operator bound, completed endpoint cancellation, finite-integral splice, "
            "complete Q_K-T or T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("certified cubic resonance normal form on the Airy event lattice", flush=True)


if __name__ == "__main__":
    main()
