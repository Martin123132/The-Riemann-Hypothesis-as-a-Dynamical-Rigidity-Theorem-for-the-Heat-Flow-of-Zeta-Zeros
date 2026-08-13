#!/usr/bin/env python3
"""Guard discrete midpoint parity from invalid interpolation promotion."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_boundary_linear_roster_pair_anchor_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "half_reflection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_kummer_reflection_branch_reduction_gate.json",
    "symmetric_residual": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_symmetric_endpoint_tail_reassembly_gate.json",
}

PRECISION = 100
A = 159_577
B = 5_122_421
L = (B - A) // 2
ALPHA_COUNT = L + 1
WITNESS_MODES = (1, 621, 622, 39_852, 39_894, 39_895)


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


def complex_record(value: acb) -> dict[str, str]:
    return {"real_ball": value.real.str(PRECISION, more=True), "imag_ball": value.imag.str(PRECISION, more=True)}


def symbolic_certificate() -> dict[str, str]:
    u = sp.symbols("u", real=True)
    endpoint_a = sp.symbols("A", real=True)
    alpha = endpoint_a + 2 * u
    canonical = alpha * sp.exp(sp.I * sp.pi * alpha**2 / 8)
    derivative = sp.factor(sp.diff(canonical, u) / sp.exp(sp.I * sp.pi * alpha**2 / 8))
    expected = 2 + sp.I * sp.pi * alpha**2 / 2
    require(sp.simplify(derivative - expected) == 0, "canonical interpolation derivative failed")
    require(sp.simplify(sp.diff(canonical, u, 2)) != 0, "canonical interpolation unexpectedly linear")

    direct_sum = ALPHA_COUNT * (A + B) // 2
    require(direct_sum == sum(range(A, B + 1, 2)), "discrete midpoint roster sum failed")
    return {
        "valid_discrete_character": "At integer u, alpha=A+2u is odd and exp(i*pi*alpha^2/8)=exp(i*pi/8).",
        "valid_direct_sum": "sum_(u=0)^L f_(1/2)(u)=exp(i*pi/8)(L+1)(A+B)/2",
        "canonical_interpolation": "f_(1/2)(u)=(A+2u)exp(i*pi(A+2u)^2/8) for every real 0<=u<=L",
        "canonical_derivative": str(derivative),
        "interpolation_guard": "Odd-square parity applies at integer roster points only. It does not turn the continuous finite-Poisson interpolation into a linear function.",
        "midpoint_mode_formula": "For m!=0, I_m(1/2)=4m(-1)^m{F[(B-4m)/2]-F[(A-4m)/2]}, with F'(q)=exp(i*pi*q^2/2).",
        "endpoint_current_at_midpoint": "The explicit E_B-E_A boundary current vanishes at x=1/2, but the finite Fresnel difference generally does not.",
    }


def fresnel(q: arb, pi: arb, imaginary_unit: acb) -> acb:
    return (1 + imaginary_unit) / 2 * (((-imaginary_unit * pi / 4).exp() * (pi / 2).sqrt() * q).erf())


def midpoint_mode(mode: int, pi: arb, imaginary_unit: acb) -> acb:
    m = arb(mode)
    sign = -1 if abs(mode) % 2 else 1
    q_a = (arb(A) - 4 * m) / 2
    q_b = (arb(B) - 4 * m) / 2
    return 4 * m * sign * (fresnel(q_b, pi, imaginary_unit) - fresnel(q_a, pi, imaginary_unit))


def witness_certificate() -> list[dict[str, Any]]:
    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    imaginary_unit = acb(0, 1)
    rows = []
    for mode in WITNESS_MODES:
        positive = midpoint_mode(mode, pi, imaginary_unit)
        negative = midpoint_mode(-mode, pi, imaginary_unit)
        pair = positive + negative
        require(not pair.real.contains(0) or not pair.imag.contains(0), f"midpoint pair counterexample vanished at mode {mode}")
        rows.append(
            {
                "mode": mode,
                "I_plus_ball": complex_record(positive),
                "I_minus_ball": complex_record(negative),
                "pair_ball": complex_record(pair),
                "pair_absolute_ball": abs(pair).str(PRECISION, more=True),
            }
        )
    require(abs(midpoint_mode(39_894, pi, imaginary_unit) + midpoint_mode(-39_894, pi, imaginary_unit)) > arb("69000"), "fold-edge counterexample scale drift")
    return rows


def render_note(artifact: dict[str, Any]) -> str:
    rows = artifact["witnesses"]
    fold = next(row for row in rows if row["mode"] == 39_894)
    return f"""# Midpoint discrete-parity interpolation guard

Date: 2026-08-13

Status: countermodel gate; not a proof of the half-domain completion

For each integer roster point `u`, `alpha=A+2u` is odd, so the direct
equation-(9) values at `x=1/2` do simplify:

```text
exp(i*pi*alpha^2/8)=exp(i*pi/8),
sum_(u=0)^L f_(1/2)(u)=exp(i*pi/8)(L+1)(A+B)/2.        (MG1)
```

That identity is valid and remains useful for checking the direct source
sum.  It does **not** imply that the canonical finite-Poisson interpolation
is linear.  Between integer points,

```text
f_(1/2)(u)=(A+2u)exp(i*pi(A+2u)^2/8),

f_(1/2)'(u)=exp(i*pi(A+2u)^2/8)
 [2+i*pi(A+2u)^2/2].                                   (MG2)
```

The quadratic chirp remains present for real noninteger `u`.  Consequently,
for nonzero integer `m`, the actual coefficient is

```text
I_m(1/2)=4m(-1)^m
 [F((B-4m)/2)-F((A-4m)/2)],                            (MG3)
```

where `F'(q)=exp(i*pi*q^2/2)`.  The explicit endpoint exponential in the
integration-by-parts formula vanishes because `A` and `B` have the same odd
character, but the Fresnel difference in (MG3) does not.

Actual-source interval witnesses disprove `I_m+I_-m=0`.  At `m=1`,

```text
|I_1+I_-1|={rows[0]['pair_absolute_ball']},
```

and at the fold edge `m=39894`,

```text
|I_39894+I_-39894|={fold['pair_absolute_ball']}>69000. (MG4)
```

The large midpoint value in (MG4) is not itself an error bound; the outer
`x` integral remains oscillatory.  Its role is to block an invalid argument:
discrete odd-square parity cannot be pushed through the continuous Poisson
interpolation coefficient.  Any half-domain completion must retain the
canonical chirp and use its phase, not replace it by a linear trapezoidal
interpolant.

Pi provenance: all `pi` factors come from the equation-(9) Kummer phase and
Fourier character.  No fitted constant is used.

Proof boundary: direct midpoint roster identity and explicit counterexamples
to coefficient-wise pair cancellation only.  No interior completion bound or
nonlinear B crossing theorem is proved.  No theorem establishing `T_upper`,
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion follows here.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    artifact = {
        "kind": STEM,
        "status": "discrete_midpoint_parity_valid_but_linear_interpolation_and_pair_cancellation_promotion_rejected",
        "passed": True,
        "scope": {"source_endpoints": [A, B], "L": L, "alpha_count": ALPHA_COUNT, "x": "1/2"},
        "symbolic_certificate": symbolic_certificate(),
        "witnesses": witness_certificate(),
        "decision": {
            "direct_discrete_midpoint_roster_simplifies": True,
            "canonical_continuous_interpolation_is_linear": False,
            "all_nonzero_symmetric_Poisson_pairs_cancel_at_midpoint": False,
            "one_sided_nonzero_midpoint_coefficients_are_purely_imaginary": False,
            "discrete_odd_square_character_may_be_applied_between_roster_points": False,
            "prior_linear_interpolation_anchor_rejected": True,
            "interior_half_domain_completion_bounded": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {name: {"path": relative(path), "sha256": file_hash(path)} for name, path in DEPENDENCIES.items()},
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {"elapsed_seconds": round(time.time() - started, 3), "workers": 1, "flint_threads": 1, "process_priority": priority},
        "next_obligation": "Retain the canonical midpoint chirp in the paired half-domain completion; derive a phase-adapted contour or paired Morse representation rather than using discrete parity inside continuous Fourier coefficients.",
        "proof_boundary": "Direct midpoint roster identity and explicit interpolation/pair-cancellation guard only. No interior completion bound, B crossing theorem, A splice, complete T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified midpoint interpolation guard: discrete parity does not cancel Poisson pairs", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
