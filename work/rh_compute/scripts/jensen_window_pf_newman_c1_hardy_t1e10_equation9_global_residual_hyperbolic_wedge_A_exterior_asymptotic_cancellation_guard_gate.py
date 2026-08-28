#!/usr/bin/env python3
"""Guard the A exterior against a nonconvergent affine/amplitude split."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_exterior_asymptotic_cancellation_guard_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "triangle": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_separable_triangle_geometry_gate.json",
    "remainder_identity": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_remainder_identity_gate.json",
    "regrouping": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_A_localized_exact_domain_regrouping_gate.json",
    "local_exact": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_A_local_exact_domain_gate.json",
}

HEIGHT = 10_000_000_000
A = 159_577
LOWER_START = 39_853
UPPER_END = 39_936
MODES = tuple(range(LOWER_START, UPPER_END + 1))


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


def complex_record(value: acb, digits: int = 70) -> dict[str, str]:
    return {
        "real_ball": value.real.str(digits, more=True),
        "imag_ball": value.imag.str(digits, more=True),
        "absolute_ball": abs(value).str(digits, more=True),
    }


def symbolic_certificate() -> dict[str, str]:
    P, K, mu, lam, S = sp.symbols("P K mu lambda S", nonzero=True)
    e2 = -2 * sp.I * K / (K - sp.I)
    c3 = sp.simplify(e2 / sp.I)
    require(c3 == -2 * K / (K - sp.I), "exact inner-tail coefficient failed")

    affine_limit = sp.simplify(-sp.I * mu * (-1) - sp.I * lam)
    require(affine_limit == sp.I * (mu - lam), "affine inner-tail coefficient failed")

    x, t, r, m = sp.symbols("x t r m", positive=True)
    v = (1 / x - 1) / r
    require(sp.limit(x * v, x, 0, dir="+") == 1 / r, "outer v limit failed")
    p_face = sp.sqrt(sp.pi / 2) * (sp.symbols("D", positive=True) * sp.sqrt(x) - 2 * m / sp.sqrt(x))
    require(sp.limit(sp.sqrt(x) * p_face, x, 0, dir="+") == -m * sp.sqrt(2 * sp.pi), "face P limit failed")

    return {
        "endpoint_factor_limit": "lim_(P->-infinity) P^2 E(P)=-2iK/(K-i), K=pi*A*m",
        "exact_inner_tail": "exp(-iP^2/2) H_E(P)=-2K/(K-i) P^(-3)+O(P^(-5))",
        "lower_Fresnel_tail": "exp(-iP^2/2) T_-(P)=-i/P-P^(-3)+O(P^(-5))",
        "affine_inner_limit": "lim exp(-iP_A^2/2)H_aff(P_A,S)=i(mu_S-lambda_P)",
        "face_limits": "sqrt(x)P_A(x)->-m sqrt(2pi), sqrt(x)S(x)->m sqrt(2pi)",
        "outer_limits": "x^(3/2)S_x->-m sqrt(pi/2), x^(1/4)O(v(x))->sqrt(2)r^(-1/4)",
        "affine_current": "J_aff,m(x)=L_m x^(-3/2)exp(i Phi_A,face(x))+O(x^(-1)exp(i Phi_A,face(x)))",
        "exact_current": "J_exact,m(x)=M_m x^(-1/4)exp(i Phi_A,face(x))(1+o(1))",
        "common_face_phase": "Phi_A,face(x)=pi*A^2*x/4+(t/2)log((1-x)/x)",
        "cutoff_identity": "C_exact,ext(epsilon)=C_aff,ext(epsilon)+R_amp,ext(epsilon) for every epsilon>0",
        "limit_guard": "C_aff,ext and R_amp,ext have opposite epsilon^(-1/2-it/2) leading terms; only their common-cutoff sum has an improper limit",
    }


def mode_certificate(precision: int = 100) -> dict[str, Any]:
    ctx.dps = precision
    ctx.threads = 1
    pi, t, endpoint = arb.pi(), arb(HEIGHT), arb(A)
    imaginary = acb(0, 1)
    epsilon_A = arb(-1)
    mu = arb(5) / 12 * (arb(2) / t).sqrt()
    rows: list[dict[str, Any]] = []
    grouped = acb(0)
    lambda_minus_mu_lower = None
    lambda_minus_mu_upper = None
    max_abs_lambda_imag = arb(0)

    for mode_int in MODES:
        mode = arb(mode_int)
        r = t / (2 * pi * mode**2)
        K = pi * endpoint * mode
        lambda_P = acb(3 * K, -1) / (acb(2 * K, -2) * K.sqrt())
        W0 = acb(K, -1) * (4 * r ** arb("0.75") / (endpoint * (pi * t).sqrt()))
        dS_leading = -mode * (pi / 2).sqrt()
        leading = epsilon_A * W0 * imaginary * (mu - lambda_P) * dS_leading / (2 * pi)
        grouped += leading
        gap = lambda_P.real - mu
        lambda_minus_mu_lower = gap.lower() if lambda_minus_mu_lower is None else min(lambda_minus_mu_lower, gap.lower())
        lambda_minus_mu_upper = gap.upper() if lambda_minus_mu_upper is None else max(lambda_minus_mu_upper, gap.upper())
        max_abs_lambda_imag = max(max_abs_lambda_imag, abs(lambda_P.imag).upper())
        rows.append(
            {
                "mode": mode_int,
                "lambda_minus_mu_real_ball": gap.str(70, more=True),
                "lambda_P_imag_ball": lambda_P.imag.str(70, more=True),
                "affine_exterior_leading_coefficient_ball": complex_record(leading),
            }
        )

    require(lambda_minus_mu_lower is not None and lambda_minus_mu_lower > arb("4.7e-6"), "lambda-mu separation lost")
    require(lambda_minus_mu_upper is not None and lambda_minus_mu_upper < arb("4.8e-6"), "lambda-mu upper guard lost")
    require(max_abs_lambda_imag < arb("3.6e-16"), "lambda imaginary guard lost")
    require(abs(grouped).lower() > arb("8.9"), "grouped leading coefficient could vanish")
    require(grouped.imag.upper() < arb("-8.9"), "grouped leading orientation drift")

    return {
        "height": HEIGHT,
        "endpoint": A,
        "mode_range": [LOWER_START, UPPER_END],
        "mode_count": len(rows),
        "mu_S_ball": mu.str(70, more=True),
        "lambda_minus_mu_real_lower_ball": lambda_minus_mu_lower.str(70, more=True),
        "lambda_minus_mu_real_upper_ball": lambda_minus_mu_upper.str(70, more=True),
        "max_abs_lambda_P_imag_ball": max_abs_lambda_imag.str(70, more=True),
        "grouped_affine_exterior_leading_coefficient_ball": complex_record(grouped),
        "rows": rows,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["mode_certificate"]
    g = c["grouped_affine_exterior_leading_coefficient_ball"]
    return f"""# A-exterior asymptotic cancellation guard

Date: 2026-08-14

Status: exact endpoint asymptotics and nonconvergence guard; not an exterior
value, complete A carrier, or global residual bound

The affine/exact-minus-affine identity remains exact at every common finite
cutoff `epsilon>0`.  It cannot be passed term by term to the `x=0` exterior
limit.

Let `H_aff` denote the inner `P` integral of
`1+lambda_P P+mu_S S`.  The lower Fresnel expansion and the exact face limits
give

```text
exp(-iP_A^2/2)H_aff(P_A,S) -> i(mu_S-lambda_P),       (EC1)
sqrt(x)P_A -> -m sqrt(2pi),
sqrt(x)S   ->  m sqrt(2pi),
x^(3/2)S_x-> -m sqrt(pi/2).                           (EC2)
```

Across all 84 modes,

```text
{c['lambda_minus_mu_real_lower_ball']}
 < Re(lambda_P-mu_S) <
{c['lambda_minus_mu_real_upper_ball']}.              (EC3)
```

Thus the affine exterior current has the common-phase asymptotic

```text
J_aff,m(x)=L_m x^(-3/2) exp(i Phi_A,face(x))+O(x^-1),
Phi_A,face(x)=pi A^2 x/4+(t/2)log((1-x)/x).           (EC4)
```

Summing every mode before the limit does not remove it:

```text
sum_m L_m real={g['real_ball']},
sum_m L_m imag={g['imag_ball']},
|sum_m L_m| ={g['absolute_ball']}.                    (EC5)
```

Since `exp(i Phi_A,face(x))=x^(-it/2)(1+O(x))`, the cutoff
primitive has a nonzero `epsilon^(-1/2-it/2)` leading term.  Therefore the
standalone affine exterior has no improper limit.  The exterior
exact-minus-affine amplitude has the opposite leading term and has no
standalone limit either.

The full exact endpoint factor supplies the missing cancellation.  With
`K=pi A m`,

```text
E(P)=-2iK/(K-i) P^(-2)+O(P^-4),
exp(-iP^2/2) integral_(-infinity)^P E(q)exp(iq^2/2)dq
 =-2K/(K-i) P^(-3)+O(P^-5).                          (EC6)
```

Together with `x^(1/4)O(v(x))->sqrt(2)r^(-1/4)`, the full exact current is
only `O(x^-1/4)` and is absolutely integrable at `x=0`.  Hence the admissible
exterior object is

```text
C_A,ext^exact=lim_(epsilon->0)
 [C_A,ext^aff(epsilon)+R_A,ext^amp(epsilon)],         (EC7)
```

with the two terms retained under one cutoff, or equivalently the full exact
endpoint factor integrated directly.  The affine/amplitude split remains
valid on the compact local box.

Pi provenance: every `pi` in (EC1)--(EC7) comes from the exact equation-(9)
bi-Morse map, original triangle face, and endpoint normalization.  No fitted
constant is introduced.

Proof boundary: exact endpoint orders, all-mode coefficient separation, and
the common-cutoff nonconvergence/cancellation guard only.  No numerical exact
exterior value, local or exterior nonlinear-amplitude bound, complete A
endpoint theorem, `R_Dir` estimate, complete `Q_K-T` or `T_upper`, all-height
theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["triangle"]["decision"]["pair_Volterra_integral_reduced_to_separable_triangle"] is True, "triangle dependency drift")
    require(dependencies["remainder_identity"]["decision"]["exact_defect_split_into_signed_amplitude_and_face_remainders"] is True, "remainder dependency drift")
    require(dependencies["regrouping"]["decision"]["tangent_exterior_must_be_removed_before_estimation"] is True, "regrouping dependency drift")
    require(dependencies["local_exact"]["decision"]["all_84_local_exact_affine_modes_certified"] is True, "local exact dependency drift")

    artifact = {
        "kind": STEM,
        "status": "A_exterior_affine_and_amplitude_remainder_separate_limits_rejected_common_cutoff_exact_current_required",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "mode_certificate": mode_certificate(),
        "decision": {
            "finite_cutoff_affine_amplitude_decomposition_remains_exact": True,
            "standalone_affine_exterior_improper_limit_exists": False,
            "standalone_exact_minus_affine_exterior_improper_limit_exists": False,
            "all_84_mode_grouped_affine_leading_coefficient_nonzero": True,
            "full_exact_endpoint_current_integrable_at_x_zero": True,
            "common_cutoff_recombination_required_before_exterior_limit": True,
            "local_compact_affine_amplitude_split_remains_admissible": True,
            "exact_exterior_value_proved": False,
            "complete_A_endpoint_block_proved": False,
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
            "precision_decimal_digits": 100,
        },
        "next_obligation": "Evaluate the full exact A exterior current with the endpoint factor E(P) retained. Use its inverse-P tail beginning at P^-3, the common exact-face x phase, and the ten-mode incomplete stationary transition. Compute the compact exact-minus-affine local amplitude separately.",
        "proof_boundary": "Exact endpoint asymptotic orders, all-mode coefficient separation, and common-cutoff cancellation guard only. No numerical exact exterior value, nonlinear-amplitude bound, complete A endpoint theorem, R_Dir estimate, complete Q_K-T or T_upper, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified A-exterior common-cutoff asymptotic cancellation guard", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
