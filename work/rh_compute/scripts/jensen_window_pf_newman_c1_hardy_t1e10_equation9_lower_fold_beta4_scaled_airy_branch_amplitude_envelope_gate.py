#!/usr/bin/env python3
"""Certify exact scaled Airy-branch amplitude envelopes on every fold cell."""

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

from flint import arb, acb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_scaled_airy_branch_amplitude_envelope_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "Airy_reduction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_canonical_airy_derivative_reduction_gate.json",
    "Hankel_split": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_canonical_exact_hankel_branch_factorization_gate.json",
    "event_selection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_hankel_carrier_event_selection_gate.json",
}

C = 159_577
Y = 64
MAX_ODD = 797
PRECISION = 105
DLMF_URL = "https://dlmf.nist.gov/9.8"


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


def airy_envelopes(x: arb) -> tuple[arb, arb, arb]:
    i = acb(0, 1)
    ai, ai_prime, bi, bi_prime = acb(-x).airy()
    w = ai - i * bi
    w_x = -ai_prime + i * bi_prime
    rho = w_x - i * x.sqrt() * w
    return abs(w), abs(w_x), abs(rho)


def direct_envelope_atlas(x_min: arb, x_max: arb) -> tuple[arb, arb, arb, dict[str, int]]:
    w_max = arb(0)
    low_wx_max = arb(0)
    rho_max = arb(0)

    def unit_panel(index: int) -> arb:
        return arb(arb(2 * index + 1) / 2, arb(1) / 2)

    # Near zero, equal log-X panels resolve the finite branch limits.
    low_panels = 8_000
    log_min = x_min.log()
    for index in range(low_panels):
        u = log_min * (1 - unit_panel(index) / low_panels)
        x = u.exp()
        w_value, wx_value, rho_value = airy_envelopes(x)
        if w_value.upper() > w_max.upper():
            w_max = w_value
        if wx_value.upper() > low_wx_max.upper():
            low_wx_max = wx_value
        if rho_value.upper() > rho_max.upper():
            rho_max = rho_value

    # Beyond X=1, equal cubic-phase panels keep every oscillation resolved.
    phase_panels = 120_000
    xi_min = arb(2) / 3
    xi_max = arb(2) * x_max ** (arb(3) / 2) / 3
    for index in range(phase_panels):
        xi = xi_min + (xi_max - xi_min) * unit_panel(index) / phase_panels
        x = (arb(3) * xi / 2) ** (arb(2) / 3)
        w_value, _, rho_value = airy_envelopes(x)
        if w_value.upper() > w_max.upper():
            w_max = w_value
        if rho_value.upper() > rho_max.upper():
            rho_max = rho_value

    return w_max, low_wx_max, rho_max, {
        "log_X_panels": low_panels,
        "uniform_cubic_phase_panels": phase_panels,
        "total_panels": low_panels + phase_panels,
    }


def certificate() -> dict[str, str]:
    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    c = arb(C)
    beta = (pi * c**2 / 8) ** (arb(1) / 3)
    eps = beta**-2
    x_min = pi / (16 * beta)
    d_max = arb(MAX_ODD) * beta / c
    x_max = d_max**2 + x_min + Y
    lambda_abs = arb(116)
    y = arb(Y)

    # Direct triangle bounds for the four rational beta^-4 weights.
    u1 = (13 * lambda_abs + 3 * y) / 60
    v1 = (8 * lambda_abs**2 + 4 * lambda_abs * y + 3 * y**2) / 60
    u2 = (
        448 * lambda_abs**5 + 280 * lambda_abs**2 * y**3 + 4565 * lambda_abs**2
        + 105 * lambda_abs * y**4 + 30 * lambda_abs * y + 63 * y**5 + 405 * y**2
    ) / 50400
    v2 = (40 * lambda_abs**3 + 20 * lambda_abs**2 * y + lambda_abs * y**2 + 9 * y**3 + 27) / 1680
    uy1 = arb(3) / 60
    vy1 = (4 * lambda_abs + 6 * y) / 60
    uy2 = (840 * lambda_abs**2 * y**2 + 420 * lambda_abs * y**3 + 30 * lambda_abs + 315 * y**4 + 810 * y) / 50400
    vy2 = (20 * lambda_abs**2 + 2 * lambda_abs * y + 27 * y**2) / 1680
    u_defect = eps * u1 + eps**2 * u2
    v_bound = eps * v1 + eps**2 * v2
    uy_bound = eps * uy1 + eps**2 * uy2
    vy_bound = eps * vy1 + eps**2 * vy2
    require(u_defect < arb("1.6e-5"), "U defect envelope failed")
    require(v_bound < arb("5.38e-4"), "V envelope failed")
    require(uy_bound < arb("7.1e-8"), "U_y envelope failed")
    require(vy_bound < arb("3.05e-6"), "V_y envelope failed")

    m_bound, n_low_bound, rho_bound, atlas = direct_envelope_atlas(x_min, x_max)
    # DLMF 9.8.21 with one retained correction and its signed next-term
    # remainder gives N(-X)^2 <= sqrt(X)/pi*(1+7/(32 X^3)) for X>=1.
    n_high_squared = x_max.sqrt() / pi * (1 + arb(7) / (32 * x_max**3))
    n_at_one_squared = (1 + arb(7) / 32) / pi
    require(n_high_squared > n_at_one_squared, "high-X derivative envelope endpoint order failed")
    n_bound = n_high_squared.sqrt()
    require(n_bound > n_low_bound, "low-X derivative atlas exceeds analytic high-X envelope")
    require(rho_bound < arb("0.512"), "exact phase-residual endpoint envelope failed")
    require(m_bound < arb("0.711"), "Airy modulus envelope failed")
    require(n_bound < arb("2.07"), "Airy derivative modulus envelope failed")

    s_bound = ((1 + u_defect) * m_bound + v_bound * n_bound) / 2
    sy_bound = (
        uy_bound * m_bound
        + vy_bound * n_bound
        + (1 + u_defect + x_max.sqrt() * v_bound) * rho_bound
    ) / 2
    require(s_bound < arb("0.356"), "scaled branch amplitude bound failed")
    require(sy_bound < arb("0.258"), "scaled branch derivative bound failed")

    return {
        "beta_ball": beta.str(PRECISION, more=True),
        "epsilon_ball": eps.str(PRECISION, more=True),
        "X_min_ball": x_min.str(PRECISION, more=True),
        "X_max_ball": x_max.str(PRECISION, more=True),
        "lambda_absolute_envelope": "116",
        "U_minus_1_absolute_bound_ball": u_defect.str(PRECISION, more=True),
        "V_absolute_bound_ball": v_bound.str(PRECISION, more=True),
        "U_y_absolute_bound_ball": uy_bound.str(PRECISION, more=True),
        "V_y_absolute_bound_ball": vy_bound.str(PRECISION, more=True),
        "W_modulus_bound_ball": m_bound.str(PRECISION, more=True),
        "W_X_modulus_bound_ball": n_bound.str(PRECISION, more=True),
        "W_X_low_X_direct_atlas_bound_ball": n_low_bound.str(PRECISION, more=True),
        "W_X_high_X_squared_asymptotic_bound_ball": n_high_squared.str(PRECISION, more=True),
        "extracted_phase_residual_bound_ball": rho_bound.str(PRECISION, more=True),
        "scaled_branch_amplitude_bound_ball": s_bound.str(PRECISION, more=True),
        "scaled_branch_y_derivative_bound_ball": sy_bound.str(PRECISION, more=True),
        "direct_interval_atlas": atlas,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Beta^-4 scaled Airy-branch amplitude envelope

Date: 2026-08-13
Status: rigorous saved-height envelope; not a proof of the ordinary splice

For branch sign `sigma in {{+1,-1}}`, set

```text
W_sigma(X)=Ai(-X)-i*sigma*Bi(-X),
W_sigma,X=d_X W_sigma,
xi=2X^(3/2)/3,

S_sigma(X,y)=e^(-i*sigma[xi-pi/4])
  [U W_sigma+V W_sigma,X]/2.                           (AE1)
```

By the exact Airy--Hankel connection in Section 11.378, `S_sigma` is the
complete beta-minus-four branch amplitude after extracting only the cubic
carrier `exp(i*sigma[xi-pi/4])`.  There is no asymptotic replacement in
(AE1).

On all 399 event cells and `0<=y<=64`, Section 11.379 gives

```text
{c['X_min_ball']} <=X<={c['X_max_ball']},
|lambda|<116.                                           (AE2)
```

Direct rational bounds for (11.377.6) give

```text
|U-1|<1.6e-5,       |V|<5.38e-4,
|U_y|<7.1e-8,       |V_y|<3.05e-6.                    (AE3)
```

For negative Airy argument, the DLMF modulus definitions identify

```text
|W_sigma|=M(-X),       |W_sigma,X|=N(-X).              (AE4)
```

A deterministic Arb atlas uses 8,000 equal `log X` panels below `X=1` and
120,000 equal cubic-phase panels above it.  For the derivative modulus on
`X>=1`, DLMF 9.8.21 and its signed next-term remainder give

```text
N(-X)^2<=sqrt(X)/pi [1+7/(32X^3)].                    (AE5)
```

The right side has no interior maximum and its upper endpoint dominates at
the saved `X_max`.  Combining this with the direct low-`X` atlas gives,
without a numerical monotonicity assumption,

```text
|W_sigma|<{c['W_modulus_bound_ball']}<0.711,
|W_sigma,X|<{c['W_X_modulus_bound_ball']}<2.07.        (AE6)
```

The derivative after extracting the cubic carrier uses

```text
rho_sigma=W_sigma,X-i*sigma*sqrt(X)W_sigma.            (AE7)
```

The same interval atlas encloses the derivative after extracting the cubic
carrier,

```text
|rho_sigma|<{c['extracted_phase_residual_bound_ball']}<0.512. (AE8)
```

Since `W_sigma,XX=-X W_sigma`, differentiation of (AE1) is exact:

```text
S_sigma,y=e^(-i*sigma[xi-pi/4])/2
 {{U_y W_sigma+V_y W_sigma,X
   +(U-i*sigma*sqrt(X)V)rho_sigma}}.                   (AE9)
```

Combining (AE3), (AE5), and (AE7) proves uniformly

```text
|S_sigma|<{c['scaled_branch_amplitude_bound_ball']}<0.356,
|S_sigma,y|<{c['scaled_branch_y_derivative_bound_ball']}<0.258. (AE10)
```

Pi provenance: the extracted `pi/4` is the forced Hankel phase, and the
`pi` inside `beta` and the event-cell endpoints comes from the Kummer/Fourier
normalization.  The Airy modulus definitions are recorded from `{DLMF_URL}`;
the finite inequalities in (AE5)--(AE10) combine its signed asymptotic
remainder with direct Arb enclosures.

This gate controls the exact canonical branch amplitude on the compact fold
strip.  It does not yet compare that amplitude with the exact finite-height
logistic/Gamma carrier, integrate the opposite branch over every owned
corridor, close `Q_K-T` or `T_upper`, prove a height-uniform theorem,
`Lambda<=0`, PF-infinity, RH, or a prize-level result.
"""


def main() -> int:
    started = time.perf_counter()
    priority = set_low_priority()
    require(CHECKER.is_file(), "missing independent checker")
    dependencies = {}
    for name, path in DEPENDENCIES.items():
        require(path.is_file(), f"missing dependency: {name}")
        dependencies[name] = json.loads(path.read_text(encoding="utf-8"))
    require(dependencies["Airy_reduction"]["decision"]["all_derivatives_reduced_to_Ai_and_first_derivative"] is True, "Airy reduction drift")
    require(dependencies["Hankel_split"]["decision"]["exact_Hankel_branch_factorization_proved"] is True, "Hankel split drift")
    require(dependencies["event_selection"]["decision"]["all_399_event_cells_strictly_positive_X"] is True, "event selection drift")

    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "exact_scaled_Airy_branch_amplitude_and_first_y_derivative_uniformly_bounded_on_all_399_fold_cells",
        "passed": True,
        "certificate": certificate(),
        "external_theorem": {
            "source": "NIST Digital Library of Mathematical Functions, Section 9.8",
            "url": DLMF_URL,
            "used_statements": ["Airy modulus definitions in 9.8(i)", "derivative-modulus expansion 9.8.21", "signed next-term remainder statement following 9.8.23"],
        },
        "decision": {
            "exact_scaled_branch_amplitude_bound_proved": True,
            "exact_scaled_branch_first_y_derivative_bound_proved": True,
            "all_399_event_cells_covered": True,
            "ordinary_logistic_Gamma_amplitude_identification_proved": False,
            "opposite_branch_corridor_integral_bound_proved": False,
            "complete_Q_K_minus_T_bound_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "next_action": "Use the certified |S| and |S_y| bounds with the opposite extracted-carrier phase gap to close its compact fold-strip integral, then derive the selected-branch amplitude bridge to the exact logistic/Gamma carrier.",
        "proof_boundary": "Uniform exact canonical branch-amplitude and first-y-derivative envelopes on the saved 399 fold cells only. No finite-height logistic/Gamma amplitude identification, opposite-branch corridor integral, ordinary splice, Q_K-T bound, complete T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
        "dependencies": {name: {"path": relative(path), "sha256": file_hash(path)} for name, path in DEPENDENCIES.items()},
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "resource_policy": {"workers": 1, "process_priority": priority},
        "runtime": {"elapsed_seconds": round(time.perf_counter() - started, 3)},
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    artifact["sources"]["note"] = {"path": relative(NOTE), "sha256": file_hash(NOTE)}
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print("certified exact scaled Airy branch envelopes: |S|<0.356, |S_y|<0.258 on all 399 cells", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
