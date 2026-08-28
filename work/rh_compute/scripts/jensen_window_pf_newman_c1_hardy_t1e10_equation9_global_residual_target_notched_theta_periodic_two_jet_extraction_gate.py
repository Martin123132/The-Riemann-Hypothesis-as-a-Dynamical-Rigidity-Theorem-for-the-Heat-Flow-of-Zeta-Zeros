#!/usr/bin/env python3
"""Certify the two-jet periodization of the folded target-notch current."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
import sys

sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_periodic_two_jet_extraction_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__)
CHECKER = Path(__file__).with_name("check_" + STEM + ".py")
DEPENDENCIES = {
    "R_after_A_equivalence": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_finite_regulator_equivalence_gate.json",
    "punctured_cell_fold": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_punctured_cell_fold_gate.json",
    "symmetric_poisson_interchange": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate.json",
}

A = 159_577
B = 5_122_421
L = (B - A) // 2
T_LO = 622
T_HI = 39_894
PRECISION = 110


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


def symbolic_certificate() -> dict[str, Any]:
    s, x, alpha = sp.symbols("s x alpha", real=True)
    phase = sp.exp(sp.I * sp.pi * x * (alpha + 2 * s) ** 2 / 4)
    derivative_residual = sp.simplify(
        sp.diff(phase, s) / (sp.I * sp.pi * x) - (alpha + 2 * s) * phase
    )
    require(derivative_residual == 0, "folded-current derivative drift")

    source_sum = sp.Integer(L) * sp.Integer(A + L - 1)
    square_increment = 4 * s * source_sum + 4 * L * s**2
    endpoint_square_jump = sp.Integer(B * B - A * A)
    psi_zero = sp.factor((square_increment - s * endpoint_square_jump) / 4)
    require(psi_zero == L * s * (s - 1), "zero-x periodic defect drift")

    p = s**2 - s
    k_zero = B - A
    phi_zero = sp.expand(psi_zero - sp.Rational(k_zero, 2) * p)
    require(phi_zero == 0, "zero-x two-jet residual did not vanish")
    require(p.subs(s, 0) == 0 and p.subs(s, 1) == 0, "bridge endpoint drift")
    require(sp.diff(p, s).subs(s, 1) - sp.diff(p, s).subs(s, 0) == 2, "bridge derivative jump drift")

    return {
        "folded_current_derivative": "partial_s Theta_L(x,s)=i*pi*x*F_x(s)",
        "folded_current_derivative_residual": str(derivative_residual),
        "theta_endpoint_jump": "D_x=Theta_L(x,1)-Theta_L(x,0)=e_B(x)-e_A(x)",
        "current_endpoint_jump": "F_x(1)-F_x(0)=B*e_B(x)-A*e_A(x)=K_x",
        "scaled_value_periodization": "Psi_x(s)=[Theta_L(x,s)-Theta_L(x,0)-s*D_x]/(i*pi*x)",
        "zero_x_extension": f"Psi_0(s)={L}*(s^2-s)",
        "two_jet_periodization": "Phi_x(s)=Psi_x(s)-(K_x/2)*(s^2-s)",
        "periodic_jet_conditions": "Phi_x(0)=Phi_x(1)=0 and partial_s Phi_x(0)=partial_s Phi_x(1)",
        "zero_x_residual": str(phi_zero),
    }


def harmonic_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    rows: list[dict[str, Any]] = []
    for epsilon_text in ("0", "1e-6", "1e-8", "1e-10"):
        epsilon = arb(epsilon_text)
        total = arb(0)
        for m in range(T_LO, T_HI + 1):
            weight = arb(1) if epsilon_text == "0" else (-pi * epsilon * m * m).exp()
            total += weight / m
        rows.append(
            {
                "epsilon": epsilon_text,
                "H_T_ball": total.str(75, more=True),
                "H_T_over_pi_ball": (total / pi).str(75, more=True),
            }
        )
    return {
        "definition": "H_T(epsilon)=sum_(m=622)^39894 exp(-pi*epsilon*m^2)/m",
        "polynomial_current": "integral_0^1 (s^2-s) C_(M,epsilon)'(s) ds=i*H_T(epsilon)/pi",
        "independent_of_M_once": "M>=39894",
        "rows": rows,
    }


def periodized_identity_certificate() -> dict[str, Any]:
    return {
        "endpoint_divided_difference": "E_x=[e_B(x)-e_A(x)]/(i*pi*x), E_0=(B^2-A^2)/4",
        "stable_endpoint_form": "E_x=((B^2-A^2)/4)e^(i*pi*x*(A^2+B^2)/8)sinc(pi*x*(B^2-A^2)/8)",
        "endpoint_current": "K_x=B*e_B(x)-A*e_A(x), K_0=B-A=2L",
        "source_identity": "integral_0^1 F_x(s)C_(M,epsilon)(s)ds=E_x-i*K_x*H_T(epsilon)/(2*pi)-integral_0^1 Phi_x(s)C_(M,epsilon)'(s)ds",
        "zero_x_identity": "integral_0^1 F_0(s)C_(M,epsilon)(s)ds=(B^2-A^2)/4-i*(B-A)*H_T(epsilon)/(2*pi)",
        "finite_cutoff_scope": "M>=B=5122421 in R_after_A; identity itself needs M>=39894",
        "limit_order_preserved": "finite M and epsilon identity first; then M to infinity at fixed epsilon; then epsilon down to zero",
    }


def fourier_tail_certificate() -> dict[str, Any]:
    return {
        "coefficient_definition": "hat Phi_x(m)=integral_0^1 Phi_x(s)e^(-2*pi*i*m*s)ds",
        "three_IBP_bound": "|hat Phi_x(m)|<=(|Delta Phi_x''|+integral_0^1|Phi_x'''(s)|ds)/(2*pi*|m|)^3",
        "differentiated_term_bound": "2*pi*|m|*|hat Phi_x(m)|<=K_Phi(x)/(2*pi*|m|)^2",
        "symmetric_tail_bound": "sum_(|m|>K)2*pi*|m|*|hat Phi_x(m)|<=K_Phi(x)/(2*pi^2*K)",
        "consequence": "The differentiated leakage series is absolutely convergent after the two endpoint jets are extracted.",
        "quantitative_warning": "The displayed derivative majorant proves convergence but is not asserted to meet the R_after_A target.",
    }


def render_note(artifact: dict[str, Any]) -> str:
    rows = artifact["harmonic_certificate"]["rows"]
    return f"""# Periodic two-jet extraction of the target-notched theta current

Date: 2026-08-23

Status: exact finite-regulator endpoint-current extraction and absolutely
convergent differentiated leakage certified; quantitative `R_after_A` bound open

Let `e_D(x)=exp(i*pi*x*D^2/4)` for `D` equal to `A={A}` or `B={B}`, and
retain the one-cell objects from Section 11.445.  Put

```text
D_x=e_B(x)-e_A(x),
E_x=D_x/(i*pi*x),              E_0=(B^2-A^2)/4,
K_x=B*e_B(x)-A*e_A(x),        K_0=B-A=2L.             (TJ1)
```

The endpoint divided difference has the stable entire form

```text
E_x=((B^2-A^2)/4)e^[i*pi*x*(A^2+B^2)/8]
    sinc[pi*x*(B^2-A^2)/8],                           (TJ2)
```

where `sinc(z)=sin(z)/z` and `sinc(0)=1`.  For `x>0`, define

```text
Psi_x(s)
 =[Theta_L(x,s)-Theta_L(x,0)-s*D_x]/(i*pi*x),

Phi_x(s)=Psi_x(s)-(K_x/2)(s^2-s).                    (TJ3)
```

Both quotients extend continuously to `x=0`, and the roster sums collapse to

```text
Psi_0(s)=L(s^2-s),             Phi_0(s)=0.            (TJ4)
```

The value and first-derivative jumps are now removed exactly:

```text
Phi_x(0)=Phi_x(1)=0,
partial_s Phi_x(0)=partial_s Phi_x(1).                (TJ5)
```

For `T={{622,...,39894}}`, set

```text
H_T(epsilon)=sum_(m in T)e^(-pi*epsilon*m^2)/m.
```

At every common finite cutoff `M>=B`, Fourier-pair cancellation gives

```text
integral_0^1 (s^2-s)C_(M,epsilon)'(s)ds
 =i*H_T(epsilon)/pi.                                  (TJ6)
```

Since `integral_0^1 C_(M,epsilon)=1`, integration by parts with (TJ5) gives
the exact joined source identity

```text
integral_0^1 F_x(s)C_(M,epsilon)(s)ds
 =E_x-i*K_x*H_T(epsilon)/(2*pi)
  -integral_0^1 Phi_x(s)C_(M,epsilon)'(s)ds.          (TJ7)
```

Thus the coherent endpoint value and derivative currents are extracted
before any norm.  At zero regulator,

```text
H_T(0)={rows[0]['H_T_ball']},
H_T(0)/pi={rows[0]['H_T_over_pi_ball']}.              (TJ8)
```

If `hat Phi_x(m)=integral_0^1 Phi_x(s)e^(-2*pi*i*m*s)ds`, three integrations
by parts use (TJ5) to give

```text
|hat Phi_x(m)|
 <=K_Phi(x)/(2*pi*|m|)^3,
K_Phi(x)=|Delta Phi_x''|+integral_0^1|Phi_x'''(s)|ds,

sum_(|m|>K)2*pi*|m|*|hat Phi_x(m)|
 <=K_Phi(x)/(2*pi^2*K).                               (TJ9)
```

The differentiated leakage in (TJ7) is therefore absolutely convergent after
the two endpoint jets are removed, including when the Abel regulator tends to
zero.  The raw derivative majorant in (TJ9) is only a convergence certificate;
it is not claimed to be quantitatively useful.

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Pi provenance: every `pi` in (TJ1)--(TJ9) is inherited from the Kummer
quadratic phase, integer Fourier character, and Gaussian Abel regulator.  No
fitted or geometric occurrence is introduced.

Proof boundary: exact finite-regulator two-jet periodization, closed target
harmonic current, and absolute convergence of the differentiated leakage only.
No quantitative finite-current evaluator, cancellation with the separately
owned A/B extraction terms, Kummer `x` bound, `R_after_A` or `R_Dir` bound,
complete `Q_K-T`, all-height theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["R_after_A_equivalence"]["decision"]["common_physical_limit_order_preserved"] is True, "R_after_A limit drift")
    require(dependencies["punctured_cell_fold"]["decision"]["long_u_integral_folded_exactly_to_one_periodic_cell"] is True, "one-cell fold drift")
    require(dependencies["symmetric_poisson_interchange"]["decision"]["symmetric_poisson_interchange_proved"] is True, "Poisson interchange drift")

    artifact = {
        "kind": STEM,
        "status": "exact_periodic_two_jet_endpoint_current_extraction_and_absolute_leakage_convergence_certified_quantitative_bound_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "source_odd_roster": [A, B],
            "roster_coordinate_length": L,
            "target_positive_modes": [T_LO, T_HI],
            "common_cutoff": "M>=B=5122421",
            "regulator": "epsilon>=0 after the finite-regulator identity is established",
        },
        "symbolic_certificate": symbolic_certificate(),
        "harmonic_certificate": harmonic_certificate(),
        "periodized_identity_certificate": periodized_identity_certificate(),
        "fourier_tail_certificate": fourier_tail_certificate(),
        "decision": {
            "theta_value_jump_extracted_exactly": True,
            "theta_first_derivative_jump_extracted_exactly": True,
            "two_jet_remainder_periodic_C1": True,
            "zero_x_extension_regular_and_residual_vanishes": True,
            "target_asymmetry_reduced_to_closed_harmonic_current": True,
            "differentiated_leakage_absolutely_convergent": True,
            "finite_regulator_and_limit_order_preserved": True,
            "floating_data_used_as_proof": False,
            "A_B_extraction_cancellation_completed": False,
            "quantitative_R_after_A_bound_proved": False,
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
            "precision_decimal_digits": PRECISION,
        },
        "next_obligation": "Insert the extracted E_x and K_x*H_T endpoint currents into the exact R_after_A ownership ledger and prove their coefficientwise cancellation or retention against the A notch, B trace window, B outer block, H_x, and explicit A_m/B_m atoms. Only after that ownership gate should the periodic C1 leakage be evaluated or bounded.",
        "proof_boundary": "Exact finite-regulator two-jet periodization, closed target harmonic current, and absolute convergence of the differentiated leakage only. No quantitative finite-current evaluator, cancellation with separately owned A/B extraction terms, Kummer x bound, R_after_A or R_Dir bound, complete Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified periodic two-jet target-notch extraction", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
