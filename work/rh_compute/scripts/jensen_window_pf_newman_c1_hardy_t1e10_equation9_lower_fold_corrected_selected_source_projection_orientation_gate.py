#!/usr/bin/env python3
"""Fix the source normalization and real-projection sign of W_corr.

This gate deliberately treats W_corr as a weighted leading-defect statistic,
not as the unweighted selector strip or the isolated fold-owned residual.
"""

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

from flint import arb, acb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_corrected_selected_source_projection_orientation_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "joint_normal_form": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_grouped_characteristic_airy_boundary_normal_form_gate.json",
    "common_profile_operator": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_finite_t_common_profile_operator_gate.json",
    "corrected_selected_splice": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_corrected_selected_defect_splice_gate.json",
    "half_reflection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_kummer_reflection_branch_reduction_gate.json",
    "paired_target_residual": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_paired_target_residual_normal_form_gate.json",
    "event_pairing": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_event_ordered_398_pairing_multiplier_gate.json",
}

C = 159_577
PRECISION = 80


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


def symbolic_certificate() -> dict[str, str | int]:
    pi = sp.pi
    raw_integral = sp.symbols("I")
    canonical = raw_integral / (2 * pi)
    physical_half = sp.sqrt(2) * raw_integral / pi
    ratio = sp.simplify(physical_half / canonical)
    require(ratio == 2 * sp.sqrt(2), "half-source/canonical normalization drift")

    parity_integer = (C * C - 1) // 8
    require(C % 8 == 1, "selector endpoint is not 1 modulo 8")
    require((C * C - 1) % 8 == 0 and parity_integer % 2 == 0, "source carrier parity failed")

    w_real, w_imag = sp.symbols("W_real W_imag", real=True)
    phase = sp.pi * (sp.Integer(C) ** 2 - 1) / 8
    projected = sp.expand_complex(
        2 * sp.re(sp.exp(sp.I * phase) * 2 * sp.sqrt(2) * (w_real + sp.I * w_imag))
    )
    require(sp.simplify(projected - 4 * sp.sqrt(2) * w_real) == 0, "real source projection drift")
    return {
        "canonical_coefficient": "G_m=(1/(2pi)) Integral_0^Y J_m(y)dy",
        "physical_half_mode": "K_m^(half)=(sqrt(2)/pi) Integral_0^Y J_m(y)dy=2sqrt(2)G_m",
        "common_carrier": "exp(i beta^3), beta^3=pi*C^2/8",
        "half_domain_projection": "2 Re[exp(-i*pi/8) K_half]",
        "parity_integer": parity_integer,
        "parity_identity": "exp(i*pi*(C^2-1)/8)=1",
        "corrected_selected_projection": "R_corr=4sqrt(2) Re W_corr",
        "scope": "This is the non-event-zero corrected selected leading-defect component, not the complete fold block or Q_K-T residual.",
    }


def numerical_certificate(splice: dict[str, Any]) -> dict[str, str]:
    ctx.dps = PRECISION
    source_factor = 4 * arb(2).sqrt()
    certificate = splice["certificate"]
    canonical_real_upper = arb(certificate["corrected_selected_uniform_real_part_upper"])
    canonical_modulus = arb(certificate["corrected_selected_uniform_modulus_bound"])
    source_real_upper = (source_factor * canonical_real_upper).upper()
    source_modulus = (source_factor * canonical_modulus).upper()
    require(source_real_upper < -arb("1.263e-4"), "source-projected corrected real part is not below -1.263e-4")
    require(source_modulus < arb("2.314e-4"), "source-projected corrected modulus exceeds 2.314e-4")
    return {
        "source_projection_factor_4sqrt2": source_factor.str(PRECISION, more=True),
        "canonical_corrected_real_part_upper": canonical_real_upper.str(PRECISION, more=True),
        "canonical_corrected_modulus_bound": canonical_modulus.str(PRECISION, more=True),
        "source_projected_corrected_real_part_upper": source_real_upper.str(PRECISION, more=True),
        "source_projected_corrected_modulus_bound": source_modulus.str(PRECISION, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Corrected selected source-projection orientation

Date: 2026-08-13

Status: exact normalization and conditional signed projection theorem; the
embedding into the isolated fold-owned paired-residual block remains open

Let `J_m(y)` denote the raw transformed fold integral with its common source
carrier suppressed.  The common-profile Fourier normalization is

```text
G_m=(1/(2pi)) Integral_0^Y J_m(y)dy.                 (SO1)
```

The exact joint Kummer normal form carries amplitude `sqrt(2)/pi`, so the
original half-domain source mode is

```text
K_m^(half)=(sqrt(2)/pi)Integral_0^Y J_m(y)dy
          =2sqrt(2)G_m.                              (SO2)
```

Restoring the suppressed carrier and applying the exact half-domain
reflection gives

```text
2 Re[e^(-i*pi/8)e^(i*beta^3) 2sqrt(2)W_corr].       (SO3)
```

Here `beta^3=pi*C^2/8`, `C={C}=1 mod 8`, and
`(C^2-1)/8={artifact['symbolic_certificate']['parity_integer']}` is even.
Therefore the phase in (SO3) is exactly one and

```text
R_corr=4sqrt(2) Re W_corr.                           (SO4)
```

Using the signed theorem of Section 11.402 proves uniformly on the ordinary
top corridor

```text
R_corr<{c['source_projected_corrected_real_part_upper']}<-1.263e-4,
|R_corr|<{c['source_projected_corrected_modulus_bound']}<2.314e-4. (SO5)
```

Equation (SO5) fixes the physical scale and sign of the weighted corrected
selected leading-defect statistic on the augmented roster.  It does not
identify that statistic with the unweighted selector strip, or prove that it
embeds unchanged into the isolated `39695..39894` fold-owned paired-residual block:
the common-profile roster is `39696..40094`, and its completed source
remainder reaches the zero, negative, and outer-positive sectors.  An exact
allocation/embedding identity must be proved before (SO5) may be called a
signed component of `Q_K-T`.

Pi provenance is explicit in (SO1)--(SO4): inverse Airy Fourier
normalization, the equation-(9) Kummer amplitude, odd-square reflection, and
`beta^3=pi*C^2/8`.  No fitted constant is used.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No weighted-statistic-to-raw-strip identity, selector-to-paired-residual
embedding, complete fold-owned residual, source/initial-data cancellation,
complete `Q_K-T` or `T_upper`, all-corridor
or height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    for name, dependency in dependencies.items():
        require(dependency.get("passed") is True, f"dependency failed: {name}")
    require(dependencies["half_reflection"]["decision"]["full_source_main_equals_twice_real_half_integral"] is True, "half reflection drift")
    require(dependencies["paired_target_residual"]["decision"]["finite_paired_half_domain_target_residual_identity_proved"] is True, "paired residual drift")
    require(dependencies["corrected_selected_splice"]["decision"]["corrected_selected_real_part_uniformly_negative"] is True, "corrected sign drift")

    symbolic = symbolic_certificate()
    certificate = numerical_certificate(dependencies["corrected_selected_splice"])
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "corrected_selected_source_projection_scale_and_negative_sign_certified_embedding_open",
        "passed": True,
        "symbolic_certificate": symbolic,
        "certificate": certificate,
        "decision": {
            "canonical_to_physical_half_mode_factor_is_2sqrt2": True,
            "source_carrier_times_half_projection_phase_is_exactly_one": True,
            "corrected_selected_augmented_strip_projection_has_negative_real_sign": True,
            "corrected_selected_weighted_leading_defect_projection_has_negative_real_sign": True,
            "corrected_selected_statistic_equals_unweighted_selector_strip": False,
            "corrected_selected_component_embedded_in_fold_owned_paired_residual": False,
            "complete_fold_owned_residual_bounded": False,
            "complete_Q_K_minus_T_bounded": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)} for name, path in DEPENDENCIES.items()
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
            "Relate the weighted corrected statistic to the unweighted selector strip with an exact companion "
            "decomposition, retaining the common source carrier, before attempting a fold-owned sign argument."
        ),
        "proof_boundary": (
            "Exact source normalization and negative real projection for the corrected selected weighted leading-defect "
            "statistic only. No identification with the unweighted selector strip, embedding into the isolated fold-owned "
            "paired residual, complete fold-owned residual, source/initial-data "
            "cancellation, complete Q_K-T or T_upper, all-corridor or height-uniform theorem, Lambda<=0, PF-infinity, RH, "
            "or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("certified weighted leading-defect projection: real part < -1.263e-4; embedding remains open", flush=True)


if __name__ == "__main__":
    main()
