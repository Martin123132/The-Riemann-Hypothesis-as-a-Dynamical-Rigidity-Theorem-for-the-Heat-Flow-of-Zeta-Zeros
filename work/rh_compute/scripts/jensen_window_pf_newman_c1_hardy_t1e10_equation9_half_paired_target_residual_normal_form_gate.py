#!/usr/bin/env python3
"""Certify the exact paired half-domain target-residual normal form."""

from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path
from typing import Any

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_paired_target_residual_normal_form_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "symmetric_poisson": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate.json",
    "symmetric_residual": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_symmetric_endpoint_tail_reassembly_gate.json",
    "half_reflection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_kummer_reflection_branch_reduction_gate.json",
    "Gamma_bulk": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_global_morse_gamma_bulk_gate.json",
    "coverage_ledger": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_ordinary_mode_coverage_ledger_gate.json",
    "B_tangent_barrier": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_kummer_B_crossing_tangent_profile_barrier_gate.json",
    "midpoint_guard": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_boundary_linear_roster_pair_anchor_gate.json",
}

LOW_END = 621
TARGET_START = 622
ORDINARY_END = 39_694
FOLD_START = 39_695
TARGET_END = 39_894
OUTER_START = 39_895


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
    q_real, q_imag, tau_real, tau_imag = sp.symbols(
        "q_real q_imag tau_real tau_imag", real=True
    )
    c = sp.cos(sp.pi / 8)
    s = sp.sin(sp.pi / 8)

    # omega=exp(-i*pi/8), omega_inv=exp(i*pi/8).  Multiplying the
    # complex target carrier by omega_inv allocates it inside Q_half without
    # asserting that it is a modewise half of the full-line Gamma bulk.
    source_projection = 2 * (c * q_real + s * q_imag)
    target_projection = 2 * tau_real
    rotated_target_real = c * (c * tau_real - s * tau_imag) + s * (
        s * tau_real + c * tau_imag
    )
    residual_projection = 2 * (c * q_real + s * q_imag - rotated_target_real)
    phase_identity = sp.simplify(
        sp.trigsimp(source_projection - target_projection - residual_projection).rewrite(sp.sqrt)
    )
    require(
        phase_identity == 0,
        "phase-allocated target identity failed",
    )

    m_cut = sp.symbols("M", integer=True, positive=True)
    count = (
        1
        + LOW_END
        + (ORDINARY_END - TARGET_START + 1)
        + (TARGET_END - FOLD_START + 1)
        + (m_cut - OUTER_START + 1)
    )
    require(sp.simplify(count - (m_cut + 1)) == 0, "nonnegative pair partition failed")

    half, zero, low, ordinary, fold, outer, target = sp.symbols(
        "H K_0 L O F X Tau"
    )
    q_half = half + zero + low + ordinary + fold + outer
    grouped = half + zero + low + (ordinary + fold - target) + outer
    require(sp.expand(q_half - target - grouped) == 0, "grouped residual failed")

    delta_f, delta_f1, remainder_plus, remainder_minus, mode = sp.symbols(
        "Delta_f Delta_f1 R_plus R_minus m", nonzero=True
    )
    denominator = 2 * sp.pi * sp.I * mode
    positive = -delta_f / denominator - delta_f1 / denominator**2 + remainder_plus / denominator**2
    negative = positive.subs({mode: -mode, remainder_plus: remainder_minus})
    pair = sp.simplify(positive + negative)
    require(sp.simplify(sp.diff(pair, delta_f)) == 0, "one-over-m current survived pairing")

    return {
        "half_mode": "K_m(t)=integral_0^(1/2) W_t(x) I_m(x) dx",
        "half_current": "K_H(t)=integral_0^(1/2) W_t(x) H_x dx",
        "finite_half_sum": "Q_(half,M)=K_H+K_0+sum_(m=1)^M[K_m+K_(-m)]",
        "physical_projection": "Q_main,M=2 Re{exp(-i*pi/8) Q_(half,M)}",
        "target_carrier": "If T_target=2 Re sum_(m=622)^39894 tau_m, put tauhat_m=exp(i*pi/8)tau_m.",
        "target_allocation_identity": "Q_main,M-T_target=2 Re{exp(-i*pi/8)[Q_(half,M)-sum tauhat_m]}.",
        "pair_identity": "K_m+K_(-m)=integral_0^(1/2)W_t(x)[I_m(x)+I_(-m)(x)]dx",
        "paired_ibp": "I_m+I_(-m)=[-2 Delta f_x'+R_m+R_(-m)]/(2*pi*i*m)^2",
        "allocation_guard": "tauhat_m is an algebraic target allocation inside the real projection, not a claim that the full-line Gamma bulk splits modewise into equal x-halves.",
        "midpoint_guard": "The pair is kept intact through the x integral; it is not set to zero at x=1/2.",
    }


def partition() -> dict[str, Any]:
    return {
        "cutoff_condition": f"integer M>={OUTER_START}",
        "endpoint_half_current": "K_H",
        "zero_mode": {"range": [0, 0], "count": 1},
        "low_B_side_pairs": {"positive_range": [1, LOW_END], "count": LOW_END},
        "ordinary_target_pairs": {
            "positive_range": [TARGET_START, ORDINARY_END],
            "count": ORDINARY_END - TARGET_START + 1,
            "target_subtraction": "tauhat_m for every mode",
        },
        "fold_owned_target_pairs": {
            "positive_range": [FOLD_START, TARGET_END],
            "count": TARGET_END - FOLD_START + 1,
            "target_subtraction": "tauhat_m for every mode",
        },
        "outer_pairs": {
            "positive_range": [OUTER_START, "M"],
            "count": f"M-{OUTER_START - 1}",
        },
        "count_identity": f"1+{LOW_END}+{ORDINARY_END - TARGET_START + 1}+{TARGET_END - FOLD_START + 1}+(M-{OUTER_START - 1})=M+1",
        "negative_modes": "Every negative mode -m occurs exactly once inside the pair K_m+K_(-m).",
        "target_count": TARGET_END - TARGET_START + 1,
    }


def render_note() -> str:
    return f"""# Exact paired half-domain target residual

Date: 2026-08-13

Status: exact-lemma certificate; not a proof of quantitative pair-block bounds

After the odd-alpha reflection, define

```text
K_m(t)=integral_0^(1/2) W_t(x)I_m(x)dx,
K_H(t)=integral_0^(1/2) W_t(x)H_x dx.                 (PR1)
```

For every symmetric cutoff `M>={OUTER_START}` the source half-sum is

```text
Q_half,M=K_H+K_0+sum_(m=1)^M (K_m+K_-m),
Q_main,M=2 Re[exp(-i*pi/8)Q_half,M].                  (PR2)
```

Let `tau_m` be the already normalized complex classical carrier, so that

```text
T_target=2 Re sum_(m={TARGET_START})^{TARGET_END} tau_m,
tauhat_m=exp(i*pi/8)tau_m.                            (PR3)
```

Pure finite algebra now gives the cancellation-preserving residual

```text
Q_main,M-T_target=2 Re exp(-i*pi/8) {{
 K_H+K_0
 +sum_(m=1)^{LOW_END}                 (K_m+K_-m)
 +sum_(m={TARGET_START})^{ORDINARY_END}(K_m+K_-m-tauhat_m)
 +sum_(m={FOLD_START})^{TARGET_END}   (K_m+K_-m-tauhat_m)
 +sum_(m={OUTER_START})^M             (K_m+K_-m) }}.   (PR4)
```

The four ranges contain respectively `{LOW_END}`, `{ORDINARY_END - TARGET_START + 1}`,
`{TARGET_END - FOLD_START + 1}`, and `M-{OUTER_START - 1}` positive modes.  Together with
`K_0` they contain exactly `M+1` nonnegative indices, while every negative
index occurs exactly once in its symmetric pair.  The target is subtracted
exactly once over its `{TARGET_END - TARGET_START + 1}` modes.

This form preserves the endpoint-current cancellation:

```text
K_m+K_-m=integral_0^(1/2)W_t(x)(I_m+I_-m)dx,
I_m+I_-m=[-2 Delta f_x'+R_m+R_-m]/(2*pi*i*m)^2.       (PR5)
```

The `1/m` current is absent before any absolute value.  The proved symmetric
Poisson interchange supplies the limit `M->infinity`, so (PR4) is also the
exact limiting residual when its outer block is understood symmetrically.

Two guards are essential.  First, `tauhat_m` in (PR3)--(PR4) is only an
algebraic allocation of the real target projection.  It does not assert that
the full-line Gamma bulk has equal modewise `x<1/2` and `x>1/2` pieces.
Second, (PR5) must retain the canonical quadratic chirp: discrete odd-square
parity at `x=1/2` does not make `I_m+I_-m` vanish there.

The tangent B profile therefore cannot be promoted by itself.  Its grouped
`q`-density must be combined with the low pairs and with the nonpositive and
outer-positive completion in (PR4).  The ordinary block remains separate
from the 200 fold-owned target pairs, matching the certified ownership
ledger.

Pi provenance: `pi/8` comes from exact odd-square Kummer reflection; all
other `pi` factors come from the equation-(9) phase or Fourier-Poisson
characters.  No fitted or geometric constant is inserted.

Proof boundary: exact finite and limiting paired residual normal form only.
No quantitative block estimate, nonlinear B-crossing theorem, A-fold splice,
complete `T_upper`, height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(
        dependencies["symmetric_poisson"]["decision"]["paired_one_over_m_endpoint_current_cancels"] is True,
        "pair-cancellation dependency drift",
    )
    require(
        dependencies["symmetric_residual"]["decision"]["symmetric_residual_limit_proved"] is True,
        "symmetric-limit dependency drift",
    )
    require(
        dependencies["half_reflection"]["decision"]["full_source_main_equals_twice_real_half_integral"] is True,
        "half-reflection dependency drift",
    )
    require(
        dependencies["Gamma_bulk"]["decision"]["ordinary_full_line_bulk_closed_at_saved_height"] is True,
        "Gamma-bulk dependency drift",
    )
    require(
        dependencies["midpoint_guard"]["decision"]["canonical_continuous_interpolation_is_linear"] is False,
        "midpoint interpolation guard drift",
    )
    require(
        dependencies["midpoint_guard"]["decision"]["all_nonzero_symmetric_Poisson_pairs_cancel_at_midpoint"] is False,
        "midpoint pair guard drift",
    )
    coverage = dependencies["coverage_ledger"]["mode_coverage"]["proposed_disjoint_target_ownership"]
    require(
        coverage["ordinary_core"]["range"] == [TARGET_START, ORDINARY_END],
        "ordinary ownership drift",
    )
    require(
        coverage["fold_owned_target"]["range"] == [FOLD_START, TARGET_END],
        "fold ownership drift",
    )

    artifact = {
        "kind": STEM,
        "status": "exact_paired_half_domain_target_residual_normal_form_proved_quantitative_bounds_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "symmetric_cutoff": f"M>={OUTER_START}",
            "target_modes": [TARGET_START, TARGET_END],
            "ordinary_owned_modes": [TARGET_START, ORDINARY_END],
            "fold_owned_modes": [FOLD_START, TARGET_END],
        },
        "symbolic_certificate": symbolic_certificate(),
        "partition": partition(),
        "decision": {
            "finite_paired_half_domain_target_residual_identity_proved": True,
            "symmetric_limit_of_paired_normal_form_proved": True,
            "every_negative_mode_retained_exactly_once": True,
            "zero_mode_and_endpoint_half_current_retained": True,
            "target_carrier_subtracted_exactly_once": True,
            "target_allocation_is_modewise_half_Gamma_bulk_claim": False,
            "midpoint_discrete_parity_cancels_Poisson_pairs": False,
            "B_tangent_profile_alone_closes_residual": False,
            "paired_blocks_quantitatively_bounded": False,
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
            "sympy_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Derive a phase-adapted enclosure for the coupled low/B-crossing and nonpositive/outer pair blocks in (PR4), retaining the q-density current and the canonical interpolation chirp; then splice the 200 fold-owned pairs to the certified atlas.",
        "proof_boundary": "Exact finite and symmetric-limit paired half-domain target-residual normal form only. No quantitative pair-block bound, nonlinear B-crossing theorem, A-fold splice, complete T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(), encoding="utf-8")
    print("certified exact paired half-domain target residual normal form", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
