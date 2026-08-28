#!/usr/bin/env python3
"""Isolate the exact Q_K-to-Hardy-upper truncation bridge and its target."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_QK_exact_hardy_truncation_bridge_target_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
PAPER = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_Brereton_2607.15310.pdf"

DEPENDENCIES = {
    "normalization_repair": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_domain_Kummer_normalization_repair_gate.json",
    "ownership_ledger": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_physical_transform_ownership_ledger_gate.json",
    "classical_upper_split": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_classical_upper_main_remainder_split_gate.json",
    "gamma_insertion": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_uniform_gamma_target_insertion_gate.json",
    "A_transition": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_complete_endpoint_projector_assembly_gate.json",
}

PRECISION = 120
HEIGHT = 10_000_000_000
PUBLISHED_QKT_LOWER = "-0.0663408"
PUBLISHED_QKT_UPPER = "-0.00727437"
DISPLAY_DELTA_LOWER = "-0.0663407986927308950"
DISPLAY_DELTA_UPPER = "-0.0072743686927308951"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {relative(path)}")
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


def central_row(artifact: dict[str, Any]) -> dict[str, Any]:
    rows = [row for row in artifact["output_rows"] if row["output_label"] == "t+0.00"]
    require(len(rows) == 1, "central output row drift")
    require(rows[0]["target_t"] == "10000000000.00000000000000000", "central height drift")
    return rows[0]


def direct_kummer_scout(alpha_values: list[int], dps: int = 50) -> list[dict[str, Any]]:
    ctx.dps = dps
    ctx.threads = 1
    t = arb(HEIGHT)
    a = acb(arb("0.75"), -t / 2)
    b = acb(arb("0.75"), t / 2)
    beta = a.gamma() * b.gamma() / acb(arb("1.5")).gamma()
    phase = acb(0, -arb.pi() / 8).exp()
    c_t = (arb.pi() / (32 * t)) ** arb("0.25")
    rows: list[dict[str, Any]] = []
    for alpha in alpha_values:
        z = acb(0, arb.pi() * alpha * alpha / 4)
        kummer = z.hypgeom_1f1(a, arb("1.5"))
        projected = c_t * alpha * (phase * beta * kummer).real
        upper = abs(projected).upper()
        require(projected.contains(0), f"direct Kummer scout unexpectedly excludes zero at alpha={alpha}")
        require(upper > arb("1e100"), f"direct Kummer scout lost catastrophic enclosure at alpha={alpha}")
        rows.append(
            {
                "alpha": alpha,
                "corrected_label_ball": projected.str(14, more=True),
                "contains_zero": True,
                "absolute_upper_ball": upper.str(12, more=True),
                "usable_for_O_0p1_target": False,
            }
        )
    return rows


def render_note(artifact: dict[str, Any]) -> str:
    target = artifact["bridge_target_certificate"]
    return f"""# Exact Kummer-to-Hardy truncation bridge target

Date: 2026-08-27

Status: exact target isolated; not a proof of the bridge value

Use the corrected half-domain normalization from the repair gate and define

```text
U=Z-L,
rho_RS=U-T,
Delta_KU=Q_K-U.                                      (KB1)
```

Here `U` is the independently enclosed complementary Hardy upper component,
`T` is the exact classical upper main, and `Q_K` is the corrected finite
equation-(9) Kummer roster.  These definitions give the exact identities

```text
Q_K-T=Delta_KU+rho_RS,
R_KGamma=Delta_KU+rho_RS-(G-T),
J_Z=Delta_KU+rho_RS-(G-T)-A_transition.               (KB2)
```

At `t=10^10`, Arb gives

```text
rho_RS={target['rho_RS_ball']}.                        (KB3)
```

Transporting the already certified sufficient `Q_K-T` corridor through
(KB2) yields the rigorous displayed target

```text
{target['published_Delta_KU_sufficient_corridor']}    (KB4)
```

The full outward-rounded endpoints retained in the artifact are

```text
lower: {target['published_Delta_KU_safe_lower_ball']}
upper: {target['published_Delta_KU_safe_upper_ball']}
```

Thus the stronger route does not ask us to prove `Q_K` is almost equal to
`U`.  It requires a definite negative finite-roster bias: even `Delta_KU=0`
misses the upper wall by more than `{target['zero_miss_margin_lower_ball']}`.

The saved hybrid value is deliberately not substituted for `Q_K`.  As a
route diagnostic only, its hybrid-minus-exact-upper ball lies above the
required upper wall by at least
`{target['source_hybrid_miss_margin_lower_ball']}`.  This reinforces the need
for an exact truncation/continuation theorem rather than telemetry transfer.

The direct fixed-precision Arb call to `1F1` is also rejected as a production
evaluator.  All three corrected-label scouts contain zero with astronomical
enclosures; the lower endpoint already has

```text
{artifact['direct_Kummer_conditioning_scout'][0]['corrected_label_ball']}.
```

This is a diagnostic of the naive representation, not an impossibility
theorem for scaled Kummer functions, steepest descent, or a common contour.
Standard contiguous relations shift Kummer parameters, whereas the odd-label
roster keeps `(a,b)` fixed and changes
`z_alpha=i*pi*alpha^2/4`; moreover
`z_(alpha+2)-z_alpha=i*pi(alpha+1)`.  The exponential phase repeats on the
odd lattice, but `1F1(a;b;z)` is not periodic in `z`.

The selected analytic object is therefore one exact continuation bridge for
`Delta_KU`: lower omitted labels, the analytically continued upper tail,
the exact `H(t)` normalization, and any lower transition term must be
assembled with signs intact before taking absolute values.

Proof boundary: exact saved-height identities and sufficient interval target,
plus fixed-precision conditioning diagnostics only.  No enclosure of
`Delta_KU`, `Q_K`, or `J_Z`, no non-A bound, all-height theorem, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    ctx.dps = PRECISION
    ctx.threads = 1
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "dependency failure")
    require(
        dependencies["normalization_repair"]["decision"]["half_domain_physical_functional_certified"]
        is True,
        "normalization repair drift",
    )

    upper_row = central_row(dependencies["classical_upper_split"])
    classical_main = arb(upper_row["high_precision"]["classical_upper_main_ball"])
    exact_upper = arb(upper_row["high_precision"]["exact_complementary_upper_ball"])
    rho = exact_upper - classical_main
    stored_rho = arb(upper_row["exact_remaining_correction_ball"])
    require(rho.overlaps(stored_rho), "central rho_RS drift")

    ownership = dependencies["ownership_ledger"]
    corridor = ownership["corridor_certificate"]
    qkt_lower = arb(corridor["exact_Q_K_minus_T_safe_lower_ball"])
    qkt_upper = arb(corridor["exact_Q_K_minus_T_safe_upper_ball"])
    delta_exact_lower = qkt_lower.upper() - rho.lower()
    delta_exact_upper = qkt_upper.lower() - rho.upper()
    require(delta_exact_lower < delta_exact_upper, "exact Delta_KU corridor is empty")

    published_qkt_lower = arb(PUBLISHED_QKT_LOWER)
    published_qkt_upper = arb(PUBLISHED_QKT_UPPER)
    delta_published_lower = published_qkt_lower - rho.lower()
    delta_published_upper = published_qkt_upper - rho.upper()
    require(delta_published_lower > delta_exact_lower, "published Delta_KU lower endpoint is not conservative")
    require(delta_published_upper < delta_exact_upper, "published Delta_KU upper endpoint is not conservative")
    require(delta_published_upper < 0, "Delta_KU target no longer excludes zero")
    display_delta_lower = arb(DISPLAY_DELTA_LOWER)
    display_delta_upper = arb(DISPLAY_DELTA_UPPER)
    require(display_delta_lower > delta_published_lower, "displayed Delta_KU lower endpoint is not conservative")
    require(display_delta_upper < delta_published_upper, "displayed Delta_KU upper endpoint is not conservative")
    require(display_delta_lower < display_delta_upper, "displayed Delta_KU corridor is empty")

    gamma_ball = arb(
        dependencies["gamma_insertion"]["correction_certificate"][
            "Gamma_to_classical_physical_absolute_bound_ball"
        ]
    )
    gamma_uncertainty = arb(0, abs(gamma_ball).upper())
    a_transition = arb(
        dependencies["A_transition"]["certificate"]["components"][
            "projector_completed_exact_A_transition_ball"
        ]
    )
    base = rho - gamma_uncertainty - a_transition
    j_lower = arb(corridor["exact_JZ_safe_lower_ball"])
    j_upper = arb(corridor["exact_JZ_safe_upper_ball"])
    delta_j_lower = j_lower.upper() - base.lower()
    delta_j_upper = j_upper.lower() - base.upper()
    require(delta_published_lower > delta_j_lower, "published Delta_KU lower misses J_Z route")
    require(delta_published_upper < delta_j_upper, "published Delta_KU upper misses J_Z route")

    source_hybrid_gap = arb(upper_row["hybrid_minus_exact_upper_ball"])
    require(source_hybrid_gap.lower() > delta_published_upper, "source hybrid no longer lies above bridge target")
    source_hybrid_miss = source_hybrid_gap.lower() - delta_published_upper
    zero_miss = -delta_published_upper
    width = delta_published_upper - delta_published_lower

    artifact = {
        "kind": STEM,
        "status": "exact_QK_minus_exact_Hardy_upper_bridge_isolated_with_rigorous_sufficient_signed_target",
        "passed": True,
        "scope": {
            "height": HEIGHT,
            "source_labels": [159_577, 5_122_421, 2],
            "source_label_count": 2_481_423,
            "precision_decimal_digits": PRECISION,
        },
        "identity_certificate": {
            "definitions": "U=Z-L; rho_RS=U-T; Delta_KU=Q_K-U",
            "classical_bridge": "Q_K-T=Delta_KU+rho_RS",
            "Gamma_bridge": "R_KGamma=Q_K-G=Delta_KU+rho_RS-(G-T)",
            "joined_bridge": "J_Z=Delta_KU+rho_RS-(G-T)-A_transition",
            "normalization_guard": "Q_K uses the corrected c_t full-Kummer label, equivalently 2c_t Re of the half-domain integral.",
        },
        "bridge_target_certificate": {
            "classical_upper_main_ball": classical_main.str(110, more=True),
            "exact_complementary_upper_ball": exact_upper.str(110, more=True),
            "rho_RS_ball": rho.str(110, more=True),
            "stored_rho_RS_ball": stored_rho.str(110, more=True),
            "exact_Q_K_minus_T_safe_lower_ball": qkt_lower.str(100, more=True),
            "exact_Q_K_minus_T_safe_upper_ball": qkt_upper.str(100, more=True),
            "exact_Delta_KU_safe_lower_ball": delta_exact_lower.str(100, more=True),
            "exact_Delta_KU_safe_upper_ball": delta_exact_upper.str(100, more=True),
            "published_Q_K_minus_T_sufficient_corridor": f"{PUBLISHED_QKT_LOWER}<Q_K-T<{PUBLISHED_QKT_UPPER}",
            "published_Delta_KU_safe_lower_ball": delta_published_lower.str(100, more=True),
            "published_Delta_KU_safe_upper_ball": delta_published_upper.str(100, more=True),
            "published_Delta_KU_sufficient_corridor": f"{DISPLAY_DELTA_LOWER}<Delta_KU<{DISPLAY_DELTA_UPPER}",
            "Delta_KU_target_width_ball": width.str(100, more=True),
            "U_minus_G_minus_A_transition_ball": base.str(100, more=True),
            "JZ_derived_Delta_KU_safe_lower_ball": delta_j_lower.str(100, more=True),
            "JZ_derived_Delta_KU_safe_upper_ball": delta_j_upper.str(100, more=True),
            "zero_miss_margin_lower_ball": zero_miss.lower().str(80, more=True),
            "source_hybrid_minus_exact_upper_ball_forbidden_as_QK_proxy": source_hybrid_gap.str(100, more=True),
            "source_hybrid_miss_margin_lower_ball": source_hybrid_miss.lower().str(80, more=True),
        },
        "argument_lattice_audit": {
            "fixed_parameters": "The hypergeometric parameters a=3/4-it/2 and c=3/2 remain fixed over the roster.",
            "argument": "z_alpha=i*pi*alpha^2/4",
            "odd_step": "z_(alpha+2)-z_alpha=i*pi(alpha+1)=2*pi*i*((alpha+1)/2)",
            "phase_fact": "exp(z_(alpha+2))=exp(z_alpha) for odd alpha",
            "nonperiodicity_guard": "1F1(a;b;z) is not periodic in z, so exponential phase repetition is not a Kummer-value recurrence.",
            "standard_contiguous_route": "Standard Kummer contiguous relations shift the hypergeometric parameters (usually denoted a and b in DLMF); they do not directly advance this fixed-parameter z roster.",
        },
        "direct_Kummer_conditioning_scout": direct_kummer_scout([159_577, 2_641_000, 5_122_421]),
        "decision": {
            "single_exact_bridge_Delta_KU_isolated": True,
            "rigorous_sufficient_Delta_KU_corridor_nonempty": True,
            "Delta_KU_zero_excluded_by_sufficient_target": True,
            "source_hybrid_telemetry_permitted_as_QK_proxy": False,
            "standard_parameter_contiguous_recurrence_selected": False,
            "naive_direct_fixed_precision_1F1_roster_evaluator_selected": False,
            "exact_truncation_continuation_bridge_selected": True,
            "Delta_KU_enclosed": False,
            "J_Z_enclosed": False,
            "non_A_bound_proved": False,
            "rh_implication": False,
        },
        "references": {
            "paper": {"path": relative(PAPER), "sha256": file_hash(PAPER), "relevant_equations": [4, 7, 9]},
            "DLMF_Kummer_recurrences": "https://dlmf.nist.gov/13.3",
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
        "next_obligation": "Derive an exact, sign-checked representation of Delta_KU from the finite equation-(9) roster: reconstruct the exact H(t) factor from equation (4), the omitted lower-label contribution, the analytically continued upper tail/Euler-Maclaurin remainder, and any lower transition term. Assemble these against the exact lower Hardy component before absolute values, then prove the displayed negative Delta_KU corridor. Do not use source-hybrid telemetry or the superseded factor-two Kummer label.",
        "proof_boundary": "Exact saved-height bridge identities and sufficient interval target, argument-lattice route audit, and fixed-precision conditioning diagnostics only. No Delta_KU, Q_K, or J_Z enclosure, non-A bound, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("isolated exact Q_K-to-Hardy-upper bridge and rigorous signed target", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
