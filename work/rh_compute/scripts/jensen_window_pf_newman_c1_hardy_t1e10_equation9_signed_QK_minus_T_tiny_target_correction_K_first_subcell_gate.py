#!/usr/bin/env python3
"""Certify the tiny target-correction Hardy operator on the first subcell."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
for directory in (SCRIPT_DIR, VENDOR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))

from flint import acb, arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_height_derivative_identity_gate as joined
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_lower_finite_cell_K_first_subcell_scout as lower_k


HEIGHT = 10_000_000_000
RADIUS = "0.0001"
STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "tiny_target_correction_K_first_subcell_gate"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
CHECKER = SCRIPT_DIR / f"check_{STEM}.py"
DEPENDENCIES = {
    "height_derivative": ROOT / "work/rh_compute/results" / (
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
        "joined_height_derivative_identity_gate.json"
    ),
    "first_height_subcell": ROOT / "work/rh_compute/results" / (
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
        "joined_first_nonzero_height_subcell_gate.json"
    ),
    "event_atlas": ROOT / "work/rh_compute/results" / (
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_R_after_A_"
        "local_height_event_atlas_gate.json"
    ),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def complex_record(value: acb, digits: int = 70) -> dict[str, str]:
    return {
        "real_ball": value.real.str(digits, more=True),
        "imag_ball": value.imag.str(digits, more=True),
        "absolute_ball": abs(value).str(digits, more=True),
        "real_radius": value.real.rad().str(35, more=True),
        "imag_radius": value.imag.rad().str(35, more=True),
    }


def real_record(value: arb, digits: int = 70) -> dict[str, str]:
    return {
        "ball": value.str(digits, more=True),
        "lower": value.lower().str(50, more=True),
        "upper": value.upper().str(50, more=True),
    }


def correction_certificate(
    t: arb, *, reverse: bool, use_duplication_h: bool
) -> dict[str, Any]:
    """Sum the simplified correction mode by mode before interval widening."""

    pi = arb.pi()
    imaginary = acb(0, 1)
    theta, theta_prime = joined.theta_and_derivative(t)
    direct_h = joined.direct_log_H_derivative(t)
    duplication_h = joined.duplication_log_H_derivative(t)
    require(direct_h.overlaps(duplication_h), "H logarithmic derivatives miss")
    h = duplication_h if use_duplication_h else direct_h
    H, H_transport = lower_k.stable_H_box(t)
    thermal = (-2 * pi * t).exp()
    C_G = (1 + thermal) ** (-arb(1) / 2)
    C_G_prime = C_G * pi * thermal / (1 + thermal)
    alpha = C_G - H
    beta = C_G_prime - h * C_G

    D = acb(0)
    D_prime = acb(0)
    phase_covariant_target = acb(0)
    compact = acb(0)
    modes = (
        range(joined.TARGET_END, joined.TARGET_START - 1, -1)
        if reverse
        else range(joined.TARGET_START, joined.TARGET_END + 1)
    )
    for mode in modes:
        log_mode = arb(mode).log()
        term = (-imaginary * t * log_mode).exp() / arb(mode).sqrt()
        D += term
        D_prime += -imaginary * log_mode * term
        phase_term = imaginary * (theta_prime - log_mode) * term
        phase_covariant_target += phase_term
        compact += (beta + imaginary * alpha * (theta_prime - log_mode)) * term

    assembled = alpha * phase_covariant_target + beta * D
    H_prime = h * H
    raw_product_rule = (
        (C_G_prime - H_prime) * D
        + alpha * D_prime
        + (imaginary * theta_prime - h) * alpha * D
    )
    require(compact.overlaps(assembled), "compact/assembled correction miss")
    require(compact.overlaps(raw_product_rule), "compact/raw product-rule correction miss")
    projection = joined.hardy_projection(theta, compact) / H
    require(abs(projection).upper() < arb("1e-15"), "tiny correction scale widened")
    return {
        "height_ball": t.str(65, more=True),
        "precision_bits": ctx.prec,
        "summation_order": "descending" if reverse else "ascending",
        "H_log_derivative_formula": "duplication" if use_duplication_h else "direct",
        "mode_range": [joined.TARGET_START, joined.TARGET_END],
        "mode_count": joined.TARGET_END - joined.TARGET_START + 1,
        "theta": real_record(theta),
        "theta_prime": real_record(theta_prime),
        "H": real_record(H),
        "H_transport": H_transport,
        "H_log_prime": real_record(h),
        "C_G": real_record(C_G),
        "C_G_prime": real_record(C_G_prime),
        "alpha_C_G_minus_H": real_record(alpha),
        "beta_C_G_prime_minus_h_C_G": real_record(beta),
        "D_T": complex_record(D),
        "D_T_prime": complex_record(D_prime),
        "D_T_prime_plus_i_theta_prime_D_T": complex_record(
            phase_covariant_target
        ),
        "compact_weighted_mode_sum_K_corr": complex_record(compact),
        "assembled_simplified_K_corr": complex_record(assembled),
        "raw_product_rule_K_corr": complex_record(raw_product_rule),
        "all_three_forms_overlap": True,
        "Hardy_projection_over_H": real_record(projection),
        "identities": {
            "value_correction": "P_corr=(C_G-H)D_T",
            "raw_operator": (
                "L_t[P_corr]=(C_G'-H')D_T+(C_G-H)D_T'"
                "+(i theta'-H'/H)(C_G-H)D_T"
            ),
            "cancelled_operator": (
                "L_t[P_corr]=(C_G-H)(D_T'+i theta'D_T)"
                "+(C_G'-(H'/H)C_G)D_T"
            ),
            "mode_weight": (
                "[C_G'-(H'/H)C_G+i(C_G-H)(theta'-log m)]m^(-1/2-it)"
            ),
            "cancelled_terms": "-H'D_T+(H'/H)H D_T=0",
        },
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Tiny target-correction K on the first height subcell

Date: 2026-08-28

Status: rigorous exact target-correction Hardy-operator interval; complete
`K_T` assembly remains open

The A-free packet contains the finite correction

```text
P_corr=(C_G-H)D_T,
D_T=sum_(m=622)^39894 m^(-1/2-it),
C_G=(1+exp(-2*pi*t))^(-1/2).                       (TC1)
```

For `h=H'/H` and `L_t=d_t+i theta'-h`, the unsimplified product rule is

```text
L_t[P_corr]=(C_G'-H')D_T+(C_G-H)D_T'
             +(i theta'-h)(C_G-H)D_T.              (TC2)
```

The two apparent `H'D_T` terms cancel exactly, leaving

```text
L_t[P_corr]
 =(C_G-H)(D_T'+i theta'D_T)+(C_G'-h C_G)D_T

 =sum_(m=622)^39894
  [C_G'-h C_G+i(C_G-H)(theta'-log m)]m^(-1/2-it). (TC3)
```

The production evaluator sums the final line of (TC3) mode by mode on
`I_1=[10^10-10^-4,10^10+10^-4]`, so it never forms the separately widened
`D_T'` and `i theta'D_T` pieces.  Arb obtains

```text
K_corr: real {c['compact_weighted_mode_sum_K_corr']['real_ball']}
        imag {c['compact_weighted_mode_sum_K_corr']['imag_ball']},

Hardy_t[K_corr]/H={c['Hardy_projection_over_H']['ball']}.       (TC4)
```

The compact weighted sum overlaps both the simplified two-term assembly and
the raw three-term product rule.  The independent checker raises precision,
reverses all 39,273 modes, switches to the Gamma-duplication formula for
`H'/H`, and repeats all three forms.  A changed six-mode analytic model also
compares (TC3) with a direct numerical derivative of `(C-H)D`.

Pi provenance: `pi` comes only from the exact thermal Gamma coefficient,
Riemann--Siegel phase, and inherited Fourier/Gamma normalization.  No fitted
geometric constant supplies `pi`.

Proof boundary: (TC1)--(TC4) certify only the exact tiny target-correction
operator on `I_1`.  Final assembly with the transition--upper-arc,
lower-plus-ordinary, and positive-real-tail packets remains open.  No wider
`Q_K-T` sign interval, full event-cell theorem, wall handoff, all-height
theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    priority = joined.set_low_priority()
    require(priority == "below_normal_one_cpu", f"resource cap unavailable: {priority}")
    ctx.prec = 384
    ctx.threads = 1
    require(CHECKER.is_file(), f"missing checker: {CHECKER}")
    dependencies: dict[str, dict[str, Any]] = {}
    for name, path in DEPENDENCIES.items():
        payload = load_json(path)
        require(payload.get("passed") is True, f"dependency did not pass: {name}")
        dependencies[name] = {
            "path": relative(path),
            "sha256": file_hash(path),
            "artifact": payload,
        }
    require(
        dependencies["height_derivative"]["artifact"]["decision"]
        ["exact_joined_first_height_derivative_proved"]
        is True,
        "height derivative identity drift",
    )
    require(
        dependencies["first_height_subcell"]["artifact"]["decision"]
        ["first_nonzero_radius_Q_K_minus_T_sign_theorem_proved"]
        is True,
        "first height subcell drift",
    )
    t = arb(arb(HEIGHT), arb(RADIUS))
    certificate = correction_certificate(
        t, reverse=False, use_duplication_h=False
    )
    artifact = {
        "kind": "rh_c1_hardy_tiny_target_correction_K_first_height_subcell_gate",
        "date": "2026-08-28",
        "status": "tiny_target_correction_Hardy_operator_first_nonzero_height_subcell_certified",
        "passed": True,
        "resource_mode": priority,
        "certificate": certificate,
        "decision": {
            "tiny_target_correction_K_interval_proved": True,
            "tiny_target_correction_all_product_rule_forms_overlap": True,
            "complete_K_T_interval_proved": False,
            "wider_Q_K_minus_T_sign_interval_proved": False,
            "maximal_event_cell_sign_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": row["path"], "sha256": row["sha256"]}
            for name, row in dependencies.items()
        },
        "sources": {
            relative(Path(__file__).resolve()): file_hash(Path(__file__).resolve()),
            relative(CHECKER): file_hash(CHECKER),
        },
        "next_obligation": (
            "Assemble the transition--upper-arc, lower-plus-ordinary, positive-real-tail, "
            "and tiny-correction Hardy packets before a final norm and independently check "
            "complete K_T on I_1."
        ),
        "proof_boundary": (
            "Exact tiny target-correction Hardy-operator interval only on "
            "|t-10^10|<=10^-4. No complete K_T assembly, wider sign interval, full "
            "event-cell theorem, wall handoff, all-height theorem, Lambda<=0, PF-infinity, "
            "RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    projection = certificate["Hardy_projection_over_H"]["ball"]
    print(
        "certified tiny target-correction K first height subcell; "
        f"Hardy_contribution={projection}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
