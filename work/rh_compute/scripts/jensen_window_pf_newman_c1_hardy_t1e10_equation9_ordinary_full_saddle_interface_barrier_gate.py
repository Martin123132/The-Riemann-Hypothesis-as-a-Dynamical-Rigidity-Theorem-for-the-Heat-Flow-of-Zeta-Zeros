#!/usr/bin/env python3
"""Certify that the bare full-saddle carrier cannot be the ordinary/fold splice."""

from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path
from typing import Any

from flint import acb, arb, ctx


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_full_saddle_interface_barrier_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "coverage_ledger": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_ordinary_mode_coverage_ledger_gate.json",
    "fresnel_partition": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fresnel_endpoint_normal_form_gate.json",
    "universal_morse": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_universal_logistic_morse_characteristic_fold_reduction_gate.json",
}

PRECISION = 90
T = 10_000_000_000
A = 159_577
B = 5_122_421
MODES = (39_694, 39_695, 39_852, 39_853)
TARGET_NORMALIZED = arb("0.000019")
TARGET_PHYSICAL = arb("0.0000086")


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


def fresnel(z: arb, pi: arb, imaginary_unit: acb) -> acb:
    scale = (-imaginary_unit * pi / 4).exp() * (pi / 2).sqrt()
    return (imaginary_unit * pi / 4).exp() / arb(2).sqrt() * (scale * z).erf()


def saddle_current_row(mode: int) -> dict[str, Any]:
    pi = arb.pi()
    imaginary_unit = acb(0, 1)
    m = arb(mode)
    t = arb(T)
    x = 2 * pi * m**2 / (t + 2 * pi * m**2)
    q_a = (x / 2).sqrt() * (arb(A) - 2 * m / x)
    q_b = (x / 2).sqrt() * (arb(B) - 2 * m / x)

    i0 = fresnel(q_b, pi, imaginary_unit) - fresnel(q_a, pi, imaginary_unit)
    i1 = ((imaginary_unit * pi * q_b**2 / 2).exp() - (imaginary_unit * pi * q_a**2 / 2).exp()) / (imaginary_unit * pi)
    full_line_i0 = 1 + imaginary_unit
    ratio = i0 / full_line_i0 + x.sqrt() * i1 / (m * arb(2).sqrt() * full_line_i0)
    ratio_defect = abs(ratio - 1)
    normalized = ratio_defect / m.sqrt()
    physical = arb(2).sqrt() / pi * normalized

    return {
        "mode": mode,
        "x_saddle_ball": x.str(PRECISION, more=True),
        "q_A_saddle_ball": q_a.str(PRECISION, more=True),
        "q_B_saddle_ball": q_b.str(PRECISION, more=True),
        "incomplete_grouped_current_ratio_ball": complex_record(ratio),
        "ratio_defect_absolute_ball": ratio_defect.str(PRECISION, more=True),
        "carrier_scaled_normalized_defect_ball": normalized.str(PRECISION, more=True),
        "carrier_scaled_physical_defect_ball": physical.str(PRECISION, more=True),
        "normalized_target_multiple_ball": (normalized / TARGET_NORMALIZED).str(PRECISION, more=True),
        "physical_target_multiple_ball": (physical / TARGET_PHYSICAL).str(PRECISION, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    rows = {row["mode"]: row for row in artifact["certificate"]["rows"]}
    edge = rows[39_695]
    structural = rows[39_852]
    return f"""# Bare full-saddle interface barrier

Date: 2026-08-13

Status: exact leading-current obstruction certified; not a proof of a full mode-error lower bound

At the reduced Morse saddle `x_m`, retain the complete grouped alpha current

```text
G_m=m*sqrt(2)*x_m^(-3/2) I_0+x_m^(-1) I_1,
I_0=integral_(q_A)^(q_B) exp(i*pi*q^2/2)dq,
I_1=integral_(q_A)^(q_B) q exp(i*pi*q^2/2)dq.
```

The full-line current is
`G_m^full=m*sqrt(2)*x_m^(-3/2)(1+i)`.  Therefore the exact frozen-saddle
ratio is

```text
R_m=I_0/(1+i)+sqrt(x_m) I_1/[m sqrt(2)(1+i)].
```

At the proposed atlas edge `m=39695`,

```text
q_A*={edge['q_A_saddle_ball']},
|R_m-1|={edge['ratio_defect_absolute_ball']},
|R_m-1|/sqrt(m)={edge['carrier_scaled_normalized_defect_ball']}.
```

That normalized leading defect is more than
`{edge['normalized_target_multiple_ball']}` times the local `1.9e-5` target.
At the true lower-interior edge `m=39852`,

```text
q_A*={structural['q_A_saddle_ball']},
|R_m-1|={structural['ratio_defect_absolute_ball']},
|R_m-1|/sqrt(m)={structural['carrier_scaled_normalized_defect_ball']},
physical={structural['carrier_scaled_physical_defect_ball']}.
```

The mode is a half-saddle to leading order: the normalized defect exceeds
the local target by `{structural['normalized_target_multiple_ball']}` and the
physical defect exceeds its target by `{structural['physical_target_multiple_ball']}`.

This rigorously rules out a **termwise splice whose ordinary overlap object is
only the bare full-line carrier**.  It does not lower-bound the complete
`x`-integrated exact-mode error, because amplitude variation and neighboring
modes may cancel.  The admissible next object must retain the incomplete
Fresnel plus endpoint current through the ordinary Morse step, or aggregate
the signed endpoint defects before taking absolute values.

Pi provenance: every occurrence of `pi` comes from the equation-(9) Kummer
quadratic phase, its Fourier-Poisson dual mode, and the already derived
equation-(9) normalization.  No geometric fit is introduced here.

Proof boundary: exact grouped alpha-current ratios at four saved-height
interface witnesses only.  No full `x`-integral lower bound, uniform ordinary
Morse remainder, complete `T_upper`, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    ctx.dps = PRECISION
    ctx.threads = 1

    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    rows = [saddle_current_row(mode) for mode in MODES]
    by_mode = {row["mode"]: row for row in rows}

    require(arb(by_mode[39_695]["ratio_defect_absolute_ball"]) > arb("0.20"), "atlas-edge ratio barrier failed")
    require(arb(by_mode[39_695]["carrier_scaled_normalized_defect_ball"]) > 55 * TARGET_NORMALIZED, "atlas-edge target barrier failed")
    require(arb(by_mode[39_852]["ratio_defect_absolute_ball"]) > arb("0.49"), "structural half-saddle barrier failed")
    require(arb(by_mode[39_852]["carrier_scaled_normalized_defect_ball"]) > 131 * TARGET_NORMALIZED, "structural normalized target barrier failed")
    require(arb(by_mode[39_852]["carrier_scaled_physical_defect_ball"]) > 130 * TARGET_PHYSICAL, "structural physical target barrier failed")

    artifact = {
        "kind": STEM,
        "status": "bare_full_saddle_termwise_splice_rigorously_rejected_at_both_overlap_edges",
        "passed": True,
        "exact_identity": {
            "grouped_current": "G_m=m*sqrt(2)*x_m^(-3/2) I0+x_m^(-1) I1",
            "moments": "I0=F(q_B)-F(q_A), I1=[exp(i*pi*q_B^2/2)-exp(i*pi*q_A^2/2)]/(i*pi)",
            "full_line_current": "G_m_full=m*sqrt(2)*x_m^(-3/2)(1+i)",
            "ratio": "R_m=I0/(1+i)+sqrt(x_m)I1/[m*sqrt(2)(1+i)]",
        },
        "certificate": {"height": T, "alpha_endpoints": [A, B], "modes": list(MODES), "rows": rows},
        "decision": {
            "bare_full_line_carrier_is_valid_termwise_overlap_object": False,
            "incomplete_fresnel_current_must_be_retained_or_signed_aggregated": True,
            "complete_x_integrated_mode_error_lower_bound_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {name: {"path": relative(path), "sha256": file_hash(path)} for name, path in DEPENDENCIES.items()},
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {"elapsed_seconds": round(time.time() - started, 3), "workers": 1, "flint_threads": 1, "process_priority": priority},
        "next_obligation": "Define the ordinary chart with the incomplete Fresnel plus endpoint current retained, derive its exact global Morse amplitude, and bound amplitude variation after summing the signed 158-mode overlap; do not compare each interface mode directly with the bare full-line carrier.",
        "proof_boundary": "Exact frozen-saddle grouped alpha-current barrier at four saved-height modes only. No lower bound for the complete x-integrated mode error, uniform ordinary remainder, complete T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }

    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "certified bare full-saddle interface barrier: "
        f"m=39695 > {arb(by_mode[39_695]['normalized_target_multiple_ball']).lower()} targets, "
        f"m=39852 > {arb(by_mode[39_852]['normalized_target_multiple_ball']).lower()} targets",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
