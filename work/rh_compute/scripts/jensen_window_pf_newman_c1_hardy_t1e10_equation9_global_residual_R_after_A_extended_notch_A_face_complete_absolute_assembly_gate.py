#!/usr/bin/env python3
"""Assemble the complete translated 42-mode A-face absolute bound."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
from flint import arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_complete_absolute_assembly_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "pre_endpoint_core": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_six_current_core_absolute_peano_interval_gate.json",
    "endpoint_layer": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_endpoint_eight_cell_peano_interval_gate.json",
}
A = 159_577
T = 10_000_000_000
ASSEMBLY_THRESHOLD = "0.00364"
R_AFTER_A_TARGET = "0.0368147039947"


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


def assemble(core: dict[str, Any], endpoint: dict[str, Any]) -> dict[str, str]:
    ctx.dps = 120
    core_bound = arb(core["interval_certificate"]["physical_complete_pre_endpoint_core_absolute_bound"])
    endpoint_bound = arb(endpoint["interval_certificate"]["physical_complete_endpoint_bound"])
    core_x16 = arb(core["interval_certificate"]["exact_x_16_ball"])
    endpoint_x16 = arb(endpoint["interval_certificate"]["exact_x_16_ball"])
    require(core_x16.overlaps(endpoint_x16), "x_16 handoff balls do not overlap")
    total = core_bound + endpoint_bound
    remaining = arb(R_AFTER_A_TARGET) - total
    require(total.upper() < arb(ASSEMBLY_THRESHOLD), "complete A-face assembly threshold failed")
    require(remaining.lower() > arb("0.03317"), "remaining R_after_A allowance drift")
    return {
        "pre_endpoint_core_bound": core_bound.str(80, more=True),
        "endpoint_layer_bound": endpoint_bound.str(80, more=True),
        "complete_translated_A_face_absolute_bound": total.str(80, more=True),
        "complete_A_face_threshold": ASSEMBLY_THRESHOLD,
        "saved_R_after_A_target": R_AFTER_A_TARGET,
        "remaining_R_after_A_allowance_ball": remaining.str(80, more=True),
        "remaining_R_after_A_allowance_lower_bound": remaining.lower().str(70, more=True),
        "common_x_16_handoff_ball": core_x16.str(80, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    cert = artifact["assembly_certificate"]
    return f"""# Complete translated A-face absolute assembly

Date: 2026-08-23

Status: rigorous saved-height local assembly; not a proof of `R_after_A` or RH

The exact handoff `x_16` partitions the translated 42-mode A-face channel as

```text
I_(A,42)=I_(A,42)^core+I_(A,42)^end,

core: 0<x<=x_16,       endpoint: x_16<=x<=1/2.      (AA1)
```

Both dependency gates use the same equation-(9) physical normalization and
their certified `x_16` balls overlap.  Their absolute bounds are

```text
|I_(A,42)^core| <= {cert['pre_endpoint_core_bound']},
|I_(A,42)^end|  <= {cert['endpoint_layer_bound']}.  (AA2)
```

Therefore the triangle inequality, with no oscillatory cancellation, gives

```text
|I_(A,42)| <= {cert['complete_translated_A_face_absolute_bound']}
             < {ASSEMBLY_THRESHOLD}.                (AA3)
```

Against the inherited joined target `R_after_A<{R_AFTER_A_TARGET}`, this
leaves a certified allowance greater than `0.03317` for the other post-A
channels.  That subtraction is budget bookkeeping only; it is not a proof
that those channels fit the allowance.

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Pi provenance: no new `pi` is introduced in (AA1)--(AA3).  Both summands
inherit their occurrences from the canonical Fresnel hierarchy and the same
equation-(9) normalization.

Proof boundary: (AA3) proves only the complete translated 42-mode A-face
channel at `t=10^10`.  The zero sector, lower edge, B-owned terms, negative
bulk, and every other post-A remainder remain outside this assembly.  No
joined `R_after_A`, `R_Dir`, `Q_K-T`, all-height theorem, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "an A-face dependency is not passed")
    core = dependencies["pre_endpoint_core"]
    endpoint = dependencies["endpoint_layer"]
    require(core["scope"]["height"] == T and endpoint["scope"]["height"] == T, "height mismatch")
    require(core["scope"]["A"] == A and endpoint["scope"]["A"] == A, "A mismatch")
    require(core["decision"]["complete_pre_endpoint_A_face_core_absolute_bound_below_1_point_2e_minus_5_proved"] is True, "core theorem drift")
    require(endpoint["decision"]["complete_A_face_endpoint_layer_absolute_bound_below_0_point_0037_proved"] is True, "endpoint theorem drift")
    require(core["decision"]["independent_replay_completed"] is True, "core replay missing")
    require(endpoint["decision"]["independent_replay_completed"] is True, "endpoint replay missing")
    require(CHECKER.is_file(), "independent checker missing")
    artifact = {
        "kind": STEM,
        "status": "rigorous_complete_translated_42_mode_A_face_absolute_bound_below_0_point_00364_at_t_1e10",
        "passed": True,
        "scope": {"height": T, "A": A, "x_region": "0<x<=1/2", "modes": [39_895, 39_936]},
        "exact_assembly": {
            "partition": "I_(A,42)=I_core+I_endpoint at the common x_16 handoff",
            "orientation": "E_42=-D_42 on both pieces; absolute values remove the common sign",
            "normalization": "both dependencies use 2*(pi/(32T))^(1/4) and the same equation-(9) weight",
        },
        "assembly_certificate": assemble(core, endpoint),
        "decision": {
            "common_x16_handoff_verified": True,
            "common_physical_normalization_verified": True,
            "complete_translated_42_mode_A_face_absolute_bound_below_0_point_00364_proved": True,
            "oscillatory_cancellation_used": False,
            "remaining_R_after_A_allowance_above_0_point_03317_bookkept": True,
            "other_post_A_channels_bounded": False,
            "R_after_A_bound_proved": False,
            "R_Dir_bound_proved": False,
            "QK_minus_T_bound_proved": False,
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
        "runtime": {"elapsed_seconds": round(time.time() - started, 3), "workers": 1},
        "next_obligation": "Return to the exact post-A ownership ledger and bound the channels outside this translated 42-mode A-face assembly against the remaining >0.03317 allowance. Do not count any A-owned term twice.",
        "proof_boundary": "Rigorous complete translated 42-mode A-face absolute bound below 0.00364 at t=10^10 only. No bound for the other post-A channels, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("assembled complete translated 42-mode A-face bound below 0.00364", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
