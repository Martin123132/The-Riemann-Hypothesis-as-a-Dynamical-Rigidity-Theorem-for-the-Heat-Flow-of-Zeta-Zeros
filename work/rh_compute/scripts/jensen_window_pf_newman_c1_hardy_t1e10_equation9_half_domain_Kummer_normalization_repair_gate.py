#!/usr/bin/env python3
"""Repair the half-domain/full-Kummer normalization in the physical ledger."""

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

from mpmath import mp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_domain_Kummer_normalization_repair_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "half_reflection": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_kummer_reflection_branch_reduction_gate.json",
    "source_reassembly": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_folded_abel_source_roster_reassembly_gate.json",
    "sparse_pilot": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_source_owned_carrier_sparse_pilot_gate.json",
    "superseded_ownership_ledger": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_physical_transform_ownership_ledger_gate.json",
}


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


def numerical_reflection_check(t_text: str = "7.25", alpha_value: int = 5, dps: int = 80) -> dict[str, str]:
    mp.dps = dps
    t = mp.mpf(t_text)
    alpha = mp.mpf(alpha_value)
    a = mp.mpf(3) / 4 - mp.j * t / 2
    b = mp.mpf(3) / 4 + mp.j * t / 2
    z = mp.j * mp.pi * alpha**2 / 4
    integrand = lambda x: alpha * x ** (a - 1) * (1 - x) ** (b - 1) * mp.exp(z * x)
    half = mp.quad(integrand, [0, mp.mpf("0.125"), mp.mpf("0.5")])
    full = mp.quad(
        integrand,
        [0, mp.mpf("0.125"), mp.mpf("0.5"), mp.mpf("0.75"), mp.mpf("0.875"), 1],
    )
    closed = alpha * mp.beta(a, b) * mp.hyp1f1(a, mp.mpf("1.5"), z)
    rotation = mp.exp(-mp.j * mp.pi / 8)
    full_error = abs(full - closed) / abs(closed)
    reflection_error = abs(rotation * full - 2 * mp.re(rotation * half)) / abs(full)
    c_t = (mp.pi / (32 * t)) ** mp.mpf("0.25")
    corrected = c_t * mp.re(rotation * closed)
    half_projection = 2 * c_t * mp.re(rotation * half)
    projection_error = abs(corrected - half_projection) / abs(corrected)
    old_literal = 2 * c_t * mp.re(rotation * closed)
    old_ratio = old_literal / corrected
    require(full_error < mp.mpf("1e-55"), "Euler-Kummer full-integral check failed")
    require(reflection_error < mp.mpf("1e-55"), "odd-label half-reflection check failed")
    require(projection_error < mp.mpf("1e-55"), "corrected physical projection check failed")
    require(abs(old_ratio - 2) < mp.mpf("1e-70"), "factor-two witness drift")
    return {
        "t": t_text,
        "alpha": str(alpha_value),
        "precision_decimal_digits": str(dps),
        "Euler_Kummer_relative_discrepancy": mp.nstr(full_error, 40),
        "half_reflection_relative_discrepancy": mp.nstr(reflection_error, 40),
        "corrected_projection_relative_discrepancy": mp.nstr(projection_error, 40),
        "superseded_literal_to_corrected_ratio": mp.nstr(old_ratio, 40),
    }


def render_note(artifact: dict[str, Any]) -> str:
    check = artifact["numerical_surrogate"]
    return f"""# Half-domain Kummer normalization repair

Date: 2026-08-27

Status: exact normalization repair certified; no numerical enclosure of `Q_K`

The source roster and every joined carrier in the half-Kummer branch live on
`0 <= x <= 1/2`.  The literal full-domain functional printed in Section
11.478 retained the half-domain factor `2 Re` while changing the integration
limit to `1`; read literally, that doubles every full Kummer label.

For an odd label `alpha`, odd-square parity gives

```text
K_alpha(1-x)=exp(i*pi/4) conjugate(K_alpha(x)),

exp(-i*pi/8) integral_0^1 K_alpha(x) dx
 =2 Re[exp(-i*pi/8) integral_0^(1/2) K_alpha(x) dx].   (NR1)
```

Hence the corrected physical functional and full-Kummer label are

```text
c_t=(pi/(32t))^(1/4),

P_t[F]=2c_t Re[exp(-i*pi/8) integral_0^(1/2) W_t(x)F(x) dx],

K_t(alpha)=c_t Re[exp(-i*pi/8) alpha B(a,b)
                   1F1(a;3/2;i*pi*alpha^2/4)].        (NR2)
```

The superseded line had `2c_t` in the last formula.  Its exact ratio to the
corrected label is `2`.  An independent altered-height surrogate verifies the
Euler-Kummer identity, odd-label reflection, and the corrected projection.
The builder discrepancies are respectively
`{check['Euler_Kummer_relative_discrepancy']}`,
`{check['half_reflection_relative_discrepancy']}`, and
`{check['corrected_projection_relative_discrepancy']}`.

This repair does not rescale the A, B, Gamma, or source-carrier certificates:
their executable scopes already use the half interval.  It corrects the
source-transform description and the direct full-Kummer formula.  Therefore
the identities

```text
J_Z=Q_K-G-A_transition,
R_nonA=J_Z-E_Btr,win-E_outer-I_(A,42)
```

and their scalar sufficient corridors remain unchanged in the corrected
half-domain normalization.

Proof boundary: exact domain/factor repair and deterministic surrogate checks
only.  No value or interval for `Q_K`, `J_Z`, `Q_K-T`, or the non-A residual,
no all-height theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion
is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "dependency failure")
    require(
        dependencies["half_reflection"]["decision"]["full_source_main_equals_twice_real_half_integral"]
        is True,
        "half-reflection dependency drift",
    )
    require(
        dependencies["source_reassembly"]["scope"]["half_x_interval"] == ["0", "1/2"],
        "source-reassembly domain drift",
    )
    require(
        dependencies["sparse_pilot"]["scope"]["x_interval"] == ["0", "1/2"],
        "sparse-pilot domain drift",
    )
    old_source = dependencies["superseded_ownership_ledger"]["source_transform_certificate"]
    require("integral_0^1" in old_source["physical_functional"], "superseded domain witness missing")
    require(old_source["one_source_label"].startswith("K_t(alpha)=2("), "superseded factor-two witness missing")

    artifact = {
        "kind": STEM,
        "status": "exact_half_domain_physical_projection_and_full_Kummer_label_normalization_repaired",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "source_odd_roster": [159_577, 5_122_421],
            "source_label_count": 2_481_423,
            "physical_x_interval": ["0", "1/2"],
        },
        "superseded_literal": {
            "physical_functional": old_source["physical_functional"],
            "one_source_label": old_source["one_source_label"],
            "defect": "The half-domain factor 2 Re was applied to a full-domain integral, doubling the full Kummer label.",
        },
        "corrected_certificate": {
            "c_t": "c_t=(pi/(32t))^(1/4)",
            "physical_functional": "P_t[F]=2c_t Re[e^(-i*pi/8) integral_0^(1/2) W_t(x)F(x)dx]",
            "odd_label_reflection": "K_alpha(1-x)=e^(i*pi/4) conjugate(K_alpha(x))",
            "half_full_identity": "e^(-i*pi/8) integral_0^1 K_alpha=2 Re[e^(-i*pi/8) integral_0^(1/2) K_alpha]",
            "Euler_Kummer_identity": "integral_0^1 x^(a-1)(1-x)^(b-1)e^(zx)dx=B(a,b) 1F1(a;3/2;z)",
            "one_source_label": "K_t(alpha)=c_t Re[e^(-i*pi/8) alpha B(a,b) 1F1(a;3/2;i*pi*alpha^2/4)]",
            "complete_source_transform": "Q_K=sum_(alpha=A,A+2,...,B)K_t(alpha)=P_t[sum_(n=0)^L f_x(n)]",
            "joined_identity": "J_Z=Q_K-G-A_transition",
        },
        "numerical_surrogate": numerical_reflection_check(),
        "decision": {
            "literal_full_domain_2Re_functional_rejected": True,
            "half_domain_physical_functional_certified": True,
            "full_Kummer_label_prefactor_corrected_from_2c_to_c": True,
            "joined_J_Z_identity_retained_in_corrected_normalization": True,
            "scalar_corridor_arithmetic_changed": False,
            "Q_K_numerically_enclosed": False,
            "non_A_bound_proved": False,
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
            "process_priority": priority,
        },
        "next_obligation": "Re-express the corrected finite Kummer roster relative to the exact complementary Hardy upper component and isolate one signed truncation/continuation bridge. Do not use the superseded factor-two label or source-hybrid telemetry.",
        "proof_boundary": "Exact half-domain/full-Kummer normalization repair only. No numerical Q_K or J_Z enclosure, non-A bound, Q_K-T theorem, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified half-domain/full-Kummer normalization repair", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
