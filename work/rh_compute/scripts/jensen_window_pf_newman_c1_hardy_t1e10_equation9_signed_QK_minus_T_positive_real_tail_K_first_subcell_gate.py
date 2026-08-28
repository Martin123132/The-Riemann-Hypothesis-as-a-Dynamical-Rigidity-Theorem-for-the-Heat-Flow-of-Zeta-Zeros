#!/usr/bin/env python3
"""Certify the positive-real-tail Hardy operator on the first height subcell."""

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

from flint import arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_height_derivative_identity_gate as joined
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_lower_finite_cell_K_first_subcell_scout as lower_k


HEIGHT = 10_000_000_000
RADIUS = "0.0001"
U = arb("39936.5")
Q_MINUS = arb("79788.5")
LABEL_COUNT = 2_481_423
STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "positive_real_tail_K_first_subcell_gate"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
CHECKER = SCRIPT_DIR / f"check_{STEM}.py"
DEPENDENCIES = {
    "quarter_disk": ROOT / "work/rh_compute/results" / (
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_"
        "quarter_disk_vertical_arc_tail_gate.json"
    ),
    "height_derivative": ROOT / "work/rh_compute/results" / (
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
        "joined_height_derivative_identity_gate.json"
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


def real_record(value: arb, digits: int = 80) -> dict[str, str]:
    return {
        "ball": value.str(digits, more=True),
        "lower": value.lower().str(55, more=True),
        "upper": value.upper().str(55, more=True),
    }


def tail_certificate(t: arb, *, use_duplication_h: bool) -> dict[str, Any]:
    """Return the elementary modulus bound for L_t[E_s T_U]."""

    pi = arb.pi()
    two_pi = 2 * pi
    R = two_pi * U
    theta, theta_prime = joined.theta_and_derivative(t)
    direct_h = joined.direct_log_H_derivative(t)
    duplication_h = joined.duplication_log_H_derivative(t)
    require(direct_h.overlaps(duplication_h), "H logarithmic derivatives miss")
    h = duplication_h if use_duplication_h else direct_h
    H, H_transport = lower_k.stable_H_box(t)
    require(H.lower() > 0, "stable H lost positivity")

    real_weight = abs(pi / 2 - h).upper()
    endpoint_phase_weight = abs(theta_prime - U.log()).upper()
    endpoint_weight = (real_weight + endpoint_phase_weight).upper()
    first_log_moment = (1 / (R * Q_MINUS)).upper()
    weighted_moment = (
        endpoint_weight / Q_MINUS + 1 / (R * Q_MINUS**2)
    ).upper()
    denominator = 1 - (-R).exp()
    require(denominator.lower() > 0, "tail geometric denominator lost positivity")

    log_prefactor = pi * t / 2 - two_pi.log() / 2
    log_modulus_bound = (
        log_prefactor
        - R.log() / 2
        - Q_MINUS * R
        - denominator.log()
        + weighted_moment.log()
    ).upper()
    log10_modulus_bound = (log_modulus_bound / arb(10).log()).upper()
    modulus_bound = log_modulus_bound.exp().upper()
    projection_log_bound = (
        log_modulus_bound + (2 / H.lower()).log()
    ).upper()
    projection_bound = projection_log_bound.exp().upper()
    require(log10_modulus_bound < arb("-1000000000"), "tail derivative is not negligible")
    return {
        "height_ball": t.str(70, more=True),
        "U": U.str(40, more=True),
        "R": R.str(70, more=True),
        "q_minus": Q_MINUS.str(40, more=True),
        "label_count": LABEL_COUNT,
        "H_log_derivative_formula": "duplication" if use_duplication_h else "direct",
        "theta": real_record(theta),
        "theta_prime": real_record(theta_prime),
        "H_log_prime": real_record(h),
        "H": real_record(H),
        "H_transport": H_transport,
        "real_weight_upper": real_record(real_weight),
        "endpoint_phase_weight_upper": real_record(endpoint_phase_weight),
        "endpoint_weight_upper": real_record(endpoint_weight),
        "first_log_moment_factor_upper": real_record(first_log_moment),
        "weighted_exponential_moment_upper": real_record(weighted_moment),
        "log_modulus_bound": real_record(log_modulus_bound),
        "log10_modulus_bound": real_record(log10_modulus_bound),
        "modulus_bound": real_record(modulus_bound),
        "Hardy_projection_over_H_absolute_bound": real_record(projection_bound),
        "Hardy_projection_over_H_log10_absolute_bound": real_record(
            projection_log_bound / arb(10).log()
        ),
        "identities": {
            "packet": "P_tail(t)=E_s(t)*integral_R^infinity F_s(u,t)du",
            "operator": (
                "L_t[P_tail]=E_s*integral_R^infinity "
                "[pi/2-H'/H+i(theta'-log(u/(2*pi)))]F_s(u,t)du"
            ),
            "finite_difference_bound": (
                "0<Delta_q(u)<=exp(-q_minus*u)/(1-exp(-R)), u>=R"
            ),
            "log_weight_bound": (
                "|theta'-log(u/(2*pi))|<=|theta'-log(U)|+log(u/R)"
            ),
            "moment_bound": (
                "log(u/R)<= (u-R)/R and u^(-1/2)<=R^(-1/2)"
            ),
        },
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Positive-real-tail K on the first height subcell

Date: 2026-08-28

Status: rigorous positive-real-tail Hardy-operator bound; complete `K_T`
is not yet assembled

For `s=1/2+it`, `R=2*pi*U`, and

```text
F_s(u,t)=u^(-s) exp(i*u^2/(4*pi)) Delta_q(u),
P_tail(t)=E_s(t) integral_R^infinity F_s(u,t)du,
E_s=exp(pi*t/2+i*pi/4)(2*pi)^(s-1),
```

the fixed-roster event atlas keeps `R`, `q_-`, and the finite label count
constant on `I_1=[10^10-10^-4,10^10+10^-4]`.  Therefore differentiation
under the absolutely convergent positive-real tail gives exactly

```text
(d_t+i theta'-H'/H)P_tail
 =E_s integral_R^infinity
  [pi/2-H'/H+i(theta'-log(u/(2*pi)))]F_s(u,t)du.     (PT1)
```

For real `u>=R`, the finite geometric difference is positive and

```text
0<Delta_q(u)<=exp(-q_- u)/(1-exp(-R)).              (PT2)
```

Writing `u=Rv`, the logarithmic weight is bounded without discarding its
endpoint cancellation:

```text
|theta'-log(u/(2*pi))|
 <=|theta'-log U|+log v,
log v<=v-1,                u^(-1/2)<=R^(-1/2).      (PT3)
```

Equations (PT1)--(PT3) yield the elementary closed bound

```text
|L_t[P_tail]|
 <= |E_s| R^(-1/2) exp(-q_-R)/(1-exp(-R))
    *[A/q_-+1/(R q_-^2)],
A=|pi/2-H'/H|+|theta'-log U|.                      (PT4)
```

Arb uses the direct and duplication formulas for `H'/H` as an overlap guard
and stable positive `H` transport.  It obtains

```text
log10 |L_t[P_tail]| <= {c['log10_modulus_bound']['upper']},
log10 |Hardy_t[L_t[P_tail]]/H|
                      <= {c['Hardy_projection_over_H_log10_absolute_bound']['upper']}.
                                                               (PT5)
```

The checker raises precision, splits `I_1` into two changed half-boxes, uses
the duplication formula for `H'/H`, and reconstructs (PT4).  Its altered
five-label packet independently differentiates `E_s T_R` numerically and
agrees with the weighted integral, catching the `pi/2`, `log(2*pi)`, and
`-log u` signs.

Pi provenance: `pi` comes only from the Mellin--Fresnel kernel, the exact
quarter-disk scale `R=2*pi*U`, the Gamma normalization in `E_s`, and the
Riemann--Siegel phase.  No fitted geometric constant supplies `pi`.

Proof boundary: (PT1)--(PT5) certify only the positive-real-tail contribution
on `I_1`.  The exact tiny target correction and final assembly with the
transition--upper-arc and lower-plus-ordinary packets remain open.  No
complete `K_T`, wider `Q_K-T` sign interval, event-cell theorem, wall handoff,
all-height theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion
is proved.
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
        dependencies["quarter_disk"]["artifact"]["production"]["Y"] == "39936.5",
        "quarter-disk U drift",
    )
    require(
        dependencies["event_atlas"]["artifact"]["decision"]
        ["maximal_same_roster_cell_certified"]
        is True,
        "event-atlas roster drift",
    )
    t = arb(arb(HEIGHT), arb(RADIUS))
    certificate = tail_certificate(t, use_duplication_h=False)
    old_tail_log10 = arb(
        dependencies["quarter_disk"]["artifact"]["arb_certificate"]
        ["error_budgets"]["positive_real_tail_log10_upper"]["ball"]
    )
    new_tail_log10 = arb(certificate["log10_modulus_bound"]["ball"])
    require(new_tail_log10.upper() - old_tail_log10.upper() < 1, "derivative tail factor drift")
    artifact = {
        "kind": "rh_c1_hardy_positive_real_tail_K_first_height_subcell_gate",
        "date": "2026-08-28",
        "status": "positive_real_tail_Hardy_operator_first_nonzero_height_subcell_certified",
        "passed": True,
        "resource_mode": priority,
        "certificate": certificate,
        "comparison_to_value_tail": {
            "saved_value_tail_log10_upper": real_record(old_tail_log10),
            "derivative_minus_value_log10_upper": real_record(
                (new_tail_log10.upper() - old_tail_log10.upper()).upper()
            ),
        },
        "decision": {
            "positive_real_tail_K_interval_proved": True,
            "positive_real_tail_Hardy_contribution_negligible": True,
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
            "Certify the exact tiny target-correction Hardy operator, then assemble it with "
            "the positive tail, transition--upper-arc, and lower-plus-ordinary packets into "
            "a complete independently checked K_T enclosure on I_1."
        ),
        "proof_boundary": (
            "Positive-real-tail Hardy-operator modulus only on |t-10^10|<=10^-4. "
            "No tiny target correction, complete K_T, wider sign interval, full event-cell "
            "theorem, wall handoff, all-height theorem, Lambda<=0, PF-infinity, RH, or "
            "prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "certified positive-real-tail K first height subcell; "
        f"log10_projection_bound={certificate['Hardy_projection_over_H_log10_absolute_bound']['upper']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
