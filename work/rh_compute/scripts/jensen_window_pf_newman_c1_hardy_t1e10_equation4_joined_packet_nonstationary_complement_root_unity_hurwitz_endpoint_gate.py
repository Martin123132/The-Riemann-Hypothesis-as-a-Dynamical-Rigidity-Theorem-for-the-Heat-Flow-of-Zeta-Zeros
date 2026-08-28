#!/usr/bin/env python3
"""Certify root-of-unity Hurwitz sums for the complementary endpoint tails."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
for directory in (SCRIPT_DIR, VENDOR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_"
    "nonstationary_complement_root_unity_hurwitz_endpoint_gate"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

COLLAR_RESULT = ROOT / (
    "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_"
    "joined_packet_upper_complement_rational_240_label_collar_gate.json"
)
COMPLEMENT_RESULT = ROOT / (
    "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_"
    "joined_packet_ordinary_complementary_half_lattice_gate.json"
)
LOWER_RESULT = ROOT / (
    "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_"
    "joined_packet_lower_cell_reversed_roster_endpoint_gate.json"
)

PRECISION_BITS = 384
MAX_POWER = 16
HEIGHT = arb("10000000000")
L = arb("621.5")
B = arb("621.5625")
A = arb("39852.5")
Q_MINUS = arb("79788.5")
Q_PLUS = arb("2561211.5")
Q_REMOTE = Q_PLUS + 240
Q_LOW = Q_MINUS - 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            process.cpu_affinity([process.cpu_affinity()[0]])
            return "below_normal_one_cpu"
        process.nice(10)
        process.cpu_affinity([process.cpu_affinity()[0]])
        return "nice_10_one_cpu"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def arb_record(value: arb, digits: int = 72) -> dict[str, Any]:
    return {
        "ball": value.str(digits, more=True),
        "lower_float": float(value.lower()),
        "mid_float": float(value.mid()),
        "upper_float": float(value.upper()),
    }


def acb_record(value: acb, digits: int = 72) -> dict[str, str]:
    return {
        "real_ball": value.real.str(digits, more=True),
        "imag_ball": value.imag.str(digits, more=True),
        "absolute_ball": abs(value).str(digits, more=True),
    }


def contains_zero(value: acb) -> bool:
    return value.real.contains(0) and value.imag.contains(0)


def q_curve(y: arb, p: arb) -> arb:
    return y + p / y


def root_of_unity(order: int, numerator: int) -> acb:
    require(order > 1, "nontrivial root order required")
    require(0 < numerator < order, "root numerator outside reduced range")
    return (acb(0, 1) * 2 * arb.pi() * numerator / order).exp()


def endpoint_phase(q0: arb, y: arb) -> acb:
    return (acb(0, 1) * 2 * arb.pi() * q0 * y).exp()


def finite_root_sum(
    delta: arb,
    order: int,
    numerator: int,
    cycles: int,
    power: int,
) -> acb:
    """Return sum_(n=0)^(hC-1) omega^n/(n+delta)^power by residues."""

    require(delta.lower() > 0, "finite shifted denominator is not positive")
    require(power >= 1 and cycles >= 1, "invalid finite root-sum parameters")
    h = arb(order)
    omega = root_of_unity(order, numerator)
    value = acb(0)
    for residue in range(order):
        argument = (delta + residue) / h
        if power == 1:
            block = (argument + cycles).digamma() - argument.digamma()
        else:
            exponent = arb(power)
            block = exponent.zeta(argument) - exponent.zeta(argument + cycles)
        value += omega**residue * block
    return value / h**power


def infinite_root_sum(delta: arb, order: int, numerator: int, power: int) -> acb:
    """Return the root-of-unity Lerch sum, with the k=1 pole canceled."""

    require(delta.lower() > 0, "infinite shifted denominator is not positive")
    require(power >= 1, "positive denominator power required")
    h = arb(order)
    omega = root_of_unity(order, numerator)
    value = acb(0)
    for residue in range(order):
        argument = (delta + residue) / h
        if power == 1:
            value -= omega**residue * argument.digamma()
        else:
            value += omega**residue * arb(power).zeta(argument)
    return value / h**power


def direct_finite_root_sum(
    delta: arb,
    order: int,
    numerator: int,
    count: int,
    power: int,
) -> acb:
    require(count > 0, "positive finite count required")
    omega = root_of_unity(order, numerator)
    weight = acb(1)
    value = acb(0)
    for index in range(count):
        value += weight / (delta + index) ** power
        weight *= omega
    return value


def alternating_eta(delta: arb, power: int) -> arb:
    require(delta.lower() > 0 and power >= 1, "invalid alternating eta parameters")
    if power == 1:
        return (((delta + 1) / 2).digamma() - (delta / 2).digamma()) / 2
    exponent = arb(power)
    return (
        exponent.zeta(delta / 2) - exponent.zeta((delta + 1) / 2)
    ) / arb(2) ** power


def finite_endpoint_rows(
    y: arb,
    q0: arb,
    delta: arb,
    order: int,
    numerator: int,
    cycles: int,
) -> list[dict[str, Any]]:
    phase = endpoint_phase(q0, y)
    count = order * cycles
    rows: list[dict[str, Any]] = []
    for power in range(1, MAX_POWER + 1):
        formula = finite_root_sum(delta, order, numerator, cycles, power)
        direct = direct_finite_root_sum(delta, order, numerator, count, power)
        discrepancy = formula - direct
        require(formula.overlaps(direct), f"finite formula/direct miss at k={power}")
        require(contains_zero(discrepancy), f"finite discrepancy excludes zero at k={power}")
        rows.append(
            {
                "power": power,
                "unphased_formula_ball": acb_record(formula),
                "phase_weighted_sum_ball": acb_record(phase * formula),
                "direct_sum_ball": acb_record(direct),
                "formula_minus_direct_ball": acb_record(discrepancy),
            }
        )
    return rows


def infinite_endpoint_rows(
    y: arb,
    q0: arb,
    delta: arb,
    lower_orientation: bool,
) -> list[dict[str, Any]]:
    omega = root_of_unity(2, 1)
    phase = endpoint_phase(q0, y)
    rows: list[dict[str, Any]] = []
    for power in range(1, MAX_POWER + 1):
        formula = infinite_root_sum(delta, 2, 1, power)
        explicit = acb(alternating_eta(delta, power))
        shifted = infinite_root_sum(delta + 1, 2, 1, power)
        recurrence = formula - omega * shifted - acb(delta ** (-power))
        require(formula.overlaps(explicit), f"alternating formula miss at k={power}")
        require(contains_zero(formula - explicit), f"eta discrepancy excludes zero at k={power}")
        require(contains_zero(recurrence), f"Lerch recurrence excludes zero at k={power}")
        orientation = -1 if lower_orientation and power % 2 else 1
        rows.append(
            {
                "power": power,
                "unphased_formula_ball": acb_record(formula),
                "explicit_alternating_eta_ball": arb_record(explicit.real),
                "phase_weighted_sum_ball": acb_record(phase * orientation * formula),
                "lower_orientation_sign": orientation,
                "functional_recurrence_discrepancy_ball": acb_record(recurrence),
            }
        )
    return rows


def endpoint_record(
    *,
    y: arb,
    q0: arb,
    delta: arb,
    order: int,
    numerator: int,
    mode: str,
    cycles: int | None = None,
    lower_orientation: bool = False,
) -> dict[str, Any]:
    require(delta.lower() > 0, f"nonstationary gap lost at y={y}")
    if mode == "finite":
        require(cycles is not None, "finite endpoint cycles missing")
        rows = finite_endpoint_rows(y, q0, delta, order, numerator, cycles)
    elif mode == "infinite":
        rows = infinite_endpoint_rows(y, q0, delta, lower_orientation)
    else:
        raise RuntimeError(f"unknown endpoint mode: {mode}")
    return {
        "endpoint_y": y.str(40, more=True),
        "base_label_q0": q0.str(40, more=True),
        "positive_shift_delta_ball": arb_record(delta),
        "root_order": order,
        "root_numerator": numerator,
        "root_ball": acb_record(root_of_unity(order, numerator)),
        "common_label_phase_ball": acb_record(endpoint_phase(q0, y)),
        "mode": mode,
        "cycles": cycles,
        "rows": rows,
    }


def build_certificate() -> dict[str, Any]:
    pi = arb.pi()
    p = HEIGHT / (2 * pi)
    q_l = q_curve(L, p)
    q_b = q_curve(B, p)
    q_a = q_curve(A, p)

    rb_b_delta = Q_PLUS - q_b
    rb_a_delta = Q_PLUS - q_a
    rl_l_delta = Q_REMOTE - q_l
    rl_a_delta = Q_REMOTE - q_a
    lower_l_delta = q_l - Q_LOW
    lower_a_delta = q_a - Q_LOW
    for name, delta in (
        ("R_B at B", rb_b_delta),
        ("R_B at a", rb_a_delta),
        ("R_L at L", rl_l_delta),
        ("R_L at a", rl_a_delta),
        ("lower tail at L", lower_l_delta),
        ("lower tail at a", lower_a_delta),
    ):
        require(delta.lower() > 0, f"{name} acquired a stationary point")

    root16_sum = sum((root_of_unity(16, 9) ** r for r in range(16)), acb(0))
    root2_sum = acb(1) + root_of_unity(2, 1)
    require(contains_zero(root16_sum), "order-16 k=1 pole coefficient does not cancel")
    require(contains_zero(root2_sum), "alternating k=1 pole coefficient does not cancel")

    families = {
        "R_B_finite_240": {
            "definition": (
                "sum_(n=0)^239 exp(2*pi*i*(q_++n)*y)/(q_++n-Q(y))^k"
            ),
            "endpoint_B": endpoint_record(
                y=B,
                q0=Q_PLUS,
                delta=rb_b_delta,
                order=16,
                numerator=9,
                mode="finite",
                cycles=15,
            ),
            "endpoint_a": endpoint_record(
                y=A,
                q0=Q_PLUS,
                delta=rb_a_delta,
                order=2,
                numerator=1,
                mode="finite",
                cycles=120,
            ),
        },
        "R_L_infinite_upper": {
            "definition": (
                "sum_(n>=0) exp(2*pi*i*(q_remote+n)*y)/(q_remote+n-Q(y))^k"
            ),
            "endpoint_L": endpoint_record(
                y=L,
                q0=Q_REMOTE,
                delta=rl_l_delta,
                order=2,
                numerator=1,
                mode="infinite",
            ),
            "endpoint_a": endpoint_record(
                y=A,
                q0=Q_REMOTE,
                delta=rl_a_delta,
                order=2,
                numerator=1,
                mode="infinite",
            ),
        },
        "R_lower_infinite": {
            "definition": (
                "sum_(n>=0) exp(2*pi*i*(q_low-n)*y)/(q_low-n-Q(y))^k"
            ),
            "endpoint_L": endpoint_record(
                y=L,
                q0=Q_LOW,
                delta=lower_l_delta,
                order=2,
                numerator=1,
                mode="infinite",
                lower_orientation=True,
            ),
            "endpoint_a": endpoint_record(
                y=A,
                q0=Q_LOW,
                delta=lower_a_delta,
                order=2,
                numerator=1,
                mode="infinite",
                lower_orientation=True,
            ),
        },
    }

    return {
        "passed": True,
        "precision_bits": PRECISION_BITS,
        "maximum_denominator_power": MAX_POWER,
        "t": HEIGHT.str(50, more=True),
        "p_equals_t_over_2pi_ball": arb_record(p),
        "Q_of_L_ball": arb_record(q_l),
        "Q_of_B_ball": arb_record(q_b),
        "Q_of_a_ball": arb_record(q_a),
        "q_minus": Q_MINUS.str(30, more=True),
        "q_plus": Q_PLUS.str(30, more=True),
        "q_remote": Q_REMOTE.str(30, more=True),
        "q_low": Q_LOW.str(30, more=True),
        "k1_pole_cancellation": {
            "order_16_exact_identity": "sum_(r=0)^15 exp(2*pi*i*9*r/16)=0",
            "order_16_numeric_zero_ball": acb_record(root16_sum),
            "order_2_exact_identity": "1+exp(pi*i)=0",
            "order_2_numeric_zero_ball": acb_record(root2_sum),
            "meaning": (
                "the common Hurwitz zeta pole at k=1 has coefficient zero; "
                "the finite value is the displayed weighted digamma sum"
            ),
        },
        "finite_formula": (
            "F_fin=h^(-k)sum_r omega^r[zeta(k,(delta+r)/h)-"
            "zeta(k,(delta+r)/h+C)] for k>1, with the zeta difference "
            "replaced by psi((delta+r)/h+C)-psi((delta+r)/h) for k=1"
        ),
        "infinite_formula": (
            "F_inf=h^(-k)sum_r omega^r zeta(k,(delta+r)/h) for k>1; "
            "F_inf=-h^(-1)sum_r omega^r psi((delta+r)/h) for k=1"
        ),
        "functional_recurrence": (
            "F_inf(k,delta,omega)-omega*F_inf(k,delta+1,omega)=delta^(-k)"
        ),
        "families": families,
        "repeated_integration_scope": (
            "all phase-weighted denominator sums for k=1..16 are available; "
            "this supports endpoint terms through seven derivatives and the "
            "coefficient powers in an eight-round recurrence, but no integrated "
            "h_8 remainder is bounded here"
        ),
    }


def note_text(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    families = c["families"]
    rb_b = families["R_B_finite_240"]["endpoint_B"]
    rl_l = families["R_L_infinite_upper"]["endpoint_L"]
    lower_a = families["R_lower_infinite"]["endpoint_a"]
    return f"""# Root-of-unity Hurwitz endpoint algebra for the nonstationary complements

Date: 2026-08-27

Status: every phase-weighted endpoint denominator sum through power
`k={MAX_POWER}` is certified for `R_B`, `R_L`, and the lower complementary
tail; no integrated remainder or complete tail enclosure is claimed.

Put

```text
p=t/(2*pi),              Q(y)=y+p/y,
D_q(y)=q-Q(y),           phi_q'(y)=2*pi*D_q(y).     (EH1)
```

For a nontrivial root of unity `omega^h=1` and `M=hC`, define

```text
F_fin(k;delta,omega,M)=sum_(n=0)^(M-1) omega^n/(n+delta)^k.
```

Grouping by residue class gives, for `k>1`,

```text
F_fin=h^(-k) sum_(r=0)^(h-1) omega^r
 [zeta(k,(delta+r)/h)-zeta(k,(delta+r)/h+C)].       (EH2)
```

At `k=1` the bracket in (EH2) is replaced by
`psi((delta+r)/h+C)-psi((delta+r)/h)`.  Direct Arb summation of all 240
terms independently overlaps (EH2) at both `B` and `a` for every
`1<=k<={MAX_POWER}`.

For the infinite root-of-unity sum,

```text
F_inf=h^(-k)sum_r omega^r zeta(k,(delta+r)/h),       k>1,
F_inf=-h^(-1)sum_r omega^r psi((delta+r)/h),         k=1.       (EH3)
```

The apparent `k=1` Hurwitz poles cancel exactly because

```text
sum_(r=0)^15 exp(2*pi*i*9*r/16)=0,
1+exp(pi*i)=0.                                      (EH4)
```

Thus `pi` here comes from the inherited Fresnel/Fourier phase and
`p=t/(2*pi)`.  The order 16 is forced by the exact fractional endpoint
`B=621+9/16`; it is not a fitted circle constant.

The production shifts remain strictly positive.  Representative certified
values are

```text
delta_RB(B)    = {rb_b['positive_shift_delta_ball']['ball']},
delta_RL(L)    = {rl_l['positive_shift_delta_ball']['ball']},
delta_lower(a) = {lower_a['positive_shift_delta_ball']['ball']}.
```

At the half-integer endpoints the infinite sums are alternating and satisfy

```text
F_inf(k,delta,-1)+F_inf(k,delta+1,-1)=delta^(-k).   (EH5)
```

Arb interval checks of (EH5), and independent agreement with the explicit
alternating eta/digamma form, pass for every power through {MAX_POWER}.  The
lower tail additionally carries the exact orientation factor `(-1)^k` because
`D=-(delta+n)`.

Starting from `h_0=y^(-1/2)` and

```text
h_(n+1)=d/dy[h_n/D],                                (EH6)
```

the certified table supplies all endpoint denominator sums needed through
the first eight integration rounds.  It does not bound the integral containing
`h_8`; that is the next quantitative obligation.

Proof boundary: exact root-of-unity/Hurwitz endpoint identities, explicit
`k=1` pole cancellation, production positive-gap checks, and rigorous endpoint
sum balls through power 16 only.  No integrated remainder, `R_B`, `R_L`, lower
complementary tail, complete ordinary or joined packet, `J_Z`, `D_K`, non-A,
all-height, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    resource_mode = set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = PRECISION_BITS
    ctx.threads = 1

    dependencies = {
        "closed_240_collar": COLLAR_RESULT,
        "complementary_half_lattice": COMPLEMENT_RESULT,
        "lower_cell": LOWER_RESULT,
    }
    for path in (*dependencies.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    for name, path in dependencies.items():
        require(load_json(path).get("passed") is True, f"dependency not passed: {name}")

    collar = load_json(COLLAR_RESULT)
    require(collar["certificate"]["collar_label_count"] == 240, "collar count drift")
    require(collar["certificate"]["root_of_unity_order_at_split"] == 16, "root order drift")
    require(collar["decision"]["upper_remainders_enclosed"] is False, "tail boundary drift")

    certificate = build_certificate()
    artifact = {
        "kind": "rh_c1_hardy_equation4_nonstationary_complement_root_unity_hurwitz_endpoint_gate",
        "date": "2026-08-27",
        "status": "production_nonstationary_complement_endpoint_sums_through_power_16_certified",
        "passed": True,
        "resource_mode": resource_mode,
        "certificate": certificate,
        "decision": {
            "production_endpoint_denominator_sums_certified": True,
            "k1_hurwitz_poles_cancel_exactly": True,
            "finite_240_formula_matches_direct_sum": True,
            "infinite_alternating_recurrences_certified": True,
            "integrated_remainders_bounded": False,
            "R_B_enclosed": False,
            "R_L_enclosed": False,
            "lower_complementary_tail_enclosed": False,
            "complete_joined_packet_enclosed": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in dependencies.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": (
            "Generate the exact h_n coefficient recurrence for the three tail orientations, "
            "combine these endpoint sums through a cancellation-preserving eight-round "
            "integration-by-parts assembly, and prove interval majorants for the three "
            "integrated h_8 remainders before enclosing any tail."
        ),
        "proof_boundary": (
            "Rigorous endpoint algebra through denominator power 16 only. No integrated "
            "remainder, complementary tail, complete ordinary or joined packet, J_Z, D_K, "
            "non-A, all-height, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved."
        ),
    }
    atomic_write(RESULT, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE, note_text(artifact))
    print("certified root-of-unity Hurwitz endpoint sums", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
