#!/usr/bin/env python3
"""Independently check the positive-real-tail Hardy-operator gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
for directory in (SCRIPT_DIR, VENDOR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))

from flint import arb, ctx
import mpmath as mp

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_positive_real_tail_K_first_subcell_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finite_difference(u: mp.mpf, q: mp.mpf, count: int) -> mp.mpf:
    return mp.exp(-q * u) * (-mp.expm1(-count * u)) / (-mp.expm1(-u))


def altered_direct_fixture() -> tuple[mp.mpf, mp.mpf]:
    """Differentiate an altered tail packet independently of the closed bound."""

    def run(dps: int) -> mp.mpc:
        with mp.workdps(dps):
            t0 = mp.mpf("3.25")
            U = mp.mpf("3.5")
            R = 2 * mp.pi * U
            q = mp.mpf("2.5")
            count = 5
            hardy_shift = mp.mpc("-0.2", "0.7")
            cuts = (R, mp.mpf(26), mp.mpf(34), mp.mpf(46), mp.mpf(64), mp.mpf(90))

            def E(t: mp.mpf) -> mp.mpc:
                s = mp.mpc("0.5", t)
                return mp.exp(mp.pi * t / 2 + 0.25j * mp.pi) * (2 * mp.pi) ** (s - 1)

            def tail(t: mp.mpf) -> mp.mpc:
                s = mp.mpc("0.5", t)
                return mp.quad(
                    lambda u: u ** (-s)
                    * mp.exp(1j * u * u / (4 * mp.pi))
                    * finite_difference(u, q, count),
                    cuts,
                )

            def packet(t: mp.mpf) -> mp.mpc:
                return E(t) * tail(t)

            numerical_operator = mp.diff(packet, t0) + hardy_shift * packet(t0)
            s0 = mp.mpc("0.5", t0)
            weighted = E(t0) * mp.quad(
                lambda u: (
                    mp.pi / 2
                    - mp.mpf("0.2")
                    + 1j * (mp.mpf("0.7") - mp.log(u / (2 * mp.pi)))
                )
                * u ** (-s0)
                * mp.exp(1j * u * u / (4 * mp.pi))
                * finite_difference(u, q, count),
                cuts,
            )
            return +(numerical_operator - weighted)

    discrepancy_80 = run(80)
    discrepancy_110 = run(110)
    drift = abs(discrepancy_110 - discrepancy_80)
    discrepancy = abs(discrepancy_110)
    require(discrepancy < mp.mpf("1e-55"), "altered tail operator identity failed")
    require(drift < mp.mpf("1e-70"), "altered tail quadrature drift")
    return discrepancy, drift


def main() -> int:
    priority = gate.joined.set_low_priority()
    require(priority == "below_normal_one_cpu", f"resource cap unavailable: {priority}")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact did not pass")
    require(
        artifact["status"]
        == "positive_real_tail_Hardy_operator_first_nonzero_height_subcell_certified",
        "status drift",
    )
    for row in artifact["dependencies"].values():
        path = ROOT / row["path"]
        require(path.is_file() and file_hash(path) == row["sha256"], "dependency hash drift")
    for relative, expected in artifact["sources"].items():
        path = ROOT / relative
        require(path.is_file() and file_hash(path) == expected, "source hash drift")
    decision = artifact["decision"]
    require(decision["positive_real_tail_K_interval_proved"] is True, "tail claim lost")
    require(
        decision["positive_real_tail_Hardy_contribution_negligible"] is True,
        "tail scale claim lost",
    )
    for key in (
        "complete_K_T_interval_proved",
        "wider_Q_K_minus_T_sign_interval_proved",
        "maximal_event_cell_sign_proved",
        "rh_implication",
    ):
        require(decision[key] is False, f"overclaim at {key}")

    ctx.prec = 512
    ctx.threads = 1
    production_log = arb(artifact["certificate"]["log10_modulus_bound"]["ball"])
    halves = (
        arb(arb(gate.HEIGHT) - arb("0.00005"), arb("0.00005")),
        arb(arb(gate.HEIGHT) + arb("0.00005"), arb("0.00005")),
    )
    replay_logs: list[arb] = []
    for t in halves:
        replay = gate.tail_certificate(t, use_duplication_h=True)
        value = arb(replay["log10_modulus_bound"]["ball"])
        require(value.upper() < arb("-1000000000"), "half-box tail scale lost")
        require(value.upper() <= production_log.upper(), "half-box bound exceeds production")
        replay_logs.append(value)

    discrepancy, drift = altered_direct_fixture()
    note = " ".join(gate.NOTE.read_text(encoding="utf-8").split()).lower()
    for token in (
        "differentiation under the absolutely convergent positive-real tail",
        "log v<=v-1",
        "no fitted geometric constant",
        "exact tiny target correction and final assembly",
        "no complete `k_t`",
    ):
        require(token in note, f"note token missing: {token}")
    print(
        "independently checked positive-real-tail K first height subcell; "
        f"half_log10={[value.upper().str(18, more=True) for value in replay_logs]}, "
        f"altered_discrepancy={mp.nstr(discrepancy, 8)}, drift={mp.nstr(drift, 8)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
