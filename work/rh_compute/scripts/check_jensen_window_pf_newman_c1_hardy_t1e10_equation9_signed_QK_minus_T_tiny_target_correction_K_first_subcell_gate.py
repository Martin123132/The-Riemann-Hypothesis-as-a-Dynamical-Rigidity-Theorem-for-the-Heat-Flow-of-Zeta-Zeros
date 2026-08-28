#!/usr/bin/env python3
"""Independently check the tiny target-correction Hardy operator."""

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

from flint import acb, arb, ctx
import mpmath as mp

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_tiny_target_correction_K_first_subcell_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def complex_from(record: dict) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def altered_product_rule_fixture() -> tuple[mp.mpf, mp.mpf]:
    """Compare the cancelled formula with a numerical product derivative."""

    def run(dps: int) -> mp.mpc:
        with mp.workdps(dps):
            t0 = mp.mpf("3.25")
            modes = tuple(range(2, 8))

            def theta(t: mp.mpf) -> mp.mpf:
                return mp.mpf("0.4") * t + mp.mpf("0.01") * t * t

            def theta_prime(t: mp.mpf) -> mp.mpf:
                return mp.mpf("0.4") + mp.mpf("0.02") * t

            def C(t: mp.mpf) -> mp.mpf:
                return 1 + mp.mpf("0.02") * mp.sin(t)

            def C_prime(t: mp.mpf) -> mp.mpf:
                return mp.mpf("0.02") * mp.cos(t)

            def H(t: mp.mpf) -> mp.mpf:
                return mp.exp(mp.mpf("0.03") * t)

            def h(t: mp.mpf) -> mp.mpf:
                return mp.mpf("0.03")

            def D(t: mp.mpf) -> mp.mpc:
                return mp.fsum(m ** (-mp.mpf("0.5") - 1j * t) for m in modes)

            def P(t: mp.mpf) -> mp.mpc:
                return (C(t) - H(t)) * D(t)

            numerical = mp.diff(P, t0) + (1j * theta_prime(t0) - h(t0)) * P(t0)
            D0 = D(t0)
            covariant_D = mp.fsum(
                1j
                * (theta_prime(t0) - mp.log(m))
                * m ** (-mp.mpf("0.5") - 1j * t0)
                for m in modes
            )
            compact = (
                (C(t0) - H(t0)) * covariant_D
                + (C_prime(t0) - h(t0) * C(t0)) * D0
            )
            return +(numerical - compact)

    first = run(85)
    second = run(115)
    discrepancy = abs(second)
    drift = abs(second - first)
    require(discrepancy < mp.mpf("1e-70"), "altered correction identity failed")
    require(drift < mp.mpf("1e-80"), "altered correction precision drift")
    return discrepancy, drift


def main() -> int:
    priority = gate.joined.set_low_priority()
    require(priority == "below_normal_one_cpu", f"resource cap unavailable: {priority}")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact did not pass")
    require(
        artifact["status"]
        == "tiny_target_correction_Hardy_operator_first_nonzero_height_subcell_certified",
        "status drift",
    )
    for row in artifact["dependencies"].values():
        path = ROOT / row["path"]
        require(path.is_file() and file_hash(path) == row["sha256"], "dependency hash drift")
    for relative, expected in artifact["sources"].items():
        path = ROOT / relative
        require(path.is_file() and file_hash(path) == expected, "source hash drift")
    decision = artifact["decision"]
    require(decision["tiny_target_correction_K_interval_proved"] is True, "claim lost")
    require(
        decision["tiny_target_correction_all_product_rule_forms_overlap"] is True,
        "product-rule overlap lost",
    )
    for key in (
        "complete_K_T_interval_proved",
        "wider_Q_K_minus_T_sign_interval_proved",
        "maximal_event_cell_sign_proved",
        "rh_implication",
    ):
        require(decision[key] is False, f"overclaim at {key}")

    ctx.prec = 448
    ctx.threads = 1
    replay = gate.correction_certificate(
        arb(arb(gate.HEIGHT), arb(gate.RADIUS)),
        reverse=True,
        use_duplication_h=True,
    )
    production_K = complex_from(
        artifact["certificate"]["compact_weighted_mode_sum_K_corr"]
    )
    replay_K = complex_from(replay["compact_weighted_mode_sum_K_corr"])
    production_projection = arb(
        artifact["certificate"]["Hardy_projection_over_H"]["ball"]
    )
    replay_projection = arb(replay["Hardy_projection_over_H"]["ball"])
    require(production_K.overlaps(replay_K), "reversed correction K misses production")
    require(
        production_projection.overlaps(replay_projection),
        "reversed correction projection misses production",
    )
    require(replay["all_three_forms_overlap"] is True, "replay product forms miss")
    require(abs(replay_projection).upper() < arb("1e-15"), "replay scale widened")
    discrepancy, drift = altered_product_rule_fixture()

    note = " ".join(gate.NOTE.read_text(encoding="utf-8").split()).lower()
    for token in (
        "the two apparent `h'd_t` terms cancel exactly",
        "never forms the separately widened",
        "no fitted geometric constant",
        "final assembly",
        "no wider `q_k-t` sign interval",
    ):
        require(token in note, f"note token missing: {token}")
    print(
        "independently checked tiny target-correction K first height subcell; "
        f"production={production_projection.str(18, more=True)}, "
        f"replay={replay_projection.str(18, more=True)}, "
        f"altered_discrepancy={mp.nstr(discrepancy, 8)}, drift={mp.nstr(drift, 8)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
