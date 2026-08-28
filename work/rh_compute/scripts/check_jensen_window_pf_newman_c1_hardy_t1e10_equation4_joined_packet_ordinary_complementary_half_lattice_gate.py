#!/usr/bin/env python3
"""Check the ordinary complementary half-lattice reduction and census."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_ordinary_complementary_half_lattice_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = ROOT / "work" / "rh_compute" / "scripts" / f"{STEM}.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def set_low_priority() -> None:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            process.cpu_affinity([process.cpu_affinity()[0]])
        else:
            process.nice(10)
            process.cpu_affinity([process.cpu_affinity()[0]])
    except Exception:
        pass


def load_builder():
    spec = importlib.util.spec_from_file_location("ordinary_complement_builder", BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load builder")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def altered_poisson_witness() -> tuple[mp.mpf, mp.mpf]:
    """Check finite roster plus explicit complements against a wide full lattice."""

    mp.mp.dps = 70
    t = mp.mpf("7.25")
    mode = 4
    q_first = mp.mpf("1.5")
    q_last = mp.mpf("6.5")

    def coefficient(q: mp.mpf) -> mp.mpc:
        relative = q - mode

        def integrand(u: mp.mpf) -> mp.mpc:
            return (
                (1 + u / mode) ** (-mp.mpf("0.5") - 1j * t)
                * mp.exp(-1j * mp.pi * u**2 + 2j * mp.pi * relative * u)
            )

        return mp.quad(integrand, [-mp.mpf("0.5"), 0, mp.mpf("0.5")])

    finite = mp.fsum(coefficient(mp.mpf(n) + mp.mpf("0.5")) for n in range(1, 7))
    cutoff_a, cutoff_b = 80, 120
    full_a = mp.fsum(coefficient(mp.mpf(n) + mp.mpf("0.5")) for n in range(-cutoff_a, cutoff_a + 1))
    full_b = mp.fsum(coefficient(mp.mpf(n) + mp.mpf("0.5")) for n in range(-cutoff_b, cutoff_b + 1))
    lower = mp.fsum(coefficient(mp.mpf(n) + mp.mpf("0.5")) for n in range(-cutoff_b, 1))
    upper = mp.fsum(coefficient(mp.mpf(n) + mp.mpf("0.5")) for n in range(7, cutoff_b + 1))
    complement_residual = abs((finite + lower + upper) - full_b)
    require(q_first == mp.mpf("1.5") and q_last == mp.mpf("6.5"), "altered roster drift")
    return max(abs(full_b - 1), abs(full_b - full_a)), complement_residual


def main() -> int:
    set_low_priority()
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True and NOTE.is_file(), "saved gate or note missing")
    for dependency in artifact["dependencies"].values():
        path = ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], "dependency hash drift")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")

    b = load_builder()
    b.ctx.prec = 448
    b.ctx.threads = 1
    replay = b.build_certificate()
    saved = artifact["certificate"]
    require(replay["upper_stationary_label_count"] == 230, "stationary census drift")
    for key in (
        "upper_stationary_span_ball",
        "upper_first_nonstationary_gap_ball",
        "lower_nearest_gap_ball",
        "ordinary_Gamma_defect_packet_log10_upper_ball",
    ):
        require(
            b.arb(replay[key]["ball"]).overlaps(b.arb(saved[key]["ball"])),
            f"precision replay misses saved {key}",
        )

    poisson_error, complement_residual = altered_poisson_witness()
    require(poisson_error < mp.mpf("0.01"), "altered full half-lattice has not converged to one")
    require(complement_residual < mp.mpf("1e-60"), "altered finite/complement partition failed")

    print(
        "checked complementary half-lattice reduction; "
        f"altered Poisson error {mp.nstr(poisson_error, 8)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
