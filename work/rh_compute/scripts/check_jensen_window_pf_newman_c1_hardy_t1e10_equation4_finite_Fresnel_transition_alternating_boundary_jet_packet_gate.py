#!/usr/bin/env python3
"""Independently check the alternating boundary-jet transition packet."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_transition_alternating_boundary_jet_packet_gate"
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
        else:
            process.nice(10)
    except Exception:
        pass


def load_builder():
    spec = importlib.util.spec_from_file_location("transition_boundary_jet_builder", BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load builder module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    set_low_priority()
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True and NOTE.is_file(), "saved gate or note missing")
    require(
        artifact["decision"]["grouped_Fresnel_Gamma_A_transition_packet_enclosed"] is True,
        "transition decision drift",
    )
    for dependency in artifact["dependencies"].values():
        path = ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], "dependency hash drift")
    for source in artifact["sources"].values():
        path = ROOT / source["path"]
        require(path.is_file() and file_hash(path) == source["sha256"], "source hash drift")

    b = load_builder()
    b.ctx.dps = 125
    b.ctx.threads = 1
    replay = b.transition_source_certificate(
        b.arb(b.HEIGHT),
        b.SOURCE_ENDPOINT,
        b.SOURCE_COUNT,
        b.CELL_LEFT,
        b.CELL_RIGHT,
        14,
        256,
        48,
        36,
        "1e-42",
    )
    saved_source = b.complex_from_record(
        artifact["source_certificate"]["complete_transition_source_integral_ball"]
    )
    replay_source = b.complex_from_record(replay["complete_transition_source_integral_ball"])
    require(replay_source.overlaps(saved_source), "higher-order transition source replay misses saved ball")
    require(
        b.arb(replay["later_label_remainder_bound_ball"]).upper()
        < b.arb(artifact["source_certificate"]["later_label_remainder_bound_ball"]).upper(),
        "higher-order remainder did not tighten",
    )

    dependencies = {
        name: json.loads((ROOT / row["path"]).read_text(encoding="utf-8"))
        for name, row in artifact["dependencies"].items()
    }
    joined = b.joined_packet_certificate(replay, dependencies)
    saved_joined = b.complex_from_record(artifact["joined_packet"]["joined_transition_packet_ball"])
    replay_joined = b.complex_from_record(joined["joined_transition_packet_ball"])
    require(replay_joined.overlaps(saved_joined), "higher-order joined packet misses saved ball")

    # Altered finite-roster identity: direct integration versus the jet formula.
    pi = b.arb.pi()
    altered_t = pi * b.arb(31) ** 2 / 8
    altered = b.transition_source_certificate(
        altered_t,
        31,
        9,
        b.arb("6.5"),
        b.arb("8.5"),
        8,
        96,
        20,
        30,
        "1e-38",
    )
    p = b.centered_parameters(altered_t, 31)
    direct = b.acb(0)
    for label_index in range(9):
        def integrand(y, analytic, j=label_index):
            return b.centered_first_label_value(y, p, 30, analytic) * (
                b.acb(0, 1) * 2 * pi * j * y
            ).exp()

        direct += b.acb.integral(
            integrand,
            b.arb("6.5"),
            b.arb("8.5"),
            abs_tol=b.arb("1e-38"),
            rel_tol=b.arb("1e-38"),
            eval_limit=300_000,
            depth_limit=40,
        )
    altered_ball = b.complex_from_record(altered["complete_transition_source_integral_ball"])
    require(direct.overlaps(altered_ball), "altered direct finite roster misses jet collapse")

    print("independently checked alternating boundary-jet transition packet", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
