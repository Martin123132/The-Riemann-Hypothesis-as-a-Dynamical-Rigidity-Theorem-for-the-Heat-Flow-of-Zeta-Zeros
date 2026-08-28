#!/usr/bin/env python3
"""Independently check the exact Q_K-to-Hardy-upper bridge target."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

from flint import acb, arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_QK_exact_hardy_truncation_bridge_target_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def altered_direct_kummer_guard() -> str:
    ctx.dps = 65
    ctx.threads = 1
    t = arb(gate.HEIGHT)
    a = acb(arb("0.75"), -t / 2)
    b = acb(arb("0.75"), t / 2)
    beta = a.gamma() * b.gamma() / acb(arb("1.5")).gamma()
    phase = acb(0, -arb.pi() / 8).exp()
    c_t = (arb.pi() / (32 * t)) ** arb("0.25")
    alpha = 159_579
    z = acb(0, arb.pi() * alpha * alpha / 4)
    value = c_t * alpha * (phase * beta * z.hypgeom_1f1(a, arb("1.5"))).real
    require(value.contains(0), "altered direct Kummer ball unexpectedly excludes zero")
    require(abs(value).upper() > arb("1e100"), "altered direct Kummer guard lost catastrophic scale")
    return value.str(12, more=True)


def main() -> int:
    require(gate.RESULT.is_file(), "missing bridge-target result")
    require(gate.NOTE.is_file(), "missing bridge-target note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == gate.STEM, "kind drift")
    require(artifact["status"] == "exact_QK_minus_exact_Hardy_upper_bridge_isolated_with_rigorous_sufficient_signed_target", "status drift")
    require(artifact.get("passed") is True, "bridge-target gate failed")

    split_path = REPO_ROOT / artifact["dependencies"]["classical_upper_split"]["path"]
    split = json.loads(split_path.read_text(encoding="utf-8"))
    row = [item for item in split["output_rows"] if item["output_label"] == "t+0.00"]
    require(len(row) == 1, "central split row drift")
    row = row[0]
    ctx.dps = 140
    exact_upper = arb(row["high_precision"]["exact_complementary_upper_ball"])
    classical_main = arb(row["high_precision"]["classical_upper_main_ball"])
    rho = exact_upper - classical_main
    target = artifact["bridge_target_certificate"]
    require(rho.overlaps(arb(target["rho_RS_ball"])), "fresh rho misses saved ball")
    lower = arb(gate.PUBLISHED_QKT_LOWER) - rho.lower()
    upper = arb(gate.PUBLISHED_QKT_UPPER) - rho.upper()
    require(lower.overlaps(arb(target["published_Delta_KU_safe_lower_ball"])), "fresh lower target misses saved ball")
    require(upper.overlaps(arb(target["published_Delta_KU_safe_upper_ball"])), "fresh upper target misses saved ball")
    require(upper < 0, "fresh Delta_KU target no longer excludes zero")
    require((upper - lower).overlaps(arb(target["Delta_KU_target_width_ball"])), "fresh target width misses saved ball")
    require(arb(gate.DISPLAY_DELTA_LOWER) > lower, "display lower endpoint is not conservative")
    require(arb(gate.DISPLAY_DELTA_UPPER) < upper, "display upper endpoint is not conservative")
    require(
        target["published_Delta_KU_sufficient_corridor"]
        == f"{gate.DISPLAY_DELTA_LOWER}<Delta_KU<{gate.DISPLAY_DELTA_UPPER}",
        "displayed Delta_KU corridor drift",
    )
    source_hybrid = arb(row["hybrid_minus_exact_upper_ball"])
    require(source_hybrid.lower() > upper, "fresh source-hybrid guard drift")

    identities = artifact["identity_certificate"]
    require(identities["classical_bridge"] == "Q_K-T=Delta_KU+rho_RS", "classical bridge drift")
    require(identities["joined_bridge"] == "J_Z=Delta_KU+rho_RS-(G-T)-A_transition", "joined bridge drift")
    lattice = artifact["argument_lattice_audit"]
    require("2*pi*i" in lattice["odd_step"], "odd argument step drift")
    require("not periodic" in lattice["nonperiodicity_guard"], "Kummer periodicity guard missing")

    decision = artifact["decision"]
    require(decision["single_exact_bridge_Delta_KU_isolated"] is True, "single bridge not isolated")
    require(decision["source_hybrid_telemetry_permitted_as_QK_proxy"] is False, "telemetry guard lost")
    require(decision["naive_direct_fixed_precision_1F1_roster_evaluator_selected"] is False, "naive direct evaluator selected")
    require(decision["exact_truncation_continuation_bridge_selected"] is True, "exact continuation route not selected")
    require(decision["Delta_KU_enclosed"] is False, "Delta_KU overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")

    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file {path}")
            require(file_hash(path) == record["sha256"], f"hash drift {path}")
    paper = artifact["references"]["paper"]
    paper_path = REPO_ROOT / paper["path"]
    require(file_hash(paper_path) == paper["sha256"], "paper hash drift")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("not a proof of the bridge value", "definite negative", "not substituted", "not periodic", "No enclosure"):
        require(token in note, f"note boundary token missing: {token}")
    altered = altered_direct_kummer_guard()
    print(f"independently checked exact Q_K-to-Hardy bridge target; altered direct ball {altered}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
