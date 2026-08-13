#!/usr/bin/env python3
"""Independently validate endpoint-current moving-roster reassembly."""

from __future__ import annotations

import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

from flint import arb, acb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_endpoint_roster_reassembly_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    priority = gate.set_low_priority()
    require(gate.RESULT.is_file(), "missing endpoint-reassembly result")
    ctx.dps = gate.PRECISION
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "endpoint-reassembly artifact does not pass")
    require(artifact["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash mismatch")
    require(artifact["sources"]["checker"]["sha256"] == gate.file_hash(gate.CHECKER), "checker hash mismatch")
    for name, path in (
        ("primitive_gate", gate.PRIMITIVE_GATE),
        ("roster_gate", gate.ROSTER_GATE),
        ("poisson_gate", gate.POISSON_GATE),
    ):
        require(artifact["dependencies"][name]["sha256"] == gate.file_hash(path), f"{name} hash mismatch")

    c = artifact["certified_reassembly"]
    require(arb(c["h_over_sigma_ball"]).overlaps(arb.pi()), "saved h/sigma mismatch")
    require(c["B_minus_C"] == gate.B - gate.C, "endpoint gap mismatch")
    require(c["B_minus_C_over_4"] == (gate.B - gate.C) // 4, "quarter gap mismatch")
    require((gate.B - gate.C) % 4 == 0 and ((gate.B - gate.C) // 4) % 2 == 1, "independent parity check failed")
    require(c["full_D84_y_B_exact"] == -84, "upper Dirichlet character mismatch")
    require(arb(c["full_D84_y64_magnitude_ball"]) < arb("0.612"), "y=64 kernel decision failed")
    saved_y64 = acb(arb(c["full_D84_y64"]["real_ball"]), arb(c["full_D84_y64"]["imag_ball"]))
    pi = arb.pi()
    endpoint_c = arb(gate.C)
    beta = (pi * endpoint_c**2 / 8) ** (arb(1) / 3)
    h = 4 * beta / endpoint_c
    i = acb(0, 1)
    direct = sum(
        (
            (-i * h * (arb(mode) - endpoint_c / 4) * gate.Y0).exp()
            for mode in range(gate.MODE_LO, gate.MODE_HI + 1)
        ),
        acb(0),
    )
    require(saved_y64.overlaps(direct), "saved y=64 kernel does not overlap independent sum")
    require(c["fixed_roster_cells_checked"] == 169, "cell count mismatch")
    require(arb(c["maximum_y64_cell_reassembly_error_ball"]).contains(0), "cell reassembly error excludes zero")

    decisions = artifact["decision"]
    require(decisions["artificial_y64_current_cancels_by_telescope"] is True, "y=64 telescope missing")
    require(decisions["upper_source_current_bounded"] is False, "upper-current boundary corrupted")
    require(decisions["grouped_bulk_cell_integrals_proved"] is False, "bulk-integral boundary corrupted")
    require(decisions["rh_implication"] is False, "RH boundary corrupted")
    print(f"validated endpoint reassembly independently: 169 cells, D84(yB)=-84, |D84(64)|<0.612; priority={priority}")


if __name__ == "__main__":
    main()
