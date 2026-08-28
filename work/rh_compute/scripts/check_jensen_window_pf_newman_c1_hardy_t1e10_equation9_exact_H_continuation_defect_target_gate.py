#!/usr/bin/env python3
"""Independently check exact H and the continuation-defect target."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_exact_H_continuation_defect_target_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_log_H(precision: int = 160) -> arb:
    ctx.dps = precision
    ctx.threads = 1
    t = arb(gate.HEIGHT)
    pi = arb.pi()
    log_thermal = (1 + (-2 * pi * t).exp()).log()
    log_cosh = pi * t - arb(2).log() + log_thermal
    return (
        arb("0.75") * arb(2).log()
        + arb("0.25") * t.log()
        + arb("0.5") * pi.log()
        - arb("0.75") * pi * t
        + arb("0.5") * log_cosh
        - log_thermal
        - acb(arb("0.75"), t / 2).lgamma().real
    )


def main() -> int:
    require(gate.RESULT.is_file(), "missing exact-H result")
    require(gate.NOTE.is_file(), "missing exact-H note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == gate.STEM, "kind drift")
    require(artifact["status"] == "exact_H_normalization_and_lower_adjusted_continuation_defect_target_certified_A73_nonpromotion_guarded", "status drift")
    require(artifact.get("passed") is True, "exact-H gate failed")

    exact_h = artifact["exact_H_certificate"]
    fresh_log_h = independent_log_H()
    require(fresh_log_h.overlaps(arb(exact_h["high_precision_log_H_direct_ball"])), "fresh H misses direct formula")
    require(fresh_log_h.overlaps(arb(exact_h["high_precision_log_H_duplication_ball"])), "fresh H misses duplication formula")
    fresh_h = fresh_log_h.exp()
    require((fresh_h - 1).lower() > 0, "fresh H-1 lost positive sign")
    scaled_h = (fresh_h - 1) * 32 * arb(gate.HEIGHT) ** 2
    require(scaled_h.lower() > 1, "fresh H next-order sign guard failed")
    require(scaled_h.upper() < arb("1.00000000000000000001"), "fresh H asymptotic scale guard failed")

    bridge_path = REPO_ROOT / artifact["dependencies"]["bridge_target"]["path"]
    component_path = REPO_ROOT / artifact["dependencies"]["component_split"]["path"]
    bridge = json.loads(bridge_path.read_text(encoding="utf-8"))["bridge_target_certificate"]
    component = json.loads(component_path.read_text(encoding="utf-8"))
    row = [item for item in component["output_rows"] if item["output_label"] == "t+0.00"]
    require(len(row) == 1, "fresh central component row drift")
    row = row[0]
    lower = arb(row["high_precision"]["exact_lower_component_ball"])
    upper = arb(row["high_precision"]["exact_complementary_upper_ball"])
    delta_lower = arb(bridge["published_Delta_KU_safe_lower_ball"])
    delta_upper = arb(bridge["published_Delta_KU_safe_upper_ball"])
    one_minus_h_u = (1 - fresh_h) * upper
    d_lower = (one_minus_h_u - fresh_h * delta_upper).upper()
    d_upper = (one_minus_h_u - fresh_h * delta_lower).lower()
    target = artifact["continuation_defect_target"]
    require(d_lower.overlaps(arb(target["D_K_safe_lower_ball"])), "fresh D_K lower target misses saved ball")
    require(d_upper.overlaps(arb(target["D_K_safe_upper_ball"])), "fresh D_K upper target misses saved ball")
    require(arb(gate.DISPLAY_DEFECT_LOWER) > d_lower, "displayed D_K lower endpoint is not conservative")
    require(arb(gate.DISPLAY_DEFECT_UPPER) < d_upper, "displayed D_K upper endpoint is not conservative")
    c_lower = d_lower + lower.upper()
    c_upper = d_upper + lower.lower()
    require(c_lower.overlaps(arb(target["C_K_safe_lower_ball"])), "fresh C_K lower target misses saved ball")
    require(c_upper.overlaps(arb(target["C_K_safe_upper_ball"])), "fresh C_K upper target misses saved ball")

    identity = artifact["continuation_defect_identity"]
    require(identity["forward"] == "D_K=(1-H)U-H*Delta_KU", "forward defect identity drift")
    require(identity["inverse"] == "Delta_KU=((1-H)U-D_K)/H", "inverse defect identity drift")
    audit = artifact["paper_exactness_audit"]
    require("approximation" in audit["A64"], "A64 exactness guard missing")
    require("leading order" in audit["A71"], "A71 exactness guard missing")
    require("approximation sign" in audit["A73"], "A73 exactness guard missing")
    require("O(1/t)" in audit["A74"], "A74 exactness guard missing")

    decision = artifact["decision"]
    require(decision["exact_H_derived_from_equation4_prefactor"] is True, "exact H derivation lost")
    require(decision["H_normalization_can_explain_required_bridge_scale"] is False, "H scale was overpromoted")
    require(
        decision["equation4_infinite_Kummer_series_equals_Z_certified_here"] is False,
        "equation-(4) exactness overpromoted",
    )
    require(decision["A73_permitted_as_exact_continuation_bridge"] is False, "A73 nonpromotion guard lost")
    require(decision["exact_equation4_or_Riemann_Siegel_contour_required"] is True, "exact contour route not selected")
    require(decision["D_K_enclosed"] is False, "D_K overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")

    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file {path}")
            require(file_hash(path) == record["sha256"], f"hash drift {path}")
    for record in artifact["references"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"reference hash drift {path}")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("continuation theorem open", "nineteen orders", "positive", "cannot currently certify", "No enclosure"):
        require(token in note, f"note boundary token missing: {token}")
    print(
        "independently checked exact H and continuation-defect target; "
        f"H-1={(fresh_h - 1).str(30, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
