#!/usr/bin/env python3
"""Independently replay the unequal-truncation saddle-map route audit."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import sys

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_unequal_truncation_saddle_map_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fresh_roots(t: mp.mpf, label: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    delta = mp.sqrt(label * label - 8 * t / mp.pi)
    return (label - delta) / 4, (label + delta) / 4


def fresh_map(t: mp.mpf, y: mp.mpf) -> tuple[mp.mpf, mp.mpf, mp.mpf, mp.mpf]:
    p = t / (2 * mp.pi)
    partner = p / y
    alpha = max(y, partner)
    beta = min(y, partner)
    return alpha, beta, mp.sqrt(alpha / beta), 2 * (y + partner)


def main() -> int:
    priority = gate.set_low_priority()
    require(priority == "below_normal_one_cpu", f"checker resource cap unavailable: {priority}")
    require(gate.RESULT.is_file(), "missing saddle-map result")
    require(gate.NOTE.is_file(), "missing saddle-map note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("kind") == gate.STEM, "kind drift")
    require(artifact.get("status") == "unequal_truncation_saddle_map_and_route_split_certified", "status drift")
    require(artifact.get("passed") is True, "gate did not pass")
    require(digest(gate.DEPENDENCY) == artifact["dependency"]["sha256"], "dependency hash drift")
    require(digest(gate.BUILDER) == artifact["source_hashes"]["builder"], "builder hash drift")
    require(digest(Path(__file__)) == artifact["source_hashes"]["checker"], "checker hash drift")

    audit = artifact["literature_audit"]
    require(audit["direct_identity_with_joined_heat_flow_packet_proved"] is False, "unsafe J_Z identification")
    require(audit["explicit_numerical_error_constants_supplied_by_paper"] is False, "implicit constants promoted")
    require(audit["theorem_1_5_scales_are_certified_bounds_here"] is False, "shape promoted to bound")
    require(audit["theorem_2_1_scales_are_certified_bounds_here"] is False, "shape promoted to bound")

    mp.mp.dps = 100
    t = mp.mpf(artifact["scope"]["height"])
    p = t / (2 * mp.pi)
    sqrt_p = mp.sqrt(p)
    c = mp.power(t, mp.mpf(1) / 6)
    lambda_star = (c + mp.sqrt(c * c - 4)) / 2
    y_star = sqrt_p / lambda_star
    split = mp.mpf(math.ceil(float(y_star - mp.mpf("0.5")))) + mp.mpf("0.5")
    require(split == mp.mpf(artifact["actual_height_route_split"]["first_half_integer_not_below_crossover"]), "split drift")

    def ratio(y: mp.mpf) -> mp.mpf:
        _, _, lam, _ = fresh_map(t, y)
        return mp.power(lam + 1 / lam, 3) / mp.sqrt(t)

    require(ratio(split) < 1, "fresh split is not contracting")
    require(ratio(split - 1) > 1, "fresh prior boundary is not noncontracting")
    require(abs(y_star - mp.mpf(artifact["actual_height_route_split"]["y_crossover"])) < mp.mpf("1e-40"), "crossover drift")

    endpoint = artifact["actual_height_route_split"]["endpoint_cells"]
    central = artifact["actual_height_route_split"]["central_cells"]
    require((endpoint["first"], endpoint["last"], endpoint["count"]) == (622, 860, 239), "endpoint cells drift")
    require((central["first"], central["last"], central["count"]) == (861, 39936, 39076), "central cells drift")
    require(endpoint["count"] + central["count"] == artifact["actual_height_route_split"]["owned_cell_count"], "cell total drift")

    for row in artifact["exact_parameter_map"]["label_rows"]:
        label = mp.mpf(row["odd_label"])
        lower, upper = fresh_roots(t, label)
        alpha, beta, lam, mapped_label = fresh_map(t, lower)
        require(abs(lower * upper - p) < mp.mpf("1e-70") * p, "fresh root product failed")
        require(abs(lower + upper - label / 2) < mp.mpf("1e-70") * label, "fresh root sum failed")
        require(abs(mapped_label - label) < mp.mpf("1e-70") * label, "fresh label map failed")
        require(abs(lam - mp.mpf(row["lambda"])) < mp.mpf("1e-40"), "fresh lambda drift")
        require(alpha >= beta >= 1, "invalid unequal truncation ordering")

    # Different exact rational-root fixtures from the builder.
    for beta, alpha, expected_label in (
        (mp.mpf("1.25"), mp.mpf("4.25"), mp.mpf("11")),
        (mp.mpf("2.5"), mp.mpf("5"), mp.mpf("15")),
        (mp.mpf("3.75"), mp.mpf("6.75"), mp.mpf("21")),
    ):
        altered_t = 2 * mp.pi * alpha * beta
        lower, upper = fresh_roots(altered_t, expected_label)
        a2, b2, _, mapped_label = fresh_map(altered_t, beta)
        require(abs(lower - beta) < mp.mpf("1e-90"), "altered lower root failed")
        require(abs(upper - alpha) < mp.mpf("1e-90"), "altered upper root failed")
        require(abs(a2 - alpha) < mp.mpf("1e-90") and abs(b2 - beta) < mp.mpf("1e-90"), "altered map failed")
        require(abs(mapped_label - expected_label) < mp.mpf("1e-90"), "altered label failed")

    decision = artifact["decision"]
    require(decision["J_Z_enclosed"] is False, "J_Z overclaim")
    require(decision["D_K_enclosed"] is False, "D_K overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    require("No identity between O'Sullivan" in artifact["proof_boundary"], "proof boundary drift")
    print("independently checked unequal-truncation saddle map and actual-height route split")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
