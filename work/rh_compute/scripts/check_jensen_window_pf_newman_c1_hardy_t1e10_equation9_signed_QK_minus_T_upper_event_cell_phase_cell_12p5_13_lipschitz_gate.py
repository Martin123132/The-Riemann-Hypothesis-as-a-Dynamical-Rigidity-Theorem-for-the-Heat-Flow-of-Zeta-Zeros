#!/usr/bin/env python3
"""Independently check the promoted 12.5--13 phase-cell sign gate."""

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

from flint import arb

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_upper_event_cell_phase_cell_12p5_13_lipschitz_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_hash(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def absolute(path_text: str) -> Path:
    return ROOT / Path(path_text)


def main() -> int:
    require(gate.RESULT.is_file(), f"missing result: {gate.RESULT}")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "phase-cell artifact is not passing")
    require(
        artifact.get("status")
        == "rigorous_two_configuration_closed_phase_cell_sign_certificate",
        "phase-cell status drift",
    )
    require(
        artifact.get("checker") == gate.relative(Path(__file__).resolve()),
        "checker path drift",
    )

    for path_text, expected in artifact["sources"].items():
        path = absolute(path_text)
        require(path.is_file(), f"missing source: {path_text}")
        require(file_hash(path) == expected, f"source hash drift: {path_text}")
    for path_text, expected in artifact["dependency_artifacts"].items():
        path = absolute(path_text)
        require(path.is_file(), f"missing dependency artifact: {path_text}")
        require(file_hash(path) == expected, f"dependency hash drift: {path_text}")

    note = absolute(artifact["proof_note"])
    require(note.is_file(), "missing proof note")
    require(file_hash(note) == artifact["proof_note_sha256"], "proof-note hash drift")
    note_text = note.read_text(encoding="utf-8")
    for token in (
        "F_t(y)=F_a(y) exp(-i h log y)",
        "Q(t) = Q(a) + integral_a^t Q'(u) du",
        "Q_K(t)-T(t) < 0",
        "does not assume or prove uniqueness",
    ):
        require(token in note_text, f"proof note lost token: {token}")

    stored = artifact["certificate"]
    rebuilt = gate.assemble_certificate()
    require(
        canonical_hash(stored) == canonical_hash(rebuilt),
        "rebuilt certificate differs from stored certificate",
    )
    require(
        stored["closed_phase_cell"] == [gate.FULL_LEFT, gate.FULL_RIGHT],
        "closed phase-cell endpoints drift",
    )
    require(
        stored["left_closed_half"][-1] == stored["right_closed_half"][0] == gate.ANCHOR,
        "half-cell union misses anchor",
    )
    require(stored["strictly_inside_fixed_roster_event_cell"] is True, "event-cell flag lost")
    require(
        stored["Q_K_minus_T_strictly_negative_on_closed_phase_cell"] is True,
        "closed-cell sign flag lost",
    )

    for variant, panels, precision in (
        ("production", gate.PRODUCTION_PANELS, gate.PRODUCTION_PANEL_PRECISION_BITS),
        ("independent", gate.INDEPENDENT_PANELS, gate.INDEPENDENT_PANEL_PRECISION_BITS),
    ):
        configuration = stored[variant]
        manifest = configuration["panel_manifest"]
        require(manifest["count"] == panels, f"{variant} panel count drift")
        require(manifest["first_index"] == 0, f"{variant} first panel drift")
        require(manifest["last_index"] == panels - 1, f"{variant} last panel drift")
        require(
            configuration["panel_precision_bits"] == precision,
            f"{variant} panel precision drift",
        )
        anchor_value = gate.real(configuration["anchor_point_row"]["Q_K_minus_T"])
        distance = arb(stored["maximum_anchor_distance"])
        threshold = arb(configuration["mean_value_threshold"])
        require(anchor_value.upper() < 0, f"{variant} anchor lost negativity")
        for side in ("left_half", "right_half"):
            half = configuration[side]
            derivative = gate.real(half["model"]["Q_K_minus_T_derivative"])
            modulus = abs(derivative).upper()
            saved_modulus = arb(half["derivative_modulus_upper"])
            require(saved_modulus.contains(modulus), f"{variant} {side} modulus drift")
            require(modulus < threshold, f"{variant} {side} misses threshold")
            transported_upper = (anchor_value.upper() + distance * modulus).upper()
            require(
                arb(half["transported_Q_K_minus_T_upper"]).contains(transported_upper),
                f"{variant} {side} transported upper drift",
            )
            require(transported_upper < 0, f"{variant} {side} lost negativity")
            require(
                half["strictly_negative_on_closed_half"] is True,
                f"{variant} {side} sign flag lost",
            )

    require(stored["production_independent_anchor_overlap"] is True, "anchor overlap flag lost")
    require(
        stored["production_independent_left_derivative_overlap"] is True,
        "left derivative overlap flag lost",
    )
    require(
        stored["production_independent_right_derivative_overlap"] is True,
        "right derivative overlap flag lost",
    )
    print(
        "validated phase-cell 12.5--13 Lipschitz gate: "
        "2 configurations, 2 half-cells, strict Q_K-T negativity",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
