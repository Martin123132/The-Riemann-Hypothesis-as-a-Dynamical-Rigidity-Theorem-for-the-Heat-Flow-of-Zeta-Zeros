#!/usr/bin/env python3
"""Independently check the 399-mode branch reflection-pairing gate."""

from __future__ import annotations

import cmath
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_opposite_branch_roster_reflection_pairing_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(Path(__file__).name.removeprefix("check_"))
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finite_ratio(numerator_multiple: int, x: float) -> float:
    denominator = cmath.sin(2 * x)
    if abs(denominator) < 1e-12:
        # The samples below use x=0 only among the removable points.
        require(abs(x) < 1e-15, "unexpected removable-zero sample")
        return numerator_multiple / 2
    return (cmath.sin(numerator_multiple * x) / denominator).real


def direct_old(x: float, c_plus: complex) -> complex:
    c_minus = c_plus.conjugate()
    negative = sum(cmath.exp(1j * (4 * j + 1) * x) for j in range(200))
    positive = sum(cmath.exp(-1j * (4 * j + 3) * x) for j in range(199))
    return c_plus * negative + c_minus * positive


def paired_old(x: float, c_plus: complex) -> complex:
    c_minus = c_plus.conjugate()
    ratio = finite_ratio(398, x)
    return cmath.exp(1j * 797 * x) * c_plus + cmath.exp(-1j * x) * ratio * (
        cmath.exp(1j * 398 * x) * c_plus + cmath.exp(-1j * 398 * x) * c_minus
    )


def direct_new(x: float, c_plus: complex) -> complex:
    c_minus = c_plus.conjugate()
    negative = sum(cmath.exp(1j * (4 * j + 3) * x) for j in range(199))
    positive = sum(cmath.exp(-1j * (4 * j + 1) * x) for j in range(200))
    return c_plus * negative + c_minus * positive


def direct_boundary(x: float, c_plus: complex) -> complex:
    c_minus = c_plus.conjugate()
    negative = sum(cmath.exp(1j * (4 * j + 1) * x) for j in range(199))
    positive = sum(cmath.exp(-1j * (4 * j + 3) * x) for j in range(200))
    return c_plus * negative + c_minus * positive


def paired_new(x: float, c_plus: complex) -> complex:
    c_minus = c_plus.conjugate()
    ratio = finite_ratio(398, x)
    return cmath.exp(-1j * 797 * x) * c_minus + cmath.exp(1j * x) * ratio * (
        cmath.exp(1j * 398 * x) * c_plus + cmath.exp(-1j * 398 * x) * c_minus
    )


def paired_boundary(x: float, c_plus: complex) -> complex:
    c_minus = c_plus.conjugate()
    ratio = finite_ratio(398, x)
    return cmath.exp(-1j * 799 * x) * c_minus + cmath.exp(-1j * x) * ratio * (
        cmath.exp(1j * 398 * x) * c_plus + cmath.exp(-1j * 398 * x) * c_minus
    )


def main() -> int:
    require(RESULT.is_file() and NOTE.is_file(), "missing reflection-pairing artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact["passed"] is True, "artifact is not passed")
    require(artifact["decision"]["grouped_opposite_branch_integral_bound_proved"] is False, "integral bound overpromoted")
    require(artifact["decision"]["edge_plus_endpoint_half_current_cancellation_proved"] is False, "endpoint cancellation overpromoted")
    require(artifact["decision"]["edge_exchange_alone_reassembles_selector_jump"] is False, "mode-set exchange overpromoted")

    old_modes = range(39_695, 40_094)
    new_modes = range(39_696, 40_095)
    old_labels = [4 * mode - 159_577 for mode in old_modes]
    new_labels = [4 * mode - 159_579 for mode in new_modes]
    boundary_labels = [4 * mode - 159_577 for mode in new_modes]
    require(old_labels == list(range(-797, 796, 4)), "independent old roster reconstruction failed")
    require(new_labels == list(range(-795, 798, 4)), "independent adjacent roster reconstruction failed")
    require(boundary_labels == list(range(-793, 800, 4)), "independent boundary roster reconstruction failed")
    require(set(old_modes) - set(new_modes) == {39_695}, "independent dropped-edge check failed")
    require(set(new_modes) - set(old_modes) == {40_094}, "independent added-edge check failed")

    samples = [0.0, 0.0017, 0.017, 0.071, 0.233, 0.619, 1.117]
    amplitudes = [0.31 - 0.17j, -0.08 + 0.22j]
    for x in samples:
        for c_plus in amplitudes:
            old_error = abs(direct_old(x, c_plus) - paired_old(x, c_plus))
            new_error = abs(direct_new(x, c_plus) - paired_new(x, c_plus))
            boundary_error = abs(direct_boundary(x, c_plus) - paired_boundary(x, c_plus))
            scale = 1 + 399 * abs(c_plus)
            require(old_error < 2e-11 * scale, f"old numeric pairing failed at x={x}")
            require(new_error < 2e-11 * scale, f"adjacent numeric pairing failed at x={x}")
            require(boundary_error < 2e-11 * scale, f"boundary numeric pairing failed at x={x}")

    sources = artifact["sources"]
    require(sources["builder"]["sha256"] == file_hash(BUILDER), "builder hash drift")
    require(sources["checker"]["sha256"] == file_hash(CHECKER), "checker hash drift")
    require(sources["note"]["sha256"] == file_hash(NOTE), "note hash drift")
    note = NOTE.read_text(encoding="utf-8")
    for marker in ("199 conjugate pairs plus", "labels are therefore `-793,-789,...,799`", "only a finite mode-set identity", "does not bound"):
        require(marker in note, f"note marker missing: {marker}")

    print("independently checked 399-mode branch pairing: 3 rosters, 199 pairs, 7 phase samples", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
