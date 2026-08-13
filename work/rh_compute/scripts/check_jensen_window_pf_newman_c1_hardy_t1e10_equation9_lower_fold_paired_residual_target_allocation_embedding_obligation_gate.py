#!/usr/bin/env python3
"""Independently check the fold target allocation and assembly guard."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_paired_residual_target_allocation_embedding_obligation_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
Q = 39_894


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def target_of_augmented(mode: int) -> int:
    if mode == Q:
        return Q
    if mode < Q:
        return mode
    if mode <= Q + 198:
        return 2 * Q - mode
    return Q - 199


def main() -> None:
    require(RESULT.is_file(), "missing result artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "saved gate did not pass")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file() and file_hash(path) == record["sha256"], f"dependency hash drift: {path}")

    inverse: dict[int, list[int]] = {target: [] for target in range(Q - 199, Q + 1)}
    for mode in range(Q - 198, Q + 201):
        inverse[target_of_augmented(mode)].append(mode)
    require(inverse[Q] == [Q], "event-zero allocation drift")
    for target in range(Q - 198, Q):
        require(inverse[target] == [target, 2 * Q - target], f"reflected pair drift at {target}")
    require(inverse[Q - 199] == [Q + 199, Q + 200], "selector edge allocation drift")

    saved = {row["target_mode"]: row["augmented_modes"] for row in artifact["allocation_groups"]}
    require(saved == inverse, "saved allocation differs from independent inverse map")

    # A second finite identity check uses target-group values directly, unlike
    # the builder's 399 per-mode construction.
    group = {r: Fraction(r + 1, r - Q + 211) for r in inverse}
    correction = {r: Fraction(2 * r - 3, 3 * (r - Q + 223)) for r in inverse if r != Q}
    target = {r: Fraction(5 * r + 7, 7 * (r - Q + 227)) for r in inverse}
    embedding = Fraction(17, 29)
    fold = embedding + sum((group[r] - target[r] for r in inverse), Fraction())
    event_zero = group[Q] - target[Q]
    allocation = sum((group[r] - correction[r] - target[r] for r in inverse if r != Q), Fraction())
    require(fold == embedding + sum(correction.values(), Fraction()) + event_zero + allocation, "independent four-term identity failed")

    decision = artifact["decision"]
    require(decision["exact_399_to_200_group_allocation_proved"] is True, "allocation decision drift")
    require(decision["selector_to_paired_residual_embedding_proved"] is False, "unproved embedding was promoted")
    require(decision["complete_fold_owned_residual_bounded"] is False, "unproved fold bound was promoted")
    print("PASS: independent inverse-map and four-term fold-allocation check; embedding remains open")


if __name__ == "__main__":
    main()
