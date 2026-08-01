#!/usr/bin/env python3
"""Check the strong-log-concave local quartic countermodel."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402

import jensen_window_pf_strong_logconcave_local_quartic_countermodel as model  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_strong_logconcave_local_quartic_countermodel.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_strong_logconcave_local_quartic_countermodel.md"
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    note = args.note.read_text(encoding="utf-8")
    issues: list[str] = []

    if payload.get("kind") != "jensen_window_pf_strong_logconcave_local_quartic_countermodel":
        issues.append("artifact kind changed")
    if payload.get("status") != "exact interval local quartic countermodel gate":
        issues.append("artifact status changed")
    rows = payload.get("rows", [])
    expected_ids = [f"slcq_{index:02d}_" for index in range(1, 11)]
    if len(rows) != 10:
        issues.append(f"expected 10 rows, found {len(rows)}")
    for prefix, row in zip(expected_ids, rows):
        if not row.get("id", "").startswith(prefix):
            issues.append(f"row order/id changed at {prefix}")

    summary = payload.get("summary", {})
    expected_summary = {
        "rows": 10,
        "positive_moments": 6,
        "local_ratio_wall_contractions": 4,
        "positive_monotone_gaps": 3,
        "strict_cubic_margins": 3,
        "negative_quartic_frontiers": 1,
        "full_support_approximation_theorems": 1,
        "open_xi_handoffs": 1,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(f"summary {key} changed")

    flint.ctx.prec = 640
    moments = model.normalized_moments()
    contractions = [
        moments[k + 1] * moments[k - 1] / moments[k] ** 2
        for k in range(1, 5)
    ]
    lower_walls = [
        flint.arb(1) / 3,
        flint.arb(3) / 5,
        flint.arb(5) / 7,
        flint.arb(7) / 9,
    ]
    gaps = [contractions[k + 1] - contractions[k] for k in range(3)]
    cubics = [
        model.cubic_frontier(contractions[k], contractions[k + 1])
        for k in range(3)
    ]
    quartic_q = model.quartic_frontier(*contractions[:3])
    if not all(bool(value > 0) for value in moments):
        issues.append("recomputed moment positivity failed")
    if not all(
        bool(contractions[k] > lower_walls[k] and contractions[k] < 1)
        for k in range(4)
    ):
        issues.append("recomputed local ratio wall failed")
    if not all(bool(value > 0) for value in gaps):
        issues.append("recomputed local monotonicity failed")
    if not all(bool(value < 0) for value in cubics):
        issues.append("recomputed cubic frontier failed")
    if not bool(quartic_q < 0):
        issues.append("recomputed quartic frontier failed")

    stored = payload.get("diagnostics", {})
    for entry, recomputed in zip(stored.get("moments", []), moments):
        if not flint.arb(entry["ball"]).overlaps(recomputed):
            issues.append(f"stored moment ball {entry.get('k')} does not overlap recomputation")
    for entry, recomputed in zip(stored.get("contractions", []), contractions):
        if not flint.arb(entry["ball"]).overlaps(recomputed):
            issues.append(
                f"stored contraction ball {entry.get('k')} does not overlap recomputation"
            )
    for entry, recomputed in zip(stored.get("monotone_gaps", []), gaps):
        if not flint.arb(entry["ball"]).overlaps(recomputed):
            issues.append(
                f"stored monotone gap {entry.get('from_k')} does not overlap recomputation"
            )
    for entry, recomputed in zip(stored.get("cubic_frontiers", []), cubics):
        if not flint.arb(entry["ball"]).overlaps(recomputed):
            issues.append(
                f"stored cubic frontier {entry.get('shift')} does not overlap recomputation"
            )
    stored_q = stored.get("quartic_Q", {}).get("ball")
    if not stored_q or not flint.arb(stored_q).overlaps(quartic_q):
        issues.append("stored quartic frontier does not overlap recomputation")

    w, x, y, z = sp.symbols("w x y z")
    quartic = 1 + 4 * w + 6 * x * w**2 + 4 * x**2 * y * w**3 + x**3 * y**2 * z * w**4
    discriminant = sp.factor(sp.discriminant(quartic, w))
    expected_q = (
        x**3 * y**4 * z**3
        - 12 * x**2 * y**3 * z**2
        - 18 * x**2 * y**2 * z**2
        + 54 * x**2 * y**2 * z
        - 27 * x**2 * y**2
        + 54 * x * y**2 * z**2
        - 6 * x * y**2 * z
        - 180 * x * y * z
        + 108 * x * y
        + 81 * x * z
        - 54 * x
        - 27 * y**2 * z**2
        + 108 * y * z
        - 64 * y
        - 54 * z
        + 36
    )
    if sp.expand(discriminant - 256 * x**6 * y**2 * expected_q) != 0:
        issues.append("independent quartic discriminant factorization failed")

    required_note_markers = [
        "1/5`-strongly log-concave",
        "x_1<x_2<x_3<x_4",
        "three consecutive cubic frontier balls",
        "negative discriminant",
        "Full-Support Guard",
        "dominated convergence",
        "not a Newman heat trajectory",
        "outputs/jensen_window_pf_quartic_double_root_threshold_lemma.md",
    ]
    for marker in required_note_markers:
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    for source in payload.get("sources", []):
        if not (REPO_ROOT / source).exists():
            issues.append(f"source missing: {source}")

    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated Jensen-window PF strong-log-concave local quartic "
        "countermodel: 10 rows, 0 issues, 6 positive Mellin moments, "
        "4 local ratio-wall contractions, 3 monotone gaps, 3 strict cubic "
        "margins, 1 negative quartic frontier, 1 full-support approximation "
        "theorem, 1 Xi-specific handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
