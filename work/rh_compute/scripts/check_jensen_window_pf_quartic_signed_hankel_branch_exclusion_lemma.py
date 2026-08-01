#!/usr/bin/env python3
"""Check the quartic signed-Hankel branch-exclusion lemma."""

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


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_quartic_signed_hankel_branch_exclusion_lemma.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_quartic_signed_hankel_branch_exclusion_lemma.md"
)
COUNTERMODEL = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_strong_logconcave_local_quartic_countermodel.json"
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    note = args.note.read_text(encoding="utf-8")
    issues: list[str] = []

    if payload.get("kind") != "jensen_window_pf_quartic_signed_hankel_branch_exclusion_lemma":
        issues.append("artifact kind changed")
    if payload.get("status") != "exact Xi quartic branch-exclusion lemma":
        issues.append("artifact status changed")
    rows = payload.get("rows", [])
    if len(rows) != 10:
        issues.append(f"expected 10 rows, found {len(rows)}")
    for index, row in enumerate(rows, start=1):
        if not row.get("id", "").startswith(f"qshb_{index:02d}_"):
            issues.append(f"row order/id changed at {index}")
    expected_summary = {
        "rows": 10,
        "exact_identities": 4,
        "imported_hankel_theorems": 1,
        "excluded_boundary_strata": 2,
        "reduced_outer_thresholds": 1,
        "countermodel_separations": 1,
        "open_outer_thresholds": 1,
    }
    for key, expected in expected_summary.items():
        if payload.get("summary", {}).get(key) != expected:
            issues.append(f"summary {key} changed")

    x, y, z, a, b, c, p = sp.symbols("x y z a b c p")
    normalized = [1, 1, x, x**2 * y, x**3 * y**2 * z]
    hankel = sp.factor(
        sp.det(sp.Matrix([[normalized[i + j] for j in range(3)] for i in range(3)]))
    )
    A = -3 * a**2 + 8 * a + p
    B = -a**2 + 2 * a + p
    boundary = {
        x: A / 6,
        y: 18 * a * B / A**2,
        z: 2 * p * A / (3 * B**2),
    }
    curvature = 3 * a**2 - 4 * a + p
    if sp.expand(sp.factor(hankel.subs(boundary)) + curvature**3 / 216) != 0:
        issues.append("independent boundary Hankel factorization failed")
    if sp.expand(
        (a - b) * (a - c)
        .subs(c, 4 - 2 * a - b)
        - curvature.subs(p, b * (4 - 2 * a - b))
    ) != 0:
        issues.append("independent root-curvature identity failed")

    counter = json.loads(COUNTERMODEL.read_text(encoding="utf-8"))
    moments = [
        flint.arb(entry["ball"]) for entry in counter["diagnostics"]["moments"][:5]
    ]
    a0, a1, a2, a3, a4 = moments
    counter_hankel = (
        a0 * (a2 * a4 - a3**2)
        - a1 * (a1 * a4 - a2 * a3)
        + a2 * (a1 * a3 - a2**2)
    )
    if not bool(counter_hankel > 0):
        issues.append("countermodel separation sign failed")
    stored_counter = rows[8].get("diagnostics", {}).get("countermodel_hankel_ball")
    if not stored_counter or not flint.arb(stored_counter).overlaps(counter_hankel):
        issues.append("stored countermodel Hankel ball changed")

    source_markers = {
        "outputs/jensen_window_pf_reciprocal_defect_compound_order3_gate.md": [
            "D_(3,n)<0 iff C_n",
            "C_n(lambda)>0 for every n>=0 and finite lambda>=-100",
        ],
        "outputs/jensen_window_pf_compound_order3_forward_invariance_certificate.md": [
            "C_n(lambda)>0 for every n>=0 and finite lambda>=-100",
            "D_(3,n)(0)<0 for every n>=0",
        ],
        "outputs/jensen_window_pf_quartic_double_root_threshold_lemma.md": [
            "(3*a^2-4*a+p)*(u-U(a,p))<=0",
        ],
    }
    for source, markers in source_markers.items():
        path = REPO_ROOT / source
        if not path.exists():
            issues.append(f"source missing: {source}")
            continue
        text = path.read_text(encoding="utf-8")
        for marker in markers:
            if marker not in text:
                issues.append(f"source marker missing in {source}: {marker}")

    required_note_markers = [
        "det[B_(i+j)]=-C^3/216",
        "C=3*a^2-4*a+p=(a-b)*(a-c)",
        "middle-root branch `C<0`",
        "triple-root",
        "u<=U(a,p)",
        "Countermodel Separation",
        "D_3=",
        "does not prove `u<=U`",
    ]
    for marker in required_note_markers:
        if marker not in note:
            issues.append(f"note marker missing: {marker}")

    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated Jensen-window PF quartic signed-Hankel branch-exclusion "
        "lemma: 10 rows, 0 issues, 4 exact identities, 1 imported Xi Hankel "
        "theorem, 2 excluded boundary strata, 1 reduced outer threshold, "
        "1 countermodel separation, 1 open outer threshold"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
