#!/usr/bin/env python3
"""Compose the global order-nine first-summand curvature theorem."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_first_summand_curvature_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order9_first_summand_curvature_certificate.md"
)
LOWER_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_localized_lower_bridge_certificate.json"
)
COMPACT_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_nested_curvature_compact_certificate.json"
)
FINITE_RAY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_nested_curvature_finite_ray_certificate.json"
)
ASYMPTOTIC_RAY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_nested_curvature_asymptotic_ray_certificate.json"
)

LOWER_THEOREM = "w_1''(t)<=4200/t^2 for every real 1250<=t<=5700"
COMPACT_THEOREM = "w_1''(t)<=4200/t^2 for every real 5700<=t<=V'(2)"
FINITE_RAY_THEOREM = (
    "w_1''(t)<=4200/t^2 for every saddle mode 2<=u<=20"
)
ASYMPTOTIC_THEOREM = "t^2*w_1''(t)<500 for every mode u>=20"
UPPER_COMPOSITION = (
    "[5700<=t<=V'(2)] union [2<=u<=20] union [u>=20]"
)
GLOBAL_THEOREM = "w_1''(t)<=4200/t^2 for every real t>=1250"


@dataclass(frozen=True)
class CurvatureRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def source_record(path: Path, artifact: dict, theorem: str) -> dict:
    return {
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "sha256": sha256(path),
        "kind": artifact["kind"],
        "status": artifact["status"],
        "theorem": theorem,
    }


def validate_sources() -> dict:
    lower = load_json(LOWER_SOURCE)
    compact = load_json(COMPACT_SOURCE)
    finite = load_json(FINITE_RAY_SOURCE)
    asymptotic = load_json(ASYMPTOTIC_RAY_SOURCE)
    if lower.get("theorem") != LOWER_THEOREM:
        raise RuntimeError("localized lower-bridge theorem changed")
    if lower.get("summary", {}).get("curvature_theorems") != 1:
        raise RuntimeError("localized lower bridge is not closed")
    if compact.get("compact", {}).get("theorem") != COMPACT_THEOREM:
        raise RuntimeError("compact curvature theorem changed")
    if compact.get("summary", {}).get("compact_curvature_theorems") != 1:
        raise RuntimeError("compact curvature range is not closed")
    if finite.get("finite_ray", {}).get("theorem") != FINITE_RAY_THEOREM:
        raise RuntimeError("finite-ray curvature theorem changed")
    if finite.get("summary", {}).get("finite_ray_theorems") != 1:
        raise RuntimeError("finite ray is not closed")
    asymptotic_row = next(
        (
            row
            for row in asymptotic.get("rows", [])
            if row.get("id") == "co9ncarc_03_dimensionless_interval"
        ),
        None,
    )
    handoff_row = next(
        (
            row
            for row in asymptotic.get("rows", [])
            if row.get("id") == "co9ncarc_04_global_handoff"
        ),
        None,
    )
    if asymptotic_row is None or asymptotic_row.get("formula") != ASYMPTOTIC_THEOREM:
        raise RuntimeError("asymptotic curvature theorem changed")
    if handoff_row is None or handoff_row.get("formula") != UPPER_COMPOSITION:
        raise RuntimeError("upper-range composition changed")
    if handoff_row.get("readiness") != "ready_to_apply":
        raise RuntimeError("upper-range composition is not ready")
    scaled = Decimal(
        asymptotic.get("dimensionless_interval", {}).get(
            "scaled_curvature_upper", "Infinity"
        )
    )
    if scaled >= Decimal(500):
        raise RuntimeError("asymptotic scaled curvature is not below 500")
    return {
        "sources": [
            source_record(LOWER_SOURCE, lower, LOWER_THEOREM),
            source_record(COMPACT_SOURCE, compact, COMPACT_THEOREM),
            source_record(FINITE_RAY_SOURCE, finite, FINITE_RAY_THEOREM),
            source_record(
                ASYMPTOTIC_RAY_SOURCE,
                asymptotic,
                ASYMPTOTIC_THEOREM,
            ),
        ],
        "lower_largest_scaled_upper": lower["coverage"][
            "largest_scaled_curvature_upper"
        ],
        "compact_largest_scaled_upper": compact["compact"][
            "largest_scaled_curvature_upper"
        ],
        "finite_ray_largest_scaled_upper": finite["finite_ray"][
            "largest_scaled_curvature_upper"
        ],
        "asymptotic_scaled_upper": asymptotic["dimensionless_interval"][
            "scaled_curvature_upper"
        ],
        "upper_composition": handoff_row["formula"],
    }


def build_artifact() -> dict:
    sources = validate_sources()
    largest = max(
        (
            sources["lower_largest_scaled_upper"],
            sources["compact_largest_scaled_upper"],
            sources["finite_ray_largest_scaled_upper"],
            sources["asymptotic_scaled_upper"],
        ),
        key=Decimal,
    )
    rows = [
        CurvatureRow(
            "co9fscc_01_lower_bridge",
            "interval_theorem",
            "ready_to_apply",
            "The localized continuation closes the finite lower real-t bridge.",
            LOWER_THEOREM,
            "Exact real-t range 1250 through 5700.",
        ),
        CurvatureRow(
            "co9fscc_02_compact",
            "interval_theorem",
            "ready_to_apply",
            "The compact saddle certificate joins t=5700 to mode two.",
            COMPACT_THEOREM,
            "Compact saddle range only.",
        ),
        CurvatureRow(
            "co9fscc_03_rays",
            "interval_analytic_theorem",
            "ready_to_apply",
            "Finite and asymptotic saddle certificates cover every mode u>=2.",
            FINITE_RAY_THEOREM + "; " + ASYMPTOTIC_THEOREM,
            "Uses the proved monotone saddle parametrization t=V'(u).",
        ),
        CurvatureRow(
            "co9fscc_04_global_composition",
            "exact_theorem_composition",
            "ready_to_apply",
            "The lower bridge and upper saddle ranges cover the complete half-line.",
            GLOBAL_THEOREM,
            "First Newman summand at lambda=-100 only.",
            {"largest_scaled_curvature_upper": largest},
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order9_first_summand_curvature_certificate",
        "date": "2026-07-14",
        "status": "rigorous global order-nine first-summand curvature theorem on t>=1250",
        "proof_boundary": (
            "This artifact composes the continuous first-summand curvature "
            "theorem only. It does not by itself prove full-kernel endpoint "
            "entry, heat-forward invariance, PF-infinity, RH, or Lambda<=0."
        ),
        "source_contract": sources,
        "theorem": GLOBAL_THEOREM,
        "largest_scaled_curvature_upper": largest,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": len(rows),
            "open_rows": 0,
            "lower_bridge_theorems": 1,
            "upper_range_compositions": 1,
            "global_first_summand_curvature_theorems": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order9_first_summand_curvature_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order9_first_summand_curvature_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    source = artifact["source_contract"]
    lines = [
        "# Order-nine global first-summand curvature certificate",
        "",
        "Date: 2026-07-14",
        "",
        "Status: rigorous global first-summand theorem on `t>=1250`; This is not a proof of full order-nine heat invariance, PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "## Composition",
        "",
        f"- `{LOWER_THEOREM}`",
        f"- `{COMPACT_THEOREM}`",
        f"- `{FINITE_RAY_THEOREM}`",
        f"- `{ASYMPTOTIC_THEOREM}`",
        "",
        "The upper artifact checks the exact handoff",
        "",
        f"`{source['upper_composition']}`.",
        "",
        "Together with the localized bridge this proves",
        "",
        f"`{artifact['theorem']}`.",
        "",
        "Largest scaled curvature upper across the four source ranges:",
        f"`{artifact['largest_scaled_curvature_upper']}<4200`.",
        "",
        "## Boundary",
        "",
        "The next composition applies the existing `<550/k^2` full-kernel",
        "transfer and the completed finite splice. That endpoint theorem is",
        "kept separate so this artifact remains a pure first-summand statement.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    print(
        "wrote global order-nine first-summand curvature certificate: "
        f"largest scaled upper {artifact['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
