#!/usr/bin/env python3
"""Compose the global order-eleven first-summand curvature theorem."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
LOWER = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_sparse_h23_lower_bridge_certificate.json"
)
COMPACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_compact_adaptive_h23_certificate.json"
)
FINITE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_nested_curvature_finite_ray_certificate.json"
)
ASYMPTOTIC = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_nested_curvature_asymptotic_ray_certificate.json"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_first_summand_curvature_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order11_first_summand_curvature_certificate.md"
)
GENERATOR_PATH = (
    "work/rh_compute/scripts/"
    "jensen_window_pf_compound_order11_first_summand_curvature_certificate.py"
)
CHECKER_PATH = (
    "work/rh_compute/scripts/"
    "check_jensen_window_pf_compound_order11_first_summand_curvature_certificate.py"
)
GLOBAL_CURVATURE_CONSTANT = 6000
LOWER_THEOREM = "y_1''(t)<=6000/t^2 for every real 1252<=t<=5700"
COMPACT_THEOREM = "y_1''(t)<=6000/t^2 for every real 5700<=t<=38020"
FINITE_THEOREM = (
    "y_1''(t)<=6000/t^2 for every saddle mode 2001/1000<=u<=20"
)
ASYMPTOTIC_THEOREM = "t^2*y_1''(t)<2000 for every mode u>=20"
GLOBAL_THEOREM = "y_1''(t)<=6000/t^2 for every real t>=1252"


@dataclass(frozen=True)
class CompositionRow:
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


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT).as_posix()


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def ready_formula(artifact: dict, row_id: str) -> str:
    rows = [row for row in artifact.get("rows", []) if row.get("id") == row_id]
    if len(rows) != 1 or rows[0].get("readiness") != "ready_to_apply":
        raise RuntimeError(f"source row is not ready: {row_id}")
    return rows[0]["formula"]


def source_record(path: Path, artifact: dict, theorem: str) -> dict:
    return {
        "path": relative(path),
        "sha256": sha256(path),
        "kind": artifact["kind"],
        "status": artifact["status"],
        "theorem": theorem,
    }


def validate_sources() -> tuple[dict, dict, dict, dict, list[dict]]:
    lower = load(LOWER)
    compact = load(COMPACT)
    finite = load(FINITE)
    asymptotic = load(ASYMPTOTIC)
    if (
        lower.get("kind")
        != "jensen_window_pf_compound_order11_sparse_h23_lower_bridge_certificate"
        or lower.get("status")
        != "rigorous order-eleven first-summand curvature theorem on 1252<=t<=5700"
        or lower.get("theorem") != LOWER_THEOREM
        or lower.get("summary", {}).get("lower_bridge_theorems") != 1
        or lower.get("summary", {}).get("global_first_summand_theorems") != 0
    ):
        raise RuntimeError("invalid lower-bridge source")
    if (
        compact.get("kind")
        != "jensen_window_pf_compound_order11_compact_adaptive_h23_certificate"
        or compact.get("status")
        != (
            "rigorous order-eleven first-summand curvature theorem on "
            "5700<=t<=38020"
        )
        or compact.get("theorem") != COMPACT_THEOREM
        or compact.get("summary", {}).get("segments") != 2020
        or compact.get("summary", {}).get("quarter_blocks") != 129280
        or compact.get("summary", {}).get("compact_first_summand_theorems") != 1
        or compact.get("summary", {}).get("global_first_summand_theorems") != 0
        or compact.get("summary", {}).get("full_kernel_theorems") != 0
    ):
        raise RuntimeError("invalid compact-bridge source")
    if (
        finite.get("kind")
        != "jensen_window_pf_compound_order11_nested_curvature_finite_ray_certificate"
        or finite.get("status")
        != (
            "rigorous order-eleven first-summand curvature theorem on "
            "2001/1000<=u<=20"
        )
        or finite.get("theorem") != FINITE_THEOREM
        or finite.get("finite_ray", {}).get("mode_range")
        != ["2001/1000", "20"]
        or finite.get("finite_ray", {}).get("all_blocks_passed") is not True
        or finite.get("summary", {}).get("finite_ray_theorems") != 1
    ):
        raise RuntimeError("invalid finite-ray source")
    asymptotic_formula = ready_formula(
        asymptotic,
        "co11ncarc_03_dimensionless_box",
    )
    if (
        asymptotic.get("kind")
        != "jensen_window_pf_compound_order11_nested_curvature_asymptotic_ray_certificate"
        or asymptotic.get("status")
        != "rigorous order-eleven first-summand curvature theorem on u>=20"
        or asymptotic_formula != ASYMPTOTIC_THEOREM
        or asymptotic.get("summary", {}).get("asymptotic_ray_theorems") != 1
    ):
        raise RuntimeError("invalid asymptotic-ray source")
    sources = [
        source_record(path, artifact, theorem)
        for path, artifact, theorem in (
            (LOWER, lower, LOWER_THEOREM),
            (COMPACT, compact, COMPACT_THEOREM),
            (FINITE, finite, FINITE_THEOREM),
            (ASYMPTOTIC, asymptotic, ASYMPTOTIC_THEOREM),
        )
    ]
    return lower, compact, finite, asymptotic, sources


def build_artifact() -> dict:
    lower, compact, finite, asymptotic, sources = validate_sources()
    scaled_uppers = {
        "lower": lower["summary"]["largest_scaled_curvature_upper"],
        "compact": compact["summary"]["largest_scaled_curvature_upper"],
        "finite_ray": finite["finite_ray"]["largest_scaled_curvature_upper"],
        "asymptotic_ray": asymptotic["dimensionless_interval"][
            "scaled_curvature_upper"
        ],
    }
    if any(
        Decimal(value) >= Decimal(GLOBAL_CURVATURE_CONSTANT)
        for value in scaled_uppers.values()
    ):
        raise RuntimeError("a source range does not fit below the 6000 global cap")
    largest_name = max(scaled_uppers, key=lambda name: Decimal(scaled_uppers[name]))
    transition_upper = compact["summary"]["saddle_transition_upper"]
    if Decimal(transition_upper) >= Decimal(38020):
        raise RuntimeError("compact bridge does not overlap the finite saddle ray")
    monotonicity = compact["source_contract"]["saddle_monotonicity"]
    if (
        monotonicity.get("formula")
        != "V''>0 for u>=1/100; V'<0 for 0<u<=1/100"
    ):
        raise RuntimeError("compact source does not bind saddle monotonicity")
    rows = [
        CompositionRow(
            "co11fscc_01_lower_bridge",
            "interval_theorem",
            "ready_to_apply",
            "The sparse-H23 continuation proves the lower real-t interval.",
            LOWER_THEOREM,
            "First Newman summand on the finite lower interval.",
        ),
        CompositionRow(
            "co11fscc_02_compact_bridge",
            "interval_theorem",
            "ready_to_apply",
            "The adaptive-H23 compact certificate extends the real-t theorem beyond the u=2.001 handoff.",
            COMPACT_THEOREM,
            "First Newman summand on the compact interval.",
            {"V_prime_2_001_upper": transition_upper},
        ),
        CompositionRow(
            "co11fscc_03_saddle_rays",
            "interval_analytic_theorem",
            "ready_to_apply",
            "Finite and asymptotic saddle certificates cover every mode u>=2.001.",
            "[2001/1000<=u<=20] union [u>=20]",
            "Uses the hash-bound strictly increasing saddle parametrization t=V'(u).",
        ),
        CompositionRow(
            "co11fscc_04_overlap",
            "exact_overlap",
            "ready_to_apply",
            "The compact real-t interval overlaps the beginning of the complete saddle ray.",
            "V'(2001/1000)<38020",
            "Joins first-summand curvature theorems only.",
        ),
        CompositionRow(
            "co11fscc_05_global_composition",
            "exact_theorem_composition",
            "ready_to_apply",
            "The lower, compact, finite-ray, and asymptotic ranges cover the complete half-line.",
            GLOBAL_THEOREM,
            "First Newman summand at lambda=-100 only; no full-kernel promotion.",
            {
                "scaled_uppers": scaled_uppers,
                "largest_source": largest_name,
                "largest_scaled_curvature_upper": scaled_uppers[largest_name],
            },
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order11_first_summand_curvature_certificate",
        "date": "2026-07-22",
        "status": (
            "rigorous global order-eleven first-summand curvature theorem "
            "on t>=1252"
        ),
        "proof_boundary": (
            "This artifact composes the continuous first-summand theorem at "
            "lambda=-100 only. It does not prove full-Newman-kernel entry, "
            "heat-forward invariance, all-shift order eleven, PF-infinity, "
            "Lambda<=0, or RH."
        ),
        "theorem": GLOBAL_THEOREM,
        "source_contract": {
            "sources": sources,
            "lower_compact_join": "t=5700",
            "compact_saddle_overlap": (
                f"V'(2001/1000)<={transition_upper}<38020"
            ),
            "finite_asymptotic_join": "u=20",
            "saddle_monotonicity": monotonicity,
            "scaled_uppers": scaled_uppers,
        },
        "largest_scaled_curvature_upper": scaled_uppers[largest_name],
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": len(rows),
            "open_rows": 0,
            "global_first_summand_curvature_theorems": 1,
            "full_kernel_theorems": 0,
            "heat_forward_theorems": 0,
            "rh_claims": 0,
        },
        "generator": GENERATOR_PATH,
        "checker": CHECKER_PATH,
    }


def write_note(path: Path, artifact: dict) -> None:
    source = artifact["source_contract"]
    lines = [
        "# Order-Eleven Global First-Summand Curvature Certificate",
        "",
        "Date: 2026-07-22",
        "",
        f"Status: **{artifact['status']}**. This is not a proof of RH or `Lambda <= 0`.",
        "",
        "## Theorem",
        "",
        f"`{artifact['theorem']}`.",
        "",
        "The four exact source ranges are:",
        "",
        f"- `{LOWER_THEOREM}`",
        f"- `{COMPACT_THEOREM}`",
        f"- `{FINITE_THEOREM}`",
        f"- `{ASYMPTOTIC_THEOREM}`",
        "",
        (
            "The compact-saddle overlap is "
            f"`{source['compact_saddle_overlap']}`."
        ),
        (
            "The largest scaled upper across the four sources is "
            f"`{artifact['largest_scaled_curvature_upper']}<6000`."
        ),
        "",
        "## Boundary",
        "",
        artifact["proof_boundary"],
        "",
        "This is not a proof of RH.",
        "",
        "## Reproduce",
        "",
        "```powershell",
        f"python {GENERATOR_PATH}",
        f"python {CHECKER_PATH}",
        "```",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


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
        "wrote global order-eleven first-summand curvature theorem: "
        f"{artifact['theorem']}; largest scaled upper "
        f"{artifact['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
