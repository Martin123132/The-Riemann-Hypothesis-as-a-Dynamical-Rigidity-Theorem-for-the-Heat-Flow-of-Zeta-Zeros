#!/usr/bin/env python3
"""Compose the global order-ten first-summand curvature theorem."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
LOWER = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_localized_lower_bridge_certificate.json"
)
COMPACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_nested_curvature_compact_certificate.json"
)
FINITE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_nested_curvature_finite_ray_certificate.json"
)
ASYMPTOTIC = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_nested_curvature_asymptotic_ray_certificate.json"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_first_summand_curvature_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order10_first_summand_curvature_certificate.md"
)
GENERATOR_PATH = (
    "work/rh_compute/scripts/"
    "jensen_window_pf_compound_order10_first_summand_curvature_certificate.py"
)
CHECKER_PATH = (
    "work/rh_compute/scripts/"
    "check_jensen_window_pf_compound_order10_first_summand_curvature_certificate.py"
)
GLOBAL_CURVATURE_CONSTANT = 4200
GLOBAL_THEOREM = "z_1''(t)<=4200/t^2 for every real t>=1251"


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
    matches = [row for row in artifact.get("rows", []) if row.get("id") == row_id]
    if len(matches) != 1 or matches[0].get("readiness") != "ready_to_apply":
        raise RuntimeError(f"source row is not ready: {row_id}")
    return matches[0]["formula"]


def validate_sources() -> tuple[dict, dict, dict, dict, list[dict]]:
    lower = load(LOWER)
    compact = load(COMPACT)
    finite = load(FINITE)
    asymptotic = load(ASYMPTOTIC)
    if (
        lower.get("kind")
        != "jensen_window_pf_compound_order10_localized_lower_bridge_certificate"
        or lower.get("status")
        != "rigorous order-ten first-summand curvature theorem on 1251<=t<=5700"
        or lower.get("theorem")
        != "z_1''(t)<=4200/t^2 for every real 1251<=t<=5700"
    ):
        raise RuntimeError("invalid lower-bridge source")
    if (
        compact.get("kind")
        != "jensen_window_pf_compound_order10_nested_curvature_compact_certificate"
        or compact.get("status")
        != "rigorous order-ten first-summand curvature theorem on 5700<=t<=38020"
        or compact.get("theorem")
        != "z_1''(t)<=4200/t^2 for every real 5700<=t<=38020"
        or compact.get("summary", {}).get("compact_first_summand_theorems") != 1
        or compact.get("summary", {}).get("full_kernel_theorems") != 0
    ):
        raise RuntimeError("invalid compact-bridge source")
    finite_formula = ready_formula(finite, "co10ncfr_02_dimensionless_cover")
    if (
        finite.get("kind")
        != "jensen_window_pf_compound_order10_nested_curvature_finite_ray_certificate"
        or finite.get("status")
        != "rigorous order-ten first-summand curvature theorem on 2001/1000<=u<=20"
        or finite_formula
        != "z_1''(t)<=5500/t^2 for every saddle mode 2001/1000<=u<=20"
        or finite.get("finite_ray", {}).get("all_blocks_passed") is not True
    ):
        raise RuntimeError("invalid finite-ray source")
    asymptotic_formula = ready_formula(asymptotic, "co10ncarc_03_dimensionless_box")
    composition_formula = ready_formula(asymptotic, "co10ncarc_04_range_composition")
    if (
        asymptotic.get("kind")
        != "jensen_window_pf_compound_order10_nested_curvature_asymptotic_ray_certificate"
        or asymptotic.get("status")
        != "rigorous order-ten first-summand curvature theorem on u>=20"
        or asymptotic_formula != "t^2*z_1''(t)<1000 for every mode u>=20"
        or composition_formula != "[2001/1000<=u<=20] union [u>=20]"
    ):
        raise RuntimeError("invalid asymptotic-ray source")
    sources = [
        {
            "kind": artifact["kind"],
            "path": relative(path),
            "sha256": sha256(path),
            "status": artifact["status"],
        }
        for path, artifact in (
            (LOWER, lower),
            (COMPACT, compact),
            (FINITE, finite),
            (ASYMPTOTIC, asymptotic),
        )
    ]
    return lower, compact, finite, asymptotic, sources


def lower_scaled_upper(lower: dict) -> str:
    rows = [
        row
        for row in lower["rows"]
        if row.get("id") == "co10llbc_05_composition"
    ]
    return rows[0]["diagnostics"]["global_maximum"]["scaled_curvature_upper"]


def build_artifact() -> dict:
    lower, compact, finite, asymptotic, sources = validate_sources()
    scaled_uppers = {
        "lower": lower_scaled_upper(lower),
        "compact": compact["summary"]["largest_scaled_curvature_upper"],
        "finite_ray": finite["finite_ray"]["largest_scaled_curvature_upper"],
        "asymptotic_ray": asymptotic["dimensionless_interval"][
            "scaled_curvature_upper"
        ],
    }
    largest_name = max(scaled_uppers, key=lambda name: Decimal(scaled_uppers[name]))
    if any(
        Decimal(value) >= Decimal(GLOBAL_CURVATURE_CONSTANT)
        for value in scaled_uppers.values()
    ):
        raise RuntimeError("a source range does not fit below the 4200 global cap")
    transition_upper = compact["summary"]["saddle_transition_upper"]
    if Decimal(transition_upper) >= Decimal(38020):
        raise RuntimeError("compact bridge does not overlap the finite saddle ray")
    rows = [
        CompositionRow(
            "co10fscc_01_lower_bridge",
            "interval_theorem",
            "ready_to_apply",
            "The localized continuation proves the lower real-t interval.",
            "z_1''(t)<=4200/t^2 for every real 1251<=t<=5700",
            "First Newman summand on the finite lower interval.",
        ),
        CompositionRow(
            "co10fscc_02_compact_bridge",
            "interval_theorem",
            "ready_to_apply",
            "The compact saddle certificate extends the real-t theorem beyond the u=2.001 handoff.",
            "z_1''(t)<=4200/t^2 for every real 5700<=t<=38020",
            "First Newman summand on the compact interval.",
            {"V_prime_2_001_upper": transition_upper},
        ),
        CompositionRow(
            "co10fscc_03_saddle_rays",
            "interval_analytic_theorem",
            "ready_to_apply",
            "Finite and asymptotic saddle certificates cover every mode u>=2.001.",
            "[2001/1000<=u<=20] union [u>=20]",
            "Uses the strictly increasing saddle parametrization t=V'(u).",
        ),
        CompositionRow(
            "co10fscc_04_overlap",
            "exact_overlap",
            "ready_to_apply",
            "The real-t compact interval overlaps the beginning of the complete saddle ray.",
            "V'(2001/1000)<38020",
            "Joins two first-summand curvature theorems only.",
        ),
        CompositionRow(
            "co10fscc_05_global_composition",
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
        "kind": "jensen_window_pf_compound_order10_first_summand_curvature_certificate",
        "date": "2026-07-16",
        "status": "rigorous global order-ten first-summand curvature theorem on t>=1251",
        "proof_boundary": (
            "This artifact composes the continuous first-summand theorem at "
            "lambda=-100 only. It does not prove a full-Newman-kernel ceiling, "
            "endpoint entry, heat-forward invariance, the Jensen hierarchy, "
            "Lambda<=0, or RH."
        ),
        "theorem": GLOBAL_THEOREM,
        "source_contract": {
            "sources": sources,
            "lower_compact_join": "5700",
            "compact_saddle_overlap": f"V'(2001/1000)<={transition_upper}<38020",
            "finite_asymptotic_join": "u=20",
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
            "endpoint_entry_theorems": 0,
            "rh_claims": 0,
        },
        "generator": GENERATOR_PATH,
        "checker": CHECKER_PATH,
    }


def write_note(path: Path, artifact: dict) -> None:
    lines = [
        "# Order-ten global first-summand curvature certificate",
        "",
        "Date: 2026-07-17",
        "",
        f"Status: **{artifact['status']}**. This certificate is not a proof of PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "## Theorem",
        "",
        f"`{artifact['theorem']}`.",
        "",
        f"The largest source upper bound is `{artifact['largest_scaled_curvature_upper']}`.",
        "",
        "## Boundary",
        "",
        artifact["proof_boundary"],
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
    artifact = build_artifact()
    DEFAULT_OUT.parent.mkdir(parents=True, exist_ok=True)
    DEFAULT_OUT.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(DEFAULT_NOTE, artifact)
    print(
        "wrote global order-ten first-summand curvature theorem: "
        f"{artifact['theorem']}; largest scaled upper "
        f"{artifact['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
