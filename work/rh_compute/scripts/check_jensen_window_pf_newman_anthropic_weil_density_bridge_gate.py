#!/usr/bin/env python3
"""Independently check the Anthropic Weil/rank-trace density audit."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re

import mpmath as mp
from pypdf import PdfReader
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_anthropic_weil_density_bridge_gate"
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE = REPO_ROOT / "outputs" / f"{KIND}.md"
EXPECTED_PDFS = {
    REPO_ROOT / "work/rh_compute/external/anthropic_zeta23_20260810/claude_more_than_two_thirds_zeta_zeros.pdf": "6792988e6cd0e17690621ce898abd5d534f98407741bc7cb14bbe7d07c77d72f",
    REPO_ROOT / "work/rh_compute/external/anthropic_zeta23_20260810/anthropic_informal_note_67_percent.pdf": "45e0330ad37965e5531fa1f4f11e5bebcae147a5237a3e5b3d029efa7ddf759d",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compact_pdf_text(path: Path) -> str:
    text = "\n".join(page.extract_text() or "" for page in PdfReader(path).pages)
    return re.sub(r"[^a-z0-9]+", "", text.lower())


def main() -> None:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == KIND, "kind drift")
    require(artifact["passed"] is True, "artifact did not pass")
    decision = artifact["decision"]
    require(decision["external_statement_source_audited"] is True, "source audit drift")
    require(decision["rank_trace_core_checked"] is True, "rank-trace drift")
    require(decision["exceptional_cap_derived"] is True, "exceptional-cap drift")
    require(decision["full_external_analytic_proof_independently_reproved"] is False, "analytic-proof overpromotion")
    require(decision["community_acceptance_treated_as_settled"] is False, "acceptance overpromotion")
    require(decision["density_amplification_proved"] is False, "density overpromotion")
    require(decision["heat_dependent_Weil_inertia_transport_proved"] is False, "inertia overpromotion")
    require(decision["direct_RH_bridge"] is False, "RH overpromotion")
    require(decision["main_Airy_route_changed"] is False, "Airy route drift")
    require(len(artifact["rows"]) == 12, "row count drift")

    for path, expected in EXPECTED_PDFS.items():
        require(file_hash(path) == expected, f"PDF hash drift: {path}")
    paper = compact_pdf_text(next(iter(EXPECTED_PDFS)))
    require("ranktraceinequality" in paper, "rank-trace source token missing")
    require("rhitselfisoutofreachofthemechanism" in paper, "method-limit source token missing")

    theta = 1 / sp.sqrt(2)
    c1 = sp.sqrt(2) * sp.tan(theta) / (1 + theta * sp.tan(theta))
    first = sp.simplify(2 - 1 / c1)
    second = sp.Rational(3, 2) - sp.cot(theta) / sp.sqrt(2)
    require(sp.simplify(sp.trigsimp(first - second)) == 0, "constant identity failed")

    mp.mp.dps = 120
    th = 1 / mp.sqrt(2)
    cc1 = mp.sqrt(2) * mp.tan(th) / (1 + th * mp.tan(th))
    critical = 2 - 1 / cc1
    exceptional = 1 - critical
    pair = exceptional / 2
    saved = artifact["constant_certificate"]
    storage_tolerance = mp.mpf("1e-78")
    require(abs(mp.mpf(saved["simple_on_line_lower_constant"]) - critical) < storage_tolerance, "saved critical constant drift")
    require(abs(mp.mpf(saved["exceptional_fraction_upper_cap"]) - exceptional) < storage_tolerance, "saved exceptional cap drift")
    require(abs(mp.mpf(saved["off_line_pair_fraction_threshold"]) - pair) < storage_tolerance, "saved pair threshold drift")

    c, r, b = sp.symbols("c r b", positive=True)
    p_trace = c * r / 2
    q_trace = c * b
    lhs = c**2 * r / 4 + c**2 * b
    rhs = c * p_trace - c**2 * r / 4 + 2 * c * q_trace - c**2 * b
    require(sp.expand(lhs - rhs) == 0, "rank-trace equality model failed")

    t = sp.symbols("t", positive=True)
    n_main = t * sp.log(t / (2 * sp.pi)) / (2 * sp.pi)
    for descendants in (1, sp.log(t), sp.sqrt(t), t):
        require(sp.limit(descendants / n_main, t, sp.oo) == 0, f"density scale failed: {descendants}")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path}")

    note = NOTE.read_text(encoding="utf-8")
    note_compact = re.sub(r"\s+", " ", note)
    for token in (
        "Required Density Amplification",
        "An infinite descendant chain is nowhere near sufficient",
        "Calling Jensen-degree growth, shift growth, or local contact index a density",
        "Do not pivot the main proof programme",
        "Pi Provenance",
    ):
        require(token in note_compact, f"note token missing: {token}")

    print(
        "validated Anthropic Weil density-bridge audit independently: "
        f"c*={mp.nstr(critical, 16)}, delta*={mp.nstr(exceptional, 16)}, "
        "12 rows, density bridge open, Airy route unchanged"
    )


if __name__ == "__main__":
    main()
