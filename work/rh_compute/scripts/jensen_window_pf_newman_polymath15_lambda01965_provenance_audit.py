#!/usr/bin/env python3
"""Build the published/candidate Newman upper-bound provenance audit."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_polymath15_lambda01965_provenance_audit.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_lambda01965_provenance_audit.md"
)

POLYMATH_URL = "https://arxiv.org/abs/1904.12438"
PLATT_TRUDGIAN_URL = "https://arxiv.org/abs/2004.09765"
RODGERS_TAO_URL = "https://arxiv.org/abs/1801.05914"
ZENODO_URL = "https://doi.org/10.5281/zenodo.20724170"
ZENODO_API_URL = "https://zenodo.org/api/records/20724170"

PT_HEIGHT = 3_000_175_332_800
CANDIDATE_X = 6_000_000_185_827
CANDIDATE_T0 = Fraction(177, 1000)
CANDIDATE_Y0_SQUARED = Fraction(39, 1000)
PUBLISHED_UPPER = Fraction(1, 5)
CANDIDATE_UPPER = Fraction(393, 2000)

ARCHIVE_MD5 = "da7ab6a1bd1e43bb26d080704927cba1"
PAPER_MD5 = "345dbbbb94fc1f96fe01edec0ad58669"
MANIFEST_ENTRIES = 505
MANIFEST_BYTES = 213_197_118


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def fraction_record(value: Fraction) -> dict:
    return {
        "exact": f"{value.numerator}/{value.denominator}",
        "decimal": f"{float(value):.15f}",
    }


PORTABLE_RUNS: tuple[tuple[str, float], ...] = (
    ("stretch_binding", 2.80),
    ("stretch_binding_secondline", 1.87),
    ("seam_exponent", 2.14),
    ("seam_exponent_secondline", 1.83),
    ("seam_kappa", 1.69),
    ("seam_kappa_secondline", 1.66),
    ("seam_ytransfer", 1.94),
    ("seam_ytransfer_secondline", 2.14),
    ("assembly_1965", 1.76),
    ("assembly_secondline", 1.57),
    ("site_glue", 0.35),
    ("site_glue_secondline", 0.25),
    ("tail_bound_box", 2.65),
    ("tail_box_anchor", 8.53),
    ("tail_bound_e3box", 62.03),
    ("tail_e3box_secondline", 71.43),
    ("tbox_blockuniform_crossgate", 30.86),
    ("tbox_blockuniform_anchor197", 11.93),
    ("tbox_blockuniform_anchor19", 32.64),
    ("y_reduction", 7.64),
    ("criterion_theorem", 1.38),
    ("error_terms_audit", 22.09),
)


def build_exact() -> dict:
    candidate_from_row = CANDIDATE_T0 + CANDIDATE_Y0_SQUARED / 2
    if candidate_from_row != CANDIDATE_UPPER:
        raise RuntimeError("candidate row arithmetic drifted")
    x_half = Fraction(CANDIDATE_X, 2)
    height_margin = Fraction(PT_HEIGHT) - x_half
    if height_margin != Fraction(350_479_773, 2):
        raise RuntimeError("Platt-Trudgian height margin drifted")
    improvement = PUBLISHED_UPPER - CANDIDATE_UPPER
    if improvement != Fraction(7, 2000):
        raise RuntimeError("candidate improvement drifted")

    return {
        "published_state": {
            "rodgers_tao_lower": fraction_record(Fraction(0)),
            "polymath15_upper": fraction_record(Fraction(11, 50)),
            "platt_trudgian_upper": fraction_record(PUBLISHED_UPPER),
            "platt_trudgian_height": PT_HEIGHT,
            "status": "peer-reviewed published input",
        },
        "candidate_record": {
            "title": (
                "A certified unconditional upper bound Lambda <= 0.1965 "
                "for the de Bruijn-Newman constant"
            ),
            "creator": "Mosaic Intelligence",
            "publication_date": "2026-06-17",
            "resource_type": "Preprint",
            "doi": "10.5281/zenodo.20724170",
            "paper_banner": "DRAFT v0.3 - not for distribution",
            "readme_status": "PRE-FREEZE DRAFT",
            "ai_disclosure": "fully AI-derived body of results",
            "x": CANDIDATE_X,
            "x_over_2": fraction_record(x_half),
            "pt_height_margin": fraction_record(height_margin),
            "t0": fraction_record(CANDIDATE_T0),
            "y0_squared": fraction_record(CANDIDATE_Y0_SQUARED),
            "claimed_upper": fraction_record(candidate_from_row),
            "improvement_over_published_0p2": fraction_record(improvement),
            "relative_improvement": fraction_record(improvement / PUBLISHED_UPPER),
            "status": "unrefereed candidate record",
        },
        "archive_integrity_audit": {
            "archive_bytes": 195_546_832,
            "archive_md5_expected": ARCHIVE_MD5,
            "archive_md5_observed": ARCHIVE_MD5,
            "paper_bytes": 228_760,
            "paper_md5_expected": PAPER_MD5,
            "zip_entries": 572,
            "manifest_entries": MANIFEST_ENTRIES,
            "manifest_bytes_hashed": MANIFEST_BYTES,
            "manifest_bad": 0,
            "manifest_missing": 0,
            "manifest_elapsed_seconds": 2.73,
        },
        "portable_rerun_audit": {
            "date": "2026-07-17",
            "python": "3.13",
            "mpmath": "1.3.0",
            "sympy": "1.14.0",
            "invocations": [
                {"name": name, "exit_code": 0, "seconds": seconds}
                for name, seconds in PORTABLE_RUNS
            ],
            "invocation_count": len(PORTABLE_RUNS),
            "pass_count": len(PORTABLE_RUNS),
            "fail_count": 0,
            "measured_seconds": 271.18,
            "scope": (
                "Direct portable entry points for every certified1965 package, "
                "the three tbox modes, and the shared criterion/error spines."
            ),
        },
        "unreproduced_boundary": {
            "compiled_packages_not_live_rerun": [
                "record/record_package_197",
                "certified1965/grid_full",
                "certified1965/grid_tbox",
            ],
            "reason": "This Windows audit had no compatible FLINT/Arb build.",
            "additional_open_checks": [
                "independent line-by-line theorem-to-code correspondence review",
                "full wrapper-level 31-verifier clean-copy run",
                "peer review or independent expert review of the preprint",
            ],
        },
        "interval_effect": {
            "published_interval": "0 <= Lambda <= 1/5",
            "candidate_interval_if_accepted": "0 <= Lambda <= 393/2000",
            "qualitative_effect": (
                "The candidate narrows the positive interval by 7/2000 but "
                "does not change the sign question Lambda<=0."
            ),
        },
    }


def build_artifact() -> dict:
    exact = build_exact()
    rows = [
        GateRow(
            id="np15lpa_01_published_lower",
            role="published_input",
            readiness="available_exact",
            claim="Rodgers and Tao prove the de Bruijn-Newman constant is nonnegative.",
            formula="Lambda>=0",
            proof_boundary="Published lower bound only; it does not prove RH.",
        ),
        GateRow(
            id="np15lpa_02_published_upper",
            role="published_input",
            readiness="available_exact",
            claim="Platt and Trudgian's exact verified height yields their published Corollary 2.",
            formula="Lambda<=1/5",
            proof_boundary="The established positive upper bound is not Lambda<=0.",
            diagnostics=exact["published_state"],
        ),
        GateRow(
            id="np15lpa_03_candidate_source",
            role="source_provenance",
            readiness="diagnostic_validated",
            claim="A June 2026 Zenodo preprint claims the sharper bound Lambda<=393/2000.",
            formula="177/1000+(39/1000)/2=393/2000",
            proof_boundary="The source labels itself a pre-freeze draft and is not peer reviewed.",
            diagnostics=exact["candidate_record"],
        ),
        GateRow(
            id="np15lpa_04_archive_integrity",
            role="finite_reproducibility",
            readiness="finite_validated",
            claim="The deposited archive and every file named by its internal SHA-256 manifest passed integrity checks.",
            formula="505/505 manifest entries; bad=0; missing=0",
            proof_boundary="Byte integrity proves provenance, not mathematical correctness.",
            diagnostics=exact["archive_integrity_audit"],
        ),
        GateRow(
            id="np15lpa_05_exact_row",
            role="exact_arithmetic",
            readiness="available_exact",
            claim="The candidate parameter row and Platt-Trudgian height margin are exact.",
            formula=(
                "Lambda_row=393/2000; T_PT-X/2=350479773/2; "
                "1/5-393/2000=7/2000"
            ),
            proof_boundary="Citation arithmetic alone does not establish the analytic certificate legs.",
        ),
        GateRow(
            id="np15lpa_06_portable_rerun",
            role="finite_reproducibility",
            readiness="finite_validated",
            claim="Twenty-two direct portable certificate invocations exited zero on this machine.",
            formula="22 pass, 0 fail, measured 271.18 seconds",
            proof_boundary="A finite software rerun is not peer review or a formal proof of source correspondence.",
            diagnostics=exact["portable_rerun_audit"],
        ),
        GateRow(
            id="np15lpa_07_compiled_boundary",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="Three FLINT/Arb producer packages were integrity-checked but not live-rerun here.",
            formula=json.dumps(exact["unreproduced_boundary"], sort_keys=True),
            proof_boundary="Do not describe this local audit as a full independent reproduction.",
            diagnostics=exact["unreproduced_boundary"],
        ),
        GateRow(
            id="np15lpa_08_review_boundary",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="The candidate remains an unrefereed preprint despite its unusually strong certificate package.",
            formula="reproduced candidate record != established peer-reviewed record",
            proof_boundary="Expert theorem-to-code review and publication status remain open.",
        ),
        GateRow(
            id="np15lpa_09_interval_effect",
            role="exact_consequence",
            readiness="available_exact",
            claim="Acceptance of the candidate would narrow only the positive Newman interval.",
            formula="[0,1/5] -> [0,393/2000]",
            proof_boundary="No qualitative step toward excluding every Lambda>0 follows.",
            diagnostics=exact["interval_effect"],
        ),
        GateRow(
            id="np15lpa_10_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The proof programme still needs a theorem excluding a positive finite multiple-zero boundary.",
            formula="exclude every Lambda in (0,393/2000] if the candidate is accepted",
            proof_boundary="The interval shrink does not prove strict Laguerre positivity, Lambda<=0, or RH.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_polymath15_lambda01965_provenance_audit",
        "date": "2026-07-17",
        "status": (
            "published Lambda<=0.2 provenance plus a reproduced but unrefereed "
            "Lambda<=0.1965 candidate audit; not a proof of Lambda<=0 or RH"
        ),
        "proof_boundary": (
            "This artifact distinguishes the peer-reviewed interval 0<=Lambda<=0.2 "
            "from the June 2026 candidate Lambda<=0.1965. It verifies exact row "
            "arithmetic, archive integrity, and 22 portable certificate invocations. "
            "It does not claim a full FLINT/Arb rerun, independent theorem-to-code "
            "correspondence, peer review, exclusion of a positive Newman boundary, "
            "Lambda<=0, or RH."
        ),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "sources": [
            POLYMATH_URL,
            PLATT_TRUDGIAN_URL,
            RODGERS_TAO_URL,
            ZENODO_URL,
            ZENODO_API_URL,
            "outputs/jensen_window_pf_newman_positive_boundary_attainment_lemma.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    published = exact["published_state"]
    candidate = exact["candidate_record"]
    integrity = exact["archive_integrity_audit"]
    rerun = exact["portable_rerun_audit"]
    boundary = exact["unreproduced_boundary"]
    lines = [
        "# Jensen-Window PF Newman Polymath-15 Lambda 0.1965 Provenance Audit",
        "",
        "Date: 2026-07-17",
        "",
        "Status: published-bound provenance and candidate-certificate audit. The",
        "candidate is reproduced in part but remains unrefereed. This is not a",
        "proof of `Lambda <= 0` or RH.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_newman_polymath15_lambda01965_provenance_audit.json",
        "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_lambda01965_provenance_audit.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_lambda01965_provenance_audit.py",
        "```",
        "",
        "## Established Record",
        "",
        "Rodgers-Tao prove `Lambda >= 0`. Polymath15 proved `Lambda <= 0.22`,",
        "and Platt-Trudgian Theorem 1 rigorously verifies RH through the exact",
        f"height `{published['platt_trudgian_height']}`. Their Corollary 2 gives",
        "",
        "```text",
        "0 <= Lambda <= 1/5 = 0.2",
        "```",
        "",
        "This is the peer-reviewed published interval used by the corpus.",
        "",
        "## Candidate Record",
        "",
        f"The Zenodo deposit `{candidate['doi']}`, dated",
        f"`{candidate['publication_date']}`, claims",
        "",
        "```text",
        "t0=177/1000",
        "y0^2=39/1000",
        "t0+y0^2/2=393/2000=0.1965",
        f"X={candidate['x']}",
        f"T_PT-X/2={candidate['pt_height_margin']['exact']}",
        "```",
        "",
        "The PDF calls itself `DRAFT v0.3 - not for distribution`; the bundle",
        "README says `PRE-FREEZE DRAFT`, and the deposit is a preprint. The paper",
        "also discloses that the result and certificate campaign are AI-derived.",
        "",
        "## Integrity And Rerun",
        "",
        "```text",
        f"archive MD5={integrity['archive_md5_observed']}",
        f"manifest entries={integrity['manifest_entries']}",
        f"manifest bytes hashed={integrity['manifest_bytes_hashed']}",
        f"manifest bad={integrity['manifest_bad']}",
        f"manifest missing={integrity['manifest_missing']}",
        f"portable invocations={rerun['invocation_count']}",
        f"portable pass={rerun['pass_count']}",
        f"portable fail={rerun['fail_count']}",
        f"portable measured seconds={rerun['measured_seconds']}",
        "```",
        "",
        "The local pass covered both assembly lines, both binding lines, all",
        "three normalization seams and their second lines, both site-glue lines,",
        "both Euler-3 tail lines, the three block-uniform finite-box modes, the",
        "y-reduction theorem, and the criterion/error-term spines.",
        "",
        "## Nonpromotion Boundary",
        "",
        "The following compiled producer packages were hash-checked but not",
        "live-rerun on this Windows machine:",
        "",
        "```text",
        *boundary["compiled_packages_not_live_rerun"],
        "```",
        "",
        "Nor has this audit supplied an independent line-by-line review linking",
        "every implementation formula to Polymath15. The correct corpus label is",
        "therefore `reproduced unrefereed candidate`, not established record.",
        "",
        "## Consequence",
        "",
        "If accepted, the candidate changes only",
        "",
        "```text",
        "[0,1/5] -> [0,393/2000]",
        "absolute shrink=7/2000",
        "relative shrink=7/400",
        "```",
        "",
        "It does not alter the qualitative target: rule out every positive",
        "boundary collision. No `Lambda <= 0` or RH conclusion is supplied.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman Polymath-15 Lambda 0.1965 provenance audit: "
        f"{len(artifact['rows'])} rows, 1 published interval, "
        f"{len(PORTABLE_RUNS)} portable passes, 3 compiled-package boundaries, "
        "1 candidate nonpromotion gate"
    )


if __name__ == "__main__":
    main()
