#!/usr/bin/env python3
"""Build the Hardy MGS local-defect obligation ledger and bypass handoff."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_hardy_local_defect_obligation_ledger_gate"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{KIND}.md"
SOURCE = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/resumable"
    / "zeta14cubicmult_resumable.f90"
)
HANDOFF_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_hardy_modewise_derivative_handoff_gate.json"
)
BUILDER = Path(__file__).resolve()
CHECKER = (
    REPO_ROOT
    / "work/rh_compute/scripts"
    / "check_jensen_window_pf_newman_c1_hardy_local_defect_obligation_ledger_gate.py"
)

ARXIV_ID = "2607.15310v1"
ARXIV_ABS_URL = "https://arxiv.org/abs/2607.15310"
ARXIV_PDF_URL = "https://arxiv.org/pdf/2607.15310"
UPSTREAM_REPOSITORY = "https://github.com/dml2391/Hardy-function-fastcodes"
UPSTREAM_COMMIT = "2e16dac3206b707052c3ac4cacdf3d1a2325e636"

SOURCE_MARKERS = (
    "ip=3",
    "I have found ip~20 to be adequate.",
    "ignore psi(2,z) correction here, which is included in MAPLE version",
    "call erf(z,6,cr1)",
    "tol=1e-7",
    "do i=jbot,ip1",
    "do i=L(nit),L(nit)-ip1+1,-1",
    "ecor=0.0",
    "t5=t5+(SUM0+SUM1)",
    "qq=conjg(t1+t2+t4)+t3+t5",
    "call q(m,k,ip,qq)",
    "csum=c1*csum+qq",
    "if (abs(smalle).gt.0.001)",
    "csum=csum+exp(fn)/(denn+i*(rminorcor-i*rminorcor2))",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative_path(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT)).replace("\\", "/")


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_marker_locations() -> dict[str, list[int]]:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    locations: dict[str, list[int]] = {}
    for marker in SOURCE_MARKERS:
        folded = marker.casefold()
        hits = [
            index + 1
            for index, line in enumerate(lines)
            if folded in line.casefold()
        ]
        require(hits, f"source marker missing: {marker}")
        locations[marker] = hits
    require(
        len(locations["call erf(z,6,cr1)"]) == 2,
        "expected two six-term complex-erf call sites",
    )
    return locations


def paper_audit() -> dict:
    return {
        "arxiv_id": ARXIV_ID,
        "abstract_url": ARXIV_ABS_URL,
        "pdf_url": ARXIV_PDF_URL,
        "audited_locations": [
            {
                "pages": "22-23",
                "equations": ["64", "67a", "69"],
                "finding": (
                    "The Euler-Maclaurin/Fourier reduction is exact through equation "
                    "(69); the paper explicitly says estimation starts with its "
                    "constituent integrals."
                ),
            },
            {
                "pages": "23-25",
                "equations": ["71", "72", "74", "75", "76", "77", "78"],
                "finding": (
                    "The saddle contour is replaced by a quadratic phase, the two "
                    "vertical endpoint phases discard order-three-and-higher terms, "
                    "and their quadratic exponentials are truncated after first order."
                ),
            },
            {
                "pages": "25-26",
                "equations": ["79", "80", "81"],
                "finding": (
                    "W1 combines exact digamma sums with endpoint/saddle corrections. "
                    "Equation (81) is formally summed over all relevant n, while the "
                    "paper uses the empirical cutoff P=3 in practice."
                ),
            },
            {
                "pages": "28-30",
                "equations": ["89", "93", "94", "96", "97"],
                "finding": (
                    "W2, W3 or W4, W5, and the endpoint half-sum complete the local "
                    "reciprocity approximation."
                ),
            },
            {
                "pages": "37-38",
                "equations": ["MGS step h", "MGS step i"],
                "finding": (
                    "The kernel is evaluated directly and the recurrence is iterated "
                    "back to the progenitor."
                ),
            },
            {
                "pages": "40-42",
                "equations": ["120", "121", "122", "123"],
                "finding": (
                    "Equation (120) defines each local defect as an exact finite-sum "
                    "difference; equations (121)-(123) propagate those defects, but the "
                    "paper leaves their higher-order maximum unbounded."
                ),
            },
        ],
        "source_conclusion": (
            "The missing theorem is not a single floating-point tolerance. It is a "
            "uniform bound for the exact local defect after analytic truncation, "
            "endpoint cutoff, hierarchy modelling, and numerical evaluation are all "
            "accounted for."
        ),
    }


def source_component_map(locations: dict[str, list[int]]) -> list[dict]:
    return [
        {
            "paper_component": "W1",
            "paper_equation": "81",
            "fortran_component": "t5",
            "source_lines": sorted(
                set(
                    locations["do i=jbot,ip1"]
                    + locations["do i=L(nit),L(nit)-ip1+1,-1"]
                    + locations["t5=t5+(SUM0+SUM1)"]
                )
            ),
            "finding": "The two ip-truncated endpoint loops and digamma base are paper W1.",
        },
        {
            "paper_component": "W2",
            "paper_equation": "89",
            "fortran_component": "t2",
            "source_lines": locations["call erf(z,6,cr1)"][:1],
            "finding": "The first six-term complex-erf correction is paper W2.",
        },
        {
            "paper_component": "W3_or_W4",
            "paper_equation": "93_or_94",
            "fortran_component": "t4",
            "source_lines": locations["call erf(z,6,cr1)"][1:],
            "finding": "The sign-selected second complex-erf correction is paper W3 or W4.",
        },
        {
            "paper_component": "W5",
            "paper_equation": "97",
            "fortran_component": "t1",
            "source_lines": locations[
                "ignore psi(2,z) correction here, which is included in MAPLE version"
            ],
            "finding": "The non-saddle digamma difference is paper W5.",
        },
        {
            "paper_component": "endpoint_half_sum",
            "paper_equation": "96",
            "fortran_component": "t3",
            "source_lines": locations["qq=conjg(t1+t2+t4)+t3+t5"],
            "finding": "t3 is the exact half-sum before floating-point evaluation error.",
        },
        {
            "paper_component": "CW_assembly_and_recurrence",
            "paper_equation": "120",
            "fortran_component": "qq_and_csum",
            "source_lines": sorted(
                set(
                    locations["qq=conjg(t1+t2+t4)+t3+t5"]
                    + locations["call q(m,k,ip,qq)"]
                    + locations["csum=c1*csum+qq"]
                )
            ),
            "finding": "q assembles CW plus the endpoint term and mgausssum applies it once per level.",
        },
    ]


def error_components() -> list[dict]:
    return [
        {
            "id": "eld_w1_saddle_phase",
            "scope": "local_MGS",
            "origin": "equations_72_to_74",
            "remainder": "order-three-and-higher phase on the diagonal saddle contour",
            "readiness": "open_analytic_bound",
            "certificate": "A contour-domain majorant retaining every Phi_q term.",
        },
        {
            "id": "eld_w1_vertical_phase",
            "scope": "local_MGS",
            "origin": "equations_75_to_78",
            "remainder": "discarded higher phase plus exponential Taylor remainder on both vertical legs",
            "readiness": "open_analytic_bound",
            "certificate": "Two explicit integral remainders with endpoint-uniform constants.",
        },
        {
            "id": "eld_w1_index_tail",
            "scope": "local_MGS",
            "origin": "equation_81_and_ip_3",
            "remainder": "all lower and upper endpoint correction indices omitted beyond ip",
            "readiness": "open_tail_bound",
            "certificate": "A summable majorant uniform in Phi, L, and every admissible branch.",
        },
        {
            "id": "eld_w1_saddle_location",
            "scope": "local_MGS",
            "origin": "series_saddle_or_Newton_tol_1e_minus_7",
            "remainder": "phase error induced by approximate saddle locations",
            "readiness": "open_root_enclosure",
            "certificate": "Interval root isolation and a residual-to-location theorem, including multiplicity separation.",
        },
        {
            "id": "eld_digamma_evaluation",
            "scope": "local_MGS_numeric",
            "origin": "PSI_truncated_series_and_asymptotic",
            "remainder": "finite psi series, asymptotic tail, decimal constants, and pole distance",
            "readiness": "open_special_function_enclosure",
            "certificate": "Directed-rounding digamma balls on every source argument domain.",
        },
        {
            "id": "eld_complex_erf_evaluation",
            "scope": "local_MGS_numeric",
            "origin": "ERF_with_N_6",
            "remainder": "piecewise six-term power, Taylor, or asymptotic approximation",
            "readiness": "open_special_function_enclosure",
            "certificate": "Region-by-region analytic remainder bounds and selector margins.",
        },
        {
            "id": "eld_other_W_analytic",
            "scope": "local_MGS",
            "origin": "W2_W3_W4_endpoint_approximations",
            "remainder": "non-W1 endpoint contour linearization and sign-selected special case",
            "readiness": "open_analytic_bound",
            "certificate": "Explicit bounds for equations (84)-(94), not only W1.",
        },
        {
            "id": "eld_hierarchy_scheme",
            "scope": "local_MGS_model",
            "origin": "scheme_phifactors_and_smalle_0_001",
            "remainder": "truncated transformed coefficients and order-raising decision",
            "readiness": "open_model_remainder",
            "certificate": "An exact transformed-level specification plus a tail bound beyond m1.",
        },
        {
            "id": "eld_kernel_model",
            "scope": "local_MGS_model",
            "origin": "finite_kernel_with_amplitude_denominator",
            "remainder": "difference between the finite source kernel model and the exact transformed child quantity",
            "readiness": "open_model_identity",
            "certificate": "A formal level adapter proving what finite quantity the denominator-weighted kernel represents.",
        },
        {
            "id": "eld_roundoff_accumulation",
            "scope": "local_MGS_numeric",
            "origin": "quad_precision_straight_line_recurrence",
            "remainder": "elementary, phase, summation, and recurrence roundoff",
            "readiness": "feasible_interval_enclosure",
            "certificate": "Directed-rounding complex balls with a recorded operation graph.",
        },
        {
            "id": "eld_hardy_representation",
            "scope": "whole_Hardy_evaluator",
            "origin": "hybrid_formula_truncation_EM_and_transition",
            "remainder": "representation error outside the MGS local defect",
            "readiness": "open_separate_bound",
            "certificate": "Explicit Euler-Maclaurin, lower-cutoff, transition, and Riemann-Siegel-tail bounds.",
        },
        {
            "id": "eld_shared_multi_shift",
            "scope": "whole_shifted_batch",
            "origin": "central_hierarchy_reused_on_0_01_grid",
            "remainder": "shift dependence suppressed by the shared hierarchy parameterization",
            "readiness": "open_C6_or_piecewise_bound",
            "certificate": "Differentiate a fixed branch through order six or enclose its analytic surrogate on a disk.",
        },
        {
            "id": "eld_transition_set",
            "scope": "whole_shifted_batch",
            "origin": "floor_round_sign_and_piecewise_special_function_selectors",
            "remainder": "branch jumps or selector changes across the physical interval/disk",
            "readiness": "mandatory_exclusion_or_partition",
            "certificate": "Strict interval margins or an explicit partition with jump-mode vectors.",
        },
    ]


RationalComplex = tuple[Fraction, Fraction]


def c_add(left: RationalComplex, right: RationalComplex) -> RationalComplex:
    return left[0] + right[0], left[1] + right[1]


def c_mul(left: RationalComplex, right: RationalComplex) -> RationalComplex:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def c_record(value: RationalComplex) -> dict[str, str]:
    return {
        "real": f"{value[0].numerator}/{value[0].denominator}",
        "imag": f"{value[1].numerator}/{value[1].denominator}",
    }


def exact_accumulation_fixture() -> dict:
    multipliers: list[RationalComplex] = [
        (Fraction(2, 3), Fraction(1, 3)),
        (Fraction(-1, 5), Fraction(2, 5)),
        (Fraction(3, 7), Fraction(-1, 7)),
    ]
    defects: list[RationalComplex] = [
        (Fraction(1, 11), Fraction(-2, 11)),
        (Fraction(3, 13), Fraction(1, 13)),
        (Fraction(-1, 17), Fraction(4, 17)),
        (Fraction(2, 19), Fraction(-3, 19)),
    ]
    recursive = defects[-1]
    for index in range(len(multipliers) - 1, -1, -1):
        recursive = c_add(defects[index], c_mul(multipliers[index], recursive))

    expanded: RationalComplex = (Fraction(0), Fraction(0))
    product: RationalComplex = (Fraction(1), Fraction(0))
    terms = []
    for index, defect in enumerate(defects):
        term = c_mul(product, defect)
        terms.append(c_record(term))
        expanded = c_add(expanded, term)
        if index < len(multipliers):
            product = c_mul(product, multipliers[index])
    require(recursive == expanded, "exact recurrence expansion failed")

    shell_level = 2
    shell_error: RationalComplex = (Fraction(0), Fraction(0))
    product = (Fraction(1), Fraction(0))
    for index in range(shell_level):
        shell_error = c_add(shell_error, c_mul(product, defects[index]))
        product = c_mul(product, multipliers[index])

    return {
        "multipliers": [c_record(value) for value in multipliers],
        "defects": [c_record(value) for value in defects],
        "expanded_terms": terms,
        "recursive_root_error": c_record(recursive),
        "expanded_root_error": c_record(expanded),
        "direct_shell_level": shell_level,
        "direct_shell_root_error": c_record(shell_error),
        "identity": "e_k=epsilon_k+a_k*e_(k+1)",
        "expanded_identity": "e_0=sum_(k=0)^r (prod_(j=0)^(k-1) a_j) epsilon_k",
        "shell_identity": (
            "If level m is evaluated directly, epsilon_m through epsilon_r are absent "
            "and e_0=sum_(k=0)^(m-1)(prod_(j=0)^(k-1)a_j)epsilon_k."
        ),
    }


def exact_defect_routes() -> dict:
    return {
        "paper_level_defect": (
            "epsilon_l=S_l-[a_l S_(l+1)+q_l], exactly as equation (120); "
            "both S terms are finite sums once a transformed-level adapter is fixed."
        ),
        "direct_defect_certificate": (
            "Enclose the two finite level sums and the source recurrence term directly. "
            "This captures W1, W2-W5, coefficient-transform, and numerical defects without "
            "requiring their cancellations to be bounded separately."
        ),
        "exact_shell_bypass": (
            "Directly evaluate one or more levels above the kernel and begin recurrence "
            "there. The skipped local defects vanish from the exact accumulation identity."
        ),
        "cost_guard": (
            "No complexity claim is made until the physical chain roster records the "
            "direct-shell lengths. A shell is admissible only if its measured finite cost "
            "fits the bounded evaluator budget."
        ),
        "adapter_obligation": (
            "The source's denominator-weighted kernel and raised-order coefficients must "
            "first be identified with a precise finite target; otherwise a direct defect "
            "could compare the wrong two quantities."
        ),
        "accumulation_fixture": exact_accumulation_fixture(),
    }


def build_rows() -> list[dict]:
    return [
        {"id": "held_01_source_pin", "role": "source_audit", "readiness": "source_pinned", "claim": "The paper version, upstream commit, and local derivative hash are fixed.", "proof_boundary": "Provenance only."},
        {"id": "held_02_exact_until_69", "role": "exact_reduction", "readiness": "available_exact", "claim": "The paper states that its reduction is exact through equation (69).", "proof_boundary": "No later contour approximation is promoted."},
        {"id": "held_03_w1_anatomy", "role": "source_audit", "readiness": "component_mapped", "claim": "Equation (81) is split into saddle, vertical-leg, endpoint-index, root, and digamma obligations.", "proof_boundary": "The component bounds remain open."},
        {"id": "held_04_fortran_mapping", "role": "source_audit", "readiness": "component_mapped", "claim": "Fortran t5,t2,t4,t1,t3 map to W1,W2,W3/4,W5, and the endpoint half-sum.", "proof_boundary": "Variable-name mapping only."},
        {"id": "held_05_cutoff_conflict", "role": "nonpromotion_guard", "readiness": "open_tail_bound", "claim": "The executable fixes ip=3 although its q routine comments say ip about 20 was found adequate.", "proof_boundary": "Neither empirical choice is a theorem."},
        {"id": "held_06_special_functions", "role": "nonpromotion_guard", "readiness": "open_enclosures", "claim": "Six-term complex-erf and truncated digamma routines require explicit remainder and selector margins.", "proof_boundary": "Quad precision does not imply rigorous special-function accuracy."},
        {"id": "held_07_scheme_kernel", "role": "nonpromotion_guard", "readiness": "open_adapter", "claim": "Hierarchy truncation and the denominator-weighted kernel need a precise transformed-level target.", "proof_boundary": "Direct summation is exact only after that target is fixed."},
        {"id": "held_08_exact_defect", "role": "exact_lemma", "readiness": "available_exact", "claim": "Equation (120) permits direct finite certification of each local recurrence defect.", "proof_boundary": "No physical defect enclosure has yet been computed."},
        {"id": "held_09_accumulation", "role": "exact_lemma", "readiness": "available_exact", "claim": "Local defects propagate by an exact affine recurrence and weighted finite sum.", "proof_boundary": "Useful only with certified local inputs."},
        {"id": "held_10_exact_shell", "role": "algorithmic_bypass", "readiness": "available_conditional", "claim": "Direct evaluation of a level above the kernel removes every skipped descendant defect from the accumulation.", "proof_boundary": "Physical shell lengths and cost remain unmeasured."},
        {"id": "held_11_termwise_route", "role": "theorem_route", "readiness": "open", "claim": "One route bounds each W1 contour and index-tail remainder with explicit constants.", "proof_boundary": "No such constants are supplied here."},
        {"id": "held_12_real_C6_route", "role": "theorem_route", "readiness": "open", "claim": "A fixed real branch can be differentiated and enclosed through order six.", "proof_boundary": "All branch and selector margins must be strict."},
        {"id": "held_13_disk_route", "role": "theorem_route", "readiness": "open", "claim": "A holomorphic mathematical surrogate can be enclosed on a complex disk and transferred by Cauchy.", "proof_boundary": "The raw source ABS/floor selectors are not themselves a holomorphic program."},
        {"id": "held_14_representation_split", "role": "scope_guard", "readiness": "mandatory", "claim": "MGS local defects and the outer Hardy representation remainder are separate budgets.", "proof_boundary": "Closing one does not close the other."},
        {"id": "held_15_next_experiment", "role": "handoff", "readiness": "ready_to_implement", "claim": "Instrument one accepted low-height chain, record every level, and compare termwise W1 versus direct-defect versus exact-shell routes.", "proof_boundary": "A low-height scout remains diagnostic until converted to uniform intervals."},
        {"id": "held_16_physical_block", "role": "nonpromotion_guard", "readiness": "not_ready_to_apply", "claim": "Physical kernel contraction remains blocked pending the external C6/disk or piecewise error theorem.", "proof_boundary": "No physical carrier value, sign theorem, Lambda<=0, or RH conclusion follows."},
    ]


def build_artifact() -> dict:
    handoff = json.loads(HANDOFF_RESULT.read_text(encoding="utf-8"))
    require(
        handoff.get("kind", "").endswith("hardy_modewise_derivative_handoff_gate"),
        "modewise derivative handoff kind mismatch",
    )
    locations = source_marker_locations()
    components = error_components()
    rows = build_rows()
    return {
        "kind": KIND,
        "status": "exact_defect_reduction_with_componentwise_open_bounds",
        "paper_audit": paper_audit(),
        "source_audit": {
            "repository": UPSTREAM_REPOSITORY,
            "upstream_commit": UPSTREAM_COMMIT,
            "path": relative_path(SOURCE),
            "sha256": file_hash(SOURCE),
            "marker_lines": locations,
            "paper_to_fortran_map": source_component_map(locations),
            "critical_findings": [
                "The executable sets ip=3 while the q-routine comment says ip about 20 was found adequate.",
                "The source explicitly omits a psi(2,z) correction present in a Maple version.",
                "The source sets ecor=0 although the surrounding comment says the correction may be significant.",
                "Complex erf is evaluated with six retained terms and saddle Newton iteration uses tol=1e-7.",
            ],
        },
        "error_components": components,
        "exact_defect_routes": exact_defect_routes(),
        "derivative_handoff": {
            "source_path": relative_path(HANDOFF_RESULT),
            "source_sha256": file_hash(HANDOFF_RESULT),
            "required_output": (
                "A branch-stable bound for D_0,...,D_6 on [T-0.07,T+0.07], "
                "one holomorphic disk bound with radius greater than 0.07, or a "
                "transition-aware piecewise replacement."
            ),
        },
        "rows": rows,
        "sources": {
            "builder": {"path": relative_path(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative_path(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_executable_target": (
            "Add read-only chain telemetry to one accepted low-height fixture. For every "
            "level record L, m1, Phi, x, fracL, ip, q components, direct finite level "
            "values, and the exact local defect at increasing precision. Then test direct "
            "evaluation of the first parent shell. Do not alter the accepted evaluator "
            "output path during this scout."
        ),
        "proof_boundary": (
            "This gate proves the paper-to-source component map, the exact affine defect "
            "accumulation identity, and the exact-shell bypass implication. It does not "
            "prove any W1 tail or contour bound, special-function enclosure, transformed-"
            "level adapter, evaluator C6 or disk theorem, transition exclusion, outer Hardy "
            "representation bound, physical carrier value, interval kernel approximation, "
            "determinant sign or bound, complete-current inequality, all-q transport, "
            "contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
        "success": (
            "built Newman C1 Hardy local-defect obligation ledger gate: 16 rows, 6 "
            "paper-to-source mappings, 13 error components, 2 exact defect identities, "
            "1 exact-shell bypass, 2 conditional derivative routes and 0 external error bounds"
        ),
    }


def render_note(artifact: dict) -> str:
    source = artifact["source_audit"]
    lines = [
        "# Newman C1 Hardy Local-Defect Obligation Ledger Gate",
        "",
        "Date: 2026-08-06",
        "",
        f"Status: `{artifact['status']}`; this is an exact reduction and source audit, not a proof of the missing external error theorem.",
        "",
        "## What equation (81) contains",
        "",
        f"Primary source: [{ARXIV_ID}]({ARXIV_ABS_URL}).  The paper says its reduction is exact through equation (69).  Approximation begins when each finite integral is replaced by the contour/saddle formulas (71)--(78).  Equation (81) then packages the resulting `W1` correction: an exact digamma base plus diagonal-saddle, lower-endpoint, and upper-endpoint approximations.",
        "",
        "The paper says the remaining endpoint-index terms can in principle be summed over all `n`, but uses `P=3` because they appear to decay quickly.  Pages 41--42 then identify `W1` as the usual largest local-error source and leave the higher-order bound unformulated.",
        "",
        "## Paper-to-source map",
        "",
        "| paper term | equation | Fortran term | role |",
        "|---|---:|---|---|",
    ]
    for row in source["paper_to_fortran_map"]:
        lines.append(
            f"| {row['paper_component']} | {row['paper_equation']} | "
            f"`{row['fortran_component']}` | {row['finding']} |"
        )
    lines.extend(
        [
            "",
            f"Pinned derivative source: `{source['path']}` (`{source['sha256']}`).",
            "",
            "The executable fixes `ip=3`, while the same `q` routine says `ip~20` was found adequate.  It also uses six-term complex-erf approximations, a truncated `PSI` implementation, saddle Newton tolerance `1e-7`, explicitly omits a `psi(2,z)` correction present in a Maple version, and sets the commented `ecor` correction to zero.  These are obligations, not accusations of numerical failure: each needs a bound before the result can be called rigorous.",
            "",
            "## Error ledger",
            "",
            "| id | scope | unresolved remainder | required certificate |",
            "|---|---|---|---|",
        ]
    )
    for row in artifact["error_components"]:
        lines.append(
            f"| `{row['id']}` | {row['scope']} | {row['remainder']} | {row['certificate']} |"
        )
    routes = artifact["exact_defect_routes"]
    lines.extend(
        [
            "",
            "## Exact-defect route",
            "",
            "Equation (120) already provides a second route that avoids proving every internal cancellation separately:",
            "",
            "```text",
            "epsilon_l = S_l - [a_l S_(l+1) + q_l].",
            "",
            "e_l = epsilon_l + a_l e_(l+1),",
            "",
            "e_0 = sum_k (product_(j<k) a_j) epsilon_k.",
            "```",
            "",
            routes["direct_defect_certificate"],
            "",
            "There is also an algorithmic bypass.  If level `m` is evaluated directly, recurrence begins from that exact level and every descendant defect `epsilon_m,epsilon_(m+1),...` disappears from the root-error identity.  In particular, directly summing the first parent above the kernel removes the very defect the paper says is usually largest.  The physical chain lengths must be measured before claiming this is cheap.",
            "",
            "The immediate caveat is structural: " + routes["adapter_obligation"],
            "",
            "## Next target",
            "",
            artifact["next_executable_target"],
            "",
            "That scout should compare three routes on identical saved chains: termwise `W1` bounds, direct finite local-defect enclosure, and one-parent exact-shell replacement.  The winner is the route that yields a useful uniform constant with the least new machinery, not the route that merely matches samples best.",
            "",
            artifact["proof_boundary"],
            "",
            "```text",
            artifact["success"],
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    atomic_write(args.out, json.dumps(artifact, indent=2) + "\n")
    atomic_write(args.note, render_note(artifact))
    print(artifact["success"])


if __name__ == "__main__":
    main()
