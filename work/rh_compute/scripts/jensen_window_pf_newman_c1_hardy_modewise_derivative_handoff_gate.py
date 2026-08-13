#!/usr/bin/env python3
"""Build the exact shifted-Hardy derivative-to-mode handoff and source audit."""

from __future__ import annotations

import argparse
from decimal import Decimal, getcontext
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path


getcontext().prec = 70

REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_hardy_modewise_derivative_handoff_gate"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{KIND}.md"
RESUMABLE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/resumable"
    / "zeta14cubicmult_resumable.f90"
)
JOINT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_hardy_joint_correlated_error_ledger_gate.json"
)
BUILDER = Path(__file__).resolve()
CHECKER = (
    REPO_ROOT
    / "work/rh_compute/scripts"
    / "check_jensen_window_pf_newman_c1_hardy_modewise_derivative_handoff_gate.py"
)

ARXIV_ID = "2607.15310v1"
ARXIV_ABS_URL = "https://arxiv.org/abs/2607.15310"
ARXIV_PDF_URL = "https://arxiv.org/pdf/2607.15310"
UPSTREAM_REPOSITORY = "https://github.com/dml2391/Hardy-function-fastcodes"
UPSTREAM_COMMIT = "2e16dac3206b707052c3ac4cacdf3d1a2325e636"
UPSTREAM_MULTI_SOURCE_SHA256 = (
    "99283214c26ac41fbbb72ba9b058659989af5580972601a8ae2be95c11c7c5e2"
)

SOURCE_MARKERS = (
    "Numbercalc<=8 is small enough to ensure the accuracy",
    "is no hard and fast rule about this",
    "tar(i)=t+(i-numbercalc)*0.01",
    "yphasear(i)=yphase+(i-numbercalc)*0.01",
    "rsphasear(i)=rsphase+(i-numbercalc)",
    "aar(i)=a+4*(i-numbercalc)*0.01/(p*a)",
    "parameters phi1,2,3 are identical for all",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative_path(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT)).replace("\\", "/")


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def decimal_fraction(value: Fraction) -> Decimal:
    return Decimal(value.numerator) / Decimal(value.denominator)


def primitive_integer_vector(values: list[Fraction]) -> list[int]:
    common = 1
    for value in values:
        common = math.lcm(common, value.denominator)
    integers = [value.numerator * (common // value.denominator) for value in values]
    divisor = 0
    for value in integers:
        divisor = math.gcd(divisor, abs(value))
    require(divisor > 0, "zero polynomial in Gram-Schmidt")
    integers = [value // divisor for value in integers]
    if integers[-1] < 0:
        integers = [-value for value in integers]
    return integers


def inner(left: list[Fraction], right: list[Fraction]) -> Fraction:
    return sum((a * b for a, b in zip(left, right, strict=True)), Fraction(0))


def exact_orthogonal_modes(max_degree: int = 5) -> list[list[int]]:
    grid = [Fraction(index, 100) for index in range(-7, 8)]
    modes: list[list[int]] = []
    for degree in range(max_degree + 1):
        vector = [value**degree for value in grid]
        for mode in modes:
            mode_fraction = [Fraction(value) for value in mode]
            scale = inner(vector, mode_fraction) / inner(mode_fraction, mode_fraction)
            vector = [
                value - scale * basis
                for value, basis in zip(vector, mode_fraction, strict=True)
            ]
        modes.append(primitive_integer_vector(vector))
    return modes


def algebraic_factor_record(numerator: Fraction, norm_squared: int) -> dict:
    require(numerator >= 0, "algebraic factor numerator must be nonnegative")
    decimal = decimal_fraction(numerator) / Decimal(norm_squared).sqrt()
    return {
        "exact": f"({fraction_text(numerator)})/sqrt({norm_squared})",
        "numerator": fraction_text(numerator),
        "norm_squared": str(norm_squared),
        "decimal": format(decimal, ".40E"),
    }


def build_exact_mode_handoff() -> dict:
    grid = [Fraction(index, 100) for index in range(-7, 8)]
    modes = exact_orthogonal_modes()
    mode_rows = []
    for degree, mode in enumerate(modes):
        norm_squared = sum(value * value for value in mode)
        derivative_weights = {}
        mode_fraction = [Fraction(value) for value in mode]
        for derivative in range(6):
            moment = inner(mode_fraction, [value**derivative for value in grid])
            if derivative < degree:
                require(moment == 0, "lower-degree orthogonality failed")
                continue
            if moment == 0:
                continue
            numerator = abs(moment) / math.factorial(derivative)
            derivative_weights[str(derivative)] = algebraic_factor_record(
                numerator, norm_squared
            )
        remainder_numerator = sum(
            Fraction(abs(coefficient)) * abs(value) ** 6
            for coefficient, value in zip(mode, grid, strict=True)
        ) / math.factorial(6)
        mode_rows.append(
            {
                "degree": degree,
                "integer_vector": mode,
                "norm_squared": str(norm_squared),
                "derivative_weights": derivative_weights,
                "sixth_derivative_remainder_weight": algebraic_factor_record(
                    remainder_numerator, norm_squared
                ),
            }
        )

    tail_squared = sum(value**12 for value in grid) / math.factorial(6) ** 2
    tail_decimal = decimal_fraction(tail_squared).sqrt()
    return {
        "shift_grid": [fraction_text(value) for value in grid],
        "maximum_absolute_shift": "7/100",
        "orthogonal_modes": mode_rows,
        "rough_tail_weight": {
            "exact_squared": fraction_text(tail_squared),
            "decimal": format(tail_decimal, ".40E"),
        },
        "hypothesis": (
            "For r=0,...,6, epsilon is C6 on [T-7/100,T+7/100] and "
            "sup |d^r epsilon/du^r| <= D_r(T)."
        ),
        "low_mode_conclusion": (
            "A_k(T)=sum_(r=k)^5 w_(k,r)D_r(T)+w_(k,6,rem)D_6(T) "
            "bounds |<epsilon_T,q_k>| for 0<=k<=5."
        ),
        "rough_tail_conclusion": (
            "R_6(T)=sqrt(sum_j |delta_j|^12)D_6(T)/6! bounds "
            "||(I-P_5)epsilon_T||_2."
        ),
        "analytic_disk_sufficient_condition": (
            "If the fixed-branch error is holomorphic on |z-T|<=rho with "
            "rho>7/100 and sup |epsilon(z)|<=M_rho(T), Cauchy gives "
            "D_r(T)<=r!*M_rho(T)/rho^r for 0<=r<=6."
        ),
    }


def source_line_numbers(text: str) -> dict[str, list[int]]:
    lines = text.splitlines()
    locations = {}
    for marker in SOURCE_MARKERS:
        hits = [index + 1 for index, line in enumerate(lines) if marker in line]
        require(hits, f"multi-evaluator source marker missing: {marker}")
        locations[marker] = hits
    return locations


def build_source_audit() -> dict:
    source_text = RESUMABLE_SOURCE.read_text(encoding="utf-8")
    locations = source_line_numbers(source_text)
    return {
        "paper": {
            "arxiv_id": ARXIV_ID,
            "title": (
                "A computational algorithm for the Hardy function Z(t), utilising "
                "sub-sequences of generalised cubic Gauss sums"
            ),
            "authors": ["David Lewis", "Ashley Brereton"],
            "submitted": "2026-07-15",
            "pages": 66,
            "abstract_url": ARXIV_ABS_URL,
            "pdf_url": ARXIV_PDF_URL,
            "audited_locations": [
                {
                    "pages": "9-11",
                    "equations": ["24", "25", "27", "29"],
                    "finding": (
                        "epsilon_t is a prescribed relative tolerance and the truncated "
                        "integral representation is stated with an O(epsilon_t) remainder."
                    ),
                },
                {
                    "pages": "41-42",
                    "equations": ["120", "121", "122", "123"],
                    "finding": (
                        "The MGS progenitor error is reduced to max_l |epsilon_l|, but "
                        "the paper says a precise upper bound needs an estimate for that "
                        "quantity and does not pursue one for the higher-order sums."
                    ),
                },
                {
                    "page": "44",
                    "equations": ["126", "127"],
                    "finding": (
                        "The hybrid Hardy representation is explicitly described as "
                        "non-exact, with both representation and Gaussian-sum errors."
                    ),
                },
                {
                    "pages": "49-50",
                    "table": "4",
                    "finding": (
                        "Accuracy support is empirical; relative-error statistics exclude "
                        "points with |Z(t)|<=1/2 and the paper states there is no guarantee "
                        "that constituent cubic-sum errors do not combine adversely."
                    ),
                },
            ],
            "audit_conclusion": (
                "The paper supplies neither an explicit higher-order local-iteration "
                "constant nor a uniform six-derivative or complex-neighborhood bound for "
                "the fifteen-value shifted evaluator."
            ),
        },
        "code": {
            "repository": UPSTREAM_REPOSITORY,
            "upstream_commit": UPSTREAM_COMMIT,
            "upstream_multi_source_sha256": UPSTREAM_MULTI_SOURCE_SHA256,
            "audited_derivative_path": relative_path(RESUMABLE_SOURCE),
            "audited_derivative_sha256": file_hash(RESUMABLE_SOURCE),
            "marker_lines": locations,
            "finding": (
                "The fifteen-value implementation shares central hierarchy parameters "
                "and varies nearby heights and phases on the fixed 0.01 grid. Its own "
                "accuracy comment says there is no hard-and-fast rule; it contains no "
                "machine-checkable derivative remainder contract."
            ),
        },
    }


def build_rows() -> list[dict]:
    return [
        {
            "id": "hmdh_01_primary_source",
            "role": "source_audit",
            "readiness": "source_pinned",
            "claim": "The arXiv v1 paper and pinned GPL-3.0 multi-evaluator are identified exactly.",
            "proof_boundary": "Bibliographic and source provenance only.",
        },
        {
            "id": "hmdh_02_paper_tolerance",
            "role": "source_audit",
            "readiness": "asymptotic_only",
            "claim": "The paper uses a prescribed epsilon_t tolerance and O(epsilon_t) relative-error notation.",
            "proof_boundary": "Big-O notation without the needed explicit shifted absolute constant is not a certificate.",
        },
        {
            "id": "hmdh_03_mgs_reduction",
            "role": "source_audit",
            "readiness": "conditional_bound",
            "claim": "Equations 120-123 reduce MGS error to a largest local iteration error.",
            "proof_boundary": "The multiplier is useful only after max_l |epsilon_l| is explicitly bounded.",
        },
        {
            "id": "hmdh_04_missing_iteration_constant",
            "role": "nonpromotion_guard",
            "readiness": "open",
            "claim": "The paper explicitly leaves the precise higher-order local-iteration bound unformulated.",
            "proof_boundary": "Sample accuracy cannot replace this missing theorem.",
        },
        {
            "id": "hmdh_05_hybrid_error",
            "role": "source_audit",
            "readiness": "nonexact_representation",
            "claim": "The hybrid Hardy formula has representation and Gaussian-sum errors of comparable stated order.",
            "proof_boundary": "Neither component has the joint shifted absolute enclosure required here.",
        },
        {
            "id": "hmdh_06_sample_boundary",
            "role": "nonpromotion_guard",
            "readiness": "diagnostic_only",
            "claim": "The paper's final accuracy evidence is sample-based and its relative statistics omit the near-zero regime.",
            "proof_boundary": "Finite samples do not prove physical-height persistence or correlation.",
        },
        {
            "id": "hmdh_07_multi_shift_structure",
            "role": "source_audit",
            "readiness": "structurally_favourable",
            "claim": "The code shares a central hierarchy across fifteen nearby values, so a fixed-branch joint expansion is a viable target.",
            "proof_boundary": "The code comment itself supplies no rigorous error or branch-domain bound.",
        },
        {
            "id": "hmdh_08_exact_grid_modes",
            "role": "exact_lemma",
            "readiness": "available_exact",
            "claim": "Integer Gram polynomials give an exact orthogonal basis through degree five on the fifteen-point 0.01 grid.",
            "proof_boundary": "Finite linear algebra only.",
        },
        {
            "id": "hmdh_09_derivative_handoff",
            "role": "conditional_theorem",
            "readiness": "available_conditional",
            "claim": "A uniform C6 error bound contracts explicitly into all six low modes and the orthogonal rough tail.",
            "proof_boundary": "The derivative bounds D_0 through D_6 remain unproved for the evaluator.",
        },
        {
            "id": "hmdh_10_analytic_disk",
            "role": "alternative_handoff",
            "readiness": "available_conditional",
            "claim": "One fixed-branch complex-disk error bound implies the seven real derivative bounds by Cauchy's estimate.",
            "proof_boundary": "Holomorphy, branch stability, and the disk supremum must all be certified.",
        },
        {
            "id": "hmdh_11_transition_guard",
            "role": "nonpromotion_guard",
            "readiness": "mandatory",
            "claim": "If the physical interval crosses a combinatorial or analytic transition, it must be partitioned and jump terms bounded before mode contraction.",
            "proof_boundary": "Observed smoothness at two calibration heights does not establish transition-free physical behavior.",
        },
        {
            "id": "hmdh_12_open_target",
            "role": "handoff",
            "readiness": "not_ready_to_apply",
            "claim": "Prove either the fixed-branch C6 bounds, one complex-disk bound, or a transition-aware piecewise substitute at physical height.",
            "proof_boundary": "No physical carrier observation is admissible before one route is closed.",
        },
    ]


def build_artifact() -> dict:
    joint = json.loads(JOINT_RESULT.read_text(encoding="utf-8"))
    require(joint.get("kind", "").endswith("joint_correlated_error_ledger_gate"), "joint error ledger kind mismatch")
    sources = {
        "joint_error_ledger": {
            "path": relative_path(JOINT_RESULT),
            "sha256": file_hash(JOINT_RESULT),
        },
        "resumable_multi_evaluator": {
            "path": relative_path(RESUMABLE_SOURCE),
            "sha256": file_hash(RESUMABLE_SOURCE),
        },
        "builder": {"path": relative_path(BUILDER), "sha256": file_hash(BUILDER)},
        "checker": {"path": relative_path(CHECKER), "sha256": file_hash(CHECKER)},
    }
    rows = build_rows()
    return {
        "kind": KIND,
        "status": "exact_conditional_handoff_with_open_external_error_theorem",
        "source_audit": build_source_audit(),
        "exact_mode_handoff": build_exact_mode_handoff(),
        "rows": rows,
        "sources": sources,
        "open_theorem": (
            "At the physical central height, certify a branch-stable C6 absolute-error "
            "envelope on [T-0.07,T+0.07], a stronger holomorphic-disk envelope, or an "
            "explicit transition-aware piecewise replacement; then apply the exact mode "
            "weights before contracting with the twenty-four physical kernels."
        ),
        "proof_boundary": (
            "This gate proves exact finite derivative-to-mode algebra and audits the "
            "published/source-code error claims. It does not prove any evaluator derivative "
            "bound, complex-disk bound, transition exclusion, external absolute error "
            "constant, physical carrier value, interval enclosure, retained observation, "
            "determinant sign or bound, complete-current inequality, all-q transport, "
            "contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
        "success": (
            "built Newman C1 Hardy modewise derivative handoff gate: 12 rows, 6 exact "
            "low modes, 7 derivative orders, 1 exact rough-tail constant, 4 primary-paper "
            "locations, 7 source markers, 2 conditional routes and 1 open external error theorem"
        ),
    }


def render_note(artifact: dict) -> str:
    handoff = artifact["exact_mode_handoff"]
    paper = artifact["source_audit"]["paper"]
    code = artifact["source_audit"]["code"]
    lines = [
        "# Newman C1 Hardy Modewise Derivative Handoff Gate",
        "",
        "Date: 2026-08-05",
        "",
        f"Status: `{artifact['status']}`; not a proof of the external error theorem, RH, or a prize-level conclusion.",
        "",
        "## Source audit",
        "",
        f"Primary paper: [{paper['arxiv_id']}]({paper['abstract_url']}) by "
        + ", ".join(paper["authors"])
        + ".",
        "",
        "The paper's equations (120)--(123) reduce the recursive MGS error to the largest "
        "local iteration error.  On pages 41--42 it then states that a precise upper bound "
        "would require an estimate of that term and does not formulate one for the "
        "higher-order sums.  Page 44 separately says the hybrid Hardy representation is "
        "not exact.  Table 4 is sample evidence; its relative-error statistics omit "
        "`|Z(t)|<=1/2`, and the surrounding text says adverse combination of cubic-sum "
        "errors is not guaranteed away.",
        "",
        "The pinned fifteen-value source is more structurally useful than an independent "
        "fifteen-run model: it shares central hierarchy parameters and varies the nearby "
        "heights and phases on one `0.01` grid.  But its own comment says there is no "
        "hard-and-fast rule behind the multi-value accuracy claim.  Source derivative hash: "
        f"`{code['audited_derivative_sha256']}`.",
        "",
        "## Exact handoff",
        "",
        "Put `delta_j=j/100`, `-7<=j<=7`, and let `q_k` be the normalized integer Gram "
        "polynomial of degree `k` listed in the artifact.  Assume on the whole shifted "
        "interval",
        "",
        "```text",
        "sup_(|u-T|<=7/100) |epsilon^(r)(u)| <= D_r(T),  0<=r<=6.",
        "```",
        "",
        "Taylor's theorem and exact discrete orthogonality give",
        "",
        "```text",
        "|<epsilon_T,q_k>|",
        " <= sum_(r=k)^5 w_(k,r) D_r(T) + w_(k,6,rem) D_6(T),",
        "",
        "||(I-P_5)epsilon_T||_2",
        " <= [sqrt(sum_j |delta_j|^12)/6!] D_6(T).",
        "```",
        "",
        "The weights are exact algebraic numbers `rational/sqrt(integer)`; no fitted "
        "decay law enters.  The rough-tail multiplier is "
        f"`{handoff['rough_tail_weight']['decimal']}` and its square is exactly "
        f"`{handoff['rough_tail_weight']['exact_squared']}`.",
        "",
        "| mode | nonzero center-derivative weights | sixth-order remainder weight |",
        "|---:|---|---:|",
    ]
    for row in handoff["orthogonal_modes"]:
        weights = ", ".join(
            f"D_{order}: {record['decimal']}"
            for order, record in row["derivative_weights"].items()
        )
        lines.append(
            f"| {row['degree']} | {weights} | "
            f"{row['sixth_derivative_remainder_weight']['decimal']} |"
        )
    lines.extend(
        [
            "",
            "A stronger alternative is one fixed-branch holomorphic disk certificate: if "
            "`|epsilon(z)|<=M_rho(T)` on `|z-T|<=rho` with `rho>0.07`, Cauchy's "
            "estimate supplies `D_r<=r! M_rho/rho^r`.  If any combinatorial or analytic "
            "transition crosses the physical interval, the interval must instead be "
            "partitioned and its jump terms entered as additional exact mode vectors.",
            "",
            "## Result",
            "",
            "The scalar two-height route remains retired.  The missing input is now exact: "
            "prove one branch-stable `C^6` envelope, one stronger complex-disk envelope, or "
            "a transition-aware piecewise substitute for the external fifteen-value error. "
            "Only then may the saved physical coefficient vectors be applied.",
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
