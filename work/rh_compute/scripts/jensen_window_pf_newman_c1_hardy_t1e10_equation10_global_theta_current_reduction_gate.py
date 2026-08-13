#!/usr/bin/env python3
"""Certify the equation-(10) pair current and global incomplete-theta reduction."""

from __future__ import annotations

import hashlib
import json
import os
from decimal import Decimal
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from pypdf import PdfReader
import sympy as sp


PAPER = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_Brereton_2607.15310.pdf"
SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/telemetry/zeta14cubicmult_telemetry.f90"
CHECKPOINT = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/accumulation_weights/t1e10_crossblock/run/checkpoint.jsonl"
SOURCE_LEDGER = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_source_aligned_hybrid_error_ledger_gate.json"
BLOCK_ZERO = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_block_zero_source_formula_identity_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation10_global_theta_current_reduction_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation10_global_theta_current_reduction_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

EXPECTED_BLOCKS = 36
EXPECTED_OUTPUTS = 15
M2 = 159577
FINAL_ALPHA = 5122421
BLOCK_ZERO_ALPHA_COUNT = 101

PAPER_PAGES = (5, 6, 18, 21, 45)
PAPER_TOKENS = (
    "starting point of this work is the exact formulation",
    "analysis is analogous for the sin",
    "Dropping this subdominant term from (10)",
    "(62)",
    "hybrid representation (126-127)",
    "is not exact",
)

SOURCE_ASSERTIONS = (
    {
        "id": "block_zero_roster_start",
        "lines": (1358, 1385),
        "tokens": ("subroutine start(a,t6,yphase,zsum,M2)", "M2=floor(a,dp1)+1", "M2=floor(a,dp1)+2"),
    },
    {
        "id": "block_zero_direct_alpha_loop",
        "lines": (677, 701),
        "tokens": ("call start(a,t6,yphase,transit,M2)", "N1=M2+ib", "do ial=M2,N1,2", "call ter(ial,1/aar(i),yphasear(i),tar(i),term)"),
    },
    {
        "id": "paired_collection_partition",
        "lines": (778, 850),
        "tokens": ("MTM=(MT-1)/2", "rae=RN1+2.0+real(MT,dp)", "RN4=2.0*(real(MT,dp)+1)*(real(jsum-1,dp))", "RN1=RN1startofblock + (2*(real(MT,dp)+1.0)*real(aenums,dp))"),
    },
    {
        "id": "retained_equation62_gaussian_formula",
        "lines": (1471, 1565),
        "tokens": ("amp(i)=sqrt(8/(aar(i)*s))", "if (MTM.lt.K) then", "call computarrs", "call mgausssum", "pcsum(i)=amp(i)*(real(c1*pcgsum(1)+c2*pcgsum(2)))"),
    },
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def paper_audit() -> dict[str, Any]:
    reader = PdfReader(PAPER)
    page_text = {
        str(page): " ".join((reader.pages[page - 1].extract_text() or "").split())
        for page in PAPER_PAGES
    }
    combined = " ".join(page_text[str(page)] for page in PAPER_PAGES)
    missing = [token for token in PAPER_TOKENS if token not in combined]
    require(not missing, f"paper audit missing tokens: {missing}")
    return {
        "pdf_pages": list(PAPER_PAGES),
        "equations": [7, 9, 10, 50, 51, 62, 126, 127],
        "required_tokens": list(PAPER_TOKENS),
        "missing_token_count": 0,
        "exact_kummer_integral_precedes_asymptotic_grouping": True,
        "equation10_contains_alpha_cosine_and_i_j_sine_pair": True,
        "paper_drops_sine_term_before_equation51": True,
        "equation62_and_hybrid_are_not_exact": True,
        "page_text_sha256": hashlib.sha256(combined.encode("utf-8")).hexdigest(),
    }


def source_audit() -> list[dict[str, Any]]:
    lines = SOURCE.read_text(encoding="utf-8", errors="replace").splitlines()
    rows: list[dict[str, Any]] = []
    for assertion in SOURCE_ASSERTIONS:
        start, end = assertion["lines"]
        excerpt = "\n".join(lines[start - 1 : end])
        missing = [token for token in assertion["tokens"] if token not in excerpt]
        require(not missing, f"source audit {assertion['id']} missing tokens: {missing}")
        rows.append(
            {
                "id": assertion["id"],
                "line_start": start,
                "line_end": end,
                "required_tokens": list(assertion["tokens"]),
                "missing_token_count": 0,
                "excerpt_sha256": hashlib.sha256(excerpt.encode("utf-8")).hexdigest(),
            }
        )
    return rows


def load_stage_zero() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with CHECKPOINT.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                row = json.loads(line, parse_float=str)
                if int(row["stage"]) == 0:
                    rows.append(row)
    require(len(rows) == EXPECTED_BLOCKS, "stage-zero block count drift")
    require([int(row["unit"]) for row in rows] == list(range(EXPECTED_BLOCKS)), "stage-zero order drift")
    return rows


def symbolic_certificates() -> dict[str, Any]:
    alpha, j, x, y = sp.symbols("alpha j x y", real=True, nonzero=True)
    u = sp.pi * alpha * j * x / 2
    common = sp.exp(sp.I * sp.pi * (alpha**2 + j**2) * x / 4)
    expanded_pair = (alpha + j) * sp.exp(sp.I * u) + (alpha - j) * sp.exp(-sp.I * u)
    trig_pair = 2 * (alpha * sp.cos(u) + sp.I * j * sp.sin(u))
    pair_residual = sp.simplify(sp.expand_complex(expanded_pair - trig_pair))
    require(pair_residual == 0, "equation-(10) pair identity failed")

    current_potential = common * sp.cos(u)
    full_pair = 2 * common * (alpha * sp.cos(u) + sp.I * j * sp.sin(u))
    current_residual = sp.simplify(4 * sp.diff(current_potential, alpha) / (sp.I * sp.pi * x) - full_pair)
    require(current_residual == 0, "pivot-current derivative identity failed")

    theta_term = sp.exp(sp.I * sp.pi * alpha**2 * x / 4 + sp.I * sp.pi * alpha * y / 2)
    theta_residual = sp.simplify(2 * sp.diff(theta_term, y) / (sp.I * sp.pi) - alpha * theta_term)
    require(theta_residual == 0, "incomplete-theta derivative identity failed")

    retained = 2 * common * alpha * sp.cos(u)
    omitted = 2 * common * sp.I * j * sp.sin(u)
    completion_residual = sp.simplify(full_pair - retained - omitted)
    require(completion_residual == 0, "retained-plus-omitted completion failed")

    # A continuous derivative is not automatically a step-two finite telescope.
    witness_sum = 2 * 1 + 2 * 3
    witness_endpoint_quotient = (5**2 - 1**2) // 2
    require(witness_sum != witness_endpoint_quotient, "discrete-telescope guard witness collapsed")
    return {
        "equation10_pair_identity": "(alpha+j)e^(iu)+(alpha-j)e^(-iu)=2[alpha cos(u)+i j sin(u)]",
        "equation10_pair_polynomial_trigonometric_residual": str(pair_residual),
        "pivot_current_identity": "P(alpha,j;x)=4/(i*pi*x)*d_alpha[e^(i*pi*(alpha^2+j^2)*x/4) cos(pi*alpha*j*x/2)]",
        "pivot_current_symbolic_residual": str(current_residual),
        "retained_cosine_piece": "C=2 alpha e^(i*pi*(alpha^2+j^2)*x/4) cos(pi*alpha*j*x/2)",
        "omitted_sine_completion": "S=2 i j e^(i*pi*(alpha^2+j^2)*x/4) sin(pi*alpha*j*x/2)",
        "completion_identity": "P=C+S",
        "completion_symbolic_residual": str(completion_residual),
        "incomplete_theta_definition": "Theta_[A,B](x,y)=sum_(A<=alpha<=B, alpha odd) exp(i*pi*alpha^2*x/4+i*pi*alpha*y/2)",
        "global_theta_derivative_identity": "sum alpha exp(i*pi*alpha^2*x/4)=2/(i*pi)*d_y Theta_[A,B](x,y)|_(y=0)",
        "theta_term_symbolic_residual": str(theta_residual),
        "discrete_telescope_guard": {
            "statement": "sum of continuous pivot derivatives on a step-two lattice is not an endpoint difference without a finite Poisson, Abel, or Euler--Maclaurin formula",
            "witness_potential": "G(alpha)=alpha^2 on alpha in {1,3}",
            "sum_of_derivatives": witness_sum,
            "naive_endpoint_quotient": witness_endpoint_quotient,
            "naive_identity_holds": False,
        },
    }


def roster_certificate(rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    block_rows: list[dict[str, Any]] = [
        {
            "block": 0,
            "alpha_start": M2,
            "alpha_end": int(Decimal(rows[0]["rn1"])),
            "alpha_count": BLOCK_ZERO_ALPHA_COUNT,
            "pair_collection_count": 0,
            "gaussian_phase_sum_count": 0,
            "coverage_identity_holds": True,
        }
    ]
    require(int(Decimal(rows[0]["rn1"])) == M2 + 2 * (BLOCK_ZERO_ALPHA_COUNT - 1), "block-zero roster drift")
    pair_collections = 0
    paired_alpha_count = 0
    recurrence_count = 0
    for block in range(1, EXPECTED_BLOCKS):
        previous = int(Decimal(rows[block - 1]["rn1"]))
        endpoint = int(Decimal(rows[block]["rn1"]))
        mt = int(rows[block]["mt"])
        aenums = int(rows[block]["aenums"])
        require(mt > 0 and mt % 2 == 1, f"MT parity drift at block {block}")
        require(aenums > 0, f"empty paired block {block}")
        expected_endpoint = previous + 2 * (mt + 1) * aenums
        require(endpoint == expected_endpoint, f"RN1 recurrence failed at block {block}")
        alpha_count = (endpoint - previous) // 2
        require(alpha_count == (mt + 1) * aenums, f"alpha count failed at block {block}")
        first_pivot = previous + 2 + mt
        require(first_pivot % 2 == 0, f"pivot parity failed at block {block}")
        require(first_pivot - mt == previous + 2, f"first collection gap at block {block}")
        last_pivot = first_pivot + 2 * (mt + 1) * (aenums - 1)
        require(last_pivot + mt == endpoint, f"last collection endpoint failed at block {block}")
        block_rows.append(
            {
                "block": block,
                "previous_endpoint": previous,
                "alpha_start": previous + 2,
                "alpha_end": endpoint,
                "alpha_step": 2,
                "mt": mt,
                "odd_offsets_per_side": (mt + 1) // 2,
                "pair_collection_count": aenums,
                "gaussian_phase_sum_count": 2 * aenums,
                "first_even_pivot": first_pivot,
                "last_even_pivot": last_pivot,
                "alpha_count": alpha_count,
                "rn1_recurrence_identity_holds": True,
                "within_block_collection_coverage_has_no_gap_or_overlap": True,
            }
        )
        recurrence_count += 1
        pair_collections += aenums
        paired_alpha_count += alpha_count

    require(int(Decimal(rows[-1]["rn1"])) == FINAL_ALPHA, "final alpha endpoint drift")
    total_alpha_count = BLOCK_ZERO_ALPHA_COUNT + paired_alpha_count
    contiguous_formula_count = (FINAL_ALPHA - M2) // 2 + 1
    require(total_alpha_count == contiguous_formula_count, "global alpha roster does not close")
    return block_rows, {
        "block_count": EXPECTED_BLOCKS,
        "block_zero_direct_alpha_count": BLOCK_ZERO_ALPHA_COUNT,
        "paired_alpha_count": paired_alpha_count,
        "total_alpha_count": total_alpha_count,
        "contiguous_formula_count": contiguous_formula_count,
        "pair_collection_count": pair_collections,
        "gaussian_phase_sum_count": 2 * pair_collections,
        "rn1_recurrence_certificate_count": recurrence_count,
        "first_alpha": M2,
        "last_alpha": FINAL_ALPHA,
        "alpha_parity": "odd",
        "alpha_step": 2,
        "global_roster_has_no_gap_or_overlap": True,
    }


def render_note(artifact: dict[str, Any]) -> str:
    agg = artifact["aggregate"]
    finite = artifact["finite_atlas_context"]
    return f"""# Equation-(10) global theta-current reduction

Date: 2026-08-10

Status: exact algebraic and roster reduction validated; not a proof of the required uniform estimate

The paper's exact paired collection in equation (10) contains both

```text
alpha_E cos(pi*j*alpha_E*x/2)
```

and

```text
i*j sin(pi*j*alpha_E*x/2).
```

Before equation (51), the paper explicitly drops the second term and applies
the saddle approximation (50).  Equation (62), equations (126)--(127), and
the source Gaussian-sum path inherit that non-exact step.

For `u=pi*alpha_E*j*x/2`, the complete pair has the exact identities

```text
(alpha_E+j)e^(iu)+(alpha_E-j)e^(-iu)
  =2[alpha_E cos(u)+i*j sin(u)],

P(alpha_E,j;x)
  =4/(i*pi*x) d/d(alpha_E)
     [e^(i*pi*(alpha_E^2+j^2)*x/4) cos(pi*alpha_E*j*x/2)].
```

Thus the omitted sine term is not an unrelated correction: it completes the
retained cosine term into a pivot derivative current.

The source checkpoint proves that its 36 blocks partition the contiguous odd
roster

```text
159577, 159579, ..., 5122421
```

with no gap or overlap.  It contains {agg['total_alpha_count']:,} alpha terms:
{agg['block_zero_direct_alpha_count']} direct terms followed by
{agg['pair_collection_count']:,} paired collections representing
{agg['paired_alpha_count']:,} terms.  The two phase branches require
{agg['gaussian_phase_sum_count']:,} Gaussian phase sums.

Define the finite incomplete theta sum

```text
Theta_[A,B](x,y)
 =sum_(A<=alpha<=B, alpha odd)
    exp(i*pi*alpha^2*x/4+i*pi*alpha*y/2).
```

Then, exactly,

```text
sum alpha exp(i*pi*alpha^2*x/4)
 =2/(i*pi) d/dy Theta_[A,B](x,y)|_(y=0).
```

This is the first global coordinate that preserves every alpha and both
members of every equation-(10) pair before asymptotic approximation.  It is
therefore a viable route for the signed source-aligned error.

The finite atlas explains why this matters: its independent block triangle is
{finite['minimum_independent_block_triangle_ball']} to
{finite['maximum_independent_block_triangle_ball']}, while only
{finite['minimum_retained_fraction_ball']} to
{finite['maximum_retained_fraction_ball']} survives in the signed sum.

There is an important guard.  A continuous pivot derivative sampled on a
step-two lattice is not automatically an endpoint telescope.  The next
theorem must apply an endpoint-complete finite Poisson, Abel, or
Euler--Maclaurin transform to the incomplete theta derivative, carry the
transition and global `REM` terms, and identify its dual saddle range with the
classical `T_upper` sum.  Only then can the source saddle truncation and
Gaussian-evaluator defects be inserted.

Proof boundary: exact symbolic identities, paper/source audit, and finite
roster coverage only.  No finite-Poisson remainder, source-error bound,
height-uniform Hardy theorem, `Lambda<=0`, RH, or prize-level conclusion is
proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    source_ledger = load_json(SOURCE_LEDGER)
    block_zero = load_json(BLOCK_ZERO)
    require(source_ledger.get("passed") is True, "source ledger is not passed")
    require(block_zero.get("passed") is True, "block-zero gate is not passed")
    require(source_ledger["scope"]["diagnostic_midpoint_used"] is False, "midpoint entered source route")
    require(block_zero["aggregate"]["alpha_terms_per_output"] == BLOCK_ZERO_ALPHA_COUNT, "block-zero alpha count drift")

    symbolic = symbolic_certificates()
    paper = paper_audit()
    source = source_audit()
    block_rows, roster = roster_certificate(load_stage_zero())
    finite = {
        key: source_ledger["aggregate"][key]
        for key in (
            "minimum_independent_block_triangle_ball",
            "maximum_independent_block_triangle_ball",
            "minimum_retained_fraction_ball",
            "maximum_retained_fraction_ball",
            "minimum_cancellation_fraction_ball",
            "maximum_cancellation_fraction_ball",
        )
    }
    aggregate = {
        **roster,
        "paper_audit_missing_token_count": paper["missing_token_count"],
        "source_audit_count": len(source),
        "source_audit_missing_token_count": sum(row["missing_token_count"] for row in source),
        "symbolic_zero_residual_count": sum(
            symbolic[key] == "0"
            for key in (
                "equation10_pair_polynomial_trigonometric_residual",
                "pivot_current_symbolic_residual",
                "completion_symbolic_residual",
                "theta_term_symbolic_residual",
            )
        ),
    }
    require(aggregate["symbolic_zero_residual_count"] == 4, "symbolic aggregate drift")

    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation10_global_theta_current_reduction_gate",
        "status": "exact_equation10_pair_current_and_contiguous_incomplete_theta_derivative_reduction_complete",
        "passed": True,
        "scope": {
            "height_center": "1e10",
            "source_blocks": "0-35",
            "alpha_roster": f"{M2}..{FINAL_ALPHA} odd",
            "diagnostic_midpoint_used": False,
        },
        "paper_audit": paper,
        "source_audit": source,
        "symbolic_reduction": symbolic,
        "block_rows": block_rows,
        "aggregate": aggregate,
        "finite_atlas_context": finite,
        "route_decision": {
            "equation10_sine_term_is_part_of_exact_pair": True,
            "paper_drops_sine_term_before_equation51": True,
            "omitted_sine_term_completes_pivot_derivative_current": True,
            "global_incomplete_theta_derivative_reduction_is_exact": True,
            "source_roster_is_one_contiguous_odd_interval": True,
            "continuous_derivative_alone_is_discrete_telescope": False,
            "blockwise_absolute_route_is_viable_at_observed_scale": False,
            "endpoint_complete_finite_poisson_or_abel_route_is_structurally_viable": True,
            "height_uniform_error_bound_proved": False,
            "rh_implication": False,
        },
        "next_obligation": "Derive an endpoint-complete finite Poisson/Abel transform for the weighted incomplete-theta derivative inside the exact equation-(9) Kummer integral, identify the dual stationary range with T_upper, and expose transition, REM, endpoint, saddle-truncation, and source-evaluator remainders with explicit constants before taking absolute values.",
        "proof_boundary": "Exact symbolic pair/current identities and finite source-roster coverage only. No endpoint-complete Poisson remainder, source-aligned height-uniform error theorem, Lambda<=0, RH, or prize-level conclusion is proved.",
        "dependencies": {
            "paper": {"path": relative(PAPER), "sha256": file_hash(PAPER)},
            "source": {"path": relative(SOURCE), "sha256": file_hash(SOURCE)},
            "source_checkpoint": {"path": relative(CHECKPOINT), "sha256": file_hash(CHECKPOINT)},
            "source_aligned_ledger": {"path": relative(SOURCE_LEDGER), "sha256": file_hash(SOURCE_LEDGER)},
            "block_zero_identity": {"path": relative(BLOCK_ZERO), "sha256": file_hash(BLOCK_ZERO)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "sympy_threads": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built equation-(10) global theta-current reduction: "
        f"blocks={EXPECTED_BLOCKS}, alpha={aggregate['total_alpha_count']}, pairs={aggregate['pair_collection_count']}"
    )


if __name__ == "__main__":
    main()
