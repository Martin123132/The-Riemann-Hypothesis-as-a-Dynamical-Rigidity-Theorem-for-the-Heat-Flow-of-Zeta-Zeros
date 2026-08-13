#!/usr/bin/env python3
"""Build a fail-closed dependency ledger from source outputs to Hardy Z."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path
import time
from typing import Any

from pypdf import PdfReader


REPO_ROOT = Path(__file__).resolve().parents[3]
SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
PAPER = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_Brereton_2607.15310.pdf"
FIXTURE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/accumulation_weights/t1e10_crossblock/fixture_result.json"
CHECKPOINT = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/accumulation_weights/t1e10_crossblock/run/checkpoint.jsonl"
INPUTS = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/accumulation_weights/t1e10_crossblock/run/inputs3.nml"
RUN_OUTPUT = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/accumulation_weights/t1e10_crossblock/run/run.output.txt"
ATLAS = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_structural_atlas_gate.json"
INNER_BOUND = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_rigorous_arb_source_q_dependency_cut_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_outer_hardy_representation_dependency_ledger_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_outer_hardy_representation_dependency_ledger_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

EXPECTED_OUTPUTS = 15
EXPECTED_BLOCKS = (20, 35)
EXPECTED_CALLS = 6784
EXPECTED_WEIGHTS = 50880


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def load_last_jsonl(path: Path) -> dict[str, Any]:
    last = ""
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                last = line
    require(bool(last), "checkpoint is empty")
    # Preserve the serialized binary128 decimal payload instead of first
    # rounding it through Python's binary64 float parser.
    return json.loads(last, parse_float=str)


SOURCE_ASSERTIONS = (
    {
        "id": "source_pi_and_phase_constants",
        "lines": [323, 335],
        "tokens": ["p=4*ATAN(C)", "tpp=2*p", "tpm=-tpp"],
    },
    {
        "id": "fifteen_point_grid_and_shifted_parameters",
        "lines": [489, 602],
        "tokens": [
            "nchalf=2*numbercalc-1",
            "tar(i)=t+(i-numbercalc)*0.01",
            "rsphasear(i)=rsphase+(i-numbercalc)*(log(0.25*a)-1/(48*t*t))*0.01",
            "aar(i)=a+4*(i-numbercalc)*0.01/(p*a)",
        ],
    },
    {
        "id": "source_hybrid_split",
        "lines": [649, 662],
        "tokens": ["ZP(t)+RSremainder(t).", "zsum", "rszsum"],
    },
    {
        "id": "block_zero_transition_and_direct_terms",
        "lines": [668, 695],
        "tokens": ["call start(a,t6,yphase,transit,M2)", "call ter(ial,1/aar(i),yphasear(i),tar(i),term)", "zsum(i)=zsum(i)*sqrt(8/aar(i))"],
    },
    {
        "id": "cubic_blocks_one_through_totblock",
        "lines": [775, 840],
        "tokens": ["do iblock=rh_block_start,totblock", "call alphasum(MTM,rae,pcsum)", "zsum=zsum+zsum2"],
    },
    {
        "id": "lower_rs_sum_and_leading_endpoint_correction",
        "lines": [911, 1025],
        "tokens": ["NC=aint(CE,dp1)", "CALL pari_calc", "rszsum(i)=2*rszsum(i)+((-1)**(RSN-1))*CE*C"],
    },
    {
        "id": "final_source_addition",
        "lines": [1035, 1058],
        "tokens": ["zsum(i)=zsum(i)+rszsum(i)", "rh_stage_complete"],
    },
    {
        "id": "fixed_64_point_transition_quadrature",
        "lines": [1345, 1432],
        "tokens": ["real(dp)      :: abb(64),wee(64),ri", "call gauleg(bot,top,abb,wee)", "do i=1,64"],
    },
    {
        "id": "direct_or_recursive_cubic_sums",
        "lines": [1437, 1542],
        "tokens": ["if (MTM.lt.K) then", "call computarrs", "call mgausssum", "pcsum(i)=amp(i)*(real(c1*pcgsum(1)+c2*pcgsum(2)))"],
    },
    {
        "id": "pari_lower_rs_arithmetic",
        "lines": [2564, 2667],
        "tokens": ["u =  gmul(t, glog", "u =  gcos", "rzsum(j) = rzsum(j) + rtodbl(w)"],
    },
    {
        "id": "central_exact_formula_phases_before_rounding",
        "lines": [2970, 3035],
        "tokens": ["glngamma(u,prec)", "rsphase_pari = gmod(theta_pari, v)", "a_pari      = gsqrt(v ,prec)"],
    },
)


def audit_source() -> list[dict[str, Any]]:
    source = SOURCE.read_text(encoding="utf-8", errors="replace")
    lines = source.splitlines()
    rows: list[dict[str, Any]] = []
    for assertion in SOURCE_ASSERTIONS:
        start, end = assertion["lines"]
        excerpt = "\n".join(lines[start - 1 : end])
        missing = [token for token in assertion["tokens"] if token not in excerpt]
        require(not missing, f"source assertion {assertion['id']} missing tokens: {missing}")
        rows.append(
            {
                "id": assertion["id"],
                "line_start": start,
                "line_end": end,
                "required_tokens": assertion["tokens"],
                "missing_token_count": 0,
                "excerpt_sha256": hashlib.sha256(excerpt.encode("utf-8")).hexdigest(),
            }
        )
    return rows


def audit_paper() -> dict[str, Any]:
    reader = PdfReader(PAPER)
    page_text = {
        "44": " ".join((reader.pages[43].extract_text() or "").split()),
        "45": " ".join((reader.pages[44].extract_text() or "").split()),
    }
    combined = page_text["44"] + " " + page_text["45"]
    tokens = (
        "(126)",
        "(127)",
        "hybrid representation (126-127)",
        "is not exact",
        "relative error",
        "further errors of a similar size",
    )
    missing = [token for token in tokens if token not in combined]
    require(not missing, f"paper audit missing tokens: {missing}")
    return {
        "pdf_pages": [44, 45],
        "equations": [126, 127],
        "required_tokens": list(tokens),
        "missing_token_count": 0,
        "paper_declares_hybrid_not_exact": True,
        "paper_supplies_explicit_uniform_constant_for_all_big_o_terms_here": False,
        "page_text_sha256": hashlib.sha256(combined.encode("utf-8")).hexdigest(),
    }


def dependency_ledger(inner: dict[str, Any]) -> list[dict[str, Any]]:
    bound = inner["aggregate"]["maximum_source_q_free_complete_majorant_upper"]
    return [
        {
            "id": "hardy_z_definition",
            "kind": "normalization",
            "state": "exact",
            "statement": "Z(t)=exp(i*theta(t))*zeta(1/2+i*t), with theta(t)=Im(log Gamma(1/4+i*t/2))-t*log(pi)/2.",
            "scope": "Exact mathematical target definition for real t.",
        },
        {
            "id": "xi_hardy_normalization",
            "kind": "normalization",
            "state": "exact",
            "statement": "Xi(t)=-((t^2+1/4)/2)*abs(pi^(-s/2)*Gamma(s/2))*Z(t), s=1/2+i*t.",
            "scope": "Exact xi-to-Hardy normalization under the stated xi convention.",
        },
        {
            "id": "critical_line_conjugation_and_reality",
            "kind": "conjugation",
            "state": "exact",
            "statement": "The functional equation and conjugation symmetry make the exact Hardy target real for real t.",
            "scope": "Exact target identity; it does not identify the source estimator with that target.",
        },
        {
            "id": "source_observer_equivalence",
            "kind": "source_arithmetic",
            "state": "exact",
            "statement": "The blocks-20--35 observer preserves the saved checkpoint state and all fifteen displayed source values exactly.",
            "scope": "Exact source reproduction only, not mathematical accuracy of those values.",
        },
        {
            "id": "central_phase_scale_rounding",
            "kind": "source_arithmetic",
            "state": "open",
            "statement": "PARI evaluates central theta and a before conversion to binary128, but no end-to-end interval encloses their conversion and downstream use.",
            "required_upgrade": "Enclose central PARI-to-binary128 conversion and every downstream phase/amplitude operation, or bypass the source values.",
        },
        {
            "id": "shifted_theta_and_a_linearization",
            "kind": "truncation",
            "state": "open",
            "statement": "The fourteen noncentral outputs use first-order updates for rsphase and a rather than independent exact evaluations.",
            "required_upgrade": "Bound the Taylor remainders uniformly on |delta|<=0.07 or recompute exact interval phases and scales at every target.",
        },
        {
            "id": "hybrid_equations_126_127",
            "kind": "normalization",
            "state": "open",
            "statement": "The paper explicitly classifies the hybrid representation as non-exact and leaves relative O(epsilon_t) and comparable Gaussian-sum errors without a usable constant here.",
            "required_upgrade": "Derive a constant-bearing equality plus remainder from an exact Hardy or Riemann-Siegel identity.",
        },
        {
            "id": "omitted_H_t_normalizer",
            "kind": "normalization",
            "state": "open",
            "statement": "Equation (127) carries H(t), while the source accumulation does not expose a certified H(t) correction column.",
            "required_upgrade": "Prove and interval-evaluate the normalization factor and its interaction with every retained term.",
        },
        {
            "id": "lower_alpha_euler_maclaurin_and_exponential_tail",
            "kind": "truncation",
            "state": "open",
            "statement": "The exact integral representation's Euler-Maclaurin and lower-alpha truncation remainders have no imported enclosure in the saved certificate.",
            "required_upgrade": "State and prove explicit remainder bounds in the actual source parameter domain.",
        },
        {
            "id": "block_zero_transition_quadrature",
            "kind": "omitted_block",
            "state": "open",
            "statement": "Block zero uses a fixed 64-point Gauss-Legendre transition quadrature with no certified quadrature and contour remainder in the current ledger.",
            "required_upgrade": "Enclose the transition integral, contour truncation, and 64-point quadrature error.",
        },
        {
            "id": "block_zero_direct_terms",
            "kind": "omitted_block",
            "state": "open",
            "statement": "The block-zero direct ter calls and their scaling are outside the blocks-20--35 weight fixture and corrected-model bound.",
            "required_upgrade": "Interval-evaluate the complete block-zero roster and source rounding at all fifteen targets.",
        },
        {
            "id": "blocks_1_19",
            "kind": "omitted_block",
            "state": "open",
            "statement": "Cubic blocks 1--19 are absent from the current cross-block fixture and internal error ledger.",
            "required_upgrade": "Extend the exact observer, corrected formula, local enclosures, and transport through blocks 1--19.",
        },
        {
            "id": "blocks_20_35_corrected_internal_model",
            "kind": "source_arithmetic",
            "state": "enclosed",
            "statement": f"The corrected internal blocks-20--35 cell-plus-tail column is bounded by {bound} across the saved outputs.",
            "scope": "Finite saved-height internal corrected-model discrepancy only; it is not an outer Hardy error.",
        },
        {
            "id": "upper_alpha_portcullis_cutoff",
            "kind": "truncation",
            "state": "open",
            "statement": "The high-alpha cutoff and its matching to the lower Riemann-Siegel cutoff have no complete constant-bearing remainder enclosure.",
            "required_upgrade": "Prove exact partition coverage and bound all cutoff, floor, and discarded-tail terms.",
        },
        {
            "id": "lower_rs_finite_sum_arithmetic",
            "kind": "source_arithmetic",
            "state": "open",
            "statement": "The n<=621 Riemann-Siegel sum is source-reproduced but not interval-enclosed through phase conversion, rtodbl, binary128 accumulation, and output addition.",
            "required_upgrade": "Re-evaluate the finite lower sum with rigorous phases and outward-rounded arithmetic.",
        },
        {
            "id": "omitted_rs_corrections_and_remainder",
            "kind": "omitted_block",
            "state": "open",
            "statement": "The source visibly adds the leading endpoint correction only; higher correction terms and the final classical Riemann-Siegel remainder are not enclosed.",
            "required_upgrade": "Choose an explicit correction order and prove a numerical remainder bound valid at every target height.",
        },
        {
            "id": "final_source_accumulation_rounding",
            "kind": "source_arithmetic",
            "state": "open",
            "statement": "Roundoff for block zero, blocks 1--19, the lower RS sum, and the final zsum+rszsum addition is outside the selected internal columns.",
            "required_upgrade": "Provide a complete operation-level interval ledger or bypass the legacy accumulation with rigorous evaluation.",
        },
        {
            "id": "assembled_outer_hardy_remainder",
            "kind": "truncation",
            "state": "open",
            "statement": "No theorem currently assembles every normalization, truncation, omitted-block, and arithmetic column into |Z_source-Z_exact|.",
            "required_upgrade": "Sum nonoverlapping explicit columns into a complete outer Hardy remainder and verify it at all fifteen outputs.",
        },
    ]


def render_note(artifact: dict[str, Any]) -> str:
    counts = artifact["aggregate"]["dependency_state_counts"]
    return f"""# Outer-Hardy representation dependency ledger

Date: 2026-08-09

Status: diagnostic validated; not a proof; Hardy-Z and Xi certification remain blocked

The fifteen saved source outputs have been traced from their exact target
definitions through the hybrid Fortran path.  The audit distinguishes exact
mathematical identities, finite enclosed internal columns, and open outer
columns.  It contains {counts['exact']} exact rows, {counts['enclosed']}
enclosed internal row, and {counts['open']} open rows.

The decisive boundary is upstream of the recurrence work.  The source follows
the paper's hybrid equations (126)--(127), and the paper explicitly says that
representation is not exact.  Its relative `O(epsilon_t)` and comparable
Gaussian-sum errors are not supplied here with constants that can be summed
into a Hardy-Z remainder.  Therefore the input `et=0.005` is a requested
accuracy parameter, not a proved error theorem.

The current bound

```text
{artifact['critical_findings']['current_internal_bound_upper']}
```

remains valid for the finite saved-height corrected internal blocks-20--35
column.  It is deliberately not relabelled as `|Z_source-Z_exact|`.

The missing outer ledger includes shifted-phase and scale Taylor remainders,
the `H(t)` normalization, lower-alpha and high-alpha truncations, block zero,
blocks 1--19, the finite lower Riemann-Siegel arithmetic, omitted higher
Riemann-Siegel corrections and final remainder, and complete source rounding.
Every one of the fifteen output rows is therefore marked
`hardy_z_certified=false` and `xi_certified=false`.

The next finite falsification gate is an independent Arb Hardy-Z calibration
at the same fifteen heights.  That finite calibration is not evidence for RH;
it tests whether the saved estimator is even inside its requested tolerance.
The proof route is separate: start from an exact classical Riemann-Siegel
identity (or an exact Hardy integral), partition it exactly, and attach an
explicit constant-bearing remainder to every transformation.
"""


def main() -> int:
    started = time.time()
    source_audit = audit_source()
    paper_audit = audit_paper()
    fixture = load_json(FIXTURE)
    atlas = load_json(ATLAS)
    inner = load_json(INNER_BOUND)
    checkpoint = load_last_jsonl(CHECKPOINT)

    require(fixture["equivalence"]["displayed_hardy_values_exact"], "source observer equivalence failed")
    require(int(fixture["weights"]["row_count"]) == EXPECTED_WEIGHTS, "weight roster drift")
    require((int(fixture["weights"]["first_block"]), int(fixture["weights"]["last_block"])) == EXPECTED_BLOCKS, "fixture block roster drift")
    require(int(atlas["aggregate"]["call_count"]) == EXPECTED_CALLS, "atlas call roster drift")
    require(int(atlas["aggregate"]["recursive_call_count"]) == 1414, "recursive roster drift")
    require(int(atlas["aggregate"]["direct_call_count"]) == 5370, "direct roster drift")
    require(int(checkpoint["stage"]) == 2, "final checkpoint stage drift")
    require(int(checkpoint["numbercalc"]) == 8 and int(checkpoint["nchalf"]) == EXPECTED_OUTPUTS, "output roster drift")
    require(int(checkpoint["totblock"]) == 35 and int(checkpoint["nc"]) == 621, "outer source metadata drift")
    require(len(checkpoint["zsum"]) == EXPECTED_OUTPUTS, "final source value roster drift")
    require(len(inner["output_rows"]) == EXPECTED_OUTPUTS, "inner-bound output roster drift")

    dependencies = dependency_ledger(inner)
    states = {state: sum(row["state"] == state for row in dependencies) for state in ("exact", "enclosed", "open")}
    kinds = sorted({row["kind"] for row in dependencies})
    require(set(kinds) >= {"normalization", "conjugation", "truncation", "omitted_block", "source_arithmetic"}, "required dependency kinds missing")
    require(states["open"] > 0, "fail-closed audit unexpectedly found no open outer columns")

    open_ids = [row["id"] for row in dependencies if row["state"] == "open"]
    base_t = Decimal(str(checkpoint["t"]))
    output_rows: list[dict[str, Any]] = []
    for index, (source_value, inner_row) in enumerate(zip(checkpoint["zsum"], inner["output_rows"]), start=1):
        offset = Decimal(index - 8) / Decimal(100)
        expected_label = f"t{offset:+.2f}"
        require(inner_row["output_label"] == expected_label, f"output label drift at {index}")
        target_t = base_t + offset
        output_rows.append(
            {
                "output_index": index,
                "output_label": expected_label,
                "target_t": format(target_t, "f"),
                "source_estimator_value": str(source_value),
                "current_internal_corrected_model_upper": inner_row["source_q_free_complete_majorant_upper"],
                "exact_hardy_target": "Z(t_j)=exp(i*theta(t_j))*zeta(1/2+i*t_j)",
                "exact_xi_target": "Xi(t_j)=-((t_j^2+1/4)/2)*abs(pi^(-s_j/2)*Gamma(s_j/2))*Z(t_j)",
                "open_dependency_ids": open_ids,
                "hardy_z_certified": False,
                "xi_certified": False,
                "requested_tolerance_certified": False,
            }
        )

    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_outer_hardy_representation_dependency_ledger_gate",
        "status": "outer_hardy_dependency_audit_validated_with_explicit_open_remainder_columns",
        "passed": True,
        "first_gate": {
            "name": "complete_outer_hardy_remainder_present",
            "satisfied": False,
            "hardy_certification_claim_permitted": False,
            "reason": "At least one open outer column remains, including the assembled Hardy remainder.",
        },
        "scope": {
            "height_center": "1e10",
            "output_count": EXPECTED_OUTPUTS,
            "offset_min": "-0.07",
            "offset_max": "+0.07",
            "source_blocks_in_internal_fixture": [20, 35],
            "source_calls_in_internal_fixture": EXPECTED_CALLS,
            "outer_weight_rows": EXPECTED_WEIGHTS,
        },
        "critical_findings": {
            "paper_representation_exact": False,
            "paper_big_o_terms_have_imported_usable_constants": False,
            "source_parameter_et": str(checkpoint["et"]),
            "source_parameter_et_is_proved_hardy_error_bound": False,
            "source_shifted_theta_exactly_recomputed": False,
            "source_shifted_a_exactly_recomputed": False,
            "source_visibly_retains_only_leading_rs_endpoint_correction": True,
            "current_internal_bound_upper": inner["aggregate"]["maximum_source_q_free_complete_majorant_upper"],
            "current_internal_bound_is_complete_hardy_error": False,
        },
        "paper_audit": paper_audit,
        "source_audit": source_audit,
        "dependency_ledger": dependencies,
        "output_rows": output_rows,
        "aggregate": {
            "dependency_count": len(dependencies),
            "dependency_state_counts": states,
            "dependency_kinds": kinds,
            "outputs_hardy_z_certified": sum(row["hardy_z_certified"] for row in output_rows),
            "outputs_xi_certified": sum(row["xi_certified"] for row in output_rows),
            "outputs_with_explicit_open_dependencies": sum(bool(row["open_dependency_ids"]) for row in output_rows),
            "open_dependency_count": len(open_ids),
        },
        "next_route": {
            "finite_falsification": "Independently enclose Hardy Z at all fifteen saved heights with Arb and compare against the source estimator. This calibrates the finite run only and is not evidence for RH by height exhaustion.",
            "proof_route": "Start from an exact classical Riemann-Siegel or exact Hardy identity, make the lower/upper partition exact, and prove constant-bearing remainders for every hybrid transformation before importing the corrected cubic machinery.",
            "priority": "Run the independent fifteen-point Arb calibration first; then choose the exact classical representation and build its explicit outer remainder ledger.",
        },
        "runtime": {"elapsed_seconds": round(time.time() - started, 3), "workers": 1},
        "dependencies": {
            "accepted_source": {"path": relative(SOURCE), "sha256": file_hash(SOURCE)},
            "paper": {"path": relative(PAPER), "sha256": file_hash(PAPER)},
            "observer_fixture": {"path": relative(FIXTURE), "sha256": file_hash(FIXTURE)},
            "final_checkpoint": {"path": relative(CHECKPOINT), "sha256": file_hash(CHECKPOINT)},
            "fixture_inputs": {"path": relative(INPUTS), "sha256": file_hash(INPUTS)},
            "fixture_output": {"path": relative(RUN_OUTPUT), "sha256": file_hash(RUN_OUTPUT)},
            "crossblock_atlas": {"path": relative(ATLAS), "sha256": file_hash(ATLAS)},
            "internal_bound": {"path": relative(INNER_BOUND), "sha256": file_hash(INNER_BOUND)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": "Independently enclose the exact Hardy Z targets at the fifteen saved heights, then derive a complete constant-bearing outer representation remainder from an exact identity.",
        "proof_boundary": "Dependency classification and fail-closed diagnostic only. It proves no complete source-to-Hardy error, no height-uniform theorem, no Lambda<=0 statement, no PF-infinity statement, no RH, and no prize-level conclusion.",
    }
    require(artifact["aggregate"]["outputs_hardy_z_certified"] == 0, "Hardy outputs were accidentally promoted")
    require(artifact["aggregate"]["outputs_xi_certified"] == 0, "Xi outputs were accidentally promoted")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "validated outer-Hardy dependency audit: "
        f"15 outputs, exact={states['exact']}, enclosed={states['enclosed']}, open={states['open']}, Hardy-certified=0"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
