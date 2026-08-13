#!/usr/bin/env python3
"""Transport the complete local endpoint theorem through block-20 weights."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import acb, arb, ctx


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_transported_endpoint_budget_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_transported_endpoint_budget_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_block20_transported_endpoint_budget_gate.py"
WEIGHT_FIXTURE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/accumulation_weights/t1e10_block20/fixture_result.json"
WEIGHTS = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/accumulation_weights/t1e10_block20/run/block20_accumulation_weights.txt"
MAJORANT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_w4_contour_endpoint_hurwitz_majorant_gate.json"
DECOMPOSITION = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_full_gate.json"
DIRECT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_direct_kernel_gate.json"
T2_MISMATCH = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_t2_lplus1_mismatch_gate.json"
REQUESTED_ERROR_SCALE = Fraction(5, 1000)
PRECISION_DIGITS = 120


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def upper(value: arb) -> Fraction:
    return cells.bound_fraction(value.upper())


def lower(value: arb) -> Fraction:
    return cells.bound_fraction(value.lower())


def upper_abs(value: arb | acb) -> Fraction:
    return upper(abs(value))


def exact_real(payload: str) -> arb:
    return cells.fraction_ball(cells.binary128_fraction(payload))


def load_weights() -> dict[tuple[int, int, int], dict[str, Any]]:
    weights: dict[tuple[int, int, int], dict[str, Any]] = {}
    for line_number, line in enumerate(WEIGHTS.read_text(encoding="ascii").splitlines(), 1):
        fields = line.split()
        require(len(fields) == 10, f"weight field-count drift at line {line_number}")
        sum_index, output_index = int(fields[0]), int(fields[1])
        amplitude_fraction = cells.binary128_fraction(fields[3])
        require(amplitude_fraction > 0, f"nonpositive amplitude at line {line_number}")
        for branch, offset in ((1, 6), (2, 8)):
            phase = acb(exact_real(fields[offset]), exact_real(fields[offset + 1]))
            key = (sum_index, output_index, branch)
            require(key not in weights, f"duplicate weight key {key}")
            weights[key] = {
                "amplitude_fraction": amplitude_fraction,
                "amplitude": cells.fraction_ball(amplitude_fraction),
                "phase": phase,
                "phase_abs_upper": upper_abs(phase),
            }
    require(len(weights) == 212 * 15 * 2, "expanded weight roster drift")
    return weights


def oriented(value: acb, conjugate: bool) -> acb:
    return value.conjugate() if conjugate else value


def projected(weight: dict[str, Any], value: acb) -> arb:
    return weight["amplitude"] * (weight["phase"] * value).real


def output_label(index: int) -> str:
    offset = index - 8
    return f"t{offset / 100:+.2f}"


def compute() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    ctx.dps = PRECISION_DIGITS
    ctx.threads = 1
    weights = load_weights()
    recursive = native_q.load_recursive_chains()
    majorant = json.loads(MAJORANT.read_text(encoding="utf-8"))
    decomposition = json.loads(DECOMPOSITION.read_text(encoding="utf-8"))
    direct = json.loads(DIRECT.read_text(encoding="utf-8"))
    mismatch = json.loads(T2_MISMATCH.read_text(encoding="utf-8"))

    majorant_rows = {int(row["chain"]): row for row in majorant["rows"]}
    correction_rows = {int(row["chain"]): row for row in decomposition["rows"]}
    direct_rows = {int(row["chain"]): row for row in direct["rows"]}
    mismatch_rows = {int(row["chain"]): row for row in mismatch["rows"]}
    require(set(majorant_rows) == set(recursive), "recursive majorant roster mismatch")
    require(set(correction_rows) == set(recursive), "recursive correction roster mismatch")
    require(set(mismatch_rows) == set(recursive), "source/paper shift roster mismatch")
    require(set(direct_rows).isdisjoint(recursive), "direct/recursive roster overlap")
    require(set(direct_rows) | set(recursive) == set(range(1, 425)), "424-chain roster gap")

    rows: list[dict[str, Any]] = []
    for output_index in range(1, 16):
        endpoint_majorant = Fraction(0)
        direct_majorant = Fraction(0)
        endpoint_signed = arb(0)
        direct_signed = arb(0)
        source_to_paper_signed = arb(0)
        source_to_exact_signed = arb(0)

        for chain in sorted(recursive):
            header = recursive[chain]["header"]
            step = recursive[chain]["steps"][1]
            sum_index = int(header["sum_index"])
            branch = int(header["branch"])
            weight = weights[(sum_index, output_index, branch)]
            local_bound = Fraction(majorant_rows[chain]["complete_majorant_upper"])
            endpoint_majorant += (
                weight["amplitude_fraction"] * weight["phase_abs_upper"] * local_bound
            )

            correction = native_q.complex_from_record(correction_rows[chain]["exact_correction"])
            deletion = native_q.complex_from_record(mismatch_rows[chain]["deletion_correction"])
            conjugate = bool(step["conjugate"])
            oriented_correction = oriented(correction, conjugate)
            oriented_deletion = oriented(deletion, conjugate)
            endpoint_signed += projected(weight, oriented_correction)
            source_to_paper_signed += projected(weight, oriented_deletion)
            source_to_exact_signed += projected(
                weight, oriented_deletion + oriented_correction
            )

        for chain in sorted(direct_rows):
            row = direct_rows[chain]
            weight = weights[(int(row["sum_index"]), output_index, int(row["branch"]))]
            local_bound = Fraction(row["total_true_gap_abs_upper"])
            direct_majorant += (
                weight["amplitude_fraction"] * weight["phase_abs_upper"] * local_bound
            )
            direct_gap = native_q.complex_from_record(row["total_true_gap"])
            direct_signed += projected(weight, direct_gap)

        combined_majorant = endpoint_majorant + direct_majorant
        endpoint_actual_upper = upper_abs(endpoint_signed)
        direct_actual_upper = upper_abs(direct_signed)
        paper_shift_upper = upper_abs(source_to_paper_signed)
        source_to_exact_upper = upper_abs(source_to_exact_signed)
        signed_combined = endpoint_signed + direct_signed
        signed_combined_upper = upper_abs(signed_combined)
        cancellation_ratio = (
            endpoint_majorant / endpoint_actual_upper
            if endpoint_actual_upper > 0
            else Fraction(0)
        )
        rows.append(
            {
                "output_index": output_index,
                "output_label": output_label(output_index),
                "recursive_call_count": len(recursive),
                "direct_call_count": len(direct_rows),
                "transported_endpoint_majorant_upper": decimal(endpoint_majorant),
                "transported_direct_majorant_upper": decimal(direct_majorant),
                "combined_local_majorant_upper": decimal(combined_majorant),
                "endpoint_majorant_to_requested_scale_ratio": decimal(
                    endpoint_majorant / REQUESTED_ERROR_SCALE
                ),
                "endpoint_majorant_within_requested_scale": endpoint_majorant <= REQUESTED_ERROR_SCALE,
                "signed_endpoint_correction": point_balls.arb_record(endpoint_signed, 50),
                "signed_endpoint_correction_abs_upper": decimal(endpoint_actual_upper),
                "endpoint_majorant_to_signed_correction_ratio": decimal(cancellation_ratio),
                "signed_direct_true_gap": point_balls.arb_record(direct_signed, 50),
                "signed_direct_true_gap_abs_upper": decimal(direct_actual_upper),
                "signed_endpoint_plus_direct": point_balls.arb_record(signed_combined, 50),
                "signed_endpoint_plus_direct_abs_upper": decimal(signed_combined_upper),
                "signed_source_to_paper_shift": point_balls.arb_record(source_to_paper_signed, 50),
                "signed_source_to_paper_shift_abs_upper": decimal(paper_shift_upper),
                "signed_source_to_exact_net_shift": point_balls.arb_record(source_to_exact_signed, 50),
                "signed_source_to_exact_net_shift_abs_upper": decimal(source_to_exact_upper),
            }
        )

    def maximum(field: str) -> tuple[Fraction, int]:
        return max((Fraction(row[field]), int(row["output_index"])) for row in rows)

    def minimum(field: str) -> tuple[Fraction, int]:
        return min((Fraction(row[field]), int(row["output_index"])) for row in rows)

    max_endpoint = maximum("transported_endpoint_majorant_upper")
    min_endpoint = minimum("transported_endpoint_majorant_upper")
    max_ratio = maximum("endpoint_majorant_to_requested_scale_ratio")
    max_signed = maximum("signed_endpoint_correction_abs_upper")
    max_direct = maximum("transported_direct_majorant_upper")
    max_paper_shift = maximum("signed_source_to_paper_shift_abs_upper")
    max_net_shift = maximum("signed_source_to_exact_net_shift_abs_upper")
    min_cancellation = minimum("endpoint_majorant_to_signed_correction_ratio")
    aggregate = {
        "output_count": len(rows),
        "recursive_call_count": len(recursive),
        "direct_call_count": len(direct_rows),
        "total_call_count": len(recursive) + len(direct_rows),
        "requested_error_scale": decimal(REQUESTED_ERROR_SCALE),
        "outputs_with_endpoint_majorant_within_requested_scale": sum(
            bool(row["endpoint_majorant_within_requested_scale"]) for row in rows
        ),
        "minimum_transported_endpoint_majorant_upper": decimal(min_endpoint[0]),
        "minimum_transported_endpoint_majorant_witness": min_endpoint[1],
        "maximum_transported_endpoint_majorant_upper": decimal(max_endpoint[0]),
        "maximum_transported_endpoint_majorant_witness": max_endpoint[1],
        "maximum_endpoint_majorant_to_requested_scale_ratio": decimal(max_ratio[0]),
        "maximum_endpoint_majorant_to_requested_scale_ratio_witness": max_ratio[1],
        "maximum_signed_endpoint_correction_abs_upper": decimal(max_signed[0]),
        "maximum_signed_endpoint_correction_witness": max_signed[1],
        "minimum_endpoint_majorant_to_signed_correction_ratio": decimal(min_cancellation[0]),
        "minimum_endpoint_majorant_to_signed_correction_ratio_witness": min_cancellation[1],
        "maximum_transported_direct_majorant_upper": decimal(max_direct[0]),
        "maximum_transported_direct_majorant_witness": max_direct[1],
        "maximum_signed_source_to_paper_shift_abs_upper": decimal(max_paper_shift[0]),
        "maximum_signed_source_to_paper_shift_witness": max_paper_shift[1],
        "maximum_signed_source_to_exact_net_shift_abs_upper": decimal(max_net_shift[0]),
        "maximum_signed_source_to_exact_net_shift_witness": max_net_shift[1],
    }
    return rows, aggregate


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    table = "\n".join(
        f"{row['output_label']:>7}  {row['transported_endpoint_majorant_upper']:>18}  "
        f"{row['signed_endpoint_correction_abs_upper']:>18}  "
        f"{row['endpoint_majorant_to_requested_scale_ratio']:>12}"
        for row in artifact["rows"]
    )
    return f"""# Hardy block-20 transported endpoint budget gate

Date: 2026-08-07

Status: rigorous finite source-weight propagation audit; not a proof of the complete evaluator or RH

## Exact transport map

The source-equivalent observer records the binary128 amplitude and the two
complex phase factors for all `212 x 15` outer accumulation rows.  Chain
`2j-1` uses branch one and chain `2j` branch two.  A recursive native-q
correction `Delta q` reaches output row `i` as

```text
amp_(j,i) Re[c_(j,i,b) O_(j,b)(Delta q)],             (1)
```

where `O` is either identity or complex conjugation exactly as recorded by
the recurrence.  The subtract-one operation is affine and does not alter an
error.  The fifty `MIT=1` calls have no endpoint correction and retain their
separate direct-kernel audit.

## Absolute budget

Transporting the complete selector-gap-uniform local endpoint majorants gives

```text
minimum over 15 outputs  <= {aggregate['minimum_transported_endpoint_majorant_upper']}
maximum over 15 outputs  <= {aggregate['maximum_transported_endpoint_majorant_upper']}
requested input scale       {aggregate['requested_error_scale']}
maximum budget/scale ratio   {aggregate['maximum_endpoint_majorant_to_requested_scale_ratio']}
outputs closing absolutely   {aggregate['outputs_with_endpoint_majorant_within_requested_scale']} / 15
```

The direct-kernel contribution is negligible on this scale:

```text
maximum transported direct majorant <= {aggregate['maximum_transported_direct_majorant_upper']}.
```

```text
 output       endpoint bound   observed signed size   bound/0.005
{table}
```

The observed signed correction is not used to shrink the theorem.  It is an
independent diagnostic showing how much cancellation the finite roster has;
the smallest bound-to-signed ratio over the 15 outputs is
`{aggregate['minimum_endpoint_majorant_to_signed_correction_ratio']}`.  The
failed Hurwitz budget can be repaired either by a sharper absolute local
theorem or by a signed cross-call estimate.  Because only a factor 1.716 is
needed, the next gate first tests whether retaining the positive homotopy
decay closes the absolute route.

The deterministic source-only `1/(L+1)` deletion is transported separately;
it is not relabelled as uncertainty.  This finite block audit does not bound
coefficient transport, Legendre truncation, cross-block accumulation, the
outer Hardy representation, or any height-uniform RH criterion.
"""


def main() -> int:
    dependencies = (WEIGHT_FIXTURE, WEIGHTS, MAJORANT, DECOMPOSITION, DIRECT, T2_MISMATCH, CHECKER)
    for path in dependencies:
        require(path.is_file(), f"missing transported-budget dependency: {path}")
    fixture = json.loads(WEIGHT_FIXTURE.read_text(encoding="utf-8"))
    require(fixture["weights"]["row_count"] == 3180, "weight fixture row-count drift")
    require(fixture["equivalence"]["displayed_hardy_values_exact"], "weight fixture lost source equivalence")
    rows, aggregate = compute()
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_transported_endpoint_budget_gate",
        "status": "rigorous_424_call_source_weight_transport_absolute_endpoint_budget_exceeds_requested_scale",
        "scope": {
            "block": 20,
            "output_count": 15,
            "recursive_calls": 374,
            "direct_calls": 50,
            "weight_rows": 3180,
            "arb_precision_decimal_digits": PRECISION_DIGITS,
        },
        "transport_identity": "amp_(j,i)*Re(c_(j,i,b)*O_(j,b)(Delta_q))",
        "rows": rows,
        "aggregate": aggregate,
        "dependencies": {
            path.stem: {"path": relative(path), "sha256": file_hash(path)}
            for path in dependencies[:-1]
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": (
            "First retain the positive linear, quadratic, and cubic homotopy decay discarded by the "
            "Hurwitz moments. If that sharper absolute theorem does not close, derive a signed cross-call "
            "phase estimate. Keep coefficient transport, Legendre truncation, cross-block accumulation, "
            "and the outer Hardy remainder as separate budgets."
        ),
        "proof_boundary": (
            "Rigorous finite block-20 transport through exact source-observed binary128 amplitudes and phases. "
            "The local endpoint majorants are rigorous, but the observed signed cancellation is diagnostic and "
            "is not promoted to a uniform theorem. No coefficient-neighborhood, Legendre, cross-block, outer "
            "Hardy, Lambda <= 0, PF-infinity, RH, or prize-level conclusion follows."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print(
        "built transported endpoint budget: "
        f"max {aggregate['maximum_transported_endpoint_majorant_upper']}, "
        f"ratio {aggregate['maximum_endpoint_majorant_to_requested_scale_ratio']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
