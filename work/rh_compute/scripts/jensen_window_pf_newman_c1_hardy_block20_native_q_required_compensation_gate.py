#!/usr/bin/env python3
"""Recover the native q compensation required by every block-20 recurrence."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import acb, arb, ctx


SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
TELEMETRY = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10_block20_full/enabled/chain_telemetry.jsonl"
)
SPECIAL = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_special_function_residual_gate.json"
REAL_ERFC = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_t5_real_erfc_residual_gate.json"
FULL_CORRECTED = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_t5_erfc_corrected_recurrence_defect_gate.json"
ATLAS = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_discrete_selector_atlas_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate.py"
PRECISIONS = (180, 260)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def fraction_decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def dyadic_fraction(payload: list[int]) -> Fraction:
    require(len(payload) == 2, "dyadic payload drift")
    mantissa, exponent = int(payload[0]), int(payload[1])
    return Fraction(mantissa * 2**exponent, 1) if exponent >= 0 else Fraction(mantissa, 2 ** (-exponent))


def interval_from_record(record: dict[str, Any]) -> arb:
    lower = dyadic_fraction(record["lower_dyadic"])
    upper = dyadic_fraction(record["upper_dyadic"])
    require(lower <= upper, "stored interval order drift")
    center = (lower + upper) / 2
    radius = (upper - lower) / 2
    return cells.fraction_ball(center) if radius == 0 else cells.interval_ball(center, radius)


def complex_from_record(record: dict[str, Any]) -> acb:
    return acb(interval_from_record(record["real"]), interval_from_record(record["imag"]))


def magnitude_upper(value: arb | acb) -> Fraction:
    return cells.bound_fraction(abs(value).upper())


def magnitude_lower(value: acb) -> Fraction:
    return cells.bound_fraction(abs(value).lower())


def complex_record(value: acb) -> dict[str, Any]:
    return point_balls.acb_record(value, 55)


def load_recursive_chains() -> dict[int, dict[str, Any]]:
    chains: dict[int, dict[str, Any]] = {}
    for line in TELEMETRY.read_text(encoding="utf-8").splitlines():
        record = json.loads(line)
        chain = int(record.get("chain", 0))
        if record["type"] == "chain":
            chains[chain] = {"header": record, "levels": {}, "steps": {}, "q_terms": {}}
        elif record["type"] == "level":
            chains[chain]["levels"][int(record["level"])] = record
        elif record["type"] == "recurrence":
            chains[chain]["steps"][int(record["nit"])] = record
        elif record["type"] == "q_terms":
            chains[chain]["q_terms"][int(record["nit"])] = record
    recursive = {chain: data for chain, data in chains.items() if int(data["header"]["mit"]) == 2}
    require(len(recursive) == 374, "native-q roster drift")
    return recursive


def inverse_affine(step: dict[str, Any], value: acb) -> acb:
    result = value + 1 if step["subtract_one"] else value
    return result.conjugate() if step["conjugate"] else result


def linear_orientation(step: dict[str, Any], value: acb) -> acb:
    return value.conjugate() if step["conjugate"] else value


def corrected_step(step: dict[str, Any], child: acb, q_value: acb) -> acb:
    value = point_balls.binary128_complex(step["multiplier_hex"]) * child + q_value
    if step["conjugate"]:
        value = value.conjugate()
    if step["subtract_one"]:
        value -= 1
    return value


def evaluate(
    data: dict[str, Any],
    special_record: dict[str, Any],
    erfc_record: dict[str, Any],
    dps: int,
) -> dict[str, acb]:
    ctx.dps = dps
    ctx.threads = 1
    levels = data["levels"]
    step = data["steps"][1]
    require(set(levels) == {1, 2} and set(data["steps"]) == {1}, "native-q chain shape drift")
    require(
        bool(levels[1]["conjugate"]) == bool(step["conjugate"])
        and bool(levels[1]["subtract_one"]) == bool(step["subtract_one"]),
        "parent/step orientation drift",
    )
    tpm = point_balls.binary128_ball(data["header"]["tpm_hex"])
    parent_raw = point_balls.direct_sum(levels[1], tpm)
    parent_adapted = point_balls.adapt(levels[1], parent_raw)
    child_raw = point_balls.direct_sum(levels[2], tpm)
    child_adapted = point_balls.adapt(levels[2], child_raw)
    multiplier = point_balls.binary128_complex(step["multiplier_hex"])
    source_q = point_balls.binary128_complex(step["qq_hex"])
    special_delta = complex_from_record(special_record)
    erfc_delta = complex_from_record(erfc_record)
    full_q = source_q + special_delta + erfc_delta

    required_q_direct = parent_raw - multiplier * child_adapted
    required_q_inverse = inverse_affine(step, parent_adapted) - multiplier * child_adapted
    required_compensation = required_q_direct - full_q
    full_model = corrected_step(step, child_adapted, full_q)
    full_defect = parent_adapted - full_model
    oriented_compensation = linear_orientation(step, required_compensation)
    return {
        "parent_raw": parent_raw,
        "parent_adapted": parent_adapted,
        "child_adapted": child_adapted,
        "full_q": full_q,
        "required_q_direct": required_q_direct,
        "required_q_inverse": required_q_inverse,
        "required_q_cross_gap": required_q_direct - required_q_inverse,
        "required_compensation": required_compensation,
        "full_defect": full_defect,
        "oriented_compensation": oriented_compensation,
        "orientation_identity_gap": full_defect - oriented_compensation,
    }


def sign_class(value: arb) -> str:
    lower, upper = cells.bounds(value)
    if lower > 0:
        return "positive"
    if upper < 0:
        return "negative"
    return "contains_zero"


def source_locations() -> dict[str, int]:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    tokens = {
        "multiplier": "c1=exp(pp)/(EPI4*sqrt(x))",
        "q_call": "call q(m,k,ip,qq)",
        "native_recurrence": "csum=c1*csum+qq",
        "conjugation": "csum=conjg(csum)",
        "subtract_one": "csum=csum-1.0",
        "q_combination": "qq=conjg(t1+t2+t4)+t3+t5",
    }
    locations: dict[str, int] = {}
    for name, token in tokens.items():
        matches = [number for number, line in enumerate(lines, 1) if token in line]
        require(matches, f"native-q source token missing: {token}")
        locations[name] = matches[0]
    return locations


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 native-q required-compensation gate

Date: 2026-08-06

Status: rigorous finite native-q inversion and required-compensation target; not a proof of the missing remainder

## Exact affine inversion

Let `P_raw` be the independently summed parent, `C` the source-adapted child,
`M` the saved multiplier, and `O` the saved conjugation/subtract-one map.  The
source recurrence has the form

```text
P_adapted = O(P_raw),
model(q)  = O(M C + q).
```

Therefore the unique additive q value required for exact finite equality is

```text
q_need = P_raw - M C.                                  (1)
```

The same value is independently recovered as
`O^{-1}(P_adapted)-M C`.  If

```text
q_full = q_source + delta_q_PSI/complex-ERF
                  + delta_q_intrinsic-real-erfc,
```

then the native correction still required by the finite sums is

```text
R_q = q_need - q_full.                                 (2)
```

Subtract-one cancels from the difference.  Consequently the already admitted
fully corrected recurrence defect is exactly `R_q` when the parent orientation
is direct and `conj(R_q)` when it is conjugated.

## Rigorous result

Every construction in (1)--(2) is evaluated at 180 and 260 decimal digits from
the exact saved binary128 source normalization.  All 374 precision pairs
overlap; both q-need constructions overlap; all orientation identities contain
zero; and all reconstructed defects overlap the prior fully corrected balls.

```text
minimum |R_q|                           >= {aggregate['minimum_required_compensation_magnitude_lower']}
maximum |R_q|                           <= {aggregate['maximum_required_compensation_magnitude_upper']}
median  |R_q|                           <= {aggregate['median_required_compensation_magnitude_upper']}
maximum |R_q|/|q_full|                  <= {aggregate['maximum_relative_to_full_q_upper']}
maximum q-need cross-construction gap   <= {aggregate['maximum_required_q_cross_gap_abs_upper']}
maximum oriented-defect identity gap    <= {aggregate['maximum_orientation_identity_gap_abs_upper']}
nonzero native compensation balls        = {aggregate['required_compensations_excluding_zero']} / 374
```

This converts the observed recurrence discrepancy into the exact finite
complex target that any omitted endpoint-saddle, finite-`ip`,
Euler--Maclaurin, or compensation term must supply.  It does not show that any
particular omitted source term equals `R_q`, does not assign a sign or formula
to the remainder, and does not provide a height-uniform bound.  It is not a
proof of an outer Hardy estimate, `Lambda<=0`, PF-infinity, RH, or a prize-level
theorem.
"""


def main() -> int:
    for path in (SOURCE, TELEMETRY, SPECIAL, REAL_ERFC, FULL_CORRECTED, ATLAS, CHECKER):
        require(path.is_file(), f"missing native-q dependency: {path}")
    special = json.loads(SPECIAL.read_text(encoding="utf-8"))
    real_erfc = json.loads(REAL_ERFC.read_text(encoding="utf-8"))
    full_corrected = json.loads(FULL_CORRECTED.read_text(encoding="utf-8"))
    special_records = {
        row["chain"]: row["correlated_replacement"]["total_q_delta"]
        for row in special["rows"]
    }
    erfc_records = {
        row["chain"]: row["correlated_replacement"]["total_q_delta"]
        for row in real_erfc["rows"]
    }
    prior_records = {row["chain"]: row["full_corrected_defect"] for row in full_corrected["rows"]}
    atlas = json.loads(ATLAS.read_text(encoding="utf-8"))
    expected = {row["chain"] for row in atlas["rows"] if row["mit"] == 2}
    recursive = load_recursive_chains()
    require(
        set(recursive) == set(special_records) == set(erfc_records) == set(prior_records) == expected,
        "native-q roster mismatch",
    )

    rows: list[dict[str, Any]] = []
    maxima: dict[str, tuple[Fraction, int]] = {}
    minima: dict[str, tuple[Fraction, int]] = {}
    magnitude_uppers: list[Fraction] = []
    orientation_histogram: Counter[str] = Counter()
    sign_histogram: Counter[str] = Counter()
    prior_overlap_count = 0

    def observe_max(name: str, value: Fraction, chain: int) -> None:
        if name not in maxima or value > maxima[name][0]:
            maxima[name] = (value, chain)

    def observe_min(name: str, value: Fraction, chain: int) -> None:
        if name not in minima or value < minima[name][0]:
            minima[name] = (value, chain)

    for chain in sorted(recursive):
        low = evaluate(recursive[chain], special_records[chain], erfc_records[chain], PRECISIONS[0])
        high = evaluate(recursive[chain], special_records[chain], erfc_records[chain], PRECISIONS[1])
        for name in (
            "required_q_direct",
            "required_q_inverse",
            "required_compensation",
            "full_defect",
            "oriented_compensation",
            "required_q_cross_gap",
            "orientation_identity_gap",
        ):
            require(low[name].overlaps(high[name]), f"chain {chain} {name} precision nonoverlap")
        require(
            high["required_q_direct"].overlaps(high["required_q_inverse"]),
            f"chain {chain} q-need construction nonoverlap",
        )
        require(
            high["required_q_cross_gap"].contains(0),
            f"chain {chain} q-need cross-gap excludes zero",
        )
        require(
            high["orientation_identity_gap"].contains(0),
            f"chain {chain} oriented compensation identity excludes zero",
        )
        prior = complex_from_record(prior_records[chain])
        require(high["full_defect"].overlaps(prior), f"chain {chain} prior full-defect nonoverlap")
        prior_overlap_count += 1

        compensation_abs = abs(high["required_compensation"])
        compensation_lower = cells.bound_fraction(compensation_abs.lower())
        compensation_upper = cells.bound_fraction(compensation_abs.upper())
        full_q_lower = magnitude_lower(high["full_q"])
        require(full_q_lower > 0, f"chain {chain} full q contains zero")
        relative_upper = compensation_upper / full_q_lower
        cross_gap_upper = magnitude_upper(high["required_q_cross_gap"])
        identity_gap_upper = magnitude_upper(high["orientation_identity_gap"])
        magnitude_uppers.append(compensation_upper)
        observe_min("compensation", compensation_lower, chain)
        observe_max("compensation", compensation_upper, chain)
        observe_max("relative", relative_upper, chain)
        observe_max("cross_gap", cross_gap_upper, chain)
        observe_max("identity_gap", identity_gap_upper, chain)

        step = recursive[chain]["steps"][1]
        orientation = (
            ("conjugated" if step["conjugate"] else "direct")
            + ("_subtract_one" if step["subtract_one"] else "")
        )
        orientation_histogram[orientation] += 1
        real_sign = sign_class(high["required_compensation"].real)
        imag_sign = sign_class(high["required_compensation"].imag)
        sign_histogram[f"real_{real_sign}__imag_{imag_sign}"] += 1
        header = recursive[chain]["header"]
        rows.append(
            {
                "chain": chain,
                "sum_index": int(header["sum_index"]),
                "branch": int(header["branch"]),
                "child_length": int(header["kernel_length"]),
                "orientation": orientation,
                "source_q_hex": recursive[chain]["steps"][1]["qq_hex"],
                "full_corrected_q": complex_record(high["full_q"]),
                "required_q_direct": complex_record(high["required_q_direct"]),
                "required_q_inverse": complex_record(high["required_q_inverse"]),
                "required_native_q_compensation": complex_record(high["required_compensation"]),
                "required_native_q_compensation_magnitude_lower": fraction_decimal(compensation_lower),
                "required_native_q_compensation_magnitude_upper": fraction_decimal(compensation_upper),
                "required_native_q_compensation_real_sign": real_sign,
                "required_native_q_compensation_imag_sign": imag_sign,
                "required_native_q_compensation_relative_to_full_q_upper": fraction_decimal(relative_upper),
                "required_q_cross_gap_abs_upper": fraction_decimal(cross_gap_upper),
                "oriented_compensation": complex_record(high["oriented_compensation"]),
                "reconstructed_full_defect": complex_record(high["full_defect"]),
                "orientation_identity_gap_abs_upper": fraction_decimal(identity_gap_upper),
                "required_compensation_excludes_zero": not high["required_compensation"].contains(0),
                "precision_overlap": True,
                "prior_full_corrected_defect_overlap": True,
            }
        )

    ordered = sorted(magnitude_uppers)
    exclude_zero = sum(row["required_compensation_excludes_zero"] for row in rows)
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate",
        "status": "rigorous_374_call_native_q_required_compensation_inversion_enclosed",
        "scope": {
            "block": 20,
            "recursive_chain_count": 374,
            "precision_ladder_decimal_digits": list(PRECISIONS),
            "precision_overlap_count": 374,
            "prior_full_corrected_defect_overlap_count": prior_overlap_count,
            "parent_upper_index": 104,
            "child_upper_indices": [1, 2],
        },
        "identity": {
            "source_model": "P_adapted=O(P_raw), model(q)=O(M*C_adapted+q)",
            "required_q": "q_need=P_raw-M*C_adapted=O^{-1}(P_adapted)-M*C_adapted",
            "full_corrected_q": "q_full=q_source+delta_q_PSI_complex_ERF+delta_q_intrinsic_real_erfc",
            "required_compensation": "R_q=q_need-q_full",
            "defect_relation": "D_full=R_q for direct parent orientation and D_full=conj(R_q) for conjugated parent orientation; subtract-one cancels",
        },
        "normalization": {
            "phase": "Exact saved binary128 source tpm is used in independent parent and child sums.",
            "boundary": "This is the native source-normalized finite q target, not a true-2*pi or height-uniform replacement theorem.",
        },
        "aggregate": {
            "orientation_histogram": dict(orientation_histogram),
            "component_sign_histogram": dict(sign_histogram),
            "minimum_required_compensation_magnitude_lower": fraction_decimal(minima["compensation"][0]),
            "minimum_required_compensation_witness": minima["compensation"][1],
            "maximum_required_compensation_magnitude_upper": fraction_decimal(maxima["compensation"][0]),
            "maximum_required_compensation_witness": maxima["compensation"][1],
            "median_required_compensation_magnitude_upper": fraction_decimal(ordered[len(ordered) // 2]),
            "p95_required_compensation_magnitude_upper": fraction_decimal(ordered[(95 * (len(ordered) - 1)) // 100]),
            "maximum_relative_to_full_q_upper": fraction_decimal(maxima["relative"][0]),
            "maximum_relative_to_full_q_witness": maxima["relative"][1],
            "maximum_required_q_cross_gap_abs_upper": fraction_decimal(maxima["cross_gap"][0]),
            "maximum_required_q_cross_gap_witness": maxima["cross_gap"][1],
            "maximum_orientation_identity_gap_abs_upper": fraction_decimal(maxima["identity_gap"][0]),
            "maximum_orientation_identity_gap_witness": maxima["identity_gap"][1],
            "required_compensations_excluding_zero": exclude_zero,
        },
        "rows": rows,
        "source": {"path": relative(SOURCE), "sha256": file_hash(SOURCE), "locations": source_locations()},
        "sources": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "special_residual": {"path": relative(SPECIAL), "sha256": file_hash(SPECIAL)},
            "t5_real_erfc_residual": {"path": relative(REAL_ERFC), "sha256": file_hash(REAL_ERFC)},
            "full_corrected_defect": {"path": relative(FULL_CORRECTED), "sha256": file_hash(FULL_CORRECTED)},
            "selector_atlas": {"path": relative(ATLAS), "sha256": file_hash(ATLAS)},
            "point_ball_builder": {
                "path": relative(Path(point_balls.__file__).resolve()),
                "sha256": file_hash(Path(point_balls.__file__).resolve()),
            },
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "python_flint": {"version": point_balls.flint.__version__, "flint_version": point_balls.flint.__FLINT_VERSION__},
        },
        "next_handoff": {
            "analytic_target": "Express or bound R_q from the source-matched endpoint-saddle and Euler-Maclaurin remainder, including finite-ip and disabled correction terms.",
            "diagnostic_rule": "Candidate compensation laws must be derived first or tested on a predeclared selector-stratified holdout; the observed R_q roster is not itself a proof bound.",
            "uniform_target": "After a local identity or bound is established, derive explicit height and selector dependence before block accumulation."
        },
        "proof_boundary": (
            "Rigorous finite native-q inversion and required-compensation balls for 374 recursive calls of "
            "one low-height block only. The result specifies the correction any exact local recurrence would "
            "need, but it does not identify that correction with an omitted endpoint term, derive a remainder "
            "formula, prove a uniform bound, control outer accumulation, establish Lambda <= 0, PF-infinity, "
            "RH, or a prize-level conclusion."
        ),
    }
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(
        "built Hardy block-20 native-q compensation gate: "
        f"374 calls, {exclude_zero} nonzero targets, "
        f"range {artifact['aggregate']['minimum_required_compensation_magnitude_lower']} to "
        f"{artifact['aggregate']['maximum_required_compensation_magnitude_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
