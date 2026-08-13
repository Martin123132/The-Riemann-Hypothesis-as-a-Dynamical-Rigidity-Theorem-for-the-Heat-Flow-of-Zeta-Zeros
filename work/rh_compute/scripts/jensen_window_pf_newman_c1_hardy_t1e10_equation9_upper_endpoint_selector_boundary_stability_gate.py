#!/usr/bin/env python3
"""Certify the Hardy source upper endpoint across the next odd-selector boundary."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import re
import struct
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

from flint import arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_strip_poisson_reassembly_gate as selector_gate


SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
FIXTURE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/t1e10/zeta14cubicmult_et005_output.txt"
SELECTOR_GATE = selector_gate.RESULT
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_upper_endpoint_selector_boundary_stability_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_upper_endpoint_selector_boundary_stability_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 100
C = 159577
C_NEXT = C + 2
SAVED_HEIGHT = 10_000_000_000
SAVED_START = C
SAVED_ENDPOINT = 5_122_421
BOUNDARY_START = C_NEXT
BOUNDARY_ENDPOINT = 5_122_423
ET_NUMERATOR = 5
ET_DENOMINATOR = 1000
GUARD_LOWER_OFFSET = -8000
GUARD_UPPER_OFFSET = 200000

EXPECTED_GSUMS = [111, 123, 138, 153, 171, 190, 212]
EXPECTED_RAW_LENGTH_FLOORS = [
    1, 2, 3, 5, 8, 11, 15, 46, 49, 54, 61, 71, 82, 94, 108, 123, 141, 161,
    184, 209, 239, 272, 310, 352, 401, 457, 520, 591, 673, 766, 871, 991,
    1127, 1283, 1459,
]
EXPECTED_ODD_LENGTHS = [
    1, 1, 3, 5, 7, 11, 15, 21, 31, 45, 61, 71, 81, 93, 107, 123, 141, 161,
    183, 209, 239, 271, 309, 351, 401, 457, 519, 591, 673, 765, 871, 991,
    1127, 1283, 1459,
]
EXPECTED_SAVED_ENDPOINTS = [
    159777, 160221, 160713, 161817, 163653, 166389, 170949, 177733, 187061,
    200629, 220133, 246421, 276949, 311717, 351573, 397365, 449941, 510149,
    578837, 656853, 745893, 847653, 962981, 1094421, 1243669, 1414117,
    1608309, 1828789, 2079797, 2365573, 2690357, 3060085, 3480693, 3958965,
    4503381, 5122421,
]

# Roots of (target-sqrt(8t/pi))*t^(1/6)=float32(3.2).  The gate verifies
# the signs on both ends of each bracket; these decimals are not trusted data.
LOWER_TRANSITION_ROOT_CENTER = "10000002368.45285511579661582259313224965527527695873791980491842836971"
UPPER_TRANSITION_ROOT_CENTER = "10000253032.91689675017597978943749343258911823964641104679983845264314"
ROOT_BRACKET_RADIUS = "1e-24"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_binary_float(value: float) -> arb:
    numerator, denominator = value.as_integer_ratio()
    return arb(numerator) / denominator


def float32(value: float) -> float:
    return struct.unpack(">f", struct.pack(">f", value))[0]


def source_default_real(value: float) -> arb:
    return exact_binary_float(float32(value))


def arb_pow(base: arb, exponent: arb | int) -> arb:
    return (base.log() * exponent).exp()


def interval_with_offsets(center: arb, lower: int, upper: int) -> arb:
    require(lower < upper and (lower + upper) % 2 == 0, "guard offsets must have an integral midpoint")
    midpoint = (lower + upper) // 2
    radius = (upper - lower) // 2
    return center + midpoint + arb(0, radius)


def source_contract() -> dict[str, Any]:
    text = SOURCE.read_text(encoding="utf-8")
    required = (
        "call start(a,t6,yphase,transit,M2)",
        "N1=M2+ib",
        "if (g.lt.3.2) then",
        "M2=M2+2",
        "do iblock=rh_block_start,totblock",
        "aenums=floor((X**iblock)*ib/2.0,dp1)",
        "MT=floor(((((et/p)**rmc5)*a)**rmc2)*((e-1.0)**rmc3)/(2.0**rmc4),dp1)",
        "MT=floor((et**rmc1)*t6*(1/b-b)/(sp*(2.0**rmc6)),dp1)",
        "RN1=RN1startofblock + (2*(real(MT,dp)+1.0)*real(aenums,dp))",
        "raecutoff=RN1",
    )
    missing = [snippet for snippet in required if snippet not in text]
    require(not missing, f"source recurrence contract drift: {missing}")
    return {
        "required_snippet_count": len(required),
        "all_required_snippets_present": True,
        "pinned_upstream_commit": "2e16dac3206b707052c3ac4cacdf3d1a2325e636",
        "source_default_real_3_2": source_default_real(3.2).str(PRECISION, more=True),
        "source_default_real_3_2_hex": float32(3.2).hex(),
    }


def parse_fixture_endpoints() -> list[int]:
    text = FIXTURE.read_text(encoding="utf-8")
    values = re.findall(r"alpha value at end of this block=\s+([0-9]+)(?:\.0)?", text)
    return [int(value) for value in values]


def source_totblock_raw(height: arb) -> arb:
    et = arb(ET_NUMERATOR) / ET_DENOMINATOR
    # The unsuffixed LOG argument and decimal constants are default real in
    # the admitted gfortran build.  Model those values exactly as binary32.
    log_005_default = source_default_real(math.log(float32(0.005)))
    ratio = et.log() / log_005_default
    return arb_pow(ratio, arb(1) / 4) * (source_default_real(3.545) * height.log()) - source_default_real(61.8)


def recurrence_interval(height: arb, start: int, fifth: arb) -> dict[str, Any]:
    p = arb.pi()
    et = arb(ET_NUMERATOR) / ET_DENOMINATOR
    a = (8 * height / p).sqrt()
    y = arb(111) / 100
    ib_raw = arb_pow(height, fifth)
    require(ib_raw.lower() >= 100 and ib_raw.upper() < 101, "block-zero length floor changed")
    ib = 200
    total_raw = source_totblock_raw(height)
    require(total_raw.upper() < 25, "automatic block count left its max-floor plateau")
    total_blocks = 35
    ht = et.log() / height.log()
    x = ((fifth - ht) / 4).exp()

    endpoint = start + ib
    odd_length = 0
    gsum_count = 0
    rows: list[dict[str, Any]] = []
    floor_margins: list[arb] = [ib_raw - 100, 101 - ib_raw]
    branch_margins: list[arb] = []

    for block in range(1, total_blocks + 1):
        inner = block <= 7
        if inner:
            branch_gap = y * a - endpoint
            require(branch_gap.lower() > 0, f"inner/outer branch changed at block {block}")
            branch_margins.append(branch_gap)
            gsum_raw = arb_pow(x, block) * ib / 2
            gsum_count = EXPECTED_GSUMS[block - 1]
            require(
                gsum_raw.lower() >= gsum_count and gsum_raw.upper() < gsum_count + 1,
                f"Gauss-sum count floor changed at block {block}",
            )
            floor_margins.extend((gsum_raw - gsum_count, gsum_count + 1 - gsum_raw))
            e = arb(endpoint) / a
            raw_length = arb_pow(arb_pow(et / p, arb(1) / 2) * a, arb(1) / 2) * arb_pow(e - 1, arb(5) / 8)
        else:
            branch_gap = arb(endpoint) - y * a
            require(branch_gap.lower() > 0, f"inner/outer branch changed at block {block}")
            branch_margins.append(branch_gap)
            b = a / (2 * endpoint)
            raw_length = (
                arb_pow(et, arb(1) / 4)
                * arb_pow(height, arb(1) / 4)
                * (1 / b - b)
                / (p.sqrt() * arb_pow(arb(2), arb(7) / 8))
            )

        raw_floor = EXPECTED_RAW_LENGTH_FLOORS[block - 1]
        require(
            raw_length.lower() >= raw_floor and raw_length.upper() < raw_floor + 1,
            f"Gauss-sum length floor changed at block {block}",
        )
        floor_margins.extend((raw_length - raw_floor, raw_floor + 1 - raw_length))

        candidate = raw_floor
        capped = False
        if not inner and candidate > 1.5 * odd_length:
            candidate = math.floor(1.5 * odd_length)
            capped = True
        if candidate % 2 == 0:
            candidate -= 1
        odd_length = candidate
        require(odd_length == EXPECTED_ODD_LENGTHS[block - 1], f"odd length changed at block {block}")

        endpoint += 2 * (odd_length + 1) * gsum_count
        rows.append(
            {
                "block": block,
                "branch": "inner" if inner else "outer",
                "gsum_count": gsum_count,
                "raw_length_floor": raw_floor,
                "odd_length": odd_length,
                "cap_applied": capped,
                "endpoint": endpoint,
                "raw_length_ball": raw_length.str(22, more=True),
            }
        )

    minimum_floor_margin = min(floor_margins, key=lambda value: value.lower())
    minimum_branch_margin = min(branch_margins, key=lambda value: value.lower())
    return {
        "total_blocks": total_blocks,
        "block_zero_length": ib,
        "start": start,
        "block_zero_endpoint": start + ib,
        "final_endpoint": endpoint,
        "rows": rows,
        "minimum_floor_margin_lower_bound": arb(minimum_floor_margin.lower()).str(22, more=True),
        "minimum_branch_margin_lower_bound": arb(minimum_branch_margin.lower()).str(22, more=True),
        "source_totblock_raw_ball": total_raw.str(22, more=True),
    }


def transition_residual(height: arb, target: int, threshold: arb) -> arb:
    return (target - (8 * height / arb.pi()).sqrt()) * arb_pow(height, arb(1) / 6) - threshold


def transition_bracket(center_text: str, target: int, threshold: arb) -> tuple[arb, arb, arb]:
    center = arb(center_text)
    radius = arb(ROOT_BRACKET_RADIUS)
    lower = center - radius
    upper = center + radius
    lower_residual = transition_residual(lower, target, threshold)
    upper_residual = transition_residual(upper, target, threshold)
    require(lower_residual.lower() > 0, f"lower root sign failed for target {target}")
    require(upper_residual.upper() < 0, f"upper root sign failed for target {target}")
    return center + arb(0, radius), lower_residual, upper_residual


def certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    p = arb.pi()
    threshold = source_default_real(3.2)
    selector_height = p * C**2 / 8
    guard_height = interval_with_offsets(selector_height, GUARD_LOWER_OFFSET, GUARD_UPPER_OFFSET)
    left_guard = interval_with_offsets(selector_height, GUARD_LOWER_OFFSET, 0)
    right_guard = interval_with_offsets(selector_height, 0, GUARD_UPPER_OFFSET)

    source_fifth = exact_binary_float(0.2)
    intended_fifth = arb(1) / 5
    saved_model = recurrence_interval(arb(SAVED_HEIGHT), SAVED_START, source_fifth)
    require(saved_model["final_endpoint"] == SAVED_ENDPOINT, "saved endpoint model mismatch")
    fixture_endpoints = parse_fixture_endpoints()
    require(fixture_endpoints == EXPECTED_SAVED_ENDPOINTS, "pinned fixture endpoint sequence changed")
    require([saved_model["block_zero_endpoint"]] + [row["endpoint"] for row in saved_model["rows"]] == fixture_endpoints, "source recurrence does not reproduce fixture")

    left_a = (8 * left_guard / p).sqrt()
    right_a = (8 * right_guard / p).sqrt()
    require(left_a.lower() > C - 1, "left guard crosses the preceding odd cell")
    require(right_a.upper() < C_NEXT + 2, "right guard crosses the following odd cell")
    left_transition = (C - left_a) * arb_pow(left_guard, arb(1) / 6)
    next_transition = (C_NEXT + 2 - right_a) * arb_pow(right_guard, arb(1) / 6)
    require(left_transition.upper() < threshold.lower(), "left guard exits the active source transition")
    require(next_transition.lower() > threshold.upper(), "right guard enters the next source transition")

    source_model = recurrence_interval(guard_height, BOUNDARY_START, source_fifth)
    intended_model = recurrence_interval(guard_height, BOUNDARY_START, intended_fifth)
    require(source_model["final_endpoint"] == BOUNDARY_ENDPOINT, "source-literal endpoint changed on guard")
    require(intended_model["final_endpoint"] == BOUNDARY_ENDPOINT, "intended-real endpoint changed on guard")
    require(
        [(row["gsum_count"], row["odd_length"], row["endpoint"]) for row in source_model["rows"]]
        == [(row["gsum_count"], row["odd_length"], row["endpoint"]) for row in intended_model["rows"]],
        "source and intended recurrence rosters differ",
    )

    lower_root, lower_left, lower_right = transition_bracket(LOWER_TRANSITION_ROOT_CENTER, C, threshold)
    upper_root, upper_left, upper_right = transition_bracket(UPPER_TRANSITION_ROOT_CENTER, C_NEXT, threshold)
    lower_distance = selector_height - lower_root
    upper_distance = upper_root - selector_height
    require(lower_distance.lower() > abs(GUARD_LOWER_OFFSET), "lower guard is not inside transition plateau")
    require(upper_distance.lower() > GUARD_UPPER_OFFSET, "upper guard is not inside transition plateau")

    old_terms_at_boundary = (BOUNDARY_ENDPOINT - C) // 2 + 1
    new_terms_at_boundary = (BOUNDARY_ENDPOINT - C_NEXT) // 2 + 1
    require(old_terms_at_boundary - new_terms_at_boundary == 1, "boundary fixed-B selector count mismatch")

    return {
        "saved_height": SAVED_HEIGHT,
        "saved_effective_start": SAVED_START,
        "saved_upper_endpoint": SAVED_ENDPOINT,
        "saved_fixture_endpoint_count_including_block_zero": len(fixture_endpoints),
        "saved_fixture_reproduced_exactly": True,
        "selector_C": C,
        "next_selector_C": C_NEXT,
        "selector_boundary_ball": selector_height.str(PRECISION, more=True),
        "boundary_effective_start_from_left": BOUNDARY_START,
        "boundary_effective_start_at_and_from_right": BOUNDARY_START,
        "boundary_upper_endpoint": BOUNDARY_ENDPOINT,
        "upper_endpoint_jump_at_selector_boundary": False,
        "upper_endpoint_change_precedes_selector_boundary": True,
        "upper_endpoint_change_size": BOUNDARY_ENDPOINT - SAVED_ENDPOINT,
        "endpoint_increment_above_effective_start": BOUNDARY_ENDPOINT - BOUNDARY_START,
        "source_default_real_transition_threshold": threshold.str(PRECISION, more=True),
        "source_transition_entry_ball": lower_root.str(PRECISION, more=True),
        "source_transition_entry_distance_below_selector_ball": lower_distance.str(PRECISION, more=True),
        "next_source_transition_entry_ball": upper_root.str(PRECISION, more=True),
        "next_source_transition_entry_distance_above_selector_ball": upper_distance.str(PRECISION, more=True),
        "transition_root_signs": {
            "lower_bracket_left_residual_ball": lower_left.str(22, more=True),
            "lower_bracket_right_residual_ball": lower_right.str(22, more=True),
            "upper_bracket_left_residual_ball": upper_left.str(22, more=True),
            "upper_bracket_right_residual_ball": upper_right.str(22, more=True),
        },
        "certified_height_guard_offsets": [GUARD_LOWER_OFFSET, GUARD_UPPER_OFFSET],
        "all_35_recurrence_rows_constant_on_guard": True,
        "source_and_intended_literal_rosters_agree_on_guard": True,
        "minimum_source_model_floor_margin_lower_bound": source_model["minimum_floor_margin_lower_bound"],
        "minimum_source_model_branch_margin_lower_bound": source_model["minimum_branch_margin_lower_bound"],
        "boundary_fixed_B_old_selector_term_count": old_terms_at_boundary,
        "boundary_fixed_B_new_selector_term_count": new_terms_at_boundary,
        "boundary_fixed_B_removed_term_count": 1,
        "boundary_rows": source_model["rows"],
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certified_upper_endpoint"]
    return f"""# Upper-endpoint stability at the odd-selector boundary

Date: 2026-08-12

Status: interval certificate for the pinned Hardy source endpoint recurrence;
not a proof of the complete source transition splice or RH

The saved output at `t=10^10` starts Block 0 at the effective odd integer
`{c['saved_effective_start']}` and ends Block 35 at

```text
B_saved={c['saved_upper_endpoint']}.                              (UE1)
```

The endpoint is not chosen by a separate cutoff.  After Block 0 the pinned
source iterates

```text
G_j=floor(X^j ib/2)                       (inner branch),
M_j=the source floor/cap/oddized Gauss-sum length,
B_j=B_(j-1)+2(M_j+1)G_j.                  (UE2)
```

The extracted recurrence reproduces all
`{c['saved_fixture_endpoint_count_including_block_zero']}` printed endpoints,
including (UE1), exactly.

Let

```text
t*=pi*{c['selector_C']}^2/8
  ={c['selector_boundary_ball']}.
```

The source `start` routine uses the unsuffixed default-real test `g<3.2`.
Under the admitted build that threshold is

```text
{c['source_default_real_transition_threshold']}.                  (UE3)
```

As the height approaches `t*` from below, this transition advances the
effective start from `{c['selector_C']}` to `{c['next_selector_C']}` before
the exact odd-selector boundary.  Its entry height is enclosed by

```text
{c['source_transition_entry_ball']},
```

which is

```text
{c['source_transition_entry_distance_below_selector_ball']}
```

below `t*`.  At and just above `t*`, the ordinary next-odd rule already gives
the same effective start `{c['next_selector_C']}`.  Consequently no start
index changes at `t*` itself.

Every automatic block count, inner/outer branch, Gauss-sum count, length
floor, cap, and oddization is constant on the Arb-certified interval

```text
t*{c['certified_height_guard_offsets'][0]:+d}
 <= t <=
t*{c['certified_height_guard_offsets'][1]:+d}.                    (UE4)
```

The conservative minimum raw floor slack is greater than
`{c['minimum_source_model_floor_margin_lower_bound']}`, and the conservative
minimum branch slack is greater than
`{c['minimum_source_model_branch_margin_lower_bound']}`.  Both the source-literal
and intended-real recurrences produce the same 35 rows.  Their common endpoint
is

```text
B_boundary={c['boundary_upper_endpoint']}.                         (UE5)
```

Thus `B` changed by `+2` during the earlier transition layer, but it does not
change simultaneously with `C -> C+2` at `t*`.  The next analogous transition
begins `{c['next_source_transition_entry_distance_above_selector_ball']}`
above `t*`.

For the exact fixed-`B` selector identity at the boundary, using
`B={c['boundary_upper_endpoint']}` gives `{c['boundary_fixed_B_old_selector_term_count']}`
formal old-selector lattice terms and `{c['boundary_fixed_B_new_selector_term_count']}`
new-selector terms, again differing by one.  Hence there is no additional
upper strip at the selector height.

Proof boundary: this closes only the question whether the pinned upper
endpoint jumps at the same odd-selector boundary.  It does not yet reassemble
the earlier coupled event consisting of the lower transition term, the direct
roster shift, and the added upper lattice point.  It also does not prove the
399-mode quantitative remainder, complete `T_upper`, `Lambda<=0`, RH, or a
prize-level conclusion.
"""


def main() -> None:
    started = time.perf_counter()
    priority = selector_gate.event_gate.window_gate.cell_gate.ode_gate.set_low_priority()
    for dependency in (SOURCE, FIXTURE, SELECTOR_GATE, CHECKER):
        require(dependency.is_file(), f"missing dependency: {dependency}")
    ctx.dps = PRECISION
    contract = source_contract()
    certified = certificate()
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_upper_endpoint_selector_boundary_stability_gate",
        "status": "source_upper_endpoint_recurrence_and_selector_boundary_stability_certified",
        "passed": True,
        "source_contract": contract,
        "certified_upper_endpoint": certified,
        "decision": {
            "saved_B_5122421_reproduced": True,
            "boundary_B_equals_5122423": True,
            "upper_endpoint_changes_at_selector_boundary": False,
            "upper_endpoint_transition_obstruction_at_selector_boundary_closed": True,
            "earlier_coupled_start_transition_reassembled": False,
            "complete_source_selector_splice_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            "selector_strip_gate": {"path": relative(SELECTOR_GATE), "sha256": file_hash(SELECTOR_GATE)},
            "pinned_source": {"path": relative(SOURCE), "sha256": file_hash(SOURCE)},
            "saved_fixture": {"path": relative(FIXTURE), "sha256": file_hash(FIXTURE)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
        "next_obligation": (
            "At the earlier source transition-entry height, reassemble the lower transition term, the +2 direct-roster shift, "
            "and the simultaneously added upper lattice point; then instantiate the quantitative 399-mode selector chart at B=5122423."
        ),
        "proof_boundary": (
            "Pinned source-rule endpoint recurrence and a finite selector-boundary height interval only. No coupled transition-term "
            "reassembly, quantitative 399-mode bound, complete source splice, T_upper, Lambda<=0, RH, or prize-level conclusion is proved."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("built upper-endpoint stability gate: B=5122423 is constant through the C-to-C+2 selector boundary")


if __name__ == "__main__":
    main()
