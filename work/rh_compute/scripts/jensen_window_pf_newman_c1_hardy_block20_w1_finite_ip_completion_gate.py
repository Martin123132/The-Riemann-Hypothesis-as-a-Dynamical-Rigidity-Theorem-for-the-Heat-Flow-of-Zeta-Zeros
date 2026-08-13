#!/usr/bin/env python3
"""Enclose every displayed W1 endpoint term omitted by the source ip clipping."""

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

import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_block20_t5_real_erfc_residual_gate as t5_erfc
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import acb, arb, ctx


SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
TELEMETRY = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10_block20_full/enabled/chain_telemetry.jsonl"
)
PROBE_RESULT = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/t5_component_probe/t1e10_block20/probe_result.json"
)
PROBE_OUTPUT = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/t5_component_probe/t1e10_block20/probe_output.txt"
)
ATLAS = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_discrete_selector_atlas_gate.json"
NATIVE_Q = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_w1_finite_ip_completion_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_w1_finite_ip_completion_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_w1_finite_ip_completion_gate.py"
PAPER_URL = "https://arxiv.org/abs/2607.15310"
PRECISIONS = (180, 260)
I = acb(0, 1)
SOURCE_SQRT_TWO_BINARY32_BITS = 0x3FB504F3


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def exact_arb(payload: str) -> arb:
    return cells.fraction_ball(cells.binary128_fraction(payload))


def exact_acb(payload: list[str]) -> acb:
    require(len(payload) == 2, "complex binary128 payload drift")
    return acb(exact_arb(payload[0]), exact_arb(payload[1]))


def fraction_decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def magnitude_upper(value: arb | acb) -> Fraction:
    return cells.bound_fraction(abs(value).upper())


def magnitude_lower(value: acb) -> Fraction:
    return cells.bound_fraction(abs(value).lower())


def ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def binary32_fraction(bits: int) -> Fraction:
    sign = -1 if bits >> 31 else 1
    exponent = (bits >> 23) & 0xFF
    fraction = bits & ((1 << 23) - 1)
    require(0 < exponent < 0xFF, "only finite normalized binary32 values are admitted")
    significand = (1 << 23) + fraction
    power = exponent - 127 - 23
    value = Fraction(sign * significand, 1)
    return value * (2**power) if power >= 0 else value / (2 ** (-power))


def load_recursive_chains() -> dict[int, dict[str, Any]]:
    chains: dict[int, dict[str, Any]] = {}
    for line in TELEMETRY.read_text(encoding="utf-8").splitlines():
        record = json.loads(line)
        chain = int(record.get("chain", 0))
        if record["type"] == "chain":
            chains[chain] = {"header": record, "levels": {}, "q_terms": {}}
        elif record["type"] == "level":
            chains[chain]["levels"][int(record["level"])] = record
        elif record["type"] == "q_terms":
            chains[chain]["q_terms"][int(record["nit"])] = record
    recursive = {chain: data for chain, data in chains.items() if int(data["header"]["mit"]) == 2}
    require(len(recursive) == 374, "finite-ip recursive roster drift")
    return recursive


def source_locations() -> dict[str, int]:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    tokens = {
        "global_ip": "ip=3",
        "ip1_seed": "ip1=ip",
        "ip1_clip_test": "if (ip1.gt.(0.5*L(nit))) then",
        "ip1_clip_value": "ip1=max1(0.2*L(nit),1.0)",
        "lower_loop": "do i=jbot,ip1",
        "upper_loop": "do i=L(nit),L(nit)-ip1+1,-1",
        "lower_default_real_sqrt_two": "sav1=sp*sx*c*sqrt(2.0)",
        "upper_default_real_sqrt_two": "sav1=sp*sx*xbeta*sqrt(2.0)",
        "t5_finish": "t5=t5+(SUM0+SUM1)",
    }
    result: dict[str, int] = {}
    for name, token in tokens.items():
        matches = [number for number, line in enumerate(lines, 1) if token in line]
        require(len(matches) == 1, f"finite-ip source token drift for {name}: {len(matches)}")
        result[name] = matches[0]
    return result


def endpoint_term(
    region: str,
    index: int,
    parent_length: int,
    child_length: int,
    a1: arb,
    a2: arb,
    a3: arb,
    xr: arb,
    frac_l: arb,
    endpoint: acb,
    p: arb,
    sqrt_p: arb,
    tpp: arb,
    tpm: arb,
    epi4: acb,
    sqrt_two: arb,
) -> tuple[acb, acb, acb]:
    con1 = arb(index) - a1
    require(not con1.contains(0), f"finite-ip zero lower denominator at index {index}")
    c = con1 / xr - 3 * a3 * con1**2 / xr**3
    sx = xr.sqrt()
    gc = arb(index) * c - a1 * c - a2 * c**2 - a3 * c**3
    phase = (I * tpp * gc).exp()
    cr2 = I * a2 / p

    if region == "lower":
        erfc_value = (sqrt_p * sx * c * sqrt_two).erfc()
        exponential = (tpm * con1 * c).exp()
        cr3 = 1 / con1 + cr2 * (1 + tpp * c * con1 * (1 + p * con1 * c)) / con1**3
        algebraic = I * cr2 / (tpp * con1**3)
        full = -phase * erfc_value / (2 * sx * epi4) - I * (exponential * cr3 - cr2 / con1**3) / tpp
    elif region == "upper":
        xbeta = arb(parent_length) - c
        distance = arb(child_length) + frac_l - arb(index)
        require(not distance.contains(0), f"finite-ip zero upper denominator at index {index}")
        erfc_value = (sqrt_p * sx * xbeta * sqrt_two).erfc()
        exponential = (tpm * xbeta * distance).exp()
        cr3 = 1 / distance + cr2 * (
            1 + tpp * xbeta * distance * (1 + p * distance * xbeta)
        ) / distance**3
        algebraic = I * endpoint * cr2 / (tpp * distance**3)
        full = -erfc_value * phase / (2 * sx * epi4) - I * endpoint * (
            exponential * cr3 - cr2 / distance**3
        ) / tpp
    else:
        raise RuntimeError(f"unknown finite-ip endpoint region: {region}")
    return full, algebraic, full - algebraic


def evaluate(
    data: dict[str, Any],
    atlas_row: dict[str, Any],
    probe_row: dict[str, Any],
    constants: dict[str, str],
    native_record: dict[str, Any],
    dps: int,
) -> dict[str, Any]:
    ctx.dps = dps
    ctx.threads = 1
    header = data["header"]
    parent = data["levels"][1]
    child = data["levels"][2]
    parent_length = int(parent["length"])
    child_length = int(child["length"])
    require(parent_length == 104 and child_length in (1, 2), "finite-ip chain length drift")
    require(int(header["ip"]) == 3, "finite-ip source cutoff drift")
    require(int(parent["degree"]) == 3, "finite-ip parent degree drift")

    a1, a2, a3 = (exact_arb(value) for value in parent["coefficients_hex"])
    xr = exact_arb(parent["xr_hex"])
    frac_l = exact_arb(parent["frac_length_hex"])
    endpoint = exact_acb(data["q_terms"][1]["endpoint_hex"])
    p = exact_arb(constants["pi"])
    sqrt_p = exact_arb(constants["sqrt_pi"])
    tpp = exact_arb(constants["two_pi"])
    tpm = exact_arb(constants["minus_two_pi"])
    epi4 = acb(exact_arb(constants["epi4_real"]), exact_arb(constants["epi4_imag"]))
    sqrt_two = cells.fraction_ball(binary32_fraction(SOURCE_SQRT_TWO_BINARY32_BITS))
    require(cells.binary128_fraction(parent["xr_hex"]) == 2 * cells.binary128_fraction(parent["coefficients_hex"][1]), "xr != 2*a2")
    require(cells.binary128_fraction(constants["two_pi"]) == 2 * cells.binary128_fraction(constants["pi"]), "tpp != 2*p")
    require(cells.binary128_fraction(constants["minus_two_pi"]) == -cells.binary128_fraction(constants["two_pi"]), "tpm != -tpp")

    jbot = ceil_fraction(cells.binary128_fraction(parent["coefficients_hex"][0]))
    ip1 = int(atlas_row["q_cell"]["q_branches"]["ip1"])
    require(jbot == int(atlas_row["q_cell"]["q_branches"]["jbot"]), "finite-ip jbot drift")
    require(ip1 == 1, "finite-ip clipped value drift")

    retained_lower = list(range(jbot, ip1 + 1))
    retained_upper = list(range(child_length, child_length - ip1, -1))
    full_indices = list(range(jbot, child_length + 1))
    missing_lower = [index for index in full_indices if index not in retained_lower]
    missing_upper = [index for index in full_indices if index not in retained_upper]

    retained_lower_value = acb(0)
    retained_upper_value = acb(0)
    for index in retained_lower:
        retained_lower_value += endpoint_term(
            "lower", index, parent_length, child_length, a1, a2, a3, xr, frac_l,
            endpoint, p, sqrt_p, tpp, tpm, epi4, sqrt_two,
        )[0]
    for index in retained_upper:
        retained_upper_value += endpoint_term(
            "upper", index, parent_length, child_length, a1, a2, a3, xr, frac_l,
            endpoint, p, sqrt_p, tpp, tpm, epi4, sqrt_two,
        )[0]

    missing_value = acb(0)
    algebraic_value = acb(0)
    decaying_value = acb(0)
    term_rows: list[dict[str, Any]] = []
    for region, indices in (("lower", missing_lower), ("upper", missing_upper)):
        for index in indices:
            full, algebraic, decaying = endpoint_term(
                region, index, parent_length, child_length, a1, a2, a3, xr, frac_l,
                endpoint, p, sqrt_p, tpp, tpm, epi4, sqrt_two,
            )
            missing_value += full
            algebraic_value += algebraic
            decaying_value += decaying
            term_rows.append(
                {
                    "region": region,
                    "index": index,
                    "displayed_w1_term": point_balls.acb_record(full, 55),
                    "psi2_equivalent_algebraic_part": point_balls.acb_record(algebraic, 55),
                    "erfc_exponential_part": point_balls.acb_record(decaying, 55),
                }
            )

    required = native_q.complex_from_record(native_record["required_native_q_compensation"])
    residual_after = required - missing_value
    source_lower = exact_acb(probe_row["lower_sum_hex"])
    source_upper = exact_acb(probe_row["upper_sum_hex"])
    return {
        "jbot": jbot,
        "ip1": ip1,
        "child_length": child_length,
        "retained_lower_indices": retained_lower,
        "retained_upper_indices": retained_upper,
        "missing_lower_indices": missing_lower,
        "missing_upper_indices": missing_upper,
        "retained_lower_value": retained_lower_value,
        "retained_upper_value": retained_upper_value,
        "retained_lower_source_gap": retained_lower_value - source_lower,
        "retained_upper_source_gap": retained_upper_value - source_upper,
        "missing_value": missing_value,
        "algebraic_value": algebraic_value,
        "decaying_value": decaying_value,
        "decomposition_gap": missing_value - algebraic_value - decaying_value,
        "required_compensation": required,
        "residual_after": residual_after,
        "term_rows": term_rows,
    }


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 displayed-W1 finite-ip completion gate

Date: 2026-08-06

Status: rigorous evaluation of every displayed equation-(81) endpoint term omitted by the source index clipping; not an exact-contour theorem

## Finite index set

The paper is exact through equation (69).  Its `W1` approximation in equation
(81) contains lower- and upper-endpoint correction terms over the finite
secondary index set.  The paper notes that all terms could be summed, but uses
a practical endpoint cutoff.  The accepted source sets `ip=3`; when the
transformed length is 1 or 2 it clips this to `ip1=1`.

For each of the 374 recursive block-20 calls this gate constructs the complete
finite set `jbot <= n <= L`, subtracts the indices actually visited by each
source loop, and evaluates every remaining displayed `W1` term with Arb at 180
and 260 decimal digits.  There is no extrapolation beyond the finite sum.  To
isolate the source's index clipping, the two erfc arguments retain the actual
default-real `sqrt(2.0)` value, binary32 `0x3FB504F3`; its replacement by the
mathematical square root is a separate correction gate.

```text
calls whose two source loops already cover every displayed index = {aggregate['fully_covered_call_count']}
calls with at least one omitted displayed endpoint term          = {aggregate['calls_with_missing_indices']}
omitted lower endpoint terms                                    = {aggregate['missing_lower_term_count']}
omitted upper endpoint terms                                    = {aggregate['missing_upper_term_count']}
total omitted displayed terms                                   = {aggregate['missing_term_count']}
```

Each omitted term is split exactly into its algebraic inverse-cube part and
its erfc/exponential part.  The finite inverse-cube sums are the terms that can
equivalently be written as differences of order-two polygamma values; no
unavailable Maple expression is assumed.

## Comparison with the required native correction

Let `T_ip` be the sum of all omitted displayed terms and let `R_q` be the exact
native correction required by the independently summed parent and child.  The
remaining target is

```text
R_after = R_q - T_ip.
```

The finite roster gives

```text
maximum |T_ip|                              <= {aggregate['maximum_completion_magnitude_upper']}
maximum inverse-cube part                   <= {aggregate['maximum_algebraic_part_magnitude_upper']}
maximum erfc/exponential part               <= {aggregate['maximum_decaying_part_magnitude_upper']}
minimum |R_after|                           >= {aggregate['minimum_residual_after_magnitude_lower']}
maximum |R_after|                           <= {aggregate['maximum_residual_after_magnitude_upper']}
calls with R_after excluding zero              {aggregate['residual_after_excluding_zero_count']} / 374
rigorously improved / worsened / unresolved    {aggregate['improved_call_count']} / {aggregate['worsened_call_count']} / {aggregate['indeterminate_call_count']}
```

The retained-term replay gap is at most
`{aggregate['maximum_retained_source_replay_gap_abs_upper']}`.  This checks the
formula transcription against the independently admitted source component
probe; it is not used as an analytic contour-error bound.

## Boundary

This gate decides only the finite index truncation inside the already
approximate displayed `W1` formula.  It does not turn equations (72), (75), or
(77) into equalities, bound omitted higher saddle phase, validate an unknown
Maple implementation, control `W2`--`W4`, establish a height-uniform recurrence,
or prove `Lambda<=0`, PF-infinity, RH, or a prize-level theorem.
"""


def main() -> int:
    for path in (SOURCE, TELEMETRY, PROBE_RESULT, PROBE_OUTPUT, ATLAS, NATIVE_Q, CHECKER):
        require(path.is_file(), f"missing finite-ip dependency: {path}")
    probe_artifact = json.loads(PROBE_RESULT.read_text(encoding="utf-8"))
    require(probe_artifact["validation"]["saved_t5_bit_match_count"] == 374, "t5 probe not admitted")
    constants, probe = t5_erfc.load_probe()
    recursive = load_recursive_chains()
    atlas_artifact = json.loads(ATLAS.read_text(encoding="utf-8"))
    atlas = {row["chain"]: row for row in atlas_artifact["rows"] if row["mit"] == 2}
    native_artifact = json.loads(NATIVE_Q.read_text(encoding="utf-8"))
    native = {row["chain"]: row for row in native_artifact["rows"]}
    require(set(recursive) == set(atlas) == set(native) == set(probe), "finite-ip roster mismatch")

    rows: list[dict[str, Any]] = []
    maxima: dict[str, tuple[Fraction, int]] = {}
    minima: dict[str, tuple[Fraction, int]] = {}
    completion_uppers: list[Fraction] = []
    child_jbot_histogram: Counter[str] = Counter()
    missing_lower_count = 0
    missing_upper_count = 0
    calls_with_missing = 0
    fully_covered = 0
    residual_nonzero = 0
    improved = 0
    worsened = 0
    indeterminate = 0

    def observe_max(name: str, value: Fraction, chain: int) -> None:
        if name not in maxima or value > maxima[name][0]:
            maxima[name] = (value, chain)

    def observe_min(name: str, value: Fraction, chain: int) -> None:
        if name not in minima or value < minima[name][0]:
            minima[name] = (value, chain)

    for chain in sorted(recursive):
        low = evaluate(recursive[chain], atlas[chain], probe[chain], constants, native[chain], PRECISIONS[0])
        high = evaluate(recursive[chain], atlas[chain], probe[chain], constants, native[chain], PRECISIONS[1])
        for name in (
            "retained_lower_value", "retained_upper_value", "missing_value", "algebraic_value",
            "decaying_value", "decomposition_gap", "required_compensation", "residual_after",
        ):
            require(low[name].overlaps(high[name]), f"chain {chain} {name} precision nonoverlap")
        require(high["decomposition_gap"].contains(0), f"chain {chain} finite-ip decomposition drift")
        require(low["term_rows"] and high["term_rows"] or not low["term_rows"] and not high["term_rows"], f"chain {chain} term roster precision drift")
        require(low["missing_lower_indices"] == high["missing_lower_indices"], f"chain {chain} lower roster drift")
        require(low["missing_upper_indices"] == high["missing_upper_indices"], f"chain {chain} upper roster drift")

        retained_gap = max(
            magnitude_upper(high["retained_lower_source_gap"]),
            magnitude_upper(high["retained_upper_source_gap"]),
        )
        require(retained_gap < Fraction(1, 10**20), f"chain {chain} retained endpoint transcription mismatch")
        completion_abs = abs(high["missing_value"])
        completion_lower = cells.bound_fraction(completion_abs.lower())
        completion_upper = cells.bound_fraction(completion_abs.upper())
        algebraic_upper = magnitude_upper(high["algebraic_value"])
        decaying_upper = magnitude_upper(high["decaying_value"])
        residual_before_abs = abs(high["required_compensation"])
        residual_after_abs = abs(high["residual_after"])
        before_lower = cells.bound_fraction(residual_before_abs.lower())
        before_upper = cells.bound_fraction(residual_before_abs.upper())
        after_lower = cells.bound_fraction(residual_after_abs.lower())
        after_upper = cells.bound_fraction(residual_after_abs.upper())
        relative_upper = completion_upper / before_lower

        if after_upper < before_lower:
            classification = "improved"
            improved += 1
        elif after_lower > before_upper:
            classification = "worsened"
            worsened += 1
        else:
            classification = "indeterminate"
            indeterminate += 1
        residual_excludes_zero = not high["residual_after"].contains(0)
        residual_nonzero += residual_excludes_zero

        lower_count = len(high["missing_lower_indices"])
        upper_count = len(high["missing_upper_indices"])
        missing_lower_count += lower_count
        missing_upper_count += upper_count
        calls_with_missing += bool(lower_count + upper_count)
        fully_covered += not bool(lower_count + upper_count)
        child_jbot_histogram[f"L{high['child_length']}_jbot{high['jbot']}"] += 1
        completion_uppers.append(completion_upper)
        observe_max("completion", completion_upper, chain)
        observe_max("algebraic", algebraic_upper, chain)
        observe_max("decaying", decaying_upper, chain)
        observe_max("relative", relative_upper, chain)
        observe_max("after", after_upper, chain)
        observe_max("replay", retained_gap, chain)
        observe_min("after", after_lower, chain)

        rows.append(
            {
                "chain": chain,
                "sum_index": int(recursive[chain]["header"]["sum_index"]),
                "branch": int(recursive[chain]["header"]["branch"]),
                "child_length": high["child_length"],
                "jbot": high["jbot"],
                "ip1": high["ip1"],
                "retained_lower_indices": high["retained_lower_indices"],
                "retained_upper_indices": high["retained_upper_indices"],
                "missing_lower_indices": high["missing_lower_indices"],
                "missing_upper_indices": high["missing_upper_indices"],
                "omitted_displayed_w1_completion": point_balls.acb_record(high["missing_value"], 55),
                "psi2_equivalent_algebraic_part": point_balls.acb_record(high["algebraic_value"], 55),
                "erfc_exponential_part": point_balls.acb_record(high["decaying_value"], 55),
                "completion_magnitude_lower": fraction_decimal(completion_lower),
                "completion_magnitude_upper": fraction_decimal(completion_upper),
                "completion_relative_to_required_q_upper": fraction_decimal(relative_upper),
                "required_native_q_compensation": point_balls.acb_record(high["required_compensation"], 55),
                "residual_after_completion": point_balls.acb_record(high["residual_after"], 55),
                "residual_after_magnitude_lower": fraction_decimal(after_lower),
                "residual_after_magnitude_upper": fraction_decimal(after_upper),
                "residual_after_excludes_zero": residual_excludes_zero,
                "magnitude_classification": classification,
                "retained_source_replay_gap_abs_upper": fraction_decimal(retained_gap),
                "precision_overlap": True,
                "terms": high["term_rows"],
            }
        )

    ordered = sorted(completion_uppers)
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_w1_finite_ip_completion_gate",
        "status": "rigorous_displayed_w1_finite_index_completion_enclosed",
        "scope": {
            "block": 20,
            "recursive_call_count": 374,
            "parent_length": 104,
            "child_lengths": [1, 2],
            "source_ip": 3,
            "effective_ip1": 1,
            "precision_ladder_decimal_digits": list(PRECISIONS),
            "precision_overlap_count": 374,
            "source_sqrt_two_binary32_bits": f"{SOURCE_SQRT_TWO_BINARY32_BITS:08X}",
            "source_sqrt_two_exact": {
                "numerator": binary32_fraction(SOURCE_SQRT_TWO_BINARY32_BITS).numerator,
                "denominator": binary32_fraction(SOURCE_SQRT_TWO_BINARY32_BITS).denominator,
                "decimal": fraction_decimal(binary32_fraction(SOURCE_SQRT_TWO_BINARY32_BITS)),
            },
        },
        "identity": {
            "complete_index_set": "jbot<=n<=L=floor(L+fracL)",
            "lower_source_set": "jbot<=n<=ip1",
            "upper_source_set": "L-ip1+1<=n<=L",
            "completion": "T_ip=sum of displayed equation-(81) lower and upper terms on the respective finite set differences",
            "residual_after": "R_after=R_q-T_ip",
            "algebraic_split": "Each displayed term is split exactly into its inverse-cube part and erfc/exponential remainder.",
            "literal_kind_boundary": "The source-equivalent erfc arguments retain default-real binary32 sqrt(2.0); exact mathematical sqrt(2) is audited separately.",
        },
        "aggregate": {
            "child_length_jbot_histogram": dict(child_jbot_histogram),
            "fully_covered_call_count": fully_covered,
            "calls_with_missing_indices": calls_with_missing,
            "missing_lower_term_count": missing_lower_count,
            "missing_upper_term_count": missing_upper_count,
            "missing_term_count": missing_lower_count + missing_upper_count,
            "maximum_completion_magnitude_upper": fraction_decimal(maxima["completion"][0]),
            "maximum_completion_witness": maxima["completion"][1],
            "median_completion_magnitude_upper": fraction_decimal(ordered[len(ordered) // 2]),
            "maximum_algebraic_part_magnitude_upper": fraction_decimal(maxima["algebraic"][0]),
            "maximum_decaying_part_magnitude_upper": fraction_decimal(maxima["decaying"][0]),
            "maximum_completion_relative_to_required_q_upper": fraction_decimal(maxima["relative"][0]),
            "minimum_residual_after_magnitude_lower": fraction_decimal(minima["after"][0]),
            "minimum_residual_after_witness": minima["after"][1],
            "maximum_residual_after_magnitude_upper": fraction_decimal(maxima["after"][0]),
            "maximum_residual_after_witness": maxima["after"][1],
            "residual_after_excluding_zero_count": residual_nonzero,
            "improved_call_count": improved,
            "worsened_call_count": worsened,
            "indeterminate_call_count": indeterminate,
            "maximum_retained_source_replay_gap_abs_upper": fraction_decimal(maxima["replay"][0]),
            "maximum_retained_source_replay_gap_witness": maxima["replay"][1],
        },
        "rows": rows,
        "source": {"path": relative(SOURCE), "sha256": file_hash(SOURCE), "locations": source_locations()},
        "paper": {
            "url": PAPER_URL,
            "exact_through": "equation (69)",
            "displayed_component": "W1, equation (81)",
            "paper_cutoff_statement": "All latter terms may be summed; the paper states they are typically negligible away from P=3 endpoints.",
        },
        "sources": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "t5_probe_result": {"path": relative(PROBE_RESULT), "sha256": file_hash(PROBE_RESULT)},
            "t5_probe_output": {"path": relative(PROBE_OUTPUT), "sha256": file_hash(PROBE_OUTPUT)},
            "selector_atlas": {"path": relative(ATLAS), "sha256": file_hash(ATLAS)},
            "native_q_target": {"path": relative(NATIVE_Q), "sha256": file_hash(NATIVE_Q)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "python_flint": {"version": point_balls.flint.__version__, "flint_version": point_balls.flint.__FLINT_VERSION__},
        },
        "next_handoff": {
            "if_residual_survives": "Derive or enclose the difference between the exact equation-(69) contour integrals and the displayed saddle/vertical-leg approximations, rather than extending ip further.",
            "maple_boundary": "The inverse-cube finite sums are psi(2)-equivalent, but this gate does not assert identity with unavailable Maple code.",
        },
        "proof_boundary": (
            "Rigorous finite evaluation of the displayed equation-(81) index completion at 374 exact source points. "
            "The displayed formula is already an approximation, so this does not bound the contour, saddle-phase, "
            "vertical-leg, W2-W4, hierarchy, outer Hardy, or height-uniform remainders and does not establish "
            "Lambda <= 0, PF-infinity, RH, or a prize-level theorem."
        ),
    }
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(
        "built Hardy block-20 displayed-W1 finite-ip completion gate: "
        f"{fully_covered} fully covered, {missing_lower_count + missing_upper_count} omitted terms, "
        f"{residual_nonzero} residuals exclude zero"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
