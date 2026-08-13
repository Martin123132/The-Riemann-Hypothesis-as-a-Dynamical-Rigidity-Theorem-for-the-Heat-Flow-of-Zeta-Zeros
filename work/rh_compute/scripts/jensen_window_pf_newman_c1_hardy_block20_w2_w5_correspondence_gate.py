#!/usr/bin/env python3
"""Audit published W2--W5 against the block-20 source decomposition."""

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
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import acb, arb, ctx


SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
PAPER = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_Brereton_2607.15310.pdf"
TELEMETRY = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10_block20_full/enabled/chain_telemetry.jsonl"
)
PROBE_RESULT = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/special_function_probe/t1e10_block20/probe_result.json"
)
PROBE_OUTPUT = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/special_function_probe/t1e10_block20/probe_output.txt"
)
SPECIAL = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_special_function_residual_gate.json"
EQ69 = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_exact_eq69_saddle_family_gate.json"
T2 = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_t2_lplus1_mismatch_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_w2_w5_correspondence_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_w2_w5_correspondence_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_w2_w5_correspondence_gate.py"
PRECISIONS = (90, 150)
I = acb(0, 1)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def upper_abs(value: arb | acb) -> Fraction:
    return cells.bound_fraction(abs(value).upper())


def lower_abs(value: acb) -> Fraction:
    return cells.bound_fraction(abs(value).lower())


def set_process_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover - platform fallback
        return f"unavailable:{type(exc).__name__}"


def source_locations() -> dict[str, int]:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    tokens = {
        "endpoint_half_sum": "t3=0.5*(1.0+endpoint)",
        "w5_internal_t1": "t1=(0.0,1.0)*(conjg(endpoint)*con2-con3)/tpp",
        "w2_internal_t2": "t2=(0.0,1.0)*conjg(endpoint)*cr3/tpp",
        "w34_positive_endpoint": "t4=-(0.0,1.0)*conjg(endpoint)/(tpp*(L(nit)+fracL(nit)))+cr3",
        "w34_negative_endpoint": "t4=(0.0,1.0)*(phicoeff(1,nit)/(con1)-1.0)/tpp+cr3",
        "assembly": "qq=conjg(t1+t2+t4)+t3+t5",
    }
    locations: dict[str, int] = {}
    for name, token in tokens.items():
        matches = [number for number, line in enumerate(lines, 1) if token in line]
        require(len(matches) == 1, f"W2--W5 source token drift for {name}: {len(matches)}")
        locations[name] = matches[0]
    return locations


def exact_fraction_ball(payload: str) -> arb:
    return cells.fraction_ball(cells.binary128_fraction(payload))


def paper_components(data: dict[str, Any], dps: int) -> dict[str, Any]:
    ctx.dps = dps
    ctx.threads = 1
    parent = data["levels"][1]
    child = data["levels"][2]
    q_terms = data["q_terms"][1]
    require(int(parent["degree"]) == 3 and int(parent["length"]) == 104, "W2--W5 parent cell drift")

    coefficient_fractions = [cells.binary128_fraction(value) for value in parent["coefficients_hex"]]
    a1_fraction, a2_fraction, a3_fraction = coefficient_fractions
    require(a2_fraction > 0, "W2--W5 nonpositive Phi2")
    a1, a2, a3 = (cells.fraction_ball(value) for value in coefficient_fractions)
    parent_length = int(parent["length"])
    xi_fraction = a1_fraction + 2 * a2_fraction * parent_length + 3 * a3_fraction * parent_length**2
    child_length = int(child["length"])
    require(xi_fraction.numerator // xi_fraction.denominator == child_length, "W2--W5 exact xi selector drift")
    frac_fraction = xi_fraction - child_length
    require(Fraction(0) < frac_fraction < Fraction(1), "W2--W5 exact xi boundary drift")
    source_frac_fraction = cells.binary128_fraction(parent["frac_length_hex"])
    require(Fraction(0) < source_frac_fraction < Fraction(1), "W2--W5 stored frac(xi) boundary drift")

    xi = cells.fraction_ball(xi_fraction)
    source_xi = cells.fraction_ball(child_length + source_frac_fraction)
    delta = cells.fraction_ball(1 - frac_fraction)
    xr = 2 * a2
    sqrt_xr = xr.sqrt()
    pi = arb.pi()
    sqrt_two = arb(2).sqrt()
    epi_minus = acb(1 / sqrt_two, -1 / sqrt_two)
    phase_at_endpoint = a1 * parent_length + a2 * parent_length**2 + a3 * parent_length**3
    source_endpoint_target = acb(0, -2 * pi * phase_at_endpoint).exp()
    paper_internal_endpoint = source_endpoint_target.conjugate()

    def fresnel_boundary(argument: arb) -> acb:
        z = epi_minus * (pi / xr).sqrt() * argument
        erfc_value = 1 - z.erf()
        oscillation = acb(0, -pi * argument**2 / xr).exp()
        return I * epi_minus * oscillation * erfc_value / (2 * sqrt_xr)

    w2_internal = I * paper_internal_endpoint * (
        epi_minus
        * acb(0, -pi * delta**2 / xr).exp()
        * (1 - (epi_minus * (pi / xr).sqrt() * delta).erf())
        / (2 * sqrt_xr)
        - 1 / (2 * pi * delta)
    )
    w2 = w2_internal.conjugate()
    if a1_fraction < 0:
        selector = "W3"
        argument = a1 + 1
        w34_internal = fresnel_boundary(argument) - I / (2 * pi * argument)
    else:
        require(a1_fraction > 0, "W2--W5 Phi1 selector on zero")
        selector = "W4"
        argument = a1
        w34_internal = fresnel_boundary(argument) - I * paper_internal_endpoint / (2 * pi * xi)
    w34 = w34_internal.conjugate()

    w5_internal = I / (2 * pi) * (
        paper_internal_endpoint * ((xi + 1).digamma() - delta.digamma())
        - ((a1 + 1).digamma() - (child_length + 1 - a1).digamma())
    )
    w5 = w5_internal.conjugate()
    half_sum = ((1 + paper_internal_endpoint) / 2).conjugate()
    aggregate = w2 + w34 + w5 + half_sum
    source_endpoint = point_balls.binary128_complex(q_terms["endpoint_hex"])
    return {
        "selector": selector,
        "xi": xi,
        "source_xi": source_xi,
        "xi_transport": xi - source_xi,
        "delta": delta,
        "endpoint": source_endpoint_target,
        "paper_internal_endpoint": paper_internal_endpoint,
        "source_endpoint": source_endpoint,
        "endpoint_transport": source_endpoint_target - source_endpoint,
        "w2": w2,
        "w34": w34,
        "w5": w5,
        "half_sum": half_sum,
        "aggregate": aggregate,
    }


def evaluate(
    data: dict[str, Any],
    eq69_row: dict[str, Any],
    t2_row: dict[str, Any],
    dps: int,
) -> dict[str, Any]:
    ctx.dps = dps
    ctx.threads = 1
    paper = paper_components(data, dps)
    parent = data["levels"][1]
    q_terms = data["q_terms"][1]
    child_length = int(data["levels"][2]["length"])
    source_two_pi = -point_balls.binary128_ball(data["header"]["tpm_hex"])
    source_q = point_balls.binary128_complex(q_terms["qq_hex"])
    source_t1_q = point_balls.binary128_complex(q_terms["t1_hex"]).conjugate()
    source_t2_q = point_balls.binary128_complex(q_terms["t2_hex"]).conjugate()
    source_t3_q = point_balls.binary128_complex(q_terms["t3_hex"])
    source_t4_q = point_balls.binary128_complex(q_terms["t4_hex"]).conjugate()
    source_t5 = point_balls.binary128_complex(q_terms["t5_hex"])
    source_endpoint = point_balls.binary128_complex(q_terms["endpoint_hex"])
    q_extra = I * source_endpoint / (source_two_pi * (child_length + 1))
    source_t2_corrected_q = source_t2_q - q_extra
    source_component_sum = source_t1_q + source_t2_corrected_q + source_t3_q + source_t4_q
    source_nonsaddle = source_q - source_t5 - q_extra
    assembly_gap = source_nonsaddle - source_component_sum

    component_deltas = {
        "w2": paper["w2"] - source_t2_corrected_q,
        "w34": paper["w34"] - source_t4_q,
        "w5": paper["w5"] - source_t1_q,
        "half_sum": paper["half_sum"] - source_t3_q,
    }
    component_delta_sum = sum(component_deltas.values(), acb(0))
    nonsaddle_transport = paper["aggregate"] - source_nonsaddle
    component_transport_gap = nonsaddle_transport - (component_delta_sum - assembly_gap)

    tpm = point_balls.binary128_ball(data["header"]["tpm_hex"])
    parent_source = point_balls.direct_sum(parent, tpm)
    parent_paper = point_balls.direct_sum(parent, -2 * arb.pi())
    source_integral_sum = native_q.complex_from_record(eq69_row["source_phase_integral_sum"])
    paper_integral_sum = native_q.complex_from_record(eq69_row["paper_phase_integral_sum"])
    parent_transport = parent_paper - parent_source
    integral_transport = paper_integral_sum - source_integral_sum
    prior_residual = native_q.complex_from_record(t2_row["residual_after_deletion"])
    residual_transport_construction = (
        prior_residual + parent_transport - integral_transport - nonsaddle_transport
    )
    residual_direct_construction = parent_paper - paper_integral_sum - paper["aggregate"]
    residual_identity_gap = residual_direct_construction - residual_transport_construction

    return {
        **paper,
        "source_t1_q": source_t1_q,
        "source_t2_corrected_q": source_t2_corrected_q,
        "source_t3_q": source_t3_q,
        "source_t4_q": source_t4_q,
        "source_nonsaddle": source_nonsaddle,
        "source_component_sum": source_component_sum,
        "assembly_gap": assembly_gap,
        "component_deltas": component_deltas,
        "component_delta_sum": component_delta_sum,
        "nonsaddle_transport": nonsaddle_transport,
        "component_transport_gap": component_transport_gap,
        "parent_transport": parent_transport,
        "integral_transport": integral_transport,
        "prior_residual": prior_residual,
        "residual_transport_construction": residual_transport_construction,
        "residual_direct_construction": residual_direct_construction,
        "residual_identity_gap": residual_identity_gap,
    }


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 W2--W5 source/paper correspondence gate

Date: 2026-08-07

Status: exact orientation ledger and rigorous finite published-component evaluation; not a proof of an exact nonsaddle-contour theorem or RH

## Term ledger

The fixed paper is Lewis--Brereton, equations (89), (93), (94), (96), and
(97).  The source evaluates the printed W terms internally with the conjugate
of its saved endpoint, then the final `qq` assembly conjugates `t1`, `t2`, and
`t4` back into the source's `exp(-2*pi*i*f)` recurrence orientation:

```text
paper endpoint half-sum        (1+E)/2       source t3
paper W2                       equation 89   source conj(t2), after deleting the source-only L+1 term
paper W3 if Phi1<0             equation 93   source conj(t4), negative branch
paper W4 if Phi1>0             equation 94   source conj(t4), positive branch
paper W5                       equation 97   source conj(t1)
```

The source emits `qq=conjg(t1+t2+t4)+t3+t5`; hence the conjugations above are
part of the implementation orientation, not an inferred sign adjustment.
The previously certified `t2` deletion removes only
`i*E/[tpp*(L+1)]` from final `qq`.

## Published components

For each exact binary128 cubic parent, this gate reconstructs

```text
xi       = Phi1 + 2 Phi2 N + 3 Phi3 N^2,
delta    = ceil(xi)-xi,
E_src    = exp(-2*pi*i*f(N)),
Q_paper  = source-oriented conjugate of [(1+conj(E_src))/2 + W2 + (W3 or W4) + W5].
```

Every `erfc` and `digamma` in W2--W5 is evaluated by Arb at 90 and 150
decimal digits.  The selector is derived from exact dyadic `Phi1`.  Exact
dyadic evaluation of `xi=f'(N)` reproduces the saved child-length selector;
its last-bit difference from the Fortran-rounded `fracL` is retained as an
explicit transport term.  The roster contains {aggregate['w3_call_count']}
W3 calls and {aggregate['w4_call_count']} W4 calls.

Let `P_pi` be the directly summed parent with mathematical `2*pi`, and let
`I69_pi` be the independently certified mathematical-pi equation-(69)
integral family.  The complete published finite residual is

```text
R_paper = P_pi - I69_pi - Q_paper.                     (1)
```

It is independently reconstructed from the prior source-normalized residual
after the `L+1` deletion by transporting the parent, integral family, and
nonsaddle aggregate to mathematical pi.  Both constructions overlap on all
374 calls.

## Finite result

```text
minimum |R_paper|                         >= {aggregate['minimum_paper_residual_magnitude_lower']}
maximum |R_paper|                         <= {aggregate['maximum_paper_residual_magnitude_upper']}
improved versus post-L+1 source residual     {aggregate['improved_call_count']} / 374
worsened                                     {aggregate['worsened_call_count']} / 374
indeterminate                                {aggregate['indeterminate_call_count']} / 374
R_paper balls excluding zero                 {aggregate['paper_residual_excluding_zero_count']} / 374
maximum |W2 source-to-paper change|        <= {aggregate['maximum_w2_component_transport_upper']}
maximum |W3/W4 source-to-paper change|     <= {aggregate['maximum_w34_component_transport_upper']}
maximum |W5 source-to-paper change|        <= {aggregate['maximum_w5_component_transport_upper']}
maximum |half-sum source-to-paper change|  <= {aggregate['maximum_half_sum_component_transport_upper']}
maximum |complete nonsaddle transport|     <= {aggregate['maximum_nonsaddle_transport_upper']}
```

Thus the exact published W2--W5 special-function formulas are now separated
from the source's finite PSI/ERF implementations and from the earlier
source-only denominator.  Any residual in (1) is not licensed as an error
bound: it is the finite discrepancy left by the contour approximations used
to obtain equations (84), (85), (87), (90), and (91), together with any
remaining normalization issue not already excluded by the two-construction
identity.

## Pi provenance

Pi enters only through the Fourier character `e(x)=exp(2*pi*i*x)` used in
the paper.  The source instead stores binary128 `p=4*atan(1)` and
`tpp=2*p`.  This gate does not insert pi from a circle, polygon, or unrelated
geometric assumption: it evaluates both normalizations and explicitly
transports the parent sum, endpoint phase, equation-(69) integrals, and all
W2--W5 factors.

## Boundary

This is a rigorous finite exact-input audit of the published approximation
on 374 block-20 calls.  It does not yet integrate the exact infinite
nonsaddle families (67b)--(67c), prove a height-uniform contour remainder,
control all recurrence levels or the outer Hardy representation, or prove a
determinant/current sign, `Lambda<=0`, PF-infinity, RH, or a prize-level
conclusion.
"""


def main() -> int:
    dependencies = (SOURCE, PAPER, TELEMETRY, PROBE_RESULT, PROBE_OUTPUT, SPECIAL, EQ69, T2, CHECKER)
    for path in dependencies:
        require(path.is_file(), f"missing W2--W5 dependency: {path}")
    probe = json.loads(PROBE_RESULT.read_text(encoding="utf-8"))
    require(
        probe["status"] == "source_derived_374_call_probe_with_exact_t1_t2_t4_replay"
        and probe["validation"]["saved_component_bit_match_count"] == 1122,
        "W2--W5 source component replay not admitted",
    )
    special = json.loads(SPECIAL.read_text(encoding="utf-8"))
    require(special["scope"]["saved_t1_t2_t4_bit_replay_count"] == 1122, "W2--W5 special gate drift")
    eq69 = json.loads(EQ69.read_text(encoding="utf-8"))
    t2 = json.loads(T2.read_text(encoding="utf-8"))
    require(eq69["status"] == "rigorous_finite_exact_equation_69_saddle_family_enclosed", "W2--W5 eq69 drift")
    require(t2["status"] == "exact_source_paper_t2_lplus1_mismatch_isolated_and_finitely_corrected", "W2--W5 t2 drift")
    eq69_rows = {int(row["chain"]): row for row in eq69["rows"]}
    t2_rows = {int(row["chain"]): row for row in t2["rows"]}
    recursive = native_q.load_recursive_chains()
    require(set(recursive) == set(eq69_rows) == set(t2_rows), "W2--W5 roster mismatch")

    priority = set_process_priority()
    rows: list[dict[str, Any]] = []
    maxima: dict[str, tuple[Fraction, int]] = {}
    minima: dict[str, tuple[Fraction, int]] = {}
    selectors: Counter[str] = Counter()
    improved = worsened = indeterminate = residual_nonzero = 0

    def observe_max(name: str, value: Fraction, chain: int) -> None:
        if name not in maxima or value > maxima[name][0]:
            maxima[name] = (value, chain)

    def observe_min(name: str, value: Fraction, chain: int) -> None:
        if name not in minima or value < minima[name][0]:
            minima[name] = (value, chain)

    compared = (
        "endpoint",
        "xi_transport",
        "w2",
        "w34",
        "w5",
        "half_sum",
        "aggregate",
        "nonsaddle_transport",
        "parent_transport",
        "integral_transport",
        "residual_transport_construction",
        "residual_direct_construction",
        "residual_identity_gap",
        "component_transport_gap",
    )
    for chain in sorted(recursive):
        low = evaluate(recursive[chain], eq69_rows[chain], t2_rows[chain], PRECISIONS[0])
        high = evaluate(recursive[chain], eq69_rows[chain], t2_rows[chain], PRECISIONS[1])
        require(low["selector"] == high["selector"], f"chain {chain} W2--W5 selector precision drift")
        for name in compared:
            require(low[name].overlaps(high[name]), f"chain {chain} W2--W5 {name} precision nonoverlap")
        require(high["component_transport_gap"].contains(0), f"chain {chain} W2--W5 component identity drift")
        require(high["residual_identity_gap"].contains(0), f"chain {chain} W2--W5 residual identity drift")
        require(upper_abs(high["assembly_gap"]) < Fraction(1, 10**30), f"chain {chain} W2--W5 source assembly drift")

        selectors[high["selector"]] += 1
        prior_abs = abs(high["prior_residual"])
        residual_abs = abs(high["residual_direct_construction"])
        prior_lower = cells.bound_fraction(prior_abs.lower())
        prior_upper = cells.bound_fraction(prior_abs.upper())
        residual_lower = cells.bound_fraction(residual_abs.lower())
        residual_upper = cells.bound_fraction(residual_abs.upper())
        if residual_upper < prior_lower:
            classification = "improved"
            improved += 1
        elif residual_lower > prior_upper:
            classification = "worsened"
            worsened += 1
        else:
            classification = "indeterminate"
            indeterminate += 1
        excludes_zero = not high["residual_direct_construction"].contains(0)
        residual_nonzero += excludes_zero

        observe_min("residual", residual_lower, chain)
        observe_max("residual", residual_upper, chain)
        observe_max("w2", upper_abs(high["component_deltas"]["w2"]), chain)
        observe_max("w34", upper_abs(high["component_deltas"]["w34"]), chain)
        observe_max("w5", upper_abs(high["component_deltas"]["w5"]), chain)
        observe_max("half", upper_abs(high["component_deltas"]["half_sum"]), chain)
        observe_max("nonsaddle", upper_abs(high["nonsaddle_transport"]), chain)
        observe_max("endpoint", upper_abs(high["endpoint_transport"]), chain)
        observe_max("xi", upper_abs(high["xi_transport"]), chain)
        observe_max("parent", upper_abs(high["parent_transport"]), chain)
        observe_max("integral", upper_abs(high["integral_transport"]), chain)
        observe_max("assembly", upper_abs(high["assembly_gap"]), chain)
        observe_max("component_gap", upper_abs(high["component_transport_gap"]), chain)
        observe_max("residual_gap", upper_abs(high["residual_identity_gap"]), chain)

        rows.append(
            {
                "chain": chain,
                "sum_index": int(recursive[chain]["header"]["sum_index"]),
                "branch": int(recursive[chain]["header"]["branch"]),
                "child_length": int(recursive[chain]["levels"][2]["length"]),
                "selector": high["selector"],
                "paper_w2": point_balls.acb_record(high["w2"], 50),
                "paper_w3_or_w4": point_balls.acb_record(high["w34"], 50),
                "paper_w5": point_balls.acb_record(high["w5"], 50),
                "paper_half_sum": point_balls.acb_record(high["half_sum"], 50),
                "paper_nonsaddle_aggregate": point_balls.acb_record(high["aggregate"], 50),
                "source_to_paper_component_changes": {
                    name: point_balls.acb_record(value, 45) for name, value in high["component_deltas"].items()
                },
                "complete_nonsaddle_transport": point_balls.acb_record(high["nonsaddle_transport"], 50),
                "source_component_assembly_gap_upper": decimal(upper_abs(high["assembly_gap"])),
                "component_transport_identity_gap_upper": decimal(upper_abs(high["component_transport_gap"])),
                "parent_source_to_pi_transport_upper": decimal(upper_abs(high["parent_transport"])),
                "eq69_source_to_pi_transport_upper": decimal(upper_abs(high["integral_transport"])),
                "endpoint_source_to_pi_transport_upper": decimal(upper_abs(high["endpoint_transport"])),
                "xi_source_rounding_transport_upper": decimal(upper_abs(high["xi_transport"])),
                "prior_post_lplus1_residual": point_balls.acb_record(high["prior_residual"], 50),
                "paper_residual": point_balls.acb_record(high["residual_direct_construction"], 50),
                "paper_residual_magnitude_lower": decimal(residual_lower),
                "paper_residual_magnitude_upper": decimal(residual_upper),
                "paper_residual_excludes_zero": excludes_zero,
                "magnitude_classification": classification,
                "residual_two_construction_gap_upper": decimal(upper_abs(high["residual_identity_gap"])),
                "precision_overlap": True,
            }
        )

    require(sum(selectors.values()) == 374 and set(selectors) == {"W3", "W4"}, "W2--W5 selector roster drift")
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_w2_w5_correspondence_gate",
        "status": "published_w2_w5_orientation_and_finite_component_residual_rigorously_enclosed",
        "scope": {
            "block": 20,
            "recursive_call_count": 374,
            "precision_ladder_decimal_digits": list(PRECISIONS),
            "precision_overlap_count": 374,
            "process_priority": priority,
        },
        "correspondence": {
            "endpoint_half_sum": {"paper": "equation (96)", "source": "t3"},
            "W2": {"paper": "equation (89)", "source": "conj(t2) after deleting i*E/[tpp*(L+1)]"},
            "W3": {"paper": "equation (93), Phi1<0", "source": "conj(t4), negative branch"},
            "W4": {"paper": "equation (94), Phi1>0", "source": "conj(t4), positive branch"},
            "W5": {"paper": "equation (97)", "source": "conj(t1)"},
            "source_assembly": "qq=conjg(t1+t2+t4)+t3+t5",
        },
        "identities": {
            "paper_nonsaddle": "Q_paper_source=conj((1+E_paper)/2+W2+(W3 or W4)+W5), E_paper=conj(E_source)",
            "paper_residual": "R_paper=P_pi-I69_pi-Q_paper",
            "transport_residual": "R_paper=R_drop+(P_pi-P_src)-(I69_pi-I69_src)-(Q_paper-Q_src_drop)",
            "component_transport": "Q_paper-Q_src_drop=sum(component paper-source changes)-source assembly gap",
        },
        "aggregate": {
            "w3_call_count": selectors["W3"],
            "w4_call_count": selectors["W4"],
            "minimum_paper_residual_magnitude_lower": decimal(minima["residual"][0]),
            "minimum_paper_residual_witness": minima["residual"][1],
            "maximum_paper_residual_magnitude_upper": decimal(maxima["residual"][0]),
            "maximum_paper_residual_witness": maxima["residual"][1],
            "improved_call_count": improved,
            "worsened_call_count": worsened,
            "indeterminate_call_count": indeterminate,
            "paper_residual_excluding_zero_count": residual_nonzero,
            "maximum_w2_component_transport_upper": decimal(maxima["w2"][0]),
            "maximum_w2_component_transport_witness": maxima["w2"][1],
            "maximum_w34_component_transport_upper": decimal(maxima["w34"][0]),
            "maximum_w34_component_transport_witness": maxima["w34"][1],
            "maximum_w5_component_transport_upper": decimal(maxima["w5"][0]),
            "maximum_w5_component_transport_witness": maxima["w5"][1],
            "maximum_half_sum_component_transport_upper": decimal(maxima["half"][0]),
            "maximum_half_sum_component_transport_witness": maxima["half"][1],
            "maximum_nonsaddle_transport_upper": decimal(maxima["nonsaddle"][0]),
            "maximum_nonsaddle_transport_witness": maxima["nonsaddle"][1],
            "maximum_endpoint_source_to_pi_transport_upper": decimal(maxima["endpoint"][0]),
            "maximum_xi_source_rounding_transport_upper": decimal(maxima["xi"][0]),
            "maximum_parent_source_to_pi_transport_upper": decimal(maxima["parent"][0]),
            "maximum_eq69_source_to_pi_transport_upper": decimal(maxima["integral"][0]),
            "maximum_source_component_assembly_gap_upper": decimal(maxima["assembly"][0]),
            "maximum_component_transport_identity_gap_upper": decimal(maxima["component_gap"][0]),
            "maximum_residual_two_construction_gap_upper": decimal(maxima["residual_gap"][0]),
        },
        "rows": rows,
        "paper": {
            "path": relative(PAPER),
            "sha256": file_hash(PAPER),
            "pages": [29, 30],
            "equations": [89, 93, 94, 96, 97],
        },
        "source": {"path": relative(SOURCE), "sha256": file_hash(SOURCE), "locations": source_locations()},
        "sources": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "probe_result": {"path": relative(PROBE_RESULT), "sha256": file_hash(PROBE_RESULT)},
            "probe_output": {"path": relative(PROBE_OUTPUT), "sha256": file_hash(PROBE_OUTPUT)},
            "special_function_gate": {"path": relative(SPECIAL), "sha256": file_hash(SPECIAL)},
            "equation_69_gate": {"path": relative(EQ69), "sha256": file_hash(EQ69)},
            "t2_lplus1_gate": {"path": relative(T2), "sha256": file_hash(T2)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "pi_provenance": {
            "paper": "pi is the Fourier normalization in e(x)=exp(2*pi*i*x)",
            "source": "binary128 p=4*atan(1), tpp=2*p, tpm=-tpp",
            "transported_objects": ["parent sum", "endpoint phase", "equation-(69) integral family", "W2--W5"],
        },
        "next_target": {
            "target": "Derive and enclose the complete exact equation-(67b)--(67c) nonsaddle aggregate before the endpoint approximations in equations (84)--(91).",
            "reason": "The published W2--W5 evaluation is now isolated from source transcription and special-function implementation error; its surviving finite residual is the contour-remainder target, not an admissible bound.",
        },
        "proof_boundary": (
            "This is a 374-call exact-input audit of the published W2--W5 approximation. It does not integrate the exact "
            "infinite nonsaddle families, prove a height-uniform contour remainder, control the outer Hardy representation, "
            "or prove Lambda<=0, PF-infinity, RH, or a prize-level conclusion."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print(
        "built Hardy block-20 W2--W5 correspondence gate: "
        f"{selectors['W3']} W3, {selectors['W4']} W4, {residual_nonzero}/374 residuals exclude zero"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
