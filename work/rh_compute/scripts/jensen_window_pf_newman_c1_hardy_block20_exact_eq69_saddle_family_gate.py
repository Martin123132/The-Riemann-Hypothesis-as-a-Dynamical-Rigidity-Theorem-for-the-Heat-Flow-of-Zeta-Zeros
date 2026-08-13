"""Enclose the exact equation-(69) saddle-family integrals on block 20."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import flint
from flint import acb, arb, ctx

import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells


SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
TELEMETRY = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10_block20_full/enabled/chain_telemetry.jsonl"
)
PAPER_PDF = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_Brereton_2607.15310.pdf"
NATIVE_Q = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_exact_eq69_saddle_family_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_exact_eq69_saddle_family_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_exact_eq69_saddle_family_gate.py"
PAPER_URL = "https://arxiv.org/abs/2607.15310"
PRECISIONS = (70, 110)
I = acb(0, 1)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def exact_arb(payload: str) -> arb:
    return point_balls.binary128_ball(payload)


def exact_acb(payload: list[str]) -> acb:
    return point_balls.binary128_complex(payload)


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def upper_abs(value: arb | acb) -> Fraction:
    return cells.bound_fraction(abs(value).upper())


def lower_abs(value: acb) -> Fraction:
    return cells.bound_fraction(abs(value).lower())


def ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def source_locations() -> dict[str, int]:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    tokens = {
        "recurrence_multiplier": "c1=exp(pp)/(EPI4*sqrt(x))",
        "q_call": "call q(m,k,ip,qq)",
        "recurrence_update": "csum=c1*csum+qq",
        "t5_start": "!      Next is t5.",
        "t5_finish": "t5=t5+(SUM0+SUM1)",
        "q_assembly": "qq=conjg(t1+t2+t4)+t3+t5",
    }
    result: dict[str, int] = {}
    for name, token in tokens.items():
        matches = [number for number, line in enumerate(lines, 1) if token in line]
        require(len(matches) == 1, f"equation-(69) source token drift for {name}: {len(matches)}")
        result[name] = matches[0]
    return result


def set_low_priority() -> str:
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


def saddle_integral_sum(data: dict[str, Any], dps: int, mathematical_pi: bool) -> dict[str, Any]:
    ctx.dps = dps
    ctx.threads = 1
    parent = data["levels"][1]
    child = data["levels"][2]
    require(int(parent["degree"]) == 3, "equation-(69) parent degree drift")
    require(int(parent["length"]) == 104, "equation-(69) parent length drift")
    require(int(child["length"]) in (1, 2), "equation-(69) child length drift")

    a1, a2, a3 = (exact_arb(value) for value in parent["coefficients_hex"])
    a1_fraction = cells.binary128_fraction(parent["coefficients_hex"][0])
    a2_fraction = cells.binary128_fraction(parent["coefficients_hex"][1])
    a3_fraction = cells.binary128_fraction(parent["coefficients_hex"][2])
    require(Fraction(-1, 2) < a1_fraction < Fraction(1, 2), "equation-(69) a1 cell drift")
    parent_length = int(parent["length"])
    require(a2_fraction > 0, "equation-(69) positive quadratic cell drift")
    require(
        2 * a2_fraction > 0 and 2 * a2_fraction + 6 * a3_fraction * parent_length > 0,
        "equation-(69) curvature-sign cell drift",
    )

    child_length = int(child["length"])
    jbot = ceil_fraction(a1_fraction)
    require(jbot in (0, 1), "equation-(69) starting index drift")
    require(bool(child["subtract_one"]) == (jbot == 1), "equation-(69) child-start adapter drift")

    source_two_pi = -exact_arb(data["header"]["tpm_hex"])
    phase = 2 * arb.pi() if mathematical_pi else source_two_pi
    tolerance = arb(f"1e-{dps - 20}")
    total = acb(0)
    integrals: list[acb] = []
    minimum_accuracy = 10**9
    for index in range(jbot, child_length + 1):
        n = arb(index)

        def integrand(y: acb, _analytic: bool) -> acb:
            g = n * y - a1 * y - a2 * y**2 - a3 * y**3
            return (I * phase * g).exp()

        value = acb.integral(
            integrand,
            arb(0),
            arb(parent_length),
            rel_tol=tolerance,
            abs_tol=tolerance,
            eval_limit=200000,
            depth_limit=30,
        )
        require(value.is_finite(), f"equation-(69) nonfinite integral at n={index}")
        minimum_accuracy = min(minimum_accuracy, int(value.rel_accuracy_bits()))
        total += value
        integrals.append(value)

    return {
        "jbot": jbot,
        "child_length": child_length,
        "indices": list(range(jbot, child_length + 1)),
        "integrals": integrals,
        "sum": total,
        "minimum_relative_accuracy_bits": minimum_accuracy,
        "source_two_pi": source_two_pi,
    }


def phase_transport_bound(data: dict[str, Any], source_two_pi: arb, dps: int) -> arb:
    ctx.dps = dps
    parent = data["levels"][1]
    child = data["levels"][2]
    a1 = exact_arb(parent["coefficients_hex"][0])
    a2 = exact_arb(parent["coefficients_hex"][1])
    a3 = exact_arb(parent["coefficients_hex"][2])
    a1_fraction = cells.binary128_fraction(parent["coefficients_hex"][0])
    jbot = ceil_fraction(a1_fraction)
    length = arb(int(parent["length"]))
    coefficient_bound = arb(0)
    for index in range(jbot, int(child["length"]) + 1):
        coefficient_bound += (
            abs(arb(index) - a1) * length**2 / 2
            + abs(a2) * length**3 / 3
            + abs(a3) * length**4 / 4
        )
    return abs(source_two_pi - 2 * arb.pi()) * coefficient_bound


def evaluate(data: dict[str, Any], dps: int) -> dict[str, Any]:
    ctx.dps = dps
    ctx.threads = 1
    step = data["steps"][1]
    parent = data["levels"][1]
    child = data["levels"][2]
    source = saddle_integral_sum(data, dps, mathematical_pi=False)
    mathematical = saddle_integral_sum(data, dps, mathematical_pi=True)
    require(source["indices"] == mathematical["indices"], "equation-(69) phase roster drift")

    tpm = exact_arb(data["header"]["tpm_hex"])
    parent_raw = point_balls.direct_sum(parent, tpm)
    child_raw = point_balls.direct_sum(child, tpm)
    child_adapted = point_balls.adapt(child, child_raw)
    multiplier = exact_acb(step["multiplier_hex"])
    source_q = exact_acb(step["qq_hex"])
    source_t5 = exact_acb(data["q_terms"][1]["t5_hex"])
    saddle_model = multiplier * child_adapted
    q_need = parent_raw - saddle_model
    source_required = q_need - source_q

    exact_w1_source_phase = source["sum"] - saddle_model
    exact_w1_paper_phase = mathematical["sum"] - saddle_model
    source_phase_replacement = exact_w1_source_phase - source_t5
    paper_phase_replacement = exact_w1_paper_phase - source_t5
    residual_after_source = source_required - source_phase_replacement
    residual_after_paper = source_required - paper_phase_replacement
    cancellation_residual = parent_raw - source["sum"] - (source_q - source_t5)
    cancellation_gap = residual_after_source - cancellation_residual
    transport = phase_transport_bound(data, source["source_two_pi"], dps)

    return {
        "indices": source["indices"],
        "source_integrals": source["integrals"],
        "paper_integrals": mathematical["integrals"],
        "source_integral_sum": source["sum"],
        "paper_integral_sum": mathematical["sum"],
        "phase_transport_bound": transport,
        "phase_transport_gap": mathematical["sum"] - source["sum"],
        "minimum_integral_accuracy_bits": min(
            source["minimum_relative_accuracy_bits"], mathematical["minimum_relative_accuracy_bits"]
        ),
        "parent_raw": parent_raw,
        "child_adapted": child_adapted,
        "saddle_model": saddle_model,
        "source_q": source_q,
        "source_t5": source_t5,
        "q_need": q_need,
        "source_required": source_required,
        "exact_w1_source_phase": exact_w1_source_phase,
        "exact_w1_paper_phase": exact_w1_paper_phase,
        "source_phase_replacement": source_phase_replacement,
        "paper_phase_replacement": paper_phase_replacement,
        "residual_after_source": residual_after_source,
        "residual_after_paper": residual_after_paper,
        "cancellation_residual": cancellation_residual,
        "cancellation_gap": cancellation_gap,
    }


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 exact equation-(69) saddle-family gate

Date: 2026-08-06

Status: rigorous finite source-normalized integration of the exact equation-(69) saddle family; not a proof of a height-uniform recurrence estimate or RH

## Exact finite object

For every recursive block-20 call, the parent has length `N=104`, the
transformed length is `L=1` or `2`, and the exact starting index is
`jbot=ceil(Phi1)` in `{{0,1}}`.  With

```text
g_n(y)=n y-Phi1 y-Phi2 y^2-Phi3 y^3,
```

this gate encloses every entire-function integral

```text
I_n=integral_0^104 exp(i*tpp*g_n(y)) dy,
jbot <= n <= L.                                        (1)
```

There are {aggregate['exact_integral_count']} integrals over the 374 calls.
Each is evaluated independently at 70 and 110 decimal digits by Arb's
certified complex integrator, with one active FLINT thread.  The two precision
levels overlap for every integral and every assembled quantity.

Let `M*C` be the saved recurrence multiplier times the independently summed
and exactly adapted child, and let `t5_src` be the emitted source component.
The source-relative exact saddle-family correction and replacement are

```text
W69_exact = sum_n I_n-M*C,
Delta69   = W69_exact-t5_src.                           (2)
```

No fitted compensation is used.  If `P` is the independently summed parent
and `q_src` the emitted full correction, the residual after replacing `t5`
has the two exactly equivalent constructions

```text
R_after = [P-M*C-q_src]-Delta69
        = P-sum_n I_n-(q_src-t5_src).                   (3)
```

The second line cancels the transformed saddle model completely.  The
checker requires the interval difference between the two constructions to
contain zero on all 374 calls.

## Result

```text
maximum |Delta69|                         <= {aggregate['maximum_exact_w1_replacement_magnitude_upper']}
maximum |Delta69|/|source native target|  <= {aggregate['maximum_replacement_relative_to_source_required_upper']}
minimum |R_after|                         >= {aggregate['minimum_residual_after_magnitude_lower']}
maximum |R_after|                         <= {aggregate['maximum_residual_after_magnitude_upper']}
improved / worsened / unresolved calls       {aggregate['improved_call_count']} / {aggregate['worsened_call_count']} / {aggregate['indeterminate_call_count']}
R_after balls excluding zero                 {aggregate['residual_after_excluding_zero_count']} / 374
```

Thus the complete finite equation-(69) saddle family, not merely the displayed
finite-`ip` approximation, can be replaced exactly on this roster.  It changes
the discrepancy materially on some calls but does not close any of the 374
recurrences.  The surviving local defect belongs to the other Euler--Maclaurin
pieces represented by `t1`--`t4` (paper W2--W5 and the endpoint half-sum), or
to their source-model identification, rather than to W1 alone.

## Pi provenance

The source-normalized integral uses the exact binary128 constant produced by
`tpp=8*atan(1)`.  The paper integral uses Arb's mathematical `2*pi`.  The
elementary inequality `|exp(iu)-exp(iv)|<=|u-v|`, integrated against an
explicit polynomial majorant for `|g_n|`, proves

```text
maximum source-tpp to mathematical-2*pi integral shift
  <= {aggregate['maximum_phase_transport_bound_upper']}.
```

Every observed interval shift lies inside that independent bound.  Pi here is
the Fourier normalization in `exp(2*pi*i*x)`; no circle or polygon is inserted
into the recurrence.

## Boundary

This is a rigorous finite exact-point replacement of the saddle-bearing
equation-(69) family on one saved block.  It is not a uniform asymptotic bound,
does not yet integrate the nonsaddle families (67b)--(67c), does not certify
W2--W5 or the outer Hardy representation, and proves no determinant/current
sign, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
"""


def main() -> int:
    for path in (SOURCE, TELEMETRY, PAPER_PDF, NATIVE_Q, CHECKER):
        require(path.is_file(), f"missing exact equation-(69) dependency: {path}")
    priority = set_low_priority()
    started = time.monotonic()
    recursive = native_q.load_recursive_chains()
    native_artifact = json.loads(NATIVE_Q.read_text(encoding="utf-8"))
    native_rows = {int(row["chain"]): row for row in native_artifact["rows"]}
    require(set(recursive) == set(native_rows), "equation-(69) native roster drift")

    rows: list[dict[str, Any]] = []
    histogram: Counter[str] = Counter()
    exact_integral_count = 0
    precision_overlap_count = 0
    residual_nonzero_count = 0
    improved = 0
    worsened = 0
    indeterminate = 0
    minimum_accuracy_bits = 10**9
    maxima: dict[str, tuple[Fraction, int]] = {}
    minima: dict[str, tuple[Fraction, int]] = {}

    def observe_max(name: str, value: Fraction, chain: int) -> None:
        if name not in maxima or value > maxima[name][0]:
            maxima[name] = (value, chain)

    def observe_min(name: str, value: Fraction, chain: int) -> None:
        if name not in minima or value < minima[name][0]:
            minima[name] = (value, chain)

    compared_names = (
        "source_integral_sum",
        "paper_integral_sum",
        "phase_transport_gap",
        "parent_raw",
        "child_adapted",
        "saddle_model",
        "source_q",
        "source_t5",
        "q_need",
        "source_required",
        "exact_w1_source_phase",
        "exact_w1_paper_phase",
        "source_phase_replacement",
        "paper_phase_replacement",
        "residual_after_source",
        "residual_after_paper",
        "cancellation_residual",
        "cancellation_gap",
    )

    for chain in sorted(recursive):
        low = evaluate(recursive[chain], PRECISIONS[0])
        high = evaluate(recursive[chain], PRECISIONS[1])
        require(low["indices"] == high["indices"], f"chain {chain} equation-(69) index drift")
        require(len(low["source_integrals"]) == len(high["source_integrals"]), f"chain {chain} integral count drift")
        for index, (low_value, high_value) in enumerate(zip(low["source_integrals"], high["source_integrals"])):
            require(low_value.overlaps(high_value), f"chain {chain} source integral {index} precision nonoverlap")
        for index, (low_value, high_value) in enumerate(zip(low["paper_integrals"], high["paper_integrals"])):
            require(low_value.overlaps(high_value), f"chain {chain} paper integral {index} precision nonoverlap")
        for name in compared_names:
            require(low[name].overlaps(high[name]), f"chain {chain} {name} precision nonoverlap")
        require(high["cancellation_gap"].contains(0), f"chain {chain} cancellation identity drift")
        require(upper_abs(high["phase_transport_gap"]) <= cells.bound_fraction(high["phase_transport_bound"].upper()), f"chain {chain} pi transport bound failed")
        precision_overlap_count += 1
        exact_integral_count += len(high["indices"])
        minimum_accuracy_bits = min(minimum_accuracy_bits, high["minimum_integral_accuracy_bits"])

        source_required_abs = abs(high["source_required"])
        replacement_abs = abs(high["source_phase_replacement"])
        residual_abs = abs(high["residual_after_source"])
        before_lower = cells.bound_fraction(source_required_abs.lower())
        before_upper = cells.bound_fraction(source_required_abs.upper())
        replacement_upper = cells.bound_fraction(replacement_abs.upper())
        after_lower = cells.bound_fraction(residual_abs.lower())
        after_upper = cells.bound_fraction(residual_abs.upper())
        relative_upper = replacement_upper / before_lower
        phase_bound_upper = cells.bound_fraction(high["phase_transport_bound"].upper())
        cancellation_gap_upper = upper_abs(high["cancellation_gap"])

        if after_upper < before_lower:
            classification = "improved"
            improved += 1
        elif after_lower > before_upper:
            classification = "worsened"
            worsened += 1
        else:
            classification = "indeterminate"
            indeterminate += 1
        residual_excludes_zero = not high["residual_after_source"].contains(0)
        residual_nonzero_count += residual_excludes_zero
        parent = recursive[chain]["levels"][1]
        child = recursive[chain]["levels"][2]
        histogram[f"L{child['length']}_jbot{high['indices'][0]}_conj{int(bool(child['conjugate']))}"] += 1

        observe_max("replacement", replacement_upper, chain)
        observe_max("relative", relative_upper, chain)
        observe_max("after", after_upper, chain)
        observe_max("phase", phase_bound_upper, chain)
        observe_max("cancellation", cancellation_gap_upper, chain)
        observe_min("after", after_lower, chain)

        rows.append(
            {
                "chain": chain,
                "sum_index": int(recursive[chain]["header"]["sum_index"]),
                "branch": int(recursive[chain]["header"]["branch"]),
                "parent_length": int(parent["length"]),
                "child_length": int(child["length"]),
                "child_conjugate": bool(child["conjugate"]),
                "child_subtract_one": bool(child["subtract_one"]),
                "indices": high["indices"],
                "source_phase_integral_sum": point_balls.acb_record(high["source_integral_sum"], 55),
                "paper_phase_integral_sum": point_balls.acb_record(high["paper_integral_sum"], 55),
                "source_to_paper_phase_transport_bound_upper": decimal(phase_bound_upper),
                "source_to_paper_observed_gap": point_balls.acb_record(high["phase_transport_gap"], 55),
                "exact_source_phase_w1": point_balls.acb_record(high["exact_w1_source_phase"], 55),
                "source_t5": point_balls.acb_record(high["source_t5"], 55),
                "exact_w1_replacement": point_balls.acb_record(high["source_phase_replacement"], 55),
                "exact_w1_replacement_magnitude_upper": decimal(replacement_upper),
                "replacement_relative_to_source_required_upper": decimal(relative_upper),
                "source_native_q_requirement": point_balls.acb_record(high["source_required"], 55),
                "residual_after_exact_w1": point_balls.acb_record(high["residual_after_source"], 55),
                "residual_after_magnitude_lower": decimal(after_lower),
                "residual_after_magnitude_upper": decimal(after_upper),
                "residual_after_excludes_zero": residual_excludes_zero,
                "magnitude_classification": classification,
                "cancellation_identity_gap_upper": decimal(cancellation_gap_upper),
                "minimum_integral_relative_accuracy_bits": high["minimum_integral_accuracy_bits"],
                "precision_overlap": True,
            }
        )

    elapsed = time.monotonic() - started
    require(exact_integral_count == 544, "equation-(69) aggregate integral count drift")
    require(precision_overlap_count == 374, "equation-(69) precision overlap drift")
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_exact_eq69_saddle_family_gate",
        "status": "rigorous_finite_exact_equation_69_saddle_family_enclosed",
        "scope": {
            "block": 20,
            "recursive_call_count": 374,
            "parent_length": 104,
            "child_lengths": [1, 2],
            "precision_ladder_decimal_digits": list(PRECISIONS),
            "precision_overlap_count": precision_overlap_count,
            "exact_integral_count": exact_integral_count,
            "minimum_integral_relative_accuracy_bits": minimum_accuracy_bits,
        },
        "identity": {
            "integral": "I_n=int_0^N exp(i*tpp*(n*y-sum_q Phi_q*y^q))dy",
            "source_relative_exact_w1": "W69_exact=sum_n I_n-M*C_adapted",
            "replacement": "Delta69=W69_exact-t5_source",
            "residual_first_construction": "R_after=(P-M*C_adapted-q_source)-Delta69",
            "residual_cancellation_construction": "R_after=P-sum_n I_n-(q_source-t5_source)",
            "starting_index": "jbot=ceil(Phi1), equal to 0 or 1 on the admitted selector cells",
            "pi_transport": "|I_math-I_source| <= |2*pi-tpp| sum_n int_0^N |g_n(y)|dy, using a coefficientwise polynomial majorant",
        },
        "aggregate": {
            "index_orientation_histogram": dict(histogram),
            "exact_integral_count": exact_integral_count,
            "precision_overlap_count": precision_overlap_count,
            "minimum_integral_relative_accuracy_bits": minimum_accuracy_bits,
            "maximum_exact_w1_replacement_magnitude_upper": decimal(maxima["replacement"][0]),
            "maximum_exact_w1_replacement_witness": maxima["replacement"][1],
            "maximum_replacement_relative_to_source_required_upper": decimal(maxima["relative"][0]),
            "maximum_replacement_relative_witness": maxima["relative"][1],
            "minimum_residual_after_magnitude_lower": decimal(minima["after"][0]),
            "minimum_residual_after_witness": minima["after"][1],
            "maximum_residual_after_magnitude_upper": decimal(maxima["after"][0]),
            "maximum_residual_after_witness": maxima["after"][1],
            "residual_after_excluding_zero_count": residual_nonzero_count,
            "improved_call_count": improved,
            "worsened_call_count": worsened,
            "indeterminate_call_count": indeterminate,
            "maximum_phase_transport_bound_upper": decimal(maxima["phase"][0]),
            "maximum_phase_transport_witness": maxima["phase"][1],
            "maximum_cancellation_identity_gap_upper": decimal(maxima["cancellation"][0]),
            "maximum_cancellation_identity_gap_witness": maxima["cancellation"][1],
        },
        "rows": rows,
        "paper": {
            "url": PAPER_URL,
            "pdf": {"path": relative(PAPER_PDF), "sha256": file_hash(PAPER_PDF)},
            "equations": {
                "exact_saddle_family": "(69)",
                "contour_decomposition": "(71)",
                "diagonal_saddle_approximation": "(72)--(74)",
                "vertical_leg_approximations": "(75)--(78)",
                "assembled_w1": "(80)--(81)",
            },
        },
        "source": {"path": relative(SOURCE), "sha256": file_hash(SOURCE), "locations": source_locations()},
        "sources": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "native_q_target": {"path": relative(NATIVE_Q), "sha256": file_hash(NATIVE_Q)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "python_flint": {"version": flint.__version__, "flint_version": flint.__FLINT_VERSION__},
        },
        "runtime": {
            "worker_count": 1,
            "flint_thread_count": 1,
            "priority": priority,
            "elapsed_seconds": elapsed,
        },
        "next_handoff": {
            "target": "Apply the same exact finite-integral strategy to the nonsaddle equation-(67b)--(67c) families represented by t1, t2, and t4, while keeping t3 exact.",
            "reason": "All residuals survive the exact t5/W1 replacement, so extending the W1 saddle expansion is no longer the leading finite-point target.",
        },
        "proof_boundary": (
            "Rigorous finite exact-point integration of the equation-(69) saddle family on 374 saved block-20 calls. "
            "This is not a height-uniform recurrence theorem, does not certify the nonsaddle equation-(67b)--(67c) "
            "families, W2-W5, outer Hardy representation, Lambda <= 0, PF-infinity, RH, or a prize-level conclusion."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print(
        "built Hardy block-20 exact equation-(69) saddle-family gate: "
        f"{exact_integral_count} integrals, {improved} improved, "
        f"{residual_nonzero_count} residuals exclude zero"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
