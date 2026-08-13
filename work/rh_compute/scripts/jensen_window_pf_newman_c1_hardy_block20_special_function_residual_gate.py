#!/usr/bin/env python3
"""Enclose source PSI/ERF residuals and their correlated q correction."""

from __future__ import annotations

from collections import Counter
from decimal import Decimal
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

import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import acb, arb, ctx


SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
PROBE_RESULT = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/special_function_probe/t1e10_block20/probe_result.json"
)
PROBE_OUTPUT = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/special_function_probe/t1e10_block20/probe_output.txt"
)
TELEMETRY = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10_block20_full/enabled/chain_telemetry.jsonl"
)
ATLAS = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_discrete_selector_atlas_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_special_function_residual_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_special_function_residual_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_special_function_residual_gate.py"


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


def interval_record(value: arb) -> dict[str, Any]:
    lower, upper = cells.bounds(value)
    return {
        "lower": fraction_decimal(lower),
        "upper": fraction_decimal(upper),
        "lower_dyadic": [lower.numerator, -(lower.denominator.bit_length() - 1)],
        "upper_dyadic": [upper.numerator, -(upper.denominator.bit_length() - 1)],
    }


def complex_interval_record(value: acb) -> dict[str, dict[str, Any]]:
    return {"real": interval_record(value.real), "imag": interval_record(value.imag)}


def magnitude_upper(value: arb | acb) -> Fraction:
    return cells.bound_fraction(abs(value).upper())


def psi_branch(argument: Fraction) -> str:
    xx = abs(argument)
    tol = Fraction(1, 10**12)
    for name, center in (("special_1", 1), ("special_2", 2), ("special_3", 3)):
        if abs(xx - center) < tol:
            return name + ("_reflected" if argument < 0 else "")
    if abs(xx - Fraction(1, 2)) < tol:
        return "special_half" + ("_reflected" if argument < 0 else "")
    if xx < 4:
        if xx < Fraction(1, 2):
            base = "small_center_0"
        elif xx <= Fraction(3, 2):
            base = "small_center_1"
        elif xx <= Fraction(5, 2):
            base = "small_center_2"
        elif xx <= Fraction(7, 2):
            base = "small_center_3"
        else:
            base = "small_center_4"
    else:
        base = "large_asymptotic"
    return base + ("_reflected" if argument < 0 else "")


def erf_branch(argument: Fraction) -> str:
    value = abs(argument)
    if value < Fraction(1, 10**14):
        return "tiny_linear"
    intervals = (
        (Fraction(4, 5), Fraction(9, 8), "taylor_1"),
        (Fraction(9, 8), Fraction(11, 8), "taylor_1_25"),
        (Fraction(11, 8), Fraction(13, 8), "taylor_1_5"),
        (Fraction(13, 8), Fraction(15, 8), "taylor_1_75"),
        (Fraction(15, 8), Fraction(17, 8), "taylor_2"),
        (Fraction(17, 8), Fraction(19, 8), "taylor_2_25"),
        (Fraction(19, 8), Fraction(21, 8), "taylor_2_5"),
    )
    for lower, upper, name in intervals:
        if lower <= value < upper:
            return name
    return "large_asymptotic" if value >= Fraction(21, 8) else "small_power"


def source_locations() -> dict[str, int]:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    tokens = {
        "pi_definition": "p=4*ATAN(C)",
        "two_pi_definition": "tpp=2*p",
        "sqrt_pi_definition": "sp=sqrt(p)",
        "epi4_definition": "epi4=(1.0,1.0)/sqrt(2*C)",
        "psi_first_pair": "call psi(con2,ps2)",
        "psi_second_pair": "call psi(phicoeff(1,nit)+1.0,ps2)",
        "psi_endpoint_pair": "call psi(con1,ps1)",
        "erf_first": "call erf(z,6,cr1)",
        "q_combination": "qq=conjg(t1+t2+t4)+t3+t5",
    }
    locations: dict[str, int] = {}
    for name, token in tokens.items():
        matches = [number for number, line in enumerate(lines, 1) if token in line]
        require(matches, f"source residual token missing: {token}")
        locations[name] = matches[0]
    return locations


def load_telemetry() -> dict[int, dict[str, Any]]:
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
    return {chain: data for chain, data in chains.items() if int(data["header"]["mit"]) == 2}


def load_probe() -> tuple[dict[str, str], dict[int, dict[str, Any]]]:
    lines = [line.split() for line in PROBE_OUTPUT.read_text(encoding="ascii").splitlines() if line.strip()]
    require(lines and lines[0][0] == "#constants" and len(lines[0]) == 7, "probe header drift")
    constants = {
        "pi": lines[0][1],
        "sqrt_pi": lines[0][2],
        "two_pi": lines[0][3],
        "minus_two_pi": lines[0][4],
        "epi4_real": lines[0][5],
        "epi4_imag": lines[0][6],
    }
    rows: dict[int, dict[str, Any]] = {}
    for fields in lines[1:]:
        require(len(fields) == 29, "probe residual field count drift")
        chain = int(fields[0])
        rows[chain] = {
            "psi": [
                {"argument_hex": fields[1 + 2 * index], "source_hex": fields[2 + 2 * index]}
                for index in range(6)
            ],
            "erf": [
                {"argument_hex": fields[13], "source_hex": fields[14:16]},
                {"argument_hex": fields[16], "source_hex": fields[17:19]},
            ],
            "erf_weights": [fields[19:21], fields[21:23]],
            "replayed_terms": {"t1_hex": fields[23:25], "t2_hex": fields[25:27], "t4_hex": fields[27:29]},
        }
    require(len(rows) == 374, "probe residual row count drift")
    return constants, rows


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 rigorous special-function residual gate

Date: 2026-08-06

Status: rigorous finite-point PSI/ERF residuals and correlated q replacement; not a proof or uniform q theorem

## Admitted source observations

The source-derived one-CPU probe copied the accepted initialization and the
actual `PSI` and `ERF` routines.  At all 374 recursive block-20 calls it replayed
the saved binary128 `t1`, `t2`, and `t4` components exactly:

```text
374 recursive calls,
2244 real PSI evaluations,
 748 complex ERF evaluations,
1122 / 1122 saved component pairs reproduced bit-for-bit.
```

The rigorous references are Arb `digamma(x)` and
`erf((1-i)x/sqrt(2))`.  The latter is exactly the source target
`erf(exp(-i*pi/4)x)` because `exp(-i*pi/4)=(1-i)/sqrt(2)`.  The evaluator's
ordinary constant is explicitly `p=4*atan(1)`, followed by `tpp=2*p`,
`sp=sqrt(p)`, and `EPI4=(1+i)/sqrt(2)`; no unexplained pi is inserted.

## Correlated PSI cancellation

Number the six source calls in execution order by `P1,...,P6`, let `E` be the
saved endpoint and `T=tpp`.  Direct substitution into the source formulas gives

```text
conj(t1)_PSI + (t5)_PSI
  = (i/T) * ( E*(P6-P1) + (P3-P5) ).
```

The `P2=PSI(1-fracL)` and `P4=PSI(L+1-a1)` residuals cancel exactly.  The gate
therefore propagates the four surviving PSI residuals and the two conjugated
ERF residuals as one complex ball for each call.

## Finite physical result

```text
maximum PSI residual             <= {aggregate['maximum_psi_residual_abs_upper']}
maximum residual difference P6-P1 <= {aggregate['maximum_psi_pair_P6_P1_abs_upper']}
maximum residual difference P3-P5 <= {aggregate['maximum_psi_pair_P3_P5_abs_upper']}
maximum complex ERF residual     <= {aggregate['maximum_erf_residual_abs_upper']}
maximum correlated PSI q shift  <= {aggregate['maximum_psi_q_delta_abs_upper']}
maximum correlated ERF q shift  <= {aggregate['maximum_erf_q_delta_abs_upper']}
maximum total special q shift   <= {aggregate['maximum_total_q_delta_abs_upper']}
median total special q shift    <= {aggregate['median_total_q_delta_abs_upper']}
calls with total shift > 1e-4  = {aggregate['q_delta_threshold_counts']['greater_than_1e-4']}
```

Every bound is evaluated at the exact saved binary128 physical argument and
contains the corresponding rigorous special-function value.  Because every
recursive chain has one q step and its final source transformations are affine
conjugation/subtraction, the same magnitude bounds the local parent-state change
from this special-function replacement.

This is not yet the full q error.  The real intrinsic `erfc` calls in `t5`, the
Euler-Maclaurin/saddle truncation, binary128 roundoff outside the observed
special-function outputs, recurrence and block accumulation, the Legendre tail,
and the outer Hardy representation remain separate obligations.  The result is
finite block-20 evidence, not a height-uniform theorem or RH.  No prize-level conclusion follows.
"""


def main() -> int:
    for path in (SOURCE, PROBE_RESULT, PROBE_OUTPUT, TELEMETRY, ATLAS, CHECKER):
        require(path.is_file(), f"missing special residual dependency: {path}")
    probe_artifact = json.loads(PROBE_RESULT.read_text(encoding="utf-8"))
    require(probe_artifact["validation"]["saved_component_bit_match_count"] == 1122, "probe not admitted")
    atlas = json.loads(ATLAS.read_text(encoding="utf-8"))
    expected_chains = {row["chain"] for row in atlas["rows"] if row["mit"] == 2}
    telemetry = load_telemetry()
    constants, probe = load_probe()
    require(set(telemetry) == set(probe) == expected_chains, "recursive residual roster mismatch")

    ctx.dps = 180
    ctx.threads = 1
    sqrt_two = arb(2).sqrt()
    imaginary_unit = acb(0, 1)
    tpp = exact_arb(constants["two_pi"])
    psi_histogram: Counter[str] = Counter()
    erf_histogram: Counter[str] = Counter()
    maxima: dict[str, tuple[Fraction, dict[str, Any]]] = {}
    q_delta_uppers: list[Fraction] = []
    rows: list[dict[str, Any]] = []

    def observe(name: str, value: Fraction, witness: dict[str, Any]) -> None:
        if name not in maxima or value > maxima[name][0]:
            maxima[name] = (value, witness)

    for chain in sorted(probe):
        data = telemetry[chain]
        header = data["header"]
        q_terms = data["q_terms"][1]
        probe_row = probe[chain]
        psi_rows: list[dict[str, Any]] = []
        psi_residuals: list[arb] = []
        for call_index, call in enumerate(probe_row["psi"], 1):
            argument_fraction = cells.binary128_fraction(call["argument_hex"])
            argument = cells.fraction_ball(argument_fraction)
            source_value = exact_arb(call["source_hex"])
            true_value = argument.digamma()
            residual = true_value - source_value
            residual_upper = magnitude_upper(residual)
            branch = psi_branch(argument_fraction)
            psi_histogram[branch] += 1
            observe("psi", residual_upper, {"chain": chain, "call": call_index, "branch": branch})
            psi_residuals.append(residual)
            psi_rows.append(
                {
                    "call": call_index,
                    "argument_hex": call["argument_hex"],
                    "source_hex": call["source_hex"],
                    "source_branch": branch,
                    "residual": interval_record(residual),
                    "residual_abs_upper": fraction_decimal(residual_upper),
                }
            )

        erf_rows: list[dict[str, Any]] = []
        erf_residuals: list[acb] = []
        for call_index, call in enumerate(probe_row["erf"], 1):
            argument_fraction = cells.binary128_fraction(call["argument_hex"])
            argument = cells.fraction_ball(argument_fraction)
            source_value = exact_acb(call["source_hex"])
            true_value = acb(argument / sqrt_two, -argument / sqrt_two).erf()
            residual = true_value - source_value
            residual_upper = magnitude_upper(residual)
            branch = erf_branch(argument_fraction)
            erf_histogram[branch] += 1
            observe("erf", residual_upper, {"chain": chain, "call": call_index, "branch": branch})
            erf_residuals.append(residual)
            erf_rows.append(
                {
                    "call": call_index,
                    "argument_hex": call["argument_hex"],
                    "source_hex": call["source_hex"],
                    "source_branch": branch,
                    "residual": complex_interval_record(residual),
                    "residual_abs_upper": fraction_decimal(residual_upper),
                }
            )

        endpoint = exact_acb(q_terms["endpoint_hex"])
        psi_pair_P6_P1 = psi_residuals[5] - psi_residuals[0]
        psi_pair_P3_P5 = psi_residuals[2] - psi_residuals[4]
        psi_pair_P6_P1_upper = magnitude_upper(psi_pair_P6_P1)
        psi_pair_P3_P5_upper = magnitude_upper(psi_pair_P3_P5)
        observe("psi_pair_P6_P1", psi_pair_P6_P1_upper, {"chain": chain})
        observe("psi_pair_P3_P5", psi_pair_P3_P5_upper, {"chain": chain})
        psi_delta = imaginary_unit / tpp * (
            endpoint * psi_pair_P6_P1 + psi_pair_P3_P5
        )
        weights = [exact_acb(payload) for payload in probe_row["erf_weights"]]
        erf_delta = (
            weights[0] * erf_residuals[0].conjugate()
            + weights[1] * erf_residuals[1].conjugate()
        )
        total_delta = psi_delta + erf_delta
        source_q = exact_acb(q_terms["qq_hex"])
        replaced_q = source_q + total_delta
        psi_delta_upper = magnitude_upper(psi_delta)
        erf_delta_upper = magnitude_upper(erf_delta)
        total_delta_upper = magnitude_upper(total_delta)
        q_delta_uppers.append(total_delta_upper)
        observe("psi_q", psi_delta_upper, {"chain": chain})
        observe("erf_q", erf_delta_upper, {"chain": chain})
        observe("total_q", total_delta_upper, {"chain": chain})

        rows.append(
            {
                "chain": chain,
                "sum_index": int(header["sum_index"]),
                "branch": int(header["branch"]),
                "child_length": int(header["kernel_length"]),
                "psi": psi_rows,
                "erf": erf_rows,
                "correlated_replacement": {
                    "psi_residual_pair_P6_P1": interval_record(psi_pair_P6_P1),
                    "psi_residual_pair_P6_P1_abs_upper": fraction_decimal(psi_pair_P6_P1_upper),
                    "psi_residual_pair_P3_P5": interval_record(psi_pair_P3_P5),
                    "psi_residual_pair_P3_P5_abs_upper": fraction_decimal(psi_pair_P3_P5_upper),
                    "psi_delta": complex_interval_record(psi_delta),
                    "psi_delta_abs_upper": fraction_decimal(psi_delta_upper),
                    "erf_delta": complex_interval_record(erf_delta),
                    "erf_delta_abs_upper": fraction_decimal(erf_delta_upper),
                    "total_q_delta": complex_interval_record(total_delta),
                    "total_q_delta_abs_upper": fraction_decimal(total_delta_upper),
                    "source_q_hex": q_terms["qq_hex"],
                    "special_replaced_q": complex_interval_record(replaced_q),
                },
            }
        )

    require(sum(psi_histogram.values()) == 2244, "PSI residual count drift")
    require(sum(erf_histogram.values()) == 748, "ERF residual count drift")
    ordered_q_deltas = sorted(q_delta_uppers)
    thresholds = {
        "greater_than_1e-3": Fraction(1, 10**3),
        "greater_than_1e-4": Fraction(1, 10**4),
        "greater_than_1e-5": Fraction(1, 10**5),
        "greater_than_1e-6": Fraction(1, 10**6),
        "greater_than_1e-8": Fraction(1, 10**8),
    }
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_special_function_residual_gate",
        "status": "rigorous_finite_point_special_residuals_and_correlated_q_replacement_enclosed",
        "scope": {
            "block": 20,
            "recursive_chain_count": 374,
            "psi_evaluation_count": 2244,
            "erf_evaluation_count": 748,
            "saved_t1_t2_t4_bit_replay_count": 1122,
            "precision_decimal_digits": 180,
        },
        "pi_provenance": {
            "source": "p=4*atan(1), tpp=2*p, sp=sqrt(p), EPI4=(1+i)/sqrt(2)",
            "rigorous_erf_target": "erf(exp(-i*pi/4)*x)=erf((1-i)*x/sqrt(2))",
            "implementation": "The algebraic (1-i)/sqrt(2) form is evaluated directly by Arb; no free geometric pi is introduced.",
        },
        "psi_correlation_identity": {
            "call_order": [
                "P1=PSI(L+fracL+1)",
                "P2=PSI(1-fracL)",
                "P3=PSI(a1+1)",
                "P4=PSI(L+1-a1)",
                "P5=PSI(jbot-a1)",
                "P6=PSI(jbot-(L+fracL))",
            ],
            "identity": "conj(t1)_PSI+(t5)_PSI=(i/tpp)*(endpoint*(P6-P1)+(P3-P5))",
            "exactly_cancelled_calls": ["P2", "P4"],
            "surviving_calls": ["P1", "P3", "P5", "P6"],
        },
        "erf_correlation_identity": {
            "identity": "delta_q_ERF=B1*conj(delta_ERF1)+B2*conj(delta_ERF2)",
            "weights": "B1 and B2 are exact binary128 source weights emitted by the admitted probe.",
        },
        "aggregate": {
            "psi_source_branch_histogram": dict(sorted(psi_histogram.items())),
            "erf_source_branch_histogram": dict(sorted(erf_histogram.items())),
            "maximum_psi_residual_abs_upper": fraction_decimal(maxima["psi"][0]),
            "maximum_psi_residual_witness": maxima["psi"][1],
            "maximum_psi_pair_P6_P1_abs_upper": fraction_decimal(maxima["psi_pair_P6_P1"][0]),
            "maximum_psi_pair_P6_P1_witness": maxima["psi_pair_P6_P1"][1],
            "maximum_psi_pair_P3_P5_abs_upper": fraction_decimal(maxima["psi_pair_P3_P5"][0]),
            "maximum_psi_pair_P3_P5_witness": maxima["psi_pair_P3_P5"][1],
            "maximum_erf_residual_abs_upper": fraction_decimal(maxima["erf"][0]),
            "maximum_erf_residual_witness": maxima["erf"][1],
            "maximum_psi_q_delta_abs_upper": fraction_decimal(maxima["psi_q"][0]),
            "maximum_psi_q_delta_witness": maxima["psi_q"][1],
            "maximum_erf_q_delta_abs_upper": fraction_decimal(maxima["erf_q"][0]),
            "maximum_erf_q_delta_witness": maxima["erf_q"][1],
            "maximum_total_q_delta_abs_upper": fraction_decimal(maxima["total_q"][0]),
            "maximum_total_q_delta_witness": maxima["total_q"][1],
            "median_total_q_delta_abs_upper": fraction_decimal(ordered_q_deltas[len(ordered_q_deltas) // 2]),
            "p95_total_q_delta_abs_upper": fraction_decimal(
                ordered_q_deltas[(95 * (len(ordered_q_deltas) - 1)) // 100]
            ),
            "q_delta_threshold_counts": {
                name: sum(value > threshold for value in q_delta_uppers)
                for name, threshold in thresholds.items()
            },
            "local_parent_state_shift_abs_upper": fraction_decimal(maxima["total_q"][0]),
        },
        "rows": rows,
        "source": {"path": relative(SOURCE), "sha256": file_hash(SOURCE), "locations": source_locations()},
        "sources": {
            "probe_result": {"path": relative(PROBE_RESULT), "sha256": file_hash(PROBE_RESULT)},
            "probe_output": {"path": relative(PROBE_OUTPUT), "sha256": file_hash(PROBE_OUTPUT)},
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "selector_atlas": {"path": relative(ATLAS), "sha256": file_hash(ATLAS)},
            "local_cell_builder": {
                "path": relative(Path(cells.__file__).resolve()),
                "sha256": file_hash(Path(cells.__file__).resolve()),
            },
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "python_flint": {"version": cells.flint.__version__, "flint_version": cells.flint.__FLINT_VERSION__},
        },
        "next_handoff": {
            "first": "Audit the real intrinsic erfc calls and endpoint saddle sums inside t5 at the same 374 points.",
            "second": "Separate mathematical q truncation error from observed binary128 roundoff and then propagate the local balls through the finite block sum.",
            "direct_route": "Audit the 50 MIT=1 direct kernels independently; they contain no q correction.",
        },
        "proof_boundary": (
            "Rigorous exact-point source PSI/ERF residuals and their correlated special-function replacement "
            "inside q for 374 finite block-20 calls only. This does not bound t5 intrinsic-erfc/saddle error, "
            "the full q truncation, recurrence or block accumulation, the Legendre tail, the outer Hardy "
            "representation, a height-uniform theorem, RH, or a prize-level conclusion."
        ),
    }
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(
        "built Hardy block-20 special-function residual gate: "
        f"2244 PSI, 748 ERF, max correlated q shift {artifact['aggregate']['maximum_total_q_delta_abs_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
