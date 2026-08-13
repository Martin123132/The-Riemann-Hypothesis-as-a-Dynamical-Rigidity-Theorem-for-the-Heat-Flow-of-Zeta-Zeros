#!/usr/bin/env python3
"""Enclose all 50 direct MIT=1 block-20 kernels with true and source 2*pi."""

from __future__ import annotations

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

import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import acb, arb, ctx


SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
TELEMETRY = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10_block20_full/enabled/chain_telemetry.jsonl"
)
ATLAS = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_discrete_selector_atlas_gate.json"
INITIAL_ADAPTER = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_initial_level_adapter_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_direct_kernel_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_direct_kernel_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_direct_kernel_gate.py"
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


def magnitude_upper(value: arb | acb) -> Fraction:
    return cells.bound_fraction(abs(value).upper())


def complex_record(value: acb) -> dict[str, Any]:
    return point_balls.acb_record(value, 55)


def decimal_quantum(value: Decimal) -> Fraction:
    exponent = value.as_tuple().exponent
    return Fraction(10**exponent, 1) if exponent >= 0 else Fraction(1, 10 ** (-exponent))


def printed_decimal_ball(payload: str) -> arb:
    value = Decimal(payload)
    center = Fraction(value)
    radius = decimal_quantum(value)
    return cells.interval_ball(center, radius)


def printed_complex_ball(payload: list[str]) -> acb:
    require(len(payload) == 2, "printed complex payload drift")
    return acb(printed_decimal_ball(payload[0]), printed_decimal_ball(payload[1]))


def load_direct_chains() -> dict[int, dict[str, Any]]:
    chains: dict[int, dict[str, Any]] = {}
    for line in TELEMETRY.read_text(encoding="utf-8").splitlines():
        record = json.loads(line)
        chain = int(record.get("chain", 0))
        if record["type"] == "chain":
            chains[chain] = {"header": record, "levels": {}, "final": None}
        elif record["type"] == "level":
            chains[chain]["levels"][int(record["level"])] = record
        elif record["type"] == "chain_end":
            chains[chain]["final"] = record
    direct = {chain: data for chain, data in chains.items() if int(data["header"]["mit"]) == 1}
    require(len(direct) == 50, "direct-kernel roster drift")
    return direct


def polynomial(level: dict[str, Any], index: int) -> arb:
    coefficients = [point_balls.binary128_ball(value) for value in level["coefficients_hex"]]
    value = arb(0)
    for coefficient in reversed(coefficients):
        value = (value + coefficient) * index
    return value


def evaluate(data: dict[str, Any], dps: int) -> dict[str, acb | arb]:
    ctx.dps = dps
    ctx.threads = 1
    header = data["header"]
    level = data["levels"][1]
    require(set(data["levels"]) == {1}, "MIT1 level shape drift")
    require(int(level["length"]) == 104 and int(level["degree"]) == 3, "direct kernel shape drift")
    require(level["subtract_one"] is False, "unexpected MIT1 subtract-one route")

    source_tpm = point_balls.binary128_ball(header["tpm_hex"])
    true_tpm = -2 * arb.pi()
    source_kernel = point_balls.adapt(level, point_balls.direct_sum(level, source_tpm))
    true_kernel = point_balls.adapt(level, point_balls.direct_sum(level, true_tpm))
    logged_kernel = printed_complex_ball(data["final"]["final_state"])
    source_roundoff_gap = source_kernel - logged_kernel
    phase_normalization_gap = true_kernel - source_kernel
    total_true_gap = true_kernel - logged_kernel

    phase_lipschitz = arb(0)
    delta_tpm = true_tpm - source_tpm
    for index in range(int(level["length"]) + 1):
        phase_lipschitz += abs(delta_tpm * polynomial(level, index))

    return {
        "source_tpm": source_tpm,
        "true_tpm": true_tpm,
        "source_kernel": source_kernel,
        "true_kernel": true_kernel,
        "logged_kernel": logged_kernel,
        "source_roundoff_gap": source_roundoff_gap,
        "phase_normalization_gap": phase_normalization_gap,
        "total_true_gap": total_true_gap,
        "phase_lipschitz": phase_lipschitz,
    }


def source_locations() -> dict[str, int]:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    tokens = {
        "pi_definition": "p=4*ATAN(C)",
        "tpm_definition": "tpm=-tpp",
        "kernel_loop": "do i=0,L(MIT-1)",
        "kernel_phase": "fn=(0.0,1.0)*tpm*i*(phicoeff(1,MIT)+i*(phicoeff(2,MIT)+i*phicoeff(3,MIT)))",
        "direct_accumulation": "csum=csum+exp(fn)",
        "orientation": "csum=conjg(csum)",
    }
    locations: dict[str, int] = {}
    for name, token in tokens.items():
        matches = [number for number, line in enumerate(lines, 1) if token in line]
        require(matches, f"direct-kernel source token missing: {token}")
        locations[name] = matches[-1] if name == "orientation" else matches[0]
    return locations


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 direct MIT=1 kernel gate

Date: 2026-08-06

Status: rigorous finite exact-point audit of all 50 direct kernels; not a proof or outer Hardy theorem

## Direct route

For `MIT=1` the source performs no recurrence and no `q` correction.  It sums

```text
S_hat = sum_(n=0)^104 exp(i*tpm_hat*(a1*n+a2*n^2+a3*n^3)),
```

then applies only the saved conjugation orientation.  All 50 scattered direct
calls are evaluated from their exact binary128 coefficients at 180 and 260
decimal digits.  Every low/high precision pair overlaps.  The saved direct
state was printed with 40 digits after the decimal rather than a hex payload;
each component is therefore enclosed by one full last-decimal unit before
comparison.

The source constant has explicit provenance: `p=4*atan(1)`, `tpp=2*p`, and
`tpm_hat=-tpp`.  The gate separately evaluates the mathematical normalization
`tpm=-2*pi` using Arb's standard pi ball.  No arbitrary circle or polygon is
introduced.

## Finite bounds

```text
maximum exact-source expression versus logged state <= {aggregate['maximum_source_roundoff_gap_abs_upper']}
maximum true-2*pi versus source-tpm kernel shift     <= {aggregate['maximum_phase_normalization_gap_abs_upper']}
maximum true-2*pi versus logged state gap            <= {aggregate['maximum_total_true_gap_abs_upper']}
maximum analytic phase-Lipschitz budget              <= {aggregate['maximum_phase_lipschitz_upper']}
minimum rigorous true-kernel magnitude               >= {aggregate['minimum_true_kernel_magnitude_lower']}
```

For real phases, `|exp(iu)-exp(iv)|<=|u-v|`; summing this inequality gives the
recorded phase-Lipschitz budget and independently dominates every observed
normalization shift.  The direct source computation is therefore enclosed at
all 50 saved points without borrowing the recursive `q` model.

This does not yet transport coefficient errors from a continuous physical
height interval, bound the outer Hardy representation or the cross-call block
sum, or address the remaining recursive `t5` saddle error.  It proves no
height-uniform theorem, determinant/current sign, `Lambda<=0`, PF-infinity, RH,
or prize-level conclusion.
"""


def main() -> int:
    for path in (SOURCE, TELEMETRY, ATLAS, INITIAL_ADAPTER, CHECKER):
        require(path.is_file(), f"missing direct-kernel dependency: {path}")
    atlas = json.loads(ATLAS.read_text(encoding="utf-8"))
    expected = {row["chain"] for row in atlas["rows"] if row["mit"] == 1}
    direct = load_direct_chains()
    require(set(direct) == expected, "direct-kernel atlas roster mismatch")

    rows: list[dict[str, Any]] = []
    maxima: dict[str, tuple[Fraction, int]] = {}
    minimum_magnitude: tuple[Fraction, int] | None = None

    def observe(name: str, value: Fraction, chain: int) -> None:
        if name not in maxima or value > maxima[name][0]:
            maxima[name] = (value, chain)

    for chain in sorted(direct):
        low = evaluate(direct[chain], PRECISIONS[0])
        high = evaluate(direct[chain], PRECISIONS[1])
        for name in (
            "source_kernel",
            "true_kernel",
            "source_roundoff_gap",
            "phase_normalization_gap",
            "total_true_gap",
        ):
            require(low[name].overlaps(high[name]), f"chain {chain} {name} precision nonoverlap")
        require(low["phase_lipschitz"].overlaps(high["phase_lipschitz"]), f"chain {chain} Lipschitz nonoverlap")

        source_roundoff_upper = magnitude_upper(high["source_roundoff_gap"])
        phase_gap_upper = magnitude_upper(high["phase_normalization_gap"])
        total_gap_upper = magnitude_upper(high["total_true_gap"])
        lipschitz_upper = magnitude_upper(high["phase_lipschitz"])
        true_magnitude_lower = cells.bound_fraction(abs(high["true_kernel"]).lower())
        require(phase_gap_upper <= lipschitz_upper, f"chain {chain} phase Lipschitz failure")
        observe("source_roundoff", source_roundoff_upper, chain)
        observe("phase_gap", phase_gap_upper, chain)
        observe("total_gap", total_gap_upper, chain)
        observe("lipschitz", lipschitz_upper, chain)
        if minimum_magnitude is None or true_magnitude_lower < minimum_magnitude[0]:
            minimum_magnitude = (true_magnitude_lower, chain)

        header = direct[chain]["header"]
        level = direct[chain]["levels"][1]
        rows.append(
            {
                "chain": chain,
                "sum_index": int(header["sum_index"]),
                "branch": int(header["branch"]),
                "length_upper_index": 104,
                "term_count": 105,
                "conjugated": bool(level["conjugate"]),
                "source_roundoff_gap": complex_record(high["source_roundoff_gap"]),
                "source_roundoff_gap_abs_upper": fraction_decimal(source_roundoff_upper),
                "phase_normalization_gap": complex_record(high["phase_normalization_gap"]),
                "phase_normalization_gap_abs_upper": fraction_decimal(phase_gap_upper),
                "total_true_gap": complex_record(high["total_true_gap"]),
                "total_true_gap_abs_upper": fraction_decimal(total_gap_upper),
                "phase_lipschitz_upper": fraction_decimal(lipschitz_upper),
                "true_kernel": complex_record(high["true_kernel"]),
                "true_kernel_magnitude_lower": fraction_decimal(true_magnitude_lower),
                "precision_overlap": True,
            }
        )

    require(minimum_magnitude is not None, "empty direct-kernel magnitude roster")
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_direct_kernel_gate",
        "status": "rigorous_50_call_true_two_pi_direct_kernel_audit",
        "scope": {
            "block": 20,
            "direct_chain_count": 50,
            "upper_index": 104,
            "terms_per_kernel": 105,
            "precision_ladder_decimal_digits": list(PRECISIONS),
            "precision_overlap_count": 50,
            "logged_state_contract": "Each ES50.40E4 component is enclosed by one full last-decimal unit.",
        },
        "phase_definition": {
            "source": "p=4*atan(1); tpp=2*p; tpm_hat=-tpp",
            "mathematical": "tpm=-2*pi",
            "finite_sum": "sum_(n=0)^104 exp(i*tpm*(a1*n+a2*n^2+a3*n^3))",
            "analytic_budget": "sum_n |(tpm-tpm_hat)*(a1*n+a2*n^2+a3*n^3)|",
        },
        "aggregate": {
            "maximum_source_roundoff_gap_abs_upper": fraction_decimal(maxima["source_roundoff"][0]),
            "maximum_source_roundoff_gap_witness": maxima["source_roundoff"][1],
            "maximum_phase_normalization_gap_abs_upper": fraction_decimal(maxima["phase_gap"][0]),
            "maximum_phase_normalization_gap_witness": maxima["phase_gap"][1],
            "maximum_total_true_gap_abs_upper": fraction_decimal(maxima["total_gap"][0]),
            "maximum_total_true_gap_witness": maxima["total_gap"][1],
            "maximum_phase_lipschitz_upper": fraction_decimal(maxima["lipschitz"][0]),
            "maximum_phase_lipschitz_witness": maxima["lipschitz"][1],
            "minimum_true_kernel_magnitude_lower": fraction_decimal(minimum_magnitude[0]),
            "minimum_true_kernel_magnitude_witness": minimum_magnitude[1],
            "all_phase_gaps_within_lipschitz_budget": True,
        },
        "rows": rows,
        "source": {"path": relative(SOURCE), "sha256": file_hash(SOURCE), "locations": source_locations()},
        "sources": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "selector_atlas": {"path": relative(ATLAS), "sha256": file_hash(ATLAS)},
            "initial_adapter": {"path": relative(INITIAL_ADAPTER), "sha256": file_hash(INITIAL_ADAPTER)},
            "point_ball_builder": {
                "path": relative(Path(point_balls.__file__).resolve()),
                "sha256": file_hash(Path(point_balls.__file__).resolve()),
            },
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "python_flint": {"version": point_balls.flint.__version__, "flint_version": point_balls.flint.__FLINT_VERSION__},
        },
        "next_handoff": {
            "recursive": "Complete the real-erfc, saddle-sum, and mathematical truncation audit inside t5 for the 374 MIT=2 calls.",
            "accumulation": "After the recursive q audit, combine all 424 certified local values through the finite block sum while keeping the outer Hardy remainder separate.",
        },
        "proof_boundary": (
            "Rigorous exact-point finite-sum audit for the 50 direct block-20 calls only. It does not "
            "transport continuous coefficient error, bound the recursive t5/q remainder, control block "
            "accumulation or the outer Hardy representation, establish height uniformity, Lambda <= 0, "
            "PF-infinity, RH, or a prize-level conclusion."
        ),
    }
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(
        "built Hardy block-20 direct-kernel gate: "
        f"50 kernels, max true/logged gap {artifact['aggregate']['maximum_total_true_gap_abs_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
