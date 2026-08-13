#!/usr/bin/env python3
"""Build exact-binary128 Arb point-ball certificates for saved Hardy defects."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_FLINT_ROOT = Path(
    r"C:\Users\ollet\Documents\Codex\third_party\python_flint_0_8_0"
)
FLINT_ROOT = Path(os.environ.get("RH_PYTHON_FLINT_ROOT", DEFAULT_FLINT_ROOT))
sys.path.insert(0, str(FLINT_ROOT))

try:
    import flint
    from flint import acb, arb, ctx
except ImportError as exc:  # pragma: no cover - explicit environment failure
    raise RuntimeError(f"python-flint runtime unavailable at {FLINT_ROOT}") from exc


TELEMETRY = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10/enabled/chain_telemetry.jsonl"
)
SCOUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_chain_telemetry_scout.json"
)
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_point_ball_defect_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_point_ball_defect_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_point_ball_defect_gate.py"
PRECISION_LADDER_DPS = (96, 160)
HEX128 = re.compile(r"[0-9A-F]{32}")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tree_hash(root: Path) -> tuple[str, int, int]:
    digest = hashlib.sha256()
    count = 0
    size = 0
    for path in sorted(root.rglob("*")):
        if not path.is_file() or "__pycache__" in path.parts or path.suffix == ".pyc":
            continue
        relative_path = path.relative_to(root).as_posix().encode("utf-8")
        payload = path.read_bytes()
        digest.update(relative_path)
        digest.update(b"\0")
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
        count += 1
        size += len(payload)
    return digest.hexdigest(), count, size


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    records = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
    require(records, "empty telemetry input")
    return records


def organize(records: list[dict[str, Any]]) -> dict[int, dict[str, Any]]:
    chains: dict[int, dict[str, Any]] = {}
    for record in records:
        record_type = record["type"]
        if record_type == "chain":
            chain_id = int(record["chain"])
            chains[chain_id] = {"header": record, "levels": {}, "steps": {}}
        elif record_type == "level":
            chains[int(record["chain"])]["levels"][int(record["level"])] = record
        elif record_type == "recurrence":
            chains[int(record["chain"])]["steps"][int(record["nit"])] = record
    require(chains, "no telemetry chains")
    return chains


def binary128_ball(payload: str) -> arb:
    require(isinstance(payload, str) and HEX128.fullmatch(payload) is not None, "invalid binary128 payload")
    bits = int(payload, 16)
    sign = -1 if bits >> 127 else 1
    exponent = (bits >> 112) & 0x7FFF
    fraction = bits & ((1 << 112) - 1)
    require(exponent != 0x7FFF, "non-finite binary128 input")
    if exponent == 0:
        if fraction == 0:
            return arb(0)
        mantissa = fraction
        binary_exponent = 1 - 16383 - 112
    else:
        mantissa = (1 << 112) + fraction
        binary_exponent = exponent - 16383 - 112
    return arb((sign * mantissa, binary_exponent))


def binary128_complex(payload: list[str]) -> acb:
    require(isinstance(payload, list) and len(payload) == 2, "invalid complex binary128 payload")
    return acb(binary128_ball(payload[0]), binary128_ball(payload[1]))


def direct_sum(level: dict[str, Any], tpm: arb) -> acb:
    coefficients = [binary128_ball(value) for value in level["coefficients_hex"]]
    require(len(coefficients) == int(level["degree"]), "coefficient/degree mismatch")
    length = int(level["length"])
    require(0 <= length <= 100000, "direct point-ball length outside scout bound")
    total = acb(0)
    for index in range(length + 1):
        phase = arb(0)
        for coefficient in reversed(coefficients):
            phase = (phase + coefficient) * index
        total += acb(arb(0), tpm * phase).exp()
    return total


def adapt(level: dict[str, Any], value: acb) -> acb:
    if level["conjugate"]:
        value = value.conjugate()
    if level["subtract_one"]:
        value -= 1
    return value


def apply_step(step: dict[str, Any], child: acb) -> acb:
    value = binary128_complex(step["multiplier_hex"]) * child
    value += binary128_complex(step["qq_hex"])
    if step["conjugate"]:
        value = value.conjugate()
    if step["subtract_one"]:
        value -= 1
    return value


def dyadic(value: arb) -> list[int]:
    mantissa, exponent = value.man_exp()
    return [int(mantissa), int(exponent)]


def arb_record(value: arb, digits: int = 45) -> dict[str, Any]:
    lower = value.lower()
    upper = value.upper()
    return {
        "display": value.str(digits),
        "lower_dyadic": dyadic(lower),
        "upper_dyadic": dyadic(upper),
        "lower_decimal": lower.str(digits, radius=False),
        "upper_decimal": upper.str(digits, radius=False),
    }


def acb_record(value: acb, digits: int = 40) -> dict[str, Any]:
    return {
        "display": value.str(digits),
        "real": arb_record(value.real, digits),
        "imag": arb_record(value.imag, digits),
        "abs": arb_record(abs(value), digits),
        "contains_zero": value.contains(0),
        "relative_accuracy_bits": int(value.rel_accuracy_bits()),
    }


def radius_upper(value: acb) -> arb:
    real_radius = value.real.rad()
    imag_radius = value.imag.rad()
    return real_radius if float(real_radius) >= float(imag_radius) else imag_radius


def compute_chain(data: dict[str, Any], dps: int) -> dict[str, Any]:
    ctx.dps = dps
    ctx.threads = 1
    header = data["header"]
    levels = data["levels"]
    steps = data["steps"]
    require(int(header["mit"]) == 2, "point-ball gate currently expects the saved MIT=2 scout")
    require(set(levels) == {1, 2} and set(steps) == {1}, "saved one-step chain shape drift")

    tpm = binary128_ball(header["tpm_hex"])
    parent_raw = direct_sum(levels[1], tpm)
    child_raw = direct_sum(levels[2], tpm)
    parent = adapt(levels[1], parent_raw)
    child = adapt(levels[2], child_raw)
    model = apply_step(steps[1], child)
    defect = parent - model
    logged_defect = binary128_complex(steps[1]["local_defect_hex"])
    source_roundoff_gap = defect - logged_defect

    return {
        "parent_raw": parent_raw,
        "child_raw": child_raw,
        "parent": parent,
        "child": child,
        "model": model,
        "defect": defect,
        "logged_defect": logged_defect,
        "source_roundoff_gap": source_roundoff_gap,
    }


def choose(values: list[arb], *, maximum: bool) -> arb:
    require(values, "empty Arb aggregate")
    return (max if maximum else min)(values, key=lambda value: float(value.mid()))


def summarize(values: list[arb]) -> dict[str, Any]:
    require(values, "empty Arb summary")
    ordered = sorted(values, key=lambda value: float(value.mid()))
    return {
        "min": arb_record(ordered[0]),
        "median": arb_record(ordered[(len(ordered) - 1) // 2]),
        "p90": arb_record(ordered[round((len(ordered) - 1) * 0.9)]),
        "max": arb_record(ordered[-1]),
    }


def build_artifact() -> dict[str, Any]:
    require(FLINT_ROOT.is_dir(), f"python-flint runtime missing: {FLINT_ROOT}")
    require(flint.__version__ == "0.8.0", "python-flint version drift")
    require(getattr(flint, "__FLINT_VERSION__", None) == "3.3.1", "FLINT version drift")
    require(TELEMETRY.is_file() and SCOUT.is_file(), "point-ball input missing")
    require(CHECKER.is_file(), "point-ball checker missing")

    scout = json.loads(SCOUT.read_text(encoding="utf-8"))
    require(scout["aggregate"]["chain_count"] == 64, "scout chain count drift")
    chains = organize(load_jsonl(TELEMETRY))
    require(sorted(chains) == list(range(1, 65)), "telemetry chain sequence drift")

    low_results = {
        chain_id: compute_chain(chains[chain_id], PRECISION_LADDER_DPS[0])
        for chain_id in sorted(chains)
    }
    high_results = {
        chain_id: compute_chain(chains[chain_id], PRECISION_LADDER_DPS[1])
        for chain_id in sorted(chains)
    }

    rows: list[dict[str, Any]] = []
    defect_abs_values: list[arb] = []
    gap_abs_values: list[arb] = []
    gap_ratio_values: list[arb] = []
    high_radii: list[arb] = []
    low_radii: list[arb] = []
    overlap_count = 0
    zero_exclusion_count = 0

    for chain_id in sorted(chains):
        header = chains[chain_id]["header"]
        low = low_results[chain_id]
        high = high_results[chain_id]
        overlaps = high["defect"].overlaps(low["defect"])
        overlap_count += int(overlaps)
        excludes_zero = not high["defect"].contains(0)
        zero_exclusion_count += int(excludes_zero)

        defect_abs = abs(high["defect"])
        gap_abs = abs(high["source_roundoff_gap"])
        gap_ratio = gap_abs / defect_abs
        high_radius = radius_upper(high["defect"])
        low_radius = radius_upper(low["defect"])
        defect_abs_values.append(defect_abs)
        gap_abs_values.append(gap_abs)
        gap_ratio_values.append(gap_ratio)
        high_radii.append(high_radius)
        low_radii.append(low_radius)

        rows.append(
            {
                "chain": chain_id,
                "block": int(header["block"]),
                "sum_index": int(header["sum_index"]),
                "branch": int(header["branch"]),
                "initial_length": int(header["initial_length"]),
                "kernel_length": int(header["kernel_length"]),
                "tpm_hex": header["tpm_hex"],
                "defect": acb_record(high["defect"]),
                "logged_binary128_defect": acb_record(high["logged_defect"]),
                "source_roundoff_gap": acb_record(high["source_roundoff_gap"]),
                "source_roundoff_gap_to_defect_ratio": arb_record(gap_ratio),
                "high_precision_component_radius": arb_record(high_radius),
                "low_precision_component_radius": arb_record(low_radius),
                "precision_balls_overlap": overlaps,
                "defect_excludes_zero": excludes_zero,
            }
        )

    require(overlap_count == len(rows), "precision ladder ball overlap failure")
    require(zero_exclusion_count == len(rows), "saved point defect contains zero")
    require(
        all(float(high.mid()) < float(low.mid()) for high, low in zip(high_radii, low_radii)),
        "higher precision did not shrink every defect radius",
    )

    tpm_payloads = {chains[chain_id]["header"]["tpm_hex"] for chain_id in chains}
    require(len(tpm_payloads) == 1, "tpm binary payload drift across chains")
    tpm_hex = next(iter(tpm_payloads))
    runtime_sha256, runtime_files, runtime_bytes = tree_hash(FLINT_ROOT)

    aggregate = {
        "chain_count": len(rows),
        "branch_counts": {
            str(branch): sum(row["branch"] == branch for row in rows) for branch in (1, 2)
        },
        "precision_overlap_count": overlap_count,
        "defect_zero_exclusion_count": zero_exclusion_count,
        "defect_abs": summarize(defect_abs_values),
        "source_roundoff_gap_abs": summarize(gap_abs_values),
        "source_roundoff_gap_to_defect_ratio": summarize(gap_ratio_values),
        "high_precision_component_radius": summarize(high_radii),
        "low_precision_component_radius": summarize(low_radii),
        "minimum_defect_abs_lower_bound": arb_record(
            choose([value.lower() for value in defect_abs_values], maximum=False)
        ),
        "maximum_source_roundoff_gap_upper_bound": arb_record(
            choose([value.upper() for value in gap_abs_values], maximum=True)
        ),
        "maximum_gap_to_defect_ratio_upper_bound": arb_record(
            choose([value.upper() for value in gap_ratio_values], maximum=True)
        ),
    }

    return {
        "kind": "jensen_window_pf_newman_c1_hardy_point_ball_defect_gate",
        "status": "exact_binary128_point_ball_certificate",
        "runtime": {
            "python_flint_version": flint.__version__,
            "flint_version": flint.__FLINT_VERSION__,
            "flint_release": flint.__FLINT_RELEASE__,
            "root": str(FLINT_ROOT),
            "tree_sha256": runtime_sha256,
            "file_count": runtime_files,
            "byte_count": runtime_bytes,
            "threads": 1,
            "precision_ladder_dps": list(PRECISION_LADDER_DPS),
        },
        "input_exactness": {
            "format": "IEEE-754 binary128 payload encoded as 32 uppercase hexadecimal digits",
            "decoder": "sign bit 127; 15-bit exponent with bias 16383; 112 explicit fraction bits",
            "arb_constructor": "arb((signed_integer_mantissa, binary_exponent))",
            "decimal_strings_used_for_arithmetic": False,
            "tpm_hex": tpm_hex,
            "tpm_exact_dyadic": dyadic(binary128_ball(tpm_hex)),
            "source_phase_definition": "p=4*ATAN(1); tpp=2*p; tpm=-tpp",
        },
        "finite_target": {
            "raw_level_sum": "sum_(k=0)^L exp(i*tpm*sum_(j=1)^m phi_j*k^j)",
            "level_adapter": "Apply the logged source conjugation flag, then the logged real subtraction flag.",
            "local_defect": "adapted_parent - adapt_k(multiplier*adapted_child + qq)",
            "q_input": "The exact logged binary128 source qq; no analytic W1-W5 truth enclosure is claimed.",
        },
        "aggregate": aggregate,
        "rows": rows,
        "route_conclusion": {
            "finding": (
                "All 64 exact-input source-convention finite defect balls exclude zero. "
                "The discrepancy is therefore not explained by binary128 operation rounding."
            ),
            "preferred_next_route": "interval_direct_local_defect",
            "next_requirement": (
                "Promote coefficients, multiplier, and analytic qq from exact point payloads "
                "to outward-rounded functions on a real cell or fixed complex disk, while "
                "proving the transformed-level adapter and every selector margin."
            ),
        },
        "sources": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "scout": {"path": relative(SCOUT), "sha256": file_hash(SCOUT)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "proof_boundary": (
            "This gate rigorously encloses a source-convention finite expression at 64 saved "
            "exact binary128 input points and proves those point defect balls exclude zero. "
            "It does not prove that the source-convention transformed levels equal the paper's "
            "exact equation-120 objects, enclose the analytic W1-W5 corrections, control any "
            "real interval or complex disk, exclude selector or hierarchy transitions, bound "
            "the outer Hardy representation error, establish physical-height complexity, or "
            "imply a carrier value, Lambda<=0, PF-infinity, RH, or a prize-level conclusion."
        ),
    }


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return "\n".join(
        [
            "# Hardy Exact-Binary128 Point-Ball Defect Gate",
            "",
            "Date: 2026-08-06",
            "Status: rigorous saved-point source-expression certificate; not a proof and no uniform analytic error theorem",
            "",
            "## Exact Inputs",
            "",
            "The telemetry records each proof-relevant `real(kind=16)` value as a 32-hex-digit payload. The gate decodes each payload as an exact IEEE-754 binary128 dyadic rational and constructs the corresponding exact Arb point. Decimal telemetry is not used for arithmetic.",
            "",
            "The phase constant remains source-derived: `p=4*ATAN(1)`, `tpp=2*p`, `tpm=-tpp`. Its saved binary payload is decoded directly; no separate value of pi is inserted.",
            "",
            "## Certificate",
            "",
            f"At `{artifact['runtime']['precision_ladder_dps'][0]}` and `{artifact['runtime']['precision_ladder_dps'][1]}` decimal digits, all `{aggregate['precision_overlap_count']}` low/high precision defect balls overlap, and all `{aggregate['defect_zero_exclusion_count']}` high-precision balls exclude zero.",
            "",
            f"The smallest rigorous defect-magnitude lower endpoint is `{aggregate['minimum_defect_abs_lower_bound']['lower_decimal']}`. The largest upper endpoint for the gap between exact-real Arb evaluation and the logged binary128 computation is `{aggregate['maximum_source_roundoff_gap_upper_bound']['upper_decimal']}`. The largest gap-to-defect ratio upper endpoint is `{aggregate['maximum_gap_to_defect_ratio_upper_bound']['upper_decimal']}`.",
            "",
            "This establishes that the observed source-convention local discrepancy is genuine at the saved points and is not a binary128 accumulation artifact. It strengthens the case for enclosing the whole finite defect directly, where the W1-W5 cancellation is retained.",
            "",
            "## Next Target",
            "",
            artifact["route_conclusion"]["next_requirement"],
            "",
            "## Proof Boundary",
            "",
            artifact["proof_boundary"],
            "",
        ]
    )


def main() -> int:
    artifact = build_artifact()
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    result_tmp = RESULT.with_suffix(".json.tmp")
    note_tmp = NOTE.with_suffix(".md.tmp")
    result_tmp.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    note_tmp.write_text(render_note(artifact), encoding="utf-8")
    result_tmp.replace(RESULT)
    note_tmp.replace(NOTE)
    print(
        "built Hardy point-ball defect gate: 64 exact binary128 chains, "
        "64 zero exclusions, 2 precision levels and 0 uniform bounds"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
