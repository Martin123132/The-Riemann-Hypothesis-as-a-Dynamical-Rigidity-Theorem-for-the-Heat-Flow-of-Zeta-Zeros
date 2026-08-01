#!/usr/bin/env python3
"""Independently validate the stable-remainder Arb matrix and cache."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR_DIR = REPO_ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR_DIR) not in sys.path:
    sys.path.insert(0, str(VENDOR_DIR))

import flint
from flint import acb, arb
import sympy as sp


STEM = (
    "jensen_window_pf_newman_theta_stable_remainder_arb_quadratic_matrix"
)
SCHEMA = "newman_theta_stable_remainder_arb_quadratic_matrix_v1"
DEFAULT_ARTIFACT = (
    REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
)
TAIL_GATE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_forward_remainder_tail_gate.json"
)
MODULAR_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_modular_blend_gate.json"
)
N_VALUES = tuple(range(4, 11))
M_ORDER = 9
FIRST_OMITTED = 13
PRECISION_BITS = 192


def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def specs() -> list[dict]:
    rows: list[dict] = []
    for b in range(M_ORDER + 1):
        rows.append(
            {
                "component": "d0",
                "p": 0,
                "b": b,
                "q": M_ORDER - b,
                "coefficient": comb(M_ORDER, b),
            }
        )
        rows.append(
            {
                "component": "d1_u",
                "p": 1,
                "b": b,
                "q": M_ORDER - b,
                "coefficient": comb(M_ORDER, b),
            }
        )
    for b in range(M_ORDER):
        rows.append(
            {
                "component": "d1_lower",
                "p": 0,
                "b": b,
                "q": M_ORDER - 1 - b,
                "coefficient": M_ORDER * comb(M_ORDER - 1, b),
            }
        )
    return rows


def task_id(task: dict) -> str:
    if task["kind"] == "mass":
        return f"mass:p{task['p']}:b{task['b']}"
    return (
        f"quadratic:N{task['N']}:p{task['p']}:"
        f"b{task['b']}:q{task['q']}"
    )


def expected_tasks() -> list[dict]:
    rows = specs()
    weights = sorted({(row["p"], row["b"]) for row in rows})
    tasks = [
        {"kind": "mass", "p": p, "b": b} for p, b in weights
    ]
    for retained in N_VALUES:
        for row in rows:
            tasks.append(
                {
                    "kind": "quadratic",
                    "N": retained,
                    "p": row["p"],
                    "b": row["b"],
                    "q": row["q"],
                }
            )
    return tasks


def read_cache(
    path: Path, expected_config: str, issues: list[str]
) -> tuple[list[dict], dict[str, dict]]:
    records: list[dict] = []
    by_id: dict[str, dict] = {}
    previous = "0" * 64
    allowed = {task_id(task) for task in expected_tasks()}
    raw = path.read_bytes()
    if raw and not raw.endswith(b"\n"):
        issues.append("cache does not end in an fsynced line terminator")
    for line_number, line in enumerate(
        raw.decode("utf-8").splitlines(), start=1
    ):
        try:
            record = json.loads(line)
        except Exception as exc:
            issues.append(
                f"cache JSON failed at line {line_number}: {exc}"
            )
            continue
        stored_hash = record.get("row_sha256")
        unhashed = dict(record)
        unhashed.pop("row_sha256", None)
        actual_hash = sha256(
            canonical_json(unhashed).encode("utf-8")
        ).hexdigest()
        if stored_hash != actual_hash:
            issues.append(f"cache row hash mismatch at line {line_number}")
        if record.get("sequence") != line_number:
            issues.append(f"cache sequence mismatch at line {line_number}")
        if record.get("schema") != SCHEMA:
            issues.append(f"cache schema mismatch at line {line_number}")
        if record.get("config_sha256") != expected_config:
            issues.append(f"cache config mismatch at line {line_number}")
        if record.get("previous_row_sha256") != previous:
            issues.append(f"cache chain mismatch at line {line_number}")
        previous = stored_hash or ""

        task = record.get("task", {})
        try:
            identifier = task_id(task)
        except Exception as exc:
            issues.append(f"bad cache task at line {line_number}: {exc}")
            continue
        if identifier not in allowed:
            issues.append(f"unexpected cache task: {identifier}")
        if identifier in by_id:
            issues.append(f"duplicate cache task: {identifier}")
        by_id[identifier] = record
        records.append(record)

        try:
            lower = arb(record["real_lower"])
            upper = arb(record["real_upper"])
            enclosure = arb(record["real_enclosure"])
            imaginary = arb(record["imaginary_enclosure"])
            if lower > upper:
                issues.append(f"reversed real endpoints: {identifier}")
            if not enclosure.is_finite():
                issues.append(f"nonfinite diagnostic enclosure: {identifier}")
            if not imaginary.contains(0):
                issues.append(f"imaginary enclosure misses zero: {identifier}")
            if upper <= 0:
                issues.append(f"nonpositive integral upper: {identifier}")
            if task.get("kind") == "mass" and lower <= 0:
                issues.append(f"mass lower is not positive: {identifier}")
        except Exception as exc:
            issues.append(f"cache enclosure parse failed for {identifier}: {exc}")
    return records, by_id


def coefficient_norms() -> dict[int, int]:
    x = sp.symbols("X", nonnegative=True)
    current = 2 * x - 3
    norms: dict[int, int] = {}
    for q in range(M_ORDER + 1):
        polynomial = sp.Poly(sp.expand(current), x)
        norms[q] = int(
            sum(abs(value) for value in polynomial.all_coeffs())
        )
        current = sp.expand(
            (5 - 4 * x) * polynomial.as_expr()
            + 4 * x * sp.diff(polynomial.as_expr(), x)
        )
    return norms


def tail_constants() -> dict[int, arb]:
    norms = coefficient_norms()
    pi = arb.pi()
    constants: dict[int, arb] = {}
    for q in range(M_ORDER + 1):
        d = 2 * q + 4
        rho = (
            arb(d) / FIRST_OMITTED
            - pi * (2 * FIRST_OMITTED + 1)
        ).exp()
        bound = (
            norms[q]
            * pi ** (q + 2)
            * FIRST_OMITTED**d
            * (-pi * FIRST_OMITTED**2).exp()
            / (1 - rho)
        )
        constants[q] = arb(bound.upper())
    return constants


def positive_upper(record: dict) -> arb:
    value = arb(record["real_upper"])
    if value.lower() <= 0:
        raise RuntimeError(
            f"nonpositive cached upper for {task_id(record['task'])}"
        )
    return value


def serialized(value: arb) -> dict:
    return {
        "upper_accumulator_enclosure": value.str(40, more=True),
        "upper": value.upper().str(40),
        "relative_accuracy_bits": value.rel_accuracy_bits(),
    }


def reconstruct_row(
    retained: int, by_id: dict[str, dict]
) -> dict | None:
    required = [
        task_id(
            {
                "kind": "quadratic",
                "N": retained,
                "p": row["p"],
                "b": row["b"],
                "q": row["q"],
            }
        )
        for row in specs()
    ]
    if any(identifier not in by_id for identifier in required):
        return None
    constants = tail_constants()
    finite = {
        "d0": arb(0),
        "d1_u": arb(0),
        "d1_lower": arb(0),
    }
    tail = {
        "d0": arb(0),
        "d1_u": arb(0),
        "d1_lower": arb(0),
    }
    for spec in specs():
        mass = positive_upper(
            by_id[
                task_id(
                    {
                        "kind": "mass",
                        "p": spec["p"],
                        "b": spec["b"],
                    }
                )
            ]
        )
        quadratic = positive_upper(
            by_id[
                task_id(
                    {
                        "kind": "quadratic",
                        "N": retained,
                        "p": spec["p"],
                        "b": spec["b"],
                        "q": spec["q"],
                    }
                )
            ]
        )
        coefficient = spec["coefficient"]
        finite[spec["component"]] += coefficient * (
            mass * quadratic
        ).sqrt()
        tail[spec["component"]] += (
            coefficient * constants[spec["q"]] * mass
        )
    total = {key: finite[key] + tail[key] for key in finite}
    finite_d1 = finite["d1_u"] + finite["d1_lower"]
    tail_d1 = tail["d1_u"] + tail["d1_lower"]
    total_d1 = total["d1_u"] + total["d1_lower"]
    return {
        "N": retained,
        "finite_stable_d0": serialized(finite["d0"]),
        "finite_stable_d1_u": serialized(finite["d1_u"]),
        "finite_stable_d1_lower": serialized(finite["d1_lower"]),
        "finite_stable_d1": serialized(finite_d1),
        "forward_arithmetic_tail_d0": serialized(tail["d0"]),
        "forward_arithmetic_tail_d1_u": serialized(tail["d1_u"]),
        "forward_arithmetic_tail_d1_lower": serialized(
            tail["d1_lower"]
        ),
        "forward_arithmetic_tail_d1": serialized(tail_d1),
        "compact_full_d0": serialized(total["d0"]),
        "compact_full_d1_u": serialized(total["d1_u"]),
        "compact_full_d1_lower": serialized(total["d1_lower"]),
        "compact_full_d1": serialized(total_d1),
        "quadratic_entries": len(specs()),
    }


def validate(path: Path) -> list[str]:
    flint.ctx.prec = PRECISION_BITS
    issues: list[str] = []
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("kind") != STEM:
        issues.append("artifact kind mismatch")

    config = artifact.get("config", {})
    config_hash = sha256(
        canonical_json(config).encode("utf-8")
    ).hexdigest()
    if artifact.get("config_sha256") != config_hash:
        issues.append("artifact config hash mismatch")
    if config.get("schema") != SCHEMA:
        issues.append("config schema mismatch")
    if config.get("precision_bits") != PRECISION_BITS:
        issues.append("precision mismatch")
    if config.get("n_values") != list(N_VALUES):
        issues.append("retained N values mismatch")
    if config.get("forward_cap") != FIRST_OMITTED - 1:
        issues.append("forward cap mismatch")
    if config.get("first_omitted") != FIRST_OMITTED:
        issues.append("first omitted index mismatch")
    expected_ids = [task_id(task) for task in expected_tasks()]
    if config.get("task_ids") != expected_ids:
        issues.append("configured task list mismatch")
    if config.get("tail_gate_sha256") != file_hash(TAIL_GATE):
        issues.append("tail-gate config hash mismatch")
    if config.get("modular_source_sha256") != file_hash(MODULAR_SOURCE):
        issues.append("modular-source config hash mismatch")

    sources = artifact.get("sources", {})
    if sources.get("tail_gate_sha256") != file_hash(TAIL_GATE):
        issues.append("tail-gate source hash mismatch")
    if sources.get("modular_partition_sha256") != file_hash(
        MODULAR_SOURCE
    ):
        issues.append("modular source hash mismatch")

    cache_info = artifact.get("cache", {})
    cache_path = REPO_ROOT / cache_info.get("path", "")
    try:
        cache_path.resolve().relative_to(REPO_ROOT.resolve())
    except Exception:
        issues.append("cache path escapes repository")
        return issues
    if not cache_path.exists():
        issues.append("cache file is missing")
        return issues
    records, by_id = read_cache(cache_path, config_hash, issues)
    if cache_info.get("records") != len(records):
        issues.append("cache record count mismatch")
    if cache_info.get("expected_records") != len(expected_ids):
        issues.append("cache expected count mismatch")
    if cache_info.get("last_row_sha256") != (
        records[-1]["row_sha256"] if records else None
    ):
        issues.append("cache terminal hash mismatch")
    missing = [
        identifier for identifier in expected_ids if identifier not in by_id
    ]
    if cache_info.get("missing_count") != len(missing):
        issues.append("cache missing count mismatch")
    if cache_info.get("missing_task_ids") != missing:
        issues.append("cache missing task list mismatch")
    if cache_info.get("complete") != (not missing):
        issues.append("cache completion flag mismatch")

    rebuilt = [
        row
        for retained in N_VALUES
        if (row := reconstruct_row(retained, by_id)) is not None
    ]
    stored = artifact.get("compact_rows", [])
    if stored != rebuilt:
        issues.append("independent compact aggregation mismatch")

    identity = artifact.get("stable_remainder_definition", {})
    if identity.get("identity") != "r_N(u)=f_N(u)+tau(u)":
        issues.append("stable remainder identity mismatch")
    if "erfc(3*sinh(4u))/2" not in identity.get(
        "switch_tail", ""
    ):
        issues.append("stable switch-tail marker missing")
    if "phi_n(u)-b_n(u)" not in identity.get(
        "retained_defect", ""
    ):
        issues.append("retained-defect marker missing")

    omega, plus, minus = sp.symbols(
        "omega plus minus", real=True
    )
    block = omega * plus + (1 - omega) * minus
    if sp.expand((plus - block) - (1 - omega) * (plus - minus)) != 0:
        issues.append("retained-defect algebra failed")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "full arithmetic remainder",
        "u>11/5",
        "first-jet",
        "transition-cell",
        "Lambda<=0",
        "RH",
        "Clay-prize",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary guard missing: {marker}")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    args = parser.parse_args()
    issues = validate(args.artifact)
    if issues:
        for issue in issues:
            print(f"ISSUE: {issue}")
        raise SystemExit(1)
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    print(
        "validated Newman theta stable-remainder Arb matrix: "
        f"{artifact['cache']['records']}/"
        f"{artifact['cache']['expected_records']} cache rows, "
        f"{len(artifact['compact_rows'])} complete retained counts, "
        "0 issues"
    )


if __name__ == "__main__":
    main()
