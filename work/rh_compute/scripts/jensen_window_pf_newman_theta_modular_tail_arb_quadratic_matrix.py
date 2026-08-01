#!/usr/bin/env python3
"""Build a resumable Arb matrix for compact modular-tail derivative budgets."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import comb, factorial
import os
from pathlib import Path
import sys
import time

try:
    import psutil

    psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
except Exception:
    pass

REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR_DIR = REPO_ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR_DIR) not in sys.path:
    sys.path.insert(0, str(VENDOR_DIR))

import flint
from flint import acb, arb
import sympy as sp


STEM = "jensen_window_pf_newman_theta_modular_tail_arb_quadratic_matrix"
RESULT_DIR = REPO_ROOT / "work" / "rh_compute" / "results"
DEFAULT_CACHE = RESULT_DIR / f"{STEM}.jsonl"
DEFAULT_OUT = RESULT_DIR / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
ENVELOPE_SOURCE = (
    RESULT_DIR
    / "jensen_window_pf_newman_theta_modular_tail_derivative_envelope_gate.json"
)
PILOT_SOURCE = (
    RESULT_DIR
    / "jensen_window_pf_newman_theta_modular_tail_arb_quadratic_pilot.json"
)

SCHEMA = "newman_theta_modular_tail_arb_quadratic_matrix_v1"
DATE = "2026-07-24"
PRECISION_BITS = 128
T_CAP = "0.2"
U_STOP = "2.2"
N_CAP = 24
N_VALUES = tuple(range(4, 11))
M_ORDER = 9
INTEGRATION_OPTIONS = {
    "deg_limit": 80,
    "eval_limit": 100_000,
    "depth_limit": 30,
    "use_heap": True,
}


def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def matrix_specs() -> list[dict]:
    specs: list[dict] = []
    for b in range(M_ORDER + 1):
        specs.append(
            {
                "component": "d0",
                "p": 0,
                "b": b,
                "q": M_ORDER - b,
                "coefficient": comb(M_ORDER, b),
            }
        )
        specs.append(
            {
                "component": "d1_u",
                "p": 1,
                "b": b,
                "q": M_ORDER - b,
                "coefficient": comb(M_ORDER, b),
            }
        )
    for b in range(M_ORDER):
        specs.append(
            {
                "component": "d1_lower",
                "p": 0,
                "b": b,
                "q": M_ORDER - 1 - b,
                "coefficient": M_ORDER * comb(M_ORDER - 1, b),
            }
        )
    return specs


def task_id(task: dict) -> str:
    if task["kind"] == "mass":
        return f"mass:p{task['p']}:b{task['b']}"
    return (
        f"quadratic:N{task['N']}:p{task['p']}:"
        f"b{task['b']}:q{task['q']}"
    )


def expected_tasks() -> list[dict]:
    specs = matrix_specs()
    weights = sorted({(row["p"], row["b"]) for row in specs})
    tasks = [
        {"kind": "mass", "p": p, "b": b}
        for p, b in weights
    ]
    for retained in N_VALUES:
        for row in specs:
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


def config_payload() -> dict:
    return {
        "schema": SCHEMA,
        "precision_bits": PRECISION_BITS,
        "t_cap_exact": "1/5",
        "u_interval_exact": ["0", "11/5"],
        "n_cap": N_CAP,
        "n_values": list(N_VALUES),
        "m_order": M_ORDER,
        "integration_options": INTEGRATION_OPTIONS,
        "envelope_source_sha256": file_hash(ENVELOPE_SOURCE),
        "pilot_source_sha256": file_hash(PILOT_SOURCE),
        "task_ids": [task_id(task) for task in expected_tasks()],
    }


def config_hash() -> str:
    return sha256(canonical_json(config_payload()).encode("utf-8")).hexdigest()


def load_cache(path: Path) -> tuple[list[dict], dict[str, dict]]:
    if not path.exists():
        return [], {}
    records: list[dict] = []
    by_id: dict[str, dict] = {}
    previous = "0" * 64
    expected_config = config_hash()
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.endswith("\n"):
                raise RuntimeError(
                    f"cache line {line_number} is not fsynced/terminated"
                )
            record = json.loads(line)
            stored_hash = record.pop("row_sha256", None)
            actual_hash = sha256(
                canonical_json(record).encode("utf-8")
            ).hexdigest()
            record["row_sha256"] = stored_hash
            if stored_hash != actual_hash:
                raise RuntimeError(f"cache hash mismatch at line {line_number}")
            if record.get("sequence") != line_number:
                raise RuntimeError(f"cache sequence mismatch at line {line_number}")
            if record.get("previous_row_sha256") != previous:
                raise RuntimeError(f"cache chain mismatch at line {line_number}")
            if record.get("config_sha256") != expected_config:
                raise RuntimeError(f"cache config mismatch at line {line_number}")
            identifier = task_id(record["task"])
            if identifier in by_id:
                raise RuntimeError(f"duplicate cache task: {identifier}")
            records.append(record)
            by_id[identifier] = record
            previous = stored_hash
    return records, by_id


def append_record(
    path: Path,
    records: list[dict],
    task: dict,
    result: acb,
) -> dict:
    real = result.real
    record = {
        "schema": SCHEMA,
        "sequence": len(records) + 1,
        "config_sha256": config_hash(),
        "task": task,
        "real_enclosure": real.str(45, more=True),
        "real_lower": real.lower().str(45),
        "real_upper": real.upper().str(45),
        "imaginary_enclosure": result.imag.str(25, more=True),
        "relative_accuracy_bits": real.rel_accuracy_bits(),
        "previous_row_sha256": (
            records[-1]["row_sha256"] if records else "0" * 64
        ),
    }
    record["row_sha256"] = sha256(
        canonical_json(record).encode("utf-8")
    ).hexdigest()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(canonical_json(record) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
    records.append(record)
    return record


def heat_polynomial(order: int, z: acb) -> acb:
    t_cap = acb(T_CAP)
    value = acb(0)
    for ell in range(order // 2 + 1):
        value += (
            factorial(order)
            * (2 * z) ** (order - 2 * ell)
            * t_cap ** (order - ell)
            / (factorial(ell) * factorial(order - 2 * ell))
        )
    return value


def weight(z: acb, p: int, b: int) -> acb:
    return (
        z**p
        * (acb(T_CAP) * z * z).exp()
        * heat_polynomial(b, z)
    )


def build_derivative_functions() -> list:
    u, n = sp.symbols("u n", positive=True, real=True)

    def phi(z: sp.Expr) -> sp.Expr:
        return (
            2 * sp.pi**2 * n**4 * sp.exp(9 * z)
            - 3 * sp.pi * n**2 * sp.exp(5 * z)
        ) * sp.exp(-sp.pi * n**2 * sp.exp(4 * z))

    omega = (1 + sp.erf(3 * sp.sinh(4 * u))) / 2
    current = omega * phi(u) + (1 - omega) * phi(-u)
    expressions: list[sp.Expr] = []
    for _ in range(M_ORDER + 1):
        expressions.append(current)
        current = sp.diff(current, u)
    functions = {
        "exp": lambda z: z.exp(),
        "sinh": lambda z: z.sinh(),
        "cosh": lambda z: z.cosh(),
        "erf": lambda z: z.erf(),
        "sqrt": lambda z: z.sqrt(),
        "pi": acb.pi(),
    }
    return [
        sp.lambdify(
            (u, n),
            expression,
            modules=[functions],
            cse=True,
            docstring_limit=0,
        )
        for expression in expressions
    ]


def integrate_task(task: dict, derivatives: list | None) -> acb:
    if task["kind"] == "mass":
        integrand = lambda z, _analytic: weight(z, task["p"], task["b"])
    else:
        if derivatives is None:
            raise RuntimeError("quadratic task requires derivative functions")

        def integrand(z: acb, _analytic: bool) -> acb:
            remainder = sum(
                (
                    derivatives[task["q"]](z, index)
                    for index in range(task["N"] + 1, N_CAP + 1)
                ),
                acb(0),
            )
            return weight(z, task["p"], task["b"]) * remainder * remainder

    options = dict(INTEGRATION_OPTIONS)
    for key in ("abs_tol", "rel_tol"):
        value = options.get(key)
        if isinstance(value, str) and value.startswith("2^-"):
            options[key] = arb(2) ** (-int(value[3:]))

    result = acb.integral(
        integrand,
        arb(0),
        arb(U_STOP),
        **options,
    )
    if not result.is_finite():
        raise RuntimeError(f"non-finite integral for {task_id(task)}")
    if task["kind"] == "mass" and result.real.lower() <= 0:
        raise RuntimeError(
            f"weight mass was not separated for {task_id(task)}: {result}"
        )
    if task["kind"] == "quadratic" and result.real.upper() <= 0:
        raise RuntimeError(
            f"quadratic upper bound was not positive for "
            f"{task_id(task)}: {result}"
        )
    if not result.imag.contains(0):
        raise RuntimeError(
            f"imaginary enclosure excludes zero for {task_id(task)}"
        )
    return result


def run_tasks(
    cache: Path,
    only_n: set[int] | None,
    max_new_tasks: int | None,
    max_seconds: float | None,
    stop_file: Path | None,
) -> tuple[list[dict], dict[str, dict], int]:
    flint.ctx.prec = PRECISION_BITS
    records, by_id = load_cache(cache)
    pending: list[dict] = []
    for task in expected_tasks():
        identifier = task_id(task)
        if identifier in by_id:
            continue
        if (
            task["kind"] == "quadratic"
            and only_n is not None
            and task["N"] not in only_n
        ):
            continue
        pending.append(task)

    started = time.monotonic()
    derivatives: list | None = None
    completed = 0
    for task in pending:
        if stop_file is not None and stop_file.exists():
            print(
                f"stop file observed before {task_id(task)}; parking",
                flush=True,
            )
            break
        if max_new_tasks is not None and completed >= max_new_tasks:
            break
        if (
            max_seconds is not None
            and completed > 0
            and time.monotonic() - started >= max_seconds
        ):
            break
        if task["kind"] == "quadratic" and derivatives is None:
            derivatives = build_derivative_functions()
        task_started = time.monotonic()
        result = integrate_task(task, derivatives)
        record = append_record(cache, records, task, result)
        by_id[task_id(task)] = record
        completed += 1
        print(
            f"completed {record['sequence']}/{len(expected_tasks())} "
            f"{task_id(task)} in {time.monotonic() - task_started:.3f}s",
            flush=True,
        )
    return records, by_id, completed


def parse_positive_upper(record: dict) -> arb:
    value = arb(record["real_upper"])
    if value.lower() <= 0:
        raise RuntimeError(
            f"cached upper endpoint is not positive: "
            f"{task_id(record['task'])}"
        )
    return value


def aggregate_n(retained: int, by_id: dict[str, dict]) -> dict | None:
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
        for row in matrix_specs()
    ]
    if any(identifier not in by_id for identifier in required):
        return None

    component_balls = {
        "d0": arb(0),
        "d1_u": arb(0),
        "d1_lower": arb(0),
    }
    for spec in matrix_specs():
        mass_id = task_id(
            {"kind": "mass", "p": spec["p"], "b": spec["b"]}
        )
        quadratic_id = task_id(
            {
                "kind": "quadratic",
                "N": retained,
                "p": spec["p"],
                "b": spec["b"],
                "q": spec["q"],
            }
        )
        mass = parse_positive_upper(by_id[mass_id])
        quadratic = parse_positive_upper(by_id[quadratic_id])
        component_balls[spec["component"]] += (
            spec["coefficient"] * (mass * quadratic).sqrt()
        )
    d1 = component_balls["d1_u"] + component_balls["d1_lower"]

    def serialize(value: arb) -> dict:
        return {
            "upper_accumulator_enclosure": value.str(40, more=True),
            "upper": value.upper().str(40),
            "relative_accuracy_bits": value.rel_accuracy_bits(),
        }

    return {
        "N": retained,
        "finite_compact_d0": serialize(component_balls["d0"]),
        "finite_compact_d1_u": serialize(component_balls["d1_u"]),
        "finite_compact_d1_lower": serialize(
            component_balls["d1_lower"]
        ),
        "finite_compact_d1": serialize(d1),
        "quadratic_entries": len(matrix_specs()),
    }


def build_artifact(
    cache: Path,
    records: list[dict],
    by_id: dict[str, dict],
) -> dict:
    tasks = expected_tasks()
    complete_rows = [
        row
        for retained in N_VALUES
        if (row := aggregate_n(retained, by_id)) is not None
    ]
    missing = [
        task_id(task)
        for task in tasks
        if task_id(task) not in by_id
    ]
    complete = not missing
    return {
        "kind": STEM,
        "date": DATE,
        "status": (
            "complete finite compact Arb matrix"
            if complete
            else "partial resumable finite compact Arb matrix"
        ),
        "config": config_payload(),
        "config_sha256": config_hash(),
        "cache": {
            "path": str(cache.relative_to(REPO_ROOT)).replace("\\", "/"),
            "records": len(records),
            "expected_records": len(tasks),
            "last_row_sha256": (
                records[-1]["row_sha256"] if records else None
            ),
            "complete": complete,
            "missing_count": len(missing),
            "missing_task_ids": missing,
        },
        "matrix_definition": {
            "order": M_ORDER,
            "entry_count_per_N": len(matrix_specs()),
            "mass_count": len(
                {(row["p"], row["b"]) for row in matrix_specs()}
            ),
            "d0": (
                "sum_(b=0)^9 binom(9,b) sqrt(M_(0,b)Q_(N,0,b,9-b))"
            ),
            "d1_u": (
                "sum_(b=0)^9 binom(9,b) sqrt(M_(1,b)Q_(N,1,b,9-b))"
            ),
            "d1_lower": (
                "9*sum_(b=0)^8 binom(8,b) "
                "sqrt(M_(0,b)Q_(N,0,b,8-b))"
            ),
            "mass": (
                "M_(p,b)=integral_0^(11/5) "
                "u^p exp(u^2/5)P_b(1/5,u)du"
            ),
            "quadratic": (
                "Q_(N,p,b,q)=integral_0^(11/5) "
                "u^p exp(u^2/5)P_b(1/5,u)"
                "[D^q sum_(n=N+1)^24 b_n(u)]^2du"
            ),
        },
        "finite_compact_rows": complete_rows,
        "proof_boundary": (
            "Every stored ball is a directed-rounding certificate for the "
            "finite arithmetic block n=N+1..24 and compact interval "
            "0<=u<=11/5. The artifact omits n>=25, u>11/5, retained J/J' "
            "lower balls, transition-cell coverage, Lambda<=0, RH, and any "
            "Clay-prize conclusion."
        ),
        "sources": {
            "derivative_envelope": str(
                ENVELOPE_SOURCE.relative_to(REPO_ROOT)
            ).replace("\\", "/"),
            "pilot": str(PILOT_SOURCE.relative_to(REPO_ROOT)).replace(
                "\\", "/"
            ),
        },
        "versions": {
            "python_flint": getattr(flint, "__version__", "unknown"),
            "sympy": sp.__version__,
        },
    }


def render_note(artifact: dict) -> str:
    cache = artifact["cache"]
    rows = artifact["finite_compact_rows"]
    title = (
        "# Newman Theta Modular-Tail Arb Quadratic Refinement"
        if "refinement" in artifact["kind"]
        else "# Newman Theta Modular-Tail Arb Quadratic Matrix"
    )
    lines = [
        title,
        "",
        f"Date: {DATE}",
        "",
        "Status: "
        + (
            "complete rigorous finite compact matrix."
            if cache["complete"]
            else "partial resumable rigorous finite compact matrix."
        ),
        "This is not an infinite-tail certificate and not a proof of",
        "`Lambda<=0`, RH, or a Clay-prize result.",
        "",
        "## Matrix",
        "",
        "For `m=9`, weighted Cauchy-Schwarz reduces the compact pieces of",
        "`d0` and `d1` to 29 analytic quadratic integrals per retained `N`.",
        "The cache is append-only, hash-chained, and fsynced after every",
        "integral.",
        "",
        "```text",
        artifact["matrix_definition"]["mass"],
        artifact["matrix_definition"]["quadratic"],
        artifact["matrix_definition"]["d0"],
        artifact["matrix_definition"]["d1_u"],
        artifact["matrix_definition"]["d1_lower"],
        "```",
        "",
        "## Progress",
        "",
        f"Cache rows: `{cache['records']}/{cache['expected_records']}`.",
        f"Complete retained counts: `{[row['N'] for row in rows]}`.",
        f"Missing tasks: `{cache['missing_count']}`.",
        "",
    ]
    if rows:
        lines.extend(
            [
                "## Certified Compact Budgets",
                "",
                "| N | d0 upper | d1 upper |",
                "|---:|---:|---:|",
            ]
        )
        for row in rows:
            lines.append(
                f"| {row['N']} | "
                f"`{row['finite_compact_d0']['upper']}` | "
                f"`{row['finite_compact_d1']['upper']}` |"
            )
        lines.append("")
    lines.extend(
        [
            "## Non-Promotion Guard",
            "",
            "Every row covers only `n=N+1..24` and `0<=u<=11/5`.",
            "The arithmetic tail `n>=25`, outer interval `u>11/5`, retained",
            "first-jet lower balls, and transition-cell theorem remain open.",
            "",
        ]
    )
    return "\n".join(lines)


def write_artifact(artifact: dict, out: Path, note: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    note.write_text(render_note(artifact), encoding="utf-8")


def parse_only_n(text: str | None) -> set[int] | None:
    if text is None:
        return None
    values = {int(item.strip()) for item in text.split(",") if item.strip()}
    invalid = sorted(values.difference(N_VALUES))
    if invalid:
        raise ValueError(f"only-n values outside {N_VALUES}: {invalid}")
    return values


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument(
        "--only-n",
        help="Comma-separated retained N values; mass tasks always run.",
    )
    parser.add_argument("--max-new-tasks", type=int)
    parser.add_argument("--max-seconds", type=float)
    parser.add_argument(
        "--stop-file",
        type=Path,
        help="Park between fsynced tasks when this file exists.",
    )
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()

    only_n = parse_only_n(args.only_n)
    records, by_id, completed = run_tasks(
        args.cache,
        only_n,
        args.max_new_tasks,
        args.max_seconds,
        args.stop_file,
    )
    artifact = build_artifact(args.cache, records, by_id)
    write_artifact(artifact, args.out, args.note)
    print(
        f"matrix progress: new={completed} "
        f"cache={artifact['cache']['records']}/"
        f"{artifact['cache']['expected_records']} "
        f"complete_N={[row['N'] for row in artifact['finite_compact_rows']]}",
        flush=True,
    )
    if args.require_complete and not artifact["cache"]["complete"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
