#!/usr/bin/env python3
"""Certify the Q_31 slab with seven retained even modular blocks."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import sys
import time

import psutil


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_theta_compact_transversality_interval_certificate as compact  # noqa: E402
import jensen_window_pf_newman_theta_stable_remainder_outer_tail_gate as outer  # noqa: E402


STEM = "jensen_window_pf_newman_theta_modular_retained_q31_interval_certificate"
RESULT_DIR = REPO_ROOT / "work" / "rh_compute" / "results"
DEFAULT_CACHE = RESULT_DIR / f"{STEM}.jsonl"
DEFAULT_OUT = RESULT_DIR / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
FULL_BUDGET_SOURCE = (
    RESULT_DIR
    / "jensen_window_pf_newman_theta_full_derivative_budget_certificate.json"
)
MODULAR_SOURCE = (
    RESULT_DIR / "jensen_window_pf_newman_theta_modular_blend_gate.json"
)
OUTER_SOURCE = (
    RESULT_DIR
    / "jensen_window_pf_newman_theta_stable_remainder_outer_tail_gate.json"
)
COMPACT_SOURCE = compact.DEFAULT_OUT
DIAGONAL_SOURCE = (
    RESULT_DIR
    / "jensen_window_pf_newman_positive_boundary_diagonal_exhaustion_gate.json"
)
WINDING_SOURCE = (
    RESULT_DIR / "jensen_window_pf_newman_first_jet_winding_gate.json"
)
BOUNDARY_SOURCE = (
    RESULT_DIR
    / "jensen_window_pf_newman_positive_boundary_attainment_lemma.json"
)

SCHEMA = "newman_theta_modular_retained_q31_interval_v1"
DATE = "2026-07-24"
PRECISION_BITS = 192
RETAINED_N = 7
DIAGONAL_J = 31
TIME_LOWER = Fraction(1, 5 * DIAGONAL_J)
TIME_UPPER = Fraction(1, 5)
X_LOWER = Fraction(38)
X_UPPER = Fraction(38 + DIAGONAL_J)
INITIAL_TIME_STEP = Fraction(1, 100)
INITIAL_X_STEP = Fraction(1, 2)
MAX_DEPTH = 7
ENDPOINT_DIGITS = 55
INTEGRATION_CUTOFF = Fraction(11, 5)
INTEGRATION_OPTIONS = {
    "deg_limit": 80,
    "eval_limit": 250_000,
    "depth_limit": 50,
    "use_heap": True,
    "abs_tol": "2^-150",
    "rel_tol": "2^-145",
}
CPU_PARK_THRESHOLD = 85.0
CPU_CONSECUTIVE_PARK_SAMPLES = 2


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def fraction_text(value: Fraction) -> str:
    return str(value)


def task_id(task: dict) -> str:
    return (
        f"t{task['t_low']}..{task['t_high']}:"
        f"x{task['x_low']}..{task['x_high']}"
    )


def initial_tasks() -> list[dict]:
    return [
        {
            "t_low": fraction_text(t_low),
            "t_high": fraction_text(t_high),
            "x_low": fraction_text(x_low),
            "x_high": fraction_text(x_high),
        }
        for t_low, t_high in compact.fraction_range(
            TIME_LOWER, TIME_UPPER, INITIAL_TIME_STEP
        )
        for x_low, x_high in compact.fraction_range(
            X_LOWER, X_UPPER, INITIAL_X_STEP
        )
    ]


def load_full_budget() -> tuple[compact.arb, compact.arb, dict]:
    source = json.loads(FULL_BUDGET_SOURCE.read_text(encoding="utf-8"))
    rows = {
        int(row["N"]): row for row in source.get("budgets", [])
    }
    if RETAINED_N not in rows:
        raise RuntimeError(f"full budget missing N={RETAINED_N}")
    row = rows[RETAINED_N]
    d0 = compact.arb(row["full_d0_m9_t1_5"]["upper"])
    d1 = compact.arb(row["full_d1_m9_t1_5"]["upper"])
    if d0.lower() <= 0 or d1.lower() <= 0:
        raise RuntimeError("full derivative budgets are not positive")
    return d0, d1, row


def cutoff_tail_bounds() -> list[compact.arb]:
    compact.flint.ctx.prec = 256
    polynomials = outer.polynomial_rows()
    switches = outer.switch_derivative_rows()
    bounds: list[compact.arb] = []
    for moment in range(10):
        forward, _ = outer.forward_outer_bound(
            0, moment, 0, 0, polynomials
        )
        defect, _ = outer.defect_outer_bound(
            RETAINED_N,
            moment,
            0,
            0,
            polynomials,
            switches,
        )
        bound = forward + defect
        if bound.lower() <= 0:
            raise RuntimeError(
                f"cutoff tail bound is not positive at moment {moment}"
            )
        bounds.append(compact.arb(bound.upper().str(80)))
    compact.flint.ctx.prec = PRECISION_BITS
    return [compact.arb(value.str(75, more=True)) for value in bounds]


def integration_options() -> dict:
    options = dict(INTEGRATION_OPTIONS)
    for key in ("abs_tol", "rel_tol"):
        value = options[key]
        options[key] = compact.arb(2) ** (-int(value[3:]))
    return options


class ModularRetainedCertifier(compact.CompactCertifier):
    def __init__(self) -> None:
        compact.flint.ctx.prec = PRECISION_BITS
        self.pi = compact.acb.pi()
        self.zero = compact.acb(0)
        self.cutoff = compact.acb(
            compact.fraction_decimal(INTEGRATION_CUTOFF)
        )
        self.options = integration_options()
        self.cutoff_tails = cutoff_tail_bounds()
        self.moment_cache: dict[tuple[Fraction, int], compact.arb] = {}
        self.evaluations = 0
        self.d0, self.d1, self.budget_row = load_full_budget()

    def phi(self, u: compact.acb, index: int) -> compact.acb:
        n2 = index * index
        x = self.pi * n2 * (4 * u).exp()
        return (
            self.pi
            * n2
            * (5 * u).exp()
            * (2 * x - 3)
            * (-x).exp()
        )

    def retained_kernel(self, u: compact.acb) -> compact.acb:
        forward = sum(
            (self.phi(u, index) for index in range(1, RETAINED_N + 1)),
            compact.acb(0),
        )
        reflected = sum(
            (
                self.phi(-u, index)
                for index in range(1, RETAINED_N + 1)
            ),
            compact.acb(0),
        )
        switch_tail = (
            3 * (4 * u).sinh()
        ).erfc() / 2
        return forward - switch_tail * (forward - reflected)

    def tail_complex(self, moment: int) -> compact.acb:
        radius = self.cutoff_tails[moment].upper()
        ball = compact.arb(0, radius.str(80))
        return compact.acb(ball, ball)

    def oscillatory_integral(
        self, moment: int, time_value: Fraction, frequency: Fraction
    ) -> compact.acb:
        time_ball = compact.acb(compact.fraction_decimal(time_value))
        frequency_ball = compact.acb(
            compact.fraction_decimal(frequency)
        )

        def integrand(
            u: compact.acb, _analytic: bool
        ) -> compact.acb:
            return (
                u**moment
                * (time_ball * u * u).exp()
                * self.retained_kernel(u)
                * (compact.acb(0, 1) * frequency_ball * u).exp()
            )

        retained = compact.acb.integral(
            integrand,
            self.zero,
            self.cutoff,
            **self.options,
        )
        return retained + self.tail_complex(moment)

    def moment_upper(
        self, moment: int, time_value: Fraction
    ) -> compact.arb:
        key = (time_value, moment)
        cached = self.moment_cache.get(key)
        if cached is not None:
            return cached
        time_ball = compact.acb(compact.fraction_decimal(time_value))

        def integrand(
            u: compact.acb, _analytic: bool
        ) -> compact.acb:
            return (
                u**moment
                * (time_ball * u * u).exp()
                * self.retained_kernel(u)
            )

        retained = compact.acb.integral(
            integrand,
            self.zero,
            self.cutoff,
            **self.options,
        ).real
        value = retained + self.cutoff_tails[moment].upper()
        if value.lower() <= 0:
            raise RuntimeError(
                f"retained moment is not positive at moment {moment}"
            )
        self.moment_cache[key] = value
        return value

    def certify_modular_box(
        self,
        time_low: Fraction,
        time_high: Fraction,
        x_low: Fraction,
        x_high: Fraction,
        depth: int,
    ) -> dict:
        self.evaluations += 1
        time_center = (time_low + time_high) / 2
        time_radius = (time_high - time_low) / 2
        x_center = (x_low + x_high) / 2
        x_radius = (x_high - x_low) / 2
        max_moment = (
            1
            + 2 * compact.TAYLOR_TIME_ORDER
            + compact.TAYLOR_X_ORDER
        )
        integrals = [
            self.oscillatory_integral(
                moment, time_center, x_center
            )
            for moment in range(max_moment + 1)
        ]
        h_box = self.transform_box(
            0,
            time_center,
            time_radius,
            x_center,
            x_radius,
            integrals,
        )
        h_prime_box = self.transform_box(
            1,
            time_center,
            time_radius,
            x_center,
            x_radius,
            integrals,
        )
        x_box = compact.interval_ball(x_center, x_radius)
        x_three = x_box**3
        x_four = x_box**4
        j_box = 16 * x_four * h_box
        j_prime_box = (
            64 * x_three * h_box + 16 * x_four * h_prime_box
        )

        x_floor = compact.arb(compact.fraction_decimal(x_low))
        tail_value = 16 * self.d0 / x_floor**5
        tail_derivative = (
            64 * self.d0 / x_floor**6
            + 16 * self.d1 / x_floor**5
        )
        tail_value_upper = tail_value.upper()
        tail_derivative_upper = tail_derivative.upper()
        value_lower = j_box.abs_lower()
        derivative_lower = j_prime_box.abs_lower()
        value_ratio = value_lower / tail_value_upper
        derivative_ratio = derivative_lower / tail_derivative_upper
        value_pass = value_lower > tail_value_upper
        derivative_pass = derivative_lower > tail_derivative_upper
        record = {
            "depth": depth,
            "t_low": str(time_low),
            "t_high": str(time_high),
            "x_low": str(x_low),
            "x_high": str(x_high),
            "j_retained_lower": j_box.lower().str(ENDPOINT_DIGITS),
            "j_retained_upper": j_box.upper().str(ENDPOINT_DIGITS),
            "j_retained_prime_lower": j_prime_box.lower().str(
                ENDPOINT_DIGITS
            ),
            "j_retained_prime_upper": j_prime_box.upper().str(
                ENDPOINT_DIGITS
            ),
            "tail_value_upper": tail_value_upper.str(ENDPOINT_DIGITS),
            "tail_derivative_upper": tail_derivative_upper.str(
                ENDPOINT_DIGITS
            ),
            "value_ratio_lower": value_ratio.lower().str(
                ENDPOINT_DIGITS
            ),
            "derivative_ratio_lower": derivative_ratio.lower().str(
                ENDPOINT_DIGITS
            ),
        }
        if value_pass or derivative_pass:
            if value_ratio > derivative_ratio:
                branch = "value"
                ratio = value_ratio
            else:
                branch = "derivative"
                ratio = derivative_ratio
            return {
                **record,
                "certified": True,
                "branch": branch,
                "certified_ratio_lower": ratio.lower().str(
                    ENDPOINT_DIGITS
                ),
                "full_value_negative": bool(
                    j_box.upper() < -tail_value_upper
                ),
                "full_value_positive": bool(
                    j_box.lower() > tail_value_upper
                ),
            }
        return {
            **record,
            "certified": False,
            "branch": None,
            "certified_ratio_lower": max(
                value_ratio, derivative_ratio
            ).lower().str(ENDPOINT_DIGITS),
            "full_value_negative": False,
            "full_value_positive": False,
        }


def split_box(
    bounds: tuple[Fraction, Fraction, Fraction, Fraction],
) -> list[tuple[Fraction, Fraction, Fraction, Fraction]]:
    time_low, time_high, x_low, x_high = bounds
    normalized_time = (
        (time_high - time_low) / INITIAL_TIME_STEP
    )
    normalized_x = (x_high - x_low) / INITIAL_X_STEP
    if normalized_x >= normalized_time:
        middle = (x_low + x_high) / 2
        return [
            (time_low, time_high, x_low, middle),
            (time_low, time_high, middle, x_high),
        ]
    middle = (time_low + time_high) / 2
    return [
        (time_low, middle, x_low, x_high),
        (middle, time_high, x_low, x_high),
    ]


def certify_initial_task(
    certifier: ModularRetainedCertifier, task: dict
) -> dict:
    initial = (
        Fraction(task["t_low"]),
        Fraction(task["t_high"]),
        Fraction(task["x_low"]),
        Fraction(task["x_high"]),
    )
    queue = [(initial, 0)]
    evaluations: list[dict] = []
    certified_leaves: list[dict] = []
    unresolved: list[dict] = []
    subdivisions = 0
    cursor = 0
    started = time.monotonic()
    while cursor < len(queue):
        bounds, depth = queue[cursor]
        cursor += 1
        record = certifier.certify_modular_box(*bounds, depth)
        evaluations.append(record)
        if record["certified"]:
            certified_leaves.append(record)
        elif depth < MAX_DEPTH:
            queue.extend(
                (child, depth + 1) for child in split_box(bounds)
            )
            subdivisions += 1
        else:
            unresolved.append(record)
    ratios = [
        compact.arb(row["certified_ratio_lower"])
        for row in certified_leaves
    ]
    minimum = (
        min(ratios, key=lambda value: float(value.lower()))
        if ratios
        else None
    )
    return {
        "status": "certified" if not unresolved else "unresolved",
        "initial_bounds": task,
        "evaluated_boxes": len(evaluations),
        "certified_leaf_boxes": len(certified_leaves),
        "subdivisions": subdivisions,
        "unresolved_boxes": len(unresolved),
        "maximum_depth": max(
            (row["depth"] for row in evaluations), default=0
        ),
        "minimum_certified_ratio_lower": (
            minimum.lower().str(ENDPOINT_DIGITS)
            if minimum is not None
            else None
        ),
        "branch_counts": {
            "value": sum(
                row.get("branch") == "value"
                for row in certified_leaves
            ),
            "derivative": sum(
                row.get("branch") == "derivative"
                for row in certified_leaves
            ),
        },
        "full_value_negative_leaves": sum(
            row["full_value_negative"] for row in certified_leaves
        ),
        "full_value_positive_leaves": sum(
            row["full_value_positive"] for row in certified_leaves
        ),
        "elapsed_seconds": round(time.monotonic() - started, 6),
        "evaluations": evaluations,
    }


def config_payload() -> dict:
    return {
        "schema": SCHEMA,
        "implementation_sha256": file_hash(Path(__file__)),
        "precision_bits": PRECISION_BITS,
        "retained_n": RETAINED_N,
        "diagonal_j": DIAGONAL_J,
        "domain": {
            "t": [str(TIME_LOWER), str(TIME_UPPER)],
            "x": [str(X_LOWER), str(X_UPPER)],
        },
        "initial_time_step": str(INITIAL_TIME_STEP),
        "initial_x_step": str(INITIAL_X_STEP),
        "max_depth": MAX_DEPTH,
        "integration_cutoff": str(INTEGRATION_CUTOFF),
        "integration_options": INTEGRATION_OPTIONS,
        "taylor_time_order": compact.TAYLOR_TIME_ORDER,
        "taylor_x_order": compact.TAYLOR_X_ORDER,
        "full_budget_sha256": file_hash(FULL_BUDGET_SOURCE),
        "modular_source_sha256": file_hash(MODULAR_SOURCE),
        "outer_source_sha256": file_hash(OUTER_SOURCE),
        "task_ids": [task_id(task) for task in initial_tasks()],
    }


def config_hash() -> str:
    return sha256(
        canonical_json(config_payload()).encode("utf-8")
    ).hexdigest()


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
                    f"cache line {line_number} lacks a terminator"
                )
            record = json.loads(line)
            stored_hash = record.pop("row_sha256", None)
            actual_hash = sha256(
                canonical_json(record).encode("utf-8")
            ).hexdigest()
            record["row_sha256"] = stored_hash
            if stored_hash != actual_hash:
                raise RuntimeError(
                    f"cache hash mismatch at line {line_number}"
                )
            if record.get("sequence") != line_number:
                raise RuntimeError(
                    f"cache sequence mismatch at line {line_number}"
                )
            if record.get("previous_row_sha256") != previous:
                raise RuntimeError(
                    f"cache chain mismatch at line {line_number}"
                )
            if record.get("config_sha256") != expected_config:
                raise RuntimeError(
                    f"cache config mismatch at line {line_number}"
                )
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
    result: dict,
) -> dict:
    record = {
        "schema": SCHEMA,
        "sequence": len(records) + 1,
        "config_sha256": config_hash(),
        "task": task,
        "result": result,
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


def run_tasks(
    cache: Path,
    max_new_tasks: int | None,
    max_seconds: float | None,
    stop_file: Path | None,
) -> tuple[list[dict], dict[str, dict], int, list[float], bool]:
    records, by_id = load_cache(cache)
    pending = [
        task
        for task in initial_tasks()
        if task_id(task) not in by_id
    ]
    certifier: ModularRetainedCertifier | None = None
    completed = 0
    cpu_samples: list[float] = []
    consecutive_high = 0
    resource_parked = False
    started = time.monotonic()
    for task in pending:
        if stop_file is not None and stop_file.exists():
            break
        if max_new_tasks is not None and completed >= max_new_tasks:
            break
        if (
            max_seconds is not None
            and completed > 0
            and time.monotonic() - started >= max_seconds
        ):
            break
        if certifier is None:
            certifier = ModularRetainedCertifier()
        task_started = time.monotonic()
        result = certify_initial_task(certifier, task)
        record = append_record(cache, records, task, result)
        by_id[task_id(task)] = record
        completed += 1
        print(
            f"completed {record['sequence']}/{len(initial_tasks())} "
            f"{task_id(task)} status={result['status']} "
            f"eval={result['evaluated_boxes']} "
            f"in {time.monotonic() - task_started:.3f}s",
            flush=True,
        )
        sample = float(psutil.cpu_percent(interval=0.2))
        cpu_samples.append(sample)
        if sample > CPU_PARK_THRESHOLD:
            consecutive_high += 1
        else:
            consecutive_high = 0
        if consecutive_high >= CPU_CONSECUTIVE_PARK_SAMPLES:
            resource_parked = True
            break
    return records, by_id, completed, cpu_samples, resource_parked


def source_audit() -> dict:
    paths = {
        "full_budget": FULL_BUDGET_SOURCE,
        "modular_partition": MODULAR_SOURCE,
        "outer_tail": OUTER_SOURCE,
        "compact_core": COMPACT_SOURCE,
        "diagonal_exhaustion": DIAGONAL_SOURCE,
        "winding": WINDING_SOURCE,
        "boundary": BOUNDARY_SOURCE,
    }
    sources = {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in paths.items()
    }
    text = {key: json.dumps(value) for key, value in sources.items()}
    markers = {
        "full_budget": ("full_d0_m9_t1_5", "N=4..10"),
        "modular_partition": (
            "sum_(n>=1)b_n(u)=Phi(u)",
            "b_n(u)>0",
        ),
        "outer_tail": ("u>11/5",),
        "compact_core": ("|x|<=38",),
        "diagonal_exhaustion": (
            "R_j=38+j",
            "delta_j=1/(5j)",
        ),
        "winding": (
            "Q_j=[1/(5j),1/4]x[0,38+j]",
            "wind(Z(partial Q_j),0)=0",
        ),
        "boundary": ("Lambda<=1/5",),
    }
    for key, required in markers.items():
        for marker in required:
            if marker not in text[key]:
                raise RuntimeError(
                    f"{key} source marker missing: {marker}"
                )
    return {
        key: {
            "kind": sources[key]["kind"],
            "path": str(paths[key].relative_to(REPO_ROOT)).replace(
                "\\", "/"
            ),
            "sha256": file_hash(paths[key]),
        }
        for key in paths
    }


def build_artifact(
    cache: Path,
    records: list[dict],
    by_id: dict[str, dict],
    cpu_samples: list[float],
    resource_parked: bool,
) -> dict:
    tasks = initial_tasks()
    missing = [
        task_id(task) for task in tasks if task_id(task) not in by_id
    ]
    results = [record["result"] for record in records]
    unresolved = [
        result for result in results if result["status"] != "certified"
    ]
    complete = not missing
    theorem_ready = complete and not unresolved
    certified_leaves = [
        evaluation
        for result in results
        for evaluation in result["evaluations"]
        if evaluation["certified"]
    ]
    ratios = [
        compact.arb(row["certified_ratio_lower"])
        for row in certified_leaves
    ]
    minimum = (
        min(ratios, key=lambda value: float(value.lower()))
        if ratios
        else None
    )
    rows = [
        GateRow(
            id="ntmrq31_01_even_modular_retained_sum",
            role="exact_input",
            readiness="proved",
            claim=(
                "The retained transform uses the first seven positive even "
                "modular blocks."
            ),
            formula="S_(7,t)=sum_(n=1)^7 B_(n,t)",
            proof_boundary="No tail is discarded without a certified bar.",
        ),
        GateRow(
            id="ntmrq31_02_full_m9_tail_bars",
            role="rigorous_input",
            readiness="proved",
            claim=(
                "The full N=7 m=9 derivative budgets bound J and J' errors "
                "uniformly on the slab."
            ),
            formula=(
                "|J-J_7|<=16*d0_7/x^5; "
                "|J'-J_7'|<=64*d0_7/x^6+16*d1_7/x^5"
            ),
            proof_boundary="Uses the independently checked full budget.",
        ),
        GateRow(
            id="ntmrq31_03_retained_cutoff_tail",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "Every retained moment through order nine has an explicit "
                "analytic tail beyond u=11/5."
            ),
            formula=(
                "tail_m<=forward_all_m+sum_(n<=7)defect_(n,m)"
            ),
            proof_boundary="No sampled improper-integral cutoff is used.",
        ),
        GateRow(
            id="ntmrq31_04_interval_slab",
            role="directed_rounding_certificate",
            readiness=("proved" if theorem_ready else "not_ready_to_apply"),
            claim=(
                "Every Taylor leaf in the Q_31 transition slab satisfies "
                "the direct retained value-or-derivative disjunction."
            ),
            formula=(
                "[1/155,1/5]x[38,69] is covered by certified leaves"
            ),
            proof_boundary=(
                "Claim activates only when the cache is complete with zero "
                "unresolved initial boxes."
            ),
        ),
        GateRow(
            id="ntmrq31_05_q31_no_contact",
            role="finite_diagonal_theorem",
            readiness=("proved" if theorem_ready else "not_ready_to_apply"),
            claim=(
                "The compact core, certified slab, and t>1/5 boundary "
                "theorem exclude every first-jet contact in Q_31."
            ),
            formula="Q_31=[1/155,1/4]x[0,69]",
            proof_boundary="A finite diagonal theorem is not cofinal.",
        ),
        GateRow(
            id="ntmrq31_06_q31_winding",
            role="finite_diagonal_theorem",
            readiness=("proved" if theorem_ready else "not_ready_to_apply"),
            claim=(
                "Boundary nonvanishing and the contact-index theorem give "
                "zero first-jet winding on Q_31."
            ),
            formula="wind((H+iH_x)(partial Q_31),0)=0",
            proof_boundary="A finite winding theorem is not Lambda<=0.",
        ),
        GateRow(
            id="ntmrq31_07_cofinal_handoff",
            role="proof_search_target",
            readiness="not_ready_to_apply",
            claim=(
                "Replace the finite retained count and finite rectangle by "
                "a cofinal adaptive-N transition theorem."
            ),
            formula=(
                "N=N(x,t) unbounded; direct C1 disjunction on every "
                "x>38, 0<t<=1/5"
            ),
            proof_boundary=(
                "This remains the missing Newman direction."
            ),
        ),
    ]
    return {
        "kind": STEM,
        "date": DATE,
        "status": (
            "complete rigorous Q_31 modular-retained interval theorem"
            if theorem_ready
            else "partial resumable Q_31 modular-retained interval certificate"
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
            "resource_parked_last_run": resource_parked,
            "cpu_samples_last_run": cpu_samples,
        },
        "summary": {
            "theorem_ready": theorem_ready,
            "certified_initial_boxes": sum(
                result["status"] == "certified" for result in results
            ),
            "unresolved_initial_boxes": len(unresolved),
            "evaluated_boxes": sum(
                result["evaluated_boxes"] for result in results
            ),
            "certified_leaf_boxes": len(certified_leaves),
            "subdivisions": sum(
                result["subdivisions"] for result in results
            ),
            "branch_counts": {
                "value": sum(
                    result["branch_counts"]["value"] for result in results
                ),
                "derivative": sum(
                    result["branch_counts"]["derivative"]
                    for result in results
                ),
            },
            "minimum_certified_ratio_lower": (
                minimum.lower().str(ENDPOINT_DIGITS)
                if minimum is not None
                else None
            ),
        },
        "source_audit": source_audit(),
        "rows": [asdict(row) for row in rows],
        "proof_boundary": (
            "Even when complete, this artifact proves only the finite "
            "Q_31 rectangle. It does not supply an unbounded adaptive-N "
            "transition theorem, prove Lambda<=0 or RH, or provide a "
            "Clay-prize conclusion."
        ),
    }


def render_note(artifact: dict) -> str:
    cache = artifact["cache"]
    summary = artifact["summary"]
    lines = [
        "# Newman Theta Modular-Retained Q31 Interval Certificate",
        "",
        f"Date: {DATE}",
        "",
        "Status: "
        + (
            "complete rigorous finite Q31 theorem."
            if summary["theorem_ready"]
            else "partial resumable interval certificate."
        ),
        "This is not a cofinal theorem and not a proof of `Lambda<=0`,",
        "RH, or a Clay-prize result.",
        "",
        "## Domain",
        "",
        "```text",
        "retained modular blocks: N=7",
        "transition slab: [1/155,1/5] x [38,69]",
        "diagonal rectangle: Q_31=[1/155,1/4] x [0,69]",
        "```",
        "",
        "## Progress",
        "",
        f"Initial boxes: `{cache['records']}/{cache['expected_records']}`.",
        f"Certified initial boxes: `{summary['certified_initial_boxes']}`.",
        f"Unresolved initial boxes: `{summary['unresolved_initial_boxes']}`.",
        f"Evaluated Taylor boxes: `{summary['evaluated_boxes']}`.",
        f"Subdivisions: `{summary['subdivisions']}`.",
        "Minimum certified ratio: "
        f"`{summary['minimum_certified_ratio_lower']}`.",
        "",
        "## Proof Boundary",
        "",
        "A complete cache proves only this finite rectangle. The surviving",
        "obligation is a cofinal adaptive-N transition theorem for every",
        "`x>38` and `0<t<=1/5`.",
        "",
    ]
    return "\n".join(lines)


def write_artifact(artifact: dict, out: Path, note: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    note.write_text(render_note(artifact), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--max-new-tasks", type=int)
    parser.add_argument("--max-seconds", type=float)
    parser.add_argument("--stop-file", type=Path)
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()

    priority_lowered = compact.request_below_normal_priority()
    records, by_id, completed, cpu_samples, resource_parked = run_tasks(
        args.cache,
        args.max_new_tasks,
        args.max_seconds,
        args.stop_file,
    )
    artifact = build_artifact(
        args.cache,
        records,
        by_id,
        cpu_samples,
        resource_parked,
    )
    artifact["below_normal_priority_applied"] = priority_lowered
    write_artifact(artifact, args.out, args.note)
    print(
        f"Q31 progress: new={completed} "
        f"cache={artifact['cache']['records']}/"
        f"{artifact['cache']['expected_records']} "
        f"unresolved={artifact['summary']['unresolved_initial_boxes']} "
        f"resource_parked={resource_parked}",
        flush=True,
    )
    if args.require_complete and not artifact["summary"]["theorem_ready"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
