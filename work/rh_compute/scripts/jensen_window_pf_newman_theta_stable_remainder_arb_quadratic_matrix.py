#!/usr/bin/env python3
"""Build the rigorous compact matrix from the stable modular remainder."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_theta_modular_tail_arb_quadratic_matrix as matrix


matrix.STEM = (
    "jensen_window_pf_newman_theta_stable_remainder_arb_quadratic_matrix"
)
matrix.SCHEMA = "newman_theta_stable_remainder_arb_quadratic_matrix_v1"
matrix.PRECISION_BITS = 192
matrix.N_VALUES = tuple(range(4, 11))
matrix.N_CAP = 12
matrix.INTEGRATION_OPTIONS = {
    "deg_limit": 100,
    "eval_limit": 200_000,
    "depth_limit": 35,
    "use_heap": True,
    "abs_tol": "2^-320",
    "rel_tol": "2^-120",
}
matrix.DEFAULT_CACHE = matrix.RESULT_DIR / f"{matrix.STEM}.jsonl"
matrix.DEFAULT_OUT = matrix.RESULT_DIR / f"{matrix.STEM}.json"
matrix.DEFAULT_NOTE = matrix.REPO_ROOT / "outputs" / f"{matrix.STEM}.md"

TAIL_GATE = (
    matrix.RESULT_DIR
    / "jensen_window_pf_newman_theta_forward_remainder_tail_gate.json"
)
MODULAR_SOURCE = (
    matrix.RESULT_DIR
    / "jensen_window_pf_newman_theta_modular_blend_gate.json"
)
FORWARD_CAP = 12
FIRST_OMITTED = 13


def config_payload() -> dict:
    return {
        "schema": matrix.SCHEMA,
        "precision_bits": matrix.PRECISION_BITS,
        "t_cap_exact": "1/5",
        "u_interval_exact": ["0", "11/5"],
        "forward_cap": FORWARD_CAP,
        "first_omitted": FIRST_OMITTED,
        "n_values": list(matrix.N_VALUES),
        "m_order": matrix.M_ORDER,
        "integration_options": matrix.INTEGRATION_OPTIONS,
        "tail_gate_sha256": matrix.file_hash(TAIL_GATE),
        "modular_source_sha256": matrix.file_hash(MODULAR_SOURCE),
        "task_ids": [
            matrix.task_id(task) for task in matrix.expected_tasks()
        ],
    }


def build_derivative_functions() -> dict[str, list]:
    u, n = matrix.sp.symbols("u n", positive=True, real=True)

    def phi(z: matrix.sp.Expr) -> matrix.sp.Expr:
        return (
            2 * matrix.sp.pi**2 * n**4 * matrix.sp.exp(9 * z)
            - 3 * matrix.sp.pi * n**2 * matrix.sp.exp(5 * z)
        ) * matrix.sp.exp(
            -matrix.sp.pi * n**2 * matrix.sp.exp(4 * z)
        )

    switch_tail = matrix.sp.erfc(3 * matrix.sp.sinh(4 * u)) / 2
    forward = phi(u)
    defect = switch_tail * (phi(u) - phi(-u))
    forward_expressions: list[matrix.sp.Expr] = []
    defect_expressions: list[matrix.sp.Expr] = []
    for _ in range(matrix.M_ORDER + 1):
        forward_expressions.append(forward)
        defect_expressions.append(defect)
        forward = matrix.sp.diff(forward, u)
        defect = matrix.sp.diff(defect, u)

    functions = {
        "exp": lambda z: z.exp(),
        "sinh": lambda z: z.sinh(),
        "cosh": lambda z: z.cosh(),
        "erf": lambda z: z.erf(),
        "erfc": lambda z: z.erfc(),
        "sqrt": lambda z: z.sqrt(),
        "pi": matrix.acb.pi(),
    }

    def compile_all(expressions: list[matrix.sp.Expr]) -> list:
        return [
            matrix.sp.lambdify(
                (u, n),
                expression,
                modules=[functions],
                cse=True,
                docstring_limit=0,
            )
            for expression in expressions
        ]

    return {
        "forward": compile_all(forward_expressions),
        "defect": compile_all(defect_expressions),
    }


def integrate_task(task: dict, derivatives: dict[str, list] | None):
    if task["kind"] == "mass":
        integrand = lambda z, _analytic: matrix.weight(
            z, task["p"], task["b"]
        )
    else:
        if derivatives is None:
            raise RuntimeError("quadratic task requires derivative functions")
        q = task["q"]
        retained = task["N"]

        def integrand(z: matrix.acb, _analytic: bool) -> matrix.acb:
            finite_remainder = sum(
                (
                    derivatives["forward"][q](z, index)
                    for index in range(retained + 1, FORWARD_CAP + 1)
                ),
                matrix.acb(0),
            )
            finite_remainder += sum(
                (
                    derivatives["defect"][q](z, index)
                    for index in range(1, retained + 1)
                ),
                matrix.acb(0),
            )
            return (
                matrix.weight(z, task["p"], task["b"])
                * finite_remainder
                * finite_remainder
            )

    options = dict(matrix.INTEGRATION_OPTIONS)
    for key in ("abs_tol", "rel_tol"):
        value = options.get(key)
        if isinstance(value, str) and value.startswith("2^-"):
            options[key] = matrix.arb(2) ** (-int(value[3:]))

    result = matrix.acb.integral(
        integrand,
        matrix.arb(0),
        matrix.arb(matrix.U_STOP),
        **options,
    )
    if not result.is_finite():
        raise RuntimeError(
            f"non-finite integral for {matrix.task_id(task)}"
        )
    if task["kind"] == "mass" and result.real.lower() <= 0:
        raise RuntimeError(
            "weight mass was not separated for "
            f"{matrix.task_id(task)}: {result}"
        )
    if task["kind"] == "quadratic" and result.real.upper() <= 0:
        raise RuntimeError(
            "quadratic upper bound was not positive for "
            f"{matrix.task_id(task)}: {result}"
        )
    if not result.imag.contains(0):
        raise RuntimeError(
            "imaginary enclosure excludes zero for "
            f"{matrix.task_id(task)}"
        )
    return result


def load_tail_constants() -> dict[int, matrix.arb]:
    source = json.loads(TAIL_GATE.read_text(encoding="utf-8"))
    constants: dict[int, matrix.arb] = {}
    pi = matrix.arb.pi()
    for row in source["tail_bounds"]:
        q = int(row["q"])
        d = 2 * q + 4
        rho = (
            matrix.arb(d) / FIRST_OMITTED
            - pi * (2 * FIRST_OMITTED + 1)
        ).exp()
        bound = (
            matrix.arb(int(source["polynomials"][q]["A_q"]))
            * pi ** (q + 2)
            * FIRST_OMITTED**d
            * (-pi * FIRST_OMITTED**2).exp()
            / (1 - rho)
        )
        if not bound.is_finite() or bound.lower() <= 0:
            raise RuntimeError(f"invalid forward-tail bound at q={q}")
        constants[q] = matrix.arb(bound.upper())
    if sorted(constants) != list(range(matrix.M_ORDER + 1)):
        raise RuntimeError("forward-tail constants are incomplete")
    return constants


def serialize_ball(value: matrix.arb) -> dict:
    return {
        "upper_accumulator_enclosure": value.str(40, more=True),
        "upper": value.upper().str(40),
        "relative_accuracy_bits": value.rel_accuracy_bits(),
    }


def aggregate_n(retained: int, by_id: dict[str, dict]) -> dict | None:
    required = [
        matrix.task_id(
            {
                "kind": "quadratic",
                "N": retained,
                "p": row["p"],
                "b": row["b"],
                "q": row["q"],
            }
        )
        for row in matrix.matrix_specs()
    ]
    if any(identifier not in by_id for identifier in required):
        return None

    tail_constants = load_tail_constants()
    finite = {
        "d0": matrix.arb(0),
        "d1_u": matrix.arb(0),
        "d1_lower": matrix.arb(0),
    }
    arithmetic_tail = {
        "d0": matrix.arb(0),
        "d1_u": matrix.arb(0),
        "d1_lower": matrix.arb(0),
    }
    for spec in matrix.matrix_specs():
        mass_id = matrix.task_id(
            {"kind": "mass", "p": spec["p"], "b": spec["b"]}
        )
        quadratic_id = matrix.task_id(
            {
                "kind": "quadratic",
                "N": retained,
                "p": spec["p"],
                "b": spec["b"],
                "q": spec["q"],
            }
        )
        mass = matrix.parse_positive_upper(by_id[mass_id])
        quadratic = matrix.parse_positive_upper(by_id[quadratic_id])
        coefficient = spec["coefficient"]
        finite[spec["component"]] += coefficient * (
            mass * quadratic
        ).sqrt()
        arithmetic_tail[spec["component"]] += (
            coefficient * tail_constants[spec["q"]] * mass
        )

    total = {
        key: finite[key] + arithmetic_tail[key] for key in finite
    }
    finite_d1 = finite["d1_u"] + finite["d1_lower"]
    arithmetic_d1 = (
        arithmetic_tail["d1_u"] + arithmetic_tail["d1_lower"]
    )
    total_d1 = total["d1_u"] + total["d1_lower"]
    return {
        "N": retained,
        "finite_stable_d0": serialize_ball(finite["d0"]),
        "finite_stable_d1_u": serialize_ball(finite["d1_u"]),
        "finite_stable_d1_lower": serialize_ball(finite["d1_lower"]),
        "finite_stable_d1": serialize_ball(finite_d1),
        "forward_arithmetic_tail_d0": serialize_ball(
            arithmetic_tail["d0"]
        ),
        "forward_arithmetic_tail_d1_u": serialize_ball(
            arithmetic_tail["d1_u"]
        ),
        "forward_arithmetic_tail_d1_lower": serialize_ball(
            arithmetic_tail["d1_lower"]
        ),
        "forward_arithmetic_tail_d1": serialize_ball(arithmetic_d1),
        "compact_full_d0": serialize_ball(total["d0"]),
        "compact_full_d1_u": serialize_ball(total["d1_u"]),
        "compact_full_d1_lower": serialize_ball(total["d1_lower"]),
        "compact_full_d1": serialize_ball(total_d1),
        "quadratic_entries": len(matrix.matrix_specs()),
    }


def build_artifact(
    cache: Path,
    records: list[dict],
    by_id: dict[str, dict],
) -> dict:
    tasks = matrix.expected_tasks()
    complete_rows = [
        row
        for retained in matrix.N_VALUES
        if (row := aggregate_n(retained, by_id)) is not None
    ]
    missing = [
        matrix.task_id(task)
        for task in tasks
        if matrix.task_id(task) not in by_id
    ]
    complete = not missing
    return {
        "kind": matrix.STEM,
        "date": matrix.DATE,
        "status": (
            "complete full-arithmetic compact Arb matrix"
            if complete
            else "partial resumable stable-remainder compact Arb matrix"
        ),
        "config": config_payload(),
        "config_sha256": matrix.config_hash(),
        "cache": {
            "path": str(cache.relative_to(matrix.REPO_ROOT)).replace(
                "\\", "/"
            ),
            "records": len(records),
            "expected_records": len(tasks),
            "last_row_sha256": (
                records[-1]["row_sha256"] if records else None
            ),
            "complete": complete,
            "missing_count": len(missing),
            "missing_task_ids": missing,
        },
        "stable_remainder_definition": {
            "switch_tail": (
                "s(u)=erfc(3*sinh(4u))/2=1-omega(u)"
            ),
            "retained_defect": (
                "delta_n(u)=s(u)*(phi_n(u)-phi_n(-u))=phi_n(u)-b_n(u)"
            ),
            "finite_part": (
                "f_N(u)=sum_(n=N+1)^12 phi_n(u)"
                "+sum_(n=1)^N delta_n(u)"
            ),
            "forward_tail": (
                "tau(u)=sum_(n>=13)phi_n(u)"
            ),
            "identity": "r_N(u)=f_N(u)+tau(u)",
        },
        "matrix_definition": {
            "order": matrix.M_ORDER,
            "entry_count_per_N": len(matrix.matrix_specs()),
            "mass_count": len(
                {
                    (row["p"], row["b"])
                    for row in matrix.matrix_specs()
                }
            ),
            "mass": (
                "M_(p,b)=integral_0^(11/5) "
                "u^p exp(u^2/5)P_b(1/5,u)du"
            ),
            "quadratic": (
                "Q_(N,p,b,q)=integral_0^(11/5) "
                "u^p exp(u^2/5)P_b(1/5,u)[D^q f_N(u)]^2du"
            ),
            "entry_bound": (
                "integral W|D^q r_N|<="
                "sqrt(M_(p,b)Q_(N,p,b,q))+B_q*M_(p,b)"
            ),
            "d0": (
                "sum_(b=0)^9 binom(9,b) entry_(N,0,b,9-b)"
            ),
            "d1_u": (
                "sum_(b=0)^9 binom(9,b) entry_(N,1,b,9-b)"
            ),
            "d1_lower": (
                "9*sum_(b=0)^8 binom(8,b) entry_(N,0,b,8-b)"
            ),
        },
        "compact_rows": complete_rows,
        # Retain the base runner's print-key contract without changing the
        # stable artifact's reader-facing name.
        "finite_compact_rows": complete_rows,
        "proof_boundary": (
            "Every completed row is a directed-rounding certificate for "
            "the full arithmetic remainder n>N on 0<=u<=11/5, including "
            "the explicit n>=13 correction. It omits u>11/5, retained "
            "first-jet J/J' lower balls, transition-cell coverage, "
            "Lambda<=0, RH, and any "
            "Clay-prize conclusion."
        ),
        "sources": {
            "tail_gate": str(
                TAIL_GATE.relative_to(matrix.REPO_ROOT)
            ).replace("\\", "/"),
            "tail_gate_sha256": matrix.file_hash(TAIL_GATE),
            "modular_partition": str(
                MODULAR_SOURCE.relative_to(matrix.REPO_ROOT)
            ).replace("\\", "/"),
            "modular_partition_sha256": matrix.file_hash(MODULAR_SOURCE),
        },
        "versions": {
            "python_flint": getattr(matrix.flint, "__version__", "unknown"),
            "sympy": matrix.sp.__version__,
        },
    }


def render_note(artifact: dict) -> str:
    cache = artifact["cache"]
    rows = artifact["compact_rows"]
    lines = [
        "# Newman Theta Stable-Remainder Arb Quadratic Matrix",
        "",
        f"Date: {matrix.DATE}",
        "",
        "Status: "
        + (
            "complete rigorous full-arithmetic compact matrix."
            if cache["complete"]
            else "partial resumable rigorous compact matrix."
        ),
        "This is not an outer-tail certificate and not a proof of",
        "`Lambda<=0`, RH, or a Clay-prize result.",
        "",
        "## Stable Identity",
        "",
        "```text",
        artifact["stable_remainder_definition"]["switch_tail"],
        artifact["stable_remainder_definition"]["retained_defect"],
        artifact["stable_remainder_definition"]["finite_part"],
        artifact["stable_remainder_definition"]["forward_tail"],
        artifact["stable_remainder_definition"]["identity"],
        "```",
        "",
        "The finite integrand evaluates the small switch defect directly,",
        "rather than subtracting two independently accumulated full sums.",
        "",
        "## Matrix",
        "",
        "```text",
        artifact["matrix_definition"]["mass"],
        artifact["matrix_definition"]["quadratic"],
        artifact["matrix_definition"]["entry_bound"],
        artifact["matrix_definition"]["d0"],
        artifact["matrix_definition"]["d1_u"],
        artifact["matrix_definition"]["d1_lower"],
        "```",
        "",
        "The cache is append-only, hash-chained, and fsynced after every",
        "Arb integral.",
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
                "| N | full d0 upper | full d1 upper |",
                "|---:|---:|---:|",
            ]
        )
        for row in rows:
            lines.append(
                f"| {row['N']} | "
                f"`{row['compact_full_d0']['upper']}` | "
                f"`{row['compact_full_d1']['upper']}` |"
            )
        lines.append("")
    lines.extend(
        [
            "## Non-Promotion Guard",
            "",
            "Completed rows cover the entire arithmetic remainder on",
            "`0<=u<=11/5`. The outer interval `u>11/5`, retained first-jet",
            "lower separation, and transition-cell theorem remain open.",
            "",
        ]
    )
    return "\n".join(lines)


matrix.config_payload = config_payload
matrix.build_derivative_functions = build_derivative_functions
matrix.integrate_task = integrate_task
matrix.aggregate_n = aggregate_n
matrix.build_artifact = build_artifact
matrix.render_note = render_note


if __name__ == "__main__":
    raise SystemExit(matrix.main())
