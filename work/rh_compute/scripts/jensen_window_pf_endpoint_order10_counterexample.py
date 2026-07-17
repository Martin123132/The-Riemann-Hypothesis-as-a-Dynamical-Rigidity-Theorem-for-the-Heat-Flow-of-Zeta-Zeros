#!/usr/bin/env python3
"""Certify the first-open-order endpoint counterexample at signed order ten."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order9_m100_finite_splice_certificate as splice  # noqa: E402
import jensen_window_pf_compound_order9_m100_prefix_certificate as prefix  # noqa: E402


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_endpoint_order10_counterexample.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_endpoint_order10_counterexample.md"
SOURCE_PATHS = splice.SOURCE_PATHS
MAX_COEFFICIENT_INDEX = 1258
PREFIX_LAST_N = 1240
DIRECT_LAST_N = 4
PRECISION_BITS = 4096


@dataclass(frozen=True)
class CounterexampleRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def orientation(order: int) -> int:
    return (-1) ** (order * (order - 1) // 2)


def sign_class(value: flint.arb) -> str:
    if bool(value > 0 and not value.contains(0)):
        return "positive"
    if bool(value < 0 and not value.contains(0)):
        return "negative"
    return "inconclusive"


def ball_text(value: flint.arb, digits: int = 100) -> str:
    return value.str(digits).replace("e", "E")


def source_record(path: Path, precedence: int, rows: int, index_range: list[int]) -> dict:
    return {
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "precedence": precedence,
        "rows": rows,
        "index_range": index_range,
        "sha256": sha256(path),
    }


def load_coefficients() -> tuple[dict[int, flint.arb], list[dict]]:
    values, diagnostics = prefix.merged_coefficients(
        SOURCE_PATHS, MAX_COEFFICIENT_INDEX
    )
    records = []
    for precedence, (path, diagnostic) in enumerate(zip(SOURCE_PATHS, diagnostics)):
        records.append(
            source_record(
                path,
                precedence,
                int(diagnostic["rows"]),
                [int(value) for value in diagnostic["index_range"]],
            )
        )
    return values, records


def direct_signed_hankel(
    values: dict[int, flint.arb], order: int, shift: int
) -> flint.arb:
    matrix = flint.arb_mat(
        [[values[shift + row + col] for col in range(order)] for row in range(order)]
    )
    return orientation(order) * matrix.det()


def direct_deep_schur(
    values: dict[int, flint.arb], order: int, width: int
) -> flint.arb:
    h = {index: value / values[0] for index, value in values.items()}
    matrix = flint.arb_mat(
        [
            [h.get(width - row + col, flint.arb(0)) for col in range(order)]
            for row in range(order)
        ]
    )
    return matrix.det()


def stable_order10_scan(values: dict[int, flint.arb]) -> dict:
    contractions = {
        index: values[index - 1] * values[index + 1] / values[index] ** 2
        for index in range(1, MAX_COEFFICIENT_INDEX)
    }
    defects = {index: 1 - value for index, value in contractions.items()}
    gaps = {
        index: (
            defects[index + 2] ** 2
            - contractions[index + 2] ** 2
            * defects[index + 1]
            * defects[index + 3]
        )
        for index in range(MAX_COEFFICIENT_INDEX - 3)
    }
    order4_margins = {
        index: (
            gaps[index + 1] ** 2
            - contractions[index + 3] ** 3 * gaps[index] * gaps[index + 2]
        )
        for index in range(MAX_COEFFICIENT_INDEX - 5)
    }
    stable_h5 = {
        index: (
            defects[index + 3]
            * defects[index + 5]
            * order4_margins[index + 1] ** 2
            - contractions[index + 4] ** 4
            * defects[index + 4] ** 2
            * order4_margins[index]
            * order4_margins[index + 2]
        )
        for index in range(MAX_COEFFICIENT_INDEX - 7)
    }
    for name, family in (
        ("contraction", contractions),
        ("defect", defects),
        ("order3_gap", gaps),
        ("order4_margin", order4_margins),
        ("order5_stable_margin", stable_h5),
    ):
        if not all(sign_class(value) == "positive" for value in family.values()):
            raise RuntimeError(f"{name} family lost strict positivity")

    def h5_value(index: int) -> flint.arb:
        ratio = values[index + 1] / values[index]
        scale = (
            values[index] ** 5
            * ratio**20
            * contractions[index + 1] ** 15
            * contractions[index + 2] ** 10
            * contractions[index + 3] ** 5
            / (
                defects[index + 3]
                * defects[index + 4] ** 2
                * defects[index + 5]
                * gaps[index + 2]
            )
        )
        return scale * stable_h5[index]

    def h4_value(index: int) -> flint.arb:
        ratio = values[index + 1] / values[index]
        return (
            values[index] ** 4
            * ratio**12
            * contractions[index + 1] ** 8
            * contractions[index + 2] ** 4
            * order4_margins[index]
            / defects[index + 3]
        )

    h5_values = {
        index: h5_value(index)
        for index in range(MAX_COEFFICIENT_INDEX - 7)
    }
    h4_values = {
        index: h4_value(index)
        for index in range(MAX_COEFFICIENT_INDEX - 7)
    }
    q6_values: dict[int, flint.arb] = {}
    for index in range(MAX_COEFFICIENT_INDEX - 9):
        relative = h5_values[index + 1] ** 2 / (
            h5_values[index] * h5_values[index + 2]
        ) - 1
        q6_values[index] = (
            h5_values[index] * h5_values[index + 2] * relative / h4_values[index + 2]
        )
    q7_values: dict[int, flint.arb] = {}
    for index in range(MAX_COEFFICIENT_INDEX - 11):
        relative = q6_values[index + 1] ** 2 / (
            q6_values[index] * q6_values[index + 2]
        ) - 1
        q7_values[index] = (
            q6_values[index] * q6_values[index + 2] * relative / h5_values[index + 2]
        )
    q8_values: dict[int, flint.arb] = {}
    for index in range(MAX_COEFFICIENT_INDEX - 13):
        relative = q7_values[index + 1] ** 2 / (
            q7_values[index] * q7_values[index + 2]
        ) - 1
        q8_values[index] = (
            q7_values[index] * q7_values[index + 2] * relative / q6_values[index + 2]
        )
    q9_values: dict[int, flint.arb] = {}
    for index in range(MAX_COEFFICIENT_INDEX - 15):
        relative = q8_values[index + 1] ** 2 / (
            q8_values[index] * q8_values[index + 2]
        ) - 1
        if sign_class(relative) != "positive":
            raise RuntimeError(f"signed order-nine input failed at n={index}")
        q9_values[index] = (
            q8_values[index] * q8_values[index + 2] * relative / q7_values[index + 2]
        )

    relative_q9 = {
        index: q9_values[index + 1] ** 2
        / (q9_values[index] * q9_values[index + 2])
        - 1
        for index in range(PREFIX_LAST_N + 1)
    }
    classifications = {
        kind: [
            index
            for index, value in relative_q9.items()
            if sign_class(value) == kind
        ]
        for kind in ("positive", "negative", "inconclusive")
    }
    if classifications["negative"] != [0, 1, 2, 3]:
        raise RuntimeError(
            f"unexpected order-ten negative set: {classifications['negative']}"
        )
    if classifications["inconclusive"]:
        raise RuntimeError("order-ten prefix contains inconclusive margins")
    if classifications["positive"] != list(range(4, PREFIX_LAST_N + 1)):
        raise RuntimeError("order-ten positive range is not contiguous from n=4")

    negative_rows = [
        {
            "n": index,
            "N": index + 9,
            "shape": f"(({index + 9}^10))",
            "relative_Q9_margin_ball": ball_text(relative_q9[index]),
            "classification": "Q_(10,n)<0",
        }
        for index in classifications["negative"]
    ]
    selected_positive_rows = [
        {
            "n": index,
            "N": index + 9,
            "shape": f"(({index + 9}^10))",
            "relative_Q9_margin_ball": ball_text(relative_q9[index]),
            "classification": "Q_(10,n)>0",
        }
        for index in (4, PREFIX_LAST_N)
    ]
    return {
        "lambda": "-100",
        "n_range": [0, PREFIX_LAST_N],
        "coefficient_range": [0, MAX_COEFFICIENT_INDEX],
        "precision_bits": PRECISION_BITS,
        "positive_rows": len(classifications["positive"]),
        "negative_rows": len(classifications["negative"]),
        "inconclusive_rows": len(classifications["inconclusive"]),
        "negative_indices": classifications["negative"],
        "positive_range": [4, PREFIX_LAST_N],
        "negative_records": negative_rows,
        "selected_positive_records": selected_positive_rows,
        "relative_margins": relative_q9,
    }


def direct_audit(
    values: dict[int, flint.arb], stable: dict
) -> list[dict]:
    rows = []
    stable_margins = stable["relative_margins"]
    for shift in range(DIRECT_LAST_N + 1):
        q8 = direct_signed_hankel(values, 8, shift + 2)
        q9_left = direct_signed_hankel(values, 9, shift)
        q9_center = direct_signed_hankel(values, 9, shift + 1)
        q9_right = direct_signed_hankel(values, 9, shift + 2)
        q10 = direct_signed_hankel(values, 10, shift)
        width = shift + 9
        schur = direct_deep_schur(values, 10, width)
        relative = q9_center**2 / (q9_left * q9_right) - 1
        toda_rhs = q9_center**2 - q9_left * q9_right
        toda_lhs = q10 * q8
        toda_residual = toda_lhs - toda_rhs
        schur_residual = schur - q10 / values[0] ** 10
        expected = "negative" if shift <= 3 else "positive"
        checks = {
            "q8_positive": sign_class(q8) == "positive",
            "q9_left_positive": sign_class(q9_left) == "positive",
            "q9_center_positive": sign_class(q9_center) == "positive",
            "q9_right_positive": sign_class(q9_right) == "positive",
            "q10_expected_sign": sign_class(q10) == expected,
            "schur_expected_sign": sign_class(schur) == expected,
            "relative_expected_sign": sign_class(relative) == expected,
            "toda_residual_contains_zero": bool(toda_residual.contains(0)),
            "schur_residual_contains_zero": bool(schur_residual.contains(0)),
            "stable_margin_overlaps_direct": bool(
                stable_margins[shift].overlaps(relative)
            ),
        }
        if not all(checks.values()):
            raise RuntimeError(f"direct order-ten audit failed at n={shift}: {checks}")
        rows.append(
            {
                "n": shift,
                "N": width,
                "shape": f"(({width}^10))",
                "expected_sign": expected,
                "Q10_ball": ball_text(q10),
                "deep_schur_ball": ball_text(schur),
                "relative_Q9_margin_ball": ball_text(relative),
                "toda_lhs_ball": ball_text(toda_lhs),
                "toda_rhs_ball": ball_text(toda_rhs),
                "toda_residual_ball": ball_text(toda_residual),
                "schur_residual_ball": ball_text(schur_residual),
                "checks": checks,
            }
        )
    return rows


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    values, sources = load_coefficients()
    stable = stable_order10_scan(values)
    direct = direct_audit(values, stable)
    relative_margins = stable.pop("relative_margins")
    exact = {
        "orientation": "epsilon_10=(-1)^45=-1, Q_(10,n)=-H_(10,n)",
        "toda_identity": (
            "Q_(10,n)*Q_(8,n+2)=Q_(9,n+1)^2-"
            "Q_(9,n)*Q_(9,n+2)"
        ),
        "relative_margin": (
            "L_n=Q_(9,n+1)^2/(Q_(9,n)*Q_(9,n+2))-1"
        ),
        "sign_equivalence": (
            "inside the positive order-eight/order-nine cone, "
            "Q_(10,n)>0 iff L_n>0"
        ),
        "schur_coordinate": (
            "Q_(10,n)(-100)=A_0(-100)^10*s_(((n+9)^10))(h)"
        ),
        "counterexample": (
            "Q_(10,n)(-100)<0 and s_(((n+9)^10))(h)<0 for n=0,1,2,3"
        ),
        "finite_positive_range": (
            "Q_(10,n)(-100)>0 for every 4<=n<=1240"
        ),
        "rejected_hierarchy": (
            "s_((N^m))(h)>0 for every m>=10,N>=m-1 is false"
        ),
        "surviving_bridge": (
            "seek a weaker Xi/Phi-specific route to Jensen hyperbolicity that "
            "does not assume all-shift signed-Hankel positivity"
        ),
    }
    rows = [
        CounterexampleRow(
            "eoc10_01_orientation",
            "exact_identity",
            "ready_to_apply",
            "Signed order ten has negative Hankel orientation.",
            exact["orientation"],
            "Exact parity only.",
        ),
        CounterexampleRow(
            "eoc10_02_interval_input",
            "interval_input",
            "ready_to_apply",
            "Merged retained-integral Arb balls cover every coefficient used by the direct and stable audits.",
            "A_k(-100) rigorously enclosed for 0<=k<=1258",
            "Finite endpoint coefficient range only.",
            {"sources": sources},
        ),
        CounterexampleRow(
            "eoc10_03_toda_coordinate",
            "exact_identity",
            "ready_to_apply",
            "The order-ten sign is exactly the width-log-concavity sign of the positive order-nine row.",
            exact["toda_identity"] + "; " + exact["relative_margin"],
            "Exact condensation inside the already proved lower cone.",
        ),
        CounterexampleRow(
            "eoc10_04_direct_counterexample",
            "interval_counterexample",
            "ready_to_apply",
            "Five direct determinant audits certify four negative order-ten signs followed by one positive sign.",
            exact["counterexample"],
            "Direct endpoint shifts n=0..4 only.",
            {"rows": direct},
        ),
        CounterexampleRow(
            "eoc10_05_schur_counterexample",
            "interval_counterexample",
            "ready_to_apply",
            "The negative signed determinants are required deep rectangular Schur values, not shallow boundary shapes.",
            exact["schur_coordinate"] + "; " + exact["counterexample"],
            "Four deep shapes at order ten.",
        ),
        CounterexampleRow(
            "eoc10_06_finite_sign_map",
            "interval_certificate",
            "ready_to_apply",
            "Cancellation-preserving condensation classifies every cached order-ten row without an inconclusive interval.",
            exact["finite_positive_range"],
            "Finite n<=1240 sign map; no all-tail theorem.",
            {key: value for key, value in stable.items()},
        ),
        CounterexampleRow(
            "eoc10_07_rejected_hierarchy",
            "forbidden_promotion",
            "ready_to_apply",
            "The proposed all-order deep endpoint hierarchy fails at its first open order.",
            exact["rejected_hierarchy"],
            "Rejects this sufficient-condition route, not RH or Jensen hyperbolicity itself.",
        ),
        CounterexampleRow(
            "eoc10_08_rerouted_bridge",
            "open_handoff",
            "not_ready_to_apply",
            "The proof programme must use a weaker determinant family or a direct Xi/Phi-specific Jensen mechanism.",
            exact["surviving_bridge"],
            "New theorem target; no bridge is proved here.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_endpoint_order10_counterexample",
        "date": "2026-07-16",
        "status": (
            "rigorous first-open-order endpoint counterexample rejecting the "
            "all-order deep rectangular hierarchy"
        ),
        "proof_boundary": (
            "This artifact proves four negative signed order-ten endpoint "
            "determinants and rejects the proposed all-order deep rectangular "
            "antecedent. It does not disprove Jensen hyperbolicity, RH, or "
            "Lambda<=0, and it does not supply the replacement Xi/Phi bridge."
        ),
        "sources": sources,
        "exact": exact,
        "direct_audit": direct,
        "finite_scan": stable,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "direct_determinant_checks": len(direct),
            "negative_deep_rectangles": 4,
            "first_positive_shift": 4,
            "scanned_order10_rows": PREFIX_LAST_N + 1,
            "positive_scanned_rows": stable["positive_rows"],
            "negative_scanned_rows": stable["negative_rows"],
            "inconclusive_scanned_rows": stable["inconclusive_rows"],
            "rejected_all_order_endpoint_hierarchies": 1,
            "surviving_xi_specific_bridge_targets": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_endpoint_order10_counterexample.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_endpoint_order10_counterexample.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    direct = artifact["direct_audit"]
    scan = artifact["finite_scan"]
    lines = [
        "# Endpoint Order-Ten Counterexample",
        "",
        "Date: 2026-07-16",
        "",
        "Status: rigorous first-open-order counterexample. This rejects the",
        "proposed all-order deep rectangular endpoint hierarchy; it is not a",
        "counterexample to RH or Jensen hyperbolicity. This is not a proof of",
        "RH or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_endpoint_order10_counterexample.json",
        "python work/rh_compute/scripts/jensen_window_pf_endpoint_order10_counterexample.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_endpoint_order10_counterexample.py",
        "```",
        "",
        "## Exact Coordinate",
        "",
        "```text",
        exact["orientation"],
        exact["toda_identity"],
        exact["relative_margin"],
        exact["sign_equivalence"],
        exact["schur_coordinate"],
        "```",
        "",
        "Orders eight and nine are already positive at every endpoint shift.",
        "Therefore the sign of `L_n` is exactly the sign of `Q_(10,n)`.",
        "",
        "## Rigorous Failure",
        "",
        "Direct `4096`-bit Arb determinants give:",
        "",
        "```text",
    ]
    for row in direct:
        lines.append(
            f"n={row['n']}, N={row['N']}, {row['shape']}: "
            f"L_n={row['relative_Q9_margin_ball']}; sign={row['expected_sign']}"
        )
    lines.extend(
        [
            "```",
            "",
            "Thus the four required deep rectangles",
            "`(9^10),(10^10),(11^10),(12^10)` are strictly negative. The",
            "next rectangle `(13^10)` is strictly positive. These are deep",
            "shapes: their smallest part satisfies `N>=10-1`.",
            "",
            "For every direct row, the checker independently verifies the raw",
            "Hankel determinant, the Jacobi-Trudi determinant, the Toda identity,",
            "and overlap with the stable condensation coordinate.",
            "",
            "## Finite Sign Map",
            "",
            "The cancellation-preserving chain through `Q_9` classifies",
            f"`{artifact['summary']['scanned_order10_rows']}` endpoint rows:",
            "",
            "```text",
            f"negative: n={scan['negative_indices']}",
            f"positive: {scan['positive_range'][0]}<=n<={scan['positive_range'][1]}",
            f"inconclusive: {scan['inconclusive_rows']}",
            "```",
            "",
            "The finite positive stretch and the fixed-order eventual-tail theorem",
            "do not repair the all-shift statement: one negative required rectangle",
            "is enough to refute it.",
            "",
            "## Consequence",
            "",
            "```text",
            exact["rejected_hierarchy"],
            "```",
            "",
            "The endpoint-to-heat theorem and Toda identity remain valid exact",
            "conditional statements, but their all-order positivity antecedent is",
            "false for the actual endpoint sequence. The programme must no longer",
            "present completion of the deep cone as the remaining theorem.",
            "",
            "What survives is the separate problem:",
            "",
            "```text",
            exact["surviving_bridge"],
            "```",
            "",
            "This is a route correction, not negative evidence for RH itself.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote endpoint order-ten counterexample: "
        f"{summary['rows']} rows, "
        f"{summary['direct_determinant_checks']} direct checks, "
        f"{summary['negative_deep_rectangles']} negative deep rectangles, "
        f"{summary['positive_scanned_rows']} positive scanned rows, "
        f"{summary['inconclusive_scanned_rows']} inconclusive"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
