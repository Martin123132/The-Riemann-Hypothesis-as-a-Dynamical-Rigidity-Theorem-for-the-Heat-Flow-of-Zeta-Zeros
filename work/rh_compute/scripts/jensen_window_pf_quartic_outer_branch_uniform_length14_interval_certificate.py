#!/usr/bin/env python3
"""Certify a uniform length-14 obstruction for one outer quartic contact."""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


PRECISION_BITS = 256
VARIABLE_COUNT = 8
ROOT_Y6_UPPER = Fraction(
    9865549779980282617079853707566678613195505108205132800,
    26486057908576767522100028185953787483471843406523808729,
)
DEFAULT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_quartic_outer_branch_uniform_length14_interval_events_endpoint.jsonl"
)
DEFAULT_PROGRESS = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_quartic_outer_branch_uniform_length14_interval_progress_endpoint.json"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_quartic_outer_branch_uniform_length14_interval_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_quartic_outer_branch_uniform_length14_interval_certificate.md"
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


class Dual:
    """Arb interval value with an eight-component forward derivative."""

    __slots__ = ("value", "gradient")

    def __init__(self, value, gradient=None):
        self.value = value if isinstance(value, flint.arb) else arb_exact(value)
        self.gradient = (
            gradient
            if gradient is not None
            else [flint.arb(0) for _ in range(VARIABLE_COUNT)]
        )

    def __add__(self, other):
        other = other if isinstance(other, Dual) else Dual(other)
        return Dual(
            self.value + other.value,
            [
                left + right
                for left, right in zip(self.gradient, other.gradient)
            ],
        )

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.value, [-entry for entry in self.gradient])

    def __sub__(self, other):
        return self + (-(other if isinstance(other, Dual) else Dual(other)))

    def __rsub__(self, other):
        return Dual(other) - self

    def __mul__(self, other):
        other = other if isinstance(other, Dual) else Dual(other)
        return Dual(
            self.value * other.value,
            [
                left * other.value + self.value * right
                for left, right in zip(self.gradient, other.gradient)
            ],
        )

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = other if isinstance(other, Dual) else Dual(other)
        denominator = other.value**2
        return Dual(
            self.value / other.value,
            [
                (left * other.value - self.value * right) / denominator
                for left, right in zip(self.gradient, other.gradient)
            ],
        )

    def __rtruediv__(self, other):
        return Dual(other) / self

    def __pow__(self, exponent: int):
        return Dual(
            self.value**exponent,
            [
                flint.arb(exponent)
                * self.value ** (exponent - 1)
                * entry
                for entry in self.gradient
            ],
        )


def lower_process_priority() -> None:
    try:
        import psutil

        psutil.Process(os.getpid()).nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:
        pass


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def arb_exact(value: Fraction | int) -> flint.arb:
    if isinstance(value, int):
        return flint.arb(value)
    return flint.arb(value.numerator) / flint.arb(value.denominator)


def arb_hull(lower: Fraction, upper: Fraction) -> flint.arb:
    return flint.arb.union(arb_exact(lower), arb_exact(upper))


def strictly_positive(value: flint.arb) -> bool:
    return value.lower() > 0


def strictly_negative(value: flint.arb) -> bool:
    return value.upper() < 0


def serialize_arb(value: flint.arb) -> str:
    """Preserve the midpoint when the radius nearly cancels it."""
    return value.str(n=80, radius=True, more=True)


def base_data() -> tuple[dict[int, Fraction], Fraction, Fraction]:
    contractions = {
        2: Fraction(24154639, 25000000),
        3: Fraction(567331181410000, 583446585220321),
        4: Fraction(909681023616163, 930852655574884),
        5: Fraction(2453, 2500),
    }
    defects = {index: 1 - value for index, value in contractions.items()}
    gap_1 = (
        defects[3] ** 2
        - contractions[3] ** 2 * defects[2] * defects[4]
    )
    gap_2 = (
        defects[4] ** 2
        - contractions[4] ** 2 * defects[3] * defects[5]
    )
    return contractions, gap_1, gap_2


def derive_root_y6_upper() -> Fraction:
    contractions, gap_1, gap_2 = base_data()
    defects = {
        index: 1 - value for index, value in contractions.items()
    }
    cap_6 = gap_2**2 / (contractions[4] ** 3 * gap_1)
    return (
        defects[5] ** 2
        - Fraction(11, 13)
        * defects[5]
        * contractions[5] ** 2
        * defects[4]
    ) / cap_6


def exponent_table() -> dict[int, dict[str, int]]:
    """Return cancellation-free monomial exponents for G_1 through G_11."""
    table: dict[int, defaultdict[str, int]] = {
        1: defaultdict(int, {"g1": 1}),
        2: defaultdict(int, {"g2": 1}),
    }
    for index in range(3, 12):
        row: defaultdict[str, int] = defaultdict(int)
        for key, exponent in table[index - 1].items():
            row[key] += 2 * exponent
        for key, exponent in table[index - 2].items():
            row[key] -= exponent
        row[f"y{index + 3}"] += 1
        row[f"x{index + 1}"] -= 3
        table[index] = row
    return {
        index: {
            key: exponent
            for key, exponent in sorted(row.items())
            if exponent
        }
        for index, row in table.items()
    }


EXPONENTS = exponent_table()


def root_box() -> tuple[tuple[Fraction, Fraction], ...]:
    return ((Fraction(0), ROOT_Y6_UPPER),) + tuple(
        (Fraction(0), Fraction(1)) for _ in range(7)
    )


def decode_box(path: str) -> tuple[tuple[Fraction, Fraction], ...]:
    box = [list(interval) for interval in root_box()]
    for token in path:
        encoded = int(token, 16)
        dimension = encoded // 2
        side = encoded % 2
        lower, upper = box[dimension]
        midpoint = (lower + upper) / 2
        box[dimension] = (
            [lower, midpoint] if side == 0 else [midpoint, upper]
        )
    return tuple((lower, upper) for lower, upper in box)


def child_path(path: str, dimension: int, side: int) -> str:
    return path + format(2 * dimension + side, "x")


def evaluate_factored(
    variables: dict[int, Dual],
) -> Dual:
    fixed, gap_1_fraction, gap_2_fraction = base_data()
    contractions = {index: Dual(value) for index, value in fixed.items()}
    defects = {
        index: Dual(1) - value for index, value in contractions.items()
    }
    gap_1 = Dual(gap_1_fraction)
    gap_2 = Dual(gap_2_fraction)

    def monomial(index: int, omit_y: int | None = None) -> Dual:
        value = Dual(1)
        for key, exponent in EXPONENTS[index].items():
            if key == "g1":
                factor = gap_1
            elif key == "g2":
                factor = gap_2
            elif key.startswith("x"):
                factor = contractions[int(key[1:])]
            else:
                parameter_index = int(key[1:])
                factor = (
                    Dual(1)
                    if parameter_index == omit_y
                    else variables[parameter_index]
                )
            value *= factor**exponent
        return value

    for index in range(6, 14):
        current_gap = monomial(index - 3)
        defects[index] = (
            defects[index - 1] ** 2 - current_gap
        ) / (
            contractions[index - 1] ** 2 * defects[index - 2]
        )
        contractions[index] = Dual(1) - defects[index]

    cap_14 = monomial(11, omit_y=14)
    repeated_gap = (
        defects[13] ** 2
        - contractions[13] ** 2 * defects[12] * defects[13]
    )
    return cap_14 - repeated_gap


def evaluate_box(
    box: tuple[tuple[Fraction, Fraction], ...],
) -> Dual:
    variables: dict[int, Dual] = {}
    for offset, (lower, upper) in enumerate(box):
        gradient = [flint.arb(0) for _ in range(VARIABLE_COUNT)]
        gradient[offset] = flint.arb(1)
        variables[offset + 6] = Dual(arb_hull(lower, upper), gradient)
    return evaluate_factored(variables)


def evaluate_point(point: list[Fraction]) -> flint.arb:
    variables = {
        offset + 6: Dual(arb_exact(value))
        for offset, value in enumerate(point)
    }
    return evaluate_factored(variables).value


def assess_box(
    box: tuple[tuple[Fraction, Fraction], ...],
) -> dict[str, object] | None:
    """Return a rigorous monotone-face mean-value upper enclosure."""
    evaluated = evaluate_box(box)
    if not evaluated.value.is_finite() or any(
        not entry.is_finite() for entry in evaluated.gradient
    ):
        return None

    anchor: list[Fraction] = []
    uncertain: list[int] = []
    contributions: list[flint.arb] = []
    monotone_signs: list[int] = []
    for dimension, (derivative, interval) in enumerate(
        zip(evaluated.gradient, box)
    ):
        lower, upper = interval
        if strictly_positive(derivative):
            anchor.append(upper)
            monotone_signs.append(1)
            contributions.append(flint.arb(0))
        elif strictly_negative(derivative):
            anchor.append(lower)
            monotone_signs.append(-1)
            contributions.append(flint.arb(0))
        else:
            anchor.append((lower + upper) / 2)
            monotone_signs.append(0)
            uncertain.append(dimension)
            contributions.append(
                abs(derivative) * arb_exact((upper - lower) / 2)
            )

    center_value = evaluate_point(anchor)
    mean_upper = center_value
    for contribution in contributions:
        mean_upper += contribution
    if not mean_upper.is_finite():
        return None
    contribution_upper = [
        float(entry.upper()) for entry in contributions
    ]
    return {
        "mean_upper": mean_upper,
        "center_value": center_value,
        "monotone_signs": monotone_signs,
        "uncertain": uncertain,
        "contribution_upper": contribution_upper,
    }


def choose_split(
    box: tuple[tuple[Fraction, Fraction], ...],
    assessment: dict[str, object] | None,
) -> tuple[int, str]:
    if assessment is None:
        normalized_widths = [
            float(
                (upper - lower)
                / (ROOT_Y6_UPPER if dimension == 0 else Fraction(1))
            )
            for dimension, (lower, upper) in enumerate(box)
        ]
        return (
            max(range(VARIABLE_COUNT), key=normalized_widths.__getitem__),
            "nonfinite_interval",
        )
    contributions = assessment["contribution_upper"]
    if max(contributions) > 0:
        return (
            max(range(VARIABLE_COUNT), key=contributions.__getitem__),
            "largest_uncertain_mean_value_contribution",
        )
    normalized_widths = [
        float(
            (upper - lower)
            / (ROOT_Y6_UPPER if dimension == 0 else Fraction(1))
        )
        for dimension, (lower, upper) in enumerate(box)
    ]
    return (
        max(range(VARIABLE_COUNT), key=normalized_widths.__getitem__),
        "unresolved_monotone_face_rounding",
    )


def event_body(
    index: int,
    path: str,
    action: str,
    **fields,
) -> dict[str, object]:
    return {"index": index, "path": path, "action": action, **fields}


def event_hash(previous: str, body: dict[str, object]) -> str:
    encoded = json.dumps(
        body, sort_keys=True, separators=(",", ":")
    ).encode("ascii")
    return hashlib.sha256(previous.encode("ascii") + encoded).hexdigest()


def load_events(
    cache: Path,
) -> tuple[list[str], dict[str, object], str]:
    stack = [""]
    previous = "0" * 64
    stats: dict[str, object] = {
        "events": 0,
        "splits": 0,
        "leaves": 0,
        "max_depth": 0,
        "worst_leaf_upper": -math.inf,
        "nonfinite_splits": 0,
    }
    if not cache.exists():
        return stack, stats, previous
    with cache.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.endswith("\n"):
                raise RuntimeError(
                    f"incomplete event-log line {line_number}; preserve "
                    "the cache and repair only the terminal partial line"
                )
            event = json.loads(line)
            if event.get("index") != stats["events"]:
                raise RuntimeError("event index discontinuity")
            if event.get("prev_sha256") != previous:
                raise RuntimeError("event hash-chain predecessor changed")
            body = {
                key: value
                for key, value in event.items()
                if key not in {"prev_sha256", "event_sha256"}
            }
            observed_hash = event_hash(previous, body)
            if event.get("event_sha256") != observed_hash:
                raise RuntimeError("event hash-chain digest changed")
            if not stack or event.get("path") != stack.pop():
                raise RuntimeError("event path is not the next DFS box")
            path = str(event["path"])
            stats["max_depth"] = max(stats["max_depth"], len(path))
            if event.get("action") == "split":
                dimension = int(event["dimension"])
                stack.append(child_path(path, dimension, 1))
                stack.append(child_path(path, dimension, 0))
                stats["splits"] += 1
                if event.get("reason") == "nonfinite_interval":
                    stats["nonfinite_splits"] += 1
            elif event.get("action") == "leaf":
                stats["leaves"] += 1
                upper_endpoint = flint.arb(
                    event["mean_upper_upper_endpoint"]
                )
                if not strictly_negative(upper_endpoint):
                    raise RuntimeError(
                        "serialized leaf upper endpoint is not strictly negative"
                    )
                upper = float(upper_endpoint.upper())
                stats["worst_leaf_upper"] = max(
                    stats["worst_leaf_upper"], upper
                )
            else:
                raise RuntimeError("unknown event action")
            stats["events"] += 1
            previous = observed_hash
    return stack, stats, previous


def append_event(
    handle,
    previous: str,
    body: dict[str, object],
) -> str:
    digest = event_hash(previous, body)
    event = {
        **body,
        "prev_sha256": previous,
        "event_sha256": digest,
    }
    handle.write(json.dumps(event, sort_keys=True) + "\n")
    return digest


def exact_corner() -> dict[str, object]:
    if derive_root_y6_upper() != ROOT_Y6_UPPER:
        raise RuntimeError("hard-coded y_6 wall changed from its exact derivation")
    fixed, gap_1, gap_2 = base_data()
    contractions = dict(fixed)
    defects = {index: 1 - value for index, value in contractions.items()}
    variables = {
        6: ROOT_Y6_UPPER,
        **{index: Fraction(1) for index in range(7, 14)},
    }

    def monomial(index: int, omit_y: int | None = None) -> Fraction:
        value = Fraction(1)
        for key, exponent in EXPONENTS[index].items():
            if key == "g1":
                factor = gap_1
            elif key == "g2":
                factor = gap_2
            elif key.startswith("x"):
                factor = contractions[int(key[1:])]
            else:
                parameter_index = int(key[1:])
                factor = (
                    Fraction(1)
                    if parameter_index == omit_y
                    else variables[parameter_index]
                )
            value *= factor**exponent
        return value

    for index in range(6, 14):
        current_gap = monomial(index - 3)
        defects[index] = (
            defects[index - 1] ** 2 - current_gap
        ) / (
            contractions[index - 1] ** 2 * defects[index - 2]
        )
        contractions[index] = 1 - defects[index]
    cap_14 = monomial(11, omit_y=14)
    repeated = (
        defects[13] ** 2
        - contractions[13] ** 2 * defects[12] * defects[13]
    )
    delta = cap_14 - repeated
    encoded = f"{delta.numerator}/{delta.denominator}".encode("ascii")
    return {
        "y6_upper": (
            f"{ROOT_Y6_UPPER.numerator}/{ROOT_Y6_UPPER.denominator}"
        ),
        "x6": (
            f"{contractions[6].numerator}/{contractions[6].denominator}"
        ),
        "delta14_decimal": str(
            flint.arb(delta.numerator) / flint.arb(delta.denominator)
        ),
        "delta14_sign": -1 if delta < 0 else 1 if delta > 0 else 0,
        "delta14_numerator_digits": len(str(abs(delta.numerator))),
        "delta14_denominator_digits": len(str(delta.denominator)),
        "delta14_sha256": hashlib.sha256(encoded).hexdigest(),
    }


def write_progress(
    path: Path,
    cache: Path,
    stack: list[str],
    stats: dict[str, object],
    final_chain: str,
    elapsed: float,
) -> None:
    payload = {
        "kind": (
            "jensen_window_pf_quartic_outer_branch_"
            "uniform_length14_interval_progress"
        ),
        "date": "2026-07-25",
        "status": "complete" if not stack else "resumable_prefix",
        "precision_bits": PRECISION_BITS,
        "event_log": cache.relative_to(REPO_ROOT).as_posix(),
        "event_log_bytes": cache.stat().st_size if cache.exists() else 0,
        "event_log_sha256": file_sha256(cache) if cache.exists() else None,
        "final_event_sha256": final_chain,
        "pending_boxes": len(stack),
        "next_path": stack[-1] if stack else None,
        "elapsed_seconds_this_run": elapsed,
        **stats,
        "resume_command": (
            "python work/rh_compute/scripts/"
            "jensen_window_pf_quartic_outer_branch_"
            "uniform_length14_interval_certificate.py --resume"
        ),
        "proof_boundary": (
            "A nonempty pending stack is an incomplete interval tree and "
            "proves no uniform obstruction."
        ),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def build_payload(cache: Path, progress: Path) -> dict:
    stack, stats, final_chain = load_events(cache)
    if stack:
        raise RuntimeError("cannot promote an incomplete interval tree")
    exact = exact_corner()
    rows = [
        GateRow(
            "qoul14_01_contact",
            "provenance",
            "proved_exact",
            "The theorem fixes the exact outward quartic contact through x_5.",
            "x_2,...,x_5 fixed",
            "One contact only; not the complete outer-contact family.",
        ),
        GateRow(
            "qoul14_02_scaled_wall",
            "exact_reduction",
            "proved_exact",
            "The first scaled-defect inequality is exactly the upper wall y_6<Y_6.",
            "13*d_6>11*d_5 iff y_6<Y_6",
            "Uses the positive order-four cap C_6.",
            {"Y_6": exact["y6_upper"], "boundary_x6": exact["x6"]},
        ),
        GateRow(
            "qoul14_03_factored_gaps",
            "exact_identity",
            "proved_exact",
            "Every corridor gap has a cancellation-free positive monomial representation.",
            "G_j=y_(j+3)*G_(j-1)^2/(x_(j+1)^3*G_(j-2))",
            "Identity on the open corridor; the factored form extends to its closure.",
            {"exponent_table": EXPONENTS},
        ),
        GateRow(
            "qoul14_04_domain",
            "exact_reduction",
            "proved_exact",
            "Every strict signed continuation satisfying the first scaled wall lies in one closed parameter box.",
            "0<=y_6<=Y_6; 0<=y_7,...,y_13<=1",
            "The closed box is used only as a proving superset of the strict corridor.",
        ),
        GateRow(
            "qoul14_05_interval_tree",
            "interval_certificate",
            "interval_validated",
            "A hash-chained dyadic Arb tree covers the complete parameter box.",
            "union certified leaves = [0,Y_6]x[0,1]^7",
            "Coverage follows from deterministic binary replay of every split.",
            {
                "events": stats["events"],
                "splits": stats["splits"],
                "leaves": stats["leaves"],
                "max_depth": stats["max_depth"],
            },
        ),
        GateRow(
            "qoul14_06_mean_value",
            "interval_certificate",
            "interval_validated",
            "Every leaf has a strictly negative monotone-face mean-value upper enclosure for Delta_14.",
            "Delta_14(c)+sum_i sup_B|partial_i Delta_14|*r_i<0",
            "Arb precision and the exact dyadic endpoints are manifest-bound.",
            {"worst_leaf_upper": stats["worst_leaf_upper"]},
        ),
        GateRow(
            "qoul14_07_corner",
            "exact_certificate",
            "proved_exact",
            "The optimizer-indicated closed boundary corner is exactly negative.",
            "Delta_14(Y_6,1,...,1)<0",
            "The corner is explanatory; the interval tree proves the whole box.",
            exact,
        ),
        GateRow(
            "qoul14_08_obstruction",
            "scope_gate",
            "proved_exact",
            "Every strict signed continuation of this contact through x_13 satisfying the first scaled step has Delta_14<0.",
            "Delta_14<0 on the admissible corridor",
            "This forbids an increasing x_14 with the next order-four sign for this contact only.",
        ),
        GateRow(
            "qoul14_09_boundary",
            "scope_gate",
            "proved_exact",
            "The result is uniform over tails from one fixed contact but not over all outer quartic contacts.",
            "one-contact uniformity != Xi quartic closure",
            "No PF-infinity, Lambda<=0, RH, or Clay-prize conclusion.",
        ),
    ]
    return {
        "kind": (
            "jensen_window_pf_quartic_outer_branch_"
            "uniform_length14_interval_certificate"
        ),
        "date": "2026-07-25",
        "status": "rigorous one-contact uniform length-fourteen obstruction",
        "proof_boundary": (
            "This interval-backed theorem proves that every strict signed "
            "order-three/order-four continuation through x_13 from one "
            "specified outward quartic contact, if it satisfies the first "
            "scaled-defect step, has Delta_14<0. It is uniform over the "
            "tail parameters for that contact. It does not cover every "
            "outer contact, prove u<=U for Xi, PF-infinity, Lambda<=0, or RH."
        ),
        "precision_bits": PRECISION_BITS,
        "root_box": {
            "y_6": ["0", exact["y6_upper"]],
            **{f"y_{index}": ["0", "1"] for index in range(7, 14)},
        },
        "exact": exact,
        "interval_certificate": {
            **stats,
            "event_log": cache.relative_to(REPO_ROOT).as_posix(),
            "event_log_bytes": cache.stat().st_size,
            "event_log_sha256": file_sha256(cache),
            "final_event_sha256": final_chain,
            "progress_manifest": progress.relative_to(REPO_ROOT).as_posix(),
            "progress_manifest_sha256": file_sha256(progress),
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "fixed_outer_contacts": 1,
            "tail_parameters": 8,
            "scaled_walls_used": 1,
            "certified_leaves": stats["leaves"],
            "uniform_length14_obstructions_for_fixed_contact": 1,
            "uniform_all_contact_theorems": 0,
        },
    }


def render_note(payload: dict) -> str:
    certificate = payload["interval_certificate"]
    exact = payload["exact"]
    return "\n".join(
        [
            "# Quartic Outer-Branch Uniform Length-14 Interval Certificate",
            "",
            "Date: 2026-07-25",
            "",
            "Status: rigorous one-contact uniform length-fourteen obstruction;",
            "not a proof of the Xi quartic threshold, `Lambda <= 0`, or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_quartic_outer_branch_uniform_length14_interval_certificate.json",
            certificate["event_log"],
            "python work/rh_compute/scripts/jensen_window_pf_quartic_outer_branch_uniform_length14_interval_certificate.py --resume",
            "python work/rh_compute/scripts/check_jensen_window_pf_quartic_outer_branch_uniform_length14_interval_certificate.py",
            "```",
            "",
            "## Exact Reduction",
            "",
            "For the fixed outward contact through `x_5`, the first",
            "scaled-defect inequality is exactly",
            "",
            "```text",
            "13*d_6>11*d_5  iff  0<y_6<Y_6,",
            f"Y_6={exact['y6_upper']}.",
            "```",
            "",
            "All later signed order-three/order-four corridors have",
            "`0<y_k<1`. Cancellation of the gap recurrence gives positive",
            "monomials in the `y_k` and inverse powers of earlier `x_k`, so",
            "the rational recurrence extends continuously to the closed box",
            "",
            "```text",
            "[0,Y_6] x [0,1]^7.",
            "```",
            "",
            "## Rigorous Cover",
            "",
            "The append-only event log is a complete binary partition of",
            "that box. On each leaf, 256-bit Arb automatic differentiation",
            "fixes every coordinate having a sign-definite derivative at its",
            "maximizing face and applies a mean-value enclosure to the",
            "remaining coordinates.",
            "",
            "```text",
            f"events: {certificate['events']}",
            f"splits: {certificate['splits']}",
            f"certified leaves: {certificate['leaves']}",
            f"maximum depth: {certificate['max_depth']}",
            f"largest certified upper endpoint: {certificate['worst_leaf_upper']}",
            "pending boxes: 0",
            "```",
            "",
            "Every upper endpoint is strictly negative. Therefore",
            "`Delta_14<0` throughout the closed proving box.",
            "",
            "The exact optimizer-indicated corner independently gives",
            "",
            "```text",
            f"Delta_14(Y_6,1,...,1)={exact['delta14_decimal']}<0.",
            "```",
            "",
            "## Scope",
            "",
            "This is uniform over every tail from this one fixed outer",
            "contact satisfying the first scaled step. It upgrades two",
            "fixed-tail examples to a complete same-contact obstruction at",
            "length 14. It is not uniform over the complete `(a,p,u)` outer",
            "contact family and does not prove the Xi quartic threshold,",
            "PF-infinity, `Lambda <= 0`, or RH.",
            "",
        ]
    )


def run(
    cache: Path,
    progress: Path,
    out: Path,
    note: Path,
    max_seconds: float,
    max_events: int,
    fsync_every: int,
) -> bool:
    flint.ctx.prec = PRECISION_BITS
    lower_process_priority()
    cache.parent.mkdir(parents=True, exist_ok=True)
    stack, stats, previous = load_events(cache)
    start = time.monotonic()
    written = 0
    with cache.open("a", encoding="utf-8") as handle:
        while (
            stack
            and written < max_events
            and time.monotonic() - start < max_seconds
        ):
            path = stack.pop()
            box = decode_box(path)
            assessment = assess_box(box)
            if (
                assessment is not None
                and strictly_negative(assessment["mean_upper"])
            ):
                body = event_body(
                    int(stats["events"]),
                    path,
                    "leaf",
                    mean_upper_ball=serialize_arb(assessment["mean_upper"]),
                    mean_upper_upper_endpoint=serialize_arb(
                        assessment["mean_upper"].upper()
                    ),
                    center_value_ball=serialize_arb(assessment["center_value"]),
                    monotone_signs=assessment["monotone_signs"],
                )
                stats["leaves"] += 1
                stats["worst_leaf_upper"] = max(
                    float(stats["worst_leaf_upper"]),
                    float(assessment["mean_upper"].upper()),
                )
            else:
                dimension, reason = choose_split(box, assessment)
                body = event_body(
                    int(stats["events"]),
                    path,
                    "split",
                    dimension=dimension,
                    reason=reason,
                )
                stack.append(child_path(path, dimension, 1))
                stack.append(child_path(path, dimension, 0))
                stats["splits"] += 1
                if reason == "nonfinite_interval":
                    stats["nonfinite_splits"] += 1
            previous = append_event(handle, previous, body)
            stats["events"] += 1
            stats["max_depth"] = max(stats["max_depth"], len(path))
            written += 1
            if written % fsync_every == 0:
                handle.flush()
                os.fsync(handle.fileno())
                write_progress(
                    progress,
                    cache,
                    stack,
                    stats,
                    previous,
                    time.monotonic() - start,
                )
        handle.flush()
        os.fsync(handle.fileno())

    elapsed = time.monotonic() - start
    write_progress(
        progress, cache, stack, stats, previous, elapsed
    )
    if stack:
        print(
            "parked uniform length-14 interval tree: "
            f"{stats['events']} events, {stats['leaves']} leaves, "
            f"{len(stack)} pending boxes, next path {stack[-1]!r}"
        )
        return False

    payload = build_payload(cache, progress)
    out.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    note.write_text(render_note(payload), encoding="utf-8")
    print(
        "completed uniform length-14 interval tree: "
        f"{stats['events']} events, {stats['leaves']} certified leaves, "
        f"maximum depth {stats['max_depth']}"
    )
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--progress", type=Path, default=DEFAULT_PROGRESS)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--max-seconds", type=float, default=900)
    parser.add_argument("--max-events", type=int, default=200000)
    parser.add_argument("--fsync-every", type=int, default=100)
    args = parser.parse_args()
    if not args.resume and args.cache.exists() and args.cache.stat().st_size:
        raise RuntimeError("event cache exists; pass --resume")
    run(
        args.cache,
        args.progress,
        args.out,
        args.note,
        args.max_seconds,
        args.max_events,
        args.fsync_every,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
