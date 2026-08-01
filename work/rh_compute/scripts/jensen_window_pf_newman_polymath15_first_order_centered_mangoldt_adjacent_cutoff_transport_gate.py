#!/usr/bin/env python3
"""Build the centered Mangoldt adjacent-cutoff transport gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "mangoldt_adjacent_cutoff_transport_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "contact_centering": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "mangoldt_contact_centering_gate.json"
    ),
    "adjacent_recurrence": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "adjacent_saddle_recurrence.json"
    ),
    "chart_stability": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "adjacent_chart_stability_certificate.json"
    ),
}

ComplexQ = tuple[Fraction, Fraction]
LogVector = dict[int, ComplexQ]


@dataclass(frozen=True)
class GateRow:
    id: str
    claim: str
    readiness: str
    exact_statement: str
    implication: str
    boundary: str


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def complex_text(value: ComplexQ) -> dict[str, str]:
    return {
        "real": fraction_text(value[0]),
        "imag": fraction_text(value[1]),
    }


def add(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return left[0] + right[0], left[1] + right[1]


def negate(value: ComplexQ) -> ComplexQ:
    return -value[0], -value[1]


def scale(value: ComplexQ, scalar: Fraction) -> ComplexQ:
    return value[0] * scalar, value[1] * scalar


def multiply(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def vector_add(left: LogVector, right: LogVector) -> LogVector:
    result = dict(left)
    for prime, value in right.items():
        result[prime] = add(
            result.get(prime, (Fraction(0), Fraction(0))),
            value,
        )
        if result[prime] == (Fraction(0), Fraction(0)):
            del result[prime]
    return result


def vector_negate(value: LogVector) -> LogVector:
    return {prime: negate(coefficient) for prime, coefficient in value.items()}


def vector_subtract(left: LogVector, right: LogVector) -> LogVector:
    return vector_add(left, vector_negate(right))


def add_log_vector(
    target: LogVector,
    source: dict[int, int],
    value: ComplexQ,
    scalar: Fraction = Fraction(1),
) -> None:
    for prime, exponent in source.items():
        contribution = scale(value, scalar * exponent)
        target[prime] = add(
            target.get(prime, (Fraction(0), Fraction(0))),
            contribution,
        )
        if target[prime] == (Fraction(0), Fraction(0)):
            del target[prime]


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8", newline="\n")
    temporary.replace(path)


def source_audit() -> dict:
    payloads: dict[str, dict] = {}
    hashes: dict[str, str] = {}
    for key, path in SOURCES.items():
        if not path.is_file():
            raise RuntimeError(f"missing source {key}: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))
        hashes[key] = file_hash(path)

    texts = {
        key: json.dumps(payload, sort_keys=True)
        for key, payload in payloads.items()
    }
    required = {
        "contact_centering": (
            "Delta Z_Lambda=log(N+1)z_(N+1)",
            "mathfrak B_N",
            "balanced",
        ),
        "adjacent_recurrence": (
            "Q_N=f_(N+1)+kappa_N*J_a",
            "Delta A_a=Re[-s_*'*ell*f_(N+1)",
        ),
        "chart_stability": (
            "|Delta X|<1100*exp(-5L/4)",
            "|Delta U|<2500*exp(-7L/4)",
        ),
    }
    # The first marker is the current checkpoint obligation rather than an
    # upstream theorem string; only require the established centering pieces.
    contact_markers = ("mathfrak B_N", "balanced_hyperbola")
    for marker in contact_markers:
        if marker not in texts["contact_centering"]:
            raise RuntimeError(f"centering source marker missing: {marker}")
    for key in ("adjacent_recurrence", "chart_stability"):
        for marker in required[key]:
            if marker not in texts[key]:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {"source_sha256": hashes}


def symbolic_audit() -> dict[str, str]:
    delta_e, delta_g, z_n = sp.symbols("delta_e delta_g z_n")
    s_prime = sp.symbols("s_prime")
    log_a, log_n = sp.symbols("log_a log_n")
    delta_endpoint = delta_g - s_prime * log_a * delta_e
    delta_mangoldt = log_n * z_n
    delta_value = z_n + delta_e
    delta_centered = sp.expand(
        delta_endpoint
        - s_prime * delta_mangoldt
        + s_prime * log_a * delta_value
    )
    expected = delta_g + s_prime * (log_a - log_n) * z_n
    if sp.expand(delta_centered - expected) != 0:
        raise RuntimeError("centered adjacent cancellation failed")

    return {
        "anchored_component_jumps": (
            "At one fixed (t,x), eta, |f_1|, s_*', and L_a=log(a) are "
            "cutoff independent. For n=N+1, "
            "Delta e=kappa_NJ_a/|f_1|, "
            "Delta g=kappa_N[J_(a,x)+mu_aJ_a]/|f_1|, "
            "Delta W_0=z_n+Delta e, and "
            "Delta Z_Lambda=log(n)z_n."
        ),
        "centered_endpoint_jump": (
            "For E_N=g_N-s_*'L_a e_N, "
            "Delta E_N=Delta g-s_*'L_a Delta e."
        ),
        "centered_jet_jump": (
            "The endpoint-value jump cancels exactly: "
            "Delta mathfrak B_N=Delta E_N-s_*'Delta Z_Lambda"
            "+s_*'L_aDelta W_0"
            "=Delta g+s_*'log(a/n)z_n."
        ),
        "physical_recurrence_match": (
            "With ell=log(n/a), z_n=f_n/|f_1|, "
            "Delta mathfrak B_N="
            "{-s_*'ell f_n+kappa_N[J_(a,x)+mu_aJ_a]}/|f_1|"
            "=Delta W_A. Thus Re(Delta mathfrak B_N)"
            "=Delta A_a/|f_1|."
        ),
        "real_projection_bound": (
            "The certified adjacent theorem gives "
            "|Re(Delta W_0)|<2200exp(-5L/4) and "
            "|Re(Delta mathfrak B_N)|<10000exp(-7L/4)."
        ),
        "complex_boundary": (
            "No complex absolute-value bound for Delta W_0 or "
            "Delta mathfrak B_N at the displayed real-projection scales "
            "is proved. The exact complex cancellation must precede "
            "real projection and absolute values."
        ),
        "ordinary_split_jump": (
            "If n=N+1 is not a square, D=floor(sqrt(N))="
            "floor(sqrt(n)); the complete square is unchanged and the "
            "wing gains exactly the factor pairs of n. Hence "
            "Delta Z_square=0 and Delta Z_wing=log(n)z_n."
        ),
        "square_split_jump": (
            "If n=r^2, D changes from r-1 to r. Put "
            "R_r=sum_(d<r)[Lambda(d)+Lambda(r)]z_(dr), "
            "D_r=Lambda(r)z_(r^2), and "
            "A_n=sum_(d|n,d<r)[Lambda(d)+Lambda(n/d)]z_n. Then "
            "Delta Z_square=R_r+D_r and "
            "Delta Z_wing=A_n-R_r."
        ),
        "square_recombination": (
            "At n=r^2 the transferred row R_r cancels, and "
            "D_r+A_n=log(n)z_n by the symmetric Mangoldt divisor "
            "identity. Perfect-square changes introduce no residual term."
        ),
        "route_decision": (
            "The centered jet and balanced hyperbola are now compatible "
            "with adjacent cutoffs. Any Type-I/II or Vaughan estimate must "
            "retain the ordinary factor-pair insertion, the perfect-square "
            "row transfer, the recurrent endpoint slope, and real-projection "
            "semantics before absolute values."
        ),
    }


def factor_vector(number: int) -> dict[int, int]:
    return {
        int(prime): int(exponent)
        for prime, exponent in sp.factorint(number).items()
    }


def lambda_vector(number: int) -> dict[int, int]:
    factors = factor_vector(number)
    if number <= 1 or len(factors) != 1:
        return {}
    return {next(iter(factors)): 1}


def deterministic_carrier(number: int) -> ComplexQ:
    return (
        Fraction((7 * number) % 19 - 9, number + 5),
        Fraction((11 * number) % 23 - 11, number + 7),
    )


def balanced_parts(cutoff: int) -> tuple[LogVector, LogVector]:
    split = math.isqrt(cutoff)
    square: LogVector = {}
    wing: LogVector = {}
    for divisor in range(1, split + 1):
        for other in range(1, split + 1):
            value = deterministic_carrier(divisor * other)
            add_log_vector(
                square,
                lambda_vector(divisor),
                value,
                Fraction(1, 2),
            )
            add_log_vector(
                square,
                lambda_vector(other),
                value,
                Fraction(1, 2),
            )
        for other in range(split + 1, cutoff // divisor + 1):
            value = deterministic_carrier(divisor * other)
            add_log_vector(wing, lambda_vector(divisor), value)
            add_log_vector(wing, lambda_vector(other), value)
    return square, wing


def direct_boundary(number: int) -> LogVector:
    result: LogVector = {}
    add_log_vector(
        result,
        factor_vector(number),
        deterministic_carrier(number),
    )
    return result


def square_transfer(number: int) -> tuple[LogVector, LogVector, LogVector]:
    root = math.isqrt(number)
    if root * root != number:
        raise ValueError("square transfer requires a perfect square")
    moved: LogVector = {}
    diagonal: LogVector = {}
    boundary: LogVector = {}
    for divisor in range(1, root):
        value = deterministic_carrier(divisor * root)
        add_log_vector(moved, lambda_vector(divisor), value)
        add_log_vector(moved, lambda_vector(root), value)
    add_log_vector(
        diagonal,
        lambda_vector(root),
        deterministic_carrier(number),
    )
    for divisor in sp.divisors(number):
        divisor = int(divisor)
        if divisor >= root:
            continue
        other = number // divisor
        value = deterministic_carrier(number)
        add_log_vector(boundary, lambda_vector(divisor), value)
        add_log_vector(boundary, lambda_vector(other), value)
    return moved, diagonal, boundary


def transition_audit(first_cutoff: int = 2, last_cutoff: int = 80) -> list[dict]:
    rows: list[dict] = []
    for cutoff in range(first_cutoff, last_cutoff + 1):
        next_cutoff = cutoff + 1
        old_square, old_wing = balanced_parts(cutoff)
        new_square, new_wing = balanced_parts(next_cutoff)
        delta_square = vector_subtract(new_square, old_square)
        delta_wing = vector_subtract(new_wing, old_wing)
        total = vector_add(delta_square, delta_wing)
        expected = direct_boundary(next_cutoff)
        if total != expected:
            raise RuntimeError(
                f"total cutoff transition failed at N={cutoff}"
            )

        is_square = math.isqrt(next_cutoff) ** 2 == next_cutoff
        if not is_square:
            if delta_square or delta_wing != expected:
                raise RuntimeError(
                    f"ordinary cutoff transition failed at N={cutoff}"
                )
            row_transfer_coordinates = 0
        else:
            moved, diagonal, boundary = square_transfer(next_cutoff)
            expected_square = vector_add(moved, diagonal)
            expected_wing = vector_subtract(boundary, moved)
            if delta_square != expected_square:
                raise RuntimeError(
                    f"square central transition failed at N={cutoff}"
                )
            if delta_wing != expected_wing:
                raise RuntimeError(
                    f"square wing transition failed at N={cutoff}"
                )
            if vector_add(diagonal, boundary) != expected:
                raise RuntimeError(
                    f"square boundary recombination failed at N={cutoff}"
                )
            row_transfer_coordinates = len(moved)

        rows.append(
            {
                "N": cutoff,
                "n": next_cutoff,
                "transition": "perfect_square" if is_square else "ordinary",
                "D_before": math.isqrt(cutoff),
                "D_after": math.isqrt(next_cutoff),
                "delta_square_coordinates": len(delta_square),
                "delta_wing_coordinates": len(delta_wing),
                "row_transfer_coordinates": row_transfer_coordinates,
                "boundary_coordinates": len(expected),
                "mismatches": 0,
            }
        )
    return rows


def rational_jump_witness() -> dict:
    kappa = (Fraction(2, 7), Fraction(-1, 5))
    j_value = (Fraction(-3, 11), Fraction(4, 13))
    j_x = (Fraction(5, 17), Fraction(-2, 19))
    mu = (Fraction(1, 23), Fraction(3, 29))
    f_n = (Fraction(-4, 31), Fraction(5, 37))
    modulus = Fraction(7, 3)
    s_prime = (Fraction(1, 17), Fraction(-1, 2))
    log_a = Fraction(7, 6)
    log_n = Fraction(5, 6)
    ell = log_n - log_a

    delta_e = scale(multiply(kappa, j_value), 1 / modulus)
    delta_g = scale(
        multiply(kappa, add(j_x, multiply(mu, j_value))),
        1 / modulus,
    )
    z_n = scale(f_n, 1 / modulus)
    delta_z_lambda = scale(z_n, log_n)
    delta_w_0 = add(z_n, delta_e)
    delta_e_centered = add(
        delta_g,
        negate(multiply(s_prime, scale(delta_e, log_a))),
    )
    delta_centered = add(
        delta_e_centered,
        add(
            negate(multiply(s_prime, delta_z_lambda)),
            multiply(s_prime, scale(delta_w_0, log_a)),
        ),
    )
    reduced = add(
        delta_g,
        multiply(s_prime, scale(z_n, log_a - log_n)),
    )
    physical = scale(
        add(
            negate(multiply(s_prime, scale(f_n, ell))),
            multiply(kappa, add(j_x, multiply(mu, j_value))),
        ),
        1 / modulus,
    )
    if delta_centered != reduced or reduced != physical:
        raise RuntimeError("rational centered jump witness failed")

    return {
        "kappa_N": complex_text(kappa),
        "J_a": complex_text(j_value),
        "J_a_x": complex_text(j_x),
        "mu_a": complex_text(mu),
        "f_n": complex_text(f_n),
        "abs_f_1": fraction_text(modulus),
        "s_prime": complex_text(s_prime),
        "log_a": fraction_text(log_a),
        "log_n": fraction_text(log_n),
        "ell": fraction_text(ell),
        "Delta_e": complex_text(delta_e),
        "Delta_g": complex_text(delta_g),
        "z_n": complex_text(z_n),
        "Delta_Z_Lambda": complex_text(delta_z_lambda),
        "Delta_W_0": complex_text(delta_w_0),
        "Delta_E_N": complex_text(delta_e_centered),
        "Delta_mathfrak_B_direct": complex_text(delta_centered),
        "Delta_mathfrak_B_reduced": complex_text(reduced),
        "Delta_mathfrak_B_physical": complex_text(physical),
        "interpretation": (
            "Exact rational-complex recurrence witness only; it is not an "
            "Xi cutoff state or a signed lower-bound witness."
        ),
    }


def build_rows(exact: dict[str, str]) -> list[GateRow]:
    return [
        GateRow(
            "mactg_01_components",
            "The anchored endpoint, carrier value, and Mangoldt moment have exact adjacent jumps.",
            "proved",
            exact["anchored_component_jumps"],
            "All components use the same cutoff-independent anchor.",
            "No sign is attached to an individual jump.",
        ),
        GateRow(
            "mactg_02_endpoint",
            "The centered endpoint jump retains the endpoint-value term.",
            "proved",
            exact["centered_endpoint_jump"],
            "The term is kept until recombination.",
            "Bounding it separately would lose cancellation.",
        ),
        GateRow(
            "mactg_03_cancellation",
            "The endpoint-value jump cancels exactly in the centered jet.",
            "proved",
            exact["centered_jet_jump"],
            "Only the endpoint-slope and entering-carrier slope remain.",
            "This is an identity, not a lower bound.",
        ),
        GateRow(
            "mactg_04_recurrence",
            "The centered jump is exactly the established physical adjacent recurrence.",
            "proved",
            exact["physical_recurrence_match"],
            "No new chart-dependent complex observable is introduced.",
            "The Xi signed scalar remains open.",
        ),
        GateRow(
            "mactg_05_projection",
            "The existing adjacent real-projection theorem transfers to the centered variables.",
            "proved",
            exact["real_projection_bound"],
            "The jump is negligible at the contact scale.",
            "Only the real projection is certified at this scale.",
        ),
        GateRow(
            "mactg_06_complex_guard",
            "Complex absolute-value promotion across adjacent charts remains invalid.",
            "guard_validated",
            exact["complex_boundary"],
            "Cancellation is performed before projection.",
            "No complex small-jump theorem is asserted.",
        ),
        GateRow(
            "mactg_07_ordinary",
            "A non-square cutoff increment changes only the hyperbolic wing.",
            "proved",
            exact["ordinary_split_jump"],
            "The new factor pairs reproduce log(n)z_n exactly.",
            "Floor boundaries are retained.",
        ),
        GateRow(
            "mactg_08_square",
            "A perfect-square increment has an explicit row-transfer decomposition.",
            "proved",
            exact["square_split_jump"],
            "The central and wing changes are separately exact.",
            "Neither change is a signed estimate.",
        ),
        GateRow(
            "mactg_09_square_rejoin",
            "The perfect-square row transfer cancels under recombination.",
            "proved",
            exact["square_recombination"],
            "No residual floor term survives.",
            "This is arithmetic bookkeeping, not contact coercivity.",
        ),
        GateRow(
            "mactg_10_transition_audit",
            "Ordinary and perfect-square transitions pass an exact finite coefficient audit.",
            "finite_certificate",
            (
                "Every N=2,...,80 transition is reconstructed with "
                "deterministic rational-complex carriers and prime-log "
                "coefficient vectors."
            ),
            "The formulas are checked independently of floating-point logs.",
            "Finite audit support for exact identities only.",
        ),
        GateRow(
            "mactg_11_band",
            "Adjacent transport does not consume the existing Abel-gap scale.",
            "proved",
            exact["real_projection_bound"],
            "The certified exp(-7L/4) jump is subordinate to exp(-5L/4).",
            "No fixed-chart Abel gap is proved.",
        ),
        GateRow(
            "mactg_12_route",
            "The centered Type-I/II coordinate is now cutoff compatible.",
            "open",
            exact["route_decision"],
            "The next step may formulate one endpoint-composed estimate.",
            (
                "No Abel gap, winding cap, contact exclusion, Lambda<=0, "
                "PF-infinity, RH, or prize-level conclusion is proved."
            ),
        ),
    ]


def build_note(payload: dict) -> str:
    exact = payload["exact"]
    summary = payload["summary"]
    return f"""# Centered Mangoldt Adjacent-Cutoff Transport Gate

Date: 2026-07-30

Status: exact componentwise adjacent transport and square/wing transition audit;
not a proof artifact. This is not a signed centered-jet lower bound, Abel
gap, contact exclusion, `Lambda<=0`, PF-infinity, RH, or a Clay-prize result.

```text
work/rh_compute/results/{STEM}.json
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```

## Pi Provenance

The `pi` in the endpoint recurrence and `a^2=x/(4*pi)+t/16` remains the
completed-zeta and Riemann-Siegel normalization. The Mangoldt cutoff,
integer square root, and row transfer introduce no new `pi`.

## Component Jumps

```text
{exact["anchored_component_jumps"]}

{exact["centered_endpoint_jump"]}

{exact["centered_jet_jump"]}
```

The two appearances of `Delta e` cancel before an absolute value is taken.

```text
{exact["physical_recurrence_match"]}
```

Thus the new centered coordinate has exactly the old physical adjacent
jet; it is not a second transport problem.

## Projection Semantics

```text
{exact["real_projection_bound"]}

{exact["complex_boundary"]}
```

The real trace is the proof-facing quantity. A small complex adjacent
jump is neither needed nor claimed.

## Ordinary Cutoffs

```text
{exact["ordinary_split_jump"]}
```

## Perfect Squares

```text
{exact["square_split_jump"]}

{exact["square_recombination"]}
```

The square/wing boundary moves at a perfect square, but the movement is
an internal transfer. It leaves no extra term in `Delta Z_Lambda`.

## Exact Audit

The independent audit checks `{summary["cutoff_transitions"]}` successive
cutoffs, comprising `{summary["ordinary_transitions"]}` ordinary and
`{summary["perfect_square_transitions"]}` perfect-square transitions.
It uses rational-complex carriers and formal prime-log coefficient
vectors, with `{summary["transition_mismatches"]}` mismatches.

## Route Decision

```text
{exact["route_decision"]}
```

The next estimate may now be written on the balanced hyperbola without
silently changing the object at a cutoff. It must still be tested at
`W_0=0`, `mathsf_X=0`, `H_a=0`, `q=1`, prime edges, ordinary transitions,
perfect-square transitions, and adjacent real projection.

## Boundary

This gate proves the centered component jumps, exact endpoint-value
cancellation, recurrence match, inherited real-projection bound, ordinary
wing insertion, perfect-square row transfer and cancellation, and finite
coefficient audit. It proves no signed lower bound, Abel-scalar gap,
horizontal successor winding cap, contact exclusion, Q209, cofinal
descendant theorem, `Lambda<=0`, PF-infinity, RH, or prize-level result.
"""


def build_payload() -> dict:
    exact = symbolic_audit()
    transitions = transition_audit()
    ordinary = sum(row["transition"] == "ordinary" for row in transitions)
    squares = sum(
        row["transition"] == "perfect_square" for row in transitions
    )
    rows = build_rows(exact)
    return {
        "kind": "mangoldt_adjacent_cutoff_transport_gate",
        "date": "2026-07-30",
        "status": (
            "exact centered adjacent transport and balanced-cutoff audit; "
            "signed lower bound remains open"
        ),
        "source_audit": source_audit(),
        "exact": exact,
        "witnesses": {
            "rational_jump": rational_jump_witness(),
            "cutoff_transitions": transitions,
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "centered_endpoint_cancellations": 1,
            "cutoff_transitions": len(transitions),
            "ordinary_transitions": ordinary,
            "perfect_square_transitions": squares,
            "transition_mismatches": sum(
                int(row["mismatches"]) for row in transitions
            ),
            "proved_abel_gaps": 0,
            "proved_successor_winding_bounds": 0,
        },
        "proof_boundary": (
            "This gate proves exact centered component jumps, endpoint-value "
            "cancellation, physical recurrence matching, inherited adjacent "
            "real-projection control, and ordinary/perfect-square balanced "
            "hyperbola transport. It proves no signed lower bound, Abel gap, "
            "winding cap, contact exclusion, Lambda<=0, PF-infinity, RH, or "
            "prize-level conclusion."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=RESULT)
    parser.add_argument("--note", type=Path, default=NOTE)
    args = parser.parse_args()

    payload = build_payload()
    atomic_write(
        args.result,
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
    )
    atomic_write(args.note, build_note(payload))
    summary = payload["summary"]
    print(
        "built centered Mangoldt cutoff transport gate: "
        f"{summary['rows']} rows, "
        f"{summary['cutoff_transitions']} transitions, "
        f"{summary['perfect_square_transitions']} square transitions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
