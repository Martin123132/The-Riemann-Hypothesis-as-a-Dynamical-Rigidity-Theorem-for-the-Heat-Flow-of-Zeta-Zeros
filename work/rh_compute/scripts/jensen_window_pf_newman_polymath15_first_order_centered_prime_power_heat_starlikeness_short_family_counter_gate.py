#!/usr/bin/env python3
"""Build exact short small-prime heat-starlikeness counter-witnesses."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_short_family_counter_gate.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_short_family_counter_gate.md"
)
SOURCES = {
    "base": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "prime_power_heat_starlikeness_base_certificate.json"
    ),
    "geometric": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "geometric_prime_power_phase_monotonicity_gate.json"
    ),
    "occupation": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "complete_prime_power_chain_occupation_gate.json"
    ),
}

Q_STAR = Fraction(9999, 10000)
NEGATIVE_TARGET = Fraction(1, 100)
Pair = tuple[Fraction, Fraction]
WITNESSES = (
    (2, 2, Fraction(-99, 100)),
    (2, 3, Fraction(-1, 2)),
    (2, 4, Fraction(0)),
    (2, 5, Fraction(1, 3)),
    (2, 6, Fraction(-1, 2)),
    (2, 7, Fraction(-1, 4)),
    (3, 2, Fraction(-99, 100)),
    (3, 3, Fraction(-1, 2)),
)
BOXES = {
    2: (Fraction(47, 50), Fraction(22, 25)),
    3: (Fraction(17, 20), Fraction(73, 100)),
}


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


def pair_text(value: Pair) -> dict[str, str]:
    return {
        "rational": fraction_text(value[0]),
        "sqrt_prime_coefficient": fraction_text(value[1]),
    }


def pair_add(left: Pair, right: Pair) -> Pair:
    return left[0] + right[0], left[1] + right[1]


def pair_scale(value: Pair, scalar: Fraction) -> Pair:
    return value[0] * scalar, value[1] * scalar


def pair_multiply(left: Pair, right: Pair, prime: int) -> Pair:
    return (
        left[0] * right[0] + prime * left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def coefficient_pair(
    prime: int,
    block_length: int,
    index: int,
    q: Fraction,
    s: Fraction,
) -> Pair:
    exponent = index * (2 * block_length - 2 - index)
    scale = q**exponent * s**index
    if index % 2 == 0:
        return scale / prime ** (index // 2), Fraction(0)
    return (
        Fraction(0),
        scale / prime ** ((index + 1) // 2),
    )


def chebyshev_values(maximum: int, x_value: Fraction) -> list[Fraction]:
    values = [Fraction(1)]
    if maximum >= 1:
        values.append(x_value)
    for _ in range(2, maximum + 1):
        values.append(2 * x_value * values[-1] - values[-2])
    return values


def current_pair(
    prime: int,
    block_length: int,
    q: Fraction,
    s: Fraction,
    x_value: Fraction,
) -> Pair:
    coefficients = [
        coefficient_pair(prime, block_length, index, q, s)
        for index in range(block_length)
    ]
    chebyshev = chebyshev_values(block_length - 1, x_value)
    result: Pair = (Fraction(0), Fraction(0))
    for index, coefficient in enumerate(coefficients):
        result = pair_add(
            result,
            pair_scale(
                pair_multiply(coefficient, coefficient, prime),
                Fraction(index + 1),
            ),
        )
    for left in range(block_length):
        for right in range(left + 1, block_length):
            result = pair_add(
                result,
                pair_scale(
                    pair_multiply(
                        coefficients[left],
                        coefficients[right],
                        prime,
                    ),
                    Fraction(left + right + 2)
                    * chebyshev[right - left],
                ),
            )
    return result


def decimal_pair(value: Pair, prime: int) -> str:
    with localcontext() as context:
        context.prec = 45
        root = Decimal(prime).sqrt()
        result = (
            Decimal(value[0].numerator) / Decimal(value[0].denominator)
            + (
                Decimal(value[1].numerator)
                / Decimal(value[1].denominator)
            )
            * root
        )
        return format(result, ".35e")


def negative_guard(value: Pair, prime: int) -> dict:
    rational, radical = value
    if radical == 0:
        if rational >= 0:
            raise RuntimeError("rational witness is not negative")
        return {
            "kind": "negative_rational",
            "positive_margin": fraction_text(-rational),
        }
    if radical < 0 and rational <= 0:
        return {
            "kind": "componentwise_negative",
            "negative_rational": fraction_text(rational),
            "negative_radical": fraction_text(radical),
        }
    if radical < 0 and rational > 0:
        square_margin = prime * radical**2 - rational**2
        if square_margin <= 0:
            raise RuntimeError("radical square guard is not positive")
        return {
            "kind": "negative_by_square",
            "positive_square_margin": fraction_text(square_margin),
            "side_conditions": "rational>0; radical<0",
        }
    if radical > 0 and rational < 0:
        square_margin = rational**2 - prime * radical**2
        if square_margin <= 0:
            raise RuntimeError("reverse radical square guard is not positive")
        return {
            "kind": "negative_by_reverse_square",
            "positive_square_margin": fraction_text(square_margin),
            "side_conditions": "rational<0; radical>0",
        }
    raise RuntimeError("witness is not negative")


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_witnesses() -> list[dict]:
    rows = []
    for prime, block_length, x_value in WITNESSES:
        q_lower, s_lower = BOXES[prime]
        if not q_lower < Q_STAR < 1:
            raise RuntimeError("q witness is outside its source box")
        if not s_lower < Q_STAR < 1:
            raise RuntimeError("s witness is outside its source box")
        if not -1 < x_value < 1:
            raise RuntimeError("angular witness is not interior")
        value = current_pair(
            prime, block_length, Q_STAR, Q_STAR, x_value
        )
        shifted = pair_add(value, (NEGATIVE_TARGET, Fraction(0)))
        guard = negative_guard(shifted, prime)
        boundary_value = current_pair(
            prime,
            block_length,
            Fraction(1),
            Fraction(1),
            x_value,
        )
        negative_guard(boundary_value, prime)
        rows.append(
            {
                "prime": prime,
                "block_length": block_length,
                "q": fraction_text(Q_STAR),
                "s": fraction_text(Q_STAR),
                "x_ang": fraction_text(x_value),
                "current_pair": pair_text(value),
                "current_decimal": decimal_pair(value, prime),
                "target": "J_0<-1/100",
                "shifted_pair": pair_text(shifted),
                "shifted_negative_guard": guard,
                "boundary_q_s_one_pair": pair_text(boundary_value),
                "boundary_q_s_one_decimal": decimal_pair(
                    boundary_value, prime
                ),
                "physical_linkage": (
                    "q=s=9999/10000; "
                    "u=log(p)/2; "
                    "t=-4log(q)/(log(p))^2>0"
                ),
            }
        )
    return rows


def build_rows() -> list[GateRow]:
    return [
        GateRow(
            "pphssf_01_current",
            "exact short-block current",
            "proved",
            "J_0=|Q|^2+Re(H*conj(Q)) with exact heat coefficients.",
            "At rational q,s,x it lies in Q(sqrt(p)).",
            "No division by Q or nonvanishing assumption is used.",
        ),
        GateRow(
            "pphssf_02_domain",
            "interior linked parameter",
            "proved",
            "q=s=9999/10000 is strictly inside both source boxes.",
            "It obeys the physical linkage at u=log(p)/2 and t>0.",
            "This is an ideal block parameter, not an asserted Xi phase hit.",
        ),
        GateRow(
            "pphssf_03_dyadic",
            "all six short dyadic lengths fail",
            "proved",
            "For p=2 and every 2<=M<=7, one rational x has J_0<-1/100.",
            "Uniform short-dyadic block heat-starlikeness is false.",
            "The M=8 positive theorem remains unchanged.",
        ),
        GateRow(
            "pphssf_04_ternary",
            "both short ternary lengths fail",
            "proved",
            "For p=3 and M=2,3, one rational x has J_0<-1/100.",
            "Uniform short-ternary block heat-starlikeness is false.",
            "The M=4 positive theorem remains unchanged.",
        ),
        GateRow(
            "pphssf_05_exact_sign",
            "radical sign certification",
            "proved",
            "Every shifted value J_0+1/100 is exactly negative.",
            "Signs use rational component guards or integer-equivalent square comparisons.",
            "Floating-point values are display-only.",
        ),
        GateRow(
            "pphssf_06_sharpness",
            "small-prime thresholds are method-sharp",
            "proved",
            "All shorter complete lengths fail before the certified M=8 and M=4 bases.",
            "The blockwise positivity thresholds cannot be lowered.",
            "This is sharpness for this current theorem, not for RH.",
        ),
        GateRow(
            "pphssf_07_continuity",
            "open negative neighborhoods",
            "proved",
            "J_0 is a polynomial in q,s,x over Q(sqrt(p)).",
            "Each strict interior witness has a nonempty negative neighborhood.",
            "The obstruction is not an isolated boundary artifact.",
        ),
        GateRow(
            "pphssf_08_route",
            "joined-phase handoff",
            "proved",
            "Short chains must remain inside the endpoint-complete joined current.",
            "Rejoin p-free bases, singletons, endpoints, and cutoffs before a sign estimate.",
            "No coefficient-only sign or blockwise winding sum is promoted.",
        ),
        GateRow(
            "pphssf_09_pi",
            "constant provenance",
            "proved",
            "x_ang=cos(theta) uses the ordinary 2*pi angular period.",
            "No new circle constant enters the counter-witnesses.",
            "Saddle pi remains inherited from completed-zeta normalization.",
        ),
        GateRow(
            "pphssf_10_boundary",
            "proof boundary",
            "open",
            "The impossible universal short-block route is closed.",
            "The joined short-chain, singleton, p-free, endpoint, Abel, and winding theorem remains open.",
            "No contact exclusion, Lambda<=0, PF-infinity, RH, or prize proof is claimed.",
        ),
    ]


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCES.items()
    }


def build_payload() -> dict:
    sources = load_sources()
    witnesses = build_witnesses()
    rows = build_rows()
    return {
        "schema_version": 1,
        "kind": "prime_power_heat_starlikeness_short_family_counter_gate",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "title": "Short small-prime heat-starlikeness counter-gate",
        "source_audit": {
            "source_sha256": {
                key: file_hash(path) for key, path in SOURCES.items()
            },
            "source_kinds": {
                key: payload.get("kind") for key, payload in sources.items()
            },
        },
        "assumptions": {
            "q_star": fraction_text(Q_STAR),
            "s_star": fraction_text(Q_STAR),
            "negative_target": "1/100",
            "dyadic_box": {
                "q": ["47/50", "1"],
                "s": ["22/25", "1"],
            },
            "ternary_box": {
                "q": ["17/20", "1"],
                "s": ["73/100", "1"],
            },
            "physical_linkage": (
                "q=exp[-t(log p)^2/4]; "
                "s=exp[-t(log p)u/2]; "
                "q=s at u=log(p)/2"
            ),
        },
        "rows": [asdict(row) for row in rows],
        "witnesses": witnesses,
        "exact": {
            "current": (
                "J_0=sum_k(k+1)c_k^2+"
                "sum_(k<l)(k+l+2)c_kc_lT_(l-k)(x)"
            ),
            "coefficient": (
                "c_k=p^(-k/2)q^[k(2M-2-k)]s^k"
            ),
            "uniform_counter_margin": (
                "all eight witnesses satisfy J_0<-1/100"
            ),
            "dyadic_failure": "p=2,M=2..7",
            "ternary_failure": "p=3,M=2..3",
            "positive_thresholds_inherited": "p=2,M>=8; p=3,M>=4",
            "continuity": (
                "strict interior polynomial witnesses imply "
                "nonempty relative-open negative neighborhoods"
            ),
        },
        "summary": {
            "rows": len(rows),
            "exact_interior_negative_witnesses": len(witnesses),
            "short_dyadic_lengths_rejected": 6,
            "short_ternary_lengths_rejected": 2,
            "uniform_short_block_positive_families": 0,
            "joined_phase_handoffs": 1,
            "joined_abel_gaps": 0,
            "successor_winding_bounds": 0,
        },
        "proof_boundary": (
            "The exact witnesses reject universal blockwise "
            "heat-starlikeness for p=2,M=2..7 and p=3,M=2..3. "
            "They do not assert that an Xi ray attains a witness phase, "
            "and they do not prove a joined p-free, singleton, endpoint, "
            "cutoff, Abel-gap, or winding theorem, contact exclusion, "
            "Lambda<=0, PF-infinity, RH, or a prize-level result."
        ),
    }


def render_note(payload: dict) -> str:
    lines = [
        "# Short Small-Prime Heat-Starlikeness Counter-Gate",
        "",
        "Date: 2026-07-30",
        "",
        "Status: exact countermodel gate. This is not a proof of a joined",
        "Xi current theorem, `Lambda<=0`, PF-infinity, or RH.",
        "",
        "```text",
        str(RESULT.relative_to(REPO_ROOT)).replace("\\", "/"),
        "python "
        + str(Path(__file__).relative_to(REPO_ROOT)).replace("\\", "/"),
        "python work/rh_compute/scripts/check_"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "prime_power_heat_starlikeness_short_family_counter_gate.py",
        "```",
        "",
        "## Outcome",
        "",
        "The single linked interior point",
        "",
        "```text",
        "q=s=9999/10000",
        "u=log(p)/2",
        "t=-4log(q)/(log(p))^2>0",
        "```",
        "",
        "lies inside both certified small-prime boxes. At the rational",
        "angular coordinates below, exact arithmetic proves `J_0<-1/100`.",
        "",
        "| p | M | x=cos(theta) | exact J_0 at q=s=1 | J_0 at q=s=9999/10000 |",
        "|---:|---:|---:|---|---:|",
    ]
    for row in payload["witnesses"]:
        boundary = row["boundary_q_s_one_pair"]
        pair = (
            f"{boundary['rational']} + "
            f"({boundary['sqrt_prime_coefficient']})sqrt({row['prime']})"
        )
        lines.append(
            f"| {row['prime']} | {row['block_length']} | "
            f"{row['x_ang']} | `{pair}` | "
            f"{row['current_decimal']} |"
        )
    lines.extend(
        [
            "",
            "Every displayed decimal is secondary. The checker reconstructs",
            "the exact pair in `Q(sqrt(p))`, shifts it by `1/100`, and",
            "proves negativity by rational signs or an exact square",
            "comparison.",
            "",
            "## Consequence",
            "",
            "```text",
            "p=2,M=2..7: universal blockwise J_0>0 is false",
            "p=3,M=2..3: universal blockwise J_0>0 is false",
            "```",
            "",
            "Together with the positive base and all-length gates, this",
            "shows that `M=8` for p=2 and `M=4` for p=3 are sharp",
            "thresholds for the present blockwise heat-current theorem.",
            "The short rows are no longer an unfinished Bernstein search:",
            "the desired uniform statement is false.",
            "",
            "The strict witnesses are interior in q, s, and x, so",
            "polynomial continuity gives nonempty negative neighborhoods.",
            "They are also on the exact q-s heat linkage. This still does",
            "not assert that a physical Xi ray attains any listed angular",
            "phase.",
            "",
            "## Route",
            "",
            "Short chains must be retained inside the fully rejoined",
            "endpoint current with their p-free base phases, singleton",
            "chains, recurrent endpoint, and cutoff terms. A sum of",
            "independent blockwise positivity or winding claims cannot",
            "close that joined object.",
            "",
            "The angular variable uses the ordinary `2*pi` period of",
            "`exp(i theta)`. Any `pi` in the saddle scale remains inherited",
            "from completed-zeta normalization.",
            "",
            "## Boundary",
            "",
            payload["proof_boundary"],
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    payload = build_payload()
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    NOTE.write_text(render_note(payload), encoding="utf-8")
    print(
        "built short small-prime heat-starlikeness counter-gate: "
        "8 exact interior witnesses, common margin J_0<-1/100"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
