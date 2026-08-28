#!/usr/bin/env python3
"""Enclose the three nonstationary complementary tails after eight IBP rounds."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_"
    "nonstationary_complement_eight_round_remainder_gate"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

ENDPOINT_RESULT = ROOT / (
    "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_"
    "joined_packet_nonstationary_complement_root_unity_hurwitz_endpoint_gate.json"
)
ENDPOINT_BUILDER = ROOT / (
    "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation4_"
    "joined_packet_nonstationary_complement_root_unity_hurwitz_endpoint_gate.py"
)
ENDPOINT_CHECKER = ROOT / (
    "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation4_"
    "joined_packet_nonstationary_complement_root_unity_hurwitz_endpoint_gate.py"
)

PRECISION_BITS = 384
ORDER = 8
SLABS = 8192
CLUSTER_POWER = 4
HEIGHT = arb("10000000000")
L = arb("621.5")
B = arb("621.5625")
A = arb("39852.5")
Q_PLUS = arb("2561211.5")
Q_REMOTE = Q_PLUS + 240
Q_LOW = arb("79787.5")


CoefficientTable = list[dict[int, dict[tuple[int, int], Fraction]]]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            process.cpu_affinity([process.cpu_affinity()[0]])
            return "below_normal_one_cpu"
        process.nice(10)
        process.cpu_affinity([process.cpu_affinity()[0]])
        return "nice_10_one_cpu"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def arb_from_fraction(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def arb_record(value: arb, digits: int = 72) -> dict[str, Any]:
    return {
        "ball": value.str(digits, more=True),
        "lower_float": float(value.lower()),
        "mid_float": float(value.mid()),
        "upper_float": float(value.upper()),
    }


def acb_record(value: acb, digits: int = 72) -> dict[str, str]:
    return {
        "real_ball": value.real.str(digits, more=True),
        "imag_ball": value.imag.str(digits, more=True),
        "absolute_ball": abs(value).str(digits, more=True),
    }


def saved_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def add_complex_error(value: acb, radius: arb) -> acb:
    return acb(arb(value.real, radius), arb(value.imag, radius))


def q_curve(y: arb, p: arb) -> arb:
    return y + p / y


def coefficient_table(order: int) -> CoefficientTable:
    """Exact monomials for h_(n+1)=(h_n/D)' with D=q-y-p/y."""

    require(order >= 1, "positive recurrence order required")
    table: CoefficientTable = [{0: {(0, -1): Fraction(1)}}]
    for _ in range(order):
        next_level: dict[int, dict[tuple[int, int], Fraction]] = {}
        for denominator_power, terms in table[-1].items():
            for (p_power, exponent2), coefficient in terms.items():
                updates = (
                    (
                        denominator_power + 1,
                        (p_power, exponent2 - 2),
                        coefficient * Fraction(exponent2, 2),
                    ),
                    (
                        denominator_power + 2,
                        (p_power, exponent2),
                        coefficient * (denominator_power + 1),
                    ),
                    (
                        denominator_power + 2,
                        (p_power + 1, exponent2 - 4),
                        -coefficient * (denominator_power + 1),
                    ),
                )
                for new_power, monomial, value in updates:
                    destination = next_level.setdefault(new_power, {})
                    destination[monomial] = destination.get(monomial, Fraction(0)) + value
        table.append(
            {
                power: {monomial: value for monomial, value in terms.items() if value}
                for power, terms in next_level.items()
            }
        )
    return table


def half_power(base: arb, exponent2: int) -> arb:
    return base.sqrt() ** exponent2


def evaluate_coefficient(
    terms: dict[tuple[int, int], Fraction], p: arb, y: arb
) -> arb:
    value = arb(0)
    for (p_power, exponent2), coefficient in terms.items():
        value += (
            arb_from_fraction(coefficient)
            * p**p_power
            * half_power(y, exponent2)
        )
    return value


def coefficient_absolute_upper(
    terms: dict[tuple[int, int], Fraction], p: arb, y_left: arb
) -> arb:
    value = arb(0)
    for (p_power, exponent2), coefficient in terms.items():
        require(exponent2 < 0, "coefficient monotonicity guard failed")
        value += (
            abs(arb_from_fraction(coefficient))
            * p**p_power
            * half_power(y_left, exponent2)
        )
    return value.upper()


def endpoint_carrier(t: arb, y: arb) -> acb:
    return (acb(0, -t * y.log() - arb.pi() * y * y)).exp()


def endpoint_level_sum(
    level: dict[int, dict[tuple[int, int], Fraction]],
    p: arb,
    t: arb,
    endpoint: dict[str, Any],
) -> acb:
    y = arb(endpoint["endpoint_y"])
    rows = {row["power"]: row for row in endpoint["rows"]}
    value = acb(0)
    for denominator_power, terms in level.items():
        required_power = denominator_power + 1
        require(required_power in rows, f"missing endpoint denominator power {required_power}")
        coefficient = evaluate_coefficient(terms, p, y)
        endpoint_sum = saved_complex(rows[required_power]["phase_weighted_sum_ball"])
        value += coefficient * endpoint_sum
    return endpoint_carrier(t, y) * value


def endpoint_expansion(
    table: CoefficientTable,
    p: arb,
    t: arb,
    left_endpoint: dict[str, Any],
    right_endpoint: dict[str, Any],
    order: int,
) -> tuple[acb, list[dict[str, Any]]]:
    partial = acb(0)
    terms: list[dict[str, Any]] = []
    imaginary_scale = acb(0, 1) * 2 * arb.pi()
    for n in range(order):
        left_value = endpoint_level_sum(table[n], p, t, left_endpoint)
        right_value = endpoint_level_sum(table[n], p, t, right_endpoint)
        factor = (-1) ** n * imaginary_scale ** (-(n + 1))
        term = factor * (right_value - left_value)
        partial += term
        terms.append(
            {
                "round": n,
                "left_boundary_ball": acb_record(left_value, 60),
                "right_boundary_ball": acb_record(right_value, 60),
                "signed_term_ball": acb_record(term, 60),
            }
        )
    return partial, terms


def clustered_edge(
    left: arb,
    right: arb,
    index: int,
    slabs: int,
    side: str,
) -> arb:
    denominator = slabs**CLUSTER_POWER
    width = right - left
    if side == "left":
        return left + width * arb(index**CLUSTER_POWER) / denominator
    if side == "right":
        return right - width * arb((slabs - index) ** CLUSTER_POWER) / denominator
    raise RuntimeError(f"unknown cluster side: {side}")


def infinite_denominator_bound(delta: arb, power: int) -> arb:
    require(delta.lower() > 0 and power > 1, "invalid infinite denominator bound")
    return (delta ** (-power) + delta ** (1 - power) / (power - 1)).upper()


def denominator_bound(delta: arb, power: int, finite_count: int | None) -> arb:
    infinite = infinite_denominator_bound(delta, power)
    if finite_count is None:
        return infinite
    finite = (arb(finite_count) * delta ** (-power)).upper()
    return finite if finite < infinite else infinite


def integrated_remainder_bound(
    level: dict[int, dict[tuple[int, int], Fraction]],
    p: arb,
    left: arb,
    right: arb,
    q0: arb,
    orientation: str,
    finite_count: int | None,
    slabs: int,
) -> tuple[arb, dict[str, Any]]:
    require(right < p.sqrt(), "production interval crossed the Q minimum")
    cluster_side = "left" if orientation == "upper" else "right"
    integral = arb(0)
    minimum_gap: arb | None = None
    for index in range(slabs):
        slab_left = clustered_edge(left, right, index, slabs, cluster_side)
        slab_right = clustered_edge(left, right, index + 1, slabs, cluster_side)
        width = slab_right - slab_left
        if orientation == "upper":
            delta = (q0 - q_curve(slab_left, p)).lower()
        elif orientation == "lower":
            delta = (q_curve(slab_right, p) - q0).lower()
        else:
            raise RuntimeError(f"unknown tail orientation: {orientation}")
        require(delta > 0, f"nonstationary gap lost in slab {index}")
        if minimum_gap is None or delta < minimum_gap:
            minimum_gap = delta
        slab_majorant = arb(0)
        for power, terms in level.items():
            coefficient = coefficient_absolute_upper(terms, p, slab_left)
            slab_majorant += coefficient * denominator_bound(delta, power, finite_count)
        integral += width * slab_majorant
    require(minimum_gap is not None and integral.is_finite(), "remainder census failed")
    scaled = (integral / (2 * arb.pi()) ** ORDER).upper()
    return scaled, {
        "slabs": slabs,
        "cluster_side": cluster_side,
        "cluster_power": CLUSTER_POWER,
        "minimum_denominator_gap_ball": arb_record(minimum_gap),
        "unscaled_integral_majorant_ball": arb_record(integral.upper()),
        "scaled_remainder_bound_ball": arb_record(scaled),
    }


def family_certificate(
    *,
    name: str,
    table: CoefficientTable,
    p: arb,
    endpoints: dict[str, Any],
    left_key: str,
    right_key: str,
    left: arb,
    right: arb,
    q0: arb,
    orientation: str,
    finite_count: int | None,
) -> dict[str, Any]:
    partial, endpoint_terms = endpoint_expansion(
        table, p, HEIGHT, endpoints[left_key], endpoints[right_key], ORDER
    )
    remainder, remainder_census = integrated_remainder_bound(
        table[ORDER], p, left, right, q0, orientation, finite_count, SLABS
    )
    enclosure = add_complex_error(partial, remainder)
    return {
        "name": name,
        "orientation": orientation,
        "finite_label_count": finite_count,
        "interval": [left.str(30, more=True), right.str(30, more=True)],
        "integration_rounds": ORDER,
        "endpoint_terms": endpoint_terms,
        "endpoint_partial_sum_ball": acb_record(partial),
        "remainder_census": remainder_census,
        "tail_integral_ball": acb_record(enclosure),
    }


def build_certificate(endpoint_artifact: dict[str, Any]) -> dict[str, Any]:
    p = HEIGHT / (2 * arb.pi())
    table = coefficient_table(ORDER)
    families = endpoint_artifact["certificate"]["families"]

    def scaled_detuning(delta: arb, y: arb) -> arb:
        q_derivative = 1 - p / y**2
        require(not q_derivative.contains(0), "endpoint quadratic scale vanished")
        return delta / abs(q_derivative).sqrt()

    rb = family_certificate(
        name="R_B",
        table=table,
        p=p,
        endpoints=families["R_B_finite_240"],
        left_key="endpoint_B",
        right_key="endpoint_a",
        left=B,
        right=A,
        q0=Q_PLUS,
        orientation="upper",
        finite_count=240,
    )
    rl = family_certificate(
        name="R_L",
        table=table,
        p=p,
        endpoints=families["R_L_infinite_upper"],
        left_key="endpoint_L",
        right_key="endpoint_a",
        left=L,
        right=A,
        q0=Q_REMOTE,
        orientation="upper",
        finite_count=None,
    )
    lower = family_certificate(
        name="R_lower",
        table=table,
        p=p,
        endpoints=families["R_lower_infinite"],
        left_key="endpoint_L",
        right_key="endpoint_a",
        left=L,
        right=A,
        q0=Q_LOW,
        orientation="lower",
        finite_count=None,
    )

    return {
        "passed": True,
        "precision_bits": PRECISION_BITS,
        "integration_rounds": ORDER,
        "remainder_slabs_per_family": SLABS,
        "cluster_power": CLUSTER_POWER,
        "coefficient_level_term_counts": [
            {
                "level": level,
                "denominator_powers": sorted(table[level]),
                "exact_monomial_count": sum(len(terms) for terms in table[level].values()),
            }
            for level in range(ORDER + 1)
        ],
        "recurrence": "h_0=y^(-1/2), h_(n+1)=d/dy[h_n/(q-y-p/y)]",
        "ibp_identity": (
            "I=sum_(n=0)^(N-1)(-1)^n(i2*pi)^(-n-1)"
            "[exp(i*phi)h_n/D]_left^right+(-1)^N(i2*pi)^(-N)"
            "integral exp(i*phi)h_N"
        ),
        "families": {"R_B": rb, "R_L": rl, "R_lower": lower},
        "endpoint_scaled_detuning": {
            "definition": "chi=delta/sqrt(abs(Q'(y))), Q'(y)=1-p/y^2",
            "R_B_at_B_ball": arb_record(scaled_detuning(Q_PLUS - q_curve(B, p), B)),
            "R_L_at_L_ball": arb_record(scaled_detuning(Q_REMOTE - q_curve(L, p), L)),
            "R_lower_at_a_ball": arb_record(scaled_detuning(q_curve(A, p) - Q_LOW, A)),
            "interpretation": (
                "R_B and R_L begin inside an exterior endpoint-transition scale; "
                "ordinary repeated integration by parts is asymptotic there"
            ),
        },
        "abel_limit_justification": (
            "derive at finite symmetric label cutoff; endpoint k=1 sums converge in the "
            "certified root-of-unity Abel convention, while every N=8 integrated "
            "remainder label sum is absolutely dominated by powers k>=8"
        ),
    }


def note_text(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    rb = c["families"]["R_B"]
    rl = c["families"]["R_L"]
    lower = c["families"]["R_lower"]
    return f"""# Eight-round enclosures of the nonstationary complementary tails

Date: 2026-08-27

Status: `R_B`, `R_L`, and the lower complementary tail are each enclosed as
rigorous complex balls after eight integration-by-parts rounds, but the balls
are numerically useless at the `D_K` scale.  This certifies rejection of the
unsplit absolute eight-round route.  This is not a proof of the complete
packet, `D_K`, or RH.

For each label put

```text
D_q(y)=q-y-p/y,             phi_q'(y)=2*pi*D_q(y),
h_0(y)=y^(-1/2),            h_(n+1)=d/dy[h_n/D_q].  (ER1)
```

Then the exact finite-cutoff identity is

```text
I_q=sum_(n=0)^7 (-1)^n(i2*pi)^(-n-1)
    [exp(i*phi_q)h_n/D_q]_left^right
    +(i2*pi)^(-8) integral exp(i*phi_q)h_8.         (ER2)
```

Every `h_n` is generated as an exact rational monomial table in
`p`, half-integer powers of `y`, and powers of `D_q^(-1)`.  The level-eight
table has {c['coefficient_level_term_counts'][8]['exact_monomial_count']}
exact monomials over denominator powers 8 through 16.  Thus its summed
integral is absolutely convergent even though the first endpoint sum uses the
root-of-unity Abel convention.

The endpoint terms use the certified Section 11.500 balls before any norm.
For the integrated remainder, 8192 rational slabs are fourth-power clustered
at the smallest-gap endpoint.  On each slab, negative powers of `y` are
maximized at its left edge and the denominator roster is bounded by

```text
sum_(m>=0)(delta+m)^(-k)
 <=delta^(-k)+delta^(1-k)/(k-1),                   (ER3)
```

with the better finite-count bound used for `R_B`.  This gives

```text
R_B endpoint partial = {rb['endpoint_partial_sum_ball']['real_ball']}
                     + i {rb['endpoint_partial_sum_ball']['imag_ball']},
R_B remainder <= {rb['remainder_census']['scaled_remainder_bound_ball']['ball']},
R_B ball = {rb['tail_integral_ball']['real_ball']}
         + i {rb['tail_integral_ball']['imag_ball']};

R_L endpoint partial = {rl['endpoint_partial_sum_ball']['real_ball']}
                     + i {rl['endpoint_partial_sum_ball']['imag_ball']},
R_L remainder <= {rl['remainder_census']['scaled_remainder_bound_ball']['ball']},
R_L ball = {rl['tail_integral_ball']['real_ball']}
         + i {rl['tail_integral_ball']['imag_ball']};

R_lower endpoint partial = {lower['endpoint_partial_sum_ball']['real_ball']}
                          + i {lower['endpoint_partial_sum_ball']['imag_ball']},
R_lower remainder <= {lower['remainder_census']['scaled_remainder_bound_ball']['ball']},
R_lower ball = {lower['tail_integral_ball']['real_ball']}
             + i {lower['tail_integral_ball']['imag_ball']}.
```

The obstruction is structural, not a failed sign check.  With
`chi=delta/sqrt(abs(Q'(y)))`, the three endpoint detunings are

```text
chi_RB(B)    = {c['endpoint_scaled_detuning']['R_B_at_B_ball']['ball']},
chi_RL(L)    = {c['endpoint_scaled_detuning']['R_L_at_L_ball']['ball']},
chi_lower(a) = {c['endpoint_scaled_detuning']['R_lower_at_a_ball']['ball']}.
```

Thus `R_B` and especially `R_L` begin within the exterior quadratic
endpoint-transition scale.  Their enormous high-order endpoint terms are the
expected asymptotic failure of ordinary integration by parts near an excluded
saddle.  They require a grouped Fresnel/erfc endpoint layer before the remote
tail is integrated by parts.  The lower tail has a different geometry and
should be sharpened separately rather than hidden inside the same failed bound.

Pi provenance: every `pi` comes from the inherited Fresnel/Fourier phase,
the derivative `phi'=2*pi*D`, or `p=t/(2*pi)`.  No fitted circle constant is
used.

Proof boundary: rigorous complex enclosures of the three isolated
nonstationary complementary tails at `t=10^10` only.  They have not yet been
joined with `P_240`, the lower cell, transition packet, upper arc, natural A
lift, or Gamma defect.  No complete ordinary or joined packet, `J_Z`, `D_K`,
non-A, all-height, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is
proved.
"""


def main() -> int:
    resource_mode = set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = PRECISION_BITS
    ctx.threads = 1
    dependencies = {
        "endpoint_result": ENDPOINT_RESULT,
        "endpoint_builder": ENDPOINT_BUILDER,
        "endpoint_checker": ENDPOINT_CHECKER,
    }
    for path in (*dependencies.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    endpoint_artifact = load_json(ENDPOINT_RESULT)
    require(endpoint_artifact.get("passed") is True, "endpoint gate not passed")
    require(
        endpoint_artifact["certificate"]["maximum_denominator_power"] >= 2 * ORDER,
        "endpoint power table is too short",
    )
    require(
        endpoint_artifact["decision"]["integrated_remainders_bounded"] is False,
        "endpoint proof boundary drift",
    )

    certificate = build_certificate(endpoint_artifact)
    artifact = {
        "kind": "rh_c1_hardy_equation4_nonstationary_complement_eight_round_remainder_gate",
        "date": "2026-08-27",
        "status": "three_tail_balls_certified_unsplit_eight_round_absolute_route_rejected",
        "passed": True,
        "resource_mode": resource_mode,
        "certificate": certificate,
        "decision": {
            "exact_eight_round_coefficient_recurrence_certified": True,
            "R_B_complex_ball_certified": True,
            "R_L_complex_ball_certified": True,
            "lower_complementary_tail_complex_ball_certified": True,
            "endpoint_carriers_retained_before_norm": True,
            "eight_round_absolute_bounds_useful_at_D_K_scale": False,
            "unsplit_eight_round_absolute_route_rejected": True,
            "complete_ordinary_packet_enclosed": False,
            "complete_joined_packet_enclosed": False,
            "D_K_corridor_certified": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in dependencies.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": (
            "Replace the near-endpoint portions of R_B and R_L by a grouped uniform "
            "Fresnel/erfc exterior-saddle layer, then apply integration by parts only after "
            "a certified scaled-detuning cutoff. Sharpen the lower complementary tail in its "
            "own geometry. Join pieces only after each new remainder closes at the D_K scale."
        ),
        "proof_boundary": (
            "Three isolated but D_K-scale-useless production tail balls and a certified route "
            "rejection only. No complete ordinary or joined packet, "
            "J_Z, D_K corridor, non-A, all-height, Lambda<=0, PF-infinity, RH, or prize-level "
            "conclusion is proved."
        ),
    }
    atomic_write(RESULT, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE, note_text(artifact))
    print("certified eight-round complementary-tail enclosures", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
