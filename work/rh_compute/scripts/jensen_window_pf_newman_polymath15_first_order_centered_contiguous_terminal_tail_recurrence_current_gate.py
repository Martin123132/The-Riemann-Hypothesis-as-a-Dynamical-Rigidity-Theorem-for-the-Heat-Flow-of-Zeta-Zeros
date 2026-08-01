#!/usr/bin/env python3
"""Build the contiguous terminal-tail recurrence/current theorem."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "contiguous_terminal_tail_recurrence_current_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
PRECISION_BITS = 192
SINC_TERMS = 30
SINC_TAIL_RADIUS = "1e-50"
E_BOXES = 16
Z_BOXES = 64
SOURCE_PATHS = {
    "endpoint_normalization": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "complex_endpoint_source_normalization_gate.json"
    ),
    "adjacent_recurrence": (
        REPO_ROOT
        / "work/rh_compute/scripts/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "adjacent_saddle_recurrence.py"
    ),
    "edge_current": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "real_edge_projective_current_gate.json"
    ),
    "pair_guard": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "edge_anchored_near_terminal_pair_current_guard.json"
    ),
    "cumulative_handoff": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "residual_placement_c1_cumulative_handoff_gate.json"
    ),
}


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def aq(value: Fraction | int) -> flint.arb:
    value = Fraction(value)
    return flint.arb(flint.fmpq(value.numerator, value.denominator))


def interval_ball(left: Fraction, right: Fraction) -> flint.arb:
    lower = aq(left)
    upper = aq(right)
    return (lower + upper) / 2 + flint.arb(0, (upper - lower) / 2)


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def load_sources() -> dict[str, dict | str]:
    missing = [str(path) for path in SOURCE_PATHS.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    payloads: dict[str, dict | str] = {}
    for key, path in SOURCE_PATHS.items():
        payloads[key] = (
            json.loads(path.read_text(encoding="utf-8"))
            if path.suffix == ".json"
            else path.read_text(encoding="utf-8")
        )
    return payloads


def source_audit(payloads: dict[str, dict | str]) -> dict:
    endpoint = payloads["endpoint_normalization"]
    if not isinstance(endpoint, dict) or "C_0(p)={exp" not in endpoint.get(
        "exact", {}
    ).get("source_normalization", {}).get("c0_formula", ""):
        raise RuntimeError("complex endpoint normalization drifted")

    recurrence = payloads["adjacent_recurrence"]
    if not isinstance(recurrence, str) or "C0(p+2)+C0(p)=exp" not in recurrence:
        raise RuntimeError("C0 recurrence source drifted")

    edge = payloads["edge_current"]
    if not isinstance(edge, dict) or "K_edge<-3/8000" not in edge.get(
        "interval", {}
    ).get("current_symbol_margin", ""):
        raise RuntimeError("edge current margin drifted")

    pair = payloads["pair_guard"]
    if not isinstance(pair, dict) or "P_3(-5/8)>1/10" not in pair.get(
        "interval_certificate", {}
    ).get("positive_margin", ""):
        raise RuntimeError("near-terminal pair guard drifted")

    handoff = payloads["cumulative_handoff"]
    if not isinstance(handoff, dict) or "K(V_J)=-R_diag+C_cross" not in handoff.get(
        "exact", {}
    ).get("diagonal_reserve", ""):
        raise RuntimeError("cumulative handoff drifted")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCE_PATHS.items()
        },
        "published_source": {
            "citation": (
                "D.H.J. Polymath, Effective approximation of heat flow "
                "evolution of the Riemann xi function, and a new upper "
                "bound for the de Bruijn-Newman constant, arXiv:1904.12438"
            ),
            "location": "equation (53), C_0, and the adjacent endpoint recurrence",
            "url": "https://arxiv.org/abs/1904.12438",
        },
    }


def symbolic_certificate() -> dict:
    p, k, x, y = sp.symbols("p k X Y", real=True)
    pi = sp.pi
    z = x + sp.I * y
    z_p = sp.I * 2 * pi * k * z
    z_pp = -sp.I * pi * z - 4 * pi**2 * k**2 * z
    d = k * y / 2
    n = y / (16 * pi) - k**2 * x / 4
    if sp.simplify(sp.re(z_p) + 4 * pi * d) != 0:
        raise RuntimeError("first carrier derivative failed")
    if sp.simplify(sp.re(z_pp) - 16 * pi**2 * n) != 0:
        raise RuntimeError("second carrier derivative failed")

    c, cp, cpp = sp.symbols("C C_p C_pp", real=True)
    current = c * cpp / (16 * pi**2) - (-cp / (4 * pi)) ** 2
    target = (c * cpp - cp**2) / (16 * pi**2)
    if sp.simplify(current - target) != 0:
        raise RuntimeError("grouped current identity failed")

    m = sp.symbols("m", integer=True, nonnegative=True)
    phi = p**2 / 2 - p + sp.Rational(3, 8)
    shifted_phase = p**2 / 2 - (2 * m + 1) * p + sp.Rational(3, 8)
    exponent_gap = sp.expand(-phi + 2 * m * p + shifted_phase)
    if exponent_gap != 0:
        raise RuntimeError("shifted recurrence phase failed")

    return {
        "endpoint_recurrence": "C_0(p+2)+C_0(p)=R(p)",
        "iterated_recurrence": (
            "C_0(p)=sum_(m=0)^M(-1)^m R(p-2m-2)"
            "+(-1)^(M+1)C_0(p-2M-2)"
        ),
        "phase_match": (
            "Re[Q(p)r(p)^m]=(-1)^m Re[R(p-2m-2)], "
            "Q=exp[-pi*i*(p^2/2-p+3/8)], r=-exp(2*pi*i*p)"
        ),
        "grouped_trace": (
            "C_M(p)=Re[C_0(p)-Q*sum_(m=0)^M r^m]="
            "(-1)^(M+1)Re C_0(p-2M-2)"
        ),
        "grouped_jets": (
            "D_M=-C_M'/(4*pi), N_M=C_M''/(16*pi^2)"
        ),
        "grouped_current": (
            "P_M=[C_M*C_M''-(C_M')^2]/(16*pi^2)"
        ),
    }


def coefficient_lists() -> tuple[list[Fraction], list[Fraction], list[Fraction]]:
    c0 = [Fraction((-1) ** n, math.factorial(2 * n + 1)) for n in range(SINC_TERMS)]
    c1 = [
        Fraction((-1) ** n, (2 * n + 1) * math.factorial(2 * n - 1))
        for n in range(1, SINC_TERMS)
    ]
    c2 = [
        Fraction((-1) ** n, (2 * n + 1) * math.factorial(2 * n - 2))
        for n in range(1, SINC_TERMS)
    ]
    return c0, c1, c2


SINC_COEFFICIENTS = coefficient_lists()


def verify_sinc_tail_bound() -> dict:
    bounds = {}
    for derivative in range(3):
        n = SINC_TERMS
        first = Fraction(
            2 ** (2 * n - derivative),
            (2 * n + 1) * math.factorial(2 * n - derivative),
        )
        exact_ratio = Fraction(
            4 * (2 * n + 1),
            (2 * n + 3) * (2 * n + 2 - derivative) * (2 * n + 1 - derivative),
        )
        ratio = Fraction(1, 800)
        if not exact_ratio < ratio:
            raise RuntimeError("sinc successive-tail ratio bound failed")
        tail = first / (1 - ratio)
        if not tail < Fraction(1, 10**50):
            raise RuntimeError("sinc derivative tail bound failed")
        bounds[str(derivative)] = f"{tail.numerator}/{tail.denominator}<1e-50"
    return {
        "argument_bound": (
            "|pi*e*z|<=pi/5<2 and "
            "|pi*z*(1+e^2*z/2)|<=13*pi/25<2"
        ),
        "terms": SINC_TERMS,
        "ratio_bound": "successive absolute tail ratio <1/800",
        "derivative_tail_bounds": bounds,
    }


def polynomial_horner(value: flint.arb, coefficients: list[Fraction]) -> flint.arb:
    result = flint.arb(0)
    for coefficient in reversed(coefficients):
        result = result * value + aq(coefficient)
    return result


def sinc_jet(value: flint.arb) -> tuple[flint.arb, flint.arb, flint.arb]:
    square = value * value
    c0, c1, c2 = SINC_COEFFICIENTS
    error = flint.arb(0, SINC_TAIL_RADIUS)
    return (
        polynomial_horner(square, c0) + error,
        value * polynomial_horner(square, c1) + error,
        polynomial_horner(square, c2) + error,
    )


def compact_curvature(e: flint.arb, z: flint.arb) -> flint.arb:
    pi = flint.arb.pi()
    e2 = e * e
    weight = 1 + e2 * z / 2
    numerator_argument = pi * z * weight
    denominator_argument = pi * e * z
    f, f1, f2 = sinc_jet(numerator_argument)
    g, g1, g2 = sinc_jet(denominator_argument)

    weight_z = e2 / 2
    numerator_z = pi * (1 + e2 * z)
    numerator_zz = pi * e2
    denominator_z = pi * e

    f_z = f1 * numerator_z
    f_zz = f2 * numerator_z**2 + f1 * numerator_zz
    g_z = g1 * denominator_z
    g_zz = g2 * denominator_z**2
    inv_g = 1 / g
    inv_g_z = -g_z / g**2
    inv_g_zz = 2 * g_z**2 / g**3 - g_zz / g**2

    q = weight * f * inv_g
    q_z = (
        weight_z * f * inv_g
        + weight * f_z * inv_g
        + weight * f * inv_g_z
    )
    q_zz = (
        2 * weight_z * f_z * inv_g
        + 2 * weight_z * f * inv_g_z
        + weight * f_zz * inv_g
        + 2 * weight * f_z * inv_g_z
        + weight * f * inv_g_zz
    )
    return q_z**2 - q * q_zz


def compact_box_certificate() -> dict:
    flint.ctx.prec = PRECISION_BITS
    tail = verify_sinc_tail_bound()
    minimum: flint.arb | None = None
    location = ""
    rows = []
    for i in range(E_BOXES):
        e_left = Fraction(2 * i, 5 * E_BOXES)
        e_right = Fraction(2 * (i + 1), 5 * E_BOXES)
        e = interval_ball(e_left, e_right)
        strip_minimum: flint.arb | None = None
        for j in range(Z_BOXES):
            z_left = Fraction(-1, 2) + Fraction(j, Z_BOXES)
            z_right = Fraction(-1, 2) + Fraction(j + 1, Z_BOXES)
            z = interval_ball(z_left, z_right)
            curvature = compact_curvature(e, z)
            lower = curvature.lower()
            if not lower > 0:
                raise RuntimeError(
                    f"unresolved compact box e={e_left}:{e_right}, "
                    f"z={z_left}:{z_right}: {curvature}"
                )
            if strip_minimum is None or lower < strip_minimum:
                strip_minimum = lower
            if minimum is None or lower < minimum:
                minimum = lower
                location = f"e=[{e_left},{e_right}], z=[{z_left},{z_right}]"
        assert strip_minimum is not None
        rows.append(
            {
                "e_interval": [str(e_left), str(e_right)],
                "minimum_lower": strip_minimum.str(40),
            }
        )
    assert minimum is not None
    if not minimum > aq(Fraction(7, 5)):
        raise RuntimeError("compact curvature margin 7/5 failed")
    return {
        "precision_bits": PRECISION_BITS,
        "domain": "0<=e<=2/5, -1/2<=z<=1/2",
        "boxes": E_BOXES * Z_BOXES,
        "e_boxes": E_BOXES,
        "z_boxes": Z_BOXES,
        "normalized_curvature": "q_z^2-q*q_zz",
        "strict_margin": "q_z^2-q*q_zz>7/5",
        "minimum_lower": minimum.str(70),
        "minimum_box": location,
        "sinc_tail": tail,
        "e_strip_minima": rows,
    }


def exact_payload() -> dict:
    return {
        "cofinal_path": (
            "Fix p in [-1,1], theta=(1-p)/2, q=2tL^2=1, "
            "a=N+theta, and M>=0; take N->infinity with M fixed."
        ),
        "tail_definition": (
            "Group the retained endpoint-terminal edge with the contiguous "
            "carriers n=N-1,...,N-M."
        ),
        "real_c0": (
            "f(t)=Re C_0(t)=cos(pi*(t^2/2+3/8))/(2*cos(pi*t)), "
            "with analytic continuation at half-integers"
        ),
        "trace_shift": (
            "For K=M+1, C_M(p)=(-1)^K*f(p-2K); since f is even, "
            "its curvature is evaluated at t=2K-p>=1."
        ),
        "curvature": (
            "D(t)=f'(t)^2-f(t)f''(t), P_M(p)=-D(2M+2-p)/(16*pi^2)."
        ),
        "base_cell": (
            "For M=0, t in [1,3], the imported 3072-box edge theorem "
            "gives D(t)>3/50."
        ),
        "away_bound": (
            "For t>=3 and |cos(pi*t)|>6/(5t), direct phase minimization "
            "gives 4*cos(pi*t)^4*D >= "
            "pi^2*t^2*cos(pi*t)^2-pi^2/2-"
            "(pi/2)*sqrt(pi^2+cos(pi*t)^4). Using pi>3, pi<22/7, "
            "and t^2*cos(pi*t)^2>36/25 proves D>4/5."
        ),
        "near_reduction": (
            "Otherwise choose the nearest alpha in Z+1/2, write y=t-alpha, "
            "e=1/alpha, z=alpha*y. Then 0<e<=2/5 and |z|<1/2. With "
            "q_e(z)=(1+e^2*z/2)*sinc(pi*z*(1+e^2*z/2))/sinc(pi*e*z), "
            "f(alpha+y)=+/-q_e(z)/(2e) and "
            "D=(q_z^2-q*q_zz)/(4e^4)."
        ),
        "near_domain_bound": (
            "arcsin(x)<=x/sqrt(1-x^2), pi>3, and t>=3 give "
            "|z|<=2/sqrt(21)+4/189<1/2."
        ),
        "theorem": (
            "For every integer M>=0 and every p in [-1,1], D(2M+2-p)>3/50 "
            "and the fixed-M contiguous terminal-tail leading symbol "
            "satisfies P_M(p)<-3/8000."
        ),
        "endpoint_diagnostic": (
            "P_M(1)=P_0(1)-M*(M+1)/16."
        ),
        "pi_provenance": (
            "Every pi is inherited from the completed-zeta/Riemann-Siegel "
            "normalization and C_0 recurrence. No circle, polygon, fitted "
            "period, or plotted symmetry introduces a new constant."
        ),
        "surviving_target": (
            "Prove a finite-height estimate uniform for a growing terminal "
            "prefix and join it to the nonterminal Abel/cross-current bulk, "
            "or prove the endpoint-complete division-free Abel gap directly."
        ),
    }


def build_rows(exact: dict, symbolic: dict, interval: dict) -> list[GateRow]:
    return [
        GateRow("cttr_01_path", "cofinal path", "proved", "The tail limit stays on the physical q=1 boundary.", exact["cofinal_path"], "M is fixed before N tends to infinity."),
        GateRow("cttr_02_tail", "canonical grouping", "proved", "The tested block is a contiguous terminal prefix.", exact["tail_definition"], "Arbitrary sparse carrier subsets are not covered."),
        GateRow("cttr_03_recurrence", "exact recurrence", "proved", "The C_0 recurrence iterates without absolute values.", symbolic["iterated_recurrence"], "This is an algebraic induction in M."),
        GateRow("cttr_04_phase", "phase matching", "proved", "Every limiting carrier phase matches one shifted recurrence term after real projection.", symbolic["phase_match"], "The equality is for real traces; the complex terms are conjugate."),
        GateRow("cttr_05_trace", "tail collapse", "proved", "The complete finite terminal prefix collapses to one shifted C_0 trace.", symbolic["grouped_trace"], "No termwise current sign is asserted."),
        GateRow("cttr_06_jets", "differentiated identity", "proved", "The grouped first and second jets are derivatives of the collapsed trace.", symbolic["grouped_jets"], "Analytic continuation covers removable half-integers."),
        GateRow("cttr_07_current", "projective current", "proved", "The grouped current is the negative Laguerre curvature of the real trace.", symbolic["grouped_current"], "The division-free polynomial is primary at trace zeros."),
        GateRow("cttr_08_base", "base cell", "interval proved", "The existing edge theorem supplies the first shifted cell.", exact["base_cell"], "This imports its 192-bit removable-chart certificate."),
        GateRow("cttr_09_away", "analytic tail", "proved", "Direct phase minimization gives a quantitative margin away from denominator removals.", exact["away_bound"], "The bound is division-free and includes trace zeros."),
        GateRow("cttr_10_near", "compact reduction", "proved", "All remaining half-integer neighborhoods reduce to one compact sinc family.", exact["near_reduction"], "The sign +/- cancels from the quadratic curvature."),
        GateRow("cttr_11_domain", "domain bound", "proved", "The physical near-pole coordinates lie strictly inside the certified compact box.", exact["near_domain_bound"], "The interval proof deliberately certifies a larger closed box."),
        GateRow("cttr_12_compact", "Arb certificate", "interval proved", "The normalized compact curvature has a strict rational margin.", interval["strict_margin"], "The sinc derivative tails are enclosed explicitly.", interval),
        GateRow("cttr_13_theorem", "asymptotic theorem", "proved", "Every fixed finite contiguous terminal prefix preserves the original strict edge margin.", exact["theorem"], "This does not permit M to grow with N."),
        GateRow("cttr_14_endpoint", "exact diagnostic", "proved", "The endpoint pattern seen in the scout follows exactly.", exact["endpoint_diagnostic"], "This diagnostic is not used as a substitute for the full-cell proof."),
        GateRow("cttr_15_pi", "constant provenance", "proved", "No unexplained pi is inserted.", exact["pi_provenance"], "The normalization is source-derived."),
        GateRow("cttr_16_handoff", "open handoff", "open", "The result identifies the next aggregate interface rather than closing it.", exact["surviving_target"], "No growing-tail estimate or full Xi cross-current theorem is claimed."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    symbolic = payload["symbolic_certificate"]
    interval = payload["interval_certificate"]
    return f"""# Contiguous Terminal-Tail Recurrence and Current Gate

Date: 2026-07-31

Status: exact recurrence and Arb-uniform fixed-M leading current theorem. This
is not a proof of a growing-tail estimate, complete aggregate bound,
`Lambda<=0`, or RH.

## Exact Collapse

{exact['cofinal_path']}

{symbolic['iterated_recurrence']}

{symbolic['phase_match']}

Therefore

```text
{symbolic['grouped_trace']}
{symbolic['grouped_jets']}
{symbolic['grouped_current']}
```

This explains why the isolated `N-3` carrier can reverse the edge current
while the canonical edge plus `N-1,N-2,N-3` prefix remains clockwise.

## Uniform Curvature

{exact['real_c0']}

{exact['base_cell']}

{exact['away_bound']}

{exact['near_reduction']}

{exact['near_domain_bound']}

A {interval['precision_bits']}-bit Arb cover of {interval['boxes']} rational
boxes proves

```text
{interval['strict_margin']}
minimum lower enclosure: {interval['minimum_lower']}
```

Since `e<=2/5`, this gives `D>875/64` in every near-removal box. The
direct phase minimum gives `D>4/5` away from those boxes, while the imported
base cell gives `D>3/50`. Hence the global fixed-prefix margin is
`P_M<-3/8000`.

The 30-term sinc derivative series includes an explicit tail ball smaller
than `1e-50`; the argument bound is below `2` on the whole compact box.

## Theorem and Boundary

{exact['theorem']}

The exact endpoint check is `{exact['endpoint_diagnostic']}`.

{exact['pi_provenance']}

This is not yet a complete aggregate theorem. The near-terminal limits are
proved with `M` fixed before `N -> infinity`; they do not justify a prefix
whose length grows with height. {exact['surviving_target']}

No growing-prefix finite-height bound, complete cross-current estimate, Abel
gap, winding cap, contact exclusion, `Q209`, cofinal descendant theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is claimed.

## Reproduce

```text
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```
"""


def build_payload() -> dict:
    payloads = load_sources()
    audit = source_audit(payloads)
    symbolic = symbolic_certificate()
    interval = compact_box_certificate()
    exact = exact_payload()
    rows = build_rows(exact, symbolic, interval)
    return {
        "kind": STEM,
        "date": "2026-07-31",
        "status": (
            "exact finite contiguous-tail recurrence and Arb-uniform fixed-M "
            "leading current theorem; no growing-tail finite-height estimate, "
            "aggregate cross-current bound, Abel gap, contact exclusion, "
            "Lambda<=0, or RH"
        ),
        "proof_boundary": (
            "This proves the exact real-trace recurrence collapse, its first "
            "two derivatives, and strict clockwise q=1 leading current for "
            "every fixed finite contiguous terminal prefix. It proves no "
            "estimate uniform when M grows with N, no finite-height grouped "
            "tail sign, complete Xi aggregate cross-current estimate, Abel "
            "gap, winding cap, contact exclusion, Q209, cofinal descendant "
            "theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
        "source_audit": audit,
        "symbolic_certificate": symbolic,
        "interval_certificate": interval,
        "exact": exact,
        "counts": {
            "rows": len(rows),
            "exact_endpoint_recurrences": 1,
            "exact_tail_collapses": 1,
            "grouped_jet_identities": 3,
            "arb_compact_boxes": interval["boxes"],
            "fixed_finite_tail_current_theorems": 1,
            "growing_tail_estimates": 0,
            "complete_aggregate_cross_current_bounds": 0,
            "abel_gaps": 0,
            "contact_exclusions": 0,
        },
        "rows": [asdict(row) for row in rows],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    atomic_write(args.out, json.dumps(payload, indent=2) + "\n")
    atomic_write(args.note, render_note(payload))
    print(
        "wrote contiguous terminal-tail recurrence/current gate: "
        f"{payload['counts']['rows']} rows, "
        f"{payload['interval_certificate']['boxes']} Arb boxes, "
        "1 all-fixed-M current theorem, 0 growing-tail estimates"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
