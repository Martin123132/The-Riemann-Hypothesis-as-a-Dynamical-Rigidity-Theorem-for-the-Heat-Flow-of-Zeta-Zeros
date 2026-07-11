#!/usr/bin/env python3
"""Build an exact quartic boundary-flow obstruction to cubic-cone promotion."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_quartic_boundary_flow_obstruction.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_quartic_boundary_flow_obstruction.md"


@dataclass(frozen=True)
class ObstructionRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def rational_text(value: sp.Expr) -> str:
    return str(sp.factor(value))


def build_exact_data() -> dict:
    w = sp.symbols("w")
    x, y, z, u, k = sp.symbols("x y z u k")
    quartic = 1 + 4 * w + 6 * x * w**2 + 4 * x**2 * y * w**3 + x**3 * y**2 * z * w**4
    discriminant = sp.factor(sp.discriminant(quartic, w))
    Q = sp.factor(discriminant / (256 * x**6 * y**2))
    if sp.simplify(discriminant - 256 * x**6 * y**2 * Q) != 0:
        raise RuntimeError("quartic discriminant normalization failed")

    derivative_discriminant = sp.factor(sp.discriminant(sp.diff(quartic, w), w))
    cubic_frontier_yz = y**2 * z**2 - 6 * y * z + 4 * y + 4 * z - 3
    if sp.simplify(
        derivative_discriminant + 6912 * x**6 * y**2 * cubic_frontier_yz
    ) != 0:
        raise RuntimeError("quartic derivative cubic factorization failed")

    x_prime = 2 * ((2 * k - 1) * (1 - x) - (2 * k + 3) * x * (1 - y))
    y_prime = 2 * y * ((2 * k + 1) * (1 - y) - (2 * k + 5) * y * (1 - z))
    z_prime = 2 * y * z * ((2 * k + 3) * (1 - z) - (2 * k + 7) * z * (1 - u))
    Q_prime = sp.factor(
        sp.diff(Q, x) * x_prime + sp.diff(Q, y) * y_prime + sp.diff(Q, z) * z_prime
    )
    if sp.degree(Q_prime, u) != 1:
        raise RuntimeError("quartic frontier derivative is not linear in the next contraction")

    a = sp.Rational(13, 20)
    b = sp.Rational(21, 50)
    c = sp.Rational(57, 25)
    if 2 * a + b + c != 4:
        raise RuntimeError("countermodel roots are not normalized")
    factored_quartic = sp.expand((1 + a * w) ** 2 * (1 + b * w) * (1 + c * w))
    e2 = a**2 + 2 * a * b + 2 * a * c + b * c
    e3 = a**2 * b + a**2 * c + 2 * a * b * c
    e4 = a**2 * b * c
    xv = sp.factor(e2 / 6)
    yv = sp.factor((e3 / 4) / xv**2)
    zv = sp.factor(e4 / (xv**3 * yv**2))
    uv = zv
    normalized_quartic = sp.expand(quartic.subs({x: xv, y: yv, z: zv}))
    if sp.expand(normalized_quartic - factored_quartic) != 0:
        raise RuntimeError("countermodel quartic factorization failed")

    cubic_frontier = lambda left, right: sp.factor(
        left**2 * right**2 - 6 * left * right + 4 * left + 4 * right - 3
    )
    cubic_margins = {
        "F_x_y": cubic_frontier(xv, yv),
        "F_y_z": cubic_frontier(yv, zv),
        "F_z_u": cubic_frontier(zv, uv),
    }
    if not all(value < 0 for value in cubic_margins.values()):
        raise RuntimeError("countermodel does not lie strictly inside every local cubic cone")

    order_margins = {
        "y_minus_x": sp.factor(yv - xv),
        "z_minus_y": sp.factor(zv - yv),
        "u_minus_z": sp.factor(uv - zv),
        "one_minus_u": sp.factor(1 - uv),
    }
    if not all(value >= 0 for value in order_margins.values()):
        raise RuntimeError("countermodel contractions are not increasing")
    lower_margins = {
        "x_minus_1_over_3": sp.factor(xv - sp.Rational(1, 3)),
        "y_minus_3_over_5": sp.factor(yv - sp.Rational(3, 5)),
        "z_minus_5_over_7": sp.factor(zv - sp.Rational(5, 7)),
        "u_minus_7_over_9": sp.factor(uv - sp.Rational(7, 9)),
    }
    if not all(value > 0 for value in lower_margins.values()):
        raise RuntimeError("countermodel misses a pointwise ratio wall")

    Q_value = sp.factor(Q.subs({x: xv, y: yv, z: zv}))
    Q_prime_value = sp.factor(Q_prime.subs({x: xv, y: yv, z: zv, u: uv, k: 1}))
    if Q_value != 0 or not Q_prime_value < 0:
        raise RuntimeError("countermodel is not an outward quartic boundary point")

    return {
        "quartic": str(quartic),
        "quartic_discriminant_factorization": "Disc(J_4)=256*x^6*y^2*Q(x,y,z)",
        "Q": str(Q),
        "derivative_discriminant_factorization": (
            "Disc(partial_w J_4)=-6912*x^6*y^2*F(y,z)"
        ),
        "cubic_frontier": "F(s,t)=s^2*t^2-6*s*t+4*s+4*t-3",
        "heat_vector_over_r_k": {
            "x_prime": str(x_prime),
            "y_prime": str(y_prime),
            "z_prime": str(z_prime),
        },
        "Q_prime_linear_in_u": True,
        "countermodel": {
            "k": 1,
            "a": rational_text(a),
            "b": rational_text(b),
            "c": rational_text(c),
            "factorization": "(1+13*w/20)^2*(1+21*w/50)*(1+57*w/25)",
            "x": rational_text(xv),
            "y": rational_text(yv),
            "z": rational_text(zv),
            "u": rational_text(uv),
            "order_margins": {name: rational_text(value) for name, value in order_margins.items()},
            "pointwise_lower_margins": {
                name: rational_text(value) for name, value in lower_margins.items()
            },
            "cubic_frontiers": {
                name: rational_text(value) for name, value in cubic_margins.items()
            },
            "Q": rational_text(Q_value),
            "Q_prime_over_r_1": rational_text(Q_prime_value),
        },
    }


def build_payload() -> dict:
    exact = build_exact_data()
    counter = exact["countermodel"]
    rows = [
        ObstructionRow(
            id="qbfo_01_normalized_quartic",
            role="exact_identity",
            readiness="available_exact",
            claim="Every shifted quartic Jensen window is represented by three consecutive contraction coordinates after positive rescaling.",
            formula="J_4=1+4*w+6*x*w^2+4*x^2*y*w^3+x^3*y^2*z*w^4",
            proof_boundary="Degree-4 normalization only.",
        ),
        ObstructionRow(
            id="qbfo_02_quartic_discriminant",
            role="exact_identity",
            readiness="available_exact",
            claim="The quartic discriminant has a positive monomial prefactor and a 16-term frontier polynomial Q.",
            formula=exact["quartic_discriminant_factorization"],
            proof_boundary="Quartic discriminant identity only.",
            diagnostics={"Q": exact["Q"]},
        ),
        ObstructionRow(
            id="qbfo_03_derivative_cubic_frontier",
            role="exact_identity",
            readiness="available_exact",
            claim="The derivative cubic is strictly hyperbolic exactly on the shifted cubic frontier sign F(y,z)<0.",
            formula=exact["derivative_discriminant_factorization"],
            proof_boundary="Derivative-cubic identity only.",
        ),
        ObstructionRow(
            id="qbfo_04_heat_frontier_derivative",
            role="exact_identity",
            readiness="available_exact",
            claim="The heat derivative of Q is explicit and linear in the next contraction u.",
            formula="Q'=Q_x*x'+Q_y*y'+Q_z*z', degree_u(Q')=1",
            proof_boundary="Local quartic boundary vector field only.",
            diagnostics=exact["heat_vector_over_r_k"],
        ),
        ObstructionRow(
            id="qbfo_05_hyperbolic_boundary_point",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="The rational example is a normalized hyperbolic quartic with one double negative root.",
            formula=counter["factorization"],
            proof_boundary="Abstract quartic boundary point, not a zeta coefficient window.",
        ),
        ObstructionRow(
            id="qbfo_06_ratio_cone_membership",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="Its four local contractions satisfy the k=1 pointwise walls and are nondecreasing.",
            formula="1/3<x<y<z=u<1 and 3/5<y, 5/7<z, 7/9<u",
            proof_boundary="Local ratio-cone membership only.",
            diagnostics={
                "coordinates": {key: counter[key] for key in ("x", "y", "z", "u")},
                "order_margins": counter["order_margins"],
                "pointwise_lower_margins": counter["pointwise_lower_margins"],
            },
        ),
        ObstructionRow(
            id="qbfo_07_strict_cubic_membership",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="Every neighboring cubic inequality is strict at the quartic boundary point.",
            formula="F(x,y)<0, F(y,z)<0, F(z,u)<0",
            proof_boundary="Three local cubic tests only.",
            diagnostics=counter["cubic_frontiers"],
        ),
        ObstructionRow(
            id="qbfo_08_outward_quartic_flow",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="The heat vector points strictly from the hyperbolic quartic boundary into negative discriminant.",
            formula=f"Q=0, Q'/r_1={counter['Q_prime_over_r_1']}<0",
            proof_boundary="Local abstract boundary-exit countermodel only.",
        ),
        ObstructionRow(
            id="qbfo_09_promotion_guard",
            role="forbidden_promotion",
            readiness="guard_validated",
            claim="The propagated ratio cone and all shifted cubic inequalities do not by themselves imply quartic heat invariance.",
            formula="ratio cone + cubic cone + local heat ODE != quartic invariant",
            proof_boundary="Blocks one promotion route; it does not show the actual zeta quartics fail.",
        ),
        ObstructionRow(
            id="qbfo_10_coupled_invariant_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="Degree four requires an additional coupled frontier condition beyond consecutive cubic increments.",
            formula="construct a Q-compatible condition involving x,y,z,u and propagate it",
            proof_boundary="Open quartic and all-degree bridge; not PF-infinity, RH, or Lambda <= 0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_quartic_boundary_flow_obstruction",
        "date": "2026-07-10",
        "status": "exact quartic boundary-flow countermodel and cubic-promotion guard",
        "proof_boundary": (
            "This artifact proves that the ratio cone plus all local cubic inequalities is not a sufficient quartic heat invariant. The countermodel is abstract and local; it does not show failure of the actual zeta trajectory, prove or disprove degree-4 zeta hyperbolicity, PF-infinity, RH, or Lambda <= 0."
        ),
        "sources": [
            "outputs/jensen_window_pf_cubic_forward_uniform_tail_certificate.md",
            "outputs/jensen_window_pf_heat_flow_jensen_hierarchy_lemma.md",
            "outputs/jensen_window_pf_bridge_target.md",
        ],
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    counter = exact["countermodel"]
    return "\n".join(
        [
            "# Jensen-Window PF Quartic Boundary-Flow Obstruction",
            "",
            "Date: 2026-07-10",
            "",
            "Status: exact quartic boundary-flow countermodel. This blocks promotion",
            "of the cubic cone but is not a counterexample to the zeta trajectory and",
            "is not a proof or disproof of RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_quartic_boundary_flow_obstruction.json",
            "python work/rh_compute/scripts/jensen_window_pf_quartic_boundary_flow_obstruction.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_quartic_boundary_flow_obstruction.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF quartic boundary-flow obstruction: 10 rows, 0 issues, 4 exact quartic identities, 1 hyperbolic boundary point, 4 positive ratio margins, 3 strict cubic margins, 1 negative quartic derivative, 1 blocked promotion, 1 open coupled-invariant handoff",
            "```",
            "",
            "## Quartic Frontier",
            "",
            "After positive rescaling, a shifted quartic window is",
            "",
            "```text",
            "J_4=1+4*w+6*x*w^2+4*x^2*y*w^3+x^3*y^2*z*w^4.",
            exact["quartic_discriminant_factorization"],
            "```",
            "",
            "The derivative cubic satisfies",
            "",
            "```text",
            exact["derivative_discriminant_factorization"],
            exact["cubic_frontier"],
            "```",
            "",
            "The exact heat vector for `(x,y,z)` makes `Q'` linear in the next",
            "contraction `u`.",
            "",
            "## Exact Boundary Countermodel",
            "",
            "Take `k=1` and",
            "",
            "```text",
            counter["factorization"],
            f"x={counter['x']}",
            f"y={counter['y']}",
            f"z=u={counter['z']}",
            "```",
            "",
            "The displayed factorization gives four negative roots with one double",
            "root, so `Q=0` is a genuine hyperbolic quartic boundary point. The",
            "contractions are nondecreasing and satisfy every k=1 pointwise wall.",
            "All neighboring cubic tests are strict:",
            "",
            "```text",
            f"F(x,y)={counter['cubic_frontiers']['F_x_y']}<0",
            f"F(y,z)={counter['cubic_frontiers']['F_y_z']}<0",
            f"F(z,u)={counter['cubic_frontiers']['F_z_u']}<0",
            "```",
            "",
            "Nevertheless, direct substitution into the heat vector gives",
            "",
            "```text",
            f"Q=0, Q'/r_1={counter['Q_prime_over_r_1']}<0.",
            "```",
            "",
            "Because `Disc(J_4)=256*x^6*y^2*Q`, the discriminant immediately becomes",
            "negative. Thus the ratio cone plus all cubic inequalities is not a",
            "quartic invariant.",
            "",
            "## Consequence",
            "",
            "Degree four needs a new coupled condition involving the quartic frontier",
            "and the next contraction. This abstract guard says nothing adverse about",
            "the actual zeta coefficient trajectory; it only rejects an insufficient",
            "proof shortcut. PF-infinity, RH, and `Lambda <= 0` remain open.",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Jensen-window PF quartic boundary-flow obstruction: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
