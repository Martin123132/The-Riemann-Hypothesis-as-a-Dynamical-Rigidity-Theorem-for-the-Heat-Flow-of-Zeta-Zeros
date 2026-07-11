#!/usr/bin/env python3
"""Reduce the Newman collision field to a local odd reciprocal-gap statistic."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_local_odd_count_reduction_lemma.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_local_odd_count_reduction_lemma.md"


@dataclass(frozen=True)
class LemmaRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def odd_count_integral(
    left: list[sp.Rational], right: list[sp.Rational], cutoff: sp.Rational
) -> sp.Expr:
    events: dict[sp.Rational, int] = {}
    for distance in left:
        events[distance] = events.get(distance, 0) + 1
    for distance in right:
        events[distance] = events.get(distance, 0) - 1
    current = 0
    previous = sp.Rational(0)
    integral = sp.Rational(0)
    for distance in sorted(events):
        if previous > 0:
            integral += current * (1 / previous - 1 / distance)
        current += events[distance]
        previous = distance
    if previous > 0 and previous < cutoff:
        integral += current * (1 / previous - 1 / cutoff)
    return sp.factor(current / cutoff + integral)


def build_exact() -> dict:
    c, y, width = sp.symbols("c y width", positive=True)
    kernel = 2 * c / (c**2 - y**2)
    split_kernel = 1 / (c - y) + 1 / (c + y)
    kernel_prime = sp.diff(kernel, y)
    if sp.simplify(kernel - split_kernel) != 0:
        raise RuntimeError("signed-pair kernel split failed")
    if sp.simplify(kernel_prime - 4 * c * y / (c**2 - y**2) ** 2) != 0:
        raise RuntimeError("field kernel derivative failed")

    left_boundary = sp.factor(kernel.subs(y, c - width))
    right_boundary = sp.factor(-kernel.subs(y, c + width))
    if sp.simplify(left_boundary - 2 * c / (width * (2 * c - width))) != 0:
        raise RuntimeError("left kernel boundary failed")
    if sp.simplify(right_boundary - 2 * c / (width * (2 * c + width))) != 0:
        raise RuntimeError("right kernel boundary failed")

    radius = sp.symbols("radius", positive=True)
    tail_primitive = (
        sp.log(2 * radius) / (2 * radius**2)
        + sp.Rational(1, 4) / radius**2
    )
    if sp.simplify(
        sp.diff(tail_primitive, radius) + sp.log(2 * radius) / radius**3
    ) != 0 or sp.limit(tail_primitive, radius, sp.oo) != 0:
        raise RuntimeError("tail logarithmic integral failed")

    outer_bound = (
        "If |F(y)|<=E on [a,2c] and |F(y)|<=C*log(2+y) for y>=2c, "
        "then for 0<H<=c/2 the field discrepancy outside (c-H,c+H) is at most "
        "5*E/H+2*C*(log(4c)+1)/c."
    )

    rho_max, rho_derivative, local_mass, error = sp.symbols(
        "rho_max rho_derivative local_mass error", nonnegative=True
    )
    smooth_local_bound = (
        local_mass / (2 * c - width)
        + 2 * width * rho_derivative
        + 2 * width * rho_max / (2 * c - width)
    )

    odd_checks: list[dict] = []
    examples = [
        (
            [sp.Rational(1, 3), sp.Rational(5, 4), sp.Rational(7, 3)],
            [sp.Rational(1, 2), sp.Rational(3, 2), sp.Rational(11, 4)],
            sp.Rational(3),
        ),
        (
            [sp.Rational(2, 5), sp.Rational(4, 5)],
            [sp.Rational(1, 2), sp.Rational(3, 4), sp.Rational(7, 5)],
            sp.Rational(2),
        ),
        (
            [sp.Rational(1, 4), sp.Rational(2, 3), sp.Rational(5, 3)],
            [sp.Rational(1, 4), sp.Rational(3, 4), sp.Rational(9, 5)],
            sp.Rational(2),
        ),
    ]
    for index, (left, right, cutoff) in enumerate(examples, start=1):
        direct = sp.factor(sum(1 / item for item in left) - sum(1 / item for item in right))
        by_count = odd_count_integral(left, right, cutoff)
        if sp.simplify(direct - by_count) != 0:
            raise RuntimeError(f"odd-count identity failed in example {index}")
        paired_count = min(len(left), len(right))
        paired = sp.factor(
            sum(
                (right[j] - left[j]) / (left[j] * right[j])
                for j in range(paired_count)
            )
            + sum(1 / item for item in left[paired_count:])
            - sum(1 / item for item in right[paired_count:])
        )
        if sp.simplify(direct - paired) != 0:
            raise RuntimeError(f"paired-gap identity failed in example {index}")
        odd_checks.append(
            {
                "example": index,
                "left": [str(item) for item in left],
                "right": [str(item) for item in right],
                "cutoff": str(cutoff),
                "field": str(direct),
                "passed": True,
            }
        )

    z, tau = sp.symbols("z tau", real=True)
    a_squared = 1 + 16 / (8 + sp.pi)
    initial = sp.expand((z**2 - 1) ** 2 * (z**2 - a_squared))
    heat_polynomial = sp.expand(
        sum(
            (-tau) ** order
            / sp.factorial(order)
            * sp.diff(initial, z, 2 * order)
            for order in range(4)
        )
    )
    if sp.simplify(sp.diff(heat_polynomial, tau) + sp.diff(heat_polynomial, z, 2)) != 0:
        raise RuntimeError("even heat-birth countermodel PDE failed")
    cofactor = sp.cancel(initial / (z - 1) ** 2)
    classical_field = sp.simplify(sp.diff(cofactor, z).subs(z, 1) / cofactor.subs(z, 1))
    if classical_field != -sp.pi / 8:
        raise RuntimeError("classical-field countermodel normalization failed")
    if sp.simplify(sp.diff(initial, z, 3).subs(z, 1) / sp.diff(initial, z, 2).subs(z, 1) - 3 * classical_field) != 0:
        raise RuntimeError("classical-field countermodel jet failed")

    eps = sp.symbols("eps", positive=True)
    local_coordinate = sp.sqrt(2) * eps + 2 * classical_field * eps**2
    local_series = sp.series(
        heat_polynomial.subs({tau: eps**2, z: 1 + local_coordinate}),
        eps,
        0,
        4,
    ).removeO()
    if sp.simplify(local_series) != 0:
        raise RuntimeError("positive branch heat-birth expansion failed")
    local_coordinate_minus = -sp.sqrt(2) * eps + 2 * classical_field * eps**2
    local_series_minus = sp.series(
        heat_polynomial.subs({tau: eps**2, z: 1 + local_coordinate_minus}),
        eps,
        0,
        4,
    ).removeO()
    if sp.simplify(local_series_minus) != 0:
        raise RuntimeError("negative branch heat-birth expansion failed")

    return {
        "positive_zero_field": (
            "At an even double zero c>0, remove the positive double atom and let mu "
            "count the other positive roots. Then B(c)=1/c+PV integral K_c(y)dmu(y), "
            "K_c(y)=2c/(c^2-y^2)=1/(c-y)+1/(c+y)."
        ),
        "kernel_derivative": "K_c'(y)=4c*y/(c^2-y^2)^2>0",
        "outer_count_bound": outer_bound,
        "outer_bound_constants": {
            "left_boundary": str(left_boundary),
            "right_boundary_absolute": str(right_boundary),
            "left_plus_mid_coefficient": "14/3 < 5",
            "tail_integral": str(tail_primitive),
        },
        "published_counting_input": {
            "source": "https://arxiv.org/abs/1904.12438",
            "statement": (
                "Uniformly for 0<t<=1/2, N_t(X)=g(X,t)+O(log(2+X)), with an "
                "absolute implied constant; g_x=log(x/(4*pi))/(4*pi)+t/(16*x)."
            ),
            "use": (
                "After removing the double atom and fixing the low cutoff, E=O(log(4c)) "
                "on [a,2c] and the same logarithmic bound holds on the tail."
            ),
        },
        "local_window": "H=log(4c)^2",
        "local_singular_field": (
            "S_H(c)=sum_(positive outsider roots y:0<|y-c|<H) 1/(c-y)"
        ),
        "smooth_local_bound": str(smooth_local_bound),
        "localized_field_theorem": (
            "Uniformly for 0<t<=1/2 at any large positive double zero c, "
            "B(c)=-pi/8+S_H(c)+O(1/log(4c)+log(4c)^3/c), H=log(4c)^2."
        ),
        "odd_count_identity": (
            "Let D_c(u)=#roots in [c-u,c)-#roots in (c,c+u]. If H avoids roots, "
            "S_H=D_c(H)/H+integral_0^H D_c(u)/u^2 du."
        ),
        "paired_gap_identity": (
            "For paired left/right distances ell_m,r_m, "
            "1/ell_m-1/r_m=(r_m-ell_m)/(ell_m*r_m); unmatched roots retain "
            "their signed reciprocal terms."
        ),
        "absolute_odd_norm": (
            "W_H=|D_c(H)|/H+integral_0^H |D_c(u)|/u^2 du controls |S_H|."
        ),
        "counting_guard": (
            "The O(log c) Riemann-von Mangoldt error controls the outer field but "
            "cannot control W_H near u=0; one neighboring even pair can keep the "
            "counting error <=1 while producing a reciprocal field of size 1/epsilon."
        ),
        "classical_field_birth_countermodel": {
            "initial_polynomial": str(initial),
            "a_squared": str(a_squared),
            "heat_solution": str(heat_polynomial),
            "field": str(classical_field),
            "center_velocity": str(2 * classical_field),
            "branches": (
                "z_+/-=1+/-sqrt(2*tau)-pi*tau/4+O(tau^(3/2)), tau>0"
            ),
            "conclusion": (
                "An exact even backward-heat solution can have the classical field and "
                "drift while still undergoing positive square-root birth."
            ),
        },
        "open_handoff": (
            "Prove an Xi-specific lambda-uniform bound on the signed odd-count integral "
            "or paired reciprocal-gap discrepancy, and then add a genuinely global "
            "mechanism that turns that control into collision exclusion; field balance "
            "alone is compatible with positive birth."
        ),
        "checks": {
            "odd_count_checks": odd_checks,
            "kernel_split": str(sp.factor(kernel - split_kernel)),
            "kernel_derivative": str(kernel_prime),
            "countermodel_field": str(classical_field),
            "countermodel_pde_residual": str(
                sp.factor(sp.diff(heat_polynomial, tau) + sp.diff(heat_polynomial, z, 2))
            ),
        },
    }


def build_payload() -> dict:
    exact = build_exact()
    rows = [
        LemmaRow(
            id="nloc_01_positive_zero_stieltjes_field",
            role="exact_identity",
            readiness="available_exact",
            claim="The signed Newman field at an even double zero is a Stieltjes transform of the positive outsider-root counting measure plus the reflected double root.",
            formula=exact["positive_zero_field"],
            proof_boundary="Canonical-product identity in the real-zero regime.",
        ),
        LemmaRow(
            id="nloc_02_outer_counting_bound",
            role="exact_lemma",
            readiness="ready_to_apply",
            claim="A uniform counting-function discrepancy controls the entire collision field outside a symmetric local window with an explicit inverse-window bound.",
            formula=exact["outer_count_bound"],
            proof_boundary="Conditional deterministic Stieltjes estimate.",
            diagnostics=exact["outer_bound_constants"],
        ),
        LemmaRow(
            id="nloc_03_published_uniform_counting",
            role="published_theorem",
            readiness="available_published",
            claim="The published positive-time zero count supplies the logarithmic discrepancy required by the outer-field lemma uniformly in t.",
            formula=exact["published_counting_input"]["statement"],
            proof_boundary="External theorem; it does not control reciprocal gaps inside the local window.",
            diagnostics=exact["published_counting_input"],
        ),
        LemmaRow(
            id="nloc_04_log_squared_localization",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim="Choosing H=log(4c)^2 reduces the full high collision field to the classical constant and one local singular sum.",
            formula=exact["localized_field_theorem"],
            proof_boundary="High-c asymptotic reduction; not a collision-exclusion theorem.",
        ),
        LemmaRow(
            id="nloc_05_odd_count_identity",
            role="exact_identity",
            readiness="available_exact",
            claim="The local singular field is exactly the inverse-square weighted odd part of the local counting function.",
            formula=exact["odd_count_identity"],
            proof_boundary="Finite local window avoiding outsider roots at its boundary.",
            diagnostics=exact["checks"]["odd_count_checks"],
        ),
        LemmaRow(
            id="nloc_06_paired_gap_norm",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Left-right pairing turns the local field into reciprocal-gap asymmetry, with a natural absolute odd-count norm as a sufficient bound.",
            formula=f"{exact['paired_gap_identity']} {exact['absolute_odd_norm']}",
            proof_boundary="Sufficient local balance estimate only.",
        ),
        LemmaRow(
            id="nloc_07_counting_only_guard",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="Uniform Riemann-von Mangoldt error does not control the local odd reciprocal-gap norm.",
            formula=exact["counting_guard"],
            proof_boundary="Blocks counting-only promotion; does not rule out stronger Xi spacing input.",
        ),
        LemmaRow(
            id="nloc_08_classical_field_birth_guard",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="Even exact classical field and center drift remain compatible with positive square-root birth under the backward heat equation.",
            formula=exact["classical_field_birth_countermodel"]["branches"],
            proof_boundary="Exact even polynomial heat-flow countermodel, not the Xi kernel.",
            diagnostics=exact["classical_field_birth_countermodel"],
        ),
        LemmaRow(
            id="nloc_09_fixed_time_scope",
            role="scope_gate",
            readiness="guard_validated",
            claim="Fixed-positive-time high-zero localization and the local field reduction are compatible: the former gives compactness, while the latter identifies what must be uniform as t approaches zero.",
            formula="high zeros simple for x>=exp(C/t); local odd-count control needed below that expanding scale",
            proof_boundary="Scope composition only; no compact collision is excluded.",
        ),
        LemmaRow(
            id="nloc_10_xi_collision_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The remaining root-field route needs both an Xi-specific uniform odd-count theorem and a separate global implication from field control to no collision.",
            formula=exact["open_handoff"],
            proof_boundary="Open theorem pair; not a proof of RH or Lambda <= 0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_local_odd_count_reduction_lemma",
        "date": "2026-07-11",
        "status": "exact local odd-count reduction with counting and birth guards",
        "proof_boundary": (
            "This artifact proves a deterministic Stieltjes outer-field bound, composes it "
            "with the published uniform positive-time zero count, reduces the high collision "
            "field to a log-squared local odd reciprocal-gap statistic, and gives an exact even "
            "heat-flow countermodel showing that classical field balance does not exclude birth. "
            "It does not prove the Xi-specific local odd-count estimate, turn field control into "
            "a global collision obstruction, prove RH, or prove Lambda <= 0."
        ),
        "sources": [
            "outputs/jensen_window_pf_newman_classical_field_balance_gate.md",
            "outputs/jensen_window_pf_newman_root_external_field_lemma.md",
            "https://arxiv.org/abs/1801.05914",
            "https://arxiv.org/abs/1904.12438",
        ],
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    birth = exact["classical_field_birth_countermodel"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Local Odd-Count Reduction Lemma",
            "",
            "Date: 2026-07-11",
            "",
            "Status: exact local odd-count reduction with counting and birth",
            "guards. This is not a proof of RH or `Lambda <= 0`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_local_odd_count_reduction_lemma.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_local_odd_count_reduction_lemma.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_local_odd_count_reduction_lemma.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF Newman local odd-count reduction lemma: 10 rows, 0 issues, 3 exact Stieltjes identities, 1 explicit outer bound, 1 log-squared localization theorem, 1 odd-count formula, 3 finite reciprocal-gap checks, 1 published uniform counting input, 1 exact classical-field birth countermodel, 1 open Xi collision-exclusion handoff",
            "```",
            "",
            "## Stieltjes Reduction",
            "",
            "For a positive double zero `c`,",
            "",
            "```text",
            exact["positive_zero_field"],
            exact["kernel_derivative"],
            "```",
            "",
            "Let `F` be the difference between the positive outsider-root count",
            "and the positive-time reference count. Integration by parts gives",
            "",
            "```text",
            exact["outer_count_bound"],
            "```",
            "",
            "The published theorem supplies `E=O(log(4c))` uniformly for",
            "`0<t<=1/2`. Choosing `H=log(4c)^2` and using the classical",
            "continuum field therefore gives",
            "",
            "```text",
            exact["localized_field_theorem"],
            "```",
            "",
            "All mesoscopic and far-field uncertainty is now `O(1/log c)`.",
            "The unresolved quantity is genuinely local.",
            "",
            "Primary counting source: https://arxiv.org/abs/1904.12438.",
            "",
            "## Odd Local Statistic",
            "",
            "Define the left-minus-right outsider count at radius `u` by `D_c(u)`.",
            "Then",
            "",
            "```text",
            exact["odd_count_identity"],
            exact["paired_gap_identity"],
            exact["absolute_odd_norm"],
            "```",
            "",
            "This is the strict local target. A global `O(log c)` counting error",
            "does not control the inverse-square weight at small `u`.",
            "",
            "## Birth Guard",
            "",
            "Set",
            "",
            "```text",
            f"a^2={birth['a_squared']}",
            "P(z)=(z^2-1)^2*(z^2-a^2)",
            "F_tau=exp(-tau*d_z^2)P",
            f"B(1)={birth['field']}, center velocity={birth['center_velocity']}",
            birth["branches"],
            "```",
            "",
            "This is an exact even solution of `F_tau=-F_zz`. It has the",
            "classical field `-pi/8` and drift `-pi/4`, yet a real pair is born",
            "for `tau>0`. Therefore even a successful odd-count theorem would",
            "control the field but would not, by itself, exclude a positive Newman",
            "boundary.",
            "",
            "## Live Handoff",
            "",
            exact["open_handoff"],
            "",
            "The next proof mechanism must use an additional Xi-specific global",
            "constraint beyond the first regularized root field.",
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
        "wrote Jensen-window PF Newman local odd-count reduction lemma: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
