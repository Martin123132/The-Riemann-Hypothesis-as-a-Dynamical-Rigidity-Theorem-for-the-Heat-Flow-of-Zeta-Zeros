#!/usr/bin/env python3
"""Build the arbitrary-regular-field first-loss heat-flow nonclosure gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_regular_field_first_loss_nonclosure_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "first_xi_jet": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_arbitrary_multiplicity_first_xi_jet_layer_gate.json",
    "positive_boundary": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_positive_boundary_attainment_lemma.json",
    "first_jet_winding": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_first_jet_winding_gate.json",
}

AUDIT_FIELDS = (
    sp.Rational(-5, 2),
    sp.Rational(-7, 4),
    sp.Rational(-1, 3),
    sp.Rational(0),
    sp.Rational(2, 5),
    sp.Rational(1, 2),
    sp.Rational(4, 3),
    sp.Rational(11, 6),
    sp.Rational(3),
    sp.Rational(7, 2),
)
MIN_MULTIPLICITY = 2
MAX_AUDIT_MULTIPLICITY = 16


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    certificate: str
    proof_boundary: str


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def expression_hash(expression: sp.Expr) -> str:
    return hashlib.sha256(sp.srepr(sp.expand(expression)).encode("utf-8")).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def load_sources() -> dict[str, dict]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        require(path.exists(), f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))

    first_xi = payloads["first_xi_jet"].get("symbolic_certificate", {})
    require(
        "ell=rho U'(rho)/U(rho)" in first_xi.get("full_expansion", ""),
        "first-Xi regular-field contract drift",
    )
    require(
        "positive Newman boundary" in payloads["positive_boundary"]
        .get("exact", {})
        .get("attainment_theorem", "")
        or "If Lambda>0" in payloads["positive_boundary"]
        .get("exact", {})
        .get("attainment_theorem", ""),
        "positive-boundary contract drift",
    )
    require(
        "floor(m/2)" in payloads["first_jet_winding"]
        .get("exact", {})
        .get("multiplicity_index", ""),
        "first-jet multiplicity contract drift",
    )
    return payloads


def source_audit() -> dict[str, dict[str, str]]:
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def heat_monomial(m: int, z: sp.Symbol, direction: int) -> sp.Expr:
    """Return exp(direction*D_z^2) z^m for direction in {-1,+1}."""
    require(direction in (-1, 1), "heat direction must be -1 or +1")
    return sp.expand(
        sum(
            direction**order
            * sp.diff(z**m, z, 2 * order)
            / sp.factorial(order)
            for order in range(m // 2 + 1)
        )
    )


def evolve_polynomial(polynomial: sp.Expr, x: sp.Symbol, tau: sp.Rational) -> sp.Expr:
    degree = sp.Poly(polynomial, x).degree()
    return sp.expand(
        sum(
            (-tau) ** order
            * sp.diff(polynomial, x, 2 * order)
            / sp.factorial(order)
            for order in range(degree // 2 + 1)
        )
    )


def model_parameters(field: sp.Expr) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    positive_contribution = sp.expand(field**2 + 2)
    negative_magnitude = sp.expand(field**2 - field + 2)
    inner_scale = sp.cancel(1 - 1 / positive_contribution)
    outer_scale = sp.cancel(1 + 1 / negative_magnitude)
    return positive_contribution, negative_magnitude, inner_scale, outer_scale


def build_symbolic_certificate() -> dict[str, str]:
    ell, s = sp.symbols("ell s", real=True)
    c = sp.symbols("c", positive=True)
    positive, negative, alpha, beta = model_parameters(ell)
    rho = -c**2
    rho_inner = -alpha * c**2
    rho_outer = -beta * c**2
    inner_field = sp.factor(rho / (rho - rho_inner))
    outer_field = sp.factor(rho / (rho - rho_outer))
    unit = (s - rho_inner) * (s - rho_outer)
    unit_field = sp.factor(rho * sp.diff(unit, s).subs(s, rho) / unit.subs(s, rho))

    require(sp.simplify(inner_field - positive) == 0, "inner-field identity drift")
    require(sp.simplify(outer_field + negative) == 0, "outer-field identity drift")
    require(sp.simplify(inner_field + outer_field - ell) == 0, "arbitrary-field identity drift")
    require(sp.simplify(unit_field - ell) == 0, "unit logarithmic derivative drift")

    return {
        "parameters": (
            "A=ell^2+2, B=ell^2-ell+2=(ell-1/2)^2+7/4, "
            "alpha=1-1/A, beta=1+1/B"
        ),
        "ordering": "A>1, B>0, and therefore 0<alpha<1<beta for every real ell.",
        "product": (
            "F_(m,ell,c)(s)=(s+c^2)^m(s+alpha*c^2)(s+beta*c^2), "
            "m>=2, c>0"
        ),
        "regular_field": (
            "At rho=-c^2, rho*U'(rho)/U(rho)=A-B=ell; the inner and outer "
            "simple zeros contribute A and -B exactly."
        ),
        "coefficient_geometry": (
            "All factors s+a*c^2 have a>0, so F has only negative zeros and "
            "strictly positive coefficients."
        ),
        "even_lift": (
            "H_(m,ell,c)(x)=F_(m,ell,c)(-x^2) is even and real-rooted, with "
            "multiplicity m at each of +/-c."
        ),
        "signed_regular_field": (
            "If H(x)=(x-c)^m V(x), then V'(c)/V(c)=m/(2c)-2c*U'(rho)/U(rho)="
            "(m+4ell)/(2c)."
        ),
        "forward_preserver": (
            "For tau>=0, exp(-tau*D_x^2)=lim_(N->infinity)"
            "[(1-sqrt(tau/N)D_x)(1+sqrt(tau/N)D_x)]^N; each real first-order "
            "factor preserves real-rootedness."
        ),
        "forward_cluster": (
            "tau^(-m/2)H_tau(c+sqrt(tau)z)/V(c) -> "
            "P_m^-(z)=exp(-D_z^2)z^m=2^(m/2)He_m(z/sqrt(2))."
        ),
        "common_drift": (
            "P_(m+1)^-=zP_m^--2(P_m^-)' implies every forward cluster root "
            "has x_j(tau)=c+sqrt(tau)h_j+2[V'(c)/V(c)]tau+O(tau^(3/2))."
        ),
        "backward_cluster": (
            "For sigma=-tau>0, sigma^(-m/2)H_(-sigma)(c+sqrt(sigma)z)/V(c) "
            "-> P_m^+(z)=exp(D_z^2)z^m. P_m^+ has no real zero for even m "
            "and only z=0 for odd m."
        ),
        "boundary": (
            "H_tau is real-rooted for every tau>=0 and is not real-rooted for "
            "all sufficiently small tau<0; heat monotonicity makes its exact "
            "Newman threshold tau=0."
        ),
        "nonclosure": (
            "The first-loss heat condition, evenness, finite genus-zero negative-root "
            "geometry, and positive s-coefficients impose no bound or sign on ell."
        ),
        "xi_handoff": (
            "Any bound on the hypothetical Xi field must use structure absent from "
            "these models, such as the theta/Fourier source, global zero placement, "
            "or zero first-jet boundary winding."
        ),
        "uniform_handoff": (
            "The cofinal Jensen route separately still needs a remainder estimate "
            "uniform in the selected degree sequence."
        ),
    }


def build_model_audit() -> list[dict]:
    s, x = sp.symbols("s x")
    rows: list[dict] = []
    c_values = (sp.Rational(1), sp.Rational(3, 2), sp.Rational(2), sp.Rational(5, 3), sp.Rational(4, 3))
    for index, field in enumerate(AUDIT_FIELDS):
        m = 2 + index % 5
        c = c_values[index % len(c_values)]
        positive, negative, alpha, beta = model_parameters(field)
        rho = -c**2
        unit = (s + alpha * c**2) * (s + beta * c**2)
        source = sp.expand((s + c**2) ** m * unit)
        lifted = sp.expand(source.subs(s, -x**2))
        computed_field = sp.factor(rho * sp.diff(unit, s).subs(s, rho) / unit.subs(s, rho))
        regular = sp.cancel(lifted / (x - c) ** m)
        signed_field = sp.factor(sp.diff(regular, x).subs(x, c) / regular.subs(x, c))
        expected_signed = sp.factor((m + 4 * field) / (2 * c))
        coefficients = sp.Poly(source, s).all_coeffs()

        require(positive > 1 and negative > 0, f"parameter positivity drift field={field}")
        require(0 < alpha < 1 < beta, f"root ordering drift field={field}")
        require(computed_field == field, f"regular field drift field={field}")
        require(signed_field == expected_signed, f"signed field drift field={field}")
        require(all(coefficient > 0 for coefficient in coefficients), f"coefficient sign drift field={field}")

        tau = sp.Rational(1, 40 + 3 * index)
        evolved = evolve_polynomial(lifted, x, tau)
        degree = int(sp.Poly(evolved, x).degree())
        real_roots = int(sp.Poly(evolved, x).count_roots(-sp.oo, sp.oo))
        require(real_roots == degree, f"positive-time root audit drift field={field}")

        rows.append(
            {
                "field": str(field),
                "multiplicity": m,
                "c": str(c),
                "A": str(positive),
                "B": str(negative),
                "alpha": str(alpha),
                "beta": str(beta),
                "computed_field": str(computed_field),
                "signed_field": str(signed_field),
                "coefficient_count": len(coefficients),
                "source_sha256": expression_hash(source),
                "positive_time": str(tau),
                "positive_time_degree": degree,
                "positive_time_real_roots": real_roots,
                "positive_time_sha256": expression_hash(evolved),
            }
        )
    return rows


def build_hermite_audit() -> list[dict]:
    z = sp.symbols("z")
    rows: list[dict] = []
    for m in range(MIN_MULTIPLICITY, MAX_AUDIT_MULTIPLICITY + 1):
        forward = heat_monomial(m, z, -1)
        backward = heat_monomial(m, z, 1)
        next_forward = heat_monomial(m + 1, z, -1)
        recurrence = sp.expand(next_forward - z * forward + 2 * sp.diff(forward, z))
        forward_real = int(sp.Poly(forward, z).count_roots(-sp.oo, sp.oo))
        backward_real = int(sp.Poly(backward, z).count_roots(-sp.oo, sp.oo))
        require(recurrence == 0, f"Hermite recurrence drift m={m}")
        require(forward_real == m, f"forward Hermite root count drift m={m}")
        require(backward_real == m % 2, f"backward nonreal count drift m={m}")
        rows.append(
            {
                "multiplicity": m,
                "forward_sha256": expression_hash(forward),
                "backward_sha256": expression_hash(backward),
                "forward_real_roots": forward_real,
                "backward_real_roots": backward_real,
                "common_drift_recurrence": "P_(m+1)^-=zP_m^--2(P_m^-)'",
            }
        )
    return rows


def build_rows(certificate: dict[str, str]) -> list[GateRow]:
    return [
        GateRow("rfn_01_sources", "source chain", "proved", "The first-Xi-jet, positive-boundary, and first-jet contracts are hash-pinned.", "Three parent artifacts are audited before construction.", "No parent theorem is strengthened."),
        GateRow("rfn_02_parameters", "arbitrary-field parameters", "proved", "The two auxiliary root scales exist for every real ell.", certificate["parameters"] + "; " + certificate["ordering"], "This is exact real algebra."),
        GateRow("rfn_03_product", "negative-root product", "proved", "The model is a finite genus-zero Laguerre-Polya product.", certificate["product"], "It is a model source, not Xi itself."),
        GateRow("rfn_04_coefficients", "coefficient geometry", "proved", "Every s-coordinate coefficient is strictly positive.", certificate["coefficient_geometry"], "Coefficient positivity alone is only a shared structural feature."),
        GateRow("rfn_05_field", "regular-field realization", "proved", "Every prescribed real ell is attained exactly.", certificate["regular_field"], "The field is at the selected multiplicity-m root."),
        GateRow("rfn_06_even", "even lift", "proved", "The signed-variable model is even and entirely real-rooted at contact.", certificate["even_lift"], "The lift has two symmetric contacts at +/-c."),
        GateRow("rfn_07_conversion", "signed field conversion", "proved", "The s-field gives an arbitrary signed-variable regular field.", certificate["signed_regular_field"], "This is a local logarithmic derivative identity."),
        GateRow("rfn_08_preserver", "forward heat flow", "proved", "Every nonnegative heat time remains real-rooted.", certificate["forward_preserver"], "The proof is for the finite polynomial model."),
        GateRow("rfn_09_forward", "forward cluster", "proved", "The contact opens into the universal real Hermite cluster.", certificate["forward_cluster"], "Local parabolic scaling only."),
        GateRow("rfn_10_drift", "regular-field drift", "proved", "The field translates every leading cluster root without changing hyperbolicity.", certificate["common_drift"], "The O(tau) formula is local at one contact."),
        GateRow("rfn_11_backward", "backward cluster", "proved", "A nonreal pair appears immediately on the negative-time side.", certificate["backward_cluster"], "Root continuity is applied at fixed finite m."),
        GateRow("rfn_12_boundary", "exact first loss", "proved", "Each model has exact Newman-style threshold zero.", certificate["boundary"], "This is a polynomial heat-flow threshold, not the Xi Newman constant."),
        GateRow("rfn_13_multiplicity", "all-multiplicity audit", "guard_validated", "The construction and local split work for every m>=2.", "Exact Hermite audits cover m=2 through 16; the closed formulas prove arbitrary fixed m.", "Finite audit range checks implementation only."),
        GateRow("rfn_14_nonclosure", "local no-go theorem", "proved", "First-loss dynamics supplies no universal sign or bound for ell.", certificate["nonclosure"], "Xi-specific global restrictions are not refuted."),
        GateRow("rfn_15_xi", "Xi-specific field handoff", "open", "A useful field bound must use actual Xi global structure.", certificate["xi_handoff"], "No such Xi bound is proved here."),
        GateRow("rfn_16_uniform", "degree-uniform handoff", "open", "The Jensen collision route still needs uniform analytic error control.", certificate["uniform_handoff"], "No degree-uniform remainder, collision exclusion, Lambda<=0, or RH is proved."),
    ]


def render_note(artifact: dict) -> str:
    certificate = artifact["symbolic_certificate"]
    counts = artifact["counts"]
    return f"""# Newman Regular-Field First-Loss Nonclosure Gate

Date: 2026-08-03

Status: exact arbitrary-field first-loss countermodel theorem; actual Xi field and degree-uniform Jensen remainder remain open; not a proof of RH.

## Construction

Fix any real `ell`, any `c>0`, and any integer `m>=2`. Set

```text
A=ell^2+2,
B=ell^2-ell+2=(ell-1/2)^2+7/4,
alpha=1-1/A,
beta=1+1/B.
```

Then `A>1`, `B>0`, and `0<alpha<1<beta`. Define

```text
F_(m,ell,c)(s)=(s+c^2)^m(s+alpha*c^2)(s+beta*c^2).
```

Every zero is negative and every coefficient is strictly positive. At the multiplicity-`m` zero `rho=-c^2`, with `F=(s-rho)^m U`, the two simple roots contribute

```text
rho/(rho+alpha*c^2)=A,
rho/(rho+beta*c^2)=-B,
rho U'(rho)/U(rho)=A-B=ell.
```

Thus this family realizes every finite real regular field, not only both signs.

## Exact Heat Boundary

Lift the product to the even real-rooted polynomial

```text
H_(m,ell,c)(x)=F_(m,ell,c)(-x^2).
```

It has multiplicity `m` at both `c` and `-c`. For `tau>=0`,

```text
exp(-tau D_x^2)
 =lim_(N->infinity)[(1-sqrt(tau/N)D_x)(1+sqrt(tau/N)D_x)]^N.
```

Every first-order factor preserves real-rootedness, so `H_tau` is real-rooted for every nonnegative time. Near `c`, the two parabolic limits are

```text
tau^(-m/2) H_tau(c+sqrt(tau)z)/V(c)
 -> exp(-D_z^2)z^m=2^(m/2)He_m(z/sqrt(2)),

sigma^(-m/2) H_(-sigma)(c+sqrt(sigma)z)/V(c)
 -> exp(D_z^2)z^m.
```

The first polynomial has `m` simple real roots. The second has no real root for even `m` and only the root zero for odd `m`, so a nonreal pair appears for every sufficiently small negative time. Heat monotonicity therefore makes zero the exact Newman-style threshold of every model.

## What The Field Does

If `H(x)=(x-c)^m V(x)`, then

```text
B_c=V'(c)/V(c)=(m+4ell)/(2c).
```

Writing `P_m^-=exp(-D_z^2)z^m`, the exact Appell-Hermite recurrence

```text
P_(m+1)^-=z P_m^--2(P_m^-)'
```

shows that every member of the forward cluster has

```text
x_j(tau)=c+sqrt(tau)h_j+2B_c tau+O(tau^(3/2)).
```

So `ell` changes the common drift but cannot change the real-versus-nonreal side of the contact. For the multiplicity-three first-Xi-jet case, `B_c=(3+4ell)/(2c)`.

## Route Decision

This is a genuine nonclosure theorem for the local route. First-loss heat dynamics, evenness, finite genus-zero negative-root geometry, and positive `s`-coefficients do not impose any sign, finite interval, or one-sided bound on `ell`. The earlier two sign examples were not exceptional; every real value occurs inside an exact Newman-boundary model.

The result does not say the actual Xi field is arbitrary. Any closing estimate must use structure absent here: the theta/Fourier source, global Xi zero placement, or zero first-jet boundary winding. The cofinal Jensen route also still requires a remainder uniform in degree.

## Audit

The builder certifies {counts['field_models']} exact rational field models and {counts['multiplicities']} Hermite multiplicities. It checks {counts['positive_time_samples']} exact positive-time evolved polynomials, all with every root real. The independent checker reconstructs separate rational models and the forward/backward local polynomials.

```text
python work/rh_compute/scripts/jensen_window_pf_newman_regular_field_first_loss_nonclosure_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_regular_field_first_loss_nonclosure_gate.py
```

Current result:

```text
validated regular-field first-loss nonclosure gate: 16 rows, 3 sources, 10 arbitrary-field models, 15 multiplicities, 10 positive-time samples, 15 forward Hermite checks, 15 backward nonreal checks, 2 open handoffs, 0 Xi field bounds, 0 degree-uniform bounds
```

No `pi` enters this construction. Its constants are rational functions of the prescribed field; any `pi` in the Xi programme remains tied to the separate Fourier/theta normalization.
"""


def main() -> None:
    load_sources()
    certificate = build_symbolic_certificate()
    models = build_model_audit()
    hermite = build_hermite_audit()
    rows = build_rows(certificate)
    artifact = {
        "kind": KIND,
        "date": "2026-08-03",
        "status": "exact arbitrary-field first-loss heat-flow nonclosure theorem",
        "proof_boundary": (
            "This gate proves that every real regular field occurs in an even, "
            "negative-rooted positive-coefficient finite genus-zero model with an "
            "exact Newman-style first-loss heat boundary, for every multiplicity "
            "m>=2. It proves no bound for the actual Xi field, no degree-uniform "
            "Jensen remainder, no Xi collision exclusion, Lambda<=0, RH, or "
            "prize-level conclusion."
        ),
        "source_audit": source_audit(),
        "symbolic_certificate": certificate,
        "model_audit": models,
        "hermite_audit": hermite,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "sources": len(SOURCE_PATHS),
            "field_models": len(models),
            "multiplicities": len(hermite),
            "positive_coefficient_checks": len(models),
            "negative_root_checks": len(models),
            "signed_field_checks": len(models),
            "positive_time_samples": len(models),
            "forward_hermite_checks": len(hermite),
            "backward_nonreal_checks": len(hermite),
            "open_handoffs": 2,
            "actual_xi_field_bounds": 0,
            "degree_uniform_bounds": 0,
            "rh_conclusions": 0,
        },
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    counts = artifact["counts"]
    print(
        "built regular-field first-loss nonclosure gate: "
        f"{counts['rows']} rows, {counts['sources']} sources, "
        f"{counts['field_models']} arbitrary-field models, "
        f"{counts['multiplicities']} multiplicities, "
        f"{counts['positive_time_samples']} positive-time samples, "
        f"{counts['forward_hermite_checks']} forward Hermite checks, "
        f"{counts['backward_nonreal_checks']} backward nonreal checks, "
        f"{counts['open_handoffs']} open handoffs, "
        f"{counts['actual_xi_field_bounds']} Xi field bounds, "
        f"{counts['degree_uniform_bounds']} degree-uniform bounds"
    )


if __name__ == "__main__":
    main()
