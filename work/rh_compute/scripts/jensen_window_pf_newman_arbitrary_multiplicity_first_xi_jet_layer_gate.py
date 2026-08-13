#!/usr/bin/env python3
"""Build the complete first-Xi-jet correction to the cubic Jensen layer."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_arbitrary_multiplicity_first_xi_jet_layer_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"
MIN_MULTIPLICITY = 3
MAX_AUDIT_MULTIPLICITY = 16

SOURCE_PATHS = {
    "tertiary_layer": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_arbitrary_multiplicity_tertiary_universal_layer_gate.json",
    "positive_boundary": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_positive_boundary_attainment_lemma.json",
    "first_jet": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_first_jet_winding_gate.json",
}


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

    tertiary = payloads["tertiary_layer"].get("symbolic_certificate", {})
    require(
        tertiary.get("tertiary_operator")
        == "Q_(m,n)=((r-y/2)*partial_y^2+(2n+1)*partial_y/4)P_m.",
        "tertiary operator drift",
    )
    require(
        "rho*[U'(rho)/U(rho)]*P_(m+1)" in tertiary.get("unit_delay", ""),
        "tertiary unit handoff drift",
    )
    require(
        "If Lambda>0" in payloads["positive_boundary"].get("exact", {}).get("attainment_theorem", ""),
        "positive-boundary attainment drift",
    )
    require(
        "floor(m/2)" in payloads["first_jet"].get("exact", {}).get("multiplicity_index", ""),
        "first-jet multiplicity-index drift",
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


def secondary_polynomial(m: int, y: sp.Symbol, q: sp.Expr) -> sp.Expr:
    result = sp.Integer(0)
    for heat_order in range(m // 2 + 1):
        for cubic_order in range((m - 2 * heat_order) // 3 + 1):
            degree = m - 2 * heat_order - 3 * cubic_order
            result += (
                sp.factorial(m)
                * q**heat_order
                * y**degree
                / (
                    sp.factorial(heat_order)
                    * sp.factorial(cubic_order)
                    * 12**cubic_order
                    * sp.factorial(degree)
                )
            )
    return sp.expand(result)


def tertiary_operator(
    polynomial: sp.Expr,
    n: int | sp.Expr,
    y: sp.Symbol,
    r: sp.Expr,
) -> sp.Expr:
    shift = sp.sympify(2 * n + 1)
    return sp.expand(
        (r - y / 2) * sp.diff(polynomial, y, 2)
        + shift * sp.diff(polynomial, y) / 4
    )


def quaternary_background(
    polynomial: sp.Expr,
    n: int | sp.Expr,
    y: sp.Symbol,
    q: sp.Expr,
    r: sp.Expr,
    v: sp.Expr,
) -> sp.Expr:
    tertiary = tertiary_operator(polynomial, n, y, r)
    shift = sp.sympify(2 * n + 1)
    return sp.expand(
        v * sp.diff(polynomial, y, 2)
        + tertiary_operator(tertiary, n, y, r) / 2
        + q * shift * sp.diff(polynomial, y) / 2
        + (q * y + shift / 8) * sp.diff(polynomial, y, 2)
        + (q**2 + r / 2) * sp.diff(polynomial, y, 3)
        + q * sp.diff(polynomial, y, 4) / 4
        + sp.diff(polynomial, y, 5) / 80
    )


def complete_quaternary(
    m: int,
    n: int | sp.Expr,
    y: sp.Symbol,
    q: sp.Expr,
    r: sp.Expr,
    v: sp.Expr,
    rho: sp.Expr,
    u1: sp.Expr,
) -> sp.Expr:
    primary = secondary_polynomial(m, y, q)
    return sp.expand(
        quaternary_background(primary, n, y, q, r, v)
        + rho * u1 * secondary_polynomial(m + 1, y, q)
    )


def build_operator_audit() -> list[dict]:
    y, q, r, v = sp.symbols("y q r v")
    rows: list[dict] = []
    for m in range(MIN_MULTIPLICITY, MAX_AUDIT_MULTIPLICITY + 1):
        primary = secondary_polynomial(m, y, q)
        require(
            sp.expand(sp.diff(primary, q) - sp.diff(primary, y, 2)) == 0,
            f"heat Appell identity drift m={m}",
        )
        require(
            sp.expand(
                sp.diff(secondary_polynomial(m + 1, y, q), y) - (m + 1) * primary
            )
            == 0,
            f"degree Appell identity drift m={m}",
        )
        shift_rows = []
        for n in (0, 3):
            tertiary = tertiary_operator(primary, n, y, r)
            background = quaternary_background(primary, n, y, q, r, v)
            shift_rows.append(
                {
                    "shift": n,
                    "tertiary_sha256": expression_hash(tertiary),
                    "quaternary_background_sha256": expression_hash(background),
                    "quaternary_degree": sp.Poly(background, y).degree(),
                    "quaternary_terms": len(sp.Poly(background, y, q, r, v).terms()),
                }
            )
        rows.append(
            {
                "multiplicity": m,
                "primary_sha256": expression_hash(primary),
                "primary_degree": sp.Poly(primary, y).degree(),
                "appell_heat": "partial_q P_m=partial_y^2 P_m",
                "appell_degree": "partial_y P_(m+1)=(m+1)P_m",
                "shift_audit": shift_rows,
            }
        )
    return rows


def exact_finite_model(m: int, n: int) -> dict:
    eps, x, z, y = sp.symbols("eps x z y")
    rho, q, r, v, u1 = sp.symbols("rho q r v u1")
    delta = (
        rho * eps**6 / 8
        + rho * q * eps**8 / 4
        + rho * r * eps**10 / 4
        + rho * v * eps**12 / 4
    )
    source = (x - rho) ** m * (1 + u1 * (x - rho))

    def generator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            4 * x * sp.diff(polynomial, x, 2)
            + (4 * n + 2) * sp.diff(polynomial, x)
        )

    evolved = sp.Integer(0)
    iterate = source
    for order in range(m + 2):
        evolved += delta**order * iterate / sp.factorial(order)
        iterate = generator(iterate)

    transformed = sp.Integer(0)
    for (degree,), coefficient in sp.Poly(sp.expand(evolved), x).terms():
        ratio = sp.prod(1 - index * eps**6 for index in range(degree))
        transformed += coefficient * ratio * z**degree
    local = sp.expand(transformed.subs(z, rho + rho * y * eps**4))

    lead = sp.expand(local.coeff(eps, 4 * m) / rho**m)
    odd_one = sp.expand(local.coeff(eps, 4 * m + 1) / rho**m)
    tertiary = sp.expand(local.coeff(eps, 4 * m + 2) / rho**m)
    odd_three = sp.expand(local.coeff(eps, 4 * m + 3) / rho**m)
    quaternary = sp.expand(local.coeff(eps, 4 * m + 4) / rho**m)
    odd_five = sp.expand(local.coeff(eps, 4 * m + 5) / rho**m)

    target_lead = secondary_polynomial(m, y, q)
    target_tertiary = tertiary_operator(target_lead, n, y, r)
    target_quaternary = complete_quaternary(m, n, y, q, r, v, rho, u1)
    require(sp.expand(lead - target_lead) == 0, f"lead drift m={m}, n={n}")
    require(odd_one == 0, f"epsilon-one drift m={m}, n={n}")
    require(sp.expand(tertiary - target_tertiary) == 0, f"tertiary drift m={m}, n={n}")
    require(odd_three == 0, f"epsilon-three drift m={m}, n={n}")
    require(sp.expand(quaternary - target_quaternary) == 0, f"quaternary drift m={m}, n={n}")
    require(odd_five == 0, f"epsilon-five drift m={m}, n={n}")

    unit_coefficient = sp.expand(quaternary.coeff(u1))
    unit_target = sp.expand(rho * secondary_polynomial(m + 1, y, q))
    require(sp.expand(unit_coefficient - unit_target) == 0, f"unit drift m={m}, n={n}")
    return {
        "multiplicity": m,
        "shift": n,
        "lead_sha256": expression_hash(lead),
        "tertiary_sha256": expression_hash(tertiary),
        "quaternary_sha256": expression_hash(quaternary),
        "unit_sha256": expression_hash(unit_coefficient),
        "epsilon_one_zero": True,
        "epsilon_three_zero": True,
        "epsilon_five_zero": True,
        "matched": True,
    }


def build_exact_finite_audit() -> list[dict]:
    return [exact_finite_model(m, n) for n in (0, 3) for m in range(3, 9)]


def contact_correction() -> dict[str, str]:
    y, q, r, v, n, ell, t = sp.symbols("y q r v n ell t")
    p3 = secondary_polynomial(3, y, q)
    p4 = secondary_polynomial(4, y, q)
    q3 = -2 ** sp.Rational(-7, 3)
    y3 = 2 ** sp.Rational(-2, 3)
    r3 = y3 / 2
    s3 = (1 - 2 * n) / 4
    tertiary = tertiary_operator(p3, n, y, r)
    quaternary = quaternary_background(p3, n, y, q, r, v) + ell * p4

    value_equation = sp.simplify(
        s3**2 * sp.diff(p3, y, 2) / 2
        + s3 * sp.diff(tertiary, y)
        + quaternary
    ).subs({q: q3, y: y3, r: r3})
    v3 = sp.solve(sp.simplify(value_equation), v)[0]
    derivative_equation = sp.simplify(
        t * sp.diff(p3, y, 2)
        + s3**2 * sp.diff(p3, y, 3) / 2
        + s3 * sp.diff(tertiary, y, 2)
        + sp.diff(quaternary, y)
    ).subs({q: q3, y: y3, r: r3, v: v3})
    t3 = sp.solve(sp.simplify(derivative_equation), t)[0]
    require(sp.simplify(v3 + ell / 4 + 3 * n / 8 + sp.Rational(1, 4)) == 0, "m3 v drift")
    require(
        sp.simplify(t3 - 2 ** sp.Rational(2, 3) * (2 * ell + n + 4) / 8) == 0,
        "m3 root drift",
    )
    return {
        "local_logarithmic_jet": "ell=rho*U'(rho)/U(rho)",
        "m3_heat_coefficient": "v_3=-(2*ell+3*n+2)/8",
        "m3_root_coefficient": "t_3=2^(-7/3)*(2*ell+n+4)",
        "m3_heat_expansion": "lambda_D=lambda_*+rho/(8D)-rho*2^(-13/3)D^(-4/3)+rho*2^(-11/3)D^(-5/3)-rho*(2*ell+3*n+2)D^(-2)/32+o(D^(-2))",
        "m3_root_expansion": "w_D=rho/D+rho*2^(-2/3)D^(-5/3)+rho*(1-2*n)D^(-2)/4+rho*2^(-7/3)*(2*ell+n+4)D^(-7/3)+o(D^(-7/3))",
        "newman_shift_zero": "For rho=-c^2 and n=0: t_D=Lambda-c^2/(8D)+c^2*2^(-13/3)D^(-4/3)-c^2*2^(-11/3)D^(-5/3)+c^2*(ell+1)D^(-2)/16+o(D^(-2)).",
    }


def canonical_product_audit() -> dict:
    s = sp.symbols("s")
    rho = sp.Integer(-1)
    models = []
    for other_root in (sp.Integer(-4), sp.Rational(-1, 4)):
        source = sp.expand((s - rho) ** 3 * (s - other_root))
        unit = s - other_root
        logarithmic_jet = sp.simplify(sp.diff(unit, s).subs(s, rho) / unit.subs(s, rho))
        derivative_jet = sp.simplify(
            sp.diff(source, s, 4).subs(s, rho)
            / (4 * sp.diff(source, s, 3).subs(s, rho))
        )
        require(sp.simplify(logarithmic_jet - derivative_jet) == 0, "derivative jet identity drift")
        ell = sp.simplify(rho * logarithmic_jet)
        models.append(
            {
                "multiple_root": str(rho),
                "other_root": str(other_root),
                "polynomial": str(source),
                "unit_logarithmic_derivative": str(logarithmic_jet),
                "ell": str(ell),
                "all_roots_negative": True,
                "all_coefficients_positive": all(coefficient > 0 for coefficient in sp.Poly(source, s).all_coeffs()),
            }
        )
    require(models[0]["ell"] == "-1/3", "negative ell witness drift")
    require(models[1]["ell"] == "4/3", "positive ell witness drift")
    return {
        "derivative_identity": "U'(rho)/U(rho)=F^(m+1)(rho)/((m+1)F^(m)(rho))",
        "genus_zero_identity": "If F(s)=F(0)prod_j(1-s/rho_j) with multiplicities mu_j, then ell_k=rho_k*sum_(j!=k)mu_j/(rho_k-rho_j).",
        "models": models,
        "verdict": "Negative-rootedness and positive coefficients do not fix the sign of ell at an interior multiple root.",
    }


def symbolic_certificate() -> dict[str, str]:
    return {
        "extended_scaling": "Set epsilon=D^(-1/6), lambda=lambda_*+rho*epsilon^6/8+(rho*q/4)*epsilon^8+(rho*r/4)*epsilon^10+(rho*v/4)*epsilon^12, and evaluate at (rho+rho*y*epsilon^4)/D.",
        "full_expansion": "After normalization by rho^m U(rho) epsilon^(4m), the cluster is P_m+epsilon^2 A_(n,r)P_m+epsilon^4{R_(n,q,r,v)P_m+ell P_(m+1)}+O(epsilon^6), where ell=rho U'(rho)/U(rho).",
        "tertiary": "A_(n,r)=(r-y/2)*partial_y^2+(2n+1)*partial_y/4.",
        "quaternary": "R P=v*P''+(A^2 P)/2+q*(2n+1)P'/2+(q*y+(2n+1)/8)P''+(q^2+r/2)P'''+q*P''''/4+P'''''/80.",
        "monomial_frame": "Conjugating the remainder by exp(-q*partial_y^2-partial_y^3/12) gives [q*(2n+1)partial_y/2+(q*y+(2n+1)/8)partial_y^2+(r/2-q^2)partial_y^3+partial_y^5/80]y^m; conjugating y back adds 2q partial_y+partial_y^2/4 and yields R.",
        "unit_entry": "The complete source-specific contribution at relative epsilon^4 is ell P_(m+1), with ell=rho U'(rho)/U(rho); U'' first enters later.",
        "m3": "At q_3=-2^(-7/3), y_3=2^(-2/3), r_3=2^(-5/3), and s_3=(1-2n)/4, the epsilon^4 contact equations give v_3=-(2ell+3n+2)/8 and t_3=2^(-7/3)(2ell+n+4).",
        "sign_guard": "The genus-zero product writes ell as a signed neighboring-root sum. Exact all-negative-root, positive-coefficient polynomial models realize ell=-1/3 and ell=4/3, so local Laguerre-Polya geometry alone supplies no sign.",
        "actual_xi_handoff": "For a hypothetical multiple zero of H_Lambda, ell is exact but depends on the location of that zero relative to the remaining H_Lambda zeros. Positive-boundary attainment does not bound this regular field.",
        "uniform_guard": "The fixed-m local expansion remains noncontradictory. Collision exclusion still requires a remainder uniform in the selected Jensen degree or a separate global first-jet winding theorem.",
    }


def build_rows(certificate: dict[str, str], contact: dict[str, str], product: dict) -> list[GateRow]:
    rows = [
        GateRow("fxj_01_sources", "source chain", "proved", "Three parent contracts are current and hash-pinned.", "The tertiary, positive-boundary, and first-jet gates are audited.", "No parent theorem is strengthened."),
        GateRow("fxj_02_scaling", "first-jet scaling", "proved", "The next heat parameter enters at D^(-2).", certificate["extended_scaling"], "This is a local selected-sequence expansion."),
        GateRow("fxj_03_expansion", "complete epsilon-four layer", "proved", "The full relative epsilon^4 coefficient is explicit.", certificate["full_expansion"], "The remainder is fixed-m, not degree-uniform."),
        GateRow("fxj_04_tertiary", "parent transport", "proved", "The epsilon^2 transport is retained exactly.", certificate["tertiary"], "This row imports rather than reproves the parent layer."),
        GateRow("fxj_05_operator", "quaternary operator", "proved", "The universal epsilon^4 background is a finite fifth-order operator.", certificate["quaternary"], "The displayed A^2 is operator composition."),
        GateRow("fxj_06_normal", "normal-order certificate", "proved", "The operator follows from one finite Weyl conjugation.", certificate["monomial_frame"], "This is an identity for every fixed multiplicity."),
        GateRow("fxj_07_unit", "first Xi unit", "proved", "The complete source-specific term at this order is isolated.", certificate["unit_entry"], "Its coefficient is exact but has no universal sign."),
        GateRow("fxj_08_appell", "Appell closure", "proved", "All terms remain in the P_m derivative family.", "partial_q P_m=partial_y^2 P_m and partial_y P_(m+1)=(m+1)P_m.", "Appell closure is not hyperbolicity preservation."),
        GateRow("fxj_09_finite", "exact finite audit", "guard_validated", "Exact radial-heat/Jensen models reproduce the whole layer.", "Models include free rho,q,r,v,u1 at shifts zero and three.", "Finite audits check algebra, not degree-uniformity."),
        GateRow("fxj_10_contact", "contact perturbation", "proved", "The epsilon^4 value and derivative equations are retained.", contact["m3_heat_coefficient"] + "; " + contact["m3_root_coefficient"], "The explicit solution is for m=3."),
        GateRow("fxj_11_newman", "Newman refinement", "proved", "The conditional unshifted collision time is explicit through D^(-2).", contact["newman_shift_zero"], "This refines rather than excludes the collision."),
        GateRow("fxj_12_product", "canonical product", "proved", "The local Xi jet is a regular neighboring-zero field.", product["genus_zero_identity"], "The identity requires the genus-zero s-coordinate product."),
        GateRow("fxj_13_countermodels", "sign countermodels", "guard_validated", "Local negative-root geometry does not determine the sign of ell.", product["verdict"], "The models are proof-safety witnesses, not Xi data."),
        GateRow("fxj_14_actual_xi", "actual Xi field", "open", "The hypothetical boundary field ell must be bounded from Xi structure.", certificate["actual_xi_handoff"], "Existence of a finite multiple zero does not locate its neighboring-root sum."),
        GateRow("fxj_15_uniform", "uniform handoff", "open", "A selected-degree-uniform remainder is still required.", certificate["uniform_guard"], "Fixed-m asymptotics do not exclude the cofinal sequence."),
        GateRow("fxj_16_boundary", "proof boundary", "guard_validated", "This gate is not a proof of RH.", "It proves the complete first-local-jet operator and exact m=3 displacement.", "No Xi field bound, uniform collision exclusion, Lambda<=0, RH, or prize-level conclusion is proved."),
    ]
    require(len(rows) == 16, "gate row count drifted")
    require(sum(row.readiness == "open" for row in rows) == 2, "open row count drifted")
    return rows


def render_note(artifact: dict) -> str:
    cert = artifact["symbolic_certificate"]
    contact = artifact["contact_correction"]
    product = artifact["canonical_product_audit"]
    counts = artifact["counts"]
    return f"""# Newman Arbitrary-Multiplicity First Xi-Jet Layer Gate

Date: 2026-08-03

Status: complete fixed-m relative `epsilon^4` operator and exact multiplicity-three Xi-jet displacement; actual Xi field bound and degree-uniform exclusion open; not a proof of RH.

## Complete Expansion

{cert['extended_scaling']}

{cert['full_expansion']}

The tertiary operator is

```text
{cert['tertiary']}
```

The complete universal background is

```text
{cert['quaternary']}
```

{cert['monomial_frame']}

## First Local Xi Jet

{cert['unit_entry']}

For multiplicity three, put `ell=rho U'(rho)/U(rho)`. The next contact coefficients are

```text
{contact['m3_heat_coefficient']}
{contact['m3_root_coefficient']}
{contact['m3_heat_expansion']}
{contact['m3_root_expansion']}
```

At `n=0` and `rho=-c^2`,

```text
{contact['newman_shift_zero']}
```

## Canonical-Product Field

{product['derivative_identity']}

{product['genus_zero_identity']}

{product['verdict']} The two exact witnesses have `ell=-1/3` and `ell=4/3` while retaining only negative roots and positive coefficients.

Therefore the new local term moves the collision but does not contradict it without additional Xi information.

## Exact Audit

The builder checks {counts['operator_checks']} symbolic shift/operator instances over {counts['multiplicities']} multiplicities and {counts['exact_finite_models']} exact finite radial-heat/Jensen models with free `rho,q,r,v,u1`. The independent checker specializes every finite model to separate rational data and reconstructs the operator without importing this builder.

## Live Handoff

Bound the actual regular field `ell` at every hypothetical positive-boundary multiple zero using the Xi canonical product, zero dynamics, or the Fourier kernel. In parallel, derive a local-uniform remainder with constants uniform along the selected degree sequence. The independent first-jet winding theorem remains the determinant-free alternative.

## Pi Provenance

No pi enters this local operator calculation. The rational `1/80` is the fifth-order finite-Jensen/heat normal-order coefficient. Pi enters only in separate Fourier normalizations of the Xi kernel.

## Proof Boundary

{cert['uniform_guard']} This gate proves no Xi regular-field sign, degree-uniform collision exclusion, all-degree Jensen hyperbolicity, PF-infinity, `Lambda<=0`, RH, or prize-level conclusion.
"""


def main() -> None:
    load_sources()
    operator_rows = build_operator_audit()
    finite_rows = build_exact_finite_audit()
    contact = contact_correction()
    product = canonical_product_audit()
    certificate = symbolic_certificate()
    rows = build_rows(certificate, contact, product)
    artifact = {
        "kind": KIND,
        "date": "2026-08-03",
        "status": "complete fixed-m first Xi-jet layer; actual Xi field and uniform exclusion open",
        "proof_boundary": rows[-1].proof_boundary,
        "source_audit": source_audit(),
        "symbolic_certificate": certificate,
        "contact_correction": contact,
        "canonical_product_audit": product,
        "operator_audit": operator_rows,
        "exact_finite_model_audit": finite_rows,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "sources": len(SOURCE_PATHS),
            "multiplicities": len(operator_rows),
            "operator_checks": 2 * len(operator_rows),
            "exact_finite_models": len(finite_rows),
            "shift_values": 2,
            "appell_identities": 2 * len(operator_rows),
            "full_xi_jet_checks": len(finite_rows),
            "exact_m3_contact_corrections": 2,
            "canonical_product_identities": 2,
            "sign_countermodels": len(product["models"]),
            "open_handoffs": sum(row.readiness == "open" for row in rows),
            "degree_uniform_bounds": 0,
            "rh_conclusions": 0,
        },
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    counts = artifact["counts"]
    print(
        "built first Xi-jet Jensen layer gate: "
        f"{counts['rows']} rows, {counts['sources']} sources, "
        f"{counts['multiplicities']} multiplicities, "
        f"{counts['operator_checks']} operator checks, "
        f"{counts['exact_finite_models']} exact finite models, "
        f"{counts['appell_identities']} Appell identities, "
        f"{counts['full_xi_jet_checks']} full Xi-jet checks, "
        f"{counts['exact_m3_contact_corrections']} exact m3 corrections, "
        f"{counts['sign_countermodels']} sign countermodels, "
        f"{counts['open_handoffs']} open handoffs, "
        f"{counts['degree_uniform_bounds']} degree-uniform bounds"
    )


if __name__ == "__main__":
    main()
