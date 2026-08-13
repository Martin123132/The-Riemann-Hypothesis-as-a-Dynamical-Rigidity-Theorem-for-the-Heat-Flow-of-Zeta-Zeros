#!/usr/bin/env python3
"""Build the arbitrary-multiplicity secondary cubic Jensen layer gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_arbitrary_multiplicity_secondary_cubic_layer_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"
MIN_MULTIPLICITY = 3
MAX_AUDIT_MULTIPLICITY = 16

SOURCE_PATHS = {
    "primary_layer": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_arbitrary_multiplicity_jensen_boundary_layer_gate.json",
    "double_layer": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_scaled_double_zero_boundary_layer_lemma.json",
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

    primary = payloads["primary_layer"].get("symbolic_certificate", {})
    require(
        primary.get("coefficient_formula")
        == "K_(m,a)(eta)=m!*sum_(q=0)^floor(m/2) a^q*eta^(m-2q)/(q!*(m-2q)!).",
        "primary arbitrary-multiplicity layer drift",
    )
    require(
        primary.get("collision_scale", "").startswith(
            "The layer changes type only at a=0, hence tau=rho/8."
        ),
        "primary collision scale drift",
    )
    require(
        payloads["double_layer"].get("exact", {}).get("second_correction")
        == "J_D(z/D)=F-z^2*F''/(2*D)+(z^3*F'''/3+z^4*F''''/8)/D^2+O_K(D^(-3))",
        "finite-Jensen second correction drift",
    )
    require(
        payloads["positive_boundary"].get("exact", {}).get("attainment_theorem", "").startswith(
            "If Lambda>0, H_Lambda has a finite real multiple zero c"
        ),
        "positive-boundary source drift",
    )
    require(
        payloads["first_jet"].get("exact", {}).get("multiplicity_index")
        == "At a real zero of spatial multiplicity m>=2 in any nonzero real-analytic solution F_t=-F_xx, the common-zero point is isolated and its local Brouwer index in the standard (x,t) orientation is +floor(m/2)",
        "first-jet multiplicity source drift",
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


def secondary_polynomial(m: int, y: sp.Symbol, q: sp.Symbol | sp.Rational) -> sp.Expr:
    result = sp.Integer(0)
    for r in range(m // 2 + 1):
        for s in range((m - 2 * r) // 3 + 1):
            degree = m - 2 * r - 3 * s
            result += (
                sp.factorial(m)
                * q**r
                * y**degree
                / (
                    sp.factorial(r)
                    * sp.factorial(s)
                    * 12**s
                    * sp.factorial(degree)
                )
            )
    return sp.expand(result)


def build_coefficient_audit() -> tuple[list[dict], int]:
    y, q = sp.symbols("y q")
    rows: list[dict] = []
    checks = 0
    for m in range(MIN_MULTIPLICITY, MAX_AUDIT_MULTIPLICITY + 1):
        polynomial = secondary_polynomial(m, y, q)
        coefficients = []
        for r in range(m // 2 + 1):
            for s in range((m - 2 * r) // 3 + 1):
                degree = m - 2 * r - 3 * s
                target = sp.factorial(m) / (
                    sp.factorial(r)
                    * sp.factorial(s)
                    * 12**s
                    * sp.factorial(degree)
                )
                observed = polynomial.coeff(y, degree).coeff(q, r)
                require(
                    observed == target,
                    f"secondary coefficient drift at m={m}, r={r}, s={s}",
                )
                coefficients.append(
                    {
                        "heat_index": r,
                        "cubic_index": s,
                        "y_power": degree,
                        "coefficient": str(target),
                    }
                )
                checks += 1

        require(polynomial.coeff(y, m - 1) == 0, f"y^(m-1) drift at m={m}")
        require(
            polynomial.subs(q, 0).coeff(y, m - 2) == 0,
            f"q=0 y^(m-2) drift at m={m}",
        )
        require(
            polynomial.subs(q, 0).coeff(y, m - 3)
            == sp.factorial(m) / (12 * sp.factorial(m - 3)),
            f"q=0 cubic coefficient drift at m={m}",
        )
        rows.append(
            {
                "multiplicity": m,
                "term_count": len(coefficients),
                "polynomial": str(polynomial),
                "coefficients": coefficients,
                "q_zero_hyperbolic": False,
                "threshold_location": "q_m<0",
            }
        )
    return rows, checks


def finite_jensen_transform(polynomial: sp.Expr, x: sp.Symbol, z: sp.Symbol, eps: sp.Symbol) -> sp.Expr:
    transformed = sp.Integer(0)
    for (degree,), coefficient in sp.Poly(sp.expand(polynomial), x).terms():
        falling_ratio = sp.prod(1 - index * eps**6 for index in range(degree))
        transformed += coefficient * falling_ratio * z**degree
    return sp.expand(transformed)


def exact_secondary_model(m: int, n: int) -> dict:
    eps, x, z, y = sp.symbols("eps x z y")
    rho, q, u1 = sp.symbols("rho q u1")
    kappa = rho * q / 4
    delta = rho * eps**6 / 8 + kappa * eps**8
    source = (x - rho) ** m * (1 + u1 * (x - rho))

    def generator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            4 * x * sp.diff(polynomial, x, 2)
            + (4 * n + 2) * sp.diff(polynomial, x)
        )

    heat_polynomial = sp.Integer(0)
    iterate = source
    for order in range(m + 2):
        heat_polynomial += delta**order * iterate / sp.factorial(order)
        iterate = generator(iterate)

    transformed = finite_jensen_transform(heat_polynomial, x, z, eps)
    local = sp.expand(transformed.subs(z, rho + rho * y * eps**4))
    for power in range(4 * m):
        require(
            sp.expand(local.coeff(eps, power)) == 0,
            f"uncancelled lower valuation at m={m}, n={n}, eps^{power}",
        )
    exact_limit = sp.expand(local.coeff(eps, 4 * m))
    target = sp.expand(rho**m * secondary_polynomial(m, y, q))
    require(
        sp.expand(exact_limit - target) == 0,
        f"secondary exact-model drift at m={m}, n={n}",
    )
    require(not exact_limit.has(u1), f"local-unit leakage at m={m}, n={n}")
    return {
        "multiplicity": m,
        "shift": n,
        "limit": str(exact_limit),
        "target": str(target),
        "local_unit_derivative_absent": True,
        "matched": True,
    }


def build_exact_model_audit() -> list[dict]:
    return [exact_secondary_model(m, n) for n in (0, 3) for m in range(3, 9)]


def exact_cubic_threshold() -> dict[str, str]:
    y, q = sp.symbols("y q")
    polynomial = secondary_polynomial(3, y, q)
    require(polynomial == y**3 + 6 * q * y + sp.Rational(1, 2), "m=3 polynomial drift")
    discriminant = sp.factor(sp.discriminant(polynomial, y))
    require(
        sp.expand(discriminant + sp.Rational(27, 4) * (128 * q**3 + 1)) == 0,
        "m=3 discriminant drift",
    )
    y_star = 2 ** sp.Rational(-2, 3)
    q_star = -2 ** sp.Rational(-7, 3)
    require(sp.simplify(polynomial.subs({y: y_star, q: q_star})) == 0, "m=3 root drift")
    require(
        sp.simplify(sp.diff(polynomial, y).subs({y: y_star, q: q_star})) == 0,
        "m=3 derivative drift",
    )
    return {
        "polynomial": "P_3(y,q)=y^3+6*q*y+1/2",
        "discriminant": "Disc_y P_3=-(27/4)*(128*q^3+1)",
        "threshold": "q_3=-2^(-7/3)",
        "double_root": "y_3=2^(-2/3)",
        "parameter_correction": "(rho*q_3/4)*D^(-4/3)=-rho*2^(-13/3)*D^(-4/3)",
    }


def symbolic_certificate() -> dict[str, str]:
    return {
        "secondary_scaling": "For epsilon=D^(-1/6), set lambda=lambda_*+rho*epsilon^6/8+(rho*q/4)*epsilon^8 and evaluate J_(D,n,lambda) at (rho+rho*y*epsilon^4)/D.",
        "secondary_limit": "If F_(n,lambda_*)(z)=(z-rho)^m U(z), rho<0, m>=3, and U(rho)!=0, then epsilon^(-4m) J_(D,n,lambda)((rho+rho*y*epsilon^4)/D)/(rho^m U(rho)) tends locally uniformly to P_m(y,q)=exp(q*partial_y^2+partial_y^3/12)y^m.",
        "coefficient_formula": "P_m(y,q)=m!*sum_(2r+3s<=m) q^r*y^(m-2r-3s)/(r!*s!*12^s*(m-2r-3s)!).",
        "cubic_origin": "At tau=rho/8 the order-D^(-1/2) combined heat/Jensen logarithm has surviving cubic coefficient rho^3/12; the eta=D^(-1/6)*rho*y zoom exponentiates all powers of that cubic term.",
        "universality": "The fixed shift n and the logarithmic derivatives of U enter only below this secondary limit. The first Xi-specific local unit derivative is delayed beyond the universal cubic layer.",
        "hyperbolicity_ray": "For s>=0, P_m(y,q-s)=exp(-s*partial_y^2)P_m(y,q), and backward polynomial heat preserves real-rootedness. Hence the hyperbolicity set in q is a closed lower ray.",
        "negative_threshold": "For q tending to -infinity, the sqrt(-2q) rescaling converges to the simple-root probabilists' Hermite polynomial He_m. At q=0 the y^(m-1) and y^(m-2) coefficients vanish while y^(m-3) is nonzero, which is impossible for an all-real root multiset. Therefore the ray is (-infinity,q_m] for a finite q_m<0.",
        "collision_refinement": "For m>=3, local root continuity gives a finite Jensen collision sequence with lambda_D=lambda_*+rho/(8D)+(rho*q_m/4)D^(-4/3)+o(D^(-4/3)) and w_D=rho/D+rho*y_m*D^(-5/3)+o(D^(-5/3)), where y_m is a multiple root of P_m(.,q_m).",
        "newman_refinement": "Under Lambda>0, rho=-c^2, so the D^(-4/3) parameter correction is positive because rho<0 and q_m<0. The secondary layer refines, but does not contradict, the fixed-shift collision sequence.",
        "double_zero_guard": "For m=2 the cubic operator annihilates y^2, so this D^(-4/3) layer is absent; the existing D^(-2) external-field correction remains the correct double-zero refinement.",
        "first_jet_guard": "Every local contact still has positive first-jet index floor(mu/2). The universal cubic layer determines local cluster collisions but supplies no boundary winding cancellation or Xi edge estimate.",
        "open_target": "The attempted first correction is universal rather than Xi-specific. The next local calculation is the D^(-1/3) correction to P_m; the local unit U first appears at D^(-2/3). A contradiction still requires a degree-uniform Xi estimate or the independent global first-jet boundary theorem.",
    }


def build_rows(certificate: dict[str, str], m3: dict[str, str]) -> list[GateRow]:
    rows = [
        GateRow("scb_01_sources", "source chain", "proved", "Four parent contracts are current and hash-pinned.", "The primary arbitrary-multiplicity layer, double-zero correction, positive-boundary theorem, and first-jet index gate are audited.", "No parent conclusion is strengthened."),
        GateRow("scb_02_scaling", "critical scaling", "proved", "Multiplicity at least three has a second critical zoom.", certificate["secondary_scaling"], "The double-zero case has a different next scale."),
        GateRow("scb_03_limit", "secondary limit", "proved", "The complete cubic secondary polynomial is locally universal.", certificate["secondary_limit"], "The convergence is local in the secondary coordinate."),
        GateRow("scb_04_coefficients", "exact coefficients", "proved", "Every coefficient of the secondary polynomial is explicit.", certificate["coefficient_formula"], "The formula is symbolic in m; the finite table is an audit."),
        GateRow("scb_05_cubic_origin", "graded operator", "proved", "The rho^3/12 cubic residual survives the cancelled Gaussian layer.", certificate["cubic_origin"], "Terms of higher graded weight do not enter the secondary limit."),
        GateRow("scb_06_universality", "Xi guard", "guard_validated", "Neither the fixed shift nor the local Xi unit enters this layer.", certificate["universality"], "No Xi-specific sign can be extracted at this order."),
        GateRow("scb_07_heat_ray", "hyperbolicity monotonicity", "proved", "Secondary hyperbolicity is a lower ray in q.", certificate["hyperbolicity_ray"], "This is polynomial backward-heat preservation, not an all-degree Xi theorem."),
        GateRow("scb_08_threshold", "finite threshold", "proved", "Every m>=3 has a finite negative universal threshold q_m.", certificate["negative_threshold"], "No closed formula for general q_m is asserted."),
        GateRow("scb_09_collision", "collision refinement", "proved", "The finite collision has a D^(-4/3) parameter correction and D^(-5/3) root correction.", certificate["collision_refinement"], "A multiple root y_m is selected from the universal threshold polynomial."),
        GateRow("scb_10_m3", "cubic exact model", "proved", "The multiplicity-three threshold and contact are exact radicals.", f"{m3['polynomial']}; {m3['discriminant']}; {m3['threshold']}; {m3['double_root']}.", "This exact constant is specific to m=3."),
        GateRow("scb_11_newman", "positive-boundary consequence", "proved", "The universal correction approaches Lambda from above the primary rho/(8D) estimate.", certificate["newman_refinement"], "This remains conditional on Lambda>0 and is not a contradiction."),
        GateRow("scb_12_double", "multiplicity split", "guard_validated", "Multiplicity two retains its existing D^(-2) refinement.", certificate["double_zero_guard"], "Do not apply the cubic secondary scale to m=2."),
        GateRow("scb_13_first_jet", "route comparison", "guard_validated", "Local positive contact charge does not close global first-jet winding.", certificate["first_jet_guard"], "The global Xi edge theorem remains open."),
        GateRow("scb_14_tertiary", "next local order", "open", "The first possible post-cubic correction must be derived without assuming an Xi sign.", certificate["open_target"], "The D^(-1/3) term may still be universal."),
        GateRow("scb_15_uniform", "degree-uniform handoff", "open", "Local asymptotics still require degree-uniform Xi control to become a contradiction.", "Bounded-degree certificates and fixed-order Hankel tails do not exclude the refined sequence.", "No uniform cofinal estimate is supplied."),
        GateRow("scb_16_boundary", "proof boundary", "guard_validated", "This gate is not a proof of RH.", "It refines arbitrary-multiplicity finite Jensen collisions and proves that the first refined layer is universal.", "No collision exclusion, all-degree hyperbolicity, PF-infinity, Lambda<=0, RH, or prize-level conclusion is proved."),
    ]
    require(len(rows) == 16, "gate row count drifted")
    require(sum(row.readiness == "open" for row in rows) == 2, "open row count drifted")
    return rows


def render_note(artifact: dict) -> str:
    cert = artifact["symbolic_certificate"]
    m3 = artifact["multiplicity_three"]
    counts = artifact["counts"]
    return f"""# Newman Arbitrary-Multiplicity Secondary Cubic Layer Gate

Date: 2026-08-03

Status: exact universal secondary cubic Jensen layer for multiplicity at least three; Xi-specific and degree-uniform handoffs open; not a proof of RH.

## Secondary Scaling

{cert['secondary_scaling']}

{cert['secondary_limit']}

Equivalently,

```text
{cert['coefficient_formula']}
```

The exact audit checks {counts['universal_coefficient_checks']} symbolic coefficients for multiplicities `{MIN_MULTIPLICITY}` through `{MAX_AUDIT_MULTIPLICITY}` and {counts['exact_finite_models']} finite heat/Jensen models with nonconstant local unit at shifts zero and three.

## Why A Cubic Appears

{cert['cubic_origin']}

The primary Gaussian heat and finite-Jensen corrections cancel at `tau=rho/8`. On the smaller `D^(-1/6)` root coordinate, every power of the residual cubic derivative has the same graded size, so the cubic must be exponentiated rather than truncated.

## Universal Threshold

{cert['hyperbolicity_ray']}

{cert['negative_threshold']}

Write that threshold as `q_m`. Then

```text
{cert['collision_refinement']}
```

## Exact Multiplicity Three

```text
{m3['polynomial']}
{m3['discriminant']}
{m3['threshold']}
{m3['double_root']}
{m3['parameter_correction']}
```

## Newman Consequence

{cert['newman_refinement']}

This sharpens the location of the forced finite collision. It does not remove it.

## Xi And First-Jet Guards

{cert['universality']}

{cert['double_zero_guard']}

{cert['first_jet_guard']}

## Live Handoff

{cert['open_target']}

## Pi Provenance

No pi enters this universal local calculation. The coefficient `1/12` is the difference between the cubic finite-Jensen logarithm and the cubic commutator generated by radial heat at `tau=rho/8`; it is not a circular constant.

## Proof Boundary

This gate proves the secondary cubic local limit, its finite negative hyperbolicity threshold, the refined collision scales, and the exact `m=3` constant. It proves no Xi-specific exclusion, degree-uniform estimate, all-degree Jensen hyperbolicity, PF-infinity, `Lambda<=0`, RH, or prize-level conclusion.
"""


def main() -> None:
    load_sources()
    coefficient_rows, coefficient_checks = build_coefficient_audit()
    exact_models = build_exact_model_audit()
    m3 = exact_cubic_threshold()
    certificate = symbolic_certificate()
    rows = build_rows(certificate, m3)
    artifact = {
        "kind": KIND,
        "date": "2026-08-03",
        "status": "exact universal secondary cubic Jensen layer; Xi-specific and degree-uniform handoffs open",
        "proof_boundary": rows[-1].proof_boundary,
        "source_audit": source_audit(),
        "symbolic_certificate": certificate,
        "multiplicity_three": m3,
        "coefficient_audit": coefficient_rows,
        "exact_finite_model_audit": exact_models,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "sources": len(SOURCE_PATHS),
            "multiplicities": MAX_AUDIT_MULTIPLICITY - MIN_MULTIPLICITY + 1,
            "universal_coefficient_checks": coefficient_checks,
            "exact_finite_models": len(exact_models),
            "fixed_shifts": 2,
            "open_handoffs": sum(row.readiness == "open" for row in rows),
            "xi_specific_signs": 0,
            "degree_uniform_bounds": 0,
            "rh_conclusions": 0,
        },
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    counts = artifact["counts"]
    print(
        "built secondary cubic Jensen layer gate: "
        f"{counts['rows']} rows, {counts['sources']} sources, "
        f"{counts['multiplicities']} multiplicities, "
        f"{counts['universal_coefficient_checks']} universal coefficients, "
        f"{counts['exact_finite_models']} exact finite models, "
        f"{counts['fixed_shifts']} fixed shifts, "
        f"{counts['open_handoffs']} open handoffs, "
        f"{counts['xi_specific_signs']} Xi-specific signs, "
        f"{counts['degree_uniform_bounds']} degree-uniform bounds"
    )


if __name__ == "__main__":
    main()
