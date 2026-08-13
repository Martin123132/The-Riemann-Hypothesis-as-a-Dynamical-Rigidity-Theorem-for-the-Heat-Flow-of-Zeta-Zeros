#!/usr/bin/env python3
"""Build the arbitrary-multiplicity cofinal Jensen boundary-layer gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_arbitrary_multiplicity_jensen_boundary_layer_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"
MIN_MULTIPLICITY = 2
MAX_AUDIT_MULTIPLICITY = 16

SOURCE_PATHS = {
    "positive_boundary": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_positive_boundary_attainment_lemma.json",
    "double_layer": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_scaled_double_zero_boundary_layer_lemma.json",
    "cofinal_scaling": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_cofinal_scaling_limit_equivalence_gate.json",
    "polar_cascade": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_polar_heat_collision_cascade_lemma.json",
    "degree361": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_real_zero_band_degree361_sector_certificate.json",
    "fixed_collision": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_first_boundary_eventual_hankel_escape_gate.json",
    "radial_heat": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_edrei_heat_flow_boundary_gate.json",
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


def row_by_id(payload: dict, row_id: str) -> dict:
    matches = [row for row in payload.get("rows", []) if row.get("id") == row_id]
    require(len(matches) == 1, f"row drift: {row_id}")
    return matches[0]


def load_sources() -> dict[str, dict]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        require(path.exists(), f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))

    require(
        payloads["positive_boundary"].get("exact", {}).get("attainment_theorem")
        == "If Lambda>0, H_Lambda has a finite real multiple zero c with |c|<=X_Lambda=exp(2*C/Lambda). Thus a positive Newman boundary cannot be realized solely by zero collisions escaping to infinity.",
        "positive-boundary attainment drift",
    )
    require(
        row_by_id(payloads["positive_boundary"], "npba_07_hermite_cluster_split").get("readiness")
        == "available_published",
        "arbitrary-multiplicity source drift",
    )
    require(
        payloads["double_layer"].get("exact", {}).get("boundary_layer")
        == "If lambda=lambda_*+tau/D and z=rho+eta/sqrt(D), then 2*D*J_D(z/D)/F''_*(rho)->eta^2+8*rho*tau-rho^2.",
        "double-zero layer drift",
    )
    require(
        payloads["double_layer"].get("exact", {}).get("heat_pde")
        == "partial_lambda F_n=(4*n+2)*partial_z F_n+4*z*partial_z^2 F_n",
        "coefficient heat PDE drift",
    )
    require(
        payloads["cofinal_scaling"].get("exact", {}).get("equivalence")
        == "cofinal degree hyperbolicity at fixed n <=> F_n is Laguerre-Polya",
        "cofinal fixed-shift equivalence drift",
    )
    require(
        payloads["polar_cascade"].get("exact", {}).get("degree_escape")
        == "At a non-exponential-polynomial LP boundary, the least nonhyperbolic Jensen degree on the bad side tends to infinity as lambda approaches the boundary.",
        "polar degree escape drift",
    )
    require(
        row_by_id(payloads["degree361"], "nrzb361_10_degree361").get("readiness")
        == "ready_to_apply",
        "degree-361 theorem drift",
    )
    require(
        payloads["fixed_collision"].get("counts", {}).get("uniform_in_order_thresholds") == 0,
        "fixed-collision quantifier guard drift",
    )
    require(
        payloads["radial_heat"].get("exact", {}).get("radial_heat")
        == "partial_lambda F_lambda=(4*z*partial_z^2+2*partial_z)F_lambda",
        "radial heat identity drift",
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


def universal_polynomial(m: int, eta: sp.Symbol, a: sp.Symbol | sp.Rational) -> sp.Expr:
    return sp.expand(
        sum(
            sp.factorial(m)
            * a**q
            * eta ** (m - 2 * q)
            / (sp.factorial(q) * sp.factorial(m - 2 * q))
            for q in range(m // 2 + 1)
        )
    )


def build_multiplicity_audit() -> tuple[list[dict], dict[str, int]]:
    eta = sp.symbols("eta")
    a = sp.symbols("a")
    rows: list[dict] = []
    coefficient_checks = 0
    hermite_checks = 0
    imaginary_checks = 0
    orientation_checks = 0

    for m in range(MIN_MULTIPLICITY, MAX_AUDIT_MULTIPLICITY + 1):
        polynomial = universal_polynomial(m, eta, a)
        coefficients = []
        for q in range(m // 2 + 1):
            coefficient = sp.factorial(m) // (
                sp.factorial(q) * sp.factorial(m - 2 * q)
            )
            require(
                sp.expand(polynomial).coeff(eta, m - 2 * q).coeff(a, q)
                == coefficient,
                f"universal coefficient drift at m={m}, q={q}",
            )
            coefficients.append(
                {
                    "q": q,
                    "eta_power": m - 2 * q,
                    "integer_coefficient": int(coefficient),
                }
            )
            coefficient_checks += 1

        hermite = sp.hermite_prob(m, eta)
        require(
            sp.expand(polynomial.subs(a, -sp.Rational(1, 2)) - hermite) == 0,
            f"Hermite good-side identity drift at m={m}",
        )
        hermite_checks += 1

        imaginary_form = sp.expand(
            sp.I**m * sp.hermite_prob(m, eta / sp.I)
        )
        require(
            sp.expand(polynomial.subs(a, sp.Rational(1, 2)) - imaginary_form) == 0,
            f"imaginary-side identity drift at m={m}",
        )
        imaginary_checks += 1

        rho, epsilon = sp.symbols("rho epsilon", negative=True)
        threshold = rho / 8
        a_good = sp.expand(4 * rho * (threshold + epsilon) - rho**2 / 2)
        a_bad = sp.expand(4 * rho * (threshold - epsilon) - rho**2 / 2)
        require(a_good == 4 * epsilon * rho, f"good-side orientation drift at m={m}")
        require(a_bad == -4 * epsilon * rho, f"bad-side orientation drift at m={m}")
        orientation_checks += 1

        rows.append(
            {
                "multiplicity": m,
                "term_count": len(coefficients),
                "coefficients": coefficients,
                "universal_polynomial": str(polynomial),
                "good_side_identity": f"He_{m}(eta) at a=-1/2",
                "good_side_real_simple_roots": m,
                "bad_side_nonreal_roots": m - (m % 2),
                "bad_side_real_roots": m % 2,
                "collision_parameter_scale": "tau=rho/8",
                "finite_collision_scale": (
                    "lambda_j=lambda_*+rho/(8D_j)+o(D_j^-1) along a collision sequence"
                ),
            }
        )

    return rows, {
        "universal_coefficient_checks": coefficient_checks,
        "hermite_good_side_checks": hermite_checks,
        "imaginary_bad_side_checks": imaginary_checks,
        "orientation_checks": orientation_checks,
    }


def build_exact_polynomial_layer_audit() -> list[dict]:
    h, x, z, eta, rho, tau = sp.symbols("h x z eta rho tau")
    rows: list[dict] = []
    for n in (0, 3):
        for m in range(2, 9):
            source = (x - rho) ** m

            def heat_generator(polynomial: sp.Expr) -> sp.Expr:
                return sp.expand(
                    4 * x * sp.diff(polynomial, x, 2)
                    + (4 * n + 2) * sp.diff(polynomial, x)
                )

            heat_polynomial = sp.Integer(0)
            iterate = source
            for order in range(m + 1):
                heat_polynomial += (tau * h**2) ** order * iterate / sp.factorial(order)
                iterate = heat_generator(iterate)

            scaled_jensen = sp.Integer(0)
            for (degree,), coefficient in sp.Poly(sp.expand(heat_polynomial), x).terms():
                falling_ratio = sp.prod(1 - index * h**2 for index in range(degree))
                scaled_jensen += coefficient * falling_ratio * z**degree

            local_expression = sp.expand(
                scaled_jensen.subs(z, rho + eta * h) / h**m
            )
            exact_limit = sp.expand(sp.limit(local_expression, h, 0))
            a = 4 * rho * tau - rho**2 / 2
            target = universal_polynomial(m, eta, a)
            require(
                sp.expand(exact_limit - target) == 0,
                f"exact polynomial layer drift at n={n}, m={m}",
            )
            rows.append(
                {
                    "shift": n,
                    "multiplicity": m,
                    "limit": str(exact_limit),
                    "target": str(target),
                    "matched": True,
                }
            )
    return rows


def symbolic_certificate() -> dict[str, str]:
    return {
        "newman_map": "F_(0,t)(s)=2*H_t(i*sqrt(s)); a finite real multiplicity-m zero c of H_t maps to the finite negative multiplicity-m zero rho=-c^2 of F_(0,t).",
        "jensen_operator": "For T_DF(z)=J_(D,n,lambda)(z/D), T_DF(z)=[(1+(z/D)*partial_x)^D F(x)]_(x=0). After translation to x=z, its logarithm is -z^2*partial_x^2/(2D)+sum_(r>=3)(-1)^(r+1)z^r*partial_x^r/(rD^(r-1)).",
        "jensen_local_limit": "Under z=rho+eta/sqrt(D), the quadratic operator tends to exp(-(rho^2/2)*partial_eta^2); every r>=3 logarithmic term is O(D^(1-r/2)) and vanishes on compact eta sets.",
        "heat_local_limit": "For lambda=lambda_*+tau/D, the coefficient PDE tends under the same scaling to exp(4*rho*tau*partial_eta^2); the drift term is O(D^-1/2).",
        "combined_layer": "If F_*(z)=(z-rho)^m U(z), U(rho)!=0, then D^(m/2)J_(D,n,lambda_*+tau/D)((rho+eta/sqrt(D))/D)/U(rho) tends locally uniformly to K_(m,a)(eta)=exp(a*partial_eta^2)eta^m, a=4*rho*tau-rho^2/2.",
        "coefficient_formula": "K_(m,a)(eta)=m!*sum_(q=0)^floor(m/2) a^q*eta^(m-2q)/(q!*(m-2q)!).",
        "good_side": "For a=-sigma^2/2<0, K_(m,a)(eta)=sigma^m*He_m(eta/sigma), so all m roots are real and simple.",
        "bad_side": "For a=beta^2/2>0, K_(m,a)(eta)=(i*beta)^m*He_m(eta/(i*beta)); for m>=2 it has a nonreal conjugate pair (and only the root zero is real when m is odd).",
        "collision_scale": "The layer changes type only at a=0, hence tau=rho/8. Root continuity and a diagonal compactness choice give a sequence D_j->infinity of finite Jensen collisions with lambda_j=lambda_*+rho/(8D_j)+o(D_j^-1) and polynomial root w_j=rho/D_j+o(D_j^(-3/2)).",
        "positive_boundary": "If Lambda>0, positive-boundary attainment supplies finite c and m>=2. At fixed shift n=0, rho=-c^2, so a sequence D_j->infinity has collisions t_j=Lambda-c^2/(8D_j)+o(D_j^-1). No unbounded-shift sequence is needed.",
        "fixed_order_guard": "The certified Section 11.195 tail is on [-100,0], not the positive Newman interval. Even if extended separately at every fixed order, a finite double collision of degree D uses determinant order k=D; more generally bounded collision multiplicity gives k=D-m_D+2 tending to infinity. Fixed-order eventual thresholds therefore do not control this layer without a new order-uniform theorem.",
        "open_target": "Prove an Xi-specific estimate uniform in degree throughout the cofinal rho/(8D) layer, or use the independent first-jet boundary-winding route. Bounded-degree certificates and fixed-order Hankel tails cannot be promoted to this target.",
    }


def build_gate_rows(certificate: dict[str, str], counts: dict[str, int]) -> list[GateRow]:
    rows = [
        GateRow("amb_01_sources", "source chain", "proved", "Seven parent artifacts are current and hash-pinned.", "Positive-boundary attainment, the double layer, cofinal scaling, polar escape, degree 361, fixed-collision, and radial-heat inputs are audited.", "No source theorem is strengthened silently."),
        GateRow("amb_02_fixed_shift", "Newman coordinate map", "proved", "Every hypothetical positive Newman boundary supplies a finite negative multiple zero at fixed Jensen shift zero.", certificate["newman_map"], "This is conditional on Lambda>0."),
        GateRow("amb_03_jensen_operator", "Jensen operator", "proved", "The scaled Jensen transform has an exact binomial differential-operator representation.", certificate["jensen_operator"], "The logarithmic series is used asymptotically on compact analytic germs."),
        GateRow("amb_04_jensen_limit", "Jensen local limit", "proved", "Only the quadratic Jensen operator survives in the parabolic root scale.", certificate["jensen_local_limit"], "A compact-set analytic remainder is required; no degree-uniform global norm is claimed."),
        GateRow("amb_05_heat_limit", "heat local limit", "proved", "The radial coefficient heat flow has a constant-coefficient parabolic limit at rho.", certificate["heat_local_limit"], "The fixed shift n affects only lower-order drift."),
        GateRow("amb_06_combined", "universal layer", "proved", "Every finite multiplicity has one universal Jensen/heat boundary layer.", certificate["combined_layer"], "The zero rho must be finite, negative, and of finite multiplicity."),
        GateRow("amb_07_coefficients", "universal polynomial", "proved", "The complete arbitrary-multiplicity polynomial is explicit.", certificate["coefficient_formula"], f"The finite audit checks multiplicities {MIN_MULTIPLICITY} through {MAX_AUDIT_MULTIPLICITY}; the formula is symbolic in m."),
        GateRow("amb_08_good_side", "Hermite splitting", "proved", "The good side has m simple real local roots.", certificate["good_side"], "This is a local cluster statement, not global polynomial hyperbolicity."),
        GateRow("amb_09_bad_side", "nonreal splitting", "proved", "The bad side has a nonreal local pair for every m>=2.", certificate["bad_side"], "For odd m one central local root remains real."),
        GateRow("amb_10_collision", "collision scale", "proved", "Finite-degree collisions converge to the universal rho/(8D) parameter layer.", certificate["collision_scale"], "Only the double-zero parent supplies the explicit D^-2 correction."),
        GateRow("amb_11_boundary", "positive-boundary consequence", "proved", "A positive Newman boundary forces a fixed-shift unbounded-degree Jensen collision sequence for every finite boundary multiplicity.", certificate["positive_boundary"], "This reduction is not a contradiction."),
        GateRow("amb_12_degree361", "bounded-degree guard", "guard_validated", "The existing all-shift theorem through degree 361 cannot meet a degree-divergent sequence.", "The degree-361 certificate remains exact and useful, but every cofinal collision sequence eventually lies above its cutoff.", "No finite cutoff is promoted to a cofinal theorem."),
        GateRow("amb_13_hankel", "fixed-order guard", "guard_validated", "The eventual-Hankel escape theorem and the cofinal boundary layer have different quantifiers.", certificate["fixed_order_guard"], "No threshold uniform in determinant order is available."),
        GateRow("amb_14_route", "uniform cofinal handoff", "open", "The live Jensen-side target is an Xi-specific estimate uniform in degree on the rho/(8D) layer.", certificate["open_target"], "The direct first-jet boundary-degree route remains an independent alternative."),
        GateRow("amb_15_proof_boundary", "proof boundary", "guard_validated", "This gate is not a proof of RH.", "It removes the double-zero assumption from the cofinal Jensen reduction and fixes the relevant shift at n=0.", "No cofinal exclusion, all-degree hyperbolicity, PF-infinity, Lambda<=0, RH, or prize-level conclusion is proved."),
    ]
    require(len(rows) == 15, "gate row count drifted")
    require(sum(row.readiness == "open" for row in rows) == 1, "open row drifted")
    return rows


def render_note(artifact: dict) -> str:
    cert = artifact["symbolic_certificate"]
    counts = artifact["counts"]
    return f"""# Newman Arbitrary-Multiplicity Jensen Boundary-Layer Gate

Date: 2026-08-03

Status: exact arbitrary-multiplicity fixed-shift cofinal Jensen boundary layer; one degree-uniform Xi handoff open; not a proof of RH.

## Positive-Boundary Coordinate

{cert['newman_map']}

Positive-boundary attainment therefore removes spatial escape for the contradiction `Lambda>0`. The relevant Jensen sequence can be taken at the single shift `n=0`; only its degree must diverge.

## Exact Jensen Operator

{cert['jensen_operator']}

{cert['jensen_local_limit']}

## Heat Scaling

{cert['heat_local_limit']}

## Universal Multiplicity Layer

{cert['combined_layer']}

Equivalently,

```text
{cert['coefficient_formula']}
```

The machine audit checks every coefficient, Hermite specialization, imaginary specialization, and parameter orientation for multiplicities `{MIN_MULTIPLICITY}` through `{MAX_AUDIT_MULTIPLICITY}`: {counts['universal_coefficient_checks']} coefficient checks and {counts['hermite_good_side_checks']} checks on each side. It also computes {counts['exact_polynomial_layer_checks']} exact finite-Jensen heat models at shifts zero and three and recovers the same limit independently.

## Root Geometry

{cert['good_side']}

{cert['bad_side']}

{cert['collision_scale']}

For `m=2`, `K_(2,a)=eta^2+2a=eta^2+8rho*tau-rho^2`, exactly recovering the existing double-zero gate. The present theorem removes the unproved assumption that a hypothetical Newman boundary zero is double.

## Newman Consequence

{cert['positive_boundary']}

Thus the prize-relevant Jensen escape is fixed shift and unbounded degree. A separate large-shift compactness theorem is not needed for this contradiction route.

## Quantifier Guard

{cert['fixed_order_guard']}

The all-shift degree-361 sector theorem is bounded-degree evidence. The fixed-order eventual-Hankel theorem excludes attained finite collisions when its upper-shift hypotheses hold, but supplies no determinant threshold uniform in the growing order used here.

## Live Handoff

{cert['open_target']}

## Pi Provenance

The universal layer introduces no pi. The value `rho=-c^2` comes from the algebraic map `s=-z^2`; pi appears only in separate Xi kernel normalizations and estimates.

## Proof Boundary

This gate proves the arbitrary-multiplicity local limit and the conditional fixed-shift cofinal collision sequence under `Lambda>0`. It proves no degree-uniform exclusion, all-degree Jensen hyperbolicity, PF-infinity, `Lambda<=0`, RH, or prize-level conclusion.
"""


def main() -> int:
    load_sources()
    multiplicity_rows, audit_counts = build_multiplicity_audit()
    exact_polynomial_rows = build_exact_polynomial_layer_audit()
    certificate = symbolic_certificate()
    counts = {
        "gate_rows": 15,
        "source_artifacts": len(SOURCE_PATHS),
        "multiplicity_rows": len(multiplicity_rows),
        "exact_polynomial_layer_checks": len(exact_polynomial_rows),
        **audit_counts,
        "fixed_jensen_shifts": 1,
        "uniform_in_degree_bounds": 0,
        "rh_conclusions": 0,
    }
    gate_rows = build_gate_rows(certificate, counts)
    artifact = {
        "kind": KIND,
        "date": "2026-08-03",
        "status": "exact arbitrary-multiplicity fixed-shift cofinal Jensen boundary layer with one degree-uniform Xi handoff open",
        "counts": counts,
        "source_audit": source_audit(),
        "symbolic_certificate": certificate,
        "multiplicity_audit": multiplicity_rows,
        "exact_polynomial_layer_audit": exact_polynomial_rows,
        "rows": [asdict(row) for row in gate_rows],
        "proof_boundary": "This proves the arbitrary-multiplicity parabolic Jensen layer and, conditionally on Lambda>0, a fixed-shift degree-divergent collision sequence with t_j=Lambda-c^2/(8D_j)+o(D_j^-1). It provides no degree-uniform exclusion, all-degree Jensen hyperbolicity, PF-infinity, Lambda<=0, RH, or prize-level conclusion.",
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    print(
        "built arbitrary-multiplicity Jensen boundary-layer gate: "
        f"{counts['gate_rows']} rows, {counts['source_artifacts']} sources, "
        f"{counts['multiplicity_rows']} multiplicities, "
        f"{counts['universal_coefficient_checks']} universal coefficients, "
        f"{counts['hermite_good_side_checks']} Hermite checks, "
        f"{counts['imaginary_bad_side_checks']} imaginary-side checks, "
        f"{counts['orientation_checks']} orientation checks, "
        f"{counts['exact_polynomial_layer_checks']} exact polynomial layers, "
        f"{counts['fixed_jensen_shifts']} fixed shift, "
        f"{counts['uniform_in_degree_bounds']} uniform-degree bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
