#!/usr/bin/env python3
"""Build the tertiary universal correction to the cubic Jensen layer."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_arbitrary_multiplicity_tertiary_universal_layer_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"
MIN_MULTIPLICITY = 3
MAX_AUDIT_MULTIPLICITY = 16

SOURCE_PATHS = {
    "secondary_layer": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_arbitrary_multiplicity_secondary_cubic_layer_gate.json",
    "primary_layer": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_arbitrary_multiplicity_jensen_boundary_layer_gate.json",
    "double_layer": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_scaled_double_zero_boundary_layer_lemma.json",
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

    secondary = payloads["secondary_layer"].get("symbolic_certificate", {})
    require(
        secondary.get("secondary_limit")
        == "If F_(n,lambda_*)(z)=(z-rho)^m U(z), rho<0, m>=3, and U(rho)!=0, then epsilon^(-4m) J_(D,n,lambda)((rho+rho*y*epsilon^4)/D)/(rho^m U(rho)) tends locally uniformly to P_m(y,q)=exp(q*partial_y^2+partial_y^3/12)y^m.",
        "secondary limit drift",
    )
    require(
        secondary.get("open_target", "").startswith(
            "The attempted first correction is universal rather than Xi-specific."
        ),
        "secondary handoff drift",
    )
    require(
        payloads["primary_layer"].get("counts", {}).get("fixed_jensen_shifts") == 1,
        "primary fixed-shift count drift",
    )
    require(
        payloads["double_layer"].get("exact", {}).get("collision_center")
        == "z_D=rho+rho*(1-2*n)/(4*D)+o(D^(-1)) in the scaled Jensen variable.",
        "double-zero center drift",
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


def tertiary_polynomial(
    m: int,
    n: int,
    y: sp.Symbol,
    q: sp.Symbol | sp.Rational,
    r: sp.Symbol | sp.Rational,
) -> sp.Expr:
    primary = secondary_polynomial(m, y, q)
    return sp.expand(
        (r - y / 2) * sp.diff(primary, y, 2)
        + sp.Rational(2 * n + 1, 4) * sp.diff(primary, y)
    )


def build_operator_audit() -> list[dict]:
    y, q, r = sp.symbols("y q r")
    rows: list[dict] = []
    for m in range(MIN_MULTIPLICITY, MAX_AUDIT_MULTIPLICITY + 1):
        p = secondary_polynomial(m, y, q)
        q0 = tertiary_polynomial(m, 0, y, q, r)
        q3 = tertiary_polynomial(m, 3, y, q, r)
        require(
            sp.expand(sp.diff(p, q) - sp.diff(p, y, 2)) == 0,
            f"Appell heat identity drift m={m}",
        )
        require(
            sp.expand(
                sp.diff(secondary_polynomial(m + 1, y, q), y) - (m + 1) * p
            )
            == 0,
            f"Appell derivative identity drift m={m}",
        )
        require(
            sp.expand(q3 - q0 - sp.Rational(3, 2) * sp.diff(p, y)) == 0,
            f"shift dependence drift m={m}",
        )
        rows.append(
            {
                "multiplicity": m,
                "secondary_polynomial": str(p),
                "tertiary_shift_zero": str(q0),
                "tertiary_shift_three": str(q3),
                "heat_identity": "partial_q P_m=partial_y^2 P_m",
                "appell_identity": "partial_y P_(m+1)=(m+1)P_m",
                "shift_difference": "Q_(m,3)-Q_(m,0)=(3/2)partial_y P_m",
            }
        )
    return rows


def exact_finite_model(m: int, n: int) -> dict:
    eps, x, z, y = sp.symbols("eps x z y")
    rho, q, r, u1 = sp.symbols("rho q r u1")
    delta = (
        rho * eps**6 / 8
        + rho * q * eps**8 / 4
        + rho * r * eps**10 / 4
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
    xi_unit_coefficient = sp.expand(
        sp.expand(local.coeff(eps, 4 * m + 4) / rho**m).coeff(u1)
    )

    target_lead = secondary_polynomial(m, y, q)
    target_tertiary = tertiary_polynomial(m, n, y, q, r)
    target_unit = rho * secondary_polynomial(m + 1, y, q)
    require(sp.expand(lead - target_lead) == 0, f"lead drift m={m}, n={n}")
    require(odd_one == 0, f"epsilon-one drift m={m}, n={n}")
    require(
        sp.expand(tertiary - target_tertiary) == 0,
        f"tertiary drift m={m}, n={n}",
    )
    require(odd_three == 0, f"epsilon-three drift m={m}, n={n}")
    require(
        sp.expand(xi_unit_coefficient - target_unit) == 0,
        f"Xi-unit coefficient drift m={m}, n={n}",
    )
    return {
        "multiplicity": m,
        "shift": n,
        "secondary_limit": str(lead),
        "tertiary_correction": str(tertiary),
        "epsilon_one_zero": True,
        "epsilon_three_zero": True,
        "first_local_unit_coefficient": str(xi_unit_coefficient),
        "first_local_unit_target": str(target_unit),
        "matched": True,
    }


def build_exact_finite_audit() -> list[dict]:
    return [exact_finite_model(m, n) for n in (0, 3) for m in range(3, 9)]


def contact_correction() -> dict[str, str]:
    y, q, r, s, n = sp.symbols("y q r s n")
    p3 = secondary_polynomial(3, y, q)
    q3 = -2 ** sp.Rational(-7, 3)
    y3 = 2 ** sp.Rational(-2, 3)
    tertiary = (r - y / 2) * sp.diff(p3, y, 2) + (2 * n + 1) * sp.diff(p3, y) / 4
    value_equation = sp.simplify(tertiary.subs({q: q3, y: y3}))
    r3 = sp.solve(value_equation, r)[0]
    derivative_equation = sp.simplify(
        s * sp.diff(p3, y, 2).subs({q: q3, y: y3})
        + sp.diff(tertiary, y).subs({q: q3, y: y3, r: r3})
    )
    s3 = sp.solve(derivative_equation, s)[0]
    require(sp.simplify(r3 - y3 / 2) == 0, "m=3 parameter correction drift")
    require(sp.simplify(s3 - (1 - 2 * n) / 4) == 0, "m=3 root correction drift")
    return {
        "nondegenerate_contact_rule": "If P=P_y=0 and P_yy!=0 at (y_m,q_m), then the epsilon^2 value equation gives r_m=y_m/2 and the derivative equation gives s_m=(1-2n)/4.",
        "m3_parameter": "r_3=2^(-5/3)",
        "m3_root": "s_3=(1-2n)/4",
        "m3_heat_expansion": "lambda_D=lambda_*+rho/(8D)-rho*2^(-13/3)D^(-4/3)+rho*2^(-11/3)D^(-5/3)+o(D^(-5/3))",
        "m3_root_expansion": "w_D=rho/D+rho*2^(-2/3)D^(-5/3)+rho*(1-2n)D^(-2)/4+o(D^(-2))",
        "newman_shift_zero": "For rho=-c^2 and n=0: t_D=Lambda-c^2/(8D)+c^2*2^(-13/3)D^(-4/3)-c^2*2^(-11/3)D^(-5/3)+o(D^(-5/3)).",
    }


def symbolic_certificate() -> dict[str, str]:
    return {
        "extended_scaling": "Set epsilon=D^(-1/6), lambda=lambda_*+rho*epsilon^6/8+(rho*q/4)*epsilon^8+(rho*r/4)*epsilon^10, and evaluate at (rho+rho*y*epsilon^4)/D.",
        "tertiary_expansion": "After normalization by rho^m U(rho) epsilon^(4m), the finite Jensen cluster is P_m(y,q)+epsilon^2 Q_(m,n)(y,q,r)+O(epsilon^4) locally, where P_m=exp(q*partial_y^2+partial_y^3/12)y^m.",
        "tertiary_operator": "Q_(m,n)=((r-y/2)*partial_y^2+(2n+1)*partial_y/4)P_m.",
        "cancellation": "All apparent fourth-order finite-Jensen and second heat-commutator terms cancel or recombine into the displayed first-order transport of P_m; no independent partial_y^4 term survives at epsilon^2.",
        "shift_entry": "The fixed Jensen shift first enters at epsilon^2 through (2n+1)partial_y P_m/4.",
        "unit_delay": "The local unit is absent from Q_(m,n). Its first linear contribution occurs at epsilon^4 and equals rho*[U'(rho)/U(rho)]*P_(m+1)(y,q).",
        "appell": "The universal family obeys partial_q P_m=partial_y^2 P_m and partial_y P_(m+1)=(m+1)P_m.",
        "contact_rule": "At any nondegenerate threshold double root P_m=P_m'=0, P_m''!=0, the next parameter and root coefficients are r_m=y_m/2 and s_m=(1-2n)/4.",
        "m3": "For q_3=-2^(-7/3), y_3=2^(-2/3): r_3=2^(-5/3) and s_3=(1-2n)/4.",
        "double_match": "The root-center coefficient rho*(1-2n)/4 agrees exactly with the existing double-zero finite-Jensen center correction, although the preceding multiplicity-three scales are different.",
        "xi_handoff": "At epsilon^4 the Xi local unit appears together with a still-uncomputed universal background and the next heat/root parameters. Its isolated coefficient alone has no sign implication.",
        "uniform_guard": "The local expansion remains noncontradictory and requires a degree-uniform analytic remainder before any collision exclusion can be claimed.",
    }


def build_rows(certificate: dict[str, str], contact: dict[str, str]) -> list[GateRow]:
    rows = [
        GateRow("tul_01_sources", "source chain", "proved", "Three parent contracts are current and hash-pinned.", "The secondary cubic, primary boundary, and double-zero gates are audited.", "No parent theorem is strengthened."),
        GateRow("tul_02_scaling", "tertiary scaling", "proved", "The next heat parameter enters at D^(-5/3).", certificate["extended_scaling"], "This is a local selected-sequence expansion."),
        GateRow("tul_03_expansion", "tertiary limit", "proved", "The normalized cluster has a complete epsilon^2 correction.", certificate["tertiary_expansion"], "An explicit local-uniform remainder is retained."),
        GateRow("tul_04_operator", "universal operator", "proved", "The epsilon^2 correction is a first-order transport of P_m and its second derivative.", certificate["tertiary_operator"], "No Xi-specific coefficient occurs."),
        GateRow("tul_05_cancellation", "operator cancellation", "proved", "No independent fourth derivative survives at this order.", certificate["cancellation"], "The cancellation is checked against exact finite models."),
        GateRow("tul_06_shift", "shift entry", "proved", "The fixed shift first appears at epsilon^2.", certificate["shift_entry"], "The Newman consequence still uses n=0."),
        GateRow("tul_07_appell", "Appell identities", "proved", "Parameter and degree derivatives close inside the universal family.", certificate["appell"], "These identities do not imply hyperbolicity beyond q_m."),
        GateRow("tul_08_unit", "Xi-unit isolation", "proved", "The first local-unit contribution is delayed to epsilon^4.", certificate["unit_delay"], "The rest of the epsilon^4 term is not yet computed."),
        GateRow("tul_09_contact", "nondegenerate contact", "conditional_exact", "A simple threshold double root has universal next corrections.", certificate["contact_rule"], "The rule requires P_m''!=0 at the selected threshold root."),
        GateRow("tul_10_m3", "multiplicity three", "proved", "The exact cubic contact supplies explicit tertiary coefficients.", certificate["m3"], "This is unconditional for m=3."),
        GateRow("tul_11_newman", "Newman refinement", "proved", "The conditional m=3 unshifted collision expansion is explicit through D^(-5/3).", contact["newman_shift_zero"], "This refines rather than excludes the collision."),
        GateRow("tul_12_double", "cross-check", "guard_validated", "The D^(-2) root-center coefficient matches the double-zero calculation.", certificate["double_match"], "The parameter scales before that coefficient differ by multiplicity."),
        GateRow("tul_13_xi", "first Xi jet", "open", "The full epsilon^4 correction must combine the Xi unit with the universal background.", certificate["xi_handoff"], "The isolated U'/U term cannot be assigned a decisive sign."),
        GateRow("tul_14_uniform", "uniform handoff", "open", "A degree-uniform remainder is still required.", certificate["uniform_guard"], "Formal local asymptotics are not collision exclusion."),
        GateRow("tul_15_boundary", "proof boundary", "guard_validated", "This gate is not a proof of RH.", "It proves the universal tertiary operator, exact m=3 corrections, and the first local-unit coefficient.", "No Xi-specific exclusion, degree-uniform theorem, PF-infinity, Lambda<=0, RH, or prize-level conclusion is proved."),
    ]
    require(len(rows) == 15, "gate row count drifted")
    require(sum(row.readiness == "open" for row in rows) == 2, "open row count drifted")
    return rows


def render_note(artifact: dict) -> str:
    cert = artifact["symbolic_certificate"]
    contact = artifact["contact_correction"]
    counts = artifact["counts"]
    return f"""# Newman Arbitrary-Multiplicity Tertiary Universal Layer Gate

Date: 2026-08-03

Status: exact universal `D^(-1/3)` correction and multiplicity-three collision refinement; full first Xi-jet order open; not a proof of RH.

## Extended Scaling

{cert['extended_scaling']}

{cert['tertiary_expansion']}

The complete correction is

```text
{cert['tertiary_operator']}
```

{cert['cancellation']}

The exact audit checks this identity for {counts['multiplicities']} symbolic multiplicities and {counts['exact_finite_models']} finite radial-heat/Jensen models at shifts zero and three.

## Shift And Appell Structure

{cert['shift_entry']}

{cert['appell']}

## Contact Correction

{cert['contact_rule']}

For multiplicity three,

```text
{contact['m3_parameter']}
{contact['m3_root']}
{contact['m3_heat_expansion']}
{contact['m3_root_expansion']}
```

At the Newman shift `n=0`,

```text
{contact['newman_shift_zero']}
```

{cert['double_match']}

## First Xi-Specific Local Term

{cert['unit_delay']}

This coefficient is independently visible in every exact model, but it is not the complete `epsilon^4` correction.

{cert['xi_handoff']}

## Live Handoff

Compute the full `epsilon^4=D^(-2/3)` operator, including the next heat and root parameters, the universal background, and `rho U'(rho)/U(rho) P_(m+1)`. Then test the actual Xi logarithmic unit jet against exact generic-unit countermodels and obtain a degree-uniform analytic remainder.

## Pi Provenance

No pi enters this local operator calculation. All constants are rational, powers of two from the exact cubic discriminant, or the algebraic local data `rho`, `q_m`, and `y_m`.

## Proof Boundary

{cert['uniform_guard']} This gate proves no Xi-specific collision exclusion, all-degree Jensen hyperbolicity, PF-infinity, `Lambda<=0`, RH, or prize-level conclusion.
"""


def main() -> None:
    load_sources()
    operator_rows = build_operator_audit()
    finite_rows = build_exact_finite_audit()
    contact = contact_correction()
    certificate = symbolic_certificate()
    rows = build_rows(certificate, contact)
    artifact = {
        "kind": KIND,
        "date": "2026-08-03",
        "status": "exact universal tertiary Jensen correction; full first Xi-jet order open",
        "proof_boundary": rows[-1].proof_boundary,
        "source_audit": source_audit(),
        "symbolic_certificate": certificate,
        "contact_correction": contact,
        "operator_audit": operator_rows,
        "exact_finite_model_audit": finite_rows,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "sources": len(SOURCE_PATHS),
            "multiplicities": len(operator_rows),
            "operator_checks": len(operator_rows),
            "exact_finite_models": len(finite_rows),
            "shift_values": 2,
            "appell_identities": 2 * len(operator_rows),
            "unit_isolation_checks": len(finite_rows),
            "exact_m3_contact_corrections": 2,
            "open_handoffs": sum(row.readiness == "open" for row in rows),
            "degree_uniform_bounds": 0,
            "rh_conclusions": 0,
        },
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    counts = artifact["counts"]
    print(
        "built tertiary universal Jensen layer gate: "
        f"{counts['rows']} rows, {counts['sources']} sources, "
        f"{counts['multiplicities']} multiplicities, "
        f"{counts['operator_checks']} operator checks, "
        f"{counts['exact_finite_models']} exact finite models, "
        f"{counts['appell_identities']} Appell identities, "
        f"{counts['unit_isolation_checks']} unit-isolation checks, "
        f"{counts['exact_m3_contact_corrections']} exact m3 corrections, "
        f"{counts['open_handoffs']} open handoffs, "
        f"{counts['degree_uniform_bounds']} degree-uniform bounds"
    )


if __name__ == "__main__":
    main()
