#!/usr/bin/env python3
"""Build the five-current Mangoldt normal-form gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "contiguous_terminal_tail_anchored_five_current_mangoldt_normal_form_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "five_current": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_current_bulk_closure_"
        "reduction.json"
    ),
    "mangoldt_abel": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "mangoldt_abel_contact_gate.json"
    ),
    "contact_centering": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "mangoldt_contact_centering_gate.json"
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


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCE_PATHS.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCE_PATHS.items()
    }


def source_audit(payloads: dict[str, dict]) -> dict:
    five = payloads["five_current"]
    if five.get("counts", {}).get("five_current_algebraic_closures") != 1:
        raise RuntimeError("five-current closure source drifted")
    if five.get("counts", {}).get("signed_bulk_bounds") != 0:
        raise RuntimeError("five-current proof boundary drifted")

    mangoldt = payloads["mangoldt_abel"].get("exact", {})
    if "1+d_n=a_0+a_1 log(n)+a_2 log(n)^2" not in mangoldt.get(
        "correction_polynomial", ""
    ):
        raise RuntimeError("quadratic correction source drifted")
    if "H_j=sum_(dm<=N)Lambda(d)" not in mangoldt.get(
        "mangoldt_moments", ""
    ):
        raise RuntimeError("Mangoldt moment source drifted")

    centering = json.dumps(payloads["contact_centering"], sort_keys=True)
    for marker in ("balanced", "Z_Lambda", "dm<=N"):
        if marker not in centering:
            raise RuntimeError(f"contact-centering marker missing: {marker}")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCE_PATHS.items()
        },
        "source_kinds": {
            key: payload.get("kind", "") for key, payload in payloads.items()
        },
        "imported_signed_type_ii_bounds": 0,
    }


def symbolic_certificate() -> dict:
    ell = sp.symbols("ell", real=True)
    rho_1, rho_2, rho_1_x, rho_2_x = sp.symbols(
        "rho_1 rho_2 rho_1_x rho_2_x"
    )
    q_0 = sp.symbols("q_0")
    correction = 1 + rho_1 * ell + rho_2 * ell**2
    correction_x = rho_1_x * ell + rho_2_x * ell**2
    delta = correction_x / correction
    if sp.simplify(delta * correction * q_0 - correction_x * q_0) != 0:
        raise RuntimeError("correction-current cancellation failed")

    g = sp.symbols("G_0:5")
    corrected = {
        "H_0": g[0] + rho_1 * g[1] + rho_2 * g[2],
        "H_1": g[1] + rho_1 * g[2] + rho_2 * g[3],
        "H_2": g[2] + rho_1 * g[3] + rho_2 * g[4],
        "D_0": rho_1_x * g[1] + rho_2_x * g[2],
        "D_1": rho_1_x * g[2] + rho_2_x * g[3],
    }
    matrix = sp.Matrix(
        [
            [1, rho_1, rho_2, 0, 0],
            [0, 1, rho_1, rho_2, 0],
            [0, 0, 1, rho_1, rho_2],
            [0, rho_1_x, rho_2_x, 0, 0],
            [0, 0, rho_1_x, rho_2_x, 0],
        ]
    )
    vector = sp.Matrix(g)
    transformed = matrix * vector
    if any(
        sp.expand(transformed[index] - corrected[key]) != 0
        for index, key in enumerate(("H_0", "H_1", "H_2", "D_0", "D_1"))
    ):
        raise RuntimeError("five-current correction transform failed")

    transform_determinant = sp.factor(matrix.det())
    expected_determinant = sp.factor(
        rho_2
        * (
            -rho_1 * rho_1_x * rho_2_x
            + rho_2 * rho_1_x**2
            + rho_2_x**2
        )
    )
    if sp.expand(transform_determinant - expected_determinant) != 0:
        raise RuntimeError("correction transform determinant failed")

    log_p = sp.symbols("log_p", positive=True, real=True)
    negative_minors: dict[str, str] = {}
    for order in range(1, 5):
        off_diagonal = log_p**order / 2
        determinant = sp.factor(-off_diagonal**2)
        if determinant != -log_p ** (2 * order) / 4:
            raise RuntimeError(f"prime-edge minor failed at order {order}")
        negative_minors[str(order)] = f"-log(p)^{2 * order}/4<0"

    return {
        "quadratic_correction": (
            "c_n=(1+d_n)/(1+d_1)=1+rho_1*log(n)+rho_2*log(n)^2"
        ),
        "correction_current": (
            "delta_n=partial_x log(c_n), so delta_n*c_n="
            "rho_(1,x)*log(n)+rho_(2,x)*log(n)^2"
        ),
        "correction_free_moments": (
            "G_r^(B)=(eta/S_a)sum_(n=1)^B log(n)^r*"
            "exp[t log(n)^2/4-s_*log(n)], 0<=r<=4"
        ),
        "five_current_transform": (
            "H_0=G_0+rho_1G_1+rho_2G_2; "
            "H_1=G_1+rho_1G_2+rho_2G_3; "
            "H_2=G_2+rho_1G_3+rho_2G_4; "
            "D_0=rho_(1,x)G_1+rho_(2,x)G_2; "
            "D_1=rho_(1,x)G_2+rho_(2,x)G_3"
        ),
        "transform_determinant": str(transform_determinant),
        "mangoldt_moments": (
            "For 1<=r<=4, G_r^(B)=(eta/S_a)sum_(dm<=B)"
            "Lambda(d)log(dm)^(r-1)q_(dm)^(0)="
            "(eta/(2S_a))sum_(dm<=B)[Lambda(d)+Lambda(m)]"
            "log(dm)^(r-1)q_(dm)^(0)"
        ),
        "heat_factorization": (
            "q_(dm)^(0)=q_d^(0)q_m^(0)"
            "exp[(t/2)log(d)log(m)]"
        ),
        "negative_prime_edge_minors": negative_minors,
    }


def arithmetic_certificate() -> dict:
    checked_cutoffs = (17, 31, 64, 97)
    divisor_rows = 0
    symmetric_rows = 0
    balanced_rows = 0

    for cutoff in checked_cutoffs:
        for n in range(2, cutoff + 1):
            factors = sp.factorint(n)
            primes = tuple(sorted(factors))
            symbols = {prime: sp.Symbol(f"L_{prime}") for prime in primes}
            log_n = sum(factors[prime] * symbols[prime] for prime in primes)

            def mangoldt(integer: int) -> sp.Expr:
                factorization = sp.factorint(integer)
                if len(factorization) != 1:
                    return sp.Integer(0)
                prime = next(iter(factorization))
                return symbols.get(prime, sp.Symbol(f"L_{prime}"))

            divisors = sp.divisors(n)
            divisor_sum = sum(mangoldt(divisor) for divisor in divisors)
            if sp.expand(divisor_sum - log_n) != 0:
                raise RuntimeError(f"Mangoldt divisor identity failed at {n}")
            divisor_rows += 1

            symmetric_sum = sum(
                (mangoldt(divisor) + mangoldt(n // divisor)) / 2
                for divisor in divisors
            )
            if sp.expand(symmetric_sum - log_n) != 0:
                raise RuntimeError(f"symmetric Mangoldt identity failed at {n}")
            symmetric_rows += 1

        root = int(cutoff**0.5)
        all_pairs = {
            (d, m)
            for d in range(1, cutoff + 1)
            for m in range(1, cutoff // d + 1)
        }
        square = {(d, m) for d, m in all_pairs if d <= root and m <= root}
        wings = {
            (d, m)
            for d, m in all_pairs
            if (d <= root < m) or (m <= root < d)
        }
        if square | wings != all_pairs or square & wings:
            raise RuntimeError(f"balanced hyperbola failed at {cutoff}")
        balanced_rows += 1

    return {
        "cutoffs": list(checked_cutoffs),
        "divisor_rows": divisor_rows,
        "symmetric_rows": symmetric_rows,
        "balanced_hyperbola_rows": balanced_rows,
        "max_moment_order": 4,
    }


def exact_payload() -> dict:
    return {
        "domain": (
            "Use the fixed physical q=1 terminal/bulk partition with "
            "B=N-M-1. All sums below stop at B and share the same eta/S_a "
            "normalizer as the certified tail."
        ),
        "coefficient_definitions": (
            "From 1+d_n=a_0+a_1 log(n)+a_2 log(n)^2 and a_0=1+d_1, "
            "put rho_1=a_1/a_0 and rho_2=a_2/a_0. Their x derivatives "
            "are exact quotient derivatives; no correction remainder is "
            "discarded."
        ),
        "five_moment_closure": (
            "The corrected five-current vector needed by Phi_B is an exact "
            "linear image of G_0,G_1,G_2,G_3,G_4. Thus the correction does "
            "not create an uncontrolled sixth current."
        ),
        "balanced_hyperbola": (
            "For D=floor(sqrt(B)), each symmetric Mangoldt moment splits "
            "exactly into the square d,m<=D and the two wings with exactly "
            "one variable <=D. The endpoint-terminal block is composed with "
            "the resulting moments before any absolute values."
        ),
        "type_ii_target": (
            "Insert the five-current transform into the exact Phi_B quadratic "
            "form and prove Phi_B<=1/400, preferably <=1/800, by a signed "
            "endpoint-coupled Type-I/II or Vaughan estimate for G_1,...,G_4 "
            "together with the value moment G_0."
        ),
        "indefiniteness_guard": (
            "For every 1<=r<=4 and any prime sqrt(B)<p<=B, the {1,p} "
            "principal minor of the real transpose-symmetric Mangoldt kernel "
            "is [[0,log(p)^r/2],[log(p)^r/2,0]] and has negative determinant. "
            "Transpose symmetry and entrywise positivity therefore do not "
            "supply the signed Phi_B bound."
        ),
        "pi_provenance": (
            "No pi is introduced by the correction transform, Mangoldt "
            "identity, or balanced hyperbola. Existing pi factors retain "
            "their completed-zeta/Riemann-Siegel origin."
        ),
        "proof_boundary": (
            "This proves an exact five-moment correction-free normal form, "
            "Mangoldt representations through order four, balanced hyperbola "
            "organization, and an indefiniteness guard. It proves no Type-I/"
            "II cancellation estimate, signed bound on Phi_B, bulk aggregate "
            "sign closure, complete retained or Xi-level current theorem, "
            "Abel gap, winding cap, contact exclusion, Q209, Lambda<=0, "
            "PF-infinity, RH, or prize-level conclusion."
        ),
    }


def build_rows(
    exact: dict, symbolic: dict, arithmetic: dict
) -> list[GateRow]:
    return [
        GateRow("ttmg_01_domain", "fixed partition", "proved", "The normal form uses the certified tail's exact bulk cutoff.", exact["domain"], "No adjacent cutoff is crossed."),
        GateRow("ttmg_02_correction", "quadratic correction", "proved", "The retained first Dirichlet correction is exactly quadratic in log(n).", symbolic["quadratic_correction"], "No asymptotic correction truncation is made."),
        GateRow("ttmg_03_current", "correction current", "proved", "Multiplication by delta_n cancels the correction denominator exactly.", symbolic["correction_current"], "rho_(j,x) is retained."),
        GateRow("ttmg_04_base", "correction-free moments", "proved", "Five correction-free logarithmic moments are sufficient.", symbolic["correction_free_moments"], "The common anchor and scale are preserved."),
        GateRow("ttmg_05_transform", "five-current transform", "proved", "The corrected current vector is an exact linear image of G_0 through G_4.", symbolic["five_current_transform"], "Invertibility is neither assumed nor required.", {"determinant": symbolic["transform_determinant"]}),
        GateRow("ttmg_06_mangoldt", "Mangoldt normal form", "proved", "Every positive-order base moment has a symmetric divisor-pair representation.", symbolic["mangoldt_moments"], "The form is complex bilinear, not Hermitian.", arithmetic),
        GateRow("ttmg_07_heat", "heat factorization", "proved", "The multiplicative pair retains the exact heat cross coupling.", symbolic["heat_factorization"], "No Gaussian fibre sign is inferred."),
        GateRow("ttmg_08_hyperbola", "balanced hyperbola", "proved", "The Type-I/II square and wings form an exact partition.", exact["balanced_hyperbola"], "This is organization, not cancellation.", arithmetic),
        GateRow("ttmg_09_join", "endpoint-coupled join", "proved", "The Mangoldt moments feed the exact terminal-tail Phi_B form.", exact["five_moment_closure"], "The endpoint is not estimated separately."),
        GateRow("ttmg_10_guard", "indefiniteness guard", "proved", "The symmetric kernels through order four are indefinite on the retained hyperbola.", exact["indefiniteness_guard"], "Positivity promotion is forbidden.", symbolic["negative_prime_edge_minors"]),
        GateRow("ttmg_11_target", "Type-I/II target", "open", "A signed endpoint-coupled cancellation theorem is now the exact arithmetic obligation.", exact["type_ii_target"], "No such estimate is proved here."),
        GateRow("ttmg_12_boundary", "proof boundary", "open", "The normal form narrows but does not close the bulk theorem.", exact["proof_boundary"], "RH remains open."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    symbolic = payload["symbolic_certificate"]
    arithmetic = payload["arithmetic_certificate"]
    return f"""# Terminal-Tail Five-Current Mangoldt Normal Form

Date: 2026-07-31

Status: exact correction-free and Mangoldt normal form with `0 signed Type-I/II
bounds`. This is not a proof artifact; `Phi_B` and RH remain open.

## Exact Correction Transform

{exact['domain']}

```text
{symbolic['quadratic_correction']}
{symbolic['correction_current']}
{symbolic['correction_free_moments']}
```

{exact['coefficient_definitions']}

The corrected five-current vector is

```text
{symbolic['five_current_transform']}
```

{exact['five_moment_closure']}

## Mangoldt Form Through Order Four

```text
{symbolic['mangoldt_moments']}
{symbolic['heat_factorization']}
```

The formal divisor audit checked `{arithmetic['divisor_rows']}` coefficient
rows and `{arithmetic['symmetric_rows']}` symmetric rows through cutoff
`{max(arithmetic['cutoffs'])}`.

{exact['balanced_hyperbola']}

## Exact Guard

{exact['indefiniteness_guard']}

The four determinants are

```text
{json.dumps(symbolic['negative_prime_edge_minors'], indent=2)}
```

Thus the old visual symmetry has a precise role as a Type-I/II coordinate,
but it is not a positivity theorem.

## Open Arithmetic Target

```text
{exact['type_ii_target']}
```

## Boundary

{exact['pi_provenance']}

{exact['proof_boundary']}

## Reproduce

```text
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```
"""


def build_payload() -> dict:
    payloads = load_sources()
    source = source_audit(payloads)
    symbolic = symbolic_certificate()
    arithmetic = arithmetic_certificate()
    exact = exact_payload()
    rows = build_rows(exact, symbolic, arithmetic)
    return {
        "kind": STEM,
        "date": "2026-07-31",
        "status": (
            "exact terminal-tail five-current correction-free Mangoldt "
            "normal form; Type-I/II Phi_B bound open; no bulk aggregate sign "
            "closure, Xi-level current theorem, Lambda<=0, or RH"
        ),
        "proof_boundary": exact["proof_boundary"],
        "source_audit": source,
        "symbolic_certificate": symbolic,
        "arithmetic_certificate": arithmetic,
        "exact": exact,
        "counts": {
            "rows": len(rows),
            "exact_correction_transforms": 1,
            "corrected_complex_currents": 5,
            "correction_free_complex_moments": 5,
            "mangoldt_moment_orders": 4,
            "balanced_hyperbola_schemas": 1,
            "negative_prime_edge_minor_families": 4,
            "signed_type_ii_targets": 1,
            "signed_type_ii_bounds": 0,
            "phi_b_bounds": 0,
            "xi_level_current_theorems": 0,
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
        "wrote five-current Mangoldt normal form: "
        f"{payload['counts']['rows']} rows, 5 correction-free moments, "
        "4 Mangoldt orders, 0 signed Type-I/II bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
