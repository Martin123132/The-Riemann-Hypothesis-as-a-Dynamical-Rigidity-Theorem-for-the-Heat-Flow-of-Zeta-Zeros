#!/usr/bin/env python3
"""Build the terminal-tail-anchored five-current bulk closure."""

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
    "contiguous_terminal_tail_anchored_five_current_bulk_closure_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "tail_abel_join": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_abel_cross_current_reduction.json"
    ),
    "joined_dyadic": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "joined_dyadic_odd_prefix_first_jet_reduction.json"
    ),
    "normalized_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "normalized_prefix_phase_flux_reduction.json"
    ),
}


@dataclass(frozen=True)
class ClosureRow:
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
    join = payloads["tail_abel_join"]
    if join.get("counts", {}).get("signed_bulk_targets") != 1:
        raise RuntimeError("tail-Abel signed target drifted")
    if "Phi_B=P_B+P_cross" not in join.get(
        "symbolic_certificate", {}
    ).get("signed_bulk_remainder", ""):
        raise RuntimeError("tail-Abel Phi_B formula drifted")

    dyadic = payloads["joined_dyadic"].get("exact", {})
    if "H_(r,x)=-s_*'*H_(r+1)+D_r" not in dyadic.get(
        "five_current_ladder", ""
    ):
        raise RuntimeError("joined-dyadic five-current ladder drifted")
    if "q_(n,x)/q_n=-s_*'*log(n)+delta_n" not in dyadic.get(
        "coefficient_chain", ""
    ):
        raise RuntimeError("joined-dyadic coefficient current drifted")

    prefix_text = json.dumps(payloads["normalized_prefix"], sort_keys=True)
    for marker in ("gamma_n=q_(n,x)/q_n", "G_(k,x)=sum_(n=1)^k"):
        if marker not in prefix_text:
            raise RuntimeError(f"normalized-prefix marker missing: {marker}")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCE_PATHS.items()
        },
        "source_kinds": {
            key: payload.get("kind", "") for key, payload in payloads.items()
        },
        "imported_signed_bulk_bounds": 0,
    }


def symbolic_certificate() -> dict:
    # Verify the moment ladder directly on three generic carriers.
    kappa_0, s_prime = sp.symbols("kappa_0 s_prime")
    ell = sp.symbols("ell_1:4", real=True)
    z = sp.symbols("z_1:4")
    delta = sp.symbols("delta_1:4")
    z_x = tuple(
        (kappa_0 - s_prime * ell[index] + delta[index]) * z[index]
        for index in range(3)
    )

    def moment(power: int) -> sp.Expr:
        return sum(ell[index] ** power * z[index] for index in range(3))

    def correction(power: int) -> sp.Expr:
        return sum(
            ell[index] ** power * delta[index] * z[index]
            for index in range(3)
        )

    def differentiated_moment(power: int) -> sp.Expr:
        return sum(
            ell[index] ** power * z_x[index] for index in range(3)
        )

    h_0, h_1, h_2 = moment(0), moment(1), moment(2)
    d_0, d_1 = correction(0), correction(1)
    h_0_x = differentiated_moment(0)
    h_1_x = differentiated_moment(1)
    if sp.expand(h_0_x - (kappa_0 * h_0 - s_prime * h_1 + d_0)) != 0:
        raise RuntimeError("H_0 moment ladder failed")
    if sp.expand(h_1_x - (kappa_0 * h_1 - s_prime * h_2 + d_1)) != 0:
        raise RuntimeError("H_1 moment ladder failed")

    log_n = sp.symbols("log_N", real=True)
    r_b = log_n * h_0 - h_1
    r_b_x_direct = log_n * h_0_x - h_1_x
    r_b_x_closed = (
        kappa_0 * r_b
        - s_prime * (log_n * h_1 - h_2)
        + log_n * d_0
        - d_1
    )
    if sp.expand(r_b_x_direct - r_b_x_closed) != 0:
        raise RuntimeError("R_B five-current closure failed")

    # Change to the terminal-centered logarithmic-distance basis.  This is
    # algebraically equivalent but avoids subtracting two log(N)-scale sums.
    k_0 = h_0
    k_1 = log_n * h_0 - h_1
    k_2 = log_n**2 * h_0 - 2 * log_n * h_1 + h_2
    e_0 = d_0
    e_1 = log_n * d_0 - d_1
    kappa_n = kappa_0 - s_prime * log_n
    k_0_x_closed = kappa_n * k_0 + s_prime * k_1 + e_0
    k_1_x_closed = kappa_n * k_1 + s_prime * k_2 + e_1
    if sp.expand(h_0_x - k_0_x_closed) != 0:
        raise RuntimeError("centered K_0 ladder failed")
    if sp.expand(r_b_x_direct - k_1_x_closed) != 0:
        raise RuntimeError("centered K_1 ladder failed")

    # Audit the real/imaginary component version consumed by A_B.
    kr, ki, c, b = sp.symbols("kappa_R kappa_I c b", real=True)
    h0r, h0i, h1r, h1i = sp.symbols(
        "H0_R H0_I H1_R H1_I", real=True
    )
    d0r, d0i = sp.symbols("D0_R D0_I", real=True)
    f_x_r = kr * h0r - ki * h0i - c * h1r + b * h1i + d0r
    f_x_i = kr * h0i + ki * h0r - c * h1i - b * h1r + d0i
    complex_product = (kr + sp.I * ki) * (h0r + sp.I * h0i)
    complex_product -= (c + sp.I * b) * (h1r + sp.I * h1i)
    complex_product += d0r + sp.I * d0i
    if sp.simplify(sp.re(sp.expand_complex(complex_product)) - f_x_r) != 0:
        raise RuntimeError("F_B real derivative closure failed")
    if sp.simplify(sp.im(sp.expand_complex(complex_product)) - f_x_i) != 0:
        raise RuntimeError("F_B imaginary derivative closure failed")

    return {
        "normalized_coefficient_current": (
            "Let hat(z_n)=z_n/S_a and write eta_x/eta=i*omega_eta. "
            "With kappa_0=i*omega_eta-partial_x log(S_a), the exact "
            "coefficient law is partial_x hat(z_n)="
            "[kappa_0-s_*'*log(n)+delta_n]hat(z_n)."
        ),
        "five_moments": (
            "H_j^(B)=sum_(n=1)^B(log n)^j*hat(z_n), j=0,1,2; "
            "D_j^(B)=sum_(n=1)^B(log n)^j*delta_n*hat(z_n), j=0,1"
        ),
        "moment_ladder": (
            "H_(0,x)=kappa_0*H_0-s_*'*H_1+D_0; "
            "H_(1,x)=kappa_0*H_1-s_*'*H_2+D_1"
        ),
        "value_closure": (
            "hat(F_B)=H_0, partial_x hat(F_B)="
            "kappa_0*H_0-s_*'*H_1+D_0"
        ),
        "abel_moment_closure": (
            "hat(R_B)=log(N)*H_0-H_1; "
            "partial_x hat(R_B)=kappa_0*hat(R_B)"
            "-s_*'[log(N)*H_1-H_2]+log(N)*D_0-D_1"
        ),
        "centered_basis": (
            "lambda_n=log(N/n), K_j=sum_(n=1)^B lambda_n^j*hat(z_n) "
            "for j=0,1,2, and E_j=sum_(n=1)^B lambda_n^j*delta_n*"
            "hat(z_n) for j=0,1; K_0=H_0, K_1=log(N)H_0-H_1, "
            "K_2=log(N)^2H_0-2log(N)H_1+H_2, E_0=D_0, "
            "E_1=log(N)D_0-D_1"
        ),
        "centered_ladder": (
            "With chi_N=kappa_0-s_*'log(N), hat(F_B)=K_0, "
            "hat(R_B)=K_1, partial_x K_0=chi_N*K_0+s_*'*K_1+E_0, "
            "and partial_x K_1=chi_N*K_1+s_*'*K_2+E_1"
        ),
        "phi_substitution": (
            "Insert the displayed hat(F_B),partial_x hat(F_B),hat(R_B),"
            "partial_x hat(R_B) into "
            "A_B=Re(s_*'R_B)-b*u_N*Im(F_B), its exact first derivative, "
            "and Phi_B=h^-2[(X_T+X_B)A_(B,x)+X_B*A_(T,x)"
            "-(A_T+A_B)X_(B,x)-A_B*X_(T,x)]."
        ),
    }


def exact_payload() -> dict:
    return {
        "domain": (
            "Use the fixed q=1, fixed-N terminal/bulk partition of the "
            "tail-anchored Abel reduction, with B=N-M-1 and no cutoff "
            "change during x differentiation."
        ),
        "algebraic_closure": (
            "The full nonterminal input to Phi_B is determined exactly by "
            "the five complex currents H_0^(B),H_1^(B),H_2^(B),D_0^(B),"
            "D_1^(B), equivalently the centered currents K_0,K_1,K_2,E_0,"
            "E_1, the common scalar coefficients, and the four already "
            "certified tail coordinates. No individual prefix is needed "
            "after those moments have been formed."
        ),
        "dyadic_rejoin": (
            "The joined dyadic odd-prefix identities remain exact with their "
            "summation cutoff replaced by B: unique dyadic valuation gives "
            "H_0,H_1,H_2 and D_0,D_1 as heat-shifted odd-prefix sums, "
            "retaining every external odd phase and correction current."
        ),
        "route_gain": (
            "The remaining theorem can be sought either as a direct signed "
            "bound on Phi_B or as a quadratic inequality for this five-"
            "current vector after exact dyadic or p-adic rejoining."
        ),
        "sharp_target": (
            "Prove Phi_B<=1/400, preferably Phi_B<=1/800, for the actual "
            "physical five-current vector and endpoint-terminal tail."
        ),
        "nonpromotion": (
            "Finite-dimensional closure supplies sufficient data, not a "
            "sign. The existing negative interior currents, amplitude "
            "contraction, angular ordering, and separate prime-power block "
            "windings do not imply the required quadratic inequality without "
            "the joined Xi correlations."
        ),
        "pi_provenance": (
            "The logarithmic moments and dyadic valuation introduce no pi. "
            "All pi factors remain inherited from the completed-zeta and "
            "Riemann-Siegel normalization already present in the tail jets."
        ),
        "proof_boundary": (
            "This proves an exact five-current algebraic closure only. It "
            "proves no signed bound on Phi_B, bulk aggregate sign closure, "
            "complete retained or Xi-level current theorem, Abel gap, winding "
            "cap, contact exclusion, Q209, Lambda<=0, PF-infinity, RH, or "
            "prize-level conclusion."
        ),
    }


def build_rows(exact: dict, symbolic: dict) -> list[ClosureRow]:
    return [
        ClosureRow("ttfc_01_domain", "fixed-chart domain", "proved", "The five-current closure uses the same fixed terminal/bulk partition.", exact["domain"], "Cutoff transport remains separate."),
        ClosureRow("ttfc_02_coefficient", "coefficient current", "proved", "Common anchor and normalizer currents are isolated in one scalar kappa_0.", symbolic["normalized_coefficient_current"], "The correction current delta_n is retained exactly."),
        ClosureRow("ttfc_03_moments", "moment definition", "proved", "Only three logarithmic moments and two correction moments are required.", symbolic["five_moments"], "These are complex signed sums, not absolute moments."),
        ClosureRow("ttfc_04_ladder", "moment first jet", "proved", "The first two derivative moments close on the five-current vector.", symbolic["moment_ladder"], "No H_3 or D_2 enters the projective current."),
        ClosureRow("ttfc_05_value", "bulk value closure", "proved", "The normalized bulk value and derivative use H_0,H_1,D_0.", symbolic["value_closure"], "The common real normalization is already included."),
        ClosureRow("ttfc_06_abel", "bulk centered-moment closure", "proved", "The Abel moment and derivative use all five currents and no more.", symbolic["abel_moment_closure"] + "; " + symbolic["centered_ladder"], "The centered basis removes artificial log(N)-scale cancellation."),
        ClosureRow("ttfc_07_phi", "tail-joined substitution", "proved", "The signed remainder Phi_B is an explicit quadratic expression in the five currents and tail jets.", symbolic["phi_substitution"], "No sign is assigned."),
        ClosureRow("ttfc_08_closure", "finite-dimensional closure", "proved", "The nonterminal join no longer requires a stored prefix path.", exact["algebraic_closure"], "Computing or bounding the moments is still arithmetic work."),
        ClosureRow("ttfc_09_dyadic", "multiplicative rejoin", "proved", "The existing joined dyadic representation applies at cutoff B.", exact["dyadic_rejoin"], "Separate block winding is not promoted."),
        ClosureRow("ttfc_10_target", "signed target", "open", "The exact five-current vector must satisfy the terminal reserve inequality.", exact["sharp_target"], "No bulk sign closure is claimed."),
        ClosureRow("ttfc_11_guard", "nonpromotion guard", "proved", "Previously proved local structure does not by itself sign the joined quadratic form.", exact["nonpromotion"], "The actual joined Xi correlation is essential."),
        ClosureRow("ttfc_12_handoff", "analytic handoff", "open", "Seek a direct or multiplicatively rejoined five-current inequality.", exact["route_gain"], exact["proof_boundary"]),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    symbolic = payload["symbolic_certificate"]
    return f"""# Terminal-Tail-Anchored Five-Current Bulk Closure

Date: 2026-07-31

Status: exact algebraic closure with `0 signed bulk bounds`. This is not a
proof artifact; `Phi_B` and RH remain open.

## Coefficient Current

{exact['domain']}

```text
{symbolic['normalized_coefficient_current']}
```

Define the five truncated complex currents

```text
{symbolic['five_moments']}
```

Their required first derivatives satisfy

```text
{symbolic['moment_ladder']}
```

## Value and Abel Moment

```text
{symbolic['value_closure']}
{symbolic['abel_moment_closure']}
```

The `log(N)` coefficient is constant because the calculation stays inside one
fixed-`N` chart. These formulas are equalities of complex currents and retain
all cancellation.

The equivalent terminal-centered basis is better conditioned:

```text
{symbolic['centered_basis']}
{symbolic['centered_ladder']}
```

## Terminal Join

```text
{symbolic['phi_substitution']}
```

{exact['algebraic_closure']}

This is stronger bookkeeping than the `O(B)` Abel-prefix form: once the five
moments are supplied, evaluating `Phi_B` takes constant additional work.

## Multiplicative Handoff

{exact['dyadic_rejoin']}

{exact['route_gain']}

The sharp open inequality is

```text
{exact['sharp_target']}
```

{exact['nonpromotion']}

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
    exact = exact_payload()
    rows = build_rows(exact, symbolic)
    return {
        "kind": STEM,
        "date": "2026-07-31",
        "status": (
            "exact terminal-tail-anchored five-current bulk closure; "
            "Phi_B bound open; no bulk aggregate sign closure, Xi-level "
            "current theorem, Lambda<=0, or RH"
        ),
        "proof_boundary": exact["proof_boundary"],
        "source_audit": source,
        "symbolic_certificate": symbolic,
        "exact": exact,
        "counts": {
            "rows": len(rows),
            "normalized_coefficient_currents": 1,
            "complex_moments": 5,
            "moment_derivative_identities": 3,
            "centered_basis_transforms": 1,
            "centered_derivative_identities": 2,
            "five_current_algebraic_closures": 1,
            "multiplicative_rejoins": 1,
            "signed_bulk_targets": 1,
            "signed_bulk_bounds": 0,
            "bulk_aggregate_sign_closures": 0,
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
        "wrote terminal-tail five-current bulk closure: "
        f"{payload['counts']['rows']} rows, 5 complex currents, "
        "1 algebraic closure, 0 signed bulk bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
