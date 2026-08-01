#!/usr/bin/env python3
"""Build the explicit five-moment real quadratic reduction for Phi_B."""

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
    "contiguous_terminal_tail_anchored_five_moment_explicit_phi_"
    "quadratic_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "abel_join": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_abel_cross_current_"
        "reduction.json"
    ),
    "five_current": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_current_bulk_closure_"
        "reduction.json"
    ),
    "mangoldt_normal_form": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_current_mangoldt_"
        "normal_form_gate.json"
    ),
}


@dataclass(frozen=True)
class ReductionRow:
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
    abel = payloads["abel_join"]
    if abel.get("counts", {}).get("signed_bulk_targets") != 1:
        raise RuntimeError("Abel-join target drifted")
    if abel.get("counts", {}).get("bulk_aggregate_closures") != 0:
        raise RuntimeError("Abel-join proof boundary drifted")
    if "Phi_B=P_B+P_cross" not in abel.get(
        "symbolic_certificate", {}
    ).get("signed_bulk_remainder", ""):
        raise RuntimeError("Abel-join Phi_B identity drifted")

    five = payloads["five_current"]
    if five.get("counts", {}).get("complex_moments") != 5:
        raise RuntimeError("five-current source drifted")
    if five.get("counts", {}).get("signed_bulk_bounds") != 0:
        raise RuntimeError("five-current proof boundary drifted")

    mangoldt = payloads["mangoldt_normal_form"]
    if mangoldt.get("counts", {}).get("correction_free_complex_moments") != 5:
        raise RuntimeError("Mangoldt five-moment source drifted")
    if mangoldt.get("counts", {}).get("signed_type_ii_bounds") != 0:
        raise RuntimeError("Mangoldt proof boundary drifted")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCE_PATHS.items()
        },
        "source_kinds": {
            key: payload.get("kind", "") for key, payload in payloads.items()
        },
        "imported_terminal_reserve": "P_T<-1/400",
        "imported_phi_bounds": 0,
    }


def _real(expr: sp.Expr) -> sp.Expr:
    return sp.expand(expr).as_real_imag()[0]


def _imag(expr: sp.Expr) -> sp.Expr:
    return sp.expand(expr).as_real_imag()[1]


def _dot(left: sp.Matrix, right: sp.Matrix) -> sp.Expr:
    return sp.expand((left.T * right)[0])


def _real_row(vector: sp.Matrix) -> sp.Matrix:
    return sp.Matrix(
        [_real(entry) for entry in vector]
        + [-_imag(entry) for entry in vector]
    )


def symbolic_certificate() -> dict:
    log_n = sp.symbols("L_N", real=True)
    rho_1_r, rho_1_i, rho_2_r, rho_2_i = sp.symbols(
        "rho_1_R rho_1_I rho_2_R rho_2_I", real=True
    )
    rho_1_x_r, rho_1_x_i, rho_2_x_r, rho_2_x_i = sp.symbols(
        "rho_1_x_R rho_1_x_I rho_2_x_R rho_2_x_I", real=True
    )
    c, b, c_x, b_x = sp.symbols("c b c_x b_x", real=True)
    kappa_r, kappa_i = sp.symbols("kappa_R kappa_I", real=True)
    u_n, u_n_x = sp.symbols("u_N u_N_x", real=True)

    rho_1 = rho_1_r + sp.I * rho_1_i
    rho_2 = rho_2_r + sp.I * rho_2_i
    rho_1_x = rho_1_x_r + sp.I * rho_1_x_i
    rho_2_x = rho_2_x_r + sp.I * rho_2_x_i
    s_prime = c + sp.I * b
    s_second = c_x + sp.I * b_x
    kappa_n = kappa_r + sp.I * kappa_i

    v_0 = sp.Matrix([1, rho_1, rho_2, 0, 0])
    v_1 = sp.Matrix(
        [
            log_n,
            log_n * rho_1 - 1,
            log_n * rho_2 - rho_1,
            -rho_2,
            0,
        ]
    )
    v_2 = sp.Matrix(
        [
            log_n**2,
            log_n**2 * rho_1 - 2 * log_n,
            log_n**2 * rho_2 - 2 * log_n * rho_1 + 1,
            rho_1 - 2 * log_n * rho_2,
            rho_2,
        ]
    )
    e_0 = sp.Matrix([0, rho_1_x, rho_2_x, 0, 0])
    e_1 = sp.Matrix(
        [
            0,
            log_n * rho_1_x,
            log_n * rho_2_x - rho_1_x,
            -rho_2_x,
            0,
        ]
    )

    q_0 = kappa_n * v_0 + s_prime * v_1 + e_0
    q_1 = kappa_n * v_1 + s_prime * v_2 + e_1
    a_0 = s_prime * v_1 + sp.I * b * u_n * v_0
    n_0 = (
        s_second * v_1
        + s_prime * q_1
        + sp.I * (b_x * u_n + b * u_n_x) * v_0
        + sp.I * b * u_n * q_0
    )

    g_r = sp.symbols("G0_R:5", real=True)
    g_i = sp.symbols("G0_I:5", real=True)
    g_complex = sp.Matrix(
        [g_r[index] + sp.I * g_i[index] for index in range(5)]
    )
    g_real = sp.Matrix((*g_r, *g_i))

    k_0 = _dot(v_0, g_complex)
    k_1 = _dot(v_1, g_complex)
    k_2 = _dot(v_2, g_complex)
    correction_0 = _dot(e_0, g_complex)
    correction_1 = _dot(e_1, g_complex)
    k_0_x = _dot(q_0, g_complex)
    k_1_x = _dot(q_1, g_complex)

    if sp.simplify(
        k_0_x
        - (kappa_n * k_0 + s_prime * k_1 + correction_0)
    ) != 0:
        raise RuntimeError("K_0 derivative coefficient vector failed")
    if sp.simplify(
        k_1_x
        - (kappa_n * k_1 + s_prime * k_2 + correction_1)
    ) != 0:
        raise RuntimeError("K_1 derivative coefficient vector failed")

    x_bulk = _real(k_0)
    x_bulk_x = _real(k_0_x)
    a_bulk_direct = _real(s_prime * k_1) - b * u_n * _imag(k_0)
    a_bulk = _real(_dot(a_0, g_complex))
    if sp.simplify(a_bulk_direct - a_bulk) != 0:
        raise RuntimeError("bulk scalar coefficient vector failed")

    a_bulk_x_direct = (
        c_x * _real(k_1)
        + c * _real(k_1_x)
        - b_x * _imag(k_1)
        - b * _imag(k_1_x)
        - (b_x * u_n + b * u_n_x) * _imag(k_0)
        - b * u_n * _imag(k_0_x)
    )
    a_bulk_x = _real(_dot(n_0, g_complex))
    if sp.simplify(a_bulk_x_direct - a_bulk_x) != 0:
        raise RuntimeError("bulk scalar derivative coefficient vector failed")

    v_row = _real_row(v_0)
    q_row = _real_row(q_0)
    a_row = _real_row(a_0)
    n_row = _real_row(n_0)
    for name, direct, row in (
        ("V", x_bulk, v_row),
        ("Q", x_bulk_x, q_row),
        ("A", a_bulk, a_row),
        ("N", a_bulk_x, n_row),
    ):
        if sp.simplify(direct - _dot(row, g_real)) != 0:
            raise RuntimeError(f"real observation row failed: {name}")

    x_t, a_t, x_t_x, a_t_x = sp.symbols(
        "X_T A_T X_T_x A_T_x", real=True
    )
    tail_current = x_t * a_t_x - a_t * x_t_x
    retained_current = (
        (x_t + x_bulk) * (a_t_x + a_bulk_x)
        - (a_t + a_bulk) * (x_t_x + x_bulk_x)
    )
    phi_numerator = sp.expand(retained_current - tail_current)
    phi_observation = sp.expand(
        (x_t + x_bulk) * a_bulk_x
        + x_bulk * a_t_x
        - (a_t + a_bulk) * x_bulk_x
        - a_bulk * x_t_x
    )
    if sp.simplify(phi_numerator - phi_observation) != 0:
        raise RuntimeError("Phi_B observation formula failed")

    matrix = sp.expand(
        (
            v_row * n_row.T
            + n_row * v_row.T
            - a_row * q_row.T
            - q_row * a_row.T
        )
        / 2
    )
    linear = (
        x_t * n_row
        + a_t_x * v_row
        - a_t * q_row
        - x_t_x * a_row
    )
    quadratic = sp.expand(_dot(g_real, matrix * g_real) + _dot(linear, g_real))
    if sp.simplify(phi_numerator - quadratic) != 0:
        raise RuntimeError("explicit real quadratic form failed")
    if matrix != matrix.T:
        raise RuntimeError("real quadratic matrix is not symmetric")

    factor = sp.Matrix.hstack(v_row, n_row, a_row, q_row)
    core = sp.Matrix(
        [
            [0, sp.Rational(1, 2), 0, 0],
            [sp.Rational(1, 2), 0, 0, 0],
            [0, 0, 0, -sp.Rational(1, 2)],
            [0, 0, -sp.Rational(1, 2), 0],
        ]
    )
    if sp.simplify(matrix - factor * core * factor.T) != sp.zeros(10):
        raise RuntimeError("rank-four quadratic factorization failed")

    witness = {
        log_n: 2,
        rho_1_r: 0,
        rho_1_i: 0,
        rho_2_r: 0,
        rho_2_i: 0,
        rho_1_x_r: 0,
        rho_1_x_i: 0,
        rho_2_x_r: 0,
        rho_2_x_i: 0,
        c: 0,
        b: -sp.Rational(1, 2),
        c_x: 0,
        b_x: 0,
        kappa_r: 1,
        kappa_i: 0,
        u_n: 1,
        u_n_x: 1,
    }
    observation = sp.Matrix.vstack(
        v_row.T, q_row.T, a_row.T, n_row.T
    )
    witness_observation = observation.subs(witness)
    witness_minor = witness_observation[:, [0, 1, 5, 6]].det()
    if witness_observation.rank() != 4 or witness_minor != sp.Rational(5, 16):
        raise RuntimeError("generic rank-four witness failed")

    total_value_fibre = sp.expand(
        phi_observation.subs(x_t, -x_bulk)
    )
    expected_total_value_fibre = sp.expand(
        x_bulk * a_t_x
        - (a_t + a_bulk) * x_bulk_x
        - a_bulk * x_t_x
    )
    if sp.simplify(total_value_fibre - expected_total_value_fibre) != 0:
        raise RuntimeError("total-real-value fibre failed")

    double_contact = sp.expand(
        total_value_fibre.subs(a_t, -a_bulk)
    )
    negative_tail_current = sp.expand(
        -(x_t * a_t_x - a_t * x_t_x).subs(
            {x_t: -x_bulk, a_t: -a_bulk}
        )
    )
    if sp.simplify(double_contact - negative_tail_current) != 0:
        raise RuntimeError("double-contact obstruction failed")

    return {
        "base_coefficient_vectors": (
            "Let L=log(N), G=(G_0,...,G_4)^T. Then "
            "v_0=(1,rho_1,rho_2,0,0), "
            "v_1=(L,Lrho_1-1,Lrho_2-rho_1,-rho_2,0), "
            "v_2=(L^2,L^2rho_1-2L,L^2rho_2-2Lrho_1+1,"
            "rho_1-2Lrho_2,rho_2), "
            "e_0=(0,rho_(1,x),rho_(2,x),0,0), and "
            "e_1=(0,Lrho_(1,x),Lrho_(2,x)-rho_(1,x),"
            "-rho_(2,x),0)."
        ),
        "derivative_coefficient_vectors": (
            "q_0=chi_N*v_0+s_*'*v_1+e_0; "
            "q_1=chi_N*v_1+s_*'*v_2+e_1; "
            "a_0=s_*'*v_1+i*b*u_N*v_0; "
            "n_0=s_*''*v_1+s_*'*q_1+i*(b_x*u_N+b*u_(N,x))*v_0"
            "+i*b*u_N*q_0."
        ),
        "four_real_observables": (
            "V=Re(v_0 dot G)=X_B, Q=Re(q_0 dot G)=X_(B,x), "
            "A=Re(a_0 dot G)=A_B, and "
            "N=Re(n_0 dot G)=A_(B,x)."
        ),
        "explicit_phi": (
            "h^2*Phi_B=(X_T+V)N+V*A_(T,x)-(A_T+A)Q-A*X_(T,x)."
        ),
        "real_quadratic_form": (
            "For g=(Re G_0,...,Re G_4,Im G_0,...,Im G_4)^T and "
            "R(z)=(Re z_0,...,Re z_4,-Im z_0,...,-Im z_4)^T, put "
            "v=R(v_0), q=R(q_0), a=R(a_0), n=R(n_0). Then "
            "h^2*Phi_B=g^T M g+ell_T^T g, "
            "M=(v n^T+n v^T-a q^T-q a^T)/2 and "
            "ell_T=X_T n+A_(T,x)v-A_T q-X_(T,x)a."
        ),
        "rank_factorization": (
            "M=[v n a q] J [v n a q]^T with "
            "J=(1/2)*diag_offdiag(+1,-1); hence rank(M)<=4. "
            "The observation map g->(V,Q,A,N) has a real kernel of "
            "dimension at least 6, and Phi_B vanishes on that kernel."
        ),
        "generic_inertia_guard": (
            "An exact algebraic specialization L=2,rho_1=rho_2="
            "rho_(1,x)=rho_(2,x)=0,s_*'=-i/2,s_*''=0,chi_N=1,"
            "u_N=u_(N,x)=1 has observation rank 4 and minor 5/16. "
            "By the displayed congruence its M has inertia (2,2,6). "
            "This is a generic nonpromotion guard, not a physical Xi state."
        ),
        "total_value_fibre": (
            "On X_T+V=0, h^2*Phi_B=V*A_(T,x)-(A_T+A)Q-"
            "A*X_(T,x), with no division by the total value."
        ),
        "double_contact_obstruction": (
            "If also A_T+A=0, then Phi_B=-P_T>1/400. Therefore the "
            "desired bound must source-specifically exclude a simultaneous "
            "retained value/centered-scalar contact; it cannot assume it away."
        ),
        "witness_minor": "5/16",
        "generic_inertia": [2, 2, 6],
    }


def exact_payload() -> dict:
    return {
        "domain": (
            "Fix the q=1, L>=50, fixed-N chart and admissible growing "
            "terminal block of the source reductions. Put B=N-M-1 and "
            "differentiate with N,M,B fixed. All G_r stop at B and use the "
            "same eta/S_a normalization as the tail."
        ),
        "coefficient_meaning": (
            "The five base vectors convert the correction-free moments to "
            "K_0,K_1,K_2,E_0,E_1. The four derived vectors then give the "
            "bulk value, value derivative, centered scalar, and centered-"
            "scalar derivative without storing a prefix path."
        ),
        "dimension_reduction": (
            "Although five complex moments have ten real components, Phi_B "
            "depends on them only through four real observations. The exact "
            "quadratic matrix has rank at most four and its tail-linear term "
            "lies in the same four-row span."
        ),
        "sharp_arithmetic_target": (
            "For the actual Xi moment vector g_B and certified tail jets, "
            "prove g_B^T M g_B+ell_T^T g_B<=h^2/400, preferably <=h^2/800."
        ),
        "type_ii_handoff": (
            "Insert the exact balanced Mangoldt or Vaughan decomposition "
            "into the four observations V,Q,A,N, or equivalently expand the "
            "rank-four form as its original two-carrier kernel. Keep the "
            "endpoint, complete square, both wings, and all transpose and "
            "Hermitian cross terms before taking absolute values."
        ),
        "falsification_program": (
            "Any proposed signed criterion must survive X_T+V=0, W_0=0, "
            "zero observation rows, p=+-1, removable p=+-1/2, prime edges, "
            "square and cube transports, n=64, q=1, and adjacent cutoffs. "
            "The simultaneous value/scalar contact is an obstruction to be "
            "excluded by the theorem, not an admissible denominator."
        ),
        "proof_boundary": (
            "This proves an exact explicit rank-four real quadratic normal "
            "form, fibre identities, and a generic indefiniteness guard. It "
            "proves no upper bound for Phi_B, endpoint-coupled Type-I/II or "
            "Vaughan estimate, contact exclusion, retained aggregate sign, "
            "Xi residual transfer, Q209, cofinal descendant theorem, "
            "Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }


def build_rows(exact: dict, symbolic: dict) -> list[ReductionRow]:
    return [
        ReductionRow("ttqp_01_domain", "fixed-chart domain", "proved", "The explicit form uses the certified terminal/bulk partition.", exact["domain"], "No cutoff transport is folded into an x derivative."),
        ReductionRow("ttqp_02_base", "base coefficient vectors", "proved", "Five explicit vectors map G_0 through G_4 to K_0,K_1,K_2,E_0,E_1.", symbolic["base_coefficient_vectors"], "All correction coefficients and derivatives remain complex and exact."),
        ReductionRow("ttqp_03_derivatives", "derived coefficient vectors", "proved", "Four derived vectors include every required first derivative.", symbolic["derivative_coefficient_vectors"], "No second correction current is introduced."),
        ReductionRow("ttqp_04_observations", "four real observations", "proved", "The ten real moment coordinates enter Phi_B through four real functionals.", symbolic["four_real_observables"], exact["coefficient_meaning"]),
        ReductionRow("ttqp_05_phi", "explicit Phi identity", "proved", "The tail/bulk join is one displayed polynomial in the four observations.", symbolic["explicit_phi"], "No value or scalar is divided out."),
        ReductionRow("ttqp_06_quadratic", "real quadratic form", "proved", "Phi_B is one explicit symmetric quadratic form plus a tail-linear term.", symbolic["real_quadratic_form"], "This is an equality, not a sign estimate."),
        ReductionRow("ttqp_07_rank", "rank factorization", "proved", "The moment quadratic has rank at most four.", symbolic["rank_factorization"], exact["dimension_reduction"]),
        ReductionRow("ttqp_08_inertia", "generic inertia guard", "proved", "A coefficient-blind semidefinite promotion is impossible.", symbolic["generic_inertia_guard"], "The witness is algebraic and is not asserted to be attained by Xi.", {"minor": symbolic["witness_minor"], "inertia": symbolic["generic_inertia"]}),
        ReductionRow("ttqp_09_kernel", "observation kernel", "proved", "Every common null direction of the four observations is invisible to Phi_B.", symbolic["rank_factorization"], "This does not bound the visible four-dimensional quotient."),
        ReductionRow("ttqp_10_value_fibre", "total real-value fibre", "proved", "The formula remains division-free when the retained real value vanishes.", symbolic["total_value_fibre"], "The crossing orientation remains unsigned."),
        ReductionRow("ttqp_11_contact", "double-contact obstruction", "proved", "A simultaneous retained value/scalar contact violates the desired Phi_B bound.", symbolic["double_contact_obstruction"], "The actual source must exclude this fibre."),
        ReductionRow("ttqp_12_target", "signed arithmetic target", "open", "The remaining theorem is one explicit rank-four inequality.", exact["sharp_arithmetic_target"], "No Type-I/II constant is supplied."),
        ReductionRow("ttqp_13_handoff", "Type-I/II handoff", "open", "The four observations must be estimated jointly in their endpoint-composed form.", exact["type_ii_handoff"], exact["proof_boundary"]),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    symbolic = payload["symbolic_certificate"]
    return f"""# Terminal-Tail Explicit Five-Moment Phi Quadratic Reduction

Date: 2026-07-31

Status: exact rank-four quadratic reduction with `0 Phi_B bounds`. This is
not a proof artifact; the signed arithmetic estimate and RH remain open.

## Five-Moment Coefficients

{exact['domain']}

```text
{symbolic['base_coefficient_vectors']}

{symbolic['derivative_coefficient_vectors']}
```

{exact['coefficient_meaning']}

Define the four real observations by

```text
{symbolic['four_real_observables']}
```

## Exact Joined Form

```text
{symbolic['explicit_phi']}

{symbolic['real_quadratic_form']}
```

{exact['dimension_reduction']}

The exact factorization is

```text
{symbolic['rank_factorization']}
```

## Indefiniteness And Fibres

{symbolic['generic_inertia_guard']}

This guard forbids a generic positivity or negativity promotion. It does not
replace analysis of the actual Xi moment vector.

```text
{symbolic['total_value_fibre']}

{symbolic['double_contact_obstruction']}
```

The second identity explains why the zero fibre cannot be discarded: the
desired inequality must itself prevent the simultaneous contact.

## Open Arithmetic Target

```text
{exact['sharp_arithmetic_target']}
```

{exact['type_ii_handoff']}

{exact['falsification_program']}

## Boundary

No new `pi` occurs in the coefficient vectors or realification. Existing
`pi` factors retain their completed-zeta/Riemann--Siegel origin.

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
            "exact terminal-tail five-moment explicit rank-four Phi_B "
            "quadratic reduction; signed Phi_B bound open; no retained "
            "aggregate sign, Xi-level current theorem, Lambda<=0, or RH"
        ),
        "proof_boundary": exact["proof_boundary"],
        "source_audit": source,
        "symbolic_certificate": symbolic,
        "exact": exact,
        "counts": {
            "rows": len(rows),
            "complex_correction_free_moments": 5,
            "real_moment_coordinates": 10,
            "complex_base_coefficient_vectors": 5,
            "real_bulk_observations": 4,
            "explicit_phi_quadratic_forms": 1,
            "quadratic_rank_upper_bound": 4,
            "generic_positive_directions": 2,
            "generic_negative_directions": 2,
            "generic_null_directions": 6,
            "division_free_value_fibres": 1,
            "double_contact_obstructions": 1,
            "signed_phi_targets": 1,
            "signed_phi_bounds": 0,
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
        "wrote explicit five-moment Phi_B quadratic reduction: "
        f"{payload['counts']['rows']} rows, rank <=4, "
        "1 double-contact obstruction, 0 signed Phi_B bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
