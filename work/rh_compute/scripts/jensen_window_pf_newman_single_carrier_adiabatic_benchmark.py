#!/usr/bin/env python3
"""Build an exact single-carrier benchmark for the adiabatic route."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_single_carrier_adiabatic_benchmark"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{STEM}.md"
SUCCESSOR_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma.json"
)
WRONSKIAN_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_critical_component_wronskian_gate.json"
)
DATE = "2026-07-25"


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def rows() -> list[dict]:
    return [
        {
            "id": "single_carrier_first_jet",
            "kind": "exact_identity",
            "claim": (
                "For E=A*exp(i*theta), A>0, r=E_x/E=u+i*v, "
                "and J=2*Re(E), the scaled first jet is "
                "V_ell=2*A*B*(cos(theta),sin(theta))^T."
            ),
            "formula": (
                "B=[[1,0],[u/ell,-v/ell]], ell>0"
            ),
        },
        {
            "id": "phase_matrix_determinant",
            "kind": "exact_identity",
            "claim": (
                "The carrier first-jet map is invertible exactly when "
                "its phase speed v is nonzero."
            ),
            "formula": (
                "det(B)=-v/ell; "
                "||B||_F^2=1+(u^2+v^2)/ell^2"
            ),
        },
        {
            "id": "carrier_clearance_floor",
            "kind": "exact_inequality",
            "claim": (
                "The determinant/Frobenius singular-value bound gives "
                "a phase-uniform first-jet clearance."
            ),
            "formula": (
                "||V_ell||_2 >= "
                "2*A*|v|/sqrt(ell^2+u^2+v^2)"
            ),
        },
        {
            "id": "carrier_heat_jet",
            "kind": "exact_identity",
            "claim": (
                "If E_t=-E_xx, then J_t=-J_xx and the logarithmic "
                "second and third x-jets determine the heat displacement."
            ),
            "formula": (
                "E_xx/E=r^2+r_x; "
                "E_xxx/E=r^3+3*r*r_x+r_xx"
            ),
        },
        {
            "id": "uniform_log_jet_hypotheses",
            "kind": "conditional_definition",
            "claim": (
                "Assume |v|>=c*ell, |r|<=C1*ell, "
                "|r_x|<=C2*ell^2, and |r_xx|<=C3*ell^3, "
                "with c>0 and fixed nonnegative C1,C2,C3."
            ),
            "formula": (
                "B2=C1^2+C2; "
                "B3=C1^3+3*C1*C2+C3"
            ),
        },
        {
            "id": "adiabatic_condition_number",
            "kind": "exact_conditional_inequality",
            "claim": (
                "Under the uniform logarithmic-jet hypotheses, the "
                "heat-jet condition number grows at most quadratically "
                "in the local frequency scale."
            ),
            "formula": (
                "||partial_t V_ell||_2/||V_ell||_2 "
                "<=K*ell^2; "
                "K=sqrt(1+C1^2)*sqrt(B2^2+B3^2)/c"
            ),
        },
        {
            "id": "cofinal_shell_cost",
            "kind": "exact_limit",
            "claim": (
                "For R_j=j+38, delta_j=1/(5*j*(j+1)), and "
                "ell_j=max(1,log(R_j/(4*pi))), the hypothetical "
                "single-carrier collar cost tends to zero."
            ),
            "formula": (
                "delta_j*ell_j^2 -> 0 as j->infinity"
            ),
        },
        {
            "id": "gronwall_transport",
            "kind": "exact_conditional_inequality",
            "claim": (
                "A pointwise condition-number bound on a collar gives "
                "multiplicative rather than additive transport and "
                "cannot create an origin contact from a nonzero old jet."
            ),
            "formula": (
                "||V(t,x)||_2 >= "
                "exp(-K*ell(x)^2*|t-t_j|)*||V(t_j,x)||_2"
            ),
        },
        {
            "id": "four_carrier_obstruction",
            "kind": "exact_countermodel",
            "claim": (
                "Positive amplitudes and distinct ordered component "
                "phase speeds do not imply a phase-speed floor for "
                "their sum."
            ),
            "formula": (
                "E_toy=exp(-3ix)-exp(-4ix)+i*exp(-ix)"
                "-(i/2)*exp(-2ix); "
                "E_toy(0)=i/2, E_toy'(0)=i, "
                "Re(E_toy)''(0)=7"
            ),
            "consequence": (
                "Re(E_toy)(0)=Re(E_toy'(0))=0 and "
                "Im(E_toy'(0)/E_toy(0))=0."
            ),
        },
        {
            "id": "xi_crossing_small_ball_handoff",
            "kind": "open_xi_target",
            "claim": (
                "The actual corrected Riemann-Siegel main needs a "
                "crossing-restricted arithmetic small-ball theorem, "
                "not a generic componentwise phase-speed argument."
            ),
            "formula": (
                "X^2+(U/L)^2 > "
                "8000000*exp(-3*L/2) on the residual critical layer"
            ),
            "status": "open and not proved",
        },
    ]


def render_note(artifact: dict) -> str:
    return f"""# Newman Single-Carrier Adiabatic Benchmark

Date: {DATE}

Status: exact conditional thought-experiment theorem and rejection guard; not a proof of `Lambda<=0` or RH.

## Hypothetical Carrier

Let

```text
E=A*exp(i*theta), A>0,
r=E_x/E=u+i*v,
J=2*Re(E),
V_ell=(J,J_x/ell), ell>0.
```

Then

```text
V_ell=2*A [[1,0],[u/ell,-v/ell]]
              (cos(theta),sin(theta))^T.          (1)
```

The matrix in (1) has determinant `-v/ell` and squared
Frobenius norm `1+(u^2+v^2)/ell^2`. Therefore

```text
||V_ell||_2
 >=2*A*|v|/sqrt(ell^2+u^2+v^2).                 (2)
```

This is the clean hypothetical geometry: a nonzero phase
speed supplies a phase-independent first-jet floor.

## Heat Transport

If the carrier itself obeys `E_t=-E_xx`, then

```text
E_xx/E=r^2+r_x,
E_xxx/E=r^3+3*r*r_x+r_xx.                       (3)
```

Assume fixed constants `c>0`, `C1,C2,C3>=0` satisfy

```text
|v|>=c*ell,
|r|<=C1*ell,
|r_x|<=C2*ell^2,
|r_xx|<=C3*ell^3.                               (4)
```

Put

```text
B2=C1^2+C2,
B3=C1^3+3*C1*C2+C3,
K=sqrt(1+C1^2)*sqrt(B2^2+B3^2)/c.
```

Equations (2)-(4) give the exact conditional bound

```text
||partial_t V_ell||_2/||V_ell||_2<=K*ell^2.     (5)
```

Gronwall then transports a nonzero old jet multiplicatively:

```text
||V_ell(t,x)||_2
 >=exp(-K*ell(x)^2*|t-t_j|)*||V_ell(t_j,x)||_2. (6)
```

For the cofinal shells,

```text
delta_j=1/(5*j*(j+1)),
R_j=j+38,
ell_j=max(1,log(R_j/(4*pi))),
delta_j*ell_j^2 -> 0.                           (7)
```

Thus the ideal single-carrier thought experiment really does
make the successor increasingly adiabatic.

## Why It Is Not Xi Yet

The exact four-carrier model

```text
E_toy=exp(-3ix)-exp(-4ix)+i*exp(-ix)
      -(i/2)*exp(-2ix)
```

has positive amplitudes and distinct ordered phase speeds, but

```text
E_toy(0)=i/2,
E_toy'(0)=i,
Re(E_toy)(0)=Re(E_toy'(0))=0,
Re(E_toy)''(0)=7,
Im(E_toy'(0)/E_toy(0))=0.                       (8)
```

So summing individually good carriers can destroy the phase-speed
floor exactly at a real double zero. Generic positive amplitudes,
ordered speeds, and omitted-tail domination cannot establish (4).

## Xi Handoff

For the corrected Riemann-Siegel main `E=X+iY` and
`E'=U+iV`, the live theorem remains the crossing-restricted
arithmetic estimate

```text
X^2+(U/L)^2
 >8000000*exp(-3*L/2)                            (9)
```

through the residual `tL->0` layer, or an equivalent
Xi-specific restriction that excludes the negative Wronskian
direction on the crossing set. Equation (7) explains why such a
theorem would make the adiabatic successor work; it does not prove
that theorem.

## Proof Boundary

The carrier algebra, singular-value floor, logarithmic-jet
condition number, cofinal `delta_j*ell_j^2` limit, Gronwall
transport, and four-carrier countermodel are exact. The uniform
carrier hypotheses are not proved for Xi. The crossing small-ball
estimate, all-`j` successor theorem, `Lambda<=0`, RH, PF-infinity,
and the Clay prize remain open.

Machine-audited files:

```text
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```
"""


def main() -> int:
    artifact = {
        "kind": STEM,
        "date": DATE,
        "status": (
            "exact conditional single-carrier adiabatic benchmark "
            "with a four-carrier rejection guard and one open Xi "
            "small-ball handoff"
        ),
        "proof_boundary": (
            "The single-carrier identities, conditional O(ell^2) "
            "condition number, cofinal limit, Gronwall transport, "
            "and four-carrier countermodel are exact. The uniform "
            "carrier hypotheses and corrected Xi crossing small-ball "
            "estimate are not proved. This does not prove an all-j "
            "successor theorem, Lambda<=0, RH, or the Clay prize."
        ),
        "source_sha256": {
            "successor_lemma": file_hash(SUCCESSOR_RESULT),
            "critical_component_wronskian_gate": file_hash(
                WRONSKIAN_RESULT
            ),
        },
        "builder_sha256": file_hash(Path(__file__)),
        "rows": rows(),
        "exact_summary": {
            "row_count": 10,
            "exact_benchmark_reductions": 8,
            "four_carrier_obstructions": 1,
            "open_xi_small_ball_handoffs": 1,
            "cofinal_cost": "delta_j*ell_j^2->0",
            "generic_multi_carrier_promotion": False,
            "all_j_xi_theorem": False,
        },
        "next_exact_obligation": (
            "Prove an Xi-specific crossing-restricted lower bound for "
            "X^2+(U/L)^2 in the corrected tL->0 layer, or an "
            "equivalent phase-critical-value avoidance theorem, and "
            "combine it with explicit higher-jet and right-strip "
            "bounds."
        ),
    }
    RESULT_PATH.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    NOTE_PATH.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman single-carrier adiabatic benchmark: "
        "10 rows, 8 exact benchmark reductions, "
        "1 four-carrier obstruction, 1 open Xi small-ball handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
