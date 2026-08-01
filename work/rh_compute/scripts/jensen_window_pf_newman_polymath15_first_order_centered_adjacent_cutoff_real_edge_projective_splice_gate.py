#!/usr/bin/env python3
"""Build the adjacent-cutoff real-edge projective splice theorem."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "adjacent_cutoff_real_edge_projective_splice_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "critical_ray_edge": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "critical_ray_finite_height_real_edge_gate.json"
    ),
    "adjacent_recurrence": (
        REPO_ROOT
        / "outputs/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "adjacent_saddle_recurrence.md"
    ),
    "real_edge_source": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "real_edge_projective_current_gate.json"
    ),
}
H_DENOMINATOR = 72_000_000_000


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


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def source_audit() -> dict:
    missing = [str(path) for path in SOURCE_PATHS.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))

    source_text = {
        key: path.read_text(encoding="utf-8")
        for key, path in SOURCE_PATHS.items()
    }
    markers = {
        "critical_ray_edge": (
            "c=A+h*e_0, |e_0|<100",
            "d_edge=h*B+h^2*e_1, |e_1|<250",
            "a^2 J_edge/S_a^2<-3749/10000000<0",
        ),
        "adjacent_recurrence": (
            "kappa_(N+1)=-kappa_N",
            "g_(N+1)-g_N=kappa_N*J_a",
        ),
        "real_edge_source": (
            "C_edge=Re(V)",
            "B=-A'/(4pi)",
        ),
    }
    for key, required in markers.items():
        text = source_text[key]
        for marker in required:
            if marker not in text:
                raise RuntimeError(f"{key} source marker missing: {marker}")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCE_PATHS.items()
        },
        "published_source": {
            "citation": (
                "D.H.J. Polymath, Effective approximation of heat flow "
                "evolution of the Riemann xi function, and a new upper "
                "bound for the de Bruijn-Newman constant, arXiv:1904.12438"
            ),
            "location": "equation (53), C_0, first correction, and Proposition 6.2",
            "url": "https://arxiv.org/abs/1904.12438",
        },
    }


def symbolic_certificate() -> dict:
    p = sp.symbols("p", real=True)
    pi = sp.pi
    theta = (1 - p) / 2
    psi = pi * (p**2 / 2 + sp.Rational(3, 8))
    delta = pi * (p**2 / 2 - p + sp.Rational(3, 8))
    f = sp.cos(psi) / (2 * sp.cos(pi * p))
    a_source = sp.simplify(f - sp.cos(delta))
    b_source = sp.simplify(
        -sp.diff(f, p) / (4 * pi) + theta * sp.sin(delta) / 2
    )

    radical_plus = sp.sqrt(2 + sp.sqrt(2))
    radical_minus = sp.sqrt(2 - sp.sqrt(2))
    expected_a = -radical_plus / 4
    expected_b_minus = -3 * radical_minus / 16
    expected_b_plus = -radical_minus / 16
    endpoint_values = {
        "A_minus": sp.simplify(a_source.subs(p, -1)),
        "A_plus": sp.simplify(a_source.subs(p, 1)),
        "B_minus": sp.simplify(b_source.subs(p, -1)),
        "B_plus": sp.simplify(b_source.subs(p, 1)),
    }
    expected = {
        "A_minus": expected_a,
        "A_plus": expected_a,
        "B_minus": expected_b_minus,
        "B_plus": expected_b_plus,
    }
    for key, value in endpoint_values.items():
        if sp.simplify(value - expected[key]) != 0:
            raise RuntimeError(f"cutoff endpoint identity failed: {key}")

    leading_wedge = sp.simplify(
        endpoint_values["A_minus"] * endpoint_values["B_plus"]
        - endpoint_values["B_minus"] * endpoint_values["A_plus"]
    )
    if sp.simplify(leading_wedge + sp.sqrt(2) / 32) != 0:
        raise RuntimeError("leading cutoff wedge failed")

    h = sp.symbols("h", positive=True, real=True)
    e0m, e0p, e1m, e1p = sp.symbols(
        "e0_minus e0_plus e1_minus e1_plus", real=True
    )
    am = endpoint_values["A_minus"]
    ap = endpoint_values["A_plus"]
    bm = endpoint_values["B_minus"]
    bp = endpoint_values["B_plus"]
    cm = am + h * e0m
    cp = ap + h * e0p
    dm = h * bm + h**2 * e1m
    dp = h * bp + h**2 * e1p
    wedge = sp.expand(cm * dp - dm * cp)
    expected_wedge = sp.expand(
        h * am * (bp - bm)
        + h**2 * (am * (e1p - e1m) + e0m * bp - bm * e0p)
        + h**3 * (e0m * e1p - e1m * e0p)
    )
    if sp.simplify(wedge - expected_wedge) != 0:
        raise RuntimeError("finite cutoff wedge expansion failed")

    s = sp.symbols("s", real=True)
    cs = (1 - s) * cm + s * cp
    ds = (1 - s) * dm + s * dp
    join_current = sp.expand(cs * sp.diff(ds, s) - ds * sp.diff(cs, s))
    if sp.simplify(join_current - wedge) != 0:
        raise RuntimeError("affine projective-current identity failed")

    return {
        "endpoint_values": {
            "A_minus": "-sqrt(2+sqrt(2))/4",
            "A_plus": "-sqrt(2+sqrt(2))/4",
            "B_minus": "-3sqrt(2-sqrt(2))/16",
            "B_plus": "-sqrt(2-sqrt(2))/16",
        },
        "leading_wedge": (
            "A_minus*B_plus-B_minus*A_plus=-sqrt(2)/32"
        ),
        "finite_wedge": (
            "Delta_cut=c_minus*d_plus-d_minus*c_plus="
            "-sqrt(2)h/32+h^2[A(e1_plus-e1_minus)"
            "+e0_minus B_plus-B_minus e0_plus]"
            "+h^3[e0_minus e1_plus-e1_minus e0_plus]"
        ),
        "affine_join": (
            "v_s=(1-s)v_minus+s v_plus; "
            "c_s*d_(s)-d_s*c_(s)=Delta_cut"
        ),
    }


def arithmetic_certificate() -> dict:
    h0 = Fraction(1, H_DENOMINATOR)
    second_order = Fraction(500) * h0
    third_order = Fraction(50_000) * h0**2
    total_error = second_order + third_order
    if not total_error < Fraction(1, 600):
        raise RuntimeError("cutoff wedge error budget failed")
    if not Fraction(100) * h0 < Fraction(1, 12):
        raise RuntimeError("cutoff nonvanishing budget failed")

    return {
        "h_upper": fraction_text(h0),
        "endpoint_bounds": (
            "A<-1/3, |A|<1/2, |B_minus|<1, |B_plus|<1"
        ),
        "defect_bounds": "|e0_minus|,|e0_plus|<100; |e1_minus|,|e1_plus|<250",
        "second_order_coefficient": 500,
        "third_order_coefficient": 50_000,
        "error_at_h_upper": fraction_text(total_error),
        "error_cap": "1/600",
        "radical_margin": "sqrt(2)/32>1/24",
        "normalized_wedge": "Delta_cut/h<-1/25",
        "signed_wedge": "Delta_cut<-h/25<0",
        "nonvanishing": "c_minus<-1/4, c_plus<-1/4, c_s<-1/4",
    }


def exact_payload() -> dict:
    return {
        "domain": (
            "L>=50, 0<tL<=25, h=1/a<1/72000000000; at an integer "
            "cutoff a=m compare p_minus=-1 in the N=m-1 chart with "
            "p_plus=1 in the N=m chart"
        ),
        "physical_alignment": (
            "S_minus=kappa_(m-1)T_0, S_plus=kappa_mT_0=-S_minus; "
            "v_minus=(C_minus/S_minus,D_minus/S_minus), "
            "v_plus=(C_plus/S_plus,D_plus/S_plus). Real nonzero rescaling "
            "does not change the projective points, and c_minus,c_plus<0 "
            "select their common short lift."
        ),
        "finite_jets": (
            "c_sigma=A_sigma+h e0_sigma, |e0_sigma|<100; "
            "d_sigma=h B_sigma+h^2 e1_sigma, |e1_sigma|<250, "
            "sigma in {minus,plus}"
        ),
        "orientation": (
            "For v_s=(c_s,d_s), partial_s arg(c_s+i d_s)="
            "Delta_cut/(c_s^2+d_s^2)<0"
        ),
        "pi_provenance": (
            "The endpoint A and B formulas inherit pi from the completed-zeta "
            "and Riemann-Siegel C_0 normalization. Evaluation at p=+-1 "
            "reduces that source phase to square radicals; no geometric "
            "circle, polygon, or fitted constant is introduced."
        ),
    }


def build_rows(exact: dict, symbolic: dict, arithmetic: dict) -> list[GateRow]:
    return [
        GateRow(
            "acres_01_domain",
            "effective_domain",
            "ready_to_apply",
            "Each canonical adjacent cutoff lies in one common high-height box.",
            exact["domain"],
            "The statement concerns the retained first-order edge model.",
        ),
        GateRow(
            "acres_02_alignment",
            "projective_normalization",
            "ready_to_apply",
            "The alternating endpoint normalizer is harmless in real projective space.",
            exact["physical_alignment"],
            "The aligned normalized rows, not the antipodal physical representatives, define the short join.",
        ),
        GateRow(
            "acres_03_endpoint_values",
            "exact_source_identity",
            "certified",
            "The two cutoff endpoint values agree and their first slopes have an exact gap.",
            "; ".join(f"{key}={value}" for key, value in symbolic["endpoint_values"].items()),
            "These are exact evaluations of the published C_0 source formulas.",
        ),
        GateRow(
            "acres_04_leading_wedge",
            "exact_orientation",
            "certified",
            "The limiting cutoff connector is strictly clockwise.",
            symbolic["leading_wedge"],
            "No numerical sign fitting is used.",
        ),
        GateRow(
            "acres_05_finite_jets",
            "finite_height_transfer",
            "certified",
            "Only the value and centered-slope edge jets are needed at a cutoff connector.",
            exact["finite_jets"],
            "The x-second jet is needed for within-cell current, not for the join-parameter current.",
        ),
        GateRow(
            "acres_06_wedge_expansion",
            "division_free_perturbation",
            "certified",
            "The finite connector determinant has an exact polynomial remainder.",
            symbolic["finite_wedge"],
            "No division by either endpoint row occurs.",
        ),
        GateRow(
            "acres_07_error_budget",
            "rational_majorant",
            "certified",
            "The O(h^2) connector error cannot consume the radical margin.",
            "500h+50000h^2<1/600 and sqrt(2)/32>1/24",
            "All finite-height constants are inherited from the critical-ray gate.",
            diagnostics=arithmetic,
        ),
        GateRow(
            "acres_08_signed_wedge",
            "effective_sign_theorem",
            "certified",
            "Every retained-model adjacent edge connector has strict negative projective current.",
            arithmetic["signed_wedge"],
            "This signs the connector, not the omitted Xi/source remainder.",
        ),
        GateRow(
            "acres_09_nonvanishing",
            "zero_fibre_guard",
            "certified",
            "The affine connector remains in one nonzero projective chart.",
            arithmetic["nonvanishing"],
            "The join cannot pass through the origin.",
        ),
        GateRow(
            "acres_10_affine_current",
            "exact_homotopy",
            "certified",
            "The affine join has constant signed projective current.",
            symbolic["affine_join"] + "; " + exact["orientation"],
            "The join parameter s is distinct from the physical x derivative.",
        ),
        GateRow(
            "acres_11_q_ge_1",
            "domain_corollary",
            "ready_to_apply",
            "The result closes every q>=1 cutoff splice in the retained critical-ray model.",
            "q=2tL^2>=1 is a subset of L>=50, 0<tL<=25",
            "Insertion into the cumulative contact scalar remains separate.",
        ),
        GateRow(
            "acres_12_boundary",
            "nonpromotion_guard",
            "guard_validated",
            "The model splice is not an Xi-level or global contact theorem.",
            "Higher-order Xi/source error and the cumulative signed-minor estimate remain open.",
            "No Xi-level edge sign, Abel gap, winding cap, contact exclusion, Lambda<=0, or RH is claimed.",
        ),
    ]


def render_note(payload: dict) -> str:
    symbolic = payload["symbolic_certificate"]
    arithmetic = payload["arithmetic_certificate"]
    exact = payload["exact"]
    values = symbolic["endpoint_values"]
    return f"""# Adjacent-Cutoff Real-Edge Projective Splice Gate

Date: 2026-07-31

Status: effective signed adjacent-cutoff splice for the retained first-order
real-edge model. This is not a proof of an Xi-level edge theorem, an Abel
gap, `Lambda <= 0`, or RH. It is not an Xi-level result.

```text
work/rh_compute/results/{STEM}.json
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```

## Cutoff Charts

{exact["domain"]}.

The adjacent recurrence gives `kappa_m=-kappa_(m-1)`. Therefore

```text
{exact["physical_alignment"]}
```

Projective points are unchanged by either real nonzero normalizer. The
proved bound `c_minus,c_plus<-1/4` aligns the two representatives in the
same half-plane and selects the short projective connector.

## Exact Endpoint Wedge

At the old endpoint `p_minus=-1` and new endpoint `p_plus=1`,

```text
A_minus={values["A_minus"]}
A_plus ={values["A_plus"]}
B_minus={values["B_minus"]}
B_plus ={values["B_plus"]}

{symbolic["leading_wedge"]}.
```

These radicals come from exact evaluation of the completed-zeta and
Riemann-Siegel `C_0` phases at `p=+-1`. No fitted circle or polygon defines
`pi` or the radicals.

## Finite-Height Budget

The critical-ray theorem supplies, independently at both endpoints,

```text
{exact["finite_jets"]}.
```

Direct expansion gives

```text
{symbolic["finite_wedge"]}.
```

Using `|A|<1/2`, `|B_minus|,|B_plus|<1`, and the saved defect caps,
the two remainder coefficients are at most `500` and `50000`. Since

```text
h<{arithmetic["h_upper"]},
500h+50000h^2<{arithmetic["error_cap"]},
sqrt(2)/32>1/24,
```

the strict cutoff sign is

```text
{arithmetic["normalized_wedge"]},
{arithmetic["signed_wedge"]}.
```

## Signed Join

For the affine connector

```text
{symbolic["affine_join"]},
{arithmetic["nonvanishing"]}.
```

Hence

```text
{exact["orientation"]}.
```

The apparent missing `rho_xx` or second adjacent x-jet is not needed for
this connector: the within-cell x-current uses `d_(edge,x)`, whereas the
cutoff homotopy current differentiates `(c_s,d_s)` with respect to `s` and
uses only the two certified row values.

## Consequence and Boundary

Together with the full critical-ray within-cell sign, this closes the
pointwise cells and every adjacent cutoff by strict clockwise projective
arcs for the retained first-order edge model, including every `q>=1` point
in the stated critical range.

The omitted higher-order Xi/source remainder has not been transferred to
this projective connector, and the edge block has not been inserted into
the complete cumulative contact/minor scalar. No Xi-level edge sign,
cumulative-minor estimate, Abel-scalar gap, successor winding cap, contact
exclusion, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion follows.
"""


def build_artifact() -> dict:
    exact = exact_payload()
    symbolic = symbolic_certificate()
    arithmetic = arithmetic_certificate()
    rows = build_rows(exact, symbolic, arithmetic)
    return {
        "kind": STEM,
        "date": "2026-07-31",
        "status": (
            "effective finite-height signed adjacent-cutoff projective "
            "splice for the retained first-order real-edge model"
        ),
        "proof_boundary": (
            "This proves the adjacent-cutoff signed splice only for the "
            "retained first-order real-edge model. It does not bound the "
            "higher-order Xi/source contribution, prove an Xi-level edge "
            "sign, cumulative-minor estimate, Abel-scalar gap, successor "
            "winding cap, contact exclusion, Lambda<=0, PF-infinity, RH, "
            "or a prize-level conclusion."
        ),
        "source_audit": source_audit(),
        "exact": exact,
        "symbolic_certificate": symbolic,
        "arithmetic_certificate": arithmetic,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "exact_cutoff_endpoint_values": 4,
            "exact_leading_wedges": 1,
            "finite_connector_jet_families": 2,
            "rational_remainder_budgets": 1,
            "nonvanishing_affine_joins": 1,
            "signed_affine_splices": 1,
            "uniform_q_ge_1_model_splices": 1,
            "required_second_adjacent_x_jets": 0,
            "xi_level_edge_splices": 0,
            "abel_gaps": 0,
            "winding_bounds": 0,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_artifact()
    atomic_write(args.output, json.dumps(payload, indent=2) + "\n")
    atomic_write(args.note, render_note(payload))
    print(
        "built adjacent-cutoff real-edge projective splice gate: "
        "12 rows, 4 endpoint values, 1 exact wedge, "
        "1 nonvanishing join, 1 signed q>=1 model splice"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
