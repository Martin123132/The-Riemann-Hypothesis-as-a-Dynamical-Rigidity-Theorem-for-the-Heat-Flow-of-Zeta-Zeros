#!/usr/bin/env python3
"""Build the edge-anchored near-terminal pair-current guard."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402
import jensen_window_pf_newman_polymath15_first_order_centered_real_edge_projective_current_gate as edge_gate  # noqa: E402


STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "edge_anchored_near_terminal_pair_current_guard"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
PRECISION_BITS = 192
SOURCE_PATHS = {
    "edge_symbol": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "real_edge_projective_current_gate.json"
    ),
    "finite_edge": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "critical_ray_finite_height_real_edge_gate.json"
    ),
    "near_terminal": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "endpoint_relative_phase_current_recurrence_gate.json"
    ),
    "interior_current": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "interior_projective_current_gate.json"
    ),
    "cumulative_handoff": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "residual_placement_c1_cumulative_handoff_gate.json"
    ),
    "edge_builder": (
        REPO_ROOT
        / "work/rh_compute/scripts/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "real_edge_projective_current_gate.py"
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


def load_sources() -> dict[str, dict | str]:
    missing = [str(path) for path in SOURCE_PATHS.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    payloads: dict[str, dict | str] = {}
    for key, path in SOURCE_PATHS.items():
        payloads[key] = (
            json.loads(path.read_text(encoding="utf-8"))
            if path.suffix == ".json"
            else path.read_text(encoding="utf-8")
        )
    return payloads


def source_audit(payloads: dict[str, dict | str]) -> dict:
    edge = payloads["edge_symbol"]
    assert isinstance(edge, dict)
    if "K_edge(p)=[A*A''-(A')^2]/(16pi^2)" not in edge.get(
        "exact", {}
    ).get("asymptotic_jets", {}).get("current_limit", ""):
        raise RuntimeError("edge symbol drifted")

    near = payloads["near_terminal"]
    assert isinstance(near, dict)
    if "a*u_(N-m)->m+theta" not in near.get("exact", {}).get(
        "near_terminal_asymptotic", {}
    ).get("source_limits", ""):
        raise RuntimeError("near-terminal source limits drifted")

    interior = payloads["interior_current"]
    assert isinstance(interior, dict)
    if "J_n<=-(5/64)" not in interior.get("exact", {}).get(
        "negative_definite_theorem", ""
    ):
        raise RuntimeError("interior current theorem drifted")

    handoff = payloads["cumulative_handoff"]
    assert isinstance(handoff, dict)
    if "K(V_J)=-R_diag+C_cross" not in handoff.get("exact", {}).get(
        "diagonal_reserve", ""
    ):
        raise RuntimeError("cumulative handoff drifted")

    builder = payloads["edge_builder"]
    assert isinstance(builder, str)
    for marker in ("def direct_a_jet", "current_symbol_margin"):
        if marker not in builder:
            raise RuntimeError(f"edge builder marker missing: {marker}")

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
            "location": "equation (53), the critical saddle, and the first correction",
            "url": "https://arxiv.org/abs/1904.12438",
        },
    }


def symbolic_certificate() -> dict:
    h = sp.symbols("h", positive=True)
    theta = sp.symbols("theta", real=True)
    m = sp.symbols("m", integer=True, positive=True)
    log_ratio = sp.log(1 - h * (theta + m)) - sp.log(1 - h * theta)
    reduced_phase = (
        2 * sp.pi * log_ratio / h**2
        + 2 * sp.pi * m * (1 / h - theta)
    )
    phase_limit = sp.simplify(sp.limit(reduced_phase, h, 0, dir="+"))
    expected_phase = -sp.pi * (4 * m * theta + m**2)
    if sp.simplify(phase_limit - expected_phase) != 0:
        raise RuntimeError("near-terminal phase limit failed")

    a, b, mm, x, y, k = sp.symbols("A B M X Y k", real=True)
    d = k * y / 2
    n = y / (16 * sp.pi) - k**2 * x / 4
    edge_current = a * mm - b**2
    interior_current = x * n - d**2
    cross_current = a * n + x * mm - 2 * b * d
    pair_current = (a + x) * (mm + n) - (b + d) ** 2
    if sp.simplify(
        pair_current - edge_current - interior_current - cross_current
    ) != 0:
        raise RuntimeError("edge-carrier polarization failed")

    return {
        "phase_limit": (
            "z_(N-m)/z_N -> R_m(theta)=(-1)^m*exp(-4*pi*i*m*theta)"
        ),
        "reduced_phase_limit": (
            "2*pi*h^-2[log(1-h(theta+m))-log(1-h theta)]"
            "+2*pi*m(h^-1-theta) -> -pi(4m theta+m^2)"
        ),
        "carrier_source": (
            "Z_m=-Q(p)R_m(theta)=X_m+iY_m, "
            "theta=(1-p)/2, k=m+theta"
        ),
        "carrier_jets": (
            "c_m/S_a->X_m; a*d_m/S_a->D_m=kY_m/2; "
            "a(c_(m,x)-r_Sc_m)/S_a->D_m; "
            "a^2(d_(m,x)-r_Sd_m)/S_a->"
            "N_m=Y_m/(16*pi)-k^2X_m/4"
        ),
        "edge_jets": (
            "c_e/S_a->A; a*d_e/S_a->B; "
            "a(c_(e,x)-r_Sc_e)/S_a->B; "
            "a^2(d_(e,x)-r_Sd_e)/S_a->M"
        ),
        "cross_symbol": "C_(e,m)=A*N_m+X_m*M-2*B*D_m",
        "pair_symbol": (
            "P_m(p)=(A+X_m)(M+N_m)-(B+D_m)^2="
            "K_edge+K_m+C_(e,m)"
        ),
    }


def evaluate_pair(point: Fraction, m: int) -> dict:
    p = edge_gate.aq(point)
    pi = flint.arb.pi()
    theta = (1 - p) / 2
    q_phase = p * p / 2 - p + edge_gate.aq(Fraction(3, 8))
    q = flint.acb(0, -pi * q_phase).exp()
    relative = (-1) ** m * flint.acb(0, -4 * pi * m * theta).exp()
    z = -q * relative

    a, a_p, a_pp = edge_gate.direct_a_jet(p)
    b = -a_p / (4 * pi)
    mm = a_pp / (16 * pi**2)
    k = m + theta
    x, y = z.real, z.imag
    d = k * y / 2
    n = y / (16 * pi) - k**2 * x / 4
    edge_current = a * mm - b**2
    interior_current = x * n - d**2
    cross_current = a * n + x * mm - 2 * b * d
    pair_current = (a + x) * (mm + n) - (b + d) ** 2

    return {
        "p": f"{point.numerator}/{point.denominator}",
        "m": m,
        "theta": str((1 - point) / 2),
        "edge_current_ball": edge_current.str(70, more=True),
        "interior_current_ball": interior_current.str(70, more=True),
        "cross_current_ball": cross_current.str(70, more=True),
        "pair_current_ball": pair_current.str(70, more=True),
        "pair_lower": pair_current.lower().str(70),
        "pair_upper": pair_current.upper().str(70),
    }


def interval_certificate() -> dict:
    flint.ctx.prec = PRECISION_BITS
    positive = evaluate_pair(Fraction(-5, 8), 3)
    negative = evaluate_pair(Fraction(0), 3)
    if not flint.arb(positive["pair_lower"]) > edge_gate.aq(Fraction(1, 10)):
        raise RuntimeError("positive pair-current witness failed")
    if not flint.arb(negative["pair_upper"]) < edge_gate.aq(Fraction(-2)):
        raise RuntimeError("negative pair-current witness failed")
    return {
        "precision_bits": PRECISION_BITS,
        "positive_margin": "P_3(-5/8)>1/10",
        "negative_margin": "P_3(0)<-2",
        "positive_witness": positive,
        "negative_witness": negative,
    }


def exact_payload() -> dict:
    return {
        "path": (
            "Fix p in [-1,1], theta=(1-p)/2, a=N+theta, "
            "q=2tL^2=1, and n=N-m with fixed integer m>=1; let N->infinity."
        ),
        "q": (
            "Q(p)=exp[-pi*i(p^2/2-p+3/8)] and "
            "S_a=kappa_N*T_0"
        ),
        "sign_reversal": (
            "For m=3, P_3(-5/8)>1/10 while P_3(0)<-2. "
            "By convergence, the retained edge plus the N-3 carrier has "
            "both projective-current orientations on sufficiently high "
            "physical q=1 cells."
        ),
        "route_guard": (
            "No proof may absorb the interior atoms into the signed edge one "
            "at a time by asserting K(v_edge+v_n)<0 or "
            "C_(e,n)<-K_edge-K_n uniformly."
        ),
        "surviving_target": (
            "The complete Xi-specific signed cross aggregate C_cross<R_diag, "
            "or the endpoint-complete division-free Abel gap, remains open. "
            "Other carriers may provide essential cancellation."
        ),
        "pi_provenance": (
            "The pi in Q, the saddle phase, and the near-terminal ratio is "
            "inherited from the completed-zeta and Riemann-Siegel "
            "normalization. It is not fitted from a circle, polygon, or plot."
        ),
    }


def build_rows(
    exact: dict, symbolic: dict, interval: dict
) -> list[GateRow]:
    return [
        GateRow(
            "eant_01_path",
            "cofinal path",
            "proved",
            "The near-terminal test stays on the physical q=1 boundary.",
            exact["path"],
            "The test concerns fixed m and fixed p as N tends to infinity.",
        ),
        GateRow(
            "eant_02_phase_ratio",
            "exact asymptotic",
            "proved",
            "The N-m carrier has an explicit phase relative to the terminal carrier.",
            symbolic["phase_limit"],
            "The integer phase exp(-2*pi*i*m*N) is exactly one.",
        ),
        GateRow(
            "eant_03_phase_derivation",
            "symbolic audit",
            "proved",
            "The phase ratio follows from the saddle logarithm before projection.",
            symbolic["reduced_phase_limit"],
            "Heat and first-correction terms vanish at this leading scale.",
        ),
        GateRow(
            "eant_04_carrier_source",
            "leading source",
            "proved",
            "The limiting carrier is a unit complex source with explicit real coordinates.",
            symbolic["carrier_source"],
            "No argument branch is selected.",
        ),
        GateRow(
            "eant_05_carrier_jets",
            "four-jet asymptotic",
            "proved",
            "All four radial-subtracted carrier jets have closed leading symbols.",
            symbolic["carrier_jets"],
            "The Y/(16*pi) term from u_x is retained.",
        ),
        GateRow(
            "eant_06_edge_jets",
            "inherited theorem",
            "proved",
            "The retained edge supplies its four already-certified leading jets.",
            symbolic["edge_jets"],
            "This does not alter the finite-height edge sign theorem.",
        ),
        GateRow(
            "eant_07_polarization",
            "exact algebra",
            "proved",
            "The edge-carrier cross and pair symbols retain every mixed term.",
            symbolic["cross_symbol"] + "; " + symbolic["pair_symbol"],
            "No diagonal current is counted twice or discarded.",
        ),
        GateRow(
            "eant_08_positive",
            "interval witness",
            "certified",
            "The edge plus N-3 carrier is counterclockwise at p=-5/8 in the leading model.",
            interval["positive_margin"],
            "This is a 192-bit Arb point certificate.",
            diagnostics=interval["positive_witness"],
        ),
        GateRow(
            "eant_09_negative",
            "interval witness",
            "certified",
            "The same pair is clockwise at p=0 in the leading model.",
            interval["negative_margin"],
            "This is a 192-bit Arb point certificate.",
            diagnostics=interval["negative_witness"],
        ),
        GateRow(
            "eant_10_reversal",
            "asymptotic guard",
            "proved",
            "The physical near-terminal pair has both signs for sufficiently high q=1 cells.",
            exact["sign_reversal"],
            "No finite threshold is needed to reject a uniform theorem.",
        ),
        GateRow(
            "eant_11_route",
            "route decision",
            "proved",
            "Termwise edge absorption is rejected; complete aggregate cancellation remains viable.",
            exact["route_guard"] + " " + exact["surviving_target"],
            "This does not falsify the complete Xi cross-current or Abel theorem.",
        ),
        GateRow(
            "eant_12_boundary",
            "proof boundary",
            "proved",
            "The guard narrows the route without claiming contact exclusion.",
            exact["pi_provenance"],
            "This artifact does not prove an aggregate cross bound, Abel gap, winding cap, contact exclusion, Lambda<=0, PF-infinity, RH, or a prize conclusion.",
        ),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    symbolic = payload["symbolic_certificate"]
    interval = payload["interval_certificate"]
    return f"""# Edge-Anchored Near-Terminal Pair-Current Guard

Date: 2026-07-31

Status: asymptotic one-carrier sign-reversal guard. This is not a proof of a
complete aggregate bound, Abel gap, `Lambda<=0`, or RH.

## Path

```text
{exact['path']}
```

The exact saddle expansion gives

```text
{symbolic['reduced_phase_limit']}
{symbolic['phase_limit']}
```

Hence

```text
{symbolic['carrier_source']}
```

## Four Jets

The near-terminal carrier has

```text
{symbolic['carrier_jets']}
```

and the retained endpoint-terminal edge has

```text
{symbolic['edge_jets']}
```

## Pair Current

Exact polarization gives

```text
{symbolic['cross_symbol']}
{symbolic['pair_symbol']}
```

The 192-bit Arb certificates are

```text
{interval['positive_margin']}
{interval['negative_margin']}
```

Thus the same edge plus `N-3` carrier has both current orientations in the
cofinal `q=1` leading model.  Convergence of all four displayed jets makes
the reversal physical for sufficiently high cells.

## Route Decision

```text
{exact['route_guard']}
```

This rejects termwise edge absorption.  It does not reject the complete
aggregate because the remaining carriers may supply essential signed
cancellation.

```text
{exact['surviving_target']}
```

{exact['pi_provenance']}

This artifact does not prove an aggregate cross-current bound, Abel gap,
winding cap, contact exclusion, `Lambda<=0`, PF-infinity, RH, or a prize-level
conclusion.
"""


def build_payload() -> dict:
    payloads = load_sources()
    symbolic = symbolic_certificate()
    interval = interval_certificate()
    exact = exact_payload()
    rows = build_rows(exact, symbolic, interval)
    return {
        "kind": STEM,
        "date": "2026-07-31",
        "status": (
            "asymptotic edge-anchored near-terminal pair-current sign-reversal "
            "guard with two Arb witnesses"
        ),
        "proof_boundary": (
            "This gate proves the fixed-m q=1 carrier phase and four-jet "
            "limits, the exact edge-carrier polarization, and opposite signs "
            "for the edge plus N-3 leading pair. It rejects uniform termwise "
            "edge absorption only. It proves no complete Xi aggregate "
            "cross-current bound, Abel gap, winding cap, contact exclusion, "
            "Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
        "source_audit": source_audit(payloads),
        "exact": exact,
        "symbolic_certificate": symbolic,
        "interval_certificate": interval,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "near_terminal_phase_limits": 1,
            "carrier_jet_limits": 4,
            "edge_jet_limits": 4,
            "exact_pair_polarizations": 1,
            "arb_point_witnesses": 2,
            "positive_pair_symbols": 1,
            "negative_pair_symbols": 1,
            "uniform_termwise_edge_absorptions": 0,
            "complete_aggregate_counterexamples": 0,
            "abel_gaps": 0,
            "contact_exclusions": 0,
        },
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    atomic_write(args.out, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    atomic_write(args.note, render_note(payload))
    print(
        "built edge-anchored near-terminal pair-current guard: "
        "12 rows, 4 carrier jets, 2 Arb witnesses, "
        "1 asymptotic sign reversal"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
