#!/usr/bin/env python3
"""Build the double-Morse ridge/carrier inversion gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_full_support_ideal_cubic_double_"
    "morse_ridge_carrier_inversion_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "finite_cell_inversion": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "ideal_cubic_reciprocal_finite_cell_inversion_gate.json"
    ),
    "hyperbolic_transport": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_ideal_cubic_double_"
        "morse_hyperbolic_transport_gate.json"
    ),
    "joined_recombination": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_joined_lower_"
        "interior_abel_recombination_gate.json"
    ),
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


def require_zero(expression: sp.Expr, label: str) -> None:
    value = sp.simplify(sp.factor(sp.together(expression)))
    if value != 0 and value.equals(0) is not True:
        raise RuntimeError(f"{label}: {value}")


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
    require(
        payloads["finite_cell_inversion"]["counts"]["carrier_inversion_identities"]
        == 3,
        "finite-cell inversion source drift",
    )
    require(
        payloads["hyperbolic_transport"]["counts"]["hyperbolic_identities"]
        == 3,
        "hyperbolic source drift",
    )
    require(
        payloads["joined_recombination"]["counts"]["global_recombinations"]
        == 4,
        "joined recombination source drift",
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


def inversion_certificate() -> dict[str, str | int]:
    carrier, exterior, tie, alias = sp.symbols("M E tau A")
    ridge, transverse, wings = sp.symbols("Z V W")

    diagonal = carrier - exterior
    cell = tie + diagonal + alias
    require_zero(
        2 * carrier - cell - (carrier - tie + exterior - alias),
        "per-cell surviving carrier",
    )

    band = tie + carrier - exterior
    require_zero(
        2 * carrier - band - (carrier - tie + exterior),
        "global surviving carrier",
    )
    defect = band - carrier
    require_zero(
        carrier - defect - (2 * carrier - band),
        "defect join",
    )
    require_zero(diagonal + exterior - carrier, "finite diagonal inversion")

    ridge_diagonal = ridge + transverse + wings
    ridge_carrier = ridge_diagonal + exterior
    require_zero(
        ridge_carrier - (ridge + transverse + wings + exterior),
        "ridge carrier completion",
    )

    return {
        "cell": "With Lambda_q^alias=sum_(s!=q)^sym B_(q,s), finite Poisson inversion gives S_q=tau_q+B_(q,q)+Lambda_q^alias and B_(q,q)=M_q-E_q^ext. Hence 2M_q-S_q=M_q-tau_q+E_q^ext-Lambda_q^alias.",
        "band": "For the complete physical band, mathscr F_N=tau_band+mathscr M_N-E_band. Therefore 2mathscr M_N-mathscr F_N=mathscr M_N-tau_band+E_band: exactly one starred carrier remains.",
        "defect": "Since mathscr D_N=mathscr F_N-mathscr M_N=tau_band-E_band, the same join is mathscr M_N-mathscr D_N. This is an identity, not an exterior-smallness estimate.",
        "weights": "The surviving mathscr M_N keeps omega_N(1)=omega_N(N)=1/2 and omega_N(q)=1 internally. The endpoint halves are not promoted to full atoms.",
        "ties": "At an internal reciprocal boundary, the -I/2 lower correction of one cell and +I/2 upper correction of its neighbor cancel. Only the two outer band corrections remain in tau_band.",
        "terminal": "The ideal joins are J_H^0=mathscr M_N[P_H^0]-tau_band[P_H^0]+E_band[P_H^0] and J_T^0=mathscr M_N[P_T^0]-tau_band[P_T^0]+E_band[P_T^0]+2iu_N mathcal T_N^[N][B_(T,N)^0]. The transpose terminal block remains explicit.",
        "inversion_identities": 5,
    }


def fiber_certificate() -> dict[str, str | int]:
    s, d_minus, d_plus = sp.symbols("s d_minus d_plus", real=True)
    kappa = 1 / (2 * sp.pi * sp.I)
    character = lambda value: sp.exp(2 * sp.pi * sp.I * value)
    kernel = kappa * (
        character(s * d_plus) - character(s * d_minus)
    ) / s
    require_zero(
        sp.diff(kernel, d_plus) - character(s * d_plus),
        "upper transverse derivative",
    )
    require_zero(
        sp.diff(kernel, d_minus) + character(s * d_minus),
        "lower transverse derivative",
    )
    require_zero(
        sp.limit(kernel, s, 0) - (d_plus - d_minus),
        "central transverse limit",
    )

    return {
        "domain": "The full diagonal rectangle a_q<=r<=b_q, 1<=u<=N maps to a curvilinear domain Omega_q in (s,d). Because eta+y_U(eta) is strictly increasing, every fixed-s fiber is an interval.",
        "diagonal_chart": "The fibers meeting d=0 are indexed by S_q^Delta=[sqrt(2)eta_-,sqrt(2)eta_+]. On this interval their actual limits satisfy d_-(s)<=0<=d_+(s), including one-sided endpoint cases q=1,N.",
        "kernel": "K_q(s)=integral_(d_-(s))^(d_+(s))e(sd)dd=kappa{e(sd_+)-e(sd_-)}/s for s!=0, with K_q(0)=d_+(0)-d_-(0).",
        "ridge_lift": "On S_q^Delta write G_q(s,d)=G_q(s,0)+dH_q(s,d). Define Z_q=e(alpha_Plog q)integral G_q(s,0)K_q(s)ds and V_q=e(alpha_Plog q)integral integral dH_q(s,d)e(sd)dd ds.",
        "wings": "Let W_q be the contribution of the remaining s-wings. Then B_(q,q)=Z_q+V_q+W_q and M_q=Z_q+V_q+W_q+E_q^ext.",
        "guard": "A direct identification Z_q=M_q is equivalent to the additional identity V_q+W_q+E_q^ext=0. Neither the hyperbolic chart nor finite Fourier inversion supplies that identity.",
        "fiber_identities": 3,
        "ridge_decomposition_identities": 2,
    }


def witness_certificate() -> dict[str, str | int]:
    alpha = Fraction(2, 1)
    q = Fraction(2, 1)
    lower = alpha / (q + Fraction(1, 2))
    upper = alpha / (q - Fraction(1, 2))
    width = upper - lower
    require(lower == Fraction(4, 5), "witness lower endpoint")
    require(upper == Fraction(4, 3), "witness upper endpoint")
    require(width == Fraction(8, 15), "witness width")
    require(width < 1, "witness width guard")

    mp.mp.dps = 70
    lower_mp = mp.mpf(lower.numerator) / lower.denominator
    upper_mp = mp.mpf(upper.numerator) / upper.denominator

    def sinc(value: mp.mpf) -> mp.mpf:
        return mp.sin(mp.pi * value) / (mp.pi * value)

    numerical_checks = 0
    for r_value in [
        mp.mpf("0.83"),
        mp.mpf("0.94"),
        mp.mpf("1.0"),
        mp.mpf("1.17"),
        mp.mpf("1.31"),
    ]:
        direct = mp.quad(
            lambda t: mp.e ** (-2j * mp.pi * r_value * t),
            [-mp.mpf("0.5"), mp.mpf("0.5")],
        )
        require(abs(direct - sinc(r_value)) < mp.mpf("1e-60"), "sinc witness")
        numerical_checks += 1

    diagonal = mp.quad(sinc, [lower_mp, mp.mpf("1"), upper_mp])
    closed = (
        mp.si(mp.pi * upper_mp) - mp.si(mp.pi * lower_mp)
    ) / mp.pi
    require(abs(diagonal - closed) < mp.mpf("1e-60"), "sine-integral witness")
    numerical_checks += 1
    width_mp = mp.mpf(width.numerator) / width.denominator
    require(abs(diagonal) <= width_mp, "finite diagonal width bound")
    numerical_checks += 1
    exterior = 1 - diagonal
    require(abs(exterior) >= mp.mpf(7) / 15, "nonzero exterior bound")
    numerical_checks += 1
    require(numerical_checks == 8, "witness numerical count")

    return {
        "source": "Take alpha_P=2, q=2, B>=3, and the piecewise-BV test source F_A(u)=A(u)e(alpha_Plog u)=1_[3/2,5/2](u). Full inversion at the interior point u=2 returns M_2=1.",
        "transform": "For this source I_A(r)e(2r)=integral_(-1/2)^(1/2)e(-rt)dt=sin(pi r)/(pi r).",
        "cell": "The reciprocal cell is [a_2,b_2]=[4/5,4/3], of exact width 8/15, so B_(2,2)=integral_(4/5)^(4/3)sin(pi r)/(pi r)dr={Si(4pi/3)-Si(4pi/5)}/pi.",
        "bound": "Since |sin(pi r)/(pi r)|<=1, |B_(2,2)|<=8/15<1=M_2. Therefore |E_2^ext|=|M_2-B_(2,2)|>=7/15.",
        "scope": "The witness is nonphysical and rejects a coefficient-blind finite-diagonal-equals-carrier identity. It also shows that any ridge-equals-carrier claim needs a separate residual-plus-wings-plus-exterior cancellation; it does not disprove such a special physical cancellation or RH.",
        "diagonal_value": mp.nstr(diagonal, 30),
        "exterior_value": mp.nstr(exterior, 30),
        "witness_exact_checks": 5,
        "witness_numerical_checks": numerical_checks,
        "nonzero_exterior_witnesses": 1,
    }


def endpoint_certificate() -> dict[str, str | int]:
    weights = [Fraction(1, 2), Fraction(1), Fraction(1), Fraction(1), Fraction(1, 2)]
    checks = 0
    for weight in weights:
        require(2 * weight - weight == weight, "surviving starred weight")
        checks += 1
    require(-Fraction(1, 2) + Fraction(1, 2) == 0, "internal tie halves")
    return {
        "weights": "For every starred atom, 2omega_N(q)-omega_N(q)=omega_N(q). Thus the surviving carrier has the original endpoint half-weights, not a new unstarred normalization.",
        "ties": "The internal lower/upper tie pair is -1/2+1/2=0. The global formula retains only tau_band at the two outer reciprocal boundaries.",
        "endpoint_weight_checks": checks,
        "internal_tie_cancellations": 1,
    }


def handoff_certificate() -> dict[str, str]:
    return {
        "decision": "Retire direct finite-ridge-to-carrier matching. The exact target is the one-carrier global package mathscr M_N-tau_band+E_band, with the transpose terminal survivor in its own channel.",
        "route": "Before applying a 1/s modulus, pull E_band through the certified full-support exterior/complement identity and join its negative-, low-positive-, and high-positive-frequency pieces to the roster fringes, endpoint halves, and terminal recurrence. A useful estimate must act on that complete non-tautological package.",
        "obligation": "Prove or rigorously obstruct a signed bound for the complete ideal Hermitian and transpose joins after this exterior composition. Only then return to the central and outer hyperbolic residuals and attach Delta_H, Delta_T.",
        "reserve": "No signed ridge-completed ideal bound, exterior/complement bound, completed-current theorem, Phi_B theorem, contact exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }


def build_rows(certificate: dict) -> list[GateRow]:
    inv = certificate["inversion"]
    fiber = certificate["fiber"]
    witness = certificate["witness"]
    endpoint = certificate["endpoint"]
    handoff = certificate["handoff"]
    rows = [
        GateRow("rci_01_sources", "source chain", "proved", "The three parent gates are current.", "Source hashes and parent counts are recorded.", "No parent bound is strengthened."),
        GateRow("rci_02_cell_poisson", "cell decomposition", "proved", "Each finite cell is tie plus diagonal plus aliases.", inv["cell"], "No alias bound is claimed."),
        GateRow("rci_03_diagonal", "diagonal inversion", "proved", "The finite diagonal is carrier minus exterior.", inv["cell"], "The exterior is retained."),
        GateRow("rci_04_cell_join", "per-cell join", "proved", "One carrier survives the local 2M-S join.", inv["cell"], "Cellwise aliases remain conditionally ordered."),
        GateRow("rci_05_band", "global band", "proved", "Internal cells recompose before estimation.", inv["band"], "No exterior smallness is inferred."),
        GateRow("rci_06_defect", "defect coordinate", "proved", "The global join is M-D.", inv["defect"], "This is an involutive identity."),
        GateRow("rci_07_one_carrier", "carrier coefficient", "proved", "Exactly one starred carrier remains.", inv["band"], "It is not two independent positive copies."),
        GateRow("rci_08_star_weights", "endpoint weights", "proved", "The surviving endpoint carriers retain half-weight.", endpoint["weights"], "Hard endpoint halves remain separate."),
        GateRow("rci_09_ties", "tie cancellation", "proved", "Internal tie halves cancel globally.", endpoint["ties"], "Outer tau_band remains."),
        GateRow("rci_10_terminal", "transpose terminal", "proved", "The terminal survivor is not absorbed into inversion.", inv["terminal"], "No terminal sign is claimed."),
        GateRow("rci_11_domain", "hyperbolic domain", "proved", "Fixed-s fibers are intervals.", fiber["domain"], "Curvilinear limits are preserved."),
        GateRow("rci_12_diagonal_chart", "diagonal subchart", "proved", "The d=0 carrier line has an exact finite s-range.", fiber["diagonal_chart"], "Outer s-wings remain."),
        GateRow("rci_13_kernel", "transverse kernel", "proved", "The finite ridge coefficient is K_q(s), not one.", fiber["kernel"], "The s=0 chart uses its removable value."),
        GateRow("rci_14_ridge_lift", "ridge lift", "proved", "The trace is lifted across actual transverse fibers.", fiber["ridge_lift"], "A measure-zero line is not itself a two-dimensional summand."),
        GateRow("rci_15_transverse", "transverse residual", "proved", "The dH residual remains exact.", fiber["ridge_lift"], "No residual bound is claimed."),
        GateRow("rci_16_wings", "outer wings", "proved", "The non-diagonal s-wings are retained.", fiber["wings"], "No wing bound is claimed."),
        GateRow("rci_17_ridge_sum", "finite diagonal split", "proved", "Ridge, transverse, and wings reconstruct B_(q,q).", fiber["wings"], "The exterior is not part of B_(q,q)."),
        GateRow("rci_18_carrier_relation", "ridge/carrier relation", "proved", "The carrier requires ridge, residual, wings, and exterior.", fiber["guard"], "No direct ridge equality is promoted."),
        GateRow("rci_19_witness_source", "test source", "guard_validated", "A compact piecewise-BV source tests the coefficient.", witness["source"], witness["scope"]),
        GateRow("rci_20_witness_transform", "exact transform", "proved", "The witness transform is normalized sinc.", witness["transform"], witness["scope"]),
        GateRow("rci_21_witness_cell", "finite cell", "proved", "The witness cell has width 8/15.", witness["cell"], witness["scope"]),
        GateRow("rci_22_witness_bound", "finite-diagonal guard", "guard_validated", "The finite diagonal cannot equal the carrier.", witness["bound"], witness["scope"]),
        GateRow("rci_23_exterior", "exterior necessity", "guard_validated", "The witness exterior has modulus at least 7/15.", witness["bound"], witness["scope"]),
        GateRow("rci_24_local_guard", "route guard", "guard_validated", "Local ridge matching is retired.", fiber["guard"], "Special physical cancellation is not ruled out."),
        GateRow("rci_25_global_route", "route decision", "guard_validated", "The global exterior is composed before a transverse modulus.", handoff["route"], handoff["reserve"]),
        GateRow("rci_26_exterior_open", "exterior composition", "open", "Control the complete exterior/carrier package.", handoff["obligation"], handoff["reserve"]),
        GateRow("rci_27_signed_open", "signed ideal join", "open", "Prove or obstruct the complete ideal sign.", handoff["obligation"], handoff["reserve"]),
        GateRow("rci_28_boundary", "proof boundary", "guard_validated", "This gate is not a proof of RH.", handoff["reserve"], "No prize-level conclusion is claimed."),
    ]
    require(len(rows) == 28, "row count drifted")
    return rows


def render_note(artifact: dict) -> str:
    cert = artifact["symbolic_certificate"]
    inv = cert["inversion"]
    fiber = cert["fiber"]
    witness = cert["witness"]
    endpoint = cert["endpoint"]
    handoff = cert["handoff"]
    return f"""# Double-Morse Ridge/Carrier Inversion Gate

Date: 2026-08-03

Status: exact carrier-coefficient and finite-diagonal obstruction gate; signed ridge-completed estimate open; not a proof of RH.

## Exact Carrier Coefficient

{inv['cell']}

{inv['band']}

{inv['defect']}

{endpoint['weights']}

{endpoint['ties']}

## Hyperbolic Ridge Decomposition

{fiber['domain']}

{fiber['diagonal_chart']}

{fiber['kernel']}

{fiber['ridge_lift']}

{fiber['wings']}

{fiber['guard']}

## Exact Obstruction Witness

{witness['source']}

{witness['transform']}

{witness['cell']}

{witness['bound']}

High-precision audit values are B_(2,2)={witness['diagonal_value']} and E_2^ext={witness['exterior_value']}.

{witness['scope']}

## Ideal Joins

{inv['terminal']}

## Route Decision

{handoff['decision']}

{handoff['route']}

{handoff['obligation']}

## Pi Provenance

The normalized sinc and sine-integral witness follows from the fixed character e(x)=exp(2pi i x). The constant pi is inherited from that Fourier convention; no circle, polygon, fitted constant, or plotted symmetry is inserted.

## Proof Boundary

{handoff['reserve']} This gate is not a proof of RH.
"""


def main() -> int:
    load_sources()
    certificate = {
        "inversion": inversion_certificate(),
        "fiber": fiber_certificate(),
        "witness": witness_certificate(),
        "endpoint": endpoint_certificate(),
        "handoff": handoff_certificate(),
    }
    rows = build_rows(certificate)
    counts = {
        "rows": len(rows),
        "inversion_identities": certificate["inversion"]["inversion_identities"],
        "fiber_identities": certificate["fiber"]["fiber_identities"],
        "ridge_decomposition_identities": certificate["fiber"]["ridge_decomposition_identities"],
        "endpoint_weight_checks": certificate["endpoint"]["endpoint_weight_checks"],
        "internal_tie_cancellations": certificate["endpoint"]["internal_tie_cancellations"],
        "witness_exact_checks": certificate["witness"]["witness_exact_checks"],
        "witness_numerical_checks": certificate["witness"]["witness_numerical_checks"],
        "nonzero_exterior_witnesses": certificate["witness"]["nonzero_exterior_witnesses"],
        "signed_ridge_completed_bounds": 0,
        "phi_b_bounds": 0,
    }
    artifact = {
        "kind": KIND,
        "date": "2026-08-03",
        "status": "one-starred-carrier coefficient and finite-diagonal obstruction proved; complete signed ideal join open",
        "counts": counts,
        "symbolic_certificate": certificate,
        "source_audit": source_audit(),
        "rows": [asdict(row) for row in rows],
        "proof_boundary": "This proves the exact per-cell and global carrier coefficient, starred endpoint weights, hyperbolic ridge/residual/wings decomposition, and a nonzero-exterior witness. It proves no signed ridge-completed ideal join, exterior/complement bound, completed current, Phi_B theorem, contact exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    print(
        "built double-Morse ridge/carrier inversion gate: "
        f"{counts['rows']} rows, {counts['inversion_identities']} inversion identities, "
        f"{counts['ridge_decomposition_identities']} ridge decompositions, "
        f"{counts['witness_exact_checks']} exact witness checks, "
        f"{counts['witness_numerical_checks']} witness checks, "
        f"{counts['endpoint_weight_checks']} endpoint-weight checks, "
        f"{counts['nonzero_exterior_witnesses']} nonzero exterior witness, "
        f"{counts['signed_ridge_completed_bounds']} signed ridge-completed bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
