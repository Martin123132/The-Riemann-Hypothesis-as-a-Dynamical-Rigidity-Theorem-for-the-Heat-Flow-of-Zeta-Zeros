#!/usr/bin/env python3
"""Build the endpoint odd-fibre correlation feasibility guard."""

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
    "endpoint_odd_fibre_correlation_feasibility_guard"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "first_pivot_guard": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "endpoint_first_pivot_odd_small_ball_guard.json"
    ),
    "absolute_phase_anchor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "absolute_phase_anchor_reduction.json"
    ),
    "joined_odd_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "joined_dyadic_odd_prefix_first_jet_reduction.json"
    ),
    "endpoint_lift": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_endpoint_holomorphic_lift.json"
    ),
    "carrier_kernel": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "carrier_kernel_abel_prefix_reduction.json"
    ),
}


@dataclass(frozen=True)
class FeasibilityRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def source_hashes() -> dict[str, str]:
    missing = [str(path) for path in SOURCES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {key: file_hash(path) for key, path in SOURCES.items()}


def source_audit() -> dict[str, str]:
    payloads = {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCES.items()
    }
    markers = {
        "first_pivot_guard": (
            "Dbar(-R_N,rho_K)",
            "d_1 remains inside O_N",
        ),
        "absolute_phase_anchor": (
            "f_1=phi*(1+d_1)",
            "g_0=(-1)^N*beta*(T_0+i)*H_a/|M_t(s)|",
        ),
        "joined_odd_prefix": (
            "H_0,H_1,H_2,D_0,D_1",
            "M_k=floor(N/2^k)",
        ),
        "endpoint_lift": (
            "M_0(i*T)*U(T)*exp(pi*i/8)=",
            "rapidly oscillating Stirling and Riemann-Siegel phases",
        ),
        "carrier_kernel": (
            "mathcal_C_N",
            "|mathcal_C_N|>A_L+epsilon_term",
        ),
    }
    for key, required in markers.items():
        text = json.dumps(payloads[key], sort_keys=True)
        for marker in required:
            if marker not in text:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {
        key: str(payload.get("kind", ""))
        for key, payload in payloads.items()
    }


def correlation_audit() -> dict[str, str]:
    x_odd, y_odd, gamma, time_0 = sp.symbols(
        "X_odd Y_odd gamma T_0", real=True
    )
    subtotal = x_odd + sp.I * y_odd
    endpoint = gamma * (time_0 + sp.I)
    norm = sp.expand_complex(
        (subtotal + endpoint) * sp.conjugate(subtotal + endpoint)
    )
    expected = (
        (x_odd + gamma * time_0) ** 2
        + (y_odd + gamma) ** 2
    )
    if sp.simplify(norm - expected) != 0:
        raise RuntimeError("odd-fibre endpoint norm audit failed")

    correlation = sp.expand_complex(
        sp.re(subtotal * sp.conjugate(endpoint))
    )
    expected_correlation = gamma * (
        time_0 * x_odd + y_odd
    )
    if sp.simplify(correlation - expected_correlation) != 0:
        raise RuntimeError("odd-fibre endpoint correlation audit failed")

    scale = sp.sqrt(time_0**2 + 1)
    parallel = (time_0 * x_odd + y_odd) / scale
    perpendicular = (-x_odd + time_0 * y_odd) / scale
    rotated = (
        (parallel + gamma * scale) ** 2 + perpendicular**2
    )
    if sp.simplify(rotated - expected) != 0:
        raise RuntimeError("endpoint parallel-coordinate audit failed")

    return {
        "norm": (
            "|S_odd+g_0|^2=(X_odd+gamma_N*T_0)^2"
            "+(Y_odd+gamma_N)^2."
        ),
        "correlation": (
            "Re(S_odd*conj(g_0))="
            "gamma_N*(T_0*X_odd+Y_odd)."
        ),
        "rotation": (
            "With D=sqrt(T_0^2+1), "
            "P_odd=(T_0X_odd+Y_odd)/D and "
            "Q_odd=(-X_odd+T_0Y_odd)/D, "
            "|S_odd+g_0|^2=(P_odd+gamma_ND)^2+Q_odd^2."
        ),
    }


def moment_null_audit() -> dict:
    ell, h = sp.symbols("ell h", real=True)
    weights = [
        sp.Integer(1),
        sp.Integer(-3),
        sp.Integer(3),
        sp.Integer(-1),
    ]
    nodes = [ell + index * h for index in range(4)]
    moments = [
        sp.expand(
            sum(
                weight * node**order
                for weight, node in zip(weights, nodes)
            )
        )
        for order in range(4)
    ]
    if moments[:3] != [0, 0, 0]:
        raise RuntimeError("finite-difference moment nullspace failed")
    if sp.simplify(moments[3] + 6 * h**3) != 0:
        raise RuntimeError("finite-difference third moment drifted")
    return {
        "support": "n=3,6,12,24",
        "log_nodes": "log(3),log(3)+h,log(3)+2h,log(3)+3h",
        "weights": "1,-3,3,-1",
        "moments": (
            "sum_(k=0)^3 w_k(log(3)+kh)^r=0 for r=0,1,2; "
            "the r=3 sum is -6h^3."
        ),
        "odd_fibre_change": (
            "Only n=3 is odd, so the k=0 odd fibre changes by lambda."
        ),
        "terminal_guard": (
            "For N>=32 the terminal power 2^K is outside "
            "{3,6,12,24}, so its coefficient is unchanged."
        ),
    }


def linked_jet_audit() -> dict:
    w = sp.symbols("w")
    lam = sp.symbols("lambda", real=True)
    polynomial = 1 + lam * (w - 1) ** 3
    if sp.simplify(polynomial.subs(w, 1) - 1) != 0:
        raise RuntimeError("linked value guard failed")
    if sp.simplify(sp.diff(polynomial, w).subs(w, 1)) != 0:
        raise RuntimeError("linked first-jet guard failed")
    if sp.simplify(sp.diff(polynomial, w, 2).subs(w, 1)) != 0:
        raise RuntimeError("linked second-jet guard failed")
    constant = sp.expand(polynomial).coeff(w, 0)
    leading = sp.expand(polynomial).coeff(w, 3)
    pivot = sp.expand(constant**2 - leading**2)
    if pivot.subs(lam, sp.Rational(1, 4)) != sp.Rational(1, 2):
        raise RuntimeError("positive linked-jet pivot guard failed")
    if pivot.subs(lam, 1) != -1:
        raise RuntimeError("negative linked-jet pivot guard failed")
    return {
        "family": "F_lambda(w)=1+lambda*(w-1)^3",
        "linked_two_jet": (
            "F_lambda(1)=1, F_lambda'(1)=F_lambda''(1)=0 "
            "for every real lambda."
        ),
        "pivot": "Delta_3=|1-lambda|^2-|lambda|^2=1-2lambda.",
        "positive_case": "lambda=1/4 gives Delta_3=1/2.",
        "negative_case": "lambda=1 gives Delta_3=-1.",
    }


def build_exact() -> dict:
    correlation = correlation_audit()
    null_family = moment_null_audit()
    linked_jet = linked_jet_audit()
    return {
        "pi_provenance": (
            "The pi in T_0, beta, and the endpoint phase cancellation "
            "comes from the completed-zeta normalization and the "
            "Riemann-Siegel saddle. The phase restoration uses only the "
            "unit factor phi=M_t(s)/|M_t(s)|. The odd/even split, finite "
            "difference guard, and correlation rotation introduce no new pi."
        ),
        "physical_coefficients": (
            "Put phi=M_t(s)/|M_t(s)| and "
            "f_n=phi*exp[t*log(n)^2/4-s_*log(n)]*(1+d_n). "
            "Then f_1=phi*(1+d_1)."
        ),
        "physical_odd_fibre": (
            "For S_odd=sum_(m<=N,m odd)f_m, the odd numerator satisfies "
            "S_odd=phi*O_N."
        ),
        "physical_endpoint": (
            "The phase-cancelled recurrent endpoint is "
            "g_0=(-1)^N*beta*(T_0+i)*H_a/|M_t(s)|=phi*R_N. "
            "Define the real scalar "
            "gamma_N=(-1)^N*beta*H_a/|M_t(s)|, so g_0=gamma_N*(T_0+i)."
        ),
        "physical_terminal": (
            "For p_K=2^K, "
            "f_(p_K)=phi*exp(i*omega_*K*h)*C_K and "
            "|f_(p_K)|=|C_K|=rho_K."
        ),
        "physical_pivot": (
            "Delta_K={|S_odd+g_0|^2-|f_(p_K)|^2}/|1+d_1|^2. "
            "Thus the first Schur pivot is positive iff the physical "
            "odd-index subtotal plus endpoint has modulus larger than the "
            "terminal power-of-two carrier."
        ),
        "correlation": correlation["correlation"],
        "cartesian_pivot": (
            correlation["norm"]
            + " Hence the pivot numerator is "
            "(X_odd+gamma_N*T_0)^2+(Y_odd+gamma_N)^2"
            "-|f_(p_K)|^2."
        ),
        "endpoint_rotation": correlation["rotation"],
        "zero_endpoint_case": (
            "If H_a=0, then gamma_N=g_0=0 and the target reduces to "
            "|S_odd|>|f_(p_K)|. Endpoint phase cancellation supplies no "
            "odd-fibre lower bound in this case."
        ),
        "total_value_decomposition": (
            "With S_even=sum_(n<=N,n even)f_n, "
            "E_[1]=S_odd+S_even+g_0 and "
            "S_odd+g_0=E_[1]-S_even. Existing control of the total "
            "real value, centered scalar, or adjacent total projections "
            "does not by itself bound this difference."
        ),
        "moment_null_family": (
            "In the correction-free joined coefficient model with N>=32, "
            "perturb q_3,q_6,q_12,q_24 by "
            "lambda*(1,-3,3,-1). "
            + null_family["moments"]
            + " D_0=D_1=0 remains unchanged, the endpoint and terminal "
            "power-of-two carrier are fixed, but "
            + null_family["odd_fibre_change"]
        ),
        "moment_null_consequence": (
            "The complete five-current data H_0,H_1,H_2,D_0,D_1, "
            "together with a fixed endpoint and fixed terminal carrier, "
            "do not determine the first pivot in the unrestricted "
            "correction-free joined class. Choosing lambda to cancel the "
            "odd-endpoint subtotal makes the pivot numerator negative; "
            "choosing |lambda| sufficiently large makes it positive."
        ),
        "linked_two_jet_guard": (
            linked_jet["family"]
            + " satisfies "
            + linked_jet["linked_two_jet"]
            + " Yet "
            + linked_jet["positive_case"]
            + " while "
            + linked_jet["negative_case"]
            + " Thus even a linked value and two z derivatives do not "
            "determine a global Schur pivot."
        ),
        "nonpromotion_boundary": (
            "The moment-null and linked-two-jet families use unrestricted "
            "coefficient perturbations. They prove insufficiency of the "
            "currently retained observables, not failure of the actual Xi "
            "coefficient family. Positivity, heat-shift, correction, cutoff, "
            "and endpoint relations may help only through a new proved "
            "Xi-specific inequality."
        ),
        "feasibility_verdict": (
            "The source audit removes the phase mystery: after restoring "
            "phi, the endpoint is the real scalar gamma_N times T_0+i, and "
            "the missing correlation is exactly "
            "gamma_N*(T_0X_odd+Y_odd). No current theorem controls that "
            "odd-fibre projection at the terminal-carrier scale. The "
            "five-current closure and the retained linked two-jet are "
            "algebraically insufficient without additional Xi coefficient "
            "structure. Full Schur disk stability is therefore downgraded "
            "from the primary proof route to a falsification coordinate."
        ),
        "retained_schur_role": (
            "Retain Delta_K for a later one-worker physical scout or for "
            "testing any source-derived odd-fibre inequality. Do not expand "
            "the second pivot unless a uniform actual first-pivot theorem "
            "is proved. A selected positive table cannot reverse this "
            "route decision."
        ),
        "signed_primary_target": (
            "Return to the actual linked point. On q>=1, prove "
            "|mathsf_X|<=delta_L implies "
            "|mathcal_C_N|>A_L+epsilon_term, including W_0=0, and prove "
            "that crossings mathsf_X=0 with mathcal_C_N>0 number less than "
            "one successor turn after complete boundary composition. This "
            "uses the endpoint-complete first jet without demanding a "
            "coefficient fibre dominate on every free disk point."
        ),
        "q_lt_1_target": (
            "The q=2tL^2<1 multiplicity-compatible parabolic/Hermite chart "
            "and finite connectors remain separate. This feasibility audit "
            "does not close them."
        ),
        "proof_boundary": (
            "The physical rephasing, odd-fibre and endpoint subtotal, "
            "terminal-carrier identity, Cartesian and rotated correlation "
            "forms, total/even decomposition, correction-free five-current "
            "null family, and linked-two-jet guard are exact. The verdict "
            "is a no-go for the current Schur inputs, not a proof that no "
            "future Xi-specific first-pivot theorem exists. No uniform "
            "actual Xi pivot, physical pivot failure, signed prefix lower "
            "bound, crossing count, strict successor flux theorem, q<1 "
            "closure, contact exclusion, Lambda<=0, PF-infinity, RH proof, "
            "or Clay-prize conclusion is asserted."
        ),
        "diagnostics": {
            "moment_null": null_family,
            "linked_two_jet": linked_jet,
        },
    }


def build_rows(exact: dict) -> list[FeasibilityRow]:
    return [
        FeasibilityRow(
            "eofcfg_00_pi_provenance",
            "definition_provenance",
            "certified",
            "The correlation audit introduces no unexplained pi.",
            exact["pi_provenance"],
            "The odd-fibre geometry is independent of circle-area formulas.",
        ),
        FeasibilityRow(
            "eofcfg_01_physical_coefficients",
            "exact_rephasing",
            "ready_to_apply",
            "The absolute normalizer phase restores the physical carriers.",
            exact["physical_coefficients"],
            "The actual d_n corrections remain.",
        ),
        FeasibilityRow(
            "eofcfg_02_odd_fibre",
            "exact_rephasing",
            "certified",
            "The first Schur constant fibre is the physical odd-index subtotal.",
            exact["physical_odd_fibre"],
            "It is not the full corrected main value.",
        ),
        FeasibilityRow(
            "eofcfg_03_endpoint",
            "exact_phase_reduction",
            "certified",
            "The recurrent endpoint is a real scalar times T_0+i.",
            exact["physical_endpoint"],
            "The scalar may vanish or change sign with H_a and parity.",
        ),
        FeasibilityRow(
            "eofcfg_04_terminal",
            "exact_rephasing",
            "certified",
            "The pivot radius is the physical terminal power-of-two carrier modulus.",
            exact["physical_terminal"],
            "Its d_(2^K) correction remains.",
        ),
        FeasibilityRow(
            "eofcfg_05_physical_pivot",
            "exact_equivalence",
            "certified",
            "The first pivot compares one physical subtotal with one physical carrier.",
            exact["physical_pivot"],
            "No lower bound follows from the identity alone.",
        ),
        FeasibilityRow(
            "eofcfg_06_correlation",
            "exact_identity",
            "certified",
            "The endpoint correlation is one explicit odd-fibre projection.",
            exact["correlation"],
            "Its sign is not fixed by endpoint phase cancellation.",
        ),
        FeasibilityRow(
            "eofcfg_07_cartesian_pivot",
            "exact_identity",
            "ready_to_apply",
            "The pivot is a translated Euclidean radius inequality.",
            exact["cartesian_pivot"],
            "Both odd-fibre coordinates remain uncontrolled.",
        ),
        FeasibilityRow(
            "eofcfg_08_endpoint_rotation",
            "exact_identity",
            "ready_to_apply",
            "Endpoint-parallel and perpendicular coordinates isolate the missing bound.",
            exact["endpoint_rotation"] + " " + exact["zero_endpoint_case"],
            "The zero-endpoint branch cannot be removed.",
        ),
        FeasibilityRow(
            "eofcfg_09_total_decomposition",
            "exact_decomposition",
            "certified",
            "The odd-endpoint subtotal is total value minus the even subtotal.",
            exact["total_value_decomposition"],
            "Existing total projection bounds do not isolate the even subtotal.",
        ),
        FeasibilityRow(
            "eofcfg_10_moment_null_family",
            "countermodel",
            "guard_validated",
            "A four-term dyadic finite difference hides an arbitrary odd-fibre change.",
            exact["moment_null_family"],
            "The family is correction-free and unrestricted, not actual Xi.",
            exact["diagnostics"]["moment_null"],
        ),
        FeasibilityRow(
            "eofcfg_11_five_current_guard",
            "countermodel_consequence",
            "guard_validated",
            "The five retained currents do not determine the first pivot generically.",
            exact["moment_null_consequence"],
            "New Xi coefficient structure could still close the gap.",
        ),
        FeasibilityRow(
            "eofcfg_12_linked_two_jet_guard",
            "countermodel",
            "guard_validated",
            "A linked value and two local derivatives do not determine a Schur pivot.",
            exact["linked_two_jet_guard"],
            "This is a generic polynomial guard.",
            exact["diagnostics"]["linked_two_jet"],
        ),
        FeasibilityRow(
            "eofcfg_13_nonpromotion",
            "nonpromotion_guard",
            "guard_validated",
            "The exact guards are not promoted to an Xi counterexample.",
            exact["nonpromotion_boundary"],
            "Only insufficiency of current inputs is asserted.",
        ),
        FeasibilityRow(
            "eofcfg_14_feasibility_verdict",
            "route_decision",
            "guard_validated",
            "The current corpus supplies no viable first-pivot correlation theorem.",
            exact["feasibility_verdict"],
            "A genuinely new Xi odd-fibre theorem could reopen the route.",
        ),
        FeasibilityRow(
            "eofcfg_15_retained_schur_role",
            "conditional_route",
            "not_ready_to_apply",
            "Schur-Cohn remains an exact falsification coordinate.",
            exact["retained_schur_role"],
            "No later pivot is expanded.",
        ),
        FeasibilityRow(
            "eofcfg_16_signed_primary_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The primary q>=1 target returns to endpoint-complete signed transport.",
            exact["signed_primary_target"],
            "The prefix lower bound and signed count remain open.",
        ),
        FeasibilityRow(
            "eofcfg_17_q_lt_1_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The multiplicity-compatible small-q chart remains separate.",
            exact["q_lt_1_target"],
            "No small-q or connector theorem follows.",
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    rows = build_rows(exact)
    return {
        "kind": STEM,
        "date": "2026-07-28",
        "status": (
            "exact physical odd-fibre endpoint-correlation normal form, "
            "five-current null-family guard, and linked-two-jet guard; "
            "full Schur disk stability is downgraded to a falsification "
            "coordinate and the signed linked-point theorem remains open"
        ),
        "proof_boundary": exact["proof_boundary"],
        "sources": {
            key: str(path.relative_to(REPO_ROOT))
            for key, path in SOURCES.items()
        },
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return "\n".join(
        [
            "# Newman Endpoint Odd-Fibre Correlation Feasibility Guard",
            "",
            "Date: 2026-07-28",
            "",
            "Status: exact correlation audit and current-input no-go;",
            "not a proof of an actual Xi pivot, `Lambda<=0`, PF-infinity,",
            "RH, or a Clay-prize result.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Pi Provenance",
            "",
            exact["pi_provenance"],
            "",
            "## Physical Rephasing",
            "",
            "```text",
            exact["physical_coefficients"],
            exact["physical_odd_fibre"],
            exact["physical_endpoint"],
            exact["physical_terminal"],
            "```",
            "",
            "## First Pivot",
            "",
            "```text",
            exact["physical_pivot"],
            exact["correlation"],
            exact["cartesian_pivot"],
            exact["endpoint_rotation"],
            "```",
            "",
            exact["zero_endpoint_case"],
            "",
            "## Total Versus Odd Fibre",
            "",
            "```text",
            exact["total_value_decomposition"],
            "```",
            "",
            "## Five-Current Null Family",
            "",
            "```text",
            exact["moment_null_family"],
            exact["moment_null_consequence"],
            "```",
            "",
            "## Linked Two-Jet Guard",
            "",
            "```text",
            exact["linked_two_jet_guard"],
            "```",
            "",
            exact["nonpromotion_boundary"],
            "",
            "## Feasibility Verdict",
            "",
            exact["feasibility_verdict"],
            "",
            "## Retained Schur Role",
            "",
            exact["retained_schur_role"],
            "",
            "## Primary Target",
            "",
            "```text",
            exact["signed_primary_target"],
            exact["q_lt_1_target"],
            "```",
            "",
            "## Boundary",
            "",
            exact["proof_boundary"],
            "",
        ]
    )


def write_outputs(artifact: dict, output: Path, note: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    note.write_text(render_note(artifact), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    write_outputs(artifact, args.output, args.note)
    print(
        "built Newman endpoint odd-fibre correlation feasibility guard: "
        "18 rows, 1 exact physical odd-fibre pivot, "
        "1 endpoint-projection normal form, "
        "1 five-current null-family guard, 1 linked-two-jet guard, "
        "1 Schur route downgrade, 2 open routes"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
