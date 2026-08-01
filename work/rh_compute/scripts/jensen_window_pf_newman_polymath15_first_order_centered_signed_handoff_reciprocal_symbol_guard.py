#!/usr/bin/env python3
"""Build the signed zeta-handoff reciprocal-symbol no-gain guard."""

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
    "signed_handoff_reciprocal_symbol_guard"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"

POLYMATH_SOURCE_URL = "https://arxiv.org/abs/1904.12438"
EXPONENT_PAIR_SOURCE_URL = "https://arxiv.org/abs/2306.05599"

R_STAR = Fraction(125_662, 155_153)
R_STAR_DUAL = 2 - R_STAR
C_STAR = Fraction(4_911_678_521, 1_933_561_194)
C_TWO = Fraction(2, 1)
C_TWO_DEFICIT = Fraction(3_133_668_399, 48_144_906_818)

SOURCES = {
    "oscillatory_handoff": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "oscillatory_zeta_handoff_theorem.json"
    ),
    "reciprocal_self_duality": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "reciprocal_saddle_self_duality_gate.json"
    ),
    "normalizer_phase": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "reciprocal_normalizer_phase_reinforcement_guard.json"
    ),
    "first_order_remainder": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_first_order_global_remainder_certificate.json"
    ),
    "adjacent_recurrence": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "adjacent_saddle_recurrence.json"
    ),
}


@dataclass(frozen=True)
class SymbolRow:
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
        "oscillatory_handoff": (
            "D_k=sum_(n<=N)exp((t/4)log(n)^2)",
            "sum_(n>N)log(n)^k*n^(-s_*)",
            "4911678521/1933561194",
        ),
        "reciprocal_self_duality": (
            "u_nu=a_omega^2/nu",
            "stationary factor u_nu/a_omega",
            "r_*=125662/155153",
        ),
        "normalizer_phase": (
            "|R-1|<8/x",
            "signed weighted-minus-unweighted tail kernel",
        ),
        "first_order_remainder": (
            "T|d|<100",
            "C_1 parity",
            "200000*L*exp(-5L/4)",
        ),
        "adjacent_recurrence": (
            "last saddle and endpoint cannot be bounded separately",
            "adjacent-saddle recurrence != contact exclusion",
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


def fraction_record(value: Fraction) -> dict[str, str]:
    return {
        "exact": str(value),
        "decimal": f"{float(value):.15f}",
    }


def symbolic_audit() -> dict:
    heat_time, saddle_log, dual_offset, sigma = sp.symbols(
        "t A z sigma", real=True
    )
    primal_log = saddle_log - dual_offset
    dual_log = saddle_log + dual_offset
    q_primal = heat_time * primal_log**2 / 4
    q_dual = heat_time * dual_log**2 / 4
    defect = 2 * sigma - 1 - heat_time * saddle_log

    weighted_left = sp.expand(q_primal + (2 * sigma - 1) * dual_offset)
    weighted_right = sp.expand(q_dual + defect * dual_offset)
    if sp.simplify(weighted_left - weighted_right) != 0:
        raise RuntimeError("signed reciprocal weighted identity failed")

    signed_kernel = (
        sp.exp((2 * sigma - 1) * dual_offset)
        * (sp.exp(q_primal) - 1)
    )
    weighted_kernel = sp.exp(q_dual + defect * dual_offset)
    normalized_symbol = sp.simplify(signed_kernel / weighted_kernel)
    expected_symbol = 1 - sp.exp(-q_primal)
    if sp.simplify(normalized_symbol - expected_symbol) != 0:
        raise RuntimeError("normalized signed symbol failed")

    log_nu, sigma_x, saddle_log_x = sp.symbols(
        "m sigma_x A_x", real=True
    )
    local_primal_log = 2 * saddle_log - log_nu
    local_dual_offset = log_nu - saddle_log
    local_q = heat_time * local_primal_log**2 / 4
    local_log_kernel = (
        (2 * sigma - 1) * local_dual_offset
        + sp.log(sp.exp(local_q) - 1)
    )
    directional = sp.diff(local_log_kernel, sigma) * sigma_x
    directional += sp.diff(local_log_kernel, saddle_log) * saddle_log_x
    expected_directional = (
        2 * sigma_x * local_dual_offset
        - (2 * sigma - 1) * saddle_log_x
        + heat_time
        * local_primal_log
        * saddle_log_x
        * sp.exp(local_q)
        / (sp.exp(local_q) - 1)
    )
    if sp.simplify(directional - expected_directional) != 0:
        raise RuntimeError("signed kernel first derivative failed")

    symbol_log_derivative = sp.simplify(
        sp.diff(sp.log(1 - sp.exp(-local_q)), saddle_log)
        * saddle_log_x
    )
    expected_symbol_derivative = (
        heat_time
        * local_primal_log
        * saddle_log_x
        / (sp.exp(local_q) - 1)
    )
    if sp.simplify(
        symbol_log_derivative - expected_symbol_derivative
    ) != 0:
        raise RuntimeError("normalized symbol derivative failed")

    weighted_dual_exponent = C_TWO * R_STAR_DUAL**2 / 8
    ordinary_dual_exponent = C_TWO * (1 - R_STAR) / 2
    relative_gap_c2 = C_TWO * R_STAR**2 / 8
    relative_gap_cstar = C_STAR * R_STAR**2 / 8
    if weighted_dual_exponent - ordinary_dual_exponent != relative_gap_c2:
        raise RuntimeError("active dual exponent decomposition failed")
    if R_STAR_DUAL - 1 != 1 - R_STAR:
        raise RuntimeError("cutoff separation exponent failed")

    endpoint_separation_c2 = (
        C_TWO_DEFICIT + Fraction(1, 2) + C_TWO / 8
    )

    return {
        "weighted_exponent_identity": (
            "q(u_nu)+(2*sigma-1)log(nu/a_omega)"
            "=q(nu)+Delta*log(nu/a_omega)"
        ),
        "normalized_symbol_identity": (
            "K_signed/[w(nu)*(nu/a_omega)^Delta]"
            "=1-exp(-q(u_nu))"
        ),
        "kernel_log_derivative": (
            "K_x/K=2*sigma_x*log(nu/a_omega)"
            "-(2*sigma-1)A_x"
            "+t*log(u_nu)*A_x*exp(q)/(exp(q)-1)"
        ),
        "symbol_log_derivative": (
            "partial_x log(1-exp(-q))"
            "=t*log(u_nu)*A_x/(exp(q)-1)"
        ),
        "constants": {
            "r_star": fraction_record(R_STAR),
            "r_star_dual": fraction_record(R_STAR_DUAL),
            "c_star": fraction_record(C_STAR),
            "c_two_deficit": fraction_record(C_TWO_DEFICIT),
            "c2_weighted_dual_exponent": fraction_record(
                weighted_dual_exponent
            ),
            "c2_unweighted_dual_exponent": fraction_record(
                ordinary_dual_exponent
            ),
            "c2_relative_gap": fraction_record(relative_gap_c2),
            "cstar_relative_gap": fraction_record(relative_gap_cstar),
            "cutoff_separation_exponent": fraction_record(1 - R_STAR),
            "endpoint_separation_at_c2": fraction_record(
                endpoint_separation_c2
            ),
        },
    }


def build_exact() -> dict:
    symbolic = symbolic_audit()
    constants = symbolic["constants"]
    return {
        "pi_provenance": (
            "The 2*pi in a_omega^2=omega/(2*pi) comes from the standard "
            "Poisson phase exp(-2*pi*i*nu*u), and the pi/4 in the global "
            "factor is the negative-curvature stationary signature. These "
            "are the same completed-zeta and Riemann-Siegel constants used "
            "by the Polymath-15 approximation; no new geometric pi is "
            "introduced."
        ),
        "signed_handoff": (
            "For w_t(u)=exp[t*log(u)^2/4], the exact zeta-handoff "
            "difference is S_N=sum_(n>=1)F_N(n)n^(-s_*), where "
            "F_N(u)=w_t(u)*1_(u<=N)-1. Equivalently, "
            "S_N=sum_(n<=N)(w_t(n)-1)n^(-s_*)"
            "-sum_(n>N)n^(-s_*)."
        ),
        "stationary_term": (
            "Write s_*=sigma-i*omega, a_omega^2=omega/(2*pi), "
            "y=nu/a_omega, and u_nu=a_omega/y=a_omega^2/nu. "
            "For a positive Poisson mode nu, stationary phase gives the "
            "leading term Q(omega)*K_N(nu)*nu^(-conj(s_*)), where "
            "Q(omega)=exp(i[omega*(2log(a_omega)-1)-pi/4]) and "
            "K_N(nu)=y^(2sigma-1)F_N(u_nu)."
        ),
        "interior_signed_symbol": (
            "On the interior reciprocal side nu>a_omega^2/N, so "
            "u_nu<N and F_N(u_nu)=w_t(u_nu)-1. Hence "
            "K_N=y^(2sigma-1)(w_t(u_nu)-1)"
            "=w_t(nu)*y^Delta*[1-w_t(u_nu)^(-1)], "
            "Delta=2sigma-1-t*log(a_omega). This factorization is exact."
        ),
        "normalized_symbol": (
            "After division by the dominant weighted reciprocal saddle, "
            "the signed weighted-minus-unweighted symbol is exactly "
            "Sigma_N(nu)=1-exp[-t*log(u_nu)^2/4]. It lies strictly between "
            "zero and one for t>0 and u_nu!=1. Thus it has neither a sign "
            "change nor a zero on the active positive-radius block."
        ),
        "active_scale": (
            "For u_nu=N^r and c=tL with log(N)=L/2+o(1), "
            "exp[-t*log(u_nu)^2/4]=N^(-c*r^2/8+o(1)). At "
            "r=r_*=125662/155153, Sigma_N=1-N^(-c*r_*^2/8+o(1)). "
            "The relative defect therefore tends to zero while the signed "
            "symbol tends to one."
        ),
        "c2_audit": (
            "At c=2 the weighted reciprocal saddle at radius "
            "2-r_* has N-exponent "
            f"{constants['c2_weighted_dual_exponent']['exact']}, while "
            "the subtracted ordinary reciprocal term has exponent "
            f"{constants['c2_unweighted_dual_exponent']['exact']}. "
            "Their exact gap is "
            f"{constants['c2_relative_gap']['exact']}="
            "r_*^2/4=0.163993860281979..., so the subtraction is a "
            "relative power-small correction. It supplies zero exponent "
            "gain against the existing positive deficit d_2."
        ),
        "cstar_audit": (
            "At c=c_* the relative exponent gap is "
            f"{constants['cstar_relative_gap']['exact']}="
            "0.208290568620833..., again with the weighted reciprocal "
            "term dominant. The conclusion is not an artifact of testing "
            "only c=2."
        ),
        "kernel_derivative": (
            "On a fixed-t, fixed-N chart and for fixed mode nu, put "
            "A=log(a_omega), U=log(u_nu)=2A-log(nu), "
            "z=log(nu/a_omega), q=tU^2/4, and "
            "A_x=omega_x/(2omega). Then "
            "K_x/K=2sigma_x*z-(2sigma-1)A_x"
            "+tU*A_x*exp(q)/(exp(q)-1). This is the exact first-x "
            "derivative of the interior leading amplitude."
        ),
        "full_term_derivative": (
            "For T_nu=Q(omega)K_N(nu)nu^(-conj(s_*)), "
            "T_(nu,x)/T_nu=-sigma_x*U-(2sigma-1)A_x"
            "+tU*A_x*exp(q)/(exp(q)-1)+i*omega_x*U. "
            "The signed factor alone obeys "
            "partial_x log(Sigma_N)=tU*A_x/(exp(q)-1), which is "
            "power-small on the active block. The first jet therefore "
            "inherits the same no-gain verdict."
        ),
        "cutoff_guard": (
            "The reciprocal cutoff is nu=a_omega^2/N, asymptotically "
            "radius one. The active dual radius is "
            "2-r_*=184644/155153 and its exact power separation from the "
            "cutoff is 1-r_*=29491/155153=0.190076891842246.... "
            "The floor seam and adjacent recurrence cannot be identified "
            "with this interior block."
        ),
        "correction_endpoint_guard": (
            "The certified first correction has T|d_n|<100 with "
            "T asymptotic to pi*N^2, so it is relative O(N^-2). The "
            "Polymath-15 endpoint and one adjacent block have scale "
            "N^(-1/2-c/8+o(1)); at c=2 this is N^(-3/4+o(1)), separated "
            "from the current N^(d_2+o(1)) active-block envelope by the "
            f"exact exponent {constants['endpoint_separation_at_c2']['exact']}"
            "=0.815088263870695.... First derivatives add logarithms only. "
            "These terms remain mandatory at the cutoff, but they cannot "
            "supply a uniform leading-symbol cancellation at r_*."
        ),
        "route_falsification": (
            "The proposed signed reciprocal mechanism is falsified at its "
            "first required test. On the active block the exact normalized "
            "symbol tends to one, its x derivative is a power-small "
            "perturbation, the normalizer is conjugate-locked, and the "
            "endpoint/correction scales do not reverse the leading symbol. "
            "No N^(-d_2-eta) pairwise gain is produced."
        ),
        "surviving_route": (
            "Retire signed primal/tail pairwise cancellation as a route to "
            "lower c_*. The reciprocal coordinate may still reorganize a "
            "dual exponential sum, but any improvement now requires a new "
            "arithmetic estimate inside that dual sum, such as bilinear or "
            "additive-energy cancellation. The primary proof programme "
            "returns to the direct Xi Abel-phase/contact theorem unless such "
            "an input is derived independently."
        ),
        "proof_boundary": (
            "This artifact proves the exact continuous signed reciprocal "
            "symbol, its fixed-chart first-x derivative, the active rational "
            "exponent audit, and a no-gain verdict for signed pairwise "
            "cancellation. It does not prove a discrete endpoint-complete "
            "Poisson/B-process remainder theorem or rule out new arithmetic "
            "cancellation inside the dual sum. The effective L_epsilon, "
            "Abel-scalar gap, inner degree theorem, contact exclusion, final "
            "Newman conclusion, PF-infinity theorem, RH proof, and "
            "prize-level conclusion remain open."
        ),
        "diagnostics": symbolic,
    }


def build_rows(exact: dict) -> list[SymbolRow]:
    constants = exact["diagnostics"]["constants"]
    return [
        SymbolRow(
            "shrs_00_pi_provenance",
            "source_provenance",
            "available_exact",
            "The signed reciprocal symbol introduces no unexplained pi.",
            exact["pi_provenance"],
            "Standard completed-zeta, Poisson, and stationary-phase constants.",
        ),
        SymbolRow(
            "shrs_01_signed_handoff",
            "exact_reindexing",
            "available_exact",
            "The weighted-minus-unweighted handoff is one piecewise sum.",
            exact["signed_handoff"],
            "Absolutely literal where the zeta series converges; elsewhere a blockwise oscillatory identity.",
        ),
        SymbolRow(
            "shrs_02_stationary_term",
            "exact_stationary_algebra",
            "available_exact",
            "The reciprocal leading term has one scalar symbol.",
            exact["stationary_term"],
            "Continuous leading stationary term, not a discrete remainder theorem.",
        ),
        SymbolRow(
            "shrs_03_interior_symbol",
            "exact_algebraic_lemma",
            "available_exact",
            "The interior signed kernel factors through the weighted saddle.",
            exact["interior_signed_symbol"],
            "Exact away from the cutoff indicator seam.",
            exact["diagnostics"],
        ),
        SymbolRow(
            "shrs_04_normalized_symbol",
            "exact_algebraic_lemma",
            "available_exact",
            "The normalized signed symbol is positive and nonvanishing.",
            exact["normalized_symbol"],
            "For t>0 and the active block u_nu>1.",
            exact["diagnostics"],
        ),
        SymbolRow(
            "shrs_05_active_scale",
            "proved_asymptotic_estimate",
            "available_exact",
            "At r_* the signed symbol tends to one.",
            exact["active_scale"],
            "Uses log(N)=L/2+o(1) from the checked saddle geometry.",
            constants,
        ),
        SymbolRow(
            "shrs_06_c2_audit",
            "exact_rational_certificate",
            "available_exact",
            "At c=2 the ordinary subtraction is power-smaller.",
            exact["c2_audit"],
            "Exact rational exponent comparison.",
            constants,
        ),
        SymbolRow(
            "shrs_07_cstar_audit",
            "exact_rational_certificate",
            "available_exact",
            "The same no-gain sign persists at c_*.",
            exact["cstar_audit"],
            "Exact rational exponent comparison.",
            constants,
        ),
        SymbolRow(
            "shrs_08_kernel_derivative",
            "exact_first_jet_identity",
            "available_exact",
            "The interior signed amplitude has an exact first-x derivative.",
            exact["kernel_derivative"],
            "Fixed t, fixed N, and fixed reciprocal mode.",
            exact["diagnostics"],
        ),
        SymbolRow(
            "shrs_09_full_term_derivative",
            "exact_first_jet_identity",
            "available_exact",
            "The normalized signed correction is power-small in the first jet.",
            exact["full_term_derivative"],
            "The derivative does not supply a hidden leading cancellation.",
        ),
        SymbolRow(
            "shrs_10_cutoff_guard",
            "nonpromotion_guard",
            "guard_validated",
            "The active reciprocal block is power-separated from the seam.",
            exact["cutoff_guard"],
            "The adjacent recurrence remains mandatory at radius one.",
            constants,
        ),
        SymbolRow(
            "shrs_11_endpoint_guard",
            "nonpromotion_guard",
            "guard_validated",
            "Corrections and endpoints cannot reverse the active leading symbol uniformly.",
            exact["correction_endpoint_guard"],
            "Exceptional small block values are not excluded; only leading-symbol cancellation is rejected.",
            constants,
        ),
        SymbolRow(
            "shrs_12_route_decision",
            "theorem_search_decision",
            "available_exact",
            "Signed reciprocal pairwise cancellation is retired.",
            exact["route_falsification"] + " " + exact["surviving_route"],
            "A genuinely new dual-sum arithmetic estimate is not ruled out.",
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    return {
        "kind": STEM,
        "date": "2026-07-28",
        "status": (
            "exact signed reciprocal leading-symbol and first-jet "
            "no-gain guard; signed pairwise cancellation is retired"
        ),
        "sources": {
            key: str(path.relative_to(REPO_ROOT))
            for key, path in SOURCES.items()
        },
        "source_urls": {
            "polymath15_primary": POLYMATH_SOURCE_URL,
            "exponent_pairs_primary": EXPONENT_PAIR_SOURCE_URL,
        },
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "exact": exact,
        "rows": [asdict(row) for row in build_rows(exact)],
        "proof_boundary": exact["proof_boundary"],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    constants = exact["diagnostics"]["constants"]
    return f"""# Newman Signed-Handoff Reciprocal-Symbol Guard

Date: 2026-07-28

Status: exact continuous leading-symbol and first-jet no-gain guard.
This is not a proof of a discrete B-process theorem, the Newman
boundary, PF-infinity, RH, or a prize-level result.

```text
work/rh_compute/results/{STEM}.json
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```

Primary sources:

```text
{POLYMATH_SOURCE_URL}
{EXPONENT_PAIR_SOURCE_URL}
```

## Exact Signed Handoff

```text
{exact["signed_handoff"]}
```

## Pi Provenance

{exact["pi_provenance"]}

## Reciprocal Symbol

```text
{exact["stationary_term"]}
{exact["interior_signed_symbol"]}
{exact["normalized_symbol"]}
```

## Active Radius

```text
{exact["active_scale"]}
{exact["c2_audit"]}
{exact["cstar_audit"]}
```

```json
{json.dumps(constants, indent=2, sort_keys=True)}
```

## First X Derivative

```text
{exact["kernel_derivative"]}
{exact["full_term_derivative"]}
```

## Endpoint And Cutoff

```text
{exact["cutoff_guard"]}
{exact["correction_endpoint_guard"]}
```

## Route Decision

```text
{exact["route_falsification"]}
{exact["surviving_route"]}
```

## Boundary

{exact["proof_boundary"]}
"""


def write_artifact(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    write_artifact(args.out, artifact)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman signed-handoff reciprocal-symbol guard: "
        "13 rows, 3 exact symbol identities, 2 first-jet identities, "
        "2 rational exponent audits, 2 endpoint/cutoff guards, "
        "1 retired pairwise route"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
