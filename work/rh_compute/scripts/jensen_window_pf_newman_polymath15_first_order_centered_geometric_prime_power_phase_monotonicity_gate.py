#!/usr/bin/env python3
"""Build the geometric prime-power phase-monotonicity gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_centered_geometric_prime_power_"
    "phase_monotonicity_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "prime_power_heat": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_prime_power_heat_block_composition_guard.json"
    ),
    "complete_chain": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_complete_prime_power_chain_occupation_gate.json"
    ),
    "ray_flow": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_ray_bottom_logarithmic_flow_reduction.json"
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


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCES.items()
    }


def audit_sources(payloads: dict[str, dict]) -> dict:
    heat = payloads["prime_power_heat"].get("exact", {})
    if "Q_R(z)=sum_(k=0)^R a_k z^k" not in heat.get(
        "uncorrected_block_winding", ""
    ):
        raise RuntimeError("prime-power one-block source drifted")
    if "winding-three sum" not in heat.get("composition_guard", ""):
        raise RuntimeError("two-block composition guard drifted")

    chain = payloads["complete_chain"]
    if "Q_(m,p)(w)=sum_(k=0)^R_m" not in chain.get(
        "chain", {}
    ).get("chain_polynomial", ""):
        raise RuntimeError("complete-chain polynomial drifted")
    if "external phase B_m" not in chain.get(
        "threshold_guard", {}
    ).get("external_phase_flip", ""):
        raise RuntimeError("external-phase guard drifted")

    ray = payloads["ray_flow"].get("exact", {})
    if "mathsf_X*(L*D_j mathcal_C_N-mathcal_C_N)" not in (
        ray.get("ray_proxy", {}).get("defect_flux", "")
    ):
        raise RuntimeError("ray signed defect drifted")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCES.items()
        },
        "source_kinds": {
            key: payload.get("kind") for key, payload in payloads.items()
        },
    }


def build_rows() -> list[GateRow]:
    return [
        GateRow(
            "gppm_01_model",
            "correction-free chain",
            "proved",
            "Isolate one complete geometric prime-power block.",
            "a=p^(-1/2), P_M(z)=z*sum_(k=0)^(M-1)(a*z)^k.",
            "This is the t=0 correction-free internal chain model.",
        ),
        GateRow(
            "gppm_02_factorization",
            "exact algebra",
            "proved",
            "The finite block has one rational geometric representation.",
            "P_M(z)=z*[1-(a*z)^M]/(1-a*z).",
            "Both denominator and numerator are nonzero on |z|=1.",
        ),
        GateRow(
            "gppm_03_block_winding",
            "topological baseline",
            "proved",
            "Every geometric block is unit-circle nonzero with winding one.",
            "wind(P_M(e^(i*theta)),0)=1.",
            "Winding one does not imply monotone phase or one upward crossing.",
        ),
        GateRow(
            "gppm_04_kernel",
            "exact phase derivative",
            "proved",
            "Differentiate the two geometric factors without an argument branch.",
            "K_m(r,theta)=m*(r^2-r*cos(m*theta))/"
            "(1-2*r*cos(m*theta)+r^2).",
            "K_m is the angular current of 1-r*exp(i*m*theta).",
        ),
        GateRow(
            "gppm_05_phase_speed",
            "exact phase derivative",
            "proved",
            "The block phase speed is an explicit difference of two kernels.",
            "Theta_M'=1+K_M(a^M,theta)-K_1(a,theta).",
            "The external constant phase does not change this derivative.",
        ),
        GateRow(
            "gppm_06_kernel_bounds",
            "exact inequality",
            "proved",
            "Each geometric-factor current has sharp elementary bounds.",
            "-m*r/(1-r)<=K_m(r,theta)<=m*r/(1+r).",
            "The endpoint values cos(m*theta)=1 and -1 attain the bounds.",
        ),
        GateRow(
            "gppm_07_monotone_margin",
            "sufficient theorem",
            "proved",
            "A positive scalar margin forces strict one-turn phase.",
            "Theta_M'>=mu_M(a):=1/(1+a)-M*a^M/(1-a^M).",
            "The bound is sufficient and may be negative for a monotone block.",
        ),
        GateRow(
            "gppm_08_projection_consequence",
            "one-block theorem",
            "conditional",
            "When mu_M(a)>0, every constant-phase real projection has one upcrossing.",
            "For all alpha, Re(exp(i*alpha)P_M(e^(i*theta))) "
            "has two simple zeros and one upward zero on [0,2*pi).",
            "This is internal to one geometric block.",
        ),
        GateRow(
            "gppm_09_length_monotonicity",
            "threshold arithmetic",
            "proved",
            "The tail penalty decreases once a is below M/(M+1).",
            "g_M(a)=M*a^M/(1-a^M); "
            "a<=M/(M+1) => g_(M+1)(a)<g_M(a).",
            "This propagates one checked base margin to all longer blocks.",
        ),
        GateRow(
            "gppm_10_large_primes",
            "prime-family theorem",
            "proved",
            "Every correction-free p>=5 block has strict phase monotonicity.",
            "M=1 is exact; for M>=2, "
            "mu_M(p^(-1/2))>=mu_2(5^(-1/2))"
            "=(3-sqrt(5))/4>0.",
            "Heat and d_n perturbations are not included.",
        ),
        GateRow(
            "gppm_11_ternary_threshold",
            "prime-family theorem",
            "proved",
            "Every correction-free ternary block of length at least four is monotone.",
            "mu_M(3^(-1/2))>=mu_4(3^(-1/2))"
            "=(2-sqrt(3))/2>0 for M>=4.",
            "Lengths two and three are outside this theorem.",
        ),
        GateRow(
            "gppm_12_dyadic_threshold",
            "prime-family theorem",
            "proved",
            "Every correction-free dyadic block of length at least eight is monotone.",
            "mu_M(2^(-1/2))>=mu_8(2^(-1/2))"
            "=22/15-sqrt(2)>0 for M>=8.",
            "Lengths two through seven are outside this theorem.",
        ),
        GateRow(
            "gppm_13_complete_short_guard",
            "exact recrossing guard",
            "proved",
            "A complete two-level chain can still have four simple crossings.",
            "For a>1/2 and alpha=pi/2, X=-sin(theta)-a*sin(2theta)"
            "=-sin(theta)*(1+2a*cos(theta)).",
            "This applies to p=2 and p=3; it is not an Xi phase-attainment claim.",
        ),
        GateRow(
            "gppm_14_guard_count",
            "first-jet guard",
            "proved",
            "The complete short-chain guard has two upward crossings.",
            "zeros: 0, pi, arccos(-1/(2a)), "
            "2*pi-arccos(-1/(2a)); wind(X+i*X_theta)=-2.",
            "The coefficient block itself remains nonzero and winds one.",
        ),
        GateRow(
            "gppm_15_completion_verdict",
            "nonpromotion guard",
            "proved",
            "Multiplicative completion alone does not imply first-jet one-turn behavior.",
            "block winding=1 but projected first-jet winding=-2 in the guard.",
            "External phase and phase speed must be controlled quantitatively.",
        ),
        GateRow(
            "gppm_16_perturbation_target",
            "Xi extension",
            "open",
            "Use the positive long-block margin as a C1 perturbation budget.",
            "|Theta_actual'-Theta_geometric'|<mu_M(p^(-1/2)).",
            "The heat quadratic, d_n phases, and moving base phase must be retained.",
        ),
        GateRow(
            "gppm_17_rejoin_guard",
            "nonpromotion guard",
            "proved",
            "Even monotone one-turn blocks can produce a multi-turn sum.",
            "The existing exact two-block model sums two winding-one blocks "
            "to winding three.",
            "P-free blocks and the recurrent endpoint must be rejoined first.",
        ),
        GateRow(
            "gppm_18_handoff",
            "route decision",
            "open",
            "Split future chain analysis into robust long blocks and explicit short blocks.",
            "Prove C1 stability on the monotone families; keep short p=2,3 "
            "blocks, singleton chains, p-free phases, cutoffs, and endpoint "
            "inside the joined ray defect.",
            "No joined Abel gap or successor winding bound is proved.",
        ),
    ]


def build_payload() -> dict:
    payloads = load_sources()
    source_audit = audit_sources(payloads)
    rows = build_rows()
    proof_boundary = (
        "This artifact proves the correction-free geometric block phase "
        "formula, elementary monotonicity margin, uniform one-turn families "
        "p>=5, p=3 with M>=4, and p=2 with M>=8, and an exact complete "
        "two-level recrossing guard for p=2,3. It does not prove C1 stability "
        "under the heat quadratic or d_n corrections, a joined p-free or "
        "endpoint theorem, the Xi Abel gap, H_j<3*pi/2, contact exclusion, "
        "Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion."
    )
    return {
        "kind": STEM,
        "date": "2026-07-29",
        "status": (
            "exact_geometric_phase_monotonicity_thresholds_and_"
            "complete_short_chain_recrossing_guard"
        ),
        "proof_boundary": proof_boundary,
        "source_audit": source_audit,
        "exact": {
            "geometric_block": {
                "definition": (
                    "a=p^(-1/2), "
                    "P_M(z)=z*sum_(k=0)^(M-1)(a*z)^k"
                ),
                "factorization": (
                    "P_M(z)=z*[1-(a*z)^M]/(1-a*z)"
                ),
                "unit_circle": (
                    "P_M is nonzero on |z|=1 and has winding one"
                ),
            },
            "phase_derivative": {
                "kernel": (
                    "K_m(r,theta)=m*(r^2-r*cos(m*theta))/"
                    "(1-2*r*cos(m*theta)+r^2)"
                ),
                "identity": (
                    "Theta_M'=1+K_M(a^M,theta)-K_1(a,theta)"
                ),
                "kernel_bounds": (
                    "-m*r/(1-r)<=K_m(r,theta)<=m*r/(1+r)"
                ),
                "margin": (
                    "Theta_M'>=mu_M(a):="
                    "1/(1+a)-M*a^M/(1-a^M)"
                ),
                "consequence": (
                    "mu_M(a)>0 implies strict phase increase by 2*pi "
                    "and exactly one upward crossing for every constant "
                    "external phase"
                ),
            },
            "thresholds": {
                "tail": "g_M(a)=M*a^M/(1-a^M)",
                "tail_monotonicity": (
                    "a<=M/(M+1) => g_(M+1)(a)<g_M(a)"
                ),
                "large_primes": (
                    "p>=5: M=1 is exact and M>=2 has "
                    "mu_M>=mu_2(5^(-1/2))=(3-sqrt(5))/4>0"
                ),
                "ternary": (
                    "p=3, M>=4: mu_M>=(2-sqrt(3))/2>0"
                ),
                "dyadic": (
                    "p=2, M>=8: mu_M>=22/15-sqrt(2)>0"
                ),
            },
            "short_chain_guard": {
                "domain": "1/2<a<1, in particular p=2 or p=3",
                "projection": (
                    "alpha=pi/2; X_a(theta)=-sin(theta)-a*sin(2theta)"
                    "=-sin(theta)*(1+2a*cos(theta))"
                ),
                "zeros": (
                    "0, pi, theta_a=arccos(-1/(2a)), 2*pi-theta_a"
                ),
                "derivatives": (
                    "X_a'(0)=-(1+2a)<0, X_a'(pi)=1-2a<0, "
                    "X_a'(theta_a)=X_a'(2*pi-theta_a)="
                    "2a*sin(theta_a)^2>0"
                ),
                "counts": (
                    "two upward and two downward crossings on [0,2*pi)"
                ),
                "first_jet_winding": "wind(X_a+i*X_a')=-2",
                "block_winding": "wind(P_2)=1",
            },
            "route": {
                "long_block_target": (
                    "|Theta_actual'-Theta_geometric'|"
                    "<mu_M(p^(-1/2))"
                ),
                "short_blocks": (
                    "Keep p=2 lengths 2..7 and p=3 lengths 2..3 "
                    "inside an explicit finite block ledger; M=1 is exact"
                ),
                "rejoin": (
                    "Rejoin all p-free bases, singleton chains, the recurrent "
                    "endpoint, and cutoff homotopies before estimating the "
                    "ray defect or successor winding"
                ),
                "open_handoff": (
                    "Derive a C1 heat-and-correction perturbation bound for "
                    "the positive-margin long families. Do not sum blockwise "
                    "winding bounds: insert the long/short decomposition into "
                    "mathsf_X*(L*D_j mathcal_C_N-mathcal_C_N), preserving "
                    "external phases and the recurrent endpoint."
                ),
            },
            "nonpromotion": {
                "completion": (
                    "A complete contracted prime-power chain can have "
                    "projected first-jet winding -2 although the complex "
                    "coefficient block is zero-free with winding one"
                ),
                "composition": (
                    "The sum of individually winding-one blocks can have "
                    "winding three"
                ),
            },
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_phase_derivative_identities": 1,
            "uniform_monotone_prime_families": 3,
            "complete_short_chain_recrossing_guards": 1,
            "nonpromotion_guards": 2,
            "open_C1_perturbation_targets": 1,
            "proved_joined_abel_gaps": 0,
            "proved_successor_winding_bounds": 0,
        },
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    lines = [
        "# Geometric Prime-Power Phase-Monotonicity Gate",
        "",
        "Date: 2026-07-29",
        "",
        "Status: exact correction-free one-block thresholds and complete",
        "short-chain recrossing guard. This is not a proof of the Xi Abel",
        "gap, a successor winding bound, or RH.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"python work/rh_compute/scripts/{STEM}.py",
        f"python work/rh_compute/scripts/check_{STEM}.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        (
            "validated geometric prime-power phase-monotonicity gate: "
            "18 rows, 1 phase-derivative identity, 3 uniform monotone "
            "prime families, 1 complete short-chain recrossing guard, "
            "2 nonpromotion guards, 1 open C1 perturbation target, "
            "0 joined Abel gaps, 0 successor winding bounds"
        ),
        "```",
        "",
        "## Exact Geometric Block",
        "",
        "For the correction-free `t=0` internal chain put",
        "",
        "```text",
        exact["geometric_block"]["definition"],
        exact["geometric_block"]["factorization"],
        exact["geometric_block"]["unit_circle"],
        "```",
        "",
        "The block is zero-free on the unit circle, but zero-free winding",
        "does not by itself control phase backtracking.",
        "",
        "## Phase-Speed Margin",
        "",
        "With `Theta_M(theta)=arg P_M(exp(i*theta))`,",
        "",
        "```text",
        exact["phase_derivative"]["kernel"],
        exact["phase_derivative"]["identity"],
        exact["phase_derivative"]["kernel_bounds"],
        exact["phase_derivative"]["margin"],
        "```",
        "",
        exact["phase_derivative"]["consequence"] + ".",
        "",
        "The tail penalty decreases by",
        "",
        "```text",
        exact["thresholds"]["tail_monotonicity"],
        "```",
        "",
        "so three uniform families follow:",
        "",
        "```text",
        exact["thresholds"]["large_primes"],
        exact["thresholds"]["ternary"],
        exact["thresholds"]["dyadic"],
        "```",
        "",
        "These are exact correction-free block theorems.",
        "",
        "## Complete Short-Chain Guard",
        "",
        "For `a>1/2`, which includes `p=2,3`, choose the constant external",
        "phase `alpha=pi/2`. Then",
        "",
        "```text",
        exact["short_chain_guard"]["projection"],
        exact["short_chain_guard"]["zeros"],
        exact["short_chain_guard"]["derivatives"],
        exact["short_chain_guard"]["counts"],
        exact["short_chain_guard"]["first_jet_winding"],
        exact["short_chain_guard"]["block_winding"],
        "```",
        "",
        "Thus a complete two-level prime-power block may be complex-zero-free",
        "with winding one while its real first jet winds `-2`. This is a",
        "generic external-phase guard, not a claim that the Xi trajectory",
        "attains that phase.",
        "",
        "## Surviving Route",
        "",
        "The positive margins support the open perturbation target",
        "",
        "```text",
        exact["route"]["long_block_target"],
        "```",
        "",
        "with the heat quadratic, coefficient corrections, and moving base",
        "phase retained. Short dyadic and ternary blocks remain explicit:",
        "",
        "```text",
        exact["route"]["short_blocks"],
        "```",
        "",
        exact["route"]["rejoin"] + ".",
        "",
        "Open handoff:",
        "",
        exact["route"]["open_handoff"],
        "",
        "The `pi` in the phase period is the standard period of",
        "`exp(i*theta)`. The prime-power chain introduces no new value of",
        "`pi`; the Xi saddle continues to use the completed-zeta constant.",
        "",
        "## Boundary",
        "",
        payload["proof_boundary"],
        "",
    ]
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "built geometric prime-power phase-monotonicity gate: "
        "18 rows, 1 phase-derivative identity, 3 uniform monotone "
        "prime families, 1 complete short-chain recrossing guard, "
        "2 nonpromotion guards, 1 open C1 perturbation target, "
        "0 joined Abel gaps, 0 successor winding bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
