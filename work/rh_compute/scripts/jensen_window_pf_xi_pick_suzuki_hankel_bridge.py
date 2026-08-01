#!/usr/bin/env python3
"""Build the exact Xi Pick/Suzuki arithmetic-Hankel bridge and guards."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_xi_pick_suzuki_hankel_bridge.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_xi_pick_suzuki_hankel_bridge.md"
)


@dataclass(frozen=True)
class BridgeRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def polynomial_guard() -> dict:
    p = Fraction(1)
    q = Fraction(2)
    a = Fraction(11, 10)
    b = Fraction(7, 5)

    alpha_real = p * p - q * q
    alpha_imag = 2 * p * q
    z_real = a * a - b * b
    z_imag = 2 * a * b
    disk_defect = (
        (z_real - alpha_real) ** 2 + z_imag**2 - alpha_imag**2
    )
    lower_denominator = (
        (z_real - alpha_real) ** 2 + (z_imag - alpha_imag) ** 2
    )
    upper_denominator = (
        (z_real - alpha_real) ** 2 + (z_imag + alpha_imag) ** 2
    )
    anti_pick = (
        (z_imag - alpha_imag) / lower_denominator
        + (z_imag + alpha_imag) / upper_denominator
    )

    roots = (
        (Fraction(1), Fraction(2)),
        (Fraction(1), Fraction(-2)),
        (Fraction(-1), Fraction(2)),
        (Fraction(-1), Fraction(-2)),
    )
    horizontal_log_derivative = sum(
        (a - root_real)
        / ((a - root_real) ** 2 + (b - root_imag) ** 2)
        for root_real, root_imag in roots
    )

    expected = {
        "alpha_real": Fraction(-3),
        "alpha_imag": Fraction(4),
        "z_real": Fraction(-3, 4),
        "z_imag": Fraction(77, 25),
        "disk_defect": Fraction(-14511, 10000),
        "lower_denominator": Fraction(59089, 10000),
        "upper_denominator": Fraction(551889, 10000),
    }
    observed = {
        "alpha_real": alpha_real,
        "alpha_imag": alpha_imag,
        "z_real": z_real,
        "z_imag": z_imag,
        "disk_defect": disk_defect,
        "lower_denominator": lower_denominator,
        "upper_denominator": upper_denominator,
    }
    if observed != expected:
        raise RuntimeError(f"unexpected polynomial witness: {observed}")
    if horizontal_log_derivative <= 0:
        raise RuntimeError("horizontal log derivative is not positive")
    if anti_pick >= 0:
        raise RuntimeError("polynomial witness did not violate the Pick sign")

    return {
        "F_star": "F_*(z)=z^2+6*z+25",
        "M_star": "M_*(w)=F_*(w^2)=w^4+6*w^2+25",
        "witness_w": "11/10+(7/5)*i",
        "witness_z": f"{fraction_text(z_real)}+({fraction_text(z_imag)})*i",
        "alpha": f"{fraction_text(alpha_real)}+({fraction_text(alpha_imag)})*i",
        "disk_defect": fraction_text(disk_defect),
        "lower_denominator": fraction_text(lower_denominator),
        "upper_denominator": fraction_text(upper_denominator),
        "horizontal_log_derivative": fraction_text(horizontal_log_derivative),
        "anti_pick_log_derivative": fraction_text(anti_pick),
        "conclusion": (
            "Horizontal modulus growth in Re(w)>1 does not imply the "
            "squared-variable Pick sign."
        ),
    }


def exact_statements() -> dict:
    return {
        "normalization": (
            "H_0(x)=xi((1+i*x)/2)/8; "
            "M(w)=integral_R Phi(u)cosh(u*w)du=2*H_0(i*w)="
            "xi((1+w)/2)/4; F(z)=M(sqrt(z))"
        ),
        "xi_directional_pick": (
            "For s=sigma+i*tau, delta=sigma-1/2>0, "
            "w=2*(delta+i*tau), and xi(s)!=0: "
            "2*|w|^2*P_Phi(w^2)/|M(w)|^2="
            "tau*Re(xi'(s)/xi(s))-delta*Im(xi'(s)/xi(s))"
        ),
        "hyperbolic_monotonicity": (
            "For X=delta^2-tau^2 and Y=2*delta*tau: "
            "partial_Y log|xi(1/2+delta+i*tau)| at fixed X="
            "[tau*Re(xi'/xi)-delta*Im(xi'/xi)]/"
            "[2*(delta^2+tau^2)]"
        ),
        "horizontal_growth_theorem": (
            "Sondow-Dumitrescu: Re(xi'(s)/xi(s))>0 in every open "
            "right half-plane lying strictly to the right of all xi zeros; "
            "in particular sigma>1"
        ),
        "zero_disk_identity": (
            "If F(z)=(z-(A+i*B))*(z-(A-i*B)) and z=x+i*y, y>0, then "
            "-Im(F'(z)/F(z))="
            "2*y*((x-A)^2+y^2-B^2)/"
            "([((x-A)^2+(y-B)^2)]*[((x-A)^2+(y+B)^2)])"
        ),
        "horizontal_countermodel": (
            "F_*(z)=z^2+6*z+25 and M_*(w)=F_*(w^2) have "
            "zeros +/-1+/-2*i; |M_*(w)| is strictly increasing in "
            "Re(w) for Re(w)>1, but P_*(w^2)<0 at "
            "w=11/10+(7/5)*i"
        ),
        "suzuki_family": (
            "E_(omega,nu)(r)=xi(1/2+omega-i*r)^nu and "
            "Theta_(omega,nu)(r)=E#(r)/E(r)="
            "[xi(1/2-omega-i*r)/xi(1/2+omega-i*r)]^nu"
        ),
        "arithmetic_hankel_operator": (
            "K_(omega,nu)=Fourier_inverse(Theta_(omega,nu)); "
            "(K_(omega,nu)[t]f)(x)=1_(x<t)*integral_(-infinity)^t "
            "K_(omega,nu)(x+y)f(y)dy"
        ),
        "local_hamiltonian": (
            "For nu*omega>1, Suzuki constructs tau>0 with "
            "det(I+/-K[t])!=0 on 0<=t<tau and "
            "H(t)=diag(1/gamma(t),gamma(t)), "
            "gamma(t)=[det(I+K[t])/det(I-K[t])]^2"
        ),
        "boundary_unitarity_guard": (
            "Theta_*(r)=(r+i)/(r-i) has |Theta_*(u)|=1 for every real u "
            "but has a pole at r=i; moving a high-contour inverse Fourier "
            "integral to the real axis acquires a nonzero residue proportional "
            "to exp(x). Boundary unimodularity alone does not make the "
            "high-contour Hankel operator unitary or contractive"
        ),
        "suzuki_equivalence": (
            "RH iff there exist omega_n decreasing to 0 and integers nu_n with "
            "nu_n*omega_n>1 such that det(I+/-K_(omega_n,nu_n)[t])!=0 "
            "for every t>=0 and J_(omega_n,nu_n)(t;r,r)->0 as t->infinity "
            "for every r in the upper half-plane"
        ),
        "fredholm_hankel_series": (
            "det(I+/-K[t])=sum_(m>=0)(+/-1)^m/m!*"
            "integral_((-infinity,t)^m) "
            "det(K(x_i+x_j))_(i,j=1)^m dx_1...dx_m"
        ),
        "rank_one_tn_guard": (
            "The constant Hankel kernel K(x+y)=1 on L2(0,1) is totally "
            "nonnegative and rank one, but its nonzero eigenvalue is 1, so "
            "det(I-K)=0 and det(I+K)=2"
        ),
        "spectral_target": (
            "The finite-t spectral target is pointwise-in-t strict "
            "contractivity ||K_(omega_n,nu_n)[t]||<1 for every finite t on "
            "a cofinal omega_n->0 sequence; no positive gap uniform in t is "
            "required or possible in the desired HB case"
        ),
        "determinant_only_reduction": (
            "A strictly decreasing cofinal omega_n->0 family with "
            "det(I+/-K_(omega_n,nu_n)[t])!=0 for every n and t already "
            "implies RH: the compatible finite forms extend to a bounded "
            "causal multiplier, and any off-line zero would generate "
            "distinct shifted zeros accumulating at itself"
        ),
        "terminal_reduction": (
            "For one fixed pair the terminal limit is not a direct "
            "consequence of contractivity. Within the cofinal RH criterion, "
            "however, determinant nonvanishing first implies RH, after which "
            "Suzuki's Theorem 2.3 supplies J(t;r,r)->0"
        ),
    }


def build_payload() -> dict:
    exact = exact_statements()
    rows = [
        BridgeRow(
            "xpsh_01_phi_xi_normalization",
            "exact_identity",
            "available_exact",
            "The Phi transform is exactly the completed xi function in the square-root coordinate.",
            exact["normalization"],
            "The constant factor does not affect logarithmic derivatives or zero locations.",
        ),
        BridgeRow(
            "xpsh_02_xi_directional_pick",
            "exact_identity",
            "available_exact",
            "The global Pick sign is a directional logarithmic-derivative inequality for xi.",
            exact["xi_directional_pick"],
            "This identity is evaluated away from zeros; strict global positivity itself excludes them.",
        ),
        BridgeRow(
            "xpsh_03_hyperbolic_modulus_monotonicity",
            "exact_equivalence",
            "ready_to_apply",
            "The Pick target is monotonicity of |xi| along the first-quadrant hyperbolas X=constant.",
            exact["hyperbolic_monotonicity"],
            "It is endpoint-equivalent and is not proved here.",
        ),
        BridgeRow(
            "xpsh_04_sondow_horizontal_growth",
            "published_exact_theorem",
            "available_exact",
            "Classical horizontal xi-modulus growth controls only Re(xi'/xi).",
            exact["horizontal_growth_theorem"],
            "It does not control the second directional term -delta*Im(xi'/xi).",
        ),
        BridgeRow(
            "xpsh_05_zero_disk_identity",
            "exact_identity",
            "available_exact",
            "A conjugate off-axis zero pair contributes with a sharp forbidden-disk sign.",
            exact["zero_disk_identity"],
            "Using actual xi zeros in this formula is diagnostic, not a noncircular proof.",
        ),
        BridgeRow(
            "xpsh_06_horizontal_growth_countermodel",
            "exact_countermodel",
            "guard_validated",
            "Horizontal modulus growth strictly to the right of every zero does not imply the Pick sign.",
            exact["horizontal_countermodel"],
            "The witness is a finite polynomial, not xi.",
        ),
        BridgeRow(
            "xpsh_07_suzuki_hb_family",
            "published_exact_definition",
            "available_exact",
            "Suzuki's shifted completed-xi family gives an arithmetic Hermite-Biehler route.",
            exact["suzuki_family"],
            "Membership in the Hermite-Biehler class for all omega>0 is equivalent to RH, not assumed here.",
        ),
        BridgeRow(
            "xpsh_08_arithmetic_hankel_operator",
            "published_exact_construction",
            "available_exact",
            "The inverse Fourier ratio produces a self-adjoint arithmetic Hankel operator.",
            exact["arithmetic_hankel_operator"],
            "Its kernel includes both gamma-factor and arithmetic Dirichlet-coefficient data.",
        ),
        BridgeRow(
            "xpsh_09_local_positive_hamiltonian",
            "published_exact_theorem",
            "available_exact",
            "For nu*omega>1 the operator and positive diagonal Hamiltonian exist unconditionally on a nonzero initial interval.",
            exact["local_hamiltonian"],
            "Local existence does not prove extension to every t or the terminal condition.",
        ),
        BridgeRow(
            "xpsh_16_boundary_unimodularity_guard",
            "exact_countermodel",
            "guard_validated",
            "Boundary modulus one does not justify replacing Suzuki's high-contour kernel by a unitary real-boundary multiplier.",
            exact["boundary_unitarity_guard"],
            "Contour deformation is valid only after excluding intervening poles; using innerness here would be circular.",
        ),
        BridgeRow(
            "xpsh_10_suzuki_global_equivalence",
            "published_exact_equivalence",
            "ready_to_apply",
            "Suzuki reduces RH to global Fredholm nonvanishing plus terminal canonical-kernel collapse on a cofinal omega sequence.",
            exact["suzuki_equivalence"],
            "This is a precise equivalent target, not a proof that either global condition holds.",
        ),
        BridgeRow(
            "xpsh_11_fredholm_hankel_hierarchy",
            "exact_identity",
            "available_exact",
            "The global determinant gate is an integrated all-order Hankel-determinant hierarchy.",
            exact["fredholm_hankel_series"],
            "Termwise determinant signs alone do not control the alternating I-K series.",
        ),
        BridgeRow(
            "xpsh_12_rank_one_tn_spectral_guard",
            "exact_countermodel",
            "guard_validated",
            "Total nonnegativity of the Hankel kernel does not exclude eigenvalue +1.",
            exact["rank_one_tn_guard"],
            "A strict contraction or another quantitative exclusion of +/-1 is still required.",
        ),
        BridgeRow(
            "xpsh_13_strict_spectral_target",
            "sufficient_theorem_target",
            "not_ready_to_apply",
            "Pointwise strict contractivity at every finite t would close Suzuki's Fredholm nonvanishing gate.",
            exact["spectral_target"],
            "No such all-t theorem is proved for the arithmetic kernels as omega tends to zero.",
        ),
        BridgeRow(
            "xpsh_14_determinant_only_reduction",
            "exact_cross_artifact_reduction",
            "internally_audited",
            "On a cofinal family, all-time determinant nonvanishing already implies RH and makes the terminal premise redundant.",
            exact["determinant_only_reduction"],
            exact["terminal_reduction"],
        ),
        BridgeRow(
            "xpsh_15_xi_pick_suzuki_handoff",
            "open_structural_handoff",
            "not_ready_to_apply",
            "The abstract positive-resolvent idea is replaced by one explicit arithmetic Hankel spectral programme.",
            (
                "Prove xpsh_13 on one cofinal omega_n->0 "
                "sequence without assuming HB, zero reality, RH, or Lambda<=0"
            ),
            "The Pick endpoint, Suzuki criterion, and signed-Hankel evidence remain distinct until the all-time spectral gate is proved.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_xi_pick_suzuki_hankel_bridge",
        "date": "2026-07-23",
        "status": (
            "exact Phi-to-xi directional Pick reduction, published Suzuki "
            "arithmetic-Hankel equivalence, and exact theorem-mismatch guards"
        ),
        "exact": exact,
        "polynomial_guard": polynomial_guard(),
        "rows": [asdict(row) for row in rows],
        "sources": [
            {
                "title": (
                    "Masatoshi Suzuki, Hamiltonians arising from L-functions "
                    "in the Selberg class, JFA 281 (2021), 109116"
                ),
                "url": "https://doi.org/10.1016/j.jfa.2021.109116",
                "relevant_result": "Theorems 2.1-2.4",
            },
            {
                "title": (
                    "Jonathan Sondow and Cristian Dumitrescu, A monotonicity "
                    "property of Riemann's xi function"
                ),
                "url": "https://arxiv.org/abs/1005.1104",
                "relevant_result": "Theorem 1",
            },
            {
                "title": "Local Phi Pick-kernel endpoint equivalence",
                "path": "outputs/jensen_window_pf_phi_pick_kernel_target.md",
            },
            {
                "title": "Local exact normalization and endpoint ledger",
                "path": "outputs/formal_core.md",
            },
            {
                "title": "Suzuki determinant-only reduction",
                "path": (
                    "outputs/"
                    "jensen_window_pf_suzuki_determinant_only_reduction.md"
                ),
            },
        ],
        "proof_boundary": (
            "This artifact proves the exact Phi-to-xi normalization, the "
            "directional and hyperbolic forms of the Pick target, and three exact "
            "countermodel guards. It records Suzuki's published unconditional "
            "local construction and exact global equivalence. It also records "
            "the internally audited determinant-only cofinal reduction. It "
            "does not prove global Fredholm nonvanishing, the Phi Pick sign, "
            "LP+, PF-infinity, Jensen hyperbolicity for zeta, RH, or Lambda<=0."
        ),
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    guard = payload["polynomial_guard"]
    return "\n".join(
        [
            "# Jensen-Window PF Xi Pick/Suzuki Hankel Bridge",
            "",
            "Date: 2026-07-23",
            "",
            "Status: exact Phi-to-xi directional Pick reduction, published Suzuki",
            "arithmetic-Hankel equivalence, and exact theorem-mismatch guards.",
            "This is not a proof of LP+, PF-infinity, RH, or `Lambda <= 0`.",
            "",
            "Artifact kind: `jensen_window_pf_xi_pick_suzuki_hankel_bridge`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_xi_pick_suzuki_hankel_bridge.json",
            "python work/rh_compute/scripts/jensen_window_pf_xi_pick_suzuki_hankel_bridge.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_xi_pick_suzuki_hankel_bridge.py",
            "```",
            "",
            "## Exact Xi Coordinate",
            "",
            "The workspace normalization gives",
            "",
            "```text",
            exact["normalization"],
            "```",
            "",
            "Put `s=sigma+i*tau`, `delta=sigma-1/2>0`, and",
            "`w=2*(delta+i*tau)`. Away from zeros, the Pick numerator becomes",
            "",
            "```text",
            exact["xi_directional_pick"],
            "```",
            "",
            "For `X=delta^2-tau^2` and `Y=2*delta*tau`, the same quantity is",
            "",
            "```text",
            exact["hyperbolic_monotonicity"],
            "```",
            "",
            "Thus the Pick endpoint is strict increase of `|xi|` along every",
            "first-quadrant branch `X=constant` as `Y` increases. This differs",
            "from ordinary horizontal increase in `sigma`.",
            "",
            "## Horizontal-Growth Guard",
            "",
            "The Sondow-Dumitrescu theorem gives",
            "",
            "```text",
            exact["horizontal_growth_theorem"],
            "```",
            "",
            "but the Pick direction also contains `-delta*Im(xi'/xi)`. The",
            "distinction is forced by an exact polynomial.",
            "",
            "For a conjugate zero pair,",
            "",
            "```text",
            exact["zero_disk_identity"],
            "```",
            "",
            "Take",
            "",
            "```text",
            guard["F_star"],
            guard["M_star"],
            "zeros(M_*)={+/-1+/-2*i}",
            f"w={guard['witness_w']}",
            f"z=w^2={guard['witness_z']}",
            f"(x+3)^2+y^2-16={guard['disk_defect']}<0",
            f"partial_Re(w) log|M_*(w)|={guard['horizontal_log_derivative']}>0",
            f"-Im(F_*'(z)/F_*(z))={guard['anti_pick_log_derivative']}<0",
            "```",
            "",
            "Every zero of `M_*` has real part at most one, so the horizontal",
            "logarithmic derivative is positive throughout `Re(w)>1`. The Pick",
            "sign still fails at the displayed point. Horizontal xi-modulus",
            "monotonicity is therefore theorem-mismatched to the required",
            "hyperbolic direction.",
            "",
            "## Suzuki Arithmetic Hankel Route",
            "",
            "For the completed zeta function define",
            "",
            "```text",
            exact["suzuki_family"],
            exact["arithmetic_hankel_operator"],
            "```",
            "",
            "For zeta, whose Selberg-class degree is one, Suzuki proves",
            "unconditionally that `nu*omega>1` gives an initial interval on",
            "which",
            "",
            "```text",
            exact["local_hamiltonian"],
            "```",
            "",
            "The global theorem is the exact equivalence",
            "",
            "```text",
            exact["suzuki_equivalence"],
            "```",
            "",
            "Suzuki states two explicit obligations for an arithmetic operator",
            "already constructed without RH. The later causal-multiplier audit",
            "sharpens their cofinal logic:",
            "",
            "```text",
            exact["determinant_only_reduction"],
            "```",
            "",
            "The real-boundary identity `|Theta(u)|=1` does not supply a",
            "shortcut to the spectral gate:",
            "",
            "```text",
            exact["boundary_unitarity_guard"],
            "```",
            "",
            "The Suzuki kernel is defined from a high horizontal contour. Moving",
            "that contour to the real line without residue already requires the",
            "relevant pole-free analyticity, so a boundary-unitary argument would",
            "be circular.",
            "",
            "## Signed-Hankel Contact",
            "",
            "The classical Fredholm series is",
            "",
            "```text",
            exact["fredholm_hankel_series"],
            "```",
            "",
            "so each coefficient is an integral of a continuum Hankel determinant.",
            "This is a genuine contact with the signed-Hankel programme. It does",
            "not make total nonnegativity sufficient:",
            "",
            "```text",
            exact["rank_one_tn_guard"],
            "```",
            "",
            "The surviving first gate is therefore quantitative:",
            "",
            "```text",
            exact["spectral_target"],
            "```",
            "",
            "Here `for every finite t` is pointwise in `t`; it does not mean",
            "one epsilon-sized gap valid uniformly as `t` tends to infinity.",
            "The exact truncation-path distinction is developed in",
            "`outputs/jensen_window_pf_suzuki_spectral_frontier.md`.",
            "",
            "For the terminal premise, the exact distinction is",
            "",
            "```text",
            exact["terminal_reduction"],
            "```",
            "",
            "The full proof is in",
            "`outputs/jensen_window_pf_suzuki_determinant_only_reduction.md`.",
            "Finite `t` grids, finite-rank discretizations, local Hamiltonian",
            "positivity, or a canonical system built only after assuming the",
            "Hermite-Biehler property still cannot prove the all-time gate.",
            "",
            "## Sources",
            "",
            "- Masatoshi Suzuki, `Hamiltonians arising from L-functions in the Selberg class`, Theorems 2.1-2.4: https://doi.org/10.1016/j.jfa.2021.109116",
            "- Jonathan Sondow and Cristian Dumitrescu, `A monotonicity property of Riemann's xi function and a reformulation of the Riemann Hypothesis`: https://arxiv.org/abs/1005.1104",
            "- `outputs/jensen_window_pf_phi_pick_kernel_target.md`",
            "- `outputs/formal_core.md`",
            "- `outputs/jensen_window_pf_suzuki_determinant_only_reduction.md`",
            "",
            "## Proof Boundary",
            "",
            payload["proof_boundary"],
            "",
        ]
    )


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
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Jensen-window PF Xi Pick/Suzuki Hankel bridge: "
        "16 rows, 7 exact coordinate/guard identities, 4 published "
        "operator steps, 3 exact countermodels, 1 open global gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
