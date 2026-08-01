#!/usr/bin/env python3
"""Build the prime-power heat-starlikeness base certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path

from prime_power_heat_bernstein import build_base_certificates


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_centered_prime_power_heat_"
    "starlikeness_base_certificate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "phase_flow": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_prime_power_"
        "logarithmic_phase_flow_gate.json"
    ),
    "complete_chain": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_complete_prime_power_"
        "chain_occupation_gate.json"
    ),
    "geometric_phase": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_geometric_prime_power_"
        "phase_monotonicity_gate.json"
    ),
    "direct_projection": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_direct_projection_"
        "regime_reduction.json"
    ),
    "normalized_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_normalized_prefix_"
        "phase_flux_reduction.json"
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
    phase = payloads["phase_flow"].get("exact", {})
    if "J_heat=|Q_heat|^2" not in phase.get(
        "replacement_target", {}
    ).get("current", ""):
        raise RuntimeError("phase-flow replacement target drifted")

    chain = payloads["complete_chain"].get("chain", {})
    if "Q_(m,p)(w)=sum_(k=0)^R_m" not in chain.get(
        "chain_polynomial", ""
    ):
        raise RuntimeError("complete-chain polynomial drifted")

    geometric = payloads["geometric_phase"].get("exact", {})
    if "22/15-sqrt(2)" not in geometric.get("thresholds", {}).get(
        "dyadic", ""
    ):
        raise RuntimeError("geometric dyadic threshold drifted")

    projection = payloads["direct_projection"].get("exact", {})
    if "0<=delta_a<1/(4x)" not in projection.get(
        "saddle_bounds", ""
    ):
        raise RuntimeError("saddle offset bound drifted")

    prefix = payloads["normalized_prefix"].get("exact", {})
    if "|epsilon_n|<8446/x^2" not in prefix.get(
        "q_ge_1_inward", ""
    ):
        raise RuntimeError("coefficient current bound drifted")

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
            "pphsb_01_domain",
            "fixed-chart domain",
            "proved",
            "Work on one complete prime-power chain in the certified ray regime.",
            "L>=50, 0<=t<=1/2, 0<=tL<=25; "
            "n_k=m*p^k, 0<=k<=M-1.",
            "The certificate concerns only complete long-block base lengths.",
        ),
        GateRow(
            "pphsb_02_ideal_coordinates",
            "exact heat coordinates",
            "proved",
            "Removing delta_a and d_n leaves a three-variable positive coefficient family.",
            "c_k=p^(-k/2)*q^[k(2M-2-k)]*s^k; "
            "q=exp(-t(log p)^2/4), s=exp(-t(log p)u/2).",
            "The omitted factors are restored quantitatively below.",
        ),
        GateRow(
            "pphsb_03_current",
            "exact angular current",
            "proved",
            "The shifted phase is increasing exactly when its division-free current is positive.",
            "Q=sum_k c_k z^k; H=sum_k k*c_k*z^k; "
            "J_0=|Q|^2+Re(H*conj(Q)).",
            "No nonvanishing assumption is needed to state J_0>0.",
        ),
        GateRow(
            "pphsb_04_chebyshev",
            "exact polynomialization",
            "proved",
            "The circle current is a real polynomial in x=cos(theta).",
            "J_0=sum_k (k+1)c_k^2"
            "+sum_(k<l)(k+l+2)c_k*c_l*T_(l-k)(x).",
            "T_d is the Chebyshev polynomial defined by cos(d theta).",
        ),
        GateRow(
            "pphsb_05_bernstein",
            "exact positivity method",
            "proved",
            "A positive lower bound for every Bernstein coefficient on a box proves positivity throughout that box.",
            "P=sum_I b_I prod_j B_(I_j,n_j)(y_j) "
            "and B_(i,n)>=0, sum_i B_(i,n)=1.",
            "All affine substitutions and de Casteljau splits use Fractions.",
        ),
        GateRow(
            "pphsb_06_radicals",
            "exact algebraic enclosure",
            "proved",
            "Every coefficient is enclosed as r+v*sqrt(p) with an outward rational square-root interval.",
            "sqrt(p) in [floor(10^90 sqrt(p))/10^90, "
            "(floor(10^90 sqrt(p))+1)/10^90].",
            "The sign of v selects the rigorous lower endpoint.",
        ),
        GateRow(
            "pphsb_07_dyadic_box",
            "parameter containment",
            "proved",
            "Every ideal dyadic base block lies in a rational certificate box.",
            "(p,M)=(2,8): q in [47/50,1], "
            "s in [22/25,1], x in [-1,1].",
            "The bounds use log(2)<7/10 and N>=2^25.",
        ),
        GateRow(
            "pphsb_08_dyadic_margin",
            "exact Bernstein certificate",
            "proved",
            "The ideal dyadic base current is uniformly positive.",
            "J_0>1/300 for (p,M)=(2,8).",
            "Six dyadic x-leaves cover [-1,1] exactly.",
        ),
        GateRow(
            "pphsb_09_ternary_box",
            "parameter containment",
            "proved",
            "Every ideal ternary base block lies in a rational certificate box.",
            "(p,M)=(3,4): q in [17/20,1], "
            "s in [73/100,1], x in [-1,1].",
            "The bounds use log(3)<11/10 and N>=2^25.",
        ),
        GateRow(
            "pphsb_10_ternary_margin",
            "exact Bernstein certificate",
            "proved",
            "The ideal ternary base current is uniformly positive.",
            "J_0>1/25 for (p,M)=(3,4).",
            "Two half-interval x-leaves cover [-1,1] exactly.",
        ),
        GateRow(
            "pphsb_11_large_prime",
            "analytic two-level family",
            "proved",
            "Every ideal two-level block with p>=5 has one uniform current margin.",
            "Q=1+r*z, 0<=r<=1/sqrt(5); "
            "min_theta J_0=(1-r)(1-2r)>1/20.",
            "This is an analytic family theorem, not a p=5-only sample.",
        ),
        GateRow(
            "pphsb_12_quinary_audit",
            "exact Bernstein cross-check",
            "proved",
            "The worst-prime p=5 two-level box independently reproduces the analytic margin.",
            "(p,M)=(5,2): q in [18/25,1], "
            "s in [13/25,1], J_0>1/20.",
            "The analytic p>=5 argument is the promoted result.",
        ),
        GateRow(
            "pphsb_13_offset",
            "saddle-offset transfer",
            "proved",
            "The positive saddle offset is isolated as a multiplicative perturbation.",
            "exp[t*k*log(p)*delta_a/2], "
            "0<=delta_a<1/(4x).",
            "It is not hidden inside the lower endpoint for s.",
        ),
        GateRow(
            "pphsb_14_coefficient_error",
            "actual coefficient transfer",
            "proved",
            "Every normalized actual coefficient differs relatively from its ideal coefficient by less than r_x.",
            "|R_k-1|<r_x=18000/x; "
            "R_k=exp[t*k*log(p)*delta_a/2]"
            "*(1+d_(m*p^k))/(1+d_m).",
            "The estimate covers M=8, M=4, and all p>=5 with M=2.",
        ),
        GateRow(
            "pphsb_15_current_error",
            "current stability",
            "proved",
            "The normalized frozen current changes by less than 900*r_x.",
            "|J_ang-J_0|<576*r_x+288*r_x^2<900*r_x.",
            "Only M<=8 and 0<c_k<=1 are used.",
        ),
        GateRow(
            "pphsb_16_ray_error",
            "fixed-ray transfer",
            "proved",
            "Radial drift and correction current keep the total ray defect below one thousandth.",
            "|J_ray-J_0|<16200448/x+6486528/x^2<1/1000.",
            "Here x>12*2^50, rho<1/(2x), and |epsilon_n|<8446/x^2.",
        ),
        GateRow(
            "pphsb_17_actual_margins",
            "actual ray positivity",
            "proved",
            "All three base families retain strict fixed-ray phase-current positivity.",
            "p=2,M=8: J_ray>7/3000; "
            "p=3,M=4: J_ray>39/1000; "
            "p>=5,M=2: J_ray>49/1000.",
            "The currents are normalized by the nonzero k=0 coefficient.",
        ),
        GateRow(
            "pphsb_18_boundary",
            "proof boundary",
            "open",
            "Base-family positivity does not yet cover longer complete chains or joined p-free bases.",
            "Open: (2,M>8), (3,M>4), (p>=5,M>2), "
            "short/singleton chains, endpoint joins, and the Xi Abel gap.",
            "No RH conclusion or successor winding bound is claimed.",
        ),
    ]


def build_payload() -> dict:
    source_payloads = load_sources()
    source_audit = audit_sources(source_payloads)
    rows = build_rows()
    certificates = [
        certificate.to_dict()
        for certificate in build_base_certificates()
    ]
    return {
        "kind": STEM,
        "date": "2026-07-30",
        "status": (
            "proved base-family correction-free Bernstein positivity "
            "and actual fixed-ray perturbation transfer"
        ),
        "proof_boundary": (
            "This certificate proves strict actual fixed-ray phase-current "
            "positivity only for the threshold base families "
            "(p,M)=(2,8), (3,4), and p>=5 with M=2. It does not prove "
            "length propagation, short or singleton chain control, joined "
            "p-free base compatibility, recurrent-endpoint composition, "
            "a Xi Abel gap, a successor winding bound, or RH."
        ),
        "source_audit": source_audit,
        "exact": {
            "ideal_family": {
                "coefficients": (
                    "c_k=p^(-k/2)*q^[k*(2M-2-k)]*s^k"
                ),
                "parameters": (
                    "q=exp(-t*(log p)^2/4); "
                    "s=exp(-t*(log p)*u/2); "
                    "0<=u<log(p)+1/N"
                ),
                "polynomial": "Q(z)=sum_(k=0)^(M-1)c_k*z^k",
                "moment": "H(z)=sum_(k=0)^(M-1)k*c_k*z^k",
                "current": "J_0=|Q|^2+Re(H*conj(Q))",
                "chebyshev_expansion": (
                    "J_0=sum_k(k+1)c_k^2"
                    "+sum_(k<l)(k+l+2)c_k*c_l*T_(l-k)(x), "
                    "x=cos(theta)"
                ),
            },
            "bernstein_method": {
                "principle": (
                    "On [0,1]^d, every tensor Bernstein basis function "
                    "is nonnegative and the basis sums to one; therefore "
                    "P>=min_I b_I."
                ),
                "arithmetic": (
                    "Exact Fraction affine transforms, exact Fraction "
                    "power-to-Bernstein conversion, exact de Casteljau "
                    "x-splits, and 90-digit outward rational sqrt bounds."
                ),
                "certificates": certificates,
            },
            "parameter_boxes": {
                "dyadic": (
                    "log(2)<7/10; q>=exp[-(log 2)^2/8]>47/50; "
                    "s>=exp[-log(2)(log(2)+1/N)/4]>22/25"
                ),
                "ternary": (
                    "log(3)<11/10; q>=exp[-(log 3)^2/8]>17/20; "
                    "s>=exp[-log(3)(log(3)+1/N)/4]>73/100"
                ),
                "quinary_audit": (
                    "log(5)<161/100; q>=exp[-(log 5)^2/8]>18/25; "
                    "s>=exp[-log(5)(log(5)+1/N)/4]>13/25"
                ),
                "cutoff": (
                    "a^2=exp(L)+t/16, N=floor(a), and L>=50 "
                    "imply N>=2^25"
                ),
            },
            "large_prime_family": {
                "reduction": (
                    "For M=2, Q=1+r*z with "
                    "0<=r=p^(-1/2)*q*s<=1/sqrt(5)"
                ),
                "current": (
                    "J_0=1+2r^2+3r*cos(theta)"
                ),
                "minimum": (
                    "min_theta J_0=(1-r)(1-2r)"
                    ">=7/5-3/sqrt(5)>1/20"
                ),
            },
            "actual_transfer": {
                "normalized_factor": (
                    "R_k=exp[t*k*log(p)*delta_a/2]"
                    "*(1+d_(m*p^k))/(1+d_m)"
                ),
                "coefficient_error": (
                    "|R_k-1|<r_x=18000/x"
                ),
                "frozen_current_error": (
                    "|J_ang-J_0|<576*r_x+288*r_x^2<900*r_x"
                ),
                "ray_terms": (
                    "J_ray=J_ang-rho*Im(H*conj(Q))"
                    "+Im(E*conj(Q))/(log(p)*(-b))"
                ),
                "ray_error": (
                    "|J_ray-J_0|"
                    "<16200448/x+6486528/x^2<1/1000"
                ),
                "margins": {
                    "dyadic": "J_ray>1/300-1/1000=7/3000",
                    "ternary": "J_ray>1/25-1/1000=39/1000",
                    "large_prime": "J_ray>1/20-1/1000=49/1000",
                },
            },
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_bernstein_base_certificates": 3,
            "analytic_large_prime_base_families": 1,
            "actual_ray_positive_base_families": 3,
            "length_propagation_theorems": 0,
            "proved_joined_abel_gaps": 0,
            "proved_successor_winding_bounds": 0,
        },
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    certificates = exact["bernstein_method"]["certificates"]
    lines = [
        "# Prime-Power Heat-Starlikeness Base Certificate",
        "",
        "Date: 2026-07-30",
        "",
        "Status: exact correction-free base-family positivity and a",
        "rigorous transfer to the actual fixed-ray current. This certificate",
        "does not promote RH or Lambda<=0.",
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
            "validated prime-power heat-starlikeness base certificate: "
            "18 rows, 3 exact Bernstein base certificates, "
            "1 analytic p>=5 two-level family, "
            "3 actual ray-positive base families, "
            "0 length-propagation theorems, 0 joined Abel gaps, "
            "0 successor winding bounds"
        ),
        "```",
        "",
        "## Exact Current",
        "",
        "After separating the saddle offset and normalized correction,",
        "",
        "```text",
        exact["ideal_family"]["coefficients"],
        exact["ideal_family"]["parameters"],
        exact["ideal_family"]["polynomial"],
        exact["ideal_family"]["moment"],
        exact["ideal_family"]["current"],
        exact["ideal_family"]["chebyshev_expansion"],
        "```",
        "",
        "This identity comes directly from pairing the diagonal and",
        "off-diagonal terms of `|Q|^2+Re(H*conj(Q))`. It replaces",
        "the false absolute geometric-C1 target by a one-sided current.",
        "",
        "## Exact Bernstein Certificates",
        "",
        exact["bernstein_method"]["principle"],
        "",
        exact["bernstein_method"]["arithmetic"],
        "",
    ]
    for certificate in certificates:
        lines.extend(
            [
                "```text",
                (
                    f"p={certificate['prime']}, M={certificate['block_length']}, "
                    f"q>={certificate['q_lower']}, "
                    f"s>={certificate['s_lower']}"
                ),
                (
                    f"power degrees={certificate['power_degrees']}, "
                    f"Bernstein shape={certificate['bernstein_shape']}"
                ),
                (
                    "x leaves="
                    + repr(certificate["x_leaf_paths"])
                ),
                (
                    "leaf coefficient lower diagnostics="
                    + repr(certificate["leaf_minimum_lowers"])
                ),
                (
                    "global coefficient lower diagnostic="
                    + certificate["global_minimum_lower"]
                ),
                (
                    "proved rational current lower="
                    + certificate["claimed_rational_lower"]
                ),
                "```",
                "",
            ]
        )
    lines.extend(
        [
            "The displayed decimals are outward-rounded lower diagnostics.",
            "The checker decides every inequality using the underlying exact",
            "Fractions and outward rational square-root intervals.",
            "",
            "## Parameter Containment",
            "",
            "```text",
            exact["parameter_boxes"]["cutoff"],
            exact["parameter_boxes"]["dyadic"],
            exact["parameter_boxes"]["ternary"],
            exact["parameter_boxes"]["quinary_audit"],
            "```",
            "",
            "The variable `s` contains `u` only. The small positive",
            "`delta_a` factor is deliberately restored as a perturbation;",
            "this avoids silently assuming `u-delta_a>=0`.",
            "",
            "## Large Primes",
            "",
            "For every `p>=5` with a two-level complete chain,",
            "",
            "```text",
            exact["large_prime_family"]["reduction"],
            exact["large_prime_family"]["current"],
            exact["large_prime_family"]["minimum"],
            "```",
            "",
            "Thus the quinary Bernstein box is an exact cross-check, while",
            "the promoted result covers the whole `p>=5, M=2` family.",
            "",
            "## Actual Fixed-Ray Transfer",
            "",
            "The normalized actual coefficient factor is",
            "",
            "```text",
            exact["actual_transfer"]["normalized_factor"],
            exact["actual_transfer"]["coefficient_error"],
            exact["actual_transfer"]["frozen_current_error"],
            exact["actual_transfer"]["ray_terms"],
            exact["actual_transfer"]["ray_error"],
            "```",
            "",
            "Consequently,",
            "",
            "```text",
            exact["actual_transfer"]["margins"]["dyadic"],
            exact["actual_transfer"]["margins"]["ternary"],
            exact["actual_transfer"]["margins"]["large_prime"],
            "```",
            "",
            "The `pi` entering `x=4*pi*exp(L)` is the usual circle",
            "constant inherited from the completed-zeta saddle",
            "normalization. The angular identity uses the same ordinary",
            "`2*pi` period of `exp(i*theta)`; no new polygonal constant is",
            "inserted.",
            "",
            "## Boundary",
            "",
            payload["proof_boundary"],
            "",
        ]
    )
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
        "built prime-power heat-starlikeness base certificate: "
        "18 rows, 3 exact Bernstein base certificates, "
        "1 analytic p>=5 two-level family, "
        "3 actual ray-positive base families, "
        "0 length-propagation theorems, 0 joined Abel gaps, "
        "0 successor winding bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
