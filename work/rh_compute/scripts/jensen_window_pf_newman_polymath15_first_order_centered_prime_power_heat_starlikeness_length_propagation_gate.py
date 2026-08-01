#!/usr/bin/env python3
"""Build the prime-power heat-starlikeness length-propagation gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
from itertools import product
from pathlib import Path

from prime_power_heat_bernstein import certify_family


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_length_propagation_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
BASE_STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_base_certificate"
)
PHASE_STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_logarithmic_phase_flow_gate"
)
SOURCES = {
    "base_certificate": (
        REPO_ROOT / "work/rh_compute/results" / f"{BASE_STEM}.json"
    ),
    "phase_flow": (
        REPO_ROOT / "work/rh_compute/results" / f"{PHASE_STEM}.json"
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


def dyadic_m9_certificate() -> dict:
    paths = tuple(
        "".join(directions)
        for directions in product("LR", repeat=5)
    )
    return certify_family(
        prime=2,
        block_length=9,
        q_lower=Fraction(47, 50),
        s_lower=Fraction(22, 25),
        x_leaf_paths=paths,
        claimed_lower=Fraction(1, 20),
    ).to_dict()


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCES.items()
    }


def audit_sources(payloads: dict[str, dict]) -> dict:
    base = payloads["base_certificate"]
    if base.get("summary", {}).get(
        "exact_bernstein_base_certificates"
    ) != 3:
        raise RuntimeError("base-certificate summary drifted")
    phase = payloads["phase_flow"].get("exact", {})
    if "J_heat=|Q_heat|^2" not in phase.get(
        "replacement_target", {}
    ).get("current", ""):
        raise RuntimeError("phase-flow target drifted")
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
            "pphslp_01_domain",
            "complete-chain domain",
            "proved",
            "Work with the normalized correction-free heat coefficients on one complete prime-power chain.",
            "c_(k,M)=p^(-k/2)q^[k(2M-2-k)]s^k; 0<q,s<=1.",
            "Short, singleton, p-free, and endpoint joins are separate.",
        ),
        GateRow(
            "pphslp_02_append",
            "append-one polynomial identity",
            "proved",
            "Increasing the block length contracts the old polynomial and appends one monomial.",
            "Q_(M+1)(z)=Q_M(q^2*z)+d_M*z^M; d_M=p^(-M/2)s^M q^(M^2).",
            "This identity alone does not sign the appended cross current.",
        ),
        GateRow(
            "pphslp_03_current_recurrence",
            "append-one current identity",
            "proved",
            "The new current is the contracted old current plus one exact diagonal and one explicit cross term.",
            "J_(M+1)=J_M(q^2*z)+(M+1)d_M^2+d_M sum_(k<M)(M+k+2)c_(k,M)q^(2k)cos((M-k)theta).",
            "The cross term is not discarded or bounded by absolute C1 closeness.",
        ),
        GateRow(
            "pphslp_04_contracted_old",
            "interior positivity transfer",
            "proved",
            "A starlike old block remains positive after the radial contraction z -> q^2 z.",
            "Re(1+zQ_M'/Q_M)>0 on |z|=1 implies the same on |z|<=1 by harmonic minimum.",
            "A separate estimate is still required for the appended monomial.",
        ),
        GateRow(
            "pphslp_05_dyadic_m9_box",
            "next dyadic length",
            "proved",
            "Every ideal dyadic M=9 chain remains in the base rational parameter box.",
            "p=2,M=9: q in [47/50,1], s in [22/25,1], x=cos(theta) in [-1,1].",
            "This is one finite extension, not an all-M dyadic theorem.",
        ),
        GateRow(
            "pphslp_06_dyadic_m9_bernstein",
            "exact next-length certificate",
            "proved",
            "Thirty-two exact angular leaves certify the entire M=9 dyadic box.",
            "p=2,M=9: J_0>1/20.",
            "The exact Fraction certificate is rebuilt by the checker.",
        ),
        GateRow(
            "pphslp_07_fourier",
            "current Fourier sequence",
            "proved",
            "The ideal current is a cosine polynomial with positive explicit coefficients.",
            "J_0=A_0/2+sum_(d=1)^(M-1)A_d cos(d theta).",
            "Here A_0=2 sum(k+1)c_k^2 and A_d=sum_k(2k+d+2)c_k c_(k+d).",
        ),
        GateRow(
            "pphslp_08_ratios",
            "large-prime ratio geometry",
            "proved",
            "For p>=5 the adjacent ratios increase with k and stay below one half.",
            "r_k=c_(k+1)/c_k=(s/sqrt(p))q^(2M-3-2k); 0<r_k<=r_(k+1)<=1/sqrt(5)<1/2.",
            "The inequality is special to the p>=5 branch.",
        ),
        GateRow(
            "pphslp_09_monotone",
            "Fourier monotonicity",
            "proved",
            "The large-prime Fourier coefficients strictly decrease.",
            "A_0>A_1>...>A_(M-1)>A_M=0.",
            "Termwise comparison uses r_k<=1/sqrt(5).",
        ),
        GateRow(
            "pphslp_10_interior_convexity",
            "interior second differences",
            "proved",
            "Every common summand in A_d-2A_(d+1)+A_(d+2) is positive.",
            "A-2(A+1)r+(A+2)rr_+ >=(1-r)[A-(A+2)r]>0.",
            "This uses r_+>=r and r<1/2<=A/(A+2).",
        ),
        GateRow(
            "pphslp_11_boundary_reduction",
            "terminal second differences",
            "proved",
            "The two terminal summands reduce to one elementary family.",
            "G_(A,d)=A-2(A+1)aq+(A+2)a^2 q^(2d+2), a=s/sqrt(p), A>=d+2.",
            "No asymptotic or sampled inequality is used.",
        ),
        GateRow(
            "pphslp_12_boundary_minimization",
            "terminal parameter reduction",
            "proved",
            "G increases with A and decreases with a on the large-prime box.",
            "G_(A,d)>=g_d(q):=d+2-2(d+3)q/sqrt(5)+(d+4)q^(2d+2)/5.",
            "The derivative signs use aq<1/2 and (d+4)/sqrt(5)<d+3.",
        ),
        GateRow(
            "pphslp_13_finite_boundary",
            "small-d terminal certificate",
            "proved",
            "For 0<=d<=6 every degree-(2d+2) Bernstein coefficient of g_d is positive.",
            "2/sqrt(5)<9/10 gives b_j>0 for j<2d+2 and b_(2d+2)>(3d+1)/10.",
            "This is an exact rational Bernstein argument.",
        ),
        GateRow(
            "pphslp_14_tail_boundary",
            "large-d terminal certificate",
            "proved",
            "For d>=7 the positive last term can be dropped.",
            "g_d(q)>d+2-(9/10)(d+3)=(d-7)/10>=0.",
            "Strictness at d=7 comes from 2/sqrt(5)<9/10.",
        ),
        GateRow(
            "pphslp_15_fejer",
            "Fejer-kernel decomposition",
            "proved",
            "Convexity converts the full current into a nonnegative sum of Fejer kernels.",
            "J_0=(1/2)sum_(d=0)^(M-1)(d+1)Delta^2 A_d F_d(theta); F_d>=0.",
            "A_M=A_(M+1)=0 and F_0=1.",
        ),
        GateRow(
            "pphslp_16_first_difference",
            "uniform strict margin",
            "proved",
            "The first Fourier second difference has a uniform rational margin.",
            "Delta^2 A_0>=14/5-6/sqrt(5)>1/10.",
            "The final comparison is the exact inequality 729>720.",
        ),
        GateRow(
            "pphslp_17_large_prime_ideal",
            "all-length ideal theorem",
            "proved",
            "Every correction-free complete p>=5 chain is uniformly heat-starlike at every length.",
            "p>=5, M>=2: J_0>(1/2)Delta^2 A_0>1/20.",
            "This closes the ideal p>=5 length obligation.",
        ),
        GateRow(
            "pphslp_18_uniform_transfer",
            "all-length perturbation budget",
            "proved",
            "Geometric coefficient sums replace the former M<=8 counting bound.",
            "sum c_k<20/11; sum k c_k<180/121; |J_ang-J_0|<21r_x, r_x=18000/x.",
            "The saddle-offset estimate uses k log(p)<=log(N), so it is uniform in M.",
        ),
        GateRow(
            "pphslp_19_actual_margins",
            "actual fixed-ray propagation",
            "proved",
            "The dyadic M=9 family and every p>=5 complete length retain strict actual fixed-ray current.",
            "p=2,M=9 and p>=5,M>=2: |J_ray-J_0|<1/1000, hence J_ray>49/1000.",
            "The p=2 estimate uses the separate ratio bound 1/sqrt(2)<3/4.",
        ),
        GateRow(
            "pphslp_20_boundary",
            "proof boundary",
            "open",
            "The remaining complete-chain length obligations are dyadic M>=10 and ternary M>=5.",
            "Open: p=2,M>=10; p=3,M>=5; short/singleton chains; p-free and endpoint joins; Xi Abel gap.",
            "No successor winding, Lambda<=0, PF-infinity, RH, or prize conclusion is promoted.",
        ),
    ]


def build_payload() -> dict:
    sources = load_sources()
    rows = build_rows()
    certificate = dyadic_m9_certificate()
    if not certificate["positive"]:
        raise RuntimeError("dyadic M=9 certificate failed")
    return {
        "kind": STEM,
        "date": "2026-07-30",
        "status": (
            "proved one dyadic next-length certificate and all-length "
            "p>=5 ideal and actual fixed-ray propagation"
        ),
        "proof_boundary": (
            "This gate proves the exact append-one identity, J_0>1/20 "
            "and J_ray>49/1000 for the dyadic M=9 family, and "
            "J_0>1/20 and J_ray>49/1000 for every complete p>=5 "
            "chain of every length M>=2. It does not prove p=2,M>=10, "
            "p=3,M>=5, short or singleton chains, p-free or endpoint "
            "joins, a Xi Abel gap, a successor winding bound, "
            "Lambda<=0, PF-infinity, or RH."
        ),
        "source_audit": audit_sources(sources),
        "exact": {
            "append_one": {
                "coefficient": (
                    "c_(k,M+1)=q^(2k)c_(k,M), 0<=k<M"
                ),
                "new_term": (
                    "d_M=c_(M,M+1)=p^(-M/2)s^M q^(M^2)"
                ),
                "polynomial": (
                    "Q_(M+1)(z)=Q_M(q^2 z)+d_M z^M"
                ),
                "current": (
                    "J_(M+1)=J_M(q^2 z)+(M+1)d_M^2"
                    "+d_M sum_(k=0)^(M-1)(M+k+2)"
                    "c_(k,M)q^(2k)cos((M-k)theta)"
                ),
            },
            "dyadic_m9": {
                "certificate": certificate,
                "ideal_margin": "J_0>1/20",
                "actual_margin": "J_ray>49/1000",
            },
            "large_prime_fejer": {
                "fourier": (
                    "A_0=2sum_(k=0)^(M-1)(k+1)c_k^2; "
                    "A_d=sum_(k=0)^(M-1-d)(2k+d+2)c_kc_(k+d)"
                ),
                "ratio": (
                    "r_k=(s/sqrt(p))q^(2M-3-2k); "
                    "0<r_k<=r_(k+1)<=1/sqrt(5)<1/2"
                ),
                "kernel": (
                    "F_d(theta)=1+2sum_(m=1)^d"
                    "(1-m/(d+1))cos(mtheta)>=0"
                ),
                "decomposition": (
                    "J_0=(1/2)sum_(d=0)^(M-1)"
                    "(d+1)(A_d-2A_(d+1)+A_(d+2))F_d"
                ),
                "first_margin": (
                    "Delta^2 A_0>=14/5-6/sqrt(5)>1/10"
                ),
                "theorem": "p>=5,M>=2: J_0>1/20",
            },
            "uniform_actual_transfer": {
                "coefficient_error": "|R_k-1|<r_x=18000/x",
                "large_prime_sums": (
                    "sum c_k<20/11; sum k c_k<180/121"
                ),
                "large_prime_frozen_error": (
                    "|J_ang-J_0|<21r_x=378000/x"
                ),
                "large_prime_ray_error": (
                    "|J_ray-J_0|<378006/x+405408/x^2<1/1000"
                ),
                "dyadic_m9_ray_error": (
                    "|J_ray-J_0|<3456096/x+1621632/x^2<1/1000"
                ),
                "margins": (
                    "p=2,M=9 and p>=5,M>=2: J_ray>49/1000"
                ),
            },
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "append_one_identities": 1,
            "exact_dyadic_next_length_certificates": 1,
            "analytic_all_length_large_prime_families": 1,
            "actual_ray_positive_propagated_families": 2,
            "open_small_prime_length_families": 2,
            "proved_joined_abel_gaps": 0,
            "proved_successor_winding_bounds": 0,
        },
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    certificate = exact["dyadic_m9"]["certificate"]
    return "\n".join(
        [
            "# Prime-Power Heat-Starlikeness Length Propagation",
            "",
            "Date: 2026-07-30",
            "",
            "Status: one exact dyadic next-length certificate and a",
            "uniform all-length theorem for every complete p>=5 chain.",
            "This does not promote RH or Lambda<=0.",
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
                "validated prime-power heat-starlikeness length "
                "propagation gate: 20 rows, 1 append-one identity, "
                "1 exact dyadic next-length certificate, "
                "1 analytic all-length p>=5 family, "
                "2 actual ray-positive propagated families, "
                "2 open small-prime length families, "
                "0 joined Abel gaps, 0 successor winding bounds"
            ),
            "```",
            "",
            "## Append-One Identity",
            "",
            "```text",
            exact["append_one"]["coefficient"],
            exact["append_one"]["new_term"],
            exact["append_one"]["polynomial"],
            exact["append_one"]["current"],
            "```",
            "",
            "The old block is evaluated at `q^2 z`; the final displayed",
            "sum is the entire new cross current. No absolute geometric",
            "`C1` approximation is inserted.",
            "",
            "## Dyadic M=9",
            "",
            "The exact rational tensor-Bernstein certificate gives",
            "",
            "```text",
            "p=2, M=9, q>=47/50, s>=22/25",
            f"power degrees={certificate['power_degrees']}",
            f"Bernstein shape={certificate['bernstein_shape']}",
            f"angular leaves={len(certificate['x_leaf_paths'])}",
            (
                "global coefficient lower diagnostic="
                + certificate["global_minimum_lower"]
            ),
            "J_0>1/20",
            "```",
            "",
            "Thus the first unresolved dyadic extension is positive; this",
            "is not yet an induction for `M>=10`.",
            "",
            "## All Large-Prime Lengths",
            "",
            "Define the Fourier coefficients by",
            "",
            "```text",
            exact["large_prime_fejer"]["fourier"],
            exact["large_prime_fejer"]["ratio"],
            "```",
            "",
            "The ratio sequence is increasing and bounded by",
            "`1/sqrt(5)<1/2`. Direct grouping proves that `A_d` is",
            "decreasing and convex after adjoining `A_M=A_(M+1)=0`.",
            "The only terminal expression is",
            "",
            "```text",
            "G_(A,d)=A-2(A+1)aq+(A+2)a^2q^(2d+2)",
            "G_(A,d)>=d+2-2(d+3)q/sqrt(5)+(d+4)q^(2d+2)/5.",
            "```",
            "",
            "For `0<=d<=6`, all ordinary degree-`2d+2` Bernstein",
            "coefficients are positive using `2/sqrt(5)<9/10`. For",
            "`d>=7`, dropping the last positive term leaves",
            "`(d-7)/10>=0`.",
            "",
            "Fejer's exact identity is",
            "",
            "```text",
            exact["large_prime_fejer"]["kernel"],
            exact["large_prime_fejer"]["decomposition"],
            exact["large_prime_fejer"]["first_margin"],
            exact["large_prime_fejer"]["theorem"],
            "```",
            "",
            "The circle constant does not enter this argument. The",
            "`cos(d theta)` terms have the ordinary `2*pi` angular period.",
            "",
            "## Actual Fixed-Ray Transfer",
            "",
            "Because `k log(p)<=log(N)`, the relative coefficient estimate",
            "is uniform in chain length. Geometric sums replace the old",
            "`M<=8` count:",
            "",
            "```text",
            exact["uniform_actual_transfer"]["coefficient_error"],
            exact["uniform_actual_transfer"]["large_prime_sums"],
            exact["uniform_actual_transfer"]["large_prime_frozen_error"],
            exact["uniform_actual_transfer"]["large_prime_ray_error"],
            exact["uniform_actual_transfer"]["dyadic_m9_ray_error"],
            exact["uniform_actual_transfer"]["margins"],
            "```",
            "",
            "## Boundary",
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
    arguments = build_parser().parse_args()
    payload = build_payload()
    arguments.out.parent.mkdir(parents=True, exist_ok=True)
    arguments.note.parent.mkdir(parents=True, exist_ok=True)
    arguments.out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    arguments.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "built prime-power heat-starlikeness length propagation gate: "
        "20 rows, 1 append-one identity, "
        "1 exact dyadic next-length certificate, "
        "1 analytic all-length p>=5 family, "
        "2 actual ray-positive propagated families, "
        "2 open small-prime length families, "
        "0 joined Abel gaps, 0 successor winding bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
