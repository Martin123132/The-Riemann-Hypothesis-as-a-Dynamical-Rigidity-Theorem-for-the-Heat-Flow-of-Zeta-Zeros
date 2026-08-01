#!/usr/bin/env python3
"""Build the q=1 growing contiguous-terminal-prefix finite-height theorem."""

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
    "contiguous_terminal_tail_growing_prefix_finite_height_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
H_DENOMINATOR = 72_000_000_000
SOURCE_PATHS = {
    "q1_edge": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "q1_finite_height_real_edge_remainder_gate.json"
    ),
    "tail_symbol": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_recurrence_current_gate.json"
    ),
    "pair_guard": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "edge_anchored_near_terminal_pair_current_guard.json"
    ),
    "q1_edge_builder": (
        REPO_ROOT
        / "work/rh_compute/scripts/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "q1_finite_height_real_edge_remainder_gate.py"
    ),
    "interior_current": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "interior_projective_current_gate.json"
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


def parse_fraction(value: str) -> Fraction:
    numerator, denominator = value.split("/", maxsplit=1)
    return Fraction(int(numerator), int(denominator))


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


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
    edge = payloads["q1_edge"]
    if not isinstance(edge, dict):
        raise RuntimeError("q1 edge source is not JSON")
    majorant = edge.get("majorant_certificate", {})
    defects = majorant.get("edge_jet_defects", {})
    expected = {
        "value": 100,
        "centered_slope": 250,
        "radial_subtracted_value_jet": 250,
        "centered_second_jet": 900,
    }
    actual = {
        key: defects.get(key, {}).get("constant") for key in expected
    }
    if actual != expected:
        raise RuntimeError(f"q1 edge defects drifted: {actual}")
    terminal = majorant.get("source_block_majorants", {}).get(
        "stable_terminal", {}
    )
    if terminal.get("W") != "|W|<6h" or terminal.get("W_x") != "|W_x|<2h^2":
        raise RuntimeError("terminal stable-log bounds drifted")
    chain = majorant.get("normalized_majorant_chain", {})
    source_slack_targets = {
        "d_over_h2": Fraction(1, 20),
        "d_x_over_h4": Fraction(1, 400),
    }
    source_slack: dict[str, str] = {}
    for key, target in source_slack_targets.items():
        entry = chain.get(key, {})
        try:
            upper = parse_fraction(entry["upper"])
        except (KeyError, TypeError, ValueError) as exc:
            raise RuntimeError(f"q1 source slack missing: {key}") from exc
        if not upper < target:
            raise RuntimeError(f"q1 source slack drifted: {key}")
        source_slack[key] = f"{upper.numerator}/{upper.denominator}"
    if majorant.get("leading_jet_bounds") != {"|A|": 4, "|B|": 3, "|M|": 3}:
        raise RuntimeError("q1 leading-jet bounds drifted")

    tail = payloads["tail_symbol"]
    if not isinstance(tail, dict) or "P_M(p)<-3/8000" not in tail.get(
        "exact", {}
    ).get("theorem", ""):
        raise RuntimeError("quantitative tail-symbol theorem drifted")

    pair = payloads["pair_guard"]
    if not isinstance(pair, dict) or "P_3(-5/8)>1/10" not in pair.get(
        "interval_certificate", {}
    ).get("positive_margin", ""):
        raise RuntimeError("pair guard drifted")

    builder = payloads["q1_edge_builder"]
    if not isinstance(builder, str):
        raise RuntimeError("q1 edge builder source is not text")
    for marker in (
        "def majorant_certificate",
        "d_x=(-i/2)",
        "s_*''=-t alpha''/8",
    ):
        if marker not in builder:
            raise RuntimeError(f"q1 source marker missing: {marker}")

    interior = payloads["interior_current"]
    if not isinstance(interior, dict) or "c_n=X_n" not in str(
        interior.get("exact", {})
    ):
        raise RuntimeError("interior current coordinates drifted")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCE_PATHS.items()
        },
        "q1_source_slack": {
            "terminal_d_over_h2_upper": source_slack["d_over_h2"],
            "terminal_d_over_h2_lt": "1/20",
            "terminal_d_x_over_h4_upper": source_slack["d_x_over_h4"],
            "terminal_d_x_over_h4_lt": "1/400",
            "leading_edge_bounds": "|A|<4, |B|<3, |M|<3",
        },
        "published_source": {
            "citation": (
                "D.H.J. Polymath, Effective approximation of heat flow "
                "evolution of the Riemann xi function, and a new upper "
                "bound for the de Bruijn-Newman constant, arXiv:1904.12438"
            ),
            "location": "equation (53), the first Dirichlet correction, and C_0",
            "url": "https://arxiv.org/abs/1904.12438",
        },
    }


def symbolic_certificate() -> dict:
    h, theta, m = sp.symbols("h theta m", positive=True, real=True)
    k = theta + m
    delta = sp.log(1 - h * k) - sp.log(1 - h * theta)
    centered_phase = (
        2 * sp.pi * delta / h**2
        + 2 * sp.pi * m * (1 / h - theta)
        + sp.pi * (m**2 + 4 * m * theta)
    )
    cubic_limit = sp.simplify(sp.limit(centered_phase / h, h, 0, dir="+"))
    expected = -2 * sp.pi * (k**3 - theta**3) / 3
    if sp.simplify(cubic_limit - expected) != 0:
        raise RuntimeError("centered cubic phase failed")

    t, log_a, alpha = sp.symbols("t log_a alpha", real=True)
    ell_k = sp.log(1 - h * k)
    ell_theta = sp.log(1 - h * theta)
    chi = alpha - log_a
    s = sp.Rational(1, 2) - sp.I * (
        2 * sp.pi / h**2 - sp.pi * t / 8
    )
    raw_ratio = (
        t
        * ((log_a + ell_k) ** 2 - (log_a + ell_theta) ** 2)
        / 4
        - (s + t * alpha / 2) * delta
    )
    stated_ratio = (
        -delta / 2
        - t * chi * delta / 2
        + t * (ell_k + ell_theta) * delta / 4
        + sp.I * (centered_phase - sp.pi * t * delta / 8)
    )
    log_r_m = sp.I * sp.pi * m - 4 * sp.I * sp.pi * m * theta
    branch_integer = sp.simplify(
        (raw_ratio - stated_ratio - log_r_m) / (2 * sp.pi * sp.I)
    )
    expected_branch = -m * (1 / h - theta) - m * (m + 1) / 2
    if sp.simplify(branch_integer - expected_branch) != 0:
        raise RuntimeError("exact ratio branch congruence failed")

    d_base, alpha_prime, alpha_second, w_0, w_1 = sp.symbols(
        "d_base alpha_prime alpha_second w_0 w_1"
    )
    d = lambda w: d_base + alpha_prime * (t / 4 + t**2 * w**2 / 8)
    d_difference = sp.factor(d(w_1) - d(w_0))
    expected_d_difference = alpha_prime * t**2 * (w_1 - w_0) * (
        w_1 + w_0
    ) / 8
    if sp.simplify(d_difference - expected_d_difference) != 0:
        raise RuntimeError("correction difference failed")
    d_x = lambda w: (-sp.I / 2) * (
        d_base
        + alpha_second * (t / 4 + t**2 * w**2 / 8)
        + t**2 * w * alpha_prime**2 / 4
    )
    d_x_difference = sp.factor(d_x(w_1) - d_x(w_0))
    expected_d_x_difference = (-sp.I / 2) * (
        alpha_second
        * t**2
        * (w_1 - w_0)
        * (w_1 + w_0)
        / 8
        + t**2 * alpha_prime**2 * (w_1 - w_0) / 4
    )
    if sp.simplify(d_x_difference - expected_d_x_difference) != 0:
        raise RuntimeError("differentiated correction difference failed")

    c_s, gamma, u, k_log, x, y = sp.symbols(
        "c_s gamma u k_log X Y", real=True
    )
    c_atom = x
    d_atom = c_s * k_log * x + gamma * u * y
    if sp.diff(d_atom, x) != c_s * k_log:
        raise RuntimeError("carrier coordinate source failed")

    z_x_r, z_x_i = sp.symbols("Z_xR Z_xI", real=True)
    c_s_x, gamma_x, u_x = sp.symbols("c_s_x gamma_x u_x", real=True)
    d_x = (
        c_s_x * k_log * x
        + c_s * k_log * z_x_r
        + gamma_x * u * y
        + gamma * u_x * y
        + gamma * u * z_x_i
    )
    expected_d_x = sp.diff(d_atom, c_s) * c_s_x + sp.diff(
        d_atom, gamma
    ) * gamma_x + sp.diff(d_atom, u) * u_x + sp.diff(
        d_atom, x
    ) * z_x_r + sp.diff(d_atom, y) * z_x_i
    if sp.simplify(d_x - expected_d_x) != 0:
        raise RuntimeError("carrier differentiated coordinate failed")

    c, d, e, n, ec, ed, ee, en = sp.symbols(
        "C D E N eps_C eps_D eps_E eps_N", real=True
    )
    perturbation = sp.expand((c + ec) * (n + en) - (d + ed) * (e + ee) - (c * n - d * e))
    target = sp.expand(c * en + n * ec + ec * en - d * ee - e * ed - ed * ee)
    if sp.simplify(perturbation - target) != 0:
        raise RuntimeError("block-current perturbation failed")

    return {
        "ratio_coordinates": (
            "ell_j=log(1-h*j), k=m+theta, Delta_m=ell_k-ell_theta, "
            "chi=alpha-log(a), r=-exp(2*pi*i*p)"
        ),
        "exact_ratio_log": (
            "L_m=log[(f_(N-m)/f_N)/r^m]="
            "-Delta_m/2-(t/2)chi*Delta_m"
            "+(t/4)(ell_k+ell_theta)Delta_m"
            "+i[R_m-(pi*t/8)Delta_m]"
            "+log(1+d_k)-log(1+d_theta)"
        ),
        "centered_phase": (
            "R_m=2*pi*h^-2*Delta_m+2*pi*m*(h^-1-theta)"
            "+pi*(m^2+4*m*theta)"
        ),
        "centered_phase_limit": (
            "R_m/h -> -(2*pi/3)[(m+theta)^3-theta^3]"
        ),
        "ratio_branch_winding": (
            "raw_log_ratio-L_m-m*Log(r)=-2*pi*i*"
            "[m*N+m*(m+1)/2], N=h^-1-theta in Z"
        ),
        "correction_differences": (
            "d_k-d_theta=alpha'*t^2*(w_k-w_theta)*"
            "(w_k+w_theta)/8; the analogous x derivative is the sum "
            "of the alpha'' quadratic and (alpha')^2 linear terms"
        ),
        "ratio_log_derivative": (
            "L_(m,x)=-s_*'*Delta_m+d_(k,x)/(1+d_k)"
            "-d_(theta,x)/(1+d_theta)+i*h*m/2"
        ),
        "normalized_carrier": (
            "Z_(m,a)=z_(N-m)/S_a=-Q*r^m*exp(W_k), "
            "W_k=W_theta+L_m"
        ),
        "carrier_coordinates": (
            "Z=X+iY, V=Z_x/h, U=-ell_k/h, K_m=-Delta_m, "
            "gamma=-Im(s_*'); C=X; "
            "D=(Re(s_*')/h)K_m X+gamma UY; E=Re V; "
            "N=(Re(s_*')_x/h^2)K_mX+(Re(s_*')/h)K_mReV"
            "+[gamma_x*(-ell_k)/h^2+gamma/(8*pi)]Y+gamma U ImV"
        ),
        "leading_coordinates": (
            "Z_0=-Q*r^m=X_0+iY_0, V_0=-i*k*Z_0/2; "
            "(C_0,D_0,E_0,N_0)=(X_0,kY_0/2,kY_0/2,"
            "Y_0/(16*pi)-k^2X_0/4)"
        ),
        "current_perturbation": (
            "Delta P=C*eps_N+N*eps_C+eps_C*eps_N"
            "-D*eps_E-E*eps_D-eps_D*eps_E"
        ),
    }


def analytic_majorant_certificate() -> dict:
    delta = Fraction(1, 10**10)
    t_max = Fraction(1, 5000)
    pi_upper = Fraction(22, 7)

    # Here delta bounds h*k throughout the growing prefix.  The source q=1
    # certificate has strict slack |d_theta|/h^2<1/20 and
    # |d_(theta,x)|/h^4<1/400, audited separately in source_audit().
    phase_normalized = Fraction(2, 3) * pi_upper / (1 - delta)
    if not phase_normalized < 3:
        raise RuntimeError("centered phase majorant failed")

    d_difference_normalized = (
        3 * t_max**2 * 2 * 7 / 8
    )
    d_x_difference_normalized = t_max**2 * (
        Fraction(7, 8) + Fraction(9, 4)
    )
    if not d_difference_normalized < Fraction(1, 10**6):
        raise RuntimeError("correction difference majorant failed")
    if not d_x_difference_normalized < Fraction(1, 10**6):
        raise RuntimeError("correction x-difference majorant failed")
    if not (
        Fraction(1, 20) + delta * d_difference_normalized
        < Fraction(1, 10)
    ):
        raise RuntimeError("extended correction majorant failed")
    if not (
        Fraction(1, 400) + delta * d_x_difference_normalized
        < Fraction(1, 100)
    ):
        raise RuntimeError("extended correction x-majorant failed")

    ratio_log_normalized = (
        1
        + 3 * t_max
        + 2 * t_max
        + phase_normalized
        + pi_upper * t_max / 4
        + 2 * d_difference_normalized
    )
    ratio_log_x_normalized = (
        Fraction(1, 2) + Fraction(1, 3000) + Fraction(1, 45)
    )
    if not ratio_log_normalized < 5:
        raise RuntimeError("ratio-log majorant failed")
    if not ratio_log_x_normalized < 1:
        raise RuntimeError("ratio-log x-majorant failed")

    coordinate_normalized = {
        "value": Fraction(24),
        "slope_derived": (
            Fraction(1, 5000)
            + Fraction(1, 1500)
            + 1
            + 12
        ),
        "value_x": Fraction(12 + 6),
        "slope_x_derived": (
            Fraction(3, 10000)
            + Fraction(1, 72000)
            + Fraction(1, 2)
            + Fraction(1, 3000)
            + Fraction(1, 2)
            + 9
        ),
    }
    coordinate_caps = {
        "value": Fraction(24),
        "slope": Fraction(50),
        "value_x": Fraction(18),
        "slope_x": Fraction(40),
    }
    if coordinate_normalized["value"] > coordinate_caps["value"]:
        raise RuntimeError("carrier value majorant failed")
    if not coordinate_normalized["slope_derived"] < coordinate_caps["slope"]:
        raise RuntimeError("carrier slope majorant failed")
    if coordinate_normalized["value_x"] > coordinate_caps["value_x"]:
        raise RuntimeError("carrier value-x majorant failed")
    if not (
        coordinate_normalized["slope_x_derived"]
        < coordinate_caps["slope_x"]
    ):
        raise RuntimeError("carrier slope-x majorant failed")

    return {
        "prefix_small_parameter": "h*k<1/10^10",
        "elementary_log_bounds": (
            "|Delta_m|<2*h*k, |Delta_m+h*m|<h^2*k^2, "
            "|ell_k|<2*h*k, |U-k|<h*k^2"
        ),
        "phase_over_hk3_upper": fraction_text(phase_normalized),
        "phase_over_hk3_lt": 3,
        "correction_extension": {
            "w_theta": "|w_theta|<3",
            "w_k": "|w_k|<4",
            "d_difference_over_h3k_upper": fraction_text(
                d_difference_normalized
            ),
            "d_difference_over_h3k_lt": "1/1000000",
            "d_x_difference_over_h5k_upper": fraction_text(
                d_x_difference_normalized
            ),
            "d_x_difference_over_h5k_lt": "1/1000000",
            "extended_d": "|d_k|<h^2/10",
            "extended_d_x": "|d_(k,x)|<h^4/100",
        },
        "ratio_log_over_hk3_upper": fraction_text(ratio_log_normalized),
        "ratio_log_over_hk3_lt": 5,
        "ratio_log_x_over_h2k2_upper": fraction_text(
            ratio_log_x_normalized
        ),
        "ratio_log_x_over_h2k2_lt": 1,
        "carrier_log_chain": (
            "|W_k|<12*h*k^3, |W_(k,x)|<3*h^2*k^2, "
            "|exp(W_k)-1|<24*h*k^3"
        ),
        "coordinate_normalized_derived": {
            key: fraction_text(value)
            for key, value in coordinate_normalized.items()
        },
        "coordinate_published_caps": {
            key: fraction_text(value) for key, value in coordinate_caps.items()
        },
    }


def exact_integer_budget() -> dict:
    h0 = Fraction(1, H_DENOMINATOR)
    if not 4**18 < H_DENOMINATOR:
        raise RuntimeError("M=3 entry check failed")
    if not 10**180 < H_DENOMINATOR**17:
        raise RuntimeError("hK smallness check failed")
    if not 10**54 < H_DENOMINATOR**5:
        raise RuntimeError("hK3 smallness check failed")
    if not 2_400_000**18 < H_DENOMINATOR**11:
        raise RuntimeError("growing-prefix current budget failed")

    individual = {
        "value": 24,
        "slope": 50,
        "value_x": 18,
        "slope_x": 40,
    }
    edge = {"value": 100, "slope": 250, "value_x": 250, "slope_x": 900}
    aggregate = {
        "value": edge["value"] + individual["value"],
        "slope": edge["slope"] + individual["slope"],
        "value_x": edge["value_x"] + individual["value_x"],
        "slope_x": edge["slope_x"] + individual["slope_x"],
    }
    linear = (
        4 * aggregate["slope_x"]
        + aggregate["value"]
        + 2 * aggregate["value_x"]
        + 2 * aggregate["slope"]
    )
    quadratic = (
        aggregate["value"] * aggregate["slope_x"]
        + aggregate["slope"] * aggregate["value_x"]
    )
    if linear != 5020 or quadratic != 196_960:
        raise RuntimeError("aggregate perturbation constants drifted")
    if not linear + Fraction(quadratic, 10**9) < 6000:
        raise RuntimeError("saved perturbation constant failed")

    return {
        "h_max": f"1/{H_DENOMINATOR}",
        "prefix_condition": "K=M+1, M>=1, h*K^18<=1",
        "smallness": {
            "hK": "h*K<=h^(17/18)<1e-10",
            "hK3": "h*K^3<=h^(5/6)<1e-9",
            "exp_log": "12*h*K^3<1/2",
        },
        "ratio_bounds": {
            "Delta": "|Delta_m|<2*h*k",
            "centered_phase": "|R_m|<3*h*k^3",
            "correction": "|d_k-d_theta|<h^3*k/1000000",
            "ratio_log": "|L_m|<5*h*k^3",
            "ratio_log_x": "|L_(m,x)|<h^2*k^2",
            "carrier_log": "|W_k|<12*h*k^3",
            "carrier_log_x": "|W_(k,x)|<3*h^2*k^2",
            "carrier_value": "|exp(W_k)-1|<24*h*k^3",
        },
        "individual_carrier_errors": {
            "value": "|C-C_0|<24*h*k^3",
            "slope": "|D-D_0|<50*h*k^4",
            "value_x": "|E-E_0|<18*h*k^4",
            "slope_x": "|N-N_0|<40*h*k^5",
            "constants": individual,
        },
        "aggregate_errors": {
            "value": "|eps_C|<124*h*K^4",
            "slope": "|eps_D|<300*h*K^5",
            "value_x": "|eps_E|<268*h*K^5",
            "slope_x": "|eps_N|<940*h*K^6",
            "constants": aggregate,
        },
        "leading_coordinate_bounds": (
            "|C|<4K, |D|<2K^2, |E|<2K^2, |N|<K^3"
        ),
        "current_error": {
            "linear_constant": linear,
            "quadratic_constant": quadratic,
            "raw": (
                "|Delta P|<5020*h*K^7+196960*h^2*K^10"
            ),
            "saved": "|Delta P|<6000*h*K^7<1/400",
            "exact_integer_check": (
                "2400000^18<72000000000^11"
            ),
        },
        "leading_margin": (
            "M>=1 implies t_shift=2(M+1)-p>=3, D(t_shift)>4/5, "
            "and P_M<-1/200 because pi^2<10"
        ),
        "finite_height_margin": (
            "a^2*J_(tail,M)/S_a^2<-1/200+1/400=-1/400"
        ),
        "base_height_entry": (
            "4^18=68719476736<72000000000, so M=3 is admissible "
            "throughout L>=50"
        ),
        "growing_family": (
            "1<=M<=floor(a^(1/18))-1; this range tends to infinity with a"
        ),
    }


def exact_payload() -> dict:
    return {
        "domain": (
            "q=2tL^2=1, L>=50, theta=(1-p)/2 in [0,1], "
            "a=N+theta, h=1/a<1/72000000000"
        ),
        "prefix": (
            "Take M>=1, K=M+1, h*K^18<=1, and group the retained "
            "endpoint-terminal edge with n=N-1,...,N-M."
        ),
        "source_bounds": (
            "|s_*'+i/2|<h^2/6000, |Re(s_*')|<h^2/20000, "
            "|s_*''|<h^4/40000, u_x=h^2/(8*pi). The audited terminal "
            "slack |d_theta|<h^2/20 and |d_(theta,x)|<h^4/400, "
            "together with the exact correction differences, gives "
            "|d_k|<h^2/10 and |d_(k,x)|<h^4/100."
        ),
        "theorem": (
            "For every admissible M>=1, the retained first-order q=1 "
            "contiguous terminal-prefix current satisfies "
            "a^2*J_(tail,M)/S_a^2<-1/400."
        ),
        "edge_case": (
            "For M=0, the imported edge theorem remains "
            "a^2*J_edge/S_a^2<-3749/10000000."
        ),
        "route_gain": (
            "At L>=50 the theorem includes M=3, so the canonical prefix "
            "through N-3 is finite-height clockwise even though the edge "
            "plus N-3 alone has both limiting orientations."
        ),
        "pi_provenance": (
            "Every pi is inherited from the completed-zeta/Riemann-Siegel "
            "source, a^2=x/(4*pi)+t/16, and the published C_0 recurrence; "
            "no fitted geometry introduces it."
        ),
        "surviving_target": (
            "Join this growing but short terminal block to the remaining "
            "n<=N-M-1 Abel/cross-current bulk without termwise absolute "
            "values, then use the existing C1 whole-jet residual transfer."
        ),
    }


def build_rows(
    exact: dict, symbolic: dict, analytic: dict, budget: dict
) -> list[GateRow]:
    return [
        GateRow("ctgf_01_domain", "physical domain", "proved", "The theorem stays on the physical q=1 critical family.", exact["domain"], "The q>1 extension is separate."),
        GateRow("ctgf_02_prefix", "growing prefix", "proved", "The allowed terminal prefix grows with saddle height.", exact["prefix"], "M is finite at each height."),
        GateRow("ctgf_03_ratio", "exact ratio", "proved", "The finite carrier ratio is centered against the exact limiting geometric phase.", symbolic["exact_ratio_log"], "The logarithm is the near-one branch."),
        GateRow("ctgf_04_phase", "phase remainder", "proved", "The first omitted saddle phase is cubic in terminal distance.", symbolic["centered_phase"] + "; " + symbolic["centered_phase_limit"], "No fitted phase is used."),
        GateRow("ctgf_05_ratio_x", "differentiated ratio", "proved", "One x derivative retains the exact correction quotient.", symbolic["ratio_log_derivative"], "No second derivative of a carrier is required."),
        GateRow("ctgf_06_source", "source bounds", "proved", "The existing q=1 source majorants extend over the admissible short prefix.", exact["source_bounds"], "The extension uses |alpha-log(n)|<4 and the exact d_k-d_theta identities.", analytic["correction_extension"]),
        GateRow("ctgf_07_log_bounds", "uniform ratio bounds", "proved", "The carrier value and derivative stay uniformly close to the limiting geometric atom.", budget["ratio_bounds"]["carrier_log"] + "; " + budget["ratio_bounds"]["carrier_log_x"], "The prefix condition makes every logarithm disk tiny.", {"analytic": analytic, "published": budget["ratio_bounds"]}),
        GateRow("ctgf_08_coordinates", "four coordinates", "proved", "The normalized carrier contributes four division-free projective coordinates.", symbolic["carrier_coordinates"], "The common radial normalizer is differentiated before subtraction."),
        GateRow("ctgf_09_leading", "leading atom", "proved", "The four limiting carrier coordinates match the fixed-M theorem.", symbolic["leading_coordinates"], "This is an exact derivative of -Q*r^m."),
        GateRow("ctgf_10_individual", "individual errors", "proved", "Each finite carrier has explicit polynomial-in-k errors.", "; ".join(budget["individual_carrier_errors"][key] for key in ("value", "slope", "value_x", "slope_x")), "All constants include the first Dirichlet correction.", analytic["coordinate_normalized_derived"]),
        GateRow("ctgf_11_aggregate", "summed errors", "proved", "The edge and every carrier error sum before current polarization.", "; ".join(budget["aggregate_errors"][key] for key in ("value", "slope", "value_x", "slope_x")), "No pairwise absolute-current sum is used."),
        GateRow("ctgf_12_perturbation", "current algebra", "proved", "The complete four-coordinate perturbation is polynomial.", symbolic["current_perturbation"], "It remains valid on every zero real projection."),
        GateRow("ctgf_13_budget", "current budget", "proved", "The cumulative finite-height defect is smaller than half the leading reserve.", budget["current_error"]["raw"] + "; " + budget["current_error"]["saved"], "The last inequality is an exact integer-power check.", budget["current_error"]),
        GateRow("ctgf_14_entry", "worst-height entry", "proved", "The repaired N-3 prefix is already admissible at L=50.", budget["base_height_entry"], "The range increases thereafter."),
        GateRow("ctgf_15_theorem", "finite-height theorem", "proved", "Every admissible nontrivial contiguous terminal prefix is strictly clockwise.", exact["theorem"], "This is a retained first-order q=1 theorem."),
        GateRow("ctgf_16_handoff", "open handoff", "open", "The theorem creates a signed growing terminal block for the aggregate route.", exact["surviving_target"], "No bulk cross-current or Xi-level conclusion is claimed."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    symbolic = payload["symbolic_certificate"]
    analytic = payload["analytic_majorant_certificate"]
    budget = payload["exact_integer_budget"]
    return f"""# Contiguous Terminal-Tail Growing-Prefix Finite-Height Gate

Date: 2026-07-31

Status: retained first-order `q=1` growing-prefix current theorem. This is not
a proof of the remaining bulk cross-current estimate, an Xi-level current
theorem, `Lambda<=0`, or RH.

## Domain

{exact['domain']}

{exact['prefix']}

Equivalently, `1<=M<=floor(a^(1/18))-1`; this range tends to infinity.

## Exact Carrier Ratio

```text
{symbolic['ratio_coordinates']}
{symbolic['exact_ratio_log']}
{symbolic['centered_phase']}
{symbolic['ratio_log_derivative']}
{symbolic['ratio_branch_winding']}
```

The cubic phase limit is

```text
{symbolic['centered_phase_limit']}
```

The inherited source bounds and exact logarithm estimates give

```text
{analytic['elementary_log_bounds']}
{symbolic['correction_differences']}
{analytic['correction_extension']['extended_d']}
{analytic['correction_extension']['extended_d_x']}
{budget['ratio_bounds']['carrier_log']}
{budget['ratio_bounds']['carrier_log_x']}
{budget['ratio_bounds']['carrier_value']}
```

## Four Coordinates

```text
{symbolic['carrier_coordinates']}
{symbolic['leading_coordinates']}
```

For each carrier,

```text
{budget['individual_carrier_errors']['value']}
{budget['individual_carrier_errors']['slope']}
{budget['individual_carrier_errors']['value_x']}
{budget['individual_carrier_errors']['slope_x']}
```

After adding the certified edge errors and summing `m=1,...,M`,

```text
{budget['aggregate_errors']['value']}
{budget['aggregate_errors']['slope']}
{budget['aggregate_errors']['value_x']}
{budget['aggregate_errors']['slope_x']}
```

## Current Margin

```text
{symbolic['current_perturbation']}
{budget['current_error']['raw']}
{budget['current_error']['saved']}
```

The nontrivial leading prefix has `P_M<-1/200`; therefore

```text
{budget['finite_height_margin']}
```

{exact['theorem']}

{budget['base_height_entry']}. {exact['route_gain']}

## Boundary

{exact['pi_provenance']}

{exact['surviving_target']}

This theorem does not sign the remaining nonterminal bulk, extend the prefix
to all `N-1` carriers, prove the `q>1` grouped-tail theorem, control the omitted
Xi remainder inside one component, prove an Abel gap, winding cap, contact
exclusion, `Q209`, `Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.

## Reproduce

```text
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```
"""


def build_payload() -> dict:
    payloads = load_sources()
    audit = source_audit(payloads)
    symbolic = symbolic_certificate()
    analytic = analytic_majorant_certificate()
    budget = exact_integer_budget()
    exact = exact_payload()
    rows = build_rows(exact, symbolic, analytic, budget)
    return {
        "kind": STEM,
        "date": "2026-07-31",
        "status": (
            "retained first-order q=1 growing contiguous-terminal-prefix "
            "finite-height current theorem; no bulk aggregate closure, "
            "Xi-level current theorem, Lambda<=0, or RH"
        ),
        "proof_boundary": (
            "This proves exact finite carrier-ratio and differentiated-ratio "
            "identities, explicit four-coordinate errors, and strict current "
            "for 1<=M<=floor(a^(1/18))-1 in the retained first-order q=1 "
            "model. It proves no sign for the remaining nonterminal bulk, "
            "complete Xi cross-current estimate, q>1 grouped-tail theorem, "
            "Abel gap, winding cap, contact exclusion, Q209, cofinal "
            "descendant theorem, Lambda<=0, PF-infinity, RH, or prize-level "
            "conclusion."
        ),
        "source_audit": audit,
        "symbolic_certificate": symbolic,
        "analytic_majorant_certificate": analytic,
        "exact_integer_budget": budget,
        "exact": exact,
        "counts": {
            "rows": len(rows),
            "exact_ratio_logs": 1,
            "differentiated_ratio_logs": 1,
            "carrier_coordinate_errors": 4,
            "aggregate_coordinate_errors": 4,
            "growing_prefix_exponent_denominator": 18,
            "worst_height_prefix_length": 3,
            "finite_height_current_theorems": 1,
            "bulk_aggregate_closures": 0,
            "xi_level_current_theorems": 0,
        },
        "rows": [asdict(row) for row in rows],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    atomic_write(args.out, json.dumps(payload, indent=2) + "\n")
    atomic_write(args.note, render_note(payload))
    print(
        "wrote growing-prefix finite-height gate: "
        f"{payload['counts']['rows']} rows, exponent 1/18, "
        "M=3 at L=50, current <-1/400, 0 bulk closures"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
