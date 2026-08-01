#!/usr/bin/env python3
"""Build a rigorous Arb pilot for quadratic modular-tail integrals."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path
import sys

try:
    import psutil

    psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
except Exception:
    pass

REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR_DIR = REPO_ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR_DIR) not in sys.path:
    sys.path.insert(0, str(VENDOR_DIR))

import flint
from flint import acb, arb
import sympy as sp


STEM = "jensen_window_pf_newman_theta_modular_tail_arb_quadratic_pilot"
DEFAULT_OUT = REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
ENVELOPE_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_modular_tail_derivative_envelope_gate.json"
)
PRECISION_BITS = 96
DERIVATIVE_ORDER = 9
N_START = 7
N_STOP = 16
T_CAP = "0.2"
U_STOP = "2.2"
INTEGRATION_OPTIONS = {
    "deg_limit": 60,
    "eval_limit": 50_000,
    "depth_limit": 25,
    "use_heap": True,
}


@dataclass(frozen=True)
class PilotRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def build_derivative_function():
    u, n = sp.symbols("u n", positive=True, real=True)

    def phi(z: sp.Expr) -> sp.Expr:
        return (
            2 * sp.pi**2 * n**4 * sp.exp(9 * z)
            - 3 * sp.pi * n**2 * sp.exp(5 * z)
        ) * sp.exp(-sp.pi * n**2 * sp.exp(4 * z))

    omega = (1 + sp.erf(3 * sp.sinh(4 * u))) / 2
    block = omega * phi(u) + (1 - omega) * phi(-u)
    derivative = sp.diff(block, u, DERIVATIVE_ORDER)
    functions = {
        "exp": lambda z: z.exp(),
        "sinh": lambda z: z.sinh(),
        "cosh": lambda z: z.cosh(),
        "erf": lambda z: z.erf(),
        "sqrt": lambda z: z.sqrt(),
        "pi": acb.pi(),
    }
    return sp.lambdify(
        (u, n),
        derivative,
        modules=[functions],
        cse=True,
        docstring_limit=0,
    )


def compute_pilot() -> dict:
    flint.ctx.prec = PRECISION_BITS
    evaluate = build_derivative_function()

    def finite_remainder(z: acb) -> acb:
        return sum(
            (evaluate(z, index) for index in range(N_START, N_STOP + 1)),
            acb(0),
        )

    midpoint_value = finite_remainder(acb("0.5"))
    if not midpoint_value.is_finite():
        raise RuntimeError("finite-tail midpoint evaluation was not finite")

    t_cap = acb(T_CAP)

    def integrand(z: acb, _analytic: bool) -> acb:
        value = finite_remainder(z)
        return (t_cap * z * z).exp() * value * value

    result = acb.integral(
        integrand,
        arb(0),
        arb(U_STOP),
        **INTEGRATION_OPTIONS,
    )
    if not result.is_finite():
        raise RuntimeError("Arb integration returned a non-finite ball")
    if result.real.lower() <= 0:
        raise RuntimeError("quadratic integral was not certified positive")
    if not result.imag.contains(0):
        raise RuntimeError("quadratic integral imaginary ball excludes zero")

    return {
        "quantity": (
            "integral_0^(11/5) exp(u^2/5)*"
            "[D^9 sum_(n=7)^16 b_n(u)]^2 du"
        ),
        "real_enclosure": result.real.str(35, more=True),
        "real_lower": result.real.lower().str(35),
        "real_upper": result.real.upper().str(35),
        "imaginary_enclosure": result.imag.str(20, more=True),
        "relative_accuracy_bits": result.real.rel_accuracy_bits(),
        "precision_bits": PRECISION_BITS,
        "derivative_order": DERIVATIVE_ORDER,
        "n_start": N_START,
        "n_stop": N_STOP,
        "t_cap_exact": "1/5",
        "u_interval_exact": ["0", "11/5"],
        "integration_options": INTEGRATION_OPTIONS,
        "integrand_entire": True,
    }


def build_artifact() -> dict:
    envelope = json.loads(ENVELOPE_SOURCE.read_text(encoding="utf-8"))
    if envelope.get("kind") != (
        "jensen_window_pf_newman_theta_modular_tail_derivative_envelope_gate"
    ):
        raise RuntimeError("unexpected derivative-envelope source")
    pilot = compute_pilot()
    rows = [
        PilotRow(
            id="ntmtaqp_01_weighted_cauchy_schwarz",
            role="exact_lemma",
            readiness="proved",
            claim=(
                "Every nonnegative weighted L1 tail term can be reduced "
                "exactly to a weight mass and an analytic quadratic integral."
            ),
            formula=(
                "int_I W|f| <= sqrt(int_I W)*sqrt(int_I W*f^2), W>=0"
            ),
            proof_boundary=(
                "This is Cauchy-Schwarz applied to sqrt(W) and "
                "sqrt(W)|f|; it does not estimate either factor."
            ),
        ),
        PilotRow(
            id="ntmtaqp_02_heat_envelope_application",
            role="exact_reduction",
            readiness="proved",
            claim=(
                "The positive heat-polynomial weights in E_(N,k) admit "
                "the weighted quadratic reduction term by term."
            ),
            formula=(
                "W_b(u)=exp(Tu^2)binom(k,b)P_b(T,u)>=0 for u>=0"
            ),
            proof_boundary=(
                "The reduction leaves certified bounds for each weight mass "
                "and exact infinite-tail quadratic integral to be supplied."
            ),
        ),
        PilotRow(
            id="ntmtaqp_03_arb_finite_compact_pilot",
            role="rigorous_numerical_diagnostic",
            readiness="certified_for_stated_finite_integral",
            claim=(
                "Arb certifies one order-9 finite-tail compact quadratic "
                "integral with a narrow positive enclosure."
            ),
            formula=pilot["quantity"],
            proof_boundary=(
                "The enclosure covers only n=7..16 and u in [0,11/5]."
            ),
            diagnostics=pilot,
        ),
        PilotRow(
            id="ntmtaqp_04_omitted_tail_guard",
            role="non_promotion_guard",
            readiness="open",
            claim=(
                "The pilot is not an enclosure for the exact r_6 tail because "
                "n>=17 and u>11/5 are not yet included."
            ),
            formula="r_6=sum_(n=7)^16 b_n + sum_(n>=17)b_n",
            proof_boundary=(
                "Both the arithmetic remainder and outer-u integral require "
                "independent analytic bounds before any C1 budget is certified."
            ),
        ),
        PilotRow(
            id="ntmtaqp_05_matrix_handoff",
            role="proof_search_target",
            readiness="not_ready_to_apply",
            claim=(
                "Complete the finite (b,k,N,t) quadratic matrix, append "
                "analytic n and u tails, and compare the resulting d0,d1 "
                "balls with retained J and J' lower balls."
            ),
            formula=(
                "certify d_(N,j,m)(T) < retained first-jet separation margin"
            ),
            proof_boundary=(
                "No retained lower separation, strict Laguerre conclusion, "
                "Lambda<=0, RH, or Clay-prize conclusion follows here."
            ),
        ),
    ]
    return {
        "kind": STEM,
        "date": "2026-07-24",
        "status": "rigorous_finite_compact_pilot_only",
        "source": {
            "derivative_envelope": str(
                ENVELOPE_SOURCE.relative_to(REPO_ROOT)
            ).replace("\\", "/"),
            "derivative_envelope_sha256": file_hash(ENVELOPE_SOURCE),
            "derivative_envelope_rows": len(envelope.get("rows", [])),
        },
        "weighted_reduction": {
            "inequality": (
                "int_I W|f| <= "
                "(int_I W)^(1/2)(int_I W*f^2)^(1/2)"
            ),
            "hypotheses": ["I measurable", "W>=0", "f real on I"],
            "application": (
                "W_b=exp(Tu^2)binom(k,b)P_b(T,u), "
                "f=r_N^(k-b)"
            ),
            "analytic_advantage": (
                "W_b*f^2 is analytic for finite theta sums, avoiding "
                "nonanalytic absolute values inside Arb integration."
            ),
        },
        "pilot": pilot,
        "rows": [asdict(row) for row in rows],
        "proof_boundary": (
            "This artifact proves the weighted Cauchy-Schwarz reduction and "
            "certifies one explicitly stated finite compact integral. It does "
            "not cover the infinite arithmetic tail, the outer-u tail, the "
            "full derivative matrix, retained first-jet lower separation, "
            "Lambda<=0, RH, or a Clay-prize proof."
        ),
        "versions": {
            "python_flint": getattr(flint, "__version__", "unknown"),
            "sympy": sp.__version__,
        },
    }


def render_note(artifact: dict) -> str:
    pilot = artifact["pilot"]
    return "\n".join(
        [
            "# Newman Theta Modular-Tail Arb Quadratic Pilot",
            "",
            "Date: 2026-07-24",
            "",
            "Status: exact weighted reduction plus one rigorous finite compact",
            "Arb integral. This is not an exact infinite-tail bound and not a",
            "proof of `Lambda<=0`, RH, or a Clay-prize result.",
            "",
            "```text",
            "work/rh_compute/results/"
            "jensen_window_pf_newman_theta_modular_tail_arb_quadratic_pilot.json",
            "python work/rh_compute/scripts/"
            "jensen_window_pf_newman_theta_modular_tail_arb_quadratic_pilot.py",
            "python work/rh_compute/scripts/"
            "check_jensen_window_pf_newman_theta_modular_tail_arb_quadratic_pilot.py",
            "```",
            "",
            "## Exact Weighted Reduction",
            "",
            "For a real function `f` and nonnegative weight `W`,",
            "",
            "```text",
            artifact["weighted_reduction"]["inequality"],
            "```",
            "",
            "This is exactly Cauchy-Schwarz. Applied term by term to the",
            "positive heat-envelope weights",
            "`W_b=exp(Tu^2)binom(k,b)P_b(T,u)`, it replaces each nonsmooth",
            "absolute-value integral by a weight mass and an analytic",
            "quadratic integral.",
            "",
            "## Rigorous Finite Compact Pilot",
            "",
            "At 96-bit Arb precision, the certified quantity is",
            "",
            "```text",
            pilot["quantity"],
            f"in {pilot['real_enclosure']}",
            "```",
            "",
            f"The enclosure has {pilot['relative_accuracy_bits']} relative",
            "accuracy bits. The integrand is entire and the Arb callback",
            "therefore has no branch-cut qualification.",
            "",
            "## Non-Promotion Guard",
            "",
            "This is only the finite arithmetic block `n=7..16` on",
            "`u in [0,11/5]`. It omits `n>=17`, `u>11/5`, the remaining",
            "heat-polynomial weights and derivative orders, and all retained",
            "first-jet lower bounds.",
            "",
            "## Next Certificate",
            "",
            "Build the finite `(b,k,N,t)` quadratic matrix, attach analytic",
            "bounds for the omitted arithmetic and outer-u tails, and compare",
            "the resulting `d0,d1` balls with rigorous retained `J,J'` lower",
            "balls.",
            "",
            "No strict Laguerre conclusion, `Lambda<=0`, RH, or Clay-prize",
            "conclusion is claimed.",
            "",
        ]
    )


def write_artifact(artifact: dict, out: Path, note: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    note.write_text(render_note(artifact), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    write_artifact(artifact, args.out, args.note)
    print(
        "wrote Newman theta modular-tail Arb quadratic pilot: "
        f"{len(artifact['rows'])} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
