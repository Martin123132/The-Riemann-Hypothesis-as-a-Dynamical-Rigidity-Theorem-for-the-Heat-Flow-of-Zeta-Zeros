#!/usr/bin/env python3
"""Independently check the beta^-4 canonical Airy reduction."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))

import mpmath as mp
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_canonical_airy_derivative_reduction_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(RESULT.is_file() and NOTE.is_file(), "missing result or note")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == STEM and artifact["passed"] is True, "gate identity drift")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file() and file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file() and file_hash(path) == record["sha256"], f"source hash drift: {path}")

    z, y, x, lam, eps = sp.symbols("z y X lambda epsilon", real=True)
    i = sp.I
    a1 = y / 2 - 3 * z**2 / 4
    r1 = -2 * z**5 / 15 + z**3 * y / 3 - z * y**2 / 4
    a2 = 13 * z**4 / 32 - 3 * z**2 * y / 8
    r2 = 17 * z**7 / 315 - 2 * z**5 * y / 15 + z**3 * y**2 / 12
    multiplier = sp.Poly(sp.expand(1 + eps * (a1 + i * r1) + eps**2 * (a2 + i * r2 + i * a1 * r1 - r1**2 / 2)), z)
    require(multiplier.degree() == 10, "independent multiplier degree drift")

    # Explicit derivatives from A''=-X A, independently frozen through k=10.
    derivative_pairs = [
        (1, 0), (0, 1), (-x, 0), (-1, -x), (x**2, -2),
        (4*x, x**2), (4-x**3, 6*x), (-9*x**2, 10-x**3),
        (x**4-28*x, -12*x**2), (16*x**3-28, x**4-52*x),
        (-x**5+100*x**2, 20*x**3-80),
    ]
    u = sp.Integer(0)
    v = sp.Integer(0)
    for k, (p_k, q_k) in enumerate(derivative_pairs):
        coefficient = multiplier.coeff_monomial(z**k)
        u += coefficient * i**k * p_k
        v += coefficient * i**k * q_k

    expected_lambda = {
        "U1": -(13 * lam + 3 * y) / 60,
        "V1": (8 * lam**2 - 4 * lam * y + 3 * y**2) / 60,
        "U2": -(448 * lam**5 + 280 * lam**2 * y**3 + 4565 * lam**2 - 105 * lam * y**4 + 30 * lam * y + 63 * y**5 - 405 * y**2) / 50400,
        "V2": (40 * lam**3 - 20 * lam**2 * y + lam * y**2 - 9 * y**3 + 27) / 1680,
    }
    stored = {
        name: sp.sympify(value.replace("lambda", "lam"), locals={"lam": lam, "y": y})
        for name, value in artifact["reduced_weights_lambda_y"].items()
    }
    for name, expected in expected_lambda.items():
        require(sp.expand(stored[name] - expected) == 0, f"stored {name} drift")
    u_expected = 1 + eps * expected_lambda["U1"] + eps**2 * expected_lambda["U2"]
    v_expected = eps * expected_lambda["V1"] + eps**2 * expected_lambda["V2"]
    require(sp.expand(u.subs(x, lam + y) - u_expected) == 0, "independent U contraction failed")
    require(sp.expand(v.subs(x, lam + y) - v_expected) == 0, "independent V contraction failed")

    mp.mp.dps = 70
    for sample in (mp.mpf("0.7"), mp.mpf("5.25"), mp.mpf("41.0")):
        airy = lambda value: mp.airyai(-value)
        a0 = airy(sample)
        a_x = mp.diff(airy, sample, 1)
        for k, (p_k, q_k) in enumerate(derivative_pairs):
            direct = mp.diff(airy, sample, k)
            exact_sample = sp.Rational(str(sample))
            p_value = mp.mpf(str(sp.N(sp.sympify(p_k).subs(x, exact_sample), 70)))
            q_value = mp.mpf(str(sp.N(sp.sympify(q_k).subs(x, exact_sample), 70)))
            reduced = p_value * a0 + q_value * a_x
            require(abs(direct - reduced) < mp.mpf("1e-50") * max(1, abs(direct)), f"numeric Airy derivative drift at k={k}")

    decision = artifact["decision"]
    require(decision["beta_minus_4_two_dimensional_model_reduced_exactly"] is True, "exact reduction flag lost")
    require(decision["ordinary_carrier_branch_match_proved"] is False, "ordinary bridge overclaim")
    require(decision["complete_Q_K_minus_T_bound_proved"] is False, "endpoint residual overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    note = NOTE.read_text(encoding="utf-8")
    for token in ("symmetric Abel limits", "2*pi", "only two scalar carrier channels", "does not perform that asymptotic match", "no RH"):
        require(token in note, f"note boundary token missing: {token}")
    print(
        "independently checked beta^-4 Airy reduction: degree 10, two carrier channels, 3 numeric ODE samples, 0 ordinary-carrier claims",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
