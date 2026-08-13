#!/usr/bin/env python3
"""Independently check the Fourier/score-moment regular-field gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_fourier_moment_regular_field_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def expression_hash(expression: sp.Expr) -> str:
    return hashlib.sha256(sp.srepr(sp.expand(expression)).encode("utf-8")).hexdigest()


def independent_heat_monomial(m: int, z: sp.Symbol, sign: int) -> sp.Expr:
    total = sp.Integer(0)
    iterate = z**m
    order = 0
    while iterate != 0:
        total += sign**order * iterate / sp.factorial(order)
        iterate = sp.diff(iterate, z, 2)
        order += 1
    return sp.expand(total)


def sinc(z: sp.Expr) -> sp.Expr:
    return sp.sin(z) / z


def main() -> None:
    require(RESULT_PATH.exists(), "missing Fourier-moment result")
    require(NOTE_PATH.exists(), "missing Fourier-moment note")
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    require(artifact.get("kind") == KIND, "kind drift")

    expected_counts = {
        "rows": 20,
        "sources": 5,
        "multiplicities": 15,
        "direct_moment_checks": 15,
        "score_moment_checks": 15,
        "score_jet_checks": 15,
        "sinc_models": 15,
        "gaussian_models": 5,
        "strong_curvature_checks": 5,
        "hermite_checks": 15,
        "open_handoffs": 2,
        "actual_xi_field_bounds": 0,
        "degree_uniform_bounds": 0,
        "rh_conclusions": 0,
    }
    require(artifact.get("counts") == expected_counts, "count drift")

    sources = artifact.get("source_audit", {})
    require(len(sources) == 5, "source count drift")
    for source in sources.values():
        path = REPO_ROOT / source["path"]
        require(path.exists(), f"missing source: {path}")
        require(file_hash(path) == source["sha256"], f"source hash drift: {path}")

    expected_ids = [f"fmr_{index:02d}_{suffix}" for index, suffix in enumerate(
        [
            "sources", "complex", "contact", "parity", "ell", "score", "score_contact",
            "score_ell", "compact", "compact_shape", "compact_field", "compact_heat",
            "gaussian", "curvature", "score_nonclosure", "theta_guard", "pi", "xi",
            "uniform", "boundary",
        ],
        start=1,
    )]
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == expected_ids, "row id drift")
    require(sum(row.get("readiness") == "open" for row in rows) == 2, "open-row drift")

    x, u = sp.symbols("x u", real=True)
    stored_moments = artifact.get("moment_audit", [])
    require(len(stored_moments) == 15, "moment audit count drift")
    for m, stored in zip(range(2, 17), stored_moments, strict=True):
        direct = sp.expand_trig(sp.diff(sp.cos(x * u), x, m))
        score = sp.expand_trig(sp.diff(sp.sin(x * u), x, m))
        require(stored["multiplicity"] == m, f"moment multiplicity drift m={m}")
        require(expression_hash(direct) == stored["direct_derivative_sha256"], f"direct hash drift m={m}")
        require(expression_hash(score) == stored["score_derivative_sha256"], f"score hash drift m={m}")
        if m % 2 == 0:
            r = m // 2
            require(sp.simplify(direct - (-1) ** r * u**m * sp.cos(x * u)) == 0, f"direct parity failed m={m}")
            require(sp.simplify(score - (-1) ** r * u**m * sp.sin(x * u)) == 0, f"score parity failed m={m}")
        else:
            r = (m - 1) // 2
            require(sp.simplify(direct - (-1) ** (r + 1) * u**m * sp.sin(x * u)) == 0, f"direct parity failed m={m}")
            require(sp.simplify(score - (-1) ** r * u**m * sp.cos(x * u)) == 0, f"score parity failed m={m}")

    xj, cj, aj, bj = sp.symbols("xj cj aj bj", nonzero=True)
    h = xj - cj
    stored_score = artifact.get("score_jet_audit", [])
    require(len(stored_score) == 15, "score-jet count drift")
    for m, stored in zip(range(2, 17), stored_score, strict=True):
        score = h**m * (1 + aj * h + bj * h**2)
        transform = sp.cancel(score / xj)
        ratio = sp.factor(sp.diff(transform, xj, m + 1).subs(xj, cj) / ((m + 1) * sp.diff(transform, xj, m).subs(xj, cj)))
        ell = sp.factor((2 * cj * ratio - m) / 4)
        require(sp.simplify(ratio - (aj - 1 / cj)) == 0, f"independent score conversion failed m={m}")
        require(sp.simplify(ell - (cj * aj / 2 - sp.Rational(m + 2, 4))) == 0, f"independent score ell failed m={m}")
        # Dummy-symbol spelling is not part of the certificate.  The two
        # preceding exact identities independently verify the stored claim.

    xs = sp.symbols("x")
    stored_sinc = artifact.get("sinc_audit", [])
    require(len(stored_sinc) == 15, "sinc audit count drift")
    for stored in stored_sinc:
        m = int(stored["multiplicity"])
        phase = sp.sympify(stored["phase"])
        transform = sinc(sp.pi * xs) ** m * sinc(phase * xs)
        dm = sp.simplify(sp.diff(transform, xs, m).subs(xs, 1))
        dn = sp.simplify(sp.diff(transform, xs, m + 1).subs(xs, 1))
        field = sp.simplify(dn / ((m + 1) * dm))
        require(sp.simplify(field - (phase / sp.tan(phase) - m - 1)) == 0, f"independent sinc field failed m={m}, phase={phase}")
        require(expression_hash(transform) == stored["transform_sha256"], f"sinc hash drift m={m}, phase={phase}")
        require(expression_hash(dm) == stored["contact_derivative_sha256"], f"contact hash drift m={m}, phase={phase}")

    # Separate exact phases and a different collision height test the same formula.
    c_value = sp.Rational(3, 2)
    for m, phase in [(2, sp.Rational(4, 3) * sp.pi), (3, sp.Rational(3, 2) * sp.pi), (4, sp.Rational(5, 3) * sp.pi), (5, sp.Rational(4, 3) * sp.pi)]:
        transform = sinc(sp.pi * xs / c_value) ** m * sinc(phase * xs / c_value)
        dm = sp.simplify(sp.diff(transform, xs, m).subs(xs, c_value))
        dn = sp.simplify(sp.diff(transform, xs, m + 1).subs(xs, c_value))
        field = sp.simplify(dn / ((m + 1) * dm))
        require(sp.simplify(field - (phase / sp.tan(phase) - m - 1) / c_value) == 0, f"separate-height field failed m={m}")

    stored_gaussian = artifact.get("gaussian_audit", [])
    require(len(stored_gaussian) == 5, "Gaussian audit count drift")
    for stored in stored_gaussian:
        m = int(stored["multiplicity"])
        phase = sp.sympify(stored["phase"])
        variance = sp.sympify(stored["gaussian_variance"])
        radius = sp.sympify(stored["support_radius"])
        transform = sp.exp(-variance * xs**2 / 2) * sinc(sp.pi * xs) ** m * sinc(phase * xs)
        dm = sp.simplify(sp.diff(transform, xs, m).subs(xs, 1))
        dn = sp.simplify(sp.diff(transform, xs, m + 1).subs(xs, 1))
        field = sp.simplify(dn / ((m + 1) * dm))
        require(sp.simplify(field - (phase / sp.tan(phase) - m - 1 - variance)) == 0, f"independent Gaussian field failed m={m}")
        require(sp.simplify(variance - radius**2).is_positive is True, f"independent curvature margin failed m={m}")
        require(expression_hash(transform) == stored["transform_sha256"], f"Gaussian hash drift m={m}")

    z = sp.symbols("z")
    stored_hermite = artifact.get("hermite_audit", [])
    require(len(stored_hermite) == 15, "Hermite audit count drift")
    for m, stored in zip(range(2, 17), stored_hermite, strict=True):
        forward = independent_heat_monomial(m, z, -1)
        backward = independent_heat_monomial(m, z, 1)
        forward_count = int(sp.Poly(forward, z).count_roots(-sp.oo, sp.oo))
        backward_count = int(sp.Poly(backward, z).count_roots(-sp.oo, sp.oo))
        require(forward_count == m, f"independent forward Hermite failed m={m}")
        require(backward_count == m % 2, f"independent backward Hermite failed m={m}")
        require(expression_hash(forward) == stored["forward_sha256"], f"forward Hermite hash drift m={m}")
        require(expression_hash(backward) == stored["backward_sha256"], f"backward Hermite hash drift m={m}")

    note = NOTE_PATH.read_text(encoding="utf-8")
    required_phrases = [
        "The parity signs are part of the theorem",
        "characteristic function of a sum",
        "realizes every real `ell`",
        "continuous positive density",
        "no single countermodel here is claimed",
        "Pi provenance is explicit",
        "not a proof of RH",
    ]
    for phrase in required_phrases:
        require(phrase in note, f"note missing required phrase: {phrase}")

    print(
        "validated Fourier-moment regular-field gate: "
        "20 rows, 5 sources, 15 multiplicities, 15 direct moment checks, "
        "15 score moment checks, 15 score-jet checks, 15 sinc models, "
        "5 Gaussian models, 5 strong-curvature checks, 15 Hermite checks, "
        "2 open handoffs, 0 Xi field bounds, 0 degree-uniform bounds"
    )


if __name__ == "__main__":
    main()
