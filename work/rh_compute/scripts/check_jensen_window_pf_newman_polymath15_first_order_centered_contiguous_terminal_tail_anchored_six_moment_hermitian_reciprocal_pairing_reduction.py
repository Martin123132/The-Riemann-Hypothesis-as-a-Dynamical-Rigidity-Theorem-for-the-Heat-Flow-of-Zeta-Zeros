#!/usr/bin/env python3
"""Validate the Hermitian reciprocal swap and near-diagonal reduction."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = REPO_ROOT / "work/rh_compute/scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_hermitian_reciprocal_pairing_reduction as gate  # noqa: E402


STEM = gate.STEM
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    f"hrpr_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "domain",
            "covariance",
            "pair",
            "derivative",
            "diagonal",
            "split",
            "coordinates",
            "stationary",
            "leading",
            "offset",
            "log",
            "lattice",
            "mixed",
            "cutoff",
            "guards",
            "handoff",
            "boundary",
        ),
        start=1,
    )
]
EXPECTED_COUNTS = {
    "rows": 17,
    "hermitian_swap_pair_identities": 1,
    "exact_full_diagonal_nulls": 1,
    "leading_real_sine_channels": 1,
    "reciprocal_lattice_defect_bounds": 1,
    "unpaired_interior_dual_main_modes": 0,
    "sign_guards": 2,
    "imaginary_resonance_guards": 1,
    "poisson_remainder_bounds": 0,
    "hard_cutoff_boundary_bounds": 0,
    "imaginary_correction_bounds": 0,
    "signed_offset_sum_bounds": 0,
    "signed_flow_bounds": 0,
    "phi_b_bounds": 0,
}
SUMMARY = (
    "validated Hermitian reciprocal pairing reduction: 17 rows, "
    "1 swap identity, 1 lattice-defect bound, "
    "0 Poisson remainder bounds, 0 signed flow bounds"
)


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path, issues: list[str]) -> dict:
    if not path.is_file():
        issues.append(f"missing result: {path}")
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"invalid result: {exc}")
        return {}


def require_zero(expression: sp.Expr, label: str, issues: list[str]) -> None:
    if sp.simplify(sp.expand_trig(sp.expand_complex(expression))) != 0:
        issues.append(label)


def independent_source_audit(payload: dict, issues: list[str]) -> None:
    saved = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in gate.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
        elif saved.get(key) != file_hash(path):
            issues.append(f"source hash drifted: {key}")
    try:
        sources = {
            key: json.loads(path.read_text(encoding="utf-8"))
            for key, path in gate.SOURCE_PATHS.items()
        }
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"independent source load failed: {exc}")
        return
    if sources["flow_matrix_phase"].get("counts", {}).get(
        "hermitian_reciprocal_diagonal_nulls"
    ) != 1:
        issues.append("independent diagonal-null source audit failed")
    if sources["saddle_flow"].get("counts", {}).get("max_moment_order") != 5:
        issues.append("independent saddle-flow source audit failed")
    if sources["two_carrier"].get("counts", {}).get("leading_q1_kernels") != 2:
        issues.append("independent leading-kernel source audit failed")


def independent_swap_audit(issues: list[str]) -> None:
    alpha, mu, real_part, imag_part = sp.symbols("a m R I", real=True)
    angle = 2 * sp.pi * alpha * mu
    z = (real_part + sp.I * imag_part) * sp.exp(sp.I * angle)
    pair = sp.I * mu * z - sp.I * mu * sp.conjugate(z)
    expected = -2 * mu * (
        real_part * sp.sin(angle) + imag_part * sp.cos(angle)
    )
    require_zero(pair - expected, "independent swap pairing failed", issues)
    undifferentiated = z + sp.conjugate(z)
    require_zero(
        sp.diff(undifferentiated, alpha) / (2 * sp.pi) - expected,
        "independent phase derivative failed",
        issues,
    )
    require_zero(pair.subs(mu, 0), "independent diagonal null failed", issues)


def independent_reciprocal_audit(issues: list[str]) -> None:
    a, k, ell = sp.symbols("a k l", positive=True)
    alpha = a**2
    n_star, m_star = alpha / k, alpha / ell
    mu = sp.expand_log(sp.log(ell / k), force=True)
    tau = sp.expand_log(sp.log(k * ell / alpha), force=True)
    u_n = sp.expand_log(sp.log(a / n_star), force=True)
    u_m = sp.expand_log(sp.log(a / m_star), force=True)
    require_zero(u_n - u_m + mu, "independent u difference failed", issues)
    require_zero(u_n + u_m - tau, "independent u sum failed", issues)
    hessian = sp.diag(-k**2 / alpha, ell**2 / alpha)
    factor = alpha / (k * ell)
    require_zero(
        factor**2 * (-hessian.det()) - 1,
        "independent stationary factor failed",
        issues,
    )


def independent_offset_audit(issues: list[str]) -> None:
    alpha, k, d = sp.symbols("alpha k d", positive=True)
    phase = alpha * sp.log(1 + d / k)
    require_zero(
        sp.diff(phase, d) - alpha / (k + d),
        "independent offset first derivative failed",
        issues,
    )
    require_zero(
        sp.diff(phase, d, 2) + alpha / (k + d) ** 2,
        "independent offset second derivative failed",
        issues,
    )
    x, t = sp.symbols("x t", positive=True)
    require_zero(
        x - sp.log(1 + x) - sp.integrate(t / (1 + t), (t, 0, x)),
        "independent log-defect integral failed",
        issues,
    )


def finite_lattice_defect_audit(issues: list[str]) -> None:
    mp.mp.dps = 80
    cases = (
        (mp.mpf(101) / 3, 7, 1),
        (mp.mpf(157) / 5, 11, 3),
        (mp.mpf(1000) + mp.mpf(1) / 7, 31, 5),
        (mp.mpf(9973) / 11, 83, 9),
    )
    for alpha, k, d in cases:
        reciprocal = alpha / k
        rho = abs(reciprocal - mp.nint(reciprocal))
        phase = alpha * mp.log(1 + mp.mpf(d) / k)
        distance = abs(phase - mp.nint(phase))
        bound = d * rho + alpha * d**2 / (2 * k**2)
        if distance > bound + mp.mpf("1e-70"):
            issues.append(f"finite lattice-defect bound failed: {alpha},{k},{d}")
        if abs(mp.sin(2 * mp.pi * phase)) > 2 * mp.pi * distance + mp.mpf("1e-70"):
            issues.append(f"finite sine-distance bound failed: {alpha},{k},{d}")


def independent_guard_audit(issues: list[str]) -> None:
    def pair(alpha: sp.Expr, mu: sp.Expr, real: sp.Expr, imag: sp.Expr) -> sp.Expr:
        return sp.simplify(
            -mu
            * (
                real * sp.sin(2 * sp.pi * alpha * mu)
                + imag * sp.cos(2 * sp.pi * alpha * mu)
            )
        )

    if pair(1, sp.Rational(1, 4), 1, 0) != -sp.Rational(1, 4):
        issues.append("independent negative sign guard failed")
    if pair(3, sp.Rational(1, 4), 1, 0) != sp.Rational(1, 4):
        issues.append("independent positive sign guard failed")
    if pair(2, sp.Rational(1, 2), 0, 1) != -sp.Rational(1, 2):
        issues.append("independent imaginary resonance guard failed")


def structural_audit(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("kind drifted")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")
    rows = payload.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids drifted")
    if len({row.get("id") for row in rows}) != len(rows):
        issues.append("duplicate row ids")
    for key in (
        "poisson_remainder_bounds",
        "hard_cutoff_boundary_bounds",
        "imaginary_correction_bounds",
        "signed_offset_sum_bounds",
        "signed_flow_bounds",
        "phi_b_bounds",
    ):
        if payload.get("counts", {}).get(key) != 0:
            issues.append(f"{key} was overpromoted")
    boundary = payload.get("proof_boundary", "")
    for marker in (
        "no discrete Poisson or B-process formula",
        "imaginary-correction bound",
        "signed flow estimate",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")


def note_audit(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for marker in (
        "# Hermitian Reciprocal Pairing Reduction",
        "P_H(k,l)=-mu",
        "C_(k,l)*mu*sin(2*pi*a^2*mu)",
        "dist(alpha*log(1+d/k),Z)",
        "imaginary correction",
        "0 Poisson remainder bounds",
        "this is not a proof",
        f"python work/rh_compute/scripts/check_{STEM}.py",
    ):
        if marker not in text:
            issues.append(f"note marker missing: {marker}")


def main() -> int:
    issues: list[str] = []
    payload = load_json(RESULT, issues)
    if payload:
        independent_source_audit(payload, issues)
        independent_swap_audit(issues)
        independent_reciprocal_audit(issues)
        independent_offset_audit(issues)
        finite_lattice_defect_audit(issues)
        independent_guard_audit(issues)
        structural_audit(payload, issues)
    note_audit(issues)
    if issues:
        for issue in issues:
            print(f"FAIL: {issue}")
        return 1
    print(SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
