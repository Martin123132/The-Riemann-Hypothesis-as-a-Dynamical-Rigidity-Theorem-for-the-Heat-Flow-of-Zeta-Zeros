#!/usr/bin/env python3
"""Independently check the scaled Airy-branch amplitude envelope."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))

from flint import arb, acb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_scaled_airy_branch_amplitude_envelope_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def envelopes(x: arb) -> tuple[arb, arb, arb]:
    i = acb(0, 1)
    ai, ai_prime, bi, bi_prime = acb(-x).airy()
    w = ai - i * bi
    wx = -ai_prime + i * bi_prime
    return abs(w), abs(wx), abs(wx - i * x.sqrt() * w)


def independent_atlas(xmin: arb, xmax: arb) -> tuple[arb, arb, arb]:
    w_max = arb(0)
    low_wx_max = arb(0)
    rho_max = arb(0)

    def unit_panel(index: int) -> arb:
        return arb(arb(2 * index + 1) / 2, arb(1) / 2)

    log_panels = 9_000
    log_min = xmin.log()
    for index in range(log_panels):
        u = log_min * (1 - unit_panel(index) / log_panels)
        w_value, wx_value, rho_value = envelopes(u.exp())
        if w_value.upper() > w_max.upper(): w_max = w_value
        if wx_value.upper() > low_wx_max.upper(): low_wx_max = wx_value
        if rho_value.upper() > rho_max.upper(): rho_max = rho_value
    phase_panels = 135_000
    xi0 = arb(2) / 3
    xi1 = arb(2) * xmax ** (arb(3) / 2) / 3
    for index in range(phase_panels):
        xi = xi0 + (xi1 - xi0) * unit_panel(index) / phase_panels
        x = (arb(3) * xi / 2) ** (arb(2) / 3)
        w_value, _, rho_value = envelopes(x)
        if w_value.upper() > w_max.upper(): w_max = w_value
        if rho_value.upper() > rho_max.upper(): rho_max = rho_value
    return w_max, low_wx_max, rho_max


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

    ctx.dps = 130
    ctx.threads = 1
    pi = arb.pi()
    c0 = arb(159_577)
    beta = (pi * c0**2 / 8) ** (arb(1) / 3)
    eps = beta**-2
    xmin = pi / (16 * beta)
    xmax = (arb(797) * beta / c0) ** 2 + xmin + 64
    m_bound, n_low_bound, rho_bound = independent_atlas(xmin, xmax)
    n_squared = xmax.sqrt()/pi*(1+arb(7)/(32*xmax**3))
    require(n_squared > (1+arb(7)/32)/pi, "independent high-X endpoint order failed")
    n_bound = n_squared.sqrt()
    require(n_bound > n_low_bound, "independent low/high derivative envelope mismatch")

    # Independently expand the stored rational weights using a coarser
    # lambda envelope than the exact event maximum.
    l = arb(116); y = arb(64)
    ud = eps * (13*l+3*y)/60 + eps**2 * (448*l**5+280*l**2*y**3+4565*l**2+105*l*y**4+30*l*y+63*y**5+405*y**2)/50400
    vb = eps * (8*l**2+4*l*y+3*y**2)/60 + eps**2 * (40*l**3+20*l**2*y+l*y**2+9*y**3+27)/1680
    uy = eps*arb(3)/60 + eps**2*(840*l**2*y**2+420*l*y**3+30*l+315*y**4+810*y)/50400
    vy = eps*(4*l+6*y)/60 + eps**2*(20*l**2+2*l*y+27*y**2)/1680
    s_bound = ((1+ud)*m_bound+vb*n_bound)/2
    sy_bound = (uy*m_bound+vy*n_bound+(1+ud+xmax.sqrt()*vb)*rho_bound)/2
    require(s_bound < arb("0.356") and sy_bound < arb("0.258"), "independent amplitude bounds failed")

    stored = artifact["certificate"]
    for key, value in (
        ("W_modulus_bound_ball", m_bound),
        ("W_X_modulus_bound_ball", n_bound),
        ("extracted_phase_residual_bound_ball", rho_bound),
        ("scaled_branch_amplitude_bound_ball", s_bound),
        ("scaled_branch_y_derivative_bound_ball", sy_bound),
    ):
        require(arb(stored[key]).overlaps(value), f"stored envelope drift: {key}")

    decision = artifact["decision"]
    require(decision["exact_scaled_branch_amplitude_bound_proved"] is True, "amplitude bound flag lost")
    require(decision["exact_scaled_branch_first_y_derivative_bound_proved"] is True, "derivative bound flag lost")
    require(decision["ordinary_logistic_Gamma_amplitude_identification_proved"] is False, "ordinary bridge overclaim")
    require(decision["opposite_branch_corridor_integral_bound_proved"] is False, "integral overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    note = NOTE.read_text(encoding="utf-8")
    for token in ("There is no asymptotic replacement", "DLMF 9.8.21", "|S_sigma|", "does not yet compare", "PF-infinity, RH"):
        require(token in note, f"note boundary token missing: {token}")
    print("independently checked scaled Airy amplitude envelopes: endpoint moduli, rational weights, |S|<0.356, |S_y|<0.258", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
