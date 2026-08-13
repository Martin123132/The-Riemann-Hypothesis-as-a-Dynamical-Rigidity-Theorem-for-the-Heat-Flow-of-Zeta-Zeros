#!/usr/bin/env python3
"""Separate the exact endpoint roster identity from the source transition approximation."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import struct
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import mpmath as mp

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_upper_endpoint_selector_boundary_stability_gate as endpoint_gate


SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
PAPER_2015 = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_2015_1502.06903.pdf"
PAPER_2026 = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_Brereton_2607.15310.pdf"
ENDPOINT_GATE = endpoint_gate.RESULT
SELECTOR_GATE = endpoint_gate.SELECTOR_GATE
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_source_transition_endpoint_scheduling_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_source_transition_endpoint_scheduling_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 110
C = 159_577
B_MINUS = 5_122_421
C_PLUS = C + 2
B_PLUS = B_MINUS + 2
ROOT_CENTER = endpoint_gate.LOWER_TRANSITION_ROOT_CENTER
ROOT_RADIUS = endpoint_gate.ROOT_BRACKET_RADIUS
RS_CUTOFF = 621


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def float32(value: float) -> float:
    return struct.unpack(">f", struct.pack(">f", value))[0]


def threshold_mpf() -> mp.mpf:
    numerator, denominator = float32(3.2).as_integer_ratio()
    return mp.mpf(numerator) / denominator


def source_contract() -> dict[str, Any]:
    text = SOURCE.read_text(encoding="utf-8")
    required = (
        "if (g.lt.3.2) then",
        "M2=M2+2",
        "call gauleg(bot,top,abb,wee)",
        "zsum=(M2-2)*ri/(sqrt(2*a)*(2**(1.5)))",
        "do ial=M2,N1,2",
        "zsum(i)=zsum(i)*sqrt(8/aar(i))",
        "RN1=RN1startofblock + (2*(real(MT,dp)+1.0)*real(aenums,dp))",
        "raecutoff=RN1",
        "CE=raecutoff*(1-sqrt(1-(a/(raecutoff*C))**2))/CE",
        "NC=aint(CE,dp1)",
    )
    missing = [snippet for snippet in required if snippet not in text]
    require(not missing, f"source transition contract drift: {missing}")
    return {
        "required_snippet_count": len(required),
        "all_required_snippets_present": True,
        "source_default_real_3_2_hex": float32(3.2).hex(),
        "source_default_real_3_2": str(threshold_mpf()),
    }


def event_interval_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    center = arb(ROOT_CENTER)
    radius = arb(ROOT_RADIUS)
    height = center + arb(0, radius)
    threshold = endpoint_gate.source_default_real(3.2)
    left = endpoint_gate.transition_residual(center - radius, C, threshold)
    right = endpoint_gate.transition_residual(center + radius, C, threshold)
    require(left.lower() > 0 and right.upper() < 0, "transition root bracket lost its signs")

    a = (8 * height / arb.pi()).sqrt()
    t6 = (height.log() / 6).exp()
    g = (C - a) * t6
    require((g - threshold).contains(0), "transition threshold is not enclosed")

    def rs_root(endpoint: int) -> arb:
        return (endpoint - (arb(endpoint) ** 2 - a**2).sqrt()) / 4

    n_minus_old = rs_root(B_MINUS)
    n_minus_new = rs_root(B_PLUS)
    for value in (n_minus_old, n_minus_new):
        require(value.lower() > RS_CUTOFF and value.upper() < RS_CUTOFF + 1, "RS cutoff floor changed")

    def direct_term(alpha: int) -> arb:
        s = arb(alpha) / a
        b = (s**2 - 1).sqrt()
        phase = height * (b * (b - s) + (s + b).log()) + height + arb.pi() / 8
        return (8 / a).sqrt() * phase.cos() / b.sqrt()

    lower_direct = direct_term(C)
    added_upper_direct = direct_term(B_PLUS)
    require(lower_direct.upper() < 0, "lower direct diagnostic lost its sign")
    require(added_upper_direct.lower() > 0, "upper direct diagnostic lost its sign")

    pre_count = (B_MINUS - C) // 2 + 1
    post_count = (B_PLUS - C_PLUS) // 2 + 1
    fixed_count = (B_PLUS - C) // 2 + 1
    require(pre_count == post_count and fixed_count == pre_count + 1, "roster counts do not match")

    return {
        "transition_entry_height_ball": height.str(PRECISION, more=True),
        "transition_residual_left_ball": left.str(28, more=True),
        "transition_residual_right_ball": right.str(28, more=True),
        "threshold_contact_ball": (g - threshold).str(28, more=True),
        "pre_effective_start": C,
        "pre_endpoint": B_MINUS,
        "post_effective_start": C_PLUS,
        "post_endpoint": B_PLUS,
        "pre_roster_count": pre_count,
        "post_roster_count": post_count,
        "fixed_chart_start": C,
        "fixed_chart_endpoint": B_PLUS,
        "fixed_chart_roster_count": fixed_count,
        "removed_lower_alpha": C,
        "added_upper_alpha": B_PLUS,
        "rs_lower_root_pre_ball": n_minus_old.str(35, more=True),
        "rs_lower_root_post_ball": n_minus_new.str(35, more=True),
        "rs_cutoff_pre": RS_CUTOFF,
        "rs_cutoff_post": RS_CUTOFF,
        "lower_direct_term_ball": lower_direct.str(35, more=True),
        "added_upper_direct_term_ball": added_upper_direct.str(35, more=True),
    }


def transition_quadrature(order: int) -> mp.mpf:
    """High-precision replica of the source's first 64-node transition branch."""
    mp.mp.dps = 90
    t = mp.mpf(ROOT_CENTER)
    p = mp.pi
    a = mp.sqrt(8 * t / p)
    t6 = t ** (mp.mpf(1) / 6)
    t3 = t6**2
    g = (C - a) * t6
    yphase = mp.fmod(t + p / 8, 2 * p)
    cc = (2 * p) ** mp.mpf("0.25")
    bot = -2 * cc * mp.sqrt(g) / t3
    aa = cc * (mp.sqrt(2 * g) * t3**2 / 4 - t3 * g * cc) / 2
    top = 4 / mp.sqrt(aa)
    wb = t3 * ((t3 / 4 - mp.mpf("1.5") * cc * mp.sqrt(2 * g)) * t3 + mp.mpf("6.375") * cc**2 * g) / 6
    bb = mp.e ** (1j * p / 4) * wb
    nodes, weights = mp.gauss_quadrature(order, "legendre")
    midpoint = (top + bot) / 2
    half_width = (top - bot) / 2
    integral = mp.mpc(0)
    for node, weight in zip(nodes, weights):
        u = midpoint + half_width * node
        integral += half_width * weight * mp.e ** (-(u**2) * (aa + u * bb))
    phase = yphase + mp.sqrt(2 * p) * g * t3 / 2 - mp.sqrt(p) * cc * g ** mp.mpf("1.5") / 3
    ri = 2 * mp.re(mp.e ** (1j * phase) * integral)
    source_default_power = mp.mpf(str(float32(2.0 ** 1.5)))
    raw = C * ri / (mp.sqrt(2 * a) * source_default_power)
    return raw * mp.sqrt(8 / a)


def direct_term_point(alpha: int) -> mp.mpf:
    mp.mp.dps = 90
    t = mp.mpf(ROOT_CENTER)
    p = mp.pi
    a = mp.sqrt(8 * t / p)
    s = mp.mpf(alpha) / a
    b = mp.sqrt(s**2 - 1)
    phase = t * (b * (b - s) + mp.log(s + b)) + t + p / 8
    return mp.sqrt(8 / a) * mp.cos(phase) / mp.sqrt(b)


def numerical_diagnostic() -> dict[str, str]:
    transition_64 = transition_quadrature(64)
    transition_96 = transition_quadrature(96)
    lower = direct_term_point(C)
    upper = direct_term_point(B_PLUS)
    defect = transition_64 - lower
    coupled = defect + upper
    require(abs(transition_64 - transition_96) < mp.mpf("1e-40"), "transition quadrature did not stabilize")
    require(abs(defect) > mp.mpf("0.03"), "transition defect diagnostic unexpectedly small")
    require(abs(coupled) > mp.mpf("0.03"), "coupled source jump diagnostic unexpectedly small")
    return {
        "classification": "high_precision_point_diagnostic_not_interval_proof",
        "transition_64": mp.nstr(transition_64, 55),
        "transition_96": mp.nstr(transition_96, 55),
        "direct_lower_C": mp.nstr(lower, 55),
        "direct_added_B_plus": mp.nstr(upper, 55),
        "transition_defect_T_minus_F_C": mp.nstr(defect, 55),
        "coupled_source_change": mp.nstr(coupled, 55),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    d = artifact["numerical_diagnostic"]
    return f"""# Source transition and endpoint scheduling gate

Date: 2026-08-12

Status: exact finite-roster decomposition and interval cutoff certificate;
the source transition quadrature remains an approximation and is not promoted
to a theorem-level identity

At the earlier source transition entry,

```text
t_tr={c['transition_entry_height_ball']}.                         (SE1)
```

the pinned implementation changes its ideal odd roster from

```text
R_-={{C,C+2,...,B_-}},      C={C}, B_-={B_MINUS},
R_+={{C+2,C+4,...,B_+}},             B_+={B_PLUS}.               (SE2)
```

Both rosters contain `{c['pre_roster_count']}` points.  For any assigned
single-alpha summand `F_t(alpha)`, finite cancellation gives the exact identity

```text
sum_(alpha in R_+) F_t(alpha)-sum_(alpha in R_-) F_t(alpha)
  =F_t(B_+)-F_t(C).                                               (SE3)
```

If `T_C(t)` denotes the source transition contribution, the complete idealized
source scheduling change is therefore

```text
Delta_source=[T_C-F_t(C)]+F_t(B_+).                              (SE4)
```

This is a decomposition, not a cancellation theorem.  The first bracket is
the transition approximation defect and the second term is the newly scheduled
upper endpoint.  Interval arithmetic also gives

```text
N_-(B_-)={c['rs_lower_root_pre_ball']},
N_-(B_+)={c['rs_lower_root_post_ball']},                          (SE5)
```

so the Riemann-Siegel floor is `{c['rs_cutoff_pre']}` on both sides.  No RS
term absorbs `F_t(B_+)`.

The 64-node source-form quadrature was independently repeated at higher
precision.  As a diagnostic only,

```text
T_C={d['transition_64']},
F_t(C)={d['direct_lower_C']},
F_t(B_+)={d['direct_added_B_plus']},
Delta_source={d['coupled_source_change']}.                        (SE6)
```

The 64- and 96-node values agree beyond 40 decimal places, but (SE6) is not an
interval enclosure of the complete compiled source and is not used as proof.

The primary-source status agrees with this separation.  The 2015 paper calls
the hybrid an approximation, describes direct numerical integration in the
transition region, and notes unproved numerical correlations.  The 2026 paper
states that its hybrid representation is not exact and again sends part of the
transition range to numerical integration.

For theorem-level height transport, fix the analytic endpoint at
`B_fix={B_PLUS}`.  The fixed chart `{C},{C + 2},...,{B_PLUS}` has
`{c['fixed_chart_roster_count']}` points and obeys

```text
source-pre  = fixed-chart-F_t(B_+),
source-post = fixed-chart+[T_C-F_t(C)].                           (SE7)
```

Thus the arbitrary implementation threshold and moving endpoint belong in an
explicit source-approximation ledger.  They must not define the exact
`C -> C+2` selector theorem.  The exact selector change remains the one-cell
Poisson identity at `t*=pi C^2/8`, with fixed `B_fix={B_PLUS}`.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

This gate does not bound the complete source approximation error, the
quantitative 399-mode chart, `T_upper`, `Lambda<=0`, RH, or a prize-level
conclusion.
"""


def main() -> None:
    started = time.perf_counter()
    priority = endpoint_gate.selector_gate.event_gate.window_gate.cell_gate.ode_gate.set_low_priority()
    for dependency in (SOURCE, PAPER_2015, PAPER_2026, ENDPOINT_GATE, SELECTOR_GATE, CHECKER):
        require(dependency.is_file(), f"missing dependency: {dependency}")
    contract = source_contract()
    certified = event_interval_certificate()
    diagnostic = numerical_diagnostic()
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_source_transition_endpoint_scheduling_gate",
        "status": "exact_coupled_roster_decomposition_and_source_approximation_boundary_certified",
        "passed": True,
        "source_contract": contract,
        "certificate": certified,
        "exact_identities": {
            "ideal_roster_change": "sum(R_plus)-sum(R_minus)=F_t(B_plus)-F_t(C)",
            "source_scheduling_change": "Delta_source=(T_C-F_t(C))+F_t(B_plus)",
            "fixed_chart_pre": "source_pre=fixed_chart-F_t(B_plus)",
            "fixed_chart_post": "source_post=fixed_chart+(T_C-F_t(C))",
        },
        "numerical_diagnostic": diagnostic,
        "paper_audit": {
            "lewis_2015": {
                "pages": [4, 5, 17, 23, 49, 50, 51],
                "finding": "hybrid and transition formulae are approximations; part of the transition zone is assigned to direct numerical integration; some observed transition correlations are not proved",
            },
            "lewis_brereton_2026": {
                "pages": [45, 46],
                "finding": "the proposed Gaussian-sum hybrid is explicitly not exact and part of the transition term is assigned to numerical integration",
            },
        },
        "decision": {
            "earlier_coupled_source_event_decomposed": True,
            "transition_is_exact_handoff_identity": False,
            "moving_source_endpoint_admissible_as_theorem_selector": False,
            "fixed_analytic_endpoint_for_next_chart": B_PLUS,
            "rs_cutoff_changes_at_event": False,
            "complete_source_error_bounded": False,
            "rh_implication": False,
        },
        "dependencies": {
            "upper_endpoint_gate": {"path": relative(ENDPOINT_GATE), "sha256": file_hash(ENDPOINT_GATE)},
            "selector_strip_gate": {"path": relative(SELECTOR_GATE), "sha256": file_hash(SELECTOR_GATE)},
            "pinned_source": {"path": relative(SOURCE), "sha256": file_hash(SOURCE)},
            "paper_2015": {"path": relative(PAPER_2015), "sha256": file_hash(PAPER_2015)},
            "paper_2026": {"path": relative(PAPER_2026), "sha256": file_hash(PAPER_2026)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
        "next_obligation": (
            "Instantiate the quantitative 399-mode selector chart at fixed B=5122423, keeping the exact symmetric one-cell Poisson "
            "splice separate from the source approximation-error ledger."
        ),
        "proof_boundary": (
            "Exact finite-roster scheduling algebra, transition-root interval, and unchanged RS floor only. The 64-node transition "
            "value is diagnostic. No complete source-error bound, quantitative 399-mode theorem, T_upper, Lambda<=0, RH, or prize-level conclusion."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("built source-transition endpoint scheduling gate: exact roster decomposition, fixed B=5122423")


if __name__ == "__main__":
    main()
