#!/usr/bin/env python3
"""Quarantine the heuristic A21 continuation and expose the exact contour target."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_A21_exactness_provenance_and_contour_jump_target_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
PAPER_2015 = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_2015_1502.06903.pdf"
PAPER_2026 = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_Brereton_2607.15310.pdf"

DEPENDENCIES = {
    "exact_H_defect_target": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_exact_H_continuation_defect_target_gate.json",
    "normalization_repair": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_domain_Kummer_normalization_repair_gate.json",
}

PRECISIONS = (90, 130)
SAMPLE_Q = ("-2", "-0.5", "-0.125", "0.125", "0.5", "2")
PARTIAL_TERMS = 24
BOUNDARY_LAYER_N = (8, 16, 32, 64)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {relative(path)}")
    return json.loads(path.read_text(encoding="utf-8"))


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def contour_values(q_value: arb, sign: int, terms: int) -> dict[str, acb | arb]:
    """Return the exact secant, a valid/invalid geometric partial sum, and residual.

    sign=+1 is the A13 expansion and sign=-1 is the A18 expansion.
    """

    pi = arb.pi()
    root_two = arb(2).sqrt()
    w = acb(pi * q_value / root_two, -pi * q_value / root_two)
    secant = 1 / w.cos()
    partial = acb(0)
    for k in range(terms):
        alpha = 2 * k + 1
        partial += 2 * (-1 if k % 2 else 1) * (acb(0, sign * alpha) * w).exp()
    ratio = -(acb(0, 2 * sign) * w).exp()
    residual = secant - partial
    residual_formula = secant * ratio**terms
    return {
        "w": w,
        "secant": secant,
        "partial": partial,
        "ratio": ratio,
        "residual": residual,
        "residual_formula": residual_formula,
    }


def row(q_text: str, precision: int) -> dict[str, Any]:
    ctx.dps = precision
    ctx.threads = 1
    q_value = arb(q_text)
    sign = 1 if q_value.upper() < 0 else -1
    values = contour_values(q_value, sign, PARTIAL_TERMS)
    require(values["residual"].overlaps(values["residual_formula"]), "finite residual identity failed")
    ratio_abs = abs(values["ratio"])
    require(ratio_abs.upper() < 1, "selected geometric expansion is not contractive")
    return {
        "q": q_text,
        "valid_expansion": "A13" if sign == 1 else "A18",
        "terms": PARTIAL_TERMS,
        "secant_real_ball": values["secant"].real.str(70, more=True),
        "secant_imag_ball": values["secant"].imag.str(70, more=True),
        "residual_real_ball": values["residual"].real.str(70, more=True),
        "residual_imag_ball": values["residual"].imag.str(70, more=True),
        "closed_residual_real_ball": values["residual_formula"].real.str(70, more=True),
        "closed_residual_imag_ball": values["residual_formula"].imag.str(70, more=True),
        "ratio_absolute_ball": ratio_abs.str(70, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    target = artifact["contour_defect_target"]
    boundary = artifact["nonuniform_boundary_layer_certificate"]
    return f"""# Equation-(4) A21 provenance and exact contour-defect target

Date: 2026-08-27

Status: provenance gap isolated; exact split-contour theorem open

The finite odd-label transform `Q_K` remains an exact finite mathematical
object, and the repaired Euler-Kummer identity still maps the complete finite
source roster to it.  What is not currently licensed is the further statement
that the infinite Kummer construction is an exact representation of Hardy's
`Z(t)`.

The 2026 paper calls equation (4) an exact formulation and sends the reader to
the 2015 Appendix A for its derivation.  That appendix starts from the exact
Riemann-Siegel contour

```text
z=1/2+q exp(-i*pi/4),     -infinity<q<infinity,       (CP1)
```

but uses different geometric expansions on the two open half-contours.  Put

```text
w=pi q exp(-i*pi/4),     1/sin(pi z)=1/cos(w).
```

Then the convergent identities are

```text
q<0:  1/cos(w)=2 sum_(k>=0)(-1)^k exp(i(2k+1)w),    (CP2)
q>0:  1/cos(w)=2 sum_(k>=0)(-1)^k exp(-i(2k+1)w).   (CP3)
```

These are the source paper's A13 and A18 expansions.  A13 is not convergent
on `q>0`, A18 is not convergent on `q<0`, and neither ordinary series converges
at `q=0`.

For `N` terms the exact residuals are

```text
q<0: R_N^-(q)=sec(w)[-exp(2iw)]^N,
q>0: R_N^+(q)=sec(w)[-exp(-2iw)]^N.                  (CP4)
```

The Arb rows in the machine artifact verify (CP4) directly at six points and
two precisions.  More importantly, convergence is not uniform at the contour
crossing.  On `q=-u/N`,

```text
|[-exp(2iw)]^N|=exp(-sqrt(2) pi u),                 (CP5)
```

which is independent of `N`.  For `u=1` the certified common value is
`{boundary['expected_boundary_factor_ball']}`.  Pointwise decay away from
`q=0` therefore does not by itself justify exchanging the whole-contour
integral, infinite sum, and limit.

Appendix A acknowledges this issue.  It says A13 is only convergent on its own
half-contour, reports that matching the two termwise asymptotic constructions
fails, introduces A21 as a heuristic continuation of A13 across the entire
`q` range, and sets the convergence question aside.  A33 is subsequently
conditional on suitable convergence criteria.  Later numerical and
asymptotic agreement is evidence for the proposed construction, but is not a
proof of the missing interchange or boundary term.

The rigorous replacement is now precise.  Starting from the exact contour,
split at `q=0`, retain A13 only on `q<0` and A18 only on `q>0`, introduce finite
cutoffs before interchanging sum and integral, and analyze the `q=O(1/N)`
boundary layer.  The chart change at `q=-1/sqrt(2)`, where `Re(z)=0`, must also
retain its branch convention.  After matching the finite roster, the resulting
projected complement must be proved to equal

```text
C_K=Z-H Q_K,       D_K=C_K-L=U-H Q_K.               (CP6)
```

At the saved height the sufficient scalar target is

```text
{target['displayed_D_K_sufficient_corridor']}.       (CP7)
```

The source of every `pi` in (CP1)--(CP5) is explicit: it comes from the
Riemann-Siegel denominator `sin(pi z)` and the fixed contour rotation
`exp(-i*pi/4)`, not from an inserted geometric circle constant.

Proof boundary: this gate certifies the two valid geometric half-contour
identities, their exact finite residuals, the nonuniform crossing layer, the
source-provenance nonpromotion, and the inherited sufficient `D_K` target.
It does not derive the contour defect, prove that a single boundary layer is
its only contribution, enclose `C_K`, `D_K`, `Delta_KU`, `Q_K`, or `J_Z`, or
prove a non-A bound, an all-height theorem, `Lambda<=0`, PF-infinity, RH, or a
prize-level conclusion.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "dependency failure")

    low_rows = [row(q_text, PRECISIONS[0]) for q_text in SAMPLE_Q]
    high_rows = [row(q_text, PRECISIONS[1]) for q_text in SAMPLE_Q]
    for low, high in zip(low_rows, high_rows, strict=True):
        require(low["q"] == high["q"], "precision-row order drift")
        require(arb(low["secant_real_ball"]).overlaps(arb(high["secant_real_ball"])), "secant real precision miss")
        require(arb(low["secant_imag_ball"]).overlaps(arb(high["secant_imag_ball"])), "secant imag precision miss")
        require(arb(low["residual_real_ball"]).overlaps(arb(high["residual_real_ball"])), "residual real precision miss")
        require(arb(low["residual_imag_ball"]).overlaps(arb(high["residual_imag_ball"])), "residual imag precision miss")

    ctx.dps = PRECISIONS[1]
    ctx.threads = 1
    pi = arb.pi()
    root_two = arb(2).sqrt()
    expected_boundary = (-root_two * pi).exp()
    boundary_rows = []
    for terms in BOUNDARY_LAYER_N:
        q_value = -arb(1) / terms
        values = contour_values(q_value, 1, terms)
        factor = abs(values["ratio"]) ** terms
        require(factor.overlaps(expected_boundary), "q=-1/N boundary factor drift")
        boundary_rows.append(
            {
                "terms": terms,
                "q": q_value.str(50, more=True),
                "omitted_geometric_factor_absolute_ball": factor.str(80, more=True),
            }
        )

    invalid_a13 = contour_values(arb("0.25"), 1, 3)
    invalid_a18 = contour_values(arb("-0.25"), -1, 3)
    require(abs(invalid_a13["ratio"]).lower() > 1, "A13 invalid-side growth guard failed")
    require(abs(invalid_a18["ratio"]).lower() > 1, "A18 invalid-side growth guard failed")

    exact_h = dependencies["exact_H_defect_target"]
    defect_target = exact_h["continuation_defect_target"]
    normalization = dependencies["normalization_repair"]
    require(
        normalization["decision"]["joined_J_Z_identity_retained_in_corrected_normalization"] is True,
        "corrected finite source transform was not retained",
    )

    artifact = {
        "kind": STEM,
        "status": "A21_heuristic_continuation_quarantined_exact_split_contour_boundary_layer_and_DK_target_certified",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "precisions_decimal_digits": list(PRECISIONS),
            "sample_q": list(SAMPLE_Q),
            "partial_terms": PARTIAL_TERMS,
        },
        "source_provenance_audit": {
            "Lewis_2015_pdf_pages_1_based": {
                "37": "A13 is stated to converge only on q<0; A18 replaces it on q>0; direct matching of the termwise asymptotic constructions fails.",
                "38": "A21 extends the A13 term integrals across the whole q range heuristically; correctness and convergence are explicitly set aside.",
                "41": "A33 is asserted only provided suitable convergence criteria hold.",
                "42": "A33 is described as proposed and supported by later asymptotic/numerical analysis.",
                "52": "The later Euler-Maclaurin discussion retains an unproved all-orders cancellation assertion.",
                "85": "The appendix describes its numerical agreement as indicative rather than a formal proof.",
            },
            "Lewis_Brereton_2026_pdf_pages_1_based": {
                "3": "Equation (4) is presented and its derivation is delegated to the 2015 source.",
                "5": "Equation (4) is called an exact formulation before equation (7) is obtained through Euler's integral.",
            },
            "logical_conclusion": "The 2026 exactness label does not replace the unresolved A21 whole-contour continuation/interchange proof in its cited derivation.",
        },
        "half_contour_identities": {
            "contour": "z=1/2+q*exp(-i*pi/4)",
            "w": "w=pi*q*exp(-i*pi/4)",
            "denominator": "1/sin(pi*z)=1/cos(w)",
            "A13_valid_q_negative": "2*sum_(k>=0)(-1)^k*exp(i*(2k+1)*w)",
            "A18_valid_q_positive": "2*sum_(k>=0)(-1)^k*exp(-i*(2k+1)*w)",
            "A13_finite_residual": "sec(w)*[-exp(2*i*w)]^N",
            "A18_finite_residual": "sec(w)*[-exp(-2*i*w)]^N",
            "q_zero_ordinary_series_convergent": False,
            "high_precision_rows": high_rows,
            "precision_ladder_overlap": True,
        },
        "invalid_side_guard": {
            "A13_at_q_positive_ratio_absolute_ball": abs(invalid_a13["ratio"]).str(80, more=True),
            "A18_at_q_negative_ratio_absolute_ball": abs(invalid_a18["ratio"]).str(80, more=True),
            "term_test_fails_on_each_invalid_half_contour": True,
            "q_zero_term_absolute_value": "2",
        },
        "nonuniform_boundary_layer_certificate": {
            "scaling": "q=-u/N on A13 (and q=+u/N on A18)",
            "exact_factor": "|ratio|^N=exp(-sqrt(2)*pi*u)",
            "u": 1,
            "expected_boundary_factor_ball": expected_boundary.str(100, more=True),
            "rows": boundary_rows,
            "uniform_convergence_through_q_zero": False,
        },
        "exact_contour_obligation": {
            "step_1": "Start from the exact Riemann-Siegel contour and fix all logarithm and power branches.",
            "step_2": "Split q<0 and q>0 before expanding the secant denominator; use A13 and A18 only in their valid domains.",
            "step_3": "Introduce finite contour, series, and puncture cutoffs; justify every interchange away from q=0 by the exact finite residuals.",
            "step_4": "Rescale q=u/N at the crossing and derive the surviving boundary-layer or jump functional without assuming it vanishes.",
            "step_5": "Track the Re(z)=0 chart point q=-1/sqrt(2) and retain the chosen branch through every contour deformation.",
            "step_6": "Match the corrected finite Kummer roster H*Q_K and derive the complete complement C_K=Z-H*Q_K.",
            "step_7": "Prove the inherited saved-height D_K corridor before any RH-level promotion.",
        },
        "contour_defect_target": {
            "C_K_definition": "C_K=Z-H*Q_K",
            "D_K_definition": "D_K=C_K-L=U-H*Q_K",
            "D_K_safe_lower_ball": defect_target["D_K_safe_lower_ball"],
            "D_K_safe_upper_ball": defect_target["D_K_safe_upper_ball"],
            "displayed_D_K_sufficient_corridor": defect_target["displayed_D_K_sufficient_corridor"],
            "boundary_layer_alone_asserted_to_equal_D_K": False,
        },
        "decision": {
            "finite_QK_definition_exact": True,
            "finite_source_to_QK_transform_identity_retained": True,
            "equation4_exactness_imported_from_2015": False,
            "A21_permitted_as_exact_whole_contour_interchange": False,
            "A33_permitted_as_unconditional_identity": False,
            "later_asymptotic_or_numerical_agreement_closes_A21": False,
            "exact_split_contour_route_selected": True,
            "q_zero_boundary_layer_requires_analysis": True,
            "C_K_enclosed": False,
            "D_K_enclosed": False,
            "non_A_bound_proved": False,
            "rh_implication": False,
        },
        "references": {
            "Lewis_2015": {
                "path": relative(PAPER_2015),
                "sha256": file_hash(PAPER_2015),
                "relevant_equations": ["A1", "A12", "A13", "A17", "A18", "A20", "A21", "A33", "A40", "A64", "A73"],
            },
            "Lewis_Brereton_2026": {
                "path": relative(PAPER_2026),
                "sha256": file_hash(PAPER_2026),
                "relevant_equations": [4, 7, 8, 9],
            },
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "arb_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Derive the exact split-contour finite-N identity with puncture and contour cutoffs, rescale the q=O(1/N) crossing layer, and prove whether the complete boundary/branch contribution yields C_K=Z-HQ_K with 0.0072743687<D_K<0.0663407986.",
        "proof_boundary": "Exact half-contour secant identities, finite residual formulas, a nonuniform q=O(1/N) crossing certificate, source-provenance nonpromotion, and the inherited saved-height D_K target only. No derivation or enclosure of the contour defect, C_K, D_K, Delta_KU, Q_K, J_Z, non-A bound, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified A21 nonpromotion and exact split-contour boundary-layer target", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
