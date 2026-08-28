#!/usr/bin/env python3
"""Derive the algebraic H(t) prefactor and the defect target behind Delta_KU."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_exact_H_continuation_defect_target_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
PAPER_2015 = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_2015_1502.06903.pdf"
PAPER_2026 = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_Brereton_2607.15310.pdf"

DEPENDENCIES = {
    "bridge_target": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_QK_exact_hardy_truncation_bridge_target_gate.json",
    "normalization_repair": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_domain_Kummer_normalization_repair_gate.json",
    "component_split": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_lower_rs_hybrid_upper_component_split_gate.json",
}

HEIGHT = 10_000_000_000
PRECISIONS = (100, 140)
DISPLAY_DEFECT_LOWER = "0.0072743687"
DISPLAY_DEFECT_UPPER = "0.0663407986"


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


def central_row(artifact: dict[str, Any]) -> dict[str, Any]:
    rows = [row for row in artifact["output_rows"] if row["output_label"] == "t+0.00"]
    require(len(rows) == 1, "central component row drift")
    return rows[0]


def exact_log_H_direct(t_text: str, precision: int) -> arb:
    ctx.dps = precision
    ctx.threads = 1
    t = arb(t_text)
    pi = arb.pi()
    log_two = arb(2).log()
    log_thermal = (1 + (-2 * pi * t).exp()).log()
    return (
        arb("1.25") * log_two
        + arb("0.25") * t.log()
        + arb("1.5") * pi.log()
        - arb("0.75") * pi * t
        - log_thermal
        - acb(arb("0.25"), t / 2).lgamma().real
        - 2 * acb(arb("0.75"), t / 2).lgamma().real
    )


def exact_log_H_duplication(t_text: str, precision: int) -> arb:
    ctx.dps = precision
    ctx.threads = 1
    t = arb(t_text)
    pi = arb.pi()
    log_two = arb(2).log()
    log_thermal = (1 + (-2 * pi * t).exp()).log()
    # log(cosh(pi*t)) = pi*t-log(2)+log(1+exp(-2*pi*t)), evaluated stably.
    log_cosh = pi * t - log_two + log_thermal
    return (
        arb("0.75") * log_two
        + arb("0.25") * t.log()
        + arb("0.5") * pi.log()
        - arb("0.75") * pi * t
        + arb("0.5") * log_cosh
        - log_thermal
        - acb(arb("0.75"), t / 2).lgamma().real
    )


def render_note(artifact: dict[str, Any]) -> str:
    exact_h = artifact["exact_H_certificate"]
    target = artifact["continuation_defect_target"]
    return f"""# Exact H(t) and continuation-defect target

Date: 2026-08-27

Status: exact normalization and target certified; continuation theorem open

Starting from the Kummer prefactor stated in equation (4) and Euler's
integral, the normalization multiplying the correctly scaled finite Kummer
main simplifies exactly to

```text
H(t)=2^(5/4)t^(1/4)pi^(3/2)exp(-3pi t/4)
     /[(1+exp(-2pi t)) |Gamma(1/4+it/2)|
       |Gamma(3/4+it/2)|^2].                          (HC1)
```

The duplication formula gives the independent equivalent expression

```text
H(t)=2^(3/4)t^(1/4)sqrt(pi)exp(-3pi t/4)sqrt(cosh(pi t))
     /[(1+exp(-2pi t)) |Gamma(3/4+it/2)|].            (HC2)
```

Arb evaluates both logarithmic forms at two precisions.  At `t=10^10`,

```text
log H={exact_h['high_precision_log_H_direct_ball']},
H-1={exact_h['high_precision_H_minus_one_ball']}.       (HC3)
```

The exactness in (HC1)--(HC3) is algebraic exactness of the displayed
special-function prefactor.  It does not certify that the proposed
equation-(4) infinite Kummer series equals Hardy's `Z(t)`; that contour and
interchange theorem remains a separate obligation.

In particular, the entire exact `H-1` correction to the complementary upper
component is below `{exact_h['H_minus_one_times_U_absolute_upper_ball']}`.
It is about nineteen orders of magnitude too small to supply the bridge
target.

Define the exact finite-roster continuation complement and its lower-adjusted
defect by

```text
C_K=Z-H Q_K,
D_K=C_K-L=U-H Q_K.                                    (HC4)
```

Since `Q_K=U+Delta_KU`, this gives the exact reversible identities

```text
D_K=(1-H)U-H Delta_KU,
Delta_KU=((1-H)U-D_K)/H.                              (HC5)
```

Transporting the certified `Delta_KU` target through (HC5) yields the simple
conservative target

```text
{target['displayed_D_K_sufficient_corridor']}.         (HC6)
```

The complete interval endpoints and the corresponding target for `C_K` are
stored in the machine artifact.  Thus a successful continuation theorem must
show that the exact complement `Z-HQ_K` differs from the independently
certified lower Hardy component by a **positive** amount between about
`0.00727` and `0.06634`.

The paper's Appendix-A Euler-Maclaurin route cannot currently certify (HC6):

```text
A64: asymptotic Kummer term with 1+O(Lambda_alpha^-1);
A65: transition approximation with O(g^2) and a stated restricted range;
A71: derivative formula explicitly only to leading order;
A73: Z(t) is written with an approximation sign;
A74: the tail integral is a Laplace estimate with 1+O(1/t).              (HC7)
```

The surrounding text also labels part of the analytic-continuation
cancellation as an unproven assertion and speculative.  Consequently A73 is
a route hint, not an exact identity or a constant-bearing remainder theorem.
The next derivation must begin from the exact equation-(4)/Riemann-Siegel
contour and produce `C_K` with every omitted lower label, upper continuation,
and transition contribution explicitly owned.

Proof boundary: exact `H(t)` normalization, exact defect identities, and a
saved-height sufficient defect target only.  No enclosure of `D_K`, `C_K`,
`Delta_KU`, `Q_K`, or `J_Z`, no non-A bound, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "dependency failure")
    require(
        dependencies["normalization_repair"]["decision"]["full_Kummer_label_prefactor_corrected_from_2c_to_c"]
        is True,
        "normalization repair drift",
    )

    low_log_direct = exact_log_H_direct(str(HEIGHT), PRECISIONS[0])
    low_log_duplication = exact_log_H_duplication(str(HEIGHT), PRECISIONS[0])
    high_log_direct = exact_log_H_direct(str(HEIGHT), PRECISIONS[1])
    high_log_duplication = exact_log_H_duplication(str(HEIGHT), PRECISIONS[1])
    require(low_log_direct.overlaps(high_log_direct), "direct H precision ladder missed")
    require(low_log_duplication.overlaps(high_log_duplication), "duplication H precision ladder missed")
    require(high_log_direct.overlaps(high_log_duplication), "exact H formulas do not overlap")
    high_H = high_log_direct.exp()
    require((high_H - 1).lower() > 0, "H-1 lost positive sign")

    component = central_row(dependencies["component_split"])
    lower = arb(component["high_precision"]["exact_lower_component_ball"])
    exact_upper = arb(component["high_precision"]["exact_complementary_upper_ball"])
    bridge = dependencies["bridge_target"]["bridge_target_certificate"]
    require(exact_upper.overlaps(arb(bridge["exact_complementary_upper_ball"])), "exact upper dependency drift")
    delta_lower = arb(bridge["published_Delta_KU_safe_lower_ball"])
    delta_upper = arb(bridge["published_Delta_KU_safe_upper_ball"])

    h_correction = (high_H - 1) * exact_upper
    require(abs(h_correction).upper() < arb("7e-22"), "H correction scale guard failed")
    one_minus_H_U = (1 - high_H) * exact_upper
    d_lower_boundary = one_minus_H_U - high_H * delta_upper
    d_upper_boundary = one_minus_H_U - high_H * delta_lower
    d_safe_lower = d_lower_boundary.upper()
    d_safe_upper = d_upper_boundary.lower()
    require(d_safe_lower < d_safe_upper, "continuation-defect corridor is empty")
    display_lower = arb(DISPLAY_DEFECT_LOWER)
    display_upper = arb(DISPLAY_DEFECT_UPPER)
    require(display_lower > d_safe_lower, "displayed D_K lower endpoint is not conservative")
    require(display_upper < d_safe_upper, "displayed D_K upper endpoint is not conservative")
    require(display_lower < display_upper, "displayed D_K corridor is empty")

    c_safe_lower = d_safe_lower + lower.upper()
    c_safe_upper = d_safe_upper + lower.lower()
    require(c_safe_lower < c_safe_upper, "C_K corridor is empty")

    artifact = {
        "kind": STEM,
        "status": "exact_H_normalization_and_lower_adjusted_continuation_defect_target_certified_A73_nonpromotion_guarded",
        "passed": True,
        "scope": {
            "height": HEIGHT,
            "precisions_decimal_digits": list(PRECISIONS),
            "source_label_count": 2_481_423,
        },
        "exact_H_certificate": {
            "direct_formula": "H=2^(5/4)t^(1/4)pi^(3/2)e^(-3pi*t/4)/[(1+e^(-2pi*t))|Gamma(1/4+it/2)||Gamma(3/4+it/2)|^2]",
            "duplication_formula": "H=2^(3/4)t^(1/4)sqrt(pi)e^(-3pi*t/4)sqrt(cosh(pi*t))/[(1+e^(-2pi*t))|Gamma(3/4+it/2)|]",
            "low_precision_log_H_direct_ball": low_log_direct.str(90, more=True),
            "low_precision_log_H_duplication_ball": low_log_duplication.str(90, more=True),
            "high_precision_log_H_direct_ball": high_log_direct.str(120, more=True),
            "high_precision_log_H_duplication_ball": high_log_duplication.str(120, more=True),
            "high_precision_H_ball": high_H.str(120, more=True),
            "high_precision_H_minus_one_ball": (high_H - 1).str(120, more=True),
            "scaled_32_t2_H_minus_one_ball": ((high_H - 1) * 32 * arb(HEIGHT) ** 2).str(100, more=True),
            "H_minus_one_times_U_ball": h_correction.str(110, more=True),
            "H_minus_one_times_U_absolute_upper_ball": abs(h_correction).upper().str(100, more=True),
            "formula_overlap": True,
            "precision_ladder_overlap": True,
        },
        "continuation_defect_identity": {
            "finite_continuation_complement": "C_K=Z-H*Q_K",
            "lower_adjusted_defect": "D_K=C_K-L=U-H*Q_K",
            "forward": "D_K=(1-H)U-H*Delta_KU",
            "inverse": "Delta_KU=((1-H)U-D_K)/H",
        },
        "continuation_defect_target": {
            "exact_lower_Hardy_component_ball": lower.str(110, more=True),
            "exact_complementary_upper_ball": exact_upper.str(110, more=True),
            "Delta_KU_safe_lower_ball": delta_lower.str(100, more=True),
            "Delta_KU_safe_upper_ball": delta_upper.str(100, more=True),
            "D_K_safe_lower_ball": d_safe_lower.str(100, more=True),
            "D_K_safe_upper_ball": d_safe_upper.str(100, more=True),
            "displayed_D_K_sufficient_corridor": f"{DISPLAY_DEFECT_LOWER}<D_K<{DISPLAY_DEFECT_UPPER}",
            "C_K_safe_lower_ball": c_safe_lower.str(100, more=True),
            "C_K_safe_upper_ball": c_safe_upper.str(100, more=True),
            "interpretation": "The exact finite-roster continuation complement C_K=Z-HQ_K must exceed the exact lower Hardy component L by the positive defect D_K.",
        },
        "paper_exactness_audit": {
            "A64": "asymptotic approximation with multiplicative 1+O(Lambda_alpha^-1)",
            "A65": "transition approximation includes O(g^2) and the text states a restricted useful range",
            "A71": "odd derivative estimate stated to leading order",
            "A73": "Euler-Maclaurin reconstruction of Z(t) retains an approximation sign and uses A64 terms",
            "A74": "tail integral is estimated by Laplace's method with 1+O(1/t)",
            "analytic_continuation_text": "The appendix calls one cancellation assertion unproven and discusses the continuation speculatively.",
            "nonpromotion": "Neither A73 nor the 2026 equation-(9) REM label supplies the exact constant-bearing C_K identity required here.",
        },
        "decision": {
            "exact_H_derived_from_equation4_prefactor": True,
            "exact_H_independent_duplication_check_passed": True,
            "equation4_infinite_Kummer_series_equals_Z_certified_here": False,
            "H_normalization_can_explain_required_bridge_scale": False,
            "single_lower_adjusted_continuation_defect_D_K_isolated": True,
            "rigorous_sufficient_D_K_corridor_nonempty": True,
            "A73_permitted_as_exact_continuation_bridge": False,
            "exact_equation4_or_Riemann_Siegel_contour_required": True,
            "D_K_enclosed": False,
            "Delta_KU_enclosed": False,
            "non_A_bound_proved": False,
            "rh_implication": False,
        },
        "references": {
            "Lewis_2015": {"path": relative(PAPER_2015), "sha256": file_hash(PAPER_2015), "relevant_equations": ["A64", "A65", "A71", "A73", "A74"]},
            "Lewis_Brereton_2026": {"path": relative(PAPER_2026), "sha256": file_hash(PAPER_2026), "relevant_equations": [4, 7, 8, 9]},
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
        "next_obligation": "Derive C_K=Z-HQ_K from the exact equation-(4) or Riemann-Siegel contour before asymptotic expansion. Partition it into an exact lower-label contour, exact upper analytic-continuation contour, and exact transition term, prove the displayed positive D_K corridor with all signs joined, and reject any import of A64/A73 big-O terms without explicit uniform constants.",
        "proof_boundary": "Exact algebraic H-prefactor normalization, exact continuation-defect definitions, paper exactness audit, and a saved-height sufficient D_K target only. The equation-(4) infinite Kummer series is not certified here as an exact identity for Z. No D_K, C_K, Delta_KU, Q_K, or J_Z enclosure, non-A bound, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified exact H and isolated the lower-adjusted continuation-defect target", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
