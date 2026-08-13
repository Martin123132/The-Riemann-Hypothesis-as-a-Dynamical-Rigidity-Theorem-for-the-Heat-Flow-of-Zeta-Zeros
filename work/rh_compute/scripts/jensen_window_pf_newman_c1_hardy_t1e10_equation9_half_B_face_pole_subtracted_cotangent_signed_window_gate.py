#!/usr/bin/env python3
"""Certify the signed physical projection of the first B cotangent current."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_pole_subtracted_cotangent_signed_window_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
CACHE = REPO_ROOT / f"work/rh_compute/results/cache/{STEM}_rows.jsonl"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "cotangent": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_pole_subtracted_cotangent_current_gate.json",
    "B_window": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_gaussian_window_fresnel_remainder_gate.json",
    "paired_residual": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_paired_target_residual_normal_form_gate.json",
    "Morse_normalization": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_universal_logistic_morse_characteristic_fold_reduction_gate.json",
}

PRECISION = 75
TOLERANCE = "1e-15"
FORMULA_VERSION = "stable-local-cot-v2-precision75"
T = 10_000_000_000
B = 5_122_421
WINDOW_XI = 70
LOCAL_MODES = (621, 622)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
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


def complex_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": value.real.str(PRECISION, more=True),
        "imag_ball": value.imag.str(PRECISION, more=True),
    }


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def symbolic_certificate() -> dict[str, str]:
    x, z, endpoint, mode, t = sp.symbols("x z B m t", positive=True, real=True)
    endpoint_phase = sp.pi * mode**2 / z + sp.pi * endpoint**2 * z / 4
    outer_phase = -sp.pi * mode**2 / x + t * sp.log((1 - x) / x) / 2
    face_phase = sp.simplify((endpoint_phase + outer_phase).subs(z, x))
    expected = sp.pi * endpoint**2 * x / 4 + t * sp.log((1 - x) / x) / 2
    require(sp.simplify(face_phase - expected) == 0, "mode did not cancel on the B face")

    c, k = sp.symbols("c k", positive=True, real=True)
    delta = c - k
    stable = (
        (sp.pi * sp.cot(sp.pi * delta) - 1 / delta) / (8 * c)
        - 1 / (8 * c * (c + k))
        - 1 / (8 * c**2)
    )
    original = sp.pi * sp.cot(sp.pi * c) / (8 * c) - 1 / (8 * c**2) - 1 / (4 * (c**2 - k**2))
    # cot(pi(c-k))=cot(pi c) for integer k.
    integer_shifted = original.xreplace({sp.cot(sp.pi * c): sp.cot(sp.pi * delta)})
    require(sp.simplify(stable - integer_shifted) == 0, "stable one-pole removal failed")

    return {
        "face_phase": "Phi_B(x)=pi*B^2*x/4+(t/2)log((1-x)/x)",
        "mode_cancellation": "phi_(m,B)(x)+psi_m(x)=Phi_B(x) exactly",
        "pole_subtracted_sum": "S_hat(c)=pi*cot(pi*c)/(8c)-1/(8c^2)-1/[4(c^2-621^2)]-1/[4(c^2-622^2)]",
        "first_current": "C_hat(x)=-2(2+i*pi*B^2*x)S_hat(Bx/2)/pi^2",
        "stable_local_form": "After choosing k=621 or 622 and delta=c-k, combine pi*cot(pi*delta)/(8c) with the kth pole before evaluation.",
        "regular_cot_series": "pi*cot(pi*delta)-1/delta=-2*sum_(n>=1)zeta(2n)delta^(2n-1)",
        "series_tail_bound": "For |delta|=r<1, the tail after n=8 is at most 4*r^17/(1-r^2), using zeta(2n)<2.",
        "physical_projection": "2*(pi/(32t))^(1/4) Re[e^(-i*pi/8) integral W_t(x)e^(i*Phi_B(x))C_hat(x)dx]",
    }


class SignedWindow:
    def __init__(self) -> None:
        self.t = arb(T)
        self.endpoint = arb(B)
        self.pi = arb.pi()
        self.i = acb(0, 1)
        self.x0 = (1 - (1 - 8 * self.t / (self.pi * self.endpoint**2)).sqrt()) / 2
        self.hessian = self.t * (1 - 2 * self.x0) / (2 * self.x0**2 * (1 - self.x0) ** 2)
        self.root_hessian = self.hessian.sqrt()
        self.log0 = ((1 - self.x0) / self.x0).log()
        self.phase0 = self.pi * self.endpoint**2 * self.x0 / 4 + self.t * self.log0 / 2
        self.normalizer = (self.pi / (32 * self.t)) ** (arb(1) / 4)

    def regular_cot(self, delta: acb) -> acb:
        radius = abs(delta).upper()
        if radius < arb("0.2"):
            value = acb(0)
            for index in range(1, 9):
                value -= 2 * arb(2 * index).zeta() * delta ** (2 * index - 1)
            remainder = 4 * radius**17 / (1 - radius**2)
            return value + acb(arb(0, remainder), arb(0, remainder))
        return self.pi * (self.pi * delta).cot() - 1 / delta

    def stable_sum(self, c: acb, local_mode: int) -> acb:
        other = 1243 - local_mode
        delta = c - local_mode
        return (
            self.regular_cot(delta) / (8 * c)
            - 1 / (8 * c * (c + local_mode))
            - 1 / (8 * c**2)
            - 1 / (4 * (c**2 - other**2))
        )

    def integrand(self, local_mode: int):
        def function(xi: acb, _: bool) -> acb:
            x = self.x0 + xi / self.root_hessian
            c = self.endpoint * x / 2
            current = -2 * (2 + self.i * self.pi * self.endpoint**2 * x) * self.stable_sum(c, local_mode) / self.pi**2
            phase = (
                self.pi * self.endpoint**2 * (x - self.x0) / 4
                + self.t * (((1 - x) / x).log() - self.log0) / 2
            )
            weight = (x * (1 - x)) ** (-arb(1) / 4)
            return weight * current * (self.i * phase).exp() / self.root_hessian

        return function

    def cuts(self) -> list[arb]:
        cuts = [arb(-WINDOW_XI + index) for index in range(2 * WINDOW_XI + 1)]
        cuts.append(self.root_hessian * (2 * arb("621.5") / self.endpoint - self.x0))
        return sorted(cuts, key=float)


def load_cache(cuts: list[arb]) -> dict[int, dict[str, Any]]:
    rows: dict[int, dict[str, Any]] = {}
    if not CACHE.is_file():
        return rows
    for line in CACHE.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        index = row.get("panel_index")
        if (
            row.get("formula_version") == FORMULA_VERSION
            and isinstance(index, int)
            and 0 <= index < len(cuts) - 1
            and row.get("xi_left") == str(cuts[index])
            and row.get("xi_right") == str(cuts[index + 1])
        ):
            rows[index] = row
    return rows


def append_cache(row: dict[str, Any]) -> None:
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    with CACHE.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(row, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def integrate(window: SignedWindow) -> tuple[acb, list[dict[str, Any]]]:
    cuts = window.cuts()
    rows = load_cache(cuts)
    for index, (left, right) in enumerate(zip(cuts, cuts[1:])):
        if index in rows:
            continue
        center = (left + right) / 2
        c_center = window.endpoint * (window.x0 + center / window.root_hessian) / 2
        local_mode = 621 if c_center < arb("621.5") else 622
        started = time.time()
        value = acb.integral(
            window.integrand(local_mode),
            left,
            right,
            abs_tol=arb(TOLERANCE),
            rel_tol=arb(TOLERANCE),
            eval_limit=250_000,
            depth_limit=40,
        )
        row = {
            "formula_version": FORMULA_VERSION,
            "panel_index": index,
            "xi_left": str(left),
            "xi_right": str(right),
            "local_removal_chart": local_mode,
            "value": complex_record(value),
            "elapsed_seconds": round(time.time() - started, 3),
        }
        append_cache(row)
        rows[index] = row
        print(f"cached signed B cotangent panel {index + 1}/{len(cuts) - 1}", flush=True)
    require(len(rows) == len(cuts) - 1, "incomplete signed-window cache")
    ordered = [rows[index] for index in range(len(cuts) - 1)]
    return sum((parse_complex(row["value"]) for row in ordered), acb(0)), ordered


def numerical_certificate(window: SignedWindow) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    reduced, rows = integrate(window)
    rotation = (window.i * window.phase0 - window.i * window.pi / 8).exp()
    physical = 2 * window.normalizer * (rotation * reduced).real
    normalized_magnitude = 2 * window.normalizer * abs(reduced)
    require(arb("4.791e-6") < physical < arb("4.792e-6"), "signed physical projection escaped enclosure")
    require(arb("2.070e-5") < normalized_magnitude < arb("2.071e-5"), "normalized complex magnitude drift")
    require(abs(physical) < arb("4.792e-6"), "signed current exceeds allocated window bound")
    return {
        "height": T,
        "endpoint": B,
        "window_xi": [-WINDOW_XI, WINDOW_XI],
        "panel_count": len(rows),
        "precision_decimal_digits": PRECISION,
        "quadrature_tolerance_per_panel": TOLERANCE,
        "x_trace_ball": window.x0.str(PRECISION, more=True),
        "trace_hessian_ball": window.hessian.str(PRECISION, more=True),
        "reduced_window_integral_ball": complex_record(reduced),
        "physically_normalized_signed_projection_ball": physical.str(PRECISION, more=True),
        "physically_normalized_signed_projection_absolute_upper_bound": "4.792e-6",
        "physically_normalized_complex_magnitude_ball": normalized_magnitude.str(PRECISION, more=True),
        "reference_physical_target": "8.6e-6",
        "absolute_target_fraction_upper_bound": "0.558",
    }, rows


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["numerical_certificate"]
    return f"""# Signed B cotangent current on the Gaussian face window

Date: 2026-08-13

Status: cancellation-preserving saved-height interval quadrature for one
boundary current; not a proof of the complete B-face estimate

After modes 621 and 622 are removed exactly, Section 11.365 gives the
analytic first face current

```text
C_hat(x)=-2(2+i*pi*B^2*x)S_hat(Bx/2)/pi^2.           (SC1)
```

All retained modes have the same exact trace phase because

```text
phi_(m,B)(x)+psi_m(x)
 =Phi_B(x)=pi*B^2*x/4+(t/2)log((1-x)/x).             (SC2)
```

Let `x_B` be the lower stationary point of `Phi_B`, let
`H_B=Phi_B''(x_B)`, and put `xi=sqrt(H_B)(x-x_B)`.  The quantity certified
here is the signed equation-(9) projection

```text
E_C=2*(pi/(32t))^(1/4) Re{{e^(-i*pi/8)
 integral_(|xi|<=70) W_t(x)e^(i*Phi_B(x))C_hat(x)dx}}. (SC3)
```

The apparent poles at `c=621,622` are never evaluated by subtracting large
singular balls.  On the chart centered at an integer `k`, the code combines
the kth pole with cotangent algebraically and evaluates the regular function

```text
pi*cot(pi*delta)-1/delta
 =-2 sum_(n>=1) zeta(2n)delta^(2n-1), delta=c-k.      (SC4)
```

For `|delta|<0.2`, eight terms are used and the omitted series is enclosed
by `4|delta|^17/(1-|delta|^2)`.  Elsewhere the exact cotangent form is used.
The 141 one-xi panels are append-only and fsynced, so the calculation is
resumable.

Arb outward rounding yields

```text
reduced complex integral
 ={c['reduced_window_integral_ball']},

E_C={c['physically_normalized_signed_projection_ball']},
|E_C|<4.792e-6.                                       (SC5)
```

This is below `0.558` of the reference physical target `8.6e-6`.  The
corresponding normalized complex magnitude is

```text
{c['physically_normalized_complex_magnitude_ball']},  (SC6)
```

which is greater than `2.07e-5`.  Thus (SC5) genuinely uses the prescribed
real projection and trace oscillation.  Taking the complex modulus or a
pointwise absolute value would not close.

This gate does not yet add the exact local 621/622 replacements, the next
two face-current terms, the nonlocal Fresnel remainder, or the outside-window
tails.  Their signed sum must be enclosed before (SC5) can be charged to a
complete B-face budget.

Pi provenance: all occurrences come from the equation-(9) face phase,
Fourier cotangent identity, and paper normalization.  The pi/8 rotation is
the exact odd-square half-Kummer reflection.  No fitted constant is used.

Proof boundary: the first pole-subtracted B face current on `|xi|<=70` at
the saved height only.  It proves no complete B estimate, no A-fold splice,
no complete paired residual or `T_upper`, no height-uniform theorem, no
`Lambda<=0`, no PF-infinity statement, no RH, and no prize-level conclusion.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    ctx.dps = PRECISION
    ctx.threads = 1
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["cotangent"]["decision"]["local_621_622_poles_removed_exactly"] is True, "pole-removal dependency drift")
    require(dependencies["B_window"]["decision"]["xi_70_window_contains_exactly_B_crossings_621_622"] is True, "window dependency drift")
    require(dependencies["paired_residual"]["decision"]["finite_paired_half_domain_target_residual_identity_proved"] is True, "paired residual dependency drift")
    require(dependencies["Morse_normalization"]["decision"]["classical_inverse_sqrt_mode_carrier_recovered"] is True, "normalization dependency drift")

    window = SignedWindow()
    certificate, rows = numerical_certificate(window)
    artifact = {
        "kind": STEM,
        "status": "first_pole_subtracted_B_face_current_signed_window_projection_below_4p792e_minus_6",
        "passed": True,
        "scope": {"height": T, "endpoint": B, "window_xi": [-WINDOW_XI, WINDOW_XI]},
        "symbolic_certificate": symbolic_certificate(),
        "numerical_certificate": certificate,
        "decision": {
            "first_pole_subtracted_current_integrated_before_absolute_value": True,
            "removable_local_poles_evaluated_stably": True,
            "signed_physical_window_projection_below_4p792e_minus_6": True,
            "complex_magnitude_bound_closes_reference_target": False,
            "exact_local_621_622_replacements_added": False,
            "higher_face_currents_added": False,
            "outside_window_tail_bounded": False,
            "complete_B_face_estimate_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "row_cache": {"path": relative(CACHE), "sha256": file_hash(CACHE)},
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "flint_threads": 1,
            "process_priority": priority,
            "cached_rows": len(rows),
        },
        "next_obligation": "In the same signed projection, add the exact 621/622 replacement on their outer-safe arcs, the first three asymptotic local/nonlocal currents on normal-safe arcs, and the certified Fresnel remainders. Then bound the two |xi|>70 B trace tails without separating the paired endpoint completion.",
        "proof_boundary": "Only the first pole-subtracted B face current on the saved-height |xi|<=70 window. No local replacement, higher-current, outside-window, complete B, A-fold, complete paired residual, complete T_upper, height-uniform, Lambda<=0, PF-infinity, RH, or prize-level theorem is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified signed B cotangent window projection below 4.792e-6", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
