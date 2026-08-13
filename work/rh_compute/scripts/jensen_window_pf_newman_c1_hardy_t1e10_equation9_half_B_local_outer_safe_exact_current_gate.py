#!/usr/bin/env python3
"""Certify the cancellation-safe local B current for modes 621 and 622."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_local_outer_safe_exact_current_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
CACHE = REPO_ROOT / f"work/rh_compute/results/cache/{STEM}_rows.jsonl"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "step_tail": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_positive_B_crossing_step_tail_normal_form_gate.json",
    "B_window": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_gaussian_window_fresnel_remainder_gate.json",
    "local_normal_remainder": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_local_normal_safe_fresnel_remainder_gate.json",
    "signed_background": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_pole_subtracted_cotangent_signed_window_gate.json",
    "pair_triangle": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_separable_triangle_geometry_gate.json",
    "truncation_dictionary": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_fresnel_boundary_truncation_dictionary_gate.json",
}

PRECISION = 60
TOLERANCE = "2e-10"
PANEL_WIDTH_XI = "0.5"
FORMULA_VERSION = "local-B-exact-positive-direct-fresnel-v1-dps60-half-xi"
T = 10_000_000_000
B = 5_122_421
WINDOW_XI = 70
MODES = (621, 622)
SPLIT_X = {
    621: (24_238_523_659, 100_000_000_000_000),
    622: (48_583_287_419, 200_000_000_000_000),
}


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
    x, endpoint, mode = sp.symbols("x B m", positive=True, real=True)
    q = sp.sqrt(x / 2) * (endpoint - 2 * mode / x)
    completion = sp.simplify(-mode**2 / x - endpoint**2 * x / 4 + q**2 / 2)
    require(sp.simplify(completion + endpoint * mode) == 0, "face carrier completion failed")

    c = endpoint * x / 2
    positive_terms = (
        endpoint / (sp.I * sp.pi * (endpoint * x - 2 * mode))
        - mode / (2 * sp.pi**2 * (c - mode) ** 3)
        + 3 * sp.I * mode * x / (4 * sp.pi**3 * (c - mode) ** 5)
    )
    negative_terms = (
        endpoint / (sp.I * sp.pi * (endpoint * x + 2 * mode))
        + mode / (2 * sp.pi**2 * (c + mode) ** 3)
        - 3 * sp.I * mode * x / (4 * sp.pi**3 * (c + mode) ** 5)
    )
    require(sp.simplify(positive_terms + negative_terms).has(mode), "paired direct truncation lost mode dependence")
    return {
        "positive_residual": "R_B,m^+=P_B,m+1_(m<=621)P_bulk,m=U_B,m+[1_(m<=621)-1_(q_B<0)]P_bulk,m",
        "exact_positive_face_coefficient": "C_exact,+,m=exp(-i*pi*q_B^2/2)R_B,m^+/x, using B odd",
        "direct_positive_three_terms": "B/[i*pi(Bx-2m)]-m/[2*pi^2(c-m)^3]+3i*m*x/[4*pi^3(c-m)^5]",
        "direct_negative_three_terms": "B/[i*pi(Bx+2m)]+m/[2*pi^2(c+m)^3]-3i*m*x/[4*pi^3(c+m)^5]",
        "outer_safe_splice": "Use C_exact,+ plus the negative three-term tail where the positive q derivative is small; use both direct three-term tails on the normal-safe arc.",
        "crossing_guard": "The exact U_B jump and target-step jump cancel before integration; panels are cut at q_B=0.",
        "physical_projection": "2*(pi/(32t))^(1/4) Re[e^(-i*pi/8) integral W_t e^(i*Phi_B) C_local dx]",
    }


class LocalWindow:
    def __init__(self) -> None:
        self.t = arb(T)
        self.endpoint = arb(B)
        self.pi = arb.pi()
        self.i = acb(0, 1)
        self.full_fresnel = acb(1, 1) / 2
        self.x0 = (1 - (1 - 8 * self.t / (self.pi * self.endpoint**2)).sqrt()) / 2
        self.hessian = self.t * (1 - 2 * self.x0) / (2 * self.x0**2 * (1 - self.x0) ** 2)
        self.root_hessian = self.hessian.sqrt()
        self.log0 = ((1 - self.x0) / self.x0).log()
        self.phase0 = self.pi * self.endpoint**2 * self.x0 / 4 + self.t * self.log0 / 2
        self.normalizer = (self.pi / (32 * self.t)) ** (arb(1) / 4)

    def crossing_xi(self, mode: int) -> arb:
        return self.root_hessian * (2 * arb(mode) / self.endpoint - self.x0)

    def split_x(self, mode: int) -> arb:
        numerator, denominator = SPLIT_X[mode]
        return arb(numerator) / denominator

    def split_xi(self, mode: int) -> arb:
        return self.root_hessian * (self.split_x(mode) - self.x0)

    def cuts(self, mode: int) -> list[arb]:
        half_steps = 2 * WINDOW_XI * 2
        cuts = [arb(-WINDOW_XI) + arb(index) / 2 for index in range(half_steps + 1)]
        cuts.append(self.crossing_xi(mode))
        split = self.split_xi(mode)
        if -arb(WINDOW_XI) < split < arb(WINDOW_XI):
            cuts.append(split)
        return sorted(cuts, key=float)

    def regimes(self, mode: int, left: arb, right: arb) -> tuple[str, int]:
        center = (left + right) / 2
        x = self.x0 + center / self.root_hessian
        split = self.split_x(mode)
        normal = (mode == 621 and x <= split) or (mode == 622 and x >= split)
        q = (x / 2).sqrt() * (self.endpoint - 2 * arb(mode) / x)
        side = -1 if q < 0 else 1
        return ("normal_safe" if normal else "outer_safe"), side

    def direct_terms(self, mode: int, x: acb) -> tuple[acb, acb]:
        m = arb(mode)
        c = self.endpoint * x / 2
        positive = (
            self.endpoint / (self.i * self.pi * (self.endpoint * x - 2 * m))
            - m / (2 * self.pi**2 * (c - m) ** 3)
            + 3 * self.i * m * x / (4 * self.pi**3 * (c - m) ** 5)
        )
        negative = (
            self.endpoint / (self.i * self.pi * (self.endpoint * x + 2 * m))
            + m / (2 * self.pi**2 * (c + m) ** 3)
            - 3 * self.i * m * x / (4 * self.pi**3 * (c + m) ** 5)
        )
        return positive, negative

    def exact_positive(self, mode: int, x: acb, side: int) -> acb:
        m = arb(mode)
        q = (x / 2).sqrt() * (self.endpoint - 2 * m / x)
        radius = side * q
        exponential = (self.i * self.pi * q**2 / 2).exp()
        argument = (-self.i * self.pi / 4).exp() * (self.pi / 2).sqrt() * radius
        tail = self.full_fresnel * argument.erfc()
        amplitude = m * (2 / x).sqrt()
        sign_adapted = exponential / (self.i * self.pi) - side * amplitude * tail
        step_coefficient = (1 if mode <= 621 else 0) - (1 if side < 0 else 0)
        step = arb(step_coefficient) * amplitude * (1 + self.i)
        return (sign_adapted + step) / (x * exponential)

    def integrand(self, mode: int, regime: str, side: int):
        def function(xi: acb, _: bool) -> acb:
            x = self.x0 + xi / self.root_hessian
            positive, negative = self.direct_terms(mode, x)
            if regime == "normal_safe":
                coefficient = positive + negative
            else:
                coefficient = self.exact_positive(mode, x, side) + negative
            phase = (
                self.pi * self.endpoint**2 * (x - self.x0) / 4
                + self.t * (((1 - x) / x).log() - self.log0) / 2
            )
            weight = (x * (1 - x)) ** (-arb(1) / 4)
            return weight * (self.i * phase).exp() * coefficient / self.root_hessian

        return function


def cache_key(mode: int, index: int) -> str:
    return f"{mode}:{index}"


def load_cache(window: LocalWindow) -> dict[str, dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {}
    if not CACHE.is_file():
        return rows
    cuts_by_mode = {mode: window.cuts(mode) for mode in MODES}
    for line in CACHE.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        mode, index = row.get("mode"), row.get("panel_index")
        if mode not in cuts_by_mode or not isinstance(index, int):
            continue
        cuts = cuts_by_mode[mode]
        if (
            row.get("formula_version") == FORMULA_VERSION
            and 0 <= index < len(cuts) - 1
            and row.get("xi_left") == str(cuts[index])
            and row.get("xi_right") == str(cuts[index + 1])
        ):
            rows[cache_key(mode, index)] = row
    return rows


def append_cache(row: dict[str, Any]) -> None:
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    with CACHE.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(row, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def integrate(window: LocalWindow) -> tuple[dict[int, acb], list[dict[str, Any]]]:
    rows = load_cache(window)
    total_panels = sum(len(window.cuts(mode)) - 1 for mode in MODES)
    completed = len(rows)
    for mode in MODES:
        cuts = window.cuts(mode)
        for index, (left, right) in enumerate(zip(cuts, cuts[1:])):
            key = cache_key(mode, index)
            if key in rows:
                continue
            regime, side = window.regimes(mode, left, right)
            started = time.time()
            value = acb.integral(
                window.integrand(mode, regime, side),
                left,
                right,
                abs_tol=arb(TOLERANCE),
                rel_tol=arb(TOLERANCE),
                eval_limit=250_000,
                depth_limit=40,
            )
            require("nan" not in str(value).lower(), f"nonconvergent local panel mode={mode} index={index}")
            row = {
                "formula_version": FORMULA_VERSION,
                "mode": mode,
                "panel_index": index,
                "xi_left": str(left),
                "xi_right": str(right),
                "regime": regime,
                "q_side": side,
                "value": complex_record(value),
                "elapsed_seconds": round(time.time() - started, 3),
            }
            append_cache(row)
            rows[key] = row
            completed += 1
            if completed % 10 == 0 or completed == total_panels:
                print(f"cached local B panels {completed}/{total_panels}", flush=True)

    require(len(rows) == total_panels, "incomplete local-current cache")
    ordered: list[dict[str, Any]] = []
    totals: dict[int, acb] = {}
    for mode in MODES:
        cuts = window.cuts(mode)
        mode_rows = [rows[cache_key(mode, index)] for index in range(len(cuts) - 1)]
        ordered.extend(mode_rows)
        totals[mode] = sum((parse_complex(row["value"]) for row in mode_rows), acb(0))
    return totals, ordered


def negative_tail_remainder(window: LocalWindow) -> arb:
    x_low = window.x0 - arb(WINDOW_XI) / window.root_hessian
    x_high = window.x0 + arb(WINDOW_XI) / window.root_hessian
    density = arb(0)
    for mode in MODES:
        m = arb(mode)
        density += 15 * m * x_high**2 / (4 * window.pi**4 * (window.endpoint * x_low / 2 + m) ** 7)
    physical = (
        2
        * window.normalizer
        * (x_high - x_low)
        * (x_low * (1 - x_low)) ** (-arb(1) / 4)
        * density
    )
    require(physical < arb("8.5e-36"), "negative local tail remainder exceeds 8.5e-36")
    return physical


def numerical_certificate(window: LocalWindow, local_remainder: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    totals, rows = integrate(window)
    reduced = totals[621] + totals[622]
    rotation = (window.i * window.phase0 - window.i * window.pi / 8).exp()
    physical_by_mode = {mode: 2 * window.normalizer * (rotation * value).real for mode, value in totals.items()}
    principal = physical_by_mode[621] + physical_by_mode[622]
    positive_remainder = arb(local_remainder["combined_physically_normalized_remainder_ball"])
    negative_remainder = negative_tail_remainder(window)
    exact_upper = principal + positive_remainder + negative_remainder

    require(arb("-1.39e-4") < principal < arb("-1.37e-4"), "local principal projection escaped scout bracket")
    require(exact_upper < arb("-1.369e-4"), "exact local B-window upper bound did not remain negative")
    return {
        "height": T,
        "endpoint": B,
        "window_xi": [-WINDOW_XI, WINDOW_XI],
        "panel_width_xi": PANEL_WIDTH_XI,
        "panel_count": len(rows),
        "precision_decimal_digits": PRECISION,
        "quadrature_tolerance_per_panel": TOLERANCE,
        "crossing_xi_balls": {str(mode): window.crossing_xi(mode).str(PRECISION, more=True) for mode in MODES},
        "balance_split_xi_balls": {str(mode): window.split_xi(mode).str(PRECISION, more=True) for mode in MODES},
        "mode_reduced_integral_balls": {str(mode): complex_record(value) for mode, value in totals.items()},
        "mode_physical_signed_projection_balls": {str(mode): value.str(PRECISION, more=True) for mode, value in physical_by_mode.items()},
        "combined_principal_physical_signed_projection_ball": principal.str(PRECISION, more=True),
        "positive_normal_remainder_absolute_allowance_ball": positive_remainder.str(PRECISION, more=True),
        "negative_partner_remainder_absolute_allowance_ball": negative_remainder.str(PRECISION, more=True),
        "exact_local_B_window_physical_upper_ball": exact_upper.str(PRECISION, more=True),
        "exact_local_B_window_physical_upper_bound": "-1.369e-4",
    }, rows


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    return f"""# Exact outer-safe local B current for modes 621 and 622

Date: 2026-08-13

Status: cancellation-safe saved-height interval certificate; not a complete
B-face estimate

For a positive mode, retain the exact continuous target-step/tail current

```text
R_B,m^+=U_B,m+[1_(m<=621)-1_(q_B<0)]P_bulk,m.        (LC1)
```

Because `B` is odd, quadratic completion gives its exact face coefficient

```text
C_exact,+,m=exp(-i*pi*q_B^2/2)R_B,m^+/x.             (LC2)
```

The jump of `U_B` at `q_B=0` is cancelled by the step in LC1 before any
integral or absolute value is taken.  The negative partner is uniformly far
from its crossing and is represented by its first three direct Fresnel
terms.  On a positive normal-safe arc those same three terms are used for
both signs; on the complementary arc LC2 replaces only the positive terms.

The mode-621 balance point lies at
`xi={c['balance_split_xi_balls']['621']}`, below the Gaussian window, so
mode 621 uses the exact outer-safe representation throughout `|xi|<=70`.
Mode 622 switches to its normal-safe direct expansion at
`xi={c['balance_split_xi_balls']['622']}`.  Every panel is also cut at the
exact mode crossing:

```text
xi_621={c['crossing_xi_balls']['621']},
xi_622={c['crossing_xi_balls']['622']}.               (LC3)
```

An append-only `{c['panel_count']}`-panel Arb quadrature evaluates the two
mode sum before projection.  It gives

```text
mode 621: {c['mode_physical_signed_projection_balls']['621']},
mode 622: {c['mode_physical_signed_projection_balls']['622']},

principal sum:
{c['combined_principal_physical_signed_projection_ball']}.          (LC4)
```

The existing positive normal-safe remainder allowance is
`{c['positive_normal_remainder_absolute_allowance_ball']}`.  The negative
partners contribute at most
`{c['negative_partner_remainder_absolute_allowance_ball']}` beyond their
three displayed terms.  Adding both adversely still proves

```text
E_local,B <= {c['exact_local_B_window_physical_upper_ball']}
           < -1.369e-4.                              (LC5)
```

Thus the local pair is not a positive budget cost: its exact physical
projection is rigorously negative by more than `1.369e-4`.  This sign is
obtained only after the crossing step and Fresnel tail remain joined.

Pi provenance: all `pi` factors come from the equation-(9) Fresnel and
outer phases, odd-endpoint completion, or the paper normalization.  No
fitted constant is used.

Proof boundary: modes 621 and 622 on the saved-height `|xi|<=70` B window
only, with a deliberately overinclusive positive normal-remainder
allowance.  No outside-window grouped tail, complete B estimate, A-fold
splice, complete `T_upper`, height-uniform theorem, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    ctx.dps = PRECISION
    ctx.threads = 1
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["step_tail"]["decision"]["tail_and_step_jumps_cancel_exactly"] is True, "step-tail dependency drift")
    require(dependencies["B_window"]["decision"]["xi_70_window_contains_exactly_B_crossings_621_622"] is True, "B-window dependency drift")
    require(dependencies["local_normal_remainder"]["decision"]["combined_physical_remainder_below_1p97e_minus_10"] is True, "local-remainder dependency drift")
    require(dependencies["signed_background"]["decision"]["first_pole_subtracted_current_integrated_before_absolute_value"] is True, "signed-orientation dependency drift")
    require(dependencies["pair_triangle"]["decision"]["pair_Volterra_integral_reduced_to_separable_triangle"] is True, "pair-triangle dependency drift")
    require(dependencies["truncation_dictionary"]["decision"]["exact_rational_dictionary_correction_derived"] is True, "dictionary dependency drift")

    window = LocalWindow()
    certificate, rows = numerical_certificate(
        window,
        dependencies["local_normal_remainder"]["interval_certificate"],
    )
    artifact = {
        "kind": STEM,
        "status": "exact_local_621_622_B_window_physical_projection_below_minus_1p369e_minus_4",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": certificate,
        "decision": {
            "target_step_and_Fresnel_tail_kept_joined": True,
            "panels_cut_at_both_B_crossings": True,
            "mode_621_balance_split_below_window": True,
            "mode_622_normal_outer_splice_certified": True,
            "negative_partner_tail_remainder_bounded": True,
            "exact_local_621_622_window_projection_below_minus_1p369e_minus_4": True,
            "outside_window_grouped_trace_tails_bounded": False,
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
        "next_obligation": "Join this signed negative local theorem to the pole-subtracted nonlocal B-window currents and all certified remainders, then bound the two grouped |xi|>70 B trace tails without separating endpoint completion sectors.",
        "proof_boundary": "Only the exact cancellation-safe modes 621 and 622 on the saved-height |xi|<=70 B window. No outside-window tail, complete B estimate, A-fold splice, complete paired residual, complete T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level theorem is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified exact local B-window projection below -1.369e-4", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
