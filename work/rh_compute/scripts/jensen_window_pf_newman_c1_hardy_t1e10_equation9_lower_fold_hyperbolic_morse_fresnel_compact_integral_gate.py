#!/usr/bin/env python3
"""Certify the signed compact Morse--Fresnel/Airy overlap integral."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
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

from flint import acb, arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_compact_integral_gate"
OVERLAP_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_overlap_core_gate.json"
LOGISTIC_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_fresnel_retaining_logistic_transform_gate.json"
AIRY_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_grouped_characteristic_airy_boundary_normal_form_gate.json"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
CACHE = REPO_ROOT / f"work/rh_compute/results/{STEM}_cache.jsonl"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 90
TOLERANCE = "1e-25"
C = 159_577
MODE = 39_894
Y_CUTOFF = 64
CORE_U = 200
PANELS_PER_HALF = 32
SERIES_DEGREE = 24
SERIES_RADIUS = "0.01"
FORMULA_VERSION = "hyperbolic_compact_difference_v1"
TARGET_ABSOLUTE = arb("0.00049")
TARGET_PHYSICAL = arb("0.000221")

ROOTS = {
    "lower_face": (
        "-0.002825759532961365076401160417159114247825864097320687324089683704073666283224381896263434854721777473",
        "0.002831092859949456681693112196296051440370287003489608130887631666653078012463105775307871867951718149",
    ),
    "event": (
        "-0.002825759532933649445131464450396983743915033772561613810770360267229064685468657752846686488110707596",
        "0.002831092859921636330917475537313707362085208216652466697578982972472672141177079747032833042459072342",
    ),
    "upper_face": (
        "-0.002825759532905933813862584261780894515755158231241335126296618978417123373928547289289614360234877368",
        "0.002831092859893815980142658768797893714991216035497558963491765267230760723646577968637109603399278758",
    ),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def append_cache(row: dict[str, Any]) -> None:
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    with CACHE.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(row, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def load_cache() -> dict[tuple[str, int], dict[str, Any]]:
    rows: dict[tuple[str, int], dict[str, Any]] = {}
    if not CACHE.is_file():
        return rows
    for line in CACHE.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        if row.get("formula_version") != FORMULA_VERSION:
            continue
        if "nan" in row["value"]["real_ball"].lower() or "nan" in row["value"]["imag_ball"].lower():
            continue
        rows[(row["face"], int(row["panel_index"]))] = row
    return rows


def symbolic_certificate() -> dict[str, str]:
    a, b, y = sp.symbols("a b y", nonzero=True)
    h0 = sp.Function("H0")(a, b, y)
    endpoint = sp.exp(sp.I * (a * y**2 + b * y))
    h1 = -b * h0 / (2 * a) - sp.I * (endpoint - 1) / (2 * a)
    require(
        sp.simplify(2 * a * h1 + b * h0 + sp.I * (endpoint - 1)) == 0,
        "normal first-moment identity failed",
    )
    q, q0, u, rho, x, Y = sp.symbols("q q0 u rho x Y")
    exact = ((q + rho * Y) ** 2 - q0**2 - u**2) / 2
    expanded = (q**2 - q0**2 - u**2) / 2 + rho * q * Y + x * Y**2
    require(sp.simplify(sp.expand(exact - expanded).subs(rho**2, 2 * x)) == 0, "exact Y-phase expansion failed")
    return {
        "exact_inner_phase": "B_t(u)+sqrt(2*x)*Q_C*Y+x*Y^2",
        "canonical_inner_phase": "B_0(u)-u*Y+Y^2/2",
        "normal_integral": "H0(a,b;Y)=int_0^Y exp(i*(a*y^2+b*y))dy",
        "normal_first_moment": "H1=-b*H0/(2*a)-i*(exp(i*(a*Y^2+b*Y))-1)/(2*a)",
        "signed_compact_difference": "int_[|u|<=200] [exp(i*Delta_c)A_t exp(i*Phi_t)-exp(i*Phi_0)]dYdu",
    }


class CompactFace:
    def __init__(self, label: str, offset_sign: int) -> None:
        self.label = label
        self.i = acb(0, 1)
        self.pi = arb.pi()
        self.c = arb(C)
        self.m = arb(MODE)
        self.tstar = self.pi * self.c**2 / 8
        self.beta = self.tstar ** (arb(1) / 3)
        self.tau = self.pi * self.m * (self.c - 2 * self.m)
        self.t = self.tau + arb(offset_sign) * self.pi / 16
        self.ymax = arb(Y_CUTOFF) / (2 * self.beta).sqrt()
        self.r = self.t / (2 * self.pi * self.m**2)
        self.d0 = 1 + self.r
        self.a0 = (self.tau - self.t) / (self.pi * self.m)
        self.b0 = self.t / (self.pi * self.m)
        self.qzero = self.a0 * (self.pi / (2 * self.d0)).sqrt()
        self.ell = 2 * (self.t - self.tau) / (self.c * self.pi.sqrt())
        self.f0_constant = (self.i * self.pi / 4).exp() * (self.pi / 2).sqrt()
        self.f0_rotation = (-self.i * self.pi / 4).exp() / arb(2).sqrt()
        self.carrier_mismatch = self._carrier_mismatch()
        self.carrier = (self.i * self.carrier_mismatch).exp()
        self.delta_lo = arb(ROOTS[label][0])
        self.delta_hi = arb(ROOTS[label][1])

    def _carrier_mismatch(self) -> arb:
        z0 = self.beta * self.r.log() / 2
        eta = self.tstar / self.t
        exact = self.t * (z0 / self.beta - eta * (z0 / self.beta).tanh())
        d = -self.beta / self.c
        lam = (self.tstar - self.t) / self.beta
        canonical = lam * d - d**3 / 3
        return exact - canonical

    def f0(self, q: acb) -> acb:
        return self.f0_constant * (self.f0_rotation * q).erf()

    def normal_integral(self, a: acb, b: acb) -> acb:
        root = (2 * a).sqrt()
        q0 = b / root
        q1 = root * self.ymax + b / root
        return (-self.i * b**2 / (4 * a)).exp() / root * (self.f0(q1) - self.f0(q0))

    def first_normal_moment(self, a: acb, b: acb, h0: acb) -> acb:
        endpoint = (self.i * (a * self.ymax**2 + b * self.ymax)).exp()
        return -b * h0 / (2 * a) - self.i * (endpoint - 1) / (2 * a)

    def k_series(self, delta: acb) -> acb:
        total = acb(0)
        for power in range(SERIES_DEGREE, -1, -1):
            coefficient = arb(1 if power % 2 == 0 else -1) / (power + 2)
            total = total * delta + coefficient
        radius = arb(SERIES_RADIUS)
        tail = radius ** (SERIES_DEGREE + 1) / ((SERIES_DEGREE + 3) * (1 - radius))
        return total + acb(arb(0, tail), arb(0, tail))

    def components(self, delta: acb) -> tuple[acb, acb, acb]:
        k = self.k_series(delta)
        denominator = self.d0 + self.r * delta
        x = 1 / denominator
        u = delta * (self.t * k).sqrt()
        qc = (self.a0 - self.b0 * delta) * (self.pi / (2 * denominator)).sqrt()

        inside = (
            -self.pi * self.d0 * self.a0 * self.b0
            - self.pi * (self.a0 * self.a0) * self.r / 2
            + self.t * self.d0 * delta * (self.r - k * denominator)
        )
        exact_boundary = delta * inside / (2 * self.d0 * denominator)
        rho = (2 * x).sqrt()
        g = (1 + delta) * (2 * k).sqrt()
        base = self.r ** (arb(1) / 4) * (1 + delta) ** (-arb(1) / 4) * g * rho
        epsilon = (self.r * x / self.t).sqrt()

        h0_exact = self.normal_integral(x, rho * qc)
        h1_exact = self.first_normal_moment(x, rho * qc, h0_exact)
        exact_inner = base * ((1 + qc * epsilon) * h0_exact + rho * epsilon * h1_exact)
        exact = (self.i * exact_boundary).exp() * exact_inner

        canonical_boundary = self.ell * u + u**2 / (2 * self.c) + u**3 / (3 * self.c * self.pi.sqrt())
        canonical = (self.i * canonical_boundary).exp() * self.normal_integral(acb(arb(1) / 2), -u)
        jacobian = (self.t / 2).sqrt() / g
        return self.carrier * exact * jacobian, canonical * jacobian, u

    def difference_integrand(self, delta: acb, analytic: bool) -> acb:
        exact, canonical, _ = self.components(delta)
        return exact - canonical

    def real_u(self, delta: arb) -> arb:
        value = self.k_series(acb(delta))
        require(value.imag.contains(0), "real k-series acquired an imaginary component")
        return delta * (self.t * value.real).sqrt()


def panel_bounds(face: CompactFace, index: int) -> tuple[arb, arb]:
    if index < PANELS_PER_HALF:
        left = face.delta_lo + (0 - face.delta_lo) * index / PANELS_PER_HALF
        right = face.delta_lo + (0 - face.delta_lo) * (index + 1) / PANELS_PER_HALF
        return left, right
    local = index - PANELS_PER_HALF
    left = face.delta_hi * local / PANELS_PER_HALF
    right = face.delta_hi * (local + 1) / PANELS_PER_HALF
    return left, right


def integrate_panel(face: CompactFace, index: int) -> acb:
    left, right = panel_bounds(face, index)
    return acb.integral(
        face.difference_integrand,
        left,
        right,
        abs_tol=arb(TOLERANCE),
        rel_tol=arb(TOLERANCE),
        eval_limit=400_000,
        depth_limit=45,
    )


def fill_cache(faces: list[CompactFace], rows: dict[tuple[str, int], dict[str, Any]]) -> None:
    panel_count = 2 * PANELS_PER_HALF
    for face in faces:
        for index in range(panel_count):
            key = (face.label, index)
            if key in rows:
                continue
            started = time.perf_counter()
            value = integrate_panel(face, index)
            left, right = panel_bounds(face, index)
            row = {
                "formula_version": FORMULA_VERSION,
                "face": face.label,
                "panel_index": index,
                "delta_left": str(left),
                "delta_right": str(right),
                "value": complex_record(value),
                "precision_decimal_digits": PRECISION,
                "tolerance": TOLERANCE,
                "elapsed_seconds": round(time.perf_counter() - started, 3),
            }
            append_cache(row)
            rows[key] = row
            print(f"cached {face.label} panel {index + 1}/{panel_count}: {row['elapsed_seconds']:.3f}s", flush=True)


def certify_face(face: CompactFace, rows: dict[tuple[str, int], dict[str, Any]]) -> dict[str, Any]:
    panel_count = 2 * PANELS_PER_HALF
    require(all((face.label, index) in rows for index in range(panel_count)), f"incomplete {face.label} cache")
    total = sum((parse_complex(rows[(face.label, index)]["value"]) for index in range(panel_count)), acb(0))

    u_lo = face.real_u(face.delta_lo)
    u_hi = face.real_u(face.delta_hi)
    endpoint_error = abs(u_lo + CORE_U).upper() + abs(u_hi - CORE_U).upper()
    sliver_error = 3 * endpoint_error
    target = acb(arb(total.real, sliver_error), arb(total.imag, sliver_error))
    physical = arb(2).sqrt() / face.pi * target
    require(abs(target) < TARGET_ABSOLUTE, f"{face.label} compact difference exceeds 0.00049")
    require(abs(physical) < TARGET_PHYSICAL, f"{face.label} physical compact difference exceeds 0.000221")
    require(abs(face.carrier_mismatch) < arb("2e-12"), f"{face.label} carrier mismatch lost its bound")
    require(endpoint_error < arb("2e-49"), f"{face.label} endpoint coordinate error too large")
    return {
        "label": face.label,
        "height_ball": face.t.str(PRECISION, more=True),
        "delta_lower_endpoint": str(face.delta_lo),
        "delta_upper_endpoint": str(face.delta_hi),
        "lower_u_ball": u_lo.str(PRECISION, more=True),
        "upper_u_ball": u_hi.str(PRECISION, more=True),
        "endpoint_coordinate_error_ball": endpoint_error.str(PRECISION, more=True),
        "endpoint_sliver_error_ball": sliver_error.str(PRECISION, more=True),
        "carrier_phase_mismatch_ball": face.carrier_mismatch.str(PRECISION, more=True),
        "normalized_compact_difference_ball": complex_record(target),
        "normalized_compact_difference_absolute_ball": abs(target).str(PRECISION, more=True),
        "physical_compact_difference_ball": complex_record(physical),
        "physical_compact_difference_absolute_ball": abs(physical).str(PRECISION, more=True),
        "panels": panel_count,
    }


def render_note(artifact: dict[str, Any]) -> str:
    rows = artifact["certified_faces"]
    lines = []
    for row in rows:
        lines.append(
            f"{row['label']}: |Delta I_core|={row['normalized_compact_difference_absolute_ball']}, "
            f"physical={row['physical_compact_difference_absolute_ball']}"
        )
    face_text = "\n".join(lines)
    return f"""# Signed hyperbolic Morse--Fresnel compact integral

Date: 2026-08-12

Status: cancellation-preserving compact integral certified at three heights;
this is not a proof of the outer tails or complete one-mode join

For mode `39894`, retain the common carrier and paired endpoint/Fresnel
current.  On `|u|<=200` and `0<=Y<=64/sqrt(2 beta)`, compare

```text
Delta I_core(t)=int int [exp(i Delta_c) A_t(u,Y) exp(i Phi_t(u,Y))
                         -exp(i Phi_0(u,Y))] dY du.       (CI1)
```

The `Y` phase is exactly quadratic and `A_t` is exactly affine in `Y`.
Writing

```text
H_0(a,b;Y)=int_0^Y exp(i[a y^2+b y])dy,
H_1(a,b;Y)=-b H_0/(2a)
            -i[exp(i[aY^2+bY])-1]/(2a),                (CI2)
```

evaluates the inner integral in closed Fresnel form.  The remaining outer
integral is performed directly on the signed difference in `delta=v-1`,
with `du/delta=sqrt(t/2)/[(1+delta)sqrt(2k(delta))]` and

```text
k(delta)=[delta-log(1+delta)]/delta^2.                 (CI3)
```

A degree-{SERIES_DEGREE} convergent series with an explicit complex geometric
tail resolves the removable value at `delta=0`.  Append-only interval panels
give

```text
{face_text}
```

uniformly at `t=tau-pi/16`, `tau`, and `tau+pi/16`.  In particular,

```text
|Delta I_core| < 0.00049,
(sqrt(2)/pi)|Delta I_core| < 0.000221.                 (CI4)
```

The displayed `pi` comes from the original Kummer quadratic phase and
Fourier--Poisson character; `sqrt(2)/pi` is the exact coordinate/Jacobian
factor already derived in the overlap chart.

## Boundary

This gate certifies only the signed compact rectangle for one prototype mode
at three heights.  It does not bound `|u|>200`, complete the one-mode
Airy/logistic integral join, propagate through 399 events, establish complete
`T_upper`, or prove `Lambda<=0`, RH, or a prize-level conclusion.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuild-cache", action="store_true")
    args = parser.parse_args()
    ctx.dps = PRECISION
    priority = set_low_priority()
    require(OVERLAP_GATE.is_file(), "missing overlap-core dependency")
    require(LOGISTIC_GATE.is_file(), "missing logistic-transform dependency")
    require(AIRY_GATE.is_file(), "missing Airy-normal-form dependency")
    require(CHECKER.is_file(), "missing independent checker")
    if args.rebuild_cache and CACHE.exists():
        CACHE.unlink()

    symbolic = symbolic_certificate()
    faces = [CompactFace("lower_face", -1), CompactFace("event", 0), CompactFace("upper_face", 1)]
    rows = load_cache()
    fill_cache(faces, rows)
    certified = [certify_face(face, rows) for face in faces]

    artifact = {
        "kind": STEM,
        "date": "2026-08-12",
        "status": "signed_compact_hyperbolic_Morse_Fresnel_overlap_integral_certified",
        "precision_decimal_digits": PRECISION,
        "resource_policy": {
            "worker_cap": 1,
            "active_workers": 1,
            "process_priority": priority,
            "panels_per_face": 2 * PANELS_PER_HALF,
            "append_only_cache": True,
        },
        "symbolic_certificate": symbolic,
        "parameters": {
            "C": C,
            "mode": MODE,
            "core_u": CORE_U,
            "y_cutoff": Y_CUTOFF,
            "series_degree": SERIES_DEGREE,
            "series_radius": SERIES_RADIUS,
        },
        "certified_faces": certified,
        "claims": {
            "inner_Fresnel_integral_exact": True,
            "paired_current_retained": True,
            "common_carrier_retained": True,
            "signed_compact_difference_integrated_before_absolute_values": True,
            "normalized_compact_difference_below_0_00049": True,
            "physical_compact_difference_below_0_000221": True,
            "outer_u_tail_proved": False,
            "complete_one_mode_join_proved": False,
            "RH_proved": False,
        },
        "proof_boundary": (
            "Rigorous signed compact integral for one mode at three heights only. No |u|>200 tail, "
            "complete one-mode join, 399-event propagation, complete T_upper, Lambda<=0, RH, or prize-level conclusion is proved."
        ),
        "next_action": (
            "Derive compatible left and right |u|>200 contour or integration-by-parts bounds in the same retained-Fresnel chart, "
            "including the cutoff boundary terms, before testing event-index propagation."
        ),
        "dependencies": {
            "overlap_gate": {"path": relative(OVERLAP_GATE), "sha256": file_hash(OVERLAP_GATE)},
            "logistic_gate": {"path": relative(LOGISTIC_GATE), "sha256": file_hash(LOGISTIC_GATE)},
            "Airy_gate": {"path": relative(AIRY_GATE), "sha256": file_hash(AIRY_GATE)},
        },
        "artifacts": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "cache": {"path": relative(CACHE), "sha256": file_hash(CACHE)},
        },
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    artifact["artifacts"]["note"] = {"path": relative(NOTE), "sha256": file_hash(NOTE)}
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    maximum = max(arb(row["normalized_compact_difference_absolute_ball"]) for row in certified)
    print(
        f"certified signed compact overlap integral: faces=3, max |difference|={maximum}; priority={priority}",
        flush=True,
    )


if __name__ == "__main__":
    main()
