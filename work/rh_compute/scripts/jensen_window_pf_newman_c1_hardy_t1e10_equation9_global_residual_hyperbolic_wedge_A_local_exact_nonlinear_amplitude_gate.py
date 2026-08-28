#!/usr/bin/env python3
"""Certify the compact exact-minus-affine transformed amplitude on A."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_local_exact_nonlinear_amplitude_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
CACHE = REPO_ROOT / f"work/rh_compute/results/{STEM}_rows.jsonl"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "factorization": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_remainder_identity_gate.json",
    "local_affine": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_A_local_exact_domain_gate.json",
    "orientation": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_corner_projector_orientation_gate.json",
    "erfc_reduction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_exact_exterior_common_phase_reduction_gate.json",
    "exact_exterior": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_exact_exterior_rational_common_phase_contour_gate.json",
}

HEIGHT = 10_000_000_000
ENDPOINT = 159_577
LOWER_START = 39_853
TARGET_END = 39_894
UPPER_END = 39_936
MODES = tuple(range(LOWER_START, UPPER_END + 1))
Y_MAX_TEXT = "0.0037"
PRODUCTION_PRECISION = 60
PRODUCTION_SERIES_TERMS = 24
PRODUCTION_TOLERANCE = "1e-15"
CACHE_VERSION = 1
ROW_FORMULA_VERSION = "local-exact-minus-affine-erfc-center-ode-disk-v1"


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


def complex_record(value: acb, digits: int = 55) -> dict[str, str]:
    return {
        "real_ball": value.real.str(digits, more=True),
        "imag_ball": value.imag.str(digits, more=True),
        "absolute_ball": abs(value).str(digits, more=True),
    }


def config_signature(precision: int, series_terms: int, tolerance: str) -> str:
    payload = {
        "cache_version": CACHE_VERSION,
        "row_formula_version": ROW_FORMULA_VERSION,
        "height": HEIGHT,
        "endpoint": ENDPOINT,
        "mode_range": [LOWER_START, UPPER_END],
        "y_max": Y_MAX_TEXT,
        "precision": precision,
        "series_terms": series_terms,
        "tolerance": tolerance,
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("ascii")).hexdigest()


def symbolic_normalization_certificate() -> dict[str, str]:
    x, mode, endpoint = sp.symbols("x m A", positive=True, real=True)
    K = sp.pi * endpoint * mode
    u = endpoint * x / (2 * mode)
    P = sp.sqrt(sp.pi / 2) * (endpoint * sp.sqrt(x) - 2 * mode / sp.sqrt(x))
    endpoint_factor = 2 * u * (K * u - sp.I) / ((u + 1) * (K - sp.I))
    normalization = 2 * sp.sqrt(2) * sp.I * (K - sp.I) / (sp.sqrt(sp.pi) * endpoint)
    driver = x ** (-sp.Rational(1, 2)) * (2 + sp.I * sp.pi * endpoint**2 * x)
    phase = sp.pi * mode**2 / x + sp.pi * endpoint**2 * x / 4
    require(
        sp.simplify(normalization * endpoint_factor * sp.diff(P, x) - driver) == 0,
        "endpoint-current normalization identity failed",
    )
    require(
        sp.simplify(P**2 / 2 - phase + sp.pi * endpoint * mode) == 0,
        "endpoint-current phase identity failed",
    )
    return {
        "endpoint_current": "I_E(P_A(x))=(-1)^(A*m) H_m(x)/C_E",
        "normalization": "C_E=2*sqrt(2)*i*(pi*A*m-i)/(sqrt(pi)*A)",
        "driver_identity": "C_E*E(P_A(x))*dP_A/dx=x^(-1/2)(2+i*pi*A^2*x)",
        "phase_identity": "P_A(x)^2/2=pi*m^2/x+pi*A^2*x/4-pi*A*m",
        "odd_endpoint_parity": "A=159577 is odd, hence (-1)^(A*m)=(-1)^m",
    }


def q_and_derivative(y: acb, series_terms: int) -> tuple[acb, acb]:
    """Enclose q(y)=2(y-log(1+y))/y^2 and its derivative at y=0 too."""
    q = acb(0)
    q_prime = acb(0)
    power = acb(1)
    derivative_power = acb(1)
    for k in range(series_terms):
        coefficient = arb(2 * (-1 if k % 2 else 1)) / arb(k + 2)
        q += coefficient * power
        if k:
            q_prime += arb(k) * coefficient * derivative_power
            derivative_power *= y
        power *= y
    radius = abs(y).upper()
    require(radius < arb("0.01"), "logistic series used outside its certified disk")
    q_tail = arb(2) * radius**series_terms / (arb(series_terms + 2) * (1 - radius))
    q_prime_tail = arb(2) * radius ** (series_terms - 1) / (1 - radius)
    q += acb(arb(0, q_tail), arb(0, q_tail))
    q_prime += acb(arb(0, q_prime_tail), arb(0, q_prime_tail))
    return q, q_prime


def mode_constants(mode_int: int) -> dict[str, Any]:
    pi, t, endpoint = arb.pi(), arb(HEIGHT), arb(ENDPOINT)
    mode = arb(mode_int)
    r = t / (2 * pi * mode**2)
    K = pi * endpoint * mode
    return {
        "mode_int": mode_int,
        "mode": mode,
        "r": r,
        "K": K,
        "lambda_P": acb(3 * K, -1) / (acb(2 * K, -2) * K.sqrt()),
        "mu_S": arb(5) / 12 * (arb(2) / t).sqrt(),
        "y_half": 2 * pi * mode**2 / t - 1,
        "W0": acb(K, -1) * (4 * r ** arb("0.75") / (endpoint * (pi * t).sqrt())),
        "carrier": (
            acb(0, 1)
            * (t * ((t / (2 * pi)).log() - 1) / 2 - t * mode.log())
        ).exp(),
    }


def exact_endpoint_current(x: acb, constants: dict[str, Any], analytic: bool) -> acb:
    """Return int_-infinity^P_A E(P) exp(iP^2/2) dP by exact erfc reduction."""
    pi = arb.pi()
    imaginary = acb(0, 1)
    endpoint = arb(ENDPOINT)
    mode = constants["mode"]
    rotation = (-imaginary * pi / 4).exp()
    X = x.sqrt(analytic=analytic)
    alpha = rotation * pi.sqrt() * endpoint / 2
    beta = rotation * pi.sqrt() * mode
    phase = pi * mode**2 / x + pi * endpoint**2 * x / 4
    w_minus = beta / X - alpha * X
    w_plus = beta / X + alpha * X
    parity = -1 if constants["mode_int"] % 2 else 1
    H = 4 * X * (imaginary * phase).exp()
    H -= 2 * pi.sqrt() * beta * parity * (w_minus.erfc() + w_plus.erfc())
    normalization = 2 * arb(2).sqrt() * imaginary * (constants["K"] - imaginary) / (pi.sqrt() * endpoint)
    return parity * H / normalization


def centered_current_enclosure(center_value: acb, radius: arb, derivative_bound: arb) -> acb:
    variation = radius * derivative_bound
    return center_value + acb(arb(0, variation), arb(0, variation))


def local_geometry(y: acb, constants: dict[str, Any], series_terms: int, analytic: bool) -> dict[str, acb]:
    pi, t = arb.pi(), arb(HEIGHT)
    q, q_prime = q_and_derivative(y, series_terms)
    sqrt_q = q.sqrt(analytic=analytic)
    v = 1 + y
    S = (t / 2).sqrt() * y * sqrt_q
    dS_dy = (t / 2).sqrt() * (sqrt_q + y * q_prime / (2 * sqrt_q))
    outer_factor = v ** arb("0.75") * sqrt_q
    g = 1 + constants["r"] * v
    sqrt_g = g.sqrt(analytic=analytic)
    P_exact = (pi / 2).sqrt() * (arb(ENDPOINT) / sqrt_g - 2 * constants["mode"] * sqrt_g)
    dP_dy = -(pi / 2).sqrt() * constants["r"] * (
        arb(ENDPOINT) / (2 * g * sqrt_g) + constants["mode"] / sqrt_g
    )
    return {
        "v": v,
        "S": S,
        "dS_dy": dS_dy,
        "outer_factor": outer_factor,
        "g": g,
        "x": 1 / g,
        "P": P_exact,
        "dP_dy": dP_dy,
    }


def exact_minus_affine_integrand(constants: dict[str, Any], series_terms: int):
    pi, t = arb.pi(), arb(HEIGHT)
    imaginary = acb(0, 1)
    root_two = arb(2).sqrt()
    fresnel_constant = (pi / 2).sqrt() * (imaginary * pi / 4).exp()
    erf_rotation = (-imaginary * pi / 4).exp() / root_two

    def integrand(y: acb, analytic: bool) -> acb:
        if abs(y).upper() >= arb("0.01"):
            if analytic:
                return acb("nan")
            raise RuntimeError("quadrature point escaped the certified logistic disk")
        geometry = local_geometry(y, constants, series_terms, analytic)
        S = geometry["S"]
        P_exact = geometry["P"]
        phase_P = (imaginary * P_exact**2 / 2).exp()

        if analytic:
            center = acb(y.real.mid(), y.imag.mid())
            radius = abs(y - center).upper()
            center_geometry = local_geometry(center, constants, series_terms, False)
            exact_center = exact_endpoint_current(center_geometry["x"], constants, False)
            u = arb(ENDPOINT) * geometry["x"] / (2 * constants["mode"])
            endpoint_factor = 2 * u * (constants["K"] * u - imaginary) / (
                (u + 1) * (constants["K"] - imaginary)
            )
            exact_derivative_bound = abs(endpoint_factor * phase_P * geometry["dP_dy"]).upper()
            exact_current = centered_current_enclosure(exact_center, radius, exact_derivative_bound)

            center_P = center_geometry["P"]
            J_center = fresnel_constant * (1 + (erf_rotation * center_P).erf())
            fresnel_derivative_bound = abs(phase_P * geometry["dP_dy"]).upper()
            J_minus = centered_current_enclosure(J_center, radius, fresnel_derivative_bound)
        else:
            exact_current = exact_endpoint_current(geometry["x"], constants, False)
            J_minus = fresnel_constant * (1 + (erf_rotation * P_exact).erf())

        exact_inner = geometry["outer_factor"] * exact_current
        affine_inner = (1 + constants["mu_S"] * S) * J_minus
        affine_inner -= imaginary * constants["lambda_P"] * phase_P
        return (
            (-imaginary * S**2 / 2).exp()
            * (exact_inner - affine_inner)
            * geometry["dS_dy"]
            / (2 * pi)
        )

    return integrand


def evaluate_mode(mode_int: int, precision: int, series_terms: int, tolerance: str) -> dict[str, Any]:
    started = time.time()
    ctx.dps = precision
    ctx.threads = 1
    constants = mode_constants(mode_int)
    y_max = arb(Y_MAX_TEXT)
    require(constants["y_half"].upper() < y_max, "local box does not contain half boundary")
    integrand = exact_minus_affine_integrand(constants, series_terms)
    points = [constants["y_half"]]
    if constants["y_half"].upper() < 0:
        points.append(arb(0))
    points.append(y_max)
    pieces: list[acb] = []
    for left, right in zip(points, points[1:]):
        pieces.append(
            acb.integral(
                integrand,
                left,
                right,
                abs_tol=arb(tolerance),
                rel_tol=arb(tolerance),
                eval_limit=400_000,
                depth_limit=48,
            )
        )
    canonical = sum(pieces, acb(0))
    normalizer = (arb.pi() / (32 * arb(HEIGHT))) ** (arb(1) / 4)
    paper_rotation = (-acb(0, 1) * arb.pi() / 8).exp()
    raw_endpoint = -constants["carrier"] * constants["W0"] * canonical
    physical = 2 * normalizer * (paper_rotation * raw_endpoint).real
    require(canonical.is_finite() and physical.is_finite(), "nonfinite local nonlinear result")
    return {
        "mode": mode_int,
        "target_indicator": 1 if mode_int <= TARGET_END else 0,
        "precision_decimal_digits": precision,
        "series_terms": series_terms,
        "tolerance": tolerance,
        "y_half_ball": constants["y_half"].str(55, more=True),
        "y_upper": Y_MAX_TEXT,
        "pieces": [complex_record(piece) for piece in pieces],
        "canonical_local_exact_minus_affine_ball": complex_record(canonical),
        "raw_A_endpoint_local_exact_minus_affine_ball": complex_record(raw_endpoint),
        "physical_local_exact_minus_affine_ball": physical.str(55, more=True),
        "elapsed_seconds": round(time.time() - started, 3),
    }


def load_cache(signature: str) -> dict[int, dict[str, Any]]:
    rows: dict[int, dict[str, Any]] = {}
    if not CACHE.is_file():
        return rows
    for line in CACHE.read_text(encoding="utf-8").splitlines():
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            continue
        if record.get("signature") != signature:
            continue
        row = record.get("row")
        if isinstance(row, dict) and row.get("mode") in MODES:
            rows[int(row["mode"])] = row
    return rows


def append_cache(signature: str, row: dict[str, Any]) -> None:
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    with CACHE.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps({"signature": signature, "row": row}, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def record_to_acb(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def aggregate_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    physical_signed = arb(0)
    physical_target = arb(0)
    physical_outer = arb(0)
    physical_termwise = arb(0)
    raw_signed = acb(0)
    for row in rows:
        physical = arb(row["physical_local_exact_minus_affine_ball"])
        raw = record_to_acb(row["raw_A_endpoint_local_exact_minus_affine_ball"])
        physical_signed += physical
        physical_termwise += abs(physical)
        raw_signed += raw
        if row["target_indicator"]:
            physical_target += physical
        else:
            physical_outer += physical
    require(physical_signed.is_finite() and raw_signed.is_finite(), "nonfinite local nonlinear aggregate")
    require(physical_termwise.lower() > 0, "invalid termwise local nonlinear sum")
    return {
        "target_side_physical_local_exact_minus_affine_ball": physical_target.str(65, more=True),
        "outer_side_physical_local_exact_minus_affine_ball": physical_outer.str(65, more=True),
        "complete_physical_local_exact_minus_affine_ball": physical_signed.str(65, more=True),
        "termwise_absolute_physical_local_exact_minus_affine_sum_ball": physical_termwise.str(65, more=True),
        "signed_to_termwise_absolute_ratio_ball": (abs(physical_signed) / physical_termwise).str(55, more=True),
        "complete_raw_A_endpoint_local_exact_minus_affine_ball": complex_record(raw_signed, 65),
    }


def render_note(artifact: dict[str, Any]) -> str:
    sums = artifact["certificate"]["signed_sums"]
    return f"""# Compact exact-minus-affine transformed amplitude on A

Date: 2026-08-14

Status: rigorous compact nonlinear-amplitude value on all 84 exact A domains;
not yet a complete A endpoint block or global residual bound

The exact transformed amplitude factors as

```text
A_exact(P,S)=E(P)O(S),
E=2u(pi A m u-i)/[(u+1)(pi A m-i)],
O=v^(3/4)sigma/(v-1).                                (NA1)
```

For `x=1/[1+r(1+y)]`, the endpoint-face variable obeys

```text
P_A=sqrt(pi/2)[A sqrt(x)-2m/sqrt(x)],
I_E(P_A)=integral_(-infinity)^P_A E(P)exp(iP^2/2)dP
        =(-1)^m H_m(x)/C_E,
C_E=2sqrt(2)i(pi A m-i)/(sqrt(pi)A).                 (NA2)
```

Equation (NA2) follows from the exact identities
`C_E E(P_A)dP_A/dx=x^(-1/2)(2+i pi A^2 x)` and
`P_A^2/2=pi m^2/x+pi A^2x/4-pi A m`; `A=159577` is odd.
The previously certified inverse-Gaussian `erfc` primitive supplies `H_m`.

Thus the compact two-dimensional amplitude defect reduces exactly to

```text
R_A,m^amp=(2pi)^(-1) integral_(y_half)^({Y_MAX_TEXT}) exp(-iS(y)^2/2)
 [O(S(y))I_E(P_A(y))-I_aff(P_A(y),S(y))] S'(y)dy.    (NA3)
```

For rigorous complex-disk callbacks, each Fresnel current is evaluated at
the disk center and its variation is enclosed by the supremum of its exact
elementary ODE derivative times the disk radius.  This avoids a numerically
pathological interval-`erfc` cancellation without changing the function or
weakening the enclosure.

After restoring the A orientation, common carrier, paper rotation, and
equation-(9) normalization, deterministic summation before norms gives

```text
target-side compact nonlinear = {sums['target_side_physical_local_exact_minus_affine_ball']},
outer-side compact nonlinear  = {sums['outer_side_physical_local_exact_minus_affine_ball']},
complete compact nonlinear    = {sums['complete_physical_local_exact_minus_affine_ball']},
sum of modewise moduli        = {sums['termwise_absolute_physical_local_exact_minus_affine_sum_ball']},
signed/modulus ratio          = {sums['signed_to_termwise_absolute_ratio_ball']}. (NA4)
```

The production calculation is resumable through a modewise fsynced JSONL
cache.  The independent checker recomputes every mode at higher precision,
with more logistic-series terms and a tighter quadrature tolerance.

Pi provenance: every `pi` in (NA1)--(NA4) comes from the original equation-(9)
triangle, exact bi-Morse/Fresnel maps, odd-endpoint carrier, and paper
normalization.  No circle, polygon, fitted constant, or plotted symmetry is
inserted.

Proof boundary: exact-minus-affine transformed-amplitude value on the compact
exact A domain `y_half<=y<={Y_MAX_TEXT}` for modes 39853..39936 at `t=10^10`
only.  The affine compact value and full exact exterior are separate certified
inputs.  Their complete A assembly, the positive full-line projector jump,
`R_Dir`, complete `Q_K-T` or `T_upper`, all-height theorem, `Lambda<=0`,
PF-infinity, RH, and any prize-level conclusion remain unproved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    ctx.dps = PRODUCTION_PRECISION
    ctx.threads = 1
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["factorization"]["decision"]["exact_transformed_amplitude_factorized"] is True, "factorization dependency drift")
    require(dependencies["local_affine"]["decision"]["all_84_local_exact_affine_modes_certified"] is True, "local affine dependency drift")
    require(dependencies["orientation"]["decision"]["A_endpoint_sign_identified_as_epsilon_A_minus_one"] is True, "orientation dependency drift")
    require(dependencies["erfc_reduction"]["decision"]["exact_endpoint_driver_erfc_primitive_proved"] is True, "erfc dependency drift")
    require(dependencies["exact_exterior"]["decision"]["full_exact_exterior_value_after_endpoint_replacement_certified"] is True, "exterior dependency drift")
    symbolic = symbolic_normalization_certificate()

    signature = config_signature(PRODUCTION_PRECISION, PRODUCTION_SERIES_TERMS, PRODUCTION_TOLERANCE)
    cached = load_cache(signature)
    computed_now = 0
    for index, mode in enumerate(MODES, start=1):
        if mode not in cached:
            row = evaluate_mode(mode, PRODUCTION_PRECISION, PRODUCTION_SERIES_TERMS, PRODUCTION_TOLERANCE)
            append_cache(signature, row)
            cached[mode] = row
            computed_now += 1
        if index % 8 == 0:
            print(f"local exact nonlinear rows {index}/{len(MODES)}", flush=True)
    rows = [cached[mode] for mode in MODES]
    signed_sums = aggregate_rows(rows)
    artifact = {
        "kind": STEM,
        "status": "compact_exact_minus_affine_transformed_amplitude_certified_for_all_84_A_modes_complete_A_assembly_open",
        "passed": True,
        "symbolic_certificate": symbolic,
        "certificate": {
            "height": HEIGHT,
            "endpoint": ENDPOINT,
            "mode_range": [LOWER_START, UPPER_END],
            "mode_count": len(rows),
            "target_side_count": sum(row["target_indicator"] for row in rows),
            "outer_side_count": sum(1 - row["target_indicator"] for row in rows),
            "y_upper": Y_MAX_TEXT,
            "precision_decimal_digits": PRODUCTION_PRECISION,
            "series_terms": PRODUCTION_SERIES_TERMS,
            "tolerance": PRODUCTION_TOLERANCE,
            "analytic_current_enclosure": "center value plus disk radius times exact ODE-derivative supremum",
            "rows": rows,
            "signed_sums": signed_sums,
        },
        "decision": {
            "exact_endpoint_current_normalization_proved": True,
            "raw_interval_erfc_disk_cancellation_avoided_by_exact_ODE_enclosure": True,
            "compact_exact_minus_affine_transformed_amplitude_certified": True,
            "all_84_local_nonlinear_modes_certified": True,
            "signed_local_nonlinear_sum_certified_before_norms": True,
            "exact_exterior_recomputed_here": False,
            "complete_A_endpoint_block_proved": False,
            "R_Dir_bound_proved": False,
            "rh_implication": False,
        },
        "cache": {
            "path": relative(CACHE),
            "sha256": file_hash(CACHE),
            "signature": signature,
            "row_formula_version": ROW_FORMULA_VERSION,
            "rows_reused": len(rows) - computed_now,
            "rows_computed_now": computed_now,
            "fsynced_after_each_row": True,
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
        "next_obligation": "Assemble the certified compact affine value, compact exact-minus-affine value, and full exact exterior with the positive full-line projector jump and orientation factors; independently check the complete A endpoint identity before using it in R_Dir.",
        "proof_boundary": "Exact-minus-affine transformed-amplitude value on the compact exact A domain y_half<=y<=0.0037 for modes 39853..39936 at t=10^10 only. No complete A endpoint block, R_Dir estimate, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified compact exact-minus-affine transformed amplitude on A", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
