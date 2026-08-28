#!/usr/bin/env python3
"""Certify the exact affine curved-face strip on a finite A Airy box."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_A_local_curved_face_strip_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
CACHE = REPO_ROOT / f"work/rh_compute/results/{STEM}_rows.jsonl"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "roster": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_A_transition_roster_gate.json",
    "affine_remainder": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_remainder_identity_gate.json",
    "orientation": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_corner_projector_orientation_gate.json",
}

HEIGHT = 10_000_000_000
A = 159_577
LOWER_START = 39_853
TARGET_END = 39_894
UPPER_END = 39_936
MODES = tuple(range(LOWER_START, UPPER_END + 1))
Y_MAX_TEXT = "0.0037"
SERIES_TERMS = 24
PRECISION = 60
ABS_TOL_TEXT = "1e-14"
CACHE_VERSION = 1


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


def config_signature(
    precision: int,
    series_terms: int,
    abs_tol_text: str,
) -> str:
    payload = {
        "cache_version": CACHE_VERSION,
        "height": HEIGHT,
        "endpoint": A,
        "mode_range": [LOWER_START, UPPER_END],
        "y_max": Y_MAX_TEXT,
        "precision": precision,
        "series_terms": series_terms,
        "abs_tol": abs_tol_text,
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("ascii")).hexdigest()


def q_and_derivative(y: acb, series_terms: int) -> tuple[acb, acb, arb, arb]:
    """Evaluate q(y)=2[y-log(1+y)]/y^2 by a certified power series."""
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
    q_tail = arb(2) * radius**series_terms / (
        arb(series_terms + 2) * (1 - radius)
    )
    q_prime_tail = arb(2) * radius ** (series_terms - 1) / (1 - radius)
    q += acb(arb(0, q_tail), arb(0, q_tail))
    q_prime += acb(arb(0, q_prime_tail), arb(0, q_prime_tail))
    return q, q_prime, q_tail, q_prime_tail


def mode_constants(mode_int: int) -> dict[str, Any]:
    pi, t, endpoint = arb.pi(), arb(HEIGHT), arb(A)
    mode = arb(mode_int)
    r = t / (2 * pi * mode**2)
    denominator = 2 * pi * mode**2 + t
    alpha = 2 * mode + t / (pi * mode)
    a = (pi * mode / alpha).sqrt() * (endpoint - alpha)
    rho = -(2 * t).sqrt() * (pi * endpoint * mode + denominator) / (
        2 * denominator ** arb("1.5")
    )
    K = pi * endpoint * mode
    lambda_P = acb(3 * K, -1) / (acb(2 * K, -2) * K.sqrt())
    mu_S = arb(5) / 12 * (arb(2) / t).sqrt()
    y_half = 2 * pi * mode**2 / t - 1
    W0 = acb(K, -1) * (4 * r ** arb("0.75") / (endpoint * (pi * t).sqrt()))
    phase = t * ((t / (2 * pi)).log() - 1) / 2 - t * mode.log()
    return {
        "mode": mode,
        "r": r,
        "a": a,
        "rho": rho,
        "lambda_P": lambda_P,
        "mu_S": mu_S,
        "y_half": y_half,
        "W0": W0,
        "carrier": (acb(0, 1) * phase).exp(),
    }


def face_integrand(constants: dict[str, Any], series_terms: int):
    pi, t, endpoint = arb.pi(), arb(HEIGHT), arb(A)
    imaginary = acb(0, 1)
    root_two = arb(2).sqrt()
    J0_scale = (pi / 2).sqrt() * (imaginary * pi / 4).exp()
    J0_rotation = (-imaginary * pi / 4).exp() / root_two

    def J0(P: acb) -> acb:
        return J0_scale * (J0_rotation * P).erf()

    def integrand(y: acb, analytic: bool) -> acb:
        # Arb deliberately probes complex disks wider than the integration
        # segment.  The power-series tail certificate is only asserted on
        # |y|<0.01, so reject wider analyticity probes with a non-finite ball;
        # the integrator will then subdivide instead of accepting an
        # uncertified extrapolation.
        if abs(y).upper() >= arb("0.01"):
            if analytic:
                return acb("nan")
            raise RuntimeError("real quadrature point escaped the certified disk")
        q, q_prime, _, _ = q_and_derivative(y, series_terms)
        sqrt_q = q.sqrt(analytic=analytic)
        S = (t / 2).sqrt() * y * sqrt_q
        dS_dy = (t / 2).sqrt() * (
            sqrt_q + y * q_prime / (2 * sqrt_q)
        )
        g = 1 + constants["r"] * (1 + y)
        sqrt_g = g.sqrt(analytic=analytic)
        P_exact = (pi / 2).sqrt() * (
            endpoint / sqrt_g - 2 * constants["mode"] * sqrt_g
        )
        P_linear = constants["a"] + constants["rho"] * S
        affine_zero = 1 + constants["mu_S"] * S
        inner = affine_zero * (J0(P_exact) - J0(P_linear))
        inner -= imaginary * constants["lambda_P"] * (
            (imaginary * P_exact**2 / 2).exp()
            - (imaginary * P_linear**2 / 2).exp()
        )
        return (
            (-imaginary * S**2 / 2).exp()
            * inner
            * dS_dy
            / (2 * pi)
        )

    return integrand


def evaluate_mode(
    mode_int: int,
    precision: int,
    series_terms: int,
    abs_tol_text: str,
) -> dict[str, Any]:
    ctx.dps = precision
    ctx.threads = 1
    pi, t = arb.pi(), arb(HEIGHT)
    imaginary = acb(0, 1)
    constants = mode_constants(mode_int)
    y_max = arb(Y_MAX_TEXT)
    require(constants["y_half"].upper() < y_max, "Airy box does not contain half boundary")
    integrand = face_integrand(constants, series_terms)
    points = [constants["y_half"]]
    if constants["y_half"].upper() < 0:
        points.append(arb(0))
    points.append(y_max)
    local_strip = acb(0)
    for left, right in zip(points, points[1:]):
        local_strip += acb.integral(
            integrand,
            left,
            right,
            abs_tol=arb(abs_tol_text),
            rel_tol=arb(abs_tol_text),
            eval_limit=300_000,
            depth_limit=48,
        )

    normalizer = (pi / (32 * t)) ** (arb(1) / 4)
    paper_rotation = (-imaginary * pi / 8).exp()
    raw_endpoint = -constants["carrier"] * constants["W0"] * local_strip
    physical = 2 * normalizer * (paper_rotation * raw_endpoint).real
    q_at_radius = q_and_derivative(acb(y_max), series_terms)
    return {
        "mode": mode_int,
        "target_indicator": 1 if mode_int <= TARGET_END else 0,
        "y_half_ball": constants["y_half"].str(50, more=True),
        "y_upper": Y_MAX_TEXT,
        "local_strip_canonical_ball": complex_record(local_strip),
        "local_strip_raw_A_endpoint_ball": complex_record(raw_endpoint),
        "physical_local_curved_face_strip_ball": physical.str(55, more=True),
        "q_series_tail_at_radius_bound": q_at_radius[2].str(40, more=True),
        "q_prime_series_tail_at_radius_bound": q_at_radius[3].str(40, more=True),
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
    record = {"signature": signature, "row": row}
    with CACHE.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def aggregate_rows(rows: list[dict[str, Any]]) -> dict[str, str]:
    signed = arb(0)
    target_side = arb(0)
    outer_side = arb(0)
    termwise = arb(0)
    for row in rows:
        value = arb(row["physical_local_curved_face_strip_ball"])
        signed += value
        termwise += abs(value)
        if row["target_indicator"]:
            target_side += value
        else:
            outer_side += value
    ratio = abs(signed) / termwise
    require(signed < arb("-0.009"), "local face-strip signed scale drift")
    require(signed > arb("-0.012"), "local face-strip signed scale exceeded audit interval")
    require(termwise > arb("0.03"), "local face-strip termwise scale drift")
    require(ratio < arb("0.35"), "local face-strip cancellation ratio drift")
    return {
        "target_side_physical_local_curved_face_strip_ball": target_side.str(65, more=True),
        "outer_side_physical_local_curved_face_strip_ball": outer_side.str(65, more=True),
        "complete_physical_local_curved_face_strip_ball": signed.str(65, more=True),
        "termwise_absolute_physical_local_curved_face_strip_sum_ball": termwise.str(65, more=True),
        "signed_to_termwise_absolute_ratio_ball": ratio.str(55, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    s = c["signed_sums"]
    return f"""# Exact local curved-face strip on the A Airy box

Date: 2026-08-13

Status: rigorous finite-box affine curved-face correction on all 84 A modes;
not a complete curved-face, amplitude, A-endpoint, or global residual bound

Put `y=v-1`.  The logistic Morse coordinate is analytic at `y=0` in the
form

```text
sigma=y sqrt(q(y)),
q(y)=2[y-log(1+y)]/y^2
    =2 sum_(k>=0)(-1)^k y^k/(k+2),
S=sqrt(t/2)sigma.                                     (CF1)
```

The evaluator uses `{SERIES_TERMS}` terms.  For `|y|<=r<1`, the omitted
tails are bounded directly by

```text
|q-q_N| <=2r^N/[(N+2)(1-r)],
|q'-q_N'|<=2r^(N-1)/(1-r).                            (CF2)
```

At `r={Y_MAX_TEXT}` their saved bounds are respectively

```text
{c['q_series_tail_at_radius_bound']},
{c['q_prime_series_tail_at_radius_bound']}.            (CF3)
```

This removes the apparent `0/0` at the null point without a floating-point
branch choice.  In the same variable the exact A face is

```text
P_A(y)=sqrt(pi/2){{A/[1+r_m(1+y)]^(1/2)
                   -2m[1+r_m(1+y)]^(1/2)}}.           (CF4)
```

For `P_lin=a+rho S`, the exact affine strip integral on
`y_half<=y<={Y_MAX_TEXT}` is evaluated through the analytic primitives

```text
J_0(P)=integral_0^P exp(iu^2/2)du,
J_1(P)=-i exp(iP^2/2),                                (CF5)

R_face,loc=(2pi)^(-1) integral e^(-iS^2/2)
 {{(1+mu S)[J_0(P_A)-J_0(P_lin)]
   +lambda[J_1(P_A)-J_1(P_lin)]}}dS.                  (CF6)
```

All 84 half-boundaries lie inside this box; its upper edge corresponds to
`S` a little above `261`.  Restoring the exact A sign, common carrier,
paper rotation, and normalization gives

```text
target-side local strip={s['target_side_physical_local_curved_face_strip_ball']},
outer-side local strip ={s['outer_side_physical_local_curved_face_strip_ball']},
complete signed strip  ={s['complete_physical_local_curved_face_strip_ball']},
sum of modewise moduli ={s['termwise_absolute_physical_local_curved_face_strip_sum_ball']},
signed/modulus ratio   ={s['signed_to_termwise_absolute_ratio_ball']}. (CF7)
```

The exact local face correction is therefore not a perturbation at the
`1.4e-4` target scale.  It has to be retained as part of the A Airy carrier,
not charged as an error to the affine tangent wedge.  Equation (CF7) does
not include `y>{Y_MAX_TEXT}`, the exact-minus-affine transformed amplitude,
or any of the global projector companions, so it is not an `R_Dir` estimate.

The production run is resumable mode by mode in an fsynced JSONL cache.  An
independent higher-precision replay recomputes every mode with more series
terms and tighter quadrature tolerance.

Pi provenance: every `pi` in (CF1)--(CF7) comes from equation (9), the exact
bi-Morse transformation, Fresnel primitives, odd-endpoint phase carrier,
and paper normalization.  No fitted constant is introduced.

Proof boundary: the exact affine curved-face strip on the finite real box
`y_half<=y<={Y_MAX_TEXT}` for modes 39853..39936 only.  No exterior face
strip, transformed-amplitude remainder, complete A endpoint theorem,
`R_Dir` estimate, complete `Q_K-T` or `T_upper`, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    # Cached rows must be parsed and aggregated at the same precision as
    # freshly evaluated rows.  Without this initialization, a fully resumed
    # run would silently use python-flint's process default precision.
    ctx.dps = PRECISION
    ctx.threads = 1
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["roster"]["decision"]["complete_84_mode_A_transition_roster_evaluated"] is True, "roster dependency drift")
    require(dependencies["affine_remainder"]["decision"]["exact_defect_split_into_signed_amplitude_and_face_remainders"] is True, "face-remainder dependency drift")
    require(dependencies["orientation"]["decision"]["A_endpoint_sign_identified_as_epsilon_A_minus_one"] is True, "orientation dependency drift")

    signature = config_signature(PRECISION, SERIES_TERMS, ABS_TOL_TEXT)
    cached = load_cache(signature)
    computed_now = 0
    for index, mode in enumerate(MODES, start=1):
        if mode not in cached:
            row = evaluate_mode(mode, PRECISION, SERIES_TERMS, ABS_TOL_TEXT)
            append_cache(signature, row)
            cached[mode] = row
            computed_now += 1
        if index % 8 == 0:
            print(f"local curved-face rows {index}/{len(MODES)}", flush=True)
    rows = [cached[mode] for mode in MODES]
    signed_sums = aggregate_rows(rows)
    q_radius = q_and_derivative(acb(arb(Y_MAX_TEXT)), SERIES_TERMS)
    certificate = {
        "height": HEIGHT,
        "endpoint": A,
        "mode_range": [LOWER_START, UPPER_END],
        "mode_count": len(rows),
        "target_side_count": sum(row["target_indicator"] for row in rows),
        "outer_side_count": sum(1 - row["target_indicator"] for row in rows),
        "y_upper": Y_MAX_TEXT,
        "series_terms": SERIES_TERMS,
        "q_series_tail_at_radius_bound": q_radius[2].str(45, more=True),
        "q_prime_series_tail_at_radius_bound": q_radius[3].str(45, more=True),
        "rows": rows,
        "signed_sums": signed_sums,
    }
    artifact = {
        "kind": STEM,
        "status": "exact_affine_curved_face_strip_on_finite_A_Airy_box_certified_for_all_84_modes_exterior_and_amplitude_open",
        "passed": True,
        "certificate": certificate,
        "decision": {
            "logistic_null_coordinate_regularized_by_certified_power_series": True,
            "all_84_half_boundaries_contained_in_finite_Airy_box": True,
            "exact_curved_face_and_tangent_used_without_face_Taylor_truncation": True,
            "affine_inner_strip_reduced_to_exact_Fresnel_primitives": True,
            "signed_local_curved_face_strip_summed_before_norms": True,
            "local_curved_face_strip_is_small_at_R_Dir_target_scale": False,
            "exterior_curved_face_strip_bound_proved": False,
            "transformed_amplitude_remainder_bound_proved": False,
            "complete_A_endpoint_block_proved": False,
            "R_Dir_bound_proved": False,
            "rh_implication": False,
        },
        "cache": {
            "path": relative(CACHE),
            "sha256": file_hash(CACHE),
            "signature": signature,
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
            "precision_decimal_digits": PRECISION,
            "absolute_and_relative_tolerance": ABS_TOL_TEXT,
        },
        "next_obligation": "Promote the local curved-face strip into the A Airy carrier rather than treating it as error. Bound the y>0.0037 face-strip complement by phase-adapted integration by parts and certify the exact-minus-affine transformed-amplitude integral on the same finite box, preserving the 84-mode signed carrier.",
        "proof_boundary": "Exact affine curved-face strip on y_half<=y<=0.0037 for modes 39853..39936 only. No exterior face strip, transformed-amplitude remainder, complete A endpoint theorem, R_Dir estimate, complete Q_K-T or T_upper, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified exact local affine A curved-face strip", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
