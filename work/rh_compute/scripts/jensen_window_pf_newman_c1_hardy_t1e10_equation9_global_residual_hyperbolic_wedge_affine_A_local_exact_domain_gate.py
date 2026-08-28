#!/usr/bin/env python3
"""Certify the localized exact-domain affine A carrier on the finite box."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_A_local_exact_domain_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
CACHE = REPO_ROOT / f"work/rh_compute/results/{STEM}_rows.jsonl"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "regrouping": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_A_localized_exact_domain_regrouping_gate.json",
    "local_face": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_A_local_curved_face_strip_gate.json",
    "orientation": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_corner_projector_orientation_gate.json",
}

HEIGHT = 10_000_000_000
A = 159_577
LOWER_START = 39_853
TARGET_END = 39_894
UPPER_END = 39_936
MODES = tuple(range(LOWER_START, UPPER_END + 1))
Y0_TEXT = "0.0037"
SERIES_TERMS = 24
PRECISION = 60
ABS_TOL_TEXT = "1e-13"
CACHE_VERSION = 1
ROW_FORMULA_VERSION = "exact-affine-domain-Jminus-v1"


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


def config_signature(precision: int, series_terms: int, abs_tol_text: str) -> str:
    payload = {
        "cache_version": CACHE_VERSION,
        "row_formula_version": ROW_FORMULA_VERSION,
        "height": HEIGHT,
        "endpoint": A,
        "mode_range": [LOWER_START, UPPER_END],
        "y_max": Y0_TEXT,
        "precision": precision,
        "series_terms": series_terms,
        "abs_tol": abs_tol_text,
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("ascii")).hexdigest()


def q_and_derivative(y: acb, series_terms: int) -> tuple[acb, acb, arb, arb]:
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
    return q, q_prime, q_tail, q_prime_tail


def mode_constants(mode_int: int) -> dict[str, Any]:
    pi, t, endpoint = arb.pi(), arb(HEIGHT), arb(A)
    mode = arb(mode_int)
    r = t / (2 * pi * mode**2)
    denominator = 2 * pi * mode**2 + t
    K = pi * endpoint * mode
    lambda_P = acb(3 * K, -1) / (acb(2 * K, -2) * K.sqrt())
    mu_S = arb(5) / 12 * (arb(2) / t).sqrt()
    y_half = 2 * pi * mode**2 / t - 1
    W0 = acb(K, -1) * (4 * r ** arb("0.75") / (endpoint * (pi * t).sqrt()))
    phase = t * ((t / (2 * pi)).log() - 1) / 2 - t * mode.log()
    return {
        "mode": mode,
        "r": r,
        "lambda_P": lambda_P,
        "mu_S": mu_S,
        "y_half": y_half,
        "W0": W0,
        "carrier": (acb(0, 1) * phase).exp(),
    }


def exact_local_integrand(constants: dict[str, Any], series_terms: int):
    pi, t, endpoint = arb.pi(), arb(HEIGHT), arb(A)
    imaginary = acb(0, 1)
    root_two = arb(2).sqrt()
    fresnel_constant = (pi / 2).sqrt() * (imaginary * pi / 4).exp()
    erf_rotation = (-imaginary * pi / 4).exp() / root_two

    def integrand(y: acb, analytic: bool) -> acb:
        if abs(y).upper() >= arb("0.01"):
            if analytic:
                return acb("nan")
            raise RuntimeError("real quadrature point escaped the certified disk")
        q, q_prime, _, _ = q_and_derivative(y, series_terms)
        sqrt_q = q.sqrt(analytic=analytic)
        S = (t / 2).sqrt() * y * sqrt_q
        dS_dy = (t / 2).sqrt() * (sqrt_q + y * q_prime / (2 * sqrt_q))
        g = 1 + constants["r"] * (1 + y)
        sqrt_g = g.sqrt(analytic=analytic)
        P_exact = (pi / 2).sqrt() * (endpoint / sqrt_g - 2 * constants["mode"] * sqrt_g)
        J0 = fresnel_constant * (erf_rotation * P_exact).erf()
        J_minus = fresnel_constant + J0
        inner = (1 + constants["mu_S"] * S) * J_minus
        inner -= imaginary * constants["lambda_P"] * (imaginary * P_exact**2 / 2).exp()
        return (-imaginary * S**2 / 2).exp() * inner * dS_dy / (2 * pi)

    return integrand


def evaluate_mode(mode_int: int, precision: int, series_terms: int, abs_tol_text: str) -> dict[str, Any]:
    ctx.dps = precision
    ctx.threads = 1
    pi, t = arb.pi(), arb(HEIGHT)
    imaginary = acb(0, 1)
    constants = mode_constants(mode_int)
    y0 = arb(Y0_TEXT)
    require(constants["y_half"].upper() < y0, "local box does not contain half boundary")
    integrand = exact_local_integrand(constants, series_terms)
    points = [constants["y_half"]]
    if constants["y_half"].upper() < 0:
        points.append(arb(0))
    points.append(y0)
    local_exact = acb(0)
    for left, right in zip(points, points[1:]):
        local_exact += acb.integral(
            integrand,
            left,
            right,
            abs_tol=arb(abs_tol_text),
            rel_tol=arb(abs_tol_text),
            eval_limit=400_000,
            depth_limit=50,
        )
    normalizer = (pi / (32 * t)) ** (arb(1) / 4)
    paper_rotation = (-imaginary * pi / 8).exp()
    raw_endpoint = -constants["carrier"] * constants["W0"] * local_exact
    physical = 2 * normalizer * (paper_rotation * raw_endpoint).real
    return {
        "mode": mode_int,
        "target_indicator": 1 if mode_int <= TARGET_END else 0,
        "y_half_ball": constants["y_half"].str(50, more=True),
        "y_upper": Y0_TEXT,
        "canonical_local_exact_affine_ball": complex_record(local_exact),
        "raw_A_endpoint_local_exact_affine_ball": complex_record(raw_endpoint),
        "physical_local_exact_affine_ball": physical.str(55, more=True),
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


def aggregate_rows(rows: list[dict[str, Any]]) -> dict[str, str]:
    signed = arb(0)
    target = arb(0)
    outer = arb(0)
    termwise = arb(0)
    for row in rows:
        value = arb(row["physical_local_exact_affine_ball"])
        signed += value
        termwise += abs(value)
        if row["target_indicator"]:
            target += value
        else:
            outer += value
    require(signed.is_finite() and target.is_finite() and outer.is_finite(), "nonfinite signed local exact sum")
    require(termwise.is_finite() and termwise.lower() > 0, "invalid local exact termwise sum")
    ratio = abs(signed) / termwise
    return {
        "target_side_physical_local_exact_affine_ball": target.str(65, more=True),
        "outer_side_physical_local_exact_affine_ball": outer.str(65, more=True),
        "complete_physical_local_exact_affine_ball": signed.str(65, more=True),
        "termwise_absolute_physical_local_exact_affine_sum_ball": termwise.str(65, more=True),
        "signed_to_termwise_absolute_ratio_ball": ratio.str(55, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    s = c["signed_sums"]
    return f"""# Local exact-domain affine A carrier

Date: 2026-08-13

Status: rigorous finite-box exact-domain affine carrier on all 84 A modes;
not an exterior, nonlinear-amplitude, complete A, or global residual bound

For `P_A(S)` the exact A face, define the lower Fresnel tail

```text
T_-(P)=integral_(-infinity)^P exp(iu^2/2)du
 =sqrt(pi/2)exp(i*pi/4)
  [1+erf(exp(-i*pi/4)P/sqrt(2))].                    (LE1)
```

The affine amplitude `1+lambda_P P+mu_S S` can then be integrated exactly
in `P` over the exact domain.  On `y_half<=y<={Y0_TEXT}` this gives

```text
C_A,loc^exact-aff=(2pi)^(-1) integral exp(-iS^2/2)
 {{(1+mu_S S)T_-(P_A(S))-i lambda_P exp(iP_A(S)^2/2)}}dS. (LE2)
```

This is the localized object required by Section 11.434.  It is equal to
local tangent plus local face, and contains no tangent exterior or artificial
tangent stationary point.

Restoring `epsilon_A=-1`, the common odd-endpoint carrier, paper rotation,
and equation-(9) normalization, deterministic summation before norms gives

```text
target-side local exact affine={s['target_side_physical_local_exact_affine_ball']},
outer-side local exact affine ={s['outer_side_physical_local_exact_affine_ball']},
complete local exact affine   ={s['complete_physical_local_exact_affine_ball']},
sum of modewise moduli        ={s['termwise_absolute_physical_local_exact_affine_sum_ball']},
signed/modulus ratio          ={s['signed_to_termwise_absolute_ratio_ball']}. (LE3)
```

The production calculation is resumable through a modewise fsynced JSONL
cache.  An independent replay uses higher precision, more logistic-series
terms, and tighter quadrature tolerance.

Equation (LE3) is only the affine-amplitude contribution on the finite exact
domain.  It excludes the exact exterior `y>{Y0_TEXT}`, the exact-minus-affine
transformed amplitude even inside the box, the positive full-line projector
jump, and all other global sectors.  It is therefore not an `R_Dir` estimate.

Pi provenance: every `pi` in (LE1)--(LE3) comes from the exact Fresnel
primitive, equation-(9) bi-Morse phase, odd-endpoint carrier, and paper
normalization.  No fitted constant is introduced.

Proof boundary: exact affine-amplitude integral on the finite exact A domain
`y_half<=y<={Y0_TEXT}` for modes 39853..39936 only.  No exact exterior affine
current, exact-minus-affine transformed-amplitude bound, complete A endpoint
theorem, `R_Dir` estimate, complete `Q_K-T` or `T_upper`, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    ctx.dps = PRECISION
    ctx.threads = 1
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["regrouping"]["decision"]["tangent_exterior_must_be_removed_before_estimation"] is True, "regrouping dependency drift")
    require(dependencies["local_face"]["decision"]["exact_curved_face_and_tangent_used_without_face_Taylor_truncation"] is True, "local-face dependency drift")
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
            print(f"local exact affine rows {index}/{len(MODES)}", flush=True)
    rows = [cached[mode] for mode in MODES]
    signed_sums = aggregate_rows(rows)
    artifact = {
        "kind": STEM,
        "status": "finite_box_exact_domain_affine_A_carrier_certified_for_all_84_modes_exterior_and_nonlinear_amplitude_open",
        "passed": True,
        "certificate": {
            "height": HEIGHT,
            "endpoint": A,
            "mode_range": [LOWER_START, UPPER_END],
            "mode_count": len(rows),
            "target_side_count": sum(row["target_indicator"] for row in rows),
            "outer_side_count": sum(1 - row["target_indicator"] for row in rows),
            "y_upper": Y0_TEXT,
            "series_terms": SERIES_TERMS,
            "rows": rows,
            "signed_sums": signed_sums,
        },
        "decision": {
            "lower_Fresnel_tail_implemented_by_exact_erf_primitive": True,
            "affine_P_moment_integrated_exactly": True,
            "local_exact_domain_used_instead_of_global_tangent_plus_local_face": True,
            "all_84_local_exact_affine_modes_certified": True,
            "signed_local_exact_affine_sum_certified_before_norms": True,
            "tangent_exterior_included_in_local_carrier": False,
            "exact_exterior_affine_current_bound_proved": False,
            "transformed_amplitude_remainder_bound_proved": False,
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
            "precision_decimal_digits": PRECISION,
            "absolute_and_relative_tolerance": ABS_TOL_TEXT,
        },
        "next_obligation": "Evaluate the exact exterior affine current using the inverse-P Fresnel expansion and common exact-face x phase, with an incomplete stationary transition for modes 39927..39936. Separately certify the exact-minus-affine transformed amplitude on the local box before assembling the complete A carrier.",
        "proof_boundary": "Exact affine-amplitude integral on the finite exact A domain y_half<=y<=0.0037 for modes 39853..39936 only. No exact exterior affine current, exact-minus-affine transformed-amplitude bound, complete A endpoint theorem, R_Dir estimate, complete Q_K-T or T_upper, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified finite-box local exact-domain affine A carrier", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
