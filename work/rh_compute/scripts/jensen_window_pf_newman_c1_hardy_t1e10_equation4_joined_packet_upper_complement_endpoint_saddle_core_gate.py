#!/usr/bin/env python3
"""Certify the grouped 230-label upper-complement endpoint-saddle core."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[3]
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_endpoint_saddle_core_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
COMPLEMENT = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_ordinary_complementary_half_lattice_gate.json"

PRECISION_BITS = 384
PANELS = 32
TOLERANCE = arb("1e-45")
HEIGHT = arb("10000000000")
LOWER = arb("621.5")
SPLIT = arb("621.5625")
Q_PLUS = arb("2561211.5")
LABEL_COUNT = 230


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8", newline="\n")
    os.replace(temporary, path)


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            process.cpu_affinity([process.cpu_affinity()[0]])
            return "below_normal_one_cpu"
        process.nice(10)
        process.cpu_affinity([process.cpu_affinity()[0]])
        return "nice_10_one_cpu"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def arb_record(value: arb, digits: int = 80) -> dict[str, Any]:
    return {
        "ball": value.str(digits, more=True),
        "lower_float": float(value.lower()),
        "mid_float": float(value.mid()),
        "upper_float": float(value.upper()),
    }


def acb_record(value: acb, digits: int = 80) -> dict[str, str]:
    return {
        "real_ball": value.real.str(digits, more=True),
        "imag_ball": value.imag.str(digits, more=True),
        "absolute_ball": abs(value).str(digits, more=True),
    }


def lower_saddle(q: arb, p: arb) -> arb:
    discriminant = (q * q - 4 * p).sqrt()
    return 2 * p / (q + discriminant)


def grouped_integrand_factory(
    *, t: arb, lower: arb, q0: arb, count: int
) -> Callable[[acb, bool], acb]:
    pi = arb.pi()
    imaginary = acb(0, 1)
    endpoint_phase = -t * lower.log() - pi * lower * lower + 2 * pi * q0 * lower
    endpoint_carrier = lower ** arb("-0.5") * (imaginary * endpoint_phase).exp()
    linear = 2 * pi * (q0 - lower) - t / lower

    def integrand(x: acb, analytic: bool) -> acb:
        w = x / lower
        root = (1 + w).sqrt(analytic=analytic)
        log_remainder = (1 + w).log() - w
        relative_phase = linear * x - t * log_remainder - pi * x * x
        step = (2 * pi * imaginary * x).exp()
        geometric = (1 - step**count) / (1 + step)
        return endpoint_carrier / root * (imaginary * relative_phase).exp() * geometric

    return integrand


def integrate_panels(
    integrand: Callable[[acb, bool], acb], width: arb, panels: int
) -> tuple[acb, list[dict[str, str]]]:
    total = acb(0)
    rows: list[dict[str, str]] = []
    for index in range(panels):
        left = width * index / panels
        right = width * (index + 1) / panels
        value = acb.integral(
            integrand,
            left,
            right,
            abs_tol=TOLERANCE,
            rel_tol=TOLERANCE,
            eval_limit=300_000,
            depth_limit=40,
        )
        require(value.is_finite(), f"nonfinite panel {index}")
        total += value
        rows.append(
            {
                "index": str(index),
                "left": left.str(35, more=True),
                "right": right.str(35, more=True),
                "absolute_ball": abs(value).str(35, more=True),
            }
        )
    return total, rows


def build_certificate() -> dict[str, Any]:
    pi = arb.pi()
    p = HEIGHT / (2 * pi)
    q_last = Q_PLUS + LABEL_COUNT - 1
    root_first = lower_saddle(Q_PLUS, p)
    root_last = lower_saddle(q_last, p)
    q_s_lower = LOWER + p / LOWER
    q_s_split = SPLIT + p / SPLIT

    require(LABEL_COUNT % 2 == 0, "endpoint cancellation needs an even roster")
    require(root_last > LOWER, "last stationary root escaped below the endpoint")
    require(root_first < SPLIT, "first stationary root escaped above the rational split")
    require(Q_PLUS - q_s_split > 0, "stationary phase remains beyond the rational split")
    require(q_s_lower - q_last > 0, "last roster label is not stationary in the core")

    width = SPLIT - LOWER
    integrand = grouped_integrand_factory(
        t=HEIGHT,
        lower=LOWER,
        q0=Q_PLUS,
        count=LABEL_COUNT,
    )
    integral, panels = integrate_panels(integrand, width, PANELS)
    require(integral.is_finite(), "grouped endpoint-saddle integral is not finite")

    return {
        "passed": True,
        "precision_bits": PRECISION_BITS,
        "panels": PANELS,
        "absolute_tolerance": TOLERANCE.str(20, more=True),
        "t": HEIGHT.str(50, more=True),
        "p_equals_t_over_2pi_ball": arb_record(p),
        "lower_endpoint": LOWER.str(30, more=True),
        "rational_split": SPLIT.str(30, more=True),
        "split_width": width.str(30, more=True),
        "q_range": [Q_PLUS.str(30, more=True), q_last.str(30, more=True)],
        "label_count": LABEL_COUNT,
        "endpoint_dirichlet_value_exact": "sum_(j=0)^229 (-1)^j=0",
        "q_s_at_lower_ball": arb_record(q_s_lower),
        "q_s_at_split_ball": arb_record(q_s_split),
        "last_label_stationary_margin_at_lower_ball": arb_record(q_s_lower - q_last),
        "first_label_nonstationary_margin_at_split_ball": arb_record(Q_PLUS - q_s_split),
        "first_label_saddle_root_ball": arb_record(root_first),
        "last_label_saddle_root_ball": arb_record(root_last),
        "first_root_to_split_margin_ball": arb_record(SPLIT - root_first),
        "last_root_from_lower_margin_ball": arb_record(root_last - LOWER),
        "grouped_endpoint_saddle_core_ball": acb_record(integral),
        "panel_rows": panels,
        "exact_decomposition": (
            "T_upper=P_core+R_upper, where P_core=sum_(q=q_+)^(q_++229) "
            "integral_L^B y^(-s)e^(-i*pi*y^2+2*pi*i*q*y)dy, B=L+1/16, "
            "and every phase in R_upper is nonstationary on its retained interval"
        ),
    }


def note_text(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Grouped upper-complement endpoint-saddle core

Date: 2026-08-27

Status: the complete 230-label stationary core is rigorously enclosed on a
rational endpoint interval; its nonstationary continuation remains open.

The complementary half-lattice gate leaves exactly the upper labels

```text
q=q_+,q_++1,...,q_++229,       q_+=2561211.5.
```

All of their lower saddle roots lie in the rational interval

```text
L=621.5 <= y <= B=621.5625=L+1/16.
```

The extreme roots are `{c['last_label_saddle_root_ball']['ball']}` and
`{c['first_label_saddle_root_ball']['ball']}`.  The first-label root remains
below `B` by `{c['first_root_to_split_margin_ball']['ball']}`, while at `B`
the nearest phase is already nonstationary by
`{c['first_label_nonstationary_margin_at_split_ball']['ball']}`.

Because the roster length is even and `L` is a half integer, its grouped
Dirichlet amplitude has the exact endpoint zero

```text
sum_(j=0)^229 exp(2*pi*i*j*L)=sum_(j=0)^229 (-1)^j=0.
```

With `x=y-L`, the builder evaluates the complete finite block through the
single stable geometric amplitude

```text
(1-exp(2*pi*i*230*x))/(1+exp(2*pi*i*x))
```

and the endpoint-centered logarithmic phase.  Thirty-two Arb panels at 384
bits give

```text
P_core = {c['grouped_endpoint_saddle_core_ball']['real_ball']}
       + i {c['grouped_endpoint_saddle_core_ball']['imag_ball']},

|P_core| = {c['grouped_endpoint_saddle_core_ball']['absolute_ball']}.
```

This yields the exact split `T_upper=P_core+R_upper`.  Every phase retained in
`R_upper` is nonstationary on its own interval: the 230 labels continue from
`B`, and labels from `q_++230` upward start at `L`.

Pi provenance: every `pi` comes from the inherited Fresnel phase, exact
half-integer Fourier spacing, or `p=t/(2*pi)`.  The rational split is
`L+1/16`; no fitted geometric constant is used.

Proof boundary: exact 230-label endpoint cancellation and rational split,
production saddle census, and rigorous grouped complex enclosure of `P_core`
only.  No quantitative enclosure of `R_upper`, the lower complementary tail,
complete ordinary or joined packet, `J_Z`, or `D_K` is proved, nor is any
non-A, all-height, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
"""


def main() -> int:
    resource_mode = set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = PRECISION_BITS
    ctx.threads = 1
    require(COMPLEMENT.is_file(), f"missing dependency: {COMPLEMENT}")
    require(load_json(COMPLEMENT).get("passed") is True, "complementary gate not passed")
    require(CHECKER.is_file(), f"missing checker: {CHECKER}")

    certificate = build_certificate()
    artifact = {
        "kind": "rh_c1_hardy_equation4_upper_complement_endpoint_saddle_core_gate",
        "date": "2026-08-27",
        "status": "grouped_230_label_endpoint_saddle_core_complex_ball_certified",
        "passed": True,
        "resource_mode": resource_mode,
        "certificate": certificate,
        "decision": {
            "all_230_stationary_roots_captured": True,
            "grouped_endpoint_amplitude_vanishes_exactly": True,
            "stationary_core_complex_ball_certified": True,
            "remaining_upper_piece_nonstationary": True,
            "remaining_upper_piece_enclosed": False,
            "complete_ordinary_packet_enclosed": False,
        },
        "sources": {
            relative(COMPLEMENT): file_hash(COMPLEMENT),
            relative(BUILDER): file_hash(BUILDER),
            relative(CHECKER): file_hash(CHECKER),
        },
        "proof_boundary": (
            "Exact endpoint cancellation, rational split and saddle census, plus a rigorous "
            "complex enclosure of the grouped 230-label stationary core only. No remaining "
            "complementary-tail, ordinary-packet, joined-packet, J_Z, D_K, non-A, all-height, "
            "Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved."
        ),
    }
    atomic_write(RESULT, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE, note_text(artifact))
    print("certified grouped 230-label endpoint-saddle core")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
