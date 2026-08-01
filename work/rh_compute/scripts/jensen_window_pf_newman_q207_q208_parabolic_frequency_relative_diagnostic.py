#!/usr/bin/env python3
"""Compare rigorous relative-transport costs on the Q207-Q208 collar."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402
from flint import arb  # noqa: E402


STEM = "jensen_window_pf_newman_q207_q208_parabolic_frequency_relative_diagnostic"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
DATE = "2026-07-26"

RESULT_DIR = REPO_ROOT / "work/rh_compute/results"
FORWARD = (
    RESULT_DIR
    / "jensen_window_pf_newman_q207_q208_"
    "forward_adiabatic_successor_certificate.json"
)
COARSE = (
    RESULT_DIR
    / "jensen_window_pf_newman_q207_q208_"
    "adiabatic_bottom_collar_certificate.jsonl"
)
REFINED = (
    RESULT_DIR
    / "jensen_window_pf_newman_q207_q208_"
    "adiabatic_bottom_collar_refined_tail_certificate.jsonl"
)
SCALE_REDUCTION = (
    RESULT_DIR
    / "jensen_window_pf_newman_"
    "ray_aligned_parabolic_frequency_reduction.json"
)

SOURCE_FILES = {
    "forward_successor": FORWARD,
    "coarse_transport_cache": COARSE,
    "refined_transport_cache": REFINED,
    "scale_reduction": SCALE_REDUCTION,
}

OLD_TIME = Fraction(1, 1035)
NEW_TIME = Fraction(1, 1040)
TIME_STEP = OLD_TIME - NEW_TIME
CORE_END = Fraction(38)
COARSE_END = Fraction(379, 2)
OLD_RADIUS = Fraction(245)
SCALE_NAMES = ("unit", "frequency", "parabolic", "blended")
PRECISION_BITS = 256
DISPLAY_DIGITS = 55


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def source_hashes() -> dict[str, str]:
    missing = [str(path) for path in SOURCE_FILES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing source artifacts: {missing}")
    return {name: file_hash(path) for name, path in SOURCE_FILES.items()}


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_jsonl(path: Path) -> list[dict]:
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def fraction_arb(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def serialize(value) -> str:
    if not value.is_finite():
        raise RuntimeError(f"nonfinite interval: {value}")
    return value.str(DISPLAY_DIGITS, more=True)


def symmetric_padding(radius) -> arb:
    return arb(0, radius.str(DISPLAY_DIGITS + 20))


def logarithm_at(x: Fraction) -> arb:
    return (fraction_arb(x) / (4 * arb.pi())).log()


def scale_bounds(
    name: str,
    x_low: Fraction,
    x_high: Fraction,
) -> tuple[object, object, object]:
    t_old = fraction_arb(OLD_TIME)
    t_new = fraction_arb(NEW_TIME)
    ell_low = logarithm_at(x_low)
    ell_high = logarithm_at(x_high)

    if name == "unit":
        return arb(1).lower(), arb(1).upper(), arb(0).upper()
    if name == "frequency":
        return (
            (1 / ell_high).lower(),
            (1 / ell_low).upper(),
            arb(0).upper(),
        )
    if name == "parabolic":
        return (
            (2 * t_new).sqrt().lower(),
            (2 * t_old).sqrt().upper(),
            (1 / (2 * t_new).sqrt()).upper(),
        )
    if name == "blended":
        def scale(time: arb, ell: arb) -> arb:
            return (2 * time / (1 + 2 * time * ell**2)).sqrt()

        def scale_time(time: arb, ell: arb) -> arb:
            q = 2 * time * ell**2
            return scale(time, ell) / (2 * time * (1 + q))

        return (
            scale(t_new, ell_high).lower(),
            scale(t_old, ell_low).upper(),
            scale_time(t_new, ell_low).upper(),
        )
    raise ValueError(f"unknown scale {name}")


def transport_lookup() -> dict[tuple[str, Fraction, Fraction], dict]:
    lookup: dict[tuple[str, Fraction, Fraction], dict] = {}
    for source, path, low, high in (
        ("coarse", COARSE, CORE_END, COARSE_END),
        ("refined", REFINED, COARSE_END, OLD_RADIUS),
    ):
        for record in load_jsonl(path):
            x_low = Fraction(record["task"]["x_low"])
            x_high = Fraction(record["task"]["x_high"])
            if low <= x_low and x_high <= high:
                if record["result"].get("status") != "certified":
                    raise RuntimeError(
                        f"uncertified source panel [{x_low},{x_high}]"
                    )
                key = (source, x_low, x_high)
                if key in lookup:
                    raise RuntimeError(f"duplicate transport row {key}")
                lookup[key] = record["result"]["transport"]
    return lookup


def panel_scale_cost(
    forward: dict,
    transport: dict,
    scale_name: str,
) -> dict:
    x_low = Fraction(forward["x_low"])
    x_high = Fraction(forward["x_high"])
    f_old = arb(forward["q207_old_cell_f"])
    fx_old = arb(forward["q207_old_cell_f_prime"])
    f_t_upper = arb(transport["f_t_abs_upper"]).upper()
    fx_t_upper = arb(transport["f_xt_abs_upper"]).upper()
    delta = fraction_arb(TIME_STEP)

    f_collar = f_old + symmetric_padding((delta * f_t_upper).upper())
    fx_collar = fx_old + symmetric_padding((delta * fx_t_upper).upper())
    s_min, s_max, s_t_max = scale_bounds(scale_name, x_low, x_high)

    f_lower = f_collar.abs_lower()
    fx_lower = fx_collar.abs_lower()
    fx_upper = fx_collar.abs_upper().upper()
    denominator = (
        f_lower**2 + arb(s_min) ** 2 * fx_lower**2
    ).sqrt().lower()
    if denominator <= 0:
        raise RuntimeError(
            f"{scale_name} denominator unresolved on [{x_low},{x_high}]"
        )

    second_upper = s_t_max * fx_upper + s_max * fx_t_upper
    numerator = (
        arb(f_t_upper) ** 2 + arb(second_upper) ** 2
    ).sqrt().upper()
    kappa = (arb(numerator) / denominator).upper()
    cost = (delta * kappa).upper()
    return {
        "denominator_lower": serialize(arb(denominator)),
        "time_derivative_upper": serialize(arb(numerator)),
        "kappa_upper": serialize(arb(kappa)),
        "one_step_log_cost_upper": serialize(arb(cost)),
    }


def build_diagnostics() -> tuple[list[dict], dict]:
    forward_artifact = load_json(FORWARD)
    forward_rows = forward_artifact.get("records", [])
    if len(forward_rows) != 525:
        raise RuntimeError(f"expected 525 forward rows, found {len(forward_rows)}")
    lookup = transport_lookup()
    records: list[dict] = []
    maxima: dict[str, tuple[object, int]] = {}
    minima: dict[str, tuple[object, int]] = {}

    for index, forward in enumerate(forward_rows):
        source = forward["source"]
        x_low = Fraction(forward["x_low"])
        x_high = Fraction(forward["x_high"])
        key = (source, x_low, x_high)
        if key not in lookup:
            raise RuntimeError(f"missing transport source for {key}")
        scale_costs = {
            name: panel_scale_cost(forward, lookup[key], name)
            for name in SCALE_NAMES
        }
        record = {
            "sequence": index + 1,
            "source": source,
            "x_low": str(x_low),
            "x_high": str(x_high),
            "scale_costs": scale_costs,
        }
        records.append(record)
        for name in SCALE_NAMES:
            cost = arb(
                scale_costs[name]["one_step_log_cost_upper"]
            ).upper()
            denominator = arb(
                scale_costs[name]["denominator_lower"]
            ).lower()
            if name not in maxima or cost > maxima[name][0]:
                maxima[name] = (cost, index)
            if name not in minima or denominator < minima[name][0]:
                minima[name] = (denominator, index)

    if len(lookup) != len(records):
        raise RuntimeError(
            f"transport/forward row mismatch: {len(lookup)} vs {len(records)}"
        )

    scale_summary: dict[str, dict] = {}
    for name in SCALE_NAMES:
        max_value, max_index = maxima[name]
        min_value, min_index = minima[name]
        max_row = records[max_index]
        min_row = records[min_index]
        scale_summary[name] = {
            "maximum_one_step_log_cost_upper": serialize(arb(max_value)),
            "maximum_cost_panel": [
                max_row["x_low"],
                max_row["x_high"],
            ],
            "minimum_collar_denominator_lower": serialize(arb(min_value)),
            "minimum_denominator_panel": [
                min_row["x_low"],
                min_row["x_high"],
            ],
            "all_525_denominators_positive": True,
            "finite_relative_certificate": True,
        }

    blended_max = arb(
        scale_summary["blended"]["maximum_one_step_log_cost_upper"]
    ).upper()
    parabolic_max = arb(
        scale_summary["parabolic"]["maximum_one_step_log_cost_upper"]
    ).upper()
    frequency_max = arb(
        scale_summary["frequency"]["maximum_one_step_log_cost_upper"]
    ).upper()
    summary = {
        "panels": len(records),
        "cover": [str(CORE_END), str(OLD_RADIUS)],
        "old_time": str(OLD_TIME),
        "new_time": str(NEW_TIME),
        "time_step": str(TIME_STEP),
        "scales": scale_summary,
        "blended_vs_parabolic_max_cost_ratio_upper": serialize(
            arb(blended_max / parabolic_max)
        ),
        "blended_vs_frequency_max_cost_ratio_upper": serialize(
            arb(blended_max / frequency_max)
        ),
        "rigorous_finite_relative_certificates": len(SCALE_NAMES),
        "all_scales_all_panels_finite": True,
        "all_j_promotions": 0,
    }
    return records, summary


def render_note(artifact: dict) -> str:
    summary = artifact["summary"]
    table_lines = [
        "| scale | max one-step relative cost | max panel | "
        "min collar denominator | min panel |",
        "|---|---:|---|---:|---|",
    ]
    for name in SCALE_NAMES:
        row = summary["scales"][name]
        table_lines.append(
            f"| {name} | {row['maximum_one_step_log_cost_upper']} | "
            f"[{row['maximum_cost_panel'][0]},{row['maximum_cost_panel'][1]}] | "
            f"{row['minimum_collar_denominator_lower']} | "
            f"[{row['minimum_denominator_panel'][0]},"
            f"{row['minimum_denominator_panel'][1]}] |"
        )
    table = "\n".join(table_lines)
    return f"""# Q207-Q208 Parabolic-Frequency Relative Diagnostic

Date: {DATE}

Status: rigorous finite source-reusing relative-transport diagnostic on the
first proved forward successor. This is not an all-stage estimate and not a
proof of `Lambda<=0` or RH.

## Source And Scope

The diagnostic reuses the 525 already certified transport panels covering

```text
[1/1040,1/1035]x[38,245].
```

It performs no new Xi quadrature. The old-edge boxes are for
`F=16(1+x^4)H`, and the stored collar bounds enclose `F_t` and `F_(xt)`.

For a positive scale `s`, put `V_s=(F,sF_x)`. On each panel the stored
componentwise transport bounds give rigorous collar boxes for `F` and
`F_x`, followed by

```text
||partial_t V_s||
 <=sqrt(|F_t|^2+(|s_t||F_x|+s|F_(xt)|)^2),

||V_s||
 >=sqrt(dist(F,0)^2+s_min^2 dist(F_x,0)^2).         (1)
```

Every denominator in (1) is strictly positive.

## Four Scales

```text
unit:       s=1,
frequency:  s=1/L,
parabolic:  s=sqrt(2t),
blended:    s=(L^2+(2t)^(-1))^(-1/2).
```

{table}

The ratio of the blended maximum cost to the parabolic maximum cost is
`{summary["blended_vs_parabolic_max_cost_ratio_upper"]}`. Its ratio to the
frequency maximum cost is
`{summary["blended_vs_frequency_max_cost_ratio_upper"]}`.

## Interpretation

At Q207-Q208, `2tL^2` is small, so the blended scale should be close to the
parabolic scale. This finite audit checks that expectation using rigorous
old-cell and heat-jet enclosures. A finite one-step relative certificate is
not an asymptotic Xi bound; it only tests the conditioning of the proposed
coordinate on the first proved successor.

## Proof Boundary

All {summary["panels"]} finite panel inequalities and all four scale
comparisons are rigorous consequences of stored certified interval data.
They prove no Q209 stage, no uniform old-collar relative theorem, no
ray-aligned all-stage successor, no `Lambda<=0`, no RH, and no Clay-prize
result.

Machine-audited files:

```text
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```
"""


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    records, summary = build_diagnostics()
    return {
        "kind": STEM,
        "date": DATE,
        "status": (
            "rigorous finite Q207-Q208 four-scale relative-transport "
            "diagnostic with zero all-j promotions"
        ),
        "proof_boundary": (
            "This source-reusing interval diagnostic proves finite relative "
            "bounds for four scales on the Q207-Q208 collar only. It proves "
            "no Q209 stage, all-j relative estimate, Lambda<=0, or RH."
        ),
        "contract": {
            "precision_bits": PRECISION_BITS,
            "old_time": str(OLD_TIME),
            "new_time": str(NEW_TIME),
            "time_step": str(TIME_STEP),
            "cover": [str(CORE_END), str(OLD_RADIUS)],
            "scales": list(SCALE_NAMES),
            "source_sha256": source_hashes(),
        },
        "summary": summary,
        "records": records,
        "builder_sha256": file_hash(Path(__file__)),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    summary = artifact["summary"]
    print(
        "built Q207-Q208 parabolic-frequency relative diagnostic: "
        f"{summary['panels']} panels, "
        f"{summary['rigorous_finite_relative_certificates']} scales, "
        f"{summary['all_j_promotions']} all-j promotions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
