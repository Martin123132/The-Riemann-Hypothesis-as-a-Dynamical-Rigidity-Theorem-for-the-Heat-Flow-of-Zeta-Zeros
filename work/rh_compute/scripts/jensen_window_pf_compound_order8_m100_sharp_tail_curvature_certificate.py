#!/usr/bin/env python3
"""Sharpen the completed order-eight curvature tail for the order-nine bridge."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as order4_compact  # noqa: E402
import jensen_window_pf_compound_order6_nested_curvature_compact_certificate as order6_compact  # noqa: E402
import jensen_window_pf_compound_order7_nested_curvature_compact_certificate as order7_compact  # noqa: E402
import jensen_window_pf_compound_order8_nested_curvature_compact_certificate as order8_compact  # noqa: E402
import jensen_window_pf_compound_order8_nested_curvature_compact_h13_h14_cache as h13_h14_cache  # noqa: E402
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order8_m100_sharp_tail_curvature_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order8_m100_sharp_tail_curvature_certificate.md"
)
FINITE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order8_nested_curvature_finite_ray_certificate.json"
)
ASYMPTOTIC_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order8_nested_curvature_asymptotic_ray_certificate.json"
)
BRIDGE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order8_first_summand_curvature_bridge.json"
)
TARGET_T = 1249
SCALED_TARGET = Decimal(2500)
INITIAL_CENTRAL_TILE_COUNT = 1000


@dataclass(frozen=True)
class SharpRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def target_tile(records: list[dict]) -> int:
    target = flint.arb(TARGET_T)
    for index, record in enumerate(records):
        left = order4_compact.interval_from_text(record["t_left"])
        right = order4_compact.interval_from_text(record["t_right"])
        if bool(right > target):
            if not bool(left < target):
                raise RuntimeError("sharp order-eight target tile does not straddle t=1249")
            return index
    raise RuntimeError("compact H cache does not reach t=1249")


def certify_sharp_block(
    base: list[dict],
    high: list[dict],
    top: list[dict],
    ultra: list[dict],
    left: int,
    right: int,
) -> list[dict]:
    result = order8_compact.compact_block(
        base, high, top, ultra, left, right
    )
    upper = result.get("scaled_curvature_upper")
    if result.get("passed") and upper is not None and Decimal(upper) < SCALED_TARGET:
        return [result]
    if right - left <= 1:
        raise RuntimeError(
            "sharp order-eight compact failed on irreducible tile "
            f"{left}: upper={upper}, failure={result.get('failure')}"
        )
    midpoint = (left + right) // 2
    return certify_sharp_block(
        base, high, top, ultra, left, midpoint
    ) + certify_sharp_block(base, high, top, ultra, midpoint, right)


def sharp_compact_certificate(
    base: list[dict],
    high: list[dict],
    top: list[dict],
    ultra: list[dict],
) -> dict:
    flint.ctx.prec = order4_compact.DEFAULT_PRECISION_BITS
    central_left = target_tile(base)
    central_right = order4_compact.mode_index(Fraction(2))
    accepted: list[dict] = []
    cursor = central_left
    while cursor < central_right:
        endpoint = min(
            cursor + INITIAL_CENTRAL_TILE_COUNT, central_right
        )
        accepted.extend(
            certify_sharp_block(
                base, high, top, ultra, cursor, endpoint
            )
        )
        cursor = endpoint
    for previous, current in zip(accepted, accepted[1:]):
        if previous["central_mode"][1] != current["central_mode"][0]:
            raise RuntimeError("sharp order-eight compact cover has a mode gap")
    if accepted[-1]["central_mode"][1] != "2":
        raise RuntimeError("sharp order-eight compact cover does not reach mode u=2")

    start_mode = Fraction(base[central_left]["mode_left"])
    start_t = potential_jet_arb(arb_rational(start_mode), 1)[1]
    first_right_t = order4_compact.interval_from_text(
        base[central_left]["t_right"]
    )
    if not bool(start_t < TARGET_T < first_right_t):
        raise RuntimeError("first sharp block does not cover t=1249")

    largest = max(
        accepted, key=lambda row: Decimal(row["scaled_curvature_upper"])
    )
    weakest = min(accepted, key=lambda row: Decimal(row["U_lower"]))
    return {
        "target_t": TARGET_T,
        "scaled_target": int(SCALED_TARGET),
        "mode_range": [str(start_mode), "2"],
        "start_t_ball": start_t.str(40).replace("e", "E"),
        "end_t_ball": potential_jet_arb(flint.arb(2), 1)[1]
        .str(40)
        .replace("e", "E"),
        "blocks": accepted,
        "block_count": len(accepted),
        "all_blocks_passed": True,
        "largest_scaled_curvature_upper": largest["scaled_curvature_upper"],
        "largest_scaled_mode": largest["central_mode"],
        "weakest_U_lower": weakest["U_lower"],
        "weakest_U_mode": weakest["central_mode"],
        "theorem": "s_1''(t)<=2500/t^2 for every real 1249<=t<=V'(2)",
    }


def validate_external_sources(
    extension_path: Path,
    collar_path: Path,
) -> dict:
    finite = load_json(FINITE_SOURCE)
    asymptotic = load_json(ASYMPTOTIC_SOURCE)
    bridge = load_json(BRIDGE_SOURCE)
    finite_ray = finite.get("finite_ray", {})
    asymptotic_interval = asymptotic.get("dimensionless_interval", {})
    if (
        finite.get("kind")
        != "jensen_window_pf_compound_order8_nested_curvature_finite_ray_certificate"
        or finite_ray.get("all_blocks_passed") is not True
        or finite_ray.get("ray_blocks") != 17999
        or Decimal(finite_ray.get("largest_scaled_curvature_upper", "Infinity"))
        >= SCALED_TARGET
    ):
        raise RuntimeError("order-eight finite-ray source changed")
    if (
        asymptotic.get("kind")
        != "jensen_window_pf_compound_order8_nested_curvature_asymptotic_ray_certificate"
        or asymptotic_interval.get("scaled_target") != 200
        or Decimal(asymptotic_interval.get("scaled_curvature_upper", "Infinity"))
        >= Decimal(200)
    ):
        raise RuntimeError("order-eight asymptotic source changed")
    transfer = bridge.get("exact", {}).get("full_transfer")
    if transfer != (
        "|W_k-W_k^(1)|<=C_(k-1)+2*C_k+C_(k+1)<190/k^2, k>=1250"
    ):
        raise RuntimeError("order-eight full-kernel transfer changed")
    cache_contract = order8_compact.source_contract(
        extension_path, collar_path
    )
    cache_contract.update(
        {
            "finite_source": FINITE_SOURCE.relative_to(REPO_ROOT).as_posix(),
            "finite_source_sha256": order8_compact.sha256(FINITE_SOURCE),
            "finite_scaled_upper": finite_ray["largest_scaled_curvature_upper"],
            "asymptotic_source": ASYMPTOTIC_SOURCE.relative_to(REPO_ROOT).as_posix(),
            "asymptotic_source_sha256": order8_compact.sha256(ASYMPTOTIC_SOURCE),
            "asymptotic_scaled_upper": asymptotic_interval[
                "scaled_curvature_upper"
            ],
            "bridge_source": BRIDGE_SOURCE.relative_to(REPO_ROOT).as_posix(),
            "bridge_source_sha256": order8_compact.sha256(BRIDGE_SOURCE),
            "full_transfer": transfer,
        }
    )
    return cache_contract


def exact_consequences() -> dict:
    return {
        "compact_theorem": "s_1''(t)<=2500/t^2 for every real 1249<=t<=V'(2)",
        "finite_ray_rebound": (
            "t^2*s_1''(t)<=355.867<2500 for every saddle mode 2<=u<=20"
        ),
        "asymptotic_rebound": "t^2*s_1''(t)<200<2500 for u>=20",
        "global_first_curvature": (
            "s_1''(t)<=2500/t^2 for every real t>=1249"
        ),
        "tent_transfer": (
            "W_k^(1)<=2500*[-log(1-1/k^2)]<2501/k^2, k>=1250"
        ),
        "tent_comparison": (
            "-log(1-x)<x/(1-x) and 2500/(k^2-1)<2501/k^2 "
            "because k^2>2501"
        ),
        "full_transfer": (
            "|W_k-W_k^(1)|<190/k^2, k>=1250"
        ),
        "sharp_full_ceiling": (
            "W_k<2501/k^2+190/k^2=2691/k^2, k>=1250"
        ),
    }


def build_artifact(
    base: list[dict],
    high: list[dict],
    top: list[dict],
    ultra: list[dict],
    extension_path: Path,
    collar_path: Path,
) -> dict:
    sources = validate_external_sources(extension_path, collar_path)
    compact = sharp_compact_certificate(base, high, top, ultra)
    exact = exact_consequences()
    rows = [
        SharpRow(
            "co8mstc_01_aligned_cache",
            "interval_input",
            "ready_to_apply",
            "Aligned compact caches enclose H derivatives through order fourteen on a strict common collar.",
            "H^(2),...,H^(14) on t+-6",
            "First Newman summand at lambda=-100 only.",
        ),
        SharpRow(
            "co8mstc_02_sharp_compact",
            "interval_theorem",
            "ready_to_apply",
            "Adaptive subdivision recovers the order-eight curvature reserve discarded by the broad compact cover.",
            exact["compact_theorem"],
            "Outward-rounded common-collar interval arithmetic.",
            {
                "blocks": compact["block_count"],
                "largest_scaled_upper": compact[
                    "largest_scaled_curvature_upper"
                ],
            },
        ),
        SharpRow(
            "co8mstc_03_existing_rays",
            "interval_composition",
            "ready_to_apply",
            "The completed finite and asymptotic ray bounds are already far below 2500.",
            exact["finite_ray_rebound"] + "; " + exact["asymptotic_rebound"],
            "Reuses completed order-eight ray certificates.",
        ),
        SharpRow(
            "co8mstc_04_global_first",
            "analytic_theorem",
            "ready_to_apply",
            "The compact and two ray regimes cover the complete sharp first-summand tail.",
            exact["global_first_curvature"],
            "Lambda=-100 first summand only.",
        ),
        SharpRow(
            "co8mstc_05_tent",
            "exact_discrete_transfer",
            "ready_to_apply",
            "The tent identity transfers the sharp continuous theorem to integer centered curvature.",
            exact["tent_transfer"] + "; " + exact["tent_comparison"],
            "Exact logarithmic comparison after the interval theorem.",
        ),
        SharpRow(
            "co8mstc_06_full_ceiling",
            "analytic_theorem",
            "ready_to_apply",
            "The completed first/full transfer gives a much sharper complete-kernel order-eight ceiling.",
            exact["sharp_full_ceiling"],
            "Complete theta kernel at lambda=-100 for integer k>=1250.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order8_m100_sharp_tail_curvature_certificate",
        "date": "2026-07-13",
        "status": "rigorous sharp order-eight tail curvature theorem for the order-nine bridge",
        "proof_boundary": (
            "This artifact sharpens the completed order-eight curvature tail at "
            "lambda=-100. It does not prove order-nine entry, PF-infinity, RH, "
            "or Lambda<=0."
        ),
        "source_contract": sources,
        "parameters": {
            "target_t": TARGET_T,
            "scaled_target": int(SCALED_TARGET),
            "initial_central_tile_count": INITIAL_CENTRAL_TILE_COUNT,
            "precision_bits": order4_compact.DEFAULT_PRECISION_BITS,
        },
        "compact": compact,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": len(rows),
            "open_rows": 0,
            "compact_blocks": compact["block_count"],
            "sharp_compact_theorems": 1,
            "global_sharp_curvature_theorems": 1,
            "sharp_full_ceiling_theorems": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order8_m100_sharp_tail_curvature_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order8_m100_sharp_tail_curvature_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    compact = artifact["compact"]
    exact = artifact["exact"]
    lines = [
        "# Order-Eight Sharp Tail Curvature Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous order-eight sharp tail theorem at `lambda=-100`.",
        "This is not a proof of order-nine entry, PF-infinity, RH, or",
        "`Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order8_m100_sharp_tail_curvature_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order8_m100_sharp_tail_curvature_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order8_m100_sharp_tail_curvature_certificate.py",
        "```",
        "",
        "## Recovered Reserve",
        "",
        "The same aligned H2-H14 caches are adaptively repartitioned from",
        "`t=1249` to saddle mode two. Every common-collar block proves all",
        "five stable-log arguments positive before bounding the curvature.",
        "",
        "```text",
        exact["compact_theorem"],
        f"all {compact['block_count']} sharp compact blocks pass,",
        "largest scaled upper="
        + compact["largest_scaled_curvature_upper"]
        + "<2500,",
        exact["finite_ray_rebound"],
        exact["asymptotic_rebound"],
        exact["global_first_curvature"],
        "```",
        "",
        "## Complete Kernel",
        "",
        "```text",
        exact["tent_transfer"],
        exact["tent_comparison"],
        exact["full_transfer"],
        exact["sharp_full_ceiling"],
        "```",
        "",
        "The earlier 4300 ceiling was sufficient for order eight but too loose",
        "for the sixth stable-log floor. The certified 2691 ceiling preserves",
        "the reserve needed by the order-nine first/full bridge.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order8_m100_entry_certificate.md",
        "outputs/jensen_window_pf_negative_lambda_first_summand_power10_rebalanced_dominance_extension.md",
        "outputs/jensen_window_pf_compound_order9_m100_tail_curvature_reduction.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--right-collar", type=Path, default=order8_compact.DEFAULT_RIGHT_COLLAR
    )
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()

    tasks = order8_compact.deterministic_tasks()
    base = order4_compact.load_cache(order4_compact.DEFAULT_CACHE, tasks)
    high = order6_compact.load_extension_cache(
        order6_compact.DEFAULT_EXTENSION_CACHE, tasks
    )
    top = order7_compact.load_extension_cache(
        order7_compact.DEFAULT_EXTENSION_CACHE, tasks
    )
    ultra = h13_h14_cache.load_cache(h13_h14_cache.DEFAULT_CACHE, tasks)
    collar = order8_compact.build_right_collar(args.right_collar)
    base, high, top, ultra = order8_compact.append_right_collar(
        base, high, top, ultra, collar
    )
    artifact = build_artifact(
        base,
        high,
        top,
        ultra,
        h13_h14_cache.DEFAULT_CACHE,
        args.right_collar,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    print(
        "wrote order-eight sharp tail curvature: "
        f"{artifact['summary']['compact_blocks']} compact blocks, "
        f"scaled upper {artifact['compact']['largest_scaled_curvature_upper']}, "
        "1 sharp full-kernel ceiling"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
