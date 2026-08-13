#!/usr/bin/env python3
"""Check the fourth-summand interval continuation certificate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import importlib
import json
import os
from pathlib import Path
import sys

for variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(variable, "1")

import psutil


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

from flint import arb  # noqa: E402


STEM = (
    "jensen_window_pf_newman_theta_fourth_summand_"
    "tangency_interval_continuation_certificate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
SOURCE = SCRIPT_DIR / f"{STEM}.py"
PARENT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_fourth_summand_tangency_homotopy_scout.json"
)


class CheckFailure(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CheckFailure(message)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def lower(text: str) -> arb:
    return arb(text).lower()


def upper(text: str) -> arb:
    return arb(text).upper()


def check_positive_ball(record: dict, label: str) -> None:
    require(lower(record["ball"]) > 0, f"{label} is not strictly positive")
    require(lower(record["lower"]) > 0, f"{label} stored lower endpoint is not positive")


def check_negative_ball(record: dict, label: str) -> None:
    require(upper(record["ball"]) < 0, f"{label} is not strictly negative")
    require(upper(record["upper"]) < 0, f"{label} stored upper endpoint is not negative")


def static_audit(payload: dict) -> dict:
    require(payload["kind"] == STEM, "result kind mismatch")
    require(payload["schema_version"] == 1, "schema version mismatch")
    require(
        payload["status"] == "rigorous interval continuation certificate complete",
        "certificate is not complete",
    )
    require(payload["source_sha256"] == sha256_path(SOURCE), "source hash mismatch")
    require(payload["parent"]["result_sha256"] == sha256_path(PARENT), "parent hash mismatch")

    resource = payload["resource_policy"]
    require(resource["mode"] == "daytime", "resource mode mismatch")
    require(resource["active_compute_workers"] == 1, "worker count is not one")
    require(resource["thread_caps"] == 1, "thread cap is not one")
    require(resource["below_normal_priority_applied"] is True, "below-normal priority missing")
    require(resource["resource_parked"] is False, "certificate was resource-parked")

    config = payload["interval_configuration"]
    require(config["precision_bits"] >= 192, "insufficient Arb precision")
    require(config["taylor_time_order"] >= 10, "insufficient time Taylor order")
    require(config["taylor_x_order"] >= 28, "insufficient x Taylor order")
    require(config["maximum_moment"] <= 60, "tail audit moment range exceeded")
    require(payload["tail_audit"]["component_bound"] == "integral tail < 10^-1800", "tail audit changed")

    continuation = payload["continuation"]
    charts = continuation["charts"]
    require(continuation["lambda_coverage"] == ["0", "0.22"], "coverage mismatch")
    require(continuation["expected_charts"] == 11, "expected chart count mismatch")
    require(continuation["certified_charts"] == 11, "not all charts were certified")
    require(continuation["unresolved_charts"] == 0, "unresolved charts remain")
    require(len(charts) == 11, "stored chart count mismatch")
    require(continuation["all_krawczyk_inclusions_passed"] is True, "Krawczyk summary failed")
    require(continuation["branch_chain_connected"] is True, "branch chain is not connected")

    minimum_krawczyk_margin: arb | None = None
    for index, chart in enumerate(charts):
        require(chart["index"] == index, f"chart {index} index mismatch")
        require(
            Fraction(chart["lambda"]["low"]) == Fraction(index, 50),
            f"chart {index} lower lambda mismatch",
        )
        require(
            Fraction(chart["lambda"]["high"]) == Fraction(index + 1, 50),
            f"chart {index} upper lambda mismatch",
        )
        require(chart["krawczyk"]["passed"] is True, f"chart {index} Krawczyk flag failed")
        for axis, margin in enumerate(chart["krawczyk"]["strict_interior_margins"]):
            check_positive_ball(margin, f"chart {index} Krawczyk margin {axis}")
            value = lower(margin["ball"])
            if minimum_krawczyk_margin is None or value < minimum_krawczyk_margin:
                minimum_krawczyk_margin = value
        signs = chart["signs"]
        require(signs["H4_strictly_positive"] is True, f"chart {index} H4 flag failed")
        require(signs["F_xx_strictly_negative"] is True, f"chart {index} Fxx flag failed")
        require(signs["dt_dlambda_strictly_negative"] is True, f"chart {index} slope flag failed")
        check_positive_ball(signs["H4"], f"chart {index} H4")
        check_negative_ball(signs["F_xx"], f"chart {index} Fxx")
        for derivative, pair in enumerate(chart["taylor_remainders"]):
            for family in ("weighted", "fourth"):
                check_positive_ball(
                    {
                        "ball": pair[family]["total_remainder_upper"],
                        "lower": pair[family]["total_remainder_upper"],
                    },
                    f"chart {index} derivative {derivative} {family} remainder",
                )

    connectors = continuation["connectors"]
    require(len(connectors) == 10, "connector count mismatch")
    minimum_connector_margin: arb | None = None
    for index, connector in enumerate(connectors):
        require(connector["left_chart"] == index, f"connector {index} left index mismatch")
        require(connector["right_chart"] == index + 1, f"connector {index} right index mismatch")
        require(connector["passed"] is True, f"connector {index} failed")
        for name, margin in connector["strict_containment_margins_in_right_chart"].items():
            check_positive_ball(margin, f"connector {index} margin {name}")
            value = lower(margin["ball"])
            if minimum_connector_margin is None or value < minimum_connector_margin:
                minimum_connector_margin = value

    endpoints = continuation["endpoints"]
    require(endpoints["initial_time_strictly_positive"] is True, "initial sign flag failed")
    require(endpoints["final_time_strictly_negative"] is True, "final sign flag failed")
    check_positive_ball(endpoints["initial"]["time"], "initial endpoint time")
    check_negative_ball(endpoints["final"]["time"], "final endpoint time")

    crossing = payload["boundary_crossing"]
    require(crossing["passed"] is True, "crossing Krawczyk failed")
    require(crossing["branch_match_passed"] is True, "crossing branch match failed")
    require(crossing["matched_to_continuation_chart"] == 10, "crossing chart mismatch")
    require(crossing["H4_strictly_positive"] is True, "crossing H4 sign failed")
    require(crossing["F_xx_strictly_negative"] is True, "crossing Fxx sign failed")
    require(crossing["dt_dlambda_strictly_negative"] is True, "crossing slope sign failed")
    for axis, margin in enumerate(crossing["strict_interior_margins"]):
        check_positive_ball(margin, f"crossing Krawczyk margin {axis}")
    require(
        Fraction("0.2") < Fraction(crossing["lambda_center"]) < Fraction("0.22"),
        "crossing center lies outside the final chart",
    )

    proof_boundary = payload["theorem_decision"]["proof_boundary"]
    for phrase in ("does not count all contacts", "n>=5", "RH", "Clay"):
        require(phrase in proof_boundary, f"proof boundary lost phrase: {phrase}")
    return {
        "charts": len(charts),
        "connectors": len(connectors),
        "minimum_krawczyk_margin_lower": str(minimum_krawczyk_margin),
        "minimum_connector_margin_lower": str(minimum_connector_margin),
    }


def independent_anchor_audit() -> dict:
    core = importlib.import_module(STEM)
    core.PRECISION_BITS = 224
    core.ABS_TOL = "1e-38"
    core.TAYLOR_TIME_ORDER = 11
    core.TAYLOR_X_ORDER = 30
    core.flint.ctx.prec = 224
    parent, rows = core.load_parent_rows()
    monitor = core.CpuMonitor()
    anchors: list[dict] = []
    final_chart = None
    for index in (0, 5, 10):
        chart = core.ContinuationChart(index, rows[index], rows[index + 1])
        record = chart.build()
        anchors.append(
            {
                "index": index,
                "time_margin": record["krawczyk"]["strict_interior_margins"][0]["lower"],
                "x_margin": record["krawczyk"]["strict_interior_margins"][1]["lower"],
                "H4_lower": record["signs"]["H4"]["lower"],
                "F_xx_upper": record["signs"]["F_xx"]["upper"],
            }
        )
        if index == 10:
            final_chart = chart
        monitor.sample()
    require(final_chart is not None, "independent final chart was not built")
    crossing = core.certify_crossing(parent, final_chart)
    monitor.sample()
    return {
        "precision_bits": 224,
        "absolute_integration_tolerance": "1e-38",
        "taylor_time_order": 11,
        "taylor_x_order": 30,
        "anchors": anchors,
        "crossing_lambda": crossing["lambda_enclosure"]["ball"],
        "crossing_x": crossing["x_enclosure"]["ball"],
        "runtime_cpu_percent": monitor.samples,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=RESULT)
    parser.add_argument("--baseline-seconds", type=int, default=5)
    parser.add_argument("--baseline-max-percent", type=float, default=60.0)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    core = importlib.import_module(STEM)
    priority_lowered = core.request_below_normal_priority()
    baseline = [float(psutil.cpu_percent(interval=1.0)) for _ in range(args.baseline_seconds)]
    baseline_mean = sum(baseline) / len(baseline)
    print(
        "interval continuation checker baseline CPU: "
        f"samples={','.join(f'{value:.1f}' for value in baseline)}, mean={baseline_mean:.1f}%",
        flush=True,
    )
    if baseline_mean > args.baseline_max_percent:
        print("interval continuation checker deferred: daytime baseline is already busy")
        return 2
    require(priority_lowered, "checker could not apply below-normal priority")
    payload = json.loads(args.result.read_text(encoding="utf-8"))
    static = static_audit(payload)
    independent = independent_anchor_audit()
    print(
        "checked fourth-summand interval continuation: "
        f"charts={static['charts']}, connectors={static['connectors']}, "
        f"independent_anchors={len(independent['anchors'])}, crossing=pass, 0 issues"
    )
    print(json.dumps({"static": static, "independent": independent}, indent=2))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except CheckFailure as error:
        print(f"CHECK FAILED: {error}", file=sys.stderr)
        sys.exit(1)
