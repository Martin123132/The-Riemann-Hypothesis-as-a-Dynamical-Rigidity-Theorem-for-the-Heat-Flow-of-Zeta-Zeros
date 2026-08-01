#!/usr/bin/env python3
"""Stress-test order-eleven compact-to-lambda-zero composition boundaries."""

from __future__ import annotations

import copy
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order11_first_summand_curvature_certificate as global_certificate  # noqa: E402
import jensen_window_pf_compound_order11_lambda0_completion_certificate as completion  # noqa: E402
import jensen_window_pf_compound_order11_m100_entry_certificate as endpoint  # noqa: E402


def expect_runtime(label: str, function, marker: str) -> None:
    try:
        function()
    except RuntimeError as exc:
        if marker not in str(exc):
            raise AssertionError(f"{label}: wrong rejection: {exc}") from exc
    else:
        raise AssertionError(f"{label}: mutation was accepted")


def compact_fixture(*, segments: int = 2020) -> dict:
    return {
        "kind": "jensen_window_pf_compound_order11_compact_adaptive_h23_certificate",
        "status": (
            "rigorous order-eleven first-summand curvature theorem on "
            "5700<=t<=38020"
        ),
        "theorem": global_certificate.COMPACT_THEOREM,
        "summary": {
            "segments": segments,
            "quarter_blocks": segments * 64,
            "compact_first_summand_theorems": 1,
            "global_first_summand_theorems": 0,
            "full_kernel_theorems": 0,
            "largest_scaled_curvature_upper": "2123",
            "saddle_transition_upper": "38019.7",
        },
        "source_contract": {
            "saddle_monotonicity": {
                "formula": "V''>0 for u>=1/100; V'<0 for 0<u<=1/100",
            },
        },
    }


def global_curvature_fixture(*, heat_forward: int = 0) -> dict:
    return {
        "kind": (
            "jensen_window_pf_compound_order11_"
            "first_summand_curvature_certificate"
        ),
        "status": (
            "rigorous global order-eleven first-summand curvature theorem "
            "on t>=1252"
        ),
        "theorem": endpoint.ORDER11_CONTINUOUS,
        "summary": {
            "global_first_summand_curvature_theorems": 1,
            "open_rows": 0,
            "full_kernel_theorems": 0,
            "heat_forward_theorems": heat_forward,
            "rh_claims": 0,
        },
    }


def endpoint_fixture(*, orders_above_11: int = 0) -> dict:
    return {
        "kind": "jensen_window_pf_compound_order11_m100_entry_certificate",
        "status": "rigorous all-shift signed order-eleven entry at lambda=-100",
        "exact": {"global_endpoint": completion.ENDPOINT},
        "summary": {
            "global_m100_order11_entry_theorems": 1,
            "open_rows": 0,
            "heat_interval_theorems": 0,
            "orders_above_11": orders_above_11,
            "rh_claims": 0,
        },
    }


def test_global_composition() -> int:
    original_load = global_certificate.load
    original_sha = global_certificate.sha256
    try:
        invalid = compact_fixture(segments=2019)
        global_certificate.load = lambda path: (
            invalid if path == global_certificate.COMPACT else original_load(path)
        )
        global_certificate.sha256 = lambda path: (
            "0" * 64
            if path == global_certificate.COMPACT
            else original_sha(path)
        )
        expect_runtime(
            "partial compact promotion",
            global_certificate.validate_sources,
            "invalid compact-bridge source",
        )
        valid = compact_fixture()
        global_certificate.load = lambda path: (
            valid if path == global_certificate.COMPACT else original_load(path)
        )
        artifact = global_certificate.build_artifact()
        if (
            artifact["summary"]["global_first_summand_curvature_theorems"] != 1
            or artifact["summary"]["full_kernel_theorems"] != 0
            or artifact["summary"]["rh_claims"] != 0
        ):
            raise AssertionError("global composition claim boundary changed")
    finally:
        global_certificate.load = original_load
        global_certificate.sha256 = original_sha
    return 2


def test_endpoint_composition() -> int:
    original_load = endpoint.load_json
    original_sha = endpoint.sha256
    try:
        invalid = global_curvature_fixture(heat_forward=1)
        endpoint.load_json = lambda path: (
            invalid if path == endpoint.ORDER11_CURVATURE else original_load(path)
        )
        endpoint.sha256 = lambda path: (
            "1" * 64 if path == endpoint.ORDER11_CURVATURE else original_sha(path)
        )
        expect_runtime(
            "endpoint circularity",
            endpoint.validate_sources,
            "global order-eleven curvature source is not closed",
        )
        valid = global_curvature_fixture()
        endpoint.load_json = lambda path: (
            valid if path == endpoint.ORDER11_CURVATURE else original_load(path)
        )
        artifact = endpoint.build_artifact()
        if (
            artifact["summary"]["global_m100_order11_entry_theorems"] != 1
            or artifact["summary"]["circularity_guards"] != 1
            or artifact["summary"]["orders_above_11"] != 0
            or artifact["summary"]["rh_claims"] != 0
        ):
            raise AssertionError("endpoint composition claim boundary changed")
    finally:
        endpoint.load_json = original_load
        endpoint.sha256 = original_sha
    return 2


def test_lambda0_composition() -> int:
    original_load = completion.load_json
    original_sha = completion.sha256
    try:
        invalid = endpoint_fixture(orders_above_11=1)
        completion.load_json = lambda path: (
            invalid if path == completion.ENDPOINT_SOURCE else original_load(path)
        )
        completion.sha256 = lambda path: (
            "2" * 64 if path == completion.ENDPOINT_SOURCE else original_sha(path)
        )
        expect_runtime(
            "higher-order leakage",
            completion.validate_sources,
            "order-eleven endpoint source is not closed",
        )
        valid = endpoint_fixture()
        bad_transfer = copy.deepcopy(original_load(completion.TRANSFER_SOURCE))
        bad_transfer["exact"]["nonpromotion_countermodel"][
            "conclusion"
        ] = "countermodel removed"
        completion.load_json = lambda path: (
            valid
            if path == completion.ENDPOINT_SOURCE
            else bad_transfer
            if path == completion.TRANSFER_SOURCE
            else original_load(path)
        )
        expect_runtime(
            "countermodel removal",
            completion.validate_sources,
            "fixed-order arbitrary-column transfer changed",
        )
        completion.load_json = lambda path: (
            valid if path == completion.ENDPOINT_SOURCE else original_load(path)
        )
        artifact = completion.build_artifact()
        if (
            artifact["summary"]["all_shift_order11_lambda0_theorems"] != 1
            or artifact["summary"]["countermodel_nonpromotion_guards"] != 1
            or artifact["summary"]["orders_above_11"] != 0
            or artifact["summary"]["pf_infinity_theorems"] != 0
            or artifact["summary"]["rh_claims"] != 0
        ):
            raise AssertionError("lambda-zero composition claim boundary changed")
    finally:
        completion.load_json = original_load
        completion.sha256 = original_sha
    return 3


def main() -> int:
    gates = (
        test_global_composition()
        + test_endpoint_composition()
        + test_lambda0_composition()
    )
    print(
        "validated order-eleven completion composition gates: "
        f"{gates} positive/negative gates, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
