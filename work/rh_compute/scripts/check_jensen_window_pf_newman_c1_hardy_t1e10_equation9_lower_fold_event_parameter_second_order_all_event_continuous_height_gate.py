#!/usr/bin/env python3
"""Independently check the all-event continuous-height beta^-4 atlas."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))

import flint
from flint import acb, arb


BUILDER = Path(__file__).with_name(Path(__file__).name.removeprefix("check_"))
CHECK_MOMENT_ORDER = 8


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_builder():
    spec = importlib.util.spec_from_file_location("all_event_height", BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load height builder")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def main() -> None:
    height = load_builder()
    require(height.RESULT.is_file() and height.CACHE.is_file() and height.NOTE.is_file(), "missing height artifacts")
    result = json.loads(height.RESULT.read_text(encoding="utf-8"))
    rows = result["certificate"]["all_rows"]
    require(len(rows) == height.EVENT_COUNT, "wrong height-row count")
    require([row["event_index"] for row in rows] == list(range(height.EVENT_COUNT)), "height roster is incomplete")
    require(result["artifacts"]["builder"]["sha256"] == file_hash(height.BUILDER), "builder hash mismatch")
    require(result["artifacts"]["checker"]["sha256"] == file_hash(Path(__file__)), "checker hash mismatch")
    require(result["artifacts"]["row_cache"]["sha256"] == file_hash(height.CACHE), "cache hash mismatch")
    require(result["artifacts"]["note"]["sha256"] == file_hash(height.NOTE), "note hash mismatch")
    require(
        result["dependencies"]["all_event_center_atlas"]["sha256"] == file_hash(height.ATLAS_RESULT),
        "center-atlas dependency hash mismatch",
    )

    normalized_max = max(rows, key=lambda row: float(arb(row["uniform_normalized_absolute_upper"])))
    physical_max = max(rows, key=lambda row: float(arb(row["uniform_physical_absolute_upper"])))
    require(normalized_max["event_index"] == result["certificate"]["maximum_normalized_event"], "normalized maximum mismatch")
    require(physical_max["event_index"] == result["certificate"]["maximum_physical_event"], "physical maximum mismatch")
    for row in rows:
        index = row["event_index"]
        require(row["config_version"] == height.CONFIG_VERSION, f"config mismatch at event {index}")
        require(row["mode"] == (159_577 - 1) // 4 - index // 2 if index % 2 == 0 else row["mode"] == (159_577 - 1) // 4 + (index + 1) // 2, f"mode mismatch at event {index}")
        require(arb(row["uniform_normalized_absolute_upper"]) < arb("0.000019"), f"normalized target failed at event {index}")
        require(arb(row["uniform_physical_absolute_upper"]) < arb("0.0000086"), f"physical target failed at event {index}")
        require(len(row["moments"]) == height.MOMENT_ORDER + 1, f"moment count mismatch at event {index}")

    atlas, second, first, module, face_class, center_rows = height.load_sources("checker")
    module.SERIES_ORDER = 70
    module.INNER_ORDER = 96
    module.QUADRATIC_ORDER = 30
    flint.ctx.cap = module.SERIES_ORDER + 8
    critical = sorted({0, 1, 2, 198, 199, 200, 397, 398, normalized_max["event_index"]})
    checked = []
    for index in critical:
        mode = atlas.mode_for_event(index)
        face = first.event_face(face_class, module, index, mode)
        moments, mass, _ = height.certified_moments(atlas, module, face, CHECK_MOMENT_ORDER)
        stored = rows[index]
        for degree in range(height.MOMENT_ORDER + 1):
            recorded = parse_complex(stored["moments"][degree]["moment_ball"])
            require(moments[degree].real.overlaps(recorded.real), f"moment {degree} real mismatch at event {index}")
            require(moments[degree].imag.overlaps(recorded.imag), f"moment {degree} imag mismatch at event {index}")

        kappa = face.pi / (16 * face.beta)
        coefficients = [
            (face.i * kappa) ** degree / arb(math.factorial(degree)) * moment
            for degree, moment in enumerate(moments)
        ]
        z_max = (module.HORIZONTAL_END**2 + module.CONTOUR_HEIGHT**2).sqrt()
        size = kappa * z_max
        remainder = mass * size.exp() * size ** (CHECK_MOMENT_ORDER + 1) / arb(math.factorial(CHECK_MOMENT_ORDER + 1))
        far = kappa.exp() * second.second_order_far_bound(first, module, face)
        fresh_bound = sum((abs(value).upper() for value in coefficients), arb(0)) + remainder + far
        require(fresh_bound < height.TARGET_NORMALIZED, f"fresh normalized target failed at event {index}")

        endpoint_overlap = []
        for theta in (-1, 1):
            reconstructed = module.add_error(height.polynomial_value(coefficients, theta), remainder + far)
            direct = height.integrate_height_event(atlas, second, first, module, face_class, index, theta)
            require(reconstructed.real.overlaps(direct.real), f"theta {theta} real mismatch at event {index}")
            require(reconstructed.imag.overlaps(direct.imag), f"theta {theta} imag mismatch at event {index}")
            endpoint_overlap.append(theta)
        checked.append({"event_index": index, "fresh_uniform_bound": fresh_bound.str(80, more=True)})
        print(f"checked event {index:03d} uniform<={fresh_bound}", flush=True)

    print(
        f"validated continuous-height atlas: {len(rows)} cells, max event {normalized_max['event_index']}, "
        f"{len(checked)} nonmatching moment and endpoint reconstructions",
        flush=True,
    )


if __name__ == "__main__":
    main()
