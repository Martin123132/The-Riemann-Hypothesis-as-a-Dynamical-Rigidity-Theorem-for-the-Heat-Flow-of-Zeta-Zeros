#!/usr/bin/env python3
"""Independently check the all-event beta^-4 normal-form atlas."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))

from flint import acb, arb


BUILDER = Path(__file__).with_name(Path(__file__).name.removeprefix("check_"))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_builder():
    spec = importlib.util.spec_from_file_location("all_event_atlas", BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load all-event builder")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def main() -> None:
    atlas = load_builder()
    require(atlas.RESULT.is_file() and atlas.CACHE.is_file() and atlas.NOTE.is_file(), "missing all-event artifacts")
    result = json.loads(atlas.RESULT.read_text(encoding="utf-8"))
    rows = result["certificate"]["all_rows"]
    require(len(rows) == atlas.EVENT_COUNT, "wrong event count")
    require([row["event_index"] for row in rows] == list(range(atlas.EVENT_COUNT)), "event roster is incomplete")
    require(result["artifacts"]["builder"]["sha256"] == file_hash(atlas.BUILDER), "builder hash mismatch")
    require(result["artifacts"]["checker"]["sha256"] == file_hash(Path(__file__)), "checker hash mismatch")
    require(result["artifacts"]["row_cache"]["sha256"] == file_hash(atlas.CACHE), "cache hash mismatch")
    require(result["artifacts"]["note"]["sha256"] == file_hash(atlas.NOTE), "note hash mismatch")
    require(
        result["dependencies"]["extreme_second_order_gate"]["sha256"] == file_hash(atlas.SECOND_RESULT),
        "extreme-gate dependency hash mismatch",
    )

    normalized_max = arb(0)
    normalized_max_index = -1
    negative_max = (-1, arb(0))
    positive_max = (-1, arb(0))
    for row in rows:
        index = row["event_index"]
        mode = atlas.mode_for_event(index)
        require(row["mode"] == mode, f"mode mismatch at event {index}")
        require(row["signed_odd"] == atlas.signed_odd_for_event(index), f"signed odd mismatch at event {index}")
        require(row["config_version"] == atlas.CONFIG_VERSION, f"config mismatch at event {index}")
        value = arb(row["normalized_second_order_absolute_upper"])
        physical = arb(row["physical_second_order_absolute_upper"])
        require(value < arb("0.000019"), f"normalized target failed at event {index}")
        require(physical < arb("0.0000086"), f"physical target failed at event {index}")
        if value > normalized_max:
            normalized_max = value
            normalized_max_index = index
        branch = negative_max if index % 2 == 0 else positive_max
        if value > branch[1]:
            if index % 2 == 0:
                negative_max = (index, value)
            else:
                positive_max = (index, value)

    require(normalized_max_index == result["certificate"]["maximum_normalized_event"], "maximum event mismatch")
    require(
        result["certificate"]["maximum_normalized_absolute_upper"]
        == rows[normalized_max_index]["normalized_second_order_absolute_upper"],
        "maximum upper-bound row mismatch",
    )

    second = atlas.load_module(atlas.SECOND_BUILDER, "checker_second_order_source")
    first = second.load_first()
    first.set_low_priority()
    import flint

    flint.ctx.dps = 80
    module = first.load_full_module()
    atlas.configure(module, checker=True)
    face_class = atlas.make_atlas_face(second, first, module)
    critical = sorted({0, 1, 198, 199, 200, negative_max[0], positive_max[0], normalized_max_index, 397, 398})
    checked = []
    for index in critical:
        fresh = atlas.integrate_event(second, first, module, face_class, index)
        stored = rows[index]
        fresh_ball = parse_complex(fresh["normalized_second_order_difference_ball"])
        stored_ball = parse_complex(stored["normalized_second_order_difference_ball"])
        require(fresh_ball.real.overlaps(stored_ball.real), f"real-ball mismatch at event {index}")
        require(fresh_ball.imag.overlaps(stored_ball.imag), f"imag-ball mismatch at event {index}")
        checked.append({
            "event_index": index,
            "fresh_normalized_absolute_upper": fresh["normalized_second_order_absolute_upper"],
        })
        print(f"checked event {index:03d} |D2|<={fresh['normalized_second_order_absolute_upper']}", flush=True)

    print(
        f"validated all-event atlas: {len(rows)} rows, max event {normalized_max_index}, "
        f"{len(checked)} independent nonmatching reconstructions",
        flush=True,
    )


if __name__ == "__main__":
    main()
