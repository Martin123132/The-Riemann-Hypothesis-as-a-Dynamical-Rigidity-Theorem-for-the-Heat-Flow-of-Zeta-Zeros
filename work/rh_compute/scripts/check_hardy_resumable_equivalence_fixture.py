#!/usr/bin/env python3
"""Independently validate resumable Hardy stop/resume fixture artifacts."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import re


REPO_ROOT = Path(__file__).resolve().parents[3]
EXTERNAL_ROOT = REPO_ROOT / "work/rh_compute/external/hardy_fastcode"
FIXTURE_ROOT = EXTERNAL_ROOT / "fixtures/resume_equivalence"
RESULT_PATH = FIXTURE_ROOT / "fixture_result.json"
RESUMABLE_ROOT = EXTERNAL_ROOT / "resumable"


CASES = {
    "t1e10": {
        "accepted": EXTERNAL_ROOT
        / "fixtures/t1e10/zeta14cubicmult_et005_output.txt",
        "totblock": 35,
        "nc": 621,
        "rs_units": 3,
        "park_count": 4,
        "negative_count": 2,
    },
    "t1e12": {
        "accepted": EXTERNAL_ROOT
        / "fixtures/t1e12/zeta14cubicmult_et005_output.txt",
        "totblock": 46,
        "nc": 2034,
        "rs_units": 8,
        "park_count": 2,
        "negative_count": 0,
    },
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_jsonl(path: Path) -> list[dict]:
    records = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        require(line.strip(), f"blank journal line {line_number}: {path}")
        records.append(json.loads(line, parse_float=Decimal))
    return records


def validate_records(records: list[dict], case: dict) -> dict:
    expected_units = (
        [(0, unit) for unit in range(case["totblock"] + 1)]
        + [(1, unit) for unit in range(1, case["rs_units"] + 1)]
        + [(2, 0)]
    )
    units = [(int(record["stage"]), int(record["unit"])) for record in records]
    require(units == expected_units, "checkpoint stage/unit sequence drift")
    run_ids = {record["run_id"] for record in records}
    require(len(run_ids) == 1, "checkpoint run-id drift")
    run_id = next(iter(run_ids))
    require(re.fullmatch(r"[0-9a-f]{64}", run_id) is not None, "run-id format drift")

    for record in records:
        require(record["format"] == "RH_HARDY_CHECKPOINT_V1", "format drift")
        require(record["comm_size"] == 1, "MPI size drift")
        require(record["numbercalc"] == 8 and record["nchalf"] == 15, "grid-size drift")
        require(record["totblock"] == case["totblock"], "block-count drift")
        require(record["rs_chunk_size"] == 257, "tail chunk-size drift")
        round_trip_digits = math.ceil(
            1 + record["real_digits"] * math.log10(record["real_radix"])
        )
        require(
            record["serialized_significand_digits"] >= round_trip_digits,
            "serialized precision is below the exact round-trip threshold",
        )
        for field in ("t", "et", "rn1"):
            value = record[field]
            require(isinstance(value, Decimal) and value.is_finite(), f"{field} finiteness drift")
        for field in ("zsum", "rszsumtot", "rszsum"):
            values = record[field]
            require(len(values) == 15, f"{field} channel count drift")
            require(all(value.is_finite() for value in values), f"{field} finiteness drift")

    tail = records[-2]
    require(tail["stage"] == 1 and tail["unit"] == case["rs_units"], "tail owner drift")
    require(tail["nc"] == case["nc"], "RS cutoff drift")
    require(tail["nmax"] == case["rs_units"] - 1, "RS complete-chunk count drift")
    return {"record_count": len(records), "run_id": run_id}


def parse_snapshot(path: Path) -> dict:
    lines = path.read_text(encoding="utf-8").splitlines()
    require(len(lines) == 9, "snapshot record count drift")
    require(lines[0].strip() == "RH_HARDY_CHECKPOINT_V1", "snapshot magic drift")
    require(re.fullmatch(r"[0-9a-f]{64}", lines[1].strip()) is not None, "snapshot run-id drift")
    integers = [int(value) for value in lines[2].split()]
    require(len(integers) == 6, "snapshot primary metadata width")
    reals = [Decimal(value) for value in lines[3].split()]
    require(len(reals) == 3, "snapshot real metadata width")
    recurrence = [int(value) for value in lines[4].split()]
    require(len(recurrence) == 7, "snapshot recurrence metadata width")
    arrays = [[Decimal(value) for value in lines[index].split()] for index in (5, 6, 7)]
    require(all(len(values) == 15 for values in arrays), "snapshot channel width")
    checksums = [Decimal(value) for value in lines[8].split()]
    require(len(checksums) == 3, "snapshot checksum width")
    require(all(value.is_finite() for values in arrays for value in values), "snapshot non-finite value")
    return {
        "run_id": lines[1].strip(),
        "stage": integers[0],
        "unit": integers[1],
        "comm_size": integers[2],
        "numbercalc": integers[3],
        "nchalf": integers[4],
        "totblock": integers[5],
        "t": reals[0],
        "et": reals[1],
        "rn1": reals[2],
        "mt": recurrence[0],
        "mto": recurrence[1],
        "aenums": recurrence[2],
        "mmax": recurrence[3],
        "nc": recurrence[4],
        "nmax": recurrence[5],
        "rs_chunk_size": recurrence[6],
        "zsum": arrays[0],
        "rszsumtot": arrays[1],
        "rszsum": arrays[2],
    }


def parse_grand_totals(path: Path) -> list[Decimal]:
    values = re.findall(
        r"Grand total of Hardy function Z\(t[^=]*=\s*([-+0-9.Ee]+)",
        path.read_text(encoding="utf-8"),
    )
    require(len(values) == 15, f"Hardy output count drift: {path}")
    return [Decimal(value) for value in values]


def compare_snapshot_to_final(snapshot: dict, final_record: dict) -> None:
    for field in (
        "run_id",
        "stage",
        "unit",
        "comm_size",
        "numbercalc",
        "nchalf",
        "totblock",
        "t",
        "et",
        "rn1",
        "mt",
        "mto",
        "aenums",
        "mmax",
        "nc",
        "nmax",
        "rs_chunk_size",
        "zsum",
        "rszsumtot",
        "rszsum",
    ):
        require(snapshot[field] == final_record[field], f"snapshot/journal drift: {field}")


def main() -> int:
    payload = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    require(payload["kind"] == "hardy_resumable_equivalence_fixture", "result kind drift")
    require(payload["image"] == "rh-hardy-fastcode:bookworm", "container image drift")
    require(
        re.fullmatch(r"sha256:[0-9a-f]{64}", payload["image_id"]) is not None,
        "container image-id drift",
    )
    require(payload["mpi_ranks"] == 1 and payload["cpu_limit"] == 1, "worker policy drift")
    require(
        payload["module_sha256"] == file_hash(RESUMABLE_ROOT / "rh_hardy_checkpoint.f90"),
        "checkpoint module hash drift",
    )
    require(
        payload["source_sha256"]
        == file_hash(RESUMABLE_ROOT / "zeta14cubicmult_resumable.f90"),
        "resumable source hash drift",
    )
    by_case = {item["case"]: item for item in payload["cases"]}
    require(set(by_case) == set(CASES), "fixture case roster drift")

    total_records = 0
    for name, expected in CASES.items():
        item = by_case[name]
        root = FIXTURE_ROOT / name
        uninterrupted_root = root / "uninterrupted"
        resumed_root = root / "resumed"
        uninterrupted_journal_path = uninterrupted_root / "checkpoint.jsonl"
        resumed_journal_path = resumed_root / "checkpoint.jsonl"
        uninterrupted = load_jsonl(uninterrupted_journal_path)
        resumed = load_jsonl(resumed_journal_path)
        left_summary = validate_records(uninterrupted, expected)
        right_summary = validate_records(resumed, expected)
        require(uninterrupted == resumed, f"full-precision resume drift: {name}")
        require(left_summary == right_summary, f"journal summary drift: {name}")
        total_records += 2 * left_summary["record_count"]

        left_snapshot = parse_snapshot(uninterrupted_root / "checkpoint.dat")
        right_snapshot = parse_snapshot(resumed_root / "checkpoint.dat")
        compare_snapshot_to_final(left_snapshot, uninterrupted[-1])
        compare_snapshot_to_final(right_snapshot, resumed[-1])
        require(left_snapshot == right_snapshot, f"final snapshot drift: {name}")

        left_output = parse_grand_totals(uninterrupted_root / "uninterrupted.output.txt")
        right_output = parse_grand_totals(resumed_root / "resumed_complete.output.txt")
        accepted_output = parse_grand_totals(expected["accepted"])
        require(left_output == right_output, f"final displayed output drift: {name}")
        accepted_error = max(abs(left - right) for left, right in zip(left_output, accepted_output))
        require(
            Decimal(item["accepted_output_max_abs_error"]) == accepted_error,
            f"accepted-output error drift: {name}",
        )
        require(item["journal_exact_match"] and item["final_output_exact_match"], "result flag drift")
        require(len(item["resumed_runs"]) == expected["park_count"] + 1, "resume run count drift")
        require(
            sum(run["returncode"] == 75 for run in item["resumed_runs"])
            == expected["park_count"],
            "parked return-code count drift",
        )
        require(item["resumed_runs"][-1]["returncode"] == 0, "completion return-code drift")
        require(
            len(item["negative_resume_tests"]) == expected["negative_count"],
            "negative-test count drift",
        )
        require(
            all(run["returncode"] == 2 for run in item["negative_resume_tests"]),
            "negative-test return-code drift",
        )

        hashes = item["hashes"]
        require(hashes["uninterrupted_journal"] == file_hash(uninterrupted_journal_path), "journal hash drift")
        require(hashes["resumed_journal"] == file_hash(resumed_journal_path), "resumed hash drift")
        for negative_name in ("reject_run_id", "reject_chunk_size"):
            negative_root = root / negative_name
            if negative_root.exists():
                require(
                    file_hash(negative_root / "checkpoint.dat")
                    == file_hash(resumed_root / "checkpoint.dat"),
                    f"negative snapshot mutation: {negative_name}",
                )
                require(
                    file_hash(negative_root / "checkpoint.jsonl")
                    == file_hash(resumed_root / "checkpoint.jsonl"),
                    f"negative journal mutation: {negative_name}",
                )

    print(
        "validated resumable Hardy equivalence fixtures: "
        f"2 heights, {total_records} journal records, exact stop/resume state, "
        "6 parked invocations, 2 fail-closed resumes"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
