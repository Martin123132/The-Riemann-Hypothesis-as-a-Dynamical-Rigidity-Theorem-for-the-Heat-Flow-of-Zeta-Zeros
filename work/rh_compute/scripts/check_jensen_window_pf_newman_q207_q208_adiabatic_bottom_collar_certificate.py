#!/usr/bin/env python3
"""Check the Q207-to-Q208 adiabatic bottom-collar certificate."""

from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import jensen_window_pf_newman_q207_q208_adiabatic_bottom_collar_certificate as collar


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT_PATH = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_q207_q208_adiabatic_bottom_collar_certificate.json"
)
NOTE_PATH = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_newman_q207_q208_adiabatic_bottom_collar_certificate.md"
)
compact = collar.compact
HISTORICAL_SUCCESSOR_SHA256 = (
    "c9451ba387d3f0e3f17a914719fc1ec899cfd41bf1667beb19073d1b89c29e90"
)


def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def validate_snapshot_cache_record(
    record: dict,
    task: dict,
    previous_hash: str,
    contract_sha256: str,
) -> None:
    if record.get("schema") != collar.SCHEMA:
        raise RuntimeError("cache schema mismatch")
    if record.get("contract_sha256") != contract_sha256:
        raise RuntimeError("cache contract hash mismatch")
    if record.get("sequence") != task["sequence"]:
        raise RuntimeError("cache sequence mismatch")
    if record.get("task") != task:
        raise RuntimeError("cache task mismatch")
    if record.get("previous_hash") != previous_hash:
        raise RuntimeError("cache hash chain is broken")
    claimed = record.get("record_sha256")
    replay = dict(record)
    replay.pop("record_sha256", None)
    actual = sha256(canonical_json(replay).encode()).hexdigest()
    if claimed != actual:
        raise RuntimeError("cache record hash mismatch")


def main() -> None:
    compact.flint.ctx.prec = collar.bridge.PRECISION_BITS
    issues: list[str] = []
    if not RESULT_PATH.exists():
        raise SystemExit("missing collar result")
    if not NOTE_PATH.exists():
        issues.append("missing collar note")
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    if artifact.get("kind") != collar.STEM:
        issues.append("kind mismatch")
    stored_contract = artifact.get("contract", {})
    current_contract = collar.contract_payload()
    normalized_contract = deepcopy(current_contract)
    stored_successor_hash = stored_contract.get(
        "source_sha256", {}
    ).get("successor_lemma")
    if stored_successor_hash not in {
        file_hash(collar.SUCCESSOR_RESULT),
        HISTORICAL_SUCCESSOR_SHA256,
    }:
        issues.append("successor source snapshot is unrecognized")
    normalized_contract["source_sha256"]["successor_lemma"] = (
        stored_successor_hash
    )
    if stored_contract != normalized_contract:
        issues.append("contract payload mismatch")
    snapshot_contract_hash = sha256(
        canonical_json(stored_contract).encode()
    ).hexdigest()
    if artifact.get("contract_sha256") != snapshot_contract_hash:
        issues.append("contract hash mismatch")
    for label, expected in collar.source_hashes().items():
        observed = stored_contract["source_sha256"].get(label)
        if label == "successor_lemma":
            continue
        if observed != expected:
            issues.append(f"source hash mismatch: {label}")
    successor = json.loads(
        collar.SUCCESSOR_RESULT.read_text(encoding="utf-8")
    )
    if successor.get("kind") != (
        "jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma"
    ):
        issues.append("current successor semantic source kind drifted")
    if successor.get("exact", {}).get("transport", {}).get(
        "strict_gate"
    ) != "delta_j*M_(j,k)<d_(j,k)":
        issues.append("current successor semantic transport gate drifted")

    tasks = collar.panel_tasks()
    references = collar.load_reference_leaves()
    records = artifact.get("records", [])
    previous_hash = "GENESIS"
    if len(records) > len(tasks):
        issues.append("too many cache records")
    for index, record in enumerate(records):
        task = tasks[index]
        try:
            validate_snapshot_cache_record(
                record,
                task,
                previous_hash,
                snapshot_contract_hash,
            )
        except RuntimeError as exc:
            issues.append(f"cache row {index + 1}: {exc}")
            break
        previous_hash = record["record_sha256"]
        result = record["result"]
        reference = references[index]
        if (
            result.get("x_low") != task["x_low"]
            or result.get("x_high") != task["x_high"]
        ):
            issues.append(f"panel {index + 1} x-domain mismatch")

        reference_f = compact.arb(reference["full_f"])
        reference_fp = compact.arb(reference["full_f_prime"])
        reference_distance = (
            reference_f.abs_lower() ** 2
            + reference_fp.abs_lower() ** 2
        ).sqrt().lower()
        stored_distance = compact.arb(
            result["reference_phase_cell"]["distance_lower"]
        )
        if not stored_distance.overlaps(reference_distance):
            issues.append(
                f"panel {index + 1} reference distance mismatch"
            )

        h_xx = compact.arb(result["full"]["h_xx"])
        h_xxx = compact.arb(result["full"]["h_xxx"])
        x_low = Fraction(task["x_low"])
        x_high = Fraction(task["x_high"])
        x_box = collar.bridge.interval_ball(x_low, x_high)
        f_t = -16 * (1 + x_box**4) * h_xx
        f_xt = (
            -64 * x_box**3 * h_xx
            - 16 * (1 + x_box**4) * h_xxx
        )
        stored_f_t = compact.arb(result["transport"]["f_t"])
        stored_f_xt = compact.arb(result["transport"]["f_xt"])
        if not stored_f_t.overlaps(f_t):
            issues.append(f"panel {index + 1} F_t formula mismatch")
        if not stored_f_xt.overlaps(f_xt):
            issues.append(f"panel {index + 1} F_xt formula mismatch")
        norm_upper = (
            stored_f_t.abs_upper().upper() ** 2
            + stored_f_xt.abs_upper().upper() ** 2
        ).sqrt().upper()
        displacement = (
            collar.fraction_to_arb(collar.TIME_STEP) * norm_upper
        ).upper()
        ratio = (displacement / reference_distance).upper()
        stored_transport = result["transport"]
        if not compact.arb(
            stored_transport["jet_time_derivative_norm_upper"]
        ).overlaps(norm_upper):
            issues.append(f"panel {index + 1} norm mismatch")
        if not compact.arb(
            stored_transport["displacement_upper"]
        ).overlaps(displacement):
            issues.append(f"panel {index + 1} displacement mismatch")
        if not compact.arb(
            stored_transport[
                "displacement_to_cell_distance_ratio_upper"
            ]
        ).overlaps(ratio):
            issues.append(f"panel {index + 1} ratio mismatch")
        certified = bool(ratio < 1)
        expected_status = "certified" if certified else "failed"
        if result.get("status") != expected_status:
            issues.append(f"panel {index + 1} status mismatch")

    summary = artifact.get("summary", {})
    certified_count = sum(
        record["result"]["status"] == "certified"
        for record in records
    )
    failed_count = len(records) - certified_count
    complete = (
        len(records) == len(tasks) and failed_count == 0
    )
    if summary.get("tasks_total") != len(tasks):
        issues.append("summary task total mismatch")
    if summary.get("tasks_completed") != len(records):
        issues.append("summary completed count mismatch")
    if summary.get("certified_panels") != certified_count:
        issues.append("summary certified count mismatch")
    if summary.get("failed_panels") != failed_count:
        issues.append("summary failed count mismatch")
    if summary.get("complete_collar_certificate") != complete:
        issues.append("summary completion mismatch")
    if complete and "complete" not in artifact.get("status", ""):
        issues.append("complete artifact status mismatch")
    if not complete and "partial" not in artifact.get("status", ""):
        issues.append("partial artifact status mismatch")

    boundary = artifact.get("proof_boundary", "")
    for phrase in (
        "one finite successor",
        "No Q208-to-Q209 transport bound",
        "no uniform all-j collar estimate",
        "no Lambda<=0",
        "no RH proof",
    ):
        if phrase not in boundary:
            issues.append(f"proof boundary lacks {phrase!r}")

    cache_path = REPO_ROOT / artifact["cache"]
    if artifact.get("cache_sha256") != file_hash(cache_path):
        issues.append("cache file hash mismatch")
    if artifact.get("cache_last_record_sha256") != previous_hash:
        issues.append("cache tail hash mismatch")

    if issues:
        raise SystemExit("\n".join(issues))
    print(
        "validated Q207-Q208 adiabatic bottom collar: "
        f"{len(records)}/{len(tasks)} panels, "
        f"{certified_count} certified, {failed_count} failed, "
        f"complete={complete}, 0 all-j promotions"
    )


if __name__ == "__main__":
    main()
