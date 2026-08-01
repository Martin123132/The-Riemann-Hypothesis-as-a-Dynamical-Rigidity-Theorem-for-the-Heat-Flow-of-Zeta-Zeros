#!/usr/bin/env python3
"""Check the hybrid refined Q207-Q208 bottom-collar certificate."""

from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import jensen_window_pf_newman_q207_q208_adiabatic_bottom_collar_refined_tail_certificate as refined


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT_PATH = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_q207_q208_adiabatic_"
    "bottom_collar_refined_tail_certificate.json"
)
NOTE_PATH = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_newman_q207_q208_adiabatic_"
    "bottom_collar_refined_tail_certificate.md"
)
compact = refined.compact
coarse = refined.coarse
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
    schema: str,
    contract_sha256: str,
) -> None:
    if record.get("schema") != schema:
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


def normalize_successor_snapshot(
    stored_contract: dict,
    current_contract: dict,
    issues: list[str],
    context: str,
) -> tuple[dict, str]:
    normalized = deepcopy(current_contract)
    stored_successor_hash = stored_contract.get(
        "source_sha256", {}
    ).get("successor_lemma")
    if stored_successor_hash not in {
        file_hash(coarse.SUCCESSOR_RESULT),
        HISTORICAL_SUCCESSOR_SHA256,
    }:
        issues.append(f"{context} successor source snapshot is unrecognized")
    normalized["source_sha256"]["successor_lemma"] = (
        stored_successor_hash
    )
    snapshot_hash = sha256(
        canonical_json(stored_contract).encode()
    ).hexdigest()
    return normalized, snapshot_hash


def load_coarse_prefix_snapshot(issues: list[str]) -> list[dict]:
    artifact = json.loads(
        refined.COARSE_RESULT.read_text(encoding="utf-8")
    )
    stored_contract = artifact.get("contract", {})
    normalized, snapshot_hash = normalize_successor_snapshot(
        stored_contract,
        coarse.contract_payload(),
        issues,
        "coarse prefix",
    )
    if stored_contract != normalized:
        issues.append("coarse prefix contract payload mismatch")
    if artifact.get("contract_sha256") != snapshot_hash:
        issues.append("coarse prefix contract hash mismatch")
    for label, expected in coarse.source_hashes().items():
        if label == "successor_lemma":
            continue
        if stored_contract.get("source_sha256", {}).get(label) != expected:
            issues.append(f"coarse prefix source hash mismatch: {label}")

    records = artifact.get("records", [])
    prefix = records[:379]
    if len(prefix) != 379:
        issues.append("coarse prefix is incomplete")
        return prefix
    tasks = coarse.panel_tasks()[:379]
    previous_hash = "GENESIS"
    for index, (record, task) in enumerate(
        zip(prefix, tasks, strict=True),
        start=1,
    ):
        try:
            validate_snapshot_cache_record(
                record,
                task,
                previous_hash,
                coarse.SCHEMA,
                snapshot_hash,
            )
        except RuntimeError as exc:
            issues.append(f"coarse prefix row {index}: {exc}")
            break
        previous_hash = record["record_sha256"]
        if record["result"]["status"] != "certified":
            issues.append(f"coarse prefix row {index} is not certified")
    if Fraction(prefix[-1]["task"]["x_high"]) != refined.PREFIX_END:
        issues.append("coarse prefix endpoint mismatch")
    return prefix


def recompute_ratio(result: dict, task: dict, reference: dict):
    reference_f = compact.arb(reference["full_f"])
    reference_fp = compact.arb(reference["full_f_prime"])
    distance = (
        reference_f.abs_lower() ** 2
        + reference_fp.abs_lower() ** 2
    ).sqrt().lower()
    stored_f_t = compact.arb(result["transport"]["f_t"])
    stored_f_xt = compact.arb(result["transport"]["f_xt"])
    norm = (
        stored_f_t.abs_upper().upper() ** 2
        + stored_f_xt.abs_upper().upper() ** 2
    ).sqrt().upper()
    displacement = (
        coarse.fraction_to_arb(coarse.TIME_STEP) * norm
    ).upper()
    ratio = (displacement / distance).upper()

    h_xx = compact.arb(result["full"]["h_xx"])
    h_xxx = compact.arb(result["full"]["h_xxx"])
    x_box = coarse.bridge.interval_ball(
        Fraction(task["x_low"]),
        Fraction(task["x_high"]),
    )
    formula_f_t = -16 * (1 + x_box**4) * h_xx
    formula_f_xt = (
        -64 * x_box**3 * h_xx
        - 16 * (1 + x_box**4) * h_xxx
    )
    if not stored_f_t.overlaps(formula_f_t):
        raise AssertionError("F_t formula mismatch")
    if not stored_f_xt.overlaps(formula_f_xt):
        raise AssertionError("F_xt formula mismatch")
    return ratio


def main() -> None:
    compact.flint.ctx.prec = coarse.bridge.PRECISION_BITS
    issues: list[str] = []
    if not RESULT_PATH.exists():
        raise SystemExit("missing refined collar result")
    if not NOTE_PATH.exists():
        issues.append("missing refined collar note")
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    if artifact.get("kind") != refined.STEM:
        issues.append("kind mismatch")
    stored_contract = artifact.get("contract", {})
    normalized_contract, snapshot_contract_hash = (
        normalize_successor_snapshot(
            stored_contract,
            refined.contract_payload(),
            issues,
            "refined collar",
        )
    )
    if stored_contract != normalized_contract:
        issues.append("contract payload mismatch")
    if artifact.get("contract_sha256") != snapshot_contract_hash:
        issues.append("contract hash mismatch")
    for label, expected in refined.source_hashes().items():
        if label == "successor_lemma":
            continue
        if stored_contract.get("source_sha256", {}).get(label) != expected:
            issues.append(f"source hash mismatch: {label}")
    successor = json.loads(
        coarse.SUCCESSOR_RESULT.read_text(encoding="utf-8")
    )
    if successor.get("kind") != (
        "jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma"
    ):
        issues.append("current successor semantic source kind drifted")
    if successor.get("exact", {}).get("transport", {}).get(
        "strict_gate"
    ) != "delta_j*M_(j,k)<d_(j,k)":
        issues.append("current successor semantic transport gate drifted")

    prefix = load_coarse_prefix_snapshot(issues)
    if len(prefix) != 379:
        issues.append("coarse prefix count mismatch")
    prefix_ratios = [refined.ratio_ball(row) for row in prefix]
    if any(not bool(ratio < 1) for ratio in prefix_ratios):
        issues.append("coarse prefix contains a failed ratio")

    tasks = refined.refined_tasks()
    references = refined.parent_references()
    records = artifact.get("refined_records", [])
    previous_hash = "GENESIS"
    checked_ratios = list(prefix_ratios)
    if len(records) > len(tasks):
        issues.append("too many refined cache records")
    for index, record in enumerate(records):
        if index >= len(tasks):
            break
        task = tasks[index]
        try:
            validate_snapshot_cache_record(
                record,
                task,
                previous_hash,
                refined.SCHEMA,
                snapshot_contract_hash,
            )
        except RuntimeError as exc:
            issues.append(f"refined row {index + 1}: {exc}")
            break
        previous_hash = record["record_sha256"]
        try:
            ratio = recompute_ratio(
                record["result"],
                task,
                references[index],
            )
        except AssertionError as exc:
            issues.append(f"refined row {index + 1}: {exc}")
            continue
        stored_ratio = compact.arb(
            record["result"]["transport"][
                "displacement_to_cell_distance_ratio_upper"
            ]
        )
        if not stored_ratio.overlaps(ratio):
            issues.append(f"refined row {index + 1}: ratio mismatch")
        expected_status = (
            "certified" if bool(ratio < 1) else "failed"
        )
        if record["result"]["status"] != expected_status:
            issues.append(f"refined row {index + 1}: status mismatch")
        checked_ratios.append(stored_ratio)

    summary = artifact.get("summary", {})
    failed = sum(
        record["result"]["status"] != "certified"
        for record in records
    )
    complete = (
        len(prefix) == 379
        and len(records) == len(tasks)
        and failed == 0
    )
    expected_summary = {
        "coarse_prefix_panels": len(prefix),
        "refined_tasks_total": len(tasks),
        "refined_tasks_completed": len(records),
        "refined_certified_panels": len(records) - failed,
        "refined_failed_panels": failed,
        "combined_panel_count": len(prefix) + len(records),
        "complete_collar_certificate": complete,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(f"summary mismatch: {key}")
    if records:
        observed_max = max(checked_ratios)
        claimed_max = compact.arb(
            summary["maximum_transport_ratio_upper"]
        )
        if not claimed_max.overlaps(observed_max):
            issues.append("summary maximum ratio mismatch")

    cache_path = REPO_ROOT / artifact["cache"]
    if artifact.get("cache_sha256") != file_hash(cache_path):
        issues.append("refined cache hash mismatch")
    if artifact.get("cache_last_record_sha256") != previous_hash:
        issues.append("refined cache tail mismatch")
    boundary = artifact.get("proof_boundary", "")
    for phrase in (
        "one finite successor",
        "No Q208-to-Q209 transport budget",
        "no uniform all-j collar estimate",
        "no Lambda<=0",
        "no RH proof",
    ):
        if phrase not in boundary:
            issues.append(f"proof boundary lacks {phrase!r}")

    if issues:
        raise SystemExit("\n".join(issues))
    print(
        "validated hybrid Q207-Q208 adiabatic collar: "
        f"379 coarse + {len(records)}/{len(tasks)} refined panels, "
        f"{failed} failures, complete={complete}, "
        "0 all-j promotions"
    )


if __name__ == "__main__":
    main()
