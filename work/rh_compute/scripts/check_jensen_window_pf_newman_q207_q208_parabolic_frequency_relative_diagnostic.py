#!/usr/bin/env python3
"""Check the Q207-Q208 parabolic-frequency relative diagnostic."""

from __future__ import annotations

import importlib.util
import json
from hashlib import sha256
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402
from flint import arb  # noqa: E402


STEM = "jensen_window_pf_newman_q207_q208_parabolic_frequency_relative_diagnostic"
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
BUILDER = REPO_ROOT / "work/rh_compute/scripts" / f"{STEM}.py"


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_builder():
    spec = importlib.util.spec_from_file_location(
        "relative_scale_builder", BUILDER
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load builder")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    if not RESULT.is_file() or not NOTE.is_file() or not BUILDER.is_file():
        raise FileNotFoundError("builder output, note, or builder is missing")
    flint.ctx.prec = 256
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    builder = load_builder()

    expected_status = (
        "rigorous finite Q207-Q208 four-scale relative-transport "
        "diagnostic with zero all-j promotions"
    )
    if artifact.get("kind") != STEM:
        raise RuntimeError("artifact kind mismatch")
    if artifact.get("date") != "2026-07-26":
        raise RuntimeError("artifact date mismatch")
    if artifact.get("status") != expected_status:
        raise RuntimeError("artifact status mismatch")
    if artifact.get("builder_sha256") != file_hash(BUILDER):
        raise RuntimeError("builder hash mismatch")

    contract = artifact.get("contract", {})
    if contract.get("precision_bits") != 256:
        raise RuntimeError("precision contract drifted")
    if contract.get("old_time") != "1/1035":
        raise RuntimeError("old time drifted")
    if contract.get("new_time") != "1/1040":
        raise RuntimeError("new time drifted")
    if contract.get("time_step") != "1/215280":
        raise RuntimeError("time step drifted")
    if contract.get("cover") != ["38", "245"]:
        raise RuntimeError("cover drifted")
    if contract.get("scales") != [
        "unit",
        "frequency",
        "parabolic",
        "blended",
    ]:
        raise RuntimeError("scale ordering drifted")
    if contract.get("source_sha256") != builder.source_hashes():
        raise RuntimeError("source hashes drifted")

    records = artifact.get("records", [])
    if len(records) != 525:
        raise RuntimeError(f"expected 525 records, found {len(records)}")
    if records[0]["x_low"] != "38" or records[-1]["x_high"] != "245":
        raise RuntimeError("record cover drifted")
    previous = "38"
    for sequence, row in enumerate(records, start=1):
        if row.get("sequence") != sequence:
            raise RuntimeError(f"sequence drift at {sequence}")
        if row.get("x_low") != previous:
            raise RuntimeError(f"cover gap before row {sequence}")
        previous = row["x_high"]
        costs = row.get("scale_costs", {})
        if list(costs) != [
            "blended",
            "frequency",
            "parabolic",
            "unit",
        ]:
            raise RuntimeError(f"stored scale map drift at row {sequence}")
        for name, values in costs.items():
            denominator = arb(values["denominator_lower"]).lower()
            derivative = arb(values["time_derivative_upper"]).upper()
            kappa = arb(values["kappa_upper"]).upper()
            cost = arb(values["one_step_log_cost_upper"]).upper()
            if denominator <= 0:
                raise RuntimeError(
                    f"nonpositive {name} denominator at row {sequence}"
                )
            if derivative < 0 or kappa < 0 or cost < 0:
                raise RuntimeError(
                    f"negative {name} upper bound at row {sequence}"
                )
    if previous != "245":
        raise RuntimeError("record cover does not terminate at 245")

    summary = artifact.get("summary", {})
    if summary.get("panels") != 525:
        raise RuntimeError("panel summary drifted")
    if summary.get("rigorous_finite_relative_certificates") != 4:
        raise RuntimeError("scale-certificate count drifted")
    if not summary.get("all_scales_all_panels_finite"):
        raise RuntimeError("finite-scale summary failed")
    if summary.get("all_j_promotions") != 0:
        raise RuntimeError("finite diagnostic was promoted to all-j")
    for name in ("unit", "frequency", "parabolic", "blended"):
        scale = summary.get("scales", {}).get(name, {})
        if not scale.get("all_525_denominators_positive"):
            raise RuntimeError(f"{name} denominator summary failed")
        if not scale.get("finite_relative_certificate"):
            raise RuntimeError(f"{name} relative summary failed")

    rebuilt_records, rebuilt_summary = builder.build_diagnostics()
    if rebuilt_summary != summary:
        raise RuntimeError("independently rebuilt summary drifted")
    if rebuilt_records[0] != records[0]:
        raise RuntimeError("first rebuilt record drifted")
    if rebuilt_records[262] != records[262]:
        raise RuntimeError("middle rebuilt record drifted")
    if rebuilt_records[-1] != records[-1]:
        raise RuntimeError("last rebuilt record drifted")

    note = NOTE.read_text(encoding="utf-8")
    required = (
        "# Q207-Q208 Parabolic-Frequency Relative Diagnostic",
        "It performs no new Xi quadrature.",
        "Every denominator in (1) is strictly positive.",
        "A finite one-step relative certificate is",
        "not an asymptotic Xi bound",
        "They prove no Q209 stage",
    )
    missing = [text for text in required if text not in note]
    if missing:
        raise RuntimeError(f"note contract missing: {missing}")

    print(
        "validated Q207-Q208 parabolic-frequency relative diagnostic: "
        "525 panels, 4 rigorous finite scale certificates, "
        "0 all-j promotions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
