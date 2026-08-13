#!/usr/bin/env python3
"""Validate the exact first-level MGS parity/orientation adapter gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_initial_level_adapter_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_initial_level_adapter_gate.md"
BUILDER = REPO_ROOT / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_initial_level_adapter_gate.py"
CHECKER = Path(__file__).resolve()
SOURCE = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
)
SOURCE_SHA256 = "0fb64f090c21b185f27edc1c9194254cf471e74bfd6af3b88cb7f0cb40ec0c3d"
EXPECTED_TPM_HEX = "C001921FB54442D18469898CC51701B8"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fraction(record: dict[str, int]) -> Fraction:
    value = Fraction(int(record["numerator"]), int(record["denominator"]))
    require(value.denominator > 0, "invalid rational record")
    return value


def dyadic(value: list[int]) -> Fraction:
    require(isinstance(value, list) and len(value) == 2, "invalid dyadic record")
    mantissa, exponent = int(value[0]), int(value[1])
    return Fraction(mantissa * (2**exponent), 1) if exponent >= 0 else Fraction(mantissa, 2 ** (-exponent))


def main() -> int:
    for path in (RESULT, NOTE, BUILDER, CHECKER, SOURCE):
        require(path.is_file(), f"missing initial-adapter artifact: {path}")
    require(file_hash(SOURCE) == SOURCE_SHA256, "accepted source hash drift")

    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_initial_level_adapter_gate",
        "initial-adapter kind drift",
    )
    require(
        artifact["status"] == "exact_mathematical_adapter_with_open_binary_phase_budget",
        "initial-adapter status drift",
    )
    require(artifact["source"]["sha256"] == SOURCE_SHA256, "result source hash drift")
    require(len(artifact["source"]["locations"]) == 7, "source location inventory drift")

    lemma = artifact["exact_lemma"]
    require("ix*k^2+sigma*k" in lemma["phase_integer"], "odd parity formula drift")
    require("k^2+k and k^2-k are even" in lemma["parity_reason"], "parity proof drift")
    require("mathematical 2*pi" in lemma["conclusion"], "mathematical periodicity scope drift")

    saved = artifact["saved_fixture"]
    require(saved["chain_count"] == 64, "saved chain count drift")
    require(saved["branch_counts"] == {"1": 32, "2": 32}, "saved branch count drift")
    require(saved["ix_histogram"] == {"0": 64}, "saved ix histogram drift")
    require(saved["coefficient_map_exact_count"] == 64, "coefficient map failure")
    require(saved["integer_phase_checks"] == 6720, "saved parity check count drift")
    require(saved["zero_integer_phase_shift_count"] == 64, "saved zero-shift count drift")
    require([row["chain"] for row in saved["rows"]] == list(range(1, 65)), "saved row sequence gap")
    for row in saved["rows"]:
        require(row["ix"] == 0, f"chain {row['chain']} ix drift")
        require(row["source_conjugates"] is True, f"chain {row['chain']} orientation drift")
        require(row["mathematical_adapter_exact"] is True, f"chain {row['chain']} adapter failure")
        require(row["binary_tpp_periodicity_bound_zero"] is True, f"chain {row['chain']} binary phase drift")
        require(fraction(row["maximum_integer_phase_shift"]) == 0, f"chain {row['chain']} nonzero phase shift")

    synthetic = artifact["synthetic_parity_coverage"]
    require(synthetic["fixture_count"] == 36, "synthetic fixture count drift")
    require(synthetic["ix_values"] == list(range(-4, 5)), "synthetic ix coverage drift")
    require(synthetic["a1_signs"] == [-1, 1], "synthetic a1 signs drift")
    require(synthetic["xr_signs"] == [-1, 1], "synthetic xr signs drift")
    require(synthetic["checked_indices_per_fixture"] == 65, "synthetic index count drift")
    require(len(synthetic["rows"]) == 36, "synthetic row count drift")
    require(all(fraction(row["maximum_integer_phase_shift"]).denominator == 1 for row in synthetic["rows"]), "synthetic phase nonintegrality")

    binary = artifact["binary_phase_normalization"]
    require(binary["tpm_hex"] == EXPECTED_TPM_HEX, "binary tpm payload drift")
    delta_lower = dyadic(binary["delta_abs_ball"]["lower_dyadic"])
    delta_upper = dyadic(binary["delta_abs_ball"]["upper_dyadic"])
    require(Fraction(1, 10**35) < delta_lower <= delta_upper < Fraction(1, 10**33), "binary phase delta scale drift")
    require("sum_(k=0)^N" in binary["generic_bound"], "generic phase bound drift")
    require("zero because ix=0" in binary["saved_bound"], "saved phase conclusion drift")
    require("Record or bound ix and N" in binary["physical_requirement"], "physical phase requirement drift")

    conclusion = artifact["route_conclusion"]
    require("mathematical first-level" in conclusion["closed"], "closed adapter scope drift")
    require("All 64 rows have ix=0" in conclusion["saved_source_case"], "saved source scope drift")
    require("Nonzero-ix physical branches" in conclusion["open"], "open phase budget drift")
    require("level-one-to-level-two" in conclusion["next_target"], "next adapter target drift")

    sources = artifact["sources"]
    for key in ("telemetry", "point_ball"):
        path = REPO_ROOT / sources[key]["path"]
        require(file_hash(path) == sources[key]["sha256"], f"{key} dependency hash drift")
    require(file_hash(BUILDER) == sources["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == sources["checker"]["sha256"], "checker hash drift")

    for token in (
        "Date: 2026-08-06",
        "Status: exact mathematical level-one adapter; not a proof and physical nonzero-ix binary phase budget open",
        "`k^2+k` and `k^2-k` are even",
        "`6720` integer-phase checks",
        "every saved row has `ix=0`",
        "mathematical periodicity alone is not a bit-level source bound",
        "does not prove the deeper transformed-child adapter",
        "prize-level conclusion",
    ):
        require(token in note, f"initial-adapter note token missing: {token}")

    print(
        "validated Hardy initial-level adapter gate: 64 exact source maps, "
        "6720 saved parity checks, 36 synthetic sign/parity fixtures and 1 open binary phase budget"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
