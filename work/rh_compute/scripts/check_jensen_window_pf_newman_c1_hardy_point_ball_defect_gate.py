#!/usr/bin/env python3
"""Validate the exact-binary128 Arb point-ball Hardy defect gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_point_ball_defect_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_point_ball_defect_gate.md"
BUILDER = REPO_ROOT / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_point_ball_defect_gate.py"
CHECKER = Path(__file__).resolve()
ACCEPTED = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
)
ACCEPTED_SHA256 = "0fb64f090c21b185f27edc1c9194254cf471e74bfd6af3b88cb7f0cb40ec0c3d"
EXPECTED_TPM_HEX = "C001921FB54442D18469898CC51701B8"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tree_hash(root: Path) -> tuple[str, int, int]:
    digest = hashlib.sha256()
    count = 0
    size = 0
    for path in sorted(root.rglob("*")):
        if not path.is_file() or "__pycache__" in path.parts or path.suffix == ".pyc":
            continue
        relative_path = path.relative_to(root).as_posix().encode("utf-8")
        payload = path.read_bytes()
        digest.update(relative_path)
        digest.update(b"\0")
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
        count += 1
        size += len(payload)
    return digest.hexdigest(), count, size


def dyadic(value: list[int]) -> Fraction:
    require(isinstance(value, list) and len(value) == 2, "invalid dyadic record")
    mantissa, exponent = int(value[0]), int(value[1])
    return Fraction(mantissa * (2**exponent), 1) if exponent >= 0 else Fraction(mantissa, 2 ** (-exponent))


def lower(record: dict) -> Fraction:
    return dyadic(record["lower_dyadic"])


def upper(record: dict) -> Fraction:
    return dyadic(record["upper_dyadic"])


def main() -> int:
    for path in (RESULT, NOTE, BUILDER, CHECKER, ACCEPTED):
        require(path.is_file(), f"missing point-ball gate artifact: {path}")
    require(file_hash(ACCEPTED) == ACCEPTED_SHA256, "accepted evaluator hash drift")

    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_point_ball_defect_gate",
        "point-ball kind drift",
    )
    require(artifact["status"] == "exact_binary128_point_ball_certificate", "status drift")

    runtime = artifact["runtime"]
    require(runtime["python_flint_version"] == "0.8.0", "python-flint version drift")
    require(runtime["flint_version"] == "3.3.1", "FLINT version drift")
    require(runtime["flint_release"] == 30301, "FLINT release drift")
    require(runtime["threads"] == 1, "Arb thread count drift")
    require(runtime["precision_ladder_dps"] == [96, 160], "Arb precision ladder drift")
    runtime_root = Path(runtime["root"])
    require(runtime_root.is_dir(), "pinned python-flint runtime missing")
    runtime_hash, runtime_files, runtime_bytes = tree_hash(runtime_root)
    require(runtime_hash == runtime["tree_sha256"], "python-flint runtime tree hash drift")
    require(runtime_files == runtime["file_count"], "runtime file count drift")
    require(runtime_bytes == runtime["byte_count"], "runtime byte count drift")

    exactness = artifact["input_exactness"]
    require(exactness["decimal_strings_used_for_arithmetic"] is False, "decimal arithmetic reintroduced")
    require(exactness["tpm_hex"] == EXPECTED_TPM_HEX, "source tpm payload drift")
    require(exactness["source_phase_definition"] == "p=4*ATAN(1); tpp=2*p; tpm=-tpp", "pi provenance drift")
    require(dyadic(exactness["tpm_exact_dyadic"]) < 0, "decoded tpm sign drift")

    aggregate = artifact["aggregate"]
    require(aggregate["chain_count"] == 64, "point-ball chain count drift")
    require(aggregate["branch_counts"] == {"1": 32, "2": 32}, "branch count drift")
    require(aggregate["precision_overlap_count"] == 64, "precision overlap failure")
    require(aggregate["defect_zero_exclusion_count"] == 64, "zero exclusion failure")
    require(
        lower(aggregate["minimum_defect_abs_lower_bound"]) > Fraction(3, 100),
        "minimum rigorous defect bound lost",
    )
    require(
        upper(aggregate["maximum_source_roundoff_gap_upper_bound"]) < Fraction(1, 10**27),
        "binary128 roundoff gap unexpectedly large",
    )
    require(
        upper(aggregate["maximum_gap_to_defect_ratio_upper_bound"]) < Fraction(1, 10**20),
        "roundoff gap is no longer negligible relative to defect",
    )

    rows = artifact["rows"]
    require([row["chain"] for row in rows] == list(range(1, 65)), "point-ball row sequence gap")
    for row in rows:
        require(row["tpm_hex"] == EXPECTED_TPM_HEX, f"chain {row['chain']} tpm drift")
        require(row["precision_balls_overlap"] is True, f"chain {row['chain']} precision mismatch")
        require(row["defect_excludes_zero"] is True, f"chain {row['chain']} contains zero")
        require(row["defect"]["contains_zero"] is False, f"chain {row['chain']} defect flag drift")
        require(lower(row["defect"]["abs"]) > 0, f"chain {row['chain']} nonpositive defect lower bound")
        require(
            upper(row["high_precision_component_radius"])
            < lower(row["low_precision_component_radius"]),
            f"chain {row['chain']} radius did not shrink",
        )
        require(
            upper(row["source_roundoff_gap_to_defect_ratio"]) < Fraction(1, 10**20),
            f"chain {row['chain']} source roundoff gap too large",
        )

    finite_target = artifact["finite_target"]
    require("exp(i*tpm" in finite_target["raw_level_sum"], "finite phase formula drift")
    require("exact logged binary128 source qq" in finite_target["q_input"], "q scope drift")
    conclusion = artifact["route_conclusion"]
    require(conclusion["preferred_next_route"] == "interval_direct_local_defect", "route conclusion drift")
    require("not explained by binary128 operation rounding" in conclusion["finding"], "roundoff finding drift")

    sources = artifact["sources"]
    require(file_hash(Path(REPO_ROOT / sources["telemetry"]["path"])) == sources["telemetry"]["sha256"], "telemetry source hash drift")
    require(file_hash(Path(REPO_ROOT / sources["scout"]["path"])) == sources["scout"]["sha256"], "scout source hash drift")
    require(file_hash(BUILDER) == sources["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == sources["checker"]["sha256"], "checker hash drift")

    for token in (
        "Date: 2026-08-06",
        "Status: rigorous saved-point source-expression certificate; not a proof and no uniform analytic error theorem",
        "32-hex-digit payload",
        "no separate value of pi is inserted",
        "all `64` high-precision balls exclude zero",
        "does not prove that the source-convention transformed levels equal the paper's",
        "prize-level conclusion",
    ):
        require(token in note, f"point-ball note token missing: {token}")

    print(
        "validated Hardy point-ball defect gate: 64 exact binary128 chains, "
        "64 zero exclusions, 2 Arb precisions, 1 pinned runtime and 0 uniform bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
