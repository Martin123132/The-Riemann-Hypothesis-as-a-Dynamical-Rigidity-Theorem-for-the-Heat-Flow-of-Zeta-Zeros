#!/usr/bin/env python3
"""Validate the logarithmic endpoint-ray reduction gate."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_reduction_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_reduction_gate.md"
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_reduction_gate.py"
)
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    for path in (RESULT, NOTE, BUILDER, CHECKER):
        require(path.is_file(), f"missing log-ray artifact: {path}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_reduction_gate",
        "log-ray kind drift",
    )
    require(
        artifact["status"] == "exact_logarithmic_endpoint_ray_reduction_with_nonrigorous_orientation_diagnostics",
        "log-ray status drift",
    )
    scope = artifact["scope"]
    aggregate = artifact["aggregate"]
    rows = artifact["rows"]
    diagnostics = artifact["diagnostic_witnesses"]
    require(scope["recursive_call_count"] == aggregate["call_count"] == len(rows) == 374, "log-ray row count drift")
    require(scope["diagnostic_witness_count"] == len(diagnostics) == 4, "log-ray diagnostic count drift")
    require(scope["diagnostic_decimal_digits"] >= 40, "log-ray diagnostic precision drift")
    require(aggregate["length_one_count"] + aggregate["length_two_count"] == 374, "log-ray length count drift")
    require(len(aggregate["sign_combination_counts"]) == 4, "log-ray sign branch count drift")
    require(all(count > 0 for count in aggregate["sign_combination_counts"].values()), "log-ray empty sign branch")
    require("-Log(1-exp(2*pi*i*u))" in artifact["kernel_identity"], "log-ray kernel identity drift")
    require("P^-=I69^-" in artifact["source_oriented_identity"]["complete"], "log-ray complete identity drift")
    require("integrable logarithmic" in artifact["analytic_contract"]["origin"], "log-ray origin contract drift")

    chain_roster = [row["chain"] for row in rows]
    require(chain_roster == sorted(chain_roster) and len(set(chain_roster)) == 374, "log-ray row order drift")
    for row in rows:
        require(row["transformed_length"] in (1, 2), "log-ray row length drift")
        require(row["first_tail_index"] == row["transformed_length"] + 1, "log-ray row tail index drift")
        require(row["first_saddle_index"] in (0, 1), "log-ray row saddle index drift")
        require(row["b_angle"] in ("pi/2", "5*pi/6"), "log-ray row b angle drift")
        require(row["c_angle"] in ("pi/2", "pi/6"), "log-ray row c angle drift")
    diagnostic_signs = {(row["phi1_sign"], row["phi3_sign"]) for row in diagnostics}
    require(len(diagnostic_signs) == 4, "log-ray diagnostic sign coverage drift")
    for row in diagnostics:
        require(not row["rigorous"], "log-ray diagnostic mislabeled rigorous")
        require(Decimal(row["absolute_gap"]) < Decimal("1e-12"), "log-ray diagnostic gap drift")

    require(artifact["paper"]["equations"] == [67, 68, 69, 82, 83, 88, 90, 92, 95, 96], "log-ray paper equation drift")
    for group in (artifact["paper"], *artifact["sources"].values()):
        path = REPO_ROOT / group["path"]
        require(path.is_file(), f"log-ray dependency missing: {group['path']}")
        require(file_hash(path) == group["sha256"], f"log-ray dependency hash drift: {group['path']}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "log-ray builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "log-ray checker hash drift")
    for token in (
        "Status: exact analytic reduction with non-rigorous branch diagnostics",
        "K_m(u)",
        "P^-=I69^-+Q_exact^-",
        "four endpoint-ray",
        "not a proof",
    ):
        require(token in note, f"log-ray note token missing: {token}")
    require("four mpmath checks are non-rigorous" in artifact["proof_boundary"], "log-ray boundary drift")
    print("validated logarithmic endpoint-ray reduction gate: 374 calls, 4 diagnostic branches")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
