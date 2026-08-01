#!/usr/bin/env python3
"""Discharge the order-eleven curvature target and compose lambda=-100 entry."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULTS = REPO_ROOT / "work/rh_compute/results"
ORDER10_CURVATURE = (
    RESULTS / "jensen_window_pf_compound_order10_first_summand_curvature_certificate.json"
)
ORDER11_CURVATURE = (
    RESULTS / "jensen_window_pf_compound_order11_first_summand_curvature_certificate.json"
)
BRIDGE_TARGET = RESULTS / "jensen_window_pf_compound_order11_curvature_bridge_target.json"
FINITE_PREFIX = RESULTS / "jensen_window_pf_compound_order11_m100_prefix_certificate.json"
DEFAULT_OUT = RESULTS / "jensen_window_pf_compound_order11_m100_entry_certificate.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_compound_order11_m100_entry_certificate.md"
GENERATOR_PATH = (
    "work/rh_compute/scripts/"
    "jensen_window_pf_compound_order11_m100_entry_certificate.py"
)
CHECKER_PATH = (
    "work/rh_compute/scripts/"
    "check_jensen_window_pf_compound_order11_m100_entry_certificate.py"
)

ORDER10_CONTINUOUS = "z_1''(t)<=4200/t^2 for every real t>=1251"
ORDER11_CONTINUOUS = "y_1''(t)<=6000/t^2 for every real t>=1252"
FIRST_DISCRETE = (
    "y_1''(t)<=6000/t^2 => Y_k^(1)<=6000*[-log(1-1/k^2)]"
    "<6001/k^2, k>=1253"
)
FULL_TRANSFER = "|Y_k-Y_k^(1)|<37/k^2 for every integer k>=1253"
FULL_CEILING = (
    "[z_1''(t)<=4200/t^2 on t>=1251 and y_1''(t)<=6000/t^2 on "
    "t>=1252] => Y_k<6001/k^2+37/k^2=6038/k^2<6100/k^2, k>=1253"
)
CONDITIONAL_TAIL = (
    "[z_1''(t)<=4200/t^2 on t>=1251 and y_1''(t)<=6000/t^2 on "
    "t>=1252] => Q_(11,n)(-100)>0 for every n>=1243"
)
TAIL_THEOREM = "Q_(11,n)(-100)>0 for every n>=1243"
FINITE_THEOREM = "Q_(11,n)(-100)>0 for every 0<=n<=1242"
GLOBAL_ENDPOINT = "Q_(11,n)(-100)>0 for every integer n>=0"


@dataclass(frozen=True)
class EntryRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def source_record(path: Path, artifact: dict) -> dict:
    return {
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "sha256": sha256(path),
        "kind": artifact["kind"],
        "status": artifact["status"],
    }


def validate_sources() -> dict:
    order10 = load_json(ORDER10_CURVATURE)
    order11 = load_json(ORDER11_CURVATURE)
    bridge = load_json(BRIDGE_TARGET)
    finite = load_json(FINITE_PREFIX)
    if (
        order10.get("status")
        != "rigorous global order-ten first-summand curvature theorem on t>=1251"
        or order10.get("theorem") != ORDER10_CONTINUOUS
        or order10.get("summary", {}).get(
            "global_first_summand_curvature_theorems"
        )
        != 1
        or order10.get("summary", {}).get("open_rows") != 0
    ):
        raise RuntimeError("global order-ten curvature source is not closed")
    if (
        order11.get("kind")
        != "jensen_window_pf_compound_order11_first_summand_curvature_certificate"
        or order11.get("status")
        != (
            "rigorous global order-eleven first-summand curvature theorem "
            "on t>=1252"
        )
        or order11.get("theorem") != ORDER11_CONTINUOUS
        or order11.get("summary", {}).get(
            "global_first_summand_curvature_theorems"
        )
        != 1
        or order11.get("summary", {}).get("open_rows") != 0
        or order11.get("summary", {}).get("full_kernel_theorems") != 0
        or order11.get("summary", {}).get("heat_forward_theorems") != 0
        or order11.get("summary", {}).get("rh_claims") != 0
    ):
        raise RuntimeError("global order-eleven curvature source is not closed")
    if (
        bridge.get("kind")
        != "jensen_window_pf_compound_order11_curvature_bridge_target"
        or bridge.get("status")
        != "exact conditional order-eleven curvature bridge with one open continuum target"
        or bridge.get("summary", {}).get("open_continuum_targets") != 1
        or bridge.get("summary", {}).get("open_rows") != 1
    ):
        raise RuntimeError("order-eleven conditional bridge identity changed")
    exact = bridge.get("exact", {})
    expected = {
        "continuous_target": ORDER11_CONTINUOUS,
        "tent_transfer": FIRST_DISCRETE,
        "full_transfer": FULL_TRANSFER,
        "conditional_full_ceiling": FULL_CEILING,
        "conditional_endpoint_tail": CONDITIONAL_TAIL,
        "finite_prefix_theorem": FINITE_THEOREM,
    }
    for key, value in expected.items():
        if exact.get(key) != value:
            raise RuntimeError(f"order-eleven bridge field changed: {key}")
    open_rows = [
        row
        for row in bridge.get("rows", [])
        if row.get("readiness") != "ready_to_apply"
    ]
    if (
        len(open_rows) != 1
        or open_rows[0].get("id") != "co11cbt_05_first_target"
        or open_rows[0].get("formula") != ORDER11_CONTINUOUS
    ):
        raise RuntimeError("conditional bridge has an unexpected open obligation")
    finite_exact = finite.get("exact", {})
    if (
        finite.get("kind")
        != "jensen_window_pf_compound_order11_m100_prefix_certificate"
        or finite.get("status")
        != "rigorous lambda=-100 signed order-eleven prefix through n=1242"
        or finite_exact.get("finite_prefix") != FINITE_THEOREM
        or finite.get("summary", {}).get("finite_prefix_theorems") != 1
        or finite.get("summary", {}).get("positive_Q11_rows") != 1243
        or finite.get("summary", {}).get("negative_Q11_rows") != 0
    ):
        raise RuntimeError("finite order-eleven endpoint source is not closed")
    return {
        "sources": [
            source_record(path, artifact)
            for path, artifact in (
                (ORDER10_CURVATURE, order10),
                (ORDER11_CURVATURE, order11),
                (BRIDGE_TARGET, bridge),
                (FINITE_PREFIX, finite),
            )
        ],
        "discharged_target": open_rows[0]["id"],
        "transfer_scaled": exact["power_envelope"]["transfer_scaled_decimal"],
    }


def build_artifact() -> dict:
    contract = validate_sources()
    rows = [
        EntryRow(
            "co11m100ec_01_order10_curvature",
            "theorem_input",
            "ready_to_apply",
            "The inherited order-ten first-summand theorem supplies the common X-floor premise.",
            ORDER10_CONTINUOUS,
            "First Newman summand at lambda=-100.",
        ),
        EntryRow(
            "co11m100ec_02_order11_curvature",
            "theorem_input",
            "ready_to_apply",
            "The global order-eleven first-summand theorem discharges the bridge's only open target.",
            ORDER11_CONTINUOUS,
            "Continuous first-summand theorem at lambda=-100.",
        ),
        EntryRow(
            "co11m100ec_03_discrete_transfer",
            "exact_theorem_application",
            "ready_to_apply",
            "The tent identity and rational perturbation envelope give the full discrete curvature ceiling.",
            FIRST_DISCRETE + "; " + FULL_TRANSFER + "; " + FULL_CEILING,
            "Full Newman kernel at lambda=-100 for integer k>=1253.",
        ),
        EntryRow(
            "co11m100ec_04_circularity_guard",
            "circularity_guard",
            "ready_to_apply",
            "The discharged continuum premise contains no full-kernel, heat-forward, or RH theorem.",
            "source full_kernel_theorems=0; heat_forward_theorems=0; rh_claims=0",
            "Prevents the endpoint conclusion from being assumed inside its first-summand premise.",
        ),
        EntryRow(
            "co11m100ec_05_analytic_tail",
            "analytic_theorem",
            "ready_to_apply",
            "The full curvature ceiling proves every order-eleven endpoint sign in the analytic tail.",
            TAIL_THEOREM,
            "Endpoint shifts n>=1243.",
        ),
        EntryRow(
            "co11m100ec_06_finite_prefix",
            "interval_theorem",
            "ready_to_apply",
            "The direct Arb and positive-cone certificate supplies every preceding endpoint shift.",
            FINITE_THEOREM,
            "Endpoint shifts 0<=n<=1242.",
        ),
        EntryRow(
            "co11m100ec_07_global_endpoint",
            "exact_theorem_composition",
            "ready_to_apply",
            "The finite prefix and analytic tail meet without an index gap.",
            GLOBAL_ENDPOINT,
            "Fixed order eleven and lambda=-100 only.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order11_m100_entry_certificate",
        "date": "2026-07-22",
        "status": "rigorous all-shift signed order-eleven entry at lambda=-100",
        "proof_boundary": (
            "This artifact proves fixed order-eleven contiguous signed Hankel "
            "positivity at lambda=-100. It does not prove the heat interval, "
            "all orders, PF-infinity, Lambda<=0, or RH."
        ),
        "source_contract": contract,
        "exact": {
            "order10_continuous": ORDER10_CONTINUOUS,
            "order11_continuous": ORDER11_CONTINUOUS,
            "first_discrete_ceiling": FIRST_DISCRETE,
            "full_kernel_transfer": FULL_TRANSFER,
            "full_curvature_ceiling": FULL_CEILING,
            "analytic_tail": TAIL_THEOREM,
            "finite_prefix": FINITE_THEOREM,
            "global_endpoint": GLOBAL_ENDPOINT,
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": len(rows),
            "open_rows": 0,
            "discharged_continuum_targets": 1,
            "circularity_guards": 1,
            "full_kernel_transfer_applications": 1,
            "analytic_tail_theorems": 1,
            "finite_prefix_theorems": 1,
            "global_m100_order11_entry_theorems": 1,
            "heat_interval_theorems": 0,
            "orders_above_11": 0,
            "rh_claims": 0,
        },
        "generator": GENERATOR_PATH,
        "checker": CHECKER_PATH,
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Order-Eleven Lambda=-100 Entry Certificate",
        "",
        "Date: 2026-07-22",
        "",
        f"Status: **{artifact['status']}**. This is not a proof of RH or `Lambda <= 0`.",
        "",
        "## Composition",
        "",
        f"`{exact['order10_continuous']}`",
        f"`{exact['order11_continuous']}`",
        f"`{exact['analytic_tail']}`",
        f"`{exact['finite_prefix']}`",
        "",
        f"Together these prove `{exact['global_endpoint']}`.",
        "",
        "## Boundary",
        "",
        artifact["proof_boundary"],
        "",
        "This is not a proof of RH.",
        "",
        "## Reproduce",
        "",
        "```powershell",
        f"python {GENERATOR_PATH}",
        f"python {CHECKER_PATH}",
        "```",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    print(f"wrote order-eleven endpoint theorem: {GLOBAL_ENDPOINT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
