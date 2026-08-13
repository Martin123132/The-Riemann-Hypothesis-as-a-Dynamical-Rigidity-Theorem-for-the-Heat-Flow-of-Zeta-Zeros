#!/usr/bin/env python3
"""Independently check the R=9 connector quadrature gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
for candidate in (REPO_ROOT / "work/rh_compute/vendor", SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

from flint import arb, acb, ctx

from check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_canonical_core_quadrature_gate import IndependentCore


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_R9_connector_quadrature_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_R9_connector_quadrature_gate.md"
FORMULA_VERSION = "R9_connector_degree17_tanh_v1"
PANELS_PER_SIDE = 20


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def overlap(left: acb, right: acb) -> bool:
    return (left.real-right.real).contains(arb(0)) and (left.imag-right.imag).contains(arb(0))


def bounds(index: int) -> tuple[arb, arb]:
    left = arb(-9)+arb(index)/4 if index<PANELS_PER_SIDE else arb(4)+arb(index-PANELS_PER_SIDE)/4
    return left,left+arb(1)/4


def main() -> None:
    require(RESULT.is_file() and NOTE.is_file(), "missing result or note")
    artifact=json.loads(RESULT.read_text(encoding="utf-8")); require(artifact.get("passed") is True,"gate not passed")
    for record in artifact["dependencies"].values():
        path=REPO_ROOT/record["path"]; require(path.is_file() and file_hash(path)==record["sha256"],f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path=REPO_ROOT/record["path"]; require(path.is_file() and file_hash(path)==record["sha256"],f"source hash drift: {path}")
    cache_path=REPO_ROOT/artifact["dependencies"]["cache"]["path"]
    rows=[json.loads(line) for line in cache_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    rows=[r for r in rows if r.get("formula_version")==FORMULA_VERSION]
    require(len(rows)==80,"connector cache row count drift")
    by={(r["integral_kind"],int(r["panel_index"])):r for r in rows}; require(len(by)==80,"duplicate connector key")
    canonical=sum((parse_complex(by[("canonical",i)]["value"]) for i in range(40)),acb(0))
    finite=sum((parse_complex(by[("finite_t_degree17",i)]["value"]) for i in range(40)),acb(0))
    values=artifact["certified_values"]
    require(overlap(canonical,parse_complex(values["canonical_connector_ball"])),"canonical connector sum drift")
    require(overlap(finite,parse_complex(values["finite_t_connector_surrogate_ball"])),"finite connector sum drift")

    ctx.dps=40; core=IndependentCore()
    for kind,index,function in (("canonical",0,core.canonical),("finite_t_degree17",39,core.finite)):
        left,right=bounds(index)
        fresh=acb.integral(function,left,right,abs_tol=arb("1e-11"),rel_tol=arb("1e-11"),eval_limit=500000,depth_limit=50)
        require(overlap(fresh,parse_complex(by[(kind,index)]["value"])),f"fresh connector panel drift: {kind} {index}")
    tail=arb(values["true_R9_outgoing_tail_relative_ball"]); correction=arb(values["finite_t_R9_relative_correction_ball"])
    require(tail<arb("0.01") and correction<arb("0.01"),"R9 threshold failed")
    require(artifact["decision"]["true_R9_tail_direct_ray_bound_proved"] is False,"ray proof boundary drift")
    require(artifact["decision"]["y_gt_64_join_proved"] is False,"outer join boundary drift")
    note=NOTE.read_text(encoding="utf-8")
    for token in ("[-9,-4] union [4,9]","I_tail","not a proof"):
        require(token in note,f"note token missing: {token}")
    print("validated R9 connector independently: 80 cache rows, true tail and finite-t correction below one percent")


if __name__ == "__main__":
    main()
