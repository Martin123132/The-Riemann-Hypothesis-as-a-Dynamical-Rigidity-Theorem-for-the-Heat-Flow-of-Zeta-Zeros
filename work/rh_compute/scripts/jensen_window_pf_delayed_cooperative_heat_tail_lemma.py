#!/usr/bin/env python3
"""Record the shifted cooperative heat-flow descent lemma."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULTS = REPO_ROOT / "work/rh_compute/results"
HEAT_SOURCE = RESULTS / "jensen_window_pf_all_order_endpoint_heat_reduction.json"
ORDER10_SOURCE = RESULTS / "jensen_window_pf_compound_order10_m100_tail_curvature_reduction.json"
DEFAULT_OUT = RESULTS / "jensen_window_pf_delayed_cooperative_heat_tail_lemma.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_delayed_cooperative_heat_tail_lemma.md"

FIXED_TAIL = (
    "for every fixed m exists N_m such that for every n>=N_m and "
    "-100<=lambda<=0, Q_(m,n)(lambda)>0"
)
COOPERATIVE_FLOW = (
    "Q_(m,n)'=a_(m,n)*Q_(m,n+1)+b_(m,n)*Q_(m,n), "
    "a_(m,n)=c_(m,n)*Q_(m-1,n)/Q_(m-1,n+1)>0, "
    "b_(m,n)=c_(m,n)/(c_(m,n)-4)*(log Q_(m-1,n+1))'"
)
VARIATION = (
    "Q_(m,n)(lambda)=E_(m,n)(lambda)*(Q_(m,n)(-100)+"
    "integral_(-100)^lambda E_(m,n)(s)^(-1)*a_(m,n)(s)*"
    "Q_(m,n+1)(s)ds), E_(m,n)>0"
)
SHIFTED_LEMMA = (
    "[Q_(m-1,n)(lambda)>0 for every n>=n0 on -100<=lambda<=0, "
    "the fixed-order m eventual tail holds, and Q_(m,n)(-100)>0 "
    "for every n>=n0] => [Q_(m,n)(lambda)>0 for every n>=n0 "
    "and -100<=lambda<=0]"
)
ORDER10_HANDOFF = (
    "[Q_(10,n)(-100)>0 for every integer n>=4] implies "
    "Q_(10,n)(lambda)>0 for every n>=4 and -100<=lambda<=0"
)
ORDER11_HANDOFF = (
    "[Q_(10,n)(lambda)>0 for every n>=4 and -100<=lambda<=0 and "
    "Q_(11,n)(-100)>0 for every n>=4] implies "
    "Q_(11,n)(lambda)>0 for every n>=4 and -100<=lambda<=0"
)


@dataclass(frozen=True)
class LemmaRow:
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


def source_record(path: Path, payload: dict) -> dict:
    return {
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "sha256": sha256(path),
        "kind": payload.get("kind"),
        "status": payload.get("status"),
    }


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def validate_sources() -> list[dict]:
    heat = load_json(HEAT_SOURCE)
    order10 = load_json(ORDER10_SOURCE)
    exact = heat.get("exact", {})
    expected = {
        "fixed_order_tail": FIXED_TAIL,
        "cooperative_flow": COOPERATIVE_FLOW,
        "variation_of_constants": VARIATION,
        "known_base": (
            "Q_(m,n)(lambda)>0 for 1<=m<=9, every n>=0, and "
            "-100<=lambda<=0"
        ),
    }
    for key, value in expected.items():
        if exact.get(key) != value:
            raise RuntimeError(f"heat-flow source contract changed at {key}")
    summary = heat.get("summary", {})
    if summary.get("all_fixed_order_tail_theorems") != 1:
        raise RuntimeError("fixed-order eventual-tail theorem is not closed")
    if summary.get("all_order_cooperative_recursions") != 1:
        raise RuntimeError("cooperative heat recursion is not closed")
    if summary.get("completed_base_order") != 9:
        raise RuntimeError("completed lower heat base changed")
    if order10.get("exact", {}).get("conditional_heat_handoff") != ORDER10_HANDOFF:
        raise RuntimeError("order-ten delayed handoff contract changed")
    return [source_record(HEAT_SOURCE, heat), source_record(ORDER10_SOURCE, order10)]


def build_artifact() -> dict:
    sources = validate_sources()
    rows = [
        LemmaRow(
            "dcht_01_eventual_tail",
            "theorem_input",
            "ready_to_apply",
            "Every fixed signed Hankel order has a uniform positive heat tail.",
            FIXED_TAIL,
            "The threshold N_m may depend on m but not on lambda.",
        ),
        LemmaRow(
            "dcht_02_cooperative_equation",
            "exact_identity",
            "ready_to_apply",
            "A positive lower layer makes the next-shift coupling strictly positive.",
            COOPERATIVE_FLOW,
            "Requires Q_(m-1,n)>0 at the relevant shifts throughout the interval.",
        ),
        LemmaRow(
            "dcht_03_scalar_solution",
            "exact_ode_solution",
            "ready_to_apply",
            "Variation of constants preserves strict positivity when the next shift is positive.",
            VARIATION,
            "E is the positive integrating factor for the scalar equation at fixed n.",
        ),
        LemmaRow(
            "dcht_04_finite_descent",
            "induction_lemma",
            "ready_to_apply",
            "Choose N=max(n0,N_m) and descend from n=N-1 to n=n0; no smaller shift is used.",
            "tail n>=N; induction order N-1,N-2,...,n0",
            "The descent has finitely many steps even though the shift chain is infinite.",
        ),
        LemmaRow(
            "dcht_05_shifted_lemma",
            "theorem",
            "ready_to_apply",
            "Endpoint positivity on a delayed shift ray propagates on that same ray.",
            SHIFTED_LEMMA,
            "No sign hypothesis or conclusion is made for n<n0.",
        ),
        LemmaRow(
            "dcht_06_order10_specialization",
            "theorem_specialization",
            "ready_to_apply",
            "The completed order-nine heat layer specializes the shifted lemma to m=10 and n0=4.",
            ORDER10_HANDOFF,
            "Conditional only on the delayed endpoint premise at lambda=-100.",
        ),
        LemmaRow(
            "dcht_07_order11_specialization",
            "theorem_specialization",
            "ready_to_apply",
            "Once the shifted order-ten heat ray is available, the same descent specializes to m=11 and n0=4.",
            ORDER11_HANDOFF,
            "Conditional on both the shifted order-ten heat ray and the order-eleven endpoint ray.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_delayed_cooperative_heat_tail_lemma",
        "date": "2026-07-16",
        "status": "exact shifted cooperative heat-flow descent lemma",
        "proof_boundary": (
            "This proves a conditional propagation lemma and its order-ten and "
            "order-eleven n>=4 specializations. It does not prove either delayed "
            "endpoint premise, any sign for n=0,1,2,3 at negative lambda, "
            "PF-infinity, RH, or Lambda<=0."
        ),
        "sources": sources,
        "exact": {
            "fixed_order_tail": FIXED_TAIL,
            "cooperative_flow": COOPERATIVE_FLOW,
            "variation_of_constants": VARIATION,
            "finite_descent": "N=max(n0,N_m); prove n=N-1 down to n0",
            "shifted_single_layer_implication": SHIFTED_LEMMA,
            "order10_delayed_handoff": ORDER10_HANDOFF,
            "order11_delayed_handoff": ORDER11_HANDOFF,
            "unused_prefix": "the proof uses no Q_(m,n) with n<n0",
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": len(rows),
            "open_rows": 0,
            "shifted_heat_lemmas": 1,
            "order10_n4_specializations": 1,
            "order11_n4_specializations": 1,
            "endpoint_premises_proved": 0,
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_delayed_cooperative_heat_tail_lemma.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_delayed_cooperative_heat_tail_lemma.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Delayed Cooperative Heat-Tail Lemma",
        "",
        "Date: 2026-07-16",
        "",
        "Status: exact shifted heat-flow propagation lemma. This is not a proof of RH or `Lambda <= 0`.",
        "",
        "## Lemma",
        "",
        "Fix an order `m` and a starting shift `n0`. Assume the lower layer is",
        "positive throughout the heat interval, the fixed-order eventual tail",
        "holds, and `Q_(m,n)(-100)>0` for every `n>=n0`.",
        "",
        "Choose `N=max(n0,N_m)`. The eventual-tail theorem supplies positivity",
        "for every `n>=N`. At `n=N-1`, variation of constants has a positive",
        "initial term and a positive forcing term from `Q_(m,N)`. Repeat at",
        "`n=N-2,...,n0`. This is a finite descending induction.",
        "",
        "```text",
        exact["variation_of_constants"],
        exact["shifted_single_layer_implication"],
        "```",
        "",
        "The proof uses no `Q_(m,n)` with `n<n0`.",
        "",
        "## Order Ten",
        "",
        "The completed order-nine layer gives `a_(10,n)>0`, so setting `n0=4` gives",
        "",
        "```text",
        exact["order10_delayed_handoff"],
        "```",
        "",
        "This remains conditional on proving the delayed endpoint premise. It",
        "does not assign signs to shifts `0,1,2,3` at negative lambda and is not",
        "PF-infinity, RH, or `Lambda<=0`.",
        "",
        "## Order Eleven",
        "",
        "Applying the same theorem at `m=11,n0=4` gives",
        "",
        "```text",
        exact["order11_delayed_handoff"],
        "```",
        "",
        "This second specialization is conditional on the order-ten heat ray",
        "and the order-eleven endpoint ray; it does not assert either premise.",
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
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print(
        "wrote shifted cooperative heat-tail lemma and order-ten/order-eleven "
        "n>=4 specializations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
