#!/usr/bin/env python3
"""Audit Q31 margin geometry and define a saddle-normalized handoff."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR_DIR = REPO_ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR_DIR) not in sys.path:
    sys.path.insert(0, str(VENDOR_DIR))

from flint import arb, ctx


STEM = "jensen_window_pf_newman_theta_q31_margin_geometry_audit"
DEFAULT_OUT = REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
Q31_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_modular_retained_q31_interval_certificate.json"
)
Q31_CACHE = Q31_SOURCE.with_suffix(".jsonl")
ARBITRARY_N_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_arbitrary_n_stable_remainder_gate.json"
)
SADDLE_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_modular_blend_adaptive_saddle_gate.json"
)
DATE = "2026-07-24"
SCHEMA = "newman_theta_modular_retained_q31_interval_v1"
ZERO_HASH = "0" * 64
WEAKEST_COUNT = 12


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def fraction_arb(text: str) -> arb:
    value = Fraction(text)
    return arb(value.numerator) / value.denominator


def lower_float(value: arb) -> float:
    return float(value.lower())


def endpoint_abs_lower(lower_text: str, upper_text: str) -> arb:
    lower = arb(lower_text)
    upper = arb(upper_text)
    if lower.lower() > 0:
        return lower.lower()
    if upper.upper() < 0:
        return -upper.upper()
    return arb(0)


def sign_pattern(lower_text: str, upper_text: str) -> str:
    lower = arb(lower_text)
    upper = arb(upper_text)
    if lower.lower() > 0:
        return "positive"
    if upper.upper() < 0:
        return "negative"
    return "contains_zero"


def read_validated_cache(config_hash: str) -> list[dict]:
    raw = Q31_CACHE.read_bytes()
    if raw and not raw.endswith(b"\n"):
        raise RuntimeError("Q31 cache lacks final line terminator")
    records: list[dict] = []
    previous = ZERO_HASH
    for line_number, line in enumerate(
        raw.decode("utf-8").splitlines(), start=1
    ):
        record = json.loads(line)
        stored_hash = record.get("row_sha256")
        payload = dict(record)
        payload.pop("row_sha256", None)
        actual_hash = sha256(
            canonical_json(payload).encode("utf-8")
        ).hexdigest()
        if stored_hash != actual_hash:
            raise RuntimeError(f"Q31 cache hash drift at line {line_number}")
        if record.get("sequence") != line_number:
            raise RuntimeError(
                f"Q31 cache sequence drift at line {line_number}"
            )
        if record.get("schema") != SCHEMA:
            raise RuntimeError(
                f"Q31 cache schema drift at line {line_number}"
            )
        if record.get("config_sha256") != config_hash:
            raise RuntimeError(
                f"Q31 cache config drift at line {line_number}"
            )
        if record.get("previous_row_sha256") != previous:
            raise RuntimeError(
                f"Q31 cache chain drift at line {line_number}"
            )
        previous = stored_hash
        records.append(record)
    return records


def leaf_diagnostic(sequence: int, leaf: dict) -> dict:
    value_lower = endpoint_abs_lower(
        leaf["j_retained_lower"], leaf["j_retained_upper"]
    )
    derivative_lower = endpoint_abs_lower(
        leaf["j_retained_prime_lower"],
        leaf["j_retained_prime_upper"],
    )
    x_floor = fraction_arb(leaf["x_low"])
    saddle_derivative_scale = 4 * (arb.pi() * x_floor).sqrt()
    scaled_derivative_lower = saddle_derivative_scale * derivative_lower
    if value_lower.lower() >= scaled_derivative_lower.lower():
        saddle_branch = "value"
        saddle_jet_lower = value_lower
    else:
        saddle_branch = "saddle_derivative"
        saddle_jet_lower = scaled_derivative_lower
    unscaled_jet_lower = (
        value_lower
        if value_lower.lower() >= derivative_lower.lower()
        else derivative_lower
    )
    return {
        "source_sequence": sequence,
        "depth": leaf["depth"],
        "branch": leaf["branch"],
        "t_low": leaf["t_low"],
        "t_high": leaf["t_high"],
        "x_low": leaf["x_low"],
        "x_high": leaf["x_high"],
        "value_sign_pattern": sign_pattern(
            leaf["j_retained_lower"], leaf["j_retained_upper"]
        ),
        "derivative_sign_pattern": sign_pattern(
            leaf["j_retained_prime_lower"],
            leaf["j_retained_prime_upper"],
        ),
        "j_retained_lower": leaf["j_retained_lower"],
        "j_retained_upper": leaf["j_retained_upper"],
        "j_retained_prime_lower": leaf["j_retained_prime_lower"],
        "j_retained_prime_upper": leaf["j_retained_prime_upper"],
        "tail_value_upper": leaf["tail_value_upper"],
        "tail_derivative_upper": leaf["tail_derivative_upper"],
        "certified_ratio_lower": leaf["certified_ratio_lower"],
        "value_abs_lower": value_lower.str(55),
        "derivative_abs_lower": derivative_lower.str(55),
        "unscaled_first_jet_lower": unscaled_jet_lower.str(55),
        "saddle_derivative_scale_lower": (
            saddle_derivative_scale.lower().str(55)
        ),
        "saddle_derivative_lower": (
            scaled_derivative_lower.lower().str(55)
        ),
        "saddle_first_jet_branch": saddle_branch,
        "saddle_first_jet_lower": saddle_jet_lower.lower().str(55),
    }


def build_geometry(records: list[dict], source: dict) -> dict:
    ctx.prec = 192
    diagnostics: list[dict] = []
    branch_counts = {"value": 0, "derivative": 0}
    value_patterns = {"positive": 0, "negative": 0, "contains_zero": 0}
    derivative_patterns = {
        "positive": 0,
        "negative": 0,
        "contains_zero": 0,
    }
    for record in records:
        sequence = record["sequence"]
        if record["result"].get("unresolved_boxes") != 0:
            raise RuntimeError(
                f"unresolved Q31 box in source sequence {sequence}"
            )
        for leaf in record["result"]["evaluations"]:
            if not leaf.get("certified"):
                # Subdivided parent evaluations are retained in the cache.
                continue
            diagnostic = leaf_diagnostic(sequence, leaf)
            diagnostics.append(diagnostic)
            branch_counts[diagnostic["branch"]] += 1
            value_patterns[diagnostic["value_sign_pattern"]] += 1
            derivative_patterns[diagnostic["derivative_sign_pattern"]] += 1

    ratio_sorted = sorted(
        diagnostics,
        key=lambda row: lower_float(arb(row["certified_ratio_lower"])),
    )
    unscaled_sorted = sorted(
        diagnostics,
        key=lambda row: lower_float(arb(row["unscaled_first_jet_lower"])),
    )
    saddle_sorted = sorted(
        diagnostics,
        key=lambda row: lower_float(arb(row["saddle_first_jet_lower"])),
    )

    coordinate_marker = 16 * arb.pi()
    marker_rows = [
        row
        for row in diagnostics
        if fraction_arb(row["x_low"]).lower() <= coordinate_marker.upper()
        and fraction_arb(row["x_high"]).upper() >= coordinate_marker.lower()
    ]
    marker_saddle_min = min(
        marker_rows,
        key=lambda row: lower_float(arb(row["saddle_first_jet_lower"])),
    )

    weakest = ratio_sorted[0]
    if weakest["value_sign_pattern"] != "contains_zero":
        raise RuntimeError("weakest Q31 ratio no longer crosses retained zero")
    if weakest["derivative_sign_pattern"] != "positive":
        raise RuntimeError("weakest Q31 ratio lost positive derivative")

    summary = source["summary"]
    if len(diagnostics) != summary["certified_leaf_boxes"]:
        raise RuntimeError("Q31 leaf count disagrees with source summary")
    if branch_counts != summary["branch_counts"]:
        raise RuntimeError("Q31 branch counts disagree with source summary")

    saddle_low = (arb(38) / (4 * arb.pi())).sqrt()
    saddle_high = (arb(69) / (4 * arb.pi())).sqrt()
    return {
        "records": len(records),
        "leaves": len(diagnostics),
        "branch_counts": branch_counts,
        "value_sign_patterns": value_patterns,
        "derivative_sign_patterns": derivative_patterns,
        "weakest_certified_ratio_leaves": ratio_sorted[:WEAKEST_COUNT],
        "weakest_unscaled_first_jet_leaf": unscaled_sorted[0],
        "weakest_saddle_first_jet_leaf": saddle_sorted[0],
        "minimum_unscaled_first_jet_lower": (
            unscaled_sorted[0]["unscaled_first_jet_lower"]
        ),
        "minimum_saddle_first_jet_lower": (
            saddle_sorted[0]["saddle_first_jet_lower"]
        ),
        "saddle_coordinate": {
            "definition": "s=sqrt(x/(4*pi))",
            "derivative_conversion": "partial_s=4*sqrt(pi*x)*partial_x",
            "domain_lower": saddle_low.lower().str(55),
            "domain_upper": saddle_high.upper().str(55),
            "first_integer_marker": "s=2, x=16*pi",
            "integer_marker_leaf_count": len(marker_rows),
            "integer_marker_minimum_saddle_first_jet_leaf": (
                marker_saddle_min
            ),
            "exact_transition_surfaces": (
                "n_*(a,t,x)=k for a in {5,9}; x=4*pi*k^2 is only "
                "the leading asymptotic marker"
            ),
        },
        "adaptive_count_audit": {
            "retained_count_in_Q31": 7,
            "kappa_one_count_at_x38": 16,
            "kappa_one_count_at_x69": 25,
            "exact_endpoint_checks": [
                "15^4 < 39^3 < 16^4",
                "24^4 < 70^3 < 25^4",
            ],
            "conclusion": (
                "The fixed N=7 Q31 certificate does not sample even the "
                "kappa=1 absolute-tail count ceil((1+x)^(3/4))."
            ),
        },
    }


def build_artifact() -> dict:
    source = json.loads(Q31_SOURCE.read_text(encoding="utf-8"))
    if source.get("kind") != (
        "jensen_window_pf_newman_theta_modular_retained_q31_interval_certificate"
    ):
        raise RuntimeError("Q31 source kind drifted")
    if not source.get("summary", {}).get("theorem_ready"):
        raise RuntimeError("Q31 source is not theorem-ready")
    arbitrary = json.loads(
        ARBITRARY_N_SOURCE.read_text(encoding="utf-8")
    )
    if arbitrary.get("kind") != (
        "jensen_window_pf_newman_theta_arbitrary_n_stable_remainder_gate"
    ):
        raise RuntimeError("arbitrary-N source kind drifted")
    saddle = json.loads(SADDLE_SOURCE.read_text(encoding="utf-8"))
    if saddle.get("kind") != (
        "jensen_window_pf_newman_theta_modular_blend_adaptive_saddle_gate"
    ):
        raise RuntimeError("adaptive-saddle source kind drifted")
    records = read_validated_cache(source["config_sha256"])
    geometry = build_geometry(records, source)

    candidate = (
        "s=sqrt(x/(4*pi)); "
        "G_N(t,x)=max(|J_N(t,x)|,"
        "|partial_s J_N(t,x)|)/A_N(t,x)="
        "max(|J_N|,4*sqrt(pi*x)*|J_N'|)/A_N; "
        "E_N=max(E_(0,N),4*sqrt(pi*x)*E_(1,N))/A_N; "
        "prove G_N>E_N on cells split at every exact "
        "n_*(a,t,x)=k surface for a in {5,9} and at every "
        "adaptive-count jump; x=4*pi*k^2 is only the leading marker"
    )
    rows = [
        GateRow(
            id="ntq31mga_01_cache_provenance",
            role="finite_certificate_audit",
            readiness="proved",
            claim=(
                "The complete Q31 hash chain and every certified leaf were "
                "replayed as immutable inputs to this geometry audit."
            ),
            formula=(
                "1240 source records; 1566 certified leaves; "
                "previous-row and canonical-row hashes verified"
            ),
            proof_boundary=(
                "This reuses directed endpoint enclosures; it does not "
                "recompute the underlying Arb integrals."
            ),
        ),
        GateRow(
            id="ntq31mga_02_weakest_leaf_transversality",
            role="finite_exact_diagnostic",
            readiness="proved",
            claim=(
                "The weakest reported Q31 ratio is a retained zero-crossing "
                "box certified by a strictly positive derivative."
            ),
            formula=(
                "0 in J_7(B), inf J_7'(B)>0, "
                "certified_ratio>1.454773029e20"
            ),
            proof_boundary=(
                "The huge ratio is tail-budget relative, not an intrinsic "
                "cofinal lower-margin theorem."
            ),
            diagnostics=geometry["weakest_certified_ratio_leaves"][0],
        ),
        GateRow(
            id="ntq31mga_03_saddle_coordinate_jet",
            role="exact_identity",
            readiness="proved",
            claim=(
                "Changing from frequency x to saddle-count coordinate s "
                "rescales the derivative by one unit of saddle motion."
            ),
            formula=(
                "s=sqrt(x/(4*pi)); x=4*pi*s^2; "
                "partial_s=8*pi*s*partial_x="
                "4*sqrt(pi*x)*partial_x"
            ),
            proof_boundary=(
                "This is a coordinate identity, not a lower bound."
            ),
        ),
        GateRow(
            id="ntq31mga_04_finite_saddle_jet_floor",
            role="finite_exact_diagnostic",
            readiness="proved",
            claim=(
                "The Q31 endpoint enclosures give a positive finite floor "
                "for max(|J_7|,|partial_s J_7|)."
            ),
            formula=(
                "min_B max(abs_lower(J_7(B)),"
                "4*sqrt(pi*x_low(B))*abs_lower(J_7'(B)))"
            ),
            proof_boundary=(
                "This floor is only on the fixed-N=7 finite Q31 cover and "
                "does not extrapolate to x>69."
            ),
            diagnostics=geometry["weakest_saddle_first_jet_leaf"],
        ),
        GateRow(
            id="ntq31mga_05_raw_ratio_extrapolation_guard",
            role="rejected_inference",
            readiness="proved",
            claim=(
                "The Q31 minimum tail ratio is not a valid empirical proxy "
                "for an intrinsic cofinal transversality constant."
            ),
            formula=(
                "ratio=max(abs_lower(J_7)/E_0,"
                "abs_lower(J_7')/E_1); E_0,E_1 are about 10^-23"
            ),
            proof_boundary=(
                "Changing the tail budget changes this ratio without "
                "changing the retained zero slope."
            ),
        ),
        GateRow(
            id="ntq31mga_06_adaptive_count_noncoverage",
            role="finite_exact_diagnostic",
            readiness="proved",
            claim=(
                "Q31 does not sample the proposed x^(3/4) adaptive-count "
                "regime, even at kappa=1."
            ),
            formula=(
                "N_Q31=7; ceil(39^(3/4))=16; "
                "ceil(70^(3/4))=25"
            ),
            proof_boundary=(
                "Therefore Q31 can anchor a cofinal cover but cannot "
                "validate its count scaling."
            ),
            diagnostics=geometry["adaptive_count_audit"],
        ),
        GateRow(
            id="ntq31mga_07_saddle_normalized_candidate",
            role="proof_search_target",
            readiness="not_ready_to_apply",
            claim=(
                "Prove a saddle-normalized retained first-jet margin against "
                "the arbitrary-N full error bars on an exact-saddle-aligned "
                "cofinal partition."
            ),
            formula=candidate,
            proof_boundary=(
                "A_N must be an explicit positive Xi/Riemann-Siegel "
                "amplitude, the switch-defect and full d0/d1 errors must be "
                "closed for arbitrary N, and the cells must terminate. No "
                "such theorem, Lambda<=0, RH, or Clay-prize conclusion is "
                "claimed here."
            ),
        ),
    ]
    return {
        "kind": STEM,
        "date": DATE,
        "status": "rigorous finite geometry audit and cofinal candidate",
        "source_audit": {
            "q31_artifact_sha256": file_hash(Q31_SOURCE),
            "q31_cache_sha256": file_hash(Q31_CACHE),
            "q31_config_sha256": source["config_sha256"],
            "q31_last_row_sha256": records[-1]["row_sha256"],
            "arbitrary_n_artifact_sha256": file_hash(ARBITRARY_N_SOURCE),
            "adaptive_saddle_artifact_sha256": file_hash(SADDLE_SOURCE),
        },
        "geometry": geometry,
        "candidate": candidate,
        "rows": [asdict(row) for row in rows],
        "proof_boundary": (
            "This artifact rigorously derives finite geometry from the "
            "completed Q31 endpoint enclosures and proposes a precise "
            "saddle-normalized cofinal target. It does not recompute Q31 "
            "integrals, prove an arbitrary-N switch-defect budget, provide "
            "an explicit amplitude lower theorem, terminate the cover for "
            "x>69, prove Lambda<=0 or RH, or supply a Clay-prize proof."
        ),
    }


def render_note(artifact: dict) -> str:
    geometry = artifact["geometry"]
    weakest = geometry["weakest_certified_ratio_leaves"][0]
    saddle_weakest = geometry["weakest_saddle_first_jet_leaf"]
    lines = [
        "# Newman Theta Q31 Margin Geometry Audit",
        "",
        f"Date: {DATE}",
        "",
        "Status: rigorous finite audit and an unproved cofinal candidate.",
        "This is not an extrapolation of Q31 and not a proof of `Lambda<=0`,",
        "RH, or a Clay-prize result.",
        "",
        "## What The Weakest Box Means",
        "",
        f"The complete cache has `{geometry['records']}` records and",
        f"`{geometry['leaves']}` certified leaves. The weakest reported",
        "tail-relative ratio occurs on",
        "",
        "```text",
        f"t in [{weakest['t_low']},{weakest['t_high']}]",
        f"x in [{weakest['x_low']},{weakest['x_high']}]",
        f"J_7 lower endpoint = {weakest['j_retained_lower']}",
        f"J_7 upper endpoint = {weakest['j_retained_upper']}",
        f"J_7' lower endpoint = {weakest['j_retained_prime_lower']}",
        f"J_7' upper endpoint = {weakest['j_retained_prime_upper']}",
        f"ratio lower = {weakest['certified_ratio_lower']}",
        "```",
        "",
        "The retained value interval crosses zero, while the retained",
        "derivative is strictly positive. The ratio is enormous because the",
        "full-tail derivative bar is around `10^-23`; it is not evidence for",
        "an intrinsic margin of size `10^20`.",
        "",
        "## Saddle Coordinate",
        "",
        "Set",
        "",
        "```text",
        "s=sqrt(x/(4*pi)),  x=4*pi*s^2,",
        "partial_s=4*sqrt(pi*x)*partial_x.",
        "```",
        "",
        "This compares a value with the change produced by moving one unit",
        "in saddle count. On the fixed Q31 cover, endpoint enclosures give",
        "",
        "```text",
        "min max(|J_7|_lower,4*sqrt(pi*x_low)*|J_7'|_lower)",
        f"  = {geometry['minimum_saddle_first_jet_lower']}",
        f"weakest box: t=[{saddle_weakest['t_low']},",
        f"{saddle_weakest['t_high']}], x=[{saddle_weakest['x_low']},",
        f"{saddle_weakest['x_high']}].",
        "```",
        "",
        "That is a finite diagnostic only.",
        "",
        "## Why Q31 Cannot Set The Cofinal Scaling",
        "",
        "Q31 keeps `N=7`. Even the `kappa=1` absolute-tail count",
        "`ceil((1+x)^(3/4))` runs from 16 to 25 across this slab. Thus Q31",
        "does not sample the adaptive-count regime; it can only anchor it.",
        "",
        "## Candidate Gate",
        "",
        "```text",
        artifact["candidate"],
        "```",
        "",
        "The amplitude `A_N` must come from an explicit Xi-specific",
        "saddle or corrected Riemann-Siegel theorem. The partition must split",
        "at every exact `n_*(a,t,x)=k` surface for `a=5,9` and every",
        "adaptive-count jump. The simpler `x=4*pi*k^2` locations are only",
        "leading markers. The partition must terminate",
        "after a proved asymptotic region. Until those pieces exist, this is",
        "a precise proof-search target rather than a theorem.",
        "",
    ]
    return "\n".join(lines)


def write_artifact(artifact: dict, out: Path, note: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    note.write_text(render_note(artifact), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    write_artifact(artifact, args.out, args.note)
    geometry = artifact["geometry"]
    print(
        "wrote Q31 margin geometry audit: "
        f"{geometry['records']} records, {geometry['leaves']} leaves, "
        f"{len(artifact['rows'])} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
