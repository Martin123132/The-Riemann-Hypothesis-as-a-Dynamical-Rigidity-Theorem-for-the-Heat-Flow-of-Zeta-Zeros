#!/usr/bin/env python3
"""Independently check the saved-height signed one-sided closure gate."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
for directory in (SCRIPT_DIR, VENDOR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx

import check_jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_eight_round_remainder_gate as interval_check


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_R_after_A_one_sided_closure_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing artifact: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def rebuild_b_trace(deps: dict[str, dict[str, Any]]) -> arb:
    local = deps["B_local"]["interval_certificate"]
    window = deps["B_window"]["interval_certificate"]
    principal = arb(local["combined_principal_physical_signed_projection_ball"])
    signed_first = arb(window["signed_nonlocal_first_current_ball"])
    # Reverse the production order deliberately.
    allowances = [
        arb(window["Fresnel_boundary_dictionary_absolute_allowance_ball"]).upper(),
        arb(window["nonlocal_post_third_current_absolute_allowance_ball"]).upper(),
        arb(window["nonlocal_higher_current_absolute_allowance_ball"]).upper(),
        arb(local["negative_partner_remainder_absolute_allowance_ball"]).upper(),
        arb(local["positive_normal_remainder_absolute_allowance_ball"]).upper(),
    ]
    radius = sum(allowances, arb(0))
    return signed_first + arb(0, radius) + principal


def main() -> int:
    resource_mode = interval_check.set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = 512
    ctx.threads = 1
    for path in (RESULT, NOTE, BUILDER, CHECKER):
        require(path.is_file(), f"missing artifact: {path}")

    artifact = load_json(RESULT)
    require(artifact.get("passed") is True, "production gate not passed")
    decision = artifact["decision"]
    require(decision["original_saved_height_one_sided_route_closed"] is True, "closure decision drift")
    require(decision["Q_K_minus_T_negative_at_saved_height"] is True, "Q_K-T sign drift")
    require(decision["all_height_transport_theorem_proved"] is False, "all-height overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")

    for source in artifact["sources"].values():
        path = ROOT / source["path"]
        require(path.is_file() and file_hash(path) == source["sha256"], f"source hash drift: {path}")
    deps: dict[str, dict[str, Any]] = {}
    for name, source in artifact["dependencies"].items():
        path = ROOT / source["path"]
        require(path.is_file() and file_hash(path) == source["sha256"], f"dependency hash drift: {path}")
        deps[name] = load_json(path)
    require(all(dep.get("passed") is True for dep in deps.values()), "dependency pass drift")

    coefficients = artifact["certificate"]["coefficient_cancellation"]
    require(coefficients["A_face_coefficient_after_collection"] == 0, "A-face did not cancel")
    require(coefficients["A_transition_coefficient_after_collection"] == 0, "A-transition did not cancel")
    require(
        coefficients["common_exact_result"] == "R_after_A=J_Z-E_Btr,win-E_outer",
        "collected exact identity drift",
    )

    saved = deps["signed_JZ"]["certificate"]
    J_Z = arb(saved["J_Z_ball"]["ball"])
    R_KGamma = arb(saved["R_KGamma_ball"]["ball"])
    QKT_saved = arb(saved["Q_K_minus_T_ball"]["ball"])
    B_trace = rebuild_b_trace(deps)
    require(
        B_trace.overlaps(arb(artifact["certificate"]["B_trace_ball"]["ball"])),
        "changed-order B trace misses production",
    )

    outer_bound = arb(
        deps["B_outer"]["interval_certificate"]["physical_outer_bound_ball"]
    ).upper()
    E_outer = arb(0, outer_bound)
    A_transition = arb(
        deps["A_transition"]["certificate"]["components"][
            "projector_completed_exact_A_transition_ball"
        ]
    )

    # Use different groupings from the production builder.
    R_after_A = (J_Z - E_outer) - B_trace
    R_Dir_route_one = R_after_A + A_transition
    R_Dir_route_two = (R_KGamma - E_outer) - B_trace
    require(R_Dir_route_one.overlaps(R_Dir_route_two), "independent R_Dir routes miss")

    gamma_bound = abs(
        arb(
            deps["Gamma_insertion"]["correction_certificate"][
                "Gamma_to_classical_physical_absolute_bound_ball"
            ]
        )
    ).upper()
    G_minus_T = arb(0, gamma_bound)
    QKT_direct = G_minus_T + R_KGamma
    QKT_rejoined = G_minus_T + R_Dir_route_one + E_outer + B_trace
    require(QKT_direct.overlaps(QKT_saved), "direct Q_K-T misses saved enclosure")
    require(QKT_rejoined.overlaps(QKT_saved), "rejoined Q_K-T misses saved enclosure")

    for name, value in (
        ("R_after_A_ball", R_after_A),
        ("R_Dir_ball", R_Dir_route_one),
        ("R_Dir_alternate_ball", R_Dir_route_two),
        ("Q_K_minus_T_ball", QKT_direct),
        ("Q_K_minus_T_rejoined_ball", QKT_rejoined),
    ):
        require(
            value.overlaps(arb(artifact["certificate"][name]["ball"])),
            f"changed-order {name} misses production",
        )

    require(R_after_A.upper() < arb("0.0368147039947"), "working R_after_A target failed")
    require(R_after_A.upper() < arb("0.0368061039947"), "strict R_after_A target failed")
    require(R_Dir_route_one.upper() < arb("0.00013197919999999995"), "R_Dir target failed")
    require(R_Dir_route_one.upper() < 0, "R_Dir negativity failed")
    require(QKT_direct.upper() < 0 and QKT_rejoined.upper() < 0, "Q_K-T negativity failed")

    note = " ".join(NOTE.read_text(encoding="utf-8").split())
    for fragment in (
        "R_after_A=J_Z-E_Btr,win-E_outer",
        "A-face channel exactly",
        "No sign or quadrature estimate",
        "Q_K-T",
        "single saved height",
        "does not prove the stronger absolute non-A bound",
        "not an RH proof",
    ):
        require(fragment.lower() in note.lower(), f"missing note fragment: {fragment}")

    print(
        "independently checked saved-height signed one-sided closure and Q_K-T negativity; "
        f"R_after_A={R_after_A.str(12, more=True)}, Q_K-T={QKT_direct.str(12, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
